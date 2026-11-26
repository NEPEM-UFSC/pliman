
pkgload::load_all("d:/Desktop/pliman")

vid_path <- "D:/Downloads/pliman4/toras2.mp4"
model_path <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"
lib_path <- pliman_onnx_library_path()

v_reader <- pliman:::open_native_video_reader_cpp(normalizePath(vid_path, winslash = "/"))
info <- pliman:::get_video_reader_info_cpp(v_reader)
src_w <- info$width
src_h <- info$height
bm_stream <- raw(3L * src_w * src_h)
attr(bm_stream, "dim") <- c(3L, src_w, src_h)

out_vid <- "D:/Downloads/pliman4/scratch_test_out.mp4"
v_writer <- pliman:::open_native_video_writer_cpp(out_vid, src_w, src_h, 30.0)

n_test <- 50
cat(sprintf("Profiling %d frames of video_detect_dl pipeline on %dx%d video...\n", n_test, src_w, src_h))

t_read <- 0
t_crop <- 0
t_preproc <- 0
t_onnx <- 0
t_track <- 0
t_draw <- 0
t_write <- 0

# Warm up ONNX session on GPU
dummy_tensor <- pliman:::preprocess_yolo_cpp(bm_stream, target_size = 640L)
pliman:::run_yolo_cpp(dummy_tensor, orig_w = 640, orig_h = 640, model_path = model_path, lib_path = lib_path, use_gpu = TRUE)

track_state <- list(next_id = 1L, total_count = 0L, tracks = list())

roi_coords <- c(as.integer(0.5 * src_w), 1L, as.integer(src_w), as.integer(src_h))
line_coords <- c(as.integer(0.5 * src_w), 1L, as.integer(0.5 * src_w), as.integer(src_h))

t_total_start <- Sys.time()

for (i in 1:n_test) {
  # 1. Read
  t0 <- Sys.time()
  has_frame <- pliman:::read_native_video_frame_cpp(v_reader, bm_stream)
  if (!has_frame) break
  t_read <- t_read + as.numeric(Sys.time() - t0, units = "secs")

  # 2. Crop
  t0 <- Sys.time()
  margin_x <- round(0.08 * src_w)
  margin_y <- round(0.08 * src_h)
  crop_x1 <- max(0L, as.integer(roi_coords[1] - margin_x))
  crop_y1 <- max(0L, as.integer(roi_coords[2] - margin_y))
  crop_x2 <- min(src_w, as.integer(roi_coords[3] + margin_x))
  crop_y2 <- min(src_h, as.integer(roi_coords[4] + margin_y))
  crop_coords <- c(crop_x1, crop_y1, crop_x2, crop_y2)
  bm_feed <- pliman:::crop_bgr_cpp(bm_stream, crop_x1, crop_y1, crop_x2, crop_y2)
  feed_dims <- dim(bm_feed)
  feed_w <- feed_dims[2]
  feed_h <- feed_dims[3]
  t_crop <- t_crop + as.numeric(Sys.time() - t0, units = "secs")

  # 3. Preprocess
  t0 <- Sys.time()
  tensor <- pliman:::preprocess_yolo_cpp(bm_feed, target_size = 640L)
  t_preproc <- t_preproc + as.numeric(Sys.time() - t0, units = "secs")

  # 4. ONNX Inference
  t0 <- Sys.time()
  raw_res <- pliman:::run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = feed_w,
    orig_h = feed_h,
    conf_threshold = 0.5,
    iou_threshold = 0.45,
    model_path = model_path,
    lib_path = lib_path,
    num_threads = 0L,
    use_gpu = TRUE,
    device_id = -1L,
    need_masks = FALSE
  )
  t_onnx <- t_onnx + as.numeric(Sys.time() - t0, units = "secs")

  # Shift box coordinates
  if (nrow(raw_res$boxes) > 0) {
    raw_res$boxes[, 1] <- raw_res$boxes[, 1] + crop_coords[1]
    raw_res$boxes[, 2] <- raw_res$boxes[, 2] + crop_coords[2]
    raw_res$boxes[, 3] <- raw_res$boxes[, 3] + crop_coords[1]
    raw_res$boxes[, 4] <- raw_res$boxes[, 4] + crop_coords[2]
  }

  num_boxes <- nrow(raw_res$boxes)
  df_b <- if (num_boxes > 0) {
    data.frame(
      xmin = raw_res$boxes[, 1], ymin = raw_res$boxes[, 2],
      xmax = raw_res$boxes[, 3], ymax = raw_res$boxes[, 4],
      label = rep("wood_log", num_boxes), score = raw_res$scores,
      class_id = raw_res$class_ids, stringsAsFactors = FALSE
    )
  } else {
    data.frame(xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
               label = character(0), score = numeric(0), class_id = integer(0))
  }

  # 5. Tracker
  t0 <- Sys.time()
  if (nrow(df_b) > 0) {
    track_res <- pliman:::update_tracker_cpp(
      xmin = df_b$xmin, ymin = df_b$ymin, xmax = df_b$xmax, ymax = df_b$ymax,
      labels = df_b$label, scores = df_b$score, state = track_state,
      count_line = line_coords, roi = roi_coords, max_dist = 150.0, min_iou = 0.25,
      max_lost = 15L, max_history = 30L
    )
    track_state <- track_res$state
    track_ids_pass <- track_res$track_ids
  } else {
    track_ids_pass <- NULL
  }
  t_track <- t_track + as.numeric(Sys.time() - t0, units = "secs")

  # 6. Draw annotations
  t0 <- Sys.time()
  pliman:::draw_yolo_detections_bgr_cpp(
    bm = bm_stream,
    xmin = df_b$xmin, ymin = df_b$ymin, xmax = df_b$xmax, ymax = df_b$ymax,
    labels = df_b$label, scores = df_b$score,
    colors = rep("#00CC66", nrow(df_b)), lwd = 2L, hud_text = "50 FPS",
    font_scale = 1.0, keypoints = NULL, kpt_threshold = 0.3, kpt_radius = 4L,
    draw_skeleton = FALSE, draw_boxes = TRUE, track_ids = track_ids_pass,
    flash = NULL, history_x = NULL, history_y = NULL, roi = roi_coords,
    roi_label = "ROI", count_line = line_coords, count_line_label = "",
    counted = NULL, hide_outside_roi = TRUE, hide_counted = FALSE,
    flash_counted = FALSE, show_text = FALSE, show_conf = FALSE,
    show_class = TRUE, show_id = FALSE, badge = TRUE, pad_scale = 1.0,
    mask_labels = NULL, mask_alpha = 0.4, mask_offset_x = 0L, mask_offset_y = 0L,
    draw_masks = FALSE, mask_ids = NULL, hud_pos = "top-left", hud_layout = "vertical",
    font_scale_hud = 1.0, font_scale_roi = 1.0
  )
  t_draw <- t_draw + as.numeric(Sys.time() - t0, units = "secs")

  # 7. Write frame
  t0 <- Sys.time()
  pliman:::write_native_video_frame_cpp(v_writer, bm_stream)
  t_write <- t_write + as.numeric(Sys.time() - t0, units = "secs")
}

t_total <- as.numeric(Sys.time() - t_total_start, units = "secs")
pliman:::close_native_video_reader_cpp(v_reader)
pliman:::close_native_video_writer_cpp(v_writer)
unlink(out_vid)

cat(sprintf("\n=== Breakdown for %d frames (Total: %.3f s | %.1f FPS | %.2f ms/frame) ===\n",
            n_test, t_total, n_test / t_total, (t_total / n_test) * 1000))
cat(sprintf("1. Read Video Frame:        %.3f s (%.1f%%) | %.2f ms/frame\n", t_read, (t_read/t_total)*100, (t_read/n_test)*1000))
cat(sprintf("2. ROI Crop:                %.3f s (%.1f%%) | %.2f ms/frame\n", t_crop, (t_crop/t_total)*100, (t_crop/n_test)*1000))
cat(sprintf("3. Preprocess YOLO (640):   %.3f s (%.1f%%) | %.2f ms/frame\n", t_preproc, (t_preproc/t_total)*100, (t_preproc/n_test)*1000))
cat(sprintf("4. ONNX Inference (GPU):    %.3f s (%.1f%%) | %.2f ms/frame\n", t_onnx, (t_onnx/t_total)*100, (t_onnx/n_test)*1000))
cat(sprintf("5. Tracker (C++):           %.3f s (%.1f%%) | %.2f ms/frame\n", t_track, (t_track/t_total)*100, (t_track/n_test)*1000))
cat(sprintf("6. Draw Annotations (C++):  %.3f s (%.1f%%) | %.2f ms/frame\n", t_draw, (t_draw/t_total)*100, (t_draw/n_test)*1000))
cat(sprintf("7. Write Video Frame:       %.3f s (%.1f%%) | %.2f ms/frame\n", t_write, (t_write/t_total)*100, (t_write/n_test)*1000))
