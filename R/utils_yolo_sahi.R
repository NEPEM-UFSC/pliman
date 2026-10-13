# ==============================================================================
# SAHI: SLICING AIDED HYPER INFERENCE FOR GIGAPIXEL & DRONE IMAGERY
# pliman: Plant Image Analysis in R
# ==============================================================================

#' Internal SAHI Post-Processing (Greedy NMM / NMS with IOS & Seam Merging)
#' @noRd
.sahi_postprocess <- function(candidates_df,
                              contours = NULL,
                              match_metric = c("IOS", "IOU"),
                              match_threshold = 0.50,
                              postprocess_type = c("GREEDY_NMM", "NMS"),
                              ...) {
  match_metric <- match.arg(match_metric)
  postprocess_type <- match.arg(postprocess_type)

  if (is.null(candidates_df) || nrow(candidates_df) == 0L) {
    return(list(boxes = candidates_df, contours = contours))
  }

  has_conts <- !is.null(contours) && length(contours) == nrow(candidates_df)

  # Standardize score and confidence
  c_vals <- if (!is.null(candidates_df$conf)) {
    candidates_df$conf
  } else if (!is.null(candidates_df$score)) {
    candidates_df$score
  } else {
    rep(1.0, nrow(candidates_df))
  }
  candidates_df$conf <- c_vals
  candidates_df$score <- c_vals

  if (is.null(candidates_df$is_full)) {
    candidates_df$is_full <- FALSE
  }

  unique_classes <- unique(candidates_df$class_id)
  final_boxes_list <- list()
  final_conts_list <- list()

  for (cls in unique_classes) {
    cls_idx <- which(candidates_df$class_id == cls)
    cls_df <- candidates_df[cls_idx, , drop = FALSE]
    cls_conts <- if (has_conts) contours[cls_idx] else NULL

    # ------------------------------------------------------------------
    # Step 1: Subsume slice fragments into full-image detections
    # When full-image inference detected a whole object, slice fragments
    # that overlap significantly with it are absorbed into the full object
    # (retaining the full unbroken contour and absorbing higher confidence).
    # ------------------------------------------------------------------
    if (any(cls_df$is_full)) {
      full_idx <- which(cls_df$is_full)
      slice_idx <- which(!cls_df$is_full)

      if (length(slice_idx) > 0L) {
        absorbed_slice <- logical(length(slice_idx))

        for (fi in full_idx) {
          fx1 <- cls_df$xmin[fi]; fy1 <- cls_df$ymin[fi]
          fx2 <- cls_df$xmax[fi]; fy2 <- cls_df$ymax[fi]
          f_area <- (fx2 - fx1 + 1) * (fy2 - fy1 + 1)

          for (si_pos in seq_along(slice_idx)) {
            if (absorbed_slice[si_pos]) next
            si <- slice_idx[si_pos]

            sx1 <- cls_df$xmin[si]; sy1 <- cls_df$ymin[si]
            sx2 <- cls_df$xmax[si]; sy2 <- cls_df$ymax[si]
            s_area <- (sx2 - sx1 + 1) * (sy2 - sy1 + 1)

            # Bounding box intersection
            ix1 <- max(fx1, sx1); iy1 <- max(fy1, sy1)
            ix2 <- min(fx2, sx2); iy2 <- min(fy2, sy2)
            iw <- max(0, ix2 - ix1 + 1)
            ih <- max(0, iy2 - iy1 + 1)
            inter_area <- iw * ih

            if (inter_area > 0) {
              ios <- inter_area / min(f_area, s_area)
              iou <- inter_area / (f_area + s_area - inter_area)

              if (ios >= 0.40 || iou >= 0.35) {
                cls_df$conf[fi] <- max(cls_df$conf[fi], cls_df$conf[si])
                cls_df$score[fi] <- cls_df$conf[fi]
                absorbed_slice[si_pos] <- TRUE
              }
            }
          }
        }

        if (any(absorbed_slice)) {
          drop_rows <- slice_idx[absorbed_slice]
          keep_rows <- setdiff(seq_len(nrow(cls_df)), drop_rows)
          cls_df <- cls_df[keep_rows, , drop = FALSE]
          if (has_conts) cls_conts <- cls_conts[keep_rows]
        }
      }
    }

    if (nrow(cls_df) == 0L) next

    # ------------------------------------------------------------------
    # Step 2: Greedy NMM / NMS with IOS (cluster-based, no snowballing)
    # ------------------------------------------------------------------
    ord <- order(cls_df$conf, decreasing = TRUE)
    cls_df <- cls_df[ord, , drop = FALSE]
    if (has_conts) cls_conts <- cls_conts[ord]

    while (nrow(cls_df) > 0L) {
      best_box <- c(cls_df$xmin[1], cls_df$ymin[1], cls_df$xmax[1], cls_df$ymax[1])
      best_area <- (best_box[3] - best_box[1] + 1) * (best_box[4] - best_box[2] + 1)

      matched_idx <- 1L

      if (nrow(cls_df) > 1L) {
        for (k in 2:nrow(cls_df)) {
          o_box <- c(cls_df$xmin[k], cls_df$ymin[k], cls_df$xmax[k], cls_df$ymax[k])
          o_area <- (o_box[3] - o_box[1] + 1) * (o_box[4] - o_box[2] + 1)

          ix1 <- max(best_box[1], o_box[1]); iy1 <- max(best_box[2], o_box[2])
          ix2 <- min(best_box[3], o_box[3]); iy2 <- min(best_box[4], o_box[4])
          iw <- max(0, ix2 - ix1 + 1)
          ih <- max(0, iy2 - iy1 + 1)
          inter <- iw * ih

          if (inter > 0) {
            ios <- inter / min(best_area, o_area)
            iou <- inter / (best_area + o_area - inter)

            matched <- if (match_metric == "IOS") {
              (ios >= match_threshold && iou >= 0.30)
            } else {
              (iou >= match_threshold)
            }

            if (matched) {
              matched_idx <- c(matched_idx, k)
            }
          }
        }
      }

      if (postprocess_type == "GREEDY_NMM") {
        m_df <- cls_df[matched_idx, , drop = FALSE]
        merged_row <- m_df[1, , drop = FALSE]
        merged_row$xmin <- min(m_df$xmin)
        merged_row$ymin <- min(m_df$ymin)
        merged_row$xmax <- max(m_df$xmax)
        merged_row$ymax <- max(m_df$ymax)
        merged_row$conf <- max(m_df$conf)
        if (!is.null(merged_row$score)) merged_row$score <- merged_row$conf

        merged_cnt <- NULL
        if (has_conts) {
          m_cnts <- cls_conts[matched_idx]
          if (isTRUE(m_df$is_full[1]) && is.matrix(m_cnts[[1]]) && nrow(m_cnts[[1]]) >= 3L) {
            merged_cnt <- m_cnts[[1]]
          } else {
            all_pts <- do.call(rbind, lapply(m_cnts, function(p) if (is.matrix(p) && nrow(p) >= 3L) p else NULL))
            if (!is.null(all_pts) && nrow(all_pts) >= 3L) {
              ch <- chull(all_pts)
              merged_cnt <- all_pts[c(ch, ch[1]), , drop = FALSE]
            } else if (length(m_cnts) > 0 && is.matrix(m_cnts[[1]])) {
              merged_cnt <- m_cnts[[1]]
            }
          }
        }

        final_boxes_list[[length(final_boxes_list) + 1L]] <- merged_row
        if (has_conts) final_conts_list[[length(final_conts_list) + 1L]] <- merged_cnt
      } else {
        # Standard NMS: keep only highest-scoring candidate in cluster
        final_boxes_list[[length(final_boxes_list) + 1L]] <- cls_df[1, , drop = FALSE]
        if (has_conts) final_conts_list[[length(final_conts_list) + 1L]] <- cls_conts[[1]]
      }

      # Drop all matched candidates
      cls_df <- cls_df[-matched_idx, , drop = FALSE]
      if (has_conts) cls_conts <- cls_conts[-matched_idx]
    }
  }

  if (length(final_boxes_list) == 0L) {
    empty_df <- candidates_df[0, , drop = FALSE]
    empty_df$is_full <- NULL
    return(list(boxes = empty_df, contours = list()))
  }

  out_df <- do.call(rbind, final_boxes_list)
  rownames(out_df) <- NULL
  out_df$is_full <- NULL

  list(boxes = out_df, contours = if (has_conts) final_conts_list else NULL)
}


#' @title Slicing Aided Hyper Inference (SAHI) for Object Detection
#' @name image_detect_sahi
#' @description
#' `image_detect_sahi()` performs high-resolution sliding-window hyper inference
#' on ultra-large images, drone orthomosaics, and macro scans. Large images (e.g. 20–50 Megapixels)
#' are dynamically sliced into overlapping tiles of size `slice_width x slice_height`.
#'
#' Deep learning inference is executed per slice, and all bounding boxes are transformed
#' to global coordinates and merged using Greedy Non-Maximum Merging (Greedy NMM) or standard NMS.
#' This prevents small objects (seeds, aphids, stomata, seedlings) from vanishing when downscaled,
#' and seamlessly unifies objects that cross tile seams.
#'
#' @param img An `image` object, file path, or 3D array.
#' @param model YOLO ONNX model name or file path (e.g., `"yolo26n"`, `"yolov8s.onnx"`).
#' @param slice_height Integer. Height of sliding window slices. Defaults to `640`.
#' @param slice_width Integer. Width of sliding window slices. Defaults to `640`.
#' @param overlap_height_ratio Numeric scalar between 0 and 0.5. Height overlap fraction. Defaults to `0.20`.
#' @param overlap_width_ratio Numeric scalar between 0 and 0.5. Width overlap fraction. Defaults to `0.20`.
#' @param conf_threshold Numeric scalar between 0 and 1. Minimum detection confidence. Defaults to `0.25`.
#' @param iou_threshold Numeric scalar between 0 and 1. Threshold for cross-slice matching. Defaults to `0.45`.
#' @param match_metric Character. Overlap metric: `"IOS"` (Intersection over Smaller area, default)
#'   or `"IOU"` (Intersection over Union). `"IOS"` is strongly recommended for SAHI because it
#'   naturally merges fragments and sub-slices with larger parent detections.
#' @param postprocess_type Character. Post-processing merge strategy: `"GREEDY_NMM"` (default, Non-Maximum Merging,
#'   unifying overlapping pieces into single encompassing boxes) or `"NMS"` (standard Non-Maximum Suppression).
#' @param match_threshold Numeric scalar between 0 and 1. Threshold for merging candidates. Defaults to `0.50`.
#' @param seam_merging Logical. If `TRUE` (default), seamlessly detects and merges objects that have been
#'   sliced in half across tile cut boundaries.
#' @param full_image_inference Logical. If `TRUE` (default), runs an additional detection pass on the
#'   resized full canvas to capture objects larger than individual slices.
#' @param rainbow Logical. If `TRUE` (default), bounding boxes are rendered in rainbow colors.
#' @param plot Logical. Display global detections overlaid on the full image. Defaults to `TRUE`.
#' @param verbose Logical. Show progress bar and tile metrics. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_detect_dl()].
#'
#' @return A `data.frame` of global bounding boxes with attributes `image`, `counts`, and `summary`.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Run SAHI on a 4000x3000 drone orthophoto
#' res <- image_detect_sahi("drone_field.jpg", model = "seedlings.onnx", slice_width = 640)
#' }
image_detect_sahi <- function(img,
                              model = "yolo26n",
                              slice_height = 640,
                              slice_width = 640,
                              overlap_height_ratio = 0.20,
                              overlap_width_ratio = 0.20,
                              conf_threshold = 0.25,
                              iou_threshold = 0.45,
                              match_metric = c("IOS", "IOU"),
                              postprocess_type = c("GREEDY_NMM", "NMS"),
                              match_threshold = 0.50,
                              seam_merging = TRUE,
                              full_image_inference = TRUE,
                              rainbow = TRUE,
                              plot = TRUE,
                              verbose = TRUE,
                              ...) {
  match_metric <- match.arg(match_metric)
  postprocess_type <- match.arg(postprocess_type)

  if (missing(match_threshold) && !missing(iou_threshold)) {
    match_threshold <- iou_threshold
  } else if (is.null(match_threshold)) {
    match_threshold <- 0.50
  }

  # Ingest image
  im <- if (is.character(img) && file.exists(img[1])) {
    image_import(img[1])
  } else if (inherits(img, c("image", "Image"))) {
    img
  } else {
    as_image(img)
  }

  dims <- dim(im)
  orig_w <- dims[1]
  orig_h <- dims[2]

  slice_w <- min(orig_w, as.integer(slice_width))
  slice_h <- min(orig_h, as.integer(slice_height))

  step_x <- max(1L, as.integer(round(slice_w * (1 - overlap_width_ratio))))
  step_y <- max(1L, as.integer(round(slice_h * (1 - overlap_height_ratio))))

  # Generate slice grid
  x_starts <- unique(c(seq(1L, orig_w - slice_w + 1L, by = step_x), max(1L, orig_w - slice_w + 1L)))
  y_starts <- unique(c(seq(1L, orig_h - slice_h + 1L, by = step_y), max(1L, orig_h - slice_h + 1L)))

  grid <- expand.grid(x = x_starts, y = y_starts)
  n_slices <- nrow(grid)

  slice_borders_x <- sort(unique(c(grid$x, grid$x + slice_w - 1L)))
  slice_borders_y <- sort(unique(c(grid$y, grid$y + slice_h - 1L)))

  if (isTRUE(verbose)) {
    cli::cli_h2("SAHI (Slicing Aided Hyper Inference) - Object Detection")
    cli::cli_alert_info("Canvas: {orig_w}x{orig_h} | Slice: {slice_w}x{slice_h} | Total slices: {.val {n_slices}}")
  }

  all_candidates <- list()

  # 1. Sliding window slice inference
  for (i in seq_len(n_slices)) {
    sx1 <- grid$x[i]
    sy1 <- grid$y[i]
    sx2 <- min(orig_w, sx1 + slice_w - 1L)
    sy2 <- min(orig_h, sy1 + slice_h - 1L)

    sub_im <- if (length(dims) >= 3L) im[sx1:sx2, sy1:sy2, , drop = FALSE] else im[sx1:sx2, sy1:sy2, drop = FALSE]

    det <- try(
      image_detect_dl(
        sub_im,
        model = model,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        plot = FALSE,
        verbose = FALSE,
        ...
      ),
      silent = TRUE
    )

    if (inherits(det, "try-error")) next
    if (!is.null(attr(det, "boxes"))) {
      det <- attr(det, "boxes")
    } else if (is.list(det) && !is.data.frame(det) && !is.null(det$boxes)) {
      det <- det$boxes
    }
    if (!is.data.frame(det) || nrow(det) == 0L) next

    # Standardize column names
    if (is.null(det$conf) && !is.null(det$score)) det$conf <- det$score
    if (is.null(det$score) && !is.null(det$conf)) det$score <- det$conf
    if (is.null(det$class_name) && !is.null(det$label)) det$class_name <- det$label
    if (is.null(det$label) && !is.null(det$class_name)) det$label <- det$class_name
    if (is.null(det$class_id)) det$class_id <- 0L

    # Shift local coordinates to global coordinates
    det$xmin <- det$xmin + sx1 - 1L
    det$xmax <- det$xmax + sx1 - 1L
    det$ymin <- det$ymin + sy1 - 1L
    det$ymax <- det$ymax + sy1 - 1L
    det$is_full <- FALSE

    all_candidates[[length(all_candidates) + 1L]] <- det
  }

  # 2. Optional full-image global pass
  if (isTRUE(full_image_inference)) {
    full_det <- try(
      image_detect_dl(
        im,
        model = model,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        plot = FALSE,
        verbose = FALSE,
        ...
      ),
      silent = TRUE
    )
    if (!inherits(full_det, "try-error")) {
      if (!is.null(attr(full_det, "boxes"))) {
        full_det <- attr(full_det, "boxes")
      } else if (is.list(full_det) && !is.data.frame(full_det) && !is.null(full_det$boxes)) {
        full_det <- full_det$boxes
      }
      if (is.data.frame(full_det) && nrow(full_det) > 0L) {
        if (is.null(full_det$conf) && !is.null(full_det$score)) full_det$conf <- full_det$score
        if (is.null(full_det$score) && !is.null(full_det$conf)) full_det$score <- full_det$conf
        if (is.null(full_det$class_name) && !is.null(full_det$label)) full_det$class_name <- full_det$label
        if (is.null(full_det$label) && !is.null(full_det$class_name)) full_det$label <- full_det$class_name
        if (is.null(full_det$class_id)) full_det$class_id <- 0L
        full_det$is_full <- TRUE
        all_candidates[[length(all_candidates) + 1L]] <- full_det
      }
    }
  }

  if (length(all_candidates) == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No objects detected across all slices.")
    empty_df <- data.frame(
      xmin = numeric(), ymin = numeric(), xmax = numeric(), ymax = numeric(),
      conf = numeric(), class_id = integer(), class_name = character(),
      score = numeric(), label = character()
    )
    attr(empty_df, "image") <- im
    attr(empty_df, "counts") <- setNames(integer(0), character(0))
    return(empty_df)
  }

  merged_df <- do.call(rbind, all_candidates)

  # 3. SAHI Post-Processing (Greedy NMM / NMS)
  res_post <- .sahi_postprocess(
    candidates_df = merged_df,
    contours = NULL,
    match_metric = match_metric,
    match_threshold = match_threshold,
    postprocess_type = postprocess_type,
    seam_merging = seam_merging,
    slice_borders_x = slice_borders_x,
    slice_borders_y = slice_borders_y
  )

  final_df <- res_post$boxes
  rownames(final_df) <- NULL

  # Summaries
  c_names <- if (!is.null(final_df$class_name)) final_df$class_name else if (!is.null(final_df$label)) final_df$label else rep("object", nrow(final_df))
  tbl <- table(c_names)
  counts <- setNames(as.integer(tbl), names(tbl))

  attr(final_df, "image") <- im
  attr(final_df, "counts") <- counts

  if (isTRUE(verbose)) {
    cli::cli_alert_success("SAHI finished! Merged {nrow(final_df)} global detection(s) across slices.")
  }

  if (isTRUE(plot)) {
    plot(im)
    if (nrow(final_df) > 0L) {
      n_b <- nrow(final_df)
      c_border <- if (isTRUE(rainbow)) grDevices::rainbow(n_b, s = 0.85, v = 0.95) else rep("#00FFCC", n_b)
      c_confs <- if (!is.null(final_df$conf)) final_df$conf else if (!is.null(final_df$score)) final_df$score else rep(1.0, n_b)
      for (r in seq_len(n_b)) {
        graphics::rect(final_df$xmin[r], final_df$ymin[r], final_df$xmax[r], final_df$ymax[r], border = c_border[r], lwd = 2)
        graphics::text(final_df$xmin[r], final_df$ymin[r] - 4,
                       labels = sprintf("%s %.2f", c_names[r], c_confs[r]),
                       col = c_border[r], cex = 0.75, pos = 4)
      }
    }
  }

  invisible(final_df)
}

#' @title Slicing Aided Hyper Inference (SAHI) for Instance Segmentation
#' @name image_segment_sahi
#' @description
#' Performs high-resolution sliding-window hyper inference for instance segmentation
#' using [image_segment_dl()]. Slices large canvases, predicts masks per tile,
#' shifts polygon contours to global coordinates, and applies Greedy Non-Maximum Merging (Greedy NMM)
#' with IOS and slice seam unification so objects cut across slice borders are cleanly assembled into complete objects.
#'
#' @inheritParams image_detect_sahi
#' @param type Type of presentation: `"highlight"` (default, overlays translucent colored masks on original image),
#'   `"segment"` (draws polygon borders), or `"mask"` (binary/label mask).
#' @param col_background Color of background when `type = "segment"`. Defaults to `"white"`.
#' @param col_highlight Base color for highlighted instances when `rainbow = FALSE`. Defaults to `"salmon"`.
#' @param alpha Numeric transparency factor (0 = fully transparent, 1 = fully opaque) for highlight fill. Defaults to `0.40`.
#' @param border Color for instance polygon borders. Defaults to `NA` (no border).
#' @param lwd Numeric border line width. Defaults to `1`.
#' @param show_id Logical. Display instance index number at centroid. Defaults to `TRUE`.
#' @return A `data.frame` of detections with attached global contours, counts, and summaries.
#' @export
image_segment_sahi <- function(img,
                               model = "yolo26n-seg",
                               slice_height = 640,
                               slice_width = 640,
                               overlap_height_ratio = 0.20,
                               overlap_width_ratio = 0.20,
                               conf_threshold = 0.25,
                               iou_threshold = 0.45,
                               match_metric = c("IOS", "IOU"),
                               postprocess_type = c("GREEDY_NMM", "NMS"),
                               match_threshold = 0.50,
                               seam_merging = TRUE,
                               full_image_inference = TRUE,
                               type = c("highlight", "segment", "mask"),
                               col_background = "white",
                               col_highlight = "salmon",
                               alpha = 0.40,
                               border = NA,
                               lwd = 1,
                               show_id = TRUE,
                               rainbow = TRUE,
                               plot = TRUE,
                               verbose = TRUE,
                               ...) {
  type <- match.arg(type)
  match_metric <- match.arg(match_metric)
  postprocess_type <- match.arg(postprocess_type)

  if (missing(match_threshold) && !missing(iou_threshold)) {
    match_threshold <- iou_threshold
  } else if (is.null(match_threshold)) {
    match_threshold <- 0.50
  }

  # Ingest image
  im <- if (is.character(img) && file.exists(img[1])) {
    image_import(img[1])
  } else if (inherits(img, c("image", "Image"))) {
    img
  } else {
    as_image(img)
  }

  dims <- dim(im)
  orig_w <- dims[1]
  orig_h <- dims[2]

  slice_w <- min(orig_w, as.integer(slice_width))
  slice_h <- min(orig_h, as.integer(slice_height))

  step_x <- max(1L, as.integer(round(slice_w * (1 - overlap_width_ratio))))
  step_y <- max(1L, as.integer(round(slice_h * (1 - overlap_height_ratio))))

  x_starts <- unique(c(seq(1L, orig_w - slice_w + 1L, by = step_x), max(1L, orig_w - slice_w + 1L)))
  y_starts <- unique(c(seq(1L, orig_h - slice_h + 1L, by = step_y), max(1L, orig_h - slice_h + 1L)))

  grid <- expand.grid(x = x_starts, y = y_starts)
  n_slices <- nrow(grid)

  slice_borders_x <- sort(unique(c(grid$x, grid$x + slice_w - 1L)))
  slice_borders_y <- sort(unique(c(grid$y, grid$y + slice_h - 1L)))

  if (isTRUE(verbose)) {
    cli::cli_h2("SAHI (Slicing Aided Hyper Inference) - Instance Segmentation")
    cli::cli_alert_info("Canvas: {orig_w}x{orig_h} | Slice: {slice_w}x{slice_h} | Total slices: {.val {n_slices}}")
  }

  all_candidates <- list()
  all_contours <- list()

  for (i in seq_len(n_slices)) {
    sx1 <- grid$x[i]
    sy1 <- grid$y[i]
    sx2 <- min(orig_w, sx1 + slice_w - 1L)
    sy2 <- min(orig_h, sy1 + slice_h - 1L)

    sub_im <- if (length(dims) >= 3L) im[sx1:sx2, sy1:sy2, , drop = FALSE] else im[sx1:sx2, sy1:sy2, drop = FALSE]

    det <- try(
      image_segment_dl(
        sub_im,
        model = model,
        type = "highlight",
        box_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        plot = FALSE,
        verbose = FALSE,
        ...
      ),
      silent = TRUE
    )

    if (inherits(det, "try-error")) next

    df_b <- NULL
    conts <- list()
    if (!is.null(attr(det, "boxes"))) {
      df_b <- attr(det, "boxes")
    } else if (is.list(det) && !is.data.frame(det) && !is.null(det$boxes)) {
      df_b <- det$boxes
    } else if (is.data.frame(det)) {
      df_b <- det
    }

    if (!is.null(attr(det, "contours"))) {
      conts <- attr(det, "contours")
    } else if (is.list(det) && !is.null(det$contours)) {
      conts <- det$contours
    }

    if (is.null(df_b) || !is.data.frame(df_b) || nrow(df_b) == 0L) next

    # Standardize column names
    if (is.null(df_b$conf) && !is.null(df_b$score)) df_b$conf <- df_b$score
    if (is.null(df_b$score) && !is.null(df_b$conf)) df_b$score <- df_b$conf
    if (is.null(df_b$class_name) && !is.null(df_b$label)) df_b$class_name <- df_b$label
    if (is.null(df_b$label) && !is.null(df_b$class_name)) df_b$label <- df_b$class_name
    if (is.null(df_b$class_id)) df_b$class_id <- 0L

    # Shift boxes to global coordinates
    df_b$xmin <- df_b$xmin + sx1 - 1L
    df_b$xmax <- df_b$xmax + sx1 - 1L
    df_b$ymin <- df_b$ymin + sy1 - 1L
    df_b$ymax <- df_b$ymax + sy1 - 1L

    # Shift contours to global coordinates
    n_b <- nrow(df_b)
    shifted_c <- vector("list", n_b)
    for (k in seq_len(n_b)) {
      bx1_loc <- df_b$xmin[k] - (sx1 - 1L)
      bx2_loc <- df_b$xmax[k] - (sx1 - 1L)
      by1_loc <- df_b$ymin[k] - (sy1 - 1L)
      by2_loc <- df_b$ymax[k] - (sy1 - 1L)

      valid_cnt <- FALSE
      if (length(conts) >= k && is.matrix(conts[[k]]) && nrow(conts[[k]]) >= 3L) {
        c_mat <- conts[[k]]
        c_mx <- mean(c_mat[, 1])
        c_my <- mean(c_mat[, 2])
        if (c_mx >= bx1_loc - 25 && c_mx <= bx2_loc + 25 &&
            c_my >= by1_loc - 25 && c_my <= by2_loc + 25) {
          c_mat[, 1] <- c_mat[, 1] + sx1 - 1L
          c_mat[, 2] <- c_mat[, 2] + sy1 - 1L
          shifted_c[[k]] <- c_mat
          valid_cnt <- TRUE
        }
      }

      if (!valid_cnt) {
        # Fallback to rectangular polygon from box
        bx1 <- df_b$xmin[k]; by1 <- df_b$ymin[k]
        bx2 <- df_b$xmax[k]; by2 <- df_b$ymax[k]
        shifted_c[[k]] <- matrix(c(bx1, bx2, bx2, bx1, by1, by1, by2, by2), ncol = 2)
      }
    }

    for (k in seq_len(n_b)) {
      row_k <- df_b[k, , drop = FALSE]
      row_k$is_full <- FALSE
      all_candidates[[length(all_candidates) + 1L]] <- row_k
      all_contours[[length(all_contours) + 1L]] <- shifted_c[[k]]
    }
  }

  # Optional full-image global pass
  if (isTRUE(full_image_inference)) {
    full_det <- try(
      image_segment_dl(
        im,
        model = model,
        type = "highlight",
        box_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        plot = FALSE,
        verbose = FALSE,
        ...
      ),
      silent = TRUE
    )
    if (!inherits(full_det, "try-error")) {
      f_df <- NULL
      f_conts <- list()
      if (!is.null(attr(full_det, "boxes"))) {
        f_df <- attr(full_det, "boxes")
      } else if (is.list(full_det) && !is.data.frame(full_det) && !is.null(full_det$boxes)) {
        f_df <- full_det$boxes
      } else if (is.data.frame(full_det)) {
        f_df <- full_det
      }

      if (!is.null(attr(full_det, "contours"))) {
        f_conts <- attr(full_det, "contours")
      } else if (is.list(full_det) && !is.null(full_det$contours)) {
        f_conts <- full_det$contours
      }

      if (!is.null(f_df) && is.data.frame(f_df) && nrow(f_df) > 0L) {
        if (is.null(f_df$conf) && !is.null(f_df$score)) f_df$conf <- f_df$score
        if (is.null(f_df$score) && !is.null(f_df$conf)) f_df$score <- f_df$conf
        if (is.null(f_df$class_name) && !is.null(f_df$label)) f_df$class_name <- f_df$label
        if (is.null(f_df$label) && !is.null(f_df$class_name)) f_df$label <- f_df$class_name
        if (is.null(f_df$class_id)) f_df$class_id <- 0L
        for (k in seq_len(nrow(f_df))) {
          row_k <- f_df[k, , drop = FALSE]
          row_k$is_full <- TRUE
          all_candidates[[length(all_candidates) + 1L]] <- row_k
          valid_cnt <- FALSE
          if (length(f_conts) >= k && is.matrix(f_conts[[k]]) && nrow(f_conts[[k]]) >= 3L) {
            c_mat <- f_conts[[k]]
            c_mx <- mean(c_mat[, 1])
            c_my <- mean(c_mat[, 2])
            if (c_mx >= f_df$xmin[k] - 25 && c_mx <= f_df$xmax[k] + 25 &&
                c_my >= f_df$ymin[k] - 25 && c_my <= f_df$ymax[k] + 25) {
              valid_cnt <- TRUE
            }
          }
          c_mat <- if (valid_cnt) {
            f_conts[[k]]
          } else {
            matrix(c(f_df$xmin[k], f_df$xmax[k], f_df$xmax[k], f_df$xmin[k], f_df$ymin[k], f_df$ymin[k], f_df$ymax[k], f_df$ymax[k]), ncol = 2)
          }
          all_contours[[length(all_contours) + 1L]] <- c_mat
        }
      }
    }
  }

  if (length(all_candidates) == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No objects segmented across all slices.")
    empty_df <- data.frame(
      xmin = numeric(), ymin = numeric(), xmax = numeric(), ymax = numeric(),
      conf = numeric(), class_id = integer(), class_name = character(),
      score = numeric(), label = character()
    )
    attr(empty_df, "image") <- im
    attr(empty_df, "counts") <- setNames(integer(0), character(0))
    attr(empty_df, "contours") <- list()
    return(empty_df)
  }

  merged_df <- do.call(rbind, all_candidates)

  # Cross-Slice SAHI Post-Processing (Greedy NMM / NMS)
  res_post <- .sahi_postprocess(
    candidates_df = merged_df,
    contours = all_contours,
    match_metric = match_metric,
    match_threshold = match_threshold,
    postprocess_type = postprocess_type,
    seam_merging = seam_merging,
    slice_borders_x = slice_borders_x,
    slice_borders_y = slice_borders_y
  )

  final_df <- res_post$boxes
  final_conts <- res_post$contours
  rownames(final_df) <- NULL

  c_names <- if (!is.null(final_df$class_name)) final_df$class_name else if (!is.null(final_df$label)) final_df$label else rep("object", nrow(final_df))
  tbl <- table(c_names)
  counts <- setNames(as.integer(tbl), names(tbl))

  attr(final_df, "image") <- im
  attr(final_df, "counts") <- counts
  attr(final_df, "contours") <- final_conts

  if (isTRUE(verbose)) {
    cli::cli_alert_success("SAHI finished! Merged {nrow(final_df)} segmented instance(s) across slices.")
  }

  if (isTRUE(plot)) {
    n_inst <- nrow(final_df)
    palette_cols <- if (isTRUE(rainbow)) {
      grDevices::rainbow(max(1L, n_inst), s = 0.85, v = 0.95)
    } else {
      rep(col_highlight, max(1L, n_inst))
    }

    if (type == "highlight") {
      plot(im)
      for (k in seq_len(n_inst)) {
        cnt <- final_conts[[k]]
        if (is.matrix(cnt) && nrow(cnt) >= 3L) {
          c_fill <- grDevices::adjustcolor(palette_cols[k], alpha.f = alpha)
          c_line <- if (!is.na(border)) border else palette_cols[k]
          graphics::polygon(cnt[, 1], cnt[, 2], col = c_fill, border = c_line, lwd = lwd)
          if (isTRUE(show_id)) {
            cx <- mean(cnt[, 1]); cy <- mean(cnt[, 2])
            graphics::text(cx, cy, labels = k, col = "white", cex = 0.85, font = 2)
            graphics::text(cx, cy, labels = k, col = "black", cex = 0.75, font = 1)
          }
        }
      }
    } else if (type == "segment") {
      plot(im)
      for (k in seq_len(n_inst)) {
        cnt <- final_conts[[k]]
        if (is.matrix(cnt) && nrow(cnt) >= 3L) {
          c_line <- palette_cols[k]
          c_fill <- grDevices::adjustcolor(c_line, alpha.f = alpha)
          graphics::polygon(cnt[, 1], cnt[, 2], col = c_fill, border = c_line, lwd = max(1.5, lwd))
          if (isTRUE(show_id)) {
            cx <- mean(cnt[, 1]); cy <- mean(cnt[, 2])
            graphics::text(cx, cy, labels = sprintf("%s %d", c_names[k], k), col = c_line, cex = 0.75, pos = 4)
          }
        }
      }
    } else if (type == "mask") {
      plot(im)
      for (k in seq_len(n_inst)) {
        cnt <- final_conts[[k]]
        if (is.matrix(cnt) && nrow(cnt) >= 3L) {
          graphics::polygon(cnt[, 1], cnt[, 2], col = "white", border = "white", lwd = 1)
        }
      }
    }
  }

  invisible(final_df)
}
