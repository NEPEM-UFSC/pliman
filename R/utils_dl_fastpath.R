# ==============================================================================
# DEEP LEARNING FAST PATHS
# pliman: Plant Image Analysis in R
#
# Internal helpers that keep per-image / per-tile / per-frame inference work in C++ or
# restricted to small windows, shared by image, SAHI and video pipelines.
# ==============================================================================

# Computes one SAM instance mask inside a window around its box.
#
# `prob` is the low-resolution probability map returned by the SAM decoder for the whole
# image. The window (box + 25% margin, at least 16 px, + morphology radius) is resampled
# with `resize_roi_bilinear_cpp()`, which yields exactly the values of the full-resolution
# bilinear resize inside the window, so the mask is the same as the former full-image
# processing except for SAM pixels far outside the prompt box. Cost is proportional to
# the window area instead of the full image.
# Returns `list(x = rx1:rx2, y = ry1:ry2, m = <logical window mask>)` (or NULL) so the
# caller can paint `labels[x, y][m] <- k` in place without copying the full label matrix.
.instance_window_mask <- function(prob, box, orig_w, orig_h, threshold,
                                  fill_hull = TRUE, filter = 0, erode = 0, dilate = 0,
                                  opening = 0, closing = 0, min_area = 0) {
  morph_pad <- 2L * as.integer(max(0, filter, erode, dilate, opening, closing)) + 2L
  bw_k <- max(1, box[3] - box[1])
  bh_k <- max(1, box[4] - box[2])
  px <- as.integer(ceiling(max(16, 0.25 * bw_k))) + morph_pad
  py <- as.integer(ceiling(max(16, 0.25 * bh_k))) + morph_pad
  rx1 <- as.integer(max(1L, floor(box[1]) - px))
  rx2 <- as.integer(min(orig_w, ceiling(box[3]) + px))
  ry1 <- as.integer(max(1L, floor(box[2]) - py))
  ry2 <- as.integer(min(orig_h, ceiling(box[4]) + py))
  if (rx2 < rx1 || ry2 < ry1) return(NULL)
  rw <- rx2 - rx1 + 1L
  rh <- ry2 - ry1 + 1L
  if (!is.matrix(prob)) prob <- as.matrix(prob)
  storage.mode(prob) <- "double"
  mask_prob <- resize_roi_bilinear_cpp(prob, as.integer(orig_w), as.integer(orig_h),
                                       rx1 - 1L, ry1 - 1L, rw, rh)
  m <- mask_prob >= threshold

  # box +/- 2 px, in window coordinates
  x1 <- max(1L, floor(box[1]) - 2L) - rx1 + 1L
  x2 <- min(orig_w, ceiling(box[3]) + 2L) - rx1 + 1L
  y1 <- max(1L, floor(box[2]) - 2L) - ry1 + 1L
  y2 <- min(orig_h, ceiling(box[4]) + 2L) - ry1 + 1L

  if (isTRUE(fill_hull) && x2 >= x1 && y2 >= y1) {
    sub_m <- m[x1:x2, y1:y2, drop = FALSE]
    m[x1:x2, y1:y2] <- (fill_holes_cpp(sub_m) > 0)
  }
  if (filter > 0) {
    m_arr <- as.array(m)
    dim(m_arr) <- c(rw, rh, 1L)
    f_res <- median_filter_binary_cpp(m_arr, rw, rh, 1L, as.integer(filter))
    m <- (f_res[, , 1] > 0)
  }
  if (erode > 0) m <- (erode_cpp(m, raio = as.integer(erode)) > 0)
  if (dilate > 0) m <- (dilate_cpp(m, raio = as.integer(dilate)) > 0)
  if (opening > 0) {
    m <- (erode_cpp(m, raio = as.integer(opening)) > 0)
    m <- (dilate_cpp(m, raio = as.integer(opening)) > 0)
  }
  if (closing > 0) {
    m <- (dilate_cpp(m, raio = as.integer(closing)) > 0)
    m <- (erode_cpp(m, raio = as.integer(closing)) > 0)
  }
  if (min_area > 0) {
    lbls <- bwlabel_cpp(m)
    if (max(lbls) > 0) {
      tb <- tabulate(lbls[lbls > 0])
      keep_ids <- which(tb >= min_area)
      m <- matrix(lbls %in% keep_ids, nrow = rw, ncol = rh)
    }
  }

  list(x = rx1:rx2, y = ry1:ry2, m = m)
}

# Resolves everything `run_yolo_cpp()` needs once per call (model file, ONNX Runtime
# library, class names, inference size), mirroring the logic of `.run_yolo()`. Used by
# the SAHI fast path so per-tile inference does not repeat model lookups.
.resolve_yolo_runtime <- function(model, dir = pliman_model_dir(), labels = NULL, imgsz = NULL,
                                  threads = 0, engine = c("gpu", "cpu"), device_id = -1) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path(engine = engine)
  if (is.null(lib_path) || !file.exists(lib_path)) return(NULL)
  model_file <- if (file.exists(model)) {
    normalizePath(model, winslash = "/")
  } else {
    hit <- .find_local_model(model, dir = dir)
    if (!is.null(hit)) hit else pliman_download_model(model = model, dir = dir)
  }
  shp <- tryCatch(
    get_onnx_model_io_shape_cpp(model_file, lib_path, as.integer(threads), use_gpu, as.integer(device_id)),
    error = function(e) list(dynamic = TRUE, expected_w = 640L, expected_h = 640L)
  )
  target_size <- if (!is.null(imgsz)) {
    as.integer(imgsz[1])
  } else {
    m_sz <- .get_onnx_imgsz(model_file)
    if (!is.null(m_sz) && length(m_sz) > 0L) as.integer(m_sz[1]) else as.integer(shp$expected_w)
  }
  if (is.na(target_size) || target_size <= 0L) target_size <- 640L
  if (!isTRUE(shp$dynamic) && target_size != shp$expected_w) target_size <- as.integer(shp$expected_w)
  class_names <- if (!is.null(labels)) {
    labels
  } else {
    cl <- .get_onnx_classes(model_file)
    if (!is.null(cl) && length(cl) > 0L) cl else .coco_classes
  }
  list(model_file = model_file, lib_path = lib_path, use_gpu = use_gpu, threads = as.integer(threads),
       device_id = as.integer(device_id), target_size = as.integer(target_size),
       class_names = class_names, batch_dynamic = isTRUE(shp$batch_dynamic))
}

# TRUE when SAHI tiles can bypass `image_detect_dl()` and run as one C++ call.
.sahi_fast_eligible <- function(model, dots) {
  if (!isTRUE(getOption("pliman.sahi_fast", TRUE))) return(FALSE)
  if (!is.character(model) || length(model) != 1L) return(FALSE)
  if (grepl("world|-cls|_cls|-sem|_sem|depth|obb", tolower(basename(model)))) return(FALSE)
  if (isTRUE(dots$return_features)) return(FALSE)
  if (!is.null(dots$background) || isTRUE(dots$filter_nested)) return(FALSE)
  TRUE
}

# Runs a SAHI tile grid (1-based tile origins from `grid`) plus, optionally, the full
# image through `run_yolo_tiles_cpp()` and returns the candidate data.frame expected by
# `.sahi_postprocess()` (global pixel coordinates, `is_full` flag).
.sahi_tiles_candidates <- function(raw_im, grid, slice_w, slice_h, rt, conf_threshold,
                                   iou_threshold, full_image_inference = TRUE) {
  if (!is.raw(raw_im) && !is.double(raw_im)) storage.mode(raw_im) <- "double"
  res <- run_yolo_tiles_cpp(
    img = raw_im,
    tile_x = as.integer(grid$x - 1L),
    tile_y = as.integer(grid$y - 1L),
    tile_w = as.integer(slice_w),
    tile_h = as.integer(slice_h),
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    model_path = rt$model_file,
    lib_path = rt$lib_path,
    num_threads = rt$threads,
    use_gpu = rt$use_gpu,
    device_id = rt$device_id,
    target_size = rt$target_size,
    include_full = isTRUE(full_image_inference)
  )
  n <- nrow(res$boxes)
  if (n == 0L) return(NULL)
  # Same candidate order as the per-slice path: slices first (in grid order), then the
  # full-image pass; ids restart at 1 within each slice, as image_detect_dl() returns them.
  tix <- as.integer(res$tile_index)
  ord <- order(tix == 0L, tix, seq_len(n))
  tix <- tix[ord]
  bx <- res$boxes[ord, , drop = FALSE]
  c_ids <- as.integer(res$class_ids)[ord]
  cn <- rt$class_names
  lbls <- ifelse(c_ids >= 0 & c_ids < length(cn), cn[c_ids + 1L], as.character(c_ids))
  sc <- round(res$scores[ord], 4)
  data.frame(
    id = stats::ave(seq_len(n), tix, FUN = seq_along),
    xmin = bx[, 1], ymin = bx[, 2],
    xmax = bx[, 3], ymax = bx[, 4],
    label = lbls, score = sc, class_id = c_ids,
    cx = round((bx[, 1] + bx[, 3]) / 2, 1),
    cy = round((bx[, 2] + bx[, 4]) / 2, 1),
    conf = sc, class_name = lbls,
    is_full = tix == 0L,
    stringsAsFactors = FALSE
  )
}

# Column-wise row-binding of many small data.frames with identical columns (e.g. one
# detection table per video frame). `do.call(rbind, lst)` on thousands of data.frames is
# quadratic-ish in R; this concatenates each column once. Falls back to rbind() when the
# tables do not share the same column layout or contain classed columns.
.rbind_df_fast <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (length(lst) == 0L) return(data.frame())
  lst <- lst[vapply(lst, function(d) is.data.frame(d) && nrow(d) > 0L, logical(1))]
  if (length(lst) == 0L) return(data.frame())
  if (length(lst) == 1L) return(lst[[1L]])
  cols <- names(lst[[1L]])
  same <- all(vapply(lst, function(d) identical(names(d), cols), logical(1)))
  # classed columns (factor, Date, POSIXct, ...) would lose their class with unlist()
  has_class <- any(vapply(lst[[1L]], function(x) !is.null(attr(x, "class")), logical(1)))
  if (!same || has_class) return(do.call(rbind, lst))
  out <- lapply(cols, function(cn) {
    v <- lapply(lst, `[[`, cn)
    if (any(vapply(v, is.list, logical(1)))) return(do.call(c, v))
    unlist(v, use.names = FALSE)
  })
  names(out) <- cols
  out <- as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
  names(out) <- cols
  out
}

# Builds a minimal enclosure for batch worker functions so that serialising the worker
# to mirai daemons does not drag the caller's frame (the full list of input images,
# intermediate results, ...) along with every chunk.
.slim_worker_env <- function(...) {
  list2env(list(...), parent = asNamespace("pliman"))
}

# Tunes ONNX settings for parallel batch execution:
# * CPU engine with `threads = 0` (ORT default = all cores) would oversubscribe the CPU
#   N-fold with N workers, so each worker gets `physical_cores %/% workers` threads.
# * GPU engine: every worker opens its own GPU session (VRAM x N, contention on a single
#   device), so concurrency is capped at `getOption("pliman.gpu_workers", 2L)`.
.tune_parallel_dl <- function(args, workers, verbose = TRUE) {
  eng <- if (!is.null(args$engine)) as.character(args$engine[1]) else "gpu"
  gpu_mode <- identical(eng, "gpu") && identical(tolower(Sys.info()[["sysname"]]), "windows")
  if (gpu_mode) {
    cap <- max(1L, as.integer(getOption("pliman.gpu_workers", 2L)))
    if (workers > cap) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning(c(
          "GPU engine: limiting concurrent workers to {.val {cap}} (each worker loads its own GPU session).",
          "i" = "Change with {.code options(pliman.gpu_workers = n)} or use {.code engine = \"cpu\"} for many workers."
        ))
      }
      workers <- cap
    }
  } else if (is.null(args$threads) || isTRUE(all(args$threads == 0))) {
    cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    if (is.na(cores) || cores < 1L) cores <- max(1L, parallel::detectCores())
    args$threads <- max(1L, as.integer(cores) %/% max(1L, as.integer(workers)))
  }
  list(args = args, workers = as.integer(workers), gpu = gpu_mode)
}
