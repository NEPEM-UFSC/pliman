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
                              agnostic_nms = FALSE,
                              ...) {
  match_metric <- match.arg(match_metric)
  postprocess_type <- match.arg(postprocess_type)

  if (is.null(candidates_df) || nrow(candidates_df) == 0L) {
    return(list(boxes = candidates_df, contours = contours))
  }

  has_conts <- !is.null(contours) && length(contours) == nrow(candidates_df)

  # Standardize score and confidence
  c_vals <- if (!is.null(candidates_df$conf)) {
    as.numeric(candidates_df$conf)
  } else if (!is.null(candidates_df$score)) {
    as.numeric(candidates_df$score)
  } else {
    rep(1.0, nrow(candidates_df))
  }
  candidates_df$conf <- c_vals
  candidates_df$score <- c_vals

  if (is.null(candidates_df$is_full)) {
    candidates_df$is_full <- FALSE
  }

  class_groups <- if (isTRUE(agnostic_nms)) {
    list(all = seq_len(nrow(candidates_df)))
  } else {
    split(seq_len(nrow(candidates_df)), candidates_df$class_id)
  }

  final_boxes_list <- list()
  final_conts_list <- list()
  is_nmm <- (postprocess_type == "GREEDY_NMM")
  is_ios <- (match_metric == "IOS")

  for (cls_idx in class_groups) {
    cls_df <- candidates_df[cls_idx, , drop = FALSE]
    cls_conts <- if (has_conts) contours[cls_idx] else NULL

    x1 <- as.numeric(cls_df$xmin)
    y1 <- as.numeric(cls_df$ymin)
    x2 <- as.numeric(cls_df$xmax)
    y2 <- as.numeric(cls_df$ymax)
    cf <- as.numeric(cls_df$conf)
    isf <- as.logical(cls_df$is_full)

    # ------------------------------------------------------------------
    # Step 1: Vectorized subsumption of slice fragments into full-image detections
    # ------------------------------------------------------------------
    if (any(isf)) {
      full_idx <- which(isf)
      slice_idx <- which(!isf)
      if (length(slice_idx) > 0L) {
        sx1 <- x1[slice_idx]; sy1 <- y1[slice_idx]
        sx2 <- x2[slice_idx]; sy2 <- y2[slice_idx]
        scf <- cf[slice_idx]
        s_area <- (sx2 - sx1 + 1.0) * (sy2 - sy1 + 1.0)
        absorbed <- logical(length(slice_idx))

        for (fi in full_idx) {
          rem_s <- which(!absorbed)
          if (length(rem_s) == 0L) break
          fx1 <- x1[fi]; fy1 <- y1[fi]; fx2 <- x2[fi]; fy2 <- y2[fi]
          rem_s <- rem_s[sx2[rem_s] >= fx1 & sx1[rem_s] <= fx2 & sy2[rem_s] >= fy1 & sy1[rem_s] <= fy2]
          if (length(rem_s) == 0L) next
          f_area <- (fx2 - fx1 + 1.0) * (fy2 - fy1 + 1.0)
          iw <- pmin(fx2, sx2[rem_s]) - pmax(fx1, sx1[rem_s]) + 1.0
          ih <- pmin(fy2, sy2[rem_s]) - pmax(fy1, sy1[rem_s]) + 1.0
          inter <- iw * ih
          sa <- s_area[rem_s]
          ios <- inter / pmin(f_area, sa)
          iou <- inter / (f_area + sa - inter)
          hit <- rem_s[ios >= 0.40 | iou >= 0.35]
          if (length(hit) > 0L) {
            max_s <- max(scf[hit])
            if (max_s > cf[fi]) cf[fi] <- max_s
            absorbed[hit] <- TRUE
          }
        }
        cls_df$conf <- cf
        cls_df$score <- cf
        if (any(absorbed)) {
          keep_rows <- rep(TRUE, nrow(cls_df))
          keep_rows[slice_idx[absorbed]] <- FALSE
          cls_df <- cls_df[keep_rows, , drop = FALSE]
          if (has_conts) cls_conts <- cls_conts[keep_rows]
          x1 <- x1[keep_rows]; y1 <- y1[keep_rows]
          x2 <- x2[keep_rows]; y2 <- y2[keep_rows]
          cf <- cf[keep_rows]; isf <- isf[keep_rows]
        }
      }
    }

    n_rem <- nrow(cls_df)
    if (n_rem == 0L) next

    # ------------------------------------------------------------------
    # Step 2: Vectorized Greedy NMM / NMS with IOS (cluster-based)
    # ------------------------------------------------------------------
    ord <- order(cf, decreasing = TRUE)
    cls_df <- cls_df[ord, , drop = FALSE]
    if (has_conts) cls_conts <- cls_conts[ord]
    x1 <- x1[ord]; y1 <- y1[ord]; x2 <- x2[ord]; y2 <- y2[ord]
    cf <- cf[ord]; isf <- isf[ord]
    areas <- (x2 - x1 + 1.0) * (y2 - y1 + 1.0)

    alive <- rep(TRUE, n_rem)
    keep_seeds <- integer(n_rem)
    out_x1 <- numeric(n_rem); out_y1 <- numeric(n_rem)
    out_x2 <- numeric(n_rem); out_y2 <- numeric(n_rem)
    out_cf <- numeric(n_rem)
    n_out <- 0L

    for (seed in seq_len(n_rem)) {
      if (!alive[seed]) next
      alive[seed] <- FALSE
      n_out <- n_out + 1L
      keep_seeds[n_out] <- seed

      bx1 <- x1[seed]; by1 <- y1[seed]; bx2 <- x2[seed]; by2 <- y2[seed]
      b_area <- areas[seed]

      cand <- which(alive)
      if (length(cand) > 0L) {
        cand <- cand[x2[cand] >= bx1 & x1[cand] <= bx2]
        if (length(cand) > 0L) {
          cand <- cand[y2[cand] >= by1 & y1[cand] <= by2]
        }
      }

      if (length(cand) > 0L) {
        iw <- pmin(bx2, x2[cand]) - pmax(bx1, x1[cand]) + 1.0
        ih <- pmin(by2, y2[cand]) - pmax(by1, y1[cand]) + 1.0
        inter <- iw * ih
        ca <- areas[cand]
        iou <- inter / (b_area + ca - inter)
        matched <- if (is_ios) {
          ios <- inter / pmin(b_area, ca)
          cand[ios >= match_threshold & iou >= 0.30]
        } else {
          cand[iou >= match_threshold]
        }
      } else {
        matched <- integer(0)
      }

      if (length(matched) > 0L) {
        alive[matched] <- FALSE
        if (is_nmm) {
          m_idx <- c(seed, matched)
          out_x1[n_out] <- min(x1[m_idx])
          out_y1[n_out] <- min(y1[m_idx])
          out_x2[n_out] <- max(x2[m_idx])
          out_y2[n_out] <- max(y2[m_idx])
          out_cf[n_out] <- max(cf[m_idx])
          if (has_conts) {
            m_cnts <- cls_conts[m_idx]
            merged_cnt <- NULL
            if (isTRUE(isf[seed]) && is.matrix(m_cnts[[1]]) && nrow(m_cnts[[1]]) >= 3L) {
              merged_cnt <- m_cnts[[1]]
            } else {
              all_pts <- do.call(rbind, lapply(m_cnts, function(p) if (is.matrix(p) && nrow(p) >= 3L) p else NULL))
              if (!is.null(all_pts) && nrow(all_pts) >= 3L) {
                ch <- chull(all_pts)
                merged_cnt <- all_pts[c(ch, ch[1]), , drop = FALSE]
              } else if (length(m_cnts) > 0L && is.matrix(m_cnts[[1]])) {
                merged_cnt <- m_cnts[[1]]
              }
            }
            final_conts_list[[length(final_conts_list) + 1L]] <- merged_cnt
          }
        } else {
          out_x1[n_out] <- bx1; out_y1[n_out] <- by1
          out_x2[n_out] <- bx2; out_y2[n_out] <- by2
          out_cf[n_out] <- cf[seed]
          if (has_conts) final_conts_list[[length(final_conts_list) + 1L]] <- cls_conts[[seed]]
        }
      } else {
        out_x1[n_out] <- bx1; out_y1[n_out] <- by1
        out_x2[n_out] <- bx2; out_y2[n_out] <- by2
        out_cf[n_out] <- cf[seed]
        if (has_conts) final_conts_list[[length(final_conts_list) + 1L]] <- cls_conts[[seed]]
      }
    }

    if (n_out > 0L) {
      idx_k <- keep_seeds[seq_len(n_out)]
      res_cls <- cls_df[idx_k, , drop = FALSE]
      if (is_nmm) {
        res_cls$xmin <- out_x1[seq_len(n_out)]
        res_cls$ymin <- out_y1[seq_len(n_out)]
        res_cls$xmax <- out_x2[seq_len(n_out)]
        res_cls$ymax <- out_y2[seq_len(n_out)]
        res_cls$conf <- out_cf[seq_len(n_out)]
        if (!is.null(res_cls$score)) res_cls$score <- res_cls$conf
      }
      final_boxes_list[[length(final_boxes_list) + 1L]] <- res_cls
    }
  }

  if (length(final_boxes_list) == 0L) {
    empty_df <- candidates_df[0, , drop = FALSE]
    empty_df$is_full <- NULL
    return(list(boxes = empty_df, contours = list()))
  }

  out_df <- if (length(final_boxes_list) == 1L) final_boxes_list[[1L]] else do.call(rbind, final_boxes_list)
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
#' @param classes Optional character or integer vector of class names or IDs to keep (e.g., `"car"` or `c("car", "truck")`).
#'   Detections belonging to other classes are discarded before post-processing. Defaults to `NULL` (all classes kept).
#' @param agnostic_nms Logical. If `TRUE`, bounding boxes and masks from different classes compete during post-processing suppression. Defaults to `FALSE`.
#' @param seam_merging Logical. If `TRUE` (default), seamlessly detects and merges objects that have been
#'   sliced in half across tile cut boundaries.
#' @param full_image_inference Logical. If `TRUE` (default), runs an additional detection pass on the
#'   resized full canvas to capture objects larger than individual slices.
#' @param cex Numeric character expansion for labels and badges. Defaults to `1.0`.
#' @param bbox Logical. If `TRUE` (default), draws bounding boxes.
#' @param show_box Logical. Explicitly control whether bounding boxes are drawn. Defaults to `TRUE` (or value of `bbox`). Supported aliases: `box`, `boxes`, `show_boxes`.
#' @param show_point Logical. If `TRUE`, draws points at the centroid of each detected object. Defaults to `FALSE`. Supported aliases: `points`, `point`, `show_points`.
#' @param point_size Numeric. Size (cex) of centroid points when `show_point = TRUE`. Defaults to `1.5`. Supported aliases: `marker_size`, `point_cex`.
#' @param point_pch Integer. Symbol code for centroid points. Defaults to `16` (filled circle). Supported alias: `marker_pch`.
#' @param point_col Color for centroid points. Defaults to `NULL` (uses matching box/instance color). Supported alias: `marker_col`.
#' @param show_text Logical. Display text labels on bounding boxes. Defaults to `TRUE`.
#' @param show_conf Logical. Display confidence score on labels/badges. Defaults to `TRUE`.
#' @param show_class Logical. Display class name on labels/badges. Defaults to `TRUE`.
#' @param show_id Logical. Display instance ID number. Defaults to `FALSE`.
#' @param badge Logical. Display rounded badge background behind text labels. Defaults to `TRUE`.
#' @param rainbow Logical. If `TRUE` (default), bounding boxes are rendered in rainbow colors.
#' @param lower_size Optional numeric value to filter out small detected objects.
#'   When `size_mode = "relative"` (default), objects with a bounding box area
#'   less than `lower_size` times the mean area of all candidate boxes are removed
#'   (e.g., `lower_size = 0.2` removes objects smaller than 20% of the average object size).
#'   When `size_mode = "absolute"`, objects with an area (in pixels) below `lower_size` are removed.
#'   Defaults to `NULL` (no lower size filter).
#' @param upper_size Optional numeric value to filter out abnormally large detected objects.
#'   When `size_mode = "relative"` (default), objects with a bounding box area
#'   greater than `upper_size` times the mean area of all candidate boxes are removed
#'   (e.g., `upper_size = 2.0` removes objects larger than 200% of the average object size).
#'   When `size_mode = "absolute"`, objects with an area (in pixels) above `upper_size` are removed.
#'   Defaults to `NULL` (no upper size filter).
#' @param size_mode Character string specifying how `lower_size` and `upper_size` are interpreted.
#'   `"relative"` (default) treats values as proportions/multiples of the mean bounding box area.
#'   `"absolute"` treats values as area in pixels.
#' @param pad Numeric padding around bounding box labels or badges. Defaults to `1.0`.
#' @param plot Logical. Display global detections overlaid on the full image. Defaults to `TRUE`.
#' @param verbose Logical. Show progress bar and tile metrics. Defaults to `TRUE`.
#' @param imgsz Target input image resolution for YOLO inference (e.g., \code{640}, \code{960}, or \code{1024}). Default is \code{NULL}.
#' @param ... Additional arguments passed to [image_detect_dl()].
#'
#' @return A `data.frame` of global bounding boxes with attributes `image`, `counts`, and `summary`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Run SAHI on a 4000x3000 drone orthophoto
#' res <- image_detect_sahi("drone_field.jpg", model = "seedlings.onnx", slice_width = 640)
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
                              classes = NULL,
                              agnostic_nms = FALSE,
                              seam_merging = TRUE,
                              full_image_inference = TRUE,
                              lower_size = NULL,
                              upper_size = NULL,
                              size_mode = c("relative", "absolute"),
                              cex = 1.0,
                              pad = 1.0,
                              bbox = TRUE,
                              show_box = NULL,
                              show_point = FALSE,
                              point_size = 1.5,
                              point_pch = 16,
                              point_col = NULL,
                              show_text = TRUE,
                              show_conf = TRUE,
                              show_class = TRUE,
                              show_id = FALSE,
                              badge = TRUE,
                              rainbow = TRUE,
                              plot = TRUE,
                              verbose = TRUE,
                              imgsz = NULL,
                              ...) {
  dots <- list(...)
  if (is.null(imgsz)) {
    if ("imgz" %in% names(dots)) imgsz <- dots$imgz
    if ("img_size" %in% names(dots)) imgsz <- dots$img_size
    if ("image_size" %in% names(dots)) imgsz <- dots$image_size
  }
  if ("lower_size" %in% names(dots)) lower_size <- dots$lower_size
  if ("upper_size" %in% names(dots)) upper_size <- dots$upper_size
  if ("size_mode" %in% names(dots)) size_mode <- dots$size_mode
  if ("show_points" %in% names(dots)) show_point <- isTRUE(dots$show_points)
  if ("show_point" %in% names(dots)) show_point <- isTRUE(dots$show_point)
  if ("points" %in% names(dots)) show_point <- isTRUE(dots$points)
  if ("point" %in% names(dots)) show_point <- isTRUE(dots$point)
  if ("show_boxes" %in% names(dots)) show_box <- isTRUE(dots$show_boxes)
  if ("box" %in% names(dots)) show_box <- isTRUE(dots$box)
  if ("boxes" %in% names(dots)) show_box <- isTRUE(dots$boxes)
  if (is.null(show_box)) {
    show_box <- isTRUE(bbox)
  }
  if ("marker_size" %in% names(dots)) point_size <- dots$marker_size
  if ("point_cex" %in% names(dots)) point_size <- dots$point_cex
  if ("marker_col" %in% names(dots)) point_col <- dots$marker_col
  if ("marker_pch" %in% names(dots)) point_pch <- dots$marker_pch
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)
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
    cli::cli_progress_step("SAHI (Slicing Aided Hyper Inference) - Object Detection")
    cli::cli_alert_info("Canvas: {orig_w}x{orig_h} | Slice: {slice_w}x{slice_h} | Total slices: {.val {n_slices}}")
  }

  all_candidates <- list()
  raw_im <- if (isS4(im) && methods::.hasSlot(im, ".Data")) methods::slot(im, ".Data") else im

  # 1. Sliding window slice inference
  for (i in seq_len(n_slices)) {
    sx1 <- grid$x[i]
    sy1 <- grid$y[i]
    sx2 <- min(orig_w, sx1 + slice_w - 1L)
    sy2 <- min(orig_h, sy1 + slice_h - 1L)

    sub_im <- if (length(dims) >= 3L) raw_im[sx1:sx2, sy1:sy2, , drop = FALSE] else raw_im[sx1:sx2, sy1:sy2, drop = FALSE]

    det <- try(
      image_detect_dl(
        sub_im,
        model = model,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        classes = classes,
        plot = FALSE,
        verbose = FALSE,
        imgsz = imgsz,
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
        classes = classes,
        plot = FALSE,
        verbose = FALSE,
        save_image = FALSE,
        imgsz = imgsz,
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

  # Filter by target classes if specified
  if (!is.null(classes) && nrow(merged_df) > 0L) {
    c_match <- tolower(as.character(classes))
    idx_keep <- tolower(merged_df$class_name) %in% c_match |
                tolower(merged_df$label) %in% c_match
    if (is.numeric(classes)) {
      idx_keep <- idx_keep | merged_df$class_id %in% classes
    }
    merged_df <- merged_df[idx_keep, , drop = FALSE]
  }

  if (nrow(merged_df) == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No objects matched the requested class(es).")
    empty_df <- data.frame(
      xmin = numeric(), ymin = numeric(), xmax = numeric(), ymax = numeric(),
      conf = numeric(), class_id = integer(), class_name = character(),
      score = numeric(), label = character()
    )
    attr(empty_df, "image") <- im
    attr(empty_df, "counts") <- setNames(integer(0), character(0))
    return(empty_df)
  }

  # 3. SAHI Post-Processing (Greedy NMM / NMS)
  if (isTRUE(verbose)) cli::cli_progress_step("SAHI Post-Processing")
  res_post <- .sahi_postprocess(
    candidates_df = merged_df,
    contours = NULL,
    match_metric = match_metric,
    match_threshold = match_threshold,
    postprocess_type = postprocess_type,
    agnostic_nms = agnostic_nms,
    seam_merging = seam_merging,
    slice_borders_x = slice_borders_x,
    slice_borders_y = slice_borders_y
  )

  final_df <- res_post$boxes
  rownames(final_df) <- NULL

  if ((!is.null(lower_size) || !is.null(upper_size)) && nrow(final_df) > 0L) {
    size_mode <- match.arg(size_mode)
    f_res <- .filter_by_size(
      boxes = final_df,
      lower_size = lower_size,
      upper_size = upper_size,
      size_mode = size_mode
    )
    final_df <- f_res$boxes
  }

  # Summaries
  c_names <- if (!is.null(final_df$class_name)) final_df$class_name else if (!is.null(final_df$label)) final_df$label else rep("object", nrow(final_df))
  tbl <- table(c_names)
  counts <- setNames(as.integer(tbl), names(tbl))

  if (nrow(final_df) > 0L) {
    if (is.null(final_df$id)) final_df$id <- seq_len(nrow(final_df))
    if (is.null(final_df$label)) final_df$label <- c_names
    if (is.null(final_df$score)) final_df$score <- if (!is.null(final_df$conf)) final_df$conf else 1.0
  }

  attr(final_df, "image") <- im
  attr(final_df, "counts") <- counts

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
    cli::cli_alert_success("SAHI finished! Merged {nrow(final_df)} global detection(s) across slices.")
  }

  if (isTRUE(plot)) {
    plot(im)
    if (nrow(final_df) > 0L) {
      n_b <- nrow(final_df)
      palette_cols <- if (isTRUE(rainbow)) grDevices::rainbow(n_b, s = 0.85, v = 0.95) else rep("#00FFCC", n_b)

      .plot_yolo_bboxes(
        boxes = final_df,
        palette_colors = palette_cols,
        lwd = 2,
        cex = cex,
        pad = pad,
        show_text = show_text,
        show_conf = show_conf,
        show_class = show_class,
        show_id = show_id,
        show_box = show_box,
        show_point = show_point,
        point_size = point_size,
        point_pch = point_pch,
        point_col = point_col,
        badge = badge
      )
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
#' @param classes Optional character or integer vector of class names or IDs to keep (e.g., `"car"` or `c("car", "truck")`).
#'   Detections belonging to other classes are discarded before post-processing. Defaults to `NULL` (all classes kept).
#' @param agnostic_nms Logical. If `TRUE`, bounding boxes and masks from different classes compete during post-processing suppression. Defaults to `FALSE`.
#' @param bbox Logical. If `TRUE`, draws bounding boxes and class badges around segmented instances. Defaults to `FALSE`.
#' @param show_box Logical. Explicitly control whether bounding boxes are drawn. Defaults to `NULL` (or value of `bbox`). Supported aliases: `box`, `boxes`, `show_boxes`.
#' @param show_point Logical. If `TRUE`, draws points at the centroid of each segmented object. Defaults to `FALSE`. Supported aliases: `points`, `point`, `show_points`.
#' @param point_size Numeric. Size (cex) of centroid points when `show_point = TRUE`. Defaults to `1.5`. Supported aliases: `marker_size`, `point_cex`.
#' @param point_pch Integer. Symbol code for centroid points. Defaults to `16` (filled circle). Supported alias: `marker_pch`.
#' @param point_col Color for centroid points. Defaults to `NULL` (uses matching instance color). Supported alias: `marker_col`.
#' @param mask Logical. If `TRUE` (default), draws instance segmentation masks.
#' @param cex Numeric character expansion for labels and badges. Defaults to `1.0`.
#' @param pad Numeric padding factor for badges. Defaults to `1.0`.
#' @param show_text Logical. Display text labels on bounding boxes. Defaults to `TRUE`.
#' @param show_conf Logical. Display confidence score on labels/badges. Defaults to `TRUE`.
#' @param show_class Logical. Display class name on labels/badges. Defaults to `TRUE`.
#' @param show_id Logical. Display instance ID number. Defaults to `FALSE`.
#' @param badge Logical. Display rounded badge background behind text labels. Defaults to `TRUE`.
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
                               classes = NULL,
                               agnostic_nms = FALSE,
                               seam_merging = TRUE,
                               full_image_inference = TRUE,
                               type = c("highlight", "segment", "mask"),
                               col_background = "white",
                               col_highlight = "salmon",
                               alpha = 0.40,
                               border = NA,
                               lwd = 1,
                               bbox = FALSE,
                               show_box = NULL,
                               show_point = FALSE,
                               point_size = 1.5,
                               point_pch = 16,
                               point_col = NULL,
                               mask = TRUE,
                               lower_size = NULL,
                               upper_size = NULL,
                               size_mode = c("relative", "absolute"),
                               cex = 1.0,
                               pad = 1.0,
                               show_text = TRUE,
                               show_conf = TRUE,
                               show_class = TRUE,
                               show_id = FALSE,
                               badge = TRUE,
                               rainbow = TRUE,
                               plot = TRUE,
                               verbose = TRUE,
                               imgsz = NULL,
                               ...) {
  dots <- list(...)
  if (is.null(imgsz)) {
    if ("imgz" %in% names(dots)) imgsz <- dots$imgz
    if ("img_size" %in% names(dots)) imgsz <- dots$img_size
    if ("image_size" %in% names(dots)) imgsz <- dots$image_size
  }
  if ("lower_size" %in% names(dots)) lower_size <- dots$lower_size
  if ("upper_size" %in% names(dots)) upper_size <- dots$upper_size
  if ("size_mode" %in% names(dots)) size_mode <- dots$size_mode
  if ("show_points" %in% names(dots)) show_point <- isTRUE(dots$show_points)
  if ("show_point" %in% names(dots)) show_point <- isTRUE(dots$show_point)
  if ("points" %in% names(dots)) show_point <- isTRUE(dots$points)
  if ("point" %in% names(dots)) show_point <- isTRUE(dots$point)
  if ("show_boxes" %in% names(dots)) show_box <- isTRUE(dots$show_boxes)
  if ("box" %in% names(dots)) show_box <- isTRUE(dots$box)
  if ("boxes" %in% names(dots)) show_box <- isTRUE(dots$boxes)
  if (is.null(show_box)) {
    show_box <- if (isTRUE(show_point)) FALSE else isTRUE(bbox)
  }
  if ("marker_size" %in% names(dots)) point_size <- dots$marker_size
  if ("point_cex" %in% names(dots)) point_size <- dots$point_cex
  if ("marker_col" %in% names(dots)) point_col <- dots$marker_col
  if ("marker_pch" %in% names(dots)) point_pch <- dots$marker_pch
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)
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
  raw_im <- if (isS4(im) && methods::.hasSlot(im, ".Data")) methods::slot(im, ".Data") else im

  for (i in seq_len(n_slices)) {
    sx1 <- grid$x[i]
    sy1 <- grid$y[i]
    sx2 <- min(orig_w, sx1 + slice_w - 1L)
    sy2 <- min(orig_h, sy1 + slice_h - 1L)

    sub_im <- if (length(dims) >= 3L) raw_im[sx1:sx2, sy1:sy2, , drop = FALSE] else raw_im[sx1:sx2, sy1:sy2, drop = FALSE]

    det <- try(
      image_segment_dl(
        sub_im,
        model = model,
        type = "highlight",
        box_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        classes = classes,
        plot = FALSE,
        verbose = FALSE,
        imgsz = imgsz,
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

    df_b$is_full <- FALSE
    all_candidates[[length(all_candidates) + 1L]] <- df_b
    all_contours <- c(all_contours, shifted_c)
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
        classes = classes,
        plot = FALSE,
        verbose = FALSE,
        imgsz = imgsz,
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
        n_f <- nrow(f_df)
        f_shifted_c <- vector("list", n_f)
        for (k in seq_len(n_f)) {
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
          f_shifted_c[[k]] <- if (valid_cnt) {
            f_conts[[k]]
          } else {
            matrix(c(f_df$xmin[k], f_df$xmax[k], f_df$xmax[k], f_df$xmin[k], f_df$ymin[k], f_df$ymin[k], f_df$ymax[k], f_df$ymax[k]), ncol = 2)
          }
        }
        f_df$is_full <- TRUE
        all_candidates[[length(all_candidates) + 1L]] <- f_df
        all_contours <- c(all_contours, f_shifted_c)
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

  # Filter by target classes if specified
  if (!is.null(classes) && nrow(merged_df) > 0L) {
    c_match <- tolower(as.character(classes))
    idx_keep <- tolower(merged_df$class_name) %in% c_match |
                tolower(merged_df$label) %in% c_match
    if (is.numeric(classes)) {
      idx_keep <- idx_keep | merged_df$class_id %in% classes
    }
    merged_df <- merged_df[idx_keep, , drop = FALSE]
    all_contours <- all_contours[idx_keep]
  }

  if (nrow(merged_df) == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No objects matched the requested class(es).")
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

  # Cross-Slice SAHI Post-Processing (Greedy NMM / NMS)
  res_post <- .sahi_postprocess(
    candidates_df = merged_df,
    contours = all_contours,
    match_metric = match_metric,
    match_threshold = match_threshold,
    postprocess_type = postprocess_type,
    agnostic_nms = agnostic_nms,
    seam_merging = seam_merging,
    slice_borders_x = slice_borders_x,
    slice_borders_y = slice_borders_y
  )

  final_df <- res_post$boxes
  final_conts <- res_post$contours
  rownames(final_df) <- NULL

  if ((!is.null(lower_size) || !is.null(upper_size)) && nrow(final_df) > 0L) {
    size_mode <- match.arg(size_mode)
    f_res <- .filter_by_size(
      boxes = final_df,
      contours = final_conts,
      lower_size = lower_size,
      upper_size = upper_size,
      size_mode = size_mode
    )
    final_df <- f_res$boxes
    final_conts <- f_res$contours
  }

  c_names <- if (!is.null(final_df$class_name)) final_df$class_name else if (!is.null(final_df$label)) final_df$label else rep("object", nrow(final_df))
  tbl <- table(c_names)
  counts <- setNames(as.integer(tbl), names(tbl))

  if (is.null(final_df$id)) final_df$id <- seq_len(nrow(final_df))
  if (is.null(final_df$label)) final_df$label <- c_names
  if (is.null(final_df$score)) {
    final_df$score <- if (!is.null(final_df$conf)) final_df$conf else rep(1.0, nrow(final_df))
  }

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

    plot(im)

    if (type == "highlight") {
      if (isTRUE(mask)) {
        for (k in seq_len(n_inst)) {
          cnt <- final_conts[[k]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            c_fill <- grDevices::adjustcolor(palette_cols[k], alpha.f = alpha)
            c_line <- if (!is.na(border)) border else palette_cols[k]
            graphics::polygon(cnt[, 1], cnt[, 2], col = c_fill, border = c_line, lwd = lwd)
          }
        }
      }
    } else if (type == "segment") {
      if (isTRUE(mask)) {
        for (k in seq_len(n_inst)) {
          cnt <- final_conts[[k]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            c_line <- palette_cols[k]
            c_fill <- grDevices::adjustcolor(c_line, alpha.f = alpha)
            graphics::polygon(cnt[, 1], cnt[, 2], col = c_fill, border = c_line, lwd = max(1.5, lwd))
          }
        }
      }
    } else if (type == "mask") {
      if (isTRUE(mask)) {
        for (k in seq_len(n_inst)) {
          cnt <- final_conts[[k]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            graphics::polygon(cnt[, 1], cnt[, 2], col = "white", border = "white", lwd = 1)
          }
        }
      }
    }

    if (isTRUE(show_box) || isTRUE(show_point) || isTRUE(bbox) || !isTRUE(mask)) {
      .plot_yolo_bboxes(
        boxes = final_df,
        palette_colors = palette_cols,
        lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd,
        cex = cex,
        pad = pad,
        show_text = show_text,
        show_conf = show_conf,
        show_class = show_class,
        show_id = show_id,
        show_box = show_box,
        show_point = show_point,
        point_size = point_size,
        point_pch = point_pch,
        point_col = point_col,
        badge = badge
      )
    } else if (isTRUE(show_id)) {
      for (k in seq_len(n_inst)) {
        cnt <- final_conts[[k]]
        if (is.matrix(cnt) && nrow(cnt) >= 3L) {
          cx <- mean(cnt[, 1]); cy <- mean(cnt[, 2])
          graphics::points(cx, cy, pch = 21, bg = "black", col = "white", cex = 2.2)
          graphics::text(cx, cy, labels = as.character(k), col = "white", font = 2, cex = 0.9)
        }
      }
    }
  }

  invisible(final_df)
}
