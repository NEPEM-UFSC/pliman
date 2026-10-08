# ==============================================================================
# SAHI: SLICING AIDED HYPER INFERENCE FOR GIGAPIXEL & DRONE IMAGERY
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Slicing Aided Hyper Inference (SAHI) for Object Detection
#' @name image_detect_sahi
#' @description
#' `image_detect_sahi()` performs high-resolution sliding-window hyper inference
#' on ultra-large images, drone orthomosaics, and macro scans. Large images (e.g. 20–50 Megapixels)
#' are dynamically sliced into overlapping tiles of size `slice_width x slice_height`.
#'
#' Deep learning inference is executed per slice, and all bounding boxes are transformed
#' to global coordinates and merged using Global Cross-Slice Non-Maximum Suppression (NMS).
#' This prevents small objects (seeds, aphids, stomata, seedlings) from vanishing when downscaled.
#'
#' @param img An `image` object, file path, or 3D array.
#' @param model YOLO ONNX model name or file path (e.g., `"yolo26n"`, `"yolov8s.onnx"`).
#' @param slice_height Integer. Height of sliding window slices. Defaults to `640`.
#' @param slice_width Integer. Width of sliding window slices. Defaults to `640`.
#' @param overlap_height_ratio Numeric scalar between 0 and 0.5. Height overlap fraction. Defaults to `0.20`.
#' @param overlap_width_ratio Numeric scalar between 0 and 0.5. Width overlap fraction. Defaults to `0.20`.
#' @param conf_threshold Numeric scalar between 0 and 1. Minimum detection confidence. Defaults to `0.25`.
#' @param iou_threshold Numeric scalar between 0 and 1. IoU threshold for cross-slice NMS. Defaults to `0.45`.
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
                              full_image_inference = TRUE,
                              rainbow = TRUE,
                              plot = TRUE,
                              verbose = TRUE,
                              ...) {
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
    if (is.list(det) && !is.data.frame(det) && !is.null(det$boxes)) {
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
      if (is.list(full_det) && !is.data.frame(full_det) && !is.null(full_det$boxes)) {
        full_det <- full_det$boxes
      }
      if (is.data.frame(full_det) && nrow(full_det) > 0L) {
        if (is.null(full_det$conf) && !is.null(full_det$score)) full_det$conf <- full_det$score
        if (is.null(full_det$score) && !is.null(full_det$conf)) full_det$score <- full_det$conf
        if (is.null(full_det$class_name) && !is.null(full_det$label)) full_det$class_name <- full_det$label
        if (is.null(full_det$label) && !is.null(full_det$class_name)) full_det$label <- full_det$class_name
        if (is.null(full_det$class_id)) full_det$class_id <- 0L
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

  # 3. Global Cross-Slice Non-Maximum Suppression (NMS)
  keep_indices <- integer()
  unique_classes <- unique(merged_df$class_id)

  for (cls in unique_classes) {
    idx_cls <- which(merged_df$class_id == cls)
    sub_df <- merged_df[idx_cls, , drop = FALSE]

    # Order by confidence / score descending safely
    c_vals <- if (!is.null(sub_df$conf)) sub_df$conf else if (!is.null(sub_df$score)) sub_df$score else rep(1.0, nrow(sub_df))
    ord <- order(c_vals, decreasing = TRUE)
    idx_cls <- idx_cls[ord]
    sub_df <- sub_df[ord, , drop = FALSE]

    while (nrow(sub_df) > 0L) {
      best_idx <- idx_cls[1]
      keep_indices <- c(keep_indices, best_idx)

      if (nrow(sub_df) == 1L) break

      # Compute IoU between best and rest
      bx1 <- sub_df$xmin[1]; by1 <- sub_df$ymin[1]
      bx2 <- sub_df$xmax[1]; by2 <- sub_df$ymax[1]
      b_area <- (bx2 - bx1 + 1) * (by2 - by1 + 1)

      rest_x1 <- sub_df$xmin[-1]; rest_y1 <- sub_df$ymin[-1]
      rest_x2 <- sub_df$xmax[-1]; rest_y2 <- sub_df$ymax[-1]
      rest_area <- (rest_x2 - rest_x1 + 1) * (rest_y2 - rest_y1 + 1)

      ix1 <- pmax(bx1, rest_x1); iy1 <- pmax(by1, rest_y1)
      ix2 <- pmin(bx2, rest_x2); iy2 <- pmin(by2, rest_y2)

      iw <- pmax(0, ix2 - ix1 + 1)
      ih <- pmax(0, iy2 - iy1 + 1)
      inter_area <- iw * ih

      union_area <- b_area + rest_area - inter_area
      iou <- inter_area / pmax(1, union_area)

      keep_rest <- which(iou < iou_threshold)
      if (length(keep_rest) == 0L) break

      sub_df <- sub_df[-1, , drop = FALSE][keep_rest, , drop = FALSE]
      idx_cls <- idx_cls[-1][keep_rest]
    }
  }

  final_df <- merged_df[sort(keep_indices), , drop = FALSE]
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
#' shifts polygon contours to global coordinates, and applies cross-slice NMS.
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
        conf_threshold = conf_threshold,
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
    if (is.list(det) && !is.data.frame(det)) {
      if (!is.null(det$boxes)) df_b <- det$boxes
      if (!is.null(det$contours)) conts <- det$contours
    } else if (is.data.frame(det)) {
      df_b <- det
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
      if (length(conts) >= k && is.matrix(conts[[k]]) && nrow(conts[[k]]) >= 3L) {
        c_mat <- conts[[k]]
        c_mat[, 1] <- c_mat[, 1] + sx1 - 1L
        c_mat[, 2] <- c_mat[, 2] + sy1 - 1L
        shifted_c[[k]] <- c_mat
      } else {
        # Fallback to rectangular polygon from box
        bx1 <- df_b$xmin[k]; by1 <- df_b$ymin[k]
        bx2 <- df_b$xmax[k]; by2 <- df_b$ymax[k]
        shifted_c[[k]] <- matrix(c(bx1, bx2, bx2, bx1, by1, by1, by2, by2), ncol = 2)
      }
    }

    for (k in seq_len(n_b)) {
      all_candidates[[length(all_candidates) + 1L]] <- df_b[k, , drop = FALSE]
      all_contours[[length(all_contours) + 1L]] <- shifted_c[[k]]
    }
  }

  # Optional full-image global pass
  if (isTRUE(full_image_inference)) {
    full_det <- try(
      image_segment_dl(
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
      f_df <- NULL
      f_conts <- list()
      if (is.list(full_det) && !is.data.frame(full_det)) {
        if (!is.null(full_det$boxes)) f_df <- full_det$boxes
        if (!is.null(full_det$contours)) f_conts <- full_det$contours
      } else if (is.data.frame(full_det)) {
        f_df <- full_det
      }
      if (!is.null(f_df) && is.data.frame(f_df) && nrow(f_df) > 0L) {
        if (is.null(f_df$conf) && !is.null(f_df$score)) f_df$conf <- f_df$score
        if (is.null(f_df$score) && !is.null(f_df$conf)) f_df$score <- f_df$conf
        if (is.null(f_df$class_name) && !is.null(f_df$label)) f_df$class_name <- f_df$label
        if (is.null(f_df$label) && !is.null(f_df$class_name)) f_df$label <- f_df$class_name
        if (is.null(f_df$class_id)) f_df$class_id <- 0L
        for (k in seq_len(nrow(f_df))) {
          all_candidates[[length(all_candidates) + 1L]] <- f_df[k, , drop = FALSE]
          c_mat <- if (length(f_conts) >= k && is.matrix(f_conts[[k]])) {
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

  # Cross-Slice NMS on merged boxes and contours
  keep_indices <- integer()
  unique_classes <- unique(merged_df$class_id)

  for (cls in unique_classes) {
    idx_cls <- which(merged_df$class_id == cls)
    sub_df <- merged_df[idx_cls, , drop = FALSE]

    c_vals <- if (!is.null(sub_df$conf)) sub_df$conf else if (!is.null(sub_df$score)) sub_df$score else rep(1.0, nrow(sub_df))
    ord <- order(c_vals, decreasing = TRUE)
    idx_cls <- idx_cls[ord]
    sub_df <- sub_df[ord, , drop = FALSE]

    while (nrow(sub_df) > 0L) {
      best_idx <- idx_cls[1]
      keep_indices <- c(keep_indices, best_idx)

      if (nrow(sub_df) == 1L) break

      bx1 <- sub_df$xmin[1]; by1 <- sub_df$ymin[1]
      bx2 <- sub_df$xmax[1]; by2 <- sub_df$ymax[1]
      b_area <- (bx2 - bx1 + 1) * (by2 - by1 + 1)

      rest_x1 <- sub_df$xmin[-1]; rest_y1 <- sub_df$ymin[-1]
      rest_x2 <- sub_df$xmax[-1]; rest_y2 <- sub_df$ymax[-1]
      rest_area <- (rest_x2 - rest_x1 + 1) * (rest_y2 - rest_y1 + 1)

      ix1 <- pmax(bx1, rest_x1); iy1 <- pmax(by1, rest_y1)
      ix2 <- pmin(bx2, rest_x2); iy2 <- pmin(by2, rest_y2)

      iw <- pmax(0, ix2 - ix1 + 1)
      ih <- pmax(0, iy2 - iy1 + 1)
      inter_area <- iw * ih

      union_area <- b_area + rest_area - inter_area
      iou <- inter_area / pmax(1, union_area)

      keep_rest <- which(iou < iou_threshold)
      if (length(keep_rest) == 0L) break

      sub_df <- sub_df[-1, , drop = FALSE][keep_rest, , drop = FALSE]
      idx_cls <- idx_cls[-1][keep_rest]
    }
  }

  sorted_keep <- sort(keep_indices)
  final_df <- merged_df[sorted_keep, , drop = FALSE]
  final_conts <- all_contours[sorted_keep]
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
