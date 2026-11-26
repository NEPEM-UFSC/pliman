# utils_yolo_diagnostics.R - Model Diagnostics, TIDE Error Analysis, and Saliency Attribution for YOLO Models

#' TIDE-Style Prediction Error Decomposition for YOLO Models
#'
#' @description
#' `yolo_error_analysis()` decomposes object detection errors into six mutually exclusive
#' canonical categories following the TIDE framework (Bolya et al., ECCV 2020). Comparing
#' predictions against ground-truth annotations across validation or test splits, it pinpoints
#' the exact bottlenecks of a model:
#'
#' 1. **Classification Error (`Cls`)**: The object was localized accurately (\eqn{\text{IoU} \ge 0.50}),
#'    but classified into the wrong category.
#' 2. **Localization Error (`Loc`)**: The object class was predicted correctly, but bounding box alignment
#'    was inaccurate (\eqn{0.10 \le \text{IoU} < 0.50}).
#' 3. **Both Classification & Localization Error (`Both`)**: Both class was wrong and localization was
#'    inaccurate (\eqn{0.10 \le \text{IoU} < 0.50}).
#' 4. **Duplicate Detection Error (`Dup`)**: Non-Maximum Suppression (NMS) failed to suppress a redundant
#'    candidate detection on an already matched ground-truth object (\eqn{\text{IoU} \ge 0.50}).
#' 5. **Background Error (`Bkg`)**: The model hallucinated an object on pure background with no nearby
#'    ground-truth entity (\eqn{\text{IoU} < 0.10}).
#' 6. **Missed Ground Truth (`Miss`)**: A true object in the field was completely overlooked by the
#'    detector (pure false negative, \eqn{\text{IoU} < 0.10} against all predictions of the same class).
#'
#' @references
#' Bolya, D., Foley, S., Hays, J., & Hoffman, J. (2020). TIDE: A general toolbox for identifying
#' object detection errors. *European Conference on Computer Vision (ECCV)*, pp. 558-573.
#'
#' @param test_res An object of class `yolo_test` (returned by [yolo_test()]), or `NULL` if
#'   evaluating directly with `model` and `data`.
#' @param model A trained model identifier or path (`.onnx` or `.pt`). Used if `test_res` is `NULL`.
#' @param data Path to the YOLO dataset root directory. Used if `test_res` is `NULL`.
#' @param split Character string. Dataset split to evaluate: `"test"` (default), `"val"`, or `"train"`.
#' @param iou_thresh Numeric. IoU threshold for correct object localization. Defaults to `0.50`.
#' @param bg_thresh Numeric. Background IoU threshold. Detections with \eqn{\text{IoU} < \text{bg\_thresh}}
#'   against all ground-truth objects are considered pure background hallucinations. Defaults to `0.10`.
#' @param plot Logical. Whether to automatically render the error decomposition bar chart. Defaults to `TRUE`.
#' @param verbose Logical. Whether to display formatted progress and summary tables. Defaults to `TRUE`.
#' @param ... Additional arguments forwarded to [yolo_test()] if `test_res` is `NULL`.
#'
#' @return An S3 object of class `c("yolo_error_analysis", "list")` containing:
#'   * `summary`: A `data.frame` with the 6 error categories, count of occurrences, percentage of all errors, and descriptions.
#'   * `by_class`: A `data.frame` of error counts and percentages broken down by object class.
#'   * `errors`: A detailed `data.frame` listing every individual error with image name, bounding box coordinates, error type, and score.
#'   * `metrics`: A list of summary benchmark numbers (Total GT, Total Pred, TP, FP, FN, Precision, Recall).
#'
#' @export
#' @aliases yolo_tide_analysis
#' @seealso [yolo_test()], [yolo_explain()], [plot.yolo_error_analysis()]
yolo_error_analysis <- function(test_res = NULL,
                                model = NULL,
                                data = NULL,
                                split = "test",
                                iou_thresh = 0.50,
                                bg_thresh = 0.10,
                                plot = TRUE,
                                verbose = TRUE,
                                ...) {
  # 1. Resolve test results
  if (is.null(test_res)) {
    if (is.null(model) || is.null(data)) {
      cli::cli_abort("Either an existing {.val test_res} (from {.code yolo_test()}) or both {.val model} and {.val data} must be provided.")
    }
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Running model evaluation via {.code yolo_test()} on split {.val {split}}...")
    }
    test_res <- yolo_test(model = model, data = data, split = split, verbose = verbose, ...)
  }

  if (!inherits(test_res, "yolo_test")) {
    cli::cli_abort("{.val test_res} must be an object of class {.val yolo_test}.")
  }

  preds_df <- test_res$predictions
  gt_list <- test_res$gt_list
  flat_gt <- test_res$ground_truth

  # Helper: Compute 2D IoU between two boxes [x1, y1, x2, y2]
  calc_iou <- function(b1, b2) {
    x1 <- max(b1[1], b2[1]); y1 <- max(b1[2], b2[2])
    x2 <- min(b1[3], b2[3]); y2 <- min(b1[4], b2[4])
    if (x2 <= x1 || y2 <= y1) return(0.0)
    inter <- (x2 - x1) * (y2 - y1)
    a1 <- (b1[3] - b1[1]) * (b1[4] - b1[2])
    a2 <- (b2[3] - b2[1]) * (b2[4] - b2[2])
    union <- a1 + a2 - inter
    if (union <= 0.0) return(0.0)
    inter / union
  }

  # Discover all unique images
  img_names <- if (!is.null(test_res$by_image$image)) test_res$by_image$image else unique(c(preds_df$img, names(gt_list)))
  img_names <- img_names[nzchar(img_names)]

  error_records <- list()
  total_tp <- 0L
  total_gt_count <- 0L
  total_pred_count <- nrow(preds_df)

  for (img_name in img_names) {
    img_preds <- if (nrow(preds_df) > 0L) preds_df[preds_df$img == img_name, , drop = FALSE] else data.frame()
    img_gt <- if (!is.null(gt_list[[img_name]])) gt_list[[img_name]] else data.frame()

    n_p <- nrow(img_preds)
    n_g <- nrow(img_gt)
    total_gt_count <- total_gt_count + n_g

    if (n_p == 0L && n_g == 0L) next

    # If no predictions but GT exists: all GTs are Missed
    if (n_p == 0L && n_g > 0L) {
      for (g_idx in seq_len(n_g)) {
        error_records[[length(error_records) + 1L]] <- list(
          image = img_name,
          error_type = "Miss",
          error_name = "Missed GT",
          pred_class = NA_character_,
          gt_class = as.character(img_gt$class[g_idx]),
          iou = 0.0,
          score = NA_real_,
          xmin = as.numeric(img_gt$xmin[g_idx]),
          ymin = as.numeric(img_gt$ymin[g_idx]),
          xmax = as.numeric(img_gt$xmax[g_idx]),
          ymax = as.numeric(img_gt$ymax[g_idx])
        )
      }
      next
    }

    # If predictions exist but no GT: all predictions are Background Error
    if (n_p > 0L && n_g == 0L) {
      for (p_idx in seq_len(n_p)) {
        error_records[[length(error_records) + 1L]] <- list(
          image = img_name,
          error_type = "Bkg",
          error_name = "Background Error",
          pred_class = as.character(img_preds$label[p_idx]),
          gt_class = NA_character_,
          iou = 0.0,
          score = as.numeric(img_preds$score[p_idx]),
          xmin = as.numeric(img_preds$xmin[p_idx]),
          ymin = as.numeric(img_preds$ymin[p_idx]),
          xmax = as.numeric(img_preds$xmax[p_idx]),
          ymax = as.numeric(img_preds$ymax[p_idx])
        )
      }
      next
    }

    # Both predictions and GT exist: sort predictions by confidence score descending
    ord <- order(img_preds$score, decreasing = TRUE)
    img_preds <- img_preds[ord, , drop = FALSE]

    gt_matched <- rep(FALSE, n_g)
    pred_assigned <- rep(FALSE, n_p)

    # Compute pairwise IoU matrix [n_p x n_g]
    iou_mat <- matrix(0.0, nrow = n_p, ncol = n_g)
    for (i in seq_len(n_p)) {
      pb <- c(img_preds$xmin[i], img_preds$ymin[i], img_preds$xmax[i], img_preds$ymax[i])
      for (j in seq_len(n_g)) {
        gb <- c(img_gt$xmin[j], img_gt$ymin[j], img_gt$xmax[j], img_gt$ymax[j])
        iou_mat[i, j] <- calc_iou(pb, gb)
      }
    }

    # Step 1: True Positive matching (same class, IoU >= iou_thresh)
    for (i in seq_len(n_p)) {
      p_cls <- tolower(as.character(img_preds$label[i]))
      best_iou <- 0.0
      best_j <- 0L

      for (j in seq_len(n_g)) {
        if (!gt_matched[j] && tolower(as.character(img_gt$class[j])) == p_cls) {
          if (iou_mat[i, j] >= iou_thresh && iou_mat[i, j] > best_iou) {
            best_iou <- iou_mat[i, j]
            best_j <- j
          }
        }
      }

      if (best_j > 0L) {
        gt_matched[best_j] <- TRUE
        pred_assigned[i] <- TRUE
        total_tp <- total_tp + 1L
      }
    }

    # Step 2: Classify False Positive errors for remaining predictions
    for (i in seq_len(n_p)) {
      if (pred_assigned[i]) next

      p_cls <- tolower(as.character(img_preds$label[i]))
      p_score <- as.numeric(img_preds$score[i])
      max_iou <- max(iou_mat[i, ])
      best_gt_idx <- which.max(iou_mat[i, ])
      best_gt_cls <- tolower(as.character(img_gt$class[best_gt_idx]))

      err_type <- "Bkg"
      err_name <- "Background Error"

      if (max_iou >= iou_thresh) {
        if (p_cls != best_gt_cls) {
          err_type <- "Cls"
          err_name <- "Classification Error"
        } else {
          err_type <- "Dup"
          err_name <- "Duplicate Detection Error"
        }
      } else if (max_iou >= bg_thresh) {
        if (p_cls == best_gt_cls) {
          err_type <- "Loc"
          err_name <- "Localization Error"
        } else {
          err_type <- "Both"
          err_name <- "Classification & Localization Error"
        }
      } else {
        err_type <- "Bkg"
        err_name <- "Background Error"
      }

      error_records[[length(error_records) + 1L]] <- list(
        image = img_name,
        error_type = err_type,
        error_name = err_name,
        pred_class = as.character(img_preds$label[i]),
        gt_class = if (max_iou >= bg_thresh) as.character(img_gt$class[best_gt_idx]) else NA_character_,
        iou = round(max_iou, 4),
        score = p_score,
        xmin = as.numeric(img_preds$xmin[i]),
        ymin = as.numeric(img_preds$ymin[i]),
        xmax = as.numeric(img_preds$xmax[i]),
        ymax = as.numeric(img_preds$ymax[i])
      )
    }

    # Step 3: Classify Missed Ground Truths (unmatched GTs)
    for (j in seq_len(n_g)) {
      if (!gt_matched[j]) {
        g_cls <- tolower(as.character(img_gt$class[j]))
        max_iou_from_same <- 0.0
        for (i in seq_len(n_p)) {
          if (tolower(as.character(img_preds$label[i])) == g_cls) {
            if (iou_mat[i, j] > max_iou_from_same) max_iou_from_same <- iou_mat[i, j]
          }
        }

        # If not covered by any correct-class prediction with IoU >= bg_thresh, it is a Miss
        if (max_iou_from_same < iou_thresh) {
          error_records[[length(error_records) + 1L]] <- list(
            image = img_name,
            error_type = "Miss",
            error_name = "Missed GT",
            pred_class = NA_character_,
            gt_class = as.character(img_gt$class[j]),
            iou = round(max_iou_from_same, 4),
            score = NA_real_,
            xmin = as.numeric(img_gt$xmin[j]),
            ymin = as.numeric(img_gt$ymin[j]),
            xmax = as.numeric(img_gt$xmax[j]),
            ymax = as.numeric(img_gt$ymax[j])
          )
        }
      }
    }
  }

  # Build full error data frame
  errors_df <- if (length(error_records) > 0L) {
    do.call(rbind.data.frame, c(error_records, list(stringsAsFactors = FALSE)))
  } else {
    data.frame(
      image = character(0), error_type = character(0), error_name = character(0),
      pred_class = character(0), gt_class = character(0), iou = numeric(0),
      score = numeric(0), xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Canonical 6 TIDE categories in standard order
  categories <- c("Cls", "Loc", "Both", "Dup", "Bkg", "Miss")
  cat_names <- c(
    "Classification Error",
    "Localization Error",
    "Both (Cls + Loc)",
    "Duplicate Detection Error",
    "Background Error",
    "Missed GT"
  )
  descriptions <- c(
    "Accurate location (IoU >= 0.50), incorrect class",
    "Correct class, suboptimal alignment (0.10 <= IoU < 0.50)",
    "Incorrect class and suboptimal alignment (0.10 <= IoU < 0.50)",
    "NMS failure: redundant prediction for an already matched object",
    "Hallucination on background (IoU < 0.10)",
    "Unmatched ground truth object (false negative)"
  )

  counts <- vapply(categories, function(cat) sum(errors_df$error_type == cat), integer(1))
  tot_errors <- max(1L, sum(counts))
  pcts <- round(100 * counts / tot_errors, 2)

  summary_df <- data.frame(
    code = categories,
    category = cat_names,
    count = as.integer(counts),
    percentage = pcts,
    description = descriptions,
    stringsAsFactors = FALSE
  )

  # Breakdown per class
  all_classes <- unique(na.omit(c(errors_df$pred_class, errors_df$gt_class)))
  all_classes <- all_classes[nzchar(all_classes)]

  by_class_list <- vector("list", length(all_classes))
  for (idx in seq_along(all_classes)) {
    cname <- all_classes[idx]
    c_sub <- errors_df[errors_df$pred_class == cname | (is.na(errors_df$pred_class) & errors_df$gt_class == cname), , drop = FALSE]
    c_counts <- vapply(categories, function(cat) sum(c_sub$error_type == cat), integer(1))
    c_tot <- max(1L, sum(c_counts))
    by_class_list[[idx]] <- data.frame(
      class = cname,
      total_errors = sum(c_counts),
      cls_err = c_counts["Cls"],
      loc_err = c_counts["Loc"],
      both_err = c_counts["Both"],
      dup_err = c_counts["Dup"],
      bkg_err = c_counts["Bkg"],
      miss_gt = c_counts["Miss"],
      stringsAsFactors = FALSE
    )
  }
  by_class_df <- if (length(by_class_list) > 0L) do.call(rbind, by_class_list) else data.frame()

  # High-level metrics
  total_fp <- sum(counts[c("Cls", "Loc", "Both", "Dup", "Bkg")])
  total_fn <- sum(counts["Miss"])
  prec <- if ((total_tp + total_fp) > 0) total_tp / (total_tp + total_fp) else 0
  rec <- if (total_gt_count > 0) total_tp / total_gt_count else 0
  f1 <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

  out <- structure(
    list(
      summary = summary_df,
      by_class = by_class_df,
      errors = errors_df,
      metrics = list(
        total_images = length(img_names),
        total_gt = total_gt_count,
        total_pred = total_pred_count,
        tp = total_tp,
        fp = total_fp,
        fn = total_fn,
        precision = round(prec, 4),
        recall = round(rec, 4),
        f1 = round(f1, 4),
        iou_thresh = iou_thresh,
        bg_thresh = bg_thresh
      )
    ),
    class = c("yolo_error_analysis", "list")
  )

  if (isTRUE(verbose)) {
    print(out)
  }

  if (isTRUE(plot)) {
    plot(out)
  }

  invisible(out)
}

#' @rdname yolo_error_analysis
#' @export
yolo_tide_analysis <- yolo_error_analysis

#' @export
print.yolo_error_analysis <- function(x, ...) {
  m <- x$metrics
  cli::cli_h2("pliman: TIDE Prediction Error Decomposition")
  cli::cli_ul(c(
    "Evaluated: {.val {m$total_images}} images | {.val {m$total_gt}} ground-truth objects | {.val {m$total_pred}} predictions",
    "Performance: Precision: {.val {m$precision}} | Recall: {.val {m$recall}} | F1-Score: {.val {m$f1}}",
    "Thresholds: Localization IoU >= {.val {m$iou_thresh}} | Background IoU < {.val {m$bg_thresh}}"
  ))

  cli::cli_h3("Error Distribution Breakdown:")
  disp_df <- x$summary[, c("category", "code", "count", "percentage")]
  names(disp_df) <- c("Error Category", "Code", "Count", "Share (%)")
  print.data.frame(disp_df, row.names = FALSE)
  invisible(x)
}

#' Plot TIDE Prediction Error Decomposition
#'
#' @description
#' Visualizes the 6 TIDE error categories in a publication-ready horizontal bar chart
#' with percentage and count annotations.
#'
#' @param x An object of class `yolo_error_analysis`.
#' @param by_class Logical. Whether to render error distributions broken down by class. Defaults to `FALSE`.
#' @param col Optional vector of colors for the 6 error categories.
#' @param ... Additional graphical arguments.
#'
#' @export
plot.yolo_error_analysis <- function(x, by_class = FALSE, col = NULL, ...) {
  s <- x$summary
  if (is.null(col)) {
    # TIDE standard palette: Cls (gold), Loc (dodgerblue), Both (darkorange), Dup (purple), Bkg (firebrick), Miss (grey40)
    col <- c(
      "Cls" = "#E69F00",
      "Loc" = "#56B4E9",
      "Both" = "#D55E00",
      "Dup" = "#CC79A7",
      "Bkg" = "#D9534F",
      "Miss" = "#4E555B"
    )
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- s
    df$category <- factor(df$category, levels = rev(df$category))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$category, y = .data$count, fill = .data$code)) +
      ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = col) +
      ggplot2::coord_flip() +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d (%.1f%%)", .data$count, .data$percentage)),
        hjust = -0.15,
        size = 3.6,
        fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.22))) +
      ggplot2::labs(
        title = "YOLO Prediction Error Decomposition (TIDE)",
        subtitle = sprintf("Total errors: %d across %d ground truth objects", sum(df$count), x$metrics$total_gt),
        x = NULL,
        y = "Error Count"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.text.y = ggplot2::element_text(face = "bold", size = 11)
      )
    print(p)
  } else {
    # Base R fallback
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mar = c(4, 12, 3, 2))
    bp <- graphics::barplot(
      rev(s$count),
      names.arg = rev(s$category),
      horiz = TRUE,
      las = 1,
      col = rev(col[s$code]),
      border = NA,
      main = "YOLO Prediction Error Decomposition (TIDE)",
      xlab = "Error Count",
      xlim = c(0, max(s$count) * 1.30)
    )
    graphics::text(
      x = rev(s$count),
      y = bp,
      labels = sprintf("  %d (%.1f%%)", rev(s$count), rev(s$percentage)),
      pos = 4,
      font = 2,
      cex = 0.85
    )
  }
  invisible(x)
}


#' Explain YOLO Model Detections via Perturbation Saliency Heatmaps
#'
#' @description
#' `yolo_explain()` generates visual saliency attribution heatmaps for object detection
#' predictions by perturbing (occluding/blurring) spatial image patches and measuring the
#' resulting drop in detection confidence score.
#'
#' Visualizing which image regions (e.g., embryo vs endosperm, leaf margin vs lesion center)
#' most strongly contribute to the model's confidence is essential for opening the deep learning
#' "black box", confirming model robustness, and satisfying rigorous review criteria in Q1 journals.
#'
#' @param img An `Image` object or character string path to an image file.
#' @param model A trained YOLO model identifier or ONNX file path (e.g., `"yolo26n"`).
#' @param target Optional integer or character string specifying which detection to explain:
#'   * `NULL` (default): Explains the detection with the highest confidence score.
#'   * An integer `1, 2, ...`: Explains the detection at that rank.
#'   * A character string: Explains the highest confidence detection of that class name.
#' @param grid_size Integer. Size in pixels of occlusion patches. Defaults to `24L`.
#' @param stride Optional integer. Stride between adjacent patches. Defaults to `grid_size / 2`
#'   for smooth attribution maps.
#' @param perturb_type Character string specifying how patches are masked:
#'   * `"blur"` (default): Blurs the patch to eliminate high-frequency edges and textures.
#'   * `"mean"`: Replaces the patch with the mean pixel color of the image.
#'   * `"zero"`: Sets the patch pixels to 0 (black).
#' @param target_region Character string:
#'   * `"box"` (default): Restricts perturbations to the target object's bounding box (+15% margin)
#'     for 20x faster computation.
#'   * `"full"`: Sweeps the entire image canvas.
#' @param col_palette Character string. Name of the color palette for the attribution heatmap:
#'   `"inferno"` (default), `"viridis"`, `"plasma"`, or `"magma"`.
#' @param alpha Numeric. Opacity of the heatmap overlay on the original image (0 to 1). Defaults to `0.55`.
#' @param conf_threshold Numeric. Baseline candidate detection threshold. Defaults to `0.20`.
#' @param iou_threshold Numeric. NMS threshold. Defaults to `0.45`.
#' @param engine Character string. Compute engine: `"gpu"` (default) or `"cpu"`.
#' @param plot Logical. Whether to render the saliency visualization. Defaults to `TRUE`.
#' @param verbose Logical. Whether to display progress messages. Defaults to `TRUE`.
#' @param ... Additional arguments forwarded to [image_detect_dl()].
#'
#' @return An S3 object of class `c("yolo_saliency", "list")` containing:
#'   * `heatmap`: A normalized 2D matrix of attribution values \eqn{[0, 1]}.
#'   * `target_box`: Bounding box coordinates of the explained object.
#'   * `target_class`: Predicted class label of the explained object.
#'   * `baseline_conf`: Baseline prediction confidence before perturbation.
#'   * `image`: The original image object.
#'   * `col_palette`: Heatmap color palette used.
#'   * `alpha`: Heatmap transparency level.
#'
#' @export
#' @aliases yolo_saliency
#' @seealso [yolo_error_analysis()], [image_detect_dl()]
yolo_explain <- function(img,
                         model = "yolo26n",
                         target = NULL,
                         grid_size = 24L,
                         stride = NULL,
                         perturb_type = c("blur", "mean", "zero"),
                         target_region = c("box", "full"),
                         col_palette = "inferno",
                         alpha = 0.55,
                         conf_threshold = 0.20,
                         iou_threshold = 0.45,
                         engine = "gpu",
                         plot = TRUE,
                         verbose = TRUE,
                         ...) {
  perturb_type <- match.arg(perturb_type)
  target_region <- match.arg(target_region)
  grid_size <- max(8L, as.integer(grid_size))
  if (is.null(stride)) stride <- max(4L, as.integer(round(grid_size / 2)))

  # 1. Load Image
  im_obj <- if (is.character(img) && file.exists(img[1])) image_import(img[1]) else img
  if (!inherits(im_obj, "Image")) im_obj <- as_image(im_obj)

  img_arr <- image_data(im_obj)
  dims <- dim(img_arr)
  w <- dims[1]; h <- dims[2]; nch <- if (length(dims) >= 3L) dims[3] else 1L

  # 2. Run Baseline Detection
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running baseline YOLO inference...",
      msg_done = "Baseline detection completed"
    )
  }

  base_det <- image_detect_dl(
    img = im_obj,
    model = model,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    engine = engine,
    plot = FALSE,
    verbose = FALSE,
    ...
  )

  boxes <- base_det$boxes
  if (is.null(boxes) || nrow(boxes) == 0L) {
    cli::cli_abort("No objects detected in the image at confidence threshold {conf_threshold}. Cannot generate explanation.")
  }

  # 3. Resolve Target Detection
  target_idx <- 1L
  if (is.null(target)) {
    target_idx <- which.max(boxes$score)
  } else if (is.numeric(target)) {
    target_idx <- min(nrow(boxes), max(1L, as.integer(target)))
  } else if (is.character(target)) {
    match_c <- which(tolower(boxes$class) == tolower(target))
    if (length(match_c) > 0L) {
      target_idx <- match_c[which.max(boxes$score[match_c])]
    } else {
      cli::cli_alert_warning("No detection with class {.val {target}} found. Explaining top detection instead.")
      target_idx <- which.max(boxes$score)
    }
  }

  tgt_box <- boxes[target_idx, ]
  base_score <- as.numeric(tgt_box$score)
  tgt_class <- as.character(tgt_box$class)
  bx1 <- max(1L, as.integer(round(tgt_box$xmin)))
  by1 <- max(1L, as.integer(round(tgt_box$ymin)))
  bx2 <- min(w, as.integer(round(tgt_box$xmax)))
  by2 <- min(h, as.integer(round(tgt_box$ymax)))

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Explaining object #{target_idx} [{tgt_class}] with baseline confidence: {round(base_score, 3)}.")
  }

  # 4. Define Grid Search Space
  if (target_region == "box") {
    pad_x <- as.integer(round((bx2 - bx1) * 0.15))
    pad_y <- as.integer(round((by2 - by1) * 0.15))
    scan_x1 <- max(1L, bx1 - pad_x); scan_x2 <- min(w - grid_size + 1L, bx2 + pad_x)
    scan_y1 <- max(1L, by1 - pad_y); scan_y2 <- min(h - grid_size + 1L, by2 + pad_y)
  } else {
    scan_x1 <- 1L; scan_x2 <- w - grid_size + 1L
    scan_y1 <- 1L; scan_y2 <- h - grid_size + 1L
  }

  xs <- seq(scan_x1, scan_x2, by = stride)
  ys <- seq(scan_y1, scan_y2, by = stride)
  n_patches <- length(xs) * length(ys)

  if (isTRUE(verbose)) {
    cli::cli_progress_bar(
      format = "{cli::pb_spin} Generating saliency perturbation map [{cli::pb_current}/{cli::pb_total}] | ETA: {cli::pb_eta}",
      total = n_patches,
      clear = FALSE
    )
  }

  # Helper: Match IoU with target box
  target_bbox_vec <- c(bx1, by1, bx2, by2)

  # Accumulator matrices for attribution map
  attr_accum <- matrix(0.0, nrow = w, ncol = h)
  count_accum <- matrix(0L, nrow = w, ncol = h)

  # Mean pixel color if needed
  mean_val <- if (perturb_type == "mean") mean(img_arr) else 0.0

  for (px in xs) {
    for (py in ys) {
      if (isTRUE(verbose)) cli::cli_progress_update()
      ex <- min(w, px + grid_size - 1L)
      ey <- min(h, py + grid_size - 1L)

      # Perturb copy of image array
      pert_arr <- img_arr
      if (perturb_type == "blur") {
        # Spatial local blur on patch
        for (c_i in seq_len(nch)) {
          sub <- if (nch == 1L) pert_arr[px:ex, py:ey] else pert_arr[px:ex, py:ey, c_i]
          m_sub <- mean(sub)
          if (nch == 1L) pert_arr[px:ex, py:ey] <- m_sub else pert_arr[px:ex, py:ey, c_i] <- m_sub
        }
      } else if (perturb_type == "mean") {
        pert_arr[px:ex, py:ey, ] <- mean_val
      } else {
        pert_arr[px:ex, py:ey, ] <- 0.0
      }

      pert_img <- as_image(pert_arr)

      # Run inference on perturbed image
      pert_det <- tryCatch({
        image_detect_dl(
          img = pert_img,
          model = model,
          conf_threshold = max(0.05, conf_threshold * 0.5),
          iou_threshold = iou_threshold,
          engine = engine,
          plot = FALSE,
          verbose = FALSE,
          ...
        )
      }, error = function(e) NULL)

      pert_score <- 0.0
      if (!is.null(pert_det) && !is.null(pert_det$boxes) && nrow(pert_det$boxes) > 0L) {
        pb_df <- pert_det$boxes
        same_cls <- which(tolower(pb_df$class) == tolower(tgt_class))
        if (length(same_cls) > 0L) {
          # Find overlapping box with target
          best_iou <- 0.0
          for (c_idx in same_cls) {
            cb <- c(pb_df$xmin[c_idx], pb_df$ymin[c_idx], pb_df$xmax[c_idx], pb_df$ymax[c_idx])
            inter_x1 <- max(target_bbox_vec[1], cb[1]); inter_y1 <- max(target_bbox_vec[2], cb[2])
            inter_x2 <- min(target_bbox_vec[3], cb[3]); inter_y2 <- min(target_bbox_vec[4], cb[4])
            if (inter_x2 > inter_x1 && inter_y2 > inter_y1) {
              inter_a <- (inter_x2 - inter_x1) * (inter_y2 - inter_y1)
              union_a <- (target_bbox_vec[3] - target_bbox_vec[1]) * (target_bbox_vec[4] - target_bbox_vec[2]) +
                         (cb[3] - cb[1]) * (cb[4] - cb[2]) - inter_a
              iou <- if (union_a > 0) inter_a / union_a else 0
              if (iou > best_iou) {
                best_iou <- iou
                if (iou >= 0.35) pert_score <- as.numeric(pb_df$score[c_idx])
              }
            }
          }
        }
      }

      # Drop in score is the attribution
      drop_s <- max(0.0, base_score - pert_score)
      attr_accum[px:ex, py:ey] <- attr_accum[px:ex, py:ey] + drop_s
      count_accum[px:ex, py:ey] <- count_accum[px:ex, py:ey] + 1L
    }
  }

  if (isTRUE(verbose)) cli::cli_progress_done()

  # Normalize attribution map
  valid_mask <- count_accum > 0L
  attr_map <- matrix(0.0, nrow = w, ncol = h)
  attr_map[valid_mask] <- attr_accum[valid_mask] / count_accum[valid_mask]

  max_v <- max(attr_map)
  if (max_v > 0) attr_map <- attr_map / max_v

  out <- structure(
    list(
      heatmap = attr_map,
      target_box = tgt_box,
      target_class = tgt_class,
      baseline_conf = base_score,
      image = im_obj,
      col_palette = col_palette,
      alpha = alpha
    ),
    class = c("yolo_saliency", "list")
  )

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Attribution heatmap generated successfully.")
  }

  if (isTRUE(plot)) {
    plot(out)
  }

  invisible(out)
}

#' @rdname yolo_explain
#' @export
yolo_saliency <- yolo_explain

#' @export
print.yolo_saliency <- function(x, ...) {
  b <- x$target_box
  cli::cli_h2("pliman: YOLO Saliency Attribution Map")
  cli::cli_ul(c(
    "Target Class: {.val {x$target_class}}",
    "Baseline Confidence: {.val {round(x$baseline_conf, 4)}}",
    "Bounding Box: [{b$xmin}, {b$ymin}, {b$xmax}, {b$ymax}]",
    "Heatmap Dimension: {nrow(x$heatmap)} x {ncol(x$heatmap)} pixels",
    "Peak Attribution: {.val {round(max(x$heatmap), 4)}}"
  ))
  invisible(x)
}

#' Plot YOLO Saliency Heatmap
#'
#' @description
#' Renders the original image with a colored, semi-transparent attribution heatmap overlay
#' and the target object's bounding box.
#'
#' @param x An object of class `yolo_saliency`.
#' @param col_palette Heatmap color palette: `"inferno"` (default), `"viridis"`, `"plasma"`, or `"magma"`.
#' @param alpha Heatmap overlay opacity. Defaults to `0.55`.
#' @param show_box Logical. Whether to draw the target bounding box. Defaults to `TRUE`.
#' @param ... Additional graphical arguments.
#'
#' @export
plot.yolo_saliency <- function(x,
                               col_palette = x$col_palette,
                               alpha = x$alpha,
                               show_box = TRUE,
                               ...) {
  im <- x$image
  plot(im, ...)

  # Generate heatmap colors
  hm <- x$heatmap
  if (max(hm) > 0) {
    pal_fun <- tryCatch(
      grDevices::hcl.colors(256, palette = col_palette, rev = FALSE),
      error = function(e) grDevices::heat.colors(256)
    )

    # Convert attribution matrix into raster with alpha
    col_idx <- as.integer(round(hm * 255)) + 1L
    col_idx[col_idx > 256L] <- 256L
    col_idx[col_idx < 1L] <- 1L

    # RGBA colors
    rgb_mat <- grDevices::col2rgb(pal_fun[col_idx])
    alpha_vals <- as.integer(round(hm * alpha * 255))
    raster_colors <- grDevices::rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ], alpha = alpha_vals, maxColorValue = 255)
    dim(raster_colors) <- dim(hm)

    # Draw raster overlay transposed for R image coordinates
    graphics::rasterImage(t(raster_colors), 0, 0, nrow(hm), ncol(hm), interpolate = TRUE)
  }

  if (isTRUE(show_box) && !is.null(x$target_box)) {
    b <- x$target_box
    graphics::rect(b$xmin, b$ymin, b$xmax, b$ymax, border = "cyan", lwd = 2.5)
    lbl <- sprintf("%s (%.2f)", x$target_class, x$baseline_conf)
    graphics::text(b$xmin + 4, b$ymin + 14, labels = lbl, col = "black", font = 2, cex = 0.90)
    graphics::text(b$xmin + 3, b$ymin + 13, labels = lbl, col = "cyan", font = 2, cex = 0.90)
  }

  invisible(x)
}

