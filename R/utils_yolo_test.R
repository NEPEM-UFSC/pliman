# utils_yolo_test.R - Comprehensive YOLO Model Testing and Validation on Test Sets in pliman

#' Evaluate and Validate a Trained YOLO Model on a Test Dataset
#'
#' @description
#' `yolo_model_test()` (and its alias `yolo_test()`) evaluates an object detection or
#' instance segmentation YOLO model on an independent test dataset (or validation split).
#' It applies batch inference using [image_detect_dl_batch()] (or [image_segment_dl_batch()]),
#' matches predictions against ground-truth bounding boxes using greedy IoU bipartite matching,
#' and computes key computer vision and counting metrics:
#'
#' * **Detection Performance**: Precision, Recall, F1-Score, mAP@50, mAP@75, mAP@50:95, and mean IoU.
#' * **Counting Validation**: Ground Truth vs. Predicted counts per image, Mean Absolute Error (MAE),
#'   Root Mean Squared Error (RMSE), and coefficient of determination (\eqn{R^2}).
#' * **Confusion Matrix**: Multi-class and single-class confusion matrix, including background False Positives and False Negatives.
#' * **Diagnostic Curves**: Precision-Recall (PR) curve, F1-Confidence curve, and IoU score distribution.
#'
#' @param model A trained model identifier. Can be:
#'   * An object of class `yolo_train` (from [yolo_train()]).
#'   * Path to a trained `.onnx` or `.pt` weights file (e.g., `"best.pt"`, `"best.onnx"`).
#'   * A pre-trained model name (e.g., `"yolo11s"`, `"yolo26n"`).
#'   If a `.pt` weights file is provided, [yolo_export()] is automatically invoked if the corresponding
#'   `.onnx` file is not already found.
#' @param data Character string or `yolo_dataset` object. Path to the YOLO dataset root directory
#'   (containing `images/` and `labels/`, or `data.yaml`). Default is `"yolo_dataset"`.
#' @param split Character string. Dataset split to evaluate: `"test"` (default), `"val"`, or `"train"`.
#'   If `"test"` is selected but the directory contains no images, gracefully falls back to `"val"`.
#' @param task Character string. Computer vision task: `"detect"` (default) or `"segment"`.
#'   If `NULL`, auto-detected from the model name or dataset annotations.
#' @param conf_threshold Numeric. Confidence score threshold for candidate detections (default `0.25`).
#' @param iou_threshold Numeric. IoU overlap threshold for matching predictions to ground-truth objects (default `0.5`).
#' @param iou_range Numeric vector. Sequence of IoU thresholds used to calculate mAP50:95.
#'   Default is `seq(0.5, 0.95, by = 0.05)`.
#' @param max_det Integer. Maximum number of detections permitted per image. Defaults to `500L`
#'   to prevent truncating dense seed/capsule populations.
#' @param classes Optional character or numeric vector of class names or IDs to evaluate. If `NULL` (default),
#'   all classes present in the dataset are evaluated.
#' @param parallel Logical. Whether to run batch inference in parallel using [mirai]. Defaults to `FALSE`.
#' @param workers Optional integer. Number of parallel worker processes. Default `NULL` (auto-selected).
#' @param verbose Logical. Whether to show progress bars, timing, and formatted CLI summaries. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_detect_dl_batch()] or [image_segment_dl_batch()].
#'
#' @return An S3 object of class `c("yolo_test", "list")` containing:
#'   * `metrics`: A single-row `data.frame` with overall benchmark statistics:
#'     `total_images`, `total_gt`, `total_pred`, `tp`, `fp`, `fn`, `precision`, `recall`, `f1`,
#'     `mean_iou`, `map50`, `map75`, `map50_95`, `mae`, `rmse`, `r2`, `slope`, `intercept`.
#'   * `by_class`: A `data.frame` of per-class detection metrics (Support, Precision, Recall, F1, AP50).
#'   * `by_image`: A `data.frame` comparing ground-truth counts vs. predicted counts, differences,
#'     TP, FP, and FN per image.
#'   * `confusion_matrix`: A confusion matrix `data.frame` of counts (Ground Truth vs. Predicted, including Background).
#'   * `confusion_matrix_norm`: Normalized confusion matrix (by actual ground truth row).
#'   * `pr_curve`: A `data.frame` containing recall and precision coordinates for plotting PR curves.
#'   * `f1_curve`: A `data.frame` containing confidence thresholds and corresponding F1 scores.
#'   * `iou_distribution`: A numeric vector of exact IoU overlap values for all matched True Positives.
#'   * `predictions`: Consolidated `data.frame` of raw predictions from [image_detect_dl_batch()].
#'   * `info`: A list of execution metadata (model, split, dataset path, timestamp, task).
#'
#' @aliases yolo_test
#' @seealso [plot.yolo_test()], [yolo_train()], [image_detect_dl_batch()]
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Test model on the independent test split
#' test_res <- yolo_model_test(
#'   model = "capsulas_yolo11s/weights/best.pt",
#'   data = "dataset_capsulas",
#'   split = "test"
#' )
#'
#' # 2. View summary metrics
#' print(test_res)
#'
#' # 3. Generate diagnostic dashboard plots
#' plot(test_res, which = "all")
#' plot(test_res, which = "counts")
#' plot(test_res, which = "confusion")
yolo_model_test <- function(model,
                            data = "yolo_dataset",
                            split = c("test", "val", "train"),
                            task = NULL,
                            conf_threshold = 0.25,
                            iou_threshold = 0.5,
                            iou_range = seq(0.5, 0.95, by = 0.05),
                            max_det = 500L,
                            classes = NULL,
                            parallel = FALSE,
                            verbose = TRUE,
                            ...) {
  split <- match.arg(split)
  t0 <- Sys.time()

  # 1. Resolve Dataset Directory and Structure
  data_dir <- if (inherits(data, "yolo_dataset")) {
    data$dir
  } else if (is.character(data) && length(data) == 1L) {
    if (file.exists(data) && !dir.exists(data) && grepl("\\.ya?ml$", data, ignore.case = TRUE)) {
      dirname(normalizePath(data, winslash = "/"))
    } else {
      normalizePath(data, winslash = "/", mustWork = FALSE)
    }
  } else {
    cli::cli_abort("Invalid {.arg data} argument. Expected a directory path, a data.yaml file, or a {.code yolo_dataset} object.")
  }

  if (!dir.exists(data_dir)) {
    cli::cli_abort("Dataset directory not found at {.path {data_dir}}.")
  }

  # Auto-detect task if not provided
  if (is.null(task)) {
    task <- tryCatch(.yolo_detect_task(data_dir), error = function(e) "detect")
    if (is.character(model) && grepl("-seg", model[1], ignore.case = TRUE)) {
      task <- "segment"
    }
  }

  # Resolve candidate image folder for split
  cand_img_dirs <- c(
    file.path(data_dir, "images", split),
    file.path(data_dir, split, "images"),
    file.path(data_dir, split)
  )
  img_dir <- NULL
  for (d in cand_img_dirs) {
    if (dir.exists(d)) {
      img_dir <- d
      break
    }
  }

  # If test split was requested but directory is missing or empty, check val split
  img_exts <- "\\.(jpe?g|png|bmp|webp|tiff?)$"
  img_files <- if (!is.null(img_dir)) list.files(img_dir, pattern = img_exts, full.names = TRUE, ignore.case = TRUE) else character(0)

  if (length(img_files) == 0L && split == "test") {
    cand_val_dirs <- c(
      file.path(data_dir, "images", "val"),
      file.path(data_dir, "val", "images"),
      file.path(data_dir, "val")
    )
    for (vd in cand_val_dirs) {
      if (dir.exists(vd)) {
        val_imgs <- list.files(vd, pattern = img_exts, full.names = TRUE, ignore.case = TRUE)
        if (length(val_imgs) > 0L) {
          if (isTRUE(verbose)) {
            cli::cli_alert_warning("No images found in {.val test} split. Falling back to {.val val} split ({length(val_imgs)} images).")
          }
          img_dir <- vd
          img_files <- val_imgs
          split <- "val"
          break
        }
      }
    }
  }

  if (length(img_files) == 0L) {
    cli::cli_abort("No valid images found for split {.val {split}} in {.path {data_dir}}.")
  }

  # Resolve corresponding labels directory
  cand_lbl_dirs <- c(
    file.path(data_dir, "labels", split),
    file.path(data_dir, split, "labels")
  )
  lbl_dir <- NULL
  for (ld in cand_lbl_dirs) {
    if (dir.exists(ld)) {
      lbl_dir <- ld
      break
    }
  }

  # Class names resolution
  class_names <- tryCatch(.yolo_get_classes(data_dir), error = function(e) character(0))
  if (length(class_names) == 0L) {
    class_names <- "object"
  }

  # 2. Resolve Model Weights
  model_path <- model
  if (inherits(model, "yolo_train")) {
    if (!is.null(model$weights$best) && file.exists(model$weights$best)) {
      model_path <- model$weights$best
    } else if (!is.null(model$model)) {
      model_path <- model$model
    }
  }

  # Handle PyTorch (.pt) weights conversion to ONNX if needed
  if (is.character(model_path) && length(model_path) == 1L && grepl("\\.pt$", model_path, ignore.case = TRUE)) {
    model_norm <- normalizePath(model_path, winslash = "/", mustWork = FALSE)
    cand_onnx <- sub("\\.pt$", ".onnx", model_norm, ignore.case = TRUE)
    if (file.exists(cand_onnx)) {
      model_path <- cand_onnx
    } else {
      local_onnx <- file.path(pliman_model_dir(), paste0(tools::file_path_sans_ext(basename(model_path)), ".onnx"))
      if (file.exists(local_onnx)) {
        model_path <- local_onnx
      } else {
        if (file.exists(model_norm)) {
          if (isTRUE(verbose)) {
            cli::cli_alert_info("Exporting PyTorch weights {.file {basename(model_path)}} to standalone ONNX for high-speed inference...")
          }
          model_path <- yolo_export(model = model_norm, imgsz = 640)
        }
      }
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_rule(
      left = cli::col_cyan("pliman: YOLO Model Independent Test Evaluation"),
      right = cli::col_cyan(format(t0, "%Y-%m-%d %H:%M"))
    )
    cli::cli_bullets(c(
      "*" = sprintf("Model: {.val %s}", if (is.character(model_path)) basename(model_path) else "in-memory model"),
      "*" = sprintf("Dataset: {.path %s} | Split: {.val %s} ({length(img_files)} images)", basename(data_dir), split),
      "*" = sprintf("Task: {.val %s} | Classes: {.val %s}", task, paste(class_names, collapse = ", ")),
      "*" = sprintf("Thresholds: Conf \u2265 %.2f | IoU \u2265 %.2f | Max Detections = %d", conf_threshold, iou_threshold, max_det)
    ))
  }

  # 3. Read and Parse Ground Truth Annotations
  all_dims <- tryCatch(fast_image_dimensions_cpp(img_files), error = function(e) NULL)
  gt_list <- vector("list", length(img_files))
  names(gt_list) <- basename(img_files)

  for (i in seq_along(img_files)) {
    im_path <- img_files[i]
    stem <- tools::file_path_sans_ext(basename(im_path))

    lbl_file <- if (!is.null(lbl_dir)) file.path(lbl_dir, paste0(stem, ".txt")) else NULL
    if (is.null(lbl_file) || !file.exists(lbl_file)) {
      # Fallback check in same directory as image or alternate label path
      alt_lbl <- file.path(dirname(im_path), paste0(stem, ".txt"))
      lbl_file <- if (file.exists(alt_lbl)) alt_lbl else NULL
    }

    if (is.null(lbl_file) || !file.exists(lbl_file)) {
      # Negative background sample (0 objects)
      gt_list[[i]] <- data.frame(
        class_id = integer(0),
        class = character(0),
        xmin = numeric(0),
        ymin = numeric(0),
        xmax = numeric(0),
        ymax = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    lines <- tryCatch(readLines(lbl_file, warn = FALSE), error = function(e) character(0))
    lines <- lines[nzchar(trimws(lines))]

    if (length(lines) == 0L) {
      gt_list[[i]] <- data.frame(
        class_id = integer(0),
        class = character(0),
        xmin = numeric(0),
        ymin = numeric(0),
        xmax = numeric(0),
        ymax = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    # Read image dimension via fast C++ header parser (zero pixel loading)
    if (!is.null(all_dims) && all_dims[i, 1] > 0 && all_dims[i, 2] > 0) {
      im_w <- as.numeric(all_dims[i, 1])
      im_h <- as.numeric(all_dims[i, 2])
    } else {
      im_w <- 640
      im_h <- 640
    }

    parsed_rows <- vector("list", length(lines))
    valid_k <- 0L

    for (ln in lines) {
      parts <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
      if (length(parts) < 5L || any(is.na(parts))) next

      cid <- as.integer(parts[1])
      cname <- if (cid >= 0 && cid < length(class_names)) class_names[cid + 1L] else paste0("class_", cid)

      if (length(parts) == 5L) {
        # Standard Bounding Box: xc, yc, w, h (normalized 0-1)
        xc <- parts[2] * im_w
        yc <- parts[3] * im_h
        bw <- parts[4] * im_w
        bh <- parts[5] * im_h
        xmin <- max(0, xc - bw / 2)
        ymin <- max(0, yc - bh / 2)
        xmax <- min(im_w, xc + bw / 2)
        ymax <- min(im_h, yc + bh / 2)
      } else {
        # Polygon Segmentation Mask: x1, y1, x2, y2, ...
        coords <- parts[-1]
        xs <- coords[seq(1, length(coords), 2)] * im_w
        ys <- coords[seq(2, length(coords), 2)] * im_h
        xmin <- max(0, min(xs))
        ymin <- max(0, min(ys))
        xmax <- min(im_w, max(xs))
        ymax <- min(im_h, max(ys))
      }

      valid_k <- valid_k + 1L
      parsed_rows[[valid_k]] <- list(
        class_id = cid,
        class = cname,
        xmin = xmin,
        ymin = ymin,
        xmax = xmax,
        ymax = ymax
      )
    }

    if (valid_k > 0L) {
      gt_list[[i]] <- do.call(rbind.data.frame, c(parsed_rows[seq_len(valid_k)], list(stringsAsFactors = FALSE)))
    } else {
      gt_list[[i]] <- data.frame(
        class_id = integer(0),
        class = character(0),
        xmin = numeric(0),
        ymin = numeric(0),
        xmax = numeric(0),
        ymax = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }

  # 4. Run Batch Inference on Test Images
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Executing batch inference across {length(img_files)} images...")
  }

  detect_fn <- if (task == "segment") image_segment_dl_batch else image_detect_dl_batch

  batch_out <- detect_fn(
    data = img_files,
    model = model_path,
    labels = if (length(class_names) > 0L) class_names else NULL,
    conf_threshold = conf_threshold,
    iou_threshold = 0.45,
    parallel = parallel,
    verbose = verbose,
    max_det = max_det,
    ...
  )

  preds_df <- batch_out$predictions
  if (is.null(preds_df) || nrow(preds_df) == 0L) {
    preds_df <- data.frame(
      img = character(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      label = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Normalize predicted column names
  if (!"label" %in% names(preds_df) && "class" %in% names(preds_df)) {
    preds_df$label <- preds_df$class
  }
  if (!"score" %in% names(preds_df) && "confidence" %in% names(preds_df)) {
    preds_df$score <- preds_df$confidence
  }

  if (length(class_names) == 1L && nrow(preds_df) > 0L) {
    preds_df$label <- class_names[1]
  }

  # Filter by requested classes if specified
  if (!is.null(classes)) {
    cls_keep <- tolower(as.character(classes))
    preds_df <- preds_df[tolower(preds_df$label) %in% cls_keep, , drop = FALSE]
  }

  # 5. High-Speed Vectorized C++ Bipartite Matching and COCO Evaluation
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Computing benchmark metrics, mAP, confusion matrix, and PR curves via C++...")
  }

  img_names <- basename(img_files)
  num_images <- length(img_names)
  all_classes_unique <- unique(c(class_names, unique(preds_df$label)))
  all_classes_unique <- all_classes_unique[nzchar(all_classes_unique)]
  num_classes <- length(all_classes_unique)

  cls_map <- setNames(seq_along(all_classes_unique) - 1L, tolower(all_classes_unique))

  # Flatten GT objects for fast C++ pass
  gt_rows_list <- vector("list", num_images)
  for (i in seq_len(num_images)) {
    g_df <- gt_list[[img_names[i]]]
    if (!is.null(g_df) && nrow(g_df) > 0L) {
      c_ids <- cls_map[tolower(g_df$class)]
      c_ids[is.na(c_ids)] <- 0L
      gt_rows_list[[i]] <- data.frame(
        img_id = rep(i - 1L, nrow(g_df)),
        cls_id = as.integer(c_ids),
        x1 = as.numeric(g_df$xmin),
        y1 = as.numeric(g_df$ymin),
        x2 = as.numeric(g_df$xmax),
        y2 = as.numeric(g_df$ymax),
        stringsAsFactors = FALSE
      )
    }
  }
  flat_gt <- do.call(rbind, gt_rows_list)

  if (is.null(flat_gt) || nrow(flat_gt) == 0L) {
    gt_img_id <- integer(0)
    gt_cls_id <- integer(0)
    gt_x1 <- numeric(0); gt_y1 <- numeric(0); gt_x2 <- numeric(0); gt_y2 <- numeric(0)
  } else {
    gt_img_id <- flat_gt$img_id
    gt_cls_id <- flat_gt$cls_id
    gt_x1 <- flat_gt$x1; gt_y1 <- flat_gt$y1; gt_x2 <- flat_gt$x2; gt_y2 <- flat_gt$y2
  }

  # Prepare Predictions for C++
  if (nrow(preds_df) > 0L) {
    p_img_idx <- match(preds_df$img, img_names) - 1L
    valid_p <- which(!is.na(p_img_idx))
    pred_img_id <- p_img_idx[valid_p]
    p_cids <- cls_map[tolower(preds_df$label[valid_p])]
    p_cids[is.na(p_cids)] <- 0L
    pred_cls_id <- as.integer(p_cids)
    pred_score <- as.numeric(preds_df$score[valid_p])
    pred_x1 <- as.numeric(preds_df$xmin[valid_p])
    pred_y1 <- as.numeric(preds_df$ymin[valid_p])
    pred_x2 <- as.numeric(preds_df$xmax[valid_p])
    pred_y2 <- as.numeric(preds_df$ymax[valid_p])
  } else {
    pred_img_id <- integer(0)
    pred_cls_id <- integer(0)
    pred_score <- numeric(0)
    pred_x1 <- numeric(0); pred_y1 <- numeric(0); pred_x2 <- numeric(0); pred_y2 <- numeric(0)
  }

  eval_res <- fast_coco_eval_cpp(
    pred_img_id = pred_img_id,
    pred_cls_id = pred_cls_id,
    pred_score = pred_score,
    pred_x1 = pred_x1,
    pred_y1 = pred_y1,
    pred_x2 = pred_x2,
    pred_y2 = pred_y2,
    gt_img_id = gt_img_id,
    gt_cls_id = gt_cls_id,
    gt_x1 = gt_x1,
    gt_y1 = gt_y1,
    gt_x2 = gt_x2,
    gt_y2 = gt_y2,
    num_images = num_images,
    num_classes = num_classes,
    iou_thresholds = as.numeric(iou_range),
    eval_iou_th = as.numeric(iou_threshold)
  )

  # Build per-image statistics table
  gt_counts <- vapply(gt_list, nrow, integer(1))
  pred_counts <- if (nrow(preds_df) > 0L) as.integer(table(factor(preds_df$img, levels = img_names))) else rep(0L, num_images)
  tp_vec <- eval_res$tp_per_img
  fp_vec <- eval_res$fp_per_img
  fn_vec <- eval_res$fn_per_img
  mean_ious <- eval_res$mean_iou_per_img

  diff_vec <- pred_counts - gt_counts
  prec_vec <- ifelse((tp_vec + fp_vec) > 0, tp_vec / (tp_vec + fp_vec), NA_real_)
  rec_vec <- ifelse(gt_counts > 0, tp_vec / gt_counts, NA_real_)
  f1_vec <- ifelse(!is.na(prec_vec) & !is.na(rec_vec) & (prec_vec + rec_vec) > 0,
                   2 * prec_vec * rec_vec / (prec_vec + rec_vec), NA_real_)

  by_image_df <- data.frame(
    image = img_names,
    gt_count = gt_counts,
    pred_count = pred_counts,
    diff = diff_vec,
    abs_diff = abs(diff_vec),
    tp = tp_vec,
    fp = fp_vec,
    fn = fn_vec,
    precision = round(prec_vec, 4),
    recall = round(rec_vec, 4),
    f1 = round(f1_vec, 4),
    mean_iou = ifelse(mean_ious > 0, round(mean_ious, 4), NA_real_),
    stringsAsFactors = FALSE
  )

  # Global summary metrics
  tot_tp <- sum(tp_vec)
  tot_fp <- sum(fp_vec)
  tot_fn <- sum(fn_vec)
  tot_gt <- sum(gt_counts)
  tot_pred <- sum(pred_counts)

  glob_prec <- if ((tot_tp + tot_fp) > 0) tot_tp / (tot_tp + tot_fp) else 0
  glob_rec <- if (tot_gt > 0) tot_tp / tot_gt else 0
  glob_f1 <- if ((glob_prec + glob_rec) > 0) 2 * (glob_prec * glob_rec) / (glob_prec + glob_rec) else 0
  all_matched_ious <- eval_res$matched_ious
  mean_tp_iou <- if (length(all_matched_ious) > 0) mean(all_matched_ious) else 0

  # Counting regression statistics
  mae_val <- mean(abs(diff_vec))
  rmse_val <- sqrt(mean(diff_vec^2))
  r2_val <- if (stats::sd(gt_counts) > 0 && stats::sd(pred_counts) > 0) {
    stats::cor(gt_counts, pred_counts)^2
  } else {
    1.0
  }
  fit_lm <- stats::lm(pred_counts ~ gt_counts)
  slope_val <- as.numeric(stats::coef(fit_lm)[2])
  intercept_val <- as.numeric(stats::coef(fit_lm)[1])

  map50_val <- eval_res$map50
  map75_val <- eval_res$map75
  map50_95_val <- eval_res$map50_95

  # PR Curve & F1 Curve
  pr_curve_df <- data.frame(
    recall = eval_res$pr_recall,
    precision = eval_res$pr_precision,
    confidence = eval_res$pr_confidence,
    stringsAsFactors = FALSE
  )

  f1_curve_df <- data.frame(
    confidence = eval_res$f1_confidence,
    f1 = eval_res$f1_values,
    stringsAsFactors = FALSE
  )

  # Confusion matrix
  cm_labels <- c(all_classes_unique, "background")
  cm_mat <- eval_res$confusion_matrix
  dimnames(cm_mat) <- list(Actual = cm_labels, Predicted = cm_labels)
  cm_df <- as.data.frame(cm_mat)
  cm_row_sums <- rowSums(cm_mat)
  cm_norm_mat <- cm_mat / pmax(1, cm_row_sums)
  cm_norm_df <- as.data.frame(round(cm_norm_mat, 4))

  # Per-Class Breakdown
  by_class_list <- vector("list", num_classes)
  for (k in seq_len(num_classes)) {
    cname <- all_classes_unique[k]
    gt_c <- if (nrow(flat_gt) > 0L) sum(flat_gt$cls_id == (k - 1L)) else 0L
    pred_c <- if (nrow(preds_df) > 0L) sum(pred_cls_id == (k - 1L)) else 0L
    tp_c <- cm_mat[cname, cname]
    fp_c <- sum(cm_mat[, cname]) - tp_c
    fn_c <- sum(cm_mat[cname, ]) - tp_c

    prec_c <- if ((tp_c + fp_c) > 0) tp_c / (tp_c + fp_c) else 0
    rec_c <- if (gt_c > 0) tp_c / gt_c else 0
    f1_c <- if ((prec_c + rec_c) > 0) 2 * (prec_c * rec_c) / (prec_c + rec_c) else 0

    by_class_list[[k]] <- list(
      class = cname,
      support = gt_c,
      pred = pred_c,
      tp = tp_c,
      fp = fp_c,
      fn = fn_c,
      precision = round(prec_c, 4),
      recall = round(rec_c, 4),
      f1 = round(f1_c, 4)
    )
  }
  by_class_df <- do.call(rbind.data.frame, c(by_class_list, list(stringsAsFactors = FALSE)))

  # Consolidated single-row metrics data frame
  metrics_df <- data.frame(
    split = split,
    images = length(img_files),
    total_gt = tot_gt,
    total_pred = tot_pred,
    tp = tot_tp,
    fp = tot_fp,
    fn = tot_fn,
    precision = round(glob_prec, 4),
    recall = round(glob_rec, 4),
    f1 = round(glob_f1, 4),
    mean_iou = round(mean_tp_iou, 4),
    map50 = round(map50_val, 4),
    map75 = round(map75_val, 4),
    map50_95 = round(map50_95_val, 4),
    mae = round(mae_val, 3),
    rmse = round(rmse_val, 3),
    r2 = round(r2_val, 4),
    slope = round(slope_val, 3),
    intercept = round(intercept_val, 3),
    stringsAsFactors = FALSE
  )

  res <- structure(
    list(
      metrics = metrics_df,
      by_class = by_class_df,
      by_image = by_image_df,
      confusion_matrix = cm_df,
      confusion_matrix_norm = cm_norm_df,
      pr_curve = pr_curve_df,
      f1_curve = f1_curve_df,
      iou_distribution = all_matched_ious,
      predictions = preds_df,
      info = list(
        model = if (is.character(model_path)) basename(model_path) else "in-memory",
        data = data_dir,
        split = split,
        task = task,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        max_det = max_det,
        n_images = length(img_files),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    ),
    class = c("yolo_test", "list")
  )

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Testing completed in {.val {round(difftime(Sys.time(), t0, units = 'secs'), 1)}s}!")
    print(res)
  }

  invisible(res)
}

#' @rdname yolo_model_test
#' @export
yolo_test <- yolo_model_test


# ==============================================================================
# S3 METHODS: print, summary, and plot for yolo_test
# ==============================================================================

#' @export
#' @export
print.yolo_test <- function(x, digits = 4, class_details = TRUE, ...) {
  m <- x$metrics
  inf <- x$info
  fmt_d <- paste0("%.", digits, "f")
  fmt_pct <- paste0("%.", max(1, digits - 2), "f%%")

  cli::cli_h1("YOLO Test Evaluation Benchmark")
  cli::cli_bullets(c(
    "*" = sprintf("Model: {.val %s} | Split: {.val %s} [{toupper(inf$task)}]", inf$model, inf$split),
    "*" = sprintf("Dataset Scale: %s test images | Ground Truth: %s objects | Predictions: %s objects",
                  format(m$images, big.mark = ","),
                  format(m$total_gt, big.mark = ","),
                  format(m$total_pred, big.mark = ",")),
    "*" = sprintf("Operating Thresholds: Conf \u2265 %.2f | IoU Match \u2265 %.2f",
                  inf$conf_threshold, inf$iou_threshold)
  ))

  cli::cli_h2("1. Computer Vision Detection Accuracy")
  cli::cli_bullets(c(
    "v" = sprintf(paste0("Precision (P)   : {.strong ", fmt_pct, "} (TP = %s | FP = %s)"),
                  m$precision * 100, format(m$tp, big.mark = ","), format(m$fp, big.mark = ",")),
    "v" = sprintf(paste0("Recall (R)      : {.strong ", fmt_pct, "} (Found %s of %s ground truth instances; FN = %s)"),
                  m$recall * 100, format(m$tp, big.mark = ","), format(m$total_gt, big.mark = ","), format(m$fn, big.mark = ",")),
    "v" = sprintf(paste0("F1-Score        : {.strong ", fmt_d, "}"), m$f1),
    "v" = sprintf(paste0("mAP @ 0.50      : {.strong ", fmt_d, "} (", fmt_pct, ")"), m$map50, m$map50 * 100),
    "v" = sprintf(paste0("mAP @ 0.75      : {.strong ", fmt_d, "} (", fmt_pct, ")"), m$map75, m$map75 * 100),
    "v" = sprintf(paste0("mAP @ 0.50:0.95 : {.strong ", fmt_d, "} (", fmt_pct, ")"), m$map50_95, m$map50_95 * 100),
    "v" = sprintf(paste0("Mean Match IoU  : {.strong ", fmt_d, "} (Quality of overlap for True Positives)"), m$mean_iou)
  ))

  cli::cli_h2("2. Seed & Object Counting Fidelity")
  cli::cli_bullets(c(
    "v" = sprintf(paste0("Coefficient of Determination (R\u00b2) : {.strong ", fmt_d, "}"), m$r2),
    "v" = sprintf(paste0("Mean Absolute Error (MAE)           : {.strong ", fmt_d, "} objects/image"), m$mae),
    "v" = sprintf(paste0("Root Mean Squared Error (RMSE)       : {.strong ", fmt_d, "} objects/image"), m$rmse),
    "v" = sprintf(paste0("Calibration Line                     : Pred = ", fmt_d, " %s ", fmt_d, " \u00d7 GT"),
                  m$intercept, if (m$intercept >= 0) "+" else "-", abs(m$slope))
  ))

  if (isTRUE(class_details) && !is.null(x$by_class) && nrow(x$by_class) > 1L) {
    cli::cli_h2("3. Per-Class Performance Breakdown")
    print(x$by_class, row.names = FALSE)
  }

  cli::cli_alert_info("Run {.run plot(x, which = 'all')} to visualize Confusion Matrix, Regression Line, PR Curve, and IoU distribution.")
  invisible(x)
}

#' @export
summary.yolo_test <- function(object, digits = 4, n_images = 10, ...) {
  m_round <- object$metrics
  num_cols_m <- vapply(m_round, is.numeric, logical(1))
  m_round[num_cols_m] <- lapply(m_round[num_cols_m], round, digits = digits)

  cls_round <- object$by_class
  if (!is.null(cls_round)) {
    num_cols_c <- vapply(cls_round, is.numeric, logical(1))
    cls_round[num_cols_c] <- lapply(cls_round[num_cols_c], round, digits = digits)
  }

  img_df <- object$by_image
  if (!is.null(img_df)) {
    num_cols_i <- vapply(img_df, is.numeric, logical(1))
    img_df[num_cols_i] <- lapply(img_df[num_cols_i], round, digits = digits)
    if (is.finite(n_images) && n_images > 0) {
      img_df <- utils::head(img_df, n_images)
    }
  }

  structure(
    list(
      metrics = m_round,
      by_class = cls_round,
      by_image = img_df,
      confusion_matrix = object$confusion_matrix,
      confusion_matrix_norm = object$confusion_matrix_norm,
      info = object$info,
      digits = digits
    ),
    class = "summary.yolo_test"
  )
}

#' @export
print.summary.yolo_test <- function(x, ...) {
  cli::cli_h1("Summary: YOLO Test Evaluation")
  cli::cli_text("Model: {.strong {x$info$model}} | Task: {.strong {toupper(x$info$task)}} | Images: {.strong {x$info$n_images}}")
  cli::cli_h2("Global Metrics")
  print(x$metrics, row.names = FALSE)
  if (!is.null(x$by_class) && nrow(x$by_class) > 1L) {
    cli::cli_h2("Per-Class Metrics")
    print(x$by_class, row.names = FALSE)
  }
  if (!is.null(x$by_image) && nrow(x$by_image) > 0L) {
    cli::cli_h2(sprintf("Sample Image Predictions (First %d rows)", nrow(x$by_image)))
    print(x$by_image, row.names = FALSE)
  }
  invisible(x)
}

#' Plot Diagnostic and Validation Graphs for YOLO Test Results
#'
#' @description
#' Generates publication-ready diagnostic plots for an evaluated YOLO model on the test dataset:
#' * `which = "all"`: A 4-panel dashboard containing Confusion Matrix, Count Validation Scatter,
#'   Precision-Recall Curve, and IoU Distribution.
#' * `which = "confusion"`: Formatted heat map of the Confusion Matrix (normalized or counts).
#' * `which = "counts"`: Ground Truth vs. Predicted Count regression scatter plot with 1:1 diagonal and \eqn{R^2}.
#' * `which = "pr"`: Precision-Recall curve with shaded area under curve.
#' * `which = "iou"`: Histogram and density distribution of True Positive IoU overlap values.
#' * `which = "f1"`: F1-score across confidence thresholds curve.
#'
#' You can also pass multiple choices to `which` (e.g. `which = c("confusion", "counts")`).
#'
#' @param x An object of class `yolo_test` (produced by [yolo_model_test()]).
#' @param which Character vector specifying which plot(s) to render: `"all"` (default),
#'   or any combination of `"confusion"`, `"counts"`, `"pr"`, `"iou"`, and `"f1"`.
#' @param ncol Optional integer. Number of columns in the plot grid. If `NULL` (default),
#'   automatically computed based on the number of plots.
#' @param nrow Optional integer. Number of rows in the plot grid. If `NULL` (default),
#'   automatically computed based on the number of plots and `ncol`.
#' @param digits Integer. Number of decimal digits for displaying statistics in plot titles,
#'   legends, regression equations, and confusion matrix cells (default `3`).
#' @param normalized Logical. If `TRUE` (default), plots normalized confusion matrix frequencies (0 to 1).
#' @param bg Character string. Background color for plots (default `"#fafafa"`).
#' @param col_cm Character vector of colors for the diagonal of the confusion matrix (correct detections).
#'   Defaults to an emerald green palette.
#' @param col_cm_error Character vector of colors for off-diagonal cells (misclassifications/errors).
#'   Defaults to a rose red palette.
#' @param col_points Character string. Fill color for scatter points in count validation (default `"#4f46e5"`).
#' @param col_line Character string. Color for linear regression line (default `"#2563eb"`).
#' @param col_ideal Character string. Color for 1:1 identity dashed line (default `"#94a3b8"`).
#' @param col_pr Character string. Color for Precision-Recall curve (default `"#059669"`).
#' @param col_iou Character string. Color for IoU distribution histogram bars (default `"#3b82f6"`).
#' @param col_f1 Character string. Color for F1-Confidence curve (default `"#7c3aed"`).
#' @param grid Logical. Whether to display subtle background grid lines (default `TRUE`).
#' @param col_grid Character string. Grid line color (default `"#e2e8f0"`).
#' @param mar Numeric vector of length 4. Margins for subplots (default `c(4.2, 4.5, 3.2, 1.5)`).
#' @param mgp Numeric vector of length 3. Margin line for axis title, labels and line (default `c(2.2, 0.7, 0)`).
#' @param cex.main Numeric. Magnification for main titles (default `1.15`).
#' @param cex.lab Numeric. Magnification for x and y labels (default `1.0`).
#' @param cex.axis Numeric. Magnification for axis tick labels (default `0.85`).
#' @param font.main Integer. Font style for titles (default `2`, bold).
#' @param save Logical. If `TRUE`, saves the generated plot to an image or document file (default `FALSE`).
#' @param file_name Optional character string. Custom file name or path to save the plot (e.g. `"test_dashboard.png"`, `"fig_counts.pdf"`).
#'   If `NULL` and `save = TRUE`, a default descriptive name is generated based on the model and selected plots.
#' @param width Numeric. Width of the output plot in inches (or pixels if > 50). If `NULL`,
#'   automatically computed based on the number of subplot columns.
#' @param height Numeric. Height of the output plot in inches (or pixels if > 50). If `NULL`,
#'   automatically computed based on the number of subplot rows.
#' @param dpi Numeric. Resolution in dots per inch (DPI) for bitmap devices (default `300`).
#' @param ... Additional arguments passed to base graphics functions.
#'
#' @return Invisibly returns the `yolo_test` object.
#' @export
plot.yolo_test <- function(x,
                           which = c("all", "confusion", "counts", "pr", "iou", "f1"),
                           ncol = NULL,
                           nrow = NULL,
                           digits = 3,
                           normalized = TRUE,
                           bg = "white",
                           col_cm = c("#f0fdf4", "#86efac", "#22c55e", "#15803d"),
                           col_cm_error = c("#fff1f2", "#fecdd3", "#f43f5e", "#be123c"),
                           col_points = "#4f46e5",
                           col_line = "#2563eb",
                           col_ideal = "#94a3b8",
                           col_pr = "#059669",
                           col_iou = "#3b82f6",
                           col_f1 = "#7c3aed",
                           grid = TRUE,
                           col_grid = "#e2e8f0",
                           mar = c(4.2, 5.2, 3.2, 1.5),
                           mgp = c(2.2, 0.7, 0),
                           cex.main = 1.15,
                           cex.lab = 1.0,
                           cex.axis = 0.85,
                           font.main = 2,
                           save = FALSE,
                           file_name = NULL,
                           width = NULL,
                           height = NULL,
                           dpi = 300,
                           ...) {
  if (!inherits(x, "yolo_test")) {
    cli::cli_abort("Object must be of class {.code yolo_test} (from {.fn yolo_model_test}).")
  }

  # Normalize aliases and synonyms for which
  norm_which <- function(w) {
    w <- tolower(trimws(w))
    w <- sub("[-_ ]+", "_", w)
    vapply(w, function(item) {
      if (item %in% c("all", "dashboard", "full")) return("all")
      if (item %in% c("counts", "counting", "count", "regression", "scatter", "correlation")) return("counts")
      if (item %in% c("confusion", "confusion_matrix", "cm", "matrix", "confusion_matrix_norm")) return("confusion")
      if (item %in% c("pr", "pr_curve", "precision_recall", "auc")) return("pr")
      if (item %in% c("iou", "ious", "iou_distribution", "hist", "histogram", "distribution")) return("iou")
      if (item %in% c("f1", "f1_curve", "f1_score", "f1_confidence", "confidence")) return("f1")
      item
    }, character(1), USE.NAMES = FALSE)
  }

  which_norm <- norm_which(which)
  valid_plots <- c("confusion", "counts", "pr", "iou", "f1")
  if ("all" %in% which_norm) {
    selected_plots <- c("confusion", "counts", "pr", "iou")
  } else {
    selected_plots <- match.arg(which_norm, valid_plots, several.ok = TRUE)
  }

  n_plots <- length(selected_plots)
  if (n_plots == 0L) {
    cli::cli_abort("No valid plots selected in {.arg which}.")
  }

  # Layout calculation
  if (n_plots == 1L) {
    num_col <- if (!is.null(ncol)) ncol else 1L
    num_row <- if (!is.null(nrow)) nrow else 1L
  } else {
    if (is.null(ncol) && is.null(nrow)) {
      if (n_plots == 2L) {
        num_col <- 2L; num_row <- 1L
      } else if (n_plots <= 4L) {
        num_col <- 2L; num_row <- 2L
      } else {
        num_col <- 3L; num_row <- ceiling(n_plots / 3L)
      }
    } else if (!is.null(ncol) && is.null(nrow)) {
      num_col <- as.integer(ncol)
      num_row <- ceiling(n_plots / num_col)
    } else if (is.null(ncol) && !is.null(nrow)) {
      num_row <- as.integer(nrow)
      num_col <- ceiling(n_plots / num_row)
    } else {
      num_col <- as.integer(ncol)
      num_row <- as.integer(nrow)
    }
  }

  # Setup saving device if requested
  model_tag <- if (!is.null(x$info$model) && nzchar(x$info$model) && x$info$model != "in-memory") {
    tools::file_path_sans_ext(basename(x$info$model))
  } else if (!is.null(x$info$data) && nzchar(x$info$data)) {
    basename(x$info$data)
  } else {
    "yolo_test"
  }

  tag_which <- if (identical(selected_plots, c("confusion", "counts", "pr", "iou")) || "all" %in% which_norm) {
    "dashboard"
  } else {
    paste(selected_plots, collapse = "_")
  }
  default_file <- paste0(model_tag, "_test_", tag_which, ".png")

  def_w <- if (num_col == 1L) 7 else if (num_col == 2L) 10 else max(9, num_col * 4.2)
  def_h <- if (num_row == 1L) 6 else if (num_row == 2L) 8.5 else max(7, num_row * 3.8)

  dev_info <- .open_plot_device(
    save = save,
    file_name = file_name,
    width = width,
    height = height,
    dpi = dpi,
    default_name = default_file,
    default_w = def_w,
    default_h = def_h
  )
  if (!is.null(dev_info)) {
    on.exit(.close_plot_device(dev_info), add = TRUE)
  }

  if (is.null(dev_info)) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
  }
  graphics::par(mfrow = c(num_row, num_col), mar = mar, mgp = mgp, bg = bg)

  add_grid <- function() {
    if (isTRUE(grid)) {
      graphics::grid(nx = NULL, ny = NULL, col = col_grid, lty = 1, lwd = 0.8)
    }
  }

  fmt_d <- paste0("%.", digits, "f")

  # --- PLOT 1: Confusion Matrix ---
  plot_cm <- function() {
    cm <- if (isTRUE(normalized)) x$confusion_matrix_norm else x$confusion_matrix
    cm_cnt <- x$confusion_matrix
    if (is.null(cm) || nrow(cm) == 0L) return(invisible(NULL))

    nr <- nrow(cm)
    nc <- ncol(cm)
    row_labels <- rownames(cm)
    col_labels <- colnames(cm)

    cex_ax <- if (nc > 4) cex.axis * 0.85 else cex.axis

    # Dynamically measure width of row labels in margin lines
    # so 'Actual Ground Truth' never collides with labels or gets clipped
    lab_w_in <- max(graphics::strwidth(row_labels, units = "inches", cex = cex_ax))
    csi_in <- max(0.01, graphics::par("csi"))
    lab_lines <- lab_w_in / csi_in

    # Calculate needed left margin and exact ylab position
    ylab_line <- lab_lines + 1.0 * cex.lab + 0.2
    needed_left_mar <- max(mar[2], ylab_line + 1.2 * cex.lab + 0.5)

    cur_mar <- graphics::par("mar")
    if (cur_mar[2] < needed_left_mar) {
      graphics::par(mar = c(cur_mar[1], needed_left_mar, cur_mar[3], cur_mar[4]))
    }

    greens <- grDevices::colorRampPalette(col_cm)(100)
    reds   <- grDevices::colorRampPalette(col_cm_error)(100)

    graphics::plot(1, type = "n", xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
                   axes = FALSE, xlab = "Predicted Class", ylab = "",
                   main = if (isTRUE(normalized)) "Confusion Matrix (Normalized)" else "Confusion Matrix (Counts)",
                   font.main = font.main, cex.main = cex.main, cex.lab = cex.lab, ...)
    graphics::title(ylab = "Actual Ground Truth", line = ylab_line, font.lab = font.main, cex.lab = cex.lab)

    mat_val <- as.matrix(cm)
    mat_cnt <- as.matrix(cm_cnt)

    for (r in seq_len(nr)) {
      for (c in seq_len(nc)) {
        val <- mat_val[r, c]
        cnt <- mat_cnt[r, c]
        y_pos <- nr - r + 1
        x_pos <- c
        is_diag <- (r == c)

        if (is_diag) {
          int_val <- if (is.na(val) || val <= 0) 0 else max(0.15, min(1, val))
          idx <- max(1, min(100, round(int_val * 99) + 1))
          cell_col <- greens[idx]
          txt_col <- if (int_val > 0.6) "#ffffff" else "#14532d"
        } else {
          int_val <- if (is.na(val) || val <= 0) 0 else max(0.12, min(1, val * 1.5))
          idx <- max(1, min(100, round(int_val * 99) + 1))
          cell_col <- if (cnt == 0) "#ffffff" else reds[idx]
          txt_col <- if (int_val > 0.6) "#ffffff" else "#881337"
        }

        graphics::rect(x_pos - 0.5, y_pos - 0.5, x_pos + 0.5, y_pos + 0.5,
                       col = cell_col, border = "#cbd5e1", lwd = 0.8)

        lbl_text <- if (isTRUE(normalized)) {
          sprintf(paste0(fmt_d, "\n(%d)"), val, cnt)
        } else {
          sprintf(paste0("%d\n(", fmt_d, "%%)"), cnt, val * 100)
        }
        graphics::text(x_pos, y_pos, labels = lbl_text, col = txt_col, font = 2, cex = 0.82)
      }
    }
    cex_ax <- if (nc > 4) cex.axis * 0.85 else cex.axis
    graphics::axis(1, at = seq_len(nc), labels = col_labels, font = 2, tick = FALSE, line = -0.5, cex.axis = cex_ax)
    graphics::axis(2, at = rev(seq_len(nr)), labels = row_labels, font = 2, tick = FALSE, line = -0.5, las = 2, cex.axis = cex_ax)
  }

  # --- PLOT 2: Ground Truth vs Predicted Scatter Plot ---
  plot_counts <- function() {
    df <- x$by_image
    m <- x$metrics
    if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))

    max_v <- max(c(df$gt_count, df$pred_count), na.rm = TRUE) * 1.08
    if (max_v <= 0) max_v <- 10

    main_title <- sprintf(paste0("Counting Validation (R\u00b2 = ", fmt_d, " | MAE = ", fmt_d, ")"), m$r2, m$mae)

    graphics::plot(df$gt_count, df$pred_count, type = "n",
                   xlim = c(0, max_v), ylim = c(0, max_v),
                   xlab = "Ground Truth Count", ylab = "Predicted Count",
                   main = main_title, font.main = font.main, cex.main = cex.main,
                   cex.lab = cex.lab, cex.axis = cex.axis, ...)
    add_grid()

    # 1:1 Identity dashed line (Ideal)
    graphics::abline(a = 0, b = 1, col = col_ideal, lty = 2, lwd = 1.8)

    # Linear regression fit line
    graphics::abline(a = m$intercept, b = m$slope, col = col_line, lwd = 2.4)

    # Semi-transparent point circles
    graphics::points(df$gt_count, df$pred_count, pch = 21,
                     bg = grDevices::adjustcolor(col_points, alpha.f = 0.65),
                     col = "#ffffff", cex = 1.35, lwd = 0.8)

    # Legend & formula badge
    sign_str <- if (m$intercept >= 0) "+" else "-"
    eq_lbl <- sprintf(paste0("Fit: y = ", fmt_d, " %s ", fmt_d, "x\nR\u00b2 = ", fmt_d, " | RMSE = ", fmt_d),
                      m$intercept, if (m$intercept >= 0) "+" else "-", abs(m$slope), m$r2, m$rmse)
    graphics::legend("topleft", legend = c("1:1 Perfect Line", eq_lbl),
                     col = c(col_ideal, col_line), lty = c(2, 1), lwd = c(1.8, 2.4),
                     bty = "n", cex = 0.82)
  }

  # --- PLOT 3: Precision-Recall Curve ---
  plot_pr <- function() {
    pr <- x$pr_curve
    m <- x$metrics
    if (is.null(pr) || nrow(pr) == 0L) return(invisible(NULL))

    main_title <- sprintf(paste0("Precision-Recall Curve (mAP50 = ", fmt_d, ")"), m$map50)

    graphics::plot(pr$recall, pr$precision, type = "n",
                   xlim = c(0, 1), ylim = c(0, 1.05),
                   xlab = "Recall", ylab = "Precision",
                   main = main_title, font.main = font.main, cex.main = cex.main,
                   cex.lab = cex.lab, cex.axis = cex.axis, ...)
    add_grid()

    graphics::polygon(c(pr$recall, 1, 0), c(pr$precision, 0, 0),
                      col = grDevices::adjustcolor(col_pr, alpha.f = 0.18), border = NA)
    graphics::lines(pr$recall, pr$precision, col = col_pr, lwd = 2.6)

    leg_text <- sprintf(paste0("mAP@0.50: ", fmt_d, "\nmAP@0.50:0.95: ", fmt_d), m$map50, m$map50_95)
    graphics::legend("bottomleft", legend = leg_text, col = col_pr, lwd = 2.6, bty = "n", cex = 0.85)
  }

  # --- PLOT 4: IoU Overlap Quality Distribution ---
  plot_iou <- function() {
    ious <- x$iou_distribution
    m <- x$metrics
    if (is.null(ious) || length(ious) == 0L) return(invisible(NULL))

    clamped_ious <- pmax(0.40, pmin(1.0, ious))
    h <- graphics::hist(clamped_ious, breaks = seq(0.4, 1.0, by = 0.05), plot = FALSE)
    mean_iou_val <- if (!is.null(m$mean_iou) && !is.na(m$mean_iou)) m$mean_iou else mean(ious, na.rm = TRUE)
    main_title <- sprintf(paste0("Match Quality (Mean IoU = ", fmt_d, ")"), mean_iou_val)

    graphics::plot(h, col = grDevices::adjustcolor(col_iou, alpha.f = 0.6),
                   border = grDevices::adjustcolor(col_iou, alpha.f = 1.0),
                   xlab = "Intersection over Union (IoU)",
                   main = main_title, font.main = font.main, cex.main = cex.main,
                   cex.lab = cex.lab, cex.axis = cex.axis, ...)
    add_grid()
    graphics::plot(h, col = grDevices::adjustcolor(col_iou, alpha.f = 0.6),
                   border = grDevices::adjustcolor(col_iou, alpha.f = 1.0), add = TRUE)
    graphics::abline(v = mean_iou_val, col = "#dc2626", lwd = 2.2, lty = 2)
    lbl_pos <- if (isTRUE(mean_iou_val > 0.82)) 2 else 4
    lbl_str <- if (lbl_pos == 2) sprintf(paste0("Mean: ", fmt_d, " "), mean_iou_val) else sprintf(paste0(" Mean: ", fmt_d), mean_iou_val)
    graphics::text(mean_iou_val, max(h$counts) * 0.9,
                   labels = lbl_str,
                   pos = lbl_pos, col = "#dc2626", font = 2, cex = 0.85)
  }

  # --- PLOT 5: F1 vs Confidence Curve ---
  plot_f1 <- function() {
    f1_df <- x$f1_curve
    m <- x$metrics
    if (is.null(f1_df) || nrow(f1_df) == 0L) return(invisible(NULL))

    best_idx <- which.max(f1_df$f1)
    best_conf <- f1_df$confidence[best_idx]
    best_f1 <- f1_df$f1[best_idx]

    main_title <- sprintf(paste0("F1-Confidence Curve (Max F1 = ", fmt_d, " at \u03c4 = %.2f)"), best_f1, best_conf)

    graphics::plot(f1_df$confidence, f1_df$f1, type = "n",
                   xlim = c(0, 1), ylim = c(0, 1.05),
                   xlab = "Confidence Threshold", ylab = "F1-Score",
                   main = main_title, font.main = font.main, cex.main = cex.main,
                   cex.lab = cex.lab, cex.axis = cex.axis, ...)
    add_grid()

    graphics::polygon(c(f1_df$confidence, 1, 0), c(f1_df$f1, 0, 0),
                      col = grDevices::adjustcolor(col_f1, alpha.f = 0.18), border = NA)
    graphics::lines(f1_df$confidence, f1_df$f1, col = col_f1, lwd = 2.6)
    graphics::points(best_conf, best_f1, pch = 21, bg = col_f1, col = "#ffffff", cex = 1.5)
    graphics::text(best_conf, best_f1, labels = sprintf(paste0("  Max F1: ", fmt_d), best_f1),
                   pos = 4, font = 2, col = col_f1, cex = 0.85)
  }

  # --- Dispatch loop for selected plots ---
  for (plt in selected_plots) {
    if (plt == "confusion") plot_cm()
    else if (plt == "counts") plot_counts()
    else if (plt == "pr") plot_pr()
    else if (plt == "iou") plot_iou()
    else if (plt == "f1") plot_f1()
  }

  invisible(x)
}
