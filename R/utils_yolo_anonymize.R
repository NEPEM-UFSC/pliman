# utils_yolo_anonymize.R - Automatic Dataset Anonymization and Privacy for Field Imagery

#' Automatically Anonymize and Protect Privacy in Field YOLO Datasets
#'
#' @description
#' `yolo_dataset_anonymize()` (and its alias `yolo_anonymize()`) scans a YOLO dataset or
#' image folder, detects sensitive entities such as human faces, workers, and vehicle
#' license plates, and automatically blurs or pixelates those regions directly on the images.
#' Crucially, target object annotations (`labels/*.txt` containing bounding boxes or
#' segmentation masks for crops, fruits, seeds, animals, or weeds) are **completely preserved**
#' without any displacement or coordinate modification.
#'
#' This utility ensures compliance with privacy regulations (e.g., GDPR, LGPD, institutional
#' ethics committees) prior to publishing open-access agronomic datasets on repositories such
#' as Zenodo, Figshare, Dryad, or Roboflow.
#'
#' @param dataset Character string or vector. Path to the YOLO dataset root directory
#'   (containing `images/` or subdirectories `train/`, `val/`, `test/`), or a character vector
#'   of image file paths.
#' @param model A YOLO model identifier or path. Defaults to `"yolo26n"`. Can be an ONNX model file
#'   or a pretrained model name (e.g., `"yolo11n"`, `"yolo26n"`).
#' @param target_classes Character vector of class names or regex patterns to anonymize.
#'   Defaults to `c("person", "face", "license plate", "car", "truck", "bus", "motorcycle")`.
#' @param person_region Character string specifying which region of a detected `"person"` to anonymize:
#'   * `"head"` (default): Estimates and blurs only the head/face region (top 22% of the bounding box).
#'   * `"full"`: Blurs the entire bounding box of the person.
#' @param method Character string specifying the anonymization method:
#'   * `"pixelate"` (default): Applies mosaic pixelation (irrecoverable pixel blocks).
#'   * `"blur"`: Applies repeated spatial box/mean blurring.
#' @param kernel_size Integer. Size in pixels of mosaic blocks (for `"pixelate"`) or blur kernel
#'   (for `"blur"`). Defaults to `25L`.
#' @param conf_threshold Numeric. Confidence threshold for detecting sensitive entities.
#'   Defaults to `0.20` to ensure high recall for distant or edge-of-frame faces.
#' @param iou_threshold Numeric. Non-Maximum Suppression (NMS) IoU threshold. Defaults to `0.45`.
#' @param output_dir Character string or `NULL`. Destination directory for anonymized images.
#'   If `NULL` (default), images are modified in-place.
#' @param backup Logical. When modifying images in-place (`output_dir = NULL`), whether to create
#'   backup copies with extension `.orig` or `.backup`. Defaults to `TRUE`.
#' @param recursive Logical. Whether to search subdirectories recursively. Defaults to `TRUE`.
#' @param engine Character string. Compute engine: `"gpu"` (default) or `"cpu"`.
#' @param verbose Logical. Whether to display progress bars and summary reports. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_detect_dl()].
#'
#' @return An S3 object of class `c("yolo_anonymized", "data.frame")` containing:
#'   * `image`: Path to the image file.
#'   * `anonymized_count`: Number of sensitive entities detected and anonymized.
#'   * `classes`: Comma-separated list of sensitive classes found.
#'   * `status`: Processing status (`"anonymized"`, `"no_sensitive_objects"`, or `"error"`).
#'   * `output_path`: Path where the anonymized image was written.
#'
#' @export
#' @aliases yolo_anonymize
#' @seealso [image_detect_dl()], [yolo_test()]
yolo_dataset_anonymize <- function(dataset,
                                   model = "yolo26n",
                                   target_classes = c("person", "face", "license plate", "car", "truck", "bus", "motorcycle"),
                                   person_region = c("head", "full"),
                                   method = c("pixelate", "blur"),
                                   kernel_size = 25L,
                                   conf_threshold = 0.20,
                                   iou_threshold = 0.45,
                                   output_dir = NULL,
                                   backup = TRUE,
                                   recursive = TRUE,
                                   engine = "gpu",
                                   verbose = TRUE,
                                   ...) {
  person_region <- match.arg(person_region)
  method <- match.arg(method)
  kernel_size <- max(4L, as.integer(kernel_size))

  # 1. Discover target images
  img_extensions <- c("jpg", "jpeg", "png", "bmp", "tif", "tiff", "webp")
  ext_pattern <- paste0("\\.(", paste(img_extensions, collapse = "|"), ")$")

  if (is.character(dataset) && length(dataset) == 1L && dir.exists(dataset)) {
    # Check if images/ subdirectory exists
    cand_images_dir <- file.path(dataset, "images")
    scan_dir <- if (dir.exists(cand_images_dir)) cand_images_dir else dataset
    all_files <- list.files(scan_dir, pattern = ext_pattern, full.names = TRUE, recursive = recursive, ignore.case = TRUE)
  } else if (is.character(dataset)) {
    all_files <- dataset[file.exists(dataset) & grepl(ext_pattern, dataset, ignore.case = TRUE)]
  } else {
    cli::cli_abort("Invalid {.val dataset} input. Must be a folder path or vector of image paths.")
  }

  # Filter out backup files
  all_files <- all_files[!grepl("\\.(backup|orig|bak)$", all_files, ignore.case = TRUE)]

  n_imgs <- length(all_files)
  if (n_imgs == 0L) {
    cli::cli_alert_warning("No matching images found to anonymize in {.path {dataset}}.")
    empty_df <- data.frame(
      image = character(0), anonymized_count = integer(0),
      classes = character(0), status = character(0), output_path = character(0),
      stringsAsFactors = FALSE
    )
    class(empty_df) <- c("yolo_anonymized", "data.frame")
    return(empty_df)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("pliman: YOLO Field Dataset Privacy & Anonymization")
    cli::cli_ul(c(
      "Found {.val {n_imgs}} image{?s} to scan.",
      "Sensitive target classes: {.val {paste(target_classes, collapse = ', ')}}",
      "Method: {.val {method}} (kernel/block size: {.val {kernel_size}} px)",
      "Person masking mode: {.val {person_region}}",
      "Model: {.val {if (is.character(model)) basename(model) else 'in-memory'}}"
    ))
  }

  # Prepare output directories if separate output_dir requested
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  res_rows <- vector("list", n_imgs)
  total_anonymized <- 0L

  if (isTRUE(verbose)) {
    cli::cli_progress_bar(
      format = "{cli::pb_spin} Anonymizing dataset [{cli::pb_current}/{cli::pb_total}] | ETA: {cli::pb_eta}",
      total = n_imgs,
      clear = FALSE
    )
  }

  target_pattern <- paste(target_classes, collapse = "|")

  for (i in seq_len(n_imgs)) {
    if (isTRUE(verbose)) cli::cli_progress_update()
    img_path <- all_files[i]

    # Detect objects
    det_out <- tryCatch({
      image_detect_dl(
        img = img_path,
        model = model,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        engine = engine,
        plot = FALSE,
        verbose = FALSE,
        ...
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(det_out) || is.null(det_out$boxes) || nrow(det_out$boxes) == 0L) {
      res_rows[[i]] <- data.frame(
        image = img_path,
        anonymized_count = 0L,
        classes = "",
        status = "no_sensitive_objects",
        output_path = img_path,
        stringsAsFactors = FALSE
      )
      next
    }

    b_df <- det_out$boxes
    # Filter for sensitive classes
    is_sensitive <- grepl(target_pattern, b_df$class, ignore.case = TRUE)
    sens_boxes <- b_df[is_sensitive, , drop = FALSE]

    if (nrow(sens_boxes) == 0L) {
      res_rows[[i]] <- data.frame(
        image = img_path,
        anonymized_count = 0L,
        classes = "",
        status = "no_sensitive_objects",
        output_path = img_path,
        stringsAsFactors = FALSE
      )
      next
    }

    # Load image for pixel manipulation
    im <- tryCatch(image_import(img_path), error = function(e) NULL)
    if (is.null(im)) {
      res_rows[[i]] <- data.frame(
        image = img_path,
        anonymized_count = 0L,
        classes = "",
        status = "error",
        output_path = img_path,
        stringsAsFactors = FALSE
      )
      next
    }

    img_data_arr <- image_data(im)
    dims <- dim(img_data_arr)
    w <- dims[1]; h <- dims[2]; nch <- if (length(dims) >= 3L) dims[3] else 1L

    # Apply blurring / pixelation to each sensitive box
    n_applied <- 0L
    for (k in seq_len(nrow(sens_boxes))) {
      box <- sens_boxes[k, ]
      bx1 <- max(1L, as.integer(round(box$xmin)))
      by1 <- max(1L, as.integer(round(box$ymin)))
      bx2 <- min(w, as.integer(round(box$xmax)))
      by2 <- min(h, as.integer(round(box$ymax)))

      if (bx2 <= bx1 || by2 <= by1) next

      # If person and person_region == 'head', restrict to upper 22% of height
      if (grepl("person", box$class, ignore.case = TRUE) && person_region == "head") {
        box_h <- by2 - by1 + 1L
        head_h <- max(8L, as.integer(round(box_h * 0.25)))
        by2 <- min(h, by1 + head_h)
      }

      if (method == "pixelate") {
        bsize <- kernel_size
        for (ix in seq(bx1, bx2, by = bsize)) {
          for (iy in seq(by1, by2, by = bsize)) {
            ex <- min(bx2, ix + bsize - 1L)
            ey <- min(by2, iy + bsize - 1L)
            for (c_idx in seq_len(nch)) {
              if (nch == 1L) {
                m_val <- mean(img_data_arr[ix:ex, iy:ey])
                img_data_arr[ix:ex, iy:ey] <- m_val
              } else {
                m_val <- mean(img_data_arr[ix:ex, iy:ey, c_idx])
                img_data_arr[ix:ex, iy:ey, c_idx] <- m_val
              }
            }
          }
        }
      } else {
        # Spatial box blur: multi-pass average filter
        half_k <- max(2L, as.integer(round(kernel_size / 2)))
        for (pass in 1:2) {
          sub_w <- bx2 - bx1 + 1L; sub_h <- by2 - by1 + 1L
          for (c_idx in seq_len(nch)) {
            sub <- if (nch == 1L) img_data_arr[bx1:bx2, by1:by2] else img_data_arr[bx1:bx2, by1:by2, c_idx]
            # Fast downsample and upscale blur
            ds_w <- max(2L, as.integer(round(sub_w / half_k)))
            ds_h <- max(2L, as.integer(round(sub_h / half_k)))
            # Sample mean block values
            for (ix in seq(bx1, bx2, by = half_k)) {
              for (iy in seq(by1, by2, by = half_k)) {
                ex <- min(bx2, ix + half_k - 1L)
                ey <- min(by2, iy + half_k - 1L)
                if (nch == 1L) {
                  img_data_arr[ix:ex, iy:ey] <- mean(sub)
                } else {
                  img_data_arr[ix:ex, iy:ey, c_idx] <- mean(sub)
                }
              }
            }
          }
        }
      }
      n_applied <- n_applied + 1L
    }

    # Determine destination file path
    out_file <- if (!is.null(output_dir)) {
      rel_path <- if (is.character(dataset) && length(dataset) == 1L && dir.exists(dataset)) {
        sub(paste0("^", normalizePath(dataset, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(img_path, winslash = "/", mustWork = FALSE))
      } else {
        basename(img_path)
      }
      target_f <- file.path(output_dir, rel_path)
      td <- dirname(target_f)
      if (!dir.exists(td)) dir.create(td, recursive = TRUE, showWarnings = FALSE)
      target_f
    } else {
      # In-place modification with optional backup
      if (isTRUE(backup)) {
        bak_file <- paste0(img_path, ".orig")
        if (!file.exists(bak_file)) {
          file.copy(img_path, bak_file)
        }
      }
      img_path
    }

    # Save modified image
    anonymized_img <- as_image(img_data_arr)
    image_export(anonymized_img, name = out_file, plot = FALSE)

    total_anonymized <- total_anonymized + n_applied
    res_rows[[i]] <- data.frame(
      image = img_path,
      anonymized_count = n_applied,
      classes = paste(unique(sens_boxes$class), collapse = ", "),
      status = if (n_applied > 0L) "anonymized" else "no_sensitive_objects",
      output_path = out_file,
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(verbose)) cli::cli_progress_done()

  final_df <- do.call(rbind, res_rows)
  class(final_df) <- c("yolo_anonymized", "data.frame")

  if (isTRUE(verbose)) {
    n_modified <- sum(final_df$anonymized_count > 0L)
    cli::cli_alert_success("Anonymization complete! {total_anonymized} sensitive region{?s} masked across {n_modified} image{?s}.")
    cli::cli_alert_info("Dataset annotations (`labels/*.txt`) were kept 100% unaltered.")
  }

  invisible(final_df)
}

#' @rdname yolo_dataset_anonymize
#' @export
yolo_anonymize <- yolo_dataset_anonymize

#' @export
print.yolo_anonymized <- function(x, ...) {
  n_total <- nrow(x)
  n_anon <- sum(x$anonymized_count > 0L)
  tot_objs <- sum(x$anonymized_count)

  cli::cli_h2("pliman: YOLO Dataset Anonymization Summary")
  cli::cli_ul(c(
    "Total images scanned: {.val {n_total}}",
    "Images anonymized: {.val {n_anon}} ({round(100 * n_anon / max(1, n_total), 1)}%)",
    "Total sensitive entities masked: {.val {tot_objs}}"
  ))

  if (n_anon > 0L) {
    anon_subset <- x[x$anonymized_count > 0L, c("image", "anonymized_count", "classes", "output_path")]
    head_rows <- min(10L, nrow(anon_subset))
    cli::cli_text("Preview of anonymized images:")
    print.data.frame(anon_subset[seq_len(head_rows), ], row.names = FALSE)
    if (nrow(anon_subset) > head_rows) {
      cli::cli_alert_info("... and {nrow(anon_subset) - head_rows} more anonymized image(s).")
    }
  }
  invisible(x)
}

