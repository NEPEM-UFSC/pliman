# utils_yolo_train.R — Automated YOLO Dataset Export, Preview, and Training in pliman
#
# ==============================================================================
# YOLO DATASET & AUTO-ANNOTATION PIPELINE
# ==============================================================================

#' @title Export Dataset in YOLO Format for Detection or Segmentation
#' @name yolo_dataset_export
#' @description
#' Automatically creates a fully structured YOLO dataset (`images/`, `labels/`,
#' and `data.yaml`) for object detection or instance segmentation. Annotations can
#' be extracted automatically from [analyze_objects()], deep learning foundation models
#' (such as Grounded-SAM or PerSAM), or custom object contours.
#'
#' @details
#' The exported directory has the standard structure expected by YOLOv8, YOLO11, and YOLO26:
#' ```
#' dataset_dir/
#' ├── data.yaml
#' ├── images/
#' │   ├── train/
#' │   └── val/
#' └── labels/
#'     ├── train/
#'     └── val/
#' ```
#'
#' For segmentation (`task = "segment"`), each line in the `.txt` label file contains:
#' \code{class_id x1 y1 x2 y2 x3 y3 ... xn yn}
#' with polygon coordinates normalized between 0 and 1.
#'
#' For detection (`task = "detect"`), each line contains:
#' \code{class_id x_center y_center width height}
#'
#' @param img An `image` object, a list of `image` objects, or a directory path
#'   containing images.
#' @param annotations Annotations source. Can be:
#'   * An object returned by [analyze_objects()] (for a single image).
#'   * A list of objects returned by [analyze_objects()].
#'   * `"auto"` (default): automatically calls [analyze_objects()] on each image
#'     using the parameters supplied in `...`.
#'   * A list of polygon matrices or bounding boxes.
#' @param dir Output directory path for the YOLO dataset. Defaults to `"yolo_dataset"`.
#' @param task Task type: `"segment"` (default) for polygon instance segmentation
#'   or `"detect"` for bounding-box detection.
#' @param class_names Character vector of class names (e.g. `c("grain")` or `c("seed", "pod")`).
#'   Defaults to `c("object")`.
#' @param class_id Numeric class ID assigned to the annotations. Defaults to `0L`.
#' @param train_prop Proportion of images assigned to the training split (default `0.8`).
#'   The remaining images are assigned to the validation split (`val`).
#' @param n_points Target number of polygon vertices per object for segmentation
#'   (default `40`). Points are sampled evenly along the contour to keep label files
#'   compact and clean.
#' @param min_area Minimum object area (in pixels) to include in the dataset. Defaults to `15`.
#' @param overwrite Logical. If `TRUE` (default), overwrites existing dataset files in `dir`.
#' @param ... Additional arguments passed directly to [analyze_objects()] when
#'   `annotations = "auto"` (e.g. `index = "B"`, `watershed = TRUE`, `tolerance = 1`).
#'
#' @return An invisible list containing:
#'   * `dir`: Root directory of the dataset.
#'   * `yaml_file`: Full path to the generated `data.yaml`.
#'   * `train_images`: Vector of paths to training images.
#'   * `val_images`: Vector of paths to validation images.
#'   * `n_objects`: Total number of annotated objects across the dataset.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#' img <- image_import("cafe.jpeg")
#'
#' # Option A: Auto-annotate directly using analyze_objects parameters
#' yolo_dataset_export(
#'   img,
#'   dir = "dataset_cafe",
#'   task = "segment",
#'   class_names = "grao_cafe",
#'   index = "B",
#'   watershed = TRUE,
#'   tolerance = 1
#' )
#'
#' # Option B: Pass pre-computed analyze_objects result
#' res <- analyze_objects(img, index = "B", watershed = TRUE, tolerance = 1)
#' yolo_dataset_export(
#'   img,
#'   annotations = res,
#'   dir = "dataset_cafe",
#'   class_names = "grao_cafe"
#' )
#' }
yolo_dataset_export <- function(img,
                                annotations = "auto",
                                dir = "yolo_dataset",
                                task = c("segment", "detect"),
                                class_names = c("object"),
                                class_id = 0L,
                                train_prop = 0.8,
                                n_points = 40L,
                                min_area = 15,
                                overwrite = TRUE,
                                ...) {
  task <- match.arg(task)
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)

  # Prepare directory structure
  train_img_dir <- file.path(dir, "images", "train")
  val_img_dir   <- file.path(dir, "images", "val")
  train_lbl_dir <- file.path(dir, "labels", "train")
  val_lbl_dir   <- file.path(dir, "labels", "val")

  for (d in c(train_img_dir, val_img_dir, train_lbl_dir, val_lbl_dir)) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Resolve images input
  img_list <- list()
  img_names <- character()

  if (is.character(img) && length(img) == 1L && dir.exists(img)) {
    # Directory path
    files <- list.files(img, pattern = "\\.(jpg|jpeg|png|tif|tiff|bmp)$", ignore.case = TRUE, full.names = TRUE)
    if (length(files) == 0L) {
      cli::cli_abort("No valid images found in directory {.path {img}}.")
    }
    img_names <- tools::file_path_sans_ext(basename(files))
    img_list <- lapply(files, image_import)
  } else if (is.character(img)) {
    # Vector of file paths
    img_names <- tools::file_path_sans_ext(basename(img))
    img_list <- lapply(img, image_import)
  } else if (inherits(img, "image") || inherits(img, "Image")) {
    # Single image
    img_list <- list(img)
    img_names <- if (!is.null(names(img)) && nzchar(names(img)[1])) names(img)[1] else "img_1"
  } else if (is.list(img)) {
    # List of images
    img_list <- img
    img_names <- names(img)
    if (is.null(img_names) || any(!nzchar(img_names))) {
      img_names <- paste0("img_", seq_along(img))
    }
  } else {
    cli::cli_abort("Invalid {.arg img} argument. Must be an image, list of images, or directory path.")
  }

  num_imgs <- length(img_list)
  if (num_imgs == 0L) {
    cli::cli_abort("No images to export.")
  }

  # Train / Val Split
  if (num_imgs == 1L) {
    is_train <- c(TRUE)
  } else {
    n_train <- max(1L, round(num_imgs * train_prop))
    if (n_train >= num_imgs && num_imgs > 1L) n_train <- num_imgs - 1L
    set.seed(42)
    train_idx <- sample(seq_len(num_imgs), size = n_train)
    is_train <- seq_len(num_imgs) %in% train_idx
  }

  total_objects_annotated <- 0L
  train_images_saved <- character()
  val_images_saved   <- character()

  cli::cli_progress_bar("Exporting YOLO dataset", total = num_imgs)

  for (i in seq_len(num_imgs)) {
    curr_img <- img_list[[i]]
    base_name <- img_names[i]
    split_name <- if (is_train[i]) "train" else "val"

    dest_img_dir <- if (is_train[i]) train_img_dir else val_img_dir
    dest_lbl_dir <- if (is_train[i]) train_lbl_dir else val_lbl_dir

    dest_img_path <- file.path(dest_img_dir, paste0(base_name, ".jpg"))
    dest_lbl_path <- file.path(dest_lbl_dir, paste0(base_name, ".txt"))

    # Save Image
    image_export(curr_img, name = paste0(base_name, ".jpg"), subfolder = dest_img_dir)

    if (is_train[i]) {
      train_images_saved <- c(train_images_saved, dest_img_path)
    } else {
      val_images_saved <- c(val_images_saved, dest_img_path)
    }

    # Extract contours
    conts <- list()

    if (identical(annotations, "auto")) {
      # Automatically run analyze_objects on the image
      res_obj <- analyze_objects(curr_img, plot = FALSE, ...)
      if (!is.null(res_obj[["contours"]])) {
        conts <- res_obj[["contours"]]
      }
    } else if (inherits(annotations, "anal_obj")) {
      # Single pre-computed analyze_objects result
      if (!is.null(annotations[["contours"]])) {
        conts <- annotations[["contours"]]
      }
    } else if (is.list(annotations) && length(annotations) == num_imgs && inherits(annotations[[i]], "anal_obj")) {
      # List of analyze_objects results
      if (!is.null(annotations[[i]][["contours"]])) {
        conts <- annotations[[i]][["contours"]]
      }
    } else if (is.list(annotations) && all(sapply(annotations, is.matrix))) {
      # Custom list of contour matrices
      conts <- annotations
    } else if (is.function(annotations)) {
      # Custom extractor function
      res_custom <- annotations(curr_img)
      if (inherits(res_custom, "anal_obj") && !is.null(res_custom[["contours"]])) {
        conts <- res_custom[["contours"]]
      } else if (is.list(res_custom)) {
        conts <- res_custom
      }
    }

    # Image dimensions in pliman: [width, height, channels]
    dims <- dim(curr_img)
    img_w <- dims[1]
    img_h <- dims[2]

    txt_lines <- character()

    if (length(conts) > 0L) {
      for (k in seq_along(conts)) {
        c_mat <- conts[[k]]
        if (is.null(c_mat) || !is.matrix(c_mat) && !is.data.frame(c_mat)) next
        c_mat <- as.matrix(c_mat)
        if (nrow(c_mat) < 4L) next

        # Filter by bounding box area approximation
        x_vals <- c_mat[, 1]
        y_vals <- c_mat[, 2]
        w_obj <- max(x_vals) - min(x_vals)
        h_obj <- max(y_vals) - min(y_vals)
        if ((w_obj * h_obj) < min_area) next

        # Resample / Decimate contour vertices for clean lightweight annotation
        n_cur <- nrow(c_mat)
        target_n <- min(as.integer(n_points), n_cur)
        if (n_cur > target_n) {
          sampled_idx <- round(seq(1, n_cur, length.out = target_n + 1L)[-(target_n + 1L)])
          c_mat <- c_mat[sampled_idx, , drop = FALSE]
        }

        # Normalize coordinates between 0 and 1
        x_norm <- pmin(pmax(c_mat[, 1] / img_w, 0.0), 1.0)
        y_norm <- pmin(pmax(c_mat[, 2] / img_h, 0.0), 1.0)

        cid <- if (length(class_id) == length(conts)) class_id[k] else class_id[1]

        if (task == "segment") {
          # YOLO-seg format: class_id x1 y1 x2 y2 ... xn yn
          interleaved <- as.vector(rbind(x_norm, y_norm))
          line_entry <- paste(c(cid, round(interleaved, 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        } else {
          # YOLO-det format: class_id cx cy w h
          xmin <- min(x_norm); xmax <- max(x_norm)
          ymin <- min(y_norm); ymax <- max(y_norm)
          cx <- (xmin + xmax) / 2.0
          cy <- (ymin + ymax) / 2.0
          bw <- xmax - xmin
          bh <- ymax - ymin
          line_entry <- paste(c(cid, round(c(cx, cy, bw, bh), 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      }
    }

    # Write text labels
    writeLines(txt_lines, dest_lbl_path)
    total_objects_annotated <- total_objects_annotated + length(txt_lines)

    cli::cli_progress_update()
  }

  cli::cli_progress_done()

  # Generate data.yaml
  yaml_path <- file.path(dir, "data.yaml")
  yaml_content <- c(
    paste0("path: ", normalizePath(dir, winslash = "/")),
    "train: images/train",
    "val: images/val",
    "",
    "names:"
  )
  for (idx in seq_along(class_names)) {
    yaml_content <- c(yaml_content, paste0("  ", idx - 1L, ": ", class_names[idx]))
  }

  writeLines(yaml_content, yaml_path)

  cli::cli_alert_success(
    "YOLO dataset successfully generated in {.path {dir}} ({total_objects_annotated} objects across {num_imgs} images)."
  )

  invisible(list(
    dir = dir,
    yaml_file = yaml_path,
    train_images = train_images_saved,
    val_images = val_images_saved,
    n_objects = total_objects_annotated
  ))
}


#' @title Preview Annotations in a YOLO Dataset
#' @name yolo_dataset_preview
#' @description
#' Loads sample images and their matching `.txt` annotation labels from a YOLO
#' dataset directory and displays them with overlaid segmentation polygons or bounding boxes.
#'
#' @param dir Path to the YOLO dataset directory (containing `images/` and `labels/`).
#' @param split Split to inspect: `"train"` (default) or `"val"`.
#' @param n Number of sample images to display (default `4`).
#' @param col Border color for contours or bounding boxes. Defaults to `"red"`.
#' @param lwd Line width for plotting. Defaults to `2`.
#' @param fill Optional fill color for polygons (e.g. `"#FF000033"` for translucent red).
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#' yolo_dataset_preview("dataset_cafe", split = "train", n = 1)
#' }
yolo_dataset_preview <- function(dir = "yolo_dataset",
                                 split = c("train", "val"),
                                 n = 4L,
                                 col = "red",
                                 lwd = 2,
                                 fill = NULL) {
  split <- match.arg(split)
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)

  img_dir <- file.path(dir, "images", split)
  lbl_dir <- file.path(dir, "labels", split)

  if (!dir.exists(img_dir) || !dir.exists(lbl_dir)) {
    cli::cli_abort("Could not find dataset directories for split {.val {split}} in {.path {dir}}.")
  }

  img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png)$", ignore.case = TRUE, full.names = TRUE)
  if (length(img_files) == 0L) {
    cli::cli_abort("No images found in {.path {img_dir}}.")
  }

  n_samples <- min(as.integer(n), length(img_files))
  sample_files <- img_files[seq_len(n_samples)]

  for (img_path in sample_files) {
    base_name <- tools::file_path_sans_ext(basename(img_path))
    lbl_path  <- file.path(lbl_dir, paste0(base_name, ".txt"))

    img <- image_import(img_path)
    plot(img)

    dims <- dim(img)
    img_w <- dims[1]
    img_h <- dims[2]

    if (file.exists(lbl_path)) {
      lines <- readLines(lbl_path, warn = FALSE)
      lines <- lines[nzchar(trimws(lines))]

      for (ln in lines) {
        vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
        if (length(vals) < 5L) next

        cls_id <- vals[1]
        coords <- vals[-1]

        if (length(coords) == 4L) {
          # Detection format: cx cy w h
          cx <- coords[1] * img_w
          cy <- coords[2] * img_h
          bw <- coords[3] * img_w
          bh <- coords[4] * img_h
          xmin <- cx - bw / 2
          xmax <- cx + bw / 2
          ymin <- cy - bh / 2
          ymax <- cy + bh / 2
          graphics::rect(xmin, ymin, xmax, ymax, border = col, lwd = lwd)
        } else if (length(coords) >= 6L && length(coords) %% 2 == 0) {
          # Segmentation polygon format: x1 y1 x2 y2 ...
          xs <- coords[seq(1, length(coords), by = 2L)] * img_w
          ys <- coords[seq(2, length(coords), by = 2L)] * img_h
          graphics::polygon(xs, ys, border = col, col = fill, lwd = lwd)
        }
      }
    }
  }

  invisible(sample_files)
}


#' @title Train or Fine-Tune a YOLO Model from R
#' @name yolo_train
#' @description
#' Automatically trains or fine-tunes a YOLO model (YOLOv8, YOLO11, or YOLO26)
#' directly from R on an exported dataset, and converts the resulting best checkpoint
#' into a standalone ONNX model ready for instant C++ execution in [image_segment_dl()]
#' or [image_detect_dl()].
#'
#' @details
#' This function provides a completely turnkey experience:
#' 1. Validates the `data.yaml` dataset structure.
#' 2. Verifies the Python / Ultralytics runtime. If missing, it can automatically install
#'    the lightweight dependencies (`pip install ultralytics onnx`).
#' 3. Runs the training job with hardware acceleration (GPU / DirectML / CPU).
#' 4. Automatically exports the resulting weights (`best.pt`) directly to `.onnx`.
#' 5. Copies the final `.onnx` model into `output_dir` (default: `pliman_model_dir()`).
#'
#' @param data Path to the `data.yaml` file (e.g. `"dataset_cafe/data.yaml"`).
#' @param model Base pre-trained model architecture to fine-tune. Defaults to `"yolo11n-seg"`.
#'   Can be `"yolo11n-seg"`, `"yolo11s-seg"`, `"yolo26n-seg"`, `"yolo11n"`, or path to a `.pt` file.
#' @param epochs Number of training epochs (default `50`).
#' @param imgsz Input image resolution for training (default `640`).
#' @param batch Batch size (default `16`).
#' @param device Device to train on: `"auto"` (default), `"cpu"`, `"0"` (GPU 0), etc.
#' @param output_dir Target directory to save the final exported `.onnx` model.
#'   Defaults to [pliman_model_dir()].
#' @param output_name Filename of the exported model. If `NULL` (default),
#'   an informative name is automatically created (e.g. `"yolo_custom_seg.onnx"`).
#' @param patience Early stopping patience in epochs. Defaults to `15`.
#' @param workers Number of worker threads for data loading. Defaults to `4`.
#' @param auto_install Logical. If `TRUE`, automatically installs `ultralytics` and `onnx`
#'   via `pip` if they are not yet available on the system. Defaults to `TRUE`.
#' @param ... Additional arguments passed to the training command (e.g. `lr0 = 0.01`).
#'
#' @return A list with:
#'   * `model_file`: Path to the exported `.onnx` model ready for `pliman`.
#'   * `results_dir`: Directory containing training plots, metrics, and logs.
#'   * `epochs`: Number of epochs trained.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Train a custom coffee grain segmentation model
#' model_res <- yolo_train(
#'   data = "dataset_cafe/data.yaml",
#'   model = "yolo11n-seg",
#'   epochs = 30,
#'   output_name = "cafe_seg.onnx"
#' )
#'
#' # Run inference with your new model directly in pliman!
#' img <- image_import("cafe.jpeg")
#' seg <- image_segment_dl(img, model = model_res$model_file, type = "highlight")
#' }
yolo_train <- function(data = "yolo_dataset/data.yaml",
                       model = "yolo11n-seg",
                       epochs = 50,
                       imgsz = 640,
                       batch = 16,
                       device = "auto",
                       output_dir = pliman_model_dir(),
                       output_name = NULL,
                       patience = 15,
                       workers = 4,
                       auto_install = TRUE,
                       ...) {
  data <- normalizePath(data, winslash = "/", mustWork = FALSE)
  if (!file.exists(data)) {
    cli::cli_abort("Could not find dataset configuration file at {.path {data}}.")
  }

  output_dir <- pliman_model_dir(output_dir)

  # Check python availability
  py_exec <- Sys.which("python")
  if (!nzchar(py_exec)) {
    py_exec <- Sys.which("python3")
  }
  if (!nzchar(py_exec)) {
    cli::cli_abort(
      c(
        "x" = "Python is required for model weight optimization (backpropagation).",
        "i" = "Please install Python 3.9+ from https://www.python.org or your package manager."
      )
    )
  }

  # Check if ultralytics and onnx are installed
  check_cmd <- paste0('"', py_exec, '" -c "import ultralytics, onnx; print(\'OK\')"')
  check_res <- suppressWarnings(system(check_cmd, intern = TRUE, ignore.stderr = TRUE))

  if (!identical(check_res, "OK")) {
    if (isTRUE(auto_install)) {
      cli::cli_alert_info("Installing {.val ultralytics} and {.val onnx} via pip...")
      install_status <- system(paste0('"', py_exec, '" -m pip install ultralytics onnx'))
      if (install_status != 0L) {
        cli::cli_abort("Failed to automatically install ultralytics. Run {.code pip install ultralytics onnx} in your terminal.")
      }
      cli::cli_alert_success("Ultralytics and ONNX successfully installed!")
    } else {
      cli::cli_abort("Ultralytics is not installed. Run {.code pip install ultralytics onnx} or set {.code auto_install = TRUE}.")
    }
  }

  # Normalize model name / path
  base_model <- model
  if (!grepl("\\.(pt|onnx|yaml)$", base_model)) {
    base_model <- paste0(base_model, ".pt")
  }

  # Build Python script to run training and immediate ONNX export
  run_dir <- file.path(tempdir(), "yolo_train_run")
  if (!dir.exists(run_dir)) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  }

  py_script_path <- file.path(run_dir, "train_and_export.py")

  # Detect task (segment or detect) from model name
  is_seg <- grepl("seg", base_model, ignore.case = TRUE)
  task_mode <- if (is_seg) "segment" else "detect"

  default_output_name <- if (is_seg) "custom_yolo_seg.onnx" else "custom_yolo_det.onnx"
  if (is.null(output_name) || !nzchar(output_name)) {
    output_name <- default_output_name
  }
  if (!grepl("\\.onnx$", output_name, ignore.case = TRUE)) {
    output_name <- paste0(output_name, ".onnx")
  }

  dest_onnx <- file.path(output_dir, output_name)

  py_code <- paste0(
    "import sys\n",
    "from ultralytics import YOLO\n\n",
    "print('--> Starting YOLO training [task=", task_mode, "]...')\n",
    "model = YOLO('", base_model, "')\n",
    "results = model.train(\n",
    "    data='", data, "',\n",
    "    epochs=", as.integer(epochs), ",\n",
    "    imgsz=", as.integer(imgsz), ",\n",
    "    batch=", as.integer(batch), ",\n",
    "    device='", device, "',\n",
    "    patience=", as.integer(patience), ",\n",
    "    workers=", as.integer(workers), ",\n",
    "    project='", normalizePath(run_dir, winslash = "/"), "',\n",
    "    name='train_run',\n",
    "    exist_ok=True\n",
    ")\n\n",
    "print('--> Exporting best model to ONNX format...')\n",
    "best_pt = '", normalizePath(run_dir, winslash = "/"), "/train_run/weights/best.pt'\n",
    "best_model = YOLO(best_pt)\n",
    "exported_path = best_model.export(format='onnx', imgsz=", as.integer(imgsz), ")\n",
    "print('--> ONNX exported to:', exported_path)\n"
  )

  writeLines(py_code, py_script_path)

  cli::cli_alert_info("Launching YOLO training ({epochs} epochs, imgsz={imgsz})...")

  start_time <- Sys.time()
  train_cmd <- paste0('"', py_exec, '" "', py_script_path, '"')
  train_status <- system(train_cmd)
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)

  if (train_status != 0L) {
    cli::cli_abort("Training process failed with exit code {train_status}.")
  }

  # Find exported ONNX file
  cand_onnx <- file.path(run_dir, "train_run", "weights", "best.onnx")
  if (!file.exists(cand_onnx)) {
    # Search in run_dir
    found <- list.files(run_dir, pattern = "\\.onnx$", recursive = TRUE, full.names = TRUE)
    if (length(found) > 0L) {
      cand_onnx <- found[1]
    }
  }

  if (!file.exists(cand_onnx)) {
    cli::cli_abort("Could not find the exported ONNX model in {.path {run_dir}}.")
  }

  # Copy to destination in models directory
  file.copy(cand_onnx, dest_onnx, overwrite = TRUE)

  cli::cli_alert_success("Model trained ({elapsed}s) and saved to {.path {dest_onnx}}!")

  invisible(list(
    model_file = dest_onnx,
    results_dir = file.path(run_dir, "train_run"),
    epochs = epochs,
    elapsed_sec = elapsed
  ))
}
