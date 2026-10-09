# utils_yolo_train.R - Automated YOLO Dataset Export, Preview, and Training in pliman
#
# ==============================================================================
# YOLO DATASET & AUTO-ANNOTATION PIPELINE
# ==============================================================================

# Helper: Draw real-time visual feedback for interactive PerSAM labeling
.draw_persam_feedback <- function(img,
                                  boxes = NULL,
                                  contours = NULL,
                                  similarity_map = NULL,
                                  task = "detect",
                                  class_name = "object",
                                  show_similarity = FALSE,
                                  sim_threshold = 0.5,
                                  img_idx = 1L,
                                  total_imgs = 1L) {
  n_boxes <- if (!is.null(boxes) && is.data.frame(boxes)) nrow(boxes) else 0L
  n_conts <- if (!is.null(contours) && is.list(contours)) length(contours) else 0L
  n_objs  <- max(n_boxes, n_conts)

  # Rainbow palette for distinguishing each object individually with matching mask + bbox
  pal <- if (n_objs > 0L) grDevices::rainbow(n_objs, s = 0.85, v = 0.95) else character()

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (isTRUE(show_similarity) && !is.null(similarity_map)) {
    graphics::par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
  } else {
    graphics::par(mar = c(2, 2, 3, 1))
  }

  # Left: Image with overlaid detections
  plot(img)

  # 1. Draw segmented masks (polygons) with matching rainbow color
  if (identical(task, "segment") && n_conts > 0L) {
    for (k in seq_len(n_conts)) {
      cm <- contours[[k]]
      if (is.null(cm) || nrow(cm) < 3) next
      col_k <- pal[((k - 1L) %% length(pal)) + 1L]
      col_fill <- grDevices::adjustcolor(col_k, alpha.f = 0.35)
      graphics::polygon(cm[, 1], cm[, 2], col = col_fill, border = col_k, lwd = 2)
    }
  }

  # 2. Draw bounding boxes and tags with the EXACT same rainbow color per object
  if (n_boxes > 0L) {
    for (b in seq_len(n_boxes)) {
      if (identical(task, "segment") && !is.null(contours)) {
        if (b > length(contours) || is.null(contours[[b]]) || (is.matrix(contours[[b]]) && nrow(contours[[b]]) < 3)) {
          next
        }
      }
      bx <- boxes[b, ]
      col_b <- pal[((b - 1L) %% length(pal)) + 1L]
      graphics::rect(bx$xmin, bx$ymin, bx$xmax, bx$ymax, border = col_b, lwd = 2)

      score_txt <- if (!is.null(bx$score) && !is.na(bx$score)) sprintf("%.2f", bx$score) else ""
      lbl_txt <- if (!is.null(bx$label) && bx$label != "exemplar") bx$label else class_name
      tag <- if (nzchar(score_txt)) paste0(lbl_txt, ": ", score_txt) else lbl_txt

      tw <- graphics::strwidth(tag, cex = 0.7)
      th <- graphics::strheight(tag, cex = 0.7)
      y_top <- max(th + 6, bx$ymin)

      rgb_b <- grDevices::col2rgb(col_b)
      lum <- (0.299 * rgb_b[1] + 0.587 * rgb_b[2] + 0.114 * rgb_b[3]) / 255.0
      txt_col <- if (lum > 0.65) "#000000" else "#FFFFFF"
      bg_tag <- grDevices::adjustcolor(col_b, alpha.f = 0.90)

      graphics::rect(bx$xmin, y_top - th - 6, bx$xmin + tw + 6, y_top,
                     col = bg_tag, border = NA)
      graphics::text(bx$xmin + 3, y_top - 3, labels = tag,
                     col = txt_col, cex = 0.7, font = 2, adj = c(0, 1))
    }
  }

  title_str <- sprintf("Image [%d/%d] - %d %s detected (sim >= %.2f)",
                       img_idx, total_imgs, n_boxes, class_name, sim_threshold)
  graphics::title(main = title_str, font.main = 2, col.main = "#006622", cex.main = 1.05)

  # Right: Feature Similarity Map ("Matriz de Coincidencia")
  if (isTRUE(show_similarity) && !is.null(similarity_map)) {
    plot(similarity_map)
    graphics::title(main = "PerSAM Cosine Similarity Map (Matriz de Coincid\u00eancia)",
                    font.main = 2, col.main = "#003366", cex.main = 1.05)
  }
}

#' @title Export Dataset in YOLO Format for Detection or Segmentation
#' @name yolo_dataset_export
#' @description
#' Automatically generates a structured YOLO dataset (`images/`, `labels/`,
#' and `data.yaml`) for object detection (`task = "detect"`) or instance
#' segmentation (`task = "segment"`). Annotations can be produced automatically
#' via deep learning foundation models:
#' \itemize{
#'   \item \strong{One-shot visual exemplar auto-labeling} with \code{model = "persam"} (Roboflow-style 1-click annotation).
#'   \item \strong{Zero-shot text prompt detection} with \code{model = "grounded-sam"} (e.g. \code{prompt = "people"}).
#'   \item \strong{Assisted labeling / pseudo-labeling} with pre-trained or fine-tuned YOLO models (\code{model = "custom.onnx"}).
#'   \item \strong{Classical image analysis} with [analyze_objects()] (\code{model = NULL}).
#' }
#'
#' @details
#' \subsection{Automated and Semi-Automated Annotation Modes}{
#' \itemize{
#'   \item \strong{PerSAM One-Shot Exemplar Auto-Labeling (\code{model = "persam"}):}
#'     The user clicks on 1 or more exemplar objects in the plot window. PerSAM extracts
#'     visual embeddings from SAM 2.1, calculates the spatial cosine feature similarity map
#'     ("matriz de coincidencia"), and identifies all matching instances across the image.
#'     When \code{interactive = TRUE}, an interactive confirmation menu allows accepting detections (<Enter>),
#'     retrying clicks (\code{'r'}), or fine-tuning the similarity threshold (\code{'t 0.55'}) in real-time.
#'   \item \strong{Grounded-SAM Text Prompting (\code{model = "grounded-sam"}):}
#'     Accepts natural language text descriptions (e.g. \code{prompt = "people"}, \code{prompt = "coffee grain"},
#'     or \code{class_names = c("person", "bicycle")}) to locate and segment all matching objects.
#'   \item \strong{YOLO Pre-Labeling / Active Learning (\code{model = "<model>.onnx"}):}
#'     Applies an existing trained YOLO model to unannotated images to produce draft labels
#'     for rapid dataset expansion.
#'   \item \strong{Classical Image Analysis (\code{model = NULL}):}
#'     Extracts contours via [analyze_objects()] using thresholding, color indexes, and watershed algorithms.
#' }
#' }
#'
#' The exported directory adheres to the standard structure expected by YOLOv8, YOLO11, and YOLO26:
#' ```
#' dataset_dir/
#' |-- data.yaml
#' |-- images/
#' |   |-- train/
#' |   \-- val/
#' \-- labels/
#'     |-- train/
#'     \-- val/
#' ```
#'
#' @param img An `image` object, a list of `image` objects, a file path, or a directory path
#'   containing images to annotate.
#' @param model Deep learning model to use for automated or semi-automated labeling:
#'   * `"persam"`: One-shot visual exemplar segmentation and detection (Roboflow-style 1-click auto-labeling).
#'   * `"grounded-sam"`: Zero-shot text-prompted instance segmentation & detection.
#'   * Path to a trained YOLO `.onnx` model (e.g. `"capsulas_det.onnx"`): Assisted labeling / pseudo-labeling.
#'   * `"ben2"`, `"u2netp"`, `"birefnet-lite"`: Salient foreground object auto-segmentation.
#'   * `NULL` (default): Uses classical computer vision segmentation via [analyze_objects()].
#' @param annotations Annotations source when `model = NULL`. Can be:
#'   * `"auto"` (default): automatically calls [analyze_objects()] on each image.
#'   * An object or list of objects returned by [analyze_objects()].
#'   * A custom list of contour matrices or bounding boxes.
#' @param dir Output directory path for the YOLO dataset. Defaults to `"yolo_dataset"`.
#' @param task Task type: `"detect"` (default) for bounding-box detection or `"segment"`
#'   for polygon instance segmentation.
#' @param class_names Character vector of class names (e.g. `c("capsula")` or `c("grain", "defect")`).
#'   Defaults to `c("object")`.
#' @param class_id Numeric class ID assigned to the annotations. Defaults to `0L`.
#' @param prompt Prompt specification for deep learning models:
#'   * For `model = "persam"`: Exemplar coordinates (`c(x, y)` or matrix). If `NULL` and `interactive = TRUE`,
#'     prompts for clicks on the image.
#'   * For `model = "grounded-sam"`: Character text prompt (e.g. `"people"`, `"grain"`, `"leaf"`).
#'     Defaults to `class_names` if not supplied.
#' @param prompt_mode Exemplar prompting mode when labeling multiple images with PerSAM:
#'   * `"each"` (default for interactive): Prompts for exemplar clicks on each image.
#'   * `"first"`: Prompts for exemplar clicks only on the first image, applying the points across all images.
#' @param interactive Logical. If `TRUE`, activates the interactive visual labeling workflow in the plot window,
#'   allowing the user to click exemplars, review instant detections with confidence tags, and confirm or adjust.
#'   Defaults to `TRUE` if `model = "persam"`, `prompt = NULL`, and R is running interactively.
#' @param sim_threshold Numeric cosine similarity threshold in `[0, 1]` for `model = "persam"`. Defaults to `0.50`.
#' @param min_dist Minimum distance (in pixels) between detected objects/peaks for PerSAM exemplar segmentation. Defaults to `20`.
#' @param box_threshold Bounding box confidence threshold for Grounded-SAM or YOLO models. Defaults to `0.25`.
#' @param text_threshold Text confidence threshold for Grounded-SAM. Defaults to `0.25`.
#' @param iou_threshold Non-Maximum Suppression (NMS) IoU threshold. Defaults to `0.50`.
#' @param feat_res Feature map resolution for PerSAM (`64`, `256` (default), `"sr"`, or `1024`).
#' @param superres_map Logical. If `TRUE`, applies Real-ESRGAN 4x super-resolution on the PerSAM similarity map
#'   with bilateral semantic background gating. Default is `FALSE`.
#' @param max_objects Maximum number of objects detected per image. Defaults to `300`.
#' @param filter_edge Logical. If `TRUE` (default), filters out spurious border-cut artifacts touching image edges.
#' @param show_similarity Logical. If `TRUE`, displays a side-by-side 2-panel preview showing the detected
#'   objects on the left and the PerSAM Cosine Similarity Map ("Matriz de Coincidencia") on the right. Defaults to `FALSE`.
#' @param train_prop Proportion of images assigned to the training split (default `0.8`).
#'   The remaining images are assigned to the validation split (`val`).
#' @param n_points Target number of polygon vertices per object for segmentation
#'   (default `40`). Points are sampled evenly along the contour to keep label files
#'   compact and clean.
#' @param min_area Minimum object area (in pixels) to include in the dataset. Defaults to `15`.
#' @param engine Computing device engine: `"cpu"` (default) or `"gpu"`.
#' @param device_id GPU device ID (default `-1`).
#' @param overwrite Logical. If `TRUE` (default), overwrites existing dataset files in `dir`.
#' @param verbose Logical. If `TRUE` (default), displays live progress information.
#' @param ... Additional arguments passed directly to [analyze_objects()] when
#'   `model = NULL` and `annotations = "auto"` (e.g. `index = "B"`, `watershed = TRUE`, `tolerance = 1`).
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
#'
#' # 1. Semi-Automated 1-Click Labeling with PerSAM:
#' # Click 1 capsule on screen -> PerSAM finds all 40+ capsules and writes YOLO boxes!
#' yolo_dataset_export(
#'   "capsulas_raw/",
#'   model = "persam",
#'   task = "detect",
#'   class_names = "capsula",
#'   dir = "dataset_capsulas"
#' )
#'
#' # 2. Zero-Shot Text Prompt Labeling with Grounded-SAM:
#' # Annotate every person in photos automatically:
#' yolo_dataset_export(
#'   "crowd_photos/",
#'   model = "grounded-sam",
#'   prompt = "people",
#'   class_names = "person",
#'   dir = "dataset_people"
#' )
#'
#' # 3. Classical Auto-Annotation with analyze_objects:
#' yolo_dataset_export(
#'   "cafe_grains/",
#'   task = "segment",
#'   class_names = "grao_cafe",
#'   index = "B",
#'   watershed = TRUE
#' )
#' }
yolo_dataset_export <- function(img,
                                model = NULL,
                                annotations = "auto",
                                dir = "yolo_dataset",
                                task = c("detect", "segment"),
                                class_names = c("object"),
                                class_id = 0L,
                                prompt = NULL,
                                prompt_mode = c("each", "first"),
                                interactive = is.null(prompt) && !is.null(model) && identical(tolower(model), "persam") && interactive(),
                                sim_threshold = 0.50,
                                min_dist = 20,
                                box_threshold = 0.25,
                                text_threshold = 0.25,
                                iou_threshold = 0.50,
                                feat_res = 256,
                                max_objects = NULL,
                                train_prop = 0.8,
                                n_points = 40L,
                                min_area = 15,
                                filter_edge = TRUE,
                                show_similarity = FALSE,
                                superres_map = FALSE,
                                engine = c("gpu", "cpu"),
                                device_id = -1,
                                overwrite = TRUE,
                                verbose = TRUE,
                                ...) {
  task <- match.arg(task)
  engine <- match.arg(engine)
  prompt_mode <- match.arg(prompt_mode)
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)

  # Auto-detect Grounded-SAM when character text prompt is given
  if (is.null(model) && is.character(prompt) && length(prompt) >= 1L && !all(prompt %in% c("center", "box", "exemplar"))) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Text prompt detected ({.val {paste(prompt, collapse = ', ')}}). Automatically activating {.val grounded-sam} model.")
    }
    model <- "grounded-sam"
  }

  # Auto-derive class_names from text prompt if left at default "object"
  if (identical(class_names, "object") && is.character(prompt) && length(prompt) >= 1L) {
    raw_tokens <- unlist(strsplit(prompt, "[,;.]|\\band\\b", perl = TRUE))
    clean_tokens <- trimws(raw_tokens)
    clean_tokens <- clean_tokens[nzchar(clean_tokens)]
    if (length(clean_tokens) > 0L) {
      class_names <- unique(clean_tokens)
    }
  }

  m_lower <- if (!is.null(model) && is.character(model)) tolower(model[1]) else ""
  is_persam <- m_lower %in% c("persam", "per-sam", "sam-exemplar")
  is_grounded_sam <- m_lower %in% c("grounded-sam", "grounded_sam")
  is_yolo <- !is_persam && !is_grounded_sam && nzchar(m_lower) && (
    grepl("yolo", m_lower) || grepl("\\.(onnx|pt)$", model[1]) || file.exists(model[1])
  )
  is_salient <- m_lower %in% c("ben2", "u2netp", "u2net", "birefnet-lite", "rmbg-1.4", "rmbg-2.0", "withoutbg", "silueta")

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
    files <- list.files(img, pattern = "\\.(jpg|jpeg|png|tif|tiff|bmp)$", ignore.case = TRUE, full.names = TRUE)
    if (length(files) == 0L) {
      cli::cli_abort("No valid images found in directory {.path {img}}.")
    }
    img_names <- tools::file_path_sans_ext(basename(files))
    img_list <- lapply(files, image_import)
  } else if (is.character(img)) {
    img_names <- tools::file_path_sans_ext(basename(img))
    img_list <- lapply(img, image_import)
  } else if (inherits(img, "image") || inherits(img, "Image")) {
    img_list <- list(img)
    img_names <- if (!is.null(names(img)) && nzchar(names(img)[1])) names(img)[1] else "img_1"
  } else if (is.list(img)) {
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

  if (isTRUE(verbose)) {
    mode_str <- if (is_persam) {
      if (isTRUE(interactive)) "PerSAM (Interactive 1-click Exemplar)" else "PerSAM (One-shot Exemplar Batch)"
    } else if (is_grounded_sam) {
      paste0("Grounded-SAM (Text Prompt: '", if (is.character(prompt)) paste(prompt, collapse = ", ") else class_names[1], "')")
    } else if (is_yolo) {
      paste0("YOLO Pseudo-Labeling (", basename(model[1]), ")")
    } else if (is_salient) {
      paste0("Salient Object Segmentation (", model[1], ")")
    } else {
      "Classical Image Analysis (analyze_objects)"
    }

    cli::cli_h2("YOLO Dataset Generation & Auto-Annotation")
    cli::cli_ul(c(
      "Annotation engine: {.val {mode_str}}",
      "Task: {.val {toupper(task)}}",
      "Classes: {.val {class_names}}",
      "Images to process: {num_imgs} ({sum(is_train)} train / {sum(!is_train)} val)",
      "Target dataset directory: {.path {dir}}"
    ))
  }

  cached_pts <- NULL
  if (!is.null(prompt) && (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt))) {
    cached_pts <- prompt
  }

  stop_early <- FALSE

  for (i in seq_len(num_imgs)) {
    curr_img <- img_list[[i]]
    base_name <- img_names[i]
    split_name <- if (is_train[i]) "train" else "val"

    dest_img_dir <- if (is_train[i]) train_img_dir else val_img_dir
    dest_lbl_dir <- if (is_train[i]) train_lbl_dir else val_lbl_dir

    dest_img_path <- file.path(dest_img_dir, paste0(base_name, ".jpg"))
    dest_lbl_path <- file.path(dest_lbl_dir, paste0(base_name, ".txt"))

    dims <- dim(curr_img)
    img_w <- dims[1]
    img_h <- dims[2]

    bx_df <- NULL
    conts_list <- list()

    # --------------------------------------------------------------------------
    # MODEL 1: PerSAM (One-shot Visual Exemplar)
    # --------------------------------------------------------------------------
    if (is_persam) {
      if (isTRUE(interactive)) {
        accepted <- FALSE
        curr_sim_thresh <- sim_threshold
        curr_min_dist <- min_dist

        while (!accepted) {
          if (is.null(cached_pts) || (identical(prompt_mode, "each") && !accepted)) {
            plot(curr_img)
            title(main = sprintf("Image [%d/%d]: Click on exemplar '%s' (Right-click/Esc to detect)",
                                 i, num_imgs, class_names[1]),
                  col.main = "#0066cc", font.main = 2)
            pts <- tryCatch(graphics::locator(n = 512, type = "p", col = "#00F0FF", pch = 19), error = function(e) NULL)
            if (is.null(pts) || length(pts$x) == 0L) {
              cli::cli_alert_warning("No point clicked on image {.val {base_name}} (<Esc> pressed).")
              cat("\n")
              cli::cli_rule(left = "Options for unannotated image")
              cli::cli_bullets(c(
                "*" = "{.bold [s] / skip}: Skip this image and proceed to the next",
                "*" = "{.bold [a] / abort}: Abort process and save dataset with annotations created so far",
                "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image"
              ))
              ans_esc <- tolower(trimws(readline("Select action [s = skip / a = abort / r = retry, default: s]: ")))

              if (ans_esc == "a" || ans_esc == "abort" || ans_esc == "abortar" || ans_esc == "q" || ans_esc == "quit") {
                cli::cli_alert_info("Interactive labeling aborted by user. Finalizing and saving dataset...")
                stop_early <- TRUE
                break
              } else if (ans_esc == "r" || ans_esc == "retry" || ans_esc == "tentar") {
                cached_pts <- NULL
                next
              } else {
                cli::cli_alert_info("Image {.val {base_name}} skipped.")
                break
              }
            }
            graphics::points(pts$x, pts$y, pch = 3, col = "#FF0055", lwd = 2, cex = 1.4)
            cached_pts <- cbind(pts$x, pts$y)
          }

          ex_coords <- cached_pts

          if (isTRUE(verbose)) {
            cli::cli_alert_info("Computing PerSAM similarity map and instance segmentations (sim_threshold = {curr_sim_thresh}, min_dist = {curr_min_dist})...")
          }

          p_res <- .run_persam(
            mat = image_data(curr_img),
            exemplar_points = ex_coords,
            sim_threshold = curr_sim_thresh,
            min_dist = curr_min_dist,
            iou_threshold = iou_threshold,
            max_objects = max_objects,
            feat_res = feat_res,
            superres_map = superres_map,
            engine = engine,
            device_id = device_id,
            verbose = FALSE,
            mask = (task == "segment")
          )

          bx_df <- p_res$boxes
          conts_list <- p_res$contours

          # Filter edge artifacts
          if (isTRUE(filter_edge) && !is.null(bx_df) && nrow(bx_df) > 0L) {
            bw_px <- bx_df$xmax - bx_df$xmin
            bh_px <- bx_df$ymax - bx_df$ymin
            aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
            is_border <- (bx_df$xmin <= 4 | bx_df$xmax >= (img_w - 4) | bx_df$ymin <= 4 | bx_df$ymax >= (img_h - 4))
            keep <- !(is_border & (aspect > 3.5 | bw_px < 15 | bh_px < 15))
            bx_df <- bx_df[keep, , drop = FALSE]
            if (length(conts_list) > 0L) {
              conts_list <- conts_list[keep[seq_along(conts_list)]]
            }
          }

          n_found <- if (!is.null(bx_df)) nrow(bx_df) else 0L

          .draw_persam_feedback(
            img = curr_img,
            boxes = bx_df,
            contours = conts_list,
            similarity_map = p_res$similarity_map,
            task = task,
            class_name = class_names[1],
            show_similarity = show_similarity,
            sim_threshold = curr_sim_thresh,
            img_idx = i,
            total_imgs = num_imgs
          )

          cli::cli_alert_success("PerSAM detected {.bold {n_found}} object(s) with similarity >= {curr_sim_thresh}.")

          cat("\n")
          cli::cli_rule(left = "Interactive Labeling Options")
          cli::cli_bullets(c(
            "*" = "{.bold [Enter] / y}: Accept detections and save label",
            "*" = "{.bold r}: Retry / re-click exemplars for this image",
            "*" = "{.bold t <val>}: Adjust similarity threshold (e.g. {.code t 0.75} or {.code t 0.60})",
            "*" = "{.bold d <val>}: Adjust minimum distance between objects (e.g. {.code d 30} or {.code d 40})",
            "*" = "{.bold s}: Skip this image",
            "*" = "{.bold q}: Quit interactive labeling and save current dataset"
          ))

          ans <- tolower(trimws(readline("Select action [Enter to accept]: ")))

          if (ans == "" || ans == "y" || ans == "yes") {
            accepted <- TRUE
            if (identical(prompt_mode, "each")) {
              cached_pts <- NULL
            }
          } else if (ans == "r" || ans == "retry") {
            cached_pts <- NULL
          } else if (startsWith(ans, "t")) {
            new_val <- as.numeric(sub("^t\\s*", "", ans))
            if (!is.na(new_val) && new_val > 0 && new_val < 1) {
              curr_sim_thresh <- new_val
              cli::cli_alert_info("Updated threshold to {curr_sim_thresh}. Re-evaluating...")
            } else {
              cli::cli_alert_warning("Invalid threshold value. Must be between 0 and 1 (e.g. 't 0.55').")
            }
          } else if (startsWith(ans, "d")) {
            new_d <- as.numeric(sub("^d\\s*", "", ans))
            if (!is.na(new_d) && new_d > 0) {
              curr_min_dist <- new_d
              cli::cli_alert_info("Updated min_dist to {curr_min_dist}px. Re-evaluating...")
            } else {
              cli::cli_alert_warning("Invalid min_dist value. Must be positive (e.g. 'd 30').")
            }
          } else if (ans == "s" || ans == "skip") {
            cli::cli_alert_info("Skipped image {base_name}.")
            if (identical(prompt_mode, "each")) cached_pts <- NULL
            break
          } else if (ans == "q" || ans == "quit") {
            cli::cli_alert_info("Quitting interactive labeling. Saving dataset with existing annotations...")
            stop_early <- TRUE
            break
          }
        }

        if (stop_early) break
        if (!accepted) next

      } else {
        # Non-interactive / Batch mode
        if (is.null(cached_pts)) {
          cli::cli_abort("In non-interactive mode with {.val persam}, please provide {.arg prompt} (exemplar coordinates).")
        }

        p_res <- .run_persam(
          mat = image_data(curr_img),
          exemplar_points = cached_pts,
          sim_threshold = sim_threshold,
          min_dist = min_dist,
          iou_threshold = iou_threshold,
          max_objects = max_objects,
          feat_res = feat_res,
          superres_map = superres_map,
          engine = engine,
          device_id = device_id,
          verbose = FALSE,
          mask = (task == "segment")
        )

        bx_df <- p_res$boxes
        conts_list <- p_res$contours

        if (isTRUE(filter_edge) && !is.null(bx_df) && nrow(bx_df) > 0L) {
          bw_px <- bx_df$xmax - bx_df$xmin
          bh_px <- bx_df$ymax - bx_df$ymin
          aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
          is_border <- (bx_df$xmin <= 4 | bx_df$xmax >= (img_w - 4) | bx_df$ymin <= 4 | bx_df$ymax >= (img_h - 4))
          keep <- !(is_border & (aspect > 3.5 | bw_px < 15 | bh_px < 15))
          bx_df <- bx_df[keep, , drop = FALSE]
          if (length(conts_list) > 0L) {
            conts_list <- conts_list[keep[seq_along(conts_list)]]
          }
        }
      }

    # --------------------------------------------------------------------------
    # MODEL 2: Grounded-SAM (Open-Vocabulary Text Prompt)
    # --------------------------------------------------------------------------
    } else if (is_grounded_sam) {
      txt_prompt <- if (!is.null(prompt) && is.character(prompt)) prompt else paste(class_names, collapse = ", ")
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Detecting objects with Grounded-SAM (prompt: {.val {txt_prompt}})...")
      }
      gs_res <- .run_grounded_sam(
        mat = image_data(curr_img),
        prompt = txt_prompt,
        box_threshold = box_threshold,
        text_threshold = text_threshold,
        iou_threshold = iou_threshold,
        engine = engine,
        device_id = device_id,
        verbose = FALSE,
        mask = (task == "segment")
      )
      bx_df <- gs_res$boxes
      conts_list <- gs_res$contours

      if (isTRUE(filter_edge) && !is.null(bx_df) && nrow(bx_df) > 0L) {
        bw_px <- bx_df$xmax - bx_df$xmin
        bh_px <- bx_df$ymax - bx_df$ymin
        aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
        is_border <- (bx_df$xmin <= 4 | bx_df$xmax >= (img_w - 4) | bx_df$ymin <= 4 | bx_df$ymax >= (img_h - 4))
        keep <- !(is_border & (aspect > 4.0 | bw_px < 10 | bh_px < 10))
        bx_df <- bx_df[keep, , drop = FALSE]
        if (length(conts_list) > 0L) {
          conts_list <- conts_list[keep[seq_along(conts_list)]]
        }
      }

    # --------------------------------------------------------------------------
    # MODEL 3: Existing YOLO Model (ONNX Assisted Pre-Labeling)
    # --------------------------------------------------------------------------
    } else if (is_yolo) {
      if (identical(task, "detect")) {
        y_res <- image_detect_dl(
          curr_img,
          model = model[1],
          conf_threshold = box_threshold,
          iou_threshold = iou_threshold,
          plot = FALSE,
          verbose = FALSE,
          engine = engine,
          device_id = device_id
        )
        bx_df <- as.data.frame(y_res)
      } else {
        y_res <- image_segment_dl(
          curr_img,
          model = model[1],
          type = "highlight",
          box_threshold = box_threshold,
          iou_threshold = iou_threshold,
          plot = FALSE,
          verbose = FALSE,
          engine = engine,
          device_id = device_id
        )
        bx_df <- if (!is.null(y_res$boxes)) as.data.frame(y_res$boxes) else data.frame()
        conts_list <- if (!is.null(y_res$contours)) y_res$contours else list()
      }

    # --------------------------------------------------------------------------
    # MODEL 4: Salient Object Background Models (ben2, u2netp, etc.)
    # --------------------------------------------------------------------------
    } else if (is_salient) {
      s_res <- image_segment_dl(
        curr_img,
        model = model[1],
        type = "mask",
        plot = FALSE,
        verbose = FALSE,
        engine = engine,
        device_id = device_id
      )
      lbl_mat <- bwlabel_cpp(image_data(s_res) > 0)
      conts_list <- if (max(lbl_mat) > 0) contour(lbl_mat) else list()

    # --------------------------------------------------------------------------
    # MODEL 5: Classical Image Analysis (analyze_objects)
    # --------------------------------------------------------------------------
    } else {
      if (identical(annotations, "auto")) {
        res_obj <- analyze_objects(curr_img, plot = FALSE, ...)
        if (!is.null(res_obj[["contours"]])) {
          conts_list <- res_obj[["contours"]]
        }
      } else if (inherits(annotations, "anal_obj")) {
        if (!is.null(annotations[["contours"]])) {
          conts_list <- annotations[["contours"]]
        }
      } else if (is.list(annotations) && length(annotations) == num_imgs && inherits(annotations[[i]], "anal_obj")) {
        if (!is.null(annotations[[i]][["contours"]])) {
          conts_list <- annotations[[i]][["contours"]]
        }
      } else if (is.list(annotations) && all(sapply(annotations, is.matrix))) {
        conts_list <- annotations
      } else if (is.function(annotations)) {
        res_custom <- annotations(curr_img)
        if (inherits(res_custom, "anal_obj") && !is.null(res_custom[["contours"]])) {
          conts_list <- res_custom[["contours"]]
        } else if (is.list(res_custom)) {
          conts_list <- res_custom
        }
      }
    }

    # --------------------------------------------------------------------------
    # FORMAT EXPORT INTO YOLO LABELS (.txt)
    # --------------------------------------------------------------------------
    txt_lines <- character()

    if (identical(task, "detect")) {
      # Detection format: class_id cx cy w h
      if (!is.null(bx_df) && nrow(bx_df) > 0L) {
        xmin <- bx_df$xmin; xmax <- bx_df$xmax
        ymin <- bx_df$ymin; ymax <- bx_df$ymax
        cx <- pmin(pmax(((xmin + xmax) / 2.0) / img_w, 0.0), 1.0)
        cy <- pmin(pmax(((ymin + ymax) / 2.0) / img_h, 0.0), 1.0)
        bw <- pmin(pmax((xmax - xmin) / img_w, 0.0), 1.0)
        bh <- pmin(pmax((ymax - ymin) / img_h, 0.0), 1.0)

        for (k in seq_along(cx)) {
          if ((bw[k] * img_w * bh[k] * img_h) < min_area) next
          cid <- if (!is.null(bx_df$label) && length(bx_df$label) >= k) {
            lbl <- tolower(trimws(bx_df$label[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) m_idx - 1L else class_id[1]
          } else {
            class_id[1]
          }
          line_entry <- paste(c(cid, round(c(cx[k], cy[k], bw[k], bh[k]), 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      } else if (length(conts_list) > 0L) {
        for (k in seq_along(conts_list)) {
          c_mat <- conts_list[[k]]
          if (is.null(c_mat) || nrow(c_mat) < 3L) next
          c_mat <- as.matrix(c_mat)
          x_vals <- c_mat[, 1]; y_vals <- c_mat[, 2]
          w_obj <- max(x_vals) - min(x_vals); h_obj <- max(y_vals) - min(y_vals)
          if ((w_obj * h_obj) < min_area) next

          cx <- ((min(x_vals) + max(x_vals)) / 2.0) / img_w
          cy <- ((min(y_vals) + max(y_vals)) / 2.0) / img_h
          bw <- w_obj / img_w
          bh <- h_obj / img_h
          cid <- if (length(class_id) == length(conts_list)) class_id[k] else class_id[1]
          line_entry <- paste(c(cid, round(c(cx, cy, bw, bh), 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      }

    } else {
      # Segmentation format: class_id x1 y1 x2 y2 ... xn yn
      if (length(conts_list) > 0L) {
        for (k in seq_along(conts_list)) {
          c_mat <- conts_list[[k]]
          if (is.null(c_mat) || nrow(c_mat) < 4L) next
          c_mat <- as.matrix(c_mat)

          x_vals <- c_mat[, 1]; y_vals <- c_mat[, 2]
          w_obj <- max(x_vals) - min(x_vals); h_obj <- max(y_vals) - min(y_vals)
          if ((w_obj * h_obj) < min_area) next

          n_cur <- nrow(c_mat)
          target_n <- min(as.integer(n_points), n_cur)
          if (n_cur > target_n) {
            sampled_idx <- round(seq(1, n_cur, length.out = target_n + 1L)[-(target_n + 1L)])
            c_mat <- c_mat[sampled_idx, , drop = FALSE]
          }

          x_norm <- pmin(pmax(c_mat[, 1] / img_w, 0.0), 1.0)
          y_norm <- pmin(pmax(c_mat[, 2] / img_h, 0.0), 1.0)

          cid <- if (!is.null(bx_df$label) && length(bx_df$label) >= k) {
            lbl <- tolower(trimws(bx_df$label[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) m_idx - 1L else class_id[1]
          } else if (length(class_id) == length(conts_list)) {
            class_id[k]
          } else {
            class_id[1]
          }

          interleaved <- as.vector(rbind(x_norm, y_norm))
          line_entry <- paste(c(cid, round(interleaved, 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      } else if (!is.null(bx_df) && nrow(bx_df) > 0L) {
        # Fallback: convert 4 box corners to polygon
        for (k in seq_len(nrow(bx_df))) {
          bx <- bx_df[k, ]
          if (((bx$xmax - bx$xmin) * (bx$ymax - bx$ymin)) < min_area) next
          x1 <- pmin(pmax(bx$xmin / img_w, 0.0), 1.0); x2 <- pmin(pmax(bx$xmax / img_w, 0.0), 1.0)
          y1 <- pmin(pmax(bx$ymin / img_h, 0.0), 1.0); y2 <- pmin(pmax(bx$ymax / img_h, 0.0), 1.0)
          cid <- if (!is.null(bx$label)) {
            lbl <- tolower(trimws(bx$label))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) m_idx - 1L else class_id[1]
          } else {
            class_id[1]
          }
          poly_pts <- c(x1, y1, x2, y1, x2, y2, x1, y2)
          line_entry <- paste(c(cid, round(poly_pts, 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      }
    }

    # Save Image to train/val
    image_export(curr_img, name = paste0(base_name, ".jpg"), subfolder = dest_img_dir)

    if (is_train[i]) {
      train_images_saved <- c(train_images_saved, dest_img_path)
    } else {
      val_images_saved <- c(val_images_saved, dest_img_path)
    }

    # Write text labels
    writeLines(txt_lines, dest_lbl_path)
    total_objects_annotated <- total_objects_annotated + length(txt_lines)

    if (isTRUE(verbose)) {
      cli::cli_alert_info("Annotated {length(txt_lines)} object(s) in {.file {base_name}.jpg} [{split_name}].")
    }
  }

  # Ensure val split has at least 1 image (required by YOLO validator)
  if (length(val_images_saved) == 0L && length(train_images_saved) > 0L) {
    src_img <- train_images_saved[1]
    src_lbl <- file.path(train_lbl_dir, paste0(tools::file_path_sans_ext(basename(src_img)), ".txt"))
    val_img <- file.path(val_img_dir, basename(src_img))
    val_lbl <- file.path(val_lbl_dir, basename(src_lbl))
    file.copy(src_img, val_img, overwrite = TRUE)
    if (file.exists(src_lbl)) file.copy(src_lbl, val_lbl, overwrite = TRUE)
    val_images_saved <- c(val_images_saved, val_img)
  }

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

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "YOLO dataset successfully generated in {.path {dir}} ({total_objects_annotated} objects across {length(train_images_saved) + length(val_images_saved)} images)."
    )
  }

  out <- list(
    dir = dir,
    yaml_file = yaml_path,
    train_images = train_images_saved,
    val_images = val_images_saved,
    n_objects = total_objects_annotated
  )
  class(out) <- "yolo_dataset"
  invisible(out)
}


#' @title Preview Annotations in a YOLO Dataset
#' @name yolo_dataset_preview
#' @description
#' Loads sample images and their matching `.txt` annotation labels from a YOLO
#' dataset directory and displays them in a multi-panel grid with overlaid segmentation polygons
#' or bounding boxes. By default, randomly samples `n` images from the requested split.
#'
#' @param dir Character or `yolo_dataset` object. Path to the dataset directory (or list from [yolo_dataset_export()] / [yolo_dataset_split()]).
#' @param split Character. The subset to preview: `"train"` (default), `"val"`, or `"test"`.
#' @param pattern Optional regex pattern to filter specific image names (e.g. `"toras5\\.png"`). Defaults to `NULL`.
#' @param n Integer. Number of sample images to preview. Defaults to `4L`.
#' @param random Logical. If `TRUE` (default), samples random images from the split. If `FALSE`, picks the first `n` images.
#' @param seed Optional integer seed for reproducible random selection.
#' @param rainbow Logical. If `TRUE` (default), assigns distinct rainbow colors to each object within the image (matching review previews). If `FALSE`, colors are assigned by class.
#' @param col Border color for bounding boxes/polygons. If provided, overrides `rainbow` and default class palette.
#' @param lwd Numeric. Line width for borders. Defaults to `1.5`.
#' @param fill Fill color for polygons with transparency (e.g. `"#00FF6633"`). If `NULL`, no fill.
#' @param show_labels Logical. If `TRUE` (default), displays class name labels.
#' @param label_bg Logical. If `TRUE`, draws a subtle background badge behind each label. If `FALSE` (default), draws clean, discreet text matching review previews without blocking objects.
#' @param cex Numeric. Text size for class name labels. Defaults to `0.6`.
#' @param title Logical or character. If `TRUE` (default), displays panel titles with filename and object count.
#'
#' @return An invisible list of the sampled image file paths.
#' @export
#' @examples
#' \dontrun{
#' dataset <- yolo_dataset_export(img, dir = "dataset_cafe", model = "persam")
#' yolo_dataset_preview(dataset, n = 4, random = TRUE)
#' plot(dataset, n = 6)
#' }
yolo_dataset_preview <- function(dir = "yolo_dataset",
                                 split = c("train", "val", "test"),
                                 pattern = NULL,
                                 n = 4L,
                                 random = TRUE,
                                 seed = NULL,
                                 mfrow = NULL,
                                 rainbow = TRUE,
                                 col = NULL,
                                 lwd = 1.5,
                                 fill = NULL,
                                 show_labels = TRUE,
                                 label_bg = FALSE,
                                 cex = 0.6,
                                 title = TRUE) {
  split <- match.arg(split)

  if (is.null(dir)) {
    cli::cli_abort("The {.arg dir} argument cannot be NULL. Provide a directory path or a {.code yolo_dataset} object.")
  }

  if (is.list(dir)) {
    if (!is.null(dir$dir) && is.character(dir$dir) && nzchar(dir$dir[1])) {
      dir <- dir$dir[1]
    } else if (!is.null(dir$yaml_file) && is.character(dir$yaml_file) && nzchar(dir$yaml_file[1])) {
      dir <- dirname(dir$yaml_file[1])
    } else if (length(dir) > 0L && is.character(dir[[1]]) && nzchar(dir[[1]])) {
      dir <- dir[[1]]
    } else {
      cli::cli_abort("Could not extract a valid directory path from the provided {.code yolo_dataset} object.")
    }
  }

  if (!is.character(dir) || length(dir) == 0L || !nzchar(dir[1])) {
    cli::cli_abort("Invalid {.arg dir} argument. Expected a character directory path or a {.code yolo_dataset} object.")
  }

  raw_dir <- dir[1]
  dir <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)

  if (!dir.exists(dir)) {
    candidates <- c(
      file.path(getwd(), raw_dir),
      file.path("D:/Downloads/pliman_dl", raw_dir),
      file.path("D:/Downloads/pliman_dl", basename(raw_dir)),
      file.path(tempdir(), raw_dir),
      file.path(pliman_model_dir(), raw_dir)
    )
    for (cand in candidates) {
      if (dir.exists(cand)) {
        dir <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }

  # Split aliases (support valid, validation, training, testing)
  split_aliases <- switch(
    split,
    "train" = c("train", "training"),
    "val"   = c("val", "valid", "validation"),
    "test"  = c("test", "testing"),
    split
  )

  img_dir <- NULL
  lbl_dir <- NULL

  # 1. Try reading paths and class names from data.yaml if present
  yaml_file <- file.path(dir, "data.yaml")
  class_map <- character()
  yaml_img_dir <- NULL

  if (file.exists(yaml_file)) {
    y_lines <- readLines(yaml_file, warn = FALSE)
    name_idx <- grep("^names:", y_lines)
    if (length(name_idx) > 0L) {
      inline_val <- sub("^names:\\s*", "", y_lines[name_idx])
      if (grepl("^\\[.*\\]$", trimws(inline_val))) {
        raw_items <- gsub("[\\[\\]'\"\\s]", "", inline_val)
        class_map <- strsplit(raw_items, ",")[[1]]
      } else {
        for (l in y_lines[(name_idx + 1L):length(y_lines)]) {
          if (!grepl("^\\s*\\d+:", l)) break
          val <- trimws(sub("^\\s*\\d+:\\s*['\"]?", "", l))
          val <- sub("['\"]?\\s*$", "", val)
          class_map <- c(class_map, val)
        }
      }
    }

    # Check for split line in yaml (e.g. train: ..., val: ..., test: ...)
    for (s in split_aliases) {
      split_line <- grep(paste0("^\\s*", s, "\\s*:"), y_lines, value = TRUE)
      if (length(split_line) > 0L) {
        raw_p <- trimws(sub(paste0("^\\s*", s, "\\s*:\\s*"), "", split_line[1]))
        raw_p <- gsub("['\"]", "", raw_p)
        clean_p <- sub("^(\\.\\.[\\\\/])+", "", raw_p)
        for (cand in c(file.path(dir, raw_p), file.path(dir, clean_p), raw_p)) {
          if (dir.exists(cand)) {
            yaml_img_dir <- normalizePath(cand, winslash = "/", mustWork = FALSE)
            break
          }
        }
        if (!is.null(yaml_img_dir)) break
      }
    }
  }

  if (!is.null(yaml_img_dir)) {
    img_dir <- yaml_img_dir
    cand_lbl1 <- sub("/images$", "/labels", img_dir)
    cand_lbl2 <- sub("/images/", "/labels/", img_dir)
    for (c_lbl in c(cand_lbl1, cand_lbl2)) {
      if (dir.exists(c_lbl)) {
        lbl_dir <- c_lbl
        break
      }
    }
  }

  # 2. Check candidate directory structures (Roboflow & standard pliman)
  if (is.null(img_dir) || is.null(lbl_dir)) {
    for (s in split_aliases) {
      # Roboflow / Ultralytics: dir/<split>/images & dir/<split>/labels
      c_img_rf <- file.path(dir, s, "images")
      c_lbl_rf <- file.path(dir, s, "labels")
      if (dir.exists(c_img_rf) && dir.exists(c_lbl_rf)) {
        img_dir <- c_img_rf
        lbl_dir <- c_lbl_rf
        break
      }

      # Standard pliman: dir/images/<split> & dir/labels/<split>
      c_img_pl <- file.path(dir, "images", s)
      c_lbl_pl <- file.path(dir, "labels", s)
      if (dir.exists(c_img_pl) && dir.exists(c_lbl_pl)) {
        img_dir <- c_img_pl
        lbl_dir <- c_lbl_pl
        break
      }
    }
  }

  # 3. Check flat / unpartitioned directory structures
  if (is.null(img_dir) || is.null(lbl_dir)) {
    c_img_flat <- file.path(dir, "images")
    c_lbl_flat <- file.path(dir, "labels")
    if (dir.exists(c_img_flat) && dir.exists(c_lbl_flat)) {
      img_dir <- c_img_flat
      lbl_dir <- c_lbl_flat
    }
  }

  if (is.null(img_dir) || !dir.exists(img_dir) || is.null(lbl_dir) || !dir.exists(lbl_dir)) {
    cli::cli_abort(c(
      "Could not find dataset directories for split {.val {split}} in {.path {dir}}.",
      "i" = "Expected either:",
      "*" = "{.path {file.path(dir, split, 'images')}} and {.path {file.path(dir, split, 'labels')}} (Roboflow structure)",
      "*" = "{.path {file.path(dir, 'images', split)}} and {.path {file.path(dir, 'labels', split)}} (Standard pliman structure)"
    ))
  }

  default_palette <- c("#00FF66", "#00F0FF", "#FF0055", "#FFD700", "#9D00FF", "#FF8800", "#00B4D8", "#E63946")

  img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", ignore.case = TRUE, full.names = TRUE)
  if (length(img_files) == 0L) {
    cli::cli_abort("No images found in {.path {img_dir}}.")
  }

  if (!is.null(pattern) && is.character(pattern) && nzchar(pattern[1])) {
    img_files <- img_files[grepl(pattern[1], basename(img_files), ignore.case = TRUE)]
    if (length(img_files) == 0L) {
      cli::cli_abort("No images matching pattern {.val {pattern[1]}} found in {.path {img_dir}}.")
    }
  }

  all_stems <- tools::file_path_sans_ext(basename(img_files))
  if (any(duplicated(all_stems))) {
    dup_names <- unique(all_stems[duplicated(all_stems)])
    cli::cli_alert_warning("Found {length(dup_names)} image stem collision(s) in split {.val {split}} (e.g. {.val {head(dup_names, 3)}}). Multiple images with different extensions share the same base name, causing them to compete for the same {.file .txt} label file!")
  }

  n_samples <- min(as.integer(n), length(img_files))
  if (isTRUE(random)) {
    if (!is.null(seed)) set.seed(seed)
    sample_files <- sample(img_files, size = n_samples)
  } else {
    sample_files <- img_files[seq_len(n_samples)]
  }

  if (n_samples > 1L) {
    if (is.null(mfrow)) {
      nc <- ceiling(sqrt(n_samples))
      nr <- ceiling(n_samples / nc)
      mfrow <- c(nr, nc)
    }
    op <- graphics::par(mfrow = mfrow, mar = c(1, 1, if (isTRUE(title)) 2.2 else 1, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  for (img_path in sample_files) {
    base_name <- tools::file_path_sans_ext(basename(img_path))
    lbl_path  <- file.path(lbl_dir, paste0(base_name, ".txt"))

    same_stem_imgs <- img_files[all_stems == base_name]
    if (length(same_stem_imgs) > 1L) {
      cli::cli_alert_warning("Image {.file {basename(img_path)}} shares base name with {.file {basename(setdiff(same_stem_imgs, img_path))}} on label {.file {paste0(base_name, '.txt')}}!")
    }

    img <- image_import(img_path)

    lines <- if (file.exists(lbl_path)) readLines(lbl_path, warn = FALSE) else character(0)
    lines <- lines[nzchar(trimws(lines))]
    n_objs <- length(lines)

    p_title <- if (isTRUE(title)) {
      if (n_objs == 0L) {
        sprintf("%s [background]", basename(img_path))
      } else {
        sprintf("%s (%d obj%s)", basename(img_path), n_objs, ifelse(n_objs == 1, "", "s"))
      }
    } else if (is.character(title)) {
      title
    } else {
      NULL
    }

    plot(img, main = p_title)

    dims <- dim(img)
    img_w <- dims[1]
    img_h <- dims[2]

    if (n_objs > 0L) {
      box_colors <- if (isTRUE(rainbow) && is.null(col)) {
        grDevices::rainbow(n_objs, s = 0.85, v = 0.95)
      } else {
        NULL
      }

      obj_i <- 0L
      for (ln in lines) {
        vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
        if (length(vals) < 5L) next
        obj_i <- obj_i + 1L

        cls_id <- as.integer(vals[1])
        coords <- vals[-1]

        draw_col <- if (!is.null(col)) {
          if (length(col) == 1L) col else col[(cls_id %% length(col)) + 1L]
        } else if (!is.null(box_colors)) {
          box_colors[min(obj_i, length(box_colors))]
        } else {
          default_palette[(cls_id %% length(default_palette)) + 1L]
        }

        lbl_text <- if (length(class_map) > cls_id) class_map[cls_id + 1L] else as.character(cls_id)

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
          graphics::rect(xmin, ymin, xmax, ymax, border = draw_col, lwd = lwd)

          if (isTRUE(show_labels)) {
            if (isTRUE(label_bg)) {
              th <- graphics::strheight(lbl_text, cex = cex) * 1.2
              tw <- graphics::strwidth(lbl_text, cex = cex) * 1.2
              lbl_y1 <- max(0, ymin - th)
              graphics::rect(xmin, lbl_y1, xmin + tw, ymin, col = grDevices::adjustcolor(draw_col, alpha.f = 0.6), border = NA)
              graphics::text(xmin + tw / 2, ymin - th / 2, labels = lbl_text, col = "black", cex = cex, font = 2)
            } else {
              graphics::text(xmin, ymin - 3, labels = lbl_text, col = draw_col, cex = cex, pos = 4, font = 2)
            }
          }
        } else if (length(coords) >= 6L && length(coords) %% 2 == 0) {
          # Segmentation polygon format: x1 y1 x2 y2 ...
          xs <- coords[seq(1, length(coords), by = 2L)] * img_w
          ys <- coords[seq(2, length(coords), by = 2L)] * img_h
          poly_fill <- if (!is.null(fill)) fill else if (isTRUE(rainbow)) grDevices::adjustcolor(draw_col, alpha.f = 0.35) else NULL
          graphics::polygon(xs, ys, border = draw_col, col = poly_fill, lwd = lwd)

          if (isTRUE(show_labels)) {
            min_x <- min(xs); min_y <- min(ys)
            if (isTRUE(label_bg)) {
              th <- graphics::strheight(lbl_text, cex = cex) * 1.2
              tw <- graphics::strwidth(lbl_text, cex = cex) * 1.2
              lbl_y1 <- max(0, min_y - th)
              graphics::rect(min_x, lbl_y1, min_x + tw, min_y, col = grDevices::adjustcolor(draw_col, alpha.f = 0.6), border = NA)
              graphics::text(min_x + tw / 2, min_y - th / 2, labels = lbl_text, col = "black", cex = cex, font = 2)
            } else {
              graphics::text(min_x, min_y - 3, labels = lbl_text, col = draw_col, cex = cex, pos = 4, font = 2)
            }
          }
        }
      }
    }
  }

  invisible(sample_files)
}

#' @export
plot.yolo_dataset <- function(x, ...) {
  yolo_dataset_preview(x, ...)
}

#' @export
print.yolo_dataset <- function(x, ...) {
  cli::cli_h3("YOLO Dataset (pliman)")
  cli::cli_bullets(c(
    "*" = "Directory: {.path {x$dir}}",
    "*" = "YAML config: {.path {x$yaml_file}}",
    "*" = "Train images: {length(x$train_images)}",
    "*" = "Validation images: {length(x$val_images)}",
    "*" = "Total objects: {x$n_objects}"
  ))
  cli::cli_alert_info("Use {.code plot(x)} or {.code yolo_dataset_preview(x)} to preview annotated images.")
  invisible(x)
}


#' @title Add Negative Background Images to a YOLO Dataset
#' @name yolo_dataset_add_background
#' @description
#' Adds negative background images (images without objects and matching empty `.txt`
#' label files) to an existing YOLO dataset. This is an official Ultralytics best practice
#' to drastically reduce false positive detections. Supports local image folders,
#' vector of image file paths, an Ultralytics Platform `.ndjson` export file, or a direct URL.
#'
#' @section Downloading Datasets from Ultralytics Hub:
#' Ultralytics provides a rich repository of ready-to-use open datasets at
#' \url{https://docs.ultralytics.com/datasets/} and \url{https://platform.ultralytics.com/}.
#' To obtain an `.ndjson` file to use here:
#' \enumerate{
#'   \item Visit \url{https://docs.ultralytics.com/datasets/} or log in to the Ultralytics Platform at \url{https://platform.ultralytics.com/}.
#'   \item Browse or search for your desired dataset (e.g. African Wildlife, COCO8, SKU-110k, VisDrone, VOC).
#'   \item Open the dataset page and navigate to the **Versions** or **Overview** tab.
#'   \item Click the **Download** (Export) icon in the top header.
#'   \item Select **NDJSON** (Newline Delimited JSON) format and download the file (e.g. \code{"african-wildlife.ndjson"}).
#'   \item Pass the downloaded file path directly to \code{yolo_dataset_add_background(source = "african-wildlife.ndjson")}!
#' }
#'
#' @param dir Character or `yolo_dataset`. Path to the YOLO dataset directory (e.g. `"yolo_toras"`).
#' @param source Character. Either a path to an Ultralytics `.ndjson` file, a web URL to an `.ndjson`, a directory with images, or a vector of image paths.
#' @param n Integer. Number of background images to add. Defaults to `20L`.
#' @param train_prop Numeric. Proportion of background images assigned to the training split. Defaults to `0.80`.
#' @param prefix Character. Filename prefix for added background images. Defaults to `"bg_"`.
#' @param verbose Logical. Whether to show progress messages. Defaults to `TRUE`.
#'
#' @return Invisible list with summary of added files.
#' @export
#' @examples
#' \dontrun{
#' # Using a local NDJSON downloaded from Ultralytics Hub:
#' yolo_dataset_add_background(
#'   dir = "yolo_toras",
#'   source = "D:/Downloads/pliman_dl/african-wildlife.ndjson",
#'   n = 20
#' )
#' }
yolo_dataset_add_background <- function(dir = "yolo_dataset",
                                        source,
                                        n = 20L,
                                        train_prop = 0.80,
                                        prefix = "bg_",
                                        verbose = TRUE) {
  if (is.null(dir)) {
    cli::cli_abort("The {.arg dir} argument cannot be NULL.")
  }

  if (is.list(dir)) {
    if (!is.null(dir$dir) && is.character(dir$dir) && nzchar(dir$dir[1])) {
      dir <- dir$dir[1]
    } else if (!is.null(dir$yaml_file) && is.character(dir$yaml_file) && nzchar(dir$yaml_file[1])) {
      dir <- dirname(dir$yaml_file[1])
    } else if (length(dir) > 0L && is.character(dir[[1]]) && nzchar(dir[[1]])) {
      dir <- dir[[1]]
    }
  }

  raw_dir <- dir[1]
  dir <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)

  if (!dir.exists(dir)) {
    candidates <- c(
      file.path(getwd(), raw_dir),
      file.path("D:/Downloads/pliman_dl", raw_dir),
      file.path("D:/Downloads/pliman_dl", basename(raw_dir)),
      file.path(tempdir(), raw_dir),
      file.path(pliman_model_dir(), raw_dir)
    )
    for (cand in candidates) {
      if (dir.exists(cand)) {
        dir <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }

  train_img_dir <- file.path(dir, "images", "train")
  val_img_dir   <- file.path(dir, "images", "val")
  train_lbl_dir <- file.path(dir, "labels", "train")
  val_lbl_dir   <- file.path(dir, "labels", "val")

  for (d in c(train_img_dir, val_img_dir, train_lbl_dir, val_lbl_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  n <- as.integer(n)
  if (n <= 0L) cli::cli_abort("{.arg n} must be at least 1.")

  n_train <- max(1L, round(n * train_prop))
  if (n_train >= n && n > 1L) n_train <- n - 1L
  n_val <- n - n_train

  added_train <- character()
  added_val <- character()

  # Check if source is a URL
  if (is.character(source) && length(source) == 1L && grepl("^https?://", source, ignore.case = TRUE) && grepl("\\.ndjson$", source, ignore.case = TRUE)) {
    if (isTRUE(verbose)) cli::cli_alert_info("Downloading NDJSON file from URL {.url {source}}...")
    tmp_nd <- file.path(tempdir(), "source_download.ndjson")
    utils::download.file(source, tmp_nd, mode = "wb", quiet = !verbose)
    source <- tmp_nd
  }

  is_ndjson <- is.character(source) && length(source) == 1L && grepl("\\.ndjson$", source, ignore.case = TRUE) && file.exists(source)
  is_zip_url <- is.character(source) && length(source) == 1L && grepl("^https?://.*\\.zip$", source, ignore.case = TRUE)
  is_local_zip <- is.character(source) && length(source) == 1L && grepl("\\.zip$", source, ignore.case = TRUE) && file.exists(source)
  is_slug_name <- is.character(source) && length(source) == 1L && !file.exists(source) && !grepl("[/\\\\]", source) &&
                  grepl("^[a-z0-9_-]+$", source, ignore.case = TRUE)

  # Helper to extract background samples from a ZIP archive
  extract_from_zip <- function(zip_path) {
    z_list <- utils::unzip(zip_path, list = TRUE)
    img_entries <- z_list$Name[grepl("\\.(jpg|jpeg|png|bmp|webp)$", z_list$Name, ignore.case = TRUE)]
    if (length(img_entries) == 0L) {
      cli::cli_abort("No image files found inside archive {.path {zip_path}}.")
    }
    n_sample <- min(n, length(img_entries))
    n_tr <- max(1L, round(n_sample * train_prop))
    if (n_tr >= n_sample && n_sample > 1L) n_tr <- n_sample - 1L
    n_v <- n_sample - n_tr

    set.seed(42)
    selected_entries <- sample(img_entries, size = n_sample)

    if (isTRUE(verbose)) {
      cli::cli_alert_info("Extracting {n_sample} background images ({n_tr} train / {n_v} val)...")
      pb <- cli::cli_progress_bar(
        name = "Extracting background images",
        total = n_sample,
        format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
      )
    }

    tmp_dir <- file.path(tempdir(), paste0("pliman_bg_", as.integer(stats::runif(1, 1000, 999999))))
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(zip_path, files = selected_entries, exdir = tmp_dir)

    tr_added <- character()
    val_added <- character()

    for (k in seq_along(selected_entries)) {
      is_tr <- (k <= n_tr)
      dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
      dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir

      src_f <- file.path(tmp_dir, selected_entries[k])
      ext <- tools::file_ext(selected_entries[k])
      if (!nzchar(ext)) ext <- "jpg"
      base_fn <- sprintf("%s%04d", prefix, k)
      dest_img <- file.path(dest_img_dir, paste0(base_fn, ".", ext))
      dest_lbl <- file.path(dest_lbl_dir, paste0(base_fn, ".txt"))

      if (file.exists(src_f)) {
        file.copy(src_f, dest_img, overwrite = TRUE)
        file.create(dest_lbl)
        if (is_tr) tr_added <- c(tr_added, dest_img) else val_added <- c(val_added, dest_img)
      }
      if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
    }
    if (isTRUE(verbose)) cli::cli_progress_done(id = pb)
    unlink(tmp_dir, recursive = TRUE, force = TRUE)

    list(train = tr_added, val = val_added)
  }

  if (is_local_zip) {
    if (isTRUE(verbose)) cli::cli_alert_info("Reading local ZIP archive {.path {source}}...")
    res_z <- extract_from_zip(source)
    added_train <- res_z$train
    added_val   <- res_z$val

  } else if (is_zip_url) {
    if (isTRUE(verbose)) cli::cli_alert_info("Downloading background dataset from {.url {source}}...")
    tmp_z <- file.path(tempdir(), paste0("bg_download_", as.integer(stats::runif(1, 1000, 999999)), ".zip"))
    if (requireNamespace("curl", quietly = TRUE)) {
      curl::curl_download(source, tmp_z, quiet = !verbose)
    } else {
      utils::download.file(source, tmp_z, mode = "wb", quiet = !verbose)
    }
    res_z <- extract_from_zip(tmp_z)
    added_train <- res_z$train
    added_val   <- res_z$val
    unlink(tmp_z, force = TRUE)

  } else if (is_slug_name) {
    slug <- tolower(gsub("[^a-z0-9]+", "-", source))
    slug <- gsub("^-+|-+$", "", slug)
    zip_url <- sprintf("https://github.com/ultralytics/assets/releases/download/v0.0.0/%s.zip", slug)
    if (isTRUE(verbose)) cli::cli_alert_info("Downloading official Ultralytics dataset archive ({slug}.zip)...")
    tmp_z <- file.path(tempdir(), paste0(slug, "_", as.integer(stats::runif(1, 1000, 999999)), ".zip"))
    ok_dl <- tryCatch({
      if (requireNamespace("curl", quietly = TRUE)) {
        curl::curl_download(zip_url, tmp_z, quiet = !verbose)
      } else {
        utils::download.file(zip_url, tmp_z, mode = "wb", quiet = !verbose)
      }
      TRUE
    }, error = function(e) FALSE)

    if (isTRUE(ok_dl) && file.exists(tmp_z) && file.size(tmp_z) > 1000) {
      res_z <- extract_from_zip(tmp_z)
      added_train <- res_z$train
      added_val   <- res_z$val
      unlink(tmp_z, force = TRUE)
    } else {
      cli::cli_abort("Failed to download Ultralytics dataset {.val {source}} from {.url {zip_url}}.")
    }

  } else if (is_ndjson) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Reading Ultralytics NDJSON file {.path {source}}...")
    }
    lines <- readLines(source, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    meta <- if (length(lines) > 0L) tryCatch(jsonlite::fromJSON(lines[1]), error = function(e) list()) else list()

    # Filter image lines with URL
    records <- list()
    for (ln in lines[-1]) {
      if (grepl('"type"\\s*:\\s*"image"', ln) && grepl('"url"\\s*:', ln)) {
        rec <- tryCatch(jsonlite::fromJSON(ln), error = function(e) NULL)
        if (!is.null(rec) && !is.null(rec$url) && nzchar(rec$url)) {
          records[[length(records) + 1L]] <- rec
        }
      }
    }
    if (length(records) == 0L) {
      cli::cli_abort("No image entries found in {.path {source}}.")
    }

    n_sample <- min(n, length(records))
    indices <- round(seq(1, length(records), length.out = n_sample))
    selected <- records[indices]

    # Pre-test first URL to detect expired signed URLs
    first_url <- selected[[1]]$url
    first_ok <- FALSE
    test_h <- tryCatch({
      if (requireNamespace("curl", quietly = TRUE)) curl::curl_fetch_memory(first_url) else NULL
    }, error = function(e) NULL)
    if (!is.null(test_h) && !is.null(test_h$status_code) && test_h$status_code == 200L) {
      first_ok <- TRUE
    }

    if (!first_ok) {
      # Extract dataset slug to check official persistent GitHub archive
      slug <- NULL
      if (!is.null(meta[["url"]]) && nzchar(meta[["url"]])) {
        slug <- basename(meta[["url"]])
      } else if (!is.null(meta[["name"]]) && nzchar(meta[["name"]])) {
        slug <- tolower(gsub("[^a-z0-9]+", "-", meta[["name"]]))
        slug <- gsub("^-+|-+$", "", slug)
      } else {
        slug <- tolower(gsub("\\.ndjson$", "", basename(source), ignore.case = TRUE))
        slug <- tolower(gsub("[^a-z0-9]+", "-", slug))
        slug <- gsub("^-+|-+$", "", slug)
      }

      fallback_url <- sprintf("https://github.com/ultralytics/assets/releases/download/v0.0.0/%s.zip", slug)
      fb_check <- tryCatch({
        if (requireNamespace("curl", quietly = TRUE)) curl::curl_fetch_memory(fallback_url) else NULL
      }, error = function(e) NULL)

      if (!is.null(fb_check) && !is.null(fb_check$status_code) && fb_check$status_code == 200L) {
        if (isTRUE(verbose)) {
          cli::cli_alert_warning("Signed CDN URLs in {.path {basename(source)}} have expired (HTTP 403 Forbidden).")
          cli::cli_alert_info("Automatically falling back to persistent Ultralytics release archive ({slug}.zip)...")
        }
        tmp_z <- file.path(tempdir(), paste0(slug, "_", as.integer(stats::runif(1, 1000, 999999)), ".zip"))
        if (requireNamespace("curl", quietly = TRUE)) {
          curl::curl_download(fallback_url, tmp_z, quiet = !verbose)
        } else {
          utils::download.file(fallback_url, tmp_z, mode = "wb", quiet = !verbose)
        }
        res_z <- extract_from_zip(tmp_z)
        added_train <- res_z$train
        added_val   <- res_z$val
        unlink(tmp_z, force = TRUE)
      } else {
        cli::cli_abort(c(
          "x" = "Failed to download background images: signed URLs in {.path {basename(source)}} have expired (HTTP 403 Forbidden).",
          "i" = "Ultralytics HUB generates pre-signed URLs that expire after 14 days.",
          "*" = "To fix: re-export a fresh NDJSON from Ultralytics HUB, pass a dataset keyword like {.code source = \"african-wildlife\"}, or provide a local folder of background images."
        ))
      }

    } else {
      # Direct NDJSON URLs are valid and accessible
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Downloading {n_sample} background images ({n_train} train / {n_val} val)...")
        pb <- cli::cli_progress_bar(
          name = "Downloading background images",
          total = n_sample,
          format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
        )
      }

      for (k in seq_along(selected)) {
        rec <- selected[[k]]
        is_tr <- (k <= n_train)
        dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
        dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir

        base_fn <- sprintf("%s%04d", prefix, k)
        dest_img <- file.path(dest_img_dir, paste0(base_fn, ".jpg"))
        dest_lbl <- file.path(dest_lbl_dir, paste0(base_fn, ".txt"))

        ok <- tryCatch({
          if (requireNamespace("curl", quietly = TRUE)) {
            curl::curl_download(rec$url, dest_img, quiet = TRUE)
          } else {
            utils::download.file(rec$url, dest_img, mode = "wb", quiet = TRUE)
          }
          file.create(dest_lbl)
          TRUE
        }, error = function(e) FALSE)

        if (ok) {
          if (is_tr) added_train <- c(added_train, dest_img) else added_val <- c(added_val, dest_img)
        }
        if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
      }
      if (isTRUE(verbose)) cli::cli_progress_done(id = pb)
    }

  } else {
    # Source is local files or folder
    img_files <- character()
    if (is.character(source) && length(source) == 1L && dir.exists(source)) {
      img_files <- list.files(source, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", ignore.case = TRUE, full.names = TRUE)
    } else if (is.character(source)) {
      img_files <- source[file.exists(source)]
    }

    if (length(img_files) == 0L) {
      cli::cli_abort("No valid image files found in {.arg source}.")
    }

    n_sample <- min(n, length(img_files))
    set.seed(42)
    selected_files <- sample(img_files, size = n_sample)

    if (isTRUE(verbose)) {
      cli::cli_alert_info("Adding {n_sample} background images ({n_train} train / {n_val} val)...")
      pb <- cli::cli_progress_bar(
        name = "Copying background images",
        total = n_sample,
        format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
      )
    }

    for (k in seq_along(selected_files)) {
      is_tr <- (k <= n_train)
      dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
      dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir

      ext <- tools::file_ext(selected_files[k])
      if (!nzchar(ext)) ext <- "jpg"
      base_fn <- sprintf("%s%04d", prefix, k)
      dest_img <- file.path(dest_img_dir, paste0(base_fn, ".", ext))
      dest_lbl <- file.path(dest_lbl_dir, paste0(base_fn, ".txt"))

      file.copy(selected_files[k], dest_img, overwrite = TRUE)
      file.create(dest_lbl)

      if (is_tr) added_train <- c(added_train, dest_img) else added_val <- c(added_val, dest_img)
      if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
    }
    if (isTRUE(verbose)) cli::cli_progress_done(id = pb)
  }

  tot_added <- length(added_train) + length(added_val)
  if (tot_added == 0L) {
    cli::cli_abort("No background images could be added to {.path {dir}}.")
  }

  # Invalidate stale YOLO label caches
  unlink(file.path(dir, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir, "labels", "test.cache"), force = TRUE)

  tot_tr <- length(list.files(train_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", ignore.case = TRUE))
  tot_val <- length(list.files(val_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", ignore.case = TRUE))

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "Successfully added {tot_added} background images to {.path {dir}}."
    )
    cli::cli_bullets(c(
      "*" = "Train images: {tot_tr} (+{length(added_train)} background)",
      "*" = "Validation images: {tot_val} (+{length(added_val)} background)",
      "v" = "Matching empty .txt label files created for all background samples.",
      "v" = "Stale YOLO dataset caches cleared."
    ))
  }

  invisible(list(
    train_added = added_train,
    val_added = added_val,
    total_train = tot_tr,
    total_val = tot_val
  ))
}


#' @title Organize and Split Images into a YOLO Classification Dataset
#' @name yolo_dataset_classify
#' @aliases yolo_dataset_cls
#' @description
#' `yolo_dataset_classify()` organizes raw images structured in class subfolders (or custom
#' lists of files) into a standardized YOLO classification dataset with `train/`, `val/`,
#' and optionally `test/` splits per class. The output directory is 100% compliant with Ultralytics
#' YOLOv8, YOLO11, and YOLO26, and can be passed directly to [yolo_train()] for deep learning training.
#'
#' @details
#' The resulting directory structure conforms to the official YOLO classification format:
#' ```
#' out_dir/
#' ├── train/
#' │   ├── class_a/
#' │   ├── class_b/
#' │   └── ...
#' ├── val/
#' │   ├── class_a/
#' │   ├── class_b/
#' │   └── ...
#' └── test/ (optional)
#'     ├── class_a/
#'     ├── class_b/
#'     └── ...
#' ```
#'
#' @param src_dir Character path to the source root directory containing class subdirectories
#'   (e.g., `"flower_photos"` containing `"daisy"`, `"dandelion"`, `"roses"`, `"sunflowers"`, `"tulips"`).
#'   Non-directory files (such as `LICENSE.txt`) and hidden folders (starting with `.`) are
#'   automatically ignored.
#' @param out_dir Output directory path for the YOLO classification dataset. Defaults to `"yolo_cls_dataset"`.
#' @param split Numeric vector specifying the train / val / test split ratios.
#'   Can be of length 2 (e.g. `c(train = 0.80, val = 0.20)`) or length 3
#'   (e.g. `c(train = 0.70, val = 0.20, test = 0.10)`). Ratios are automatically normalized to sum to 1.
#' @param classes Optional character vector of specific class folder names to include. If `NULL` (default),
#'   all detected subdirectories in `src_dir` are included.
#' @param ext Character vector of valid image file extensions (default: `c("jpg", "jpeg", "png", "bmp", "tif", "tiff", "webp")`).
#' @param copy Logical. If `TRUE` (default), image files are copied to `out_dir`. If `FALSE`,
#'   files are moved using [file.rename()].
#' @param balance Logical or integer. If `TRUE`, balances the dataset by randomly sampling
#'   each class down to the count of the smallest class. If an integer `N`, samples at most `N` images per class.
#'   Default is `FALSE`.
#' @param min_images Integer threshold for displaying a warning if a class contains very few images
#'   (default: `10L`). Ultralytics recommends >= 30-100 images per class for robust generalization.
#' @param seed Optional integer for reproducible random shuffling and train/val/test partitioning (default: 42).
#' @param overwrite Logical. If `TRUE`, overwrites `out_dir` if it already exists (default: `FALSE`).
#' @param verbose Logical. If `TRUE` (default), displays progress and an informative summary table.
#'
#' @return An object of class `c("yolo_dataset_cls", "yolo_dataset", "list")` containing:
#' * `dir`: Absolute path to the created dataset directory (ready to pass to `yolo_train(data = ...)`).
#' * `classes`: Character vector of class names.
#' * `splits`: Character vector of split names created (`"train"`, `"val"`, and optionally `"test"`).
#' * `summary`: `data.frame` with the number of images per class in each split and totals.
#' * `total_images`: Total count of processed images.
#' * `split_ratio`: Named numeric vector of split ratios used.
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Organize flower photos into train (70%), val (20%), test (10%)
#' ds <- yolo_dataset_classify(
#'   src_dir = "D:/Downloads/pliman_dl/flower_photos",
#'   out_dir = "dataset_flores_yolo",
#'   split = c(train = 0.70, val = 0.20, test = 0.10),
#'   min_images = 10
#' )
#'
#' # Inspect summary table
#' ds
#'
#' # Train a YOLO classification model directly with pliman!
#' model <- yolo_train(
#'   data = ds,
#'   model = "yolo11n-cls",
#'   epochs = 50
#' )
#' }
yolo_dataset_classify <- function(src_dir,
                                  out_dir = "yolo_cls_dataset",
                                  split = c(train = 0.70, val = 0.20, test = 0.10),
                                  classes = NULL,
                                  ext = c("jpg", "jpeg", "png", "bmp", "tif", "tiff", "webp"),
                                  copy = TRUE,
                                  balance = FALSE,
                                  min_images = 10L,
                                  seed = 42,
                                  overwrite = FALSE,
                                  verbose = TRUE) {
  if (missing(src_dir) || is.null(src_dir) || !nzchar(src_dir[1])) {
    cli::cli_abort("The {.arg src_dir} argument must be provided (path to the directory containing class subfolders).")
  }

  src_dir <- normalizePath(src_dir[1], winslash = "/", mustWork = FALSE)
  if (!dir.exists(src_dir)) {
    cli::cli_abort("Source directory {.path {src_dir}} does not exist.")
  }

  # Parse splits
  if (!is.numeric(split) || length(split) < 2L || length(split) > 3L) {
    cli::cli_abort("{.arg split} must be a numeric vector of length 2 (train, val) or length 3 (train, val, test).")
  }
  if (any(split <= 0)) {
    cli::cli_abort("All values in {.arg split} must be strictly positive.")
  }

  split_names <- names(split)
  if (is.null(split_names) || any(!nzchar(split_names))) {
    split_names <- if (length(split) == 2L) c("train", "val") else c("train", "val", "test")
  }
  split_norm <- split / sum(split)
  names(split_norm) <- split_names

  # Identify class subdirectories
  all_subdirs <- list.dirs(src_dir, full.names = FALSE, recursive = FALSE)
  all_subdirs <- all_subdirs[nzchar(all_subdirs) & !startsWith(all_subdirs, ".")]

  if (length(all_subdirs) == 0L) {
    cli::cli_abort("No subdirectories found in {.path {src_dir}}. Each class must have its own folder.")
  }

  if (!is.null(classes)) {
    cls_found <- intersect(classes, all_subdirs)
    if (length(cls_found) == 0L) {
      cli::cli_abort("None of the specified {.arg classes} were found in {.path {src_dir}}.")
    }
    all_subdirs <- cls_found
  }

  # Setup output directory
  out_dir_norm <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)
  if (dir.exists(out_dir_norm)) {
    existing_files <- list.files(out_dir_norm, all.files = TRUE, no.. = TRUE)
    if (length(existing_files) > 0L) {
      if (isTRUE(overwrite)) {
        unlink(out_dir_norm, recursive = TRUE, force = TRUE)
        dir.create(out_dir_norm, recursive = TRUE, showWarnings = FALSE)
      } else {
        cli::cli_abort("Output directory {.path {out_dir_norm}} already exists and is not empty. Set {.code overwrite = TRUE} to replace it.")
      }
    }
  } else {
    dir.create(out_dir_norm, recursive = TRUE, showWarnings = FALSE)
  }

  ext_pat <- paste0("\\.(", paste(tolower(ext), collapse = "|"), ")$")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Discover files per class
  class_files <- list()
  for (cls in all_subdirs) {
    c_path <- file.path(src_dir, cls)
    f_list <- list.files(c_path, pattern = ext_pat, full.names = TRUE, ignore.case = TRUE)
    if (length(f_list) > 0L) {
      class_files[[cls]] <- f_list
    }
  }

  if (length(class_files) < 2L) {
    cli::cli_abort("Found fewer than 2 valid classes with images in {.path {src_dir}}. Classification requires at least 2 classes.")
  }

  active_classes <- names(class_files)

  # Handle class balancing if requested
  if (isTRUE(balance)) {
    min_cnt <- min(vapply(class_files, length, integer(1)))
    for (cls in active_classes) {
      class_files[[cls]] <- sample(class_files[[cls]], size = min_cnt)
    }
  } else if (is.numeric(balance) && length(balance) == 1L && balance > 0) {
    bal_cnt <- as.integer(balance)
    for (cls in active_classes) {
      if (length(class_files[[cls]]) > bal_cnt) {
        class_files[[cls]] <- sample(class_files[[cls]], size = bal_cnt)
      }
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Classification Dataset Generator (pliman)")
    cli::cli_alert_info("Source directory: {.path {src_dir}}")
    cli::cli_alert_info("Destination directory: {.path {out_dir_norm}}")
    cli::cli_alert_info("Found {length(active_classes)} classes: {.val {active_classes}}")
    cli::cli_alert_info(
      "Splits: {paste(paste0(split_names, ' (', round(split_norm * 100, 1), '%)'), collapse = ' | ')}"
    )
    pb <- cli::cli_progress_bar(
      name = if (isTRUE(copy)) "Copying images" else "Moving images",
      total = sum(vapply(class_files, length, integer(1))),
      format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
    )
  }

  summary_rows <- list()
  total_moved <- 0L

  for (cls in active_classes) {
    fl <- sample(class_files[[cls]]) # shuffle
    n_tot <- length(fl)

    # Compute partition indices
    if (length(split_norm) == 2L) {
      n_tr <- max(1L, round(n_tot * split_norm[1]))
      n_tr <- min(n_tr, n_tot - 1L)
      part_list <- list(
        train = fl[seq_len(n_tr)],
        val   = fl[(n_tr + 1L):n_tot]
      )
      names(part_list) <- split_names
    } else {
      n_tr <- max(1L, round(n_tot * split_norm[1]))
      n_val <- max(1L, round(n_tot * split_norm[2]))
      if (n_tr + n_val >= n_tot) {
        n_tr <- max(1L, n_tot - 2L)
        n_val <- 1L
      }
      part_list <- list(
        train = fl[seq_len(n_tr)],
        val   = fl[(n_tr + 1L):(n_tr + n_val)],
        test  = fl[(n_tr + n_val + 1L):n_tot]
      )
      names(part_list) <- split_names
    }

    row_data <- list(class = cls)
    for (s_name in split_names) {
      dest_dir <- file.path(out_dir_norm, s_name, cls)
      if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

      cur_files <- part_list[[s_name]]
      row_data[[s_name]] <- length(cur_files)

      if (length(cur_files) > 0L) {
        dest_paths <- file.path(dest_dir, basename(cur_files))
        if (isTRUE(copy)) {
          file.copy(cur_files, dest_paths, overwrite = TRUE)
        } else {
          file.rename(cur_files, dest_paths)
        }
        total_moved <- total_moved + length(cur_files)
        if (isTRUE(verbose)) cli::cli_progress_update(id = pb, inc = length(cur_files))
      }
    }
    row_data$total <- n_tot
    summary_rows[[cls]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
  }

  if (isTRUE(verbose)) cli::cli_progress_done(id = pb)

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  # Add Total row
  tot_row <- list(class = "TOTAL")
  for (nm in setdiff(names(summary_df), "class")) {
    tot_row[[nm]] <- sum(summary_df[[nm]])
  }
  summary_df_full <- rbind(summary_df, as.data.frame(tot_row, stringsAsFactors = FALSE))

  out <- list(
    dir = out_dir_norm,
    classes = active_classes,
    splits = split_names,
    summary = summary_df_full,
    total_images = total_moved,
    split_ratio = split_norm,
    min_images = min_images
  )
  class(out) <- c("yolo_dataset_cls", "yolo_dataset", "list")

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "YOLO classification dataset successfully created in {.path {out_dir_norm}} ({total_moved} images across {length(active_classes)} classes)."
    )
    cli::cli_rule(left = "Dataset Summary")
    print(summary_df_full, row.names = FALSE)
    cli::cli_rule()

    # Health Check 1: Classes with very few total images
    low_classes <- summary_df[summary_df$total < min_images, ]
    if (nrow(low_classes) > 0L) {
      for (i in seq_len(nrow(low_classes))) {
        cli::cli_alert_warning(
          "Class {.val {low_classes$class[i]}} has only {low_classes$total[i]} image(s) (< {min_images}). Model accuracy for this class may be limited. (Recommended: >= 30-100 images/class)."
        )
      }
    }

    # Health Check 2: Validation split with critically few images
    if ("val" %in% names(summary_df)) {
      low_val <- summary_df[summary_df$val < 3L, ]
      if (nrow(low_val) > 0L) {
        for (i in seq_len(nrow(low_val))) {
          cli::cli_alert_warning(
            "Validation split for class {.val {low_val$class[i]}} has only {low_val$val[i]} image(s). Validation metrics may be unstable."
          )
        }
      }
    }

    # Health Check 3: Severe class imbalance
    if (nrow(summary_df) >= 2L) {
      max_row <- summary_df[which.max(summary_df$total), ]
      min_row <- summary_df[which.min(summary_df$total), ]
      imb_ratio <- max_row$total / max(1L, min_row$total)
      if (imb_ratio >= 3.0) {
        cli::cli_alert_info(
          "Class imbalance detected: {.val {max_row$class}} ({max_row$total} imgs) has {round(imb_ratio, 1)}x more samples than {.val {min_row$class}} ({min_row$total} imgs). Tip: Use {.code balance = TRUE} to balance sample sizes."
        )
      }
    }

    cli::cli_alert_info("Ready to train! Run: {.code yolo_train(data = '{out_dir_norm}', model = 'yolo11n-cls', epochs = 50)}")
  }

  invisible(out)
}

#' @export
yolo_dataset_cls <- yolo_dataset_classify

#' @export
print.yolo_dataset_cls <- function(x, ...) {
  cli::cli_h3("YOLO Classification Dataset (pliman)")
  cli::cli_bullets(c(
    "*" = "Directory: {.path {x$dir}}",
    "*" = "Classes ({length(x$classes)}): {.val {x$classes}}",
    "*" = "Splits: {paste(x$splits, collapse = ', ')}",
    "*" = "Total images: {x$total_images}"
  ))
  if (!is.null(x$summary)) {
    cli::cli_rule(left = "Split Breakdown")
    print(x$summary, row.names = FALSE)
    cli::cli_rule()

    # Alert if any class has low sample counts
    sub_df <- x$summary[x$summary$class != "TOTAL", , drop = FALSE]
    min_thresh <- if (!is.null(x$min_images)) x$min_images else 10L
    low_classes <- sub_df[sub_df$total < min_thresh, ]
    if (nrow(low_classes) > 0L) {
      cli::cli_alert_warning(
        "Low sample size in class(es): {paste0(low_classes$class, ' (', low_classes$total, ')', collapse = ', ')}."
      )
    }
  }
  cli::cli_alert_info("Train with: {.code yolo_train(data = '{x$dir}', model = 'yolo11n-cls', epochs = 50)}")
  invisible(x)
}

#' @export
plot.yolo_dataset_cls <- function(x, n_per_class = 2, max_images = 12, ...) {
  train_dir <- file.path(x$dir, "train")
  if (!dir.exists(train_dir)) {
    train_dir <- file.path(x$dir, x$splits[1])
  }

  sample_imgs <- list()
  sample_labels <- character()

  for (cls in x$classes) {
    cls_folder <- file.path(train_dir, cls)
    if (dir.exists(cls_folder)) {
      files <- list.files(cls_folder, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
      if (length(files) > 0L) {
        pick <- files[seq_len(min(length(files), n_per_class))]
        for (f in pick) {
          sample_imgs[[length(sample_imgs) + 1L]] <- f
          sample_labels <- c(sample_labels, cls)
        }
      }
    }
  }

  if (length(sample_imgs) == 0L) {
    cli::cli_warn("No images found to display.")
    return(invisible())
  }

  if (length(sample_imgs) > max_images) {
    sample_imgs <- sample_imgs[seq_len(max_images)]
    sample_labels <- sample_labels[seq_len(max_images)]
  }

  imgs_loaded <- lapply(sample_imgs, function(f) try(image_import(f), silent = TRUE))
  valid_idx <- vapply(imgs_loaded, function(im) !inherits(im, "try-error") && is_image(im), logical(1))

  if (!any(valid_idx)) {
    cli::cli_warn("Could not load sample images for plotting.")
    return(invisible())
  }

  imgs_loaded <- imgs_loaded[valid_idx]
  sample_labels <- sample_labels[valid_idx]

  n_plots <- length(imgs_loaded)
  nc <- ceiling(sqrt(n_plots))
  nr <- ceiling(n_plots / nc)

  op <- graphics::par(mfrow = c(nr, nc), mar = c(1, 1, 2.5, 1))
  on.exit(graphics::par(op), add = TRUE)

  for (i in seq_along(imgs_loaded)) {
    plot(imgs_loaded[[i]], main = sample_labels[i], axes = FALSE)
  }
  invisible(x)
}


#' @title Import an Ultralytics NDJSON Dataset and Export to YOLO Structure
#' @name yolo_dataset_from_ndjson
#' @aliases yolo_dataset_import_ndjson
#' @description
#' Reads an Ultralytics Hub dataset export (`.ndjson` format), parses dataset metadata,
#' classes, bounding boxes or segmentation polygons, downloads images before their
#' signed URLs expire, and writes a complete standard YOLO dataset folder
#' with `images/`, `labels/`, and `data.yaml`. Optionally packages the entire dataset
#' into a `.zip` archive for permanent offline storage and distribution.
#'
#' @section Downloading Datasets from Ultralytics Hub:
#' Ultralytics provides a rich repository of ready-to-use open datasets at
#' \url{https://docs.ultralytics.com/datasets/} and \url{https://platform.ultralytics.com/}.
#' To obtain an `.ndjson` file to use here:
#' \enumerate{
#'   \item Visit \url{https://docs.ultralytics.com/datasets/} or log in to the Ultralytics Platform at \url{https://platform.ultralytics.com/}.
#'   \item Browse or search for your desired dataset (e.g. African Wildlife, COCO8, SKU-110k, VisDrone, VOC).
#'   \item Open the dataset page and navigate to the **Versions** or **Overview** tab.
#'   \item Click the **Download** (Export) icon in the top header.
#'   \item Select **NDJSON** (Newline Delimited JSON) format and download the file (e.g. \code{"african-wildlife.ndjson"}).
#'   \item Run \code{yolo_dataset_from_ndjson(file = "african-wildlife.ndjson", zip = TRUE)} to download all
#'     images and generate the complete YOLO dataset before the cloud URLs expire!
#' }
#'
#' @param file Character. Path to the `.ndjson` file (e.g. `"african-wildlife.ndjson"`), or a direct URL to an `.ndjson` file.
#' @param dir Character or `NULL`. Target directory for the exported YOLO dataset. If `NULL` (default),
#'   creates a directory named `yolo_<dataset_name>` in the same folder as `file`.
#' @param splits Character vector. Which splits to process (`"train"`, `"val"`, `"test"`).
#'   Defaults to all available splits in the file.
#' @param train_prop Numeric or `NULL`. If specified (e.g. `0.80`), re-splits the images randomly into
#'   train and validation subsets instead of using the original splits from the `.ndjson`. Defaults to `NULL`.
#' @param max_images Integer or `NULL`. Maximum number of images to download. Defaults to `NULL` (all images).
#' @param download_images Logical. Whether to download image files. Defaults to `TRUE`.
#' @param parallel Logical. Whether to download images concurrently using `curl::multi_download`. Defaults to `TRUE`.
#' @param zip Logical. Whether to create a `.zip` archive of the exported dataset directory. Defaults to `FALSE`.
#' @param overwrite Logical. Whether to overwrite existing files. Defaults to `FALSE` (enables resuming).
#' @param verbose Logical. Whether to show progress bar and informative messages. Defaults to `TRUE`.
#'
#' @return An object of class `yolo_dataset` with dataset directory paths and classes.
#' @export
#' @examples
#' \dontrun{
#' # Download and export the entire dataset permanently before URLs expire:
#' ds <- yolo_dataset_from_ndjson(
#'   file = "D:/Downloads/pliman_dl/african-wildlife.ndjson",
#'   zip = TRUE
#' )
#' yolo_dataset_preview(ds)
#' }
yolo_dataset_from_ndjson <- function(file,
                                     dir = NULL,
                                     splits = c("train", "val", "test"),
                                     train_prop = NULL,
                                     max_images = NULL,
                                     download_images = TRUE,
                                     parallel = TRUE,
                                     zip = FALSE,
                                     overwrite = FALSE,
                                     verbose = TRUE) {
  # Handle remote URL
  if (is.character(file) && length(file) == 1L && grepl("^https?://", file, ignore.case = TRUE)) {
    if (isTRUE(verbose)) cli::cli_alert_info("Downloading NDJSON file from URL {.url {file}}...")
    tmp_nd <- file.path(tempdir(), "file_download.ndjson")
    utils::download.file(file, tmp_nd, mode = "wb", quiet = !verbose)
    file <- tmp_nd
  }

  if (!file.exists(file)) cli::cli_abort("File {.path {file}} not found.")

  if (isTRUE(verbose)) cli::cli_alert_info("Reading NDJSON dataset {.path {file}}...")
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) cli::cli_abort("NDJSON file is empty or invalid.")

  meta <- tryCatch(jsonlite::fromJSON(lines[1]), error = function(e) list())
  dataset_name <- if (!is.null(meta$name) && nzchar(meta$name)) meta$name else "custom_dataset"
  clean_name <- tolower(gsub("[^A-Za-z0-9]+", "_", dataset_name))
  clean_name <- gsub("^_|_$", "", clean_name)

  if (is.null(dir)) {
    dir <- file.path(dirname(file), paste0("yolo_", clean_name))
  }
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)

  task <- if (!is.null(meta$task) && nzchar(meta$task)) meta$task else "detect"
  class_names <- meta$class_names
  names_vec <- if (is.list(class_names)) unlist(class_names) else as.character(class_names)
  if (is.null(names(names_vec)) || !any(nzchar(names(names_vec)))) {
    names(names_vec) <- as.character(seq_along(names_vec) - 1L)
  }

  records <- list()
  for (i in 2:length(lines)) {
    rec <- tryCatch(jsonlite::fromJSON(lines[i]), error = function(e) NULL)
    if (!is.null(rec) && !is.null(rec$url) && nzchar(rec$url)) {
      sp <- if (!is.null(rec$split) && nzchar(rec$split)) rec$split else "train"
      if (sp %in% splits) {
        rec$split <- sp
        records[[length(records) + 1L]] <- rec
      }
    }
  }

  if (length(records) == 0L) {
    cli::cli_abort("No valid image entries found for splits: {paste(splits, collapse = ', ')}.")
  }

  if (!is.null(max_images)) {
    records <- records[seq_len(min(as.integer(max_images), length(records)))]
  }

  # Optional re-splitting into train/val
  if (is.numeric(train_prop) && length(train_prop) == 1L && train_prop > 0 && train_prop < 1) {
    n_tot <- length(records)
    n_tr <- max(1L, round(n_tot * train_prop))
    set.seed(42)
    shuf_idx <- sample(n_tot)
    for (k in seq_along(records)) {
      records[[k]]$split <- if (shuf_idx[k] <= n_tr) "train" else "val"
    }
  }

  active_splits <- unique(sapply(records, function(r) r$split))
  for (s in active_splits) {
    dir.create(file.path(dir, "images", s), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir, "labels", s), recursive = TRUE, showWarnings = FALSE)
  }

  files <- sapply(records, function(r) {
    fn <- if (!is.null(r$file) && nzchar(r$file)) r$file else basename(strsplit(r$url, "\\?")[[1]][1])
    gsub("[ ()]", "_", fn)
  })
  splits_vec <- sapply(records, function(r) r$split)
  urls_vec <- sapply(records, function(r) r$url)

  dest_imgs <- file.path(dir, "images", splits_vec, files)
  dest_lbls <- file.path(dir, "labels", splits_vec, paste0(tools::file_path_sans_ext(files), ".txt"))

  if (isTRUE(download_images)) {
    need_dl <- if (isTRUE(overwrite)) rep(TRUE, length(dest_imgs)) else !file.exists(dest_imgs)
    n_need <- sum(need_dl)
    if (n_need > 0L) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Downloading {n_need} images to {.path {dir}}...")
        pb <- cli::cli_progress_bar(
          name = "Downloading images",
          total = n_need,
          format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
        )
      }
      dl_idx <- which(need_dl)
      if (isTRUE(parallel) && requireNamespace("curl", quietly = TRUE)) {
        batch_size <- 25L
        batches <- split(dl_idx, ceiling(seq_along(dl_idx) / batch_size))
        for (b in batches) {
          curl::multi_download(urls_vec[b], dest_imgs[b], progress = FALSE)
          if (isTRUE(verbose)) cli::cli_progress_update(id = pb, inc = length(b))
        }
      } else {
        for (idx in dl_idx) {
          tryCatch({
            if (requireNamespace("curl", quietly = TRUE)) {
              curl::curl_download(urls_vec[idx], dest_imgs[idx], quiet = TRUE)
            } else {
              utils::download.file(urls_vec[idx], dest_imgs[idx], mode = "wb", quiet = TRUE)
            }
          }, error = function(e) NULL)
          if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
        }
      }
      if (isTRUE(verbose)) cli::cli_progress_done(id = pb)
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_progress_step("Writing YOLO annotation label files...", msg_done = "Labels created")
  }
  tot_ann <- 0L
  for (i in seq_along(records)) {
    r <- records[[i]]
    lbl_f <- dest_lbls[i]
    lbl_lines <- character()
    if (!is.null(r$annotations$boxes) && length(r$annotations$boxes) > 0) {
      bx <- as.matrix(r$annotations$boxes)
      tot_ann <- tot_ann + nrow(bx)
      lbl_lines <- apply(bx, 1, function(row) {
        sprintf("%d %.6f %.6f %.6f %.6f", as.integer(row[1]), row[2], row[3], row[4], row[5])
      })
    } else if (!is.null(r$annotations$segments) && length(r$annotations$segments) > 0) {
      seg <- r$annotations$segments
      if (is.list(seg)) {
        for (s_item in seg) {
          tot_ann <- tot_ann + 1L
          lbl_lines <- c(lbl_lines, paste(as.character(s_item), collapse = " "))
        }
      } else if (is.matrix(seg)) {
        tot_ann <- tot_ann + nrow(seg)
        lbl_lines <- apply(seg, 1, function(row) paste(row, collapse = " "))
      }
    }
    if (length(lbl_lines) > 0) writeLines(lbl_lines, lbl_f) else file.create(lbl_f)
  }
  if (isTRUE(verbose)) cli::cli_progress_done()

  yaml_file <- file.path(dir, "data.yaml")
  yaml_lines <- c(
    paste0("path: ", dir),
    if (dir.exists(file.path(dir, "images", "train"))) "train: images/train" else NULL,
    if (dir.exists(file.path(dir, "images", "val"))) "val: images/val" else NULL,
    if (dir.exists(file.path(dir, "images", "test"))) "test: images/test" else NULL,
    "names:"
  )
  for (nm_idx in names(names_vec)) {
    yaml_lines <- c(yaml_lines, sprintf("  %s: %s", nm_idx, names_vec[[nm_idx]]))
  }
  writeLines(yaml_lines, yaml_file)

  zip_file <- NULL
  if (isTRUE(zip)) {
    if (isTRUE(verbose)) cli::cli_progress_step("Creating ZIP package...", msg_done = "ZIP created")
    zip_file <- paste0(dir, ".zip")
    orig_wd <- getwd()
    setwd(dirname(dir))
    utils::zip(zipfile = basename(zip_file), files = basename(dir))
    setwd(orig_wd)
    if (isTRUE(verbose)) cli::cli_progress_done()
  }

  tr_imgs <- list.files(file.path(dir, "images", "train"), full.names = TRUE)
  v_imgs  <- list.files(file.path(dir, "images", "val"), full.names = TRUE)
  te_imgs <- list.files(file.path(dir, "images", "test"), full.names = TRUE)

  res <- structure(
    list(
      dir = dir,
      yaml_file = yaml_file,
      zip_file = zip_file,
      task = task,
      train_images = tr_imgs,
      val_images = v_imgs,
      test_images = te_imgs,
      classes = names_vec,
      n_objects = tot_ann
    ),
    class = "yolo_dataset"
  )

  if (isTRUE(verbose)) {
    cli::cli_alert_success("YOLO dataset successfully created at {.path {dir}}!")
    cli::cli_bullets(c(
      "*" = "Task: {task}",
      "*" = "Classes: {paste(names_vec, collapse = ', ')}",
      "*" = "Train images: {length(tr_imgs)}",
      "*" = "Validation images: {length(v_imgs)}",
      "*" = if (length(te_imgs) > 0) "Test images: {length(te_imgs)}" else NULL,
      "*" = "Total objects: {tot_ann}",
      "v" = if (!is.null(zip_file)) "ZIP package created: {.path {zip_file}}" else NULL
    ))
    cli::cli_alert_info("Images and annotations are now permanently stored locally (offline, safe from URL expiration).")
  }

  invisible(res)
}


#' @rdname yolo_dataset_from_ndjson
#' @export
yolo_dataset_import_ndjson <- yolo_dataset_from_ndjson


# ==============================================================================
# GPU, PYTHON & PRE-TRAINED WEIGHTS HELPERS
# ==============================================================================

# Internal helper to detect NVIDIA GPU presence, model, and driver
.detect_nvidia_gpu <- function() {
  smi <- Sys.which("nvidia-smi")
  if (!nzchar(smi) && .Platform$OS.type == "windows") {
    cand <- file.path(Sys.getenv("SystemRoot", "C:/Windows"), "System32", "nvidia-smi.exe")
    if (file.exists(cand)) smi <- cand
  }
  if (nzchar(smi)) {
    out <- suppressWarnings(tryCatch(
      system2(smi, args = c("--query-gpu=name,driver_version", "--format=csv,noheader"), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    ))
    if (length(out) > 0L && nzchar(trimws(out[1]))) {
      parts <- strsplit(trimws(out[1]), ",")[[1]]
      gpu_name <- trimws(parts[1])
      driver_ver <- if (length(parts) > 1L) trimws(parts[2]) else ""
      return(list(has_gpu = TRUE, name = gpu_name, driver = driver_ver))
    }
  }
  list(has_gpu = FALSE, name = NULL, driver = NULL)
}

# Internal helper to locate an available Python executable
.find_python_exec <- function() {
  # 1. System PATH
  for (cmd in c("python", "python3", "py")) {
    p <- Sys.which(cmd)
    if (nzchar(p)) return(normalizePath(p, winslash = "/"))
  }

  # 2. Windows standard paths
  if (.Platform$OS.type == "windows") {
    cand_dirs <- c(
      file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "Python", "Python312", "python.exe"),
      file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "Python", "Python311", "python.exe"),
      file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "Python", "Python310", "python.exe"),
      file.path(Sys.getenv("ProgramFiles"), "Python312", "python.exe"),
      file.path(Sys.getenv("ProgramFiles"), "Python311", "python.exe"),
      file.path(Sys.getenv("LOCALAPPDATA"), "Microsoft", "WindowsApps", "python.exe")
    )
    for (cand in cand_dirs) {
      if (file.exists(cand)) return(normalizePath(cand, winslash = "/"))
    }
  } else {
    # 3. Unix standard paths (Linux / macOS)
    cand_dirs <- c(
      "/usr/bin/python3",
      "/usr/local/bin/python3",
      "/opt/homebrew/bin/python3",
      path.expand("~/.pyenv/shims/python3"),
      path.expand("~/miniconda3/bin/python3"),
      path.expand("~/anaconda3/bin/python3")
    )
    for (cand in cand_dirs) {
      if (file.exists(cand)) return(normalizePath(cand, winslash = "/"))
    }
  }

  character()
}

# Internal helper to verify or install Python seamlessly
.ensure_python <- function(auto_install = TRUE) {
  py_exec <- .find_python_exec()
  if (length(py_exec) > 0L && nzchar(py_exec)) {
    return(py_exec)
  }

  cli::cli_rule("{.strong pliman: YOLO Training Environment}")
  cli::cli_alert_warning("Python is required for model backpropagation training, but was not found.")
  cli::cli_alert_info("Inference with pre-trained models in pliman is 100% C++ (zero-Python required).")

  do_install <- FALSE
  if (interactive()) {
    ans <- utils::askYesNo("Would you like pliman to automatically download and install Python now?", default = TRUE)
    do_install <- isTRUE(ans)
  } else {
    do_install <- isTRUE(auto_install)
  }

  if (!do_install) {
    cli::cli_abort(
      c(
        "x" = "Python is required for training.",
        "i" = "Install Python 3.9+ from https://www.python.org or run with {.code auto_install = TRUE}."
      )
    )
  }

  cli::cli_alert_info("Starting automatic Python installation...")

  if (.Platform$OS.type == "windows") {
    winget <- Sys.which("winget")
    winget_done <- FALSE
    if (nzchar(winget)) {
      cli::cli_alert_info("Installing Python 3.12 via Windows Package Manager (winget)...")
      res <- suppressWarnings(system2(winget, args = c("install", "--id", "Python.Python.3.12", "-e", "--silent", "--accept-package-agreements", "--accept-source-agreements")))
      winget_done <- (res == 0L)
    }

    if (!winget_done) {
      cli::cli_alert_info("Downloading Python 3.12 installer from python.org...")
      installer_url <- "https://www.python.org/ftp/python/3.12.8/python-3.12.8-amd64.exe"
      dest_exe <- file.path(tempdir(), "python_setup.exe")
      dl_ok <- .download_url_robust(installer_url, dest_exe, min_size = 10000000)
      if (dl_ok) {
        cli::cli_alert_info("Running silent Python installer...")
        system2(dest_exe, args = c("/quiet", "InstallAllUsers=0", "PrependPath=1", "Include_pip=1", "Include_test=0"), wait = TRUE)
        unlink(dest_exe)
      }
    }
  } else {
    cli::cli_alert_info("On Linux/macOS, please ensure Python 3.9+ is installed (e.g. via apt-get, dnf, or brew install python3).")
  }

  py_exec <- .find_python_exec()
  if (length(py_exec) == 0L || !nzchar(py_exec)) {
    cli::cli_abort(
      c(
        "x" = "Automatic Python installation could not be completed.",
        "i" = "Please install Python 3.9+ from https://www.python.org and ensure it is added to your PATH."
      )
    )
  }

  cli::cli_alert_success("Python installed and located at {.path {py_exec}}")
  py_exec
}

# Internal helper to configure PyTorch (CUDA / CPU), Ultralytics, and ONNX
.setup_python_yolo_env <- function(py_exec, auto_install = TRUE) {
  gpu_info <- .detect_nvidia_gpu()

  chk_py <- paste0(
    "import sys\n",
    "try:\n",
    "    import torch\n",
    "    t_ver = torch.__version__\n",
    "    cuda_ok = torch.cuda.is_available()\n",
    "except ImportError:\n",
    "    t_ver = 'none'\n",
    "    cuda_ok = False\n",
    "try:\n",
    "    import ultralytics\n",
    "    u_ok = True\n",
    "except ImportError:\n",
    "    u_ok = False\n",
    "try:\n",
    "    import onnx\n",
    "    o_ok = True\n",
    "except ImportError:\n",
    "    o_ok = False\n",
    "print(f'{t_ver}|{cuda_ok}|{u_ok}|{o_ok}')\n"
  )

  status_out <- suppressWarnings(tryCatch(
    system2(py_exec, args = c("-c", shQuote(chk_py)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  ))

  t_ver <- "none"
  cuda_ok <- FALSE
  u_ok <- FALSE
  o_ok <- FALSE

  if (length(status_out) > 0L && grepl("\\|", status_out[length(status_out)])) {
    parts <- strsplit(trimws(status_out[length(status_out)]), "\\|")[[1]]
    t_ver <- parts[1]
    cuda_ok <- identical(tolower(parts[2]), "true")
    u_ok <- identical(tolower(parts[3]), "true")
    o_ok <- identical(tolower(parts[4]), "true")
  }

  needs_torch_cuda_upgrade <- (gpu_info$has_gpu && (!cuda_ok || grepl("cpu", t_ver, ignore.case = TRUE) || t_ver == "none"))
  needs_torch_install <- (t_ver == "none")

  if (needs_torch_cuda_upgrade || needs_torch_install || !u_ok || !o_ok) {
    if (!isTRUE(auto_install)) {
      cli::cli_abort("Required Python dependencies are missing. Set {.code auto_install = TRUE} or install them via pip.")
    }

    if (gpu_info$has_gpu && (needs_torch_cuda_upgrade || needs_torch_install)) {
      cli::cli_alert_success("NVIDIA GPU detected: {.strong {gpu_info$name}} (CUDA capable) \U0001f4a5")
      cli::cli_alert_info("Configuring PyTorch with CUDA 12.4 acceleration support...")
      system2(py_exec, args = c("-m", "pip", "install", "--upgrade", "torch", "torchvision", "--index-url", "https://download.pytorch.org/whl/cu124"))
      cuda_ok <- TRUE
    } else if (needs_torch_install) {
      cli::cli_alert_info("No NVIDIA GPU detected. Installing CPU-optimized PyTorch...")
      system2(py_exec, args = c("-m", "pip", "install", "--upgrade", "torch", "torchvision"))
    }

    if (!u_ok) {
      cli::cli_alert_info("Installing {.pkg ultralytics}...")
      system2(py_exec, args = c("-m", "pip", "install", "ultralytics"))
    }

    if (!o_ok) {
      cli::cli_alert_info("Installing {.pkg onnx} and {.pkg onnxslim}...")
      system2(py_exec, args = c("-m", "pip", "install", "onnx", "onnxslim"))
    }

    cli::cli_alert_success("Python training dependencies verified!")
  }

  list(gpu_info = gpu_info, cuda_available = cuda_ok)
}

# Internal helper to locate or download pre-trained .pt model weights from NEPEM GitHub
.resolve_yolo_pt <- function(model, output_dir = pliman_model_dir()) {
  if (file.exists(model[1])) {
    return(normalizePath(model[1], winslash = "/"))
  }

  pt_file <- model[1]
  if (!grepl("\\.(pt|yaml)$", pt_file, ignore.case = TRUE)) {
    resolved <- .resolve_model_name(pt_file)
    pt_file <- paste0(resolved, ".pt")
  } else if (grepl("\\.pt$", pt_file, ignore.case = TRUE)) {
    stem <- sub("\\.pt$", "", pt_file, ignore.case = TRUE)
    resolved <- .resolve_model_name(stem)
    pt_file <- paste0(resolved, ".pt")
  }

  cands <- c(
    file.path(output_dir, pt_file),
    file.path("D:/Desktop/models", pt_file),
    file.path(pliman_model_dir(), pt_file)
  )
  for (cand in unique(cands)) {
    if (file.exists(cand) && file.info(cand)$size > 100000) {
      return(normalizePath(cand, winslash = "/"))
    }
  }

  dest_file <- file.path(output_dir, pt_file)
  nepem_base <- pliman_models_base_url()
  nepem_url <- paste0(nepem_base, pt_file)

  cli::cli_alert_info("Pre-trained weights {.val {pt_file}} not found locally.")
  cli::cli_alert_info("Downloading from NEPEM repository ({.url {nepem_url}})...")

  dl_success <- .download_url_robust(nepem_url, dest_file, min_size = 100000)

  if (!dl_success) {
    cli::cli_alert_warning("Download from NEPEM failed. Retrying with Ultralytics assets...")
    ultra_tag <- if (grepl("yolo26", pt_file, ignore.case = TRUE)) {
      "v8.4.0"
    } else if (grepl("yolo11", pt_file, ignore.case = TRUE)) {
      "v8.3.0"
    } else {
      "v8.2.0"
    }
    actual_ultra_file <- sub("^yolo26", "yolo11", pt_file, ignore.case = TRUE)
    ultra_url <- paste0("https://github.com/ultralytics/assets/releases/download/", ultra_tag, "/", actual_ultra_file)
    dl_success <- .download_url_robust(ultra_url, dest_file, min_size = 100000)
  }

  if (!dl_success || !file.exists(dest_file)) {
    cli::cli_abort(
      c(
        "x" = "Failed to download model weights {.val {pt_file}}.",
        "i" = "Please check your internet connection or place {.file {pt_file}} in {.path {output_dir}}."
      )
    )
  }

  cli::cli_alert_success("Model weights {.file {pt_file}} ready!")
  normalizePath(dest_file, winslash = "/")
}


# ==============================================================================
# YOLO TRAINING & EXPORT PIPELINE
# ==============================================================================

#' @title Export a Trained PyTorch YOLO Checkpoint (.pt) to ONNX
#' @name yolo_export
#' @description
#' Converts a trained PyTorch YOLO checkpoint (`best.pt` or any custom `.pt` model)
#' into a standalone ONNX model ready for instant C++ execution in [image_segment_dl()]
#' or [image_detect_dl()].
#'
#' @param model Path to the `.pt` weights file.
#' @param output_dir Target directory to save the exported `.onnx` file. Defaults to [pliman_model_dir()].
#' @param output_name Filename of the exported model. If `NULL` (default), derived from `model`.
#' @param imgsz Input image resolution (default `640`).
#' @param auto_install Logical. If `TRUE` (default), automatically installs missing export libraries.
#'
#' @return The full path to the exported `.onnx` model.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#' yolo_export("best.pt", output_name = "custom_model.onnx")
#' }
yolo_export <- function(model,
                        output_dir = pliman_model_dir(),
                        output_name = NULL,
                        imgsz = 640,
                        auto_install = TRUE) {
  if (!file.exists(model)) {
    cli::cli_abort("Model checkpoint file not found at {.path {model}}.")
  }
  model <- normalizePath(model, winslash = "/")
  output_dir <- pliman_model_dir(output_dir)

  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)

  if (is.null(output_name) || !nzchar(output_name)) {
    output_name <- paste0(tools::file_path_sans_ext(basename(model)), ".onnx")
  }
  if (!grepl("\\.onnx$", output_name, ignore.case = TRUE)) {
    output_name <- paste0(output_name, ".onnx")
  }
  dest_onnx <- file.path(output_dir, output_name)

  run_dir <- file.path(tempdir(), "yolo_export_run")
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  py_export_script <- file.path(run_dir, "export.py")
  py_code <- paste0(
    "import sys\n",
    "from ultralytics import YOLO\n",
    "print('--> Loading weights from: ", model, "')\n",
    "m = YOLO('", model, "')\n",
    "try:\n",
    "    out = m.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=False)\n",
    "except Exception as e:\n",
    "    print(f'Standard export notice: {e}. Retrying with opset=17...')\n",
    "    out = m.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=False, opset=17)\n",
    "print('--> Exported successfully:', out)\n"
  )
  writeLines(py_code, py_export_script)

  cli::cli_alert_info("Exporting {.file {basename(model)}} to standalone ONNX format...")
  status <- system2(py_exec, args = shQuote(py_export_script))
  if (status != 0L) {
    cli::cli_abort("ONNX export failed with exit code {status}.")
  }

  cand_onnx <- file.path(dirname(model), paste0(tools::file_path_sans_ext(basename(model)), ".onnx"))
  if (!file.exists(cand_onnx)) {
    found <- list.files(dirname(model), pattern = "\\.onnx$", full.names = TRUE)
    if (length(found) > 0L) cand_onnx <- found[1]
  }

  if (!file.exists(cand_onnx)) {
    cli::cli_abort("Could not locate exported ONNX file near {.path {model}}.")
  }

  file.copy(cand_onnx, dest_onnx, overwrite = TRUE)
  cli::cli_alert_success("Exported ONNX saved to: {.file {dest_onnx}}")
  invisible(dest_onnx)
}

# Helper function to format seconds into human-readable strings (s, m s, h m s)
.format_hms <- function(s) {
  if (is.null(s) || !is.numeric(s) || is.na(s) || s < 0) return("0s")
  if (s < 60) {
    sprintf("%.1fs", s)
  } else if (s < 3600) {
    mins <- floor(s / 60)
    secs <- round(s %% 60, 1)
    sprintf("%dm %04.1fs", mins, secs)
  } else {
    hours <- floor(s / 3600)
    rem <- s %% 3600
    mins <- floor(rem / 60)
    secs <- round(rem %% 60)
    sprintf("%dh %02dm %02ds", hours, mins, secs)
  }
}

# Helper to serialize R value to Python dictionary literal
.r_to_py_val <- function(v) {
  if (is.null(v)) return("None")
  if (is.logical(v)) return(if (isTRUE(v)) "True" else "False")
  if (is.numeric(v)) {
    if (length(v) == 1L) return(as.character(v))
    return(paste0("[", paste(v, collapse = ", "), "]"))
  }
  if (is.character(v)) {
    if (length(v) == 1L) return(paste0("'", gsub("'", "\\\\'", v), "'"))
    return(paste0("[", paste0("'", gsub("'", "\\\\'", v), "'", collapse = ", "), "]"))
  }
  return(paste0("'", as.character(v), "'"))
}

#' @title Train and Fine-Tune YOLO Models for Segmentation or Detection
#' @name yolo_train
#' @description
#' High-level training interface for Ultralytics YOLO models (YOLO26, YOLO11, YOLOv8)
#' in `pliman`. Supports instance segmentation (`task = "segment"`) and object
#' detection (`task = "detect"`), live formatted training progress, early stopping,
#' automatic `.pt` and `.onnx` export, and rich S3 diagnostic visualization (`plot()`).
#'
#' @details
#' \subsection{Training Lifecycle and Pipeline}{
#' \enumerate{
#'   \item \strong{Dataset Structure Verification}: Inspects bounding box and polygon coordinates in `data.yaml`
#'     to confirm task configuration, automatically reconciling task mismatches (e.g. bounding box dataset with segmentation model).
#'   \item \strong{Hardware Acceleration}: Auto-detects NVIDIA CUDA GPUs with cuDNN acceleration, falling back gracefully
#'     to multi-threaded CPU processing.
#'   \item \strong{Pre-trained Weights}: Automatically verifies or downloads official pre-trained YOLO weights from NEPEM repositories.
#'   \item \strong{Training Execution}: Runs Ultralytics training with real-time CLI feedback, early stopping patience,
#'     and automated learning rate scheduling.
#'   \item \strong{Pure Numerical Evaluation Export}: Extracts tabular CSV metrics (\code{confusion_matrix.csv},
#'     \code{confusion_matrix_normalized.csv}, \code{pr_curve.csv}, \code{f1_curve.csv}, and \code{results.csv})
#'     while suppressing low-resolution matplotlib raster images.
#'   \item \strong{Dual Model Export}: Saves both PyTorch \code{.pt} weights (enabling future fine-tuning or transfer learning)
#'     and optimized \code{.onnx} models (for pure C++ inference via ONNX Runtime in \code{pliman} without requiring Python).
#' }
#' }
#'
#' \subsection{Loss Functions and Mathematical Formulations}{
#' YOLO models optimize a composite multi-task objective balancing bounding box regression, classification,
#' boundary distribution focal loss, and (for instance segmentation) pixel-wise mask prediction:
#' \deqn{\mathcal{L}_{\text{total}} = \lambda_{\text{box}} \mathcal{L}_{\text{box}} + \lambda_{\text{cls}} \mathcal{L}_{\text{cls}} + \lambda_{\text{dfl}} \mathcal{L}_{\text{dfl}} + \lambda_{\text{seg}} \mathcal{L}_{\text{seg}}}{L_total = lambda_box * L_box + lambda_cls * L_cls + lambda_dfl * L_dfl + lambda_seg * L_seg}
#'
#' \strong{1. Complete IoU (CIoU) Bounding Box Loss (\eqn{\mathcal{L}_{\text{box}}}):}
#' Bounding box regression uses the Complete IoU loss, which penalizes bounding box overlap error, normalized distance
#' between center coordinates, and aspect ratio discrepancy:
#' \deqn{\mathcal{L}_{\text{CIoU}} = 1 - \text{IoU} + \frac{\rho^2(b, b^{\text{gt}})}{c^2} + \alpha v}{L_CIoU = 1 - IoU + (distance^2 / c^2) + alpha * v}
#' where \eqn{\rho(b, b^{\text{gt}})} is the Euclidean distance between the center points of predicted box \eqn{b} and ground truth box \eqn{b^{\text{gt}}},
#' \eqn{c} is the diagonal length of the smallest enclosing bounding box covering both boxes, and \eqn{v} quantifies aspect ratio consistency:
#' \deqn{v = \frac{4}{\pi^2} \left( \arctan\frac{w^{\text{gt}}}{h^{\text{gt}}} - \arctan\frac{w}{h} \right)^2}{v = (4 / pi^2) * (atan(w_gt / h_gt) - atan(w / h))^2}
#' \deqn{\alpha = \frac{v}{(1 - \text{IoU}) + v}}{alpha = v / ((1 - IoU) + v)}
#'
#' \strong{2. Classification Loss (\eqn{\mathcal{L}_{\text{cls}}}):}
#' Class probabilities are optimized using Binary Cross-Entropy (BCE) with sigmoid activation across all classes:
#' \deqn{\mathcal{L}_{\text{cls}} = -\frac{1}{N_c} \sum_{c=1}^{N_c} \left[ y_c \log(\hat{y}_c) + (1 - y_c) \log(1 - \hat{y}_c) \right]}{L_cls = -sum[ y * log(p) + (1 - y) * log(1 - p) ]}
#' where \eqn{y_c \in \{0, 1\}} represents ground truth class membership and \eqn{\hat{y}_c \in [0, 1]} is the predicted class probability.
#'
#' \strong{3. Instance Segmentation Mask Loss (\eqn{\mathcal{L}_{\text{seg}}}):}
#' For segmentation models, a pixel-wise binary cross-entropy loss compares the predicted prototype mask \eqn{\hat{M}} with the ground truth binary mask \eqn{M} across mask domain \eqn{\Omega}:
#' \deqn{\mathcal{L}_{\text{seg}} = -\frac{1}{|\Omega|} \sum_{(x,y) \in \Omega} \left[ M(x,y) \log \hat{M}(x,y) + (1 - M(x,y)) \log(1 - \hat{M}(x,y)) \right]}{L_seg = -(1 / |Omega|) * sum[ M * log(M_hat) + (1 - M) * log(1 - M_hat) ]}
#'
#' \strong{4. Distribution Focal Loss (\eqn{\mathcal{L}_{\text{dfl}}}):}
#' DFL models continuous bounding box coordinates \eqn{y} as probability distributions over discrete bins \eqn{\{y_i, y_{i+1}\}} around the target edges:
#' \deqn{\mathcal{L}_{\text{dfl}}(S_i, S_{i+1}) = -\left[ (y_{i+1} - y) \log(S_i) + (y - y_i) \log(S_{i+1}) \right]}{L_dfl(S_i, S_{i+1}) = -[ (y_{i+1} - y) * log(S_i) + (y - y_i) * log(S_{i+1}) ]}
#' where \eqn{S_i} and \eqn{S_{i+1}} are softmax probabilities of adjacent integer anchor bins satisfying \eqn{y_i \le y \le y_{i+1}}.
#' This enables flexible sub-pixel boundary localization even under heavy occlusion or blurred object boundaries.
#'
#' \strong{5. Cosine Annealing Learning Rate Schedule:}
#' During the warm-up phase (\eqn{t < T_{\text{warm}}}), the learning rate increases linearly from 0 to \eqn{\eta_0}.
#' For remaining epochs (\eqn{t \ge T_{\text{warm}}}), the learning rate decays according to a cosine annealing schedule:
#' \deqn{\eta_t = \eta_{\text{final}} + \frac{1}{2} (\eta_0 - \eta_{\text{final}}) \left( 1 + \cos\left( \frac{t - T_{\text{warm}}}{T_{\text{max}} - T_{\text{warm}}} \pi \right) \right)}{eta_t = eta_final + 0.5 * (eta_0 - eta_final) * (1 + cos((t - T_warm) / (T_max - T_warm) * pi))}
#' where \eqn{\eta_0 = \text{lr0}} and \eqn{\eta_{\text{final}} = \text{lr0} \times \text{lrf}}.
#' }
#'
#' @param data Path to the `data.yaml` file (e.g. `"dataset_cafe/data.yaml"`).
#' @param model Base pre-trained model architecture to fine-tune. Defaults to `"yolo26n-seg"`.
#'   Can be `"yolo26n-seg"`, `"yolo26s-seg"`, `"yolo11n-seg"`, `"yolo26n"`, or path to a `.pt` file.
#' @param epochs Number of training epochs (default `50`).
#' @param imgsz Input image resolution for training (default `640`).
#' @param batch Batch size (default `16`).
#' @param device Device to train on: `"auto"` (default), `"0"` (GPU 0), `"cpu"`, etc.
#' @param output_dir Target directory to save the final exported `.onnx` model.
#'   Defaults to [pliman_model_dir()].
#' @param output_name Filename of the exported model. If `NULL` (default),
#'   an informative name is automatically created (e.g. `"yolo_custom_seg.onnx"`).
#' @param results_dir Directory where training plots, confusion matrices, and `results.csv`
#'   are permanently preserved. If `NULL` (default), creates a folder named `<output_name>_results`
#'   inside `output_dir`.
#' @param export_pt Logical. If `TRUE` (default), also exports the PyTorch `.pt`
#'   weights (e.g., `model.pt`) alongside the `.onnx` model, allowing future
#'   re-training or fine-tuning on more complex datasets.
#' @param patience Early stopping patience in epochs. Defaults to `15`.
#' @param lr0 Initial learning rate (default `0.01`).
#' @param lrf Final learning rate fraction (default `0.01`).
#' @param optimizer Optimizer to use (`"auto"`, `"SGD"`, `"Adam"`, `"AdamW"`, `"RMSProp"`). Default is `"auto"`.
#' @param close_mosaic Disable mosaic augmentation in the last N epochs. Default is `10`.
#' @param amp Automatic Mixed Precision (AMP). Set to `FALSE` to disable AMP if GPU checks fail. Default is `NULL` (auto-detected).
#' @param plots Logical. If `FALSE` (default), suppresses generation of low-resolution matplotlib raster image files from YOLO; native vector plots are generated directly in R from exported CSV data.
#' @param workers Number of worker threads for data loading. Defaults to `0` on Windows and `4` on Linux/macOS.
#' @param auto_install Logical. If `TRUE`, automatically installs missing dependencies via `pip`. Defaults to `TRUE`.
#' @param ... Additional arguments passed directly to `model.train()` in Ultralytics
#'   (e.g. `weight_decay = 0.0005`, `momentum = 0.937`, `mosaic = 0.5`, `mixup = 0.1`, `degrees = 10`, `box = 7.5`, `cls = 0.5`).
#'
#' @return An object of class `yolo_train` containing:
#'   * `model_file`: Path to the exported `.onnx` model ready for `pliman`.
#'   * `pt_file`: Path to the saved `.pt` PyTorch weights (if `export_pt = TRUE`).
#'   * `results_dir`: Directory containing training plots, metrics, and logs.
#'   * `plots`: Named list of file paths to generated diagnostic plots.
#'   * `metrics`: Data frame of per-epoch loss and mAP progressions (`results.csv`).
#'   * `confusion_matrix`: Data frame of the final confusion matrix (counts).
#'   * `confusion_matrix_norm`: Data frame of the normalized confusion matrix.
#'   * `pr_curve`: Data frame of the Precision-Recall curve.
#'   * `f1_curve`: Data frame of the F1 vs Confidence curve.
#'   * `best_metrics`: Data frame row of the best performing epoch.
#'   * `epochs`: Number of epochs trained.
#'   * `elapsed_sec`: Elapsed training time in seconds.
#'   * `elapsed_formatted`: Formatted training time (e.g. `"3h 35m 14s"`).
#'   * `task`: Task type (`"segment"` or `"detect"`).
#'   * `call`: Matched function call.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Train a custom coffee grain segmentation model
#' model_res <- yolo_train(
#'   data = "dataset_cafe/data.yaml",
#'   model = "yolo26n-seg",
#'   epochs = 30,
#'   output_name = "cafe_seg.onnx"
#' )
#'
#' # Inspect training summary and plot diagnostic dashboard
#' print(model_res)
#' plot(model_res)
#'
#' # Plot specific evaluation curves:
#' plot(model_res, which = "confusion_matrix")
#' plot(model_res, which = "pr_curve")
#' plot(model_res, which = "metrics")
#'
#' # Run inference with your new model directly in pliman!
#' img <- image_import("cafe.jpeg")
#' seg <- image_segment_dl(img, model = model_res$model_file, type = "highlight")
#' }
yolo_train <- function(data = "yolo_dataset/data.yaml",
                       model = "yolo26n-seg",
                       epochs = 50,
                       imgsz = 640,
                       batch = 16,
                       device = "auto",
                       output_dir = pliman_model_dir(),
                       output_name = NULL,
                       results_dir = NULL,
                       export_pt = TRUE,
                       patience = 15,
                       lr0 = 0.01,
                       lrf = 0.01,
                       optimizer = "auto",
                       close_mosaic = 10,
                       amp = NULL,
                       plots = FALSE,
                       workers = if (.Platform$OS.type == "windows") 0L else 4L,
                       auto_install = TRUE,
                       ...) {
  if (is.list(data)) {
    if (!is.null(data$data_yaml) && is.character(data$data_yaml) && nzchar(data$data_yaml[1])) {
      data <- data$data_yaml[1]
    } else if (!is.null(data$yaml_file) && is.character(data$yaml_file) && nzchar(data$yaml_file[1])) {
      data <- data$yaml_file[1]
    } else if (!is.null(data$dir) && is.character(data$dir) && nzchar(data$dir[1])) {
      data <- data$dir[1]
    }
  }

  data <- normalizePath(data, winslash = "/", mustWork = FALSE)
  if (!file.exists(data) && !dir.exists(data)) {
    cli::cli_abort("Could not find dataset configuration file or directory at {.path {data}}.")
  }

  output_dir <- pliman_model_dir(output_dir)

  # 1. Ensure Python runtime & hardware environment
  py_exec <- .ensure_python(auto_install = auto_install)
  env_info <- .setup_python_yolo_env(py_exec, auto_install = auto_install)

  # 2. Inspect dataset structure to prevent task mismatch
  lbl_dir <- if (dir.exists(file.path(data, "labels"))) {
    file.path(data, "labels")
  } else {
    file.path(dirname(data), "labels")
  }
  lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  dataset_task <- "unknown"
  if (length(lbl_files) > 0L) {
    for (lf in lbl_files[seq_len(min(10L, length(lbl_files)))]) {
      lines <- readLines(lf, warn = FALSE)
      lines <- lines[nzchar(trimws(lines))]
      for (ln in lines) {
        vals <- strsplit(trimws(ln), "\\s+")[[1]]
        if (length(vals) == 5L) {
          dataset_task <- "detect"
          break
        } else if (length(vals) >= 6L) {
          dataset_task <- "segment"
          break
        }
      }
      if (dataset_task != "unknown") break
    }
  } else {
    # Check if dataset is a classification directory (contains train/ and val/ class folders)
    train_dir <- if (dir.exists(file.path(data, "train"))) file.path(data, "train") else file.path(dirname(data), "train")
    val_dir   <- if (dir.exists(file.path(data, "val"))) file.path(data, "val") else file.path(dirname(data), "val")
    if (dir.exists(train_dir) && dir.exists(val_dir)) {
      tr_sub <- list.dirs(train_dir, recursive = FALSE, full.names = FALSE)
      if (length(tr_sub) >= 2L) {
        dataset_task <- "classify"
      }
    }
  }

  if (identical(dataset_task, "classify") && !grepl("cls", model, ignore.case = TRUE)) {
    cls_model <- if (grepl("yolo26", model, ignore.case = TRUE)) "yolo26n-cls" else "yolo11n-cls"
    cli::cli_alert_warning(
      "The dataset at {.path {data}} is structured for image classification ({.val classify}), but model {.val {model}} was requested."
    )
    cli::cli_alert_info(
      "Automatically adjusting base model to classification ({.val {cls_model}})."
    )
    model <- cls_model
    if (missing(imgsz) || imgsz == 640) imgsz <- 224
  } else if (identical(dataset_task, "detect") && grepl("seg", model, ignore.case = TRUE)) {
    det_model <- sub("[-_]seg", "", model, ignore.case = TRUE)
    cli::cli_alert_warning(
      "The dataset at {.path {data}} contains bounding boxes ({.val detect}), but a segmentation model ({.val {model}}) was requested."
    )
    cli::cli_alert_info(
      "Automatically adjusting base model to detection ({.val {det_model}}) to avoid training failure."
    )
    cli::cli_alert_info(
      "Tip: To train instance segmentation (YOLO-seg), re-export your dataset with {.code task = 'segment'} in {.fn yolo_dataset_export}."
    )
    model <- det_model
    if (!is.null(output_name) && grepl("_seg", output_name, ignore.case = TRUE)) {
      output_name <- sub("_seg", "_det", output_name, ignore.case = TRUE)
    }
  }

  # 3. Resolve or download base pre-trained model (.pt)
  base_model <- .resolve_yolo_pt(model, output_dir = output_dir)

  # 4. Detect task (segment, obb, pose, classify, or detect)
  is_seg  <- grepl("seg", base_model, ignore.case = TRUE)
  is_obb  <- grepl("obb", base_model, ignore.case = TRUE)
  is_pose <- grepl("pose", base_model, ignore.case = TRUE)
  is_cls  <- grepl("cls", base_model, ignore.case = TRUE)
  task_mode <- if (is_seg) {
    "segment"
  } else if (is_obb) {
    "obb"
  } else if (is_pose) {
    "pose"
  } else if (is_cls) {
    "classify"
  } else {
    "detect"
  }

  data_folder_name <- if (dir.exists(data)) basename(data) else basename(dirname(data))
  if (is.null(output_name) || !nzchar(output_name)) {
    prefix <- if (nzchar(data_folder_name) && data_folder_name != ".") data_folder_name else "custom_yolo"
    suffix <- if (is_seg) "_seg" else if (is_obb) "_obb" else if (is_pose) "_pose" else if (is_cls) "_cls" else "_det"
    output_name <- paste0(prefix, suffix)
  }
  base_stem <- sub("\\.(onnx|pt)$", "", output_name, ignore.case = TRUE)
  dest_onnx <- file.path(output_dir, paste0(base_stem, ".onnx"))
  dest_pt <- file.path(output_dir, paste0(base_stem, ".pt"))

  if (is.null(results_dir) || !nzchar(results_dir)) {
    results_dir <- file.path(output_dir, paste0(base_stem, "_results"))
  }
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 4. Display training configuration banner
  hw_desc <- if (isTRUE(env_info$cuda_available)) {
    paste0("NVIDIA GPU (", env_info$gpu_info$name, ") \U0001f680")
  } else {
    "Multi-threaded CPU \U0001f4bb"
  }

  cli::cli_h2("YOLO Model Training (pliman)")
  cli::cli_ul(c(
    "Dataset: {.path {data}}",
    "Base model: {.val {basename(base_model)}}",
    "Hardware: {hw_desc}",
    "Epochs: {epochs} (patience: {patience})",
    "Batch size: {batch}",
    "Resolution: {imgsz}x{imgsz} px",
    "Results folder: {.path {results_dir}}"
  ))

  # 5. Build Python training & export script
  run_dir <- file.path(tempdir(), "yolo_train_run")
  if (!dir.exists(run_dir)) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  }

  py_script_path <- file.path(run_dir, "train_and_export.py")

  # Hyperparameters and training options
  train_args <- list(
    data = data,
    epochs = as.integer(epochs),
    imgsz = as.integer(imgsz),
    batch = as.integer(batch),
    patience = as.integer(patience),
    workers = as.integer(workers),
    lr0 = as.numeric(lr0),
    lrf = as.numeric(lrf),
    optimizer = optimizer,
    close_mosaic = as.integer(close_mosaic),
    plots = isTRUE(plots)
  )
  if (!is.null(amp)) {
    train_args$amp <- isTRUE(amp)
  }

  dots <- list(...)
  for (nm in names(dots)) {
    if (nzchar(nm)) {
      train_args[[nm]] <- dots[[nm]]
    }
  }

  kwargs_py_lines <- vapply(names(train_args), function(k) {
    paste0("        '", k, "': ", .r_to_py_val(train_args[[k]]))
  }, character(1))
  kwargs_py_code <- paste(kwargs_py_lines, collapse = ",\n")

  base_model_norm <- normalizePath(base_model, winslash = "/", mustWork = FALSE)
  run_dir_norm <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  is_seg_py <- if (is_seg) "True" else "False"
  is_cls_py <- if (is_cls) "True" else "False"
  dev_str <- as.character(device)

  py_code <- paste0(
    "import sys, logging, time, torch\n",
    "import ultralytics.utils\n",
    "import ultralytics.engine.trainer\n",
    "import ultralytics.engine.validator\n\n",
    "class SilentTQDM:\n",
    "    def __init__(self, iterable=None, *args, **kwargs):\n",
    "        self.iterable = iterable if iterable is not None else []\n",
    "    def __iter__(self):\n",
    "        return iter(self.iterable)\n",
    "    def __enter__(self):\n",
    "        return self\n",
    "    def __exit__(self, *args):\n",
    "        pass\n",
    "    def update(self, *args, **kwargs):\n",
    "        pass\n",
    "    def set_description(self, *args, **kwargs):\n",
    "        pass\n",
    "    def set_postfix(self, *args, **kwargs):\n",
    "        pass\n",
    "    def close(self):\n",
    "        pass\n",
    "    def refresh(self):\n",
    "        pass\n",
    "    def clear(self):\n",
    "        pass\n",
    "    @staticmethod\n",
    "    def write(*args, **kwargs):\n",
    "        pass\n\n",
    "ultralytics.utils.TQDM = SilentTQDM\n",
    "ultralytics.engine.trainer.TQDM = SilentTQDM\n",
    "ultralytics.engine.validator.TQDM = SilentTQDM\n",
    "ultralytics.utils.LOGGER.setLevel(logging.WARNING)\n\n",
    "from ultralytics import YOLO\n\n",
    "best_fitness = -1.0\n",
    "best_epoch = 1\n",
    "best_map50 = 0.0\n",
    "t_start = None\n",
    "ep_start = None\n",
    "is_seg = ", is_seg_py, "\n",
    "is_cls = ", is_cls_py, "\n\n",
    "def format_time(s):\n",
    "    if s is None or s < 0:\n",
    "        return '0s'\n",
    "    s = round(float(s), 1)\n",
    "    if s < 60:\n",
    "        return f'{s:.1f}s'\n",
    "    elif s < 3600:\n",
    "        m = int(s // 60)\n",
    "        rem_s = round(s - (m * 60), 1)\n",
    "        if rem_s >= 59.95:\n",
    "            m += 1\n",
    "            rem_s = 0.0\n",
    "        return f'{m}m {rem_s:04.1f}s'\n",
    "    else:\n",
    "        tot = int(round(s))\n",
    "        h = int(tot // 3600)\n",
    "        rem = tot - (h * 3600)\n",
    "        m = int(rem // 60)\n",
    "        rem_s = int(rem - (m * 60))\n",
    "        return f'{h}h {m:02d}m {rem_s:02d}s'\n\n",
    "def on_train_start(trainer):\n",
    "    global t_start\n",
    "    t_start = time.time()\n",
    "    print('')\n",
    "    if is_cls:\n",
    "        print('  Epoch    GPU Mem       Loss     Top1-Acc     Top5-Acc          Time')\n",
    "        print('  ' + '-' * 70)\n",
    "    elif is_seg:\n",
    "        print('  Epoch    GPU Mem     Box Loss     Seg Loss     Cls Loss       mAP50     mAP50-95          Time')\n",
    "        print('  ' + '-' * 88)\n",
    "    else:\n",
    "        print('  Epoch    GPU Mem     Box Loss     Cls Loss       mAP50     mAP50-95          Time')\n",
    "        print('  ' + '-' * 76)\n",
    "    sys.stdout.flush()\n\n",
    "def on_train_epoch_start(trainer):\n",
    "    global ep_start\n",
    "    ep_start = time.time()\n\n",
    "def on_fit_epoch_end(trainer):\n",
    "    global best_fitness, best_epoch, best_map50\n",
    "    ep = trainer.epoch + 1\n",
    "    total = trainer.epochs\n",
    "    if ep > total:\n",
    "        return\n\n",
    "    if torch.cuda.is_available():\n",
    "        gpu_mem = f'{torch.cuda.memory_reserved() / (1024**3):.1f} GB'\n",
    "    else:\n",
    "        gpu_mem = 'CPU'\n\n",
    "    tloss = getattr(trainer, 'tloss', {})\n",
    "    metrics = getattr(trainer, 'metrics', {})\n",
    "    if is_cls:\n",
    "        c_loss = float(tloss[0]) if isinstance(tloss, (list, tuple)) and len(tloss) > 0 else float(getattr(trainer, 'loss', 0.0))\n",
    "        top1 = float(metrics.get('metrics/accuracy_top1', 0.0))\n",
    "        top5 = float(metrics.get('metrics/accuracy_top5', 0.0))\n",
    "        fitness = top1\n",
    "        m50 = top1\n",
    "    else:\n",
    "        if isinstance(tloss, dict):\n",
    "            b_loss = float(tloss.get('box_loss', 0.0))\n",
    "            s_loss = float(tloss.get('seg_loss', 0.0))\n",
    "            c_loss = float(tloss.get('cls_loss', 0.0))\n",
    "        elif isinstance(tloss, (list, tuple)):\n",
    "            b_loss = float(tloss[0]) if len(tloss) > 0 else 0.0\n",
    "            s_loss = float(tloss[1]) if len(tloss) > 2 else 0.0\n",
    "            c_loss = float(tloss[-1]) if len(tloss) > 1 else 0.0\n",
    "        else:\n",
    "            b_loss, s_loss, c_loss = 0.0, 0.0, 0.0\n",
    "        m50 = float(metrics.get('metrics/mAP50(B)', metrics.get('metrics/mAP50(M)', 0.0)))\n",
    "        m95 = float(metrics.get('metrics/mAP50-95(B)', metrics.get('metrics/mAP50-95(M)', 0.0)))\n",
    "        fitness = float(getattr(trainer, 'fitness', 0.0))\n\n",
    "    is_best = False\n",
    "    if fitness > best_fitness:\n",
    "        best_fitness = fitness\n",
    "        best_epoch = ep\n",
    "        best_map50 = m50\n",
    "        is_best = True\n\n",
    "    ep_t = format_time(time.time() - ep_start) if ep_start else ''\n",
    "    tag = '  * (best)' if is_best else ''\n\n",
    "    if is_cls:\n",
    "        line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}   {c_loss:>8.4f}     {top1:>7.3f}      {top5:>7.3f}    {ep_t:>10}{tag}'\n",
    "    elif is_seg:\n",
    "        line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}     {b_loss:>8.4f}     {s_loss:>8.4f}     {c_loss:>8.4f}     {m50:>7.3f}     {m95:>7.3f}    {ep_t:>10}{tag}'\n",
    "    else:\n",
    "        line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}     {b_loss:>8.4f}     {c_loss:>8.4f}     {m50:>7.3f}     {m95:>7.3f}    {ep_t:>10}{tag}'\n\n",
    "    print(line)\n",
    "    sys.stdout.flush()\n\n",
    "def on_train_end(trainer):\n",
    "    if is_cls:\n",
    "        print('  ' + '-' * 70)\n",
    "        elapsed_str = format_time(time.time() - t_start)\n",
    "        print(f'  [OK] Training finished in {elapsed_str} (Best Accuracy Top-1: {best_fitness:.3f} at epoch {best_epoch})')\n",
    "    else:\n",
    "        print('  ' + '-' * (88 if is_seg else 76))\n",
    "        elapsed_str = format_time(time.time() - t_start)\n",
    "        print(f'  [OK] Training finished in {elapsed_str} (Best mAP50: {best_map50:.3f} at epoch {best_epoch})')\n",
    "    print('')\n",
    "    sys.stdout.flush()\n\n",
    "def main():\n",
    "    dev = '", dev_str, "'\n",
    "    if dev == 'auto':\n",
    "        dev = 0 if torch.cuda.is_available() else 'cpu'\n",
    "    elif dev.isdigit():\n",
    "        dev = int(dev)\n\n",
    "    model = YOLO(r'", base_model_norm, "')\n",
    "    model.add_callback('on_train_start', on_train_start)\n",
    "    model.add_callback('on_train_epoch_start', on_train_epoch_start)\n",
    "    model.add_callback('on_fit_epoch_end', on_fit_epoch_end)\n",
    "    model.add_callback('on_train_end', on_train_end)\n\n",
    "    train_kwargs = {\n",
    kwargs_py_code, "\n",
    "    }\n",
    "    train_kwargs['device'] = dev\n",
    "    train_kwargs['project'] = r'", run_dir_norm, "'\n",
    "    train_kwargs['name'] = 'train_run'\n",
    "    train_kwargs['exist_ok'] = True\n\n",
    "    results = model.train(**train_kwargs)\n\n",
    "    best_pt = r'", run_dir_norm, "/train_run/weights/best.pt'\n",
    "    best_model = YOLO(best_pt)\n",
    "    try:\n",
    "        exported_path = best_model.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=False)\n",
    "    except Exception as e:\n",
    "        print(f'Standard export notice: {e}. Retrying with opset=17...')\n",
    "        exported_path = best_model.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=False, opset=17)\n\n",
    "    import csv, os, numpy as np\n",
    "    from pathlib import Path\n",
    "    target_dir = r'", run_dir_norm, "/train_run'\n",
    "    try:\n",
    "        val_res = best_model.val(plots=True, workers=0)\n",
    "        cm = getattr(val_res, 'confusion_matrix', None)\n",
    "        if cm is not None and hasattr(cm, 'matrix') and cm.matrix is not None:\n",
    "            matrix = np.array(cm.matrix)\n",
    "            n_rows, n_cols = matrix.shape\n",
    "            raw_names = [best_model.names[i] for i in range(len(best_model.names))]\n",
    "            col_names = (raw_names + ['background']) if n_cols > len(raw_names) else raw_names[:n_cols]\n",
    "            row_names = raw_names[:n_rows]\n",
    "            cm_csv = os.path.join(target_dir, 'confusion_matrix.csv')\n",
    "            with open(cm_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                writer = csv.writer(f)\n",
    "                writer.writerow(['true_class'] + col_names)\n",
    "                for i, row in enumerate(matrix):\n",
    "                    r_lbl = row_names[i] if i < len(row_names) else f'class_{i}'\n",
    "                    writer.writerow([r_lbl] + [int(val) for val in row])\n",
    "            row_sums = matrix.sum(axis=1, keepdims=True)\n",
    "            norm_matrix = np.divide(matrix, row_sums, out=np.zeros_like(matrix, dtype=float), where=row_sums > 0)\n",
    "            norm_csv = os.path.join(target_dir, 'confusion_matrix_normalized.csv')\n",
    "            with open(norm_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                writer = csv.writer(f)\n",
    "                writer.writerow(['true_class'] + col_names)\n",
    "                for i, row in enumerate(norm_matrix):\n",
    "                    r_lbl = row_names[i] if i < len(row_names) else f'class_{i}'\n",
    "                    writer.writerow([r_lbl] + [round(float(val), 4) for val in row])\n",
    "        box = getattr(val_res, 'box', None)\n",
    "        if box is not None and hasattr(box, 'curves_results') and len(box.curves_results) >= 1:\n",
    "            cr_pr = box.curves_results[0]\n",
    "            x_rec = cr_pr[0]\n",
    "            y_prec = cr_pr[1]\n",
    "            cls_cols = [best_model.names[i] for i in range(y_prec.shape[0])]\n",
    "            pr_csv = os.path.join(target_dir, 'pr_curve.csv')\n",
    "            with open(pr_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                writer = csv.writer(f)\n",
    "                writer.writerow(['recall'] + cls_cols)\n",
    "                for j in range(len(x_rec)):\n",
    "                    row_vals = [round(float(x_rec[j]), 4)] + [round(float(y_prec[c, j]), 4) for c in range(y_prec.shape[0])]\n",
    "                    writer.writerow(row_vals)\n",
    "            if len(box.curves_results) >= 2:\n",
    "                cr_f1 = box.curves_results[1]\n",
    "                x_conf = cr_f1[0]\n",
    "                y_f1 = cr_f1[1]\n",
    "                f1_csv = os.path.join(target_dir, 'f1_curve.csv')\n",
    "                with open(f1_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                    writer = csv.writer(f)\n",
    "                    writer.writerow(['confidence'] + cls_cols)\n",
    "                    for j in range(len(x_conf)):\n",
    "                        row_vals = [round(float(x_conf[j]), 4)] + [round(float(y_f1[c, j]), 4) for c in range(y_f1.shape[0])]\n",
    "                        writer.writerow(row_vals)\n",
    "    except Exception as e:\n",
    "        print(f'Notice during metrics export: {e}')\n\n",
    "    for ext in ('*.png', '*.jpg', '*.jpeg'):\n",
    "        for p in Path(r'", run_dir_norm, "').rglob(ext):\n",
    "            try:\n",
    "                p.unlink()\n",
    "            except:\n",
    "                pass\n\n",
    "if __name__ == '__main__':\n",
    "    main()\n"
  )

  writeLines(py_code, py_script_path)

  cli::cli_alert_info("Training progress:")

  start_time <- Sys.time()
  train_status <- system2(py_exec, args = shQuote(py_script_path))
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)

  if (train_status != 0L) {
    cli::cli_abort("Training process failed with exit code {train_status}.")
  }

  # 6. Preserve all training run artifacts and diagnostic plots in results_dir
  train_run_dir <- file.path(run_dir, "train_run")
  if (dir.exists(train_run_dir)) {
    contents <- list.files(train_run_dir, full.names = TRUE)
    for (item in contents) {
      file.copy(item, results_dir, recursive = TRUE, overwrite = TRUE)
    }
  }

  # Ensure no YOLO matplotlib raster image files remain if plots is FALSE
  if (!isTRUE(plots)) {
    unlink(list.files(results_dir, pattern = "\\.(png|jpg|jpeg)$", recursive = TRUE, full.names = TRUE))
  }

  # 7. Locate exported ONNX model
  cand_onnx <- file.path(results_dir, "weights", "best.onnx")
  if (!file.exists(cand_onnx)) {
    cand_onnx <- file.path(run_dir, "train_run", "weights", "best.onnx")
  }
  if (!file.exists(cand_onnx)) {
    found <- list.files(run_dir, pattern = "\\.onnx$", recursive = TRUE, full.names = TRUE)
    if (length(found) > 0L) {
      cand_onnx <- found[1]
    }
  }
  if (!file.exists(cand_onnx)) {
    found_res <- list.files(results_dir, pattern = "\\.onnx$", recursive = TRUE, full.names = TRUE)
    if (length(found_res) > 0L) {
      cand_onnx <- found_res[1]
    }
  }

  if (!file.exists(cand_onnx)) {
    cli::cli_abort("Could not find the exported ONNX model in {.path {run_dir}}.")
  }

  file.copy(cand_onnx, dest_onnx, overwrite = TRUE)

  # 8. Locate saved PyTorch .pt weights
  cand_pt <- file.path(results_dir, "weights", "best.pt")
  if (!file.exists(cand_pt)) {
    cand_pt <- file.path(run_dir, "train_run", "weights", "best.pt")
  }
  if (!file.exists(cand_pt)) {
    found_pt <- list.files(run_dir, pattern = "best\\.pt$", recursive = TRUE, full.names = TRUE)
    if (length(found_pt) > 0L) {
      cand_pt <- found_pt[1]
    }
  }

  saved_pt <- FALSE
  if (isTRUE(export_pt) && file.exists(cand_pt)) {
    file.copy(cand_pt, dest_pt, overwrite = TRUE)
    saved_pt <- TRUE
  }

  # 9. Catalog diagnostic plots and parse metrics
  find_plot <- function(patterns) {
    for (pat in patterns) {
      matches <- list.files(results_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
      if (length(matches) > 0L) return(matches[1])
    }
    NULL
  }

  plots_list <- list(
    results = find_plot(c("^results\\.png$", "^results\\.jpg$")),
    confusion_matrix_norm = find_plot(c("^confusion_matrix_normalized\\.png$", "^confusion_matrix_normalized\\.jpg$")),
    confusion_matrix = find_plot(c("^confusion_matrix\\.png$", "^confusion_matrix\\.jpg$")),
    pr_curve = find_plot(c("^(Mask|Box)?PR_curve\\.png$", "^(Mask|Box)?PR_curve\\.jpg$")),
    f1_curve = find_plot(c("^(Mask|Box)?F1_curve\\.png$", "^(Mask|Box)?F1_curve\\.jpg$")),
    p_curve = find_plot(c("^(Mask|Box)?P_curve\\.png$", "^(Mask|Box)?P_curve\\.jpg$")),
    r_curve = find_plot(c("^(Mask|Box)?R_curve\\.png$", "^(Mask|Box)?R_curve\\.jpg$")),
    labels = find_plot(c("^labels\\.(jpg|png)$")),
    labels_correlogram = find_plot(c("^labels_correlogram\\.(jpg|png)$"))
  )

  read_csv_safe <- function(filename) {
    p <- file.path(results_dir, filename)
    if (file.exists(p)) {
      tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
    } else {
      NULL
    }
  }

  cm_df <- read_csv_safe("confusion_matrix.csv")
  cm_norm_df <- read_csv_safe("confusion_matrix_normalized.csv")
  pr_df <- read_csv_safe("pr_curve.csv")
  f1_df <- read_csv_safe("f1_curve.csv")

  results_csv <- file.path(results_dir, "results.csv")
  metrics_df <- NULL
  best_metrics <- NULL
  if (file.exists(results_csv)) {
    tryCatch({
      metrics_df <- utils::read.csv(results_csv, strip.white = TRUE, check.names = FALSE)
      names(metrics_df) <- trimws(names(metrics_df))
      map50_col <- grep("mAP50", names(metrics_df), value = TRUE)
      map50_col <- map50_col[!grepl("95", map50_col)]
      if (length(map50_col) > 0L) {
        best_idx <- which.max(metrics_df[[map50_col[1]]])
        if (length(best_idx) > 0L) {
          best_metrics <- metrics_df[best_idx, , drop = FALSE]
        }
      }
    }, error = function(e) NULL)
  }

  elapsed_formatted <- .format_hms(elapsed)
  cli::cli_alert_success("Model trained ({elapsed_formatted}) and saved successfully! \U0001f389")
  cli::cli_ul(c(
    paste0("ONNX: {.file ", dest_onnx, "} (fast, zero-Python C++ inference in pliman)"),
    if (saved_pt) paste0("PyTorch: {.file ", dest_pt, "} (weights for future re-training / fine-tuning)") else character(0),
    paste0("Diagnostics & evaluations: {.path ", results_dir, "}")
  ))

  cli::cli_alert_info("Inference in pliman is 100% C++ (pure ONNX Runtime, zero-Python required):")
  if (is_seg) {
    cli::cli_code(paste0(
      'img <- image_import("your_image.jpg")\n',
      'res <- image_segment_dl(img, model = "', basename(dest_onnx), '")\n',
      'plot(res)'
    ))
  } else {
    cli::cli_code(paste0(
      'img <- image_import("your_image.jpg")\n',
      'res <- image_detect_dl(img, model = "', basename(dest_onnx), '")\n',
      'plot(res)'
    ))
  }

  cli::cli_alert_info("To inspect diagnostic plots (confusion matrix, PR curve, loss curves):")
  cli::cli_code(paste0(
    'plot(model_res)\n',
    'plot(model_res, which = "confusion_matrix")\n',
    'plot(model_res, which = "pr_curve")'
  ))

  res <- structure(
    list(
      model_file = dest_onnx,
      pt_file = if (saved_pt) dest_pt else NULL,
      results_dir = results_dir,
      plots = plots_list,
      metrics = metrics_df,
      best_metrics = best_metrics,
      confusion_matrix = cm_df,
      confusion_matrix_norm = cm_norm_df,
      pr_curve = pr_df,
      f1_curve = f1_df,
      epochs = epochs,
      elapsed_sec = elapsed,
      elapsed_formatted = elapsed_formatted,
      task = task_mode,
      base_model = basename(base_model),
      call = match.call()
    ),
    class = c("yolo_train", "list")
  )

  invisible(res)
}

#' @title Load Training Results and Diagnostics for a YOLO Model
#' @name yolo_results
#' @description
#' Reconstructs a `yolo_train` S3 object from an existing model name, `.onnx`/`.pt` file,
#' or results directory without needing to retrain. The resulting object can be directly
#' printed with [print.yolo_train()] and visualized with [plot.yolo_train()] or [yolo_plot()].
#'
#' @details
#' When a YOLO model is trained via [yolo_train()], all evaluation tables, weights, and logs
#' are preserved in `<output_name>_results`. [yolo_results()] reconstructs the full training
#' session state by parsing:
#' \itemize{
#'   \item \code{results.csv}: Epoch-by-epoch loss progressions, precision, recall, mAP50, and mAP50-95.
#'   \item \code{confusion_matrix.csv}: Complete confusion matrix counts including the background class.
#'   \item \code{confusion_matrix_normalized.csv}: Normalized classification rates.
#'   \item \code{pr_curve.csv}: 1,000-point Precision-Recall evaluation curve per class.
#'   \item \code{f1_curve.csv}: 1,000-point F1 vs. Confidence evaluation curve per class.
#'   \item Model weights: Automatically locates exported \code{.onnx} and \code{.pt} files.
#' }
#' This enables instant, offline performance inspection, publication-ready vector plotting,
#' and threshold selection for any trained model without needing Python or GPU re-execution.
#'
#' @param model Name or path to a trained YOLO model (e.g. `"capsulas_det"`, `"capsulas_det.onnx"`).
#' @param results_dir Optional direct path to the `<model>_results` folder. If `NULL` (default),
#'   the folder is searched in `output_dir`, in `pliman_model_dir()`, or alongside the model file.
#' @param output_dir Directory where pliman models and results are stored. Defaults to [pliman_model_dir()].
#'
#' @return An object of class `c("yolo_train", "list")`.
#' @export
#' @seealso [plot.yolo_train()], [yolo_plot()], [yolo_train()]
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Load diagnostics from a previously trained model
#' res <- yolo_results("capsulas_det")
#'
#' # Inspect training metrics and plot diagnostic dashboard
#' print(res)
#' plot(res)
#' plot(res, which = "confusion_matrix")
#' plot(res, which = "pr_curve")
#' }
yolo_results <- function(model = NULL,
                         results_dir = NULL,
                         output_dir = pliman_model_dir()) {
  output_dir <- pliman_model_dir(output_dir)

  if (is.null(model) && is.null(results_dir)) {
    cli::cli_abort("Please provide either {.arg model} or {.arg results_dir}.")
  }

  base_stem <- if (!is.null(model)) {
    sub("\\.(onnx|pt)$", "", basename(model), ignore.case = TRUE)
  } else {
    sub("_results$", "", basename(results_dir), ignore.case = TRUE)
  }

  # Search for results_dir if not explicitly given
  if (is.null(results_dir) || !dir.exists(results_dir)) {
    stem_clean <- sub("^yolo[-_]", "", base_stem, ignore.case = TRUE)
    cand_dirs <- c(
      if (!is.null(model)) file.path(dirname(model), paste0(base_stem, "_results")),
      if (!is.null(model)) file.path(dirname(model), paste0(stem_clean, "_results")),
      file.path(output_dir, paste0(base_stem, "_results")),
      file.path(output_dir, paste0(stem_clean, "_results")),
      file.path(pliman_model_dir(), paste0(base_stem, "_results")),
      file.path(pliman_model_dir(), paste0(stem_clean, "_results")),
      file.path("D:/Desktop/models", paste0(base_stem, "_results")),
      list.files(output_dir, pattern = paste0("^", stem_clean, ".*_results$"), full.names = TRUE),
      list.files(pliman_model_dir(), pattern = paste0("^", stem_clean, ".*_results$"), full.names = TRUE)
    )
    found_dir <- unique(cand_dirs[dir.exists(cand_dirs)])
    if (length(found_dir) > 0L) {
      results_dir <- found_dir[1]
      base_stem <- sub("_results$", "", basename(results_dir), ignore.case = TRUE)
    } else {
      avail <- list.files(output_dir, pattern = "_results$")
      avail_names <- sub("_results$", "", avail)
      hint_msg <- if (length(avail_names) > 0) {
        paste0(" Available trained models: ", paste0("'", avail_names, "'", collapse = ", "), ".")
      } else {
        ""
      }
      cli::cli_abort("Could not find results directory for model {.val {base_stem}}.{hint_msg}")
    }
  }

  # Find onnx and pt models
  dest_onnx <- file.path(output_dir, paste0(base_stem, ".onnx"))
  if (!file.exists(dest_onnx)) {
    cand_onnx <- c(
      file.path(results_dir, "weights", "best.onnx"),
      file.path(dirname(results_dir), paste0(base_stem, ".onnx")),
      list.files(results_dir, pattern = "\\.onnx$", recursive = TRUE, full.names = TRUE)
    )
    cand_onnx <- cand_onnx[file.exists(cand_onnx)]
    if (length(cand_onnx) > 0L) dest_onnx <- cand_onnx[1]
  }

  dest_pt <- file.path(output_dir, paste0(base_stem, ".pt"))
  if (!file.exists(dest_pt)) {
    cand_pt <- c(
      file.path(results_dir, "weights", "best.pt"),
      file.path(dirname(results_dir), paste0(base_stem, ".pt")),
      list.files(results_dir, pattern = "best\\.pt$", recursive = TRUE, full.names = TRUE)
    )
    cand_pt <- cand_pt[file.exists(cand_pt)]
    if (length(cand_pt) > 0L) dest_pt <- cand_pt[1] else dest_pt <- NULL
  }

  find_plot <- function(patterns) {
    for (pat in patterns) {
      matches <- list.files(results_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
      if (length(matches) > 0L) return(matches[1])
    }
    NULL
  }

  plots_list <- list(
    results = find_plot(c("^results\\.png$", "^results\\.jpg$")),
    confusion_matrix_norm = find_plot(c("^confusion_matrix_normalized\\.png$", "^confusion_matrix_normalized\\.jpg$")),
    confusion_matrix = find_plot(c("^confusion_matrix\\.png$", "^confusion_matrix\\.jpg$")),
    pr_curve = find_plot(c("^(Mask|Box)?PR_curve\\.png$", "^(Mask|Box)?PR_curve\\.jpg$")),
    f1_curve = find_plot(c("^(Mask|Box)?F1_curve\\.png$", "^(Mask|Box)?F1_curve\\.jpg$")),
    p_curve = find_plot(c("^(Mask|Box)?P_curve\\.png$", "^(Mask|Box)?P_curve\\.jpg$")),
    r_curve = find_plot(c("^(Mask|Box)?R_curve\\.png$", "^(Mask|Box)?R_curve\\.jpg$")),
    labels = find_plot(c("^labels\\.(jpg|png)$")),
    labels_correlogram = find_plot(c("^labels_correlogram\\.(jpg|png)$"))
  )

  read_csv_safe <- function(filename) {
    p <- file.path(results_dir, filename)
    if (file.exists(p)) {
      tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
    } else {
      NULL
    }
  }

  cm_df <- read_csv_safe("confusion_matrix.csv")
  cm_norm_df <- read_csv_safe("confusion_matrix_normalized.csv")
  pr_df <- read_csv_safe("pr_curve.csv")
  f1_df <- read_csv_safe("f1_curve.csv")

  results_csv <- file.path(results_dir, "results.csv")
  metrics_df <- NULL
  best_metrics <- NULL
  epochs <- 0L
  elapsed_sec <- 0
  if (file.exists(results_csv)) {
    tryCatch({
      metrics_df <- utils::read.csv(results_csv, strip.white = TRUE, check.names = FALSE)
      names(metrics_df) <- trimws(names(metrics_df))
      epochs <- nrow(metrics_df)
      if ("time" %in% names(metrics_df)) {
        elapsed_sec <- round(max(metrics_df$time, na.rm = TRUE), 1)
      }
      map50_col <- grep("mAP50", names(metrics_df), value = TRUE)
      # Exclude mAP50-95 from mAP50 candidate
      map50_col <- map50_col[!grepl("95", map50_col)]
      if (length(map50_col) > 0L) {
        best_idx <- which.max(metrics_df[[map50_col[1]]])
        if (length(best_idx) > 0L) {
          best_metrics <- metrics_df[best_idx, , drop = FALSE]
        }
      }
    }, error = function(e) NULL)
  }

  is_seg <- grepl("seg", base_stem, ignore.case = TRUE)

  structure(
    list(
      model_file = dest_onnx,
      pt_file = dest_pt,
      results_dir = results_dir,
      plots = plots_list,
      metrics = metrics_df,
      best_metrics = best_metrics,
      confusion_matrix = cm_df,
      confusion_matrix_norm = cm_norm_df,
      pr_curve = pr_df,
      f1_curve = f1_df,
      epochs = epochs,
      elapsed_sec = elapsed_sec,
      elapsed_formatted = .format_hms(elapsed_sec),
      task = if (is_seg) "segment" else "detect",
      base_model = base_stem,
      call = match.call()
    ),
    class = c("yolo_train", "list")
  )
}

#' @title Print YOLO Training Summary
#' @description Displays a concise summary of YOLO model training results, performance, and saved artifacts.
#' @param x An object of class `yolo_train`.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.yolo_train <- function(x, ...) {
  cli::cli_h2("YOLO Model Training Summary (pliman)")
  cli::cli_ul(c(
    paste0("Task: {.val ", toupper(x$task), "}"),
    paste0("Base architecture: {.val ", x$base_model, "}"),
    paste0("Epochs trained: ", x$epochs),
    paste0("Training duration: ", x$elapsed_formatted)
  ))

  if (!is.null(x$best_metrics) && nrow(x$best_metrics) > 0L) {
    map50_col <- grep("mAP50", names(x$best_metrics), value = TRUE)
    map50_col <- map50_col[!grepl("95", map50_col)]
    map95_col <- grep("mAP50-95|mAP50.95", names(x$best_metrics), value = TRUE)
    ep_val <- x$best_metrics$epoch
    m50_val <- if (length(map50_col) > 0L) round(x$best_metrics[[map50_col[1]]], 3) else NULL
    m95_val <- if (length(map95_col) > 0L) round(x$best_metrics[[map95_col[1]]], 3) else NULL

    summary_txt <- paste0("Best epoch: ", ep_val)
    if (!is.null(m50_val)) {
      summary_txt <- paste0(summary_txt, " (mAP50: ", m50_val, if (!is.null(m95_val)) paste0(", mAP50-95: ", m95_val) else "", ")")
    }
    cli::cli_alert_success(summary_txt)
  }

  cli::cli_h3("Model Artifacts")
  cli::cli_ul(c(
    paste0("ONNX model: {.file ", x$model_file, "}"),
    if (!is.null(x$pt_file)) paste0("PyTorch weights: {.file ", x$pt_file, "}") else NULL,
    paste0("Diagnostics directory: {.path ", x$results_dir, "}")
  ))

  diag_items <- c(
    if (!is.null(x$metrics)) paste0("Loss & mAP curves (", nrow(x$metrics), " epochs)"),
    if (!is.null(x$confusion_matrix) || !is.null(x$confusion_matrix_norm)) "Confusion Matrix (native vector)",
    if (!is.null(x$pr_curve)) "Precision-Recall curve (native vector)",
    if (!is.null(x$f1_curve)) "F1-Confidence curve (native vector)"
  )
  if (length(diag_items) > 0L) {
    cli::cli_alert_info("Diagnostic evaluations: {.val {diag_items}}")
    cli::cli_alert_info("Use {.code plot(x)} for the 4-panel dashboard, or {.code plot(x, which = 'confusion_matrix')}, {.code plot(x, which = 'pr_curve')}, {.code plot(x, which = 'metrics')}.")
  }

  invisible(x)
}

# Helper to generate high-resolution native R vector plots for YOLO metrics
.plot_yolo_metrics_native <- function(df, task = "detect", bg = "#fafafa", ...) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    cli::cli_alert_warning("No metrics data available to plot.")
    return(invisible(NULL))
  }

  clean_cols <- gsub("[ /()-]", ".", tolower(names(df)))

  get_col <- function(patterns) {
    for (pat in patterns) {
      idx <- grep(pat, clean_cols, perl = TRUE)
      if (length(idx) > 0L) return(names(df)[idx[1]])
    }
    NULL
  }

  ep <- if ("epoch" %in% names(df)) df$epoch else seq_len(nrow(df))

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  c_seg_train <- get_col(c("train.*seg.*loss", "seg.*loss"))
  c_seg_val <- get_col(c("val.*seg.*loss"))
  has_seg <- !is.null(c_seg_train) || identical(task, "segment")

  c_box_train <- get_col(c("train.*box.*loss", "box.*loss"))
  c_box_val <- get_col(c("val.*box.*loss"))
  c_cls_train <- get_col(c("train.*cls.*loss", "cls.*loss"))
  c_cls_val <- get_col(c("val.*cls.*loss"))

  c_prec_b <- get_col(c("precision.*b", "precision"))
  c_rec_b <- get_col(c("recall.*b", "recall"))
  c_map50_b <- get_col(c("map50.*b", "map50(?!.*95)"))
  c_map95_b <- get_col(c("map50.*95.*b", "map50.*95"))

  c_map50_m <- get_col(c("map50.*m"))
  c_map95_m <- get_col(c("map50.*95.*m"))

  c_lr <- get_col(c("lr.*pg0", "^lr"))

  graphics::par(mfrow = c(2, 3), mar = c(3.8, 3.8, 2.5, 1.2), mgp = c(2.2, 0.7, 0), bg = bg)

  add_grid <- function() {
    graphics::grid(nx = NULL, ny = NULL, col = "#e5e7eb", lty = 1, lwd = 0.8)
  }

  # 1. Box Loss
  if (!is.null(c_box_train)) {
    tr_box <- df[[c_box_train]]
    vl_box <- if (!is.null(c_box_val)) df[[c_box_val]] else NULL
    y_max <- max(c(tr_box, vl_box), na.rm = TRUE) * 1.05
    graphics::plot(ep, tr_box, type = "n", xlab = "Epoch", ylab = "Box Loss",
                   main = "Box Loss (Train vs Val)", ylim = c(0, y_max), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, tr_box, col = "#2563EB", lwd = 2.2)
    if (!is.null(vl_box)) {
      graphics::lines(ep, vl_box, col = "#DC2626", lwd = 2.2, lty = 2)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    }
  }

  # 2. Classification or Segmentation Loss
  if (has_seg && !is.null(c_seg_train)) {
    tr_seg <- df[[c_seg_train]]
    vl_seg <- if (!is.null(c_seg_val)) df[[c_seg_val]] else NULL
    y_max <- max(c(tr_seg, vl_seg), na.rm = TRUE) * 1.05
    graphics::plot(ep, tr_seg, type = "n", xlab = "Epoch", ylab = "Mask Loss",
                   main = "Segmentation Loss", ylim = c(0, y_max), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, tr_seg, col = "#2563EB", lwd = 2.2)
    if (!is.null(vl_seg)) {
      graphics::lines(ep, vl_seg, col = "#DC2626", lwd = 2.2, lty = 2)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    }
  } else if (!is.null(c_cls_train)) {
    tr_cls <- df[[c_cls_train]]
    vl_cls <- if (!is.null(c_cls_val)) df[[c_cls_val]] else NULL
    y_max <- max(c(tr_cls, vl_cls), na.rm = TRUE) * 1.05
    graphics::plot(ep, tr_cls, type = "n", xlab = "Epoch", ylab = "Class Loss",
                   main = "Classification Loss", ylim = c(0, y_max), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, tr_cls, col = "#2563EB", lwd = 2.2)
    if (!is.null(vl_cls)) {
      graphics::lines(ep, vl_cls, col = "#DC2626", lwd = 2.2, lty = 2)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    }
  }

  # 3. Precision & Recall
  if (!is.null(c_prec_b) && !is.null(c_rec_b)) {
    p_val <- df[[c_prec_b]]
    r_val <- df[[c_rec_b]]
    graphics::plot(ep, p_val, type = "n", xlab = "Epoch", ylab = "Score",
                   main = "Precision & Recall", ylim = c(0, 1.05), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, p_val, col = "#059669", lwd = 2.2)
    graphics::lines(ep, r_val, col = "#D97706", lwd = 2.2, lty = 2)
    graphics::legend("bottomright", legend = c("Precision", "Recall"), col = c("#059669", "#D97706"),
                     lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
  }

  # 4. Validation mAP
  target_map50 <- if (has_seg && !is.null(c_map50_m)) c_map50_m else c_map50_b
  target_map95 <- if (has_seg && !is.null(c_map95_m)) c_map95_m else c_map95_b
  map_title <- if (has_seg && !is.null(c_map50_m)) "Mask Validation mAP" else "Validation mAP"

  if (!is.null(target_map50)) {
    m50 <- df[[target_map50]]
    m95 <- if (!is.null(target_map95)) df[[target_map95]] else NULL
    graphics::plot(ep, m50, type = "n", xlab = "Epoch", ylab = "mAP",
                   main = map_title, ylim = c(0, 1.05), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, m50, col = "#16A34A", lwd = 2.5)
    if (!is.null(m95)) graphics::lines(ep, m95, col = "#EA580C", lwd = 2.5, lty = 2)
    best_ep <- which.max(m50)
    if (length(best_ep) > 0L && !is.na(best_ep)) {
      graphics::points(ep[best_ep], m50[best_ep], pch = 21, bg = "#16A34A", col = "white", cex = 1.5)
      text_pos <- if (m50[best_ep] > 0.88) 1 else if (best_ep > length(ep) / 2) 2 else 4
      graphics::text(ep[best_ep], m50[best_ep],
                     labels = sprintf(" \u2605 Best: %.3f (ep %d)", m50[best_ep], ep[best_ep]),
                     pos = text_pos, font = 2, cex = 0.85, col = "#15803D")
    }
    leg_labels <- if (!is.null(m95)) c("mAP50", "mAP50-95") else "mAP50"
    leg_cols <- if (!is.null(m95)) c("#16A34A", "#EA580C") else "#16A34A"
    graphics::legend("bottomright", legend = leg_labels, col = leg_cols,
                     lty = if (!is.null(m95)) c(1, 2) else 1, lwd = 2.5, bty = "n", cex = 0.85)
  }

  # 5. F1 Score or Box mAP
  if (has_seg && !is.null(c_map50_b)) {
    m50_b <- df[[c_map50_b]]
    m95_b <- if (!is.null(c_map95_b)) df[[c_map95_b]] else NULL
    graphics::plot(ep, m50_b, type = "n", xlab = "Epoch", ylab = "mAP",
                   main = "Box Validation mAP", ylim = c(0, 1.05), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, m50_b, col = "#0284C7", lwd = 2.5)
    if (!is.null(m95_b)) graphics::lines(ep, m95_b, col = "#D97706", lwd = 2.5, lty = 2)
    best_b <- which.max(m50_b)
    if (length(best_b) > 0L && !is.na(best_b)) {
      graphics::points(ep[best_b], m50_b[best_b], pch = 21, bg = "#0284C7", col = "white", cex = 1.5)
      text_pos_b <- if (m50_b[best_b] > 0.88) 1 else if (best_b > length(ep) / 2) 2 else 4
      graphics::text(ep[best_b], m50_b[best_b],
                     labels = sprintf(" \u2605 Best: %.3f (ep %d)", m50_b[best_b], ep[best_b]),
                     pos = text_pos_b, font = 2, cex = 0.85, col = "#0369A1")
    }
    leg_labels <- if (!is.null(m95_b)) c("Box mAP50", "Box mAP50-95") else "Box mAP50"
    leg_cols <- if (!is.null(m95_b)) c("#0284C7", "#D97706") else "#0284C7"
    graphics::legend("bottomright", legend = leg_labels, col = leg_cols,
                     lty = if (!is.null(m95_b)) c(1, 2) else 1, lwd = 2.5, bty = "n", cex = 0.85)
  } else if (!is.null(c_prec_b) && !is.null(c_rec_b)) {
    p_val <- df[[c_prec_b]]
    r_val <- df[[c_rec_b]]
    f1 <- 2 * (p_val * r_val) / (p_val + r_val + 1e-6)
    graphics::plot(ep, f1, type = "n", xlab = "Epoch", ylab = "F1 Score",
                   main = "F1 Score Progression", ylim = c(0, 1.05), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, f1, col = "#7C3AED", lwd = 2.2)
    best_f1 <- which.max(f1)
    if (length(best_f1) > 0L && !is.na(best_f1)) {
      graphics::points(ep[best_f1], f1[best_f1], pch = 21, bg = "#7C3AED", col = "white", cex = 1.5)
      text_pos_f1 <- if (f1[best_f1] > 0.88) 1 else if (best_f1 > length(ep) / 2) 2 else 4
      graphics::text(ep[best_f1], f1[best_f1],
                     labels = sprintf(" Max: %.3f", f1[best_f1]),
                     pos = text_pos_f1, font = 2, cex = 0.85, col = "#6D28D9")
    }
  }

  # 6. Learning Rate Schedule
  if (!is.null(c_lr)) {
    lr_val <- df[[c_lr]]
    graphics::plot(ep, lr_val, type = "n", xlab = "Epoch", ylab = "Learning Rate",
                   main = "Learning Rate Schedule", ylim = c(0, max(lr_val, na.rm = TRUE) * 1.05),
                   cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, lr_val, col = "#64748B", lwd = 2.2)
  }

  invisible(df)
}

# Helper: Native R vector Confusion Matrix heatmap
.plot_yolo_confusion_matrix_native <- function(cm_norm = NULL, cm_counts = NULL, normalized = TRUE, bg = "#fafafa", ...) {
  if (is.null(cm_norm) && is.null(cm_counts)) {
    cli::cli_alert_warning("No confusion matrix data available to plot.")
    return(invisible(NULL))
  }

  # Helper to remove columns that are entirely NA (e.g. background column in classification)
  clean_cm_df <- function(df) {
    if (is.null(df) || !is.data.frame(df) || ncol(df) <= 2L) return(df)
    na_cols <- vapply(df[, -1, drop = FALSE], function(col) all(is.na(col)), logical(1))
    if (any(na_cols)) {
      df <- df[, c(TRUE, !na_cols), drop = FALSE]
    }
    df
  }

  if (!is.null(cm_norm)) cm_norm <- clean_cm_df(cm_norm)
  if (!is.null(cm_counts)) cm_counts <- clean_cm_df(cm_counts)

  if (is.null(cm_norm) && !is.null(cm_counts)) {
    cm_norm <- cm_counts
    mat_tmp <- as.matrix(cm_counts[, -1, drop = FALSE])
    rsums <- rowSums(mat_tmp, na.rm = TRUE)
    rsums[rsums == 0] <- 1
    cm_norm[, -1] <- round(mat_tmp / rsums, 4)
  }

  if (is.null(cm_counts) && !is.null(cm_norm)) {
    cm_counts <- cm_norm
  }

  row_labels <- as.character(cm_norm[[1]])
  col_labels <- names(cm_norm)[-1]
  mat_norm <- as.matrix(cm_norm[, -1, drop = FALSE])
  mat_counts <- as.matrix(cm_counts[, -1, drop = FALSE])

  mat_norm[is.na(mat_norm)] <- 0
  mat_counts[is.na(mat_counts)] <- 0

  nr <- nrow(mat_norm)
  nc <- ncol(mat_norm)

  op <- graphics::par(no.readonly = TRUE)
  if (identical(graphics::par("mfrow"), c(1L, 1L))) {
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = c(5, 6.5, 3.5, 2), bg = bg)
  }

  main_title <- if (isTRUE(normalized)) "Confusion Matrix (Normalized)" else "Confusion Matrix (Counts)"

  graphics::plot(1, type = "n", xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
                 axes = FALSE, xlab = "Predicted Class", ylab = "",
                 main = main_title, font.main = 2, cex.main = 1.15)
  graphics::title(ylab = "True Class", line = 4.2, font.lab = 1)

  blues <- grDevices::colorRampPalette(c("#f8fafc", "#bfdbfe", "#3b82f6", "#1d4ed8"))(100)

  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      val <- mat_norm[r, c]
      cnt <- mat_counts[r, c]
      y_pos <- nr - r + 1
      x_pos <- c

      col_idx <- if (is.na(val)) 1 else max(1, min(100, round(val * 99) + 1))
      bg_col <- blues[col_idx]

      graphics::rect(x_pos - 0.48, y_pos - 0.48, x_pos + 0.48, y_pos + 0.48,
                     col = bg_col, border = "#cbd5e1", lwd = 1.2)

      txt_col <- if (!is.na(val) && val > 0.55) "#ffffff" else "#0f172a"
      if (isTRUE(normalized)) {
        graphics::text(x_pos, y_pos + 0.08, labels = sprintf("%.2f", val), col = txt_col, font = 2, cex = 1.05)
        graphics::text(x_pos, y_pos - 0.18, labels = sprintf("(%d)", cnt), col = txt_col, cex = 0.82)
      } else {
        graphics::text(x_pos, y_pos + 0.08, labels = sprintf("%d", cnt), col = txt_col, font = 2, cex = 1.05)
        graphics::text(x_pos, y_pos - 0.18, labels = sprintf("(%.1f%%)", val * 100), col = txt_col, cex = 0.82)
      }
    }
  }

  cex_axis <- if (nc > 6) 0.75 else 0.9
  graphics::axis(1, at = seq_len(nc), labels = col_labels, font = 2, tick = FALSE, line = -0.5, cex.axis = cex_axis)
  graphics::axis(2, at = rev(seq_len(nr)), labels = row_labels, font = 2, tick = FALSE, line = -0.5, las = 2, cex.axis = cex_axis)

  invisible(cm_norm)
}

# Helper: Native R vector Precision-Recall Curve plot
.plot_yolo_pr_curve_native <- function(pr_df, mAP50 = NULL, bg = "#fafafa", ...) {
  if (is.null(pr_df) || !is.data.frame(pr_df) || nrow(pr_df) == 0L) {
    cli::cli_alert_warning("No Precision-Recall curve data available to plot.")
    return(invisible(NULL))
  }

  rec <- pr_df[[1]]
  cls_cols <- names(pr_df)[-1]

  op <- graphics::par(no.readonly = TRUE)
  if (identical(graphics::par("mfrow"), c(1L, 1L))) {
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = c(4.2, 4.2, 3, 1.5), mgp = c(2.2, 0.7, 0), bg = bg)
  }

  graphics::plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1.05),
                 xlab = "Recall", ylab = "Precision",
                 main = "Precision-Recall Curve", font.main = 2, cex.main = 1.15)
  graphics::grid(nx = NULL, ny = NULL, col = "#e2e8f0", lty = 1, lwd = 0.8)

  pal <- c("#2563eb", "#059669", "#d97706", "#dc2626", "#7c3aed", "#0891b2")

  leg_labels <- character(length(cls_cols))
  for (ci in seq_along(cls_cols)) {
    cn <- cls_cols[ci]
    prec <- pr_df[[cn]]
    col_k <- pal[((ci - 1) %% length(pal)) + 1]

    auc_val <- sum(diff(rec) * (prec[-1] + prec[-length(prec)]) / 2, na.rm = TRUE)
    auc_txt <- if (!is.null(mAP50) && length(cls_cols) == 1L) {
      sprintf("%.3f", mAP50)
    } else {
      sprintf("%.3f", auc_val)
    }
    leg_labels[ci] <- paste0(cn, " (mAP50: ", auc_txt, ")")

    graphics::polygon(c(rec, 1, 0), c(prec, 0, 0),
                      col = grDevices::adjustcolor(col_k, alpha.f = 0.15), border = NA)
    graphics::lines(rec, prec, col = col_k, lwd = 2.5)
  }

  graphics::legend("bottomleft", legend = leg_labels,
                   col = pal[seq_along(cls_cols)], lwd = 2.5, bty = "n", cex = 0.85)

  invisible(pr_df)
}

# Helper: Native R vector F1-Confidence Curve plot
.plot_yolo_f1_curve_native <- function(f1_df, bg = "#fafafa", ...) {
  if (is.null(f1_df) || !is.data.frame(f1_df) || nrow(f1_df) == 0L) {
    cli::cli_alert_warning("No F1-Confidence curve data available to plot.")
    return(invisible(NULL))
  }

  conf <- f1_df[[1]]
  cls_cols <- names(f1_df)[-1]

  op <- graphics::par(no.readonly = TRUE)
  if (identical(graphics::par("mfrow"), c(1L, 1L))) {
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = c(4.2, 4.2, 3, 1.5), mgp = c(2.2, 0.7, 0), bg = bg)
  }

  graphics::plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1.05),
                 xlab = "Confidence", ylab = "F1 Score",
                 main = "F1-Confidence Curve", font.main = 2, cex.main = 1.15)
  graphics::grid(nx = NULL, ny = NULL, col = "#e2e8f0", lty = 1, lwd = 0.8)

  pal <- c("#7c3aed", "#2563eb", "#059669", "#d97706", "#dc2626", "#0891b2")

  leg_labels <- character(length(cls_cols))
  for (ci in seq_along(cls_cols)) {
    cn <- cls_cols[ci]
    f1_val <- f1_df[[cn]]
    col_k <- pal[((ci - 1) %% length(pal)) + 1]

    best_idx <- which.max(f1_val)
    max_f1 <- if (length(best_idx) > 0L) f1_val[best_idx] else NA
    leg_labels[ci] <- paste0(cn, " (Max F1: ", round(max_f1, 3), " at conf ", round(conf[best_idx], 2), ")")

    graphics::lines(conf, f1_val, col = col_k, lwd = 2.5)
    if (!is.na(max_f1)) {
      graphics::points(conf[best_idx], max_f1, pch = 21, bg = col_k, col = "white", cex = 1.4)
    }
  }

  graphics::legend("bottomleft", legend = leg_labels,
                   col = pal[seq_along(cls_cols)], lwd = 2.5, bty = "n", cex = 0.85)

  invisible(f1_df)
}

# Helper: Native R 2x2 comprehensive diagnostic dashboard
.plot_yolo_all_dashboard_native <- function(x, bg = "#fafafa", ...) {
  df <- x$metrics
  has_metrics <- !is.null(df) && is.data.frame(df) && nrow(df) > 0L
  has_cm <- !is.null(x$confusion_matrix_norm) || !is.null(x$confusion_matrix)
  has_pr <- !is.null(x$pr_curve)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 5, 3, 1.5), mgp = c(2.2, 0.7, 0), bg = bg)

  add_grid <- function() {
    graphics::grid(nx = NULL, ny = NULL, col = "#e2e8f0", lty = 1, lwd = 0.8)
  }

  clean_cols <- if (has_metrics) gsub("[ /()-]", ".", tolower(names(df))) else character(0)
  get_col <- function(patterns) {
    for (pat in patterns) {
      idx <- grep(pat, clean_cols, perl = TRUE)
      if (length(idx) > 0L) return(names(df)[idx[1]])
    }
    NULL
  }

  ep <- if (has_metrics && "epoch" %in% names(df)) df$epoch else if (has_metrics) seq_len(nrow(df)) else NULL

  is_cls_task <- identical(x$task, "classify") ||
    (has_metrics && any(grepl("accuracy_top", names(df))))

  # Panel 1: Loss curves
  if (has_metrics) {
    if (is_cls_task) {
      c_tr_loss <- get_col(c("train.*loss", "^loss"))
      c_vl_loss <- get_col(c("val.*loss"))
      tr_loss <- if (!is.null(c_tr_loss)) df[[c_tr_loss]] else NULL
      vl_loss <- if (!is.null(c_vl_loss)) df[[c_vl_loss]] else NULL

      all_losses <- c(tr_loss, vl_loss)
      y_max <- if (length(all_losses) > 0L) max(all_losses, na.rm = TRUE) * 1.05 else 1

      graphics::plot(ep, tr_loss, type = "n", xlab = "Epoch", ylab = "Loss",
                     main = "Training & Validation Loss", ylim = c(0, y_max),
                     cex.main = 1.15, font.main = 2)
      add_grid()
      if (!is.null(tr_loss)) graphics::lines(ep, tr_loss, col = "#2563eb", lwd = 2.2)
      if (!is.null(vl_loss)) graphics::lines(ep, vl_loss, col = "#dc2626", lwd = 2.2, lty = 2)
      graphics::legend("topright", legend = c("Train Loss", "Val Loss"),
                       col = c("#2563eb", "#dc2626"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    } else {
      c_tr_box <- get_col(c("train.*box.*loss", "box.*loss"))
      c_vl_box <- get_col("val.*box.*loss")
      c_tr_cls <- get_col(c("train.*cls.*loss", "cls.*loss", "train.*seg.*loss"))
      c_vl_cls <- get_col(c("val.*cls.*loss", "val.*seg.*loss"))

      tr_box <- if (!is.null(c_tr_box)) df[[c_tr_box]] else NULL
      vl_box <- if (!is.null(c_vl_box)) df[[c_vl_box]] else NULL
      tr_cls <- if (!is.null(c_tr_cls)) df[[c_tr_cls]] else NULL
      vl_cls <- if (!is.null(c_vl_cls)) df[[c_vl_cls]] else NULL

      all_losses <- c(tr_box, vl_box, tr_cls, vl_cls)
      y_max <- if (length(all_losses) > 0L) max(all_losses, na.rm = TRUE) * 1.05 else 1

      graphics::plot(ep, tr_box, type = "n", xlab = "Epoch", ylab = "Loss",
                     main = "Training & Validation Loss", ylim = c(0, y_max),
                     cex.main = 1.15, font.main = 2)
      add_grid()
      if (!is.null(tr_box)) graphics::lines(ep, tr_box, col = "#2563eb", lwd = 2.2)
      if (!is.null(vl_box)) graphics::lines(ep, vl_box, col = "#dc2626", lwd = 2.2, lty = 2)
      if (!is.null(tr_cls)) graphics::lines(ep, tr_cls, col = "#059669", lwd = 2)
      if (!is.null(vl_cls)) graphics::lines(ep, vl_cls, col = "#d97706", lwd = 2, lty = 2)

      leg_names <- c(
        if (!is.null(tr_box)) "Box (Train)",
        if (!is.null(vl_box)) "Box (Val)",
        if (!is.null(tr_cls)) "Class (Train)",
        if (!is.null(vl_cls)) "Class (Val)"
      )
      leg_cols <- c(
        if (!is.null(tr_box)) "#2563eb",
        if (!is.null(vl_box)) "#dc2626",
        if (!is.null(tr_cls)) "#059669",
        if (!is.null(vl_cls)) "#d97706"
      )
      leg_ltys <- c(
        if (!is.null(tr_box)) 1,
        if (!is.null(vl_box)) 2,
        if (!is.null(tr_cls)) 1,
        if (!is.null(vl_cls)) 2
      )
      if (length(leg_names) > 0L) {
        graphics::legend("topright", legend = leg_names, col = leg_cols,
                         lty = leg_ltys, lwd = 2, bty = "n", cex = 0.8)
      }
    }
  }

  # Panel 2: Accuracy Progression (Classification) or mAP Progression (Detection/Segmentation)
  best_m50 <- NULL
  if (has_metrics) {
    if (is_cls_task) {
      c_top1 <- get_col(c("accuracy_top1", "top1", "accuracy"))
      c_top5 <- get_col(c("accuracy_top5", "top5"))
      top1_v <- if (!is.null(c_top1)) df[[c_top1]] else NULL
      top5_v <- if (!is.null(c_top5)) df[[c_top5]] else NULL

      if (!is.null(top1_v)) {
        graphics::plot(ep, top1_v, type = "n", xlab = "Epoch", ylab = "Accuracy",
                       main = "Validation Accuracy Progression", ylim = c(0, 1.05),
                       cex.main = 1.15, font.main = 2)
        add_grid()
        graphics::lines(ep, top1_v, col = "#16a34a", lwd = 2.5)
        if (!is.null(top5_v)) graphics::lines(ep, top5_v, col = "#ea580c", lwd = 2.5, lty = 2)

        best_ep <- which.max(top1_v)
        if (length(best_ep) > 0L && !is.na(best_ep)) {
          graphics::points(ep[best_ep], top1_v[best_ep], pch = 21, bg = "#16a34a", col = "white", cex = 1.5)
          text_pos <- if (top1_v[best_ep] > 0.88) 1 else 4
          graphics::text(ep[best_ep], top1_v[best_ep],
                         labels = sprintf(" \u2605 Best: %.1f%% (ep %d)", top1_v[best_ep] * 100, ep[best_ep]),
                         pos = text_pos, font = 2, cex = 0.85, col = "#15803d")
        }
        leg_acc_names <- c("Top-1 Acc", if (!is.null(top5_v)) "Top-5 Acc")
        leg_acc_cols <- c("#16a34a", if (!is.null(top5_v)) "#ea580c")
        graphics::legend("bottomright", legend = leg_acc_names, col = leg_acc_cols,
                         lty = if (!is.null(top5_v)) c(1, 2) else 1, lwd = 2.5, bty = "n", cex = 0.85)
      }
    } else {
      c_m50 <- get_col(c("map50.*b", "map50(?!.*95)", "map50"))
      c_m95 <- get_col(c("map50.*95.*b", "map50.*95"))
      m50 <- if (!is.null(c_m50)) df[[c_m50]] else NULL
      m95 <- if (!is.null(c_m95)) df[[c_m95]] else NULL

      if (!is.null(m50)) {
        graphics::plot(ep, m50, type = "n", xlab = "Epoch", ylab = "mAP",
                       main = "Validation mAP Progression", ylim = c(0, 1.05),
                       cex.main = 1.15, font.main = 2)
        add_grid()
        graphics::lines(ep, m50, col = "#16a34a", lwd = 2.5)
        if (!is.null(m95)) graphics::lines(ep, m95, col = "#ea580c", lwd = 2.5, lty = 2)
        best_ep <- which.max(m50)
        if (length(best_ep) > 0L && !is.na(best_ep)) {
          best_m50 <- m50[best_ep]
          graphics::points(ep[best_ep], m50[best_ep], pch = 21, bg = "#16a34a", col = "white", cex = 1.5)
          text_pos <- if (m50[best_ep] > 0.88) 1 else 4
          graphics::text(ep[best_ep], m50[best_ep],
                         labels = sprintf(" \u2605 Best: %.3f (ep %d)", m50[best_ep], ep[best_ep]),
                         pos = text_pos, font = 2, cex = 0.85, col = "#15803d")
        }
        leg_m_names <- c("mAP50", if (!is.null(m95)) "mAP50-95")
        leg_m_cols <- c("#16a34a", if (!is.null(m95)) "#ea580c")
        graphics::legend("bottomright", legend = leg_m_names, col = leg_m_cols,
                         lty = if (!is.null(m95)) c(1, 2) else 1, lwd = 2.5, bty = "n", cex = 0.85)
      }
    }
  }

  # Panel 3: Confusion Matrix
  if (has_cm) {
    .plot_yolo_confusion_matrix_native(x$confusion_matrix_norm, x$confusion_matrix, normalized = TRUE, bg = bg, ...)
  } else if (has_metrics) {
    c_prec <- get_col(c("precision.*b", "precision"))
    c_rec <- get_col(c("recall.*b", "recall"))
    if (!is.null(c_prec) && !is.null(c_rec)) {
      graphics::plot(ep, df[[c_prec]], type = "n", xlab = "Epoch", ylab = "Score",
                     main = "Precision & Recall Progression", ylim = c(0, 1.05), cex.main = 1.15, font.main = 2)
      add_grid()
      graphics::lines(ep, df[[c_prec]], col = "#059669", lwd = 2.2)
      graphics::lines(ep, df[[c_rec]], col = "#d97706", lwd = 2.2, lty = 2)
      graphics::legend("bottomright", legend = c("Precision", "Recall"), col = c("#059669", "#d97706"),
                       lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    }
  }

  # Panel 4: Precision-Recall Curve (or Learning Rate Schedule)
  if (has_pr) {
    .plot_yolo_pr_curve_native(x$pr_curve, mAP50 = best_m50, bg = bg, ...)
  } else if (!is.null(x$f1_curve)) {
    .plot_yolo_f1_curve_native(x$f1_curve, bg = bg, ...)
  } else if (has_metrics) {
    c_lr <- get_col(c("lr.*pg0", "^lr"))
    if (!is.null(c_lr)) {
      lr_v <- df[[c_lr]]
      graphics::plot(ep, lr_v, type = "n", xlab = "Epoch", ylab = "Learning Rate",
                     main = "Learning Rate Schedule", ylim = c(0, max(lr_v, na.rm = TRUE) * 1.05),
                     cex.main = 1.15, font.main = 2)
      add_grid()
      graphics::lines(ep, lr_v, col = "#64748b", lwd = 2.2)
    }
  }

  invisible(x)
}

#' @title Plot Diagnostic Evaluations and Performance Curves for YOLO Models
#' @name plot.yolo_train
#' @description
#' Generates high-resolution, publication-ready base R vector graphics for evaluating trained YOLO models.
#' Supports multi-panel overview dashboards, confusion matrices (raw counts and normalized rates),
#' precision-recall (PR) curves, confidence vs. F1-score curves, and epoch-wise loss/metric progressions.
#'
#' @details
#' \subsection{Diagnostic Evaluation Statistics and Mathematical Formulas}{
#'
#' \strong{1. Spatial Overlap (Intersection over Union, IoU):}
#' Given a predicted bounding box or segmentation mask \eqn{\mathcal{B}_{\text{pred}}} and a ground truth
#' annotation \eqn{\mathcal{B}_{\text{true}}}, the spatial overlap is quantified by the Jaccard index:
#' \deqn{\text{IoU} = \frac{|\mathcal{B}_{\text{pred}} \cap \mathcal{B}_{\text{true}}|}{|\mathcal{B}_{\text{pred}} \cup \mathcal{B}_{\text{true}}|}}{IoU = area(overlap) / area(union)}
#' A detection is designated as a True Positive (TP) if \eqn{\text{IoU} \ge \text{threshold}} (typically \eqn{0.50})
#' and the predicted class matches ground truth.
#'
#' \strong{2. Confusion Matrix with Background Class:}
#' Object detection models must detect both the presence and location of objects. Hence, detections can fail by predicting
#' an incorrect class, by hallucinating an object on the background, or by failing to detect a real object.
#' The evaluation matrix \eqn{M} of dimensions \eqn{(N_c + 1) \times (N_c + 1)} includes the background class:
#' \itemize{
#'   \item \strong{Cell \eqn{(i, j)} for \eqn{i, j \le N_c}:} Ground truth object of class \eqn{i} predicted as class \eqn{j}.
#'     Diagonal entries \eqn{M(i, i)} represent \strong{True Positives (TP)}.
#'   \item \strong{Last column \eqn{(i, N_c + 1)}:} Ground truth object of class \eqn{i} completely missed by the model
#'     (misclassified as background, i.e., \strong{False Negatives, FN}).
#'   \item \strong{Last row \eqn{(N_c + 1, j)}:} Background region falsely predicted as an object of class \eqn{j}
#'     (i.e., \strong{False Positives, FP}).
#' }
#' The row-normalized confusion matrix computes the empirical conditional probability \eqn{P(\hat{Y} = j \mid Y = i)}:
#' \deqn{M_{\text{norm}}(i, j) = P(\hat{Y} = j \mid Y = i) = \frac{M(i, j)}{\sum_{k=1}^{N_c + 1} M(i, k)}}{M_norm(i, j) = M(i, j) / sum_k M(i, k)}
#' The diagonal entry \eqn{M_{\text{norm}}(i, i)} corresponds directly to the per-class Sensitivity (Recall).
#' Conversely, the column-normalized entry represents the empirical precision \eqn{P(Y = i \mid \hat{Y} = j)}:
#' \deqn{M_{\text{prec}}(i, j) = P(Y = i \mid \hat{Y} = j) = \frac{M(i, j)}{\sum_{k=1}^{N_c + 1} M(k, j)}}{M_prec(i, j) = M(i, j) / sum_k M(k, j)}
#'
#' \strong{3. Precision and Recall:}
#' At any given confidence threshold \eqn{\tau \in [0, 1]}:
#' \deqn{\text{Precision}(\tau) = \frac{\text{TP}(\tau)}{\text{TP}(\tau) + \text{FP}(\tau)}}{Precision(tau) = TP(tau) / (TP(tau) + FP(tau))}
#' \deqn{\text{Recall}(\tau) = \frac{\text{TP}(\tau)}{\text{TP}(\tau) + \text{FN}(\tau)}}{Recall(tau) = TP(tau) / (TP(tau) + FN(tau))}
#' Precision measures the accuracy of positive predictions (minimizing false alarms), while Recall measures the completeness of coverage (minimizing missed objects).
#'
#' \strong{4. F1-Score vs. Confidence Curve:}
#' The F1-score is the harmonic mean of Precision and Recall as a function of the confidence threshold \eqn{\tau}:
#' \deqn{F_1(\tau) = 2 \times \frac{\text{Precision}(\tau) \times \text{Recall}(\tau)}{\text{Precision}(\tau) + \text{Recall}(\tau)} = \frac{2\,\text{TP}(\tau)}{2\,\text{TP}(\tau) + \text{FP}(\tau) + \text{FN}(\tau)}}{F1(tau) = 2 * (Precision * Recall) / (Precision + Recall)}
#' The peak of this curve identifies the optimal operating threshold \eqn{\tau^*} that maximizes the harmonic trade-off:
#' \deqn{\tau^* = \arg\max_{\tau} F_1(\tau)}{tau* = argmax F1(tau)}
#'
#' \strong{5. Precision-Recall (PR) Curve and Average Precision (AP):}
#' As the confidence threshold \eqn{\tau} ranges from 1 down to 0, recall monotonically increases from 0 to 1, while precision generally declines.
#' The interpolated precision \eqn{p_{\text{interp}}(r)} at recall level \eqn{r} is the maximum precision achieved at any recall \eqn{\tilde{r} \ge r}:
#' \deqn{p_{\text{interp}}(r) = \max_{\tilde{r} \ge r} p(\tilde{r})}{p_interp(r) = max_{r' >= r} p(r')}
#' Average Precision (AP) is the integral area under this interpolated curve, numerically approximated via the 101-point COCO evaluation grid:
#' \deqn{\text{AP} = \int_{0}^{1} p_{\text{interp}}(r)\,dr \approx \frac{1}{101} \sum_{k=0}^{100} p_{\text{interp}}\left(\frac{k}{100}\right)}{AP = integral_0^1 p_interp(r) dr = (1 / 101) * sum_{k=0}^100 p_interp(k / 100)}
#'
#' \strong{6. Mean Average Precision (mAP):}
#' \itemize{
#'   \item \strong{mAP@50 (mAP50):} Mean AP across all \eqn{N_c} target classes evaluated at an IoU threshold of \eqn{0.50}:
#'     \deqn{\text{mAP50} = \frac{1}{N_c} \sum_{c=1}^{N_c} \text{AP}_c^{(0.50)}}{mAP50 = (1 / Nc) * sum(AP_c)}
#'   \item \strong{mAP@50-95 (mAP50-95):} The primary metric of the COCO benchmark, computed as the unweighted mean of mAP across 10 IoU thresholds from \eqn{0.50} to \eqn{0.95} in steps of \eqn{0.05}:
#'     \deqn{\text{mAP50-95} = \frac{1}{10 \cdot N_c} \sum_{k=0}^{9} \sum_{c=1}^{N_c} \text{AP}_c^{(0.50 + 0.05k)}}{mAP50-95 = mean(mAP at IoU 0.50, 0.55, ..., 0.95)}
#' }
#' }
#'
#' \subsection{Interpretation Guide for Diagnostic Plots}{
#' \itemize{
#'   \item \strong{\code{which = "all"}}: Displays a 2x2 publication-ready dashboard combining:
#'     (1) Training & Validation Loss curves (detecting divergence, underfitting, or overfitting),
#'     (2) Validation mAP progression across epochs (with the best epoch highlighted with a star),
#'     (3) Normalized Confusion Matrix (identifying frequent misclassifications and background leakage), and
#'     (4) Precision-Recall curve with shaded area under the curve (quantifying overall detector precision).
#'   \item \strong{\code{which = "confusion_matrix"}}: Displays raw detection counts across all classes and background.
#'     High counts in the bottom row indicate excessive false positive background detections; high counts in the rightmost column indicate missed objects.
#'   \item \strong{\code{which = "confusion_matrix_norm"}}: Displays normalized conditional probabilities. The diagonal values represent per-class sensitivity (percentage of true objects correctly detected).
#'   \item \strong{\code{which = "pr_curve"}}: A steep curve maintaining precision close to 1.0 across high recall indicates a highly robust model with minimal false positives.
#'   \item \strong{\code{which = "f1_curve"}}: Used to choose the ideal confidence threshold for downstream inference in \code{\link{image_detect_dl}} or \code{\link{image_segment_dl}}.
#'   \item \strong{\code{which = "metrics"}}: Displays a 2x3 panel of detailed epoch-by-epoch loss curves (box, class, mask), precision, recall, and learning rate schedules.
#' }
#' }
#'
#' @param x An object of class `yolo_train` or character path/name of a model.
#' @param which Which plot to display. One of:
#'   * `"all"` (default): A 2x2 high-resolution diagnostic dashboard with Loss curves, Validation mAP progression, Confusion Matrix, and Precision-Recall curve.
#'   * `"metrics"`: High-resolution native R vector curves plotting loss, mAP, precision, and recall progression from `results.csv`.
#'   * `"confusion_matrix_norm"`: Normalized confusion matrix heatmap (counts + rates).
#'   * `"confusion_matrix"`: Raw count confusion matrix heatmap.
#'   * `"pr_curve"`: Precision-Recall (PR) curve with shaded area under the curve and mAP50.
#'   * `"f1_curve"`: F1 score vs. Confidence curve.
#'   * `"results"`: Alias for `"metrics"`.
#'   * `"labels"`: Dataset label distribution and bounding box correlogram (if available).
#' @param native Logical. When `TRUE` (default), renders crisp, resizable native R base graphic vector plots
#'   directly from metrics and CSV evaluation tables rather than displaying raster image files.
#' @param bg Background color for native plots. Defaults to `"#fafafa"`.
#' @param mar Plot margins passed to base graphics. Defaults to `0.5`.
#' @param ... Additional arguments passed to base graphics.
#' @return Invisibly returns `x` or plotted data frame.
#' @export
#' @seealso [yolo_plot()], [yolo_results()], [yolo_train()]
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Load diagnostics from a previously trained model
#' res <- yolo_results("capsulas_det")
#'
#' # Inspect training summary and plot 2x2 diagnostic dashboard
#' print(res)
#' plot(res)
#'
#' # Plot individual evaluation curves in high-resolution vector format:
#' plot(res, which = "confusion_matrix")
#' plot(res, which = "confusion_matrix_norm")
#' plot(res, which = "pr_curve")
#' plot(res, which = "f1_curve")
#' plot(res, which = "metrics")
#' }
#' @export
plot.yolo_train <- function(x,
                            which = c("all", "metrics", "confusion_matrix_norm",
                                      "confusion_matrix", "pr_curve", "f1_curve",
                                      "results", "labels"),
                            native = TRUE,
                            bg = "#fafafa",
                            mar = 0.5,
                            ...) {
  if (is.character(x)) {
    x <- yolo_results(x)
  }
  if (!inherits(x, "yolo_train")) {
    cli::cli_abort("Object must be of class {.code yolo_train} or a model name/path.")
  }

  which <- match.arg(which)

  # Dynamically lazy-load CSV evaluation data if missing from an existing / older in-memory object
  if (!is.null(x$results_dir) && dir.exists(x$results_dir)) {
    if (is.null(x$confusion_matrix)) {
      p <- file.path(x$results_dir, "confusion_matrix.csv")
      if (file.exists(p)) {
        x$confusion_matrix <- tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
      }
    }
    if (is.null(x$confusion_matrix_norm)) {
      p <- file.path(x$results_dir, "confusion_matrix_normalized.csv")
      if (file.exists(p)) {
        x$confusion_matrix_norm <- tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
      }
    }
    if (is.null(x$pr_curve)) {
      p <- file.path(x$results_dir, "pr_curve.csv")
      if (file.exists(p)) {
        x$pr_curve <- tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
      }
    }
    if (is.null(x$f1_curve)) {
      p <- file.path(x$results_dir, "f1_curve.csv")
      if (file.exists(p)) {
        x$f1_curve <- tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
      }
    }
    if (is.null(x$metrics)) {
      p <- file.path(x$results_dir, "results.csv")
      if (file.exists(p)) {
        x$metrics <- tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
      }
    }
  }

  show_single <- function(path, title = "") {
    if (is.null(path) || !file.exists(path)) {
      cli::cli_alert_warning("Plot data/file for {.val {which}} was not found in {.path {x$results_dir}}.")
      return(invisible(NULL))
    }
    img <- image_import(path)
    plot(img, main = title, ...)
    invisible(img)
  }

  if (which == "all") {
    if (isTRUE(native) && (!is.null(x$metrics) || !is.null(x$confusion_matrix_norm) || !is.null(x$pr_curve))) {
      return(invisible(.plot_yolo_all_dashboard_native(x, bg = bg, ...)))
    }
    # Fallback to raster plots if native not requested
    dashboard <- list(
      "Training Curves" = x$plots$results,
      "Confusion Matrix" = if (!is.null(x$plots$confusion_matrix_norm)) x$plots$confusion_matrix_norm else x$plots$confusion_matrix,
      "PR Curve" = x$plots$pr_curve,
      "F1 Curve" = x$plots$f1_curve
    )
    dashboard <- dashboard[!vapply(dashboard, function(p) is.null(p) || !file.exists(p), logical(1))]
    if (length(dashboard) == 0L) {
      cli::cli_alert_warning("No diagnostic plots found in {.path {x$results_dir}}.")
      return(invisible(NULL))
    }
    if (length(dashboard) == 1L) {
      return(show_single(dashboard[[1L]], names(dashboard)[1L]))
    }
    imgs <- lapply(dashboard, image_import)
    image_combine(imgs, labels = names(dashboard), mar = mar, ...)
    return(invisible(imgs))
  }

  if (which %in% c("metrics", "results")) {
    if (isTRUE(native) && !is.null(x$metrics) && nrow(x$metrics) > 0L) {
      return(invisible(.plot_yolo_metrics_native(x$metrics, task = x$task, bg = bg, ...)))
    }
    return(show_single(x$plots$results, "Training & Validation Progress"))
  }

  if (which == "confusion_matrix_norm") {
    if (isTRUE(native) && (!is.null(x$confusion_matrix_norm) || !is.null(x$confusion_matrix))) {
      return(invisible(.plot_yolo_confusion_matrix_native(x$confusion_matrix_norm, x$confusion_matrix, normalized = TRUE, bg = bg, ...)))
    }
    tgt <- if (!is.null(x$plots$confusion_matrix_norm)) x$plots$confusion_matrix_norm else x$plots$confusion_matrix
    return(show_single(tgt, "Normalized Confusion Matrix"))
  }

  if (which == "confusion_matrix") {
    if (isTRUE(native) && (!is.null(x$confusion_matrix) || !is.null(x$confusion_matrix_norm))) {
      return(invisible(.plot_yolo_confusion_matrix_native(x$confusion_matrix_norm, x$confusion_matrix, normalized = FALSE, bg = bg, ...)))
    }
    tgt <- if (!is.null(x$plots$confusion_matrix)) x$plots$confusion_matrix else x$plots$confusion_matrix_norm
    return(show_single(tgt, "Confusion Matrix"))
  }

  if (which == "pr_curve") {
    is_cls <- identical(x$task, "classify") ||
      (!is.null(x$metrics) && any(grepl("accuracy_top", names(x$metrics))))
    if (is_cls) {
      cli::cli_alert_info("Precision-Recall (PR) curves apply to Object Detection/Segmentation tasks and are not generated for Image Classification.")
      cli::cli_alert_info("Displaying Confusion Matrix instead.")
      return(invisible(.plot_yolo_confusion_matrix_native(x$confusion_matrix_norm, x$confusion_matrix, normalized = TRUE, bg = bg, ...)))
    }
    if (isTRUE(native) && !is.null(x$pr_curve)) {
      best_m50 <- NULL
      if (!is.null(x$best_metrics) && nrow(x$best_metrics) > 0L) {
        map_col <- grep("mAP50", names(x$best_metrics), value = TRUE)
        map_col <- map_col[!grepl("95", map_col)]
        if (length(map_col) > 0L) best_m50 <- x$best_metrics[[map_col[1]]]
      }
      return(invisible(.plot_yolo_pr_curve_native(x$pr_curve, mAP50 = best_m50, bg = bg, ...)))
    }
    return(show_single(x$plots$pr_curve, "Precision-Recall Curve"))
  }

  if (which == "f1_curve") {
    is_cls <- identical(x$task, "classify") ||
      (!is.null(x$metrics) && any(grepl("accuracy_top", names(x$metrics))))
    if (is_cls) {
      cli::cli_alert_info("F1-Confidence curves apply to Object Detection/Segmentation tasks and are not generated for Image Classification.")
      cli::cli_alert_info("Displaying Confusion Matrix instead.")
      return(invisible(.plot_yolo_confusion_matrix_native(x$confusion_matrix_norm, x$confusion_matrix, normalized = TRUE, bg = bg, ...)))
    }
    if (isTRUE(native) && !is.null(x$f1_curve)) {
      return(invisible(.plot_yolo_f1_curve_native(x$f1_curve, bg = bg, ...)))
    }
    return(show_single(x$plots$f1_curve, "F1-Confidence Curve"))
  }

  if (which == "labels") {
    tgt <- if (!is.null(x$plots$labels)) x$plots$labels else x$plots$labels_correlogram
    return(show_single(tgt, "Dataset Labels & Correlogram"))
  }

  invisible(x)
}

#' @title Plot Diagnostics for YOLO Training Objects or Model Names
#' @name yolo_plot
#' @description
#' Convenient wrapper function for [plot.yolo_train()]. Accepts either an existing
#' `yolo_train` S3 object, a model name string (e.g. `"capsulas_det"`), or a direct path
#' to a results folder.
#'
#' @details
#' This function simplifies diagnostic inspections by allowing users to directly pass model names
#' (e.g. `yolo_plot("my_model")`) without first having to call [yolo_results()].
#'
#' See [plot.yolo_train()] for full mathematical formulas, statistical definitions of IoU, Confusion Matrix,
#' Precision, Recall, F1-Score, PR Curve, mAP50, and mAP50-95, and plot interpretation guidelines.
#'
#' @param x An object of class `yolo_train`, or a character string specifying a model name or results directory.
#' @param ... Arguments passed directly to [plot.yolo_train()] (e.g. `which`, `native`, `bg`).
#'
#' @return Invisibly returns the `yolo_train` object or plotted data frame.
#' @export
#' @seealso [plot.yolo_train()], [yolo_results()], [yolo_train()]
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Plot full 2x2 vector diagnostic dashboard directly by model name
#' yolo_plot("capsulas_det")
#'
#' # Plot individual evaluation curves in high-resolution vector format:
#' yolo_plot("capsulas_det", which = "confusion_matrix")
#' yolo_plot("capsulas_det", which = "confusion_matrix_norm")
#' yolo_plot("capsulas_det", which = "pr_curve")
#' yolo_plot("capsulas_det", which = "f1_curve")
#' yolo_plot("capsulas_det", which = "metrics")
#' }
yolo_plot <- function(x, ...) {
  if (is.character(x)) {
    x <- yolo_results(x)
  }
  plot(x, ...)
}


