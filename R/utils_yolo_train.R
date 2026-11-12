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
.yolo_crop_to_objects <- function(img, boxes = NULL, polygons = NULL, margin = 0.05) {
  if (!is_image(img) && !is.array(img) && !is.matrix(img)) {
    return(list(img = img, boxes = boxes, polygons = polygons, cropped = FALSE))
  }
  dims <- dim(img)
  orig_w <- dims[1]
  orig_h <- dims[2]
  obj_x <- numeric(0)
  obj_y <- numeric(0)
  if (!is.null(boxes) && is.data.frame(boxes) && nrow(boxes) > 0L) {
    if (all(c("xmin", "ymin", "xmax", "ymax") %in% names(boxes))) {
      obj_x <- c(obj_x, boxes$xmin, boxes$xmax)
      obj_y <- c(obj_y, boxes$ymin, boxes$ymax)
    }
  }
  if (!is.null(polygons) && is.list(polygons) && length(polygons) > 0L) {
    for (k in seq_along(polygons)) {
      p <- polygons[[k]]
      if (is.null(p)) next
      if ((is.matrix(p) || is.data.frame(p)) && nrow(p) >= 3L) {
        p_mat <- as.matrix(p)
        obj_x <- c(obj_x, p_mat[, 1])
        obj_y <- c(obj_y, p_mat[, 2])
      } else if (is.vector(p) && length(p) >= 6L) {
        obj_x <- c(obj_x, p[c(TRUE, FALSE)] * orig_w)
        obj_y <- c(obj_y, p[c(FALSE, TRUE)] * orig_h)
      }
    }
  }
  obj_x <- obj_x[!is.na(obj_x)]
  obj_y <- obj_y[!is.na(obj_y)]
  if (length(obj_x) == 0L || length(obj_y) == 0L) {
    return(list(img = img, boxes = boxes, polygons = polygons, cropped = FALSE))
  }
  min_ox <- min(obj_x)
  max_ox <- max(obj_x)
  min_oy <- min(obj_y)
  max_oy <- max(obj_y)
  bw_span <- max(1.0, max_ox - min_ox)
  bh_span <- max(1.0, max_oy - min_oy)
  cmarg <- if (is.numeric(margin) && length(margin) >= 1L && !is.na(margin[1])) as.numeric(margin[1]) else 0.05
  pad_x <- max(1.0, bw_span * cmarg)
  pad_y <- max(1.0, bh_span * cmarg)
  crop_xmin <- max(1L, as.integer(floor(min_ox - pad_x)))
  crop_xmax <- min(as.integer(orig_w), as.integer(ceiling(max_ox + pad_x)))
  crop_ymin <- max(1L, as.integer(floor(min_oy - pad_y)))
  crop_ymax <- min(as.integer(orig_h), as.integer(ceiling(max_oy + pad_y)))
  if (crop_xmax <= crop_xmin || crop_ymax <= crop_ymin) {
    return(list(img = img, boxes = boxes, polygons = polygons, cropped = FALSE))
  }
  if (crop_xmin == 1L && crop_xmax == orig_w && crop_ymin == 1L && crop_ymax == orig_h) {
    return(list(img = img, boxes = boxes, polygons = polygons, cropped = FALSE))
  }
  cropped_img <- image_crop(img, width = crop_xmin:crop_xmax, height = crop_ymin:crop_ymax, plot = FALSE)
  new_w <- dim(cropped_img)[1]
  new_h <- dim(cropped_img)[2]
  new_boxes <- boxes
  if (!is.null(new_boxes) && is.data.frame(new_boxes) && nrow(new_boxes) > 0L) {
    if (all(c("xmin", "ymin", "xmax", "ymax") %in% names(boxes))) {
      new_boxes$xmin <- pmax(0.0, pmin(as.numeric(new_w), new_boxes$xmin - crop_xmin + 1.0))
      new_boxes$xmax <- pmax(0.0, pmin(as.numeric(new_w), new_boxes$xmax - crop_xmin + 1.0))
      new_boxes$ymin <- pmax(0.0, pmin(as.numeric(new_h), new_boxes$ymin - crop_ymin + 1.0))
      new_boxes$ymax <- pmax(0.0, pmin(as.numeric(new_h), new_boxes$ymax - crop_ymin + 1.0))
    }
  }
  new_polygons <- polygons
  if (!is.null(new_polygons) && is.list(new_polygons) && length(new_polygons) > 0L) {
    for (k in seq_along(new_polygons)) {
      p <- new_polygons[[k]]
      if (is.null(p)) next
      if ((is.matrix(p) || is.data.frame(p)) && nrow(p) >= 3L) {
        p_mat <- as.matrix(p)
        p_mat[, 1] <- pmax(0.0, pmin(as.numeric(new_w), p_mat[, 1] - crop_xmin + 1.0))
        p_mat[, 2] <- pmax(0.0, pmin(as.numeric(new_h), p_mat[, 2] - crop_ymin + 1.0))
        new_polygons[[k]] <- p_mat
      } else if (is.vector(p) && length(p) >= 6L) {
        px_orig <- p[c(TRUE, FALSE)] * orig_w
        py_orig <- p[c(FALSE, TRUE)] * orig_h
        px_new <- pmax(0.0, pmin(1.0, (px_orig - crop_xmin + 1.0) / new_w))
        py_new <- pmax(0.0, pmin(1.0, (py_orig - crop_ymin + 1.0) / new_h))
        new_polygons[[k]] <- as.vector(rbind(px_new, py_new))
      }
    }
  }
  list(
    img = cropped_img,
    boxes = new_boxes,
    polygons = new_polygons,
    crop_bounds = c(xmin = crop_xmin, xmax = crop_xmax, ymin = crop_ymin, ymax = crop_ymax),
    cropped = TRUE
  )
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
#'     for rapid dataset expansion. If a YOLO detection model is passed with \code{task = "segment"},
#'     \code{pliman} automatically pairs the detected bounding boxes with Meta SAM 2.1 to extract
#'     precise sub-pixel instance segmentation polygons.
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
#' @param dir Output directory path for the YOLO dataset. If `NULL` (default when `img` is a directory path),
#'   the dataset structure (`images/`, `labels/`, `data.yaml`) is generated directly inside `img`.
#'   Specify a directory name only if you wish to create a separate dataset folder.
#' @param task Task type: `"detect"` (default) for bounding-box detection or `"segment"`
#'   for polygon instance segmentation.
#' @param class_names Character vector of class names (e.g. `c("capsula")` or `c("grain", "defect")`).
#'   Defaults to `c("object")`.
#' @param class_id Numeric class ID or vector of IDs. Defaults to `NULL`, which automatically
#'   assigns sequential IDs (0, 1, ...) matching `class_names`.
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
#' @param filter_edge Logical. Whether to automatically filter out spurious border-cut artifacts touching image edges (default: `FALSE`).
#' @param crop_to_objects Logical. If `TRUE`, crops the image to the bounding box containing all detected objects plus a margin (default: `FALSE`).
#' @param crop_margin Numeric. Margin fraction added around the object bounding box when `crop_to_objects = TRUE` (default: `0.05`, i.e., 5%).
#' @param box_margin Safety margin (padding proportion) added around detected bounding boxes before segmenting with Meta SAM 2.1 (e.g. `0.05` expands 5% on each side). Prevents masks from cutting off object borders. Defaults to `0.05` (5%).
#' @param sahi Logical. If `TRUE`, applies Slicing-Aided Hyper Inference (SAHI) during YOLO pre-labeling,
#'   dividing high-resolution images into overlapping slices to detect small and crowded objects with high precision. Defaults to `FALSE`.
#' @param slice_height,slice_width Integer dimensions (in pixels) for SAHI slices. Defaults to `640`.
#' @param overlap_height_ratio,overlap_width_ratio Overlap ratio between adjacent slices along height and width. Defaults to `0.20` (20%).
#' @param match_metric Match metric used for merging slice predictions across boundaries: `"IOS"` (Intersection over Smaller, default) or `"IOU"`.
#' @param match_threshold Threshold for match metric to merge overlapping slice predictions. Defaults to `0.50`.
#' @param full_image_inference Logical. If `TRUE` (default), runs an additional full-image inference pass alongside slice predictions to ensure large objects spanning across multiple slices are captured.
#' @param sahi_if Conditional SAHI trigger. Can be `NULL` (default, disabled), `TRUE` (activates SAHI if average relative box area <= 1% of total image area), or a numeric threshold in `(0, 1]` (e.g. `0.005` for 0.5% or `0.02` for 2%). When specified, standard inference runs first; if the objects are small relative to the canvas (or if no objects are detected when `sahi_if_empty = TRUE`), SAHI is automatically triggered to recover fine details and small instances.
#' @param sahi_metric Character metric to aggregate relative bounding box areas when checking `sahi_if`. One of `"mean"` (default), `"median"`, or `"min"`.
#' @param sahi_if_empty Logical. If `TRUE` (default), automatically triggers SAHI when standard inference finds 0 objects on an image larger than the slice dimensions.
#' @param show_similarity Logical. If `TRUE`, displays a side-by-side 2-panel preview showing the detected
#'   objects on the left and the PerSAM Cosine Similarity Map ("Matriz de Coincidencia") on the right. Defaults to `FALSE`.
#' @param split_ratio Numeric vector of length 3 for `c(train, val, test)` split.
#'   Defaults to `c(0.70, 0.20, 0.10)`. If a vector of length 2 is provided,
#'   `c(train, val)` is assumed. If a single number is given, it is treated as the training proportion.
#' @param n_points Target number of polygon vertices per object for segmentation
#'   (default `40`). Points are sampled evenly along the contour to keep label files
#'   compact and clean.
#' @param min_area Minimum object area (in pixels) to include in the dataset. Defaults to `15`.
#' @param review Logical or integer. If `TRUE` (default), writes annotated review images to a `review/` subfolder. An integer can be provided (e.g. `review = 100`) to randomly sample that number of images across the dataset for review previews. Defaults to `TRUE`.
#'   subfolder in `dir` for visual inspection of all images. An integer can be provided (e.g. `review = 30`)
#'   to limit the number of review images generated.
#' @param rainbow Logical. If `TRUE` (default), plots instances in distinct rainbow colors during interactive review.
#' @param show_bbox Logical. Whether to draw bounding boxes on review previews. Defaults to `TRUE`. Supported alias: `bbox`.
#' @param show_polygon Logical. Whether to draw segmentation polygons on review previews. Defaults to `TRUE`. Supported alias: `polygon`.
#' @param show_points Logical. If `TRUE`, draws points at the centroid of each object using the same color as the box. Defaults to `FALSE`. Supported alias: `points`.
#' @param show_labels Logical. Whether to display class name on labels/badges during review previews. Defaults to `TRUE`. Supported aliases: `show_class`, `labels`.
#' @param show_conf Logical. Whether to display confidence score on labels/badges during review previews. Defaults to `TRUE`. Supported alias: `show_scores`.
#' @param bbox,polygon,points,labels Backward-compatible aliases for `show_bbox`, `show_polygon`, `show_points`, and `show_labels`.
#' @param badge Logical. Whether to draw a filled background badge behind labels during interactive review. Defaults to `TRUE`.
#' @param engine Computing device engine: `"cpu"` (default) or `"gpu"`.
#' @param device_id GPU device ID (default `-1`).
#' @param lower_size Optional numeric value to filter out small detected objects.
#'   When `size_mode = "relative"` (default), objects with a bounding box area
#'   less than `lower_size` times the mean area of all candidate boxes in that image are removed
#'   (e.g., `lower_size = 0.2` removes small debris, dust, or grains smaller than 20% of the average object size).
#'   When `size_mode = "absolute"`, objects with an area (in pixels) below `lower_size` are removed.
#'   Defaults to `NULL` (no lower size filter).
#' @param upper_size Optional numeric value to filter out abnormally large detected objects.
#'   When `size_mode = "relative"` (default), objects with a bounding box area
#'   greater than `upper_size` times the mean area of all candidate boxes in that image are removed
#'   (e.g., `upper_size = 2.0` removes objects larger than 200% of the average object size).
#'   When `size_mode = "absolute"`, objects with an area (in pixels) above `upper_size` are removed.
#'   Defaults to `NULL` (no upper size filter).
#' @param size_mode Character string specifying how `lower_size` and `upper_size` are interpreted.
#'   `"relative"` (default) treats values as proportions/multiples of the mean bounding box area.
#'   `"absolute"` treats values as area in pixels.
#' @param move Logical. If `TRUE` (default), moves the source images directly into the dataset's `images/` split directories instead of copying them.
#' @param overwrite Logical. If `TRUE` (default), overwrites existing dataset files in `dir`.
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs the dataset export as an RStudio Background Job in the Jobs tab, keeping the current R session free. Defaults to `FALSE`.
#' @param verbose Logical. If `TRUE` (default), displays live progress information.
#' @param ... Additional arguments passed directly to [analyze_objects()] when
#'   `model = NULL` and `annotations = "auto"` (e.g. `index = "B"`, `watershed = TRUE`, `tolerance = 1`).
#'
#' @return An invisible list containing:
#'   * `dir`: Root directory of the dataset.
#'   * `yaml_file`: Full path to the generated `data.yaml`.
#'   * `train_images`: Vector of paths to training images.
#'   * `val_images`: Vector of paths to validation images.
#'   * `test_images`: Vector of paths to test images (if test split was generated).
#'   * `review_dir`: Path to the review directory containing visual previews (if `review = TRUE`).
#'   * `n_objects`: Total number of annotated objects across the dataset.
#' @export
#' @examplesIf interactive()
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
yolo_dataset_export <- function(img,
                                model = NULL,
                                annotations = "auto",
                                dir = NULL,
                                task = c("detect", "segment"),
                                class_names = c("object"),
                                class_id = NULL,
                                prompt = NULL,
                                prompt_mode = c("first", "each"),
                                interactive = is.null(prompt) && !is.null(model) && is.character(model) && identical(tolower(model[1]), "persam") && base::interactive(),
                                sim_threshold = 0.50,
                                min_dist = 20,
                                box_threshold = 0.25,
                                text_threshold = 0.25,
                                iou_threshold = 0.50,
                                feat_res = 256,
                                max_objects = NULL,
                                split_ratio = c(0.70, 0.20, 0.10),
                                n_points = 40L,
                                min_area = 15,
                                filter_edge = FALSE,
                                crop_to_objects = FALSE,
                                crop_margin = 0.05,
                                box_margin = 0.05,
                                sahi = FALSE,
                                slice_height = 640,
                                slice_width = 640,
                                overlap_height_ratio = 0.20,
                                overlap_width_ratio = 0.20,
                                match_metric = c("IOS", "IOU"),
                                match_threshold = 0.50,
                                full_image_inference = TRUE,
                                sahi_if = NULL,
                                sahi_metric = c("mean", "median", "min"),
                                sahi_if_empty = TRUE,
                                lower_size = NULL,
                                upper_size = NULL,
                                size_mode = c("relative", "absolute"),
                                show_similarity = FALSE,
                                superres_map = FALSE,
                                review = TRUE,
                                rainbow = TRUE,
                                show_bbox = TRUE,
                                show_polygon = TRUE,
                                show_points = FALSE,
                                show_labels = TRUE,
                                show_conf = TRUE,
                                bbox = show_bbox,
                                polygon = show_polygon,
                                points = show_points,
                                labels = show_labels,
                                badge = TRUE,
                                engine = c("gpu", "cpu"),
                                device_id = -1,
                                move = TRUE,
                                overwrite = TRUE,
                                as_job = FALSE,
                                verbose = TRUE,
                                ...) {
  dots <- list(...)
  if ("as_job" %in% names(dots)) as_job <- isTRUE(dots$as_job)

  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with dataset export in current session.")
    } else {
      mc <- match.call(expand.dots = TRUE)
      eval_env <- parent.frame()
      evaluated_args <- list()
      arg_names <- setdiff(names(mc), c("", "as_job"))
      for (arg_name in arg_names) {
        arg_val <- tryCatch(eval(mc[[arg_name]], envir = eval_env), error = function(e) mc[[arg_name]])
        evaluated_args[[arg_name]] <- arg_val
      }
      evaluated_args$as_job <- FALSE

      if (isTRUE(evaluated_args$interactive) && is.null(evaluated_args$prompt) &&
          !is.null(evaluated_args$model) && is.character(evaluated_args$model) &&
          identical(tolower(evaluated_args$model[1]), "persam")) {
        cli::cli_abort("Interactive prompt selection ({.code interactive = TRUE}) cannot be run as an RStudio background job. Please provide {.arg prompt} coordinates or run with {.code as_job = FALSE}.")
      }
      evaluated_args$interactive <- FALSE

      job_dir <- file.path(tempdir(), "yolo_jobs")
      if (!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("yolo_export_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_export_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (!is.null(evaluated_args$dir) && is.character(evaluated_args$dir) && nzchar(evaluated_args$dir[1])) {
          basename(evaluated_args$dir[1])
        } else if (!is.null(evaluated_args$img) && is.character(evaluated_args$img) && nzchar(evaluated_args$img[1])) {
          basename(evaluated_args$img[1])
        } else {
          "dataset"
        }
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Dataset Export: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_export()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Dataset Export - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_export, job_info$args)\n",
        "target_dir <- if (!is.null(res$dir)) res$dir else job_info$args$dir\n",
        "if (!is.null(target_dir) && dir.exists(target_dir)) {\n",
        "  save_path <- file.path(target_dir, 'dataset_export.rds')\n",
        "  saveRDS(res, save_path)\n",
        "  cat('\\n[OK] Dataset S3 object successfully saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Dataset Export Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Dataset Export (pliman)")
      cli::cli_alert_success("Dataset export has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target: {.val {target_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track live export progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  if (!missing(bbox)) show_bbox <- bbox
  if (!missing(polygon)) show_polygon <- polygon
  if (!missing(points)) show_points <- points
  if (!missing(labels)) show_labels <- labels
  if ("bbox" %in% names(dots)) show_bbox <- isTRUE(dots$bbox)
  if ("show_bbox" %in% names(dots)) show_bbox <- isTRUE(dots$show_bbox)
  if ("show_boxes" %in% names(dots)) show_bbox <- isTRUE(dots$show_boxes)
  if ("polygon" %in% names(dots)) show_polygon <- isTRUE(dots$polygon)
  if ("show_polygon" %in% names(dots)) show_polygon <- isTRUE(dots$show_polygon)
  if ("points" %in% names(dots)) show_points <- isTRUE(dots$points)
  if ("point" %in% names(dots)) show_points <- isTRUE(dots$point)
  if ("show_points" %in% names(dots)) show_points <- isTRUE(dots$show_points)
  if ("show_point" %in% names(dots)) show_points <- isTRUE(dots$show_point)
  if ("show_class" %in% names(dots)) show_labels <- isTRUE(dots$show_class)
  if ("show_labels" %in% names(dots)) show_labels <- isTRUE(dots$show_labels)
  if ("labels" %in% names(dots)) show_labels <- isTRUE(dots$labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("show_conf" %in% names(dots)) show_conf <- isTRUE(dots$show_conf)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)
  show_class <- show_labels
  if ("crop_to_objects" %in% names(dots)) crop_to_objects <- isTRUE(dots$crop_to_objects)
  if ("crop_margin" %in% names(dots)) crop_margin <- as.numeric(dots$crop_margin[1])
  if (is.numeric(crop_to_objects) && length(crop_to_objects) == 1L && crop_to_objects > 0) {
    crop_margin <- as.numeric(crop_to_objects)
    crop_to_objects <- TRUE
  }
  if ("box_margin" %in% names(dots)) box_margin <- as.numeric(dots$box_margin[1])
  if ("margin" %in% names(dots)) box_margin <- as.numeric(dots$margin[1])
  if ("sahi" %in% names(dots)) sahi <- isTRUE(dots$sahi)
  if ("sahi_if" %in% names(dots)) sahi_if <- dots$sahi_if
  sahi_metric <- match.arg(sahi_metric)
  if ("sahi_metric" %in% names(dots)) sahi_metric <- dots$sahi_metric[1]
  if ("sahi_if_empty" %in% names(dots)) sahi_if_empty <- isTRUE(dots$sahi_if_empty)
  if ("slice_height" %in% names(dots)) slice_height <- as.integer(dots$slice_height[1])
  if ("slice_width" %in% names(dots)) slice_width <- as.integer(dots$slice_width[1])
  if ("overlap_height_ratio" %in% names(dots)) overlap_height_ratio <- as.numeric(dots$overlap_height_ratio[1])
  if ("overlap_width_ratio" %in% names(dots)) overlap_width_ratio <- as.numeric(dots$overlap_width_ratio[1])
  match_metric <- match.arg(match_metric)
  if ("match_metric" %in% names(dots)) match_metric <- dots$match_metric[1]
  if ("match_threshold" %in% names(dots)) match_threshold <- as.numeric(dots$match_threshold[1])
  if ("full_image_inference" %in% names(dots)) full_image_inference <- isTRUE(dots$full_image_inference)
  dilate_val <- if ("dilate" %in% names(dots)) as.integer(dots$dilate[1]) else 0L
  erode_val  <- if ("erode" %in% names(dots)) as.integer(dots$erode[1]) else 0L
  if ("lower_size" %in% names(dots)) lower_size <- dots$lower_size
  if ("upper_size" %in% names(dots)) upper_size <- dots$upper_size
  if ("size_mode" %in% names(dots)) size_mode <- dots$size_mode
  if (!is.null(size_mode)) size_mode <- match.arg(size_mode, c("relative", "absolute"))
  if ("train_prop" %in% names(dots)) {
    tp <- as.numeric(dots$train_prop[1])
    split_ratio <- c(tp, 1.0 - tp, 0.0)
  }
  extra_args <- dots[setdiff(names(dots), c(
    "train_prop", "show_labels", "show_scores", "show_class", "show_conf", "label_box",
    "crop_to_objects", "crop_margin", "box_margin", "margin", "show_bbox", "bbox", "show_boxes",
    "show_polygon", "polygon", "show_points", "points", "point", "show_point", "labels",
    "lower_size", "upper_size", "size_mode", "sahi", "sahi_if", "sahi_metric", "sahi_if_empty", "slice_height", "slice_width",
    "overlap_height_ratio", "overlap_width_ratio", "match_metric", "match_threshold", "full_image_inference",
    "as_job"
  ))]
  split_ratio <- as.numeric(split_ratio)
  if (length(split_ratio) == 1L) {
    split_ratio <- c(split_ratio, max(0, 1 - split_ratio), 0)
  } else if (length(split_ratio) == 2L) {
    split_ratio <- c(split_ratio, max(0, 1 - sum(split_ratio)))
  }
  split_ratio <- split_ratio / sum(split_ratio)
  splits <- c("train", "val")
  if (split_ratio[3] > 0) splits <- c(splits, "test")
  task <- match.arg(task)
  engine <- match.arg(engine)
  prompt_mode <- match.arg(prompt_mode)
  if (is.null(dir) || !nzchar(dir)) {
    dir <- if (is.character(img) && length(img) == 1L && dir.exists(img)) img else "yolo_dataset"
  }
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
  if (is.null(class_id) || length(class_id) == 0L) {
    class_id <- seq_along(class_names) - 1L
  } else if (length(class_id) == 1L && length(class_names) > 1L) {
    class_id <- seq_along(class_names) - 1L
  } else if (length(class_id) < length(class_names)) {
    class_id <- c(class_id, seq(from = max(class_id) + 1L, length.out = length(class_names) - length(class_id)))
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
  test_img_dir  <- file.path(dir, "images", "test")
  train_lbl_dir <- file.path(dir, "labels", "train")
  val_lbl_dir   <- file.path(dir, "labels", "val")
  test_lbl_dir  <- file.path(dir, "labels", "test")
  for (sp in splits) {
    dir.create(file.path(dir, "images", sp), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir, "labels", sp), recursive = TRUE, showWarnings = FALSE)
  }
  # Resolve images input
  img_list <- list()
  img_names <- character()
  img_src_files <- character()
  input_orig_dir <- NULL
  if (is.character(img) && length(img) == 1L && dir.exists(img)) {
    files <- list.files(img, pattern = "\\.(jpg|jpeg|png|tif|tiff|bmp)$", ignore.case = TRUE, full.names = TRUE, recursive = FALSE)
    files <- files[!grepl("^review_", basename(files))]
    if (length(files) == 0L) {
      files <- list.files(img, pattern = "\\.(jpg|jpeg|png|tif|tiff|bmp)$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
      files <- files[!grepl("(_aug\\d+|_tile_\\d+_\\d+|^review_)", basename(files))]
    }
    if (length(files) == 0L) {
      cli::cli_abort("No valid images found in directory {.path {img}}.")
    }
    img_names <- tools::file_path_sans_ext(basename(files))
    img_list <- as.list(files)
    img_src_files <- normalizePath(files, winslash = "/", mustWork = FALSE)
    input_orig_dir <- normalizePath(img, winslash = "/", mustWork = FALSE)
  } else if (is.character(img)) {
    img_names <- tools::file_path_sans_ext(basename(img))
    img_list <- as.list(img)
    img_src_files <- normalizePath(img, winslash = "/", mustWork = FALSE)
    input_orig_dir <- normalizePath(dirname(img[1]), winslash = "/", mustWork = FALSE)
  } else if (inherits(img, "image") || inherits(img, "Image")) {
    img_list <- list(img)
    img_names <- if (!is.null(names(img)) && nzchar(names(img)[1])) names(img)[1] else "img_1"
    img_src_files <- rep(NA_character_, 1L)
  } else if (is.list(img)) {
    img_list <- img
    img_names <- names(img)
    if (is.null(img_names) || any(!nzchar(img_names))) {
      img_names <- paste0("img_", seq_along(img))
    }
    img_src_files <- rep(NA_character_, length(img_list))
  } else {
    cli::cli_abort("Invalid {.arg img} argument. Must be an image, list of images, or directory path.")
  }
  move <- isTRUE(move)
  logged_orig_paths <- character()
  logged_dest_imgs  <- character()
  logged_dest_lbls  <- character()
  logged_rev_imgs   <- character()
  num_imgs <- length(img_list)
  do_review <- isTRUE(review) || (is.numeric(review) && length(review) == 1L && review[1] > 0)
  review_indices <- integer(0)
  review_dir <- file.path(dir, "review")
  if (do_review) {
    if (is.numeric(review) && length(review) == 1L && review[1] > 0) {
      n_rev <- min(num_imgs, as.integer(review[1]))
      review_indices <- sort(sample(seq_len(num_imgs), n_rev))
    } else {
      review_indices <- seq_len(num_imgs)
    }
    dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (num_imgs == 0L) {
    cli::cli_abort("No images to export.")
  }
  # Train / Val / Test Split
  if (num_imgs == 1L) {
    split_assign <- c("train")
  } else if (num_imgs == 2L) {
    split_assign <- c("train", "val")
  } else {
    set.seed(42)
    shuffled_idx <- sample(seq_len(num_imgs))
    n_tr  <- max(1L, round(num_imgs * split_ratio[1]))
    n_val <- max(1L, round(num_imgs * split_ratio[2]))
    split_assign <- character(num_imgs)
    split_assign[shuffled_idx[seq_len(min(num_imgs, n_tr))]] <- "train"
    val_indices <- shuffled_idx[(n_tr + 1L):min(num_imgs, n_tr + n_val)]
    split_assign[val_indices] <- "val"
    if (num_imgs > (n_tr + n_val) && split_ratio[3] > 0) {
      test_indices <- shuffled_idx[(n_tr + n_val + 1L):num_imgs]
      split_assign[test_indices] <- "test"
    } else {
      split_assign[split_assign == ""] <- "val"
    }
  }
  total_objects_annotated <- 0L
  train_images_saved <- character()
  val_images_saved   <- character()
  test_images_saved  <- character()
  if (isTRUE(verbose)) {
    mode_str <- if (is_persam) {
      if (isTRUE(interactive)) "PerSAM (Interactive 1-click Exemplar)" else "PerSAM (One-shot Exemplar Batch)"
    } else if (is_grounded_sam) {
      paste0("Grounded-SAM (Text Prompt: '", if (is.character(prompt)) paste(prompt, collapse = ", ") else class_names[1], "')")
    } else if (is_yolo) {
      if (identical(task, "segment") && !grepl("[-_]seg", tolower(model[1]))) {
        paste0("YOLO Detection + Meta SAM 2.1 Refinement (", basename(model[1]), ")")
      } else {
        paste0("YOLO Pseudo-Labeling (", basename(model[1]), ")")
      }
    } else if (is_salient) {
      paste0("Salient Object Segmentation (", model[1], ")")
    } else {
      "Classical Image Analysis (analyze_objects)"
    }
    split_counts_str <- paste0(
      sum(split_assign == "train"), " train / ",
      sum(split_assign == "val"), " val",
      if (any(split_assign == "test")) paste0(" / ", sum(split_assign == "test"), " test") else ""
    )
    cli::cli_h2("YOLO Dataset Generation & Auto-Annotation")
    cli::cli_ul(c(
      "Annotation engine: {.val {mode_str}}",
      "Task: {.val {toupper(task)}}",
      "Classes: {.val {class_names}}",
      "Images to process: {num_imgs} ({split_counts_str})",
      "Target dataset directory: {.path {dir}}"
    ))
  }
  accumulated_prototypes_list <- vector("list", length(class_names))
  names(accumulated_prototypes_list) <- class_names
  cached_pts_list <- vector("list", length(class_names))
  names(cached_pts_list) <- class_names
  if (!is.null(prompt)) {
    if (is.list(prompt) && !is.data.frame(prompt)) {
      for (cl_i in seq_along(class_names)) {
        c_name <- class_names[cl_i]
        if (c_name %in% names(prompt)) {
          cached_pts_list[[cl_i]] <- prompt[[c_name]]
        } else if (length(prompt) >= cl_i) {
          cached_pts_list[[cl_i]] <- prompt[[cl_i]]
        }
      }
    } else if (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt)) {
      cached_pts_list[[1]] <- prompt
    }
  }
  stop_early <- FALSE
  show_pb <- isTRUE(verbose) && !(is_persam && isTRUE(interactive))
  pb <- if (show_pb) {
    cli::cli_progress_bar(
      format = "{cli::pb_spin} Annotating images [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      total = num_imgs,
      clear = FALSE
    )
  } else NULL
  for (i in seq_len(num_imgs)) {
    curr_item <- img_list[[i]]
    curr_img <- if (is.character(curr_item)) image_import(curr_item) else curr_item
    base_name <- img_names[i]
    split_name <- split_assign[i]
    dest_img_dir <- file.path(dir, "images", split_name)
    dest_lbl_dir <- file.path(dir, "labels", split_name)
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
      img_boxes_list <- list()
      img_polys_list <- list()
      for (cl_idx in seq_along(class_names)) {
        curr_cl_name <- class_names[cl_idx]
        curr_cl_id <- class_id[cl_idx]
        cached_pts <- cached_pts_list[[cl_idx]]
        if (isTRUE(interactive)) {
          accepted <- FALSE
          curr_sim_thresh <- sim_threshold
          curr_min_dist <- min_dist
          while (!accepted) {
            # Prompt for exemplars on image 1 OR if prompt_mode == "each" OR if no prototypes yet
            need_prompt <- (is.null(accumulated_prototypes_list[[cl_idx]]) && is.null(cached_pts)) ||
                           (identical(prompt_mode, "each") && !accepted)
            curr_pts <- NULL
            if (need_prompt) {
              if (isTRUE(verbose)) {
                if (is.null(accumulated_prototypes_list[[cl_idx]])) {
                  cli::cli_alert_info("Image {i}/{num_imgs} ({base_name}) - Class {cl_idx}/{length(class_names)} ({.val {curr_cl_name}}): Click on 2 to 5 exemplar instances in the plot window.")
                  cli::cli_alert_info("Press <Esc> or right-click when finished (or to skip this class if absent in this image).")
                } else {
                  cli::cli_alert_info("Image {i}/{num_imgs} ({base_name}) - Class {cl_idx}/{length(class_names)} ({.val {curr_cl_name}}): Click exemplar instances, or press <Esc> to autolabel using learned prototypes.")
                }
              }
              plot(curr_img)
              # Overlay any previously annotated classes on this image
              if (length(img_boxes_list) > 0L) {
                prev_bx_df <- do.call(rbind, img_boxes_list)
                prev_conts <- unlist(img_polys_list, recursive = FALSE)
                .render_overlay_annotations(
                  boxes = prev_bx_df,
                  polygons = prev_conts,
                  cw = img_w,
                  ch = img_h,
                  alpha = 0.25,
                  lwd = 1.5
                )
              }
              title_str <- sprintf("Image [%d/%d] (%s) - Class [%d/%d] ('%s'): Click exemplars (<Esc> if none/done)",
                                   i, num_imgs, base_name, cl_idx, length(class_names), curr_cl_name)
              graphics::title(main = title_str, col.main = "#0066cc", font.main = 2, cex.main = 1.05)
              pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.3), error = function(e) NULL)
              if (is.null(pts) || length(pts$x) == 0L) {
                if (is.null(accumulated_prototypes_list[[cl_idx]]) && is.null(cached_pts)) {
                  cli::cli_alert_warning("No point clicked for class {.val {curr_cl_name}} on image {.val {base_name}} (<Esc> pressed).")
                  cli::cli_bullets(c(
                    "*" = paste0(cli::col_green("[s] / skip"), ": Skip class {.val {curr_cl_name}} for this image (default)"),
                    "*" = paste0(cli::col_green("[r] / retry"), ": Retry clicking exemplar objects for this class"),
                    "*" = paste0(cli::col_green("[q] / quit"), ": Abort dataset export")
                  ))
                  ans_esc <- tryCatch(
                    tolower(trimws(readline("Select action [s = skip / r = retry / q = quit, default: s]: "))),
                    interrupt = function(e) "s"
                  )
                  if (ans_esc %in% c("q", "quit", "abort")) {
                    cli::cli_alert_info("Dataset export cancelled by user.")
                    stop_early <- TRUE
                    break
                  } else if (ans_esc %in% c("r", "retry")) {
                    next
                  } else {
                    accepted <- TRUE
                    break
                  }
                } else {
                  if (identical(prompt_mode, "each")) {
                    cli::cli_alert_info("No point clicked for class {.val {curr_cl_name}} (<Esc> pressed).")
                    cli::cli_bullets(c(
                      "*" = paste0(cli::col_green("[a] / autolabel"), ": Autolabel from visual memory (default)"),
                      "*" = paste0(cli::col_green("[s] / skip"), ": Skip class {.val {curr_cl_name}} for this image"),
                      "*" = paste0(cli::col_green("[r] / retry"), ": Retry clicking exemplar objects")
                    ))
                    ans_esc <- tryCatch(
                      tolower(trimws(readline("Select action [a = autolabel / s = skip / r = retry, default: a]: "))),
                      interrupt = function(e) "a"
                    )
                    if (ans_esc %in% c("s", "skip")) {
                      accepted <- TRUE
                      break
                    } else if (ans_esc %in% c("r", "retry")) {
                      next
                    } else {
                      curr_pts <- NULL
                    }
                  } else {
                    curr_pts <- NULL
                  }
                }
              } else {
                graphics::points(pts$x, pts$y, col = "#000000", pch = 1, cex = 1.3, lwd = 1.5)
                curr_pts <- cbind(x = pts$x, y = pts$y)
              }
            } else if (!is.null(cached_pts) && i == 1L) {
              curr_pts <- cached_pts
            } else {
              curr_pts <- NULL
            }
            if (isTRUE(stop_early)) break
            if (isTRUE(verbose)) {
              if (!is.null(accumulated_prototypes_list[[cl_idx]]) && is.null(curr_pts)) {
                n_proto <- length(accumulated_prototypes_list[[cl_idx]]$semantic)
                cli::cli_alert_info("Autolabeling {.val {curr_cl_name}} using {n_proto} accumulated prototype{?s}...")
              }
              cli::cli_progress_step(
                msg = "Running SAM 2.1 inference on {.val {base_name}} [{curr_cl_name}] [{toupper(engine)}]... ",
                msg_done = "SAM 2.1 inference complete for {.val {curr_cl_name}}"
              )
            }
            c_num_arr <- image_data(curr_img, type = "numeric")
            p_res <- tryCatch({
              .run_persam(
                mat = c_num_arr,
                exemplar_points = curr_pts,
                precomputed_prototypes = accumulated_prototypes_list[[cl_idx]],
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
            }, error = function(e) {
              cli::cli_alert_warning("Exemplar matching error on image {i} for class '{curr_cl_name}': {e$message}")
              NULL
            })
            if (isTRUE(verbose)) {
              cli::cli_progress_done()
            }
            r_boxes <- data.frame()
            r_polys <- list()
            if (!is.null(p_res) && !is.null(p_res$boxes) && nrow(p_res$boxes) > 0L) {
              r_boxes <- data.frame(
                xmin = as.numeric(p_res$boxes$xmin),
                ymin = as.numeric(p_res$boxes$ymin),
                xmax = as.numeric(p_res$boxes$xmax),
                ymax = as.numeric(p_res$boxes$ymax),
                conf = round(as.numeric(p_res$boxes$score), 3),
                class_id = curr_cl_id,
                class_name = curr_cl_name,
                label = curr_cl_name,
                stringsAsFactors = FALSE
              )
              # Filter edge artifacts
              if (isTRUE(filter_edge) && nrow(r_boxes) > 0L) {
                bw_px <- r_boxes$xmax - r_boxes$xmin
                bh_px <- r_boxes$ymax - r_boxes$ymin
                aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
                is_border <- (r_boxes$xmin <= 4 | r_boxes$xmax >= (img_w - 4) | r_boxes$ymin <= 4 | r_boxes$ymax >= (img_h - 4))
                keep <- !(is_border & (aspect > 3.5 | bw_px < 15 | bh_px < 15))
                r_boxes <- r_boxes[keep, , drop = FALSE]
                if (!is.null(p_res$contours) && length(p_res$contours) > 0L) {
                  p_res$contours <- p_res$contours[keep[seq_along(p_res$contours)]]
                }
              }
              if (!is.null(p_res$contours)) {
                for (ci in seq_along(p_res$contours)) {
                  cnt <- p_res$contours[[ci]]
                  if (is.matrix(cnt) && nrow(cnt) >= 3L) {
                    gx <- cnt[, 1] / img_w; gy <- cnt[, 2] / img_h
                    if (length(gx) > 32L) {
                      idx <- round(seq(1, length(gx), length.out = 32))
                      gx <- gx[idx]; gy <- gy[idx]
                    }
                    r_polys[[ci]] <- as.vector(rbind(gx, gy))
                  }
                }
              }
            }
            prev_bxs <- if (length(img_boxes_list) > 0L) do.call(rbind, img_boxes_list) else NULL
            prev_pls <- if (length(img_polys_list) > 0L) unlist(img_polys_list, recursive = FALSE) else NULL
            next_cl_name <- if (cl_idx < length(class_names)) class_names[cl_idx + 1L] else NULL
            # Interactive Review Step
            rev <- .interactive_fewshot_review_step(
              curr_im = curr_img,
              r_boxes = r_boxes,
              r_polys = r_polys,
              curr_points = curr_pts,
              c_num_arr = c_num_arr,
              cw = img_w,
              ch = img_h,
              pre_proto = accumulated_prototypes_list[[cl_idx]],
              latest_proto = if (!is.null(p_res) && !is.null(p_res$prototypes)) p_res$prototypes else NULL,
              task = task,
              label = curr_cl_name,
              class_id = curr_cl_id,
              conf_threshold = curr_sim_thresh,
              radius = NULL,
              min_dist = curr_min_dist,
              max_objects = max_objects,
              method = "sam",
              engine = engine,
              device_id = device_id,
              rainbow = rainbow,
              show_class = show_class,
              show_conf = show_conf,
              badge = badge,
              img_name = base_name,
              img_idx = i,
              total_imgs = num_imgs,
              class_idx = cl_idx,
              total_classes = length(class_names),
              next_class_name = next_cl_name,
              prior_boxes = prev_bxs,
              prior_polygons = prev_pls
            )
            if (identical(rev$action, "resample")) {
              accumulated_prototypes_list[[cl_idx]] <- NULL
              cached_pts_list[[cl_idx]] <- NULL
              curr_pts <- NULL
              next
            } else if (identical(rev$action, "stop")) {
              if (nrow(rev$boxes) > 0L) {
                img_boxes_list[[length(img_boxes_list) + 1L]] <- rev$boxes
                img_polys_list[[length(img_polys_list) + 1L]] <- rev$polygons
              }
              accepted <- TRUE
              stop_early <- TRUE
              break
            } else {
              if (nrow(rev$boxes) > 0L) {
                img_boxes_list[[length(img_boxes_list) + 1L]] <- rev$boxes
                img_polys_list[[length(img_polys_list) + 1L]] <- rev$polygons
              }
              accepted <- TRUE
              new_protos <- if (!is.null(rev$prototypes)) rev$prototypes else if (!is.null(p_res$prototypes)) p_res$prototypes else NULL
              if (!is.null(new_protos)) {
                new_protos <- .compact_prototypes(new_protos, max_protos = 8L)
                accumulated_prototypes_list[[cl_idx]] <- new_protos
                if (isTRUE(verbose)) {
                  cli::cli_alert_success("Visual prototype memory for {.val {curr_cl_name}}: {length(accumulated_prototypes_list[[cl_idx]]$semantic)} prototype{?s} active.")
                }
              }
            }
          } # while (!accepted)
          if (isTRUE(stop_early)) break
        } else {
          # Non-interactive / batch mode for class cl_idx
          c_num_arr <- image_data(curr_img, type = "numeric")
          p_res <- NULL
          if (i == 1L && is.null(accumulated_prototypes_list[[cl_idx]])) {
            if (is.null(cached_pts)) {
              cli::cli_abort("In non-interactive mode with {.val persam}, please provide {.arg prompt} for class {.val {curr_cl_name}}.")
            }
            p_res <- .run_persam(
              mat = c_num_arr,
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
            if (!is.null(p_res$prototypes)) accumulated_prototypes_list[[cl_idx]] <- p_res$prototypes
          } else {
            p_res <- .run_persam(
              mat = c_num_arr,
              exemplar_points = NULL,
              precomputed_prototypes = accumulated_prototypes_list[[cl_idx]],
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
          }
          if (!is.null(p_res) && !is.null(p_res$boxes) && nrow(p_res$boxes) > 0L) {
            cl_bx <- data.frame(
              xmin = as.numeric(p_res$boxes$xmin),
              ymin = as.numeric(p_res$boxes$ymin),
              xmax = as.numeric(p_res$boxes$xmax),
              ymax = as.numeric(p_res$boxes$ymax),
              conf = round(as.numeric(p_res$boxes$score), 3),
              class_id = curr_cl_id,
              class_name = curr_cl_name,
              label = curr_cl_name,
              stringsAsFactors = FALSE
            )
            cl_pl <- list()
            if (!is.null(p_res$contours)) {
              for (ci in seq_along(p_res$contours)) {
                cnt <- p_res$contours[[ci]]
                if (is.matrix(cnt) && nrow(cnt) >= 3L) {
                  gx <- cnt[, 1] / img_w; gy <- cnt[, 2] / img_h
                  if (length(gx) > 32L) {
                    idx <- round(seq(1, length(gx), length.out = 32))
                    gx <- gx[idx]; gy <- gy[idx]
                  }
                  cl_pl[[ci]] <- as.vector(rbind(gx, gy))
                }
              }
            }
            if (isTRUE(filter_edge) && nrow(cl_bx) > 0L) {
              bw_px <- cl_bx$xmax - cl_bx$xmin
              bh_px <- cl_bx$ymax - cl_bx$ymin
              aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
              is_border <- (cl_bx$xmin <= 4 | cl_bx$xmax >= (img_w - 4) | cl_bx$ymin <= 4 | cl_bx$ymax >= (img_h - 4))
              keep <- !(is_border & (aspect > 3.5 | bw_px < 15 | bh_px < 15))
              cl_bx <- cl_bx[keep, , drop = FALSE]
              if (length(cl_pl) > 0L) cl_pl <- cl_pl[keep[seq_along(cl_pl)]]
            }
            if (nrow(cl_bx) > 0L) {
              img_boxes_list[[length(img_boxes_list) + 1L]] <- cl_bx
              img_polys_list[[length(img_polys_list) + 1L]] <- cl_pl
            }
          }
        }
      } # for (cl_idx in seq_along(class_names))
      # Combine annotations across all classes for image i
      if (length(img_boxes_list) > 0L) {
        bx_df <- do.call(rbind, img_boxes_list)
        rownames(bx_df) <- NULL
      } else {
        bx_df <- data.frame()
      }
      conts_list <- unlist(img_polys_list, recursive = FALSE)
      if (is.null(conts_list)) conts_list <- list()
      if (isTRUE(stop_early)) break
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
    # MODEL 3: Existing YOLO Model (ONNX Assisted Pre-Labeling & Auto-SAM)
    # --------------------------------------------------------------------------
    } else if (is_yolo) {
      use_cond_sahi <- !isTRUE(sahi) && !is.null(sahi_if) && !isFALSE(sahi_if)
      can_slice <- (img_w > slice_width || img_h > slice_height)
      thresh_sahi <- if (isTRUE(sahi_if)) 0.01 else if (is.numeric(sahi_if)) as.numeric(sahi_if[1]) else 0.01
      if (is.numeric(thresh_sahi) && thresh_sahi > 1) thresh_sahi <- thresh_sahi / 100

      .eval_sahi_trigger <- function(df_boxes) {
        if (!use_cond_sahi || !can_slice) return(list(trigger = FALSE, reason = ""))
        if (is.null(df_boxes) || nrow(df_boxes) == 0L) {
          if (isTRUE(sahi_if_empty)) {
            return(list(trigger = TRUE, reason = "0 objects detected at full image scale"))
          }
          return(list(trigger = FALSE, reason = ""))
        }
        if (!("xmax" %in% names(df_boxes)) || !("xmin" %in% names(df_boxes))) {
          return(list(trigger = FALSE, reason = ""))
        }
        b_areas <- (as.numeric(df_boxes$xmax) - as.numeric(df_boxes$xmin)) * 
                   (as.numeric(df_boxes$ymax) - as.numeric(df_boxes$ymin))
        rel_areas <- b_areas / (as.numeric(img_w) * as.numeric(img_h))
        stat_val <- switch(sahi_metric,
                           "mean" = mean(rel_areas, na.rm = TRUE),
                           "median" = stats::median(rel_areas, na.rm = TRUE),
                           "min" = min(rel_areas, na.rm = TRUE),
                           mean(rel_areas, na.rm = TRUE))
        if (is.finite(stat_val) && stat_val <= thresh_sahi) {
          msg <- sprintf("relative %s box area (%.4f%%) <= threshold (%.4f%%)",
                         sahi_metric, stat_val * 100, thresh_sahi * 100)
          return(list(trigger = TRUE, reason = msg))
        }
        list(trigger = FALSE, reason = "")
      }

      if (identical(task, "detect")) {
        pad_val <- if (is.numeric(box_margin) && length(box_margin) >= 1L && !is.na(box_margin[1]) && box_margin[1] > 0) box_margin[1] else 0
        if (isTRUE(sahi)) {
          y_res <- image_detect_sahi(
            curr_img,
            model = model[1],
            slice_height = slice_height,
            slice_width = slice_width,
            overlap_height_ratio = overlap_height_ratio,
            overlap_width_ratio = overlap_width_ratio,
            conf_threshold = box_threshold,
            iou_threshold = iou_threshold,
            match_metric = match_metric,
            match_threshold = match_threshold,
            full_image_inference = full_image_inference,
            pad = pad_val,
            plot = FALSE,
            verbose = FALSE,
            engine = engine,
            device_id = device_id
          )
          bx_df <- as.data.frame(y_res)
        } else {
          y_res <- image_detect_dl(
            curr_img,
            model = model[1],
            conf_threshold = box_threshold,
            iou_threshold = iou_threshold,
            pad = pad_val,
            plot = FALSE,
            verbose = FALSE,
            engine = engine,
            device_id = device_id
          )
          bx_df <- as.data.frame(y_res)
          chk <- .eval_sahi_trigger(bx_df)
          if (isTRUE(chk$trigger)) {
            if (isTRUE(verbose)) {
              cli::cli_alert_info("Triggering SAHI inference for {.file {base_name}} ({chk$reason})...")
            }
            y_res <- image_detect_sahi(
              curr_img,
              model = model[1],
              slice_height = slice_height,
              slice_width = slice_width,
              overlap_height_ratio = overlap_height_ratio,
              overlap_width_ratio = overlap_width_ratio,
              conf_threshold = box_threshold,
              iou_threshold = iou_threshold,
              match_metric = match_metric,
              match_threshold = match_threshold,
              full_image_inference = full_image_inference,
              pad = pad_val,
              plot = FALSE,
              verbose = FALSE,
              engine = engine,
              device_id = device_id
            )
            bx_df <- as.data.frame(y_res)
          }
        }
      } else {
        # Task is segment!
        is_seg_model <- grepl("[-_]seg", tolower(model[1]))
        y_res <- NULL
        if (is_seg_model) {
          y_res <- tryCatch({
            if (isTRUE(sahi)) {
              image_segment_sahi(
                curr_img,
                model = model[1],
                slice_height = slice_height,
                slice_width = slice_width,
                overlap_height_ratio = overlap_height_ratio,
                overlap_width_ratio = overlap_width_ratio,
                conf_threshold = box_threshold,
                iou_threshold = iou_threshold,
                match_metric = match_metric,
                match_threshold = match_threshold,
                full_image_inference = full_image_inference,
                type = "highlight",
                plot = FALSE,
                verbose = FALSE,
                engine = engine,
                device_id = device_id
              )
            } else {
              init_res <- image_segment_dl(
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
              df_chk <- if (!is.null(init_res$results)) as.data.frame(init_res$results) else if (!is.null(init_res$boxes)) as.data.frame(init_res$boxes) else as.data.frame(init_res)
              chk <- .eval_sahi_trigger(df_chk)
              if (isTRUE(chk$trigger)) {
                if (isTRUE(verbose)) {
                  cli::cli_alert_info("Triggering SAHI segmentation for {.file {base_name}} ({chk$reason})...")
                }
                image_segment_sahi(
                  curr_img,
                  model = model[1],
                  slice_height = slice_height,
                  slice_width = slice_width,
                  overlap_height_ratio = overlap_height_ratio,
                  overlap_width_ratio = overlap_width_ratio,
                  conf_threshold = box_threshold,
                  iou_threshold = iou_threshold,
                  match_metric = match_metric,
                  match_threshold = match_threshold,
                  full_image_inference = full_image_inference,
                  type = "highlight",
                  plot = FALSE,
                  verbose = FALSE,
                  engine = engine,
                  device_id = device_id
                )
              } else {
                init_res
              }
            }
          }, error = function(e) NULL)
        }
        conts_cand <- if (!is.null(y_res)) {
          if (!is.null(y_res$contours) && length(y_res$contours) > 0L) {
            y_res$contours
          } else if (!is.null(attr(y_res, "contours")) && length(attr(y_res, "contours")) > 0L) {
            attr(y_res, "contours")
          } else list()
        } else list()

        if (length(conts_cand) > 0L) {
          bx_df <- if (!is.null(y_res$results)) as.data.frame(y_res$results) else if (!is.null(y_res$boxes)) as.data.frame(y_res$boxes) else as.data.frame(y_res)
          conts_list <- conts_cand
        } else {
          # Model is a detection model (or custom detection weights).
          # Automatically run object detection + Meta SAM 2.1 polygon refinement!
          if (i == 1L && isTRUE(verbose)) {
            bm_info <- if (is.numeric(box_margin) && box_margin[1] > 0) {
              paste0(" (margin: ", round(box_margin[1] * (if (box_margin[1] <= 1) 100 else 1)), "%)")
            } else ""
            sahi_info <- if (isTRUE(sahi)) " with SAHI" else ""
            cli::cli_alert_info("YOLO detection model detected with {.code task = 'segment'}. Automatically activating {.strong Meta SAM 2.1} to segment bounding boxes{sahi_info} into polygon masks{bm_info}!")
          }
          if (isTRUE(sahi)) {
            det_res <- image_detect_sahi(
              curr_img,
              model = model[1],
              slice_height = slice_height,
              slice_width = slice_width,
              overlap_height_ratio = overlap_height_ratio,
              overlap_width_ratio = overlap_width_ratio,
              conf_threshold = box_threshold,
              iou_threshold = iou_threshold,
              match_metric = match_metric,
              match_threshold = match_threshold,
              full_image_inference = full_image_inference,
              plot = FALSE,
              verbose = FALSE,
              engine = engine,
              device_id = device_id
            )
          } else {
            det_res <- image_detect_dl(
              curr_img,
              model = model[1],
              conf_threshold = box_threshold,
              iou_threshold = iou_threshold,
              plot = FALSE,
              verbose = FALSE,
              engine = engine,
              device_id = device_id
            )
            df_chk <- as.data.frame(det_res)
            chk <- .eval_sahi_trigger(df_chk)
            if (isTRUE(chk$trigger)) {
              if (isTRUE(verbose)) {
                cli::cli_alert_info("Triggering SAHI detection for {.file {base_name}} ({chk$reason})...")
              }
              det_res <- image_detect_sahi(
                curr_img,
                model = model[1],
                slice_height = slice_height,
                slice_width = slice_width,
                overlap_height_ratio = overlap_height_ratio,
                overlap_width_ratio = overlap_width_ratio,
                conf_threshold = box_threshold,
                iou_threshold = iou_threshold,
                match_metric = match_metric,
                match_threshold = match_threshold,
                full_image_inference = full_image_inference,
                plot = FALSE,
                verbose = FALSE,
                engine = engine,
                device_id = device_id
              )
            }
          }
          bx_df <- as.data.frame(det_res)
          if (!is.null(bx_df) && nrow(bx_df) > 0L) {
            sam_res <- .sam_segment_boxes(
              curr_img,
              boxes = bx_df,
              threshold = 0.5,
              n_points = n_points,
              min_area = min_area,
              box_margin = box_margin,
              dilate = dilate_val,
              erode = erode_val,
              engine = engine,
              device_id = device_id,
              verbose = FALSE
            )
            conts_list <- sam_res$contours
          } else {
            conts_list <- list()
          }
        }
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
        res_obj <- if (length(extra_args) > 0L) {
          do.call(analyze_objects, c(list(img = curr_img, plot = FALSE), extra_args))
        } else {
          analyze_objects(curr_img, plot = FALSE)
        }
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
      } else if (is.data.frame(annotations) && all(c("xmin", "ymin", "xmax", "ymax") %in% names(annotations))) {
        bx_df <- annotations
        if (identical(task, "segment") && nrow(bx_df) > 0L) {
          sam_res <- .sam_segment_boxes(curr_img, boxes = bx_df, threshold = 0.5, n_points = n_points, min_area = min_area, box_margin = box_margin, dilate = dilate_val, erode = erode_val, engine = engine, device_id = device_id, verbose = FALSE)
          conts_list <- sam_res$contours
        }
      } else if (is.list(annotations) && !is.null(annotations$boxes) && is.data.frame(annotations$boxes)) {
        bx_df <- annotations$boxes
        if (identical(task, "segment") && nrow(bx_df) > 0L) {
          sam_res <- .sam_segment_boxes(curr_img, boxes = bx_df, threshold = 0.5, n_points = n_points, min_area = min_area, box_margin = box_margin, dilate = dilate_val, erode = erode_val, engine = engine, device_id = device_id, verbose = FALSE)
          conts_list <- sam_res$contours
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
    # OPTIONAL SIZE FILTERING
    # --------------------------------------------------------------------------
    if ((!is.null(lower_size) || !is.null(upper_size)) && !is.null(bx_df) && nrow(bx_df) > 0L) {
      size_m_use <- if (!is.null(size_mode)) size_mode[1] else "relative"
      f_res <- .filter_by_size(
        boxes = bx_df,
        contours = conts_list,
        lower_size = lower_size,
        upper_size = upper_size,
        size_mode = size_m_use
      )
      bx_df <- f_res$boxes
      if (!is.null(conts_list) && length(conts_list) > 0L) {
        conts_list <- f_res$contours
      }
    }

    # --------------------------------------------------------------------------
    # OPTIONAL CROP TO OBJECT BOUNDING BOXES
    # --------------------------------------------------------------------------
    if (isTRUE(crop_to_objects)) {
      cmarg <- if (is.numeric(crop_margin) && length(crop_margin) >= 1L && !is.na(crop_margin[1])) as.numeric(crop_margin[1]) else 0.05
      c_res <- .yolo_crop_to_objects(curr_img, boxes = bx_df, polygons = conts_list, margin = cmarg)
      if (isTRUE(c_res$cropped)) {
        curr_img <- c_res$img
        bx_df <- c_res$boxes
        conts_list <- c_res$polygons
        img_w <- dim(curr_img)[1]
        img_h <- dim(curr_img)[2]
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
          cid <- if (!is.null(bx_df$class_id) && length(bx_df$class_id) >= k && !is.na(bx_df$class_id[k])) {
            as.integer(bx_df$class_id[k])
          } else if (!is.null(bx_df$class_name) && length(bx_df$class_name) >= k && !is.na(bx_df$class_name[k])) {
            lbl <- tolower(trimws(bx_df$class_name[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
          } else if (!is.null(bx_df$label) && length(bx_df$label) >= k && !is.na(bx_df$label[k])) {
            lbl <- tolower(trimws(bx_df$label[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
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
          if (is.null(c_mat)) next
          if (is.vector(c_mat) && length(c_mat) >= 6L && !is.matrix(c_mat)) {
            cid <- if (!is.null(bx_df$class_id) && length(bx_df$class_id) >= k && !is.na(bx_df$class_id[k])) {
              as.integer(bx_df$class_id[k])
            } else if (!is.null(bx_df$class_name) && length(bx_df$class_name) >= k && !is.na(bx_df$class_name[k])) {
              lbl <- tolower(trimws(bx_df$class_name[k]))
              m_idx <- match(lbl, tolower(trimws(class_names)))
              if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
            } else if (!is.null(bx_df$label) && length(bx_df$label) >= k && !is.na(bx_df$label[k])) {
              lbl <- tolower(trimws(bx_df$label[k]))
              m_idx <- match(lbl, tolower(trimws(class_names)))
              if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
            } else if (length(class_id) == length(conts_list)) {
              class_id[k]
            } else {
              class_id[1]
            }
            line_entry <- paste(c(cid, sprintf("%.6f", c_mat)), collapse = " ")
            txt_lines <- c(txt_lines, line_entry)
            next
          }
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
          cid <- if (!is.null(bx_df$class_id) && length(bx_df$class_id) >= k && !is.na(bx_df$class_id[k])) {
            as.integer(bx_df$class_id[k])
          } else if (!is.null(bx_df$class_name) && length(bx_df$class_name) >= k && !is.na(bx_df$class_name[k])) {
            lbl <- tolower(trimws(bx_df$class_name[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
          } else if (!is.null(bx_df$label) && length(bx_df$label) >= k && !is.na(bx_df$label[k])) {
            lbl <- tolower(trimws(bx_df$label[k]))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
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
        # Fallback: try SAM 2.1 or convert 4 box corners to polygon
        if (identical(task, "segment") && length(conts_list) == 0L) {
          sam_res <- tryCatch(.sam_segment_boxes(curr_img, boxes = bx_df, threshold = 0.5, n_points = n_points, min_area = min_area, box_margin = box_margin, dilate = dilate_val, erode = erode_val, engine = engine, device_id = device_id, verbose = FALSE), error = function(e) NULL)
          if (!is.null(sam_res) && length(sam_res$contours) > 0L) {
            conts_list <- sam_res$contours
          }
        }
        for (k in seq_len(nrow(bx_df))) {
          bx <- bx_df[k, ]
          if (((bx$xmax - bx$xmin) * (bx$ymax - bx$ymin)) < min_area) next
          x1 <- pmin(pmax(bx$xmin / img_w, 0.0), 1.0); x2 <- pmin(pmax(bx$xmax / img_w, 0.0), 1.0)
          y1 <- pmin(pmax(bx$ymin / img_h, 0.0), 1.0); y2 <- pmin(pmax(bx$ymax / img_h, 0.0), 1.0)
          cid <- if (!is.null(bx$class_id) && !is.na(bx$class_id)) {
            as.integer(bx$class_id)
          } else if (!is.null(bx$class_name) && !is.na(bx$class_name)) {
            lbl <- tolower(trimws(bx$class_name))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
          } else if (!is.null(bx$label) && !is.na(bx$label)) {
            lbl <- tolower(trimws(bx$label))
            m_idx <- match(lbl, tolower(trimws(class_names)))
            if (!is.na(m_idx)) class_id[m_idx] else class_id[1]
          } else {
            class_id[1]
          }
          poly_pts <- c(x1, y1, x2, y1, x2, y2, x1, y2)
          line_entry <- paste(c(cid, round(poly_pts, 5L)), collapse = " ")
          txt_lines <- c(txt_lines, line_entry)
        }
      }
    }
    orig_src <- if (length(img_src_files) >= i) img_src_files[i] else NA_character_
    logged_orig_paths <- c(logged_orig_paths, orig_src)
    logged_dest_imgs <- c(logged_dest_imgs, dest_img_path)
    logged_dest_lbls <- c(logged_dest_lbls, dest_lbl_path)
    # Save or move image to split folder
    if (isTRUE(move) && !is.na(orig_src) && file.exists(orig_src)) {
      if (isTRUE(crop_to_objects)) {
        # Cropped image is a derived sub-region: never delete the raw original image
        image_export(curr_img, name = paste0(base_name, ".jpg"), subfolder = dest_img_dir)
      } else {
        if (normalizePath(orig_src, winslash = "/", mustWork = FALSE) != normalizePath(dest_img_path, winslash = "/", mustWork = FALSE)) {
          mv_ok <- file.rename(orig_src, dest_img_path)
          if (!mv_ok) {
            image_export(curr_img, name = paste0(base_name, ".jpg"), subfolder = dest_img_dir)
            unlink(orig_src, force = TRUE)
          }
        }
      }
    } else {
      if (is.na(orig_src) || normalizePath(orig_src, winslash = "/", mustWork = FALSE) != normalizePath(dest_img_path, winslash = "/", mustWork = FALSE)) {
        if (!is.na(orig_src) && !isTRUE(crop_to_objects) && tolower(tools::file_ext(orig_src)) %in% c("jpg", "jpeg") && file.exists(orig_src)) {
          file.copy(orig_src, dest_img_path, overwrite = TRUE)
        } else {
          image_export(curr_img, name = paste0(base_name, ".jpg"), subfolder = dest_img_dir)
        }
      }
    }
    if (split_name == "train") {
      train_images_saved <- c(train_images_saved, dest_img_path)
    } else if (split_name == "val") {
      val_images_saved <- c(val_images_saved, dest_img_path)
    } else {
      test_images_saved <- c(test_images_saved, dest_img_path)
    }
    # Write text labels
    writeLines(txt_lines, dest_lbl_path)
    total_objects_annotated <- total_objects_annotated + length(txt_lines)
    if (!is.null(pb)) {
      cli::cli_progress_update(id = pb)
    } else if (isTRUE(verbose)) {
      cli::cli_alert_info("Annotated {length(txt_lines)} object(s) in {.file {base_name}.jpg} [{split_name}].")
    }
    # Optional review image rendering
    cur_rev <- NA_character_
    if (do_review && (i %in% review_indices)) {
      tryCatch({
        rev_out <- file.path(review_dir, paste0("review_", base_name, ".jpg"))
        out_w <- 1200L
        out_h <- max(100L, as.integer(round(1200 * img_h / img_w)))
        dev_opened <- FALSE
        tryCatch({
          grDevices::jpeg(rev_out, width = out_w, height = out_h, quality = 92, res = 120)
          dev_opened <- TRUE
          cur_rev <- rev_out
        }, error = function(e) {
          tryCatch({
            alt_out <- file.path(review_dir, paste0("review_", i, "_", base_name, ".jpg"))
            grDevices::jpeg(alt_out, width = out_w, height = out_h, quality = 92, res = 120)
            dev_opened <<- TRUE
            cur_rev <<- alt_out
          }, error = function(e2) {
            tryCatch({
              grDevices::png(rev_out, width = out_w, height = out_h, res = 120)
              dev_opened <<- TRUE
              cur_rev <<- rev_out
            }, error = function(e3) NULL)
          })
        })
        if (dev_opened) {
          plot(curr_img)
          n_objs <- length(txt_lines)
          if (isTRUE(rainbow) && n_objs > 0L) {
            rev_border <- grDevices::rainbow(n_objs, s = 0.85, v = 0.95)
            rev_poly <- grDevices::adjustcolor(rev_border, alpha.f = 0.40)
          } else {
            rev_border <- rep("#00FFCC", max(1L, n_objs))
            rev_poly <- rep("#00FFCC55", max(1L, n_objs))
          }
          if (n_objs > 0L) {
            parsed_objs <- lapply(txt_lines, function(ln) {
              vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
              if (length(vals) < 5L) return(NULL)
              cls_id <- as.integer(vals[1])
              coords <- vals[-1]
              if (length(coords) == 4L) {
                cx <- coords[1] * img_w
                cy <- coords[2] * img_h
                bw <- coords[3] * img_w
                bh <- coords[4] * img_h
                list(
                  type = "box",
                  cls_id = cls_id,
                  xmin = max(0, cx - bw / 2),
                  xmax = min(img_w, cx + bw / 2),
                  ymin = max(0, cy - bh / 2),
                  ymax = min(img_h, cy + bh / 2)
                )
              } else if (length(coords) >= 6L && length(coords) %% 2 == 0) {
                px <- coords[seq(1, length(coords), by = 2L)] * img_w
                py <- coords[seq(2, length(coords), by = 2L)] * img_h
                list(
                  type = "poly",
                  cls_id = cls_id,
                  px = px,
                  py = py,
                  xmin = max(0, min(px)),
                  xmax = min(img_w, max(px)),
                  ymin = max(0, min(py)),
                  ymax = min(img_h, max(py))
                )
              } else {
                NULL
              }
            })
            parsed_objs <- parsed_objs[!sapply(parsed_objs, is.null)]
            # Draw contours / boxes / points
            for (r in seq_along(parsed_objs)) {
              obj <- parsed_objs[[r]]
              col_b <- rev_border[((r - 1L) %% length(rev_border)) + 1L]
              col_p <- rev_poly[((r - 1L) %% length(rev_poly)) + 1L]
              if (identical(obj$type, "poly")) {
                if (isTRUE(show_polygon)) {
                  graphics::polygon(obj$px, obj$py, col = col_p, border = col_b, lwd = 1.5)
                }
                if (isTRUE(show_bbox)) {
                  graphics::rect(obj$xmin, obj$ymin, obj$xmax, obj$ymax, border = col_b, lwd = 1.5, lty = if (isTRUE(show_polygon)) 2 else 1)
                }
                if (isTRUE(show_points)) {
                  graphics::points(mean(obj$px), mean(obj$py), col = col_b, bg = col_b, pch = 19, cex = 1.2)
                }
              } else {
                if (isTRUE(show_bbox)) {
                  graphics::rect(obj$xmin, obj$ymin, obj$xmax, obj$ymax, border = col_b, lwd = 2)
                }
                if (isTRUE(show_points)) {
                  cx <- (obj$xmin + obj$xmax) / 2
                  cy <- (obj$ymin + obj$ymax) / 2
                  graphics::points(cx, cy, col = col_b, bg = col_b, pch = 19, cex = 1.2)
                }
              }
            }
            # Draw labels / badges
            if ((isTRUE(show_labels) || isTRUE(show_conf)) && (isTRUE(show_bbox) || isTRUE(show_polygon) || isTRUE(show_points))) {
              for (r in seq_along(parsed_objs)) {
                obj <- parsed_objs[[r]]
                col_b <- rev_border[((r - 1L) %% length(rev_border)) + 1L]
                lbl_parts <- character(0)
                if (isTRUE(show_labels)) {
                  m_idx <- match(obj$cls_id, class_id)
                  cur_lbl <- if (!is.na(m_idx) && m_idx <= length(class_names)) class_names[m_idx] else if (length(class_names) > obj$cls_id) class_names[obj$cls_id + 1L] else as.character(obj$cls_id)
                  if (nzchar(cur_lbl)) lbl_parts <- c(lbl_parts, cur_lbl)
                }
                if (isTRUE(show_conf)) {
                  if (!is.null(bx_df) && "conf" %in% names(bx_df) && length(bx_df$conf) >= r && !is.na(bx_df$conf[r])) {
                    lbl_parts <- c(lbl_parts, sprintf("%.2f", bx_df$conf[r]))
                  }
                }
                tag <- paste(lbl_parts, collapse = " ")
                if (nzchar(tag)) {
                  tw <- graphics::strwidth(tag, cex = 0.8)
                  th <- graphics::strheight(tag, cex = 0.8)
                  y_top <- max(th + 6, obj$ymin)
                  if (isTRUE(badge)) {
                    bg_tag <- grDevices::adjustcolor(col_b, alpha.f = 0.90)
                    graphics::rect(obj$xmin, y_top - th - 6, obj$xmin + tw + 6, y_top, col = bg_tag, border = NA)
                    graphics::text(obj$xmin + 3, y_top - 3, labels = tag, col = "#FFFFFF", cex = 0.8, font = 2, adj = c(0, 1))
                  } else {
                    graphics::text(obj$xmin + 3, y_top - 3, labels = tag, col = col_b, cex = 0.8, font = 2, adj = c(0, 1))
                  }
                }
              }
            }
          }
          graphics::title(
            sprintf("%s: %d %s", paste0(base_name, ".jpg"), length(txt_lines), if (identical(task, "segment")) "segmented" else "detected"),
            col.main = if (isTRUE(rainbow)) "#222222" else "#00FFCC",
            cex.main = 1.0
          )
          grDevices::dev.off()
        }
      }, error = function(e) {
        if (grDevices::dev.cur() > 1L) {
          tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        }
      })
    }
    logged_rev_imgs <- c(logged_rev_imgs, cur_rev)
    if (is.character(curr_item)) {
      rm(curr_img)
    }
  }
  if (!is.null(pb)) {
    cli::cli_progress_done(id = pb)
  }
  # Ensure val split has at least 1 image (required by YOLO validator)
  if (length(val_images_saved) == 0L && length(train_images_saved) > 1L) {
    src_img <- train_images_saved[length(train_images_saved)]
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
    "val: images/val"
  )
  if (length(test_images_saved) > 0L) {
    yaml_content <- c(yaml_content, "test: images/test")
  }
  yaml_content <- c(yaml_content, "", "names:")
  for (idx in seq_along(class_names)) {
    cid <- if (length(class_id) >= idx) class_id[idx] else (idx - 1L)
    yaml_content <- c(yaml_content, paste0("  ", cid, ": ", class_names[idx]))
  }
  writeLines(yaml_content, yaml_path)
  # Save undo history manifest (.yolo_history.rds)
  hist_file <- file.path(dir, ".yolo_history.rds")
  existing_hist <- if (file.exists(hist_file)) tryCatch(readRDS(hist_file), error = function(e) list()) else list()
  if (!is.list(existing_hist)) existing_hist <- list()
  history_entry <- list(
    id = paste0("export_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    timestamp = Sys.time(),
    operation = "export",
    input_dir = input_orig_dir,
    output_dir = dir,
    moved = isTRUE(move),
    images = data.frame(
      orig = logged_orig_paths,
      dest = logged_dest_imgs,
      lbl = logged_dest_lbls,
      rev = logged_rev_imgs,
      stringsAsFactors = FALSE
    ),
    yaml_file = yaml_path,
    review_dir = if (do_review) review_dir else NULL
  )
  existing_hist[[length(existing_hist) + 1L]] <- history_entry
  saveRDS(existing_hist, hist_file)
  if (isTRUE(verbose)) {
    tot_saved <- length(train_images_saved) + length(val_images_saved) + length(test_images_saved)
    cli::cli_alert_success(
      "YOLO dataset successfully generated in {.path {dir}} ({total_objects_annotated} objects across {tot_saved} images)."
    )
    if (do_review) {
      cli::cli_alert_info("Review previews written to: {.path {review_dir}}")
    }
  }
  out <- list(
    dir = dir,
    yaml_file = yaml_path,
    train_images = train_images_saved,
    val_images = val_images_saved,
    test_images = test_images_saved,
    review_dir = if (do_review) review_dir else NULL,
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
#' @param show_bbox Logical. Whether to draw bounding boxes and annotations. Defaults to `TRUE`. Supported alias: `bbox`.
#' @param show_polygon Logical. For segmentation datasets, whether to draw polygons. Defaults to `TRUE`. Supported alias: `polygon`.
#' @param show_points Logical. If `TRUE`, draws points at the centroid of each bounding box or polygon using the same color as the box. Defaults to `FALSE`. Supported alias: `points`.
#' @param show_labels Logical. If `TRUE` (default), displays class name labels. Supported aliases: `labels`, `show_class`.
#' @param bbox,polygon,points,labels Backward-compatible aliases for `show_bbox`, `show_polygon`, `show_points`, and `show_labels`.
#' @param n Integer. Number of sample images to preview. Defaults to `4L`. Supported aliases: `n_images`, `num_images`.
#' @param which Optional integer vector of image indices (e.g. `1:4`) or character vector of
#'   image names or stems (e.g. `c("img_3_aug1", "img_553_aug2")`) to plot specific images.
#' @param images Alias for `which`.
#' @param only_augmented Logical. If `TRUE`, only previews augmented images (containing `_aug`).
#'   Defaults to `TRUE` if `dir` is an augmented dataset (`yolo_dataset_aug`), and `FALSE` otherwise.
#' @param only_original Logical. If `TRUE`, only previews original pristine images (excluding synthetic
#'   augmentations `_aug[0-9]+`, synthetic backgrounds `synth_bg_`, external backgrounds `bg_`, and scrambled backgrounds `scramble_`).
#'   Defaults to `FALSE` (or `NULL`).
#' @param pattern Optional regex pattern to filter specific image names (e.g. `"toras5\\.png"`). Defaults to `NULL`.
#' @param random Logical. If `TRUE` (default when `which` is `NULL`), samples random images from the split. If `FALSE`, picks the first `n` images.
#' @param seed Optional integer seed for reproducible random selection.
#' @param mfrow Numeric vector of length 2 giving `c(rows, cols)` for the multi-panel layout.
#' @param rainbow Logical. If `TRUE` (default), assigns distinct rainbow colors to each object within the image (matching review previews). If `FALSE`, colors are assigned by class.
#' @param col Border color for bounding boxes/polygons. If provided, overrides `rainbow` and default class palette.
#' @param lwd Numeric. Line width for borders. Defaults to `1.5`.
#' @param fill Fill color for polygons with transparency (e.g. `"#00FF6633"`). If `NULL`, no fill.
#' @param label_bg Logical. If `TRUE`, draws a subtle background badge behind each label. If `FALSE` (default), draws clean, discreet text matching review previews without blocking objects.
#' @param cex Numeric. Text size for class name labels. Defaults to `0.6`.
#' @param title Logical or character. If `TRUE` (default), displays panel titles with filename and object count.
#' @param save Logical. If `TRUE`, saves the generated plot to a file. Defaults to `FALSE`.
#' @param file_name Optional character string. Custom file name or path to save the plot (e.g. `"dataset_preview.png"`).
#' @param width Numeric. Width of the output plot. If `<= 50`, interpreted as inches; if `> 50`, interpreted as pixels. Defaults to `max(7, nc * 3.5)`.
#' @param height Numeric. Height of the output plot. If `<= 50`, interpreted as inches; if `> 50`, interpreted as pixels. Defaults to `max(6, nr * 3.5)`.
#' @param dpi Numeric. Resolution in dots per inch (DPI) for raster formats (PNG, JPEG, TIFF). Defaults to `300`.
#' @param ... Additional arguments.
#'
#' @return An invisible list of the sampled image file paths.
#' @export
#' @examplesIf interactive()
#' dataset <- yolo_dataset_export(img, dir = "dataset_cafe", model = "persam")
#' yolo_dataset_preview(dataset, n = 4, random = TRUE)
#' plot(dataset, n = 6)
#' plot(dataset, n = 4, save = TRUE, file_name = "preview.png")
yolo_dataset_preview <- function(dir = "yolo_dataset",
                                 split = c("train", "val", "test"),
                                 show_bbox = TRUE,
                                 show_polygon = TRUE,
                                 show_points = FALSE,
                                 show_labels = TRUE,
                                 bbox = show_bbox,
                                 polygon = show_polygon,
                                 points = show_points,
                                 labels = show_labels,
                                 n = 4L,
                                 which = NULL,
                                 images = NULL,
                                 only_augmented = NULL,
                                 only_original = NULL,
                                 pattern = NULL,
                                 random = is.null(which) && is.null(images),
                                 seed = NULL,
                                 mfrow = NULL,
                                 rainbow = TRUE,
                                 col = NULL,
                                 lwd = 1.5,
                                 fill = NULL,
                                 label_bg = FALSE,
                                 cex = 0.6,
                                 title = TRUE,
                                 save = FALSE,
                                 file_name = NULL,
                                 width = NULL,
                                 height = NULL,
                                 dpi = 300,
                                 ...) {
  dots <- list(...)
  if (!missing(bbox)) show_bbox <- bbox
  if (!missing(polygon)) show_polygon <- polygon
  if (!missing(points)) show_points <- points
  if (!missing(labels)) show_labels <- labels
  if (!is.null(dots$show_bbox)) show_bbox <- dots$show_bbox
  if (!is.null(dots$bbox)) show_bbox <- dots$bbox
  if (!is.null(dots$show_boxes)) show_bbox <- dots$show_boxes
  if (!is.null(dots$show_polygon)) show_polygon <- dots$show_polygon
  if (!is.null(dots$polygon)) show_polygon <- dots$polygon
  if (!is.null(dots$show_points)) show_points <- dots$show_points
  if (!is.null(dots$points)) show_points <- dots$points
  if (!is.null(dots$show_point)) show_points <- dots$show_point
  if (!is.null(dots$point)) show_points <- dots$point
  if (!is.null(dots$show_labels)) show_labels <- dots$show_labels
  if (!is.null(dots$labels)) show_labels <- dots$labels
  if (!is.null(dots$show_class)) show_labels <- dots$show_class
  if (!is.null(dots$n_images)) n <- dots$n_images
  if (!is.null(dots$num_images)) n <- dots$num_images
  if (!is.null(dots$n_samples)) n <- dots$n_samples
  if (!is.null(dots$img) && is.null(which) && is.null(images)) which <- dots$img
  if (!is.null(images) && is.null(which)) which <- images
  if (!is.null(dots$only_augmented) && is.null(only_augmented)) only_augmented <- dots$only_augmented
  if (!is.null(dots$only_original) && is.null(only_original)) only_original <- dots$only_original
  point_pch <- if (!is.null(dots$point_pch)) dots$point_pch else if (!is.null(dots$pch)) dots$pch else 19
  point_cex <- if (!is.null(dots$point_cex)) dots$point_cex else if (!is.null(dots$cex_point)) dots$cex_point else if (!is.null(dots$point_size)) dots$point_size else 1.2
  split <- match.arg(split)
  if (is.null(dir)) {
    cli::cli_abort("The {.arg dir} argument cannot be NULL. Provide a directory path or a {.code yolo_dataset} object.")
  }
  dir_obj <- dir
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
    if (any(grepl("(kpt_shape|task:\\s*pose)", y_lines, ignore.case = TRUE))) {
      return(yolo_dataset_pose_preview(
        dir = dir,
        split = split,
        show_bbox = show_bbox,
        n = n,
        which = which,
        images = images,
        pattern = pattern,
        random = random,
        seed = seed,
        mfrow = mfrow,
        rainbow = rainbow,
        col = col,
        line_lwd = lwd,
        title = title,
        save = save,
        file_name = file_name,
        width = width,
        height = height,
        dpi = dpi,
        ...
      ))
    }
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
    cand_lbl1 <- sub("([/\\\\])images([/\\\\]|$)", "\\1labels\\2", img_dir)
    cand_lbl2 <- sub("([/\\\\])images([/\\\\])", "\\1labels\\2", img_dir)
    for (c_lbl in c(cand_lbl1, cand_lbl2)) {
      if (c_lbl != img_dir && dir.exists(c_lbl)) {
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
  img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", ignore.case = TRUE, full.names = TRUE)
  if (length(img_files) == 0L) {
    cli::cli_abort("No images found in {.path {img_dir}}.")
  }
  is_aug_context <- inherits(dir_obj, "yolo_dataset_aug") ||
                    (is.list(dir_obj) && isTRUE(dir_obj$is_augmented)) ||
                    (is.character(raw_dir) && grepl("_aug", raw_dir[1]))
  if (isTRUE(only_original)) {
    if (isTRUE(only_augmented)) {
      cli::cli_abort("Cannot specify both {.arg only_original = TRUE} and {.arg only_augmented = TRUE}.")
    }
    only_augmented <- FALSE
  } else if (is.null(only_augmented)) {
    only_augmented <- is_aug_context && is.null(which) && is.null(pattern)
  }
  if (isTRUE(only_augmented)) {
    aug_matches <- grepl("_aug[0-9]+", basename(img_files))
    if (any(aug_matches)) {
      img_files <- img_files[aug_matches]
    }
  } else if (isTRUE(only_original)) {
    # Exclude offline augmented variants, synthetic/external backgrounds, and scrambles
    synth_patt <- "(_aug[0-9]+|^synth_bg_|^bg_|^scramble_)"
    orig_matches <- !grepl(synth_patt, basename(img_files), ignore.case = TRUE)
    if (any(orig_matches)) {
      img_files <- img_files[orig_matches]
    } else {
      cli::cli_abort("No original images found in split {.val {split}} of {.path {img_dir}}.")
    }
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
  norm_num <- function(x) {
    x <- sub("^(.*?[a-zA-Z_-])0+([0-9]+)$", "\\1\\2", x)
    x <- sub("^0+([0-9]+)$", "\\1", x)
    x
  }
  if (!is.null(which)) {
    if (is.numeric(which)) {
      w_idx <- as.integer(which)
      valid_idx <- w_idx[w_idx >= 1L & w_idx <= length(img_files)]
      if (length(valid_idx) == 0L) {
        cli::cli_abort("Specified index {.val {which}} is out of bounds. Available images in split {.val {split}}: 1 to {length(img_files)}.")
      }
      sample_files <- img_files[valid_idx]
    } else if (is.character(which)) {
      selected_files <- character(0)
      for (w in which) {
        w_raw <- trimws(as.character(w))
        w_clean <- sub("[./\\]+$", "", w_raw)
        if (!nzchar(w_clean)) w_clean <- w_raw
        target_base <- basename(w_clean)
        target_stem <- tools::file_path_sans_ext(target_base)

        # 1. Exact match on basename in current split
        m <- which(tolower(basename(img_files)) == tolower(target_base))

        # 2. Exact match on stem in current split
        if (length(m) == 0L) {
          m <- which(tolower(all_stems) == tolower(target_stem))
        }

        # 3. Numeric zero-padding equivalence (e.g. img_3 matches img_003 or vice-versa)
        if (length(m) == 0L) {
          m <- which(tolower(norm_num(all_stems)) == tolower(norm_num(target_stem)))
        }

        # 4. Glob wildcard ONLY if * or ? explicitly provided
        if (length(m) == 0L && grepl("[*?]", w_raw)) {
          rx_base <- utils::glob2rx(target_base)
          m <- grep(rx_base, basename(img_files), ignore.case = TRUE)
          if (length(m) == 0L) {
            rx_stem <- utils::glob2rx(target_stem)
            m <- grep(rx_stem, all_stems, ignore.case = TRUE)
          }
        }

        if (length(m) > 0L) {
          selected_files <- c(selected_files, img_files[m])
          next
        }

        # 5. Search in other splits if not found in current split
        all_dataset_imgs <- list.files(dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
        all_dataset_imgs <- all_dataset_imgs[!grepl("[/\\\\]review[/\\\\]|^review_", basename(all_dataset_imgs))]
        other_files <- all_dataset_imgs[!normalizePath(all_dataset_imgs, winslash = "/", mustWork = FALSE) %in% normalizePath(img_files, winslash = "/", mustWork = FALSE)]

        if (length(other_files) > 0L) {
          other_stems <- tools::file_path_sans_ext(basename(other_files))
          om <- which(tolower(basename(other_files)) == tolower(target_base))
          if (length(om) == 0L) {
            om <- which(tolower(other_stems) == tolower(target_stem))
          }
          if (length(om) == 0L) {
            om <- which(tolower(norm_num(other_stems)) == tolower(norm_num(target_stem)))
          }
          if (length(om) == 0L && grepl("[*?]", w_raw)) {
            rx_base <- utils::glob2rx(target_base)
            om <- grep(rx_base, basename(other_files), ignore.case = TRUE)
            if (length(om) == 0L) {
              rx_stem <- utils::glob2rx(target_stem)
              om <- grep(rx_stem, other_stems, ignore.case = TRUE)
            }
          }
          if (length(om) > 0L) {
            found_paths <- other_files[om]
            for (fp in found_paths) {
              found_split <- if (grepl("[/\\\\]train[/\\\\]?", fp)) "train" else if (grepl("[/\\\\](val|valid)[/\\\\]?", fp)) "val" else if (grepl("[/\\\\]test[/\\\\]?", fp)) "test" else "another split"
              cli::cli_alert_info("Image {.file {basename(fp)}} was found in split {.val {found_split}} instead of {.val {split}}.")
            }
            selected_files <- c(selected_files, found_paths)
            next
          }
        }

        cli::cli_abort("No image matching {.val {w}} found in dataset {.path {dir}} (searched split {.val {split}} and all other splits).")
      }
      sample_files <- unique(selected_files)
    } else {
      cli::cli_abort("{.arg which} must be an integer index vector or character vector of image names.")
    }
    n_samples <- length(sample_files)
  } else {
    n_samples <- min(as.integer(n), length(img_files))
    if (isTRUE(random)) {
      if (!is.null(seed)) set.seed(seed)
      sample_files <- sample(img_files, size = n_samples)
    } else {
      sample_files <- img_files[seq_len(n_samples)]
    }
  }
  if (n_samples > 1L) {
    if (is.null(mfrow)) {
      if (n_samples > 16L) {
        cli::cli_alert_warning("Showing the first 16 of {n_samples} matched images to fit graphics device margins. Specify {.arg mfrow} to display more.")
        sample_files <- sample_files[seq_len(16L)]
        n_samples <- 16L
      }
      nc <- ceiling(sqrt(n_samples))
      nr <- ceiling(n_samples / nc)
      mfrow <- c(nr, nc)
    }
  } else {
    mfrow <- c(1L, 1L)
  }

  ds_tag <- if (is.character(dir)) basename(dir) else "yolo_dataset"
  default_file <- paste0(ds_tag, "_preview.png")
  def_w <- max(7, mfrow[2] * 3.5)
  def_h <- max(6, mfrow[1] * 3.5)

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

  if (n_samples > 1L) {
    op <- graphics::par(mfrow = mfrow, mar = c(1, 1, if (isTRUE(title)) 2.2 else 1, 1))
    on.exit(graphics::par(op), add = TRUE)
  } else {
    op <- graphics::par(mar = c(1, 1, if (isTRUE(title)) 2.2 else 1, 1))
    on.exit(graphics::par(op), add = TRUE)
  }
  for (img_path in sample_files) {
    base_name <- tools::file_path_sans_ext(basename(img_path))
    alt_lbl <- sub("([/\\\\])images([/\\\\])", "\\1labels\\2", img_path)
    alt_lbl <- paste0(tools::file_path_sans_ext(alt_lbl), ".txt")
    if (file.exists(alt_lbl)) {
      lbl_path <- alt_lbl
    } else {
      lbl_path  <- file.path(lbl_dir, paste0(base_name, ".txt"))
      if (!file.exists(lbl_path)) {
        for (cand_parent in c(dir, dirname(img_dir), dirname(dirname(img_dir)))) {
          alt_lbl2 <- file.path(cand_parent, "labels", split, paste0(base_name, ".txt"))
          if (file.exists(alt_lbl2)) {
            lbl_path <- alt_lbl2
            break
          }
          alt_lbl3 <- file.path(cand_parent, "labels", paste0(base_name, ".txt"))
          if (file.exists(alt_lbl3)) {
            lbl_path <- alt_lbl3
            break
          }
        }
      }
    }
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
          if (isTRUE(show_bbox)) {
            graphics::rect(xmin, ymin, xmax, ymax, border = draw_col, lwd = lwd, xpd = TRUE)
          }
          if (isTRUE(show_points)) {
            pt_col <- if (!is.null(dots$col_point)) dots$col_point else if (!is.null(dots$point_col)) dots$point_col else draw_col
            graphics::points(cx, cy, col = pt_col, bg = pt_col, pch = point_pch, cex = point_cex, xpd = TRUE)
          }
          if (isTRUE(show_labels) && (isTRUE(show_bbox) || isTRUE(show_points))) {
            if (isTRUE(label_bg)) {
              th <- graphics::strheight(lbl_text, cex = cex) * 1.2
              tw <- graphics::strwidth(lbl_text, cex = cex) * 1.2
              lbl_y1 <- max(0, ymin - th)
              graphics::rect(xmin, lbl_y1, xmin + tw, ymin, col = grDevices::adjustcolor(draw_col, alpha.f = 0.6), border = NA, xpd = TRUE)
              graphics::text(xmin + tw / 2, ymin - th / 2, labels = lbl_text, col = "black", cex = cex, font = 2, xpd = TRUE)
            } else {
              graphics::text(xmin, ymin - 3, labels = lbl_text, col = draw_col, cex = cex, pos = 4, font = 2, xpd = TRUE)
            }
          }
        } else if (length(coords) >= 6L && length(coords) %% 2 == 0) {
          # Segmentation polygon format: x1 y1 x2 y2 ...
          xs <- coords[seq(1, length(coords), by = 2L)] * img_w
          ys <- coords[seq(2, length(coords), by = 2L)] * img_h
          cx <- mean(xs)
          cy <- mean(ys)
          if (isTRUE(show_polygon)) {
            poly_fill <- if (!is.null(fill)) fill else if (isTRUE(rainbow)) grDevices::adjustcolor(draw_col, alpha.f = 0.35) else NULL
            graphics::polygon(xs, ys, border = draw_col, col = poly_fill, lwd = lwd, xpd = TRUE)
          }
          if (isTRUE(show_bbox)) {
            min_x <- min(xs); max_x <- max(xs)
            min_y <- min(ys); max_y <- max(ys)
            graphics::rect(min_x, min_y, max_x, max_y, border = draw_col, lwd = lwd, lty = if (isTRUE(show_polygon)) 2 else 1, xpd = TRUE)
          }
          if (isTRUE(show_points)) {
            pt_col <- if (!is.null(dots$col_point)) dots$col_point else if (!is.null(dots$point_col)) dots$point_col else draw_col
            graphics::points(cx, cy, col = pt_col, bg = pt_col, pch = point_pch, cex = point_cex, xpd = TRUE)
          }
          if (isTRUE(show_labels) && (isTRUE(show_bbox) || isTRUE(show_polygon) || isTRUE(show_points))) {
            min_x <- min(xs); min_y <- min(ys)
            if (isTRUE(label_bg)) {
              th <- graphics::strheight(lbl_text, cex = cex) * 1.2
              tw <- graphics::strwidth(lbl_text, cex = cex) * 1.2
              lbl_y1 <- max(0, min_y - th)
              graphics::rect(min_x, lbl_y1, min_x + tw, min_y, col = grDevices::adjustcolor(draw_col, alpha.f = 0.6), border = NA, xpd = TRUE)
              graphics::text(min_x + tw / 2, min_y - th / 2, labels = lbl_text, col = "black", cex = cex, font = 2, xpd = TRUE)
            } else {
              graphics::text(min_x, min_y - 3, labels = lbl_text, col = draw_col, cex = cex, pos = 4, font = 2, xpd = TRUE)
            }
          }
        }
      }
    }
  }
  invisible(sample_files)
}


#' @title Interactively or Programmatically Review YOLO Dataset Annotations
#' @name yolo_dataset_review
#' @aliases yolo_review
#' @description
#' Revies annotations for a specific image in a YOLO dataset. Supports interactive inspection,
#' removing false positive objects (by clicking on boxes/masks, center dots, or by index),
#' adding missed objects, re-labeling from scratch or appending detections using deep learning models
#' (e.g., PerSAM one-shot visual exemplars, YOLO ONNX/PT models, Grounded-SAM text prompts, or classical
#' segmentation via [analyze_objects()]), and synchronizing updated labels with visual review previews.
#'
#' @details
#' `yolo_dataset_review()` follows the exact same engine and argument conventions as [yolo_dataset_export()].
#' It can be used in two complementary modes:
#' \itemize{
#'   \item \strong{Interactive GUI}: Launches an interactive graphics session where the user can click on
#'     unwanted objects to delete them (`[r]`), click opposite corners to add new bounding boxes (`[a]`),
#'     delete objects by index directly in console (`[d <index>]`), or re-run exemplar matching with PerSAM.
#'   \item \strong{Programmatic Scripting}: Specify `remove = c(2, 5)` to drop objects, `keep = c(1, 3)` to retain
#'     specific instances, or pass `model = "persam"` / `model = "best.onnx"` to re-annotate automatically.
#' }
#' When changes are confirmed, the original `.txt` file is safely backed up as `<stem>.txt.bak` (if `backup = TRUE`),
#' coordinates are written in normalized YOLO format, and the review image in `review/review_<stem>.jpg` is updated.
#'
#' @param dir Character or `yolo_dataset` object. Path to the YOLO dataset root directory (default: `"yolo_dataset"`).
#' @param image Character vector, integer vector, or `NULL`. The image(s) to edit. Can be:
#'   \itemize{
#'     \item A character vector of image paths, filenames, or stems (e.g., `imgs` from [list.files()], `"img_001.jpg"`, or `c("ds1_001.jpg", "ds2_002.jpg")`).
#'     \item A numeric vector of 1-based indices (e.g., `image = 1` or `image = 1:5`).
#'     \item `"all"`: Reviews all images in the dataset sequentially.
#'     \item `NULL` (default): lists available images and prompts the user in interactive sessions.
#'   }
#' @param model Character or `NULL`. Model to use for re-labeling or object extraction:
#'   \itemize{
#'     \item `"persam"` or `"sam-exemplar"`: One-shot interactive exemplar-guided segmentation/detection using SAM 2.1.
#'     \item `"grounded-sam"`: Open-vocabulary text prompt detection using Grounding DINO + SAM.
#'     \item An ONNX or PT model file path (e.g., `"yolo11n.onnx"`, `"best.pt"`, or standard YOLO model name).
#'     \item `NULL` (default): edits existing annotations without running a new model (or runs classical [analyze_objects()] if segmentation arguments are passed).
#'   }
#' @param task Character. The YOLO task: `"detect"` (bounding boxes) or `"segment"` (polygons).
#'   Auto-detected from dataset if available.
#' @param mode Character. How model predictions interact with existing labels:
#'   \itemize{
#'     \item `"replace"` (default when `model` is provided): Replaces existing annotations with new predictions.
#'     \item `"append"`: Keeps existing annotations and appends new model predictions.
#'     \item `"edit"` (default when `model` is `NULL`): Starts from existing annotations and allows interactive or programmatic additions/deletions.
#'   }
#' @param class_names Character vector. Class names for the dataset. If `NULL` (default),
#'   automatically inferred from `data.yaml` or existing annotations.
#' @param class_id Integer vector. Numeric class IDs corresponding to `class_names`.
#' @param prompt Character, list, or matrix. Prompt for model inference (text prompt for Grounded-SAM,
#'   or exemplar coordinates for PerSAM).
#' @param interactive Logical. Whether to activate the interactive graphical editor (using `graphics::locator()`)
#'   to click-to-delete unwanted objects and click-to-add missed objects.
#'   Defaults to `base::interactive()`.
#' @param remove Integer vector or `"all"`. Specific object indices to remove programmatically (e.g., `remove = c(2, 5)`).
#' @param keep Integer vector. Specific object indices to retain, removing all others.
#' @param delete_image Logical. If `TRUE`, permanently deletes the specified image(s), matching `.txt` label file(s),
#'   and review preview(s) from the dataset. Can also be triggered interactively by typing `rm`. Defaults to `FALSE`. Supported alias: `rm`.
#' @param box_threshold Numeric. Confidence threshold for bounding box detection (default: `0.25`).
#' @param text_threshold Numeric. Text threshold for Grounded-SAM (default: `0.25`).
#' @param sim_threshold Numeric. Cosine similarity threshold for PerSAM exemplar matching (default: `0.50`).
#' @param min_dist Numeric. Minimum pixel distance between object centers for PerSAM (default: `20`).
#' @param iou_threshold Numeric. NMS IoU threshold for overlapping detections (default: `0.50`).
#' @param min_area Numeric. Minimum pixel area to keep an object (default: `15`).
#' @param filter_edge Logical. Whether to automatically filter out border artifacts touching image edges (default: `FALSE`).
#' @param backup Logical. Whether to create a backup of the original `.txt` label file as `<stem>.txt.bak` before saving (default: `TRUE`).
#' @param review Logical. Whether to update the visual review image in `review/review_<stem>.jpg` if present (default: `TRUE`).
#' @param rainbow Logical. If `TRUE` (default), uses distinct rainbow palette colors for each object.
#' @param show_bbox Logical. Whether to display bounding box rectangles (default: `TRUE`).
#' @param show_polygon Logical. Whether to display segmentation polygons (default: `TRUE`).
#' @param show_points Logical. Whether to display solid colored dots at object centers (default: `FALSE`).
#' @param show_labels Logical. Whether to display class names and confidence badges (default: `TRUE`).
#' @param show_indices Logical. Whether to display object index numbers (`#1`, `#2`, ...) for easy identification (default: `TRUE`).
#' @param point_pch Integer. Plotting symbol for center points (default: `19`, solid circle).
#' @param point_cex Numeric. Size of center points (default: `1.6`).
#' @param engine Character. Execution engine: `"gpu"` (default) or `"cpu"`.
#' @param device_id Integer. GPU device ID (default: `-1` for auto).
#' @param verbose Logical. Whether to display formatted progress and messages via `cli` (default: `TRUE`).
#' @param ... Additional arguments passed to [analyze_objects()] or model inference.
#'
#' @return An invisible list containing:
#'   \item{image_path}{Full path to the image file}
#'   \item{label_path}{Full path to the updated `.txt` label file}
#'   \item{boxes}{`data.frame` of updated bounding boxes and coordinates}
#'   \item{polygons}{`list` of updated segmentation polygons}
#'   \item{n_objects}{Number of saved objects}
#'   \item{backup_path}{Path to backup file (if created)}
#'
#' @export
#' @examples
#' \dontrun{
#' # 1. Interactive editing: click to remove broken/unwanted seeds, press [c] to save
#' yolo_dataset_review("meu_dataset", image = "img_005.jpg")
#'
#' # 2. Relabel image with PerSAM (click 1 exemplar object to re-segment all instances)
#' yolo_dataset_review("meu_dataset", image = "img_005.jpg", model = "persam")
#'
#' # 3. Relabel using a custom trained ONNX model
#' yolo_dataset_review("meu_dataset", image = "img_005.jpg", model = "yolo11n.onnx")
#'
#' # 4. Programmatic removal of objects #2 and #5 without opening GUI
#' yolo_dataset_review("meu_dataset", image = "img_005.jpg", remove = c(2, 5))
#' }
yolo_dataset_review <- function(dir = "yolo_dataset",
                              image = NULL,
                              model = NULL,
                              task = c("detect", "segment"),
                              mode = c("edit", "replace", "append"),
                              class_names = NULL,
                              class_id = NULL,
                              prompt = NULL,
                              interactive = base::interactive(),
                              remove = NULL,
                              keep = NULL,
                              delete_image = FALSE,
                              box_threshold = 0.25,
                              text_threshold = 0.25,
                              sim_threshold = 0.50,
                              min_dist = 20,
                              iou_threshold = 0.50,
                              min_area = 15,
                              filter_edge = FALSE,
                              backup = TRUE,
                              review = TRUE,
                              rainbow = TRUE,
                              show_bbox = TRUE,
                              show_polygon = TRUE,
                              show_points = FALSE,
                              show_labels = TRUE,
                              show_indices = TRUE,
                              point_pch = 19,
                              point_cex = 1.6,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              verbose = TRUE,
                              ...) {
  dots <- list(...)
  if (!is.null(dots$rm) && missing(delete_image)) delete_image <- dots$rm
  if (!is.null(dots$delete_image) && missing(delete_image)) delete_image <- dots$delete_image
  engine <- match.arg(engine)

  # Resolve dataset root
  dir_path <- if (inherits(dir, "yolo_dataset") || (is.list(dir) && !is.null(dir$dir))) {
    dir$dir
  } else if (is.character(dir) && length(dir) == 1L) {
    dir
  } else {
    cli::cli_abort("{.arg dir} must be a dataset root directory or a {.code yolo_dataset} object.")
  }
  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")
  }

  # Auto-detect task if not explicitly supplied
  if (missing(task)) {
    detected_task <- tryCatch(.yolo_detect_task(dir_path), error = function(e) "detect")
    task <- if (identical(detected_task, "segment")) "segment" else "detect"
  } else {
    task <- match.arg(task)
  }

  # Mode defaults
  if (missing(mode)) {
    mode <- if (!is.null(model)) "replace" else "edit"
  } else {
    mode <- match.arg(mode)
  }

  # Infer classes from data.yaml if NULL
  if (is.null(class_names)) {
    for (cand_root in c(dir_path, dirname(dir_path), dirname(dirname(dir_path)))) {
      cl <- tryCatch(.yolo_get_classes(cand_root), error = function(e) NULL)
      if (!is.null(cl) && length(cl) > 0L) {
        class_names <- cl
        break
      }
    }
    if (is.null(class_names) || length(class_names) == 0L) {
      class_names <- "object"
    }
  }
  if (is.null(class_id)) {
    class_id <- seq_along(class_names) - 1L
  }

  # Locate all image files in dataset
  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  img_files <- list.files(dir_path, pattern = img_pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  img_files <- img_files[!grepl("[/\\\\]review[/\\\\]|^review_", basename(img_files))]

  if (length(img_files) == 0L) {
    cli::cli_abort("No image files found in dataset {.path {dir_path}}.")
  }

  # Natural alphanumeric sorting (e.g. img_1, img_2, ..., img_10, img_100)
  num_part <- suppressWarnings(as.numeric(gsub(".*?([0-9]+)[^0-9]*$", "\\1", basename(img_files))))
  prefix_part <- gsub("[0-9]+[^0-9]*$", "", basename(img_files))
  img_files <- img_files[order(prefix_part, ifelse(is.na(num_part), 0, num_part), img_files)]

  all_stems <- tools::file_path_sans_ext(basename(img_files))
  all_rel <- gsub(paste0("^", gsub("([/\\\\])", "[/\\\\]", dir_path), "[/\\\\]?"), "", img_files)

  norm_num <- function(x) {
    x <- sub("^(.*?[a-zA-Z_-])0+([0-9]+)$", "\\1\\2", x)
    x <- sub("^0+([0-9]+)$", "\\1", x)
    x
  }

  # Resolve target images / indices
  if (is.null(image) || (is.character(image) && length(image) == 1L && identical(tolower(image), "all"))) {
    is_sequential_mode <- TRUE
    subset_indices <- seq_along(img_files)
    if (isTRUE(interactive)) {
      cli::cli_h3("Available images in dataset ({length(img_files)} images):")
      preview_list <- head(all_rel, 10L)
      for (k in seq_along(preview_list)) {
        cli::cli_bullets(c(" " = paste0("[", k, "] {.file ", preview_list[k], "}")))
      }
      if (length(img_files) > 10L) {
        cli::cli_alert_info("... and {length(img_files) - 10L} more images.")
      }
      ans_img <- trimws(readline("Select starting image (enter number or filename, default: 1, 'q' to cancel): "))
      if (tolower(ans_img) %in% c("q", "quit", "abort", "cancel")) {
        cli::cli_alert_info("Review session cancelled.")
        return(invisible(NULL))
      }
      if (nzchar(ans_img)) {
        start_idx <- 1L
        if (grepl("^\\d+$", ans_img)) {
          si <- as.integer(ans_img)
          if (si >= 1L && si <= length(img_files)) start_idx <- si
        } else {
          clean_ans <- sub("[./\\]+$", "", ans_img)
          target_base <- basename(clean_ans)
          target_stem <- tools::file_path_sans_ext(target_base)
          m_idx <- which(tolower(all_stems) == tolower(target_stem) |
                           tolower(basename(img_files)) == tolower(target_base) |
                           tolower(all_rel) == tolower(clean_ans))
          if (length(m_idx) == 0L) {
            m_idx <- which(tolower(norm_num(all_stems)) == tolower(norm_num(target_stem)))
          }
          if (length(m_idx) > 0L) start_idx <- m_idx[1]
        }
        if (start_idx > 1L) {
          subset_indices <- subset_indices[start_idx:length(subset_indices)]
        }
      }
    }
  } else if (is.numeric(image)) {
    if (any(image < 1L | image > length(img_files))) {
      bad <- image[image < 1L | image > length(img_files)]
      cli::cli_abort("Image index out of range [1, {length(img_files)}]: {.val {bad}}.")
    }
    subset_indices <- as.integer(image)
    is_sequential_mode <- (length(subset_indices) > 1L)
  } else if (is.character(image)) {
    norm_img_files <- normalizePath(img_files, winslash = "/", mustWork = FALSE)
    matched_indices <- integer(0)
    unmatched <- character(0)

    for (im_item in image) {
      if (!nzchar(trimws(im_item))) next
      item_clean <- trimws(im_item)

      # 1. Match directory path on disk
      if (dir.exists(item_clean)) {
        dir_norm <- normalizePath(item_clean, winslash = "/", mustWork = FALSE)
        dir_match <- which(startsWith(norm_img_files, paste0(dir_norm, "/")) | norm_img_files == dir_norm)
        if (length(dir_match) > 0L) {
          matched_indices <- c(matched_indices, dir_match)
          next
        }
      }

      # 2. Match dataset split/subfolder (e.g., "train", "val", "test", "valid")
      clean_sub <- tolower(sub("[./\\]+$", "", item_clean))
      rel_segments <- strsplit(all_rel, "[/\\\\]")
      split_match <- which(vapply(rel_segments, function(parts) {
        any(tolower(parts) == clean_sub)
      }, logical(1)))
      if (length(split_match) > 0L) {
        matched_indices <- c(matched_indices, split_match)
        next
      }

      idx <- integer(0)
      if (file.exists(item_clean)) {
        target_norm <- normalizePath(item_clean, winslash = "/", mustWork = FALSE)
        idx <- which(norm_img_files == target_norm)
      }
      if (length(idx) == 0L) {
        clean_img <- sub("[./\\]+$", "", item_clean)
        target_base <- basename(clean_img)
        target_stem <- tools::file_path_sans_ext(target_base)
        idx <- which(tolower(all_stems) == tolower(target_stem) |
                       tolower(basename(img_files)) == tolower(target_base) |
                       tolower(all_rel) == tolower(clean_img))
      }
      if (length(idx) == 0L) {
        idx <- which(tolower(norm_num(all_stems)) == tolower(norm_num(target_stem)))
      }
      if (length(idx) > 0L) {
        matched_indices <- c(matched_indices, idx[1])
      } else {
        unmatched <- c(unmatched, im_item)
      }
    }

    if (length(matched_indices) == 0L) {
      cli::cli_abort("None of the specified images could be found in dataset {.path {dir_path}}.")
    }
    if (length(unmatched) > 0L && isTRUE(verbose)) {
      cli::cli_alert_warning("{length(unmatched)} of {length(image)} image(s) not found in dataset: {.val {head(unmatched, 3)}}")
    }

    subset_indices <- unique(matched_indices)
    is_sequential_mode <- (length(subset_indices) > 1L)
  } else {
    cli::cli_abort("Invalid {.arg image} parameter. Must be NULL, 'all', numeric indices, or a character vector of image paths/stems.")
  }

  step_i <- 1L
  total_subset <- length(subset_indices)
  total_imgs <- length(img_files)
  reviewed_results <- list()
  session_saved <- 0L
  active_class_idx <- 1L

  repeat {
    if (step_i < 1L || step_i > total_subset) {
      if (isTRUE(verbose)) {
        if (total_subset < total_imgs) {
          cli::cli_alert_success("Finished reviewing all {total_subset} selected image{?s}! ({session_saved} image{?s} updated).")
        } else {
          cli::cli_alert_success("Finished reviewing all {total_imgs} images! ({session_saved} image{?s} updated).")
        }
      }
      break
    }

    curr_idx <- subset_indices[step_i]
    target_img_path <- normalizePath(img_files[curr_idx], winslash = "/", mustWork = FALSE)
    base_stem <- tools::file_path_sans_ext(basename(target_img_path))
    img_dir <- dirname(target_img_path)

    # Determine split
    split_name <- if (grepl("[/\\\\]train[/\\\\]?", img_dir)) {
      "train"
    } else if (grepl("[/\\\\]val[/\\\\]?", img_dir)) {
      "val"
    } else if (grepl("[/\\\\]test[/\\\\]?", img_dir)) {
      "test"
    } else {
      "images"
    }

    # Resolve target label path
    cand_lbl_same <- file.path(img_dir, paste0(base_stem, ".txt"))
    cand_lbl_yolo <- file.path(sub("((/|\\\\))images((/|\\\\)|$)", "\\1labels\\3", img_dir), paste0(base_stem, ".txt"))
    dest_lbl_path <- if (file.exists(cand_lbl_same)) {
      cand_lbl_same
    } else if (file.exists(cand_lbl_yolo)) {
      cand_lbl_yolo
    } else if (grepl("((/|\\\\))images((/|\\\\)|$)", img_dir)) {
      dir.create(dirname(cand_lbl_yolo), recursive = TRUE, showWarnings = FALSE)
      cand_lbl_yolo
    } else {
      cand_lbl_same
    }

    # Handle programmatic image deletion (delete_image = TRUE or rm = TRUE)
    if (isTRUE(delete_image)) {
      unlink(target_img_path, force = TRUE)
      unlink(dest_lbl_path, force = TRUE)
      unlink(paste0(dest_lbl_path, ".bak"), force = TRUE)
      rev_cand <- file.path(dir_path, "review", paste0("review_", base_stem, ".jpg"))
      unlink(rev_cand, force = TRUE)
      unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
      unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
      unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
      unlink(file.path(dir_path, "train", "labels.cache"), force = TRUE)
      unlink(file.path(dir_path, "val", "labels.cache"), force = TRUE)
      unlink(file.path(dir_path, "test", "labels.cache"), force = TRUE)

      if (isTRUE(verbose)) {
        cli::cli_alert_success("Permanently removed {.file {basename(target_img_path)}} and its label from dataset.")
      }

      reviewed_results[[length(reviewed_results) + 1L]] <- list(
        image_path = target_img_path,
        label_path = dest_lbl_path,
        boxes = data.frame(),
        polygons = list(),
        n_objects = 0L,
        backup_path = NULL,
        status = "deleted"
      )

      step_i <- step_i + 1L
      next
    }

    # Load image
    curr_img <- image_import(target_img_path)
    dims <- dim(curr_img)
    img_w <- dims[1]
    img_h <- dims[2]

    # Parse existing labels
    existing_boxes <- data.frame(
      index = integer(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      cx = numeric(0),
      cy = numeric(0),
      class_id = integer(0),
      class_name = character(0),
      conf = numeric(0),
      stringsAsFactors = FALSE
    )
    existing_polys <- list()

    if (file.exists(dest_lbl_path)) {
      raw_lines <- readLines(dest_lbl_path, warn = FALSE)
      raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
      if (length(raw_lines) > 0L) {
        b_rows <- list()
        p_rows <- list()
        for (li in seq_along(raw_lines)) {
          vals <- as.numeric(strsplit(trimws(raw_lines[li]), "\\s+")[[1]])
          if (length(vals) < 5L) next
          cid <- as.integer(vals[1])
          coords <- vals[-1]
          cname <- if (length(class_names) > cid) class_names[cid + 1L] else as.character(cid)
          if (length(coords) == 4L) {
            cx <- coords[1] * img_w
            cy <- coords[2] * img_h
            bw <- coords[3] * img_w
            bh <- coords[4] * img_h
            b_rows[[length(b_rows) + 1L]] <- data.frame(
              index = length(b_rows) + 1L,
              xmin = cx - bw / 2,
              ymin = cy - bh / 2,
              xmax = cx + bw / 2,
              ymax = cy + bh / 2,
              cx = cx,
              cy = cy,
              class_id = cid,
              class_name = cname,
              conf = 1.0,
              stringsAsFactors = FALSE
            )
            p_rows[[length(p_rows) + 1L]] <- NULL
          } else {
            px <- coords[c(TRUE, FALSE)] * img_w
            py <- coords[c(FALSE, TRUE)] * img_h
            b_rows[[length(b_rows) + 1L]] <- data.frame(
              index = length(b_rows) + 1L,
              xmin = min(px),
              ymin = min(py),
              xmax = max(px),
              ymax = max(py),
              cx = mean(px),
              cy = mean(py),
              class_id = cid,
              class_name = cname,
              conf = 1.0,
              stringsAsFactors = FALSE
            )
            p_rows[[length(p_rows) + 1L]] <- coords
          }
        }
        if (length(b_rows) > 0L) existing_boxes <- do.call(rbind, b_rows)
        existing_polys <- p_rows
      }
    }

    boxes <- existing_boxes
    polygons <- existing_polys

    # ----------------------------------------------------------------------------
    # MODEL RE-ANNOTATION STEP (if model is supplied)
    # ----------------------------------------------------------------------------
    if (!is.null(model)) {
      m_lower <- tolower(model[1])
      is_persam <- m_lower %in% c("persam", "per-sam", "sam-exemplar")
      is_grounded_sam <- m_lower %in% c("grounded-sam", "grounded_sam")
      is_yolo <- !is_persam && !is_grounded_sam && nzchar(m_lower) && (
        grepl("yolo", m_lower) || grepl("\\.(onnx|pt)$", model[1]) || file.exists(model[1])
      )
      is_salient <- m_lower %in% c("ben2", "u2netp", "u2net", "birefnet-lite", "rmbg-1.4", "rmbg-2.0", "withoutbg", "silueta")

      pred_boxes <- data.frame()
      pred_polys <- list()

      if (is_persam) {
        if (is.null(prompt) && isTRUE(interactive)) {
          plot(
            curr_img,
            main = sprintf("PerSAM: Click 1 to 5 exemplar object(s) on '%s' (<Esc> when done)", base_stem),
            col.main = "#0066cc",
            font.main = 2,
            cex.main = 1.05,
            mar = c(1.5, 1.5, 3.8, 1.5)
          )
          cli::cli_alert_info("Click on exemplar instances of {.val {class_names[1]}} in the plot window (<Esc> when finished).")
          pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.4), error = function(e) NULL)
          if (!is.null(pts) && length(pts$x) > 0L) {
            prompt <- cbind(x = pts$x, y = pts$y)
          }
        }
        c_num_arr <- image_data(curr_img, type = "numeric")
        if (isTRUE(verbose)) {
          cli::cli_progress_step(msg = "Running PerSAM / SAM 2.1 matching on {.file {base_stem}}...",
                                 msg_done = "PerSAM matching complete")
        }
        p_res <- tryCatch({
          .run_persam(
            mat = c_num_arr,
            exemplar_points = prompt,
            sim_threshold = sim_threshold,
            min_dist = min_dist,
            iou_threshold = iou_threshold,
            feat_res = 256L,
            engine = engine,
            device_id = device_id,
            verbose = FALSE,
            mask = (task == "segment")
          )
        }, error = function(e) {
          cli::cli_alert_warning("PerSAM inference error: {e$message}")
          NULL
        })
        if (isTRUE(verbose)) cli::cli_progress_done()

        if (!is.null(p_res) && !is.null(p_res$boxes) && nrow(p_res$boxes) > 0L) {
          pred_boxes <- data.frame(
            index = seq_len(nrow(p_res$boxes)),
            xmin = as.numeric(p_res$boxes$xmin),
            ymin = as.numeric(p_res$boxes$ymin),
            xmax = as.numeric(p_res$boxes$xmax),
            ymax = as.numeric(p_res$boxes$ymax),
            cx = (as.numeric(p_res$boxes$xmin) + as.numeric(p_res$boxes$xmax)) / 2,
            cy = (as.numeric(p_res$boxes$ymin) + as.numeric(p_res$boxes$ymax)) / 2,
            class_id = class_id[1],
            class_name = class_names[1],
            conf = round(as.numeric(p_res$boxes$score), 3),
            stringsAsFactors = FALSE
          )
          if (!is.null(p_res$contours)) {
            for (ci in seq_along(p_res$contours)) {
              cnt <- p_res$contours[[ci]]
              if (is.matrix(cnt) && nrow(cnt) >= 3L) {
                gx <- cnt[, 1] / img_w; gy <- cnt[, 2] / img_h
                if (length(gx) > 32L) {
                  s_idx <- round(seq(1, length(gx), length.out = 32))
                  gx <- gx[s_idx]; gy <- gy[s_idx]
                }
                pred_polys[[ci]] <- as.vector(rbind(gx, gy))
              } else {
                pred_polys[[ci]] <- NULL
              }
            }
          }
        }
      } else if (is_grounded_sam) {
        t_prompt <- if (is.character(prompt) && length(prompt) >= 1L) prompt[1] else class_names[1]
        if (isTRUE(verbose)) cli::cli_alert_info("Running Grounded-SAM with text prompt: {.val {t_prompt}}...")
        gs_res <- image_segment_dl(
          curr_img,
          model = "grounded-sam",
          prompt = t_prompt,
          box_threshold = box_threshold,
          text_threshold = text_threshold,
          iou_threshold = iou_threshold,
          plot = FALSE,
          verbose = FALSE,
          engine = engine,
          device_id = device_id
        )
        if (!is.null(gs_res$boxes) && nrow(gs_res$boxes) > 0L) {
          pred_boxes <- data.frame(
            index = seq_len(nrow(gs_res$boxes)),
            xmin = as.numeric(gs_res$boxes$xmin),
            ymin = as.numeric(gs_res$boxes$ymin),
            xmax = as.numeric(gs_res$boxes$xmax),
            ymax = as.numeric(gs_res$boxes$ymax),
            cx = (as.numeric(gs_res$boxes$xmin) + as.numeric(gs_res$boxes$xmax)) / 2,
            cy = (as.numeric(gs_res$boxes$ymin) + as.numeric(gs_res$boxes$ymax)) / 2,
            class_id = class_id[1],
            class_name = class_names[1],
            conf = round(as.numeric(gs_res$boxes$score), 3),
            stringsAsFactors = FALSE
          )
          pred_polys <- if (!is.null(gs_res$contours)) gs_res$contours else list()
        }
      } else if (is_yolo) {
        if (isTRUE(verbose)) cli::cli_alert_info("Running YOLO inference ({basename(model[1])}) on {.file {base_stem}}...")
        use_sahi <- isTRUE(dots$sahi)
        use_cond_sahi <- !use_sahi && !is.null(dots$sahi_if) && !isFALSE(dots$sahi_if)
        thresh_sahi <- if (isTRUE(dots$sahi_if)) 0.01 else if (is.numeric(dots$sahi_if)) as.numeric(dots$sahi_if[1]) else 0.01
        if (is.numeric(thresh_sahi) && thresh_sahi > 1) thresh_sahi <- thresh_sahi / 100
        sahi_metric <- if (!is.null(dots$sahi_metric)) dots$sahi_metric[1] else "mean"
        sahi_if_empty <- if (!is.null(dots$sahi_if_empty)) isTRUE(dots$sahi_if_empty) else TRUE
        slice_h <- if (!is.null(dots$slice_height)) as.integer(dots$slice_height[1]) else 640L
        slice_w <- if (!is.null(dots$slice_width)) as.integer(dots$slice_width[1]) else 640L
        can_slice <- (img_w > slice_w || img_h > slice_h)

        .eval_review_sahi <- function(df_boxes) {
          if (!use_cond_sahi || !can_slice) return(list(trigger = FALSE, reason = ""))
          if (is.null(df_boxes) || nrow(df_boxes) == 0L) {
            if (isTRUE(sahi_if_empty)) return(list(trigger = TRUE, reason = "0 objects detected at full image scale"))
            return(list(trigger = FALSE, reason = ""))
          }
          if (!("xmax" %in% names(df_boxes)) || !("xmin" %in% names(df_boxes))) return(list(trigger = FALSE, reason = ""))
          b_areas <- (as.numeric(df_boxes$xmax) - as.numeric(df_boxes$xmin)) * (as.numeric(df_boxes$ymax) - as.numeric(df_boxes$ymin))
          rel_areas <- b_areas / (as.numeric(img_w) * as.numeric(img_h))
          stat_val <- switch(sahi_metric,
                             "mean" = mean(rel_areas, na.rm = TRUE),
                             "median" = stats::median(rel_areas, na.rm = TRUE),
                             "min" = min(rel_areas, na.rm = TRUE),
                             mean(rel_areas, na.rm = TRUE))
          if (is.finite(stat_val) && stat_val <= thresh_sahi) {
            msg <- sprintf("relative %s box area (%.4f%%) <= threshold (%.4f%%)",
                           sahi_metric, stat_val * 100, thresh_sahi * 100)
            return(list(trigger = TRUE, reason = msg))
          }
          list(trigger = FALSE, reason = "")
        }

        if (identical(task, "detect")) {
          if (use_sahi) {
            y_res <- image_detect_sahi(
              curr_img,
              model = model[1],
              slice_height = slice_h,
              slice_width = slice_w,
              overlap_height_ratio = if (!is.null(dots$overlap_height_ratio)) as.numeric(dots$overlap_height_ratio[1]) else 0.20,
              overlap_width_ratio = if (!is.null(dots$overlap_width_ratio)) as.numeric(dots$overlap_width_ratio[1]) else 0.20,
              conf_threshold = box_threshold,
              iou_threshold = iou_threshold,
              match_metric = if (!is.null(dots$match_metric)) dots$match_metric[1] else "IOS",
              match_threshold = if (!is.null(dots$match_threshold)) as.numeric(dots$match_threshold[1]) else 0.50,
              full_image_inference = if (!is.null(dots$full_image_inference)) isTRUE(dots$full_image_inference) else TRUE,
              plot = FALSE,
              verbose = FALSE,
              engine = engine,
              device_id = device_id
            )
          } else {
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
            df_y_chk <- as.data.frame(y_res)
            chk <- .eval_review_sahi(df_y_chk)
            if (isTRUE(chk$trigger)) {
              if (isTRUE(verbose)) cli::cli_alert_info("Triggering SAHI inference for {.file {base_stem}} ({chk$reason})...")
              y_res <- image_detect_sahi(
                curr_img,
                model = model[1],
                slice_height = slice_h,
                slice_width = slice_w,
                overlap_height_ratio = if (!is.null(dots$overlap_height_ratio)) as.numeric(dots$overlap_height_ratio[1]) else 0.20,
                overlap_width_ratio = if (!is.null(dots$overlap_width_ratio)) as.numeric(dots$overlap_width_ratio[1]) else 0.20,
                conf_threshold = box_threshold,
                iou_threshold = iou_threshold,
                match_metric = if (!is.null(dots$match_metric)) dots$match_metric[1] else "IOS",
                match_threshold = if (!is.null(dots$match_threshold)) as.numeric(dots$match_threshold[1]) else 0.50,
                full_image_inference = if (!is.null(dots$full_image_inference)) isTRUE(dots$full_image_inference) else TRUE,
                plot = FALSE,
                verbose = FALSE,
                engine = engine,
                device_id = device_id
              )
            }
          }
          df_y <- as.data.frame(y_res)
          if (nrow(df_y) > 0L) {
            pred_boxes <- data.frame(
              index = seq_len(nrow(df_y)),
              xmin = as.numeric(df_y$xmin),
              ymin = as.numeric(df_y$ymin),
              xmax = as.numeric(df_y$xmax),
              ymax = as.numeric(df_y$ymax),
              cx = (as.numeric(df_y$xmin) + as.numeric(df_y$xmax)) / 2,
              cy = (as.numeric(df_y$ymin) + as.numeric(df_y$ymax)) / 2,
              class_id = if (!is.null(df_y$class_id)) as.integer(df_y$class_id) else class_id[1],
              class_name = if (!is.null(df_y$class_name)) as.character(df_y$class_name) else class_names[1],
              conf = if (!is.null(df_y$conf)) round(as.numeric(df_y$conf), 3) else 1.0,
              stringsAsFactors = FALSE
            )
          }
        } else {
          if (use_sahi) {
            y_res <- image_segment_sahi(
              curr_img,
              model = model[1],
              slice_height = slice_h,
              slice_width = slice_w,
              overlap_height_ratio = if (!is.null(dots$overlap_height_ratio)) as.numeric(dots$overlap_height_ratio[1]) else 0.20,
              overlap_width_ratio = if (!is.null(dots$overlap_width_ratio)) as.numeric(dots$overlap_width_ratio[1]) else 0.20,
              conf_threshold = box_threshold,
              iou_threshold = iou_threshold,
              match_metric = if (!is.null(dots$match_metric)) dots$match_metric[1] else "IOS",
              match_threshold = if (!is.null(dots$match_threshold)) as.numeric(dots$match_threshold[1]) else 0.50,
              full_image_inference = if (!is.null(dots$full_image_inference)) isTRUE(dots$full_image_inference) else TRUE,
              type = "highlight",
              plot = FALSE,
              verbose = FALSE,
              engine = engine,
              device_id = device_id
            )
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
            y_boxes_chk <- if (!is.null(y_res$results)) y_res$results else if (!is.null(y_res$boxes)) y_res$boxes else as.data.frame(y_res)
            chk <- .eval_review_sahi(y_boxes_chk)
            if (isTRUE(chk$trigger)) {
              if (isTRUE(verbose)) cli::cli_alert_info("Triggering SAHI segmentation for {.file {base_stem}} ({chk$reason})...")
              y_res <- image_segment_sahi(
                curr_img,
                model = model[1],
                slice_height = slice_h,
                slice_width = slice_w,
                overlap_height_ratio = if (!is.null(dots$overlap_height_ratio)) as.numeric(dots$overlap_height_ratio[1]) else 0.20,
                overlap_width_ratio = if (!is.null(dots$overlap_width_ratio)) as.numeric(dots$overlap_width_ratio[1]) else 0.20,
                conf_threshold = box_threshold,
                iou_threshold = iou_threshold,
                match_metric = if (!is.null(dots$match_metric)) dots$match_metric[1] else "IOS",
                match_threshold = if (!is.null(dots$match_threshold)) as.numeric(dots$match_threshold[1]) else 0.50,
                full_image_inference = if (!is.null(dots$full_image_inference)) isTRUE(dots$full_image_inference) else TRUE,
                type = "highlight",
                plot = FALSE,
                verbose = FALSE,
                engine = engine,
                device_id = device_id
              )
            }
          }
          y_boxes <- if (!is.null(y_res$results)) y_res$results else if (!is.null(y_res$boxes)) y_res$boxes else as.data.frame(y_res)
          if (!is.null(y_boxes) && nrow(y_boxes) > 0L) {
            df_b <- as.data.frame(y_boxes)
            pred_boxes <- data.frame(
              index = seq_len(nrow(df_b)),
              xmin = as.numeric(df_b$xmin),
              ymin = as.numeric(df_b$ymin),
              xmax = as.numeric(df_b$xmax),
              ymax = as.numeric(df_b$ymax),
              cx = (as.numeric(df_b$xmin) + as.numeric(df_b$xmax)) / 2,
              cy = (as.numeric(df_b$ymin) + as.numeric(df_b$ymax)) / 2,
              class_id = if (!is.null(df_b$class_id)) as.integer(df_b$class_id) else class_id[1],
              class_name = if (!is.null(df_b$class_name)) as.character(df_b$class_name) else class_names[1],
              conf = if (!is.null(df_b$score)) round(as.numeric(df_b$score), 3) else 1.0,
              stringsAsFactors = FALSE
            )
            pred_polys <- if (!is.null(y_res$contours)) y_res$contours else if (!is.null(attr(y_res, "contours"))) attr(y_res, "contours") else list()
          }
        }
      } else if (is_salient) {
        s_res <- image_segment_dl(curr_img, model = model[1], type = "mask", plot = FALSE, verbose = FALSE,
                                  engine = engine, device_id = device_id)
        lbl_mat <- bwlabel_cpp(image_data(s_res) > 0)
        if (max(lbl_mat) > 0) {
          cnts <- contour(lbl_mat)
          b_list <- list()
          for (ci in seq_along(cnts)) {
            c_m <- cnts[[ci]]
            if (nrow(c_m) >= 3L) {
              b_list[[length(b_list) + 1L]] <- data.frame(
                index = length(b_list) + 1L,
                xmin = min(c_m[, 1]), ymin = min(c_m[, 2]),
                xmax = max(c_m[, 1]), ymax = max(c_m[, 2]),
                cx = mean(c_m[, 1]), cy = mean(c_m[, 2]),
                class_id = class_id[1], class_name = class_names[1], conf = 1.0,
                stringsAsFactors = FALSE
              )
              gx <- c_m[, 1] / img_w; gy <- c_m[, 2] / img_h
              pred_polys[[length(pred_polys) + 1L]] <- as.vector(rbind(gx, gy))
            }
          }
          if (length(b_list) > 0L) pred_boxes <- do.call(rbind, b_list)
        }
      } else {
        extra_args <- dots[setdiff(names(dots), c("show_bbox", "show_polygon", "show_points", "show_labels", "show_indices"))]
        res_obj <- if (length(extra_args) > 0L) {
          do.call(analyze_objects, c(list(img = curr_img, plot = FALSE), extra_args))
        } else {
          analyze_objects(curr_img, plot = FALSE)
        }
        if (!is.null(res_obj$results)) {
          df_r <- res_obj$results
          b_list <- list()
          for (ri in seq_len(nrow(df_r))) {
            cx <- df_r$x[ri]; cy <- df_r$y[ri]
            w_r <- sqrt(df_r$area[ri])
            b_list[[length(b_list) + 1L]] <- data.frame(
              index = length(b_list) + 1L,
              xmin = max(0, cx - w_r / 2), ymin = max(0, cy - w_r / 2),
              xmax = min(img_w, cx + w_r / 2), ymax = min(img_h, cy + w_r / 2),
              cx = cx, cy = cy,
              class_id = class_id[1], class_name = class_names[1], conf = 1.0,
              stringsAsFactors = FALSE
            )
          }
          if (length(b_list) > 0L) pred_boxes <- do.call(rbind, b_list)
          pred_polys <- if (!is.null(res_obj$contours)) res_obj$contours else list()
        }
      }

      if (mode == "replace") {
        boxes <- pred_boxes
        polygons <- pred_polys
      } else if (mode == "append") {
        if (nrow(pred_boxes) > 0L) {
          boxes <- rbind(existing_boxes, pred_boxes)
          boxes$index <- seq_len(nrow(boxes))
          polygons <- c(existing_polys, pred_polys)
        }
      }
    }

    # ----------------------------------------------------------------------------
    # PROGRAMMATIC FILTERING & REMOVALS
    # ----------------------------------------------------------------------------
    if (!is.null(remove)) {
      if (identical(remove, "all")) {
        boxes <- data.frame(
          index = integer(0), xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
          cx = numeric(0), cy = numeric(0), class_id = integer(0), class_name = character(0),
          conf = numeric(0), stringsAsFactors = FALSE
        )
        polygons <- list()
      } else if (is.numeric(remove) && nrow(boxes) > 0L) {
        del_idx <- intersect(as.integer(remove), seq_len(nrow(boxes)))
        if (length(del_idx) > 0L) {
          boxes <- boxes[-del_idx, , drop = FALSE]
          if (length(polygons) >= max(del_idx)) polygons <- polygons[-del_idx]
          if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
        }
      }
    }

    if (!is.null(keep) && is.numeric(keep) && nrow(boxes) > 0L) {
      keep_idx <- intersect(as.integer(keep), seq_len(nrow(boxes)))
      boxes <- boxes[keep_idx, , drop = FALSE]
      if (length(polygons) > 0L) polygons <- polygons[keep_idx]
      if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
    }

    # Filter micro-boxes by min_area
    if (min_area > 0 && nrow(boxes) > 0L) {
      b_areas <- (boxes$xmax - boxes$xmin) * (boxes$ymax - boxes$ymin)
      valid_a <- which(b_areas >= min_area)
      if (length(valid_a) < nrow(boxes)) {
        boxes <- boxes[valid_a, , drop = FALSE]
        if (length(polygons) > 0L) polygons <- polygons[valid_a]
        if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
      }
    }

    # Filter edge artifacts if requested
    if (isTRUE(filter_edge) && nrow(boxes) > 0L) {
      bw_px <- boxes$xmax - boxes$xmin
      bh_px <- boxes$ymax - boxes$ymin
      aspect <- pmax(bw_px, bh_px) / pmax(1.0, pmin(bw_px, bh_px))
      is_border <- (boxes$xmin <= 4 | boxes$xmax >= (img_w - 4) | boxes$ymin <= 4 | boxes$ymax >= (img_h - 4))
      keep_e <- !(is_border & (aspect > 3.5 | bw_px < 15 | bh_px < 15))
      boxes <- boxes[keep_e, , drop = FALSE]
      if (length(polygons) > 0L) polygons <- polygons[keep_e]
      if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
    }

    # ----------------------------------------------------------------------------
    # INTERACTIVE EDITOR SESSION
    # ----------------------------------------------------------------------------
    nav_action <- "next"
    skip_save <- FALSE

    if (isTRUE(interactive)) {
      .yolo_render_edit_canvas(
        img = curr_img,
        boxes = boxes,
        polygons = polygons,
        img_name = base_stem,
        split_name = split_name,
        task = task,
        rainbow = rainbow,
        show_bbox = show_bbox,
        show_polygon = show_polygon,
        show_points = show_points,
        show_labels = show_labels,
        show_indices = show_indices,
        point_pch = point_pch,
        point_cex = point_cex
      )

      repeat {
        n_cur <- nrow(boxes)
        cat("\n")
        header_title <- if (total_subset < total_imgs) {
          sprintf("YOLO Dataset Image Editor [%s | %s | %d/%d (Dataset #%d/%d) | %d obj%s]",
                  base_stem, split_name, step_i, total_subset, curr_idx, total_imgs, n_cur, if (n_cur == 1) "" else "s")
        } else {
          sprintf("YOLO Dataset Image Editor [%s | %s | Image %d of %d | %d obj%s]",
                  base_stem, split_name, curr_idx, total_imgs, n_cur, if (n_cur == 1) "" else "s")
        }
        active_cl_name <- if (active_class_idx <= length(class_names)) class_names[active_class_idx] else class_names[1]
        active_cl_id <- if (active_class_idx <= length(class_id)) class_id[active_class_idx] else class_id[1]

        cli::cli_rule(left = header_title)
        if (is_sequential_mode) {
          next_pos <- min(total_subset, step_i + 1L)
          prev_pos <- max(1L, step_i - 1L)
          cli::cli_bullets(c(
            "*" = paste0(cli::col_green("[Enter] / c / n"), ": Save labels and go to {.val NEXT} image ({step_i} -> {next_pos})"),
            "*" = paste0(cli::col_green("p"), ": Save labels and go to {.val PREVIOUS} image ({step_i} -> {prev_pos})"),
            "*" = paste0(cli::col_green("s / skip"), ": Skip this image (leave untouched) and go to {.val NEXT}"),
            "*" = paste0(cli::col_red("rm"), ": ", cli::col_red("Permanently DELETE this image and labels from dataset (with confirmation)")),
            "*" = paste0(cli::col_green("a [class]"), ": Add bounding boxes (click 2 diagonal corners per box; <Esc> when done)"),
            "*" = paste0(cli::col_green("r"), ": Click on 1 or more objects (boxes/dots/masks) to delete them (<Esc> when finished)"),
            "*" = paste0(cli::col_green("d <numbers>"), ": Delete specific objects by index (e.g. {.code d 3} or {.code d 2,5})"),
            "*" = paste0(cli::col_green("class <name/id>"), ": Change active class for new boxes (active: {.val {active_cl_name}} [id: {active_cl_id}])"),
            "*" = paste0(cli::col_green("clear"), ": Remove all objects (set image as background/negative)"),
            "*" = paste0(cli::col_green("goto <number/name>"), ": Jump directly to another image"),
            "*" = paste0(cli::col_green("q / finish"), ": Save current image and finish review session"),
            "*" = paste0(cli::col_green("abort"), ": Exit review session immediately without saving current image")
          ))
        } else {
          cli::cli_bullets(c(
            "*" = paste0(cli::col_green("[Enter] / c"), ": Save and apply changes to label file"),
            "*" = paste0(cli::col_red("rm"), ": ", cli::col_red("Permanently DELETE this image and label from dataset (with confirmation)")),
            "*" = paste0(cli::col_green("n / p"), ": Save and navigate to {.val NEXT} or {.val PREVIOUS} image"),
            "*" = paste0(cli::col_green("a [class]"), ": Add bounding boxes (click 2 diagonal corners per box; <Esc> when done)"),
            "*" = paste0(cli::col_green("r"), ": Click on 1 or more objects to delete them (<Esc> when finished)"),
            "*" = paste0(cli::col_green("d <numbers>"), ": Delete specific objects by index number (e.g. {.code d 3})"),
            "*" = paste0(cli::col_green("class <name/id>"), ": Change active class for new boxes (active: {.val {active_cl_name}} [id: {active_cl_id}])"),
            "*" = paste0(cli::col_green("clear"), ": Remove all objects"),
            "*" = paste0(cli::col_green("q / abort"), ": Cancel editing and leave label file untouched")
          ))
        }

        prompt_txt <- if (is_sequential_mode) {
          if (total_subset < total_imgs) {
            sprintf("Select action [Enter = Save & Next (%d/%d)]: ", step_i, total_subset)
          } else {
            sprintf("Select action [Enter = Save & Next (%d/%d)]: ", curr_idx, total_imgs)
          }
        } else {
          "Select action [Enter = Save]: "
        }

        ans <- tryCatch(
          tolower(trimws(readline(prompt_txt))),
          interrupt = function(e) "q"
        )

        if (ans %in% c("rm", "del_img", "delete_img", "remove_img", "remove_image", "delete_image")) {
          cli::cli_alert_warning("You are about to PERMANENTLY DELETE {.file {basename(target_img_path)}} and its annotations from the dataset!")
          conf_del <- tolower(trimws(readline("Are you sure? Type 'yes' or 'y' to confirm (any other key cancels): ")))
          if (!conf_del %in% c("y", "yes")) {
            cli::cli_alert_info("Deletion cancelled.")
            next
          }
          unlink(target_img_path, force = TRUE)
          unlink(dest_lbl_path, force = TRUE)
          unlink(paste0(dest_lbl_path, ".bak"), force = TRUE)
          rev_cand <- file.path(dir_path, "review", paste0("review_", base_stem, ".jpg"))
          unlink(rev_cand, force = TRUE)
          unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
          unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
          unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
          unlink(file.path(dir_path, "train", "labels.cache"), force = TRUE)
          unlink(file.path(dir_path, "val", "labels.cache"), force = TRUE)
          unlink(file.path(dir_path, "test", "labels.cache"), force = TRUE)

          cli::cli_alert_success("Permanently removed {.file {basename(target_img_path)}} and its label from dataset.")

          reviewed_results[[length(reviewed_results) + 1L]] <- list(
            image_path = target_img_path,
            label_path = dest_lbl_path,
            boxes = data.frame(),
            polygons = list(),
            n_objects = 0L,
            backup_path = NULL,
            status = "deleted"
          )
          skip_save <- TRUE
          nav_action <- if (is_sequential_mode) "next" else "exit"
          break
        }

        if (ans == "" || ans %in% c("c", "continue", "save", "ok", "y", "yes")) {
          nav_action <- if (is_sequential_mode) "next" else "exit"
          break
        }

        if (ans %in% c("n", "next")) {
          nav_action <- "next"
          break
        }

        if (ans %in% c("p", "prev", "previous")) {
          nav_action <- "prev"
          break
        }

        if (ans %in% c("s", "skip")) {
          skip_save <- TRUE
          nav_action <- "next"
          cli::cli_alert_info("Skipped {.file {basename(target_img_path)}} without changes.")
          break
        }

        if (ans %in% c("q", "quit", "finish", "done", "exit")) {
          nav_action <- "exit"
          break
        }

        if (ans %in% c("abort", "cancel")) {
          skip_save <- TRUE
          nav_action <- "exit"
          cli::cli_alert_info("Editing cancelled for {.file {basename(target_img_path)}}.")
          break
        }

        if (grepl("^(goto|jump)\\s+", ans) || (grepl("^\\d+$", ans) && as.integer(ans) <= total_imgs)) {
          target_str <- sub("^(goto|jump)\\s+", "", ans)
          if (grepl("^\\d+$", target_str)) {
            new_i <- as.integer(target_str)
            if (new_i >= 1L && new_i <= total_subset) {
              nav_action <- paste0("goto_", new_i)
              break
            } else if (new_i %in% subset_indices) {
              nav_action <- paste0("goto_", which(subset_indices == new_i)[1])
              break
            }
          } else {
            clean_t <- sub("[./\\]+$", "", target_str)
            target_stem_cand <- tolower(tools::file_path_sans_ext(basename(clean_t)))
            target_base_cand <- tolower(basename(clean_t))
            subset_stems <- tolower(all_stems[subset_indices])
            subset_bases <- tolower(basename(img_files[subset_indices]))
            m_s <- which(subset_stems == target_stem_cand | subset_bases == target_base_cand)
            if (length(m_s) > 0L) {
              nav_action <- paste0("goto_", m_s[1])
              break
            }
            m_t <- which(tolower(all_stems) == target_stem_cand |
                         tolower(basename(img_files)) == target_base_cand |
                         tolower(all_rel) == tolower(clean_t))
            if (length(m_t) > 0L) {
              nav_action <- paste0("goto_raw_", m_t[1])
              break
            }
          }
          cli::cli_alert_warning("Could not find image {.val {target_str}}.")
          next
        }

        if (ans %in% c("clear", "empty", "all")) {
          boxes <- data.frame(
            index = integer(0), xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
            cx = numeric(0), cy = numeric(0), class_id = integer(0), class_name = character(0),
            conf = numeric(0), stringsAsFactors = FALSE
          )
          polygons <- list()
          cli::cli_alert_success("Cleared all objects. Image is now set as negative/background.")
          .yolo_render_edit_canvas(curr_img, boxes, polygons, base_stem, split_name, task, rainbow,
                                   show_bbox, show_polygon, show_points, show_labels, show_indices, point_pch, point_cex)
          next
        }

        # Handle 'class' (change active class for new boxes)
        if (grepl("^(class|cl)\\s+", ans)) {
          cl_arg <- trimws(sub("^(class|cl)\\s+", "", ans))
          if (grepl("^\\d+$", cl_arg)) {
            c_num <- as.integer(cl_arg)
            if (c_num %in% class_id) {
              active_class_idx <- which(class_id == c_num)[1]
            } else if (c_num >= 1L && c_num <= length(class_names)) {
              active_class_idx <- c_num
            } else {
              cli::cli_alert_warning("Class ID {c_num} not found. Available IDs: {.val {class_id}}.")
            }
          } else {
            m_c <- which(tolower(class_names) == tolower(cl_arg))
            if (length(m_c) > 0L) {
              active_class_idx <- m_c[1]
            } else {
              cli::cli_alert_warning("Class {.val {cl_arg}} not found. Available classes: {.val {class_names}}.")
            }
          }
          cli::cli_alert_info("Active class set to: {.val {class_names[active_class_idx]}} (ID: {class_id[active_class_idx]}).")
          next
        } else if (ans %in% c("class", "classes")) {
          cli::cli_alert_info("Available classes:")
          for (cli_i in seq_along(class_names)) {
            mark <- if (cli_i == active_class_idx) " *" else ""
            cli::cli_bullets(c(" " = paste0("[", class_id[cli_i], "] {.val ", class_names[cli_i], "}", mark)))
          }
          ans_cl <- trimws(readline("Enter class name or ID to select: "))
          if (nzchar(ans_cl)) {
            if (grepl("^\\d+$", ans_cl)) {
              c_num <- as.integer(ans_cl)
              if (c_num %in% class_id) active_class_idx <- which(class_id == c_num)[1]
              else if (c_num >= 1L && c_num <= length(class_names)) active_class_idx <- c_num
            } else {
              m_c <- which(tolower(class_names) == tolower(ans_cl))
              if (length(m_c) > 0L) active_class_idx <- m_c[1]
            }
            cli::cli_alert_info("Active class set to: {.val {class_names[active_class_idx]}} (ID: {class_id[active_class_idx]}).")
          }
          next
        }

        # Handle 'd <indices>'
        if (grepl("^d\\s+", ans)) {
          raw_idx <- sub("^d\\s+", "", ans)
          nums <- as.integer(unlist(strsplit(raw_idx, "[,;\\s]+")))
          nums <- nums[!is.na(nums)]
          del_idx <- intersect(nums, seq_len(nrow(boxes)))
          if (length(del_idx) > 0L) {
            boxes <- boxes[-del_idx, , drop = FALSE]
            if (length(polygons) >= max(del_idx)) polygons <- polygons[-del_idx]
            if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
            cli::cli_alert_success("Deleted object(s): {.val {del_idx}}. {nrow(boxes)} remaining.")
            .yolo_render_edit_canvas(curr_img, boxes, polygons, base_stem, split_name, task, rainbow,
                                     show_bbox, show_polygon, show_points, show_labels, show_indices, point_pch, point_cex)
          } else {
            cli::cli_alert_warning("No matching object indices found.")
          }
          next
        }

        # Handle 'r' (click to delete)
        if (ans %in% c("r", "remove", "del", "delete")) {
          if (nrow(boxes) == 0L) {
            cli::cli_alert_warning("No objects remaining to remove.")
            next
          }
          cli::cli_alert_info("Click on objects (bounding boxes, center dots, or masks) to delete them.")
          cli::cli_alert_info("Press <Esc> or right-click when finished.")
          del_pts <- tryCatch(
            graphics::locator(n = 200, type = "p", pch = 4, col = "#FF0033", cex = 1.6, lwd = 2.5),
            error = function(e) NULL
          )
          if (!is.null(del_pts) && length(del_pts$x) > 0L) {
            del_indices <- integer(0)
            for (pt_i in seq_along(del_pts$x)) {
              cx <- del_pts$x[pt_i]; cy <- del_pts$y[pt_i]
              best_idx <- NULL
              best_dist <- Inf
              for (bi in seq_len(nrow(boxes))) {
                if (bi %in% del_indices) next
                bx <- boxes[bi, ]
                d_cent <- (bx$cx - cx)^2 + (bx$cy - cy)^2
                inside <- (cx >= bx$xmin && cx <= bx$xmax && cy >= bx$ymin && cy <= bx$ymax)
                if (inside) {
                  if (d_cent < best_dist) {
                    best_dist <- d_cent
                    best_idx <- bi
                  }
                } else if (d_cent < (min(img_w, img_h) * 0.08)^2 && d_cent < best_dist) {
                  best_dist <- d_cent
                  best_idx <- bi
                }
              }
              if (!is.null(best_idx)) del_indices <- c(del_indices, best_idx)
            }
            del_indices <- unique(del_indices)
            if (length(del_indices) > 0L) {
              boxes <- boxes[-del_indices, , drop = FALSE]
              if (length(polygons) >= max(del_indices)) polygons <- polygons[-del_indices]
              if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
              cli::cli_alert_success("Removed {length(del_indices)} object(s). {nrow(boxes)} remaining.")
              .yolo_render_edit_canvas(curr_img, boxes, polygons, base_stem, split_name, task, rainbow,
                                       show_bbox, show_polygon, show_points, show_labels, show_indices, point_pch, point_cex)
            } else {
              cli::cli_alert_warning("No object was close to the clicked point(s).")
            }
          }
          next
        }

        # Handle 'a' (click to add bounding boxes continuously)
        if (ans %in% c("a", "add", "new") || grepl("^(a|add)\\s+", ans)) {
          add_cl_idx <- active_class_idx
          if (grepl("^(a|add)\\s+", ans)) {
            cl_arg <- trimws(sub("^(a|add)\\s+", "", ans))
            if (grepl("^\\d+$", cl_arg)) {
              c_num <- as.integer(cl_arg)
              if (c_num %in% class_id) add_cl_idx <- which(class_id == c_num)[1]
              else if (c_num >= 1L && c_num <= length(class_names)) add_cl_idx <- c_num
            } else {
              m_c <- which(tolower(class_names) == tolower(cl_arg))
              if (length(m_c) > 0L) add_cl_idx <- m_c[1]
            }
          }
          target_cl_name <- if (add_cl_idx <= length(class_names)) class_names[add_cl_idx] else class_names[1]
          target_cl_id <- if (add_cl_idx <= length(class_id)) class_id[add_cl_idx] else class_id[1]

          cli::cli_alert_info("Adding bounding boxes for {.val {target_cl_name}} (2 diagonal clicks per box).")
          cli::cli_alert_info("Click corners continuously. Press <Esc> or right-click when finished.")

          n_added <- 0L
          repeat {
            new_pts <- tryCatch(
              graphics::locator(n = 2, type = "p", pch = 3, col = "#00FF66", cex = 1.6, lwd = 2.5),
              error = function(e) NULL
            )
            if (is.null(new_pts) || length(new_pts$x) < 2L) {
              break
            }
            x1 <- min(new_pts$x); x2 <- max(new_pts$x)
            y1 <- min(new_pts$y); y2 <- max(new_pts$y)

            # Skip accidental clicks on virtually the same point
            if (abs(x2 - x1) < 2 || abs(y2 - y1) < 2) {
              cli::cli_alert_warning("Clicked points were too close. Box skipped.")
              next
            }

            n_added <- n_added + 1L
            new_idx <- nrow(boxes) + 1L
            new_row <- data.frame(
              index = new_idx,
              xmin = x1, ymin = y1, xmax = x2, ymax = y2,
              cx = (x1 + x2) / 2, cy = (y1 + y2) / 2,
              class_id = target_cl_id,
              class_name = target_cl_name,
              conf = 1.0,
              stringsAsFactors = FALSE
            )
            boxes <- rbind(boxes, new_row)
            boxes$index <- seq_len(nrow(boxes))
            polygons[[length(polygons) + 1L]] <- NULL

            # Draw the box immediately on the screen in real-time
            graphics::rect(x1, y1, x2, y2, border = "#00FF66", lwd = 2.5, lty = 1)
            lbl_txt <- paste0("#", new_idx, " ", target_cl_name)
            graphics::text(x1, y1 - 8, labels = lbl_txt, col = "#00FF66", font = 2, cex = 0.9, adj = c(0, 1))

            cli::cli_alert_success("Box #{new_idx} added ({target_cl_name}). Click next 2 corners, or <Esc> to finish.")
          }

          if (n_added > 0L) {
            cli::cli_alert_success("Finished adding boxes! Added {n_added} new box{?es} ({nrow(boxes)} total).")
            .yolo_render_edit_canvas(curr_img, boxes, polygons, base_stem, split_name, task, rainbow,
                                     show_bbox, show_polygon, show_points, show_labels, show_indices, point_pch, point_cex)
          } else {
            cli::cli_alert_info("No boxes added.")
          }
          next
        }
      }
    }

    # ----------------------------------------------------------------------------
    # COMMIT / WRITE EDITS TO DISK
    # ----------------------------------------------------------------------------
    backup_file <- NULL
    if (!skip_save) {
      if (isTRUE(backup) && file.exists(dest_lbl_path)) {
        backup_file <- paste0(dest_lbl_path, ".bak")
        file.copy(dest_lbl_path, backup_file, overwrite = TRUE)
      }

      txt_lines <- character(0)
      n_final <- nrow(boxes)

      if (n_final > 0L) {
        if (identical(task, "segment") && length(polygons) == n_final && any(!vapply(polygons, is.null, logical(1)))) {
          # Write segmentation lines
          for (k in seq_len(n_final)) {
            cid <- boxes$class_id[k]
            poly <- polygons[[k]]
            if (!is.null(poly) && length(poly) >= 6L) {
              px <- if (max(poly) <= 1.05) poly[c(TRUE, FALSE)] else poly[c(TRUE, FALSE)] / img_w
              py <- if (max(poly) <= 1.05) poly[c(FALSE, TRUE)] else poly[c(FALSE, TRUE)] / img_h
              interleaved <- as.vector(rbind(round(px, 6L), round(py, 6L)))
              txt_lines <- c(txt_lines, paste(c(cid, sprintf("%.6f", interleaved)), collapse = " "))
            } else {
              # Box fallback to 4 corners
              bx <- boxes[k, ]
              x1 <- round(max(0, bx$xmin) / img_w, 6L); x2 <- round(min(img_w, bx$xmax) / img_w, 6L)
              y1 <- round(max(0, bx$ymin) / img_h, 6L); y2 <- round(min(img_h, bx$ymax) / img_h, 6L)
              poly_pts <- c(x1, y1, x2, y1, x2, y2, x1, y2)
              txt_lines <- c(txt_lines, paste(c(cid, sprintf("%.6f", poly_pts)), collapse = " "))
            }
          }
        } else {
          # Write detection bounding boxes (cx cy w h)
          for (k in seq_len(n_final)) {
            bx <- boxes[k, ]
            cid <- bx$class_id
            cx <- round(((bx$xmin + bx$xmax) / 2) / img_w, 6L)
            cy <- round(((bx$ymin + bx$ymax) / 2) / img_h, 6L)
            bw <- round((bx$xmax - bx$xmin) / img_w, 6L)
            bh <- round((bx$ymax - bx$ymin) / img_h, 6L)
            txt_lines <- c(txt_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, cx, cy, bw, bh))
          }
        }
      }

      # Write updated label file
      writeLines(txt_lines, dest_lbl_path)

      # Update review preview if review directory exists
      rev_dir <- file.path(dir_path, "review")
      if (isTRUE(review) && dir.exists(rev_dir)) {
        rev_out <- file.path(rev_dir, paste0("review_", base_stem, ".jpg"))
        tryCatch({
          out_w <- 1200L
          out_h <- max(100L, as.integer(round(1200 * img_h / img_w)))
          grDevices::jpeg(rev_out, width = out_w, height = out_h, quality = 92, res = 120)
          .yolo_render_edit_canvas(
            img = curr_img,
            boxes = boxes,
            polygons = polygons,
            img_name = base_stem,
            split_name = split_name,
            task = task,
            rainbow = rainbow,
            show_bbox = show_bbox,
            show_polygon = show_polygon,
            show_points = show_points,
            show_labels = show_labels,
            show_indices = show_indices,
            point_pch = point_pch,
            point_cex = point_cex
          )
          grDevices::dev.off()
        }, error = function(e) NULL)
      }

      if (isTRUE(verbose)) {
        cli::cli_alert_success("Saved {.file {basename(dest_lbl_path)}} ({n_final} object{if (n_final == 1) '' else 's'}).")
        if (!is.null(backup_file)) {
          cli::cli_alert_info("Backup saved to {.file {basename(backup_file)}}.")
        }
      }

      session_saved <- session_saved + 1L
      reviewed_results[[length(reviewed_results) + 1L]] <- list(
        image_path = target_img_path,
        label_path = dest_lbl_path,
        boxes = boxes,
        polygons = polygons,
        n_objects = n_final,
        backup_path = backup_file,
        status = "saved"
      )
    }

    # Navigation logic
    if (nav_action == "next") {
      if (step_i >= total_subset) {
        if (isTRUE(verbose)) {
          if (total_subset < total_imgs) {
            cli::cli_alert_success("Reached end of selected images ({total_subset}/{total_subset} images).")
          } else {
            cli::cli_alert_success("Reached end of dataset ({total_imgs}/{total_imgs} images).")
          }
        }
        break
      }
      step_i <- step_i + 1L
    } else if (nav_action == "prev") {
      if (step_i <= 1L) {
        cli::cli_alert_warning("Already at the first image of selection (1/{total_subset}).")
      } else {
        step_i <- step_i - 1L
      }
    } else if (grepl("^goto_raw_\\d+$", nav_action)) {
      raw_num <- as.integer(sub("^goto_raw_", "", nav_action))
      if (raw_num %in% subset_indices) {
        step_i <- which(subset_indices == raw_num)[1]
      } else {
        subset_indices <- c(subset_indices, raw_num)
        total_subset <- length(subset_indices)
        step_i <- total_subset
      }
    } else if (grepl("^goto_\\d+$", nav_action)) {
      target_num <- as.integer(sub("^goto_", "", nav_action))
      if (target_num >= 1L && target_num <= total_subset) {
        step_i <- target_num
      } else if (target_num %in% subset_indices) {
        step_i <- which(subset_indices == target_num)[1]
      } else {
        cli::cli_alert_warning("Target {target_num} is not in the current image selection.")
      }
    } else if (nav_action == "exit") {
      break
    } else {
      break
    }
  }

  if (length(reviewed_results) == 1L) {
    invisible(reviewed_results[[1]])
  } else {
    invisible(reviewed_results)
  }
}

#' @rdname yolo_dataset_review
#' @export
yolo_dataset_edit <- yolo_dataset_review

#' @rdname yolo_dataset_review
#' @export
yolo_image_review <- yolo_dataset_review

#' @rdname yolo_dataset_review
#' @export
yolo_image_edit <- yolo_dataset_review


.yolo_render_edit_canvas <- function(img,
                                     boxes,
                                     polygons,
                                     img_name = "",
                                     split_name = "",
                                     task = "detect",
                                     rainbow = TRUE,
                                     show_bbox = TRUE,
                                     show_polygon = TRUE,
                                     show_points = TRUE,
                                     show_labels = TRUE,
                                     show_indices = TRUE,
                                     point_pch = 19,
                                     point_cex = 1.6) {
  dims <- dim(img)
  img_w <- dims[1]
  img_h <- dims[2]
  n_objs <- if (is.data.frame(boxes)) nrow(boxes) else 0L

  title_str <- if (n_objs == 0L) {
    sprintf("Image: %s [%s] | 0 objects (background)", img_name, split_name)
  } else {
    sprintf("Image: %s [%s] | %d objects | Task: %s", img_name, split_name, n_objs, toupper(task))
  }
  sub_str <- "[r] remove  |  [a] add box  |  [d <numbers>] delete  |  [Enter/c] save  |  [q] cancel"

  plot(
    img,
    main = title_str,
    sub = sub_str,
    col.main = "#006622",
    font.main = 2,
    cex.main = 1.15,
    col.sub = "#555555",
    font.sub = 3,
    cex.sub = 0.85,
    mar = c(2.2, 1.5, 4.0, 1.5)
  )

  if (n_objs == 0L) {
    return(invisible(NULL))
  }

  cols_border <- if (isTRUE(rainbow)) {
    grDevices::rainbow(n_objs, s = 0.85, v = 0.95)
  } else {
    rep("#00FFCC", n_objs)
  }
  cols_poly <- grDevices::adjustcolor(cols_border, alpha.f = 0.35)

  # Draw polygons if available
  if (isTRUE(show_polygon) && length(polygons) > 0L) {
    for (i in seq_along(polygons)) {
      poly <- polygons[[i]]
      if (!is.null(poly) && length(poly) >= 6L) {
        px <- if (max(poly) <= 1.05) poly[c(TRUE, FALSE)] * img_w else poly[c(TRUE, FALSE)]
        py <- if (max(poly) <= 1.05) poly[c(FALSE, TRUE)] * img_h else poly[c(FALSE, TRUE)]
        cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
        cp <- cols_poly[((i - 1L) %% length(cols_poly)) + 1L]
        graphics::polygon(px, py, col = cp, border = cb, lwd = 2)
      }
    }
  }

  # Draw bounding boxes
  if (isTRUE(show_bbox) && nrow(boxes) > 0L) {
    for (i in seq_len(nrow(boxes))) {
      bx <- boxes[i, ]
      cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
      graphics::rect(bx$xmin, bx$ymin, bx$xmax, bx$ymax, border = cb, lwd = 2)
    }
  }

  # Draw center points (solid colored circle with dark border)
  if (isTRUE(show_points) && nrow(boxes) > 0L) {
    for (i in seq_len(nrow(boxes))) {
      bx <- boxes[i, ]
      cx <- (bx$xmin + bx$xmax) / 2
      cy <- (bx$ymin + bx$ymax) / 2
      cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
      graphics::points(cx, cy, pch = point_pch, col = cb, cex = point_cex)
      graphics::points(cx, cy, pch = 1, col = "#000000", cex = point_cex, lwd = 1.2)
    }
  }

  # Draw badges with index numbers (#1, #2...) and labels
  if (nrow(boxes) > 0L && (isTRUE(show_labels) || isTRUE(show_indices))) {
    for (i in seq_len(nrow(boxes))) {
      bx <- boxes[i, ]
      cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
      tag_parts <- character(0)
      if (isTRUE(show_indices)) tag_parts <- c(tag_parts, paste0("#", i))
      if (isTRUE(show_labels) && !is.null(bx$class_name) && nzchar(bx$class_name)) {
        tag_parts <- c(tag_parts, bx$class_name)
      }
      tag <- paste(tag_parts, collapse = " ")

      if (nzchar(tag)) {
        tw <- graphics::strwidth(tag, cex = 0.72)
        th <- graphics::strheight(tag, cex = 0.72)
        y_top <- max(th + 6, bx$ymin)
        graphics::rect(bx$xmin, y_top - th - 5, bx$xmin + tw + 6, y_top,
                       col = grDevices::adjustcolor(cb, alpha.f = 0.85), border = NA)
        graphics::text(bx$xmin + 3, y_top - 2.5, labels = tag, col = "#FFFFFF",
                       cex = 0.72, font = 2, adj = c(0, 1))
      }
    }
  }

  invisible(NULL)
}

#' @title Plot and Preview YOLO Dataset Annotations and Bounding Boxes
#' @name plot.yolo_dataset
#' @aliases plot.yolo_dataset_aug
#' @description
#' Visualizes sample images and their corresponding YOLO bounding boxes or
#' segmentation masks from a YOLO dataset or augmented dataset. Allows inspecting
#' bounding box transformations (rotations, flips) and image augmentations.
#'
#' @param x A `yolo_dataset` or `yolo_dataset_aug` object (or dataset directory path).
#' @param show_bbox Logical. Whether to draw bounding boxes and annotations. Defaults to `TRUE`. Supported alias: `bbox`.
#' @param show_polygon Logical. For segmentation datasets, whether to draw polygons. Defaults to `TRUE`. Supported alias: `polygon`.
#' @param show_points Logical. If `TRUE`, draws points at the centroid of each bounding box or polygon using the same color as the box. Defaults to `FALSE`. Supported alias: `points`.
#' @param show_labels Logical. Whether to display class label names. Defaults to `TRUE`. Supported aliases: `labels`, `show_class`.
#' @param bbox,polygon,points,labels Backward-compatible aliases for `show_bbox`, `show_polygon`, `show_points`, and `show_labels`.
#' @param n Integer. Number of images to plot. Defaults to `4L`. Supported aliases: `n_images`, `num_images`.
#' @param which Optional integer vector of image indices (e.g. `1:4`) or character vector of
#'   image names/stems (e.g. `c("img_3_aug1", "img_553_aug2")`) to plot specific images.
#' @param images Alias for `which`.
#' @param split Character. Dataset split to preview (`"train"`, `"val"`, or `"test"`). Defaults to `"train"`.
#' @param only_augmented Logical. If `TRUE`, only previews augmented images (containing `_aug`).
#'   Defaults to `TRUE` if `x` is an augmented dataset (`yolo_dataset_aug`), and `FALSE` otherwise.
#' @param only_original Logical. If `TRUE`, only previews original pristine images (excluding synthetic
#'   augmentations `_aug[0-9]+`, synthetic backgrounds `synth_bg_`, external backgrounds `bg_`, and scrambled backgrounds `scramble_`).
#'   Defaults to `FALSE` (or `NULL`).
#' @param pattern Optional regex pattern to filter image filenames.
#' @param random Logical. Whether to sample images randomly. Defaults to `TRUE` when `which` is `NULL`.
#' @param seed Optional integer seed for reproducible random selection.
#' @param mfrow Numeric vector of length 2 giving `c(rows, cols)` for the multi-panel layout.
#' @param rainbow Logical. If `TRUE`, assigns distinct rainbow colors to each object. Defaults to `TRUE`.
#' @param col Custom border color for bounding boxes/polygons.
#' @param lwd Numeric. Line width for bounding boxes. Defaults to `1.5`.
#' @param fill Fill color for segmentation polygons. Defaults to semi-transparent color.
#' @param label_bg Logical. Whether to draw a background badge behind label text. Defaults to `FALSE`.
#' @param cex Numeric. Text scale for labels. Defaults to `0.6`.
#' @param title Logical or character. Title for each panel. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [yolo_dataset_preview()].
#' @return An invisible list of the sampled image file paths.
#' @method plot yolo_dataset
#' @export
plot.yolo_dataset <- function(x,
                              show_bbox = TRUE,
                              show_polygon = TRUE,
                              show_points = FALSE,
                              show_labels = TRUE,
                              bbox = show_bbox,
                              polygon = show_polygon,
                              points = show_points,
                              labels = show_labels,
                              n = 4L,
                              which = NULL,
                              images = NULL,
                              split = c("train", "val", "test"),
                              only_augmented = NULL,
                              only_original = NULL,
                              pattern = NULL,
                              random = is.null(which) && is.null(images),
                              seed = NULL,
                              mfrow = NULL,
                              rainbow = TRUE,
                              col = NULL,
                              lwd = 1.5,
                              fill = NULL,
                              label_bg = FALSE,
                              cex = 0.6,
                              title = TRUE,
                              ...) {
  if (!missing(bbox)) show_bbox <- bbox
  if (!missing(polygon)) show_polygon <- polygon
  if (!missing(points)) show_points <- points
  if (!missing(labels)) show_labels <- labels
  yolo_dataset_preview(
    dir = x,
    split = split,
    show_bbox = show_bbox,
    show_polygon = show_polygon,
    show_points = show_points,
    show_labels = show_labels,
    n = n,
    which = which,
    images = images,
    only_augmented = only_augmented,
    only_original = only_original,
    pattern = pattern,
    random = random,
    seed = seed,
    mfrow = mfrow,
    rainbow = rainbow,
    col = col,
    lwd = lwd,
    fill = fill,
    label_bg = label_bg,
    cex = cex,
    title = title,
    ...
  )
}
#' @method plot yolo_dataset_aug
#' @export
plot.yolo_dataset_aug <- function(x, ...) {
  plot.yolo_dataset(x, ...)
}
#' @method print yolo_dataset_aug
#' @export
print.yolo_dataset_aug <- function(x, ...) {
  print.yolo_dataset(x, ...)
}
#' @export
print.yolo_dataset <- function(x, ...) {
  if (inherits(x, "yolo_dataset_aug") || isTRUE(x$is_augmented)) {
    cli::cli_h3("YOLO Augmented Dataset (pliman)")
  } else {
    cli::cli_h3("YOLO Dataset (pliman)")
  }
  bullets <- c(
    "*" = "Directory: {.path {x$dir}}",
    "*" = "YAML config: {.path {x$yaml_file}}",
    "*" = "Train images: {length(x$train_images)}",
    "*" = "Validation images: {length(x$val_images)}"
  )
  if (length(x$test_images) > 0L) {
    bullets <- c(bullets, "*" = "Test images: {length(x$test_images)}")
  }
  bullets <- c(bullets, "*" = "Total objects: {x$n_objects}")
  if (!is.null(x$total_augmented)) {
    bullets <- c(bullets, "*" = "Augmented variants: {x$total_augmented}")
  }
  if (!is.null(x$review_dir) && dir.exists(x$review_dir)) {
    bullets <- c(bullets, "*" = "Review previews: {.path {x$review_dir}}")
  }
  cli::cli_bullets(bullets)
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
#' @param as_job Logical. If `TRUE` and running inside RStudio, launches background image addition
#'   as an RStudio Background Job in the Jobs tab, keeping the current R session free. Defaults to `FALSE`.
#' @param verbose Logical. Whether to show progress messages. Defaults to `TRUE`.
#'
#' @return Invisible list with summary of added files.
#' @export
#' @examplesIf interactive()
#' # Using a local NDJSON downloaded from Ultralytics Hub:
#' yolo_dataset_add_background(
#'   dir = "yolo_toras",
#'   source = "african-wildlife.ndjson",
#'   n = 20
#' )
yolo_dataset_add_background <- function(dir = "yolo_dataset",
                                        source,
                                        n = 20L,
                                        train_prop = 0.80,
                                        prefix = "bg_",
                                        as_job = FALSE,
                                        verbose = TRUE) {
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with background image addition in current session.")
    } else {
      mc <- match.call(expand.dots = TRUE)
      eval_env <- parent.frame()
      evaluated_args <- list()
      arg_names <- setdiff(names(mc), c("", "as_job"))
      for (arg_name in arg_names) {
        arg_val <- tryCatch(eval(mc[[arg_name]], envir = eval_env), error = function(e) mc[[arg_name]])
        evaluated_args[[arg_name]] <- arg_val
      }
      evaluated_args$as_job <- FALSE

      job_dir <- file.path(tempdir(), "yolo_jobs")
      if (!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("yolo_add_bg_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_add_bg_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (inherits(dir, "yolo_dataset")) basename(dir$dir)
        else if (is.character(dir) && nzchar(dir[1])) basename(dir[1])
        else "dataset"
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Add Background: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_add_background()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Add Background - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_add_background, job_info$args)\n",
        "cat('\\n[OK] Add Background Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Add Background (pliman)")
      cli::cli_alert_success("Background image addition has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target directory: {.val {target_name}}",
        "Count: {.val {n}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }
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
  if (dir.exists(file.path(dir, "train", "images"))) {
    train_img_dir <- file.path(dir, "train", "images")
    val_img_dir   <- file.path(dir, "val", "images")
    train_lbl_dir <- file.path(dir, "train", "labels")
    val_lbl_dir   <- file.path(dir, "val", "labels")
  } else {
    train_img_dir <- file.path(dir, "images", "train")
    val_img_dir   <- file.path(dir, "images", "val")
    train_lbl_dir <- file.path(dir, "labels", "train")
    val_lbl_dir   <- file.path(dir, "labels", "val")
  }
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
  synth_patt_names <- c("synth", "synthetic", "crop_rows", "field", "plantation", "rows",
                        "soil", "earth", "dirt", "straw", "mulch", "residue",
                        "pasture", "grass", "lawn", "forest", "canopy", "woods",
                        "mountain", "sea", "water", "sky", "clouds", "wood",
                        "lab", "paper", "board", "bench", "germitest",
                        "leaf_litter", "leaves", "litter", "mulch_film", "plastic", "plastic_mulch",
                        "gravel", "concrete", "floor", "composite", "field_mixed", "mixed",
                        "waves", "wave", "gradient", "plasma", "perlin", "blobs", "texture", "random")
  is_synth <- is.character(source) && length(source) == 1L &&
              tolower(source) %in% synth_patt_names
  if (is_synth) {
    patt <- tolower(source)
    if (patt %in% c("synth", "synthetic")) patt <- "random"
    if (patt == "wave") patt <- "waves"
    return(invisible(yolo_background_synth(
      dir = dir,
      n = n,
      pattern = patt,
      train_prop = train_prop,
      prefix = prefix,
      verbose = verbose
    )))
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
    uids_zip <- tolower(gsub("-", "", uuid(length(selected_entries))))
    for (k in seq_along(selected_entries)) {
      is_tr <- (k <= n_tr)
      dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
      dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir
      src_f <- file.path(tmp_dir, selected_entries[k])
      ext <- tools::file_ext(selected_entries[k])
      if (!nzchar(ext)) ext <- "jpg"
      base_fn <- sprintf("%s%s", prefix, uids_zip[k])
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
      uids_nd <- tolower(gsub("-", "", uuid(length(selected))))
      for (k in seq_along(selected)) {
        rec <- selected[[k]]
        is_tr <- (k <= n_train)
        dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
        dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir
        base_fn <- sprintf("%s%s", prefix, uids_nd[k])
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
    uids_local <- tolower(gsub("-", "", uuid(length(selected_files))))
    for (k in seq_along(selected_files)) {
      is_tr <- (k <= n_train)
      dest_img_dir <- if (is_tr) train_img_dir else val_img_dir
      dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir
      ext <- tools::file_ext(selected_files[k])
      if (!nzchar(ext)) ext <- "jpg"
      base_fn <- sprintf("%s%s", prefix, uids_local[k])
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

#' @title Undo and Remove Added Background Images from a YOLO Dataset
#' @name yolo_dataset_add_background_undo
#' @aliases yolo_dataset_background_undo yolo_dataset_bg_undo yolo_dataset_add_bg_undo yolo_dataset_clean_background yolo_dataset_clean_bg
#' @description
#' `yolo_dataset_add_background_undo()` removes all negative background images and matching
#' empty label files (matching `^bg_` or a custom `prefix`) previously added by
#' [yolo_dataset_add_background()], reverting the dataset back to its pristine state.
#'
#' It also automatically invalidates stale YOLO cache files (`train.cache`, `val.cache`, `test.cache`, `labels.cache`).
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path. Defaults to `"."`.
#' @param prefix Character. Filename prefix of the background files to remove (default: `"bg_"`).
#'   If set to `"all"`, cleans both generic (`"^bg_"`) and synthetic (`"^synth_bg_"`) background files.
#' @param pattern Optional character regular expression. If provided, overrides `prefix` to match files.
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Add negative background images to a YOLO dataset
#' yolo_dataset_add_background("my_dataset", source = "background_folder", n = 20)
#'
#' # 2. Undo and remove added background images
#' yolo_dataset_add_background_undo("my_dataset")
yolo_dataset_add_background_undo <- function(data = ".",
                                             prefix = "bg_",
                                             pattern = NULL,
                                             verbose = TRUE) {
  dir_path <- tryCatch(.yolo_get_dir(data), error = function(e) {
    if (is.list(data)) {
      if (!is.null(data$dir) && is.character(data$dir) && nzchar(data$dir[1])) {
        data$dir[1]
      } else if (!is.null(data$yaml_file) && is.character(data$yaml_file) && nzchar(data$yaml_file[1])) {
        dirname(data$yaml_file[1])
      } else if (length(data) > 0L && is.character(data[[1]]) && nzchar(data[[1]])) {
        data[[1]]
      } else {
        "."
      }
    } else if (is.character(data) && length(data) >= 1L && nzchar(data[1])) {
      data[1]
    } else {
      "."
    }
  })
  raw_dir <- dir_path
  dir_path <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    candidates <- c(
      file.path(getwd(), raw_dir),
      file.path(tempdir(), raw_dir),
      file.path(pliman_model_dir(), raw_dir)
    )
    for (cand in candidates) {
      if (dir.exists(cand)) {
        dir_path <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }
  task <- tryCatch(.yolo_detect_task(dir_path), error = function(e) "detect")
  classes <- tryCatch(.yolo_get_classes(dir_path), error = function(e) character())

  # Resolve regex pattern for matching files
  patt_regex <- if (!is.null(pattern) && is.character(pattern) && nzchar(pattern[1])) {
    pattern[1]
  } else if (identical(tolower(prefix[1]), "all")) {
    "^(bg_|synth_bg_)"
  } else if (length(prefix) > 1L) {
    paste0("^(", paste(prefix, collapse = "|"), ")")
  } else {
    paste0("^", prefix[1])
  }

  img_ext_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  n_removed_img <- 0L
  n_removed_lbl <- 0L
  n_orig_img <- 0L

  # Search images in images/ or dir_path (covers images/train, train/images, etc.)
  img_dir <- file.path(dir_path, "images")
  search_img_dir <- if (dir.exists(img_dir)) img_dir else dir_path
  all_imgs <- list.files(search_img_dir, pattern = img_ext_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  bg_imgs <- all_imgs[grepl(patt_regex, basename(all_imgs), ignore.case = TRUE)]
  if (length(bg_imgs) > 0L) {
    unlink(bg_imgs, force = TRUE)
    n_removed_img <- length(bg_imgs)
  }
  n_orig_img <- length(all_imgs) - n_removed_img

  # Search labels in labels/, train/labels, val/labels, etc.
  lbl_dirs <- c(
    file.path(dir_path, "labels"),
    file.path(dir_path, "train", "labels"),
    file.path(dir_path, "val", "labels"),
    file.path(dir_path, "test", "labels")
  )
  lbl_dirs <- lbl_dirs[dir.exists(lbl_dirs)]
  search_lbl_dir <- if (length(lbl_dirs) > 0L) lbl_dirs else dir_path
  all_lbls <- list.files(search_lbl_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  bg_lbls <- all_lbls[grepl(patt_regex, basename(all_lbls), ignore.case = TRUE)]
  if (length(bg_lbls) > 0L) {
    unlink(bg_lbls, force = TRUE)
    n_removed_lbl <- length(bg_lbls)
  }

  # Invalidate stale cache files across all potential locations
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
  unlink(file.path(dir_path, "train", "labels.cache"), force = TRUE)
  unlink(file.path(dir_path, "val", "labels.cache"), force = TRUE)
  unlink(file.path(dir_path, "test", "labels.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No background files matching {.val {patt_regex}} found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned background files from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "*" = "Background images removed: {n_removed_img}",
        "*" = "Background empty label files removed: {n_removed_lbl}",
        "*" = "Remaining images in dataset: {n_orig_img}",
        "v" = "YOLO label caches invalidated."
      ))
    }
  }

  res <- tryCatch(.yolo_wrap_dataset(dir_path, task = task, classes = classes), error = function(e) {
    list(dir = dir_path, images_removed = n_removed_img, labels_removed = n_removed_lbl, remaining_images = n_orig_img)
  })
  invisible(res)
}

#' @export
yolo_dataset_background_undo <- yolo_dataset_add_background_undo

#' @export
yolo_dataset_bg_undo <- yolo_dataset_add_background_undo

#' @export
yolo_dataset_add_bg_undo <- yolo_dataset_add_background_undo

#' @export
yolo_dataset_clean_background <- yolo_dataset_add_background_undo

#' @export
yolo_dataset_clean_bg <- yolo_dataset_add_background_undo

#' @title Generate Synthetic Background Images for YOLO Datasets
#' @name yolo_background_synth
#' @aliases yolo_synth_background yolo_dataset_synth_background
#' @description
#' `yolo_background_synth()` creates realistic synthetic background images featuring agricultural crop rows,
#' tilled soil with aggregate crumbs, straw/mulch crop residues, dense pasture turf, forest canopies,
#' mountain relief silhouettes, sea/water surfaces with specular reflections, atmospheric sky/clouds,
#' laboratory phenotyping bench/paper with optional calibration grids, fallen leaf litter,
#' agricultural plastic mulch film, greenhouse gravel/concrete floors, multi-layer composite fields,
#' wood grain, sinusoidal waves, organic color gradients, Perlin/plasma noise, or textured surfaces.
#'
#' In YOLO detection and segmentation models, adding negative background images (images without any target
#' objects and with empty `.txt` annotation files) is a recommended practice (typically 1% to 10% of the dataset)
#' to drastically minimize False Positive (FP) detections on complex, non-target backgrounds.
#'
#' @param dir Character or `yolo_dataset` object. Target YOLO dataset directory or output folder.
#'   If `NULL` (default), generated images are returned directly as a list of pliman `Image` objects.
#'   If a YOLO dataset is supplied, images are split between `images/train/` and `images/val/` according
#'   to `train_prop`, and matching empty `.txt` label files are created automatically.
#' @param n Integer. Number of synthetic background images to generate. Defaults to `10L`.
#' @param width,height Integer. Dimensions of generated images in pixels. Defaults to `640L`.
#' @param pattern Character vector or string specifying the geometric/visual pattern(s):
#'   * `"crop_rows"` (or `"field"`): Agricultural plantation rows with perspective convergence,
#'     foliage clumps along the ridges, and tilled furrows.
#'   * `"soil"` (or `"earth"`): Tilled soil with multiscale fractal roughness and directional shaded clods/crumbs.
#'   * `"straw"` (or `"mulch"`): Crop residues and dry stalks scattered across ground substrate (no-till farming).
#'   * `"pasture"` (or `"grass"`): Dense pasture turf with macro soil-fertility patches and anisotropic grass fibers.
#'   * `"forest"` (or `"canopy"`): Aerial forest and tree crowns with 3D dome shading, canopy gaps, and chlorophyll variations.
#'   * `"mountain"`: Mountain terrain profiles with atmospheric sky, distant misty ranges, and foreground rocky ridges.
#'   * `"sea"` (or `"water"`): Water and aquatic surfaces with multi-directional wave interference and specular sunlight highlights.
#'   * `"sky"` (or `"clouds"`): Atmospheric sky gradient with soft cumulus cloud billows.
#'   * `"lab"` (or `"paper"`, `"board"`, `"bench"`, `"germitest"`): Laboratory phenotyping benchtop, paper sheets,
#'     subtle millimeter grids, and studio lighting drop-off for seed/leaf scanning.
#'   * `"leaf_litter"` (or `"leaves"`, `"litter"`): Scattered fallen leaves and plant residue on soil.
#'   * `"mulch_film"` (or `"plastic"`): Agricultural black/silver plastic mulching bed with soil margins and specular wrinkles.
#'   * `"gravel"` (or `"concrete"`, `"floor"`): Greenhouse concrete or gravel substrate with fine aggregate micro-shadows.
#'   * `"composite"` (or `"field_mixed"`): Multi-layer natural field combining soil furrows, crop rows, straw residue, and weeds.
#'   * `"wood"`: Wood grain and timber growth rings with organic waviness.
#'   * `"waves"`: Multi-frequency sinusoidal/cosine color waves simulating fluid undulations.
#'   * `"gradient"`: Linear and radial color gradients with smooth illumination transitions.
#'   * `"plasma"`: Multi-octave turbulence noise producing cloudy, organic textures.
#'   * `"blobs"`: Random Gaussian color spots simulating out-of-focus background elements.
#'   * `"texture"`: Textured substrate with high-frequency chromatic variations.
#'   * `"random"` (default): Randomly picks among all available patterns (or the subset supplied in `pattern`).
#' @param colors Preset name (`"random"`, `"agriculture"`, `"soil"`, `"crop"`, `"straw"`,
#'   `"pasture"`, `"forest"`, `"mountain"`, `"sea"`, `"sky"`, `"wood"`, `"lab"`, `"germitest"`,
#'   `"leaf_litter"`, `"plastic_film"`, `"gravel"`, `"earth"`, `"foliage"`, `"neutral"`,
#'   `"vibrant"`, `"pastel"`, `"dark"`, `"light"`), a numeric RGB vector (e.g., `c(255, 250, 250)`
#'   in 0-255 or `c(1, 0.98, 0.98)` in 0-1 for a custom tone), a matrix or list of RGB values,
#'   a standard R color name (e.g. `"white"`, `"snow"`), or a character vector of hex colors
#'   (e.g. `c("#1c120c", "#5c3a21", "#90be6d")` or `"#fffafa"`). Defaults to `"random"`, which smartly pairs palettes to patterns.
#' @param color Alias for `colors` (convenient for passing a single RGB vector or color). Defaults to `NULL`.
#' @param noise Numeric. Standard deviation of subtle Gaussian sensor noise added to
#'   synthetic imagery (0 to 0.2). Defaults to `0.01`.
#' @param blur Numeric. Standard deviation (sigma) for Gaussian smoothing filter to soften edges (0 to 5). Defaults to `0`.
#' @param vignette Numeric. Strength of natural camera lens vignetting / radial illumination falloff (0 to 1). Defaults to `0.10`.
#' @param exposure Numeric. Strength of subtle photographic exposure jitter (+/- variation in brightness). Defaults to `0.05`.
#' @param train_prop Numeric. Proportion of images assigned to the training split (`images/train/`)
#'   versus validation (`images/val/`) when saving to a YOLO dataset. Defaults to `0.80`.
#' @param prefix Character. Filename prefix for saved images. Defaults to `"synth_bg_"`.
#' @param preview Logical. If `TRUE` and in an interactive session, displays a visual
#'   preview grid of the generated images. Defaults to `FALSE`.
#' @param as_job Logical. If `TRUE` and running inside RStudio, launches the synthetic
#'   background generation as a background job in RStudio, freeing up the console. Defaults to `FALSE`.
#' @param verbose Logical. Whether to print progress messages. Defaults to `TRUE`.
#' @return If `dir = NULL`, a list of pliman `Image` objects. If `dir` is specified,
#'   an invisible list with paths of generated training and validation images: `list(train = ..., val = ...)`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Generate 12 diverse agricultural, field, and laboratory backgrounds in memory
#' bgs <- yolo_background_synth(
#'   n = 12,
#'   pattern = c("crop_rows", "soil", "straw", "pasture", "forest", "lab", "leaf_litter", "composite"),
#'   preview = TRUE
#' )
#'
#' # 2. Generate warm white backgrounds using custom RGB vector c(255, 250, 250)
#' bgs_white <- yolo_background_synth(
#'   n = 5,
#'   color = c(255, 250, 250),
#'   pattern = "solid",
#'   preview = TRUE
#' )
#'
#' # 3. Generate and inject 20 synthetic background images into a YOLO dataset
#' yolo_background_synth(dir = "my_dataset", n = 20, pattern = "random")
yolo_background_synth <- function(dir = NULL,
                                  n = 10L,
                                  width = 640L,
                                  height = 640L,
                                  pattern = c("random", "crop_rows", "soil", "straw", "pasture",
                                              "forest", "mountain", "sea", "sky", "wood",
                                              "lab", "leaf_litter", "mulch_film", "gravel", "composite",
                                              "waves", "gradient", "plasma", "blobs", "texture", "solid"),
                                  colors = c("random", "agriculture", "soil", "crop", "straw",
                                             "pasture", "forest", "mountain", "sea", "sky", "wood",
                                             "lab", "germitest", "leaf_litter", "plastic_film", "gravel",
                                             "earth", "foliage", "neutral", "vibrant", "pastel", "dark", "light"),
                                  color = NULL,
                                  noise = 0.01,
                                  blur = 0,
                                  vignette = 0.10,
                                  exposure = 0.05,
                                  train_prop = 0.80,
                                  prefix = "synth_bg_",
                                  preview = FALSE,
                                  as_job = FALSE,
                                  verbose = TRUE) {
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (is.null(dir)) {
      cli::cli_abort("{.arg dir} must be specified when using {.code as_job = TRUE}.")
    }
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with background image synthesis in current session.")
    } else {
      mc <- match.call(expand.dots = TRUE)
      eval_env <- parent.frame()
      evaluated_args <- list()
      arg_names <- setdiff(names(mc), c("", "as_job"))
      for (arg_name in arg_names) {
        arg_val <- tryCatch(eval(mc[[arg_name]], envir = eval_env), error = function(e) mc[[arg_name]])
        evaluated_args[[arg_name]] <- arg_val
      }
      evaluated_args$as_job <- FALSE

      job_dir <- file.path(tempdir(), "yolo_jobs")
      if (!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("yolo_bg_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_bg_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (inherits(dir, "yolo_dataset")) basename(dir$dir)
        else if (is.character(dir) && nzchar(dir[1])) basename(dir[1])
        else "dataset"
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Synth Background: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_background_synth()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Synthetic Background - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_background_synth, job_info$args)\n",
        "cat('\\n[OK] Synthetic Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Synthetic Background Generation (pliman)")
      cli::cli_alert_success("Background image synthesis has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target directory: {.val {target_name}}",
        "Count: {.val {n}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }
  if (!is.null(color)) {
    colors <- color
  }
  n <- as.integer(n[1])
  if (n <= 0L) cli::cli_abort("{.arg n} must be at least 1.")
  width <- as.integer(width[1])
  height <- as.integer(height[1])
  if (width < 32L || height < 32L) cli::cli_abort("Dimensions must be at least 32x32 pixels.")
  noise <- max(0, min(0.5, as.numeric(noise[1])))
  blur <- max(0, min(10, as.numeric(blur[1])))
  vignette <- max(0, min(1, as.numeric(vignette[1])))
  exposure <- max(0, min(0.5, as.numeric(exposure[1])))
  palette_presets <- list(
    agriculture  = c("#1c120c", "#3d2314", "#5c3a21", "#2d5a27", "#4f772d", "#90be6d", "#d4a373", "#e9edc9", "#ffffff"),
    soil         = c("#1f140e", "#3e2723", "#5d4037", "#8d6e63", "#a0522d", "#b97a57", "#c85a17", "#d7ccc8", "#f5deb3"),
    crop         = c("#0b2512", "#1b4332", "#2d6a4f", "#40916c", "#52b788", "#74c69d", "#90be6d", "#b7e4c7", "#e9f5ed"),
    straw        = c("#2c1b0e", "#5c3a1d", "#8c5828", "#b8860b", "#d4a017", "#e6c35c", "#faedcd", "#fff8dc"),
    pasture      = c("#0d2818", "#1e4620", "#386641", "#6a994e", "#a7c957", "#d8f3dc", "#f2e8cf"),
    forest       = c("#051f15", "#113f2d", "#1b4332", "#2d6a4f", "#385e38", "#52796f", "#2b1d0c", "#84a98c"),
    mountain     = c("#0f172a", "#1e293b", "#334155", "#475569", "#64748b", "#94a3b8", "#cbd5e1", "#7dd3fc", "#ffffff"),
    sea          = c("#021a36", "#033b66", "#0077b6", "#0096c7", "#00b4d8", "#48cae4", "#90e0ef", "#caf0f8", "#ffffff"),
    sky          = c("#0284c7", "#38bdf8", "#7dd3fc", "#bae6fd", "#e0f2fe", "#f8fafc", "#ffffff"),
    wood         = c("#2b1d0c", "#4a3219", "#6e4720", "#8c5b29", "#b07d4f", "#d4a77a", "#e8cbab"),
    lab          = c("#cbd5e1", "#e2e8f0", "#edf2f7", "#f1f5f9", "#f8fafc", "#ffffff"),
    germitest    = c("#0c2d48", "#145da0", "#2e8bc0", "#84b6d4", "#b1d4e0", "#d4eaf7"),
    leaf_litter  = c("#2d1808", "#5c3317", "#8b4513", "#a0522d", "#cd853f", "#6b8e23", "#d2b48c"),
    plastic_film = c("#09090b", "#18181b", "#27272a", "#3f3f46", "#71717a", "#cbd5e1", "#f8fafc"),
    gravel       = c("#1f2937", "#374151", "#4b5563", "#6b7280", "#9ca3af", "#d1d5db"),
    earth        = c("#2c1810", "#8b5a2b", "#a0522d", "#d2b48c", "#f5deb3", "#e8d8c8", "#5c4033", "#3d2817"),
    foliage      = c("#1b4332", "#2d6a4f", "#40916c", "#52b788", "#74c69d", "#d8f3dc", "#90be6d", "#4d7c0f"),
    neutral      = c("#111827", "#1f2937", "#374151", "#4b5563", "#6b7280", "#9ca3af", "#d1d5db", "#f3f4f6", "#ffffff"),
    vibrant      = c("#ef4444", "#f97316", "#f59e0b", "#10b981", "#06b6d4", "#3b82f6", "#8b5cf6", "#ec4899"),
    pastel       = c("#fecdd3", "#fed7aa", "#fef08a", "#bbf7d0", "#a5f3fc", "#bfdbfe", "#ddd6fe", "#fce7f3"),
    dark         = c("#09090b", "#18181b", "#27272a", "#3f3f46", "#1e293b", "#0f172a", "#020617"),
    light        = c("#f8fafc", "#f1f5f9", "#e2e8f0", "#cbd5e1", "#ffffff", "#fafafa", "#f5f5f5")
  )
  # Color aliases
  palette_presets$harvest <- palette_presets$straw
  palette_presets$water   <- palette_presets$sea
  palette_presets$clouds  <- palette_presets$sky
  palette_presets$grass   <- palette_presets$pasture
  palette_presets$paper   <- palette_presets$lab
  pattern_aliases <- c(
    solid          = "solid",
    plain          = "solid",
    uniform        = "solid",
    flat           = "solid",
    mono           = "solid",
    field          = "crop_rows",
    rows           = "crop_rows",
    plantation     = "crop_rows",
    earth          = "soil",
    ground         = "soil",
    dirt           = "soil",
    mulch          = "straw",
    residue        = "straw",
    grass          = "pasture",
    meadow         = "pasture",
    lawn           = "pasture",
    canopy         = "forest",
    trees          = "forest",
    woods          = "forest",
    water          = "sea",
    ocean          = "sea",
    lake           = "sea",
    pond           = "sea",
    clouds         = "sky",
    cloud          = "sky",
    paper          = "lab",
    board          = "lab",
    bench          = "lab",
    germitest      = "lab",
    leaves         = "leaf_litter",
    litter         = "leaf_litter",
    foliage_litter = "leaf_litter",
    plastic        = "mulch_film",
    plastic_mulch  = "mulch_film",
    film           = "mulch_film",
    concrete       = "gravel",
    floor          = "gravel",
    road           = "gravel",
    field_mixed    = "composite",
    mixed          = "composite",
    natural        = "composite",
    bark           = "wood",
    timber         = "wood"
  )
  all_known_patterns <- c("crop_rows", "soil", "straw", "pasture", "forest",
                          "mountain", "sea", "sky", "wood", "lab", "leaf_litter",
                          "mulch_film", "gravel", "composite",
                          "waves", "gradient", "plasma", "blobs", "texture", "solid")
  # Resolve candidate patterns
  user_patterns <- tolower(pattern)
  user_patterns <- vapply(user_patterns, function(p) {
    if (p %in% names(pattern_aliases)) pattern_aliases[[p]] else p
  }, character(1))
  candidate_patterns <- if ("random" %in% user_patterns) {
    all_known_patterns
  } else {
    valid_sel <- intersect(user_patterns, all_known_patterns)
    if (length(valid_sel) == 0L) all_known_patterns else valid_sel
  }
  patt_to_palettes <- list(
    crop_rows   = c("agriculture", "crop", "soil", "straw", "pasture"),
    soil        = c("soil", "earth", "agriculture"),
    straw       = c("straw", "soil", "agriculture", "earth"),
    pasture     = c("pasture", "crop", "agriculture", "foliage"),
    forest      = c("forest", "foliage", "crop", "agriculture"),
    mountain    = c("mountain", "earth", "neutral", "sky"),
    sea         = c("sea", "sky", "light"),
    sky         = c("sky", "sea", "light", "pastel"),
    wood        = c("wood", "earth", "soil", "straw"),
    lab         = c("lab", "germitest", "light", "dark"),
    leaf_litter = c("leaf_litter", "earth", "soil", "foliage"),
    mulch_film  = c("plastic_film", "dark", "soil"),
    gravel      = c("gravel", "neutral", "soil", "earth"),
    composite   = c("agriculture", "soil", "straw", "crop"),
    solid       = c("light", "lab", "neutral", "white")
  )
  hex_to_rgb <- function(hex_vec) {
    rgb_mat <- grDevices::col2rgb(hex_vec) / 255
    lapply(seq_len(ncol(rgb_mat)), function(j) rgb_mat[, j])
  }
  generate_single_bg <- function(w, h, patt_arg, col_arg) {
    patt <- if (identical(patt_arg, "random") || !patt_arg %in% all_known_patterns) {
      sample(candidate_patterns, 1)
    } else {
      patt_arg
    }
    # Resolve palette
    rgb_list <- NULL
    if (is.numeric(col_arg) && length(col_arg) == 3L) {
      # Single RGB vector e.g. c(255, 250, 250) or c(1, 0.98, 0.98)
      rgb_norm <- if (max(col_arg, na.rm = TRUE) > 1) col_arg / 255 else col_arg
      rgb_norm <- pmax(0, pmin(1, as.numeric(rgb_norm)))
      if (patt %in% c("solid", "plain", "uniform", "flat")) {
        rgb_list <- list(rgb_norm, rgb_norm)
      } else {
        rgb_list <- list(
          pmax(0, rgb_norm * 0.93),
          pmax(0, rgb_norm * 0.97),
          rgb_norm,
          pmin(1, rgb_norm * 1.01)
        )
      }
    } else if (is.list(col_arg) && length(col_arg) >= 1L && all(vapply(col_arg, function(x) is.numeric(x) && length(x) == 3L, logical(1)))) {
      rgb_list <- lapply(col_arg, function(x) {
        v <- if (max(x, na.rm = TRUE) > 1) x / 255 else x
        pmax(0, pmin(1, as.numeric(v)))
      })
      if (length(rgb_list) == 1L) rgb_list <- c(rgb_list, rgb_list)
    } else if (is.numeric(col_arg) && is.matrix(col_arg) && (ncol(col_arg) == 3L || nrow(col_arg) == 3L)) {
      if (ncol(col_arg) == 3L) {
        rgb_list <- lapply(seq_len(nrow(col_arg)), function(r) {
          v <- col_arg[r, ]
          val <- if (max(v, na.rm = TRUE) > 1) v / 255 else v
          pmax(0, pmin(1, as.numeric(val)))
        })
      } else {
        rgb_list <- lapply(seq_len(ncol(col_arg)), function(c) {
          v <- col_arg[, c]
          val <- if (max(v, na.rm = TRUE) > 1) v / 255 else v
          pmax(0, pmin(1, as.numeric(val)))
        })
      }
      if (length(rgb_list) == 1L) rgb_list <- c(rgb_list, rgb_list)
    } else if (is.character(col_arg) && length(col_arg) == 1L && !col_arg %in% names(palette_presets) && !identical(col_arg, "random")) {
      rgb_try <- tryCatch(as.numeric(grDevices::col2rgb(col_arg) / 255), error = function(e) NULL)
      if (!is.null(rgb_try) && length(rgb_try) == 3L) {
        if (patt %in% c("solid", "plain", "uniform", "flat")) {
          rgb_list <- list(rgb_try, rgb_try)
        } else {
          rgb_list <- list(
            pmax(0, rgb_try * 0.93),
            pmax(0, rgb_try * 0.97),
            rgb_try,
            pmin(1, rgb_try * 1.01)
          )
        }
      }
    }

    if (is.null(rgb_list)) {
      pal_hex <- if (is.character(col_arg) && length(col_arg) == 1L && col_arg %in% names(palette_presets)) {
        palette_presets[[col_arg]]
      } else if (identical(col_arg, "random")) {
        chosen_pal_name <- if (patt %in% names(patt_to_palettes)) {
          sample(patt_to_palettes[[patt]], 1)
        } else {
          sample(names(palette_presets), 1)
        }
        palette_presets[[chosen_pal_name]]
      } else if (is.character(col_arg) && length(col_arg) > 1L && all(col_arg %in% names(palette_presets))) {
        palette_presets[[sample(col_arg, 1)]]
      } else if (is.character(col_arg) && all(grepl("^#", col_arg))) {
        col_arg
      } else {
        hex_try <- tryCatch({
          m <- grDevices::col2rgb(col_arg)
          grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
        }, error = function(e) NULL)
        if (!is.null(hex_try) && length(hex_try) >= 1L) {
          hex_try
        } else {
          palette_presets$agriculture
        }
      }
      if (length(pal_hex) < 2L) pal_hex <- c(pal_hex, "#ffffff")
      # Preserve monotonic luminance progression for ultra-natural shading
      k_cols <- min(length(pal_hex), sample(4:6, 1))
      idx <- sort(sample(seq_along(pal_hex), k_cols))
      sampled_hex <- pal_hex[idx]
      rgb_list <- hex_to_rgb(sampled_hex)
    }
    # In EBImage / pliman: dim 1 is width (x, cols), dim 2 is height (y, rows)
    x <- matrix(seq(-1, 1, length.out = w), nrow = w, ncol = h, byrow = FALSE)
    y <- matrix(seq(-1, 1, length.out = h), nrow = w, ncol = h, byrow = TRUE)
    if (patt == "crop_rows") {
      is_perspective <- stats::runif(1) > 0.4
      meander <- 0.035 * sin(2.5 * pi * y) + 0.015 * cos(5 * pi * x)
      if (is_perspective) {
        vx <- stats::runif(1, -0.3, 0.3)
        vy <- stats::runif(1, -1.8, -1.25)
        u <- (x + meander - vx) / (y - vy)
        f_rows <- stats::runif(1, 9, 15)
      } else {
        ang <- stats::runif(1, -pi/6, pi/6)
        u <- (x + meander) * cos(ang) + y * sin(ang)
        f_rows <- stats::runif(1, 7, 13)
      }
      row_prof <- 0.5 + 0.5 * cos(2 * pi * f_rows * u)
      clumps <- matrix(stats::rnorm(w * h, 0, 0.28), w, h)
      val <- row_prof + pmax(0, row_prof - 0.25)^2 * clumps + 0.08 * y
    } else if (patt == "soil") {
      val <- matrix(0, w, h)
      for (s in c(1, 2, 4, 8, 16)) {
        ang <- stats::runif(1, 0, 2 * pi)
        u <- x * cos(ang) + y * sin(ang)
        val <- val + (1 / sqrt(s)) * sin(s * pi * u + stats::runif(1, 0, 2 * pi))
      }
      val <- val * 0.4 + matrix(stats::rnorm(w * h, 0, 0.2), w, h)
      n_clods <- sample(45:80, 1)
      sun_ang <- stats::runif(1, 0, 2 * pi)
      sun_dx <- cos(sun_ang)
      sun_dy <- sin(sun_ang)
      for (k in seq_len(n_clods)) {
        cx <- stats::runif(1, -1.05, 1.05)
        cy <- stats::runif(1, -1.05, 1.05)
        cr <- stats::runif(1, 0.03, 0.12)
        c_min <- max(1L, as.integer(floor((cx - cr + 1) * 0.5 * (w - 1))) + 1L)
        c_max <- min(w, as.integer(ceiling((cx + cr + 1) * 0.5 * (w - 1))) + 1L)
        r_min <- max(1L, as.integer(floor((cy - cr + 1) * 0.5 * (h - 1))) + 1L)
        r_max <- min(h, as.integer(ceiling((cy + cr + 1) * 0.5 * (h - 1))) + 1L)
        if (c_min <= c_max && r_min <= r_max) {
          sub_x <- x[c_min:c_max, r_min:r_max, drop = FALSE]
          sub_y <- y[c_min:c_max, r_min:r_max, drop = FALSE]
          dist_sq <- (sub_x - cx)^2 + (sub_y - cy)^2
          mask <- dist_sq < (cr^2)
          if (any(mask)) {
            dome <- sqrt(pmax(0, 1 - dist_sq / (cr^2)))
            shade <- 1 + 0.45 * (((sub_x - cx) * sun_dx + (sub_y - cy) * sun_dy) / cr)
            val[c_min:c_max, r_min:r_max] <- val[c_min:c_max, r_min:r_max] + (dome * shade * mask) * 0.65
          }
        }
      }
    } else if (patt == "straw") {
      val <- matrix(stats::rnorm(w * h, 0, 0.08), w, h)
      for (s in c(1, 2, 4)) {
        ang <- stats::runif(1, 0, 2 * pi)
        val <- val + (0.15 / s) * sin(s * pi * (x * cos(ang) + y * sin(ang)))
      }
      n_straws <- sample(180:280, 1)
      for (s in seq_len(n_straws)) {
        sx <- stats::runif(1, -1.05, 1.05)
        sy <- stats::runif(1, -1.05, 1.05)
        s_ang <- stats::runif(1, 0, pi)
        s_len <- stats::runif(1, 0.06, 0.24)
        s_w <- stats::runif(1, 0.005, 0.014)
        s_curv <- stats::runif(1, -0.3, 0.3)
        rad_box <- s_len * 1.15
        c_min <- max(1L, as.integer(floor((sx - rad_box + 1) * 0.5 * (w - 1))) + 1L)
        c_max <- min(w, as.integer(ceiling((sx + rad_box + 1) * 0.5 * (w - 1))) + 1L)
        r_min <- max(1L, as.integer(floor((sy - rad_box + 1) * 0.5 * (h - 1))) + 1L)
        r_max <- min(h, as.integer(ceiling((sy + rad_box + 1) * 0.5 * (h - 1))) + 1L)
        if (c_min <= c_max && r_min <= r_max) {
          sub_x <- x[c_min:c_max, r_min:r_max, drop = FALSE]
          sub_y <- y[c_min:c_max, r_min:r_max, drop = FALSE]
          cos_a <- cos(s_ang)
          sin_a <- sin(s_ang)
          u_s <- (sub_x - sx) * cos_a + (sub_y - sy) * sin_a
          v_s <- -(sub_x - sx) * sin_a + (sub_y - sy) * cos_a - s_curv * (u_s^2)
          needle <- exp(-(v_s^2) / (2 * s_w^2)) * pmax(0, 1 - (u_s / s_len)^2)^2
          val[c_min:c_max, r_min:r_max] <- pmax(val[c_min:c_max, r_min:r_max], needle * stats::runif(1, 0.75, 1.35))
        }
      }
    } else if (patt == "pasture") {
      val_patch <- matrix(0, w, h)
      for (s in c(1, 2, 4)) {
        ang <- stats::runif(1, 0, 2 * pi)
        u <- x * cos(ang) + y * sin(ang)
        val_patch <- val_patch + (1 / s) * sin(s * pi * u + stats::runif(1, 0, 2 * pi))
      }
      ang_g1 <- stats::runif(1, -pi/7, pi/7)
      ang_g2 <- ang_g1 + stats::runif(1, pi/5, pi/3)
      u_g1 <- x * cos(ang_g1) + y * sin(ang_g1)
      u_g2 <- x * cos(ang_g2) + y * sin(ang_g2)
      f_grass1 <- stats::runif(1, 35, 60)
      f_grass2 <- stats::runif(1, 30, 50)
      blades1 <- 0.5 + 0.5 * sin(2 * pi * f_grass1 * u_g1 + matrix(stats::rnorm(w * h, 0, 0.4), w, h))
      blades2 <- 0.5 + 0.5 * sin(2 * pi * f_grass2 * u_g2 + matrix(stats::rnorm(w * h, 0, 0.4), w, h))
      val <- 0.35 * val_patch + 0.35 * blades1 + 0.2 * blades2 + matrix(stats::rnorm(w * h, 0, 0.12), w, h)
    } else if (patt == "forest") {
      n_trees <- sample(70:120, 1)
      val <- matrix(0, w, h)
      sun_ang <- stats::runif(1, 0, 2 * pi)
      sun_dx <- cos(sun_ang)
      sun_dy <- sin(sun_ang)
      for (t in seq_len(n_trees)) {
        tx <- stats::runif(1, -1.05, 1.05)
        ty <- stats::runif(1, -1.05, 1.05)
        tr <- stats::runif(1, 0.07, 0.22)
        c_min <- max(1L, as.integer(floor((tx - tr + 1) * 0.5 * (w - 1))) + 1L)
        c_max <- min(w, as.integer(ceiling((tx + tr + 1) * 0.5 * (w - 1))) + 1L)
        r_min <- max(1L, as.integer(floor((ty - tr + 1) * 0.5 * (h - 1))) + 1L)
        r_max <- min(h, as.integer(ceiling((ty + tr + 1) * 0.5 * (h - 1))) + 1L)
        if (c_min <= c_max && r_min <= r_max) {
          sub_x <- x[c_min:c_max, r_min:r_max, drop = FALSE]
          sub_y <- y[c_min:c_max, r_min:r_max, drop = FALSE]
          dist_sq <- (sub_x - tx)^2 + (sub_y - ty)^2
          mask <- dist_sq < (tr^2)
          if (any(mask)) {
            dome <- sqrt(pmax(0, 1 - dist_sq / (tr^2)))
            shade <- 1 + 0.35 * (((sub_x - tx) * sun_dx + (sub_y - ty) * sun_dy) / tr)
            bump <- dome * shade
            val[c_min:c_max, r_min:r_max] <- pmax(val[c_min:c_max, r_min:r_max], bump * mask)
          }
        }
      }
      val <- val + matrix(stats::rnorm(w * h, 0, 0.08), w, h) * (val > 0.05)
    } else if (patt == "mountain") {
      val <- (y + 1) * 0.25
      x_vec <- seq(-1, 1, length.out = w)
      ridge1 <- rep(-0.15, w)
      for (o in 1:4) {
        ridge1 <- ridge1 + (0.2 / o) * sin(o * 2.2 * x_vec + stats::runif(1, 0, 2*pi))
      }
      r1_mat <- matrix(ridge1, nrow = w, ncol = h, byrow = FALSE)
      mask1 <- y >= r1_mat
      val[mask1] <- 0.45 + 0.25 * (y[mask1] - r1_mat[mask1]) + matrix(stats::rnorm(sum(mask1), 0, 0.05), sum(mask1), 1)
      ridge2 <- rep(0.2, w)
      for (o in 1:5) {
        ridge2 <- ridge2 + (0.28 / (o^0.9)) * sin(o * 3.1 * x_vec + stats::runif(1, 0, 2*pi))
      }
      r2_mat <- matrix(ridge2, nrow = w, ncol = h, byrow = FALSE)
      mask2 <- y >= r2_mat
      slope_tex <- 0.6 + 0.3 * (y[mask2] - r2_mat[mask2]) + 0.15 * sin(16 * x[mask2] + 8 * y[mask2])
      val[mask2] <- slope_tex + matrix(stats::rnorm(sum(mask2), 0, 0.08), sum(mask2), 1)
    } else if (patt == "sea") {
      val <- matrix(0, w, h)
      for (w_i in 1:4) {
        ang_w <- stats::runif(1, -pi/3, pi/3)
        k_w <- stats::runif(1, 2.5, 7.5)
        u_w <- x * cos(ang_w) + y * sin(ang_w)
        val <- val + (1 / w_i) * (cos(k_w * pi * u_w + stats::runif(1, 0, 2*pi)) + 0.3 * (sin(2 * k_w * pi * u_w))^2)
      }
      sparkle <- pmax(0, val - 1.1)^3 * 2.5
      val <- val + sparkle + 0.2 * y
    } else if (patt == "sky") {
      val_sky <- matrix(0, w, h)
      for (s in c(1, 2, 4, 8)) {
        ang <- stats::runif(1, 0, 2 * pi)
        u <- x * cos(ang) + y * sin(ang)
        val_sky <- val_sky + (1 / s) * sin(s * pi * u + stats::runif(1, 0, 2 * pi))
      }
      clouds <- pmax(0, val_sky - 0.2)^1.4 * 2
      val <- (y + 1) * 0.4 + clouds
    } else if (patt == "lab") {
      spot_x <- stats::runif(1, -0.3, 0.3)
      spot_y <- stats::runif(1, -0.3, 0.3)
      lighting <- 1 - 0.25 * ((x - spot_x)^2 + (y - spot_y)^2)
      paper_fiber <- matrix(stats::rnorm(w * h, 0, 0.04), w, h)
      grid_mat <- matrix(0, w, h)
      if (stats::runif(1) > 0.4) {
        grid_spacing <- sample(c(20L, 30L, 40L), 1)
        c_idx <- seq(1, w, by = grid_spacing)
        r_idx <- seq(1, h, by = grid_spacing)
        grid_mat[c_idx, ] <- grid_mat[c_idx, ] - 0.12
        grid_mat[, r_idx] <- grid_mat[, r_idx] - 0.12
      }
      val <- lighting + paper_fiber + grid_mat
    } else if (patt == "leaf_litter") {
      val <- matrix(stats::rnorm(w * h, 0, 0.06), w, h)
      n_leaves <- sample(50:90, 1)
      for (l in seq_len(n_leaves)) {
        lx <- stats::runif(1, -1.05, 1.05)
        ly <- stats::runif(1, -1.05, 1.05)
        l_ang <- stats::runif(1, 0, 2 * pi)
        l_len <- stats::runif(1, 0.08, 0.22)
        l_w   <- l_len * stats::runif(1, 0.4, 0.65)
        rad_box <- l_len * 1.1
        c_min <- max(1L, as.integer(floor((lx - rad_box + 1) * 0.5 * (w - 1))) + 1L)
        c_max <- min(w, as.integer(ceiling((lx + rad_box + 1) * 0.5 * (w - 1))) + 1L)
        r_min <- max(1L, as.integer(floor((ly - rad_box + 1) * 0.5 * (h - 1))) + 1L)
        r_max <- min(h, as.integer(ceiling((ly + rad_box + 1) * 0.5 * (h - 1))) + 1L)
        if (c_min <= c_max && r_min <= r_max) {
          sub_x <- x[c_min:c_max, r_min:r_max, drop = FALSE]
          sub_y <- y[c_min:c_max, r_min:r_max, drop = FALSE]
          cos_a <- cos(l_ang)
          sin_a <- sin(l_ang)
          u_l <- (sub_x - lx) * cos_a + (sub_y - ly) * sin_a
          v_l <- -(sub_x - lx) * sin_a + (sub_y - ly) * cos_a
          dist_norm <- (u_l / l_len)^2 + (v_l / l_w)^2
          mask <- dist_norm < 1
          if (any(mask)) {
            leaf_dome <- sqrt(pmax(0, 1 - dist_norm))
            midrib <- 1 - 0.25 * exp(-(v_l^2) / (2 * (0.005)^2))
            leaf_val <- leaf_dome * midrib * stats::runif(1, 0.6, 1.1)
            val[c_min:c_max, r_min:r_max] <- pmax(val[c_min:c_max, r_min:r_max], leaf_val * mask)
          }
        }
      }
    } else if (patt == "mulch_film") {
      bed_w <- stats::runif(1, 0.55, 0.75)
      wrinkles <- 0.15 * sin(12 * pi * y + 3 * sin(2 * pi * x)) + 0.1 * cos(6 * pi * x)
      sheen <- exp(-((x - 0.1)^2) / 0.15) * 0.4
      plastic_val <- 0.45 + wrinkles + sheen
      soil_margin <- matrix(stats::rnorm(w * h, 0, 0.2), w, h)
      for (s in c(1, 2, 4)) {
        soil_margin <- soil_margin + (0.3 / s) * sin(s * pi * (x + y))
      }
      edge_trans <- pmin(1, pmax(0, (bed_w - abs(x)) / 0.08))
      val <- plastic_val * edge_trans + soil_margin * (1 - edge_trans)
    } else if (patt == "gravel") {
      val <- matrix(stats::rnorm(w * h, 0, 0.25), w, h)
      for (s in c(2, 4, 8, 16)) {
        ang <- stats::runif(1, 0, 2 * pi)
        val <- val + (0.2 / sqrt(s)) * sin(s * pi * (x * cos(ang) + y * sin(ang)))
      }
    } else if (patt == "composite") {
      val_soil <- matrix(0, w, h)
      for (s in c(1, 2, 4)) {
        ang <- stats::runif(1, 0, 2 * pi)
        val_soil <- val_soil + (0.3 / s) * sin(s * pi * (x * cos(ang) + y * sin(ang)))
      }
      f_rows <- stats::runif(1, 6, 10)
      ang_r <- stats::runif(1, -pi/8, pi/8)
      u_r <- x * cos(ang_r) + y * sin(ang_r)
      row_prof <- 0.5 + 0.5 * cos(2 * pi * f_rows * u_r)
      val_straw <- matrix(0, w, h)
      for (s in seq_len(sample(70:120, 1))) {
        sx <- stats::runif(1, -1.05, 1.05)
        sy <- stats::runif(1, -1.05, 1.05)
        s_ang <- stats::runif(1, 0, pi)
        s_len <- stats::runif(1, 0.05, 0.18)
        rad_box <- s_len * 1.1
        c_min <- max(1L, as.integer(floor((sx - rad_box + 1) * 0.5 * (w - 1))) + 1L)
        c_max <- min(w, as.integer(ceiling((sx + rad_box + 1) * 0.5 * (w - 1))) + 1L)
        r_min <- max(1L, as.integer(floor((sy - rad_box + 1) * 0.5 * (h - 1))) + 1L)
        r_max <- min(h, as.integer(ceiling((sy + rad_box + 1) * 0.5 * (h - 1))) + 1L)
        if (c_min <= c_max && r_min <= r_max) {
          sub_x <- x[c_min:c_max, r_min:r_max, drop = FALSE]
          sub_y <- y[c_min:c_max, r_min:r_max, drop = FALSE]
          u_s <- (sub_x - sx) * cos(s_ang) + (sub_y - sy) * sin(s_ang)
          v_s <- -(sub_x - sx) * sin(s_ang) + (sub_y - sy) * cos(s_ang)
          needle <- exp(-(v_s^2) / (2 * (0.008)^2)) * pmax(0, 1 - (u_s / s_len)^2)^2
          val_straw[c_min:c_max, r_min:r_max] <- pmax(val_straw[c_min:c_max, r_min:r_max], needle)
        }
      }
      val <- 0.45 * row_prof + 0.3 * val_soil + 0.35 * val_straw + matrix(stats::rnorm(w * h, 0, 0.1), w, h)
    } else if (patt == "wood") {
      ang_w <- stats::runif(1, -pi/12, pi/12)
      u_w <- x * cos(ang_w) + y * sin(ang_w)
      v_w <- -x * sin(ang_w) + y * cos(ang_w)
      f_wood <- stats::runif(1, 15, 30)
      waviness <- 0.15 * sin(3 * pi * v_w) + 0.08 * cos(6 * pi * u_w)
      grain <- sin(2 * pi * f_wood * (u_w + waviness) + matrix(stats::rnorm(w * h, 0, 0.2), w, h))
      val <- 0.6 * (0.5 + 0.5 * grain) + 0.3 * (v_w + 1) + matrix(stats::rnorm(w * h, 0, 0.1), w, h)
    } else if (patt == "waves") {
      n_h <- sample(2:5, 1)
      val <- matrix(0, w, h)
      for (k in seq_len(n_h)) {
        ang <- stats::runif(1, 0, 2 * pi)
        freq <- stats::runif(1, 0.5, 4.5)
        phase <- stats::runif(1, 0, 2 * pi)
        u <- x * cos(ang) + y * sin(ang)
        val <- val + (1 / k) * sin(2 * pi * freq * u + phase)
      }
    } else if (patt == "gradient") {
      ang <- stats::runif(1, 0, 2 * pi)
      u <- x * cos(ang) + y * sin(ang)
      radial <- sqrt(x^2 + y^2) * stats::runif(1, 0.2, 0.9)
      val <- u + radial
    } else if (patt == "plasma") {
      val <- matrix(0, w, h)
      for (s in c(1, 2, 4, 8)) {
        ang <- stats::runif(1, 0, 2 * pi)
        u <- x * cos(ang) + y * sin(ang)
        val <- val + (1 / s) * sin(s * pi * u + stats::runif(1, 0, 2 * pi))
      }
    } else if (patt == "blobs") {
      val <- matrix(0, w, h)
      n_blobs <- sample(15:35, 1)
      for (b in seq_len(n_blobs)) {
        bx <- stats::runif(1, -1.1, 1.1)
        by <- stats::runif(1, -1.1, 1.1)
        br <- stats::runif(1, 0.15, 0.55)
        val <- val + exp(-((x - bx)^2 + (y - by)^2) / (2 * br^2))
      }
    } else if (patt %in% c("solid", "plain", "uniform", "flat")) {
      val <- matrix(0.5, w, h) + matrix(stats::rnorm(w * h, 0, 0.002), w, h)
    } else { # texture
      val <- matrix(stats::rnorm(w * h), w, h)
      ang <- stats::runif(1, 0, 2 * pi)
      val <- val * 0.35 + (x * cos(ang) + y * sin(ang)) * 0.65
    }
    rng <- range(val)
    norm_val <- if (diff(rng) > 0) (val - rng[1]) / diff(rng) else matrix(0.5, w, h)
    r_mat <- matrix(0, w, h)
    g_mat <- matrix(0, w, h)
    b_mat <- matrix(0, w, h)
    n_segs <- length(rgb_list) - 1L
    for (i in seq_len(n_segs)) {
      c1 <- rgb_list[[i]]
      c2 <- rgb_list[[i + 1L]]
      t1 <- (i - 1L) / n_segs
      t2 <- i / n_segs
      mask <- (norm_val >= t1 & norm_val <= t2)
      if (any(mask)) {
        local_t <- (norm_val[mask] - t1) / (t2 - t1)
        r_mat[mask] <- c1[1] * (1 - local_t) + c2[1] * local_t
        g_mat[mask] <- c1[2] * (1 - local_t) + c2[2] * local_t
        b_mat[mask] <- c1[3] * (1 - local_t) + c2[3] * local_t
      }
    }
    if (noise > 0) {
      r_mat <- r_mat + matrix(stats::rnorm(w * h, 0, noise), w, h)
      g_mat <- g_mat + matrix(stats::rnorm(w * h, 0, noise), w, h)
      b_mat <- b_mat + matrix(stats::rnorm(w * h, 0, noise), w, h)
    }
    if (vignette > 0) {
      # Natural camera lens vignetting falloff
      v_falloff <- 1 - vignette * (((x^2 + y^2) / 2))
      r_mat <- r_mat * v_falloff
      g_mat <- g_mat * v_falloff
      b_mat <- b_mat * v_falloff
    }
    if (exposure > 0) {
      # Subtle exposure / illumination jitter per image
      exp_factor <- stats::runif(1, 1 - exposure, 1 + exposure)
      r_mat <- r_mat * exp_factor
      g_mat <- g_mat * exp_factor
      b_mat <- b_mat * exp_factor
    }
    r_mat <- pmax(0, pmin(1, r_mat))
    g_mat <- pmax(0, pmin(1, g_mat))
    b_mat <- pmax(0, pmin(1, b_mat))
    arr <- array(c(r_mat, g_mat, b_mat), dim = c(w, h, 3))
    out_im <- as_image(arr, colormode = "Color")
    if (blur > 0) {
      out_im <- image_blur(out_im, sigma = blur)
    }
    out_im
  }
  color_desc <- if (is.numeric(colors) && length(colors) == 3L) {
    sprintf("RGB(%s)", paste(round(colors, 1), collapse = ", "))
  } else if (is.list(colors) || (is.numeric(colors) && is.matrix(colors))) {
    "Custom RGB Palette"
  } else if (is.character(colors) && length(colors) > 1L) {
    paste(colors[1:min(3L, length(colors))], collapse = ", ")
  } else {
    as.character(colors[1])
  }
  if (isTRUE(verbose)) {
    cli::cli_h2("Generating {n} Synthetic YOLO Background Image{?s}")
    cli::cli_bullets(c(
      "*" = "Dimensions: {width}x{height} pixels",
      "*" = "Pattern{?s}: {.val {candidate_patterns}}",
      "*" = "Colors: {.val {color_desc}}",
      "*" = "Noise: {noise} (blur: {blur})"
    ))
    pb <- cli::cli_progress_bar(
      name = "Synthesizing background images",
      total = n,
      format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
    )
  }
  valid_preset_names <- names(palette_presets)
  is_preset_selection <- is.character(colors) &&
    all(tolower(colors) %in% c("random", valid_preset_names)) &&
    !any(grepl("^#", colors))

  candidate_color_presets <- if (is_preset_selection) {
    p_sel <- intersect(tolower(colors), valid_preset_names)
    if (length(p_sel) == 0L || "random" %in% tolower(colors)) valid_preset_names else p_sel
  } else NULL

  is_custom_col <- !is_preset_selection && (
    is.numeric(colors) || is.list(colors) ||
    (is.character(colors) && (all(grepl("^#", colors)) || tryCatch(!is.null(grDevices::col2rgb(colors)), error = function(e) FALSE)))
  )

  uids <- tolower(gsub("-", "", uuid(n)))
  images_list <- vector("list", n)
  for (i in seq_len(n)) {
    p_curr <- sample(candidate_patterns, 1)
    c_curr <- if (is_preset_selection) {
      sample(candidate_color_presets, 1)
    } else if (is_custom_col) {
      colors
    } else if (identical(colors[1], "random")) {
      "random"
    } else {
      colors[1]
    }
    images_list[[i]] <- generate_single_bg(width, height, p_curr, c_curr)
    names(images_list)[i] <- sprintf("%s%s", prefix, uids[i])
    if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
  }
  if (isTRUE(verbose)) cli::cli_progress_done(id = pb)
  # Preview if requested
  if (isTRUE(preview)) {
    tryCatch({
      image_combine(images_list, ncol = min(5L, ceiling(sqrt(n))))
    }, error = function(e) {
      if (length(images_list) > 0L) plot(images_list[[1]])
    })
  }
  # If no directory supplied, return images in memory
  if (is.null(dir)) {
    return(images_list)
  }
  # Resolve directory (YOLO dataset or output folder)
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
  is_yolo_ds <- dir.exists(file.path(dir, "images")) ||
    file.exists(file.path(dir, "data.yaml")) ||
    dir.exists(file.path(dir, "train", "images"))

  if (dir.exists(file.path(dir, "train", "images"))) {
    train_img_dir <- file.path(dir, "train", "images")
    val_img_dir   <- if (dir.exists(file.path(dir, "val", "images"))) file.path(dir, "val", "images") else file.path(dir, "train", "images")
    train_lbl_dir <- file.path(dir, "train", "labels")
    val_lbl_dir   <- if (dir.exists(file.path(dir, "val", "labels"))) file.path(dir, "val", "labels") else file.path(dir, "train", "labels")
  } else {
    train_img_dir <- file.path(dir, "images", "train")
    val_img_dir   <- file.path(dir, "images", "val")
    train_lbl_dir <- file.path(dir, "labels", "train")
    val_lbl_dir   <- file.path(dir, "labels", "val")
  }

  if (is_yolo_ds || dir.exists(file.path(dir, "images", "train")) || dir.exists(file.path(dir, "train", "images"))) {
    for (d in c(train_img_dir, val_img_dir, train_lbl_dir, val_lbl_dir)) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
    n_train <- max(1L, round(n * train_prop))
    if (n_train >= n && n > 1L) n_train <- n - 1L
    n_val <- n - n_train
  } else {
    # Plain directory export
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    n_train <- n
    n_val <- 0L
    train_img_dir <- dir
  }
  added_train <- character()
  added_val <- character()
  for (k in seq_along(images_list)) {
    is_tr <- (k <= n_train)
    dest_dir <- if (is_tr) train_img_dir else val_img_dir
    fn_base <- names(images_list)[k]
    if (is.null(fn_base) || !nzchar(fn_base)) {
      fn_base <- sprintf("%s%s", prefix, uids[k])
    }
    dest_img <- file.path(dest_dir, paste0(fn_base, ".jpg"))
    image_export(images_list[[k]], name = paste0(fn_base, ".jpg"), subfolder = dest_dir, verbose = FALSE)
    if (is_yolo_ds || dir.exists(file.path(dir, "labels"))) {
      dest_lbl_dir <- if (is_tr) train_lbl_dir else val_lbl_dir
      dest_lbl <- file.path(dest_lbl_dir, paste0(fn_base, ".txt"))
      file.create(dest_lbl)
    }
    if (is_tr) added_train <- c(added_train, dest_img) else added_val <- c(added_val, dest_img)
  }
  # Invalidate caches
  if (dir.exists(file.path(dir, "labels"))) {
    unlink(file.path(dir, "labels", "train.cache"), force = TRUE)
    unlink(file.path(dir, "labels", "val.cache"), force = TRUE)
  }
  if (dir.exists(file.path(dir, "train", "labels"))) {
    unlink(file.path(dir, "train", "labels.cache"), force = TRUE)
    unlink(file.path(dir, "val", "labels.cache"), force = TRUE)
  }
  if (isTRUE(verbose)) {
    cli::cli_alert_success("Saved {n} synthetic background images to {.path {dir}}.")
    if (length(added_val) > 0L) {
      cli::cli_bullets(c(
        "*" = "Training background images: {length(added_train)}",
        "*" = "Validation background images: {length(added_val)}",
        "v" = "Empty label .txt files generated for all background samples."
      ))
    }
  }
  invisible(list(train = added_train, val = added_val))
}
#' @export
yolo_synth_background <- yolo_background_synth
#' @export
yolo_dataset_synth_background <- yolo_background_synth
#' @title Remove Synthetic Background Images and Labels from a YOLO Dataset
#' @name yolo_background_synth_clean
#' @aliases yolo_synth_background_clean yolo_dataset_synth_clean yolo_dataset_clean_synth yolo_clean_synth_background yolo_dataset_background_clean
#' @description
#' `yolo_background_synth_clean()` removes all synthetic background images and matching empty label files
#' (matching the pattern `^synth_bg_` or a custom `prefix`) previously generated by [yolo_background_synth()]
#' or [yolo_dataset_add_background()], reverting the dataset back to its pristine state.
#'
#' It automatically invalidates stale YOLO cache files (`train.cache`, `val.cache`, `test.cache`)
#' to ensure subsequent YOLO training or validation does not reference deleted files.
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param prefix Character. Filename prefix of the synthetic background files to remove (default: `"synth_bg_"`).
#'   If set to `"all"`, it cleans both synthetic (`"^synth_bg_"`) and generic (`"^bg_"`) background files.
#' @param pattern Optional character regular expression. If provided, overrides `prefix` to match files.
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset, or an invisible list with removal counts.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Inject synthetic backgrounds into a YOLO dataset
#' yolo_background_synth(dir = "my_dataset", n = 15, pattern = "waves")
#'
#' # 2. Clean up synthetic backgrounds and return to original state
#' yolo_background_synth_clean("my_dataset")
yolo_background_synth_clean <- function(data,
                                        prefix = "synth_bg_",
                                        pattern = NULL,
                                        verbose = TRUE) {
  dir_path <- tryCatch(.yolo_get_dir(data), error = function(e) {
    if (is.list(data)) {
      if (!is.null(data$dir) && is.character(data$dir) && nzchar(data$dir[1])) {
        data$dir[1]
      } else if (!is.null(data$yaml_file) && is.character(data$yaml_file) && nzchar(data$yaml_file[1])) {
        dirname(data$yaml_file[1])
      } else if (length(data) > 0L && is.character(data[[1]]) && nzchar(data[[1]])) {
        data[[1]]
      } else {
        "."
      }
    } else if (is.character(data) && length(data) >= 1L && nzchar(data[1])) {
      data[1]
    } else {
      "."
    }
  })
  raw_dir <- dir_path
  dir_path <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    candidates <- c(
      file.path(getwd(), raw_dir),
      file.path(tempdir(), raw_dir),
      file.path(pliman_model_dir(), raw_dir)
    )
    for (cand in candidates) {
      if (dir.exists(cand)) {
        dir_path <- normalizePath(cand, winslash = "/", mustWork = FALSE)
        break
      }
    }
  }
  task <- tryCatch(.yolo_detect_task(dir_path), error = function(e) "detect")
  classes <- tryCatch(.yolo_get_classes(dir_path), error = function(e) character())
  # Resolve regex pattern for matching files
  patt_regex <- if (!is.null(pattern) && is.character(pattern) && nzchar(pattern[1])) {
    pattern[1]
  } else if (identical(tolower(prefix[1]), "all")) {
    "^(synth_bg_|bg_)"
  } else {
    paste0("^", prefix[1])
  }
  img_ext_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  n_removed_img <- 0L
  n_removed_lbl <- 0L
  n_orig_img <- 0L
  # Search images
  img_dir <- file.path(dir_path, "images")
  search_img_dir <- if (dir.exists(img_dir)) img_dir else dir_path
  all_imgs <- list.files(search_img_dir, pattern = img_ext_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  synth_imgs <- all_imgs[grepl(patt_regex, basename(all_imgs), ignore.case = TRUE)]
  if (length(synth_imgs) > 0L) {
    unlink(synth_imgs, force = TRUE)
    n_removed_img <- length(synth_imgs)
  }
  n_orig_img <- length(all_imgs) - n_removed_img
  # Search labels
  lbl_dir <- file.path(dir_path, "labels")
  if (dir.exists(lbl_dir)) {
    all_lbls <- list.files(lbl_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    synth_lbls <- all_lbls[grepl(patt_regex, basename(all_lbls), ignore.case = TRUE)]
    if (length(synth_lbls) > 0L) {
      unlink(synth_lbls, force = TRUE)
      n_removed_lbl <- length(synth_lbls)
    }
    # Invalidate stale cache files
    unlink(file.path(lbl_dir, "train.cache"), force = TRUE)
    unlink(file.path(lbl_dir, "val.cache"), force = TRUE)
    unlink(file.path(lbl_dir, "test.cache"), force = TRUE)
  }
  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No synthetic background files matching {.val {patt_regex}} found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned synthetic background files from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "*" = "Synthetic background images removed: {n_removed_img}",
        "*" = "Background label files removed: {n_removed_lbl}",
        "*" = "Remaining images in dataset: {n_orig_img}",
        "v" = "YOLO label caches invalidated."
      ))
    }
  }
  res <- tryCatch(.yolo_wrap_dataset(dir_path, task = task, classes = classes), error = function(e) {
    list(dir = dir_path, images_removed = n_removed_img, labels_removed = n_removed_lbl, remaining_images = n_orig_img)
  })
  invisible(res)
}
#' @export
yolo_synth_background_clean <- yolo_background_synth_clean
#' @export
yolo_dataset_synth_clean <- yolo_background_synth_clean
#' @export
yolo_dataset_clean_synth <- yolo_background_synth_clean
#' @export
yolo_clean_synth_background <- yolo_background_synth_clean
#' @export
yolo_dataset_background_clean <- yolo_background_synth_clean
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
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Organize flower photos into train (70%), val (20%), test (10%)
#' ds <- yolo_dataset_classify(
#'   src_dir = "flower_photos",
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
#' @examplesIf interactive()
#' # Download and export the entire dataset permanently before URLs expire:
#' ds <- yolo_dataset_from_ndjson(
#'   file = "african-wildlife.ndjson",
#'   zip = TRUE
#' )
#' yolo_dataset_preview(ds)
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
    file.path(pliman_model_dir(), pt_file)
  )
  for (cand in unique(cands)) {
    if (file.exists(cand) && file.info(cand)$size > 100000) {
      return(normalizePath(cand, winslash = "/"))
    }
  }
  dest_file <- file.path(output_dir, pt_file)
  # Download directly from Ultralytics Hugging Face repository for YOLO26
  if (grepl("yolo26", pt_file, ignore.case = TRUE)) {
    hf_ultra_url <- paste0("https://huggingface.co/Ultralytics/YOLO26/resolve/main/", pt_file)
    cli::cli_alert_info("Pre-trained weights {.val {pt_file}} not found locally.")
    cli::cli_alert_info("Downloading from Ultralytics Hugging Face repository ({.url {hf_ultra_url}})...")
    dl_success <- .download_url_robust(hf_ultra_url, dest_file, min_size = 100000)
  } else {
    nepem_base <- pliman_models_base_url()
    nepem_url <- paste0(nepem_base, pt_file)
    cli::cli_alert_info("Pre-trained weights {.val {pt_file}} not found locally.")
    cli::cli_alert_info("Downloading pre-trained weights ({.url {nepem_url}})...")
    dl_success <- .download_url_robust(nepem_url, dest_file, min_size = 100000)
  }
  if (!dl_success) {
    cli::cli_alert_warning("Primary download failed. Retrying with Ultralytics GitHub release assets fallback...")
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
  dest_norm <- normalizePath(dest_file, winslash = "/")
  cli::cli_alert_success("Model weights {.file {pt_file}} successfully downloaded to {.path {dest_norm}}")
  return(dest_norm)
}
# ==============================================================================
# YOLO TRAINING & EXPORT PIPELINE
# ==============================================================================
# Helper to recursively copy directory contents safely
.copy_dir_contents <- function(src_dir, dest_dir) {
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  norm_src <- normalizePath(src_dir, winslash = "/", mustWork = FALSE)
  norm_dest <- normalizePath(dest_dir, winslash = "/", mustWork = FALSE)
  if (norm_src == norm_dest) return(TRUE)

  items <- list.files(src_dir, full.names = FALSE, all.files = FALSE, include.dirs = TRUE)
  for (item in items) {
    src_item <- file.path(src_dir, item)
    dest_item <- file.path(dest_dir, item)
    if (dir.exists(src_item)) {
      if (!dir.exists(dest_item)) {
        dir.create(dest_item, recursive = TRUE, showWarnings = FALSE)
      }
      .copy_dir_contents(src_item, dest_item)
    } else {
      file.copy(src_item, dest_item, overwrite = TRUE)
    }
  }
  TRUE
}

#' @title Export YOLO Model to ONNX and Pliman Model Zoo
#' @name yolo_export
#' @description
#' Exports a trained YOLO PyTorch (`.pt`) model to standalone ONNX (`.onnx`) format, or
#' imports/registers an existing `.onnx` model and its training results into the pliman model
#' repository ([pliman_model_dir()]).
#'
#' When `export_results = TRUE` (default), all diagnostic plots, loss curves, confusion matrices,
#' and performance metrics from the training run are also copied to `<output_name>_results` in
#' `output_dir`. This allows inspecting the exported model anytime with [yolo_results()] or
#' [yolo_model_details()].
#'
#' @param model A model checkpoint path (`.pt` or `.onnx`), a training results directory,
#'   a results `.zip` archive, or a `yolo_train` object from [yolo_results()] or [yolo_train()].
#' @param output_dir Target directory to save the exported model and results. Defaults to [pliman_model_dir()].
#' @param output_name Desired filename or identifier for the exported model (e.g. `"capsulas_det"` or `"capsulas_det.onnx"`).
#'   If `NULL` (default), it is inferred automatically from the model name or results folder.
#' @param export_results Logical. If `TRUE` (default), copies the full training results directory
#'   (curves, metrics, confusion matrices, weights) alongside the model so that [yolo_results()] can load them.
#' @param imgsz Integer image resolution (width and height) for ONNX export. Default is `640`.
#' @param dynamic Logical. If `TRUE` (default), exports the ONNX model with dynamic input and output
#'   axes (`[batch, 3, height, width]`), enabling inference with variable batch sizes and flexible
#'   image resolutions (e.g., 640, 960, 1024, 1280) in downstream functions like [image_detect_dl()]
#'   and [video_detect_dl()]. If `FALSE`, inputs are fixed to `imgsz x imgsz`.
#' @param simplify Logical. If `TRUE` (default), runs graph optimization via `onnxslim` / `onnx-simplifier`
#'   to prune redundant operators, fuse layers, and reduce model file size and runtime latency.
#' @param opset Integer specifying the ONNX operator set version (default is `17L`). Modern YOLO architectures
#'   typically recommend opset 17 or 18 for compatibility with advanced operators and DirectML/CUDA.
#' @param half Logical. If `TRUE`, converts model weights and graph activations to FP16 half-precision,
#'   reducing file size by approximately 50\% and accelerating GPU/NPU inference on hardware with native FP16 support.
#'   Defaults to `FALSE` (standard FP32 single precision).
#' @param auto_install Logical. If `TRUE` (default), automatically installs missing Python export libraries
#'   (`ultralytics`, `onnx`, `onnxslim`).
#'
#' @return The full path to the exported `.onnx` model (invisibly).
#' @export
#' @seealso [yolo_train()], [yolo_results()], [yolo_model_details()], [image_detect_dl()], [video_detect_dl()], [pliman_model_dir()]
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Export from a PyTorch checkpoint with dynamic axes
#' yolo_export("best.pt", output_name = "custom_model.onnx", imgsz = 640, dynamic = TRUE)
#'
#' # 2. Export with FP16 half-precision for GPU optimization
#' yolo_export("best.pt", output_name = "custom_fp16.onnx", half = TRUE)
#'
#' # 3. Register an external training folder with all results
#' yolo_export("path/to/model_capsulas", output_name = "capsulas_det")
#'
#' # 4. Directly export from a yolo_results object
#' res <- yolo_results("path/to/model_capsulas")
#' yolo_export(res, output_name = "capsulas_det")
#'
#' # 5. Inspect exported model details anytime
#' yolo_model_details("capsulas_det.onnx")
yolo_export <- function(model,
                        output_dir = pliman_model_dir(),
                        output_name = NULL,
                        export_results = TRUE,
                        imgsz = 640,
                        dynamic = TRUE,
                        simplify = TRUE,
                        opset = 17L,
                        half = FALSE,
                        auto_install = TRUE) {
  orig_results_dir <- NULL
  if (is.list(model)) {
    orig_results_dir <- model$results_dir %||% NULL
    cand_model <- model$model_file %||% model$best_onnx %||% model$best_pt %||% model$model
    if (!is.null(cand_model) && is.character(cand_model) && file.exists(cand_model[1])) {
      model <- cand_model[1]
    }
  }

  if (is.null(model) || !nzchar(as.character(model)[1])) {
    cli::cli_abort("Please specify a valid model checkpoint, results directory, or yolo_train object.")
  }

  model_raw <- as.character(model)[1]

  # Check if model is a directory containing weights/ or results
  cand_results_dir <- orig_results_dir
  if (dir.exists(model_raw)) {
    if (is.null(cand_results_dir)) cand_results_dir <- model_raw
    # Search for model weights inside directory
    inner_cands <- c(
      file.path(model_raw, "weights", "best.onnx"),
      file.path(model_raw, "weights", "best.pt"),
      file.path(model_raw, "best.onnx"),
      file.path(model_raw, "best.pt"),
      list.files(file.path(model_raw, "weights"), pattern = "\\.onnx$", full.names = TRUE),
      list.files(file.path(model_raw, "weights"), pattern = "\\.pt$", full.names = TRUE),
      list.files(model_raw, pattern = "\\.onnx$", full.names = TRUE),
      list.files(model_raw, pattern = "\\.pt$", full.names = TRUE)
    )
    inner_cands <- inner_cands[file.exists(inner_cands) & !dir.exists(inner_cands)]
    if (length(inner_cands) > 0L) {
      ox_cands <- grep("\\.onnx$", inner_cands, value = TRUE)
      model_raw <- if (length(ox_cands) > 0L) ox_cands[1] else inner_cands[1]
    } else {
      cli::cli_abort("No .onnx or .pt model weights found inside directory {.path {model_raw}}.")
    }
  }

  if (!file.exists(model_raw)) {
    cli::cli_abort("Model checkpoint file not found at {.path {model_raw}}.")
  }

  model <- normalizePath(model_raw, winslash = "/")
  output_dir <- pliman_model_dir(output_dir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Check for surrounding results directory if not already known
  if (is.null(cand_results_dir) || !dir.exists(cand_results_dir)) {
    p1 <- dirname(model)
    p2 <- dirname(p1)
    is_res_dir <- function(d) {
      if (!dir.exists(d)) return(FALSE)
      file.exists(file.path(d, "args.yaml")) ||
      file.exists(file.path(d, "results.csv")) ||
      dir.exists(file.path(d, "weights"))
    }
    if (basename(p1) == "weights" && is_res_dir(p2)) {
      cand_results_dir <- p2
    } else if (is_res_dir(p1)) {
      cand_results_dir <- p1
    } else {
      resolved <- .resolve_yolo_results_dir(model = model)
      if (!is.null(resolved) && dir.exists(resolved)) {
        cand_results_dir <- resolved
      }
    }
  }

  # Infer base_stem
  if (!is.null(output_name) && nzchar(output_name)) {
    base_stem <- sub("\\.(onnx|pt)$", "", basename(output_name), ignore.case = TRUE)
  } else {
    raw_stem <- tools::file_path_sans_ext(basename(model))
    if (raw_stem %in% c("best", "last") && !is.null(cand_results_dir)) {
      res_stem <- sub("_results$", "", basename(cand_results_dir), ignore.case = TRUE)
      if (!res_stem %in% c("weights", "train", "val", "")) {
        raw_stem <- res_stem
      }
    }
    base_stem <- raw_stem
  }

  dest_onnx <- file.path(output_dir, paste0(base_stem, ".onnx"))
  dest_pt <- file.path(output_dir, paste0(base_stem, ".pt"))
  dest_results_dir <- file.path(output_dir, paste0(base_stem, "_results"))

  # Internal helper to sync results directory
  .sync_results <- function() {
    if (isTRUE(export_results) && !is.null(cand_results_dir) && dir.exists(cand_results_dir)) {
      .copy_dir_contents(cand_results_dir, dest_results_dir)
      dest_w_dir <- file.path(dest_results_dir, "weights")
      if (!dir.exists(dest_w_dir)) dir.create(dest_w_dir, recursive = TRUE, showWarnings = FALSE)
      if (file.exists(dest_onnx)) {
        file.copy(dest_onnx, file.path(dest_w_dir, "best.onnx"), overwrite = TRUE)
      }
      # Copy best.pt if present
      cand_pt <- file.path(cand_results_dir, "weights", "best.pt")
      if (!file.exists(cand_pt)) cand_pt <- file.path(cand_results_dir, "best.pt")
      if (!file.exists(cand_pt)) cand_pt <- file.path(dirname(model), "best.pt")
      if (file.exists(cand_pt)) {
        file.copy(cand_pt, file.path(dest_w_dir, "best.pt"), overwrite = TRUE)
        file.copy(cand_pt, dest_pt, overwrite = TRUE)
      }
      cli::cli_alert_success("Training results and diagnostics exported to: {.file {dest_results_dir}}")
      cli::cli_alert_info("You can inspect the results anytime with: {.code yolo_results('{base_stem}')}")
    }
  }

  # If model is already an ONNX file, simply copy it directly to destination (no Python needed)
  if (grepl("\\.onnx$", model, ignore.case = TRUE)) {
    norm_dest <- normalizePath(dest_onnx, winslash = "/", mustWork = FALSE)
    if (model != norm_dest) {
      file.copy(model, dest_onnx, overwrite = TRUE)
      cli::cli_alert_success("Model saved to: {.file {dest_onnx}}")
    } else {
      cli::cli_alert_success("Model is already in ONNX format at: {.file {dest_onnx}}")
    }
    .sync_results()
    return(invisible(dest_onnx))
  }

  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)
  run_dir <- file.path(tempdir(), "yolo_export_run")
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  py_export_script <- file.path(run_dir, "export.py")
  py_dyn_str <- if (isTRUE(dynamic)) "True" else "False"
  py_simp_str <- if (isTRUE(simplify)) "True" else "False"
  py_half_str <- if (isTRUE(half)) "True" else "False"
  py_code <- paste0(
    "import sys\n",
    "from ultralytics import YOLO\n",
    "print('--> Loading weights from: ", model, "')\n",
    "m = YOLO('", model, "')\n",
    "try:\n",
    "    out = m.export(format='onnx', imgsz=", as.integer(imgsz),
    ", dynamic=", py_dyn_str, ", simplify=", py_simp_str,
    ", opset=", as.integer(opset), ", half=", py_half_str, ")\n",
    "except Exception as e:\n",
    "    print(f'Standard export notice: {e}. Retrying without half/simplify...')\n",
    "    out = m.export(format='onnx', imgsz=", as.integer(imgsz),
    ", dynamic=", py_dyn_str, ", opset=", as.integer(opset), ")\n",
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

  # Also copy .pt model if present
  if (file.exists(model)) {
    file.copy(model, dest_pt, overwrite = TRUE)
  }

  tryCatch(clear_onnx_sessions_cpp(), error = function(e) NULL)

  .sync_results()
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
#' @param data Path to the dataset directory (e.g. `"yolo_dataset"` or `"pasta"`), or to a dataset
#'   YAML configuration file (e.g. `"pasta/data.yaml"` or `"pasta/dataset.yaml"`). Can also be a
#'   `yolo_dataset` object returned by [yolo_dataset_export()], [yolo_dataset_autolabel()], or
#'   [yolo_dataset_tiling()]. Defaults to `"yolo_dataset"`.
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
#' @param resume Logical (`TRUE` to automatically resume from the latest checkpoint) or
#'   character path pointing to an existing `last.pt` checkpoint file. Defaults to `FALSE`.
#' @param checkpoint Logical. If `TRUE` (default), continuously saves real-time checkpoints:
#'   mirrors `last.pt` and `results.csv` after each epoch, and copies `best.pt` directly to
#'   `results_dir` and `output_dir` immediately whenever a new `*(best)` epoch is reached.
#'   This guarantees no progress is lost even if training is interrupted. Defaults to `TRUE`.
#' @param checkpoint_onnx Logical. If `TRUE`, also generates and updates a checkpoint `.onnx`
#'   model in real-time whenever a new `*(best)` epoch is achieved. Defaults to `FALSE`.
#' @param as_job Logical. If `TRUE` and running inside RStudio, launches the training
#'   as an asynchronous RStudio Background Job, immediately freeing up the R console while
#'   live progress and epoch tables stream in RStudio's "Jobs" tab. Checkpoints continue
#'   to be saved in real time, and the resulting model object is automatically saved to
#'   `model_res.rds` inside `results_dir`. Defaults to `FALSE`.
#' @param patience Early stopping patience in epochs. Defaults to `15`.
#' @param lr0 Initial learning rate (default `0.01`).
#' @param lrf Final learning rate fraction (default `0.01`).
#' @param optimizer Optimizer to use (`"auto"`, `"SGD"`, `"Adam"`, `"AdamW"`, `"RMSProp"`). Default is `"auto"`.
#' @param close_mosaic Disable mosaic augmentation in the last N epochs. Default is `10`.
#' @param amp Automatic Mixed Precision (AMP). Set to `FALSE` to disable AMP if GPU checks fail. Default is `NULL` (auto-detected).
#' @param plots Logical. If `FALSE` (default), suppresses generation of low-resolution matplotlib raster image files from YOLO; native vector plots are generated directly in R from exported CSV data.
#' @param verbose Logical. Whether to display informative CLI messages during training and export. Defaults to `TRUE`.
#' @param progress Optional logical or `NULL`. Whether to display a live progress indicator during training (defaults to `NULL`, following interactive/verbose defaults).
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
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Train a custom coffee grain segmentation model with real-time checkpointing
#' model_res <- yolo_train(
#'   data = "dataset_cafe", # or "dataset_cafe/data.yaml"
#'   model = "yolo26n-seg",
#'   epochs = 30,
#'   output_name = "cafe_seg.onnx"
#' )
#'
#' # 2. Resume training seamlessly from the last saved checkpoint if interrupted:
#' yolo_train(data = "dataset_cafe", resume = TRUE, epochs = 50)
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
yolo_train <- function(data = "yolo_dataset",
                       model = "yolo26n-seg",
                       epochs = 50,
                       imgsz = 640,
                       batch = 16,
                       device = "auto",
                       output_dir = pliman_model_dir(),
                       output_name = NULL,
                       results_dir = NULL,
                       export_pt = TRUE,
                       resume = FALSE,
                       checkpoint = TRUE,
                       checkpoint_onnx = FALSE,
                       as_job = FALSE,
                       patience = 15,
                       lr0 = 0.01,
                       lrf = 0.01,
                       optimizer = "auto",
                       close_mosaic = 10,
                       amp = NULL,
                       plots = FALSE,
                       verbose = TRUE,
                       progress = NULL,
                       workers = if (.Platform$OS.type == "windows") 0L else 4L,
                       auto_install = TRUE,
                       ...) {
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with training in current session.")
    } else {
      mc <- match.call(expand.dots = TRUE)
      eval_env <- parent.frame()
      evaluated_args <- list()
      arg_names <- setdiff(names(mc), c("", "as_job"))
      for (arg_name in arg_names) {
        arg_val <- tryCatch(eval(mc[[arg_name]], envir = eval_env), error = function(e) mc[[arg_name]])
        evaluated_args[[arg_name]] <- arg_val
      }
      evaluated_args$as_job <- FALSE
      job_dir <- file.path(tempdir(), "yolo_jobs")
      if (!dir.exists(job_dir)) {
        dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
      }
      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("yolo_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_args_", time_stamp, ".rds"))
      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)
      target_name <- if (!is.null(output_name) && nzchar(output_name)) {
        output_name
      } else if (is.character(data) && nzchar(data[1])) {
        basename(data[1])
      } else {
        "model"
      }
      j_name <- paste0("YOLO Train: ", target_name)
      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_train()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Training - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_train, job_info$args)\n",
        "if (!is.null(res$results_dir) && dir.exists(res$results_dir)) {\n",
        "  save_path <- file.path(res$results_dir, 'model_res.rds')\n",
        "  saveRDS(res, save_path)\n",
        "  cat('\\n[OK] Model results successfully saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)
      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )
      cli::cli_h2("YOLO Model Training (pliman)")
      cli::cli_alert_success("Training has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Real-time checkpoints: Enabled (continuous sync to results_dir)"
      ))
      cli::cli_alert_info("Your R console is free! You can track live epochs, losses and mAP in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }
  # 0. Robust resolution of dataset input (directory, YAML file, or yolo_dataset object)
  if (is.list(data)) {
    if (!is.null(data$data_yaml) && is.character(data$data_yaml) && nzchar(data$data_yaml[1]) && file.exists(data$data_yaml[1])) {
      data <- data$data_yaml[1]
    } else if (!is.null(data$yaml_file) && is.character(data$yaml_file) && nzchar(data$yaml_file[1]) && file.exists(data$yaml_file[1])) {
      data <- data$yaml_file[1]
    } else if (!is.null(data$dir) && is.character(data$dir) && nzchar(data$dir[1])) {
      data <- data$dir[1]
    } else if (!is.null(data$dataset) && is.list(data$dataset)) {
      if (!is.null(data$dataset$data_yaml) && file.exists(data$dataset$data_yaml[1])) {
        data <- data$dataset$data_yaml[1]
      } else if (!is.null(data$dataset$dir)) {
        data <- data$dataset$dir[1]
      }
    } else if (length(data) > 0L && is.character(data[[1]]) && nzchar(data[[1]])) {
      data <- data[[1]]
    }
  }
  raw_data <- as.character(data)[1]
  # Check relative to working directory or tempdir if needed
  if (!file.exists(raw_data) && !dir.exists(raw_data)) {
    for (cand in c(file.path(getwd(), raw_data), file.path(tempdir(), raw_data))) {
      if (file.exists(cand) || dir.exists(cand)) {
        raw_data <- cand
        break
      }
    }
  }
  # Also check if raw_data was e.g. "pasta" and "pasta/data.yaml" or "pasta/dataset.yaml" exists
  if (!file.exists(raw_data) && !dir.exists(raw_data)) {
    for (sub_yaml in c("data.yaml", "dataset.yaml")) {
      cand_sub <- file.path(raw_data, sub_yaml)
      if (file.exists(cand_sub)) {
        raw_data <- cand_sub
        break
      }
      for (base_cand in c(getwd(), tempdir())) {
        cand_sub2 <- file.path(base_cand, raw_data, sub_yaml)
        if (file.exists(cand_sub2)) {
          raw_data <- cand_sub2
          break
        }
      }
      if (file.exists(raw_data)) break
    }
  }
  data_norm <- normalizePath(raw_data, winslash = "/", mustWork = FALSE)
  if (!file.exists(data_norm) && !dir.exists(data_norm)) {
    cli::cli_abort("Could not find dataset configuration file or directory at {.path {data}}.")
  }
  # Resolve training progress display mode: 'batch' (Ultralytics live progress bars), 'epoch' (table), 'none' (silent)
  if (!is.null(progress)) {
    if (isTRUE(progress) || identical(progress, "batch")) {
      verbose_mode <- "batch"
    } else if (identical(progress, "epoch") || identical(progress, "table")) {
      verbose_mode <- "epoch"
    } else if (isFALSE(progress) || identical(progress, "none") || identical(progress, "silent")) {
      verbose_mode <- "none"
    } else {
      verbose_mode <- "batch"
    }
  } else {
    if (isTRUE(verbose) || identical(verbose, "batch")) {
      verbose_mode <- "batch"
    } else if (identical(verbose, "epoch") || identical(verbose, "table")) {
      verbose_mode <- "epoch"
    } else if (isFALSE(verbose) || identical(verbose, "none") || identical(verbose, "silent")) {
      verbose_mode <- "none"
    } else {
      verbose_mode <- "batch"
    }
  }
  # Split into data_dir (directory) and yaml_file (configuration file)
  if (dir.exists(data_norm)) {
    data_dir <- data_norm
    yaml_candidates <- c(
      file.path(data_dir, "data.yaml"),
      file.path(data_dir, "dataset.yaml")
    )
    existing_yaml <- yaml_candidates[file.exists(yaml_candidates)]
    if (length(existing_yaml) == 0L) {
      all_yamls <- list.files(data_dir, pattern = "\\.ya?ml$", full.names = TRUE)
      if (length(all_yamls) > 0L) existing_yaml <- all_yamls[1]
    }
    yaml_file <- if (length(existing_yaml) > 0L) existing_yaml[1] else NULL
  } else {
    yaml_file <- data_norm
    data_dir <- dirname(data_norm)
  }
  output_dir <- pliman_model_dir(output_dir)
  # 1. Ensure Python runtime & hardware environment
  py_exec <- .ensure_python(auto_install = auto_install)
  env_info <- .setup_python_yolo_env(py_exec, auto_install = auto_install)
  # 2. Inspect dataset structure to prevent task mismatch
  lbl_dir <- if (dir.exists(file.path(data_dir, "labels"))) {
    file.path(data_dir, "labels")
  } else if (dir.exists(file.path(dirname(data_dir), "labels"))) {
    file.path(dirname(data_dir), "labels")
  } else {
    file.path(data_dir, "labels")
  }
  lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  dataset_task <- "unknown"
  if (length(lbl_files) > 0L) {
    for (lf in lbl_files[seq_len(min(15L, length(lbl_files)))]) {
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
    train_dir <- if (dir.exists(file.path(data_dir, "train"))) file.path(data_dir, "train") else file.path(dirname(data_dir), "train")
    val_dir   <- if (dir.exists(file.path(data_dir, "val"))) file.path(data_dir, "val") else file.path(dirname(data_dir), "val")
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
      "The dataset at {.path {data_dir}} is structured for image classification ({.val classify}), but model {.val {model}} was requested."
    )
    cli::cli_alert_info(
      "Automatically adjusting base model to classification ({.val {cls_model}})."
    )
    model <- cls_model
    if (missing(imgsz) || imgsz == 640) imgsz <- 224
  } else if (identical(dataset_task, "detect") && grepl("seg", model, ignore.case = TRUE)) {
    det_model <- sub("[-_]seg", "", model, ignore.case = TRUE)
    cli::cli_alert_warning(
      "The dataset at {.path {data_dir}} contains bounding boxes ({.val detect}), but a segmentation model ({.val {model}}) was requested."
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
  # Ensure YAML file exists and has valid, synchronized paths for detection and segmentation tasks
  if (!identical(dataset_task, "classify") && !grepl("cls", model, ignore.case = TRUE)) {
    yaml_file <- .yolo_sync_data_yaml(data_dir, yaml_file = yaml_file, verbose = isTRUE(verbose))
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
  data_folder_name <- basename(data_dir)
  if (!nzchar(data_folder_name) || data_folder_name == ".") {
    data_folder_name <- basename(normalizePath(data_dir, winslash = "/", mustWork = FALSE))
  }
  if (is.null(output_name) || !nzchar(output_name)) {
    prefix <- if (nzchar(data_folder_name) && data_folder_name != ".") data_folder_name else "custom_yolo"
    suffix <- if (is_seg) "_seg" else if (is_obb) "_obb" else if (is_pose) "_pose" else if (is_cls) "_cls" else "_det"
    output_name <- paste0(prefix, suffix)
  }
  base_stem <- sub("\\.(onnx|onxx|pt)$", "", output_name, ignore.case = TRUE)
  dest_onnx <- file.path(output_dir, paste0(base_stem, ".onnx"))
  dest_pt <- file.path(output_dir, paste0(base_stem, ".pt"))
  if (is.null(results_dir) || !nzchar(results_dir)) {
    results_dir <- file.path(output_dir, paste0(base_stem, "_results"))
  }
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  }
  resume_pt_path <- NULL
  if (isTRUE(resume)) {
    candidates <- c(
      file.path(results_dir, "weights", "last.pt"),
      file.path(output_dir, paste0(base_stem, "_results"), "weights", "last.pt")
    )
    temp_candidates <- tryCatch({
      list.files(dirname(tempdir()), pattern = "^last\\.pt$", recursive = TRUE, full.names = TRUE)
    }, error = function(e) character(0))
    if (length(temp_candidates) > 0) {
      temp_candidates <- temp_candidates[grepl("train_run[\\\\/]weights[\\\\/]last\\.pt$", temp_candidates)]
      if (length(temp_candidates) > 0) {
        mtimes <- file.info(temp_candidates)$mtime
        temp_candidates <- temp_candidates[order(mtimes, decreasing = TRUE)]
        candidates <- c(candidates, temp_candidates)
      }
    }
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) {
      resume_pt_path <- normalizePath(existing[1], winslash = "/", mustWork = TRUE)
      cli::cli_alert_success("Resuming training from checkpoint: {.path {resume_pt_path}}")
    } else {
      cli::cli_abort("No checkpoint ('last.pt') found to resume from in {.path {results_dir}} or temporary session folders.")
    }
  } else if (is.character(resume) && nzchar(resume)) {
    if (!file.exists(resume)) {
      cli::cli_abort("Specified checkpoint file does not exist: {.path {resume}}")
    }
    resume_pt_path <- normalizePath(resume, winslash = "/", mustWork = TRUE)
    cli::cli_alert_success("Resuming training from checkpoint: {.path {resume_pt_path}}")
  }
  # 4. Display training configuration banner
  hw_desc <- if (isTRUE(env_info$cuda_available)) {
    paste0("NVIDIA GPU (", env_info$gpu_info$name, ") \U0001f680")
  } else {
    "Multi-threaded CPU \U0001f4bb"
  }
  dataset_display <- if (!identical(dataset_task, "classify") && !is.null(yaml_file) && file.exists(yaml_file)) yaml_file else data_dir
  banner_items <- c(
    "Dataset: {.path {dataset_display}}",
    "Base model: {.val {basename(base_model)}}",
    "Hardware: {hw_desc}",
    "Epochs: {epochs} (patience: {patience})",
    "Batch size: {batch}",
    "Resolution: {imgsz}x{imgsz} px",
    "Results folder: {.path {results_dir}}",
    paste0("Progress display: ", if (verbose_mode == "batch") "Live batch progress bar & metrics (it/s, ETA)" else if (verbose_mode == "epoch") "Compact epoch summary table" else "Silent")
  )
  if (!is.null(resume_pt_path)) {
    banner_items <- c(banner_items, paste0("Resume: Resuming from ", resume_pt_path))
  }
  if (isTRUE(checkpoint)) {
    banner_items <- c(banner_items, "Checkpoints: Real-time mirror to results_dir (weights/last.pt & best.pt)")
  }
  cli::cli_h2("YOLO Model Training (pliman)")
  cli::cli_ul(banner_items)
  # 5. Build Python training & export script
  run_dir <- file.path(tempdir(), "yolo_train_run")
  if (!dir.exists(run_dir)) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  }
  py_script_path <- file.path(run_dir, "train_and_export.py")
  # Argument passed to Ultralytics model.train()
  train_data_arg <- if (identical(dataset_task, "classify") || is_cls) {
    data_dir
  } else {
    yaml_file
  }
  # Hyperparameters and training options
  train_args <- list(
    data = train_data_arg,
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
  if (!is.null(resume_pt_path)) {
    train_args$resume <- TRUE
  }
  dots <- list(...)
  for (nm in names(dots)) {
    if (nzchar(nm) && !nm %in% c("verbose", "progress")) {
      train_args[[nm]] <- dots[[nm]]
    }
  }
  kwargs_py_lines <- vapply(names(train_args), function(k) {
    paste0("        '", k, "': ", .r_to_py_val(train_args[[k]]))
  }, character(1))
  kwargs_py_code <- paste(kwargs_py_lines, collapse = ",\n")
  base_model_norm <- normalizePath(base_model, winslash = "/", mustWork = FALSE)
  run_dir_norm <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  results_dir_norm <- normalizePath(results_dir, winslash = "/", mustWork = FALSE)
  dest_pt_norm <- normalizePath(dest_pt, winslash = "/", mustWork = FALSE)
  dest_onnx_norm <- normalizePath(dest_onnx, winslash = "/", mustWork = FALSE)
  resume_pt_norm <- if (!is.null(resume_pt_path)) normalizePath(resume_pt_path, winslash = "/", mustWork = FALSE) else ""
  checkpoint_py <- if (isTRUE(checkpoint)) "True" else "False"
  checkpoint_onnx_py <- if (isTRUE(checkpoint_onnx)) "True" else "False"
  is_seg_py <- if (is_seg) "True" else "False"
  is_cls_py <- if (is_cls) "True" else "False"
  dev_str <- as.character(device)
  verbose_mode_py <- sprintf("'%s'", verbose_mode)
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
    "verbose_mode = ", verbose_mode_py, "\n",
    "if verbose_mode != 'batch':\n",
    "    ultralytics.utils.TQDM = SilentTQDM\n",
    "    ultralytics.engine.trainer.TQDM = SilentTQDM\n",
    "    ultralytics.engine.validator.TQDM = SilentTQDM\n",
    "    logging.getLogger('ultralytics').setLevel(logging.WARNING)\n",
    "    try:\n",
    "        ultralytics.utils.LOGGER.setLevel(logging.WARNING)\n",
    "    except Exception:\n",
    "        pass\n",
    "else:\n",
    "    try:\n",
    "        ultralytics.utils.LOGGER.setLevel(logging.INFO)\n",
    "    except Exception:\n",
    "        pass\n\n",
    "from ultralytics import YOLO, settings\n",
    "try:\n",
    "    settings.update({'runs_dir': r'", run_dir_norm, "'})\n",
    "except Exception:\n",
    "    pass\n\n",
    "best_fitness = -1.0\n",
    "best_epoch = 1\n",
    "best_map50 = 0.0\n",
    "t_start = None\n",
    "ep_start = None\n",
    "is_seg = ", is_seg_py, "\n",
    "is_cls = ", is_cls_py, "\n",
    "results_dir = r'", results_dir_norm, "'\n",
    "dest_pt = r'", dest_pt_norm, "'\n",
    "dest_onnx = r'", dest_onnx_norm, "'\n",
    "checkpoint_enabled = ", checkpoint_py, "\n",
    "checkpoint_onnx = ", checkpoint_onnx_py, "\n",
    "resume_pt = r'", resume_pt_norm, "'\n\n",
    "import os, shutil\n\n",
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
    "    if verbose_mode == 'epoch':\n",
    "        print('')\n",
    "        if is_cls:\n",
    "            print('  Epoch    GPU Mem       Loss     Top1-Acc     Top5-Acc          Time')\n",
    "            print('  ' + '-' * 70)\n",
    "        elif is_seg:\n",
    "            print('  Epoch    GPU Mem     Box Loss     Seg Loss     Cls Loss       mAP50     mAP50-95          Time')\n",
    "            print('  ' + '-' * 88)\n",
    "        else:\n",
    "            print('  Epoch    GPU Mem     Box Loss     Cls Loss       mAP50     mAP50-95          Time')\n",
    "            print('  ' + '-' * 76)\n",
    "        sys.stdout.flush()\n\n",
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
    "    if verbose_mode == 'epoch':\n",
    "        if is_cls:\n",
    "            line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}   {c_loss:>8.4f}     {top1:>7.3f}      {top5:>7.3f}    {ep_t:>10}{tag}'\n",
    "        elif is_seg:\n",
    "            line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}     {b_loss:>8.4f}     {s_loss:>8.4f}     {c_loss:>8.4f}     {m50:>7.3f}     {m95:>7.3f}    {ep_t:>10}{tag}'\n",
    "        else:\n",
    "            line = f'  {ep:>3}/{total:<3}  {gpu_mem:>8}     {b_loss:>8.4f}     {c_loss:>8.4f}     {m50:>7.3f}     {m95:>7.3f}    {ep_t:>10}{tag}'\n",
    "        print(line)\n",
    "        sys.stdout.flush()\n",
    "    elif verbose_mode == 'batch' and is_best:\n",
    "        if is_cls:\n",
    "            print(f'  --> Epoch {ep}/{total}: new best model saved! (Top1-Acc: {best_fitness:.3f}){tag}')\n",
    "        else:\n",
    "            print(f'  --> Epoch {ep}/{total}: new best model saved! (mAP50: {best_map50:.3f}, mAP50-95: {m95:.3f}){tag}')\n",
    "        sys.stdout.flush()\n",
    "    if checkpoint_enabled:\n",
    "        try:\n",
    "            w_dir = os.path.join(trainer.save_dir, 'weights')\n",
    "            cw_dir = os.path.join(results_dir, 'weights')\n",
    "            os.makedirs(cw_dir, exist_ok=True)\n",
    "            f_last = os.path.join(w_dir, 'last.pt')\n",
    "            if os.path.isfile(f_last):\n",
    "                shutil.copy2(f_last, os.path.join(cw_dir, 'last.pt'))\n",
    "            f_csv = os.path.join(trainer.save_dir, 'results.csv')\n",
    "            if os.path.isfile(f_csv):\n",
    "                shutil.copy2(f_csv, os.path.join(results_dir, 'results.csv'))\n",
    "            if is_best:\n",
    "                f_best = os.path.join(w_dir, 'best.pt')\n",
    "                if os.path.isfile(f_best):\n",
    "                    shutil.copy2(f_best, os.path.join(cw_dir, 'best.pt'))\n",
    "                    try:\n",
    "                        if dest_pt:\n",
    "                            os.makedirs(os.path.dirname(dest_pt), exist_ok=True)\n",
    "                            shutil.copy2(f_best, dest_pt)\n",
    "                    except Exception:\n",
    "                        pass\n",
    "                    if checkpoint_onnx:\n",
    "                        try:\n",
    "                            ckpt_m = YOLO(f_best)\n",
    "                            ckpt_m.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=True)\n",
    "                            src_ox = os.path.splitext(f_best)[0] + '.onnx'\n",
    "                            if os.path.isfile(src_ox):\n",
    "                                shutil.copy2(src_ox, os.path.join(cw_dir, 'best.onnx'))\n",
    "                                if dest_onnx:\n",
    "                                    shutil.copy2(src_ox, dest_onnx)\n",
    "                        except Exception:\n",
    "                            pass\n",
    "        except Exception:\n",
    "            pass\n\n",
    "def on_train_end(trainer):\n",
    "    if verbose_mode != 'none':\n",
    "        elapsed_str = format_time(time.time() - t_start)\n",
    "        if verbose_mode == 'epoch':\n",
    "            print('  ' + '-' * (70 if is_cls else (88 if is_seg else 76)))\n",
    "        if is_cls:\n",
    "            print(f'  [OK] Training finished in {elapsed_str} (Best Accuracy Top-1: {best_fitness:.3f} at epoch {best_epoch})')\n",
    "        else:\n",
    "            print(f'  [OK] Training finished in {elapsed_str} (Best mAP50: {best_map50:.3f} at epoch {best_epoch})')\n",
    "        print('')\n",
    "        sys.stdout.flush()\n\n",
    "def main():\n",
    "    dev = '", dev_str, "'\n",
    "    if dev == 'auto':\n",
    "        dev = 0 if torch.cuda.is_available() else 'cpu'\n",
    "    elif dev.isdigit():\n",
    "        dev = int(dev)\n\n",
    "    if resume_pt and os.path.isfile(resume_pt):\n",
    "        model = YOLO(resume_pt)\n",
    "    else:\n",
    "        model = YOLO(r'", base_model_norm, "')\n",
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
    "    train_kwargs['exist_ok'] = True\n",
    "    train_kwargs['verbose'] = (verbose_mode == 'batch')\n\n",
    "    results = model.train(**train_kwargs)\n\n",
    "    best_pt = r'", run_dir_norm, "/train_run/weights/best.pt'\n",
    "    best_model = YOLO(best_pt)\n",
    "    try:\n",
    "        exported_path = best_model.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=True)\n",
    "    except Exception as e:\n",
    "        print(f'Standard export notice: {e}. Retrying with opset=17...')\n",
    "        exported_path = best_model.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=True, opset=17)\n\n",
    "    import csv, os, numpy as np\n",
    "    from pathlib import Path\n",
    "    target_dir = r'", run_dir_norm, "/train_run'\n",
    "    try:\n",
    "        val_res = best_model.val(project=r'", run_dir_norm, "', name='val_run', exist_ok=True, plots=True, workers=0)\n",
    "        import shutil\n",
    "        val_dir = os.path.join(r'", run_dir_norm, "', 'val_run')\n",
    "        if os.path.exists(val_dir):\n",
    "            for fname in os.listdir(val_dir):\n",
    "                s = os.path.join(val_dir, fname)\n",
    "                d = os.path.join(target_dir, fname)\n",
    "                if os.path.isfile(s) and not os.path.exists(d):\n",
    "                    try:\n",
    "                        shutil.copy2(s, d)\n",
    "                    except Exception:\n",
    "                        pass\n",
    "        cm = getattr(val_res, 'confusion_matrix', None)\n",
    "        if cm is not None and hasattr(cm, 'matrix') and cm.matrix is not None:\n",
    "            # Ultralytics stores cm.matrix as [predicted_class, true_class].\n",
    "            # Transpose to standard [true_class, predicted_class] where rows = ground truth and cols = predicted.\n",
    "            matrix = np.array(cm.matrix).T\n",
    "            n_rows, n_cols = matrix.shape\n",
    "            raw_names = [best_model.names[i] for i in range(len(best_model.names))]\n",
    "            col_names = (raw_names + ['background']) if n_cols > len(raw_names) else raw_names[:n_cols]\n",
    "            row_names = (raw_names + ['background']) if n_rows > len(raw_names) else raw_names[:n_rows]\n",
    "            cm_csv = os.path.join(target_dir, 'confusion_matrix.csv')\n",
    "            with open(cm_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                writer = csv.writer(f)\n",
    "                writer.writerow(['ground_truth'] + col_names)\n",
    "                for i, row in enumerate(matrix):\n",
    "                    r_lbl = row_names[i] if i < len(row_names) else f'class_{i}'\n",
    "                    writer.writerow([r_lbl] + [int(val) for val in row])\n",
    "            norm_matrix = np.zeros_like(matrix, dtype=float)\n",
    "            has_bg = len(row_names) > len(raw_names)\n",
    "            for r_idx in range(n_rows):\n",
    "                if has_bg and r_idx == n_rows - 1:\n",
    "                    c_sums = matrix.sum(axis=0)\n",
    "                    for c_idx in range(n_cols):\n",
    "                        norm_matrix[r_idx, c_idx] = round(float(matrix[r_idx, c_idx] / c_sums[c_idx]), 4) if (c_sums[c_idx] > 0 and c_idx < n_cols - 1) else 0.0\n",
    "                else:\n",
    "                    r_sum = float(matrix[r_idx].sum())\n",
    "                    if r_sum > 0:\n",
    "                        norm_matrix[r_idx] = np.round(matrix[r_idx] / r_sum, 4)\n",
    "            norm_csv = os.path.join(target_dir, 'confusion_matrix_normalized.csv')\n",
    "            with open(norm_csv, 'w', newline='', encoding='utf-8') as f:\n",
    "                writer = csv.writer(f)\n",
    "                writer.writerow(['ground_truth'] + col_names)\n",
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
  if (verbose_mode != "none") cli::cli_alert_info("Training progress:")
  start_time <- Sys.time()
  train_status <- system2(py_exec, args = shQuote(py_script_path))
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  if (train_status != 0L) {
    # Salvage any checkpoints generated before failure
    train_run_dir <- file.path(run_dir, "train_run")
    if (dir.exists(train_run_dir)) {
      contents <- list.files(train_run_dir, full.names = TRUE)
      for (item in contents) {
        tryCatch(file.copy(item, results_dir, recursive = TRUE, overwrite = TRUE), error = function(e) NULL)
      }
    }
    salvaged_last <- file.path(results_dir, "weights", "last.pt")
    salvaged_best <- file.path(results_dir, "weights", "best.pt")
    if (file.exists(salvaged_best)) {
      tryCatch(file.copy(salvaged_best, dest_pt, overwrite = TRUE), error = function(e) NULL)
    }
    if (file.exists(salvaged_last)) {
      cli::cli_alert_warning("Training process exited with status {train_status}, but checkpoints were successfully salvaged to {.path {results_dir}}!")
      cli::cli_alert_info("To resume training from the last saved epoch, run:")
      cli::cli_code(paste0("yolo_train(..., resume = \"", normalizePath(salvaged_last, winslash = "/", mustWork = FALSE), "\")"))
    }
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
  val_run_dir <- file.path(run_dir, "val_run")
  if (dir.exists(val_run_dir)) {
    contents <- list.files(val_run_dir, full.names = TRUE)
    for (item in contents) {
      file.copy(item, results_dir, recursive = TRUE, overwrite = TRUE)
    }
  }
  # Ensure no unwanted 'runs' directory remains in working directory or output_dir
  wd_runs <- file.path(getwd(), "runs")
  if (dir.exists(wd_runs)) {
    unlink(wd_runs, recursive = TRUE, force = TRUE)
  }
  if (!is.null(output_dir)) {
    out_runs <- file.path(output_dir, "runs")
    if (dir.exists(out_runs)) {
      unlink(out_runs, recursive = TRUE, force = TRUE)
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
  aligned_cm <- .align_yolo_confusion_matrix(cm_counts = cm_df, cm_norm = cm_norm_df)
  cm_df <- aligned_cm$cm_counts
  cm_norm_df <- aligned_cm$cm_norm
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
      dataset_dir = data_dir,
      data_yaml = yaml_file,
      call = match.call()
    ),
    class = c("yolo_train", "list")
  )
  invisible(res)
}
#' Internal helper to locate a YOLO subrun directory
#' @noRd
.find_yolo_subrun_dir <- function(dir) {
  if (is.null(dir) || !dir.exists(dir)) return(NULL)
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)

  # 1. Check if results.csv or weights/ or plots exist directly in dir
  has_csv <- file.exists(file.path(dir, "results.csv"))
  has_weights <- dir.exists(file.path(dir, "weights"))
  has_plots <- any(file.exists(file.path(dir, c("confusion_matrix.png", "results.png", "PR_curve.png", "F1_curve.png", "confusion_matrix.csv"))))

  if (has_csv || has_weights || has_plots) {
    return(dir)
  }

  # 2. Search recursively for results.csv
  c_csv <- list.files(dir, pattern = "^results\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(c_csv) > 0L) {
    return(normalizePath(dirname(c_csv[1]), winslash = "/", mustWork = FALSE))
  }

  # 3. Search recursively for weights/best.pt or weights/best.onnx
  c_weights <- list.files(dir, pattern = "^(best|last)\\.(pt|onnx)$", recursive = TRUE, full.names = TRUE)
  if (length(c_weights) > 0L) {
    w_dir <- dirname(c_weights[1])
    parent_dir <- dirname(w_dir)
    return(normalizePath(parent_dir, winslash = "/", mustWork = FALSE))
  }

  # 4. Search recursively for confusion_matrix or any training plot
  c_plots <- list.files(dir, pattern = "^(confusion_matrix|results|F1_curve|PR_curve)\\.(png|jpg|csv)$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(c_plots) > 0L) {
    return(normalizePath(dirname(c_plots[1]), winslash = "/", mustWork = FALSE))
  }

  dir
}

.resolve_yolo_results_dir <- function(model = NULL, results_dir = NULL, output_dir = NULL) {
  # 1. Direct path given via results_dir
  if (!is.null(results_dir) && nzchar(as.character(results_dir)[1])) {
    rd <- as.character(results_dir)[1]
    if (file.exists(rd) && !dir.exists(rd) && grepl("\\.zip$", rd, ignore.case = TRUE)) {
      ex_dir <- sub("\\.zip$", "", rd, ignore.case = TRUE)
      if (!dir.exists(ex_dir)) {
        dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
        utils::unzip(rd, exdir = ex_dir)
      }
      return(.find_yolo_subrun_dir(ex_dir))
    }
    if (dir.exists(rd)) {
      return(.find_yolo_subrun_dir(rd))
    }
    rd_cwd <- file.path(getwd(), rd)
    if (dir.exists(rd_cwd)) {
      return(.find_yolo_subrun_dir(rd_cwd))
    }
  }

  target <- if (!is.null(model)) as.character(model)[1] else NULL
  if (is.null(target) || !nzchar(target)) return(NULL)

  # 2. Check if target is directly an existing zip archive
  if (file.exists(target) && !dir.exists(target) && grepl("\\.zip$", target, ignore.case = TRUE)) {
    ex_dir <- sub("\\.zip$", "", target, ignore.case = TRUE)
    if (!dir.exists(ex_dir)) {
      dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
      utils::unzip(target, exdir = ex_dir)
    }
    return(.find_yolo_subrun_dir(ex_dir))
  }

  # 3. Check if target is directly an existing directory
  if (dir.exists(target)) {
    return(.find_yolo_subrun_dir(target))
  }

  # 3b. Check if target is a file in a results folder or weights subfolder
  if (file.exists(target) && !dir.exists(target)) {
    p_dir <- dirname(target)
    if (basename(p_dir) == "weights") {
      gp_dir <- dirname(p_dir)
      if (file.exists(file.path(gp_dir, "args.yaml")) || file.exists(file.path(gp_dir, "results.csv")) || dir.exists(file.path(gp_dir, "weights"))) {
        return(.find_yolo_subrun_dir(gp_dir))
      }
    } else if (file.exists(file.path(p_dir, "args.yaml")) || file.exists(file.path(p_dir, "results.csv")) || dir.exists(file.path(p_dir, "weights"))) {
      return(.find_yolo_subrun_dir(p_dir))
    }
  }

  # 4. Check in current working directory (getwd())
  cwd_path <- file.path(getwd(), target)
  if (dir.exists(cwd_path)) {
    return(.find_yolo_subrun_dir(cwd_path))
  }
  if (file.exists(paste0(cwd_path, ".zip"))) {
    return(.resolve_yolo_results_dir(model = paste0(cwd_path, ".zip")))
  }

  base_stem <- sub("\\.(onnx|pt|zip)$", "", basename(target), ignore.case = TRUE)
  base_stem <- sub("_results$", "", base_stem, ignore.case = TRUE)
  stem_clean <- sub("^yolo[-_]", "", base_stem, ignore.case = TRUE)

  # 5. Build candidate directories across working directory, output_dir, and pliman_model_dir()
  search_roots <- unique(c(
    getwd(),
    if (!is.null(output_dir)) output_dir,
    pliman_model_dir(),
    dirname(target)
  ))
  search_roots <- search_roots[dir.exists(search_roots)]

  cand_dirs <- character(0)
  for (root in search_roots) {
    cand_dirs <- c(cand_dirs,
      file.path(root, target),
      file.path(root, paste0(target, "_results")),
      file.path(root, paste0(base_stem, "_results")),
      file.path(root, paste0(stem_clean, "_results")),
      file.path(root, base_stem),
      file.path(root, stem_clean),
      file.path(root, "runs_yolo", base_stem),
      file.path(root, "runs_yolo", paste0(base_stem, "temp")),
      file.path(root, "runs", "detect", base_stem),
      file.path(root, "runs", "segment", base_stem),
      file.path(root, "runs", "detect", "train"),
      file.path(root, "runs", "detect", "val"),
      file.path(root, "runs", "segment", "train"),
      file.path(root, "runs", "segment", "val")
    )

    # Check zip files in root
    cand_zips <- c(
      file.path(root, paste0(target, ".zip")),
      file.path(root, paste0(base_stem, "_results.zip")),
      file.path(root, paste0(base_stem, ".zip"))
    )
    for (cz in cand_zips) {
      if (file.exists(cz) && !dir.exists(cz)) {
        res_z <- .resolve_yolo_results_dir(model = cz)
        if (!is.null(res_z)) return(res_z)
      }
    }

    # Pattern matching in root
    pat <- paste0("(?i)(^|_)", stem_clean, ".*(_results)?$")
    matched_dirs <- list.files(root, pattern = pat, full.names = TRUE)
    cand_dirs <- c(cand_dirs, matched_dirs[dir.exists(matched_dirs)])
  }

  cand_dirs <- unique(cand_dirs[dir.exists(cand_dirs)])
  if (length(cand_dirs) > 0L) {
    for (cd in cand_dirs) {
      subrun <- .find_yolo_subrun_dir(cd)
      if (!is.null(subrun)) return(subrun)
    }
  }

  # 6. Fuzzy search in getwd(): any folder containing stem_clean
  if (nzchar(stem_clean)) {
    fuzzy_dirs <- list.files(getwd(), pattern = paste0("(?i)", stem_clean), full.names = TRUE)
    fuzzy_dirs <- fuzzy_dirs[dir.exists(fuzzy_dirs)]
    for (fd in fuzzy_dirs) {
      subrun <- .find_yolo_subrun_dir(fd)
      if (!is.null(subrun)) return(subrun)
    }
  }

  NULL
}

#' @title Inspect YOLO Training Results and Diagnostics
#' @name yolo_results
#' @description
#' `yolo_results()` loads, inspects, and plots diagnostics, training metrics, loss curves,
#' confusion matrices, and precision-recall graphs for a trained YOLO model.
#'
#' It supports models trained locally with [yolo_train()] or [yolo_train_workflow()], as well
#' as training runs downloaded from Google Colab, cloud VMs, or external directories.
#'
#' You can pass:
#' * A model name (e.g. `"capsulas_det"`, `"model_capsulas"`).
#' * A direct path to an external results folder (e.g. `"runs/detect/train"`, `"capsulas_small_results"`, or downloaded from Colab).
#' * A direct path to a `.zip` archive containing the downloaded training results (automatically extracted).
#'
#' @param model A model name, a direct path to a results directory, or path to a results `.zip` archive.
#'   Can be omitted if `results_dir` is explicitly supplied.
#' @param results_dir Optional direct path to the training results folder. If `model` is already
#'   a directory or `.zip` path, this is inferred automatically.
#' @param output_dir Base directory to search for trained models. Defaults to [pliman_model_dir()].
#'
#' @return An object of class `c("yolo_train", "list")`.
#' @export
#' @seealso [plot.yolo_train()], [yolo_plot()], [yolo_train()]
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Inspect a model by name
#' res <- yolo_results("capsulas_det")
#'
#' # 2. Inspect a folder downloaded from Google Colab
#' res_colab <- yolo_results("capsulas_small_results")
#'
#' # 3. Inspect a .zip archive directly
#' res_zip <- yolo_results("capsulas_small_results.zip")
#'
#' # Plot diagnostic dashboard
#' plot(res)
#' plot(res, which = "confusion_matrix")
#' plot(res, which = "pr_curve")
yolo_results <- function(model = NULL,
                         results_dir = NULL,
                         output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- pliman_model_dir()
  } else {
    output_dir <- normalizePath(output_dir[1], winslash = "/", mustWork = FALSE)
  }

  if (is.list(model)) {
    if (!is.null(model$results_dir) && dir.exists(model$results_dir)) {
      results_dir <- results_dir %||% model$results_dir
    }
    model <- model$model_file %||% model$results_dir %||% model$best_onnx %||% model$best_pt %||% NULL
  }

  if (is.null(model) && is.null(results_dir)) {
    cli::cli_abort("Please provide either {.arg model} or {.arg results_dir} (e.g. a model name, results folder path, or results .zip).")
  }

  # Resolve the results folder
  found_results_dir <- .resolve_yolo_results_dir(model = model, results_dir = results_dir, output_dir = output_dir)

  if (is.null(found_results_dir) || !dir.exists(found_results_dir)) {
    target_name <- if (!is.null(model)) as.character(model)[1] else as.character(results_dir)[1]
    avail_pliman <- list.files(output_dir, pattern = "_results$")
    avail_pliman_names <- sub("_results$", "", avail_pliman)
    avail_cwd <- list.files(getwd(), pattern = "_results$")
    avail_cwd_names <- sub("_results$", "", avail_cwd)

    cli::cli_abort(c(
      "Could not find training results directory for {.val {target_name}}.",
      "i" = "Searched in: {.path {getwd()}}, {.path {output_dir}}, and {.path {pliman_model_dir()}}.",
      if (length(avail_pliman_names) > 0L) c("*" = "Available models in pliman: {.val {avail_pliman_names}}") else NULL,
      if (length(avail_cwd_names) > 0L) c("*" = "Available results folders in working directory: {.val {avail_cwd_names}}") else NULL,
      "v" = "Tip: You can pass the exact path to the downloaded Colab results folder or .zip file directly:\n  {.code yolo_results(\"path/to/results_folder\")} or {.code yolo_results(\"results.zip\")}"
    ))
  }

  results_dir <- found_results_dir
  base_stem <- sub("_results$", "", basename(results_dir), ignore.case = TRUE)
  if (base_stem %in% c("val", "train") && basename(dirname(results_dir)) %in% c("detect", "segment", "classify")) {
    base_stem <- basename(dirname(dirname(results_dir)))
  }

  # Find onnx and pt models (check results_dir first, then parent, then output_dir)
  dest_onnx <- file.path(results_dir, "weights", "best.onnx")
  if (!file.exists(dest_onnx)) {
    cand_onnx <- c(
      file.path(results_dir, "best.onnx"),
      file.path(dirname(results_dir), paste0(base_stem, ".onnx")),
      file.path(output_dir, paste0(base_stem, ".onnx")),
      list.files(results_dir, pattern = "\\.onnx$", recursive = TRUE, full.names = TRUE),
      list.files(dirname(results_dir), pattern = "\\.onnx$", full.names = TRUE)
    )
    cand_onnx <- cand_onnx[file.exists(cand_onnx)]
    dest_onnx <- if (length(cand_onnx) > 0L) cand_onnx[1] else NULL
  }

  dest_pt <- file.path(results_dir, "weights", "best.pt")
  if (!file.exists(dest_pt)) {
    cand_pt <- c(
      file.path(results_dir, "best.pt"),
      file.path(dirname(results_dir), paste0(base_stem, ".pt")),
      file.path(output_dir, paste0(base_stem, ".pt")),
      list.files(results_dir, pattern = "(best|last)\\.pt$", recursive = TRUE, full.names = TRUE),
      list.files(dirname(results_dir), pattern = "\\.pt$", full.names = TRUE)
    )
    cand_pt <- cand_pt[file.exists(cand_pt)]
    dest_pt <- if (length(cand_pt) > 0L) cand_pt[1] else NULL
  }

  find_plot <- function(patterns) {
    for (pat in patterns) {
      matches <- list.files(results_dir, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
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
      return(tryCatch(utils::read.csv(p, check.names = FALSE, strip.white = TRUE), error = function(e) NULL))
    }
    matches <- list.files(results_dir, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(matches) > 0L) {
      return(tryCatch(utils::read.csv(matches[1], check.names = FALSE, strip.white = TRUE), error = function(e) NULL))
    }
    NULL
  }

  cm_df <- read_csv_safe("confusion_matrix.csv")
  cm_norm_df <- read_csv_safe("confusion_matrix_normalized.csv")
  aligned_cm <- .align_yolo_confusion_matrix(cm_counts = cm_df, cm_norm = cm_norm_df)
  cm_df <- aligned_cm$cm_counts
  cm_norm_df <- aligned_cm$cm_norm
  pr_df <- read_csv_safe("pr_curve.csv")
  f1_df <- read_csv_safe("f1_curve.csv")

  results_csv <- file.path(results_dir, "results.csv")
  if (!file.exists(results_csv)) {
    c_csv <- list.files(results_dir, pattern = "^results\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(c_csv) > 0L) results_csv <- c_csv[1]
  }

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
      map50_col <- map50_col[!grepl("95", map50_col)]
      if (length(map50_col) > 0L) {
        best_idx <- which.max(metrics_df[[map50_col[1]]])
        if (length(best_idx) > 0L) {
          best_metrics <- metrics_df[best_idx, , drop = FALSE]
        }
      }
    }, error = function(e) NULL)
  }

  # Detect task
  task_detected <- "detect"
  args_yaml <- file.path(results_dir, "args.yaml")
  if (!file.exists(args_yaml)) {
    c_yaml <- list.files(results_dir, pattern = "^args\\.yaml$", recursive = TRUE, full.names = TRUE)
    if (length(c_yaml) > 0L) args_yaml <- c_yaml[1]
  }
  if (file.exists(args_yaml)) {
    tryCatch({
      ylines <- readLines(args_yaml, warn = FALSE)
      t_line <- grep("^task:", ylines, value = TRUE)
      if (length(t_line) > 0L) {
        val <- trimws(sub("^task:\\s*", "", t_line[1]))
        if (val %in% c("detect", "segment", "classify")) task_detected <- val
      }
    }, error = function(e) NULL)
  } else if (!is.null(metrics_df)) {
    if (any(grepl("seg", names(metrics_df), ignore.case = TRUE)) || any(grepl("mask", names(metrics_df), ignore.case = TRUE))) {
      task_detected <- "segment"
    } else if (any(grepl("acc", names(metrics_df), ignore.case = TRUE)) && !any(grepl("box", names(metrics_df), ignore.case = TRUE))) {
      task_detected <- "classify"
    }
  } else if (grepl("seg", base_stem, ignore.case = TRUE)) {
    task_detected <- "segment"
  } else if (grepl("cls|class", base_stem, ignore.case = TRUE)) {
    task_detected <- "classify"
  }

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
      task = task_detected,
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
.plot_yolo_metrics_native <- function(df, task = "detect", bg = "white", ...) {
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
      graphics::lines(ep, vl_box, col = "#DC2626", lwd = 2.2, lty = 1)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
      graphics::lines(ep, vl_seg, col = "#DC2626", lwd = 2.2, lty = 1)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
      graphics::lines(ep, vl_cls, col = "#DC2626", lwd = 2.2, lty = 1)
      graphics::legend("topright", legend = c("Train", "Val"), col = c("#2563EB", "#DC2626"),
                       lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
    graphics::lines(ep, r_val, col = "#D97706", lwd = 2.2, lty = 1)
    graphics::legend("bottomright", legend = c("Precision", "Recall"), col = c("#059669", "#D97706"),
                     lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
    if (!is.null(m95)) graphics::lines(ep, m95, col = "#EA580C", lwd = 2.5, lty = 1)
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
                     lty = if (!is.null(m95)) c(1, 1) else 1, lwd = 2.5, bty = "n", cex = 0.85)
  }
  # 5. F1 Score or Box mAP
  if (has_seg && !is.null(c_map50_b)) {
    m50_b <- df[[c_map50_b]]
    m95_b <- if (!is.null(c_map95_b)) df[[c_map95_b]] else NULL
    graphics::plot(ep, m50_b, type = "n", xlab = "Epoch", ylab = "mAP",
                   main = "Box Validation mAP", ylim = c(0, 1.05), cex.main = 1.1, font.main = 2)
    add_grid()
    graphics::lines(ep, m50_b, col = "#0284C7", lwd = 2.5)
    if (!is.null(m95_b)) graphics::lines(ep, m95_b, col = "#D97706", lwd = 2.5, lty = 1)
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
                     lty = if (!is.null(m95_b)) c(1, 1) else 1, lwd = 2.5, bty = "n", cex = 0.85)
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
# Helper: Align and normalize confusion matrix
# Ensures rows = Ground Truth and columns = Predicted Class
.align_yolo_confusion_matrix <- function(cm_counts = NULL, cm_norm = NULL) {
  if (is.null(cm_counts) && is.null(cm_norm)) return(list(cm_counts = NULL, cm_norm = NULL))

  clean_cm_df <- function(df) {
    if (is.null(df) || !is.data.frame(df) || ncol(df) <= 2L) return(df)
    na_cols <- vapply(df[, -1, drop = FALSE], function(col) all(is.na(col)), logical(1))
    if (any(na_cols)) {
      df <- df[, c(TRUE, !na_cols), drop = FALSE]
    }
    df
  }
  if (!is.null(cm_counts)) cm_counts <- clean_cm_df(cm_counts)
  if (!is.null(cm_norm)) cm_norm <- clean_cm_df(cm_norm)

  base_df <- if (!is.null(cm_counts)) cm_counts else cm_norm
  if (!is.data.frame(base_df) || ncol(base_df) <= 1L) {
    return(list(cm_counts = cm_counts, cm_norm = cm_norm))
  }

  first_col <- names(base_df)[1]
  # Older exports used 'true_class' but contained un-transposed Ultralytics matrix [pred, true]
  needs_transpose <- identical(first_col, "true_class")

  if (needs_transpose) {
    row_labels <- as.character(base_df[[1]])
    col_labels <- names(base_df)[-1]
    if (length(row_labels) == length(col_labels) &&
        grepl("^class_\\d+$", row_labels[length(row_labels)]) &&
        any(grepl("^background$", col_labels, ignore.case = TRUE))) {
      row_labels[length(row_labels)] <- "background"
    }

    if (!is.null(cm_counts)) {
      mat_c <- as.matrix(cm_counts[, -1, drop = FALSE])
      if (nrow(mat_c) == ncol(mat_c)) {
        mat_c_t <- t(mat_c)
        cm_counts <- stats::setNames(
          data.frame(ground_truth = col_labels, as.data.frame(mat_c_t), check.names = FALSE, stringsAsFactors = FALSE),
          c("ground_truth", row_labels)
        )
      }
    }
    if (!is.null(cm_norm) && is.null(cm_counts)) {
      mat_n <- as.matrix(cm_norm[, -1, drop = FALSE])
      if (nrow(mat_n) == ncol(mat_n)) {
        mat_n_t <- t(mat_n)
        cm_norm <- stats::setNames(
          data.frame(ground_truth = col_labels, as.data.frame(mat_n_t), check.names = FALSE, stringsAsFactors = FALSE),
          c("ground_truth", row_labels)
        )
      }
    }
  }

  if (!is.null(cm_counts)) names(cm_counts)[1] <- "ground_truth"
  if (!is.null(cm_norm)) names(cm_norm)[1] <- "ground_truth"

  # Recalculate normalized matrix from counts to guarantee correct row/column proportions
  if (!is.null(cm_counts)) {
    row_labels <- as.character(cm_counts[[1]])
    col_labels <- names(cm_counts)[-1]
    mat_c <- as.matrix(cm_counts[, -1, drop = FALSE])
    mat_c[is.na(mat_c)] <- 0
    nr <- nrow(mat_c)
    nc <- ncol(mat_c)
    mat_n <- matrix(0, nrow = nr, ncol = nc)
    for (r in seq_len(nr)) {
      r_name <- tolower(row_labels[r])
      if (r_name == "background" || grepl("^class_\\d+$", r_name)) {
        # False Positives on background: normalize by column sum (detections for that class)
        c_sums <- colSums(mat_c)
        for (c in seq_len(nc)) {
          mat_n[r, c] <- if (c_sums[c] > 0 && c <= length(c_sums) && tolower(col_labels[c]) != "background") {
            round(mat_c[r, c] / c_sums[c], 4)
          } else {
            0
          }
        }
      } else {
        # True class instances: normalize by row sum (Total Ground Truth Instances -> Sensitivity / Recall)
        r_sum <- sum(mat_c[r, ])
        mat_n[r, ] <- if (r_sum > 0) round(mat_c[r, ] / r_sum, 4) else 0
      }
    }
    cm_norm <- stats::setNames(
      data.frame(ground_truth = row_labels, as.data.frame(mat_n), check.names = FALSE, stringsAsFactors = FALSE),
      c("ground_truth", col_labels)
    )
  }

  list(cm_counts = cm_counts, cm_norm = cm_norm)
}

# Helper: Native R vector Confusion Matrix heatmap
.plot_yolo_confusion_matrix_native <- function(cm_norm = NULL, cm_counts = NULL, normalized = TRUE, bg = "white", ...) {
  if (is.null(cm_norm) && is.null(cm_counts)) {
    cli::cli_alert_warning("No confusion matrix data available to plot.")
    return(invisible(NULL))
  }
  aligned <- .align_yolo_confusion_matrix(cm_counts = cm_counts, cm_norm = cm_norm)
  cm_counts <- aligned$cm_counts
  cm_norm <- aligned$cm_norm

  if (is.null(cm_counts) && !is.null(cm_norm)) {
    cm_counts <- cm_norm
  }
  if (is.null(cm_norm) && !is.null(cm_counts)) {
    cm_norm <- cm_counts
  }

  row_labels <- as.character(cm_norm[[1]])
  col_labels <- names(cm_norm)[-1]
  if (length(row_labels) == length(col_labels) && grepl("^class_\\d+$", row_labels[length(row_labels)]) && any(grepl("^background$", col_labels, ignore.case = TRUE))) {
    row_labels[length(row_labels)] <- "background"
  }
  mat_norm <- as.matrix(cm_norm[, -1, drop = FALSE])
  mat_counts <- as.matrix(cm_counts[, -1, drop = FALSE])
  mat_norm[is.na(mat_norm)] <- 0
  mat_counts[is.na(mat_counts)] <- 0
  nr <- nrow(mat_norm)
  nc <- ncol(mat_norm)
  is_binary_2x2 <- (nr == 2L && nc == 2L)

  op <- graphics::par(no.readonly = TRUE)
  if (identical(graphics::par("mfrow"), c(1L, 1L))) {
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = c(5, 6.5, 3.5, 2), bg = bg)
  }
  main_title <- if (isTRUE(normalized)) "Confusion Matrix (Normalized)" else "Confusion Matrix (Counts)"
  graphics::plot(1, type = "n", xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
                 axes = FALSE, xlab = "Predicted value", ylab = "",
                 main = main_title, font.main = 2, cex.main = 1.15)
  graphics::title(ylab = "Actual value", line = 4.2, font.lab = 2)

  # Modern harmonious palettes:
  # Diagonal (TRUEs: TP, TN) -> Green tones
  # Off-diagonal (FALSEs: FP, FN, misclassifications) -> Soft Rose/Red tones
  greens <- grDevices::colorRampPalette(c("#f0fdf4", "#86efac", "#22c55e", "#15803d"))(100)
  reds   <- grDevices::colorRampPalette(c("#fff1f2", "#fecdd3", "#f43f5e", "#be123c"))(100)

  tags_2x2 <- matrix(c("TP", "FP", "FN", "TN"), nrow = 2, ncol = 2)

  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      val <- mat_norm[r, c]
      cnt <- mat_counts[r, c]
      y_pos <- nr - r + 1
      x_pos <- c
      is_diag <- (r == c)

      if (is_diag) {
        # Diagonal = TRUEs (TP / TN) -> Green
        if (is_binary_2x2 && r == 2L) {
          # TN in object detection (background correctly unpredicted)
          intensity <- 0.45
          bg_col <- greens[45]
          txt_col <- "#14532d"
        } else {
          intensity <- if (is.na(val) || val <= 0) 0 else max(0.15, min(1, val))
          col_idx <- max(1, min(100, round(intensity * 99) + 1))
          bg_col <- greens[col_idx]
          txt_col <- if (intensity > 0.6) "#ffffff" else "#14532d"
        }
      } else {
        # Off-diagonal = FALSEs (FP / FN / Misclassifications) -> Red/Pink
        if (cnt == 0) {
          bg_col <- "#ffffff"
          txt_col <- "#64748b"
        } else {
          intensity <- if (is.na(val) || val <= 0) 0.25 else max(0.2, min(1, val))
          col_idx <- max(1, min(100, round(intensity * 99) + 1))
          bg_col <- reds[col_idx]
          txt_col <- if (intensity > 0.6) "#ffffff" else "#881337"
        }
      }

      graphics::rect(x_pos - 0.48, y_pos - 0.48, x_pos + 0.48, y_pos + 0.48,
                     col = bg_col, border = "#cbd5e1", lwd = 1.3)

      if (is_binary_2x2) {
        tag <- tags_2x2[r, c]
        graphics::text(x_pos, y_pos + 0.22, labels = tag, col = txt_col, font = 2, cex = 1.1)

        cnt_lbl <- if (is_diag && r == 2L && cnt == 0) "0" else sprintf("%d", cnt)
        graphics::text(x_pos, y_pos - 0.02, labels = cnt_lbl, col = txt_col, font = 2, cex = 1.25)

        rate_lbl <- if (r == 1L) {
          sprintf("(%.1f%%)", val * 100)
        } else if (r == 2L && c == 1L) {
          sprintf("(FP: %d)", cnt)
        } else {
          "(TN)"
        }
        graphics::text(x_pos, y_pos - 0.24, labels = rate_lbl, col = txt_col, cex = 0.85)
      } else {
        if (isTRUE(normalized)) {
          graphics::text(x_pos, y_pos + 0.08, labels = sprintf("%.2f", val), col = txt_col, font = 2, cex = 1.05)
          graphics::text(x_pos, y_pos - 0.18, labels = sprintf("(%d)", cnt), col = txt_col, cex = 0.82)
        } else {
          graphics::text(x_pos, y_pos + 0.08, labels = sprintf("%d", cnt), col = txt_col, font = 2, cex = 1.05)
          graphics::text(x_pos, y_pos - 0.18, labels = sprintf("(%.1f%%)", val * 100), col = txt_col, cex = 0.82)
        }
      }
    }
  }
  cex_axis <- if (nc > 6) 0.75 else 0.95
  graphics::axis(1, at = seq_len(nc), labels = col_labels, font = 2, tick = FALSE, line = -0.5, cex.axis = cex_axis)
  graphics::axis(2, at = rev(seq_len(nr)), labels = row_labels, font = 2, tick = FALSE, line = -0.5, las = 2, cex.axis = cex_axis)
  invisible(cm_norm)
}
# Helper: Native R vector Precision-Recall Curve plot
.plot_yolo_pr_curve_native <- function(pr_df, mAP50 = NULL, bg = "white", ...) {
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
.plot_yolo_f1_curve_native <- function(f1_df, bg = "white", ...) {
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
# Helper: Native R comprehensive diagnostic dashboard
.plot_yolo_all_dashboard_native <- function(x, bg = "white", ncol = NULL, nrow = NULL, digits = 3, ...) {
  df <- x$metrics
  has_metrics <- !is.null(df) && is.data.frame(df) && nrow(df) > 0L
  has_cm <- !is.null(x$confusion_matrix_norm) || !is.null(x$confusion_matrix)
  has_pr <- !is.null(x$pr_curve)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  num_col <- if (!is.null(ncol)) as.integer(ncol) else if (!is.null(nrow)) ceiling(4L / as.integer(nrow)) else 2L
  num_row <- if (!is.null(nrow)) as.integer(nrow) else ceiling(4L / num_col)
  graphics::par(mfrow = c(num_row, num_col), mar = c(4.2, 5, 3, 1.5), mgp = c(2.2, 0.7, 0), bg = bg)
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
      if (!is.null(vl_loss)) graphics::lines(ep, vl_loss, col = "#dc2626", lwd = 2.2, lty = 1)
      graphics::legend("topright", legend = c("Train Loss", "Val Loss"),
                       col = c("#2563eb", "#dc2626"), lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
      if (!is.null(vl_box)) graphics::lines(ep, vl_box, col = "#dc2626", lwd = 2.2, lty = 1)
      if (!is.null(tr_cls)) graphics::lines(ep, tr_cls, col = "#059669", lwd = 2)
      if (!is.null(vl_cls)) graphics::lines(ep, vl_cls, col = "#d97706", lwd = 2, lty = 1)
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
        if (!is.null(vl_box)) 1,
        if (!is.null(tr_cls)) 1,
        if (!is.null(vl_cls)) 1
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
        if (!is.null(top5_v)) graphics::lines(ep, top5_v, col = "#ea580c", lwd = 2.5, lty = 1)
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
                         lty = if (!is.null(top5_v)) c(1, 1) else 1, lwd = 2.5, bty = "n", cex = 0.85)
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
        if (!is.null(m95)) graphics::lines(ep, m95, col = "#ea580c", lwd = 2.5, lty = 1)
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
                         lty = if (!is.null(m95)) c(1, 1) else 1, lwd = 2.5, bty = "n", cex = 0.85)
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
      graphics::lines(ep, df[[c_rec]], col = "#d97706", lwd = 2.2, lty = 1)
      graphics::legend("bottomright", legend = c("Precision", "Recall"), col = c("#059669", "#d97706"),
                       lty = c(1, 1), lwd = 2, bty = "n", cex = 0.85)
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
#' @param bg Background color for native plots. Defaults to `"white"`.
#' @param mar Plot margins passed to base graphics. Defaults to `0.5`.
#' @param ncol Optional integer. Number of columns in multi-panel layout.
#' @param nrow Optional integer. Number of rows in multi-panel layout.
#' @param digits Integer. Number of decimal digits for statistical displays (default `3`).
#' @param save Logical. If `TRUE`, saves the generated plot to a file. Defaults to `FALSE`.
#' @param file_name Optional character string. Custom file name or path to save the plot (e.g. `"metrics.png"`, `"dashboard.pdf"`).
#'   If not specified and `save = TRUE`, a default filename based on the model name and plotted metric is automatically generated.
#' @param width Numeric. Width of the output plot. If `<= 50`, interpreted as inches; if `> 50`, interpreted as pixels.
#'   Defaults to `10` inches for the multi-panel dashboard (`which = "all"`), `9` for metrics curves, and `7` for single curves.
#' @param height Numeric. Height of the output plot. If `<= 50`, interpreted as inches; if `> 50`, interpreted as pixels.
#'   Defaults to `8` inches for the multi-panel dashboard (`which = "all"`), and `6` for single curves.
#' @param dpi Numeric. Resolution in dots per inch (DPI) for raster formats (PNG, JPEG, TIFF). Defaults to `300`.
#' @param ... Additional arguments passed to base graphics.
#' @return Invisibly returns `x` or plotted data frame.
#' @export
#' @seealso [yolo_plot()], [yolo_results()], [yolo_train()]
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Load diagnostics from a previously trained model
#' res <- yolo_results("capsulas_det")
#'
#' # Inspect training summary and plot 2x2 diagnostic dashboard
#' print(res)
#' plot(res)
#'
#' # Save high-resolution 300 DPI publication plots:
#' plot(res, save = TRUE)
#' plot(
#'   res,
#'   which = "confusion_matrix_norm",
#'   save = TRUE,
#'   file_name = "fig_cm.png",
#'   width = 8,
#'   height = 7,
#'   dpi = 300
#' )
#' plot(res, which = "metrics", save = TRUE, file_name = "fig_metrics.pdf", width = 9, height = 6)
#'
#' # Plot individual evaluation curves in high-resolution vector format:
#' plot(res, which = "confusion_matrix")
#' plot(res, which = "confusion_matrix_norm")
#' plot(res, which = "pr_curve")
#' plot(res, which = "f1_curve")
#' plot(res, which = "metrics")
#' @export
plot.yolo_train <- function(x,
                            which = c("all", "metrics", "confusion_matrix_norm",
                                      "confusion_matrix", "pr_curve", "f1_curve",
                                      "results", "labels"),
                            native = TRUE,
                            bg = "white",
                            mar = 0.5,
                            ncol = NULL,
                            nrow = NULL,
                            digits = 3,
                            save = FALSE,
                            file_name = NULL,
                            width = NULL,
                            height = NULL,
                            dpi = 300,
                            ...) {
  if (is.character(x)) {
    x <- yolo_results(x)
  }
  if (!inherits(x, "yolo_train")) {
    cli::cli_abort("Object must be of class {.code yolo_train} or a model name/path.")
  }
  which <- match.arg(which)

  # Setup saving device if requested
  model_tag <- if (!is.null(x$model_name) && nzchar(x$model_name)) {
    x$model_name
  } else if (!is.null(x$results_dir) && nzchar(x$results_dir)) {
    basename(x$results_dir)
  } else {
    "yolo_train"
  }
  tag_which <- if (which == "all") "dashboard" else which
  default_file <- paste0(model_tag, "_", tag_which, ".png")

  def_w <- if (which == "all") 10 else if (which %in% c("metrics", "results")) 9 else 7
  def_h <- if (which == "all") 8 else 6

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
  if (!is.null(x$confusion_matrix) || !is.null(x$confusion_matrix_norm)) {
    aligned_cm <- .align_yolo_confusion_matrix(cm_counts = x$confusion_matrix, cm_norm = x$confusion_matrix_norm)
    x$confusion_matrix <- aligned_cm$cm_counts
    x$confusion_matrix_norm <- aligned_cm$cm_norm
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
      return(invisible(.plot_yolo_all_dashboard_native(x, bg = bg, ncol = ncol, nrow = nrow, digits = digits, ...)))
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
#' @examplesIf interactive()
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
yolo_plot <- function(x, ...) {
  if (is.character(x)) {
    x <- yolo_results(x)
  }
  plot(x, ...)
}
# ==============================================================================
# YOLO ACTIVE LEARNING & BOOTSTRAP TRAINING WORKFLOW
# ==============================================================================
#' @title End-to-End YOLO Active Learning & Bootstrap Training Workflow
#' @name yolo_train_workflow
#' @aliases yolo_workflow
#' @description
#' `yolo_train_workflow()` implements a complete active-learning bootstrap pipeline:
#' 1. **Few-Shot Sampling**: Selects a small subset of representative sample images
#'    (e.g., 3–4 images) from an unlabelled dataset directory.
#' 2. **Auto-Annotation via Few-Shot / Exemplar Matching**: Uses [yolo_dataset_export()] with Meta PerSAM (or Grounded-SAM)
#'    to accurately annotate the sampled images.
#' 3. **Offline Augmentation**: Optionally expands those few images with [yolo_dataset_augment()] into 20–50+ diverse images
#'    (3D perspective warping, contrast CLAHE, random cutout, noise, blur, rotations).
#' 4. **High-Resolution Tiling**: Optionally slices images into uniform tiles with [yolo_dataset_tiling()] (e.g. 640x640),
#'    optimizing small object resolution and dataset volume.
#' 5. **Automated Training**: Trains a high-performance YOLO model via [yolo_train()].
#' 6. **Model Export**: Produces ready-to-deploy `.onnx` and `.pt` model files that can be directly passed
#'    into [yolo_dataset_autolabel()] to annotate thousands of remaining images in a few seconds!
#'
#' @param data Directory path containing unlabelled images (e.g. `"dataset_complex"`), or a vector of image paths.
#' @param sample_n Number of images to sample from `data` if `sample_images` is not provided (default `4L`).
#' @param sample_images Optional vector of specific image paths or filenames to label. If `NULL` (default),
#'   samples `sample_n` images from `data`.
#' @param model_fewshot Few-shot / exemplar model to label the sample images. Defaults to `"persam"`.
#'   Can also be `"grounded-sam"`, `"ben2"`, or an existing model.
#' @param class_names Names of the classes to detect. Defaults to `"object"`.
#' @param task Target task: `"detect"` (bounding boxes) or `"segment"` (polygons). Defaults to `"detect"`.
#' @param prompt Visual exemplar prompt points or bounding box for PerSAM, or text prompt for Grounded-SAM.
#'   If `NULL` (default) and `interactive = TRUE`, opens an interactive locator window on the first sampled image.
#' @param prompt_mode Prompt reuse mode: `"first"` (default, prompt once on image 1 and propagate to all sample images)
#'   or `"each"` (prompt on each sample image).
#' @param sim_threshold Similarity threshold for PerSAM (default `0.70`).
#' @param crop_to_objects Logical. Whether to crop sampled images to object bounding boxes with margin. Defaults to `TRUE`.
#' @param crop_margin Margin proportion around objects when cropping (default `0.10`).
#' @param augment Logical. Whether to apply offline data augmentation (default `TRUE`).
#' @param augment_total Target total number of augmented training images (default `30L`).
#' @param augment_args List of custom arguments passed to [yolo_dataset_augment()].
#' @param tiling Logical. Whether to apply dataset tiling (default `TRUE`).
#' @param tile_size Tile size in pixels (default `640`).
#' @param tiling_args List of custom arguments passed to [yolo_dataset_tiling()].
#' @param train_model Base pre-trained YOLO architecture to train (default `"yolo26n"`, or `"yolo11n"`, `"yolo26s"`).
#' @param epochs Number of training epochs (default `50`).
#' @param batch Batch size (default `8L`, or `4L` on low-VRAM GPUs).
#' @param device Device to train on (`"auto"`, `0`, or `"cpu"`). Defaults to `"auto"`.
#' @param close_mosaic Number of final epochs to disable mosaic augmentation (default `15`).
#' @param output_dir Target directory to save the final exported `.onnx` and `.pt` model files.
#'   Defaults to [pliman_model_dir()].
#' @param output_name Filename for the exported model (e.g. `"model_capsulas.onnx"`). If `NULL`,
#'   an informative name based on `data` is automatically generated.
#' @param temp_dir Directory used to construct the bootstrap training dataset. If `NULL` (default),
#'   creates a temporary folder inside [tempdir()].
#' @param dataset_dir Optional character directory path to save or reuse the constructed YOLO training dataset.
#'   If `NULL` (default), uses `temp_dir` or creates a persistent directory when `keep_dataset = TRUE`.
#' @param keep_dataset Logical. If `FALSE` (default), removes the temporary bootstrap dataset after training.
#' @param as_job Logical. If `TRUE` and running inside RStudio, launches the entire
#'   workflow (bootstrapping, augmentation, tiling, and YOLO training) as an asynchronous
#'   RStudio Background Job, immediately freeing up the R console. Default is `FALSE`.
#'   If `TRUE`, keeps the folder intact.
#' @param review Logical. Whether to generate review preview images during few-shot labeling (default `FALSE`).
#' @param show_class Logical. Whether to display class names on labels/badges during review (default `FALSE`).
#' @param show_conf Logical. Whether to display confidence values on labels/badges during review (default `FALSE`).
#' @param badge Logical. Whether to draw a filled badge behind label text (default `FALSE`).
#' @param interactive Logical. Whether to allow interactive clicking/refinement (default `interactive()`).
#' @param verbose Logical. Whether to display informative progress messages (default `TRUE`).
#' @param ... Additional arguments passed to [yolo_train()] (e.g., `scale = 0.35`, `degrees = 15`).
#'
#' @return An object of class `c("yolo_workflow", "yolo_train")` containing:
#'   * `model_file`: Path to the exported `.onnx` model ready for inference.
#'   * `pt_file`: Path to the PyTorch `.pt` model weights.
#'   * `classes`: Class names.
#'   * `task`: Task type (`"detect"` or `"segment"`).
#'   * `sample_images`: Paths of the sample images used for bootstrap labeling.
#'   * `dataset_dir`: Path to the training dataset (if `keep_dataset = TRUE`).
#'   * `workflow`: Summary list of workflow parameters and metadata.
#'   * Plus all metrics, diagnostic plots, and summaries from [yolo_train()].
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Run bootstrap workflow: label 4 sample images with PerSAM, augment to 30, tile 640, train YOLO
#' wf <- yolo_train_workflow(
#'   data = "dataset_complex",
#'   sample_n = 4,
#'   model_fewshot = "persam",
#'   class_names = "cápsula",
#'   augment = TRUE,
#'   augment_total = 30,
#'   tiling = TRUE,
#'   tile_size = 640,
#'   train_model = "yolo26n",
#'   epochs = 50,
#'   output_name = "modelo_capsulas.onnx"
#' )
#'
#' # 2. Automatically propagate trained model across ALL hundreds of images in the folder!
#' ds_full <- yolo_dataset_autolabel("dataset_complex", model = wf, output_dir = "dataset_full")
yolo_train_workflow <- function(data,
                                sample_n = 4L,
                                sample_images = NULL,
                                model_fewshot = "persam",
                                class_names = "object",
                                task = c("detect", "segment"),
                                prompt = NULL,
                                prompt_mode = c("first", "each"),
                                sim_threshold = 0.70,
                                crop_to_objects = TRUE,
                                crop_margin = 0.10,
                                augment = TRUE,
                                augment_total = 30L,
                                augment_args = list(),
                                tiling = TRUE,
                                tile_size = 640L,
                                tiling_args = list(),
                                train_model = "yolo26n",
                                epochs = 50L,
                                batch = 8L,
                                device = "auto",
                                close_mosaic = 15L,
                                output_dir = pliman_model_dir(),
                                output_name = NULL,
                                temp_dir = NULL,
                                dataset_dir = NULL,
                                keep_dataset = TRUE,
                                as_job = FALSE,
                                review = FALSE,
                                show_class = FALSE,
                                show_conf = FALSE,
                                badge = FALSE,
                                interactive = is.null(prompt) && base::interactive(),
                                verbose = TRUE,
                                ...) {
  task <- match.arg(task)
  prompt_mode <- match.arg(prompt_mode)
  # 1. Resolve source unlabelled images
  if (is.character(data) && length(data) == 1L && dir.exists(data)) {
    data_dir <- normalizePath(data, winslash = "/", mustWork = TRUE)
    img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
    all_imgs <- list.files(data_dir, pattern = img_patt, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
    all_imgs <- all_imgs[!grepl("(_aug\\d+|_tile_\\d+_\\d+|^review_)", basename(all_imgs))]
    if (length(all_imgs) == 0L) {
      all_imgs <- list.files(data_dir, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      all_imgs <- all_imgs[!grepl("(_aug\\d+|_tile_\\d+_\\d+|^review_)", basename(all_imgs))]
    }
    if (length(all_imgs) == 0L) {
      cli::cli_abort("No valid images found in directory {.path {data}}.")
    }
  } else if (is.character(data) && all(file.exists(data))) {
    data_dir <- dirname(normalizePath(data[1], winslash = "/", mustWork = FALSE))
    all_imgs <- normalizePath(data, winslash = "/", mustWork = TRUE)
  } else {
    cli::cli_abort("Invalid {.arg data} argument. Must be a directory path or character vector of image files.")
  }
  # 2. Select sample images
  if (!is.null(sample_images)) {
    if (is.numeric(sample_images)) {
      n_sel <- min(length(all_imgs), as.integer(sample_images[1]))
      selected_imgs <- sample(all_imgs, n_sel)
    } else if (is.character(sample_images)) {
      all_stems <- tools::file_path_sans_ext(basename(all_imgs))
      all_bases <- basename(all_imgs)
      matched_idx <- integer(0)
      for (si in sample_images) {
        si_stem <- tools::file_path_sans_ext(basename(si))
        si_base <- basename(si)
        m <- which(all_imgs == si |
                   all_bases == si_base |
                   all_stems == si_stem |
                   tolower(all_stems) == tolower(si_stem) |
                   tolower(all_bases) == tolower(si_base))
        if (length(m) > 0L) {
          matched_idx <- c(matched_idx, m[1])
        }
      }
      matched_idx <- unique(matched_idx)
      if (length(matched_idx) == 0L) {
        cli::cli_abort(c(
          "None of the specified {.arg sample_images} were found in {.path {data_dir}}.",
          "i" = "Searched for: {.val {sample_images}}",
          "i" = "Found {length(all_imgs)} image(s) in directory (e.g. {.val {head(all_stems, 3)}})."
        ))
      }
      if (length(matched_idx) < length(sample_images)) {
        found_stems <- tolower(all_stems[matched_idx])
        missing_si <- sample_images[!tolower(tools::file_path_sans_ext(basename(sample_images))) %in% found_stems]
        if (length(missing_si) > 0L && isTRUE(verbose)) {
          cli::cli_alert_warning("{length(missing_si)} sample image(s) not found and ignored: {.val {missing_si}}")
        }
      }
      selected_imgs <- all_imgs[matched_idx]
    } else {
      cli::cli_abort("Invalid {.arg sample_images}. Must be numeric (count) or character (filenames/paths).")
    }
  } else {
    n_sel <- min(length(all_imgs), as.integer(max(1L, sample_n[1])))
    selected_imgs <- if (length(all_imgs) <= n_sel) all_imgs else sample(all_imgs, n_sel)
  }
  # 3. Setup temporary workspace
  if (is.null(temp_dir) && !is.null(dataset_dir)) {
    temp_dir <- dataset_dir
  }
  user_temp <- !is.null(temp_dir)
  if (is.null(temp_dir)) {
    d_stem <- basename(data_dir)
    if (!nzchar(d_stem) || d_stem == ".") d_stem <- "workflow"
    if (isTRUE(keep_dataset)) {
      out_base <- if (!is.null(output_dir) && dir.exists(output_dir)) output_dir else getwd()
      temp_dir <- file.path(out_base, paste0(d_stem, "_dataset"))
    } else {
      temp_dir <- file.path(tempdir(), paste0("yolo_wf_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    }
  }
  temp_dir <- normalizePath(temp_dir, winslash = "/", mustWork = FALSE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(as_job)) {
    keep_dataset <- TRUE
  }
  if (isFALSE(keep_dataset)) {
    on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  if (isTRUE(verbose)) {
    cli::cli_h1("YOLO Active Learning Workflow (pliman)")
    cli::cli_bullets(c(
      "*" = "Source image pool: {.val {length(all_imgs)}} images in {.path {data_dir}}",
      "*" = "Bootstrap sample: {.val {length(selected_imgs)}} images ({.val {paste(basename(selected_imgs), collapse = ', ')}})",
      "*" = "Few-shot annotator: {.val {model_fewshot}} (class: {.val {paste(class_names, collapse = ', ')}})",
      "*" = "Augmentation: {.val {isTRUE(augment)}} {if (isTRUE(augment)) paste0('(target: ', augment_total, ' imgs)') else ''}",
      "*" = "Tiling: {.val {isTRUE(tiling)}} {if (isTRUE(tiling)) paste0('(size: ', tile_size, 'x', tile_size, ' px)') else ''}",
      "*" = "YOLO architecture: {.val {train_model}} (epochs: {epochs}, batch: {batch})"
    ))
  }
  # Step 1: Few-shot bootstrapping
  if (isTRUE(verbose)) {
    cli::cli_h2("Step 1/4: Few-Shot Bootstrapping ({length(selected_imgs)} sample images)")
  }
  exp_args <- list(
    img = selected_imgs,
    model = model_fewshot,
    dir = temp_dir,
    task = task,
    class_names = class_names,
    prompt = prompt,
    prompt_mode = prompt_mode,
    sim_threshold = sim_threshold,
    crop_to_objects = crop_to_objects,
    crop_margin = crop_margin,
    review = review,
    show_class = show_class,
    show_conf = show_conf,
    badge = badge,
    move = FALSE,
    interactive = interactive,
    verbose = verbose
  )
  do.call(yolo_dataset_export, exp_args)
  # Step 2: Offline Augmentation
  if (isTRUE(augment)) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Step 2/4: Offline Dataset Augmentation (expanding to {augment_total} images)")
    }
    aug_params <- list(
      data = temp_dir,
      total_imgs = augment_total,
      brightness = TRUE,
      brightness_range = c(0.6, 1.2),
      contrast = TRUE,
      contrast_range = c(1, 3),
      perspective = TRUE,
      perspective_x = c(20, 40),
      perspective_y = c(20, 40),
      perspective_z = c(20, 40),
      noise = TRUE,
      noise_sd = c(0.05, 0.1),
      blur = TRUE,
      blur_sigma = c(2, 4),
      rotate = c(-45, 45),
      cutout = TRUE,
      cutout_target = "object",
      cutout_fill = "mean",
      cutout_size = c(0.4, 0.6),
      cutout_n = c(3, 5),
      max_ops = c(1, 3),
      verbose = verbose
    )
    if (is.list(augment_args) && length(augment_args) > 0L) {
      aug_params[names(augment_args)] <- augment_args
    }
    do.call(yolo_dataset_augment, aug_params)
  }
  # Step 3: Tiling
  if (isTRUE(tiling)) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Step 3/4: High-Resolution Tiling (tile size: {tile_size}x{tile_size} px)")
    }
    tile_params <- list(
      data = temp_dir,
      tile_size = tile_size,
      overlap = 0.20,
      min_coverage = 0.30,
      max_bg_ratio = 0.25,
      include_empty = FALSE,
      verbose = verbose
    )
    if (is.list(tiling_args) && length(tiling_args) > 0L) {
      tile_params[names(tiling_args)] <- tiling_args
    }
    do.call(yolo_dataset_tiling, tile_params)
  }
  # Step 4: YOLO Training
  if (isTRUE(verbose)) {
    cli::cli_h2("Step 4/4: Fast YOLO Training ({train_model})")
  }
  if (is.null(output_name) || !nzchar(output_name)) {
    data_stem <- basename(data_dir)
    if (!nzchar(data_stem) || data_stem == ".") data_stem <- "workflow"
    output_name <- paste0(data_stem, "_workflow.onnx")
  }
  train_params <- list(
    data = temp_dir,
    model = train_model,
    epochs = epochs,
    imgsz = tile_size,
    batch = batch,
    device = device,
    close_mosaic = close_mosaic,
    output_dir = output_dir,
    output_name = output_name,
    plots = FALSE,
    verbose = verbose
  )
  extra_dots <- list(...)
  if (length(extra_dots) > 0L) {
    train_params[names(extra_dots)] <- extra_dots
  }
  if (isTRUE(as_job)) {
    train_params$as_job <- TRUE
  }
  fit <- do.call(yolo_train, train_params)
  if (isTRUE(as_job)) {
    cli::cli_h2("YOLO Workflow - Training in Background")
    cli::cli_alert_success("Interactive bootstrapping and dataset generation completed!")
    cli::cli_alert_info("YOLO model training ({train_model}) was dispatched to RStudio's {.strong Jobs} tab.")
    cli::cli_alert_info("Your R console is free while the model trains. Checkpoints are automatically synced.")
    return(invisible(fit))
  }
  if (isTRUE(keep_dataset) || user_temp) {
    if (isTRUE(verbose)) cli::cli_alert_success("Bootstrap dataset preserved at: {.path {temp_dir}}")
  }
  res <- fit
  class(res) <- unique(c("yolo_workflow", class(fit)))
  res$classes <- class_names
  res$task <- task
  res$sample_images <- selected_imgs
  res$dataset_dir <- if (isTRUE(keep_dataset)) temp_dir else NULL
  res$show_class <- show_class
  res$show_conf <- show_conf
  res$badge <- badge
  res$workflow <- list(
    data = data_dir,
    model_fewshot = model_fewshot,
    sample_count = length(selected_imgs),
    augmented = isTRUE(augment),
    augment_total = if (isTRUE(augment)) augment_total else length(selected_imgs),
    tiled = isTRUE(tiling),
    tile_size = tile_size
  )
  # Ensure no unwanted 'runs' directory remains from Ultralytics
  wd_runs <- file.path(getwd(), "runs")
  if (dir.exists(wd_runs)) {
    unlink(wd_runs, recursive = TRUE, force = TRUE)
  }
  if (!is.null(output_dir)) {
    out_runs <- file.path(output_dir, "runs")
    if (dir.exists(out_runs)) {
      unlink(out_runs, recursive = TRUE, force = TRUE)
    }
  }
  if (!is.null(data_dir)) {
    d_runs <- file.path(data_dir, "runs")
    if (dir.exists(d_runs)) {
      unlink(d_runs, recursive = TRUE, force = TRUE)
    }
  }
  invisible(res)
}
#' @export
yolo_workflow <- yolo_train_workflow
#' @method print yolo_workflow
#' @export
print.yolo_workflow <- function(x, ...) {
  cli::cli_h2("YOLO Active Learning Workflow Model (pliman)")
  b_map <- if (!is.null(x$best_metrics) && "map50" %in% names(x$best_metrics)) {
    sprintf("%.4f (epoch %s)", x$best_metrics$map50[1], x$best_metrics$epoch[1])
  } else {
    "N/A"
  }
  cli::cli_bullets(c(
    "*" = "Trained model: {.file {basename(x$model_file)}} ({.path {x$model_file}})",
    "*" = "PyTorch weights: {.file {basename(x$pt_file)}} ({.path {x$pt_file}})",
    "*" = "Task: {.val {x$task}}",
    "*" = "Classes: {.val {paste(x$classes, collapse = ', ')}}",
    "*" = "Bootstrap images used: {.val {length(x$sample_images)}}",
    "*" = "Best mAP50: {.val {b_map}}",
    "v" = "Ready for batch auto-labeling!"
  ))
  cli::cli_alert_info("To auto-label all remaining images in a folder, run:")
  cli::cli_code(sprintf("ds_full <- yolo_dataset_autolabel(input_dir = '%s', model = mod, output_dir = 'dataset_auto')",
                        if (!is.null(x$workflow$data)) x$workflow$data else "my_folder"))
  invisible(x)
}
#' @method plot yolo_workflow
#' @export
plot.yolo_workflow <- function(x, ...) {
  plot.yolo_train(x, ...)
}
#' Recover and Inspect YOLO Training Checkpoints
#'
#' Scans temporary session folders, training run directories, or specified paths
#' for in-progress or interrupted YOLO training checkpoints (`last.pt` and `best.pt`).
#' Allows rescuing weights, results, and optionally exporting an in-progress checkpoint
#' to ONNX.
#'
#' @param path Character string specifying the folder or file to search or recover.
#'   If `NULL` (default), scans temporary session directories (`tempdir()` and parent)
#'   as well as the current working directory for `train_run` or `yolo_train_run` folders.
#' @param dest_dir Destination directory where recovered checkpoints and results should
#'   be copied. Defaults to `"recovered_yolo_checkpoint"`.
#' @param export_onnx Logical, whether to automatically export `best.pt` (or `last.pt`)
#'   to ONNX using the Python environment. Default is `FALSE`.
#' @param imgsz Integer image resolution for ONNX export if `export_onnx = TRUE`. Default is `640`.
#'
#' @return A list with discovered checkpoints, file paths, and metadata, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Scan for any interrupted training runs and save to ./recovered_checkpoint
#' yolo_checkpoint_recover()
#'
#' # Recover and export to ONNX directly
#' yolo_checkpoint_recover(dest_dir = "my_model", export_onnx = TRUE)
#' }
yolo_checkpoint_recover <- function(path = NULL,
                                   dest_dir = "recovered_yolo_checkpoint",
                                   export_onnx = FALSE,
                                   imgsz = 640) {
  cli::cli_h2("YOLO Checkpoint Recovery (pliman)")
  search_dirs <- character(0)
  if (!is.null(path) && nzchar(path)) {
    if (file.exists(path)) {
      if (dir.exists(path)) {
        search_dirs <- c(search_dirs, path)
      } else if (grepl("\\.pt$", path, ignore.case = TRUE)) {
        search_dirs <- c(search_dirs, dirname(path))
      }
    } else {
      cli::cli_alert_warning("Specified path does not exist: {.path {path}}")
    }
  }
  search_dirs <- c(
    search_dirs,
    file.path(tempdir(), "yolo_train_run"),
    tempdir(),
    getwd()
  )
  # Also search parent tempdir for other R sessions
  parent_temp <- dirname(tempdir())
  if (dir.exists(parent_temp)) {
    sub_rtmps <- list.dirs(parent_temp, recursive = FALSE, full.names = TRUE)
    rtmp_trains <- file.path(sub_rtmps[grepl("Rtmp", basename(sub_rtmps))], "yolo_train_run")
    search_dirs <- c(search_dirs, rtmp_trains[dir.exists(rtmp_trains)])
  }
  search_dirs <- unique(search_dirs[dir.exists(search_dirs)])
  pt_files <- character(0)
  for (sd in search_dirs) {
    found <- tryCatch({
      list.files(sd, pattern = "\\.pt$", recursive = TRUE, full.names = TRUE)
    }, error = function(e) character(0))
    pt_files <- c(pt_files, found)
  }
  pt_files <- unique(normalizePath(pt_files, winslash = "/", mustWork = FALSE))
  pt_files <- pt_files[file.exists(pt_files)]
  if (length(pt_files) == 0L) {
    cli::cli_alert_info("No YOLO checkpoint weights (`.pt`) found in scanned directories.")
    return(invisible(list(recovered = character(0))))
  }
  # Sort by modification time (most recent first)
  mtimes <- file.info(pt_files)$mtime
  pt_files <- pt_files[order(mtimes, decreasing = TRUE)]
  # Identify best and last candidates
  cand_best <- pt_files[grepl("best\\.pt$", basename(pt_files))]
  cand_last <- pt_files[grepl("last\\.pt$", basename(pt_files))]
  primary_best <- if (length(cand_best) > 0) cand_best[1] else NULL
  primary_last <- if (length(cand_last) > 0) cand_last[1] else NULL
  cli::cli_alert_success("Found {length(pt_files)} checkpoint weight file{?s}:")
  for (f in pt_files) {
    sz_mb <- round(file.info(f)$size / (1024 * 1024), 1)
    mt <- format(file.info(f)$mtime, "%Y-%m-%d %H:%M:%S")
    cli::cli_ul(c("{.file {basename(f)}} ({sz_mb} MB, modified {mt}): {.path {f}}"))
  }
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dest_weights_dir <- file.path(dest_dir, "weights")
  if (!dir.exists(dest_weights_dir)) {
    dir.create(dest_weights_dir, recursive = TRUE, showWarnings = FALSE)
  }
  recovered_paths <- character(0)
  for (f in c(primary_best, primary_last)) {
    if (!is.null(f) && file.exists(f)) {
      target <- file.path(dest_weights_dir, basename(f))
      file.copy(f, target, overwrite = TRUE)
      recovered_paths <- c(recovered_paths, target)
    }
  }
  # Look for results.csv to recover metrics
  results_csvs <- character(0)
  for (sd in search_dirs) {
    found_csv <- tryCatch({
      list.files(sd, pattern = "^results\\.csv$", recursive = TRUE, full.names = TRUE)
    }, error = function(e) character(0))
    results_csvs <- c(results_csvs, found_csv)
  }
  latest_metrics <- NULL
  if (length(results_csvs) > 0L) {
    mtimes_csv <- file.info(results_csvs)$mtime
    latest_csv <- results_csvs[order(mtimes_csv, decreasing = TRUE)][1]
    file.copy(latest_csv, file.path(dest_dir, "results.csv"), overwrite = TRUE)
    recovered_paths <- c(recovered_paths, file.path(dest_dir, "results.csv"))
    csv_data <- tryCatch(utils::read.csv(latest_csv, check.names = FALSE, strip.white = TRUE), error = function(e) NULL)
    if (!is.null(csv_data) && nrow(csv_data) > 0L) {
      latest_metrics <- csv_data[nrow(csv_data), , drop = FALSE]
      ep_col <- names(latest_metrics)[grepl("epoch", names(latest_metrics), ignore.case = TRUE)]
      curr_ep <- if (length(ep_col) > 0) latest_metrics[[ep_col[1]]] else nrow(csv_data)
      cli::cli_alert_info("Latest completed epoch recorded in results.csv: {.val {curr_ep}}")
    }
  }
  # Export ONNX if requested
  exported_onnx_path <- NULL
  if (isTRUE(export_onnx)) {
    model_to_export <- if (!is.null(primary_best)) primary_best else primary_last
    if (!is.null(model_to_export)) {
      cli::cli_alert_info("Exporting recovered checkpoint to ONNX...")
      py_exec <- .find_python_exec()
      if (length(py_exec) > 0L && nzchar(py_exec[1])) {
        out_onnx_target <- file.path(dest_dir, paste0(sub("\\.pt$", "", basename(model_to_export)), ".onnx"))
        export_py <- paste0(
          "from ultralytics import YOLO\n",
          "model = YOLO(r'", normalizePath(model_to_export, winslash = "/", mustWork = FALSE), "')\n",
          "model.export(format='onnx', imgsz=", as.integer(imgsz), ", dynamic=True)\n"
        )
        t_py <- file.path(tempdir(), "temp_recover_export.py")
        writeLines(export_py, t_py)
        res_exp <- system2(py_exec, args = shQuote(t_py))
        src_ox <- sub("\\.pt$", ".onnx", model_to_export)
        if (file.exists(src_ox)) {
          file.copy(src_ox, out_onnx_target, overwrite = TRUE)
          exported_onnx_path <- out_onnx_target
          recovered_paths <- c(recovered_paths, exported_onnx_path)
          cli::cli_alert_success("ONNX model exported: {.path {exported_onnx_path}}")
        }
      }
    }
  }
  cli::cli_alert_success("Checkpoints safely recovered to {.path {dest_dir}}")
  if (!is.null(primary_last)) {
    cli::cli_alert_info("To resume training from this checkpoint in pliman, run:")
    cli::cli_code(paste0("yolo_train(..., resume = \"", normalizePath(primary_last, winslash = "/", mustWork = FALSE), "\")"))
  }
  invisible(list(
    best = primary_best,
    last = primary_last,
    all_checkpoints = pt_files,
    dest_dir = dest_dir,
    recovered_files = recovered_paths,
    metrics = latest_metrics,
    onnx = exported_onnx_path
  ))
}
