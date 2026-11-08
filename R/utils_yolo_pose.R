# ==============================================================================
# YOLO-POSE DATASET INTERACTIVE LABELING & ANNOTATION ENGINE
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Interactive Dataset Annotation for YOLO-Pose (Keypoints & Bounding Boxes)
#' @name yolo_dataset_label_pose
#' @aliases yolo_pose_label yolo_dataset_pose
#' @description
#' Provides an ultra-efficient, ergonomic interactive labeling tool for creating and editing
#' YOLO-Pose datasets. Supports three essential labeling strategies:
#' 1. **Auto BBox from Points** (`box_mode = "auto"`): Click *only* the K keypoints; the bounding box
#'    is automatically computed with a configurable padding margin (`box_pad`). Cuts clicks in half!
#' 2. **Manual BBox + Points** (`box_mode = "manual"`): Click 2 opposite corners for the bounding box,
#'    followed by the K keypoints in sequence.
#' 3. **Visibility & Occlusion Control** (`allow_occlusion = TRUE`): Supports marking points as visible (`v=2`),
#'    occluded (`v=1`), or absent/skipped (`v=0`) with instantaneous keyboard/right-click shortcuts.
#'
#' @param dir Character or `yolo_dataset` object. Dataset root directory (containing `images/` and `labels/` or a flat directory).
#'   Default is `"."`.
#' @param image Optional character vector. Specific image name(s) or file paths to annotate. If `NULL` (default),
#'   iterates through all images in the dataset.
#' @param kpts Integer. Number of keypoints per object (e.g. `2` for log/stick endpoints, `4` for corners/seeds,
#'   `5` for facial landmarks, `17` for human pose). Default is `2L`.
#' @param kpt_names Character vector. Optional human-readable names for each keypoint (e.g. `c("top", "bottom")` or `c("start", "end")`).
#'   These names are displayed on the prompt and badge as each point is clicked.
#' @param skeleton Optional list or matrix. Pairs of 1-based keypoint indices to connect with skeleton lines
#'   (e.g., `list(c(1, 2))` or `matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE)`).
#'   If `"chain"`, `TRUE`, or `NULL` with `kpts >= 2`, points are connected sequentially (`1-2-3...`).
#'   Set to `FALSE` or `"none"` to disable skeleton lines.
#' @param box_mode Character. How bounding boxes are defined:
#'   - `"auto"` (default, Strategy 2): User clicks ONLY the K keypoints. Bounding box is auto-generated enclosing all points + `box_pad`.
#'   - `"manual"` (Strategy 1): User clicks 2 opposite corners to define the bounding box, then clicks the K keypoints.
#'   - `"prompt"`: Prompts the user per object whether to click a manual box or use auto-box.
#' @param box_pad Numeric or numeric vector of length 2. Margin (relative to box dimension or image size)
#'   to expand the bounding box when `box_mode = "auto"`. Default is `0.10` (10% padding).
#'   Ensures non-zero width/height even when points are strictly collinear (e.g. 2 endpoints of a line).
#' @param allow_occlusion Logical. If `TRUE` (default, Strategy 3), allows skipping absent keypoints (`v = 0`)
#'   or marking them as occluded (`v = 1`) via keyboard/right-click shortcuts (`s` for skip/absent, `o` for occluded).
#' @param class_names Character vector. Names of object classes (default: `"object"`).
#' @param class_id Integer. Default class index for newly added objects (default: `0L`).
#' @param split Character. Subset to annotate when YOLO folder structure is detected (`"all"`, `"train"`, `"val"`, or `"test"`). Default is `"all"`.
#' @param resume Logical. If `TRUE` (default), starts from the first unannotated image in the dataset.
#' @param update_yaml Logical. If `TRUE` (default), automatically writes/updates `kpt_shape: [kpts, 3]`
#'   and `skeleton` in `data.yaml`.
#' @param backup Logical. If `TRUE` (default), creates a backup (`<name>.txt.bak`) before saving modifications.
#' @param rainbow Logical. If `TRUE` (default), assigns distinct, harmonious colors to each keypoint.
#' @param pt_cex Numeric. Size of the keypoint circles in the plot (default: `1.8`).
#' @param line_lwd Numeric. Line width for skeleton connections and bounding box borders (default: `2.0`).
#' @param verbose Logical. If `TRUE` (default), prints progress and statistics.
#'
#' @return An invisible list containing:
#'   \item{total_images}{Total number of images inspected}
#'   \item{annotated_images}{Number of images with annotations}
#'   \item{total_objects}{Total number of pose objects labeled}
#'   \item{kpts}{Number of keypoints configured}
#'   \item{labels_dir}{Path to the labels directory}
#'
#' @export
#' @examples
#' \dontrun{
#' # Label wood logs with 2 keypoints (extremities) using auto-bounding box (only 2 clicks per log!)
#' yolo_dataset_label_pose("wood_log", kpts = 2, kpt_names = c("end1", "end2"), box_mode = "auto")
#'
#' # Label plants with 4 keypoints (base, node1, node2, apex) using manual box + occlusion control
#' yolo_dataset_label_pose("plants", kpts = 4, kpt_names = c("base", "n1", "n2", "apex"), box_mode = "manual")
#' }
yolo_dataset_label_pose <- function(dir = ".",
                                    image = NULL,
                                    kpts = 2L,
                                    kpt_names = NULL,
                                    skeleton = NULL,
                                    box_mode = c("auto", "manual", "prompt"),
                                    box_pad = 0.10,
                                    allow_occlusion = TRUE,
                                    class_names = NULL,
                                    class_id = 0L,
                                    split = c("all", "train", "val", "test"),
                                    resume = TRUE,
                                    update_yaml = TRUE,
                                    backup = TRUE,
                                    split_data = FALSE,
                                    split_ratio = c(train = 0.70, val = 0.20, test = 0.10),
                                    rainbow = TRUE,
                                    pt_cex = 1.8,
                                    line_lwd = 2.0,
                                    verbose = TRUE) {
  if (!base::interactive()) {
    cli::cli_abort("yolo_dataset_label_pose() requires an interactive R session with graphics display.")
  }

  box_mode <- match.arg(box_mode)
  split <- match.arg(split)
  kpts <- as.integer(max(1L, kpts[1]))

  # Default keypoint names if not supplied
  if (is.null(kpt_names) || length(kpt_names) == 0L) {
    if (kpts == 2L) {
      kpt_names <- c("p1", "p2")
    } else if (kpts == 4L) {
      kpt_names <- c("p1", "p2", "p3", "p4")
    } else {
      kpt_names <- paste0("kpt_", seq_len(kpts))
    }
  } else {
    kpt_names <- as.character(kpt_names)
    if (length(kpt_names) < kpts) {
      kpt_names <- c(kpt_names, paste0("kpt_", seq(length(kpt_names) + 1L, kpts)))
    }
  }

  # Resolve skeleton pairs (1-indexed)
  skeleton_mat <- .yolo_resolve_skeleton(skeleton, kpts)

  # Resolve directory
  dir_path <- if (inherits(dir, "yolo_dataset") || (is.list(dir) && !is.null(dir$dir))) {
    dir$dir
  } else if (is.character(dir) && length(dir) == 1L) {
    dir
  } else {
    cli::cli_abort("{.arg dir} must be a single directory path or a {.code yolo_dataset} object.")
  }
  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")
  }

  # Resolve class names
  if (is.null(class_names)) {
    cl_found <- tryCatch(.yolo_get_classes(dir_path), error = function(e) NULL)
    class_names <- if (!is.null(cl_found) && length(cl_found) > 0L) cl_found else "object"
  } else {
    class_names <- as.character(class_names)
  }
  class_id <- as.integer(class_id[1])

  # Locate images and corresponding labels
  img_ext_patt <- "(?i)\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  has_split_folders <- dir.exists(file.path(dir_path, "images")) || dir.exists(file.path(dir_path, "train"))

  img_files <- character(0)
  if (has_split_folders) {
    search_splits <- if (split == "all") c("train", "val", "test") else split
    for (sp in search_splits) {
      sp_dir <- file.path(dir_path, "images", sp)
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp, "images")
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp)
      if (dir.exists(sp_dir)) {
        fls <- list.files(sp_dir, pattern = img_ext_patt, full.names = TRUE, recursive = TRUE)
        img_files <- c(img_files, fls)
      }
    }
  }

  if (length(img_files) == 0L) {
    # Scan root directory directly
    img_files <- list.files(dir_path, pattern = img_ext_patt, full.names = TRUE, recursive = TRUE)
  }

  # Filter temporary files
  img_files <- img_files[!grepl("[/\\\\]\\.__tmp_", img_files)]

  if (length(img_files) == 0L) {
    cli::cli_abort("No image files found in dataset {.path {dir_path}}.")
  }

  # Filter specific images if requested
  if (!is.null(image) && length(image) > 0L) {
    img_targets <- tolower(basename(as.character(image)))
    stem_targets <- tolower(tools::file_path_sans_ext(img_targets))
    all_bases <- tolower(basename(img_files))
    all_stems <- tolower(tools::file_path_sans_ext(all_bases))
    keep_mask <- all_bases %in% img_targets | all_stems %in% stem_targets
    if (any(keep_mask)) {
      img_files <- img_files[keep_mask]
    } else {
      cli::cli_abort("None of the specified images were found in the dataset.")
    }
  }

  total_imgs <- length(img_files)

  # Synchronize data.yaml with kpt_shape and skeleton if requested
  if (isTRUE(update_yaml)) {
    .yolo_sync_pose_yaml(dir_path, kpts = kpts, skeleton = skeleton_mat, class_names = class_names)
  }

  # Palette for keypoints: fixed harmonious colors per keypoint index
  kpt_colors <- if (isTRUE(rainbow)) {
    grDevices::rainbow(kpts, s = 0.90, v = 0.95)
  } else {
    rep("#00FF66", kpts)
  }

  # Find starting index if resume = TRUE
  start_idx <- 1L
  if (isTRUE(resume)) {
    for (ri in seq_along(img_files)) {
      lbl_cand <- .yolo_find_matching_label(img_files[ri], dir_path)
      if (!file.exists(lbl_cand)) {
        start_idx <- ri
        break
      }
      lines_cand <- readLines(lbl_cand, warn = FALSE)
      lines_cand <- lines_cand[nzchar(trimws(lines_cand))]
      if (length(lines_cand) == 0L) {
        start_idx <- ri
        break
      }
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO-Pose Interactive Dataset Annotator")
    cli::cli_bullets(c(
      "*" = "Dataset: {.path {dir_path}}",
      "*" = "Images: {total_imgs} images found (starting at #{start_idx})",
      "*" = "Keypoints: {kpts} point{?s} per object ({paste(kpt_names, collapse = ', ')})",
      "*" = "Box strategy: {.val {toupper(box_mode)}} (Auto-box pad: {box_pad * 100}%)",
      "*" = "Occlusion control: {.val {if (allow_occlusion) 'Enabled (Right-click / Esc to flag)' else 'Disabled'}}"
    ))
  }

  # Main navigation loop across images
  curr_idx <- start_idx
  total_annotated <- 0L

  while (curr_idx >= 1L && curr_idx <= total_imgs) {
    target_img_path <- img_files[curr_idx]
    base_stem <- tools::file_path_sans_ext(basename(target_img_path))
    dest_lbl_path <- .yolo_find_matching_label(target_img_path, dir_path)

    # Ensure label directory exists
    lbl_dir <- dirname(dest_lbl_path)
    if (!dir.exists(lbl_dir)) dir.create(lbl_dir, recursive = TRUE, showWarnings = FALSE)

    # Backup if requested
    if (isTRUE(backup) && file.exists(dest_lbl_path) && !file.exists(paste0(dest_lbl_path, ".bak"))) {
      file.copy(dest_lbl_path, paste0(dest_lbl_path, ".bak"), overwrite = FALSE)
    }

    # Load image
    curr_img <- tryCatch(image_import(target_img_path), error = function(e) NULL)
    if (is.null(curr_img)) {
      cli::cli_alert_warning("Could not read image {.path {target_img_path}}. Skipping.")
      curr_idx <- curr_idx + 1L
      next
    }

    img_dims <- dim(curr_img)
    img_w <- img_dims[1]
    img_h <- img_dims[2]

    # Read existing pose annotations
    parsed_lbl <- .yolo_read_pose_label(dest_lbl_path, img_w = img_w, img_h = img_h, kpts = kpts, class_names = class_names)
    boxes <- parsed_lbl$boxes
    keypoints <- parsed_lbl$keypoints

    # Current image editing loop
    save_on_exit <- TRUE
    nav_action <- "next"

    # Initial render
    .yolo_render_pose_canvas(
      img = curr_img,
      boxes = boxes,
      keypoints = keypoints,
      kpt_names = kpt_names,
      skeleton = skeleton_mat,
      img_name = base_stem,
      curr_idx = curr_idx,
      total_imgs = total_imgs,
      box_mode = box_mode,
      kpt_colors = kpt_colors,
      pt_cex = pt_cex,
      line_lwd = line_lwd
    )

    repeat {
      n_objs <- nrow(boxes)
      cat("\n")
      cli::cli_h3(sprintf("YOLO-Pose [%s | %d/%d | %d obj%s | K=%d | Mode: %s]",
                          base_stem, curr_idx, total_imgs, n_objs, if (n_objs == 1) "" else "s", kpts, toupper(box_mode)))
      cli::cli_bullets(c(
        "*" = paste0(cli::col_green("[Enter] / c / n"), ": Save & go to NEXT image"),
        "*" = paste0(cli::col_green("p"), ": Save & go to PREVIOUS image"),
        "*" = paste0(cli::col_green("a"), ": Add a new pose object (", kpts, " points)"),
        "*" = paste0(cli::col_green("aa / cont"), ": Continuous rapid annotation mode"),
        "*" = paste0(cli::col_green("m"), ": Toggle box mode (auto-bbox <-> manual box)"),
        "*" = paste0(cli::col_green("e <idx>"), ": Re-annotate / edit points of object <idx>"),
        "*" = paste0(cli::col_green("d <idx>"), ": Delete object <idx>"),
        "*" = paste0(cli::col_green("r"), ": Click on an object to delete it"),
        "*" = paste0(cli::col_green("clear"), ": Clear all objects (mark as negative/background)"),
        "*" = paste0(cli::col_green("s / skip"), ": Skip image without saving changes"),
        "*" = paste0(cli::col_red("rm"), ": Permanently delete image and label from dataset"),
        "*" = paste0(cli::col_green("q / finish"), ": Save and exit session")
      ))

      prompt_txt <- sprintf("Select action [Enter = Save & Next (%d/%d)]: ", curr_idx, total_imgs)
      ans <- tryCatch(tolower(trimws(readline(prompt_txt))), interrupt = function(e) "q")

      # Navigation: Save & Next
      if (ans == "" || ans %in% c("c", "continue", "n", "next", "save", "ok", "y")) {
        nav_action <- "next"
        break
      }

      # Navigation: Save & Prev
      if (ans %in% c("p", "prev", "previous")) {
        nav_action <- "prev"
        break
      }

      # Skip image
      if (ans %in% c("s", "skip")) {
        save_on_exit <- FALSE
        nav_action <- "next"
        cli::cli_alert_info("Skipped {.file {base_stem}} without saving.")
        break
      }

      # Finish session
      if (ans %in% c("q", "quit", "finish", "done", "exit")) {
        nav_action <- "exit"
        break
      }

      # Abort session
      if (ans %in% c("abort", "cancel")) {
        save_on_exit <- FALSE
        nav_action <- "exit"
        cli::cli_alert_info("Exited session without saving {.file {base_stem}}.")
        break
      }

      # Toggle box mode
      if (ans %in% c("m", "mode", "toggle")) {
        box_mode <- if (box_mode == "auto") "manual" else "auto"
        cli::cli_alert_success("Switched box mode to: {.val {toupper(box_mode)}}.")
        .yolo_render_pose_canvas(
          img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
          skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
          box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
        )
        next
      }

      # Permanently remove image and label
      if (ans %in% c("rm", "delete_image", "del_img", "remove_image")) {
        unlink(target_img_path, force = TRUE)
        unlink(dest_lbl_path, force = TRUE)
        unlink(paste0(dest_lbl_path, ".bak"), force = TRUE)
        cli::cli_alert_success("Permanently removed {.file {basename(target_img_path)}} from dataset.")
        save_on_exit <- FALSE
        nav_action <- "next"
        break
      }

      # Clear all objects
      if (ans %in% c("clear", "empty", "reset")) {
        boxes <- data.frame(
          index = integer(0), xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
          cx = numeric(0), cy = numeric(0), class_id = integer(0), class_name = character(0), stringsAsFactors = FALSE
        )
        keypoints <- list()
        cli::cli_alert_success("Cleared all annotations for {.file {base_stem}} (negative background).")
        .yolo_render_pose_canvas(
          img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
          skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
          box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
        )
        next
      }

      # Delete specific object by index (e.g. "d 2")
      if (grepl("^(d|del|delete|remove)\\s+\\d+", ans)) {
        del_num <- as.integer(sub("^(d|del|delete|remove)\\s+", "", ans))
        if (del_num >= 1L && del_num <= nrow(boxes)) {
          boxes <- boxes[-del_num, , drop = FALSE]
          keypoints <- keypoints[-del_num]
          if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
          cli::cli_alert_success("Deleted object #{del_num}. {nrow(boxes)} remaining.")
          .yolo_render_pose_canvas(
            img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
            skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
            box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
          )
        } else {
          cli::cli_alert_warning("Invalid object index #{del_num}.")
        }
        next
      }

      # Click to delete object
      if (ans %in% c("r", "rem")) {
        if (nrow(boxes) == 0L) {
          cli::cli_alert_info("No objects to delete.")
          next
        }
        cli::cli_alert_info("Click on any keypoint or inside a bounding box to delete it (<Esc> to cancel):")
        click_del <- tryCatch(graphics::locator(n = 1, type = "p", pch = 4, col = "#FF0033", cex = 1.8, lwd = 2.5), error = function(e) NULL)
        if (!is.null(click_del) && length(click_del$x) > 0L) {
          cx <- click_del$x[1]; cy <- click_del$y[1]
          best_idx <- NULL
          best_d <- Inf
          for (bi in seq_len(nrow(boxes))) {
            bx <- boxes[bi, ]
            inside <- (cx >= bx$xmin && cx <= bx$xmax && cy >= bx$ymin && cy <= bx$ymax)
            d_mid <- (bx$cx - cx)^2 + (bx$cy - cy)^2
            if (inside && d_mid < best_d) {
              best_d <- d_mid
              best_idx <- bi
            }
          }
          if (!is.null(best_idx)) {
            boxes <- boxes[-best_idx, , drop = FALSE]
            keypoints <- keypoints[-best_idx]
            if (nrow(boxes) > 0L) boxes$index <- seq_len(nrow(boxes))
            cli::cli_alert_success("Deleted object #{best_idx}. {nrow(boxes)} remaining.")
          } else {
            cli::cli_alert_warning("No object was found near clicked position.")
          }
          .yolo_render_pose_canvas(
            img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
            skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
            box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
          )
        }
        next
      }

      # Add new pose object(s)
      if (ans %in% c("a", "add", "aa", "cont", "continuous")) {
        is_continuous <- ans %in% c("aa", "cont", "continuous")
        repeat {
          res_add <- .yolo_annotate_single_pose(
            img_w = img_w,
            img_h = img_h,
            kpts = kpts,
            kpt_names = kpt_names,
            skeleton = skeleton_mat,
            box_mode = box_mode,
            box_pad = box_pad,
            allow_occlusion = allow_occlusion,
            class_id = class_id,
            class_name = class_names[1],
            kpt_colors = kpt_colors,
            pt_cex = pt_cex,
            line_lwd = line_lwd
          )

          if (!is.null(res_add)) {
            new_box <- res_add$box
            new_box$index <- nrow(boxes) + 1L
            boxes <- rbind(boxes, new_box)
            boxes$index <- seq_len(nrow(boxes))
            keypoints[[length(keypoints) + 1L]] <- res_add$kpts
            cli::cli_alert_success("Added pose object #{nrow(boxes)} with {kpts} keypoint{?s}.")

            .yolo_render_pose_canvas(
              img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
              skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
              box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
            )
          } else {
            break
          }

          if (!is_continuous) break
        }
        next
      }

      # Re-edit an existing object's keypoints
      if (grepl("^(e|edit)\\s+\\d+", ans)) {
        e_num <- as.integer(sub("^(e|edit)\\s+", "", ans))
        if (e_num >= 1L && e_num <= nrow(boxes)) {
          cli::cli_alert_info("Re-annotating keypoints for object #{e_num}...")
          res_re <- .yolo_annotate_single_pose(
            img_w = img_w,
            img_h = img_h,
            kpts = kpts,
            kpt_names = kpt_names,
            skeleton = skeleton_mat,
            box_mode = box_mode,
            box_pad = box_pad,
            allow_occlusion = allow_occlusion,
            class_id = boxes$class_id[e_num],
            class_name = boxes$class_name[e_num],
            kpt_colors = kpt_colors,
            pt_cex = pt_cex,
            line_lwd = line_lwd
          )
          if (!is.null(res_re)) {
            boxes[e_num, ] <- res_re$box
            boxes$index[e_num] <- e_num
            keypoints[[e_num]] <- res_re$kpts
            cli::cli_alert_success("Updated object #{e_num}.")
            .yolo_render_pose_canvas(
              img = curr_img, boxes = boxes, keypoints = keypoints, kpt_names = kpt_names,
              skeleton = skeleton_mat, img_name = base_stem, curr_idx = curr_idx, total_imgs = total_imgs,
              box_mode = box_mode, kpt_colors = kpt_colors, pt_cex = pt_cex, line_lwd = line_lwd
            )
          }
        } else {
          cli::cli_alert_warning("Invalid object index #{e_num}.")
        }
        next
      }

      # Jump to specific image index (e.g. "goto 5" or "5")
      if (grepl("^(goto|jump)\\s+\\d+", ans) || (grepl("^\\d+$", ans) && as.integer(ans) <= total_imgs)) {
        t_idx <- as.integer(sub("^(goto|jump)\\s+", "", ans))
        if (t_idx >= 1L && t_idx <= total_imgs) {
          curr_idx <- t_idx
          save_on_exit <- TRUE
          nav_action <- "goto"
          break
        }
      }

      cli::cli_alert_warning("Unknown command. Type {.code a} to add object, {.code [Enter]} to save & next, or {.code q} to exit.")
    }

    # Save label if modified
    if (isTRUE(save_on_exit)) {
      .yolo_write_pose_label(dest_lbl_path, boxes = boxes, keypoints = keypoints, img_w = img_w, img_h = img_h, kpts = kpts)
      if (nrow(boxes) > 0L) total_annotated <- total_annotated + 1L
    }

    if (nav_action == "exit") {
      break
    } else if (nav_action == "prev") {
      curr_idx <- max(1L, curr_idx - 1L)
    } else if (nav_action == "next") {
      curr_idx <- curr_idx + 1L
    }
  }

  # In-place dataset partitioning if requested
  if (isTRUE(split_data) && !has_split_folders && total_annotated > 0L) {
    if (isTRUE(verbose)) cli::cli_alert_info("Organizing dataset into train/val/test split folders...")
    .yolo_split_pose_dataset_inplace(dir_path, split_ratio = split_ratio, verbose = verbose)
    has_split_folders <- TRUE
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Pose labeling session completed!")
    cli::cli_bullets(c(
      "*" = "Images annotated: {total_annotated} of {total_imgs}",
      "*" = "Keypoint format: YOLO-Pose ({kpts} keypoints per object)",
      "*" = "YAML updated: {.path {file.path(dir_path, 'data.yaml')}}",
      if (!has_split_folders) c("i" = paste0("Dataset is in a flat directory. To organize into train/val/test splits for YOLO training, run: yolo_dataset_split('", dir_path, "') or pass split_data = TRUE.")) else NULL
    ))
  }

  invisible(list(
    total_images = total_imgs,
    annotated_images = total_annotated,
    kpts = kpts,
    labels_dir = if (has_split_folders) file.path(dir_path, "labels") else dir_path
  ))
}

#' @rdname yolo_dataset_label_pose
#' @export
yolo_pose_label <- yolo_dataset_label_pose

#' @rdname yolo_dataset_label_pose
#' @export
yolo_dataset_pose <- yolo_dataset_label_pose


# ==============================================================================
# YOLO-POSE DATASET PREVIEW (yolo_dataset_pose_preview)
# ==============================================================================

#' @title Preview YOLO-Pose Dataset Annotations with Skeletons and Keypoints
#' @name yolo_dataset_pose_preview
#' @aliases yolo_pose_preview
#' @description
#' Visualizes sample images from a YOLO-Pose dataset with their bounding boxes,
#' keypoints, and skeleton lines. Allows inspecting pose annotations, occlusion flags,
#' and multi-panel layouts.
#'
#' @param dir Character or `yolo_dataset` object. Path to the dataset directory. Defaults to `"."`.
#' @param split Character. Dataset split to preview (`"train"`, `"val"`, `"test"`, or `"all"`). Defaults to `"train"`.
#' @param show_bbox Logical. Whether to draw bounding boxes around objects. Defaults to `TRUE`.
#' @param show_skeleton Logical. Whether to connect keypoints with skeleton lines. Defaults to `TRUE`.
#' @param show_points Logical. Whether to draw keypoint dots. Defaults to `TRUE`.
#' @param show_labels Logical. Whether to display class label names or object index. Defaults to `TRUE`.
#' @param show_kpt_labels Logical. Whether to display keypoint index or names next to points. Defaults to `FALSE`.
#' @param kpts Integer or `NULL`. Number of keypoints. If `NULL` (default), extracted from `data.yaml` or first label.
#' @param kpt_names Character vector. Optional human-readable names for each keypoint.
#' @param skeleton Optional skeleton connection pairs. If `NULL`, reads from `data.yaml` or chains points sequentially.
#' @param n Integer. Number of images to plot. Defaults to `4L`.
#' @param which Optional integer vector of image indices or character vector of image names/stems.
#' @param images Alias for `which`.
#' @param pattern Optional regex pattern to filter image filenames.
#' @param random Logical. Whether to sample images randomly. Defaults to `TRUE` when `which` is `NULL`.
#' @param seed Optional integer seed for reproducible random selection.
#' @param mfrow Numeric vector of length 2 giving `c(rows, cols)` for the multi-panel layout.
#' @param rainbow Logical. If `TRUE` (default), assigns distinct harmonious colors to each keypoint.
#' @param col Custom border color for bounding boxes.
#' @param pt_cex Numeric. Size of keypoint dots (default: `1.5`).
#' @param line_lwd Numeric. Line width for bounding boxes and skeleton lines (default: `2.0`).
#' @param title Logical or character. Title for each panel. Defaults to `TRUE`.
#' @param save Logical. Whether to save the preview to an image file. Defaults to `FALSE`.
#' @param file_name Character. File name to save preview if `save = TRUE`.
#' @param width,height,dpi Numeric. Dimensions and DPI for saved preview image.
#' @param ... Additional arguments passed to graphics methods.
#'
#' @return An invisible list of the sampled image file paths and parsed annotations.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Preview 4 random images with pose annotations from a dataset
#' yolo_dataset_pose_preview("milho", n = 4)
yolo_dataset_pose_preview <- function(dir = ".",
                                      split = c("train", "val", "test", "all"),
                                      show_bbox = TRUE,
                                      show_skeleton = TRUE,
                                      show_points = TRUE,
                                      show_labels = TRUE,
                                      show_kpt_labels = FALSE,
                                      kpts = NULL,
                                      kpt_names = NULL,
                                      skeleton = NULL,
                                      n = 4L,
                                      which = NULL,
                                      images = NULL,
                                      pattern = NULL,
                                      random = is.null(which) && is.null(images),
                                      seed = NULL,
                                      mfrow = NULL,
                                      rainbow = TRUE,
                                      col = NULL,
                                      pt_cex = 1.5,
                                      line_lwd = 2.0,
                                      title = TRUE,
                                      save = FALSE,
                                      file_name = NULL,
                                      width = NULL,
                                      height = NULL,
                                      dpi = 300,
                                      ...) {
  dots <- list(...)
  if (!is.null(dots$bbox)) show_bbox <- isTRUE(dots$bbox)
  if (!is.null(dots$points)) show_points <- isTRUE(dots$points)
  if (!is.null(dots$labels)) show_labels <- isTRUE(dots$labels)
  if (!is.null(dots$skeleton)) show_skeleton <- isTRUE(dots$skeleton)
  if (!is.null(images) && is.null(which)) which <- images

  split <- match.arg(split)

  dir_path <- dir
  if (is.list(dir)) {
    if (!is.null(dir$dir)) dir_path <- dir$dir[1]
    else if (!is.null(dir$yaml_file)) dir_path <- dirname(dir$yaml_file[1])
    else if (length(dir) > 0L && is.character(dir[[1]])) dir_path <- dir[[1]]
  }
  dir_path <- normalizePath(dir_path[1], winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")
  }

  # Read metadata from data.yaml
  yaml_meta <- .yolo_read_pose_yaml_meta(dir_path)
  if (is.null(kpts)) kpts <- yaml_meta$kpts
  if (is.null(skeleton)) skeleton <- yaml_meta$skeleton
  class_names <- yaml_meta$classes
  if (is.null(class_names) || length(class_names) == 0L) class_names <- "object"

  # Discover image files
  img_ext_patt <- "(?i)\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  has_split_folders <- dir.exists(file.path(dir_path, "images")) || dir.exists(file.path(dir_path, "train"))

  img_files <- character(0)
  if (has_split_folders) {
    search_splits <- if (split == "all") c("train", "val", "test") else split
    for (sp in search_splits) {
      sp_dir <- file.path(dir_path, "images", sp)
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp, "images")
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp)
      if (dir.exists(sp_dir)) {
        fls <- list.files(sp_dir, pattern = img_ext_patt, full.names = TRUE, recursive = TRUE)
        img_files <- c(img_files, fls)
      }
    }
  }

  if (length(img_files) == 0L) {
    img_files <- list.files(dir_path, pattern = img_ext_patt, full.names = TRUE, recursive = FALSE)
  }

  img_files <- img_files[!grepl("[/\\\\]\\.__tmp_", img_files)]

  if (!is.null(pattern) && nzchar(pattern)) {
    img_files <- img_files[grepl(pattern, basename(img_files), ignore.case = TRUE)]
  }

  if (length(img_files) == 0L) {
    cli::cli_abort("No image files found in {.path {dir_path}} for split {.val {split}}.")
  }

  # Subsetting
  if (!is.null(which)) {
    if (is.numeric(which)) {
      w_idx <- which[which >= 1L & which <= length(img_files)]
      if (length(w_idx) > 0L) img_files <- img_files[w_idx]
    } else if (is.character(which)) {
      w_stems <- tolower(tools::file_path_sans_ext(basename(which)))
      all_stems <- tolower(tools::file_path_sans_ext(basename(img_files)))
      m_idx <- which(all_stems %in% w_stems | tolower(basename(img_files)) %in% tolower(basename(which)))
      if (length(m_idx) > 0L) img_files <- img_files[m_idx]
    }
  } else {
    n_sample <- min(as.integer(n[1]), length(img_files))
    if (isTRUE(random)) {
      if (!is.null(seed)) set.seed(seed)
      img_files <- sample(img_files, n_sample)
    } else {
      img_files <- img_files[seq_len(n_sample)]
    }
  }

  total_to_plot <- length(img_files)
  if (total_to_plot == 0L) {
    cli::cli_abort("No matching images to preview.")
  }

  # Inspect first label to infer kpts if still NULL
  if (is.null(kpts) || is.na(kpts) || kpts < 1L) {
    for (im_f in img_files) {
      l_cand <- .yolo_find_matching_label(im_f, dir_path)
      if (file.exists(l_cand)) {
        lines_cand <- readLines(l_cand, warn = FALSE)
        lines_cand <- lines_cand[nzchar(trimws(lines_cand))]
        if (length(lines_cand) > 0L) {
          vals <- as.numeric(strsplit(trimws(lines_cand[1]), "\\s+")[[1]])
          if (length(vals) > 5L && (length(vals) - 5L) %% 3L == 0L) {
            kpts <- as.integer((length(vals) - 5L) / 3L)
            break
          }
        }
      }
    }
    if (is.null(kpts) || is.na(kpts) || kpts < 1L) kpts <- 2L
  }

  skeleton_mat <- .yolo_resolve_skeleton(skeleton, kpts)

  # Palette
  kpt_colors <- if (isTRUE(rainbow)) {
    grDevices::rainbow(kpts, s = 0.90, v = 0.95)
  } else {
    rep("#00FF66", kpts)
  }

  # Layout
  if (is.null(mfrow)) {
    if (total_to_plot == 1L) mfrow <- c(1L, 1L)
    else if (total_to_plot == 2L) mfrow <- c(1L, 2L)
    else if (total_to_plot <= 4L) mfrow <- c(2L, 2L)
    else if (total_to_plot <= 6L) mfrow <- c(2L, 3L)
    else if (total_to_plot <= 9L) mfrow <- c(3L, 3L)
    else {
      cols <- ceiling(sqrt(total_to_plot))
      rows <- ceiling(total_to_plot / cols)
      mfrow <- c(rows, cols)
    }
  }

  # Save setup if requested
  if (isTRUE(save)) {
    if (is.null(file_name) || !nzchar(file_name)) {
      file_name <- paste0("yolo_pose_preview_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    }
    w_px <- if (!is.null(width)) width else (mfrow[2] * 600)
    h_px <- if (!is.null(height)) height else (mfrow[1] * 500)
    grDevices::png(file_name, width = w_px, height = h_px, res = dpi)
    on.exit({
      if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    }, add = TRUE)
  }

  old_par <- graphics::par(c("mfrow", "mar"))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = mfrow, mar = c(1.5, 1.5, 2.5, 1.5))

  annotated_results <- vector("list", total_to_plot)

  for (idx in seq_along(img_files)) {
    im_path <- img_files[idx]
    stem <- tools::file_path_sans_ext(basename(im_path))
    lbl_path <- .yolo_find_matching_label(im_path, dir_path)

    im <- tryCatch(image_import(im_path), error = function(e) NULL)
    if (is.null(im)) next

    dims <- dim(im)
    img_w <- dims[1]; img_h <- dims[2]

    parsed <- .yolo_read_pose_label(lbl_path, img_w = img_w, img_h = img_h, kpts = kpts, class_names = class_names)
    boxes <- parsed$boxes
    keypoints <- parsed$keypoints
    annotated_results[[idx]] <- list(image = im_path, boxes = boxes, keypoints = keypoints)

    # Plot base image
    plot(im)

    n_objs <- nrow(boxes)

    # 1. Draw Skeleton Lines
    if (isTRUE(show_skeleton) && n_objs > 0L && nrow(skeleton_mat) > 0L) {
      for (i in seq_len(n_objs)) {
        k_mat <- keypoints[[i]]
        if (!is.null(k_mat) && nrow(k_mat) >= 2L) {
          for (s in seq_len(nrow(skeleton_mat))) {
            p1 <- skeleton_mat[s, 1]
            p2 <- skeleton_mat[s, 2]
            if (p1 <= nrow(k_mat) && p2 <= nrow(k_mat)) {
              v1 <- k_mat[p1, 3]
              v2 <- k_mat[p2, 3]
              if (v1 > 0L && v2 > 0L) {
                sk_col <- if (v1 == 1L || v2 == 1L) "#FFAA00" else "#00FFCC"
                sk_lty <- if (v1 == 1L || v2 == 1L) 2 else 1
                graphics::lines(
                  c(k_mat[p1, 1], k_mat[p2, 1]),
                  c(k_mat[p1, 2], k_mat[p2, 2]),
                  col = sk_col,
                  lwd = line_lwd,
                  lty = sk_lty
                )
              }
            }
          }
        }
      }
    }

    # 2. Draw Bounding Boxes
    if (isTRUE(show_bbox) && n_objs > 0L) {
      box_cols <- if (!is.null(col)) rep(col, n_objs) else if (isTRUE(rainbow)) grDevices::rainbow(n_objs, s = 0.85, v = 0.90) else rep("#00FF66", n_objs)
      for (i in seq_len(n_objs)) {
        bx <- boxes[i, ]
        bc <- box_cols[((i - 1L) %% length(box_cols)) + 1L]
        graphics::rect(bx$xmin, bx$ymin, bx$xmax, bx$ymax, border = bc, lwd = line_lwd, lty = 2)

        if (isTRUE(show_labels)) {
          tag <- if (!is.null(bx$class_name) && nzchar(bx$class_name)) paste0("#", i, " ", bx$class_name) else paste0("#", i)
          tw <- graphics::strwidth(tag, cex = 0.75)
          th <- graphics::strheight(tag, cex = 0.75)
          y_top <- max(th + 6, bx$ymin)
          graphics::rect(bx$xmin, y_top - th - 5, bx$xmin + tw + 6, y_top,
                         col = grDevices::adjustcolor(bc, alpha.f = 0.85), border = NA)
          graphics::text(bx$xmin + 3, y_top - 2.5, labels = tag, col = "#FFFFFF",
                         cex = 0.75, font = 2, adj = c(0, 1))
        }
      }
    }

    # 3. Draw Keypoint Dots
    if (isTRUE(show_points) && n_objs > 0L) {
      for (i in seq_len(n_objs)) {
        k_mat <- keypoints[[i]]
        if (!is.null(k_mat) && nrow(k_mat) > 0L) {
          for (k in seq_len(nrow(k_mat))) {
            v <- k_mat[k, 3]
            if (v == 0L) next
            kc <- kpt_colors[((k - 1L) %% length(kpt_colors)) + 1L]
            kx <- k_mat[k, 1]
            ky <- k_mat[k, 2]

            if (v == 2L) {
              graphics::points(kx, ky, pch = 21, bg = kc, col = "#000000", cex = pt_cex, lwd = 1.2)
            } else if (v == 1L) {
              graphics::points(kx, ky, pch = 21, bg = "transparent", col = kc, cex = pt_cex * 1.1, lwd = 2.0)
              graphics::points(kx, ky, pch = 3, col = kc, cex = pt_cex * 0.7, lwd = 1.5)
            }

            if (isTRUE(show_kpt_labels)) {
              kname <- if (!is.null(kpt_names) && length(kpt_names) >= k) kpt_names[k] else paste0("P", k)
              graphics::text(kx + 6, ky - 6, labels = kname, col = kc, cex = 0.70, font = 2, adj = c(0, 0))
            }
          }
        }
      }
    }

    # Title
    if (isTRUE(title)) {
      graphics::title(sprintf("%s (%d object%s)", stem, n_objs, if (n_objs == 1) "" else "s"),
                      cex.main = 1.05, font.main = 2, col.main = "#222222")
    }
  }

  invisible(annotated_results)
}

#' @rdname yolo_dataset_pose_preview
#' @export
yolo_pose_preview <- yolo_dataset_pose_preview

.yolo_split_pose_dataset_inplace <- function(dir_path, split_ratio = c(train = 0.70, val = 0.20, test = 0.10), verbose = TRUE) {
  img_ext_patt <- "(?i)\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  fls <- list.files(dir_path, pattern = img_ext_patt, full.names = TRUE, recursive = FALSE)
  if (length(fls) == 0L) return(invisible(NULL))

  split_norm <- split_ratio / sum(split_ratio)
  split_names <- if (length(split_ratio) == 2L) c("train", "val") else c("train", "val", "test")
  names(split_norm) <- split_names

  for (s in split_names) {
    dir.create(file.path(dir_path, "images", s), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir_path, "labels", s), recursive = TRUE, showWarnings = FALSE)
  }

  set.seed(42)
  shuffled <- sample(fls)
  n_tot <- length(shuffled)

  if (length(split_names) == 2L) {
    n_tr <- max(1L, round(n_tot * split_norm[1]))
    n_tr <- min(n_tr, n_tot - 1L)
    splits_list <- list(train = shuffled[seq_len(n_tr)], val = shuffled[(n_tr + 1L):n_tot])
  } else {
    n_tr <- max(1L, round(n_tot * split_norm[1]))
    n_val <- max(1L, round(n_tot * split_norm[2]))
    if (n_tr + n_val >= n_tot) { n_tr <- max(1L, n_tot - 2L); n_val <- 1L }
    splits_list <- list(
      train = shuffled[seq_len(n_tr)],
      val   = shuffled[(n_tr + 1L):(n_tr + n_val)],
      test  = shuffled[(n_tr + n_val + 1L):n_tot]
    )
  }

  for (s in split_names) {
    cur_imgs <- splits_list[[s]]
    for (im_f in cur_imgs) {
      stem <- tools::file_path_sans_ext(basename(im_f))
      dest_im <- file.path(dir_path, "images", s, basename(im_f))
      dest_lbl <- file.path(dir_path, "labels", s, paste0(stem, ".txt"))
      orig_lbl <- file.path(dir_path, paste0(stem, ".txt"))

      # Move image
      file.rename(im_f, dest_im)

      # Move label if exists
      if (file.exists(orig_lbl)) {
        file.rename(orig_lbl, dest_lbl)
      } else {
        file.create(dest_lbl)
      }
    }
  }

  # Update data.yaml
  yaml_file <- file.path(dir_path, "data.yaml")
  if (file.exists(yaml_file)) {
    lines <- readLines(yaml_file, warn = FALSE)
    tr_idx <- grep("^\\s*train\\s*:", lines)
    if (length(tr_idx) > 0L) lines[tr_idx[1]] <- "train: images/train"
    vl_idx <- grep("^\\s*val\\s*:", lines)
    if (length(vl_idx) > 0L) lines[vl_idx[1]] <- "val: images/val"
    ts_idx <- grep("^\\s*test\\s*:", lines)
    if (length(ts_idx) > 0L && "test" %in% split_names) {
      lines[ts_idx[1]] <- "test: images/test"
    } else if ("test" %in% split_names) {
      lines <- c(lines[1:vl_idx[1]], "test: images/test", lines[(vl_idx[1] + 1L):length(lines)])
    }
    writeLines(lines, yaml_file)
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Partitioned {n_tot} image/label pairs into {.path images/} and {.path labels/} ({paste(split_names, collapse = ', ')}).")
  }
  invisible(TRUE)
}

.yolo_read_pose_yaml_meta <- function(dir_path) {
  yaml_cand <- c(file.path(dir_path, "data.yaml"), file.path(dir_path, "dataset.yaml"), list.files(dir_path, pattern = "\\.ya?ml$", full.names = TRUE))
  yaml_cand <- yaml_cand[file.exists(yaml_cand)]
  res <- list(kpts = NULL, skeleton = NULL, classes = NULL)
  if (length(yaml_cand) == 0L) return(res)

  lines <- readLines(yaml_cand[1], warn = FALSE)

  k_idx <- grep("^\\s*kpt_shape\\s*:", lines)
  if (length(k_idx) > 0L) {
    m <- regmatches(lines[k_idx[1]], regexpr("\\[[0-9 ,]+\\]", lines[k_idx[1]]))
    if (length(m) > 0L) {
      nums <- as.integer(strsplit(gsub("[^0-9,]", "", m), ",")[[1]])
      if (length(nums) >= 1L && !is.na(nums[1])) res$kpts <- nums[1]
    }
  }

  s_idx <- grep("^\\s*skeleton\\s*:", lines)
  if (length(s_idx) > 0L) {
    sk_pairs <- list()
    end_s <- s_idx[1] + 1L
    while (end_s <= length(lines) && grepl("^[[:space:]]*-[[:space:]]*\\[", lines[end_s])) {
      m <- regmatches(lines[end_s], regexpr("\\[[0-9 ,]+\\]", lines[end_s]))
      if (length(m) > 0L) {
        nums <- as.integer(strsplit(gsub("[^0-9,]", "", m), ",")[[1]])
        if (length(nums) >= 2L && !any(is.na(nums[1:2]))) sk_pairs[[length(sk_pairs) + 1L]] <- nums[1:2]
      }
      end_s <- end_s + 1L
    }
    if (length(sk_pairs) > 0L) {
      res$skeleton <- do.call(rbind, sk_pairs)
    }
  }

  res$classes <- tryCatch(.yolo_get_classes(dir_path), error = function(e) NULL)
  res
}



# ==============================================================================
# INTERNAL HELPER FUNCTIONS
# ==============================================================================

.yolo_resolve_skeleton <- function(skeleton, kpts) {
  if (is.null(skeleton) || identical(skeleton, "chain") || isTRUE(skeleton)) {
    if (kpts >= 2L) {
      from <- seq_len(kpts - 1L)
      to <- from + 1L
      return(cbind(from, to))
    } else {
      return(matrix(0L, nrow = 0L, ncol = 2L))
    }
  }
  if (isFALSE(skeleton) || identical(tolower(skeleton[1]), "none")) {
    return(matrix(0L, nrow = 0L, ncol = 2L))
  }
  if (is.matrix(skeleton) && ncol(skeleton) == 2L) {
    return(skeleton)
  }
  if (is.list(skeleton)) {
    pairs <- do.call(rbind, lapply(skeleton, function(p) as.integer(p[1:2])))
    return(pairs)
  }
  matrix(0L, nrow = 0L, ncol = 2L)
}

.yolo_find_matching_label <- function(img_path, dir_path) {
  stem <- tools::file_path_sans_ext(basename(img_path))
  img_d <- dirname(img_path)

  # Check standard YOLO images/ -> labels/ structure
  if (grepl("(/|\\\\)images(/|\\\\|$)", img_d)) {
    lbl_d <- sub("((/|\\\\))images((/|\\\\)|$)", "\\1labels\\3", img_d)
    return(file.path(lbl_d, paste0(stem, ".txt")))
  }

  # Or labels directory at root
  if (dir.exists(file.path(dir_path, "labels"))) {
    return(file.path(dir_path, "labels", paste0(stem, ".txt")))
  }

  # Same folder as image
  file.path(img_d, paste0(stem, ".txt"))
}

.yolo_read_pose_label <- function(lbl_path, img_w, img_h, kpts, class_names = "object") {
  empty_res <- list(
    boxes = data.frame(
      index = integer(0), xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
      cx = numeric(0), cy = numeric(0), class_id = integer(0), class_name = character(0), stringsAsFactors = FALSE
    ),
    keypoints = list()
  )

  if (!file.exists(lbl_path)) return(empty_res)

  lines <- readLines(lbl_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0L) return(empty_res)

  boxes_list <- vector("list", length(lines))
  kpts_list <- vector("list", length(lines))

  for (i in seq_along(lines)) {
    vals <- as.numeric(strsplit(trimws(lines[i]), "\\s+")[[1]])
    if (length(vals) < 5L) next

    cid <- as.integer(vals[1])
    cx <- vals[2] * img_w
    cy <- vals[3] * img_h
    bw <- vals[4] * img_w
    bh <- vals[5] * img_h

    xmin <- cx - bw / 2
    xmax <- cx + bw / 2
    ymin <- cy - bh / 2
    ymax <- cy + bh / 2

    cname <- if (cid >= 0 && cid < length(class_names)) class_names[cid + 1L] else "object"

    boxes_list[[i]] <- data.frame(
      index = i, xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
      cx = cx, cy = cy, class_id = cid, class_name = cname, stringsAsFactors = FALSE
    )

    # Keypoints parsing (format: x1 y1 v1 x2 y2 v2 ...)
    kpt_mat <- matrix(0, nrow = kpts, ncol = 3L)
    colnames(kpt_mat) <- c("x", "y", "v")

    rem_vals <- vals[-(1:5)]
    if (length(rem_vals) >= (3L * kpts)) {
      # xyv format (Ultralytics standard)
      for (k in seq_len(kpts)) {
        kpt_mat[k, 1] <- rem_vals[(k - 1L) * 3L + 1L] * img_w
        kpt_mat[k, 2] <- rem_vals[(k - 1L) * 3L + 2L] * img_h
        kpt_mat[k, 3] <- as.integer(rem_vals[(k - 1L) * 3L + 3L])
      }
    } else if (length(rem_vals) >= (2L * kpts)) {
      # xy format (without visibility)
      for (k in seq_len(kpts)) {
        kpt_mat[k, 1] <- rem_vals[(k - 1L) * 2L + 1L] * img_w
        kpt_mat[k, 2] <- rem_vals[(k - 1L) * 2L + 2L] * img_h
        kpt_mat[k, 3] <- 2L
      }
    }
    kpts_list[[i]] <- kpt_mat
  }

  boxes_df <- do.call(rbind, boxes_list)
  if (!is.null(boxes_df) && nrow(boxes_df) > 0L) {
    boxes_df$index <- seq_len(nrow(boxes_df))
  } else {
    boxes_df <- empty_res$boxes
  }

  list(boxes = boxes_df, keypoints = kpts_list)
}

.yolo_write_pose_label <- function(lbl_path, boxes, keypoints, img_w, img_h, kpts) {
  if (nrow(boxes) == 0L) {
    # Write empty file (negative/background sample)
    writeLines(character(0), lbl_path)
    return(invisible(TRUE))
  }

  out_lines <- character(nrow(boxes))
  for (i in seq_len(nrow(boxes))) {
    cid <- boxes$class_id[i]
    cx <- round(boxes$cx[i] / img_w, 6)
    cy <- round(boxes$cy[i] / img_h, 6)
    bw <- round((boxes$xmax[i] - boxes$xmin[i]) / img_w, 6)
    bh <- round((boxes$ymax[i] - boxes$ymin[i]) / img_h, 6)

    # Clamping
    cx <- max(0.0, min(1.0, cx))
    cy <- max(0.0, min(1.0, cy))
    bw <- max(0.0001, min(1.0, bw))
    bh <- max(0.0001, min(1.0, bh))

    k_mat <- keypoints[[i]]
    k_parts <- character(kpts * 3L)
    for (k in seq_len(kpts)) {
      v <- if (!is.null(k_mat) && nrow(k_mat) >= k) as.integer(k_mat[k, 3]) else 0L
      if (v > 0L && !is.null(k_mat) && nrow(k_mat) >= k) {
        kx <- round(k_mat[k, 1] / img_w, 6)
        ky <- round(k_mat[k, 2] / img_h, 6)
        kx <- max(0.0, min(1.0, kx))
        ky <- max(0.0, min(1.0, ky))
      } else {
        kx <- 0.0
        ky <- 0.0
        v <- 0L
      }
      k_parts[(k - 1L) * 3L + 1L] <- sprintf("%.6f", kx)
      k_parts[(k - 1L) * 3L + 2L] <- sprintf("%.6f", ky)
      k_parts[(k - 1L) * 3L + 3L] <- as.character(v)
    }

    out_lines[i] <- paste(c(cid, sprintf("%.6f", cx), sprintf("%.6f", cy), sprintf("%.6f", bw), sprintf("%.6f", bh), k_parts), collapse = " ")
  }

  writeLines(out_lines, lbl_path)
  invisible(TRUE)
}

.yolo_render_pose_canvas <- function(img,
                                     boxes,
                                     keypoints,
                                     kpt_names,
                                     skeleton,
                                     img_name,
                                     curr_idx,
                                     total_imgs,
                                     box_mode,
                                     kpt_colors,
                                     pt_cex = 1.8,
                                     line_lwd = 2.0) {
  dims <- dim(img)
  img_w <- dims[1]
  img_h <- dims[2]
  n_objs <- nrow(boxes)

  title_str <- sprintf("YOLO-Pose [%s | Image %d of %d | %d object%s | Box: %s]",
                       img_name, curr_idx, total_imgs, n_objs, if (n_objs == 1) "" else "s", toupper(box_mode))
  sub_str <- "[a] add pose | [aa] rapid mode | [m] toggle box mode | [r] remove | [Enter] next | [q] exit"

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

  if (n_objs == 0L) return(invisible(NULL))

  obj_colors <- grDevices::rainbow(n_objs, s = 0.85, v = 0.95)

  for (i in seq_len(n_objs)) {
    bx <- boxes[i, ]
    oc <- obj_colors[((i - 1L) %% length(obj_colors)) + 1L]

    # Draw bounding box
    graphics::rect(bx$xmin, bx$ymin, bx$xmax, bx$ymax, border = oc, lwd = line_lwd, lty = 2)

    # Object badge
    tag <- sprintf("#%d %s", i, bx$class_name)
    tw <- graphics::strwidth(tag, cex = 0.72)
    th <- graphics::strheight(tag, cex = 0.72)
    y_badge <- max(th + 6, bx$ymin)
    graphics::rect(bx$xmin, y_badge - th - 4, bx$xmin + tw + 6, y_badge, col = grDevices::adjustcolor(oc, alpha.f = 0.85), border = NA)
    graphics::text(bx$xmin + 3, y_badge - 2, labels = tag, col = "#FFFFFF", cex = 0.70, font = 2, adj = c(0, 1))

    # Draw skeleton lines
    k_mat <- keypoints[[i]]
    if (!is.null(k_mat) && nrow(skeleton) > 0L) {
      for (sk in seq_len(nrow(skeleton))) {
        p1 <- skeleton[sk, 1]
        p2 <- skeleton[sk, 2]
        if (p1 <= nrow(k_mat) && p2 <= nrow(k_mat)) {
          v1 <- k_mat[p1, 3]
          v2 <- k_mat[p2, 3]
          if (v1 > 0L && v2 > 0L) {
            line_col <- if (v1 == 2L && v2 == 2L) "#00FFEE" else "#FFAA00"
            graphics::lines(c(k_mat[p1, 1], k_mat[p2, 1]), c(k_mat[p1, 2], k_mat[p2, 2]), col = line_col, lwd = line_lwd + 0.5)
          }
        }
      }
    }

    # Draw keypoints
    if (!is.null(k_mat)) {
      for (k in seq_len(nrow(k_mat))) {
        kx <- k_mat[k, 1]
        ky <- k_mat[k, 2]
        kv <- k_mat[k, 3]
        kc <- kpt_colors[((k - 1L) %% length(kpt_colors)) + 1L]

        if (kv == 2L) {
          # Visible: solid circle with dark border
          graphics::points(kx, ky, pch = 21, bg = kc, col = "#000000", cex = pt_cex, lwd = 1.5)
        } else if (kv == 1L) {
          # Occluded: hollow square with cross
          graphics::points(kx, ky, pch = 22, bg = "transparent", col = kc, cex = pt_cex * 1.2, lwd = 2)
          graphics::points(kx, ky, pch = 4, col = kc, cex = pt_cex * 0.9, lwd = 1.8)
        }

        if (kv > 0L) {
          k_label <- sprintf("%d", k)
          graphics::text(kx + 7, ky - 7, labels = k_label, col = "#000000", cex = 0.70, font = 2)
          graphics::text(kx + 6, ky - 8, labels = k_label, col = kc, cex = 0.70, font = 2)
        }
      }
    }
  }

  invisible(NULL)
}

.yolo_annotate_single_pose <- function(img_w,
                                       img_h,
                                       kpts,
                                       kpt_names,
                                       skeleton,
                                       box_mode,
                                       box_pad,
                                       allow_occlusion,
                                       class_id,
                                       class_name,
                                       kpt_colors,
                                       pt_cex = 1.8,
                                       line_lwd = 2.0) {
  manual_box <- (box_mode == "manual")

  # Strategy 1: Manual Bounding Box (Click 2 corners first)
  bx_min <- 0; bx_max <- 0; by_min <- 0; by_max <- 0
  if (manual_box) {
    cli::cli_alert_info("Click TWO opposite corners of the bounding box (or <Esc>/right-click to cancel):")
    b_pts <- tryCatch(
      graphics::locator(n = 2, type = "p", pch = 3, col = "#00FF66", cex = 1.6, lwd = 2.2),
      error = function(e) NULL
    )
    if (is.null(b_pts) || length(b_pts$x) < 2L) {
      cli::cli_alert_warning("Bounding box addition cancelled.")
      return(NULL)
    }
    bx_min <- max(0, min(b_pts$x))
    bx_max <- min(img_w, max(b_pts$x))
    by_min <- max(0, min(b_pts$y))
    by_max <- min(img_h, max(b_pts$y))

    # Draw temporary manual box
    graphics::rect(bx_min, by_min, bx_max, by_max, border = "#00FF66", lty = 2, lwd = line_lwd)
  }

  # Strategy 2 & 3: Keypoints Clicking with Occlusion Control
  pts_x <- numeric(kpts)
  pts_y <- numeric(kpts)
  pts_v <- integer(kpts)

  for (k in seq_len(kpts)) {
    kname <- if (!is.null(kpt_names) && length(kpt_names) >= k) kpt_names[k] else paste0("point_", k)
    kc <- kpt_colors[((k - 1L) %% length(kpt_colors)) + 1L]

    cli::cli_alert_info("Click point {k}/{kpts} [{kname}] (Left-click = Visible | Right-click/Esc = Occluded/Skip):")

    pt <- tryCatch(
      graphics::locator(n = 1, type = "n"),
      error = function(e) NULL
    )

    if (is.null(pt) || length(pt$x) == 0L) {
      # Strategy 3: Occlusion / Skip handling on Right-click or Esc
      if (isTRUE(allow_occlusion)) {
        cli::cli_alert_warning("Point {k}/{kpts} [{kname}] interrupted.")
        occ_ans <- tolower(trimws(readline("Select flag [o = Occluded (1) | s = Skip/Absent (0) | c = Cancel object]: ")))
        if (occ_ans %in% c("o", "1", "occ", "occluded")) {
          cli::cli_alert_info("Click estimated occluded location for [{kname}]:")
          pt_occ <- tryCatch(graphics::locator(n = 1, type = "n"), error = function(e) NULL)
          if (!is.null(pt_occ) && length(pt_occ$x) > 0L) {
            pts_x[k] <- pt_occ$x[1]
            pts_y[k] <- pt_occ$y[1]
            pts_v[k] <- 1L
          } else {
            pts_x[k] <- 0.0; pts_y[k] <- 0.0; pts_v[k] <- 0L
          }
        } else if (occ_ans %in% c("s", "0", "skip", "absent")) {
          pts_x[k] <- 0.0; pts_y[k] <- 0.0; pts_v[k] <- 0L
        } else {
          cli::cli_alert_warning("Object annotation cancelled.")
          return(NULL)
        }
      } else {
        cli::cli_alert_warning("Object annotation cancelled.")
        return(NULL)
      }
    } else {
      # Normal visible point
      pts_x[k] <- pt$x[1]
      pts_y[k] <- pt$y[1]
      pts_v[k] <- 2L
    }

    # Draw point immediately for feedback
    if (pts_v[k] == 2L) {
      graphics::points(pts_x[k], pts_y[k], pch = 21, bg = kc, col = "#000000", cex = pt_cex, lwd = 1.5)
      graphics::text(pts_x[k] + 7, pts_y[k] - 7, labels = as.character(k), col = kc, cex = 0.75, font = 2)
    } else if (pts_v[k] == 1L) {
      graphics::points(pts_x[k], pts_y[k], pch = 22, bg = "transparent", col = kc, cex = pt_cex * 1.2, lwd = 2)
      graphics::points(pts_x[k], pts_y[k], pch = 4, col = kc, cex = pt_cex * 0.9, lwd = 1.8)
    }

    # Draw skeleton lines dynamically as points are clicked
    if (nrow(skeleton) > 0L && pts_v[k] > 0L) {
      for (sk in seq_len(nrow(skeleton))) {
        p1 <- skeleton[sk, 1]
        p2 <- skeleton[sk, 2]
        if (p2 == k && p1 < k && pts_v[p1] > 0L) {
          graphics::lines(c(pts_x[p1], pts_x[k]), c(pts_y[p1], pts_y[k]), col = "#00FFEE", lwd = line_lwd)
        }
      }
    }
  }

  # Strategy 2: Auto Bounding Box calculation
  if (!manual_box) {
    valid_idx <- which(pts_v > 0L)
    if (length(valid_idx) == 0L) {
      bx_min <- 0; bx_max <- 10; by_min <- 0; by_max <- 10
    } else {
      min_x <- min(pts_x[valid_idx]); max_x <- max(pts_x[valid_idx])
      min_y <- min(pts_y[valid_idx]); max_y <- max(pts_y[valid_idx])
      span_x <- max_x - min_x
      span_y <- max_y - min_y

      pad_val <- if (is.numeric(box_pad) && length(box_pad) >= 1L) max(0.01, box_pad[1]) else 0.10
      pad_px_x <- max(span_x * pad_val, img_w * 0.025, 12)
      pad_px_y <- max(span_y * pad_val, img_h * 0.025, 12)

      bx_min <- max(0, min_x - pad_px_x)
      bx_max <- min(img_w, max_x + pad_px_x)
      by_min <- max(0, min_y - pad_px_y)
      by_max <- min(img_h, max_y + pad_px_y)
    }
  }

  box_df <- data.frame(
    index = 1L,
    xmin = bx_min, ymin = by_min, xmax = bx_max, ymax = by_max,
    cx = (bx_min + bx_max) / 2, cy = (by_min + by_max) / 2,
    class_id = class_id, class_name = class_name, stringsAsFactors = FALSE
  )

  kpts_mat <- cbind(x = pts_x, y = pts_y, v = pts_v)

  list(box = box_df, kpts = kpts_mat)
}

.yolo_sync_pose_yaml <- function(dir_path, kpts, skeleton, class_names) {
  yaml_file <- file.path(dir_path, "data.yaml")
  if (!file.exists(yaml_file)) {
    yaml_candidates <- c(file.path(dir_path, "dataset.yaml"), list.files(dir_path, pattern = "\\.ya?ml$", full.names = TRUE))
    if (length(yaml_candidates) > 0L && file.exists(yaml_candidates[1])) {
      yaml_file <- yaml_candidates[1]
    }
  }

  # Build skeleton lines for YAML
  sk_lines <- character(0)
  if (nrow(skeleton) > 0L) {
    sk_lines <- c("skeleton:")
    for (s in seq_len(nrow(skeleton))) {
      sk_lines <- c(sk_lines, sprintf("  - [%d, %d]", skeleton[s, 1], skeleton[s, 2]))
    }
  }

  if (file.exists(yaml_file)) {
    lines <- readLines(yaml_file, warn = FALSE)

    # Replace or add kpt_shape
    kpt_idx <- grep("^\\s*kpt_shape\\s*:", lines)
    kpt_line <- sprintf("kpt_shape: [%d, 3]", kpts)
    if (length(kpt_idx) > 0L) {
      lines[kpt_idx[1]] <- kpt_line
    } else {
      lines <- c(lines, kpt_line)
    }

    # Replace or add skeleton
    sk_idx <- grep("^\\s*skeleton\\s*:", lines)
    if (length(sk_idx) > 0L) {
      # Find end of existing skeleton block
      end_sk <- sk_idx[1]
      while ((end_sk + 1L) <= length(lines) && grepl("^[[:space:]]*-[[:space:]]*\\[", lines[end_sk + 1L])) {
        end_sk <- end_sk + 1L
      }
      lines <- c(lines[seq_len(sk_idx[1] - 1L)], sk_lines, if (end_sk < length(lines)) lines[(end_sk + 1L):length(lines)] else character(0))
    } else if (length(sk_lines) > 0L) {
      lines <- c(lines, sk_lines)
    }

    writeLines(lines, yaml_file)
  } else {
    # Create new data.yaml
    tr_rel <- if (dir.exists(file.path(dir_path, "images", "train"))) "images/train" else "train"
    vl_rel <- if (dir.exists(file.path(dir_path, "images", "val"))) "images/val" else "val"
    new_lines <- c(
      paste0("path: ", dir_path),
      paste0("train: ", tr_rel),
      paste0("val: ", vl_rel),
      "",
      "names:",
      sprintf("  0: %s", class_names[1]),
      "",
      sprintf("kpt_shape: [%d, 3]", kpts),
      sk_lines
    )
    writeLines(new_lines, yaml_file)
  }
}
