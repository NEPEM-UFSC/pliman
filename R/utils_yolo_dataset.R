# utils_yolo_dataset.R - Comprehensive, Pipe-Friendly YOLO Dataset Architecture for pliman
# Author: Tiago Olivoto <tiagoolivoto@gmail.com>

# ==============================================================================
# 1. INTERNAL HELPERS
# ==============================================================================
.yolo_get_dir <- function(data) {
  if (is.null(data)) {
    cli::cli_abort("The {.arg data} argument cannot be NULL.")
  }
  if (is.list(data)) {
    if (!is.null(data$dir) && is.character(data$dir) && nzchar(data$dir[1])) {
      return(normalizePath(data$dir[1], winslash = "/", mustWork = FALSE))
    }
    if (!is.null(data$yaml_file) && is.character(data$yaml_file) && nzchar(data$yaml_file[1])) {
      return(normalizePath(dirname(data$yaml_file[1]), winslash = "/", mustWork = FALSE))
    }
    if (!is.null(data$data_yaml) && is.character(data$data_yaml) && nzchar(data$data_yaml[1])) {
      return(normalizePath(dirname(data$data_yaml[1]), winslash = "/", mustWork = FALSE))
    }
    if (!is.null(data$dataset) && is.list(data$dataset)) {
      return(.yolo_get_dir(data$dataset))
    }
    if (length(data) > 0L && is.character(data[[1]]) && nzchar(data[[1]])) {
      return(normalizePath(data[[1]], winslash = "/", mustWork = FALSE))
    }
  }
  raw_path <- as.character(data)[1]
  if (!file.exists(raw_path) && !dir.exists(raw_path)) {
    for (cand in c(file.path(getwd(), raw_path), file.path(tempdir(), raw_path))) {
      if (file.exists(cand) || dir.exists(cand)) {
        raw_path <- cand
        break
      }
    }
  }
  if (dir.exists(raw_path)) {
    return(normalizePath(raw_path, winslash = "/", mustWork = FALSE))
  }
  if (file.exists(raw_path)) {
    return(normalizePath(dirname(raw_path), winslash = "/", mustWork = FALSE))
  }
  normalizePath(raw_path, winslash = "/", mustWork = FALSE)
}

.yolo_detect_task <- function(dir) {
  # 1. First check YAML file for explicit task / kpt_shape definition
  yaml_files <- c(file.path(dir, "data.yaml"), file.path(dir, "dataset.yaml"))
  for (yf in yaml_files) {
    if (file.exists(yf)) {
      y_lines <- readLines(yf, warn = FALSE)
      if (any(grepl("(kpt_shape|task:\\s*pose)", y_lines, ignore.case = TRUE))) return("pose")
      if (any(grepl("task:\\s*obb", y_lines, ignore.case = TRUE))) return("obb")
      if (any(grepl("task:\\s*segment", y_lines, ignore.case = TRUE))) return("segment")
      if (any(grepl("task:\\s*classify", y_lines, ignore.case = TRUE))) return("classify")
      if (any(grepl("task:\\s*detect", y_lines, ignore.case = TRUE))) return("detect")
    }
  }

  tr_dir <- file.path(dir, "train")
  val_dir <- file.path(dir, "val")
  lbl_dir <- file.path(dir, "labels")

  if (!dir.exists(lbl_dir) && dir.exists(tr_dir) && dir.exists(val_dir)) {
    tr_subs <- list.dirs(tr_dir, recursive = FALSE, full.names = FALSE)
    tr_subs <- tr_subs[nzchar(tr_subs) & !startsWith(tr_subs, ".")]
    if (length(tr_subs) >= 2L) {
      return("classify")
    }
  }

  lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(lbl_files) > 0L) {
    has_9 <- FALSE
    has_seg <- FALSE
    has_box <- FALSE
    for (lf in lbl_files[seq_len(min(20L, length(lbl_files)))]) {
      if (!file.exists(lf)) next
      lines_lf <- tryCatch(readLines(lf, warn = FALSE), error = function(e) character(0))
      lines_lf <- lines_lf[nzchar(trimws(lines_lf))]
      for (ln in lines_lf) {
        vals <- strsplit(trimws(ln), "\\s+")[[1]]
        if (length(vals) == 5L) {
          has_box <- TRUE
        } else if (length(vals) == 9L) {
          has_9 <- TRUE
        } else if (length(vals) >= 6L) {
          has_seg <- TRUE
        }
      }
    }
    if (has_9 && !has_seg && !has_box) return("obb")
    if (has_seg || has_9) return("segment")
    if (has_box) return("detect")
    return("detect")
  }

  "detect"
}

.yolo_get_classes <- function(dir) {
  yaml_path <- file.path(dir, "data.yaml")
  if (file.exists(yaml_path)) {
    lines <- readLines(yaml_path, warn = FALSE)
    in_names <- FALSE
    cls_list <- character()
    for (ln in lines) {
      t_ln <- trimws(ln)
      if (grepl("^names:", t_ln)) {
        in_names <- TRUE
        m <- regmatches(t_ln, regexpr("\\[.*\\]", t_ln))
        if (length(m) > 0L) {
          items <- gsub("[\\[\\]'\" ]", "", strsplit(m, ",")[[1]])
          return(items[nzchar(items)])
        }
        next
      }
      if (in_names) {
        if (grepl("^[a-zA-Z_]+:", t_ln) && !grepl("^\\d+:", t_ln)) {
          break
        }
        val <- sub("^[-0-9:'\" ]+", "", t_ln)
        val <- gsub("['\"]", "", trimws(val))
        if (nzchar(val)) cls_list <- c(cls_list, val)
      }
    }
    if (length(cls_list) > 0L) return(cls_list)
  }

  tr_dir <- file.path(dir, "train")
  if (dir.exists(tr_dir)) {
    subs <- list.dirs(tr_dir, recursive = FALSE, full.names = FALSE)
    subs <- subs[nzchar(subs) & !startsWith(subs, ".")]
    if (length(subs) > 0L) return(subs)
  }

  character()
}

.yolo_write_yaml <- function(out_dir, classes, train_dir = "images/train", val_dir = "images/val", test_dir = NULL, task = "detect") {
  yaml_path <- file.path(out_dir, "data.yaml")
  lines <- c(
    paste0("# YOLO dataset configuration generated by pliman"),
    paste0("path: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
    paste0("train: ", train_dir),
    paste0("val: ", val_dir)
  )
  if (!is.null(test_dir) && nzchar(test_dir) && dir.exists(file.path(out_dir, test_dir))) {
    lines <- c(lines, paste0("test: ", test_dir))
  }
  if (!is.null(task) && nzchar(task)) {
    lines <- c(lines, paste0("task: ", task))
  }
  lines <- c(lines, paste0("names:"))
  for (i in seq_along(classes)) {
    lines <- c(lines, paste0("  ", i - 1L, ": ", classes[i]))
  }
  writeLines(lines, yaml_path)
  yaml_path
}

.yolo_sync_data_yaml <- function(data_dir, yaml_file = NULL, verbose = TRUE) {
  if (is.null(data_dir) || !dir.exists(data_dir)) return(yaml_file)
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
  if (is.null(yaml_file) || !file.exists(yaml_file)) {
    yaml_candidates <- c(
      file.path(data_dir, "data.yaml"),
      file.path(data_dir, "dataset.yaml")
    )
    existing <- yaml_candidates[file.exists(yaml_candidates)]
    yaml_file <- if (length(existing) > 0L) existing[1] else file.path(data_dir, "data.yaml")
  }
  yaml_file <- normalizePath(yaml_file, winslash = "/", mustWork = FALSE)

  # Check actual folder structure on disk
  tr_exists <- dir.exists(file.path(data_dir, "images", "train")) || dir.exists(file.path(data_dir, "train", "images"))
  val_exists <- dir.exists(file.path(data_dir, "images", "val")) || dir.exists(file.path(data_dir, "val", "images"))
  test_exists <- dir.exists(file.path(data_dir, "images", "test")) || dir.exists(file.path(data_dir, "test", "images"))

  if (!file.exists(yaml_file)) {
    cl_names <- .yolo_get_classes(data_dir)
    if (is.null(cl_names) || length(cl_names) == 0L) cl_names <- c("object")
    tr_rel <- if (dir.exists(file.path(data_dir, "images", "train"))) "images/train" else "train"
    vl_rel <- if (dir.exists(file.path(data_dir, "images", "val"))) "images/val" else "val"
    ts_rel <- if (dir.exists(file.path(data_dir, "images", "test"))) "images/test" else if (test_exists) "test" else NULL
    .yolo_write_yaml(data_dir, cl_names, train_dir = tr_rel, val_dir = vl_rel, test_dir = ts_rel)
    return(yaml_file)
  }

  # File exists: inspect and check for stale / moved paths
  lines <- readLines(yaml_file, warn = FALSE)
  needs_repair <- FALSE

  path_line_idx <- grep("^\\s*path\\s*:", lines)
  train_line_idx <- grep("^\\s*train\\s*:", lines)
  val_line_idx <- grep("^\\s*val\\s*:", lines)

  # Check path: line
  curr_path <- ""
  if (length(path_line_idx) > 0L) {
    curr_path <- trimws(sub("^\\s*path\\s*:\\s*", "", lines[path_line_idx[1]]))
    curr_path <- gsub("['\"]", "", curr_path)
    if (curr_path == "." || !nzchar(curr_path)) {
      curr_path <- data_dir
      needs_repair <- TRUE
    } else {
      norm_curr <- tryCatch(normalizePath(curr_path, winslash = "/", mustWork = FALSE), error = function(e) curr_path)
      if (norm_curr != data_dir) {
        needs_repair <- TRUE
      }
    }
  } else {
    curr_path <- data_dir
    needs_repair <- TRUE
  }

  # Check if train or val line contains an outdated folder name prefix (e.g. "dataset_complex/images/val")
  if (length(train_line_idx) > 0L) {
    tr_val <- trimws(sub("^\\s*train\\s*:\\s*", "", lines[train_line_idx[1]]))
    tr_val <- gsub("['\"]", "", tr_val)
    full_tr <- if (file.exists(tr_val) || dir.exists(tr_val)) tr_val else file.path(curr_path, tr_val)
    if (!dir.exists(full_tr) && !file.exists(full_tr)) {
      needs_repair <- TRUE
    }
  } else {
    needs_repair <- TRUE
  }

  if (length(val_line_idx) > 0L) {
    vl_val <- trimws(sub("^\\s*val\\s*:\\s*", "", lines[val_line_idx[1]]))
    vl_val <- gsub("['\"]", "", vl_val)
    full_vl <- if (file.exists(vl_val) || dir.exists(vl_val)) vl_val else file.path(curr_path, vl_val)
    if (!dir.exists(full_vl) && !file.exists(full_vl)) {
      needs_repair <- TRUE
    }
  } else {
    needs_repair <- TRUE
  }

  if (isTRUE(needs_repair)) {
    cl_names <- .yolo_get_classes(data_dir)
    if (is.null(cl_names) || length(cl_names) == 0L) cl_names <- c("object")
    tr_rel <- if (dir.exists(file.path(data_dir, "images", "train"))) "images/train" else "train"
    vl_rel <- if (dir.exists(file.path(data_dir, "images", "val"))) "images/val" else "val"
    ts_rel <- if (dir.exists(file.path(data_dir, "images", "test"))) "images/test" else if (test_exists) "test" else NULL
    .yolo_write_yaml(data_dir, cl_names, train_dir = tr_rel, val_dir = vl_rel, test_dir = ts_rel)
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Healed {.file {basename(yaml_file)}}: synchronized dataset paths to {.path {data_dir}}.")
    }
  }

  return(yaml_file)
}

.check_in_place <- function(target_dir, check = TRUE) {
  if (!isTRUE(check)) return(TRUE)
  if (!interactive()) return(TRUE)

  cli::cli_alert_warning(
    "Caution: {.code in_place = TRUE} will permanently replace and overwrite existing images and labels in {.path {target_dir}}."
  )
  cli::cli_alert_danger(
    "This operation cannot be undone and original files will be lost if you do not have a backup!"
  )

  ans <- readline("Are you sure you want to proceed and overwrite existing files? [y/N]: ")
  confirmed <- tolower(trimws(ans)) %in% c("y", "yes", "s", "sim")
  if (!confirmed) {
    cli::cli_alert_info("Operation cancelled by user.")
    return(FALSE)
  }
  return(TRUE)
}

.yolo_check_in_place <- .check_in_place

.yolo_wrap_dataset <- function(dir, task = NULL, classes = NULL, n_objects = NULL) {
  dir_norm <- normalizePath(dir, winslash = "/", mustWork = FALSE)
  if (is.null(task)) task <- .yolo_detect_task(dir_norm)
  if (is.null(classes) || length(classes) == 0L) classes <- .yolo_get_classes(dir_norm)

  yaml_path <- file.path(dir_norm, "data.yaml")

  if (task == "classify") {
    splits <- c("train", "val", "test")
    splits <- splits[dir.exists(file.path(dir_norm, splits))]
    summary_rows <- list()
    for (cls in classes) {
      r <- list(class = cls)
      for (s in splits) {
        c_dir <- file.path(dir_norm, s, cls)
        r[[s]] <- if (dir.exists(c_dir)) length(list.files(c_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", ignore.case = TRUE)) else 0L
      }
      r$total <- sum(unlist(r[splits]))
      summary_rows[[cls]] <- as.data.frame(r, stringsAsFactors = FALSE)
    }
    summary_df <- if (length(summary_rows) > 0L) do.call(rbind, summary_rows) else data.frame()
    if (nrow(summary_df) > 0L) {
      tot_row <- list(class = "TOTAL")
      for (nm in setdiff(names(summary_df), "class")) tot_row[[nm]] <- sum(summary_df[[nm]])
      summary_full <- rbind(summary_df, as.data.frame(tot_row, stringsAsFactors = FALSE))
    } else {
      summary_full <- summary_df
    }

    out <- list(
      dir = dir_norm,
      classes = classes,
      splits = splits,
      task = "classify",
      summary = summary_full,
      total_images = if (nrow(summary_df) > 0L) sum(summary_df$total) else 0L
    )
    class(out) <- c("yolo_dataset_cls", "yolo_dataset", "list")
    return(out)
  }

  tr_imgs <- list.files(file.path(dir_norm, "images", "train"), pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
  val_imgs <- list.files(file.path(dir_norm, "images", "val"), pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
  test_imgs <- list.files(file.path(dir_norm, "images", "test"), pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
  if (!is.null(n_objects)) {
    n_objs <- as.integer(n_objects)
  } else {
    all_lbls <- list.files(file.path(dir_norm, "labels"), pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
    n_objs <- fast_count_yolo_objects_cpp(all_lbls)
  }

  out <- list(
    dir = dir_norm,
    yaml_file = if (file.exists(yaml_path)) yaml_path else NULL,
    data_yaml = if (file.exists(yaml_path)) yaml_path else NULL,
    classes = classes,
    task = task,
    train_images = tr_imgs,
    val_images = val_imgs,
    test_images = test_imgs,
    n_objects = n_objs
  )
  class(out) <- c("yolo_dataset", "list")
  out
}


# ==============================================================================
# 2. YOLO DATASET HEALTH AUDIT & REPAIR (yolo_dataset_audit)
# ==============================================================================
#' @title Comprehensive YOLO Dataset Health Audit and Auto-Repair
#' @name yolo_dataset_audit
#' @aliases yolo_dataset_health
#' @description
#' `yolo_dataset_audit()` conducts a deep structural, syntactic, and geometric health
#' check on any YOLO dataset (Detection, Instance Segmentation, or Image Classification).
#' It detects duplicate file stems across extensions (stem collisions, e.g. `foo.jpg` and `foo.png`
#' colliding on `foo.txt`), orphaned images/labels, corrupt files, invalid coordinates (< 0 or > 1),
#' negative bounding box dimensions, malformed polygons, class index out-of-bounds,
#' and severe class imbalance.
#'
#' With `repair = TRUE`, it automatically repairs the dataset in-place:
#' * Relocates misplaced label files when images are moved across splits (e.g. from train to val) to preserve annotations.
#' * Renames colliding images and review previews uniquely (`foo_jpg.jpg`, `foo_png.png`) and links existing `.txt` annotations to the best-matching file based on modification timestamps.
#' * Clamps out-of-bounds coordinates to `[0, 1]`.
#' * Prunes zero-area bounding boxes.
#' * Creates empty `.txt` background labels for unannotated images.
#' * Removes orphaned `.txt` label files.
#'
#' The returned object wraps the audited dataset seamlessly so it can be chained
#' directly with the pipe operator `|>` into `yolo_train()` or `yolo_dataset_augment()`.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param repair Logical. If `TRUE`, automatically fixes detected issues in-place. Defaults to `FALSE`.
#' @param min_images Integer threshold for warning if a class has few samples. Defaults to `10L`.
#' @param check_leakage Logical. Whether to check for near-duplicate / leaky images across splits using perceptual hashing. Defaults to `TRUE`.
#' @param leak_threshold Integer Hamming distance threshold for flagging near-duplicate images across splits. Defaults to `4L`.
#' @param verbose Logical. Whether to display progress and diagnostic alerts. Defaults to `TRUE`.
#'
#' @return An object of class `c("yolo_dataset_audit", "yolo_dataset", "list")` containing:
#' * `dir`: Root directory of the dataset.
#' * `is_healthy`: Logical indicator of dataset health.
#' * `issues_found`: Total count of issues detected.
#' * `issues_repaired`: Count of issues successfully fixed (when `repair = TRUE`).
#' * `details`: `data.frame` with file-by-file diagnostic records.
#' * `dataset`: Underlying `yolo_dataset` or `yolo_dataset_cls` object.
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Audit a YOLO dataset and automatically fix issues
#' ds_clean <- yolo_dataset_audit("dataset_cafe", repair = TRUE)
#'
#' # Pipe directly into training!
#' ds_clean |> yolo_train(model = "yolo26n-seg", epochs = 50)
yolo_dataset_audit <- function(data,
                               repair = FALSE,
                               min_images = 10L,
                               check_leakage = TRUE,
                               leak_threshold = 4L,
                               verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)
  n_classes <- length(classes)

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Health Audit (pliman)")
    cli::cli_alert_info("Directory: {.path {dir_path}}")
    cli::cli_alert_info("Detected task: {.val {task}} | Classes ({n_classes}): {.val {classes}}")
    if (isTRUE(repair)) cli::cli_alert_warning("Auto-repair mode: {.val ENABLED} (fixes will be applied in-place)")
  }

  issues <- list()
  repaired_count <- 0L

  # CASE A: CLASSIFICATION DATASET
  if (task == "classify") {
    splits <- c("train", "val", "test")
    splits <- splits[dir.exists(file.path(dir_path, splits))]

    img_exts <- "\\.(jpg|jpeg|png|bmp|webp)$"
    for (s in splits) {
      s_dir <- file.path(dir_path, s)
      cls_dirs <- list.dirs(s_dir, recursive = FALSE, full.names = FALSE)
      cls_dirs <- cls_dirs[nzchar(cls_dirs) & !startsWith(cls_dirs, ".")]
      for (cls in cls_dirs) {
        c_folder <- file.path(s_dir, cls)
        fls <- list.files(c_folder, pattern = img_exts, full.names = TRUE, ignore.case = TRUE)
        if (length(fls) == 0L) next
        szs <- file.size(fls)
        bad_idx <- which(is.na(szs) | szs == 0)
        if (length(bad_idx) > 0L) {
          bad_files <- fls[bad_idx]
          issues[[length(issues) + 1L]] <- data.frame(
            file = bad_files, type = "corrupt_image", message = "0-byte or unreadable image file", stringsAsFactors = FALSE
          )
          if (isTRUE(repair)) {
            unlink(bad_files)
            repaired_count <- repaired_count + length(bad_files)
          }
        }
      }
    }
  } else {
    # CASE B: DETECTION / SEGMENTATION DATASET
    img_dir <- file.path(dir_path, "images")
    lbl_dir <- file.path(dir_path, "labels")

    if (!dir.exists(img_dir)) {
      cli::cli_abort("Dataset is missing required {.path images/} directory.")
    }
    if (!dir.exists(lbl_dir)) {
      if (isTRUE(repair)) {
        dir.create(lbl_dir, recursive = TRUE, showWarnings = FALSE)
        repaired_count <- repaired_count + 1L
      } else {
        cli::cli_abort("Dataset is missing required {.path labels/} directory.")
      }
    }

    img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)

    # Check for corrupt (0-byte) images
    if (length(img_files) > 0L) {
      img_szs <- file.size(img_files)
      bad_img_idx <- which(is.na(img_szs) | img_szs == 0)
      if (length(bad_img_idx) > 0L) {
        bad_imgs <- img_files[bad_img_idx]
        issues[[length(issues) + 1L]] <- data.frame(
          file = bad_imgs, type = "corrupt_image", message = "0-byte or unreadable image file", stringsAsFactors = FALSE
        )
        if (isTRUE(repair)) {
          unlink(bad_imgs)
          repaired_count <- repaired_count + length(bad_imgs)
          img_files <- img_files[-bad_img_idx]
        }
      }
    }

    img_dir_norm <- gsub("\\\\", "/", normalizePath(img_dir, winslash = "/", mustWork = FALSE))
    lbl_dir_norm <- gsub("\\\\", "/", normalizePath(lbl_dir, winslash = "/", mustWork = FALSE))
    img_files_clean <- gsub("\\\\", "/", img_files)
    lbl_files_clean <- gsub("\\\\", "/", lbl_files)

    img_rel <- sub(paste0("^", img_dir_norm, "/?"), "", img_files_clean)
    lbl_rel <- sub(paste0("^", lbl_dir_norm, "/?"), "", lbl_files_clean)

    img_stems <- sub("\\.[a-zA-Z0-9]+$", "", img_rel)
    lbl_stems <- sub("\\.txt$", "", lbl_rel)

    # Check data.yaml integrity & valid image paths
    yaml_candidates <- c(file.path(dir_path, "data.yaml"), file.path(dir_path, "dataset.yaml"))
    yaml_exists <- yaml_candidates[file.exists(yaml_candidates)]
    if (length(yaml_exists) == 0L) {
      issues[[length(issues) + 1L]] <- data.frame(
        file = file.path(dir_path, "data.yaml"),
        type = "missing_data_yaml",
        message = "Dataset lacks a data.yaml descriptor configuration file",
        stringsAsFactors = FALSE
      )
      if (isTRUE(repair)) {
        .yolo_sync_data_yaml(dir_path, verbose = FALSE)
        repaired_count <- repaired_count + 1L
      }
    } else {
      yf <- yaml_exists[1]
      y_lines <- readLines(yf, warn = FALSE)
      p_idx <- grep("^\\s*path\\s*:", y_lines)
      v_idx <- grep("^\\s*val\\s*:", y_lines)
      t_idx <- grep("^\\s*train\\s*:", y_lines)
      p_val <- if (length(p_idx) > 0) gsub("['\"]", "", trimws(sub("^\\s*path\\s*:\\s*", "", y_lines[p_idx[1]]))) else dir_path
      if (!nzchar(p_val)) p_val <- dir_path
      broken <- FALSE
      broken_msg <- ""
      if (length(v_idx) > 0) {
        v_val <- gsub("['\"]", "", trimws(sub("^\\s*val\\s*:\\s*", "", y_lines[v_idx[1]])))
        full_v <- if (file.exists(v_val) || dir.exists(v_val)) v_val else file.path(p_val, v_val)
        if (!dir.exists(full_v) && !file.exists(full_v)) {
          broken <- TRUE
          broken_msg <- paste0("Validation path does not exist: '", full_v, "' (dataset was likely renamed or moved)")
        }
      }
      if (!broken && length(t_idx) > 0) {
        t_val <- gsub("['\"]", "", trimws(sub("^\\s*train\\s*:\\s*", "", y_lines[t_idx[1]])))
        full_t <- if (file.exists(t_val) || dir.exists(t_val)) t_val else file.path(p_val, t_val)
        if (!dir.exists(full_t) && !file.exists(full_t)) {
          broken <- TRUE
          broken_msg <- paste0("Training path does not exist: '", full_t, "' (dataset was likely renamed or moved)")
        }
      }
      if (broken) {
        issues[[length(issues) + 1L]] <- data.frame(
          file = yf,
          type = "stale_yaml_paths",
          message = paste0("data.yaml has stale paths: ", broken_msg),
          stringsAsFactors = FALSE
        )
        if (isTRUE(repair)) {
          .yolo_sync_data_yaml(dir_path, yaml_file = yf, verbose = FALSE)
          repaired_count <- repaired_count + 1L
        }
      }
    }

    # 1. Detect missing labels and orphan labels
    missing_lbl_stems <- setdiff(img_stems, lbl_stems)
    orphan_lbl_stems <- setdiff(lbl_stems, img_stems)
    n_misplaced <- 0L

    # 1a. Misplaced label files (e.g., image moved across split folders such as train -> val)
    # Check if any missing image label exists as an orphan label in another folder
    if (length(missing_lbl_stems) > 0L && length(orphan_lbl_stems) > 0L) {
      missing_base <- basename(missing_lbl_stems)
      orphan_base <- basename(orphan_lbl_stems)
      common_base <- intersect(missing_base, orphan_base)

      if (length(common_base) > 0L) {
        matched_missing <- character()
        matched_orphan <- character()

        for (cb in common_base) {
          mis_idx <- which(missing_base == cb)
          orph_idx <- which(orphan_base == cb)

          # If multiple orphan files share the same base name, prioritize non-empty and larger files
          if (length(orph_idx) > 1L) {
            orph_paths <- file.path(lbl_dir, paste0(orphan_lbl_stems[orph_idx], ".txt"))
            orph_sizes <- file.size(orph_paths)
            orph_sizes[is.na(orph_sizes)] <- 0
            orph_idx <- orph_idx[order(orph_sizes, decreasing = TRUE)]
          }

          n_pairs <- min(length(mis_idx), length(orph_idx))
          for (k in seq_len(n_pairs)) {
            matched_missing <- c(matched_missing, missing_lbl_stems[mis_idx[k]])
            matched_orphan <- c(matched_orphan, orphan_lbl_stems[orph_idx[k]])
          }
        }

        if (length(matched_missing) > 0L) {
          n_misplaced <- length(matched_missing)
          issues[[length(issues) + 1L]] <- data.frame(
            file = matched_missing,
            type = "misplaced_label",
            message = sprintf("Label found in '%s.txt' but image is in '%s'", matched_orphan, matched_missing),
            stringsAsFactors = FALSE
          )

          if (isTRUE(repair)) {
            for (k in seq_along(matched_missing)) {
              src_txt <- file.path(lbl_dir, paste0(matched_orphan[k], ".txt"))
              dest_txt <- file.path(lbl_dir, paste0(matched_missing[k], ".txt"))
              dest_dir <- dirname(dest_txt)
              if (!dir.exists(dest_dir)) {
                dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
              }
              ok <- file.rename(src_txt, dest_txt)
              if (!ok) {
                file.copy(src_txt, dest_txt, overwrite = TRUE)
                unlink(src_txt)
              }
            }
            repaired_count <- repaired_count + n_misplaced
          }

          # Remove matched files from missing and orphan lists so they are not treated as blank or deleted
          missing_lbl_stems <- setdiff(missing_lbl_stems, matched_missing)
          orphan_lbl_stems <- setdiff(orphan_lbl_stems, matched_orphan)
        }
      }
    }

    # 1b. Truly unannotated images (missing .txt anywhere in dataset)
    if (length(missing_lbl_stems) > 0L) {
      issues[[length(issues) + 1L]] <- data.frame(
        file = missing_lbl_stems, type = "missing_label", message = "Image lacks matching .txt label file", stringsAsFactors = FALSE
      )
      if (isTRUE(repair)) {
        dest_txts <- file.path(lbl_dir, paste0(missing_lbl_stems, ".txt"))
        unique_dirs <- unique(dirname(dest_txts))
        for (ud in unique_dirs) {
          if (!dir.exists(ud)) dir.create(ud, recursive = TRUE, showWarnings = FALSE)
        }
        file.create(dest_txts)
        repaired_count <- repaired_count + length(dest_txts)
      }
    }

    # 2. Truly orphan label files (no matching image anywhere in dataset)
    if (length(orphan_lbl_stems) > 0L) {
      issues[[length(issues) + 1L]] <- data.frame(
        file = orphan_lbl_stems, type = "orphan_label", message = "Label .txt exists but image is missing", stringsAsFactors = FALSE
      )
      if (isTRUE(repair)) {
        txt_paths <- file.path(lbl_dir, paste0(orphan_lbl_stems, ".txt"))
        unlink(txt_paths[file.exists(txt_paths)])
        repaired_count <- repaired_count + length(orphan_lbl_stems)
      }
    }

    # 2b. Stem collisions (multiple images sharing the exact same stem with different extensions)
    if (any(duplicated(img_stems))) {
      dup_stems <- unique(img_stems[duplicated(img_stems)])
      for (ds in dup_stems) {
        colliding_files <- img_rel[img_stems == ds]
        issues[[length(issues) + 1L]] <- data.frame(
          file = ds, type = "stem_collision",
          message = sprintf("Multiple images share base name '%s', colliding on single label file: %s",
                            basename(ds), paste(basename(colliding_files), collapse = ", ")),
          stringsAsFactors = FALSE
        )

        if (isTRUE(repair)) {
          lbl_file <- file.path(lbl_dir, paste0(ds, ".txt"))
          img_paths <- file.path(img_dir, colliding_files)
          best_idx <- 1L
          if (file.exists(lbl_file)) {
            lbl_mtime <- file.info(lbl_file)$mtime
            img_mtimes <- file.info(img_paths)$mtime
            diffs <- abs(as.numeric(difftime(img_mtimes, lbl_mtime, units = "secs")))
            diffs[is.na(diffs)] <- Inf
            best_idx <- which.min(diffs)
          }

          base_s <- basename(ds)
          for (k in seq_along(colliding_files)) {
            old_rel <- colliding_files[k]
            ext <- tools::file_ext(old_rel)
            split_dir <- dirname(old_rel)
            new_s <- paste0(base_s, "_", tolower(ext))
            new_img_name <- paste0(new_s, ".", ext)

            old_img_path <- file.path(img_dir, old_rel)
            new_img_path <- file.path(img_dir, split_dir, new_img_name)

            if (file.exists(old_img_path) && old_img_path != new_img_path) {
              file.rename(old_img_path, new_img_path)
            }

            rev_dir <- file.path(dir_path, "review")
            if (dir.exists(rev_dir)) {
              old_rev <- file.path(rev_dir, paste0("review_", basename(old_rel)))
              new_rev <- file.path(rev_dir, paste0("review_", new_img_name))
              if (file.exists(old_rev) && old_rev != new_rev) {
                file.rename(old_rev, new_rev)
              }
            }
          }

          if (file.exists(lbl_file)) {
            matched_ext <- tools::file_ext(colliding_files[best_idx])
            matched_new_lbl <- file.path(lbl_dir, dirname(colliding_files[best_idx]),
                                         paste0(base_s, "_", tolower(matched_ext), ".txt"))
            file.rename(lbl_file, matched_new_lbl)
          }

          for (k in seq_along(colliding_files)) {
            ext_k <- tools::file_ext(colliding_files[k])
            lbl_k <- file.path(lbl_dir, dirname(colliding_files[k]), paste0(base_s, "_", tolower(ext_k), ".txt"))
            if (!file.exists(lbl_k)) {
              file.create(lbl_k)
            }
          }
          repaired_count <- repaired_count + 1L
        }
      }

      if (isTRUE(repair)) {
        img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
        lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
        img_files_clean <- gsub("\\\\", "/", img_files)
        lbl_files_clean <- gsub("\\\\", "/", lbl_files)
        img_rel <- sub(paste0("^", img_dir_norm, "/?"), "", img_files_clean)
        lbl_rel <- sub(paste0("^", lbl_dir_norm, "/?"), "", lbl_files_clean)
        img_stems <- sub("\\.[a-zA-Z0-9]+$", "", img_rel)
        lbl_stems <- sub("\\.txt$", "", lbl_rel)
      }
    }

    # Refresh label files if any repairs occurred
    if (isTRUE(repair) && (n_misplaced > 0L || length(missing_lbl_stems) > 0L || length(orphan_lbl_stems) > 0L || any(duplicated(img_stems)))) {
      lbl_files <- list.files(lbl_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
    }

    # 3. Label file syntactic and coordinate checks (C++ OpenMP)
    total_objs_found <- 0L
    if (length(lbl_files) > 0L) {
      audit_res <- fast_audit_yolo_labels_cpp(
        lbl_files = lbl_files,
        n_classes = n_classes,
        repair = isTRUE(repair),
        num_threads = 0L
      )
      repaired_count <- repaired_count + audit_res$repaired_count
      total_objs_found <- audit_res$n_objects

      if (length(audit_res$files) > 0L) {
        issues[[length(issues) + 1L]] <- data.frame(
          file = audit_res$files,
          type = audit_res$types,
          message = audit_res$messages,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Check Cross-Split Data Leakage via fast C++ dHash (64-bit perceptual hash)
  if (isTRUE(check_leakage) && task != "classify") {
    img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
    all_sp <- c("train", "val", "test")
    all_im_files <- character(0)
    all_im_splits <- character(0)
    for (sp in all_sp) {
      sp_dir <- file.path(dir_path, "images", sp)
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp, "images")
      if (dir.exists(sp_dir)) {
        fls <- list.files(sp_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
        if (length(fls) > 0L) {
          all_im_files <- c(all_im_files, fls)
          all_im_splits <- c(all_im_splits, rep(sp, length(fls)))
        }
      }
    }

    has_tr <- any(all_im_splits == "train")
    has_ev <- any(all_im_splits %in% c("val", "test"))

    if (has_tr && has_ev && length(all_im_files) >= 2L) {
      leak_res <- tryCatch({
        fast_dhash_leakage_cpp(all_im_files, all_im_splits, max_distance = as.integer(leak_threshold), num_threads = 0L)
      }, error = function(e) NULL)

      if (!is.null(leak_res) && leak_res$leak_count > 0L) {
        leaks_df <- leak_res$leaks
        for (lk_i in seq_len(nrow(leaks_df))) {
          ev_f <- leaks_df$img_eval[lk_i]
          tr_f <- leaks_df$img_train[lk_i]
          sp_e <- leaks_df$split_eval[lk_i]
          dist_v <- leaks_df$hamming_dist[lk_i]
          sim_v <- leaks_df$similarity_pct[lk_i]

          issues[[length(issues) + 1L]] <- data.frame(
            file = ev_f,
            type = "data_leakage",
            message = sprintf("Near-duplicate of train/%s in %s/ (Hamming dist: %d, %.1f%% match) - artificially inflates mAP!",
                              basename(tr_f), sp_e, dist_v, sim_v),
            stringsAsFactors = FALSE
          )

          if (isTRUE(repair)) {
            unlink(ev_f, force = TRUE)
            ev_lbl <- file.path(dir_path, "labels", sp_e, paste0(tools::file_path_sans_ext(basename(ev_f)), ".txt"))
            if (file.exists(ev_lbl)) unlink(ev_lbl, force = TRUE)
            repaired_count <- repaired_count + 1L
          }
        }
      }
    }
  }

  issues_df <- if (length(issues) > 0L) do.call(rbind, issues) else data.frame(file = character(), type = character(), message = character(), stringsAsFactors = FALSE)
  is_healthy <- (nrow(issues_df) == 0L)

  if (isTRUE(verbose)) {
    if (is_healthy) {
      cli::cli_alert_success("Dataset passed all audit checks! No anomalies or corruption detected.")
    } else {
      cli::cli_alert_danger("Found {nrow(issues_df)} issue(s) during dataset audit.")
      tab_types <- table(issues_df$type)
      for (nm in names(tab_types)) {
        cli::cli_bullets(c("*" = "{nm}: {tab_types[[nm]]} occurrences"))
      }
      if (isTRUE(repair)) {
        cli::cli_alert_success("Auto-repair applied: successfully resolved {repaired_count} issue(s)!")
      } else {
        cli::cli_alert_info("Run with {.code repair = TRUE} to automatically heal these issues.")
      }
    }
  }

  wrapped_ds <- .yolo_wrap_dataset(dir_path, task = task, classes = classes, n_objects = if (task == "classify") NULL else total_objs_found)

  out <- list(
    dir = dir_path,
    yaml_file = wrapped_ds$yaml_file,
    data_yaml = wrapped_ds$data_yaml,
    classes = classes,
    task = task,
    splits = wrapped_ds$splits,
    summary = wrapped_ds$summary,
    total_images = wrapped_ds$total_images,
    train_images = wrapped_ds$train_images,
    val_images = wrapped_ds$val_images,
    test_images = wrapped_ds$test_images,
    n_objects = wrapped_ds$n_objects,
    is_healthy = is_healthy,
    issues_found = nrow(issues_df),
    issues_repaired = repaired_count,
    details = issues_df,
    dataset = wrapped_ds
  )
  class(out) <- unique(c("yolo_dataset_audit", class(wrapped_ds)))
  invisible(out)
}

#' @export
yolo_dataset_health <- yolo_dataset_audit

#' @export
print.yolo_dataset_audit <- function(x, ...) {
  cli::cli_h3("YOLO Dataset Audit Report (pliman)")
  if (x$is_healthy) {
    cli::cli_alert_success("Health Status: EXCELLENT (100% compliant, ready for training)")
  } else {
    cli::cli_alert_warning("Health Status: ISSUES DETECTED ({x$issues_found} anomalies)")
    if (x$issues_repaired > 0L) {
      cli::cli_alert_success("Issues Repaired: {x$issues_repaired}")
    }
    if (nrow(x$details) > 0L) {
      print(utils::head(x$details, 10L), row.names = FALSE)
      if (nrow(x$details) > 10L) {
        cli::cli_alert_info("... and {nrow(x$details) - 10L} more issue(s).")
      }
    }
  }
  if (!is.null(x$dataset)) {
    print(x$dataset, ...)
  }
  invisible(x)
}

#' @export
plot.yolo_dataset_audit <- function(x, ...) {
  if (!is.null(x$dataset)) {
    plot(x$dataset, ...)
  } else {
    plot(x$dir, ...)
  }
}


# ==============================================================================
# 2b. PERCEPTUAL DEDUPLICATION & CROSS-SPLIT LEAKAGE DETECTION (yolo_dataset_dedup)
# ==============================================================================
#' @title Perceptual Deduplication and Cross-Split Data Leakage Detection
#' @name yolo_dataset_dedup
#' @aliases yolo_dataset_deduplicate
#' @description
#' `yolo_dataset_dedup()` computes ultrarrápido 64-bit Perceptual Difference Hashes (dHash)
#' using native OpenMP C++ to identify near-duplicate images within and across YOLO dataset splits.
#'
#' Crucially, it detects **Cross-Split Data Leakage** (where near-identical frames or drone photos
#' leak between `train/` and `val/` or `test/`), which artificially inflates validation metrics (mAP).
#'
#' @param data A `yolo_dataset` object or directory path.
#' @param threshold Integer (0-64). Maximum Hamming distance to consider two images as near-duplicates
#'   (default: `4L`, corresponding to \eqn{\ge 93.8\%} perceptual visual similarity). Distance 0 denotes exact duplicates.
#' @param action Action to execute: `"report"` (default, displays analysis without altering files),
#'   `"quarantine"` (moves leaking or duplicate images and labels to a separate directory), or
#'   `"remove"` (deletes duplicate/leaking files from the target split).
#' @param check_leakage Logical. Whether to inspect and flag cross-split leakage (`train` vs `val`/`test`). Defaults to `TRUE`.
#' @param target_split Character specifying which split to act upon when `action = "quarantine"` or `"remove"`:
#'   `"val"` (default, cleans validation leaks preserving train), `"test"`, `"train"`, or `"both"`.
#' @param quarantine_dir Optional directory path for quarantined files. Defaults to `<dataset>/quarantine`.
#' @param num_threads Integer. Number of OpenMP threads for C++ execution (default: `0L`, auto-detect all CPU cores).
#' @param verbose Logical. Whether to show progress messages (default: `TRUE`).
#'
#' @return An invisible list containing `leaks` (data frame of cross-split leaks), `duplicates` (data frame of in-split duplicates),
#'   and `hashes` (named character vector of 64-bit hexadecimal hashes).
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Inspect dataset for duplicates and cross-split leakage
#' res <- yolo_dataset_dedup("my_dataset")
#'
#' # 2. Quarantine near-duplicate validation images leaking from train
#' yolo_dataset_dedup("my_dataset", action = "quarantine")
yolo_dataset_dedup <- function(data,
                               threshold = 4L,
                               action = c("report", "quarantine", "remove"),
                               check_leakage = TRUE,
                               target_split = c("val", "test", "train", "both"),
                               quarantine_dir = NULL,
                               num_threads = 0L,
                               verbose = TRUE) {
  action <- match.arg(action)
  target_split <- match.arg(target_split)
  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")

  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  all_splits <- c("train", "val", "test")
  all_imgs <- character(0)
  all_sp <- character(0)

  if (task == "classify") {
    for (sp in all_splits) {
      sp_dir <- file.path(dir_path, sp)
      if (dir.exists(sp_dir)) {
        fls <- list.files(sp_dir, pattern = img_pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
        if (length(fls) > 0L) {
          all_imgs <- c(all_imgs, fls)
          all_sp <- c(all_sp, rep(sp, length(fls)))
        }
      }
    }
  } else {
    for (sp in all_splits) {
      sp_dir <- file.path(dir_path, "images", sp)
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp, "images")
      if (!dir.exists(sp_dir)) sp_dir <- file.path(dir_path, sp)
      if (dir.exists(sp_dir)) {
        fls <- list.files(sp_dir, pattern = img_pattern, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
        if (length(fls) > 0L) {
          all_imgs <- c(all_imgs, fls)
          all_sp <- c(all_sp, rep(sp, length(fls)))
        }
      }
    }
    if (length(all_imgs) == 0L) {
      all_imgs <- list.files(dir_path, pattern = img_pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      all_sp <- rep("train", length(all_imgs))
    }
  }

  n_imgs <- length(all_imgs)
  if (n_imgs < 2L) {
    if (isTRUE(verbose)) cli::cli_alert_info("Dataset contains fewer than 2 images. No deduplication needed.")
    return(invisible(list(leaks = data.frame(), duplicates = data.frame(), hashes = character(0))))
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Perceptual Deduplication & Leakage Detection (pliman)")
    cli::cli_alert_info("Analyzing {n_imgs} images via 64-bit dHash (Hamming threshold <= {threshold})...")
  }

  res_cpp <- fast_dhash_leakage_cpp(
    img_files = all_imgs,
    splits = all_sp,
    max_distance = as.integer(threshold),
    num_threads = as.integer(num_threads)
  )

  leaks_df <- res_cpp$leaks
  dups_df <- res_cpp$duplicates

  if (isTRUE(verbose)) {
    if (nrow(leaks_df) == 0L) {
      cli::cli_alert_success("Cross-Split Data Leakage: NONE detected! Train, Val and Test splits are strictly isolated.")
    } else {
      cli::cli_alert_danger("Cross-Split Data Leakage: FOUND {nrow(leaks_df)} leaking pair(s) between train and eval splits!")
      print(utils::head(leaks_df, 10L), row.names = FALSE)
      if (nrow(leaks_df) > 10L) cli::cli_alert_info("... and {nrow(leaks_df) - 10L} more leaking pair(s).")
    }

    if (nrow(dups_df) == 0L) {
      cli::cli_alert_success("In-Split Duplicate Images: NONE detected.")
    } else {
      cli::cli_alert_warning("In-Split Duplicate Images: FOUND {nrow(dups_df)} duplicate or near-duplicate pair(s)!")
      print(utils::head(dups_df, 10L), row.names = FALSE)
      if (nrow(dups_df) > 10L) cli::cli_alert_info("... and {nrow(dups_df) - 10L} more duplicate pair(s).")
    }
  }

  # Execute actions: quarantine or remove
  if (action %in% c("quarantine", "remove") && (nrow(leaks_df) > 0L || nrow(dups_df) > 0L)) {
    to_act <- character(0)

    if (isTRUE(check_leakage) && nrow(leaks_df) > 0L) {
      if (target_split %in% c("val", "test")) {
        to_act <- c(to_act, leaks_df$img_eval[leaks_df$split_eval == target_split])
      } else if (target_split == "train") {
        to_act <- c(to_act, leaks_df$img_train)
      } else { # both
        to_act <- c(to_act, leaks_df$img_eval)
      }
    }

    if (nrow(dups_df) > 0L) {
      to_act <- c(to_act, dups_df$img_2)
    }

    to_act <- unique(to_act[file.exists(to_act)])

    if (length(to_act) > 0L) {
      if (action == "quarantine") {
        q_dir <- if (!is.null(quarantine_dir)) quarantine_dir else file.path(dir_path, "quarantine")
        dir.create(q_dir, recursive = TRUE, showWarnings = FALSE)
        for (im_path in to_act) {
          file.rename(im_path, file.path(q_dir, basename(im_path)))
          stem <- tools::file_path_sans_ext(basename(im_path))
          lbl_path <- file.path(dirname(dirname(im_path)), "labels", basename(dirname(im_path)), paste0(stem, ".txt"))
          if (file.exists(lbl_path)) file.rename(lbl_path, file.path(q_dir, paste0(stem, ".txt")))
        }
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Quarantined {length(to_act)} image and label file(s) to {.path {q_dir}}.")
        }
      } else if (action == "remove") {
        for (im_path in to_act) {
          unlink(im_path, force = TRUE)
          stem <- tools::file_path_sans_ext(basename(im_path))
          lbl_path <- file.path(dirname(dirname(im_path)), "labels", basename(dirname(im_path)), paste0(stem, ".txt"))
          if (file.exists(lbl_path)) unlink(lbl_path, force = TRUE)
        }
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Removed {length(to_act)} leaking/duplicate image and label file(s).")
        }
      }
      unlink(list.files(dir_path, pattern = "\\.cache$", full.names = TRUE, recursive = TRUE), force = TRUE)
    }
  }

  invisible(list(
    leaks = leaks_df,
    duplicates = dups_df,
    hashes = setNames(res_cpp$hashes, all_imgs)
  ))
}

#' @rdname yolo_dataset_dedup
#' @export
yolo_dataset_deduplicate <- yolo_dataset_dedup


# ==============================================================================
# 3. DATASET MERGING & CLASS CONFLICT RESOLUTION (yolo_dataset_merge)
# ==============================================================================
#' @title Merge Multiple YOLO Datasets with Class Conflict Resolution
#' @name yolo_dataset_merge
#' @aliases yolo_dataset_combine
#' @description
#' `yolo_dataset_merge()` combines two or more YOLO datasets (Detection, Segmentation,
#' or Classification) into a single, cohesive, conflict-free master dataset.
#' It handles filename collisions by prepending dataset prefixes, harmonizes class indices
#' across differing `data.yaml` definitions (with optional custom class renaming via `class_map`),
#' updates label files on-the-fly, and generates a unified master `data.yaml`.
#'
#' @param datasets A list or character vector of dataset directory paths or `yolo_dataset` objects.
#' @param out_dir Output directory path for the merged dataset. Defaults to `"merged_yolo_dataset"`.
#' @param class_map Optional named list or vector mapping old class names to unified names
#'   (e.g., `list("cerca" = "fence", "folha" = "leaf")`).
#' @param copy Logical. If `TRUE` (default), copies files. If `FALSE`, moves files using [file.rename()].
#' @param overwrite Logical. If `TRUE`, overwrites `out_dir` if it exists. Defaults to `FALSE`.
#' @param parallel Logical. Whether to use parallel multithreading for file copying and label remapping (default: `TRUE`).
#' @param num_threads Integer. Number of OpenMP threads to use (default: `0`, uses all available CPU cores).
#' @param verbose Logical. Whether to show progress. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` object pointing to `out_dir` (ready for `yolo_train()` or piping `|>`).
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Merge two datasets collected from different seasons
#' ds_total <- yolo_dataset_merge(
#'   datasets = list("dataset_2024", "dataset_2025"),
#'   out_dir  = "dataset_unificado"
#' )
yolo_dataset_merge <- function(datasets,
                               out_dir = "merged_yolo_dataset",
                               class_map = NULL,
                               copy = TRUE,
                               overwrite = FALSE,
                               parallel = TRUE,
                               num_threads = 0,
                               verbose = TRUE) {
  if (length(datasets) < 2L) {
    cli::cli_abort("At least 2 datasets must be provided in {.arg datasets} to merge.")
  }

  ds_dirs <- vapply(datasets, .yolo_get_dir, character(1))
  out_dir_norm <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)

  if (dir.exists(out_dir_norm)) {
    if (isTRUE(overwrite)) {
      unlink(out_dir_norm, recursive = TRUE, force = TRUE)
    } else {
      cli::cli_abort("Output directory {.path {out_dir_norm}} already exists. Use {.code overwrite = TRUE} to replace.")
    }
  }
  dir.create(out_dir_norm, recursive = TRUE, showWarnings = FALSE)

  # Determine number of threads
  threads <- if (isTRUE(parallel)) as.integer(num_threads)[1] else 1L
  if (is.na(threads) || threads < 0L) threads <- 0L

  # Detect tasks across datasets
  tasks <- vapply(ds_dirs, .yolo_detect_task, character(1))
  is_all_cls <- all(tasks == "classify")
  is_seg <- any(tasks == "segment")
  master_task <- if (is_all_cls) "classify" else if (is_seg) "segment" else "detect"

  # Gather and unify classes
  all_local_classes <- lapply(ds_dirs, .yolo_get_classes)
  master_classes <- character()

  for (i in seq_along(all_local_classes)) {
    cls_vec <- all_local_classes[[i]]
    if (!is.null(class_map) && length(class_map) > 0L) {
      for (k in seq_along(cls_vec)) {
        if (cls_vec[k] %in% names(class_map)) {
          cls_vec[k] <- as.character(class_map[[cls_vec[k]]])
        }
      }
    }
    all_local_classes[[i]] <- cls_vec
    master_classes <- union(master_classes, cls_vec)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Merging Engine (pliman)")
    cli::cli_alert_info("Merging {length(ds_dirs)} datasets into {.path {out_dir_norm}}")
    cli::cli_alert_info("Unified Task: {.val {master_task}} | Unified Classes ({length(master_classes)}): {.val {master_classes}}")
    if (threads != 1L) {
      cli::cli_alert_info("Parallel acceleration: active with OpenMP multithreading.")
    }
  }

  total_images_merged <- 0L
  total_objects_merged <- 0L

  # CASE A: CLASSIFICATION MERGE
  if (master_task == "classify") {
    splits <- c("train", "val", "test")
    for (i in seq_along(ds_dirs)) {
      src_d <- ds_dirs[i]
      prefix <- sprintf("ds%d_", i)
      loc_classes <- all_local_classes[[i]]

      for (s in splits) {
        s_src <- file.path(src_d, s)
        if (!dir.exists(s_src)) next
        for (lc in list.dirs(s_src, recursive = FALSE, full.names = FALSE)) {
          target_c <- if (!is.null(class_map) && lc %in% names(class_map)) class_map[[lc]] else lc
          dest_folder <- file.path(out_dir_norm, s, target_c)
          dir.create(dest_folder, recursive = TRUE, showWarnings = FALSE)

          fls <- list.files(file.path(s_src, lc), pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
          if (length(fls) == 0L) next
          new_fns <- paste0(prefix, basename(fls))
          dest_fls <- file.path(dest_folder, new_fns)
          parallel_copy_files_cpp(fls, dest_fls, copy = copy, num_threads = threads)
          total_images_merged <- total_images_merged + length(fls)
        }
      }
    }
    .yolo_write_yaml(out_dir_norm, classes = master_classes, task = "classify")
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Merged classification dataset written to {.path {out_dir_norm}} ({total_images_merged} images).")
    }
    return(.yolo_wrap_dataset(out_dir_norm, task = "classify", classes = master_classes))
  }

  # CASE B: DETECTION / SEGMENTATION MERGE
  splits <- c("train", "val", "test")
  for (s in splits) {
    dir.create(file.path(out_dir_norm, "images", s), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir_norm, "labels", s), recursive = TRUE, showWarnings = FALSE)
  }

  for (i in seq_along(ds_dirs)) {
    src_d <- ds_dirs[i]
    prefix <- sprintf("ds%d_", i)
    loc_classes <- all_local_classes[[i]]

    if (isTRUE(verbose)) {
      cli::cli_alert_info("Processing dataset [{i}/{length(ds_dirs)}]: {.file {basename(src_d)}}...")
    }

    # Vector for class ID mapping (0-indexed)
    id_map <- integer(length(loc_classes))
    for (k in seq_along(loc_classes)) {
      m_idx <- which(master_classes == loc_classes[k]) - 1L
      id_map[k] <- m_idx
    }

    for (s in splits) {
      src_img_dir <- file.path(src_d, "images", s)
      src_lbl_dir <- file.path(src_d, "labels", s)
      dest_img_dir <- file.path(out_dir_norm, "images", s)
      dest_lbl_dir <- file.path(out_dir_norm, "labels", s)

      if (!dir.exists(src_img_dir)) next
      img_files <- list.files(src_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
      if (length(img_files) == 0L) next

      stems <- sub("\\.[a-zA-Z0-9]+$", "", basename(img_files))
      exts <- tools::file_ext(img_files)
      new_fn <- paste0(prefix, stems, ".", exts)
      dest_im_paths <- file.path(dest_img_dir, new_fn)

      src_lbl_paths <- file.path(src_lbl_dir, paste0(stems, ".txt"))
      dest_txt_paths <- file.path(dest_lbl_dir, paste0(prefix, stems, ".txt"))

      # 1. High-speed parallel copy/move of images via OpenMP
      parallel_copy_files_cpp(img_files, dest_im_paths, copy = copy, num_threads = threads)

      # 2. High-speed parallel remapping and writing of labels via OpenMP
      n_objs <- parallel_remap_labels_cpp(src_lbl_paths, dest_txt_paths, id_map, num_threads = threads)

      total_images_merged <- total_images_merged + length(img_files)
      total_objects_merged <- total_objects_merged + n_objs
    }
  }

  .yolo_write_yaml(out_dir_norm, classes = master_classes, task = master_task)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Merged dataset written to {.path {out_dir_norm}} ({total_images_merged} images, {total_objects_merged} objects).")
  }

  .yolo_wrap_dataset(out_dir_norm, task = master_task, classes = master_classes)
}

#' @export
yolo_dataset_combine <- yolo_dataset_merge


# ------------------------------------------------------------------------------
# Internal helper: Permute target object cluster and fill background
# ------------------------------------------------------------------------------
.permute_cluster_single <- function(im,
                                    raw_lines = character(0),
                                    position = "random",
                                    fill = c("border_slices", "mean", "blur", "noise"),
                                    margin = 0.02,
                                    blend_edge = 5L) {
  fill <- if (length(fill) > 1L) fill[1] else fill
  if (!fill %in% c("border_slices", "mean", "blur", "noise")) fill <- "border_slices"

  u <- unclass(im)
  arr <- if (is.raw(u)) as.numeric(u) / 255.0 else as.numeric(u)
  dim(arr) <- dim(u)
  if (length(dim(arr)) == 2L) dim(arr) <- c(dim(arr)[1], dim(arr)[2], 1L)

  W <- dim(arr)[1]
  H <- dim(arr)[2]
  n_channels <- dim(arr)[3]

  parsed_boxes <- list()
  parsed_items <- list()

  if (length(raw_lines) > 0L) {
    for (ln in raw_lines) {
      vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
      if (length(vals) < 5L || any(is.na(vals))) next
      cid <- as.integer(vals[1])
      if (length(vals) == 5L) {
        xc <- vals[2]; yc <- vals[3]; bw <- vals[4]; bh <- vals[5]
        parsed_boxes[[length(parsed_boxes) + 1L]] <- c(xc - bw * 0.5, xc + bw * 0.5, yc - bh * 0.5, yc + bh * 0.5)
        parsed_items[[length(parsed_items) + 1L]] <- list(type = "box", cid = cid, xc = xc, yc = yc, w = bw, h = bh)
      } else if (length(vals) > 5L && (length(vals) - 5L) %% 3L == 0L) {
        xc <- vals[2]; yc <- vals[3]; bw <- vals[4]; bh <- vals[5]
        k_raw <- vals[6:length(vals)]
        parsed_boxes[[length(parsed_boxes) + 1L]] <- c(xc - bw * 0.5, xc + bw * 0.5, yc - bh * 0.5, yc + bh * 0.5)
        parsed_items[[length(parsed_items) + 1L]] <- list(type = "pose", cid = cid, xc = xc, yc = yc, w = bw, h = bh, kpts = k_raw)
      } else {
        coords <- vals[-1]
        xs <- coords[seq(1, length(coords), by = 2L)]
        ys <- coords[seq(2, length(coords), by = 2L)]
        parsed_boxes[[length(parsed_boxes) + 1L]] <- c(min(xs), max(xs), min(ys), max(ys))
        parsed_items[[length(parsed_items) + 1L]] <- list(type = "poly", cid = cid, xs = xs, ys = ys)
      }
    }
  }

  if (length(parsed_boxes) == 0L) {
    return(list(image = im, labels = raw_lines, shifted = FALSE, pos = "none"))
  }

  all_x1 <- vapply(parsed_boxes, function(b) b[1], numeric(1))
  all_x2 <- vapply(parsed_boxes, function(b) b[2], numeric(1))
  all_y1 <- vapply(parsed_boxes, function(b) b[3], numeric(1))
  all_y2 <- vapply(parsed_boxes, function(b) b[4], numeric(1))

  c_xmin <- max(0.0, min(all_x1) - margin)
  c_xmax <- min(1.0, max(all_x2) + margin)
  c_ymin <- max(0.0, min(all_y1) - margin)
  c_ymax <- min(1.0, max(all_y2) + margin)

  c_w <- c_xmax - c_xmin
  c_h <- c_ymax - c_ymin

  if (c_w >= 0.92 && c_h >= 0.92) {
    return(list(image = im, labels = raw_lines, shifted = FALSE, pos = "none"))
  }

  positions_all <- c("left", "right", "top", "bottom", "top_left", "top_right", "bottom_left", "bottom_right")
  chosen_pos <- if (length(position) > 1L) sample(position, 1L) else position
  if (chosen_pos == "random") {
    cur_cx <- (c_xmin + c_xmax) * 0.5
    cur_cy <- (c_ymin + c_ymax) * 0.5
    candidates <- character(0)
    if (cur_cx > 0.35) candidates <- c(candidates, "left", "top_left", "bottom_left")
    if (cur_cx < 0.65) candidates <- c(candidates, "right", "top_right", "bottom_right")
    if (cur_cy > 0.35) candidates <- c(candidates, "top", "top_left", "top_right")
    if (cur_cy < 0.65) candidates <- c(candidates, "bottom", "bottom_left", "bottom_right")
    if (length(candidates) == 0L) candidates <- positions_all
    chosen_pos <- sample(candidates, 1L)
  } else if (chosen_pos == "corner") {
    chosen_pos <- sample(c("top_left", "top_right", "bottom_left", "bottom_right"), 1L)
  }

  target_x1 <- switch(chosen_pos,
    "left" = 0.0,
    "top_left" = 0.0,
    "bottom_left" = 0.0,
    "right" = max(0.0, 1.0 - c_w),
    "top_right" = max(0.0, 1.0 - c_w),
    "bottom_right" = max(0.0, 1.0 - c_w),
    c_xmin
  )

  target_y1 <- switch(chosen_pos,
    "top" = 0.0,
    "top_left" = 0.0,
    "top_right" = 0.0,
    "bottom" = max(0.0, 1.0 - c_h),
    "bottom_left" = max(0.0, 1.0 - c_h),
    "bottom_right" = max(0.0, 1.0 - c_h),
    c_ymin
  )

  src_px1 <- max(1L, floor(c_xmin * W) + 1L)
  src_px2 <- min(W, ceiling(c_xmax * W))
  src_py1 <- max(1L, floor(c_ymin * H) + 1L)
  src_py2 <- min(H, ceiling(c_ymax * H))

  patch <- arr[src_px1:src_px2, src_py1:src_py2, , drop = FALSE]
  patch_w <- dim(patch)[1]
  patch_h <- dim(patch)[2]

  dst_px1 <- max(1L, min(W - patch_w + 1L, floor(target_x1 * W) + 1L))
  dst_px2 <- dst_px1 + patch_w - 1L
  dst_py1 <- max(1L, min(H - patch_h + 1L, floor(target_y1 * H) + 1L))
  dst_py2 <- dst_py1 + patch_h - 1L

  actual_shift_x <- (dst_px1 - src_px1) / W
  actual_shift_y <- (dst_py1 - src_py1) / H

  if (abs(actual_shift_x) < 0.02 && abs(actual_shift_y) < 0.02) {
    if (src_px1 > (W - src_px2)) {
      dst_px1 <- 1L
      dst_px2 <- dst_px1 + patch_w - 1L
      chosen_pos <- "left"
    } else {
      dst_px2 <- W
      dst_px1 <- dst_px2 - patch_w + 1L
      chosen_pos <- "right"
    }
    actual_shift_x <- (dst_px1 - src_px1) / W
    actual_shift_y <- (dst_py1 - src_py1) / H
  }

  buf <- max(2L, min(15L, as.integer(round(min(W, H) * 0.01))))
  border_rects <- list()
  if (src_px1 > (buf + 10L)) border_rects[[length(border_rects) + 1L]] <- list(x1 = 1L, x2 = src_px1 - buf, y1 = 1L, y2 = H)
  if ((W - src_px2) >= (buf + 10L)) border_rects[[length(border_rects) + 1L]] <- list(x1 = src_px2 + buf, x2 = W, y1 = 1L, y2 = H)
  if (src_py1 > (buf + 10L)) border_rects[[length(border_rects) + 1L]] <- list(x1 = 1L, x2 = W, y1 = 1L, y2 = src_py1 - buf)
  if ((H - src_py2) >= (buf + 10L)) border_rects[[length(border_rects) + 1L]] <- list(x1 = 1L, x2 = W, y1 = src_py2 + buf, y2 = H)

  new_arr <- arr

  if (fill == "border_slices" && length(border_rects) > 0L) {
    m_r <- mean(vapply(border_rects, function(r) mean(arr[r$x1:r$x2, r$y1:r$y2, 1L]), numeric(1)))
    m_g <- if (n_channels >= 3L) mean(vapply(border_rects, function(r) mean(arr[r$x1:r$x2, r$y1:r$y2, 2L]), numeric(1))) else m_r
    m_b <- if (n_channels >= 3L) mean(vapply(border_rects, function(r) mean(arr[r$x1:r$x2, r$y1:r$y2, 3L]), numeric(1))) else m_r

    hole_w <- src_px2 - src_px1 + 1L
    hole_h <- src_py2 - src_py1 + 1L
    noise_base <- rnorm(hole_w * hole_h, mean = 0, sd = 0.02)

    new_arr[src_px1:src_px2, src_py1:src_py2, 1L] <- pmax(0.0, pmin(1.0, m_r + noise_base))
    if (n_channels >= 3L) {
      new_arr[src_px1:src_px2, src_py1:src_py2, 2L] <- pmax(0.0, pmin(1.0, m_g + noise_base))
      new_arr[src_px1:src_px2, src_py1:src_py2, 3L] <- pmax(0.0, pmin(1.0, m_b + noise_base))
    }

    step_x <- max(16L, min(100L, hole_w %/% 3L))
    step_y <- max(16L, min(100L, hole_h %/% 3L))

    xs_grid <- seq(src_px1, src_px2, by = max(8L, step_x %/% 2L))
    ys_grid <- seq(src_py1, src_py2, by = max(8L, step_y %/% 2L))

    for (gx in xs_grid) {
      for (gy in ys_grid) {
        br <- border_rects[[sample(seq_along(border_rects), 1L)]]
        avail_w <- br$x2 - br$x1 + 1L
        avail_h <- br$y2 - br$y1 + 1L

        sw <- min(avail_w, sample(seq(step_x, min(avail_w, step_x * 2L)), 1L))
        sh <- min(avail_h, sample(seq(step_y, min(avail_h, step_y * 2L)), 1L))

        sx1 <- if (avail_w > sw) sample(br$x1:(br$x2 - sw + 1L), 1L) else br$x1
        sy1 <- if (avail_h > sh) sample(br$y1:(br$y2 - sh + 1L), 1L) else br$y1

        slice_patch <- arr[sx1:(sx1 + sw - 1L), sy1:(sy1 + sh - 1L), , drop = FALSE]

        if (runif(1) > 0.5) slice_patch <- slice_patch[rev(seq_len(sw)), , , drop = FALSE]
        if (runif(1) > 0.5) slice_patch <- slice_patch[, rev(seq_len(sh)), , drop = FALSE]

        tx2 <- min(src_px2, gx + sw - 1L)
        ty2 <- min(src_py2, gy + sh - 1L)
        cur_sw <- tx2 - gx + 1L
        cur_sh <- ty2 - gy + 1L

        if (cur_sw > 0L && cur_sh > 0L) {
          alpha <- runif(1, 0.70, 0.95)
          old_target <- new_arr[gx:tx2, gy:ty2, , drop = FALSE]
          new_target <- slice_patch[1:cur_sw, 1:cur_sh, , drop = FALSE]
          new_arr[gx:tx2, gy:ty2, ] <- alpha * new_target + (1.0 - alpha) * old_target
        }
      }
    }
  } else if (fill == "mean" || length(border_rects) == 0L) {
    for (ch in seq_len(n_channels)) {
      m_ch <- mean(arr[, , ch])
      new_arr[src_px1:src_px2, src_py1:src_py2, ch] <- m_ch + rnorm((src_px2 - src_px1 + 1L) * (src_py2 - src_py1 + 1L), 0, 0.03)
    }
  } else if (fill == "noise") {
    for (ch in seq_len(n_channels)) {
      new_arr[src_px1:src_px2, src_py1:src_py2, ch] <- runif((src_px2 - src_px1 + 1L) * (src_py2 - src_py1 + 1L))
    }
  } else if (fill == "blur") {
    for (ch in seq_len(n_channels)) {
      m_ch <- mean(arr[, , ch])
      new_arr[src_px1:src_px2, src_py1:src_py2, ch] <- m_ch
    }
  }

  if (blend_edge > 0L && blend_edge < min(patch_w, patch_h) %/% 4L) {
    bx_seq <- seq_len(blend_edge)
    for (b in bx_seq) {
      wt <- b / (blend_edge + 1.0)
      if (src_px1 + b - 1L <= src_px2) {
        new_arr[src_px1 + b - 1L, src_py1:src_py2, ] <- wt * new_arr[src_px1 + b - 1L, src_py1:src_py2, ] + (1.0 - wt) * arr[max(1L, src_px1 - 1L), src_py1:src_py2, ]
      }
      if (src_px2 - b + 1L >= src_px1) {
        new_arr[src_px2 - b + 1L, src_py1:src_py2, ] <- wt * new_arr[src_px2 - b + 1L, src_py1:src_py2, ] + (1.0 - wt) * arr[min(W, src_px2 + 1L), src_py1:src_py2, ]
      }
      if (src_py1 + b - 1L <= src_py2) {
        new_arr[src_px1:src_px2, src_py1 + b - 1L, ] <- wt * new_arr[src_px1:src_px2, src_py1 + b - 1L, ] + (1.0 - wt) * arr[src_px1:src_px2, max(1L, src_py1 - 1L), ]
      }
      if (src_py2 - b + 1L >= src_py1) {
        new_arr[src_px1:src_px2, src_py2 - b + 1L, ] <- wt * new_arr[src_px1:src_px2, src_py2 - b + 1L, ] + (1.0 - wt) * arr[src_px1:src_px2, min(H, src_py2 + 1L), ]
      }
    }
  }

  new_arr[dst_px1:dst_px2, dst_py1:dst_py2, ] <- patch

  new_lines <- character(0)
  for (item in parsed_items) {
    if (item$type == "box") {
      n_xc <- item$xc + actual_shift_x
      n_yc <- item$yc + actual_shift_y
      bx1 <- max(0.0, n_xc - item$w * 0.5)
      bx2 <- min(1.0, n_xc + item$w * 0.5)
      by1 <- max(0.0, n_yc - item$h * 0.5)
      by2 <- min(1.0, n_yc + item$h * 0.5)
      n_w <- bx2 - bx1
      n_h <- by2 - by1
      f_xc <- (bx1 + bx2) * 0.5
      f_yc <- (by1 + by2) * 0.5
      if (n_w >= 0.005 && n_h >= 0.005) {
        new_lines <- c(new_lines, sprintf("%d %.6f %.6f %.6f %.6f", item$cid, f_xc, f_yc, n_w, n_h))
      }
    } else if (item$type == "pose") {
      n_xc <- item$xc + actual_shift_x
      n_yc <- item$yc + actual_shift_y
      bx1 <- max(0.0, n_xc - item$w * 0.5)
      bx2 <- min(1.0, n_xc + item$w * 0.5)
      by1 <- max(0.0, n_yc - item$h * 0.5)
      by2 <- min(1.0, n_yc + item$h * 0.5)
      n_w <- bx2 - bx1
      n_h <- by2 - by1
      f_xc <- (bx1 + bx2) * 0.5
      f_yc <- (by1 + by2) * 0.5

      k_raw <- item$kpts
      nk <- length(k_raw) / 3L
      k_out <- character(nk * 3L)
      for (k in seq_len(nk)) {
        kx <- k_raw[(k - 1L) * 3L + 1L]
        ky <- k_raw[(k - 1L) * 3L + 2L]
        kv <- as.integer(k_raw[(k - 1L) * 3L + 3L])
        if (kv > 0L) {
          kx <- kx + actual_shift_x
          ky <- ky + actual_shift_y
        }
        k_out[(k - 1L) * 3L + 1L] <- sprintf("%.6f", max(0.0, min(1.0, kx)))
        k_out[(k - 1L) * 3L + 2L] <- sprintf("%.6f", max(0.0, min(1.0, ky)))
        k_out[(k - 1L) * 3L + 3L] <- as.character(kv)
      }
      if (n_w >= 0.005 && n_h >= 0.005) {
        new_lines <- c(new_lines, paste(c(item$cid, sprintf("%.6f", f_xc), sprintf("%.6f", f_yc), sprintf("%.6f", n_w), sprintf("%.6f", n_h), k_out), collapse = " "))
      }
    } else {
      n_xs <- item$xs + actual_shift_x
      n_ys <- item$ys + actual_shift_y
      c_xs <- pmax(0.0, pmin(1.0, n_xs))
      c_ys <- pmax(0.0, pmin(1.0, n_ys))
      pw <- max(c_xs) - min(c_xs)
      ph <- max(c_ys) - min(c_ys)
      if (pw >= 0.005 && ph >= 0.005) {
        coords_str <- sprintf("%.6f %.6f", c_xs, c_ys)
        new_lines <- c(new_lines, paste(c(item$cid, coords_str), collapse = " "))
      }
    }
  }

  out_im <- as_image(new_arr)
  list(image = out_im, labels = new_lines, shifted = TRUE, pos = chosen_pos)
}


# ==============================================================================
# 4. OFFLINE DATA AUGMENTATION (yolo_dataset_augment)
# ==============================================================================
#' @title Offline Data Augmentation with Geometric Bounding Box & Polygon Transformations
#' @name yolo_dataset_augment
#' @aliases yolo_dataset_aug
#' @description
#' `yolo_dataset_augment()` applies robust offline geometric and photometric augmentations
#' (horizontal/vertical flips, 90/180/270 rotations, brightness, contrast, noise, and blur)
#' to a YOLO dataset. Crucially, all bounding boxes and segmentation polygon vertices
#' are dynamically recomputed with mathematical precision to match the transformed images.
#'
#' Augmentations are exclusively applied to the `train/` split, strictly keeping `val/`
#' and `test/` pristine for unbiased evaluation.
#'
#' @details
#' \strong{How `copies`, `total_imgs`, and transformations work:}
#'
#' \itemize{
#'   \item \code{total_imgs}: When specified, sets the \strong{exact total number of training images}
#'     desired after augmentation. For example, if you currently have 500 original training images
#'     and specify \code{total_imgs = 1250}, the function generates exactly 750 new augmented variants
#'     (distributed evenly across training images and classes). Overrides \code{copies}.
#'   \item \code{copies = 1L} (default when \code{total_imgs = NULL}): Generates \strong{1 new augmented variant per original training image} (doubling training size).
#'   \item \code{copies = 3L}: Generates \strong{3 distinct augmented variants per original training image} (quadrupling training size).
#' }
#'
#' The parameters (\code{rotate}, \code{flip_h}, \code{flip_v}, \code{brightness},
#' \code{contrast}, \code{noise}, \code{blur}) define the \emph{pool of eligible transformations}.
#' Rather than generating a full Cartesian grid product of all possible combinations, the
#' function samples a unique combination of transformations for each copy requested.
#'
#' Already augmented images (named with suffix \code{_aug[0-9]+}) are automatically skipped,
#' preventing recursive dataset explosion if the function is executed multiple times.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param out_dir Output directory. If `NULL` (default), augments the dataset in-place.
#' @param total_imgs Integer. The total number of training images desired after data augmentation.
#'   For example, if you have 500 images and request `total_imgs = 1250`, exactly 750 new augmented
#'   variants will be generated and distributed across images. If `NULL` (default), `copies` is used.
#' @param copies Integer. The number of augmented variants generated per original training image (default: `1L`).
#'   Ignored if `total_imgs` is specified.
#' @param flip_h Logical. Allow horizontal flipping with probability 0.5 (default: `FALSE`).
#' @param flip_v Logical. Allow vertical flipping with probability 0.5 (default: `FALSE`).
#' @param rotate Vector or logical. Rotation augmentation in degrees. Can be:
#'   * A single numeric angle (e.g., `rotate = 45` or `rotate = 90` to strictly rotate by that angle),
#'   * A continuous range `c(min, max)` of angles (default: `c(-45, 45)`),
#'   * A vector of discrete angles to sample from (e.g. `c(0, 90, 180, 270)`),
#'   * `TRUE` (samples uniformly from `rotate_range`), or `FALSE`/`0` (no rotation).
#' @param rotate_range Numeric vector of length 2 giving `c(min, max)` rotation angles in degrees
#'   when `rotate = TRUE` (default: `c(-45, 45)`).
#' @param perspective Logical. Whether to apply 3D camera perspective warp augmentation (default: `FALSE`).
#'   When `TRUE`, simulates realistic pitch, yaw, and roll camera tilts.
#' @param perspective_x Numeric vector of length 2 giving `c(min, max)` pitch tilt range around X in degrees (default: `c(20, 40)`).
#'   If given as a positive range (e.g. `c(15, 30)`), samples random positive or negative tilts between 20 and 40 degrees.
#' @param perspective_y Numeric vector of length 2 giving `c(min, max)` yaw tilt range around Y in degrees (default: `c(20, 40)`).
#'   If given as a positive range (e.g. `c(15, 30)`), samples random positive or negative tilts between 20 and 40 degrees.
#' @param perspective_z Numeric vector of length 2 giving `c(min, max)` roll angle range around Z in degrees (default: `c(20, 40)`).
#'   If given as a positive range (e.g. `c(15, 30)`), samples random positive or negative rolls between 20 and 40 degrees.
#' @param perspective_scale Numeric vector of length 2 giving `c(min, max)` zoom scale range (default: `c(0.95, 1.05)`).
#' @param bg_col Background fill color for exposed canvas (default: `"white"`).
#' @param brightness Logical. Allow random brightness jitter (default: `FALSE`).
#' @param brightness_range Numeric vector of length 2 giving `c(min, max)` brightness scaling factors (default: `c(0.6, 1.2)`).
#' @param contrast Logical. Allow random contrast adjustment via adaptive histogram equalization ([image_contrast()]). Defaults to `FALSE`.
#' @param contrast_range Numeric vector of length 2 giving `c(min, max)` clip limit for [image_contrast()] (default: `c(1, 3)`).
#' @param noise Logical. Allow additive Gaussian noise (default: `FALSE`).
#' @param noise_sd Numeric vector of length 2 giving `c(min, max)` standard deviation of Gaussian noise (default: `c(0.05, 0.1)`).
#' @param blur Logical. Allow Gaussian blurring (default: `FALSE`).
#' @param blur_sigma Numeric vector of length 2 giving `c(min, max)` Gaussian blur sigma (default: `c(2, 4)`).
#' @param cutout Logical. Whether to apply Cutout / Random Erasing augmentation to simulate occlusions (default: `FALSE`).
#' @param cutout_target Target location for cutout rectangles. Either `"object"` (default, occludes parts of annotated objects) or `"image"` (places cutout rectangles randomly across the entire image canvas).
#' @param cutout_pos Positioning strategy when `cutout_target = "object"`. Either `"corner"` (default, randomly anchors the cutout at one of the four bounding box corners, simulating realistic partial occlusion of objects) or `"random"` (places cutout anywhere within the bounding box).
#' @param cutout_n Number of cutout rectangles to generate per image/object. Can be an integer or integer range `c(min, max)` (default: `c(3, 5)`).
#' @param cutout_size Numeric vector of length 2 giving `c(min, max)` relative dimensions of the cutout rectangle. When `cutout_target = "object"`, it is relative to object bounding box dimensions (default: `c(0.4, 0.6)`). When `cutout_target = "image"`, it is relative to image dimensions.
#' @param cutout_fill Fill color or pattern for cutout rectangles. Can be `"mean"` (channel-wise mean image color, default), `"black"`, `"white"`, `"noise"` (random uniform noise), or any valid R color name.
#' @param mosaic Logical or numeric. Whether to generate 4-image mosaic augmentations ([image_mosaic()]).
#'   Can be `TRUE` (applies mosaic with probability `mosaic_prob`), `FALSE` (default, disabled),
#'   or a numeric probability between 0 and 1 (e.g. `mosaic = 1.0` for 100% mosaic augmentation).
#' @param mosaic_prob Numeric between 0 and 1. Probability of applying mosaic augmentation when `mosaic = TRUE`.
#'   Defaults to `0.5`.
#' @param mosaic_overlap Fraction of overlap between adjacent tiles in the mosaic composite
#'   (default `0.20`, i.e., 20%). Passed to [image_mosaic()].
#' @param mosaic_crop Region to crop from each image for the mosaic: `"quadrant"` (default),
#'   `"random"`, or `"center"`.
#' @param mosaic_blend Blending mode for overlapping regions in the mosaic: `"overlay"` (default
#'   for object detection datasets), `"feather"`, or `"none"`.
#' @param mosaic_min_area Minimum visible area fraction (between 0 and 1) required
#'   to retain an object cut by mosaic quadrant seams (default `0.10`, i.e., 10 percent).
#' @param scramble_background Logical or numeric. Whether to synthesize negative background
#'   images via jigsaw patch scrambling ([yolo_dataset_scramble()]). If numeric (e.g. `0.05`),
#'   it specifies the proportion of negative background images to generate relative to the training set.
#'   If `TRUE`, defaults to `scramble_prop = 0.05` (5% negative background samples). Defaults to `FALSE`.
#' @param scramble_prop Numeric between 0 and 1. Proportion of negative background images to generate
#'   when `scramble_background = TRUE` (default: `0.05`).
#' @param scramble_grid Integer vector of length 1 or 2 specifying grid divisions for tile scrambling
#'   (default: `c(8, 8)`). Ignored when `scramble_mode = "pixel"`.
#' @param scramble_mode Character specifying the scrambling strategy: `"grid"` (default jigsaw tile shuffle),
#'   `"block_swap"` (local block swapping preserving partial background continuity), or `"pixel"` (global pixel shuffle).
#' @param permute Logical. Whether to apply object cluster permutation/relocation augmentation to break center bias. Defaults to `FALSE`.
#' @param permute_prob Numeric probability in `[0, 1]` of applying cluster permutation to an image when `permute = TRUE`. Defaults to `0.5`.
#' @param permute_pos Character string specifying target placement regions for relocated clusters: `"random"`, `"left"`, `"right"`, `"top"`, `"bottom"`, or `"corner"`. Defaults to `"random"`.
#' @param permute_fill Character string specifying how vacated object regions are filled: `"border_slices"`, `"mean"`, `"blur"`, or `"noise"`. Defaults to `"border_slices"`.
#' @param permute_margin Numeric relative margin fraction kept clear around image borders when relocating clusters. Defaults to `0.02`.
#' @param quality Integer between 1 and 100. JPEG compression quality for exported augmented images (default: `90`). Quality 85-90 preserves pristine pixel fidelity while encoding several times faster than 95+.
#' @param max_ops Maximum number of transformations/alterations applied to each augmented image copy (default: `c(1, 3)`). Can be an integer (e.g., `2` or `3`, meaning between 1 and `max_ops` operations are randomly selected), a range `c(min, max)` (e.g. `c(1, 3)`), or `NULL` (no limit, independent probabilities for all enabled operations).
#' @param seed Optional integer seed for reproducibility (default: 42).
#' @param verbose Logical. Whether to print informative messages. Defaults to `TRUE`.
#' @param progress Logical. Whether to display a CLI progress bar during data augmentation. Defaults to `TRUE`.
#' @param as_job Logical. If `TRUE` and running inside RStudio, launches the data augmentation
#'   process in the background via RStudio Background Jobs, leaving the main R console completely free.
#'   Defaults to `FALSE`.
#' @param ... Additional arguments (supports alias `blurr = TRUE`, `blurr_sigma`, `noise_range`).
#'
#' @return A `yolo_dataset` object pointing to the augmented dataset directory.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Default: 1 augmented copy per image (doubles training size)
#' ds_aug1 <- yolo_dataset_augment("dataset_cafe")
#'
#' # 2. Set an exact target total number of images (e.g. exactly 1250 training images)
#' ds_aug1250 <- yolo_dataset_augment("dataset_cafe", total_imgs = 1250)
#'
#' # 3. Cutout augmentation focused on objects
#' ds_cutout <- yolo_dataset_augment("dataset_cafe",
#'                                   cutout = TRUE,
#'                                   cutout_target = "object",
#'                                   cutout_size = c(0.4, 0.6))
#'
#' # 4. Mosaic augmentation (4 images stitched with 20% overlap)
#' ds_mosaic <- yolo_dataset_augment("dataset_cafe",
#'                                   mosaic = TRUE,
#'                                   mosaic_overlap = 0.20)
#'
#' # 5. Break center bias: permute/relocate target clusters to borders with authentic background fill
#' ds_perm <- yolo_dataset_augment("dataset_cafe",
#'                                 permute = TRUE,
#'                                 permute_prob = 0.5,
#'                                 permute_fill = "border_slices")
yolo_dataset_augment <- function(data,
                                 out_dir = NULL,
                                 total_imgs = NULL,
                                 copies = 1L,
                                 flip_h = FALSE,
                                 flip_v = FALSE,
                                 rotate = TRUE,
                                 rotate_range = c(-45, 45),
                                 perspective = FALSE,
                                 perspective_x = c(15, 30),
                                 perspective_y = c(15, 30),
                                 perspective_z = c(15, 30),
                                 perspective_scale = c(0.95, 1.05),
                                 bg_col = "white",
                                 brightness = FALSE,
                                 brightness_range = c(0.6, 1.2),
                                 contrast = FALSE,
                                 contrast_range = c(1, 3),
                                 noise = FALSE,
                                 noise_sd = c(0.05, 0.1),
                                 blur = FALSE,
                                 blur_sigma = c(2, 4),
                                 cutout = FALSE,
                                 cutout_target = c("object", "image"),
                                 cutout_pos = c("corner", "random"),
                                 cutout_n = c(3, 5),
                                 cutout_size = c(0.4, 0.6),
                                 cutout_fill = "mean",
                                 mosaic = FALSE,
                                 mosaic_prob = 0.5,
                                 mosaic_overlap = 0.20,
                                 mosaic_crop = c("quadrant", "random", "center"),
                                 mosaic_blend = c("none", "feather", "overlay"),
                                 mosaic_min_area = 0.10,
                                 scramble_background = FALSE,
                                 scramble_prop = 0.05,
                                 scramble_grid = c(8, 8),
                                 scramble_mode = c("grid", "block_swap", "pixel"),
                                 permute = FALSE,
                                 permute_prob = 0.5,
                                 permute_pos = c("random", "left", "right", "top", "bottom", "corner"),
                                 permute_fill = c("border_slices", "mean", "blur", "noise"),
                                 permute_margin = 0.02,
                                 quality = 90,
                                 max_ops = c(1, 3),
                                 seed = 42,
                                 verbose = TRUE,
                                 progress = TRUE,
                                 as_job = FALSE,
                                 ...) {
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with data augmentation in current session.")
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
      job_script_path <- file.path(job_dir, paste0("yolo_aug_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_aug_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        d_p <- .yolo_get_dir(data)
        basename(d_p)
      }, error = function(e) {
        if (inherits(data, "yolo_dataset")) basename(data$dir)
        else if (is.character(data) && nzchar(data[1])) basename(data[1])
        else "dataset"
      })
      j_name <- paste0("YOLO Augment: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_augment()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Data Augmentation - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_augment, job_info$args)\n",
        "target_dir <- if (!is.null(res$dir)) res$dir else job_info$args$out_dir\n",
        "if (!is.null(target_dir) && dir.exists(target_dir)) {\n",
        "  save_path <- file.path(target_dir, 'dataset_augmented.rds')\n",
        "  saveRDS(res, save_path)\n",
        "  cat('\\n[OK] Augmented dataset object successfully saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Augmentation Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Data Augmentation (pliman)")
      cli::cli_alert_success("Data augmentation has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target dataset: {.val {target_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track live augmentation progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  cutout_target <- match.arg(cutout_target)
  cutout_pos <- match.arg(cutout_pos)
  permute_pos <- if (length(permute_pos) > 1L) permute_pos[1] else permute_pos
  permute_fill <- if (length(permute_fill) > 1L) permute_fill[1] else permute_fill
  mosaic_crop <- match.arg(mosaic_crop)
  mosaic_blend <- match.arg(mosaic_blend)
  dots <- list(...)
  if ("jpeg_quality" %in% names(dots)) quality <- as.integer(dots$jpeg_quality)
  quality <- max(1L, min(100L, as.integer(quality)))
  if ("cutout_position" %in% names(dots)) cutout_pos <- dots$cutout_position
  if ("blurr" %in% names(dots)) blur <- isTRUE(dots$blurr)
  if ("blurr_sigma" %in% names(dots)) blur_sigma <- dots$blurr_sigma
  if ("noise_range" %in% names(dots)) noise_sd <- dots$noise_range
  if ("contrast_clip_limit" %in% names(dots)) contrast_range <- dots$contrast_clip_limit
  if ("max_alterations" %in% names(dots)) max_ops <- dots$max_alterations
  if ("max_transforms" %in% names(dots)) max_ops <- dots$max_transforms
  if ("mosaic" %in% names(dots)) mosaic <- dots$mosaic
  if ("mosaic_probability" %in% names(dots)) mosaic_prob <- dots$mosaic_probability
  if ("mosaic_prob" %in% names(dots)) mosaic_prob <- dots$mosaic_prob

  mosaic_p <- if (is.numeric(mosaic)) max(0.0, min(1.0, mosaic[1L])) else if (isTRUE(mosaic)) max(0.0, min(1.0, mosaic_prob[1L])) else 0.0
  has_mosaic <- (mosaic_p > 0.0)

  .wrap_aug <- function(tdir, tsk, cls, n_created = 0L) {
    w <- .yolo_wrap_dataset(tdir, task = tsk, classes = cls)
    class(w) <- unique(c("yolo_dataset_aug", class(w)))
    w$is_augmented <- TRUE
    w$total_augmented <- as.integer(n_created)
    w
  }

  .sample_one <- function(x) {
    if (length(x) == 0L) return(0)
    if (length(x) == 1L) return(x[1])
    sample(x, 1)
  }

  .get_rotate_sampler <- function(rot, r_range = c(-45, 45)) {
    if (is.null(rot) || isFALSE(rot)) {
      return(function() 0)
    }
    if (isTRUE(rot)) {
      return(function() round(runif(1, min(r_range), max(r_range))))
    }
    if (is.numeric(rot)) {
      if (length(rot) == 1L) {
        val <- as.numeric(rot)
        return(function() val)
      }
      if (length(rot) == 2L && !all(rot %in% c(0, 90, 180, 270))) {
        return(function() round(runif(1, min(rot), max(rot))))
      }
      return(function() .sample_one(rot))
    }
    function() 0
  }

  .sample_perspective_val <- function(p) {
    if (is.null(p)) return(0.0)
    p <- as.numeric(p)
    if (length(p) == 0L || all(is.na(p))) return(0.0)
    if (length(p) == 1L) {
      if (p == 0) return(0.0)
      if (p > 0) return(sample(c(-1, 1), 1L) * p)
      return(p)
    }
    p_min <- min(p, na.rm = TRUE)
    p_max <- max(p, na.rm = TRUE)
    if (p_min == 0 && p_max == 0) return(0.0)
    if (p_min >= 0) {
      mag <- if (p_min == p_max) p_min else runif(1, p_min, p_max)
      return(sample(c(-1, 1), 1L) * mag)
    }
    runif(1, p_min, p_max)
  }

  .get_perspective_sampler <- function(persp, px, py, pz, pscale) {
    if (is.null(persp) || isFALSE(persp)) {
      return(NULL)
    }
    function() {
      list(
        px = .sample_perspective_val(px),
        py = .sample_perspective_val(py),
        pz = .sample_perspective_val(pz),
        pscale = if (length(pscale) == 1L) as.numeric(pscale) else runif(1, min(pscale), max(pscale))
      )
    }
  }

  .yolo_flip_annotations <- function(raw_lines, do_h = FALSE, do_v = FALSE, flip_idx = NULL) {
    if (length(raw_lines) == 0L || (!do_h && !do_v)) return(raw_lines)
    out_lines <- character(0)
    for (ln in raw_lines) {
      vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
      if (length(vals) < 5L) next
      cid <- as.integer(vals[1])
      is_pose <- (length(vals) > 5L && (length(vals) - 5L) %% 3L == 0L)

      if (is_pose && length(vals) >= 8L) {
        xc <- vals[2]; yc <- vals[3]; w <- vals[4]; h <- vals[5]
        if (do_h) xc <- 1.0 - xc
        if (do_v) yc <- 1.0 - yc

        k_raw <- vals[6:length(vals)]
        nk <- length(k_raw) / 3L
        k_out <- character(nk * 3L)
        for (k in seq_len(nk)) {
          kx <- k_raw[(k - 1L) * 3L + 1L]
          ky <- k_raw[(k - 1L) * 3L + 2L]
          kv <- as.integer(k_raw[(k - 1L) * 3L + 3L])
          if (kv > 0L) {
            if (do_h) kx <- 1.0 - kx
            if (do_v) ky <- 1.0 - ky
          }
          k_out[(k - 1L) * 3L + 1L] <- sprintf("%.6f", max(0.0, min(1.0, kx)))
          k_out[(k - 1L) * 3L + 2L] <- sprintf("%.6f", max(0.0, min(1.0, ky)))
          k_out[(k - 1L) * 3L + 3L] <- as.character(kv)
        }
        if (do_h && !is.null(flip_idx) && length(flip_idx) == nk) {
          reordered <- character(nk * 3L)
          for (k in seq_len(nk)) {
            t_k <- flip_idx[k]
            if (t_k >= 1L && t_k <= nk) {
              reordered[(k - 1L) * 3L + 1L] <- k_out[(t_k - 1L) * 3L + 1L]
              reordered[(k - 1L) * 3L + 2L] <- k_out[(t_k - 1L) * 3L + 2L]
              reordered[(k - 1L) * 3L + 3L] <- k_out[(t_k - 1L) * 3L + 3L]
            } else {
              reordered[(k - 1L) * 3L + 1L] <- k_out[(k - 1L) * 3L + 1L]
              reordered[(k - 1L) * 3L + 2L] <- k_out[(k - 1L) * 3L + 2L]
              reordered[(k - 1L) * 3L + 3L] <- k_out[(k - 1L) * 3L + 3L]
            }
          }
          k_out <- reordered
        }
        out_lines <- c(out_lines, paste(c(cid, sprintf("%.6f", xc), sprintf("%.6f", yc), sprintf("%.6f", w), sprintf("%.6f", h), k_out), collapse = " "))
      } else if (length(vals) == 5L) {
        xc <- vals[2]; yc <- vals[3]; w <- vals[4]; h <- vals[5]
        if (do_h) xc <- 1.0 - xc
        if (do_v) yc <- 1.0 - yc
        out_lines <- c(out_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, w, h))
      } else {
        coords <- vals[-1]
        xs <- coords[seq(1, length(coords), by = 2)]
        ys <- coords[seq(2, length(coords), by = 2)]
        if (do_h) xs <- 1.0 - xs
        if (do_v) ys <- 1.0 - ys
        coord_pairs <- sprintf("%.6f %.6f", xs, ys)
        out_lines <- c(out_lines, paste(c(cid, coord_pairs), collapse = " "))
      }
    }
    out_lines
  }

  .build_rot_matrix <- function(W_in, H_in, W_out, H_out, ang) {
    rad <- -ang * pi / 180.0
    cosA <- cos(rad); sinA <- sin(rad)
    cx_src <- (W_in - 1) * 0.5
    cy_src <- (H_in - 1) * 0.5
    cx_dst <- (W_out - 1) * 0.5
    cy_dst <- (H_out - 1) * 0.5

    Tin <- matrix(c(1, 0, -cx_src, 0, 1, -cy_src, 0, 0, 1), 3, 3, byrow = TRUE)
    Rmat <- matrix(c(cosA, sinA, 0, -sinA, cosA, 0, 0, 0, 1), 3, 3, byrow = TRUE)
    Tout <- matrix(c(1, 0, cx_dst, 0, 1, cy_dst, 0, 0, 1), 3, 3, byrow = TRUE)
    Tout %*% Rmat %*% Tin
  }

  .angles_ellipse <- seq(0, 2 * pi, length.out = 33L)[-33L]
  .cos_ellipse <- cos(.angles_ellipse)
  .sin_ellipse <- sin(.angles_ellipse)

  .yolo_transform_annotations <- function(raw_lines, W_src, H_src, W_dst, H_dst, M) {
    if (length(raw_lines) == 0L) return(character(0))
    aug_lines <- character(0)

    for (ln in raw_lines) {
      vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
      if (length(vals) < 5L) next
      cid <- as.integer(vals[1])
      is_pose <- (length(vals) > 5L && (length(vals) - 5L) %% 3L == 0L)

      if (is_pose && length(vals) >= 8L) {
        # 1. Transform Bounding Box using ellipse projection
        xc <- vals[2]; yc <- vals[3]; w <- vals[4]; h <- vals[5]
        cx <- xc * W_src; cy <- yc * H_src
        rx <- (w * W_src) * 0.5; ry <- (h * H_src) * 0.5

        pts <- rbind(
          cx + rx * .cos_ellipse,
          cy + ry * .sin_ellipse,
          1.0
        )

        trans <- M %*% pts
        tw <- trans[3, ]; tw[abs(tw) < 1e-9] <- 1e-9
        valid <- tw > 0.01
        if (!any(valid)) next
        tx <- trans[1, valid] / tw[valid]
        ty <- trans[2, valid] / tw[valid]

        c_xmin <- max(0.0, min(as.numeric(W_dst), min(tx)))
        c_xmax <- max(0.0, min(as.numeric(W_dst), max(tx)))
        c_ymin <- max(0.0, min(as.numeric(H_dst), min(ty)))
        c_ymax <- max(0.0, min(as.numeric(H_dst), max(ty)))

        if ((c_xmax - c_xmin) >= 2.0 && (c_ymax - c_ymin) >= 2.0) {
          n_xc <- (c_xmin + c_xmax) / (2.0 * W_dst)
          n_yc <- (c_ymin + c_ymax) / (2.0 * H_dst)
          n_w  <- (c_xmax - c_xmin) / W_dst
          n_h  <- (c_ymax - c_ymin) / H_dst

          n_xc <- max(0.0, min(1.0, n_xc))
          n_yc <- max(0.0, min(1.0, n_yc))
          n_w  <- max(0.0, min(1.0, n_w))
          n_h  <- max(0.0, min(1.0, n_h))

          # 2. Transform Keypoints using affine matrix M
          k_raw <- vals[6:length(vals)]
          nk <- length(k_raw) / 3L
          k_out <- character(nk * 3L)

          for (k in seq_len(nk)) {
            kx <- k_raw[(k - 1L) * 3L + 1L] * W_src
            ky <- k_raw[(k - 1L) * 3L + 2L] * H_src
            kv <- as.integer(k_raw[(k - 1L) * 3L + 3L])

            if (kv > 0L) {
              pt_k <- c(kx, ky, 1.0)
              trans_k <- M %*% pt_k
              tw_k <- if (abs(trans_k[3]) < 1e-9) 1e-9 else trans_k[3]
              tx_k <- trans_k[1] / tw_k
              ty_k <- trans_k[2] / tw_k

              if (tx_k >= 0 && tx_k <= W_dst && ty_k >= 0 && ty_k <= H_dst) {
                norm_kx <- max(0.0, min(1.0, tx_k / W_dst))
                norm_ky <- max(0.0, min(1.0, ty_k / H_dst))
              } else {
                norm_kx <- 0.0; norm_ky <- 0.0; kv <- 0L # Out of bounds -> absent
              }
            } else {
              norm_kx <- 0.0; norm_ky <- 0.0; kv <- 0L
            }

            k_out[(k - 1L) * 3L + 1L] <- sprintf("%.6f", norm_kx)
            k_out[(k - 1L) * 3L + 2L] <- sprintf("%.6f", norm_ky)
            k_out[(k - 1L) * 3L + 3L] <- as.character(kv)
          }

          aug_lines <- c(aug_lines, paste(c(cid, sprintf("%.6f", n_xc), sprintf("%.6f", n_yc), sprintf("%.6f", n_w), sprintf("%.6f", n_h), k_out), collapse = " "))
        }
      } else if (length(vals) == 5L) {
        # Box: xc, yc, w, h
        # Inscribed ellipse projection (tight bounding box avoiding corner inflation)
        xc <- vals[2]; yc <- vals[3]; w <- vals[4]; h <- vals[5]
        cx <- xc * W_src
        cy <- yc * H_src
        rx <- (w * W_src) * 0.5
        ry <- (h * H_src) * 0.5

        pts <- rbind(
          cx + rx * .cos_ellipse,
          cy + ry * .sin_ellipse,
          1.0
        )

        trans <- M %*% pts
        tw <- trans[3, ]
        tw[abs(tw) < 1e-9] <- 1e-9
        valid <- tw > 0.01
        if (!any(valid)) next
        tx <- trans[1, valid] / tw[valid]
        ty <- trans[2, valid] / tw[valid]

        new_xmin <- min(tx)
        new_xmax <- max(tx)
        new_ymin <- min(ty)
        new_ymax <- max(ty)

        c_xmin <- max(0.0, min(as.numeric(W_dst), new_xmin))
        c_xmax <- max(0.0, min(as.numeric(W_dst), new_xmax))
        c_ymin <- max(0.0, min(as.numeric(H_dst), new_ymin))
        c_ymax <- max(0.0, min(as.numeric(H_dst), new_ymax))

        if ((c_xmax - c_xmin) >= 2.0 && (c_ymax - c_ymin) >= 2.0) {
          n_xc <- (c_xmin + c_xmax) / (2.0 * W_dst)
          n_yc <- (c_ymin + c_ymax) / (2.0 * H_dst)
          n_w  <- (c_xmax - c_xmin) / W_dst
          n_h  <- (c_ymax - c_ymin) / H_dst

          n_xc <- max(0.0, min(1.0, n_xc))
          n_yc <- max(0.0, min(1.0, n_yc))
          n_w  <- max(0.0, min(1.0, n_w))
          n_h  <- max(0.0, min(1.0, n_h))

          aug_lines <- c(aug_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, n_xc, n_yc, n_w, n_h))
        }
      } else {
        # Polygon: class x1 y1 x2 y2 ...
        coords <- vals[-1]
        xs <- coords[seq(1, length(coords), by = 2)] * W_src
        ys <- coords[seq(2, length(coords), by = 2)] * H_src
        n_pts <- length(xs)

        pts <- rbind(xs, ys, rep(1, n_pts))
        trans <- M %*% pts
        tw <- trans[3, ]
        tw[abs(tw) < 1e-9] <- 1e-9
        tx <- trans[1, ] / tw
        ty <- trans[2, ] / tw

        norm_x <- pmax(0.0, pmin(1.0, tx / W_dst))
        norm_y <- pmax(0.0, pmin(1.0, ty / H_dst))

        interleaved <- as.vector(rbind(norm_x, norm_y))
        aug_lines <- c(aug_lines, paste(c(cid, sprintf("%.6f", interleaved)), collapse = " "))
      }
    }

    aug_lines
  }

  .apply_cutout <- function(arr,
                            raw_lines = character(0),
                            target = "object",
                            n_cutouts = c(1, 3),
                            size_range = c(0.15, 0.40),
                            fill = "black",
                            pos = "corner") {
    W <- dim(arr)[1]
    H <- dim(arr)[2]
    n_channels <- if (length(dim(arr)) >= 3L) dim(arr)[3] else 1L

    fill_is_noise <- identical(fill, "noise")
    fill_val <- if (fill_is_noise) {
      NULL
    } else if (identical(fill, "mean")) {
      if (n_channels == 3L) {
        c(mean(arr[, , 1]), mean(arr[, , 2]), mean(arr[, , 3]))
      } else {
        mean(arr)
      }
    } else {
      tryCatch({
        rgb_vec <- grDevices::col2rgb(fill)[, 1] / 255.0
        if (n_channels == 3L) rgb_vec else mean(rgb_vec)
      }, error = function(e) {
        if (n_channels == 3L) c(0, 0, 0) else 0
      })
    }

    fill_patch <- function(arr_in, x1, y1, x2, y2) {
      x1 <- max(1L, min(W, as.integer(round(x1))))
      x2 <- max(1L, min(W, as.integer(round(x2))))
      y1 <- max(1L, min(H, as.integer(round(y1))))
      y2 <- max(1L, min(H, as.integer(round(y2))))
      if (x2 < x1 || y2 < y1) return(arr_in)

      pw <- (x2 - x1 + 1L)
      ph <- (y2 - y1 + 1L)

      if (fill_is_noise) {
        if (n_channels == 3L) {
          arr_in[x1:x2, y1:y2, 1] <- runif(pw * ph)
          arr_in[x1:x2, y1:y2, 2] <- runif(pw * ph)
          arr_in[x1:x2, y1:y2, 3] <- runif(pw * ph)
        } else {
          arr_in[x1:x2, y1:y2, 1] <- runif(pw * ph)
        }
      } else {
        if (n_channels == 3L) {
          arr_in[x1:x2, y1:y2, 1] <- fill_val[1]
          arr_in[x1:x2, y1:y2, 2] <- fill_val[2]
          arr_in[x1:x2, y1:y2, 3] <- fill_val[3]
        } else {
          arr_in[x1:x2, y1:y2, 1] <- fill_val[1]
        }
      }
      arr_in
    }

    parsed_boxes <- list()
    if (length(raw_lines) > 0L) {
      for (ln in raw_lines) {
        vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
        if (length(vals) < 5L || any(is.na(vals))) next
        if (length(vals) == 5L || (length(vals) > 5L && (length(vals) - 5L) %% 3L == 0L)) {
          xc <- vals[2]; yc <- vals[3]; bw <- vals[4]; bh <- vals[5]
          parsed_boxes[[length(parsed_boxes) + 1L]] <- c(
            (xc - bw * 0.5) * W,
            (xc + bw * 0.5) * W,
            (yc - bh * 0.5) * H,
            (yc + bh * 0.5) * H
          )
        } else if (length(vals) >= 6L) {
          coords <- vals[-1]
          px <- coords[seq(1, length(coords), by = 2L)] * W
          py <- coords[seq(2, length(coords), by = 2L)] * H
          parsed_boxes[[length(parsed_boxes) + 1L]] <- c(
            min(px), max(px), min(py), max(py)
          )
        }
      }
    }

    n_k <- if (length(n_cutouts) == 1L) {
      as.integer(n_cutouts)
    } else {
      sample(as.integer(n_cutouts[1]):as.integer(n_cutouts[length(n_cutouts)]), 1L)
    }
    n_k <- max(1L, n_k)

    s_min <- min(size_range)
    s_max <- max(size_range)

    if (target == "object" && length(parsed_boxes) > 0L) {
      chosen_indices <- sample(seq_along(parsed_boxes), size = min(n_k, length(parsed_boxes)), replace = (n_k > length(parsed_boxes)))
      for (b_idx in chosen_indices) {
        b <- parsed_boxes[[b_idx]]
        bx1 <- b[1]; bx2 <- b[2]; by1 <- b[3]; by2 <- b[4]
        bw <- max(4.0, bx2 - bx1)
        bh <- max(4.0, by2 - by1)

        s_fx <- if (s_min == s_max) s_min else runif(1, s_min, s_max)
        s_fy <- if (s_min == s_max) s_min else runif(1, s_min, s_max)
        cw <- max(2.0, min(bw * 0.85, bw * s_fx))
        ch <- max(2.0, min(bh * 0.85, bh * s_fy))

        if (identical(pos, "corner")) {
          corner <- sample(1L:4L, 1L)
          if (corner == 1L) {
            # Top-Left corner
            rx1 <- bx1; rx2 <- bx1 + cw
            ry1 <- by1; ry2 <- by1 + ch
          } else if (corner == 2L) {
            # Top-Right corner
            rx1 <- bx2 - cw; rx2 <- bx2
            ry1 <- by1; ry2 <- by1 + ch
          } else if (corner == 3L) {
            # Bottom-Left corner
            rx1 <- bx1; rx2 <- bx1 + cw
            ry1 <- by2 - ch; ry2 <- by2
          } else {
            # Bottom-Right corner
            rx1 <- bx2 - cw; rx2 <- bx2
            ry1 <- by2 - ch; ry2 <- by2
          }
        } else {
          rx1 <- runif(1, min = bx1, max = max(bx1, bx2 - cw))
          ry1 <- runif(1, min = by1, max = max(by1, by2 - ch))
          rx2 <- rx1 + cw
          ry2 <- ry1 + ch
        }

        arr <- fill_patch(arr, rx1, ry1, rx2, ry2)
      }
    } else {
      for (k in seq_len(n_k)) {
        s_f <- if (s_min == s_max) s_min else runif(1, s_min, s_max)
        cw <- max(4.0, W * s_f)
        ch <- max(4.0, H * s_f)

        rx1 <- runif(1, min = 1, max = max(1, W - cw))
        ry1 <- runif(1, min = 1, max = max(1, H - ch))
        rx2 <- rx1 + cw
        ry2 <- ry1 + ch

        arr <- fill_patch(arr, rx1, ry1, rx2, ry2)
      }
    }

    arr
  }

  sample_rot <- .get_rotate_sampler(rotate, rotate_range)
  persp_sampler <- .get_perspective_sampler(perspective, perspective_x, perspective_y, perspective_z, perspective_scale)
  is_fixed_rot <- is.numeric(rotate) && length(rotate) == 1L && rotate[1] != 0

  .sample_chosen_ops <- function() {
    eligible_ops <- character(0)
    if (isTRUE(flip_h)) eligible_ops <- c(eligible_ops, "flip_h")
    if (isTRUE(flip_v)) eligible_ops <- c(eligible_ops, "flip_v")
    if (!isFALSE(rotate) && !identical(rotate, 0)) eligible_ops <- c(eligible_ops, "rot")
    if (!is.null(persp_sampler)) eligible_ops <- c(eligible_ops, "persp")
    if (isTRUE(brightness)) eligible_ops <- c(eligible_ops, "bright")
    if (isTRUE(contrast)) eligible_ops <- c(eligible_ops, "contrast")
    if (isTRUE(noise)) eligible_ops <- c(eligible_ops, "noise")
    if (isTRUE(blur)) eligible_ops <- c(eligible_ops, "blur")
    if (isTRUE(cutout)) eligible_ops <- c(eligible_ops, "cutout")

    if (length(eligible_ops) == 0L) {
      return(list(
        do_fh = FALSE, do_fv = FALSE, do_rot = FALSE, ang = 0,
        do_persp = FALSE, persp_pars = NULL,
        do_bright = FALSE, do_contrast = FALSE, do_noise = FALSE,
        do_blur = FALSE, do_cutout = FALSE
      ))
    }

    if (!is.null(max_ops)) {
      k_min <- if (length(max_ops) >= 2L) max(1L, as.integer(min(max_ops))) else 1L
      k_max <- as.integer(max(max_ops))
      k_max <- max(k_min, k_max)
      k_min <- min(k_min, length(eligible_ops))
      k_max <- min(k_max, length(eligible_ops))
      n_ops <- if (k_min == k_max) k_min else sample(k_min:k_max, 1L)
      chosen <- sample(eligible_ops, size = n_ops, replace = FALSE)

      do_fh <- "flip_h" %in% chosen
      do_fv <- "flip_v" %in% chosen
      do_rot <- "rot" %in% chosen
      do_persp <- "persp" %in% chosen
      do_bright <- "bright" %in% chosen
      do_contrast <- "contrast" %in% chosen
      do_noise <- "noise" %in% chosen
      do_blur <- "blur" %in% chosen
      do_cutout <- "cutout" %in% chosen

      ang <- 0
      if (do_rot) {
        ang <- sample_rot()
        if (ang == 0) {
          if (is.numeric(rotate) && any(rotate != 0)) {
            non_zeros <- rotate[rotate != 0]
            ang <- .sample_one(non_zeros)
          } else {
            ang <- round(runif(1, min(rotate_range), max(rotate_range)))
            if (ang == 0) ang <- 90
          }
        }
      }

      persp_pars <- NULL
      if (do_persp && !is.null(persp_sampler)) {
        persp_pars <- persp_sampler()
        if (persp_pars$px == 0 && persp_pars$py == 0 && persp_pars$pz == 0) {
          val <- .sample_perspective_val(perspective_x)
          if (val == 0) val <- .sample_perspective_val(perspective_y)
          if (val == 0) val <- sample(c(-10, 10), 1L)
          persp_pars$px <- val
        }
      }
    } else {
      do_fh <- isTRUE(flip_h) && (runif(1) > 0.45)
      do_fv <- isTRUE(flip_v) && (runif(1) > 0.5)
      ang <- sample_rot()
      do_rot <- (ang != 0) && (is_fixed_rot || runif(1) > 0.35)
      persp_pars <- if (!is.null(persp_sampler)) persp_sampler() else NULL
      do_persp <- !is.null(persp_pars) && (runif(1) > 0.4)
      do_bright <- isTRUE(brightness) && (runif(1) > 0.35)
      do_contrast <- isTRUE(contrast) && (runif(1) > 0.4)
      do_noise <- isTRUE(noise) && (runif(1) > 0.45)
      do_blur <- isTRUE(blur) && (runif(1) > 0.45)
      do_cutout <- isTRUE(cutout) && (runif(1) > 0.4)

      # Ensure at least one eligible transformation is performed
      if (!do_fh && !do_fv && (!do_rot || ang == 0) && !do_persp && !do_bright && !do_contrast && !do_noise && !do_blur && !do_cutout) {
        forced <- .sample_one(eligible_ops)
        if (forced == "flip_h") do_fh <- TRUE
        else if (forced == "flip_v") do_fv <- TRUE
        else if (forced == "rot") {
          ang <- sample_rot()
          if (ang == 0) ang <- 90
          do_rot <- TRUE
        }
        else if (forced == "persp") {
          do_persp <- TRUE
          persp_pars <- if (!is.null(persp_sampler)) persp_sampler() else NULL
        }
        else if (forced == "bright") do_bright <- TRUE
        else if (forced == "contrast") do_contrast <- TRUE
        else if (forced == "noise") do_noise <- TRUE
        else if (forced == "blur") do_blur <- TRUE
        else if (forced == "cutout") do_cutout <- TRUE
      }
    }

    list(
      do_fh = do_fh, do_fv = do_fv, do_rot = do_rot, ang = ang,
      do_persp = do_persp, persp_pars = persp_pars,
      do_bright = do_bright, do_contrast = do_contrast, do_noise = do_noise,
      do_blur = do_blur, do_cutout = do_cutout
    )
  }

  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  flip_idx <- NULL
  yaml_src <- file.path(dir_path, "data.yaml")
  if (file.exists(yaml_src)) {
    y_lines_src <- readLines(yaml_src, warn = FALSE)
    f_idx <- grep("^\\s*flip_idx\\s*:", y_lines_src)
    if (length(f_idx) > 0L) {
      m <- regmatches(y_lines_src[f_idx[1]], regexpr("\\[[0-9 ,]+\\]", y_lines_src[f_idx[1]]))
      if (length(m) > 0L) {
        flip_idx <- as.integer(unlist(strsplit(gsub("[\\[\\] ]", "", m), ","))) + 1L
      }
    }
  }

  target_dir <- if (!is.null(out_dir)) normalizePath(out_dir[1], winslash = "/", mustWork = FALSE) else dir_path
  if (!is.null(out_dir) && target_dir != dir_path) {
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    for (item in list.files(dir_path, full.names = TRUE)) {
      file.copy(item, target_dir, recursive = TRUE, overwrite = TRUE)
    }
    yaml_f <- file.path(target_dir, "data.yaml")
    if (file.exists(yaml_f)) {
      y_lines <- readLines(yaml_f, warn = FALSE)
      path_idx <- grep("^path:", y_lines)
      if (length(path_idx) > 0L) {
        y_lines[path_idx[1]] <- paste0("path: ", target_dir)
        writeLines(y_lines, yaml_f)
      }
    }
  }

  if (!is.null(seed)) set.seed(seed)
  copies <- as.integer(max(1L, copies))

  # CASE A: CLASSIFICATION AUGMENTATION
  if (task == "classify") {
    tr_dir <- file.path(target_dir, "train")
    if (!dir.exists(tr_dir)) cli::cli_abort("Missing {.path train/} directory for classification dataset.")

    class_subdirs <- list.dirs(tr_dir, recursive = FALSE, full.names = FALSE)
    class_subdirs <- class_subdirs[nzchar(class_subdirs) & !startsWith(class_subdirs, ".")]

    cls_img_map <- list()
    for (cls in class_subdirs) {
      c_dir <- file.path(tr_dir, cls)
      fls <- list.files(c_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
      fls <- fls[!grepl("_aug\\d+", basename(fls))]
      cls_img_map[[cls]] <- fls
    }

    n_orig <- sum(vapply(cls_img_map, length, integer(1)))
    if (n_orig == 0L) {
      if (isTRUE(verbose)) cli::cli_alert_warning("No original training images found in {.path {tr_dir}}.")
      return(.wrap_aug(target_dir, "classify", classes, 0L))
    }

    if (!is.null(total_imgs)) {
      target_tot <- as.integer(total_imgs[1])
      if (target_tot <= n_orig) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Training split already contains {n_orig} images (>= total_imgs = {target_tot}). No augmented images generated.")
        }
        return(.wrap_aug(target_dir, "classify", classes, 0L))
      }
      target_new <- target_tot - n_orig
      class_targets <- vapply(class_subdirs, function(cls) {
        length(cls_img_map[[cls]]) / n_orig * target_new
      }, numeric(1))
      class_new <- floor(class_targets)
      rem_cls <- target_new - sum(class_new)
      if (rem_cls > 0L) {
        ord_rem <- order(class_targets - class_new, decreasing = TRUE)
        class_new[ord_rem[seq_len(rem_cls)]] <- class_new[ord_rem[seq_len(rem_cls)]] + 1L
      }
    }

    total_to_create <- if (!is.null(total_imgs)) sum(unlist(class_new)) else n_orig * copies
    show_pb <- isTRUE(verbose) && isTRUE(progress)
    if (show_pb) {
      pb <- cli::cli_progress_bar(
        name = "Augmenting dataset",
        total = total_to_create,
        format = "{cli::pb_spin} Augmenting dataset [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
        clear = FALSE
      )
    }

    total_created <- 0L
    for (cls in class_subdirs) {
      c_dir <- file.path(tr_dir, cls)
      fls <- cls_img_map[[cls]]
      n_c <- length(fls)
      if (n_c == 0L) next

      if (!is.null(total_imgs)) {
        c_target <- class_new[[cls]]
        c_base <- c_target %/% n_c
        c_rem <- c_target %% n_c
        img_copies <- rep(c_base, n_c)
        if (c_rem > 0L) {
          extra_idx <- sample(seq_len(n_c), c_rem)
          img_copies[extra_idx] <- img_copies[extra_idx] + 1L
        }
      } else {
        img_copies <- rep(copies, n_c)
      }

      for (i in seq_along(fls)) {
        n_cp <- img_copies[i]
        if (n_cp == 0L) next

        f <- fls[i]
        im <- try(image_import(f), silent = TRUE)
        if (inherits(im, "try-error") || !is_image(im)) next

        ext <- tools::file_ext(f)
        stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(f))

        for (cp in seq_len(n_cp)) {
          do_mos <- has_mosaic && (runif(1) < mosaic_p)
          if (do_mos && n_c >= 2L) {
            other_fls <- sample(fls, 3L, replace = (n_c < 4L))
            four_fls <- c(f, other_fls)
            cur_im <- image_mosaic(images = four_fls, overlap = mosaic_overlap,
                                   crop_region = mosaic_crop, blend = mosaic_blend,
                                   plot = FALSE, verbose = FALSE)
          } else {
            cur_im <- im
          }
          ops <- .sample_chosen_ops()
          do_fh <- ops$do_fh
          do_fv <- ops$do_fv
          do_rot <- ops$do_rot
          ang <- ops$ang
          do_persp <- ops$do_persp
          persp_pars <- ops$persp_pars
          do_bright <- ops$do_bright
          do_contrast <- ops$do_contrast
          do_noise <- ops$do_noise
          do_blur <- ops$do_blur
          do_cutout <- ops$do_cutout

          if (do_fh) {
            d <- dim(cur_im); attrs <- attributes(cur_im)
            cur_im <- cur_im[rev(seq_len(d[1])), , , drop = FALSE]
            attributes(cur_im) <- attrs
          }
          if (do_fv) {
            d <- dim(cur_im); attrs <- attributes(cur_im)
            cur_im <- cur_im[, rev(seq_len(d[2])), , drop = FALSE]
            attributes(cur_im) <- attrs
          }
          if (do_persp) {
            cur_im <- image_perspective(cur_im, x = persp_pars$px, y = persp_pars$py, z = persp_pars$pz,
                                        scale = persp_pars$pscale, bg_col = bg_col, keep_size = TRUE, plot = FALSE)
          }
          if (do_rot && ang != 0) {
            cur_im <- image_rotate(cur_im, angle = ang, bg_col = bg_col, plot = FALSE)
          }

          if (do_contrast) {
            c_lim <- if (length(contrast_range) == 1L) {
              as.numeric(contrast_range)
            } else if (min(contrast_range) < 1.0 && max(contrast_range) <= 1.5) {
              runif(1, 1.5, 3.5)
            } else {
              runif(1, min(contrast_range), max(contrast_range))
            }
            cur_im <- image_contrast(cur_im, clip_limit = c_lim, verbose = FALSE, plot = FALSE)
          }

          has_photo <- (do_blur || do_bright || do_noise || do_cutout)
          if (has_photo) {
            u <- unclass(cur_im)
            arr <- if (is.raw(u)) as.numeric(u) / 255.0 else as.numeric(u)
            dim(arr) <- dim(u)
            if (length(dim(arr)) == 2L) dim(arr) <- c(dim(arr)[1], dim(arr)[2], 1L)

            # Photometric & Blur transforms
            if (do_blur) {
              b_sig <- runif(1, min(blur_sigma), max(blur_sigma))
              arr <- gaussian_blur_cpp(arr, sigma = b_sig)
            }
            if (do_bright) {
              b_fac <- runif(1, min(brightness_range), max(brightness_range))
              arr <- image_brightness_cpp(arr, factor = b_fac, offset = 0.0)
            }
            if (do_noise) {
              sd_val <- runif(1, min(noise_sd), max(noise_sd))
              seed_val <- as.integer(sample.int(1000000L, 1L))
              arr <- image_noise_cpp(arr, type = 1L, mean = 0.0, sd = sd_val, p = 0.0, seed = seed_val)
            }
            if (do_cutout) {
              arr <- .apply_cutout(arr, raw_lines = character(0), target = "image",
                                   n_cutouts = cutout_n, size_range = cutout_size, fill = cutout_fill,
                                   pos = cutout_pos)
            }
            img_to_write <- arr
          } else {
            img_to_write <- unclass(cur_im)
          }

          applied_ops <- if (do_mos && n_c >= 2L) "mosaic" else character(0)
          if (do_fh) applied_ops <- c(applied_ops, "fliph")
          if (do_fv) applied_ops <- c(applied_ops, "flipv")
          if (do_persp) {
            px_r <- round(persp_pars$px)
            py_r <- round(persp_pars$py)
            px_str <- if (px_r < 0) paste0("_", abs(px_r)) else as.character(px_r)
            py_str <- if (py_r < 0) paste0("_", abs(py_r)) else as.character(py_r)
            applied_ops <- c(applied_ops, sprintf("persp_x%s_y%s", px_str, py_str))
          }
          if (do_rot && ang != 0) {
            ang_r <- round(ang)
            ang_str <- if (ang_r < 0) paste0("_", abs(ang_r)) else as.character(ang_r)
            applied_ops <- c(applied_ops, paste0("rot", ang_str))
          }
          if (do_bright) applied_ops <- c(applied_ops, "bright")
          if (do_contrast) applied_ops <- c(applied_ops, "contrast")
          if (do_noise) applied_ops <- c(applied_ops, "noise")
          if (do_blur) applied_ops <- c(applied_ops, "blur")
          if (do_cutout) applied_ops <- c(applied_ops, "cutout")
          ops_tag <- if (length(applied_ops) > 0L) paste0("_", paste(applied_ops, collapse = "_")) else ""
          ops_tag <- gsub("-", "_", ops_tag, fixed = TRUE)
          ops_tag <- gsub("_+", "_", ops_tag)

          new_fn <- paste0(stem, "_aug", cp, ops_tag, ".", ext)
          write_image_cpp(img_to_write, file.path(c_dir, new_fn), quality = quality)
          total_created <- total_created + 1L
          if (show_pb) cli::cli_progress_update(id = pb)
        }
      }
    }

    if (show_pb) cli::cli_progress_done(id = pb)

    if (isTRUE(verbose)) {
      tot_tr <- length(list.files(tr_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", recursive = TRUE, ignore.case = TRUE))
      cli::cli_alert_success("Generated {total_created} augmented images across training classes.")
      cli::cli_bullets(c(
        "*" = "Original training images: {n_orig}",
        "*" = "Augmented variants created: {total_created}",
        "*" = "Total training images now: {tot_tr}{if (!is.null(total_imgs)) paste0(' (target total_imgs: ', total_imgs, ')') else ''}",
        "v" = "Validation and test splits remained pristine (unaugmented)."
      ))
    }
    return(.wrap_aug(target_dir, "classify", classes, total_created))
  }

  # CASE B: DETECTION / SEGMENTATION AUGMENTATION
  tr_img_dir <- file.path(target_dir, "images", "train")
  tr_lbl_dir <- file.path(target_dir, "labels", "train")

  if (!dir.exists(tr_img_dir) || !dir.exists(tr_lbl_dir)) {
    cli::cli_abort("Missing {.path images/train/} or {.path labels/train/} directory.")
  }

  img_files <- list.files(tr_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  img_files <- img_files[!grepl("_aug\\d+", basename(img_files))]
  n_orig <- length(img_files)

  if (n_orig == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No original training images found in {.path {tr_img_dir}}.")
    return(.wrap_aug(target_dir, task, classes, 0L))
  }

  if (!is.null(total_imgs)) {
    target_tot <- as.integer(total_imgs[1])
    if (target_tot <= n_orig) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Training split already contains {n_orig} images (>= total_imgs = {target_tot}). No augmented images generated.")
      }
      return(.wrap_aug(target_dir, task, classes, 0L))
    }
    target_new <- target_tot - n_orig
    base_copies <- target_new %/% n_orig
    rem <- target_new %% n_orig
    img_copies <- rep(base_copies, n_orig)
    if (rem > 0L) {
      extra_idx <- sample(seq_len(n_orig), rem)
      img_copies[extra_idx] <- img_copies[extra_idx] + 1L
    }
  } else {
    img_copies <- rep(copies, n_orig)
  }

  total_to_create <- sum(img_copies)
  if (total_to_create == 0L) {
    return(.wrap_aug(target_dir, task, classes, 0L))
  }

  show_pb <- isTRUE(verbose) && isTRUE(progress)
  if (show_pb) {
    pb <- cli::cli_progress_bar(
      name = "Augmenting dataset",
      total = total_to_create,
      format = "{cli::pb_spin} Augmenting dataset [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      clear = FALSE
    )
  }

  total_created <- 0L

  for (i in seq_along(img_files)) {
    n_cp <- img_copies[i]
    if (n_cp == 0L) next

    im_f <- img_files[i]
    im <- try(image_import(im_f), silent = TRUE)
    if (inherits(im, "try-error") || !is_image(im)) next

    stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_f))
    ext <- tools::file_ext(im_f)
    lbl_f <- file.path(tr_lbl_dir, paste0(stem, ".txt"))
    raw_lines <- if (file.exists(lbl_f)) readLines(lbl_f, warn = FALSE) else character(0)
    raw_lines <- raw_lines[nzchar(trimws(raw_lines))]

    for (cp in seq_len(n_cp)) {
      do_mos <- has_mosaic && (runif(1) < mosaic_p)
      if (do_mos && n_orig >= 2L) {
        other_idx <- sample(seq_along(img_files), 3L, replace = (n_orig < 4L))
        four_idx <- c(i, other_idx)
        four_img_files <- img_files[four_idx]
        four_lbls <- lapply(four_idx, function(idx) {
          st_k <- sub("\\.[a-zA-Z0-9]+$", "", basename(img_files[idx]))
          lf_k <- file.path(tr_lbl_dir, paste0(st_k, ".txt"))
          if (file.exists(lf_k)) readLines(lf_k, warn = FALSE) else character(0)
        })
        mos_out <- image_mosaic(
          images = four_img_files,
          labels = four_lbls,
          overlap = mosaic_overlap,
          crop_region = mosaic_crop,
          blend = mosaic_blend,
          min_area_ratio = mosaic_min_area,
          plot = FALSE,
          verbose = FALSE
        )
        cur_im <- mos_out
        cur_lines <- attr(mos_out, "labels")
      } else {
        cur_im <- im
        cur_lines <- raw_lines
      }

      did_perm_pos <- NULL
      do_perm <- isTRUE(permute) && (runif(1) < permute_prob)
      if (do_perm && length(cur_lines) > 0L) {
        perm_res <- .permute_cluster_single(cur_im, cur_lines, position = permute_pos, fill = permute_fill, margin = permute_margin)
        if (isTRUE(perm_res$shifted)) {
          cur_im <- perm_res$image
          cur_lines <- perm_res$labels
          did_perm_pos <- perm_res$pos
        }
      }

      ops <- .sample_chosen_ops()
      do_fh <- ops$do_fh
      do_fv <- ops$do_fv
      do_rot <- ops$do_rot
      ang <- ops$ang
      do_persp <- ops$do_persp
      persp_pars <- ops$persp_pars
      do_bright <- ops$do_bright
      do_contrast <- ops$do_contrast
      do_noise <- ops$do_noise
      do_blur <- ops$do_blur
      do_cutout <- ops$do_cutout

      # 1. Flip horizontal
      if (do_fh) {
        d <- dim(cur_im); attrs <- attributes(cur_im)
        cur_im <- cur_im[rev(seq_len(d[1])), , , drop = FALSE]
        attributes(cur_im) <- attrs
        cur_lines <- .yolo_flip_annotations(cur_lines, do_h = TRUE, do_v = FALSE, flip_idx = flip_idx)
      }

      # 2. Flip vertical
      if (do_fv) {
        d <- dim(cur_im); attrs <- attributes(cur_im)
        cur_im <- cur_im[, rev(seq_len(d[2])), , drop = FALSE]
        attributes(cur_im) <- attrs
        cur_lines <- .yolo_flip_annotations(cur_lines, do_h = FALSE, do_v = TRUE, flip_idx = flip_idx)
      }

      # 3. 3D Perspective warp
      if (do_persp) {
        W_in <- dim(cur_im)[1]; H_in <- dim(cur_im)[2]
        cur_im <- image_perspective(cur_im, x = persp_pars$px, y = persp_pars$py, z = persp_pars$pz,
                                    scale = persp_pars$pscale, bg_col = bg_col, keep_size = TRUE, plot = FALSE)
        M_persp <- attr(cur_im, "homography")
        cur_lines <- .yolo_transform_annotations(cur_lines, W_in, H_in, dim(cur_im)[1], dim(cur_im)[2], M_persp)
      }

      # 4. Rotation (supports arbitrary angles, e.g. 45, 90, 180, 270)
      if (do_rot && ang != 0) {
        W_in <- dim(cur_im)[1]; H_in <- dim(cur_im)[2]
        cur_im <- image_rotate(cur_im, angle = ang, bg_col = bg_col, plot = FALSE)
        W_out <- dim(cur_im)[1]; H_out <- dim(cur_im)[2]
        M_rot <- .build_rot_matrix(W_in, H_in, W_out, H_out, ang)
        cur_lines <- .yolo_transform_annotations(cur_lines, W_in, H_in, W_out, H_out, M_rot)
      }

      if (do_contrast) {
        c_lim <- if (length(contrast_range) == 1L) {
          as.numeric(contrast_range)
        } else if (min(contrast_range) < 1.0 && max(contrast_range) <= 1.5) {
          runif(1, 1.5, 3.5)
        } else {
          runif(1, min(contrast_range), max(contrast_range))
        }
        cur_im <- image_contrast(cur_im, clip_limit = c_lim, verbose = FALSE, plot = FALSE)
      }

      # 5. Photometric & Blur transforms
      has_photo <- (do_blur || do_bright || do_noise || do_cutout)
      if (has_photo) {
        u <- unclass(cur_im)
        arr <- if (is.raw(u)) as.numeric(u) / 255.0 else as.numeric(u)
        dim(arr) <- dim(u)
        if (length(dim(arr)) == 2L) dim(arr) <- c(dim(arr)[1], dim(arr)[2], 1L)

        if (do_blur) {
          b_sig <- runif(1, min(blur_sigma), max(blur_sigma))
          arr <- gaussian_blur_cpp(arr, sigma = b_sig)
        }
        if (do_bright) {
          b_fac <- runif(1, min(brightness_range), max(brightness_range))
          arr <- image_brightness_cpp(arr, factor = b_fac, offset = 0.0)
        }
        if (do_noise) {
          sd_val <- runif(1, min(noise_sd), max(noise_sd))
          seed_val <- as.integer(sample.int(1000000L, 1L))
          arr <- image_noise_cpp(arr, type = 1L, mean = 0.0, sd = sd_val, p = 0.0, seed = seed_val)
        }
        if (do_cutout) {
          arr <- .apply_cutout(arr, raw_lines = cur_lines, target = cutout_target,
                               n_cutouts = cutout_n, size_range = cutout_size, fill = cutout_fill,
                               pos = cutout_pos)
        }
        img_to_write <- arr
      } else {
        img_to_write <- unclass(cur_im)
      }

      applied_ops <- if (do_mos && n_orig >= 2L) "mosaic" else character(0)
      if (!is.null(did_perm_pos)) applied_ops <- c(applied_ops, paste0("perm_", did_perm_pos))
      if (do_fh) applied_ops <- c(applied_ops, "fliph")
      if (do_fv) applied_ops <- c(applied_ops, "flipv")
      if (do_persp) {
        px_r <- round(persp_pars$px)
        py_r <- round(persp_pars$py)
        px_str <- if (px_r < 0) paste0("_", abs(px_r)) else as.character(px_r)
        py_str <- if (py_r < 0) paste0("_", abs(py_r)) else as.character(py_r)
        applied_ops <- c(applied_ops, sprintf("persp_x%s_y%s", px_str, py_str))
      }
      if (do_rot && ang != 0) {
        ang_r <- round(ang)
        ang_str <- if (ang_r < 0) paste0("_", abs(ang_r)) else as.character(ang_r)
        applied_ops <- c(applied_ops, paste0("rot", ang_str))
      }
      if (do_bright) applied_ops <- c(applied_ops, "bright")
      if (do_contrast) applied_ops <- c(applied_ops, "contrast")
      if (do_noise) applied_ops <- c(applied_ops, "noise")
      if (do_blur) applied_ops <- c(applied_ops, "blur")
      if (do_cutout) applied_ops <- c(applied_ops, paste0("cutout_", cutout_target))
      ops_tag <- if (length(applied_ops) > 0L) paste0("_", paste(applied_ops, collapse = "_")) else ""
      ops_tag <- gsub("-", "_", ops_tag, fixed = TRUE)
      ops_tag <- gsub("_+", "_", ops_tag)

      new_fn <- paste0(stem, "_aug", cp, ops_tag, ".", ext)
      new_txt <- paste0(stem, "_aug", cp, ops_tag, ".txt")

      write_image_cpp(img_to_write, file.path(tr_img_dir, new_fn), quality = quality)
      writeLines(cur_lines, file.path(tr_lbl_dir, new_txt))
      total_created <- total_created + 1L
      if (show_pb) cli::cli_progress_update(id = pb)
    }
  }

  if (show_pb) cli::cli_progress_done(id = pb)

  # Invalidate stale YOLO label caches
  unlink(file.path(target_dir, "labels", "train.cache"), force = TRUE)
  unlink(file.path(target_dir, "labels", "val.cache"), force = TRUE)
  unlink(file.path(target_dir, "labels", "test.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    tot_tr <- length(list.files(tr_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", ignore.case = TRUE))
    cli::cli_alert_success("Generated {total_created} augmented image and label pairs.")
    cli::cli_bullets(c(
      "*" = "Original training images: {n_orig}",
      "*" = "Augmented variants created: {total_created}",
      "*" = "Total training images now: {tot_tr}{if (!is.null(total_imgs)) paste0(' (target total_imgs: ', total_imgs, ')') else ''}",
      "v" = "Validation and test splits remained pristine (unaugmented).",
      "v" = "Stale YOLO dataset caches cleared."
    ))
  }

  # Optional negative background synthesis via puzzle scrambling
  do_scramble <- FALSE
  actual_scramble_prop <- scramble_prop
  if (is.numeric(scramble_background) && length(scramble_background) >= 1L && scramble_background[1L] > 0) {
    do_scramble <- TRUE
    actual_scramble_prop <- scramble_background[1L]
  } else if (isTRUE(scramble_background)) {
    do_scramble <- TRUE
  }

  if (do_scramble) {
    scramble_mode <- match.arg(scramble_mode)
    if (isTRUE(verbose)) cli::cli_rule("Negative Background Scrambling")
    tryCatch({
      yolo_dataset_scramble(
        data = target_dir,
        prop = actual_scramble_prop,
        grid = scramble_grid,
        mode = scramble_mode,
        quality = quality,
        seed = seed,
        as_job = FALSE,
        verbose = verbose
      )
    }, error = function(e) {
      cli::cli_alert_warning("Could not generate scrambled backgrounds: {e$message}")
    })
  }

  .wrap_aug(target_dir, task, classes, total_created)
}

#' @export
yolo_dataset_aug <- yolo_dataset_augment

#' @title Remove Augmented Images and Labels from a YOLO Dataset
#' @name yolo_dataset_augment_clean
#' @aliases yolo_dataset_aug_clean yolo_dataset_clean_aug
#' @description
#' `yolo_dataset_augment_clean()` removes all synthetic augmented image and label files
#' (matching the pattern `_aug[0-9]+`) previously generated by [yolo_dataset_augment()],
#' reverting the dataset back to its pristine original state.
#'
#' It automatically invalidates stale YOLO cache files (`train.cache`, `val.cache`, `test.cache`)
#' to ensure subsequent YOLO training or validation does not reference deleted files.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Augment dataset in-place
#' ds <- yolo_dataset_augment("dataset_complex", copies = 2)
#'
#' # 2. Clean up augmented variants and return to pristine state
#' ds_clean <- yolo_dataset_augment_clean("dataset_complex")
yolo_dataset_augment_clean <- function(data, verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"

  n_removed_img <- 0L
  n_removed_lbl <- 0L
  n_orig_img <- 0L

  if (task == "classify") {
    all_imgs <- list.files(dir_path, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    aug_imgs <- all_imgs[grepl("(_aug\\d+|_perm\\d+)", basename(all_imgs))]
    if (length(aug_imgs) > 0L) {
      unlink(aug_imgs, force = TRUE)
      n_removed_img <- length(aug_imgs)
    }
    n_orig_img <- length(all_imgs) - n_removed_img
  } else {
    # Detection or Segmentation
    img_dir <- file.path(dir_path, "images")
    lbl_dir <- file.path(dir_path, "labels")

    if (dir.exists(img_dir)) {
      all_imgs <- list.files(img_dir, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      aug_imgs <- all_imgs[grepl("(_aug\\d+|_perm\\d+)", basename(all_imgs))]
      if (length(aug_imgs) > 0L) {
        unlink(aug_imgs, force = TRUE)
        n_removed_img <- length(aug_imgs)
      }
      n_orig_img <- length(all_imgs) - n_removed_img
    }

    if (dir.exists(lbl_dir)) {
      all_lbls <- list.files(lbl_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      aug_lbls <- all_lbls[grepl("(_aug\\d+|_perm\\d+)", basename(all_lbls))]
      if (length(aug_lbls) > 0L) {
        unlink(aug_lbls, force = TRUE)
        n_removed_lbl <- length(aug_lbls)
      }
    }

    # Invalidate stale cache files
    unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
    unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
    unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
  }

  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No augmented files matching `_aug[0-9]+` found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned augmented files from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "*" = "Augmented images removed: {n_removed_img}",
        if (task != "classify") "*" = "Augmented label files removed: {n_removed_lbl}",
        "*" = "Remaining original images: {n_orig_img}",
        "v" = "YOLO label caches invalidated."
      ))
    }
  }

  invisible(.yolo_wrap_dataset(dir_path, task = task, classes = classes))
}

#' @export
yolo_dataset_aug_clean <- yolo_dataset_augment_clean

#' @export
yolo_dataset_clean_aug <- yolo_dataset_augment_clean

# ==============================================================================
# 4c. BREAK SPATIAL CENTER BIAS VIA CLUSTER PERMUTATION (yolo_dataset_permute)
# ==============================================================================
#' @title Relocate / Permute Target Object Clusters to Canvas Edges for YOLO Datasets
#' @name yolo_dataset_permute
#' @aliases yolo_dataset_relocate
#' @description
#' `yolo_dataset_permute()` (and alias `yolo_dataset_relocate()`) breaks spatial position bias
#' (Center Bias) in YOLO object detection and instance segmentation datasets by identifying
#' the bounding cluster/window of target objects (such as a stack of wood logs loaded on a truck)
#' and relocating/permuting it towards the canvas borders and corners (`"left"`, `"right"`,
#' `"top"`, `"bottom"`, `"top_left"`, `"top_right"`, `"bottom_left"`, `"bottom_right"`, or `"random"`).
#'
#' The hole left in the original object location is seamlessly filled using authentic slices/patches
#' sampled from the outer image borders (`fill = "border_slices"`), preserving authentic background
#' texture, lighting, and sensor noise while training the neural network that objects can appear
#' anywhere across the canvas.
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param out_dir Character. Directory to save the augmented dataset. If `NULL` (default),
#'   generates permuted variants in-place within the current dataset using suffix `_perm<copy>_<pos>`.
#' @param copies Integer. Number of permuted variants to generate per training image (default: `1L`).
#' @param total_imgs Integer. Target total number of training images (default: `NULL`).
#' @param positions Character vector of allowed target positions: `"random"` (default),
#'   `"left"`, `"right"`, `"top"`, `"bottom"`, `"corner"`, `"top_left"`, `"top_right"`,
#'   `"bottom_left"`, `"bottom_right"`.
#' @param fill Character specifying the hole filling method: `"border_slices"` (default, uses authentic
#'   patches from image borders), `"mean"` (mean border color + noise), `"blur"`, or `"noise"`.
#' @param margin Numeric in `[0, 0.2]`. Margin padding around the cluster bounding box (default: `0.02`).
#' @param blend_edge Integer. Pixel radius for feathered boundary blending of the filled hole (default: `5L`).
#' @param quality Integer (1-100). JPEG compression quality (default: `90`).
#' @param seed Optional integer seed for reproducibility (default: `42`).
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs as a background job (default: `FALSE`).
#' @param verbose Logical. Whether to show progress messages (default: `TRUE`).
#' @param progress Logical. Whether to show a progress bar (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the augmented dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Generate 1 permuted variant per image (randomly shifted to borders/corners)
#' ds_perm <- yolo_dataset_permute("my_dataset")
#'
#' # 2. Shift clusters specifically to the left and right borders
#' ds_sides <- yolo_dataset_permute("my_dataset", positions = c("left", "right"), copies = 2)
#'
#' # 3. Revert and clean up permuted images
#' yolo_dataset_permute_clean("my_dataset")
yolo_dataset_permute <- function(data,
                                 out_dir = NULL,
                                 copies = 1L,
                                 total_imgs = NULL,
                                 positions = c("random", "left", "right", "top", "bottom", "corner", "top_left", "top_right", "bottom_left", "bottom_right"),
                                 fill = c("border_slices", "mean", "blur", "noise"),
                                 margin = 0.02,
                                 blend_edge = 5L,
                                 quality = 90,
                                 seed = 42,
                                 as_job = FALSE,
                                 verbose = TRUE,
                                 progress = TRUE) {
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding in current session.")
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
      job_script_path <- file.path(job_dir, paste0("yolo_perm_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_perm_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        d_p <- .yolo_get_dir(data)
        basename(d_p)
      }, error = function(e) {
        if (inherits(data, "yolo_dataset")) basename(data$dir)
        else if (is.character(data) && nzchar(data[1])) basename(data[1])
        else "dataset"
      })
      j_name <- paste0("YOLO Permute: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_permute()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Cluster Permutation - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_permute, job_info$args)\n",
        "target_dir <- if (!is.null(res$dir)) res$dir else job_info$args$out_dir\n",
        "if (!is.null(target_dir) && dir.exists(target_dir)) {\n",
        "  save_path <- file.path(target_dir, 'dataset_permuted.rds')\n",
        "  saveRDS(res, save_path)\n",
        "  cat('\\n[OK] Permuted dataset object successfully saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Permutation Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Cluster Permutation (pliman)")
      cli::cli_alert_success("Cluster permutation has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target dataset: {.val {target_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track live progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  fill <- if (length(fill) > 1L) fill[1] else fill
  if (!fill %in% c("border_slices", "mean", "blur", "noise")) fill <- "border_slices"

  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  if (task == "classify") {
    cli::cli_abort("Classification datasets do not contain bounding boxes or polygons to permute. yolo_dataset_permute() applies to detection and segmentation datasets.")
  }

  target_dir <- if (!is.null(out_dir)) normalizePath(out_dir[1], winslash = "/", mustWork = FALSE) else dir_path
  if (!is.null(out_dir) && target_dir != dir_path) {
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    for (item in list.files(dir_path, full.names = TRUE)) {
      file.copy(item, target_dir, recursive = TRUE, overwrite = TRUE)
    }
    yaml_f <- file.path(target_dir, "data.yaml")
    if (file.exists(yaml_f)) {
      y_lines <- readLines(yaml_f, warn = FALSE)
      path_idx <- grep("^path:", y_lines)
      if (length(path_idx) > 0L) {
        y_lines[path_idx[1]] <- paste0("path: ", target_dir)
        writeLines(y_lines, yaml_f)
      }
    }
  }

  if (!is.null(seed)) set.seed(seed)
  copies <- as.integer(max(1L, copies))

  tr_img_dir <- file.path(target_dir, "images", "train")
  tr_lbl_dir <- file.path(target_dir, "labels", "train")

  if (!dir.exists(tr_img_dir) || !dir.exists(tr_lbl_dir)) {
    cli::cli_abort("Missing {.path images/train/} or {.path labels/train/} directory in {.path {target_dir}}.")
  }

  img_files <- list.files(tr_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  img_files <- img_files[!grepl("(_aug\\d+|_perm\\d+)", basename(img_files))]
  n_orig <- length(img_files)

  if (n_orig == 0L) {
    if (isTRUE(verbose)) cli::cli_alert_warning("No original training images found in {.path {tr_img_dir}}.")
    return(invisible(.yolo_wrap_dataset(target_dir, task = task, classes = classes)))
  }

  if (!is.null(total_imgs)) {
    target_tot <- as.integer(total_imgs[1])
    if (target_tot <= n_orig) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Training split already contains {n_orig} images (>= total_imgs = {target_tot}). No permuted images generated.")
      }
      return(invisible(.yolo_wrap_dataset(target_dir, task = task, classes = classes)))
    }
    target_new <- target_tot - n_orig
    base_copies <- target_new %/% n_orig
    rem <- target_new %% n_orig
    img_copies <- rep(base_copies, n_orig)
    if (rem > 0L) {
      extra_idx <- sample(seq_len(n_orig), rem)
      img_copies[extra_idx] <- img_copies[extra_idx] + 1L
    }
  } else {
    img_copies <- rep(copies, n_orig)
  }

  total_to_create <- sum(img_copies)
  if (total_to_create == 0L) {
    return(invisible(.yolo_wrap_dataset(target_dir, task = task, classes = classes)))
  }

  show_pb <- isTRUE(verbose) && isTRUE(progress)
  if (show_pb) {
    pb <- cli::cli_progress_bar(
      name = "Permuting object clusters",
      total = total_to_create,
      format = "{cli::pb_spin} Permuting object clusters [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      clear = FALSE
    )
  }

  total_created <- 0L

  for (i in seq_along(img_files)) {
    n_cp <- img_copies[i]
    if (n_cp == 0L) next

    im_f <- img_files[i]
    im <- try(image_import(im_f), silent = TRUE)
    if (inherits(im, "try-error") || !is_image(im)) next

    stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_f))
    ext <- tools::file_ext(im_f)
    lbl_f <- file.path(tr_lbl_dir, paste0(stem, ".txt"))
    raw_lines <- if (file.exists(lbl_f)) readLines(lbl_f, warn = FALSE) else character(0)
    raw_lines <- raw_lines[nzchar(trimws(raw_lines))]

    if (length(raw_lines) == 0L) {
      next
    }

    for (cp in seq_len(n_cp)) {
      cur_pos <- if (length(positions) > 1L) sample(positions, 1L) else positions[1]
      perm_res <- .permute_cluster_single(
        im = im,
        raw_lines = raw_lines,
        position = cur_pos,
        fill = fill,
        margin = margin,
        blend_edge = blend_edge
      )

      if (isTRUE(perm_res$shifted)) {
        pos_tag <- paste0("perm", cp, "_", perm_res$pos)
        new_fn <- paste0(stem, "_", pos_tag, ".", ext)
        new_txt <- paste0(stem, "_", pos_tag, ".txt")

        write_image_cpp(unclass(perm_res$image), file.path(tr_img_dir, new_fn), quality = quality)
        writeLines(perm_res$labels, file.path(tr_lbl_dir, new_txt))
        total_created <- total_created + 1L
      }
      if (show_pb) cli::cli_progress_update(id = pb)
    }
  }

  if (show_pb) cli::cli_progress_done(id = pb)

  # Invalidate stale label caches
  unlink(file.path(target_dir, "labels", "train.cache"), force = TRUE)
  unlink(file.path(target_dir, "labels", "val.cache"), force = TRUE)
  unlink(file.path(target_dir, "labels", "test.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    tot_tr <- length(list.files(tr_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", ignore.case = TRUE))
    cli::cli_alert_success("Generated {total_created} cluster-permuted images to break spatial position bias.")
    cli::cli_bullets(c(
      "*" = "Original training images: {n_orig}",
      "*" = "Permuted variants created: {total_created}",
      "*" = "Total training images now: {tot_tr}{if (!is.null(total_imgs)) paste0(' (target total_imgs: ', total_imgs, ')') else ''}",
      "v" = "Holes filled seamlessly using authentic {.val {fill}} background synthesis.",
      "v" = "Validation and test splits remained pristine (unaugmented)."
    ))
  }

  invisible(.yolo_wrap_dataset(target_dir, task = task, classes = classes))
}

#' @export
yolo_dataset_relocate <- yolo_dataset_permute

#' @title Clean Permuted Augmented Files from a YOLO Dataset
#' @name yolo_dataset_permute_clean
#' @aliases yolo_dataset_clean_perm
#' @description
#' `yolo_dataset_permute_clean()` removes all synthetic permuted image and label files
#' (matching the pattern `_perm[0-9]+`) previously generated by [yolo_dataset_permute()],
#' returning the training split back to its pristine original state.
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
yolo_dataset_permute_clean <- function(data, verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  n_removed_img <- 0L
  n_removed_lbl <- 0L

  img_dir <- file.path(dir_path, "images")
  lbl_dir <- file.path(dir_path, "labels")

  if (dir.exists(img_dir)) {
    all_imgs <- list.files(img_dir, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    perm_imgs <- all_imgs[grepl("_perm\\d+", basename(all_imgs))]
    if (length(perm_imgs) > 0L) {
      unlink(perm_imgs, force = TRUE)
      n_removed_img <- length(perm_imgs)
    }
  }

  if (dir.exists(lbl_dir)) {
    all_lbls <- list.files(lbl_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    perm_lbls <- all_lbls[grepl("_perm\\d+", basename(all_lbls))]
    if (length(perm_lbls) > 0L) {
      unlink(perm_lbls, force = TRUE)
      n_removed_lbl <- length(perm_lbls)
    }
  }

  # Invalidate stale cache files
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No permuted files matching `_perm[0-9]+` found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned {n_removed_img} permuted images and {n_removed_lbl} labels from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "v" = "YOLO label caches invalidated."
      ))
    }
  }

  invisible(.yolo_wrap_dataset(dir_path, task = task, classes = classes))
}

#' @export
yolo_dataset_clean_perm <- yolo_dataset_permute_clean


# ==============================================================================
# 4b. NEGATIVE BACKGROUND SYNTHESIS VIA SCRAMBLING (yolo_dataset_scramble)
# ==============================================================================
#' @title Generate Negative Background Images via Image Scrambling for YOLO Datasets
#' @name yolo_dataset_scramble
#' @aliases yolo_dataset_shuffle
#' @description
#' `yolo_dataset_scramble()` (and alias `yolo_dataset_shuffle()`) generates authentic
#' negative background images (empty label files with zero objects) directly from existing
#' images in a YOLO dataset using patch scrambling / jigsaw tile permutation.
#'
#' By dividing original images into a grid of tiles (e.g. 8x8 or 16x16), randomly shuffling
#' the tile locations, and optionally applying random 90-degree rotations and flips to individual
#' tiles, the image retains 100% of the true domain characteristics:
#' - Authentic sensor noise, ISO grain, and camera compression
#' - Exact spectral color distribution, white balance, and contrast
#' - Natural scene lighting and shadow hues
#'
#' Crucially, the puzzle scrambling completely destroys all coherent object contours,
#' aspect ratios, and spatial semantics. Providing these images with empty `.txt` labels
#' instructs the YOLO object detector that random textures and color patches without
#' structured object contours must produce zero detections, drastically suppressing false-positive
#' hallucinations in production.
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param n Integer. Exact number of scrambled background images to generate. If `NULL`,
#'   `prop` is used. Defaults to `NULL`.
#' @param prop Numeric between 0 and 1. Proportion of original images to generate as scrambled
#'   backgrounds (e.g., `0.05` for 5% negatives, following Ultralytics recommendation of 2-10%).
#'   Defaults to `0.05`.
#' @param grid Integer vector of length 1 or 2. Grid dimensions for tile slicing (e.g. `c(8, 8)`
#'   or `8`). Defaults to `c(8, 8)`. Ignored when `mode = "pixel"`.
#' @param mode Character specifying the scrambling strategy: `"grid"` (default jigsaw tile shuffle),
#'   `"block_swap"` (local block swapping preserving partial background continuity), or `"pixel"` (global pixel shuffle).
#' @param swap_prob Numeric between 0 and 1. Fraction of tiles swapped when `mode = "block_swap"` (default: `0.5`).
#' @param shuffle_pixels Logical. Whether to randomize/shuffle all pixels internally within
#'   each individual tile (default: `TRUE`). Combining tile permutation with intra-tile pixel
#'   randomization completely dissolves any recognizable sub-objects (e.g. wood logs) while
#'   strictly conserving the local background color and textural palette of each region!
#' @param rotate_tiles Logical. Whether to apply random 90-degree rotations (0°, 90°, 180°, 270°)
#'   to individual scrambled tiles (default: `TRUE`). Ignored when `mode = "pixel"`.
#' @param flip_tiles Logical. Whether to apply random horizontal/vertical flips to individual
#'   scrambled tiles (default: `TRUE`). Ignored when `mode = "pixel"`.
#' @param train_prop Numeric. Proportion of scrambled background images placed in the `train` split
#'   vs `val` split. Defaults to `0.80`.
#' @param prefix Character. Filename prefix for the generated scrambled background files (default: `"scramble_"`).
#'   Each generated image and label is named `<prefix><uuid>.<ext>` (e.g. `scramble_6dfb8eb1fe94439dac814f8b1fbe1b45.jpg`),
#'   guaranteeing unique identifiers and preventing accidental file overwrites when running repeatedly.
#' @param quality Integer (1-100). JPEG compression quality (default: `90`).
#' @param seed Optional integer seed for reproducibility (default: `42`).
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs as a background job (default: `FALSE`).
#' @param verbose Logical. Whether to show progress messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the augmented dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Tile permutation with intra-tile pixel randomization (default)
#' yolo_dataset_scramble("my_dataset", prop = 0.05, grid = c(8, 8), shuffle_pixels = TRUE)
#'
#' # 2. Local block swap preserving partial background continuity
#' yolo_dataset_scramble("my_dataset", prop = 0.05, grid = c(10, 10), mode = "block_swap")
#'
#' # 3. Full pixel shuffle (exact spectral color noise)
#' yolo_dataset_scramble("my_dataset", prop = 0.05, mode = "pixel")
#'
#' # Undo/clean scrambled backgrounds
#' yolo_dataset_scramble_undo("my_dataset")
yolo_dataset_scramble <- function(data,
                                  n = NULL,
                                  prop = 0.05,
                                  grid = c(8, 8),
                                  mode = c("grid", "block_swap", "pixel"),
                                  swap_prob = 0.5,
                                  shuffle_pixels = TRUE,
                                  rotate_tiles = TRUE,
                                  flip_tiles = TRUE,
                                  train_prop = 0.80,
                                  prefix = "scramble_",
                                  quality = 90,
                                  seed = 42,
                                  as_job = FALSE,
                                  verbose = TRUE) {
  mode <- match.arg(mode)
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with scrambling in current session.")
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
      job_script_path <- file.path(job_dir, paste0("yolo_scramble_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_scramble_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (inherits(data, "yolo_dataset")) basename(data$dir)
        else if (is.character(data) && nzchar(data[1])) basename(data[1])
        else "dataset"
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Scramble: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_scramble()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Dataset Scramble - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_scramble, job_info$args)\n",
        "cat('\\n[OK] YOLO Dataset Scramble Background Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Dataset Scramble (pliman)")
      cli::cli_alert_success("Dataset scrambling has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target dataset: {.val {target_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  # Check standard YOLO directories
  if (dir.exists(file.path(dir_path, "train", "images"))) {
    train_img_dir <- file.path(dir_path, "train", "images")
    val_img_dir   <- file.path(dir_path, "val", "images")
    train_lbl_dir <- file.path(dir_path, "train", "labels")
    val_lbl_dir   <- file.path(dir_path, "val", "labels")
  } else {
    train_img_dir <- file.path(dir_path, "images", "train")
    val_img_dir   <- file.path(dir_path, "images", "val")
    train_lbl_dir <- file.path(dir_path, "labels", "train")
    val_lbl_dir   <- file.path(dir_path, "labels", "val")
  }

  for (d in c(train_img_dir, val_img_dir, train_lbl_dir, val_lbl_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  # Ignore existing scrambles, synth backgrounds, or augments when sampling source originals
  ignore_patt <- sprintf("^(%s|synth_bg_|bg_|_aug)", prefix)

  all_tr_imgs <- list.files(train_img_dir, pattern = img_patt, full.names = TRUE, ignore.case = TRUE)
  tr_sources <- all_tr_imgs[!grepl(ignore_patt, basename(all_tr_imgs), ignore.case = TRUE)]

  all_val_imgs <- list.files(val_img_dir, pattern = img_patt, full.names = TRUE, ignore.case = TRUE)
  val_sources <- all_val_imgs[!grepl(ignore_patt, basename(all_val_imgs), ignore.case = TRUE)]

  all_sources <- c(tr_sources, val_sources)
  if (length(all_sources) == 0L) {
    all_sources <- c(all_tr_imgs, all_val_imgs)
  }
  if (length(all_sources) == 0L) {
    cli::cli_abort("No source images found in {.path {dir_path}} to scramble.")
  }

  # Calculate target counts
  if (!is.null(n)) {
    total_scramble <- as.integer(n)
  } else {
    prop <- max(0.001, min(1.0, as.numeric(prop)))
    total_scramble <- max(1L, round(length(all_sources) * prop))
  }

  if (total_scramble <= 0L) total_scramble <- 1L

  n_train <- max(1L, round(total_scramble * train_prop))
  if (n_train >= total_scramble && total_scramble > 1L) n_train <- total_scramble - 1L
  n_val <- total_scramble - n_train

  # Sample source images
  set.seed(if (!is.null(seed)) seed else sample.int(2147483647L, 1L))
  sample_tr <- if (length(tr_sources) > 0L) {
    sample(tr_sources, size = n_train, replace = (n_train > length(tr_sources)))
  } else {
    sample(all_sources, size = n_train, replace = (n_train > length(all_sources)))
  }

  sample_val <- if (n_val > 0L) {
    if (length(val_sources) > 0L) {
      sample(val_sources, size = n_val, replace = (n_val > length(val_sources)))
    } else {
      sample(all_sources, size = n_val, replace = (n_val > length(all_sources)))
    }
  } else {
    character(0)
  }

  if (length(grid) == 1L) {
    grid_r <- as.integer(grid[1L])
    grid_c <- grid_r
  } else {
    grid_r <- as.integer(grid[1L])
    grid_c <- as.integer(grid[2L])
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Generating {total_scramble} scrambled negative background images ({n_train} train / {n_val} val)...")
    pb <- cli::cli_progress_bar(
      name = "Scrambling images",
      total = total_scramble,
      format = "{cli::pb_spin} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
    )
  }

  created_train <- character()
  created_val <- character()

  uids_train <- tolower(gsub("-", "", uuid(length(sample_tr))))
  uids_val   <- if (length(sample_val) > 0L) tolower(gsub("-", "", uuid(length(sample_val)))) else character(0)

  # Helper to process a slice
  process_slice <- function(src_list, uids_list, dest_img_dir, dest_lbl_dir, split_tag) {
    added <- character()
    for (i in seq_along(src_list)) {
      src_file <- src_list[i]
      uid_i <- uids_list[i]
      src_stem <- tools::file_path_sans_ext(basename(src_file))
      ext <- tools::file_ext(src_file)
      if (!nzchar(ext)) ext <- "jpg"

      cur_seed <- if (!is.null(seed)) (seed + i * 1013L) %% 2147483647L else sample.int(2147483647L, 1L)

      img_obj <- tryCatch(read_image_cpp(src_file), error = function(e) NULL)
      if (is.null(img_obj)) {
        if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
        next
      }

      scrambled <- image_scramble_cpp(
        img_sexp = img_obj,
        grid_r = grid_r,
        grid_c = grid_c,
        rotate_tiles = isTRUE(rotate_tiles),
        flip_tiles = isTRUE(flip_tiles),
        mode = mode,
        swap_prob = as.numeric(swap_prob),
        shuffle_pixels = isTRUE(shuffle_pixels),
        seed = cur_seed
      )

      out_fn <- sprintf("%s%s.%s", prefix, uid_i, ext)
      out_lbl <- sprintf("%s%s.txt", prefix, uid_i)

      dest_img_path <- file.path(dest_img_dir, out_fn)
      dest_lbl_path <- file.path(dest_lbl_dir, out_lbl)

      write_image_cpp(scrambled, dest_img_path, quality = quality)
      file.create(dest_lbl_path) # Empty label for negative background sample

      added <- c(added, dest_img_path)
      if (isTRUE(verbose)) cli::cli_progress_update(id = pb)
    }
    added
  }

  created_train <- process_slice(sample_tr, uids_train, train_img_dir, train_lbl_dir, "train")
  if (n_val > 0L) {
    created_val <- process_slice(sample_val, uids_val, val_img_dir, val_lbl_dir, "val")
  }

  if (isTRUE(verbose)) cli::cli_progress_done(id = pb)

  # Invalidate stale YOLO label caches
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)

  tot_tr <- length(list.files(train_img_dir, pattern = img_patt, ignore.case = TRUE))
  tot_val <- length(list.files(val_img_dir, pattern = img_patt, ignore.case = TRUE))
  tot_created <- length(created_train) + length(created_val)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Generated {tot_created} scrambled negative background images ({length(created_train)} train / {length(created_val)} val).")
    bullets <- c("*" = sprintf("Scrambling mode: %s", mode))
    if (mode != "pixel") {
      bullets <- c(bullets, "*" = sprintf("Grid divisions: %dx%d (%d jigsaw tiles per image)", grid_r, grid_c, grid_r * grid_c))
      bullets <- c(bullets, "*" = sprintf("Intra-tile pixel randomization: %s", isTRUE(shuffle_pixels)))
      bullets <- c(bullets, "*" = sprintf("Random tile rotations & flips: %s", isTRUE(rotate_tiles) || isTRUE(flip_tiles)))
    }
    if (mode == "block_swap") {
      bullets <- c(bullets, "*" = sprintf("Swap probability: %.2f", swap_prob))
    }
    bullets <- c(bullets,
      "*" = sprintf("Training images now: %d (+%d scrambled)", tot_tr, length(created_train)),
      "*" = sprintf("Validation images now: %d (+%d scrambled)", tot_val, length(created_val)),
      "v" = "Matching empty .txt label files created (0 objects) to suppress false positives.",
      "v" = "Stale YOLO dataset caches cleared."
    )
    cli::cli_bullets(bullets)
  }

  invisible(.yolo_wrap_dataset(dir_path, task = task, classes = classes))
}

#' @export
yolo_dataset_shuffle <- yolo_dataset_scramble

#' @title Undo and Remove Scrambled Background Images from a YOLO Dataset
#' @name yolo_dataset_scramble_undo
#' @aliases yolo_dataset_shuffle_undo yolo_dataset_scramble_clean yolo_dataset_clean_scramble
#' @description
#' `yolo_dataset_scramble_undo()` removes all scrambled background images and matching
#' empty label files (matching `^scramble_` or a custom `prefix`) generated by
#' [yolo_dataset_scramble()], reverting the dataset back to its previous state.
#'
#' It also automatically invalidates stale YOLO cache files (`train.cache`, `val.cache`, `test.cache`).
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param prefix Character. Filename prefix of the scrambled files to remove (default: `"scramble_"`).
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' yolo_dataset_scramble_undo("my_dataset")
yolo_dataset_scramble_undo <- function(data, prefix = "scramble_", verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  patt_regex <- paste0("^", prefix)
  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"

  n_removed_img <- 0L
  n_removed_lbl <- 0L
  n_orig_img <- 0L

  img_dir <- file.path(dir_path, "images")
  search_img_dir <- if (dir.exists(img_dir)) img_dir else dir_path
  all_imgs <- list.files(search_img_dir, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  scr_imgs <- all_imgs[grepl(patt_regex, basename(all_imgs), ignore.case = TRUE)]

  if (length(scr_imgs) > 0L) {
    unlink(scr_imgs, force = TRUE)
    n_removed_img <- length(scr_imgs)
  }
  n_orig_img <- length(all_imgs) - n_removed_img

  lbl_dirs <- c(
    file.path(dir_path, "labels"),
    file.path(dir_path, "train", "labels"),
    file.path(dir_path, "val", "labels"),
    file.path(dir_path, "test", "labels")
  )
  lbl_dirs <- lbl_dirs[dir.exists(lbl_dirs)]
  search_lbl_dir <- if (length(lbl_dirs) > 0L) lbl_dirs else dir_path
  all_lbls <- list.files(search_lbl_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  scr_lbls <- all_lbls[grepl(patt_regex, basename(all_lbls), ignore.case = TRUE)]
  if (length(scr_lbls) > 0L) {
    unlink(scr_lbls, force = TRUE)
    n_removed_lbl <- length(scr_lbls)
  }
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
  unlink(file.path(dir_path, "train", "labels.cache"), force = TRUE)
  unlink(file.path(dir_path, "val", "labels.cache"), force = TRUE)
  unlink(file.path(dir_path, "test", "labels.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No scrambled files matching prefix `^{prefix}` found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned scrambled background files from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "*" = "Scrambled images removed: {n_removed_img}",
        "*" = "Scrambled empty label files removed: {n_removed_lbl}",
        "*" = "Remaining images: {n_orig_img}",
        "v" = "YOLO label caches invalidated."
      ))
    }
  }

  invisible(.yolo_wrap_dataset(dir_path, task = task, classes = classes))
}

#' @export
yolo_dataset_shuffle_undo <- yolo_dataset_scramble_undo

#' @export
yolo_dataset_scramble_clean <- yolo_dataset_scramble_undo

#' @export
yolo_dataset_clean_scramble <- yolo_dataset_scramble_undo


# ==============================================================================
# 5. CLASS FILTERING & ISOLATION (yolo_dataset_filter)
# ==============================================================================
#' @title Filter or Subset Classes in a YOLO Dataset
#' @name yolo_dataset_filter
#' @aliases yolo_dataset_subset
#' @description
#' `yolo_dataset_filter()` isolates specific classes of interest from a multi-class
#' dataset (or drops unwanted categories). It rewrites all `.txt` label files, re-indexes
#' surviving class IDs contiguously from 0 to K-1, updates `data.yaml`, and optionally
#' cleans or preserves images that now have zero annotations (as background negatives).
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param keep_classes Character vector of class names to retain. Mutually exclusive with `drop_classes`.
#' @param drop_classes Character vector of class names to exclude.
#' @param out_dir Output directory path for the filtered dataset. Defaults to `"filtered_yolo_dataset"`.
#' @param drop_empty_images Logical. If `TRUE` (default), images left with zero annotations are deleted.
#'   If `FALSE`, they are kept as negative background images with empty `.txt` labels.
#' @param copy Logical. Whether to copy files to `out_dir`. Defaults to `TRUE`.
#' @param overwrite Logical. Whether to overwrite `out_dir` if it exists. Defaults to `FALSE`.
#' @param verbose Logical. Whether to show progress. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` object pointing to the filtered dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Focus exclusively on 2 classes from a 10-class dataset
#' ds_focused <- yolo_dataset_filter(
#'   "dataset_geral",
#'   keep_classes = c("ferrugem", "mancha_angular"),
#'   out_dir      = "dataset_doencas"
#' )
yolo_dataset_filter <- function(data,
                                keep_classes = NULL,
                                drop_classes = NULL,
                                out_dir = "filtered_yolo_dataset",
                                drop_empty_images = TRUE,
                                copy = TRUE,
                                overwrite = FALSE,
                                verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  if (is.null(keep_classes) && is.null(drop_classes)) {
    cli::cli_abort("Must specify either {.arg keep_classes} or {.arg drop_classes}.")
  }

  surviving_classes <- classes
  if (!is.null(keep_classes)) {
    surviving_classes <- intersect(classes, keep_classes)
  }
  if (!is.null(drop_classes)) {
    surviving_classes <- setdiff(surviving_classes, drop_classes)
  }

  if (length(surviving_classes) == 0L) {
    cli::cli_abort("Filtering criteria resulted in 0 surviving classes.")
  }

  out_dir_norm <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)
  if (dir.exists(out_dir_norm)) {
    if (isTRUE(overwrite)) unlink(out_dir_norm, recursive = TRUE, force = TRUE)
    else cli::cli_abort("Directory {.path {out_dir_norm}} already exists. Set {.code overwrite = TRUE}.")
  }
  dir.create(out_dir_norm, recursive = TRUE, showWarnings = FALSE)

  # Mapping from old class index to new class index
  id_map <- integer()
  for (i in seq_along(classes)) {
    old_c <- classes[i]
    if (old_c %in% surviving_classes) {
      new_idx <- which(surviving_classes == old_c) - 1L
      id_map[as.character(i - 1L)] <- new_idx
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Class Filtering (pliman)")
    cli::cli_alert_info("Filtered classes: from {length(classes)} down to {length(surviving_classes)}: {.val {surviving_classes}}")
  }

  # CASE A: CLASSIFICATION
  if (task == "classify") {
    splits <- c("train", "val", "test")
    for (s in splits) {
      s_dir <- file.path(dir_path, s)
      if (!dir.exists(s_dir)) next
      for (sc in surviving_classes) {
        c_src <- file.path(s_dir, sc)
        if (!dir.exists(c_src)) next
        c_dest <- file.path(out_dir_norm, s, sc)
        dir.create(c_dest, recursive = TRUE, showWarnings = FALSE)
        fls <- list.files(c_src, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
        for (f in fls) {
          if (isTRUE(copy)) file.copy(f, file.path(c_dest, basename(f)), overwrite = TRUE) else file.rename(f, file.path(c_dest, basename(f)))
        }
      }
    }
    return(.yolo_wrap_dataset(out_dir_norm, task = "classify", classes = surviving_classes))
  }

  # CASE B: DETECTION / SEGMENTATION
  splits <- c("train", "val", "test")
  for (s in splits) {
    src_img_dir <- file.path(dir_path, "images", s)
    src_lbl_dir <- file.path(dir_path, "labels", s)
    dest_img_dir <- file.path(out_dir_norm, "images", s)
    dest_lbl_dir <- file.path(out_dir_norm, "labels", s)

    if (!dir.exists(src_img_dir)) next
    dir.create(dest_img_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(dest_lbl_dir, recursive = TRUE, showWarnings = FALSE)

    img_files <- list.files(src_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
    for (im_f in img_files) {
      stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_f))
      lbl_f <- file.path(src_lbl_dir, paste0(stem, ".txt"))
      kept_lines <- character()

      if (file.exists(lbl_f)) {
        lines <- readLines(lbl_f, warn = FALSE)
        lines <- lines[nzchar(trimws(lines))]
        for (ln in lines) {
          parts <- strsplit(trimws(ln), "\\s+")[[1]]
          old_cid <- parts[1]
          if (old_cid %in% names(id_map)) {
            parts[1] <- as.character(id_map[[old_cid]])
            kept_lines <- c(kept_lines, paste(parts, collapse = " "))
          }
        }
      }

      if (length(kept_lines) == 0L && isTRUE(drop_empty_images)) {
        next
      }

      if (isTRUE(copy)) {
        file.copy(im_f, file.path(dest_img_dir, basename(im_f)), overwrite = TRUE)
      } else {
        file.rename(im_f, file.path(dest_img_dir, basename(im_f)))
      }
      writeLines(kept_lines, file.path(dest_lbl_dir, paste0(stem, ".txt")))
    }
  }

  .yolo_write_yaml(out_dir_norm, classes = surviving_classes, task = task)
  .yolo_wrap_dataset(out_dir_norm, task = task, classes = surviving_classes)
}

#' @export
yolo_dataset_subset <- yolo_dataset_filter


# ==============================================================================
# 6. UNIVERSAL CONVERTER: CSV / COCO / PASCAL VOC TO YOLO (yolo_dataset_convert)
# ==============================================================================
#' @title Convert External Annotations (CSV, COCO, Pascal VOC) to YOLO Format
#' @name yolo_dataset_convert
#' @description
#' `yolo_dataset_convert()` converts annotations from standard computer vision formats
#' into a clean, ready-to-train YOLO dataset.
#' Supports:
#' * **CSV**: Tabular tables (`image`, `xmin`, `ymin`, `xmax`, `ymax`, `class`) commonly exported by ImageJ, Fiji, or field records.
#' * **COCO**: Standard COCO JSON (`instances_train.json` with `images`, `annotations`, `categories`).
#' * **Pascal VOC**: Directory of XML annotation files.
#'
#' Automatically computes normalized YOLO coordinates and creates `train/` and `val/` splits.
#'
#' @param input Path to the annotation file (`.csv`, `.json`) or directory containing Pascal VOC `.xml` files.
#' @param format Format type: `"csv"` (default), `"coco"`, or `"voc"`.
#' @param img_dir Directory where the raw image files are stored.
#' @param out_dir Destination directory for the YOLO dataset. Defaults to `"yolo_converted_dataset"`.
#' @param split Proportions for train and validation splits (default: `c(train = 0.80, val = 0.20)`).
#' @param copy Logical. Whether to copy image files into the dataset. Defaults to `TRUE`.
#' @param seed Random seed for split partitioning (default: 42).
#' @param overwrite Logical. Overwrite existing `out_dir`. Defaults to `FALSE`.
#' @param verbose Logical. Show progress. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` object pointing to `out_dir`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Convert CSV bounding boxes from field phenotyping
#' ds_yolo <- yolo_dataset_convert(
#'   input   = "field_labels.csv",
#'   format  = "csv",
#'   img_dir = "field_photos",
#'   out_dir = "dataset_field_yolo"
#' )
yolo_dataset_convert <- function(input,
                                 format = c("csv", "coco", "voc"),
                                 img_dir = NULL,
                                 out_dir = "yolo_converted_dataset",
                                 split = c(train = 0.80, val = 0.20),
                                 copy = TRUE,
                                 seed = 42,
                                 overwrite = FALSE,
                                 verbose = TRUE) {
  format <- match.arg(format)
  out_dir_norm <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)

  if (dir.exists(out_dir_norm)) {
    if (isTRUE(overwrite)) unlink(out_dir_norm, recursive = TRUE, force = TRUE)
    else cli::cli_abort("Directory {.path {out_dir_norm}} already exists. Use {.code overwrite = TRUE}.")
  }

  for (s in c("train", "val")) {
    dir.create(file.path(out_dir_norm, "images", s), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir_norm, "labels", s), recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  split_norm <- split / sum(split)
  master_classes <- character()
  records_by_img <- list()

  # FORMAT 1: CSV
  if (format == "csv") {
    if (!file.exists(input)) cli::cli_abort("CSV file {.path {input}} does not exist.")
    df <- utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- tolower(trimws(names(df)))

    im_col <- intersect(c("image", "img", "filename", "file", "image_name"), names(df))[1]
    if (is.na(im_col)) cli::cli_abort("Could not find image filename column in CSV.")

    cls_col <- intersect(c("class", "label", "category", "name"), names(df))[1]
    if (is.na(cls_col)) { df$class <- "object"; cls_col <- "class" }

    has_xyxy <- all(c("xmin", "ymin", "xmax", "ymax") %in% names(df))
    has_xywh <- all(c("x", "y", "w", "h") %in% names(df)) || all(c("x_center", "y_center", "width", "height") %in% names(df))

    if (!has_xyxy && !has_xywh) {
      cli::cli_abort("CSV must contain bounding box columns: (xmin, ymin, xmax, ymax) or (x, y, w, h).")
    }

    master_classes <- unique(as.character(df[[cls_col]]))

    unique_imgs <- unique(as.character(df[[im_col]]))
    for (im_name in unique_imgs) {
      sub_d <- df[df[[im_col]] == im_name, ]
      box_lines <- character()

      full_im_path <- if (!is.null(img_dir)) file.path(img_dir, im_name) else im_name
      w_px <- 1.0; h_px <- 1.0
      is_pixel_coords <- FALSE

      if (has_xyxy) {
        if (max(sub_d$xmax, na.rm = TRUE) > 1.0 || max(sub_d$ymax, na.rm = TRUE) > 1.0) is_pixel_coords <- TRUE
      } else {
        if (max(sub_d$w, sub_d$width, na.rm = TRUE) > 1.0) is_pixel_coords <- TRUE
      }

      if (is_pixel_coords && file.exists(full_im_path)) {
        im_probe <- try(image_import(full_im_path), silent = TRUE)
        if (!inherits(im_probe, "try-error") && is_image(im_probe)) {
          w_px <- dim(im_probe)[1]
          h_px <- dim(im_probe)[2]
        }
      }

      for (r in seq_len(nrow(sub_d))) {
        cid <- which(master_classes == as.character(sub_d[[cls_col]][r])) - 1L
        if (has_xyxy) {
          x1 <- sub_d$xmin[r] / w_px; y1 <- sub_d$ymin[r] / h_px
          x2 <- sub_d$xmax[r] / w_px; y2 <- sub_d$ymax[r] / h_px
          xc <- (x1 + x2) / 2; yc <- (y1 + y2) / 2
          bw <- abs(x2 - x1); bh <- abs(y2 - y1)
        } else {
          xc <- sub_d[[intersect(c("x", "x_center"), names(sub_d))[1]]][r] / w_px
          yc <- sub_d[[intersect(c("y", "y_center"), names(sub_d))[1]]][r] / h_px
          bw <- sub_d[[intersect(c("w", "width"), names(sub_d))[1]]][r] / w_px
          bh <- sub_d[[intersect(c("h", "height"), names(sub_d))[1]]][r] / h_px
        }
        xc <- max(0.001, min(0.999, xc))
        yc <- max(0.001, min(0.999, yc))
        bw <- max(0.001, min(0.999, bw))
        bh <- max(0.001, min(0.999, bh))
        box_lines <- c(box_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, bw, bh))
      }
      records_by_img[[im_name]] <- list(full_path = full_im_path, lines = box_lines)
    }
  } else if (format == "coco") {
    # FORMAT 2: COCO JSON
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg jsonlite} is required to parse COCO JSON annotations. Install with: {.code install.packages('jsonlite')}")
    }
    if (!file.exists(input)) cli::cli_abort("COCO file {.path {input}} does not exist.")

    coco_data <- jsonlite::fromJSON(input, simplifyVector = FALSE)
    cat_map <- character()
    if (!is.null(coco_data$categories)) {
      for (cat_item in coco_data$categories) {
        cat_map[as.character(cat_item$id)] <- cat_item$name
      }
    }
    master_classes <- unname(cat_map)

    img_info_map <- list()
    if (!is.null(coco_data$images)) {
      for (im_item in coco_data$images) {
        img_info_map[[as.character(im_item$id)]] <- im_item
      }
    }

    # Group annotations by image_id
    annos_by_img <- list()
    if (!is.null(coco_data$annotations)) {
      for (ann in coco_data$annotations) {
        im_id <- as.character(ann$image_id)
        annos_by_img[[im_id]] <- c(annos_by_img[[im_id]], list(ann))
      }
    }

    for (im_id in names(img_info_map)) {
      im_meta <- img_info_map[[im_id]]
      im_fn <- im_meta$file_name
      full_im_path <- if (!is.null(img_dir)) file.path(img_dir, im_fn) else im_fn

      im_w <- as.numeric(im_meta$width)
      im_h <- as.numeric(im_meta$height)

      box_lines <- character()
      if (im_id %in% names(annos_by_img)) {
        for (ann in annos_by_img[[im_id]]) {
          cat_name <- cat_map[[as.character(ann$category_id)]]
          cid <- which(master_classes == cat_name) - 1L

          # Bbox in COCO: [xmin, ymin, width, height]
          if (!is.null(ann$bbox) && length(ann$bbox) == 4L) {
            bx <- as.numeric(ann$bbox)
            xc <- (bx[1] + bx[3] / 2) / im_w
            yc <- (bx[2] + bx[4] / 2) / im_h
            bw <- bx[3] / im_w
            bh <- bx[4] / im_h
            xc <- max(0.001, min(0.999, xc))
            yc <- max(0.001, min(0.999, yc))
            bw <- max(0.001, min(0.999, bw))
            bh <- max(0.001, min(0.999, bh))
            box_lines <- c(box_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, bw, bh))
          }
        }
      }
      records_by_img[[im_fn]] <- list(full_path = full_im_path, lines = box_lines)
    }
  } else if (format == "voc") {
    # FORMAT 3: PASCAL VOC XML
    xml_files <- if (dir.exists(input)) {
      list.files(input, pattern = "\\.xml$", full.names = TRUE, ignore.case = TRUE)
    } else if (file.exists(input)) {
      input
    } else {
      cli::cli_abort("VOC path {.path {input}} does not exist.")
    }

    if (length(xml_files) == 0L) cli::cli_abort("No XML files found in {.path {input}}.")

    for (xf in xml_files) {
      txt <- paste(readLines(xf, warn = FALSE), collapse = "\n")

      # Extract filename
      fn_m <- regmatches(txt, regexec("<filename>(.*?)</filename>", txt))[[1]]
      im_fn <- if (length(fn_m) > 1L) trimws(fn_m[2]) else sub("\\.xml$", ".jpg", basename(xf))

      # Extract width and height
      w_m <- regmatches(txt, regexec("<width>(.*?)</width>", txt))[[1]]
      h_m <- regmatches(txt, regexec("<height>(.*?)</height>", txt))[[1]]
      im_w <- if (length(w_m) > 1L) as.numeric(trimws(w_m[2])) else 1.0
      im_h <- if (length(h_m) > 1L) as.numeric(trimws(h_m[2])) else 1.0

      full_im_path <- if (!is.null(img_dir)) file.path(img_dir, im_fn) else file.path(dirname(xf), im_fn)

      # Extract objects
      obj_blocks <- strsplit(txt, "<object>")[[1]][-1]
      box_lines <- character()

      for (ob in obj_blocks) {
        nm_m <- regmatches(ob, regexec("<name>(.*?)</name>", ob))[[1]]
        if (length(nm_m) <= 1L) next
        cls_name <- trimws(nm_m[2])
        if (!(cls_name %in% master_classes)) master_classes <- c(master_classes, cls_name)
        cid <- which(master_classes == cls_name) - 1L

        xmin_m <- regmatches(ob, regexec("<xmin>(.*?)</xmin>", ob))[[1]]
        ymin_m <- regmatches(ob, regexec("<ymin>(.*?)</ymin>", ob))[[1]]
        xmax_m <- regmatches(ob, regexec("<xmax>(.*?)</xmax>", ob))[[1]]
        ymax_m <- regmatches(ob, regexec("<ymax>(.*?)</ymax>", ob))[[1]]

        if (length(xmin_m) > 1L && length(ymin_m) > 1L && length(xmax_m) > 1L && length(ymax_m) > 1L) {
          x1 <- as.numeric(xmin_m[2]); y1 <- as.numeric(ymin_m[2])
          x2 <- as.numeric(xmax_m[2]); y2 <- as.numeric(ymax_m[2])
          xc <- ((x1 + x2) / 2) / im_w
          yc <- ((y1 + y2) / 2) / im_h
          bw <- abs(x2 - x1) / im_w
          bh <- abs(y2 - y1) / im_h
          xc <- max(0.001, min(0.999, xc))
          yc <- max(0.001, min(0.999, yc))
          bw <- max(0.001, min(0.999, bw))
          bh <- max(0.001, min(0.999, bh))
          box_lines <- c(box_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, bw, bh))
        }
      }
      records_by_img[[im_fn]] <- list(full_path = full_im_path, lines = box_lines)
    }
  }

  # Partition across train / val
  img_names_shuffled <- sample(names(records_by_img))
  n_tot <- length(img_names_shuffled)
  n_tr <- max(1L, round(n_tot * split_norm[1]))
  if (n_tot >= 2L && n_tr >= n_tot) n_tr <- n_tot - 1L
  train_names <- img_names_shuffled[seq_len(n_tr)]
  val_names <- if (n_tr < n_tot) img_names_shuffled[(n_tr + 1L):n_tot] else character(0)

  for (s in c("train", "val")) {
    set_imgs <- if (s == "train") train_names else val_names
    for (im_nm in set_imgs) {
      if (is.na(im_nm) || !nzchar(im_nm)) next
      rec <- records_by_img[[im_nm]]
      if (is.null(rec)) next
      stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_nm))
      dest_im <- file.path(out_dir_norm, "images", s, basename(im_nm))
      dest_txt <- file.path(out_dir_norm, "labels", s, paste0(stem, ".txt"))

      if (!is.null(rec$full_path) && file.exists(rec$full_path)) {
        if (isTRUE(copy)) file.copy(rec$full_path, dest_im, overwrite = TRUE) else file.rename(rec$full_path, dest_im)
      }
      writeLines(rec$lines, dest_txt)
    }
  }

  .yolo_write_yaml(out_dir_norm, classes = master_classes, task = "detect")

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Successfully converted {length(records_by_img)} images ({format} format) to YOLO in {.path {out_dir_norm}}.")
  }

  .yolo_wrap_dataset(out_dir_norm, task = "detect", classes = master_classes)
}


# ==============================================================================
# 6. CONVERT DETECTION DATASET TO SEGMENTATION VIA SAM 2.1 (yolo_dataset_to_segment)
# ==============================================================================
#' @title Convert a YOLO Detection Dataset to a Segmentation Dataset using SAM 2.1
#' @name yolo_dataset_to_segment
#' @aliases yolo_dataset_detect_to_segment yolo_dataset_det_to_seg
#' @description
#' `yolo_dataset_to_segment()` transforms an existing YOLO object detection dataset
#' (with bounding box labels: `class_id cx cy w h`) into a high-precision YOLO instance
#' segmentation dataset (with polygon labels: `class_id x1 y1 x2 y2 ... xn yn`) using
#' Meta AI's Segment Anything Model 2.1 (**SAM 2.1**).
#'
#' For each image, the existing bounding boxes are extracted as spatial prompts (*box prompts*)
#' and fed to SAM 2.1, which isolates object pixels and produces sub-pixel polygon contours.
#'
#' @param data A `yolo_dataset` object or directory path to an existing YOLO detection dataset.
#' @param out_dir Output directory path for the new segmentation dataset.
#'   Defaults to `paste0(basename(dir_path), "_seg")`.
#' @param threshold Numeric. Probability threshold for SAM 2.1 mask binarization (default: `0.5`).
#' @param n_points Integer. Maximum number of polygon vertices per object (default: `40L`).
#' @param min_area Numeric. Minimum polygon pixel area to retain (default: `15`).
#' @param fill_hull Logical. Whether to fill holes inside segmented masks (default: `TRUE`).
#' @param box_margin Safety margin (padding proportion) added around bounding boxes before segmenting with Meta SAM 2.1 (e.g. 0.05 expands 5% on each side). Prevents masks from cutting off object borders. Defaults to `0.05` (5%).
#' @param copy Logical. Whether to copy image files to `out_dir` (default: `TRUE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @param device_id GPU device ID (default: `-1`).
#' @param batch_size Integer. Maximum batch size of bounding boxes passed to SAM 2.1 decoder at once (default: `64L`).
#'   Processing boxes in batches leverages GPU tensor cores for maximum throughput while staying safely within DirectML memory limits.
#' @param overwrite Logical. Whether to overwrite `out_dir` if it already exists (default: `FALSE`).
#' @param resume Logical. If `TRUE`, resumes processing an existing dataset directory by skipping images that have already been converted (default: `FALSE`).
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs the conversion as an RStudio Background Job (default: `FALSE`).
#' @param verbose Logical. Whether to display progress messages (default: `TRUE`).
#' @param ... Additional arguments. Legacy parameters such as \code{parallel} or \code{workers} are gracefully ignored with an informative notice.
#'
#' @return An object of class `c("yolo_dataset", "list")` pointing to `out_dir`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Convert a bounding-box dataset into a polygon segmentation dataset
#' ds_seg <- yolo_dataset_to_segment(
#'   data    = "dataset_cafe_det",
#'   out_dir = "dataset_cafe_seg"
#' )
yolo_dataset_to_segment <- function(data,
                                    out_dir = NULL,
                                    threshold = 0.5,
                                    n_points = 40L,
                                    min_area = 15,
                                    fill_hull = TRUE,
                                    box_margin = 0.05,
                                    copy = TRUE,
                                    engine = c("cpu", "gpu"),
                                    device_id = -1,
                                    batch_size = 64L,
                                    overwrite = FALSE,
                                    resume = FALSE,
                                    as_job = FALSE,
                                    verbose = TRUE,
                                    ...) {
  dots <- list(...)
  if ("box_margin" %in% names(dots)) box_margin <- as.numeric(dots$box_margin[1])
  if ("margin" %in% names(dots)) box_margin <- as.numeric(dots$margin[1])
  if (!is.null(dots$parallel) && isTRUE(dots$parallel)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("The {.code parallel} argument has been removed because multi-process GPU execution causes DirectML VRAM contention and freezes. Processing sequentially with batched decoder for maximum speed and stability.")
    }
  }
  engine <- match.arg(engine)

  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with segmentation conversion in current session.")
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
      job_script_path <- file.path(job_dir, paste0("yolo_seg_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_seg_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (inherits(data, "yolo_dataset")) basename(data$dir)
        else if (is.character(data) && nzchar(data[1])) basename(data[1])
        else "dataset"
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO To Segment: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_to_segment()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Detection to Segmentation - RStudio Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "ds_res <- do.call(pliman::yolo_dataset_to_segment, job_info$args)\n",
        "out_dir_target <- if (!is.null(job_info$args$out_dir)) job_info$args$out_dir else ds_res$dir\n",
        "if (!is.null(out_dir_target) && dir.exists(out_dir_target)) {\n",
        "  save_path <- file.path(out_dir_target, 'dataset_seg.rds')\n",
        "  saveRDS(ds_res, save_path)\n",
        "  cat('\\n[OK] Dataset S3 object successfully saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Detection to Segmentation Job completed successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Detection to Segmentation (pliman)")
      cli::cli_alert_success("Dataset segmentation conversion has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target dataset: {.val {target_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Directory {.path {dir_path}} does not exist.")
  }
  task_curr <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)
  if (is.null(classes) || length(classes) == 0L) {
    classes <- c("object")
  }

  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- paste0(dir_path, "_seg")
  }
  out_dir_norm <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

  if (dir.exists(out_dir_norm)) {
    if (isTRUE(overwrite)) {
      unlink(out_dir_norm, recursive = TRUE, force = TRUE)
    } else if (!isTRUE(resume)) {
      cli::cli_abort("Directory {.path {out_dir_norm}} already exists. Use {.code overwrite = TRUE} or {.code resume = TRUE}.")
    }
  }

  # Discover splits
  all_splits <- c("train", "val", "test")
  found_splits <- character()
  for (s in all_splits) {
    if (dir.exists(file.path(dir_path, "images", s)) || dir.exists(file.path(dir_path, "labels", s))) {
      found_splits <- c(found_splits, s)
    }
  }

  is_flat <- FALSE
  if (length(found_splits) == 0L) {
    if (dir.exists(file.path(dir_path, "images")) || dir.exists(file.path(dir_path, "labels"))) {
      is_flat <- TRUE
      found_splits <- c("")
    } else {
      cli::cli_abort("Could not find {.path images/} or {.path labels/} directories inside {.path {dir_path}}.")
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("Converting YOLO Detection Dataset to Segmentation via Meta SAM 2.1")
    cli::cli_ul(c(
      "Source dataset: {.path {dir_path}}",
      "Current task: {.val {task_curr}}",
      "Classes: {.val {classes}}",
      "Destination: {.path {out_dir_norm}}",
      "Engine: {.val {toupper(engine)}}",
      paste0("Batch size: ", batch_size, " (decoder)"),
      if (is.numeric(box_margin) && box_margin[1] > 0) paste0("Box margin: ", round(box_margin[1] * (if (box_margin[1] <= 1) 100 else 1)), "%") else NULL,
      "Execution: Sequential (hardware-accelerated)",
      if (isTRUE(resume)) "Resume: Enabled (skipping already converted images)" else NULL
    ))
  }

  total_converted <- 0L
  total_objects <- 0L

  for (s in found_splits) {
    if (is_flat) {
      src_img_dir <- file.path(dir_path, "images")
      src_lbl_dir <- file.path(dir_path, "labels")
      dest_img_dir <- file.path(out_dir_norm, "images")
      dest_lbl_dir <- file.path(out_dir_norm, "labels")
    } else {
      src_img_dir <- file.path(dir_path, "images", s)
      src_lbl_dir <- file.path(dir_path, "labels", s)
      dest_img_dir <- file.path(out_dir_norm, "images", s)
      dest_lbl_dir <- file.path(out_dir_norm, "labels", s)
    }

    if (!dir.exists(src_img_dir)) next
    dir.create(dest_img_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(dest_lbl_dir, recursive = TRUE, showWarnings = FALSE)

    img_files <- list.files(src_img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
    if (length(img_files) == 0L) next

    sp_label <- if (nzchar(s)) paste0(" [", s, "]") else ""

    # Worker function to process one image
    convert_single_image <- function(im_f) {
      stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_f))
      lbl_f <- file.path(src_lbl_dir, paste0(stem, ".txt"))
      dest_im_path <- file.path(dest_img_dir, basename(im_f))
      dest_lbl_path <- file.path(dest_lbl_dir, paste0(stem, ".txt"))

      # Resume support: skip if label file was already generated
      if (isTRUE(resume) && file.exists(dest_lbl_path)) {
        if (!file.exists(dest_im_path)) {
          if (isTRUE(copy)) {
            file.copy(im_f, dest_im_path, overwrite = TRUE)
          } else {
            file.rename(im_f, dest_im_path)
          }
        }
        return(list(converted = 1L, objects = 0L))
      }

      # Copy or move image
      if (isTRUE(copy)) {
        file.copy(im_f, dest_im_path, overwrite = TRUE)
      } else {
        file.rename(im_f, dest_im_path)
      }

      if (!file.exists(lbl_f)) {
        writeLines(character(0), dest_lbl_path)
        return(list(converted = 1L, objects = 0L))
      }

      lines <- readLines(lbl_f, warn = FALSE)
      lines <- lines[nzchar(trimws(lines))]
      if (length(lines) == 0L) {
        writeLines(character(0), dest_lbl_path)
        return(list(converted = 1L, objects = 0L))
      }

      src_im_to_read <- if (file.exists(im_f)) im_f else dest_im_path
      im_obj <- tryCatch(image_import(src_im_to_read), error = function(e) NULL)
      if (is.null(im_obj)) {
        writeLines(lines, dest_lbl_path)
        return(list(converted = 1L, objects = 0L))
      }
      dims <- dim(im_obj)
      img_w <- dims[1]
      img_h <- dims[2]

      boxes_list <- list()
      class_id_list <- integer()
      poly_lines_existing <- character()

      for (ln in lines) {
        vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
        if (length(vals) == 5L) {
          cid <- as.integer(vals[1])
          cx <- vals[2]; cy <- vals[3]
          bw <- vals[4]; bh <- vals[5]
          xmin <- max(1.0, (cx - bw / 2.0) * img_w)
          ymin <- max(1.0, (cy - bh / 2.0) * img_h)
          xmax <- min(as.double(img_w), (cx + bw / 2.0) * img_w)
          ymax <- min(as.double(img_h), (cy + bh / 2.0) * img_h)
          boxes_list[[length(boxes_list) + 1L]] <- c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
          class_id_list <- c(class_id_list, cid)
        } else if (length(vals) >= 6L) {
          poly_lines_existing <- c(poly_lines_existing, ln)
        }
      }

      seg_lines_out <- poly_lines_existing

      if (length(boxes_list) > 0L) {
        df_bx <- as.data.frame(do.call(rbind, boxes_list))
        sam_res <- tryCatch(
          .sam_segment_boxes(
            im_obj,
            boxes = df_bx,
            threshold = threshold,
            fill_hull = fill_hull,
            min_area = min_area,
            n_points = n_points,
            batch_size = batch_size,
            box_margin = box_margin,
            engine = engine,
            device_id = device_id,
            verbose = FALSE
          ),
          error = function(e) NULL
        )

        if (!is.null(sam_res) && length(sam_res$contours) > 0L) {
          for (k in seq_along(sam_res$contours)) {
            cnt <- sam_res$contours[[k]]
            if (is.null(cnt) || nrow(cnt) < 3L) next
            cid <- class_id_list[k]

            n_cur <- nrow(cnt)
            target_n <- min(as.integer(n_points), n_cur)
            if (n_cur > target_n) {
              sampled_idx <- round(seq(1, n_cur, length.out = target_n + 1L)[-(target_n + 1L)])
              cnt <- cnt[sampled_idx, , drop = FALSE]
            }

            x_norm <- pmin(pmax(cnt[, 1] / img_w, 0.0), 1.0)
            y_norm <- pmin(pmax(cnt[, 2] / img_h, 0.0), 1.0)
            interleaved <- as.vector(rbind(x_norm, y_norm))
            line_str <- paste(c(cid, round(interleaved, 5L)), collapse = " ")
            seg_lines_out <- c(seg_lines_out, line_str)
          }
        } else {
          for (k in seq_len(nrow(df_bx))) {
            bx <- df_bx[k, ]
            cid <- class_id_list[k]
            x1 <- pmin(pmax(bx$xmin / img_w, 0.0), 1.0); x2 <- pmin(pmax(bx$xmax / img_w, 0.0), 1.0)
            y1 <- pmin(pmax(bx$ymin / img_h, 0.0), 1.0); y2 <- pmin(pmax(bx$ymax / img_h, 0.0), 1.0)
            interleaved <- c(x1, y1, x2, y1, x2, y2, x1, y2)
            seg_lines_out <- c(seg_lines_out, paste(c(cid, round(interleaved, 5L)), collapse = " "))
          }
        }
      }

      writeLines(seg_lines_out, dest_lbl_path)
      list(converted = 1L, objects = length(seg_lines_out))
    }

    # Sequential processing with progress bar and periodic garbage collection
    pb <- if (isTRUE(verbose)) cli::cli_progress_bar(format = paste0("Processing images", sp_label, " {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [{cli::pb_percent}] | ETA: {cli::pb_eta}"), total = length(img_files)) else NULL
    img_counter <- 0L
    for (im_f in img_files) {
      res_f <- convert_single_image(im_f)
      total_converted <- total_converted + res_f$converted
      total_objects <- total_objects + res_f$objects
      if (!is.null(pb)) cli::cli_progress_update(id = pb)
      img_counter <- img_counter + 1L
      if (img_counter %% 25L == 0L) {
        gc(verbose = FALSE)
      }
    }
    if (!is.null(pb)) cli::cli_progress_done(id = pb)
  }

  .yolo_write_yaml(out_dir_norm, classes = classes, task = "segment")

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "Converted {total_converted} image(s) and {total_objects} object(s) to YOLO segmentation in {.path {out_dir_norm}}."
    )
  }

  .yolo_wrap_dataset(out_dir_norm, task = "segment", classes = classes)
}

#' @rdname yolo_dataset_to_segment
#' @export
yolo_dataset_detect_to_segment <- yolo_dataset_to_segment

#' @rdname yolo_dataset_to_segment
#' @export
yolo_dataset_det_to_seg <- yolo_dataset_to_segment


# 7. DATASET PARTITIONING & RE-SPLITTING (yolo_dataset_split)
# ==============================================================================

# Internal helper: Multi-label iterative stratification (Sechidis et al. / Szymański & Kajdanowicz)
.iterative_stratify <- function(Y, proportions) {
  N <- nrow(Y)
  C <- ncol(Y)
  K <- length(proportions)
  p <- proportions / sum(proportions)

  class_totals <- colSums(Y)
  target_per_split_class <- outer(p, class_totals) # K x C matrix
  target_imgs_per_split <- round(p * N)
  diff_imgs <- N - sum(target_imgs_per_split)
  if (diff_imgs != 0) {
    target_imgs_per_split[1] <- target_imgs_per_split[1] + diff_imgs
  }

  allocated_split <- integer(N)
  unassigned <- rep(TRUE, N)

  has_objs <- which(rowSums(Y) > 0L)
  pure_bg <- which(rowSums(Y) == 0L)

  alloc_class <- matrix(0.0, nrow = K, ncol = C)
  alloc_imgs <- integer(K)

  # Allocate images with target objects iteratively
  while (any(unassigned[has_objs])) {
    unassigned_has_objs <- has_objs[unassigned[has_objs]]

    Y_un <- Y[unassigned_has_objs, , drop = FALSE]
    col_has <- colSums(Y_un > 0L)
    pos_classes <- which(col_has > 0L)
    if (length(pos_classes) == 0L) break

    # Prioritize rarest remaining class
    min_c <- pos_classes[which.min(col_has[pos_classes])]

    # Pick image containing min_c with fewest other labels
    c_imgs <- unassigned_has_objs[Y[unassigned_has_objs, min_c] > 0L]
    img_lbl_counts <- rowSums(Y[c_imgs, , drop = FALSE] > 0L)
    best_img <- c_imgs[which.min(img_lbl_counts)]

    # Best split: prioritize splits with 0 instances if class_totals >= K, else relative deficit
    splits_with_zero <- which(alloc_class[, min_c] == 0 & p > 0)
    best_split <- if (class_totals[min_c] >= K && length(splits_with_zero) > 0L) {
      splits_with_zero[which.max(p[splits_with_zero])]
    } else {
      norm_deficits <- (target_per_split_class[, min_c] - alloc_class[, min_c]) / (target_per_split_class[, min_c] + 1e-4)
      which.max(norm_deficits)
    }

    allocated_split[best_img] <- best_split
    unassigned[best_img] <- FALSE
    alloc_imgs[best_split] <- alloc_imgs[best_split] + 1L
    alloc_class[best_split, ] <- alloc_class[best_split, ] + Y[best_img, ]
  }

  # Allocate background images proportionally
  if (length(pure_bg) > 0L) {
    for (bg_img in pure_bg) {
      rem_cap <- target_imgs_per_split - alloc_imgs
      best_split <- which.max(rem_cap)
      allocated_split[bg_img] <- best_split
      alloc_imgs[best_split] <- alloc_imgs[best_split] + 1L
    }
  }

  # Fallback for any unassigned
  remaining <- which(unassigned)
  if (length(remaining) > 0L) {
    for (ri in remaining) {
      rem_cap <- target_imgs_per_split - alloc_imgs
      best_split <- which.max(rem_cap)
      allocated_split[ri] <- best_split
      alloc_imgs[best_split] <- alloc_imgs[best_split] + 1L
    }
  }

  allocated_split
}

#' @title Partition or Re-stratify an Existing YOLO Dataset
#' @name yolo_dataset_split
#' @description
#' `yolo_dataset_split()` takes any YOLO dataset (flat or previously partitioned)
#' and redistributes all images and matching `.txt` labels into new `train/`, `val/`,
#' and optionally `test/` splits. It ensures that pairs are preserved and regenerates
#' an updated `data.yaml`.
#'
#' For object detection and segmentation datasets, setting `stratify = TRUE` (default)
#' applies multi-label iterative stratification (Sechidis et al.) to ensure that rare and
#' minority classes maintain exact representative proportions across all splits, preventing
#' zero-instance validation splits.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param out_dir Output directory path for the re-split dataset. Defaults to `"split_yolo_dataset"`.
#' @param split Proportions for partitioning (default: `c(train = 0.70, val = 0.20, test = 0.10)`).
#' @param stratify Logical. If `TRUE` (default), applies multi-label iterative stratification
#'   to ensure balanced class distributions across all splits.
#' @param seed Optional random seed for reproducible shuffling (default: 42).
#' @param copy Logical. Whether to copy files. Defaults to `TRUE`.
#' @param overwrite Logical. Overwrite existing `out_dir`. Defaults to `FALSE`.
#' @param verbose Logical. Whether to show progress. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` object pointing to the newly split dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Re-stratify into 70% train, 20% validation, 10% test with iterative stratification
#' ds_split <- yolo_dataset_split("dataset_cafe", split = c(0.7, 0.2, 0.1), stratify = TRUE)
yolo_dataset_split <- function(data,
                               out_dir = "split_yolo_dataset",
                               split = c(train = 0.70, val = 0.20, test = 0.10),
                               stratify = TRUE,
                               seed = 42,
                               copy = TRUE,
                               overwrite = FALSE,
                               verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")

  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  out_dir_norm <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)
  if (dir.exists(out_dir_norm)) {
    if (isTRUE(overwrite)) unlink(out_dir_norm, recursive = TRUE, force = TRUE)
    else cli::cli_abort("Directory {.path {out_dir_norm}} already exists. Use {.code overwrite = TRUE}.")
  }
  dir.create(out_dir_norm, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(seed)) set.seed(seed)
  split_norm <- split / sum(split)
  split_names <- if (length(split) == 2L) c("train", "val") else c("train", "val", "test")
  names(split_norm) <- split_names

  # CASE A: CLASSIFICATION RE-SPLIT
  if (task == "classify") {
    all_splits <- c("train", "val", "test")
    all_splits <- all_splits[dir.exists(file.path(dir_path, all_splits))]

    for (s in split_names) {
      for (cls in classes) {
        dir.create(file.path(out_dir_norm, s, cls), recursive = TRUE, showWarnings = FALSE)
      }
    }

    total_moved <- 0L
    for (cls in classes) {
      cls_imgs <- character()
      for (s in all_splits) {
        c_folder <- file.path(dir_path, s, cls)
        if (dir.exists(c_folder)) {
          fls <- list.files(c_folder, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
          cls_imgs <- c(cls_imgs, fls)
        }
      }
      cls_imgs <- sample(cls_imgs)
      n_tot <- length(cls_imgs)
      if (n_tot == 0L) next

      if (length(split_names) == 2L) {
        n_tr <- max(1L, round(n_tot * split_norm[1]))
        n_tr <- min(n_tr, n_tot - 1L)
        splits_map <- list(train = cls_imgs[seq_len(n_tr)], val = cls_imgs[(n_tr + 1L):n_tot])
      } else {
        n_tr <- max(1L, round(n_tot * split_norm[1]))
        n_val <- max(1L, round(n_tot * split_norm[2]))
        if (n_tr + n_val >= n_tot) { n_tr <- max(1L, n_tot - 2L); n_val <- 1L }
        splits_map <- list(
          train = cls_imgs[seq_len(n_tr)],
          val   = cls_imgs[(n_tr + 1L):(n_tr + n_val)],
          test  = cls_imgs[(n_tr + n_val + 1L):n_tot]
        )
      }

      for (sn in split_names) {
        for (im_f in splits_map[[sn]]) {
          dest_f <- file.path(out_dir_norm, sn, cls, basename(im_f))
          if (isTRUE(copy)) file.copy(im_f, dest_f, overwrite = TRUE) else file.rename(im_f, dest_f)
          total_moved <- total_moved + 1L
        }
      }
    }

    if (isTRUE(verbose)) {
      cli::cli_alert_success("Re-stratified {total_moved} classification images across {paste(split_names, collapse = ', ')} splits.")
    }
    return(.yolo_wrap_dataset(out_dir_norm, task = "classify", classes = classes))
  }

  # CASE B: DETECTION / SEGMENTATION / POSE / OBB
  all_imgs <- if (dir.exists(file.path(dir_path, "images"))) {
    list.files(file.path(dir_path, "images"), pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  } else {
    list.files(dir_path, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", recursive = FALSE, full.names = TRUE, ignore.case = TRUE)
  }
  all_lbls <- if (dir.exists(file.path(dir_path, "labels"))) {
    list.files(file.path(dir_path, "labels"), pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  } else {
    list.files(dir_path, pattern = "\\.txt$", recursive = FALSE, full.names = TRUE)
  }

  for (s in split_names) {
    dir.create(file.path(out_dir_norm, "images", s), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir_norm, "labels", s), recursive = TRUE, showWarnings = FALSE)
  }

  n_tot <- length(all_imgs)
  if (n_tot == 0L) cli::cli_abort("No image files found to split in {.path {dir_path}}.")

  # Multi-label Iterative Stratification
  splits_list <- list()
  if (isTRUE(stratify) && length(classes) > 0L) {
    # Build count matrix Y: N x C
    Y <- matrix(0L, nrow = n_tot, ncol = length(classes))
    stems_imgs <- tools::file_path_sans_ext(basename(all_imgs))
    stems_lbls <- tools::file_path_sans_ext(basename(all_lbls))
    match_indices <- match(stems_imgs, stems_lbls)

    for (i in seq_len(n_tot)) {
      lbl_idx <- match_indices[i]
      if (!is.na(lbl_idx)) {
        lns <- readLines(all_lbls[lbl_idx], warn = FALSE)
        lns <- trimws(lns)[nzchar(trimws(lns))]
        for (ln in lns) {
          vals <- suppressWarnings(as.numeric(strsplit(ln, "\\s+")[[1]]))
          if (length(vals) >= 5L) {
            cid <- as.integer(vals[1])
            if (cid + 1L <= length(classes)) {
              Y[i, cid + 1L] <- Y[i, cid + 1L] + 1L
            }
          }
        }
      }
    }

    # Run iterative stratification
    split_assign <- .iterative_stratify(Y, split_norm)
    for (si in seq_along(split_names)) {
      splits_list[[split_names[si]]] <- all_imgs[split_assign == si]
    }
  } else {
    # Simple random shuffling fallback
    shuffled_imgs <- sample(all_imgs)
    if (length(split_names) == 2L) {
      n_tr <- max(1L, round(n_tot * split_norm[1]))
      n_tr <- min(n_tr, n_tot - 1L)
      splits_list <- list(train = shuffled_imgs[seq_len(n_tr)], val = shuffled_imgs[(n_tr + 1L):n_tot])
    } else {
      n_tr <- max(1L, round(n_tot * split_norm[1]))
      n_val <- max(1L, round(n_tot * split_norm[2]))
      if (n_tr + n_val >= n_tot) { n_tr <- max(1L, n_tot - 2L); n_val <- 1L }
      splits_list <- list(
        train = shuffled_imgs[seq_len(n_tr)],
        val   = shuffled_imgs[(n_tr + 1L):(n_tr + n_val)],
        test  = shuffled_imgs[(n_tr + n_val + 1L):n_tot]
      )
    }
  }

  for (s in split_names) {
    cur_imgs <- splits_list[[s]]
    for (im_f in cur_imgs) {
      stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(im_f))
      dest_im <- file.path(out_dir_norm, "images", s, basename(im_f))
      dest_lbl <- file.path(out_dir_norm, "labels", s, paste0(stem, ".txt"))

      if (isTRUE(copy)) file.copy(im_f, dest_im, overwrite = TRUE) else file.rename(im_f, dest_im)

      match_lbl <- all_lbls[sub("\\.txt$", "", basename(all_lbls)) == stem]
      if (length(match_lbl) > 0L && file.exists(match_lbl[1])) {
        if (isTRUE(copy)) file.copy(match_lbl[1], dest_lbl, overwrite = TRUE) else file.rename(match_lbl[1], dest_lbl)
      } else {
        file.create(dest_lbl)
      }
    }
  }

  .yolo_write_yaml(out_dir_norm, classes = classes, task = task)
  if (identical(task, "pose")) {
    orig_yaml <- file.path(dir_path, "data.yaml")
    if (file.exists(orig_yaml)) {
      y_orig <- readLines(orig_yaml, warn = FALSE)
      pose_extra <- character(0)
      k_idx <- grep("^\\s*kpt_shape\\s*:", y_orig)
      if (length(k_idx) > 0L) pose_extra <- c(pose_extra, y_orig[k_idx[1]])
      s_idx <- grep("^\\s*skeleton\\s*:", y_orig)
      if (length(s_idx) > 0L) {
        pose_extra <- c(pose_extra, y_orig[s_idx[1]])
        es <- s_idx[1] + 1L
        while (es <= length(y_orig) && grepl("^\\s*-\\s*\\[", y_orig[es])) {
          pose_extra <- c(pose_extra, y_orig[es])
          es <- es + 1L
        }
      }
      if (length(pose_extra) > 0L) {
        dest_yaml <- file.path(out_dir_norm, "data.yaml")
        dest_lines <- readLines(dest_yaml, warn = FALSE)
        writeLines(c(dest_lines, "", pose_extra), dest_yaml)
      }
    }
  }

  if (isTRUE(verbose)) {
    strat_text <- if (isTRUE(stratify)) " (multi-label iterative stratification)" else " (random)"
    cli::cli_alert_success("Re-split {n_tot} image and label pairs across {paste(split_names, collapse = ', ')} splits{strat_text}.")
  }

  .yolo_wrap_dataset(out_dir_norm, task = task, classes = classes)
}


# ==============================================================================
# 8. INSTANCE COPY-PASTE AUGMENTATION & CLASS BALANCING
# ==============================================================================

# Internal helper: Separable box filter for alpha feathering
.feather_mask <- function(m, radius = 2L) {
  if (radius <= 0L) return(m)
  w <- nrow(m); h <- ncol(m)
  temp <- matrix(0.0, w, h)
  for (y in seq_len(h)) {
    row_cum <- c(0, cumsum(m[, y]))
    for (x in seq_len(w)) {
      x1 <- max(1L, x - radius); x2 <- min(w, x + radius)
      temp[x, y] <- (row_cum[x2 + 1] - row_cum[x1]) / (x2 - x1 + 1)
    }
  }
  out <- matrix(0.0, w, h)
  for (x in seq_len(w)) {
    col_cum <- c(0, cumsum(temp[x, ]))
    for (y in seq_len(h)) {
      y1 <- max(1L, y - radius); y2 <- min(h, y + radius)
      out[x, y] <- (col_cum[y2 + 1] - col_cum[y1]) / (y2 - y1 + 1)
    }
  }
  out
}

# Internal helper: Build instance bank from training images and labels
.build_instance_bank <- function(dir_path, target_classes = NULL, blend_edge = 3L) {
  img_dir <- file.path(dir_path, "images", "train")
  if (!dir.exists(img_dir)) img_dir <- file.path(dir_path, "train", "images")
  if (!dir.exists(img_dir)) img_dir <- file.path(dir_path, "train")
  lbl_dir <- file.path(dir_path, "labels", "train")
  if (!dir.exists(lbl_dir)) lbl_dir <- file.path(dir_path, "train", "labels")
  if (!dir.exists(lbl_dir)) lbl_dir <- file.path(dir_path, "train")

  classes <- .yolo_get_classes(dir_path)
  fls <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  stems <- tools::file_path_sans_ext(basename(fls))
  lbl_files <- file.path(lbl_dir, paste0(stems, ".txt"))

  bank <- list()
  for (cls in classes) bank[[cls]] <- list()

  target_cids <- if (!is.null(target_classes)) {
    if (is.character(target_classes)) which(classes %in% target_classes) - 1L else as.integer(target_classes)
  } else {
    seq_along(classes) - 1L
  }

  for (i in seq_along(fls)) {
    if (!file.exists(lbl_files[i])) next
    lines <- readLines(lbl_files[i], warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0L) next

    parsed_lines <- lapply(lines, function(ln) suppressWarnings(as.numeric(strsplit(ln, "\\s+")[[1]])))
    cids_in_img <- vapply(parsed_lines, function(v) if (length(v) >= 5) as.integer(v[1]) else -1L, integer(1))
    matching_objs <- which(cids_in_img %in% target_cids)
    if (length(matching_objs) == 0L) next

    im <- try(image_import(fls[i]), silent = TRUE)
    if (inherits(im, "try-error") || !is_image(im)) next

    W <- dim(im)[1]; H <- dim(im)[2]
    u <- unclass(im)
    arr <- if (is.raw(u)) as.numeric(u) / 255.0 else as.numeric(u)
    dim(arr) <- dim(u)

    for (obj_idx in matching_objs) {
      vals <- parsed_lines[[obj_idx]]
      cid <- as.integer(vals[1])
      cls_name <- if (cid + 1L <= length(classes)) classes[cid + 1L] else as.character(cid)

      is_poly <- (length(vals) > 5L && (length(vals) - 1L) %% 2L == 0L)

      if (is_poly) {
        xs <- vals[c(TRUE, FALSE)][-1]
        ys <- vals[c(FALSE, TRUE)]
        px_x <- xs * W; px_y <- ys * H
        xmin <- max(1L, floor(min(px_x))); xmax <- min(W, ceiling(max(px_x)))
        ymin <- max(1L, floor(min(px_y))); ymax <- min(H, ceiling(max(px_y)))
        pw <- xmax - xmin + 1L; ph <- ymax - ymin + 1L
        if (pw < 4L || ph < 4L) next

        patch <- arr[xmin:xmax, ymin:ymax, , drop = FALSE]
        rel_x <- px_x - xmin + 1L; rel_y <- px_y - ymin + 1L
        mask_bin <- try(polygon_to_binary(cbind(rel_x, rel_y)), silent = TRUE)
        if (inherits(mask_bin, "try-error") || length(mask_bin) == 0L) {
          mask_bin <- matrix(TRUE, pw, ph)
        }
        alpha <- .feather_mask(mask_bin * 1.0, radius = blend_edge)
        dim(alpha) <- c(pw, ph)

        bank[[cls_name]][[length(bank[[cls_name]]) + 1L]] <- list(
          cid = cid,
          cls_name = cls_name,
          patch = patch,
          alpha = alpha,
          is_poly = TRUE,
          rel_poly_norm = cbind((rel_x - 1L) / pw, (rel_y - 1L) / ph)
        )
      } else if (length(vals) == 5L) {
        xc <- vals[2]; yc <- vals[3]; bw <- vals[4]; bh <- vals[5]
        xmin <- max(1L, floor((xc - bw / 2) * W))
        xmax <- min(W, ceiling((xc + bw / 2) * W))
        ymin <- max(1L, floor((yc - bh / 2) * H))
        ymax <- min(H, ceiling((yc + bh / 2) * H))
        pw <- xmax - xmin + 1L; ph <- ymax - ymin + 1L
        if (pw < 4L || ph < 4L) next

        patch <- arr[xmin:xmax, ymin:ymax, , drop = FALSE]
        mask_rect <- matrix(1.0, pw, ph)
        alpha <- .feather_mask(mask_rect, radius = blend_edge)
        dim(alpha) <- c(pw, ph)

        bank[[cls_name]][[length(bank[[cls_name]]) + 1L]] <- list(
          cid = cid,
          cls_name = cls_name,
          patch = patch,
          alpha = alpha,
          is_poly = FALSE,
          bw = bw,
          bh = bh
        )
      }
    }
  }

  bank
}

#' @title Simple Copy-Paste Instance Augmentation for YOLO Datasets
#' @name yolo_dataset_copypaste
#' @aliases yolo_dataset_paste
#' @description
#' `yolo_dataset_copypaste()` implements native instance-level Copy-Paste augmentation
#' (inspired by Ghiasi et al., CVPR) for YOLO object detection and instance segmentation datasets.
#'
#' It extracts real instance objects (bounding boxes or exact polygon masks) of target or
#' minority classes into an in-memory Object Bank, and intelligently pastes them onto other
#' training images with:
#' * Scale jitter (\code{scale_range})
#' * Random horizontal/vertical reflection (\code{flip_h}, \code{flip_v})
#' * Ambient luminance and contrast adaptation (\code{match_luminance})
#' * Smooth alpha-feathered edge blending (\code{blend_edge})
#' * Geometric occlusion checking and annotation updating (\code{max_iou})
#'
#' @param data A `yolo_dataset` object, directory path, or YAML file path.
#' @param out_dir Character. Directory to save the augmented dataset. If `NULL` (default),
#'   augments in-place using suffix `_copypaste<copy>`.
#' @param classes Character vector or integer IDs of target classes to extract and paste.
#'   If `NULL` (default), automatically identifies minority classes for balanced injection.
#' @param target_count Optional integer target instances per class in `train/`. If specified,
#'   continues copy-pasting until each target class achieves at least `target_count` instances.
#' @param copies Integer. Number of copy-paste augmented images generated per training image (default: `1L`).
#' @param total_imgs Optional integer. Target total number of training images after augmentation.
#' @param n_pastes Integer vector of length 1 or 2. Range of instances pasted per augmented image (default: `c(1, 3)`).
#' @param scale_range Numeric vector of length 2. Scaling factor range for pasted instances (default: `c(0.75, 1.25)`).
#' @param rotation Numeric vector of length 1 or 2. Rotation angles in degrees for pasted instances (default: `c(-20, 20)`).
#' @param flip_h Logical. Whether to apply 50% random horizontal reflection to pasted instances (default: `TRUE`).
#' @param flip_v Logical. Whether to apply random vertical reflection to pasted instances (default: `FALSE`).
#' @param match_luminance Logical. Whether to adapt local mean luminance and contrast of pasted instances
#'   to match the surrounding target image background (default: `TRUE`).
#' @param match_strength Numeric between 0 and 1. Intensity of luminance matching (default: `0.40`).
#' @param blend_edge Integer. Pixel radius for alpha-feathering instance edges (default: `3L`).
#' @param max_iou Numeric between 0 and 1. Maximum allowable occlusion ratio before removing
#'   underlying annotations (default: `0.70`).
#' @param quality Integer (1-100). JPEG compression quality for exported images (default: `90`).
#' @param clean Logical. If `TRUE`, cleans previously generated `_copypaste` augmented images and labels before generating new copies. Defaults to `FALSE`.
#' @param seed Optional integer seed for reproducibility (default: `42`).
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs as a background job (default: `FALSE`).
#' @param verbose Logical. Whether to display progress messages (default: `TRUE`).
#' @param progress Logical. Whether to show a CLI progress bar (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the augmented dataset directory.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Copy-paste minority instances onto training images
#' ds_cp <- yolo_dataset_copypaste("my_dataset")
#'
#' # 2. Specify target classes and higher paste density
#' ds_cp2 <- yolo_dataset_copypaste("my_dataset",
#'                                  classes = "rare_defect",
#'                                  n_pastes = c(2, 5),
#'                                  scale_range = c(0.8, 1.3))
#'
#' # 3. Revert and clean up copy-pasted files
#' yolo_dataset_copypaste_clean("my_dataset")
yolo_dataset_copypaste <- function(data,
                                   out_dir = NULL,
                                   classes = NULL,
                                   target_count = NULL,
                                   copies = 1L,
                                   total_imgs = NULL,
                                   n_pastes = c(1L, 3L),
                                   scale_range = c(0.75, 1.25),
                                   rotation = c(-20, 20),
                                   flip_h = TRUE,
                                   flip_v = FALSE,
                                   match_luminance = TRUE,
                                   match_strength = 0.40,
                                   blend_edge = 3L,
                                   max_iou = 0.70,
                                   quality = 90L,
                                   clean = FALSE,
                                   seed = 42,
                                   as_job = FALSE,
                                   verbose = TRUE,
                                   progress = TRUE) {
  if (isTRUE(clean)) {
    yolo_dataset_copypaste_clean(data, verbose = FALSE)
  }
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding in current session.")
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
      job_script_path <- file.path(job_dir, paste0("yolo_cp_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_cp_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        d_p <- .yolo_get_dir(data)
        basename(d_p)
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Copy-Paste: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_copypaste()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Copy-Paste Augmentation - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_copypaste, job_info$args)\n",
        "target_dir <- if (!is.null(res$dir)) res$dir else job_info$args$out_dir\n",
        "if (!is.null(target_dir) && dir.exists(target_dir)) {\n",
        "  save_path <- file.path(target_dir, 'dataset_copypasted.rds')\n",
        "  saveRDS(res, save_path)\n",
        "  cat('\\n[OK] Copy-pasted dataset object saved to:\\n  ', save_path, '\\n')\n",
        "}\n",
        "cat('\\n[OK] YOLO Copy-Paste Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )
      if (isTRUE(verbose)) {
        cli::cli_alert_success("Launched YOLO Copy-Paste augmentation as RStudio Background Job: {.val {j_name}}")
      }
      return(invisible(job_id))
    }
  }

  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")

  task <- .yolo_detect_task(dir_path)
  all_classes <- .yolo_get_classes(dir_path)
  if (!is.null(seed)) set.seed(seed)

  if (task == "classify") {
    cli::cli_alert_warning("Copy-Paste augmentation is designed for object detection and instance segmentation datasets.")
    return(.yolo_wrap_dataset(dir_path, task = task, classes = all_classes))
  }

  # Target directory handling (in-place vs clone)
  target_dir <- if (!is.null(out_dir)) normalizePath(out_dir[1], winslash = "/", mustWork = FALSE) else dir_path
  if (!is.null(out_dir) && target_dir != dir_path) {
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    subdirs <- c("images/train", "images/val", "images/test", "labels/train", "labels/val", "labels/test")
    for (sd in subdirs) {
      src_sd <- file.path(dir_path, sd)
      dst_sd <- file.path(target_dir, sd)
      if (dir.exists(src_sd)) {
        dir.create(dst_sd, recursive = TRUE, showWarnings = FALSE)
        fls <- list.files(src_sd, full.names = TRUE)
        if (length(fls) > 0L) file.copy(fls, dst_sd, overwrite = TRUE)
      }
    }
    y_src <- file.path(dir_path, "data.yaml")
    if (file.exists(y_src)) file.copy(y_src, file.path(target_dir, "data.yaml"), overwrite = TRUE)
  }

  # Locate train directories
  tr_img_dir <- file.path(target_dir, "images", "train")
  if (!dir.exists(tr_img_dir)) tr_img_dir <- file.path(target_dir, "train", "images")
  if (!dir.exists(tr_img_dir)) tr_img_dir <- file.path(target_dir, "train")
  tr_lbl_dir <- file.path(target_dir, "labels", "train")
  if (!dir.exists(tr_lbl_dir)) tr_lbl_dir <- file.path(target_dir, "train", "labels")
  if (!dir.exists(tr_lbl_dir)) tr_lbl_dir <- file.path(target_dir, "train")

  if (!dir.exists(tr_img_dir)) cli::cli_abort("Could not find training images directory in {.path {target_dir}}.")

  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  train_imgs <- list.files(tr_img_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
  train_imgs <- train_imgs[!grepl("(_copypaste\\d+|_bal\\d+)", basename(train_imgs))]

  if (length(train_imgs) == 0L) {
    cli::cli_alert_warning("No pristine training images found to augment.")
    return(.yolo_wrap_dataset(target_dir, task = task, classes = all_classes))
  }

  # Extract instance bank
  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Copy-Paste Instance Augmentation (pliman)")
    cli::cli_alert_info("Extracting instance objects from training annotations...")
  }

  bank <- .build_instance_bank(target_dir, target_classes = classes, blend_edge = blend_edge)
  bank_counts <- vapply(bank, length, integer(1))
  active_classes <- names(bank)[bank_counts > 0L]

  if (length(active_classes) == 0L) {
    cli::cli_alert_warning("No valid instance objects found for the specified classes in the training set.")
    return(.yolo_wrap_dataset(target_dir, task = task, classes = all_classes))
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Instance Bank built with {sum(bank_counts)} total instances across {length(active_classes)} classes:")
    for (ac in active_classes) {
      cli::cli_bullets(c(" " = paste0("{.val ", ac, "}: ", bank_counts[ac], " instances")))
    }
  }

  # Determine target images and copies
  n_target <- length(train_imgs)
  n_cp_vec <- if (!is.null(total_imgs)) {
    tot <- as.integer(total_imgs)
    delta <- max(0L, tot - n_target)
    if (delta == 0L) rep(0L, n_target) else {
      base_cp <- delta %/% n_target; rem <- delta %% n_target
      v <- rep(base_cp, n_target)
      if (rem > 0L) v[seq_len(rem)] <- v[seq_len(rem)] + 1L
      v
    }
  } else {
    rep(as.integer(max(1L, copies)), n_target)
  }

  total_to_create <- sum(n_cp_vec)
  if (total_to_create == 0L) {
    cli::cli_alert_info("Target image count already satisfied. No new images created.")
    return(.yolo_wrap_dataset(target_dir, task = task, classes = all_classes))
  }

  pb <- if (isTRUE(verbose) && isTRUE(progress)) {
    cli::cli_progress_bar("Copy-pasting instances", total = total_to_create)
  } else NULL

  n_pastes_range <- if (length(n_pastes) == 1L) c(n_pastes, n_pastes) else n_pastes
  scale_min <- min(scale_range); scale_max <- max(scale_range)

  total_created <- 0L

  for (i in seq_along(train_imgs)) {
    n_cp <- n_cp_vec[i]
    if (n_cp == 0L) next

    src_img <- train_imgs[i]
    stem <- tools::file_path_sans_ext(basename(src_img))
    ext <- tools::file_ext(src_img)
    src_lbl <- file.path(tr_lbl_dir, paste0(stem, ".txt"))

    im_dest <- try(image_import(src_img), silent = TRUE)
    if (inherits(im_dest, "try-error") || !is_image(im_dest)) next

    W_dest <- dim(im_dest)[1]; H_dest <- dim(im_dest)[2]
    u_dest <- unclass(im_dest)
    base_arr <- if (is.raw(u_dest)) as.numeric(u_dest) / 255.0 else as.numeric(u_dest)
    dim(base_arr) <- dim(u_dest)

    orig_lines <- if (file.exists(src_lbl)) {
      lns <- readLines(src_lbl, warn = FALSE)
      trimws(lns)[nzchar(trimws(lns))]
    } else character(0)

    for (cp in seq_len(n_cp)) {
      cur_arr <- base_arr
      cur_lines <- orig_lines

      k_pastes <- sample(n_pastes_range[1]:n_pastes_range[2], 1L)

      for (kp in seq_len(k_pastes)) {
        # Pick class (prioritize minority)
        chosen_cls <- sample(active_classes, 1L)
        cls_instances <- bank[[chosen_cls]]
        if (length(cls_instances) == 0L) next

        inst <- sample(cls_instances, 1L)[[1]]
        s <- runif(1L, scale_min, scale_max)

        pw <- nrow(inst$patch); ph <- ncol(inst$patch)
        nw <- max(4L, round(pw * s)); nh <- max(4L, round(ph * s))
        if (nw >= W_dest - 10L) nw <- W_dest - 10L
        if (nh >= H_dest - 10L) nh <- H_dest - 10L
        if (nw < 4L || nh < 4L) next

        scaled_patch <- image_resize_cpp(inst$patch, nw, nh, 0L)
        scaled_alpha <- image_resize_cpp(inst$alpha, nw, nh, 0L)
        dim(scaled_alpha) <- c(nw, nh)

        # Flips
        poly_rel <- if (inst$is_poly) inst$rel_poly_norm else NULL
        if (isTRUE(flip_h) && runif(1) > 0.5) {
          scaled_patch <- scaled_patch[rev(seq_len(nw)), , , drop = FALSE]
          scaled_alpha <- scaled_alpha[rev(seq_len(nw)), , drop = FALSE]
          if (!is.null(poly_rel)) poly_rel[, 1] <- 1.0 - poly_rel[, 1]
        }
        if (isTRUE(flip_v) && runif(1) > 0.5) {
          scaled_patch <- scaled_patch[, rev(seq_len(nh)), , drop = FALSE]
          scaled_alpha <- scaled_alpha[, rev(seq_len(nh)), drop = FALSE]
          if (!is.null(poly_rel)) poly_rel[, 2] <- 1.0 - poly_rel[, 2]
        }

        # Random placement inside canvas
        max_x0 <- W_dest - nw
        max_y0 <- H_dest - nh
        if (max_x0 < 1L || max_y0 < 1L) next
        px0 <- sample(seq_len(max_x0), 1L)
        py0 <- sample(seq_len(max_y0), 1L)

        dest_sub <- cur_arr[px0:(px0 + nw - 1L), py0:(py0 + nh - 1L), , drop = FALSE]

        # Luminance & contrast matching
        if (isTRUE(match_luminance)) {
          fg_idx <- which(scaled_alpha > 0.3)
          if (length(fg_idx) > 10L) {
            n_ch <- dim(scaled_patch)[3]
            for (ch in seq_len(n_ch)) {
              m_d <- mean(dest_sub[, , ch][fg_idx], na.rm = TRUE)
              s_d <- sd(dest_sub[, , ch][fg_idx], na.rm = TRUE)
              m_s <- mean(scaled_patch[, , ch][fg_idx], na.rm = TRUE)
              s_s <- sd(scaled_patch[, , ch][fg_idx], na.rm = TRUE)
              if (!is.na(s_s) && s_s > 0.001 && !is.na(s_d)) {
                sc <- min(1.8, max(0.6, s_d / s_s))
                matched <- (scaled_patch[, , ch] - m_s) * sc + m_d
                matched <- pmax(0.0, pmin(1.0, matched))
                scaled_patch[, , ch] <- (1.0 - match_strength) * scaled_patch[, , ch] + match_strength * matched
              }
            }
          }
        }

        # Alpha blending
        n_ch <- dim(scaled_patch)[3]
        for (ch in seq_len(n_ch)) {
          cur_arr[px0:(px0 + nw - 1L), py0:(py0 + nh - 1L), ch] <-
            dest_sub[, , ch] * (1.0 - scaled_alpha) + scaled_patch[, , ch] * scaled_alpha
        }

        # Occlusion check on existing annotations
        new_bx1 <- px0 / W_dest; new_bx2 <- (px0 + nw) / W_dest
        new_by1 <- py0 / H_dest; new_by2 <- (py0 + nh) / H_dest

        surviving_lines <- character(0)
        for (ln in cur_lines) {
          vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
          if (length(vals) < 5L) next
          if (length(vals) == 5L) {
            ox1 <- vals[2] - vals[4] / 2; ox2 <- vals[2] + vals[4] / 2
            oy1 <- vals[3] - vals[5] / 2; oy2 <- vals[3] + vals[5] / 2
            oarea <- vals[4] * vals[5]
          } else {
            oxs <- vals[c(TRUE, FALSE)][-1]; oys <- vals[c(FALSE, TRUE)]
            ox1 <- min(oxs); ox2 <- max(oxs); oy1 <- min(oys); oy2 <- max(oys)
            oarea <- (ox2 - ox1) * (oy2 - oy1)
          }

          ix1 <- max(ox1, new_bx1); ix2 <- min(ox2, new_bx2)
          iy1 <- max(oy1, new_by1); iy2 <- min(oy2, new_by2)
          if (ix2 > ix1 && iy2 > iy1) {
            iarea <- (ix2 - ix1) * (iy2 - iy1)
            ioa <- iarea / (oarea + 1e-6)
            if (ioa >= max_iou) next # occluded by pasted object
          }
          surviving_lines <- c(surviving_lines, ln)
        }

        # Add newly pasted object annotation
        if (inst$is_poly) {
          poly_px_x <- px0 + poly_rel[, 1] * nw
          poly_px_y <- py0 + poly_rel[, 2] * nh
          gx <- pmax(0.0, pmin(1.0, poly_px_x / W_dest))
          gy <- pmax(0.0, pmin(1.0, poly_px_y / H_dest))
          c_str <- sprintf("%.6f %.6f", gx, gy)
          surviving_lines <- c(surviving_lines, paste(c(inst$cid, c_str), collapse = " "))
        } else {
          n_xc <- (px0 + nw / 2.0) / W_dest
          n_yc <- (py0 + nh / 2.0) / H_dest
          n_bw <- nw / W_dest
          n_bh <- nh / H_dest
          surviving_lines <- c(surviving_lines, sprintf("%d %.6f %.6f %.6f %.6f", inst$cid, n_xc, n_yc, n_bw, n_bh))
        }

        cur_lines <- surviving_lines
      }

      # Export augmented image and label file
      new_fn <- paste0(stem, "_copypaste", cp, ".", ext)
      new_txt <- paste0(stem, "_copypaste", cp, ".txt")

      write_image_cpp(cur_arr, file.path(tr_img_dir, new_fn), quality = quality)
      writeLines(cur_lines, file.path(tr_lbl_dir, new_txt))
      total_created <- total_created + 1L

      if (!is.null(pb)) cli::cli_progress_update(id = pb)
    }
  }

  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  # Invalidate stale cache files
  unlink(file.path(tr_lbl_dir, "train.cache"), force = TRUE)
  unlink(file.path(target_dir, "train.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Generated {total_created} Copy-Paste augmented image(s) in {.path {tr_img_dir}}.")
  }

  .yolo_wrap_dataset(target_dir, task = task, classes = all_classes)
}

#' @rdname yolo_dataset_copypaste
#' @export
yolo_dataset_paste <- yolo_dataset_copypaste

#' @title Clean Copy-Paste Augmented Images and Labels from a YOLO Dataset
#' @name yolo_dataset_copypaste_clean
#' @aliases yolo_dataset_clean_copypaste
#' @description
#' `yolo_dataset_copypaste_clean()` removes all synthetic copy-paste augmented image and label
#' files (matching pattern `_copypaste\\d+`) from the dataset, reverting it to its pristine state.
#'
#' @param data A `yolo_dataset` object or directory path.
#' @param verbose Logical. Whether to print progress messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
yolo_dataset_copypaste_clean <- function(data, verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  all_imgs <- list.files(dir_path, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  cp_imgs <- all_imgs[grepl("_copypaste\\d+", basename(all_imgs))]

  n_removed_img <- 0L
  if (length(cp_imgs) > 0L) {
    unlink(cp_imgs, force = TRUE)
    n_removed_img <- length(cp_imgs)
  }

  all_lbls <- list.files(dir_path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  cp_lbls <- all_lbls[grepl("_copypaste\\d+", basename(all_lbls))]

  n_removed_lbl <- 0L
  if (length(cp_lbls) > 0L) {
    unlink(cp_lbls, force = TRUE)
    n_removed_lbl <- length(cp_lbls)
  }

  # Invalidate caches
  unlink(list.files(dir_path, pattern = "\\.cache$", full.names = TRUE, recursive = TRUE), force = TRUE)

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Copy-Paste Clean")
    cli::cli_alert_success("Removed {n_removed_img} augmented image(s) and {n_removed_lbl} label file(s).")
  }

  .yolo_wrap_dataset(dir_path, task = task, classes = classes)
}

#' @rdname yolo_dataset_copypaste_clean
#' @export
yolo_dataset_clean_copypaste <- yolo_dataset_copypaste_clean


# ==============================================================================
# 8b. CLASS BALANCING (yolo_dataset_balance)
# ==============================================================================
#' @title Balance Class Distributions via Copy-Paste, Augmentation or Downsampling
#' @name yolo_dataset_balance
#' @description
#' `yolo_dataset_balance()` balances class distributions across all YOLO dataset types
#' (classification, object detection, instance segmentation, and OBB) to eliminate majority class bias.
#'
#' * For **object detection and instance segmentation** (`task %in% c("detect", "segment", "obb")`):
#'   * `strategy = "copypaste"` (default): Uses state-of-the-art instance-level Copy-Paste augmentation
#'     to extract instances of minority classes and blend them into training images until all classes
#'     reach `target_count`.
#'   * `strategy = "undersample"`: Downsamples images containing exclusively majority classes.
#' * For **image classification** (`task == "classify"`):
#'   * `strategy = "oversample"`: Balances folder images via geometric/photometric augmentation.
#'   * `strategy = "undersample"`: Downsamples majority class folders to `target_count`.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param strategy Strategy to apply: `"copypaste"` (default for detection/segmentation),
#'   `"oversample"`, or `"undersample"`.
#' @param target_count Optional integer target count per class in `train/`. If `NULL`,
#'   uses the maximum class count for oversampling/copy-paste and the minimum count for undersampling.
#' @param classes Optional character vector or integer IDs of specific classes to balance.
#' @param clean Logical. If `TRUE` (default), cleans previously balanced files (`_bal*.ext`, `_copypaste*.ext`)
#'   before computing class distributions.
#' @param seed Optional random seed (default: `42`).
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs as a background job (default: `FALSE`).
#' @param verbose Logical. Whether to show progress messages (default: `TRUE`).
#' @param ... Additional arguments forwarded to [yolo_dataset_copypaste()].
#'
#' @return A `yolo_dataset` object with balanced class distributions.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Balance detection/segmentation dataset using instance Copy-Paste
#' ds_balanced <- yolo_dataset_balance("my_yolo_dataset")
yolo_dataset_balance <- function(data,
                                 strategy = c("copypaste", "oversample", "undersample"),
                                 target_count = NULL,
                                 classes = NULL,
                                 clean = TRUE,
                                 seed = 42,
                                 as_job = FALSE,
                                 verbose = TRUE,
                                 ...) {
  strategy <- match.arg(strategy)
  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")
  task <- .yolo_detect_task(dir_path)
  all_classes <- .yolo_get_classes(dir_path)

  if (!is.null(seed)) set.seed(seed)

  # Clean previous balancing artifacts if requested
  if (isTRUE(clean)) {
    old_bal <- list.files(dir_path, pattern = "(_bal\\d+|_copypaste\\d+)", full.names = TRUE, recursive = TRUE)
    if (length(old_bal) > 0L) {
      unlink(old_bal, force = TRUE)
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Cleaned {length(old_bal)} previously balanced image/label file(s).")
      }
    }
  }

  # ----------------------------------------------------------------------------
  # CASE 1: Classification datasets
  # ----------------------------------------------------------------------------
  if (task == "classify") {
    if (strategy == "copypaste") strategy <- "oversample"
    tr_dir <- file.path(dir_path, "train")
    if (!dir.exists(tr_dir)) cli::cli_abort("Missing {.path train/} directory.")

    counts <- integer()
    files_by_cls <- list()
    for (cls in all_classes) {
      c_folder <- file.path(tr_dir, cls)
      fls <- list.files(c_folder, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
      fls <- fls[!grepl("_bal\\d+", basename(fls))]
      files_by_cls[[cls]] <- fls
      counts[cls] <- length(fls)
    }

    if (is.null(target_count)) {
      target_count <- if (strategy == "oversample") max(counts) else min(counts)
    } else {
      target_count <- as.integer(target_count)
    }

    if (isTRUE(verbose)) {
      cli::cli_h2("YOLO Classification Class Balancing (pliman)")
      cli::cli_alert_info("Strategy: {.val {strategy}} | Target count per class: {.val {target_count}}")
    }

    for (cls in all_classes) {
      fls <- files_by_cls[[cls]]
      curr_n <- length(fls)
      c_folder <- file.path(tr_dir, cls)

      if (strategy == "undersample" && curr_n > target_count) {
        to_remove <- sample(fls, size = curr_n - target_count)
        unlink(to_remove)
      } else if (strategy == "oversample" && curr_n < target_count && curr_n > 0L) {
        needed <- target_count - curr_n
        picks <- sample(fls, size = needed, replace = TRUE)
        for (i in seq_along(picks)) {
          src_f <- picks[i]
          im <- try(image_import(src_f), silent = TRUE)
          if (inherits(im, "try-error") || !is_image(im)) next

          im_aug <- im
          if (runif(1) > 0.5) im_aug <- image_hreflect(im_aug, verbose = FALSE)
          if (runif(1) > 0.5) im_aug <- image_vreflect(im_aug, verbose = FALSE)
          ang <- sample(c(0, 90, 180, 270), 1)
          if (ang != 0) im_aug <- image_rotate(im_aug, angle = ang, verbose = FALSE)

          arr <- image_data(im_aug, type = "normalized")
          arr <- arr * runif(1, 0.90, 1.10)
          arr[arr > 1] <- 1; arr[arr < 0] <- 0

          ext <- tools::file_ext(src_f)
          stem <- sub("\\.[a-zA-Z0-9]+$", "", basename(src_f))
          new_name <- paste0(stem, "_bal", i, ".", ext)
          image_export(as_image(arr, colormode = "Color"), file.path(c_folder, new_name), verbose = FALSE)
        }
      }
    }

    if (isTRUE(verbose)) cli::cli_alert_success("Classification class balancing completed successfully!")
    return(.yolo_wrap_dataset(dir_path, task = "classify", classes = all_classes))
  }

  # ----------------------------------------------------------------------------
  # CASE 2: Object Detection, Instance Segmentation & OBB
  # ----------------------------------------------------------------------------
  tr_lbl_dir <- file.path(dir_path, "labels", "train")
  if (!dir.exists(tr_lbl_dir)) tr_lbl_dir <- file.path(dir_path, "train", "labels")
  if (!dir.exists(tr_lbl_dir)) tr_lbl_dir <- file.path(dir_path, "train")

  tr_img_dir <- file.path(dir_path, "images", "train")
  if (!dir.exists(tr_img_dir)) tr_img_dir <- file.path(dir_path, "train", "images")
  if (!dir.exists(tr_img_dir)) tr_img_dir <- file.path(dir_path, "train")

  lbl_fls <- list.files(tr_lbl_dir, pattern = "\\.txt$", full.names = TRUE)
  lbl_fls <- lbl_fls[!grepl("(_bal\\d+|_copypaste\\d+)", basename(lbl_fls))]

  # Count objects per class
  class_counts <- setNames(integer(length(all_classes)), all_classes)
  for (lf in lbl_fls) {
    lns <- readLines(lf, warn = FALSE)
    lns <- trimws(lns)[nzchar(trimws(lns))]
    for (ln in lns) {
      vals <- suppressWarnings(as.numeric(strsplit(ln, "\\s+")[[1]]))
      if (length(vals) >= 5L) {
        cid <- as.integer(vals[1])
        if (cid + 1L <= length(all_classes)) {
          cname <- all_classes[cid + 1L]
          class_counts[cname] <- class_counts[cname] + 1L
        }
      }
    }
  }

  if (is.null(target_count)) {
    pos_counts <- class_counts[class_counts > 0L]
    target_count <- if (strategy %in% c("copypaste", "oversample")) {
      if (length(class_counts) > 0L) max(class_counts) else 0L
    } else {
      if (length(pos_counts) > 0L) min(pos_counts) else 0L
    }
  } else {
    target_count <- as.integer(target_count)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Class Balancing (pliman)")
    cli::cli_alert_info("Task: {.val {task}} | Strategy: {.val {strategy}} | Target instances per class: {.val {target_count}}")
    cli::cli_alert_info("Current instance distribution:")
    for (cn in all_classes) {
      status_icon <- if (class_counts[cn] >= target_count) "v" else "x"
      cli::cli_bullets(c(" " = paste0("{.val ", cn, "}: ", class_counts[cn], " instances")))
    }
  }

  if (strategy %in% c("copypaste", "oversample")) {
    minority_classes <- names(class_counts)[class_counts < target_count]
    if (!is.null(classes)) {
      minority_classes <- intersect(minority_classes, classes)
    }

    if (length(minority_classes) == 0L) {
      if (isTRUE(verbose)) cli::cli_alert_success("All target classes already meet or exceed the target count.")
      return(.yolo_wrap_dataset(dir_path, task = task, classes = all_classes))
    }

    # Run Copy-Paste augmentation targeting minority classes
    res <- yolo_dataset_copypaste(
      data = dir_path,
      classes = minority_classes,
      target_count = target_count,
      clean = FALSE,
      seed = seed,
      as_job = as_job,
      verbose = verbose,
      ...
    )
    return(res)
  } else if (strategy == "undersample") {
    # Downsample images containing exclusively majority classes
    majority_classes <- names(class_counts)[class_counts > target_count]
    if (length(majority_classes) == 0L) {
      if (isTRUE(verbose)) cli::cli_alert_info("No majority classes exceed target count for undersampling.")
      return(.yolo_wrap_dataset(dir_path, task = task, classes = all_classes))
    }

    # Find images that contain only majority classes
    img_majority_only <- character(0)
    for (lf in lbl_fls) {
      lns <- readLines(lf, warn = FALSE)
      lns <- trimws(lns)[nzchar(trimws(lns))]
      cids <- vapply(lns, function(x) as.integer(strsplit(x, "\\s+")[[1]][1]), integer(1))
      cnames <- all_classes[cids + 1L]
      if (all(cnames %in% majority_classes)) {
        stem <- tools::file_path_sans_ext(basename(lf))
        img_candidates <- list.files(tr_img_dir, pattern = paste0("^", stem, "\\."), full.names = TRUE)
        if (length(img_candidates) > 0L) {
          img_majority_only <- c(img_majority_only, img_candidates[1])
        }
      }
    }

    if (length(img_majority_only) > 0L) {
      n_to_drop <- min(length(img_majority_only), max(class_counts) - target_count)
      if (n_to_drop > 0L) {
        dropped <- sample(img_majority_only, size = n_to_drop)
        dropped_lbls <- file.path(tr_lbl_dir, paste0(tools::file_path_sans_ext(basename(dropped)), ".txt"))
        unlink(dropped, force = TRUE)
        unlink(dropped_lbls[file.exists(dropped_lbls)], force = TRUE)
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Undersampled {length(dropped)} training image(s) containing exclusively majority classes.")
        }
      }
    }

    unlink(list.files(dir_path, pattern = "\\.cache$", full.names = TRUE, recursive = TRUE), force = TRUE)
    return(.yolo_wrap_dataset(dir_path, task = task, classes = all_classes))
  }
}
# ==============================================================================
# 9. STATISTICAL PROFILING & SUMMARY (yolo_dataset_stats)
# ==============================================================================
#' @title Compute Deep Statistical Profile for a YOLO Dataset
#' @name yolo_dataset_stats
#' @description
#' `yolo_dataset_stats()` computes granular distribution metrics for any YOLO dataset.
#' It reports sample frequencies per split, class proportion percentages, imbalance ratios,
#' and (for object detection/segmentation) bounding box area distributions.
#'
#' @param data A `yolo_dataset` object, `yolo_dataset_cls` object, or directory path.
#' @param verbose Logical. Display formatted tables in the console. Defaults to `TRUE`.
#'
#' @return A `data.frame` containing summary statistics per class.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' stats <- yolo_dataset_stats("dataset_flores_yolo")
yolo_dataset_stats <- function(data, verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  wrapped <- .yolo_wrap_dataset(dir_path, task = task, classes = classes)

  if (task == "classify") {
    df <- wrapped$summary
    if (isTRUE(verbose)) {
      cli::cli_h2("YOLO Classification Dataset Statistics (pliman)")
      print(df, row.names = FALSE)
    }
    return(invisible(df))
  }

  # Detection / Segmentation stats
  lbl_files <- list.files(file.path(dir_path, "labels"), pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  counts <- setNames(integer(length(classes)), classes)
  box_areas <- list()

  for (lf in lbl_files) {
    if (!file.exists(lf)) next
    lines <- tryCatch(readLines(lf, warn = FALSE), error = function(e) character(0))
    lines <- lines[nzchar(trimws(lines))]
    for (ln in lines) {
      vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
      if (length(vals) >= 5L) {
        cid <- as.integer(vals[1]) + 1L
        if (cid >= 1L && cid <= length(classes)) {
          c_name <- classes[cid]
          counts[c_name] <- counts[c_name] + 1L
          if (length(vals) == 5L) {
            box_areas[[c_name]] <- c(box_areas[[c_name]], vals[4] * vals[5])
          }
        }
      }
    }
  }

  df_rows <- list()
  for (cls in classes) {
    areas <- box_areas[[cls]]
    df_rows[[cls]] <- data.frame(
      class = cls,
      objects = counts[cls],
      mean_box_pct = if (length(areas) > 0L) round(mean(areas) * 100, 2) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  res_df <- do.call(rbind, df_rows)
  rownames(res_df) <- NULL

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Object Statistics (pliman)")
    print(res_df, row.names = FALSE)
  }

  invisible(res_df)
}


# ==============================================================================
# 9. TILING & SLIDING WINDOW WITH OVERLAP (yolo_dataset_tiling)
# ==============================================================================
#' @title Tiling and Sliding Window Slicing for High-Resolution YOLO Datasets
#' @name yolo_dataset_tiling
#' @aliases yolo_dataset_tile
#' @description
#' `yolo_dataset_tiling()` slices high-resolution images into uniform overlapping tiles
#' (sliding window) suitable for training YOLO object detection and instance segmentation models.
#'
#' It maps and adapts all existing YOLO bounding boxes (or polygons) to the local coordinate
#' system of each tile, applying the standard truncation/preservation criterion: keeping
#' objects whose center falls within the tile or whose intersection preserves at least 50%
#' of the original area (configurable via `min_area_ratio` and `assign_by`).
#'
#' Crucially, data splitting across `train`, `val`, and `test` is performed at the **source image**
#' level (or inherited from existing split directories), ensuring that all tiles derived from the
#' same original image remain in the same split to strictly prevent data leakage.
#'
#' Empty background tiles (containing no objects) are preserved with empty `.txt` files to provide
#' essential negative background samples for YOLO training.
#'
#' @param data A `yolo_dataset` object, or path to a YOLO dataset directory (or folder containing images and labels).
#' @param path Alias for `data`.
#' @param out_dir Output directory where the tiled dataset will be created. If `NULL` (default), tiles the dataset in-place within the same directory.
#' @param n Optional integer. Number of source images to sample for tiling (e.g. `n = 50`).
#'   If `NULL` (default) and `prop` is `NULL`, all available images are tiled.
#' @param prop Optional numeric between 0 and 1. Proportion of source images to sample for tiling
#'   (e.g. `prop = 0.20` for 20% of the dataset). Ignored if `n` is provided.
#' @param tile_size Integer or numeric vector of length 2 giving `c(width, height)`. Standard multiple of 32 (default: `640L`).
#' @param overlap Numeric fraction of overlap between adjacent tiles (default: `0.20`, i.e., 20%). Ignored if `n_tiles` is set.
#' @param n_tiles Optional integer or numeric vector of length 2 giving exact number of tiles per axis (e.g. `n_tiles = 5` or `n_tiles = c(5, 5)`). If provided, overrides `overlap`.
#' @param min_area_ratio Minimum proportion of original object area preserved when truncated by a tile boundary (default: `0.50`, i.e., 50%).
#' @param assign_by Criterion for keeping truncated annotations: `"center_or_overlap"` (default, keeps if center is inside tile OR area ratio >= `min_area_ratio`), `"center_and_overlap"` (keeps only if center is inside AND area ratio >= `min_area_ratio`), `"center"` (only if center is inside), or `"overlap"` (only if area ratio >= `min_area_ratio`).
#' @param min_box_size Minimum width/height in pixels for a cropped bounding box in the tile (default: `4`).
#' @param include_empty Logical. Whether to generate empty `.txt` files and keep background tiles containing no objects (default: `TRUE`).
#' @param max_bg_ratio Optional numeric fraction. Maximum ratio of background tiles relative to object tiles. If `NULL` (default), keeps all background tiles.
#' @param split_ratio Numeric vector of length 3 giving `c(train, val, test)` proportions (default: `c(0.70, 0.20, 0.10)`). Ignored if the source dataset is already partitioned into splits.
#' @param pad Logical. If an image is smaller than `tile_size`, whether to pad it with `bg_col` to `tile_size` (default: `TRUE`).
#' @param bg_col Background color for padding when `pad = TRUE` (default: `"black"`).
#' @param class_names Character vector of class names. If `NULL` (default), auto-detected from source `data.yaml` or set to `c("object")`.
#' @param seed Optional integer seed for reproducible splitting (default: `42`).
#' @param verbose Logical. Whether to print progress and statistics. Defaults to `TRUE`.
#' @param progress Logical. Whether to show a CLI progress bar. Defaults to `TRUE`.
#' @param ... Additional arguments.
#'
#' @return A `yolo_dataset_tiled` object (inheriting from `yolo_dataset`) pointing to the generated dataset directory.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Tile a high-resolution YOLO dataset with 20% overlap
#' tiled <- yolo_dataset_tiling("highres_dataset",
#'                              out_dir = "tiled_dataset",
#'                              tile_size = 640,
#'                              overlap = 0.20)
#' plot(tiled)
#'
#' # 2. Uniform 5x5 grid per image with 50% minimum area preservation
#' tiled_5x5 <- yolo_dataset_tiling("highres_dataset",
#'                                  out_dir = "tiled_5x5",
#'                                  tile_size = 640,
#'                                  n_tiles = c(5, 5),
#'                                  min_area_ratio = 0.50)
#'
#' # 3. Tile a subset of images (e.g. 50 images or 20% from a large 500-image dataset)
#' tiled_sub <- yolo_dataset_tiling("highres_dataset",
#'                                  out_dir = "tiled_sub",
#'                                  prop = 0.20,
#'                                  tile_size = 640)
yolo_dataset_tiling <- function(data = NULL,
                                out_dir = NULL,
                                path = NULL,
                                n = NULL,
                                prop = NULL,
                                tile_size = 640L,
                                overlap = 0.20,
                                n_tiles = NULL,
                                min_area_ratio = 0.50,
                                assign_by = c("center_or_overlap", "center_and_overlap", "center", "overlap"),
                                min_box_size = 4,
                                include_empty = TRUE,
                                max_bg_ratio = NULL,
                                split_ratio = c(0.70, 0.20, 0.10),
                                pad = TRUE,
                                bg_col = "black",
                                class_names = NULL,
                                seed = 42,
                                verbose = TRUE,
                                progress = TRUE,
                                ...) {
  if (is.null(data)) {
    if (!is.null(path)) {
      data <- path
    } else {
      cli::cli_abort("Argument {.arg data} (or {.arg path}) must be provided.")
    }
  }

  extra_args <- list(...)
  if (!is.null(extra_args$train_prop) && is.numeric(extra_args$train_prop)) {
    tp <- extra_args$train_prop[1]
    if (tp >= 1.0) {
      split_ratio <- c(1.0, 0.0, 0.0)
    } else {
      split_ratio <- c(tp, 1.0 - tp, 0.0)
    }
  }

  assign_by <- match.arg(assign_by)
  mode_int <- switch(assign_by,
    "center_or_overlap" = 0L,
    "center" = 1L,
    "overlap" = 2L,
    "center_and_overlap" = 3L
  )

  tile_w <- as.integer(tile_size[1])
  tile_h <- as.integer(if (length(tile_size) >= 2L) tile_size[2] else tile_size[1])

  if (tile_w <= 0L || tile_h <= 0L) {
    cli::cli_abort("{.arg tile_size} must be positive integer dimensions.")
  }

  if (tile_w %% 32L != 0L || tile_h %% 32L != 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_warning("tile_size ({tile_w}x{tile_h}) is not a multiple of 32. YOLO models typically perform best with multiples of 32.")
    }
  }

  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Source dataset directory {.path {dir_path}} does not exist.")
  }

  detected_task <- .yolo_detect_task(dir_path)
  if (is.null(class_names) || length(class_names) == 0L) {
    class_names <- .yolo_get_classes(dir_path)
  }
  if (is.null(class_names) || length(class_names) == 0L) {
    class_names <- c("object")
  }

  # Background color conversion
  bg_rgb <- tryCatch(grDevices::col2rgb(bg_col)[, 1], error = function(e) c(0, 0, 0))
  bg_num <- as.numeric(bg_rgb)

  # Discover images and their existing or assigned split
  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"

  # Case A: Standard YOLO directory structure with images/<split>
  splits_found <- character(0)
  if (dir.exists(file.path(dir_path, "images"))) {
    for (sp in c("train", "val", "test")) {
      if (dir.exists(file.path(dir_path, "images", sp))) {
        splits_found <- c(splits_found, sp)
      }
    }
  }

  source_records <- list()

  if (length(splits_found) > 0L) {
    for (sp in splits_found) {
      sp_im_dir <- file.path(dir_path, "images", sp)
      sp_lb_dir <- file.path(dir_path, "labels", sp)
      fls <- list.files(sp_im_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
      fls <- fls[!grepl("_tile_\\d+_\\d+", basename(fls))]
      for (f in fls) {
        stem <- tools::file_path_sans_ext(basename(f))
        lbl_cand <- file.path(sp_lb_dir, paste0(stem, ".txt"))
        source_records[[length(source_records) + 1L]] <- list(
          img_path = f,
          lbl_path = if (file.exists(lbl_cand)) lbl_cand else NULL,
          split = sp,
          stem = stem
        )
      }
    }
  } else if (dir.exists(file.path(dir_path, "images"))) {
    # Case B: images/ and labels/ without train/val/test
    sp_im_dir <- file.path(dir_path, "images")
    sp_lb_dir <- file.path(dir_path, "labels")
    fls <- list.files(sp_im_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
    fls <- fls[!grepl("_tile_\\d+_\\d+", basename(fls))]
    for (f in fls) {
      stem <- tools::file_path_sans_ext(basename(f))
      lbl_cand <- file.path(sp_lb_dir, paste0(stem, ".txt"))
      source_records[[length(source_records) + 1L]] <- list(
        img_path = f,
        lbl_path = if (file.exists(lbl_cand)) lbl_cand else NULL,
        split = NULL,
        stem = stem
      )
    }
  } else {
    # Case C: Flat folder containing images and .txt files side-by-side
    fls <- list.files(dir_path, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
    fls <- fls[!grepl("_tile_\\d+_\\d+", basename(fls))]
    for (f in fls) {
      stem <- tools::file_path_sans_ext(basename(f))
      lbl_cand <- file.path(dir_path, paste0(stem, ".txt"))
      source_records[[length(source_records) + 1L]] <- list(
        img_path = f,
        lbl_path = if (file.exists(lbl_cand)) lbl_cand else NULL,
        split = NULL,
        stem = stem
      )
    }
  }

  n_total_orig <- length(source_records)
  if (n_total_orig == 0L) {
    cli::cli_abort("No valid images found in {.path {dir_path}} matching {.val {img_pattern}}.")
  }

  # Subsampling source images via n or prop
  n_sample <- n_total_orig
  if (!is.null(n)) {
    n_val <- as.integer(n[1])
    if (n_val <= 0L) cli::cli_abort("{.arg n} must be a positive integer.")
    n_sample <- min(n_total_orig, n_val)
  } else if (!is.null(prop)) {
    p_val <- as.numeric(prop[1])
    if (p_val <= 0 || p_val > 1) cli::cli_abort("{.arg prop} must be a number between 0 and 1.")
    n_sample <- max(1L, min(n_total_orig, as.integer(round(n_total_orig * p_val))))
  }

  if (n_sample < n_total_orig) {
    if (!is.null(seed)) set.seed(seed)
    if (length(splits_found) > 0L) {
      # Stratified sampling preserving proportions of existing splits
      sp_names <- vapply(source_records, function(x) x$split, character(1))
      sp_list <- split(seq_along(source_records), sp_names)
      sp_props <- vapply(sp_list, length, integer(1)) / n_total_orig
      sp_targets <- stats::setNames(as.integer(round(n_sample * sp_props)), names(sp_list))
      if (n_sample >= length(sp_list)) {
        sp_targets[sp_targets == 0L] <- 1L
      }
      diff_c <- n_sample - sum(sp_targets)
      if (diff_c > 0L) {
        largest <- which.max(vapply(sp_list, length, integer(1)))
        sp_targets[largest] <- sp_targets[largest] + diff_c
      } else if (diff_c < 0L) {
        while (diff_c < 0L) {
          largest <- which.max(sp_targets)
          if (sp_targets[largest] <= 1L) break
          sp_targets[largest] <- sp_targets[largest] - 1L
          diff_c <- diff_c + 1L
        }
      }
      sel_indices <- integer(0)
      for (sp_i in names(sp_list)) {
        cand <- sp_list[[sp_i]]
        take_k <- min(length(cand), sp_targets[[sp_i]])
        if (take_k > 0L) {
          sel <- if (length(cand) == 1L) cand else sample(cand, take_k)
          sel_indices <- c(sel_indices, sel)
        }
      }
      source_records <- source_records[sort(sel_indices)]
    } else {
      # Flat sampling across all source images
      sel_indices <- if (n_total_orig == 1L) 1L else sample(n_total_orig, n_sample)
      source_records <- source_records[sort(sel_indices)]
    }

    if (isTRUE(verbose)) {
      sampling_desc <- if (!is.null(n)) paste0("n = ", n_sample) else paste0("prop = ", round(prop[1] * 100, 1), "% (", n_sample, " images)")
      cli::cli_alert_info("Subsampling {n_sample} of {n_total_orig} images for tiling ({sampling_desc}).")
    }
  }

  n_source_imgs <- length(source_records)

  # Assign splits at the SOURCE IMAGE level if not pre-existing (ZERO data leakage!)
  if (is.null(source_records[[1L]]$split)) {
    if (!is.null(seed)) set.seed(seed)
    shuffled_idx <- if (n_source_imgs > 1L) sample(n_source_imgs) else 1L
    sp_norm <- as.numeric(split_ratio)
    if (length(sp_norm) == 1L) sp_norm <- c(sp_norm, max(0.0, 1.0 - sp_norm), 0.0)
    else if (length(sp_norm) == 2L) sp_norm <- c(sp_norm, max(0.0, 1.0 - sum(sp_norm)))
    sp_sum <- sum(sp_norm)
    sp_norm <- if (sp_sum > 0) sp_norm / sp_sum else c(1.0, 0.0, 0.0)

    assigned_splits <- character(n_source_imgs)
    if (n_source_imgs == 1L) {
      assigned_splits[1] <- if (sp_norm[1] > 0) "train" else if (sp_norm[2] > 0) "val" else "test"
    } else {
      counts <- as.integer(round(n_source_imgs * sp_norm))
      if (counts[1] < 1L && sp_norm[1] > 0) counts[1] <- 1L
      diff <- n_source_imgs - sum(counts)
      if (diff > 0L) {
        counts[1] <- counts[1] + diff
      } else if (diff < 0L) {
        while (diff < 0L) {
          if (counts[3] > 0L) counts[3] <- counts[3] - 1L
          else if (counts[2] > 0L) counts[2] <- counts[2] - 1L
          else if (counts[1] > 1L) counts[1] <- counts[1] - 1L
          else break
          diff <- diff + 1L
        }
      }

      cur <- 1L
      if (counts[1] > 0L && cur <= n_source_imgs) {
        idx <- seq(from = cur, length.out = min(counts[1], n_source_imgs - cur + 1L))
        assigned_splits[shuffled_idx[idx]] <- "train"
        cur <- cur + length(idx)
      }
      if (counts[2] > 0L && cur <= n_source_imgs) {
        idx <- seq(from = cur, length.out = min(counts[2], n_source_imgs - cur + 1L))
        assigned_splits[shuffled_idx[idx]] <- "val"
        cur <- cur + length(idx)
      }
      if (counts[3] > 0L && cur <= n_source_imgs) {
        idx <- seq(from = cur, length.out = min(counts[3], n_source_imgs - cur + 1L))
        assigned_splits[shuffled_idx[idx]] <- "test"
        cur <- cur + length(idx)
      }
      assigned_splits[assigned_splits == ""] <- "train"
    }

    for (i in seq_len(n_source_imgs)) {
      source_records[[i]]$split <- assigned_splits[i]
    }
  }

  # Prepare destination directories
  out_dir <- if (!is.null(out_dir)) normalizePath(out_dir[1], winslash = "/", mustWork = FALSE) else dir_path
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  active_splits <- unique(vapply(source_records, function(x) x$split, character(1)))
  for (sp in active_splits) {
    dir.create(file.path(out_dir, "images", sp), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "labels", sp), recursive = TRUE, showWarnings = FALSE)
  }

  # Grid calculation helper
  .calc_axis_starts <- function(dim_size, t_size, n_t, ov) {
    if (dim_size <= t_size) return(1L)
    if (!is.null(n_t)) {
      nx <- max(1L, as.integer(n_t))
      if (nx == 1L) return(1L)
      st <- (dim_size - t_size) / (nx - 1L)
      return(as.integer(round(1.0 + (0:(nx - 1L)) * st)))
    }
    step <- max(1L, as.integer(round(t_size * (1.0 - ov))))
    nx <- as.integer(ceiling((dim_size - t_size) / step)) + 1L
    st <- (dim_size - t_size) / (nx - 1L)
    as.integer(round(1.0 + (0:(nx - 1L)) * st))
  }

  # Calculate total tiles across dataset for progress bar
  total_tiles_expected <- 0L
  img_meta <- vector("list", n_source_imgs)

  for (i in seq_len(n_source_imgs)) {
    rec <- source_records[[i]]
    w_cur <- tile_w; h_cur <- tile_h
    im_temp <- try(image_import(rec$img_path), silent = TRUE)
    if (!inherits(im_temp, "try-error") && is_image(im_temp)) {
      w_cur <- as.integer(dim(im_temp)[1])
      h_cur <- as.integer(dim(im_temp)[2])
    }

    xs <- .calc_axis_starts(w_cur, tile_w, if (!is.null(n_tiles)) n_tiles[1] else NULL, overlap)
    ys <- .calc_axis_starts(h_cur, tile_h, if (!is.null(n_tiles) && length(n_tiles) >= 2L) n_tiles[2] else if (!is.null(n_tiles)) n_tiles[1] else NULL, overlap)
    img_meta[[i]] <- list(w = w_cur, h = h_cur, xs = xs, ys = ys)
    total_tiles_expected <- total_tiles_expected + (length(xs) * length(ys))
  }

  show_pb <- isTRUE(verbose) && isTRUE(progress)
  if (show_pb) {
    pb <- cli::cli_progress_bar(
      name = "Tiling dataset",
      total = total_tiles_expected,
      format = "{cli::pb_spin} Tiling dataset [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      clear = FALSE
    )
  }

  total_tiles_created <- 0L
  tiles_with_objs <- 0L
  tiles_empty <- 0L

  for (i in seq_len(n_source_imgs)) {
    rec <- source_records[[i]]
    im <- try(image_import(rec$img_path), silent = TRUE)
    if (inherits(im, "try-error") || !is_image(im)) {
      if (show_pb) cli::cli_progress_update(id = pb, inc = length(img_meta[[i]]$xs) * length(img_meta[[i]]$ys))
      next
    }

    W <- dim(im)[1]
    H <- dim(im)[2]
    cur_split <- rec$split
    stem <- rec$stem
    ext <- tolower(tools::file_ext(rec$img_path))
    if (!nzchar(ext)) ext <- "jpg"

    # Parse annotations
    raw_lines <- if (!is.null(rec$lbl_path) && file.exists(rec$lbl_path)) {
      readLines(rec$lbl_path, warn = FALSE)
    } else {
      character(0)
    }
    raw_lines <- raw_lines[nzchar(trimws(raw_lines))]

    boxes_list <- list()
    polys_list <- list()
    poly_cids <- numeric(0)

    for (ln in raw_lines) {
      vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
      if (length(vals) < 5L || any(is.na(vals))) next
      cid <- vals[1]

      if (length(vals) == 5L) {
        # Detection: class_id xc yc w h
        xc <- vals[2]; yc <- vals[3]; bw <- vals[4]; bh <- vals[5]
        b_x1 <- (xc - bw * 0.5) * W
        b_x2 <- (xc + bw * 0.5) * W
        b_y1 <- (yc - bh * 0.5) * H
        b_y2 <- (yc + bh * 0.5) * H
        boxes_list[[length(boxes_list) + 1L]] <- c(cid, b_x1, b_y1, b_x2, b_y2)
      } else if (length(vals) >= 6L && (length(vals) - 1L) %% 2L == 0L) {
        # Segmentation: class_id x1 y1 x2 y2 ... xn yn
        coords <- vals[-1]
        px <- coords[seq(1, length(coords), by = 2L)] * W
        py <- coords[seq(2, length(coords), by = 2L)] * H
        interleaved <- as.vector(rbind(px, py))
        polys_list[[length(polys_list) + 1L]] <- interleaved
        poly_cids <- c(poly_cids, cid)
      }
    }

    boxes_mat <- if (length(boxes_list) > 0L) {
      do.call(rbind, boxes_list)
    } else {
      matrix(numeric(0), nrow = 0L, ncol = 5L)
    }

    xs <- img_meta[[i]]$xs
    ys <- img_meta[[i]]$ys

    for (r_idx in seq_along(ys)) {
      y1 <- ys[r_idx]
      for (c_idx in seq_along(xs)) {
        x1 <- xs[c_idx]

        # Crop tile using C++ OpenMP kernel
        tile_img <- crop_tile_cpp(im, x1, y1, tile_w, tile_h, bg_num)

        # Coordinate frame in continuous space [0, W]
        t_x1 <- static_cast_double <- x1 - 1.0
        t_x2 <- t_x1 + tile_w
        t_y1 <- y1 - 1.0
        t_y2 <- t_y1 + tile_h

        # Clip boxes
        res_b <- clip_yolo_boxes_cpp(
          boxes = boxes_mat,
          tile_x1 = t_x1, tile_y1 = t_y1,
          tile_x2 = t_x2, tile_y2 = t_y2,
          tile_w = tile_w, tile_h = tile_h,
          min_area_ratio = min_area_ratio,
          mode = mode_int,
          min_box_size = min_box_size
        )

        # Clip polygons (if any)
        res_p <- if (length(polys_list) > 0L) {
          clip_yolo_polygons_cpp(
            polygons = polys_list,
            class_ids = poly_cids,
            tile_x1 = t_x1, tile_y1 = t_y1,
            tile_x2 = t_x2, tile_y2 = t_y2,
            tile_w = tile_w, tile_h = tile_h,
            min_area_ratio = min_area_ratio,
            mode = mode_int,
            max_vertices = 32L
          )
        } else {
          list(class_ids = integer(0), polygons = list())
        }

        # Build local YOLO label lines
        tile_lines <- character(0)
        if (nrow(res_b) > 0L) {
          for (bi in seq_len(nrow(res_b))) {
            tile_lines <- c(tile_lines, sprintf("%d %.6f %.6f %.6f %.6f",
                                                as.integer(res_b[bi, 1]),
                                                res_b[bi, 2], res_b[bi, 3],
                                                res_b[bi, 4], res_b[bi, 5]))
          }
        }
        if (length(res_p$polygons) > 0L) {
          for (pi in seq_along(res_p$polygons)) {
            tile_lines <- c(tile_lines, paste(c(as.integer(res_p$class_ids[pi]), sprintf("%.6f", res_p$polygons[[pi]])), collapse = " "))
          }
        }

        has_objects <- (length(tile_lines) > 0L)

        # Filter empty background tiles if requested
        if (!has_objects) {
          if (!isTRUE(include_empty)) {
            if (show_pb) cli::cli_progress_update(id = pb)
            next
          }
          if (!is.null(max_bg_ratio) && is.numeric(max_bg_ratio)) {
            if (tiles_empty >= (max_bg_ratio * max(1L, tiles_with_objs))) {
              if (show_pb) cli::cli_progress_update(id = pb)
              next
            }
          }
          tiles_empty <- tiles_empty + 1L
        } else {
          tiles_with_objs <- tiles_with_objs + 1L
        }

        # Save tile image and label file
        tile_stem <- sprintf("%s_tile_%d_%d", stem, r_idx, c_idx)
        tile_fn <- paste0(tile_stem, ".", ext)
        tile_txt_fn <- paste0(tile_stem, ".txt")

        dest_img <- file.path(out_dir, "images", cur_split, tile_fn)
        dest_lbl <- file.path(out_dir, "labels", cur_split, tile_txt_fn)

        .write_image_file(tile_img, dest_img, quality = 95)
        writeLines(tile_lines, dest_lbl)

        total_tiles_created <- total_tiles_created + 1L
        if (show_pb) cli::cli_progress_update(id = pb)
      }
    }
  }

  if (show_pb) cli::cli_progress_done(id = pb)

  # Programmatically generate data.yaml
  yaml_lines <- c(paste0("path: ", out_dir))
  if (dir.exists(file.path(out_dir, "images", "train")) &&
      length(list.files(file.path(out_dir, "images", "train"))) > 0L) {
    yaml_lines <- c(yaml_lines, "train: images/train")
  }
  if (dir.exists(file.path(out_dir, "images", "val")) &&
      length(list.files(file.path(out_dir, "images", "val"))) > 0L) {
    yaml_lines <- c(yaml_lines, "val: images/val")
  } else if (dir.exists(file.path(out_dir, "images", "train"))) {
    yaml_lines <- c(yaml_lines, "val: images/train")
  }
  if (dir.exists(file.path(out_dir, "images", "test")) &&
      length(list.files(file.path(out_dir, "images", "test"))) > 0L) {
    yaml_lines <- c(yaml_lines, "test: images/test")
  }
  yaml_lines <- c(
    yaml_lines,
    "",
    paste0("nc: ", length(class_names)),
    "names:"
  )
  for (ci in seq_along(class_names)) {
    yaml_lines <- c(yaml_lines, sprintf("  %d: %s", ci - 1L, class_names[ci]))
  }
  writeLines(yaml_lines, file.path(out_dir, "data.yaml"))

  # Invalidate stale YOLO label caches
  unlink(file.path(out_dir, "labels", "train.cache"), force = TRUE)
  unlink(file.path(out_dir, "labels", "val.cache"), force = TRUE)
  unlink(file.path(out_dir, "labels", "test.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Dataset tiling complete! Created {total_tiles_created} tiles.")
    cli::cli_bullets(c(
      "*" = "Source images processed: {n_source_imgs}{if (n_source_imgs < n_total_orig) paste0(' (sampled from ', n_total_orig, ')') else ''}",
      "*" = "Tile dimensions: {tile_w} x {tile_h} pixels",
      "*" = "Tiles with objects: {tiles_with_objs} ({round(tiles_with_objs / max(1L, total_tiles_created) * 100, 1)}%)",
      "*" = "Empty background tiles: {tiles_empty} ({round(tiles_empty / max(1L, total_tiles_created) * 100, 1)}%)",
      "*" = "Classes ({length(class_names)}): {.val {paste(class_names, collapse = ', ')}}",
      "v" = "Data leakage strictly prevented: all tiles from the same source image share the same split.",
      "v" = "Config file written: {.file data.yaml}."
    ))
  }

  out_wrap <- .yolo_wrap_dataset(out_dir, task = detected_task, classes = class_names)
  class(out_wrap) <- unique(c("yolo_dataset_tiled", class(out_wrap)))
  out_wrap$is_tiled <- TRUE
  out_wrap$tile_size <- c(tile_w, tile_h)
  out_wrap$overlap <- overlap
  out_wrap$n_source_images_tiled <- n_source_imgs
  out_wrap$n_source_images_total <- n_total_orig
  out_wrap$n_tiles_total <- total_tiles_created
  out_wrap$n_tiles_with_objects <- tiles_with_objs
  out_wrap$n_tiles_background <- tiles_empty

  invisible(out_wrap)
}

#' @export
yolo_dataset_tile <- yolo_dataset_tiling

#' @method plot yolo_dataset_tiled
#' @export
plot.yolo_dataset_tiled <- function(x, ...) {
  plot.yolo_dataset(x, ...)
}

#' @method print yolo_dataset_tiled
#' @export
print.yolo_dataset_tiled <- function(x, ...) {
  cli::cli_h2("Tiled YOLO Dataset (pliman)")
  cli::cli_bullets(c(
    "*" = "Directory: {.path {x$dir}}",
    "*" = "Tile size: {.val {paste(x$tile_size, collapse = 'x')}}",
    "*" = "Total tiles created: {.val {x$n_tiles_total}} ({.val {x$n_tiles_with_objects}} with objects, {.val {x$n_tiles_background}} background)",
    "*" = "Classes: {.val {paste(x$classes, collapse = ', ')}}"
  ))
  invisible(x)
}

#' @title Remove Generated Tiles from a YOLO Dataset
#' @name yolo_dataset_tiling_clean
#' @aliases yolo_dataset_tile_clean yolo_dataset_clean_tile yolo_dataset_clean_tiling
#' @description
#' `yolo_dataset_tiling_clean()` removes all synthetic tile images and labels
#' (matching the pattern `_tile_[0-9]+_[0-9]+`) previously generated by [yolo_dataset_tiling()],
#' reverting the dataset back to its pre-tiled state.
#'
#' It automatically invalidates stale YOLO cache files (`train.cache`, `val.cache`, `test.cache`)
#' to ensure subsequent YOLO training does not reference deleted files.
#'
#' @param data A `yolo_dataset` object, or directory path.
#' @param verbose Logical. Whether to print informative messages (default: `TRUE`).
#'
#' @return An invisible `yolo_dataset` object pointing to the cleaned dataset.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Tile dataset in-place
#' ds <- yolo_dataset_tiling("dataset_complex", tile_size = 640)
#'
#' # 2. Clean up tiles and return to pre-tiled state
#' ds_clean <- yolo_dataset_tiling_clean("dataset_complex")
yolo_dataset_tiling_clean <- function(data, verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)

  img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"

  # Find all tile images and labels matching _tile_\d+_\d+
  all_imgs <- list.files(dir_path, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  tile_imgs <- all_imgs[grepl("_tile_\\d+_\\d+", basename(all_imgs))]

  all_lbls <- list.files(dir_path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  tile_lbls <- all_lbls[grepl("_tile_\\d+_\\d+", basename(all_lbls))]

  n_removed_img <- length(tile_imgs)
  n_removed_lbl <- length(tile_lbls)

  if (n_removed_img > 0L) {
    unlink(tile_imgs, force = TRUE)
  }
  if (n_removed_lbl > 0L) {
    unlink(tile_lbls, force = TRUE)
  }

  n_remaining_imgs <- length(all_imgs) - n_removed_img

  # Invalidate stale cache files
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)

  if (isTRUE(verbose)) {
    if (n_removed_img == 0L && n_removed_lbl == 0L) {
      cli::cli_alert_info("No tile files matching `_tile_[0-9]+_[0-9]+` found in {.path {dir_path}}.")
    } else {
      cli::cli_alert_success("Cleaned tile files from {.path {dir_path}}.")
      cli::cli_bullets(c(
        "*" = "Tile images removed: {n_removed_img}",
        "*" = "Tile label files removed: {n_removed_lbl}",
        "*" = "Remaining images: {n_remaining_imgs}",
        "v" = "YOLO label caches invalidated."
      ))
    }
  }

  invisible(.yolo_wrap_dataset(dir_path, task = task, classes = classes))
}

#' @rdname yolo_dataset_tiling_clean
#' @export
yolo_dataset_tile_clean <- yolo_dataset_tiling_clean

#' @rdname yolo_dataset_tiling_clean
#' @export
yolo_dataset_clean_tile <- yolo_dataset_tiling_clean

#' @rdname yolo_dataset_tiling_clean
#' @export
yolo_dataset_clean_tiling <- yolo_dataset_tiling_clean

#' @title Undo YOLO Dataset Operations (Ctrl + Z)
#' @name yolo_dataset_undo
#' @description
#' Reverses the last labeling, export, or autolabeling operation performed on a YOLO dataset,
#' functioning like a "Ctrl + Z" action.
#'
#' If source images were moved (`move = TRUE`), they are safely moved back to their original
#' source directory. If images were copied (`move = FALSE`), the copies in `images/` are removed.
#' All generated labels (`labels/`), review preview images (`review/`), and YAML configs are
#' cleanly purged. If no history manifest is found, it safely restores all images from `images/`
#' back to the dataset folder (or a custom `to` directory).
#'
#' @param data Path to the dataset directory (e.g. `"dataset_rotulado"`, `"."`) or a `yolo_dataset` object.
#' @param to Optional destination folder to return the images to (used if original location cannot be determined
#'   or to override original locations). Defaults to restoring images directly inside the dataset directory.
#' @param steps Integer. Number of operations to undo (default `1L`).
#' @param purge_empty Logical. Whether to remove empty dataset directories (`images/`, `labels/`, `review/`) after restoring images. Defaults to `TRUE`.
#' @param force Logical. If `TRUE`, forces undo even if `data` is a root/system directory or a parent container holding multiple child datasets. Defaults to `FALSE`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#'
#' @return Invisible list containing details of restored images, removed labels, and cleaned files.
#' @export
yolo_dataset_undo <- function(data = ".",
                              to = NULL,
                              steps = 1L,
                              purge_empty = TRUE,
                              force = FALSE,
                              verbose = TRUE) {
  dir_path <- tryCatch(.yolo_get_dir(data), error = function(e) {
    if (is.character(data) && length(data) == 1L) normalizePath(data, winslash = "/", mustWork = FALSE) else "."
  })
  hist_file <- file.path(dir_path, ".yolo_history.rds")

  # Safety check 1: Prevent running undo on filesystem/drive root or user system folders
  norm_dir <- tolower(normalizePath(dir_path, winslash = "/", mustWork = FALSE))
  is_root_dir <- grepl("^[a-z]:/?$", norm_dir) || norm_dir == "/" || norm_dir == ""

  sys_dirs <- tolower(normalizePath(c(
    path.expand("~"),
    path.expand("~/Downloads"),
    path.expand("~/Desktop"),
    path.expand("~/Documents"),
    path.expand("~/Pictures"),
    Sys.getenv("USERPROFILE"),
    Sys.getenv("HOME"),
    Sys.getenv("SystemRoot")
  ), winslash = "/", mustWork = FALSE))
  sys_dirs <- unique(sys_dirs[nzchar(sys_dirs)])
  is_sys_dir <- norm_dir %in% sys_dirs

  # Safety check 2: Check if dir_path is a container folder holding child datasets
  sub_dirs <- list.dirs(dir_path, recursive = FALSE, full.names = TRUE)
  child_datasets <- if (length(sub_dirs) > 0L) {
    sub_dirs[vapply(sub_dirs, function(sd) {
      file.exists(file.path(sd, "data.yaml")) ||
      file.exists(file.path(sd, "dataset.yaml")) ||
      file.exists(file.path(sd, ".yolo_history.rds")) ||
      (dir.exists(file.path(sd, "images")) && dir.exists(file.path(sd, "labels")))
    }, logical(1))]
  } else {
    character(0)
  }

  is_container_dir <- length(child_datasets) > 0L &&
                      !file.exists(file.path(dir_path, "data.yaml")) &&
                      !file.exists(hist_file)

  if (!isTRUE(force)) {
    if (is_root_dir || is_sys_dir) {
      cli::cli_abort(c(
        "Cannot perform {.fn yolo_dataset_undo} on root/system directory {.path {dir_path}}.",
        "x" = "Performing undo on a root or system folder risks moving or deleting files across your system.",
        "i" = "Please specify the exact dataset folder (e.g., {.code yolo_dataset_undo('path/to/dataset')}).",
        "i" = "If you really intend to force an undo on this directory, pass {.code force = TRUE}."
      ))
    }
    if (is_container_dir) {
      child_names <- basename(child_datasets)
      cli::cli_abort(c(
        "Cannot perform {.fn yolo_dataset_undo} on directory {.path {dir_path}}.",
        "x" = "This folder appears to be a container holding {length(child_names)} child dataset(s): {.val {child_names}}.",
        "i" = "Running undo here would pull images from all child datasets into the root folder.",
        "i" = "Please specify the specific dataset you want to undo, for example: {.code yolo_dataset_undo('{child_datasets[1]}')}.",
        "i" = "If you really intend to force an undo on this container folder, pass {.code force = TRUE}."
      ))
    }
  }

  n_restored <- 0L
  n_cleaned_copies <- 0L
  n_labels_removed <- 0L
  n_reviews_removed <- 0L
  restored_to_dirs <- character()

  if (file.exists(hist_file)) {
    hist_list <- tryCatch(readRDS(hist_file), error = function(e) list())
    if (is.list(hist_list) && length(hist_list) > 0L) {
      n_ops <- min(length(hist_list), as.integer(max(1L, steps[1])))
      ops_to_undo <- rev(tail(hist_list, n_ops))
      remaining_hist <- head(hist_list, length(hist_list) - n_ops)

      for (op in ops_to_undo) {
        img_df <- op$images
        was_moved <- isTRUE(op$moved)

        if (!is.null(img_df) && nrow(img_df) > 0L) {
          for (r in seq_len(nrow(img_df))) {
            d_img <- img_df$dest[r]
            o_img <- img_df$orig[r]
            d_lbl <- img_df$lbl[r]
            d_rev <- img_df$rev[r]

            # 1. Image restoration or deletion
            if (file.exists(d_img)) {
              if (was_moved) {
                target_dest <- if (!is.null(to) && nzchar(to)) {
                  file.path(normalizePath(to, winslash = "/", mustWork = FALSE), basename(d_img))
                } else if (!is.na(o_img) && nzchar(o_img)) {
                  o_img
                } else if (!is.null(op$input_dir) && nzchar(op$input_dir)) {
                  file.path(op$input_dir, basename(d_img))
                } else {
                  file.path(dir_path, basename(d_img))
                }

                target_dir <- dirname(target_dest)
                if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
                restored_to_dirs <- c(restored_to_dirs, target_dir)

                if (normalizePath(d_img, winslash = "/", mustWork = FALSE) != normalizePath(target_dest, winslash = "/", mustWork = FALSE)) {
                  mv_ok <- file.rename(d_img, target_dest)
                  if (!mv_ok) {
                    file.copy(d_img, target_dest, overwrite = TRUE)
                    unlink(d_img, force = TRUE)
                  }
                }
                n_restored <- n_restored + 1L
              } else {
                unlink(d_img, force = TRUE)
                n_cleaned_copies <- n_cleaned_copies + 1L
              }
            }

            # 2. Label removal
            if (!is.na(d_lbl) && file.exists(d_lbl)) {
              unlink(d_lbl, force = TRUE)
              n_labels_removed <- n_labels_removed + 1L
            }

            # 3. Review removal
            if (!is.na(d_rev) && file.exists(d_rev)) {
              unlink(d_rev, force = TRUE)
              n_reviews_removed <- n_reviews_removed + 1L
            }
          }
        }

        # Purge yaml file if created
        if (!is.null(op$yaml_file) && file.exists(op$yaml_file)) {
          unlink(op$yaml_file, force = TRUE)
        }
      }

      if (length(remaining_hist) > 0L) {
        saveRDS(remaining_hist, hist_file)
      } else {
        unlink(hist_file, force = TRUE)
      }
    }
  } else {
    # Fallback when no .yolo_history.rds exists:
    # First verify that dir_path actually looks like a YOLO dataset
    has_yolo_struct <- file.exists(file.path(dir_path, "data.yaml")) ||
                       file.exists(file.path(dir_path, "dataset.yaml")) ||
                       dir.exists(file.path(dir_path, "images")) ||
                       dir.exists(file.path(dir_path, "labels"))

    if (!has_yolo_struct && !isTRUE(force)) {
      cli::cli_abort(c(
        "Directory {.path {dir_path}} does not appear to be a valid YOLO dataset.",
        "x" = "No {.file .yolo_history.rds}, {.file data.yaml}, or standard dataset structure ({.path images/}, {.path labels/}) was found.",
        "i" = "Please verify the path or specify the exact dataset folder."
      ))
    }

    img_patt <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
    all_imgs <- character(0)
    if (dir.exists(file.path(dir_path, "images"))) {
      all_imgs <- list.files(file.path(dir_path, "images"), pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    }
    if (length(all_imgs) == 0L) {
      split_img_dirs <- file.path(dir_path, c("train", "val", "valid", "test"), "images")
      split_img_dirs <- split_img_dirs[dir.exists(split_img_dirs)]
      if (length(split_img_dirs) > 0L) {
        all_imgs <- list.files(split_img_dirs, pattern = img_patt, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      }
    }
    all_imgs <- all_imgs[!grepl("(_aug\\d+|_tile_\\d+_\\d+|^review_)", basename(all_imgs))]

    target_dir <- if (!is.null(to) && nzchar(to)) {
      normalizePath(to, winslash = "/", mustWork = FALSE)
    } else {
      dir_path
    }
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    restored_to_dirs <- c(restored_to_dirs, target_dir)

    parent_data_dir <- file.path(dirname(dir_path), "data")
    from_sibling_data <- FALSE
    if (length(all_imgs) == 0L) {
      # Rescue images if a prior undo moved images to a sibling 'data' directory
      if (dir.exists(parent_data_dir) && normalizePath(parent_data_dir, winslash = "/") != normalizePath(dir_path, winslash = "/")) {
        cand_in_data <- list.files(parent_data_dir, pattern = img_patt, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
        if (length(cand_in_data) > 0L) {
          all_imgs <- cand_in_data
          from_sibling_data <- TRUE
        }
      }
    }

    for (im_f in all_imgs) {
      target_file <- file.path(target_dir, basename(im_f))
      if (normalizePath(im_f, winslash = "/", mustWork = FALSE) != normalizePath(target_file, winslash = "/", mustWork = FALSE)) {
        mv_ok <- file.rename(im_f, target_file)
        if (!mv_ok) {
          file.copy(im_f, target_file, overwrite = TRUE)
          unlink(im_f, force = TRUE)
        }
      }
      n_restored <- n_restored + 1L
    }

    if (isTRUE(from_sibling_data) && dir.exists(parent_data_dir)) {
      rem_data <- list.files(parent_data_dir, all.files = TRUE, no.. = TRUE)
      if (length(rem_data) == 0L) unlink(parent_data_dir, recursive = TRUE, force = TRUE)
    }

    # Remove all label files (.txt)
    all_lbls <- list.files(file.path(dir_path, "labels"), pattern = "\\.txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (length(all_lbls) > 0L) {
      unlink(all_lbls, force = TRUE)
      n_labels_removed <- length(all_lbls)
    }

    # Remove review images
    rev_imgs <- list.files(file.path(dir_path, "review"), full.names = TRUE, recursive = TRUE)
    if (length(rev_imgs) > 0L) {
      unlink(rev_imgs, force = TRUE)
      n_reviews_removed <- length(rev_imgs)
    }

    # Remove yaml
    unlink(file.path(dir_path, "data.yaml"), force = TRUE)
    unlink(file.path(dir_path, "dataset.yaml"), force = TRUE)
  }

  # Invalidate label cache files and purge leftover runs directories
  unlink(file.path(dir_path, "labels", "train.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "val.cache"), force = TRUE)
  unlink(file.path(dir_path, "labels", "test.cache"), force = TRUE)
  if (dir.exists(file.path(dir_path, "runs"))) unlink(file.path(dir_path, "runs"), recursive = TRUE, force = TRUE)
  if (dir.exists(file.path(getwd(), "runs"))) unlink(file.path(getwd(), "runs"), recursive = TRUE, force = TRUE)

  # Purge empty directories
  if (isTRUE(purge_empty)) {
    sub_dirs <- c(
      file.path(dir_path, "images", c("train", "val", "test")),
      file.path(dir_path, "images"),
      file.path(dir_path, "labels", c("train", "val", "test")),
      file.path(dir_path, "labels"),
      file.path(dir_path, "review"),
      file.path(dir_path, c("train", "val", "test"))
    )
    for (sd in sub_dirs) {
      if (dir.exists(sd)) {
        fls <- list.files(sd, all.files = TRUE, no.. = TRUE)
        if (length(fls) == 0L) unlink(sd, recursive = TRUE, force = TRUE)
      }
    }
    if (dir.exists(dir_path) && normalizePath(dir_path, winslash = "/") != normalizePath(".", winslash = "/")) {
      remaining_files <- list.files(dir_path, all.files = TRUE, no.. = TRUE)
      if (length(remaining_files) == 0L) {
        unlink(dir_path, recursive = TRUE, force = TRUE)
      }
    }
  }

  res_dest <- unique(restored_to_dirs)

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Undo (Ctrl + Z)")
    bullets <- character()
    if (n_restored > 0L) {
      bullets <- c(bullets, "v" = sprintf("Restored %d image(s) to: %s", n_restored, paste(res_dest, collapse = ", ")))
    }
    if (n_cleaned_copies > 0L) {
      bullets <- c(bullets, "v" = sprintf("Removed %d copied image(s) from dataset", n_cleaned_copies))
    }
    if (n_labels_removed > 0L) {
      bullets <- c(bullets, "v" = sprintf("Removed %d label file(s) (.txt)", n_labels_removed))
    }
    if (n_reviews_removed > 0L) {
      bullets <- c(bullets, "v" = sprintf("Cleaned %d review preview image(s)", n_reviews_removed))
    }
    if (length(bullets) == 0L) {
      cli::cli_alert_info("No files were modified or pending undo in {.path {dir_path}}.")
    } else {
      cli::cli_bullets(bullets)
      cli::cli_alert_success("Operation successfully reverted!")
    }
  }

  invisible(list(
    restored_images = n_restored,
    cleaned_copies = n_cleaned_copies,
    removed_labels = n_labels_removed,
    removed_reviews = n_reviews_removed,
    target_dirs = res_dest
  ))
}

#' @rdname yolo_dataset_undo
#' @export
yolo_undo <- yolo_dataset_undo

#' @rdname yolo_dataset_undo
#' @export
yolo_autolabel_undo <- yolo_dataset_undo

#' @rdname yolo_dataset_undo
#' @export
yolo_dataset_export_undo <- yolo_dataset_undo

#' @rdname yolo_dataset_undo
#' @export
yolo_dataset_label_undo <- yolo_dataset_undo


# ==============================================================================
# Delete Backup and Temporary Files from a YOLO Dataset
# ==============================================================================

#' Delete Backup and Temporary Files from a YOLO Dataset
#'
#' @description
#' `yolo_dataset_backup_delete()` scans a YOLO dataset directory (recursively across
#' `labels/`, `images/`, and dataset root) and safely deletes backup files generated
#' during manual editing (e.g., `<stem>.txt.bak` from [yolo_dataset_review()]), previous
#' version copies (`*.old`, `*.backup`, `*~`), and history manifest files (`.yolo_history.rds`).
#'
#' This helps clean up disk space, prevents repository clutter, and ensures a pristine
#' directory structure before training or sharing the dataset.
#'
#' @param data Path to the YOLO dataset root directory (e.g. `"dataset_capsulas"`, `"."`)
#'   or an object of class `yolo_dataset`.
#' @param pattern Character string containing a regular expression of file patterns to delete.
#'   Defaults to `"(\\.(bak|old|backup)$|~|\\.yolo_history\\.rds$)"`.
#' @param history Logical. Whether to also delete undo/history manifest files (`.yolo_history.rds`).
#'   Defaults to `TRUE`.
#' @param backup_dirs Logical. Whether to remove entire backup folders if found
#'   (e.g., `"backup"`, `"backups"`, `"labels_backup"`, `"images_backup"`). Defaults to `TRUE`.
#' @param recursive Logical. Whether to search subdirectories recursively (`labels/train`,
#'   `labels/val`, etc.). Defaults to `TRUE`.
#' @param confirm Logical. Whether to ask for confirmation before deleting files in
#'   interactive sessions. Defaults to `interactive()`.
#' @param verbose Logical. Whether to show informative messages with count and reclaimed disk space.
#'   Defaults to `TRUE`.
#'
#' @return An invisible object of class `c("yolo_dataset_backup_delete", "list")` containing:
#'   * `deleted_files`: Character vector of deleted file paths.
#'   * `deleted_dirs`: Character vector of deleted backup directory paths.
#'   * `count`: Total number of deleted files.
#'   * `reclaimed_bytes`: Total number of bytes freed.
#'   * `reclaimed_mb`: Megabytes freed.
#'
#' @aliases yolo_backup_delete yolo_dataset_delete_backup yolo_delete_backup
#' @seealso [yolo_dataset_review()], [yolo_dataset_undo()]
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Delete all .bak and history files from dataset
#' yolo_dataset_backup_delete("dataset_capsulas")
#' }
yolo_dataset_backup_delete <- function(data = ".",
                                       pattern = "(\\.(bak|old|backup)$|~|\\.yolo_history\\.rds$)",
                                       history = TRUE,
                                       backup_dirs = TRUE,
                                       recursive = TRUE,
                                       confirm = interactive(),
                                       verbose = TRUE) {
  dir_path <- tryCatch(.yolo_get_dir(data), error = function(e) {
    if (is.character(data) && length(data) == 1L) normalizePath(data, winslash = "/", mustWork = FALSE) else "."
  })

  if (!dir.exists(dir_path)) {
    cli::cli_abort("Dataset directory {.path {dir_path}} does not exist.")
  }

  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = TRUE)

  # 1. Search for matching backup files
  all_candidates <- list.files(
    path = dir_path,
    pattern = pattern,
    full.names = TRUE,
    recursive = isTRUE(recursive),
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  # Explicitly check for .yolo_history.rds if history is requested
  if (isTRUE(history)) {
    hist_file <- file.path(dir_path, ".yolo_history.rds")
    if (file.exists(hist_file) && !(hist_file %in% all_candidates)) {
      all_candidates <- c(all_candidates, hist_file)
    }
  } else {
    all_candidates <- all_candidates[!grepl("\\.yolo_history\\.rds$", all_candidates, ignore.case = TRUE)]
  }

  # Normalize paths
  all_candidates <- unique(normalizePath(all_candidates, winslash = "/", mustWork = FALSE))
  all_candidates <- all_candidates[file.exists(all_candidates)]

  # 2. Search for backup directories (e.g. backup, backups, labels_backup, images_backup)
  candidate_dirs <- character(0)
  if (isTRUE(backup_dirs)) {
    all_subdirs <- list.dirs(dir_path, full.names = TRUE, recursive = isTRUE(recursive))
    backup_dir_pattern <- "(^|/)(backup|backups|labels_backup|images_backup)$"
    candidate_dirs <- all_subdirs[grepl(backup_dir_pattern, all_subdirs, ignore.case = TRUE)]
    candidate_dirs <- unique(normalizePath(candidate_dirs, winslash = "/", mustWork = FALSE))
    candidate_dirs <- candidate_dirs[dir.exists(candidate_dirs)]
  }

  total_files_count <- length(all_candidates)
  total_dirs_count <- length(candidate_dirs)

  if (total_files_count == 0L && total_dirs_count == 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_success("No backup files or directories found in {.path {basename(dir_path)}}.")
    }
    out <- list(
      deleted_files = character(0),
      deleted_dirs = character(0),
      count = 0L,
      reclaimed_bytes = 0,
      reclaimed_mb = 0
    )
    class(out) <- c("yolo_dataset_backup_delete", "list")
    return(invisible(out))
  }

  # Calculate reclaimed size
  file_sizes <- file.size(all_candidates)
  total_bytes <- sum(file_sizes, na.rm = TRUE)
  total_mb <- round(total_bytes / (1024 * 1024), 2)
  size_fmt <- if (total_bytes < 1024 * 1024) {
    sprintf("%.1f KB", total_bytes / 1024)
  } else {
    sprintf("%.2f MB", total_mb)
  }

  if (isTRUE(verbose)) {
    cli::cli_rule(left = cli::col_cyan("pliman: YOLO Dataset Backup Cleanup"))
    bullets <- c(
      "*" = sprintf("Dataset: {.path %s}", dir_path),
      "*" = sprintf("Found {.strong %d} backup file(s) ({.val %s})", total_files_count, size_fmt)
    )
    if (total_dirs_count > 0L) {
      bullets <- c(bullets, "*" = sprintf("Found {.strong %d} backup director(y/ies)", total_dirs_count))
    }
    cli::cli_bullets(bullets)
  }

  # Confirmation prompt in interactive sessions if requested
  if (isTRUE(confirm)) {
    msg <- sprintf("Permanently delete %d backup file(s) and free %s?", total_files_count, size_fmt)
    ans <- readline(prompt = paste(msg, "[y/N]: "))
    if (!grepl("^[yY](es)?$", trimws(ans))) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning("Cleanup cancelled by user. No files were deleted.")
      }
      out <- list(
        deleted_files = character(0),
        deleted_dirs = character(0),
        count = 0L,
        reclaimed_bytes = 0,
        reclaimed_mb = 0
      )
      class(out) <- c("yolo_dataset_backup_delete", "list")
      return(invisible(out))
    }
  }

  # Execute deletion of files
  deleted_files <- character(0)
  for (f in all_candidates) {
    ok <- tryCatch(unlink(f, force = TRUE) == 0L || !file.exists(f), error = function(e) FALSE)
    if (ok) {
      deleted_files <- c(deleted_files, f)
    }
  }

  # Execute deletion of backup directories
  deleted_dirs <- character(0)
  for (d in candidate_dirs) {
    ok <- tryCatch(unlink(d, recursive = TRUE, force = TRUE) == 0L || !dir.exists(d), error = function(e) FALSE)
    if (ok) {
      deleted_dirs <- c(deleted_dirs, d)
    }
  }

  if (isTRUE(verbose)) {
    dir_info <- if (length(deleted_dirs) > 0L) sprintf(" and %d director(y/ies)", length(deleted_dirs)) else ""
    cli::cli_alert_success(sprintf(
      "Successfully deleted {.strong %d} backup file(s)%s! Reclaimed {.val %s}.",
      length(deleted_files),
      dir_info,
      size_fmt
    ))
  }

  out <- list(
    deleted_files = deleted_files,
    deleted_dirs = deleted_dirs,
    count = length(deleted_files),
    reclaimed_bytes = total_bytes,
    reclaimed_mb = total_mb
  )
  class(out) <- c("yolo_dataset_backup_delete", "list")
  invisible(out)
}

#' @rdname yolo_dataset_backup_delete
#' @export
yolo_backup_delete <- yolo_dataset_backup_delete

#' @rdname yolo_dataset_backup_delete
#' @export
yolo_dataset_delete_backup <- yolo_dataset_backup_delete

#' @rdname yolo_dataset_backup_delete
#' @export
yolo_delete_backup <- yolo_dataset_backup_delete

#' @export
print.yolo_dataset_backup_delete <- function(x, ...) {
  if (x$count == 0L && length(x$deleted_dirs) == 0L) {
    cli::cli_alert_info("No backup files were deleted.")
  } else {
    bullets <- c(
      "v" = sprintf("Deleted backup files       : {.strong %d}", x$count)
    )
    if (length(x$deleted_dirs) > 0L) {
      bullets <- c(bullets, "v" = sprintf("Deleted backup directories : {.strong %d}", length(x$deleted_dirs)))
    }
    bullets <- c(bullets, "v" = sprintf("Disk space reclaimed       : {.val %s MB}", x$reclaimed_mb))
    cli::cli_bullets(bullets)
  }
  invisible(x)
}

.yolo_create_zip <- function(zipfile, source_dir, files_to_zip = NULL) {
  zipfile <- normalizePath(zipfile, winslash = "/", mustWork = FALSE)
  source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)

  if (file.exists(zipfile)) {
    unlink(zipfile, force = TRUE)
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(source_dir)

  if (is.null(files_to_zip)) {
    files_to_zip <- list.files(".", all.files = FALSE, full.names = FALSE)
  }

  success <- FALSE

  # 1. Try 'zip' package if available (pure C cross-platform, fast)
  if (requireNamespace("zip", quietly = TRUE)) {
    res <- tryCatch({
      zip::zip(zipfile = zipfile, files = files_to_zip, recurse = TRUE, mode = "cherry-pick")
      TRUE
    }, error = function(e) {
      tryCatch({
        zip::zipr(zipfile = zipfile, files = files_to_zip, recurse = TRUE)
        TRUE
      }, error = function(e2) FALSE)
    })
    if (isTRUE(res) && file.exists(zipfile) && file.info(zipfile)$size > 0) {
      success <- TRUE
    }
  }

  # 2. Try utils::zip
  if (!success) {
    res <- tryCatch({
      utils::zip(zipfile, files = files_to_zip, flags = "-r9Xq") == 0
    }, error = function(e) FALSE)
    if (isTRUE(res) && file.exists(zipfile) && file.info(zipfile)$size > 0) {
      success <- TRUE
    }
  }

  # 3. Try Windows PowerShell Compress-Archive
  if (!success && .Platform$OS.type == "windows") {
    items_str <- paste(shQuote(files_to_zip), collapse = ", ")
    ps_cmd <- sprintf("Compress-Archive -Path %s -DestinationPath %s -Force", items_str, shQuote(zipfile))
    res <- tryCatch({
      system2("powershell", args = c("-NoProfile", "-Command", ps_cmd), stdout = FALSE, stderr = FALSE) == 0
    }, error = function(e) FALSE)
    if (isTRUE(res) && file.exists(zipfile) && file.info(zipfile)$size > 0) {
      success <- TRUE
    }
  }

  # 4. Try tar
  if (!success && nzchar(Sys.which("tar"))) {
    res <- tryCatch({
      system2(Sys.which("tar"), args = c("-a", "-c", "-f", shQuote(zipfile), files_to_zip), stdout = FALSE, stderr = FALSE) == 0
    }, error = function(e) FALSE)
    if (isTRUE(res) && file.exists(zipfile) && file.info(zipfile)$size > 0) {
      success <- TRUE
    }
  }

  if (!success || !file.exists(zipfile) || file.info(zipfile)$size == 0) {
    cli::cli_abort("Failed to create ZIP archive at {.path {zipfile}}.")
  }

  invisible(zipfile)
}

#' @title Compress a YOLO Dataset into a ZIP Archive
#' @name yolo_dataset_zip
#' @aliases yolo_zip
#' @description
#' Creates a portable `.zip` archive containing a YOLO dataset. By default, the zip
#' file is named after the dataset folder (e.g. `<dataset_name>.zip`) and saved directly
#' in the same location (parent directory) as the dataset itself.
#'
#' @param data A `yolo_dataset` object or character path to a YOLO dataset directory.
#' @param dest_file Character. Optional path or filename for the resulting `.zip` file.
#'   If `NULL` (default), the file is saved as `<dataset_name>.zip` in the parent directory
#'   of the dataset. If a directory path is provided, the zip is saved inside it using the
#'   dataset's name.
#' @param include_review Logical. Whether to include the `review/` visual inspection folder
#'   if present in the dataset. Defaults to `FALSE` to produce clean and compact training archives
#'   suitable for Google Colab, Kaggle, or server deployments.
#' @param overwrite Logical. Whether to overwrite an existing zip file. Defaults to `TRUE`.
#' @param verbose Logical. If `TRUE` (default), prints status and size information.
#'
#' @return The normalized path to the created `.zip` file (invisibly).
#' @export
#' @examplesIf interactive()
#' yolo_dataset_zip("my_dataset")
yolo_dataset_zip <- function(data,
                             dest_file = NULL,
                             include_review = FALSE,
                             overwrite = TRUE,
                             verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  if (is.null(dir_path) || !nzchar(dir_path) || !dir.exists(dir_path)) {
    cli::cli_abort("Dataset directory {.path {data}} does not exist.")
  }
  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = TRUE)
  ds_name  <- basename(dir_path)
  parent_dir <- dirname(dir_path)

  # Resolve destination zip path:
  # By default: saved in the same local of the dataset with the name of the dataset
  if (is.null(dest_file) || !nzchar(dest_file[1])) {
    dest_path <- file.path(parent_dir, paste0(ds_name, ".zip"))
  } else {
    dest_arg <- dest_file[1]
    if (dir.exists(dest_arg)) {
      dest_path <- file.path(dest_arg, paste0(ds_name, ".zip"))
    } else {
      if (!grepl("\\.zip$", dest_arg, ignore.case = TRUE)) {
        dest_path <- paste0(dest_arg, ".zip")
      } else {
        dest_path <- dest_arg
      }
    }
  }
  dest_path <- normalizePath(dest_path, winslash = "/", mustWork = FALSE)

  if (file.exists(dest_path) && !isTRUE(overwrite)) {
    cli::cli_abort("Destination file {.path {dest_path}} already exists and {.arg overwrite = FALSE}.")
  }

  dest_parent <- dirname(dest_path)
  if (!dir.exists(dest_parent)) {
    dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Compression (pliman)")
    cli::cli_alert_info("Compressing dataset {.val {ds_name}} from {.path {dir_path}}...")
  }

  # Discover files and folders to include
  all_items <- list.files(dir_path, all.files = FALSE, full.names = FALSE)
  if (!isTRUE(include_review)) {
    all_items <- all_items[all_items != "review"]
  }
  if (length(all_items) == 0L) {
    cli::cli_abort("No files found to compress inside {.path {dir_path}}.")
  }

  # Build zip archive
  .yolo_create_zip(zipfile = dest_path, source_dir = dir_path, files_to_zip = all_items)

  if (file.exists(dest_path)) {
    sz_mb <- round(file.info(dest_path)$size / (1024 * 1024), 2)
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Dataset zipped successfully: {.path {dest_path}} ({sz_mb} MB)")
      cli::cli_bullets(c(
        "*" = "Dataset: {.val {ds_name}}",
        "*" = "Included items: {.val {paste(all_items, collapse = ', ')}}",
        "*" = "Destination: {.path {dest_path}}"
      ))
    }
  }

  invisible(dest_path)
}

#' @rdname yolo_dataset_zip
#' @export
yolo_zip <- yolo_dataset_zip

#' @title Decompress a YOLO Dataset or Results ZIP Archive
#' @name yolo_dataset_unzip
#' @aliases yolo_unzip
#' @description
#' Extracts a `.zip` archive containing a YOLO dataset or YOLO training results.
#' By default, extracts into a folder named after the `.zip` archive (without the `.zip`
#' extension) located in the same parent directory as the archive.
#'
#' If the extracted directory contains a valid YOLO dataset (`data.yaml` or `images/` directory),
#' the function automatically returns an active `yolo_dataset` S3 object for direct inspection,
#' training with [yolo_train()], or preview with [yolo_dataset_preview()].
#'
#' @param file Character. Path to the `.zip` archive to extract. If the `.zip` extension
#'   is omitted, the function checks for `<file>.zip` in the current working directory.
#' @param exdir Character. Directory where the files should be extracted. If `NULL` (default),
#'   the contents are extracted into a folder named after the `.zip` archive (without extension)
#'   in the same directory as the `.zip` file.
#' @param overwrite Logical. Whether to overwrite existing files during extraction. Defaults to `TRUE`.
#' @param as_dataset Logical. If `TRUE` (default) and the extracted archive is a YOLO dataset,
#'   returns a `yolo_dataset` S3 object. If `FALSE` or if the archive contains training results,
#'   returns the path to the extracted directory.
#' @param verbose Logical. Whether to print informative status messages. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` object if the extracted contents represent a YOLO dataset, or the
#'   path to the extracted folder (invisibly).
#' @export
#' @examplesIf interactive()
#' # Unzip dataset into 'my_dataset/'
#' ds <- yolo_unzip("my_dataset.zip")
#'
#' # Unzip training results from Google Colab
#' res_dir <- yolo_unzip("capsulas_small_results.zip")
#' yolo_results(res_dir)
yolo_dataset_unzip <- function(file,
                               exdir = NULL,
                               overwrite = TRUE,
                               as_dataset = TRUE,
                               verbose = TRUE) {
  if (missing(file) || is.null(file) || !nzchar(as.character(file)[1])) {
    cli::cli_abort("Please provide the path to the {.file .zip} file.")
  }
  file_arg <- as.character(file)[1]

  # Resolve zip path
  zip_path <- NULL
  cand_paths <- c(
    file_arg,
    if (!grepl("\\.zip$", file_arg, ignore.case = TRUE)) paste0(file_arg, ".zip"),
    file.path(getwd(), file_arg),
    if (!grepl("\\.zip$", file_arg, ignore.case = TRUE)) file.path(getwd(), paste0(file_arg, ".zip"))
  )
  for (cp in cand_paths) {
    if (file.exists(cp) && !dir.exists(cp)) {
      zip_path <- normalizePath(cp, winslash = "/", mustWork = TRUE)
      break
    }
  }

  if (is.null(zip_path)) {
    cli::cli_abort("ZIP file {.file {file_arg}} was not found in the current directory or specified path.")
  }

  zip_stem <- sub("\\.zip$", "", basename(zip_path), ignore.case = TRUE)
  parent_dir <- dirname(zip_path)

  # Resolve destination directory
  if (is.null(exdir) || !nzchar(exdir[1])) {
    target_dir <- file.path(parent_dir, zip_stem)
  } else {
    target_dir <- normalizePath(exdir[1], winslash = "/", mustWork = FALSE)
  }

  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO ZIP Archive Extraction (pliman)")
    cli::cli_alert_info("Extracting {.file {basename(zip_path)}} to {.path {target_dir}}...")
  }

  # Extract using utils::unzip
  unzip_res <- tryCatch({
    utils::unzip(zipfile = zip_path, exdir = target_dir, overwrite = isTRUE(overwrite))
  }, error = function(e) {
    cli::cli_abort("Failed to unzip {.file {zip_path}}: {e$message}")
  })

  # Detect if all extracted files were placed inside a single subfolder
  candidate_dirs <- c(
    target_dir,
    list.dirs(target_dir, recursive = FALSE, full.names = TRUE)
  )

  # Locate actual dataset or results directory
  resolved_dir <- target_dir
  for (cd in candidate_dirs) {
    if (file.exists(file.path(cd, "data.yaml")) ||
        (dir.exists(file.path(cd, "images")) && dir.exists(file.path(cd, "labels"))) ||
        (dir.exists(file.path(cd, "train")) && dir.exists(file.path(cd, "val")))) {
      resolved_dir <- normalizePath(cd, winslash = "/", mustWork = FALSE)
      break
    }
  }

  n_extracted <- length(unzip_res)
  sz_mb <- round(file.info(zip_path)$size / (1024 * 1024), 2)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Extraction completed successfully ({n_extracted} files extracted).")
    cli::cli_bullets(c(
      "*" = "Archive: {.file {basename(zip_path)}} ({sz_mb} MB)",
      "*" = "Extracted location: {.path {resolved_dir}}"
    ))
  }

  # If as_dataset = TRUE and it contains a valid YOLO dataset, wrap it
  if (isTRUE(as_dataset)) {
    is_dataset <- file.exists(file.path(resolved_dir, "data.yaml")) ||
                  (dir.exists(file.path(resolved_dir, "images")) && dir.exists(file.path(resolved_dir, "labels"))) ||
                  (dir.exists(file.path(resolved_dir, "train")) && !file.exists(file.path(resolved_dir, "results.csv")))
    if (is_dataset) {
      wrapped <- tryCatch(.yolo_wrap_dataset(resolved_dir), error = function(e) NULL)
      if (!is.null(wrapped)) {
        return(wrapped)
      }
    }
  }

  invisible(resolved_dir)
}

#' @rdname yolo_dataset_unzip
#' @export
yolo_unzip <- yolo_dataset_unzip


# ==============================================================================
# DATASET RESIZING & LABEL RESCALING (yolo_dataset_resize)
# ==============================================================================

.yolo_resize_annotations <- function(raw_lines, orig_w, orig_h, target_w, target_h, mode = "letterbox") {
  if (length(raw_lines) == 0L) return(character(0))
  if (mode == "stretch") {
    # In stretch mode, coordinates are normalized [0, 1] relative to width and height,
    # so they remain numerically unchanged after stretching.
    return(raw_lines)
  }

  scale_lb <- min(target_w / orig_w, target_h / orig_h)
  new_w_lb <- max(1L, round(orig_w * scale_lb))
  new_h_lb <- max(1L, round(orig_h * scale_lb))
  pad_x_lb <- as.numeric((target_w - new_w_lb) %/% 2)
  pad_y_lb <- as.numeric((target_h - new_h_lb) %/% 2)

  scale_crop <- max(target_w / orig_w, target_h / orig_h)
  crop_w <- round(orig_w * scale_crop)
  crop_h <- round(orig_h * scale_crop)
  crop_x <- as.numeric((crop_w - target_w) %/% 2)
  crop_y <- as.numeric((crop_h - target_h) %/% 2)

  res_lines <- character(0)

  for (ln in raw_lines) {
    vals <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
    if (length(vals) < 5L || any(is.na(vals))) next
    cid <- as.integer(vals[1])

    if (length(vals) == 5L) {
      # Detection Bounding Box: class cx cy w h
      cx <- vals[2]; cy <- vals[3]; bw <- vals[4]; bh <- vals[5]

      if (mode == "letterbox") {
        new_cx <- (cx * orig_w * scale_lb + pad_x_lb) / target_w
        new_cy <- (cy * orig_h * scale_lb + pad_y_lb) / target_h
        new_bw <- (bw * orig_w * scale_lb) / target_w
        new_bh <- (bh * orig_h * scale_lb) / target_h

        new_cx <- max(0.0, min(1.0, new_cx))
        new_cy <- max(0.0, min(1.0, new_cy))
        new_bw <- max(0.0, min(1.0, new_bw))
        new_bh <- max(0.0, min(1.0, new_bh))

        res_lines <- c(res_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, new_cx, new_cy, new_bw, new_bh))

      } else if (mode == "crop") {
        x_min_pix <- (cx - bw / 2.0) * orig_w * scale_crop - crop_x
        x_max_pix <- (cx + bw / 2.0) * orig_w * scale_crop - crop_x
        y_min_pix <- (cy - bh / 2.0) * orig_h * scale_crop - crop_y
        y_max_pix <- (cy + bh / 2.0) * orig_h * scale_crop - crop_y

        c_xmin <- max(0.0, min(as.numeric(target_w), x_min_pix))
        c_xmax <- max(0.0, min(as.numeric(target_w), x_max_pix))
        c_ymin <- max(0.0, min(as.numeric(target_h), y_min_pix))
        c_ymax <- max(0.0, min(as.numeric(target_h), y_max_pix))

        if ((c_xmax - c_xmin) >= 2.0 && (c_ymax - c_ymin) >= 2.0) {
          new_cx <- (c_xmin + c_xmax) / (2.0 * target_w)
          new_cy <- (c_ymin + c_ymax) / (2.0 * target_h)
          new_bw <- (c_xmax - c_xmin) / target_w
          new_bh <- (c_ymax - c_ymin) / target_h

          new_cx <- max(0.0, min(1.0, new_cx))
          new_cy <- max(0.0, min(1.0, new_cy))
          new_bw <- max(0.0, min(1.0, new_bw))
          new_bh <- max(0.0, min(1.0, new_bh))

          res_lines <- c(res_lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, new_cx, new_cy, new_bw, new_bh))
        }
      }

    } else {
      # Polygon / Segmentation / OBB: class x1 y1 x2 y2 ... xn yn
      coords <- vals[-1]
      xs <- coords[seq(1, length(coords), by = 2)]
      ys <- coords[seq(2, length(coords), by = 2)]

      if (mode == "letterbox") {
        new_xs <- (xs * orig_w * scale_lb + pad_x_lb) / target_w
        new_ys <- (ys * orig_h * scale_lb + pad_y_lb) / target_h
        new_xs <- pmax(0.0, pmin(1.0, new_xs))
        new_ys <- pmax(0.0, pmin(1.0, new_ys))

        interleaved <- as.vector(rbind(new_xs, new_ys))
        res_lines <- c(res_lines, paste(c(cid, sprintf("%.6f", interleaved)), collapse = " "))

      } else if (mode == "crop") {
        new_xs <- (xs * orig_w * scale_crop - crop_x) / target_w
        new_ys <- (ys * orig_h * scale_crop - crop_y) / target_h
        new_xs <- pmax(0.0, pmin(1.0, new_xs))
        new_ys <- pmax(0.0, pmin(1.0, new_ys))

        interleaved <- as.vector(rbind(new_xs, new_ys))
        res_lines <- c(res_lines, paste(c(cid, sprintf("%.6f", interleaved)), collapse = " "))
      }
    }
  }

  res_lines
}

.yolo_resize_single_image <- function(im, target_w, target_h, mode = "letterbox", bg_col = "black", pad_mode = "color", filter = "lanczos") {
  orig_w <- dim(im)[1]
  orig_h <- dim(im)[2]
  dims_orig <- dim(im)
  n_ch <- if (length(dims_orig) >= 3L) dims_orig[3] else 1L
  cm <- attr(im, "colormode") %||% "Color"

  if (mode == "stretch") {
    return(image_resize(im, width = target_w, height = target_h, filter = filter, verbose = FALSE, plot = FALSE))
  }

  if (mode == "letterbox") {
    scale <- min(target_w / orig_w, target_h / orig_h)
    new_w <- max(1L, round(orig_w * scale))
    new_h <- max(1L, round(orig_h * scale))
    pad_x <- (target_w - new_w) %/% 2
    pad_y <- (target_h - new_h) %/% 2

    scaled_im <- image_resize(im, width = new_w, height = new_h, filter = filter, verbose = FALSE, plot = FALSE)
    is_raw_im <- is.raw(scaled_im)

    # Replicate border pixels outward to simulate expanded natural background
    if (pad_mode == "edge") {
      map_x <- pmax(1L, pmin(new_w, seq_len(target_w) - pad_x))
      map_y <- pmax(1L, pmin(new_h, seq_len(target_h) - pad_y))
      if (n_ch >= 3L) {
        arr <- scaled_im[map_x, map_y, 1:n_ch, drop = FALSE]
        return(as_image(arr, colormode = cm))
      } else {
        arr <- scaled_im[map_x, map_y, drop = FALSE]
        return(as_image(arr, colormode = "grayscale"))
      }
    }

    # Fill padding bars with the average color along the respective image border
    if (pad_mode == "mean_edge") {
      map_x <- pmax(1L, pmin(new_w, seq_len(target_w) - pad_x))
      map_y <- pmax(1L, pmin(new_h, seq_len(target_h) - pad_y))
      if (n_ch >= 3L) {
        arr <- scaled_im[map_x, map_y, 1:n_ch, drop = FALSE]
        if (pad_y > 0L) {
          for (c in seq_len(n_ch)) {
            top_val <- mean(as.numeric(scaled_im[, 1L, c]))
            bot_val <- mean(as.numeric(scaled_im[, new_h, c]))
            if (is_raw_im) {
              arr[, seq_len(pad_y), c] <- as.raw(pmax(0L, pmin(255L, round(top_val))))
              arr[, (pad_y + new_h + 1L):target_h, c] <- as.raw(pmax(0L, pmin(255L, round(bot_val))))
            } else {
              arr[, seq_len(pad_y), c] <- top_val
              arr[, (pad_y + new_h + 1L):target_h, c] <- bot_val
            }
          }
        }
        if (pad_x > 0L) {
          for (c in seq_len(n_ch)) {
            left_val <- mean(as.numeric(scaled_im[1L, , c]))
            right_val <- mean(as.numeric(scaled_im[new_w, , c]))
            if (is_raw_im) {
              arr[seq_len(pad_x), , c] <- as.raw(pmax(0L, pmin(255L, round(left_val))))
              arr[(pad_x + new_w + 1L):target_w, , c] <- as.raw(pmax(0L, pmin(255L, round(right_val))))
            } else {
              arr[seq_len(pad_x), , c] <- left_val
              arr[(pad_x + new_w + 1L):target_w, , c] <- right_val
            }
          }
        }
        return(as_image(arr, colormode = cm))
      } else {
        arr <- scaled_im[map_x, map_y, drop = FALSE]
        if (pad_y > 0L) {
          top_val <- mean(as.numeric(scaled_im[, 1L]))
          bot_val <- mean(as.numeric(scaled_im[, new_h]))
          if (is_raw_im) {
            arr[, seq_len(pad_y)] <- as.raw(pmax(0L, pmin(255L, round(top_val))))
            arr[, (pad_y + new_h + 1L):target_h] <- as.raw(pmax(0L, pmin(255L, round(bot_val))))
          } else {
            arr[, seq_len(pad_y)] <- top_val
            arr[, (pad_y + new_h + 1L):target_h] <- bot_val
          }
        }
        if (pad_x > 0L) {
          left_val <- mean(as.numeric(scaled_im[1L, ]))
          right_val <- mean(as.numeric(scaled_im[new_w, ]))
          if (is_raw_im) {
            arr[seq_len(pad_x), ] <- as.raw(pmax(0L, pmin(255L, round(left_val))))
            arr[(pad_x + new_w + 1L):target_w, ] <- as.raw(pmax(0L, pmin(255L, round(right_val))))
          } else {
            arr[seq_len(pad_x), ] <- left_val
            arr[(pad_x + new_w + 1L):target_w, ] <- right_val
          }
        }
        return(as_image(arr, colormode = "grayscale"))
      }
    }

    # pad_mode == "color" (solid color fill)
    rgb_vals <- tryCatch(grDevices::col2rgb(bg_col) / 255.0, error = function(e) c(0, 0, 0))

    if (n_ch >= 3L) {
      if (is_raw_im) {
        arr <- array(as.raw(0), dim = c(target_w, target_h, n_ch))
        bg_raw <- as.raw(pmax(0L, pmin(255L, round(rgb_vals * 255.0))))
        arr[,,1] <- bg_raw[1]
        arr[,,2] <- bg_raw[2]
        arr[,,3] <- bg_raw[3]
        if (n_ch == 4L) arr[,,4] <- as.raw(255)
      } else {
        arr <- array(0.0, dim = c(target_w, target_h, n_ch))
        arr[,,1] <- rgb_vals[1]
        arr[,,2] <- rgb_vals[2]
        arr[,,3] <- rgb_vals[3]
        if (n_ch == 4L) arr[,,4] <- 1.0
      }
      x_idx <- (pad_x + 1L):(pad_x + new_w)
      y_idx <- (pad_y + 1L):(pad_y + new_h)
      arr[x_idx, y_idx, 1:n_ch] <- scaled_im[,,1:n_ch]
      return(as_image(arr, colormode = cm))
    } else {
      if (is_raw_im) {
        bg_raw <- as.raw(pmax(0L, pmin(255L, round(mean(rgb_vals) * 255.0))))
        arr <- array(bg_raw, dim = c(target_w, target_h))
      } else {
        arr <- array(mean(rgb_vals), dim = c(target_w, target_h))
      }
      x_idx <- (pad_x + 1L):(pad_x + new_w)
      y_idx <- (pad_y + 1L):(pad_y + new_h)
      arr[x_idx, y_idx] <- scaled_im[,]
      return(as_image(arr, colormode = "grayscale"))
    }
  }

  if (mode == "crop") {
    scale <- max(target_w / orig_w, target_h / orig_h)
    scaled_w <- max(target_w, round(orig_w * scale))
    scaled_h <- max(target_h, round(orig_h * scale))
    crop_x <- (scaled_w - target_w) %/% 2
    crop_y <- (scaled_h - target_h) %/% 2

    scaled_im <- image_resize(im, width = scaled_w, height = scaled_h, filter = filter, verbose = FALSE, plot = FALSE)
    crop_x1 <- crop_x + 1L
    crop_y1 <- crop_y + 1L

    if (n_ch >= 3L) {
      cropped_arr <- scaled_im[crop_x1:(crop_x1 + target_w - 1L), crop_y1:(crop_y1 + target_h - 1L), 1:n_ch, drop = FALSE]
      return(as_image(cropped_arr, colormode = cm))
    } else {
      cropped_arr <- scaled_im[crop_x1:(crop_x1 + target_w - 1L), crop_y1:(crop_y1 + target_h - 1L), drop = FALSE]
      return(as_image(cropped_arr, colormode = "grayscale"))
    }
  }

  image_resize(im, width = target_w, height = target_h, filter = filter, verbose = FALSE, plot = FALSE)
}


#' @title Resize Images and Adjust Annotations in a YOLO Dataset
#' @name yolo_dataset_resize
#'
#' @description
#' `yolo_dataset_resize()` resizes all images in a YOLO dataset to a target
#' resolution using \code{\link{image_resize}} while automatically recalculating
#' and adapting bounding box, segmentation polygon, or OBB annotations to maintain
#' exact alignment.
#'
#' It supports different resizing strategies:
#' \itemize{
#'   \item \code{"letterbox"} (default): Resizes the image preserving its aspect
#'   ratio to fit inside the target dimensions, then centers it on a solid-colored
#'   canvas (e.g. black or gray) of the target dimensions. Bounding box and polygon
#'   coordinates are scaled and translated to match the padded canvas. Ideal for
#'   converting non-square (e.g. 16:9, 4:3) images to square (640x640, 960x960)
#'   formats for YOLO training without object distortion.
#'   \item \code{"stretch"}: Directly scales image width and height to target
#'   dimensions. Normalized coordinates (in \code{[0, 1]}) remain numerically unchanged.
#'   \item \code{"crop"}: Resizes image preserving aspect ratio to cover the
#'   target dimensions, then center-crops to target dimensions. Coordinates
#'   are shifted and clipped accordingly.
#' }
#'
#' @param data A path to a YOLO dataset directory, a \code{yolo_dataset} object,
#'   or an object returned by \code{yolo_dataset_split()}, \code{yolo_dataset_autolabel()}, etc.
#' @param out_dir Destination directory for the resized dataset. If \code{NULL}
#'   and \code{in_place = FALSE}, defaults to \code{"<data>_resized"}.
#' @param size Target size in pixels. Can be a single integer (e.g., \code{640} or \code{960})
#'   to produce square images \code{c(size, size)}, or a length-2 integer vector
#'   \code{c(width, height)}. Default is \code{640}.
#' @param width Optional integer. Target width in pixels. If provided, overrides \code{size}.
#' @param height Optional integer. Target height in pixels. If provided, overrides \code{size}.
#' @param rel_size Optional numeric. Relative scaling factor as percentage (e.g., \code{50} for 50 percent).
#' @param mode Character string specifying the resize strategy: \code{"letterbox"} (default, scale preserving aspect ratio and pad canvas), \code{"stretch"} (distort width and height to target dimensions), or \code{"crop"} (scale preserving aspect ratio and center-crop).
#' @param bg_col Color for canvas padding in \code{"letterbox"} mode when \code{pad_mode = "color"}. Default is \code{"black"}. Can also accept \code{"edge"} or \code{"mean_edge"} as a shortcut to \code{pad_mode}.
#' @param pad_mode Character string specifying how padding areas are filled in \code{"letterbox"} mode:
#' \itemize{
#'   \item \code{"edge"} (default): Replicates boundary pixels outward (clamping / border replicate), seamlessly expanding the natural background of the scene (e.g. soil, table, lab paper).
#'   \item \code{"color"}: Fills padding bars with the solid color specified by \code{bg_col} (default: \code{"black"}).
#'   \item \code{"mean_edge"}: Fills padding bars with the average color along each adjacent image border.
#' }
#' @param filter Character string specifying the interpolation filter used by \code{\link{image_resize}}.
#'   One of \code{"lanczos"} (default), \code{"bilinear"}, or \code{"nearest"}.
#' @param quality Integer between 1 and 100. JPEG compression quality when exporting resized images. Default is \code{95}.
#' @param skip_existing Logical. If \code{TRUE} (default), checks image dimensions and skips images that are already at the target width and height (copying them directly without recompression if not in-place), saving compute time and preventing generation loss.
#' @param in_place Logical. If \code{TRUE} (default), modifies images and labels directly in the source directory.
#' @param check_in_place Logical. If \code{TRUE} (default), prompts for interactive user confirmation before modifying files in-place to prevent accidental data loss. Default is \code{TRUE}.
#' @param copy_yaml Logical. Whether to copy / generate \code{data.yaml} in the output directory. Default is \code{TRUE}.
#' @param verbose Logical. Whether to show informative messages and summaries. Default is \code{TRUE}.
#' @param progress Logical. Whether to display a progress bar. Default is \code{TRUE}.
#' @param as_job Logical. If \code{TRUE} and running inside RStudio, runs the resizing task in an RStudio Background Job. Default is \code{FALSE}.
#'
#' @return A \code{yolo_dataset} object pointing to the resized dataset directory.
#'
#' @examples
#' \dontrun{
#' # 1. Resize dataset to 640x640 with standard black letterbox padding
#' ds_resized <- yolo_dataset_resize("capsulas", size = 640)
#'
#' # 2. Resize replicating border colors to simulate expanded natural background
#' ds_exp <- yolo_dataset_resize("capsulas", size = 640, pad_mode = "edge")
#'
#' # 3. Resize to 960x960 in-place
#' yolo_dataset_resize("capsulas", size = 960, in_place = TRUE)
#' }
#'
#' @aliases yolo_resize yolo_dataset_rescale
#' @export
yolo_dataset_resize <- function(data,
                                out_dir = NULL,
                                size = 640,
                                width = NULL,
                                height = NULL,
                                rel_size = NULL,
                                mode = c("letterbox", "stretch", "crop"),
                                bg_col = "black",
                                pad_mode = c("edge", "color",  "mean_edge"),
                                filter = c("lanczos", "bilinear", "nearest"),
                                quality = 95,
                                skip_existing = TRUE,
                                in_place = TRUE,
                                check_in_place = TRUE,
                                copy_yaml = TRUE,
                                verbose = TRUE,
                                progress = TRUE,
                                as_job = FALSE) {

  mode <- match.arg(mode)
  filter <- match.arg(filter)

  # Check if bg_col was used as a shortcut for pad_mode (e.g. bg_col = "edge")
  if (is.character(bg_col) && length(bg_col) == 1L) {
    low_bg <- tolower(trimws(bg_col))
    if (low_bg %in% c("edge", "border", "replicate")) {
      pad_mode <- "edge"
    } else if (low_bg %in% c("mean_edge", "edge_mean", "mean")) {
      pad_mode <- "mean_edge"
    }
  }
  pad_mode <- match.arg(pad_mode)

  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Source dataset directory {.path {dir_path}} does not exist.")
  }

  # Safety check when modifying in-place
  is_dest_inplace <- isTRUE(in_place) || (!is.null(out_dir) && identical(
    normalizePath(dir_path, winslash = "/", mustWork = FALSE),
    normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)
  ))

  if (is_dest_inplace) {
    if (!.check_in_place(dir_path, check = check_in_place)) {
      return(invisible(NULL))
    }
  }

  # Handle RStudio Background Job dispatch
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside RStudio). Executing in foreground.")
    } else {
      evaluated_args <- list(
        data = if (inherits(data, "yolo_dataset")) data$dir else data,
        out_dir = out_dir,
        size = size,
        width = width,
        height = height,
        rel_size = rel_size,
        mode = mode,
        bg_col = bg_col,
        pad_mode = pad_mode,
        filter = filter,
        quality = quality,
        skip_existing = skip_existing,
        in_place = in_place,
        check_in_place = FALSE,
        copy_yaml = copy_yaml,
        verbose = verbose,
        progress = progress,
        as_job = FALSE
      )

      job_dir <- file.path(tempdir(), "yolo_jobs")
      if (!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("yolo_resize_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("yolo_resize_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      target_name <- tryCatch({
        if (inherits(data, "yolo_dataset")) basename(data$dir)
        else if (is.character(data) && nzchar(data[1])) basename(data[1])
        else "dataset"
      }, error = function(e) "dataset")
      j_name <- paste0("YOLO Dataset Resize: ", target_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::yolo_dataset_resize()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman YOLO Dataset Resize - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::yolo_dataset_resize, job_info$args)\n",
        "cat('\\n[OK] YOLO Dataset Resize Job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("YOLO Dataset Resize (pliman)")
      cli::cli_alert_success("Dataset resizing has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Target directory: {.val {target_name}}",
        "Mode: {.val {mode}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

  # Resolve target directory
  if (isTRUE(in_place)) {
    target_dir <- dir_path
  } else if (is.null(out_dir)) {
    target_dir <- paste0(dir_path, "_resized")
  } else {
    target_dir <- normalizePath(out_dir[1], winslash = "/", mustWork = FALSE)
  }

  is_inplace <- identical(
    normalizePath(dir_path, winslash = "/", mustWork = FALSE),
    normalizePath(target_dir, winslash = "/", mustWork = FALSE)
  )

  if (!is_inplace && !dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Resolve target dimensions
  is_rel <- !is.null(rel_size)
  if (is_rel) {
    target_w <- NULL
    target_h <- NULL
    dim_desc <- sprintf("%.1f%% relative scaling", as.numeric(rel_size[1]))
  } else {
    if (!is.null(width) || !is.null(height)) {
      if (!is.null(width) && !is.null(height)) {
        target_w <- as.integer(width[1])
        target_h <- as.integer(height[1])
      } else if (!is.null(width)) {
        target_w <- as.integer(width[1])
        target_h <- as.integer(width[1])
      } else {
        target_w <- as.integer(height[1])
        target_h <- as.integer(height[1])
      }
    } else {
      if (length(size) >= 2L) {
        target_w <- as.integer(size[1])
        target_h <- as.integer(size[2])
      } else {
        target_w <- as.integer(size[1])
        target_h <- as.integer(size[1])
      }
    }
    dim_desc <- sprintf("%d x %d pixels", target_w, target_h)
  }

  detected_task <- .yolo_detect_task(dir_path)
  class_names <- .yolo_get_classes(dir_path)

  # Discover images and corresponding label paths
  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  splits_found <- character(0)
  if (dir.exists(file.path(dir_path, "images"))) {
    for (sp in c("train", "val", "test")) {
      if (dir.exists(file.path(dir_path, "images", sp))) {
        splits_found <- c(splits_found, sp)
      }
    }
  }

  records <- list()

  if (length(splits_found) > 0L) {
    # Case 1: Standard images/<split> and labels/<split>
    for (sp in splits_found) {
      sp_im_dir <- file.path(dir_path, "images", sp)
      sp_lb_dir <- file.path(dir_path, "labels", sp)
      dest_im_dir <- file.path(target_dir, "images", sp)
      dest_lb_dir <- file.path(target_dir, "labels", sp)
      dir.create(dest_im_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(dest_lb_dir, recursive = TRUE, showWarnings = FALSE)

      fls <- list.files(sp_im_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
      for (f in fls) {
        stem <- tools::file_path_sans_ext(basename(f))
        ext <- tools::file_ext(f)
        lbl_cand <- file.path(sp_lb_dir, paste0(stem, ".txt"))
        records[[length(records) + 1L]] <- list(
          src_img = f,
          src_lbl = if (file.exists(lbl_cand)) lbl_cand else NULL,
          dest_img = file.path(dest_im_dir, paste0(stem, ".", ext)),
          dest_lbl = file.path(dest_lb_dir, paste0(stem, ".txt")),
          split = sp,
          stem = stem
        )
      }
    }
  } else if (dir.exists(file.path(dir_path, "images"))) {
    # Case 2: images/ and labels/ without sub-splits
    sp_im_dir <- file.path(dir_path, "images")
    sp_lb_dir <- file.path(dir_path, "labels")
    dest_im_dir <- file.path(target_dir, "images")
    dest_lb_dir <- file.path(target_dir, "labels")
    dir.create(dest_im_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(dest_lb_dir, recursive = TRUE, showWarnings = FALSE)

    fls <- list.files(sp_im_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
    for (f in fls) {
      stem <- tools::file_path_sans_ext(basename(f))
      ext <- tools::file_ext(f)
      lbl_cand <- file.path(sp_lb_dir, paste0(stem, ".txt"))
      records[[length(records) + 1L]] <- list(
        src_img = f,
        src_lbl = if (file.exists(lbl_cand)) lbl_cand else NULL,
        dest_img = file.path(dest_im_dir, paste0(stem, ".", ext)),
        dest_lbl = file.path(dest_lb_dir, paste0(stem, ".txt")),
        split = NULL,
        stem = stem
      )
    }
  } else if (detected_task == "classify") {
    # Case 3: Classification dataset
    cls_splits <- c("train", "val", "test")
    cls_splits <- cls_splits[dir.exists(file.path(dir_path, cls_splits))]
    if (length(cls_splits) > 0L) {
      for (sp in cls_splits) {
        subdirs <- list.dirs(file.path(dir_path, sp), full.names = TRUE, recursive = FALSE)
        for (sd in subdirs) {
          cls_name <- basename(sd)
          dest_cls_dir <- file.path(target_dir, sp, cls_name)
          dir.create(dest_cls_dir, recursive = TRUE, showWarnings = FALSE)
          fls <- list.files(sd, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
          for (f in fls) {
            records[[length(records) + 1L]] <- list(
              src_img = f,
              src_lbl = NULL,
              dest_img = file.path(dest_cls_dir, basename(f)),
              dest_lbl = NULL,
              split = sp,
              stem = tools::file_path_sans_ext(basename(f))
            )
          }
        }
      }
    } else {
      subdirs <- list.dirs(dir_path, full.names = TRUE, recursive = FALSE)
      for (sd in subdirs) {
        cls_name <- basename(sd)
        dest_cls_dir <- file.path(target_dir, cls_name)
        dir.create(dest_cls_dir, recursive = TRUE, showWarnings = FALSE)
        fls <- list.files(sd, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
        for (f in fls) {
          records[[length(records) + 1L]] <- list(
            src_img = f,
            src_lbl = NULL,
            dest_img = file.path(dest_cls_dir, basename(f)),
            dest_lbl = NULL,
            split = NULL,
            stem = tools::file_path_sans_ext(basename(f))
          )
        }
      }
    }
  } else {
    # Case 4: Flat directory containing images and .txt files side-by-side
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    fls <- list.files(dir_path, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
    for (f in fls) {
      stem <- tools::file_path_sans_ext(basename(f))
      ext <- tools::file_ext(f)
      lbl_cand <- file.path(dir_path, paste0(stem, ".txt"))
      records[[length(records) + 1L]] <- list(
        src_img = f,
        src_lbl = if (file.exists(lbl_cand)) lbl_cand else NULL,
        dest_img = file.path(target_dir, paste0(stem, ".", ext)),
        dest_lbl = file.path(target_dir, paste0(stem, ".txt")),
        split = NULL,
        stem = stem
      )
    }
  }

  n_total <- length(records)
  if (n_total == 0L) {
    cli::cli_abort("No valid images found in {.path {dir_path}} matching {.val {img_pattern}}.")
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("YOLO Dataset Resize (pliman)")
    cli::cli_bullets(c(
      "*" = "Source dataset: {.path {dir_path}}",
      "*" = "Target dataset: {.path {target_dir}}{if (is_inplace) ' (in-place)' else ''}",
      "*" = "Target resolution: {.val {dim_desc}}",
      "*" = "Strategy mode: {.val {mode}}{if (mode == 'letterbox') paste0(' (pad: ', pad_mode, if (pad_mode == 'color') paste0(' [', bg_col, ']') else '', ')') else ''}",
      "*" = "Images to resize: {.val {n_total}}"
    ))
  }

  # Step 1: Fast C++ batch dimension check across all images
  all_src_imgs <- vapply(records, function(r) r$src_img, character(1))
  dim_mat <- tryCatch(
    fast_image_dimensions_cpp(all_src_imgs),
    error = function(e) matrix(0L, nrow = n_total, ncol = 2L)
  )

  # Step 2: Determine target dimensions and identify which images need resizing
  if (!is_rel) {
    target_w_vec <- rep(target_w, n_total)
    target_h_vec <- rep(target_h, n_total)
  } else {
    target_w_vec <- pmax(1L, round(dim_mat[, 1] * rel_size / 100))
    target_h_vec <- pmax(1L, round(dim_mat[, 2] * rel_size / 100))
  }

  can_skip <- isTRUE(skip_existing) &
              (dim_mat[, 1] > 0L) & (dim_mat[, 2] > 0L) &
              (dim_mat[, 1] == target_w_vec) & (dim_mat[, 2] == target_h_vec)

  skipped_indices <- which(can_skip)
  process_indices <- which(!can_skip)

  total_skipped <- length(skipped_indices)
  total_images_processed <- 0L
  total_labels_modified <- 0L

  # Step 3: Handle skipped images (fast path: no decoding, no resizing)
  if (total_skipped > 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Fast C++ check: {total_skipped} of {n_total} image{?s} already match target resolution ({dim_desc}).")
    }
    if (!is_inplace) {
      for (idx in skipped_indices) {
        rec <- records[[idx]]
        if (normalizePath(rec$src_img, winslash = "/", mustWork = FALSE) !=
            normalizePath(rec$dest_img, winslash = "/", mustWork = FALSE)) {
          file.copy(rec$src_img, rec$dest_img, overwrite = TRUE)
        }
        if (!is.null(rec$dest_lbl)) {
          if (!is.null(rec$src_lbl) && file.exists(rec$src_lbl)) {
            if (normalizePath(rec$src_lbl, winslash = "/", mustWork = FALSE) !=
                normalizePath(rec$dest_lbl, winslash = "/", mustWork = FALSE)) {
              file.copy(rec$src_lbl, rec$dest_lbl, overwrite = TRUE)
            }
          } else {
            file.create(rec$dest_lbl)
          }
        }
      }
    }
  }

  # Step 4: Resize ONLY the images that actually need it
  n_to_process <- length(process_indices)
  if (n_to_process > 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Resizing {n_to_process} image{?s}...")
    }
    show_pb <- isTRUE(verbose) && isTRUE(progress)
    if (show_pb) {
      pb <- cli::cli_progress_bar(
        name = "Resizing images & labels",
        total = n_to_process,
        format = "{cli::pb_spin} Resizing [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
        clear = FALSE
      )
    }

    for (i in process_indices) {
      rec <- records[[i]]
      im <- try(image_import(rec$src_img), silent = TRUE)
      if (inherits(im, "try-error") || !is_image(im)) {
        if (show_pb) cli::cli_progress_update(id = pb)
        next
      }

      orig_w <- dim(im)[1]
      orig_h <- dim(im)[2]

      cur_tw <- if (is_rel) max(1L, round(orig_w * rel_size / 100)) else target_w
      cur_th <- if (is_rel) max(1L, round(orig_h * rel_size / 100)) else target_h

      # Second-chance check (in case fast C++ header reader returned 0 for this file)
      if (isTRUE(skip_existing) && orig_w == cur_tw && orig_h == cur_th) {
        total_skipped <- total_skipped + 1L
        if (!is_inplace) {
          if (normalizePath(rec$src_img, winslash = "/", mustWork = FALSE) !=
              normalizePath(rec$dest_img, winslash = "/", mustWork = FALSE)) {
            file.copy(rec$src_img, rec$dest_img, overwrite = TRUE)
          }
          if (!is.null(rec$dest_lbl)) {
            if (!is.null(rec$src_lbl) && file.exists(rec$src_lbl)) {
              if (normalizePath(rec$src_lbl, winslash = "/", mustWork = FALSE) !=
                  normalizePath(rec$dest_lbl, winslash = "/", mustWork = FALSE)) {
                file.copy(rec$src_lbl, rec$dest_lbl, overwrite = TRUE)
              }
            } else {
              file.create(rec$dest_lbl)
            }
          }
        }
        if (show_pb) cli::cli_progress_update(id = pb)
        next
      }

      # Resize image
      resized_im <- .yolo_resize_single_image(
        im = im,
        target_w = cur_tw,
        target_h = cur_th,
        mode = mode,
        bg_col = bg_col,
        pad_mode = pad_mode,
        filter = filter
      )

      # Save resized image
      if (is_inplace) {
        tmp_img_file <- file.path(dirname(rec$dest_img), paste0(".tmp_res_", basename(rec$dest_img)))
        image_export(resized_im, name = tmp_img_file, quality = quality)
        file.rename(tmp_img_file, rec$dest_img)
      } else {
        image_export(resized_im, name = rec$dest_img, quality = quality)
      }
      total_images_processed <- total_images_processed + 1L

      # Transform and save label file if applicable
      if (!is.null(rec$dest_lbl)) {
        if (!is.null(rec$src_lbl) && file.exists(rec$src_lbl)) {
          raw_lines <- readLines(rec$src_lbl, warn = FALSE)
          raw_lines <- raw_lines[nzchar(trimws(raw_lines))]

          if (length(raw_lines) > 0L) {
            new_lines <- .yolo_resize_annotations(
              raw_lines = raw_lines,
              orig_w = orig_w,
              orig_h = orig_h,
              target_w = cur_tw,
              target_h = cur_th,
              mode = mode
            )
          } else {
            new_lines <- character(0)
          }

          if (is_inplace) {
            tmp_lbl_file <- file.path(dirname(rec$dest_lbl), paste0(".tmp_res_", basename(rec$dest_lbl)))
            writeLines(new_lines, tmp_lbl_file)
            file.rename(tmp_lbl_file, rec$dest_lbl)
          } else {
            writeLines(new_lines, rec$dest_lbl)
          }
          total_labels_modified <- total_labels_modified + 1L

        } else if (!is_inplace) {
          # Negative / background sample in detection dataset: create empty label
          file.create(rec$dest_lbl)
        }
      }

      if (show_pb) cli::cli_progress_update(id = pb)
    }

    if (show_pb) cli::cli_progress_done(id = pb)
  }

  # Invalidate stale label caches so YOLO training always re-reads updated coordinates
  cache_files <- list.files(target_dir, pattern = "\\.cache$", recursive = TRUE, full.names = TRUE)
  if (length(cache_files) > 0L) {
    unlink(cache_files, force = TRUE)
  }

  # Copy / generate data.yaml if not in-place
  if (!is_inplace && isTRUE(copy_yaml)) {
    src_yaml <- file.path(dir_path, "data.yaml")
    dst_yaml <- file.path(target_dir, "data.yaml")
    if (file.exists(src_yaml)) {
      yaml_lines <- readLines(src_yaml, warn = FALSE)
      path_idx <- grep("^path\\s*:", yaml_lines)
      if (length(path_idx) > 0L) {
        yaml_lines[path_idx[1]] <- paste0("path: ", normalizePath(target_dir, winslash = "/", mustWork = FALSE))
      }
      writeLines(yaml_lines, dst_yaml)
    } else if (length(splits_found) > 0L) {
      tr_d <- if ("train" %in% splits_found) "images/train" else "images"
      val_d <- if ("val" %in% splits_found) "images/val" else tr_d
      tst_d <- if ("test" %in% splits_found) "images/test" else NULL
      .yolo_write_yaml(
        out_dir = target_dir,
        classes = class_names,
        train_dir = tr_d,
        val_dir = val_d,
        test_dir = tst_d,
        task = detected_task
      )
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Dataset resizing complete! Processed {total_images_processed} image{?s}.")
    cli::cli_bullets(c(
      "*" = "Images resized: {total_images_processed}",
      "*" = "Images skipped (already target size): {total_skipped}",
      "*" = "Label files transformed: {total_labels_modified}",
      "*" = "Output dataset directory: {.path {target_dir}}",
      "v" = "Label caches invalidated."
    ))
  }

  # Wrap and return yolo_dataset object
  out_wrap <- tryCatch(
    .yolo_wrap_dataset(target_dir, task = detected_task, classes = class_names),
    error = function(e) list(dir = target_dir, task = detected_task, classes = class_names)
  )
  class(out_wrap) <- unique(c("yolo_dataset_resized", class(out_wrap)))
  out_wrap$is_resized <- TRUE
  out_wrap$target_size <- if (is_rel) rel_size else c(target_w, target_h)
  out_wrap$resize_mode <- mode
  out_wrap$pad_mode <- pad_mode
  out_wrap$n_images_resized <- total_images_processed
  out_wrap$n_images_skipped <- total_skipped
  out_wrap$n_labels_modified <- total_labels_modified

  invisible(out_wrap)
}

#' @rdname yolo_dataset_resize
#' @export
yolo_resize <- yolo_dataset_resize

#' @rdname yolo_dataset_resize
#' @export
yolo_dataset_rescale <- yolo_dataset_resize


# ==============================================================================
# 13. COMPREHENSIVE DATASET PROFILING & DESCRIPTION (yolo_dataset_describe)
# ==============================================================================
#' @title Comprehensive YOLO Dataset Inspection, Statistics, and Profiling
#' @name yolo_dataset_describe
#' @aliases yolo_describe yolo_dataset_summary
#'
#' @description
#' `yolo_dataset_describe()` computes a comprehensive, deep statistical profile
#' of any YOLO dataset (Object Detection, Instance Segmentation, Oriented Bounding Box,
#' or Image Classification).
#'
#' It profiles:
#' \itemize{
#'   \item **Splits allocation**: Exact image counts and percentages across \code{train},
#'   \code{val}, and \code{test} splits, annotated vs. negative/background images,
#'   and mean/max object density per image.
#'   \item **Class distribution**: Granular object frequencies per class across splits,
#'   class balance/imbalance ratios, image co-occurrence frequencies, and average spatial
#'   object area percentages.
#'   \item **Image properties**: Sampled resolutions (e.g. 640x640), aspect ratio distributions,
#'   file formats, and total disk footprint in MB.
#'   \item **Key insights**: Diagnostic alerts regarding class imbalance, background sample
#'   adequacy, and split ratios.
#' }
#'
#' @param data A path to a YOLO dataset directory, a \code{yolo_dataset} object,
#'   or an object returned by \code{\link{yolo_dataset_split}}, \code{\link{yolo_dataset_autolabel}},
#'   \code{\link{yolo_dataset_resize}}, \code{\link{yolo_dataset_audit}}, etc.
#' @param by_class Logical. Whether to compute granular per-class statistics and spatial
#'   area distributions. Default is \code{TRUE}.
#' @param by_split Logical. Whether to compute metrics broken down by split (\code{train},
#'   \code{val}, \code{test}). Default is \code{TRUE}.
#' @param sample_dims Integer. Number of images to sample to profile resolutions and aspect
#'   ratios without full-decoding delays on large datasets. Default is \code{50}.
#'   Set to \code{Inf} or \code{NULL} to inspect all images.
#' @param plot Logical. If \code{TRUE}, automatically displays visual distribution plots.
#'   Default is \code{FALSE}.
#' @param verbose Logical. Whether to print formatted CLI summary tables and insights
#'   to the console. Default is \code{TRUE}.
#'
#' @return A \code{yolo_dataset_desc} object (class \code{c("yolo_dataset_desc", "list")})
#' containing structured summary tables:
#' \itemize{
#'   \item \code{dir}: Dataset path.
#'   \item \code{task}: Detected YOLO task (\code{"detect"}, \code{"segment"}, \code{"obb"}, or \code{"classify"}).
#'   \item \code{classes}: Character vector of class names.
#'   \item \code{total_images}: Total image count across all splits.
#'   \item \code{total_objects}: Total labeled objects across all classes.
#'   \item \code{total_annotated}: Count of images with at least one object.
#'   \item \code{total_background}: Count of negative/background images (0 annotations).
#'   \item \code{total_original}: Count of canonical original images.
#'   \item \code{total_augmented}: Count of offline-augmented images (\code{_aug}).
#'   \item \code{total_synthetic}: Count of synthetic background images (\code{synth_bg_}).
#'   \item \code{total_tiled}: Count of sliced/tiled crop images (\code{_tile_}).
#'   \item \code{provenance_summary}: A \code{data.frame} detailing image counts and percentages by origin (Original, Augmented, Synthetic, Tiled) across splits.
#'   \item \code{splits_summary}: A \code{data.frame} summarizing images, background ratios, and object counts per split.
#'   \item \code{classes_summary}: A \code{data.frame} detailing object frequencies, proportions, image occurrences, and spatial sizes per class.
#'   \item \code{density_summary}: A \code{data.frame} categorizing images into object count density buckets.
#'   \item \code{scale_summary}: A \code{data.frame} breaking down objects by spatial area into COCO benchmark categories (Small, Medium, Large).
#'   \item \code{formats_summary}: A \code{data.frame} of image file extensions found in the dataset.
#'   \item \code{imbalance_ratio}: Numerical ratio between most and least frequent classes.
#'   \item \code{disk_mb}: Total disk footprint in megabytes.
#'   \item \code{avg_file_kb}: Average image file size in kilobytes.
#' }
#'
#' @examples
#' \dontrun{
#' # Profile a YOLO dataset
#' desc <- yolo_dataset_describe("capsulas")
#'
#' # Generate visual summary charts
#' plot(desc)
#' plot(desc, which = "provenance")
#' plot(desc, which = "density")
#' }
#'
#' @export
yolo_dataset_describe <- function(data,
                                  by_class = TRUE,
                                  by_split = TRUE,
                                  sample_dims = 50,
                                  plot = FALSE,
                                  verbose = TRUE) {
  dir_path <- .yolo_get_dir(data)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Source dataset directory {.path {dir_path}} does not exist.")
  }

  task <- .yolo_detect_task(dir_path)
  classes <- .yolo_get_classes(dir_path)
  if (is.null(classes) || length(classes) == 0L) {
    classes <- "object"
  }
  n_classes <- length(classes)

  yaml_path <- file.path(dir_path, "data.yaml")
  has_yaml <- file.exists(yaml_path)

  # Check disk footprint using fast C++ directory recursion
  total_bytes <- fast_dir_size_cpp(dir_path)
  disk_mb <- round(total_bytes / (1024 * 1024), 2)

  img_pattern <- "\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  max_sample <- if (is.null(sample_dims) || is.infinite(sample_dims)) .Machine$integer.max else as.integer(sample_dims[1])

  # ----------------------------------------------------------------------------
  # CASE 1: Classification datasets
  # ----------------------------------------------------------------------------
  if (task == "classify") {
    cls_splits <- c("train", "val", "test")
    cls_splits <- cls_splits[dir.exists(file.path(dir_path, cls_splits))]
    is_split <- length(cls_splits) > 0L

    if (!is_split) {
      cls_splits <- "all"
    }

    class_split_counts <- matrix(0L, nrow = n_classes, ncol = length(cls_splits),
                                 dimnames = list(classes, cls_splits))
    split_img_counts <- setNames(integer(length(cls_splits)), cls_splits)
    img_dims <- list()
    img_exts <- character(0)

    prov_types <- c("Original", "Augmented", "Synthetic_BG", "Tiled")
    class_split_prov <- matrix(0L, nrow = length(prov_types), ncol = length(cls_splits),
                               dimnames = list(prov_types, cls_splits))
    all_img_exts <- character(0)

    split_records <- list()
    for (sp in cls_splits) {
      sp_dir <- if (sp == "all") dir_path else file.path(dir_path, sp)
      recs <- list()
      for (cls in classes) {
        c_dir <- file.path(sp_dir, cls)
        if (dir.exists(c_dir)) {
          fls <- list.files(c_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
          class_split_counts[cls, sp] <- length(fls)
          split_img_counts[sp] <- split_img_counts[sp] + length(fls)
          for (f in fls) {
            stem <- tools::file_path_sans_ext(basename(f))
            is_aug <- grepl("_aug\\d+", stem, ignore.case = TRUE)
            is_synth <- grepl("(^synth_bg_|^synth_|_synth_|_synth$)", stem, ignore.case = TRUE)
            is_tile <- grepl("(_tile_\\d+_\\d+|_tile_)", stem, ignore.case = TRUE)
            ptype <- if (is_aug) "Augmented" else if (is_synth) "Synthetic_BG" else if (is_tile) "Tiled" else "Original"
            class_split_prov[ptype, sp] <- class_split_prov[ptype, sp] + 1L
            all_img_exts <- c(all_img_exts, tolower(tools::file_ext(f)))
            recs[[length(recs) + 1L]] <- list(
              img = f,
              lbl = NULL,
              stem = stem,
              split = sp,
              ptype = ptype,
              class = cls,
              n_objs = 1L,
              is_background = FALSE
            )
          }
          if (length(fls) > 0L && length(img_dims) < max_sample) {
            needed <- max(1L, max_sample - length(img_dims))
            sample_sub <- utils::head(fls, needed)
            d_mat <- fast_image_dimensions_cpp(sample_sub)
            v_idx <- which(d_mat[, 1] > 0 & d_mat[, 2] > 0)
            if (length(v_idx) > 0L) {
              for (vi in v_idx) {
                img_dims[[length(img_dims) + 1L]] <- c(d_mat[vi, 1], d_mat[vi, 2])
                img_exts <- c(img_exts, tolower(tools::file_ext(sample_sub[vi])))
              }
            }
          }
        }
      }
      split_records[[sp]] <- recs
    }

    total_images <- sum(split_img_counts)
    if (total_images == 0L) {
      cli::cli_abort(c(
        "x" = "No images found in classification dataset directory {.path {dir_path}}.",
        "i" = "Please verify that class subfolders exist and contain supported image files (jpg, jpeg, png, bmp, webp, tif, tiff)."
      ))
    }

    # Splits summary
    splits_df <- data.frame(
      Split = cls_splits,
      Images = as.integer(split_img_counts),
      Pct_Images = if (total_images > 0L) round(as.integer(split_img_counts) / total_images * 100, 1) else 0.0,
      Classes = apply(class_split_counts > 0L, 2, sum),
      stringsAsFactors = FALSE
    )
    if (nrow(splits_df) > 1L) {
      splits_df <- rbind(splits_df, data.frame(
        Split = "TOTAL",
        Images = total_images,
        Pct_Images = 100.0,
        Classes = sum(apply(class_split_counts > 0L, 1, any)),
        stringsAsFactors = FALSE
      ))
    }

    # Classes summary
    tot_per_class <- rowSums(class_split_counts)
    classes_df <- data.frame(
      ID = seq_along(classes) - 1L,
      Class = classes,
      class_split_counts,
      Total = as.integer(tot_per_class),
      Pct_Total = if (total_images > 0L) round(tot_per_class / total_images * 100, 1) else 0.0,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(classes_df) <- NULL

    # Provenance summary
    tot_per_prov <- rowSums(class_split_prov)
    prov_df <- data.frame(
      Type = prov_types,
      class_split_prov,
      Total = as.integer(tot_per_prov),
      Pct = if (total_images > 0L) round(tot_per_prov / total_images * 100, 1) else 0.0,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(prov_df) <- NULL
    if (nrow(prov_df) > 1L) {
      prov_df <- rbind(prov_df, data.frame(
        Type = "TOTAL",
        as.list(colSums(class_split_prov)),
        Total = total_images,
        Pct = 100.0,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }

    # Formats summary
    fmt_tbl <- table(all_img_exts)
    formats_df <- if (length(fmt_tbl) > 0L) {
      data.frame(
        Format = toupper(names(fmt_tbl)),
        Count = as.integer(fmt_tbl),
        Pct = if (total_images > 0L) round(as.integer(fmt_tbl) / total_images * 100, 1) else rep(0.0, length(fmt_tbl)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Format = character(0),
        Count = integer(0),
        Pct = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    pos_counts <- tot_per_class[tot_per_class > 0L]
    imbalance_ratio <- if (length(pos_counts) > 1L && min(pos_counts) > 0L) {
      round(max(pos_counts) / min(pos_counts), 2)
    } else {
      1.0
    }

    res <- list(
      dir = dir_path,
      task = "classify",
      yaml_file = if (has_yaml) yaml_path else NULL,
      classes = classes,
      splits = cls_splits,
      total_images = total_images,
      total_objects = total_images,
      total_annotated = total_images,
      total_background = 0L,
      pct_background = 0.0,
      total_original = as.integer(tot_per_prov["Original"]),
      total_augmented = as.integer(tot_per_prov["Augmented"]),
      total_synthetic = as.integer(tot_per_prov["Synthetic_BG"]),
      total_tiled = as.integer(tot_per_prov["Tiled"]),
      avg_file_kb = if (total_images > 0L) round((total_bytes / 1024) / total_images, 1) else 0.0,
      mean_density = 1.0,
      median_density = 1.0,
      min_density = 1L,
      max_density = 1L,
      n_boxes = 0L,
      n_polygons = 0L,
      provenance_summary = prov_df,
      formats_summary = formats_df,
      density_summary = NULL,
      scale_summary = NULL,
      splits_summary = splits_df,
      classes_summary = classes_df,
      imbalance_ratio = imbalance_ratio,
      disk_mb = disk_mb,
      sample_img_dims = img_dims,
      sample_img_exts = img_exts,
      split_records = split_records
    )
    class(res) <- c("yolo_dataset_desc", "list")

    if (isTRUE(verbose)) print(res)
    if (isTRUE(plot)) plot(res)
    return(invisible(res))
  }

  # ----------------------------------------------------------------------------
  # CASE 2: Detection, Segmentation, or OBB
  # ----------------------------------------------------------------------------
  std_splits <- c("train", "val", "test")
  found_splits <- character(0)
  has_split_dirs <- FALSE

  # Check standard images/<split> and labels/<split>
  if (dir.exists(file.path(dir_path, "images"))) {
    for (sp in std_splits) {
      if (dir.exists(file.path(dir_path, "images", sp))) {
        found_splits <- c(found_splits, sp)
      }
    }
    if (length(found_splits) > 0L) has_split_dirs <- TRUE
  }
  # Check <split>/images
  if (!has_split_dirs) {
    for (sp in std_splits) {
      if (dir.exists(file.path(dir_path, sp, "images")) || dir.exists(file.path(dir_path, sp))) {
        found_splits <- c(found_splits, sp)
      }
    }
    if (length(found_splits) > 0L) has_split_dirs <- TRUE
  }

  if (!has_split_dirs) {
    found_splits <- "all"
  }

  prov_types <- c("Original", "Augmented", "Synthetic_BG", "Tiled")
  class_split_prov <- matrix(0L, nrow = length(prov_types), ncol = length(found_splits),
                             dimnames = list(prov_types, found_splits))

  # Discover images and labels per split
  split_records <- list()
  for (sp in found_splits) {
    im_dir <- if (sp == "all") {
      if (dir.exists(file.path(dir_path, "images"))) file.path(dir_path, "images") else dir_path
    } else {
      if (dir.exists(file.path(dir_path, "images", sp))) file.path(dir_path, "images", sp)
      else if (dir.exists(file.path(dir_path, sp, "images"))) file.path(dir_path, sp, "images")
      else file.path(dir_path, sp)
    }

    lb_dir <- if (sp == "all") {
      if (dir.exists(file.path(dir_path, "labels"))) file.path(dir_path, "labels") else dir_path
    } else {
      if (dir.exists(file.path(dir_path, "labels", sp))) file.path(dir_path, "labels", sp)
      else if (dir.exists(file.path(dir_path, sp, "labels"))) file.path(dir_path, sp, "labels")
      else file.path(dir_path, sp)
    }

    fls <- list.files(im_dir, pattern = img_pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(fls) == 0L) {
      split_records[[sp]] <- list()
      next
    }
    stems <- tools::file_path_sans_ext(basename(fls))
    lbl_cands <- file.path(lb_dir, paste0(stems, ".txt"))
    lbl_exists <- file.exists(lbl_cands)

    is_aug <- grepl("_aug\\d+", stems, ignore.case = TRUE)
    is_synth <- grepl("(^synth_bg_|^synth_|_synth_|_synth$)", stems, ignore.case = TRUE)
    is_tile <- grepl("(_tile_\\d+_\\d+|_tile_)", stems, ignore.case = TRUE)
    ptypes <- ifelse(is_aug, "Augmented", ifelse(is_synth, "Synthetic_BG", ifelse(is_tile, "Tiled", "Original")))

    tbl_pt <- table(ptypes)
    for (pt_name in names(tbl_pt)) {
      if (pt_name %in% prov_types) {
        class_split_prov[pt_name, sp] <- as.integer(tbl_pt[pt_name])
      }
    }

    recs <- vector("list", length(fls))
    for (idx in seq_along(fls)) {
      recs[[idx]] <- list(
        img = fls[idx],
        lbl = if (lbl_exists[idx]) lbl_cands[idx] else NULL,
        stem = stems[idx],
        split = sp,
        ptype = ptypes[idx]
      )
    }
    split_records[[sp]] <- recs
  }

  total_imgs_all <- sum(vapply(split_records, length, integer(1)))
  if (total_imgs_all == 0L) {
    cli::cli_abort(c(
      "x" = "No images found in dataset directory {.path {dir_path}}.",
      "i" = "Searched for extensions ({.val jpg, jpeg, png, bmp, webp, tif, tiff}) in {.path {dir_path}} and subdirectories ({.path images/}, {.path images/train/}, etc.).",
      "i" = "Please check that the path points to a valid YOLO dataset directory or export a dataset first using {.fun yolo_dataset_export}."
    ))
  }

  all_recs <- unlist(split_records, recursive = FALSE)
  flat_img_files <- vapply(all_recs, function(r) r$img, character(1))
  flat_stems <- vapply(all_recs, function(r) r$stem, character(1))
  flat_splits <- vapply(all_recs, function(r) r$split, character(1))
  flat_ptypes <- vapply(all_recs, function(r) r$ptype, character(1))
  flat_lbl_files <- vapply(all_recs, function(r) if (!is.null(r$lbl)) r$lbl else "", character(1))
  all_img_exts <- tolower(tools::file_ext(flat_img_files))

  # Fast image dimensions via STB image header inspection (zero pixel decoding)
  sample_count <- min(total_imgs_all, max_sample)
  sample_files <- flat_img_files[seq_len(sample_count)]
  dim_mat <- fast_image_dimensions_cpp(sample_files)
  valid_dims <- which(dim_mat[, 1] > 0 & dim_mat[, 2] > 0)
  sample_img_dims <- lapply(valid_dims, function(idx) c(dim_mat[idx, 1], dim_mat[idx, 2]))
  sample_img_exts <- all_img_exts[valid_dims]

  # Fast parallel OpenMP C++ execution across all label files
  cpp_res <- parse_yolo_labels_cpp(
    lbl_files = flat_lbl_files,
    img_files = flat_img_files,
    stems = flat_stems,
    splits = flat_splits,
    ptypes = flat_ptypes,
    classes = classes,
    split_keys = found_splits,
    max_anomalies = 500L,
    num_threads = 0L
  )

  class_split_objs <- cpp_res$class_split_objs
  dimnames(class_split_objs) <- list(classes, found_splits)
  class_img_counts <- setNames(as.integer(cpp_res$class_img_counts), classes)
  class_box_areas <- setNames(cpp_res$class_box_areas, classes)
  all_box_areas <- as.numeric(cpp_res$all_box_areas)
  all_box_ars <- as.numeric(cpp_res$all_box_ars)
  all_objs_per_img <- as.integer(cpp_res$file_n_objs)
  n_boxes_total <- as.integer(cpp_res$n_boxes_total)
  n_polys_total <- as.integer(cpp_res$n_polys_total)
  spatial_grid <- cpp_res$spatial_grid
  dim_grid <- cpp_res$dim_grid
  class_cooccurrence <- cpp_res$class_cooccurrence
  dimnames(class_cooccurrence) <- list(classes, classes)
  anomaly_records <- cpp_res$anomalies_records
  n_oob_total <- as.integer(cpp_res$n_oob)
  n_zero_area_total <- as.integer(cpp_res$n_zero_area)
  n_tiny_total <- as.integer(cpp_res$n_tiny)
  n_huge_total <- as.integer(cpp_res$n_huge)
  n_extreme_ar_total <- as.integer(cpp_res$n_extreme_ar)
  n_duplicate_pairs_total <- as.integer(cpp_res$n_duplicate_pairs)
  n_invalid_cids_total <- as.integer(cpp_res$n_invalid_cids)

  # Re-assign file_n_objs and is_background back to split_records & compute split statistics
  offset <- 0L
  split_stats <- list()
  for (sp in found_splits) {
    n_sp_imgs <- length(split_records[[sp]])
    sp_indices <- seq(offset + 1L, offset + n_sp_imgs)
    sp_n_objs <- all_objs_per_img[sp_indices]
    sp_is_bg <- cpp_res$file_is_bg[sp_indices]

    for (k in seq_len(n_sp_imgs)) {
      split_records[[sp]][[k]]$n_objs <- sp_n_objs[k]
      split_records[[sp]][[k]]$is_background <- sp_is_bg[k]
    }

    n_sp_annotated <- sum(!sp_is_bg)
    n_sp_background <- sum(sp_is_bg)
    tot_sp_objs <- sum(sp_n_objs)

    split_stats[[sp]] <- list(
      Split = sp,
      Images = n_sp_imgs,
      Pct_Images = if (total_imgs_all > 0L) round(n_sp_imgs / total_imgs_all * 100, 1) else 0.0,
      Annotated = n_sp_annotated,
      Background = n_sp_background,
      Pct_Bg = if (n_sp_imgs > 0L) round(n_sp_background / n_sp_imgs * 100, 1) else 0.0,
      Objects = tot_sp_objs,
      Mean_Objs = if (n_sp_imgs > 0L) round(mean(sp_n_objs), 2) else 0.0,
      Max_Objs = if (n_sp_imgs > 0L) max(sp_n_objs) else 0L
    )
    offset <- offset + n_sp_imgs
  }

  splits_df <- do.call(rbind, lapply(split_stats, as.data.frame, stringsAsFactors = FALSE))
  rownames(splits_df) <- NULL

  if (nrow(splits_df) > 1L) {
    all_objs <- sum(splits_df$Objects)
    all_ann <- sum(splits_df$Annotated)
    all_bg <- sum(splits_df$Background)
    tot_row <- data.frame(
      Split = "TOTAL",
      Images = total_imgs_all,
      Pct_Images = 100.0,
      Annotated = all_ann,
      Background = all_bg,
      Pct_Bg = if (total_imgs_all > 0L) round(all_bg / total_imgs_all * 100, 1) else 0.0,
      Objects = all_objs,
      Mean_Objs = if (total_imgs_all > 0L) round(all_objs / total_imgs_all, 2) else 0.0,
      Max_Objs = max(splits_df$Max_Objs),
      stringsAsFactors = FALSE
    )
    splits_df <- rbind(splits_df, tot_row)
  }

  tot_objs_per_class <- rowSums(class_split_objs)
  grand_total_objs <- sum(tot_objs_per_class)

  mean_areas <- vapply(classes, function(cn) {
    a <- class_box_areas[[cn]]
    if (length(a) > 0L) round(mean(a), 2) else NA_real_
  }, numeric(1))

  size_cats <- vapply(mean_areas, function(a) {
    if (is.na(a)) return("N/A")
    if (a < 0.5) return("Tiny (<0.5%)")
    if (a < 2.0) return("Small (0.5-2%)")
    if (a < 10.0) return("Medium (2-10%)")
    return("Large (>10%)")
  }, character(1))

  classes_df <- data.frame(
    ID = seq_along(classes) - 1L,
    Class = classes,
    class_split_objs,
    Total_Objs = as.integer(tot_objs_per_class),
    Pct_Objs = if (grand_total_objs > 0L) round(tot_objs_per_class / grand_total_objs * 100, 1) else 0.0,
    Images = as.integer(class_img_counts),
    Pct_Images = if (total_imgs_all > 0L) round(class_img_counts / total_imgs_all * 100, 1) else 0.0,
    Avg_Area_Pct = mean_areas,
    Size_Cat = size_cats,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(classes_df) <- NULL

  # Provenance summary table
  tot_per_prov <- rowSums(class_split_prov)
  prov_df <- data.frame(
    Type = prov_types,
    class_split_prov,
    Total = as.integer(tot_per_prov),
    Pct = if (total_imgs_all > 0L) round(tot_per_prov / total_imgs_all * 100, 1) else 0.0,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(prov_df) <- NULL
  if (nrow(prov_df) > 1L) {
    prov_df <- rbind(prov_df, data.frame(
      Type = "TOTAL",
      as.list(colSums(class_split_prov)),
      Total = total_imgs_all,
      Pct = 100.0,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  # Formats summary table
  fmt_tbl <- table(all_img_exts)
  formats_df <- if (length(fmt_tbl) > 0L) {
    data.frame(
      Format = toupper(names(fmt_tbl)),
      Count = as.integer(fmt_tbl),
      Pct = if (total_imgs_all > 0L) round(as.integer(fmt_tbl) / total_imgs_all * 100, 1) else rep(0.0, length(fmt_tbl)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Format = character(0),
      Count = integer(0),
      Pct = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Density summary table
  b_0 <- sum(all_objs_per_img == 0L)
  b_1 <- sum(all_objs_per_img == 1L)
  b_2_5 <- sum(all_objs_per_img >= 2L & all_objs_per_img <= 5L)
  b_6_15 <- sum(all_objs_per_img >= 6L & all_objs_per_img <= 15L)
  b_16_30 <- sum(all_objs_per_img >= 16L & all_objs_per_img <= 30L)
  b_gt30 <- sum(all_objs_per_img > 30L)

  density_df <- data.frame(
    Bucket = c("0 (Background)", "1 object", "2 - 5 objects", "6 - 15 objects", "16 - 30 objects", "> 30 objects"),
    Images = c(b_0, b_1, b_2_5, b_6_15, b_16_30, b_gt30),
    Pct_Images = if (total_imgs_all > 0L) round(c(b_0, b_1, b_2_5, b_6_15, b_16_30, b_gt30) / total_imgs_all * 100, 1) else 0.0,
    stringsAsFactors = FALSE
  )
  density_df <- density_df[density_df$Images > 0L, , drop = FALSE]

  # Scale summary table (COCO benchmark standard)
  tot_boxes <- length(all_box_areas)
  scale_df <- if (tot_boxes > 0L) {
    n_sm <- sum(all_box_areas < 0.5)
    n_med <- sum(all_box_areas >= 0.5 & all_box_areas < 2.0)
    n_lg <- sum(all_box_areas >= 2.0)
    data.frame(
      Scale = c("Small (< 0.5% area)", "Medium (0.5% - 2.0% area)", "Large (> 2.0% area)"),
      Objects = c(n_sm, n_med, n_lg),
      Pct_Objects = round(c(n_sm, n_med, n_lg) / tot_boxes * 100, 1),
      stringsAsFactors = FALSE
    )
  } else NULL

  pos_counts <- tot_objs_per_class[tot_objs_per_class > 0L]
  imbalance_ratio <- if (length(pos_counts) > 1L && min(pos_counts) > 0L) {
    round(max(pos_counts) / min(pos_counts), 2)
  } else {
    1.0
  }

  total_bg_count <- sum(splits_df$Background[splits_df$Split != "TOTAL"])

  out <- list(
    dir = dir_path,
    task = task,
    yaml_file = if (has_yaml) yaml_path else NULL,
    classes = classes,
    splits = found_splits,
    total_images = total_imgs_all,
    total_objects = grand_total_objs,
    total_annotated = sum(splits_df$Annotated[splits_df$Split != "TOTAL"]),
    total_background = total_bg_count,
    pct_background = if (total_imgs_all > 0L) round(total_bg_count / total_imgs_all * 100, 1) else 0.0,
    total_original = as.integer(tot_per_prov["Original"]),
    total_augmented = as.integer(tot_per_prov["Augmented"]),
    total_synthetic = as.integer(tot_per_prov["Synthetic_BG"]),
    total_tiled = as.integer(tot_per_prov["Tiled"]),
    avg_file_kb = if (total_imgs_all > 0L) round((total_bytes / 1024) / total_imgs_all, 1) else 0.0,
    mean_density = if (length(all_objs_per_img) > 0L) round(mean(all_objs_per_img), 2) else 0.0,
    median_density = if (length(all_objs_per_img) > 0L) stats::median(all_objs_per_img) else 0.0,
    min_density = if (length(all_objs_per_img) > 0L) min(all_objs_per_img) else 0L,
    max_density = if (length(all_objs_per_img) > 0L) max(all_objs_per_img) else 0L,
    n_boxes = n_boxes_total,
    n_polygons = n_polys_total,
    provenance_summary = prov_df,
    formats_summary = formats_df,
    density_summary = density_df,
    scale_summary = scale_df,
    aspect_ratio_summary = if (length(all_box_ars) > 0L) {
      ar_tall <- sum(all_box_ars < 0.75)
      ar_square <- sum(all_box_ars >= 0.75 & all_box_ars <= 1.33)
      ar_wide <- sum(all_box_ars > 1.33)
      tot_ar <- length(all_box_ars)
      data.frame(
        Category = c("Portrait / Tall (AR < 0.75)", "Square-ish (0.75 - 1.33)", "Landscape / Wide (AR > 1.33)"),
        Objects = c(ar_tall, ar_square, ar_wide),
        Pct = if (tot_ar > 0L) round(c(ar_tall, ar_square, ar_wide) / tot_ar * 100, 1) else c(0, 0, 0),
        stringsAsFactors = FALSE
      )
    } else NULL,
    spatial_density = spatial_grid,
    box_dim_density = dim_grid,
    class_cooccurrence = if (n_classes > 1L) class_cooccurrence else NULL,
    anomalies = list(
      out_of_bounds = n_oob_total,
      zero_area = n_zero_area_total,
      micro_boxes = n_tiny_total,
      huge_boxes = n_huge_total,
      extreme_aspect_ratios = n_extreme_ar_total,
      duplicate_pairs = n_duplicate_pairs_total,
      invalid_class_ids = n_invalid_cids_total,
      records = anomaly_records
    ),
    splits_summary = splits_df,
    classes_summary = classes_df,
    imbalance_ratio = imbalance_ratio,
    disk_mb = disk_mb,
    sample_img_dims = sample_img_dims,
    sample_img_exts = sample_img_exts,
    split_records = split_records
  )
  class(out) <- c("yolo_dataset_desc", "list")

  if (isTRUE(verbose)) {
    print(out)
  }
  if (isTRUE(plot)) {
    plot(out)
  }
  invisible(out)
}

#' @rdname yolo_dataset_describe
#' @export
yolo_describe <- yolo_dataset_describe

#' @rdname yolo_dataset_describe
#' @export
yolo_dataset_summary <- yolo_dataset_describe


#' @export
print.yolo_dataset_desc <- function(x, ...) {
  cli::cli_h1("YOLO Dataset Description (pliman)")

  task_tag <- switch(x$task,
    detect = "Object Detection (BBoxes)",
    segment = "Instance Segmentation (Polygons)",
    obb = "Oriented Bounding Box (OBB)",
    classify = "Image Classification (Classes)",
    x$task
  )

  cli::cli_bullets(c(
    "*" = "Directory: {.path {x$dir}}",
    "*" = "Task: {.strong {task_tag}}",
    "*" = "Config: {if (!is.null(x$yaml_file)) cli::format_inline('{.path {x$yaml_file}}') else 'No data.yaml found'}",
    "*" = "Total Images: {.val {x$total_images}} ({x$disk_mb} MB on disk{if (!is.null(x$avg_file_kb)) paste0(', avg ~', x$avg_file_kb, ' KB/img') else ''})",
    "*" = "Total Objects: {.val {x$total_objects}} across {.val {length(x$classes)}} class{?es}",
    "*" = "Provenance: {.val {x$total_original}} original | {.val {x$total_augmented}} augmented | {.val {x$total_synthetic}} synthetic | {.val {x$total_tiled}} tiled"
  ))

  if (!is.null(x$provenance_summary)) {
    cat("\n")
    cli::cli_h2("Dataset Provenance & Composition")
    print(x$provenance_summary, row.names = FALSE)
  }

  cat("\n")
  cli::cli_h2("Splits Overview")
  print(x$splits_summary, row.names = FALSE)

  cat("\n")
  cli::cli_h2("Class Distribution")
  print(x$classes_summary, row.names = FALSE)

  if (x$task != "classify" && (!is.null(x$density_summary) || !is.null(x$scale_summary))) {
    cat("\n")
    cli::cli_h2("Object Density & Scale Profile")

    if (!is.null(x$density_summary) && nrow(x$density_summary) > 0L) {
      cli::cli_alert_info("Objects per image: mean = {.val {x$mean_density}}, median = {.val {x$median_density}} (min: {.val {x$min_density}}, max: {.val {x$max_density}})")
      print(x$density_summary, row.names = FALSE)
    }

    if (!is.null(x$scale_summary) && nrow(x$scale_summary) > 0L) {
      cat("\n")
      cli::cli_alert_info("Spatial scale breakdown (COCO benchmark):")
      print(x$scale_summary, row.names = FALSE)
    }

    if (!is.null(x$n_boxes) && (x$n_boxes + x$n_polygons) > 0L) {
      tot_geom <- x$n_boxes + x$n_polygons
      cli::cli_bullets(c(
        "*" = "Annotation geometry: {.val {x$n_boxes}} Bounding Box{?es} ({round(x$n_boxes / tot_geom * 100, 1)}%) | {.val {x$n_polygons}} Polygon{?s} ({round(x$n_polygons / tot_geom * 100, 1)}%)"
      ))
    }
  }

  if (!is.null(x$formats_summary) && nrow(x$formats_summary) > 0L) {
    cat("\n")
    cli::cli_h2("Image Specifications")
    fmt_str <- paste(sprintf("%s: %d (%.1f%%)", x$formats_summary$Format, x$formats_summary$Count, x$formats_summary$Pct), collapse = ", ")
    cli::cli_bullets(c("*" = "File formats: {fmt_str}"))
    if (length(x$sample_img_dims) > 0L) {
      dim_strings <- vapply(x$sample_img_dims, function(d) sprintf("%d x %d", d[1], d[2]), character(1))
      dim_tbl <- table(dim_strings)
      top_dim <- names(dim_tbl)[which.max(dim_tbl)]
      top_dim_pct <- round(max(dim_tbl) / length(dim_strings) * 100, 1)
      cli::cli_bullets(c("*" = "Dominant resolution: {.val {top_dim}} ({top_dim_pct}% of inspected sample)"))
    }
  }

  cat("\n")
  cli::cli_h2("Key Insights & Distribution Profile")

  sp_sub <- x$splits_summary[x$splits_summary$Split != "TOTAL", ]
  ratio_str <- paste(paste0(sp_sub$Split, ": ", sp_sub$Pct_Images, "%"), collapse = " / ")

  imb_note <- if (x$imbalance_ratio == 1.0) {
    "Perfect balance"
  } else if (x$imbalance_ratio < 2.0) {
    sprintf("Well balanced (%.1fx max/min ratio)", x$imbalance_ratio)
  } else if (x$imbalance_ratio < 5.0) {
    sprintf("Moderate imbalance (%.1fx max/min ratio)", x$imbalance_ratio)
  } else {
    sprintf("High imbalance (%.1fx max/min ratio) - consider augmentation or class re-sampling", x$imbalance_ratio)
  }

  bullets <- c(
    "v" = "Split allocation: {ratio_str}",
    "v" = "Class balance: {imb_note}"
  )

  # Provenance insight
  if (!is.null(x$total_augmented) && x$total_augmented > 0L) {
    bullets <- c(bullets, "v" = "Augmentation: {.val {x$total_augmented}} images generated offline via data augmentation.")
  } else {
    bullets <- c(bullets, "i" = "Augmentation: 0 augmented images. Can expand training set via {.fun yolo_dataset_augment}.")
  }

  if (!is.null(x$total_synthetic) && x$total_synthetic > 0L) {
    bullets <- c(bullets, "v" = "Synthetic backgrounds: {.val {x$total_synthetic}} synthetic negative images included.")
  }

  tot_bg <- sum(sp_sub$Background)
  tot_imgs <- sum(sp_sub$Images)
  bg_pct <- if (tot_imgs > 0L) round(tot_bg / tot_imgs * 100, 1) else 0.0
  if (x$task != "classify") {
    if (tot_bg > 0L) {
      bullets <- c(bullets, "i" = "Background / Negative samples: {tot_bg} ({bg_pct}%) - healthy for false positive reduction.")
    } else {
      bullets <- c(bullets, "!" = "Background samples: 0 (0%) - consider adding 2-5% negative samples via {.fun yolo_dataset_add_background}.")
    }
  }

  if (length(x$sample_img_dims) > 0L) {
    dim_strings <- vapply(x$sample_img_dims, function(d) sprintf("%d x %d", d[1], d[2]), character(1))
    dim_tbl <- table(dim_strings)
    top_dim <- names(dim_tbl)[which.max(dim_tbl)]
    top_dim_pct <- round(max(dim_tbl) / length(dim_strings) * 100, 1)
    bullets <- c(bullets, "i" = "Dominant image resolution: {.val {top_dim}} ({top_dim_pct}% of inspected sample).")
  }

  cli::cli_bullets(bullets)
  invisible(x)
}


#' @export
plot.yolo_dataset_desc <- function(x,
                                  which = c("all", "classes", "splits", "provenance", "density"),
                                  col_palette = NULL,
                                  ...) {
  which <- match.arg(which)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  classes_df <- x$classes_summary
  splits_df <- x$splits_summary[x$splits_summary$Split != "TOTAL", ]

  splits <- setdiff(colnames(classes_df), c("ID", "Class", "Total_Objs", "Pct_Objs", "Images", "Pct_Images", "Avg_Area_Pct", "Size_Cat", "Total", "Pct_Total"))

  sp_cols <- c("#3b82f6", "#10b981", "#f59e0b", "#64748b", "#8b5cf6")[seq_along(splits_df$Split)]

  has_prov <- !is.null(x$provenance_summary) && nrow(x$provenance_summary[x$provenance_summary$Type != "TOTAL" & x$provenance_summary$Total > 0, , drop = FALSE]) > 1L
  has_dens <- !is.null(x$density_summary) && nrow(x$density_summary) > 1L

  if (which == "all") {
    if (has_prov || has_dens) {
      graphics::par(mfrow = c(2, 2), mar = c(4, 5, 3, 1) + 0.1)
    } else {
      graphics::par(mfrow = c(1, 2), mar = c(5, 7, 4, 2) + 0.1)
    }
  } else {
    graphics::par(mar = c(5, 7, 4, 2) + 0.1)
  }

  # Helper for class plot
  plot_classes <- function() {
    obj_counts <- classes_df$Total_Objs %||% classes_df$Total
    names(obj_counts) <- classes_df$Class
    pct_vals <- classes_df$Pct_Objs %||% classes_df$Pct_Total

    if (length(splits) > 1L) {
      mat <- as.matrix(classes_df[, splits, drop = FALSE])
      rownames(mat) <- classes_df$Class
      mat_t <- t(mat)

      bp <- graphics::barplot(
        mat_t,
        horiz = TRUE,
        las = 1,
        col = sp_cols[seq_along(splits)],
        border = NA,
        xlim = c(0, max(colSums(mat_t)) * 1.35),
        xlab = "Object Count",
        main = paste0("Class Distribution (Total: ", sum(mat_t), ")"),
        legend.text = splits,
        args.legend = list(x = "topright", bty = "n", cex = 0.8)
      )
      graphics::grid(nx = NULL, ny = NA, col = "#e2e8f0", lty = 1)
      graphics::barplot(
        mat_t,
        horiz = TRUE,
        las = 1,
        col = sp_cols[seq_along(splits)],
        border = NA,
        add = TRUE
      )
      tot_counts <- colSums(mat_t)
      graphics::text(
        x = tot_counts + max(tot_counts) * 0.03,
        y = bp,
        labels = sprintf("%d (%.1f%%)", tot_counts, pct_vals),
        adj = c(0, 0.5),
        cex = 0.8,
        font = 2
      )
    } else {
      cols <- col_palette %||% c("#3b82f6", "#10b981", "#f59e0b", "#ef4444", "#8b5cf6", "#ec4899", "#14b8a6")
      if (length(cols) < nrow(classes_df)) cols <- grDevices::rainbow(nrow(classes_df))

      bp <- graphics::barplot(
        obj_counts,
        horiz = TRUE,
        las = 1,
        col = cols[seq_along(obj_counts)],
        border = NA,
        xlim = c(0, max(obj_counts) * 1.35),
        xlab = "Object Count",
        main = paste0("Class Distribution (Total: ", sum(obj_counts), ")")
      )
      graphics::grid(nx = NULL, ny = NA, col = "#e2e8f0", lty = 1)
      graphics::barplot(
        obj_counts,
        horiz = TRUE,
        las = 1,
        col = cols[seq_along(obj_counts)],
        border = NA,
        add = TRUE
      )
      graphics::text(
        x = obj_counts + max(obj_counts) * 0.03,
        y = bp,
        labels = sprintf("%d (%.1f%%)", obj_counts, pct_vals),
        adj = c(0, 0.5),
        cex = 0.8,
        font = 2
      )
    }
  }

  # Helper for splits plot
  plot_splits <- function() {
    sp_counts <- splits_df$Images
    names(sp_counts) <- splits_df$Split
    sp_pcts <- splits_df$Pct_Images

    bp_sp <- graphics::barplot(
      sp_counts,
      col = sp_cols[seq_along(sp_counts)],
      border = NA,
      ylim = c(0, max(sp_counts) * 1.35),
      ylab = "Image Count",
      main = paste0("Dataset Splits (Total: ", sum(sp_counts), " images)"),
      las = 1
    )
    graphics::grid(nx = NA, ny = NULL, col = "#e2e8f0", lty = 1)
    graphics::barplot(
      sp_counts,
      col = sp_cols[seq_along(sp_counts)],
      border = NA,
      add = TRUE
    )
    graphics::text(
      x = bp_sp,
      y = sp_counts + max(sp_counts) * 0.05,
      labels = sprintf("%d\n(%.1f%%)", sp_counts, sp_pcts),
      cex = 0.8,
      font = 2
    )
  }

  # Helper for provenance plot
  plot_provenance <- function() {
    if (is.null(x$provenance_summary)) return(invisible())
    p_sub <- x$provenance_summary[x$provenance_summary$Type != "TOTAL" & x$provenance_summary$Total > 0, , drop = FALSE]
    if (nrow(p_sub) == 0L) return(invisible())
    p_counts <- p_sub$Total
    names(p_counts) <- p_sub$Type
    p_pcts <- p_sub$Pct
    p_cols <- c("Original" = "#3b82f6", "Augmented" = "#8b5cf6", "Synthetic_BG" = "#f59e0b", "Tiled" = "#10b981")[p_sub$Type]
    p_cols[is.na(p_cols)] <- "#64748b"

    bp_p <- graphics::barplot(
      p_counts,
      col = p_cols,
      border = NA,
      ylim = c(0, max(p_counts) * 1.35),
      ylab = "Image Count",
      main = paste0("Dataset Composition (", x$total_images, " images)"),
      las = 1
    )
    graphics::grid(nx = NA, ny = NULL, col = "#e2e8f0", lty = 1)
    graphics::barplot(p_counts, col = p_cols, border = NA, add = TRUE)
    graphics::text(
      x = bp_p,
      y = p_counts + max(p_counts) * 0.05,
      labels = sprintf("%d\n(%.1f%%)", p_counts, p_pcts),
      cex = 0.8,
      font = 2
    )
  }

  # Helper for density plot
  plot_density <- function() {
    if (is.null(x$density_summary) || nrow(x$density_summary) == 0L) return(invisible())
    d_counts <- x$density_summary$Images
    names(d_counts) <- x$density_summary$Bucket
    d_pcts <- x$density_summary$Pct_Images
    d_col <- "#0284c7"

    bp_d <- graphics::barplot(
      d_counts,
      col = d_col,
      border = NA,
      ylim = c(0, max(d_counts) * 1.35),
      ylab = "Image Count",
      main = "Object Density Distribution",
      las = 1
    )
    graphics::grid(nx = NA, ny = NULL, col = "#e2e8f0", lty = 1)
    graphics::barplot(d_counts, col = d_col, border = NA, add = TRUE)
    graphics::text(
      x = bp_d,
      y = d_counts + max(d_counts) * 0.05,
      labels = sprintf("%d\n(%.1f%%)", d_counts, d_pcts),
      cex = 0.75,
      font = 2
    )
  }

  if (which == "all") {
    plot_classes()
    plot_splits()
    if (has_prov || has_dens) {
      if (has_prov) plot_provenance() else plot_density()
      if (has_dens && has_prov) plot_density()
    }
  } else if (which == "classes") {
    plot_classes()
  } else if (which == "splits") {
    plot_splits()
  } else if (which == "provenance") {
    plot_provenance()
  } else if (which == "density") {
    plot_density()
  }

  invisible(x)
}


#' @title Standardize and Rename Image Files
#' @name image_rename
#' @aliases yolo_image_rename yolo_dataset_rename
#' @description
#' Standardizes and renames image files (and their corresponding YOLO label files, if present)
#' in a directory or YOLO dataset. Supports sequential numbering, zero-padding, custom prefix/suffix,
#' regex search-and-replace, and hyphenless UUID generation (RFC 4122 v4) to prevent collisions.
#'
#' @details
#' Renaming files is performed using a two-step collision-safe staging procedure, guaranteeing
#' that files are never accidentally overwritten even if old and new filenames overlap (e.g., renaming
#' `img_1.jpg` to `img_001.jpg` when another `img_001.jpg` already exists).
#'
#' When `rename_labels = TRUE`, matching `.txt` annotation files (either located in the same directory
#' or in a corresponding YOLO `labels/` directory, such as `labels/train` for `images/train`) are
#' automatically renamed in lockstep to keep annotations perfectly synchronized with their images.
#'
#' @param dir Character or `yolo_dataset` object. Directory containing the images to rename (default: `"."`).
#'   Can be a folder with images, a YOLO split directory (e.g. `images/train`), or a root YOLO dataset directory.
#' @param pattern Character or `NULL`. Regex pattern or file extension to filter images. If `NULL` (default),
#'   matches all standard image formats (`.jpg`, `.jpeg`, `.png`, `.bmp`, `.webp`, `.tif`, `.tiff`).
#'   Can also be a simple extension like `".jpg"` or `"jpg"`.
#' @param pattern_replace Character or `NULL`. Optional regex pattern to find in existing filenames.
#'   When specified alongside `replacement`, standardizes names by regex substitution rather than
#'   sequential numbering (e.g., `pattern_replace = "DSC_", replacement = "leaf_"`).
#' @param replacement Character or `NULL`. Replacement string used when `pattern_replace` is provided.
#' @param prefix Character. Prefix to add to the new filenames (default: `""`).
#' @param suffix Character. Optional suffix added before the file extension (default: `""`).
#' @param start_index Integer. Starting number for sequential numbering (default: `1L`).
#' @param pad Logical. If `TRUE` (default), pads numbers with leading zeros (e.g. `img_001.jpg` ... `img_381.jpg`)
#'   so that alphabetical sorting matches numerical sorting.
#' @param digits Integer or `NULL`. Number of digits for zero-padding. If `NULL` (default) and `pad = TRUE`,
#'   it is automatically determined from the total number of files (e.g., 3 digits for up to 999 files).
#' @param uuid Logical. If `TRUE` (default), generates a unique UUID for each file.
#' @param uuid_hyphen Logical. If `FALSE` (default), hyphens are removed from the UUID string
#'   (e.g., hexadecimal characters: `"e4d909c290"`).
#' @param uuid_position Character. Where to place the UUID relative to the prefix and index:
#'   `"only"` (default, e.g. `e4d909c290.jpg` or with prefix), `"suffix"` (e.g. `001_e4d909c290`),
#'   or `"prefix"` (e.g. `e4d909c290_001`).
#' @param uuid_length Integer. Length of the UUID string to keep (default: `10L`).
#'   Common values are 8, 10, 12, 16, or 32.
#' @param keep_prefix Logical or Character. If `TRUE` (default), automatically preserves characteristic prefixes
#'   and tags that identify synthetic backgrounds (`synth_bg_`, `bg_`, `background_`), negative samples (`neg_`, `empty_`),
#'   tiling slices (`_tile_0_0`), and data augmentations (`_aug1`, `_aug2_rot90`, `_mosaic`, `_cutout`).
#'   This ensures dataset provenance is preserved when standardizing names with UUIDs.
#'   Can also be `"all"` (to preserve any prefix before the first underscore), or a character vector / regex
#'   specifying custom prefixes to preserve (e.g. `keep_prefix = c("ds1_", "ds2_")` or `keep_prefix = "^(tora_|wood_)"`).
#'   Set to `FALSE` to disable preservation and rename files unconditionally.
#' @param rename_labels Logical. If `TRUE` (default), searches for matching `.txt` annotation files
#'   with the same stem (either in `dir` or in a corresponding YOLO `labels/` directory) and renames them.
#' @param recursive Logical. If `TRUE`, renames images recursively in subdirectories (e.g. `train`, `val`, `test`).
#'   Default is `FALSE`.
#' @param copy Logical. If `FALSE` (default), moves/renames files using safe staging. If `TRUE`, leaves original
#'   files intact and creates copies with the new names.
#' @param dry_run Logical. If `TRUE`, previews the planned renames without modifying any files on disk.
#'   Default is `FALSE`.
#' @param verbose Logical. If `TRUE` (default), prints progress and summary messages via `cli`.
#'
#' @return A `data.frame` containing:
#'   \item{original_name}{The original filename}
#'   \item{new_name}{The new filename}
#'   \item{original_path}{Full path to original image}
#'   \item{new_path}{Full path to new image}
#'   \item{status}{Outcome: `"renamed"`, `"copied"`, `"unchanged"`, or `"dry_run"`}
#'   \item{label_original}{Original label path, or `NA` if none found}
#'   \item{label_new}{New label path, or `NA`}
#'   Returned invisibly.
#'
#' @export
#' @examples
#' \dontrun{
#' # Standardize all jpgs in current directory to img_001.jpg, img_002.jpg, ...
#' image_rename(pattern = ".jpg", prefix = "img_")
#'
#' # Rename with a 32-character hyphenless UUID
#' image_rename(pattern = ".jpg", prefix = "sample_", uuid = TRUE)
#'
#' # Rename using only UUID without numbers: sample_e4d909c2...jpg
#' image_rename(pattern = ".jpg", prefix = "sample_", uuid = TRUE, uuid_position = "only")
#'
#' # Regex substitution: replace "IMG_" with "field_" in all image names
#' image_rename(pattern_replace = "IMG_", replacement = "field_")
#'
#' # Preview without renaming (dry-run)
#' image_rename(pattern = ".jpg", dry_run = TRUE)
#' }
image_rename <- function(dir = ".",
                         pattern = NULL,
                         pattern_replace = NULL,
                         replacement = NULL,
                         prefix = "",
                         suffix = "",
                         start_index = 1L,
                         digits = NULL,
                         pad = TRUE,
                         uuid = TRUE,
                         uuid_hyphen = FALSE,
                         uuid_position = c("only", "suffix", "prefix"),
                         uuid_length = 10L,
                         keep_prefix = TRUE,
                         rename_labels = TRUE,
                         recursive = FALSE,
                         copy = FALSE,
                         dry_run = FALSE,
                         verbose = TRUE) {
  uuid_position <- match.arg(uuid_position)

  # Extract directory path if a yolo_dataset or list is passed
  if (inherits(dir, "yolo_dataset") || (is.list(dir) && !is.null(dir$dir))) {
    dir_path <- dir$dir
  } else if (is.character(dir) && length(dir) == 1L) {
    dir_path <- dir
  } else {
    cli::cli_abort("{.arg dir} must be a single directory path or a {.code yolo_dataset} object.")
  }

  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    cli::cli_abort("Directory {.path {dir_path}} does not exist.")
  }

  # If user supplied replacement without pattern_replace, treat pattern as pattern_replace
  if (!is.null(replacement) && is.null(pattern_replace)) {
    pattern_replace <- pattern
    pattern <- NULL
  }

  # Resolve file matching pattern
  if (!is.null(pattern) && grepl("^\\.?[a-zA-Z0-9]+$", pattern)) {
    ext_clean <- sub("^\\.", "", pattern)
    file_pattern <- ext_clean
  } else if (!is.null(pattern)) {
    file_pattern <- pattern
  } else {
    file_pattern <- "(?i)\\.(jpg|jpeg|png|bmp|webp|tif|tiff)$"
  }

  # Auto-detect if images reside in a nested "images" folder (standard YOLO structure) or split folders (Roboflow structure)
  img_search_dir <- dir_path
  if (!recursive && dir.exists(file.path(dir_path, "images"))) {
    sub_imgs <- list.files(file.path(dir_path, "images"), pattern = file_pattern, full.names = TRUE, recursive = TRUE)
    if (length(sub_imgs) > 0L) {
      img_search_dir <- file.path(dir_path, "images")
      recursive <- TRUE
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Found YOLO {.path images/} subdirectory. Scanning recursively inside {.path {img_search_dir}}.")
      }
    }
  } else if (!recursive && any(dir.exists(file.path(dir_path, c("train", "val", "valid", "test"), "images")))) {
    sub_imgs <- list.files(dir_path, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
    if (length(sub_imgs) > 0L) {
      recursive <- TRUE
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Found Roboflow/split directory structure. Scanning recursively across splits inside {.path {img_search_dir}}.")
      }
    }
  }

  img_files <- list.files(img_search_dir, pattern = file_pattern, full.names = TRUE, recursive = recursive)
  # Filter out any hidden temp files from previous runs
  img_files <- img_files[!grepl("[/\\\\]\\.__tmp_ren_", img_files)]

  if (length(img_files) == 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_warning("No image files matching {.val {file_pattern}} found in {.path {img_search_dir}}.")
    }
    return(invisible(data.frame(
      original_name = character(0),
      new_name = character(0),
      original_path = character(0),
      new_path = character(0),
      status = character(0),
      label_original = character(0),
      label_new = character(0),
      stringsAsFactors = FALSE
    )))
  }

  n <- length(img_files)
  img_dirs <- dirname(img_files)
  img_basenames <- basename(img_files)
  img_exts <- tools::file_ext(img_basenames)
  img_stems <- tools::file_path_sans_ext(img_basenames)

  # Generate new stems
  if (!is.null(pattern_replace) && !is.null(replacement)) {
    new_stems <- gsub(pattern_replace, replacement, img_stems)
  } else {
    idx_seq <- seq(from = as.integer(start_index), length.out = n)
    if (isTRUE(pad)) {
      if (is.null(digits)) {
        digits <- max(nchar(as.character(max(idx_seq))), 3L)
      }
      fmt <- paste0("%0", digits, "d")
      idx_str <- sprintf(fmt, idx_seq)
    } else {
      idx_str <- as.character(idx_seq)
    }

    if (isTRUE(uuid)) {
      # Use built-in pliman uuid function with hyphen stripping
      raw_uids <- uuid(n = n, usetime = FALSE)
      if (!isTRUE(uuid_hyphen)) {
        raw_uids <- gsub("-", "", raw_uids, fixed = TRUE)
      }
      if (uuid_length < 32L) {
        raw_uids <- substr(raw_uids, 1L, uuid_length)
      }

      if (uuid_position == "suffix") {
        sep_mid <- if (nzchar(prefix) || nzchar(idx_str[1])) "_" else ""
        core_stems <- paste0(idx_str, sep_mid, raw_uids)
      } else if (uuid_position == "prefix") {
        sep_mid <- if (nzchar(prefix) || nzchar(idx_str[1])) "_" else ""
        core_stems <- paste0(raw_uids, sep_mid, idx_str)
      } else {
        # "only"
        core_stems <- raw_uids
      }
    } else {
      core_stems <- idx_str
    }

    # Handle keep_prefix (preserves synthetic background prefixes, tile/aug suffixes, or custom prefixes)
    if (!isFALSE(keep_prefix) && !is.null(keep_prefix)) {
      std_prefix_regex <- "^(synth_bg_|bg_synth_|synth_|bg_|background_|neg_|negative_|empty_)"
      tag_regex <- "(_tile_\\d+_\\d+.*|_tile_.*|_aug\\d+.*|_aug_.*|_mosaic.*|_cutout.*)$"

      if (isTRUE(keep_prefix)) {
        prefix_regex <- std_prefix_regex
      } else if (is.character(keep_prefix)) {
        if (identical(tolower(keep_prefix[1]), "all")) {
          prefix_regex <- paste0("(", std_prefix_regex, "|^([a-zA-Z0-9]+_))")
        } else {
          escaped <- vapply(keep_prefix, function(p) {
            p_esc <- gsub("([.\\^$*+?()[{\\\\|])", "\\\\\\1", p)
            paste0("^", p_esc)
          }, character(1))
          custom_patt <- paste0("(", paste(escaped, collapse = "|"), ")")
          prefix_regex <- paste0("(", std_prefix_regex, "|", custom_patt, ")")
        }
      } else {
        prefix_regex <- std_prefix_regex
      }

      pres_pfx <- rep("", n)
      if (nzchar(prefix_regex)) {
        has_pfx <- grepl(prefix_regex, img_stems, ignore.case = TRUE)
        if (any(has_pfx)) {
          m <- regexpr(prefix_regex, img_stems[has_pfx], ignore.case = TRUE)
          pres_pfx[has_pfx] <- regmatches(img_stems[has_pfx], m)
        }
      }

      pres_tag <- rep("", n)
      has_tag <- grepl(tag_regex, img_stems, ignore.case = TRUE)
      if (any(has_tag)) {
        m <- regexpr(tag_regex, img_stems[has_tag], ignore.case = TRUE)
        pres_tag[has_tag] <- regmatches(img_stems[has_tag], m)
      }

      n_pres_pfx <- sum(nzchar(pres_pfx))
      n_pres_tag <- sum(nzchar(pres_tag))
      if ((n_pres_pfx > 0L || n_pres_tag > 0L) && isTRUE(verbose)) {
        cli::cli_alert_info("Preserving provenance: {n_pres_pfx} background/prefix file{?s}, {n_pres_tag} tile/aug/slice file{?s}.")
      }

      new_stems <- character(n)
      for (i in seq_len(n)) {
        pf <- pres_pfx[i]
        tg <- pres_tag[i]
        if (nzchar(pf)) {
          new_stems[i] <- paste0(pf, prefix, core_stems[i], tg, suffix)
        } else {
          new_stems[i] <- paste0(prefix, core_stems[i], tg, suffix)
        }
      }
    } else {
      new_stems <- paste0(prefix, core_stems, suffix)
    }
  }

  new_basenames <- paste0(new_stems, ".", img_exts)
  new_img_files <- file.path(img_dirs, new_basenames)

  # Check matching YOLO annotation labels (.txt)
  lbl_orig <- rep(NA_character_, n)
  lbl_new <- rep(NA_character_, n)

  if (isTRUE(rename_labels)) {
    for (i in seq_len(n)) {
      stem_i <- img_stems[i]
      d_i <- img_dirs[i]

      # Scenario 1: Same directory as image
      cand_same <- file.path(d_i, paste0(stem_i, ".txt"))
      if (file.exists(cand_same)) {
        lbl_orig[i] <- cand_same
        lbl_new[i] <- file.path(d_i, paste0(new_stems[i], ".txt"))
      } else {
        # Scenario 2: Standard YOLO dataset structure (/images/ -> /labels/)
        if (grepl("(/|\\\\)images(/|\\\\|$)", d_i)) {
          d_lbl <- sub("((/|\\\\))images((/|\\\\)|$)", "\\1labels\\3", d_i)
          cand_yolo <- file.path(d_lbl, paste0(stem_i, ".txt"))
          if (file.exists(cand_yolo)) {
            lbl_orig[i] <- cand_yolo
            lbl_new[i] <- file.path(d_lbl, paste0(new_stems[i], ".txt"))
          }
        }
      }
    }
  }

  action_name <- if (copy) "copied" else "renamed"
  status_vec <- if (isTRUE(dry_run)) {
    rep("dry_run", n)
  } else {
    ifelse(img_files == new_img_files, "unchanged", action_name)
  }

  res_df <- data.frame(
    original_name = img_basenames,
    new_name = new_basenames,
    original_path = img_files,
    new_path = new_img_files,
    status = status_vec,
    label_original = lbl_orig,
    label_new = lbl_new,
    stringsAsFactors = FALSE
  )

  # Dry-run handling
  if (isTRUE(dry_run)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("[DRY-RUN] Preview of planned image renames (first 6 of {n}):")
      preview_n <- min(n, 6L)
      for (k in seq_len(preview_n)) {
        lbl_note <- if (!is.na(lbl_orig[k])) " [+ matching label]" else ""
        cli::cli_bullets(c(" " = paste0("{.val ", img_basenames[k], "} -> {.val ", new_basenames[k], "}", lbl_note)))
      }
      if (n > preview_n) {
        cli::cli_bullets(c(" " = "... and {n - preview_n} more files."))
      }
    }
    return(invisible(res_df))
  }

  # Execute renames / copies
  to_do <- which(img_files != new_img_files)

  if (length(to_do) > 0L) {
    if (isTRUE(copy)) {
      cp_ok <- file.copy(img_files[to_do], new_img_files[to_do], overwrite = TRUE)
      if (!all(cp_ok)) {
        cli::cli_warn("Some images could not be copied.")
      }
      has_lbl <- which(!is.na(lbl_orig) & lbl_orig != lbl_new)
      if (length(has_lbl) > 0L) {
        file.copy(lbl_orig[has_lbl], lbl_new[has_lbl], overwrite = TRUE)
      }
    } else {
      # Safe 2-step staging rename to prevent any collisions or data loss
      tmp_uuid <- tolower(gsub("-", "", uuid(1, usetime = FALSE)))
      tmp_imgs <- file.path(img_dirs[to_do], paste0(".__tmp_ren_", tmp_uuid, "_", seq_along(to_do), ".", img_exts[to_do]))

      # Stage 1: rename images to temp
      ok_s1 <- file.rename(img_files[to_do], tmp_imgs)
      if (!all(ok_s1)) {
        # Rollback partial stage 1
        rolled <- which(ok_s1)
        if (length(rolled) > 0L) {
          file.rename(tmp_imgs[rolled], img_files[to_do][rolled])
        }
        cli::cli_abort("Failed to safely stage files for renaming. Operation cancelled without data loss.")
      }

      # Stage 2: rename temp to target
      ok_s2 <- file.rename(tmp_imgs, new_img_files[to_do])
      if (!all(ok_s2)) {
        cli::cli_warn("Some files could not be moved to their final names. Please inspect directory.")
      }

      # Labels staging rename
      has_lbl <- which(!is.na(lbl_orig) & lbl_orig != lbl_new)
      if (length(has_lbl) > 0L) {
        lbl_dirs <- dirname(lbl_orig[has_lbl])
        tmp_lbls <- file.path(lbl_dirs, paste0(".__tmp_ren_lbl_", tmp_uuid, "_", seq_along(has_lbl), ".txt"))
        ok_l1 <- file.rename(lbl_orig[has_lbl], tmp_lbls)
        if (all(ok_l1)) {
          file.rename(tmp_lbls, lbl_new[has_lbl])
        } else {
          # Rollback labels if needed
          rolled_l <- which(ok_l1)
          if (length(rolled_l) > 0L) {
            file.rename(tmp_lbls[rolled_l], lbl_orig[has_lbl][rolled_l])
          }
          cli::cli_warn("Some label files could not be safely renamed.")
        }
      }
    }
  }

  if (isTRUE(verbose)) {
    n_lbls <- sum(!is.na(lbl_orig) & lbl_orig != lbl_new)
    lbl_text <- if (n_lbls > 0L) paste0(" and ", n_lbls, " corresponding label (.txt) file(s)") else ""
    if (isTRUE(copy)) {
      cli::cli_alert_success("Successfully copied {length(to_do)} image file(s){lbl_text}.")
    } else {
      cli::cli_alert_success("Successfully renamed {length(to_do)} image file(s){lbl_text}.")
    }
  }

  invisible(res_df)
}

#' @rdname image_rename
#' @export
yolo_image_rename <- image_rename

#' @rdname image_rename
#' @export
yolo_dataset_rename <- function(dir = ".",
                                pattern = NULL,
                                pattern_replace = NULL,
                                replacement = NULL,
                                prefix = "",
                                suffix = "",
                                start_index = 1L,
                                digits = NULL,
                                pad = TRUE,
                                uuid = TRUE,
                                uuid_hyphen = FALSE,
                                uuid_position = c("only", "suffix", "prefix"),
                                uuid_length = 10L,
                                keep_prefix = TRUE,
                                rename_labels = TRUE,
                                recursive = FALSE,
                                copy = FALSE,
                                dry_run = FALSE,
                                verbose = TRUE) {
  image_rename(
    dir = dir,
    pattern = pattern,
    pattern_replace = pattern_replace,
    replacement = replacement,
    prefix = prefix,
    suffix = suffix,
    start_index = start_index,
    digits = digits,
    pad = pad,
    uuid = uuid,
    uuid_hyphen = uuid_hyphen,
    uuid_position = match.arg(uuid_position),
    uuid_length = uuid_length,
    keep_prefix = keep_prefix,
    rename_labels = rename_labels,
    recursive = recursive,
    copy = copy,
    dry_run = dry_run,
    verbose = verbose
  )
}







