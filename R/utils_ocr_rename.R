#' Automated Image File Renaming and Sorting Using Deep Learning OCR
#'
#' @description
#' Scans a directory of images, detects and recognizes alphanumeric text tags
#' (e.g., plot tags, accession IDs, ear tags, seed packets), and automatically
#' renames or copies the files based on the recognized text.
#'
#' Includes advanced spatial positioning (`text_pos`) to select the relevant text
#' when multiple labels appear in the image, regular expression pattern matching,
#' conflict resolution, dry-run simulation, and full undo capabilities via
#' \code{\link{image_rename_undo}}.
#'
#' @param dir Character string. Path to the directory containing images to scan.
#'   Default is \code{"."}.
#' @param pattern Optional character string containing a regular expression to
#'   filter image files. If \code{NULL} (default), matches common image extensions:
#'   \code{"\\.(jpe?g|png|tif{1,2}|bmp|webp)$"}.
#' @param recursive Logical. Whether to scan subdirectories recursively. Default is \code{FALSE}.
#' @param text_pos Spatial or heuristic strategy used to select the target text
#'   when multiple text regions are detected in an image:
#'   \itemize{
#'     \item \code{"first"} (default): First detected text region in scan order.
#'     \item \code{"top-left"}: Closest bounding box to the top-left canvas corner.
#'     \item \code{"top-right"}: Closest bounding box to the top-right canvas corner.
#'     \item \code{"bottom-left"}: Closest bounding box to the bottom-left canvas corner.
#'     \item \code{"bottom-right"}: Closest bounding box to the bottom-right canvas corner.
#'     \item \code{"top"}: Text region positioned highest vertically (minimum y).
#'     \item \code{"bottom"}: Text region positioned lowest vertically (maximum y).
#'     \item \code{"center"}: Closest text bounding box centroid to the image center.
#'     \item \code{"largest"}: Text region with the largest bounding box area.
#'     \item \code{"confidence"}: Text region with the highest recognition confidence score.
#'   }
#' @param pattern_match Optional character string containing a regular expression.
#'   When provided, only recognized text strings matching this regex pattern are
#'   considered (e.g. \code{"^G[0-9]+"} or \code{"^[A-Z]{2,4}-[0-9]{3}"}).
#' @param prefix Character string to prepend to the new filename. Default is \code{""}.
#' @param suffix Character string to append to the new filename before the extension. Default is \code{""}.
#' @param case Character string specifying case transformation: \code{"as_is"} (default),
#'   \code{"upper"}, or \code{"lower"}.
#' @param sep Character string used as separator between prefix, recognized text, and suffix.
#'   Default is \code{"_"}.
#' @param min_chars Integer. Minimum number of characters required for a recognized text
#'   to be considered valid. Default is \code{1L}.
#' @param min_score Numeric. Minimum recognition confidence score (0 to 1) required. Default is \code{0.35}.
#' @param sanitize Logical. Whether to remove illegal filesystem characters and normalize spaces.
#'   Default is \code{TRUE}.
#' @param keep_extension Logical. Whether to preserve the original file extension. Default is \code{TRUE}.
#' @param copy Logical. If \code{TRUE}, files are copied rather than renamed/moved, preserving original files.
#'   Default is \code{FALSE}.
#' @param output_dir Optional character string specifying destination directory. If \code{NULL} (default),
#'   files are renamed/copied within their source directory.
#' @param conflict Character string specifying how to resolve filename collisions:
#'   \itemize{
#'     \item \code{"rename"} (default): Append numeric suffix (\code{"_1"}, \code{"_2"}, etc.).
#'     \item \code{"skip"}: Do not rename or overwrite; skip the file.
#'     \item \code{"overwrite"}: Overwrite existing target file.
#'   }
#' @param dry_run Logical. If \code{TRUE}, simulates the entire OCR scanning and renaming process
#'   without modifying any files on disk, returning a preview table. Default is \code{FALSE}.
#' @param save_log Logical. Whether to save an audit log file (\code{.pliman_rename_log_*.rds})
#'   for subsequent recovery with \code{\link{image_rename_undo}}. Default is \code{TRUE}.
#' @param det_model Text detection model backbone: \code{"ppocr_det"} (default),
#'   \code{"ppocr_det_server"}, or \code{"craft"}.
#' @param rec_model Text recognition model backbone. Default is \code{"ppocr_rec"}.
#' @param det_thresh Pixel binarization threshold for DBNet detection. Default is \code{0.30}.
#' @param box_thresh Bounding box confidence threshold. Default is \code{0.50}.
#' @param unclip_ratio Polygon expansion ratio. Default is \code{1.30}.
#' @param text_threshold CRAFT character threshold. Default is \code{0.50}.
#' @param link_threshold CRAFT link threshold. Default is \code{0.40}.
#' @param low_text CRAFT low text threshold. Default is \code{0.40}.
#' @param engine Compute engine: \code{"gpu"} (default) or \code{"cpu"}.
#' @param threads Number of CPU threads (0 = auto-detect).
#' @param dir_models Directory where ONNX models are stored. Default is \code{pliman_model_dir()}.
#' @param verbose Logical. Whether to show informative progress messages. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{image_ocr_dl}}.
#'
#' @return An object of class \code{c("image_rename_ocr", "data.frame")} containing the audit log
#'   of all processed files with columns:
#'   \itemize{
#'     \item \code{original_path}: Full path to original file.
#'     \item \code{original_name}: Original file basename.
#'     \item \code{detected_text}: Text selected for renaming.
#'     \item \code{score}: Recognition confidence score.
#'     \item \code{text_pos}: Selection rule applied.
#'     \item \code{new_name}: Final filename generated.
#'     \item \code{new_path}: Full destination path.
#'     \item \code{status}: Execution outcome (\code{"renamed"}, \code{"copied"}, \code{"preview"},
#'       \code{"no_text_detected"}, \code{"no_pattern_match"}, \code{"skipped_conflict"}, \code{"error"}).
#'     \item \code{action}: Performed action (\code{"rename"} or \code{"copy"}).
#'   }
#'
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examplesIf interactive()
#'   library(pliman)
#'
#'   # Dry-run preview in a folder of plot tag images
#'   preview <- image_rename_ocr(
#'     dir = "field_photos",
#'     text_pos = "top-left",
#'     pattern_match = "^G[0-9]+",
#'     dry_run = TRUE
#'   )
#'   print(preview)
#'
#'   # Apply rename safely with audit log
#'   # res <- image_rename_ocr(dir = "field_photos", text_pos = "center")
#'   # Undo if needed:
#'   # image_rename_undo(res)
image_rename_ocr <- function(dir = ".",
                             pattern = NULL,
                             recursive = FALSE,
                             text_pos = c("first", "top-left", "top-right", "bottom-left",
                                          "bottom-right", "top", "bottom", "center",
                                          "largest", "confidence"),
                             pattern_match = NULL,
                             prefix = "",
                             suffix = "",
                             case = c("as_is", "upper", "lower"),
                             sep = "_",
                             min_chars = 1L,
                             min_score = 0.35,
                             sanitize = TRUE,
                             keep_extension = TRUE,
                             copy = FALSE,
                             output_dir = NULL,
                             conflict = c("rename", "skip", "overwrite"),
                             dry_run = FALSE,
                             save_log = TRUE,
                             det_model = c("ppocr_det", "ppocr_det_server", "craft"),
                             rec_model = "ppocr_rec",
                             det_thresh = 0.30,
                             box_thresh = 0.50,
                             unclip_ratio = 1.30,
                             text_threshold = 0.50,
                             link_threshold = 0.40,
                             low_text = 0.40,
                             engine = c("gpu", "cpu"),
                             threads = 0,
                             dir_models = pliman_model_dir(),
                             verbose = TRUE,
                             ...) {
  text_pos <- match.arg(text_pos)
  case <- match.arg(case)
  conflict <- match.arg(conflict)
  engine <- match.arg(engine)
  if (is.character(det_model) && length(det_model) > 1L) {
    det_model <- det_model[1L]
  }

  if (!dir.exists(dir)) {
    stop("Directory does not exist: ", dir)
  }

  if (is.null(pattern)) {
    pattern <- "(?i)\\.(jpe?g|png|tif{1,2}|bmp|webp)$"
  }

  files <- list.files(path = dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(files) == 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_warning("No matching image files found in '{dir}'.")
    }
    empty_df <- data.frame(
      original_path = character(0), original_name = character(0),
      detected_text = character(0), score = numeric(0),
      text_pos = character(0), new_name = character(0),
      new_path = character(0), status = character(0),
      action = character(0), stringsAsFactors = FALSE
    )
    class(empty_df) <- c("image_rename_ocr", "data.frame")
    return(empty_df)
  }

  # Ensure output directory exists if provided
  if (!is.null(output_dir)) {
    if (!dry_run && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  n_files <- length(files)
  if (isTRUE(verbose)) {
    act_word <- if (dry_run) "Previewing OCR renaming for" else (if (copy) "OCR scanning & copying" else "OCR scanning & renaming")
    cli::cli_h1("{act_word} {n_files} image file{?s}")
    cli::cli_alert_info("Selection strategy: {.field {text_pos}} | Engine: {.val {toupper(engine)}}")
    if (!is.null(pattern_match)) {
      cli::cli_alert_info("Regex text filter: {.code {pattern_match}}")
    }
  }

  # Prepare helper to sanitize filename
  sanitize_name <- function(txt) {
    if (!isTRUE(sanitize)) return(txt)
    # Remove control characters and illegal filename characters
    s <- gsub('[\r\n\t\\/\\\\:\'*?"<>|]', "", txt)
    # Replace whitespace with sep
    s <- gsub("\\s+", sep, s)
    # Consolidate multiple consecutive separators
    if (nzchar(sep)) {
      pattern_dup <- paste0(regex_escape(sep), "+")
      s <- gsub(pattern_dup, sep, s)
      # Strip leading/trailing separators
      s <- gsub(paste0("^", regex_escape(sep), "|", regex_escape(sep), "$"), "", s)
    }
    s
  }

  regex_escape <- function(str) {
    gsub("([.\\^$*+?()[{\\\\|])", "\\\\\\1", str)
  }

  results <- vector("list", n_files)

  for (i in seq_len(n_files)) {
    f_path <- files[i]
    f_orig_name <- basename(f_path)
    f_dir <- if (!is.null(output_dir)) output_dir else dirname(f_path)
    f_ext <- tools::file_ext(f_orig_name)
    dot_ext <- if (keep_extension && nzchar(f_ext)) paste0(".", f_ext) else ""

    # Run OCR on single image
    ocr_res <- tryCatch({
      image_ocr_dl(
        img = f_path,
        det_model = det_model,
        rec_model = rec_model,
        det_thresh = det_thresh,
        box_thresh = box_thresh,
        unclip_ratio = unclip_ratio,
        text_threshold = text_threshold,
        link_threshold = link_threshold,
        low_text = low_text,
        engine = engine,
        threads = threads,
        dir = dir_models,
        plot = FALSE,
        verbose = FALSE,
        ...
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(ocr_res) || nrow(ocr_res) == 0L) {
      results[[i]] <- list(
        original_path = f_path,
        original_name = f_orig_name,
        detected_text = NA_character_,
        score = NA_real_,
        text_pos = text_pos,
        new_name = f_orig_name,
        new_path = f_path,
        status = if (is.null(ocr_res)) "error" else "no_text_detected",
        action = if (copy) "copy" else "rename"
      )
      next
    }

    # Filter out empty, low confidence, or short texts
    candidates <- ocr_res[nzchar(ocr_res$text) &
                          nchar(ocr_res$text) >= min_chars &
                          ocr_res$score >= min_score, , drop = FALSE]

    # Filter by regex pattern_match if specified
    if (!is.null(pattern_match) && nrow(candidates) > 0L) {
      match_idx <- grepl(pattern_match, candidates$text, perl = TRUE)
      candidates <- candidates[match_idx, , drop = FALSE]
      if (nrow(candidates) == 0L) {
        results[[i]] <- list(
          original_path = f_path,
          original_name = f_orig_name,
          detected_text = NA_character_,
          score = NA_real_,
          text_pos = text_pos,
          new_name = f_orig_name,
          new_path = f_path,
          status = "no_pattern_match",
          action = if (copy) "copy" else "rename"
        )
        next
      }
    }

    if (nrow(candidates) == 0L) {
      results[[i]] <- list(
        original_path = f_path,
        original_name = f_orig_name,
        detected_text = NA_character_,
        score = NA_real_,
        text_pos = text_pos,
        new_name = f_orig_name,
        new_path = f_path,
        status = "no_text_detected",
        action = if (copy) "copy" else "rename"
      )
      next
    }

    # Retrieve image dimensions for spatial metrics
    img_attr <- attr(ocr_res, "image")
    if (!is.null(img_attr)) {
      dims <- dim(image_data(img_attr))
      canvas_w <- dims[1]
      canvas_h <- dims[2]
    } else {
      canvas_w <- max(candidates$xmax, na.rm = TRUE)
      canvas_h <- max(candidates$ymax, na.rm = TRUE)
    }

    # Select the winning text region based on text_pos
    chosen_idx <- 1L
    if (nrow(candidates) > 1L) {
      cx <- (candidates$xmin + candidates$xmax) / 2
      cy <- (candidates$ymin + candidates$ymax) / 2
      area <- (candidates$xmax - candidates$xmin + 1) * (candidates$ymax - candidates$ymin + 1)

      chosen_idx <- switch(
        text_pos,
        "first" = 1L,
        "top-left" = which.min(candidates$xmin^2 + candidates$ymin^2),
        "top-right" = which.min((canvas_w - candidates$xmax)^2 + candidates$ymin^2),
        "bottom-left" = which.min(candidates$xmin^2 + (canvas_h - candidates$ymax)^2),
        "bottom-right" = which.min((canvas_w - candidates$xmax)^2 + (canvas_h - candidates$ymax)^2),
        "top" = which.min(candidates$ymin),
        "bottom" = which.max(candidates$ymax),
        "center" = which.min((cx - canvas_w / 2)^2 + (cy - canvas_h / 2)^2),
        "largest" = which.max(area),
        "confidence" = which.max(candidates$score),
        1L
      )
    }

    chosen_row <- candidates[chosen_idx, , drop = FALSE]
    raw_text <- chosen_row$text
    chosen_score <- chosen_row$score

    # Apply case transformation
    formatted_text <- switch(
      case,
      "upper" = toupper(raw_text),
      "lower" = tolower(raw_text),
      raw_text
    )

    clean_text <- sanitize_name(formatted_text)

    # Compose new filename parts
    name_parts <- character(0)
    if (nzchar(prefix)) name_parts <- c(name_parts, prefix)
    name_parts <- c(name_parts, clean_text)
    if (nzchar(suffix)) name_parts <- c(name_parts, suffix)

    base_new_name <- paste(name_parts, collapse = sep)
    candidate_new_name <- paste0(base_new_name, dot_ext)
    target_path <- file.path(f_dir, candidate_new_name)

    # Handle collisions/conflicts
    final_status <- if (dry_run) "preview" else (if (copy) "copied" else "renamed")
    if (file.exists(target_path) && normalizePath(target_path, mustWork = FALSE) != normalizePath(f_path, mustWork = FALSE)) {
      if (conflict == "skip") {
        results[[i]] <- list(
          original_path = f_path,
          original_name = f_orig_name,
          detected_text = raw_text,
          score = chosen_score,
          text_pos = text_pos,
          new_name = candidate_new_name,
          new_path = target_path,
          status = "skipped_conflict",
          action = if (copy) "copy" else "rename"
        )
        next
      } else if (conflict == "rename") {
        counter <- 1L
        while (file.exists(target_path)) {
          candidate_new_name <- paste0(base_new_name, sep, counter, dot_ext)
          target_path <- file.path(f_dir, candidate_new_name)
          counter <- counter + 1L
        }
      }
    }

    # Perform action on disk if not dry_run
    if (!dry_run) {
      if (copy) {
        ok <- file.copy(from = f_path, to = target_path, overwrite = (conflict == "overwrite"))
        if (!ok) final_status <- "error"
      } else {
        ok <- file.rename(from = f_path, to = target_path)
        if (!ok) final_status <- "error"
      }
    }

    results[[i]] <- list(
      original_path = f_path,
      original_name = f_orig_name,
      detected_text = raw_text,
      score = chosen_score,
      text_pos = text_pos,
      new_name = candidate_new_name,
      new_path = target_path,
      status = final_status,
      action = if (copy) "copy" else "rename"
    )
  }

  res_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
  class(res_df) <- c("image_rename_ocr", "data.frame")

  # Save audit log for recovery if requested
  if (isTRUE(save_log) && !dry_run) {
    timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_name <- paste0(".pliman_rename_log_", timestamp_str, ".rds")
    log_path <- file.path(if (!is.null(output_dir)) output_dir else dir, log_name)
    saveRDS(res_df, file = log_path)
    attr(res_df, "log_file") <- log_path
  }

  if (isTRUE(verbose)) {
    print.image_rename_ocr(res_df)
  }

  invisible(res_df)
}

#' Undo Previous OCR Renaming or Copying Operations
#'
#' @description
#' Restores original filenames by inverting an audit log produced by \code{\link{image_rename_ocr}}.
#' If files were renamed, they are moved back to their original names. If files were copied,
#' the copied instances are safely removed.
#'
#' @param log An \code{image_rename_ocr} object returned by \code{\link{image_rename_ocr}},
#'   or a character string specifying the path to a \code{.pliman_rename_log_*.rds} file.
#'   If \code{NULL} (default), searches for the most recent log file in \code{dir}.
#' @param dir Character string. Path to search for audit log files when \code{log = NULL}.
#'   Default is \code{"."}.
#' @param verbose Logical. Whether to show progress messages. Default is \code{TRUE}.
#'
#' @return A \code{data.frame} summarizing the inverted operations and their status.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examplesIf interactive()
#'   library(pliman)
#'   # Undo the most recent rename in the current folder:
#'   # image_rename_undo()
image_rename_undo <- function(log = NULL, dir = ".", verbose = TRUE) {
  if (is.null(log)) {
    log_files <- list.files(path = dir, pattern = "^\\.pliman_rename_log_.*\\.rds$", full.names = TRUE, all.files = TRUE)
    if (length(log_files) == 0L) {
      stop("No pliman rename log files found in '", dir, "'. Specify 'log' explicitly.")
    }
    # Pick the most recent log file
    file_info <- file.info(log_files)
    latest_file <- log_files[which.max(file_info$mtime)]
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Using most recent audit log: {.file {basename(latest_file)}}")
    }
    log_df <- readRDS(latest_file)
  } else if (is.character(log) && file.exists(log)) {
    log_df <- readRDS(log)
  } else if (inherits(log, "image_rename_ocr") || is.data.frame(log)) {
    log_df <- log
  } else {
    stop("Invalid 'log' argument. Must be a log file path or an 'image_rename_ocr' data frame.")
  }

  # Invert successful rename or copy actions
  invertible <- log_df[log_df$status %in% c("renamed", "copied"), , drop = FALSE]
  if (nrow(invertible) == 0L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("No reversible renaming or copying operations found in log.")
    }
    return(invisible(log_df))
  }

  undo_statuses <- character(nrow(invertible))

  for (i in rev(seq_len(nrow(invertible)))) {
    curr <- invertible[i, ]
    act <- curr$action
    orig <- curr$original_path
    cur_p <- curr$new_path

    if (act == "rename") {
      if (file.exists(cur_p)) {
        ok <- file.rename(from = cur_p, to = orig)
        undo_statuses[i] <- if (ok) "reverted" else "failed_revert"
      } else {
        undo_statuses[i] <- "target_missing"
      }
    } else if (act == "copy") {
      if (file.exists(cur_p)) {
        ok <- file.remove(cur_p)
        undo_statuses[i] <- if (ok) "copy_removed" else "failed_removal"
      } else {
        undo_statuses[i] <- "copy_missing"
      }
    }
  }

  invertible$undo_status <- undo_statuses

  if (isTRUE(verbose)) {
    n_rev <- sum(undo_statuses %in% c("reverted", "copy_removed"))
    cli::cli_alert_success("Successfully reverted {n_rev} of {nrow(invertible)} file operation{?s}.")
  }

  invisible(invertible)
}

#' @export
print.image_rename_ocr <- function(x, ...) {
  if (nrow(x) == 0L) {
    cli::cli_alert_info("No files processed.")
    return(invisible(x))
  }

  n_total <- nrow(x)
  n_renamed <- sum(x$status == "renamed")
  n_copied <- sum(x$status == "copied")
  n_preview <- sum(x$status == "preview")
  n_no_text <- sum(x$status == "no_text_detected")
  n_no_pat <- sum(x$status == "no_pattern_match")
  n_skip_conf <- sum(x$status == "skipped_conflict")
  n_err <- sum(x$status == "error")

  if (n_preview > 0L) {
    cli::cli_h2("OCR Rename Simulation (Dry Run) - {n_total} file{?s}")
  } else {
    cli::cli_h2("OCR Rename Audit - {n_total} file{?s}")
  }

  if (n_renamed > 0L) cli::cli_alert_success("Renamed: {.val {n_renamed}} file{?s}")
  if (n_copied > 0L) cli::cli_alert_success("Copied: {.val {n_copied}} file{?s}")
  if (n_preview > 0L) cli::cli_alert_info("Previewed: {.val {n_preview}} file{?s}")
  if (n_no_text > 0L) cli::cli_alert_warning("No text detected: {.val {n_no_text}} file{?s}")
  if (n_no_pat > 0L) cli::cli_alert_warning("No pattern match: {.val {n_no_pat}} file{?s}")
  if (n_skip_conf > 0L) cli::cli_alert_warning("Skipped on conflict: {.val {n_skip_conf}} file{?s}")
  if (n_err > 0L) cli::cli_alert_danger("Errors: {.val {n_err}} file{?s}")

  show_cols <- intersect(c("original_name", "detected_text", "score", "new_name", "status"), colnames(x))
  print.data.frame(x[, show_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}
