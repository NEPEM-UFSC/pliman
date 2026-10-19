# ==============================================================================
# BATCH DEEP LEARNING INFERENCE
# pliman: Plant Image Analysis in R
# ==============================================================================

#' Batch Object Detection with Deep Learning
#'
#' Runs deep learning object detection (e.g., YOLO models) on a batch of images
#' defined by a directory path, a filename pattern, a vector of file paths,
#' or a list of `Image` objects. Supports sequential or parallel processing (via `mirai`),
#' and returns a consolidated predictions dataset along with per-image and per-class summaries.
#'
#' @param data Input data. Can be:
#'   * A character string representing a directory path containing images (e.g., `"flax_grains"`).
#'   * A character pattern to match image files in `path` or the current directory (e.g., `"-"`, `"^img_"`).
#'   * A character vector of image file paths.
#'   * A list of `Image` objects or numeric arrays.
#' @param pattern Optional regex pattern to filter image filenames when `data` is a directory path.
#' @param path Optional directory path to search when `data` is a pattern. Defaults to `"."`.
#' @param parallel Logical. If `TRUE`, runs inference in parallel across images using `mirai`.
#'   Defaults to `FALSE`.
#' @param workers Integer. Number of parallel workers/daemons to launch when `parallel = TRUE`.
#'   Defaults to 60% of available CPU cores.
#' @param show_point Logical. Whether to show a point at the centroid of each detected object.
#'   Defaults to `FALSE`.
#' @param keep_images Logical. Whether to keep the full original images attached as attributes in `$results`.
#'   Defaults to `FALSE` to prevent excessive RAM consumption during large batches (saving gigabytes of memory).
#' @param save_image Logical. Whether to save annotated prediction images with bounding boxes to disk.
#'   Defaults to `FALSE`.
#' @param dir_processed Character string. Directory where annotated images should be saved when
#'   `save_image = TRUE`. Defaults to `"./proc"`.
#' @param prefix Character string. Prefix prepended to saved output filenames. Defaults to `"proc_"`.
#' @param verbose Logical. Whether to show progress bars, timing, and formatted CLI summaries.
#'   Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_detect_dl()] (e.g., `model`, `conf_threshold`,
#'   `iou_threshold`, `classes`, `labels`, `col`, `rainbow`, `lwd`, `cex`, `badge`, `engine`, `device_id`).
#'
#' @return An S3 object of class `c("dl_detect_batch", "list")` containing:
#'   * `predictions`: A consolidated `data.frame` of all detected bounding boxes across all images,
#'     with column `img` identifying the source image.
#'   * `summary`: A per-image `data.frame` summarizing total objects detected, mean confidence scores,
#'     and individual counts for each detected class.
#'   * `class_summary`: A summary `data.frame` across the entire batch with total counts, frequencies,
#'     average per image, and mean confidence scores per class.
#'   * `results`: A named list of raw individual detection outputs from [image_detect_dl()].
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Batch detection in a specific folder (Sequential)
#' res <- image_detect_dl_batch(
#'   data = "flax_grains",
#'   model = "yolo26n",
#'   conf_threshold = 0.35
#' )
#' head(res$summary)
#' head(res$predictions)
#'
#' # 2. Batch detection using filename pattern in parallel
#' res_par <- image_detect_dl_batch(
#'   data = "-",
#'   path = "flax_grains",
#'   parallel = TRUE,
#'   workers = 4
#' )
#' res_par$class_summary
image_detect_dl_batch <- function(data,
                                  pattern = NULL,
                                  path = NULL,
                                  parallel = FALSE,
                                  workers = NULL,
                                  show_point = FALSE,
                                  keep_images = FALSE,
                                  save_image = FALSE,
                                  dir_processed = "./proc",
                                  prefix = "proc_",
                                  verbose = TRUE,
                                  ...) {
  # 1. Resolve image inputs
  resolved <- .resolve_batch_input(data = data, pattern = pattern, path = path)
  items <- resolved$items
  names_items <- resolved$names
  n_items <- length(items)

  if (n_items == 0L) {
    cli::cli_abort("No valid images found for batch processing.")
  }

  dots_args <- list(...)
  dots_args$plot <- FALSE
  dots_args$verbose <- FALSE
  dots_args$show_point <- show_point
  dots_args$attach_image <- isTRUE(keep_images)

  t0 <- Sys.time()

  if (isTRUE(verbose)) {
    cli::cli_rule(
      left = cli::col_blue("pliman: Deep Learning Batch Detection"),
      right = cli::col_blue(format(t0, "%Y-%m-%d | %H:%M:%S"))
    )
    cli::cli_alert_info("Found {.val {n_items}} image{?s} to process.")
  }

  worker_detect_fn <- function(im, args) {
    do.call(pliman::image_detect_dl, c(list(img = im), args))
  }

  # 2. Execution: Parallel vs Sequential
  if (isTRUE(parallel)) {
    if (!requireNamespace("mirai", quietly = TRUE)) {
      cli::cli_abort("Package {.val mirai} is required for parallel processing. Please install it with {.code install.packages('mirai')}.")
    }

    nworkers <- if (is.null(workers)) {
      max(1L, floor(parallel::detectCores() * 0.6))
    } else {
      as.integer(workers)
    }
    eff_workers <- min(nworkers, n_items)

    st <- tryCatch(mirai::status()$daemons, error = function(e) 0L)
    existing_daemons <- if (is.numeric(st)) st else 0L

    daemon_status <- tryCatch({
      if (existing_daemons == 0L || (!is.null(workers) && existing_daemons != eff_workers)) {
        empty_prof <- tempfile()
        file.create(empty_prof)
        old_prof <- Sys.getenv("R_PROFILE_USER", unset = NA)
        Sys.setenv(R_PROFILE_USER = empty_prof)

        mirai::daemons(eff_workers)

        if (is.na(old_prof)) Sys.unsetenv("R_PROFILE_USER") else Sys.setenv(R_PROFILE_USER = old_prof)
        unlink(empty_prof)

        if (existing_daemons == 0L) {
          on.exit(mirai::daemons(0), add = TRUE)
        }
        mirai::everywhere(
          {
            .libPaths(lp)
            library(pliman)
          },
          lp = .libPaths()
        )
      }
      TRUE
    }, error = function(e) {
      cli::cli_warn("Failed to initialize parallel daemons: {e$message}. Falling back to sequential execution.")
      FALSE
    })

    if (isTRUE(daemon_status)) {
      if (isTRUE(verbose)) {
        cli::cli_rule(
          left = cli::col_blue("Parallel detection using {.val {eff_workers}} cores"),
          right = cli::col_blue(format(Sys.time(), "%H:%M:%S"))
        )
      }

      n_chunks <- if (n_items <= eff_workers) eff_workers else min(n_items, eff_workers * 3L)
      chunk_size <- ceiling(n_items / n_chunks)
      chunk_factor <- rep(seq_len(n_chunks), each = chunk_size, length.out = n_items)
      chunks <- split(items, chunk_factor)

      results_chunks <- mirai::mirai_map(
        .x = chunks,
        .f = function(chunk_items, fn, args) {
          lapply(chunk_items, function(it) {
            fn(it, args)
          })
        },
        .args = list(fn = worker_detect_fn, args = dots_args)
      )[.progress]

      # Check for worker-level errors
      failed_chunks <- sapply(results_chunks, function(res) inherits(res, "errorValue") || !is.list(res))
      if (any(failed_chunks)) {
        err_msg <- paste(sapply(which(failed_chunks), function(idx) {
          if (inherits(results_chunks[[idx]], "errorValue")) as.character(results_chunks[[idx]]) else "Worker failure"
        }), collapse = "; ")
        cli::cli_abort("Parallel detection failed: {err_msg}")
      }

      results <- do.call(c, results_chunks)
      if (!isTRUE(keep_images)) {
        for (i in seq_len(n_items)) {
          if (!is.null(results[[i]])) attr(results[[i]], "image") <- NULL
        }
        gc(verbose = FALSE)
      }
    } else {
      parallel <- FALSE
    }
  }

  if (!isTRUE(parallel)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        total = n_items,
        clear = FALSE
      )
    }
    results <- vector("list", n_items)
    for (i in seq_len(n_items)) {
      if (isTRUE(verbose)) cli::cli_progress_update()
      res_i <- worker_detect_fn(items[[i]], dots_args)
      if (!isTRUE(keep_images) && !is.null(res_i)) {
        attr(res_i, "image") <- NULL
      }
      results[[i]] <- res_i
      if (i %% 50L == 0L) gc(verbose = FALSE)
    }
    if (isTRUE(verbose)) cli::cli_progress_done()
  }

  names(results) <- names_items

  # 3. Optional saving of annotated images
  if (isTRUE(save_image)) {
    if (!dir.exists(dir_processed)) {
      dir.create(dir_processed, recursive = TRUE)
    }
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Saving annotated detection images to {.path {dir_processed}}...",
        msg_done = "Annotated images successfully saved"
      )
    }
    for (i in seq_len(n_items)) {
      res_i <- results[[i]]
      im_obj <- attr(res_i, "image")
      if (is.null(im_obj) && is.character(items[[i]])) {
        im_obj <- tryCatch(image_import(items[[i]]), error = function(e) NULL)
      }
      if (!is.null(im_obj)) {
        out_file <- file.path(dir_processed, paste0(prefix, names_items[i], if (!grepl("\\.[a-zA-Z]+$", names_items[i])) ".jpg" else ""))
        tryCatch({
          dims <- dim(image_data(im_obj))
          grDevices::png(out_file, width = dims[1], height = dims[2])
          on.exit(grDevices::dev.off(), add = TRUE)
          plot(im_obj)
          if (is.data.frame(res_i) && nrow(res_i) > 0) {
            pal <- .get_class_palette(
              labels = res_i$label,
              n = nrow(res_i),
              rainbow = if (!is.null(dots_args$rainbow)) dots_args$rainbow else TRUE,
              user_col = dots_args$col
            )
            .plot_yolo_bboxes(
              boxes = res_i,
              palette_colors = pal,
              lwd = if (!is.null(dots_args$lwd)) dots_args$lwd else 2,
              cex = if (!is.null(dots_args$cex)) dots_args$cex else 1.0,
              pad = if (!is.null(dots_args$pad)) dots_args$pad else 1.0,
              show_text = if (!is.null(dots_args$show_text)) dots_args$show_text else TRUE,
              show_conf = if (!is.null(dots_args$show_conf)) dots_args$show_conf else TRUE,
              show_class = if (!is.null(dots_args$show_class)) dots_args$show_class else TRUE,
              show_id = if (!is.null(dots_args$show_id)) dots_args$show_id else FALSE,
              show_box = dots_args$show_box,
              show_point = show_point,
              point_size = if (!is.null(dots_args$point_size)) dots_args$point_size else 1.5,
              point_pch = if (!is.null(dots_args$point_pch)) dots_args$point_pch else 16,
              point_col = dots_args$point_col,
              badge = if (!is.null(dots_args$badge)) dots_args$badge else TRUE
            )
          }
          grDevices::dev.off()
        }, error = function(e) NULL)
        rm(im_obj)
        if (i %% 50L == 0L) gc(verbose = FALSE)
      }
    }
    if (!isTRUE(keep_images)) {
      for (i in seq_len(n_items)) {
        if (!is.null(results[[i]])) attr(results[[i]], "image") <- NULL
      }
      gc(verbose = FALSE)
    }
  }

  # 4. Consolidate predictions
  pred_list <- lapply(seq_len(n_items), function(i) {
    df <- results[[i]]
    if (is.data.frame(df) && nrow(df) > 0) {
      data.frame(img = names_items[i], df, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  predictions <- do.call(rbind, pred_list)
  if (is.null(predictions)) {
    predictions <- data.frame(
      img = character(0), id = integer(0), xmin = numeric(0), ymin = numeric(0),
      xmax = numeric(0), ymax = numeric(0), label = character(0), score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # 5. Build per-image summary
  all_classes <- if (nrow(predictions) > 0) sort(unique(predictions$label)) else character(0)
  summary_list <- lapply(seq_len(n_items), function(i) {
    im_name <- names_items[i]
    df <- results[[i]]
    n_obj <- if (is.data.frame(df)) nrow(df) else 0L
    mean_sc <- if (n_obj > 0 && "score" %in% names(df)) round(mean(df$score, na.rm = TRUE), 4) else NA_real_

    row_data <- list(img = im_name, num_objects = n_obj, mean_score = mean_sc)
    if (length(all_classes) > 0) {
      if (n_obj > 0) {
        cls_tab <- table(factor(df$label, levels = all_classes))
        for (cl in all_classes) {
          row_data[[paste0("count_", cl)]] <- as.integer(cls_tab[[cl]])
        }
      } else {
        for (cl in all_classes) {
          row_data[[paste0("count_", cl)]] <- 0L
        }
      }
    }
    as.data.frame(row_data, stringsAsFactors = FALSE)
  })
  summary_df <- do.call(rbind, summary_list)

  # 6. Overall class summary
  class_summary <- if (nrow(predictions) > 0) {
    tab <- table(predictions$label)
    cls_df <- data.frame(
      label = names(tab),
      total_count = as.integer(tab),
      mean_per_img = round(as.integer(tab) / n_items, 2),
      pct = round(as.integer(tab) / nrow(predictions) * 100, 2),
      stringsAsFactors = FALSE
    )
    sc_means <- tapply(predictions$score, predictions$label, mean, na.rm = TRUE)
    cls_df$mean_score <- round(as.numeric(sc_means[cls_df$label]), 4)
    cls_df[order(cls_df$total_count, decreasing = TRUE), , drop = FALSE]
  } else {
    data.frame(
      label = character(0), total_count = integer(0),
      mean_per_img = numeric(0), pct = numeric(0), mean_score = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(class_summary) <- NULL

  elapsed_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)

  out <- structure(
    list(
      predictions = predictions,
      summary = summary_df,
      class_summary = class_summary,
      results = results,
      n_images = n_items,
      elapsed = elapsed_sec,
      parallel = parallel
    ),
    class = c("dl_detect_batch", "list")
  )

  if (isTRUE(verbose)) {
    print(out)
  }

  invisible(out)
}

#' @export
print.dl_detect_batch <- function(x, ...) {
  cli::cli_rule(
    left = cli::col_green("✔ pliman Deep Learning Batch Detection Finished"),
    right = cli::col_blue(paste0(x$elapsed, "s"))
  )
  cli::cli_alert_success(
    "Processed {.val {x$n_images}} image{?s} | Found {.val {nrow(x$predictions)}} total object{?s} across {.val {nrow(x$class_summary)}} class{?es}."
  )
  if (nrow(x$class_summary) > 0) {
    cli::cli_h3("Class Summary:")
    print(x$class_summary, row.names = FALSE)
  }
  cli::cli_alert_info("Per-image results available in {.code $summary} and box coordinates in {.code $predictions}.")
  invisible(x)
}


#' Batch Instance and Semantic Segmentation with Deep Learning
#'
#' Runs deep learning segmentation models (e.g., YOLO-seg, Grounded-SAM, U2-Net, StarDist)
#' on a batch of images defined by a directory path, a filename pattern, a vector of file paths,
#' or a list of `Image` objects. Supports sequential or parallel processing (via `mirai`),
#' optional automated morphometric feature extraction (via [dl_extract_phenotypes()]),
#' and returns a consolidated predictions dataset along with per-image summaries.
#'
#' @param data Input data. Can be:
#'   * A character string representing a directory path containing images (e.g., `"flax_grains"`).
#'   * A character pattern to match image files in `path` or the current directory (e.g., `"-"`, `"^img_"`).
#'   * A character vector of image file paths.
#'   * A list of `Image` objects or numeric arrays.
#' @param pattern Optional regex pattern to filter image filenames when `data` is a directory path.
#' @param path Optional directory path to search when `data` is a pattern. Defaults to `"."`.
#' @param parallel Logical. If `TRUE`, runs inference in parallel across images using `mirai`.
#'   Defaults to `FALSE`.
#' @param workers Integer. Number of parallel workers/daemons to launch when `parallel = TRUE`.
#'   Defaults to 60% of available CPU cores.
#' @param measure Logical. If `TRUE`, automatically extracts granular morphometric, spectral,
#'   and geometric phenotypes for every segmented instance using [dl_extract_phenotypes()]. Defaults to `TRUE`.
#' @param pixel_size Optional numeric physical size per pixel (e.g., `0.05` for 0.05 cm/px) passed to [dl_extract_phenotypes()].
#' @param show_point Logical. Whether to show a point at the centroid of each segmented object.
#'   Defaults to `FALSE`.
#' @param keep_images Logical. Whether to keep the full original images attached as attributes in `$results`.
#'   Defaults to `FALSE` to prevent excessive RAM consumption during large batches (saving gigabytes of memory).
#' @param save_image Logical. Whether to save segmented or highlight mask images to disk.
#'   Defaults to `FALSE`.
#' @param dir_processed Character string. Directory where segmented images should be saved when
#'   `save_image = TRUE`. Defaults to `"./proc"`.
#' @param prefix Character string. Prefix prepended to saved output filenames. Defaults to `"proc_"`.
#' @param verbose Logical. Whether to show progress bars, timing, and formatted CLI summaries.
#'   Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_segment_dl()] (e.g., `model`, `type`, `threshold`,
#'   `box_threshold`, `text_threshold`, `iou_threshold`, `prompt`, `col_background`, `col_highlight`,
#'   `alpha`, `bbox`, `mask`, `engine`, `device_id`).
#'
#' @return An S3 object of class `c("dl_segment_batch", "list")` containing:
#'   * `predictions`: A consolidated `data.frame` of all detected bounding boxes or instance centroids across all images.
#'   * `phenotypes`: (If `measure = TRUE` and instance segmentation was run) A consolidated `data.frame` of all
#'     morphological, spectral, and shape traits with column `img`.
#'   * `summary`: A per-image `data.frame` summarizing total objects detected, total segmented area, and average size.
#'   * `results`: A named list of raw individual segmentation outputs from [image_segment_dl()].
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Batch instance segmentation with YOLO26-seg (Parallel)
#' res <- image_segment_dl_batch(
#'   data = "flax_grains",
#'   model = "yolo26n-seg",
#'   parallel = TRUE,
#'   workers = 4,
#'   measure = TRUE
#' )
#' head(res$summary)
#' head(res$phenotypes)
image_segment_dl_batch <- function(data,
                                   pattern = NULL,
                                   path = NULL,
                                   parallel = FALSE,
                                   workers = NULL,
                                   measure = TRUE,
                                   pixel_size = NULL,
                                   show_point = FALSE,
                                   keep_images = FALSE,
                                   save_image = FALSE,
                                   dir_processed = "./proc",
                                   prefix = "proc_",
                                   verbose = TRUE,
                                   ...) {
  # 1. Resolve image inputs
  resolved <- .resolve_batch_input(data = data, pattern = pattern, path = path)
  items <- resolved$items
  names_items <- resolved$names
  n_items <- length(items)

  if (n_items == 0L) {
    cli::cli_abort("No valid images found for batch segmentation.")
  }

  dots_args <- list(...)
  dots_args$plot <- FALSE
  dots_args$verbose <- FALSE
  dots_args$show_point <- show_point
  dots_args$attach_image <- isTRUE(keep_images)

  t0 <- Sys.time()

  if (isTRUE(verbose)) {
    cli::cli_rule(
      left = cli::col_blue("pliman: Deep Learning Batch Segmentation"),
      right = cli::col_blue(format(t0, "%Y-%m-%d | %H:%M:%S"))
    )
    cli::cli_alert_info("Found {.val {n_items}} image{?s} to segment.")
  }

  worker_seg_fn <- function(im, args) {
    im_obj <- if (is.character(im) && file.exists(im[1])) {
      pliman::image_import(im[1])
    } else {
      im
    }
    do.call(pliman::image_segment_dl, c(list(img = im_obj), args))
  }

  # 2. Execution: Parallel vs Sequential
  if (isTRUE(parallel)) {
    if (!requireNamespace("mirai", quietly = TRUE)) {
      cli::cli_abort("Package {.val mirai} is required for parallel processing. Please install it with {.code install.packages('mirai')}.")
    }

    nworkers <- if (is.null(workers)) {
      max(1L, floor(parallel::detectCores() * 0.6))
    } else {
      as.integer(workers)
    }
    eff_workers <- min(nworkers, n_items)

    st <- tryCatch(mirai::status()$daemons, error = function(e) 0L)
    existing_daemons <- if (is.numeric(st)) st else 0L

    daemon_status <- tryCatch({
      if (existing_daemons == 0L || (!is.null(workers) && existing_daemons != eff_workers)) {
        empty_prof <- tempfile()
        file.create(empty_prof)
        old_prof <- Sys.getenv("R_PROFILE_USER", unset = NA)
        Sys.setenv(R_PROFILE_USER = empty_prof)

        mirai::daemons(eff_workers)

        if (is.na(old_prof)) Sys.unsetenv("R_PROFILE_USER") else Sys.setenv(R_PROFILE_USER = old_prof)
        unlink(empty_prof)

        if (existing_daemons == 0L) {
          on.exit(mirai::daemons(0), add = TRUE)
        }
        mirai::everywhere(
          {
            .libPaths(lp)
            library(pliman)
          },
          lp = .libPaths()
        )
      }
      TRUE
    }, error = function(e) {
      cli::cli_warn("Failed to initialize parallel daemons: {e$message}. Falling back to sequential execution.")
      FALSE
    })

    if (isTRUE(daemon_status)) {
      if (isTRUE(verbose)) {
        cli::cli_rule(
          left = cli::col_blue("Parallel segmentation using {.val {eff_workers}} cores"),
          right = cli::col_blue(format(Sys.time(), "%H:%M:%S"))
        )
      }

      n_chunks <- if (n_items <= eff_workers) eff_workers else min(n_items, eff_workers * 3L)
      chunk_size <- ceiling(n_items / n_chunks)
      chunk_factor <- rep(seq_len(n_chunks), each = chunk_size, length.out = n_items)
      chunks <- split(items, chunk_factor)

      results_chunks <- mirai::mirai_map(
        .x = chunks,
        .f = function(chunk_items, fn, args) {
          lapply(chunk_items, function(it) {
            fn(it, args)
          })
        },
        .args = list(fn = worker_seg_fn, args = dots_args)
      )[.progress]

      failed_chunks <- sapply(results_chunks, function(res) inherits(res, "errorValue") || !is.list(res))
      if (any(failed_chunks)) {
        err_msg <- paste(sapply(which(failed_chunks), function(idx) {
          if (inherits(results_chunks[[idx]], "errorValue")) as.character(results_chunks[[idx]]) else "Worker failure"
        }), collapse = "; ")
        cli::cli_abort("Parallel segmentation failed: {err_msg}")
      }

      results <- do.call(c, results_chunks)
      if (!isTRUE(keep_images)) {
        for (i in seq_len(n_items)) {
          if (!is.null(results[[i]])) {
            attr(results[[i]], "image") <- NULL
            attr(results[[i]], "similarity_map") <- NULL
          }
        }
        gc(verbose = FALSE)
      }
    } else {
      parallel <- FALSE
    }
  }

  if (!isTRUE(parallel)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        total = n_items,
        clear = FALSE
      )
    }
    results <- vector("list", n_items)
    for (i in seq_len(n_items)) {
      if (isTRUE(verbose)) cli::cli_progress_update()
      res_i <- worker_seg_fn(items[[i]], dots_args)
      if (!isTRUE(keep_images) && !is.null(res_i)) {
        attr(res_i, "image") <- NULL
        attr(res_i, "similarity_map") <- NULL
      }
      results[[i]] <- res_i
      if (i %% 50L == 0L) gc(verbose = FALSE)
    }
    if (isTRUE(verbose)) cli::cli_progress_done()
  }

  names(results) <- names_items

  # 3. Optional saving of segmented output images
  if (isTRUE(save_image)) {
    if (!dir.exists(dir_processed)) {
      dir.create(dir_processed, recursive = TRUE)
    }
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Saving segmented images to {.path {dir_processed}}...",
        msg_done = "Segmented images successfully saved"
      )
    }
    for (i in seq_len(n_items)) {
      res_i <- results[[i]]
      out_img <- if (inherits(res_i, c("Image", "image"))) res_i else attr(res_i, "image")
      if (is.null(out_img) && is.character(items[[i]])) {
        out_img <- tryCatch(image_import(items[[i]]), error = function(e) NULL)
      }
      if (!is.null(out_img)) {
        out_file <- file.path(dir_processed, paste0(prefix, names_items[i], if (!grepl("\\.[a-zA-Z]+$", names_items[i])) ".png" else ""))
        tryCatch({
          # If res_i is a highlight list or show_point is requested, render onto png canvas
          is_highlight_res <- !inherits(res_i, c("Image", "image")) && is.list(res_i)
          if (is_highlight_res || isTRUE(show_point)) {
            dims <- dim(image_data(out_img))
            grDevices::png(out_file, width = dims[1], height = dims[2])
            on.exit(grDevices::dev.off(), add = TRUE)
            plot(out_img)
            bx <- if (is.data.frame(res_i)) res_i else if (!is.null(attr(res_i, "boxes"))) attr(res_i, "boxes") else res_i$boxes
            num_inst <- if (is.data.frame(bx)) nrow(bx) else 0L
            if (num_inst > 0L) {
              pal <- .get_class_palette(
                labels = bx$label,
                n = num_inst,
                rainbow = if (!is.null(dots_args$rainbow)) dots_args$rainbow else TRUE,
                user_col = dots_args$col
              )
              conts <- if (!is.null(attr(res_i, "contours"))) attr(res_i, "contours") else res_i$contours
              draw_mask <- if (!is.null(dots_args$mask)) isTRUE(dots_args$mask) else TRUE
              if (draw_mask && is.list(conts) && length(conts) > 0L) {
                alpha_val <- if (!is.null(dots_args$alpha)) dots_args$alpha else 0.4
                for (k in seq_len(min(num_inst, length(conts)))) {
                  ck <- conts[[k]]
                  if (is.null(ck) || length(ck) == 0L) next
                  k_col <- pal[((k - 1) %% length(pal)) + 1]
                  c_rgb <- tryCatch(grDevices::col2rgb(k_col) / 255.0, error = function(e) c(0, 1, 0))
                  poly_col <- grDevices::rgb(c_rgb[1], c_rgb[2], c_rgb[3], alpha = alpha_val)
                  p_list <- if (is.matrix(ck)) list(ck) else if (is.list(ck)) ck else list()
                  for (p in p_list) {
                    if (is.matrix(p) && nrow(p) >= 3) {
                      graphics::polygon(p[, 1], p[, 2], col = poly_col, border = NA)
                    }
                  }
                }
              }
              .plot_yolo_bboxes(
                boxes = bx,
                palette_colors = pal,
                lwd = if (!is.null(dots_args$lwd)) dots_args$lwd else 2,
                cex = if (!is.null(dots_args$cex)) dots_args$cex else 1.0,
                pad = if (!is.null(dots_args$pad)) dots_args$pad else 1.0,
                show_text = if (!is.null(dots_args$show_text)) dots_args$show_text else TRUE,
                show_conf = if (!is.null(dots_args$show_conf)) dots_args$show_conf else TRUE,
                show_class = if (!is.null(dots_args$show_class)) dots_args$show_class else TRUE,
                show_id = if (!is.null(dots_args$show_id)) dots_args$show_id else FALSE,
                show_box = if (!is.null(dots_args$show_box)) dots_args$show_box else (if (!is.null(dots_args$bbox)) dots_args$bbox else FALSE),
                show_point = show_point,
                point_size = if (!is.null(dots_args$point_size)) dots_args$point_size else 1.5,
                point_pch = if (!is.null(dots_args$point_pch)) dots_args$point_pch else 16,
                point_col = dots_args$point_col,
                badge = if (!is.null(dots_args$badge)) dots_args$badge else TRUE
              )
            }
            grDevices::dev.off()
          } else {
            image_export(out_img, out_file)
          }
        }, error = function(e) NULL)
        rm(out_img)
        if (i %% 50L == 0L) gc(verbose = FALSE)
      }
    }
    if (!isTRUE(keep_images)) {
      for (i in seq_len(n_items)) {
        if (!is.null(results[[i]])) {
          attr(results[[i]], "image") <- NULL
          attr(results[[i]], "similarity_map") <- NULL
        }
      }
      gc(verbose = FALSE)
    }
  }

  # 4. Consolidate bounding boxes / predictions
  pred_list <- lapply(seq_len(n_items), function(i) {
    res_i <- results[[i]]
    bx <- if (is.data.frame(res_i)) res_i else attr(res_i, "boxes")
    if (is.data.frame(bx) && nrow(bx) > 0) {
      data.frame(img = names_items[i], bx, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  predictions <- do.call(rbind, pred_list)
  if (is.null(predictions)) {
    predictions <- data.frame(
      img = character(0), id = integer(0), xmin = numeric(0), ymin = numeric(0),
      xmax = numeric(0), ymax = numeric(0), label = character(0), score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # 5. Extract phenotypes if measure = TRUE
  pheno_list <- list()
  if (isTRUE(measure)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Extracting morphometric phenotypes from segmented instances...",
        msg_done = "Morphometric feature extraction complete"
      )
    }
    for (i in seq_len(n_items)) {
      res_i <- results[[i]]
      ph <- tryCatch(
        dl_extract_phenotypes(res_i, pixel_size = pixel_size, verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(ph) && is.data.frame(ph) && nrow(ph) > 0) {
        pheno_list[[i]] <- data.frame(img = names_items[i], ph, stringsAsFactors = FALSE)
      }
    }
  }
  phenotypes <- if (length(pheno_list) > 0) do.call(rbind, pheno_list) else NULL

  # 6. Build per-image summary
  summary_list <- lapply(seq_len(n_items), function(i) {
    im_name <- names_items[i]
    res_i <- results[[i]]
    bx <- if (is.data.frame(res_i)) res_i else attr(res_i, "boxes")
    n_obj <- if (is.data.frame(bx)) nrow(bx) else 0L

    ph_i <- if (length(pheno_list) >= i && !is.null(pheno_list[[i]])) pheno_list[[i]] else NULL
    tot_area <- if (!is.null(ph_i) && "area" %in% names(ph_i)) round(sum(ph_i$area, na.rm = TRUE), 2) else NA_real_
    mean_area <- if (!is.null(ph_i) && "area" %in% names(ph_i) && nrow(ph_i) > 0) round(mean(ph_i$area, na.rm = TRUE), 2) else NA_real_

    row_data <- list(img = im_name, num_objects = n_obj, total_area = tot_area, mean_area = mean_area)
    as.data.frame(row_data, stringsAsFactors = FALSE)
  })
  summary_df <- do.call(rbind, summary_list)

  elapsed_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)

  out <- structure(
    list(
      predictions = predictions,
      phenotypes = phenotypes,
      summary = summary_df,
      results = results,
      n_images = n_items,
      elapsed = elapsed_sec,
      parallel = parallel
    ),
    class = c("dl_segment_batch", "list")
  )

  if (isTRUE(verbose)) {
    print(out)
  }

  invisible(out)
}

#' @export
print.dl_segment_batch <- function(x, ...) {
  cli::cli_rule(
    left = cli::col_green("✔ pliman Deep Learning Batch Segmentation Finished"),
    right = cli::col_blue(paste0(x$elapsed, "s"))
  )
  cli::cli_alert_success(
    "Processed {.val {x$n_images}} image{?s} | Segmented {.val {nrow(x$predictions)}} object instance{?s}."
  )
  if (!is.null(x$phenotypes) && nrow(x$phenotypes) > 0) {
    cli::cli_alert_info("Extracted {.val {ncol(x$phenotypes) - 1L}} phenotypic traits per object in {.code $phenotypes}.")
  }
  cli::cli_alert_info("Per-image area and counts available in {.code $summary}.")
  invisible(x)
}


# ==============================================================================
# Helper for resolving batch image inputs
# ==============================================================================
.resolve_batch_input <- function(data, pattern = NULL, path = NULL) {
  ext_pat <- "\\.(png|jpe?g|tiff?|webp|bmp|PNG|JPE?G|TIFF?|WEBP|BMP)$"

  # 1. List of Image objects
  if (is.list(data) && !is.character(data)) {
    names_img <- names(data)
    if (is.null(names_img) || any(names_img == "")) {
      names_img <- paste0("img_", seq_along(data))
      names(data) <- names_img
    }
    return(list(type = "images", items = data, names = names_img))
  }

  # 2. Character directory or single string
  if (is.character(data) && length(data) == 1L) {
    # Directory exists?
    if (dir.exists(data)) {
      patt <- if (!is.null(pattern)) pattern else ext_pat
      fls <- list.files(path = data, pattern = patt, full.names = TRUE)
      fls <- fls[grepl(ext_pat, fls)]
      return(list(type = "files", items = fls, names = basename(fls)))
    }

    # Single existing image file?
    if (file.exists(data)) {
      return(list(type = "files", items = data, names = basename(data)))
    }

    # Otherwise, treat as pattern
    target_dir <- if (!is.null(path) && dir.exists(path)) path else "."
    fls <- list.files(path = target_dir, pattern = data, full.names = TRUE)
    fls <- fls[grepl(ext_pat, fls)]
    return(list(type = "files", items = fls, names = basename(fls)))
  }

  # 3. Vector of file paths
  if (is.character(data) && length(data) > 1L) {
    return(list(type = "files", items = data, names = basename(data)))
  }

  # 4. Pattern with path fallback
  if ((missing(data) || is.null(data)) && !is.null(pattern)) {
    target_dir <- if (!is.null(path) && dir.exists(path)) path else "."
    fls <- list.files(path = target_dir, pattern = pattern, full.names = TRUE)
    fls <- fls[grepl(ext_pat, fls)]
    return(list(type = "files", items = fls, names = basename(fls)))
  }

  cli::cli_abort("Unable to resolve {.arg data}. Please provide a folder path, file pattern, vector of file paths, or list of images.")
}
