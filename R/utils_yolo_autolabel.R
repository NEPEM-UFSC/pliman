# ==============================================================================
# FEW-SHOT EXEMPLAR ANNOTATION & AUTO-LABELING FOR YOLO DATASETS
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Interactive Few-Shot Exemplar Fitting for Object Detection and Segmentation
#' @name yolo_fewshot_fit
#' @description
#' `yolo_fewshot_fit()` implements an interactive few-shot visual learning engine powered by
#' Meta AI's Segment Anything 2.1 (SAM 2.1 / PerSAM). The user selects 3–5 sample instances
#' of an object (e.g. seeds, beans, logs, leaves, lesions, fruits) by clicking on an image.
#'
#' Using the PerSAM exemplar matching engine, the function computes visual prototypes from the
#' exemplars and segments all other matching instances across the entire image with high precision.
#'
#' The user can specify `task = "segment"` to obtain instance masks and polygonal boundaries
#' for YOLO segmentation datasets, or `task = "detect"` for oriented/axis-aligned bounding boxes.
#'
#' The resulting model object (`yolo_fewshot_model`) can be inspected, plotted, and directly
#' propagated to hundreds of unannotated images in a folder using [yolo_dataset_autolabel()].
#'
#' @param img An `image` object, file path, or matrix.
#' @param points Optional $(x, y)$ coordinates of exemplar points. If `NULL` (default),
#'   an interactive plotting window opens for the user to click on 3–5 representative objects.
#' @param label Character vector naming the target class(es). Can be a single class (e.g., `"bean"`, `"leaf"`, `"seed"`, `"tora"`)
#'   or a character vector of multiple classes (e.g., `c("leaf", "lesion")`, `c("healthy", "damaged")`). Defaults to `"object"`.
#'   When multiple classes are specified, `yolo_fewshot_fit()` interactively collects exemplars for each class in sequence,
#'   or uses the respective elements from `points` if a list is provided.
#' @param task Target task: `"segment"` (instance segmentation polygons and masks) or
#'   `"detect"` (bounding boxes). Defaults to `"segment"`.
#' @param conf_threshold Numeric scalar between 0 and 1 specifying detection confidence/similarity threshold. Defaults to `0.50` (or `0.25` for adaptive).
#' @param radius Optional estimated radius in pixels for target objects. If `NULL` (default),
#'   automatically detected from exemplar edge gradients.
#' @param min_dist Minimum distance in pixels between detected instances. Defaults to `16`.
#' @param min_area Minimum object area in pixels. If `NULL` (default), automatically calibrated
#'   from the exemplar scale (`0.20 * median_exemplar_area`) to filter out background noise.
#' @param max_area Maximum object area in pixels. If `NULL` (default), automatically calibrated
#'   (`3.80 * median_exemplar_area`).
#' @param method Segmentation method: `"sam"` (or `"persam"`, Meta Segment Anything 2.1 via ONNX)
#'   or `"adaptive"` (fast color-texture morphology). Defaults to `"sam"`.
#' @param engine Execution engine for SAM: `"cpu"` or `"gpu"`. Defaults to `"cpu"`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param rainbow Logical. If `TRUE` (default), plots instances in distinct rainbow colors..
#' @param plot Logical. Display the detection and segmentation overlay. Defaults to `TRUE`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#'
#' @return An object of class `yolo_fewshot_model` containing exemplar signatures,
#'   detected instances on the reference image, and parameters for batch propagation.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # 1. Fit exemplar model on reference photo (clicks or supplied coordinates)
#' modelo_tora <- yolo_fewshot_fit("toras/toras22.png", label = "tora", task = "segment", method = "sam")
#'
#' # 2. Inspect segmentation overlay on reference photo
#' plot(modelo_tora)
#'
#' # 3. Auto-label entire directory of 300 photos into a YOLO dataset
#' ds <- yolo_dataset_autolabel("toras_raw/", model = modelo_tora, output_dir = "dataset_toras_yolo")
#' }
yolo_fewshot_fit <- function(img,
                             points = NULL,
                             sample = NULL,
                             label = "object",
                             task = c("segment", "detect"),
                             conf_threshold = 0.50,
                             radius = NULL,
                             min_dist = 12,
                             min_area = NULL,
                             max_area = NULL,
                             max_objects = NULL,
                             method = c("sam", "persam", "adaptive"),
                             engine = c("gpu", "cpu"),
                             device_id = -1,
                             rainbow = TRUE,
                             plot = TRUE,
                             verbose = TRUE) {
  task <- match.arg(task)
  method <- match.arg(method)
  engine <- match.arg(engine)

  # Support multi-class few-shot fitting
  if (length(label) > 1L) {
    if (isTRUE(verbose)) {
      cli::cli_h1("Multi-Class Few-Shot Training (pliman)")
      cli::cli_alert_info("Fitting {length(label)} visual classes: {.val {label}}")
    }
    models <- vector("list", length(label))
    names(models) <- label
    for (i in seq_along(label)) {
      lbl <- label[i]
      if (isTRUE(verbose)) {
        cli::cli_h2("Class {i}/{length(label)}: {.val {lbl}}")
      }
      pts_i <- if (is.list(points) && lbl %in% names(points)) {
        points[[lbl]]
      } else if (is.list(points) && length(points) >= i && !is.data.frame(points)) {
        points[[i]]
      } else if (i == 1L && (is.matrix(points) || is.data.frame(points))) {
        points
      } else {
        NULL
      }
      mod_i <- yolo_fewshot_fit(
        img = img,
        points = pts_i,
        sample = sample,
        label = lbl,
        task = task,
        conf_threshold = conf_threshold,
        radius = radius,
        min_dist = min_dist,
        min_area = min_area,
        max_area = max_area,
        max_objects = max_objects,
        method = method,
        engine = engine,
        device_id = device_id,
        rainbow = rainbow,
        plot = FALSE,
        verbose = verbose
      )
      mod_i$class_id <- as.integer(i - 1L)
      if (!is.null(mod_i$detected_boxes) && nrow(mod_i$detected_boxes) > 0L) {
        mod_i$detected_boxes$class_id <- as.integer(i - 1L)
      }
      models[[i]] <- mod_i
    }

    comb_obj <- structure(
      list(
        label = label,
        classes = label,
        models = models,
        task = task,
        conf_threshold = conf_threshold,
        method = method,
        engine = engine,
        reference_image = models[[1]]$reference_image,
        reference_images = models[[1]]$reference_images
      ),
      class = c("yolo_fewshot_multimodel", "yolo_fewshot_model")
    )

    if (isTRUE(verbose)) {
      cli::cli_alert_success("Multi-class few-shot model successfully fitted for classes: {.val {label}}!")
    }

    if (isTRUE(plot)) {
      plot(comb_obj, rainbow = rainbow)
    }

    return(invisible(comb_obj))
  }

  # Resolve input image(s): directory with sample, vector of file paths, or single image
  img_list <- list()
  if (is.character(img) && length(img) == 1L && dir.exists(img[1])) {
    img_files <- list.files(img[1], pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
    if (length(img_files) == 0L) {
      cli::cli_abort("No image files found in directory {.path {img[1]}}.")
    }
    n_s <- if (!is.null(sample) && sample > 0L) min(length(img_files), as.integer(sample)) else min(3L, length(img_files))
    img_list <- as.list(sample(img_files, n_s))
  } else if (is.character(img) && length(img) > 1L) {
    if (!is.null(sample) && sample > 0L && sample < length(img)) {
      img_list <- as.list(sample(img, as.integer(sample)))
    } else {
      img_list <- as.list(img)
    }
  } else if (is.list(img) && !inherits(img, c("image", "Image"))) {
    img_list <- img
  } else {
    img_list <- list(img)
  }

  # Multi-Image Training Mode
  if (length(img_list) > 1L) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Multi-Image Few-Shot Training (pliman)")
      cli::cli_alert_info("Training visual prototypes across {length(img_list)} reference images for robust generalization.")
    }

    all_prototypes_sem <- list()
    all_prototypes_hr0 <- list()
    all_areas <- numeric()
    all_min_areas <- numeric()
    all_max_areas <- numeric()
    all_ref_boxes <- list()
    all_ref_polys <- list()
    all_ref_paths <- character()
    all_ref_imgs <- list()
    all_ref_points <- list()
    stop_early <- FALSE

    for (s in seq_along(img_list)) {
      curr_target <- img_list[[s]]
      curr_im <- if (is.character(curr_target) && file.exists(curr_target[1])) {
        image_import(curr_target[1])
      } else if (inherits(curr_target, c("image", "Image"))) {
        curr_target
      } else {
        as_image(curr_target)
      }

      c_dims <- dim(curr_im)
      cw <- c_dims[1]; ch <- c_dims[2]

      curr_points <- NULL
      if (is.list(points) && length(points) >= s) {
        curr_points <- points[[s]]
      } else if (is.matrix(points) && s == 1L) {
        curr_points <- points
      }

      if (is.null(curr_points)) {
        repeat {
          if (isTRUE(verbose)) {
            if (s == 1L || length(all_prototypes_sem) == 0L) {
              cli::cli_alert_info("Image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}): Click on 2 to 5 exemplar instances of {.val {label}} in the plot window.")
              cli::cli_alert_info("Press <Esc> or right-click when finished.")
            } else {
              cli::cli_alert_info("Image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}): Click on exemplar instances of {.val {label}}, or press <Esc> to finish or continue.")
            }
          }
          plot(curr_im)
          pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.3), error = function(e) NULL)

          if (is.null(pts) || length(pts$x) == 0L) {
            if (s == 1L && length(all_prototypes_sem) == 0L) {
              if (!interactive()) {
                cli::cli_abort("No exemplar points selected on first reference image.")
              }
              cli::cli_alert_warning("No point clicked on first image ({basename(as.character(curr_target[1]))}) (<Esc> pressed).")
              cli::cli_bullets(c(
                "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image",
                "*" = "{.bold [s] / stop}: Abort model fitting"
              ))
              ans <- tolower(trimws(readline("Select action [r = retry / s = stop, default: r]: ")))
              if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar")) {
                cli::cli_abort("Model fitting cancelled by user.")
              } else {
                next
              }
            } else {
              if (interactive()) {
                cli::cli_alert_warning("No point clicked on image {s}/{length(img_list)} ({basename(as.character(curr_target[1]))}) (<Esc> pressed).")
                cli::cli_bullets(c(
                  "*" = "{.bold [s] / stop}: Finalize model now with the {length(all_ref_paths)} reference image(s) completed so far",
                  "*" = "{.bold [c] / continue}: Continue to remaining images (auto-detect on this image)",
                  "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image"
                ))
                ans <- tolower(trimws(readline("Select action [s = stop / c = continue / r = retry, default: s]: ")))
                if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar", "")) {
                  cli::cli_alert_info("Finalizing model with the {length(all_ref_paths)} reference image(s) completed so far...")
                  stop_early <- TRUE
                  break
                } else if (ans %in% c("r", "retry", "tentar", "t")) {
                  next
                } else {
                  curr_points <- NULL
                  break
                }
              } else {
                curr_points <- NULL
                break
              }
            }
          } else {
            curr_points <- cbind(x = pts$x, y = pts$y)
            break
          }
        }

        if (isTRUE(stop_early)) {
          break
        }
      } else if (is.vector(curr_points) && length(curr_points) >= 2L && length(curr_points) %% 2 == 0L) {
        half <- length(curr_points) / 2
        curr_points <- cbind(x = curr_points[seq_len(half)], y = curr_points[(half + 1L):length(curr_points)])
      } else if (is.data.frame(curr_points)) {
        curr_points <- as.matrix(curr_points[, 1:2])
      }

      c_num_arr <- image_data(curr_im, type = "numeric")
      pre_proto <- if (is.null(curr_points) && length(all_prototypes_sem) > 0L) {
        list(
          semantic = all_prototypes_sem,
          hr0 = all_prototypes_hr0,
          avg_area = if (length(all_areas) > 0) mean(all_areas) else 100.0,
          min_area = if (length(all_min_areas) > 0) min(all_min_areas) else 10.0,
          max_area = if (length(all_max_areas) > 0) max(all_max_areas) else 500.0
        )
      } else {
        NULL
      }

      p_res <- tryCatch({
        .run_persam(
          mat = c_num_arr,
          exemplar_points = curr_points,
          precomputed_prototypes = pre_proto,
          sim_threshold = conf_threshold,
          min_dist = min_dist,
          max_objects = max_objects,
          feat_res = 256L,
          engine = engine,
          device_id = device_id,
          fill_hull = TRUE,
          verbose = FALSE
        )
      }, error = function(e) {
        if (isTRUE(verbose)) {
          cli::cli_alert_warning("Exemplar matching error on image {s}: {e$message}")
        }
        NULL
      })

      if (!is.null(p_res) && !is.null(p_res$prototypes)) {
        all_prototypes_sem <- c(all_prototypes_sem, p_res$prototypes$semantic)
        all_prototypes_hr0 <- c(all_prototypes_hr0, p_res$prototypes$hr0)
        if (!is.null(p_res$prototypes$avg_area)) {
          all_areas <- c(all_areas, as.numeric(p_res$prototypes$avg_area))
        }
        if (!is.null(p_res$prototypes$min_area)) {
          all_min_areas <- c(all_min_areas, as.numeric(p_res$prototypes$min_area))
        }
        if (!is.null(p_res$prototypes$max_area)) {
          all_max_areas <- c(all_max_areas, as.numeric(p_res$prototypes$max_area))
        }
      }

      ref_id <- if (is.character(curr_target)) normalizePath(curr_target[1], winslash = "/", mustWork = FALSE) else paste0("ref_", s)
      all_ref_paths <- c(all_ref_paths, ref_id)
      all_ref_imgs[[ref_id]] <- curr_im
      if (!is.null(curr_points)) {
        all_ref_points[[ref_id]] <- curr_points
      }

      r_boxes <- data.frame()
      r_polys <- list()

      if (!is.null(p_res) && nrow(p_res$boxes) > 0L) {
        r_boxes <- data.frame(
          xmin = as.numeric(p_res$boxes$xmin),
          ymin = as.numeric(p_res$boxes$ymin),
          xmax = as.numeric(p_res$boxes$xmax),
          ymax = as.numeric(p_res$boxes$ymax),
          conf = round(as.numeric(p_res$boxes$score), 3),
          class_id = 0L,
          class_name = label,
          stringsAsFactors = FALSE
        )
        for (ci in seq_along(p_res$contours)) {
          cnt <- p_res$contours[[ci]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            gx <- cnt[, 1] / cw; gy <- cnt[, 2] / ch
            if (length(gx) > 32L) {
              idx <- round(seq(1, length(gx), length.out = 32))
              gx <- gx[idx]; gy <- gy[idx]
            }
            r_polys[[ci]] <- as.vector(rbind(gx, gy))
          }
        }
      }
      all_ref_boxes[[ref_id]] <- r_boxes
      all_ref_polys[[ref_id]] <- r_polys
    }

    if (length(all_ref_paths) == 0L) {
      cli::cli_abort("Model fitting aborted. No reference images were annotated.")
    }

    all_areas <- as.numeric(unlist(all_areas, use.names = FALSE))
    all_min_areas <- as.numeric(unlist(all_min_areas, use.names = FALSE))
    all_max_areas <- as.numeric(unlist(all_max_areas, use.names = FALSE))
    all_areas <- all_areas[!is.na(all_areas)]
    all_min_areas <- all_min_areas[!is.na(all_min_areas)]
    all_max_areas <- all_max_areas[!is.na(all_max_areas)]

    merged_prototypes <- list(
      semantic = all_prototypes_sem,
      hr0 = all_prototypes_hr0,
      avg_area = if (length(all_areas) > 0) mean(all_areas) else 100.0,
      min_area = if (length(all_min_areas) > 0) min(all_min_areas) else 10.0,
      max_area = if (length(all_max_areas) > 0) max(all_max_areas) else 500.0
    )

    primary_id <- all_ref_paths[1]
    primary_im <- all_ref_imgs[[primary_id]]
    primary_boxes <- all_ref_boxes[[primary_id]]
    primary_polys <- all_ref_polys[[primary_id]]

    model_obj <- structure(
      list(
        label = label,
        task = task,
        target_rgb = c(0.5, 0.5, 0.5),
        rgb_sd = c(0.2, 0.2, 0.2),
        median_radius = 20,
        median_area = 1200,
        conf_threshold = conf_threshold,
        min_area = if (!is.null(min_area)) min_area else 0,
        max_area = max_area,
        exemplar_points = if (length(all_ref_points) > 0L) all_ref_points[[1]] else NULL,
        exemplar_points_map = all_ref_points,
        detected_boxes = primary_boxes,
        detected_polygons = primary_polys,
        detected_labels = NULL,
        detected_mask = NULL,
        similarity_map = NULL,
        prototypes = merged_prototypes,
        reference_image_path = primary_id,
        reference_image_paths = all_ref_paths,
        reference_boxes_map = all_ref_boxes,
        reference_polys_map = all_ref_polys,
        reference_image = primary_im,
        reference_images = all_ref_imgs,
        method = "sam",
        engine = engine
      ),
      class = "yolo_fewshot_model"
    )

    if (isTRUE(verbose)) {
      total_det <- sum(vapply(all_ref_boxes, nrow, integer(1)))
      cli::cli_alert_success("Model fitted across {length(all_ref_paths)} images! Learned {length(all_prototypes_sem)} prototypes ({total_det} reference instances segmented).")
    }

    if (isTRUE(plot) && !is.null(primary_im)) {
      plot(model_obj, rainbow = rainbow)
    }

    return(invisible(model_obj))
  }

  # Single Image Training Mode
  img_single <- img_list[[1]]
  im <- if (is.character(img_single) && file.exists(img_single[1])) {
    image_import(img_single[1])
  } else if (inherits(img_single, c("image", "Image"))) {
    img_single
  } else {
    as_image(img_single)
  }

  dims <- dim(im)
  w <- dims[1]
  h <- dims[2]

  # Interactive clicking if points not provided
  if (is.null(points)) {
    if (isTRUE(verbose)) {
      cli::cli_h2("Few-Shot Exemplar Selection (pliman)")
      cli::cli_alert_info("Click on 3 to 5 sample instances of {.val {label}} in the plot window.")
      cli::cli_alert_info("Press <Esc> or right-click when finished.")
    }
    repeat {
      plot(im)
      pts <- tryCatch(graphics::locator(n = 64, type = "p", col = "#FFFF00", pch = 19, cex = 1.3), error = function(e) NULL)
      if (is.null(pts) || length(pts$x) == 0L) {
        if (!interactive()) {
          cli::cli_abort("No exemplar points selected.")
        }
        cli::cli_alert_warning("No point clicked on image (<Esc> pressed).")
        cli::cli_bullets(c(
          "*" = "{.bold [r] / retry}: Retry clicking exemplar objects on this image",
          "*" = "{.bold [s] / stop}: Abort model fitting"
        ))
        ans <- tolower(trimws(readline("Select action [r = retry / s = stop, default: r]: ")))
        if (ans %in% c("s", "stop", "abort", "a", "q", "quit", "p", "parar")) {
          cli::cli_abort("Model fitting cancelled by user.")
        } else {
          next
        }
      } else {
        points <- cbind(x = pts$x, y = pts$y)
        break
      }
    }
  } else if (is.vector(points) && length(points) >= 2L && length(points) %% 2 == 0L) {
    half <- length(points) / 2
    points <- cbind(x = points[seq_len(half)], y = points[(half + 1L):length(points)])
  } else if (is.data.frame(points)) {
    points <- as.matrix(points[, 1:2])
  }

  n_exemplars <- nrow(points)
  if (isTRUE(verbose)) {
    cli::cli_alert_success("Collected {n_exemplars} exemplar point{?s} for class {.val {label}}.")
  }

  norm_arr <- image_data(im, type = "normalized")
  num_arr <- image_data(im, type = "numeric")

  # 1. Scale auto-calibration from exemplars
  angles <- seq(0, 2 * pi, length.out = 16)[-16]
  max_search <- min(140L, as.integer(round(min(w, h) * 0.20)))
  radii <- numeric(n_exemplars)

  for (k in seq_len(n_exemplars)) {
    px <- as.integer(round(points[k, 1]))
    py <- as.integer(round(points[k, 2]))
    c_rgb <- norm_arr[px, py, ]
    r_search <- 5:max_search
    diffs <- numeric(length(r_search))
    for (i in seq_along(r_search)) {
      r <- r_search[i]
      xs <- pmax(1L, pmin(w, as.integer(round(px + r * cos(angles)))))
      ys <- pmax(1L, pmin(h, as.integer(round(py + r * sin(angles)))))
      d_rgb <- sqrt((norm_arr[cbind(xs, ys, 1)] - c_rgb[1])^2 +
                    (norm_arr[cbind(xs, ys, 2)] - c_rgb[2])^2 +
                    (norm_arr[cbind(xs, ys, 3)] - c_rgb[3])^2)
      diffs[i] <- mean(d_rgb)
    }
    grad <- c(0, diff(diffs))
    start_idx <- min(8L, length(grad))
    edge_r <- r_search[which.max(grad[start_idx:length(grad)]) + (start_idx - 1L)]
    radii[k] <- edge_r
  }

  median_rad <- if (!is.null(radius)) as.numeric(radius) else max(8, median(radii))
  median_area <- pi * median_rad^2

  if (is.null(min_area) || identical(min_area, 25)) {
    min_area <- 0
  }

  # Extract exemplar color signatures
  exemplar_colors <- list()
  for (k in seq_len(n_exemplars)) {
    px <- as.integer(round(points[k, 1]))
    py <- as.integer(round(points[k, 2]))
    rk <- max(3L, as.integer(round(radii[k] * 0.35)))
    sub_r <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 1]
    sub_g <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 2]
    sub_b <- norm_arr[max(1L, px - rk):min(w, px + rk), max(1L, py - rk):min(h, py + rk), 3]
    exemplar_colors[[k]] <- c(mean(sub_r), mean(sub_g), mean(sub_b))
  }
  col_mat <- do.call(rbind, exemplar_colors)
  target_rgb <- colMeans(col_mat)
  rgb_sd <- pmax(apply(col_mat, 2, stats::sd), 0.05)

  # 2. Execute segmentation engine
  used_method <- method
  df_boxes <- data.frame()
  cand_polys <- list()
  det_labels <- NULL
  det_mask <- NULL
  det_sim <- NULL

  if (method %in% c("sam", "persam")) {
    used_method <- "sam"
    persam_res <- tryCatch({
      .run_persam(
        mat = num_arr,
        exemplar_points = points,
        sim_threshold = conf_threshold,
        min_dist = if (!is.null(radius)) min(min_dist, radius * 0.4) else min_dist,
        max_objects = max_objects,
        feat_res = 256L,
        engine = engine,
        device_id = device_id,
        fill_hull = TRUE,
        verbose = verbose
      )
    }, error = function(e) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning("PerSAM execution encountered: {e$message}. Falling back to adaptive morphology.")
      }
      NULL
    })

    if (!is.null(persam_res) && nrow(persam_res$boxes) > 0L) {
      raw_boxes <- persam_res$boxes
      raw_conts <- persam_res$contours
      det_labels <- persam_res$labels
      det_mask <- persam_res$mask
      det_sim <- persam_res$similarity_map

      # Filter by scale if specified
      valid_idx <- seq_len(nrow(raw_boxes))
      if (!is.null(min_area) && min_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- which(b_areas >= min_area)
      }
      if (!is.null(max_area) && max_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- intersect(valid_idx, which(b_areas <= max_area))
      }

      if (length(valid_idx) > 0L) {
        raw_boxes <- raw_boxes[valid_idx, , drop = FALSE]
        raw_conts <- raw_conts[valid_idx]

        df_boxes <- data.frame(
          xmin = as.numeric(raw_boxes$xmin),
          ymin = as.numeric(raw_boxes$ymin),
          xmax = as.numeric(raw_boxes$xmax),
          ymax = as.numeric(raw_boxes$ymax),
          conf = round(as.numeric(raw_boxes$score), 3),
          class_id = 0L,
          class_name = label,
          stringsAsFactors = FALSE
        )

        for (i in seq_along(raw_conts)) {
          cnt <- raw_conts[[i]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            gx <- cnt[, 1] / w
            gy <- cnt[, 2] / h
            if (length(gx) > 32L) {
              idx <- round(seq(1, length(gx), length.out = 32))
              gx <- gx[idx]; gy <- gy[idx]
            }
            cand_polys[[i]] <- as.vector(rbind(gx, gy))
          } else {
            cand_polys[[i]] <- NULL
          }
        }
      }
    } else {
      used_method <- "adaptive"
    }
  }

  if (identical(used_method, "adaptive") || nrow(df_boxes) == 0L) {
    used_method <- "adaptive"
    # Adaptive Gaussian Mahalanobis similarity map
    dr <- (norm_arr[, , 1] - target_rgb[1]) / rgb_sd[1]
    dg <- (norm_arr[, , 2] - target_rgb[2]) / rgb_sd[2]
    db <- (norm_arr[, , 3] - target_rgb[3]) / rgb_sd[3]
    dist_sq <- (dr^2 + dg^2 + db^2) / 3.0
    sim_raw <- exp(-0.5 * dist_sq)

    blur_sigma <- max(2, as.integer(round(median_rad * 0.12)))
    sim_im <- as_image(matrix(sim_raw, w, h), colormode = "Grayscale")
    sim_blur <- image_blur(sim_im, sigma = blur_sigma)
    sim_mat <- image_data(sim_blur, type = "normalized")

    adapt_t <- if (conf_threshold > 0.4) 0.25 else conf_threshold
    bin_mask <- sim_mat >= adapt_t
    lbl <- bwlabel_cpp(bin_mask)
    counts <- table(lbl[lbl > 0])
    valid <- as.integer(names(counts[counts >= min_area & counts <= max_area]))

    cand_b <- list()
    cand_c <- numeric()
    conts_all <- if (length(valid) > 0L) extract_contours_cpp(lbl) else list()
    cand_p <- list()

    for (obj in valid) {
      coords <- which(lbl == obj, arr.ind = TRUE)
      bx1 <- min(coords[, 1]); bx2 <- max(coords[, 1])
      by1 <- min(coords[, 2]); by2 <- max(coords[, 2])
      bw <- bx2 - bx1 + 1L; bh <- by2 - by1 + 1L
      ar <- bw / pmax(1L, bh)
      if (ar < 0.38 || ar > 2.60) next

      if (bw >= 12L && bh >= 12L) {
        sub_p <- norm_arr[bx1:bx2, by1:by2, 1]
        nr <- nrow(sub_p); nc <- ncol(sub_p)
        mr <- max(2L, as.integer(round(nr * 0.15))); mc <- max(2L, as.integer(round(nc * 0.15)))
        c_sub <- sub_p[mr:(nr - mr), mc:(nc - mc)]
        b_sub <- c(sub_p[1:mr, ], sub_p[(nr - mr + 1L):nr, ],
                   sub_p[, 1:mc], sub_p[, (nc - mc + 1L):nc])
        if (mean(c_sub) - mean(b_sub) < 0.035) next
      }

      obj_conf <- mean(sim_mat[coords])
      if (obj_conf < adapt_t) next

      poly_norm <- NULL
      if (length(conts_all) >= obj && is.matrix(conts_all[[obj]]) && nrow(conts_all[[obj]]) >= 3L) {
        poly_pts <- conts_all[[obj]]
        gx <- poly_pts[, 1] / w; gy <- poly_pts[, 2] / h
        if (length(gx) > 32L) {
          idx <- round(seq(1, length(gx), length.out = 32))
          gx <- gx[idx]; gy <- gy[idx]
        }
        poly_norm <- as.vector(rbind(gx, gy))
      }

      cand_b[[length(cand_b) + 1L]] <- c(bx1, by1, bx2, by2)
      cand_c <- c(cand_c, obj_conf)
      cand_p[[length(cand_p) + 1L]] <- poly_norm
    }

    if (length(cand_b) > 0L) {
      box_mat <- do.call(rbind, cand_b)
      if (nrow(box_mat) > 1L) {
        keep_idx <- .nms_boxes(box_mat, cand_c, iou_thresh = 0.35)
        box_mat <- box_mat[keep_idx, , drop = FALSE]
        cand_c <- cand_c[keep_idx]
        cand_p <- cand_p[keep_idx]
      }
      df_boxes <- data.frame(
        xmin = as.numeric(box_mat[, 1]),
        ymin = as.numeric(box_mat[, 2]),
        xmax = as.numeric(box_mat[, 3]),
        ymax = as.numeric(box_mat[, 4]),
        conf = round(cand_c, 3),
        class_id = 0L,
        class_name = label,
        stringsAsFactors = FALSE
      )
      cand_polys <- cand_p
    }
  }

  ref_path <- if (is.character(img) && file.exists(img[1])) normalizePath(img[1], winslash = "/", mustWork = FALSE) else NULL
  prototypes <- if (exists("persam_res") && !is.null(persam_res)) persam_res$prototypes else NULL

  ref_boxes_map <- list()
  ref_polys_map <- list()
  if (!is.null(ref_path)) {
    ref_boxes_map[[ref_path]] <- df_boxes
    ref_polys_map[[ref_path]] <- cand_polys
  }

  model_obj <- structure(
    list(
      label = label,
      task = task,
      target_rgb = target_rgb,
      rgb_sd = rgb_sd,
      median_radius = median_rad,
      median_area = median_area,
      conf_threshold = conf_threshold,
      min_area = min_area,
      max_area = max_area,
      exemplar_points = points,
      exemplar_points_map = if (!is.null(ref_path)) stats::setNames(list(points), ref_path) else list(points),
      detected_boxes = df_boxes,
      detected_polygons = cand_polys,
      detected_labels = det_labels,
      detected_mask = det_mask,
      similarity_map = det_sim,
      prototypes = prototypes,
      reference_image_path = ref_path,
      reference_image_paths = if (!is.null(ref_path)) ref_path else character(0),
      reference_boxes_map = ref_boxes_map,
      reference_polys_map = ref_polys_map,
      reference_image = im,
      reference_images = if (!is.null(ref_path)) stats::setNames(list(im), ref_path) else list(im),
      method = used_method,
      engine = engine
    ),
    class = "yolo_fewshot_model"
  )

  if (isTRUE(verbose)) {
    n_det <- if (task == "segment" && length(cand_polys) > 0L) sum(!vapply(cand_polys, is.null, logical(1))) else nrow(df_boxes)
    cli::cli_alert_success("Model fitted! {n_det} instances {if (task == 'segment') 'segmented' else 'detected'} on reference image.")
  }

  if (isTRUE(plot)) {
    plot(model_obj, rainbow = rainbow)
  }

  invisible(model_obj)
}

#' @title Combine Few-Shot Models
#' @description Combines multiple `yolo_fewshot_model` objects into a multi-class few-shot model.
#' @param ... `yolo_fewshot_model` objects to combine.
#' @return A `yolo_fewshot_multimodel` object.
#' @export
c.yolo_fewshot_model <- function(...) {
  dots <- list(...)
  dots <- dots[!vapply(dots, is.null, logical(1))]
  all_submodels <- list()
  for (item in dots) {
    if (inherits(item, "yolo_fewshot_multimodel")) {
      all_submodels <- c(all_submodels, item$models)
    } else if (inherits(item, "yolo_fewshot_model")) {
      all_submodels <- c(all_submodels, list(item))
    }
  }
  if (length(all_submodels) == 0L) {
    return(NULL)
  }
  cls_names <- character(length(all_submodels))
  for (i in seq_along(all_submodels)) {
    c_name <- all_submodels[[i]]$label[1]
    cls_names[i] <- c_name
    all_submodels[[i]]$class_id <- as.integer(i - 1L)
    if (!is.null(all_submodels[[i]]$detected_boxes) && nrow(all_submodels[[i]]$detected_boxes) > 0L) {
      all_submodels[[i]]$detected_boxes$class_id <- as.integer(i - 1L)
    }
  }
  names(all_submodels) <- cls_names
  structure(
    list(
      label = cls_names,
      classes = cls_names,
      models = all_submodels,
      task = all_submodels[[1]]$task,
      conf_threshold = all_submodels[[1]]$conf_threshold,
      method = all_submodels[[1]]$method,
      engine = all_submodels[[1]]$engine,
      reference_image = all_submodels[[1]]$reference_image,
      reference_images = all_submodels[[1]]$reference_images
    ),
    class = c("yolo_fewshot_multimodel", "yolo_fewshot_model")
  )
}

#' @title Predict Few-Shot Instances on a New Image
#' @name yolo_fewshot_predict
#' @description
#' Applies a trained `yolo_fewshot_model` to detect and segment instances in a new image.
#'
#' @param model A `yolo_fewshot_model` fitted via [yolo_fewshot_fit()].
#' @param img An `image` object or file path to evaluate.
#' @param task Task override: `"segment"` (polygons) or `"detect"` (bounding boxes). Defaults to `model$task`.
#' @param conf_threshold Numeric scalar override for detection confidence. Defaults to model threshold.
#' @param max_objects Maximum number of instances to detect.
#' @param engine Execution engine: `"cpu"` or `"gpu"`. If `NULL`, uses `model$engine`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param plot Logical. Display detection results. Defaults to `FALSE`.
#' @param rainbow Logical. If `TRUE`, plots instances in distinct rainbow colors. Defaults to `FALSE`.
#'
#' @return A `list` containing `boxes` (bounding boxes data.frame) and `polygons` (list of normalized coordinates).
#' @export
yolo_fewshot_predict <- function(model,
                                img,
                                task = model$task,
                                conf_threshold = NULL,
                                max_objects = NULL,
                                engine = "gpu",
                                device_id = -1,
                                plot = FALSE,
                                rainbow = TRUE) {
  if (is.list(model) && !inherits(model, "yolo_fewshot_model")) {
    if (all(vapply(model, function(m) inherits(m, "yolo_fewshot_model"), logical(1)))) {
      model <- do.call(c, model)
    }
  }

  if (!inherits(model, "yolo_fewshot_model")) {
    cli::cli_abort("{.arg model} must be a {.cls yolo_fewshot_model} from {.fn yolo_fewshot_fit}.")
  }

  im <- if (is.character(img) && file.exists(img[1])) {
    image_import(img[1])
  } else if (inherits(img, c("image", "Image"))) {
    img
  } else {
    as_image(img)
  }

  dims <- dim(im)
  w <- dims[1]; h <- dims[2]

  task <- if (!missing(task) && !is.null(task)) match.arg(task, c("segment", "detect")) else if (!is.null(model$task)) model$task else "segment"
  conf_t <- if (!is.null(conf_threshold)) as.numeric(conf_threshold) else model$conf_threshold

  # Handle multi-class models
  if (inherits(model, "yolo_fewshot_multimodel") || !is.null(model$models)) {
    all_b <- list()
    all_p <- list()
    for (m_idx in seq_along(model$models)) {
      sub_m <- model$models[[m_idx]]
      sub_res <- yolo_fewshot_predict(
        model = sub_m,
        img = im,
        task = task,
        conf_threshold = conf_threshold,
        max_objects = max_objects,
        engine = engine,
        device_id = device_id,
        plot = FALSE,
        rainbow = FALSE
      )
      if (!is.null(sub_res$boxes) && nrow(sub_res$boxes) > 0L) {
        all_b[[length(all_b) + 1L]] <- sub_res$boxes
        all_p <- c(all_p, sub_res$polygons)
      }
    }
    comb_boxes <- if (length(all_b) > 0L) do.call(rbind, all_b) else data.frame()
    res <- list(boxes = comb_boxes, polygons = all_p, image = im)
    if (isTRUE(plot)) {
      plot(im)
      if (nrow(res$boxes) > 0L) {
        classes_vec <- if (!is.null(model$classes)) model$classes else unique(res$boxes$class_name)
        pal <- if (isTRUE(rainbow)) grDevices::rainbow(max(1L, length(classes_vec)), s = 0.85, v = 0.95) else rep("#00FFCC", length(classes_vec))
        names(pal) <- classes_vec
        if (task == "segment" && length(res$polygons) > 0L) {
          for (i in seq_along(res$polygons)) {
            poly <- res$polygons[[i]]
            c_name <- if (i <= nrow(res$boxes)) res$boxes$class_name[i] else classes_vec[1]
            c_base <- if (c_name %in% names(pal)) pal[[c_name]] else "#00FFCC"
            if (!is.null(poly) && length(poly) >= 6L) {
              px <- poly[c(TRUE, FALSE)] * w
              py <- poly[c(FALSE, TRUE)] * h
              graphics::polygon(px, py, col = grDevices::adjustcolor(c_base, alpha.f = 0.40), border = c_base, lwd = 1.5)
            }
          }
        }
        for (r in seq_len(nrow(res$boxes))) {
          c_name <- res$boxes$class_name[r]
          c_col <- if (c_name %in% names(pal)) pal[[c_name]] else "#00FFCC"
          if (task == "detect" || length(res$polygons) == 0L) {
            graphics::rect(res$boxes$xmin[r], res$boxes$ymin[r], res$boxes$xmax[r], res$boxes$ymax[r], border = c_col, lwd = 2)
          }
          graphics::text(res$boxes$xmin[r], res$boxes$ymin[r] - 5,
                         labels = sprintf("%s %.2f", res$boxes$class_name[r], res$boxes$conf[r]),
                         col = c_col, cex = 0.8, pos = 4)
        }
      }
    }
    return(invisible(res))
  }

  # 1. Check if this is one of the exact reference images used during fit
  is_ref <- FALSE
  ref_match_id <- NULL
  if (is.character(img)) {
    norm_in <- normalizePath(img[1], winslash = "/", mustWork = FALSE)
    if (!is.null(model$reference_image_paths) && norm_in %in% model$reference_image_paths) {
      is_ref <- TRUE
      ref_match_id <- norm_in
    } else if (!is.null(model$reference_image_path) && identical(norm_in, model$reference_image_path)) {
      is_ref <- TRUE
      ref_match_id <- norm_in
    }
  }
  if (!is_ref && !is.null(model$reference_image)) {
    ref_d <- dim(model$reference_image)
    if (identical(dims, ref_d)) {
      s_idx <- seq(1, length(im), length.out = min(100L, length(im)))
      if (isTRUE(all.equal(as.numeric(im)[s_idx], as.numeric(model$reference_image)[s_idx]))) {
        is_ref <- TRUE
      }
    }
  }

  if (is_ref) {
    b_match <- if (!is.null(ref_match_id) && !is.null(model$reference_boxes_map[[ref_match_id]])) {
      model$reference_boxes_map[[ref_match_id]]
    } else {
      model$detected_boxes
    }
    p_match <- if (!is.null(ref_match_id) && !is.null(model$reference_polys_map[[ref_match_id]])) {
      model$reference_polys_map[[ref_match_id]]
    } else {
      model$detected_polygons
    }

    if (!is.null(b_match) && nrow(b_match) > 0L) {
      if (!is.null(conf_t) && !is.null(b_match$conf)) {
        keep_ref <- which(b_match$conf >= conf_t)
        b_match <- b_match[keep_ref, , drop = FALSE]
        if (length(p_match) > 0L) {
          p_match <- p_match[keep_ref]
        }
      }
      res <- list(boxes = b_match, polygons = p_match, image = im)
      if (isTRUE(plot)) {
        plot(im)
        n_inst <- if (task == "segment" && length(res$polygons) > 0L) length(res$polygons) else nrow(res$boxes)
        if (isTRUE(rainbow) && n_inst > 0L) {
          c_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
          c_poly <- grDevices::adjustcolor(c_border, alpha.f = 0.40)
        } else {
          c_border <- rep("#00FFCC", max(1L, n_inst))
          c_poly <- rep("#00FFCC55", max(1L, n_inst))
        }
        if (task == "segment" && length(res$polygons) > 0L) {
          for (i in seq_along(res$polygons)) {
            poly <- res$polygons[[i]]
            if (!is.null(poly) && length(poly) >= 6L) {
              px <- poly[c(TRUE, FALSE)] * w
              py <- poly[c(FALSE, TRUE)] * h
              graphics::polygon(px, py, col = c_poly[i], border = c_border[i], lwd = 1.5)
            }
          }
        } else if (nrow(res$boxes) > 0L) {
          df_b <- res$boxes
          for (r in seq_len(nrow(df_b))) {
            graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = c_border[r], lwd = 2)
            graphics::text(df_b$xmin[r], df_b$ymin[r] - 5, labels = sprintf("%s %.2f", df_b$class_name[r], df_b$conf[r]),
                           col = c_border[r], cex = 0.8, pos = 4)
          }
        }
      }
      return(invisible(res))
    }
  }

  df_b <- data.frame()
  cand_polys <- list()

  # 2. Run PerSAM visual prototype transfer if SAM method
  if (identical(model$method, "sam") && !is.null(model$prototypes)) {
    num_arr <- image_data(im, type = "numeric")
    persam_res <- tryCatch({
      .run_persam(
        mat = num_arr,
        exemplar_points = NULL,
        precomputed_prototypes = model$prototypes,
        sim_threshold = conf_t,
        min_dist = if (!is.null(model$median_radius)) min(12.0, max(6.0, model$median_radius * 0.25)) else 10.0,
        max_objects = if (!is.null(max_objects)) max_objects else model$max_objects,
        feat_res = 256L,
        engine = if (!is.null(engine)) engine else if (!is.null(model$engine)) model$engine else "cpu",
        device_id = if (!is.null(device_id) && device_id >= 0) device_id else if (!is.null(model$device_id)) model$device_id else -1L,
        fill_hull = TRUE,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(persam_res) && nrow(persam_res$boxes) > 0L) {
      raw_boxes <- persam_res$boxes
      raw_conts <- persam_res$contours

      valid_idx <- seq_len(nrow(raw_boxes))
      if (!is.null(model$min_area) && model$min_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- which(b_areas >= model$min_area)
      }
      if (!is.null(model$max_area) && model$max_area > 0) {
        b_areas <- (raw_boxes$xmax - raw_boxes$xmin) * (raw_boxes$ymax - raw_boxes$ymin)
        valid_idx <- intersect(valid_idx, which(b_areas <= model$max_area))
      }

      if (length(valid_idx) > 0L) {
        raw_boxes <- raw_boxes[valid_idx, , drop = FALSE]
        raw_conts <- raw_conts[valid_idx]

        df_b <- data.frame(
          xmin = as.numeric(raw_boxes$xmin),
          ymin = as.numeric(raw_boxes$ymin),
          xmax = as.numeric(raw_boxes$xmax),
          ymax = as.numeric(raw_boxes$ymax),
          conf = round(as.numeric(raw_boxes$score), 3),
          class_id = if (!is.null(model$class_id)) as.integer(model$class_id) else 0L,
          class_name = model$label[1],
          stringsAsFactors = FALSE
        )

        for (i in seq_along(raw_conts)) {
          cnt <- raw_conts[[i]]
          if (is.matrix(cnt) && nrow(cnt) >= 3L) {
            gx <- cnt[, 1] / w
            gy <- cnt[, 2] / h
            if (length(gx) > 32L) {
              idx <- round(seq(1, length(gx), length.out = 32))
              gx <- gx[idx]; gy <- gy[idx]
            }
            cand_polys[[i]] <- as.vector(rbind(gx, gy))
          } else {
            cand_polys[[i]] <- NULL
          }
        }
      }
    }
  }

  # 3. Fallback to adaptive morphology if PerSAM did not find instances or method != "sam"
  if (nrow(df_b) == 0L) {
    norm_arr <- image_data(im, type = "normalized")
    dr <- (norm_arr[, , 1] - model$target_rgb[1]) / model$rgb_sd[1]
    dg <- (norm_arr[, , 2] - model$target_rgb[2]) / model$rgb_sd[2]
    db <- (norm_arr[, , 3] - model$target_rgb[3]) / model$rgb_sd[3]
    dist_sq <- (dr^2 + dg^2 + db^2) / 3.0
    sim_raw <- exp(-0.5 * dist_sq)

    blur_sigma <- max(2, as.integer(round(model$median_radius * 0.12)))
    sim_im <- as_image(matrix(sim_raw, w, h), colormode = "Grayscale")
    sim_blur <- image_blur(sim_im, sigma = blur_sigma)
    sim_mat <- image_data(sim_blur, type = "normalized")

    adapt_t <- if (conf_t > 0.4) 0.25 else conf_t
    bin_mask <- sim_mat >= adapt_t
    lbl <- bwlabel_cpp(bin_mask)
    counts <- table(lbl[lbl > 0])
    valid <- as.integer(names(counts[counts >= model$min_area & counts <= model$max_area]))

    cand_boxes <- list()
    cand_confs <- numeric()
    conts_all <- if (length(valid) > 0L) extract_contours_cpp(lbl) else list()
    cand_p <- list()

    for (obj in valid) {
      coords <- which(lbl == obj, arr.ind = TRUE)
      bx1 <- min(coords[, 1]); bx2 <- max(coords[, 1])
      by1 <- min(coords[, 2]); by2 <- max(coords[, 2])
      bw <- bx2 - bx1 + 1L; bh <- by2 - by1 + 1L
      ar <- bw / pmax(1L, bh)
      if (ar < 0.38 || ar > 2.60) next

      obj_conf <- mean(sim_mat[coords])
      if (obj_conf < adapt_t) next

      poly_norm <- NULL
      if (length(conts_all) >= obj && is.matrix(conts_all[[obj]]) && nrow(conts_all[[obj]]) >= 3L) {
        poly_pts <- conts_all[[obj]]
        gx <- poly_pts[, 1] / w; gy <- poly_pts[, 2] / h
        if (length(gx) > 32L) {
          idx <- round(seq(1, length(gx), length.out = 32))
          gx <- gx[idx]; gy <- gy[idx]
        }
        poly_norm <- as.vector(rbind(gx, gy))
      }

      cand_boxes[[length(cand_boxes) + 1L]] <- c(bx1, by1, bx2, by2)
      cand_confs <- c(cand_confs, obj_conf)
      cand_p[[length(cand_p) + 1L]] <- poly_norm
    }

    if (length(cand_boxes) > 0L) {
      box_mat <- do.call(rbind, cand_boxes)
      if (nrow(box_mat) > 1L) {
        keep_idx <- .nms_boxes(box_mat, cand_confs, iou_thresh = 0.35)
        box_mat <- box_mat[keep_idx, , drop = FALSE]
        cand_confs <- cand_confs[keep_idx]
        cand_p <- cand_p[keep_idx]
      }
      df_b <- data.frame(
        xmin = as.numeric(box_mat[, 1]),
        ymin = as.numeric(box_mat[, 2]),
        xmax = as.numeric(box_mat[, 3]),
        ymax = as.numeric(box_mat[, 4]),
        conf = round(cand_confs, 3),
        class_id = if (!is.null(model$class_id)) as.integer(model$class_id) else 0L,
        class_name = model$label[1],
        stringsAsFactors = FALSE
      )
      cand_polys <- cand_p
    }
  }

  if (nrow(df_b) == 0L) {
    df_b <- data.frame(xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
                       conf = numeric(0), class_id = integer(0), class_name = character(0), stringsAsFactors = FALSE)
  }

  res <- list(boxes = df_b, polygons = cand_polys, image = im)

  if (isTRUE(plot)) {
    plot(im)
    n_inst <- if (task == "segment" && length(cand_polys) > 0L) length(cand_polys) else nrow(df_b)
    if (isTRUE(rainbow) && n_inst > 0L) {
      c_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
      c_poly <- grDevices::adjustcolor(c_border, alpha.f = 0.40)
    } else {
      c_border <- rep("#00FFCC", max(1L, n_inst))
      c_poly <- rep("#00FFCC55", max(1L, n_inst))
    }
    if (task == "segment" && length(cand_polys) > 0L) {
      for (i in seq_along(cand_polys)) {
        poly <- cand_polys[[i]]
        if (!is.null(poly) && length(poly) >= 6L) {
          px <- poly[c(TRUE, FALSE)] * w
          py <- poly[c(FALSE, TRUE)] * h
          graphics::polygon(px, py, col = c_poly[i], border = c_border[i], lwd = 1.5)
        }
      }
    } else if (nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = c_border[r], lwd = 2)
        graphics::text(df_b$xmin[r], df_b$ymin[r] - 5, labels = sprintf("%s %.2f", df_b$class_name[r], df_b$conf[r]),
                       col = c_border[r], cex = 0.8, pos = 4)
      }
    }
  }

  invisible(res)
}

# Vectorized Non-Maximum Suppression (NMS) helper
.nms_boxes <- function(boxes, scores, iou_thresh = 0.35) {
  if (nrow(boxes) <= 1L) return(seq_len(nrow(boxes)))
  order_idx <- order(scores, decreasing = TRUE)
  keep <- integer()

  while (length(order_idx) > 0L) {
    curr <- order_idx[1]
    keep <- c(keep, curr)
    if (length(order_idx) == 1L) break

    rest <- order_idx[-1]
    xx1 <- pmax(boxes[curr, 1], boxes[rest, 1])
    yy1 <- pmax(boxes[curr, 2], boxes[rest, 2])
    xx2 <- pmin(boxes[curr, 3], boxes[rest, 3])
    yy2 <- pmin(boxes[curr, 4], boxes[rest, 4])

    w_inter <- pmax(0, xx2 - xx1 + 1)
    h_inter <- pmax(0, yy2 - yy1 + 1)
    inter <- w_inter * h_inter

    area_curr <- (boxes[curr, 3] - boxes[curr, 1] + 1) * (boxes[curr, 4] - boxes[curr, 2] + 1)
    area_rest <- (boxes[rest, 3] - boxes[rest, 1] + 1) * (boxes[rest, 4] - boxes[rest, 2] + 1)
    iou <- inter / (area_curr + area_rest - inter)

    order_idx <- rest[iou <= iou_thresh]
  }
  return(keep)
}

#' @title Auto-Label an Entire Directory into a Production YOLO Dataset
#' @name yolo_dataset_autolabel
#' @description
#' `yolo_dataset_autolabel()` batch-processes an entire folder of unannotated images using
#' a trained `yolo_fewshot_model` (from [yolo_fewshot_fit()]). It predicts bounding boxes
#' or instance segmentation polygons, writes standard YOLO `.txt` annotation files, splits
#' the data into `train`, `val`, and `test`, and generates a `data.yaml` descriptor.
#'
#' It also creates a `review/` visual inspection folder with rendered boxes/polygons so the user
#' can quickly audit the automated labels.
#'
#' @param input_dir Directory containing raw unannotated images.
#' @param model A `yolo_fewshot_model` object fitted with [yolo_fewshot_fit()].
#' @param output_dir Target directory to write the YOLO dataset.
#' @param task Target task: `"segment"` (polygons) or `"detect"` (bounding boxes). Defaults to `model$task`.
#' @param conf_threshold Minimum confidence / similarity threshold for accepting predicted instances (e.g. 0.50). If `NULL` (default), uses `model$conf_threshold`.
#' @param split_ratio Numeric vector of length 3 for `c(train, val, test)` split. Defaults to `c(0.70, 0.20, 0.10)`.
#' @param engine Execution engine override: `"cpu"` or `"gpu"`. If `NULL`, uses `model$engine`.
#' @param device_id GPU device ID if `engine = "gpu"`. Defaults to `-1`.
#' @param review Logical. If `TRUE` (default), writes annotated review images to `review/` subfolder.
#' @param rainbow Logical. If `TRUE`, review previews are rendered with distinct rainbow colors per instance. Defaults to `FALSE`.
#' @param seed Random seed for reproducible splitting. Defaults to `42`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#'
#' @return A `yolo_dataset` S3 object pointing to the created dataset.
#' @export
yolo_dataset_autolabel <- function(input_dir,
                                   model,
                                   output_dir,
                                   task = if (!is.null(model$task)) model$task else c("segment", "detect"),
                                   conf_threshold = NULL,
                                   max_objects = NULL,
                                   split_ratio = c(0.70, 0.20, 0.10),
                                   engine = "gpu",
                                   device_id = -1,
                                   review = TRUE,
                                   rainbow = TRUE,
                                   seed = 42,
                                   verbose = TRUE) {
  task <- match.arg(task, c("segment", "detect"))
  if (is.list(model) && !inherits(model, "yolo_fewshot_model")) {
    if (all(vapply(model, function(m) inherits(m, "yolo_fewshot_model"), logical(1)))) {
      model <- do.call(c, model)
    }
  }
  if (!inherits(model, "yolo_fewshot_model")) {
    cli::cli_abort("{.arg model} must be a {.cls yolo_fewshot_model} from {.fn yolo_fewshot_fit}.")
  }

  classes <- if (!is.null(model$classes)) {
    as.character(model$classes)
  } else if (!is.null(model$label)) {
    as.character(model$label)
  } else {
    "object"
  }

  if (!dir.exists(input_dir)) {
    cli::cli_abort("Directory {.path {input_dir}} does not exist.")
  }

  fls <- list.files(input_dir, pattern = "\\.(jpg|jpeg|png|bmp|webp)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fls) == 0L) {
    cli::cli_abort("No valid images found in {.path {input_dir}}.")
  }

  if (!is.null(seed)) set.seed(seed)

  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  splits <- c("train", "val", "test")
  for (sp in splits) {
    dir.create(file.path(output_dir, "images", sp), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "labels", sp), recursive = TRUE, showWarnings = FALSE)
  }

  review_dir <- file.path(output_dir, "review")
  if (isTRUE(review)) dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  # Stratify images
  n_total <- length(fls)
  shuffled <- sample(fls)
  n_tr <- max(1L, round(n_total * split_ratio[1]))
  n_val <- max(1L, round(n_total * split_ratio[2]))

  split_assign <- character(n_total)
  split_assign[seq_len(n_tr)] <- "train"
  split_assign[(n_tr + 1L):min(n_total, n_tr + n_val)] <- "val"
  if (n_total > (n_tr + n_val)) {
    split_assign[(n_tr + n_val + 1L):n_total] <- "test"
  } else {
    split_assign[split_assign == ""] <- "val"
  }

  used_eng <- if (!is.null(engine)) engine else if (!is.null(model$engine)) model$engine else "cpu"
  conf_t <- if (!is.null(conf_threshold)) as.numeric(conf_threshold) else model$conf_threshold

  cli_pb <- NULL
  if (isTRUE(verbose)) {
    cli::cli_h2("Automated YOLO Dataset Labeling (pliman)")
    cls_disp <- paste(classes, collapse = ", ")
    cli::cli_alert_info("Class{?es}: {.val {cls_disp}} | Task: {.val {task}} | Conf: {.val {conf_t}} | Images: {.val {n_total}} | Engine: {.val {used_eng}}")
    cli_pb <- cli::cli_progress_bar(
      format = "{cli::pb_spin} Auto-labeling [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
      total = n_total,
      type = "iterator",
      clear = FALSE
    )
  }

  # Ensure unique destination stems across images (prevents label overwrite if same stem with different extensions)
  stems_raw <- tools::file_path_sans_ext(basename(fls))
  exts_raw  <- tools::file_ext(basename(fls))

  unique_stems <- character(length(fls))
  names(unique_stems) <- fls
  unique_fns <- character(length(fls))
  names(unique_fns) <- fls

  if (any(duplicated(stems_raw))) {
    dup_stems <- unique(stems_raw[duplicated(stems_raw)])
    for (idx in seq_along(fls)) {
      f <- fls[idx]
      s <- stems_raw[idx]
      e <- exts_raw[idx]
      if (s %in% dup_stems) {
        stem_clean <- paste0(s, "_", tolower(e))
        unique_stems[f] <- stem_clean
        unique_fns[f]   <- paste0(stem_clean, ".", e)
      } else {
        unique_stems[f] <- s
        unique_fns[f]   <- basename(f)
      }
    }
  } else {
    for (idx in seq_along(fls)) {
      unique_stems[fls[idx]] <- stems_raw[idx]
      unique_fns[fls[idx]]   <- basename(fls[idx])
    }
  }

  total_objects <- 0L

  for (i in seq_along(shuffled)) {
    im_path <- shuffled[i]
    sp <- split_assign[i]
    pred <- try(yolo_fewshot_predict(model, im_path, task = task, conf_threshold = conf_t, max_objects = max_objects, engine = engine, device_id = device_id, plot = FALSE), silent = TRUE)
    if (inherits(pred, "try-error")) {
      if (isTRUE(verbose) && !is.null(cli_pb)) cli::cli_progress_update(id = cli_pb)
      next
    }

    df_b <- pred$boxes
    polys <- pred$polygons

    # Filter strictly by conf_t
    if (!is.null(conf_t) && nrow(df_b) > 0L && !is.null(df_b$conf)) {
      pass_idx <- which(df_b$conf >= conf_t)
      df_b <- df_b[pass_idx, , drop = FALSE]
      if (length(polys) > 0L) {
        polys <- polys[pass_idx]
      }
    }

    im_obj <- pred$image
    dims <- dim(im_obj)
    w <- dims[1]; h <- dims[2]

    # Destination paths with unique stems
    fn_dest <- unique_fns[[im_path]]
    stem_dest <- unique_stems[[im_path]]
    dest_img <- file.path(output_dir, "images", sp, fn_dest)
    dest_lbl <- file.path(output_dir, "labels", sp, paste0(stem_dest, ".txt"))

    file.copy(im_path, dest_img, overwrite = TRUE)

    # Build YOLO label lines
    lines <- character()
    if (nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        cid <- if ("class_id" %in% names(df_b)) as.integer(df_b$class_id[r]) else 0L
        if (task == "detect") {
          xc <- ((df_b$xmin[r] + df_b$xmax[r]) / 2) / w
          yc <- ((df_b$ymin[r] + df_b$ymax[r]) / 2) / h
          bw <- (df_b$xmax[r] - df_b$xmin[r]) / w
          bh <- (df_b$ymax[r] - df_b$ymin[r]) / h
          lines <- c(lines, sprintf("%d %.6f %.6f %.6f %.6f", cid, xc, yc, bw, bh))
        } else {
          poly <- polys[[r]]
          if (!is.null(poly) && length(poly) >= 6L) {
            lines <- c(lines, paste(c(as.character(cid), sprintf("%.6f", poly)), collapse = " "))
          } else {
            # Bounding box fallback as 4-corner polygon
            xc1 <- df_b$xmin[r] / w; yc1 <- df_b$ymin[r] / h
            xc2 <- df_b$xmax[r] / w; yc2 <- df_b$ymax[r] / h
            lines <- c(lines, sprintf("%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f",
                                      cid, xc1, yc1, xc2, yc1, xc2, yc2, xc1, yc2))
          }
        }
      }
      total_objects <- total_objects + nrow(df_b)
    }

    writeLines(lines, dest_lbl)

    # Optional review image rendering (wrapped in tryCatch to prevent file locks or device errors from crashing autolabel)
    if (isTRUE(review) && i <= 30L) {
      tryCatch({
        rev_out <- file.path(review_dir, paste0("review_", fn_dest))
        out_w <- 1200L
        out_h <- max(100L, as.integer(round(1200 * h / w)))
        dev_opened <- FALSE
        tryCatch({
          grDevices::png(rev_out, width = out_w, height = out_h, res = 120)
          dev_opened <- TRUE
        }, error = function(e) {
          # Try alternative filename if primary is locked by image viewer
          alt_out <- file.path(review_dir, paste0("review_", i, "_", fn_dest))
          tryCatch({
            grDevices::png(alt_out, width = out_w, height = out_h, res = 120)
            dev_opened <<- TRUE
          }, error = function(e2) NULL)
        })

        if (dev_opened) {
          plot(im_obj)
          n_b <- nrow(df_b)
          if (isTRUE(rainbow) && n_b > 0L) {
            rev_border <- grDevices::rainbow(n_b, s = 0.85, v = 0.95)
            rev_poly <- grDevices::adjustcolor(rev_border, alpha.f = 0.40)
          } else {
            rev_border <- rep("#00FFCC", max(1L, n_b))
            rev_poly <- rep("#00FFCC55", max(1L, n_b))
          }
          if (task == "segment" && length(polys) > 0L) {
            for (r in seq_len(n_b)) {
              poly <- polys[[r]]
              if (!is.null(poly) && length(poly) >= 6L) {
                px <- poly[c(TRUE, FALSE)] * w
                py <- poly[c(FALSE, TRUE)] * h
                graphics::polygon(px, py, col = rev_poly[r], border = rev_border[r], lwd = 1.5)
              } else {
                graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = rev_border[r], lwd = 2)
              }
              cur_lbl <- if ("class_name" %in% names(df_b) && nzchar(df_b$class_name[r])) df_b$class_name[r] else classes[1]
              lbl_text <- if (!is.null(df_b$conf) && !is.na(df_b$conf[r])) {
                sprintf("%s (%.2f)", cur_lbl, df_b$conf[r])
              } else {
                cur_lbl
              }
              graphics::text(df_b$xmin[r], df_b$ymin[r] - 4,
                             labels = lbl_text,
                             col = rev_border[r], cex = 0.8, pos = 4)
            }
          } else if (n_b > 0L) {
            for (r in seq_len(n_b)) {
              graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = rev_border[r], lwd = 2)
              cur_lbl <- if ("class_name" %in% names(df_b) && nzchar(df_b$class_name[r])) df_b$class_name[r] else classes[1]
              lbl_text <- if (!is.null(df_b$conf) && !is.na(df_b$conf[r])) {
                sprintf("%s (%.2f)", cur_lbl, df_b$conf[r])
              } else {
                cur_lbl
              }
              graphics::text(df_b$xmin[r], df_b$ymin[r] - 4,
                             labels = lbl_text,
                             col = rev_border[r], cex = 0.8, pos = 4)
            }
          }
          graphics::title(sprintf("%s: %d %s", fn_dest, n_b, if (task == "segment") "segmented" else "detected"),
                          col.main = if (isTRUE(rainbow)) "#222222" else "#00FFCC",
                          cex.main = 1.0)
          grDevices::dev.off()
        }
      }, error = function(e) {
        # Ensure graphics device is closed if an unexpected error occurs during plotting
        if (grDevices::dev.cur() > 1L) {
          tryCatch(grDevices::dev.off(), error = function(e) NULL)
        }
      })
    }

    if (isTRUE(verbose) && !is.null(cli_pb)) {
      cli::cli_progress_update(id = cli_pb)
    }
  }

  if (isTRUE(verbose) && !is.null(cli_pb)) {
    cli::cli_progress_done(id = cli_pb)
  }

  # Write data.yaml
  yaml_content <- c(
    paste0("path: ", output_dir),
    "train: images/train",
    "val: images/val",
    "test: images/test",
    "",
    "names:",
    paste0("  ", seq_along(classes) - 1L, ": ", classes)
  )
  writeLines(yaml_content, file.path(output_dir, "data.yaml"))

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Dataset auto-labeled successfully!")
    cli::cli_alert_info("Total instances generated: {.val {total_objects}}")
    if (isTRUE(review)) {
      cli::cli_alert_info("Review previews written to: {.path {review_dir}}")
    }
  }

  .yolo_wrap_dataset(output_dir, task = task, classes = classes)
}

#' @title Plot a Few-Shot Exemplar Model
#' @name plot.yolo_fewshot_model
#' @description
#' Plots reference image(s) with bounding boxes, segmentation polygons, and clicked exemplar points.
#' When the model contains multiple reference images, they are arranged automatically in a multi-panel grid.
#' Setting \code{rainbow = TRUE} renders each detected instance in distinct vibrant colors.
#'
#' @param x A \code{yolo_fewshot_model} object.
#' @param task Target task: \code{"segment"} or \code{"detect"}. Defaults to \code{x$task}.
#' @param rainbow Logical. If \code{TRUE}, colors each detected instance with distinct vibrant hues. Defaults to \code{FALSE}.
#' @param show_boxes Logical. Show bounding boxes. Defaults to \code{TRUE} for detection tasks.
#' @param show_polygons Logical. Show segmentation polygons. Defaults to \code{TRUE} for segmentation tasks.
#' @param show_points Logical. Show clicked exemplar points. Defaults to \code{TRUE}.
#' @param col_border Border color when \code{rainbow = FALSE}. Defaults to \code{"#00FFCC"}.
#' @param col_poly Fill color for polygons when \code{rainbow = FALSE}. Defaults to \code{"#00FFCC55"}.
#' @param col_point Point color for exemplars. Defaults to \code{"#FFFF00"}.
#' @param lwd Line width for borders. Defaults to \code{2}.
#' @param mfrow Numeric vector of length 2 for grid layout (e.g. \code{c(2, 2)}). If \code{NULL}, calculated automatically.
#' @param which Indices or names of reference images to plot. If \code{NULL} (default), plots all reference images.
#' @param ... Additional arguments passed to graphics methods.
#'
#' @return Invisible \code{x}.
#' @export
plot.yolo_fewshot_model <- function(x,
                                    task = x$task,
                                    rainbow = TRUE,
                                    show_boxes = NULL,
                                    show_polygons = NULL,
                                    show_points = TRUE,
                                    col_border = "#00FFCC",
                                    col_poly = "#00FFCC55",
                                    col_point = "#FFFF00",
                                    lwd = 2,
                                    mfrow = NULL,
                                    which = NULL,
                                    ...) {
  if (inherits(x, "yolo_fewshot_multimodel") && !is.null(x$models)) {
    for (m in x$models) {
      plot(m, task = task, rainbow = rainbow, show_boxes = show_boxes,
           show_polygons = show_polygons, show_points = show_points,
           lwd = lwd, ...)
    }
    return(invisible(x))
  }
  imgs <- if (!is.null(x$reference_images) && length(x$reference_images) > 0L) {
    x$reference_images
  } else if (!is.null(x$reference_image)) {
    list(x$reference_image)
  } else {
    cli::cli_abort("No reference images found in the model to plot.")
  }

  if (!is.null(which)) {
    imgs <- imgs[which]
  }

  n_imgs <- length(imgs)
  if (n_imgs == 0L) {
    cli::cli_abort("Selected index in {.arg which} does not match any reference image.")
  }

  task <- if (!is.null(task)) match.arg(task, c("segment", "detect")) else if (!is.null(x$task)) x$task else "segment"
  if (is.null(show_polygons)) show_polygons <- (task == "segment")
  if (is.null(show_boxes)) show_boxes <- (task == "detect" || !show_polygons)

  if (n_imgs > 1L) {
    if (is.null(mfrow)) {
      nc <- ceiling(sqrt(n_imgs))
      nr <- ceiling(n_imgs / nc)
      mfrow <- c(nr, nc)
    }
    op <- graphics::par(mfrow = mfrow, mar = c(1, 1, 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  img_keys <- names(imgs)

  for (k in seq_along(imgs)) {
    im_k <- imgs[[k]]
    k_id <- if (!is.null(img_keys) && length(img_keys) >= k) img_keys[k] else NULL

    df_b <- if (!is.null(k_id) && !is.null(x$reference_boxes_map[[k_id]])) {
      x$reference_boxes_map[[k_id]]
    } else if (k == 1L && !is.null(x$detected_boxes)) {
      x$detected_boxes
    } else {
      data.frame()
    }

    polys <- if (!is.null(k_id) && !is.null(x$reference_polys_map[[k_id]])) {
      x$reference_polys_map[[k_id]]
    } else if (k == 1L && !is.null(x$detected_polygons)) {
      x$detected_polygons
    } else {
      list()
    }

    pts_k <- if (!is.null(k_id) && !is.null(x$exemplar_points_map[[k_id]])) {
      x$exemplar_points_map[[k_id]]
    } else if (k == 1L && !is.null(x$exemplar_points)) {
      x$exemplar_points
    } else {
      NULL
    }

    plot(im_k)
    w <- dim(im_k)[1]
    h <- dim(im_k)[2]

    n_inst <- if (task == "segment" && length(polys) > 0L) {
      length(polys)
    } else {
      nrow(df_b)
    }

    if (isTRUE(rainbow) && n_inst > 0L) {
      cols_border <- grDevices::rainbow(n_inst, s = 0.85, v = 0.95)
      cols_poly <- grDevices::adjustcolor(cols_border, alpha.f = 0.40)
    } else {
      cols_border <- rep(col_border, max(1L, n_inst))
      cols_poly <- rep(col_poly, max(1L, n_inst))
    }

    # Draw segmentation polygons if available
    if (isTRUE(show_polygons) && length(polys) > 0L) {
      for (i in seq_along(polys)) {
        poly <- polys[[i]]
        if (!is.null(poly) && length(poly) >= 6L) {
          px <- poly[c(TRUE, FALSE)] * w
          py <- poly[c(FALSE, TRUE)] * h
          cb <- cols_border[((i - 1L) %% length(cols_border)) + 1L]
          cp <- cols_poly[((i - 1L) %% length(cols_poly)) + 1L]
          graphics::polygon(px, py, col = cp, border = cb, lwd = max(1, lwd - 1))
        }
      }
    }

    # Draw exemplar points
    if (isTRUE(show_points) && !is.null(pts_k) && length(pts_k) > 0L) {
      if (is.matrix(pts_k) || is.data.frame(pts_k)) {
        graphics::points(pts_k[, 1], pts_k[, 2], col = col_point, pch = 19, cex = 1.3)
      }
    }

    # Draw detected bounding boxes
    if (isTRUE(show_boxes) && nrow(df_b) > 0L) {
      for (r in seq_len(nrow(df_b))) {
        cb <- cols_border[((r - 1L) %% length(cols_border)) + 1L]
        graphics::rect(df_b$xmin[r], df_b$ymin[r], df_b$xmax[r], df_b$ymax[r], border = cb, lwd = lwd)
        if (!is.null(df_b$conf[r])) {
          graphics::text(df_b$xmin[r], df_b$ymin[r] - 4,
                         labels = sprintf("%s (%.2f)", x$label, df_b$conf[r]),
                         col = cb, cex = 0.8, pos = 4)
        }
      }
    }

    n_det_k <- if (task == "segment" && length(polys) > 0L) {
      sum(!vapply(polys, is.null, logical(1)))
    } else {
      nrow(df_b)
    }

    if (n_imgs > 1L) {
      fn_label <- if (!is.null(k_id) && file.exists(k_id)) basename(k_id) else sprintf("Reference %d", k)
      graphics::title(sprintf("%s: %d %s", fn_label, n_det_k, if (task == "segment") "segmented" else "detected"),
                      col.main = if (isTRUE(rainbow)) "#222222" else col_border,
                      cex.main = 1.0)
    } else {
      graphics::title(sprintf("Few-Shot Model: '%s' (%d %s) [%s]",
                              x$label, n_det_k,
                              if (task == "segment") "segmented" else "detected",
                              task),
                      col.main = col_border)
    }
  }

  invisible(x)
}

#' @export
print.yolo_fewshot_model <- function(x, ...) {
  cli::cli_h2("pliman Few-Shot Exemplar Model")
  cls_txt <- if (!is.null(x$classes)) paste(x$classes, collapse = ", ") else paste(x$label, collapse = ", ")
  cli::cli_text("Target Class{?es}  : {.val {cls_txt}}")
  n_refs <- if (!is.null(x$reference_images)) length(x$reference_images) else 1L
  cli::cli_text("Reference Images : {.val {n_refs}}")
  n_det <- if (!is.null(x$reference_boxes_map)) sum(vapply(x$reference_boxes_map, nrow, integer(1))) else if (!is.null(x$detected_boxes)) nrow(x$detected_boxes) else 0L
  cli::cli_text("Detected on Refs : {.val {n_det}} instances")
  cli::cli_text("Confidence Thresh: {.val {x$conf_threshold}}")
  cli::cli_text("Task             : {.val {x$task}}")
  cli::cli_text("Engine           : {.val {x$engine}}")
  invisible(x)
}
