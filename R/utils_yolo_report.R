# ------------------------------------------------------------------------------
# YOLO DATASET INTERACTIVE AUDIT & PROFILING REPORT
# High-fidelity, self-contained HTML report with modern UI, annotated galleries,
# provenance matrices, spatial heatmaps, training advisor, and automated health checks.
# ------------------------------------------------------------------------------

#' Generate an Interactive HTML Report for a YOLO Dataset
#'
#' @description
#' `yolo_dataset_report()` conducts an end-to-end audit and visual profiling of a
#' YOLO dataset (detection, segmentation, OBB, or classification) and generates a
#' self-contained, interactive HTML report with a modern UI.
#'
#' The report features:
#' * **Executive KPI Dashboard**: Metrics for images, objects, density,
#'   provenance, and class balance, alongside an **Overall Dataset Health Score Gauge** (0-100%).
#' * **Dataset Integrity & Anomaly Inspector**: Interactive clickable anomaly pills
#'   (Duplicate Pairs, Extreme AR, Micro-Boxes, Out-of-Bounds, Zero-Area) linking
#'   directly to an interactive table and modal viewer with image thumbnails and annotations.
#' * **Dataset Provenance & Composition**: Matrix and visual breakdown of original,
#'   augmented (`_aug`), synthetic background (`synth_bg_`), and tiled (`_tile_`) images.
#' * **Splits & Class Distribution**: Interactive tables and visual allocation
#'   across `train`, `val`, and `test`.
#' * **Spatial Object Location Heatmap**: 12x12 grid of normalized object center coordinates
#'   with top and left marginal distribution bars to identify camera/sensor bias.
#' * **Bounding Box Geometry & Aspect Ratio**: 10x10 Width vs. Height distribution grid
#'   with 1:1 diagonal reference line and aspect ratio categories (Portrait, Square-ish, Landscape).
#' * **Interactive Annotated Image Gallery**: Filterable thumbnail gallery showcasing
#'   fairly sampled images across splits (`train`, `val`, `test`) and provenance categories
#'   with pixel-accurate SVG overlays, live search-by-name, split filtering,
#'   anomalies tab, object count sorting, and full-screen modal Lightbox.
#' * **Multi-Class Co-occurrence Matrix**: Pairwise co-occurrence heatmap when multiple classes exist.
#' * **AI & Training Hyperparameter Advisor**: Tailor-made recommendations for `imgsz`,
#'   `epochs`, `batch`, `max_det` (preventing truncated detections in dense scenes),
#'   and augmentations, with a copy-to-clipboard R `yolo_train()` code snippet.
#' * **Automated Health Checklist**: 10-point automated diagnostic evaluating train/val
#'   ratios, class balance, background sample presence, resolution consistency, and
#'   annotation integrity.
#' * **Export Tools**: Dark/Light theme toggle, Print/PDF, and direct JSON audit export.
#'
#' @param data A directory path to a YOLO dataset, or an existing `yolo_dataset_desc`
#'   object returned by [yolo_dataset_describe()].
#' @param output_file Character string specifying the output HTML file path.
#'   Default is `"yolo_dataset_report.html"` in the current working directory.
#' @param title Character string for the report title. Default is
#'   `"YOLO Dataset Audit & Profiling Report"`.
#' @param max_samples_per_category Integer. Maximum number of candidate images
#'   to sample and embed per category (Original, Augmented, Synthetic, Tiled,
#'   Negative/Background). Default is `8L`.
#' @param max_thumbnail_size Integer. Maximum pixel dimension (width or height)
#'   for embedded thumbnail images. Large images are resized before embedding to keep
#'   the HTML report compact, fast, and portable. Default is `640L`.
#' @param lwd Numeric. Line width (in pixels) for bounding boxes and segmentation
#'   polygon contours in image thumbnail overlays and the Lightbox modal (default: `1.5`).
#'   Can also be adjusted interactively in real time via the report's slider.
#' @param task Character specifying the YOLO task: `"auto"` (default), `"detect"`,
#'   `"segment"`, `"obb"`, or `"classify"`. When `"auto"`, the task is inferred
#'   from the dataset directory structure or `data.yaml`.
#' @param check_blur Logical. If `TRUE` (default), runs blur detection using Laplacian variance.
#' @param blur_threshold Numeric threshold for blur detection. Defaults to `100.0`.
#' @param max_blur_samples Integer. Maximum number of original images to sample and audit for blurriness (default: `500L`).
#'   For large datasets with thousands of high-resolution images, sampling avoids reading every single image from disk,
#'   dramatically speeding up report generation while retaining a reliable quality assessment. Set to `Inf` or `NULL`
#'   to audit 100% of images.
#' @param browse Logical. If `TRUE` (default when interactive), automatically opens
#'   the generated HTML report in the default web browser or RStudio Viewer.
#' @param theme Character specifying the initial color theme: `"dark"` (default,
#'   modern AI dashboard look) or `"light"` (clean paper document style).
#'   Users can also toggle themes interactively inside the HTML report.
#' @param verbose Logical. If `TRUE` (default), prints progress and summary
#'   messages to the console.
#' @param ... Additional options passed to internal report renderers (e.g., `box_lwd`).
#'
#' @return A character string containing the normalized absolute path to the
#'   generated HTML file (invisibly).
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @export
#' @rdname yolo_dataset_report
#'
#' @examples
#' \dontrun{
#' # Generate an interactive HTML report from a dataset directory
#' yolo_dataset_report("path/to/dataset", output_file = "capsulas_report.html")
#'
#' # Adjust bounding box line width to 2.0px
#' yolo_dataset_report("path/to/dataset", lwd = 2.0)
#'
#' # Generate a report from an existing describe object
#' desc <- yolo_dataset_describe("path/to/dataset")
#' yolo_dataset_report(desc, output_file = "dataset_report.html")
#' }
yolo_dataset_report <- function(data,
                                output_file = "yolo_dataset_report.html",
                                title = "YOLO Dataset Audit & Profiling Report",
                                max_samples_per_category = 8L,
                                max_thumbnail_size = 640L,
                                lwd = 1.5,
                                task = c("auto", "detect", "segment", "obb", "classify"),
                                check_blur = TRUE,
                                blur_threshold = 80,
                                max_blur_samples = 500L,
                                browse = interactive(),
                                theme = c("dark", "light"),
                                verbose = TRUE,
                                ...) {
  dots <- list(...)
  if (!is.null(dots$box_lwd)) lwd <- dots$box_lwd
  if (!is.null(dots$linewidth)) lwd <- dots$linewidth
  if (!is.null(dots$line_width)) lwd <- dots$line_width
  lwd <- as.numeric(lwd)[1]
  if (is.na(lwd) || lwd <= 0) lwd <- 1.5
  theme <- match.arg(theme)
  task <- match.arg(task)
  blur_threshold <- as.numeric(blur_threshold)[1]
  if (is.na(blur_threshold) || blur_threshold <= 0) blur_threshold <- 80


  start_time <- Sys.time()

  # 1. Obtain or reuse yolo_dataset_desc object
  if (inherits(data, "yolo_dataset_desc")) {
    desc <- data
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Reusing provided {.cls yolo_dataset_desc} statistical profile.")
    }
  } else if (is.character(data) && length(data) == 1L) {
    if (!dir.exists(data)) {
      cli::cli_abort("Directory {.path {data}} does not exist.")
    }
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Profiling dataset at {.path {data}} via {.fun yolo_dataset_describe}...")
    }
    desc <- yolo_dataset_describe(
      data = data,
      sample_dims = 30,
      verbose = FALSE,
      plot = FALSE
    )
    if (task != "auto") {
      desc$task <- task
    }
  } else {
    cli::cli_abort("{.arg data} must be either a directory path or a {.cls yolo_dataset_desc} object.")
  }

  # Build anomaly lookup by image stem
  anomaly_lookup <- list()
  if (!is.null(desc$anomalies$records) && length(desc$anomalies$records) > 0L) {
    for (ar in desc$anomalies$records) {
      st <- ar$stem
      anomaly_lookup[[st]] <- c(anomaly_lookup[[st]], list(ar))
    }
  }

  # 2. Analyze blur & focus quality across original images
  blur_audit <- list(
    checked = 0L,
    blurry = 0L,
    pct_blurry = 0.0,
    threshold = blur_threshold,
    scores = list()
  )
  blurry_records <- list()

  if (isTRUE(check_blur)) {
    all_recs <- do.call(c, unname(desc$split_records))
    if (!is.null(all_recs) && length(all_recs) > 0L) {
      orig_idx <- which(vapply(all_recs, function(r) identical(r$ptype, "Original"), logical(1)))
      target_recs <- if (length(orig_idx) > 0L) all_recs[orig_idx] else all_recs

      img_paths <- vapply(target_recs, function(r) r$img, character(1))
      img_stems <- vapply(target_recs, function(r) r$stem %||% tools::file_path_sans_ext(basename(r$img)), character(1))

      valid_mask <- file.exists(img_paths)
      if (any(valid_mask)) {
        valid_paths <- img_paths[valid_mask]
        valid_stems <- img_stems[valid_mask]
        valid_recs <- target_recs[valid_mask]
        n_valid <- length(valid_paths)

        # Apply sampling limit if dataset has many images to prevent bottleneck on disk I/O & decompression
        if (!is.null(max_blur_samples) && is.finite(max_blur_samples) && n_valid > max_blur_samples) {
          sel_idx <- round(seq(1L, n_valid, length.out = as.integer(max_blur_samples)))
          sel_paths <- valid_paths[sel_idx]
          sel_stems <- valid_stems[sel_idx]
          sel_recs <- valid_recs[sel_idx]
        } else {
          sel_paths <- valid_paths
          sel_stems <- valid_stems
          sel_recs <- valid_recs
        }

        if (isTRUE(verbose)) {
          sample_note <- if (length(sel_paths) < n_valid) sprintf(" (sampled %d from %d)", length(sel_paths), n_valid) else ""
          cli::cli_alert_info("Analyzing image sharpness & focus quality across {length(sel_paths)} original images{sample_note}...")
        }

        scores <- images_blur_cpp(sel_paths, max_size = max_thumbnail_size)
        is_blurry_vec <- scores < blur_threshold

        blur_audit$checked <- length(scores)
        blur_audit$blurry <- sum(is_blurry_vec)
        blur_audit$pct_blurry <- round((blur_audit$blurry / max(1L, blur_audit$checked)) * 100, 1)

        for (k in seq_along(sel_stems)) {
          blur_audit$scores[[sel_stems[k]]] <- scores[k]
        }
        blurry_records <- sel_recs[is_blurry_vec]
      }
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Sampling candidate images across categories, splits (train, val, test), anomalies, and blurry samples...")
  }

  # 3. Sample gallery candidates across categories, splits, anomalies, and blurry samples
  sample_records <- .yolo_sample_gallery_candidates(
    split_records = desc$split_records,
    max_per_cat = max_samples_per_category,
    anomaly_records = desc$anomalies$records,
    blurry_records = blurry_records
  )

  # 4. Process candidate images (thumbnails, base64 data URIs, SVG annotations, sharpness score)
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Generating thumbnails and SVG annotation overlays for {length(sample_records)} sample images...")
  }

  processed_cards <- list()
  for (r in sample_records) {
    stem_name <- r$stem %||% tools::file_path_sans_ext(basename(r$img))
    card_info <- .yolo_process_sample_card(
      record = r,
      classes = desc$classes,
      task = desc$task,
      max_thumbnail_size = max_thumbnail_size,
      anomalies = anomaly_lookup[[stem_name]],
      lwd = lwd,
      blur_score = blur_audit$scores[[stem_name]],
      blur_threshold = blur_threshold
    )
    if (!is.null(card_info)) {
      processed_cards[[length(processed_cards) + 1L]] <- card_info
    }
  }

  # 5. Generate automated health checklist and score
  health_audit <- .yolo_build_health_audit(desc, blur_audit = blur_audit)

  # 6. Build full HTML document
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Assembling interactive HTML report...")
  }

  html_content <- .yolo_generate_report_html(
    desc = desc,
    cards = processed_cards,
    health = health_audit,
    title = title,
    theme = theme,
    lwd = lwd,
    blur_threshold = blur_threshold
  )

  # Ensure destination directory exists
  out_dir <- dirname(output_file)
  if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  writeLines(html_content, con = output_file, useBytes = TRUE)
  abs_path <- normalizePath(output_file, mustWork = FALSE)
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Report generated in {elapsed}s: {.path {abs_path}}")
  }

  if (isTRUE(browse)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() && rstudioapi::hasFun("viewer")) {
      rstudioapi::viewer(abs_path)
    } else {
      utils::browseURL(abs_path)
    }
  }

  invisible(abs_path)
}

#' @rdname yolo_dataset_report
#' @export
yolo_report <- yolo_dataset_report


# ==============================================================================
# AUGMENTATION DEFINITIONS & RECOGNITION
# ==============================================================================
.yolo_aug_defs <- list(
  blur     = list(label = "Blur",         icon = "\U0001f32b\ufe0f", color = "#ec4899"),
  noise    = list(label = "Noise",        icon = "\u26a1", color = "#f59e0b"),
  rot      = list(label = "Rotate",       icon = "\U0001f504", color = "#6366f1"),
  persp    = list(label = "Perspective",  icon = "\U0001f4d0", color = "#8b5cf6"),
  fliph    = list(label = "Flip H",       icon = "\u2194\ufe0f", color = "#10b981"),
  flipv    = list(label = "Flip V",       icon = "\u2195\ufe0f", color = "#14b8a6"),
  bright   = list(label = "Brightness",   icon = "\u2600\ufe0f", color = "#eab308"),
  contrast = list(label = "Contrast",     icon = "\U0001f313", color = "#f97316"),
  cutout   = list(label = "Cutout",       icon = "\u2702\ufe0f", color = "#ef4444"),
  mosaic   = list(label = "Mosaic",       icon = "\U0001f9e9", color = "#06b6d4"),
  synth_bg = list(label = "Synthetic BG", icon = "\U0001f5bc\ufe0f", color = "#f59e0b"),
  tile     = list(label = "Tiled",        icon = "\U0001fa9f", color = "#10b981")
)

.yolo_extract_aug_tags <- function(stem) {
  if (is.null(stem) || !nzchar(stem)) return(character(0))
  tags <- character(0)
  s <- tolower(stem)
  if (grepl("(^|_)blur($|_|[0-9])", s)) tags <- c(tags, "blur")
  if (grepl("(^|_)noise($|_|[0-9])", s)) tags <- c(tags, "noise")
  if (grepl("(^|_)rot(_|-)?[0-9]+", s)) tags <- c(tags, "rot")
  if (grepl("(^|_)persp($|_)", s)) tags <- c(tags, "persp")
  if (grepl("(^|_)fliph($|_)", s)) tags <- c(tags, "fliph")
  if (grepl("(^|_)flipv($|_)", s)) tags <- c(tags, "flipv")
  if (grepl("(^|_)bright($|_)", s)) tags <- c(tags, "bright")
  if (grepl("(^|_)contrast($|_)", s)) tags <- c(tags, "contrast")
  if (grepl("(^|_)cutout($|_)", s)) tags <- c(tags, "cutout")
  if (grepl("(^|_)mosaic($|_)", s)) tags <- c(tags, "mosaic")
  if (grepl("(^|_)synth(_bg)?($|_)", s) || grepl("synthetic", s)) tags <- c(tags, "synth_bg")
  if (grepl("_tile_", s)) tags <- c(tags, "tile")
  unique(tags)
}


# ==============================================================================
# INTERNAL HELPER: SAMPLE CANDIDATE IMAGES PER PROVENANCE CATEGORY & SPLIT
# ==============================================================================
.yolo_sample_gallery_candidates <- function(split_records, max_per_cat = 8L, anomaly_records = NULL, blurry_records = NULL) {
  if (is.null(split_records) || length(split_records) == 0L) {
    return(list())
  }

  all_recs <- do.call(c, unname(split_records))
  if (length(all_recs) == 0L) return(list())

  cats <- c("Original", "Augmented", "Synthetic_BG", "Tiled", "Background")
  sampled <- list()
  sampled_stems <- character(0)

  get_cat <- function(r) {
    if (isTRUE(r$is_background) || identical(r$n_objs, 0L)) {
      return("Background")
    }
    pt <- r$ptype %||% "Original"
    if (pt == "Augmented") return("Augmented")
    if (pt == "Synthetic_BG") return("Synthetic_BG")
    if (pt == "Tiled") return("Tiled")
    return("Original")
  }

  splits_available <- names(split_records)
  n_splits <- length(splits_available)
  per_split_target <- max(2L, ceiling(max_per_cat / max(1L, n_splits)))

  # 1. Fair sampling across each split and each category
  for (cat_name in cats) {
    cat_chosen_stems <- character(0)
    for (sp in splits_available) {
      sp_recs <- split_records[[sp]]
      if (length(sp_recs) == 0L) next
      sp_cats <- vapply(sp_recs, get_cat, character(1))
      match_idx <- which(sp_cats == cat_name)
      if (length(match_idx) > 0L) {
        pick_n <- min(length(match_idx), per_split_target)
        for (i in seq_len(pick_n)) {
          rec <- sp_recs[[match_idx[i]]]
          s <- rec$stem %||% tools::file_path_sans_ext(basename(rec$img))
          if (!s %in% sampled_stems) {
            sampled[[length(sampled) + 1L]] <- rec
            sampled_stems <- c(sampled_stems, s)
            cat_chosen_stems <- c(cat_chosen_stems, s)
          }
        }
      }
    }
    # If category has fewer than max_per_cat, backfill from remaining available
    if (length(cat_chosen_stems) < max_per_cat) {
      all_cat_recs <- all_recs[vapply(all_recs, get_cat, character(1)) == cat_name]
      for (rec in all_cat_recs) {
        if (length(cat_chosen_stems) >= max_per_cat) break
        s <- rec$stem %||% tools::file_path_sans_ext(basename(rec$img))
        if (!s %in% sampled_stems) {
          sampled[[length(sampled) + 1L]] <- rec
          sampled_stems <- c(sampled_stems, s)
          cat_chosen_stems <- c(cat_chosen_stems, s)
        }
      }
    }
  }

  # 2. Guarantee that images containing anomalies are sampled so they can be inspected in Lightbox
  if (!is.null(anomaly_records) && length(anomaly_records) > 0L) {
    anom_stems <- unique(vapply(anomaly_records, function(a) a$stem, character(1)))
    anom_stems_to_add <- utils::head(setdiff(anom_stems, sampled_stems), 40L)
    if (length(anom_stems_to_add) > 0L) {
      all_stems <- vapply(all_recs, function(r) r$stem %||% tools::file_path_sans_ext(basename(r$img)), character(1))
      for (astem in anom_stems_to_add) {
        idx <- which(all_stems == astem)[1]
        if (!is.na(idx)) {
          rec <- all_recs[[idx]]
          sampled[[length(sampled) + 1L]] <- rec
          sampled_stems <- c(sampled_stems, astem)
        }
      }
    }
  }

  # 3. Guarantee that images detected as blurry are sampled for immediate inspection
  if (!is.null(blurry_records) && length(blurry_records) > 0L) {
    blur_stems <- unique(vapply(blurry_records, function(b) b$stem %||% tools::file_path_sans_ext(basename(b$img)), character(1)))
    blur_stems_to_add <- utils::head(setdiff(blur_stems, sampled_stems), 30L)
    if (length(blur_stems_to_add) > 0L) {
      all_stems <- vapply(all_recs, function(r) r$stem %||% tools::file_path_sans_ext(basename(r$img)), character(1))
      for (bstem in blur_stems_to_add) {
        idx <- which(all_stems == bstem)[1]
        if (!is.na(idx)) {
          rec <- all_recs[[idx]]
          sampled[[length(sampled) + 1L]] <- rec
          sampled_stems <- c(sampled_stems, bstem)
        }
      }
    }
  }

  # 4. Guarantee representative sampling across ALL augmentation transformation types
  aug_types <- names(.yolo_aug_defs)
  for (at in aug_types) {
    matching_recs <- list()
    for (r in all_recs) {
      st <- r$stem %||% tools::file_path_sans_ext(basename(r$img))
      if (at %in% .yolo_extract_aug_tags(st)) {
        matching_recs[[length(matching_recs) + 1L]] <- r
      }
    }
    if (length(matching_recs) > 0L) {
      target_n <- min(5L, length(matching_recs))
      sampled_for_type <- 0L
      for (r in matching_recs) {
        if (sampled_for_type >= target_n) break
        st <- r$stem %||% tools::file_path_sans_ext(basename(r$img))
        if (!st %in% sampled_stems) {
          sampled[[length(sampled) + 1L]] <- r
          sampled_stems <- c(sampled_stems, st)
          sampled_for_type <- sampled_for_type + 1L
        }
      }
    }
  }

  sampled
}


# ==============================================================================
# INTERNAL HELPER: PROCESS A SAMPLE IMAGE CARD (THUMBNAIL & SVG ANNOTATIONS)
# ==============================================================================
.yolo_process_sample_card <- function(record, classes, task, max_thumbnail_size = 640L, anomalies = NULL, lwd = 1.5, blur_score = NULL, blur_threshold = 80) {
  img_path <- record$img
  if (!file.exists(img_path)) return(NULL)

  im <- try(image_import(img_path), silent = TRUE)
  if (inherits(im, "try-error") || !is_image(im)) {
    return(NULL)
  }

  d <- dim(im)
  h_orig <- d[1]
  w_orig <- d[2]

  # Compute sharpness/blur score if not provided
  b_score <- if (!is.null(blur_score)) blur_score else image_blur_cpp(im, max_size = max_thumbnail_size)
  b_score <- round(as.numeric(b_score)[1], 1)
  is_blurry <- isTRUE(b_score < blur_threshold)

  if (max(h_orig, w_orig) > max_thumbnail_size) {
    scale_fac <- max_thumbnail_size / max(h_orig, w_orig)
    im <- image_resize(im, rel_size = scale_fac * 100)
  }

  tf <- tempfile(fileext = ".jpg")
  on.exit(if (file.exists(tf)) unlink(tf), add = TRUE)
  .write_image_file(im, tf, quality = 80)

  b64_str <- base64_encode(tf, file = TRUE)
  data_uri <- paste0("data:image/jpeg;base64,", b64_str)

  stem <- record$stem %||% tools::file_path_sans_ext(basename(img_path))
  split_name <- record$split %||% "all"
  is_bg <- isTRUE(record$is_background) || identical(record$n_objs, 0L)
  ptype <- record$ptype %||% "Original"

  cat_slug <- if (is_bg) "background" else tolower(gsub("[^a-zA-Z0-9]", "", ptype))
  cat_label <- if (is_bg) "NEGATIVE (BG)" else toupper(ptype)

  palette <- c(
    "#6366f1", "#10b981", "#f59e0b", "#ef4444", "#8b5cf6",
    "#ec4899", "#06b6d4", "#f97316", "#14b8a6", "#84cc16",
    "#3b82f6", "#d946ef", "#eab308", "#64748b"
  )

  objects_list <- list()
  svg_elements <- character(0)

  if (!is_bg && !is.null(record$lbl) && file.exists(record$lbl)) {
    lines <- tryCatch(readLines(record$lbl, warn = FALSE), error = function(e) character(0))
    lines <- lines[nzchar(trimws(lines))]

    for (ln_idx in seq_along(lines)) {
      parts <- suppressWarnings(as.numeric(strsplit(trimws(lines[ln_idx]), "\\s+")[[1]]))
      if (length(parts) < 5L) next

      cid <- as.integer(parts[1])
      cname <- if (cid >= 0L && (cid + 1L) <= length(classes)) classes[cid + 1L] else paste0("class_", cid)
      color <- palette[(cid %% length(palette)) + 1L]

      if (length(parts) == 5L) {
        xc <- parts[2] * 100
        yc <- parts[3] * 100
        bw <- parts[4] * 100
        bh <- parts[5] * 100
        xmin <- max(0, xc - bw / 2)
        ymin <- max(0, yc - bh / 2)

        objects_list[[length(objects_list) + 1L]] <- list(
          type = "BBox",
          class = cname,
          color = color,
          area_pct = round(parts[4] * parts[5] * 100, 2),
          coords = sprintf("xc=%.2f, yc=%.2f, w=%.2f, h=%.2f", parts[2], parts[3], parts[4], parts[5])
        )

        svg_elements <- c(svg_elements, sprintf(
          '<rect class="yolo-geo yolo-bbox" x="%.2f" y="%.2f" width="%.2f" height="%.2f" stroke="%s" stroke-width="%.2f" fill="%s" fill-opacity="0.18" rx="0.5" data-class="%s" vector-effect="non-scaling-stroke" />',
          xmin, ymin, bw, bh, color, lwd, color, cname
        ))

        lbl_w <- max(4.0, nchar(cname) * 1.5 + 1.2)
        lbl_y <- max(0.5, ymin - 3.2)
        svg_elements <- c(svg_elements, sprintf(
          '<g class="yolo-label" transform="translate(%.2f, %.2f)"><rect x="0" y="0" width="%.2f" height="3.2" rx="0.6" fill="%s" /><text x="0.8" y="2.3" font-size="2.1" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-weight="700" fill="#ffffff">%s</text></g>',
          xmin, lbl_y, lbl_w, color, cname
        ))
      } else {
        coords <- parts[-1]
        xs <- coords[seq(1, length(coords), 2)] * 100
        ys <- coords[seq(2, length(coords), 2)] * 100
        pts_str <- paste(sprintf("%.2f,%.2f", xs, ys), collapse = " ")

        bw <- max(xs) - min(xs)
        bh <- max(ys) - min(ys)
        xmin <- min(xs)
        ymin <- min(ys)

        objects_list[[length(objects_list) + 1L]] <- list(
          type = "Polygon",
          class = cname,
          color = color,
          area_pct = round((bw * bh) / 100, 2),
          coords = sprintf("%d vertices", length(xs))
        )

        svg_elements <- c(svg_elements, sprintf(
          '<polygon class="yolo-geo yolo-poly" points="%s" stroke="%s" stroke-width="%.2f" fill="%s" fill-opacity="0.22" data-class="%s" vector-effect="non-scaling-stroke" />',
          pts_str, color, lwd, color, cname
        ))

        lbl_w <- max(4.0, nchar(cname) * 1.5 + 1.2)
        lbl_y <- max(0.5, ymin - 3.2)
        svg_elements <- c(svg_elements, sprintf(
          '<g class="yolo-label" transform="translate(%.2f, %.2f)"><rect x="0" y="0" width="%.2f" height="3.2" rx="0.6" fill="%s" /><text x="0.8" y="2.3" font-size="2.1" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-weight="700" fill="#ffffff">%s</text></g>',
          xmin, lbl_y, lbl_w, color, cname
        ))
      }
    }
  } else if (!is.null(record$class)) {
    cname <- record$class
    color <- palette[1]
    objects_list[[1]] <- list(
      type = "Class",
      class = cname,
      color = color,
      area_pct = 100.0,
      coords = "Full image"
    )
  }

  svg_overlay <- if (length(svg_elements) > 0L) {
    sprintf(
      '<svg class="yolo-svg-overlay" viewBox="0 0 100 100" preserveAspectRatio="none">%s</svg>',
      paste(svg_elements, collapse = "")
    )
  } else ""

  has_anom <- !is.null(anomalies) && length(anomalies) > 0L
  anom_labels <- if (has_anom) unique(vapply(anomalies, function(a) a$type_label %||% a$type, character(1))) else character(0)

  list(
    stem = stem,
    split = split_name,
    category_slug = cat_slug,
    category_label = cat_label,
    provenance = ptype,
    is_background = is_bg,
    n_objects = length(objects_list),
    width = w_orig,
    height = h_orig,
    data_uri = data_uri,
    svg_overlay = svg_overlay,
    objects = objects_list,
    has_anomaly = has_anom,
    anomalies = anomalies %||% list(),
    anomaly_summary = paste(anom_labels, collapse = ", "),
    blur_score = b_score,
    is_blurry = is_blurry,
    aug_tags = .yolo_extract_aug_tags(stem)
  )
}


# ==============================================================================
# INTERNAL HELPER: AUTOMATED DATASET HEALTH AUDIT & SCORING
# ==============================================================================
.yolo_build_health_audit <- function(desc, blur_audit = NULL) {
  items <- list()
  anom <- desc$anomalies %||% list()

  # 1. Dataset Scale & Volume Check
  tot_imgs <- desc$total_images %||% 0L
  if (tot_imgs >= 500L) {
    items[[length(items) + 1L]] <- list(
      title = "Dataset Sample Volume",
      status = "pass",
      badge = "OPTIMAL",
      metric = sprintf("%s images", format(tot_imgs, big.mark = ",")),
      detail = "Sufficient volume for robust transfer learning and generalization across diverse scenarios."
    )
  } else if (tot_imgs >= 100L) {
    items[[length(items) + 1L]] <- list(
      title = "Dataset Sample Volume",
      status = "info",
      badge = "MODERATE",
      metric = sprintf("%s images", format(tot_imgs, big.mark = ",")),
      detail = "Sufficient for transfer learning on pre-trained YOLO backbones; offline augmentation can further enrich diversity."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Dataset Sample Volume",
      status = "warning",
      badge = "LOW VOLUME",
      metric = sprintf("%d images", tot_imgs),
      detail = "Small dataset (< 100 images). High risk of overfitting. Recommend expanding via yolo_dataset_augment() and synthetic backgrounds."
    )
  }

  # 2. Split Ratio Allocation
  sp_df <- desc$splits_summary
  if (!is.null(sp_df)) {
    sub_sp <- sp_df[sp_df$Split != "TOTAL", , drop = FALSE]
    train_row <- sub_sp[sub_sp$Split == "train", , drop = FALSE]
    val_row <- sub_sp[sub_sp$Split == "val", , drop = FALSE]
    test_row <- sub_sp[sub_sp$Split == "test", , drop = FALSE]

    if (nrow(train_row) > 0L && nrow(val_row) > 0L) {
      tr_pct <- train_row$Pct_Images[1]
      val_pct <- val_row$Pct_Images[1]
      has_test <- nrow(test_row) > 0L

      if (tr_pct >= 60 && tr_pct <= 85 && val_pct >= 10 && val_pct <= 30) {
        items[[length(items) + 1L]] <- list(
          title = "Train / Validation Allocation",
          status = "pass",
          badge = "BALANCED",
          metric = sprintf("Train: %.1f%% | Val: %.1f%%%s", tr_pct, val_pct, if (has_test) sprintf(" | Test: %.1f%%", test_row$Pct_Images[1]) else ""),
          detail = "Standard partition matches machine learning best practices (~70-80% train, ~10-20% validation)."
        )
      } else {
        items[[length(items) + 1L]] <- list(
          title = "Train / Validation Allocation",
          status = "warning",
          badge = "UNEVEN SPLIT",
          metric = sprintf("Train: %.1f%% | Val: %.1f%%", tr_pct, val_pct),
          detail = "Train/Val ratio departs from standard recommended distribution. Review partition to avoid unrepresentative validation."
        )
      }

      if (has_test) {
        items[[length(items) + 1L]] <- list(
          title = "Unbiased Test Set",
          status = "pass",
          badge = "PRESENT",
          metric = sprintf("%s test images (%.1f%%)", format(test_row$Images[1], big.mark = ","), test_row$Pct_Images[1]),
          detail = "Dedicated test set ensures objective evaluation isolated from hyperparameter tuning."
        )
      } else {
        items[[length(items) + 1L]] <- list(
          title = "Unbiased Test Set",
          status = "info",
          badge = "OPTIONAL",
          metric = "0 test images",
          detail = "No test partition found. While train/val suffice for development, a test set is recommended for final benchmark validation."
        )
      }
    } else {
      items[[length(items) + 1L]] <- list(
        title = "Train / Validation Allocation",
        status = "warning",
        badge = "NO SPLITS",
        metric = "Single partition",
        detail = "Dataset is not partitioned into train/val/test splits. Partition with yolo_dataset_split() prior to model training."
      )
    }
  }

  # 3. Class Imbalance
  imb <- desc$imbalance_ratio %||% 1.0
  n_classes <- length(desc$classes)
  if (n_classes <= 1L) {
    items[[length(items) + 1L]] <- list(
      title = "Class Distribution",
      status = "pass",
      badge = "SINGLE CLASS",
      metric = sprintf("1 class (%s)", desc$classes[1]),
      detail = "Single-class object detection / segmentation dataset. Class imbalance is inherently not applicable."
    )
  } else if (imb <= 2.0) {
    items[[length(items) + 1L]] <- list(
      title = "Class Distribution",
      status = "pass",
      badge = "WELL BALANCED",
      metric = sprintf("%.2fx Max/Min Ratio", imb),
      detail = "Even distribution across classes ensures equal gradient attention and avoids minority class neglect."
    )
  } else if (imb <= 4.0) {
    items[[length(items) + 1L]] <- list(
      title = "Class Distribution",
      status = "info",
      badge = "MODERATE SKEW",
      metric = sprintf("%.2fx Max/Min Ratio", imb),
      detail = "Moderate class frequency difference. Pre-trained YOLO models generally handle ratios up to 3-4x without loss reweighting."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Class Distribution",
      status = "warning",
      badge = "HIGH IMBALANCE",
      metric = sprintf("%.2fx Max/Min Ratio", imb),
      detail = "High class skew detected. Minority classes risk under-detection. Consider targeted augmentation on scarce classes."
    )
  }

  # 4. Background / Negative Sample Ratio
  if (desc$task != "classify") {
    bg_pct <- desc$pct_background %||% 0.0
    bg_cnt <- desc$total_background %||% 0L
    if (bg_pct >= 1.5 && bg_pct <= 12.0) {
      items[[length(items) + 1L]] <- list(
        title = "Background / False-Positive Defense",
        status = "pass",
        badge = "OPTIMAL",
        metric = sprintf("%s negative images (%.1f%%)", format(bg_cnt, big.mark = ","), bg_pct),
        detail = "Ideal proportion (recommended 2-10% negative images without objects) to drastically reduce background false positives."
      )
    } else if (bg_cnt == 0L) {
      items[[length(items) + 1L]] <- list(
        title = "Background / False-Positive Defense",
        status = "warning",
        badge = "ZERO BACKGROUND",
        metric = "0 negative images (0%)",
        detail = "Zero background images. Models trained without negative samples often trigger false positive hallucinations. Consider adding 2-5% negative images via yolo_dataset_add_background()."
      )
    } else {
      items[[length(items) + 1L]] <- list(
        title = "Background / False-Positive Defense",
        status = "info",
        badge = "ELEVATED",
        metric = sprintf("%s negative images (%.1f%%)", format(bg_cnt, big.mark = ","), bg_pct),
        detail = "Background images exceed 12%. Check that unannotated positive samples are not mistakenly treated as empty background."
      )
    }
  }

  # 5. Provenance & Offline Augmentation
  tot_aug <- desc$total_augmented %||% 0L
  tot_synth <- desc$total_synthetic %||% 0L
  tot_tile <- desc$total_tiled %||% 0L
  if ((tot_aug + tot_synth + tot_tile) > 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Dataset Augmentation Diversity",
      status = "pass",
      badge = "ENRICHED",
      metric = sprintf("%s offline augmented / synthetic / tiled images", format(tot_aug + tot_synth + tot_tile, big.mark = ",")),
      detail = "Offline data augmentation expands training variance and improves model invariance to lighting, orientation, and scale."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Dataset Augmentation Diversity",
      status = "info",
      badge = "BASELINE ONLY",
      metric = "0 offline augmented images",
      detail = "Entire dataset consists of original annotations. YOLO performs online augmentations during training, but offline expansions via yolo_dataset_augment() can double sample diversity."
    )
  }

  # 6. Resolution Consistency
  if (length(desc$sample_img_dims) > 0L) {
    dim_strings <- vapply(desc$sample_img_dims, function(d) sprintf("%d x %d", d[1], d[2]), character(1))
    uniq_dims <- unique(dim_strings)
    if (length(uniq_dims) == 1L) {
      items[[length(items) + 1L]] <- list(
        title = "Resolution Uniformity",
        status = "pass",
        badge = "PERFECT UNIFORMITY",
        metric = sprintf("All samples: %s", uniq_dims[1]),
        detail = "Uniform image dimensions accelerate batch tensor loading and avoid padding artifacts."
      )
    } else {
      top_dim <- names(sort(table(dim_strings), decreasing = TRUE))[1]
      items[[length(items) + 1L]] <- list(
        title = "Resolution Uniformity",
        status = "info",
        badge = "MIXED RESOLUTIONS",
        metric = sprintf("Dominant: %s (%d variations)", top_dim, length(uniq_dims)),
        detail = "Mixed image dimensions. Ultralytics YOLO automatically resizes and pads with letterbox during preprocessing."
      )
    }
  }

  # 7. Coordinate Boundary Integrity
  n_oob <- anom$out_of_bounds %||% 0L
  if (n_oob == 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Bounding Box Coordinate Bounds",
      status = "pass",
      badge = "PERFECT BOUNDS",
      metric = "0 out-of-bounds coordinates",
      detail = "All normalized coordinates reside strictly within the valid [0.0, 1.0] image domain."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Bounding Box Coordinate Bounds",
      status = "warning",
      badge = "COORDINATE OVERFLOW",
      metric = sprintf("%s boxes out-of-bounds", format(n_oob, big.mark = ",")),
      detail = "Detected bounding boxes exceeding [0, 1] normalized coordinates. Click anomaly badge to view affected images."
    )
  }

  # 8. Box Scale & Geometry Sanity
  n_zero <- anom$zero_area %||% 0L
  n_micro <- anom$micro_boxes %||% 0L
  if (n_zero == 0L && n_micro == 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Annotation Geometry Sanity",
      status = "pass",
      badge = "VALID GEOMETRIES",
      metric = "0 zero-area or micro-boxes",
      detail = "No collapsed degenerate boxes or sub-pixel artifacts (< 0.05% image area) detected."
    )
  } else if (n_zero > 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Annotation Geometry Sanity",
      status = "warning",
      badge = "DEGENERATE BOXES",
      metric = sprintf("%s zero-area boxes detected", format(n_zero, big.mark = ",")),
      detail = "Zero-area boxes cause division-by-zero errors in CIoU/GIoU loss computation. Click badge to view images."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Annotation Geometry Sanity",
      status = "info",
      badge = "MICRO-OBJECTS",
      metric = sprintf("%s micro-boxes (<0.05%% area)", format(n_micro, big.mark = ",")),
      detail = "Micro-boxes detected. Click badge to inspect affected images and verify recognizable target features."
    )
  }

  # 9. Duplicate & Overlapping Annotations
  n_dup <- anom$duplicate_pairs %||% 0L
  if (n_dup == 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Annotation Redundancy",
      status = "pass",
      badge = "NO DUPLICATES",
      metric = "0 duplicate bounding boxes",
      detail = "No high-overlap redundant duplicate annotation pairs (IoU > 0.85) identified."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Annotation Redundancy",
      status = "warning",
      badge = "DUPLICATES DETECTED",
      metric = sprintf("%s duplicate pairs (IoU > 0.85)", format(n_dup, big.mark = ",")),
      detail = "Duplicate overlapping boxes detected. Click badge to inspect images and remove redundant labels."
    )
  }

  # 10. Aspect Ratio Extremes
  n_ext_ar <- anom$extreme_aspect_ratios %||% 0L
  if (n_ext_ar == 0L) {
    items[[length(items) + 1L]] <- list(
      title = "Aspect Ratio Plausibility",
      status = "pass",
      badge = "NATURAL SHAPES",
      metric = "0 extreme aspect ratios",
      detail = "All bounding boxes have standard proportions (0.125 <= AR <= 8.0), matching standard anchor distributions."
    )
  } else {
    items[[length(items) + 1L]] <- list(
      title = "Aspect Ratio Plausibility",
      status = "info",
      badge = "EXTREME SLENDERNESS",
      metric = sprintf("%s boxes with AR > 8 or < 0.125", format(n_ext_ar, big.mark = ",")),
      detail = "Ultra-slender or needle-like boxes detected. Click badge to inspect whether annotations are edge slicing artifacts."
    )
  }

  # 11. Image Focus & Blur Quality Check
  if (!is.null(blur_audit) && blur_audit$checked > 0L) {
    n_blur <- blur_audit$blurry
    pct_blur <- blur_audit$pct_blurry
    b_thresh <- blur_audit$threshold

    if (n_blur == 0L) {
      items[[length(items) + 1L]] <- list(
        title = "Image Sharpness & Focus Quality",
        status = "pass",
        badge = "CRISP & SHARP",
        metric = "0% blurry images (100% sharp)",
        detail = sprintf("All %d inspected original images exhibit sharp focus and crisp edge contrast (Laplacian variance >= %.0f).", blur_audit$checked, b_thresh)
      )
    } else {
      items[[length(items) + 1L]] <- list(
        title = "Image Sharpness & Focus Quality",
        status = "warning",
        badge = "BLUR DETECTED",
        metric = sprintf("%s blurry original images (%.1f%%)", format(n_blur, big.mark = ","), pct_blur),
        detail = sprintf("%.1f%% of inspected original images appear out-of-focus or motion-blurred (Laplacian variance < %.0f). Blurry training samples weaken CNN/backbone feature extraction and increase false negatives.", pct_blur, b_thresh),
        action_html = sprintf('<button type="button" class="btn-inspect-blur" onclick="inspectBlurryImages()" style="margin-top:0.65rem; display:inline-flex; align-items:center; gap:0.4rem; padding:0.35rem 0.75rem; font-size:0.8rem; font-weight:600; color:#ef4444; background:rgba(239, 68, 68, 0.12); border:1px solid rgba(239, 68, 68, 0.35); border-radius:6px; cursor:pointer; transition:all 0.2s ease;">\U0001f50d View %d Blurry Images in Gallery</button>', n_blur)
      )
    }
  }

  # Score calculation
  n_pass <- sum(vapply(items, function(x) x$status == "pass", logical(1)))
  n_info <- sum(vapply(items, function(x) x$status == "info", logical(1)))
  n_warn <- sum(vapply(items, function(x) x$status == "warning", logical(1)))
  tot_checks <- length(items)

  pts <- (n_pass * 10) + (n_info * 7) + (n_warn * 3)
  score <- max(10L, min(100L, round(pts / (tot_checks * 10) * 100)))

  grade <- if (score >= 94) "A+" else if (score >= 85) "A" else if (score >= 75) "B" else if (score >= 60) "C" else "D"
  grade_cls <- if (score >= 85) "grade-opt" else if (score >= 70) "grade-good" else if (score >= 55) "grade-med" else "grade-warn"

  list(
    items = items,
    score = score,
    grade = grade,
    grade_class = grade_cls,
    n_pass = n_pass,
    n_info = n_info,
    n_warn = n_warn
  )
}


# ==============================================================================
# INTERNAL HELPER: ANOMALY BADGES STRIP WITH CLICKABLE INSPECTORS
# ==============================================================================
.yolo_build_anomaly_strip_html <- function(desc) {
  anom <- desc$anomalies %||% list()
  n_oob <- anom$out_of_bounds %||% 0L
  n_zero <- anom$zero_area %||% 0L
  n_micro <- anom$micro_boxes %||% 0L
  n_huge <- anom$huge_boxes %||% 0L
  n_dup <- anom$duplicate_pairs %||% 0L
  n_ext_ar <- anom$extreme_aspect_ratios %||% 0L
  n_inv_cid <- anom$invalid_class_ids %||% 0L

  tot_critical <- n_oob + n_zero + n_dup + n_inv_cid
  has_any <- (tot_critical + n_micro + n_huge + n_ext_ar) > 0L

  pills <- character(0)

  if (!has_any) {
    pills <- c(pills, '<span class="anomaly-pill pill-perfect" title="All annotations strictly adhere to YOLO format specifications"><span>\u2713</span> 100% Clean Geometry \u2022 Zero Anomalies Detected</span>')
  } else {
    if (tot_critical == 0L) {
      pills <- c(pills, '<span class="anomaly-pill pill-pass" title="Zero critical bounding box errors found"><span>\u2713</span> Zero Critical Defects</span>')
    }
    if (n_dup > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-danger pill-clickable" onclick="openAnomalyModal(\'duplicate_pair\')" title="Click to view images with duplicate pairs"><span>\U0001f504</span> %s Duplicate Pairs <span class="pill-inspect">View \u2197</span></button>', format(n_dup, big.mark = ",")))
    }
    if (n_micro > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-warn pill-clickable" onclick="openAnomalyModal(\'micro_box\')" title="Click to view images with micro-boxes"><span>\U0001f50d</span> %s Micro-Boxes (&lt;0.05%%) <span class="pill-inspect">View \u2197</span></button>', format(n_micro, big.mark = ",")))
    }
    if (n_ext_ar > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-info pill-clickable" onclick="openAnomalyModal(\'extreme_aspect_ratio\')" title="Click to view images with extreme aspect ratio"><span>\U0001f4cf</span> %s Extreme AR (&gt;8 or &lt;0.125) <span class="pill-inspect">View \u2197</span></button>', format(n_ext_ar, big.mark = ",")))
    }
    if (n_oob > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-danger pill-clickable" onclick="openAnomalyModal(\'out_of_bounds\')" title="Click to view images with out-of-bounds coordinates"><span>\u26a0\ufe0f</span> %s Out-of-Bounds <span class="pill-inspect">View \u2197</span></button>', format(n_oob, big.mark = ",")))
    }
    if (n_zero > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-danger pill-clickable" onclick="openAnomalyModal(\'zero_area\')" title="Click to view images with zero-area boxes"><span>\U0001f6ab</span> %s Zero-Area <span class="pill-inspect">View \u2197</span></button>', format(n_zero, big.mark = ",")))
    }
    if (n_inv_cid > 0L) {
      pills <- c(pills, sprintf('<button type="button" class="anomaly-pill pill-danger pill-clickable" onclick="openAnomalyModal(\'invalid_class_id\')" title="Click to view images with invalid class IDs"><span>\u2753</span> %s Invalid Class IDs <span class="pill-inspect">View \u2197</span></button>', format(n_inv_cid, big.mark = ",")))
    }
  }

  inspect_btn <- if (has_any) {
    '<button type="button" class="btn-action btn-inspect-all" onclick="openAnomalyModal(\'all\')" title="Open Interactive Anomaly Inspector"><span>\U0001f50d</span> Open Anomaly Inspector</button>'
  } else ""

  paste0(
    '<div class="anomaly-strip-container">\n',
    '  <div class="anomaly-strip-title">Integrity & Geometry:</div>\n',
    '  <div class="anomaly-pills">', paste(pills, collapse = "\n"), '</div>\n',
    '  ', inspect_btn, '\n',
    '</div>\n'
  )
}


# ==============================================================================
# INTERNAL HELPER: SPATIAL 12x12 & DIMENSION 10x10 CHARTS
# ==============================================================================
.yolo_build_spatial_charts_html <- function(desc) {
  # 1. 12x12 Spatial Center Heatmap
  sp_mat <- desc$spatial_density
  has_sp <- !is.null(sp_mat) && is.matrix(sp_mat) && nrow(sp_mat) == 12L && ncol(sp_mat) == 12L

  spatial_html <- ""
  if (has_sp) {
    max_val <- max(sp_mat, 1L)
    tot_pts <- sum(sp_mat)

    x_marg <- colSums(sp_mat)
    y_marg <- rowSums(sp_mat)
    max_x_marg <- max(x_marg, 1L)
    max_y_marg <- max(y_marg, 1L)

    top_bars <- character(12)
    for (j in 1:12) {
      pct_h <- round(x_marg[j] / max_x_marg * 100, 1)
      pct_tot <- if (tot_pts > 0) round(x_marg[j] / tot_pts * 100, 1) else 0
      top_bars[j] <- sprintf(
        '<div class="marg-col" title="X Bin %d (%.2f-%.2f): %s objects (%.1f%%)"><div class="marg-col-bar" style="height: %.1f%%;"></div></div>',
        j, (j - 1) / 12, j / 12, format(x_marg[j], big.mark = ","), pct_tot, max(pct_h, 4)
      )
    }

    left_bars <- character(12)
    for (i in 1:12) {
      pct_w <- round(y_marg[i] / max_y_marg * 100, 1)
      pct_tot <- if (tot_pts > 0) round(y_marg[i] / tot_pts * 100, 1) else 0
      left_bars[i] <- sprintf(
        '<div class="marg-row" title="Y Bin %d (%.2f-%.2f): %s objects (%.1f%%)"><div class="marg-row-bar" style="width: %.1f%%;"></div></div>',
        i, (i - 1) / 12, i / 12, format(y_marg[i], big.mark = ","), pct_tot, max(pct_w, 4)
      )
    }

    cells <- character(144)
    idx <- 1L
    for (i in 1:12) {
      for (j in 1:12) {
        val <- sp_mat[i, j]
        ratio <- if (max_val > 0) val / max_val else 0
        pct_tot <- if (tot_pts > 0) round(val / tot_pts * 100, 2) else 0

        opacity <- if (val == 0) 0.04 else min(1, 0.15 + (ratio ^ 0.65) * 0.85)
        bg_color <- if (ratio > 0.70) {
          sprintf("rgba(245, 158, 11, %.2f)", opacity)
        } else if (ratio > 0.35) {
          sprintf("rgba(14, 165, 233, %.2f)", opacity)
        } else {
          sprintf("rgba(99, 102, 241, %.2f)", opacity)
        }

        cells[idx] <- sprintf(
          '<div class="heatmap-cell" style="background-color: %s;" title="Location [%.2f-%.2f, %.2f-%.2f]&#10;Count: %s (%.2f%%)"></div>',
          bg_color,
          (j - 1) / 12, j / 12, (i - 1) / 12, i / 12,
          format(val, big.mark = ","), pct_tot
        )
        idx <- idx + 1L
      }
    }

    spatial_html <- paste0(
      '<div class="content-card">',
      '  <div class="card-head-row">',
      '    <h3 class="card-subtitle">\U0001f4cd Spatial Object Location Heatmap (12\u00d712)</h3>',
      '    <span class="info-badge" title="Normalized object center coordinates (xc, yc)">Normalized Center (xc, yc)</span>',
      '  </div>',
      '  <p class="chart-caption">Visualizes sensor & camera placement bias. Marginal distribution bars on Top (X axis) and Left (Y axis) reveal peripheral cropping, edge clustering, or center bias.</p>',
      '  <div class="spatial-heatmap-wrapper">',
      '    <div class="marg-corner"></div>',
      '    <div class="marg-top-container">', paste(top_bars, collapse = ""), '</div>',
      '    <div class="marg-left-container">', paste(left_bars, collapse = ""), '</div>',
      '    <div class="heatmap-grid-12x12">',
      '      <div class="heatmap-crosshair-v"></div>',
      '      <div class="heatmap-crosshair-h"></div>',
      paste(cells, collapse = ""),
      '    </div>',
      '  </div>',
      '  <div class="heatmap-legend">',
      '    <span>0 objects</span>',
      '    <div class="legend-gradient"></div>',
      '    <span>Peak (', format(max_val, big.mark = ","), ')</span>',
      '  </div>',
      '</div>'
    )
  }

  # 2. 10x10 Bounding Box Dimension Distribution (Width vs. Height)
  dim_mat <- desc$box_dim_density
  has_dim <- !is.null(dim_mat) && is.matrix(dim_mat) && nrow(dim_mat) == 10L && ncol(dim_mat) == 10L

  dim_html <- ""
  if (has_dim) {
    max_d_val <- max(dim_mat, 1L)
    tot_d_pts <- sum(dim_mat)

    d_cells <- character(100)
    d_idx <- 1L
    for (i in 10:1) { # height descending so row 10 (tall) is at top
      for (j in 1:10) { # width ascending (0.0 -> 1.0)
        val <- dim_mat[i, j]
        ratio <- if (max_d_val > 0) val / max_d_val else 0
        pct_tot <- if (tot_d_pts > 0) round(val / tot_d_pts * 100, 2) else 0

        opacity <- if (val == 0) 0.04 else min(1, 0.15 + (ratio ^ 0.65) * 0.85)
        bg_color <- if (ratio > 0.70) {
          sprintf("rgba(236, 72, 153, %.2f)", opacity)
        } else if (ratio > 0.35) {
          sprintf("rgba(168, 85, 247, %.2f)", opacity)
        } else {
          sprintf("rgba(99, 102, 241, %.2f)", opacity)
        }

        d_cells[d_idx] <- sprintf(
          '<div class="heatmap-cell" style="background-color: %s;" title="Width: [%.1f-%.1f], Height: [%.1f-%.1f]&#10;Count: %s (%.2f%%)"></div>',
          bg_color,
          (j - 1) / 10, j / 10, (i - 1) / 10, i / 10,
          format(val, big.mark = ","), pct_tot
        )
        d_idx <- d_idx + 1L
      }
    }

    ar_rows <- character(0)
    if (!is.null(desc$aspect_ratio_summary)) {
      for (k in seq_len(nrow(desc$aspect_ratio_summary))) {
        r <- desc$aspect_ratio_summary[k, ]
        bar_col <- if (grepl("Portrait", r$Category)) "#0ea5e9" else if (grepl("Square", r$Category)) "#10b981" else "#f59e0b"
        ar_rows <- c(ar_rows, sprintf(
          '<tr><td><strong>%s</strong></td><td>%s</td><td><div class="table-bar-cell"><div class="mini-bar-bg"><div class="mini-bar-fill" style="width: %.1f%%; background-color: %s;"></div></div><span>%.1f%%</span></div></td></tr>',
          r$Category, format(r$Objects, big.mark = ","), r$Pct, bar_col, r$Pct
        ))
      }
    }

    dim_html <- paste0(
      '<div class="content-card">',
      '  <div class="card-head-row">',
      '    <h3 class="card-subtitle">\U0001f4d0 Bounding Box Geometry & Aspect Ratio</h3>',
      '    <span class="info-badge" title="Width vs Height normalized dimensions with 1:1 aspect ratio diagonal">Width \u00d7 Height</span>',
      '  </div>',
      '  <p class="chart-caption">Normalized Width (X) vs Height (Y). The dashed diagonal marks 1:1 aspect ratio (square). Points above diagonal are tall/portrait; points below are wide/landscape.</p>',
      '  <div class="dim-heatmap-container">',
      '    <div class="dim-y-axis-label">Height (norm) \u2192</div>',
      '    <div class="dim-heatmap-inner">',
      '      <div class="heatmap-grid-10x10">',
      '        <div class="diagonal-reference-line" title="1:1 Aspect Ratio (Square Objects)"></div>',
      paste(d_cells, collapse = ""),
      '      </div>',
      '      <div class="dim-x-axis-label">Width (norm) \u2192</div>',
      '    </div>',
      '  </div>',
      '  <div style="margin-top: 1.25rem;">',
      '    <h4 style="font-size:0.85rem; font-weight:700; text-transform:uppercase; color:var(--text-muted); margin-bottom:0.5rem;">Aspect Ratio Categories</h4>',
      '    <table class="report-table" style="font-size:0.82rem;">',
      '      <thead><tr><th>Shape Category</th><th>Objects</th><th>Share</th></tr></thead>',
      '      <tbody>', paste(ar_rows, collapse = "\n"), '</tbody>',
      '    </table>',
      '  </div>',
      '</div>'
    )
  }

  paste0(
    '<div class="grid-2col" style="margin-bottom: 2rem;">',
    spatial_html,
    dim_html,
    '</div>'
  )
}


# ==============================================================================
# INTERNAL HELPER: MULTI-CLASS CO-OCCURRENCE MATRIX
# ==============================================================================
.yolo_build_cooccurrence_html <- function(desc) {
  cooc <- desc$class_cooccurrence
  if (is.null(cooc) || !is.matrix(cooc) || length(desc$classes) <= 1L) {
    return("")
  }

  clnames <- desc$classes
  n <- length(clnames)
  max_cooc <- max(cooc, 1L)

  th_headers <- paste(vapply(seq_along(clnames), function(i) {
    sprintf('<th title="%s (ID %d)"><div class="th-class-rot">%s</div></th>', clnames[i], i - 1L, clnames[i])
  }, character(1)), collapse = "")

  rows <- character(0)
  for (i in 1:n) {
    tds <- character(n)
    for (j in 1:n) {
      val <- cooc[i, j]
      is_diag <- (i == j)
      ratio <- if (max_cooc > 0) val / max_cooc else 0
      opacity <- if (val == 0) 0.05 else min(1, 0.15 + (ratio ^ 0.7) * 0.85)

      cell_bg <- if (is_diag) {
        sprintf("rgba(99, 102, 241, %.2f)", opacity)
      } else {
        sprintf("rgba(16, 185, 129, %.2f)", opacity)
      }

      diag_cls <- if (is_diag) ' cooc-diag' else ''
      tds[j] <- sprintf(
        '<td class="cooc-cell%s" style="background-color: %s;" title="%s + %s: %s images">%s</td>',
        diag_cls, cell_bg, clnames[i], clnames[j], format(val, big.mark = ","),
        if (val > 0) format(val, big.mark = ",") else "-"
      )
    }

    rows <- c(rows, sprintf(
      '<tr><td class="cooc-row-header"><strong>%s</strong></td>%s</tr>',
      clnames[i], paste(tds, collapse = "")
    ))
  }

  paste0(
    '<section class="section-block" id="cooccurrence">\n',
    '  <div class="section-header">\n',
    '    <div>\n',
    '      <h2 class="section-title">\U0001f517 Multi-Class Co-occurrence Matrix</h2>\n',
    '      <p class="section-desc">Pairwise frequencies showing how often classes appear together in the same image. Diagonal represents total images containing that class.</p>\n',
    '    </div>\n',
    '  </div>\n',
    '  <div class="content-card">\n',
    '    <div class="table-responsive">\n',
    '      <table class="cooc-table">\n',
    '        <thead><tr><th>Class</th>', th_headers, '</tr></thead>\n',
    '        <tbody>', paste(rows, collapse = "\n"), '</tbody>\n',
    '      </table>\n',
    '    </div>\n',
    '  </div>\n',
    '</section>\n\n'
  )
}


# ==============================================================================
# INTERNAL HELPER: AI & TRAINING HYPERPARAMETER ADVISOR
# ==============================================================================
.yolo_build_training_advisor <- function(desc) {
  tot_imgs <- desc$total_images %||% 0L
  max_dens <- desc$max_density %||% 0L
  task <- desc$task %||% "detect"

  dominant_dim <- c(640, 640)
  if (length(desc$sample_img_dims) > 0L) {
    dim_strings <- vapply(desc$sample_img_dims, function(d) sprintf("%d x %d", d[1], d[2]), character(1))
    top_dim_str <- names(sort(table(dim_strings), decreasing = TRUE))[1]
    parts <- as.integer(strsplit(top_dim_str, " x ")[[1]])
    if (length(parts) == 2L && !any(is.na(parts))) dominant_dim <- parts
  }
  max_native_dim <- max(dominant_dim)

  pct_small <- 0.0
  if (!is.null(desc$scale_summary)) {
    sm_row <- desc$scale_summary[grepl("Small", desc$scale_summary$Scale), , drop = FALSE]
    if (nrow(sm_row) > 0L) pct_small <- sm_row$Pct_Objects[1]
  }

  rec_imgsz <- 640L
  imgsz_rationale <- "Standard 640px input provides optimal inference speed and gradient stability."
  if (max_native_dim >= 1024L && pct_small >= 20.0) {
    rec_imgsz <- 1024L
    imgsz_rationale <- sprintf("Higher resolution 1024px preserves high-frequency details for %.1f%% small objects (<0.5%% area).", pct_small)
  } else if (max_native_dim >= 1280L && pct_small >= 10.0) {
    rec_imgsz <- 1280L
    imgsz_rationale <- sprintf("Native resolution %dx%d with small objects benefits from 1280px tensor scaling.", dominant_dim[1], dominant_dim[2])
  } else if (max_native_dim <= 480L) {
    rec_imgsz <- 480L
    imgsz_rationale <- "Dataset images are natively small; upscaling to 640px adds compute without extra detail."
  }

  rec_epochs <- 100L
  rec_patience <- 20L
  if (tot_imgs < 150L) {
    rec_epochs <- 250L
    rec_patience <- 50L
    epoch_rationale <- "Small sample volume (<150 images) requires extended fine-tuning with higher patience to converge."
  } else if (tot_imgs < 500L) {
    rec_epochs <- 180L
    rec_patience <- 35L
    epoch_rationale <- "Moderate sample volume benefits from ~180 epochs with cosine learning rate scheduling."
  } else if (tot_imgs > 5000L) {
    rec_epochs <- 80L
    rec_patience <- 15L
    epoch_rationale <- "Substantial image volume reaches validation convergence within 80 epochs."
  } else {
    epoch_rationale <- "Standard 100-150 epochs provides balanced convergence for transfer learning."
  }

  rec_batch <- if (rec_imgsz >= 1280L) 8L else if (rec_imgsz >= 1024L) 16L else 32L
  batch_rationale <- sprintf("Recommended batch=%d optimizes GPU memory utilization for %dpx tensor grids.", rec_batch, rec_imgsz)

  rec_max_det <- 300L
  max_det_rationale <- "Default 300 maximum detections per image accommodates standard scenarios."
  if (max_dens >= 200L) {
    rec_max_det <- as.integer(ceiling(max_dens * 1.4 / 50) * 50)
    rec_max_det <- max(rec_max_det, 500L)
    max_det_rationale <- sprintf("CRITICAL: Peak density reaches %d objects/image! Set max_det=%d to prevent truncating detections (default YOLO cap is 300).", max_dens, rec_max_det)
  }

  model_name <- switch(task,
    detect = "yolov8m.pt",
    segment = "yolov8m-seg.pt",
    obb = "yolov8m-obb.pt",
    classify = "yolov8m-cls.pt",
    "yolov8m.pt"
  )

  yaml_str <- desc$yaml_file %||% "data.yaml"
  yaml_clean <- if (file.exists(yaml_str)) normalizePath(yaml_str, winslash = "/", mustWork = FALSE) else "path/to/data.yaml"

  r_code <- sprintf(
'# Ultralytics YOLO Training Command via pliman
library(pliman)

# Step 1: Train model with dataset-tailored hyperparameters
model <- yolo_train(
  data = "%s",
  model = "%s",
  epochs = %d,
  imgsz = %d,
  batch = %d,
  patience = %d,
  max_det = %d,
  plots = TRUE,
  save = TRUE
)

# Step 2: Evaluate validation benchmarks
val_metrics <- yolo_val(model = model, data = "%s")

# Step 3: Run inference on new test images
preds <- yolo_predict(
  model = model,
  source = "path/to/test_images",
  conf = 0.25,
  iou = 0.45,
  max_det = %d
)',
    yaml_clean, model_name, rec_epochs, rec_imgsz, rec_batch, rec_patience, rec_max_det,
    yaml_clean, rec_max_det
  )

  r_code_escaped <- gsub("&", "&amp;", gsub("<", "&lt;", gsub(">", "&gt;", r_code)))

  paste0(
    '<section class="section-block" id="advisor">\n',
    '  <div class="section-header">\n',
    '    <div>\n',
    '      <h2 class="section-title">\U0001f9e0 AI & Training Hyperparameter Advisor</h2>\n',
    '      <p class="section-desc">Data-driven hyperparameter configurations derived from this dataset\u2019s resolution, object scale distribution, and peak density.</p>\n',
    '    </div>\n',
    '  </div>\n',
    '  <div class="advisor-grid">\n',
    '    <div class="advisor-card card-glow">\n',
    '      <div class="advisor-param-title">Optimal Input Resolution (<code>imgsz</code>)</div>\n',
    '      <div class="advisor-param-val">', rec_imgsz, ' px</div>\n',
    '      <div class="advisor-param-desc">', imgsz_rationale, '</div>\n',
    '    </div>\n',
    '    <div class="advisor-card card-glow">\n',
    '      <div class="advisor-param-title">Maximum Detections (<code>max_det</code>)</div>\n',
    '      <div class="advisor-param-val" style="color: ', if (rec_max_det > 300) "var(--color-warning)" else "var(--color-primary)", ';">', rec_max_det, '</div>\n',
    '      <div class="advisor-param-desc">', max_det_rationale, '</div>\n',
    '    </div>\n',
    '    <div class="advisor-card card-glow">\n',
    '      <div class="advisor-param-title">Training Epochs & Patience</div>\n',
    '      <div class="advisor-param-val">', rec_epochs, ' <span style="font-size:1rem; font-weight:500; color:var(--text-muted);">epochs (pat ', rec_patience, ')</span></div>\n',
    '      <div class="advisor-param-desc">', epoch_rationale, '</div>\n',
    '    </div>\n',
    '    <div class="advisor-card card-glow">\n',
    '      <div class="advisor-param-title">Recommended Batch Size</div>\n',
    '      <div class="advisor-param-val">', rec_batch, '</div>\n',
    '      <div class="advisor-param-desc">', batch_rationale, '</div>\n',
    '    </div>\n',
    '  </div>\n',
    '  <div class="code-snippet-card">\n',
    '    <div class="code-snippet-header">\n',
    '      <div class="snippet-title"><span>\U0001f680 Ready-to-Execute R Training Script</span></div>\n',
    '      <button class="btn-copy" onclick="copyTrainingCode(this)" title="Copy script to clipboard"><span>\U0001f4cb Copy Script</span></button>\n',
    '    </div>\n',
    '    <pre class="code-snippet-body"><code id="training-code-block">', r_code_escaped, '</code></pre>\n',
    '  </div>\n',
    '</section>\n\n'
  )
}


# ==============================================================================
# INTERNAL HELPER: REPORT CSS STYLESHEET
# ==============================================================================
.yolo_report_css <- function(lwd = 1.5) {
  lwd_val <- suppressWarnings(as.numeric(lwd)[1])
  if (is.na(lwd_val) || lwd_val <= 0) lwd_val <- 1.5
  paste0('
    :root {
      --yolo-box-lwd: ', sprintf("%.1fpx", lwd_val), ';
      --font-body: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      --font-heading: "Outfit", "Inter", -apple-system, sans-serif;
      --font-mono: "JetBrains Mono", Consolas, Monaco, monospace;
      --color-primary: #6366f1;
      --color-primary-hover: #4f46e5;
      --color-accent: #0ea5e9;
      --color-success: #10b981;
      --color-warning: #f59e0b;
      --color-danger: #ef4444;
      --color-purple: #8b5cf6;
      --bg-canvas: #090d16;
      --bg-surface: #111827;
      --bg-surface-elevated: #1a2234;
      --bg-card: rgba(17, 24, 39, 0.85);
      --border-subtle: rgba(255, 255, 255, 0.08);
      --border-focus: rgba(99, 102, 241, 0.5);
      --text-main: #f3f4f6;
      --text-muted: #9ca3af;
      --text-faint: #6b7280;
      --shadow-card: 0 4px 20px -2px rgba(0, 0, 0, 0.5);
      --shadow-hover: 0 10px 30px -5px rgba(99, 102, 241, 0.25);
    }
    [data-theme="light"] {
      --bg-canvas: #f8fafc;
      --bg-surface: #ffffff;
      --bg-surface-elevated: #f1f5f9;
      --bg-card: rgba(255, 255, 255, 0.95);
      --border-subtle: #e2e8f0;
      --border-focus: rgba(99, 102, 241, 0.4);
      --text-main: #0f172a;
      --text-muted: #475569;
      --text-faint: #94a3b8;
      --shadow-card: 0 4px 15px -2px rgba(0, 0, 0, 0.06);
      --shadow-hover: 0 10px 25px -5px rgba(99, 102, 241, 0.15);
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    html { scroll-behavior: smooth; }
    body {
      font-family: var(--font-body);
      background-color: var(--bg-canvas);
      color: var(--text-main);
      line-height: 1.55;
      -webkit-font-smoothing: antialiased;
      transition: background-color 0.25s ease, color 0.25s ease;
    }
    .navbar {
      position: sticky;
      top: 0;
      z-index: 100;
      background: var(--bg-card);
      backdrop-filter: blur(16px);
      -webkit-backdrop-filter: blur(16px);
      border-bottom: 1px solid var(--border-subtle);
      padding: 0.85rem 2rem;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .nav-brand {
      display: flex;
      align-items: center;
      gap: 0.85rem;
      font-family: var(--font-heading);
      font-weight: 700;
      font-size: 1.25rem;
      letter-spacing: -0.02em;
    }
    .brand-logo {
      background: linear-gradient(135deg, #6366f1 0%, #ec4899 100%);
      color: #fff;
      padding: 0.35rem 0.65rem;
      border-radius: 8px;
      font-size: 0.95rem;
      font-weight: 800;
      box-shadow: 0 2px 8px rgba(99, 102, 241, 0.4);
    }
    .nav-menu {
      display: flex;
      gap: 1.25rem;
      list-style: none;
    }
    .nav-link {
      color: var(--text-muted);
      text-decoration: none;
      font-size: 0.88rem;
      font-weight: 500;
      padding: 0.35rem 0.65rem;
      border-radius: 6px;
      transition: all 0.2s ease;
    }
    .nav-link:hover, .nav-link.active {
      color: var(--color-primary);
      background: rgba(99, 102, 241, 0.1);
    }
    .nav-actions {
      display: flex;
      align-items: center;
      gap: 0.65rem;
    }
    .btn-action {
      background: var(--bg-surface-elevated);
      color: var(--text-main);
      border: 1px solid var(--border-subtle);
      padding: 0.45rem 0.85rem;
      border-radius: 8px;
      font-size: 0.85rem;
      font-weight: 600;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      transition: all 0.2s ease;
    }
    .btn-action:hover {
      border-color: var(--color-primary);
      color: var(--color-primary);
    }
    .container {
      max-width: 1380px;
      margin: 0 auto;
      padding: 2.5rem 1.5rem 5rem;
    }
    .header-banner { margin-bottom: 2rem; }
    .dataset-title {
      font-family: var(--font-heading);
      font-size: 2.35rem;
      font-weight: 800;
      letter-spacing: -0.03em;
      margin-bottom: 0.4rem;
      background: linear-gradient(135deg, var(--text-main) 60%, var(--color-primary) 100%);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
    }
    .dataset-meta {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 1.25rem;
      font-size: 0.88rem;
      color: var(--text-muted);
    }
    .badge-task {
      font-size: 0.75rem;
      font-weight: 700;
      padding: 0.25rem 0.65rem;
      border-radius: 9999px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }
    .badge-detect { background: rgba(99, 102, 241, 0.15); color: #818cf8; border: 1px solid rgba(99, 102, 241, 0.3); }
    .badge-segment { background: rgba(16, 185, 129, 0.15); color: #34d399; border: 1px solid rgba(16, 185, 129, 0.3); }
    .badge-obb { background: rgba(245, 158, 11, 0.15); color: #fbbf24; border: 1px solid rgba(245, 158, 11, 0.3); }
    .badge-classify { background: rgba(14, 165, 233, 0.15); color: #38bdf8; border: 1px solid rgba(14, 165, 233, 0.3); }

    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 1.25rem;
      margin-bottom: 1.25rem;
    }
    .kpi-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1.35rem 1.25rem;
      box-shadow: var(--shadow-card);
      transition: all 0.2s ease;
      position: relative;
      overflow: hidden;
    }
    .kpi-card:hover { transform: translateY(-2px); border-color: var(--border-focus); }
    .kpi-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.65rem;
    }
    .kpi-title {
      font-size: 0.78rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: var(--text-muted);
    }
    .kpi-icon { font-size: 1.15rem; }
    .kpi-value {
      font-family: var(--font-heading);
      font-size: 2.1rem;
      font-weight: 800;
      line-height: 1.1;
      color: var(--text-main);
      margin-bottom: 0.35rem;
    }
    .kpi-val-sub { font-size: 1rem; font-weight: 500; color: var(--text-muted); }
    .kpi-subtitle {
      font-size: 0.8rem;
      color: var(--text-faint);
    }
    .dot-sep { margin: 0 0.35rem; opacity: 0.5; }

    /* Health Score Circular Gauge */
    .kpi-gauge-card {
      display: flex;
      align-items: center;
      gap: 1.25rem;
    }
    .gauge-wrapper {
      position: relative;
      width: 82px;
      height: 82px;
      flex-shrink: 0;
    }
    .gauge-svg {
      transform: rotate(-90deg);
      width: 100%;
      height: 100%;
    }
    .gauge-bg {
      fill: none;
      stroke: var(--border-subtle);
      stroke-width: 8;
    }
    .gauge-fill {
      fill: none;
      stroke-width: 8;
      stroke-linecap: round;
      transition: stroke-dashoffset 1s ease;
    }
    .gauge-center {
      position: absolute;
      top: 50%; left: 50%;
      transform: translate(-50%, -50%);
      text-align: center;
    }
    .gauge-score {
      font-family: var(--font-heading);
      font-size: 1.35rem;
      font-weight: 800;
      line-height: 1;
      color: var(--text-main);
    }
    .gauge-grade {
      font-size: 0.7rem;
      font-weight: 800;
      letter-spacing: 0.05em;
      text-transform: uppercase;
      margin-top: 2px;
    }
    .grade-opt { color: #10b981; stroke: #10b981; }
    .grade-good { color: #0ea5e9; stroke: #0ea5e9; }
    .grade-med { color: #f59e0b; stroke: #f59e0b; }
    .grade-warn { color: #ef4444; stroke: #ef4444; }

    /* Anomaly Strip */
    .anomaly-strip-container {
      margin-bottom: 2.5rem;
      padding: 0.85rem 1.25rem;
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      display: flex;
      align-items: center;
      gap: 1rem;
      flex-wrap: wrap;
    }
    .anomaly-strip-title {
      font-size: 0.8rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-muted);
      flex-shrink: 0;
    }
    .anomaly-pills {
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      align-items: center;
    }
    .anomaly-pill {
      display: inline-flex;
      align-items: center;
      gap: 0.35rem;
      padding: 0.3rem 0.7rem;
      border-radius: 9999px;
      font-size: 0.76rem;
      font-weight: 600;
      transition: all 0.15s ease;
      font-family: inherit;
    }
    .pill-clickable {
      cursor: pointer;
      border: 1px solid transparent;
      outline: none;
    }
    .pill-clickable:hover {
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(0,0,0,0.3);
      filter: brightness(1.15);
    }
    .pill-inspect {
      font-size: 0.65rem;
      opacity: 0.8;
      text-transform: uppercase;
      margin-left: 0.2rem;
      border-bottom: 1px dotted currentColor;
    }
    .pill-perfect { background: rgba(16, 185, 129, 0.15); color: #34d399; border: 1px solid rgba(16, 185, 129, 0.3); }
    .pill-pass { background: rgba(16, 185, 129, 0.12); color: #34d399; }
    .pill-warn { background: rgba(245, 158, 11, 0.15); color: #fbbf24; border: 1px solid rgba(245, 158, 11, 0.3); }
    .pill-danger { background: rgba(239, 68, 68, 0.15); color: #f87171; border: 1px solid rgba(239, 68, 68, 0.3); }
    .pill-info { background: rgba(14, 165, 233, 0.15); color: #38bdf8; border: 1px solid rgba(14, 165, 233, 0.3); }
    .btn-inspect-all {
      margin-left: auto;
      background: var(--bg-surface-elevated);
      border: 1px solid var(--border-subtle);
      color: var(--text-main);
      font-size: 0.8rem;
      padding: 0.35rem 0.75rem;
      border-radius: 8px;
      cursor: pointer;
    }
    .btn-inspect-all:hover {
      border-color: var(--color-primary);
      color: var(--color-primary);
    }

    .section-block {
      margin-bottom: 3.5rem;
      scroll-margin-top: 5rem;
    }
    .section-header {
      display: flex;
      justify-content: space-between;
      align-items: flex-end;
      margin-bottom: 1.25rem;
      border-bottom: 1px solid var(--border-subtle);
      padding-bottom: 0.85rem;
    }
    .section-title {
      font-family: var(--font-heading);
      font-size: 1.45rem;
      font-weight: 700;
      letter-spacing: -0.02em;
    }
    .section-desc {
      font-size: 0.85rem;
      color: var(--text-muted);
      margin-top: 0.2rem;
    }

    .prov-bar-container {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1.25rem;
      margin-bottom: 1.5rem;
    }
    .prov-bar {
      display: flex;
      height: 24px;
      border-radius: 6px;
      overflow: hidden;
      background: var(--bg-surface-elevated);
      margin-bottom: 0.85rem;
    }
    .prov-seg {
      height: 100%;
      transition: width 0.4s ease;
      position: relative;
    }
    .prov-orig { background: linear-gradient(90deg, #6366f1, #4f46e5); }
    .prov-aug { background: linear-gradient(90deg, #8b5cf6, #7c3aed); }
    .prov-synth { background: linear-gradient(90deg, #f59e0b, #d97706); }
    .prov-tile { background: linear-gradient(90deg, #10b981, #059669); }
    .prov-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 1.5rem;
      font-size: 0.82rem;
    }
    .legend-item { display: inline-flex; align-items: center; gap: 0.5rem; }
    .legend-color { width: 12px; height: 12px; border-radius: 3px; }
    .prov-orig-bg { background: #6366f1; }
    .prov-aug-bg { background: #8b5cf6; }
    .prov-synth-bg { background: #f59e0b; }
    .prov-tile-bg { background: #10b981; }

    .table-responsive {
      width: 100%;
      overflow-x: auto;
      border-radius: 12px;
      border: 1px solid var(--border-subtle);
      background: var(--bg-surface);
      box-shadow: var(--shadow-card);
    }
    .report-table {
      width: 100%;
      border-collapse: collapse;
      font-size: 0.85rem;
      text-align: left;
    }
    .report-table th {
      background: var(--bg-surface-elevated);
      color: var(--text-muted);
      font-weight: 600;
      padding: 0.75rem 1rem;
      border-bottom: 1px solid var(--border-subtle);
      text-transform: uppercase;
      font-size: 0.75rem;
      letter-spacing: 0.05em;
    }
    .report-table td {
      padding: 0.75rem 1rem;
      border-bottom: 1px solid var(--border-subtle);
      color: var(--text-main);
    }
    .report-table tr:last-child td { border-bottom: none; }
    .report-table tr:hover td { background: rgba(99, 102, 241, 0.04); }
    .row-total td {
      font-weight: 700;
      background: var(--bg-surface-elevated);
      border-top: 2px solid var(--border-subtle);
    }
    .class-id-badge {
      display: inline-block;
      width: 22px;
      height: 22px;
      line-height: 22px;
      text-align: center;
      border-radius: 5px;
      font-size: 0.75rem;
      font-weight: 700;
      color: #fff;
    }
    .table-bar-cell {
      display: flex;
      align-items: center;
      gap: 0.75rem;
    }
    .mini-bar-bg {
      flex: 1;
      height: 6px;
      background: var(--bg-surface-elevated);
      border-radius: 3px;
      overflow: hidden;
      min-width: 60px;
    }
    .mini-bar-fill {
      height: 100%;
      border-radius: 3px;
      transition: width 0.3s ease;
    }
    .tag-badge {
      display: inline-block;
      padding: 0.2rem 0.55rem;
      border-radius: 6px;
      font-size: 0.75rem;
      font-weight: 600;
    }
    .prov-badge-orig { background: rgba(99, 102, 241, 0.15); color: #818cf8; }
    .prov-badge-aug { background: rgba(139, 92, 246, 0.15); color: #a78bfa; }
    .prov-badge-synth { background: rgba(245, 158, 11, 0.15); color: #fbbf24; }
    .prov-badge-tile { background: rgba(16, 185, 129, 0.15); color: #34d399; }
    .prov-badge-tot { background: var(--bg-surface-elevated); color: var(--text-muted); }
    .size-small { background: rgba(14, 165, 233, 0.15); color: #38bdf8; }
    .size-medium { background: rgba(168, 85, 247, 0.15); color: #c084fc; }
    .size-large { background: rgba(236, 72, 153, 0.15); color: #f472b6; }

    /* Spatial Heatmap & Dimension Grid */
    .chart-caption {
      font-size: 0.8rem;
      color: var(--text-muted);
      margin-bottom: 1.25rem;
      line-height: 1.4;
    }
    .card-head-row {
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.4rem;
    }
    .info-badge {
      font-size: 0.72rem;
      font-weight: 600;
      padding: 0.2rem 0.5rem;
      border-radius: 6px;
      background: var(--bg-surface-elevated);
      color: var(--text-muted);
      border: 1px solid var(--border-subtle);
    }
    .spatial-heatmap-wrapper {
      display: grid;
      grid-template-columns: 36px 1fr;
      grid-template-rows: 36px 1fr;
      gap: 6px;
      width: 100%;
      max-width: 440px;
      margin: 0 auto;
      aspect-ratio: 1 / 1;
    }
    .marg-corner { width: 36px; height: 36px; }
    .marg-top-container {
      display: grid;
      grid-template-columns: repeat(12, 1fr);
      gap: 2px;
      align-items: end;
      height: 36px;
    }
    .marg-col {
      height: 100%;
      display: flex;
      align-items: flex-end;
      background: rgba(255,255,255,0.03);
      border-radius: 2px;
      cursor: pointer;
    }
    .marg-col-bar {
      width: 100%;
      background: var(--color-accent);
      border-radius: 2px 2px 0 0;
      transition: height 0.3s ease;
    }
    .marg-left-container {
      display: grid;
      grid-template-rows: repeat(12, 1fr);
      gap: 2px;
      justify-items: end;
      width: 36px;
    }
    .marg-row {
      width: 100%;
      display: flex;
      justify-content: flex-end;
      background: rgba(255,255,255,0.03);
      border-radius: 2px;
      cursor: pointer;
    }
    .marg-row-bar {
      height: 100%;
      background: var(--color-accent);
      border-radius: 2px 0 0 2px;
      transition: width 0.3s ease;
    }
    .heatmap-grid-12x12 {
      position: relative;
      display: grid;
      grid-template-columns: repeat(12, 1fr);
      grid-template-rows: repeat(12, 1fr);
      gap: 2px;
      background: rgba(0, 0, 0, 0.4);
      border: 1px solid var(--border-subtle);
      border-radius: 6px;
      overflow: hidden;
      aspect-ratio: 1 / 1;
    }
    .heatmap-cell {
      position: relative;
      width: 100%;
      height: 100%;
      border-radius: 1px;
      cursor: crosshair;
      transition: transform 0.1s ease, filter 0.1s ease;
    }
    .heatmap-cell:hover {
      filter: brightness(1.4) drop-shadow(0 0 4px var(--color-primary));
      z-index: 10;
      transform: scale(1.15);
    }
    .heatmap-crosshair-v {
      position: absolute;
      top: 0; bottom: 0; left: 50%;
      width: 1px;
      background: rgba(255, 255, 255, 0.2);
      pointer-events: none;
      z-index: 5;
    }
    .heatmap-crosshair-h {
      position: absolute;
      left: 0; right: 0; top: 50%;
      height: 1px;
      background: rgba(255, 255, 255, 0.2);
      pointer-events: none;
      z-index: 5;
    }
    .heatmap-legend {
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 0.75rem;
      margin-top: 1rem;
      font-size: 0.75rem;
      color: var(--text-muted);
    }
    .legend-gradient {
      width: 140px;
      height: 8px;
      border-radius: 4px;
      background: linear-gradient(to right, rgba(99, 102, 241, 0.1), rgba(99, 102, 241, 0.8), #0ea5e9, #f59e0b);
    }

    .dim-heatmap-container {
      display: flex;
      align-items: center;
      gap: 0.5rem;
      max-width: 380px;
      margin: 0 auto;
    }
    .dim-y-axis-label {
      writing-mode: vertical-rl;
      transform: rotate(180deg);
      font-size: 0.75rem;
      font-weight: 600;
      color: var(--text-muted);
      text-align: center;
    }
    .dim-heatmap-inner {
      flex: 1;
      display: flex;
      flex-direction: column;
      gap: 0.4rem;
    }
    .heatmap-grid-10x10 {
      position: relative;
      display: grid;
      grid-template-columns: repeat(10, 1fr);
      grid-template-rows: repeat(10, 1fr);
      gap: 2px;
      background: rgba(0, 0, 0, 0.4);
      border: 1px solid var(--border-subtle);
      border-radius: 6px;
      overflow: hidden;
      aspect-ratio: 1 / 1;
    }
    .diagonal-reference-line {
      position: absolute;
      top: 0; left: 0; width: 100%; height: 100%;
      background: linear-gradient(to top right, transparent calc(50% - 1px), rgba(255,255,255,0.3) 50%, transparent calc(50% + 1px));
      pointer-events: none;
      z-index: 5;
    }
    .dim-x-axis-label {
      text-align: center;
      font-size: 0.75rem;
      font-weight: 600;
      color: var(--text-muted);
    }

    /* Training Advisor & Code Snippet */
    .advisor-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 1.25rem;
      margin-bottom: 1.5rem;
    }
    .advisor-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1.25rem;
      display: flex;
      flex-direction: column;
      justify-content: space-between;
    }
    .advisor-param-title {
      font-size: 0.78rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-muted);
      margin-bottom: 0.5rem;
    }
    .advisor-param-val {
      font-family: var(--font-heading);
      font-size: 1.75rem;
      font-weight: 800;
      color: var(--color-primary);
      line-height: 1.2;
      margin-bottom: 0.4rem;
    }
    .advisor-param-desc {
      font-size: 0.78rem;
      color: var(--text-faint);
      line-height: 1.45;
    }
    .code-snippet-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      overflow: hidden;
      box-shadow: var(--shadow-card);
    }
    .code-snippet-header {
      background: var(--bg-surface-elevated);
      border-bottom: 1px solid var(--border-subtle);
      padding: 0.75rem 1.25rem;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .snippet-title {
      font-family: var(--font-heading);
      font-size: 0.9rem;
      font-weight: 700;
      color: var(--text-main);
    }
    .btn-copy {
      background: var(--color-primary);
      color: #fff;
      border: none;
      padding: 0.35rem 0.75rem;
      border-radius: 6px;
      font-size: 0.78rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.2s;
    }
    .btn-copy:hover {
      background: var(--color-primary-hover);
      transform: translateY(-1px);
    }
    .code-snippet-body {
      padding: 1.25rem;
      overflow-x: auto;
      font-family: var(--font-mono);
      font-size: 0.85rem;
      line-height: 1.6;
      color: #a5b4fc;
      background: rgba(0, 0, 0, 0.4);
      margin: 0;
    }

    /* Co-occurrence Table */
    .cooc-table {
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
      text-align: center;
    }
    .cooc-table th, .cooc-table td {
      padding: 0.65rem 0.5rem;
      border: 1px solid var(--border-subtle);
    }
    .cooc-row-header {
      text-align: left;
      background: var(--bg-surface-elevated);
      font-weight: 600;
    }
    .cooc-diag {
      font-weight: 700;
      border: 1px solid var(--color-primary) !important;
    }
    .th-class-rot {
      max-width: 90px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }

    /* Enhanced Gallery Toolbar */
    .gallery-toolbar {
      display: flex;
      flex-direction: column;
      gap: 1rem;
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1rem 1.25rem;
      margin-bottom: 1.5rem;
    }
    .gallery-toolbar-top {
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 0.75rem;
    }
    .gallery-tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
    }
    .tab-btn {
      background: var(--bg-surface-elevated);
      color: var(--text-muted);
      border: 1px solid var(--border-subtle);
      padding: 0.4rem 0.85rem;
      border-radius: 8px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.2s ease;
    }
    .tab-btn:hover { color: var(--text-main); border-color: var(--color-primary); }
    .tab-btn.active {
      background: var(--color-primary);
      color: #ffffff;
      border-color: var(--color-primary);
      box-shadow: 0 2px 8px rgba(99, 102, 241, 0.35);
    }
    .gallery-toolbar-bottom {
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 1rem;
      border-top: 1px solid var(--border-subtle);
      padding-top: 0.85rem;
    }
    .gallery-filter-group {
      display: flex;
      align-items: center;
      gap: 0.65rem;
      flex-wrap: wrap;
    }
    .gallery-search-input {
      background: var(--bg-surface-elevated);
      border: 1px solid var(--border-subtle);
      color: var(--text-main);
      padding: 0.4rem 0.85rem;
      border-radius: 8px;
      font-size: 0.82rem;
      font-family: var(--font-body);
      min-width: 180px;
      outline: none;
      transition: border-color 0.2s;
    }
    .gallery-search-input:focus { border-color: var(--color-primary); }
    .gallery-select {
      background: var(--bg-surface-elevated);
      border: 1px solid var(--border-subtle);
      color: var(--text-main);
      padding: 0.4rem 0.75rem;
      border-radius: 8px;
      font-size: 0.82rem;
      cursor: pointer;
      outline: none;
    }
    .gallery-count-badge {
      font-size: 0.78rem;
      font-weight: 600;
      color: var(--text-muted);
      padding: 0.25rem 0.6rem;
      border-radius: 6px;
      background: var(--bg-surface-elevated);
    }
    .gallery-toggles {
      display: flex;
      align-items: center;
      gap: 1.25rem;
      font-size: 0.82rem;
      color: var(--text-muted);
      font-weight: 500;
    }
    .toggle-label {
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      cursor: pointer;
      user-select: none;
    }
    .toggle-label input { accent-color: var(--color-primary); cursor: pointer; }

    .gallery-grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
      gap: 1.35rem;
    }
    .gallery-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      overflow: hidden;
      box-shadow: var(--shadow-card);
      transition: transform 0.2s ease, box-shadow 0.2s ease, border-color 0.2s ease;
      display: flex;
      flex-direction: column;
    }
    .gallery-card:hover {
      transform: translateY(-3px);
      box-shadow: var(--shadow-hover);
      border-color: var(--border-focus);
    }
    .card-media-wrapper {
      position: relative;
      width: 100%;
      aspect-ratio: 1 / 1;
      background: #000000;
      overflow: hidden;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .card-media-frame {
      position: relative;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      max-width: 100%;
      max-height: 100%;
      line-height: 0;
    }
    .card-thumb {
      max-width: 100%;
      max-height: 100%;
      width: auto;
      height: auto;
      object-fit: contain;
      display: block;
      transition: transform 0.3s ease;
    }
    .card-media-wrapper:hover .card-thumb { transform: scale(1.03); }
    .card-media-frame .yolo-svg-overlay,
    .yolo-svg-overlay {
      position: absolute;
      top: 0; left: 0;
      width: 100%; height: 100%;
      pointer-events: none;
    }
    .yolo-geo,
    .yolo-bbox,
    .yolo-poly {
      vector-effect: non-scaling-stroke;
      stroke-width: var(--yolo-box-lwd, 1.5px) !important;
    }
    .card-badges {
      position: absolute;
      top: 0.65rem; left: 0.65rem;
      display: flex;
      flex-wrap: wrap;
      gap: 0.4rem;
      pointer-events: none;
      max-width: 90%;
    }
    .card-badge {
      font-size: 0.68rem;
      font-weight: 700;
      padding: 0.2rem 0.5rem;
      border-radius: 6px;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      box-shadow: 0 2px 6px rgba(0,0,0,0.4);
    }
    .badge-original { background: #6366f1; color: #fff; }
    .badge-augmented { background: #8b5cf6; color: #fff; }
    .badge-synthetic_bg { background: #f59e0b; color: #000; }
    .badge-tiled { background: #10b981; color: #fff; }
    .badge-background { background: #64748b; color: #fff; }
    .badge-split { background: rgba(0, 0, 0, 0.65); color: #fff; border: 1px solid rgba(255,255,255,0.2); }
    .badge-anomaly { background: #ef4444 !important; color: #fff !important; }
    .badge-blur { background: #dc2626 !important; color: #fff !important; border: 1px solid rgba(255,255,255,0.25); }
    .meta-blur-warning { color: #f87171 !important; background: rgba(239, 68, 68, 0.15) !important; font-weight: 700; }
    .meta-blur-sharp { color: #10b981 !important; background: rgba(16, 185, 129, 0.12) !important; }
    .btn-inspect-blur:hover { background: rgba(239, 68, 68, 0.22) !important; border-color: #ef4444 !important; transform: translateY(-1px); }

    /* Augmentation Filter Pills */
    .gallery-toolbar-aug {
      margin-top: 0.85rem;
      padding-top: 0.85rem;
      border-top: 1px solid var(--border-subtle);
      display: flex;
      flex-direction: column;
      gap: 0.5rem;
    }
    .aug-toolbar-header {
      display: flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 0.6rem;
    }
    .aug-toolbar-label {
      font-weight: 700;
      color: var(--text-main);
      text-transform: uppercase;
      letter-spacing: 0.04em;
      font-size: 0.78rem;
    }
    .aug-toolbar-desc {
      color: var(--text-muted);
      font-size: 0.76rem;
    }
    .aug-clear-btn {
      background: rgba(239, 68, 68, 0.12);
      color: #ef4444;
      border: 1px solid rgba(239, 68, 68, 0.35);
      border-radius: 9999px;
      padding: 0.2rem 0.65rem;
      font-size: 0.72rem;
      font-weight: 700;
      cursor: pointer;
      transition: all 0.2s;
      margin-left: auto;
    }
    .aug-clear-btn:hover { background: #ef4444; color: #fff; }
    .aug-pills-container {
      display: flex;
      flex-wrap: wrap;
      gap: 0.45rem;
    }
    .aug-pill {
      background: var(--bg-surface-elevated);
      color: var(--text-main);
      border: 1px solid var(--border-subtle);
      border-radius: 9999px;
      padding: 0.32rem 0.75rem;
      font-size: 0.8rem;
      font-weight: 600;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 0.45rem;
      transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
      user-select: none;
    }
    .aug-pill:hover {
      border-color: var(--color-primary);
      transform: translateY(-1px);
    }
    .aug-pill .pill-count {
      font-size: 0.72rem;
      font-weight: 700;
      padding: 0.1rem 0.42rem;
      border-radius: 9999px;
      background: rgba(255, 255, 255, 0.08);
      color: var(--text-muted);
    }
    .aug-pill.active {
      background: var(--color-primary) !important;
      color: #ffffff !important;
      border-color: var(--color-primary) !important;
      box-shadow: 0 4px 12px rgba(99, 102, 241, 0.35);
    }
    .aug-pill.active .pill-count {
      background: rgba(0, 0, 0, 0.28) !important;
      color: #ffffff !important;
    }
    .card-aug-tag {
      font-size: 0.66rem;
      font-weight: 600;
      padding: 0.12rem 0.35rem;
      border-radius: 4px;
      background: rgba(99, 102, 241, 0.12);
      color: #818cf8;
      border: 1px solid rgba(99, 102, 241, 0.25);
      white-space: nowrap;
    }

    .card-footer {
      padding: 0.85rem 1rem;
      background: var(--bg-surface);
      border-top: 1px solid var(--border-subtle);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .card-stem {
      font-family: var(--font-mono);
      font-size: 0.78rem;
      font-weight: 500;
      color: var(--text-main);
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      max-width: 160px;
    }
    .card-meta {
      display: flex;
      align-items: center;
      gap: 0.4rem;
    }
    .meta-tag {
      font-size: 0.72rem;
      font-weight: 600;
      padding: 0.15rem 0.4rem;
      border-radius: 4px;
      background: var(--bg-surface-elevated);
      color: var(--text-muted);
    }
    .obj-tag { color: var(--color-primary); background: rgba(99, 102, 241, 0.1); }

    .hide-bboxes .yolo-geo,
    .hide-bboxes .yolo-label { display: none !important; }
    .hide-labels .yolo-label { display: none !important; }
    .hide-fills .yolo-geo { fill-opacity: 0 !important; }

    .health-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 1.25rem;
    }
    .health-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1.25rem;
      display: flex;
      gap: 1rem;
      box-shadow: var(--shadow-card);
      transition: all 0.2s ease;
    }
    .health-card:hover { transform: translateY(-2px); border-color: var(--border-focus); }
    .health-icon {
      width: 38px;
      height: 38px;
      border-radius: 10px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-size: 1.15rem;
      font-weight: 800;
      flex-shrink: 0;
    }
    .health-pass .health-icon { background: rgba(16, 185, 129, 0.15); color: #10b981; }
    .health-warning .health-icon { background: rgba(245, 158, 11, 0.15); color: #f59e0b; }
    .health-info .health-icon { background: rgba(14, 165, 233, 0.15); color: #0ea5e9; }
    .health-body { flex: 1; }
    .health-head {
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.35rem;
    }
    .health-title { font-weight: 700; font-size: 0.95rem; }
    .health-badge {
      font-size: 0.68rem;
      font-weight: 700;
      padding: 0.2rem 0.5rem;
      border-radius: 6px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .badge-health-pass { background: rgba(16, 185, 129, 0.15); color: #34d399; }
    .badge-health-warning { background: rgba(245, 158, 11, 0.15); color: #fbbf24; }
    .badge-health-info { background: rgba(14, 165, 233, 0.15); color: #38bdf8; }
    .health-metric {
      font-family: var(--font-mono);
      font-size: 0.85rem;
      color: var(--color-primary);
      margin-bottom: 0.4rem;
      font-weight: 600;
    }
    .health-detail {
      font-size: 0.8rem;
      color: var(--text-muted);
      line-height: 1.45;
    }

    /* Lightbox Modal */
    .lightbox-modal {
      display: none;
      position: fixed;
      top: 0; left: 0;
      width: 100%; height: 100%;
      background: rgba(0, 0, 0, 0.88);
      backdrop-filter: blur(8px);
      z-index: 1000;
      align-items: center;
      justify-content: center;
      padding: 2rem;
    }
    .lightbox-modal.open { display: flex; }
    .lightbox-content {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 16px;
      max-width: 1100px;
      width: 100%;
      max-height: 90vh;
      display: flex;
      overflow: hidden;
      box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.7);
      position: relative;
    }
    .lightbox-preview {
      flex: 2;
      background: #000;
      display: flex;
      align-items: center;
      justify-content: center;
      position: relative;
      overflow: hidden;
      min-height: 380px;
    }
    .lightbox-media-frame {
      position: relative;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      max-width: 100%;
      max-height: 80vh;
      line-height: 0;
    }
    .lightbox-img {
      max-width: 100%;
      max-height: 80vh;
      width: auto;
      height: auto;
      object-fit: contain;
      display: block;
    }
    .lightbox-svg-layer {
      position: absolute;
      top: 0; left: 0;
      width: 100%; height: 100%;
      pointer-events: none;
    }
    .lightbox-sidebar {
      flex: 1;
      padding: 2rem;
      overflow-y: auto;
      max-height: 80vh;
      border-left: 1px solid var(--border-subtle);
      display: flex;
      flex-direction: column;
    }
    .lightbox-close {
      position: absolute;
      top: 1rem; right: 1rem;
      width: 36px; height: 36px;
      background: var(--bg-surface-elevated);
      color: var(--text-main);
      border: 1px solid var(--border-subtle);
      border-radius: 50%;
      display: flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      font-size: 1.2rem;
      z-index: 10;
      transition: all 0.2s;
    }
    .lightbox-close:hover { background: var(--color-danger); color: #fff; }
    .lightbox-nav {
      position: absolute;
      top: 50%; transform: translateY(-50%);
      width: 44px; height: 44px;
      background: rgba(0, 0, 0, 0.65);
      color: #fff;
      border: 1px solid rgba(255,255,255,0.2);
      border-radius: 50%;
      display: flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      font-size: 1.4rem;
      transition: all 0.2s;
      user-select: none;
    }
    .lightbox-nav:hover { background: var(--color-primary); }
    .lightbox-prev { left: 1rem; }
    .lightbox-next { right: 1rem; }

    .lightbox-anomaly-box {
      background: rgba(239, 68, 68, 0.1);
      border: 1px solid rgba(239, 68, 68, 0.3);
      border-radius: 8px;
      padding: 0.75rem 1rem;
      margin-bottom: 1.25rem;
      font-size: 0.8rem;
      color: #fca5a5;
    }
    .lightbox-anomaly-box ul {
      margin-left: 1.2rem;
      margin-top: 0.35rem;
    }

    /* Lightbox Overlay Display Controls */
    .lightbox-overlay-controls {
      background: var(--bg-surface-elevated);
      border: 1px solid var(--border-subtle);
      border-radius: 10px;
      padding: 0.65rem 0.85rem;
      margin-bottom: 0.95rem;
      display: flex;
      flex-direction: column;
      gap: 0.55rem;
      box-shadow: 0 2px 6px rgba(0, 0, 0, 0.15);
    }
    .lightbox-controls-title {
      display: flex;
      align-items: center;
      justify-content: space-between;
      font-size: 0.72rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-muted);
    }
    .lb-keys-hint {
      font-size: 0.65rem;
      font-weight: 500;
      color: var(--text-muted);
      opacity: 0.85;
      text-transform: none;
      letter-spacing: normal;
    }
    .lb-keys-hint b {
      color: var(--color-primary);
      font-weight: 700;
      padding: 0 0.15rem;
    }
    .lightbox-toggle-group {
      display: flex;
      gap: 0.4rem;
      flex-wrap: wrap;
    }
    .lb-control-pill {
      display: inline-flex;
      align-items: center;
      gap: 0.35rem;
      padding: 0.35rem 0.65rem;
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 6px;
      font-size: 0.76rem;
      font-weight: 600;
      color: var(--text-main);
      cursor: pointer;
      user-select: none;
      transition: all 0.15s ease;
      flex: 1 1 auto;
      justify-content: center;
    }
    .lb-control-pill:hover {
      border-color: var(--color-primary);
      transform: translateY(-1px);
    }
    .lb-control-pill:has(input:checked),
    .lb-control-pill.active {
      background: rgba(99, 102, 241, 0.15);
      border-color: var(--color-primary);
      color: var(--color-primary);
    }
    .lb-control-pill:has(input:not(:checked)) {
      opacity: 0.65;
      color: var(--text-muted);
    }
    .lb-control-pill input {
      accent-color: var(--color-primary);
      cursor: pointer;
      margin: 0;
      width: 14px;
      height: 14px;
    }
    .lightbox-slider-line {
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding-top: 0.35rem;
      border-top: 1px solid var(--border-subtle);
    }
    .lightbox-slider-label {
      font-size: 0.72rem;
      font-weight: 600;
      color: var(--text-muted);
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }

    /* Anomaly Inspector Modal */
    .anomaly-modal-content {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 16px;
      max-width: 1050px;
      width: 100%;
      max-height: 85vh;
      display: flex;
      flex-direction: column;
      overflow: hidden;
      box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.7);
      position: relative;
      padding: 1.75rem;
    }
    .anomaly-modal-header {
      margin-bottom: 1.25rem;
    }
    .anomaly-modal-title {
      font-family: var(--font-heading);
      font-size: 1.35rem;
      font-weight: 700;
      margin-bottom: 0.3rem;
    }
    .anomaly-modal-desc {
      font-size: 0.82rem;
      color: var(--text-muted);
      margin-bottom: 1rem;
    }
    .anomaly-filter-tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin-bottom: 0.85rem;
    }
    .anomaly-thumb {
      width: 52px;
      height: 52px;
      object-fit: cover;
      border-radius: 8px;
      border: 1px solid var(--border-subtle);
      cursor: pointer;
      transition: transform 0.2s;
      display: block;
    }
    .anomaly-thumb:hover { transform: scale(1.15); border-color: var(--color-primary); }
    .anomaly-thumb-placeholder {
      width: 52px; height: 52px;
      background: var(--bg-surface-elevated);
      border-radius: 8px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-size: 1.2rem;
      cursor: pointer;
    }
    .anomaly-stem-code {
      font-family: var(--font-mono);
      font-size: 0.82rem;
      color: var(--text-main);
    }
    .anomaly-detail-text {
      font-size: 0.82rem;
      color: var(--text-muted);
      line-height: 1.4;
    }

    .grid-2col {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 1.5rem;
    }
    .content-card {
      background: var(--bg-surface);
      border: 1px solid var(--border-subtle);
      border-radius: 12px;
      padding: 1.25rem;
      box-shadow: var(--shadow-card);
    }
    .card-subtitle {
      font-size: 0.95rem;
      font-weight: 700;
      margin-bottom: 0.85rem;
      color: var(--text-main);
    }

    .report-footer {
      margin-top: 4rem;
      padding-top: 2rem;
      border-top: 1px solid var(--border-subtle);
      display: flex;
      justify-content: space-between;
      align-items: center;
      font-size: 0.85rem;
      color: var(--text-faint);
    }

    @media print {
      .navbar, .gallery-toolbar, .btn-action, .lightbox-modal, .btn-copy, .btn-inspect-all { display: none !important; }
      body { background: #fff !important; color: #000 !important; }
      .kpi-card, .table-responsive, .health-card, .advisor-card { box-shadow: none !important; border: 1px solid #ccc !important; }
    }
    @media (max-width: 900px) {
      .lightbox-content { flex-direction: column; max-height: 95vh; }
      .lightbox-sidebar { border-left: none; border-top: 1px solid var(--border-subtle); max-height: 40vh; }
      .navbar { padding: 0.85rem 1rem; }
      .nav-menu { display: none; }
    }
  ')
}


# ==============================================================================
# INTERNAL HELPER: REPORT JAVASCRIPT LOGIC
# ==============================================================================
.yolo_report_js <- function() {
  '
    function updateBoxLwd(val) {
      var pxVal = parseFloat(val);
      if (isNaN(pxVal) || pxVal <= 0) pxVal = 1.5;
      document.documentElement.style.setProperty("--yolo-box-lwd", pxVal + "px");
      var el1 = document.getElementById("box-lwd-val");
      if (el1) el1.textContent = pxVal.toFixed(1) + "px";
      var s1 = document.getElementById("box-lwd-slider");
      if (s1 && Math.abs(parseFloat(s1.value) - pxVal) > 0.01) s1.value = pxVal;
      var el2 = document.getElementById("lightboxLwdVal");
      if (el2) el2.textContent = pxVal.toFixed(1) + "px";
      var s2 = document.getElementById("lightboxLwdSlider");
      if (s2 && Math.abs(parseFloat(s2.value) - pxVal) > 0.01) s2.value = pxVal;
      try { localStorage.setItem("pliman_report_box_lwd", pxVal); } catch(e) {}
    }

    function toggleTheme() {
      var htmlEl = document.documentElement;
      var cur = htmlEl.getAttribute("data-theme") || "dark";
      var next = cur === "dark" ? "light" : "dark";
      htmlEl.setAttribute("data-theme", next);
      document.getElementById("theme-icon").textContent = next === "dark" ? "\\u{1f319}" : "\\u{2600}\\u{fe0f}";
      document.getElementById("theme-text").textContent = next === "dark" ? "Dark" : "Light";
      try { localStorage.setItem("pliman_report_theme", next); } catch(e) {}
    }

    var activeAugFilters = [];

    function toggleAugFilter(tag, btn) {
      var idx = activeAugFilters.indexOf(tag);
      if (idx !== -1) {
        activeAugFilters.splice(idx, 1);
        if (btn) btn.classList.remove("active");
      } else {
        activeAugFilters.push(tag);
        if (btn) btn.classList.add("active");
      }
      var clearBtn = document.getElementById("aug-clear-btn");
      if (clearBtn) {
        clearBtn.style.display = activeAugFilters.length > 0 ? "inline-flex" : "none";
      }
      filterGalleryCombined();
    }

    function clearAugFilters() {
      activeAugFilters = [];
      document.querySelectorAll(".aug-pill").forEach(function(b) { b.classList.remove("active"); });
      var clearBtn = document.getElementById("aug-clear-btn");
      if (clearBtn) clearBtn.style.display = "none";
      filterGalleryCombined();
    }

    function filterGalleryCombined() {
      var activeTab = document.querySelector(".tab-btn.active");
      var catFilter = activeTab ? (activeTab.getAttribute("data-category") || "all") : "all";
      var splitSelect = document.getElementById("gallery-split-select");
      var splitFilter = splitSelect ? splitSelect.value : "all";
      var searchInput = document.getElementById("gallery-search");
      var searchVal = searchInput ? searchInput.value.toLowerCase().trim() : "";

      var cards = document.querySelectorAll(".gallery-card");
      var visibleCount = 0;

      cards.forEach(function(c) {
        var cCat = c.getAttribute("data-category") || "";
        var cSplit = (c.getAttribute("data-split") || "").toLowerCase();
        var cStem = (c.getAttribute("data-stem") || "").toLowerCase();
        var hasAnom = c.getAttribute("data-has-anomaly") === "true";
        var isBlur = c.getAttribute("data-is-blurry") === "true";
        var cAugStr = (c.getAttribute("data-aug-tags") || "").trim();
        var cAugList = cAugStr ? cAugStr.split(/\\s+/) : [];

        var matchCat = false;
        if (catFilter === "all") {
          matchCat = true;
        } else if (catFilter === "anomaly") {
          matchCat = hasAnom;
        } else if (catFilter === "blurry") {
          matchCat = isBlur;
        } else {
          matchCat = (cCat === catFilter);
        }

        var matchSplit = (splitFilter === "all" || cSplit === splitFilter.toLowerCase());
        var matchSearch = (!searchVal || cStem.indexOf(searchVal) !== -1);

        var matchAug = true;
        if (activeAugFilters.length > 0) {
          for (var a = 0; a < activeAugFilters.length; a++) {
            if (cAugList.indexOf(activeAugFilters[a]) === -1) {
              matchAug = false;
              break;
            }
          }
        }

        if (matchCat && matchSplit && matchSearch && matchAug) {
          c.style.display = "flex";
          visibleCount++;
        } else {
          c.style.display = "none";
        }
      });

      var badge = document.getElementById("gallery-count-badge");
      if (badge) {
        badge.textContent = visibleCount + " shown";
      }
    }

    function filterGallery(category, btn) {
      document.querySelectorAll(".tab-btn").forEach(function(b) { b.classList.remove("active"); });
      if (btn) btn.classList.add("active");
      filterGalleryCombined();
    }

    function sortGalleryCards(order) {
      var grid = document.getElementById("main-gallery-grid");
      if (!grid) return;
      var cards = Array.from(grid.querySelectorAll(".gallery-card"));

      cards.sort(function(a, b) {
        if (order === "blur_asc") {
          var bA = parseFloat(a.getAttribute("data-blur") || "99999");
          var bB = parseFloat(b.getAttribute("data-blur") || "99999");
          if (bA < 0) bA = 99999;
          if (bB < 0) bB = 99999;
          return bA - bB;
        } else if (order === "blur_desc") {
          var bA = parseFloat(a.getAttribute("data-blur") || "-1");
          var bB = parseFloat(b.getAttribute("data-blur") || "-1");
          return bB - bA;
        } else if (order === "desc") {
          var nA = parseInt(a.getAttribute("data-objects") || "0", 10);
          var nB = parseInt(b.getAttribute("data-objects") || "0", 10);
          return nB - nA;
        } else if (order === "asc") {
          var nA = parseInt(a.getAttribute("data-objects") || "0", 10);
          var nB = parseInt(b.getAttribute("data-objects") || "0", 10);
          return nA - nB;
        } else if (order === "name") {
          var sA = a.getAttribute("data-stem") || "";
          var sB = b.getAttribute("data-stem") || "";
          return sA.localeCompare(sB);
        } else {
          var iA = parseInt(a.getAttribute("data-card-idx") || "0", 10);
          var iB = parseInt(b.getAttribute("data-card-idx") || "0", 10);
          return iA - iB;
        }
      });

      cards.forEach(function(c) { grid.appendChild(c); });
    }

    function inspectBlurryImages() {
      var gallerySec = document.getElementById("gallery");
      if (gallerySec) gallerySec.scrollIntoView({ behavior: "smooth" });
      var blurryTab = document.getElementById("tab-blurry-btn");
      if (blurryTab) filterGallery("blurry", blurryTab);
      var sortSelect = document.getElementById("gallery-sort-select");
      if (sortSelect) {
        sortSelect.value = "blur_asc";
        sortGalleryCards("blur_asc");
      }
    }

    function copyTrainingCode(btn) {
      var codeEl = document.getElementById("training-code-block");
      if (!codeEl) return;
      var txt = codeEl.textContent;
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(txt).then(function() {
          showCopied(btn);
        }).catch(function() {
          fallbackCopy(txt, btn);
        });
      } else {
        fallbackCopy(txt, btn);
      }
    }

    function fallbackCopy(text, btn) {
      var ta = document.createElement("textarea");
      ta.value = text;
      ta.style.position = "fixed";
      ta.style.left = "-9999px";
      document.body.appendChild(ta);
      ta.focus();
      ta.select();
      try {
        document.execCommand("copy");
        showCopied(btn);
      } catch (err) {}
      document.body.removeChild(ta);
    }

    function showCopied(btn) {
      var orig = btn.innerHTML;
      btn.innerHTML = "<span>\\u{2713} Copied!</span>";
      btn.style.background = "var(--color-success)";
      setTimeout(function() {
        btn.innerHTML = orig;
        btn.style.background = "";
      }, 2000);
    }

    function exportAuditJSON() {
      var auditData = window.__YOLO_AUDIT_DATA__ || {};
      var jsonStr = JSON.stringify(auditData, null, 2);
      var blob = new Blob([jsonStr], { type: "application/json" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url;
      a.download = "yolo_dataset_audit_" + (auditData.dataset_name || "report") + ".json";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }

    function setOverlayState(className, apply) {
      var grid = document.getElementById("main-gallery-grid");
      var preview = document.getElementById("lightboxPreview");
      if (grid) {
        if (apply) grid.classList.add(className);
        else grid.classList.remove(className);
      }
      if (preview) {
        if (apply) preview.classList.add(className);
        else preview.classList.remove(className);
      }

      var isChecked = !apply;
      var idMap = {
        "hide-bboxes": ["chk-bboxes", "lb-chk-bboxes"],
        "hide-labels": ["chk-labels", "lb-chk-labels"],
        "hide-fills": ["chk-fills", "lb-chk-fills"]
      };
      var ids = idMap[className] || [];
      ids.forEach(function(id) {
        var el = document.getElementById(id);
        if (el && el.checked !== isChecked) el.checked = isChecked;
      });

      var pillMap = {
        "hide-bboxes": "lb-pill-bboxes",
        "hide-labels": "lb-pill-labels",
        "hide-fills": "lb-pill-fills"
      };
      var pillId = pillMap[className];
      if (pillId) {
        var pEl = document.getElementById(pillId);
        if (pEl) pEl.classList.toggle("active", isChecked);
      }

      try {
        localStorage.setItem("pliman_report_" + className, apply ? "1" : "0");
      } catch(e) {}
    }

    function toggleOverlay(className, apply) {
      setOverlayState(className, apply);
    }

    var currentLightboxIdx = 0;
    var allCards = [];

    function initReport() {
      allCards = document.querySelectorAll(".gallery-card");
      try {
        var savedTheme = localStorage.getItem("pliman_report_theme");
        if (savedTheme && (savedTheme === "dark" || savedTheme === "light")) {
          document.documentElement.setAttribute("data-theme", savedTheme);
          document.getElementById("theme-icon").textContent = savedTheme === "dark" ? "\\u{1f319}" : "\\u{2600}\\u{fe0f}";
          document.getElementById("theme-text").textContent = savedTheme === "dark" ? "Dark" : "Light";
        }
      } catch(e) {}

      try {
        var savedLwd = localStorage.getItem("pliman_report_box_lwd");
        if (savedLwd) updateBoxLwd(savedLwd);
      } catch(e) {}

      // Restore overlay toggles from localStorage
      ["hide-bboxes", "hide-labels", "hide-fills"].forEach(function(cls) {
        try {
          var saved = localStorage.getItem("pliman_report_" + cls);
          if (saved === "1") {
            setOverlayState(cls, true);
          } else if (saved === "0") {
            setOverlayState(cls, false);
          }
        } catch(e) {}
      });

      // Count anomalies for gallery tab badge
      var anomCards = document.querySelectorAll(\'.gallery-card[data-has-anomaly="true"]\');
      var anomBadge = document.getElementById("gallery-anom-count");
      if (anomBadge) anomBadge.textContent = anomCards.length;

      // Count blurry images for gallery tab badge
      var blurCards = document.querySelectorAll(\'.gallery-card[data-is-blurry="true"]\');
      var blurBadge = document.getElementById("gallery-blurry-count");
      if (blurBadge) blurBadge.textContent = blurCards.length;
    }

    function openLightbox(cardIdx) {
      allCards = document.querySelectorAll(".gallery-card");
      currentLightboxIdx = cardIdx;
      renderLightbox();
      document.getElementById("lightboxModal").classList.add("open");
      document.addEventListener("keydown", handleLightboxKeys);
    }

    function openLightboxByStem(stem) {
      var card = document.querySelector(\'.gallery-card[data-stem="\' + stem + \'"]\');
      if (card) {
        var idx = parseInt(card.getAttribute("data-card-idx"), 10);
        closeAnomalyModal();
        openLightbox(idx);
      } else {
        alert("Image thumbnail was not sampled into the gallery. Image stem: " + stem);
      }
    }

    function closeLightbox() {
      document.getElementById("lightboxModal").classList.remove("open");
      document.removeEventListener("keydown", handleLightboxKeys);
    }

    function handleModalBackdropClick(e) {
      if (e.target.id === "lightboxModal") {
        closeLightbox();
      }
    }

    function prevLightboxImage(e) {
      if (e) e.stopPropagation();
      if (!allCards || allCards.length === 0) return;
      currentLightboxIdx = (currentLightboxIdx - 1 + allCards.length) % allCards.length;
      renderLightbox();
    }

    function nextLightboxImage(e) {
      if (e) e.stopPropagation();
      if (!allCards || allCards.length === 0) return;
      currentLightboxIdx = (currentLightboxIdx + 1) % allCards.length;
      renderLightbox();
    }

    function handleLightboxKeys(e) {
      if (e.key === "Escape") {
        closeLightbox();
        closeAnomalyModal();
        return;
      }
      if (e.target && (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA")) {
        if (e.target.type !== "checkbox" && e.target.type !== "range") return;
      }
      if (e.key === "ArrowLeft") {
        prevLightboxImage();
        return;
      }
      if (e.key === "ArrowRight") {
        nextLightboxImage();
        return;
      }
      if (e.key === "b" || e.key === "B") {
        var preview = document.getElementById("lightboxPreview");
        var isHidden = preview ? preview.classList.contains("hide-bboxes") : false;
        setOverlayState("hide-bboxes", !isHidden);
      } else if (e.key === "l" || e.key === "L") {
        var preview = document.getElementById("lightboxPreview");
        var isHidden = preview ? preview.classList.contains("hide-labels") : false;
        setOverlayState("hide-labels", !isHidden);
      } else if (e.key === "f" || e.key === "F") {
        var preview = document.getElementById("lightboxPreview");
        var isHidden = preview ? preview.classList.contains("hide-fills") : false;
        setOverlayState("hide-fills", !isHidden);
      }
    }

    function renderLightbox() {
      if (!allCards || allCards.length === 0) return;
      var card = allCards[currentLightboxIdx];
      if (!card) return;

      var img = card.querySelector(".card-thumb");
      var svg = card.querySelector(".yolo-svg-overlay");
      var hidden = card.querySelector(".hidden-card-data");

      var lbImg = document.getElementById("lightboxImg");
      var frame = document.getElementById("lightboxMediaFrame");
      lbImg.src = img.src;
      document.getElementById("lightboxSvgContainer").innerHTML = svg ? svg.outerHTML : "";

      // Sync lightbox overlay toggles with current preview classes
      var prevEl = document.getElementById("lightboxPreview");
      if (prevEl) {
        ["hide-bboxes", "hide-labels", "hide-fills"].forEach(function(cls) {
          var isHidden = prevEl.classList.contains(cls);
          var isChecked = !isHidden;
          var chkId = cls === "hide-bboxes" ? "lb-chk-bboxes" : (cls === "hide-labels" ? "lb-chk-labels" : "lb-chk-fills");
          var pillId = cls === "hide-bboxes" ? "lb-pill-bboxes" : (cls === "hide-labels" ? "lb-pill-labels" : "lb-pill-fills");
          var chk = document.getElementById(chkId);
          if (chk) chk.checked = isChecked;
          var pill = document.getElementById(pillId);
          if (pill) pill.classList.toggle("active", isChecked);
        });
      }

      var dimStr = hidden.getAttribute("data-dim") || "";
      if (dimStr.indexOf("\\u{d7}") !== -1) {
        var parts = dimStr.split("\\u{d7}");
        var w = parseFloat(parts[0]);
        var h = parseFloat(parts[1]);
        if (w > 0 && h > 0) {
          frame.style.aspectRatio = (w / h);
        } else {
          frame.style.aspectRatio = "auto";
        }
      } else {
        frame.style.aspectRatio = "auto";
      }

      document.getElementById("lightboxStem").textContent = hidden.getAttribute("data-stem");
      var catSlug = card.getAttribute("data-category");
      var catLabel = hidden.getAttribute("data-category");
      var catEl = document.getElementById("lightboxCategory");
      catEl.textContent = catLabel;
      catEl.className = "card-badge badge-" + catSlug;

      document.getElementById("lightboxSplit").textContent = hidden.getAttribute("data-split");
      document.getElementById("lightboxDim").textContent = hidden.getAttribute("data-dim");

      var blurEl = document.getElementById("lightboxBlur");
      if (blurEl) {
        var blurVal = parseFloat(card.getAttribute("data-blur") || "-1");
        var isBlur = card.getAttribute("data-is-blurry") === "true";
        if (blurVal >= 0) {
          blurEl.style.display = "inline-block";
          if (isBlur) {
            blurEl.className = "meta-tag meta-blur-warning";
            blurEl.textContent = "\\u{1f32b}\\u{fe0f} Blurry (Laplacian Var: " + blurVal.toFixed(1) + ")";
          } else {
            blurEl.className = "meta-tag meta-blur-sharp";
            blurEl.textContent = "\\u{26a1} Sharp (Laplacian Var: " + blurVal.toFixed(1) + ")";
          }
        } else {
          blurEl.style.display = "none";
        }
      }

      var lbAugContainer = document.getElementById("lightboxAugTags");
      if (lbAugContainer) {
        var augStr = (hidden.getAttribute("data-aug-tags") || "").trim();
        if (augStr) {
          var tags = augStr.split(/\\s+/);
          var augDefs = {
            blur: { label: "Blur", icon: "\\u{1f32b}\\u{fe0f}" },
            noise: { label: "Noise", icon: "\\u{26a1}" },
            rot: { label: "Rotate", icon: "\\u{1f504}" },
            persp: { label: "Perspective", icon: "\\u{1f4d0}" },
            fliph: { label: "Flip H", icon: "\\u{2194}\\u{fe0f}" },
            flipv: { label: "Flip V", icon: "\\u{2195}\\u{fe0f}" },
            bright: { label: "Brightness", icon: "\\u{2600}\\u{fe0f}" },
            contrast: { label: "Contrast", icon: "\\u{1f313}" },
            cutout: { label: "Cutout", icon: "\\u{2702}\\u{fe0f}" },
            mosaic: { label: "Mosaic", icon: "\\u{1f9e9}" },
            synth_bg: { label: "Synthetic BG", icon: "\\u{1f5bc}\\u{fe0f}" },
            tile: { label: "Tiled", icon: "\\u{1fa9f}" }
          };
          var augHtml = \'<span style="font-size:0.75rem; font-weight:700; color:var(--text-muted); text-transform:uppercase; margin-right:0.3rem; align-self:center;">Augmentations:</span>\';
          for (var k = 0; k < tags.length; k++) {
            var t = tags[k];
            var d = augDefs[t] || { label: t, icon: "\\u{1f3f7}\\u{fe0f}" };
            augHtml += \'<span class="card-aug-tag" style="font-size:0.75rem; padding:0.2rem 0.5rem;">\' + d.icon + \' \' + d.label + \'</span>\';
          }
          lbAugContainer.innerHTML = augHtml;
          lbAugContainer.style.display = "flex";
        } else {
          lbAugContainer.innerHTML = "";
          lbAugContainer.style.display = "none";
        }
      }

      // Check anomalies for this image
      var anomContainer = document.getElementById("lightboxAnomaliesBox");
      var anomData = [];
      try {
        anomData = JSON.parse(hidden.getAttribute("data-anomalies") || "[]");
      } catch(e) { anomData = []; }

      if (anomContainer) {
        if (anomData.length > 0) {
          var anomHtml = \'<div class="lightbox-anomaly-box"><strong>\\u{26a0}\\u{fe0f} Detected Integrity Issues:</strong><ul>\';
          for (var a = 0; a < anomData.length; a++) {
            anomHtml += \'<li><strong>\' + (anomData[a].type_label || anomData[a].type) + \':</strong> \' + anomData[a].detail + \'</li>\';
          }
          anomHtml += \'</ul></div>\';
          anomContainer.innerHTML = anomHtml;
          anomContainer.style.display = "block";
        } else {
          anomContainer.innerHTML = "";
          anomContainer.style.display = "none";
        }
      }

      var listContainer = document.getElementById("lightboxObjectsList");
      var objs = [];
      try {
        objs = JSON.parse(hidden.getAttribute("data-objects") || "[]");
      } catch (err) { objs = []; }

      if (objs.length === 0) {
        listContainer.innerHTML = \'<div style="color:var(--text-muted); font-size:0.85rem; padding:1rem 0;">No annotations (negative background sample).</div>\';
      } else {
        var html = \'<table class="report-table" style="font-size:0.8rem;"><thead><tr><th>#</th><th>Class</th><th>Type</th><th>Area %</th></tr></thead><tbody>\';
        for (var i = 0; i < objs.length; i++) {
          var o = objs[i];
          html += \'<tr><td>\' + (i + 1) + \'</td><td><strong style="color:\' + o.color + \'">\' + o.class + \'</strong></td><td>\' + o.type + \'</td><td>\' + o.area_pct + \'%</td></tr>\';
        }
        html += \'</tbody></table>\';
        listContainer.innerHTML = html;
      }
    }

    // Anomaly Inspector Modal Logic
    var currentAnomalyFilter = "all";

    function openAnomalyModal(filterType) {
      currentAnomalyFilter = filterType || "all";
      var modal = document.getElementById("anomalyModal");
      if (!modal) return;

      // Update active filter button in modal
      var tabs = document.querySelectorAll("#anomalyFilterTabs .tab-btn");
      tabs.forEach(function(b) {
        var t = b.getAttribute("data-type") || "all";
        if (t === currentAnomalyFilter) {
          b.classList.add("active");
        } else {
          b.classList.remove("active");
        }
      });

      renderAnomalyTable();
      modal.classList.add("open");
    }

    function closeAnomalyModal() {
      var modal = document.getElementById("anomalyModal");
      if (modal) modal.classList.remove("open");
    }

    function handleAnomalyBackdropClick(e) {
      if (e.target.id === "anomalyModal") {
        closeAnomalyModal();
      }
    }

    function filterAnomalyTable(filterType, btn) {
      currentAnomalyFilter = filterType;
      if (btn) {
        var tabs = document.querySelectorAll("#anomalyFilterTabs .tab-btn");
        tabs.forEach(function(b) { b.classList.remove("active"); });
        btn.classList.add("active");
      }
      renderAnomalyTable();
    }

    function renderAnomalyTable() {
      var tbody = document.getElementById("anomalyTableBody");
      var searchInput = document.getElementById("anomalySearchInput");
      var searchVal = searchInput ? searchInput.value.toLowerCase().trim() : "";
      var allAnoms = window.__YOLO_ANOMALIES__ || [];

      var rows = [];
      for (var i = 0; i < allAnoms.length; i++) {
        var a = allAnoms[i];
        var matchType = (currentAnomalyFilter === "all" || a.type === currentAnomalyFilter);
        var matchSearch = (!searchVal || (a.stem || "").toLowerCase().indexOf(searchVal) !== -1);

        if (matchType && matchSearch) {
          var thumbHtml = "";
          if (a.thumb_uri) {
            thumbHtml = \'<img src="\' + a.thumb_uri + \'" class="anomaly-thumb" onclick="openLightboxByStem(\\\'\' + a.stem + \'\\\')" title="Click to view in Lightbox" />\';
          } else {
            thumbHtml = \'<div class="anomaly-thumb-placeholder" onclick="openLightboxByStem(\\\'\' + a.stem + \'\\\')" title="Click to view">\\u{1f5bc}\\u{fe0f}</div>\';
          }

          var splitBadge = \'<span class="card-badge badge-split">\' + (a.split || "train").toUpperCase() + \'</span>\';
          var typeBadge = \'<span class="tag-badge \' + (a.badge_class || "pill-warn") + \'">\' + (a.type_label || a.type) + \'</span>\';

          var tr = \'<tr>\' +
            \'<td>\' + thumbHtml + \'</td>\' +
            \'<td><strong class="anomaly-stem-code">\' + a.stem + \'</strong></td>\' +
            \'<td>\' + splitBadge + \'</td>\' +
            \'<td>\' + typeBadge + \'</td>\' +
            \'<td><div class="anomaly-detail-text">\' + a.detail + \'</div></td>\' +
            \'<td style="text-align: right;"><button type="button" class="btn-action" onclick="openLightboxByStem(\\\'\' + a.stem + \'\\\')" style="font-size:0.75rem; padding:0.25rem 0.6rem;">\\u{1f50d} Lightbox</button></td>\' +
            \'</tr>\';
          rows.push(tr);
        }
      }

      if (rows.length === 0) {
        tbody.innerHTML = \'<tr><td colspan="6" style="text-align:center; padding:2rem; color:var(--text-muted);">No anomalies found matching current filter.</td></tr>\';
      } else {
        tbody.innerHTML = rows.join("");
      }

      var countEl = document.getElementById("anomalyModalCount");
      if (countEl) countEl.textContent = "(" + rows.length + " instances)";
    }

    document.addEventListener("DOMContentLoaded", initReport);
  '
}


# ==============================================================================
# INTERNAL HELPER: ASSEMBLE COMPLETE MODERN HTML REPORT
# ==============================================================================
.yolo_generate_report_html <- function(desc, cards, health, title, theme = "dark", lwd = 1.5, blur_threshold = 80) {
  tot_imgs <- format(desc$total_images %||% 0L, big.mark = ",")
  tot_objs <- format(desc$total_objects %||% 0L, big.mark = ",")
  disk_mb <- desc$disk_mb %||% 0
  avg_kb <- desc$avg_file_kb %||% 0
  tot_orig <- desc$total_original %||% 0L
  tot_aug <- desc$total_augmented %||% 0L
  tot_synth <- desc$total_synthetic %||% 0L
  tot_tile <- desc$total_tiled %||% 0L
  pct_orig <- if ((desc$total_images %||% 0) > 0) round(tot_orig / desc$total_images * 100, 1) else 0

  task_badge <- switch(desc$task,
    detect = "Object Detection (BBoxes)",
    segment = "Instance Segmentation (Polygons)",
    obb = "Oriented Bounding Box (OBB)",
    classify = "Classification",
    toupper(desc$task)
  )

  task_badge_class <- switch(desc$task,
    detect = "badge-detect",
    segment = "badge-segment",
    obb = "badge-obb",
    classify = "badge-classify",
    "badge-detect"
  )

  imb_val <- desc$imbalance_ratio %||% 1.0
  imb_label <- if (imb_val <= 2.0) "Balanced" else if (imb_val <= 4.0) "Moderate" else "High Skew"

  # Circular Health Score Gauge calculations
  h_score <- health$score %||% 85L
  h_grade <- health$grade %||% "A"
  h_grade_cls <- health$grade_class %||% "grade-opt"

  # Circumference for r=32: C = 2 * pi * 32 = 201.06
  circ <- 201.06
  dash_offset <- round(circ * (1 - (h_score / 100)), 2)

  # Executive KPI Cards HTML
  kpi_html <- paste0(
    '<div class="kpi-grid">\n',
    '  <div class="kpi-card card-glow kpi-gauge-card">\n',
    '    <div class="gauge-wrapper">\n',
    '      <svg class="gauge-svg" viewBox="0 0 80 80">\n',
    '        <circle class="gauge-bg" cx="40" cy="40" r="32"></circle>\n',
    '        <circle class="gauge-fill ', h_grade_cls, '" cx="40" cy="40" r="32" stroke-dasharray="', circ, '" stroke-dashoffset="', dash_offset, '"></circle>\n',
    '      </svg>\n',
    '      <div class="gauge-center">\n',
    '        <div class="gauge-score">', h_score, '</div>\n',
    '        <div class="gauge-grade ', h_grade_cls, '">', h_grade, '</div>\n',
    '      </div>\n',
    '    </div>\n',
    '    <div>\n',
    '      <div class="kpi-header"><span class="kpi-title">Health Score</span><span class="kpi-icon">\U0001fa7a</span></div>\n',
    '      <div style="font-weight: 700; font-size: 1.05rem; margin-bottom: 0.2rem;">Grade ', h_grade, ' (', h_score, '/100)</div>\n',
    '      <div class="kpi-subtitle">', health$n_pass, ' pass \u2022 ', health$n_info, ' info \u2022 ', health$n_warn, ' warn</div>\n',
    '    </div>\n',
    '  </div>\n',
    '  <div class="kpi-card card-glow">\n',
    '    <div class="kpi-header"><span class="kpi-title">Total Images</span><span class="kpi-icon">\U0001f5bc\ufe0f</span></div>\n',
    '    <div class="kpi-value">', tot_imgs, '</div>\n',
    '    <div class="kpi-subtitle">', disk_mb, ' MB on disk <span class="dot-sep">\u2022</span> ~', avg_kb, ' KB/img</div>\n',
    '  </div>\n',
    '  <div class="kpi-card card-glow">\n',
    '    <div class="kpi-header"><span class="kpi-title">Total Objects</span><span class="kpi-icon">\U0001f3af</span></div>\n',
    '    <div class="kpi-value">', tot_objs, '</div>\n',
    '    <div class="kpi-subtitle">Mean ', desc$mean_density %||% 0, ' obj/img <span class="dot-sep">\u2022</span> Max ', desc$max_density %||% 0, '</div>\n',
    '  </div>\n',
    '  <div class="kpi-card card-glow">\n',
    '    <div class="kpi-header"><span class="kpi-title">Original vs Generated</span><span class="kpi-icon">\U0001f9ec</span></div>\n',
    '    <div class="kpi-value">', pct_orig, '% <span class="kpi-val-sub">orig</span></div>\n',
    '    <div class="kpi-subtitle">', tot_orig, ' orig <span class="dot-sep">\u2022</span> ', tot_aug, ' aug <span class="dot-sep">\u2022</span> ', tot_synth, ' synth</div>\n',
    '  </div>\n',
    '  <div class="kpi-card card-glow">\n',
    '    <div class="kpi-header"><span class="kpi-title">Classes & Balance</span><span class="kpi-icon">\u2696\ufe0f</span></div>\n',
    '    <div class="kpi-value">', length(desc$classes), ' <span class="kpi-val-sub">classes</span></div>\n',
    '    <div class="kpi-subtitle">Imbalance: ', sprintf("%.1fx", imb_val), ' (', imb_label, ')</div>\n',
    '  </div>\n',
    '</div>\n'
  )

  # Anomaly strip
  anomaly_strip_html <- .yolo_build_anomaly_strip_html(desc)

  # Provenance visual progress bar
  p_tot <- max(1L, desc$total_images %||% 1L)
  pct_p_orig <- round(tot_orig / p_tot * 100, 1)
  pct_p_aug <- round(tot_aug / p_tot * 100, 1)
  pct_p_synth <- round(tot_synth / p_tot * 100, 1)
  pct_p_tile <- round(tot_tile / p_tot * 100, 1)

  provenance_bar_html <- paste0(
    '<div class="prov-bar-container">\n',
    '  <div class="prov-bar">\n',
    '    <div class="prov-seg prov-orig" style="width: ', pct_p_orig, '%" title="Original: ', tot_orig, ' (', pct_p_orig, '%)"></div>\n',
    '    <div class="prov-seg prov-aug" style="width: ', pct_p_aug, '%" title="Augmented: ', tot_aug, ' (', pct_p_aug, '%)"></div>\n',
    '    <div class="prov-seg prov-synth" style="width: ', pct_p_synth, '%" title="Synthetic BG: ', tot_synth, ' (', pct_p_synth, '%)"></div>\n',
    '    <div class="prov-seg prov-tile" style="width: ', pct_p_tile, '%" title="Tiled: ', tot_tile, ' (', pct_p_tile, '%)"></div>\n',
    '  </div>\n',
    '  <div class="prov-legend">\n',
    '    <span class="legend-item"><span class="legend-color prov-orig-bg"></span> Original (', tot_orig, ' \u2022 ', pct_p_orig, '%)</span>\n',
    '    <span class="legend-item"><span class="legend-color prov-aug-bg"></span> Augmented (', tot_aug, ' \u2022 ', pct_p_aug, '%)</span>\n',
    '    <span class="legend-item"><span class="legend-color prov-synth-bg"></span> Synthetic BG (', tot_synth, ' \u2022 ', pct_p_synth, '%)</span>\n',
    '    <span class="legend-item"><span class="legend-color prov-tile-bg"></span> Tiled (', tot_tile, ' \u2022 ', pct_p_tile, '%)</span>\n',
    '  </div>\n',
    '</div>\n'
  )

  # Splits Table HTML
  splits_table_html <- ""
  if (!is.null(desc$splits_summary) && nrow(desc$splits_summary) > 0L) {
    sp_df <- desc$splits_summary
    rows <- character(0)
    for (i in seq_len(nrow(sp_df))) {
      r <- sp_df[i, ]
      is_tot <- r$Split == "TOTAL"
      row_cls <- if (is_tot) ' class="row-total"' else ""
      bg_val <- if (!is.null(r$Background)) sprintf("%d (%.1f%%)", r$Background, r$Pct_Bg) else "-"
      objs_val <- if (!is.null(r$Objects)) format(r$Objects, big.mark = ",") else "-"
      mean_val <- if (!is.null(r$Mean_Objs)) sprintf("%.2f", r$Mean_Objs) else "-"
      max_val <- if (!is.null(r$Max_Objs)) as.character(r$Max_Objs) else "-"

      rows <- c(rows, sprintf(
        '<tr%s><td><strong>%s</strong></td><td>%s</td><td>%.1f%%</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
        row_cls,
        r$Split,
        format(r$Images, big.mark = ","),
        r$Pct_Images,
        bg_val,
        objs_val,
        mean_val,
        max_val
      ))
    }
    splits_table_html <- paste0(
      '<div class="table-responsive">\n',
      '  <table class="report-table">\n',
      '    <thead><tr><th>Split</th><th>Images</th><th>% of Total</th><th>Negative / BG</th><th>Total Objects</th><th>Mean Objs/Img</th><th>Max Objs</th></tr></thead>\n',
      '    <tbody>', paste(rows, collapse = "\n"), '</tbody>\n',
      '  </table>\n',
      '</div>\n'
    )
  }

  # Classes Table HTML
  classes_table_html <- ""
  if (!is.null(desc$classes_summary) && nrow(desc$classes_summary) > 0L) {
    cls_df <- desc$classes_summary
    splits_in_cls <- setdiff(colnames(cls_df), c("ID", "Class", "Total_Objs", "Pct_Objs", "Images", "Pct_Images", "Avg_Area_Pct", "Size_Cat", "Total", "Pct_Total"))
    sp_th <- paste(sprintf("<th>%s</th>", toupper(splits_in_cls)), collapse = "")

    rows <- character(0)
    palette <- c("#6366f1", "#10b981", "#f59e0b", "#ef4444", "#8b5cf6", "#ec4899", "#06b6d4", "#f97316", "#14b8a6", "#84cc16")

    for (i in seq_len(nrow(cls_df))) {
      r <- cls_df[i, ]
      cid <- r$ID
      cname <- r$Class
      c_col <- palette[(cid %% length(palette)) + 1L]
      tot_val <- r$Total_Objs %||% r$Total %||% 0L
      pct_val <- r$Pct_Objs %||% r$Pct_Total %||% 0.0

      split_tds <- ""
      if (length(splits_in_cls) > 0L) {
        split_tds <- paste(vapply(splits_in_cls, function(s) sprintf("<td>%s</td>", format(r[[s]], big.mark = ",")), character(1)), collapse = "")
      }

      size_badge <- if (!is.null(r$Size_Cat) && nzchar(r$Size_Cat)) {
        sprintf('<span class="tag-badge size-%s">%s</span>', tolower(r$Size_Cat), r$Size_Cat)
      } else "-"
      area_val <- if (!is.null(r$Avg_Area_Pct)) sprintf("%.2f%%", r$Avg_Area_Pct) else "-"

      rows <- c(rows, sprintf(
        '<tr>
          <td><span class="class-id-badge" style="background-color: %s">%d</span></td>
          <td><strong>%s</strong></td>
          %s
          <td><strong>%s</strong></td>
          <td>
            <div class="table-bar-cell">
              <div class="mini-bar-bg"><div class="mini-bar-fill" style="width: %.1f%%; background-color: %s"></div></div>
              <span>%.1f%%</span>
            </div>
          </td>
          <td>%s</td>
          <td>%s</td>
        </tr>',
        c_col, cid, cname, split_tds, format(tot_val, big.mark = ","), pct_val, c_col, pct_val, area_val, size_badge
      ))
    }

    classes_table_html <- paste0(
      '<div class="table-responsive">\n',
      '  <table class="report-table">\n',
      '    <thead><tr><th>ID</th><th>Class Name</th>', sp_th, '<th>Total Objects</th><th>Share</th><th>Avg Image Area</th><th>Scale Category</th></tr></thead>\n',
      '    <tbody>', paste(rows, collapse = "\n"), '</tbody>\n',
      '  </table>\n',
      '</div>\n'
    )
  }

  # Provenance Table HTML
  prov_table_html <- ""
  if (!is.null(desc$provenance_summary) && nrow(desc$provenance_summary) > 0L) {
    p_df <- desc$provenance_summary
    splits_in_p <- setdiff(colnames(p_df), c("Type", "Total", "Pct"))
    sp_th <- paste(sprintf("<th>%s</th>", toupper(splits_in_p)), collapse = "")

    rows <- character(0)
    for (i in seq_len(nrow(p_df))) {
      r <- p_df[i, ]
      is_tot <- r$Type == "TOTAL"
      row_cls <- if (is_tot) ' class="row-total"' else ""
      split_tds <- ""
      if (length(splits_in_p) > 0L) {
        split_tds <- paste(vapply(splits_in_p, function(s) sprintf("<td>%s</td>", format(r[[s]], big.mark = ",")), character(1)), collapse = "")
      }

      badge_cls <- switch(tolower(r$Type),
        original = "prov-badge-orig",
        augmented = "prov-badge-aug",
        synthetic_bg = "prov-badge-synth",
        tiled = "prov-badge-tile",
        "prov-badge-tot"
      )

      rows <- c(rows, sprintf(
        '<tr%s><td><span class="tag-badge %s">%s</span></td>%s<td><strong>%s</strong></td><td>%.1f%%</td></tr>',
        row_cls, badge_cls, r$Type, split_tds, format(r$Total, big.mark = ","), r$Pct
      ))
    }

    prov_table_html <- paste0(
      '<div class="table-responsive">\n',
      '  <table class="report-table">\n',
      '    <thead><tr><th>Provenance Type</th>', sp_th, '<th>Total Images</th><th>% of Dataset</th></tr></thead>\n',
      '    <tbody>', paste(rows, collapse = "\n"), '</tbody>\n',
      '  </table>\n',
      '</div>\n'
    )
  }

  # Spatial & Dimension Heatmaps
  spatial_charts_html <- .yolo_build_spatial_charts_html(desc)

  # Co-occurrence Matrix
  cooccurrence_html <- .yolo_build_cooccurrence_html(desc)

  # Training Advisor
  advisor_html <- .yolo_build_training_advisor(desc)

  # Density & Scale Breakdown HTML
  density_cards_html <- ""
  if (!is.null(desc$density_summary) || !is.null(desc$scale_summary)) {
    d_rows <- character(0)
    if (!is.null(desc$density_summary)) {
      for (i in seq_len(nrow(desc$density_summary))) {
        r <- desc$density_summary[i, ]
        d_rows <- c(d_rows, sprintf(
          '<tr><td><strong>%s</strong></td><td>%s</td><td><div class="table-bar-cell"><div class="mini-bar-bg"><div class="mini-bar-fill" style="width: %.1f%%; background-color: #0ea5e9"></div></div><span>%.1f%%</span></div></td></tr>',
          r$Bucket, format(r$Images, big.mark = ","), r$Pct_Images, r$Pct_Images
        ))
      }
    }

    s_rows <- character(0)
    if (!is.null(desc$scale_summary)) {
      for (i in seq_len(nrow(desc$scale_summary))) {
        r <- desc$scale_summary[i, ]
        s_rows <- c(s_rows, sprintf(
          '<tr><td><strong>%s</strong></td><td>%s</td><td><div class="table-bar-cell"><div class="mini-bar-bg"><div class="mini-bar-fill" style="width: %.1f%%; background-color: #a855f7"></div></div><span>%.1f%%</span></div></td></tr>',
          r$Scale, format(r$Objects, big.mark = ","), r$Pct_Objects, r$Pct_Objects
        ))
      }
    }

    density_cards_html <- paste0(
      '<div class="grid-2col">\n',
      '  <div class="content-card">\n',
      '    <h3 class="card-subtitle">Object Density Distribution (Objects per Image)</h3>\n',
      '    <table class="report-table"><thead><tr><th>Density Bracket</th><th>Images</th><th>Share</th></tr></thead><tbody>\n',
      paste(d_rows, collapse = "\n"),
      '    </tbody></table>\n',
      '  </div>\n',
      '  <div class="content-card">\n',
      '    <h3 class="card-subtitle">Spatial Object Scale (COCO Benchmark Standard)</h3>\n',
      '    <table class="report-table"><thead><tr><th>Scale Bracket</th><th>Objects</th><th>Share</th></tr></thead><tbody>\n',
      paste(s_rows, collapse = "\n"),
      '    </tbody></table>\n',
      '  </div>\n',
      '</div>\n'
    )
  }

  # Compute Augmentation Tag Frequencies across sampled cards
  aug_tag_counts <- list()
  for (c in cards) {
    tags_in_c <- c$aug_tags %||% character(0)
    if (length(tags_in_c) > 0L) {
      for (t in tags_in_c) {
        aug_tag_counts[[t]] <- (aug_tag_counts[[t]] %||% 0L) + 1L
      }
    }
  }

  aug_pills_html <- character(0)
  for (t in names(.yolo_aug_defs)) {
    cnt <- aug_tag_counts[[t]] %||% 0L
    if (cnt > 0L) {
      def <- .yolo_aug_defs[[t]]
      aug_pills_html <- c(aug_pills_html, sprintf(
        '<button type="button" class="aug-pill" data-tag="%s" onclick="toggleAugFilter(\'%s\', this)" title="Filter images containing %s augmentation">%s %s <span class="pill-count">%d</span></button>',
        t, t, def$label, def$icon, def$label, cnt
      ))
    }
  }

  gallery_aug_toolbar_html <- if (length(aug_pills_html) > 0L) {
    paste0(
      '        <div class="gallery-toolbar-aug" id="gallery-aug-toolbar">\n',
      '          <div class="aug-toolbar-header">\n',
      '            <span class="aug-toolbar-label">\U0001f9ea Augmentation Filters:</span>\n',
      '            <span class="aug-toolbar-desc">(Click pills to combine filters with AND, e.g. Blur + Noise)</span>\n',
      '            <button type="button" class="aug-clear-btn" id="aug-clear-btn" onclick="clearAugFilters()" style="display:none;">\u2715 Clear Augmentation Filters</button>\n',
      '          </div>\n',
      '          <div class="aug-pills-container">\n',
      '            ', paste(aug_pills_html, collapse = "\n            "), '\n',
      '          </div>\n',
      '        </div>\n'
    )
  } else ""

  # Build Gallery Cards HTML
  gallery_cards_html <- character(0)
  card_thumb_map <- list()

  for (idx in seq_along(cards)) {
    c <- cards[[idx]]
    badge_cls <- paste0("badge-", c$category_slug)
    card_thumb_map[[c$stem]] <- list(uri = c$data_uri, idx = idx - 1L)

    objs_json <- jsonlite::toJSON(c$objects, auto_unbox = TRUE)
    objs_json_safe <- gsub("'", "&#39;", gsub('"', "&quot;", objs_json))

    anom_json <- jsonlite::toJSON(c$anomalies, auto_unbox = TRUE)
    anom_json_safe <- gsub("'", "&#39;", gsub('"', "&quot;", anom_json))

    anomaly_badge_html <- if (isTRUE(c$has_anomaly)) {
      sprintf('<span class="card-badge badge-anomaly" title="%s">\u26a0\ufe0f ANOMALY</span>', c$anomaly_summary)
    } else ""

    blur_score_val <- if (!is.null(c$blur_score) && !is.na(c$blur_score)) as.numeric(c$blur_score) else -1
    is_blurry_val <- isTRUE(c$is_blurry)

    blur_badge_html <- if (is_blurry_val) {
      sprintf('<span class="card-badge badge-blur" title="Blur score (Laplacian var: %.1f < %.1f)">\U0001f32b\ufe0f BLUR</span>', blur_score_val, blur_threshold)
    } else ""

    blur_meta_html <- if (blur_score_val >= 0) {
      if (is_blurry_val) {
        sprintf('<span class="meta-tag meta-blur-warning" title="Laplacian Var: %.1f (Low Sharpness)">\U0001f32b\ufe0f %.0f</span>', blur_score_val, blur_score_val)
      } else {
        sprintf('<span class="meta-tag meta-blur-sharp" title="Laplacian Var: %.1f (Sharp)">\u26a1 %.0f</span>', blur_score_val, blur_score_val)
      }
    } else ""

    aug_tags_list <- c$aug_tags %||% character(0)
    aug_tags_str <- paste(aug_tags_list, collapse = " ")

    aug_chips_html <- if (length(aug_tags_list) > 0L) {
      chips <- vapply(aug_tags_list, function(t) {
        def <- .yolo_aug_defs[[t]]
        lbl <- if (!is.null(def)) paste(def$icon, def$label) else t
        sprintf('<span class="card-aug-tag" title="Augmentation: %s">%s</span>', t, lbl)
      }, character(1))
      paste(chips, collapse = "")
    } else ""

    gallery_cards_html <- c(gallery_cards_html, sprintf('
    <div class="gallery-card" data-category="%s" data-split="%s" data-stem="%s" data-objects="%d" data-card-idx="%d" data-has-anomaly="%s" data-blur="%.1f" data-is-blurry="%s" data-aug-tags="%s">
      <div class="card-media-wrapper" onclick="openLightbox(%d)">
        <div class="card-media-frame">
          <img class="card-thumb" src="%s" alt="%s" loading="lazy" />
          %s
        </div>
        <div class="card-badges">
          <span class="card-badge %s">%s</span>
          <span class="card-badge badge-split">%s</span>
          %s
          %s
        </div>
      </div>
      <div class="card-footer">
        <div class="card-stem" title="%s">%s</div>
        <div class="card-meta">
          <span class="meta-tag">%d\u00d7%d</span>
          <span class="meta-tag obj-tag">%d obj%s</span>
          %s
          %s
        </div>
      </div>
      <div class="hidden-card-data" style="display:none;" data-stem="%s" data-split="%s" data-category="%s" data-dim="%d\u00d7%d" data-objects="%s" data-anomalies="%s" data-blur="%.1f" data-is-blurry="%s" data-aug-tags="%s"></div>
    </div>
    ',
      c$category_slug,
      c$split,
      c$stem,
      c$n_objects,
      idx - 1L,
      if (isTRUE(c$has_anomaly)) "true" else "false",
      blur_score_val,
      if (is_blurry_val) "true" else "false",
      aug_tags_str,
      idx - 1L,
      c$data_uri,
      c$stem,
      c$svg_overlay,
      badge_cls,
      c$category_label,
      toupper(c$split),
      anomaly_badge_html,
      blur_badge_html,
      c$stem,
      c$stem,
      c$width, c$height,
      c$n_objects, if (c$n_objects == 1L) "" else "s",
      blur_meta_html,
      aug_chips_html,
      c$stem,
      toupper(c$split),
      c$category_label,
      c$width, c$height,
      objs_json_safe,
      anom_json_safe,
      blur_score_val,
      if (is_blurry_val) "true" else "false",
      aug_tags_str
    ))
  }

  # Health Checklist HTML
  health_items_html <- character(0)
  for (h in health$items) {
    status_cls <- paste0("health-", h$status)
    badge_cls <- paste0("badge-health-", h$status)
    icon <- switch(h$status, pass = "\u2713", warning = "!", info = "i", "\u2022")
    act_html <- if (!is.null(h$action_html) && nzchar(h$action_html)) h$action_html else ""

    health_items_html <- c(health_items_html, sprintf('
    <div class="health-card %s">
      <div class="health-icon">%s</div>
      <div class="health-body">
        <div class="health-head">
          <span class="health-title">%s</span>
          <span class="health-badge %s">%s</span>
        </div>
        <div class="health-metric">%s</div>
        <div class="health-detail">%s</div>
        %s
      </div>
    </div>
    ', status_cls, icon, h$title, badge_cls, h$badge, h$metric, h$detail, act_html))
  }

  # Build augmented anomalies JSON with thumbnail URIs
  anom_raw <- desc$anomalies$records %||% list()
  anom_records_for_ui <- list()
  for (a in anom_raw) {
    tm <- card_thumb_map[[a$stem]]
    anom_records_for_ui[[length(anom_records_for_ui) + 1L]] <- list(
      stem = a$stem,
      split = a$split,
      type = a$type,
      type_label = a$type_label %||% a$type,
      badge_class = a$badge_class %||% "pill-warn",
      detail = a$detail,
      thumb_uri = if (!is.null(tm)) tm$uri else "",
      card_idx = if (!is.null(tm)) tm$idx else -1L
    )
  }

  anomalies_json_str <- jsonlite::toJSON(anom_records_for_ui, auto_unbox = TRUE)

  generation_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  pliman_ver <- tryCatch(as.character(utils::packageVersion("pliman")), error = function(e) "devel")

  # Prepare JSON Audit Payload for direct client-side export
  audit_export_obj <- list(
    generator = paste0("pliman v", pliman_ver),
    generated_at = generation_timestamp,
    dataset_name = basename(desc$dir %||% "yolo_dataset"),
    dataset_directory = desc$dir %||% "",
    task = desc$task %||% "detect",
    health_score = health$score,
    health_grade = health$grade,
    summary = list(
      total_images = desc$total_images %||% 0L,
      total_objects = desc$total_objects %||% 0L,
      total_background = desc$total_background %||% 0L,
      pct_background = desc$pct_background %||% 0.0,
      mean_density = desc$mean_density %||% 0.0,
      max_density = desc$max_density %||% 0L,
      classes = desc$classes %||% character(0),
      imbalance_ratio = desc$imbalance_ratio %||% 1.0,
      disk_mb = desc$disk_mb %||% 0.0
    ),
    anomalies = desc$anomalies %||% list(),
    health_checks = health$items
  )
  audit_json_str <- jsonlite::toJSON(audit_export_obj, auto_unbox = TRUE, pretty = TRUE)

  # Full HTML Document Assembly
  paste0(
    '<!DOCTYPE html>\n',
    '<html lang="en" data-theme="', theme, '">\n',
    '<head>\n',
    '  <meta charset="UTF-8">\n',
    '  <meta name="viewport" content="width=device-width, initial-scale=1.0">\n',
    '  <title>', title, ' - pliman</title>\n',
    '  <link rel="preconnect" href="https://fonts.googleapis.com">\n',
    '  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n',
    '  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=Outfit:wght@400;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">\n',
    '  <style>', .yolo_report_css(lwd = lwd), '</style>\n',
    '</head>\n',
    '<body>\n',
    '  <nav class="navbar">\n',
    '    <div class="nav-brand">\n',
    '      <span class="brand-logo">pliman</span>\n',
    '      <span>YOLO Dataset Audit</span>\n',
    '    </div>\n',
    '    <ul class="nav-menu">\n',
    '      <li><a class="nav-link" href="#overview">Overview</a></li>\n',
    '      <li><a class="nav-link" href="#provenance">Provenance</a></li>\n',
    '      <li><a class="nav-link" href="#splits">Splits & Classes</a></li>\n',
    '      <li><a class="nav-link" href="#geometry">Geometry & Heatmaps</a></li>\n',
    if (length(desc$classes) > 1L) '      <li><a class="nav-link" href="#cooccurrence">Co-occurrence</a></li>\n' else '',
    '      <li><a class="nav-link" href="#gallery">Gallery</a></li>\n',
    '      <li><a class="nav-link" href="#advisor">AI Advisor</a></li>\n',
    '      <li><a class="nav-link" href="#health">Health Audit</a></li>\n',
    '    </ul>\n',
    '    <div class="nav-actions">\n',
    '      <button class="btn-action" onclick="openAnomalyModal(\'all\')" title="Inspect Dataset Anomalies">\n',
    '        <span>\U0001f50d</span> Anomaly Inspector\n',
    '      </button>\n',
    '      <button class="btn-action" onclick="exportAuditJSON()" title="Export JSON Audit Report">\n',
    '        <span>\U0001f4e5</span> JSON Audit\n',
    '      </button>\n',
    '      <button class="btn-action" onclick="toggleTheme()" title="Toggle Theme">\n',
    '        <span id="theme-icon">', if (theme == "dark") "\U0001f319" else "\u2600\ufe0f", '</span> <span id="theme-text">', if (theme == "dark") "Dark" else "Light", '</span>\n',
    '      </button>\n',
    '      <button class="btn-action" onclick="window.print()" title="Print or Export to PDF">\n',
    '        <span>\U0001f5a8\ufe0f</span> Print\n',
    '      </button>\n',
    '    </div>\n',
    '  </nav>\n\n',
    '  <div class="container">\n',
    '    <header class="header-banner" id="overview">\n',
    '      <h1 class="dataset-title">', title, '</h1>\n',
    '      <div class="dataset-meta">\n',
    '        <span class="badge-task ', task_badge_class, '">', task_badge, '</span>\n',
    '        <span>\U0001f4c1 <code>', desc$dir %||% "Dataset", '</code></span>\n',
    '        <span>\U0001f552 Generated on ', generation_timestamp, '</span>\n',
    '      </div>\n',
    '    </header>\n\n',
    '    ', kpi_html, '\n',
    '    ', anomaly_strip_html, '\n\n',
    '    <section class="section-block" id="provenance">\n',
    '      <div class="section-header">\n',
    '        <div>\n',
    '          <h2 class="section-title">\U0001f9ec Dataset Provenance & Augmentation</h2>\n',
    '          <p class="section-desc">Composition of original ground-truth images, offline augmentations, synthetic backgrounds, and tiled crops.</p>\n',
    '        </div>\n',
    '      </div>\n',
    '      ', provenance_bar_html, '\n',
    '      ', prov_table_html, '\n',
    '    </section>\n\n',
    '    <section class="section-block" id="splits">\n',
    '      <div class="section-header">\n',
    '        <div>\n',
    '          <h2 class="section-title">\u2696\ufe0f Splits Allocation & Class Balance</h2>\n',
    '          <p class="section-desc">Sample distribution across training, validation, and test partitions alongside per-class frequency.</p>\n',
    '        </div>\n',
    '      </div>\n',
    '      <div style="margin-bottom: 2rem;">\n',
    '        <h3 style="font-size: 1.05rem; font-weight: 700; margin-bottom: 0.85rem; color: var(--text-main);">Partition Allocations</h3>\n',
    '        ', splits_table_html, '\n',
    '      </div>\n',
    '      <div>\n',
    '        <h3 style="font-size: 1.05rem; font-weight: 700; margin-bottom: 0.85rem; color: var(--text-main);">Class Frequency & Geometry</h3>\n',
    '        ', classes_table_html, '\n',
    '      </div>\n',
    '    </section>\n\n',
    '    <section class="section-block" id="geometry">\n',
    '      <div class="section-header">\n',
    '        <div>\n',
    '          <h2 class="section-title">\U0001f5fa\ufe0f Spatial Distribution & Bounding Box Geometry</h2>\n',
    '          <p class="section-desc">Two-dimensional spatial density profiling across normalized image coordinates alongside aspect ratio classifications.</p>\n',
    '        </div>\n',
    '      </div>\n',
    '      ', spatial_charts_html, '\n',
    '      ', density_cards_html, '\n',
    '    </section>\n\n',
    cooccurrence_html,
    '    <section class="section-block" id="gallery">\n',
    '      <div class="section-header">\n',
    '        <div>\n',
    '          <h2 class="section-title">\U0001f5bc\ufe0f Interactive Annotated Image Gallery</h2>\n',
    '          <p class="section-desc">Pixel-accurate SVG overlays for sample images across each provenance category and splits (train, val, test) with live search and anomaly filtering.</p>\n',
    '        </div>\n',
    '      </div>\n',
    '      <div class="gallery-toolbar">\n',
    '        <div class="gallery-toolbar-top">\n',
    '          <div class="gallery-tabs">\n',
    '            <button class="tab-btn active" data-category="all" onclick="filterGallery(\'all\', this)">All Samples (', length(cards), ')</button>\n',
    '            <button class="tab-btn" data-category="original" onclick="filterGallery(\'original\', this)">Original</button>\n',
    '            <button class="tab-btn" data-category="augmented" onclick="filterGallery(\'augmented\', this)">Augmented</button>\n',
    '            <button class="tab-btn" data-category="synthetic_bg" onclick="filterGallery(\'synthetic_bg\', this)">Synthetic BG</button>\n',
    '            <button class="tab-btn" data-category="tiled" onclick="filterGallery(\'tiled\', this)">Tiled</button>\n',
    '            <button class="tab-btn" data-category="background" onclick="filterGallery(\'background\', this)">Negative / BG</button>\n',
    '            <button class="tab-btn" data-category="blurry" id="tab-blurry-btn" onclick="filterGallery(\'blurry\', this)" style="border-color: rgba(244, 63, 94, 0.4); color: #f43f5e;">\U0001f32b\ufe0f Blurry (<span id="gallery-blurry-count">0</span>)</button>\n',
    '            <button class="tab-btn" data-category="anomaly" onclick="filterGallery(\'anomaly\', this)" style="border-color: rgba(245, 158, 11, 0.4); color: #fbbf24;">\u26a0\ufe0f Anomalies (<span id="gallery-anom-count">0</span>)</button>\n',
    '          </div>\n',
    '          <span class="gallery-count-badge" id="gallery-count-badge">', length(cards), ' shown</span>\n',
    '        </div>\n',
    '        <div class="gallery-toolbar-bottom">\n',
    '          <div class="gallery-filter-group">\n',
    '            <input type="text" id="gallery-search" class="gallery-search-input" placeholder="\U0001f50d Search image stem..." oninput="filterGalleryCombined()" />\n',
    '            <select id="gallery-split-select" class="gallery-select" onchange="filterGalleryCombined()">\n',
    '              <option value="all">All Splits</option>\n',
    '              <option value="train">Train Split</option>\n',
    '              <option value="val">Val Split</option>\n',
    '              <option value="test">Test Split</option>\n',
    '            </select>\n',
    '            <select id="gallery-sort-select" class="gallery-select" onchange="sortGalleryCards(this.value)">\n',
    '              <option value="default">Default Sort</option>\n',
    '              <option value="blur_asc">Most Blurry First</option>\n',
    '              <option value="blur_desc">Sharpest First</option>\n',
    '              <option value="desc">Most Objects First</option>\n',
    '              <option value="asc">Fewest Objects First</option>\n',
    '              <option value="name">Name (A-Z)</option>\n',
    '            </select>\n',
    '          </div>\n',
    '          <div class="gallery-toggles">\n',
    '            <label class="toggle-label">\n',
    '              <input type="checkbox" id="chk-bboxes" checked onchange="toggleOverlay(\'hide-bboxes\', !this.checked)">\n',
    '              <span>BBoxes / Masks</span>\n',
    '            </label>\n',
    '            <label class="toggle-label">\n',
    '              <input type="checkbox" id="chk-labels" checked onchange="toggleOverlay(\'hide-labels\', !this.checked)">\n',
    '              <span>Labels</span>\n',
    '            </label>\n',
    '            <label class="toggle-label">\n',
    '              <input type="checkbox" id="chk-fills" checked onchange="toggleOverlay(\'hide-fills\', !this.checked)">\n',
    '              <span>Shaded Fill</span>\n',
    '            </label>\n',
    '            <div class="lwd-control" style="display:inline-flex; align-items:center; gap:0.45rem; background:var(--bg-surface-elevated); padding:0.25rem 0.6rem; border-radius:6px; border:1px solid var(--border-subtle);">\n',
    '              <span style="font-size:0.75rem; font-weight:600; color:var(--text-muted); text-transform:uppercase;">Line:</span>\n',
    '              <input type="range" id="box-lwd-slider" min="0.5" max="5.0" step="0.5" value="', sprintf("%.1f", lwd), '" style="width:70px; accent-color:var(--color-primary); cursor:pointer;" oninput="updateBoxLwd(this.value)" />\n',
    '              <span id="box-lwd-val" style="font-size:0.75rem; font-weight:700; color:var(--color-primary); min-width:2.4ch;">', sprintf("%.1fpx", lwd), '</span>\n',
    '            </div>\n',
    '          </div>\n',
    '        </div>\n',
    gallery_aug_toolbar_html,
    '      </div>\n',
    '      <div class="gallery-grid" id="main-gallery-grid">\n',
    '        ', paste(gallery_cards_html, collapse = "\n"), '\n',
    '      </div>\n',
    '    </section>\n\n',
    advisor_html,
    '    <section class="section-block" id="health">\n',
    '      <div class="section-header">\n',
    '        <div>\n',
    '          <h2 class="section-title">\U0001fa7a Automated Dataset Health Audit</h2>\n',
    '          <p class="section-desc">Comprehensive diagnostic based on computer vision and YOLO best practices for model convergence and false-positive prevention.</p>\n',
    '        </div>\n',
    '      </div>\n',
    '      <div class="health-grid">\n',
    '        ', paste(health_items_html, collapse = "\n"), '\n',
    '      </div>\n',
    '    </section>\n\n',
    '    <footer class="report-footer">\n',
    '      <div>Generated with <strong>pliman ', pliman_ver, '</strong> \u2022 Automated YOLO Vision Pipeline</div>\n',
    '      <div>', generation_timestamp, '</div>\n',
    '    </footer>\n',
    '  </div>\n\n',
    '  <!-- Lightbox Modal -->\n',
    '  <div class="lightbox-modal" id="lightboxModal" onclick="handleModalBackdropClick(event)">\n',
    '    <div class="lightbox-content">\n',
    '      <button class="lightbox-close" onclick="closeLightbox()">&times;</button>\n',
    '      <div class="lightbox-preview" id="lightboxPreview">\n',
    '        <button class="lightbox-nav lightbox-prev" onclick="prevLightboxImage(event)">&lsaquo;</button>\n',
    '        <div class="lightbox-media-frame" id="lightboxMediaFrame">\n',
    '          <img class="lightbox-img" id="lightboxImg" src="" alt="Sample Preview" />\n',
    '          <div id="lightboxSvgContainer" class="lightbox-svg-layer"></div>\n',
    '        </div>\n',
    '        <button class="lightbox-nav lightbox-next" onclick="nextLightboxImage(event)">&rsaquo;</button>\n',
    '      </div>\n',
    '      <div class="lightbox-sidebar">\n',
    '        <h3 id="lightboxStem" style="font-family:var(--font-heading); font-size:1.3rem; margin-bottom:0.5rem; word-break:break-all;"></h3>\n',
    '        <div style="display:flex; flex-wrap:wrap; gap:0.5rem; margin-bottom:0.75rem;">\n',
    '          <span class="card-badge" id="lightboxCategory"></span>\n',
    '          <span class="card-badge badge-split" id="lightboxSplit"></span>\n',
    '          <span class="meta-tag" id="lightboxDim"></span>\n',
    '          <span class="meta-tag" id="lightboxBlur" style="display:none;"></span>\n',
    '        </div>\n',
    '        <div id="lightboxAugTags" style="display:flex; flex-wrap:wrap; gap:0.4rem; margin-bottom:1rem;"></div>\n',
    '        <div class="lightbox-overlay-controls">\n',
    '          <div class="lightbox-controls-title">\n',
    '            <span>\U0001f441\ufe0f Overlay Controls</span>\n',
    '            <span class="lb-keys-hint">Shortcuts: <b>B</b> Box \u2022 <b>L</b> Label \u2022 <b>F</b> Fill</span>\n',
    '          </div>\n',
    '          <div class="lightbox-toggle-group">\n',
    '            <label class="lb-control-pill active" id="lb-pill-bboxes" title="Toggle Bounding Boxes & Masks (Shortcut: B)">\n',
    '              <input type="checkbox" id="lb-chk-bboxes" checked onchange="setOverlayState(\'hide-bboxes\', !this.checked)" />\n',
    '              <span class="lb-pill-icon">\U0001f532</span>\n',
    '              <span class="lb-pill-text">Box / Mask</span>\n',
    '            </label>\n',
    '            <label class="lb-control-pill active" id="lb-pill-labels" title="Toggle Class Labels (Shortcut: L)">\n',
    '              <input type="checkbox" id="lb-chk-labels" checked onchange="setOverlayState(\'hide-labels\', !this.checked)" />\n',
    '              <span class="lb-pill-icon">\U0001f3f7\ufe0f</span>\n',
    '              <span class="lb-pill-text">Labels</span>\n',
    '            </label>\n',
    '            <label class="lb-control-pill active" id="lb-pill-fills" title="Toggle Shaded Fill (Shortcut: F)">\n',
    '              <input type="checkbox" id="lb-chk-fills" checked onchange="setOverlayState(\'hide-fills\', !this.checked)" />\n',
    '              <span class="lb-pill-icon">\U0001f3a8</span>\n',
    '              <span class="lb-pill-text">Shaded Fill</span>\n',
    '            </label>\n',
    '          </div>\n',
    '          <div class="lightbox-slider-line">\n',
    '            <span class="lightbox-slider-label">Line Width:</span>\n',
    '            <div style="display:flex; align-items:center; gap:0.45rem;">\n',
    '              <input type="range" id="lightboxLwdSlider" min="0.5" max="5.0" step="0.5" value="', sprintf("%.1f", lwd), '" style="width:75px; accent-color:var(--color-primary); cursor:pointer;" oninput="updateBoxLwd(this.value)" />\n',
    '              <span id="lightboxLwdVal" style="font-size:0.75rem; font-weight:700; color:var(--color-primary); min-width:2.4ch;">', sprintf("%.1fpx", lwd), '</span>\n',
    '            </div>\n',
    '          </div>\n',
    '        </div>\n',
    '        <div id="lightboxAnomaliesBox" style="display:none;"></div>\n',
    '        <h4 style="font-size:0.85rem; font-weight:700; text-transform:uppercase; color:var(--text-muted); margin-bottom:0.65rem;">Annotated Objects</h4>\n',
    '        <div id="lightboxObjectsList" style="flex:1; overflow-y:auto;"></div>\n',
    '      </div>\n',
    '    </div>\n',
    '  </div>\n\n',
    '  <!-- Interactive Anomaly Inspector Modal -->\n',
    '  <div class="lightbox-modal" id="anomalyModal" onclick="handleAnomalyBackdropClick(event)">\n',
    '    <div class="anomaly-modal-content">\n',
    '      <button class="lightbox-close" onclick="closeAnomalyModal()">&times;</button>\n',
    '      <div class="anomaly-modal-header">\n',
    '        <h2 class="anomaly-modal-title">\U0001f50d Dataset Anomaly & Integrity Inspector <span id="anomalyModalCount" style="font-size:1rem; font-weight:500; color:var(--text-muted);"></span></h2>\n',
    '        <p class="anomaly-modal-desc">Inspect exact images and bounding box annotations causing integrity warnings. Click thumbnail or "Lightbox" button to inspect annotations full-screen.</p>\n',
    '        <div class="anomaly-filter-tabs" id="anomalyFilterTabs">\n',
    '          <button type="button" class="tab-btn active" data-type="all" onclick="filterAnomalyTable(\'all\', this)">All Anomalies</button>\n',
    '          <button type="button" class="tab-btn" data-type="duplicate_pair" onclick="filterAnomalyTable(\'duplicate_pair\', this)">Duplicate Pairs</button>\n',
    '          <button type="button" class="tab-btn" data-type="extreme_aspect_ratio" onclick="filterAnomalyTable(\'extreme_aspect_ratio\', this)">Extreme AR</button>\n',
    '          <button type="button" class="tab-btn" data-type="micro_box" onclick="filterAnomalyTable(\'micro_box\', this)">Micro-Boxes</button>\n',
    '          <button type="button" class="tab-btn" data-type="out_of_bounds" onclick="filterAnomalyTable(\'out_of_bounds\', this)">Out of Bounds</button>\n',
    '          <button type="button" class="tab-btn" data-type="zero_area" onclick="filterAnomalyTable(\'zero_area\', this)">Zero Area</button>\n',
    '        </div>\n',
    '        <div style="margin-top: 0.5rem;">\n',
    '          <input type="text" id="anomalySearchInput" class="gallery-search-input" placeholder="\U0001f50d Search anomaly by image stem..." oninput="renderAnomalyTable()" style="width: 100%; max-width: 320px;" />\n',
    '        </div>\n',
    '      </div>\n',
    '      <div class="anomaly-modal-body" style="flex:1; overflow-y:auto;">\n',
    '        <div class="table-responsive">\n',
    '          <table class="report-table" id="anomalyTable">\n',
    '            <thead>\n',
    '              <tr>\n',
    '                <th style="width: 65px;">Preview</th>\n',
    '                <th>Image Stem</th>\n',
    '                <th>Split</th>\n',
    '                <th>Anomaly Type</th>\n',
    '                <th>Diagnostic Detail</th>\n',
    '                <th style="text-align: right;">Action</th>\n',
    '              </tr>\n',
    '            </thead>\n',
    '            <tbody id="anomalyTableBody">\n',
    '            </tbody>\n',
    '          </table>\n',
    '        </div>\n',
    '      </div>\n',
    '    </div>\n',
    '  </div>\n\n',
    '  <script>\n',
    '    window.__YOLO_AUDIT_DATA__ = ', audit_json_str, ';\n',
    '    window.__YOLO_ANOMALIES__ = ', anomalies_json_str, ';\n',
    '    ', .yolo_report_js(), '\n',
    '  </script>\n',
    '</body>\n',
    '</html>'
  )
}
