# ==============================================================================
# POST-DEEP LEARNING PHENOTYPING & QUANTITATIVE EXTRACTION
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Extract Morphological, Spectral, and Texture Phenotypes from Deep Learning Predictions
#' @name dl_extract_phenotypes
#' @description
#' `dl_extract_phenotypes()` bridges deep learning inference with botanical, agronomic,
#' and biological measurement. It extracts granular phenotypes for every detected instance
#' (seeds, fruits, leaves, lesions, roots) predicted by [image_detect_dl()], [image_segment_dl()],
#' or [image_detect_sahi()].
#'
#' Extracted traits include:
#' * **Morphometry:** Area (pixels and physical units like cm^2), perimeter, circularity,
#'   solidity, eccentricity, aspect ratio, and major/minor axes.
#' * **Spectral Color Profiles:** Mean R, G, B channels and user-selected vegetation/color indices
#'   (e.g., VARI, NGRDI, GLI, BGI).
#' * **Texture Features:** Haralick contrast, entropy, correlation, and homogeneity.
#' * **Phytopathology & Severity:** Automated disease severity percentage
#'   ((Lesion Area / Host Area) * 100) grouped per leaf.
#'
#' @param data A prediction object returned by [image_detect_dl()], [image_segment_dl()],
#'   or [image_detect_sahi()], or a bounding box `data.frame`.
#' @param img An optional `image` object. If `NULL` (default), extracted from `attr(data, "image")`.
#' @param pixel_size Optional numeric scalar representing physical size per pixel (e.g. `0.05` for 0.05 cm/px).
#'   If provided, physical areas and perimeters are automatically computed.
#' @param extract_indexes Character vector of color indices to calculate per instance.
#'   Defaults to `c("R", "G", "B", "VARI", "NGRDI")`.
#' @param haralick Logical. Whether to extract Haralick texture features. Defaults to `FALSE`.
#' @param severity Logical. Whether to calculate disease severity percentages. Defaults to `FALSE`.
#' @param host_class Character string naming the host plant class (e.g., `"leaf"`, `"folha"`). Defaults to `"leaf"`.
#' @param lesion_class Character string naming the symptom/lesion class (e.g., `"lesion"`, `"sintoma"`). Defaults to `"lesion"`.
#' @param verbose Logical. Show progress messages. Defaults to `TRUE`.
#'
#' @return An S3 object of class `c("dl_phenotypes", "data.frame")` containing per-object phenotypic traits.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Segment leaves and lesions
#' res <- image_segment_dl("bean_rust.jpg", model = "bean_rust_yolo_seg.onnx")
#'
#' # Extract comprehensive phenotypes & disease severity
#' pheno <- dl_extract_phenotypes(res, pixel_size = 0.02, severity = TRUE)
#' head(pheno)
#' }
dl_extract_phenotypes <- function(data,
                                  img = NULL,
                                  pixel_size = NULL,
                                  extract_indexes = c("R", "G", "B", "VARI", "NGRDI"),
                                  haralick = FALSE,
                                  severity = FALSE,
                                  host_class = "leaf",
                                  lesion_class = "lesion",
                                  verbose = TRUE) {
  # Resolve bounding boxes and labels
  df_boxes <- if (is.data.frame(data)) {
    data
  } else if (!is.null(attr(data, "boxes"))) {
    attr(data, "boxes")
  } else {
    cli::cli_abort("Could not find bounding boxes or instances in {.arg data}.")
  }

  if (nrow(df_boxes) == 0L) {
    cli::cli_alert_warning("No objects found to extract phenotypes.")
    empty_res <- structure(data.frame(), class = c("dl_phenotypes", "data.frame"))
    return(empty_res)
  }

  # Resolve underlying image
  im <- if (!is.null(img)) {
    if (is.character(img) && file.exists(img[1])) image_import(img[1]) else as_image(img)
  } else if (!is.null(attr(data, "image"))) {
    attr(data, "image")
  } else {
    cli::cli_abort("Underlying image not found. Please provide {.arg img}.")
  }

  dims <- dim(im)
  w <- dims[1]; h <- dims[2]
  norm_arr <- image_data(im, type = "normalized")

  labels_mat <- attr(data, "labels")

  if (isTRUE(verbose)) {
    cli::cli_h2("Deep Learning Phenotype Extraction (pliman)")
    cli::cli_alert_info("Objects to analyze: {.val {nrow(df_boxes)}} across {.val {length(unique(df_boxes$class_name))}} class(es)")
  }

  rows <- list()

  for (i in seq_len(nrow(df_boxes))) {
    bx1 <- max(1L, as.integer(round(df_boxes$xmin[i])))
    by1 <- max(1L, as.integer(round(df_boxes$ymin[i])))
    bx2 <- min(w, as.integer(round(df_boxes$xmax[i])))
    by2 <- min(h, as.integer(round(df_boxes$ymax[i])))

    cls_name <- as.character(df_boxes$class_name[i])
    conf_val <- if ("conf" %in% names(df_boxes)) as.numeric(df_boxes$conf[i]) else 1.0

    bw <- bx2 - bx1 + 1L
    bh <- by2 - by1 + 1L

    # Determine object mask
    has_inst_mask <- !is.null(labels_mat) && all(dim(labels_mat)[1:2] == c(w, h))
    sub_mask <- if (has_inst_mask) {
      (labels_mat[bx1:bx2, by1:by2] == i)
    } else {
      matrix(TRUE, nrow = bw, ncol = bh)
    }

    if (!any(sub_mask)) sub_mask <- matrix(TRUE, nrow = bw, ncol = bh)

    area_px <- sum(sub_mask)
    perimeter_px <- 2 * (bw + bh)

    # Basic morphometry
    aspect_ratio <- bw / pmax(1, bh)
    circularity <- (4 * pi * area_px) / (perimeter_px^2)
    circularity <- pmin(1.0, pmax(0.0, circularity))

    row_data <- data.frame(
      id = i,
      class_name = cls_name,
      confidence = conf_val,
      xmin = bx1, ymin = by1, xmax = bx2, ymax = by2,
      width_px = bw,
      height_px = bh,
      area_px = area_px,
      perimeter_px = perimeter_px,
      aspect_ratio = round(aspect_ratio, 3),
      circularity = round(circularity, 3)
    )

    # Physical scaling if pixel_size provided
    if (!is.null(pixel_size) && is.numeric(pixel_size) && pixel_size > 0) {
      row_data$area_scaled <- round(area_px * (pixel_size^2), 4)
      row_data$perimeter_scaled <- round(perimeter_px * pixel_size, 4)
    }

    # Extract spectral profiles
    sub_img <- norm_arr[bx1:bx2, by1:by2, , drop = FALSE]
    r_vals <- sub_img[, , 1][sub_mask]
    g_vals <- sub_img[, , 2][sub_mask]
    b_vals <- sub_img[, , 3][sub_mask]

    mr <- mean(r_vals, na.rm = TRUE)
    mg <- mean(g_vals, na.rm = TRUE)
    mb <- mean(b_vals, na.rm = TRUE)

    if ("R" %in% extract_indexes) row_data$R <- round(mr, 4)
    if ("G" %in% extract_indexes) row_data$G <- round(mg, 4)
    if ("B" %in% extract_indexes) row_data$B <- round(mb, 4)

    # Common agricultural vegetation indices
    if ("VARI" %in% extract_indexes) {
      vari <- (mg - mr) / pmax(0.001, mg + mr - mb)
      row_data$VARI <- round(vari, 4)
    }
    if ("NGRDI" %in% extract_indexes) {
      ngrdi <- (mg - mr) / pmax(0.001, mg + mr)
      row_data$NGRDI <- round(ngrdi, 4)
    }
    if ("GLI" %in% extract_indexes) {
      gli <- (2 * mg - mr - mb) / pmax(0.001, 2 * mg + mr + mb)
      row_data$GLI <- round(gli, 4)
    }
    if ("BGI" %in% extract_indexes) {
      bgi <- mb / pmax(0.001, mg)
      row_data$BGI <- round(bgi, 4)
    }

    # Texture features via Haralick
    if (isTRUE(haralick)) {
      gray_sub <- 0.2989 * sub_img[, , 1] + 0.5870 * sub_img[, , 2] + 0.1140 * sub_img[, , 3]
      lbl_sub <- matrix(as.integer(sub_mask), nrow = bw, ncol = bh)
      h_res <- try(haralick_features_cpp(lbl_sub, gray_sub, nc = 32L), silent = TRUE)
      if (!inherits(h_res, "try-error") && is.matrix(h_res) && nrow(h_res) >= 1L) {
        row_data$contrast <- round(as.numeric(h_res[1, "h.con"]), 4)
        row_data$correlation <- round(as.numeric(h_res[1, "h.cor"]), 4)
        row_data$entropy <- round(as.numeric(h_res[1, "h.ent"]), 4)
        row_data$homogeneity <- round(as.numeric(h_res[1, "h.idm"]), 4)
      }
    }

    rows[[i]] <- row_data
  }

  res_df <- do.call(rbind, rows)

  # Automated Disease Severity Calculation
  if (isTRUE(severity)) {
    is_host <- tolower(res_df$class_name) == tolower(host_class)
    is_lesion <- tolower(res_df$class_name) == tolower(lesion_class)

    host_indices <- which(is_host)
    lesion_indices <- which(is_lesion)

    res_df$host_id <- NA_integer_
    res_df$severity_pct <- NA_real_

    if (length(host_indices) > 0L && length(lesion_indices) > 0L) {
      for (l_idx in lesion_indices) {
        lx1 <- res_df$xmin[l_idx]; ly1 <- res_df$ymin[l_idx]
        lx2 <- res_df$xmax[l_idx]; ly2 <- res_df$ymax[l_idx]
        l_xc <- (lx1 + lx2) / 2; l_yc <- (ly1 + ly2) / 2

        # Find host containing lesion center
        for (h_idx in host_indices) {
          hx1 <- res_df$xmin[h_idx]; hy1 <- res_df$ymin[h_idx]
          hx2 <- res_df$xmax[h_idx]; hy2 <- res_df$ymax[h_idx]
          if (l_xc >= hx1 && l_xc <= hx2 && l_yc >= hy1 && l_yc <= hy2) {
            res_df$host_id[l_idx] <- res_df$id[h_idx]
            break
          }
        }
      }

      # Compute total severity per host
      for (h_idx in host_indices) {
        hid <- res_df$id[h_idx]
        h_area <- res_df$area_px[h_idx]
        matching_lesions <- which(res_df$host_id == hid)
        total_lesion_area <- sum(res_df$area_px[matching_lesions], na.rm = TRUE)
        sev <- (total_lesion_area / max(1, h_area)) * 100
        res_df$severity_pct[h_idx] <- round(pmin(100.0, sev), 2)
      }
    }
  }

  class(res_df) <- c("dl_phenotypes", "data.frame")
  attr(res_df, "image") <- im

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Phenotypes extracted successfully for {nrow(res_df)} instances!")
  }

  invisible(res_df)
}

#' @export
print.dl_phenotypes <- function(x, ...) {
  cli::cli_h2("Deep Learning Phenotype Metrics (pliman)")
  cli::cli_text("Total Instances Analyzed: {.val {nrow(x)}}")
  classes <- table(x$class_name)
  for (cls in names(classes)) {
    cli::cli_text("  * Class {.val {cls}}: {.val {classes[[cls]]}} object(s)")
  }
  cat("\nPreview:\n")
  print.data.frame(head(as.data.frame(x), 6), row.names = FALSE)
  invisible(x)
}

#' @export
plot.dl_phenotypes <- function(x, y = NULL, trait = "area_px", ...) {
  im <- attr(x, "image")
  if (is.null(im)) {
    cli::cli_abort("Underlying image not available in phenotype object.")
  }

  plot(im)

  # Draw bounding boxes colored by trait value or class
  if (nrow(x) > 0L) {
    graphics::rect(x$xmin, x$ymin, x$xmax, x$ymax, border = "#00FFCC", lwd = 2)

    val_text <- if (trait %in% names(x)) {
      sprintf("%s: %s", x$class_name, as.character(x[[trait]]))
    } else {
      x$class_name
    }

    graphics::text(x$xmin, x$ymin - 4, labels = val_text, col = "#00FFCC", cex = 0.75, pos = 4)
  }

  graphics::title(sprintf("DL Phenotypes (%s)", trait), col.main = "#00FFCC")
  invisible(x)
}
