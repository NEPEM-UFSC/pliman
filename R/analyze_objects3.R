#' Analyzes objects in an image (Hyper-Optimized Engine v3)
#'
#' @description
#' `analyze_objects3()` is the highest-performance computer vision pipeline in `pliman`,
#' utilizing SIMD vectorization, division-free Otsu thresholding, 2D distance transform
#' with column linear scans, Akl-Toussaint heuristic convex hulls, and a unified zero-allocation C++ engine.
#'
#' @inheritParams analyze_objects
#' @param area_mode The method to compute object area. Default is `"contour"` (Shoelace area
#'   of the boundary polygon). If `"pixel"`, counts the exact number of labeled mask pixels
#'   belonging to the object.
#' @return An object of class `anal_obj`. See [analyze_objects()] for details.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
analyze_objects3 <- function(img,
                             foreground = NULL,
                             background = NULL,
                             pick_palettes = FALSE,
                             segment_objects = TRUE,
                             viewer = get_pliman_viewer(),
                             reference = FALSE,
                             reference_area = NULL,
                             back_fore_index = "R/(G/B)",
                             fore_ref_index = "B-R",
                             reference_img = NULL,
                             reference_larger = FALSE,
                             reference_smaller = FALSE,
                             pattern = NULL,
                             parallel = FALSE,
                             workers = NULL,
                             watershed = TRUE,
                             veins = FALSE,
                             sigma_veins = 1,
                             veins_thinning = FALSE,
                             veins_rel_erode = 0.005,
                             ab_angles = FALSE,
                             ab_angles_percentiles = c(0.25, 0.75),
                             width_at = FALSE,
                             width_at_percentiles = c(0.05, 0.25, 0.50, 0.75, 0.95),
                             haralick = FALSE,
                             har_nbins = 32,
                             har_scales = 1,
                             har_band = 1,
                             smooth = FALSE,
                             contours = TRUE,
                             pcv = FALSE,
                             pcv_niter = 100,
                             resize = FALSE,
                             trim = FALSE,
                             fill_hull = FALSE,
                             erode = FALSE,
                             dilate = FALSE,
                             opening = FALSE,
                             closing = FALSE,
                             filter = FALSE,
                             return_exact = FALSE,
                             filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                             invert = FALSE,
                             sensitivity = "medium",
                             area_mode = c("contour", "pixel"),
                             index = "NB",
                             r = 1,
                             g = 2,
                             b = 3,
                             re = 4,
                             nir = 5,
                             swir = 6,
                             object_index = NULL,
                             pixel_level_index = FALSE,
                             return_mask = FALSE,
                             efourier = FALSE,
                             nharm = 10,
                             threshold = "Otsu",
                             k = 0.1,
                             windowsize = NULL,
                             tolerance = NULL,
                             extension = NULL,
                             lower_noise = 0.10,
                             rel_size = NULL,
                             lower_size = NULL,
                             upper_size = NULL,
                             topn_lower = NULL,
                             topn_upper = NULL,
                             lower_eccent = NULL,
                             upper_eccent = NULL,
                             lower_circ = NULL,
                             upper_circ = NULL,
                             randomize = TRUE,
                             nrows = 1000,
                             plot = TRUE,
                             show_original = TRUE,
                             show_chull = FALSE,
                             show_contour = FALSE,
                             show_bbox = FALSE,
                             contour_col = "red",
                             contour_size = 1,
                             show_lw = FALSE,
                             show_background = TRUE,
                             show_segmentation = FALSE,
                             col_foreground = NULL,
                             col_background = NULL,
                             marker = "point",
                             marker_col = "red",
                             marker_size = NULL,
                             save_image = FALSE,
                             max_pixels = 2e6,
                             prefix = "proc_",
                             dir_original = NULL,
                             dir_processed = NULL,
                             verbose = TRUE){

  area_mode <- area_mode[[1]]
  if (!area_mode %in% c("contour", "pixel")) {
    cli::cli_abort("Argument {.arg area_mode} must be one of {.val contour} or {.val pixel}.")
  }

  if (!is.null(rel_size)) {
    lower_noise <- rel_size
  }
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)

  if (!missing(img) && !missing(pattern)) {
    cli::cli_abort("Only one of {.arg img} or {.arg pattern} can be used.")
  }

  diretorio_original <- if (is.null(dir_original)) "./" else (if (grepl("[/\\]", dir_original)) dir_original else paste0("./", dir_original))
  diretorio_processada <- if (is.null(dir_processed)) "./" else (if (grepl("[/\\]", dir_processed)) dir_processed else paste0("./", dir_processed))

  help_count3 <- function(img, foreground = NULL, background = NULL, pick_palettes = FALSE, resize = FALSE, trim = FALSE, fill_hull = FALSE, threshold = "Otsu",
                          erode = FALSE, dilate = FALSE, opening = FALSE, closing = FALSE, filter = FALSE, filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                          tolerance = 0.2, extension = 1, randomize = TRUE, nrows = 1000, plot = TRUE, show_original = TRUE, show_chull = FALSE,
                          show_contour = FALSE, show_bbox = FALSE, contour_col = "red", contour_size = 1, show_lw = FALSE,
                          show_background = TRUE, show_segmentation = FALSE, marker = "point", marker_col = "red", marker_size = NULL, save_image = FALSE,
                          prefix = "proc_", dir_original = NULL, dir_processed = NULL, verbose = TRUE, col_background = NULL,
                          col_foreground = NULL, lower_noise = 0.1, lower_size = NULL, upper_size = NULL, topn_lower = NULL, topn_upper = NULL,
                          lower_eccent = NULL, upper_eccent = NULL, lower_circ = NULL, upper_circ = NULL,
                          ab_angles = FALSE, ab_angles_percentiles = c(0.25, 0.75), width_at = FALSE, width_at_percentiles = c(0.05, 0.25, 0.50, 0.75, 0.95),
                          return_mask = FALSE, pcv = FALSE, pcv_niter = 100, object_index = NULL, max_pixels = 2e6, return_exact = FALSE, area_mode = "contour",
                          watershed = TRUE, reference = FALSE, reference_area = NULL, reference_larger = FALSE, reference_smaller = FALSE,
                          veins = FALSE, sigma_veins = 1, veins_thinning = FALSE, veins_rel_erode = 0.005,
                          haralick = FALSE, har_nbins = 32, har_scales = 1, har_band = 1, smooth = FALSE, contours = TRUE,
                          index = "NB", r = 1, g = 2, b = 3, re = 4, nir = 5, swir = 6, invert = FALSE, k = 0.1, windowsize = NULL,
                          pixel_level_index = FALSE, efourier = FALSE, nharm = 10, back_fore_index = "R/(G/B)", viewer = get_pliman_viewer()) {

    if (is.character(img)) {
      all_files <- sapply(list.files(diretorio_original), file_name)
      check_names_dir(img, all_files, diretorio_original)
      imag <- list.files(diretorio_original, pattern = paste0("^", img, "\\."))
      name_ori <- file_name(imag)
      extens_ori <- file_extension(imag)
      img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
    } else {
      name_ori <- as.character(substitute(img))[1]
      extens_ori <- "png"
    }

    im_dat <- if (inherits(img, "Image")) img@.Data else image_data(img)
    d_img <- dim(im_dat)
    if (length(d_img) == 2) {
      arr <- array(c(im_dat, im_dat, im_dat), dim = c(d_img[1], d_img[2], 3))
      img <- as_image(arr)
      im_dat <- img@.Data
    } else if (length(d_img) == 3 && d_img[3] == 1) {
      im_data_1 <- im_dat[,,1]
      arr <- array(c(im_data_1, im_data_1, im_data_1), dim = c(d_img[1], d_img[2], 3))
      img <- as_image(arr)
      im_dat <- img@.Data
    }

    if (trim != FALSE) {
      if (!is.numeric(trim)) cli::cli_abort("Argument {.arg trim} must be numeric.")
      img <- image_trim(img, trim)
      im_dat <- if (inherits(img, "Image")) img@.Data else image_data(img)
    }
    if (resize != FALSE) {
      if (!is.numeric(resize)) cli::cli_abort("Argument {.arg resize} must be numeric.")
      img <- image_resize(img, resize)
      im_dat <- if (inherits(img, "Image")) img@.Data else image_data(img)
    }

    if (is.null(tolerance)) tolerance <- 0.2
    if (is.null(extension)) extension <- 1

    bin_precomputed <- NULL

    if (isFALSE(reference)) {
      if (isTRUE(pick_palettes) && interactive()) {
        background <- pick_palette(img, r = 5, verbose = FALSE, palette = FALSE, plot = FALSE, col = "blue", external_device = FALSE, title = "Pick BACKGROUND colors", viewer = viewer)
        foreground <- pick_palette(img, r = 5, verbose = FALSE, palette = FALSE, plot = FALSE, col = "salmon", external_device = FALSE, title = "Pick FOREGROUND colors", viewer = viewer)
      }

      if (!is.null(foreground) && !is.null(background)) {
        if (is.character(foreground)) foreground <- image_import(foreground)
        if (is.character(background)) background <- image_import(background)
        fore_num <- image_data(foreground, type = "numeric")
        back_num <- image_data(background, type = "numeric")
        df_fore <- data.frame(CODE = "foreground", R = c(fore_num[,,1]), G = c(fore_num[,,2]), B = c(fore_num[,,3]))
        df_back <- data.frame(CODE = "background", R = c(back_num[,,1]), G = c(back_num[,,2]), B = c(back_num[,,3]))
        back_fore <- transform(rbind(df_fore[sample(1:nrow(df_fore)),][1:nrows,],
                                     df_back[sample(1:nrow(df_back)),][1:nrows,]),
                               Y = ifelse(CODE == "background", 0, 1))
        formula_glm <- as.formula(paste("Y ~ ", back_fore_index))
        modelo1 <- suppressWarnings(glm(formula_glm, family = binomial("logit"), data = back_fore))
        bin_precomputed <- predict_binary_glm(modelo1, img) == 1
      }
    }

    thresh_method <- if (is.character(threshold)) threshold[[1]] else "numeric"
    thresh_num <- if (is.numeric(threshold)) threshold[[1]] else 0.5
    win_sz <- if (is.null(windowsize)) 0L else as.integer(windowsize)
    smooth_val <- if (is.numeric(smooth)) as.integer(smooth) else 0L

    tol_val <- if (is.null(tolerance)) {
      if (is.character(sensitivity)) {
        sens_str <- tolower(trimws(sensitivity[1]))
        lvl <- switch(sens_str,
                      "low"     = 1.0,
                      "medium"  = 2.0,
                      "high"    = 3.0,
                      "extreme" = 4.0,
                      cli::cli_abort("Invalid {.arg sensitivity} '{sens_str}'."))
      } else if (is.numeric(sensitivity)) {
        lvl <- as.double(sensitivity[1])
      } else {
        lvl <- 2.0
      }
      0.20 / (2.0 ^ (lvl - 1.0))
    } else {
      as.double(tolerance)
    }

    ext_val <- if (is.null(extension)) 1L else as.integer(extension)

    # Chamada direta ao motor C++ analyze_objects3_cpp
    cpp_res <- analyze_objects3_cpp(
      img_sexp = im_dat,
      bin_sexp = if (!is.null(bin_precomputed)) bin_precomputed else NULL,
      index_str = as.character(index),
      r = as.integer(r), g = as.integer(g), b = as.integer(b),
      re = as.integer(re), nir = as.integer(nir), swir = as.integer(swir),
      threshold_method = thresh_method,
      threshold_val = thresh_num,
      k_adj = as.numeric(k),
      windowsize = win_sz,
      invert = isTRUE(invert),
      erode_sz = if (is.numeric(erode)) as.integer(erode) else 0L,
      dilate_sz = if (is.numeric(dilate)) as.integer(dilate) else 0L,
      opening_sz = if (is.numeric(opening)) as.integer(opening) else 0L,
      closing_sz = if (is.numeric(closing)) as.integer(closing) else 0L,
      filter_sz = if (is.numeric(filter)) as.integer(filter) else 0L,
      fill_hull = isTRUE(fill_hull),
      filter_order = filter_order,
      return_exact = isTRUE(return_exact),
      watershed = isTRUE(watershed),
      tolerance = tol_val,
      ext = ext_val,
      haralick = isTRUE(haralick),
      har_nbins = as.integer(har_nbins),
      har_band = if (is.numeric(har_band)) as.integer(har_band) else 1L,
      smooth = smooth_val,
      return_contours = isTRUE(contours)
    )

    nmask <- cpp_res[["labels"]]
    ocont <- cpp_res[["contours"]]
    shape <- cpp_res[["shape"]]

    valid <- which(!is.na(shape$x))
    if (length(valid) < nrow(shape)) {
      shape <- shape[valid, , drop = FALSE]
      if (length(ocont) > 0) ocont <- ocont[valid]
    }
    if (length(ocont) > 0) names(ocont) <- shape$id

    if (area_mode == "pixel") {
      pix_areas <- get_area_mask(nmask)
      shape$area <- pix_areas[shape$id]
      shape$coverage <- shape$area / (nrow(nmask) * ncol(nmask))
      shape$form_factor <- 4 * pi * shape$area / (shape$perimeter ^ 2)
      shape$rectangularity <- (shape$length * shape$width) / shape$area
      shape$solidity <- shape$area / shape$area_ch
      shape$circularity <- (shape$perimeter ^ 2) / shape$area
      shape$circularity_norm <- (shape$area * 4 * pi) / (shape$perimeter ^ 2)
    }

    if (isTRUE(haralick) && !is.null(cpp_res[["haralick"]]) && nrow(cpp_res[["haralick"]]) > 0) {
      hal <- as.data.frame(cpp_res[["haralick"]])
      shape <- cbind(shape, hal[valid, , drop = FALSE])
      colnames(shape) <- c(names_measures(), har_names())
    }

    ch <- if (!is.null(cpp_res[["chull"]]) && length(cpp_res[["chull"]]) > 0) {
      if (length(valid) < length(cpp_res[["chull"]])) cpp_res[["chull"]][valid] else cpp_res[["chull"]]
    } else if (length(ocont) > 0) {
      conv_hull(ocont)
    } else {
      NULL
    }
    if (!is.null(ch) && length(ch) > 0) names(ch) <- shape$id

    # Reference Correction
    if (isTRUE(reference)) {
      if (is.null(reference_area)) cli::cli_abort("A known {.strong area} must be declared when a reference object is used.")
      if (isFALSE(reference_larger) && isFALSE(reference_smaller)) {
        cli::cli_abort("When {.code reference = TRUE}, one of {.arg reference_larger} or {.arg reference_smaller} must be TRUE.")
      }
      r_val <- if (is.raw(img)) 255 else 1
      if (isTRUE(reference_larger)) {
        lineid <- which.max(shape$area)
        id_ref <- shape[lineid, "id"]
        ref_measures <- shape[lineid, , drop = FALSE]
        rownames(ref_measures) <- NULL
        pix_ref <- which(nmask == id_ref)
        img[,,1][pix_ref] <- r_val
        img[,,2][pix_ref] <- 0
        img[,,3][pix_ref] <- 0
        npix_ref <- shape[lineid, "area"]
        shape <- shape[-lineid, ]
        shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
      } else if (isTRUE(reference_smaller)) {
        shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
        lineid <- which.min(shape$area)
        id_ref <- shape[lineid, "id"]
        ref_measures <- shape[lineid, , drop = FALSE]
        rownames(ref_measures) <- NULL
        pix_ref <- which(nmask == id_ref)
        img[,,1][pix_ref] <- r_val
        img[,,2][pix_ref] <- 0
        img[,,3][pix_ref] <- 0
        npix_ref <- shape[lineid, "area"]
        shape <- shape[-lineid, ]
      }
      if (exists("npix_ref")) {
        px_side <- sqrt(reference_area / npix_ref)
        shape$area <- shape$area * px_side ^ 2
        shape$area_ch <- shape$area_ch * px_side ^ 2
        shape[6:18] <- apply(shape[6:18], 2, function(x) x * px_side)
      }
    }

    # Size and Shape Filters
    if (nrow(shape) > 0) {
      if (!is.null(lower_size)) {
        shape <- shape[which(shape$area > lower_size), , drop = FALSE]
      } else if (!is.null(lower_noise)) {
        shape <- shape[which(shape$area > mean(shape$area, na.rm = TRUE) * lower_noise), , drop = FALSE]
      }
      if (!is.null(upper_size) && nrow(shape) > 0) shape <- shape[which(shape$area < upper_size), , drop = FALSE]
      if (!is.null(topn_lower) && nrow(shape) > 0) shape <- shape[order(shape$area), , drop = FALSE][1:min(nrow(shape), topn_lower), , drop = FALSE]
      if (!is.null(topn_upper) && nrow(shape) > 0) shape <- shape[order(shape$area, decreasing = TRUE), , drop = FALSE][1:min(nrow(shape), topn_upper), , drop = FALSE]
      if (!is.null(lower_eccent) && nrow(shape) > 0) shape <- shape[which(shape$eccentricity > lower_eccent), , drop = FALSE]
      if (!is.null(upper_eccent) && nrow(shape) > 0) shape <- shape[which(shape$eccentricity < upper_eccent), , drop = FALSE]
      if (!is.null(lower_circ) && nrow(shape) > 0) shape <- shape[which(shape$circularity > lower_circ), , drop = FALSE]
      if (!is.null(upper_circ) && nrow(shape) > 0) shape <- shape[which(shape$circularity < upper_circ), , drop = FALSE]
    }

    if (nrow(shape) < length(ocont) && nrow(shape) > 0 && length(ocont) > 0) {
      keep_idx <- match(shape$id, names(ocont))
      ocont <- ocont[keep_idx]
      if (!is.null(ch) && length(ch) > 0) ch <- ch[keep_idx]
      filter_labels_cpp(nmask, shape$id)
    } else if (nrow(shape) == 0) {
      ocont <- list()
      ch <- list()
      nmask[] <- 0L
    }

    # Secondary Analyses
    if (isTRUE(efourier) && nrow(shape) > 0) {
      efr <- efourier(ocont, nharm = nharm)
      efer <- efourier_error(efr, plot = FALSE)$stats
      efpowwer <- efourier_power(efr, plot = FALSE)
      efpow <- efpowwer$cum_power
      min_harm <- efpowwer$min_harm
      efrn <- efourier_norm(efr)
      efr <- efourier_coefs(efr)
      names(efr)[1] <- "id"
      efrn <- efourier_coefs(efrn)
      names(efrn)[1] <- "id"
    } else {
      efr <- efrn <- efer <- efpow <- min_harm <- NULL
    }

    angles <- if (isTRUE(ab_angles) && nrow(shape) > 0) poly_apex_base_angle(ocont, ab_angles_percentiles) else NULL
    widths <- if (isTRUE(width_at) && nrow(shape) > 0) poly_width_at(ocont, width_at_percentiles) else NULL

    if (isTRUE(veins) && nrow(shape) > 0) {
      vein_result <- detect_veins_cpp(
        R_sexp = img[,,1], G_sexp = img[,,2], B_sexp = img[,,3],
        labels_sexp = as.matrix(nmask),
        sigma1 = sigma_veins, sigma2 = sigma_veins * 5,
        threshold = -1, channel = 0, erode_size = -1,
        rel_erode = veins_rel_erode, thinning = veins_thinning,
        return_map = FALSE
      )
      vein_props <- vein_result$proportion
      valid_labs <- which(!is.na(vein_props) & seq_along(vein_props) > 0)
      prop_veins <- data.frame(id = valid_labs, prop_veins = vein_props[valid_labs])
      prop_veins <- prop_veins[prop_veins$id %in% shape$id, ]
    } else {
      prop_veins <- NULL
    }

    pcv_res <- if (isTRUE(pcv) && nrow(shape) > 0) poly_pcv(ocont, niter = pcv_niter) else NULL

    # Object Indexing
    if (!is.null(object_index) && nrow(shape) > 0) {
      object_index_used <- object_index[1]
      ind_df <- get_cached_indexes_df()
      ind_formula <- ifelse(object_index %in% ind_df$Index, ind_df[match(object_index, ind_df$Index), 2], object_index)
      R_ch <- img[,,1]; G_ch <- img[,,2]; B_ch <- img[,,3]
      hsb_arrays <- rgb_to_hsb_cpp(R_ch, G_ch, B_ch)
      eval_env <- list(
        R = R_ch, G = G_ch, B = B_ch,
        h = hsb_arrays$H, s = hsb_arrays$S, b = hsb_arrays$B,
        H = hsb_arrays$H, S = hsb_arrays$S
      )
      parsed <- lapply(ind_formula, function(f) parse(text = f)[[1]])
      idx_mat <- do.call(cbind, lapply(parsed, function(expr) {
        as.vector(eval(expr, envir = list2env(eval_env, parent = baseenv())))
      }))
      valid_ids <- as.integer(shape$id)
      means_mat <- compute_index_means_cpp(idx_mat, as.vector(nmask), valid_ids)
      colnames(means_mat) <- object_index
      indexes <- data.frame(id = valid_ids, means_mat, check.names = FALSE)
      if (isTRUE(pixel_level_index)) {
        obj_rgb <- object_rgb(img, nmask)
        obj_rgb <- subset(obj_rgb, id %in% shape$id)
        obj_rgb <- cbind(obj_rgb, rgb_to_hsb(obj_rgb[, 2:4]))
        pixel_idx <- do.call(cbind, lapply(parsed, function(expr) {
          as.vector(eval(expr, envir = list2env(eval_env, parent = baseenv())))
        }))
        mask_flat <- as.vector(nmask)
        pixel_keep <- mask_flat %in% valid_ids
        obj_rgb <- cbind(obj_rgb, as.data.frame(pixel_idx[pixel_keep, , drop = FALSE], col.names = object_index))
      } else {
        obj_rgb <- NULL
      }
    } else {
      obj_rgb <- NULL
      indexes <- NULL
      object_index_used <- NULL
    }

    mask_out <- if (isTRUE(return_mask)) nmask else NULL

    n_obj  <- nrow(shape)
    a_mean <- if (n_obj > 0) round(mean(shape$area, na.rm = TRUE), 1) else 0
    a_min  <- if (n_obj > 0) round(min(shape$area, na.rm = TRUE), 1) else 0
    a_max  <- if (n_obj > 0) round(max(shape$area, na.rm = TRUE), 1) else 0
    a_sd   <- if (n_obj > 0) round(sd(shape$area, na.rm = TRUE), 1) else 0
    a_sum  <- if (n_obj > 0) round(sum(shape$area, na.rm = TRUE), 1) else 0
    cov_val <- if (n_obj > 0) sum(shape$coverage) else 0

    stats <- data.frame(
      stat = c("n", "min_area", "mean_area", "max_area", "sd_area", "sum_area", "coverage"),
      value = c(n_obj, a_min, a_mean, a_max, a_sd, a_sum, cov_val)
    )

    results <- list(
      results = shape, statistics = stats, object_rgb = obj_rgb,
      object_index = indexes, efourier = efr, efourier_norm = efrn,
      efourier_error = efer, efourier_power = efpow, efourier_minharm = min_harm,
      veins = prop_veins, angles = angles, width_at = widths, mask = mask_out,
      pcv = pcv_res, contours = ocont,
      parms = list(
        index = index,
        object_index = object_index_used,
        reference = isTRUE(reference),
        reference_area = reference_area,
        reference_larger = reference_larger,
        reference_smaller = reference_smaller,
        npix_ref = if (exists("npix_ref")) npix_ref else NULL,
        px_side = if (exists("px_side")) px_side else NULL,
        reference_measures = if (exists("ref_measures")) ref_measures else NULL
      )
    )
    class(results) <- "anal_obj"

    if (plot == TRUE || save_image == TRUE) {
      backg <- !is.null(col_background)
      col_bg <- if (is.null(col_background)) col2rgb("white") else (if (is.character(col_background)) col2rgb(col_background) else col_background)
      col_fg <- if (is.null(col_foreground)) col2rgb("gray") else (if (is.character(col_foreground)) col2rgb(col_foreground) else col_foreground)
      img_max <- max(as.numeric(image_data(img)))
      if (max(col_bg) > 1 && img_max <= 1) col_bg <- col_bg / 255
      if (max(col_fg) > 1 && img_max <= 1) col_fg <- col_fg / 255

      ID <- which(nmask != 0)
      ID2 <- which(nmask == 0)

      if (show_original == TRUE && show_segmentation == FALSE) {
        im2 <- img[,,1:3]
        if (backg) {
          im3 <- image_color_labels(nmask)
          im2[,,1][which(im3[,,1] == 0)] <- col_bg[1]
          im2[,,2][which(im3[,,2] == 0)] <- col_bg[2]
          im2[,,3][which(im3[,,3] == 0)] <- col_bg[3]
        }
      } else if (show_original == TRUE && show_segmentation == TRUE) {
        im2 <- image_color_labels(nmask)
        if (backg) {
          im2[,,1][which(im2[,,1] == 0)] <- col_bg[1]
          im2[,,2][which(im2[,,2] == 0)] <- col_bg[2]
          im2[,,3][which(im2[,,3] == 0)] <- col_bg[3]
        } else {
          im2[,,1][which(im2[,,1] == 0)] <- as.numeric(img[,,1])[which(im2[,,1] == 0)]
          im2[,,2][which(im2[,,2] == 0)] <- as.numeric(img[,,2])[which(im2[,,2] == 0)]
          im2[,,3][which(im2[,,3] == 0)] <- as.numeric(img[,,3])[which(im2[,,3] == 0)]
        }
      } else {
        if (show_segmentation == TRUE) {
          im2 <- image_color_labels(nmask)
          im2[,,1][which(im2[,,1] == 0)] <- col_bg[1]
          im2[,,2][which(im2[,,2] == 0)] <- col_bg[2]
          im2[,,3][which(im2[,,3] == 0)] <- col_bg[3]
        } else {
          im2 <- img[,,1:3]
          im2[,,1][ID] <- col_fg[1]
          im2[,,2][ID] <- col_fg[2]
          im2[,,3][ID] <- col_fg[3]
          im2[,,1][ID2] <- col_bg[1]
          im2[,,2][ID2] <- col_bg[2]
          im2[,,3][ID2] <- col_bg[3]
        }
      }

      show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
      marker_val <- ifelse(is.null(marker), "id", marker)
      marker_col_val <- ifelse(is.null(marker_col), "white", marker_col)
      marker_sz <- ifelse(is.null(marker_size), 0.75, marker_size)

      if (plot == TRUE) {
        plot(im2, max_pixels = max_pixels)
        if (nrow(shape) > 0) {
          if (isTRUE(show_contour)) plot_contour(ocont, col = contour_col, lwd = contour_size)
          if (show_bbox) plot_bbox(ocont, col = contour_col)
          if (show_mark && marker_val != "point") {
            text(shape[, 2], shape[, 3], round(shape[, marker_val], 3), col = marker_col_val, cex = marker_sz)
          } else if (show_mark && marker_val == "point") {
            points(shape[, 2], shape[, 3], col = marker_col_val, pch = 16, cex = marker_sz)
          }
          if (isTRUE(show_lw)) plot_lw(results)
          if (show_chull) plot_contour(ch, col = contour_col, lwd = contour_size)
        }
      }

      if (save_image == TRUE) {
        dim_img <- dim(im2)
        dev.new(width = dim_img[[1]], height = dim_img[[2]], noRStudioGD = TRUE, units = "px")
        plot(im2, max_pixels = max_pixels)
        if (nrow(shape) > 0) {
          if (isTRUE(show_contour)) plot_contour(ocont, col = contour_col, lwd = contour_size)
          if (show_bbox) plot_bbox(ocont, col = contour_col)
          if (show_mark && marker_val != "point") {
            text(shape[, 2], shape[, 3], round(shape[, marker_val], 3), col = marker_col_val, cex = marker_sz)
          } else if (show_mark && marker_val == "point") {
            points(shape[, 2], shape[, 3], col = marker_col_val, pch = 16, cex = marker_sz)
          }
          if (isTRUE(show_lw)) plot_lw(results)
          if (show_chull) plot_contour(ch, col = contour_col, lwd = contour_size)
        }
        fig_name <- paste0(prefix, name_ori, ".", extens_ori)
        dev.print(device = png, file = file.path(diretorio_processada, fig_name), width = dim_img[[1]], height = dim_img[[2]])
        dev.off()
      }
    }

    invisible(results)
  }

  if (missing(pattern)) {
    if (verbose) cli::cli_progress_step("{.pkg pliman} is processing the image. Please wait.")
    help_count3(
      img = img, foreground = foreground, background = background, pick_palettes = pick_palettes, resize = resize, trim = trim, fill_hull = fill_hull, threshold = threshold,
      erode = erode, dilate = dilate, opening = opening, closing = closing, filter = filter, filter_order = filter_order,
      tolerance = tolerance, extension = extension, randomize = randomize, nrows = nrows, plot = plot, show_original = show_original, show_chull = show_chull,
      show_contour = show_contour, show_bbox = show_bbox, contour_col = contour_col, contour_size = contour_size, show_lw = show_lw,
      show_background = show_background, show_segmentation = show_segmentation, marker = marker, marker_col = marker_col, marker_size = marker_size, save_image = save_image,
      prefix = prefix, dir_original = dir_original, dir_processed = dir_processed, verbose = verbose, col_background = col_background,
      col_foreground = col_foreground, lower_noise = lower_noise, lower_size = lower_size, upper_size = upper_size, topn_lower = topn_lower, topn_upper = topn_upper,
      lower_eccent = lower_eccent, upper_eccent = upper_eccent, lower_circ = lower_circ, upper_circ = upper_circ,
      ab_angles = ab_angles, ab_angles_percentiles = ab_angles_percentiles, width_at = width_at, width_at_percentiles = width_at_percentiles,
      return_mask = return_mask, pcv = pcv, pcv_niter = pcv_niter, object_index = object_index, max_pixels = max_pixels, return_exact = return_exact, area_mode = area_mode,
      watershed = watershed, reference = reference, reference_area = reference_area, reference_larger = reference_larger, reference_smaller = reference_smaller,
      veins = veins, sigma_veins = sigma_veins, veins_thinning = veins_thinning, veins_rel_erode = veins_rel_erode,
      haralick = haralick, har_nbins = har_nbins, har_scales = har_scales, har_band = har_band, smooth = smooth, contours = contours,
      index = index, r = r, g = g, b = b, re = re, nir = nir, swir = swir, invert = invert, k = k, windowsize = windowsize,
      pixel_level_index = pixel_level_index, efourier = efourier, nharm = nharm, back_fore_index = back_fore_index, viewer = viewer
    )
  } else {
    if (pattern %in% as.character(0:9)) {
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    imgpath <- file.path(getwd(), sub('./', '', diretorio_original)) |> trunc_path(max_chars = 50)

    if (length(grep(pattern, names_plant)) == 0) {
      cli::cli_abort("Pattern {.val {pattern}} not found in directory {.path {imgpath}}.")
    }

    allowed_ext <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF")
    if (!all(extensions %in% allowed_ext)) {
      cli::cli_abort("Allowed extensions are {.val {allowed_ext}}.")
    }

    old_opt <- options(cli.progress_bar_style = "bar")
    on.exit(options(old_opt), add = TRUE)

    args_list <- list(
      foreground = foreground, background = background, pick_palettes = pick_palettes, resize = resize, trim = trim, fill_hull = fill_hull, threshold = threshold,
      erode = erode, dilate = dilate, opening = opening, closing = closing, filter = filter, filter_order = filter_order,
      tolerance = tolerance, extension = extension, randomize = randomize, nrows = nrows, plot = plot, show_original = show_original, show_chull = show_chull,
      show_contour = show_contour, show_bbox = show_bbox, contour_col = contour_col, contour_size = contour_size, show_lw = show_lw,
      show_background = show_background, show_segmentation = show_segmentation, marker = marker, marker_col = marker_col, marker_size = marker_size, save_image = save_image,
      prefix = prefix, dir_original = dir_original, dir_processed = dir_processed, verbose = verbose, col_background = col_background,
      col_foreground = col_foreground, lower_noise = lower_noise, lower_size = lower_size, upper_size = upper_size, topn_lower = topn_lower, topn_upper = topn_upper,
      lower_eccent = lower_eccent, upper_eccent = upper_eccent, lower_circ = lower_circ, upper_circ = upper_circ,
      ab_angles = ab_angles, ab_angles_percentiles = ab_angles_percentiles, width_at = width_at, width_at_percentiles = width_at_percentiles,
      return_mask = return_mask, pcv = pcv, pcv_niter = pcv_niter, object_index = object_index, max_pixels = max_pixels, return_exact = return_exact, area_mode = area_mode,
      watershed = watershed, reference = reference, reference_area = reference_area, reference_larger = reference_larger, reference_smaller = reference_smaller,
      veins = veins, sigma_veins = sigma_veins, veins_thinning = veins_thinning, veins_rel_erode = veins_rel_erode,
      haralick = haralick, har_nbins = har_nbins, har_scales = har_scales, har_band = har_band, smooth = smooth, contours = contours,
      index = index, r = r, g = g, b = b, re = re, nir = nir, swir = swir, invert = invert, k = k, windowsize = windowsize,
      pixel_level_index = pixel_level_index, efourier = efourier, nharm = nharm, back_fore_index = back_fore_index, viewer = viewer
    )

    if (parallel == TRUE) {
      if (!requireNamespace("mirai", quietly = TRUE)) {
        cli::cli_abort("Package {.val mirai} is required for parallel processing.")
      }
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.3), workers)
      pkg_root <- tryCatch(rprojroot::find_package_root_file(), error = function(e) getwd())
      is_dev_mode <- file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE)

      mirai::daemons(nworkers)
      if (is_dev_mode) {
        mirai::everywhere({
          .libPaths(lp)
          pkgload::load_all(p_dir, quiet = TRUE, helpers = FALSE)
        }, lp = .libPaths(), p_dir = pkg_root)
      } else {
        mirai::everywhere({
          .libPaths(lp)
          library(pliman)
        }, lp = .libPaths())
      }
      on.exit(mirai::daemons(0), add = TRUE)

      jobs <- vector("list", length(plants))
      for (i in seq_along(plants)) {
        args_i <- c(list(img = plants[i]), args_list)
        jobs[[i]] <- mirai::mirai(do.call(analyze_objects3, args_i), args_i = args_i)
      }
      results <- vector("list", length(plants))
      for (i in seq_along(plants)) {
        results[[i]] <- mirai::call_mirai(jobs[[i]])$data
      }
      names(results) <- names_plant
      class(results) <- "anal_obj_ls"
      invisible(results)
    } else {
      results <- vector("list", length(plants))
      if (verbose) pb <- progress_bar$new(total = length(plants))
      for (i in seq_along(plants)) {
        if (verbose) pb$tick()
        args_i <- c(list(img = plants[i]), args_list)
        results[[i]] <- suppressMessages(do.call(analyze_objects3, args_i))
      }
      names(results) <- names_plant
      class(results) <- "anal_obj_ls"
      invisible(results)
    }
  }
}
