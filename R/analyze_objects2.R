.index_csv_cache_env <- new.env(parent = emptyenv())
get_cached_indexes_df <- function() {
  if (is.null(.index_csv_cache_env$df)) {
    fpath <- system.file("indexes.csv", package = "pliman", mustWork = TRUE)
    .index_csv_cache_env$df <- read.csv(file = fpath, header = TRUE, sep = ";")
  }
  .index_csv_cache_env$df
}

#' Analyzes objects in an image (Refactored C++ Engine)
#'
#' @description
#' `analyze_objects2()` provides an ultra-fast, cache-optimized implementation of
#' [analyze_objects()] by executing core computations (vegetation index, thresholding,
#' morphological filtering, watershed segmentation, contour tracing, and 35+ shape measures)
#' directly in C++.
#'
#' @inheritParams analyze_objects
#' @param area_mode The method to compute object area. Default is `"contour"` (Shoelace area
#'   of the boundary polygon). If `"pixel"`, counts the exact number of labeled mask pixels
#'   belonging to the object (useful when objects contain internal holes/islands).
#' @return An object of class `anal_obj`. See [analyze_objects()] for details.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
analyze_objects2 <- function(img,
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

  t0 <- Sys.time()
  format_time <- function(t_start) {
    elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    if (elapsed_sec < 60) {
      sprintf("%.2fs", elapsed_sec)
    } else if (elapsed_sec < 3600) {
      mins <- floor(elapsed_sec / 60)
      secs <- round(elapsed_sec %% 60, 1)
      sprintf("%dm %.1fs", mins, secs)
    } else {
      hrs <- floor(elapsed_sec / 3600)
      rem_mins <- floor((elapsed_sec %% 3600) / 60)
      sprintf("%dh %dm", hrs, rem_mins)
    }
  }

  area_mode <- area_mode[[1]]
  if (!area_mode %in% c("contour", "pixel")) {
    cli::cli_abort("Argument {.arg area_mode} must be one of {.val contour} or {.val pixel}.")
  }

  check_filter_order(filter_order, verbose, erode, dilate, opening, closing, filter, fill_hull)
  if (!is.null(rel_size)) {
    lower_noise <- rel_size
  }
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)

  if (!missing(img) && !missing(pattern)) {
    cli::cli_abort("Only one of {.arg img} or {.arg pattern} can be used.")
  }

  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <- ifelse(grepl("[/\\]", dir_original), dir_original, paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <- ifelse(grepl("[/\\]", dir_processed), dir_processed, paste0("./", dir_processed))
  }

  help_count2 <- function(img, foreground = NULL, background = NULL, pick_palettes = FALSE, resize = FALSE, trim = FALSE, fill_hull = FALSE, threshold = "Otsu",
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
                          haralick = FALSE, har_nbins = 32, har_scales = 1, har_band = 1, smooth = FALSE,
                          index = "NB", r = 1, g = 2, b = 3, re = 4, nir = 5, swir = 6, invert = FALSE, k = 0.1, windowsize = NULL,
                          pixel_level_index = FALSE, efourier = FALSE, nharm = 10, back_fore_index = "R/(G/B)", viewer = get_pliman_viewer()) {

    if (is.null(dir_original)) {
      diretorio_original <- "./"
    } else {
      diretorio_original <- ifelse(grepl("[/\\]", dir_original), dir_original, paste0("./", dir_original))
    }
    if (is.null(dir_processed)) {
      diretorio_processada <- "./"
    } else {
      diretorio_processada <- ifelse(grepl("[/\\]", dir_processed), dir_processed, paste0("./", dir_processed))
    }

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

    if(trim != FALSE){
      if(!is.numeric(trim)) cli::cli_abort("Argument {.arg trim} must be numeric.")
      img <- image_trim(img, trim)
      im_dat <- if (inherits(img, "Image")) img@.Data else image_data(img)
    }
    if(resize != FALSE){
      if(!is.numeric(resize)) cli::cli_abort("Argument {.arg resize} must be numeric.")
      img <- image_resize(img, resize)
      im_dat <- if (inherits(img, "Image")) img@.Data else image_data(img)
    }

    # Setup watershed defaults matching image_watershed
    if (is.null(tolerance)) {
      tolerance <- 0.2
    }
    if (is.null(extension)) {
      extension <- 1
    }

    bin_precomputed <- NULL

    if (isFALSE(reference)) {
      if (isTRUE(pick_palettes)) {
        if (interactive()) {
          background <- pick_palette(img, r = 5, verbose = FALSE, palette = FALSE, plot = FALSE, col = "blue", external_device = FALSE, title = "Pick BACKGROUND colors", viewer = viewer)
          foreground <- pick_palette(img, r = 5, verbose = FALSE, palette = FALSE, plot = FALSE, col = "salmon", external_device = FALSE, title = "Pick FOREGROUND colors", viewer = viewer)
        }
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

    # Extract raw data / SEXP for C++ engine
    img_data_raw <- im_dat

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
                      cli::cli_abort("Invalid {.arg sensitivity} '{sens_str}'. Must be one of: 'low', 'medium', 'high', 'extreme', or a numeric value."))
      } else if (is.numeric(sensitivity)) {
        lvl <- as.double(sensitivity[1])
      } else {
        cli::cli_abort("{.arg sensitivity} must be a character string ('low', 'medium', 'high', 'extreme') or a numeric value.")
      }
      0.20 / (2.0 ^ (lvl - 1.0))
    } else {
      as.double(tolerance)
    }

    ext_val <- if (is.null(extension)) 1L else as.integer(extension)

    # Execute fast C++ pipeline
    cpp_res <- analyze_objects_cpp(
      img_sexp = img_data_raw,
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
      smooth = smooth_val
    )

    nmask <- cpp_res[["labels"]]
    ocont <- cpp_res[["contours"]]
    shape <- cpp_res[["shape"]]

    valid <- which(!is.na(shape$x))
    if (length(valid) < nrow(shape)) {
      shape <- shape[valid, , drop = FALSE]
      ocont <- ocont[valid]
    }
    names(ocont) <- shape$id

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

    ch <- if (!is.null(cpp_res[["chull"]])) {
      if (length(valid) < length(cpp_res[["chull"]])) cpp_res[["chull"]][valid] else cpp_res[["chull"]]
    } else {
      conv_hull(ocont)
    }
    if (!is.null(ch)) names(ch) <- shape$id

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
    if (!is.null(lower_size) && !is.null(topn_lower) || !is.null(upper_size) && !is.null(topn_upper)) {
      cli::cli_abort("Only one of {.arg lower_*} or {.arg topn_*} can be used.")
    }
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

    keep_idx <- match(shape$id, names(ocont))
    ocont <- ocont[keep_idx]
    if (!is.null(ch)) ch <- ch[keep_idx]
    if (nrow(shape) > 0) {
      filter_labels_cpp(nmask, shape$id)
    } else {
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

    if (isTRUE(width_at) && nrow(shape) > 0) {
      widths <- do.call(rbind, lapply(ocont, function(x) {
        x |> poly_align(plot = FALSE) |> poly_width_at(width_at_percentiles)
      })) |> as.data.frame() |> rownames_to_column("id")
      names(widths) <- c("id", paste0("width", width_at_percentiles))
    } else {
      widths <- NULL
    }

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

    shape_markers <- if (!is.null(object_index) && nrow(shape) > 0) cbind(shape, indexes) else shape

    if (plot == TRUE || save_image == TRUE) {
      backg <- !is.null(col_background)
      if (is.null(col_background)) {
        col_bg <- col2rgb("white")
      } else if (is.character(col_background)) {
        col_bg <- col2rgb(col_background)
      } else {
        col_bg <- col_background
      }
      if (is.null(col_foreground)) {
        col_fg <- col2rgb("gray")
      } else if (is.character(col_foreground)) {
        col_fg <- col2rgb(col_foreground)
      } else {
        col_fg <- col_foreground
      }
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
      if (!isFALSE(show_mark) && marker_val != "point" && !marker_val %in% colnames(shape_markers)) {
        cli::cli_warn("Accepted {.arg marker} values are: {.val {paste(colnames(shape_markers), collapse = \", \")}}. Drawing object id instead.")
        marker_val <- "id"
      }
      marker_col_val <- ifelse(is.null(marker_col), "white", marker_col)
      marker_sz <- ifelse(is.null(marker_size), 0.75, marker_size)

      if (plot == TRUE) {
        plot(im2, max_pixels = max_pixels)
        if (nrow(shape) > 0) {
          if (isTRUE(show_contour)) plot_contour(ocont, col = contour_col, lwd = contour_size)
          if (show_bbox) plot_bbox(ocont, col = contour_col)
          if (show_mark && marker_val != "point") {
            text(shape_markers[, 2], shape_markers[, 3], round(shape_markers[, marker_val], 3), col = marker_col_val, cex = marker_sz)
          } else if (show_mark && marker_val == "point") {
            points(shape_markers[, 2], shape_markers[, 3], col = marker_col_val, pch = 16, cex = marker_sz)
          }
          if (isTRUE(show_chull)) plot_contour(ch |> poly_close(), col = "black")
          if (isTRUE(show_lw)) plot_lw(results)
        }
      }

      if (save_image == TRUE) {
        if (!dir.exists(diretorio_processada)) {
          dir.create(diretorio_processada, recursive = TRUE)
        }
        png_file <- file.path(diretorio_processada, paste0(prefix, name_ori, ".", extens_ori))
        img_d <- dim(image_data(im2))
        tot_pix <- img_d[1] * img_d[2]
        if (!is.null(max_pixels) && tot_pix > max_pixels) {
          scale_factor <- sqrt(max_pixels / tot_pix)
          png_w <- round(img_d[1] * scale_factor)
          png_h <- round(img_d[2] * scale_factor)
        } else {
          png_w <- img_d[1]
          png_h <- img_d[2]
        }
        png(png_file, width = png_w, height = png_h)
        dev_num <- dev.cur()
        on.exit(if (dev.cur() == dev_num) dev.off(), add = TRUE)

        plot(im2, max_pixels = max_pixels)
        if (nrow(shape) > 0) {
          if (isTRUE(show_contour)) plot_contour(ocont, col = contour_col, lwd = contour_size)
          if (show_bbox) plot_bbox(ocont, col = contour_col)
          if (show_mark && marker_val != "point") {
            text(shape_markers[, 2], shape_markers[, 3], round(shape_markers[, marker_val], 3), col = marker_col_val, cex = marker_sz)
          } else if (show_mark && marker_val == "point") {
            points(shape_markers[, 2], shape_markers[, 3], col = marker_col_val, pch = 16, cex = marker_sz)
          }
          if (isTRUE(show_lw)) plot_lw(results)
        }
        footer_txt <- sprintf("N: %d  |  Area (px\u00b2)  mean: %s  |  min: %s  |  max: %s", n_obj, a_mean, a_min, a_max)
        usr <- par("usr")
        text(usr[1] + diff(usr[1:2]) * 0.01, usr[3] + diff(usr[3:4]) * 0.01, footer_txt, adj = c(0, 0), cex = 1, col = "#555555", font = 1, xpd = NA)
        dev.off()
      }
    }

    invisible(results)
  }

  if (missing(pattern)) {
    if (verbose) cli::cli_progress_step("{.pkg pliman} is processing the image. Please wait.")
    help_count2(
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
      haralick = haralick, har_nbins = har_nbins, har_scales = har_scales, har_band = har_band, smooth = smooth,
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
      haralick = haralick, har_nbins = har_nbins, har_scales = har_scales, har_band = har_band, smooth = smooth,
    index = index, r = r, g = g, b = b, re = re, nir = nir, swir = swir, invert = invert, k = k, windowsize = windowsize,
      pixel_level_index = pixel_level_index, efourier = efourier, nharm = nharm, back_fore_index = back_fore_index, viewer = viewer
    )

    if (parallel == TRUE) {
      if (!requireNamespace("mirai", quietly = TRUE)) {
        cli::cli_abort("Package {.val mirai} is required for parallel processing. Please install it with {.code install.packages('mirai')}.")
      }
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.3), workers)
      pkg_root <- tryCatch(rprojroot::find_package_root_file(), error = function(e) getwd())
      is_dev_mode <- file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE)

      daemon_status <- tryCatch({
        mirai::daemons(nworkers)
        if (is_dev_mode) {
          mirai::everywhere(
            {
              .libPaths(lp)
              pkgload::load_all(p_dir, quiet = TRUE, helpers = FALSE)
            },
            lp = .libPaths(),
            p_dir = pkg_root
          )
        } else {
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
        on.exit(mirai::daemons(0), add = TRUE)
        if (verbose) {
          cli::cli_rule(
            left = cli::col_blue("Parallel processing using {nworkers} cores"),
            right = cli::col_blue("Started on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
          )
        }

        results <- mirai::mirai_map(
          .x = names_plant,
          .f = help_count2,
          .args = args_list
        )[.progress]

        failed_imgs <- sapply(results, function(res) inherits(res, "errorValue") || !is.list(res) || is.null(res[["statistics"]]))
        if (any(failed_imgs)) {
          err_msg <- paste(sapply(which(failed_imgs), function(idx) {
            if (inherits(results[[idx]], "errorValue")) as.character(results[[idx]]) else "Worker failed to return results"
          }), collapse = "; ")
          cli::cli_abort("Parallel processing failed for image(s) {.val {names_plant[failed_imgs]}}: {err_msg}")
        }
      } else {
        parallel <- FALSE
      }
    }

    if (parallel == FALSE) {
      if (verbose) {
        cli::cli_rule(
          left = cli::col_blue("Analyzing {length(names_plant)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%S')}")
        )

        cli::cli_progress_bar(
          format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
          total = length(names_plant),
          clear = FALSE
        )
      }
      results <- vector("list", length(names_plant))
      for (i in seq_along(names_plant)) {
        if (verbose) cli::cli_progress_update()
        results[[i]] <- do.call(help_count2, c(list(img = names_plant[i]), args_list))
      }
      if (verbose) {
        cli::cli_progress_done()
      }
    }

    ## bind the results
    names(results) <- names_plant
    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                st <- results[[i]][["statistics"]]
                st$img <- names(results[i])
                st[, c(3, 1, 2)]
              })
      )

    if (!is.null(object_index)) {
      obj_rgb_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["object_rgb"]]
        if (!is.null(df) && nrow(df) > 0) {
          df$img <- names(results[i])
          df[, c(ncol(df), 1:(ncol(df) - 1))]
        } else NULL
      })
      obj_rgb <- if (any(!sapply(obj_rgb_list, is.null))) do.call(rbind, obj_rgb_list) else NULL

      obj_idx_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["object_index"]]
        if (!is.null(df) && nrow(df) > 0) {
          df$img <- names(results[i])
          df[, c(ncol(df), 1:(ncol(df) - 1))]
        } else NULL
      })
      object_index <- if (any(!sapply(obj_idx_list, is.null))) do.call(rbind, obj_idx_list) else NULL
    } else {
      obj_rgb <- NULL
      object_index <- NULL
    }

    if (!isFALSE(efourier)) {
      bind_ef <- function(key, col2_name = "id") {
        lst <- lapply(seq_along(results), function(i) {
          df <- results[[i]][[key]]
          if (!is.null(df) && nrow(df) > 0) {
            df$img <- names(results[i])
            df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
            if (!is.null(col2_name) && ncol(df) >= 2) names(df)[2] <- col2_name
            df
          } else NULL
        })
        if (any(!sapply(lst, is.null))) do.call(rbind, lst) else NULL
      }
      efourier <- bind_ef("efourier")
      efourier_norm <- bind_ef("efourier_norm")
      efourier_error <- bind_ef("efourier_error")
      efourier_power <- bind_ef("efourier_power")
      efourier_minharm <- bind_ef("efourier_minharm")
    } else {
      efourier <- efourier_norm <- efourier_error <- efourier_power <- efourier_minharm <- NULL
    }

    if (isTRUE(veins)) {
      veins_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["veins"]]
        if (!is.null(df) && nrow(df) > 0) {
          df$img <- names(results[i])
          df[, c(ncol(df), 1:(ncol(df) - 1))]
        } else NULL
      })
      veins <- if (any(!sapply(veins_list, is.null))) do.call(rbind, veins_list) else NULL
    } else {
      veins <- NULL
    }

    if (isTRUE(ab_angles)) {
      angles_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["angles"]]
        if (!is.null(df) && nrow(df) > 0) {
          df$img <- names(results[i])
          df[, c(ncol(df), 1:(ncol(df) - 1))]
        } else NULL
      })
      angles <- if (any(!sapply(angles_list, is.null))) do.call(rbind, angles_list) else NULL
    } else {
      angles <- NULL
    }

    if (isTRUE(width_at)) {
      width_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["width_at"]]
        if (!is.null(df) && nrow(df) > 0) {
          df$img <- names(results[i])
          df[, c(ncol(df), 1:(ncol(df) - 1))]
        } else NULL
      })
      width_at <- if (any(!sapply(width_list, is.null))) do.call(rbind, width_list) else NULL
    } else {
      width_at <- NULL
    }

    if (isTRUE(pcv)) {
      pcv_list <- lapply(seq_along(results), function(i) {
        df <- results[[i]][["pcv"]]
        if (!is.null(df) && length(df) > 0) {
          data.frame(img = names(results[i]), pcv = df)
        } else NULL
      })
      pcv <- if (any(!sapply(pcv_list, is.null))) do.call(rbind, pcv_list) else NULL
    } else {
      pcv <- NULL
    }

    res_list <- lapply(seq_along(results), function(i) {
      df <- results[[i]][["results"]]
      if (nrow(df) > 0) {
        df$img <- names(results[i])
      } else {
        df$img <- character(0)
      }
      df
    })
    results <- do.call(rbind, res_list)

    if ("img" %in% colnames(results) && nrow(results) > 0) {
      results <- results[, c("img", setdiff(colnames(results), "img"))]
    }
    nimages <- length(unique(stats$img))
    n_img <-
      results |>
      dplyr::group_by(img) |>
      dplyr::summarise(
        n = dplyr::n(),
        area_mean = mean(area, na.rm = TRUE),
        area_min = min(area, na.rm = TRUE),
        area_max = max(area, na.rm = TRUE),
        area_sum = sum(area, na.rm = TRUE),
        area_sd = sd(area, na.rm = TRUE)
      )

    if(verbose == TRUE){
      average_n <- mean(n_img$n)
      min_n <- min(n_img$n)
      max_n <- max(n_img$n)
      average_area <- mean(n_img$area_mean)
      min_area <- min(n_img$area_max)
      max_area <- max(n_img$area_min)

      # Global statistics
      glob_stat <- cli::ansi_columns(
        paste(
          c(
            "Total objects:",
            "Total area:",
            "Overall mean area:",
            "Overall SD:",
            "Min area:",
            "Max area:"
          ),
          c(
            sum(n_img$n),
            round(sum(n_img$area_sum, na.rm = TRUE), 2),
            round(mean(results$area, na.rm = TRUE), 2),
            round(sd(results$area, na.rm = TRUE), 2),
            round(min(results$area, na.rm = TRUE), 2),
            round(max(results$area, na.rm = TRUE), 2)
          )
        ),
        width = 60,
        fill = "rows",
        align = "left",
        sep = "",
        max_cols = 2
      )
      cli::boxx(glob_stat, header = "Global statistics ") |> cat(sep = "\n")

      cross_imgstat <-
        cli::ansi_columns(
          paste(
            c(
              "Avg objects:",
              "Avg sum area:",
              "Min objects:",
              "Max objects:",
              "Avg area:",
              "Avg SD area:",
              "Min mean area:",
              "Max mean area:"
            ),
            c(
              round(mean(n_img$n), 2),
              round(mean(n_img$area_sum, na.rm = TRUE), 2),
              min(n_img$n),
              max(n_img$n),
              round(mean(n_img$area_mean, na.rm = TRUE), 2),
              round(mean(n_img$area_sd, na.rm = TRUE), 2),
              round(min(n_img$area_mean, na.rm = TRUE), 2),
              round(max(n_img$area_mean, na.rm = TRUE), 2)
            )
          ),
          width = 60,
          fill = "rows",
          align = "left",
          sep = "",
          max_cols = 2
        )

      cli::boxx(cross_imgstat,
                header = "Across-image statistics (per-image averages)",
                footer = paste0("Based on ", nimages, " images")) |>
        cat(sep = "\n")

      cli::cli_progress_done()

      time_str <- format_time(t0)
      cli::cli_rule(
        left = cli::col_green(paste0("\u2714 Processing successfully finished in ", time_str)),
        right = cli::col_blue(format(Sys.time(), format = "%Y-%m-%d | %H:%M:%S"))
      )
    }

    invisible(
      structure(
        list(statistics = n_img,
             count = stats[stats$stat == "n", c(1, 3)],
             results = results,
             obj_rgb = obj_rgb,
             object_index = object_index,
             efourier = efourier,
             efourier_norm = efourier_norm,
             efourier_error = efourier_error,
             efourier_minharm = efourier_minharm,
             veins = veins,
             angles = angles,
             width_at = width_at,
             pcv = pcv),
        class = "anal_obj_ls"
      )
    )
  }
}
