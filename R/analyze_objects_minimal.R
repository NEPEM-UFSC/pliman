#'Analyzes objects in an image
#'
#' A lighter option to [analyze_objects()]
#'
#' @export
#' @name analyze_objects_minimal
#' @inheritParams analyze_objects
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' obj <- analyze_objects(img)
#' obj$statistics
#'
#' }
#'
analyze_objects_minimal <- function(img,
                                    segment_objects = TRUE,
                                    reference = FALSE,
                                    reference_area = NULL,
                                    back_fore_index = "R/(G/B)",
                                    fore_ref_index = "B-R",
                                    reference_larger = FALSE,
                                    reference_smaller = FALSE,
                                    pattern = NULL,
                                    parallel = FALSE,
                                    workers = NULL,
                                    watershed = TRUE,
                                    fill_hull = FALSE,
                                    opening = FALSE,
                                    closing = FALSE,
                                    filter = FALSE,
                                    erode = FALSE,
                                    dilate = FALSE,
                                    invert = FALSE,
                                    object_size = "medium",
                                    index = "NB",
                                    r = 1,
                                    g = 2,
                                    b = 3,
                                    re = 4,
                                    nir = 5,
                                    threshold = "Otsu",
                                    tolerance = NULL,
                                    extension = NULL,
                                    lower_noise = 0.10,
                                    lower_size = NULL,
                                    upper_size = NULL,
                                    topn_lower = NULL,
                                    topn_upper = NULL,
                                    lower_eccent = NULL,
                                    upper_eccent = NULL,
                                    lower_circ = NULL,
                                    upper_circ = NULL,
                                    plot = TRUE,
                                    show_original = TRUE,
                                    show_contour = TRUE,
                                    contour_col = "red",
                                    contour_size = 1,
                                    col_foreground = NULL,
                                    col_background = NULL,
                                    marker = FALSE,
                                    marker_col = NULL,
                                    marker_size = NULL,
                                    save_image = FALSE,
                                    prefix = "proc_",
                                    dir_original = NULL,
                                    dir_processed = NULL,
                                    verbose = TRUE){
  check_ebi()
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)
  if (!object_size %in% c("small", "medium", "large", "elarge")) {
    cli::cli_abort("{.arg object_size} must be one of {.val small}, {.val medium}, {.val large}, or {.val elarge}.")
  }

  if (!missing(img) && !missing(pattern)) {
    cli::cli_abort("Only one of {.arg img} or {.arg pattern} can be used.")
  }

  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  help_count <-
    function(img, fill_hull, threshold, opening, closing, filter, erode, dilate, tolerance, extension,  plot,
             show_original,  marker, marker_col, marker_size,
             save_image, prefix, dir_original, dir_processed, verbose,
             col_background, col_foreground, lower_noise){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      # when reference is not used
      if(isFALSE(reference)){
        if(isTRUE(segment_objects)){
          img2 <- help_binary(img,
                              index = index,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              invert = invert,
                              fill_hull = fill_hull,
                              threshold = threshold,
                              opening = opening,
                              closing = closing,
                              filter = filter,
                              erode = erode,
                              dilate = dilate,
                              resize = FALSE)
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img2)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(img2),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(img2)
          }
        } else{
          img2 <- img[,,1]
          img2[img2@.Data == 0 | img2@.Data != 0] <- TRUE
          nmask <- EBImage::bwlabel(img2)
        }

        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
        if(isTRUE(fill_hull)){
          nmask <- EBImage::fillHull(nmask)
        }
        shape <- compute_measures_minimal(mask = nmask)
        object_contour <- shape$cont
        shape <- shape$shape

      } else{
        # when reference is used
        if(is.null(reference_area)){
          cli::cli_abort("A known area must be declared when a template is used.")

        }
        if(isFALSE(reference_larger) & isFALSE(reference_smaller)){
          # segment back and fore
          if(!isFALSE(invert)){
            invert1 <- ifelse(length(invert) == 1, invert, invert[1])
          } else{
            invert1 <- FALSE
          }
          img_bf <-
            help_binary(img,
                        threshold = threshold,
                        index = back_fore_index,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        erode = erode,
                        dilate = dilate,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        invert = invert1,
                        fill_hull = fill_hull)
          img3 <- img
          img3@.Data[,,1][which(img_bf != 1)] <- 2
          img3@.Data[,,2][which(img_bf != 1)] <- 2
          img3@.Data[,,3][which(img_bf != 1)] <- 2
          ID <-  which(img_bf == 1) # IDs for foreground
          ID2 <- which(img_bf == 0) # IDs for background
          # segment fore and ref
          if(!isFALSE(invert)){
            invert2 <- ifelse(length(invert) == 1, invert, invert[2])
          } else{
            invert2 <- FALSE
          }
          img4 <-
            help_binary(img3,
                        threshold = threshold,
                        index = fore_ref_index,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        erode = erode,
                        dilate = dilate,
                        invert = invert2)
          mask <- img_bf
          pix_ref <- which(img4 != 1)
          img@.Data[,,1][pix_ref] <- 1
          img@.Data[,,2][pix_ref] <- 0
          img@.Data[,,3][pix_ref] <- 0
          npix_ref <- length(pix_ref)
          mask[pix_ref] <- 0
          if(is.numeric(filter) & filter > 1){
            mask <- EBImage::medianFilter(mask, size = filter)
          }
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(mask),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(mask)
          }

          shape <- compute_measures_minimal(mask = nmask)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape

          # correct measures based on the area of the reference object
          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side^2
          shape[5:8] <- apply(shape[5:8], 2, function(x){
            x * px_side
          })
        } else{
          # correct the measures based on larger or smaller objects
          mask <-
            help_binary(img,
                        threshold = threshold,
                        index = index,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        erode = erode,
                        dilate = dilate,
                        invert = invert,
                        fill_hull = fill_hull)
          ID <-  which(mask == 1) # IDs for foreground
          ID2 <- which(mask == 0) # IDs for background
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(mask)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(mask),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(mask)
          }

          shape <- compute_measures_minimal(mask = nmask)

          object_contour <- shape$cont
          shape <- shape$shape

          if(isTRUE(reference_larger)){
            id_ref <- which.max(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
          } else{
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
            id_ref <- which.min(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
          }

          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side ^ 2
          shape[5:8] <- apply(shape[5:8], 2, function(x){
            x * px_side
          })
        }
      }


      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        cli::cli_abort("Only one of {.arg lower_*} or {.arg topn_*} can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$area > lower_size, ],
             shape <- shape[shape$area > mean(shape$area) * lower_noise, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$area, decreasing = TRUE),][1:topn_upper,]
      }
      if(!is.null(lower_eccent)){
        shape <- shape[shape$eccentricity > lower_eccent, ]
      }
      if(!is.null(upper_eccent)){
        shape <- shape[shape$eccentricity < upper_eccent, ]
      }
      if(!is.null(lower_circ)){
        shape <- shape[shape$circularity > lower_circ, ]
      }
      if(!is.null(upper_circ)){
        shape <- shape[shape$circularity < upper_circ, ]
      }
      object_contour <- object_contour[as.character(shape$id)]



      stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area"),
                          value = c(length(shape$area),
                                    min(shape$area),
                                    mean(shape$area),
                                    max(shape$area)))
      results <- list(results = shape,
                      statistics = stats)
      class(results) <- "anal_obj_minimal"
      if(plot == TRUE | save_image == TRUE){
        backg <- !is.null(col_background)
        # color for background
        if (is.null(col_background)){
          col_background <- col2rgb("white") / 255
        } else{
          ifelse(is.character(col_background),
                 col_background <- col2rgb(col_background) / 255,
                 col_background <- col_background / 255)
        }
        # color for lesions
        if (is.null(col_foreground)){
          col_foreground <- col2rgb("black") / 255
        } else{
          ifelse(is.character(col_foreground),
                 col_foreground <- col2rgb(col_foreground) / 255,
                 col_foreground <- col_foreground / 255)
        }

        if(show_original == TRUE){
          im2 <- img[,,1:3]
          EBImage::colorMode(im2) <- "Color"
          if(backg){
            im3 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- col_background[3]
          }
        }
        show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
        marker <- ifelse(is.null(marker), "id", marker)
        if (!isFALSE(show_mark) && marker != "point" && !marker %in% colnames(shape)) {
          cli::cli_warn(c(
            "!" = "{.arg marker} must be one of: {.val {colnames(shape)}}.",
            "i" = "Defaulting to {.val 'id'} (object ID will be drawn)."
          ))
          marker <- "id"
        }

        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
        # correct the contour
        object_contour <- lapply(object_contour, function(x){
          x + 1
        })

        if(plot == TRUE){
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2] + 1,
                   shape[, 3] + 1,
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour)  & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2] + 1,
                     shape[, 3] + 1,
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
        }

        if(save_image == TRUE){
          if(dir.exists(diretorio_processada) == FALSE){
            dir.create(diretorio_processada, recursive = TRUE)
          }
          png(paste0(diretorio_processada, "/",
                     prefix,
                     name_ori, ".",
                     extens_ori),
              width = dim(im2@.Data)[1],
              height = dim(im2@.Data)[2])
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2] + 1,
                   shape[, 3] + 1,
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2] + 1,
                     shape[, 3] + 1,
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          dev.off()
        }
      }
      invisible(results)
    }

  if(missing(pattern)){
    if(verbose){
      cli::cli_progress_step(
        msg = "Processing a single image. Please, wait.",
        msg_done = "Image {.emph Successfully} analyzed!",
        msg_failed = "Oops, something went wrong."
      )
    }
    help_count(img, fill_hull, threshold, opening, closing, filter, erode, dilate, tolerance, extension,  plot,
               show_original,  marker, marker_col, marker_size,
               save_image, prefix, dir_original, dir_processed, verbose,
               col_background, col_foreground, lower_noise)
  } else{
    if (pattern %in% as.character(0:9)) {
      old_pattern <- pattern
      pattern <- "^[0-9].*$"
      cli::cli_alert_info(
        "Numeric pattern {.val {old_pattern}} converted to {.val {pattern}}."
      )
    }

    # List files
    plants      <- list.files(path = diretorio_original, pattern = pattern)
    extensions  <- tolower(vapply(plants, tools::file_ext,   ""))
    names_plant <-        vapply(plants, tools::file_path_sans_ext, "")
    imgpath <- file.path(getwd(), sub('./', '', diretorio_original)) |> trunc_path(max_chars = 50)

    # Error if no matches
    if (length(plants) == 0) {
      cli::cli_abort(c(
        "x" = "Pattern {.val {pattern}} not found in {.path {diretorio_original}}",
        "i" = "Check your working directory: {.path {getwd()}}"
      ))
    }

    # Error on unsupported extensions
    bad_ext <- setdiff(extensions, c("png", "jpeg", "jpg", "tiff"))
    if (length(bad_ext) > 0) {
      cli::cli_abort(c(
        "x" = "Unsupported extension{?s}: {.val {unique(bad_ext)}} found.",
        "i" = "Allowed extensions are {.val png}, {.val jpeg}, {.val jpg}, {.val tiff}."
      ))
    }

    if(parallel == TRUE){
      init_time <- Sys.time()
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.3), workers)

      # Inicia workers persistentes do mirai
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0))

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Parallel processing using {nworkers} cores"),
          right = cli::col_blue("Started on  {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg       = "Processing {.val {length(names_plant)}} images found on {.path {imgpath}}. Please, wait.",
          msg_done  = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

      # Define a função para processar cada imagem
      process_image <- function(img) {
        help_count(
          img,
          fill_hull, threshold, opening, closing, filter, erode, dilate, tolerance, extension, plot,
          show_original, marker, marker_col, marker_size,
          save_image, prefix, dir_original, dir_processed, verbose,
          col_background, col_foreground, lower_noise
        )
      }

      # Executa paralelamente com barra de progresso
      results <- mirai::mirai_map(
        .x = names_plant,
        .f = process_image
      )[.progress]

    } else{
      cli::cli_rule(
        left = cli::col_blue("Analyzing {length(names_plant)} images"),
        right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%S')}")
      )
      cli::cli_alert_info("Directory: {.path {imgpath}}")
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta} | {.val {cli::pb_status}}",
        total = length(names_plant),
        clear = FALSE
      )

      results <- vector("list", length(names_plant))
      for (i in seq_along(names_plant)) {
        img_name <- names_plant[i]
        cli::cli_progress_update(status = names_plant[i])

        results[[i]] <- help_count(
          img = names_plant[i],
          fill_hull, threshold, opening, closing, filter, erode, dilate,
          tolerance, extension,  plot, show_original,  marker, marker_col, marker_size,
          save_image, prefix, dir_original, dir_processed, verbose,
          col_background, col_foreground, lower_noise
        )
      }
      cli::cli_progress_done()
    }


    ## bind the results
    names(results) <- names_plant

    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["statistics"]],
                          id =  names(results[i]))[,c(3, 1, 2)]
              })
      )

    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["results"]],
                          img =  names(results[i]))
              })
      )

    if("img" %in% colnames(results)){
      results <- results[, c(ncol(results), 1:ncol(results) - 1)]
    }
    nimages <- length(unique(stats$id))
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
      cli::boxx(glob_stat, header = "Global statistics ")  |> cat(sep = "\n")

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

      cli::cli_rule(
        left = cli::col_blue("Processing successfully finished"),
        right = cli::col_blue("on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
      )
    }

    invisible(
      structure(
        list(statistics = stats,
             count = stats[stats$stat == "n", c(1, 3)],
             results = results),
        class = "anal_obj_ls_minimal"
      )
    )
  }
}


#' @name analyze_objects_minimal
#' @inheritParams plot.analyze_objects
#' @method plot anal_obj_minimal
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'
#' img <- image_pliman("soy_green.jpg")
#' # Segment the foreground (grains) using the normalized blue index (NB, default)
#' # Shows the average value of the blue index in each object
#'
#' rgb <- analyze_objects_minimal(img)
#' # density of area
#' plot(rgb)
#'
#' # histogram of area
#' plot(rgb, type = "histogram") # or 'hist'
#' }
plot.anal_obj_minimal <- function(x,
                                  which = "measure",
                                  measure = "area",
                                  type = c("density", "histogram"),
                                  ...){
  if (!which %in% c("measure", "index")) {
    cli::cli_abort("{.arg which} must be one of {.val measure} or {.val index}.")
  }

  nam <- colnames(x$results)
  if (!measure %in% nam) {
    cli::cli_abort(c(
      "x" = "Measure {.val {measure}} not available in {.arg x}.",
      "i" = "Try one of {.val {paste(nam, collapse = \", \")}}."
    ))
  }

  temp <- x$results[[measure]]
  types <- c("density", "histogram")
  matches <- grepl(type[1], types)
  type <- types[matches]
  if(type == "histogram"){
    hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
  } else{
    density_data <- density(temp)  # Calculate the density for the column
    plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
    points(x = temp, y = rep(0, length(temp)), col = "red")
  }
}
#' @name analyze_objects_minimal
#' @export
plot.anal_obj_ls_minimal <- function(x,
                                     which = "measure",
                                     measure = "area",
                                     type = c("density", "histogram"),
                                     ...){
  if (!which %in% c("measure", "index")) {
    cli::cli_abort("{.arg which} must be one of {.val measure} or {.val index}.")
  }

  nam <- colnames(x$results)
  if (!measure %in% nam) {
    cli::cli_abort(c(
      "x" = "Measure {.val {measure}} not available in {.arg x}.",
      "i" = "Try one of {.val {paste(nam, collapse = \", \")}}."
    ))
  }
  temp <- x$results[[measure]]
  types <- c("density", "histogram")
  matches <- grepl(type[1], types)
  type <- types[matches]
  if(type == "histogram"){
    hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
  } else{
    density_data <- density(temp)  # Calculate the density for the column
    plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
    points(x = temp, y = rep(0, length(temp)), col = "red")
  }
}


