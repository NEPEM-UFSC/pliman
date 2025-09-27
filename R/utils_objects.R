#' Utilities for working with image objects
#'
#' * `object_id()` get the object identification in an image.
#' * `object_coord()` get the object coordinates and (optionally) draw a
#' bounding rectangle around multiple objects in an image.
#' * `object_contour()` returns the coordinates (`x` and `y`) for the contours
#' of each object in the image.
#' * `object_isolate()` isolates an object from an image.
#' @name utils_objects
#'
#' @inheritParams analyze_objects
#' @param img An image of class `Image` or a list of `Image` objects.
#' @param center If `TRUE` returns the object contours centered on the origin.
#' @param id
#' * For `object_coord()`, a vector (or scalar) of object `id` to compute the
#' bounding rectangle. Object ids can be obtained with [object_id()]. Set `id =
#' "all"` to compute the coordinates for all objects in the image. If `id =
#' NULL` (default) a bounding rectangle is drawn including all the objects.
#' * For `object_isolate()`, a scalar that identifies the object to be extracted.
#'
#' @param dir_original The directory containing the original images. Defaults
#'    to `NULL`, which means that the current working directory will be
#'    considered.
#' @param index The index to produce a binary image used to compute bounding
#'   rectangle coordinates. See [image_binary()] for more details.
#' @param invert Inverts the binary image, if desired. Defaults to `FALSE`.
#' @param opening,closing,filter **Morphological operations (brush size)**
#'  * `opening` performs an erosion followed by a dilation. This helps to
#'   remove small objects while preserving the shape and size of larger objects.
#'  * `closing` performs a dilatation followed by an erosion. This helps to
#'   fill small holes while preserving the shape and size of larger objects.
#'  * `filter` performs median filtering in the binary image. Provide a positive
#'  integer > 1 to indicate the size of the median filtering. Higher values are
#'  more efficient to remove noise in the background but can dramatically impact
#'  the perimeter of objects, mainly for irregular perimeters such as leaves
#'  with serrated edges.
#'
#'   Hierarchically, the operations are performed as opening > closing > filter.
#'   The value declared in each argument will define the brush size.
#'@param smooth whether the object contours should be smoothed with
#'  [poly_smooth()]. Defaults to `FALSE`. To smooth use a numeric value
#'  indicating the number of interactions used to smooth the contours.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param watershed If `TRUE` (default) performs watershed-based object
#'   detection. This will detect objects even when they are touching one other.
#'   If `FALSE`, all pixels for each connected set of foreground pixels are set
#'   to a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param threshold By default (`threshold = "Otsu"`), a threshold value based
#'   on Otsu's method is used to reduce the grayscale image to a binary image.
#'   If a numeric value is informed, this value will be used as a threshold.
#'   Inform any non-numeric value different than "Otsu" to iteratively chosen
#'   the threshold based on a raster plot showing pixel intensity of the index.
#' @param edge The number of pixels in the edge of the bounding rectangle.
#'   Defaults to `2`.
#' @param extension,tolerance,object_size Controls the watershed segmentation of
#'   objects in the image. See [analyze_objects()] for more details.
#' @param plot Shows the image with bounding rectangles? Defaults to
#'   `TRUE`.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 50% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param ...
#' * For `object_isolate()`, further arguments passed on to [object_coord()].
#' * For `object_id()`, further arguments passed on to [analyze_objects()].
#' @return
#' * `object_id()` An image of class `"Image"` containing the object's
#' identification.
#' * `object_coord()` A list with the coordinates for the bounding rectangles.
#' If `id = "all"` or a numeric vector, a list with a vector of coordinates is
#' returned.
#' * `object_isolate()` An image of class `"Image"` containing the isolated
#' object.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg")
#' # Get the object's (leaves) identification
#' object_id(img)
#'
#' # Get the coordinates and draw a bounding rectangle around leaves 1 and 3
#' object_coord(img, id = c(1, 3))
#'
#' # Isolate leaf 3
#' isolated <- object_isolate(img, id = 3)
#' plot(isolated)
#'
#' }
object_coord <- function(img,
                         id =  NULL,
                         index = "NB",
                         watershed = TRUE,
                         invert = FALSE,
                         opening = FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         edge = 2,
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         parallel = FALSE,
                         workers = NULL,
                         plot = TRUE,
                         verbose = TRUE){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.cls Image}")
    }
    if(parallel == TRUE){
      # decide number of workers
      nworkers <- ifelse(is.null(workers),
                         trunc(parallel::detectCores() * 0.5),
                         workers)

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # CLI info + progress step
      if (verbose) {
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "All object coordinates extracted.",
          msg_failed = "Object coordinate extraction failed."
        )
      }

      # run object_coord in parallel
      results <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          pliman::object_coord(
            im, id, index, invert,
            fill_hull, threshold, edge, extension, tolerance,
            object_size, plot
          )
        }
      )[.progress]

    } else{
      lapply(img, object_coord, id, index, invert, fill_hull, threshold,
             edge, extension, tolerance, object_size, plot)
    }
  } else{
    img2 <- help_binary(img,
                        index = index,
                        invert = invert,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        fill_hull = fill_hull,
                        threshold = threshold)
    if(is.null(id)){
      data_mask <- img2@.Data
      coord <- t(as.matrix(bounding_box(data_mask, edge)))
      colnames(coord) <- c("xleft", "xright", "ybottom", "ytop")
      if(plot == TRUE){
        plot(img)
        rect(xleft = coord[1],
             xright = coord[2],
             ybottom = coord[3],
             ytop = coord[4])
      }
    } else{
      if(isTRUE(watershed)){
        res <- length(img2)
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
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
      data_mask <- nmask@.Data
      ifelse(id == "all",
             ids <- 1:max(data_mask),
             ids <- id)
      list_mask <- list()
      for (i in ids) {
        temp <- data_mask
        temp[which(data_mask != i)] <- FALSE
        list_mask[[i]] <- temp
      }
      list_mask <- list_mask[ids]
      coord <- t(sapply(list_mask, bounding_box, edge))
      colnames(coord) <- c("xleft", "xright", "ybottom", "ytop")
      if(plot == TRUE){
        plot(img)
        rect(xleft = coord[,1],
             xright = coord[,2],
             ybottom = coord[,3],
             ytop = coord[,4])
      }
    }
    invisible(coord)
  }
}
#' @name utils_objects
#' @inheritParams analyze_objects
#' @export
#'

object_contour <- function(img,
                           pattern = NULL,
                           dir_original = NULL,
                           center =  FALSE,
                           index = "NB",
                           invert = FALSE,
                           opening = FALSE,
                           closing = FALSE,
                           filter = FALSE,
                           fill_hull = FALSE,
                           smooth = FALSE,
                           threshold = "Otsu",
                           watershed = TRUE,
                           extension = NULL,
                           tolerance = NULL,
                           object_size = "medium",
                           parallel = FALSE,
                           workers = NULL,
                           plot = TRUE,
                           verbose = TRUE){
  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }

  if(is.null(pattern) && inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.cls Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.5), workers)
      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # inform user
      cli::cli_alert_info("Image processing using multiple sessions ({.val {nworkers}}). Please wait.")
      cli::cli_progress_step(
        msg        = "Processing {.val {length(img)}} images in parallel...",
        msg_done   = "Object contour extraction complete.",
        msg_failed = "Object contour extraction failed."
      )

      # run object_contour in parallel
      results <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          pliman::object_contour(
            im,
            pattern       = pattern,
            dir_original  = dir_original,
            center        = center,
            index         = index,
            invert        = invert,
            opening       = opening,
            closing       = closing,
            filter        = filter,
            fill_hull     = fill_hull,
            smooth        = smooth,
            threshold     = threshold,
            watershed     = watershed,
            extension     = extension,
            tolerance     = tolerance,
            object_size   = object_size,
            plot          = plot
          )
        }
      )[.progress]
    } else{
      lapply(img, object_contour, pattern, dir_original, center, index, invert, opening, closing, filter, fill_hull, smooth, threshold,
             watershed, extension, tolerance, object_size, plot = plot)
    }
  } else{
    if(is.null(pattern)){
      img2 <- help_binary(img,
                          index = index,
                          invert = invert,
                          opening = opening,
                          closing = closing,
                          filter = filter,
                          fill_hull = fill_hull,
                          threshold = threshold)
      if(isTRUE(watershed)){
        res <- length(img2)
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
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
      contour <- EBImage::ocontour(nmask)
      contour <- lapply(contour, function(x){
        x + 1
      })

      if(isTRUE(center)){
        contour <-
          lapply(contour, function(x){
            transform(x,
                      X1 = X1 - mean(X1),
                      X2 = X2 - mean(X2))
          })
      }
      dims <- sapply(contour, function(x){dim(x)[1]})
      contour <- contour[which(dims > mean(dims * 0.1))]
      if(is.numeric(smooth) & smooth > 0){
        contour <- poly_smooth(contour, niter = smooth, plot = FALSE) |> poly_close()
      }
      if(isTRUE(plot)){
        if(isTRUE(center)){
          plot_polygon(contour)
        } else{
          plot(img)
          plot_contour(contour, col = "red")
        }
      }
      invisible(contour)
    } else{
      if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
        pattern <- "^[0-9].*$"
      }
      plants <- list.files(pattern = pattern, diretorio_original)
      extensions <- as.character(sapply(plants, file_extension))
      names_plant <- as.character(sapply(plants, file_name))
      if (length(grep(pattern, names_plant)) == 0) {
        cli::cli_abort("Pattern {.val {pattern}} not found in directory {.path {file.path(getwd(), sub('.', '', diretorio_original))}}.")
      }

      allowed_ext <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF")

      if (!all(extensions %in% allowed_ext)) {
        cli::cli_abort(c(
          "!" = "Some image files have unsupported extensions.",
          "i" = "Allowed extensions are: {.val {tolower(unique(allowed_ext))}}"
        ))
      }

      help_contour <- function(img){
        img <- image_import(img)
        img2 <- help_binary(img,
                            index = index,
                            invert = invert,
                            opening = opening,
                            closing = closing,
                            filter = filter,
                            fill_hull = fill_hull,
                            threshold = threshold)
        if(isTRUE(watershed)){
          res <- length(img2)
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
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
        contour <- EBImage::ocontour(nmask)
        if(isTRUE(center)){
          contour <-
            lapply(contour, function(x){
              transform(x,
                        X1 = X1 - mean(X1),
                        X2 = X2 - mean(X2))
            })
        }
        dims <- sapply(contour, function(x){dim(x)[1]})
        contour[which(dims > mean(dims * 0.1))]
      }


      if(parallel == TRUE){
        # decide number of workers
        nworkers <- ifelse(is.null(workers),
                           trunc(parallel::detectCores() * 0.5),
                           workers)

        # start mirai daemons
        mirai::daemons(nworkers)
        on.exit(mirai::daemons(0), add = TRUE)

        # CLI header + progress step
        if (verbose) {
          cli::cli_rule(
            left  = cli::col_blue("Parallel processing using {.val {nworkers}} cores"),
            right = cli::col_blue("Started on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
          )
          cli::cli_progress_step(
            msg        = "Processing {.val {length(plants)}} images using {.val {nworkers}} cores...",
            msg_done   = "Batch processing finished",
            msg_failed = "Oops, something went wrong."
          )
        }

        # run help_contour in parallel
        results <- mirai::mirai_map(
          .x = plants,
          .f = help_contour
        )[.progress]


      } else{
        if(verbose){
          cli::cli_rule(
            left = cli::col_blue("Analyzing {.val {length(names_plant)}} images found on {.path {diretorio_original}}."),
            right = cli::col_blue("Start on {.val {format(Sys.time(), format = '%Y-%m-%d - %H:%M:%OS0')}}")
          )
          cli::cli_progress_bar(
            format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [ETA:{cli::pb_eta}] | Current: {.val {cli::pb_status}}",
            total = length(names_plant),
            clear = FALSE
          )
        }
        results <- list()
        for (i in seq_along(plants)) {
          cli::cli_progress_update(status = plants[i])
          results[[i]] <- help_contour(img = plants[i])
        }


      }
      names(results) <- plants
      invisible(results)
    }
  }
}



#' @name utils_objects
#' @export
object_isolate <- function(img,
                           id = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           ...){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      # decide number of workers
      nworkers <- ifelse(is.null(workers),
                         trunc(parallel::detectCores() * 0.5),
                         workers)

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # CLI header + progress
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Isolating objects in {.val {length(img)}} images"),
          right = cli::col_blue("Using {.val {nworkers}} workers")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Object isolation complete.",
          msg_failed = "Object isolation failed."
        )
      }

      # run object_isolate in parallel
      results <- mirai::mirai_map(
        .x = img,
        .f = function(im){pliman::object_isolate(im, id)}
      )[.progress]

    } else{
      lapply(img, object_isolate, id, ...)
    }
  } else{
    coord <- object_coord(img,
                          id = id,
                          plot = FALSE,
                          ...)
    segmented <- img[coord[1]:coord[2],
                     coord[3]:coord[4],
                     1:3]
    invisible(segmented)
  }
}
#' @name utils_objects
#' @export
object_id <- function(img,
                      parallel = FALSE,
                      workers = NULL,
                      verbose = TRUE,
                      ...){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      # decide number of workers
      nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.5) else workers

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # CLI header + progress step
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Extracting object IDs from {.val {length(img)}} images"),
          right = cli::col_blue("Using {.val {nworkers}} workers")
        )
        cli::cli_progress_step(
          msg      = "Processing images in parallel...",
          msg_done = "Object ID extraction complete.",
          msg_failed = "Object ID extraction failed."
        )
      }

      # run object_id in parallel
      results <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          pliman::object_id(im, ...)
        }
      )[.progress]
    } else{
      lapply(img, object_id, ...)
    }
  } else{
    analyze_objects(img, verbose = FALSE, marker = "id", ...)
  }
}




#' Splits objects from an image into multiple images
#'
#' Using threshold-based segmentation, objects are first isolated from
#' background. Then, a new image is created for each single object. A list of
#' images is returned.
#'
#' @inheritParams analyze_objects
#' @param lower_size Plant images often contain dirt and dust. To prevent dust from
#'   affecting the image analysis, objects with lesser than 10% of the mean of all objects
#'   are removed. Set `lower_limit = 0` to keep all the objects.
#' @param edge The number of pixels to be added in the edge of the segmented
#'   object. Defaults to 5.
#' @param remove_bg If `TRUE`, the pixels that are not part of objects are
#'   converted to white.
#' @param ... Additional arguments passed on to [image_combine()]
#' @return A list of objects of class `Image`.
#' @export
#' @seealso [analyze_objects()], [image_binary()]
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg", plot = TRUE)
#' imgs <- object_split(img) # set to NULL to use 50% of the cores
#' }
#'
object_split <- function(img,
                         index = "NB",
                         lower_size = NULL,
                         watershed = TRUE,
                         invert = FALSE,
                         fill_hull = FALSE,
                         opening = 3,
                         closing = FALSE,
                         filter = FALSE,
                         erode = FALSE,
                         dilate = FALSE,
                         threshold = "Otsu",
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         edge = 3,
                         remove_bg = FALSE,
                         plot = TRUE,
                         verbose = TRUE,
                         ...){
  check_ebi()

  img2 <- help_binary(img,
                      opening = opening,
                      closing = closing,
                      filter = filter,
                      erode = erode,
                      dilate = dilate,
                      index = index,
                      invert = invert,
                      fill_hull = fill_hull,
                      threshold = threshold)
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

  objcts <- get_area_mask(nmask)
  av_area <- mean(objcts)
  ifelse(!is.null(lower_size),
         cutsize <- lower_size,
         cutsize <-  av_area * 0.1)
  selected <- which(objcts > cutsize)

  split_objects <- function(img, nmask){
    objects <- help_isolate_object(img[,,1], img[,,2], img[,,3], nmask, remove_bg, edge)
    lapply(seq_along(objects), function(x){
      dimx <- dim(objects[[x]][[1]])
      EBImage::Image(array(c(objects[[x]][[1]], objects[[x]][[2]], objects[[x]][[3]]), dim = c(dimx, 3)), colormode = "Color")
    })
  }
  list_objects <- split_objects(img, nmask)
  names(list_objects) <- 1:length(list_objects)
  list_objects <- list_objects[selected]
  if(isTRUE(verbose)){
    cat("==============================\n")
    cat("Summary of the procedure\n")
    cat("==============================\n")
    cat("Number of objects:", length(objcts), "\n")
    cat("Average area     :", mean(objcts), "\n")
    cat("Minimum area     :", min(objcts), "\n")
    cat("Maximum area     :", max(objcts), "\n")
    cat("Objects created  :", length(list_objects), "\n")
    cat("==============================\n")
  }
  if(isTRUE(plot)){
    image_combine(list_objects, ...)
  }
  invisible(list_objects)
}


#' Augment Images
#'
#' This function takes an image and augments it by rotating it multiple times.
#' @inheritParams analyze_objects
#' @param img An `Image` object.
#' @param pattern A regular expression pattern to select multiple images from a
#'   directory.
#' @param times The number of times to rotate the image.
#' @param type The type of output: "export" to save images or "return" to return
#'   a list of augmented images.
#' @param dir_original The directory where original images are located.
#' @param dir_processed The directory where processed images will be saved.
#' @param parallel Whether to perform image augmentation in parallel.
#' @param verbose Whether to display progress messages.
#'
#' @return If type is "export," augmented images are saved. If type is "return,"
#'   a list of augmented images is returned.
#'
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' imgs <- image_augment(img, type = "return", times = 4)
#' image_combine(imgs)
#' }
#'
image_augment <- function(img,
                          pattern = NULL,
                          times = 12,
                          type = "export",
                          dir_original = NULL,
                          dir_processed = NULL,
                          parallel = FALSE,
                          verbose = TRUE,
                          workers = NULL){
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


  if(is.null(pattern)){
    angles <- seq(0, 360, by = 360 / times)
    angles <- angles[-length(angles)]
    obj_list <- list()
    for(i in 1:times){
      top <- img@.Data[1:10,,]
      bottom <- img@.Data[(nrow(img)-10):nrow(img),,]
      left <- img@.Data[,1:10,]
      right <- img@.Data[,(ncol(img) - 10):ncol(img),]

      rval <- mean(c(c(top[,,1]), c(bottom[,,1]), c(left[,,1]), c(right[,,1])))
      gval <- mean(c(c(top[,,2]), c(bottom[,,2]), c(left[,,2]), c(right[,,2])))
      bval <- mean(c(c(top[,,3]), c(bottom[,,3]), c(left[,,3]), c(right[,,3])))

      tmp <- EBImage::rotate(img, angles[i], bg.col = rgb(rval, gval, bval))
      if(type == "export"){
        image_export(tmp,
                     name = paste0("v", sub("\\.", "_", round(angles[i], 2)), ".jpg"),
                     subfolder = diretorio_processada)
      } else{
        obj_list[[paste0("v_", sub("\\.", "_", round(angles[i], 2)), ".jpg")]] <- tmp
      }
    }
  } else{

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

    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if (length(grep(pattern, names_plant)) == 0) {
      cli::cli_abort(c(
        "!" = "Pattern {.val {pattern}} not found.",
        "x" = "Not found in directory {.path {file.path(getwd(), sub('^\\.', '', diretorio_original))}}."
      ))
    }

    if (!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))) {
      cli::cli_abort("Allowed extensions are: {.val .png}, {.val .jpeg}, {.val .jpg}, {.val .tiff}")
    }


    if(isTRUE(parallel)){

      # 1. make sure the output directory exists
      if (!dir.exists(diretorio_processada)) {
        dir.create(diretorio_processada, recursive = TRUE)
      }

      # decide number of workers
      nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.3) else workers

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)


      # CLI header + progress
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Augmenting {.val {length(plants)}} images"),
          right = cli::col_blue("Using {.val {nworkers}} workers")
        )
        cli::cli_progress_step(
          msg        = "Dispatching rotation tasks...",
          msg_done   = "All images processed.",
          msg_failed = "Something went wrong."
        )
      }

      # run in parallel
      mirai::mirai_map(
        .x = plants,
        .f = function(path) {
          tmpimg <- pliman::image_import(path, path = diretorio_original)
          angles <- seq(0, 360, length.out = times + 1)[- (times + 1)]

          for (ang in angles) {
            # compute background colour
            frame  <- tmpimg@.Data
            border <- c(
              c(frame[1:10,,1]),   c(frame[(nrow(tmpimg)-9):nrow(tmpimg),,1]),
              c(frame[,1:10,1]),   c(frame[,(ncol(tmpimg)-9):ncol(tmpimg),1])
            )
            rval <- mean(border)
            # same for G/B...
            gval <- mean(c(frame[1:10,,2], frame[(nrow(tmpimg)-9):nrow(tmpimg),,2],
                           frame[,1:10,2], frame[,(ncol(tmpimg)-9):ncol(tmpimg),2]))
            bval <- mean(c(frame[1:10,,3], frame[(nrow(tmpimg)-9):nrow(tmpimg),,3],
                           frame[,1:10,3], frame[,(ncol(tmpimg)-9):ncol(tmpimg),3]))

            rotated <- EBImage::rotate(tmpimg, ang, bg.col = grDevices::rgb(rval, gval, bval))

            if (type == "export") {
              # build a full file path
              fname <- file.path(
                diretorio_processada,
                paste0(pliman::file_name(path), "_", sub("\\.", "-", round(ang, 2)), ".jpg")
              )
              # write it out
              pliman::image_export(rotated, name = fname, subfolder = NULL)
            } else {
              # return objects if needed
              NULL
            }
          }

          NULL
        }
      )[.progress]


    } else{
      obj_list <- list()
      for(i in seq_along(plants)){

        tmpimg <- image_import(plants[[i]], path = diretorio_original)
        angles <- seq(0, 360, by = 360 / times)
        angles <- angles[-length(angles)]
        for(j in 1:times){
          top <- tmpimg@.Data[1:10,,]
          bottom <- tmpimg@.Data[(nrow(tmpimg)-10):nrow(tmpimg),,]
          left <- tmpimg@.Data[,1:10,]
          right <- tmpimg@.Data[,(ncol(tmpimg) - 10):ncol(tmpimg),]

          rval <- mean(c(c(top[,,1]), c(bottom[,,1]), c(left[,,1]), c(right[,,1])))
          gval <- mean(c(c(top[,,2]), c(bottom[,,2]), c(left[,,2]), c(right[,,2])))
          bval <- mean(c(c(top[,,3]), c(bottom[,,3]), c(left[,,3]), c(right[,,3])))

          tmp <- EBImage::rotate(tmpimg, angles[j], bg.col = rgb(rval, gval, bval))


          if(type == "export"){
            image_export(tmp,
                         name = paste0(file_name(plants[[i]]), "_", sub("\\.", "-", round(angles[j], 2)), ".jpg"),
                         subfolder = diretorio_processada)
          } else{
            obj_list[[paste0(file_name(plants[[i]]), "_", sub("\\.", "-", round(angles[j], 2)), ".jpg")]] <- tmp
          }




        }
      }
    }
  }

  if(type == "return"){
    invisible(obj_list)
  }

}


#' Export multiple objects from an image to multiple images
#'
#' Givin an image with multiple objects, `object_export()` will split the
#' objects into a list of objects using [object_split()] and then export them to
#' multiple images into the current working directory (or a subfolder). Batch
#' processing is performed by declaring a file name pattern that matches the
#' images within the working directory.
#'
#' @inheritParams object_split
#' @inheritParams utils_image
#' @inheritParams analyze_objects
#' @inheritParams image_augment
#'
#' @param pattern A pattern of file name used to identify images to be
#'   processed. For example, if `pattern = "im"` all images in the current
#'   working directory that the name matches the pattern (e.g., img1.-,
#'   image1.-, im2.-) will be imported and processed. Providing any number as
#'   pattern (e.g., `pattern = "1"`) will select images that are named as 1.-,
#'   2.-, and so on. An error will be returned if the pattern matches any file
#'   that is not supported (e.g., img1.pdf).
#' @param augment A logical indicating if exported objects should be augmented using
#'   [image_augment()]. Defaults to `FALSE`.
#'@param dir_original The directory containing the original images. Defaults to
#'  `NULL`. It can be either a full path, e.g., `"C:/Desktop/imgs"`, or a
#'  subfolder within the current working directory, e.g., `"/imgs"`.
#' @param dir_processed Optional character string indicating a subfolder within the
#'   current working directory to save the image(s). If the folder doesn't
#'   exist, it will be created.
#' @param format The format of image to be exported.
#' @param squarize Squarizes the image before the exportation? If `TRUE`,
#'   [image_square()] will be called internally.
#' @return A `NULL` object.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("potato_leaves.jpg")
#' object_export(img,
#'               remove_bg = TRUE)
#' }
object_export <- function(img,
                          pattern = NULL,
                          dir_original = NULL,
                          dir_processed = NULL,
                          format = ".jpg",
                          squarize = FALSE,
                          augment = FALSE,
                          times = 12,
                          index = "NB",
                          lower_size = NULL,
                          watershed = FALSE,
                          invert = FALSE,
                          fill_hull = FALSE,
                          opening = 3,
                          closing = FALSE,
                          filter = FALSE,
                          erode = FALSE,
                          dilate = FALSE,
                          threshold = "Otsu",
                          extension = NULL,
                          tolerance = NULL,
                          object_size = "medium",
                          edge = 20,
                          remove_bg = FALSE,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  if(is.null(pattern)){
    list_objects <- object_split(img = img,
                                 index = index,
                                 lower_size = lower_size,
                                 watershed = watershed,
                                 invert = invert,
                                 fill_hull = fill_hull,
                                 opening = opening,
                                 closing = closing,
                                 erode = erode,
                                 dilate = dilate,
                                 filter = filter,
                                 threshold = threshold,
                                 extension = extension,
                                 tolerance = tolerance,
                                 object_size = object_size,
                                 edge = edge,
                                 remove_bg = remove_bg,
                                 plot = FALSE,
                                 verbose = FALSE)
    names(list_objects) <-  leading_zeros(as.numeric(names(list_objects)), n = 4)

    if(isTRUE(augment)){
      bb <-
        lapply(seq_along(list_objects), function(x){
          image_augment(list_objects[[x]], type = "return", times = times)
        })
      names(bb) <- names(list_objects)
      unlisted <- do.call(c, bb)
      names(unlisted) <- sub("\\.", "_", names(unlisted))
      list_objects <- unlisted
    }


    a <- lapply(seq_along(list_objects), function(i){
      tmp <- list_objects[[i]]
      if(isTRUE(squarize)){
        tmp <- image_square(tmp,
                            plot = FALSE,
                            sample_left = 5,
                            sample_top = 5,
                            sample_right = 5,
                            sample_bottom = 5)
      }
      image_export(tmp,
                   name = paste0(file_name(names(list_objects[i])), ".jpg"),
                   subfolder = dir_processed)
    })
  } else{

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

    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if (length(grep(pattern, names_plant)) == 0) {
      cli::cli_abort(c(
        "!" = "Pattern {.val {pattern}} not found.",
        "x" = "Not found in directory {.path {file.path(getwd(), sub('^\\.', '', diretorio_original))}}."
      ))
    }

    if (!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))) {
      cli::cli_abort("Allowed extensions are: {.val .png}, {.val .jpeg}, {.val .jpg}, {.val .tiff}")
    }


    if(isTRUE(parallel)){

      # ensure output dir exists up front
      if (!dir.exists(diretorio_processada)) {
        dir.create(diretorio_processada, recursive = TRUE)
      }

      # decide number of workers
      nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.3) else workers

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # CLI header + progress step
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Augmenting and exporting {.val {length(plants)}} objects"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Dispatching {.val {length(plants)}} image tasks...",
          msg_done   = "All batches finished.",
          msg_failed = "Batch processing failed."
        )
      }

      # run in parallel
      mirai::mirai_map(
        .x = plants,
        .f = function(path) {
          tmpimg <- pliman::image_import(path, path = diretorio_original)

          # split into objects
          list_objects <- object_split(
            img         = tmpimg,
            index       = index,
            lower_size  = lower_size,
            watershed   = watershed,
            invert      = invert,
            fill_hull   = fill_hull,
            opening     = opening,
            closing     = closing,
            filter      = filter,
            threshold   = threshold,
            extension   = extension,
            tolerance   = tolerance,
            object_size = object_size,
            edge        = edge,
            remove_bg   = remove_bg,
            verbose     = FALSE,
            plot        = FALSE
          )
          names(list_objects) <- paste0(
            leading_zeros(as.numeric(names(list_objects)), n = 4),
            ".jpg"
          )

          # optional augmentation
          if (isTRUE(augment)) {
            list_objects <- unlist(lapply(names(list_objects), function(nm) {
              imgs <- image_augment(list_objects[[nm]], type = "return", times = times)
              setNames(imgs, sub("\\.", "_", names(imgs)))
            }), recursive = FALSE)
            names(list_objects) <- sub("\\.jpg$", "", names(list_objects))
          }
          return(list_objects)
          # export each object to an absolute path
          for (nm in names(list_objects)) {
            obj <- list_objects[[nm]]
            if (isTRUE(squarize)) {
              obj <- try(
                image_square(obj, plot = FALSE,
                             sample_left   = 5,
                             sample_top    = 5,
                             sample_right  = 5,
                             sample_bottom = 5),
                silent = TRUE
              ) %||% obj
            }
            out_file <- file.path(
              diretorio_processada,
              paste0(file_name(path), "_", nm, format)
            )
            pliman::image_export(
              obj,
              name      = out_file,
              subfolder = NULL
            )
          }

          NULL
        }
      )[.progress]

    } else{
      cli::cli_rule(
        left  = cli::col_blue("Augmenting and exporting {.val {length(plants)}} objects"),
        right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        total = length(names_plant),
        clear = FALSE
      )
      for(i in seq_along(plants)){
        tmpimg <- image_import(plants[[i]], path = diretorio_original)
        cli::cli_progress_update()

        list_objects <- object_split(img = tmpimg,
                                     index = index,
                                     lower_size = lower_size,
                                     watershed = watershed,
                                     invert = invert,
                                     fill_hull = fill_hull,
                                     opening = opening,
                                     closing = closing,
                                     filter = filter,
                                     threshold = threshold,
                                     extension = extension,
                                     tolerance = tolerance,
                                     object_size = object_size,
                                     edge = edge,
                                     remove_bg = remove_bg,
                                     verbose = FALSE,
                                     plot = FALSE)
        names(list_objects) <-  paste0(leading_zeros(as.numeric(names(list_objects)), n = 4), ".jpg")
        if(isTRUE(augment)){
          bb <-
            lapply(seq_along(list_objects), function(x){
              image_augment(list_objects[[x]], type = "return", times = times)
            })
          names(bb) <- names(list_objects)
          unlisted <- do.call(c, bb)
          names(unlisted) <- sub("\\.", "_", names(unlisted))
          list_objects <- unlisted
          names(list_objects) <- sub("jpg.", "", names(list_objects))
        }

        a <- lapply(seq_along(list_objects), function(j){
          tmp <- list_objects[[j]]
          if(isTRUE(squarize)){
            try(
              tmp <- image_square(tmp,
                                  plot = FALSE,
                                  sample_left = 5,
                                  sample_top = 5,
                                  sample_right = 5,
                                  sample_bottom = 5),
              silent = TRUE
            )

          }
          image_export(tmp,
                       name = paste0(file_name(plants[[i]]), "_", names(list_objects[j])),
                       subfolder = diretorio_processada)
        }
        )
      }
    }

  }
}



#' Extract red, green and blue values from objects
#'
#' Given an image and a matrix of labels that identify each object, the function
#' extracts the red, green, and blue values from each object.
#'
#' @param img An `Image` object
#' @param labels A mask containing the labels for each object. This can be
#'   obtained with [EBImage::bwlabel()] or [EBImage::watershed()]
#'
#' @return A data.frame with `n` rows (number of pixels for all the objects) and
#'   the following columns:
#'  * `id`: the object id;
#'  * `R`: the value for the red band;
#'  * `G`: the value for the blue band;
#'  * `B`: the value for the green band;
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' # segment the objects using the "B" (blue) band (default)
#'
#' labs <- object_label(img, watershed = TRUE)
#' rgb <- object_rgb(img, labs[[1]])
#' head(rgb)
#' }
object_rgb <- function(img, labels){
  dd <- help_get_rgb(img[,,1], img[,,2], img[,,3], labels)
  df2 <- data.frame(do.call(rbind,  lapply(dd, function(x){
    matrix(x, ncol = 4, byrow = TRUE)
  })))
  colnames(df2) <- c("id", "R", "G", "B")
  if(dim(img)[[3]] == 5){
    renir <- help_get_renir(img[,,4], img[,,5], labels)
    df3 <- data.frame(do.call(rbind,  lapply(renir, function(x){
      matrix(x, ncol = 3, byrow = TRUE)
    })))
    df2 <- cbind(df2, df3[, 2:3])
    colnames(df2) <- c("id", "R", "G", "B", "RE", "NIR")
  }
  invisible(df2)
}



#' Apply color to image objects
#'
#' The function applies the color informed in the argument `color` to segmented
#' objects in the image. The segmentation is performed using image indexes. Use
#' [image_index()] to identify the better candidate index to segment objects.
#'
#' @inheritParams image_binary
#' @param color The color to apply in the image objects. Defaults to `"blue"`.
#' @param plot Plots the modified image? Defaults to `TRUE`.
#' @param pick_palettes  Logical argument indicating wheater the user needs to
#'   pick up the color palettes for foreground and background for the image. If
#'   `TRUE` [pick_palette()] will be called internally so that the user can sample
#'   color points representing foreground and background.
#'@param foreground,background A color palette for the foregrond and background,
#'  respectively (optional).
#' @param ... Additional arguments passed on to [image_binary()].
#'
#' @return An object of class `Image`
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg")
#' img2 <- object_to_color(img, index = "G-R")
#' image_combine(img, img2)
#' }
#'
object_to_color <- function(img,
                            pick_palettes = FALSE,
                            background = NULL,
                            foreground = NULL,
                            index = "NB",
                            color = "blue",
                            plot = TRUE,
                            ...){
  if(isTRUE(pick_palettes) || (!is.null(background) & !is.null(foreground))){
    if(interactive()){
      plot(img)
      if(is.null(background) && is.null(foreground)){
        cli::cli_inform("Use the first mouse button to pick up BACKGROUND colors. Press {.kbd Esc} to exit.")

        back <- pick_palette(img,
                             r = 5,
                             verbose = FALSE,
                             palette  = FALSE,
                             plot = FALSE,
                             col = "blue",
                             external_device = FALSE,
                             title = "Use the first mouse button to pick up BACKGROUND colors. Click 'Done' to finish",
                             viewer = "base")
      } else{
        back <- background
      }
      if(is.null(background) && is.null(foreground)){
        cli::cli_inform("Use the first mouse button to pick up FOREGROUND colors. Press {.kbd Esc} to exit.")

        fore <- pick_palette(img,
                             r = 5,
                             verbose = FALSE,
                             palette  = FALSE,
                             plot = FALSE,
                             col = "salmon",
                             external_device = FALSE,
                             title = "Use the first mouse button to pick up FOREGROUND colors. Click 'Done' to finish",
                             viewer = "base")
      } else{
        fore <- foreground
      }

      original <-
        data.frame(CODE = "img",
                   R = c(img@.Data[,,1]),
                   G = c(img@.Data[,,2]),
                   B = c(img@.Data[,,3]))
      foreground <-
        data.frame(CODE = "foreground",
                   R = c(fore@.Data[,,1]),
                   G = c(fore@.Data[,,2]),
                   B = c(fore@.Data[,,3]))
      background <-
        data.frame(CODE = "background",
                   R = c(back@.Data[,,1]),
                   G = c(back@.Data[,,2]),
                   B = c(back@.Data[,,3]))
      back_fore <-
        transform(rbind(foreground[sample(1:nrow(foreground)),][1:2000,],
                        background[sample(1:nrow(background)),][1:2000,]),
                  Y = ifelse(CODE == "background", 0, 1))

      formula <- as.formula(paste("Y ~ ", "R+G+B"))

      modelo1 <- suppressWarnings(glm(formula,
                                      family = binomial("logit"),
                                      data = back_fore))
      pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
      bin <- EBImage::Image(matrix(pred1, ncol = dim(img)[[2]]))

    }
  } else{
    bin <- help_binary(img,
                       index = index,
                       ...)
  }

  pix_ref <- which(bin == 1)
  colto <- col2rgb(color) / 255
  img@.Data[,,1][pix_ref] <- colto[1]
  img@.Data[,,2][pix_ref] <- colto[2]
  img@.Data[,,3][pix_ref] <- colto[3]
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}

#' Compute Bounding Boxes from Contours
#'
#' This function calculates the bounding boxes for a given list of contours.
#'
#' @param contours A list of matrices, where each matrix contains two columns
#'   representing (x, y) coordinates of a contour.
#'
#' @return A list of bounding boxes, where each bounding box is represented as a
#'   list with `x_min`, `y_min`, `x_max`, and `y_max` values.
#'
#' @examples
#' if(interactive()){
#' contours <- list(
#'   matrix(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
#'            110, 120, 130, 140, 150, 160, 170, 180, 190, 200),
#'          ncol = 2, byrow = FALSE)
#' )
#' bbox_list <- object_bbox(contours)
#' print(bbox_list)
#' }
#'
#' @export
object_bbox <- function(contours) {
  # Ensure contours is a list
  if (!is.list(contours)) {
    cli::cli_abort("contours must be a list of coordinate matrices")
  }
  bbox_list <- lapply(contours, function(coords) {
    list(
      x_min = min(coords[, 1]),
      y_min = min(coords[, 2]),
      x_max = max(coords[, 1]),
      y_max = max(coords[, 2])
    )
  })
  return(bbox_list)
}

#' Add Bounding Boxes to an Existing Plot
#'
#' This function overlays bounding boxes onto an existing plot.
#'
#' @param bbox_list A list of bounding boxes, as returned by `object_bbox()`.
#' @param col The color for the bounding boxes. Defaults to `"red"`.
#' @return None (adds bounding boxes to an existing plot).
#'
#' @examples
#' if(interactive()){
#' plot(NA,
#'     xlim = c(0, 200),
#'     ylim = c(0, 200),
#'     asp = 1)
#' contours <- list(
#'   matrix(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
#'            110, 120, 130, 140, 150, 160, 170, 180, 190, 200),
#'          ncol = 2, byrow = FALSE)
#' )
#' bbox_list <- object_bbox(contours)
#' plot_bbox(bbox_list)
#' }
#'
#' @export
plot_bbox <- function(bbox_list, col = "red") {
  if(is.matrix(bbox_list[[1]])){
    bbox_list <- object_bbox(bbox_list)
  }
  if (!is.list(bbox_list) || length(bbox_list) == 0) {
    cli::cli_abort("bbox_list must be a non-empty list of bounding boxes.")
  }
  for (bbox in bbox_list) {
    rect(bbox$x_min, bbox$y_min, bbox$x_max, bbox$y_max, border = col, lwd = 1)
  }
}


#' @title Plot object thumbnails at (x, y) coordinates derived from image features
#'
#' @description
#' Extracts connected objects from an image, computes their features, crops each
#' object, converts the crop to a raster with alpha, and draws each thumbnail
#' centered at its corresponding `(x, y)` feature location in a ggplot.
#' Optionally overlays object IDs. Caching can be used to avoid recomputing
#' object features on repeated calls.
#' @param img An image of class `EBImage::Image` (or compatible) from which
#'   objects will be segmented and measured.
#' @param x Character scalar. Name of the feature column (returned by
#'   `get_measures(res)`) to use on the x-axis.
#' @param y Character scalar. Name of the feature column to use on the y-axis.
#' @param scale Numeric in (0, 1]. Relative thumbnail size as a fraction of the
#'   data range along each axis (larger values draw larger thumbnails).
#' @param xy_ratio Numeric scalar. Factor applied to the vertical scaling
#'   (y-axis) of thumbnails. Use values other than 1 to stretch or compress
#'   thumbnails vertically.
#' @param xlab,ylab Character scalars used as x- and y-axis labels. Defaults to
#'   `x` and `y`, respectively.
#' @param erosion,dilatation Integer (non-negative). Size of the structuring
#'   element for morphological erosion/dilatation of the segmented objects.
#' @param show_id Logical. If `TRUE`, overlays object IDs at their
#'   `x, y` locations.
#' @param color_id Character. Color used for the ID labels when
#'   `show_id = TRUE`.
#' @param size_id Numeric. Text size for the ID labels when
#'   `show_id = TRUE`.
#' @param cache Logical. If `TRUE` (default), caches results of object
#'   extraction using a simple key based on image dimensions and parameters.
#' @param verbose If `TRUE` (default), shows the progress of analysis.
#' @param ... Additional arguments forwarded to [analyze_objects()].
#' @return A list with two elements:
#' * `features` — a data frame (or tibble) with object-level features
#'   returned by `get_measures(res)`. Must contain columns named `x` and `y`.
#' * `plot` — a `ggplot` object. The thumbnail scatter plot.
#'
#' @section Scaling behavior:
#' Thumbnails are sized relative to the observed ranges in `x` and `y`.
#' If the two axes differ substantially in range, perceived thumbnail aspect
#' on the plotting device may vary. Use `xy_ratio` to adjust vertical scaling.
#' @importFrom grDevices as.raster rgb extendrange
#'
#' @export
#' @examples
#' if(interactive()){
#' img <- image_pliman("potato_leaves.jpg")
#' plot(img)
#' res <- object_scatter(
#'  img = img,
#'  index = "B",
#'  x = "area",
#'  y = "solidity",
#'  watershed = FALSE,
#'  scale = 0.5
#' )
#' res$plot
#'
#' # remove cached data
#' clear_pliman_cache()
#' }
#'

object_scatter <- function(img,
                           x,
                           y,
                           scale = 0.1,
                           xy_ratio = 1,
                           xlab = x,
                           ylab = y,
                           erosion = 2,
                           dilatation = FALSE,
                           show_id = FALSE,
                           color_id = "black",
                           size_id = 3,
                           cache = TRUE,
                           verbose = TRUE,
                           ...) {
  cache_key <- NULL
  cache_dir <- tools::R_user_dir("pliman", which = "cache")

  # ---- helpers -------------------------------------------------------------
  ebimg_to_raster <- function(img, flip_vertical = FALSE) {
    a <- EBImage::imageData(img)
    w <- dim(a)[1]; h <- dim(a)[2]
    c <- ifelse(length(dim(a)) == 3, dim(a)[3], 1)

    norm01 <- function(x) {
      r <- range(x, finite = TRUE)
      if (!is.finite(r[1]) || r[1] == r[2]) return(ifelse(is.finite(x), 0, x))
      (x - r[1]) / (r[2] - r[1])
    }

    if (c == 1) {
      g <- as.vector(norm01(a)); col <- grDevices::rgb(g, g, g, 1)
    } else {
      r <- as.vector(norm01(a[, , 1]))
      g <- as.vector(norm01(a[, , 2]))
      b <- as.vector(norm01(a[, , 3]))
      al <- if (c >= 4) as.vector(norm01(a[, , 4])) else 1
      col <- grDevices::rgb(r, g, b, alpha = al)
    }
    mat <- t(matrix(col, nrow = w, ncol = h))
    if (flip_vertical) mat <- mat[h:1, , drop = FALSE]
    grDevices::as.raster(mat)
  }

  make_key <- function(obj) {
    rw  <- serialize(obj, NULL, xdr = FALSE)
    v   <- as.integer(rw)
    mod <- .Machine$integer.max
    h <- 0
    for (i in seq_along(v)) {
      h <- (h * 131 + v[i] + i) %% mod
    }
    sprintf("%08x", as.integer(h))
  }

  # ---- caching -------------------------------------------------------------
  if (isTRUE(cache)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    key_input <- list(
      img_dim    = dim(img),
      erosion    = erosion,
      dilatation = dilatation,
      dots       = as.list(match.call(expand.dots = FALSE)$...)
    )
    cache_key  <- make_key(key_input)
  }
  cache_path <- if (isTRUE(cache)) file.path(cache_dir, paste0("object_scatter_", cache_key, ".rds")) else NULL

  if (isTRUE(cache) && length(cache_path) && file.exists(cache_path)) {
    if(verbose){
      cli::cli_progress_step("Getting cached data...")
    }
    results <- readRDS(cache_path)
    vals    <- results$vals
    rasters <- results$rasters
  } else {
    if(verbose){
      cli::cli_progress_step("Extracting object features...")
    }
    res <- analyze_objects(img,
                           ...,
                           erode = FALSE,
                           dilate = FALSE,
                           plot = FALSE,
                           show_contour = FALSE,
                           return_mask = TRUE,
                           verbose = FALSE)

    bb   <- object_bbox(res[["contours"]])
    mask <- res$mask
    if (is.numeric(erosion)) {
      mask <- image_erode(mask, size = erosion)
    }
    if (is.numeric(dilatation)) {
      mask <- image_dilate(mask, size = dilatation)
    }
    ima <- image_alpha(img, mask)

    # base R instead of purrr::map
    rasters <- lapply(bb, function(b) {
      ebimg_to_raster(ima[b$x_min:b$x_max, b$y_min:b$y_max, ])
    })

    vals <- get_measures(res)

    if (isTRUE(cache) && length(cache_path)) {
      if(verbose){
        cli::cli_progress_step(
          "Saving data in cache. You can clear cached data with {.fn clear_pliman_cache}"
        )
      }
      saveRDS(list(vals = vals, rasters = rasters), cache_path)
    }
  }

  if(verbose){
    cli::cli_progress_step("Putting objects in their positions...")
  }

  if (!(x %in% colnames(vals))) {
    cli::cli_abort(c(
      "!" = "Column {.val {x}} not found in computed features.",
      "i" = "Available columns are: {.val {colnames(vals)}}"
    ))
  }
  if (!(y %in% colnames(vals))) {
    cli::cli_abort(c(
      "!" = "Column {.val {y}} not found in computed features.",
      "i" = "Available columns are: {.val {colnames(vals)}}"
    ))
  }

  xval <- vals[[x]]
  yval <- vals[[y]]

  if (!(length(xval) == length(yval) && length(xval) == length(rasters))) {
    cli::cli_abort(c(
      "!" = "Inconsistent lengths detected.",
      "x" = "`xval` has {length(xval)} values.",
      "x" = "`yval` has {length(yval)} values.",
      "x" = "`rasters` has {length(rasters)} elements."
    ))
  }

  xlim0 <- range(xval)
  ylim0 <- range(yval)
  wpx <- vapply(rasters, ncol, integer(1))
  hpx <- vapply(rasters, nrow, integer(1))
  scale_x <- wpx / max(wpx) * scale
  scale_y <- hpx / max(hpx) * scale * xy_ratio
  xminp <- xmaxp <- yminp <- ymaxp <- numeric(length(xval))
  for (i in seq_along(xval)) {
    wi <- diff(xlim0) / 2 * scale_x[i]
    hi <- diff(ylim0) / 2 * scale_y[i]
    xminp[i] <- xval[i] - wi; xmaxp[i] <- xval[i] + wi
    yminp[i] <- yval[i] - hi; ymaxp[i] <- yval[i] + hi
  }
  xlim <- extendrange(range(xminp, xmaxp))
  ylim <- extendrange(range(yminp, ymaxp))
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(bg = "white", bty = "n")

  xt <- pretty(xlim, n = 5)
  yt <- pretty(ylim, n = 5)
  plot(NA, xlim = range(xt), ylim = range(yt), xlab = xlab, ylab = ylab, type = "n")
  abline(v = xt, col = "grey85", lty = "dotted")
  abline(h = yt, col = "grey85", lty = "dotted")
  axis(1, at = xt, labels = format(xt, trim = TRUE), lwd = 0, lwd.ticks = 1)
  axis(2, at = yt, labels = format(yt, trim = TRUE), lwd = 0, lwd.ticks = 1)
  # draw rasters
  for (i in seq_along(xval)) {
    graphics::rasterImage(
      rasters[[i]],
      xleft   = xminp[i],
      ybottom = yminp[i],
      xright  = xmaxp[i],
      ytop    = ymaxp[i],
      interpolate = TRUE
    )
  }

  if (isTRUE(show_id)) {
    graphics::text(xval, yval, labels = vals$id, col = color_id, cex = size_id / 3)
  }

  # capture a replayable plot object
  rp <- grDevices::recordPlot()
  invisible(list(features = vals, plot = rp))
}

