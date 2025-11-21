#'Combines images to a grid
#'
#'Combines several images to a grid
#' @param ... a comma-separated name of image objects or a list containing image
#'   objects.
#' @param labels A character vector with the same length of the number of
#'   objects in `...` to indicate the plot labels.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param col The color for the plot labels. Defaults to `col = "black"`.
#' @param verbose Shows the name of objects declared in `...` or a numeric
#'   sequence if a list with no names is provided. Set to `FALSE` to supress the
#'   text.
#' @importFrom stats reshape IQR quantile
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A grid with the images in `...`
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img1 <- image_pliman("sev_leaf.jpg")
#' img2 <- image_pliman("sev_leaf_nb.jpg")
#' image_combine(img1, img2)
#' }
image_combine <- function(...,
                          labels = NULL,
                          nrow = NULL,
                          ncol = NULL,
                          col = "black",
                          verbose = TRUE){
  if(is.list(c(...))){
    plots <- as.list(...)
    if(class(plots) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      plots <- lapply(plots, function(x){x[[1]]})
    }
    if(!is.null(labels)){
      names(plots) <- labels
    }
  }else{
    plots <- list(...)
    if(is.null(labels)){
      names(plots) <- unlist(strsplit(gsub("c\\(|\\)",  "", substitute(c(...))), "\\s*(\\s|,)\\s*"))[-1]
    } else{
      names(plots) <- labels
    }
  }
  num_plots <- length(plots)
  if (is.null(nrow) && is.null(ncol)){
    ncol <- ceiling(sqrt(num_plots))
    nrow <- ceiling(num_plots/ncol)
  }
  if (is.null(ncol)){
    ncol <- ceiling(num_plots/nrow)
  }
  if (is.null(nrow)){
    nrow <- ceiling(num_plots/ncol)
  }
  op <- par(mfrow = c(nrow, ncol))
  on.exit(par(op))
  ifelse(is.null(names(plots)), index <- 1:length(plots), index <- names(plots))
  for(i in 1:length(plots)){
    plot(plots[[i]])
    if(verbose == TRUE){
      dim <- image_dimension(plots[[i]], verbose = FALSE)
      text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = col)
    }
  }
}

#'Import and export images
#'
#'Import images from files and URLs and write images to files, possibly with
#'batch processing.
#' @name utils_image
#' @param img
#' * For `image_import()`, a character vector of file names or URLs.
#' * For `image_input()`, a character vector of file names or URLs or an array
#' containing the pixel intensities of an image.
#' * For `image_export()`, an Image object, an array or a list of images.
#' * For `image_pliman()`, a charactere value specifying the image example. See
#' `?pliman_images` for more details.
#' @param which logical scalar or integer vector to indicate which image are
#'   imported if a TIFF files is informed. Defaults to `1` (the first image is
#'   returned).
#' @param name An string specifying the name of the image. It can be either a
#'   character with the image name (e.g., "img1") or name and extension (e.g.,
#'   "img1.jpg"). If none file extension is provided, the image will be saved as
#'   a *.jpg file.
#' @param prefix A prefix to include in the image name when exporting a list of
#'   images. Defaults to `""`, i.e., no prefix.
#' @param extension When `image` is a list, `extension` can be used to define
#'   the extension of exported files. This will overwrite the file extensions
#'   given in `image`.
#' @param pattern A pattern of file name used to identify images to be imported.
#'   For example, if `pattern = "im"` all images in the current working
#'   directory that the name matches the pattern (e.g., img1.-, image1.-, im2.-)
#'   will be imported as a list. Providing any number as pattern (e.g., `pattern
#'   = "1"`) will select images that are named as 1.-, 2.-, and so on. An error
#'   will be returned if the pattern matches any file that is not supported
#'   (e.g., img1.pdf).
#' @param subfolder Optional character string indicating a subfolder within the
#'   current working directory to save the image(s). If the folder doesn't
#'   exist, it will be created.
#' @param path A character vector of full path names; the default corresponds to
#'   the working directory, [getwd()]. It will overwrite (if given) the path
#'   informed in `image` argument.
#' @param resize Resize the image after importation? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param plot Plots the image after importing? Defaults to `FALSE`.
#' @param nrow,ncol Passed on to [image_combine()]. The number of rows and
#'   columns to use in the composite image when `plot = TRUE`.
#' @param ...
#'  * For `image_import()` alternative arguments passed to the corresponding
#'  functions from the `jpeg`, `png`, and `tiff` packages.
#'  * For `image_input()` further arguments passed on to [EBImage::Image()].
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_import()` returns a new `Image` object.
#' * `image_export()` returns an invisible vector of file names.
#' * `image_pliman()` returns a new `Image` object with the example image
#' required. If an empty call is used, the path to the `tmp_images` directory
#' installed with the package is returned.
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' folder <- image_pliman()
#' full_path <- paste0(folder, "/sev_leaf.jpg")
#' (path <- file_dir(full_path))
#' (file <- basename(full_path))
#' image_import(img = full_path)
#' image_import(img = file, path = path)
#' }
image_import <- function(img,
                         ...,
                         which = 1,
                         pattern = NULL,
                         path = NULL,
                         resize = FALSE,
                         plot = FALSE,
                         nrow = NULL,
                         ncol = NULL){
  check_ebi()
  valid_extens <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF", "TIF", "tif", "gri", "grd")
  if(!is.null(pattern)){
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    path <- ifelse(is.null(path), getwd(), path)
    imgs <- list.files(pattern = pattern, path)
    if(length(grep(pattern, imgs)) == 0){
      cli::cli_abort("Pattern {.val {pattern}} not found in {.dir {path}}.")
    }
    extensions <- as.character(sapply(imgs, file_extension))
    all_valid <- extensions %in% valid_extens
    if (any(!all_valid)) {
      cli::cli_warn("Image{?s} {.val {imgs[!all_valid]}} of invalid format ignored.")
    }

    imgs <- paste0(path, "/", imgs[all_valid])
    list_img <-
      lapply(imgs, function(x){
        EBImage::readImage(x)
      })
    names(list_img) <- basename(imgs)
    if(isTRUE(plot)){
      image_combine(list_img, nrow = nrow, ncol = ncol)
    }
    if(resize != FALSE){
      if(!is.numeric(resize)){
        cli::cli_abort("Argument {.val resize} must be numeric.")
      }
      list_img <- image_resize(list_img, resize, parallel = FALSE)
    }
    invisible(list_img)
  } else{
    img_dir <- ifelse(is.null(path), file_dir(img), path)
    all_files <- sapply(list.files(img_dir), file_name)
    img_name <- file_name(img)
    test <- img_name %in% file_name(list.files(img_dir))
    if(!any(grepl("http", img_dir, fixed = TRUE)) & !all(test)){
      cli::cli_abort("Image {.val {img_name[which(test == FALSE)]}} not found in {.dir {img_dir[which(test == FALSE)]}}.")
    }
    fext <- file_extension(img)
    img_name <- paste0(img_dir, "/", img_name , ".", fext[length(fext)])
    if(length(img) > 1){
      ls <-
        lapply(seq_along(img_name),
               function(x){
                 fext <- file_extension(img_name[[1]])
                 if(fext[length(fext)] %in% c("tif", "TIF", "tiff", "TIFF", "gri", "grd")){
                   terra::rast(img_name[x])
                 } else{
                   EBImage::readImage(img_name[x], ...)
                 }

               })
      names(ls) <- basename(img_name)
      if(isTRUE(plot)){
        image_combine(ls, nrow = nrow, ncol = ncol)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          cli::cli_abort("Argument {.val resize} must be numeric.")
        }
        ls <- image_resize(ls, resize)
      }
      invisible(ls)
    } else{
      fext <- file_extension(img_name)
      if(fext[length(fext)] %in% c("tif", "TIF", "tiff", "TIFF", "gri", "grd")){
        img <- terra::rast(img_name)
      } else{
        img <- EBImage::readImage(img_name, ...)
      }
      if(isTRUE(plot)){
        plot(img)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          cli::cli_abort("Argument {.val resize} must be numeric.")
        }
        img <- image_resize(img, resize)
      }
      invisible(img)
    }
  }
}

#' @export
#' @name utils_image
image_export <- function(img,
                         name,
                         prefix = "",
                         extension = NULL,
                         subfolder = NULL,
                         ...){
  check_ebi()
  if(class(img) %in% c("binary_list", "index_list",
                       "img_mat_list", "palette_list")){
    img <- lapply(img, function(x){x[[1]]})
  }
  if(inherits(img, "segment_list")){
    img <- lapply(img, function(x){x[[1]][[1]]})
  }
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    name <- file_name(names(img))
    extens <- file_extension(names(img))
    if(any(sapply(extens, length)) ==  0 & is.null(extension)){
      extens <- rep("jpg", length(img))
      cli::cli_inform(c(
        "v" = "Image{?s} exported as {.val *.jpg} file{?s}."
      ))

    }
    if(!is.null(extension)){
      extens <- rep(extension, length(img))
    }
    names(img) <- paste0(name, ".", extens)
    if(!missing(subfolder)){
      dir_out <- paste0(getwd(), "/", subfolder)
      if(dir.exists(dir_out) == FALSE){
        dir.create(dir_out, recursive = TRUE)
      }
      names(img) <- paste0(dir_out, "/", prefix, name, ".", extens)
      a <-
        lapply(seq_along(img), function(i){
          EBImage::writeImage(x = img[[i]], files = names(img[i]), ...)
        })
    } else{
      a <-
        lapply(seq_along(img), function(i){
          EBImage::writeImage(x = img[[i]], files = paste0(prefix, names(img[i])), ...)
        })
    }

  } else{
    filname <- file_name(name)
    extens <- unlist(file_extension(name))
    dir_out <- file_dir(name)
    if(length(extens) ==  1){
      extens <- extens
    } else if(length(extens) ==  0 & is.null(extension)){
      extens <- "jpg"
      cli::cli_inform(c(
        "v" = "Image{?s} exported as {.val *.jpg} file{?s}."
      ))
    } else if(!is.null(extension)){
      extens <- extension
    }
    if(!missing(subfolder) & nchar(dir_out) == 2){
      dir_out <- paste0("./", subfolder)
    }
    if(dir.exists(dir_out) == FALSE){
      dir.create(dir_out, recursive = TRUE)
    }
    name <- paste0(dir_out, "/", filname, ".", extens)
    EBImage::writeImage(img, name)
  }
}

#' @export
#' @name utils_image
image_input <- function(img, ...){
  check_ebi()
  if(inherits(img, "character")){
    image_import(img, ...)
  } else if(inherits(img, "array")){
    range <- apply(img, 3, max)
    if(any(range > 1)){
      EBImage::Image(img / 255, colormode = "color")
    } else{
      EBImage::Image(img, colormode = "color")
    }
  }
}

#' @export
#' @name utils_image
image_pliman <- function(img, plot = FALSE){
  check_ebi()
  path <- system.file("tmp_images", package = "pliman")
  files <- list.files(path)
  if(!missing(img)){
    if(!img %in% files){
      cli::cli_abort(c(
        "!" = "Image not available in {.pkg pliman}.",
        "i" = "Available images: {.val {paste(files, collapse = ', ')}}"
      ))
    }
    im <- image_import(system.file(paste0("tmp_images/", img), package = "pliman"))
    if(isTRUE(plot)){
      plot(im)
    }
    invisible(im)
  } else{
    path
  }
}




##### Spatial transformations
#'Spatial transformations
#'
#' Performs image rotation and reflection
#' * `image autocrop()` Crops automatically  an image to the area of objects.
#' * `image_crop()` Crops an image to the desired area.
#' * `image_trim()` Remove pixels from the edges of an image (20 by default).
#' * `image_dimension()` Gives the dimension (width and height) of an image.
#' * `image_rotate()` Rotates the image clockwise by the given angle.
#' * `image_horizontal()` Converts (if needed) an image to a horizontal image.
#' * `image_vertical()` Converts (if needed) an image to a vertical image.
#' * `image_hreflect()` Performs horizontal reflection of the `image`.
#' * `image_vreflect()` Performs vertical reflection of the `image`.
#' * `image_resize()` Resize the `image`. See more at [EBImage::resize()].
#' * `image_contrast()` Improve contrast locally by performing adaptive
#' histogram equalization. See more at [EBImage::clahe()].
#' * `image_dilate()` Performs image dilatation. See more at [EBImage::dilate()].
#' * `image_erode()` Performs image erosion. See more at [EBImage::erode()].
#' * `image_opening()` Performs an erosion followed by a dilation. See more at
#' [EBImage::opening()].
#' * `image_closing()` Performs a dilation followed by an erosion. See more at
#' [EBImage::closing()].
#' * `image_filter()` Performs median filtering in constant time. See more at
#' [EBImage::medianFilter()].
#' * `image_blur()` Performs blurring filter of images. See more at
#' [EBImage::gblur()].
#' * `image_skeleton()` Performs image skeletonization.
#'
#'
#' @name utils_transform
#' @inheritParams image_view
#' @inheritParams analyze_objects
#' @param img An image or a list of images of class `Image`.
#' @param index The index to segment the image. See [image_index()] for more
#'   details. Defaults to `"NB"` (normalized blue).
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param show How to plot in mapview viewer, either `"rgb"` or `"index"`.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param edge
#' * for [image_autocrop()] the number of pixels in the edge of the cropped
#' image. If `edge = 0` the image will be cropped to create a bounding rectangle
#' (x and y coordinates) around the image objects.
#' * for [image_trim()], the number of pixels removed from the edges. By
#' default, 20 pixels are removed from all the edges.
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
#' @param top,bottom,left,right The number of pixels removed from `top`,
#'   `bottom`, `left`, and `right` when using [image_trim()].
#' @param angle The rotation angle in degrees.
#' @param bg_col Color used to fill the background pixels, defaults to `"white"`.
#' @param rel_size The relative size of the resized image. Defaults to 100. For
#'   example, setting `rel_size = 50` to an image of width `1280 x 720`, the new
#'   image will have a size of `640 x 360`.
#' @param width,height
#'  * For `image_resize()` the Width and height of the resized image. These arguments
#'   can be missing. In this case, the image is resized according to the
#'   relative size informed in `rel_size`.
#'  * For `image_crop()` a numeric vector indicating the pixel range (x and y,
#' respectively) that will be maintained in the cropped image, e.g., width =
#' 100:200
#' @param kern An `Image` object or an array, containing the structuring
#'   element. Defaults to a brushe generated with [EBImage::makeBrush()].
#' @param niter The number of iterations to perform in the thinning procedure.
#'   Defaults to 3. Set to `NULL` to iterate until the binary image is no longer
#'   changing.
#' @param shape A character vector indicating the shape of the brush. Can be
#'   `box`, `disc`, `diamond`, `Gaussian` or `line`. Default is `disc`.
#' @param size
#' * For `image_filter()` is the median filter radius (integer). Defaults to `3`.
#' * For `image_dilate()` and `image_erode()` is an odd number containing the
#' size of the brush in pixels. Even numbers are rounded to the next odd one.
#' The default depends on the image resolution and is computed as the image
#' resolution (megapixels) times 20.
#' @param sigma A numeric denoting the standard deviation of the Gaussian filter
#'   used for blurring. Defaults to `3`.
#' @param cache The the L2 cache size of the system CPU in kB (integer).
#'   Defaults to `512`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param plot If `TRUE` plots the modified image. Defaults to `FALSE`.
#' @param ... Additional arguments passed on to [image_binary()].
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_skeleton()` returns a binary `Image` object.
#' * All other functions returns a  modified version of `image` depending on the
#' `image_*()` function used.
#' * If `image` is a list, a list of the same length will be returned.
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'plot(img)
#'img <- image_resize(img, 50)
#'img1 <- image_rotate(img, 45)
#'img2 <- image_hreflect(img)
#'img3 <- image_vreflect(img)
#'img4 <- image_vertical(img)
#'image_combine(img1, img2, img3, img4)
#' }
image_autocrop <- function(img,
                           index = "NB",
                           edge = 5,
                           opening = 5,
                           closing = FALSE,
                           filter = FALSE,
                           invert = FALSE,
                           threshold = "Otsu",
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                         "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.4), workers)
      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # verbose header with cli
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Image processing using {nworkers} workers"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%S')}")
        )
      }

      # run image_autocrop in parallel with built-in progress
      res <- mirai::mirai_map(
        .x       = img,
        .f       = function(image) image_autocrop(image, index, edge),
        .promise = if (verbose) cli::cli_progress_update
      )[.progress]

      # final completion message
      if (verbose) {
        cli::cli_rule(
          left = cli::col_green("All {length(img)} images processed")
        )
      }
    } else{
      res <- lapply(img, image_autocrop, index, edge)
    }
    invisible(structure(res, class = "autocrop_list"))
  } else{
    conv_hull <- object_coord(img,
                              index = index,
                              id = NULL,
                              edge = edge,
                              plot = FALSE,
                              opening = opening,
                              closing = closing,
                              filter = filter,
                              invert = invert,
                              threshold = threshold)
    segmented <- img[conv_hull[1]:conv_hull[2],
                     conv_hull[3]:conv_hull[4],
                     1:3]
    if(isTRUE(plot)){
      plot(segmented)
    }
    invisible(segmented)
  }
}
#' @name utils_transform
#' @export

image_crop <- function(img,
                       width = NULL,
                       height = NULL,
                       viewer = get_pliman_viewer(),
                       downsample = NULL,
                       max_pixels = 1000000,
                       show = "rgb",
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                         "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.4), workers)
      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # optional verbose message
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Parallel processing using {nworkers} cores"),
          right = cli::col_blue("Started on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(names_plant)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

      # run image_crop in parallel using mirai
      raw <- mirai::mirai_map(
        .x = img,
        .f = function(image) {
          pliman::image_crop(
            image,
            width,
            height,
            viewer,
            downsample,
            max_pixels
          )
        }
      )[.progress]
    } else{
      res <- lapply(img, image_crop, width, height, viewer, downsample, max_pixels)
    }
    invisible(res)
  } else{
    if (!is.null(width) | !is.null(height)) {
      dim <- dim(img)[1:2]
      if (!is.null(width)  & is.null(height)) {
        height <- 1:dim[2]
      }
      if (is.null(width) & !is.null(height)) {
        width <- 1:dim[1]
      }
      if(!is.null(height) & !is.null(width)){
        width <- width
        height <- height
      }
      if (!is.numeric(width) | !is.numeric(height)) {
        cli::cli_abort("Vectors {.val width} and {.val height} must be numeric.")
      }
      img@.Data <- img@.Data[width, height, ]
    }
    if (is.null(width) & is.null(height)) {
      if(vieweropt == "base"){
        cli::cli_inform(c("i" = "Use the {cli::col_blue(cli::style_bold('left mouse button'))} to crop the image."))

        if(EBImage::numberOfFrames(img) > 2){
          plot(EBImage::Image(img[,,1:3], colormode = "Color"))
        } else if(EBImage::numberOfFrames(img) == 1){
          plot(img)
        }
        cord <- locator(type = "p", n = 2, col = "red", pch = 19)
        minw <- min(cord$x[[1]], cord$x[[2]])
        maxw <- max(cord$x[[1]], cord$x[[2]])
        minh <- min(cord$y[[1]], cord$y[[2]])
        maxh <- max(cord$y[[1]], cord$y[[2]])
        w <- round(minw, 0):round(maxw, 0)
        h <- round(minh, 0):round(maxh, 0)
      } else{
        nc <- ncol(img)
        mv <- mv_rectangle(img, show = show, downsample = downsample, max_pixels = max_pixels)
        w <- round(min(mv[,1]):max(mv[,1]))
        h <- round((min(mv[,2]))):round(max(mv[,2]))
      }
      img@.Data <- img@.Data[w, h, ]
      if(isTRUE(verbose)){
        cat(paste0("width = ", w[1], ":", w[length(w)]), "\n")
        cat(paste0("height = ", h[1], ":", h[length(h)]), "\n")
      }
    }
    if (isTRUE(plot)) {
      if(EBImage::numberOfFrames(img) > 2){
        plot(EBImage::Image(img[,,1:3], colormode = "Color"))
      } else if(EBImage::numberOfFrames(img) == 1){
        plot(img)
      }
    }
    invisible(img)
  }
}


#' @name utils_transform
#' @export
image_dimension <- function(img,
                            parallel = FALSE,
                            workers = NULL,
                            verbose = TRUE){
  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Parallel dimension extraction of {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
      }

      raw <- mirai::mirai_map(
        .x = img,
        .f = function(image) {
          pliman::image_dimension(image, verbose = FALSE)
        }
      )[.progress]

      res <- as.data.frame(do.call(rbind, raw))

    } else {
      res <- do.call(rbind, lapply(img, function(x) {
        dim <- image_dimension(x, verbose = FALSE)
        data.frame(width = dim[[1]], height = dim[[2]])
      }))
      res <- transform(res, image = rownames(res))[, c(3, 1, 2)]
      rownames(res) <- NULL
    }

    if (verbose) {
      cli::cli_rule("Image dimension summary")
      cli::cli_alert_info("Processed {.val {nrow(res)}} image(s)")
      print(res, row.names = FALSE)
    }

    invisible(res)

  } else {
    width  <- dim(img)[[1]]
    height <- dim(img)[[2]]

    if (verbose) {
      cli::cli_rule("Image dimension")
      cli::cli_text("Width : {.val {width}}")
      cli::cli_text("Height: {.val {height}}")
    }

    invisible(list(width = width, height = height))
  }
}

#' @name utils_transform
#' @export
image_rotate <- function(img,
                         angle,
                         bg_col = "white",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = TRUE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Rotating {length(img)} images in parallel"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          EBImage::rotate(im, angle, bg.col = bg_col)
        }
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left = cli::col_blue("Rotating {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images sequentially...",
          msg_done   = "Processing complete",
          msg_failed = "Sequential processing failed"
        )
      }

      res <- lapply(img, function(im) {
        EBImage::rotate(im, angle, bg.col = bg_col)
      })
    }

    if (isTRUE(plot)) {
      for (r in res) {
        if (EBImage::numberOfFrames(r) > 2) {
          plot(EBImage::Image(r[,,1:3], colormode = "Color"))
        }
      }
    }

    invisible(res)

  } else {
    rotated <- EBImage::rotate(img, angle, bg.col = bg_col)

    if (isTRUE(plot) && EBImage::numberOfFrames(rotated) > 2) {
      plot(EBImage::Image(rotated[,,1:3], colormode = "Color"))
    }

    invisible(rotated)
  }
}


#' @name utils_transform
#' @export
image_horizontal <- function(img,
                             parallel = FALSE,
                             workers = NULL,
                             verbose = TRUE,
                             plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Ensuring horizontal orientation of {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          w <- dim(im)[[1]]
          h <- dim(im)[[2]]
          if (w < h) {
            EBImage::rotate(im, 90)
          } else {
            im
          }
        }
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left = cli::col_blue("Processing images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images sequentially...",
          msg_done   = "Sequential processing complete",
          msg_failed = "Sequential processing failed"
        )
      }

      res <- lapply(img, function(im) {
        w <- dim(im)[[1]]
        h <- dim(im)[[2]]
        if (w < h) {
          EBImage::rotate(im, 90)
        } else {
          im
        }
      })
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    w <- dim(img)[[1]]
    h <- dim(img)[[2]]
    if (w < h) {
      img <- EBImage::rotate(img, 90)
    }

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_vertical <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Ensuring vertical orientation of {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          w <- dim(im)[[1]]
          h <- dim(im)[[2]]
          if (w > h) {
            EBImage::rotate(im, 90)
          } else {
            im
          }
        }
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Processing images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images sequentially...",
          msg_done   = "Sequential processing complete",
          msg_failed = "Sequential processing failed"
        )
      }

      res <- lapply(img, function(im) {
        w <- dim(im)[[1]]
        h <- dim(im)[[2]]
        if (w > h) {
          EBImage::rotate(im, 90)
        } else {
          im
        }
      })
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    w <- dim(img)[[1]]
    h <- dim(img)[[2]]
    if (w > h) {
      img <- EBImage::rotate(img, 90)
    }

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_hreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Horizontal reflection of {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg = "Processing images in parallel...",
          msg_done = "Reflection complete.",
          msg_failed = "Parallel reflection failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) EBImage::flop(im)
      )[.progress]

    } else {
      res <- lapply(img, EBImage::flop)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    img <- EBImage::flop(img)
    if (isTRUE(plot)) plot(img)
    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_vreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Vertical reflection of {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg = "Processing images in parallel...",
          msg_done = "Reflection complete.",
          msg_failed = "Parallel reflection failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) EBImage::flip(im)
      )[.progress]

    } else {
      res <- lapply(img, EBImage::flip)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    img <- EBImage::flip(img)
    if (isTRUE(plot)) plot(img)
    invisible(img)
  }
}


#' @name utils_transform
#' @export
image_resize <- function(img,
                         rel_size = 100,
                         width,
                         height,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE) {
  check_ebi()
  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    # define resize_fun fora do escopo de missing()
    resize_fun <- function(im, new_width, height) {
      EBImage::resize(im, new_width, height)
    }

    # calcula width/height antes
    res <- NULL
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Resizing {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg = "Resizing images in parallel...",
          msg_done = "Resize complete.",
          msg_failed = "Parallel resize failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im, width, height) {
          w <- dim(im)[[1]]
          new_width <- if (!missing(width)) width else w * rel_size / 100
          resize_fun(im, new_width, height)
        }
      )[.progress]

    } else {
      res <- lapply(img, function(im, width, height) {
        w <- dim(im)[[1]]
        new_width <- if (!missing(width)) width else w * rel_size / 100
        resize_fun(im, new_width, height)
      })
    }
    if(plot){
      image_combine(res)
    }
    invisible(res)

  } else {
    w <- dim(img)[[1]]
    new_width <- if (!missing(width)) width else w * rel_size / 100
    img <- EBImage::resize(img, new_width, height)
    if (isTRUE(plot)) plot(img)
    invisible(img)
  }
}



#' @name utils_transform
#' @export
image_trim <- function(img,
                       edge = NULL,
                       top = NULL,
                       bottom = NULL,
                       left = NULL,
                       right = NULL,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE) {
  check_ebi()

  # define bordas
  if (is.null(edge) && all(sapply(list(top, bottom, left, right), is.null))) {
    edge <- 20
  }
  if (is.null(edge) && !all(sapply(list(top, bottom, left, right), is.null))) {
    edge <- 0
  }

  top    <- ifelse(is.null(top), edge, top)
  bottom <- ifelse(is.null(bottom), edge, bottom)
  left   <- ifelse(is.null(left), edge, left)
  right  <- ifelse(is.null(right), edge, right)

  # processamento em lista
  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # modo paralelo
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Trimming {length(img)} images in parallel"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Trimming images...",
          msg_done   = "Trim completed.",
          msg_failed = "Parallel trimming failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          im <- im[, -c(1:top), ]
          im <- im[, -c((dim(im)[2] - bottom + 1):dim(im)[2]), ]
          im <- im[-c((dim(im)[1] - right + 1):dim(im)[1]), , ]
          im <- im[-c(1:left), , ]
          im
        }
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Trimming images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Trimming images...",
          msg_done   = "Trim completed.",
          msg_failed = "Sequential trimming failed."
        )
      }

      res <- lapply(img, function(im) {
        im <- im[, -c(1:top), ]
        im <- im[, -c((dim(im)[2] - bottom + 1):dim(im)[2]), ]
        im <- im[-c((dim(im)[1] - right + 1):dim(im)[1]), , ]
        im <- im[-c(1:left), , ]
        im
      })
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    img <- img[, -c(1:top), ]
    img <- img[, -c((dim(img)[2] - bottom + 1):dim(img)[2]), ]
    img <- img[-c((dim(img)[1] - right + 1):dim(img)[1]), , ]
    img <- img[-c(1:left), , ]

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_dilate <- function(img,
                         kern = NULL,
                         size = NULL,
                         shape = "disc",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # função auxiliar para cada imagem
    dilate_image <- function(im) {
      if (is.null(kern)) {
        d <- dim(im)
        s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
        s <- ifelse(s == 0, 2, s)
        k <- suppressWarnings(EBImage::makeBrush(s, shape = shape))
      } else {
        k <- kern
      }
      EBImage::dilate(im, k)
    }

    # modo paralelo
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Dilating {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Dilation complete.",
          msg_failed = "Parallel dilation failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = dilate_image
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Dilating {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images...",
          msg_done   = "Dilation complete.",
          msg_failed = "Sequential dilation failed."
        )
      }

      res <- lapply(img, dilate_image)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    if (is.null(kern)) {
      d <- dim(img)
      s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      s <- ifelse(s == 0, 2, s)
      kern <- suppressWarnings(EBImage::makeBrush(s, shape = shape))
    }
    img <- EBImage::dilate(img, kern)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_erode <- function(img,
                        kern = NULL,
                        size = NULL,
                        shape = "disc",
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # função auxiliar para cada imagem
    erode_image <- function(im) {
      if (is.null(kern)) {
        d <- dim(im)
        s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
        s <- ifelse(s == 0, 2, s)
        k <- suppressWarnings(EBImage::makeBrush(s, shape = shape))
      } else {
        k <- kern
      }
      EBImage::erode(im, k)
    }

    # modo paralelo
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Eroding {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Erosion complete.",
          msg_failed = "Parallel erosion failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = erode_image
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Eroding {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images...",
          msg_done   = "Erosion complete.",
          msg_failed = "Sequential erosion failed."
        )
      }

      res <- lapply(img, erode_image)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    if (is.null(kern)) {
      d <- dim(img)
      size <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    }

    img <- EBImage::erode(img, kern)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_opening <- function(img,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # função auxiliar para aplicar abertura
    opening_image <- function(im) {
      if (is.null(kern)) {
        d <- dim(im)
        s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
        s <- ifelse(s == 0, 2, s)
        k <- suppressWarnings(EBImage::makeBrush(s, shape = shape))
      } else {
        k <- kern
      }
      EBImage::opening(im, k)
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Opening {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Opening complete.",
          msg_failed = "Parallel opening failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = opening_image
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Opening {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images...",
          msg_done   = "Opening complete.",
          msg_failed = "Sequential opening failed."
        )
      }

      res <- lapply(img, opening_image)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    if (is.null(kern)) {
      d <- dim(img)
      size <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    }
    img <- EBImage::opening(img, kern)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_closing <- function(img,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # função auxiliar para aplicar fechamento
    closing_image <- function(im) {
      if (is.null(kern)) {
        d <- dim(im)
        s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
        s <- ifelse(s == 0, 2, s)
        k <- suppressWarnings(EBImage::makeBrush(s, shape = shape))
      } else {
        k <- kern
      }
      EBImage::closing(im, k)
    }

    # modo paralelo
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Closing {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Closing complete.",
          msg_failed = "Parallel closing failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = closing_image
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Closing {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images...",
          msg_done   = "Closing complete.",
          msg_failed = "Sequential closing failed."
        )
      }

      res <- lapply(img, closing_image)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    if (is.null(kern)) {
      d <- dim(img)
      size <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    }

    img <- EBImage::closing(img, kern)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_skeleton <- function(img,
                           kern = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE,
                           ...) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # função auxiliar para aplicar skeletonization
    skel_fun <- function(im) {
      if (EBImage::colorMode(im) != 0) {
        im <- help_binary(im, ..., resize = FALSE)
      }

      s <- matrix(1, nrow(im), ncol(im))
      skel <- matrix(0, nrow(im), ncol(im))
      k <- if (is.null(kern)) suppressWarnings(EBImage::makeBrush(2, shape = "diamond")) else kern

      while (max(s) == 1) {
        opened <- EBImage::opening(im, k)
        s <- im - opened
        skel <- skel | s
        im <- EBImage::erode(im, k)
      }

      EBImage::Image(skel)
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Skeletonizing {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing in parallel...",
          msg_done   = "Skeletonization complete.",
          msg_failed = "Skeletonization failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = skel_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Skeletonizing {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing...",
          msg_done   = "Skeletonization complete.",
          msg_failed = "Sequential skeletonization failed."
        )
      }

      res <- lapply(img, skel_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    if (EBImage::colorMode(img) != 0) {
      img <- help_binary(img, ..., resize = FALSE)
    }

    s <- matrix(1, nrow(img), ncol(img))
    skel <- matrix(0, nrow(img), ncol(img))
    kern <- if (is.null(kern)) suppressWarnings(EBImage::makeBrush(2, shape = "diamond")) else kern

    while (max(s) == 1) {
      opened <- EBImage::opening(img, kern)
      s <- img - opened
      skel <- skel | s
      img <- EBImage::erode(img, kern)
    }

    img <- EBImage::Image(skel)
    if (isTRUE(plot)) plot(img)
    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_thinning <- function(img,
                           niter = 3,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE,
                           ...) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    thin_fun <- function(im) {
      if (EBImage::colorMode(im) != 0) {
        im <- help_binary(im, ..., resize = FALSE)
      }

      if (is.null(niter)) {
        li <- sum(im)
        lf <- 1
        while ((li - lf) != 0) {
          li <- sum(im)
          im <- help_edge_thinning(im)
          lf <- sum(im)
        }
      } else {
        for (i in seq_len(niter)) {
          im <- help_edge_thinning(im)
        }
      }

      EBImage::Image(im)
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Thinning {.val {length(img)}} images"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing in parallel...",
          msg_done   = "Thinning complete.",
          msg_failed = "Thinning failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = thin_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Thinning {.val {length(img)}} images sequentially"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing...",
          msg_done   = "Thinning complete.",
          msg_failed = "Sequential thinning failed."
        )
      }

      res <- lapply(img, thin_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(invisible(res))

  } else {
    if (EBImage::colorMode(img) != 0) {
      img <- help_binary(img, ..., resize = FALSE)
    }

    if (is.null(niter)) {
      li <- sum(img)
      lf <- 1
      while ((li - lf) != 0) {
        li <- sum(img)
        img <- help_edge_thinning(img)
        lf <- sum(img)
      }
    } else {
      for (i in seq_len(niter)) {
        img <- help_edge_thinning(img)
      }
    }

    img <- EBImage::Image(img)
    if (isTRUE(plot)) plot(img)
    invisible(img)
  }
}



#' Perform Guo-Hall thinning on a binary image or list of binary images
#'
#' This function performs the Guo-Hall thinning algorithm (Guo and Hall, 1989)
#' on a binary image or a list of binary images.
#'
#' @param img The binary image or a list of binary images to be thinned. It can
#'   be either a single binary image of class 'Image' or a list of binary
#'   images.
#' @param parallel Logical, whether to perform thinning using multiple cores
#'   (parallel processing). If TRUE, the function will use multiple cores for
#'   processing if available. Default is FALSE.
#' @param workers Integer, the number of workers (cores) to use for parallel
#'   processing. If NULL (default), it will use 40% of available cores.
#' @param verbose Logical, whether to display progress messages during parallel
#'   processing. Default is TRUE.
#' @param plot Logical, whether to plot the thinned images. Default is FALSE.
#' @param ... Additional arguments to be passed to [image_binary()] if
#'   \code{img} is not a binary image.
#'
#' @references Guo, Z., and R.W. Hall. 1989. Parallel thinning with
#'    two-subiteration algorithms. Commun. ACM 32(3): 359–373.
#'    \doi{10.1145/62065.62074}
#' @return If \code{img} is a single binary image, the function returns the
#'   thinned binary image. If \code{img} is a list of binary images, the
#'   function returns a list containing the thinned binary images.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("potato_leaves.jpg", plot = TRUE)
#' image_thinning_guo_hall(img, index = "R", plot = TRUE)
#' }
#'
image_thinning_guo_hall <- function(img,
                                    parallel = FALSE,
                                    workers = NULL,
                                    verbose = TRUE,
                                    plot = FALSE,
                                    ...) {
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    thin_fun <- function(im) {
      if (EBImage::colorMode(im) != 0) {
        im <- help_binary(im, ..., resize = FALSE)
      }
      helper_guo_hall(im)
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Guo-Hall thinning {.val {length(img)}} images"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing in parallel...",
          msg_done   = "Thinning complete.",
          msg_failed = "Thinning failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = thin_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Guo-Hall thinning {.val {length(img)}} images (sequential)"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing...",
          msg_done   = "Thinning complete.",
          msg_failed = "Sequential thinning failed."
        )
      }

      res <- lapply(img, thin_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(invisible(res))

  } else {
    if (EBImage::colorMode(img) != 0) {
      img <- help_binary(img, ..., resize = FALSE)
    }

    thin <- helper_guo_hall(img)

    if (isTRUE(plot)) {
      plot(thin)
    }

    invisible(thin)
  }
}

#' @name utils_transform
#' @export
image_filter <- function(img,
                         size = 2,
                         cache = 512,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE) {
  check_ebi()

  if (size < 2) {
    cli::cli_abort("Using {.arg size} < 2 may crash the R session. Use 2 or more.")
  }

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    filter_fun <- function(im) {
      EBImage::medianFilter(im, size, cache)
    }
    cli::cli_rule(
      left  = cli::col_blue("Median filtering {.val {length(img)}} images"),
      right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
    )
    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_progress_step(
          msg        = "Processing in parallel...",
          msg_done   = "Filtering complete.",
          msg_failed = "Parallel filtering failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = filter_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_progress_step(
          msg        = "Processing...",
          msg_done   = "Filtering complete.",
          msg_failed = "Sequential filtering failed."
        )
      }

      res <- lapply(img, filter_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(invisible(res))

  } else {
    img <- EBImage::medianFilter(img, size, cache)
    if (isTRUE(plot)) {
      plot(img)
    }
    invisible(img)
  }
}

#' @name utils_transform
#' @export
image_blur <- function(img,
                       sigma = 3,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  check_ebi()

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    blur_fun <- function(im) {
      EBImage::gblur(im, sigma)
    }

    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Blurring {.val {length(img)}} images"),
        right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
      )
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_progress_step(
          msg        = "Processing in parallel...",
          msg_done   = "Blurring complete.",
          msg_failed = "Parallel blurring failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = blur_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_progress_step(
          msg        = "Processing...",
          msg_done   = "Blurring complete.",
          msg_failed = "Sequential blurring failed."
        )
      }

      res <- lapply(img, blur_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(invisible(res))

  } else {
    img <- EBImage::gblur(img, sigma)
    if (isTRUE(plot)) {
      plot(img)
    }
    invisible(img)
  }
}
#' @name utils_transform
#' @export
image_contrast <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {
  check_ebi()

  get_factors <- function(x) {
    factors <- vector()
    for (i in 1:x) {
      if ((x %% i) == 0) {
        factors[i] <- i
      }
    }
    factors[!is.na(factors)]
  }

  contrast_fun <- function(im) {
    img_width <- dim(im)[1]
    img_height <- dim(im)[2]

    fx <- get_factors(img_width)
    nx <- suppressWarnings(fx[max(which(fx > 1 & fx < 100))])

    fy <- get_factors(img_height)
    ny <- suppressWarnings(fy[max(which(fy > 1 & fy < 100))])

    testx <- !any(fx > 1 & fx < 100)
    if (testx) {
      while (testx) {
        img_width <- img_width + 1
        fx <- get_factors(img_width)
        testx <- !any(fx > 1 & fx < 100)
        if (any(fx) > 100) break
      }
      im <- EBImage::resize(im, w = img_width, h = img_height)
      nx <- suppressWarnings(fx[max(which(fx > 1 & fx < 100))])
    }

    testy <- !any(fy > 1 & fy < 100)
    if (testy) {
      while (testy) {
        img_height <- img_height + 1
        fy <- get_factors(img_height)
        testy <- !any(fy > 1 & fy < 100)
        if (any(fy) > 100) break
      }
      im <- EBImage::resize(im, w = img_width, h = img_height)
      ny <- suppressWarnings(fy[max(which(fy > 1 & fy < 100))])
    }

    EBImage::clahe(im, nx = nx, ny = ny, bins = 256)
  }

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    if (verbose) {
      cli::cli_rule(
        left = cli::col_blue("Contrast adjustment"),
        right = cli::col_blue("Processing {length(img)} images")
      )
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers),
                         trunc(parallel::detectCores() * 0.4),
                         workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_progress_step(
          msg        = "Running contrast enhancement in parallel...",
          msg_done   = "Contrast enhancement complete.",
          msg_failed = "Contrast enhancement failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = contrast_fun
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_progress_step(
          msg        = "Running contrast enhancement...",
          msg_done   = "Contrast enhancement complete.",
          msg_failed = "Contrast enhancement failed."
        )
      }

      res <- lapply(img, contrast_fun)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(invisible(res))

  } else {
    img <- contrast_fun(img)
    if (isTRUE(plot)) {
      plot(img)
    }
    invisible(img)
  }
}

#' Create an `Image` object of a given color
#'
#' image_create() can be used to create an `Image` object with a desired color and size.
#'
#' @param color either a color name (as listed by [grDevices::colors()]), or a hexadecimal
#'   string of the form `"#rrggbb"`.
#' @param width,heigth The width and heigth of the image in pixel units.
#' @param plot Plots the image after creating it? Defaults to `FALSE`.
#'
#' @return An object of class `Image`.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' image_create("red")
#' image_create("#009E73", width = 300, heigth = 100)
#' }

image_create <- function(color,
                         width = 200,
                         heigth = 200,
                         plot = FALSE){
  check_ebi()
  width <- as.integer(width)
  heigth <- as.integer(heigth)
  rgb <- col2rgb(color) / 255
  r <- rep(rgb[1], width*heigth)
  g <- rep(rgb[2], width*heigth)
  b <- rep(rgb[3], width*heigth)
  img <- EBImage::Image(c(r, g, b),
                        dim = c(width, heigth, 3),
                        colormode = "color")
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}

#' Creates a binary image
#'
#' Reduce a color, color near-infrared, or grayscale images to a binary image
#' using a given color channel (red, green blue) or even color indexes. The
#' Otsu's thresholding method (Otsu, 1979) is used to automatically perform
#' clustering-based image thresholding.
#' @inheritParams image_index
#' @param img An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. See the available indexes with
#'   [pliman_indexes()] and [image_index()] for more details.
#' @param threshold The theshold method to be used.
#'  * By default (`threshold = "Otsu"`), a threshold value based
#'  on Otsu's method is used to reduce the grayscale image to a binary image. If
#'  a numeric value is informed, this value will be used as a threshold.
#'
#'  * If `threshold = "adaptive"`, adaptive thresholding (Shafait et al. 2008)
#'  is used, and will depend on the `k` and `windowsize` arguments.
#'
#'  * If any non-numeric value different than `"Otsu"` and `"adaptive"` is used,
#'  an iterative section will allow you to choose the threshold based on a
#'  raster plot showing pixel intensity of the index.
#' @param k a numeric in the range 0-1. when `k` is high, local threshold
#'   values tend to be lower. when `k` is low, local threshold value tend to be
#'   higher.
#' @param windowsize windowsize controls the number of local neighborhood in
#'   adaptive thresholding. By default it is set to `1/3 * minxy`, where
#'   `minxy` is the minimum dimension of the image (in pixels).
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If `TRUE`, pixels that have R, G, and B values equals to 1 will be
#'   considered as `NA`. This may be useful to compute an image index for
#'   objects that have, for example, a white background. In such cases, the
#'   background will not be considered for the threshold computation.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value as the percentage of desired resizing. For example, if
#'   `resize = 30`, the resized image will have 30% of the size of original
#'   image.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#'
#' @param erode,dilate,opening,closing,filter **Morphological operations (brush size)**
#'  * `dilate` puts the mask over every background pixel, and sets it to
#'  foreground if any of the pixels covered by the mask is from the foreground.
#'  * `erode` puts the mask over every foreground pixel, and sets it to
#'  background if any of the pixels covered by the mask is from the background.
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
#' @param invert Inverts the binary image, if desired.
#' @param plot Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references
#' Otsu, N. 1979. Threshold selection method from gray-level histograms. IEEE
#' Trans Syst Man Cybern SMC-9(1): 62–66. \doi{10.1109/tsmc.1979.4310076}
#'
#' Shafait, F., D. Keysers, and T.M. Breuel. 2008. Efficient implementation of
#' local adaptive thresholding techniques using integral images. Document
#' Recognition and Retrieval XV. SPIE. p. 317–322 \doi{10.1117/12.767755}
#'
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing binary images. The length will depend on the number
#'   of indexes used.
#' @importFrom utils read.csv
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg")
#'image_binary(img, index = c("R, G"))
#' }
#'
image_binary <- function(img,
                         index = "R",
                         r = 1,
                         g = 2,
                         b = 3,
                         re = 4,
                         nir = 5,
                         return_class = "ebimage",
                         threshold = c("Otsu", "adaptive"),
                         k = 0.15,
                         windowsize = NULL,
                         has_white_bg = FALSE,
                         resize = FALSE,
                         fill_hull = FALSE,
                         erode = FALSE,
                         dilate = FALSE,
                         opening = FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         invert = FALSE,
                         plot = TRUE,
                         nrow = NULL,
                         ncol = NULL,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE) {
  check_ebi()
  threshold <- threshold[[1]]

  bin_img <- function(imgs) {
    if(threshold == "adaptive"){
      if(is.null(windowsize)){
        windowsize <- min(dim(imgs)) / 3
        if(windowsize %% 2 == 0) windowsize <- as.integer(windowsize + 1)
      }
      if (windowsize <= 2) {
        cli::cli_abort("{.arg windowsize} must be >= 3")
      }
      if (windowsize %% 2 == 0) windowsize <- as.integer(windowsize + 1)
      if (windowsize >= dim(imgs)[[1]] || windowsize >= dim(imgs)[[2]]) {
        windowsize <- min(dim(imgs)) / 3
      }
      if (k > 1) {
        cli::cli_abort("{.arg k} must be in [0, 1].")
      }
      imgs <- EBImage::thresh(imgs, w=windowsize, h=windowsize, offset=k)
      # imgs <- EBImage::Image(threshold_adaptive(as.matrix(imgs), k, windowsize, 0.5))
    } else {
      if(threshold == "Otsu"){
        threshold_val <- help_otsu(imgs@.Data[!is.infinite(imgs@.Data) & !is.na(imgs@.Data)])
      } else if(is.numeric(threshold)) {
        threshold_val <- threshold
      } else {
        pixels <- terra::rast(t(imgs@.Data))
        terra::plot(pixels, col = custom_palette(n = 100), axes = FALSE, asp = NA)
        threshold_val <- readline("Selected threshold: ")
      }
      imgs <- EBImage::Image(imgs < threshold_val)
    }

    if(invert) imgs <- 1 - imgs
    imgs[is.na(imgs)] <- FALSE
    if (is.numeric(erode) & erode > 0) {
      imgs <- image_erode(imgs, size = erode)
    }
    if (is.numeric(dilate) & dilate > 0) {
      imgs <- image_dilate(imgs, size = dilate)
    }
    if (is.numeric(opening) & opening > 0) {
      imgs <- image_opening(imgs, size = opening)
    }
    if (is.numeric(closing) & closing > 0) {
      imgs <- image_closing(imgs, size = closing)
    }
    if (is.numeric(filter) & filter > 1) {
      imgs <- EBImage::medianFilter(imgs, filter)
    }
    if (isTRUE(fill_hull)) {
      imgs <- EBImage::fillHull(imgs)
    }
    invisible(imgs)
  }

  process_image <- function(im) {
    imgs <- lapply(
      image_index(im, index, r, g, b, re, nir, return_class,
                  resize, re, nir, has_white_bg, plot = FALSE,
                  nrow, ncol, verbose = verbose),
      bin_img
    )
    imgs
  }

  if(is.list(img) && all(sapply(img, inherits, what = "Image"))){
    cli::cli_rule("Image binarization")
    if (isTRUE(parallel)) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * .4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)
      cli::cli_progress_step(
        msg = "Processing {.val {length(img)}} images in parallel...",
        msg_done = "Batch processing finished",
        msg_failed = "Oops, something went wrong."
      )
      res <- mirai::mirai_map(
        .x = img,
        .f = process_image
      )[.progress]
    } else {
      res <- lapply(img, process_image)
    }
    invisible(structure(res, class = "binary_list"))
  } else {
    imgs <- process_image(img)
    if (isTRUE(plot)) {
      num_plots <- length(imgs)
      if (is.null(nrow) && is.null(ncol)){
        ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
        nrow <- ceiling(num_plots/ncol)
      }
      if (is.null(ncol)) ncol <- ceiling(num_plots/nrow)
      if (is.null(nrow)) nrow <- ceiling(num_plots/ncol)
      op <- par(mfrow = c(nrow, ncol))
      on.exit(par(op))
      index_names <- names(imgs)
      for(i in seq_along(imgs)){
        plot(imgs[[i]])
        if(verbose){
          dim <- image_dimension(imgs[[i]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index_names[[i]], pos = 4, col = "red")
        }
      }
    }
    invisible(imgs)
  }
}


#' Image indexes
#'
#' `image_index()` Builds image indexes using Red, Green, Blue, Red-Edge, and
#' NIR bands. See [this
#' page](https://nepem-ufsc.github.io/pliman/articles/indexes.html) for a
#' detailed list of available indexes.
#'
#'
#' @name image_index
#' @inheritParams plot_index
#' @param img An `Image` object. Multispectral mosaics can be converted to an
#'   `Image` object using `mosaic_as_ebimage()`.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to a binary image. Use [pliman_indexes()] or the
#'   `details` section to see the available indexes. Defaults to `NULL`
#'   (normalized Red, Green, and Blue). You can also use "RGB" for RGB only,
#'   "NRGB" for normalized RGB,  "MULTISPECTRAL" for multispectral indices
#'   (provided NIR and RE bands are available) or "all" for all indexes. Users
#'   can also calculate their own index using the band names, e.g., `index =
#'   "R+B/G"`.
#' @param r,g,b,re,nir The red, green, blue, red-edge, and near-infrared bands
#'   of the image, respectively. Defaults to 1, 2, 3, 4, and 5, respectively. If
#'   a multispectral image is provided (5 bands), check the order of bands,
#'   which are frequently presented in the 'BGR' format.
#' @param return_class The class of object to be returned. If `"terra` returns a
#'   SpatRaster object with the number of layers equal to the number of indexes
#'   computed. If `"ebimage"` (default) returns a list of `Image` objects, where
#'   each element is one index computed.
#' @param resize Resize the image before processing? Defaults to `resize =
#'   FALSE`. Use `resize = 50`, which resizes the image to 50% of the original
#'   size to speed up image processing.
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If TRUE, pixels that have R, G, and B values equals to 1 will be considered
#'   as NA. This may be useful to compute an image index for objects that have,
#'   for example, a white background. In such cases, the background will not be
#'   considered for the threshold computation.
#' @param plot Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param ... Additional arguments passed on to [plot.image_index()].
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references
#' Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#'
#' Karcher, D.E., and M.D. Richardson. 2003. Quantifying Turfgrass Color Using
#' Digital Image Analysis. Crop Science 43(3): 943–951.
#' \doi{10.2135/cropsci2003.9430}
#'
#' Bannari, A., D. Morin, F. Bonn, and A.R. Huete. 1995. A review of vegetation
#' indices. Remote Sensing Reviews 13(1–2): 95–120.
#' \doi{10.1080/02757259509532298}
#'
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing Grayscale images. The length will depend on the
#'   number of indexes used.
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg")
#'image_index(img, index = c("R, NR"))
#' }
image_index <- function(img,
                        index = NULL,
                        r = 1,
                        g = 2,
                        b = 3,
                        re = 4,
                        nir = 5,
                        return_class = c("ebimage", "terra"),
                        resize = FALSE,
                        has_white_bg = FALSE,
                        plot = TRUE,
                        nrow = NULL,
                        ncol = NULL,
                        max_pixels = 100000,
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        ...){
  check_ebi()
  return_classopt <- c("terra", "ebimage")
  return_classopt <- return_classopt[pmatch(return_class[1], return_classopt)]
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.cls Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)
      cli::cli_progress_step(
        msg = "Processing {.val {length(img)}} images in parallel...",
        msg_done = "Image index extraction finished",
        msg_failed = "Something went wrong during image index extraction."
      )
      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          image_index(im, index, r, g, b, re, nir, resize, has_white_bg, plot, nrow, ncol, max_pixels)
        }
      )[.progress]

    } else{
      res <- lapply(img, image_index, index, r, g, b, re, nir, resize, has_white_bg, plot, nrow, ncol, max_pixels)
    }
    invisible(structure(res, class = "index_list"))
  } else{
    if(resize != FALSE){
      img <- image_resize(img, resize)
    }
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
    hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])
    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    }else{
      RE <- try(img@.Data[,,re], TRUE)
      NIR <- try(img@.Data[,,nir], TRUE)
      test_multi <- any(sapply(list(RE, NIR), class) == "try-error")
      if(isTRUE(test_multi)){
        all_ind <- ind$Index[!ind$Index %in% nir_ind]
      } else{
        all_ind <- ind$Index
      }
      if(index[[1]] %in% c("RGB", "NRGB", "MULTISPECTRAL", "all")){
        index <-  switch (index,
                          RGB = c("R", "G", "B"),
                          NRGB = c("NR", "NG", "NB"),
                          MULTISPECTRAL = c("NDVI", "PSRI", "GNDVI", "RVI", "NDRE", "TVI", "CVI", "EVI", "CIG", "CIRE", "DVI", "NDWI"),
                          all = all_ind
        )} else{
          if(length(index) > 1){
            index <- index
          } else{
            index <- strsplit(index, "\\s*(,)\\s*")[[1]]
          }
        }
    }

    R <- try(img@.Data[,,r], TRUE)
    G <- try(img@.Data[,,g], TRUE)
    B <- try(img@.Data[,,b], TRUE)
    test_band <- any(sapply(list(R, G, B), class) == "try-error")
    if(any(index %in% hsb_ind)){
      hsb <- rgb_to_hsb(data.frame(R = c(R), G = c(G), B = c(B)))
      h <- matrix(hsb$h, nrow = nrow(img), ncol = ncol(img))
      s <- matrix(hsb$s, nrow = nrow(img), ncol = ncol(img))
      b <- matrix(hsb$b, nrow = nrow(img), ncol = ncol(img))
    }
    if(any(index %in% nir_ind)){
      if(isTRUE(test_multi)){
        cli::cli_abort("Near-Infrared and RedeEdge bands are not available in the provided image.")
      }
    }
    if(isTRUE(test_band)){
      cli::cli_abort("At least 3 bands (RGB) are necessary to calculate indices available in pliman.")
    }
    imgs <- list()
    for(i in 1:length(index)){
      indx <- index[[i]]
      if(!indx %in% ind$Index){
        if (isTRUE(verbose)) {
          cli::cli_inform(c("i" = "Index {.val {indx}} is not available. Trying to compute your own index."))
        }

      }
      if(isTRUE(has_white_bg)){
        R[which(R == 1 & G == 1 & B == 1)] <- NA
        G[which(R == 1 & G == 1 & B == 1)] <- NA
        B[which(R == 1 & G == 1 & B == 1)] <- NA
      }

      if(indx %in% ind$Index){
        imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==indx]))))
      } else{
        imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(indx))))
      }
    }
    names(imgs) <- index
    class(imgs) <- "image_index"
    if(plot == TRUE){
      plot_index(imgs, nrow = nrow, ncol = ncol, max_pixels = max_pixels, ...)
    }
    if(return_classopt == "ebimage"){
      invisible(imgs)
    } else{
      terras <-
        terra::rast(
          lapply(1:length(imgs), function(i){
            terra::rast(t(imgs[[i]]@.Data))
          }
          )
        )
      names(terras) <- names(imgs)
      invisible(terras)
    }
  }
}


#' Plots an `image_index` object
#'
#' The S3 method `plot()` can be used to generate a raster or density plot of
#' the index values computed with `image_index()`
#'
#' @details When `type = "raster"` (default), the function calls [plot_index()]
#' to create a raster plot for each index present in `x`. If `type = "density"`,
#' a for loop is used to create a density plot for each index. Both types of
#' plots can be arranged in a grid controlled by the `ncol` and `nrow`
#' arguments.
#'
#'
#' @name image_index
#' @param x An object of class `image_index`.
#' @param type The type of plot. Use `type = "raster"` (default) to produce a
#'   raster plot showing the intensity of the pixels for each image index or
#'   `type = "density"` to produce a density plot with the pixels' intensity.
#' @param ... Additional arguments passed to [plot_index()] for customization.
#' @method plot image_index
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A `NULL` object
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' # Example for S3 method plot()
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' # compute the index
#' ind <- image_index(img, index = c("R, G, B, NGRDI"), plot = FALSE)
#' plot(ind)
#'
#' # density plot
#' plot(ind, type = "density")
#' }
#'
#
plot.image_index <- function(x,
                             type = c("raster", "density"),
                             nrow = NULL,
                             ncol = NULL,
                             ...){
  check_ebi()
  typeop <- c("raster", "density")
  typeop <- typeop[pmatch(type[1], typeop)]

  if(!typeop %in% c("raster", "density")){
    cli::cli_abort("`type` must be one of the 'raster' or 'density'. ")
  }
  if(typeop == "density"){
    mat <-
      as.data.frame(
        do.call(cbind,
                lapply(x, function(i){
                  as.vector(i)}
                ))
      )
    mat <- data.frame(mat[sample(1:nrow(mat), 70000, replace = TRUE),])
    colnames(mat) <- names(x)
    num_plots <- ncol(mat)

    if (is.null(nrow) && is.null(ncol)){
      ncols <- ceiling(sqrt(num_plots))
      nrows <- ceiling(num_plots/ncols)
    }
    if (is.null(ncol)){
      ncols <- ceiling(num_plots/nrows)
    }
    if (is.null(nrow)){
      nrows <- ceiling(num_plots/ncols)
    }
    op <- par(mfrow = c(nrows, ncols),
              mar = c(3, 2.5, 3, 3))
    on.exit(par(op))

    for (col in names(mat)) {
      density_data <- density(mat[[col]])  # Calculate the density for the column
      plot(density_data, main = col, col = "red", lwd = 2, xlab = NA, ylab = "Density")  # Create the density plot
    }

  } else{
    plot_index(x, ncol = ncol, nrow = nrow, ...)
  }
}



#' Image segmentation
#' @description
#' * `image_segment()` reduces a color, color near-infrared, or grayscale images
#' to a segmented image using a given color channel (red, green blue) or even
#' color indexes (See [image_index()] for more details). The Otsu's thresholding
#' method (Otsu, 1979) is used to automatically perform clustering-based image
#' thresholding.
#'
#' * `image_segment_iter()` Provides an iterative image segmentation, returning
#' the proportions of segmented pixels.
#'
#' @inheritParams image_binary
#' @inheritParams image_index
#' @param img An image object or a list of image objects.
#' @param index
#'  * For `image_segment()`, a character value (or a vector of characters)
#'  specifying the target mode for conversion to binary image. See the available
#'  indexes with [pliman_indexes()].  See [image_index()] for more details.
#' * For `image_segment_iter()` a character or a vector of characters with the
#' same length of `nseg`. It can be either an available index (described above)
#' or any operation involving the RGB values (e.g., `"B/R+G"`).
#' @param col_background The color of the segmented background. Defaults to
#'   `NULL` (white background).
#' @param na_background Consider the background as NA? Defaults to FALSE.
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If `TRUE`, pixels that have R, G, and B values equals to 1 will be
#'   considered as `NA`. This may be useful to compute an image index for
#'   objects that have, for example, a white background. In such cases, the
#'   background will not be considered for the threshold computation.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param invert Inverts the binary image, if desired. For
#'   `image_segmentation_iter()` use a vector with the same length of `nseg`.
#' @param plot Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param nseg The number of iterative segmentation steps to be performed.
#' @param ... Additional arguments passed on to `image_segment()`.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @export
#' @name image_segment
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_segment()` returns list containing `n` objects where `n` is the
#' number of indexes used. Each objects contains:
#'    * `image` an image with the RGB bands (layers) for the segmented object.
#'    * `mask` A mask with logical values of 0 and 1 for the segmented image.
#'
#' * `image_segment_iter()` returns a list with (1) a data frame with the
#' proportion of pixels in the segmented images and (2) the segmented images.
#'

#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg", plot = TRUE)
#'image_segment(img, index = c("R, G, B"))
#' }
#'
#'
image_segment <- function(img,
                          index = NULL,
                          r = 1,
                          g = 2,
                          b = 3,
                          re = 4,
                          nir = 5,
                          threshold = c("Otsu", "adaptive"),
                          k = 0.1,
                          windowsize = NULL,
                          col_background = NULL,
                          na_background = FALSE,
                          has_white_bg = FALSE,
                          fill_hull = FALSE,
                          erode = FALSE,
                          dilate = FALSE,
                          opening = FALSE,
                          closing = FALSE,
                          filter = FALSE,
                          invert = FALSE,
                          plot = TRUE,
                          nrow = NULL,
                          ncol = NULL,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  check_ebi()
  threshold <- threshold[[1]]
  if(inherits(img, "img_segment")){
    img <- img[[1]]
  }
  if(is.list(img)){
    if(!all(sapply(img, class)  %in% c("Image", "img_segment"))){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)

      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      cli::cli_progress_step(
        msg        = "Processing {.val {length(img)}} images in parallel...",
        msg_done   = "Image segmentation finished",
        msg_failed = "Something went wrong during image segmentation."
      )

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          image_segment(
            im, index, r, g, b, re, nir,
            threshold, k, windowsize,
            col_background, has_white_bg,
            fill_hull, erode, dilate,
            opening, closing, filter,
            invert, plot = plot, nrow, ncol
          )
        }
      )[.progress]
    } else{
      res <- lapply(img, image_segment, index, r, g, b, re, nir, threshold, k, windowsize, col_background, has_white_bg, fill_hull, erode, dilate, opening, closing, filter, invert, plot = plot, nrow, ncol)
    }
    invisible(structure(res, class = "segment_list"))
  } else{
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
    hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])
    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    }else{
      RE <- try(img@.Data[,,re], TRUE)
      NIR <- try(img@.Data[,,nir], TRUE)
      test_multi <- any(sapply(list(RE, NIR), class) == "try-error")
      if(isTRUE(test_multi)){
        all_ind <- ind$Index[!ind$Index %in% nir_ind]
      } else{
        all_ind <- ind$Index
      }
      if(index[[1]] %in% c("RGB", "NRGB", "MULTISPECTRAL", "all")){
        index <-  switch (index,
                          RGB = c("R", "G", "B"),
                          NRGB = c("NR", "NG", "NB"),
                          MULTISPECTRAL = c("NDVI", "PSRI", "GNDVI", "RVI", "NDRE", "TVI", "CVI", "EVI", "CIG", "CIRE", "DVI", "NDWI"),
                          all = all_ind
        )} else{
          index <- strsplit(index, "\\s*(,)\\s*")[[1]]
        }
    }
    imgs <- list()
    # color for background
    if (is.null(col_background)){
      col_background <- col2rgb("white") / 255
    } else{
      ifelse(is.character(col_background),
             col_background <- col2rgb(col_background) / 255,
             col_background <- col_background / 255)
    }
    for(i in 1:length(index)){
      imgmask <- img
      indx <- index[[i]]
      img2 <- help_binary(img,
                          index = indx,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          threshold = threshold,
                          k = k,
                          windowsize = windowsize,
                          has_white_bg = has_white_bg,
                          resize = FALSE,
                          fill_hull = fill_hull,
                          erode = erode,
                          dilate = dilate,
                          opening = opening,
                          closing = closing,
                          filter = filter,
                          invert = invert)
      ID <- which(img2@.Data == FALSE)
      if(!na_background){
        imgmask@.Data[,,1][ID] <- col_background[1]
        imgmask@.Data[,,2][ID] <- col_background[2]
        imgmask@.Data[,,3][ID] <- col_background[3]
        if(dim(img)[[3]] > 3){
          imgmask@.Data[,,4][ID] <- 1
          imgmask@.Data[,,5][ID] <- 1
        }
      } else{
        imgmask@.Data[,,1][ID] <- NA
        imgmask@.Data[,,2][ID] <- NA
        imgmask@.Data[,,3][ID] <- NA
        if(dim(img)[[3]] > 3){
          imgmask@.Data[,,4][ID] <- NA
          imgmask@.Data[,,5][ID] <- NA
        }
      }

      imgs[[i]] <- imgmask
    }
    names(imgs) <- index
    num_plots <- length(imgs)
    if (is.null(nrow) && is.null(ncol)){
      ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
      nrow <- ceiling(num_plots/ncol)
    }
    if (is.null(ncol)){
      ncol <- ceiling(num_plots/nrow)
    }
    if (is.null(nrow)){
      nrow <- ceiling(num_plots/ncol)
    }
    if(plot == TRUE){
      op <- par(mfrow = c(nrow, ncol))
      on.exit(par(op))
      for(i in 1:length(imgs)){
        tmps <- imgs[[i]][,,1:3]
        EBImage::colorMode(tmps) <- "Color"
        plot(tmps)
        if(verbose == TRUE){
          dim <- image_dimension(imgs[[i]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
        }
      }
    }
    if(length(imgs) == 1){
      invisible(imgs[[1]])
    } else{
      invisible(structure(imgs, class = "img_segment"))
    }
  }
}




#' @export
#' @name image_segment
image_segment_iter <- function(img,
                               nseg = 2,
                               index = NULL,
                               invert = NULL,
                               threshold = NULL,
                               k = 0.1,
                               windowsize = NULL,
                               has_white_bg = FALSE,
                               plot = TRUE,
                               verbose = TRUE,
                               nrow = NULL,
                               ncol = NULL,
                               parallel = FALSE,
                               workers = NULL,
                               ...){
  check_ebi()
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if(verbose){
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Image segmentation (iterative) finished",
          msg_failed = "Something went wrong during the segmentation process."
        )
      }
      a <- mirai::mirai_map(
        .x = img,
        .f = function(im){
          image_segment_iter(im, nseg, index, invert, threshold, has_white_bg, plot, verbose, nrow, ncol,  ...)
        }
      )[.progress]
    } else{
      a <- lapply(img, image_segment_iter, nseg, index, invert, threshold, has_white_bg, plot, verbose, nrow, ncol, ...)
    }
    results <-
      do.call(rbind, lapply(a, function(x){
        x$results
      }))
    images <-
      lapply(a, function(x){
        x$images
      })
    invisible(list(results = results,
                   images = images))
  } else{
    avali_index <- pliman_indexes()
    if(nseg == 1){
      if(is.null(invert)){
        invert <- FALSE
      } else{
        invert <- invert
      }
      if(is.null(threshold)){
        threshold <- "Otsu"
      } else{
        threshold <- threshold
      }
      if(is.null(index)){
        image_segment(img,
                      invert = invert[1],
                      index = "all",
                      has_white_bg = has_white_bg,
                      ...)
        index <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
      } else{
        index <- index[1]
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      segmented <-
        image_segment(img,
                      index = index,
                      threshold = my_thresh,
                      invert = invert[1],
                      plot = FALSE,
                      has_white_bg = has_white_bg,
                      ...)
      total <- length(img)
      segm <- length(which(segmented != 1))
      prop <- segm / total * 100
      results <- data.frame(total = total,
                            segmented = segm,
                            prop = prop)
      imgs <- list(img, segmented)
      if(verbose){
        print(results)
      }
      if(plot == TRUE){
        image_combine(imgs, ...)
      }
      invisible(list(results = results,
                     images = imgs))
    } else{
      if(is.null(index)){
        image_segment(img,
                      index = "all",
                      ...)
        indx <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
      } else{
        if(length(index) != nseg){
          cli::cli_abort("Length of {.code index} must be equal to {.code nseg}.")
        }
        indx <- index[1]
      }
      if(is.null(invert)){
        invert <- rep(FALSE, nseg)
      } else{
        invert <- invert
      }
      segmented <- list()
      total <- length(img)
      if(is.null(threshold)){
        threshold <- rep("Otsu", nseg)
      } else{
        threshold <- threshold
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      first <-
        image_segment(img,
                      index = indx,
                      invert = invert[1],
                      threshold = my_thresh[1],
                      plot = FALSE,
                      has_white_bg = has_white_bg,
                      ...)
      segmented[[1]] <- first
      for (i in 2:(nseg)) {
        if(is.null(index)){
          image_segment(first,
                        index = "all",
                        plot = TRUE,
                        has_white_bg = has_white_bg,
                        ncol = ncol,
                        nrow = nrow,
                        ...)
          indx <-
            switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                   "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                   "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                   "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
          if(is.null(indx)){
            break
          }
        } else{
          indx <- index[i]
        }
        my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[i]))),
                            as.character(threshold[i]),
                            as.numeric(threshold[i]))
        second <-
          image_segment(first,
                        index = indx,
                        threshold = my_thresh,
                        invert = invert[i],
                        plot = FALSE,
                        ...)
        segmented[[i]] <- second
        first <- second
      }
      pixels <-
        rbind(total,
              do.call(rbind,
                      lapply(segmented, function(x){
                        length(which(x != 1))
                      })
              )
        )
      rownames(pixels) <- NULL
      colnames(pixels) <- "pixels"
      prop <- NULL
      for(i in 2:nrow(pixels)){
        prop[1] <- 100
        prop[i] <- pixels[i] / pixels[i - 1] * 100
      }
      pixels <- data.frame(pixels)
      pixels$percent <- prop
      imgs <- lapply(segmented, function(x){
        x[[1]]
      })
      imgs <- c(list(img), segmented)
      names <- paste("seg", 1:length(segmented), sep = "")
      names(imgs) <- c("original", names)
      pixels <- transform(pixels, image = c("original",names))
      pixels <- pixels[,c(3, 1, 2)]
      if(verbose){
        print(pixels)
      }
      if(plot == TRUE){
        image_combine(imgs, ncol = ncol, nrow = nrow, ...)
      }
      invisible(list(results = pixels,
                     images = imgs))
    }
  }
}



#' Image segmentation using k-means clustering
#'
#' Segments image objects using clustering by the k-means clustering algorithm
#' @inheritParams image_segment
#' @inheritParams analyze_objects
#' @param img An `Image` object.
#' @param bands A numeric integer/vector indicating the RGB band used in the
#'   segmentation. Defaults to `1:3`, i.e., all the RGB bands are used.
#' @param nclasses The number of desired classes after image segmentation.
#' @param invert Invert the segmentation? Defaults to `FALSE`. If `TRUE` the
#'   binary matrix is inverted.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param plot Plot the segmented image?
#' @return A list with the following values:
#' * `image` The segmented image considering only two classes (foreground and
#' background)
#' * `clusters` The class of each pixel. For example, if `ncluster = 3`,
#' `clusters` will be a two-way matrix with values ranging from 1 to 3.
#' `masks` A list with the binary matrices showing the segmentation.
#' @export
#' @references Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A
#'   K-means clustering algorithm. Applied Statistics, 28, 100–108.
#'   \doi{10.2307/2346830}
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' img <- image_pliman("la_leaves.jpg", plot = TRUE)
#' seg <- image_segment_kmeans(img)
#' seg <- image_segment_kmeans(img, fill_hull = TRUE, invert = TRUE, filter = 10)
#' }

image_segment_kmeans <-   function (img,
                                    bands = 1:3,
                                    nclasses = 2,
                                    invert = FALSE,
                                    opening = FALSE,
                                    closing = FALSE,
                                    filter = FALSE,
                                    erode = FALSE,
                                    dilate = FALSE,
                                    fill_hull = FALSE,
                                    plot = TRUE){
  check_ebi()
  imm <- img@.Data[, , bands]
  if(length(dim(imm)) < 3){
    imb <- data.frame(B1 = image_to_mat(imm)[,3])
  } else{
    imb <- image_to_mat(imm)[, -c(1, 2)]
  }
  # rownames(imb) <- paste0("r", 1:nrow(imb))
  x <- suppressWarnings(stats::kmeans(na.omit(imb), nclasses))
  imm <- cbind(imb, 'clus'=NA)
  imm[names(x$cluster), ] <- x$cluster
  x2 <- x3 <- imm$clus
  nm <- names(sort(table(x2)))
  for (i in 1:length(nm)) {
    x3[x2 == nm[i]] <- i
  }
  m <- matrix(x3, nrow = dim(img)[1])
  LIST <- list()
  for (i in 1:length(nm)) {
    list <- list(m == i)
    LIST <- c(LIST, list)
  }
  if(isTRUE(fill_hull)){
    LIST <- lapply(LIST, EBImage::fillHull)
  }
  if(is.numeric(opening) & opening > 0){
    LIST <- lapply(LIST, image_opening, size = opening)
  }
  if(is.numeric(erode) & erode > 0){
    LIST <- lapply(LIST, image_erode, size = erode)
  }
  if(is.numeric(dilate) & dilate > 0){
    LIST <- lapply(LIST, image_dilate, size = dilate)
  }
  if(is.numeric(closing) & closing > 0){
    LIST <- lapply(LIST, image_closing, size = closing)
  }
  if(is.numeric(filter) & filter > 1){
    LIST <- lapply(LIST, EBImage::medianFilter, size = filter)
  }
  mask <- LIST[[1]]
  if(isFALSE(invert)){
    id <- which(mask == 1)
  } else{
    id <- which(mask != 1)
  }
  im2 <- img
  im2@.Data[, , 1][id] <- 1
  im2@.Data[, , 2][id] <- 1
  im2@.Data[, , 3][id] <- 1
  if(isTRUE(plot)){
    if(nclasses == 2){
      plot(im2)
    } else{
      suppressWarnings(image(m, useRaster = TRUE))
    }
  }
  invisible(list(img = im2,
                 clusters = m,
                 masks = LIST))
}


#' Image segmentation by hand
#'
#' This R code is a function that allows the user to manually segment an image based on the parameters provided. This only works in an interactive section.
#'
#' @details If the shape is "free", it allows the user to draw a perimeter to
#'   select/remove objects. If the shape is "circle", it allows the user to
#'   click on the center and edge of the circle to define the desired area. If
#'   the shape is "rectangle", it allows the user to select two points to define
#'   the area.
#'
#' @param img An `Image` object.
#' @param shape The type of shape to use. Defaults to "free". Other possible
#'   values are "circle" and "rectangle". Partial matching is allowed.
#' @param type The type of segmentation. By default (`type = "select"`) objects
#'   are selected. Use `type = "remove"` to remove the selected area from the
#'   image.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param resize By default, the segmented object is resized to fill the
#'   original image size. Use `resize = FALSE` to keep the segmented object in
#'   the original scale.
#' @param edge Number of pixels to add in the edge of the segmented object when
#'   `resize = TRUE`. Defaults to 5.
#' @param plot Plot the segmented object? Defaults to `TRUE`.
#'
#' @return A list with the segmented image and the mask used for segmentation.
#' @export
#'
#' @examples
#' if (interactive()) {
#' img <- image_pliman("la_leaves.jpg")
#' seg <- image_segment_manual(img)
#' plot(seg$mask)
#'
#' }
image_segment_manual <-  function(img,
                                  shape = c("free", "circle", "rectangle"),
                                  type = c("select", "remove"),
                                  viewer = get_pliman_viewer(),
                                  resize = TRUE,
                                  edge = 5,
                                  plot = TRUE){
  check_ebi()
  vals <- c("free", "circle", "rectangle")
  shape <- vals[[pmatch(shape[1], vals)]]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    if(shape == "free"){
      if(vieweropt == "base"){
        plot(img)
        cli::cli_inform(c(
          "i" = "Please, draw a perimeter to select/remove objects. Click {.key Esc} to finish."
        ))

        stop <- FALSE
        n <- 1e+06
        coor <- NULL
        a <- 0
        while (isFALSE(stop)) {
          if (a > 1) {
            if (nrow(coor) > 1) {
              lines(coor[(nrow(coor) - 1):nrow(coor), 1], coor[(nrow(coor) -
                                                                  1):nrow(coor), 2], col = "red")
            }
          }
          x = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
          if (is.null(x)){
            stop <- TRUE
          }
          coor <- rbind(coor, x)
          a <- a + 1
          if (a >= n) {
            stop = TRUE
          }
        }
        coor <- rbind(coor, coor[1, ])
      } else{
        coor <- mv_polygon(img)
        plot(img)

      }
    }

    if(shape == "circle"){
      if(vieweropt == "base"){
        plot(img)
        cli::cli_alert_info("Click on the center of the circle")
        cent = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
        cli::cli_alert_info("Click on the edge of the circle")
        ext = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
        radius = sqrt(sum((cent - ext)^2))
        x1 = seq(-1, 1, l = 2000)
        x2 = x1
        y1 = sqrt(1 - x1^2)
        y2 = (-1) * y1
        x = c(x1, x2) * radius + cent[1]
        y = c(y1, y2) * radius + cent[2]
      } else{
        mv <- mv_two_points(img)
        radius = sqrt(sum((c(mv$x1, mv$y1) - c(mv$x2, mv$y2))^2))
        x1 = seq(-1, 1, l = 2000)
        x2 = x1
        y1 = sqrt(1 - x1^2)
        y2 = (-1) * y1
        x = c(x1, x2) * radius + mv$x1
        y = c(y1, y2) * radius + mv$y1
      }
      coor = cbind(x, y)
      plot(img)
    }

    if(shape == "rectangle"){
      if(vieweropt == "base"){
        plot(img)
        cli::cli_alert_info("Select {.val 2} points drawing the diagonal that includes the area of interest.")
        cord <- unlist(locator(type = "p", n = 2, col = "red", pch = 19))
        coor <-
          rbind(c(cord[1], cord[3]),
                c(cord[2], cord[3]),
                c(cord[2], cord[4]),
                c(cord[1], cord[4]))
      } else{
        coor <- mv_rectangle(img)
        plot(img)
      }
    }
    mat <- NULL
    for (i in 1:(nrow(coor) - 1)) {
      c1<-  coor[i, ]
      c2 <- coor[i + 1, ]
      a <- c1[2]
      b <- (c2[2] - c1[2])/(c2[1] - c1[1])
      Xs <- round(c1[1], 0):round(c2[1], 0) - round(c1[1], 0)
      Ys <- round(a + b * Xs, 0)
      mat <- rbind(mat, cbind(Xs + round(c1[1], 0), Ys))
      lines(Xs + round(c1[1], 0), Ys, col = "red")
    }
    n = dim(img)
    imF = matrix(0, n[1], n[2])
    id = unique(mat[, 1])
    for (i in id) {
      coorr <- mat[mat[, 1] == i, ]
      imF[i, min(coorr[, 2], na.rm = T):max(coorr[, 2], na.rm = T)] = 1
    }
    mask <- EBImage::fillHull(EBImage::bwlabel(imF))
    # invisible(mask)
    if(type[1] == "select"){
      id <- mask != 1
    } else{
      id <- mask == 1
    }
    img@.Data[, , 1][id] = 1
    img@.Data[, , 2][id] = 1
    img@.Data[, , 3][id] = 1

    if(isTRUE(resize)){
      nrows <- nrow(mask)
      ncols <- ncol(mask)
      a <- apply(mask, 2, function(x) {
        any(x != 0)
      })
      col_min <- min(which(a == TRUE))
      col_min <- ifelse(col_min < 1, 1, col_min) - edge
      col_max <- max(which(a == TRUE))
      col_max <- ifelse(col_max > ncols, ncols, col_max) + edge
      b <- apply(mask, 1, function(x) {
        any(x != 0)
      })
      row_min <- min(which(b == TRUE))
      row_min <- ifelse(row_min < 1, 1, row_min) - edge
      row_max <- max(which(b == TRUE))
      row_max <- ifelse(row_max > nrows, nrows, row_max) + edge
      img <- img[row_min:row_max, col_min:col_max, 1:3]
    }
    if(isTRUE(plot)){
      plot(img)
    }
    invisible(list(img = image, mask = EBImage::Image(mask)))
  }
}


#' Convert an image to a data.frame
#'
#' Given an object image, converts it into a data frame where each row corresponds to the intensity values of each pixel in the image.
#' @param img An image object.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing three matrices (R, G, and B), and a data frame
#'   containing four columns: the name of the image in `image` and the R, G, B
#'   values.
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' dim(img)
#' mat <- image_to_mat(img)
#' dim(mat[[1]])
#' }
image_to_mat <- function(img,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE){
  check_ebi()
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      cli::cli_abort("All images must be of class {.code Image}.")
    }
    if(parallel == TRUE){
      # decide number of workers
      nworkers <- ifelse(is.null(workers),
                         trunc(parallel::detectCores() * 0.4),
                         workers)

      # start mirai daemons
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # CLI header + progress step
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Converting {.val {length(img)}} images to matrices"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Matrix conversion complete.",
          msg_failed = "Matrix conversion failed."
        )
      }

      # run image_to_mat in parallel
      res <- mirai::mirai_map(
        .x = img,
        .f = pliman::image_to_mat
      )[.progress]

    } else{
      res <- lapply(img, image_to_mat)
    }
    invisible(structure(res, class = "img_mat_list"))
  } else{
    mat <- cbind(expand.grid(Row = 1:dim(img)[1], Col = 1:dim(img)[2]))
    if(length(dim(img)) == 3){
      for (i in 1:dim(img)[3]) {
        mat <- cbind(mat, c(img[, , i]))
      }
      colnames(mat) = c("row", "col", paste0("B", 1:dim(img)[3]))
    } else{
      mat <- cbind(mat, c(img))
      colnames(mat) = c("row", "col", "B1")
    }
    invisible(mat)
  }
}


#' Create image palettes
#'
#' `image_palette()`  creates image palettes by applying the k-means algorithm
#' to the RGB values.
#' @inheritParams analyze_objects
#' @param img An image object.
#' @param npal The number of color palettes.
#' @param proportional Creates a joint palette with proportional size equal to
#'   the number of pixels in the image? Defaults to `TRUE`.
#' @param plot Plot the generated palette? Defaults to `TRUE`.
#' @param colorspace The color space to produce the clusters. Defaults to `rgb`.
#'   If `hsb`, the color space is first converted from RGB > HSB before k-means
#'   algorithm be applied.
#' @param remove_bg Remove background from the color palette? Defaults to
#'   `FALSE`.
#' @param index An image index used to remove the background, passed to
#'   [image_binary()].
#' @param parallel If TRUE processes the images asynchronously (in parallel) in
#'   separate R sessions running in the background on the same machine.
#' @param return_pal Return the color palette image? Defaults to `FALSE`.
#' @return `image_palette()` returns a list with two elements:
#' * `palette_list` A list with `npal` color palettes of class `Image`.
#' * `joint` An object of class `Image` with the color palettes
#' * `proportions` The proportion of the entire image corresponding to each color in the palette
#' * `rgbs` The average RGB value for each palette
#' @name palettes
#' @export
#' @importFrom stats na.omit
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'pal <- image_palette(img, npal = 5)
#'}
#'
#'

image_palette <- function (img,
                           pattern = NULL,
                           npal = 5,
                           proportional = TRUE,
                           colorspace = c("rgb", "hsb"),
                           remove_bg = FALSE,
                           index = "B",
                           plot = TRUE,
                           save_image = FALSE,
                           prefix = "proc_",
                           dir_original = NULL,
                           dir_processed = NULL,
                           return_pal = FALSE,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE) {
  check_ebi()
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

  help_pal <- function(img, npal, proportional, colorspace, plot, save_image, prefix){
    if(is.character(img)){
      all_files <- sapply(list.files(diretorio_original), file_name)
      imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
      name_ori <- file_name(imag)
      extens_ori <- file_extension(imag)
      img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
    } else{
      name_ori <- match.call()[[2]]
      extens_ori <- "jpg"
    }
    if (!colorspace[[1]] %in% c("rgb", "hsb")) {
      cli::cli_warn(c(
        "!" = "`colorspace` must be one of {.val 'rgb'} or {.val 'hsb'}.",
        " " = "Setting to {.val 'rgb'}."
      ))
      colorspace <- "rgb"
    }

    # remove BG if needed
    if(remove_bg){
      mask <- image_binary(img, index = "B-R", opening = 5, plot = FALSE, verbose = FALSE)[[1]]
      ID <- which(mask == FALSE)
      img@.Data[,,1][ID] <- NA
      img@.Data[,,2][ID] <- NA
      img@.Data[,,3][ID] <- NA
    }
    nc <- ncol(img)
    nr <- nrow(img)
    if(length(dim(img)) == 1){
      if(colorspace[[1]] == "hsb"){
        cli::cli_abort("HSB can only be computed with an 3 layers array (RGB).")
      } else{
        imb <- data.frame(B1 = rgb_to_hsb(img)[,3])
      }
    } else if(length(dim(img)) == 3){
      if(colorspace[[1]] == "hsb"){
        imb <- image_to_mat(img)[, -c(1, 2)]
        rownames(imb) <- paste0("r", 1:nrow(imb))
        hsb <- rgb_to_hsb(img)
        rownames(hsb) <- paste0("r", 1:nrow(hsb))
      } else{
        imb <- image_to_mat(img)[, -c(1, 2)]
        rownames(imb) <- paste0("r", 1:nrow(imb))
      }
    }
    if(any(is.na(imb[, 1]))){
      if(colorspace[[1]] == "rgb"){
        set.seed(10)
        km <- suppressWarnings(stats::kmeans(na.omit(imb), npal))
      } else{
        set.seed(10)
        km <- suppressWarnings(stats::kmeans(na.omit(hsb), npal))
      }
      imb <- cbind(imb, 'cluster'=NA)
      imb[names(km$cluster), "cluster"] <- km$cluster
    } else{
      if(colorspace[[1]] == "rgb"){
        set.seed(10)
        km <- suppressWarnings(stats::kmeans(imb, npal))
        imb$cluster <- km$cluster
      } else{
        set.seed(10)
        km <- suppressWarnings(stats::kmeans(hsb, npal))
        imb$cluster <- km$cluster
      }
    }
    if(colorspace[[1]] == "hsb"){
      props <-
        imb |>
        dplyr::bind_cols(hsb) |>
        dplyr::relocate(cluster, .before = 1) |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(
          n = dplyr::n(),
          R = mean(B1, na.rm = TRUE),
          G = mean(B2, na.rm = TRUE),
          B = mean(B3, na.rm = TRUE),
          h = mean(h, na.rm = TRUE),
          s = mean(s, na.rm = TRUE),
          b = mean(b, na.rm = TRUE)
        )
    } else{
      props <-
        imb |>
        dplyr::filter(!is.na(cluster)) |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(
          n = dplyr::n(),
          R = mean(B1, na.rm = TRUE),
          G = mean(B2, na.rm = TRUE),
          B = mean(B3, na.rm = TRUE)
        )
    }
    props <-
      props |>
      dplyr::mutate(prop = n / sum(n), .after = n) |>
      dplyr::arrange(prop) |>
      dplyr::mutate(cluster = paste0("c", 1:npal),
                    .before = 1) |>
      dplyr::ungroup()
    if(plot){

      pal_list <- list()
      pal_rgb <- list()
      for(i in 1:nrow(props)){
        R <- matrix(rep(props[[i, 4]], 10000), 100, 100)
        G <- matrix(rep(props[[i, 5]], 10000), 100, 100)
        B <- matrix(rep(props[[i, 6]], 10000), 100, 100)
        pal_list[[paste0("pal_", i)]] <- EBImage::rgbImage(R, G, B)
        pal_rgb[[paste0("pal_", i)]] <- c(R = R[1], G = G[1], B = B[1])
      }


      rownames(props) <- NULL
      if (proportional == FALSE) {
        n <- nrow(props)
        ARR <- array(NA, dim = c(100, 66 * n, 3))
        c = 1
        f = 66
        for (i in 1:n) {
          ARR[1:100, c:f, 1] <- props[[i, 4]]
          ARR[1:100, c:f, 2] <- props[[i, 5]]
          ARR[1:100, c:f, 3] <- props[[i, 6]]
          c = f + 1
          f = f + 66
        }
      }
      if (proportional == TRUE) {
        n <- nrow(props)
        ARR <- array(NA, dim = c(100, 66 * n, 3))
        nn <- round(66 * n * props$prop, 0)
        a <- 1
        b <- nn[1]
        nn <- c(nn, 0)
        for (i in 1:n) {
          ARR[1:100, a:b, 1] <- props[[i, 4]]
          ARR[1:100, a:b, 2] <- props[[i, 5]]
          ARR[1:100, a:b, 3] <- props[[i, 6]]
          a <- b + 1
          b <- b + nn[i + 1]
          if (b > (66 * n)) {
            b <- 66 * n
          }
        }
      }
      im2 <- EBImage::as.Image(ARR)
      EBImage::colorMode(im2) <- 2
      im2 <- image_resize(im2, height = ncol(img), width = nc * 0.1)
      im2@.Data[1:nrow(im2), 1:1, ] <- 0
      im2@.Data[1:1, 1:ncol(im2), ] <- 0
      im2@.Data[1:nrow(im2), ncol(im2):(ncol(im2)-1), ] <- 0
      im2@.Data[nrow(im2):(nrow(im2)-1), 1:ncol(im2), ] <- 0
      im2 <- EBImage::abind(img, im2, along = 1)
      if (plot == TRUE) {
        plot(im2)
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada, recursive = TRUE)
        }
        jpeg(paste0(diretorio_processada, "/",
                    prefix,
                    name_ori, ".",
                    extens_ori),
             width = dim(im2@.Data)[1],
             height = dim(im2@.Data)[2])
        plot(im2)
        dev.off()
      }
    } else{
      im2 <- NULL
      pal_list <- NULL
    }
    if(!return_pal){
      im2 <- NULL
      pal_list <- NULL
    }
    # gc()
    return(list(palette_list = pal_list,
                joint = im2,
                proportions = props))
  }
  if(missing(pattern)){
    if(verbose){
      cli::cli_progress_step(
        msg = "Processing a single image. Please, wait.",
        msg_done = "Image {.emph Successfully} analyzed!",
        msg_failed = "Oops, something went wrong."
      )
    }
    help_pal(img, npal, proportional, colorspace, plot, save_image, prefix)
  } else {
    # anchor purely-numeric patterns
    if (pattern %in% as.character(0:9)) {
      pattern <- "^[0-9].*$"
    }

    # list files
    plants      <- list.files(pattern = pattern, diretorio_original)
    extensions  <- tolower(vapply(plants, tools::file_ext,   ""))
    names_plant <-         vapply(plants, tools::file_path_sans_ext, "")

    # abort if no matches
    if (length(plants) == 0) {
      cli::cli_abort(c(
        "x" = "Pattern {.val {pattern}} not found in {.path {diretorio_original}}.",
        "i" = "Check working dir: {.path {getwd()}}"
      ))
    }

    # abort on bad extensions
    bad_ext <- setdiff(extensions, c("png","jpeg","jpg","tiff"))
    if (length(bad_ext)) {
      cli::cli_abort(c(
        "x" = "Unsupported extension{?s}: {.val {unique(bad_ext)}} found.",
        "i" = "Allowed: {.val png}, {.val jpeg}, {.val jpg}, {.val tiff}."
      ))
    }


    if (parallel) {
      # parallel setup
      nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.3) else workers
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = "Parallel palette extraction of {.val {length(names_plant)}} images",
          right = "Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}"
        )
        cli::cli_progress_step(
          msg        = "Dispatching batches...",
          msg_done   = "All batches complete!",
          msg_failed = "Batch failed."
        )
      }

      # run help_pal in parallel with mirai
      results <- mirai::mirai_map(
        .x = names_plant,
        .f = function(img_path) {
          help_pal(
            img          = img_path,
            npal          = npal,
            proportional  = proportional,
            colorspace    = colorspace,
            plot          = plot,
            save_image    = save_image,
            prefix        = prefix
          )
        }
      )[.progress]

    } else {
      # sequential processing
      if (verbose) {
        cli::cli_rule(
          left  = "Sequential palette extraction of {.val {length(names_plant)}} images",
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_bar(
          format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | Current: {.val {cli::pb_status}}",
          total  = length(names_plant),
          clear  = TRUE
        )
      }

      results <- vector("list", length(names_plant))
      for (i in seq_along(names_plant)) {
        if (verbose) cli::cli_progress_update(status = names_plant[i])
        results[[i]] <- help_pal(
          img          = names_plant[i],
          npal, proportional, colorspace, plot,
          save_image, prefix
        )
      }


    }

    # assemble outputs
    names(results) <- names_plant
    proportions  <- do.call(rbind, lapply(seq_along(results), function(i) {
      results[[i]]$proportions |>
        dplyr::mutate(img = names_plant[i], .before = 1)
    }))
    palette_list <- lapply(results, `[[`, "palette_list")
    joint        <- lapply(results, `[[`, "joint")
    names(joint) <- names_plant

    if (verbose) {
      cli::cli_progress_done()
      cli::cli_rule(
        left  = cli::col_green("All {.val {length(names_plant)}} images processed"),
        right = cli::col_blue("Finished on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
    }

    return(list(
      proportions  = proportions,
      palette_list = palette_list,
      joint        = joint
    ))
  }

}








#' Expands an image
#'
#' Expands an image towards the left, top, right, or bottom by sampling pixels
#' from the image edge. Users can choose how many pixels (rows or columns) are
#' sampled and how many pixels the expansion will have.
#'
#' @param img An `Image` object.
#' @param left,top,right,bottom The number of pixels to expand in the left, top,
#'   right, and bottom directions, respectively.
#' @param edge The number of pixels to expand in all directions. This can be
#'   used to avoid calling all the above arguments
#' @param sample_left,sample_top,sample_right,sample_bottom The number of pixels
#'   to sample from each side. Defaults to 20.
#' @param random Randomly sampling of the edge's pixels? Defaults to `FALSE`.
#' @param filter Apply a median filter in the sampled pixels? Defaults to
#'   `FALSE`.
#' @param plot Plots the extended image? defaults to `FALSE`.
#'
#' @return An `Image` object
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' image_expand(img, left = 200)
#' image_expand(img, right = 150, bottom = 250, filter = 5)
#' }
#'
image_expand <- function(img,
                         left = NULL,
                         top = NULL,
                         right = NULL,
                         bottom = NULL,
                         edge = NULL,
                         sample_left = 10,
                         sample_top = 10,
                         sample_right = 10,
                         sample_bottom = 10,
                         random = FALSE,
                         filter = NULL,
                         plot = TRUE){
  check_ebi()
  if(!is.null(edge)){
    left <- edge
    top <- edge
    right <- edge
    bottom <- edge
  }
  if (sample_left < 2) {
    cli::cli_warn("{.arg sample_left} must be > {.val 1}. Setting to {.val 2}.")
    sample_left <- 2
  }
  if (sample_top < 2) {
    cli::cli_warn("{.arg sample_top} must be > {.val 1}. Setting to {.val 2}.")
    sample_top <- 2
  }
  if (sample_right < 2) {
    cli::cli_warn("{.arg sample_right} must be > {.val 1}. Setting to {.val 2}.")
    sample_right <- 2
  }
  if (sample_bottom < 2) {
    cli::cli_warn("{.arg sample_bottom} must be > {.val 1}. Setting to {.val 2}.")
    sample_bottom <- 2
  }

  if(!is.null(left)){
    left_img <- img@.Data[1:sample_left,,] |> EBImage::Image(colormode = "Color")
    left_img <- EBImage::resize(left_img, w = left, h = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(left_img)
      for (i in 1:nc[1]) {
        left_img@.Data[i,,] <- left_img@.Data[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      left_img <- EBImage::medianFilter(left_img, size = filter)
    }
    img <- EBImage::abind(left_img, img, along = 1)
  }
  if(!is.null(top)){
    top_img <- img@.Data[,1:sample_top,] |> EBImage::Image(colormode = "Color")
    top_img <- EBImage::resize(top_img, w = dim(img)[1], h = top)
    if(isTRUE(random)){
      nc <- dim(top_img)
      for (i in 1:nc[2]) {
        top_img@.Data[,i,] <- top_img@.Data[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      top_img <- EBImage::medianFilter(top_img, size = filter)
    }
    img <- EBImage::abind(top_img, img, along = 2)
  }
  if(!is.null(right)){
    dimx <- dim(img)[1]
    right_img <- img@.Data[(dimx-sample_right):dimx,,] |> EBImage::Image(colormode = "Color")
    right_img <- EBImage::resize(right_img, w = right, h = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(right_img)
      for (i in 1:nc[1]) {
        right_img@.Data[i,,] <- right_img@.Data[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      right_img <- EBImage::medianFilter(right_img, size = filter)
    }
    img <- EBImage::abind(img, right_img, along = 1)
  }
  if(!is.null(bottom)){
    dimy <- dim(img)[2]
    bot_img <- img@.Data[,(dimy-sample_bottom):dimy,] |> EBImage::Image(colormode = "Color")
    bot_img <- EBImage::resize(bot_img, w = dim(img)[1], h = bottom)
    if(isTRUE(random)){
      nc <- dim(bot_img)
      for (i in 1:nc[2]) {
        bot_img@.Data[,i,] <- bot_img@.Data[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      bot_img <- EBImage::medianFilter(bot_img, size = filter)
    }
    img <- EBImage::abind(img, bot_img, along = 2)
  }
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}


#' Squares an image
#'
#' Converts a rectangular image into a square image by expanding the
#' rows/columns using [image_expand()].
#'
#' @inheritParams image_expand
#'
#' @return The modified `Image` object.
#' @param ... Further arguments passed on to [image_expand()].
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' dim(img)
#' square <- image_square(img)
#' dim(square)
#' }
image_square <- function(img, plot = TRUE, ...){
  len <- dim(img)
  n <- max(len[1], len[2])
  if (len[1] > len[2]) {
    ni1 <- ceiling((n - len[2])/2)
    if((ni1*2 + len[2]) != n){
      ni2 <- ni1 - 1
    } else{
      ni2 <- ni1
    }
    img <- image_expand(img, bottom = ni1, top = ni2, plot = FALSE, ...)
  }
  if (len[2] > len[1]) {
    ni1 <- ceiling((n - len[1])/2)
    if((ni1*2 + len[1]) != n){
      ni2 <- ni1 - 1
    } else{
      ni2 <- ni1
    }
    img <- image_expand(img, left = ni1, right = ni2, plot = FALSE, ...)
  }
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}




#' Utilities for image resolution
#'
#' Provides useful conversions between size (cm), number of pixels (px) and
#' dots per inch (dpi).
#' * [dpi_to_cm()] converts a known dpi value to centimeters.
#' * [cm_to_dpi()] converts a known centimeter values to dpi.
#' * [pixels_to_cm()] converts the number of pixels to centimeters, given a
#' known resolution (dpi).
#' * [cm_to_pixels()] converts a distance (cm) to number of pixels, given a
#' known resolution (dpi).
#' * [distance()] Computes the distance between two points in an image based on
#' the Pythagorean theorem.
#' * [dpi()] An interactive function to compute the image resolution given a
#' known distance informed by the user. See more information in the **Details**
#' section.
#' * [npixels()] returns the number of pixels of an image.
#' @details [dpi()] only run in an interactive section. To compute the image
#'   resolution (dpi) the user must use the left button mouse to create a line
#'   of known distance. This can be done, for example, using a template with
#'   known distance in the image (e.g., `la_leaves.jpg`).
#'
#' @name utils_dpi
#' @inheritParams image_view
#' @param img An image object.
#' @param dpi The image resolution in dots per inch.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param px The number of pixels.
#' @param cm The size in centimeters.
#' @return
#' * [dpi_to_cm()], [cm_to_dpi()], [pixels_to_cm()], and [cm_to_pixels()] return
#' a numeric value or a vector of numeric values if the input data is a vector.
#' * [dpi()] returns the computed dpi (dots per inch) given the known distance
#' informed in the plot.
#' @export
#' @importFrom grDevices rgb2hsv convertColor
#' @importFrom graphics locator
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' library(pliman)
#' # Convert  dots per inch to centimeter
#' dpi_to_cm(c(1, 2, 3))
#'
#' # Convert centimeters to dots per inch
#' cm_to_dpi(c(1, 2, 3))
#'
#' # Convert centimeters to number of pixels with resolution of 96 dpi.
#' cm_to_pixels(c(1, 2, 3), 96)
#'
#' # Convert number of pixels to cm with resolution of 96 dpi.
#' pixels_to_cm(c(1, 2, 3), 96)
#'
#' if(isTRUE(interactive())){
#' #### compute the dpi (dots per inch) resolution ####
#' # only works in an interactive section
#' # objects_300dpi.jpg has a known resolution of 300 dpi
#' img <- image_pliman("objects_300dpi.jpg")
#' # Higher square: 10 x 10 cm
#' # 1) Run the function dpi()
#' # 2) Use the left mouse button to create a line in the higher square
#' # 3) Declare a known distance (10 cm)
#' # 4) See the computed dpi
#' dpi(img)
#'
#'
#' img2 <- image_pliman("la_leaves.jpg")
#' # square leaf sample (2 x 2 cm)
#' dpi(img2)
#' }
dpi_to_cm <- function(dpi){
  2.54 / dpi
}
#' @name utils_dpi
#' @export
cm_to_dpi <- function(cm){
  cm / 2.54
}
#' @name utils_dpi
#' @export
pixels_to_cm <- function(px, dpi){
  px * (2.54 / dpi)
}
#' @name utils_dpi
#' @export
cm_to_pixels <- function(cm, dpi){
  cm / (2.54 / dpi)
}
#' @name utils_dpi
#' @export
npixels <- function(img){
  if(!inherits(img, "Image")){
    cli::cli_abort("Image must be of class 'Image'.")
  }
  dim <- dim(img)
  dim[[1]] * dim[[2]]
}
#' @name utils_dpi
#' @export
dpi <- function(img,
                viewer = get_pliman_viewer(),
                downsample = NULL,
                max_pixels = 1000000){
  if(isTRUE(interactive())){
    pix <- distance(img, viewer = viewer, downsample = downsample, max_pixels = max_pixels)
    known <- as.numeric(readline("known distance (cm): "))
    pix / (known / 2.54)
  }
}

#' @name utils_dpi
#' @export
distance <- function(img,
                     viewer = get_pliman_viewer(),
                     downsample = NULL,
                     max_pixels = 1000000){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(isTRUE(interactive())){
    if(vieweropt == "base"){
      plot(img)
      cli::cli_alert_info("Use the first mouse button to create a line in the plot.")
      coords <- locator(type = "l",
                        n = 2,
                        lwd = 2,
                        col = "red")
      pix <- sqrt((coords$x[1] - coords$x[2])^2 + (coords$y[1] - coords$y[2])^2)
    } else{
      coords2 <- mv_two_points(img, downsample = downsample, max_pixels = max_pixels)
      pix <- sqrt((coords2$x1 - coords2$x2)^2 + (coords2$y1 - coords2$y2)^2)
    }
    invisible(pix)
  }
}


#' Convert between colour spaces
#' @description
#'  * `rgb_to_srgb()` Transforms colors from RGB space (red/green/blue) to
#'  Standard Red Green Blue (sRGB), using a gamma correction of 2.2. The
#'  function performs the conversion by applying a gamma correction to the input
#'  RGB values (raising them to the power of 2.2) and then transforming them
#'  using a specific transformation matrix. The result is clamped to the range
#'  0-1 to ensure valid sRGB values.
#'
#'
#' * `rgb_to_hsb()` Transforms colors from RGB space (red/green/blue) to HSB
#' space (hue/saturation/brightness). The HSB values are calculated as follows
#' (see https://www.rapidtables.com/convert/color/rgb-to-hsv.html for more
#' details).
#'    - Hue: The hue is determined based on the maximum value among R, G, and B,
#'    and it ranges from 0 to 360 degrees.
#'    - Saturation: Saturation is calculated as the difference between the maximum
#'    and minimum channel values, expressed as a percentage.
#'    - Brightness: Brightness is equal to the maximum channel value, expressed as
#'    a percentage.
#'
#'
#'  * `rgb_to_lab()` Transforms colors from RGB space (red/green/blue) to CIE-LAB
#'  space, using the sRGB values. See [grDevices::convertColor()] for more
#'  details.
#'
#'
#' @param object An `Image` object, an object computed with `analyze_objects()`
#'   with a valid `object_index` argument, or a `data.frame/matrix`. For the
#'   last, a three-column data (R, G, and B, respectively) is required.
#'
#' @references
#'See [the detailed formulas here](https://www.example.com)
#'
#' @export
#' @name utils_colorspace
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A data frame with the columns of the converted color space
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' rgb_to_lab(img)
#'
#' # analyze the object and convert the pixels
#' anal <- analyze_objects(img, object_index = "B", pixel_level_index = TRUE)
#' rgb_to_lab(anal)
#' }
rgb_to_hsb <- function(object){
  if (any(class(object) %in%  c("data.frame", "matrix"))){
    hsb <-
      rgb_to_hsb_help(r = object[,1],
                      g = object[,2],
                      b = object[,3])
    colnames(hsb) <- c("h", "s", "b")
  }
  if (any(class(object)  %in% c("anal_obj", "anal_obj_ls"))){
    if(!is.null(object$object_rgb)){
      tmp <- object$object_rgb
      if ("img" %in% colnames(tmp)){
        hsb <-
          rgb_to_hsb_help(r = c(tmp[,3]),
                          g = c(tmp[,4]),
                          b = c(tmp[,5]))
        hsb <- data.frame(cbind(tmp[,1:2], hsb))
        colnames(hsb)[1:2] <- c("img", "id")
        colnames(hsb)[3:5] <- c("h", "s", "b")
      }
      hsb <-
        rgb_to_hsb_help(r = c(tmp[,2]),
                        g = c(tmp[,3]),
                        b = c(tmp[,4]))
      hsb <- data.frame(cbind(tmp[,1], hsb))
      colnames(hsb)[1] <- "id"
      colnames(hsb)[2:4] <- c("h", "s", "b")
    } else{
      cli::cli_abort(c(
        "!" = "Cannot obtain the RGB for each object since the {.arg object_index} argument was not used.",
        "i" = "Have you accidentally missed the argument {.arg pixel_level_index} = TRUE?"
      ))

    }
  }
  if (any(class(object) == "Image")){
    hsb <-
      rgb_to_hsb_help(r = c(object[,,1]),
                      g = c(object[,,2]),
                      b = c(object[,,3]))
    colnames(hsb) <- c("h", "s", "b")
  }
  invisible(data.frame(hsb))
}

#' @export
#' @name utils_colorspace
rgb_to_srgb <- function(object){
  if (any(class(object) %in%  c("data.frame", "matrix"))){
    srgb <- rgb_to_srgb_help(object[, 1:3])
    colnames(srgb) <- c("sR", "sG", "sB")
  }
  if (any(class(object)  %in% c("anal_obj", "anal_obj_ls"))){
    if(!is.null(object$object_rgb)){
      tmp <- object$object_rgb
      if ("img" %in% colnames(tmp)){
        srgb <- rgb_to_srgb_help(as.matrix(tmp[, 3:5]))
        srgb <- data.frame(cbind(tmp[,1:2], srgb))
        colnames(srgb)[1:2] <- c("img", "id")
        colnames(srgb)[3:5] <- c("sR", "sG", "sB")
      } else{
        srgb <- rgb_to_srgb_help(as.matrix(tmp[,2:4]))
        srgb <- data.frame(cbind(tmp[,1], srgb))
        colnames(srgb)[1] <- "id"
        colnames(srgb)[2:4] <- c("sR", "sG", "sB")
      }
    } else{
      cli::cli_abort(c(
        "!" = "Cannot obtain the RGB for each object since the {.arg object_index} argument was not used.",
        "i" = "Have you accidentally missed the argument {.arg pixel_level_index} = TRUE?"
      ))

    }
  }
  if (any(class(object) == "Image")){
    srgb <- rgb_to_srgb_help(cbind(c(object[,,1]), c(object[,,2]), c(object[,,3])))
    colnames(srgb) <- c("sR", "sG", "sB")
  }
  invisible(data.frame(srgb))
}


#' @export
#' @name utils_colorspace
rgb_to_lab <- function(object){
  object <- rgb_to_srgb(object)
  srgb <- data.frame(r = object[, 1],
                     g = object[, 2],
                     b = object[, 3])
  lab <- convertColor(srgb, from = "sRGB", to = "Lab")
  invisible(lab)
}




# Faster alternatives (makes only the needed)
help_segment <- function(img,
                         index = NULL,
                         r = 1,
                         g = 2,
                         b = 3,
                         re = 4,
                         nir = 5,
                         threshold = c("Otsu", "adaptive"),
                         k = 0.1,
                         windowsize = NULL,
                         col_background = NULL,
                         has_white_bg = FALSE,
                         fill_hull = FALSE,
                         opening = FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         dilate = FALSE,
                         erode = FALSE,
                         invert = FALSE){
  img2 <- help_binary(img,
                      index = index,
                      r = r,
                      g = g,
                      b = b,
                      re = re,
                      nir = nir,
                      threshold = threshold,
                      k = k,
                      windowsize = windowsize,
                      has_white_bg = has_white_bg,
                      resize = FALSE,
                      fill_hull = fill_hull,
                      opening = opening,
                      closing = closing,
                      filter = filter,
                      dilate = dilate,
                      erode = erode,
                      invert = invert)
  ID <- which(img2@.Data == FALSE)
  if(dim(img)[3] == 3){
    img@.Data[,,r][ID] <- 1
    img@.Data[,,g][ID] <- 1
    img@.Data[,,b][ID] <- 1
  } else if(dim(img)[3] == 4){
    img@.Data[,,r][ID] <- 1
    img@.Data[,,g][ID] <- 1
    img@.Data[,,b][ID] <- 1
    img@.Data[,,re][ID] <- 1
  } else{
    img@.Data[,,r][ID] <- 1
    img@.Data[,,g][ID] <- 1
    img@.Data[,,b][ID] <- 1
    img@.Data[,,re][ID] <- 1
    img@.Data[,,nir][ID] <- 1
  }
  invisible(img)
}



help_binary <- function(img,
                        index = NULL,
                        r = 1,
                        g = 2,
                        b = 3,
                        re = 4,
                        nir = 5,
                        threshold = c("Otsu", "adaptive"),
                        k = 0.15,
                        windowsize = NULL,
                        has_white_bg = FALSE,
                        resize = FALSE,
                        fill_hull = FALSE,
                        erode = FALSE,
                        dilate = FALSE,
                        opening = FALSE,
                        closing = FALSE,
                        filter = FALSE,
                        invert = FALSE){
  threshold <- threshold[[1]]

  bin_img <- function(imgs,
                      invert,
                      fill_hull,
                      threshold,
                      erode,
                      dilate,
                      opening,
                      closing,
                      filter){
    # adapted from imagerExtra  https://bit.ly/3Wp4pwv
    if(threshold == "adaptive"){
      if(is.null(windowsize)){
        windowsize <- min(dim(imgs)) / 3
        if(windowsize %% 2 == 0){
          windowsize <- as.integer(windowsize + 1)
        }
      }
      if (windowsize <= 2) {
        cli::cli_abort("{.arg windowsize} must be greater than or equal to {.val 3}.")
      }

      if (windowsize %% 2 == 0) {
        cli::cli_warn(
          "{.arg windowsize} is even ({.val {windowsize}}). It will be treated as {.val {windowsize + 1}}."
        )
        windowsize <- as.integer(windowsize + 1)
      }

      if (windowsize >= dim(imgs)[[1]] || windowsize >= dim(imgs)[[2]]) {
        cli::cli_warn(
          "{.arg windowsize} is too large. Setting to {.code min(dim(img)) / 3}."
        )
        windowsize <- min(dim(imgs)) / 3
      }

      if (k > 1) {
        cli::cli_abort("{.arg k} must be in range {.val [0, 1]}.")
      }
      imgs <- EBImage::thresh(imgs, w=windowsize, h=windowsize, offset=k)
      # imgs <- EBImage::Image(threshold_adaptive(as.matrix(imgs), k, windowsize, 0.5))
    }
    if(threshold != "adaptive"){
      if(threshold == "Otsu"){
        if(any(is.infinite(imgs)) | any(is.na(imgs))){
          threshold <- help_otsu(imgs@.Data[!is.infinite(imgs@.Data) & !is.na(imgs@.Data)])
        } else{
          threshold <- help_otsu(imgs@.Data)
        }
      } else{
        if(is.numeric(threshold)){
          threshold <- threshold
        } else{
          pixels <- terra::rast(EBImage::transpose(imgs)@.Data)
          terra::plot(pixels, col = custom_palette(n = 100),  axes = FALSE, asp = NA)
          threshold <- readline("Selected threshold: ")
        }
      }
      imgs <- EBImage::Image(imgs < threshold)
    }

    if(invert == TRUE){
      imgs <- 1 - imgs
    }

    imgs[which(is.na(imgs))] <- FALSE
    if(is.numeric(erode) & erode > 0){
      imgs <- image_erode(imgs, size = erode)
    }
    if(is.numeric(dilate) & dilate > 0){
      imgs <- image_dilate(imgs, size = dilate)
    }
    if(is.numeric(opening) & opening > 0){
      imgs <- image_opening(imgs, size = opening)
    }
    if(is.numeric(closing) & closing > 0){
      imgs <- image_closing(imgs, size = closing)
    }
    if(is.numeric(filter) & filter > 1){
      imgs <- EBImage::medianFilter(imgs, filter)
    }
    if(isTRUE(fill_hull)){
      imgs <- EBImage::fillHull(imgs)
    }
    invisible(imgs)
  }

  gray_img <- help_imageindex(img, index, r, g, b, re, nir, resize, has_white_bg)
  bin_img <- bin_img(gray_img,
                     invert,
                     fill_hull,
                     threshold,
                     erode,
                     dilate,
                     opening,
                     closing,
                     filter)
  invisible(bin_img)
}



help_imageindex <- function(img,
                            index = NULL,
                            r = 1,
                            g = 2,
                            b = 3,
                            re = 4,
                            nir = 5,
                            resize = FALSE,
                            has_white_bg = FALSE){
  if(resize != FALSE){
    img <- image_resize(img, resize)
  }
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")

  nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
  hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])

  R <- try(img@.Data[,,r], TRUE)
  G <- try(img@.Data[,,g], TRUE)
  B <- try(img@.Data[,,b], TRUE)
  RE <- try(img@.Data[,,re], TRUE)
  NIR <- try(img@.Data[,,nir], TRUE)

  if(any(index %in% hsb_ind)){
    hsb <- rgb_to_hsb(data.frame(R = c(R), G = c(G), B = c(B)))
    h <- matrix(hsb$h, nrow = nrow(img), ncol = ncol(img))
    s <- matrix(hsb$s, nrow = nrow(img), ncol = ncol(img))
    b <- matrix(hsb$b, nrow = nrow(img), ncol = ncol(img))
  }

  if(any(index %in% nir_ind)){
    test_multi <- any(sapply(list(RE, NIR), class) == "try-error")
    if(isTRUE(test_multi)){
      cli::cli_abort("Near-Infrared and RedeEdge bands are not available in the provided image.")
    }
  }
  if(isTRUE(has_white_bg)){
    R[which(R == 1 & G == 1 & B == 1)] <- NA
    G[which(R == 1 & G == 1 & B == 1)] <- NA
    B[which(R == 1 & G == 1 & B == 1)] <- NA
  }

  if(index %in% ind$Index){
    img_gray <- EBImage::Image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index]))))
  } else{
    img_gray <- EBImage::Image(eval(parse(text = as.character(index))))
  }
  invisible(img_gray)
}


#' Create an `Image` object
#'
#' This function is a simple wrapper around [EBImage::Image()].
#'
#' @param data A vector or array containing the pixel intensities of an image.
#'   If missing, the default 1x1 zero-filled array is used.
#' @param ... Additional arguments passed to [EBImage::Image()].
#' @return An `Image` object.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' img <-
#' as_image(rnorm(150 * 150 * 3),
#'          dim = c(150, 150, 3),
#'          colormode = 'Color')
#' plot(img)
#' }
as_image <- function(data, ...){
  check_ebi()
  EBImage::Image(data, ...)
}

#' Prepare images to analyze_objects_shp()
#'
#' It is a simple wrapper around [image_align()] and [image_crop()]. In this case, only the option `viewer = "base"` is used. To use `viewer = "mapview"`, please, use such functions separately.
#'
#' @param img A `Image` object
#' @inheritParams image_align
#'
#' @return An aligned and cropped `Image` object.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' img <- image_pliman("flax_leaves.jpg")
#' prepare_to_shp(img)
#' }
prepare_to_shp <- function(img,
                           align = "vertical"){
  check_ebi()
  aligned <- image_align(img, viewer = "base")
  cropped <- image_crop(aligned, viewer = "base", plot = TRUE)
  invisible(cropped)
}

#' Add Alpha Layer to an RGB Image
#'
#' This function adds an alpha (transparency) layer to an RGB image using the EBImage package.
#' The alpha layer can be specified as a single numeric value for uniform transparency
#' or as a matrix/array matching the dimensions of the image for varying transparency.
#'
#' @param img An RGB image of class `Image` from the EBImage package. The image must be in RGB format (color mode 2).
#' @param mask A numeric value or matrix/array specifying the alpha layer:
#'     * If `mask` is a single numeric value, it sets a uniform transparency level (0 for fully transparent, 1 for fully opaque).
#'     * If `mask` is a matrix or array, it must have the same dimensions as the image channels, allowing for varying transparency.
#'
#' @return An `Image` object with an added alpha layer, maintaining the RGBA format.
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' # Load the EBImage package
#' library(pliman)
#'
#' # Load a sample RGB image
#' img <- image_pliman("soybean_touch.jpg")
#'
#' # 50% transparency
#' image_alpha(img, 0.5) |> plot()
#'
#' # transparent background
#' mask <- image_binary(img, "NB")[[1]]
#' img_tb <- image_alpha(img, mask)
#' plot(img_tb)
#'
#' }
#'
#' @export
#'
image_alpha <- function(img, mask) {
  check_ebi()
  # Check if the input image is an RGB image
  if (EBImage::colorMode(img) != 2) {
    cli::cli_abort("Input image must be in RGB format.")
  }
  r <- img[,,1]
  g <- img[,,2]
  b <- img[,,3]
  EBImage::colorMode(r) <- "Grayscale"
  EBImage::colorMode(g) <- "Grayscale"
  EBImage::colorMode(b) <- "Grayscale"
  if (is.numeric(mask) && length(mask) == 1) {
    alpha_layer <- array(mask, dim = dim(r))
  }
  # If mask is a matrix or array, ensure it matches the image dimensions
  else if (all(dim(mask) == dim(r))) {
    alpha_layer <- array(as.numeric(mask), dim = dim(img)[1:2])
  } else {
    cli::cli_abort("Mask must be either a single numeric value or a matrix with the same dimensions as the image channels.")
  }
  img_with_alpha <- EBImage::combine(r, g, b, alpha_layer)
  EBImage::colorMode(img_with_alpha) <- "Color"
  return(img_with_alpha)
}

#' Label Connected Components in a Binary Image
#'
#' This function labels connected components in a binary image while allowing
#' for a specified maximum gap between pixels to still be considered part of
#' the same object.
#'
#' @param img A binary image matrix where `1` represents foreground pixels and
#'   `0` represents background pixels. This should be compatible with the
#'   EBImage package.
#' @param max_gap An integer specifying the maximum allowable gap (in pixels)
#'   between connected components to be considered as part of the same object.
#'   Default is `1`.
#'
#' @return An object of class `Image` (from the EBImage package), where each
#'   connected component is assigned a unique integer label.
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- matrix(c(
#'   1, 1, 0, 0, 0, 1, 1, 1, 0,
#'   0, 0, 0, 0, 0, 1, 0, 0, 0,
#'   1, 1, 0, 0, 1, 1, 1, 0, 0,
#'   0, 0, 0, 0, 0, 0, 0, 0, 1
#' ), nrow = 4, byrow = TRUE)
#'
#' image_label(img, max_gap = 1)
#' image_label(img, max_gap = 2)
#' image_label(img, max_gap = 3)
#' }
image_label <- function(img, max_gap = 0){
  help_label(img, max_gap = max_gap) |> EBImage::as.Image()
}



#' Smooth Contour Line Detection
#'
#'
#' @param img An `Image` object.
#' @param index A character string with the index to be used. Defaults to `"GRAY"`.
#' @param Q numeric value with the pixel quantization step
#' @return A list with the contour lines.
#' @importFrom utils tail
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' conts <- image_contour_line(img, index = "B")
#' plot(img)
#' plot_contour(conts, col = "black")
#' }
#'
image_contour_line <- function(img, index = "GRAY", Q = 2.0){
  ind <- image_index(img, index, plot = FALSE)[[1]]
  mat <- ind@.Data * 255
  contourlines <- utils_contours(mat,
                                 X = nrow(mat),
                                 Y = ncol(mat),
                                 Q = Q)
  names(contourlines) <- c("x", "y", "curvelimits", "curves", "contourpoints")
  contourlines$curvelimits <- contourlines$curvelimits+1L
  from <- contourlines$curvelimits
  to <- c(tail(contourlines$curvelimits, contourlines$curves-1L)-1L, contourlines$contourpoints)
  curve <- unlist(mapply(seq_along(from), from, to, FUN=function(contourid, from, to) rep(contourid, to-from+1L), SIMPLIFY = FALSE))
  contourlines$data <- data.frame(x = contourlines$x, y = contourlines$y, curve = curve)
  res <- split(contourlines$data, contourlines$data$curve)
  res <-
    lapply(res, function(x){
      as.matrix(x[, 1:2])
    })
  return(res)
}


#' @title Canny Edge Detector
#' @description Canny Edge Detector for Images. Adapted from \url{https://github.com/bnosac/image/tree/master/image.CannyEdges}.
#' @param img An `Image` object.
#' @param index A character string with the index to be used. Defaults to `"GRAY"`.
#' @param s sigma, the Gaussian filter variance. Defaults to 5.
#' @param low_thr lower threshold value of the algorithm. Defaults to 10.
#' @param high_thr upper threshold value of the algorithm. Defaults to 20
#' @return a list with an `Image` object with values 0 or 255, and the number of
#'   pixels which have value 255 (pixels_nonzero).
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' conts <- image_canny_edge(img, index = "B")
#' par(mfrow = c(1, 2))
#' plot(img)
#' plot(conts$edges)
#' par(mfrow = c(1, 1))
#' }
image_canny_edge <- function(img,
                             index = "GRAY",
                             s = 5,
                             low_thr = 10,
                             high_thr = 20) {
  ind <- image_index(img, index, plot = FALSE)[[1]]
  mat <- ind@.Data * 255
  res <- canny_edge_detector(mat, nrow(mat), ncol(mat), s, low_thr, high_thr, TRUE)
  res$edges <- EBImage::as.Image(res$edges)
  return(res[1:2])
}

#' @title Line Segment Detection in an Image
#' @description Detects line segments in a digital image using the Line Segment
#'   Detector (LSD), a linear-time method that controls false detections and
#'   requires no parameter tuning. Based on Burns, Hanson, and Riseman's method
#'   with an a-contrario validation approach.
#'
#' @param img An `Image` object.
#' @param index A character string with the index to be used. Defaults to `"GRAY"`.
#' @param scale A positive numeric value. Scales the input image before detection using Gaussian filtering.
#'   A value <1 downscales, >1 upscales. Default is 0.8.
#' @param sigma_scale A positive numeric value determining the Gaussian filter sigma.
#'   If scale <1, sigma = sigma_scale / scale; otherwise, sigma = sigma_scale. Default is 0.6.
#' @param quant A positive numeric value controlling gradient quantization error. Default is 2.0.
#' @param ang_th A numeric value (0-180) defining the gradient angle tolerance in degrees. Default is 22.5.
#' @param log_eps A numeric detection threshold. Larger values make detection stricter. Default is 0.0.
#' @param density_th A numeric value (0-1) defining the minimum proportion of supporting points in a rectangle. Default is 0.7.
#' @param n_bins A positive integer specifying the number of bins for pseudo-ordering gradient modulus. Default is 1024.
#' @param union Logical. If TRUE, merges close line segments. Default is FALSE.
#' @param union_min_length Numeric. Minimum segment length to merge. Default is 5.
#' @param union_max_distance Numeric. Maximum distance between segments to merge. Default is 5.
#' @param union_ang_th Numeric. Angle threshold for merging segments. Default is 7.
#' @param union_use_NFA Logical. If TRUE, uses NFA in merging. Default is FALSE.
#' @param union_log_eps Numeric. Detection threshold for merging. Default is 0.0.
#'
#' @return A list of class `lsd` containing:
#' \itemize{
#'   \item `n` - Number of detected line segments.
#'   \item `lines` - A matrix with detected segments (columns: x1, y1, x2, y2, width, p, -log_nfa).
#'   \item `pixels` - A matrix assigning each pixel to a detected segment (0 = unused pixels).
#' }
#'
#' @references
#' Grompone von Gioi, R., Jakubowicz, J., Morel, J.-M., & Randall, G. (2010).
#' LSD: A Fast Line Segment Detector with a False Detection Control.
#' IEEE Transactions on Pattern Analysis and Machine Intelligence, 32(4), 722-732.\doi{10.5201/ipol.2012.gjmr-lsd}
#'
#' @export
#' @examples
#' library(pliman)
image_line_segment <- function(img,
                               index = "GRAY",
                               scale = 0.8,
                               sigma_scale = 0.6,
                               quant = 2.0,
                               ang_th = 22.5,
                               log_eps = 0.0,
                               density_th = 0.7,
                               n_bins = 1024,
                               union = FALSE,
                               union_min_length = 5,
                               union_max_distance = 5,
                               union_ang_th = 7,
                               union_use_NFA = FALSE,
                               union_log_eps = 0.0) {
  ind <- image_index(EBImage::flop(EBImage::transpose(img)), index, plot = FALSE)[[1]]
  x <- ind@.Data * 255
  lines <- detect_line_segments(as.numeric(x),
                                X = nrow(x),
                                Y = ncol(x),
                                scale = as.numeric(scale),
                                sigma_scale = as.numeric(sigma_scale),
                                quant = as.numeric(quant),
                                ang_th = as.numeric(ang_th),
                                log_eps = as.numeric(log_eps),
                                need_to_union = as.logical(union),
                                union_use_NFA = as.logical(union_use_NFA),
                                union_ang_th = as.numeric(union_ang_th),
                                union_log_eps = as.numeric(union_log_eps),
                                length_threshold = as.numeric(union_min_length),
                                dist_threshold = as.numeric(union_max_distance))

  names(lines) <- c("lines", "pixels")
  colnames(lines$lines) <- c("x1", "y1", "x2", "y2", "width", "p", "-log_nfa")
  lines$n <- nrow(lines$lines)

  return(lines)
}

#' @title Plot Detected Line Segments
#' @description Plots the detected line segments from the output of [image_line_segment()].
#' Each segment is drawn as a red line on the existing plot.
#'
#' @param x A list returned by [image_line_segment()], containing detected line segments.
#' @param col The color of lines
#' @param lwd The width of lines. Defaults to 1
#' @return No return value. The function adds line segments to an existing plot.
#'
#' @examples
#' library(pliman)
#'
#' @export
plot_line_segment <- function(x, col = "red", lwd = 1){
  a <- lapply(seq_len(x$n), FUN=function(i){
    l <- rbind(
      x$lines[i, c("x1", "y1")],
      x$lines[i, c("x2", "y2")])
    segments(l[1, 1], l[1, 2], l[2, 1], l[2, 2], col = col, lwd = lwd)
  })
}

#' @title Extract Mean Colors from a Color Checker Card
#'
#' @description
#' This function identifies a color checker card in an image, finds its
#' four corners to correct for perspective distortion, generates a grid
#' corresponding to the color patches, and extracts the mean RGB values from
#' the center of each patch.
#'
#' **Note:** The function attempts to automatically guess appropriate values
#' for the parameters related to the card's dimensions (`nrow`, `ncol`) and
#' the sampling process (`erode`, `xpix`, `ypix`) based on the detected
#' coverage area of the card in the image. You may override these guesses
#' by providing explicit values.
#'
#' @param img An `Image` object from the `EBImage` or `pliman` package.
#' @param index The vegetation or color index string (e.g., "GRAY", "B", "R")
#'   passed to `pliman::image_binary()` to segment the card from the
#'   background. Default is **"GRAY"**.
#' @param nrow The number of rows of color patches on the card. If **`NULL`**
#'   (default), the function automatically determines this based on the card's
#'   aspect ratio (6 for 6x4, 4 for 4x6).
#' @param ncol The number of columns of color patches on the card. If **`NULL`**
#'   (default), the function automatically determines this based on the card's
#'   aspect ratio (4 for 6x4, 6 for 4x6).
#' @param erode The size (in pixels) of the erosion kernel applied to the
#'   binary mask to shrink the selection and avoid patch edges. If **`NULL`**
#'   (default), a value is automatically calculated based on the card's
#'   coverage in the image.
#' @param xpix The total width (in pixels) of the sampling rectangle at the
#'   center of each patch. If **`NULL`** (default), a value is automatically
#'   calculated based on the card's coverage in the image.
#' @param ypix The total height (in pixels) of the sampling rectangle at the
#'   center of each patch. If **`NULL`** (default), a value is automatically
#'   calculated based on the card's coverage in the image.
#' @param plot Logical. If **`TRUE`** (default), displays the original image with
#'   detected corners (red dots), sampling boxes (red rectangles), and
#'   patch IDs.
#'
#' @return a `data.frame` with the RGB values for each color chip.
#' A `data.frame` with `nrow * ncol` rows and 4 columns:
#' \itemize{
#'   \item `id`: The patch identifier (from 1 to `nrow * ncol`).
#'   \item `R`: The mean Red channel value (typically 0-1).
#'   \item `G`: The mean Green channel value (typically 0-1).
#'   \item `B`: The mean Blue channel value (typically 0-1).
#' }
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("colorcheck.jpg")
#' get_card_colors(img, xpix = 30, ypix = 30)
#' }
get_card_colors <- function(img,
                            index = "GRAY",
                            nrow = NULL,
                            ncol = NULL,
                            erode = NULL,
                            xpix = NULL,
                            ypix = NULL,
                            plot = TRUE){

  cover <- seq(0.05, 0.8, length.out = 50)
  exp_model <- function(x, a, b, c) {
    a + b * exp(-c * x)
  }
  size <- exp_model(cover, 3, 12, 40)
  erodeval  <- exp_model(cover, 30, -30, 5) * 1.8
  cex <-  seq(1.5, 0.4, length.out = 50)
  cext <-  seq(0.6, 1.5, length.out = 50)
  erodfact <- seq(0.5, 1, length.out = 50)
  bin <-
    image_binary(img,
                 index = index,
                 fill_hull = TRUE,
                 plot = FALSE)[[1]]
  #
  lab <- EBImage::bwlabel(bin)
  area <-
    lab |>
    EBImage::computeFeatures.shape()

  id_ref <- which.max(area[, 1])
  lab@.Data[which(lab != id_ref)] <- 0
  lab@.Data[which(lab == id_ref)] <- 1
  coverage <- sum(lab) / length(lab)
  configp <- which.min(abs(cover - coverage))
  suggestpix <- round((coverage * 100) * size[configp])
  erode <- ifelse(is.null(erode), erodeval[configp], erode)
  xpix <- ifelse(is.null(xpix), suggestpix, xpix)
  ypix <- ifelse(is.null(ypix), suggestpix, ypix)
  lab <- image_erode(lab, size = erode)
  oco <- EBImage::ocontour(lab)[[1]]
  hull_idx <- chull(oco)
  corners <- oco[hull_idx, ]
  order_idx <- order(atan2(corners[,2] - mean(corners[,2]),
                           corners[,1] - mean(corners[,1])))
  corners_ordered <- corners[order_idx, ]
  soma <- corners_ordered[, 1] + corners_ordered[, 2]
  diferenca <- corners_ordered[, 1] - corners_ordered[, 2]
  idx_bottom_left <- which.min(soma)
  idx_top_right <- which.max(soma)
  idx_bottom_right <- which.max(diferenca)
  idx_top_left <- which.min(diferenca)
  bottom_left <- corners_ordered[idx_bottom_left, ]
  top_right <- corners_ordered[idx_top_right, ]
  bottom_right <- corners_ordered[idx_bottom_right, ]
  top_left <- corners_ordered[idx_top_left, ]
  corners <- rbind(bottom_left, top_right, bottom_right, top_left)
  dblbr <- sqrt(sum((corners[1, ] - corners[3, ])^2) )
  dbltl <- sqrt(sum((corners[1, ] - corners[4, ])^2) )

  if (dblbr > dbltl) {
    nrow <- ifelse(is.null(nrow), 4, nrow)
    ncol <- ifelse(is.null(ncol), 6, ncol)
  } else {
    nrow <- ifelse(is.null(nrow), 6, nrow)
    ncol <- ifelse(is.null(ncol), 4, ncol)
  }
  grid_points <- matrix(NA, nrow = nrow * ncol, ncol = 2)
  colnames(grid_points) <- c("x", "y")
  for (r in 1:nrow) {
    for (c in 1:ncol) {
      c_norm <- (c - 0.5) / ncol
      r_norm <- (nrow - r + 0.5) / nrow # Close to 1 for r=1 (bottom edge), close to 0 for r=nrow (top edge)
      p_top_x <- corners[4, 1] * (1 - c_norm) + corners[2, 1] * c_norm
      p_top_y <- corners[4, 2] * (1 - c_norm) + corners[2, 2] * c_norm
      # p_bottom is on the line segment connecting BL (1) and BR (3)
      p_bottom_x <- corners[1, 1] * (1 - c_norm) + corners[3, 1] * c_norm
      p_bottom_y <- corners[1, 2] * (1 - c_norm) + corners[3, 2] * c_norm
      grid_x <- p_top_x * (1 - r_norm) + p_bottom_x * r_norm
      grid_y <- p_top_y * (1 - r_norm) + p_bottom_y * r_norm

      idx <- (r - 1) * ncol + c
      grid_points[idx, ] <- c(grid_x, grid_y)
    }
  }
  vals <-
    grid_points |>
    as.data.frame() |>
    mutate(lab = 1:(nrow * ncol))

  # criar grid para extrair cor
  coords <-
    lapply(1:nrow(vals), function(i){
      x <- vals[i, 1]
      y <- vals[i, 2]
      rx <- c(x - xpix / 2, x + xpix / 2)
      ry <- c(y - ypix / 2, y + ypix / 2)
      c(rx, ry)
    })

  rgbobs <-
    do.call(rbind, lapply(coords, function(x){
      apply(img@.Data[x[1]:x[2], x[3]:x[4], ], 3, mean)
    })) |>
    as.data.frame() |>
    dplyr::mutate(id = 1:(nrow * ncol), .before = 1)

  colnames(rgbobs) <- c("id", "R", "G", "B")
  if(plot){
    plot(img)
    coords_matrix <- do.call(rbind, coords)
    colnames(coords_matrix) <- c("x_min", "x_max", "y_min", "y_max")
    rect(
      xleft   = coords_matrix[, "x_min"],
      ybottom = coords_matrix[, "y_min"],
      xright  = coords_matrix[, "x_max"],
      ytop    = coords_matrix[, "y_max"],
      border = "red",  # Cor da borda
      lwd = 1.5       # Largura da linha
    )
    points(corners, col = "red", pch = 19, cex = cex[configp])
    text(x = vals[, 1], y = vals[, 2], labels = vals[, 3],
         cex = cext[configp])
  }
  return(rgbobs)
}


#' @title Correct Image Colors using a Color Checker
#'
#' @description
#' Calibrates the color of an image using a set of known color references (e.g.,
#' from a color checker) by implementing polynomial color correction models.
#' It works by finding a transformation matrix (K) that maps the
#' observed colors (sampled from the color card) to their known reference values.
#' This matrix K is then applied to every pixel in the image.
#'
#' The correction model is based on solving the equation \eqn{S \times K = T}, where:
#' \itemize{
#'   \item \code{S} is the source matrix of observed colors extended by polynomial terms.
#'   Dimensions are \eqn{n \times 9} for \code{model = "cubic"} or \eqn{n \times 20}
#'   for \code{model = "root_polynomial"}.
#'   \item \code{T} is the target matrix (\eqn{n \times 3}) of known reference colors (R, G, B).
#'   \item \code{K} is the transformation matrix (\eqn{P \times 3}) solved using the
#'   Moore-Penrose pseudo-inverse.
#' }
#'
#' @param img An `Image` object to be corrected.
#' @param card_colors A data frame of the *observed* colors from the color
#'   checker in the image. Typically the output from `get_card_colors()`.
#'   Must contain 'R', 'G', and 'B' columns.
#' @param known_colors A data frame of the *known* reference colors for the
#'   color checker. Must have 'R', 'G', and 'B' columns and the same
#'   number of rows as `card_colors`.
#' @param k_mat A pre-calculated transformation matrix (**K**).
#'   If provided, \code{card_colors} and \code{known_colors} are ignored, and the
#'   correction is applied directly using this matrix.
#'   This allows reusing a calculated matrix for multiple images taken in the same
#'   light conditions.
#'   \itemize{
#'     \item If \code{model = "cubic"}, \code{k_mat} must be a \eqn{9 \times 3} matrix.
#'     \item If \code{model = "root_polynomial"}, \code{k_mat} must be a \eqn{20 \times 3} matrix.
#'   }
#'   Defaults to \code{NULL}.
#' @param model The correction model to use.
#'   \itemize{
#'     \item \code{"cubic"} (default): Uses a 3rd-order polynomial with 9 terms
#'     (R, G, B, R^2, G^2, B^2, RG, RB, GB).
#'     \item \code{"root_polynomial"}: Uses a root-polynomial model with 20 terms,
#'     which usually provides higher accuracy for non-linear color distortions.
#'   }
#' @return
#' If both \code{card_colors} and \code{known_colors} are provided (or if \code{k_mat} is \code{NULL}
#' and the required color data frames are present), a list is returned:
#' \itemize{
#'   \item **img**: An \code{Image} object with the color correction applied,
#'     maintaining the original dimensions and color mode (\code{'Color'}).
#'   \item **k**: The calculated transformation matrix (**K**). Dimensions are
#'   \eqn{9 \times 3} or \eqn{20 \times 3} depending on the \code{model}.
#' }
#' If only \code{img} and a pre-calculated \code{k_mat} are provided, the function
#' returns only the corrected \code{Image} object.
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("colorcheck.jpg")
#'
#' # Standard Macbeth ColorChecker values (approximate)
#' kvals <- data.frame(
#'  id = 1:24,
#'  R = c(0.976, 0, 0.870, 0.384, 0.792, 0.752, 0.227, 0.494, 0.631, 0.960,
#'        0.764, 0.321, 0.478, 0.729, 0.325, 0.341, 0.313, 0.223, 0.615, 0.772,
#'        0.168, 0.098, 0.933, 0.439),
#'  G = c(0.949, 0.498, 0.462, 0.733, 0.776, 0.294, 0.345, 0.490, 0.615, 0.803,
#'        0.309, 0.415, 0.462, 0.101, 0.227, 0.470, 0.313, 0.572, 0.737, 0.568,
#'        0.160, 0.215, 0.619, 0.298),
#'  B = c(0.933, 0.623, 0.125, 0.650, 0.764, 0.568, 0.623, 0.682, 0.603, 0,
#'        0.372, 0.235, 0.454, 0.200, 0.415, 0.607, 0.305, 0.250, 0.211, 0.490,
#'        0.168, 0.529, 0.098, 0.235)
#'  )
#'
#' card_colors <- get_card_colors(img, xpix = 20, ypix = 30, erode = 20)
#'
#' # Correct using the default cubic model
#' corrected <- image_correction(img, card_colors, kvals, model = "root-polynomial")
#'
#' # Compare original and corrected
#' image_combine(img, corrected$img)
#' }
image_correction <- function(img,
                             card_colors = NULL,
                             known_colors = NULL,
                             k_mat = NULL,
                             model = "cubic"){
  if (!model %in% c("cubic", "root_polynomial")) {
    cli::cli_abort(
      c("Argument {.val model} must be one of {.val cubic} or {.val root_polynomial}.",
        "x" = "Value provided: {.val {model}}")
    )
  }
  if(!is.null(card_colors) & !is.null(known_colors)){
    # Check if known_colors has the correct number of rows
    if (nrow(card_colors) != nrow(known_colors)) {
      cli::cli_abort(
        c("The number of rows in {.arg card_colors} ({.val {nrow(card_colors)}}) does not equal {.arg known_colors} ({.val {nrow(known_colors)}}).",
          "i" = "Please provide the correct reference color data frame.")
      )
    }

    # Check if known_colors has R, G, B columns
    if (!all(c("R", "G", "B") %in% colnames(known_colors))) {
      cli::cli_abort(
        c("The {.arg known_colors} data frame must contain columns {.field R}, {.field G}, and {.field B}.",
          "x" = "Found columns: {.field {colnames(known_colors)}}")
      )
    }
    K <- mpinv(create_poly_matrix(card_colors, model = model)) %*% as.matrix(known_colors[, c("R", "G", "B")])
    # 4. Aplicar a correção no C++
    if (model == "cubic") {
      if(nrow(K) != 9){
        cli::cli_abort("Para o modelo 'cubic', 'k_mat' deve ter 9 linhas.")
      }
    } else if (model == "root_polynomial") {
      if(nrow(K) != 20){
        cli::cli_abort("Para o modelo 'root-polynomial', 'k_mat' deve ter 20 linhas.")
      }
    }

    return(list(
      img = EBImage::Image(correct_image_rcpp(img, K, model = model),
                           dim = dim(img),
                           colormode = 'Color'),
      k = K
    ))
  } else{
    return(EBImage::Image(correct_image_rcpp(img, k_mat, model = model),
                          dim = dim(img),
                          colormode = 'Color'))
  }
}


#' Interactive Color Correction
#'
#' @description
#' Interactively calibrates the color of an image by prompting the user to
#' sample color patches. It first calls `pliman::pick_rgb_area()` to pause and
#' prompt the user to click on `color_chips` patches in the image. These
#' interactively sampled colors are then used as the *observed* colors to
#' compute the 3rd-order polynomial transformation matrix (K), which maps the
#' observed colors to the provided `known_colors`. This matrix K is then applied
#' to every pixel in the image.
#'
#' @param img An `Image` object to be corrected.
#' @param known_colors A `data.frame` containing the target reference values.
#'   Must contain the columns `R`, `G`, and `B` (normalized between 0-1) and
#'   the number of rows must equal `color_chips`.
#' @param color_chips The number of color patches to be interactively
#'   sampled from the image. Defaults to `24`.
#'
#' @return An `Image` object with corrected colors.
#'
#' @export
#' @examples
#' if(interactive()){
#' # colorcheck image available at https://github.com/HarryCWright/PlantSizeClr
#' library(pliman)
#'
#' # known values in the following sequence (row, column)
#' # blue (6, 2)
#' # red (4, 2)
#' # yellow (3, 2)
#' # white (1, 1)
#' known_colors <- data.frame(
#'  id = c(1, 2, 3, 4),
#'  R = c(0.09803922, 0.72941176, 0.96078431, 0.97647059),
#'  G = c(0.2156863, 0.1019608, 0.8039216, 0.9490196),
#'  B = c(0.5294118, 0.2000000, 0.0000000, 0.9333333)
#' )
#' img <- image_pliman("colorcheck.jpg")
#' # draw four samples, following the order (blue, red, yellow, white)
#'  img_cor <- image_correction_pick(img, known_colors, color_chips = 4)
#'  image_combine(img, img_cor)
#' }
#'
#'
image_correction_pick <- function(img,
                                  known_colors,
                                  color_chips = 24){
  # Check if known_colors has the correct number of rows
  if (nrow(known_colors) != color_chips) {
    cli::cli_abort(
      c("The number of rows in {.arg known_colors} ({.val {nrow(known_colors)}}) does not equal {.arg color_chips} ({.val {color_chips}}).",
        "i" = "Please provide the correct reference color data frame.")
    )
  }
  if (!all(c("R", "G", "B") %in% colnames(known_colors))) {
    cli::cli_abort(
      c("The {.arg known_colors} data frame must contain columns {.field R}, {.field G}, and {.field B}.",
        "x" = "Found columns: {.field {colnames(known_colors)}}")
    )
  }
  cli::cli_progress_step(
    msg        = "Pick the {.val {color_chips}} color chips in the image that correspond to the {.arg known_colors} RGB values...",
    msg_done   = "Color sampling done",
    msg_failed = "Oops, something went wrong."
  )
  rgbimg <- pick_rgb_area(img, n = color_chips, verbose = FALSE)
  cli::cli_progress_step(
    msg        = "Performing color correction...",
    msg_done   = "Color correction finished",
    msg_failed = "Oops, something went wrong."
  )
  K <- mpinv(create_poly_matrix(rgbimg)) %*% create_poly_matrix(known_colors)
  return(EBImage::Image(correct_image_rcpp(img, K),
                        dim = dim(img),
                        colormode = 'Color'))
}

