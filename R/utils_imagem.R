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
#' @param mar,oma Margins for each panel and outer margins for the combined plot.
#'   Defaults to `mar = c(1.5, 1.5, 1.5, 1.5)` and `oma = c(0, 0, 0, 0)`. Set
#'   `mar = 0` or `mar = c(0, 0, 0, 0)` to eliminate margins and maximize the
#'   plotting area.
#' @importFrom stats reshape IQR quantile
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A grid with the images in `...`
#' @examples
#' if (interactive()) {
#' library(pliman)
#' img1 <- image_pliman("sev_leaf.jpg")
#' img2 <- image_pliman("sev_leaf_nb.jpg")
#' image_combine(img1, img2)
#'
#' # Remove margins to optimize plotting area:
#' image_combine(img1, img2, mar = 0)
#' }
image_combine <- function(...,
                          labels = NULL,
                          nrow = NULL,
                          ncol = NULL,
                          col = "black",
                          verbose = TRUE,
                          mar = c(1.5, 1.5, 1.5, 1.5),
                          oma = c(0, 0, 0, 0)) {
  dots <- list(...)
  call_args <- match.call(expand.dots = FALSE)$...

  get_img <- function(x) {
    if (inherits(x, c("image", "Image")) || is.matrix(x) || is.array(x)) {
      return(x)
    }
    if (is.list(x) && !is.null(x$img)) {
      return(x$img)
    }
    x
  }

  is_img_list <- function(lst) {
    if (!is.list(lst) || length(lst) == 0L) return(FALSE)
    first <- lst[[1L]]
    inherits(first, c("image", "Image")) || (is.list(first) && !is.null(first$img))
  }

  if (length(dots) == 1L && is_img_list(dots[[1L]])) {
    plots <- dots[[1L]]
    if (inherits(plots, c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list"))) {
      plots <- lapply(plots, function(x) x[[1L]])
    }
  } else {
    plots <- dots
  }

  plots <- lapply(plots, get_img)

  if (!is.null(labels)) {
    names(plots) <- labels
  } else if (is.null(names(plots)) && !is.null(call_args)) {
    names(plots) <- sapply(call_args, deparse)
  }

  num_plots <- length(plots)
  if (num_plots == 0L) {
    return(return(NULL))
  }

  if (is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(num_plots))
    nrow <- ceiling(num_plots / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(num_plots / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(num_plots / ncol)
  }

  if (length(mar) == 1L) {
    mar <- rep(mar, 4L)
  }
  if (length(oma) == 1L) {
    oma <- rep(oma, 4L)
  }

  op <- par(mfrow = c(nrow, ncol), mar = mar, oma = oma)
  on.exit(par(op), add = TRUE)

  idx_names <- if (is.null(names(plots))) as.character(seq_len(num_plots)) else names(plots)

  for (i in seq_len(num_plots)) {
    img <- plots[[i]]
    plot(img)
    if (isTRUE(verbose)) {
      lbl <- idx_names[i]
      if (nzchar(lbl)) {
        d <- dim(img)
        h <- if (!is.null(d) && length(d) >= 2L) d[2L] else 100
        text(0, h * 0.075, lbl, pos = 4, col = col, cex = 1.1, font = 2)
      }
    }
  }
  return(plots)
}

# Internal helper function to read a single image file (Zero-dependency C++ engine)
.read_image_file <- function(file, ...) {
  if (grepl("^https?://", file, ignore.case = TRUE)) {
    ext <- tools::file_ext(file)
    if (ext == "") ext <- "jpg"
    tmp <- tempfile(fileext = paste0(".", ext))
    on.exit(unlink(tmp), add = TRUE)
    utils::download.file(file, tmp, mode = "wb", quiet = TRUE)
    file <- tmp
  }

  ext <- tolower(tools::file_ext(file))

  # 1. Spatial Rasters (.gri, .grd)
  if (ext %in% c("gri", "grd")) {
    if (requireNamespace("terra", quietly = TRUE)) {
      r <- terra::rast(file)
      img_arr <- terra::as.array(r)
      if (max(img_arr, na.rm = TRUE) <= 1.0) {
        raw_data <- float_to_raw_cpp(img_arr)
      } else {
        raw_data <- as.raw(pmin(pmax(round(img_arr), 0), 255))
        dim(raw_data) <- dim(img_arr)
      }
      return(as_image(raw_data))
    } else {
      cli::cli_abort(c(
        "!" = "Package {.pkg terra} is required to import spatial raster files ({.val .{ext}}).",
        "i" = "Please install it with: {.code install.packages('terra')}"
      ))
    }
  }

  # 2. WebP (.webp)
  if (ext == "webp") {
    if (requireNamespace("webp", quietly = TRUE)) {
      img_arr <- webp::read_webp(file, ...)
      if (length(dim(img_arr)) == 3L && dim(img_arr)[3L] == 4L) {
        img_arr <- img_arr[, , 1:3, drop = FALSE]
      }
      raw_data <- float_to_raw_cpp(img_arr)
      return(as_image(raw_data, colormode = if (length(dim(raw_data)) == 3L) "Color" else "Grayscale"))
    } else {
      cli::cli_abort(c(
        "!" = "Package {.pkg webp} is required to import {.val .webp} image files.",
        "i" = "Please install it with: {.code install.packages('webp')}"
      ))
    }
  }

  # 3. Native C++ Engine (JPEG, PNG, TIFF, BMP, TGA, GIF, etc.) with Fallback for Complex TIFFs
  out <- tryCatch({
    read_image_cpp(file)
  }, error = function(e) {
    if (ext %in% c("tif", "tiff")) {
      if (requireNamespace("tiff", quietly = TRUE)) {
        img_arr <- tiff::readTIFF(file, native = FALSE, ...)
        if (is.list(img_arr)) img_arr <- img_arr[[1L]]
        if (length(dim(img_arr)) == 3L && dim(img_arr)[3L] == 4L) {
          img_arr <- img_arr[, , 1:3, drop = FALSE]
        }
        if (is.raw(img_arr)) {
          return(as_image(img_arr, colormode = if (length(dim(img_arr)) == 3L) "Color" else "Grayscale"))
        }
        raw_data <- float_to_raw_cpp(img_arr)
        return(as_image(raw_data, colormode = if (length(dim(raw_data)) == 3L) "Color" else "Grayscale"))
      } else if (requireNamespace("terra", quietly = TRUE)) {
        r <- terra::rast(file)
        img_arr <- terra::as.array(r)
        if (max(img_arr, na.rm = TRUE) <= 1.0) {
          raw_data <- float_to_raw_cpp(img_arr)
        } else {
          raw_data <- as.raw(pmin(pmax(round(img_arr), 0), 255))
          dim(raw_data) <- dim(img_arr)
        }
        return(as_image(raw_data, colormode = if (length(dim(raw_data)) == 3L && dim(raw_data)[3L] >= 3) "Color" else "Grayscale"))
      } else {
        cli::cli_abort(c(
          "!" = "Failed to decode complex/compressed TIFF file {.val {basename(file)}} with native engine.",
          "i" = "Package {.pkg tiff} or {.pkg terra} is required for complex TIFF files.",
          "i" = "Please install with: {.code install.packages('tiff')}"
        ))
      }
    }
    cli::cli_abort(c(
      "!" = "Failed to import image {.val {basename(file)}}.",
      "x" = conditionMessage(e)
    ))
  })

  return(out)
}

# Internal helper function to write a single image file (Zero-dependency C++ engine)
.write_image_file <- function(img, file, quality = 95, ...) {
  dir_out <- dirname(file)
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  ext <- tolower(tools::file_ext(file))

  # 1. Native C++ Engine (JPEG, PNG, TIFF, BMP, TGA)
  if (ext %in% c("jpg", "jpeg", "png", "bmp", "tif", "tiff")) {
    write_image_cpp(img, file, quality = quality)
    return(return(file))
  }

  # 2. WebP (.webp)
  if (ext == "webp") {
    if (requireNamespace("webp", quietly = TRUE)) {
      img_d <- image_data(img)
      float_data <- if (is.raw(img_d)) raw_to_float_cpp(img_d) else img_d
      webp::write_webp(float_data, target = file, quality = round(quality * 100))
      return(return(file))
    } else {
      cli::cli_abort(c(
        "!" = "Package {.pkg webp} is required to export {.val .webp} image files.",
        "i" = "Please install it with: {.code install.packages('webp')}"
      ))
    }
  }

  # Native C++ fallback
  write_image_cpp(img, file, quality = quality)
  return(return(file))
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
#'  * For `image_input()` further arguments passed on to [as_image()].
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_import()` returns a new `image` object.
#' * `image_export()` returns an return vector of file names.
#' * `image_pliman()` returns a new `image` object with the example image
#' required. If an empty call is used, the path to the `tmp_images` directory
#' installed with the package is returned.
#' @examples
#' if (interactive()) {
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
  valid_extens <- c("png", "jpeg", "jpg", "tiff", "tif", "webp", "bmp",
                    "PNG", "JPEG", "JPG", "TIFF", "TIF", "WEBP", "BMP", "gri", "grd")
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
    list_img <- lapply(imgs, function(x) .read_image_file(x, ...))
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
    return(list_img)
  } else {
    resolved_img <- sapply(seq_along(img), function(idx) {
      img_item <- img[idx]
      if (grepl("http", img_item, fixed = TRUE)) {
        return(img_item)
      }

      item_dir <- ifelse(is.null(path), dirname(img_item), path)
      item_base <- file_name(img_item)
      item_ext <- tools::file_ext(img_item)

      priority_exts <- c("jpg", "jpeg", "png", "tiff", "tif", "webp", "bmp", "gri", "grd",
                         "JPG", "JPEG", "PNG", "TIFF", "TIF", "WEBP", "BMP", "GRI", "GRD")

      if (item_ext == "" || !(item_ext %in% priority_exts)) {
        dir_files <- list.files(item_dir)
        matching_files <- dir_files[tools::file_path_sans_ext(dir_files) == item_base]
        found_matches <- matching_files[tools::file_ext(matching_files) %in% priority_exts]

        if (length(found_matches) > 0) {
          match_exts <- tools::file_ext(found_matches)
          best_idx <- min(match(match_exts, priority_exts), na.rm = TRUE)
          best_ext <- priority_exts[best_idx]
          best_file <- found_matches[tools::file_ext(found_matches) == best_ext][1]

          cli::cli_inform("Extension not provided or unmatched for {.val {item_base}}. Importing {.file {best_file}}.")
          return(paste0(item_dir, "/", best_file))
        }
      }

      if (!is.null(path)) {
        return(paste0(path, "/", basename(img_item)))
      } else {
        if (dirname(img_item) == ".") {
          return(paste0(item_dir, "/", img_item))
        }
        return(img_item)
      }
    })

    img_name <- resolved_img

    if(!any(grepl("http", img_name, fixed = TRUE))) {
      test <- file.exists(img_name)
      if (!all(test)) {
        cli::cli_abort("Image {.val {basename(img_name[which(test == FALSE)])}} not found in {.dir {dirname(img_name[which(test == FALSE)])}}.")
      }
    }

    if(length(img_name) > 1){
      ls <- lapply(seq_along(img_name), function(x) .read_image_file(img_name[x], ...))
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
      return(ls)
    } else {
      img_out <- .read_image_file(img_name, ...)
      if(isTRUE(plot)){
        plot(img_out)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          cli::cli_abort("Argument {.val resize} must be numeric.")
        }
        img_out <- image_resize(img_out, resize)
      }
      return(img_out)
    }
  }
}

#' @export
#' @name utils_image
image_export <- function(img,
                         name = NULL,
                         prefix = "",
                         extension = NULL,
                         subfolder = NULL,
                         ...){
  if(inherits(img, c("binary_list", "index_list",
                       "img_mat_list", "palette_list"))){
    img <- lapply(img, function(x){x[[1]]})
  }
  if(inherits(img, "segment_list")){
    img <- lapply(img, function(x){x[[1]][[1]]})
  }
  if(is.list(img)){
    if(!all(sapply(img, is_image))){
      cli::cli_abort("All images must be of class {.code image} or {.code Image}.")
    }

    if (!missing(name) && is.character(name) && length(name) == 1L && is.null(subfolder)) {
      if (dir.exists(name) || nchar(tools::file_ext(name)) == 0L) {
        subfolder <- name
        name <- NULL
      }
    }

    if (is.null(name) || length(name) == 0L) {
      name <- names(img)
    }
    if (is.null(name) || length(name) == 0L) {
      name <- paste0("img_", seq_along(img))
    }

    base_names <- file_name(name)
    extens <- unlist(file_extension(name))

    if ((length(extens) == 0L || all(nchar(extens) == 0L)) && is.null(extension)) {
      extens <- rep("jpg", length(img))
      n_img <- length(img)
      cli::cli_inform(c("v" = "{n_img} image{?s} exported as {.val *.jpg} file{?s}."))
    } else if (!is.null(extension)) {
      extens <- rep(extension, length(img))
    } else if (length(extens) < length(img)) {
      extens <- rep(extens, length.out = length(img))
    }

    if(!missing(subfolder) && !is.null(subfolder)){
      dir_out <- subfolder
      if(dir.exists(dir_out) == FALSE){
        dir.create(dir_out, recursive = TRUE)
      }
      out_names <- file.path(dir_out, paste0(prefix, base_names, ".", extens))
    } else {
      out_names <- paste0(prefix, base_names, ".", extens)
    }
    out_names <- gsub("/+", "/", out_names)
    out_names <- sub("^\\./", "", out_names)
    lapply(seq_along(img), function(i){
      .write_image_file(img[[i]], out_names[i], ...)
    })
    return(invisible(out_names))
  } else {
    filname <- file_name(name)
    extens <- unlist(file_extension(name))
    dir_out <- file_dir(name)

    if (is.null(filname) || nchar(filname) == 0L) {
      filname <- "img"
    }

    if(length(extens) == 1 && nchar(extens) > 0){
      extens <- extens
    } else if((length(extens) == 0 || nchar(extens) == 0) && is.null(extension)){
      extens <- "jpg"
      cli::cli_inform(c("v" = "Image exported as {.val *.jpg} file."))
    } else if(!is.null(extension)){
      extens <- extension
    }

    if(!missing(subfolder) && !is.null(subfolder)){
      dir_out <- subfolder
    } else if (dir_out == ".//" || dir_out == "./" || dir_out == ".") {
      dir_out <- ""
    }

    if(nchar(dir_out) > 0 && dir.exists(dir_out) == FALSE){
      dir.create(dir_out, recursive = TRUE)
    }

    if (nchar(prefix) > 0) {
      filname <- paste0(prefix, filname)
    }

    out_file <- if (nchar(dir_out) > 0) file.path(dir_out, paste0(filname, ".", extens)) else paste0(filname, ".", extens)
    out_file <- gsub("/+", "/", out_file)
    out_file <- sub("^\\./", "", out_file)

    .write_image_file(img, out_file, ...)
    return(invisible(out_file))
  }
}

#' @export
#' @name utils_image
image_input <- function(img, ...){
  if(inherits(img, "character")){
    image_import(img, ...)
  } else if(inherits(img, "array") || inherits(img, "matrix")){
    range <- max(img, na.rm = TRUE)
    if(range > 1){
      as_image(img / 255, colormode = "Color")
    } else {
      as_image(img, colormode = "Color")
    }
  } else if(inherits(img, c("image", "Image"))) {
    as_image(img)
  }
}

#' @export
#' @name utils_image
image_pliman <- function(img, plot = FALSE){
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
    return(im)
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
#' * `image_resize()` Resize the `image`. See more at [image_resize()].
#' * `image_contrast()` Improve contrast locally by performing adaptive
#' histogram equalization.
#' * `image_dilate()` Performs image dilatation.
#' * `image_erode()` Performs image erosion.
#' * `image_opening()` Performs an erosion followed by a dilation.
#' * `image_closing()` Performs a dilation followed by an erosion.
#' * `image_filter()` Performs median filtering in constant time.
#' * `image_blur()` Performs blurring filter of images.
#' * `image_skeleton()` Performs image skeletonization.
#'
#'
#' @name utils_transform
#' @inheritParams image_view
#' @inheritParams analyze_objects
#' @param img An image or a list of images of class `image`.
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
#' @param bg_col,bg_color Color used to fill the background pixels, defaults to `"white"`.
#' @param nx Number of contextual regions in the x-direction (tile grid) for CLAHE contrast enhancement. Defaults to 8.
#' @param ny Number of contextual regions in the y-direction (tile grid) for CLAHE contrast enhancement. Defaults to 8.
#' @param clip_limit Normalized contrast limit value for CLAHE contrast enhancement. Defaults to 3.
#' @param bins Number of histogram bins used for CLAHE contrast enhancement. Defaults to 256.
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
#' @param kern An `image` object or an array, containing the structuring
#'   element. Defaults to a brushe generated with [make_brush()].
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
#' * `image_skeleton()` returns a binary `image` object.
#' * All other functions returns a  modified version of `image` depending on the
#' `image_*()` function used.
#' * If `image` is a list, a list of the same length will be returned.
#' @examples
#' if (interactive()) {
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
                           opening = FALSE,
                           closing = FALSE,
                           filter = FALSE,
                           invert = FALSE,
                           threshold = "Otsu",
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  if(is.list(img)){
    if(inherits(img, c("binary_list", "segment_list", "index_list",
                         "img_mat_list", "palette_list"))){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, is_image))){
      cli::cli_abort("All images must be of class {.code image} or {.code Image}.")
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
    return(structure(res, class = "autocrop_list"))
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
    return(segmented)
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
    if(inherits(img, c("binary_list", "segment_list", "index_list",
                         "img_mat_list", "palette_list"))){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, is_image))){
      cli::cli_abort("All images must be of class {.code image} or {.code Image}.")
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
      if(verbose){
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }
    } else{
      res <- lapply(img, image_crop, width, height, viewer, downsample, max_pixels)
    }
    return(res)
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
      img <- if (length(dim(img)) == 3) img[width, height, , drop = FALSE] else img[width, height, drop = FALSE]
    }
    if (is.null(width) & is.null(height)) {
      if(vieweropt == "base"){
        cli::cli_inform(c("i" = "Use the {cli::col_blue(cli::style_bold('left mouse button'))} to crop the image."))

        if(length(dim(img)) >= 3 && dim(img)[3] >= 3){
          plot(img[, , 1:3])
        } else {
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
      img <- if (length(dim(img)) == 3) img[w, h, , drop = FALSE] else img[w, h, drop = FALSE]
      if(isTRUE(verbose)){
        cat(paste0("width = ", w[1], ":", w[length(w)]), "\n")
        cat(paste0("height = ", h[1], ":", h[length(h)]), "\n")
      }
    }
    if (isTRUE(plot)) {
      if(dim(img)[[3]] > 2){
        plot(as_image(img[,,1:3], colormode = "Color"))
      } else if(dim(img)[[3]] == 1){
        plot(img)
      }
    }
    return(img)
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

    return(res)

  } else {
    width  <- dim(img)[[1]]
    height <- dim(img)[[2]]

    if (verbose) {
      cli::cli_rule("Image dimension")
      cli::cli_text("Width : {.val {width}}")
      cli::cli_text("Height: {.val {height}}")
    }

    return(list(width = width, height = height))
  }
}

#' @name utils_transform
#' @export
image_rotate <- function(img,
                         angle,
                         bg_col = "white",
                         bg_color = NULL,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = TRUE) {

  if (!is.null(bg_color)) {
    bg_col <- bg_color
  }

  # -- inner helper ---------------------------------------------------------
  # Convert bg_col to a per-channel numeric vector [0, 1] for the C++ kernel.
  # Then call image_rotate_cpp (backward-mapping bilinear, OpenMP) and wrap
  .one <- function(im) {
    cm    <- attr(im, "colormode") %||% "Color"
    nch   <- if (length(dim(im)) == 3) dim(im)[3] else 1L
    bg_v  <- if (is.numeric(bg_col)) {
      bg_col
    } else {
      tryCatch(
        grDevices::col2rgb(bg_col)[, 1] / 255,
        error = function(e) rep(1, nch)
      )
    }
    arr <- image_rotate_cpp(image_data(im), angle, bg_v)
    as_image(arr, colormode = cm)
  }

  if (is.list(img) && !is_image(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }
    if (!all(sapply(img, is_image))) {
      cli::cli_abort("All images must be of class {.cls image}.")
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

      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) pliman::image_rotate(im, angle, bg_col)
      )[.progress]
      if(verbose){
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }
    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Rotating {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images sequentially...",
          msg_done   = "Processing complete",
          msg_failed = "Sequential processing failed"
        )
      }

      res <- lapply(img, .one)
    }

    if (isTRUE(plot)) {
      for (r in res) {
        plot(r)
      }
    }

    return(res)

  } else {
    rotated <- .one(img)

    if (isTRUE(plot)) {
      plot(rotated)
    }

    return(rotated)
  }
}


#' @name utils_transform
#' @export
image_horizontal <- function(img,
                             parallel = FALSE,
                             workers = NULL,
                             verbose = TRUE,
                             plot = FALSE) {


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

      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          w <- dim(im)[[1]]
          h <- dim(im)[[2]]
          if (w < h) {
            image_rotate(im, angle = 90)
          } else {
            im
          }
        }
      )[.progress]
      if(verbose){
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

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
          image_rotate(im, angle = 90)
        } else {
          im
        }
      })
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(res)

  } else {
    w <- dim(img)[[1]]
    h <- dim(img)[[2]]
    if (w < h) {
      img <- image_rotate(img, angle = 90)
    }

    if (isTRUE(plot)) {
      plot(img)
    }

    return(img)
  }
}

#' @name utils_transform
#' @export
image_vertical <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {


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
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = function(im) {
          w <- dim(im)[[1]]
          h <- dim(im)[[2]]
          if (w > h) {
            image_rotate(im, angle = 90)
          } else {
            im
          }
        }
      )[.progress]
      if(verbose){
        cli::cli_progress_step(
          msg        = "Processing {.val {length(img)}} images in parallel...",
          msg_done   = "Batch processing finished",
          msg_failed = "Oops, something went wrong."
        )
      }

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
          image_rotate(im, angle = 90)
        } else {
          im
        }
      })
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(res)

  } else {
    w <- dim(img)[[1]]
    h <- dim(img)[[2]]
    if (w > h) {
      img <- image_rotate(img, angle = 90)
    }

    if (isTRUE(plot)) {
      plot(img)
    }

    return(img)
  }
}

#' @name utils_transform
#' @export
image_hreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {


  # -- inner helper: reflect one Image object ---------------------------------
  # 1. Extract image_data() (O(1) reference - no copy until modified).
  # 2. Call the optimised C++ kernel (allocates exactly one output buffer).
  # 3. Wrap result with as_image() - O(1) S4 constructor, no data copy.
  # Total extra allocation: exactly one buffer of the image size.
  .one <- function(im) {
    as_image(
      image_hreflect_cpp(image_data(im)),
      colormode = attr(im, "colormode") %||% "Color"
    )
  }

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
        .f = function(im) pliman::image_hreflect(im)
      )[.progress]

    } else {
      res <- lapply(img, .one)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(res)

  } else {
    img <- .one(img)
    if (isTRUE(plot)) plot(img)
    return(img)
  }
}

#' @name utils_transform
#' @export
image_vreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {


  # -- inner helper: reflect one Image object ---------------------------------
  .one <- function(im) {
    as_image(
      image_vreflect_cpp(image_data(im)),
      colormode = attr(im, "colormode") %||% "Color"
    )
  }

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
        .f = function(im) pliman::image_vreflect(im)
      )[.progress]

    } else {
      res <- lapply(img, .one)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    return(res)

  } else {
    img <- .one(img)
    if (isTRUE(plot)) plot(img)
    return(img)
  }
}


#' Resize an image
#'
#' @description
#' Resizes a 2D (grayscale) or 3D (RGB) image or a list of images. It supports
#' three reconstruction filters: Lanczos-3 (highest quality), Bilinear, and
#' Nearest-Neighbor.
#'
#'
#' @param img An \code{Image} object or a list of \code{Image} objects.
#' @param rel_size The relative size of the resized image (as a percentage of the
#'   original size). Defaults to 100. For example, \code{rel_size = 50} resizes a
#'   \code{1000 x 800} image to \code{500 x 400}.
#' @param width The target width (number of rows in the output matrix) in pixels.
#'   If missing, it is computed automatically to preserve the aspect ratio based
#'   on \code{height} or \code{rel_size}.
#' @param height The target height (number of columns in the output matrix) in
#'   pixels. If missing, it is computed automatically to preserve the aspect ratio.
#' @param filter The interpolation filter to use. Options are:
#'   * \code{"lanczos"}: Lanczos-3 windowed sinc resampling (highest quality).
#'   * \code{"bilinear"}: Bilinear interpolation (fast, smooth preview).
#'   * \code{"nearest"}: Nearest-neighbor interpolation (use for binary masks or label images).
#' @param parallel Logical value indicating whether to process lists of images in parallel.
#' @param workers The maximum number of parallel processes to run when \code{parallel = TRUE}.
#' @param verbose Logical value. If \code{TRUE} (default), displays step notifications in the console.
#' @param plot Logical value. If \code{TRUE}, plots the resized image(s). Defaults to \code{FALSE}.
#'
#' @return A modified \code{Image} object or a list of modified \code{Image} objects.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#'
#' @examples
#' if (interactive()) {
#' library(pliman)
#'
#' # Load example image
#' img <- image_pliman("sev_leaf.jpg")
#'
#' # Resize to 50% using high-quality Lanczos filter
#' img_half <- image_resize(img, rel_size = 50, filter = "lanczos", plot = TRUE)
#'
#' # Resize specifying custom width and height
#' img_custom <- image_resize(img, width = 400, height = 300, filter = "bilinear")
#'
#' # Resize a list of images
#' img_list <- list(img, img)
#' resized_list <- image_resize(img_list, rel_size = 25, parallel = FALSE)
#' }
image_resize <- function(img,
                         rel_size = 100,
                         width,
                         height,
                         filter = c("lanczos", "bilinear", "nearest"),
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE) {


  filter <- match.arg(filter)
  filter_code <- switch(filter, lanczos = 0L, bilinear = 1L, nearest = 2L)

  # Internal helper: resize a single Image using the C++ backend
  .resize_one <- function(im, target_w, target_h) {
    cm <- attr(im, "colormode") %||% "Color"
    arr <- image_resize_cpp(image_data(im),
                            as.integer(target_w),
                            as.integer(target_h),
                            filter_code)
    as_image(arr, colormode = cm)
  }

  # Compute target dimensions from a single image
  .compute_dims <- function(im, w_arg, h_arg, has_w, has_h) {
    d <- dim(im)[1:2]
    if (has_w && has_h) {
      return(c(as.integer(w_arg), as.integer(h_arg)))
    }
    if (has_w && !has_h) {
      tw <- as.integer(w_arg)
      th <- as.integer(round(tw * d[2] / d[1]))
      return(c(tw, th))
    }
    if (!has_w && has_h) {
      th <- as.integer(h_arg)
      tw <- as.integer(round(th * d[1] / d[2]))
      return(c(tw, th))
    }
    # Neither explicit: use rel_size
    tw <- as.integer(round(d[1] * rel_size / 100))
    th <- as.integer(round(d[2] * rel_size / 100))
    c(tw, th)
  }

  has_w <- !missing(width)
  has_h <- !missing(height)
  w_val <- if (has_w) width else NULL
  h_val <- if (has_h) height else NULL

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All images must be of class {.cls Image}.")
    }

    resize_single <- function(im) {
      dims <- .compute_dims(im, w_val, h_val, has_w, has_h)
      .resize_one(im, dims[1], dims[2])
    }

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
        .f = resize_single
      )[.progress]

    } else {
      res <- lapply(img, resize_single)
    }
    if(plot){
      image_combine(res)
    }
    return(res)

  } else {
    dims <- .compute_dims(img, w_val, h_val, has_w, has_h)
    img <- .resize_one(img, dims[1], dims[2])
    if (isTRUE(plot)) plot(img)
    return(img)
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

    return(res)

  } else {
    img <- img[, -c(1:top), ]
    img <- img[, -c((dim(img)[2] - bottom + 1):dim(img)[2]), ]
    img <- img[-c((dim(img)[1] - right + 1):dim(img)[1]), , ]
    img <- img[-c(1:left), , ]

    if (isTRUE(plot)) {
      plot(img)
    }

    return(img)
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


  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    # funcao auxiliar para aplicar skeletonization
    skel_fun <- function(im) {
      if (attr(im, "colormode") != "Grayscale") {
        im <- help_binary(im, ..., resize = FALSE)
      }

      s <- matrix(1, nrow(im), ncol(im))
      skel <- matrix(0, nrow(im), ncol(im))
      k <- if (is.null(kern)) make_brush(2, shape = "diamond") else kern

      while (max(s) == 1) {
        opened <- image_opening(im, k)
        s <- im - opened
        skel <- skel | s
        im <- image_erode(im, k)
      }

      as_image(skel)
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

    return(res)

  } else {
    if (attr(img, "colormode") != "Grayscale") {
      img <- help_binary(img, ..., resize = FALSE)
    }

    s <- matrix(1, nrow(img), ncol(img))
    skel <- matrix(0, nrow(img), ncol(img))
    kern <- if (is.null(kern)) make_brush(2, shape = "diamond") else kern

    while (max(s) == 1) {
      opened <- image_opening(img, kern)
      s <- img - opened
      skel <- skel | s
      img <- image_erode(img, kern)
    }

    img <- as_image(skel)
    if (isTRUE(plot)) plot(img)
    return(img)
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


  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    thin_fun <- function(im) {
      if (attr(im, "colormode") != "Grayscale") {
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

      as_image(im)
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

    return(return(res))

  } else {
    if (attr(img, "colormode") != "Grayscale") {
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

    img <- as_image(img)
    if (isTRUE(plot)) plot(img)
    return(img)
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
#'    two-subiteration algorithms. Commun. ACM 32(3): 359-373.
#'    \doi{10.1145/62065.62074}
#' @return If \code{img} is a single binary image, the function returns the
#'   thinned binary image. If \code{img} is a list of binary images, the
#'   function returns a list containing the thinned binary images.
#' @export
#'
#' @examples
#' if (interactive()) {
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


  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    thin_fun <- function(im) {
      if (attr(im, "colormode") != "Grayscale") {
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

    return(return(res))

  } else {
    if (attr(img, "colormode") != "Grayscale") {
      img <- help_binary(img, ..., resize = FALSE)
    }

    thin <- helper_guo_hall(img)

    if (isTRUE(plot)) {
      plot(thin)
    }

    return(thin)
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

  size <- as.integer(size)
  if (size < 1L) {
    cli::cli_abort("Using {.arg size} < 1 may crash. Use 1 or more.")
  }

  filter_fun <- function(im) {
    if (inherits(im, "Image")) {
      arr <- as.numeric(image_data(im))
    } else {
      arr <- as.array(im)
    }

    dm  <- dim(arr)
    nr  <- dm[1L]
    nc  <- dm[2L]
    nch <- if (length(dm) == 3L) dm[3L] else 1L

    # Garante array 3D para o C++
    if (length(dm) == 2L) dim(arr) <- c(nr, nc, 1L)

    if (storage.mode(arr) == "logical") {
      res <- median_filter_binary_cpp(arr, nr, nc, nch, size)
    } else {
      res <- median_filter_cpp(arr, nr, nc, nch, size)
    }
    dim(res) <- dm

    if (inherits(im, "Image")) {
      return(as_image(res, colormode = im@colormode))
    }
    res
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
        left  = cli::col_blue("Median filtering {.val {length(img)}} images"),
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

    return(return(res))

  } else {
    res <- filter_fun(img)
    if (isTRUE(plot)) {
      plot(res)
    }
    return(res |> as_image(storage = typeof(img)))
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

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    blur_fun <- function(im) {
      help_gblur(im, sigma)
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

    return(res)

  } else {
    img <- help_gblur(img, sigma)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_contrast <- function(img,
                           nx = 8,
                           ny = 8,
                           clip_limit = 3,
                           bins = 256,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE) {

  contrast_fun <- function(im) {
    cpp_clahe(im, nx = nx, ny = ny, clip_limit = clip_limit, nbins = bins)
  }

  if (is.list(img)) {
    if (inherits(img, c("binary_list", "segment_list", "index_list",
                        "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1]])
    }

    if (!all(sapply(img, is_image))) {
      cli::cli_abort("All elements in the list must be of class {.cls image} or {.cls Image}.")
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

    return(return(res))

  } else {
    res <- contrast_fun(img)
    if (isTRUE(plot)) {
      plot(res)
    }
    return(res)
  }
}

#' Create an `image` object of a given color
#'
#' image_create() can be used to create an `image` object with a desired color and size.
#'
#' @param color either a color name (as listed by [grDevices::colors()]), or a hexadecimal
#'   string of the form `"#rrggbb"`.
#' @param width,heigth The width and heigth of the image in pixel units.
#' @param plot Plots the image after creating it? Defaults to `FALSE`.
#'
#' @return An object of class `image`.
#' @export
#'
#' @examples
#' if (interactive()) {
#' image_create("red")
#' image_create("#009E73", width = 300, heigth = 100)
#' }

image_create <- function(color,
                         width = 200,
                         heigth = 200,
                         plot = FALSE){

  width <- as.integer(width)
  heigth <- as.integer(heigth)
  rgb <- col2rgb(color) / 255
  r <- rep(rgb[1], width*heigth)
  g <- rep(rgb[2], width*heigth)
  b <- rep(rgb[3], width*heigth)
  img <- as_image(as_image(c(r, g, b),
                           dim = c(width, heigth, 3),
                           colormode = "color"))
  if(isTRUE(plot)){
    plot(img)
  }
  return(img)
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
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param max_size Maximum size of objects to keep. Larger objects are removed. Default is `NULL`.
#' @param remove_reflection Logical or numeric. If `TRUE` or numeric, removes internal reflection holes from binary objects. If numeric, specifies the neck distance threshold. Default is `FALSE`.
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
#' @param filter_order A character vector indicating the order in which morphological
#'   operations are applied. Defaults to `c("erode", "dilate", "opening", "closing", "filter", "fill_hull")`.
#' @param invert Inverts the binary image, if desired.
#' @param return_exact Logical. If `TRUE`, exact morphological reconstruction via
#'   watershed-partitioned seed matching is performed when `opening > 0`.
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
#' Trans Syst Man Cybern SMC-9(1): 62-66. \doi{10.1109/tsmc.1979.4310076}
#'
#' Shafait, F., D. Keysers, and T.M. Breuel. 2008. Efficient implementation of
#' local adaptive thresholding techniques using integral images. Document
#' Recognition and Retrieval XV. SPIE. p. 317-322 \doi{10.1117/12.767755}
#'
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing binary images. The length will depend on the number
#'   of indexes used.
#' @importFrom utils read.csv
#' @examples
#' if (interactive()) {
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
                         return_class = "image",
                         threshold = c("Otsu", "adaptive"),
                         k = 0.15,
                         windowsize = NULL,
                         has_white_bg = FALSE,
                         resize = FALSE,
                         fill_hull = FALSE,
                         max_size = NULL,
                         remove_reflection = FALSE,
                         erode = FALSE,
                         dilate = FALSE,
                         opening = FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                         invert = FALSE,
                         return_exact = FALSE,
                         plot = TRUE,
                         nrow = NULL,
                         ncol = NULL,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE) {

  # remove_reflection: fill internal holes in binary (reflections appear as holes inside grains).
  # Uses min_neck_dist to avoid filling inter-grain spaces.
  # TRUE -> min_neck_dist = 3 (default); numeric N -> min_neck_dist = N.
  rr_neck_dist <- if (!isFALSE(remove_reflection)) {
    fill_hull <- TRUE
    if (is.numeric(remove_reflection) && remove_reflection > 0) as.double(remove_reflection) else 3.0
  } else {
    0.0
  }

  check_filter_order(filter_order, verbose, erode, dilate, opening, closing, filter, fill_hull)
  threshold <- threshold[[1]]
  if (is.character(threshold)) {
    if (tolower(threshold) == "otsu") {
      threshold <- "Otsu"
    } else if (tolower(threshold) == "adaptive") {
      threshold <- "adaptive"
    }
  }

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
      binary_mat <- threshold_adaptive(image_data(imgs), k, windowsize)
    } else {
      if(threshold == "Otsu"){
        threshold_val <- help_otsu(image_data(imgs))
      } else if(is.numeric(threshold)) {
        threshold_val <- threshold
      } else {
        t_data <- image_data(imgs, type = "numeric")
        pixels <- terra::rast(t(t_data))
        terra::plot(pixels, col = custom_palette(n = 100), axes = FALSE, asp = NA)
        threshold_val <- readline("Selected threshold: ")
      }
      binary_mat <- imgs < threshold_val
    }

    if (isTRUE(invert)) {
      binary_mat <- !binary_mat
    }

    binary_mat[is.na(binary_mat)] <- FALSE

    do_exact_opening <- isTRUE(return_exact) && ((is.numeric(opening) && opening > 0) || isTRUE(opening))

    raw_bin <- unclass(image_data(binary_mat))
    if (!is.logical(raw_bin)) raw_bin <- raw_bin != 0
    dim(raw_bin) <- dim(binary_mat)[1:2]
    res_mat <- help_binary_filters_cpp(
      img_sexp = raw_bin,
      erode = if (is.numeric(erode)) as.integer(erode) else 0L,
      dilate = if (is.numeric(dilate)) as.integer(dilate) else 0L,
      opening = if (!do_exact_opening && is.numeric(opening)) as.integer(opening) else 0L,
      closing = if (is.numeric(closing)) as.integer(closing) else 0L,
      filter = if (is.numeric(filter)) as.integer(filter) else 0L,
      fill_hull = isTRUE(fill_hull),
      filter_order = filter_order,
      max_size = if (is.null(max_size)) -1.0 else as.numeric(max_size),
      min_neck_dist = rr_neck_dist
    )

    out_img <- as_image(res_mat, storage = "raw")

    if (do_exact_opening) {
      sz <- if (is.numeric(opening)) opening else NULL
      out_img <- image_opening(out_img, size = sz, return_exact = TRUE, verbose = FALSE, plot = FALSE)
    }

    # Return as a logical image S3 object
    return(out_img)
  }

  process_image <- function(im) {
    is_gray <- length(dim(im)) == 2 || (length(dim(im)) == 3 && dim(im)[3] == 1)
    if (is_gray) {
      imgs <- list(gray = bin_img(im))
    } else {
      imgs <- lapply(
        image_index(im,
                    index = index,
                    r = r,
                    g = g,
                    b = b,
                    re = re,
                    nir = nir,
                    return_class = return_class,
                    resize = resize,
                    has_white_bg = has_white_bg,
                    plot = FALSE,
                    nrow = nrow,
                    ncol = ncol,
                    verbose = verbose),
        bin_img
      )
    }
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
    return(structure(res, class = "binary_list"))
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
    return(imgs)
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
#' @param img An `image` object. Multispectral mosaics can be converted to an
#'   `image` object using `mosaic_as_ebimage()`.
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
#'   computed. If `"ebimage"` (default) returns a list of `image` objects, where
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
#' Digital Image Analysis. Crop Science 43(3): 943-951.
#' \doi{10.2135/cropsci2003.9430}
#'
#' Bannari, A., D. Morin, F. Bonn, and A.R. Huete. 1995. A review of vegetation
#' indices. Remote Sensing Reviews 13(1-2): 95-120.
#' \doi{10.1080/02757259509532298}
#'
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing Grayscale images. The length will depend on the
#'   number of indexes used.
#' @examples
#' if (interactive()) {
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
                        return_class = c("image", "terra"),
                        resize = FALSE,
                        has_white_bg = FALSE,
                        plot = TRUE,
                        nrow = NULL,
                        ncol = NULL,
                        max_pixels = 500000,
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        ...){

  return_classopt <- c("terra", "image", "ebimage")
  return_classopt <- return_classopt[pmatch(return_class[1], return_classopt)]
  if (is.na(return_classopt)) {
    return_classopt <- "image"
  }
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
    return(structure(res, class = "index_list"))
  } else{
    if(resize != FALSE){
      img <- image_resize(img, resize)
    }
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
    hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])
    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    } else {
      if(index[[1]] %in% c("RGB", "NRGB", "MULTISPECTRAL", "all")){
        index <- switch(index,
          RGB = c("R", "G", "B"),
          NRGB = c("NR", "NG", "NB"),
          MULTISPECTRAL = c("NDVI", "PSRI", "GNDVI", "RVI", "NDRE", "TVI", "CVI", "EVI", "CIG", "CIRE", "DVI", "NDWI"),
          all = ind$Index
        )
      } else {
        if(length(index) > 1){
          index <- index
        } else {
          index <- strsplit(index, "\\s*(,)\\s*")[[1]]
        }
      }
    }

    img_num <- NULL
    get_img_num <- function() {
      if (is.null(img_num)) {
        img_num <<- image_data(img, type = "numeric")
      }
      img_num
    }

    if(any(index %in% hsb_ind)){
      inum <- get_img_num()
      R <- inum[,,r]
      G <- inum[,,g]
      B <- inum[,,b]
      hsb <- rgb_to_hsb(data.frame(R = c(R), G = c(G), B = c(B)))
      h <- matrix(hsb$h, nrow = nrow(img), ncol = ncol(img))
      s <- matrix(hsb$s, nrow = nrow(img), ncol = ncol(img))
      b <- matrix(hsb$b, nrow = nrow(img), ncol = ncol(img))
    }
    if(any(index %in% nir_ind)){
      if(dim(img)[3] < max(re, nir)){
        cli::cli_abort("Near-Infrared and RedEdge bands are not available in the provided image.")
      }
    }
    if(dim(img)[3] < 3){
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
        thresh <- if (is.raw(as.vector(img))) as.raw(0xff) else 1
        dat <- image_data(img)
        white <- dat[,,r] == thresh & dat[,,g] == thresh & dat[,,b] == thresh
        if (is.raw(dat)) {
          dat <- image_data(img, type = "normalized")
        }
        dat[,,r][white] <- NA
        dat[,,g][white] <- NA
        dat[,,b][white] <- NA
        img <- as_image(dat, storage = "double")
      }

      cpp_res <- try(compute_single_index_cpp(img, indx, r, g, b, re, nir), silent = TRUE)
      if (!is.null(cpp_res) && !inherits(cpp_res, "try-error")) {
        imgs[[i]] <- cpp_res
      } else if (indx %in% ind$Index) {
        imgs[[i]] <- as_image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index) == indx]))), storage = "double")
      } else {
        imgs[[i]] <- as_image(eval(parse(text = as.character(indx))), storage = "double")
      }
    }
    names(imgs) <- index
    class(imgs) <- "image_index"
    if(plot == TRUE){
      plot_index(imgs, nrow = nrow, ncol = ncol, max_pixels = max_pixels, ...)
    }
    if(return_classopt %in% c("image", "ebimage")){
      return(imgs)
    } else{
      terras <-
        terra::rast(
          lapply(1:length(imgs), function(i){
            u <- image_data(imgs[[i]])
            class(u) <- NULL
            dim(u) <- dim(imgs[[i]])[1:2]
            if (is.raw(u)) u <- as.double(u) / 255
            terra::rast(t(u))
          })
        )
      names(terras) <- names(imgs)
      return(terras)
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
#' if (interactive()) {
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
#' if (interactive()) {
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
    return(structure(res, class = "segment_list"))
  } else{
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
    hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])
    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    }else{
      RE <- try(as.numeric(image_data(img)[,,re]), TRUE)
      NIR <- try(as.numeric(image_data(img)[,,nir]), TRUE)
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
      col_bg_norm <- c(1, 1, 1)
    } else {
      if (is.character(col_background)) {
        col_bg_norm <- col2rgb(col_background) / 255
      } else if (is.numeric(col_background) && max(col_background, na.rm = TRUE) > 1.0) {
        col_bg_norm <- col_background / 255
      } else {
        col_bg_norm <- col_background
      }
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
      ID <- which(image_data(img2) == FALSE)
      if(!na_background){
        if (is.raw(image_data(imgmask))) {
          bg_col <- as.raw(round(col_bg_norm * 255))
          bg_extra <- as.raw(255)
        } else {
          bg_col <- col_bg_norm
          bg_extra <- 1
        }
        imgmask[,,1][ID] <- bg_col[1]
        imgmask[,,2][ID] <- bg_col[2]
        imgmask[,,3][ID] <- bg_col[3]
        if(dim(img)[[3]] > 3){
          imgmask[,,4][ID] <- bg_extra
          imgmask[,,5][ID] <- bg_extra
        }
      } else{
        if (is.raw(image_data(imgmask))) {
          imgmask <- to_numeric(imgmask)
        }
        imgmask[,,1][ID] <- NA
        imgmask[,,2][ID] <- NA
        imgmask[,,3][ID] <- NA
        if(dim(img)[[3]] > 3){
          imgmask[,,4][ID] <- NA
          imgmask[,,5][ID] <- NA
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
        tmps <- as_image(tmps, colormode = "Color")
        plot(tmps)
        if(verbose == TRUE){
          dim <- image_dimension(imgs[[i]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
        }
      }
    }
    if(length(imgs) == 1){
      return(imgs[[1]])
    } else{
      return(structure(imgs, class = "img_segment"))
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
    return(list(results = results,
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
      return(list(results = results,
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
      return(list(results = pixels,
                     images = imgs))
    }
  }
}



#' Image segmentation using k-means clustering
#'
#' Segments image objects using clustering by the k-means clustering algorithm
#' @inheritParams image_segment
#' @inheritParams analyze_objects
#' @param img An `image` object.
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
#'   K-means clustering algorithm. Applied Statistics, 28, 100-108.
#'   \doi{10.2307/2346830}
#'
#' @examples
#' if (interactive()) {
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

  imm <- img[, , bands, drop = FALSE]
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
    LIST <- lapply(LIST, image_fill_hull)
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
    LIST <- lapply(LIST, image_filter, size = filter)
  }
  mask <- LIST[[1]]
  if(isFALSE(invert)){
    id <-  which(mask != 0)
  } else{
    id <- which(mask != 1)
  }
  im2 <- img
  im2[, , 1][id] <- 1
  im2[, , 2][id] <- 1
  im2[, , 3][id] <- 1
  if(isTRUE(plot)){
    if(nclasses == 2){
      plot(im2)
    } else{
      suppressWarnings(image(m, useRaster = TRUE))
    }
  }
  return(list(img = im2,
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
#' @param img An `image` object.
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
    mask <- image_fill_hull(image_bwlabel(imF))
    # return(mask)
    if(type[1] == "select"){
      id <- mask != 1
    } else{
      id <- mask == 1
    }
    img[, , 1][id] <- 1
    img[, , 2][id] <- 1
    img[, , 3][id] <- 1

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
    return(list(img = image, mask = as_image(mask)))
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
#' if (interactive()) {
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
    return(structure(res, class = "img_mat_list"))
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
    return(mat)
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
#' * `palette_list` A list with `npal` color palettes of class `image`.
#' * `joint` An object of class `image` with the color palettes
#' * `proportions` The proportion of the entire image corresponding to each color in the palette
#' * `rgbs` The average RGB value for each palette
#' @name palettes
#' @export
#' @importFrom stats na.omit
#' @examples
#' if (interactive()) {
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
                           filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                           plot = TRUE,
                           save_image = FALSE,
                           prefix = "proc_",
                           dir_original = NULL,
                           dir_processed = NULL,
                           return_pal = FALSE,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE) {

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
      mask <- image_binary(img, index = "B-R", opening = 5, filter_order = filter_order, plot = FALSE, verbose = FALSE)[[1]]
      ID <- which(mask == FALSE)
      img[,,1][ID] <- NA
      img[,,2][ID] <- NA
      img[,,3][ID] <- NA
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
        pal_list[[paste0("pal_", i)]] <- rgb_image(R, G, B)
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
      im2 <- as_image(ARR)
      im2 <- image_resize(im2, height = ncol(img), width = nc * 0.1)
      im2[1:nrow(im2), 1:1, ] <- 0
      im2[1:1, 1:ncol(im2), ] <- 0
      im2[1:nrow(im2), ncol(im2):(ncol(im2)-1), ] <- 0
      im2[nrow(im2):(nrow(im2)-1), 1:ncol(im2), ] <- 0
      im2 <- image_combine(img, im2, along = 1)
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
             width = dim(im2)[1],
             height = dim(im2)[2])
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
        if (verbose) cli::cli_progress_update()
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
#' @param img An `image` object.
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
#' @return An `image` object
#' @export
#'
#' @examples
#' if (interactive()) {
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
    left_img <- as_image(img[1:sample_left,,])
    left_img <- image_resize(left_img, width = left, height = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(left_img)
      for (i in 1:nc[1]) {
        left_img[i,,] <- left_img[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      left_img <- image_filter(left_img, size = filter)
    }
    img <- image_combine(left_img, img, along = 1)
  }
  if(!is.null(top)){
    top_img <- as_image(img[,1:sample_top,])
    top_img <- image_resize(top_img, width = dim(img)[1], height = top)
    if(isTRUE(random)){
      nc <- dim(top_img)
      for (i in 1:nc[2]) {
        top_img[,i,] <- top_img[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      top_img <- image_filter(top_img, size = filter)
    }
    img <- image_combine(top_img, img, along = 2)
  }
  if(!is.null(right)){
    dimx <- dim(img)[1]
    right_img <- as_image(img[(dimx-sample_right):dimx,,])
    right_img <- image_resize(right_img, width = right, height = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(right_img)
      for (i in 1:nc[1]) {
        right_img[i,,] <- right_img[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      right_img <- image_filter(right_img, size = filter)
    }
    img <- image_combine(img, right_img, along = 1)
  }
  if(!is.null(bottom)){
    dimy <- dim(img)[2]
    bot_img <- as_image(img[,(dimy-sample_bottom):dimy,])
    bot_img <- image_resize(bot_img, width = dim(img)[1], height = bottom)
    if(isTRUE(random)){
      nc <- dim(bot_img)
      for (i in 1:nc[2]) {
        bot_img[,i,] <- bot_img[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      bot_img <- image_filter(bot_img, size = filter)
    }
    img <- image_combine(img, bot_img, along = 2)
  }
  if(isTRUE(plot)){
    plot(img)
  }
  return(img)
}


#' Squares an image
#'
#' Converts a rectangular image into a square image by expanding the
#' rows/columns using [image_expand()].
#'
#' @inheritParams image_expand
#'
#' @return The modified `image` object.
#' @param ... Further arguments passed on to [image_expand()].
#' @export
#'
#' @examples
#' if (interactive()) {
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
  return(img)
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
    return(pix)
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
#' @param object An `image` object, an object computed with `analyze_objects()`
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
#' if (interactive()) {
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
  return(data.frame(hsb))
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
  return(data.frame(srgb))
}


#' @export
#' @name utils_colorspace
rgb_to_lab <- function(object){
  object <- rgb_to_srgb(object)
  srgb <- data.frame(r = object[, 1],
                     g = object[, 2],
                     b = object[, 3])
  lab <- convertColor(srgb, from = "sRGB", to = "Lab")
  return(lab)
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
                         na_background = FALSE,
                         has_white_bg = FALSE,
                         fill_hull = FALSE,
                         opening = FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         dilate = FALSE,
                         erode = FALSE,
                         invert = FALSE,
                         filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull")){
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
                      invert = invert,
                      filter_order = filter_order)
  ID <- which(image_data(img2) == FALSE)
  is_raw_img <- is.raw(image_data(img))
  if(!na_background && !is.null(col_background)){
    if (is_raw_img) {
      bg_val <- if (is.character(col_background)) as.raw(col2rgb(col_background)) else as.raw(round(col_background * 255))
    } else {
      bg_val <- if (is.character(col_background)) col2rgb(col_background)/255 else col_background/255
    }
    img[,,r][ID] <- bg_val[1]
    img[,,g][ID] <- bg_val[2]
    img[,,b][ID] <- bg_val[3]
  } else if (isTRUE(na_background)) {
    if (is_raw_img) {
      img <- to_numeric(img)
    }
    img[,,r][ID] <- NA
    img[,,g][ID] <- NA
    img[,,b][ID] <- NA
    if(dim(img)[3] > 3){
      img[,,4][ID] <- NA
    }
    if(dim(img)[3] > 4){
      img[,,5][ID] <- NA
    }
  } else {
    bg_one <- if (is_raw_img) as.raw(255) else 1
    if(dim(img)[3] == 3){
      img[,,r][ID] <- bg_one
      img[,,g][ID] <- bg_one
      img[,,b][ID] <- bg_one
    } else if(dim(img)[3] == 4){
      img[,,r][ID] <- bg_one
      img[,,g][ID] <- bg_one
      img[,,b][ID] <- bg_one
      img[,,re][ID] <- bg_one
    } else{
      img[,,r][ID] <- bg_one
      img[,,g][ID] <- bg_one
      img[,,b][ID] <- bg_one
      img[,,re][ID] <- bg_one
      img[,,nir][ID] <- bg_one
    }
  }
  return(img)
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
                        invert = FALSE,
                        filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                        return_exact = FALSE){
  check_filter_order(filter_order, verbose = FALSE, erode, dilate, opening, closing, filter, fill_hull)
  threshold <- threshold[[1]]

  bin_img <- function(imgs,
                      invert,
                      fill_hull,
                      threshold,
                      erode,
                      dilate,
                      opening,
                      closing,
                      filter,
                      filter_order){
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
      binary_mat <- threshold_adaptive(image_data(imgs), k, windowsize)
    }
    if(threshold != "adaptive"){
      if(threshold == "Otsu"){
        threshold <- help_otsu(image_data(imgs))
      } else{
        if(is.numeric(threshold)){
          threshold <- threshold
        } else{
          pixels <- terra::rast(t(image_data(imgs)))
          terra::plot(pixels, col = custom_palette(n = 100),  axes = FALSE, asp = NA)
          threshold <- readline("Selected threshold: ")
        }
      }
      binary_mat <- imgs < threshold
    }

    if(invert == TRUE){
      binary_mat <- !binary_mat
    }

    binary_mat[is.na(binary_mat)] <- FALSE

    do_exact_opening <- isTRUE(return_exact) && ((is.numeric(opening) && opening > 0) || isTRUE(opening))

    raw_bin <- unclass(image_data(binary_mat))
    if (!is.logical(raw_bin)) raw_bin <- raw_bin != 0
    dim(raw_bin) <- dim(binary_mat)[1:2]
    res_mat <- help_binary_filters_cpp(
      img_sexp = raw_bin,
      erode = if (is.numeric(erode)) as.integer(erode) else 0L,
      dilate = if (is.numeric(dilate)) as.integer(dilate) else 0L,
      opening = if (!do_exact_opening && is.numeric(opening)) as.integer(opening) else 0L,
      closing = if (is.numeric(closing)) as.integer(closing) else 0L,
      filter = if (is.numeric(filter)) as.integer(filter) else 0L,
      fill_hull = isTRUE(fill_hull),
      filter_order = filter_order
    )

    out_img <- as_image(res_mat, storage = "raw")
    if (do_exact_opening) {
      sz <- if (is.numeric(opening)) opening else NULL
      out_img <- image_opening(out_img, size = sz, return_exact = TRUE, verbose = FALSE, plot = FALSE)
    }

    return(out_img)
  }

  is_gray <- length(dim(img)) == 2 || (length(dim(img)) == 3 && dim(img)[3] == 1)
  if (is_gray) {
    gray_img <- img
  } else {
    gray_img <- help_imageindex(img, index, r, g, b, re, nir, resize, has_white_bg)
  }
  bin_img <- bin_img(gray_img,
                     invert,
                     fill_hull,
                     threshold,
                     erode,
                     dilate,
                     opening,
                     closing,
                     filter,
                     filter_order)
  return(bin_img)
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

  if(isTRUE(has_white_bg)){
    # Set pure-white pixels to NA in RGB bands (raw: 255; double: 1)
    thresh <- if (is.raw(as.vector(img))) as.raw(0xff) else 1
    dat <- image_data(img)
    white <- dat[,,r] == thresh & dat[,,g] == thresh & dat[,,b] == thresh
    if (is.raw(dat)) {
      dat <- image_data(img, type = "normalized")
    }
    dat[,,r][white] <- NA
    dat[,,g][white] <- NA
    dat[,,b][white] <- NA
    img <- as_image(dat, storage = "double")
  }

  # Try the fast C++ path first
  cpp_res <- try(compute_single_index_cpp(img, index, r, g, b, re, nir), silent = TRUE)
  if (!is.null(cpp_res) && !inherits(cpp_res, "try-error")) {
    return(return(cpp_res))
  }

  # Fallback: R-based eval(parse) for custom expressions or unknown indices
  ind <- read.csv(file = system.file("indexes.csv", package = "pliman", mustWork = TRUE),
                  header = TRUE, sep = ";")
  nir_ind <- as.character(ind$Index[ind$Band %in% c("MULTI")])
  hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])

  R <- try(as.numeric(image_data(img)[,,r]), TRUE)
  G <- try(as.numeric(image_data(img)[,,g]), TRUE)
  B <- try(as.numeric(image_data(img)[,,b]), TRUE)
  RE  <- try(as.numeric(image_data(img)[,,re]),  TRUE)
  NIR <- try(as.numeric(image_data(img)[,,nir]), TRUE)

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

  if(index %in% ind$Index){
    res_eval <- eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index])))
  } else{
    res_eval <- eval(parse(text = as.character(index)))
  }
  dim(res_eval) <- dim(img)[1:2]
  img_gray <- as_image(res_eval)
  return(img_gray)
}




#' Prepare images to analyze_objects_shp()
#'
#' It is a simple wrapper around [image_align()] and [image_crop()]. In this case, only the option `viewer = "base"` is used. To use `viewer = "mapview"`, please, use such functions separately.
#'
#' @param img A `image` object
#' @inheritParams image_align
#'
#' @return An aligned and cropped `image` object.
#' @export
#'
#' @examples
#' if (interactive()) {
#' img <- image_pliman("flax_leaves.jpg")
#' prepare_to_shp(img)
#' }
prepare_to_shp <- function(img,
                           align = "vertical"){

  aligned <- image_align(img, viewer = "base")
  cropped <- image_crop(aligned, viewer = "base", plot = TRUE)
  return(cropped)
}

#' Add Alpha Layer to an RGB Image
#'
#' This function adds an alpha (transparency) layer to an RGB image.
#' The alpha layer can be specified as a single numeric value for uniform transparency
#' or as a matrix/array matching the dimensions of the image for varying transparency.
#'
#' @param img An RGB image of class `image` ..
#' @param mask A numeric value or matrix/array specifying the alpha layer:
#'     * If `mask` is a single numeric value, it sets a uniform transparency level (0 for fully transparent, 1 for fully opaque).
#'     * If `mask` is a matrix or array, it must have the same dimensions as the image channels, allowing for varying transparency.
#'
#' @return An `image` object with an added alpha layer, maintaining the RGBA format.
#'
#' @examples
#' if (interactive()) {
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
  if (attr(img, "colormode") != "Color") {
    cli::cli_abort("Input image must be in RGB format.")
  }

  m_dim <- if (is.array(mask) || is.matrix(mask)) dim(mask)[1:2] else NULL
  if (is.numeric(mask) && length(mask) == 1) {
    mask_bin <- array(mask > 0, dim = dim(img)[1:2])
  } else if (!is.null(m_dim) && all(m_dim == dim(img)[1:2])) {
    mask_bin <- mask > 0
    if (length(dim(mask_bin)) == 3) {
      mask_bin <- mask_bin[,,1]
    }
  } else {
    cli::cli_abort("Mask must be either a single numeric value or a matrix with the same dimensions as the image channels.")
  }

  a <- image_data(img)
  is_raw_img <- is.raw(a)

  if (is_raw_img) {
    alpha_layer <- as.raw(ifelse(mask_bin, 255, 0))
    arr <- array(c(a[,,1], a[,,2], a[,,3], alpha_layer), dim = c(dim(img)[1:2], 4))
    as_image(arr, colormode = "Color")
  } else {
    alpha_layer <- ifelse(mask_bin, 1.0, 0.0)
    arr <- array(c(a[,,1], a[,,2], a[,,3], alpha_layer), dim = c(dim(img)[1:2], 4))
    as_image(arr, colormode = "Color")
  }
}

#' Label Connected Components in a Binary Image
#'
#' This function labels connected components in a binary image while allowing
#' for a specified maximum gap between pixels to still be considered part of
#' the same object.
#'
#' @param img A binary image matrix where `1` represents foreground pixels and
#'   `0` represents background pixels.
#' @param max_gap An integer specifying the maximum allowable gap (in pixels)
#'   between connected components to be considered as part of the same object.
#'   Default is `1`.
#'
#' @return An object of class `image`, where each
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
#' image_label(img, max_gap = 1) |> plot()
#' image_label(img, max_gap = 2)
#' image_label(img, max_gap = 3)
#' }
image_label <- function(img, max_gap = 0){
  help_label(img, max_gap = max_gap) |> as_image()
}



#' Smooth Contour Line Detection
#'
#'
#' @param img An `image` object.
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
  mat <- as.numeric(image_data(ind)) * 255
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
#' @param img An `image` object.
#' @param index A character string with the index to be used. Defaults to `"GRAY"`.
#' @param s sigma, the Gaussian filter variance. Defaults to 5.
#' @param low_thr lower threshold value of the algorithm. Defaults to 10.
#' @param high_thr upper threshold value of the algorithm. Defaults to 20
#' @return a list with an `image` object with values 0 or 255, and the number of
#'   pixels which have value 255 (pixels_nonzero).
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' conts <- image_canny_edge(img, index = "B")
#' image_combine(img, conts$edges)
#' }
image_canny_edge <- function(img,
                             index = "GRAY",
                             s = 5,
                             low_thr = 10,
                             high_thr = 20) {
  ind <- image_index(img, index, plot = FALSE)[[1]]
  d <- dim(ind)
  w <- d[1L]
  h <- d[2L]
  mat <- as.integer(round(as.numeric(image_data(ind)) * 255))
  res <- canny_edge_detector(mat, w, h, s, low_thr, high_thr, TRUE)
  res$edges <- as_image(res$edges, dim = c(w, h), colormode = "Grayscale")
  return(res[1:2])
}

#' @title Line Segment Detection in an Image
#' @description Detects line segments in a digital image using the Line Segment
#'   Detector (LSD), a linear-time method that controls false detections and
#'   requires no parameter tuning. Based on Burns, Hanson, and Riseman's method
#'   with an a-contrario validation approach.
#'
#' @param img An `image` object.
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
  ind <- image_index(image_hreflect(image_transpose(img)), index, plot = FALSE)[[1]]
  x <- as.numeric(image_data(ind)) * 255
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

#' @title Standard ColorChecker Reference Chart
#'
#' @description
#' Provides the standard reference RGB values for the 24-patch ColorChecker
#' Classic (Macbeth) chart in sRGB D65 color space.
#'
#' @param chart The reference chart standard. Currently `"colorchecker24"` (default)
#'   or `"macbeth"`.
#'
#' @return A `data.frame` with 24 rows and columns:
#' \itemize{
#'   \item `id`: Patch identifier (1 to 24).
#'   \item `code`: Patch code (`No.001` to `No.024`).
#'   \item `name`: Descriptive name of the color patch.
#'   \item `R`, `G`, `B`: Standard reference sRGB values (0-255 scale).
#'   \item `ref_R`, `ref_G`, `ref_B`: Alias reference columns.
#' }
#'
#' @export
#' @examples
#' ref <- colorchecker_ref()
#' head(ref)
colorchecker_ref <- function(chart = c("colorchecker24", "macbeth")) {
  chart <- match.arg(chart)
  data.frame(
    id = 1:24,
    code = sprintf("No.%03d", 1:24),
    name = c(
      "Dark skin", "Light skin", "Blue sky", "Foliage", "Blue flower", "Bluish green",
      "Orange", "Purplish blue", "Moderate red", "Purple", "Yellow green", "Orange yellow",
      "Blue", "Green", "Red", "Yellow", "Magenta", "Cyan",
      "White (0.04*)", "Neutral 8 (0.02*)", "Neutral 6.5 (0.42*)", "Neutral 5 (0.68*)", "Neutral 3.5 (1.00*)", "Black (1.50*)"
    ),
    R = c(115, 204, 101,  89, 141, 132, 249,  80, 222,  91, 173, 255,  44,  74, 179, 250, 191,   6, 252, 230, 200, 143, 100,  50),
    G = c( 82, 161, 134, 109, 137, 228, 118,  91,  91,  63, 232, 164,  56, 148,  42, 226,  81, 142, 252, 230, 200, 143, 100,  50),
    B = c( 69, 141, 179,  61, 194, 208,  35, 182, 125, 123,  91,  26, 142,  81,  50,  21, 160, 172, 252, 230, 200, 142, 100,  50),
    ref_R = c(115, 204, 101,  89, 141, 132, 249,  80, 222,  91, 173, 255,  44,  74, 179, 250, 191,   6, 252, 230, 200, 143, 100,  50),
    ref_G = c( 82, 161, 134, 109, 137, 228, 118,  91,  91,  63, 232, 164,  56, 148,  42, 226,  81, 142, 252, 230, 200, 143, 100,  50),
    ref_B = c( 69, 141, 179,  61, 194, 208,  35, 182, 125, 123,  91,  26, 142,  81,  50,  21, 160, 172, 252, 230, 200, 142, 100,  50),
    stringsAsFactors = FALSE
  )
}

#' @title Order Quadrilateral Corners Topologically
#'
#' @description
#' Sorts 4 quadrilateral corner vertices into standard Top-Left (TL),
#' Top-Right (TR), Bottom-Right (BR), and Bottom-Left (BL) ordering.
#'
#' @param pts A 4x2 matrix, data frame with 2 columns, list with `TL, TR, BR, BL`,
#'   or an 8-element numeric vector.
#'
#' @return A named list with `$TL`, `$TR`, `$BR`, `$BL` coordinates (numeric pairs `c(x, y)`).
#' @export
order_quad_corners <- function(pts) {
  if (is.list(pts) && all(c("TL", "TR", "BR", "BL") %in% names(pts))) {
    return(pts)
  }
  if (is.vector(pts) && length(pts) == 8) {
    pts <- matrix(pts, ncol = 2, byrow = TRUE)
  }
  if (is.data.frame(pts)) {
    pts <- as.matrix(pts[, 1:2])
  }
  if (!is.matrix(pts) || nrow(pts) < 4) {
    cli::cli_abort("Corners must be a 4x2 matrix, a list with TL, TR, BR, BL, or an 8-element vector.")
  }

  pts <- pts[1:4, 1:2, drop = FALSE]

  ord_y <- order(pts[, 2])
  top_pts <- pts[ord_y[1:2], , drop = FALSE]
  bot_pts <- pts[ord_y[3:4], , drop = FALSE]

  tl <- as.numeric(unname(top_pts[which.min(top_pts[, 1]), ]))
  tr <- as.numeric(unname(top_pts[which.max(top_pts[, 1]), ]))
  bl <- as.numeric(unname(bot_pts[which.min(bot_pts[, 1]), ]))
  br <- as.numeric(unname(bot_pts[which.max(bot_pts[, 1]), ]))

  list(TL = tl, TR = tr, BR = br, BL = bl)
}


#' @title Extract ColorChecker Palettes using Bilinear Homography and Negative Buffering
#'
#' @description
#' Automatically detects a ColorChecker chart in an image (or uses user-specified
#' corners), maps each color patch using perspective-aware bilinear homography,
#' and extracts pure nucleus color statistics using a dynamic negative buffer.
#'
#' @details
#' The extraction employs a bilinear interpolation mapping:
#' \deqn{P(u, v) = (1 - u)(1 - v)\mathbf{TL} + u(1 - v)\mathbf{TR} + uv\mathbf{BR} + (1 - u)v\mathbf{BL}}
#' Outer card margins (`margin_x`, `margin_y`) exclude the surrounding plastic/cardboard
#' frame. A dynamic negative buffer (`buffer`, default 45%) shrinks each nominal cell
#' inwards towards its center, ensuring pixel sampling is strictly restricted to the
#' homogeneous nucleus of the color chip, immune to edge glare, shadows, and plastic borders.
#'
#' @param img An `image` object (S3 pliman image or array).
#' @param reference Reference chart standard for comparison. Default `"colorchecker24"`.
#'   Can be `"colorchecker24"`, a custom `data.frame` with reference `R, G, B` columns,
#'   or `NULL` to extract observed colors without reference metrics.
#' @param nrow The number of rows of color patches. If `NULL` (default), automatically
#'   determined based on card orientation (4 for landscape, 6 for portrait for 24-patch charts).
#' @param ncol The number of columns of color patches. If `NULL` (default), automatically
#'   determined based on card orientation (6 for landscape, 4 for portrait).
#' @param buffer Negative buffer shrinkage fraction (0 to 1, default `0.45` for 45% shrink).
#'   If provided as a percentage > 1 (e.g. 45), it is automatically scaled to fraction.
#' @param margin_x Outer horizontal border margin fraction (default `0.05` / 5%).
#' @param margin_y Outer vertical border margin fraction (default `0.07` / 7%).
#' @param corners Optional quadrilateral corners. Can be `NULL` (automatic detection),
#'   a 4x2 matrix, a list with `$TL, $TR, $BR, $BL`, or an 8-element numeric vector.
#' @param index Segmentation index used for auto-detecting the card (default `"GRAY"`).
#' @param erode Erosion kernel size for segmentation mask. If `NULL` (default),
#'   automatically computed based on card coverage.
#' @param fill_hull Logical. Fill holes in the segmentation mask (default `TRUE`).
#' @param stat Primary statistic assigned to `R, G, B` columns: `"mean"` (default)
#'   or `"median"`.
#' @param auto_orient Logical. If `TRUE` (default for 24-patch ColorChecker charts),
#'   automatically detects the physical chart orientation (0, 90, 180, 270 degrees) by
#'   identifying the White (patch 19) and Black (patch 24) corners, ensuring that
#'   patch 1 (Dark skin) to patch 24 (Black) are always indexed correctly regardless
#'   of camera tilt or card rotation.
#' @param plot Logical. If `TRUE` (default), displays the image with overlay showing
#'   the detected card boundary, cell grid, negative buffer sampling polygons, and patch badges.
#' @param verbose Logical. If `TRUE` (default), prints diagnostic summary.
#' @param ... Legacy arguments (`xpix`, `ypix`, `filter_order`) supported for backward compatibility.
#'
#' @return An object of class `c("pliman_card_colors", "data.frame")` containing:
#' \itemize{
#'   \item `id`: Patch identifier (1 to `nrow * ncol`).
#'   \item `row`, `col`: Grid row and column positions.
#'   \item `center_x`, `center_y`: Pixel center coordinates.
#'   \item `R`, `G`, `B`: Sampled primary RGB values (0-255 scale).
#'   \item `R_med`, `G_med`, `B_med`: Median RGB values.
#'   \item `R_sd`, `G_sd`, `B_sd`: Standard deviations per channel.
#'   \item `n_pixels`: Number of pixels sampled in the nucleus.
#'   \item If reference comparison is active: `code`, `name`, `ref_R`, `ref_G`, `ref_B`,
#'     `Delta_R`, `Delta_G`, `Delta_B`, `Delta_E`, and `Delta_E_Calib`.
#' }
#' The object also stores attributes:
#' \itemize{
#'   \item `attr(res, "patches")`: List of patch bounding boxes and polygon vertices.
#'   \item `attr(res, "corners")`: Sorted corner coordinates (`TL, TR, BR, BL`).
#'   \item `attr(res, "wb_gains")`: Calculated White Balance gains.
#'   \item `attr(res, "ccm")`: Calculated 3x3 Color Correction Matrix.
#'   \item `attr(res, "avg_delta_raw")`: Mean Euclidean color error (\eqn{\Delta E}) before calibration.
#'   \item `attr(res, "avg_delta_calib")`: Mean \eqn{\Delta E} after CCM 3x3 calibration.
#' }
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("colorcheck.jpg")
#'
#' # Auto-detect card, auto-orient, and extract 24 patches with 45% negative buffer
#' card <- get_card_colors(img)
#' head(card)
#'
#' # View patch overlay
#' plot(card)
#' }
get_card_colors <- function(img,
                            reference = "colorchecker24",
                            nrow = NULL,
                            ncol = NULL,
                            buffer = 0.45,
                            margin_x = 0.05,
                            margin_y = 0.07,
                            corners = NULL,
                            auto_orient = TRUE,
                            index = "GRAY",
                            erode = NULL,
                            fill_hull = TRUE,
                            stat = c("mean", "median"),
                            plot = TRUE,
                            verbose = TRUE,
                            ...) {
  stat <- match.arg(stat)
  dots <- list(...)

  # Handle buffer percentage
  if (!is.null(buffer) && is.numeric(buffer)) {
    if (buffer > 1.0) buffer <- buffer / 100
    buffer <- max(0.0, min(0.95, buffer))
  } else {
    buffer <- 0.45
  }

  arr_raw <- unclass(img)
  d <- dim(arr_raw)
  w <- d[1]
  h <- d[2]

  # Harmonize image data to 0-255 scale numeric matrices
  scale_mult <- if (is.raw(arr_raw)) 1.0 else if (max(arr_raw, na.rm = TRUE) <= 1.0) 255.0 else 1.0
  r_mat <- as.numeric(arr_raw[, , 1]) * scale_mult
  g_mat <- as.numeric(arr_raw[, , 2]) * scale_mult
  b_mat <- as.numeric(arr_raw[, , 3]) * scale_mult
  dim(r_mat) <- c(w, h)
  dim(g_mat) <- c(w, h)
  dim(b_mat) <- c(w, h)

  # Determine / Detect corners
  corners_sorted <- NULL
  if (!is.null(corners)) {
    corners_sorted <- order_quad_corners(corners)
  } else {
    # Multi-candidate detection across segmentation strategies
    detect_methods <- list(
      function() {
        bin <- image_binary(img, index = index, fill_hull = fill_hull, plot = FALSE, verbose = FALSE)[[1]]
        if (!is.null(erode) && is.numeric(erode) && erode > 0) bin <- image_erode(bin, size = erode)
        bin
      },
      function() {
        # Dark plastic frame segmentation (ColorChecker frame is dark)
        gray <- (0.299 * r_mat + 0.587 * g_mat + 0.114 * b_mat) / 255
        dark_bin <- (gray < 0.35)
        tryCatch(image_fill_hull(dark_bin), error = function(e) dark_bin)
      },
      function() {
        # High-saturation color patches
        max_c <- pmax(r_mat, g_mat, b_mat)
        min_c <- pmin(r_mat, g_mat, b_mat)
        sat <- (max_c - min_c) / pmax(max_c, 1)
        sat_bin <- (sat > 0.15)
        tryCatch(image_fill_hull(sat_bin), error = function(e) sat_bin)
      }
    )

    best_candidate <- NULL
    best_score <- -1e9
    min_area <- (w * h) * 0.005

    for (method_fn in detect_methods) {
      bin_mask <- tryCatch(method_fn(), error = function(e) NULL)
      if (is.null(bin_mask)) next

      lab <- image_bwlabel(bin_mask)
      tab <- table(as.vector(lab))
      tab <- tab[names(tab) != "0"]
      if (length(tab) == 0) next

      valid_ids <- as.integer(names(tab[tab >= min_area]))
      if (length(valid_ids) == 0) next

      for (cand_id in valid_ids) {
        mask_cand <- as.matrix(lab == cand_id)
        cnts <- extract_contours_cpp(mask_cand)
        if (length(cnts) == 0) next

        raw_cnrs <- find_card_corners_cpp(cnts[[1]])
        if (is.null(raw_cnrs) || nrow(raw_cnrs) < 4) next

        cnrs <- order_quad_corners(raw_cnrs)
        c_tl_cand <- cnrs$TL; c_tr_cand <- cnrs$TR; c_br_cand <- cnrs$BR; c_bl_cand <- cnrs$BL

        # Quadrilateral Area (Shoelace formula)
        x_q <- c(c_tl_cand[1], c_tr_cand[1], c_br_cand[1], c_bl_cand[1])
        y_q <- c(c_tl_cand[2], c_tr_cand[2], c_br_cand[2], c_bl_cand[2])
        quad_area <- 0.5 * abs(sum(x_q * c(y_q[2:4], y_q[1]) - y_q * c(x_q[2:4], x_q[1])))
        if (quad_area < min_area) next

        # Side lengths & Aspect ratio
        s_top <- sqrt(sum((c_tr_cand - c_tl_cand)^2))
        s_bot <- sqrt(sum((c_br_cand - c_bl_cand)^2))
        s_left <- sqrt(sum((c_bl_cand - c_tl_cand)^2))
        s_right <- sqrt(sum((c_br_cand - c_tr_cand)^2))

        avg_w <- (s_top + s_bot) / 2
        avg_h <- (s_left + s_right) / 2
        aspect <- max(avg_w, avg_h) / max(min(avg_w, avg_h), 1)

        # Aspect ratio score: Standard ColorChecker 24 is ~1.5 (6:4)
        aspect_score <- 1 - min(1, abs(aspect - 1.5) / 1.5)

        # Rectangularity score: mask area vs quad area
        obj_area <- tab[[as.character(cand_id)]]
        rect_score <- min(1, obj_area / quad_area)

        # Corner contrast score
        center_cand <- (c_tl_cand + c_tr_cand + c_br_cand + c_bl_cand) / 4
        sample_lum_cand <- function(pt) {
          p <- pt * 0.85 + center_cand * 0.15
          x0 <- max(1, min(w, round(p[1] - 5)))
          x1 <- max(1, min(w, round(p[1] + 5)))
          y0 <- max(1, min(h, round(p[2] - 5)))
          y1 <- max(1, min(h, round(p[2] + 5)))
          sub_r <- r_mat[x0:x1, y0:y1]
          sub_g <- g_mat[x0:x1, y0:y1]
          sub_b <- b_mat[x0:x1, y0:y1]
          r_val <- if (length(sub_r) > 0) mean(sub_r) else 0
          g_val <- if (length(sub_g) > 0) mean(sub_g) else 0
          b_val <- if (length(sub_b) > 0) mean(sub_b) else 0
          0.299 * r_val + 0.587 * g_val + 0.114 * b_val
        }

        corner_lums <- sapply(list(c_tl_cand, c_tr_cand, c_br_cand, c_bl_cand), sample_lum_cand)
        lum_range <- max(corner_lums) - min(corner_lums)
        contrast_score <- min(1, lum_range / 150)

        total_score <- aspect_score * 0.3 + rect_score * 0.4 + contrast_score * 0.3

        if (total_score > best_score && rect_score > 0.55) {
          best_score <- total_score
          best_candidate <- cnrs
        }
      }
      if (!is.null(best_candidate) && best_score > 0.7) break
    }

    if (is.null(best_candidate)) {
      best_candidate <- list(
        TL = c(w * 0.1, h * 0.1),
        TR = c(w * 0.9, h * 0.1),
        BR = c(w * 0.9, h * 0.9),
        BL = c(w * 0.1, h * 0.9)
      )
    }
    corners_sorted <- best_candidate
  }

  # Check whether we should perform universal orientation optimization
  is_standard_24 <- (is.null(reference) || (is.character(reference) && reference %in% c("colorchecker24", "macbeth")) ||
                     (!is.null(reference) && is.data.frame(reference) && nrow(reference) == 24)) &&
                    (is.null(nrow) || is.null(ncol) || (nrow * ncol == 24))

  ref_target_mat <- NULL
  if (is_standard_24) {
    ref_target_df <- if (is.data.frame(reference)) reference else colorchecker_ref()
    ref_target_mat <- as.matrix(ref_target_df[, c("R", "G", "B")])
  }

  if (isTRUE(auto_orient) && !is.null(ref_target_mat)) {
    P <- list(
      corners_sorted$TL,
      corners_sorted$TR,
      corners_sorted$BR,
      corners_sorted$BL
    )

    # 12 topological walk configurations (4x6 and 6x4 across all rotations & reflections)
    configs <- list(
      # 4 rows x 6 cols
      list(nrow = 4, ncol = 6, corners = list(TL = P[[1]], TR = P[[2]], BR = P[[3]], BL = P[[4]])),
      list(nrow = 4, ncol = 6, corners = list(TL = P[[2]], TR = P[[3]], BR = P[[4]], BL = P[[1]])),
      list(nrow = 4, ncol = 6, corners = list(TL = P[[3]], TR = P[[4]], BR = P[[1]], BL = P[[2]])),
      list(nrow = 4, ncol = 6, corners = list(TL = P[[4]], TR = P[[1]], BR = P[[2]], BL = P[[3]])),
      list(nrow = 4, ncol = 6, corners = list(TL = P[[2]], TR = P[[1]], BR = P[[4]], BL = P[[3]])),
      list(nrow = 4, ncol = 6, corners = list(TL = P[[4]], TR = P[[3]], BR = P[[2]], BL = P[[1]])),

      # 6 rows x 4 cols
      list(nrow = 6, ncol = 4, corners = list(TL = P[[1]], TR = P[[2]], BR = P[[3]], BL = P[[4]])),
      list(nrow = 6, ncol = 4, corners = list(TL = P[[2]], TR = P[[3]], BR = P[[4]], BL = P[[1]])),
      list(nrow = 6, ncol = 4, corners = list(TL = P[[3]], TR = P[[4]], BR = P[[1]], BL = P[[2]])),
      list(nrow = 6, ncol = 4, corners = list(TL = P[[4]], TR = P[[1]], BR = P[[2]], BL = P[[3]])),
      list(nrow = 6, ncol = 4, corners = list(TL = P[[2]], TR = P[[1]], BR = P[[4]], BL = P[[3]])),
      list(nrow = 6, ncol = 4, corners = list(TL = P[[4]], TR = P[[3]], BR = P[[2]], BL = P[[1]]))
    )

    best_score <- Inf
    best_cfg <- configs[[1]]

    for (cfg in configs) {
      c_tl_i <- cfg$corners$TL; c_tr_i <- cfg$corners$TR; c_br_i <- cfg$corners$BR; c_bl_i <- cfg$corners$BL
      map_i <- function(u, v) {
        (1 - u) * (1 - v) * c_tl_i + u * (1 - v) * c_tr_i + u * v * c_br_i + (1 - u) * v * c_bl_i
      }

      u_s <- margin_x; u_e <- 1 - margin_x
      v_s <- margin_y; v_e <- 1 - margin_y
      pw <- (u_e - u_s) / cfg$ncol
      ph <- (v_e - v_s) / cfg$nrow
      hw <- (pw / 2) * (1 - buffer)
      hh <- (ph / 2) * (1 - buffer)

      obs_cand <- matrix(0, nrow = 24, ncol = 3)
      idx_c <- 1
      for (r_i in 1:cfg$nrow) {
        for (c_i in 1:cfg$ncol) {
          uc <- u_s + (c_i - 0.5) * pw
          vc <- v_s + (r_i - 0.5) * ph
          p_tl_c <- map_i(uc - hw, vc - hh)
          p_br_c <- map_i(uc + hw, vc + hh)
          p_tr_c <- map_i(uc + hw, vc - hh)
          p_bl_c <- map_i(uc - hw, vc + hh)
          xs <- c(p_tl_c[1], p_tr_c[1], p_br_c[1], p_bl_c[1])
          ys <- c(p_tl_c[2], p_tr_c[2], p_br_c[2], p_bl_c[2])
          x0 <- max(1, min(w, floor(min(xs)))); x1 <- max(1, min(w, ceiling(max(xs))))
          y0 <- max(1, min(h, floor(min(ys)))); y1 <- max(1, min(h, ceiling(max(ys))))
          if (x0 > x1) { t_tmp <- x0; x0 <- x1; x1 <- t_tmp }
          if (y0 > y1) { t_tmp <- y0; y0 <- y1; y1 <- t_tmp }
          sub_r <- r_mat[x0:x1, y0:y1]; sub_g <- g_mat[x0:x1, y0:y1]; sub_b <- b_mat[x0:x1, y0:y1]
          r_m <- if (length(sub_r) > 0) mean(sub_r) else 0
          g_m <- if (length(sub_g) > 0) mean(sub_g) else 0
          b_m <- if (length(sub_b) > 0) mean(sub_b) else 0
          obs_cand[idx_c, ] <- c(r_m, g_m, b_m)
          idx_c <- idx_c + 1
        }
      }
      d_e <- mean(sqrt(rowSums((obs_cand - ref_target_mat)^2)))
      if (!is.na(d_e) && d_e < best_score) {
        best_score <- d_e
        best_cfg <- cfg
      }
    }

    nrow <- best_cfg$nrow
    ncol <- best_cfg$ncol
    c_tl <- best_cfg$corners$TL
    c_tr <- best_cfg$corners$TR
    c_br <- best_cfg$corners$BR
    c_bl <- best_cfg$corners$BL
  } else {
    c_tl <- corners_sorted$TL
    c_tr <- corners_sorted$TR
    c_br <- corners_sorted$BR
    c_bl <- corners_sorted$BL

    w_top <- sqrt(sum((c_tr - c_tl)^2))
    w_bot <- sqrt(sum((c_br - c_bl)^2))
    h_left <- sqrt(sum((c_bl - c_tl)^2))
    h_right <- sqrt(sum((c_br - c_tr)^2))
    card_w <- (w_top + w_bot) / 2
    card_h <- (h_left + h_right) / 2

    if (is.null(nrow) || is.null(ncol)) {
      if (card_w >= card_h) {
        if (is.null(nrow)) nrow <- 4
        if (is.null(ncol)) ncol <- 6
      } else {
        if (is.null(nrow)) nrow <- 6
        if (is.null(ncol)) ncol <- 4
      }
    }
  }

  # Bilinear mapping function in canonical card space:
  # u in [0, 1] across cols (1..ncol)
  # v in [0, 1] across rows (1..nrow)
  map_uv <- function(u, v) {
    (1 - u) * (1 - v) * c_tl +
    u * (1 - v) * c_tr +
    u * v * c_br +
    (1 - u) * v * c_bl
  }

  u_start <- margin_x
  u_end   <- 1 - margin_x
  v_start <- margin_y
  v_end   <- 1 - margin_y

  patch_w_norm <- (u_end - u_start) / ncol
  patch_h_norm <- (v_end - v_start) / nrow

  # Check for legacy xpix / ypix overrides
  legacy_xpix <- dots$xpix
  legacy_ypix <- dots$ypix

  results <- vector("list", nrow * ncol)
  idx <- 1

  for (r_idx in 1:nrow) {
    for (c_idx in 1:ncol) {
      u_c <- u_start + (c_idx - 0.5) * patch_w_norm
      v_c <- v_start + (r_idx - 0.5) * patch_h_norm

      half_w <- (patch_w_norm / 2) * (1 - buffer)
      half_h <- (patch_h_norm / 2) * (1 - buffer)

      cell_tl <- map_uv(u_start + (c_idx - 1) * patch_w_norm, v_start + (r_idx - 1) * patch_h_norm)
      cell_tr <- map_uv(u_start + c_idx * patch_w_norm, v_start + (r_idx - 1) * patch_h_norm)
      cell_br <- map_uv(u_start + c_idx * patch_w_norm, v_start + r_idx * patch_h_norm)
      cell_bl <- map_uv(u_start + (c_idx - 1) * patch_w_norm, v_start + r_idx * patch_h_norm)

      p_tl <- map_uv(u_c - half_w, v_c - half_h)
      p_tr <- map_uv(u_c + half_w, v_c - half_h)
      p_br <- map_uv(u_c + half_w, v_c + half_h)
      p_bl <- map_uv(u_c - half_w, v_c + half_h)
      p_center <- map_uv(u_c, v_c)

      if (!is.null(legacy_xpix) && !is.null(legacy_ypix)) {
        x_min <- max(1, min(w, round(p_center[1] - legacy_xpix / 2)))
        x_max <- max(1, min(w, round(p_center[1] + legacy_xpix / 2)))
        y_min <- max(1, min(h, round(p_center[2] - legacy_ypix / 2)))
        y_max <- max(1, min(h, round(p_center[2] + legacy_ypix / 2)))
      } else {
        xs <- c(p_tl[1], p_tr[1], p_br[1], p_bl[1])
        ys <- c(p_tl[2], p_tr[2], p_br[2], p_bl[2])
        x_min <- max(1, min(w, floor(min(xs))))
        x_max <- max(1, min(w, ceiling(max(xs))))
        y_min <- max(1, min(h, floor(min(ys))))
        y_max <- max(1, min(h, ceiling(max(ys))))
      }

      if (x_min > x_max) { tmp <- x_min; x_min <- x_max; x_max <- tmp }
      if (y_min > y_max) { tmp <- y_min; y_min <- y_max; y_max <- tmp }

      sub_r <- r_mat[x_min:x_max, y_min:y_max, drop = FALSE]
      sub_g <- g_mat[x_min:x_max, y_min:y_max, drop = FALSE]
      sub_b <- b_mat[x_min:x_max, y_min:y_max, drop = FALSE]

      r_mean <- if (length(sub_r) > 0) mean(sub_r) else 0
      g_mean <- if (length(sub_g) > 0) mean(sub_g) else 0
      b_mean <- if (length(sub_b) > 0) mean(sub_b) else 0

      r_med <- if (length(sub_r) > 0) stats::median(sub_r) else 0
      g_med <- if (length(sub_g) > 0) stats::median(sub_g) else 0
      b_med <- if (length(sub_b) > 0) stats::median(sub_b) else 0

      r_sd <- if (length(sub_r) > 1) stats::sd(sub_r) else 0
      g_sd <- if (length(sub_g) > 1) stats::sd(sub_g) else 0
      b_sd <- if (length(sub_b) > 1) stats::sd(sub_b) else 0

      results[[idx]] <- list(
        id = idx,
        row = r_idx,
        col = c_idx,
        center_x = p_center[1],
        center_y = p_center[2],
        cell_x = c(cell_tl[1], cell_tr[1], cell_br[1], cell_bl[1]),
        cell_y = c(cell_tl[2], cell_tr[2], cell_br[2], cell_bl[2]),
        poly_x = c(p_tl[1], p_tr[1], p_br[1], p_bl[1]),
        poly_y = c(p_tl[2], p_tr[2], p_br[2], p_bl[2]),
        R_obs = round(r_mean, 1),
        G_obs = round(g_mean, 1),
        B_obs = round(b_mean, 1),
        R_med = round(r_med, 1),
        G_med = round(g_med, 1),
        B_med = round(b_med, 1),
        R_sd  = round(r_sd, 2),
        G_sd  = round(g_sd, 2),
        B_sd  = round(b_sd, 2),
        n_pixels = length(sub_r)
      )
      idx <- idx + 1
    }
  }

  # Build data frame
  df_obs <- do.call(rbind, lapply(results, function(x) {
    primary_r <- if (stat == "median") x$R_med else x$R_obs
    primary_g <- if (stat == "median") x$G_med else x$G_obs
    primary_b <- if (stat == "median") x$B_med else x$B_obs
    data.frame(
      id = x$id,
      row = x$row,
      col = x$col,
      center_x = round(x$center_x, 1),
      center_y = round(x$center_y, 1),
      R = primary_r,
      G = primary_g,
      B = primary_b,
      R_mean = x$R_obs,
      G_mean = x$G_obs,
      B_mean = x$B_obs,
      R_med = x$R_med,
      G_med = x$G_med,
      B_med = x$B_med,
      R_sd = x$R_sd,
      G_sd = x$G_sd,
      B_sd = x$B_sd,
      n_pixels = x$n_pixels,
      stringsAsFactors = FALSE
    )
  }))

  # Reference matching & error calculation
  ref_df <- NULL
  if (!is.null(reference)) {
    if (is.character(reference) && reference %in% c("colorchecker24", "macbeth") && nrow * ncol == 24) {
      ref_df <- colorchecker_ref(reference)
    } else if (is.data.frame(reference) && nrow(reference) == nrow(df_obs)) {
      ref_df <- reference
      if (!all(c("R", "G", "B") %in% colnames(ref_df)) && all(c("ref_R", "ref_G", "ref_B") %in% colnames(ref_df))) {
        ref_df$R <- ref_df$ref_R
        ref_df$G <- ref_df$ref_G
        ref_df$B <- ref_df$ref_B
      }
    }
  }

  wb_gains <- NULL
  ccm_mat <- NULL
  avg_delta_raw <- NULL
  avg_delta_calib <- NULL

  if (!is.null(ref_df)) {
    df_obs$code <- if (!is.null(ref_df$code)) ref_df$code else sprintf("No.%03d", df_obs$id)
    df_obs$name <- if (!is.null(ref_df$name)) ref_df$name else sprintf("Patch %d", df_obs$id)
    df_obs$ref_R <- ref_df$R
    df_obs$ref_G <- ref_df$G
    df_obs$ref_B <- ref_df$B

    df_obs$Delta_R <- round(df_obs$R - df_obs$ref_R, 1)
    df_obs$Delta_G <- round(df_obs$G - df_obs$ref_G, 1)
    df_obs$Delta_B <- round(df_obs$B - df_obs$ref_B, 1)
    df_obs$Delta_E <- round(sqrt(df_obs$Delta_R^2 + df_obs$Delta_G^2 + df_obs$Delta_B^2), 1)

    obs_mat <- as.matrix(df_obs[, c("R", "G", "B")])
    ref_mat <- as.matrix(df_obs[, c("ref_R", "ref_G", "ref_B")])

    # White balance gains from neutral patch (Patch 19 White or patch with highest luminance)
    white_idx <- if (nrow(df_obs) >= 19) 19 else which.max(rowMeans(ref_mat))
    w_obs <- obs_mat[white_idx, ]
    wb_gains <- if (all(w_obs > 0)) w_obs[2] / w_obs else c(1, 1, 1)

    # 3x3 CCM with Ridge Regularization
    ccm_mat <- tryCatch({
      solve(t(obs_mat) %*% obs_mat + diag(1e-4, 3)) %*% t(obs_mat) %*% ref_mat
    }, error = function(e) diag(1, 3))

    calib_mat <- matrix(pmax(0, pmin(255, obs_mat %*% ccm_mat)), nrow = nrow(obs_mat), ncol = 3)
    delta_e_calib <- round(sqrt(rowSums((calib_mat - ref_mat)^2)), 1)
    df_obs$Delta_E_Calib <- delta_e_calib

    avg_delta_raw <- round(mean(df_obs$Delta_E), 1)
    avg_delta_calib <- round(mean(delta_e_calib), 1)
  }

  class(df_obs) <- c("pliman_card_colors", "data.frame")
  attr(df_obs, "patches") <- results
  attr(df_obs, "corners") <- corners_sorted
  attr(df_obs, "wb_gains") <- wb_gains
  attr(df_obs, "ccm") <- ccm_mat
  attr(df_obs, "avg_delta_raw") <- avg_delta_raw
  attr(df_obs, "avg_delta_calib") <- avg_delta_calib
  attr(df_obs, "reference") <- ref_df
  attr(df_obs, "img_dim") <- c(w, h)
  attr(df_obs, "image") <- img

  if (isTRUE(plot)) {
    plot.pliman_card_colors(df_obs, img = img, type = "overlay")
  }

  if (isTRUE(verbose) && !is.null(ref_df)) {
    cli::cli_alert_success(
      "Extracted {.val {nrow(df_obs)}} patches | Raw {.field Delta E} = {.val {avg_delta_raw}} | Calibrated (CCM) {.field Delta E} = {.val {avg_delta_calib}}"
    )
  }

  return(df_obs)
}


#' @export
plot.pliman_card_colors <- function(x,
                                   img = NULL,
                                   type = c("overlay", "swatches"),
                                   ...) {
  type <- match.arg(type)
  if (is.null(img)) {
    img <- attr(x, "image")
  }
  if (is.null(img)) {
    cli::cli_abort("An {.arg img} object must be provided to plot {.fn plot.pliman_card_colors}.")
  }

  if (type == "overlay") {
    plot(img)
    patches <- attr(x, "patches")
    corners <- attr(x, "corners")

    if (!is.null(corners)) {
      c_tl <- corners$TL; c_tr <- corners$TR; c_br <- corners$BR; c_bl <- corners$BL
      graphics::polygon(
        c(c_tl[1], c_tr[1], c_br[1], c_bl[1]),
        c(c_tl[2], c_tr[2], c_br[2], c_bl[2]),
        border = "#00E5FF", lwd = 2.5, lty = "dashed"
      )
    }

    if (!is.null(patches)) {
      for (p in patches) {
        if (!is.null(p$cell_x) && !is.null(p$cell_y)) {
          graphics::polygon(p$cell_x, p$cell_y, border = grDevices::rgb(1, 1, 1, 0.4), lwd = 1)
        }
        if (!is.null(p$poly_x) && !is.null(p$poly_y)) {
          graphics::polygon(p$poly_x, p$poly_y, border = "#00FF66",
                            col = grDevices::rgb(0, 1, 0.4, 0.22), lwd = 1.8)
        }
        graphics::points(p$center_x, p$center_y, pch = 21, bg = "#000000", col = "#00FF66", cex = 2.0, lwd = 1.5)
        graphics::text(p$center_x, p$center_y, sprintf("%02d", p$id), col = "#FFFFFF", cex = 0.8, font = 2)
      }
    }
  } else if (type == "swatches") {
    n_patches <- nrow(x)
    has_ref <- all(c("ref_R", "ref_G", "ref_B") %in% colnames(x))

    op <- graphics::par(mar = c(1, 1, 2, 1), bg = "#f8fafc")
    on.exit(graphics::par(op), add = TRUE)

    n_cols <- min(6, n_patches)
    n_rows <- ceiling(n_patches / n_cols)

    graphics::plot.new()
    graphics::plot.window(xlim = c(0, n_cols), ylim = c(n_rows, 0))
    graphics::title(main = "ColorChecker Patches: Measured vs Reference Swatches", col.main = "#0f172a", font.main = 2)

    for (i in seq_len(n_patches)) {
      r_idx <- ceiling(i / n_cols)
      c_idx <- (i - 1) %% n_cols + 1

      x0 <- c_idx - 0.95; x1 <- c_idx - 0.05
      y0 <- r_idx - 0.95; y1 <- r_idx - 0.05

      graphics::rect(x0, y0, x1, y1, col = "#ffffff", border = "#cbd5e1", lwd = 1)

      obs_col <- grDevices::rgb(
        min(255, max(0, x$R[i])),
        min(255, max(0, x$G[i])),
        min(255, max(0, x$B[i])),
        maxColorValue = 255
      )

      if (has_ref) {
        ref_col <- grDevices::rgb(
          min(255, max(0, x$ref_R[i])),
          min(255, max(0, x$ref_G[i])),
          min(255, max(0, x$ref_B[i])),
          maxColorValue = 255
        )
        graphics::rect(x0 + 0.05, y0 + 0.25, x0 + 0.43, y1 - 0.05, col = ref_col, border = "#94a3b8")
        graphics::rect(x0 + 0.47, y0 + 0.25, x1 - 0.05, y1 - 0.05, col = obs_col, border = "#94a3b8")
        de_txt <- if (!is.null(x$Delta_E)) sprintf("\u0394E: %.1f", x$Delta_E[i]) else ""
        graphics::text(mean(c(x0, x1)), y0 + 0.14, sprintf("#%02d %s", x$id[i], de_txt),
                       cex = 0.7, font = 2, col = "#1e293b")
      } else {
        graphics::rect(x0 + 0.05, y0 + 0.25, x1 - 0.05, y1 - 0.05, col = obs_col, border = "#94a3b8")
        graphics::text(mean(c(x0, x1)), y0 + 0.14, sprintf("#%02d", x$id[i]),
                       cex = 0.75, font = 2, col = "#1e293b")
      }
    }
  }

  invisible(x)
}

#' @export
print.pliman_card_colors <- function(x, ...) {
  cat(cli::col_cyan(cli::style_bold("\n-- ColorChecker Card Extraction --\n")))
  cat(sprintf("Number of patches: %d\n", nrow(x)))

  wb <- attr(x, "wb_gains")
  if (!is.null(wb)) {
    cat(sprintf("White Balance Gains: R = %.3f, G = %.3f, B = %.3f\n", wb[1], wb[2], wb[3]))
  }

  raw_de <- attr(x, "avg_delta_raw")
  calib_de <- attr(x, "avg_delta_calib")
  if (!is.null(raw_de)) {
    cat(sprintf("Average Delta E (Raw):   %.2f\n", raw_de))
  }
  if (!is.null(calib_de)) {
    cat(sprintf("Average Delta E (CCM):   %.2f\n", calib_de))
  }
  cat("\n")
  print.data.frame(head(as.data.frame(x), 10), ...)
  if (nrow(x) > 10) {
    cat(cli::col_grey(sprintf("... and %d more rows (use as.data.frame() or [ ] to view all)\n", nrow(x) - 10)))
  }
  invisible(x)
}

#' @export
summary.pliman_card_colors <- function(object, ...) {
  cat(cli::col_cyan(cli::style_bold("\n-- ColorChecker Extraction Summary --\n")))
  cat(sprintf("Patches sampled: %d\n", nrow(object)))
  cat(sprintf("Sampling pixels per patch: min = %d, median = %d, max = %d\n",
              min(object$n_pixels), as.integer(stats::median(object$n_pixels)), max(object$n_pixels)))

  if ("Delta_E" %in% colnames(object)) {
    cat("\nEuclidean Color Error (Delta E):\n")
    print(summary(object$Delta_E))
  }
  if ("Delta_E_Calib" %in% colnames(object)) {
    cat("\nCalibrated Delta E (CCM 3x3):\n")
    print(summary(object$Delta_E_Calib))
  }
  invisible(object)
}


#' @title Correct Image Colors using a Color Checker
#'
#' @description
#' Calibrates the color of an image using a set of known color references (e.g.,
#' from a color checker) by implementing multiple calibration models including
#' industry-standard **CCM (Color Correction Matrix 3x3)**, Affine (4x3),
#' White Balance, Cubic (9-term polynomial), and Root-Polynomial (20-term).
#'
#' @details
#' Color calibration finds an optimal transformation matrix \eqn{\mathbf{K}} mapping observed
#' colors \eqn{\mathbf{S}} to reference target colors \eqn{\mathbf{T}}:
#' \itemize{
#'   \item \code{"ccm"} (default): Linear 3x3 Color Correction Matrix computed via Ridge-regularized
#'     least squares:
#'     \deqn{\mathbf{K} = (\mathbf{S}^T \mathbf{S} + \lambda \mathbf{I})^{-1} \mathbf{S}^T \mathbf{T}}
#'   \item \code{"affine"}: 4x3 linear transformation with translation offset (intercept) \eqn{[1, R, G, B]}.
#'   \item \code{"white_balance"}: Channel-wise gain correction derived from neutral/white patches.
#'   \item \code{"cubic"}: 9-term polynomial mapping (\eqn{R, G, B, R^2, G^2, B^2, R^3, G^3, B^3}).
#'   \item \code{"root_polynomial"}: 20-term root-polynomial mapping with cross-channel interactions.
#' }
#'
#' @param img An `image` object to be corrected.
#' @param card_colors A data frame of observed colors from the color checker
#'   (typically output directly from `get_card_colors()`). Must contain `R, G, B`
#'   (or `R_obs, G_obs, B_obs`) columns. If `NULL`, `k_mat` must be provided.
#' @param known_colors A data frame of known reference colors. Must contain `R, G, B`
#'   (or `ref_R, ref_G, ref_B`) columns. If `NULL` and `card_colors` contains reference
#'   data or `reference = "colorchecker24"`, the built-in reference chart is used automatically.
#' @param k_mat A pre-calculated transformation matrix (\eqn{\mathbf{K}}).
#'   If provided, `card_colors` and `known_colors` are optional, and the correction
#'   is applied directly to `img`. Useful for batch calibration of image series.
#' @param model The calibration model to use:
#'   \itemize{
#'     \item `"ccm"` (default): 3x3 linear Color Correction Matrix (robust, standard in computer vision).
#'     \item `"affine"`: 4x3 matrix with intercept.
#'     \item `"white_balance"`: Channel gains derived from neutral white patch.
#'     \item `"cubic"`: 9-term 3rd-order polynomial model.
#'     \item `"root_polynomial"`: 20-term root-polynomial model.
#'   }
#' @param reference Reference chart standard to use when `known_colors = NULL` (default `"colorchecker24"`).
#' @param lambda Ridge regularization factor for matrix inversion (default `1e-4`).
#' @param plot Logical. If `TRUE` (default `FALSE`), displays a side-by-side comparison
#'   of the original and calibrated image.
#' @param verbose Logical. If `TRUE` (default), prints calibration diagnostics.
#'
#' @return
#' If `card_colors` or `known_colors` are provided, returns a list containing:
#' \itemize{
#'   \item `img`: An `image` object with the color correction applied.
#'   \item `k`: The calculated transformation matrix (\eqn{\mathbf{K}}).
#'   \item `model`: The model used.
#'   \item `avg_delta_raw`: Mean Euclidean color error (\eqn{\Delta E}) before calibration.
#'   \item `avg_delta_calib`: Mean \eqn{\Delta E} after calibration.
#'   \item `wb_gains`: White balance channel gains.
#'   \item `df`: Comparison data frame with observed, calibrated, and reference colors.
#' }
#' If only `img` and a pre-calculated `k_mat` are supplied, returns the corrected `image` directly.
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("colorcheck.jpg")
#'
#' # 1. Extract card colors
#' card <- get_card_colors(img)
#'
#' # 2. Calibrate image using default CCM 3x3 model
#' calib <- image_correction(img, card, model = "ccm", plot = TRUE)
#'
#' # 3. Or using polynomial models
#' calib_poly <- image_correction(img, card, model = "root_polynomial")
#' }
image_correction <- function(img,
                             card_colors = NULL,
                             known_colors = NULL,
                             k_mat = NULL,
                             model = c("ccm", "affine", "white_balance", "cubic", "root_polynomial"),
                             reference = "colorchecker24",
                             lambda = 1e-4,
                             plot = FALSE,
                             verbose = TRUE) {
  model <- match.arg(model)

  if (is.null(card_colors) && is.null(k_mat)) {
    cli::cli_abort("Either {.arg card_colors} or a pre-calculated {.arg k_mat} must be provided.")
  }

  if (!is.null(card_colors)) {
    # Extract observed RGB matrix
    col_names <- colnames(card_colors)
    cc_r <- if ("R" %in% col_names) card_colors$R else if ("R_obs" %in% col_names) card_colors$R_obs else card_colors[[1]]
    cc_g <- if ("G" %in% col_names) card_colors$G else if ("G_obs" %in% col_names) card_colors$G_obs else card_colors[[2]]
    cc_b <- if ("B" %in% col_names) card_colors$B else if ("B_obs" %in% col_names) card_colors$B_obs else card_colors[[3]]
    cc_rgb <- cbind(as.numeric(cc_r), as.numeric(cc_g), as.numeric(cc_b))
    colnames(cc_rgb) <- c("R", "G", "B")

    # Determine reference colors
    kc_rgb <- NULL
    if (!is.null(known_colors)) {
      kc_col <- colnames(known_colors)
      kc_r <- if ("R" %in% kc_col) known_colors$R else if ("ref_R" %in% kc_col) known_colors$ref_R else known_colors[[1]]
      kc_g <- if ("G" %in% kc_col) known_colors$G else if ("ref_G" %in% kc_col) known_colors$ref_G else known_colors[[2]]
      kc_b <- if ("B" %in% kc_col) known_colors$B else if ("ref_B" %in% kc_col) known_colors$ref_B else known_colors[[3]]
      kc_rgb <- cbind(as.numeric(kc_r), as.numeric(kc_g), as.numeric(kc_b))
      colnames(kc_rgb) <- c("R", "G", "B")
    } else if (all(c("ref_R", "ref_G", "ref_B") %in% col_names)) {
      kc_rgb <- cbind(as.numeric(card_colors$ref_R), as.numeric(card_colors$ref_G), as.numeric(card_colors$ref_B))
      colnames(kc_rgb) <- c("R", "G", "B")
    } else if (!is.null(reference) && is.character(reference) && reference %in% c("colorchecker24", "macbeth") && nrow(card_colors) == 24) {
      ref_tab <- colorchecker_ref(reference)
      kc_rgb <- cbind(ref_tab$R, ref_tab$G, ref_tab$B)
      colnames(kc_rgb) <- c("R", "G", "B")
    } else {
      cli::cli_abort("Reference colors could not be determined. Please provide {.arg known_colors}.")
    }

    if (nrow(cc_rgb) != nrow(kc_rgb)) {
      cli::cli_abort(
        c("Number of rows in {.arg card_colors} ({.val {nrow(cc_rgb)}}) does not match {.arg known_colors} ({.val {nrow(kc_rgb)}}).")
      )
    }

    # Harmonize scale to 0-255
    if (max(cc_rgb, na.rm = TRUE) <= 1.0) cc_rgb <- cc_rgb * 255
    if (max(kc_rgb, na.rm = TRUE) <= 1.0) kc_rgb <- kc_rgb * 255

    # Compute transformation matrix K based on model
    wb_gains <- NULL
    K <- NULL

    if (model == "ccm") {
      K <- solve(t(cc_rgb) %*% cc_rgb + diag(lambda, 3)) %*% t(cc_rgb) %*% kc_rgb
    } else if (model == "affine") {
      s_aff <- cbind(1, cc_rgb)
      K <- solve(t(s_aff) %*% s_aff + diag(lambda, 4)) %*% t(s_aff) %*% kc_rgb
    } else if (model == "white_balance") {
      # Use patch 19 (white) or highest luminance patch
      white_idx <- if (nrow(cc_rgb) >= 19) 19 else which.max(rowMeans(kc_rgb))
      w_obs <- cc_rgb[white_idx, ]
      w_ref <- kc_rgb[white_idx, ]
      wb_gains <- if (all(w_obs > 0)) w_ref / w_obs else c(1, 1, 1)
      K <- matrix(wb_gains, nrow = 1, ncol = 3)
    } else if (model == "cubic") {
      s_cub <- create_poly_matrix(as.data.frame(cc_rgb), model = "cubic")
      K <- mpinv(s_cub) %*% kc_rgb
    } else if (model == "root_polynomial") {
      s_root <- create_poly_matrix(as.data.frame(cc_rgb), model = "root_polynomial")
      K <- mpinv(s_root) %*% kc_rgb
    }

    # Apply correction via C++
    res_raw <- correct_image_rcpp(img, K, model = model)
    res_clamped <- pmin(pmax(round(res_raw), 0), 255)
    corrected_img <- as_image(res_clamped, dim = dim(img), colormode = 'Color')

    # Compute calibrated colors and delta E
    calib_raw <- if (model == "ccm") {
      cc_rgb %*% K
    } else if (model == "affine") {
      cbind(1, cc_rgb) %*% K
    } else if (model == "white_balance") {
      sweep(cc_rgb, 2, wb_gains, "*")
    } else if (model == "cubic") {
      create_poly_matrix(as.data.frame(cc_rgb), model = "cubic") %*% K
    } else if (model == "root_polynomial") {
      create_poly_matrix(as.data.frame(cc_rgb), model = "root_polynomial") %*% K
    }
    calib_obs <- matrix(pmax(0, pmin(255, calib_raw)), nrow = nrow(cc_rgb), ncol = 3)

    delta_raw <- round(sqrt(rowSums((cc_rgb - kc_rgb)^2)), 1)
    delta_calib <- round(sqrt(rowSums((calib_obs - kc_rgb)^2)), 1)
    avg_delta_raw <- round(mean(delta_raw), 1)
    avg_delta_calib <- round(mean(delta_calib), 1)

    df_comp <- data.frame(
      id = seq_len(nrow(cc_rgb)),
      R_obs = round(cc_rgb[, 1], 1),
      G_obs = round(cc_rgb[, 2], 1),
      B_obs = round(cc_rgb[, 3], 1),
      R_calib = round(calib_obs[, 1], 1),
      G_calib = round(calib_obs[, 2], 1),
      B_calib = round(calib_obs[, 3], 1),
      ref_R = round(kc_rgb[, 1], 1),
      ref_G = round(kc_rgb[, 2], 1),
      ref_B = round(kc_rgb[, 3], 1),
      Delta_E_Raw = delta_raw,
      Delta_E_Calib = delta_calib,
      stringsAsFactors = FALSE
    )

    if (isTRUE(verbose)) {
      cli::cli_alert_success(
        "Color correction applied ({.field model} = {.val {model}}) | Raw {.field Delta E} = {.val {avg_delta_raw}} -> Calibrated {.field Delta E} = {.val {avg_delta_calib}}"
      )
    }

    if (isTRUE(plot)) {
      image_combine(img, corrected_img, labels = c("Original", paste0("Calibrated (", model, ")")))
    }

    return(list(
      img = corrected_img,
      k = K,
      model = model,
      avg_delta_raw = avg_delta_raw,
      avg_delta_calib = avg_delta_calib,
      wb_gains = wb_gains,
      df = df_comp
    ))
  } else {
    # Direct application of pre-calculated matrix K
    if (model == "cubic" && nrow(k_mat) != 9) {
      cli::cli_abort("For model 'cubic', 'k_mat' must have 9 rows.")
    } else if (model == "root_polynomial" && nrow(k_mat) != 20) {
      cli::cli_abort("For model 'root_polynomial', 'k_mat' must have 20 rows.")
    } else if (model == "ccm" && (nrow(k_mat) != 3 || ncol(k_mat) != 3)) {
      cli::cli_abort("For model 'ccm', 'k_mat' must be a 3x3 matrix.")
    } else if (model == "affine" && (nrow(k_mat) != 4 || ncol(k_mat) != 3)) {
      cli::cli_abort("For model 'affine', 'k_mat' must be a 4x3 matrix.")
    }

    res_raw <- correct_image_rcpp(img, k_mat, model = model)
    res_clamped <- pmin(pmax(round(res_raw), 0), 255)
    corrected_img <- as_image(res_clamped, dim = dim(img), colormode = 'Color')

    if (isTRUE(plot)) {
      image_combine(img, corrected_img, labels = c("Original", paste0("Calibrated (", model, ")")))
    }

    return(corrected_img)
  }
}


#' Interactive Color Correction
#'
#' @description
#' Interactively calibrates the color of an image by prompting the user to
#' sample color patches. It calls `pliman::pick_rgb_area()` to pause and
#' prompt the user to click on `color_chips` patches in the image. These
#' sampled colors are then mapped to `known_colors` using the specified `model`.
#'
#' @param img An `image` object to be corrected.
#' @param known_colors A `data.frame` containing target reference values with `R, G, B`
#'   columns. If `NULL` and `color_chips == 24`, defaults to `colorchecker_ref()`.
#' @param color_chips The number of color patches to sample interactively (default `24`).
#' @param model The calibration model to use (`"ccm"`, `"affine"`, `"white_balance"`,
#'   `"cubic"`, or `"root_polynomial"`). Default `"ccm"`.
#' @param lambda Ridge regularization factor (default `1e-4`).
#'
#' @return An `image` object with corrected colors.
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#'
#' # Sample 4 known chips (e.g. blue, red, yellow, white)
#' known_colors <- data.frame(
#'  id = c(1, 2, 3, 4),
#'  R = c(25, 186, 245, 249),
#'  G = c(55, 26, 205, 242),
#'  B = c(135, 51, 0, 238)
#' )
#' img <- image_pliman("colorcheck.jpg")
#' img_cor <- image_correction_pick(img, known_colors, color_chips = 4)
#' image_combine(img, img_cor)
#' }
image_correction_pick <- function(img,
                                  known_colors = NULL,
                                  color_chips = 24,
                                  model = c("ccm", "affine", "white_balance", "cubic", "root_polynomial"),
                                  lambda = 1e-4) {
  model <- match.arg(model)

  if (is.null(known_colors)) {
    if (color_chips == 24) {
      known_colors <- colorchecker_ref()
    } else {
      cli::cli_abort("Please provide {.arg known_colors} matching the {.val {color_chips}} sampled chips.")
    }
  }

  if (nrow(known_colors) != color_chips) {
    cli::cli_abort(
      c("The number of rows in {.arg known_colors} ({.val {nrow(known_colors)}}) does not equal {.arg color_chips} ({.val {color_chips}}).")
    )
  }

  col_k <- colnames(known_colors)
  if (!all(c("R", "G", "B") %in% col_k) && !all(c("ref_R", "ref_G", "ref_B") %in% col_k)) {
    cli::cli_abort("The {.arg known_colors} data frame must contain columns R, G, and B.")
  }

  cli::cli_progress_step(
    msg        = "Pick {.val {color_chips}} color chips in the image corresponding to {.arg known_colors}...",
    msg_done   = "Color sampling done",
    msg_failed = "Oops, something went wrong."
  )
  rgbimg <- pick_rgb_area(img, n = color_chips, verbose = FALSE)

  cli::cli_progress_step(
    msg        = "Performing color correction ({model})...",
    msg_done   = "Color correction finished",
    msg_failed = "Oops, something went wrong."
  )

  res <- image_correction(img, card_colors = rgbimg, known_colors = known_colors,
                          model = model, lambda = lambda, verbose = FALSE)
  return(res$img)
}

check_filter_order <- function(filter_order, verbose, erode, dilate, opening, closing, filter, fill_hull) {
  valid_filters <- c("erode", "dilate", "opening", "closing", "filter", "fill_hull")
  invalid <- setdiff(filter_order, valid_filters)
  if (length(invalid) > 0) {
    cli::cli_warn("The following filter{?s} {?is/are} not recognized: {.val {invalid}}.")
  }
  if (isTRUE(verbose)) {
    used <- character()
    for (op in filter_order) {
      if (op %in% valid_filters) {
        if (op == "erode" && is.numeric(erode) && erode > 0) {
          used <- c(used, "erode")
        } else if (op == "dilate" && is.numeric(dilate) && dilate > 0) {
          used <- c(used, "dilate")
        } else if (op == "opening" && is.numeric(opening) && opening > 0) {
          used <- c(used, "opening")
        } else if (op == "closing" && is.numeric(closing) && closing > 0) {
          used <- c(used, "closing")
        } else if (op == "filter" && is.numeric(filter) && filter > 1) {
          used <- c(used, "filter")
        } else if (op == "fill_hull" && isTRUE(fill_hull)) {
          used <- c(used, "fill_hull")
        }
      }
    }
    if (length(used) > 0) {
      cli::cli_alert_info("Filter order: {.val {paste(used, collapse = ' -> ')}}")
    }
  }
  return(filter_order)
}

#' Fast Binary Image Prediction from a Logistic Regression Model
#'
#' Evaluates a fitted `glm` logistic regression model directly on image channel
#' matrices without creating massive data frames or running slow `predict.glm`
#' loops over millions of pixels.
#'
#' @param model A fitted `glm` logistic regression model object.
#' @param img An `image` object.
#' @return An integer matrix containing 1 (foreground) and 0 (background).
#' @export
predict_binary_glm <- function(model, img) {
  co <- stats::coef(model)
  co[is.na(co)] <- 0

  img_num <- image_data(img, type = "numeric")
  nrow_img <- dim(img)[1]
  ncol_img <- dim(img)[2]

  R_ch <- c(img_num[,,1])
  G_ch <- c(img_num[,,2])
  B_ch <- c(img_num[,,3])

  eval_env <- list(
    R = R_ch, G = G_ch, B = B_ch,
    r = R_ch, g = G_ch, b = B_ch
  )

  intercept_name <- "(Intercept)"
  eta <- if (intercept_name %in% names(co)) co[[intercept_name]] else 0

  term_names <- setdiff(names(co), intercept_name)
  for (tn in term_names) {
    b_val <- co[[tn]]
    if (b_val == 0) next
    expr_str <- gsub(":", "*", tn)
    expr <- parse(text = expr_str)[[1]]
    x_val <- eval(expr, envir = list2env(eval_env, parent = baseenv()))
    eta <- eta + b_val * x_val
  }

  matrix(as.integer(eta >= 0), nrow = nrow_img, ncol = ncol_img)
}

#' Create an RGB Image from Red, Green, and Blue Channels
#'
#' Combines single-channel matrices, arrays, or `image` objects representing
#' Red, Green, and Blue channels into a 3D color RGB `image` object.
#'
#' @param red,green,blue Single-channel matrices, arrays, or `image` objects.
#'   Any channel can be `NULL` (will be filled with zeros).
#' @return A 3D `image` object of color mode `"Color"`.
#' @export
rgb_image <- function(red = NULL, green = NULL, blue = NULL) {
  non_null <- list(red = red, green = green, blue = blue)
  non_null <- non_null[!sapply(non_null, is.null)]

  if (length(non_null) == 0) {
    cli::cli_abort("At least one channel (red, green, or blue) must be provided.")
  }

  first_ch <- non_null[[1]]
  d <- if (is.array(first_ch) || is.matrix(first_ch)) dim(first_ch)[1:2] else NULL
  if (is.null(d)) {
    cli::cli_abort("Channels must be 2D matrices or arrays.")
  }

  w <- d[1]
  h <- d[2]

  # Check if raw mode should be used
  is_raw <- any(sapply(non_null, is.raw)) || any(sapply(non_null, function(x) storage.mode(x) == "raw"))

  # Helper to process each channel
  get_channel <- function(ch) {
    if (is.null(ch)) {
      if (is_raw) return(array(as.raw(0), dim = c(w, h)))
      else return(array(0.0, dim = c(w, h)))
    }
    dat <- image_data(ch)
    if (length(dim(dat)) == 3) dat <- dat[,,1]

    if (is_raw) {
      if (is.raw(dat)) return(dat)
      if (is.numeric(dat)) {
        mx <- suppressWarnings(max(dat, na.rm = TRUE))
        if (is.finite(mx) && mx <= 1.0 && mx >= 0 && any(dat > 0 & dat < 1)) {
          return(array(as.raw(round(dat * 255)), dim = c(w, h)))
        } else {
          return(array(as.raw(pmin(pmax(round(dat), 0), 255)), dim = c(w, h)))
        }
      }
      return(array(as.raw(dat), dim = c(w, h)))
    } else {
      if (is.raw(dat)) return(array(as.numeric(dat) / 255, dim = c(w, h)))
      return(array(as.numeric(dat), dim = c(w, h)))
    }
  }

  r_layer <- get_channel(red)
  g_layer <- get_channel(green)
  b_layer <- get_channel(blue)

  arr <- array(c(r_layer, g_layer, b_layer), dim = c(w, h, 3))
  as_image(arr, colormode = "Color")
}

#' @rdname rgb_image
#' @export
rgbImage <- rgb_image
