# ===========================================================================
# Internal helpers (not exported)
# ===========================================================================

# Validates and extracts the 2D logical matrix from a binary image.
# Generates informative errors via {cli} if the input is invalid.
# Color (3D) images are collapsed channel by channel with OR.
.check_binary <- function(img, arg = "img") {
  if (!requireNamespace("cli", quietly = TRUE))
    stop("{cli} package is required. Install it with: install.packages('cli')")

  if (!inherits(img, "Image") && !is.matrix(img) && !is.array(img))
    cli::cli_abort(c(
      "{.arg {arg}} must be an {.cls Image} object, a matrix, or an array.",
      "x" = "Provided object has class {.cls {class(img)}}."
    ))

  # Extract data array without EBImage dependency
  if (isS4(img) && inherits(img, "Image") && .hasSlot(img, ".Data")) {
    d <- img@.Data
  } else {
    d <- as.array(img)
  }

  if (storage.mode(d) != "logical")
    cli::cli_abort(c(
      "{.arg {arg}} must be a binary image ({.val logical}).",
      "x" = "{.code storage.mode({arg})} is {.val {storage.mode(d)}}, not {.val logical}.",
      "i" = "Binarize the image first, for example: {.code img > 0.5}."
    ))

  if (length(dim(d)) == 3L) {
    out <- d[, , 1L]
    for (k in seq_len(dim(d)[3L])[-1L]) out <- out | d[, , k]
    return(out)
  }
  d
}

# Converts a matrix back to an image object
.mat_to_img <- function(mat) {
  EBImage::Image(mat, colormode = "Grayscale")
}


# ===========================================================================
# Public functions
# ===========================================================================

#' Morphological erosion of a binary image
#'
#' @description
#' Performs morphological erosion of a binary image.
#'
#' @param img An `Image` object or a list of `Image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param size The size of the structuring element (radius in pixels). Default
#'   is `NULL`, which automatically determines the size based on image
#'   dimensions.
#' @param shape A character string specifying the shape of the structuring element.
#'   Either `"disc"` (default) or `"square"`.
#' @param external Logical. If `TRUE`, erosion does not consume pixels inside
#'   internal holes. If `FALSE` (default), erosion acts on all object boundaries
#'   (essential for separating touching grains).
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `Image` object or a list of `Image` objects.
#' @name utils_transform
#' @export
image_erode <- function(img,
                        size = NULL,
                        shape = "disc",
                        external = FALSE,
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
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    erode_image <- function(im) {
      d <- dim(im)
      s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      s <- ifelse(s == 0, 2, s)
      mat <- .check_binary(im)
      if (isTRUE(external)) {
        res <- erode_external_cpp(mat, raio = as.integer(s), forma = shape)
      } else {
        res <- erode_cpp(mat, raio = as.integer(s), forma = shape)
      }
      return(.mat_to_img(res))
    }

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
    d <- dim(img)
    s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
    s <- ifelse(s == 0, 2, s)
    mat <- .check_binary(img)
    if (isTRUE(external)) {
      res <- erode_external_cpp(mat, raio = as.integer(s), forma = shape)
    } else {
      res <- erode_cpp(mat, raio = as.integer(s), forma = shape)
    }
    img <- .mat_to_img(res)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' Morphological dilation of a binary image
#'
#' @description
#' Performs morphological dilation of a binary image.
#'
#' @param img An `Image` object or a list of `Image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param size The size of the structuring element (radius in pixels). Default is `NULL`,
#'   which automatically determines the size based on image dimensions.
#' @param shape A character string specifying the shape of the structuring element.
#'   Either `"disc"` (default) or `"square"`.
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `Image` object or a list of `Image` objects.
#' @name utils_transform
#' @export
image_dilate <- function(img,
                         size = NULL,
                         shape = "disc",
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
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    dilate_image <- function(im) {
      d <- dim(im)
      s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
      s <- ifelse(s == 0, 2, s)
      mat <- .check_binary(im)
      res <- dilate_cpp(mat, raio = as.integer(s), forma = shape)
      return(.mat_to_img(res))
    }

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
    d <- dim(img)
    s <- ifelse(is.null(size), round(d[[1]] * d[[2]] / 1e06 * 5, 0), size)
    s <- ifelse(s == 0, 2, s)
    mat <- .check_binary(img)
    res <- dilate_cpp(mat, raio = as.integer(s), forma = shape)
    img <- .mat_to_img(res)

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' Fill internal holes in a binary image
#'
#' @description
#' Fills internal holes in a binary image.
#'
#' @param img An `Image` object or a list of `Image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `Image` object or a list of `Image` objects.
#' @name utils_transform
#' @export
image_fill_hull <- function(img,
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
      cli::cli_abort("All elements in the list must be of class {.cls Image}.")
    }

    fill_image <- function(im) {
      mat <- .check_binary(im)
      return(.mat_to_img(fill_holes_cpp(mat)))
    }

    if (parallel) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Filling holes in {length(img)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images in parallel...",
          msg_done   = "Hole filling complete.",
          msg_failed = "Parallel hole filling failed."
        )
      }

      res <- mirai::mirai_map(
        .x = img,
        .f = fill_image
      )[.progress]

    } else {
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Filling holes in {length(img)} images sequentially"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%OS0')}")
        )
        cli::cli_progress_step(
          msg        = "Processing images...",
          msg_done   = "Hole filling complete.",
          msg_failed = "Sequential hole filling failed."
        )
      }

      res <- lapply(img, fill_image)
    }

    if (isTRUE(plot)) {
      for (r in res) plot(r)
    }

    invisible(res)

  } else {
    mat <- .check_binary(img)
    img <- .mat_to_img(fill_holes_cpp(mat))

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' Morphological opening of a binary image
#'
#' @description
#' Performs morphological opening of a binary image. Opening is defined as an
#' erosion followed by a dilation. It is useful for removing small objects or
#' noise from the image.
#'
#' @param img An `Image` object or a list of `Image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param size The size of the structuring element (radius in pixels). Default is `NULL`,
#'   which automatically determines the size based on image dimensions.
#' @param shape A character string specifying the shape of the structuring element.
#'   Either `"disc"` (default) or `"square"`.
#' @param external Logical. If `TRUE`, erosion does not consume pixels inside
#'   internal holes (passed to `image_erode()`).
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `Image` object or a list of `Image` objects.
#' @name utils_transform
#' @export
image_opening <- function(img,
                          size = NULL,
                          shape = "disc",
                          external = FALSE,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE) {
  # Opening is an erosion followed by a dilation
  if (verbose && is.list(img)) {
    cli::cli_alert_info("Starting morphological opening (erosion followed by dilation).")
  }
  img <- image_erode(img = img, size = size, shape = shape,
                     external = external, parallel = parallel, workers = workers,
                     verbose = verbose, plot = FALSE)
  img <- image_dilate(img = img, size = size, shape = shape,
                      parallel = parallel, workers = workers,
                      verbose = verbose, plot = plot)
  invisible(img)
}

#' Morphological closing of a binary image
#'
#' @description
#' Performs morphological closing of a binary image. Closing is defined as a
#' dilation followed by an erosion. It is useful for filling small holes or
#' connecting nearby objects.
#'
#' @param img An `Image` object or a list of `Image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param size The size of the structuring element (radius in pixels). Default is `NULL`,
#'   which automatically determines the size based on image dimensions.
#' @param shape A character string specifying the shape of the structuring element.
#'   Either `"disc"` (default) or `"square"`.
#' @param external Logical. If `TRUE`, erosion does not consume pixels inside
#'   internal holes (passed to `image_erode()`).
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `Image` object or a list of `Image` objects.
#' @name utils_transform
#' @export
image_closing <- function(img,
                          size = NULL,
                          shape = "disc",
                          external = FALSE,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE) {
  # Closing is a dilation followed by an erosion
  if (verbose && is.list(img)) {
    cli::cli_alert_info("Starting morphological closing (dilation followed by erosion).")
  }
  img <- image_dilate(img = img, size = size, shape = shape,
                      parallel = parallel, workers = workers,
                      verbose = verbose, plot = FALSE)
  img <- image_erode(img = img, size = size, shape = shape,
                     external = external, parallel = parallel, workers = workers,
                     verbose = verbose, plot = plot)
  invisible(img)
}

#' Extract object contours from a labeled image
#'
#' @description
#' Extracts the external contour of each object in a labeled matrix (e.g., from
#' watershed segmentation). It returns a list of coordinates for each object.
#'
#' @param labels A labeled `Image` object or a matrix containing integer labels. Typically
#'   the result of [image_watershed()] or similar segmentation.
#'
#' @return A list where each element is a matrix of (x, y) coordinates for the
#'   contour of the respective object.
#' @name contour
#' @export
#' @examples
#' library(pliman)
#' # generate a binary image
#' img <- image_pliman("soybean_touching.jpg")
#' bin <- image_binary(img, "B")[[1]]
#' # watershed segmentation
#' wat <- image_watershed(bin)
#' # extract contours
#' cont <- contour(wat)
#' plot(img)
#' plot_contour(cont, col = "red", lwd = 2)
contour <- function(labels) {
  if(storage.mode(labels) != "integer"){
    cli::cli_abort("The input must be a labeled matrix or Image object.")
  }
  res <- extract_contours_cpp(labels)
  return(res)
}
