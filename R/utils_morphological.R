# ===========================================================================
# Internal helpers (not exported)
# ===========================================================================

# Gaussian blur helper: wraps the C++ separable blur and preserves
help_gblur <- function(img, sigma) {
  mat <- image_data(img)
  if (is.null(dim(mat))) {
    stop("help_gblur: input must have dim attribute (matrix or array).")
  }
  res <- gaussian_blur_cpp(mat, sigma)
  if (is_image(img)) {
    return(as_image(res))
  }
  res
}

# Validates and extracts the 2D logical matrix from a binary image.
# Generates informative errors via {cli} if the input is invalid.
# Color (3D) images are collapsed channel by channel with OR.
.check_binary <- function(img, arg = "img") {
  if (!requireNamespace("cli", quietly = TRUE))
    stop("{cli} package is required. Install it with: install.packages('cli')")

  if (!inherits(img, "Image") && !inherits(img, "image") && !is.matrix(img) && !is.array(img))
    cli::cli_abort(c(
      "{.arg {arg}} must be an {.cls image} object, a matrix, or an array.",
      "x" = "Provided object has class {.cls {class(img)}}."
    ))

  d <- image_data(img)

  if (storage.mode(d) != "logical") {
    if (is.raw(d) || is.numeric(d)) {
      d <- d > 0
    } else {
      cli::cli_abort(c(
        "{.arg {arg}} must be a binary image ({.val logical}).",
        "x" = "{.code storage.mode({arg})} is {.val {storage.mode(d)}}, not {.val logical}.",
        "i" = "Binarize the image first, for example: {.code img > 0.5}."
      ))
    }
  }

  if (length(dim(d)) == 3L) {
    out <- d[, , 1L]
    for (k in seq_len(dim(d)[3L])[-1L]) out <- out | d[, , k]
    d <- out
  }
  d <- unclass(d)
  dim(d) <- dim(d)[1:2]
  attributes(d) <- list(dim = dim(d))
  d
}

# Converts a matrix back to an image object
.mat_to_img <- function(mat) {
  as_image(mat, colormode = "Grayscale")
}


# ===========================================================================
# Public functions
# ===========================================================================

#' Morphological erosion of a binary image
#'
#' @description
#' Performs morphological erosion of a binary image.
#'
#' @param img An `image` object or a list of `image` objects. Must be binary
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
#' @return A binary `image` object or a list of `image` objects.
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
#' @param img An `image` object or a list of `image` objects. Must be binary
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
#' @return A binary `image` object or a list of `image` objects.
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
#' @param img An `image` object or a list of `image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param max_size Maximum size of holes to fill (in pixels). Default is `NULL` (fills all holes).
#' @param min_neck_dist Minimum neck distance threshold for hole filling. Default is `0.0`.
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `image` object or a list of `image` objects.
#' @name utils_transform
#' @export
image_fill_hull <- function(img,
                            max_size = NULL,
                            min_neck_dist = 0.0,
                            parallel = FALSE,
                            workers = NULL,
                            verbose = TRUE,
                            plot = FALSE) {

  m_sz <- if (is.null(max_size)) -1.0 else as.numeric(max_size)
  n_dist <- if (is.null(min_neck_dist)) 0.0 else as.numeric(min_neck_dist)

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
      return(.mat_to_img(fill_holes_cpp(mat, max_size = m_sz, min_neck_dist = n_dist)))
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
    img <- .mat_to_img(fill_holes_cpp(mat, max_size = m_sz, min_neck_dist = n_dist))

    if (isTRUE(plot)) {
      plot(img)
    }

    invisible(img)
  }
}

#' Remove Specular Light Reflections
#'
#' Removes specular light reflections from grayscale index images prior to thresholding.
#'
#' @param img Grayscale `Image` object or numeric matrix.
#' @param radius Radius of structuring window for specular peak removal (default = 5).
#' @param has_white_bg Logical indicating whether the background is white/light (default = TRUE).
#' @return A filtered `Image` object or matrix.
#' @export
image_remove_reflection <- function(img, radius = 5, has_white_bg = TRUE) {
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    return(lapply(img, function(x) image_remove_reflection(x, radius = radius, has_white_bg = has_white_bg)))
  }
  is_img <- inherits(img, "Image") || inherits(img, "image")
  mat <- if (is_img) image_data(img) else img
  if (length(dim(mat)) == 3L) {
    mat <- mat[,,1L]
  }
  res_mat <- remove_reflection_cpp(mat, radius = as.integer(radius), has_white_bg = isTRUE(has_white_bg))
  if (is_img) {
    return(as_image(res_mat, storage = "double"))
  } else {
    return(res_mat)
  }
}

#' Morphological opening of a binary image
#'
#' @description
#' Performs morphological opening of a binary image (erosion followed by dilation).
#' When `return_exact = TRUE`, morphological reconstruction via watershed-partitioned seed matching is performed.
#' This isolates noise touching large objects into separate watershed labels, identifies surviving labels using the un-dilated eroded marker,
#' and restores 100% of the exact original boundaries and shapes of surviving objects without rounding corners or eroding tips.
#'
#' @param img An `image` object or a list of `image` objects. Must be binary
#'   (`storage.mode = "logical"`).
#' @param size The size of the structuring element (radius in pixels). Default is `NULL`,
#'   which automatically determines the size based on image dimensions.
#' @param shape A character string specifying the shape of the structuring element.
#'   Either `"disc"` (default) or `"square"`.
#' @param external Logical. If `TRUE`, erosion does not consume pixels inside
#'   internal holes (passed to `image_erode()`).
#' @param return_exact Logical. If `TRUE`, returns the exact original pixels and
#'   boundaries of objects in `img` that survived the opening. Uses internal
#'   watershed partitioning to isolate and eliminate noise touching large
#'   objects, preserving 100% of sharp corners and geometric features of
#'   surviving objects.
#' @param parallel Logical. If `TRUE`, processes the list of images in parallel.
#' @param workers Number of parallel workers. If `NULL`, defaults to 40% of
#'   available cores.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param plot Logical. If `TRUE`, plots the resulting image(s).
#'
#' @return A binary `image` object or a list of `image` objects.
#' @name utils_transform
#' @export
image_opening <- function(img,
                          size = NULL,
                          shape = "disc",
                          external = FALSE,
                          return_exact = FALSE,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE) {
  .one_opening <- function(im) {
    if (isTRUE(return_exact)) {
      orig_mask <- im > 0.5

      # 1. Erosion to obtain core marker seeds (un-dilated to prevent touching attached noise)
      eroded <- image_erode(img = im, size = size, shape = shape,
                            external = external, parallel = FALSE,
                            verbose = FALSE, plot = FALSE)
      eroded_mask <- image_data(eroded) > 0.5

      # 2. Watershed segmentation on original mask to isolate attached noise into separate label IDs
      labeled <- image_watershed(orig_mask)
      lab_mat <- image_data(labeled)

      # 3. Find surviving label IDs (watershed segments that intersect eroded_mask)
      surviving_ids <- unique(lab_mat[eroded_mask])
      surviving_ids <- surviving_ids[surviving_ids > 0]

      # 4. Reconstruct exact original mask for surviving watershed objects
      if (length(surviving_ids) == 0L) {
        exact_mat <- matrix(FALSE, nrow = nrow(lab_mat), ncol = ncol(lab_mat))
      } else {
        exact_mat <- lab_mat %in% surviving_ids
        dim(exact_mat) <- dim(lab_mat)
      }

      out_img <- as_image(exact_mat, storage = "raw")
      if (isTRUE(plot)) plot(out_img)
      return(out_img)
    } else {
      if (verbose && is.list(img)) {
        cli::cli_alert_info("Starting morphological opening (erosion followed by dilation).")
      }
      eroded <- image_erode(img = im, size = size, shape = shape,
                            external = external, parallel = FALSE,
                            verbose = verbose, plot = FALSE)
      opened <- image_dilate(img = eroded, size = size, shape = shape,
                             parallel = FALSE, verbose = verbose, plot = plot)
      if (isTRUE(plot)) plot(opened)
      return(opened)
    }
  }

  if (is.list(img)) {
    res <- lapply(img, .one_opening)
    return(res)
  } else {
    return(.one_opening(img))
  }
}

#' Remove Small Binary Objects
#'
#' @description
#' Removes small noise objects from a binary image without altering or eroding
#' the boundaries of larger objects. Filtering can be based on an absolute pixel count (`min_size`)
#' or a relative fraction of the mean object area (`rel_size`, e.g. `0.1` for objects smaller than 10% of average object area).
#'
#' @param img A binary `image` object or a list of binary `image` objects.
#' @param min_size Minimum area (in pixels) for an object to be retained.
#' @param rel_size Relative threshold as a fraction of the mean object area (e.g. `0.1` for 10% of mean area).
#' @param watershed Logical. Default `TRUE`. Applies watershed segmentation before measuring component sizes,
#'   ensuring small noise projections touching large objects are split and removed.
#' @param plot Logical. If `TRUE`, plots the filtered image.
#'
#' @return A binary `image` object or a list of binary `image` objects.
#' @name utils_transform
#' @export
image_remove_small <- function(img,
                               min_size = NULL,
                               rel_size = NULL,
                               watershed = TRUE,
                               plot = FALSE) {
  .one_remove_small <- function(im) {
    mask <- im > 0.5
    labeled <- if (isTRUE(watershed)) image_watershed(mask) else image_bwlabel(mask)
    lab_mat <- image_data(labeled)

    if (max(lab_mat) == 0L) {
      if (isTRUE(plot)) plot(im)
      return(im)
    }

    tbl <- table(lab_mat[lab_mat > 0])
    areas <- as.numeric(tbl)
    ids <- as.integer(names(tbl))

    cutoff <- 0
    if (!is.null(min_size) && is.numeric(min_size)) {
      cutoff <- max(cutoff, min_size)
    }
    if (!is.null(rel_size) && is.numeric(rel_size)) {
      mean_area <- mean(areas)
      rel_cutoff <- rel_size * mean_area
      cutoff <- max(cutoff, rel_cutoff)
    }

    if (cutoff <= 0) {
      if (isTRUE(plot)) plot(im)
      return(im)
    }

    keep_ids <- ids[areas >= cutoff]

    if (length(keep_ids) == 0L) {
      res_mat <- matrix(FALSE, nrow = nrow(lab_mat), ncol = ncol(lab_mat))
    } else {
      res_mat <- lab_mat %in% keep_ids
      dim(res_mat) <- dim(lab_mat)
    }

    res_img <- as_image(res_mat, storage = "raw")
    if (isTRUE(plot)) plot(res_img)
    return(res_img)
  }

  if (is.list(img)) {
    res <- lapply(img, .one_remove_small)
    return(res)
  } else {
    return(.one_remove_small(img))
  }
}

#' Morphological closing of a binary image
#'
#' @description
#' Performs morphological closing of a binary image. Closing is defined as a
#' dilation followed by an erosion. It is useful for filling small holes or
#' connecting nearby objects.
#'
#' @param img An `image` object or a list of `image` objects. Must be binary
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
#' @return A binary `image` object or a list of `image` objects.
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
#' @param labels A labeled `image` object or a matrix containing integer labels. Typically
#'   the result of [image_watershed()] or similar segmentation.
#'
#' @return A list where each element is a matrix of (x, y) coordinates for the
#'   contour of the respective object.
#' @name contour
#' @export
#' @examples
#' library(pliman)
#' # generate a binary image
#' img <- image_pliman("soybean_touch.jpg")
#' bin <- image_binary(img, "B")[[1]]
#' # watershed segmentation
#' wat <- image_watershed(bin)
#' # extract contours
#' cont <- contour(wat)
#' plot(img)
#' plot_contour(cont, col = "red", lwd = 2)
contour <- function(labels) {
  if (inherits(labels, c("Image", "image")) || is.matrix(labels) || is.array(labels)) {
    if (storage.mode(labels) != "integer") {
      storage.mode(labels) <- "integer"
    }
  } else if (is.logical(labels) || is.numeric(labels) || is.raw(labels)) {
    if (storage.mode(labels) != "integer") {
      storage.mode(labels) <- "integer"
    }
  } else {
    cli::cli_abort("The input must be a labeled matrix or Image object.")
  }
  res <- extract_contours_cpp(labels)
  return(res)
}
