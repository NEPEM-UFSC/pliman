# nova função


#' Intelligent Morphological Watershed (Spanning Forest + Union-Find)
#'
#' This function performs a morphological watershed segmentation on a binary
#' image. It implements a highly optimized Spanning Forest algorithm coupled
#' with a Disjoint-Set (Union-Find) data structure. This approach eliminates the
#' need for slow iterative scans ($O(N^2)$ complexity), allowing the algorithm
#' to execute in milliseconds even when resolving tens of thousands of distinct
#' objects. It also features a mechanism to merge highly connected objects,
#' effectively preventing over-segmentation.
#'
#' @param img           A binary `image` object or a matrix/array.
#' @param sensitivity   Text ("low", "medium", "high", "extreme") or Number
#'   (e.g. 2.5). Automatically defines the best parameters for optimal
#'   precision. (Overrides manual adjustment of tolerance).
#' @param tolerance     Manual Fine-Tuning. If provided, overrides `sensitivity`
#'   and directly controls the merging threshold for connected components.
#'   Values < 1.0 are interpreted as relative tolerance (percentage of the peak
#'   height). Values >= 1.0 are interpreted as absolute tolerance in pixels
#'   (linear Euclidean distance).
#' @param ext           Pixel neighborhood radius (default = 1). Specifies the
#'   extension radius used during the regional maxima calculation and watershed
#'   flooding.
#' @param min_size Minimum area (in pixels) for an object segment retained after watershed.
#' @param rel_size Relative threshold as a fraction of the mean segment area (e.g. `0.1` for 10% of mean area).
#'
#' @return A grayscale `image` object containing the labels.
#' @importFrom methods .hasSlot slot
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   bin <- image_pliman("soybean_touch.jpg") |> image_binary(index = "B")
#'   w_med <- image_watershed(bin$B)
#' }
image_watershed <- function(img,
                            sensitivity = "medium",
                            tolerance = NULL,
                            ext = 1,
                            min_size = NULL,
                            rel_size = NULL,
                            fill_hull = FALSE,
                            max_size = NULL) {
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    if (inherits(img, c("binary_list", "segment_list", "index_list", "img_mat_list", "palette_list"))) {
      img <- lapply(img, function(x) x[[1L]])
    }
    res <- lapply(img, function(x) image_watershed(x, sensitivity = sensitivity, tolerance = tolerance, ext = ext, min_size = min_size, rel_size = rel_size, fill_hull = fill_hull, max_size = max_size))
    return(res)
  }

  if (!requireNamespace("cli", quietly = TRUE))
    stop("{cli} package is required.")

  if (!inherits(img, "Image") && !inherits(img, "image") && !is.matrix(img) && !is.array(img))
    cli::cli_abort(c("{.arg img} must be an {.cls image} object, a matrix or an array."))

  mat <- image_data(img)

  if (length(dim(mat)) == 3L) {
    d2 <- mat[, , 1L]
    for (k in seq_len(dim(mat)[3L])[-1L]) d2 <- d2 | mat[, , k]
    mat <- d2
  }

  if (isTRUE(fill_hull) || !is.null(max_size)) {
    m_sz <- if (is.null(max_size)) -1.0 else as.numeric(max_size)
    mat <- fill_holes_cpp(mat > 0, max_size = m_sz)
  }

  # Intelligent Mapping
  if (is.null(tolerance)) {
    if (is.character(sensitivity)) {
      sens_str <- tolower(trimws(sensitivity[1]))
      lvl <- switch(sens_str,
                    "low"     = 1.0,
                    "medium"  = 2.0,
                    "high"    = 3.0,
                    "extreme" = 4.0,
                    cli::cli_abort("Invalid sensitivity '{sens_str}'. Use: 'low', 'medium', 'high', 'extreme'."))
    } else if (is.numeric(sensitivity)) {
      lvl <- as.double(sensitivity[1])
    } else {
      cli::cli_abort("Sensitivity must be text or numeric.")
    }

    # Scale-Invariant Relative Tolerance (Percentage drop)
    tolerancia <- 0.20 / (2 ^ (lvl - 1.0))
  } else {
    tolerancia <- as.double(tolerance)
  }

  res <- watershed_cpp(mat, tolerance = tolerancia, ext = as.integer(ext))

  # Post-Watershed Size Filtering
  if (!is.null(min_size) || !is.null(rel_size)) {
    tbl <- table(res[res > 0])
    if (length(tbl) > 0) {
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

      if (cutoff > 0) {
        keep_ids <- ids[areas >= cutoff]
        if (length(keep_ids) < length(ids)) {
          filter_labels_cpp(res, keep_ids)
        }
      }
    }
  }

  return(as_image(res, colormode = "Grayscale", storage = "integer"))
}

#' Label connected components in a binary image
#'
#' This function performs Connected Component Labeling (CCL) on a binary image.
#' It implements a highly optimized 2-pass Union-Find algorithm, providing
#' lightning-fast execution even for massive images.
#'
#' @param binary A binary `image` object or a logical matrix/array.
#'
#' @return A grayscale `image` object containing the labels, where each isolated
#' object is assigned a unique integer value. Background pixels are 0.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   bin <- image_pliman("soybean_touch.jpg") |> image_binary(index = "B")
#'   labels <- image_bwlabel(bin$B)
#' }
image_bwlabel <- function(binary) {
  if (!is.logical(binary)) {
    binary <- binary > 0
  }

  if (storage.mode(binary) != "logical") {
    cli::cli_abort("The input must be a logical matrix or binary Image.")
  }

  res <- bwlabel_cpp(binary)

  return(as_image(res, colormode = "Grayscale", storage = "integer"))
}

#' Colorize Labeled Images
#'
#' Colorizes a labeled image by allocating a different color to each object.
#' Background pixels (value 0) are colorized with black.
#'
#' @param labels A labeled image (an `image` object or a matrix containing integer labels).
#'
#' @return A color `image` object.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   bin <- image_pliman("soybean_touch.jpg") |> image_binary(index = "B")
#'   lbl <- image_bwlabel(bin$B)
#'   color_lbl <- image_color_labels(lbl)
#'   plot(color_lbl)
#' }
image_color_labels <- function(labels) {
  if (inherits(labels, "Image") || inherits(labels, "image")) {
    mat <- unclass(labels)
  } else {
    mat <- labels
  }

  res <- color_labels_cpp(mat)
  return(as_image(res, colormode = "Color"))
}

#' Distance map transform
#'
#' Computes the distance map transform of a binary image. The distance map is a
#' matrix which contains for each pixel the distance to its nearest background
#' pixel.
#'
#' @param binary A binary `image` object or logical matrix.
#' @param normalize Logical argument indicating whether the distance map should
#'   be normalized to the 0-1 range. Defaults to `FALSE`.
#'
#' @return An `image` object or an array, with pixels containing the distances
#'   to the nearest background points
#' @export
#' @examples
#' if (interactive()) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' binary <- image_binary(img, "B")[[1]]
#' wts <- dist_transform(binary)
#' range(wts)
#'}
dist_transform <- function(binary, normalize = FALSE){
  if (is.raw(binary) || is.numeric(binary)) {
    binary <- binary > 0
  }

  if (storage.mode(binary) != "logical") {
    cli::cli_abort("The input must be a logical matrix or binary Image.")
  }

  res <- help_dist_transform(binary)

  if (isTRUE(normalize)) {
    max_val <- max(res, na.rm = TRUE)
    if (max_val > 0) {
      res <- res / max_val
    }
  }

  if (inherits(binary, "Image") || inherits(binary, "image")) {
    return(as_image(res, colormode = "Grayscale"))
  }
  return(res)
}


#' Labels objects
#'
#' All pixels for each connected set of foreground (non-zero) pixels in x are
#' set to an unique increasing integer, starting from 1. Hence, max(x) gives the
#' number of connected objects in `x`.
#' @inheritParams image_binary
#' @inheritParams analyze_objects
#' @return A list with the same length of `img` containing the labeled objects.
#' @export
#'
#' @examples
#' if (interactive()) {
#' img <- image_pliman("soybean_touch.jpg")
#' # segment the objects using the "B" (blue) band.
#' object_label(img, index = "B")
#' object_label(img, index = "B", watershed = TRUE)
#' }

object_label <- function(img,
                         index = "B",
                         invert = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         k = 0.1,
                         windowsize = NULL,
                         opening =  FALSE,
                         closing = FALSE,
                         filter = FALSE,
                         erode = FALSE,
                         dilate = FALSE,
                         filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                         watershed = FALSE,
                         tolerance = NULL,
                         extension = NULL,
                         object_size = "medium",
                         plot = TRUE,
                         ncol = NULL,
                         nrow = NULL,
                         verbose = TRUE){

  img2 <- image_binary(img,
                       index = index,
                       invert = invert,
                       fill_hull = fill_hull,
                       threshold = threshold,
                       k = k,
                       windowsize = windowsize,
                       opening =  opening,
                       closing = closing,
                       filter = filter,
                       erode = erode,
                       dilate = dilate,
                       filter_order = filter_order,
                       resize = FALSE,
                       plot = FALSE)
  labels <- list()
  img2_len <- length(img2)
  for (i in 1:length(img2)){
    if(img2_len > 1){
      tmp <- img2[[i]][[1]]
    } else{
      tmp <- img2[[i]]
    }
    if(isTRUE(watershed)){
      labels[[i]] <- image_watershed(tmp, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
    } else{
      labels[[i]] <- image_bwlabel(tmp)
    }
  }
  if(plot == TRUE){
    num_plots <- length(labels)
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
    op <- par(mfrow = c(nrow, ncol))
    on.exit(par(op))
    index <- names(labels)
    for(i in 1:length(labels)){
      plot(image_color_labels(labels[[i]]))
      if(verbose == TRUE){
        dim <- image_dimension(labels[[i]], verbose = FALSE)
        text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
      }
    }
  }
  invisible(labels)
}

#' Calculate Otsu's threshold
#'
#' Given a numeric vector with the pixel's intensities, returns the threshold
#' value based on Otsu's method, which minimizes the combined intra-class
#' variance
#'
#' @param values A numeric vector with the pixel values.
#'
#' @return
#' A double (threshold value).
#'
#' @references Otsu, N. 1979. Threshold selection method from gray-level
#'   histograms. IEEE Trans Syst Man Cybern SMC-9(1): 62–66. doi:
#'   \doi{10.1109/tsmc.1979.4310076}

#' @export
#'
#' @examples
#' if (interactive()) {
#' img <- image_pliman("soybean_touch.jpg")
#' thresh <- otsu(image_data(img[,,3]))
#' plot(img[,,3] < thresh)
#' }
#'
otsu <- function(values){
  help_otsu(values)
}

