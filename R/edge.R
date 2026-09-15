#' Object edges
#'
#' Applies the Sobel-Feldman Operator to detect edges. The operator is based on
#' convolving the image with a small, separable, and integer-valued filter in
#' the horizontal and vertical directions.
#'
#' @param img An image or a list of images of class `image`.
#' @param sigma Gaussian kernel standard deviation used in the gaussian blur.
#' @param threshold The theshold method to be used.  If `threshold = "Otsu"`
#'   (default), a threshold value based on Otsu's method is used to reduce the
#'   grayscale image to a binary image. If any non-numeric value different than
#'   `"Otsu"` is used, an iterative section will allow you to choose the
#'   threshold based on a raster plot showing pixel intensity of the index.
#'   Alternatively, provide a numeric value to be used as the threshold value.
#' @param thinning Logical value indicating whether a thinning procedure should
#'   be applied to the detected edges. See [image_skeleton()]
#' @param plot Logical value indicating whether a plot should be created
#' @return A binary version of `image`.
#' @references Sobel, I., and G. Feldman. 1973. A 3×3 isotropic gradient
#'   operator for image processing. Pattern Classification and Scene Analysis:
#'   271–272.
#' @export
#'
#'
#'
#' @examples
#' if (interactive()) {
#' library(pliman)
#' img <- image_pliman("sev_leaf_nb.jpg", plot = TRUE)
#' object_edge(img)
#' }
#'
object_edge <- function(img,
                        sigma = 1,
                        threshold = "Otsu",
                        thinning = FALSE,
                        plot = TRUE){
  
  gray <- image_index(img,
                      "GRAY",
                      plot = FALSE,
                      verbose = FALSE)[[1]]
  if(!isFALSE(sigma)){
    gray <- help_gblur(gray, sigma)
  }
  edata <- sobel_help(gray)

  if (threshold == "Otsu") {
    threshold <- help_otsu(as.numeric(image_data(edata)))
  }  else {
    if (is.numeric(threshold)) {
      threshold <- threshold
    }
    else {
      pixels <- terra::rast(image_data(image_transpose(edata)))
      terra::plot(pixels, col = custom_palette(),  axes = FALSE, asp = NA)
      threshold <- readline("Selected threshold: ")
    }
  }
  edata <- as_image(edata > threshold)
  if(isTRUE(thinning)){
    edata <- image_thinning(edata, verbose = FALSE, plot = FALSE)
  }
  if(isTRUE(plot)){
    plot(edata)
  }
  invisible(edata)
}


#' Object veins detection and quantification
#'
#' Detects leaf veins using Difference of Gaussians (DoG) bandpass filtering,
#' dynamic resolution-independent boundary erosion, and per-object Otsu adaptive
#' thresholding. Optionally applies Guo-Hall skeleton thinning.
#'
#' @param img An image or a list of images of class `image`.
#' @param labels Optional matrix or `image` object containing object labels
#'   (e.g., obtained from [object_label()]). If `NULL` (default), labels are
#'   automatically computed using [object_label()].
#' @param index The image index used for object segmentation when `labels = NULL` (default `"GRAY"`).
#'   See [image_index()].
#' @param watershed Logical argument indicating whether watershed segmentation should be performed
#'   when `labels = NULL` (default `TRUE`).
#' @param sigma1 Fine Gaussian scale parameter for vein detection (default `0.75`).
#' @param sigma2 Coarse Gaussian scale parameter for vein detection. Defaults to `sigma1 * 5`.
#' @param threshold Fixed vein threshold value in `[0, 1]`. If `< 0` (default `-1`),
#'   adaptive per-object Otsu thresholding is performed.
#' @param channel Channel to extract for vein detection: `0` for green (default,
#'   best for leaves) or `1` for perceptual grayscale.
#' @param erode_size Fixed boundary erosion radius in pixels. If `< 0` (default `-1`),
#'   the erosion radius is dynamically calculated using `rel_erode`.
#' @param rel_erode Relative boundary erosion fraction of `min(width, height)` (default `0.005` = 0.5%).
#'   Trims boundary pixels to eliminate outer edge artifacts independently of image resolution.
#' @param thinning Logical value indicating whether Guo-Hall skeleton thinning should
#'   be applied to detected veins (default `FALSE`). If `TRUE`, measures vein length
#'   skeleton density instead of vein area fraction.
#' @param overlay Logical argument indicating whether detected veins should be overlaid
#'   in color onto the original image (default `TRUE`).
#' @param col Color to use for overlaying detected veins (default `"red"`).
#' @param plot Logical value indicating whether a plot should be created (default `TRUE`).
#'   If `overlay = TRUE`, plots the original image with overlaid veins.
#'
#' @return A list containing:
#'   * `proportion`: A numeric vector with vein proportions for each object ID.
#'   * `vein_map`: An `image` object representing the binary vein map (or thinned skeleton).
#'   * `overlay`: An `image` object showing the original image with veins highlighted in `col` (if `overlay = TRUE`).
#' @export
#'
#' @examples
#' if (interactive()) {
#' library(pliman)
#' img <- image_pliman("sev_leaf_nb.jpg", plot = TRUE)
#' veins <- object_veins(img, overlay = TRUE, col = "red")
#' veins$proportion
#' }
#'
object_veins <- function(img,
                         labels = NULL,
                         index = "GRAY",
                         watershed = TRUE,
                         sigma1 = 0.75,
                         sigma2 = sigma1 * 5,
                         threshold = -1,
                         channel = 0,
                         erode_size = -1,
                         rel_erode = 0.005,
                         thinning = FALSE,
                         overlay = TRUE,
                         col = "red",
                         plot = TRUE) {
  

  .overlay_veins <- function(im, vmap_data, color = "red") {
    col_rgb <- tryCatch(
      grDevices::col2rgb(color)[, 1] / 255,
      error = function(e) c(1, 0, 0)
    )
    res <- im
    mask <- (vmap_data == 1)
    n_ch <- min(3, dim(im)[3])
    for (ch in 1:n_ch) {
      arr <- res[,,ch]
      arr[mask] <- col_rgb[ch]
      res[,,ch] <- arr
    }
    res
  }

  .one_veins <- function(im, lab) {
    if (is.null(lab)) {
      lab <- object_label(im, index = index, watershed = watershed, plot = FALSE, verbose = FALSE, erode = )[[1]] |> image_data() |> as.numeric()
    } else if (inherits(lab, "Image")) {
      lab <- as.numeric(image_data(lab))
    }

    res <- detect_veins_cpp(
      R_sexp      = im[,,1],
      G_sexp      = im[,,2],
      B_sexp      = im[,,3],
      labels_sexp = lab,
      sigma1     = sigma1,
      sigma2     = sigma2,
      threshold  = threshold,
      channel    = channel,
      erode_size = erode_size,
      rel_erode  = rel_erode,
      thinning   = thinning,
      return_map = TRUE
    )

    vmap <- as_image(res$vein_map)
    over_img <- if (isTRUE(overlay)) .overlay_veins(im, res$vein_map, col) else NULL

    if (isTRUE(plot)) {
      if (isTRUE(overlay) && !is.null(over_img)) {
        plot(over_img)
      } else {
        plot(vmap)
      }
    }

    list(
      proportion  = res$proportion,
      vein_map    = vmap,
      overlay     = over_img
    )
  }

  if (is.list(img)) {
    if (!all(sapply(img, inherits, "Image"))) {
      cli::cli_abort("All elements in {.arg img} must be of class {.cls Image}.")
    }
    if (!is.null(labels) && is.list(labels)) {
      res <- lapply(seq_along(img), function(i) .one_veins(img[[i]], labels[[i]]))
    } else {
      res <- lapply(img, function(im) .one_veins(im, labels))
    }
    invisible(res)
  } else {
    res <- .one_veins(img, labels)
    invisible(res)
  }
}
