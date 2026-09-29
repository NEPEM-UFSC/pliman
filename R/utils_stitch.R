#' Native Planar Image Stitching and Panorama
#'
#' @description
#' `image_stitch()` merges multiple overlapping planar images (e.g. root scanner tiles,
#' oversized plant leaves, bench scans, or drone orthomosaic patches) into a single
#' seamless composite without requiring OpenCV, Python, or external software.
#'
#' The native pipeline performs:
#' 1. **Multi-Scale Normalized Cross-Correlation (NCC)**: Coarse-to-fine registration
#'    finding the exact sub-pixel horizontal and vertical shifts between adjacent tiles.
#' 2. **Cosine-Smooth Seam Feathering**: Eliminates visible seam lines and exposure steps
#'    by applying non-linear cosine weighted transitions across the overlap zone.
#' 3. **Dynamic Canvas Expansion**: Recursively chains two or more images into a single
#'    composite image array.
#'
#' @param images A list of `Image` objects or a character vector of file paths.
#' @param direction Stitching direction: `"horizontal"` (left-to-right), `"vertical"`
#'   (top-to-bottom), or `"auto"` (automatically chosen based on best correlation peak).
#' @param blend Logical. If `TRUE` (default), applies smooth seam feathering across
#'   overlapping image bands. If `FALSE`, places images directly with a sharp cut.
#' @param overlap_hint Expected fractional overlap between adjacent images (default: 0.2 = 20%).
#' @param plot Logical. If `TRUE`, plots the stitched result.
#' @param return_details Logical. If `TRUE`, returns a named list with the stitched `image`
#'   and registration diagnostics (shifts, NCC scores, overlap dimensions). If `FALSE`
#'   (default), returns the stitched `Image` directly.
#' @param verbose Logical. If `TRUE` (default), displays progress messages.
#'
#' @return An `Image` object (or a named list if `return_details = TRUE`).
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Stitch two overlapping leaf scans
#' imgs <- list(
#'   image_read("leaf_left.jpg"),
#'   image_read("leaf_right.jpg")
#' )
#' stitched <- image_stitch(imgs, direction = "horizontal")
#' plot(stitched)
#' }
image_stitch <- function(images,
                         direction = c("horizontal", "vertical", "auto"),
                         blend = TRUE,
                         overlap_hint = 0.2,
                         plot = FALSE,
                         return_details = FALSE,
                         verbose = TRUE) {

  direction <- match.arg(direction)

  # Validate input
  if (is.character(images)) {
    images <- lapply(images, image_import)
  }

  if (!is.list(images) || length(images) < 2) {
    cli::cli_abort("At least two images are required for {.fn image_stitch}.")
  }

  n <- length(images)
  if (verbose) {
    cli::cli_rule(left = cli::col_blue("Native Planar Image Stitching"))
    cli::cli_alert_info("Stitching {n} images along {direction} axis...")
  }

  curr_img <- images[[1]]
  shifts <- list()

  for (i in 2:n) {
    next_img <- images[[i]]

    if (verbose) {
      cli::cli_progress_step("Aligning image {i}/{n}...")
    }

    res <- stitch_pair_cpp(
      img1 = curr_img,
      img2 = next_img,
      direction = direction,
      blend = blend,
      overlap_hint = overlap_hint
    )

    shifts[[i - 1]] <- list(
      pair = c(i - 1, i),
      dx = res$dx,
      dy = res$dy,
      ncc = res$ncc,
      direction = res$direction,
      overlap_w = res$overlap_width,
      overlap_h = res$overlap_height
    )

    curr_img <- res$image
  }

  if (verbose) {
    cli::cli_alert_success("Stitching completed successfully! Canvas dimensions: {dim(curr_img)[1]} x {dim(curr_img)[2]}")
  }

  if (isTRUE(plot)) {
    plot(curr_img)
  }

  if (return_details) {
    return(list(
      image = curr_img,
      diagnostics = shifts
    ))
  } else {
    return(curr_img)
  }
}
