#' Measures Injury in Images
#'
#' The `measures_injury` function calculates the percentage of injury in images
#' by performing binary segmentation and identifying lesions. It processes
#' either a single image or a batch of images specified by a pattern in a
#' directory.
#'
#' @inheritParams analyze_objects
#' @param dir_original The directory containing the original and processed
#'   images. Defaults to NULL. In this case, the function will search for the
#'   image img in the current working directory.
#'
#' @return A numeric value representing the injury percentage for a single
#'   image, or a data frame with injury percentages for batch processing.
#'
#' @details The function processes each image by reading it, applying binary
#'   segmentation to detect lesions, filling the segmented areas, calculating
#'   the injury percentage, and optionally saving the processed image with
#'   highlighted lesions. In batch mode, it uses the provided pattern to
#'   identify images in the specified directory and can utilize parallel
#'   processing for efficiency.
#'
#'
#' @export
measure_injury <- function(img = NULL,
                           pattern = NULL,
                           index = "GRAY",
                           threshold = "Otsu",
                           invert = FALSE,
                           opening = 5,
                           closing = FALSE,
                           filter = FALSE,
                           erode = FALSE,
                           dilate = FALSE,
                           plot = TRUE,
                           dir_original = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE) {

  # Function to process a single image
  process_image <- function(image_path) {
    # Import image
    if (is.character(image_path)) {
      img <- EBImage::readImage(image_path)
    } else {
      img <- image_path
    }

    # Binary segmentation
    seg <- image_binary(img,
                        index = index,
                        threshold = threshold,
                        opening = opening,
                        closing = closing,
                        erode = erode,
                        dilate = dilate,
                        filter = filter,
                        invert = invert,
                        plot = FALSE)[[1]]

    # Fill segmentation
    segfill <- EBImage::fillHull(seg)
    lesions <- segfill - seg
    ID <- which(lesions == 1)

    # Color lesions
    img@.Data[,,1][ID] <- 165 / 255
    img@.Data[,,2][ID] <- 42 / 255
    img@.Data[,,3][ID] <- 42 / 255

    # Plot if required
    if (plot) {
      plot(img)
    }

    # Calculate injury percentage
    injury_percentage <- (sum(segfill) - sum(seg)) / sum(segfill) * 100

    return(injury_percentage)
  }

  if (is.null(img) & !is.null(pattern)) {
    # Batch processing mode
    if (is.null(dir_original)) {
      dir_original <- "./"
    }
    image_files <- list.files(dir_original, pattern = pattern, full.names = TRUE)

    if (length(image_files) == 0) {
      stop("No images found matching the pattern.")
    }

    if(parallel == TRUE){
      init_time <- Sys.time()
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.3), workers)
      future::plan(future::multisession, workers = nworkers)
      on.exit(future::plan(future::sequential))
      `%dofut%` <- doFuture::`%dofuture%`

      if(verbose == TRUE){
        message("Processing ", length(image_files), " images in multiple sessions (",nworkers, "). Please, wait.")
      }

      results <-
        foreach::foreach(i = seq_along(image_files)) %dofut%{
          process_image(image_files[i])
        }

    } else{
      init_time <- Sys.time()
      pb <- progress(max = length(image_files), style = 4)
      foo <- function(plants, ...){
        if(verbose == TRUE){
          run_progress(pb, ...)
        }
        process_image(plants)
      }
      results <-
        lapply(seq_along(image_files), function(i){
          foo(image_files[i],
              actual = i,
              text = paste("Processing image", file_name(image_files[i])))
        })
    }

    names(results) <- basename(image_files)
    results <- data.frame(do.call(rbind, results))
    results$img <- rownames(results)
    rownames(results) <- NULL
    colnames(results) <- c("injury", "img")
    return(results[, 2:1])
  } else if (!is.null(img)) {
    result <- process_image(img)
    return(result)
    if (verbose) {
      cat("Injury percentage for the image:", result, "%\n")
    }
  } else {
    stop("Either 'img' or 'pattern' must be provided.")
  }
}

