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

  # helper to process one image
  process_image <- function(image_path) {
    img_obj <- if (is.character(image_path)) {
      EBImage::readImage(image_path)
    } else {
      image_path
    }

    seg <- image_binary(
      img_obj, index, threshold,
      opening = opening, closing = closing,
      erode = erode, dilate = dilate,
      filter = filter, invert = invert,
      plot = FALSE
    )[[1]]

    segfill <- EBImage::fillHull(seg)
    lesions <- segfill - seg
    idx <- which(lesions == 1)

    img_obj@.Data[,,1][idx] <- 165/255
    img_obj@.Data[,,2][idx] <-  42/255
    img_obj@.Data[,,3][idx] <-  42/255

    if (plot){
      plot(img_obj)
    }

    pct <- (sum(segfill) - sum(seg)) / sum(segfill) * 100
    pct
  }

  # single-image mode
  if (!is.null(img)) {
    result <- process_image(img)
    if (verbose) {
      cli::cli_alert_info("Injury percentage: {.val {round(result, 2)}}%")
    }
    return(result)
  }

  # batch mode
  if (is.null(pattern)) {
    cli::cli_abort("{.arg pattern} must be provided for batch processing.")
  }
  dir_original <- dir_original %||% "./"
  files <- list.files(dir_original, pattern, full.names = TRUE)

  if (length(files) == 0) {
    cli::cli_abort(c(
      "x" = "No images found matching {.val {pattern}} in {.path {dir_original}}.",
      "i" = "Check your working directory: {.path {getwd()}}"
    ))
  }

  init_time <- Sys.time()

  if (parallel) {
    # setup
    # capture start time
    init_time <- Sys.time()

    # decide number of workers
    nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.3) else workers

    # start mirai daemons
    mirai::daemons(nworkers)
    on.exit(mirai::daemons(0), add = TRUE)

    # header + dispatch message
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Parallel injury measurement of {.val {length(files)}} images"),
        right = cli::col_blue("Started at {.val {format(init_time, '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_step(
        msg        = "Dispatching batches...",
        msg_done   = "All batches complete!",
        msg_failed = "Batch run failed."
      )
    }

    # run all images in parallel
    raw <- mirai::mirai_map(
      .x = files,
      .f = process_image
    )[.progress]

    # final completion rule
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_green("Parallel processing finished"),
        right = cli::col_blue("Ended at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
    }


  } else {
    # sequential
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Sequential injury measurement of {.val {length(files)}} images"),
        right = cli::col_blue("Started at {.val {format(init_time, '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | Current: {.val {cli::pb_status}}",
        total  = length(files),
        clear  = FALSE
      )
    }

    raw <- vector("numeric", length(files))
    for (i in seq_along(files)) {
      if (verbose) cli::cli_progress_update(status = basename(files[i]))
      raw[i] <- process_image(files[i])
    }

    if (verbose) {
      cli::cli_progress_done()
      cli::cli_rule(
        left  = cli::col_green("Sequential processing finished"),
        right = cli::col_blue("Ended at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
    }
  }

  # assemble and return
  names(raw) <- basename(files)
  df <- data.frame(
    img    = names(raw),
    injury = raw,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  df
}


