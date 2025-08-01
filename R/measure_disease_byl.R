
#' Performs plant disease measurements by leaf
#'
#' Computes the percentage of symptomatic leaf area using color palettes or RGB
#' indexes `by` each `l`eaf of an image. This allows, for example, processing
#' replicates of the same treatment  and obtaining the results for each
#' replication with a single image. To do that, leaf samples are first splitten
#' with [object_split()] and then, [measure_disease()] is applied to the list of
#' leaves.
#'
#' @inheritParams object_split
#' @inheritParams measure_disease
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory. It can be either a full path, e.g.,
#'   `"C:/Desktop/imgs"`, or a subfolder within the current working directory,
#'   e.g., `"/imgs"`.
#' @param index A character value specifying the target mode for
#'   conversion to binary to segment the leaves from background. Defaults to "B"
#'   (blue). See [image_index()] for more details. Personalized indexes can be
#'   informed as, e.g., `index = "R*G/B`.
#' @param lower_size To prevent dust from affecting object segmentation, objects
#'   with lesser than `10%` of the mean of all objects are removed. . One can
#'   set a known area or use `lower_limit = 0` to select all objects (not
#'   advised).
#' @param ... Additional arguments passed on to [measure_disease()].
#' @return
#' * A list with the following objects:
#'   - `severity` A data frame with the percentage of healthy and symptomatic
#'  areas for each leaf in the image(s).
#'   - `shape`,`statistics` If `show_features = TRUE` is used, returns the shape
#'  (area, perimeter, etc.) for each lesion and a summary statistic of the
#'  results.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#'library(pliman)
#'img <- image_pliman("mult_leaves.jpg", plot = TRUE)
#'sev <-
#'  measure_disease_byl(img = img,
#'                      index_lb = "B",
#'                      index_dh = "NGRDI",
#'                      workers = 2)
#' sev$severity
#' }
#'
#'
measure_disease_byl <- function(img,
                                index = "B",
                                index_lb = "B",
                                index_dh = "NGRDI",
                                lower_size = NULL,
                                watershed = TRUE,
                                invert = FALSE,
                                fill_hull = FALSE,
                                opening = c(10, 0),
                                closing = c(0, 0),
                                filter = c(0, 0),
                                erode = c(0, 0),
                                dilate = c(0, 0),
                                threshold = "Otsu",
                                extension = NULL,
                                tolerance = NULL,
                                object_size = "large",
                                img_healthy = NULL,
                                img_symptoms = NULL,
                                plot = TRUE,
                                save_image = FALSE,
                                dir_original = NULL,
                                dir_processed = NULL,
                                pattern = NULL,
                                parallel = FALSE,
                                workers = NULL,
                                show_features = FALSE,
                                verbose = TRUE,
                                ...){
  check_ebi()

  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste("./", sep = "")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  if(is.character(img_healthy)){
    all_files <- sapply(list.files(diretorio_original), file_name)
    imag <- list.files(diretorio_original, pattern = img_healthy)
    check_names_dir(img_healthy, all_files, "")
    name_h <- file_name(imag)
    extens <- file_extension(imag)
    img_healthy <- image_import(paste(diretorio_original, "/", name_h, ".", extens, sep = ""))
  }
  if(is.character(img_symptoms)){
    all_files <- sapply(list.files(diretorio_original), file_name)
    imag <- list.files(diretorio_original, pattern = img_symptoms)
    check_names_dir(img_symptoms, all_files, "")
    name_h <- file_name(imag)
    extens <- file_extension(imag)
    img_symptoms <- image_import(paste(diretorio_original, "/", name_h, ".", extens, sep = ""))
  }
  back <- EBImage::Image(rep(1, 100*300),dim=c(100,300,3), colormode = 'Color')

  help_byl <- function(img,
                       img_healthy,
                       img_symptoms,
                       back,
                       index_dh,
                       index_lb){
    if(is.character(img)){
      all_files <- sapply(list.files(diretorio_original), file_name)
      check_names_dir(img, all_files, diretorio_original)
      imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
      name_ori <- file_name(imag)
      extens_ori <- file_extension(imag)
      img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
    } else{
      name_ori <- match.call()[[2]]
      extens_ori <- "jpg"
    }

    splits <- object_split(img,
                           index = index,
                           watershed = watershed,
                           invert = invert[[1]],
                           fill_hull = fill_hull,
                           opening = opening[[1]],
                           closing = closing[[1]],
                           filter = filter[[1]],
                           erode = erode[[1]],
                           dilate = dilate[[1]],
                           threshold = threshold,
                           extension = extension,
                           tolerance = tolerance,
                           object_size = object_size,
                           remove_bg = TRUE,
                           plot = FALSE,
                           verbose = FALSE)
    results <- list()
    if(is.null(img_healthy)){
      results <-
        lapply(seq_along(splits),
               function(i){
                 measure_disease(splits[[i]],
                                 index_dh = index_dh,
                                 index_lb = index_lb,
                                 plot = plot,
                                 save_image = save_image,
                                 show_features = show_features,
                                 dir_processed = diretorio_processada,
                                 prefix = paste0(name_ori, "_", i),
                                 name = "",
                                 opening = opening,
                                 closing = closing,
                                 filter = filter,
                                 erode = erode,
                                 dilate = dilate,
                                 threshold = threshold,
                                 invert = ifelse(length(invert) > 1, invert[[2]], invert[[1]]),
                                 ...)
               })

    } else{
      results <-
        lapply(seq_along(splits),
               function(i){
                 measure_disease(splits[[i]],
                                 img_healthy = img_healthy,
                                 img_symptoms = img_symptoms,
                                 img_background = back,
                                 plot = plot,
                                 show_features = show_features,
                                 save_image = save_image,
                                 dir_processed = diretorio_processada,
                                 prefix = paste0(name_ori, "_", i),
                                 name = "",
                                 opening = opening,
                                 closing = closing,
                                 filter = filter,
                                 erode = erode,
                                 dilate = dilate,
                                 threshold = threshold,
                                 ...)
               })
    }
    names(results) <- paste0(name_ori, "-", 1:length(results))
    severity <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["severity"]],
                          img =  names(results[i]))[, c(3, 1:2)]
              })
      ) |>
      separate_col(img, into = c("img", "leaf"), sep = "-")
    if(isTRUE(show_features)){
      stats <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  if(is.null(results[[i]][["statistics"]])){
                    res <- data.frame(stat = "NA", value = 0)
                  } else{
                    res <- results[[i]][["statistics"]]
                  }
                  transform(res, img =  names(results[i]))[,c(3, 1, 2)]
                })
        ) |>
        separate_col(img, into = c("img", "leaf"), sep = "-")
      shape <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  if(is.null(results[[i]][["shape"]])){
                    names <- names_measures()
                    res <- data.frame(matrix(nrow = 1, ncol = length(names)))
                    res[1, ] <- 0
                    colnames(res) <- names
                  } else{
                    res <- results[[i]][["shape"]]
                  }
                  transform(res, img =  names(results[i]))
                })
        ) |>
        separate_col(img, into = c("img", "leaf"), sep = "-")

    } else{
      shape <- NULL
      stats <- NULL
    }
    invisible(
      structure(
        list(severity = severity,
             stats = stats,
             shape = shape),
        class = "plm_disease_byl")
    )
  }

  if(missing(pattern)){
    if(verbose){
      cli::cli_progress_step(
        msg = "Processing a single image. Please, wait.",
        msg_done = "Image {.emph Successfully} analyzed!",
        msg_failed = "Oops, something went wrong."
      )
    }
    results <- help_byl(img, img_healthy, img_symptoms, back, index_dh, index_lb)
  } else {
    if (pattern %in% as.character(0:9)) {
      old_pattern <- pattern
      pattern <- "^[0-9].*$"
      cli::cli_alert_info(
        "Numeric pattern {.val {old_pattern}} converted to {.val {pattern}}."
      )
    }

    # list and name files
    plants       <- list.files(pattern = pattern, diretorio_original)
    extensions   <- tolower(vapply(plants, tools::file_ext,   ""))
    names_plant  <-         vapply(plants, tools::file_path_sans_ext, "")


    if (length(plants) == 0) {
      cli::cli_abort(c(
        "x" = "Pattern {.val {pattern}} not found in directory {.path {diretorio_original}}.",
        "i" = "Check your working directory: {.path {getwd()}}"
      ))
    }


    bad <- setdiff(extensions, c("png", "jpeg", "jpg", "tiff"))
    if (length(bad)) {
      cli::cli_abort(c(
        "x" = "Unsupported file extension{?s}: {.val {unique(bad)}} found.",
        "i" = "Allowed: {.val png}, {.val jpeg}, {.val jpg}, {.val tiff}."
      ))
    }


    if (parallel) {
      # define número de workers
      nworkers <- ifelse(is.null(workers),
                         trunc(parallel::detectCores() * 0.5),
                         workers)

      # inicia os workers do mirai
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)

      # mensagens CLI
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Parallel processing of {.val {length(names_plant)}} images"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Dispatching batches...",
          msg_done   = "All batches complete!",
          msg_failed = "Batch run failed."
        )
      }

      # define função para aplicar em paralelo
      process_image <- function(img) {
        help_byl(
          img          = img,
          img_healthy  = img_healthy,
          img_symptoms = img_symptoms,
          back         = back,
          index_dh     = index_dh,
          index_lb     = index_lb
        )
      }

      # executa com barra de progresso automática
      results <- mirai::mirai_map(
        .x = names_plant,
        .f = process_image
      )[.progress]


    } else {
      # sequential header + bar
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Sequential processing of {.val {length(names_plant)}} images"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_bar(
          format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | Current: {.val {cli::pb_status}}",
          total  = length(names_plant),
          clear  = FALSE
        )
      }

      # run loop
      results <- vector("list", length(names_plant))
      for (i in seq_along(names_plant)) {
        if (verbose) cli::cli_progress_update(status = names_plant[i])
        results[[i]] <- help_byl(
          img          = names_plant[i],
          img_healthy, img_symptoms, back,
          index_dh, index_lb
        )
      }

      # done rule
      if (verbose) {
        cli::cli_progress_done()
        cli::cli_rule(
          left  = cli::col_green("Finished processing {.val {length(names_plant)}} images"),
          right = cli::col_blue("Ended at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
      }
    }
  }
  invisible(structure(
    results, class = "plm_disease_byl"
  ))
}




