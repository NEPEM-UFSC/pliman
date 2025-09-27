#' Apply a function to images
#'
#' @description
#' Most of the functions in pliman can be applied to a list of images, but this can be not
#' ideal to deal with lots of images, mainly if they have a high resolution. For
#' curiosity, a 6000 x 4000 image use nearly 570 Megabytes of RAM. So, it would be
#' impossible to deal with lots of images within R. `apply_fun_to_img()` applies a
#' function to images stored in a given directory as follows:
#'  * Create a vector of image names that contain a given pattern of name.
#'  * Import each image of such a list.
#'  * Apply a function to the imported image.
#'  * Export the mutated image to the computer.
#'
#' If parallel is set to `FALSE` (default), the images are processed sequentially, which
#' means that one image needs to be imported, processed, and exported so that the other
#' image can be processed. If parallel is set to `TRUE`,  the images are processed
#' asynchronously (in parallel) in separate R sessions (3) running in the background on
#' the same machine. It may speed up the processing time when lots of images need to be
#' processed.
#'
#' @param pattern A pattern to match the images' names.
#' @param fun A function to apply to the images.
#' @param ... Arguments passed on to `fun`.
#' @param dir_original,dir_processed The directory containing the original and processed
#'   images. Defaults to `NULL`, which means that the current working directory will be
#'   considered. **The processed image will overwrite the original image unless a
#'   prefix/suffix be used or a subfolder is informed in dir_processed argument**.
#' @param prefix,suffix A prefix and/or suffix to be included in the name of processed
#'   images. Defaults to `""`.
#' @param parallel If `TRUE` processes the images asynchronously (in parallel) in separate
#'   R sessions (3 by default) running in the background on the same machine. It may speed
#'   up the processing time, especially when pattern is used is informed.
#' @param workers A positive numeric scalar or a function specifying the number of
#'   parallel processes that can be active at the same time. Defaults to 3.
#' @param verbose Shows the progress in console? Defaults to `TRUE`.
#'
#' @return Nothing. The processed images are saved to the current working directory.
#' @export
#'
#' @examples
#' # apply_fun_to_imgs("pattern", image_resize, rel_size = 50)
apply_fun_to_imgs <- function(pattern,
                               fun,
                               ...,
                               dir_original = NULL,
                               dir_processed = NULL,
                               prefix = "",
                               suffix = "",
                               parallel = FALSE,
                               workers = 3,
                               verbose = TRUE){
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

  imgs <- list.files(pattern = pattern, diretorio_original)
  # check for no images found
  if (length(imgs) == 0) {
    cli::cli_abort(c(
      "x" = "Pattern {.val {pattern}} not found in directory {.path {diretorio_original}}.",
      "i" = "Current working directory: {.path {getwd()}}"
    ))
  }

  # check for invalid extensions
  extensions <- as.character(vapply(imgs, tools::file_ext, character(1)))
  invalid <- setdiff(tolower(extensions), c("png", "jpeg", "jpg", "tiff"))
  if (length(invalid) > 0) {
    cli::cli_abort(c(
      "x" = "Found file(s) with unsupported extension{?s}: {.val {unique(invalid)}}.",
      "i" = "Allowed extensions are {.val png}, {.val jpeg}, {.val jpg}, {.val tiff}."
    ))
  }

  help_apply <- function(img,
                         fun,
                         ...,
                         diretorio_original,
                         diretorio_processada,
                         prefix = prefix,
                         suffix = suffix){
    img_file <- image_import(img, path = diretorio_original)
    img_list <- list(img_file)
    names(img_list) <- paste0(prefix, file_name(img), suffix, ".", file_extension(img))
    res <- lapply(img_list, fun, ...)
    image_export(res, subfolder = diretorio_processada)
  }

  if (parallel) {
    # decide number of workers
    workers <- ifelse(is.null(workers),
                      ceiling(parallel::detectCores() * 0.75),
                      workers)

    # silenciar mensagens do pliman, capturar opções antigas
    opar <- options(pliman_quiet = TRUE)
    on.exit(options(opar), add = TRUE)

    # iniciar workers persistentes
    mirai::daemons(workers)
    on.exit(mirai::daemons(0), add = TRUE)

    # carregar pliman nos workers (se necessário)
    mirai::everywhere({
      library(pliman)
    })

    # mensagens CLI
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Parallel processing using {workers} cores"),
        right = cli::col_blue("Started on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )

      cli::cli_progress_step(
        msg        = "Processing {.val {length(imgs)}} images in parallel...",
        msg_done   = "All batches complete!",
        msg_failed = "Batch failed"
      )
    }
    # função que será executada em paralelo
    process_image <- function(img) {
      help_apply(
        img                  = img,
        fun                  = fun,
        ...,
        diretorio_original   = diretorio_original,
        diretorio_processada = diretorio_processada,
        prefix               = prefix,
        suffix               = suffix
      )
    }

    # executa em paralelo com barra de progresso automática
    results <- mirai::mirai_map(
      .x = imgs,
      .f = process_image
    )[.progress]

  } else {
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Sequential processing of {.val {length(imgs)}} images"),
        right = cli::col_blue("Started on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [ETA: {cli::pb_eta}] | Current: {.val {cli::pb_status}}",
        total  = length(imgs),
        clear  = TRUE
      )
    }

    results <- vector("list", length(imgs))
    for (i in seq_along(imgs)) {
      if (verbose) {
        cli::cli_progress_update(status = imgs[i])
      }
      results[[i]] <- help_apply(
        img                  = imgs[[i]],
        fun                  = fun,
        ...,
        diretorio_original   = diretorio_original,
        diretorio_processada = diretorio_processada,
        prefix               = prefix,
        suffix               = suffix
      )
    }

    if (verbose) {
      cli::cli_progress_done()
      cli::cli_rule(
        left  = cli::col_blue("Function {.fn {deparse(substitute(fun))}} successfully applied to the images"),
        right = cli::col_blue("Finished on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
    }
  }

}




