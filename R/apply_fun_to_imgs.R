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
#' @param as_job Logical. If `TRUE` and running inside RStudio, runs the image processing
#'   as an RStudio Background Job (default: `FALSE`).
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
                               verbose = TRUE,
                               as_job = FALSE){
  # 0. Support for RStudio Background Jobs (as_job = TRUE)
  if (isTRUE(as_job)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
      cli::cli_alert_warning("RStudio Job API is not available (not running inside interactive RStudio). Proceeding with image processing in current session.")
    } else {
      mc <- match.call(expand.dots = TRUE)
      eval_env <- parent.frame()
      evaluated_args <- list()
      arg_names <- setdiff(names(mc), c("", "as_job"))
      for (arg_name in arg_names) {
        arg_val <- tryCatch(eval(mc[[arg_name]], envir = eval_env), error = function(e) mc[[arg_name]])
        evaluated_args[[arg_name]] <- arg_val
      }
      evaluated_args$as_job <- FALSE

      job_dir <- file.path(tempdir(), "pliman_jobs")
      if (!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

      time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      job_script_path <- file.path(job_dir, paste0("apply_fun_job_", time_stamp, ".R"))
      job_args_rds <- file.path(job_dir, paste0("apply_fun_args_", time_stamp, ".rds"))

      saveRDS(list(args = evaluated_args, getwd = getwd()), job_args_rds)

      fun_name <- tryCatch({
        deparse(substitute(fun))
      }, error = function(e) "function")
      if (length(fun_name) > 1L || !nzchar(fun_name)) fun_name <- "function"

      pat_name <- if (is.character(pattern) && nzchar(pattern[1])) pattern[1] else "images"
      j_name <- paste0("Apply ", fun_name, ": ", pat_name)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::apply_fun_to_imgs()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman apply_fun_to_imgs() - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "do.call(pliman::apply_fun_to_imgs, job_info$args)\n",
        "cat('\\n[OK] Image processing background job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("Apply Function to Images (pliman)")
      cli::cli_alert_success("Image processing task has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Working directory: {.path {getwd()}}",
        "Pattern: {.val {pat_name}}",
        "Function: {.val {fun_name}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track live progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    }
  }

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

    if(verbose){
      cli::cli_progress_step(
        msg        = "Processing {.val {length(imgs)}} images in parallel...",
        msg_done   = "All batches complete!",
        msg_failed = "Batch failed"
      )
    }

  } else {
    if (verbose) {
      cli::cli_rule(
        left  = cli::col_blue("Sequential processing of {.val {length(imgs)}} images"),
        right = cli::col_blue("Started on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [ETA: {cli::pb_eta}]",
        total  = length(imgs),
        clear  = TRUE
      )
    }

    results <- vector("list", length(imgs))
    for (i in seq_along(imgs)) {
      if (verbose) {
        cli::cli_progress_update()
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




