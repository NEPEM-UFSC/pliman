# ==============================================================================
# DEEP LEARNING BENCHMARK AND SESSION DIAGNOSTICS
# pliman: Plant Image Analysis in R
# ==============================================================================

#' Benchmark YOLO inference speed
#'
#' @description
#' `yolo_benchmark()` measures how fast a YOLO ONNX model runs on the current machine,
#' separating the raw C++ inference path (`run_yolo_cpp()`: letterbox + ONNX Runtime +
#' decoding/NMS) from the complete R pipeline of [image_detect_dl()] (model lookup,
#' post-processing and result tables). Optionally it also times tiled inference with
#' [image_detect_sahi()].
#'
#' Warm-up iterations are run first and excluded from the timings, so the first-call
#' costs (session creation, graph optimisation and, on GPU, DirectML shader compilation)
#' do not distort the averages.
#'
#' @details
#' Speed depends heavily on the hardware and on the model input size (`imgsz`). Some
#' practical notes:
#' \itemize{
#'   \item `engine = "gpu"` uses DirectML (Windows / DirectX 12). Gains depend on the GPU;
#'     integrated GPUs may be slower than a multi-core CPU for small models.
#'   \item FP16 models ([yolo_export_fp16()]) are faster only on GPUs with native
#'     half-precision throughput (e.g. NVIDIA Tensor Cores, RTX 20xx and newer). On GPUs
#'     without fast FP16 (e.g. GTX 10xx/16xx) or on CPU, the casting overhead can make
#'     them slower than FP32.
#'   \item Dynamic INT8 models ([onnx_quantize()]) target CPU execution; on DirectML they
#'     may fall back to slower kernels.
#'   \item With `engine = "cpu"`, `threads` controls the ONNX Runtime intra-op threads
#'     (`0` = all physical cores).
#' }
#'
#' @param model Character. Model name (e.g. `"yolo26n"`) or path to a `.onnx` file.
#' @param img Image used for timing: an `image` object, a file path or `NULL` (default),
#'   which uses the example image `"soybean_touch.jpg"`.
#' @param imgsz Optional integer inference size. `NULL` uses the model's native size.
#' @param engine Execution engine: `"gpu"` (DirectML, Windows) or `"cpu"`.
#' @param threads Integer. ONNX Runtime CPU threads (`0` = runtime default).
#' @param n Integer. Number of timed iterations per stage.
#' @param warmup Integer. Untimed warm-up iterations before each stage.
#' @param device_id Integer GPU index (`-1` = adapter with the most VRAM).
#' @param sahi Logical. Also time [image_detect_sahi()] with 640 px slices.
#' @param conf_threshold Confidence threshold used during timing.
#' @param verbose Logical. Print a summary.
#'
#' @return A data.frame (class `pliman_benchmark`) with one row per stage and columns
#'   `stage`, `ms_mean`, `ms_sd`, `ms_min`, `fps` and `n_objects`. Attributes hold the
#'   model file, engine, threads, image size and the active ONNX sessions.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' yolo_benchmark("yolo26n", engine = "cpu")
#' yolo_benchmark("yolo26n", engine = "gpu", sahi = TRUE)
yolo_benchmark <- function(model = "yolo26n",
                           img = NULL,
                           imgsz = NULL,
                           engine = c("gpu", "cpu"),
                           threads = 0,
                           n = 20,
                           warmup = 3,
                           device_id = -1,
                           sahi = FALSE,
                           conf_threshold = 0.25,
                           verbose = TRUE) {
  engine <- match.arg(engine)
  n <- max(1L, as.integer(n))
  warmup <- max(0L, as.integer(warmup))
  im <- if (is.null(img)) {
    image_pliman("soybean_touch.jpg")
  } else if (is.character(img)) {
    image_import(img[1])
  } else {
    img
  }
  mat <- if (isS4(im) && methods::.hasSlot(im, ".Data")) methods::slot(im, ".Data") else unclass(im)
  if (!is.raw(mat) && !is.double(mat)) storage.mode(mat) <- "double"
  rt <- .resolve_yolo_runtime(model, imgsz = imgsz, threads = threads, engine = engine, device_id = device_id)
  if (is.null(rt)) {
    cli::cli_abort("ONNX Runtime library not found. Run {.code onnx_install()} first.")
  }
  w <- dim(mat)[1]
  h <- dim(mat)[2]

  time_stage <- function(fun) {
    for (i in seq_len(warmup)) fun()
    t <- numeric(n)
    n_obj <- NA_integer_
    for (i in seq_len(n)) {
      t0 <- proc.time()[["elapsed"]]
      r <- fun()
      t[i] <- (proc.time()[["elapsed"]] - t0) * 1000
      if (i == 1L) n_obj <- tryCatch({
        if (is.data.frame(r)) nrow(r) else if (!is.null(r$boxes)) NROW(r$boxes) else NA_integer_
      }, error = function(e) NA_integer_)
    }
    c(ms_mean = mean(t), ms_sd = if (n > 1) stats::sd(t) else 0, ms_min = min(t), n_objects = n_obj)
  }

  stages <- list()
  stages[["inference (C++)"]] <- time_stage(function() {
    run_yolo_cpp(mat, w, h, conf_threshold, 0.45, rt$model_file, rt$lib_path, rt$threads,
                 rt$use_gpu, rt$device_id, need_masks = FALSE, target_size = rt$target_size)
  })
  stages[["image_detect_dl()"]] <- time_stage(function() {
    image_detect_dl(im, model = rt$model_file, conf_threshold = conf_threshold, imgsz = rt$target_size,
                    threads = threads, engine = engine, device_id = device_id,
                    plot = FALSE, verbose = FALSE)
  })
  if (isTRUE(sahi)) {
    stages[["image_detect_sahi()"]] <- time_stage(function() {
      image_detect_sahi(im, model = rt$model_file, conf_threshold = conf_threshold,
                        slice_height = 640, slice_width = 640, imgsz = rt$target_size,
                        threads = threads, engine = engine, device_id = device_id,
                        plot = FALSE, verbose = FALSE)
    })
  }

  out <- data.frame(
    stage = names(stages),
    ms_mean = vapply(stages, function(s) s[["ms_mean"]], numeric(1)),
    ms_sd = vapply(stages, function(s) s[["ms_sd"]], numeric(1)),
    ms_min = vapply(stages, function(s) s[["ms_min"]], numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out$fps <- ifelse(out$ms_mean > 0, 1000 / out$ms_mean, NA_real_)
  out$n_objects <- vapply(stages, function(s) as.integer(s[["n_objects"]]), integer(1))
  sess <- tryCatch(onnx_session_info_cpp(), error = function(e) NULL)
  attr(out, "model") <- rt$model_file
  attr(out, "engine") <- engine
  attr(out, "gpu_active") <- if (!is.null(sess)) any(sess$gpu[sess$model == rt$model_file]) else NA
  attr(out, "threads") <- threads
  attr(out, "imgsz") <- rt$target_size
  attr(out, "image_size") <- c(width = w, height = h)
  class(out) <- c("pliman_benchmark", "data.frame")
  if (isTRUE(verbose)) print(out)
  invisible(out)
}

#' @export
print.pliman_benchmark <- function(x, digits = 2, ...) {
  gpu_txt <- if (isTRUE(attr(x, "gpu_active"))) "GPU (DirectML)" else if (identical(attr(x, "engine"), "gpu")) "CPU (GPU unavailable)" else "CPU"
  isz <- attr(x, "image_size")
  cli::cli_rule(left = "YOLO benchmark: {.file {basename(attr(x, 'model'))}}")
  cli::cli_alert_info("Device: {gpu_txt} | imgsz: {attr(x, 'imgsz')} | image: {isz[1]}x{isz[2]} | threads: {attr(x, 'threads')}")
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df[c("ms_mean", "ms_sd", "ms_min", "fps")] <- lapply(df[c("ms_mean", "ms_sd", "ms_min", "fps")], round, digits)
  print(df, row.names = FALSE)
  invisible(x)
}

#' Inspect and configure cached ONNX Runtime sessions
#'
#' @description
#' Deep-learning functions in pliman keep ONNX Runtime sessions in memory, so the model
#' is loaded (and, on GPU, compiled) only once per R session. `pliman_session_info()`
#' lists the cached sessions; `pliman_max_sessions()` reads or changes how many sessions
#' are kept (least recently used sessions are released first). Use
#' [pliman_clear_sessions()] to release all of them.
#'
#' @details
#' Pipelines such as Grounded-SAM (Grounding DINO + SAM encoder + SAM decoder) or
#' detection with `return_features = TRUE` (YOLO + CLIP/DINOv2) need several sessions at
#' once. The default limit is 8 and can also be set before loading the package with the
#' environment variable `PLIMAN_ONNX_MAX_SESSIONS` or with
#' `options(pliman.max_sessions = n)`. Lower it if RAM or GPU memory is tight.
#'
#' ONNX Runtime worker threads keep the runtime default of busy-waiting ("spinning")
#' briefly between runs, which gives the lowest latency. When several R processes share
#' the CPU (e.g. parallel batch workers), set the environment variable
#' `PLIMAN_ORT_SPIN=0` before the first inference to disable spinning.
#'
#' @param n Integer. New maximum number of cached sessions (`>= 1`). `NULL` only
#'   returns the current value.
#' @return `pliman_session_info()`: a data.frame with one row per cached session
#'   (model, whether it runs on the GPU, threads, input dimensions, warm-up and CPU
#'   fallback flags), most recently used first. `pliman_max_sessions()`: the previous
#'   limit, invisibly when `n` is given.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' res <- image_detect_dl(img, model = "yolo26n", plot = FALSE)
#' pliman_session_info()
#' pliman_max_sessions()
pliman_session_info <- function() {
  info <- tryCatch(onnx_session_info_cpp(), error = function(e) NULL)
  if (is.null(info)) {
    return(data.frame(model = character(0), gpu = logical(0), stringsAsFactors = FALSE))
  }
  info$model <- basename(info$model)
  info
}

#' @rdname pliman_session_info
#' @export
pliman_max_sessions <- function(n = NULL) {
  if (is.null(n)) return(onnx_set_max_sessions_cpp(0L))
  n <- as.integer(n)
  if (is.na(n) || n < 1L) cli::cli_abort("{.arg n} must be a positive integer.")
  invisible(onnx_set_max_sessions_cpp(n))
}
