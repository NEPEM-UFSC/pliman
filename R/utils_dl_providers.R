# ==============================================================================
# ONNX RUNTIME EXECUTION PROVIDERS, PROFILING AND SESSION START-UP CACHE
# pliman: Plant Image Analysis in R
#
# Process-wide settings read by the C++ session factory (src/onnx_inference.cpp,
# get_or_create_cached_session()). They are passed through environment variables so
# that parallel workers (mirai, callr) started afterwards inherit them.
# ==============================================================================

.gpu_provider_choices <- c("auto", "dml", "cuda", "tensorrt", "coreml", "openvino")

#' Choose the GPU execution provider used by ONNX Runtime
#'
#' @description
#' Every deep-learning function in pliman that takes `engine = "gpu"` runs its ONNX
#' models through an ONNX Runtime *execution provider*. `pliman_gpu_provider()` reads or
#' sets which provider is tried first:
#'
#' * `"auto"` (default): CUDA when the loaded ONNX Runtime library includes it,
#'   otherwise DirectML (Windows). This reproduces the previous behaviour with the
#'   libraries installed by [onnx_install()].
#' * `"dml"`: DirectML (Windows 10/11, any DirectX 12 GPU: NVIDIA, AMD, Intel).
#' * `"cuda"`: NVIDIA CUDA. Needs the CUDA build of ONNX Runtime
#'   (`onnx_install(engine = "cuda")`) **plus** the NVIDIA CUDA 12 runtime, cuBLAS and
#'   cuDNN 9 libraries on the library search path (see `cuda_path`).
#' * `"tensorrt"`: NVIDIA TensorRT on top of CUDA (also needs the TensorRT libraries).
#'   Engines are compiled the first time a model is loaded, which can take minutes;
#'   they are cached in `trt_cache_dir` so later sessions start quickly.
#' * `"coreml"`: Apple CoreML (macOS builds of ONNX Runtime).
#' * `"openvino"`: Intel OpenVINO (requires an OpenVINO build of ONNX Runtime). The
#'   device is taken from the environment variable `PLIMAN_OPENVINO_DEVICE`
#'   (default `"AUTO:GPU,CPU"`).
#'
#' A provider that the loaded library does not contain, or that fails to initialise
#' (missing driver or CUDA libraries), is skipped: pliman then tries CUDA, DirectML and
#' finally the CPU, so inference never fails just because of this setting. Use
#' [pliman_session_info()] (column `provider`) to see which provider actually serves
#' each model, and [onnx_available_providers()] to list what the library offers.
#'
#' @section Hardware dependence:
#' Speed-ups depend entirely on the GPU, its driver and, for CUDA/TensorRT, the
#' installed NVIDIA libraries. TensorRT with `trt_fp16 = TRUE` and FP16 models only
#' accelerate on GPUs with native half-precision throughput (NVIDIA Tensor Cores, RTX
#' 20xx and newer); on GPUs without it (e.g. GTX 10xx/16xx) they can be slower than FP32
#' because of casting. Always compare with [yolo_benchmark()] on your own machine.
#'
#' @param provider One of `"auto"`, `"dml"`, `"cuda"`, `"tensorrt"`, `"coreml"`,
#'   `"openvino"`. `NULL` (default) only returns the current setting.
#' @param trt_cache_dir Directory for TensorRT engine and timing caches. Defaults to
#'   `file.path(pliman_model_dir(), "trt_cache")` when `provider = "tensorrt"`.
#' @param trt_fp16 Logical. Let TensorRT build FP16 kernels. `NULL` keeps the current
#'   setting (initially `FALSE`).
#' @param cuda_path Optional character vector of directories holding the CUDA runtime,
#'   cuBLAS and cuDNN libraries (e.g. `"C:/Program Files/NVIDIA GPU Computing
#'   Toolkit/CUDA/v12.4/bin"`). They are prepended to `PATH` (Windows) or
#'   `LD_LIBRARY_PATH` for the current R process, before the first GPU session is
#'   created.
#' @return The previous provider, invisibly when a new one is set.
#' @seealso [onnx_available_providers()], [pliman_session_info()], [onnx_install()]
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' pliman_gpu_provider()             # "auto"
#' onnx_available_providers()        # what the installed runtime offers
#'
#' # NVIDIA CUDA (after onnx_install(engine = "cuda"))
#' pliman_gpu_provider("cuda", cuda_path = "C:/cuda/bin")
#' img <- image_pliman("soybean_touch.jpg")
#' res <- image_detect_dl(img, model = "yolo26n", plot = FALSE)
#' pliman_session_info()[, c("model", "provider")]
#'
#' pliman_gpu_provider("auto")       # back to the default
pliman_gpu_provider <- function(provider = NULL,
                                trt_cache_dir = NULL,
                                trt_fp16 = NULL,
                                cuda_path = NULL) {
  prev <- tolower(Sys.getenv("PLIMAN_GPU_PROVIDER", "auto"))
  if (!prev %in% .gpu_provider_choices) prev <- "auto"
  if (!is.null(cuda_path)) .prepend_library_path(cuda_path)
  if (!is.null(trt_fp16)) Sys.setenv(PLIMAN_TRT_FP16 = if (isTRUE(trt_fp16)) "1" else "0")
  if (!is.null(trt_cache_dir)) {
    dir.create(trt_cache_dir, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(PLIMAN_TRT_CACHE_DIR = normalizePath(trt_cache_dir, winslash = "/", mustWork = FALSE))
  }
  if (is.null(provider)) {
    if (is.null(trt_cache_dir) && is.null(trt_fp16) && is.null(cuda_path)) return(prev)
    return(invisible(prev))
  }
  provider <- tolower(as.character(provider[1]))
  if (provider == "directml") provider <- "dml"
  if (provider == "trt") provider <- "tensorrt"
  provider <- match.arg(provider, .gpu_provider_choices)
  if (provider == "tensorrt" && !nzchar(Sys.getenv("PLIMAN_TRT_CACHE_DIR"))) {
    d <- file.path(pliman_model_dir(), "trt_cache")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(PLIMAN_TRT_CACHE_DIR = normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  Sys.setenv(PLIMAN_GPU_PROVIDER = provider)
  # the library path depends on the provider (CUDA builds live in their own folder)
  .clear_onnx_lib_cache()
  invisible(prev)
}

# Drops the memoised ONNX Runtime library paths (see pliman_onnx_lib_path()).
.clear_onnx_lib_cache <- function() {
  keys <- ls(envir = .onnx_meta_cache, all.names = TRUE)
  keys <- keys[startsWith(keys, "onnx_lib_path::")]
  if (length(keys)) rm(list = keys, envir = .onnx_meta_cache)
  invisible(NULL)
}

.prepend_library_path <- function(dirs) {
  dirs <- normalizePath(dirs[dir.exists(dirs)], winslash = "\\", mustWork = FALSE)
  if (!length(dirs)) {
    cli::cli_alert_warning("None of the {.arg cuda_path} directories exist; library search path unchanged.")
    return(invisible(FALSE))
  }
  var <- if (.Platform$OS.type == "windows") "PATH" else "LD_LIBRARY_PATH"
  cur <- Sys.getenv(var, "")
  parts <- strsplit(cur, .Platform$path.sep, fixed = TRUE)[[1]]
  new <- c(dirs[!dirs %in% parts], parts)
  do.call(Sys.setenv, stats::setNames(list(paste(new, collapse = .Platform$path.sep)), var))
  invisible(TRUE)
}

# Library override used by pliman_onnx_lib_path(): an explicit PLIMAN_ONNX_LIB, or the
# CUDA build installed by onnx_install(engine = "cuda") when the CUDA or TensorRT
# provider is selected. The override applies to both engines on purpose: alternating
# two different runtime libraries in one process forces every cached session to be
# rebuilt, and the CUDA build also contains the CPU provider.
.onnx_lib_override <- function() {
  ev <- Sys.getenv("PLIMAN_ONNX_LIB", "")
  if (nzchar(ev) && file.exists(ev)) return(normalizePath(ev, winslash = "/"))
  prov <- tolower(Sys.getenv("PLIMAN_GPU_PROVIDER", "auto"))
  if (prov %in% c("cuda", "tensorrt")) {
    p <- file.path(pliman_onnx_dir(), "cuda", .onnx_lib_name())
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NULL
}

.onnx_lib_name <- function() {
  switch(tolower(Sys.info()[["sysname"]]),
         "windows" = "onnxruntime.dll",
         "linux"   = "libonnxruntime.so",
         "darwin"  = "libonnxruntime.dylib",
         "onnxruntime.dll")
}

# Downloads the CUDA/TensorRT-enabled ONNX Runtime build from the official GitHub
# release into pliman_onnx_dir()/cuda. Called by pliman_download_onnx(engine = "cuda").
.download_onnx_cuda <- function(version, dir, force = FALSE) {
  os <- tolower(Sys.info()[["sysname"]])
  if (!os %in% c("windows", "linux")) {
    cli::cli_abort("The CUDA build of ONNX Runtime is only published for Windows and Linux (x64).")
  }
  dest_dir <- file.path(dir, "cuda")
  lib_name <- .onnx_lib_name()
  dest_lib <- file.path(dest_dir, lib_name)
  if (!isTRUE(force) && file.exists(dest_lib) && file.info(dest_lib)$size > 1000000) {
    return(normalizePath(dest_lib, winslash = "/"))
  }
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  url <- if (os == "windows") {
    paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version,
           "/onnxruntime-win-x64-gpu-", version, ".zip")
  } else {
    paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version,
           "/onnxruntime-linux-x64-gpu-", version, ".tgz")
  }
  archive <- file.path(tempdir(), basename(url))
  old_timeout <- getOption("timeout")
  on.exit({
    options(timeout = old_timeout)
    unlink(archive)
  }, add = TRUE)
  options(timeout = max(1800, old_timeout))
  cli::cli_alert_info("Downloading ONNX Runtime {version} with CUDA/TensorRT providers (~250 MB)...")
  if (requireNamespace("curl", quietly = TRUE)) {
    curl::curl_download(url, destfile = archive, quiet = TRUE, mode = "wb")
  } else {
    utils::download.file(url, destfile = archive, mode = "wb", quiet = TRUE)
  }
  if (!file.exists(archive) || file.info(archive)$size < 1e6) {
    cli::cli_abort("Failed to download {.url {url}}.")
  }
  ex <- file.path(tempdir(), paste0("ort_cuda_", Sys.getpid()))
  dir.create(ex, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(ex, recursive = TRUE), add = TRUE)
  if (os == "windows") utils::unzip(archive, exdir = ex) else utils::untar(archive, exdir = ex)
  pat <- if (os == "windows") "\\.dll$" else "\\.so(\\.[0-9.]+)?$"
  libs <- list.files(ex, pattern = pat, recursive = TRUE, full.names = TRUE)
  libs <- libs[grepl("onnxruntime", basename(libs))]
  if (!length(libs)) cli::cli_abort("No ONNX Runtime libraries found in the downloaded archive.")
  file.copy(libs, dest_dir, overwrite = TRUE)
  if (os == "linux" && !file.exists(dest_lib)) {
    versioned <- list.files(dest_dir, pattern = "^libonnxruntime\\.so\\.[0-9]", full.names = TRUE)
    if (length(versioned)) file.copy(versioned[1], dest_lib, overwrite = TRUE)
  }
  if (!file.exists(dest_lib)) cli::cli_abort("Could not install {.file {lib_name}}.")
  cli::cli_alert_success("ONNX Runtime (CUDA/TensorRT) installed at {.path {dest_dir}}.")
  cli::cli_bullets(c(
    "i" = "Activate it with {.code pliman_gpu_provider(\"cuda\")} (or {.val tensorrt}).",
    "!" = "CUDA 12 runtime, cuBLAS and cuDNN 9 must be installed and on the library search path; pass their folder with {.arg cuda_path} if needed."
  ))
  normalizePath(dest_lib, winslash = "/")
}

# Grounding DINO variant for one call.
#
# The distributed `groundingdino-tiny.onnx` is dynamically quantized (INT8:
# MatMulInteger / DynamicQuantizeLinear / ConvInteger). Profiled on a GTX 1650:
# * DirectML runs 99.8% of the node time on the GPU (DML fuses the INT8 pattern into
#   MatMulIntegerToFloat), ~0.8 s per call vs ~8.6 s on the 12-thread CPU;
# * CUDA has no kernels for those operators: 392 MatMulInteger and 291
#   DynamicQuantizeLinear nodes fall back to the CPU (68% of the time, plus 693
#   host-to-device copies), ~3.2 s per call.
# So the FP16 export (`groundingdino-tiny_fp16.onnx`, same graph, FP32 I/O) is used
# automatically only with the CUDA/TensorRT providers, when the file is present in the
# model directory. `options(pliman.dino_model = <file name or path>)` overrides the
# choice for every engine.
.resolve_dino_model <- function(dir, use_gpu) {
  opt <- getOption("pliman.dino_model")
  if (!is.null(opt) && nzchar(opt)) {
    p <- if (file.exists(opt)) opt else file.path(dir, opt)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
    cli::cli_alert_warning("{.code options(pliman.dino_model)} = {.file {opt}} not found; using the default Grounding DINO model.")
  }
  if (isTRUE(use_gpu) && tolower(Sys.getenv("PLIMAN_GPU_PROVIDER", "auto")) %in% c("cuda", "tensorrt")) {
    p <- file.path(dir, "groundingdino-tiny_fp16.onnx")
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  file.path(dir, "groundingdino-tiny.onnx")
}

# SAM 2.1 image encoder variant for one call. Controlled by
# options(pliman.sam_encoder = ...):
#   "fp32"  (default) the distributed sam2.1.encoder.onnx
#   "int8"  sam2.1.encoder_int8.onnx (dynamic INT8; CPU sessions only)
#   "fp16"  sam2.1.encoder_fp16.onnx (FP16 weights, FP32 I/O; GPU with Tensor Cores)
#   or a file name in `dir` / a path.
# Variants are built with pliman_build_model_variant(); a missing file falls back to
# FP32. There is deliberately no automatic choice. Measured on an i7 (12 threads) +
# GTX 1650, encode only, mask IoU against FP32 over three prompts:
#
#                  FP32      INT8               FP16
#   CPU           2.31 s    2.31 s (IoU .995)  5.75 s (IoU .9999)
#   DirectML      0.36 s   38.3 s  (IoU 0!)   36.7 s  (IoU .9999)
#   CUDA          0.39 s    4.65 s (IoU .996)  1.52 s (IoU .9998)
#
# INT8 only pays off on CPUs with VNNI/AMX and returns wrong masks under DirectML, so
# it is refused for GPU sessions; FP16 needs native half-precision throughput (Tensor
# Cores), without which it is 4-100x slower.
.resolve_sam_encoder <- function(dir, use_gpu) {
  base <- file.path(dir, "sam2.1.encoder.onnx")
  opt <- getOption("pliman.sam_encoder", "fp32")
  if (is.null(opt) || !nzchar(opt[1])) return(base)
  opt <- as.character(opt[1])
  variant <- switch(tolower(opt),
    "fp32" = return(base),
    "int8" = "sam2.1.encoder_int8.onnx",
    "fp16" = "sam2.1.encoder_fp16.onnx",
    opt
  )
  if (isTRUE(use_gpu) && grepl("int8", basename(variant), ignore.case = TRUE)) {
    .warn_once("sam_enc_int8_gpu", c(
      "!" = "The INT8 SAM 2.1 encoder is not used on the GPU (DirectML returned wrong masks with it); using FP32.",
      "i" = "INT8 is meant for {.code engine = \"cpu\"} / {.code sam_device = \"cpu\"}."
    ))
    return(base)
  }
  p <- if (file.exists(variant)) variant else file.path(dir, variant)
  if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  .warn_once(paste0("sam_enc_missing_", variant), c(
    "!" = "SAM 2.1 encoder variant {.file {variant}} not found; using FP32.",
    "i" = "Build it with {.code pliman_build_model_variant()}."
  ))
  base
}

.warned_env <- new.env(parent = emptyenv())
.warn_once <- function(key, msg, .envir = parent.frame()) {
  if (isTRUE(.warned_env[[key]])) return(invisible(FALSE))
  assign(key, TRUE, envir = .warned_env)
  cli::cli_bullets(msg, .envir = .envir)
  invisible(TRUE)
}

#' List the execution providers of the installed ONNX Runtime
#'
#' Returns the execution providers compiled into the ONNX Runtime library that pliman
#' loads for `engine` (for example `"DmlExecutionProvider"`, `"CUDAExecutionProvider"`,
#' `"CPUExecutionProvider"`). A provider being listed means the library supports it; it
#' can still fail to start if the matching drivers or libraries are missing, in which
#' case pliman falls back as described in [pliman_gpu_provider()].
#'
#' @param engine `"gpu"` (default) or `"cpu"`: which library to query.
#' @return A character vector of provider names.
#' @export
#' @examplesIf interactive()
#' onnx_available_providers()
onnx_available_providers <- function(engine = c("gpu", "cpu")) {
  engine <- match.arg(engine)
  lib <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib) || !file.exists(lib)) {
    cli::cli_alert_warning("ONNX Runtime library not found. Run {.code onnx_install()} first.")
    return(character(0))
  }
  onnx_available_providers_cpp(lib)
}

#' Profile ONNX Runtime sessions
#'
#' @description
#' `pliman_onnx_profiling()` turns ONNX Runtime's built-in profiler on or off for
#' sessions created afterwards. Each session writes one JSON trace, named after the
#' model and the provider, when it is released, i.e. on [pliman_clear_sessions()], LRU
#' eviction, or at the end of the R session. `onnx_profile_summary()` reads such a
#' trace and reports, per execution provider and operator type, how many nodes ran
#' there and how much time they took.
#'
#' The main use is to find **CPU fallbacks inside a GPU session**: operators that the
#' GPU provider cannot execute run on the CPU, with a device-to-host copy around each
#' one. A model with many such nodes (for example dynamically quantized INT8 models on
#' DirectML) may be faster entirely on the CPU, or need an FP16/FP32 variant.
#'
#' @param enable Logical. `TRUE` starts profiling new sessions, `FALSE` stops.
#' @param dir Directory for the JSON traces. Defaults to a temporary directory.
#' @param file Path to a JSON trace produced by the profiler.
#' @return `pliman_onnx_profiling()`: the trace prefix (invisibly). Already cached
#'   sessions are not profiled: call [pliman_clear_sessions()] first.
#'   `onnx_profile_summary()`: a data.frame with columns `provider`, `op_type`,
#'   `nodes`, `calls`, `total_ms` and `share` (fraction of the node time), ordered by
#'   time, with attribute `"providers"` summarising time per provider.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' pliman_clear_sessions()
#' prefix <- pliman_onnx_profiling(TRUE)
#' img <- image_pliman("soybean_touch.jpg")
#' res <- image_detect_dl(img, model = "yolo26n", plot = FALSE)
#' pliman_clear_sessions()                     # writes the trace
#' pliman_onnx_profiling(FALSE)
#' f <- list.files(dirname(prefix), "\\.json$", full.names = TRUE)
#' onnx_profile_summary(f[1])
pliman_onnx_profiling <- function(enable = TRUE, dir = file.path(tempdir(), "pliman_ort_profile")) {
  if (!isTRUE(enable)) {
    Sys.unsetenv("PLIMAN_ORT_PROFILE")
    return(invisible(NULL))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(normalizePath(dir, winslash = "/", mustWork = FALSE), "ort")
  Sys.setenv(PLIMAN_ORT_PROFILE = prefix)
  invisible(prefix)
}

#' @rdname pliman_onnx_profiling
#' @export
onnx_profile_summary <- function(file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg jsonlite} is required to read ONNX Runtime profiles.")
  }
  ev <- jsonlite::fromJSON(file, simplifyVector = TRUE)
  if (!is.data.frame(ev) || !"cat" %in% names(ev)) {
    cli::cli_abort("{.file {file}} does not look like an ONNX Runtime profile.")
  }
  nodes <- ev[ev$cat == "Node" & grepl("_kernel_time$", ev$name), , drop = FALSE]
  if (!nrow(nodes)) return(data.frame())
  prov <- nodes$args$provider
  op <- nodes$args$op_name
  node_name <- sub("_kernel_time$", "", nodes$name)
  if (is.null(prov)) prov <- rep(NA_character_, nrow(nodes))
  if (is.null(op)) op <- rep(NA_character_, nrow(nodes))
  prov <- sub("ExecutionProvider$", "", prov)
  key <- paste(prov, op, sep = "\r")
  dur <- as.numeric(nodes$dur) / 1000
  agg_ms <- tapply(dur, key, sum)
  agg_calls <- tapply(dur, key, length)
  agg_nodes <- tapply(node_name, key, function(x) length(unique(x)))
  parts <- strsplit(names(agg_ms), "\r", fixed = TRUE)
  out <- data.frame(
    provider = vapply(parts, `[`, character(1), 1L),
    op_type = vapply(parts, function(p) if (length(p) > 1L) p[2] else NA_character_, character(1)),
    nodes = as.integer(agg_nodes),
    calls = as.integer(agg_calls),
    total_ms = round(as.numeric(agg_ms), 3),
    stringsAsFactors = FALSE
  )
  tot <- sum(out$total_ms)
  out$share <- if (tot > 0) round(out$total_ms / tot, 4) else 0
  out <- out[order(-out$total_ms), , drop = FALSE]
  rownames(out) <- NULL
  by_prov <- tapply(out$total_ms, out$provider, sum)
  n_prov <- tapply(out$nodes, out$provider, sum)
  attr(out, "providers") <- data.frame(
    provider = names(by_prov),
    nodes = as.integer(n_prov[names(by_prov)]),
    total_ms = round(as.numeric(by_prov), 3),
    share = if (tot > 0) round(as.numeric(by_prov) / tot, 4) else 0,
    stringsAsFactors = FALSE
  )
  out
}

#' Cache optimized ONNX graphs on disk
#'
#' @description
#' When enabled, the first CPU session of each model saves the graph produced by ONNX
#' Runtime's optimizer to `dir`, and later sessions (new R sessions, parallel workers,
#' after [pliman_clear_sessions()]) load that optimized graph with the optimizer
#' turned off. This shortens session start-up for large models such as the SAM 2.1
#' encoder, Grounding DINO or RMBG. Predictions are unchanged: the cached file is the
#' same graph the optimizer would build again.
#'
#' Cached files are keyed by model size, modification time and ONNX Runtime API version,
#' so an updated model or runtime never reuses a stale graph. They are tied to the CPU
#' that produced them (the optimizer may insert CPU-specific layouts) and are therefore
#' kept locally; GPU sessions are never cached because provider-specific graphs are not
#' portable. Models stored with external weight files (`.onnx.data`) are skipped.
#'
#' @param enable Logical. `TRUE` enables the cache, `FALSE` disables it, `NULL`
#'   (default) only reports whether it is enabled.
#' @param dir Cache directory. Defaults to `file.path(pliman_model_dir(), "ort_graph_cache")`.
#' @param clear Logical. Delete the cached graphs in `dir`.
#' @return Logical, whether the cache was enabled before the call (invisibly when
#'   changing it).
#' @export
#' @examplesIf interactive()
#' pliman_onnx_graph_cache(TRUE)
#' pliman_session_info()[, c("model", "graph_cache", "create_ms")]
pliman_onnx_graph_cache <- function(enable = NULL,
                                    dir = file.path(pliman_model_dir(), "ort_graph_cache"),
                                    clear = FALSE) {
  prev <- nzchar(Sys.getenv("PLIMAN_ORT_OPT_CACHE", ""))
  if (isTRUE(clear)) {
    old <- list.files(dir, pattern = "\\.opt\\.onnx(\\.tmp[0-9]+)?$", full.names = TRUE)
    if (length(old)) unlink(old)
    cli::cli_alert_success("Removed {length(old)} cached graph{?s} from {.path {dir}}.")
  }
  if (is.null(enable)) return(if (isTRUE(clear)) invisible(prev) else prev)
  if (isTRUE(enable)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(PLIMAN_ORT_OPT_CACHE = normalizePath(dir, winslash = "/", mustWork = FALSE))
  } else {
    Sys.unsetenv("PLIMAN_ORT_OPT_CACHE")
  }
  invisible(prev)
}
