# ==============================================================================
# DEEP LEARNING RUNTIME CACHES
# pliman: Plant Image Analysis in R
#
# Package-level lookup tables and small session caches used on the hot path of the
# ONNX inference helpers (model-name aliases, input sizes, directory listings,
# resolved model paths, GPU adapters, Python environment). All caches live in
# `.onnx_meta_cache` (R/utils_neural.R) and are emptied by `pliman_clear_sessions()`.
# ==============================================================================

# User-friendly model aliases -> canonical model names (used by `.resolve_model_name()`).
.model_aliases <- c(
  "rmbg"                   = "rmbg-1.4",
  "bria"                   = "rmbg-1.4",
  "bria-rmbg"              = "rmbg-1.4",
  "bria_rmbg"              = "rmbg-1.4",
  "rmbg1.4"                = "rmbg-1.4",
  "rmbg_1.4"               = "rmbg-1.4",
  "rmbg2"                  = "rmbg-2.0",
  "rmbg-2"                 = "rmbg-2.0",
  "rmbg_2.0"               = "rmbg-2.0",
  "rmbg2.0"                = "rmbg-2.0",
  "bria-rmbg-2.0"          = "rmbg-2.0",
  "birefnet"               = "birefnet-lite",
  "birefnet_lite"          = "birefnet-lite",
  "without_bg"             = "withoutbg",
  "without-bg"             = "withoutbg",
  "withoutbg-open-weights" = "withoutbg",
  "u2net_p"                = "u2netp",
  "ben"                    = "ben2",
  "ben2"                   = "ben2",
  "ben-2"                  = "ben2",
  "boundary-aware"         = "ben2",
  "boundary_aware"         = "ben2",
  "boundary-aware-extraction-network" = "ben2",
  "sam"                    = "sam2.1",
  "sam2"                   = "sam2.1",
  "sam-2"                  = "sam2.1",
  "sam2.1"                 = "sam2.1",
  "sam-2.1"                = "sam2.1",
  "segment-anything"       = "sam2.1",
  "segment_anything"       = "sam2.1",
  "segmentanything"        = "sam2.1",
  "grounded-sam"           = "grounded-sam",
  "grounded_sam"           = "grounded-sam",
  "groundedsam"            = "grounded-sam",
  "grounding-dino"         = "grounded-sam",
  "grounding_dino"         = "grounded-sam",
  "groundingdino"          = "grounded-sam",
  "persam"                 = "persam",
  "per-sam"                = "persam",
  "persam2"                = "persam",
  "persam2.1"              = "persam",
  "depth-anything"         = "depth-anything-v2",
  "depthanything"          = "depth-anything-v2",
  "depth_anything"         = "depth-anything-v2",
  "depth-anything-v2"      = "depth-anything-v2",
  "depthanythingv2"        = "depth-anything-v2",
  "depth_anything_v2"      = "depth-anything-v2",
  "depth"                  = "depth-anything-v2",
  "dinov2"                 = "dinov2",
  "dino-v2"                = "dinov2",
  "dino_v2"                = "dinov2",
  "dino"                   = "dinov2",
  # YOLO Object Detection
  "yolo"                   = "yolo26n",
  "yolo26"                 = "yolo26n",
  "yolo-26"                = "yolo26n",
  "yolo26n"                = "yolo26n",
  "yolo-26n"               = "yolo26n",
  "yolo26s"                = "yolo26s",
  "yolo-26s"               = "yolo26s",
  "yolo26m"                = "yolo26m",
  "yolo-26m"               = "yolo26m",
  "yolo26l"                = "yolo26l",
  "yolo-26l"               = "yolo26l",
  "yolo26x"                = "yolo26x",
  "yolo-26x"               = "yolo26x",
  "yolo11"                 = "yolo26n",
  "yolo-11"                = "yolo26n",
  "yolo11n"                = "yolo26n",
  "yolo-11n"               = "yolo26n",
  "yolov8"                 = "yolo26n",
  "yolo8"                  = "yolo26n",

  # YOLO Instance Segmentation
  "yolo-seg"               = "yolo26n-seg",
  "yoloseg"                = "yolo26n-seg",
  "yolo26-seg"             = "yolo26n-seg",
  "yolo-26-seg"            = "yolo26n-seg",
  "yolo26n-seg"            = "yolo26n-seg",
  "yolo26s-seg"            = "yolo26s-seg",
  "yolo26m-seg"            = "yolo26m-seg",
  "yolo26l-seg"            = "yolo26l-seg",
  "yolo26x-seg"            = "yolo26x-seg",
  "yolo11-seg"             = "yolo26n-seg",
  "yolo11n-seg"            = "yolo26n-seg",

  # YOLO Pose Estimation
  "yolo-pose"              = "yolo26n-pose",
  "yolopose"               = "yolo26n-pose",
  "yolo26-pose"            = "yolo26n-pose",
  "yolo-26-pose"           = "yolo26n-pose",
  "yolo26n-pose"           = "yolo26n-pose",
  "yolo26s-pose"           = "yolo26s-pose",
  "yolo26m-pose"           = "yolo26m-pose",
  "yolo26l-pose"           = "yolo26l-pose",
  "yolo26x-pose"           = "yolo26x-pose",

  # YOLO Classification
  "yolo-cls"               = "yolo26n-cls",
  "yolocls"                = "yolo26n-cls",
  "yolo26-cls"             = "yolo26n-cls",
  "yolo-26-cls"            = "yolo26n-cls",
  "yolo26n-cls"            = "yolo26n-cls",
  "yolo26s-cls"            = "yolo26s-cls",
  "yolo26m-cls"            = "yolo26m-cls",
  "yolo26l-cls"            = "yolo26l-cls",
  "yolo26x-cls"            = "yolo26x-cls",
  "yolo-26n-cls"           = "yolo26n-cls",
  "yolo-26s-cls"           = "yolo26s-cls",
  "yolo-26m-cls"           = "yolo26m-cls",
  "yolo-26l-cls"           = "yolo26l-cls",
  "yolo-26x-cls"           = "yolo26x-cls",
  "yolo11-cls"             = "yolo26n-cls",
  "yolo11n-cls"            = "yolo26n-cls",
  "yolo11s-cls"            = "yolo26s-cls",
  "yolo11m-cls"            = "yolo26m-cls",
  "yolo11l-cls"            = "yolo26l-cls",
  "yolo11x-cls"            = "yolo26x-cls",
  "yolo-11-cls"            = "yolo26n-cls",
  "yolo-11n-cls"           = "yolo26n-cls",
  "yolo-11s-cls"           = "yolo26s-cls",
  "yolo-11m-cls"           = "yolo26m-cls",
  "yolo-11l-cls"           = "yolo26l-cls",
  "yolo-11x-cls"           = "yolo26x-cls",
  "yolov8-cls"             = "yolo26n-cls",
  "yolov8n-cls"            = "yolo26n-cls",
  "yolov8s-cls"            = "yolo26s-cls",
  "yolov8m-cls"            = "yolo26m-cls",
  "yolov8l-cls"            = "yolo26l-cls",
  "yolov8x-cls"            = "yolo26x-cls",
  "yolo8-cls"              = "yolo26n-cls",
  "yolo8n-cls"             = "yolo26n-cls",
  "yolo8s-cls"             = "yolo26s-cls",
  "yolo8m-cls"             = "yolo26m-cls",
  "yolo8l-cls"             = "yolo26l-cls",
  "yolo8x-cls"             = "yolo26x-cls",
  "stardist"               = "stardist",
  "star-dist"              = "stardist",
  "stardist-dsb2018"       = "stardist",
  "star_dist"              = "stardist",
  "purescale"              = "SAFMN4x",
  "purescale-4x"           = "SAFMN4x",
  "4x-purescale"           = "SAFMN4x",
  "4x_purescale"           = "SAFMN4x",
  "SAFMN4x"                = "SAFMN4x",
  "safmn"                  = "SAFMN4x",
  "safmn-4x"               = "SAFMN4x",
  "safmn4x"                = "SAFMN4x",
  "safmn_4x"               = "SAFMN4x",
  "4x_safmn"               = "SAFMN4x",
  "4x-safmn"               = "SAFMN4x",
  "SAFMN4x.onnx"           = "SAFMN4x",
  "purescale-2x"           = "SAFMN2x",
  "2x-purescale"           = "SAFMN2x",
  "2x_purescale"           = "SAFMN2x",
  "SAFMN2x"                = "SAFMN2x",
  "safmn-2x"               = "SAFMN2x",
  "safmn2x"                = "SAFMN2x",
  "safmn_2x"               = "SAFMN2x",
  "2x_safmn"               = "SAFMN2x",
  "2x-safmn"               = "SAFMN2x",
  "SAFMN2x.onnx"           = "SAFMN2x",
  "purescale-2x-sharper"   = "SAFMN2x_sharper",
  "2x-purescale-sharper"   = "SAFMN2x_sharper",
  "2x_purescale_sharper"   = "SAFMN2x_sharper",
  "SAFMN2x_sharper"        = "SAFMN2x_sharper",
  "safmn-2x-sharper"       = "SAFMN2x_sharper",
  "safmn2x_sharper"        = "SAFMN2x_sharper",
  "safmn_2x_sharper"       = "SAFMN2x_sharper",
  "2x_safmn_sharper"       = "SAFMN2x_sharper",
  "2x-safmn-sharper"       = "SAFMN2x_sharper",
  "SAFMN2x_sharper.onnx"   = "SAFMN2x_sharper",
  "realesrgan"             = "realesrgan_x2_full",
  "real-esrgan"            = "realesrgan_x2_full",
  "realesrgan-compact"     = "realesrgan",
  "realesrgan_compact"     = "realesrgan",
  "real-esrgan-compact"    = "realesrgan",
  "realesrgan-x2"          = "realesrgan_x2_full",
  "realesrgan_x2"          = "realesrgan_x2_full",
  "real-esrgan-x2"         = "realesrgan_x2_full",
  "real_esrgan_x2"         = "realesrgan_x2_full",
  "realesrgan-x2-full"     = "realesrgan_x2_full",
  "realesrgan_x2_full"     = "realesrgan_x2_full",
  "real-esrgan-x2-full"    = "realesrgan_x2_full",
  "real_esrgan_x2_full"    = "realesrgan_x2_full",
  "realesrgan-x4"          = "realesrgan_x4_full",
  "realesrgan_x4"          = "realesrgan_x4_full",
  "real-esrgan-x4"         = "realesrgan_x4_full",
  "real_esrgan_x4"         = "realesrgan_x4_full",
  "realesrgan-x4-full"     = "realesrgan_x4_full",
  "realesrgan_x4_full"     = "realesrgan_x4_full",
  "real-esrgan-x4-full"    = "realesrgan_x4_full",
  "real_esrgan_x4_full"    = "realesrgan_x4_full",
  "superres"               = "realesrgan_x2_full",
  "super-res"              = "realesrgan_x2_full",
  "clip"                   = "clip-vit-b32",
  "clip-vit-b32"           = "clip-vit-b32",
  "clip-b32"               = "clip-vit-b32",
  "clip-vit-base"          = "clip-vit-b32",
  "clip_vit_b32"           = "clip-vit-b32",
  "yolo-world"             = "yolov8s-worldv2",
  "yoloworld"              = "yolov8s-worldv2",
  "yolo-world-v2"          = "yolov8s-worldv2",
  "yolov8s-world"          = "yolov8s-worldv2",
  "yolov8sworld"           = "yolov8s-worldv2",
  "yolov8s-worldv2"        = "yolov8s-worldv2",
  "yolov8m-world"          = "yolov8m-worldv2",
  "yolov8mworld"           = "yolov8m-worldv2",
  "yolov8m-worldv2"        = "yolov8m-worldv2",
  "yolov8l-world"          = "yolov8l-worldv2",
  "yolov8lworld"           = "yolov8l-worldv2",
  "yolov8l-worldv2"        = "yolov8l-worldv2",
  "yolov8x-world"          = "yolov8x-worldv2",
  "yolov8xworld"           = "yolov8x-worldv2",
  "yolov8x-worldv2"        = "yolov8x-worldv2",
  "yolo-nas-s"             = "yolo_nas_s",
  "yolonas-s"              = "yolo_nas_s",
  "yolonas_s"              = "yolo_nas_s",
  "yolo_nas_s"             = "yolo_nas_s",
  "yolo-nas-l"             = "yolo_nas_l",
  "yolonas-l"              = "yolo_nas_l",
  "yolonas_l"              = "yolo_nas_l",
  "yolo_nas_l"             = "yolo_nas_l",
  "yolo-nas"               = "yolo_nas_s",
  "yolonas"                = "yolo_nas_s",
  # YOLO Oriented Bounding Box (OBB)
  "yolo-obb"               = "yolo26n-obb",
  "yoloobb"                = "yolo26n-obb",
  "yolo26-obb"             = "yolo26n-obb",
  "yolo-26-obb"            = "yolo26n-obb",
  "yolo26n-obb"            = "yolo26n-obb",
  "yolo-26n-obb"           = "yolo26n-obb",
  "yolo26s-obb"            = "yolo26s-obb",
  "yolo-26s-obb"           = "yolo26s-obb",
  "yolo26m-obb"            = "yolo26m-obb",
  "yolo-26m-obb"           = "yolo26m-obb",
  "yolo26l-obb"            = "yolo26l-obb",
  "yolo-26l-obb"           = "yolo26l-obb",
  "yolo26x-obb"            = "yolo26x-obb",
  "yolo-26x-obb"           = "yolo26x-obb",
  "yolo11-obb"             = "yolo26n-obb",
  "yolo-11-obb"            = "yolo26n-obb",
  "yolo11n-obb"            = "yolo26n-obb",
  "yolo-11n-obb"           = "yolo26n-obb",
  "yolo11s-obb"            = "yolo26s-obb",
  "yolo-11s-obb"           = "yolo26s-obb",
  "yolo11m-obb"            = "yolo26m-obb",
  "yolo-11m-obb"           = "yolo26m-obb",
  "yolo11l-obb"            = "yolo26l-obb",
  "yolo-11l-obb"           = "yolo26l-obb",
  "yolo11x-obb"            = "yolo26x-obb",
  "yolo-11x-obb"           = "yolo26x-obb",
  "yolov8-obb"             = "yolo26n-obb",
  "yolo8-obb"              = "yolo26n-obb",
  "yolov8n-obb"            = "yolo26n-obb",
  "yolo8n-obb"             = "yolo26n-obb",
  "yolov8s-obb"            = "yolo26s-obb",
  "yolo8s-obb"             = "yolo26s-obb",
  "yolov8m-obb"            = "yolo26m-obb",
  "yolo8m-obb"             = "yolo26m-obb",
  "yolov8l-obb"            = "yolo26l-obb",
  "yolo8l-obb"             = "yolo26l-obb",
  "yolov8x-obb"            = "yolo26x-obb",
  "yolo8x-obb"             = "yolo26x-obb",
  "sd-turbo"               = "sd-turbo",
  "sdturbo"                = "sd-turbo",
  "sd_turbo"               = "sd-turbo",
  "lcm"                    = "sd-turbo",
  "ifcnn"                  = "ifcnn",
  "ifcnn-max"              = "ifcnn",
  "ifcnn_max"              = "ifcnn",
  "mfif"                   = "ifcnn",
  "focus-stack"            = "ifcnn",
  # Optical Character Recognition (PP-OCRv4 / DBNet)
  "ppocr"                  = "ch_PP-OCRv4_det",
  "pp-ocr"                 = "ch_PP-OCRv4_det",
  "ppocr_det"              = "ch_PP-OCRv4_det",
  "ppocr-det"              = "ch_PP-OCRv4_det",
  "ppocr_rec"              = "ch_PP-OCRv4_rec",
  "ppocr-rec"              = "ch_PP-OCRv4_rec",
  "ch_pp-ocrv4_det"        = "ch_PP-OCRv4_det",
  "ch_pp-ocrv4-det"        = "ch_PP-OCRv4_det",
  "ch_pp-ocrv4_rec"        = "ch_PP-OCRv4_rec",
  "ch_pp-ocrv4-rec"        = "ch_PP-OCRv4_rec",
  "ppocr_det_server"       = "ch_PP-OCRv4_det_server",
  "ppocr-det-server"       = "ch_PP-OCRv4_det_server",
  "ppocr_server"           = "ch_PP-OCRv4_det_server",
  "ch_pp-ocrv4_det_server" = "ch_PP-OCRv4_det_server",
  "ch_pp-ocrv4-det-server" = "ch_PP-OCRv4_det_server",
  "craft"                  = "craft",
  "craft_mlt"              = "craft",
  "craft-mlt"              = "craft",
  "trocr"                  = "trocr_small",
  "trocr_small"            = "trocr_small",
  "trocr-small"            = "trocr_small",
  "trocr_lat"              = "trocr_small_lat",
  "trocr_small_lat"        = "trocr_small_lat",
  "trocr-small-lat"        = "trocr_small_lat",
  "trocr_encoder"          = "trocr_encoder",
  "trocr-encoder"          = "trocr_encoder",
  "trocr_decoder"          = "trocr_decoder",
  "trocr-decoder"          = "trocr_decoder",
  "trocr_vocab"            = "trocr_vocab",
  "trocr-vocab"            = "trocr_vocab",
  "trocr_dict"             = "trocr_vocab",
  "trocr_encoder_lat"      = "trocr_encoder_lat",
  "trocr-encoder-lat"      = "trocr_encoder_lat",
  "trocr_decoder_lat"      = "trocr_decoder_lat",
  "trocr-decoder-lat"      = "trocr_decoder_lat",
  "trocr_vocab_lat"        = "trocr_vocab_lat",
  "trocr-vocab-lat"        = "trocr_vocab_lat",
  "trocr_dict_lat"         = "trocr_vocab_lat",
  "tocr_encoder_lat"       = "trocr_encoder_lat",
  "tocr_decoder_lat"       = "trocr_decoder_lat",
  "tocr_vocab_lat"         = "trocr_vocab_lat",
  # Image Restoration (Zero-DCE & NAFNet)
  "zero_dce"               = "zero_dce",
  "zero-dce"               = "zero_dce",
  "zerodce"                = "zero_dce",
  "lowlight"               = "zero_dce",
  "low-light"              = "zero_dce",
  "nafnet"                 = "nafnet",
  "naf-net"                = "nafnet",
  "deblur"                 = "nafnet",
  "denoise"                = "nafnet"
)

# Native input resolution of the salient-object / foundation models that go through
# `.compute_dl_mask()` and friends. Previously this was read from
# `pliman_available_models()` on every call, which rebuilds the whole model catalogue
# (and may touch the network when the asset cache is stale).
.dl_model_input_sizes <- c(
  "u2netp"                  = 320L,
  "u2net"                   = 320L,
  "silueta"                 = 320L,
  "birefnet-lite"           = 1024L,
  "isnet-general-use"       = 1024L,
  "rmbg-1.4"                = 1024L,
  "rmbg-2.0"                = 1024L,
  "withoutbg"               = 448L,
  "ben2"                    = 1024L,
  "sam2.1"                  = 1024L,
  "grounded-sam"            = 800L,
  "groundingdino-tiny"      = 800L,
  "depth-anything-v2-small" = 518L,
  "dinov2"                  = 518L,
  "dinov2-vits14"           = 518L
)

# Input size for a catalogue model: static table first, then (once per session) the
# model registry, falling back to `default`.
.dl_input_size <- function(model, dir = pliman_model_dir(), default = 320L) {
  hit <- .dl_model_input_sizes[model]
  if (!is.na(hit)) return(unname(hit))
  key <- paste0("input_size::", model)
  if (exists(key, envir = .onnx_meta_cache, inherits = FALSE)) {
    return(get(key, envir = .onnx_meta_cache, inherits = FALSE))
  }
  val <- tryCatch({
    models_df <- pliman_available_models(dir = dir)
    row <- models_df[models_df$model == model, ]
    if (nrow(row) > 0) as.integer(strsplit(as.character(row$input_size[1]), "x")[[1]][1]) else default
  }, error = function(e) default)
  if (is.na(val)) val <- default
  assign(key, val, envir = .onnx_meta_cache)
  val
}

# Directory listing of model files cached by directory modification time, so repeated
# fuzzy lookups do not call list.files() on every inference.
.cached_model_listing <- function(dir) {
  mt <- suppressWarnings(file.info(dir)$mtime)
  key <- paste0("listing::", dir)
  if (exists(key, envir = .onnx_meta_cache, inherits = FALSE)) {
    hit <- get(key, envir = .onnx_meta_cache, inherits = FALSE)
    if (identical(hit$mtime, mt)) return(hit$files)
  }
  files <- list.files(dir, pattern = "\\.(onnx|pt)$", ignore.case = TRUE, full.names = FALSE)
  assign(key, list(mtime = mt, files = files), envir = .onnx_meta_cache)
  files
}

# Negative-lookup cache for `.find_local_model()`: remembers names that did not resolve,
# keyed by the directory mtime (adding/removing a model invalidates it).
.local_model_miss <- function(dir, name, set = FALSE) {
  key <- paste0("local_miss::", dir, "::", name)
  mt <- suppressWarnings(file.info(dir)$mtime)
  if (isTRUE(set)) {
    assign(key, mt, envir = .onnx_meta_cache)
    return(invisible(TRUE))
  }
  if (!exists(key, envir = .onnx_meta_cache, inherits = FALSE)) return(FALSE)
  identical(get(key, envir = .onnx_meta_cache, inherits = FALSE), mt)
}

# Resolved-path cache for `pliman_download_model()`, so inference helpers that call it on
# every image (U2Net, SAM, Grounding DINO, StarDist, DINOv2, ...) return immediately.
.dl_path_cache_get <- function(dir, model) {
  key <- paste0("dl::", dir, "::", paste(model, collapse = "|"))
  if (!exists(key, envir = .onnx_meta_cache, inherits = FALSE)) return(NULL)
  hit <- get(key, envir = .onnx_meta_cache, inherits = FALSE)
  if (length(hit) > 0 && all(file.exists(hit))) hit else NULL
}

.dl_path_cache_set <- function(dir, model, path) {
  if (is.null(path) || length(path) == 0 || !all(file.exists(path))) return(invisible(path))
  key <- paste0("dl::", dir, "::", paste(model, collapse = "|"))
  assign(key, path, envir = .onnx_meta_cache)
  invisible(path)
}

# GPU adapter enumeration (DXGI) is cached for the session.
.gpu_info_cached <- function() {
  if (exists("gpu_info", envir = .onnx_meta_cache, inherits = FALSE)) {
    return(get("gpu_info", envir = .onnx_meta_cache, inherits = FALSE))
  }
  info <- tryCatch(pliman_gpu_info_cpp(), error = function(e) NULL)
  assign("gpu_info", info, envir = .onnx_meta_cache)
  info
}

# Python + ultralytics environment check used by YOLO-World / generative helpers.
# Spawning Python to verify imports costs seconds, so do it once per session.
.python_yolo_env_cached <- function(auto_install = TRUE) {
  if (exists("py_env_exec", envir = .onnx_meta_cache, inherits = FALSE)) {
    py <- get("py_env_exec", envir = .onnx_meta_cache, inherits = FALSE)
    if (!is.null(py) && nzchar(py) && (file.exists(py) || nzchar(Sys.which(py)))) return(py)
  }
  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)
  assign("py_env_exec", py_exec, envir = .onnx_meta_cache)
  py_exec
}
