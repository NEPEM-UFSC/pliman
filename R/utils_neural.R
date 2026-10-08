# utils_neural.R - Pre-trained Deep Learning Segmentation and Background Removal for pliman
#
# ==============================================================================
# DEEP LEARNING ARCHITECTURE & ONNX RUNTIME DOCUMENTATION
# ==============================================================================
#' @title Deep Learning Background Removal & Salient Object Segmentation in pliman
#' @name utils_neural
#' @description
#' The `pliman` neural module provides high-accuracy, zero-Python foreground segmentation
#' and background removal powered by pre-trained Deep Learning models in **ONNX**
#' (Open Neural Network Exchange) format.
#'
#' @details
#' # 1. Supported Neural Architectures
#' * **`"u2netp"` (Default for CPU)**:
#'   * **Size:** ~4.4 MB | **Input Resolution:** 320 x 320
#'   * **Architecture:** Ultra-lightweight U2-Net with nested U-structure residual blocks.
#'   * **Best For:** Fast batch processing on standard laptops and CPUs without dedicated GPU.
#'
#' * **`"ben2"` (Boundary-Aware Extraction Network)**:
#'   * **Size:** ~212.6 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** BEN2 Base (Prama LLC) with boundary-focused refinement stream.
#'   * **Best For:** Crisp sub-pixel object boundaries, fine petioles, leaf margins, roots, and hairs.
#'
#' * **`"grounded-sam"` (Zero-Shot Text-Prompted Detection & Instance Segmentation)**:
#'   * **Size:** ~314.8 MB (Grounding DINO ~194.4 MB + SAM 2.1 ~120.2 MB) | **Input Resolution:** 800 x 800 (DINO) & 1024 x 1024 (SAM)
#'   * **Architecture:** Combines Grounding DINO Tiny for open-vocabulary text-prompted bounding box detection with SAM 2.1 Hiera-Tiny for instance mask extraction.
#'   * **Best For:** Detecting and segmenting specific objects using natural language prompts (e.g. `"cow"`, `"leaf"`, `"lesion"`), with support for bounding-box-only detection (`mask = FALSE, bbox = TRUE`).
#'
#' * **`"sam2.1"` (Segment Anything 2.1 by Meta AI, alias `"sam2"`)**:
#'   * **Size:** ~120.2 MB (Encoder + Decoder) | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** Meta AI Hiera-Tiny foundation vision transformer with promptable mask decoder.
#'   * **Best For:** Zero-shot foundation segmentation with point prompts (`pick_object = TRUE`), box prompts, or center point defaults.
#'
#' * **`"sam3.1"` (Segment Anything 3.1 by Meta AI, alias `"sam3"`)**:
#'   * **Size:** ~868.1 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** Meta AI Segment Anything 3.1 with real-time concept-driven segmentation.
#'   * **Best For:** Advanced foundation segmentation across diverse biological structures.
#'
#' * **`"birefnet-lite"` (Bilateral Reference Network Lite)**:
#'   * **Size:** ~213.6 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** Bilateral Reference Network for High-Resolution Dichotomous Image Segmentation (BiRefNet Lite).
#'   * **Best For:** Extremely fine margins, hair-thin leaf serrations, complex petioles, translucent halos, and subtle lesions.
#'
#' * **`"isnet-general-use"` (High Precision Boundary Matting)**:
#'   * **Size:** ~170.4 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** IS-Net (Intermediate Supervision Network) specialized in high-resolution boundary matting.
#'   * **Best For:** Plant leaves, fine serrated margins, lesion borders, chlorotic halos, veins, and complex natural backgrounds.
#'
#' * **`"rmbg-1.4"` (BRIA AI Standard, alias `"rmbg"`)**:
#'   * **Size:** ~168.0 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** BRIA RMBG 1.4 trained on extensive curated foreground datasets.
#'   * **Best For:** General background removal under challenging lighting, shadows, reflections, and fine object edges.
#'
#' * **`"rmbg-2.0"` (BRIA AI Next-Gen Model)**:
#'   * **Size:** ~976.9 MB | **Input Resolution:** 1024 x 1024
#'   * **Architecture:** BRIA RMBG 2.0 based on a BiRefNet backbone trained on extensive curated commercial datasets.
#'   * **Best For:** Highest-fidelity boundary matting, transparent surfaces, fine plant hairs, and complex scenes.
#'
#' * **`"withoutbg"` (State-of-the-Art DepthAnythingV2 + ConvNeXt Matting)**:
#'   * **Size:** ~433.4 MB | **Input Resolution:** 448 x 448
#'   * **Architecture:** withoutBG Open Weights: DINOv3 ConvNeXt-fused U-Net guided by DepthAnythingV2 depth features.
#'   * **Best For:** Complex 3D depth-separated backgrounds, transparent materials, and fine details.
#'
#' * **`"silueta"`**:
#'   * **Size:** ~42.1 MB | **Input Resolution:** 320 x 320
#'   * **Architecture:** Compact silhouette extraction network.
#'   * **Best For:** Balanced speed-accuracy tradeoff on edge/CPU hardware.
#'
#' * **`"u2net"`**:
#'   * **Size:** ~167.8 MB | **Input Resolution:** 320 x 320
#'   * **Architecture:** Full-depth U2-Net salient object detection model.
#'   * **Best For:** Robust multi-scale salient object detection on standard resolution images.
#'
#' # 2. Zero-Python Architecture (Native C++ Runtime)
#'
#' Unlike standard Python-dependent deep learning packages, `pliman` executes ONNX models
#' natively in C++ via direct Microsoft ONNX Runtime C API dynamic bindings.
#'
#' * **Zero Python Requirement:** No Python, Conda, Virtualenv, or pip required.
#' * **Direct C++ Binding:** Connects directly to Microsoft's official ONNX Runtime C++ engine (`onnxruntime.dll` on Windows / `libonnxruntime.so` on Linux / `libonnxruntime.dylib` on macOS).
#' * **CPU & GPU Support:** Runs optimized CPU multi-threading (AVX2/AVX512) or CUDA GPU acceleration.
#'
#' # 3. Model Weight Storage & Cache
#'
#' Downloaded `.onnx` model weights are stored permanently in the user's R data directory:
#' `tools::R_user_dir("pliman", which = "data")/models`
#' (e.g., `C:/Users/<User>/AppData/Local/R/data/R/pliman/models` on Windows).
#' This follows CRAN policy and ensures models are preserved across R session restarts
#' and package upgrades.
#'
#' # 4. Quick Setup Guide
#'
#' To configure everything in a single call:
#' ```r
#' library(pliman)
#'
#' # 1. Set up the environment (downloads ONNX Runtime C++ library & default models)
#' pliman_configure_dl()
#'
#' # 2. Extract a binary mask
#' mask <- image_binary_dl(img, model = "u2netp")
#'
#' # 3. Remove background with transparent alpha channel (RGBA)
#' img_trans <- image_remove_bg_dl(img, model = "isnet-general-use")
#'
#' # 4. General segmentation
#' seg <- image_segment_dl(img, model = "u2netp")
#' ```
NULL

# ==============================================================================
# SECTION 1: MODEL MANAGEMENT & DIRECTORY MANAGEMENT
# ==============================================================================

#' Model Directory Management for pliman Deep Learning Models
#'
#' Gets or sets the local directory where pre-trained ONNX models are stored.
#' Defaults to the standard R user data directory:
#' `tools::R_user_dir("pliman", which = "data")/models`.
#'
#' @param dir Optional character string specifying a custom directory path.
#'   If provided, creates the directory if it does not exist.
#' @return A character string with the normalized absolute path to the model directory.
#' @export
#' @examples
#' pliman_model_dir()
pliman_model_dir <- function(dir = NULL) {
  if (!is.null(dir)) {
    dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(dir)
  }

  default_dir <- file.path(tools::R_user_dir("pliman", which = "data"), "models")
  default_dir <- normalizePath(default_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(default_dir)) {
    dir.create(default_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(default_dir)
}

# Internal helper to resolve user-friendly model aliases
.resolve_model_name <- function(model) {
  aliases <- c(
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
    "sam3"                   = "sam3.1",
    "sam-3"                  = "sam3.1",
    "sam3.1"                 = "sam3.1",
    "sam-3.1"                = "sam3.1",
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
    "realesrgan"             = "realesrgan-compact",
    "real-esrgan"            = "realesrgan-compact",
    "realesrgan-compact"     = "realesrgan-compact",
    "real_esrgan"            = "realesrgan-compact",
    "superres"               = "realesrgan-compact",
    "super-res"              = "realesrgan-compact",
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
    "focus-stack"            = "ifcnn"
  )
  vapply(model, function(m) {
    ml <- tolower(m)
    if (ml %in% names(aliases)) aliases[[ml]] else m
  }, character(1), USE.NAMES = FALSE)
}

.coco_classes <- c(
  "person", "bicycle", "car", "motorcycle", "airplane", "bus", "train", "truck", "boat", "traffic light",
  "fire hydrant", "stop sign", "parking meter", "bench", "bird", "cat", "dog", "horse", "sheep", "cow",
  "elephant", "bear", "zebra", "giraffe", "backpack", "umbrella", "handbag", "tie", "suitcase", "frisbee",
  "skis", "snowboard", "sports ball", "kite", "baseball bat", "baseball glove", "skateboard", "surfboard", "tennis racket", "bottle",
  "wine glass", "cup", "fork", "knife", "spoon", "bowl", "banana", "apple", "sandwich", "orange",
  "broccoli", "carrot", "hot dog", "pizza", "donut", "cake", "chair", "couch", "potted plant", "bed",
  "dining table", "toilet", "tv", "laptop", "mouse", "remote", "keyboard", "cell phone", "microwave", "oven",
  "toaster", "sink", "refrigerator", "book", "clock", "vase", "scissors", "teddy bear", "hair drier", "toothbrush"
)

.dota_classes <- c(
  "plane", "ship", "storage tank", "baseball diamond", "tennis court",
  "basketball court", "ground track field", "harbor", "bridge", "large vehicle",
  "small vehicle", "helicopter", "roundabout", "soccer ball field", "swimming pool"
)

# Helper to automatically extract embedded class names from ONNX model metadata
.get_onnx_classes <- function(model_file) {
  if (is.null(model_file) || !file.exists(model_file)) return(NULL)
  tryCatch({
    sz <- file.info(model_file)$size
    if (is.na(sz) || sz < 1024) return(NULL)
    read_sz <- min(sz, 262144L)
    con <- file(model_file, "rb")
    seek(con, max(0, sz - read_sz))
    chunk <- readBin(con, "raw", n = read_sz)
    close(con)

    txt <- rawToChar(chunk[chunk >= 32 & chunk <= 126])
    pos <- regexpr("names[^{]*\\{([^}]+)\\}", txt)
    if (pos > 0) {
      match_full <- substr(txt, pos, pos + attr(pos, "match.length") - 1)
      inner <- sub("^[^{]*\\{", "", match_full)
      inner <- sub("\\}[^}]*$", "", inner)
      matches <- gregexpr("(\\d+)\\s*:\\s*(['\"])(.*?)\\2", inner, perl = TRUE)
      if (matches[[1]][1] > 0) {
        matched_tokens <- regmatches(inner, matches)[[1]]
        indices <- as.integer(sub("^(\\d+)\\s*:.*$", "\\1", matched_tokens, perl = TRUE))
        class_names <- sub("^\\d+\\s*:\\s*['\"]", "", matched_tokens, perl = TRUE)
        class_names <- sub("['\"]\\s*$", "", class_names, perl = TRUE)
        ord <- order(indices)
        return(class_names[ord])
      }
    }
    NULL
  }, error = function(e) NULL)
}
# Helper to map colors by object class/label or rainbow per instance (persistent by ID when available)
.get_class_palette <- function(labels = NULL, n = length(labels), user_col = NULL, default_col = "#00CC66", rainbow = TRUE, ids = NULL) {
  if ((is.null(labels) || length(labels) == 0L) && n <= 0L) {
    return(character(0))
  }
  if (is.null(labels) || length(labels) == 0L) {
    chr_labels <- rep("object", max(1L, n))
  } else {
    chr_labels <- as.character(labels)
    chr_labels[is.na(chr_labels) | !nzchar(chr_labels)] <- "object"
  }

  # Curated modern vibrant palette for multi-class detection, tracking & rainbow mode
  palette_pool <- c(
    "#00CC66", "#0099FF", "#FF9900", "#FF3366", "#9933FF",
    "#00E5FF", "#FFD700", "#FF6600", "#00B4D8", "#7209B7",
    "#2EC4B6", "#E71D36", "#FF9F1C", "#4361EE", "#F72585",
    "#10B981", "#3B82F6", "#8B5CF6", "#EC4899", "#F59E0B",
    "#14B8A6", "#6366F1", "#A855F7", "#D946EF", "#EF4444",
    "#84CC16", "#06B6D4", "#F97316", "#FB7185", "#38BDF8",
    "#22C55E", "#A855F7", "#E11D48", "#EAB308", "#0284C7",
    "#16A34A", "#9333EA", "#C026D3", "#EA580C", "#0D9488"
  )

  # If rainbow = TRUE, assign a distinct color to each object instance.
  # When tracking IDs are available, the color is strictly locked to the object's ID
  # so that it NEVER changes across frames!
  if (isTRUE(rainbow)) {
    pool <- if (!is.null(user_col) && length(user_col) > 0L) user_col else palette_pool
    n_pool <- length(pool)

    if (!is.null(ids) && length(ids) == length(chr_labels)) {
      id_num <- suppressWarnings(as.integer(ids))
      cols <- vapply(seq_along(id_num), function(idx) {
        val <- id_num[idx]
        if (!is.na(val) && val > 0L) {
          pool[((val - 1L) %% n_pool) + 1L]
        } else {
          pool[((idx - 1L) %% n_pool) + 1L]
        }
      }, character(1))
      return(cols)
    }

    n_inst <- length(chr_labels)
    if (n_inst <= n_pool) {
      return(pool[seq_len(n_inst)])
    } else {
      return(grDevices::rainbow(n_inst, s = 0.85, v = 0.95))
    }
  }

  u_labels <- unique(chr_labels)
  n_cls <- length(u_labels)

  if (!is.null(user_col)) {
    if (length(user_col) == 1L) {
      cls_map <- stats::setNames(rep(user_col, n_cls), u_labels)
    } else {
      cls_map <- stats::setNames(rep_len(user_col, n_cls), u_labels)
    }
  } else if (n_cls == 1L) {
    cls_map <- stats::setNames(default_col, u_labels)
  } else if (n_cls <= length(palette_pool)) {
    cls_map <- stats::setNames(palette_pool[seq_len(n_cls)], u_labels)
  } else {
    cls_map <- stats::setNames(grDevices::rainbow(n_cls, s = 0.85, v = 0.95), u_labels)
  }

  unname(cls_map[chr_labels])
}

.coco_keypoints <- c(
  "nose", "left_eye", "right_eye", "left_ear", "right_ear",
  "left_shoulder", "right_shoulder", "left_elbow", "right_elbow",
  "left_wrist", "right_wrist", "left_hip", "right_hip",
  "left_knee", "right_knee", "left_ankle", "right_ankle"
)

# Standard COCO 17-keypoint anatomical skeleton pairs (1-indexed)
.coco_skeleton_pairs <- matrix(c(
  1, 2,  1, 3,  2, 4,  3, 5,       # facial (nose-eye-ear)
  6, 7,  6, 12, 7, 13, 12, 13,     # torso (shoulders & hips)
  6, 8,  8, 10,                    # left arm (shoulder-elbow-wrist)
  7, 9,  9, 11,                    # right arm (shoulder-elbow-wrist)
  12, 14, 14, 16,                  # left leg (hip-knee-ankle)
  13, 15, 15, 17                   # right leg (hip-knee-ankle)
), ncol = 2, byrow = TRUE)

.imagenet_classes <- c(
  'tench', 'goldfish', 'great_white_shark', 'tiger_shark', 'hammerhead', 'electric_ray', 'stingray', 'cock', 'hen', 'ostrich',
  'brambling', 'goldfinch', 'house_finch', 'junco', 'indigo_bunting', 'robin', 'bulbul', 'jay', 'magpie', 'chickadee',
  'water_ouzel', 'kite', 'bald_eagle', 'vulture', 'great_grey_owl', 'European_fire_salamander', 'common_newt', 'eft', 'spotted_salamander', 'axolotl',
  'bullfrog', 'tree_frog', 'tailed_frog', 'loggerhead', 'leatherback_turtle', 'mud_turtle', 'terrapin', 'box_turtle', 'banded_gecko', 'common_iguana',
  'American_chameleon', 'whiptail', 'agama', 'frilled_lizard', 'alligator_lizard', 'Gila_monster', 'green_lizard', 'African_chameleon', 'Komodo_dragon', 'African_crocodile',
  'American_alligator', 'triceratops', 'thunder_snake', 'ringneck_snake', 'hognose_snake', 'green_snake', 'king_snake', 'garter_snake', 'water_snake', 'vine_snake',
  'night_snake', 'boa_constrictor', 'rock_python', 'Indian_cobra', 'green_mamba', 'sea_snake', 'horned_viper', 'diamondback', 'sidewinder', 'trilobite',
  'harvestman', 'scorpion', 'black_and_gold_garden_spider', 'barn_spider', 'garden_spider', 'black_widow', 'tarantula', 'wolf_spider', 'tick', 'centipede',
  'black_grouse', 'ptarmigan', 'ruffed_grouse', 'prairie_chicken', 'peacock', 'quail', 'partridge', 'African_grey', 'macaw', 'sulphur-crested_cockatoo',
  'lorikeet', 'coucal', 'bee_eater', 'hornbill', 'hummingbird', 'jacamar', 'toucan', 'drake', 'red-breasted_merganser', 'goose',
  'black_swan', 'tusker', 'echidna', 'platypus', 'wallaby', 'koala', 'wombat', 'jellyfish', 'sea_anemone', 'brain_coral',
  'flatworm', 'nematode', 'conch', 'snail', 'slug', 'sea_slug', 'chiton', 'chambered_nautilus', 'Dungeness_crab', 'rock_crab',
  'fiddler_crab', 'king_crab', 'American_lobster', 'spiny_lobster', 'crayfish', 'hermit_crab', 'isopod', 'white_stork', 'black_stork', 'spoonbill',
  'flamingo', 'little_blue_heron', 'American_egret', 'bittern', 'crane_(bird)', 'limpkin', 'European_gallinule', 'American_coot', 'bustard', 'ruddy_turnstone',
  'red-backed_sandpiper', 'redshank', 'dowitcher', 'oystercatcher', 'pelican', 'king_penguin', 'albatross', 'grey_whale', 'killer_whale', 'dugong',
  'sea_lion', 'Chihuahua', 'Japanese_spaniel', 'Maltese_dog', 'Pekinese', 'Shih-Tzu', 'Blenheim_spaniel', 'papillon', 'toy_terrier', 'Rhodesian_ridgeback',
  'Afghan_hound', 'basset', 'beagle', 'bloodhound', 'bluetick', 'black-and-tan_coonhound', 'Walker_hound', 'English_foxhound', 'redbone', 'borzoi',
  'Irish_wolfhound', 'Italian_greyhound', 'whippet', 'Ibizan_hound', 'Norwegian_elkhound', 'otterhound', 'Saluki', 'Scottish_deerhound', 'Weimaraner', 'Staffordshire_bullterrier',
  'American_Staffordshire_terrier', 'Bedlington_terrier', 'Border_terrier', 'Kerry_blue_terrier', 'Irish_terrier', 'Norfolk_terrier', 'Norwich_terrier', 'Yorkshire_terrier', 'wire-haired_fox_terrier', 'Lakeland_terrier',
  'Sealyham_terrier', 'Airedale', 'cairn', 'Australian_terrier', 'Dandie_Dinmont', 'Boston_bull', 'miniature_schnauzer', 'giant_schnauzer', 'standard_schnauzer', 'Scotch_terrier',
  'Tibetan_terrier', 'silky_terrier', 'soft-coated_wheaten_terrier', 'West_Highland_white_terrier', 'Lhasa', 'flat-coated_retriever', 'curly-coated_retriever', 'golden_retriever', 'Labrador_retriever', 'Chesapeake_Bay_retriever',
  'German_short-haired_pointer', 'vizsla', 'English_setter', 'Irish_setter', 'Gordon_setter', 'Brittany_spaniel', 'clumber', 'English_springer', 'Welsh_springer_spaniel', 'cocker_spaniel',
  'Sussex_spaniel', 'Irish_water_spaniel', 'kuvasz', 'schipperke', 'groenendael', 'malinois', 'briard', 'kelpie', 'komondor', 'Old_English_sheepdog',
  'Shetland_sheepdog', 'collie', 'Border_collie', 'Bouvier_des_Flandres', 'Rottweiler', 'German_shepherd', 'Doberman', 'miniature_pinscher', 'Greater_Swiss_Mountain_dog', 'Bernese_mountain_dog',
  'Appenzeller', 'EntleBucher', 'boxer', 'bull_mastiff', 'Tibetan_mastiff', 'French_bulldog', 'Great_Dane', 'Saint_Bernard', 'Eskimo_dog', 'malamute',
  'Siberian_husky', 'dalmatian', 'affenpinscher', 'basenji', 'pug', 'Leonberg', 'Newfoundland', 'Great_Pyrenees', 'Samoyed', 'Pomeranian',
  'chow', 'keeshond', 'Brabancon_griffon', 'Pembroke', 'Cardigan', 'toy_poodle', 'miniature_poodle', 'standard_poodle', 'Mexican_hairless', 'timber_wolf',
  'white_wolf', 'red_wolf', 'coyote', 'dingo', 'dhole', 'African_hunting_dog', 'hyena', 'red_fox', 'kit_fox', 'Arctic_fox',
  'grey_fox', 'tabby', 'tiger_cat', 'Persian_cat', 'Siamese_cat', 'Egyptian_cat', 'cougar', 'lynx', 'leopard', 'snow_leopard',
  'jaguar', 'lion', 'tiger', 'cheetah', 'brown_bear', 'American_black_bear', 'ice_bear', 'sloth_bear', 'mongoose', 'meerkat',
  'tiger_beetle', 'ladybug', 'ground_beetle', 'long-horned_beetle', 'leaf_beetle', 'dung_beetle', 'rhinoceros_beetle', 'weevil', 'fly', 'bee',
  'ant', 'grasshopper', 'cricket', 'walking_stick', 'cockroach', 'mantis', 'cicada', 'leafhopper', 'lacewing', 'dragonfly',
  'damselfly', 'admiral', 'ringlet', 'monarch', 'cabbage_butterfly', 'sulphur_butterfly', 'lycaenid', 'starfish', 'sea_urchin', 'sea_cucumber',
  'wood_rabbit', 'hare', 'Angora', 'hamster', 'porcupine', 'fox_squirrel', 'marmot', 'beaver', 'guinea_pig', 'sorrel',
  'zebra', 'hog', 'wild_boar', 'warthog', 'hippopotamus', 'ox', 'water_buffalo', 'bison', 'ram', 'bighorn',
  'ibex', 'hartebeest', 'impala', 'gazelle', 'Arabian_camel', 'llama', 'weasel', 'mink', 'polecat', 'black-footed_ferret',
  'otter', 'skunk', 'badger', 'armadillo', 'three-toed_sloth', 'orangutan', 'gorilla', 'chimpanzee', 'gibbon', 'siamang',
  'guenon', 'patas', 'baboon', 'macaque', 'langur', 'colobus', 'proboscis_monkey', 'marmoset', 'capuchin', 'howler_monkey',
  'titi', 'spider_monkey', 'squirrel_monkey', 'Madagascar_cat', 'indri', 'Indian_elephant', 'African_elephant', 'lesser_panda', 'giant_panda', 'barracouta',
  'eel', 'coho', 'rock_beauty', 'anemone_fish', 'sturgeon', 'gar', 'lionfish', 'puffer', 'abacus', 'abaya',
  'academic_gown', 'accordion', 'acoustic_guitar', 'aircraft_carrier', 'airliner', 'airship', 'altar', 'ambulance', 'amphibian', 'analog_clock',
  'apiary', 'apron', 'ashcan', 'assault_rifle', 'backpack', 'bakery', 'balance_beam', 'balloon', 'ballpoint', 'Band_Aid',
  'banjo', 'bannister', 'barbell', 'barber_chair', 'barbershop', 'barn', 'barometer', 'barrel', 'barrow', 'baseball',
  'basketball', 'bassinet', 'bassoon', 'bathing_cap', 'bath_towel', 'bathtub', 'beach_wagon', 'beacon', 'beaker', 'bearskin',
  'beer_bottle', 'beer_glass', 'bell_cote', 'bib', 'bicycle-built-for-two', 'bikini', 'binder', 'binoculars', 'birdhouse', 'boathouse',
  'bobsled', 'bolo_tie', 'bonnet', 'bookcase', 'bookshop', 'bottlecap', 'bow', 'bow_tie', 'brass', 'brassiere',
  'breakwater', 'breastplate', 'broom', 'bucket', 'buckle', 'bulletproof_vest', 'bullet_train', 'butcher_shop', 'cab', 'caldron',
  'candle', 'cannon', 'canoe', 'can_opener', 'cardigan', 'car_mirror', 'carousel', "carpenter's_kit", 'carton', 'car_wheel',
  'cash_machine', 'cassette', 'cassette_player', 'castle', 'catamaran', 'CD_player', 'cello', 'cellular_telephone', 'chain', 'chainlink_fence',
  'chain_mail', 'chain_saw', 'chest', 'chiffonier', 'chime', 'china_cabinet', 'Christmas_stocking', 'church', 'cinema', 'cleaver',
  'cliff_dwelling', 'cloak', 'clog', 'cocktail_shaker', 'coffee_mug', 'coffeepot', 'coil', 'combination_lock', 'computer_keyboard', 'confectionery',
  'container_ship', 'convertible', 'corkscrew', 'cornet', 'cowboy_boot', 'cowboy_hat', 'cradle', 'crane_(machine)', 'crash_helmet', 'crate',
  'crib', 'Crock_Pot', 'croquet_ball', 'crutch', 'cuirass', 'dam', 'desk', 'desktop_computer', 'dial_telephone', 'diaper',
  'digital_clock', 'digital_watch', 'dining_table', 'dishrag', 'dishwasher', 'disk_brake', 'dock', 'dogsled', 'dome', 'doormat',
  'drilling_platform', 'drum', 'drumstick', 'dumbbell', 'Dutch_oven', 'electric_fan', 'electric_guitar', 'electric_locomotive', 'entertainment_center', 'envelope',
  'espresso_maker', 'face_powder', 'feather_boa', 'file', 'fireboat', 'fire_engine', 'fire_screen', 'flagpole', 'flute', 'folding_chair',
  'football_helmet', 'forklift', 'fountain', 'fountain_pen', 'four-poster', 'freight_car', 'French_horn', 'frying_pan', 'fur_coat', 'garbage_truck',
  'gasmask', 'gas_pump', 'goblet', 'go-kart', 'golf_ball', 'golfcart', 'gondola', 'gong', 'gown', 'grand_piano',
  'greenhouse', 'grille', 'grocery_store', 'guillotine', 'hair_slide', 'hair_spray', 'half_track', 'hammer', 'hamper', 'hand_blower',
  'hand-held_computer', 'handkerchief', 'hard_disc', 'harmonica', 'harp', 'harvester', 'hatchet', 'holster', 'home_theater', 'honeycomb',
  'hook', 'hoopskirt', 'horizontal_bar', 'horse_cart', 'hourglass', 'iPod', 'iron', "jack-o'-lantern", 'jean', 'jeep',
  'jersey', 'jigsaw_puzzle', 'jinrikisha', 'joystick', 'kimono', 'knee_pad', 'knot', 'lab_coat', 'ladle', 'lampshade',
  'laptop', 'lawn_mower', 'lens_cap', 'letter_opener', 'library', 'lifeboat', 'lighter', 'limousine', 'liner', 'lipstick',
  'Loafer', 'lotion', 'loudspeaker', 'loupe', 'lumbermill', 'magnetic_compass', 'mailbag', 'mailbox', 'maillot_(tights)', 'maillot_(tank_suit)',
  'manhole_cover', 'maraca', 'marimba', 'mask', 'matchstick', 'maypole', 'maze', 'measuring_cup', 'medicine_chest', 'megalith',
  'microphone', 'microwave', 'military_uniform', 'milk_can', 'minibus', 'miniskirt', 'minivan', 'missile', 'mitten', 'mixing_bowl',
  'mobile_home', 'Model_T', 'modem', 'monastery', 'monitor', 'moped', 'mortar', 'mortarboard', 'mosque', 'mosquito_net',
  'motor_scooter', 'mountain_bike', 'mountain_tent', 'mouse', 'mousetrap', 'moving_van', 'muzzle', 'nail', 'neck_brace', 'necklace',
  'nipple', 'notebook', 'obelisk', 'oboe', 'ocarina', 'odometer', 'oil_filter', 'organ', 'oscilloscope', 'overskirt',
  'oxcart', 'oxygen_mask', 'packet', 'paddle', 'paddlewheel', 'padlock', 'paintbrush', 'pajama', 'palace', 'panpipe',
  'paper_towel', 'parachute', 'parallel_bars', 'park_bench', 'parking_meter', 'passenger_car', 'patio', 'pay-phone', 'pedestal', 'pencil_box',
  'pencil_sharpener', 'perfume', 'Petri_dish', 'photocopier', 'pick', 'pickelhaube', 'picket_fence', 'pickup', 'pier', 'piggy_bank',
  'pill_bottle', 'pillow', 'ping-pong_ball', 'pinwheel', 'pirate', 'pitcher', 'plane', 'planetarium', 'plastic_bag', 'plate_rack',
  'plow', 'plunger', 'Polaroid_camera', 'pole', 'police_van', 'poncho', 'pool_table', 'pop_bottle', 'pot', "potter's_wheel",
  'power_drill', 'prayer_rug', 'printer', 'prison', 'projectile', 'projector', 'puck', 'punching_bag', 'purse', 'quill',
  'quilt', 'racer', 'racket', 'radiator', 'radio', 'radio_telescope', 'rain_barrel', 'recreational_vehicle', 'reel', 'reflex_camera',
  'refrigerator', 'remote_control', 'restaurant', 'revolver', 'rifle', 'rocking_chair', 'rotisserie', 'rubber_eraser', 'rugby_ball', 'rule',
  'running_shoe', 'safe', 'safety_pin', 'saltshaker', 'sandal', 'sarong', 'sax', 'scabbard', 'scale', 'school_bus',
  'schooner', 'scoreboard', 'screen', 'screw', 'screwdriver', 'seat_belt', 'sewing_machine', 'shield', 'shoe_shop', 'shoji',
  'shopping_basket', 'shopping_cart', 'shovel', 'shower_cap', 'shower_curtain', 'ski', 'ski_mask', 'sleeping_bag', 'slide_rule', 'sliding_door',
  'slot', 'snorkel', 'snowmobile', 'snowplow', 'soap_dispenser', 'soccer_ball', 'sock', 'solar_dish', 'sombrero', 'soup_bowl',
  'space_bar', 'space_heater', 'space_shuttle', 'spatula', 'speedboat', 'spider_web', 'spindle', 'sports_car', 'spotlight', 'stage',
  'steam_locomotive', 'steel_arch_bridge', 'steel_drum', 'stethoscope', 'stole', 'stone_wall', 'stopwatch', 'stove', 'strainer', 'streetcar',
  'stretcher', 'studio_couch', 'stupa', 'submarine', 'suit', 'sundial', 'sunglass', 'sunglasses', 'sunscreen', 'suspension_bridge',
  'swab', 'sweatshirt', 'swimming_trunks', 'swing', 'switch', 'syringe', 'table_lamp', 'tank', 'tape_player', 'teapot',
  'teddy', 'television', 'tennis_ball', 'thatch', 'theater_curtain', 'thimble', 'thresher', 'throne', 'tile_roof', 'toaster',
  'tobacco_shop', 'toilet_seat', 'torch', 'totem_pole', 'tow_truck', 'toyshop', 'tractor', 'trailer_truck', 'tray', 'trench_coat',
  'tricycle', 'trimaran', 'tripod', 'triumphal_arch', 'trolleybus', 'trombone', 'tub', 'turnstile', 'typewriter_keyboard', 'umbrella',
  'unicycle', 'upright', 'vacuum', 'vase', 'vault', 'velvet', 'vending_machine', 'vestment', 'viaduct', 'violin',
  'volleyball', 'waffle_iron', 'wall_clock', 'wallet', 'wardrobe', 'warplane', 'washbasin', 'washer', 'water_bottle', 'water_jug',
  'water_tower', 'whiskey_jug', 'whistle', 'wig', 'window_screen', 'window_shade', 'Windsor_tie', 'wine_bottle', 'wing', 'wok',
  'wooden_spoon', 'wool', 'worm_fence', 'wreck', 'yawl', 'yurt', 'web_site', 'comic_book', 'crossword_puzzle', 'street_sign',
  'traffic_light', 'book_jacket', 'menu', 'plate', 'guacamole', 'consomme', 'hot_pot', 'trifle', 'ice_cream', 'ice_lolly',
  'French_loaf', 'bagel', 'pretzel', 'cheeseburger', 'hotdog', 'mashed_potato', 'head_cabbage', 'broccoli', 'cauliflower', 'zucchini',
  'spaghetti_squash', 'acorn_squash', 'butternut_squash', 'cucumber', 'artichoke', 'bell_pepper', 'cardoon', 'mushroom', 'Granny_Smith', 'strawberry',
  'orange', 'lemon', 'fig', 'pineapple', 'banana', 'jackfruit', 'custard_apple', 'pomegranate', 'hay', 'carbonara',
  'chocolate_sauce', 'dough', 'meat_loaf', 'pizza', 'potpie', 'burrito', 'red_wine', 'espresso', 'cup', 'eggnog',
  'alp', 'bubble', 'cliff', 'coral_reef', 'geyser', 'lakeside', 'promontory', 'sandbar', 'seashore', 'valley',
  'volcano', 'ballplayer', 'groom', 'scuba_diver', 'rapeseed', 'daisy', "yellow_lady's_slipper", 'corn', 'acorn', 'hip',
  'buckeye', 'coral_fungus', 'agaric', 'gyromitra', 'stinkhorn', 'earthstar', 'hen-of-the-woods', 'bolete', 'ear', 'toilet_tissue'
)

#' Get Base URL for Pre-Trained Neural Models
#'
#' Returns the base URL for downloading pre-trained ONNX models from GitHub Releases.
#' Can be customized globally via \code{options(pliman.models_base_url = "...")}.
#'
#' @return Character URL string ending in a slash.
#' @export
pliman_models_base_url <- function() {
  getOption("pliman.models_base_url", "https://github.com/NEPEM-UFSC/models/releases/download/v1.0.0/")
}

# Internal helper for robust downloads with SSL error recovery
.download_url_robust <- function(url, dest_file, min_size = 1000) {
  dl_success <- FALSE
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    tryCatch({
      curl::curl_download(url, destfile = dest_file, quiet = FALSE, mode = "wb", handle = h)
      dl_success <- file.exists(dest_file) && file.info(dest_file)$size >= min_size
    }, error = function(e) {
      if (grepl("certificate|schannel|SEC_E|SSL", e$message, ignore.case = TRUE)) {
        cli::cli_alert_warning("SSL handshake failed ({e$message}). Retrying with SSL verification bypass...")
        h_insecure <- curl::new_handle()
        curl::handle_setopt(h_insecure, ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)
        tryCatch({
          curl::curl_download(url, destfile = dest_file, quiet = FALSE, mode = "wb", handle = h_insecure)
          dl_success <<- file.exists(dest_file) && file.info(dest_file)$size >= min_size
        }, error = function(e2) NULL)
      }
    })
  }

  if (!dl_success) {
    tryCatch({
      utils::download.file(url, destfile = dest_file, mode = "wb", quiet = FALSE)
      dl_success <- file.exists(dest_file) && file.info(dest_file)$size >= min_size
    }, error = function(e) NULL)
  }

  if (!isTRUE(dl_success) || !file.exists(dest_file) || file.info(dest_file)$size < min_size) {
    if (file.exists(dest_file)) unlink(dest_file)
    return(FALSE)
  }
  return(TRUE)
}

.pliman_assets_cache <- new.env(parent = emptyenv())

# Internal helper to query GitHub release assets with local disk & memory cache
.fetch_github_release_assets <- function(refresh = FALSE, dir = pliman_model_dir()) {
  cache_file <- file.path(dir, ".github_assets_cache.rds")

  # 1. In-memory session cache (valid for 12 hours)
  if (!isTRUE(refresh) && exists("assets", envir = .pliman_assets_cache)) {
    cache_time <- get0("time", envir = .pliman_assets_cache, ifnotfound = 0)
    if (as.numeric(difftime(Sys.time(), cache_time, units = "hours")) < 12) {
      return(get("assets", envir = .pliman_assets_cache))
    }
  }

  # 2. Disk cache (valid for 24 hours)
  if (!isTRUE(refresh) && file.exists(cache_file)) {
    f_info <- file.info(cache_file)
    age_hours <- as.numeric(difftime(Sys.time(), f_info$mtime, units = "hours"))
    if (!is.na(age_hours) && age_hours < 24) {
      cached_data <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.data.frame(cached_data) && nrow(cached_data) > 0) {
        assign("assets", cached_data, envir = .pliman_assets_cache)
        assign("time", Sys.time(), envir = .pliman_assets_cache)
        return(cached_data)
      }
    }
  }

  # 3. Query GitHub API
  api_url <- "https://api.github.com/repos/NEPEM-UFSC/models/releases/tags/v1.0.0"
  assets_df <- NULL

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    raw_res <- tryCatch({
      if (requireNamespace("curl", quietly = TRUE)) {
        h <- curl::new_handle("useragent" = "pliman-R-package")
        con <- curl::curl(api_url, handle = h)
        on.exit(tryCatch(close(con), error = function(e) NULL), add = TRUE)
        jsonlite::fromJSON(con)
      } else {
        jsonlite::fromJSON(api_url)
      }
    }, error = function(e) NULL)

    if (is.list(raw_res) && !is.null(raw_res$assets) && is.data.frame(raw_res$assets)) {
      raw_assets <- raw_res$assets
      valid_rows <- grepl("\\.(onnx|pt)$", raw_assets$name, ignore.case = TRUE)
      if (any(valid_rows)) {
        filtered <- raw_assets[valid_rows, c("name", "size", "browser_download_url")]
        filtered$size_mb <- round(as.numeric(filtered$size) / (1024 * 1024), 1)
        filtered$model <- tools::file_path_sans_ext(filtered$name)
        assets_df <- filtered
      }
    }
  }

  if (is.data.frame(assets_df) && nrow(assets_df) > 0) {
    assign("assets", assets_df, envir = .pliman_assets_cache)
    assign("time", Sys.time(), envir = .pliman_assets_cache)
    tryCatch({
      if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(assets_df, cache_file)
    }, error = function(e) invisible(NULL))
    return(assets_df)
  }

  # Fallback to existing disk cache if API failed
  if (file.exists(cache_file)) {
    cached_data <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.data.frame(cached_data) && nrow(cached_data) > 0) {
      return(cached_data)
    }
  }

  data.frame(
    name = character(0),
    size = numeric(0),
    browser_download_url = character(0),
    size_mb = numeric(0),
    model = character(0),
    stringsAsFactors = FALSE
  )
}

#' List Available Pre-Trained Neural Models
#'
#' Lists all supported background removal, object detection, segmentation, and classification models,
#' their disk sizes, input resolutions, descriptions, local download status, and models dynamically
#' discovered from GitHub release assets.
#'
#' @param dir Directory where models are stored (default: `pliman_model_dir()`).
#' @param refresh Logical. If `TRUE`, forces re-fetching the asset catalog from GitHub releases (default `FALSE`).
#' @return A data frame containing model metadata and local availability.
#' @export
#' @examples
#' \dontrun{
#' pliman_available_models()
#' pliman_available_models(refresh = TRUE)
#' }
pliman_available_models <- function(dir = pliman_model_dir(), refresh = FALSE) {
  base_url <- pliman_models_base_url()
  hf_yolo_url <- "https://huggingface.co/zwh20081/yolo26-onnx/resolve/main/"
  models_info <- list(
    list(
      name = "u2netp",
      description = "U2-Net Portable: Ultra-lightweight (~4.6 MB) and fast for CPU",
      size_mb = 4.6,
      input_size = 320,
      filename = "u2netp.onnx",
      url = paste0(base_url, "u2netp.onnx")
    ),
    list(
      name = "birefnet-lite",
      description = "BiRefNet Lite: Bilateral Reference Network for High-Resolution Dichotomous Segmentation (1024x1024)",
      size_mb = 213.6,
      input_size = 1024,
      filename = "birefnet-lite.onnx",
      url = paste0(base_url, "birefnet-lite.onnx")
    ),
    list(
      name = "isnet-general-use",
      description = "IS-Net: High-resolution (1024x1024) boundary matting (best for leaf contours & lesions)",
      size_mb = 178.6,
      input_size = 1024,
      filename = "isnet-general-use.onnx",
      url = paste0(base_url, "isnet-general-use.onnx")
    ),
    list(
      name = "rmbg-1.4",
      description = "BRIA RMBG 1.4: State-of-the-Art background removal (1024x1024, alias 'rmbg')",
      size_mb = 176.1,
      input_size = 1024,
      filename = "rmbg-1.4.onnx",
      url = paste0(base_url, "rmbg-1.4.onnx")
    ),
    list(
      name = "rmbg-2.0",
      description = "BRIA RMBG 2.0: Next-Gen BiRefNet-based background removal (1024x1024)",
      size_mb = 976.9,
      input_size = 1024,
      filename = "rmbg-2.0.onnx",
      url = paste0(base_url, "rmbg-2.0.onnx")
    ),
    list(
      name = "silueta",
      description = "Silueta: Compact model optimized for high accuracy and low latency (320x320)",
      size_mb = 44.2,
      input_size = 320,
      filename = "silueta.onnx",
      url = paste0(base_url, "silueta.onnx")
    ),
    list(
      name = "u2net",
      description = "U2-Net: Full salient object detection model (320x320)",
      size_mb = 176.3,
      input_size = 320,
      filename = "u2net.onnx",
      url = paste0(base_url, "u2net.onnx")
    ),
    list(
      name = "withoutbg",
      description = "withoutBG Open Weights: DepthAnythingV2 + ConvNeXt-fused U-Net matting (448x448)",
      size_mb = 433.4,
      input_size = 448,
      filename = "withoutbg.onnx",
      url = paste0(base_url, "withoutbg.onnx")
    ),
    list(
      name = "sam2.1",
      description = "Segment Anything Model 2.1 (Hiera-Tiny, Meta AI): Foundation image segmentation (1024x1024)",
      size_mb = 120.2,
      input_size = 1024,
      filename = "sam2.1.encoder.onnx",
      url = paste0(base_url, "sam2.1.encoder.onnx")
    ),
    list(
      name = "sam3.1",
      description = "Segment Anything Model 3.1 (Meta AI): Concept-driven foundation segmentation (1024x1024)",
      size_mb = 868.1,
      input_size = 1024,
      filename = "sam3.1.onnx",
      url = paste0(base_url, "sam3.1.onnx")
    ),
    list(
      name = "ben2",
      description = "Boundary-aware Extraction Network (BEN2 Base, Prama LLC): High-accuracy crisp boundary extraction (1024x1024)",
      size_mb = 212.6,
      input_size = 1024,
      filename = "ben2.onnx",
      url = paste0(base_url, "ben2.onnx")
    ),
    list(
      name = "grounded-sam",
      description = "Grounded-SAM: Grounding DINO Tiny + SAM 2.1 zero-shot text-prompted instance segmentation (800x800 & 1024x1024)",
      size_mb = 314.8,
      input_size = 800,
      filename = "groundingdino-tiny.onnx",
      url = paste0(base_url, "groundingdino-tiny.onnx")
    ),
    list(
      name = "persam",
      description = "PerSAM: One-shot visual exemplar instance segmentation and counting via SAM 2.1 embeddings (1024x1024)",
      size_mb = 120.2,
      input_size = 1024,
      filename = "sam2.1.encoder.onnx",
      url = paste0(base_url, "sam2.1.encoder.onnx")
    ),
    list(
      name = "depth-anything-v2",
      description = "Depth Anything V2 (Small): Monocular relative 3D depth estimation (518x518)",
      size_mb = 95.5,
      input_size = 518,
      filename = "depth-anything-v2-small.onnx",
      url = paste0(base_url, "depth-anything-v2-small.onnx")
    ),
    list(
      name = "dinov2",
      description = "DINOv2 (ViT-S/14, Meta AI): Self-supervised dense foundation vision transformer (518x518)",
      size_mb = 85.2,
      input_size = 518,
      filename = "dinov2-vits14.onnx",
      url = paste0(base_url, "dinov2-vits14.onnx")
    ),
    # YOLO26 Object Detection (Ultralytics)
    list(
      name = "yolo26n",
      description = "YOLO26 Nano Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 9.5,
      input_size = 640,
      filename = "yolo26n.onnx",
      url = paste0(hf_yolo_url, "yolo26n.onnx")
    ),
    list(
      name = "yolo26s",
      description = "YOLO26 Small Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 36.5,
      input_size = 640,
      filename = "yolo26s.onnx",
      url = paste0(hf_yolo_url, "yolo26s.onnx")
    ),
    list(
      name = "yolo26m",
      description = "YOLO26 Medium Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 78.2,
      input_size = 640,
      filename = "yolo26m.onnx",
      url = paste0(hf_yolo_url, "yolo26m.onnx")
    ),
    list(
      name = "yolo26l",
      description = "YOLO26 Large Detection: Real-time general object detection (640x640, 80 COCO classes)",
      size_mb = 95.0,
      input_size = 640,
      filename = "yolo26l.onnx",
      url = paste0(hf_yolo_url, "yolo26l.onnx")
    ),
    list(
      name = "yolo26x",
      description = "YOLO26 Extra-Large Detection: High-accuracy general object detection (640x640, 80 COCO classes)",
      size_mb = 212.9,
      input_size = 640,
      filename = "yolo26x.onnx",
      url = paste0(hf_yolo_url, "yolo26x.onnx")
    ),

    # YOLO26 Instance Segmentation (Ultralytics)
    list(
      name = "yolo26n-seg",
      description = "YOLO26 Nano Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 10.7,
      input_size = 640,
      filename = "yolo26n-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26n-seg.onnx")
    ),
    list(
      name = "yolo26s-seg",
      description = "YOLO26 Small Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 40.0,
      input_size = 640,
      filename = "yolo26s-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26s-seg.onnx")
    ),
    list(
      name = "yolo26m-seg",
      description = "YOLO26 Medium Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 90.2,
      input_size = 640,
      filename = "yolo26m-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26m-seg.onnx")
    ),
    list(
      name = "yolo26l-seg",
      description = "YOLO26 Large Segmentation: Real-time instance segmentation (640x640, 80 COCO classes)",
      size_mb = 107.1,
      input_size = 640,
      filename = "yolo26l-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26l-seg.onnx")
    ),
    list(
      name = "yolo26x-seg",
      description = "YOLO26 Extra-Large Segmentation: High-accuracy instance segmentation (640x640, 80 COCO classes)",
      size_mb = 240.0,
      input_size = 640,
      filename = "yolo26x-seg.onnx",
      url = paste0(hf_yolo_url, "yolo26x-seg.onnx")
    ),

    # YOLO26 Pose Estimation (Ultralytics)
    list(
      name = "yolo26n-pose",
      description = "YOLO26 Nano Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 11.6,
      input_size = 640,
      filename = "yolo26n-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26n-pose.onnx")
    ),
    list(
      name = "yolo26s-pose",
      description = "YOLO26 Small Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 39.9,
      input_size = 640,
      filename = "yolo26s-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26s-pose.onnx")
    ),
    list(
      name = "yolo26m-pose",
      description = "YOLO26 Medium Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 82.6,
      input_size = 640,
      filename = "yolo26m-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26m-pose.onnx")
    ),
    list(
      name = "yolo26l-pose",
      description = "YOLO26 Large Pose: Real-time human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 99.4,
      input_size = 640,
      filename = "yolo26l-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26l-pose.onnx")
    ),
    list(
      name = "yolo26x-pose",
      description = "YOLO26 Extra-Large Pose: High-accuracy human pose estimation (640x640, 17 COCO keypoints)",
      size_mb = 220.0,
      input_size = 640,
      filename = "yolo26x-pose.onnx",
      url = paste0(hf_yolo_url, "yolo26x-pose.onnx")
    ),

    # YOLO26 Image Classification (Ultralytics)
    list(
      name = "yolo26n-cls",
      description = "YOLO26 Nano Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 10.8,
      input_size = 640,
      filename = "yolo26n-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26n-cls.onnx")
    ),
    list(
      name = "yolo26s-cls",
      description = "YOLO26 Small Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 25.7,
      input_size = 640,
      filename = "yolo26s-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26s-cls.onnx")
    ),
    list(
      name = "yolo26m-cls",
      description = "YOLO26 Medium Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 44.4,
      input_size = 640,
      filename = "yolo26m-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26m-cls.onnx")
    ),
    list(
      name = "yolo26l-cls",
      description = "YOLO26 Large Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 53.9,
      input_size = 640,
      filename = "yolo26l-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26l-cls.onnx")
    ),
    list(
      name = "yolo26x-cls",
      description = "YOLO26 Extra-Large Classification: Real-time image classification (640x640, 1000 ImageNet classes)",
      size_mb = 113.2,
      input_size = 640,
      filename = "yolo26x-cls.onnx",
      url = paste0(hf_yolo_url, "yolo26x-cls.onnx")
    ),
    list(
      name = "stardist",
      description = "StarDist (DSB 2018): Star-convex polygon detection for round/overlapping objects & cells",
      size_mb = 34.8,
      input_size = 256,
      filename = "stardist-dsb2018.onnx",
      url = paste0(base_url, "stardist-dsb2018.onnx")
    ),
    list(
      name = "realesrgan-compact",
      description = "Real-ESRGAN Compact: Fast 4x generative super-resolution with edge preservation",
      size_mb = 16.7,
      input_size = 256,
      filename = "realesrgan-compact.onnx",
      url = paste0(base_url, "realesrgan-compact.onnx")
    ),
    list(
      name = "clip-vit-b32",
      description = "CLIP (ViT-B/32, OpenAI): Zero-shot image classification and cross-modal embedding (224x224)",
      size_mb = 78.4,
      input_size = 224,
      filename = "clip-vit-b32-vision.onnx",
      url = paste0(base_url, "vision_model_quantized.onnx")
    ),
    list(
      name = "yolov8s-worldv2",
      description = "YOLO-World Small (v8s): Real-time zero-shot open-vocabulary object detection (640x640)",
      size_mb = 25.9,
      input_size = 640,
      filename = "yolov8s-worldv2.pt",
      url = "https://github.com/ultralytics/assets/releases/download/v8.2.0/yolov8s-worldv2.pt"
    ),
    list(
      name = "yolov8m-worldv2",
      description = "YOLO-World Medium (v8m): Real-time zero-shot open-vocabulary object detection (640x640)",
      size_mb = 54.4,
      input_size = 640,
      filename = "yolov8m-worldv2.pt",
      url = "https://github.com/ultralytics/assets/releases/download/v8.2.0/yolov8m-worldv2.pt"
    ),
    list(
      name = "yolov8l-worldv2",
      description = "YOLO-World Large (v8l): Real-time zero-shot open-vocabulary object detection (640x640)",
      size_mb = 89.9,
      input_size = 640,
      filename = "yolov8l-worldv2.pt",
      url = "https://github.com/ultralytics/assets/releases/download/v8.2.0/yolov8l-worldv2.pt"
    ),
    list(
      name = "yolov8x-worldv2",
      description = "YOLO-World XLarge (v8x): Real-time zero-shot open-vocabulary object detection (640x640)",
      size_mb = 140.5,
      input_size = 640,
      filename = "yolov8x-worldv2.pt",
      url = "https://github.com/ultralytics/assets/releases/download/v8.2.0/yolov8x-worldv2.pt"
    ),
    list(
      name = "yolo_nas_s",
      description = "YOLO-NAS Small: Neural Architecture Search object detection (640x640)",
      size_mb = 46.6,
      input_size = 640,
      filename = "yolo_nas_s.onnx",
      url = "https://github.com/CVHub520/X-AnyLabeling/releases/download/v0.1.0/yolo_nas_s.onnx"
    ),
    list(
      name = "yolo_nas_l",
      description = "YOLO-NAS Large: Neural Architecture Search object detection (640x640)",
      size_mb = 160.4,
      input_size = 640,
      filename = "yolo_nas_l.onnx",
      url = "https://github.com/CVHub520/X-AnyLabeling/releases/download/v0.1.0/yolo_nas_l.onnx"
    ),
    # YOLO26 Oriented Bounding Box (Ultralytics)
    list(
      name = "yolo26n-obb",
      description = "YOLO26 Nano OBB: Oriented Bounding Box detection (dynamic/1024x1024, 15 DOTA classes)",
      size_mb = 11.5,
      input_size = "dynamic",
      filename = "yolo26n-obb.onnx",
      url = paste0(base_url, "yolo26n-obb.onnx")
    ),
    list(
      name = "yolo26s-obb",
      description = "YOLO26 Small OBB: Oriented Bounding Box detection (dynamic/1024x1024, 15 DOTA classes)",
      size_mb = 39.7,
      input_size = "dynamic",
      filename = "yolo26s-obb.onnx",
      url = paste0(base_url, "yolo26s-obb.onnx")
    ),
    list(
      name = "yolo26m-obb",
      description = "YOLO26 Medium OBB: Oriented Bounding Box detection (dynamic/1024x1024, 15 DOTA classes)",
      size_mb = 84.6,
      input_size = "dynamic",
      filename = "yolo26m-obb.onnx",
      url = paste0(base_url, "yolo26m-obb.onnx")
    ),
    list(
      name = "yolo26l-obb",
      description = "YOLO26 Large OBB: Oriented Bounding Box detection (dynamic/1024x1024, 15 DOTA classes)",
      size_mb = 106.1,
      input_size = "dynamic",
      filename = "yolo26l-obb.onnx",
      url = paste0(base_url, "yolo26l-obb.onnx")
    ),
    list(
      name = "yolo26x-obb",
      description = "YOLO26 Extra-Large OBB: Oriented Bounding Box detection (dynamic/1024x1024, 15 DOTA classes)",
      size_mb = 236.5,
      input_size = "dynamic",
      filename = "yolo26x-obb.onnx",
      url = paste0(base_url, "yolo26x-obb.onnx")
    ),
    list(
      name = "ifcnn",
      description = "IFCNN (Image Fusion CNN): General multi-focus and all-in-focus image fusion",
      size_mb = 0.34,
      input_size = "dynamic",
      filename = "ifcnn.onnx",
      url = paste0(base_url, "ifcnn.onnx")
    )
  )

  # Merge dynamic assets discovered from GitHub release
  gh_assets <- .fetch_github_release_assets(refresh = refresh, dir = dir)
  if (is.data.frame(gh_assets) && nrow(gh_assets) > 0) {
    existing_files <- vapply(models_info, function(m) m$filename, character(1))
    existing_names <- vapply(models_info, function(m) m$name, character(1))
    new_assets <- gh_assets[!gh_assets$name %in% existing_files, ]
    if (nrow(new_assets) > 0) {
      for (i in seq_len(nrow(new_assets))) {
        a_name <- new_assets$name[i]
        a_model <- new_assets$model[i]
        a_size <- new_assets$size_mb[i]
        a_url <- new_assets$browser_download_url[i]
        clean_name <- if (a_model %in% existing_names) a_name else a_model

        desc <- if (grepl("[-_]cls", a_name, ignore.case = TRUE)) {
          paste0("YOLO Classification Model: Real-time image classification (", a_name, ")")
        } else if (grepl("[-_]seg", a_name, ignore.case = TRUE)) {
          paste0("YOLO Instance Segmentation Model: Real-time mask segmentation (", a_name, ")")
        } else if (grepl("[-_]pose", a_name, ignore.case = TRUE)) {
          paste0("YOLO Pose Estimation Model: Keypoint & skeleton detection (", a_name, ")")
        } else if (grepl("[-_]obb", a_name, ignore.case = TRUE)) {
          paste0("YOLO Oriented Bounding Box Model: Oriented object detection (", a_name, ")")
        } else if (grepl("world", a_name, ignore.case = TRUE)) {
          paste0("YOLO-World Open-Vocabulary Model: Zero-shot object detection (", a_name, ")")
        } else if (grepl("^yolo", a_name, ignore.case = TRUE)) {
          paste0("YOLO Object Detection Model: General object detection (", a_name, ")")
        } else {
          paste0("Pre-trained neural model asset (", a_name, ")")
        }

        inp_size <- if (grepl("obb|ifcnn", a_name, ignore.case = TRUE)) "dynamic" else 640

        models_info[[length(models_info) + 1]] <- list(
          name = clean_name,
          description = desc,
          size_mb = a_size,
          input_size = inp_size,
          filename = a_name,
          url = a_url
        )
      }
    }
  }

  res <- do.call(rbind, lapply(models_info, function(m) {
    file_path <- file.path(dir, m$filename)
    downloaded <- file.exists(file_path) && file.info(file_path)$size > 1000
    data.frame(
      model = m$name,
      size_mb = m$size_mb,
      input_size = if (is.numeric(m$input_size)) paste0(m$input_size, "x", m$input_size) else as.character(m$input_size),
      downloaded = downloaded,
      description = m$description,
      filename = m$filename,
      url = m$url,
      path = ifelse(downloaded, file_path, NA_character_),
      stringsAsFactors = FALSE
    )
  }))

  return(res)
}

#' Download Pre-Trained Neural Models
#'
#' Downloads one or more pre-trained neural models into the user's model directory.
#' Supports downloading by model name, model asset name, dynamic GitHub release asset,
#' or directly from a custom URL (e.g. `https://.../model.onnx`).
#'
#' @param model Character vector naming the model(s) or direct URL(s) to download.
#'   Can be a known model name (e.g. `"yolo26n"`, `"rmbg-2.0"`), an asset filename
#'   (e.g. `"yolo26x-cls.onnx"`, `"yolo26s-obb.pt"`), `"all"`, or a direct HTTP/HTTPS URL.
#' @param dir Directory to save the model files (default: `pliman_model_dir()`).
#' @param force Logical. If `TRUE`, forces re-download even if the file exists (default `FALSE`).
#' @param refresh Logical. If `TRUE`, forces re-fetching the GitHub release asset catalog (default `FALSE`).
#' @return A character vector with the absolute path(s) to the downloaded model file(s).
#' @export
#' @examples
#' \dontrun{
#'   # Download a single model
#'   pliman_download_model("u2netp")
#'
#'   # Download directly from a custom URL
#'   pliman_download_model("https://example.com/custom_yolo.onnx")
#'
#'   # Download all models at once
#'   pliman_download_model("all")
#' }
pliman_download_model <- function(model = "u2netp",
                                  dir = pliman_model_dir(),
                                  force = FALSE,
                                  refresh = FALSE) {
  dir <- pliman_model_dir(dir)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(600, old_timeout))

  download_single <- function(m) {
    # CASE 1: Direct URL download
    if (grepl("^https?://", m)) {
      clean_url <- sub("[?#].*$", "", m)
      fname <- basename(clean_url)
      if (!nzchar(fname) || fname %in% c("/", ".")) {
        fname <- paste0("model_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".onnx")
      }
      dest_file <- file.path(dir, fname)
      if (!isTRUE(force) && file.exists(dest_file) && file.info(dest_file)$size > 1000) {
        return(normalizePath(dest_file, winslash = "/"))
      }
      cli::cli_alert_info("Downloading model directly from URL {.url {m}} to {.path {dest_file}}...")
      ok <- .download_url_robust(m, dest_file, min_size = 1000)
      if (!ok || !file.exists(dest_file) || file.info(dest_file)$size < 1000) {
        if (file.exists(dest_file)) unlink(dest_file)
        stop("Failed to download model from URL: ", m)
      }
      cli::cli_alert_success("Model {.file {fname}} successfully downloaded!")
      return(normalizePath(dest_file, winslash = "/"))
    }

    # CASE 2: Local file already exists
    if (file.exists(m)) {
      return(normalizePath(m, winslash = "/"))
    }
    if (file.exists(file.path(dir, m))) {
      return(normalizePath(file.path(dir, m), winslash = "/"))
    }
    if (file.exists(file.path(dir, paste0(m, ".onnx")))) {
      return(normalizePath(file.path(dir, paste0(m, ".onnx")), winslash = "/"))
    }
    if (file.exists(file.path(dir, paste0(m, ".pt")))) {
      return(normalizePath(file.path(dir, paste0(m, ".pt")), winslash = "/"))
    }

    # CASE 3: Model lookup in catalog
    resolved_m <- .resolve_model_name(m)
    models_df <- pliman_available_models(dir = dir, refresh = refresh)

    match_idx <- which(
      models_df$model == resolved_m |
      models_df$filename == resolved_m |
      models_df$model == m |
      models_df$filename == m
    )

    if (length(match_idx) == 0 && !isTRUE(refresh)) {
      models_df <- pliman_available_models(dir = dir, refresh = TRUE)
      match_idx <- which(
        models_df$model == resolved_m |
        models_df$filename == resolved_m |
        models_df$model == m |
        models_df$filename == m
      )
    }

    base_url <- pliman_models_base_url()

    # CASE 4: Direct probe on GitHub releases
    if (length(match_idx) == 0) {
      candidate_fnames <- unique(c(
        m,
        paste0(m, ".onnx"),
        paste0(m, ".pt"),
        resolved_m,
        paste0(resolved_m, ".onnx"),
        paste0(resolved_m, ".pt")
      ))
      candidate_fnames <- candidate_fnames[grepl("\\.(onnx|pt)$", candidate_fnames)]

      downloaded_probe <- NULL
      for (cand in candidate_fnames) {
        probe_url <- paste0(base_url, cand)
        probe_dest <- file.path(dir, cand)
        cli::cli_alert_info("Probing GitHub release asset {.val {cand}} at {.url {probe_url}}...")
        ok_probe <- .download_url_robust(probe_url, probe_dest, min_size = 1000)
        if (ok_probe && file.exists(probe_dest) && file.info(probe_dest)$size > 1000) {
          cli::cli_alert_success("Model asset {.file {cand}} found and downloaded successfully!")
          downloaded_probe <- probe_dest
          break
        } else {
          if (file.exists(probe_dest)) unlink(probe_dest)
        }
      }

      if (!is.null(downloaded_probe)) {
        return(normalizePath(downloaded_probe, winslash = "/"))
      }

      valid_options <- unique(c(models_df$model, models_df$filename))
      stop("Unknown model(s): '", m, "'.\n",
           "Could not find model locally, in catalog, or as an asset on GitHub releases (NEPEM-UFSC/models).\n",
           "Valid available options are: 'all', ",
           paste(paste0("'", utils::head(valid_options, 35), "'"), collapse = ", "),
           if (length(valid_options) > 35) paste0(", ... (", length(valid_options) - 35, " more)") else "")
    }

    row <- models_df[match_idx[1], ]
    actual_model <- row$model

    if (actual_model %in% c("sam2.1", "persam", "per-sam")) {
      enc_file <- file.path(dir, "sam2.1.encoder.onnx")
      dec_file <- file.path(dir, "sam2.1.decoder.onnx")
      if (!isTRUE(force) && file.exists(enc_file) && file.info(enc_file)$size > 10000000 &&
          file.exists(dec_file) && file.info(dec_file)$size > 1000000) {
        return(normalizePath(enc_file, winslash = "/"))
      }
      cli::cli_alert_info("Downloading SAM 2.1 encoder and decoder ({row$size_mb} MB) to {.path {dir}}...")
      enc_url <- paste0(base_url, "sam2.1.encoder.onnx")
      dec_url <- paste0(base_url, "sam2.1.decoder.onnx")
      ok_enc <- .download_url_robust(enc_url, enc_file, min_size = 10000000)
      ok_dec <- .download_url_robust(dec_url, dec_file, min_size = 1000000)
      if (!ok_enc || !ok_dec) {
        stop("Failed to download SAM 2.1 model. Please check your internet connection.")
      }
      cli::cli_alert_success("Model {.val {actual_model}} successfully downloaded!")
      return(normalizePath(enc_file, winslash = "/"))
    }

    if (actual_model %in% c("grounded-sam", "grounding-dino")) {
      model_file <- file.path(dir, "groundingdino-tiny.onnx")
      vocab_file <- file.path(dir, "vocab.txt")
      if (!isTRUE(force) && file.exists(model_file) && file.info(model_file)$size > 100000000 &&
          file.exists(vocab_file) && file.info(vocab_file)$size > 50000) {
        return(normalizePath(model_file, winslash = "/"))
      }
      cli::cli_alert_info("Downloading Grounding DINO Tiny and BERT vocab ({row$size_mb} MB) to {.path {dir}}...")
      dino_url <- paste0(base_url, "groundingdino-tiny.onnx")
      vocab_url <- paste0(base_url, "vocab.txt")
      ok_dino <- .download_url_robust(dino_url, model_file, min_size = 100000000)
      ok_vocab <- .download_url_robust(vocab_url, vocab_file, min_size = 50000)
      if (!ok_dino || !ok_vocab) {
        stop("Failed to download Grounding DINO model/vocab. Please check your internet connection.")
      }
      cli::cli_alert_success("Grounding DINO model and vocab successfully downloaded!")
      return(normalizePath(model_file, winslash = "/"))
    }

    if (actual_model %in% c("clip", "clip-vit-b32", "clip-vit-base-patch32")) {
      vis_file <- file.path(dir, "clip-vit-b32-vision.onnx")
      txt_file <- file.path(dir, "clip-vit-b32-text.onnx")
      vocab_file <- file.path(dir, "clip_vocab.json")
      merges_file <- file.path(dir, "clip_merges.txt")
      if (!isTRUE(force) && file.exists(vis_file) && file.info(vis_file)$size > 30000000 &&
          file.exists(txt_file) && file.info(txt_file)$size > 30000000 &&
          file.exists(vocab_file) && file.info(vocab_file)$size > 500000 &&
          file.exists(merges_file) && file.info(merges_file)$size > 200000) {
        return(normalizePath(vis_file, winslash = "/"))
      }
      cli::cli_alert_info("Downloading CLIP ViT-B/32 models and tokenizer files ({row$size_mb} MB) to {.path {dir}}...")
      vis_url <- "https://huggingface.co/Xenova/clip-vit-base-patch32/resolve/main/onnx/vision_model_quantized.onnx"
      txt_url <- "https://huggingface.co/Xenova/clip-vit-base-patch32/resolve/main/onnx/text_model_quantized.onnx"
      vocab_url <- "https://huggingface.co/Xenova/clip-vit-base-patch32/raw/main/vocab.json"
      merges_url <- "https://huggingface.co/Xenova/clip-vit-base-patch32/raw/main/merges.txt"
      ok_vis <- .download_url_robust(vis_url, vis_file, min_size = 30000000)
      ok_txt <- .download_url_robust(txt_url, txt_file, min_size = 30000000)
      ok_vocab <- .download_url_robust(vocab_url, vocab_file, min_size = 500000)
      ok_merges <- .download_url_robust(merges_url, merges_file, min_size = 200000)
      if (!ok_vis || !ok_txt || !ok_vocab || !ok_merges) {
        stop("Failed to download CLIP models/tokenizer. Please check your internet connection.")
      }
      cli::cli_alert_success("CLIP ViT-B/32 successfully downloaded!")
      return(normalizePath(vis_file, winslash = "/"))
    }

    if (grepl("world", actual_model)) {
      world_pt_name <- if (grepl("\\.pt$", actual_model)) actual_model else paste0(actual_model, ".pt")
      if (actual_model %in% c("yolo-world", "yoloworld", "yolov8s-world")) world_pt_name <- "yolov8s-worldv2.pt"
      if (actual_model %in% c("yolov8m-world", "yolov8mworld")) world_pt_name <- "yolov8m-worldv2.pt"
      if (actual_model %in% c("yolov8l-world", "yolov8lworld")) world_pt_name <- "yolov8l-worldv2.pt"
      if (actual_model %in% c("yolov8x-world", "yolov8xworld")) world_pt_name <- "yolov8x-worldv2.pt"
      world_file <- file.path(dir, world_pt_name)
      if (!isTRUE(force) && file.exists(world_file) && file.info(world_file)$size > 15000000) {
        return(normalizePath(world_file, winslash = "/"))
      }
      size_disp <- if (!is.null(row$size_mb) && length(row$size_mb) > 0 && !is.na(row$size_mb[1])) paste0(" (", row$size_mb[1], " MB)") else ""
      cli::cli_alert_info("Downloading YOLO-World weights{size_disp} to {.path {dir}}...")
      world_url <- if (!is.null(row$url) && length(row$url) > 0 && !is.na(row$url[1]) && nzchar(row$url[1])) {
        row$url[1]
      } else {
        paste0("https://github.com/ultralytics/assets/releases/download/v8.2.0/", world_pt_name)
      }
      ok_world <- .download_url_robust(world_url, world_file, min_size = 15000000)
      if (!ok_world) {
        stop("Failed to download YOLO-World model. Please check your internet connection.")
      }
      cli::cli_alert_success("YOLO-World successfully downloaded!")
      return(normalizePath(world_file, winslash = "/"))
    }

    filename <- if (!is.null(row$filename) && length(row$filename) > 0 && !is.na(row$filename[1]) && nzchar(row$filename[1])) {
      row$filename[1]
    } else {
      paste0(actual_model, ".onnx")
    }
    dest_file <- file.path(dir, filename)

    if (!isTRUE(force) && file.exists(dest_file) && file.info(dest_file)$size > 1000) {
      return(normalizePath(dest_file, winslash = "/"))
    }

    url <- if (!is.null(row$url) && length(row$url) > 0 && !is.na(row$url[1]) && nzchar(row$url[1])) {
      row$url[1]
    } else {
      paste0(base_url, filename)
    }

    size_disp <- if (!is.null(row$size_mb) && length(row$size_mb) > 0 && !is.na(row$size_mb[1])) paste0(" (", row$size_mb[1], " MB)") else ""
    cli::cli_alert_info("Downloading model {.val {actual_model}}{size_disp} to {.path {dir}}...")

    dl_success <- .download_url_robust(url, dest_file, min_size = 1000)

    if (!isTRUE(dl_success) || !file.exists(dest_file) || file.info(dest_file)$size < 1000) {
      if (file.exists(dest_file)) unlink(dest_file)
      stop("Failed to download model ", actual_model, ". Please check your internet connection (e.g., Wi-Fi captive portal login).")
    }

    cli::cli_alert_success("Model {.val {actual_model}} successfully downloaded!")
    return(normalizePath(dest_file, winslash = "/"))
  }

  if (identical(model, "all") || ("all" %in% model)) {
    models_df <- pliman_available_models(dir = dir, refresh = refresh)
    model <- unique(models_df$model)
  }

  res <- vapply(model, download_single, character(1), USE.NAMES = FALSE)
  if (length(res) == 1L) {
    return(res[1])
  }
  return(res)
}

# ==============================================================================
# SECTION 2: ENVIRONMENT CONFIGURATION & DIRECT ONNX RUNTIME DOWNLOADER
# ==============================================================================

#' Directory for ONNX Runtime C++ Shared Libraries
#'
#' Gets or sets the local directory where the official Microsoft ONNX Runtime
#' C++ shared libraries (`onnxruntime.dll`, `libonnxruntime.so`, or `libonnxruntime.dylib`)
#' are saved.
#'
#' @param dir Optional custom directory path.
#' @return A character string with the path to the bin directory.
#' @export
#' @examples
#' pliman_onnx_dir()
pliman_onnx_dir <- function(dir = NULL) {
  if (!is.null(dir)) {
    dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(dir)
  }

  default_dir <- file.path(tools::R_user_dir("pliman", which = "data"), "bin")
  default_dir <- normalizePath(default_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(default_dir)) {
    dir.create(default_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(default_dir)
}

#' Download Microsoft ONNX Runtime C++ Library Directly
#'
#' Downloads and installs the official Microsoft ONNX Runtime C++ shared library
#' directly from Microsoft's GitHub releases, with **zero external package dependencies**
#' and **zero Python**.
#'
#' Automatically detects the host operating system (Windows, Linux, or macOS),
#' downloads the pre-compiled CPU binary (~15-60 MB), extracts the shared library
#' (`onnxruntime.dll` on Windows, `libonnxruntime.so` on Linux, `libonnxruntime.dylib` on macOS),
#' places it in `pliman`'s binary directory, and configures the environment.
#'
#' @param version Character string with the ONNX Runtime release version (default `"1.20.1"`).
#' @param dir Directory to store the library (default: `pliman_onnx_dir()`).
#' @param force Logical. If `TRUE`, forces re-download even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`. When `"gpu"` on Windows,
#'   downloads and configures `onnxruntime-directml.dll` and `DirectML.dll` from Microsoft's NuGet CDN.
#' @return A character string with the path to the installed shared library.
#' @export
#' @examples
#' \dontrun{
#'   pliman_download_onnx()
#'   pliman_download_onnx(engine = "gpu")
#' }
pliman_download_onnx <- function(version = "1.20.1",
                                 dir = pliman_onnx_dir(),
                                 force = FALSE,
                                 engine = c("gpu", "cpu")) {
  engine <- match.arg(engine)
  os <- tolower(Sys.info()[["sysname"]])

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  if (engine == "gpu") {
    if (os != "windows") {
      cli::cli_alert_warning("GPU acceleration via DirectML is currently only supported on Windows (DirectX 12). Falling back to CPU engine.")
      return(pliman_download_onnx(version = version, dir = dir, force = force, engine = "cpu"))
    }

    dest_dml_ort <- file.path(dir, "onnxruntime-directml.dll")
    dest_dml_core <- file.path(dir, "DirectML.dll")

    if (!isTRUE(force) &&
        file.exists(dest_dml_ort) && file.info(dest_dml_ort)$size > 1000000 &&
        file.exists(dest_dml_core) && file.info(dest_dml_core)$size > 1000000) {
      return(normalizePath(dest_dml_ort, winslash = "/"))
    }

    # Fast direct download of pre-extracted DLLs from repository
    base_url <- pliman_models_base_url()
    cli::cli_alert_info("Downloading DirectML (GPU) binaries from pliman repository...")
    ok_ort <- .download_url_robust(paste0(base_url, "onnxruntime-directml.dll"), dest_dml_ort, min_size = 1000000)
    ok_dml <- .download_url_robust(paste0(base_url, "DirectML.dll"), dest_dml_core, min_size = 1000000)

    if (!isTRUE(ok_ort) || !isTRUE(ok_dml)) {
      stop("Failed to download DirectML GPU libraries from ", base_url,
           ". Please check your connection or repository availability.")
    }

    cli::cli_alert_success("Microsoft ONNX Runtime DirectML (GPU) installed at {.path {dest_dml_ort}}!")
    return(normalizePath(dest_dml_ort, winslash = "/"))
  }

  lib_name <- switch(os,
    "windows" = "onnxruntime.dll",
    "linux"   = "libonnxruntime.so",
    "darwin"  = "libonnxruntime.dylib",
    stop("Unsupported operating system: ", os)
  )

  dest_lib <- file.path(dir, lib_name)
  onnxr_lib_dir <- file.path(tools::R_user_dir("onnxr", which = "data"), "lib")
  dest_onnxr_lib <- file.path(onnxr_lib_dir, lib_name)

  if (!isTRUE(force) && file.exists(dest_lib) && file.info(dest_lib)$size > 1000000) {
    if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(dest_onnxr_lib) || file.info(dest_onnxr_lib)$size < 1000000) {
      file.copy(dest_lib, dest_onnxr_lib, overwrite = TRUE)
    }
    Sys.setenv(ORT_ROOT = dir)
    return(normalizePath(dest_lib, winslash = "/"))
  }

  if (os == "windows") {
    base_url <- pliman_models_base_url()
    cli::cli_alert_info("Downloading ONNX Runtime CPU library from pliman repository...")
    ok_cpu <- .download_url_robust(paste0(base_url, "onnxruntime.dll"), dest_lib, min_size = 1000000)
    if (isTRUE(ok_cpu)) {
      if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(dest_lib, dest_onnxr_lib, overwrite = TRUE)
      Sys.setenv(ORT_ROOT = dir)
      cli::cli_alert_success("Microsoft ONNX Runtime CPU library installed at {.path {dest_lib}}!")
      return(normalizePath(dest_lib, winslash = "/"))
    }
  }

  cli::cli_alert_info("Downloading official Microsoft ONNX Runtime C++ library (v{version}) for {.val {os}}...")

  url <- switch(os,
    "windows" = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-win-x64-", version, ".zip"),
    "linux"   = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-linux-x64-", version, ".tgz"),
    "darwin"  = paste0("https://github.com/microsoft/onnxruntime/releases/download/v", version, "/onnxruntime-osx-universal-", version, ".tgz")
  )

  temp_archive <- file.path(tempdir(), paste0("ort_archive_", version, if (os == "windows") ".zip" else ".tgz"))
  old_timeout <- getOption("timeout")
  on.exit({
    options(timeout = old_timeout)
    if (file.exists(temp_archive)) unlink(temp_archive)
  }, add = TRUE)
  options(timeout = max(600, old_timeout))

  if (requireNamespace("curl", quietly = TRUE)) {
    curl::curl_download(url, destfile = temp_archive, quiet = FALSE, mode = "wb")
  } else {
    utils::download.file(url, destfile = temp_archive, mode = "wb", quiet = FALSE)
  }

  if (!file.exists(temp_archive) || file.info(temp_archive)$size < 100000) {
    stop("Failed to download ONNX Runtime from ", url, ". Please check your internet connection.")
  }

  cli::cli_alert_info("Extracting {.val {lib_name}}...")
  extract_tmp <- file.path(tempdir(), paste0("ort_extracted_", Sys.getpid()))
  dir.create(extract_tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(extract_tmp, recursive = TRUE), add = TRUE)

  if (os == "windows") {
    utils::unzip(temp_archive, exdir = extract_tmp)
  } else {
    utils::untar(temp_archive, exdir = extract_tmp)
  }

  found_files <- list.files(extract_tmp, pattern = paste0("^", gsub("\\.", "\\\\.", lib_name), "$"), recursive = TRUE, full.names = TRUE)
  if (length(found_files) == 0) {
    found_files <- list.files(extract_tmp, pattern = "onnxruntime", recursive = TRUE, full.names = TRUE)
    found_files <- found_files[grepl(paste0("\\", tools::file_ext(lib_name), "$"), found_files)]
  }

  if (length(found_files) == 0) {
    stop("Could not locate ", lib_name, " in the downloaded archive.")
  }

  src_lib <- found_files[1]
  file.copy(src_lib, dest_lib, overwrite = TRUE)

  tryCatch({
    if (!dir.exists(onnxr_lib_dir)) dir.create(onnxr_lib_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(src_lib, dest_onnxr_lib, overwrite = TRUE)
  }, error = function(e) NULL, warning = function(w) NULL)

  Sys.setenv(ORT_ROOT = dir)
  cli::cli_alert_success("Microsoft ONNX Runtime C++ library installed at {.path {dest_lib}}!")
  return(normalizePath(dest_lib, winslash = "/"))
}

#' Proprietary ONNX Runtime Installer for pliman
#'
#' Downloads and configures the official Microsoft ONNX Runtime C++ shared library
#' directly into pliman's local bin directory, requiring **zero external packages**
#' and **zero Python**.
#'
#' @param version Character string with the ONNX Runtime release version (default `"1.20.1"`).
#' @param force Logical. If `TRUE`, forces re-download even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return A character string with the path to the installed shared library.
#' @export
#' @examples
#' \dontrun{
#'   onnx_install()
#'   onnx_install(engine = "gpu")
#' }
onnx_install <- function(version = "1.20.1", force = FALSE, engine = c("gpu", "cpu")) {
  pliman_download_onnx(version = version, force = force, engine = engine)
}

#' @rdname onnx_install
#' @export
pliman_install_onnx <- function(version = "1.20.1", force = FALSE, engine = c("gpu", "cpu")) {
  onnx_install(version = version, force = force, engine = engine)
}

#' Locate ONNX Runtime C++ Shared Library
#'
#' Finds the local path to the ONNX Runtime dynamic library (`onnxruntime.dll` on Windows,
#' `libonnxruntime.so` on Linux, `libonnxruntime.dylib` on macOS), or `onnxruntime-directml.dll`
#' when `engine = "gpu"` on Windows.
#'
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return Character path to the dynamic library or `NULL` if not found.
#' @export
pliman_onnx_lib_path <- function(engine = c("gpu", "cpu")) {
  engine <- match.arg(engine)
  os <- tolower(Sys.info()[["sysname"]])

  if (os == "windows") {
    p_gpu <- file.path(pliman_onnx_dir(), "onnxruntime-directml.dll")
    if (file.exists(p_gpu)) {
      return(normalizePath(p_gpu, winslash = "/"))
    }
    if (engine == "gpu") {
      return(NULL)
    }
  } else if (engine == "gpu") {
    cli::cli_alert_warning("DirectML GPU acceleration is only supported on Windows (DirectX 12). Falling back to 'cpu'.")
    engine <- "cpu"
  }

  lib_name <- switch(os,
    "windows" = "onnxruntime.dll",
    "linux"   = "libonnxruntime.so",
    "darwin"  = "libonnxruntime.dylib",
    "onnxruntime.dll"
  )

  p1 <- file.path(pliman_onnx_dir(), lib_name)
  if (file.exists(p1)) return(normalizePath(p1, winslash = "/"))

  ort_root <- Sys.getenv("ORT_ROOT", "")
  if (nzchar(ort_root)) {
    p2 <- file.path(ort_root, lib_name)
    if (file.exists(p2)) return(normalizePath(p2, winslash = "/"))
    p2b <- file.path(ort_root, "lib", lib_name)
    if (file.exists(p2b)) return(normalizePath(p2b, winslash = "/"))
  }

  p3 <- file.path(tools::R_user_dir("onnxr", which = "data"), "lib", lib_name)
  if (file.exists(p3)) return(normalizePath(p3, winslash = "/"))

  return(NULL)
}

pliman_onnx_library_path <- function(engine = c("gpu", "cpu")) {
  pliman_onnx_lib_path(engine = engine)
}


#' GPU Information and DirectML Status in pliman
#'
#' Queries the system's graphics hardware (via DirectX DXGI on Windows) to list
#' available GPU adapters, dedicated video memory (VRAM), and DirectML status.
#'
#' @return A data frame containing information on detected GPU adapters, or a message if no GPU is available.
#' @export
#' @examples
#' \dontrun{
#'   pliman_gpu_info()
#' }
pliman_gpu_info <- function() {
  info <- pliman_gpu_info_cpp()
  if (!isTRUE(info$available) || length(info$device_id) == 0) {
    cli::cli_alert_info("No DirectX 12 compatible GPU adapters detected or running on non-Windows OS.")
    return(invisible(NULL))
  }

  df <- data.frame(
    device_id = info$device_id,
    name = info$name,
    vram_mb = round(info$vram_mb, 1),
    is_dedicated = info$is_dedicated,
    is_default = (info$device_id == info$default_device_id),
    stringsAsFactors = FALSE
  )

  cli::cli_h1("Detected GPU Adapters (DirectX 12 / DirectML)")
  for (i in seq_len(nrow(df))) {
    prefix <- if (df$is_default[i]) "* (Default) " else "  "
    type_str <- if (df$is_dedicated[i]) "Dedicated" else "Integrated"
    cli::cli_alert_info("{prefix}[Device {df$device_id[i]}] {df$name[i]} ({type_str}, {df$vram_mb[i]} MB VRAM)")
  }

  dml_lib <- pliman_onnx_lib_path(engine = "gpu")
  if (!is.null(dml_lib) && file.exists(dml_lib)) {
    cli::cli_alert_success("DirectML GPU library: Installed at {.path {dml_lib}}")
  } else {
    cli::cli_alert_warning("DirectML GPU library: Not installed. Run {.code onnx_install(engine = 'gpu')} to configure.")
  }

  invisible(df)
}

#' Configure Complete Deep Learning Environment for pliman (100% Native C++)
#'
#' Sets up the entire zero-Python, zero-external-package Deep Learning environment
#' in a single call with interactive CLI feedback.
#'
#' This function:
#' 1. Creates and verifies the local model directory (`tools::R_user_dir("pliman", "data")/models`).
#' 2. Downloads and verifies the official Microsoft ONNX Runtime C++ shared library directly.
#' 3. Downloads the specified pre-trained models (`u2netp` and `isnet-general-use` by default).
#' 4. Runs a self-test inference using the native C++ engine to confirm 100% operational status.
#'
#' @param models Character vector of models to pre-download. Defaults to `"u2netp"`.
#'   Pass `"all"` to download all 5 models, or `NULL` to only configure the engine without downloading models.
#' @param dir Directory to store models (default: `pliman_model_dir()`).
#' @param force Logical. Force re-download and re-configuration even if already present (default `FALSE`).
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"`.
#' @return Logical `TRUE` invisibly on success.
#' @export
#' @examples
#' \dontrun{
#'   pliman_configure_dl()
#'   pliman_configure_dl(engine = "gpu")
#' }
pliman_configure_dl <- function(models = "u2netp",
                                dir = pliman_model_dir(),
                                force = FALSE,
                                engine = c("gpu", "cpu")) {
  engine <- match.arg(engine)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(600, old_timeout))

  cli::cli_h1("Configuring pliman Deep Learning Environment (100% Native C++)")

  # 1. Check and configure model directory
  cli::cli_h2("1. Model Storage Directory")
  model_path <- pliman_model_dir(dir)
  cli::cli_alert_success("Model directory: {.path {model_path}}")

  # 2. Check and configure ONNX Runtime C++ engine directly
  cli::cli_h2("2. Microsoft ONNX Runtime C++ Engine (Zero-Package Native)")
  lib_path <- onnx_install(force = force, engine = engine)
  cli::cli_alert_success("ONNX Runtime C++ library verified: {.path {lib_path}}")

  # 3. Download Requested Models
  if (!is.null(models)) {
    cli::cli_h2("3. Pre-trained Neural Models")
    if (identical(models, "all")) {
      models <- c("u2netp", "isnet-general-use", "rmbg-1.4", "silueta", "u2net")
    }

    for (m in models) {
      tryCatch({
        pliman_download_model(model = m, dir = model_path, force = force)
      }, error = function(e) {
        cli::cli_alert_danger("Failed to download model {.val {m}}: {e$message}")
      })
    }
  }

  # 4. Perform a quick verification test
  cli::cli_h2("4. Verification & Self-Test")
  test_model <- if ("u2netp" %in% models) "u2netp" else models[1]
  test_file <- file.path(model_path, paste0(test_model, ".onnx"))

  if (file.exists(test_file)) {
    cli::cli_alert_info("Running self-test with {.val {test_model}} in native C++...")
    test_ok <- tryCatch({
      dummy_tensor <- numeric(1 * 3 * 320 * 320)
      attr(dummy_tensor, "dims") <- as.integer(c(1, 3, 320, 320))
      res <- .run_onnx_inference(dummy_tensor, test_file, target_size = 320, engine = engine)
      if (is.matrix(res) && nrow(res) == 320 && ncol(res) == 320) {
        TRUE
      } else {
        FALSE
      }
    }, error = function(e) {
      cli::cli_alert_warning("Self-test error: {e$message}")
      FALSE
    })

    if (test_ok) {
      cli::cli_alert_success("Self-test passed! Native C++ neural inference is working at full speed.")
    }
  }

  cli::cli_h2("Configuration Complete!")
  cli::cli_alert_info("You can now run: {.code mask <- image_binary_dl(img, model = 'u2netp')}")
  cli::cli_alert_info("Or remove background: {.code img_tb <- image_remove_bg_dl(img, model = 'u2netp')}")
  return(invisible(TRUE))
}

# ==============================================================================
# SECTION 3: INTERNAL TENSOR PREPROCESSING & INFERENCE RUNNER
# ==============================================================================

# Bilinear interpolation resize for numeric/raw matrix via high-speed C++
.bilinear_resize_2d <- function(mat, out_w, out_h) {
  in_w <- nrow(mat)
  in_h <- ncol(mat)
  if (in_w == out_w && in_h == out_h) return(mat)
  image_resize_cpp(mat, as.integer(out_w), as.integer(out_h), filter = 1L)
}

# Fast polygon decimation for screen display (preserves full resolution in returned contours)
.decimate_poly <- function(p, max_pts = 2500L) {
  n <- nrow(p)
  if (n <= max_pts) return(p)
  step <- ceiling(n / max_pts)
  idx <- c(seq(1L, n, by = step), 1L)
  p[idx, , drop = FALSE]
}

.preprocess_nchw <- function(mat,
                             target_size = 320,
                             mean = c(0.485, 0.456, 0.406),
                             std = c(0.229, 0.224, 0.225),
                             letterbox = FALSE) {
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  # Extract R, G, B normalized [0, 1]
  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  dim(R) <- dim(G) <- dim(B) <- c(orig_w, orig_h)

  if (isTRUE(letterbox)) {
    scale <- target_size / max(orig_w, orig_h)
    new_w <- max(1L, round(orig_w * scale))
    new_h <- max(1L, round(orig_h * scale))

    R_res <- .bilinear_resize_2d(R, new_w, new_h)
    G_res <- .bilinear_resize_2d(G, new_w, new_h)
    B_res <- .bilinear_resize_2d(B, new_w, new_h)

    # Pad at top-left of canvas target_size x target_size [width, height]
    R_pad <- matrix(0.0, nrow = target_size, ncol = target_size)
    G_pad <- matrix(0.0, nrow = target_size, ncol = target_size)
    B_pad <- matrix(0.0, nrow = target_size, ncol = target_size)

    R_pad[1:new_w, 1:new_h] <- R_res
    G_pad[1:new_w, 1:new_h] <- G_res
    B_pad[1:new_w, 1:new_h] <- B_res

    R_norm <- (R_pad - mean[1]) / std[1]
    G_norm <- (G_pad - mean[2]) / std[2]
    B_norm <- (B_pad - mean[3]) / std[3]

    vec_R <- as.numeric(R_norm)
    vec_G <- as.numeric(G_norm)
    vec_B <- as.numeric(B_norm)

    flat_tensor <- c(vec_R, vec_G, vec_B)
    attr(flat_tensor, "dims") <- as.integer(c(1, 3, target_size, target_size))
    attr(flat_tensor, "letterbox_dims") <- c(new_w, new_h)
    return(flat_tensor)
  }

  # Resize each channel to target_size x target_size
  R_res <- .bilinear_resize_2d(R, target_size, target_size)
  G_res <- .bilinear_resize_2d(G, target_size, target_size)
  B_res <- .bilinear_resize_2d(B, target_size, target_size)

  # Normalize
  R_norm <- (R_res - mean[1]) / std[1]
  G_norm <- (G_res - mean[2]) / std[2]
  B_norm <- (B_res - mean[3]) / std[3]

  # In pliman, M[x, y] is column-major with x as row and y as col.
  # Reading column-by-column traverses y=1 (all x), then y=2 (all x),
  # which corresponds to the C-order NCHW row-major layout expected by ONNX.
  vec_R <- as.numeric(R_norm)
  vec_G <- as.numeric(G_norm)
  vec_B <- as.numeric(B_norm)

  # Combined NCHW flat vector [1, 3, target_size, target_size]
  flat_tensor <- c(vec_R, vec_G, vec_B)
  attr(flat_tensor, "dims") <- as.integer(c(1, 3, target_size, target_size))

  return(flat_tensor)
}

.run_onnx_inference <- function(tensor, model_path, target_size = 320, threads = 0,
                                engine = c("gpu", "cpu"), device_id = -1) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_path) || !file.exists(lib_path)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_path <- onnx_install(engine = engine)
  }

  tensor_dims <- attr(tensor, "dims")
  if (is.null(tensor_dims)) {
    tensor_dims <- as.integer(c(1, 3, target_size, target_size))
  }
  tensor_vec <- as.numeric(tensor)

  pred_mat <- run_onnx_inference_cpp(
    tensor_vec = tensor_vec,
    tensor_dims = tensor_dims,
    model_path = normalizePath(model_path, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_path, winslash = "/", mustWork = FALSE),
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  return(pred_mat)
}

# Internal core mask computation helper
.compute_dl_mask <- function(mat,
                             model = "u2netp",
                             threshold = 0.5,
                             fill_hull = TRUE,
                             filter = 0,
                             erode = 0,
                             dilate = 0,
                             opening = 0,
                             closing = 0,
                             min_area = 0,
                             invert = FALSE,
                             pick_object = FALSE,
                             prompt = NULL,
                             exemplar = FALSE,
                             threads = 0,
                             engine = c("gpu", "cpu"),
                             device_id = -1,
                             verbose = TRUE,
                             dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  model <- .resolve_model_name(model[1])
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Interactive picking of prompt point(s)
  is_text_prompt <- is.character(prompt) && length(prompt) >= 1 && !all(prompt %in% c("center", "box", "exemplar"))
  is_sam <- grepl("sam", model, ignore.case = TRUE)
  use_persam <- isTRUE(exemplar) || (model == "persam") || (is_sam && identical(prompt, "exemplar"))
  use_grounded_sam <- !use_persam && ((model == "grounded-sam") || (is_sam && is_text_prompt))

  if (use_persam) {
    p_res <- .run_persam(
      mat = mat,
      exemplar_points = if (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt)) prompt else NULL,
      threshold = threshold,
      threads = threads,
      engine = engine,
      device_id = device_id,
      fill_hull = fill_hull,
      filter = filter,
      erode = erode,
      dilate = dilate,
      opening = opening,
      closing = closing,
      min_area = min_area,
      invert = invert,
      mask = TRUE,
      verbose = verbose,
      dir = dir
    )
    if (isTRUE(verbose)) {
      .print_detection_summary(p_res$summary)
    }
    return(p_res$mask * 1.0)
  }

  if (use_grounded_sam) {
    if (is.null(prompt)) {
      cli::cli_abort("Model {.val grounded-sam} requires a text {.arg prompt} (e.g., {.code prompt = 'leaf'}).")
    }
    gs_res <- .run_grounded_sam(
      mat = mat,
      prompt = prompt,
      threshold = threshold,
      threads = threads,
      engine = engine,
      device_id = device_id,
      fill_hull = fill_hull,
      filter = filter,
      erode = erode,
      dilate = dilate,
      opening = opening,
      closing = closing,
      min_area = min_area,
      invert = invert,
      mask = TRUE,
      verbose = verbose,
      dir = dir
    )
    if (isTRUE(verbose)) {
      .print_detection_summary(gs_res$summary)
    }
    return(gs_res$mask * 1.0)
  }

  if (isTRUE(pick_object)) {
    cli::cli_alert_info("Click on the object(s) of interest in the plot window. Press <Esc> or right-click when finished.")
    plot(as_image(mat))
    pts <- tryCatch(graphics::locator(n = 512, type = "p", col = "red", pch = 19), error = function(e) NULL)
    if (!is.null(pts) && length(pts$x) > 0) {
      prompt <- cbind(pts$x, pts$y)
    } else {
      cli::cli_alert_warning("No point selected. Defaulting to center of image.")
      prompt <- c(orig_w / 2, orig_h / 2)
    }
  }

  if (is_sam) {
    if (is.null(prompt)) {
      prompt <- c(orig_w / 2, orig_h / 2)
    }

    # Ensure SAM models are downloaded
    model_file <- pliman_download_model(model = model, dir = dir)
    enc_path <- file.path(dir, "sam2.1.encoder.onnx")
    dec_path <- file.path(dir, "sam2.1.decoder.onnx")

    # SAM uses standard ImageNet normalization and 1024x1024
    tensor <- .preprocess_nchw(mat,
                               target_size = 1024L,
                               mean = c(0.485, 0.456, 0.406),
                               std = c(0.229, 0.224, 0.225),
                               letterbox = FALSE)

    # Convert prompt coordinates to 1024x1024 space
    if (is.character(prompt)) {
      if (prompt == "center") {
        pts_x <- 512.0
        pts_y <- 512.0
        pts_lbl <- 1L
      } else { # "box"
        pts_x <- c(10.0, 1014.0)
        pts_y <- c(10.0, 1014.0)
        pts_lbl <- c(2L, 3L)
      }
    } else if (is.matrix(prompt) || is.data.frame(prompt)) {
      pts_x <- (as.numeric(prompt[, 1]) / orig_w) * 1024.0
      pts_y <- (as.numeric(prompt[, 2]) / orig_h) * 1024.0
      pts_lbl <- rep(1L, length(pts_x))
    } else if (length(prompt) == 4) { # box: c(xmin, ymin, xmax, ymax)
      pts_x <- c(prompt[1] / orig_w * 1024.0, prompt[3] / orig_w * 1024.0)
      pts_y <- c(prompt[2] / orig_h * 1024.0, prompt[4] / orig_h * 1024.0)
      pts_lbl <- c(2L, 3L)
    } else { # single point: c(x, y)
      pts_x <- (as.numeric(prompt[1]) / orig_w) * 1024.0
      pts_y <- (as.numeric(prompt[2]) / orig_h) * 1024.0
      pts_lbl <- 1L
    }

    lib_path <- pliman_onnx_lib_path(engine = engine)
    if (is.null(lib_path) || !file.exists(lib_path)) {
      cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
      lib_path <- onnx_install(engine = engine)
    }

    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Running SAM 2.1 inference [{toupper(engine)}]...",
        msg_done = "SAM 2.1 inference complete"
      )
    }

    pred_raw <- run_sam2_inference_cpp(
      tensor_vec = as.numeric(tensor),
      points_x = pts_x,
      points_y = pts_y,
      point_labels = pts_lbl,
      encoder_path = normalizePath(enc_path, winslash = "/", mustWork = FALSE),
      decoder_path = normalizePath(dec_path, winslash = "/", mustWork = FALSE),
      lib_path = normalizePath(lib_path, winslash = "/", mustWork = FALSE),
      num_threads = as.integer(threads),
      use_gpu = use_gpu,
      device_id = as.integer(device_id)
    )

    mask_prob <- .bilinear_resize_2d(pred_raw, orig_w, orig_h)
  } else {
    # Ensure model is downloaded
    model_file <- pliman_download_model(model = model, dir = dir)

    # Determine target input size, normalization, and letterboxing
    is_wbg <- grepl("withoutbg", model)
    if (model %in% c("rmbg-1.4", "isnet-general-use")) {
      target_size <- 1024L
      mean_val <- c(0.5, 0.5, 0.5)
      std_val <- c(1.0, 1.0, 1.0)
      use_letterbox <- FALSE
    } else if (model == "ben2") {
      target_size <- 1024L
      mean_val <- c(0, 0, 0)
      std_val <- c(1, 1, 1)
      use_letterbox <- FALSE
    } else if (is_wbg) {
      target_size <- 448L
      mean_val <- c(0, 0, 0)
      std_val <- c(1, 1, 1)
      use_letterbox <- TRUE
    } else {
      models_df <- pliman_available_models(dir = dir)
      row <- models_df[models_df$model == model, ]
      target_size <- if (nrow(row) > 0) as.integer(strsplit(as.character(row$input_size), "x")[[1]][1]) else 320L
      mean_val <- c(0.485, 0.456, 0.406)
      std_val <- c(0.229, 0.224, 0.225)
      use_letterbox <- FALSE
    }

    # 1. Pre-process to NCHW normalized tensor (row-major flat vector)
    tensor <- .preprocess_nchw(mat,
                               target_size = target_size,
                               mean = mean_val,
                               std = std_val,
                               letterbox = use_letterbox)

    # 2. Run ONNX Inference
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Running neural inference with {.val {model}} [{toupper(engine)}]...",
        msg_done = "Inference with {.val {model}} complete"
      )
    }

    pred_raw <- tryCatch({
      .run_onnx_inference(tensor, model_file, target_size = target_size, threads = threads, engine = engine, device_id = device_id)
    }, error = function(e) {
      if (engine == "gpu" && grepl("Recursos de mem[o\u00f3]ria|out of memory|8007000E", e$message, ignore.case = TRUE)) {
        cli::cli_alert_danger("Insufficient GPU VRAM for model {.val {model}} at resolution {target_size}x{target_size}.")
        cli::cli_alert_info("Tip: Use {.code model = 'rmbg-1.4'}, {.code model = 'birefnet-lite'}, or {.code model = 'u2netp'} on GPU, or switch to {.code engine = 'cpu'}.")
      }
      stop(e)
    })

    # Normalize raw predictions to [0, 1] probability
    p_min <- min(pred_raw, na.rm = TRUE)
    p_max <- max(pred_raw, na.rm = TRUE)
    if (p_min < 0 || p_max > 1) {
      prob_map <- 1.0 / (1.0 + exp(-pred_raw))
    } else {
      prob_map <- pred_raw
    }

    # 3. Un-letterbox if applicable and resize probability map back to original image dimensions
    if (isTRUE(use_letterbox)) {
      lb_dims <- attr(tensor, "letterbox_dims")
      prob_crop <- prob_map[1:lb_dims[1], 1:lb_dims[2], drop = FALSE]
      mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
    } else {
      mask_prob <- .bilinear_resize_2d(prob_map, orig_w, orig_h)
    }
  }

  # 4. Threshold into initial binary mask
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing binary mask...",
      msg_done = "Mask post-processing complete"
    )
  }

  mask <- mask_prob >= threshold

  # 5. Morphological Post-Processing on the Mask
  # 5.1 Fill internal holes in the mask
  if (isTRUE(fill_hull)) {
    mask <- fill_holes_cpp(mask)
    mask <- (mask > 0)
  }

  # 5.2 Median filter smoothing to reduce noise/speckles
  if (filter > 0) {
    mask_arr <- as.array(mask)
    if (length(dim(mask_arr)) == 2) dim(mask_arr) <- c(orig_w, orig_h, 1L)
    f_res <- median_filter_binary_cpp(mask_arr, orig_w, orig_h, 1L, as.integer(filter))
    mask <- (f_res[, , 1] > 0)
  }

  # 5.3 Morphological erosion (shrinks mask boundaries)
  if (erode > 0) {
    mask <- erode_cpp(mask, raio = as.integer(erode))
    mask <- (mask > 0)
  }

  # 5.4 Morphological dilation (expands mask boundaries)
  if (dilate > 0) {
    mask <- dilate_cpp(mask, raio = as.integer(dilate))
    mask <- (mask > 0)
  }

  # 5.5 Morphological opening (erode then dilate)
  if (opening > 0) {
    mask <- erode_cpp(mask, raio = as.integer(opening))
    mask <- dilate_cpp(mask, raio = as.integer(opening))
    mask <- (mask > 0)
  }

  # 5.6 Morphological closing (dilate then erode)
  if (closing > 0) {
    mask <- dilate_cpp(mask, raio = as.integer(closing))
    mask <- erode_cpp(mask, raio = as.integer(closing))
    mask <- (mask > 0)
  }

  # 5.7 Filter small artifacts / retain largest object
  if (min_area <= 0) {
    lbls <- bwlabel_cpp(mask)
    if (max(lbls) > 0) {
      tbl <- table(lbls[lbls > 0])
      largest_id <- as.integer(names(which.max(tbl)))
      mask <- (lbls == largest_id)
    }
  } else {
    lbls <- bwlabel_cpp(mask)
    if (max(lbls) > 0) {
      tbl <- table(lbls[lbls > 0])
      keep_ids <- as.integer(names(tbl[tbl >= min_area]))
      mask <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
    }
  }

  # 5.8 Invert mask if requested
  if (isTRUE(invert)) {
    mask <- !mask
  }

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }

  return(mask)
}

# ==============================================================================
# SECTION 4: USER-FACING SEGMENTATION, BINARIZATION & BACKGROUND REMOVAL
# ==============================================================================

#' Deep Learning Binary Mask Extraction
#'
#' Extracts a pixel-accurate binary mask (foreground = 1, background = 0) from an image using
#' state-of-the-art Deep Learning neural networks executed in ONNX format via native C++ bindings.
#' Supports salient object detection, dichotomous image segmentation, high-resolution boundary matting,
#' foundation segmentation models (SAM 2.1 / SAM 3.1), zero-shot text-prompted instance segmentation (Grounded-SAM),
#' and one-shot visual exemplar segmentation (PerSAM).
#'
#' @param img An `image` object, 2D grayscale matrix, 3D color array (RGB/RGBA), or a `list`
#'   of images/arrays. When a list is passed, each element is automatically processed in batch,
#'   returning a corresponding list of binary masks.
#' @param model Character string specifying the pre-trained neural network architecture to use:
#'   * **General Salient Object Detection & Background Extraction:**
#'     - `"u2netp"` (default): Ultra-lightweight U2-Net Portable (~4.6 MB, 320x320 input). Highly recommended
#'       for rapid CPU batch processing on standard laptops without requiring a dedicated GPU.
#'     - `"ben2"` (or `"ben"`): Boundary-aware Extraction Network Base (~212.6 MB, 1024x1024 input). High-accuracy
#'       crisp boundary extraction, exceptional for fine leaf contours, intricate roots, thin petioles, and plant serrations.
#'     - `"rmbg-1.4"` (or `"rmbg"`): BRIA AI RMBG 1.4 (~176.1 MB, 1024x1024 input). Enterprise-grade salient foreground
#'       extractor trained on diverse commercial datasets, exceptionally robust against studio lighting, reflections, and shadows.
#'     - `"rmbg-2.0"`: BRIA AI RMBG 2.0 (~976.9 MB, 1024x1024 input). Next-generation BiRefNet-based model for maximum boundary
#'       detail, fine textures, and translucent edge matting.
#'     - `"isnet-general-use"`: Intermediate Supervision Network (~178.6 MB, 1024x1024 input). High-precision boundary matting,
#'       ideal for subtle leaf lesions, chlorotic halos, and fine botanical structures.
#'     - `"withoutbg"`: withoutBG Open Weights (~433.4 MB, 448x448 input). Combines DepthAnythingV2 3D depth features with
#'       ConvNeXt-fused U-Net matting, excelling in depth-cluttered backgrounds.
#'     - `"birefnet-lite"` (or `"birefnet"`): Bilateral Reference Network Lite (~213.6 MB, 1024x1024 input). Bilateral reference
#'       architecture specialized in dichotomous image segmentation and complex shapes.
#'     - `"silueta"`: Compact silhouette extraction model (~44.2 MB, 320x320 input). Fast inference with solid accuracy on edge/CPU hardware.
#'     - `"u2net"`: Full-depth U2-Net architecture (~176.3 MB, 320x320 input). Deep multi-scale salient object detection model.
#'   * **Zero-Shot & Interactive Foundation Models:**
#'     - `"grounded-sam"`: Zero-shot open-vocabulary text-prompted instance segmentation (~314.8 MB: Grounding DINO Tiny + SAM 2.1).
#'       Detects and segments specific concepts described in natural language (e.g., `prompt = "leaf"`, `prompt = "fruit"`, `prompt = "insect"`).
#'     - `"persam"`: Personalize Segment Anything via SAM 2.1 embeddings (~120.2 MB, 1024x1024 input). Segments all instances
#'       visually similar to one or more user-clicked exemplar objects using Multi-Prototype Cosine Similarity Max-Pooling.
#'     - `"sam2.1"` (or `"sam2"`): Meta AI Segment Anything 2.1 (Hiera-Tiny, ~120.2 MB, 1024x1024 input). Foundation vision transformer
#'       supporting point prompts, box prompts, or interactive clicking (`pick_object = TRUE`).
#'     - `"sam3.1"` (or `"sam3"`): Meta AI Segment Anything 3.1 (~868.1 MB, 1024x1024 input). High-capacity concept-driven
#'       foundation segmentation model.
#' @param threshold Numeric value in `[0, 1]` specifying the cutoff threshold applied to the model's
#'   continuous probability/logit output map. Default is `0.5`. Increasing this threshold (e.g., `0.7`) makes
#'   segmentation more conservative (avoiding background noise), whereas decreasing it (e.g., `0.3`) includes
#'   faint, thin, or translucent object edges.
#' @param fill_hull Logical. If `TRUE` (default), automatically fills internal holes and voids
#'   (such as specular glare or reflective spots) within the segmented foreground object using morphological
#'   hole reconstruction (`fill_holes_cpp()`).
#' @param filter Integer specifying the window radius for binary median filtering (`median_filter_binary_cpp()`).
#'   Default is `0` (disabled). A positive value (e.g., `2` or `3`) eliminates salt-and-pepper noise and
#'   smoothes jagged boundaries without blurring sharp object corners.
#' @param erode Integer specifying the radius for morphological erosion (`erode_cpp()`). Default is `0` (disabled).
#'   Shrinks foreground boundaries inward by the given pixel radius, useful for severing thin touching bridges
#'   between adjacent objects or removing outer halo artifacts.
#' @param dilate Integer specifying the radius for morphological dilation (`dilate_cpp()`). Default is `0` (disabled).
#'   Expands foreground boundaries outward by the given pixel radius, useful for restoring peripheral details
#'   or slightly undersized masks.
#' @param opening Integer specifying the radius for morphological opening (erosion followed by dilation).
#'   Default is `0` (disabled). Removes small stray islands and thin protrusions while preserving overall
#'   object geometry.
#' @param closing Integer specifying the radius for morphological closing (dilation followed by erosion).
#'   Default is `0` (disabled). Fuses narrow internal cracks, gaps, and small indentations along the perimeter.
#' @param min_area Integer minimum object size in pixels evaluated via connected component labeling (`bwlabel_cpp()`).
#'   If positive, discards all foreground connected components smaller than this area. If `0` (default),
#'   automatically isolates and retains only the single largest connected component in the image.
#' @param invert Logical. If `TRUE`, inverts the output binary mask such that foreground becomes `0` and
#'   background becomes `1`. Default is `FALSE`.
#' @param pick_object Logical. If `TRUE`, activates an interactive graphics prompt via `graphics::locator()`,
#'   allowing the user to directly click on the object(s) of interest in the plot window to guide SAM models.
#'   Press `<Esc>` or right-click when finished. Default is `FALSE`.
#' @param prompt Flexible prompt specification for SAM, Grounded-SAM, and PerSAM models:
#'   * `NULL` (default): Defaults to the center of the image `c(width / 2, height / 2)` for SAM models.
#'   * Character string (e.g., `"leaf"`, `"fruit"`, `"weed"`): Triggers Grounded-SAM open-vocabulary object detection.
#'   * Numeric vector `c(x, y)`: Coordinates of a single foreground point prompt or single visual exemplar.
#'   * Two-column numeric matrix or data frame `cbind(x, y)`: Specifies multiple positive point prompts or multiple exemplar prototypes.
#'   * Numeric vector of length 4 `c(xmin, ymin, xmax, ymax)`: Specifies a bounding box prompt in pixel coordinates.
#'   * `"center"`: Uses the image center coordinates `c(width / 2, height / 2)`.
#'   * `"box"`: Uses a bounding box spanning the full image canvas.
#' @param exemplar Logical. If `TRUE` (or if `model = "persam"`), activates one-shot visual exemplar segmentation (PerSAM)
#'   using SAM 2.1 embeddings. Default is `FALSE`.
#' @param threads Integer specifying the number of threads for parallel ONNX Runtime execution. Default is `0`,
#'   which automatically detects and utilizes all available logical CPU cores on the system.
#' @param engine Character specifying the backend execution engine: `"cpu"` (default) or `"gpu"`. On Windows,
#'   `"gpu"` utilizes Microsoft DirectML over DirectX 12 hardware (compatible with NVIDIA, AMD, and Intel GPUs
#'   without requiring CUDA toolkit installation).
#' @param device_id Integer GPU device adapter index. Default is `-1`, which automatically selects the high-performance
#'   discrete GPU adapter with the maximum dedicated video memory (VRAM).
#' @param verbose Logical. If `TRUE` (default), displays step-by-step progress feedback and timing indicators
#'   via `cli::cli_progress_step()`. Set to `FALSE` for silent execution in automated pipelines.
#' @param plot Logical. If `TRUE` (default), renders the resulting binary mask to the active graphics device.
#' @param dir Character string specifying the local directory where ONNX model files are stored and cached.
#'   Defaults to `pliman_model_dir()`.
#' @param ... Additional arguments passed down to internal plotting methods.
#'
#' @return A binary `image` object (Grayscale, values 0 and 1) or a list of binary `image` objects.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # Fast binarization on CPU (ultra-lightweight, 4.6 MB)
#'   mask <- image_binary_dl(img, model = "u2netp")
#'
#'   # Sub-pixel boundary segmentation for fine veins and serrations
#'   mask_ben <- image_binary_dl(img, model = "ben2")
#'
#'   # Interactive prompt: click directly on objects of interest
#'   mask_interactive <- image_binary_dl(img, model = "sam2.1", pick_object = TRUE)
#'
#'   # Open-vocabulary text-prompted binarization
#'   mask_weed <- image_binary_dl(img, model = "grounded-sam", prompt = "weed")
#'
#'   # One-shot visual exemplar segmentation (PerSAM)
#'   mask_exemplar <- image_binary_dl(img, model = "persam", prompt = c(150, 200))
#' }
image_binary_dl <- function(img,
                            model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                            threshold = 0.5,
                            fill_hull = TRUE,
                            filter = 0,
                            erode = 0,
                            dilate = 0,
                            opening = 0,
                            closing = 0,
                            min_area = 0,
                            invert = FALSE,
                            pick_object = FALSE,
                            prompt = NULL,
                            exemplar = FALSE,
                            threads = 0,
                            engine = c("gpu", "cpu"),
                            device_id = -1,
                            verbose = TRUE,
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_binary_dl(x, model = model, threshold = threshold,
                      fill_hull = fill_hull, filter = filter,
                      erode = erode, dilate = dilate,
                      opening = opening, closing = closing,
                      min_area = min_area, invert = invert,
                      pick_object = pick_object, prompt = prompt,
                      exemplar = exemplar,
                      threads = threads,
                      engine = engine, device_id = device_id,
                      verbose = verbose,
                      plot = FALSE, dir = dir)
    })
    return(res)
  }

  if (is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.null(dim(mat)) || length(dim(mat)) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  mask <- .compute_dl_mask(mat = mat,
                           model = model,
                           threshold = threshold,
                           fill_hull = fill_hull,
                           filter = filter,
                           erode = erode,
                           dilate = dilate,
                           opening = opening,
                           closing = closing,
                           min_area = min_area,
                           invert = invert,
                           pick_object = pick_object,
                           prompt = prompt,
                           exemplar = exemplar,
                           threads = threads,
                           engine = engine,
                           device_id = device_id,
                           verbose = verbose,
                           dir = dir)

  mask_img <- as_image(mask, colormode = "Grayscale")

  if (isTRUE(plot)) {
    plot(mask_img, ...)
  }

  invisible(mask_img)
}

#' Deep Learning Background Removal (Transparent / RGBA)
#'
#' Removes the background of an image using pre-trained Deep Learning models,
#' producing a 4-channel **RGBA** image with a transparent background by default,
#' or replacing the background with a solid color.
#'
#' @inheritParams image_binary_dl
#' @param transparent Logical. If `TRUE` (default), returns a 4-channel RGBA `image` object
#'   where background pixels have an alpha channel value of `0` (100% transparent).
#' @param bg_color Optional background color to fill when `transparent = FALSE`. Accepts standard
#'   R color names (e.g., `"white"`, `"black"`) or hex strings (e.g., `"#FFFFFF"`). Defaults to `"black"`.
#' @return An `image` object (4-channel RGBA if `transparent = TRUE`, 3-channel RGB otherwise) or a list of images.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # Fast background removal with transparent alpha channel (CPU friendly)
#'   img_trans <- image_remove_bg_dl(img, model = "u2netp")
#'
#'   # Razor-sharp plant boundary cutout (leaf margins, serrations, petioles)
#'   img_ben <- image_remove_bg_dl(img, model = "ben2")
#'
#'   # Robust background removal under complex studio shadows and reflections
#'   img_rmbg <- image_remove_bg_dl(img, model = "rmbg-1.4")
#'
#'   # Replace background with solid white
#'   img_white <- image_remove_bg_dl(img, model = "u2netp", transparent = FALSE, bg_color = "white")
#'
#'   # Extract specific object via natural language prompt (Grounded-SAM)
#'   img_fruit <- image_remove_bg_dl(img, model = "grounded-sam", prompt = "red fruit")
#'
#'   # Extract instances matching a visual exemplar click (PerSAM)
#'   img_exemplar <- image_remove_bg_dl(img, model = "persam", prompt = c(120, 150))
#' }
image_remove_bg_dl <- function(img,
                               model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                               threshold = 0.5,
                               fill_hull = TRUE,
                               filter = 0,
                               erode = 0,
                               dilate = 0,
                               opening = 0,
                               closing = 0,
                               min_area = 0,
                               invert = FALSE,
                               pick_object = FALSE,
                               prompt = NULL,
                               exemplar = FALSE,
                               threads = 0,
                               engine = c("gpu", "cpu"),
                               device_id = -1,
                               verbose = TRUE,
                               transparent = TRUE,
                               bg_color = "black",
                               plot = TRUE,
                               dir = pliman_model_dir(),
                               ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_remove_bg_dl(x, model = model, threshold = threshold,
                         fill_hull = fill_hull, filter = filter,
                         erode = erode, dilate = dilate,
                         opening = opening, closing = closing,
                         min_area = min_area, invert = invert,
                         pick_object = pick_object, prompt = prompt,
                         exemplar = exemplar,
                         threads = threads,
                         engine = engine, device_id = device_id,
                         verbose = verbose,
                         transparent = transparent, bg_color = bg_color,
                         plot = FALSE, dir = dir)
    })
    return(res)
  }

  if (is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  mask <- .compute_dl_mask(mat = mat,
                           model = model,
                           threshold = threshold,
                           fill_hull = fill_hull,
                           filter = filter,
                           erode = erode,
                           dilate = dilate,
                           opening = opening,
                           closing = closing,
                           min_area = min_area,
                           invert = invert,
                           pick_object = pick_object,
                           prompt = prompt,
                           exemplar = exemplar,
                           threads = threads,
                           engine = engine,
                           device_id = device_id,
                           verbose = verbose,
                           dir = dir)

  if (isTRUE(transparent)) {
    # Generate 4-channel RGBA image with alpha layer
    is_raw_img <- is.raw(mat)
    nch <- if (length(dims) >= 3) dims[3] else 1

    if (nch >= 3) {
      R <- mat[, , 1]
      G <- mat[, , 2]
      B <- mat[, , 3]
    } else {
      R <- G <- B <- mat[, , 1]
    }

    if (is_raw_img) {
      alpha_layer <- as.raw(ifelse(mask, 255, 0))
      arr <- array(c(R, G, B, alpha_layer), dim = c(orig_w, orig_h, 4))
    } else {
      alpha_layer <- ifelse(mask, 1.0, 0.0)
      arr <- array(c(R, G, B, alpha_layer), dim = c(orig_w, orig_h, 4))
    }
    out <- as_image(arr, colormode = "Color")
  } else {
    out <- as_image(mat)
    # Fill background with bg_color
    col_rgb <- tryCatch(col2rgb(bg_color) / 255.0, error = function(e) c(0, 0, 0))
    if (length(dims) == 3) {
      for (k in seq_len(dims[3])) {
        fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
        out[, , k][!mask] <- fill_val
      }
    } else {
      fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
      out[!mask] <- fill_val
    }
  }

  if (isTRUE(plot)) {
    plot(out, ...)
  }

  invisible(out)
}

# ==============================================================================
# SECTION 3.5: GROUNDED-SAM (GROUNDING DINO + SAM 2.1) INSTANCE SEGMENTATION
# ==============================================================================

# Global cache for BERT vocabulary
.pliman_vocab_cache <- new.env(parent = emptyenv())

# Native WordPiece tokenizer for BERT / Grounding DINO (Zero external dependencies)
.bert_tokenize <- function(text, vocab_path) {
  if (!file.exists(vocab_path)) {
    stop("BERT vocab file not found: ", vocab_path)
  }
  if (!exists(vocab_path, envir = .pliman_vocab_cache, inherits = FALSE)) {
    vocab_lines <- readLines(vocab_path, warn = FALSE, encoding = "UTF-8")
    vocab_env <- new.env(hash = TRUE, parent = emptyenv(), size = length(vocab_lines))
    for (i in seq_along(vocab_lines)) {
      vocab_env[[vocab_lines[i]]] <- i - 1L
    }
    assign(vocab_path, vocab_env, envir = .pliman_vocab_cache)
  } else {
    vocab_env <- get(vocab_path, envir = .pliman_vocab_cache)
  }

  if (length(text) > 1) {
    text <- paste(trimws(text), collapse = " and ")
  }
  text_clean <- tolower(trimws(text))
  # Support commas, semicolons, and normalize delimiters into natural conjunctions
  text_clean <- gsub("[,;]+", " and ", text_clean)
  # Normalize internal periods between class names (e.g. "people . book ." -> "people and book .")
  text_clean <- gsub("\\s*\\.\\s*(?=[a-zA-Z0-9])", " and ", text_clean, perl = TRUE)
  # Collapse redundant 'and's or spaces
  text_clean <- gsub("\\b(and\\s*)+", "and ", text_clean)
  text_clean <- gsub("\\s+", " ", text_clean)
  text_clean <- trimws(text_clean)
  # Ensure clean trailing period expected by Grounding DINO
  text_clean <- gsub("\\s*\\.*$", " .", text_clean)
  if (!grepl("\\.$", text_clean)) {
    text_clean <- paste0(text_clean, " .")
  }

  # Insert spaces around punctuation
  text_spaced <- gsub("([[:punct:]])", " \\1 ", text_clean)
  words <- strsplit(text_spaced, "\\s+")[[1]]
  words <- words[nzchar(words)]

  token_ids <- 101L # [CLS]
  token_words <- "[CLS]"
  unk_id <- 100L

  for (w in words) {
    start <- 1L
    w_len <- nchar(w)
    cur_ids <- integer(0)
    is_unk <- FALSE

    while (start <= w_len) {
      end <- w_len
      found_sub <- FALSE
      while (start <= end) {
        sub <- substr(w, start, end)
        if (start > 1L) sub <- paste0("##", sub)
        if (exists(sub, envir = vocab_env, inherits = FALSE)) {
          cur_ids <- c(cur_ids, as.integer(get(sub, envir = vocab_env)))
          start <- end + 1L
          found_sub <- TRUE
          break
        }
        end <- end - 1L
      }
      if (!found_sub) {
        is_unk <- TRUE
        break
      }
    }

    if (is_unk) {
      token_ids <- c(token_ids, unk_id)
      token_words <- c(token_words, w)
    } else {
      token_ids <- c(token_ids, cur_ids)
      token_words <- c(token_words, rep(w, length(cur_ids)))
    }
  }

  token_ids <- c(token_ids, 102L) # [SEP]
  token_words <- c(token_words, "[SEP]")

  return(list(
    input_ids = token_ids,
    attention_mask = rep(1L, length(token_ids)),
    token_type_ids = rep(0L, length(token_ids)),
    token_words = token_words
  ))
}

# Helper to generate rounded rectangle polygon vertices
.rounded_rect_pts <- function(xleft, ybottom, xright, ytop, r = NULL, n = 8) {
  x0 <- min(xleft, xright)
  x1 <- max(xleft, xright)
  y0 <- min(ybottom, ytop)
  y1 <- max(ybottom, ytop)

  w <- x1 - x0
  h <- y1 - y0

  if (is.null(r) || is.na(r) || r <= 0) {
    r <- min(w * 0.25, h * 0.35, 4)
  } else {
    r <- min(r, w / 2, h / 2)
  }

  theta <- seq(0, pi / 2, length.out = n)
  arc_tr <- cbind((x1 - r) + r * cos(theta), (y0 + r) - r * sin(theta))
  arc_tl <- cbind((x0 + r) - r * sin(theta), (y0 + r) - r * cos(theta))
  arc_bl <- cbind((x0 + r) - r * cos(theta), (y1 - r) + r * sin(theta))
  arc_br <- cbind((x1 - r) + r * sin(theta), (y1 - r) + r * cos(theta))

  rbind(arc_tr, arc_tl, arc_bl, arc_br)
}

# Helper to draw YOLO-style bounding boxes with class label and confidence score
.plot_yolo_bboxes <- function(boxes,
                              palette_colors,
                              lwd = 2,
                              cex = 1.0,
                              pad = 1.0,
                              show_text = TRUE,
                              show_conf = TRUE,
                              show_class = TRUE,
                              show_id = FALSE,
                              badge = TRUE,
                              cex_scale = NULL) {
  if (is.null(boxes) || !is.data.frame(boxes) || nrow(boxes) == 0) return(invisible(NULL))

  if (!is.null(cex_scale)) {
    cex <- cex_scale
  }

  pad_mult_x <- if (length(pad) >= 1 && is.numeric(pad)) max(0, pad[1]) else 1.0
  pad_mult_y <- if (length(pad) >= 2 && is.numeric(pad)) max(0, pad[2]) else pad_mult_x

  u <- graphics::par("usr")
  img_h <- abs(u[3] - u[4])
  base_th <- abs(graphics::strheight("Ag", units = "user", cex = 1))
  target_th <- max(7, min(14, img_h * 0.012))
  base_cex <- if (base_th > 0) max(0.3, min(1.0, target_th / base_th)) else 0.55
  cex_val <- max(0.01, base_cex * cex)
  draw_text <- isTRUE(show_text) && cex_val > 0.005

  num_inst <- nrow(boxes)
  for (i in seq_len(num_inst)) {
    bx <- boxes[i, ]
    k_col <- if ("color" %in% names(bx) && !is.null(bx$color) && !is.na(bx$color) && nzchar(as.character(bx$color))) {
      as.character(bx$color)
    } else if (length(palette_colors) > 0) {
      palette_colors[((i - 1) %% length(palette_colors)) + 1]
    } else {
      "#00CC66"
    }
    x1 <- bx$xmin; y1 <- bx$ymin; x2 <- bx$xmax; y2 <- bx$ymax

    # Draw box outline (rotated polygon for OBB, rectangle for standard)
    if (all(c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4") %in% names(bx))) {
      graphics::polygon(x = c(bx$x1, bx$x2, bx$x3, bx$x4), y = c(bx$y1, bx$y2, bx$y3, bx$y4), border = k_col, lwd = lwd)
    } else {
      graphics::rect(xleft = x1, ybottom = y1, xright = x2, ytop = y2, border = k_col, lwd = lwd)
    }

    if (!draw_text) next

    # Label text formatting
    lbl_parts <- character(0)
    if (isTRUE(show_id) && !is.null(bx$id) && !is.na(bx$id)) {
      lbl_parts <- c(lbl_parts, paste0("#", bx$id))
    }
    if (isTRUE(show_class) && !is.null(bx$label) && !is.na(bx$label) && nzchar(as.character(bx$label))) {
      lbl_parts <- c(lbl_parts, as.character(bx$label))
    }
    if (isTRUE(show_conf) && !is.null(bx$score) && !is.na(bx$score)) {
      lbl_parts <- c(lbl_parts, sprintf("%.2f", bx$score))
    }
    lbl_txt <- paste(lbl_parts, collapse = " ")
    if (!nzchar(lbl_txt)) next

    raw_tw <- abs(graphics::strwidth(lbl_txt, units = "user", cex = cex_val, font = 2))
    th <- abs(graphics::strheight("Ag", units = "user", cex = cex_val, font = 2))
    if (th <= 0) th <- abs(graphics::strheight(lbl_txt, units = "user", cex = cex_val, font = 2))

    # Calculate character-proportional width to prevent oversized badges across different devices/DPI
    chars <- strsplit(lbl_txt, "")[[1]]
    char_weights <- vapply(chars, function(ch) {
      if (ch %in% c(" ", ".", ",", ":", ";", "!", "|", "'", "\"", "`", "(", ")", "[", "]", "{", "}", "/", "\\", "-")) {
        0.28
      } else if (ch %in% c("i", "l", "t", "r", "f", "j", "I", "1")) {
        0.35
      } else if (ch %in% c("#", "@", "%", "&", "m", "w", "M", "W")) {
        0.75
      } else if (grepl("[A-Z]", ch)) {
        0.60
      } else if (grepl("[0-9]", ch)) {
        0.54
      } else {
        0.50
      }
    }, numeric(1))
    approx_tw <- sum(char_weights) * th

    # Snug fit: adapt tightly to text width without oversized trailing background
    tw <- if (raw_tw > 0 && approx_tw > 0) min(raw_tw, approx_tw * 1.08) else max(raw_tw, approx_tw)
    pad_x <- th * 0.18 * pad_mult_x
    pad_y <- th * 0.12 * pad_mult_y

    badge_w <- tw + 5 * pad_x
    badge_h <- th + 2 * pad_y

    outside <- (y1 - badge_h) >= min(u[3], u[4])
    by0 <- if (outside) y1 - badge_h else y1
    by1 <- if (outside) y1 else y1 + badge_h
    bx0 <- x1
    bx1 <- x1 + badge_w

    # Ensure badge does not overflow right canvas boundary
    x_max_plot <- max(u[1], u[2])
    if (bx1 > x_max_plot) {
      shift <- bx1 - x_max_plot
      bx0 <- max(min(u[1], u[2]), bx0 - shift)
      bx1 <- bx0 + badge_w
    }

    # Text contrast color (black on bright badges, white on dark badges)
    rgb_vals <- tryCatch(grDevices::col2rgb(k_col)[, 1] / 255.0, error = function(e) c(0, 1, 0))
    lum <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
    txt_col <- if (lum > 0.55) "#111111" else "#ffffff"

    if (isTRUE(badge)) {
      # YOLO rounded badge tag
      r_val <- min(badge_h * 0.28, max(0, 3 * pad_mult_y))
      pts <- .rounded_rect_pts(bx0, by0, bx1, by1, r = r_val)
      graphics::polygon(pts[, 1], pts[, 2], col = k_col, border = NA)
      graphics::text(bx0 + pad_x, (by0 + by1) / 2, labels = lbl_txt, col = txt_col, font = 2, adj = c(0, 0.5), cex = cex_val)
    } else {
      # Draw text directly without badge background
      graphics::text(bx0, (by0 + by1) / 2, labels = lbl_txt, col = k_col, font = 2, adj = c(0, 0.5), cex = cex_val)
    }
  }
}

# Helper to draw YOLO pose keypoints and skeleton connections
.plot_yolo_keypoints <- function(keypoints_list,
                                 palette_colors,
                                 kpt_threshold = 0.3,
                                 lwd = 2,
                                 kpt_radius = 4) {
  if (length(keypoints_list) == 0) return(invisible(NULL))

  # Define colors for limbs
  limb_colors <- c(
    rep("#FF4B4B", 4), # facial: red
    rep("#FFA500", 4), # torso: orange
    rep("#00CC66", 2), # left arm: green
    rep("#0099FF", 2), # right arm: blue
    rep("#9933FF", 2), # left leg: purple
    rep("#FF00FF", 2)  # right leg: magenta
  )

  for (i in seq_along(keypoints_list)) {
    kp_df <- keypoints_list[[i]]
    if (is.null(kp_df) || nrow(kp_df) < 17) next
    k_col <- palette_colors[((i - 1) %% length(palette_colors)) + 1]

    # Draw skeleton limbs
    for (b in seq_len(nrow(.coco_skeleton_pairs))) {
      p1 <- .coco_skeleton_pairs[b, 1]
      p2 <- .coco_skeleton_pairs[b, 2]
      if (kp_df$conf[p1] >= kpt_threshold && kp_df$conf[p2] >= kpt_threshold) {
        graphics::lines(
          c(kp_df$x[p1], kp_df$x[p2]),
          c(kp_df$y[p1], kp_df$y[p2]),
          col = limb_colors[min(b, length(limb_colors))],
          lwd = lwd
        )
      }
    }

    # Draw keypoint dots
    valid_kpts <- which(kp_df$conf >= kpt_threshold)
    if (length(valid_kpts) > 0) {
      graphics::points(
        kp_df$x[valid_kpts],
        kp_df$y[valid_kpts],
        pch = 21,
        bg = k_col,
        col = "white",
        cex = max(0.6, kpt_radius / 3),
        lwd = 1.5
      )
    }
  }
  invisible(NULL)
}

# Helper to format human-readable instance detection summary
.format_detection_summary <- function(labels) {
  if (length(labels) == 0) return("0 objects")
  counts <- table(factor(labels, levels = unique(labels)))
  parts <- vapply(names(counts), function(cls) {
    cnt <- counts[[cls]]
    cls_lower <- tolower(cls)
    if (cnt == 1) {
      paste0(cnt, " ", cls)
    } else {
      # Inherent or already plural nouns
      if (cls_lower %in% c("shorts", "pants", "jeans", "trousers", "glasses", "sunglasses", "scissors", "children", "people")) {
        paste0(cnt, " ", cls)
      } else if (cls_lower == "child") {
        paste0(cnt, " children")
      } else if (cls_lower == "leaf") {
        paste0(cnt, " leaves")
      } else if (cls_lower == "person") {
        paste0(cnt, " people")
      } else if (grepl("[^aeiou]y$", cls_lower)) {
        # Words ending in consonant + y (e.g., butterfly -> butterflies)
        stem <- substr(cls, 1, nchar(cls) - 1)
        paste0(cnt, " ", stem, "ies")
      } else if (grepl("(ss|sh|ch|x|z)$", cls_lower)) {
        # Words ending in -ss, -sh, -ch, -x, -z (e.g. dress -> dresses, sunglass -> sunglasses)
        paste0(cnt, " ", cls, "es")
      } else if (grepl("s$", cls_lower)) {
        # Words ending in single s (already plural, e.g. cars, shoes, cows)
        paste0(cnt, " ", cls)
      } else {
        paste0(cnt, " ", cls, "s")
      }
    }
  }, character(1), USE.NAMES = FALSE)

  if (length(parts) == 1) {
    parts[1]
  } else if (length(parts) == 2) {
    paste(parts[1], "and", parts[2])
  } else {
    paste0(paste(parts[-length(parts)], collapse = ", "), ", and ", parts[length(parts)])
  }
}

# Helper to print highlighted detection summary block
.print_detection_summary <- function(summary_str, title = "Detection Summary") {
  if (is.null(summary_str)) return(invisible(NULL))

  if (is.data.frame(summary_str)) {
    if (nrow(summary_str) == 0) return(invisible(NULL))
    lbl_col <- intersect(c("class", "label"), tolower(names(summary_str)))
    cnt_col <- intersect(c("count", "freq"), tolower(names(summary_str)))
    if (length(lbl_col) > 0 && length(cnt_col) > 0) {
      lbls <- rep(as.character(summary_str[[lbl_col[1]]]), as.integer(summary_str[[cnt_col[1]]]))
      summary_text <- .format_detection_summary(lbls)
    } else {
      summary_text <- paste(apply(summary_str, 1, paste, collapse = ": "), collapse = ", ")
    }
  } else if (is.table(summary_str)) {
    if (length(summary_str) == 0) return(invisible(NULL))
    lbls <- rep(names(summary_str), as.integer(summary_str))
    summary_text <- .format_detection_summary(lbls)
  } else if (is.character(summary_str)) {
    if (length(summary_str) == 0 || !any(nzchar(summary_str))) return(invisible(NULL))
    if (length(summary_str) > 1) {
      summary_text <- .format_detection_summary(summary_str)
    } else {
      summary_text <- summary_str
    }
  } else {
    summary_text <- as.character(summary_str)
  }

  if (is.null(summary_text) || !nzchar(summary_text) || summary_text == "0 objects") return(invisible(NULL))

  cli::cli_alert_success("pliman detected {.bold {summary_text}} in the image.")
}



# Grounded-SAM runner: detects instances from text prompt and segments each with SAM 2.1
.run_grounded_sam <- function(mat,
                              prompt,
                              threshold = 0.5,
                              box_threshold = 0.25,
                              text_threshold = 0.25,
                              iou_threshold = 0.5,
                              threads = 0,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              fill_hull = TRUE,
                              filter = 0,
                              erode = 0,
                              dilate = 0,
                              opening = 0,
                              closing = 0,
                              min_area = 0,
                              invert = FALSE,
                              mask = TRUE,
                              verbose = TRUE,
                              dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  dir <- pliman_model_dir(dir)
  lib_file <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_file) || !file.exists(lib_file)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_file <- onnx_install(engine = engine)
  }

  # Ensure Grounding DINO and SAM 2.1 (if masking) are downloaded
  pliman_download_model("grounded-sam", dir = dir)
  if (isTRUE(mask)) {
    pliman_download_model("sam2.1", dir = dir)
  }

  dino_model <- file.path(dir, "groundingdino-tiny.onnx")
  vocab_file <- file.path(dir, "vocab.txt")
  sam_enc <- file.path(dir, "sam2.1.encoder.onnx")
  sam_dec <- file.path(dir, "sam2.1.decoder.onnx")

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # 1. Preprocess for Grounding DINO (800x800, ImageNet norm)
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Tokenizing prompt and preparing image...",
      msg_done = "Prompt tokenized"
    )
  }

  dino_tensor <- .preprocess_nchw(mat, target_size = 800L,
                                  mean = c(0.485, 0.456, 0.406),
                                  std = c(0.229, 0.224, 0.225),
                                  letterbox = FALSE)

  # 2. Tokenize prompt
  tok <- .bert_tokenize(prompt, vocab_file)

  # 3. Run Grounding DINO
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Detecting objects with Grounding DINO [{toupper(engine)}]...",
      msg_done = "Grounding DINO object detection complete"
    )
  }

  dino_res <- run_grounding_dino_cpp(
    pixel_values = as.numeric(dino_tensor),
    input_ids = as.integer(tok$input_ids),
    token_type_ids = as.integer(tok$token_type_ids),
    attention_mask = as.integer(tok$attention_mask),
    model_path = normalizePath(dino_model, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
    box_threshold = box_threshold,
    text_threshold = text_threshold,
    iou_threshold = iou_threshold,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  num_boxes <- length(dino_res$scores)
  if (num_boxes == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
      cli::cli_alert_warning("Grounding DINO detected 0 instances for prompt: {.val {prompt}}.")
    }
    return(list(
      boxes = data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                         xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                         label = character(0), stringsAsFactors = FALSE),
      labels = as_image(matrix(0L, nrow = orig_w, ncol = orig_h), colormode = "Grayscale", storage = "integer"),
      mask = as_image(matrix(FALSE, nrow = orig_w, ncol = orig_h))
    ))
  }

  # Scale normalized box coordinates to pixel coordinates on original image
  raw_boxes <- dino_res$boxes
  boxes_px <- matrix(0.0, nrow = num_boxes, ncol = 4)
  boxes_px[, 1] <- pmax(1.0, raw_boxes[, 1] * orig_w)
  boxes_px[, 2] <- pmax(1.0, raw_boxes[, 2] * orig_h)
  boxes_px[, 3] <- pmin(as.double(orig_w), raw_boxes[, 3] * orig_w)
  boxes_px[, 4] <- pmin(as.double(orig_h), raw_boxes[, 4] * orig_h)

  # Map token indices to labels (avoiding conjunctions or punctuation)
  labels_str <- character(num_boxes)
  stop_tokens <- c("[cls]", "[sep]", "and", ".", ",", ";")
  for (i in seq_len(num_boxes)) {
    t_idx <- dino_res$token_indices[i] + 1L
    w_cand <- if (t_idx >= 1 && t_idx <= length(tok$token_words)) tok$token_words[t_idx] else ""
    if (tolower(w_cand) %in% stop_tokens || !nzchar(w_cand)) {
      dists <- abs(seq_along(tok$token_words) - t_idx)
      valid_mask <- !(tolower(tok$token_words) %in% stop_tokens)
      if (any(valid_mask)) {
        dists[!valid_mask] <- Inf
        w_cand <- tok$token_words[which.min(dists)]
      } else {
        w_cand <- prompt
      }
    }
    labels_str[i] <- w_cand
  }

  boxes_df <- data.frame(
    id = seq_len(num_boxes),
    xmin = round(boxes_px[, 1], 1),
    ymin = round(boxes_px[, 2], 1),
    xmax = round(boxes_px[, 3], 1),
    ymax = round(boxes_px[, 4], 1),
    score = round(dino_res$scores, 3),
    label = labels_str,
    stringsAsFactors = FALSE
  )

  counts_df <- as.data.frame(table(factor(labels_str, levels = unique(labels_str))), stringsAsFactors = FALSE)
  colnames(counts_df) <- c("label", "count")
  summary_str <- .format_detection_summary(labels_str)

  # If no boxes found, return early
  if (nrow(boxes_df) == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = counts_df,
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = NULL
    ))
  }

  # 4. If mask = FALSE, skip SAM 2.1 completely (pure object detection mode)
  if (!isTRUE(mask)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = counts_df,
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = matrix(FALSE, nrow = orig_w, ncol = orig_h),
      features = data.frame()
    ))
  }

  # 4. Run SAM 2.1 on all detected boxes
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Segmenting {num_boxes} instance{?s} with SAM 2.1 [{toupper(engine)}]...",
      msg_done = "SAM 2.1 segmented {num_boxes} instance{?s}"
    )
  }

  sam_tensor <- .preprocess_nchw(mat, target_size = 1024L,
                                mean = c(0.485, 0.456, 0.406),
                                std = c(0.229, 0.224, 0.225),
                                letterbox = FALSE)

  raw_masks <- run_sam2_instances_cpp(
    tensor_vec = as.numeric(sam_tensor),
    boxes = boxes_px,
    orig_w = as.double(orig_w),
    orig_h = as.double(orig_h),
    encoder_path = normalizePath(sam_enc, winslash = "/", mustWork = FALSE),
    decoder_path = normalizePath(sam_dec, winslash = "/", mustWork = FALSE),
    lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  # 5. Post-process masks directly into a single multi-label matrix
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing multi-label instance masks...",
      msg_done = "Instance masks ready"
    )
  }

  combined_labels <- matrix(0L, nrow = orig_w, ncol = orig_h)

  for (k in seq_len(num_boxes)) {
    prob_crop <- raw_masks[[k]]
    # Resize continuous float logits to full resolution for silky-smooth sub-pixel boundaries
    mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
    m <- mask_prob >= threshold

    bx <- boxes_px[k, ]
    x1 <- max(1L, floor(bx[1]) - 2L)
    x2 <- min(orig_w, ceiling(bx[3]) + 2L)
    y1 <- max(1L, floor(bx[2]) - 2L)
    y2 <- min(orig_h, ceiling(bx[4]) + 2L)

    if (isTRUE(fill_hull)) {
      sub_m <- m[x1:x2, y1:y2]
      m[x1:x2, y1:y2] <- (fill_holes_cpp(sub_m) > 0)
    }
    if (filter > 0) {
      m_arr <- as.array(m)
      if (length(dim(m_arr)) == 2) dim(m_arr) <- c(orig_w, orig_h, 1L)
      f_res <- median_filter_binary_cpp(m_arr, orig_w, orig_h, 1L, as.integer(filter))
      m <- (f_res[, , 1] > 0)
    }
    if (erode > 0) {
      m <- (erode_cpp(m, raio = as.integer(erode)) > 0)
    }
    if (dilate > 0) {
      m <- (dilate_cpp(m, raio = as.integer(dilate)) > 0)
    }
    if (opening > 0) {
      m <- (erode_cpp(m, raio = as.integer(opening)) > 0)
      m <- (dilate_cpp(m, raio = as.integer(opening)) > 0)
    }
    if (closing > 0) {
      m <- (dilate_cpp(m, raio = as.integer(closing)) > 0)
      m <- (erode_cpp(m, raio = as.integer(closing)) > 0)
    }
    if (min_area > 0) {
      lbls <- bwlabel_cpp(m)
      if (max(lbls) > 0) {
        tbl <- table(lbls[lbls > 0])
        keep_ids <- as.integer(names(tbl[tbl >= min_area]))
        m <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
      }
    }

    combined_labels[m] <- as.integer(k)
  }

  labels_img <- as_image(combined_labels, colormode = "Grayscale", storage = "integer")

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }
  contornos <- contour(labels_img)
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  return(list(
    boxes = boxes_df,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (combined_labels > 0L),
    features = poly_measures(contornos)
  ))
}

# PerSAM (Personalize Segment Anything) runner: segments all instances similar to visual exemplar(s)
.run_persam <- function(mat,
                        exemplar_points = NULL,
                        precomputed_prototypes = NULL,
                        sim_threshold = 0.5,
                        min_dist = 16,
                        iou_threshold = 0.5,
                        max_objects = NULL,
                        feat_res = 256,
                        threshold = 0.5,
                        superres_map = FALSE,
                        threads = 0,
                        engine = c("gpu", "cpu"),
                        device_id = -1,
                        fill_hull = TRUE,
                        filter = 0,
                        erode = 0,
                        dilate = 0,
                        opening = 0,
                        closing = 0,
                        min_area = 0,
                        invert = FALSE,
                        mask = TRUE,
                        verbose = TRUE,
                        dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  if (is.character(feat_res)) {
    feat_res <- switch(tolower(feat_res[1]),
      "fast" = 64L,
      "vit" = 64L,
      "low" = 64L,
      "medium" = 256L,
      "high" = 1024L,
      "dense" = 1024L,
      "conv" = 1024L,
      "grain" = 1024L,
      "grains" = 1024L,
      "ultra" = 1024L,
      1024L
    )
  }
  feat_res <- as.integer(feat_res)
  if (is.na(feat_res) || feat_res <= 64L) {
    feat_res <- 64L
  } else if (feat_res <= 256L) {
    feat_res <- 256L
  } else {
    feat_res <- 1024L
  }
  if (engine == "gpu" && .Platform$OS.type == "windows") {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("SAM 2.1 Vision Transformer uses CPU execution (DirectML GPU has known driver limitations with ViT attention layers).")
    }
    use_gpu <- FALSE
    engine <- "cpu"
  } else {
    use_gpu <- (engine == "gpu")
  }
  c_max_objects <- if (is.null(max_objects) || is.infinite(max_objects) || max_objects <= 0) -1L else as.integer(max_objects)
  dir <- pliman_model_dir(dir)
  lib_file <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_file) || !file.exists(lib_file)) {
    cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading Microsoft library...")
    lib_file <- onnx_install(engine = engine)
  }

  # Ensure SAM 2.1 model is downloaded
  pliman_download_model("sam2.1", dir = dir)
  sam_enc <- file.path(dir, "sam2.1.encoder.onnx")
  sam_dec <- file.path(dir, "sam2.1.decoder.onnx")

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Interactive picking of exemplar points if not supplied and no precomputed prototypes
  if (!is.null(precomputed_prototypes)) {
    ex_x <- numeric(0)
    ex_y <- numeric(0)
  } else if (is.null(exemplar_points)) {
    cli::cli_alert_info("Click on 1 or more exemplar object(s) in the plot window. Press <Esc> or right-click when finished.")
    plot(as_image(mat))
    pts <- tryCatch(graphics::locator(n = 512, type = "p", col = "cyan", pch = 19), error = function(e) NULL)
    if (!is.null(pts) && length(pts$x) > 0) {
      ex_x <- pts$x
      ex_y <- pts$y
    } else {
      cli::cli_abort("No point selected.")
    }
  } else if (is.matrix(exemplar_points) || is.data.frame(exemplar_points)) {
    ex_x <- as.numeric(exemplar_points[, 1])
    ex_y <- as.numeric(exemplar_points[, 2])
    if (length(ex_x) == 0) {
      cli::cli_abort("No exemplar point provided.")
    }
  } else if (is.numeric(exemplar_points)) {
    if (length(exemplar_points) == 2) {
      ex_x <- exemplar_points[1]
      ex_y <- exemplar_points[2]
    } else if (length(exemplar_points) > 0 && length(exemplar_points) %% 2 == 0) {
      half <- length(exemplar_points) / 2
      ex_x <- exemplar_points[seq_len(half)]
      ex_y <- exemplar_points[(half + 1L):length(exemplar_points)]
    } else {
      cli::cli_abort("Invalid {.arg exemplar_points}. Expected numeric coordinates or a matrix/data.frame.")
    }
  } else {
    cli::cli_abort("Invalid {.arg exemplar_points}. Expected numeric coordinates or a matrix/data.frame.")
  }

  # 1. Preprocess for SAM 2.1 (1024x1024, ImageNet norm)
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Computing visual embeddings with SAM 2.1 [{toupper(engine)}]...",
      msg_done = "SAM 2.1 visual embeddings ready"
    )
  }

  sam_tensor <- .preprocess_nchw(mat, target_size = 1024L,
                                mean = c(0.485, 0.456, 0.406),
                                std = c(0.229, 0.224, 0.225),
                                letterbox = FALSE)

  # 2. Run PerSAM C++ routine
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Calculating PerSAM feature similarity map & segmenting matches...",
      msg_done = "PerSAM exemplar matching complete"
    )
  }

  persam_res <- tryCatch({
    run_sam2_persam_cpp(
      tensor_vec = as.numeric(sam_tensor),
      exemplar_x = as.numeric(ex_x),
      exemplar_y = as.numeric(ex_y),
      orig_w = as.double(orig_w),
      orig_h = as.double(orig_h),
      encoder_path = normalizePath(sam_enc, winslash = "/", mustWork = FALSE),
      decoder_path = normalizePath(sam_dec, winslash = "/", mustWork = FALSE),
      lib_path = normalizePath(lib_file, winslash = "/", mustWork = FALSE),
      sim_threshold = as.double(sim_threshold),
      min_dist = as.double(min_dist),
      iou_threshold = as.double(iou_threshold),
      max_objects = c_max_objects,
      feat_res = as.integer(feat_res),
      num_threads = as.integer(threads),
      use_gpu = use_gpu,
      device_id = as.integer(device_id),
      precomputed_prototypes = precomputed_prototypes
    )
  }, error = function(e) {
    if (isTRUE(use_gpu)) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning("GPU execution failed ({e$message}). Retrying on CPU...")
      }
      run_sam2_persam_cpp(
        tensor_vec = as.numeric(sam_tensor),
        exemplar_x = as.numeric(ex_x),
        exemplar_y = as.numeric(ex_y),
        orig_w = as.double(orig_w),
        orig_h = as.double(orig_h),
        encoder_path = normalizePath(sam_enc, winslash = "/", mustWork = FALSE),
        decoder_path = normalizePath(sam_dec, winslash = "/", mustWork = FALSE),
        lib_path = normalizePath(pliman_onnx_lib_path(engine = "cpu"), winslash = "/", mustWork = FALSE),
        sim_threshold = as.double(sim_threshold),
        min_dist = as.double(min_dist),
        iou_threshold = as.double(iou_threshold),
        max_objects = c_max_objects,
        feat_res = as.integer(feat_res),
        num_threads = as.integer(threads),
        use_gpu = FALSE,
        device_id = -1L,
        precomputed_prototypes = precomputed_prototypes
      )
    } else {
      stop(e)
    }
  })

  # Process similarity map: normalize to [0, 1] as an image object (resolution governed by feat_res)
  sim_raw <- persam_res$similarity_map
  sim_min <- min(sim_raw)
  sim_max <- max(sim_raw)
  sim_norm <- sim_raw / sim_max
  sim_img <- as_image(sim_norm, colormode = "Grayscale", storage = "double")
  attr(sim_img, "raw_matrix") <- sim_norm
  attr(sim_img, "raw_range") <- c(sim_min, sim_max)

  if (isTRUE(superres_map)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Applying Real-ESRGAN neural super-resolution to similarity map...")
    }
    sim_img <- tryCatch({
      sr <- image_superres_dl(
        sim_img,
        model = "realesrgan-compact",
        engine = engine,
        device_id = device_id,
        threads = threads,
        verbose = FALSE,
        plot = FALSE,
        dir = dir
      )
      sr_d <- as.numeric(image_data(sr))
      dim(sr_d) <- dim(image_data(sr))
      sr_gray <- if (length(dim(sr_d)) >= 3 && dim(sr_d)[3] >= 3) {
        (sr_d[, , 1] + sr_d[, , 2] + sr_d[, , 3]) / 3.0
      } else {
        sr_d[, , 1]
      }
      # Re-normalize cleanly to [0, 1]
      s_min <- min(sr_gray)
      s_max <- max(sr_gray)
      if (s_max > s_min) {
        sr_gray <- (sr_gray - s_min) / (s_max - s_min)
      }
      sr_obj <- as_image(sr_gray, colormode = "Grayscale", storage = "double")
      attr(sr_obj, "raw_matrix") <- sr_gray
      attr(sr_obj, "raw_range") <- c(s_min, s_max)
      sr_obj
    }, error = function(e) {
      if (isTRUE(verbose)) {
        cli::cli_alert_warning("Super-resolution failed: {e$message}. Using native similarity map.")
      }
      sim_img
    })
  }

  num_inst <- length(persam_res$scores)
  if (num_inst == 0) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
      cli::cli_alert_warning("PerSAM found 0 instances matching the exemplar (sim_threshold = {sim_threshold}).")
    }
    return(list(
      boxes = data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                         xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                         label = character(0), stringsAsFactors = FALSE),
      counts = data.frame(label = "exemplar", count = 0L, stringsAsFactors = FALSE),
      summary = "0 objects similar to the exemplar",
      contours = list(),
      labels = as_image(matrix(0L, nrow = orig_w, ncol = orig_h), colormode = "Grayscale", storage = "integer"),
      mask = as_image(matrix(FALSE, nrow = orig_w, ncol = orig_h)),
      similarity_map = sim_img,
      prototypes = persam_res$prototypes
    ))
  }

  raw_boxes <- persam_res$boxes
  boxes_df <- data.frame(
    id = seq_len(num_inst),
    xmin = round(raw_boxes[, 1], 1),
    ymin = round(raw_boxes[, 2], 1),
    xmax = round(raw_boxes[, 3], 1),
    ymax = round(raw_boxes[, 4], 1),
    score = round(persam_res$scores, 3),
    label = "exemplar",
    stringsAsFactors = FALSE
  )

  summary_str <- paste0(num_inst, if (num_inst == 1) " object similar to the exemplar" else " objects similar to the exemplar")

  # If mask = FALSE, skip mask post-processing and contour extraction
  if (!isTRUE(mask)) {
    if (isTRUE(verbose)) {
      cli::cli_progress_done()
    }
    return(list(
      boxes = boxes_df,
      counts = data.frame(label = "exemplar", count = num_inst, stringsAsFactors = FALSE),
      summary = summary_str,
      contours = list(),
      labels = NULL,
      mask = matrix(FALSE, nrow = orig_w, ncol = orig_h),
      similarity_map = sim_img,
      features = data.frame()
    ))
  }

  # 3. Post-process masks into multi-label matrix
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Post-processing multi-label instance masks...",
      msg_done = "Instance masks ready"
    )
  }

  if (!is.null(persam_res$labels) &&
      filter == 0 && erode == 0 && dilate == 0 && opening == 0 && closing == 0 && min_area == 0) {
    combined_labels <- persam_res$labels
  } else if (!is.null(persam_res$labels)) {
    combined_labels <- persam_res$labels
    if (min_area > 0) {
      tbl <- table(combined_labels[combined_labels > 0L])
      drop_ids <- as.integer(names(tbl[tbl < min_area]))
      if (length(drop_ids) > 0) {
        combined_labels[combined_labels %in% drop_ids] <- 0L
      }
    }
    if (filter > 0 || erode > 0 || dilate > 0 || opening > 0 || closing > 0) {
      for (k in seq_len(num_inst)) {
        bx <- raw_boxes[k, ]
        x1 <- max(1L, floor(bx[1]) - 2L)
        x2 <- min(orig_w, ceiling(bx[3]) + 2L)
        y1 <- max(1L, floor(bx[2]) - 2L)
        y2 <- min(orig_h, ceiling(bx[4]) + 2L)
        sub_m <- (combined_labels[x1:x2, y1:y2] == as.integer(k))
        if (!any(sub_m)) next

        if (filter > 0) {
          m_arr <- as.array(sub_m)
          if (length(dim(m_arr)) == 2) dim(m_arr) <- c(nrow(sub_m), ncol(sub_m), 1L)
          f_res <- median_filter_binary_cpp(m_arr, nrow(sub_m), ncol(sub_m), 1L, as.integer(filter))
          sub_m <- (f_res[, , 1] > 0)
        }
        if (erode > 0) sub_m <- (erode_cpp(sub_m, raio = as.integer(erode)) > 0)
        if (dilate > 0) sub_m <- (dilate_cpp(sub_m, raio = as.integer(dilate)) > 0)
        if (opening > 0) {
          sub_m <- (erode_cpp(sub_m, raio = as.integer(opening)) > 0)
          sub_m <- (dilate_cpp(sub_m, raio = as.integer(opening)) > 0)
        }
        if (closing > 0) {
          sub_m <- (dilate_cpp(sub_m, raio = as.integer(closing)) > 0)
          sub_m <- (erode_cpp(sub_m, raio = as.integer(closing)) > 0)
        }
        patch <- combined_labels[x1:x2, y1:y2]
        patch[patch == as.integer(k)] <- 0L
        patch[sub_m] <- as.integer(k)
        combined_labels[x1:x2, y1:y2] <- patch
      }
    }
  } else {
    combined_labels <- matrix(0L, nrow = orig_w, ncol = orig_h)
    for (k in seq_len(num_inst)) {
      prob_crop <- persam_res$masks[[k]]
      mask_prob <- .bilinear_resize_2d(prob_crop, orig_w, orig_h)
      m <- mask_prob >= threshold

      bx <- raw_boxes[k, ]
      x1 <- max(1L, floor(bx[1]) - 2L)
      x2 <- min(orig_w, ceiling(bx[3]) + 2L)
      y1 <- max(1L, floor(bx[2]) - 2L)
      y2 <- min(orig_h, ceiling(bx[4]) + 2L)

      if (isTRUE(fill_hull)) {
        sub_m <- m[x1:x2, y1:y2]
        m[x1:x2, y1:y2] <- (fill_holes_cpp(sub_m) > 0)
      }
      if (filter > 0) {
        m_arr <- as.array(m)
        if (length(dim(m_arr)) == 2) dim(m_arr) <- c(orig_w, orig_h, 1L)
        f_res <- median_filter_binary_cpp(m_arr, orig_w, orig_h, 1L, as.integer(filter))
        m <- (f_res[, , 1] > 0)
      }
      if (erode > 0) {
        m <- (erode_cpp(m, raio = as.integer(erode)) > 0)
      }
      if (dilate > 0) {
        m <- (dilate_cpp(m, raio = as.integer(dilate)) > 0)
      }
      if (opening > 0) {
        m <- (erode_cpp(m, raio = as.integer(opening)) > 0)
        m <- (dilate_cpp(m, raio = as.integer(opening)) > 0)
      }
      if (closing > 0) {
        m <- (dilate_cpp(m, raio = as.integer(closing)) > 0)
        m <- (erode_cpp(m, raio = as.integer(closing)) > 0)
      }
      if (min_area > 0) {
        lbls <- bwlabel_cpp(m)
        if (max(lbls) > 0) {
          tbl <- table(lbls[lbls > 0])
          keep_ids <- as.integer(names(tbl[tbl >= min_area]))
          m <- matrix(lbls %in% keep_ids, nrow = orig_w, ncol = orig_h)
        }
      }

      combined_labels[m] <- as.integer(k)
    }
  }

  labels_img <- as_image(combined_labels, colormode = "Grayscale", storage = "integer")

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }

  contornos_raw <- contour(labels_img)
  valid_ids <- which(vapply(contornos_raw, function(p) !is.null(p) && is.matrix(p) && nrow(p) >= 3L, logical(1)))
  if (length(valid_ids) > 0L) {
    contornos <- contornos_raw[valid_ids]
    if (nrow(boxes_df) >= max(valid_ids)) {
      boxes_df <- boxes_df[valid_ids, , drop = FALSE]
      boxes_df$id <- seq_len(nrow(boxes_df))
    }
  } else {
    contornos <- list()
    boxes_df <- boxes_df[0, , drop = FALSE]
  }
  num_inst <- nrow(boxes_df)
  summary_str <- paste0(num_inst, if (num_inst == 1) " object similar to the exemplar" else " objects similar to the exemplar")
  return(list(
    boxes = boxes_df,
    counts = data.frame(label = "exemplar", count = num_inst, stringsAsFactors = FALSE),
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (combined_labels > 0L),
    similarity_map = sim_img,
    prototypes = persam_res$prototypes,
    features = tryCatch(poly_measures(contornos), error = function(e) data.frame())
  ))
}

# Internal helper to preprocess an image for YOLO models with letterboxing to 640x640
.preprocess_yolo <- function(mat, target_size = 640L) {
  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  dim(R) <- dim(G) <- dim(B) <- c(orig_w, orig_h)

  gain <- min(target_size / orig_w, target_size / orig_h)
  new_w <- max(1L, round(orig_w * gain))
  new_h <- max(1L, round(orig_h * gain))

  pad_x <- (target_size - new_w) / 2.0
  pad_y <- (target_size - new_h) / 2.0
  x1 <- floor(pad_x) + 1L
  y1 <- floor(pad_y) + 1L
  x2 <- x1 + new_w - 1L
  y2 <- y1 + new_h - 1L

  R_res <- .bilinear_resize_2d(R, new_w, new_h)
  G_res <- .bilinear_resize_2d(G, new_w, new_h)
  B_res <- .bilinear_resize_2d(B, new_w, new_h)

  fill_val <- 114.0 / 255.0
  canvas_R <- matrix(fill_val, nrow = target_size, ncol = target_size)
  canvas_G <- matrix(fill_val, nrow = target_size, ncol = target_size)
  canvas_B <- matrix(fill_val, nrow = target_size, ncol = target_size)

  canvas_R[x1:x2, y1:y2] <- R_res
  canvas_G[x1:x2, y1:y2] <- G_res
  canvas_B[x1:x2, y1:y2] <- B_res

  vec <- c(as.numeric(canvas_R), as.numeric(canvas_G), as.numeric(canvas_B))
  attr(vec, "gain") <- gain
  attr(vec, "pad_x") <- pad_x
  attr(vec, "pad_y") <- pad_y
  vec
}

# Internal helper to run YOLO detection, segmentation, or pose estimation
.run_yolo <- function(mat,
                      model = "yolo26n",
                      conf_threshold = 0.25,
                      iou_threshold = 0.45,
                      labels = NULL,
                      threads = 0,
                      engine = c("gpu", "cpu"),
                      device_id = -1,
                      fill_hull = TRUE,
                      filter = 0,
                      erode = 0,
                      dilate = 0,
                      opening = 0,
                      closing = 0,
                      min_area = 0,
                      need_masks = TRUE,
                      return_features = FALSE,
                      txt_feats = NULL,
                      dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (file.exists(model)) {
    model_file <- normalizePath(model, winslash = "/")
  } else if (file.exists(file.path(dir, model))) {
    model_file <- normalizePath(file.path(dir, model), winslash = "/")
  } else if (file.exists(file.path(dir, paste0(model, ".onnx")))) {
    model_file <- normalizePath(file.path(dir, paste0(model, ".onnx")), winslash = "/")
  } else if (file.exists(file.path("D:/Desktop/models", model))) {
    model_file <- normalizePath(file.path("D:/Desktop/models", model), winslash = "/")
  } else if (file.exists(file.path("D:/Desktop/models", paste0(model, ".onnx")))) {
    model_file <- normalizePath(file.path("D:/Desktop/models", paste0(model, ".onnx")), winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model, dir = dir)
  }

  # Detect model type from name
  model_lower <- tolower(model)
  is_obb_model  <- grepl("obb", model_lower)
  is_nas_model  <- grepl("nas", model_lower)

  tensor <- tryCatch(
    preprocess_yolo_cpp(mat, target_size = 640L),
    error = function(e) .preprocess_yolo(mat, target_size = 640L)
  )
  orig_w <- attr(tensor, "orig_w")
  if (is.null(orig_w)) orig_w <- dim(mat)[1]
  orig_h <- attr(tensor, "orig_h")
  if (is.null(orig_h)) orig_h <- dim(mat)[2]

  raw_res <- run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = orig_w,
    orig_h = orig_h,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id),
    need_masks = isTRUE(need_masks),
    return_features = isTRUE(return_features),
    txt_feats = txt_feats
  )

  if (!is.null(labels)) {
    class_names <- labels
  } else {
    model_classes <- .get_onnx_classes(model_file)
    if (!is.null(model_classes) && length(model_classes) > 0L) {
      class_names <- model_classes
    } else if (is_obb_model) {
      class_names <- .dota_classes
    } else if (is_nas_model) {
      class_names <- .coco_classes
    } else {
      class_names <- .coco_classes
    }
  }

  num_boxes <- nrow(raw_res$boxes)
  if (num_boxes > 0) {
    c_ids <- raw_res$class_ids
    lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))

    # Build base box dataframe
    df_boxes <- data.frame(
      id = seq_len(num_boxes),
      xmin = raw_res$boxes[, 1],
      ymin = raw_res$boxes[, 2],
      xmax = raw_res$boxes[, 3],
      ymax = raw_res$boxes[, 4],
      label = lbls,
      score = round(raw_res$scores, 4),
      class_id = c_ids,
      stringsAsFactors = FALSE
    )

    # Append OBB corner coordinates for oriented bounding boxes
    has_obb_corners <- is_obb_model && !is.null(raw_res$corners) &&
                       is.matrix(raw_res$corners) && ncol(raw_res$corners) == 8 &&
                       nrow(raw_res$corners) == num_boxes
    if (has_obb_corners) {
      df_boxes$x1 <- raw_res$corners[, 1]
      df_boxes$y1 <- raw_res$corners[, 2]
      df_boxes$x2 <- raw_res$corners[, 3]
      df_boxes$y2 <- raw_res$corners[, 4]
      df_boxes$x3 <- raw_res$corners[, 5]
      df_boxes$y3 <- raw_res$corners[, 6]
      df_boxes$x4 <- raw_res$corners[, 7]
      df_boxes$y4 <- raw_res$corners[, 8]
      df_boxes$angle <- raw_res$angles
    }
    counts <- table(factor(df_boxes$label, levels = unique(df_boxes$label)))
    counts_df <- data.frame(
      class = names(counts),
      count = as.integer(counts),
      stringsAsFactors = FALSE
    )
    summary_str <- .format_detection_summary(df_boxes$label)
  } else {
    df_boxes <- data.frame(
      id = integer(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      label = character(0),
      score = numeric(0),
      class_id = integer(0),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(class = character(0), count = integer(0), stringsAsFactors = FALSE)
    summary_str <- "0 objects"
  }

  lbl_mat <- raw_res$labels
  if (num_boxes > 0 && !is.null(lbl_mat) && is.matrix(lbl_mat) && any(lbl_mat > 0L)) {
    for (k in seq_len(num_boxes)) {
      bx <- df_boxes[k, ]
      x1 <- max(1L, floor(bx$xmin) - 2L)
      x2 <- min(orig_w, ceiling(bx$xmax) + 2L)
      y1 <- max(1L, floor(bx$ymin) - 2L)
      y2 <- min(orig_h, ceiling(bx$ymax) + 2L)
      sub_m <- (lbl_mat[x1:x2, y1:y2] == k)
      if (!any(sub_m)) next

      if (isTRUE(fill_hull)) {
        sub_m <- (fill_holes_cpp(sub_m) > 0)
      }
      if (filter > 0) {
        m_arr <- as.array(sub_m)
        if (length(dim(m_arr)) == 2) dim(m_arr) <- c(nrow(sub_m), ncol(sub_m), 1L)
        sub_m <- (median_filter_binary_cpp(m_arr, nrow(sub_m), ncol(sub_m), 1L, as.integer(filter))[, , 1] > 0)
      }
      if (erode > 0) {
        sub_m <- (erode_cpp(sub_m, raio = as.integer(erode)) > 0)
      }
      if (dilate > 0) {
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(dilate)) > 0)
      }
      if (opening > 0) {
        sub_m <- (erode_cpp(sub_m, raio = as.integer(opening)) > 0)
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(opening)) > 0)
      }
      if (closing > 0) {
        sub_m <- (dilate_cpp(sub_m, raio = as.integer(closing)) > 0)
        sub_m <- (erode_cpp(sub_m, raio = as.integer(closing)) > 0)
      }
      if (min_area > 0) {
        lbls <- bwlabel_cpp(sub_m)
        if (max(lbls) > 0) {
          tbl <- table(lbls[lbls > 0])
          keep_ids <- as.integer(names(tbl[tbl >= min_area]))
          sub_m <- matrix(lbls %in% keep_ids, nrow = nrow(sub_m), ncol = ncol(sub_m))
        }
      }
      lbl_mat[x1:x2, y1:y2][lbl_mat[x1:x2, y1:y2] == k] <- 0L
      lbl_mat[x1:x2, y1:y2][sub_m] <- as.integer(k)
    }
  }

  labels_img <- if (!is.null(lbl_mat) && is.matrix(lbl_mat)) {
    as_image(lbl_mat, colormode = "Grayscale", storage = "integer")
  } else NULL

  contornos <- if (num_boxes > 0 && !is.null(lbl_mat) && is.matrix(lbl_mat) && any(lbl_mat > 0L)) contour(labels_img) else list()
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  feats <- if (length(contornos) > 0) {
    tryCatch(poly_measures(contornos), error = function(e) data.frame())
  } else {
    data.frame()
  }

  has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && num_boxes > 0
  kpts_list <- if (has_kpts) {
    lapply(seq_len(num_boxes), function(k) {
      vals <- raw_res$keypoints[k, ]
      data.frame(
        id = k,
        keypoint = .coco_keypoints,
        x = vals[seq(1, 51, by = 3)],
        y = vals[seq(2, 51, by = 3)],
        conf = round(vals[seq(3, 51, by = 3)], 4),
        stringsAsFactors = FALSE
      )
    })
  } else {
    list()
  }

  heatmap_img <- if (!is.null(raw_res$heatmap)) {
    as_image(raw_res$heatmap, colormode = "Grayscale", storage = "double")
  } else NULL

  energy_img <- if (!is.null(raw_res$feature_energy)) {
    as_image(raw_res$feature_energy, colormode = "Grayscale", storage = "double")
  } else NULL

  res_out <- list(
    boxes = df_boxes,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = if (!is.null(lbl_mat) && is.matrix(lbl_mat)) (lbl_mat > 0L) else NULL,
    features = feats,
    keypoints = kpts_list,
    raw_keypoints = if (has_kpts) raw_res$keypoints else NULL
  )
  if (isTRUE(return_features)) {
    res_out$heatmap <- heatmap_img
    res_out$feature_energy <- energy_img
    res_out$feature_map <- raw_res$feature_map
    res_out$embeddings <- raw_res$embeddings
    if (nrow(df_boxes) > 0 && !is.null(raw_res$embeddings) && nrow(raw_res$embeddings) == nrow(df_boxes)) {
      attr(res_out$boxes, "embeddings") <- raw_res$embeddings
    }
  }
  res_out
}

# Internal helper to run StarDist polygon detection
.run_stardist <- function(mat,
                          model = "stardist",
                          prob_threshold = 0.5,
                          nms_threshold = 0.3,
                          threads = 0,
                          engine = c("gpu", "cpu"),
                          device_id = -1,
                          dir = pliman_model_dir()) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (file.exists(model)) {
    model_file <- normalizePath(model, winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model, dir = dir)
  }

  dims <- dim(mat)
  orig_w <- dims[1]
  orig_h <- dims[2]

  target_size <- 256L
  tensor <- .preprocess_nchw(
    mat,
    target_size = target_size,
    mean = c(0, 0, 0),
    std = c(1, 1, 1),
    letterbox = FALSE
  )

  raw_res <- run_stardist_cpp(
    tensor_vec = tensor,
    in_w = target_size,
    in_h = target_size,
    orig_w = orig_w,
    orig_h = orig_h,
    prob_threshold = prob_threshold,
    nms_threshold = nms_threshold,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  num_objs <- nrow(raw_res$boxes)
  if (num_objs > 0) {
    df_boxes <- data.frame(
      id = seq_len(num_objs),
      xmin = raw_res$boxes[, 1],
      ymin = raw_res$boxes[, 2],
      xmax = raw_res$boxes[, 3],
      ymax = raw_res$boxes[, 4],
      center_x = raw_res$centers_x,
      center_y = raw_res$centers_y,
      label = rep("object", num_objs),
      score = round(raw_res$scores, 4),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(
      class = "object",
      count = num_objs,
      stringsAsFactors = FALSE
    )
    summary_str <- .format_detection_summary(df_boxes$label)
  } else {
    df_boxes <- data.frame(
      id = integer(0),
      xmin = numeric(0),
      ymin = numeric(0),
      xmax = numeric(0),
      ymax = numeric(0),
      center_x = numeric(0),
      center_y = numeric(0),
      label = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    )
    counts_df <- data.frame(class = character(0), count = integer(0), stringsAsFactors = FALSE)
    summary_str <- "0 objects"
  }

  lbl_mat <- raw_res$labels
  labels_img <- as_image(lbl_mat, colormode = "Grayscale", storage = "integer")
  contornos <- if (num_objs > 0 && any(lbl_mat > 0L)) contour(labels_img) else list()
  if (is.list(contornos)) {
    contornos <- contornos[!vapply(contornos, is.null, logical(1))]
  }
  feats <- if (length(contornos) > 0) {
    tryCatch(poly_measures(contornos), error = function(e) data.frame())
  } else {
    data.frame()
  }

  list(
    boxes = df_boxes,
    counts = counts_df,
    summary = summary_str,
    contours = contornos,
    labels = labels_img,
    mask = (lbl_mat > 0L),
    features = feats,
    polygons_x = raw_res$polygons_x,
    polygons_y = raw_res$polygons_y,
    centers_x = raw_res$centers_x,
    centers_y = raw_res$centers_y
  )
}

#' @keywords internal
.conts_to_boxes <- function(conts, max_w, max_h) {
  if (is.null(conts) || length(conts) == 0L) {
    return(data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                      xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }
  valid <- lapply(seq_along(conts), function(i) {
    p <- conts[[i]]
    if (!is.matrix(p) || nrow(p) == 0L) return(NULL)
    data.frame(
      id = i,
      xmin = max(1L, round(min(p[, 1]))),
      ymin = max(1L, round(min(p[, 2]))),
      xmax = min(max_w, round(max(p[, 1]))),
      ymax = min(max_h, round(max(p[, 2]))),
      score = 1.0,
      label = "object",
      stringsAsFactors = FALSE
    )
  })
  valid <- valid[!vapply(valid, is.null, logical(1))]
  if (length(valid) == 0L) {
    return(data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                      xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, valid)
}

#' @keywords internal
.extract_object_features <- function(mat, boxes, feature_model = c("clip-vit-b32", "dinov2"),
                                     threads = 0, engine = "cpu", device_id = -1, dir = pliman_model_dir()) {
  feature_model <- match.arg(feature_model)
  dim_out <- if (feature_model == "dinov2") 384L else 512L
  if (is.null(boxes) || nrow(boxes) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = dim_out))
  }

  dims <- dim(mat)
  max_x <- dims[1]
  max_y <- dims[2]
  n_boxes <- nrow(boxes)

  feats_list <- lapply(seq_len(n_boxes), function(i) {
    x1 <- max(1L, min(max_x, round(boxes$xmin[i])))
    x2 <- max(1L, min(max_x, round(boxes$xmax[i])))
    y1 <- max(1L, min(max_y, round(boxes$ymin[i])))
    y2 <- max(1L, min(max_y, round(boxes$ymax[i])))
    if (x2 < x1) { tmp <- x1; x1 <- x2; x2 <- tmp }
    if (y2 < y1) { tmp <- y1; y1 <- y2; y2 <- tmp }
    if ((x2 - x1) < 2L || (y2 - y1) < 2L) {
      return(rep(0.0, dim_out))
    }

    crop <- mat[x1:x2, y1:y2, , drop = FALSE]
    if (feature_model == "dinov2") {
      res_dino <- image_features_dl(crop, model = "dinov2", return_pca = FALSE,
                                    threads = threads, engine = engine, device_id = device_id,
                                    dir = dir, plot = FALSE)
      as.numeric(res_dino$cls_token)
    } else {
      as.numeric(image_embed_dl(crop, model = "clip-vit-b32", threads = threads,
                                engine = engine, device_id = device_id, dir = dir))
    }
  })

  feat_mat <- do.call(rbind, feats_list)
  row_ids <- if (!is.null(boxes$id)) paste0("obj_", boxes$id) else paste0("obj_", seq_len(n_boxes))
  rownames(feat_mat) <- row_ids
  feat_mat
}

#' @keywords internal
.attach_features <- function(obj,
                             mat,
                             boxes = NULL,
                             conts = NULL,
                             mask_mat = NULL,
                             return_features = FALSE,
                             feature_model = c("yolo", "clip-vit-b32", "dinov2"),
                             plot_features = FALSE,
                             threads = 0,
                             engine = "cpu",
                             device_id = -1,
                             dir = pliman_model_dir(),
                             verbose = TRUE,
                             yolo_features = NULL) {
  if (!isTRUE(return_features)) {
    return(obj)
  }

  feature_model <- match.arg(feature_model)

  # 1. Native YOLO Feature Maps, Heatmap, and Instance Embeddings
  if (feature_model == "yolo" && !is.null(yolo_features)) {
    sp_map <- yolo_features$heatmap
    feat_mat <- yolo_features$embeddings

    if (is.list(obj) && !inherits(obj, c("Image", "image", "data.frame")) && !is.null(obj$boxes)) {
      obj$heatmap <- yolo_features$heatmap
      obj$feature_energy <- yolo_features$feature_energy
      obj$feature_map <- yolo_features$feature_map
      obj$embeddings <- yolo_features$embeddings
      if (nrow(obj$boxes) > 0 && !is.null(feat_mat)) {
        attr(obj$boxes, "embeddings") <- feat_mat
      }
    }
    attr(obj, "heatmap") <- yolo_features$heatmap
    attr(obj, "feature_energy") <- yolo_features$feature_energy
    attr(obj, "feature_map") <- yolo_features$feature_map
    attr(obj, "embeddings") <- yolo_features$embeddings

    if (isTRUE(plot_features) && !is.null(sp_map)) {
      plot_yolo_heatmap(mat, sp_map)
    }
    return(obj)
  }

  # If "yolo" was chosen but yolo_features is NULL (e.g. non-YOLO model), fallback to dinov2
  if (feature_model == "yolo") {
    feature_model <- "dinov2"
  }

  dims <- dim(mat)
  max_w <- dims[1]
  max_h <- dims[2]

  # Resolve boxes if missing
  if (is.null(boxes) || nrow(boxes) == 0L) {
    if (!is.null(conts) && length(conts) > 0L) {
      boxes <- .conts_to_boxes(conts, max_w = max_w, max_h = max_h)
    } else if (!is.null(mask_mat) && any(mask_mat > 0L)) {
      lbls <- bwlabel_cpp(mask_mat > 0L)
      if (max(lbls) > 0L) {
        c_list <- contour(lbls)
        boxes <- .conts_to_boxes(c_list, max_w = max_w, max_h = max_h)
      }
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Extracting deep object features and spatial map [{feature_model}]...",
      msg_done = "Deep features and spatial map ready"
    )
  }

  feat_mat <- .extract_object_features(
    mat = mat,
    boxes = boxes,
    feature_model = feature_model,
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  dino_sp <- image_features_dl(
    img = mat,
    model = "dinov2",
    return_pca = TRUE,
    interpolate = TRUE,
    threads = threads,
    engine = engine,
    device_id = device_id,
    verbose = FALSE,
    plot = FALSE,
    dir = dir
  )
  sp_map <- dino_sp$pca_image

  if (isTRUE(plot_features) && !is.null(sp_map)) {
    plot(sp_map)
  }

  if (isTRUE(verbose)) {
    cli::cli_progress_done()
  }

  if (is.list(obj) && !inherits(obj, c("Image", "image", "data.frame")) && !is.null(obj$boxes)) {
    obj$features <- feat_mat
    obj$feature_map <- sp_map
  }
  attr(obj, "features") <- feat_mat
  attr(obj, "feature_map") <- sp_map
  obj
}

#' Plot YOLO Feature Activation Heatmap
#'
#' Visualizes the 2D continuous activation heatmap or prototype feature energy
#' extracted during YOLO inference, either standalone or overlaid with alpha transparency
#' onto the original image.
#'
#' @param img An `Image` object, matrix, array, or an object returned by [image_segment_dl()]
#'   or [image_detect_dl()] with `return_features = TRUE`.
#' @param heatmap An `Image` or 2D matrix representing the activation heatmap. If `NULL`
#'   and `img` contains a `$heatmap` or `attr(img, "heatmap")`, it is extracted automatically.
#' @param col_palette Color palette name (`"inferno"`, `"magma"`, `"viridis"`, `"plasma"`,
#'   `"spectral"`), or a custom vector of colors. Default is `"inferno"`.
#' @param alpha Numeric transparency of the heatmap overlay on the original image (0 = only original image,
#'   1 = only heatmap). Default is `0.55`.
#' @param overlay Logical. If `TRUE` (default), blends the colored heatmap over the original image.
#'   If `FALSE`, plots the standalone colorized heatmap.
#' @param title Optional title for the plot.
#' @param plot Logical. Whether to display the plot immediately. Default is `TRUE`.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return An `Image` object containing the colorized / blended heatmap visualization, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("tomatoes.jpg")
#'   res <- image_segment_dl(img, model = "yolo26n-seg", return_features = TRUE)
#'   plot_yolo_heatmap(res)
#' }
plot_yolo_heatmap <- function(img,
                              heatmap = NULL,
                              col_palette = "inferno",
                              alpha = 0.55,
                              overlay = TRUE,
                              title = NULL,
                              plot = TRUE,
                              ...) {
  base_mat <- NULL
  if (is.null(heatmap)) {
    if (is.list(img) && !inherits(img, c("Image", "image")) && !is.null(img$heatmap)) {
      heatmap <- img$heatmap
    } else if (!is.null(attr(img, "heatmap"))) {
      heatmap <- attr(img, "heatmap")
    }
  }

  if (inherits(img, c("Image", "image"))) {
    base_mat <- image_data(img)
  } else if (is.array(img) && length(dim(img)) >= 2) {
    base_mat <- img
  } else if (is.list(img) && !is.null(attr(img, "image"))) {
    orig_img <- attr(img, "image")
    base_mat <- if (inherits(orig_img, c("Image", "image"))) image_data(orig_img) else orig_img
  }

  if (is.null(heatmap)) {
    stop("Could not find a valid activation heatmap. Ensure return_features = TRUE was used during inference.")
  }

  h_mat <- if (inherits(heatmap, c("Image", "image"))) image_data(heatmap) else as.matrix(heatmap)
  if (length(dim(h_mat)) == 3) h_mat <- h_mat[, , 1]

  min_val <- min(h_mat, na.rm = TRUE)
  max_val <- max(h_mat, na.rm = TRUE)
  rng <- if (max_val > min_val) (max_val - min_val) else 1.0
  norm_h <- (h_mat - min_val) / rng

  pal_colors <- switch(col_palette,
    "viridis"  = grDevices::hcl.colors(256, "Viridis"),
    "magma"    = grDevices::hcl.colors(256, "Inferno"),
    "inferno"  = grDevices::hcl.colors(256, "Inferno"),
    "plasma"   = grDevices::hcl.colors(256, "Plasma"),
    "spectral" = grDevices::hcl.colors(256, "Spectral"),
    if (is.character(col_palette) && length(col_palette) > 1) {
      grDevices::colorRampPalette(col_palette)(256)
    } else {
      grDevices::hcl.colors(256, "Inferno")
    }
  )

  rgb_mat <- grDevices::col2rgb(pal_colors) / 255.0
  idx <- pmax(1L, pmin(256L, round(norm_h * 255.0) + 1L))

  hw <- nrow(norm_h)
  hh <- ncol(norm_h)

  heat_r <- matrix(rgb_mat[1, idx], nrow = hw, ncol = hh)
  heat_g <- matrix(rgb_mat[2, idx], nrow = hw, ncol = hh)
  heat_b <- matrix(rgb_mat[3, idx], nrow = hw, ncol = hh)

  if (isTRUE(overlay) && !is.null(base_mat)) {
    bw <- dim(base_mat)[1]
    bh <- dim(base_mat)[2]
    bc <- if (length(dim(base_mat)) >= 3) dim(base_mat)[3] else 1L

    if (bw == hw && bh == hh) {
      base_norm <- if (is.raw(base_mat)) as.numeric(base_mat) / 255.0 else as.numeric(base_mat)
      dim(base_norm) <- dim(base_mat)

      br <- base_norm[, , 1]
      bg <- if (bc >= 2) base_norm[, , 2] else br
      bb <- if (bc >= 3) base_norm[, , 3] else br

      active_mask <- (norm_h > 0.05)
      blend_a <- pmin(1.0, pmax(0.0, alpha))

      out_r <- br
      out_g <- bg
      out_b <- bb

      out_r[active_mask] <- (1.0 - blend_a) * br[active_mask] + blend_a * heat_r[active_mask]
      out_g[active_mask] <- (1.0 - blend_a) * bg[active_mask] + blend_a * heat_g[active_mask]
      out_b[active_mask] <- (1.0 - blend_a) * bb[active_mask] + blend_a * heat_b[active_mask]

      out_arr <- array(0.0, dim = c(hw, hh, 3L))
      out_arr[, , 1] <- out_r
      out_arr[, , 2] <- out_g
      out_arr[, , 3] <- out_b
    } else {
      out_arr <- array(0.0, dim = c(hw, hh, 3L))
      out_arr[, , 1] <- heat_r
      out_arr[, , 2] <- heat_g
      out_arr[, , 3] <- heat_b
    }
  } else {
    out_arr <- array(0.0, dim = c(hw, hh, 3L))
    out_arr[, , 1] <- heat_r
    out_arr[, , 2] <- heat_g
    out_arr[, , 3] <- heat_b
  }

  res_img <- as_image(out_arr, colormode = "Color")
  attr(res_img, "heatmap") <- h_mat

  if (isTRUE(plot)) {
    plot(res_img, ...)
    if (!is.null(title)) {
      graphics::title(title)
    }
  }

  invisible(res_img)
}

#' Deep Learning Semantic & Instance Segmentation, Object Detection, and Visual Exemplar Counting
#'
#' Segments foreground objects (such as leaves, fruits, seeds, grains, roots,
#' animals, or plants) from complex, textured, or out-of-focus backgrounds using pre-trained
#' Deep Learning models in **ONNX** format.
#'
#' This unified function supports five major computer vision workflows:
#' 1. **Salient Foreground Segmentation & Background Removal:** Semantic cutout using architectures
#'    such as U2-Net, BEN2, BRIA RMBG 1.4/2.0, IS-Net, withoutBG, BiRefNet, and Silueta.
#' 2. **Real-Time YOLO26 Instance Segmentation (`yolo26n-seg`):** Real-time instance segmentation
#'    combining bounding box detection with 32 prototype mask channels and continuous bilinear sub-pixel
#'    interpolation. Returns a single multi-label mask image (`type = "mask"`), smooth polygon overlays
#'    (`type = "highlight"`), or background cutouts (`type = "segment"`), accompanied by per-instance
#'    morphological measurements (area, perimeter, circularity, etc.).
#' 3. **Star-Convex Object Detection & Segmentation (`stardist`):** Predicts 48 radial star-convex
#'    distance rays and object probability, specifically tailored for densely packed, touching, or convex
#'    structures such as seeds, grains, nuclei, and cells.
#' 4. **Zero-Shot Text-Prompted Detection & Instance Segmentation (Grounded-SAM):** Combines Grounding DINO
#'    with SAM 2.1 to locate, outline, and count arbitrary concepts described in natural language (e.g., `prompt = "leaf, fruit, insect"`).
#' 5. **One-Shot Visual Exemplar Segmentation & Counting (PerSAM):** Personalize Segment Anything via SAM 2.1.
#'    The user clicks on 1 or more exemplar objects (or provides coordinates), and the algorithm automatically searches,
#'    segments, and counts all visually similar objects across the entire image using Multi-Prototype Max-Pooling.
#'
#' @inheritParams image_binary_dl
#' @param type Character specifying the visual and output modality:
#'   * `"segment"` (default): Returns an `image` object where detected foreground is preserved
#'     and background is replaced by a solid color specified by `col_background`.
#'   * `"mask"`: Returns the binary mask (for semantic models) or single multi-label integer instance map
#'     (for YOLO-seg, StarDist, Grounded-SAM, and PerSAM) as an `image` object (`colormode = "Grayscale"`).
#'   * `"highlight"`: Overlays the original image with semi-transparent colored polygons covering
#'     each detected instance, ideal for visual inspection, quality control, and publication figures.
#'     When `mask = FALSE`, overlays bounding boxes directly.
#'   * `"boxes"`: Overlays bounding boxes on the original image and returns a data frame of detected boxes.
#' @param model Model architecture or path to a custom ONNX file. Built-in choices:
#'   * `"yolo26n-seg"` (or `"yolo26s-seg"`, `"yolo26m-seg"`, `"yolo26l-seg"`, `"yolo26x-seg"`): YOLO26 Instance Segmentation (real-time polygon mask & box detection, 80 COCO categories).
#'   * `"stardist"` (or `"stardist-dsb2018"`): StarDist star-convex radial polygon segmentation for overlapping grains/cells.
#'   * `"grounded-sam"`: Zero-shot open-vocabulary instance segmentation prompted by text.
#'   * `"persam"`: One-shot visual exemplar segmentation guided by reference click points.
#'   * `"u2netp"`, `"u2net"`: U2-Net salient foreground cutout (fast CPU-friendly default).
#'   * `"ben2"`: Boundary-aware Enhanced Network (BEN2) for razor-sharp leaf margins, serrations, and roots.
#'   * `"rmbg-1.4"`, `"rmbg-2.0"`: BRIA state-of-the-art commercial-grade background removal.
#'   * `"withoutbg"`, `"birefnet-lite"`, `"isnet-general-use"`, `"silueta"`: Specialized salient cutout backends.
#' @param col_background Character string or hex code specifying the solid background color when
#'   `type = "segment"`. Defaults to `"white"`. Accepts any valid R color name (e.g., `"black"`, `"transparent"`)
#'   or hexadecimal code (e.g., `"#FFFFFF"`).
#' @param col_highlight Fill color for the semi-transparent overlay polygon when `type = "highlight"`.
#'   Defaults to `"salmon"`. When multiple instances are detected, distinct categorical
#'   colors are automatically generated from a rainbow palette.
#' @param alpha Numeric transparency value in `[0, 1]` for the highlight polygon overlay when `type = "highlight"`.
#'   Default is `0.4` (40% opacity).
#' @param border Color for the border stroke of the highlight polygon when `type = "highlight"`.
#'   Default is `NA` (no border stroke).
#' @param lwd Numeric line width for the bounding boxes and highlight polygon contours. Default is `1`.
#' @param bbox Logical. If `TRUE`, draws YOLO-style bounding boxes with class labels and confidence scores
#'   around detected instances. Defaults to `FALSE`.
#' @param mask Logical. If `TRUE` (default), computes pixel-accurate instance or semantic masks.
#'   When `mask = FALSE` and `bbox = TRUE`, skips the mask decoder step and returns only the detected
#'   bounding box coordinates (fast object detection). Defaults to `TRUE`.
#' @param show_id Logical. If `TRUE`, overlays object instance IDs in bounding box
#'   labels or at their center of mass (computed via [poly_mass()]) when `type = "highlight"`
#'   and `bbox = FALSE`. Defaults to `FALSE`.
#' @param cex Size factor for bounding box label text (default `1.0`). Use smaller values (e.g., `0.2` or `0.5`)
#'   to reduce text size for dense scenes or many objects.
#' @param pad Padding multiplier controlling the background badge size around label text (default `1.0`).
#'   Can be a single value or a length-2 vector `c(pad_x, pad_y)`. Use `pad = 0` for a badge fitting tightly to the text.
#' @param show_text Logical. Whether to display text labels on bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in labels (default `TRUE`).
#' @param show_class Logical. Whether to display class names in labels (default `TRUE`).
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param rainbow Logical. If `TRUE`, assigns a distinct, unique color to each segmented object.
#'   If `FALSE` (the default), objects of the same class share the exact same color. Defaults to `FALSE`.
#' @param exemplar Logical. If `TRUE` (or if `model = "persam"`), activates one-shot visual exemplar segmentation (PerSAM)
#'   using SAM 2.1 embeddings. The user can click on 1 or more representative exemplar objects in the plot window
#'   (or pass coordinates via `prompt`), and the model will automatically compute visual feature embeddings,
#'   search the entire image for all similar instances via cosine similarity, and segment all matching objects.
#'   **Multi-Prototype Max-Pooling:** When multiple exemplars are clicked (e.g., one bright log, one dark log, one with moss),
#'   each click forms an independent 256-D prototype vector. The similarity map calculates the maximum cosine similarity
#'   against *any* prototype, preventing vector dilution and capturing wide visual diversity without prototype cancellation. Default is `FALSE`.
#' @param sim_threshold Numeric cosine similarity threshold in `[0, 1]` for PerSAM visual exemplar matching.
#'   Higher values (e.g. `0.55` - `0.65`) require closer visual resemblance, while lower values (e.g. `0.40` - `0.50`)
#'   capture more diverse, shadowed, or smaller instances. Default is `0.5`.
#' @param min_dist Numeric minimum spatial distance in pixels (within the 1024x1024 embedding space) between candidate
#'   peaks to prevent duplicate instance detections of the same object. Default is `16`.
#' @param box_threshold Numeric confidence threshold in `[0, 1]` for Grounded-SAM and YOLO detection bounding
#'   box proposals. Default is `0.25`.
#' @param text_threshold Numeric token logit confidence threshold in `[0, 1]` for Grounded-SAM text-to-box
#'   association. Default is `0.25`.
#' @param iou_threshold Numeric Non-Maximum Suppression (NMS) Intersection over Union (IoU) threshold in `[0, 1]`
#'   for deduplicating overlapping bounding box proposals. Default is `0.5`.
#' @param feat_res Spatial resolution / algorithm for exemplar matching. Options:
#'   * `64` (or `"fast"`, `"vit"`, `"low"`): Fast 64x64 ViT patch token matching (stride 16). Extremely fast
#'     and ideal for large objects (e.g. people, leaves, fruits, animals).
#'   * `256` (or `"medium"`): 256x256 multi-scale CNN feature fusion (stride 4).
#'   * `1024` (or `"dense"`, `"conv"`, `"high"`, `"grain"`): Dense Stride-1 Convolutional Exemplar matching
#'     at full 1024x1024 pixel resolution with single-pass Zero-Mean Normalized Cross-Correlation (ZNCC) and
#'     ViT semantic guidance. Essential for resolving small touching objects, seeds, and grains with sharp pixel boundaries.
#'   Defaults to `256`.
#' @param max_objects Maximum number of instances to detect and segment. Default is `300`.
#' @param superres_map Logical. If `TRUE`, applies Real-ESRGAN neural super-resolution to upscale the PerSAM similarity map for finer object boundaries and visualization. Defaults to `FALSE`.
#' @param return_features Logical. If `TRUE`, extracts deep feature representations for detected or
#'   segmented objects and the global image. Attaches `"features"` (an \eqn{N \times D} matrix of per-object
#'   deep embeddings, 512-D for CLIP or 384-D for DINOv2) and `"feature_map"` (dense spatial 3-component PCA false-color `Image`)
#'   as attributes to the returned object. Defaults to `FALSE`.
#' @param feature_model Foundation model used for feature representation when `return_features = TRUE`.
#'   Options are `"clip-vit-b32"` (default, OpenAI CLIP ViT-B/32 producing 512-D semantic embeddings) or
#'   `"dinov2"` (Meta DINOv2-S ViT producing 384-D self-supervised patch embeddings).
#' @param plot_features Logical. If `TRUE` and `return_features = TRUE`, visualizes the dense spatial PCA feature map. Defaults to `FALSE`.
#' @param show_mask Deprecated. If `TRUE`, equivalent to `type = "mask"`.
#'
#' @return
#'   * If `mask = FALSE` and `bbox = TRUE`: A data frame containing the bounding box coordinates
#'     (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `score`, `label`).
#'   * If `type = "segment"`: An `image` object containing the segmented image.
#'   * If `type = "mask"`: A single multi-label or binary `image` object (`colormode = "Grayscale"`).
#'   * If `type = "highlight"`: Invisibly returns a list with:
#'     - `"boxes"`: Data frame of detected bounding boxes and scores.
#'     - `"labels"`: Single multi-label integer mask `image` where each pixel value corresponds to instance `id`.
#'     - `"summary"`: Text summary of detected objects and counts.
#'     - `"counts"`: Data frame of instance counts per class label.
#'     - `"contours"`: List of polygon contour coordinate matrices.
#'     - `"features"`: If `return_features = TRUE`, numeric matrix of per-object deep feature embeddings.
#'     - `"feature_map"`: If `return_features = TRUE`, dense spatial 3-component PCA false-color `Image`.
#'     - `"similarity_map"`: Continuous cosine similarity `image` object (when `exemplar = TRUE`).
#'   * If `img` is a list: A list of the corresponding outputs.
#' @export
#'
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'
#'   # 1. Basic salient foreground segmentation with white background (CPU-friendly)
#'   seg <- image_segment_dl(img, model = "u2netp")
#'
#'   # 2. Razor-sharp plant boundary extraction (leaf margins, serrations, roots)
#'   seg_ben <- image_segment_dl(img, model = "ben2")
#'
#'   # 3. Semi-transparent polygon overlay for visual inspection
#'   image_segment_dl(img, model = "u2netp", type = "highlight")
#'
#'   # 4. Zero-shot text-prompted instance segmentation (Grounded-SAM)
#'   # Segments leaves and lesions, displaying bounding boxes and summary
#'   res_gs <- image_segment_dl(
#'     img,
#'     model = "grounded-sam",
#'     prompt = "plant leaf, lesion",
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#'   # Inspect counts and morphological measures
#'   print(attr(res_gs, "counts"))
#'   print(head(attr(res_gs, "features")))
#'
#'   # 5. Ultra-fast pure object detection without mask generation
#'   # (Skips SAM 2.1 completely for 10x faster execution)
#'   boxes <- image_segment_dl(
#'     img,
#'     model = "grounded-sam",
#'     prompt = "insect, leaf",
#'     mask = FALSE,
#'     bbox = TRUE
#'   )
#'
#'   # 6. One-shot visual exemplar segmentation (PerSAM):
#'   # Click on 1 or more representative items (e.g. seeds, logs, grains, cells)
#'   # to automatically find, outline, and count all matching instances:
#'   res_persam <- image_segment_dl(
#'     img,
#'     exemplar = TRUE,
#'     sim_threshold = 0.48,
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#'
#'   # 7. Multi-prototype visual exemplar segmentation with predefined coordinates:
#'   # Passing coordinates of distinct log types (light, dark, weathered)
#'   pts <- rbind(c(65, 380), c(425, 380), c(630, 380))
#'   res_logs <- image_segment_dl(
#'     img,
#'     exemplar = TRUE,
#'     prompt = pts,
#'     sim_threshold = 0.45,
#'     type = "highlight",
#'     engine = "gpu"
#'   )
#'   # 8. YOLO26 instance segmentation (real-time polygon mask & box detection)
#'   res_yolo <- image_segment_dl(
#'     img,
#'     model = "yolo26n-seg",
#'     type = "highlight",
#'     bbox = TRUE
#'   )
#' }
image_segment_dl <- function(img,
                             type = c("segment", "mask", "highlight"),
                             col_background = "white",
                             col_highlight = "salmon",
                             alpha = 0.4,
                             border = NA,
                             lwd = 1,
                             bbox = FALSE,
                             mask = TRUE,
                             cex = 1.0,
                             pad = 1.0,
                             show_text = TRUE,
                             show_conf = TRUE,
                             show_class = TRUE,
                             show_id = FALSE,
                             badge = TRUE,
                             rainbow = TRUE,
                             model = c("u2netp", "ben2", "yolo26n-seg", "stardist", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                             threshold = 0.5,
                             box_threshold = 0.25,
                             text_threshold = 0.25,
                             iou_threshold = 0.5,
                             exemplar = FALSE,
                             sim_threshold = 0.5,
                             min_dist = 16,
                             feat_res = 256,
                             max_objects = NULL,
                             superres_map = FALSE,
                             return_features = FALSE,
                             feature_model = c("yolo", "clip-vit-b32", "dinov2"),
                             plot_features = FALSE,
                             threads = 0,
                             engine = c("gpu", "cpu"),
                             device_id = -1,
                             fill_hull = TRUE,
                             filter = 0,
                             erode = 0,
                             dilate = 0,
                             opening = 0,
                             closing = 0,
                             min_area = 0,
                             invert = FALSE,
                             pick_object = FALSE,
                             prompt = NULL,
                             show_mask = FALSE,
                             verbose = TRUE,
                             plot = TRUE,
                             dir = pliman_model_dir(),
                             ...) {
  engine <- match.arg(engine)
  dots <- list(...)
  if ("label_size" %in% names(dots)) cex <- dots$label_size
  if ("text_size" %in% names(dots)) cex <- dots$text_size
  if ("font_scale" %in% names(dots)) cex <- dots$font_scale
  if ("cex_scale" %in% names(dots)) cex <- dots$cex_scale
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  if (length(cex) > 1L) {
    if (missing(pad)) pad <- cex[2]
    cex <- cex[1]
  }

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_segment_dl(x,
                       type = type,
                       col_background = col_background,
                       col_highlight = col_highlight,
                       alpha = alpha,
                       border = border,
                       lwd = lwd,
                       bbox = bbox,
                       mask = mask,
                       cex = cex,
                       pad = pad,
                       show_text = show_text,
                       show_conf = show_conf,
                       show_class = show_class,
                       show_id = show_id,
                       badge = badge,
                       rainbow = rainbow,
                       model = model,
                       threshold = threshold,
                       box_threshold = box_threshold,
                       text_threshold = text_threshold,
                       iou_threshold = iou_threshold,
                       exemplar = exemplar,
                       sim_threshold = sim_threshold,
                       min_dist = min_dist,
                       feat_res = feat_res,
                       max_objects = max_objects,
                       superres_map = superres_map,
                       return_features = return_features,
                       feature_model = feature_model,
                       plot_features = plot_features,
                       threads = threads,
                       engine = engine,
                       device_id = device_id,
                       fill_hull = fill_hull,
                       filter = filter,
                       erode = erode,
                       dilate = dilate,
                       opening = opening,
                       closing = closing,
                       min_area = min_area,
                       invert = invert,
                       pick_object = pick_object,
                       prompt = prompt,
                       show_mask = show_mask,
                       verbose = verbose,
                       plot = FALSE,
                       dir = dir,
                       ...)
    })
    return(res)
  }

  if (!isTRUE(mask) && !isTRUE(bbox)) {
    cli::cli_abort("At least one of {.arg mask} or {.arg bbox} must be TRUE.")
  }

  if (isTRUE(show_mask)) {
    type <- "mask"
  } else {
    type_str <- tolower(type[1])
    if (grepl("^high", type_str)) {
      type <- "highlight"
    } else {
      type <- match.arg(type, c("segment", "mask", "highlight"))
    }
  }

  is_url <- is.character(model) && grepl("^https?://", model[1])
  if (!is_url && is.character(model)) {
    model <- .resolve_model_name(model[1])
  }

  if (!is_url && grepl("[-_]cls", tolower(model))) {
    cli::cli_abort(c(
      "Model {.val {model}} is an image classification model, not an instance segmentation model.",
      "i" = "Please use {.fn image_classify_dl} instead."
    ))
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  # Grounded-SAM, PerSAM, YOLO-Seg, and StarDist instance segmentation
  is_text_prompt <- is.character(prompt) && length(prompt) >= 1 && !all(prompt %in% c("center", "box", "exemplar"))
  is_sam_model <- grepl("sam", model, ignore.case = TRUE)
  use_persam <- isTRUE(exemplar) || (model == "persam") || (is_sam_model && identical(prompt, "exemplar"))
  use_grounded_sam <- !use_persam && ((model == "grounded-sam") || (is_sam_model && is_text_prompt))
  non_yolo_models <- c(
    "u2netp", "u2net", "birefnet-lite", "birefnet", "isnet-general-use", "isnet",
    "rmbg-1.4", "rmbg-2.0", "rmbg", "silueta", "withoutbg", "ben2", "ben",
    "sam2.1", "sam3.1", "sam", "grounded-sam", "persam", "stardist", "star-dist",
    "depth-anything-v2", "dinov2", "realesrgan-compact",
    "yolov8s-world", "yolo-world", "yoloworld", "yolov8s-worldv2"
  )

  if (grepl("world", tolower(model[1]))) {
    if (isTRUE(bbox) && !isTRUE(mask)) {
      world_classes <- if (!is.null(prompt)) prompt else if (!is.null(list(...)$labels)) list(...)$labels else if (!is.null(list(...)$classes)) list(...)$classes else c("object")
      return(image_detect_world(
        img = mat,
        classes = world_classes,
        conf_threshold = box_threshold,
        iou_threshold = iou_threshold,
        rainbow = rainbow,
        col = if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL,
        lwd = lwd,
        plot = plot,
        verbose = verbose,
        dir = dir,
        ...
      ))
    } else {
      stop(
        "O modelo '", model[1], "' e um detector de caixas delimitadoras (bounding boxes) por texto aberto (.pt), nao um modelo de segmentacao de mascaras/poligonos.\n",
        "  - Para segmentacao de instancias com YOLO: use model = 'yolo26n-seg' (ou yolo26s-seg, yolo26l-seg, etc.)\n",
        "  - Para segmentacao guiada por texto: use model = 'grounded-sam' com prompt = 'sua classe'\n",
        "  - Para deteccao de caixas com YOLO-World: use image_detect_dl(img, model = '", model[1], "', labels = c('sua classe'))"
      )
    }
  }

  is_known_non_yolo <- model %in% non_yolo_models

  is_custom_yolo <- !is_known_non_yolo && (
    grepl("yolo", model, ignore.case = TRUE) ||
    grepl("[-_]seg", model, ignore.case = TRUE) ||
    grepl("[-_]det", model, ignore.case = TRUE)
  )
  use_yolo_seg <- !use_persam && !use_grounded_sam && !is_known_non_yolo && (
    model %in% c("yolo26n-seg", "yolo26s-seg", "yolo26m-seg", "yolo26l-seg", "yolo26x-seg",
                 "yolo11n-seg", "yolo11s-seg", "yolo11m-seg", "yolo11l-seg", "yolo11x-seg", "yolo-seg") ||
    is_custom_yolo
  )
  use_stardist <- !use_persam && !use_grounded_sam && !use_yolo_seg && (model %in% c("stardist", "star-dist"))

  if (use_persam || use_grounded_sam || use_yolo_seg || use_stardist) {
    if (use_persam) {
      if (!missing(box_threshold) && missing(sim_threshold)) {
        sim_threshold <- box_threshold
      }
      dots <- list(...)
      ex_pts <- if (is.numeric(prompt) || is.matrix(prompt) || is.data.frame(prompt)) {
        prompt
      } else if (!is.null(dots$exemplar_points)) {
        dots$exemplar_points
      } else {
        NULL
      }
      gs_res <- .run_persam(
        mat = mat,
        exemplar_points = ex_pts,
        sim_threshold = sim_threshold,
        min_dist = min_dist,
        iou_threshold = iou_threshold,
        max_objects = max_objects,
        feat_res = feat_res,
        threshold = threshold,
        superres_map = superres_map,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        invert = invert,
        mask = mask,
        verbose = verbose,
        dir = dir
      )
    } else if (use_grounded_sam) {
      gs_res <- .run_grounded_sam(
        mat = mat,
        prompt = prompt,
        threshold = threshold,
        box_threshold = box_threshold,
        text_threshold = text_threshold,
        iou_threshold = iou_threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        invert = invert,
        mask = mask,
        verbose = verbose,
        dir = dir
      )
    } else if (use_yolo_seg) {
      gs_res <- .run_yolo(
        mat = mat,
        model = model,
        conf_threshold = box_threshold,
        iou_threshold = iou_threshold,
        labels = if (!is.null(list(...)$labels)) list(...)$labels else NULL,
        threads = threads,
        engine = engine,
        device_id = device_id,
        fill_hull = fill_hull,
        filter = filter,
        erode = erode,
        dilate = dilate,
        opening = opening,
        closing = closing,
        min_area = min_area,
        need_masks = isTRUE(mask) || type == "mask" || type == "highlight",
        return_features = return_features,
        dir = dir
      )
    } else if (use_stardist) {
      gs_res <- .run_stardist(
        mat = mat,
        model = model,
        prob_threshold = threshold,
        nms_threshold = iou_threshold,
        threads = threads,
        engine = engine,
        device_id = device_id,
        dir = dir
      )
    }

    sum_title <- if (use_persam) {
      "PerSAM Detection Summary"
    } else if (use_yolo_seg) {
      "YOLO Detection Summary"
    } else if (use_stardist) {
      "StarDist Detection Summary"
    } else {
      "Grounded-SAM Detection Summary"
    }

    if (nrow(gs_res$boxes) == 0L && isTRUE(verbose)) {
      cli::cli_alert_warning("No objects were detected above the confidence threshold ({.code box_threshold = {box_threshold}}).")
      cli::cli_alert_info("Tip: If the model was trained for few epochs (e.g. 3), try {.code box_threshold = 0.01} or train for 30-50 epochs so the model learns confident scores.")
    }

    if (!isTRUE(mask) && isTRUE(bbox)) {
      if (isTRUE(plot)) {
        plot(as_image(mat), ...)
        num_inst <- nrow(gs_res$boxes)
        if (num_inst > 0) {
          user_hl <- if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL
          palette_colors <- .get_class_palette(
            labels = gs_res$boxes$label,
            n = num_inst,
            user_col = user_hl,
            default_col = if (!is.null(user_hl)) user_hl else "#00CC66",
            rainbow = rainbow
          )
          .plot_yolo_bboxes(
            boxes = gs_res$boxes,
            palette_colors = palette_colors,
            lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd,
            cex = cex,
            pad = pad,
            show_text = show_text,
            show_conf = show_conf,
            show_class = show_class,
            show_id = show_id,
            badge = badge
          )
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(gs_res$boxes, "summary") <- gs_res$summary
      attr(gs_res$boxes, "counts") <- gs_res$counts
      attr(gs_res$boxes, "image") <- as_image(mat)
      if (!is.null(gs_res$similarity_map)) attr(gs_res$boxes, "similarity_map") <- gs_res$similarity_map
      gs_res$boxes <- .attach_features(gs_res$boxes, mat = mat, boxes = gs_res$boxes,
                                       return_features = return_features, feature_model = feature_model,
                                       plot_features = plot_features, threads = threads, engine = engine,
                                       device_id = device_id, dir = dir, verbose = verbose,
                                       yolo_features = gs_res)
      return(invisible(gs_res$boxes))
    }

    if (type == "mask") {
      if (!isTRUE(mask) || is.null(gs_res$labels)) {
        cli::cli_alert_warning("type = 'mask' requested but mask = FALSE. Returning bounding boxes instead.")
        return(invisible(gs_res$boxes))
      }
      out_img <- gs_res$labels
      if (isTRUE(plot)) {
        plot(out_img, ...)
        if (isTRUE(bbox) && nrow(gs_res$boxes) > 0) {
          num_inst <- nrow(gs_res$boxes)
          user_hl <- if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL
          palette_colors <- .get_class_palette(
            labels = gs_res$boxes$label,
            n = num_inst,
            user_col = user_hl,
            default_col = if (!is.null(user_hl)) user_hl else "#00CC66",
            rainbow = rainbow
          )
          for (k in seq_len(num_inst)) {
            bx <- gs_res$boxes[k, ]
            k_col <- palette_colors[((k - 1) %% length(palette_colors)) + 1]
            graphics::rect(xleft = bx$xmin, ybottom = bx$ymin, xright = bx$xmax, ytop = bx$ymax,
                           border = k_col, lwd = lwd)
          }
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(out_img, "boxes") <- gs_res$boxes
      attr(out_img, "summary") <- gs_res$summary
      attr(out_img, "counts") <- gs_res$counts
      attr(out_img, "image") <- as_image(mat)
      if (!is.null(gs_res$similarity_map)) attr(out_img, "similarity_map") <- gs_res$similarity_map
      out_img <- .attach_features(out_img, mat = mat, boxes = gs_res$boxes,
                                  return_features = return_features, feature_model = feature_model,
                                  plot_features = plot_features, threads = threads, engine = engine,
                                  device_id = device_id, dir = dir, verbose = verbose,
                                  yolo_features = gs_res)
      return(invisible(out_img))
    }

    if (type == "segment") {
      out <- as_image(mat)
      if (isTRUE(mask) && !is.null(gs_res$mask)) {
        col_rgb <- tryCatch(grDevices::col2rgb(col_background) / 255.0, error = function(e) c(1, 1, 1))
        m_bg <- gs_res$mask
        if (length(dims) == 3) {
          for (k in seq_len(dims[3])) {
            fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
            out[, , k][!m_bg] <- fill_val
          }
        } else {
          fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
          out[!m_bg] <- fill_val
        }
      }

      if (isTRUE(plot)) {
        plot(out, ...)
        if ((isTRUE(bbox) || !isTRUE(mask)) && nrow(gs_res$boxes) > 0) {
          num_inst <- nrow(gs_res$boxes)
          user_hl <- if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL
          palette_colors <- .get_class_palette(
            labels = gs_res$boxes$label,
            n = num_inst,
            user_col = user_hl,
            default_col = if (!is.null(user_hl)) user_hl else "#00CC66",
            rainbow = rainbow
          )
          .plot_yolo_bboxes(
            boxes = gs_res$boxes,
            palette_colors = palette_colors,
            lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd,
            cex = cex,
            pad = pad,
            show_text = show_text,
            show_conf = show_conf,
            show_class = show_class,
            show_id = show_id,
            badge = badge
          )
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      attr(out, "boxes") <- gs_res$boxes
      attr(out, "labels") <- gs_res$labels
      attr(out, "summary") <- gs_res$summary
      attr(out, "counts") <- gs_res$counts
      attr(out, "image") <- as_image(mat)
      if (!is.null(gs_res$similarity_map)) attr(out, "similarity_map") <- gs_res$similarity_map
      out <- .attach_features(out, mat = mat, boxes = gs_res$boxes,
                              return_features = return_features, feature_model = feature_model,
                              plot_features = plot_features, threads = threads, engine = engine,
                              device_id = device_id, dir = dir, verbose = verbose,
                              yolo_features = gs_res)
      return(invisible(out))
    }

    if (type == "highlight") {
      num_inst <- nrow(gs_res$boxes)
      conts_all <- vector("list", num_inst)

      has_boxes <- !is.null(gs_res$boxes) && is.data.frame(gs_res$boxes) &&
                   nrow(gs_res$boxes) == num_inst && all(c("xmin", "ymin", "xmax", "ymax") %in% names(gs_res$boxes))

      if (isTRUE(mask) && !is.null(gs_res$labels) && num_inst > 0) {
        lbl_mat <- image_data(gs_res$labels)
        pad <- 2L
        for (k in seq_len(num_inst)) {
          if (has_boxes) {
            bx <- gs_res$boxes[k, ]
            x1 <- max(1L, floor(bx$xmin) - pad)
            x2 <- min(orig_w, ceiling(bx$xmax) + pad)
            y1 <- max(1L, floor(bx$ymin) - pad)
            y2 <- min(orig_h, ceiling(bx$ymax) + pad)
            sub_lbl <- lbl_mat[x1:x2, y1:y2, drop = FALSE]
            sub_m <- (sub_lbl == k)
            if (!any(sub_m)) next
            lbls <- bwlabel_cpp(sub_m)
            conts <- list()
            if (max(lbls) > 0) {
              conts <- contour(lbls)
              for (ci in seq_along(conts)) {
                if (is.matrix(conts[[ci]])) {
                  conts[[ci]][, 1] <- conts[[ci]][, 1] + (x1 - 1L)
                  conts[[ci]][, 2] <- conts[[ci]][, 2] + (y1 - 1L)
                }
              }
            }
          } else {
            inst_m <- (lbl_mat == k)
            if (!any(inst_m)) next
            lbls <- bwlabel_cpp(inst_m)
            conts <- if (max(lbls) > 0) contour(lbls) else list()
          }
          conts_all[[k]] <- conts
        }
      }

      # Compute instance centroids using poly_mass (center of mass)
      cx <- numeric(num_inst)
      cy <- numeric(num_inst)
      for (k in seq_len(num_inst)) {
        conts_k <- conts_all[[k]]
        pm <- NULL
        if (is.list(conts_k) && length(conts_k) > 0L) {
          pt_counts <- vapply(conts_k, function(p) if (is.matrix(p)) nrow(p) else 0L, integer(1))
          if (max(pt_counts) >= 3L) {
            main_p <- conts_k[[which.max(pt_counts)]]
            pm <- tryCatch(poly_mass(main_p), error = function(e) NULL)
          }
        }
        if (!is.null(pm) && is.matrix(pm) && nrow(pm) >= 1L && is.finite(pm[1, "x"]) && is.finite(pm[1, "y"])) {
          cx[k] <- pm[1, "x"]
          cy[k] <- pm[1, "y"]
        } else if (has_boxes && k <= nrow(gs_res$boxes)) {
          cx[k] <- (gs_res$boxes$xmin[k] + gs_res$boxes$xmax[k]) / 2
          cy[k] <- (gs_res$boxes$ymin[k] + gs_res$boxes$ymax[k]) / 2
        } else {
          cx[k] <- 0
          cy[k] <- 0
        }
      }

      if (isTRUE(plot)) {
        if (isTRUE(verbose)) {
          cli::cli_progress_step(
            msg = if (isTRUE(mask)) "Rendering highlight overlay..." else "Rendering detection overlay...",
            msg_done = if (isTRUE(mask)) "Highlight overlay ready" else "Detection overlay ready"
          )
        }
        plot(as_image(mat), ...)
        if (num_inst > 0) {
          user_hl <- if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL
          palette_colors <- .get_class_palette(
            labels = gs_res$boxes$label,
            n = num_inst,
            user_col = user_hl,
            default_col = if (!is.null(user_hl)) user_hl else "#00CC66",
            rainbow = rainbow
          )

          if (isTRUE(mask) && !is.null(gs_res$labels)) {
            for (k in seq_len(num_inst)) {
              conts <- conts_all[[k]]
              if (is.null(conts) || length(conts) == 0) next

              k_col <- palette_colors[((k - 1) %% length(palette_colors)) + 1]
              c_rgb <- tryCatch(grDevices::col2rgb(k_col) / 255.0, error = function(e) c(0, 1, 0))
              poly_col <- grDevices::rgb(c_rgb[1], c_rgb[2], c_rgb[3], alpha = alpha)
              cur_border <- if (is.na(border)) NA else k_col

              for (p in conts) {
                if (is.matrix(p) && nrow(p) >= 3) {
                  p_plot <- .decimate_poly(p)
                  graphics::polygon(p_plot[, 1], p_plot[, 2], col = poly_col, border = cur_border, lwd = lwd)
                }
              }
            }
          }

          if (isTRUE(bbox) || !isTRUE(mask)) {
            .plot_yolo_bboxes(
              boxes = gs_res$boxes,
              palette_colors = palette_colors,
              lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd,
              cex = cex,
              pad = pad,
              show_text = show_text,
              show_conf = show_conf,
              show_class = show_class,
              show_id = show_id,
              badge = badge
            )
          } else if (isTRUE(show_id)) {
            graphics::points(cx, cy, pch = 21, bg = "black", col = "white", cex = 2.2)
            graphics::text(cx, cy, labels = as.character(seq_len(num_inst)), col = "white", font = 2, cex = 0.9)
          }
        }
        if (isTRUE(verbose)) {
          cli::cli_progress_done()
        }
      }
      if (isTRUE(verbose)) {
        .print_detection_summary(gs_res$summary, title = sum_title)
      }
      flat_conts <- unlist(conts_all, recursive = FALSE)
      if (is.null(flat_conts)) flat_conts <- list()
      attr(gs_res, "contours") <- flat_conts
      attr(gs_res, "centroids") <- data.frame(id = seq_len(num_inst), x = cx, y = cy)
      if (has_boxes && is.data.frame(gs_res$boxes) && nrow(gs_res$boxes) == num_inst) {
        gs_res$boxes$cx <- round(cx, 1)
        gs_res$boxes$cy <- round(cy, 1)
      }
      attr(gs_res, "summary") <- gs_res$summary
      attr(gs_res, "counts") <- gs_res$counts
      attr(gs_res, "image") <- as_image(mat)
      if (!is.null(gs_res$similarity_map)) attr(gs_res, "similarity_map") <- gs_res$similarity_map
      gs_res <- .attach_features(gs_res, mat = mat, boxes = gs_res$boxes, conts = flat_conts,
                                 return_features = return_features, feature_model = feature_model,
                                 plot_features = plot_features, threads = threads, engine = engine,
                                 device_id = device_id, dir = dir, verbose = verbose,
                                 yolo_features = gs_res)
      return(invisible(gs_res))
    }
  }

  # 1. Compute refined mask with all morphological operations applied
  mask_mat <- .compute_dl_mask(mat = mat,
                               model = model,
                               threshold = threshold,
                               fill_hull = fill_hull,
                               filter = filter,
                               erode = erode,
                               dilate = dilate,
                               opening = opening,
                               closing = closing,
                               min_area = min_area,
                               invert = invert,
                               pick_object = pick_object,
                               prompt = prompt,
                               exemplar = exemplar,
                               threads = threads,
                               engine = engine,
                               device_id = device_id,
                               verbose = verbose,
                               dir = dir)

  if (!isTRUE(mask) && isTRUE(bbox)) {
    lbls <- bwlabel_cpp(mask_mat > 0)
    conts <- if (max(lbls) > 0) contour(lbls) else list()
    num_inst <- length(conts)
    if (num_inst > 0) {
      boxes_list <- lapply(seq_along(conts), function(i) {
        p <- conts[[i]]
        data.frame(
          id = i,
          xmin = round(min(p[, 1]), 1),
          ymin = round(min(p[, 2]), 1),
          xmax = round(max(p[, 1]), 1),
          ymax = round(max(p[, 2]), 1),
          score = 1.0,
          label = "object",
          stringsAsFactors = FALSE
        )
      })
      boxes_df <- do.call(rbind, boxes_list)
    } else {
      boxes_df <- data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                             xmax = numeric(0), ymax = numeric(0), score = numeric(0),
                             label = character(0), stringsAsFactors = FALSE)
    }
    if (isTRUE(plot)) {
      plot(as_image(mat), ...)
      if (nrow(boxes_df) > 0) {
        user_hl <- if (!missing(col_highlight)) col_highlight else if (!is.null(list(...)$col)) list(...)$col else NULL
        palette_colors <- .get_class_palette(
          labels = boxes_df$label,
          n = nrow(boxes_df),
          user_col = user_hl,
          default_col = if (!is.null(user_hl)) user_hl else "#00CC66",
          rainbow = rainbow
        )
        .plot_yolo_bboxes(
          boxes = boxes_df,
          palette_colors = palette_colors,
          lwd = if (is.null(lwd) || is.na(lwd) || lwd <= 1) 2 else lwd,
          cex = cex,
          pad = pad,
          show_text = show_text,
          show_conf = show_conf,
          show_class = show_class,
          show_id = show_id,
          badge = badge
        )
      }
    }
    boxes_df <- .attach_features(boxes_df, mat = mat, boxes = boxes_df, conts = conts,
                                 return_features = return_features, feature_model = feature_model,
                                 plot_features = plot_features, threads = threads, engine = engine,
                                 device_id = device_id, dir = dir, verbose = verbose)
    return(invisible(boxes_df))
  }

  # 2. Output according to selected type
  if (type == "mask") {
    mask_img <- as_image(mask_mat, colormode = "Grayscale")
    if (isTRUE(plot)) {
      plot(mask_img, ...)
    }
    mask_img <- .attach_features(mask_img, mat = mat, mask_mat = mask_mat,
                                 return_features = return_features, feature_model = feature_model,
                                 plot_features = plot_features, threads = threads, engine = engine,
                                 device_id = device_id, dir = dir, verbose = verbose)
    return(invisible(mask_img))
  }

  if (type == "segment") {
    out <- as_image(mat)
    col_rgb <- tryCatch(grDevices::col2rgb(col_background) / 255.0, error = function(e) c(1, 1, 1))
    if (length(dims) == 3) {
      for (k in seq_len(dims[3])) {
        fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[k] * 255)) else col_rgb[k]
        out[, , k][!mask_mat] <- fill_val
      }
    } else {
      fill_val <- if (is.raw(mat)) as.raw(round(col_rgb[1] * 255)) else col_rgb[1]
      out[!mask_mat] <- fill_val
    }

    if (isTRUE(plot)) {
      plot(out, ...)
    }
    out <- .attach_features(out, mat = mat, mask_mat = mask_mat,
                            return_features = return_features, feature_model = feature_model,
                            plot_features = plot_features, threads = threads, engine = engine,
                            device_id = device_id, dir = dir, verbose = verbose)
    return(invisible(out))
  }

  if (type == "highlight") {
    lbls <- bwlabel_cpp(mask_mat > 0)
    conts <- if (max(lbls) > 0) contour(lbls) else list()

    if (isTRUE(plot)) {
      if (isTRUE(verbose)) {
        cli::cli_progress_step(
          msg = "Rendering highlight overlay...",
          msg_done = "Highlight overlay ready"
        )
      }
      plot(as_image(mat), ...)
      if (isTRUE(rainbow) && length(conts) > 0) {
        palette_pool <- c(
          "#00CC66", "#0099FF", "#FF9900", "#FF3366", "#9933FF",
          "#00E5FF", "#FFD700", "#FF6600", "#00B4D8", "#7209B7",
          "#2EC4B6", "#E71D36", "#FF9F1C", "#4361EE", "#F72585"
        )
        pal <- if (length(conts) <= length(palette_pool)) palette_pool[seq_along(conts)] else grDevices::rainbow(length(conts), s = 0.85, v = 0.95)
      } else {
        pal <- rep(col_highlight, length(conts))
      }
      for (ci in seq_along(conts)) {
        p <- conts[[ci]]
        k_col <- pal[ci]
        c_rgb <- tryCatch(grDevices::col2rgb(k_col) / 255.0, error = function(e) c(0, 1, 0))
        poly_col <- grDevices::rgb(c_rgb[1], c_rgb[2], c_rgb[3], alpha = alpha)
        cur_border <- if (is.na(border)) NA else k_col
        if (is.matrix(p) && nrow(p) >= 3) {
          p_plot <- .decimate_poly(p)
          graphics::polygon(p_plot[, 1], p_plot[, 2], col = poly_col, border = cur_border, lwd = lwd)
        }
        if (isTRUE(bbox) && is.matrix(p) && nrow(p) >= 1) {
          graphics::rect(xleft = min(p[, 1]), ybottom = min(p[, 2]),
                         xright = max(p[, 1]), ytop = max(p[, 2]),
                         border = k_col, lwd = lwd)
        }
      }
      if (isTRUE(verbose)) {
        cli::cli_progress_done()
      }
    }
    conts <- .attach_features(conts, mat = mat, conts = conts, mask_mat = mask_mat,
                              return_features = return_features, feature_model = feature_model,
                              plot_features = plot_features, threads = threads, engine = engine,
                              device_id = device_id, dir = dir, verbose = verbose)
    return(invisible(conts))
  }
}

#' Clear Cached ONNX Runtime Sessions
#'
#' Releases all in-memory compiled neural network sessions (e.g., Grounded-SAM,
#' SAM 2.1, U2-Net) and frees allocated RAM.
#'
#' @return Invisible \code{NULL}.
#' @export
pliman_clear_sessions <- function() {
  clear_onnx_sessions_cpp()
  if (exists(".pliman_vocab_cache")) {
    rm(list = ls(envir = .pliman_vocab_cache), envir = .pliman_vocab_cache)
  }
  invisible(NULL)
}

#' @rdname image_remove_bg_dl
#' @export
image_remove_bg <- function(img,
                            model = c("u2netp", "ben2", "grounded-sam", "sam2.1", "persam", "rmbg-1.4", "rmbg-2.0", "withoutbg", "birefnet-lite", "sam3.1", "isnet-general-use", "silueta", "u2net"),
                            threshold = 0.5,
                            fill_hull = TRUE,
                            filter = 0,
                            erode = 0,
                            dilate = 0,
                            opening = 0,
                            closing = 0,
                            min_area = 0,
                            invert = FALSE,
                            pick_object = FALSE,
                            prompt = NULL,
                            exemplar = FALSE,
                            threads = 0,
                            engine = c("gpu", "cpu"),
                            device_id = -1,
                            verbose = TRUE,
                            transparent = TRUE,
                            bg_color = "black",
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  image_remove_bg_dl(img = img,
                     model = model,
                     threshold = threshold,
                     fill_hull = fill_hull,
                     filter = filter,
                     erode = erode,
                     dilate = dilate,
                     opening = opening,
                     closing = closing,
                     min_area = min_area,
                     invert = invert,
                     pick_object = pick_object,
                     prompt = prompt,
                     exemplar = exemplar,
                     threads = threads,
                     engine = engine,
                     device_id = device_id,
                     verbose = verbose,
                     transparent = transparent,
                     bg_color = bg_color,
                     plot = plot,
                     dir = dir,
                     ...)
}

#' Estimate Relative 3D Depth with Depth Anything V2
#'
#' Estimates a metric-aligned relative 3D depth map from any 2D RGB image using
#' the state-of-the-art Depth Anything V2 foundation model running purely via ONNX Runtime.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"depth-anything-v2"`.
#' @param col_palette Color palette used to visualize the depth map when `plot = TRUE`.
#'   Options include `"viridis"`, `"magma"`, `"inferno"`, `"plasma"`, `"spectral"`, or a custom vector of colors.
#' @param invert Logical. If `TRUE`, inverts the depth map so near objects are darker. Defaults to `FALSE`.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the resulting depth heatmap.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return An `Image` object representing the colorized depth map, with the raw normalized
#'   depth matrix stored in the `"depth"` attribute (`attr(res, "depth")`).
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'   depth <- image_depth_dl(img)
#' }
image_depth_dl <- function(img,
                           model = "depth-anything-v2",
                           col_palette = "magma",
                           invert = FALSE,
                           threads = 0,
                           engine = c("gpu", "cpu"),
                           device_id = -1,
                           verbose = TRUE,
                           plot = TRUE,
                           dir = pliman_model_dir(),
                           ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_depth_dl(x, model = model, col_palette = col_palette, invert = invert,
                     threads = threads, engine = engine, device_id = device_id,
                     verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "depth-anything-v2"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running Depth Anything V2 [{toupper(engine)}]...",
      msg_done = "Depth estimation complete"
    )
  }

  in_size <- 518L
  tensor <- .preprocess_nchw(
    mat,
    target_size = in_size,
    mean = c(0.485, 0.456, 0.406),
    std = c(0.229, 0.224, 0.225),
    letterbox = FALSE
  )

  depth_mat <- run_depth_anything_cpp(
    tensor_vec = tensor,
    in_w = in_size,
    in_h = in_size,
    orig_w = orig_w,
    orig_h = orig_h,
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  min_d <- min(depth_mat, na.rm = TRUE)
  max_d <- max(depth_mat, na.rm = TRUE)
  rng <- if (max_d > min_d) (max_d - min_d) else 1.0
  norm_d <- (depth_mat - min_d) / rng
  if (isTRUE(invert)) norm_d <- 1.0 - norm_d

  pal_colors <- switch(col_palette,
    "viridis"  = grDevices::hcl.colors(256, "Viridis"),
    "magma"    = grDevices::hcl.colors(256, "Inferno"),
    "inferno"  = grDevices::hcl.colors(256, "Inferno"),
    "plasma"   = grDevices::hcl.colors(256, "Plasma"),
    "spectral" = grDevices::hcl.colors(256, "Spectral"),
    if (is.character(col_palette) && length(col_palette) > 1) {
      grDevices::colorRampPalette(col_palette)(256)
    } else {
      grDevices::hcl.colors(256, "Viridis")
    }
  )

  rgb_mat <- grDevices::col2rgb(pal_colors) / 255.0
  idx <- pmax(1L, pmin(256L, round(norm_d * 255.0) + 1L))

  out_arr <- array(0.0, dim = c(orig_w, orig_h, 3L))
  out_arr[, , 1] <- matrix(rgb_mat[1, idx], nrow = orig_w, ncol = orig_h)
  out_arr[, , 2] <- matrix(rgb_mat[2, idx], nrow = orig_w, ncol = orig_h)
  out_arr[, , 3] <- matrix(rgb_mat[3, idx], nrow = orig_w, ncol = orig_h)

  res_img <- as_image(out_arr)
  attr(res_img, "depth") <- depth_mat

  if (isTRUE(plot)) {
    plot(res_img, ...)
  }

  invisible(res_img)
}

#' Extract Foundation Features and Semantic Visualization with DINOv2
#'
#' Extracts dense patch embeddings and global representation vectors from DINOv2
#' (Vision Transformer ViT-S/14), and computes top-3 PCA false-color semantic segmentation
#' without any labels or supervised training.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"dinov2"`.
#' @param patch_size Size of ViT patches in pixels (default 14).
#' @param return_pca Logical. If `TRUE`, computes top-3 principal components of patch embeddings.
#' @param interpolate Logical. If `TRUE`, interpolates the PCA map from patch grid (37x37) to original image resolution.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the PCA false-color semantic image.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list containing:
#'   * `cls_token`: 384-dimensional global semantic representation vector.
#'   * `pca_image`: An `Image` object representing the 3-component PCA false-color semantic map.
#'   * `pca_r`, `pca_g`, `pca_b`: Matrices of normalized PCA coordinates.
#'   * `num_patches`: Total number of spatial patches.
#'   * `embed_dim`: Dimensionality of patch embeddings.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'   feats <- image_features_dl(img)
#' }
image_features_dl <- function(img,
                              model = "dinov2",
                              patch_size = 14,
                              return_pca = TRUE,
                              interpolate = TRUE,
                              threads = 0,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_features_dl(x, model = model, patch_size = patch_size, return_pca = return_pca,
                        interpolate = interpolate, threads = threads, engine = engine,
                        device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "dinov2"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running DINOv2 foundation ViT [{toupper(engine)}]...",
      msg_done = "DINOv2 feature extraction complete"
    )
  }

  in_size <- 518L
  tensor <- .preprocess_nchw(
    mat,
    target_size = in_size,
    mean = c(0.485, 0.456, 0.406),
    std = c(0.229, 0.224, 0.225),
    letterbox = FALSE
  )

  dino_res <- run_dinov2_cpp(
    tensor_vec = tensor,
    in_w = in_size,
    in_h = in_size,
    patch_size = as.integer(patch_size),
    return_pca = isTRUE(return_pca),
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  pca_img <- NULL
  if (isTRUE(return_pca)) {
    wp <- dino_res$wp
    hp <- dino_res$hp
    r_mat <- dino_res$pca_r
    g_mat <- dino_res$pca_g
    b_mat <- dino_res$pca_b

    if (isTRUE(interpolate)) {
      r_mat <- .bilinear_resize_2d(r_mat, orig_w, orig_h)
      g_mat <- .bilinear_resize_2d(g_mat, orig_w, orig_h)
      b_mat <- .bilinear_resize_2d(b_mat, orig_w, orig_h)
      out_w <- orig_w
      out_h <- orig_h
    } else {
      out_w <- wp
      out_h <- hp
    }

    arr <- array(0.0, dim = c(out_w, out_h, 3L))
    arr[, , 1] <- r_mat
    arr[, , 2] <- g_mat
    arr[, , 3] <- b_mat
    pca_img <- as_image(arr)

    if (isTRUE(plot)) {
      plot(pca_img, ...)
    }
  }

  out <- list(
    cls_token = dino_res$cls_token,
    pca_image = pca_img,
    pca_r = dino_res$pca_r,
    pca_g = dino_res$pca_g,
    pca_b = dino_res$pca_b,
    num_patches = dino_res$num_patches,
    embed_dim = dino_res$embed_dim
  )
  invisible(out)
}

#' Universal Object Detection with YOLO26
#'
#' Performs real-time object detection using YOLO26 (or custom YOLO ONNX models)
#' with bounding box regression and End-to-End detection.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n"`.
#'   Options include `"yolo26n"`, `"yolo26s"`, `"yolo26m"`, `"yolo26l"`, `"yolo26x"`.
#' @param conf_threshold Minimum confidence score for candidate boxes (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param labels Optional character vector of class names (defaults to standard 80 COCO classes).
#' @param col Optional color or palette for bounding boxes.
#' @param rainbow Logical. If `TRUE`, assigns a distinct, unique color to each detected object.
#'   If `FALSE` (the default), objects of the same class share the exact same color. Defaults to `FALSE`.
#' @param lwd Line width for bounding boxes (default 2).
#' @param cex Size factor for label text (default `1.0`). Use smaller values (e.g., `cex = 0.2` or `0.5`)
#'   to reduce text size for dense scenes or many objects.
#' @param pad Padding multiplier controlling the area/background badge size around label text (default `1.0`).
#'   Can be a single value or a length-2 vector `c(pad_x, pad_y)`. Use `pad = 0` or smaller values for a compact area.
#' @param show_text Logical. Whether to display text labels on bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in labels (default `TRUE`).
#' @param show_class Logical. Whether to display class names in labels (default `TRUE`).
#' @param show_id Logical. Whether to display object IDs (`#1`, `#2`, ...) in labels (default `FALSE`).
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param return_features Logical. If `TRUE`, extracts deep feature representations for detected
#'   bounding boxes and the global image. Attaches `"features"` (an \eqn{N \times D} matrix of per-object
#'   deep embeddings, 512-D for CLIP or 384-D for DINOv2) and `"feature_map"` (dense spatial 3-component PCA false-color `Image`)
#'   as attributes to the returned object. Defaults to `FALSE`.
#' @param feature_model Foundation model used for feature representation when `return_features = TRUE`.
#'   Options are `"clip-vit-b32"` (default, OpenAI CLIP ViT-B/32 producing 512-D semantic embeddings) or
#'   `"dinov2"` (Meta DINOv2-S ViT producing 384-D self-supervised patch embeddings).
#' @param plot_features Logical. If `TRUE` and `return_features = TRUE`, visualizes the dense spatial PCA feature map. Defaults to `FALSE`.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the original image with bounding boxes.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A data frame containing detected bounding boxes (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`),
#'   with `summary` and `counts` stored as attributes. If `return_features = TRUE`, `features` (an \eqn{N \times D} matrix)
#'   and `feature_map` (dense spatial PCA `Image`) are also attached as attributes.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("objects.png")
#'   boxes <- image_detect_dl(img)
#'   # Rainbow coloring per box:
#'   boxes_rb <- image_detect_dl(img, rainbow = TRUE)
#' }
image_detect_dl <- function(img,
                            model = "yolo26n",
                            conf_threshold = 0.25,
                            iou_threshold = 0.45,
                            labels = NULL,
                            col = NULL,
                            rainbow = TRUE,
                            lwd = 2,
                            cex = 1.0,
                            pad = 1.0,
                            show_text = TRUE,
                            show_conf = TRUE,
                            show_class = TRUE,
                            show_id = FALSE,
                            badge = TRUE,
                            return_features = FALSE,
                            feature_model = c("yolo", "clip-vit-b32", "dinov2"),
                            plot_features = FALSE,
                            threads = 0,
                            engine = c("gpu", "cpu"),
                            device_id = -1,
                            verbose = TRUE,
                            plot = TRUE,
                            dir = pliman_model_dir(),
                            ...) {
  dots <- list(...)
  if ("label_size" %in% names(dots)) cex <- dots$label_size
  if ("text_size" %in% names(dots)) cex <- dots$text_size
  if ("font_scale" %in% names(dots)) cex <- dots$font_scale
  if ("cex_scale" %in% names(dots)) cex <- dots$cex_scale
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  if (length(cex) > 1L) {
    if (missing(pad)) pad <- cex[2]
    cex <- cex[1]
  }

  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_detect_dl(x, model = model, conf_threshold = conf_threshold, iou_threshold = iou_threshold,
                      labels = labels, col = col, rainbow = rainbow, lwd = lwd,
                      cex = cex, pad = pad, show_text = show_text, show_conf = show_conf,
                      show_class = show_class, show_id = show_id, badge = badge,
                      return_features = return_features, feature_model = feature_model, plot_features = plot_features,
                      threads = threads, engine = engine, device_id = device_id, verbose = verbose, plot = FALSE,
                      dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  is_url <- is.character(model) && grepl("^https?://", model[1])
  model_str <- if (is_url) {
    model[1]
  } else if (is.character(model)) {
    .resolve_model_name(model[1])
  } else {
    "yolo26n"
  }

  if (!is_url && grepl("[-_]cls", tolower(model_str))) {
    cli::cli_abort(c(
      "Model {.val {model[1]}} is an image classification model, not an object detector.",
      "i" = "Please use {.fn image_classify_dl} instead."
    ))
  }

  # Route all YOLO-World variants to image_detect_world
  is_world_model <- grepl("world", tolower(model_str))
  if (is_world_model) {
    world_classes <- if (!is.null(labels)) labels else if (!is.null(list(...)$classes)) list(...)$classes else c("object")
    # Determine which YOLO-World variant to use
    world_model <- model_str  # e.g. "yolov8s-worldv2" or "yolov8l-worldv2"
    return(image_detect_world(
      img = mat,
      classes = world_classes,
      conf_threshold = conf_threshold,
      iou_threshold = iou_threshold,
      rainbow = rainbow,
      col = col,
      lwd = lwd,
      cex = cex,
      pad = pad,
      show_text = show_text,
      show_conf = show_conf,
      show_class = show_class,
      show_id = show_id,
      badge = badge,
      return_features = return_features,
      feature_model = feature_model,
      plot_features = plot_features,
      plot = plot,
      verbose = verbose,
      dir = dir,
      model = world_model,
      ...
    ))
  }

  if (isTRUE(verbose)) {
    pb <- cli::cli_process_start(
      "Running YOLO object detection [{toupper(engine)}]...",
      on_exit = "failed"
    )
  }

  res <- .run_yolo(
    mat = mat,
    model = model_str,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    labels = labels,
    threads = threads,
    engine = engine,
    device_id = device_id,
    need_masks = FALSE,
    return_features = return_features,
    dir = dir
  )

  if (isTRUE(verbose)) {
    sum_txt <- if (!is.null(res$summary) && nzchar(res$summary) && res$summary != "0 objects") {
      paste0("pliman detected {.bold ", res$summary, "} in the image.")
    } else {
      "No objects detected in the image."
    }
    cli::cli_process_done(id = pb, msg_done = sum_txt)
  }

  df_boxes <- res$boxes

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    num_inst <- nrow(df_boxes)
    if (num_inst > 0) {
      palette_colors <- .get_class_palette(
        labels = df_boxes$label,
        n = num_inst,
        user_col = col,
        default_col = if (!is.null(col)) col else "#00CC66",
        rainbow = rainbow
      )
      .plot_yolo_bboxes(
        boxes = df_boxes,
        palette_colors = palette_colors,
        lwd = lwd,
        cex = cex,
        pad = pad,
        show_text = show_text,
        show_conf = show_conf,
        show_class = show_class,
        show_id = show_id,
        badge = badge
      )
    }
  }

  df_boxes <- .attach_features(
    obj = df_boxes,
    mat = mat,
    boxes = df_boxes,
    return_features = return_features,
    feature_model = feature_model,
    plot_features = plot_features,
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir,
    verbose = verbose,
    yolo_features = res
  )

  attr(df_boxes, "summary") <- res$summary
  attr(df_boxes, "counts") <- res$counts
  attr(df_boxes, "image") <- as_image(mat)
  if (isTRUE(return_features)) {
    attr(df_boxes, "heatmap") <- res$heatmap
    attr(df_boxes, "feature_energy") <- res$feature_energy
    attr(df_boxes, "feature_map") <- res$feature_map
    attr(df_boxes, "embeddings") <- res$embeddings
  }
  invisible(df_boxes)
}

#' List Connected Video Cameras
#'
#' Scans and lists all video capture devices (webcams, USB cameras, capture cards)
#' available on the system, indicating device IDs and native camera names.
#'
#' @param details Logical. If `TRUE`, also displays supported resolutions, framerates,
#'   and pixel formats for each device. Default is `FALSE`.
#' @return A data frame containing camera `id` and friendly `name`, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   # List available cameras
#'   list_cameras()
#'
#'   # View detailed resolutions and formats
#'   list_cameras(details = TRUE)
#' }
list_cameras <- function(details = FALSE) {
  if (.Platform$OS.type == "windows" && has_native_camera_cpp()) {
    cams <- tryCatch(list_native_cameras_cpp(), error = function(e) list(id = integer(0), name = character(0)))
    df <- data.frame(
      id = as.integer(cams$id),
      name = as.character(cams$name),
      stringsAsFactors = FALSE
    )
    if (nrow(df) == 0) {
      cli::cli_alert_warning("No video cameras detected on this system.")
      return(invisible(df))
    }

    cli::cli_h1("Connected Video Cameras")
    for (i in seq_len(nrow(df))) {
      cid <- df$id[i]
      cname <- df$name[i]
      formats <- tryCatch(list_camera_formats_cpp(cid), error = function(e) NULL)
      res_info <- ""
      if (!is.null(formats) && nrow(formats) > 0) {
        unique_res <- unique(paste0(formats$width, "x", formats$height))
        top_res <- unique_res[1]
        res_info <- sprintf(" [Native: %s | %d modes]", top_res, length(unique_res))
      }
      cli::cli_alert_info("Camera [{cid}]: {.strong {cname}}{res_info}")
      if (isTRUE(details) && !is.null(formats) && nrow(formats) > 0) {
        fmt_summary <- unique(data.frame(
          resolution = paste0(formats$width, "x", formats$height),
          fps = formats$fps,
          format = formats$format,
          stringsAsFactors = FALSE
        ))
        for (j in seq_len(min(8L, nrow(fmt_summary)))) {
          cli::cli_text("    - {fmt_summary$resolution[j]} @ {round(fmt_summary$fps[j])} FPS ({fmt_summary$format[j]})")
        }
        if (nrow(fmt_summary) > 8L) {
          cli::cli_text("    - ... ({nrow(fmt_summary) - 8L} more modes)")
        }
      }
    }
    cli::cli_rule()
    cli::cli_text("To use a camera in {.fun video_detect_dl}, pass:")
    for (i in seq_len(nrow(df))) {
      cli::cli_text("  {.code video_detect_dl(video = {df$id[i]})} or {.code video_detect_dl(video = '{df$name[i]}')}")
    }
    invisible(df)
  } else {
    cli::cli_alert_info("Camera listing is currently supported on Windows native drivers.")
    invisible(data.frame(id = integer(0), name = character(0)))
  }
}

#' Real-Time Video Object Detection with YOLO (ONNX Runtime)
#'
#' Performs real-time object detection and tracking on video files (.mp4, .avi, .mov)
#' or live camera feeds (webcam) using YOLO models via native C++ ONNX Runtime with
#' zero Python dependencies. When `return_data = TRUE`, it computes an extensive ecological
#' and temporal community profile (richness, abundance, Shannon diversity, Simpson index,
#' Pielou evenness, species accumulation curves, and class persistence rates).
#'
#' @details
#' \subsection{Ecological and Analytical Framework}{
#' In automated phenotyping, precision agriculture, ecological monitoring, and industrial inspection,
#' a video stream represents a continuous temporal transect composed of \eqn{F} sequential frames
#' recorded at a frame rate of \eqn{\text{FPS}}{FPS} frames per second. In each frame \eqn{t \in \{1, \dots, F\}},
#' deep neural network inference detects candidate bounding boxes belonging to \eqn{S_t} recognized classes.
#'
#' By aggregating detections over space and time, `video_detect_dl` translates raw computer vision
#' detections into standardized ecological metrics, quantifying instantaneous community state,
#' temporal persistence, class dominance, and cumulative discovery rates.
#' }
#'
#' \subsection{Abundance Metrics and Density Dynamics}{
#' Abundance characterizes the numerical density of detected objects:
#' \itemize{
#'   \item \strong{Instantaneous Abundance (\eqn{N_t})}: The total count of object instances detected in frame \eqn{t}:
#'     \deqn{N_t = \sum_{c=1}^S n_{t, c}}{N_t = \sum_{c=1}^S n_{t, c}}
#'     where \eqn{n_{t, c}} is the count of detected objects belonging to class \eqn{c} in frame \eqn{t}, and \eqn{S} is the total number of candidate classes.
#'   \item \strong{Total Video Abundance (\eqn{N})}: The cumulative sum of all detections identified across the processed footage:
#'     \deqn{N = \sum_{t=1}^F N_t = \sum_{c=1}^S n_c}{N = \sum_{t=1}^F N_t = \sum_{c=1}^S n_c}
#'     where \eqn{n_c = \sum_{t=1}^F n_{t, c}} is the total detection count for class \eqn{c}.
#'   \item \strong{Average Abundance per Frame (\eqn{\bar{N}})}: Quantifies the baseline object density per frame:
#'     \deqn{\bar{N} = \frac{N}{F}}{mean_N = N / F}
#'   \item \strong{Peak Abundance (\eqn{N_{\text{peak}}})}: The maximum instantaneous count observed in any single frame:
#'     \deqn{N_{\text{peak}} = \max_{1 \le t \le F} N_t}{N_peak = max(N_t)}
#'     The exact frame index \eqn{t_{\text{peak}} = \operatorname{argmax}_t N_t} and corresponding timestamp (in seconds)
#'     \eqn{T_{\text{peak}} = (t_{\text{peak}} - 1) / \text{FPS}} identify the moment of maximum congestion or cluster density.
#' }
#' }
#'
#' \subsection{Richness and Species Accumulation Dynamics}{
#' Richness quantifies class variety and categorical complexity over time:
#' \itemize{
#'   \item \strong{Instantaneous Richness (\eqn{S_t})}: The number of distinct classes co-occurring simultaneously in frame \eqn{t}:
#'     \deqn{S_t = \sum_{c=1}^S \mathbb{I}(n_{t, c} > 0)}{S_t = \sum_{c=1}^S I(n_{t, c} > 0)}
#'     where \eqn{\mathbb{I}(\cdot)} is the indicator function equal to 1 if class \eqn{c} has at least one detection in frame \eqn{t}, and 0 otherwise.
#'   \item \strong{Total Video Richness (\eqn{S})}: The total count of unique classes recognized across the entire video:
#'     \deqn{S = |\{c : n_c > 0\}|}{S = count of unique classes detected}
#'   \item \strong{Cumulative Richness / Species Accumulation Curve (Collector's Curve)}:
#'     \deqn{S_{\text{cum}}(t) = |\{c : \exists t' \le t \text{ such that } n_{t', c} > 0\}|}{S_cum(t) = count of distinct classes discovered up to frame t}
#'     The Collector's Curve plots the progressive discovery of novel classes over time. A plateau indicates that sampling effort
#'     (video duration) is sufficient to capture the full assemblage of classes present in the environment.
#'   \item \strong{Cumulative Detections (\eqn{N_{\text{cum}}(t)})}:
#'     \deqn{N_{\text{cum}}(t) = \sum_{t'=1}^t N_{t'}}{N_cum(t) = \sum_{t'=1}^t N_{t'}}
#'     Tracking the cumulative volume of detections provides insights into detection throughput and arrival rates.
#' }
#' }
#'
#' \subsection{Per-Class Ecological Metrics and Constancy}{
#' For each class \eqn{c}, `video_detect_dl` extracts detailed persistence and distribution metrics:
#' \itemize{
#'   \item \strong{Total Detections (\eqn{n_c})}: Total bounding box instances identified for class \eqn{c}.
#'   \item \strong{Relative Abundance (\eqn{p_c})}: The percentage contribution of class \eqn{c} to total video detections:
#'     \deqn{p_c = \frac{n_c}{N} \times 100\%}{p_c = (n_c / N) * 100\%}
#'   \item \strong{Frames Present (\eqn{F_c})}: The number of frames in which class \eqn{c} was detected at least once:
#'     \deqn{F_c = \sum_{t=1}^F \mathbb{I}(n_{t, c} > 0)}{F_c = \sum_{t=1}^F I(n_{t, c} > 0)}
#'   \item \strong{Occurrence Rate / Constancy (\eqn{O_c})}: The temporal persistence of class \eqn{c} across the video footage:
#'     \deqn{O_c = \frac{F_c}{F} \times 100\%}{O_c = (F_c / F) * 100\%}
#'     In ecological literature, species with \eqn{O_c \ge 50\%} are classified as constant, \eqn{25\% \le O_c < 50\%} as accessory, and \eqn{O_c < 25\%} as accidental.
#'   \item \strong{Mean Abundance When Present}: The expected number of individuals per frame given that the class is present:
#'     \deqn{\bar{n}_{c, \text{present}} = \frac{n_c}{\max(1, F_c)}}{mean_n_present = n_c / F_c}
#'   \item \strong{Mean Abundance Overall}: The unconditional expected number of individuals per frame across all processed frames:
#'     \deqn{\bar{n}_{c, \text{overall}} = \frac{n_c}{F}}{mean_n_overall = n_c / F}
#'   \item \strong{Maximum per Frame (\eqn{\max_t n_{t, c}})}: Peak instantaneous count of class \eqn{c}.
#'   \item \strong{Mean Confidence (\eqn{\bar{s}_c})}: The arithmetic mean of detection confidence scores for class \eqn{c}:
#'     \deqn{\bar{s}_c = \frac{1}{n_c} \sum_{i=1}^{n_c} \text{score}_i}{mean_conf = (1 / n_c) * \sum score_i}
#'   \item \strong{Temporal Span}: First observed timestamp (\eqn{t_{\text{first}, c}}) and last observed timestamp (\eqn{t_{\text{last}, c}}) in seconds.
#' }
#' }
#'
#' \subsection{Diversity, Dominance, and Evenness Indices}{
#' Summary diversity indices synthesize class variety and equitability into scalar coefficients:
#' \itemize{
#'   \item \strong{Shannon-Wiener Diversity Index (\eqn{H'})}:
#'     \deqn{H' = -\sum_{c=1}^S p_c \ln(p_c)}{H' = - \sum_{c=1}^S p_c * ln(p_c)}
#'     where \eqn{p_c = n_c / N}. \eqn{H'} quantifies the uncertainty in predicting the class identity of a randomly sampled detection.
#'     Higher values indicate higher diversity and a more equitable distribution across classes.
#'     In addition to the global video-level index, \eqn{H'_t} is computed dynamically for each frame in the `$temporal` data frame.
#'   \item \strong{Simpson's Diversity Index (\eqn{1 - D}) and Dominance (\eqn{D})}:
#'     \deqn{D = \sum_{c=1}^S p_c^2, \qquad 1 - D = 1 - \sum_{c=1}^S p_c^2}{D = \sum_{c=1}^S p_c^2, 1 - D = 1 - \sum_{c=1}^S p_c^2}
#'     \eqn{D} represents the probability that two randomly chosen detections belong to the same class (dominance),
#'     while \eqn{1 - D} (Gini-Simpson index) represents the probability that they belong to different classes.
#'     Values range from 0 (monoculture / single class dominating entirely) to \eqn{1 - 1/S} (maximal diversity).
#'   \item \strong{Pielou's Evenness Index (\eqn{J'})}:
#'     \deqn{J' = \frac{H'}{\ln(S)}}{J' = H' / ln(S)}
#'     for \eqn{S > 1} (with \eqn{J' = 1.0} if \eqn{S = 1}, and \eqn{J' = 0.0} if \eqn{S = 0}).
#'     \eqn{J'} normalizes Shannon diversity against its theoretical maximum (\eqn{H'_{\text{max}} = \ln(S)}),
#'     providing an equitability index bounded between 0 and 1. A value of 1.0 indicates that all classes are detected with equal frequency.
#'   \item \strong{Berger-Parker Dominance Index (\eqn{d})}:
#'     \deqn{d = \max_{1 \le c \le S} (p_c)}{d = max(p_c)}
#'     Quantifies the proportional dominance of the most abundant single class in the assemblage.
#' }
#' }
#'
#' \subsection{Multi-Object Tracking and Unique Entity Counting}{
#' When `track = TRUE`, a spatial-temporal Kalman filter tracking algorithm (SORT/ByteTrack)
#' associates detections across consecutive frames using Intersection over Union (IoU) cost matrices:
#' \itemize{
#'   \item \strong{Unique Objects (\eqn{U})}: Persistent IDs (\verb{#1}, \verb{#2}, ...) track physical entities over time,
#'     preventing double-counting when calculating true entity abundance:
#'     \deqn{U = |\{ \text{distinct persistent IDs} > 0 \}|}{U = count of distinct persistent IDs}
#'     For each class \eqn{c}, \eqn{U_c} records the unique tracked individuals assigned to that class.
#'   \item \strong{Cumulative Unique Objects Curve (\eqn{U_{\text{cum}}(t)})}:
#'     Tracks the progressive entry of novel unique physical entities into the monitored frame.
#'   \item \strong{Virtual Tripwire Counting (\code{count_line})}:
#'     When a virtual line is provided as \code{c(x1, y1, x2, y2)} in normalized coordinates \eqn{[0, 1]},
#'     the trajectory segments between consecutive centroids \eqn{\mathbf{p}_{t-1}} and \eqn{\mathbf{p}_t} of tracked objects
#'     are evaluated for intersection with the virtual line segment.
#'     Crossing events are registered precisely when the signed 2D cross-product changes orientation,
#'     cataloguing the event frame, timestamp, persistent ID, and class label in `$crossings`.
#' }
#' }
#'
#' \subsection{Visualization and S3 Methods}{
#' The returned object has class `c("pliman_video_detect", "list")` and provides built-in S3 methods:
#' \itemize{
#'   \item `plot(x, type = ...)`: Renders clean, publication-ready graphics. Supported `type` options:
#'     \itemize{
#'       \item `"temporal"` (default): Two-panel stacked time series displaying total abundance \eqn{N(t)} and richness \eqn{S(t)}.
#'       \item `"abundance"`: Multi-class abundance trajectories over time.
#'       \item `"richness"`: Step-function curve of class richness over time.
#'       \item `"cumulative"`: Two-panel plot of the Species Accumulation Curve (Collector's Curve) and cumulative detection volume.
#'       \item `"class"`: Horizontal barplot displaying total detections and relative abundance percentages for each class.
#'       \item `"all"`: 2x2 multi-panel overview combining abundance, richness, class distribution, and the Collector's Curve.
#'     }
#'   \item `summary(x)`: Formats and displays an overview table with video metadata, ecological metrics, and class distribution.
#'   \item `print(x)`: Displays an informative console summary with key video metrics and class rankings.
#'   \item `as.data.frame(x)`, `head(x)`, `tail(x)`, `dim(x)`, `nrow(x)`, `ncol(x)`, `x$column`, `x[i, j]`:
#'     Provide transparent backward compatibility with data frame operations on raw detections.
#' }
#' }
#'
#' @param video Input video source. Either a camera index (e.g. `0`, `1`), a camera name
#'   or substring (e.g. `"EMEET"`, `"HD User Facing"`), `"camera"` (default camera `0`),
#'   or a path to a video file (e.g., `"field_drone.mp4"`).
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n"`.
#'   Pre-trained options include `"yolo26n"`, `"yolo26s"`, `"yolov8n"`, `"yolov8s"`,
#'   `"yolov11n"`, etc., or a full path to a custom user-trained `.onnx` model.
#' @param conf_threshold Minimum confidence threshold for candidate detections (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param labels Optional character vector of class names. If `NULL` (default), labels are
#'   automatically detected from model metadata or fallback to standard classes.
#' @param classes Optional character vector of specific class names to filter and retain.
#'   If specified, only detections matching these classes will be kept.
#' @param detect Logical. Whether to run object detection model inference. If `FALSE`,
#'   the camera/video stream runs in pass-through mode without neural network inference, allowing
#'   fast real-time streaming and recording. Defaults to `TRUE`.
#' @param output Optional output file path (e.g. `"annotated_video.mp4"`) to save the processed video.
#' @param save_video Logical. Whether to save the processed video. Automatically `TRUE` if `output` is specified.
#' @param return_data Logical. Whether to compute and return the detection dataset and ecological summary statistics (default `TRUE`). If `FALSE`, returns `invisible(NULL)`.
#' @param track Logical. Whether to enable multi-object tracking and assign persistent IDs (`#1`, `#2`, etc.) across frames. Defaults to `FALSE`.
#' @param trail Logical. Whether to draw historical trajectory motion trail lines behind tracked objects when `track = TRUE`. Defaults to `TRUE`. Set to `FALSE` to maintain persistent IDs while hiding the trail.
#' @param cex Size factor for label text (default `1.0`). Use smaller values (e.g., `cex = 0.2` or `0.5`)
#'   to reduce text size for dense scenes or many objects.
#' @param pad Padding multiplier controlling the area/background badge size around label text (default `1.0`).
#'   Can be a single value or a length-2 vector `c(pad_x, pad_y)`. Use `pad = 0` or smaller values for a compact area.
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param text_size Numeric font scale for bounding box text labels and badges (alias for `cex`, default `1.0`).
#' @param show_text Logical. Whether to display text badges above bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in badges (default `TRUE`).
#' @param show_class Logical. Whether to display class label names in badges (default `TRUE`).
#' @param show_id Logical. Whether to display persistent tracking IDs (`#1`, `#2`) in badges when `track = TRUE` (default `FALSE`).
#' @param count_line Optional counting tripwire line specified as `c(x1, y1, x2, y2)` in normalized coordinates (`0` to `1`). Objects crossing this virtual line are counted.
#' @param count_line_label Optional label displayed along the counting line (default `NULL`, no text).
#' @param roi Optional Region of Interest coordinates `c(xmin, ymin, xmax, ymax)` in normalized coordinates (`0` to `1`).
#' @param roi_mode Either `"crop"` (infer only inside ROI) or `"filter"` (detect everywhere, retain only detections inside ROI).
#' @param roi_label Label displayed on the ROI box (default `"ROI / MONITORED ZONE"`).
#' @param hide_outside_roi Logical. Whether to remove / suppress bounding boxes and visual annotations for objects outside the ROI zone. Defaults to `TRUE` when `roi` is provided.
#' @param hide_counted Logical. Whether to remove / suppress bounding boxes for objects that have already crossed the virtual counting line. Defaults to `TRUE` when counting / tracking is active.
#' @param flash_counted Logical. Whether to show a brief neon green confirmation flash (~12 frames) when an object crosses the counting line before removing its bounding box. Defaults to `TRUE`.
#' @param fullscreen Logical. Whether to display the video feed in full screen mode (default `FALSE`).
#' @param window_size Optional numeric vector `c(width, height)` for the display window.
#' @param window_scale Scaling factor for the display window (default 1.0).
#' @param resolution Optional camera capture resolution as a numeric vector `c(width, height)`,
#'   e.g. `c(1920, 1080)` or `c(1280, 720)`. Default is `NULL` (uses native sensor resolution).
#' @param crop Optional frame cropping. Set to `"1:1"`, `"square"`, or `TRUE` for an automatic 1:1 centered square crop (e.g. 1080x1080 or 720x720), or `c(width, height)` for a custom centered rectangular crop. Default is `NULL`.
#' @param backend Camera and display backend: `"native"` (default, uses native Windows Media Foundation / Win32 GDI with zero external dependencies).
#' @param skeleton Logical. Whether to render the 17-keypoint anatomical skeleton when running pose estimation models (e.g. `"yolo26n-pose"`). Defaults to `TRUE`.
#' @param kpt_threshold Minimum confidence threshold for rendering individual keypoints and skeleton limbs (default 0.3).
#' @param kpt_radius Integer radius of keypoint joint markers in pixels (default 4).
#' @param bbox Logical. Whether to render bounding boxes around detected objects or people (default `TRUE`).
#' @param mask Logical. Whether to render semi-transparent instance segmentation masks when using segmentation models (default `TRUE`).
#' @param alpha Opacity/transparency of instance segmentation masks (default `0.4`, range `0` to `1`).
#' @param show Logical. Whether to display the video feed with bounding boxes in real time.
#'   Default is `interactive()`.
#' @param show_fps Logical. Whether to overlay a real-time FPS and detection counter badge (default `TRUE`).
#' @param hud_pos Screen position of the telemetry HUD card: `"top-right"` (default, avoids top-left ROI zone),
#'   `"top-left"`, `"bottom-left"`, `"bottom-right"`, `"top"`, or `"bottom"`.
#' @param hud_layout Layout style of the telemetry HUD: `"vertical"` (default modern multi-row glassmorphic card)
#'   or `"horizontal"` (sleek capsule pill bar).
#' @param infer_every Integer stride for neural network inference (default `1L`, running inference on every frame).
#'   Setting e.g. `infer_every = 2L` (or using alias `stride = 2` / `skip_frames = 1`) runs the heavy neural network
#'   on odd frames and smoothly propagates tracks and detections on intermediate frames, cutting compute time in half
#'   and maximizing real-time video display smoothness up to the camera's physical frame rate.
#' @param max_frames Optional maximum number of frames to process. Default is `NULL` (all frames).
#' @param fps Target framerate for video output. If `NULL` (default), matches the source video or camera.
#' @param col Optional color or palette for bounding boxes. If `NULL`, distinct colors are assigned.
#' @param lwd Line width for bounding boxes (default 2).
#' @param rainbow Logical. If `TRUE`, assigns distinct colors per detected instance.
#' @param engine Execution engine: `"gpu"` (DirectML, default) or `"cpu"`.
#' @param device_id GPU device ID (default -1 for auto).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param dir Directory where pre-trained models are stored (default `pliman_model_dir()`).
#' @param verbose Logical. Whether to show progress messages.
#' @param ... Additional arguments passed to graphics plotting functions or aliases
#'   (`camera`, `cam`, `device`, `infer`, `inference`, `save`, `return_df`, `return_detections`).
#' @return If `return_data = TRUE`, an object of class `c("pliman_video_detect", "list")` containing:
#' \itemize{
#'   \item `summary`: A `data.frame` of video-level metadata and ecological diversity metrics:
#'     `video_source`, `model`, `engine`, `total_frames`, `fps`, `duration_sec`,
#'     `total_detections` (\eqn{N}), `total_richness` (\eqn{S}), `unique_objects` (\eqn{U}),
#'     `avg_abundance_frame` (\eqn{\bar{N}}), `max_abundance_frame` (\eqn{N_{\text{peak}}}),
#'     `peak_frame` (\eqn{t_{\text{peak}}}), `peak_time_sec` (\eqn{T_{\text{peak}}}),
#'     `shannon_index` (\eqn{H'}), `simpson_index` (\eqn{1-D}), `pielou_evenness` (\eqn{J'}),
#'     and optionally `total_counted` (when `count_line` is used).
#'   \item `by_class`: A `data.frame` sorted by total detections detailing per-class metrics:
#'     `class`, `total_detections` (\eqn{n_c}), `relative_abundance` (\eqn{p_c}, \%),
#'     `unique_objects` (\eqn{U_c}), `n_frames_present` (\eqn{F_c}), `occurrence_rate` (\eqn{O_c}, \%),
#'     `mean_abundance_when_present` (\eqn{\bar{n}_{c, \text{present}}}), `mean_abundance_overall` (\eqn{\bar{n}_{c, \text{overall}}}),
#'     `max_per_frame`, `mean_confidence` (\eqn{\bar{s}_c}), `first_seen_sec` (\eqn{t_{\text{first}, c}}),
#'     `last_seen_sec` (\eqn{t_{\text{last}, c}}), and optionally `total_counted`.
#'   \item `unique_objects`: A `data.frame` aggregating every unique tracked object across all observed frames:
#'     `id`, `label`, `n_frames`, `first_frame`, `last_frame`, `duration_sec`, `mean_score`,
#'     and (when using instance segmentation models) detailed morphological shape features:
#'     `area` (in pixels), `area_sd`, `perimeter`, `radius_mean`, `length`, `width`, `circularity`, `eccentricity`, and `asp_ratio`.
#'   \item `temporal`: A `data.frame` of frame-by-frame time series at resolution \eqn{t \in \{1, \dots, F\}}:
#'     `frame`, `timestamp` (in seconds), instantaneous `abundance` (\eqn{N_t}), instantaneous `richness` (\eqn{S_t}),
#'     instantaneous `shannon` diversity (\eqn{H'_t}), and individual count columns for every detected class.
#'   \item `cumulative`: A `data.frame` of cumulative progression:
#'     `frame`, `timestamp`, `cumulative_detections` (\eqn{N_{\text{cum}}(t)}),
#'     `cumulative_richness` (Species Accumulation Curve / Collector's Curve, \eqn{S_{\text{cum}}(t)}),
#'     and optionally `cumulative_unique_objects` (\eqn{U_{\text{cum}}(t)}).
#'   \item `diversity`: A `list` containing global ecological diversity indices:
#'     `richness_S` (\eqn{S}), `total_abundance_N` (\eqn{N}), `shannon_H` (\eqn{H'}),
#'     `simpson_1_minus_D` (\eqn{1 - D}), `simpson_D` (\eqn{D}), `pielou_J` (\eqn{J'}),
#'     `berger_parker_d` (\eqn{d}), `mean_richness_per_frame`, and `mean_abundance_per_frame`.
#'   \item `detections`: A `data.frame` of all raw candidate detections across all frames:
#'     `frame`, `timestamp`, `id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`, and `class_id`
#'     (plus shape features when running segmentation models).
#'   \item `crossings`: (Optional) A `data.frame` of virtual tripwire crossing events when `count_line` is specified.
#'   \item `keypoints`: (Optional) A `data.frame` of 17 keypoint anatomical landmarks when running pose estimation models.
#'   \item `roi`, `roi_mode`, `count_line`, `output`: Additional metadata reflecting execution configuration.
#' }
#' If `return_data = FALSE`, returns `invisible(NULL)`.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'
#'   # 1. Video detection and ecological diversity analysis
#'   dets <- video_detect_dl(
#'     video = "field_drone.mp4",
#'     model = "yolo26n",
#'     conf_threshold = 0.5,
#'     return_data = TRUE,
#'     output = "annotated_field.mp4"
#'   )
#'
#'   # Summary table of ecological metrics
#'   summary(dets)
#'
#'   # Plot community temporal dynamics (abundance & richness)
#'   plot(dets, type = "temporal")
#'
#'   # Plot Species Accumulation Curve (Collector's Curve)
#'   plot(dets, type = "cumulative")
#'
#'   # Plot multi-class abundance trajectories
#'   plot(dets, type = "abundance")
#'
#'   # Multi-panel 2x2 diagnostic overview
#'   plot(dets, type = "all")
#'
#'   # Access per-class breakdown
#'   dets$by_class
#'
#'   # 2. Multi-object tracking with virtual tripwire line counting
#'   tracked_dets <- video_detect_dl(
#'     video = "conveyor_belt.mp4",
#'     model = "seed_weights.onnx",
#'     track = TRUE,
#'     count_line = c(0.5, 0, 0.5, 1), # vertical counting tripwire
#'     return_data = TRUE
#'   )
#' }
video_detect_dl <- function(video = 0,
                            model = "yolo26n",
                            conf_threshold = 0.25,
                            iou_threshold = 0.45,
                            labels = NULL,
                            classes = NULL,
                            detect = TRUE,
                            output = NULL,
                            save_video = !is.null(output),
                            return_data = TRUE,
                            track = FALSE,
                            trail = TRUE,
                            cex = 1.0,
                            pad = 1.0,
                            badge = TRUE,
                            show_text = TRUE,
                            show_conf = TRUE,
                            show_class = TRUE,
                            show_id = FALSE,
                            text_size = 1.0,
                            count_line = NULL,
                            count_line_label = NULL,
                            roi = NULL,
                            roi_mode = c("crop", "filter"),
                            roi_label = "ROI / MONITORED ZONE",
                            hide_outside_roi = TRUE,
                            hide_counted = TRUE,
                            flash_counted = TRUE,
                            fullscreen = FALSE,
                            window_size = NULL,
                            window_scale = 1.0,
                            resolution = NULL,
                            crop = NULL,
                            backend = c("native", "auto"),
                            skeleton = TRUE,
                            kpt_threshold = 0.3,
                            kpt_radius = 4,
                            bbox = TRUE,
                            mask = TRUE,
                            alpha = 0.4,
                            show = NULL,
                            show_fps = TRUE,
                            hud_pos = c("top-right", "top-left", "bottom-left", "bottom-right", "top", "bottom"),
                            hud_layout = c("vertical", "horizontal"),
                            infer_every = 1L,
                            max_frames = NULL,
                            fps = NULL,
                            col = NULL,
                            lwd = 2,
                            rainbow = TRUE,
                            engine = c("gpu", "cpu"),
                            device_id = -1,
                            threads = 0,
                            dir = pliman_model_dir(),
                            verbose = TRUE,
                            ...) {
  if (!missing(cex) && missing(text_size)) {
    text_size <- cex
  } else if (!missing(text_size) && missing(cex)) {
    cex <- text_size
  }
  dots <- list(...)
  if ("camera" %in% names(dots)) video <- dots$camera
  if ("cam" %in% names(dots)) video <- dots$cam
  if ("device" %in% names(dots)) video <- dots$device
  if ("square" %in% names(dots) && isTRUE(dots$square)) crop <- "1:1"
  if ("aspect_ratio" %in% names(dots) && dots$aspect_ratio %in% c("1:1", "square")) crop <- "1:1"
  if ("trail" %in% names(dots)) trail <- dots$trail
  if ("draw_trail" %in% names(dots)) trail <- dots$draw_trail
  if ("remove_counted" %in% names(dots)) hide_counted <- isTRUE(dots$remove_counted)
  if ("bbox_counted" %in% names(dots)) hide_counted <- !isTRUE(dots$bbox_counted)
  if ("show_counted" %in% names(dots)) hide_counted <- !isTRUE(dots$show_counted)
  if ("bbox_outside_roi" %in% names(dots)) hide_outside_roi <- !isTRUE(dots$bbox_outside_roi)
  if ("show_outside_roi" %in% names(dots)) hide_outside_roi <- !isTRUE(dots$show_outside_roi)
  if ("remove_outside_roi" %in% names(dots)) hide_outside_roi <- isTRUE(dots$remove_outside_roi)
  if ("track_trail" %in% names(dots)) trail <- dots$track_trail
  if ("cex" %in% names(dots)) { text_size <- dots$cex; cex <- dots$cex }
  if ("text_size" %in% names(dots)) { text_size <- dots$text_size; cex <- dots$text_size }
  if ("font_scale" %in% names(dots)) { text_size <- dots$font_scale; cex <- dots$font_scale }
  if ("label_size" %in% names(dots)) { text_size <- dots$label_size; cex <- dots$label_size }
  if ("cex_scale" %in% names(dots)) { text_size <- dots$cex_scale; cex <- dots$cex_scale }
  if ("pad" %in% names(dots)) pad <- dots$pad
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("badge" %in% names(dots)) badge <- dots$badge
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("labels" %in% names(dots) && is.logical(dots$labels)) show_class <- isTRUE(dots$labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("scores" %in% names(dots) && is.logical(dots$scores)) show_conf <- isTRUE(dots$scores)
  if ("conf" %in% names(dots) && is.logical(dots$conf)) show_conf <- isTRUE(dots$conf)
  if ("show_text" %in% names(dots)) show_text <- isTRUE(dots$show_text)
  if ("show_id" %in% names(dots)) show_id <- isTRUE(dots$show_id)
  if ("mask" %in% names(dots)) mask <- dots$mask
  if ("masks" %in% names(dots)) mask <- dots$masks
  if ("draw_mask" %in% names(dots)) mask <- dots$draw_mask
  if ("draw_masks" %in% names(dots)) mask <- dots$draw_masks
  if ("alpha" %in% names(dots)) alpha <- dots$alpha
  if ("mask_alpha" %in% names(dots)) alpha <- dots$mask_alpha
  if ("infer" %in% names(dots)) detect <- dots$infer
  if ("inference" %in% names(dots)) detect <- dots$inference
  if ("save" %in% names(dots)) save_video <- dots$save
  if ("return_df" %in% names(dots)) return_data <- dots$return_df
  if ("return_detections" %in% names(dots)) return_data <- dots$return_detections

  detect <- isTRUE(detect)
  save_video <- isTRUE(save_video)
  return_data <- isTRUE(return_data)
  trail <- isTRUE(trail)
  mask <- isTRUE(mask)
  alpha <- as.numeric(alpha)
  if (is.na(alpha) || alpha < 0) alpha <- 0.4
  if (alpha > 1.0) alpha <- 1.0
  src_fps <- if (!is.null(fps) && is.numeric(fps) && fps > 0) as.numeric(fps) else 30.0

  if ("hud_pos" %in% names(dots)) hud_pos <- dots$hud_pos
  if ("pos" %in% names(dots)) hud_pos <- dots$pos
  if ("position" %in% names(dots)) hud_pos <- dots$position
  if ("hud_position" %in% names(dots)) hud_pos <- dots$hud_position
  if ("hud_layout" %in% names(dots)) hud_layout <- dots$hud_layout
  if ("layout" %in% names(dots)) hud_layout <- dots$layout
  if ("hud_style" %in% names(dots)) hud_layout <- dots$hud_style
  if ("style" %in% names(dots)) hud_layout <- dots$style

  valid_pos <- c("top-right", "top-left", "bottom-left", "bottom-right", "top", "bottom")
  hud_pos <- if (is.character(hud_pos) && length(hud_pos) > 0) {
    p_match <- match.arg(tolower(hud_pos[1]), valid_pos)
    p_match
  } else "top-right"

  valid_layout <- c("vertical", "horizontal")
  hud_layout <- if (is.character(hud_layout) && length(hud_layout) > 0) {
    l_match <- match.arg(tolower(hud_layout[1]), valid_layout)
    l_match
  } else "vertical"

  if ("infer_every" %in% names(dots)) infer_every <- dots$infer_every
  if ("infer_each" %in% names(dots)) infer_every <- dots$infer_each
  if ("infer_step" %in% names(dots)) infer_every <- dots$infer_step
  if ("step" %in% names(dots)) infer_every <- dots$step
  if ("stride" %in% names(dots)) infer_every <- dots$stride
  if ("skip_frames" %in% names(dots)) infer_every <- as.integer(dots$skip_frames) + 1L
  infer_every <- max(1L, as.integer(infer_every))

  text_size_num <- suppressWarnings(as.numeric(text_size))
  if (is.na(text_size_num)) text_size_num <- 1.0
  if (isFALSE(show_text) || text_size_num <= 0.05) {
    show_text <- FALSE
    f_scale <- 0.0
  } else {
    f_scale <- text_size_num
  }
  pad_num <- if (is.numeric(pad) && length(pad) > 0) pad[1] else 1.0
  pad_num <- max(0.0, as.numeric(pad_num))
  badge <- isTRUE(badge)

  # For small text (<= 0.3), hide redundant class name and confidence score by default
  # unless the user explicitly requested them in arguments / dots
  if (!"show_conf" %in% names(dots) && !"show_scores" %in% names(dots) && !"scores" %in% names(dots) && text_size_num <= 0.3) {
    show_conf <- FALSE
  }
  if (!"show_class" %in% names(dots) && !"show_labels" %in% names(dots) && text_size_num <= 0.3 && (isTRUE(track) || !is.null(count_line))) {
    show_class <- FALSE
  }
  track_max_dist <- if ("max_dist" %in% names(dots)) as.numeric(dots$max_dist) else 150.0
  track_min_iou <- if ("min_iou" %in% names(dots)) as.numeric(dots$min_iou) else 0.25
  track_max_lost <- if ("max_lost" %in% names(dots)) as.integer(dots$max_lost) else 15L
  track_max_history <- if ("max_history" %in% names(dots)) as.integer(dots$max_history) else 30L
  roi_mode <- match.arg(roi_mode)
  backend <- match.arg(backend)
  if (backend == "auto") backend <- "native"

  # Safely wipe any dangling ghost progress bars from previously interrupted sessions
  try({
    cli_ns <- asNamespace("cli")
    if (exists("clienv", envir = cli_ns, inherits = FALSE)) {
      assign("progress", list(), envir = cli_ns$clienv)
    }
  }, silent = TRUE)

  on.exit({
    try(cli::cli_progress_cleanup(), silent = TRUE)
    try({
      cli_ns <- asNamespace("cli")
      if (exists("clienv", envir = cli_ns, inherits = FALSE)) {
        assign("progress", list(), envir = cli_ns$clienv)
      }
    }, silent = TRUE)
  }, add = TRUE)

  if (save_video && is.null(output)) {
    output <- if (detect) "annotated_video.mp4" else "recorded_video.mp4"
  }
  if (!save_video) {
    output <- NULL
  }

  out_video_file <- if (save_video && !is.null(output)) {
    normalizePath(output, winslash = "/", mustWork = FALSE)
  } else {
    ""
  }

  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (detect) {
    # Model resolution: custom ONNX file or pre-trained download
    if (file.exists(model[1])) {
      model_file <- normalizePath(model[1], winslash = "/")
    } else if (file.exists(file.path(dir, model[1]))) {
      model_file <- normalizePath(file.path(dir, model[1]), winslash = "/")
    } else if (file.exists(file.path(dir, paste0(model[1], ".onnx")))) {
      model_file <- normalizePath(file.path(dir, paste0(model[1], ".onnx")), winslash = "/")
    } else if (file.exists(file.path("D:/Desktop/models", model[1]))) {
      model_file <- normalizePath(file.path("D:/Desktop/models", model[1]), winslash = "/")
    } else if (file.exists(file.path("D:/Desktop/models", paste0(model[1], ".onnx")))) {
      model_file <- normalizePath(file.path("D:/Desktop/models", paste0(model[1], ".onnx")), winslash = "/")
    } else {
      model_file <- pliman_download_model(model = model[1], dir = dir)
    }

    # Label resolution: custom vector -> ONNX embedded metadata -> COCO fallback
    if (!is.null(labels)) {
      class_names <- labels
    } else {
      model_classes <- .get_onnx_classes(model_file)
      if (!is.null(model_classes) && length(model_classes) > 0L) {
        class_names <- model_classes
      } else {
        class_names <- .coco_classes
      }
    }
  } else {
    model_file <- "None"
    class_names <- character(0)
  }

  has_video_ext <- is.character(video) && tolower(tools::file_ext(video[1])) %in% c("mp4", "avi", "mov", "mkv", "wmv", "flv", "webm", "m4v")
  is_camera <- is.numeric(video) || (is.character(video) && !has_video_ext && (
    tolower(video[1]) %in% c("0", "1", "2", "3", "4", "5", "camera", "live", "webcam") ||
    !file.exists(video[1])
  ))
  if (is.null(show)) {
    show <- if (is_camera) interactive() else FALSE
  }
  if (!is_camera && isTRUE(show) && isTRUE(verbose)) {
    cli::cli_alert_info("Interactive frame display enabled (rendering plots to R graphics device reduces throughput).")
  }

  # Initialize object tracking state
  track_state <- list(
    next_id = 1L,
    total_count = 0L,
    tracks = list(),
    counts_by_class = list(),
    crossing_events = list()
  )
  total_count_val <- 0L
  in_zone_val <- 0L

  # ============================================================================
  # MODE A: LIVE WEBCAM STREAMING (NATIVE C++ MEDIA FOUNDATION)
  # ============================================================================
  if (is_camera) {
    if (.Platform$OS.type != "windows" || !has_native_camera_cpp()) {
      cli::cli_abort(c(
        "x" = "Live camera streaming is currently supported natively on Windows (via Windows Media Foundation).",
        "i" = "To analyze video on Linux/macOS, please pass a video file path (e.g., {.path 'sample.mp4'})."
      ))
    }

    if (isTRUE(verbose)) {
      if (detect) {
        cli::cli_alert_info("Starting real-time camera inference with model {.val {basename(model_file)}} [{toupper(engine)} | native C++ engine]")
      } else {
        cli::cli_alert_info("Starting real-time camera stream [pass-through mode, inference disabled | native C++ engine]")
      }
      if (isTRUE(fullscreen)) {
        cli::cli_alert_info("Display mode: Fullscreen enabled.")
      }
      if (!is.null(roi)) {
        cli::cli_alert_info("Region of Interest (ROI) active: mode = {.val {roi_mode}}.")
      }
      if (isTRUE(track) || !is.null(count_line)) {
        cli::cli_alert_info("Real-time object tracking and counting tripwire enabled.")
      }
      cli::cli_alert_info("Press {.kbd ESC} or close the video window to stop.")
    }

    if (save_video) {
      out_video_file <- normalizePath(output, winslash = "/", mustWork = FALSE)
      cam_dir <- file.path(tempdir(), paste0("cam_yolo_", as.integer(stats::runif(1, 1000, 9999))))
      if (!dir.exists(cam_dir)) dir.create(cam_dir, recursive = TRUE, showWarnings = FALSE)
      on.exit(unlink(cam_dir, recursive = TRUE, force = TRUE), add = TRUE)
      recorded_frames <- list()
    }

    frame_records <- list()
    frame_keypoints <- list()
    frame_count <- 0L
    t_start <- Sys.time()
    last_t <- Sys.time()
    fps_smooth <- 0.0

    cam_id <- 0L
    cam_name <- "Default Camera"
      avail_cams <- tryCatch(list_native_cameras_cpp(), error = function(e) list(id = integer(0), name = character(0)))
      if (length(avail_cams$id) > 0) {
        if (is.numeric(video)) {
          idx_req <- as.integer(video[1])
          match_idx <- which(avail_cams$id == idx_req)
          if (length(match_idx) > 0) {
            cam_id <- avail_cams$id[match_idx[1]]
            cam_name <- avail_cams$name[match_idx[1]]
          } else {
            cli::cli_alert_warning("Camera index {idx_req} not found. Available: {paste0('[', avail_cams$id, '] ', avail_cams$name, collapse = ', ')}. Defaulting to [{avail_cams$id[1]}].")
            cam_id <- avail_cams$id[1]
            cam_name <- avail_cams$name[1]
          }
        } else if (is.character(video)) {
          v_num <- suppressWarnings(as.integer(video[1]))
          if (!is.na(v_num) && v_num %in% avail_cams$id) {
            cam_id <- v_num
            cam_name <- avail_cams$name[which(avail_cams$id == v_num)[1]]
          } else {
            matched <- grep(video[1], avail_cams$name, ignore.case = TRUE)
            if (length(matched) > 0) {
              cam_id <- avail_cams$id[matched[1]]
              cam_name <- avail_cams$name[matched[1]]
            } else {
              cli::cli_alert_warning("No camera matching '{video[1]}' found. Available: {paste0('[', avail_cams$id, '] ', avail_cams$name, collapse = ', ')}. Defaulting to [{avail_cams$id[1]}].")
              cam_id <- avail_cams$id[1]
              cam_name <- avail_cams$name[1]
            }
          }
        }
      }

      target_w <- 0L
      target_h <- 0L
      if (!is.null(resolution) && length(resolution) >= 2) {
        target_w <- as.integer(resolution[1])
        target_h <- as.integer(resolution[2])
      }

      cam <- tryCatch({
        open_native_camera_cpp(cam_id, target_w, target_h)
      }, error = function(e) {
        cli::cli_abort(c(
          "x" = "Failed to open native camera {cam_id} ({cam_name}): {e$message}",
          "i" = "Please check if your webcam is connected and not in use by another application."
        ))
      })
      on.exit(close_native_camera_cpp(cam), add = TRUE)

      cam_dims <- get_camera_dims_cpp(cam)
      cam_w <- cam_dims$width
      cam_h <- cam_dims$height

      # Detect if square aspect ratio or custom crop requested
      is_square <- isTRUE(crop) ||
                   (is.character(crop) && tolower(crop[1]) %in% c("1:1", "square")) ||
                   (!is.null(resolution) && length(resolution) == 2 && resolution[1] == resolution[2]) ||
                   (!is.null(resolution) && length(resolution) == 1 && is.numeric(resolution))

      crop_x <- 0L
      crop_y <- 0L
      crop_w <- 0L
      crop_h <- 0L
      orig_w <- cam_w
      orig_h <- cam_h

      if (is_square) {
        req_size <- if (!is.null(resolution)) as.integer(resolution[1]) else min(cam_w, cam_h)
        sq_size <- min(req_size, cam_w, cam_h)
        crop_w <- as.integer(sq_size)
        crop_h <- as.integer(sq_size)
        crop_x <- as.integer(round((cam_w - sq_size) / 2))
        crop_y <- as.integer(round((cam_h - sq_size) / 2))
        orig_w <- crop_w
        orig_h <- crop_h
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Square frame active: 1:1 aspect ratio ({orig_w}x{orig_h}) centered from {cam_w}x{cam_h} camera feed.")
        }
      } else if (is.numeric(crop) && length(crop) >= 2) {
        req_w <- min(as.integer(crop[1]), cam_w)
        req_h <- min(as.integer(crop[2]), cam_h)
        crop_w <- req_w
        crop_h <- req_h
        crop_x <- as.integer(round((cam_w - req_w) / 2))
        crop_y <- as.integer(round((cam_h - req_h) / 2))
        orig_w <- crop_w
        orig_h <- crop_h
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Custom cropped frame: {orig_w}x{orig_h} centered from {cam_w}x{cam_h} camera feed.")
        }
      } else {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Active camera: [{cam_id}] {.strong {cam_name}} ({orig_w}x{orig_h} | {round(orig_w/orig_h, 2)}:1 FOV)")
        }
      }

      bm <- raw(3L * orig_w * orig_h)
      attr(bm, "dim") <- c(3L, orig_w, orig_h)

      win <- NULL
      if (isTRUE(show)) {
        w_w <- if (!is.null(window_size)) window_size[1] else round(orig_w * window_scale)
        w_h <- if (!is.null(window_size)) window_size[2] else round(orig_h * window_scale)
        win_title <- if (detect) {
          sprintf("pliman - YOLO [%s | %s] (ESC to exit, Space to pause, F for fullscreen)", toupper(engine), basename(model_file))
        } else {
          sprintf("pliman - Camera Stream [%s] (ESC to exit, Space to pause, F for fullscreen)", cam_name)
        }
        win <- create_native_window_cpp(title = win_title, width = as.integer(w_w), height = as.integer(w_h), fullscreen = isTRUE(fullscreen))
        on.exit(close_native_window_cpp(win), add = TRUE)
      }

      last_df_b <- NULL
      last_raw_res <- NULL
      last_has_kpts <- FALSE
      last_has_masks <- FALSE
      last_kpt_mat_pass <- NULL
      last_track_ids_pass <- NULL
      last_flash_pass <- NULL
      last_history_x_pass <- NULL
      last_history_y_pass <- NULL
      last_total_count_val <- 0L
      last_in_zone_val <- 0L
      last_use_crop <- FALSE

      while (TRUE) {
        if (!is.null(max_frames) && frame_count >= max_frames) {
          break
        }

        has_frame <- grab_native_frame_cpp(cam, bm, crop_x1 = crop_x, crop_y1 = crop_y, crop_w = crop_w, crop_h = crop_h)
        if (!has_frame) {
          break
        }

        t_cur <- Sys.time()
        dt <- as.numeric(t_cur - last_t, units = "secs")
        last_t <- t_cur
        cur_fps <- if (dt > 0) (1.0 / dt) else 30.0
        fps_smooth <- if (fps_smooth == 0) cur_fps else (0.85 * fps_smooth + 0.15 * cur_fps)

        frame_count <- frame_count + 1L
        t_elapsed <- as.numeric(t_cur - t_start, units = "secs")

        # ROI coordinates calculation
        roi_coords <- NULL
        if (!is.null(roi) && length(roi) >= 4) {
          rx1 <- if (roi[1] <= 1.0) round(roi[1] * orig_w) else round(roi[1])
          ry1 <- if (roi[2] <= 1.0) round(roi[2] * orig_h) else round(roi[2])
          rx2 <- if (roi[3] <= 1.0) round(roi[3] * orig_w) else round(roi[3])
          ry2 <- if (roi[4] <= 1.0) round(roi[4] * orig_h) else round(roi[4])
          rx1 <- max(0L, min(orig_w - 2L, as.integer(rx1)))
          ry1 <- max(0L, min(orig_h - 2L, as.integer(ry1)))
          rx2 <- max(rx1 + 2L, min(orig_w, as.integer(rx2)))
          ry2 <- max(ry1 + 2L, min(orig_h, as.integer(ry2)))
          roi_coords <- c(rx1, ry1, rx2, ry2)
        }

        # Virtual count line coordinates calculation
        line_coords <- NULL
        if (!is.null(count_line)) {
          if (length(count_line) == 1) {
            pos_x <- if (count_line <= 1.0) round(count_line * orig_w) else round(count_line)
            line_coords <- c(as.numeric(pos_x), 0.0, as.numeric(pos_x), as.numeric(orig_h))
          } else if (length(count_line) >= 4) {
            lx1 <- if (count_line[1] <= 1.0) round(count_line[1] * orig_w) else round(count_line[1])
            ly1 <- if (count_line[2] <= 1.0) round(count_line[2] * orig_h) else round(count_line[2])
            lx2 <- if (count_line[3] <= 1.0) round(count_line[3] * orig_w) else round(count_line[3])
            ly2 <- if (count_line[4] <= 1.0) round(count_line[4] * orig_h) else round(count_line[4])
            line_coords <- c(as.numeric(lx1), as.numeric(ly1), as.numeric(lx2), as.numeric(ly2))
          }
        }

        df_b <- data.frame(
          frame = integer(0), timestamp = numeric(0), id = integer(0),
          xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
          label = character(0), score = numeric(0), class_id = integer(0),
          stringsAsFactors = FALSE
        )
        has_kpts <- FALSE
        has_masks <- FALSE
        kpt_mat_pass <- NULL
        track_ids_pass <- NULL
        flash_pass <- NULL
        history_x_pass <- NULL
        history_y_pass <- NULL

        use_crop <- FALSE
        should_infer <- (frame_count == 1L || (frame_count %% infer_every == 0L))

        if (detect && should_infer) {
          use_crop <- (roi_mode == "crop" && !is.null(roi_coords))
          if (use_crop) {
            margin_x <- round(0.08 * orig_w)
            margin_y <- round(0.08 * orig_h)
            crop_x1 <- max(0L, as.integer(roi_coords[1] - margin_x))
            crop_y1 <- max(0L, as.integer(roi_coords[2] - margin_y))
            crop_x2 <- min(orig_w, as.integer(roi_coords[3] + margin_x))
            crop_y2 <- min(orig_h, as.integer(roi_coords[4] + margin_y))
            crop_coords <- c(crop_x1, crop_y1, crop_x2, crop_y2)
            bm_feed <- crop_bgr_cpp(bm, crop_x1, crop_y1, crop_x2, crop_y2)
          } else {
            crop_coords <- c(0L, 0L, orig_w, orig_h)
            bm_feed <- bm
          }
          feed_dims <- dim(bm_feed)
          feed_w <- feed_dims[2]
          feed_h <- feed_dims[3]

          tensor <- tryCatch(
            preprocess_yolo_cpp(bm_feed, target_size = 640L),
            error = function(e) preprocess_yolo_cpp(bm, target_size = 640L)
          )

          raw_res <- run_yolo_cpp(
            tensor_vec = tensor,
            orig_w = feed_w,
            orig_h = feed_h,
            conf_threshold = conf_threshold,
            iou_threshold = iou_threshold,
            model_path = model_file,
            lib_path = lib_path,
            num_threads = as.integer(threads),
            use_gpu = use_gpu,
            device_id = as.integer(device_id),
            need_masks = isTRUE(mask)
          )

          has_masks <- isTRUE(mask) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && any(raw_res$labels > 0L)
          num_boxes_orig <- nrow(raw_res$boxes)
          orig_idx <- seq_len(num_boxes_orig)

          # Shift box coordinates back to full frame if crop was used
          if (use_crop && nrow(raw_res$boxes) > 0) {
            raw_res$boxes[, 1] <- raw_res$boxes[, 1] + crop_coords[1]
            raw_res$boxes[, 2] <- raw_res$boxes[, 2] + crop_coords[2]
            raw_res$boxes[, 3] <- raw_res$boxes[, 3] + crop_coords[1]
            raw_res$boxes[, 4] <- raw_res$boxes[, 4] + crop_coords[2]
            if (!is.null(raw_res$keypoints) && ncol(raw_res$keypoints) >= 51) {
              raw_res$keypoints[, seq(1, 51, by = 3)] <- raw_res$keypoints[, seq(1, 51, by = 3)] + crop_coords[1]
              raw_res$keypoints[, seq(2, 51, by = 3)] <- raw_res$keypoints[, seq(2, 51, by = 3)] + crop_coords[2]
            }
          } else if (roi_mode == "filter" && is.null(count_line) && !is.null(roi_coords) && nrow(raw_res$boxes) > 0) {
            cxs <- (raw_res$boxes[, 1] + raw_res$boxes[, 3]) / 2
            cys <- (raw_res$boxes[, 2] + raw_res$boxes[, 4]) / 2
            keep_roi <- which(cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            raw_res$boxes <- raw_res$boxes[keep_roi, , drop = FALSE]
            raw_res$scores <- raw_res$scores[keep_roi]
            raw_res$class_ids <- raw_res$class_ids[keep_roi]
            orig_idx <- orig_idx[keep_roi]
            if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
              raw_res$keypoints <- raw_res$keypoints[keep_roi, , drop = FALSE]
            }
          }

          num_boxes <- nrow(raw_res$boxes)
          if (num_boxes > 0) {
            c_ids <- raw_res$class_ids
            lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))
            df_b <- data.frame(
              frame = frame_count,
              timestamp = round(t_elapsed, 3),
              id = seq_len(num_boxes),
              xmin = raw_res$boxes[, 1],
              ymin = raw_res$boxes[, 2],
              xmax = raw_res$boxes[, 3],
              ymax = raw_res$boxes[, 4],
              label = lbls,
              score = round(raw_res$scores, 4),
              class_id = c_ids,
              orig_idx = orig_idx,
              stringsAsFactors = FALSE
            )
            if (!is.null(classes)) {
              keep_cls <- which(df_b$label %in% classes)
              df_b <- df_b[keep_cls, , drop = FALSE]
              if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
                raw_res$keypoints <- raw_res$keypoints[keep_cls, , drop = FALSE]
              }
            }

            if ((isTRUE(track) || !is.null(count_line))) {
              if (nrow(df_b) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = df_b$xmin,
                  ymin = df_b$ymin,
                  xmax = df_b$xmax,
                  ymax = df_b$ymax,
                  labels = df_b$label,
                  scores = df_b$score,
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history
                )
                track_state <- track_res$state
                df_b$id <- track_res$track_ids
                df_b$counted <- track_res$counted
                if (!is.null(track_res$in_zone_vec)) {
                  df_b$in_zone <- track_res$in_zone_vec
                } else if (!is.null(roi_coords)) {
                  cxs <- (df_b$xmin + df_b$xmax) / 2
                  cys <- (df_b$ymin + df_b$ymax) / 2
                  df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                }
                track_ids_pass <- track_res$track_ids
                flash_pass <- track_res$flash
                history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
                total_count_val <- track_res$total_count
                in_zone_val <- track_res$in_zone
              } else if (!is.null(track_state$tracks) && length(track_state$tracks) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = numeric(0),
                  ymin = numeric(0),
                  xmax = numeric(0),
                  ymax = numeric(0),
                  labels = character(0),
                  scores = numeric(0),
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history
                )
                track_state <- track_res$state
                total_count_val <- track_res$total_count
                in_zone_val <- track_res$in_zone
              }
            }

            if (!("in_zone" %in% names(df_b)) && !is.null(roi_coords) && nrow(df_b) > 0) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
            }

            has_masks <- isTRUE(mask) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && any(raw_res$labels > 0L)
            if (has_masks && nrow(df_b) > 0) {
              exact_area <- tabulate(raw_res$labels, nbins = num_boxes_orig)
              conts <- extract_contours_cpp(raw_res$labels)
              meas <- poly_measures_minimal_cpp(conts)
              m_idx <- df_b$orig_idx
              df_b$area <- as.numeric(exact_area[m_idx])
              df_b$perimeter <- round(as.numeric(meas$perimeter[m_idx]), 1)
              df_b$radius_mean <- round((as.numeric(meas$maj_axis[m_idx]) + as.numeric(meas$min_axis[m_idx])) / 2, 2)
              df_b$length <- round(as.numeric(meas$length[m_idx]), 1)
              df_b$width <- round(as.numeric(meas$width[m_idx]), 1)
              df_b$circularity <- round(as.numeric(meas$circularity_norm[m_idx]), 4)
              df_b$eccentricity <- round(as.numeric(meas$eccentricity[m_idx]), 4)
              df_b$asp_ratio <- round(as.numeric(meas$length[m_idx]) / pmax(0.001, as.numeric(meas$width[m_idx])), 3)
            }

            has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && nrow(raw_res$keypoints) == nrow(df_b)
            if (has_kpts) {
              kpt_mat_pass <- raw_res$keypoints
            }
          }

          last_df_b <- df_b
          last_raw_res <- raw_res
          last_has_kpts <- has_kpts
          last_has_masks <- has_masks
          last_kpt_mat_pass <- kpt_mat_pass
          last_track_ids_pass <- track_ids_pass
          last_flash_pass <- flash_pass
          last_history_x_pass <- history_x_pass
          last_history_y_pass <- history_y_pass
          last_total_count_val <- total_count_val
          last_in_zone_val <- in_zone_val
          last_use_crop <- use_crop
        } else if (detect && !should_infer && !is.null(last_df_b)) {
          df_b <- last_df_b
          raw_res <- last_raw_res
          has_kpts <- last_has_kpts
          has_masks <- last_has_masks
          kpt_mat_pass <- last_kpt_mat_pass
          track_ids_pass <- last_track_ids_pass
          flash_pass <- last_flash_pass
          history_x_pass <- last_history_x_pass
          history_y_pass <- last_history_y_pass
          total_count_val <- last_total_count_val
          in_zone_val <- last_in_zone_val
          use_crop <- last_use_crop
          if (nrow(df_b) > 0) {
            df_b$frame <- frame_count
            df_b$timestamp <- round(t_elapsed, 3)
          }
        }

        palette_colors <- if (detect && nrow(df_b) > 0) {
          .get_class_palette(
            labels = df_b$label,
            n = nrow(df_b),
            user_col = col,
            default_col = if (!is.null(col)) col else "#00CC66",
            rainbow = rainbow,
            ids = if ("id" %in% names(df_b)) df_b$id else NULL
          )
        } else {
          character(0)
        }

        if (detect && return_data) {
          if (nrow(df_b) > 0) df_b$color <- palette_colors
          frame_records[[frame_count]] <- df_b
          if (has_kpts) {
            kpts_list_frame <- lapply(seq_len(nrow(raw_res$keypoints)), function(k) {
              vals <- raw_res$keypoints[k, ]
              data.frame(
                frame = frame_count,
                id = df_b$id[k],
                keypoint = .coco_keypoints,
                x = vals[seq(1, 51, by = 3)],
                y = vals[seq(2, 51, by = 3)],
                conf = round(vals[seq(3, 51, by = 3)], 4),
                stringsAsFactors = FALSE
              )
            })
            frame_keypoints[[frame_count]] <- do.call(rbind, kpts_list_frame)
          }
        }

        hud_msg <- if (isTRUE(show_fps)) {
          if (isTRUE(track) || !is.null(count_line)) {
            sprintf("%.1f FPS | IN ZONE: %d | COUNT: %d", fps_smooth, in_zone_val, total_count_val)
          } else if (detect) {
            sprintf("%.1f FPS | %d Dets", fps_smooth, nrow(df_b))
          } else {
            sprintf("%.1f FPS", fps_smooth)
          }
        } else ""

        if ((detect && nrow(df_b) > 0) || nzchar(hud_msg) || !is.null(roi_coords) || !is.null(line_coords)) {
          draw_yolo_detections_bgr_cpp(
            bm = bm,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            labels = if (detect && nrow(df_b) > 0) df_b$label else character(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            colors = palette_colors,
            lwd = as.integer(lwd),
            hud_text = hud_msg,
            font_scale = as.numeric(f_scale),
            keypoints = kpt_mat_pass,
            kpt_threshold = as.numeric(kpt_threshold),
            kpt_radius = as.integer(kpt_radius),
            draw_skeleton = isTRUE(skeleton),
            draw_boxes = isTRUE(bbox),
            track_ids = track_ids_pass,
            flash = flash_pass,
            history_x = history_x_pass,
            history_y = history_y_pass,
            roi = roi_coords,
            roi_label = roi_label,
            count_line = line_coords,
            count_line_label = if (!is.null(count_line_label)) as.character(count_line_label) else "",
            counted = if ("counted" %in% names(df_b)) df_b$counted else NULL,
            hide_outside_roi = isTRUE(hide_outside_roi),
            hide_counted = isTRUE(hide_counted),
            flash_counted = isTRUE(flash_counted),
            show_text = isTRUE(show_text),
            show_conf = isTRUE(show_conf),
            show_class = isTRUE(show_class),
            show_id = isTRUE(show_id),
            badge = badge,
            pad_scale = pad_num,
            mask_labels = if (isTRUE(mask) && has_masks) raw_res$labels else NULL,
            mask_alpha = as.numeric(alpha),
            mask_offset_x = if (use_crop) as.integer(crop_coords[1]) else 0L,
            mask_offset_y = if (use_crop) as.integer(crop_coords[2]) else 0L,
            draw_masks = isTRUE(mask),
            mask_ids = if (has_masks && nrow(df_b) > 0) as.integer(df_b$orig_idx) else NULL,
            hud_pos = as.character(hud_pos),
            hud_layout = as.character(hud_layout)
          )
        }

        if (save_video) {
          f_rec_path <- file.path(cam_dir, sprintf("rec_%06d.jpg", frame_count))
          save_frame_bgr_cpp(bm, f_rec_path, quality = 90L)
          recorded_frames[[frame_count]] <- f_rec_path
        }

        if (isTRUE(show)) {
          show_native_frame_cpp(win, bm, preserve_aspect = TRUE)
          status <- poll_native_window_cpp(win)
          if (status == 1L) {
            break
          }
          while (status == 2L) {
            Sys.sleep(0.05)
            status <- poll_native_window_cpp(win)
            if (status == 1L) break
          }
          if (status == 1L) break
        }
      }

    all_dets <- if (length(frame_records) > 0) do.call(rbind, frame_records) else data.frame()
    all_kpts <- if (length(frame_keypoints) > 0) do.call(rbind, frame_keypoints) else NULL
    total_time <- as.numeric(Sys.time() - t_start, units = "secs")
    avg_fps <- if (total_time > 0 && frame_count > 0) (frame_count / total_time) else fps_smooth

    if (save_video && length(recorded_frames) > 0) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Encoding recorded camera session to {.path {out_video_file}}...")
      }
      rec_fps <- if (!is.null(fps)) fps else (if (avg_fps > 0) round(avg_fps, 1) else 24.0)
      av::av_encode_video(
        input = unlist(recorded_frames),
        output = out_video_file,
        framerate = rec_fps,
        verbose = FALSE
      )
      if (isTRUE(verbose)) cli::cli_alert_success("Camera video saved: {.path {out_video_file}}")
    }

  # ============================================================================
  # MODE B: VIDEO FILE PROCESSING (AV / FFMPEG)
  # ============================================================================
  } else {
    if (!file.exists(video)) {
      cli::cli_abort("Video file not found: {.path {video}}")
    }
    ext <- tolower(tools::file_ext(video))
    if (ext %in% c("jpg", "jpeg", "png", "bmp", "webp", "tif", "tiff")) {
      cli::cli_abort(c(
        "x" = "File {.path {video}} is a static image, not a video file.",
        "i" = "For static image detection or pose estimation, please use {.code image_detect_dl()} or {.code image_pose_dl()}."
      ))
    }
    if (!requireNamespace("av", quietly = TRUE)) {
      cli::cli_abort(c(
        "!" = "Package {.pkg av} is required for video file inference.",
        "i" = "Please install it with: {.code install.packages('av')}"
      ))
    }

    v_info <- av::av_video_info(video)
    src_w <- v_info$video$width
    src_h <- v_info$video$height
    orig_fps <- v_info$video$framerate
    if (is.null(orig_fps) || is.na(orig_fps) || orig_fps <= 0) orig_fps <- 30.0
    src_fps <- if (!is.null(fps)) fps else orig_fps
    if (is.null(src_fps) || is.na(src_fps) || src_fps <= 0) src_fps <- 30.0

    # Destination directory for decoded frames
    run_dir <- file.path(tempdir(), paste0("vid_yolo_", as.integer(stats::runif(1, 1000, 9999))))
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)

    if (isTRUE(verbose)) {
      cli::cli_progress_step("Decoding video frames with FFmpeg...", msg_done = "Video frames decoded")
    }
    filter_fps <- if (length(fps) && !is.null(fps)) paste0("fps=fps=", fps) else NULL
    vfilter <- if (length(filter_fps)) filter_fps else "null"
    if (!is.null(max_frames) && is.numeric(max_frames) && max_frames > 0) {
      lim_n <- as.integer(max_frames) + 2L
      sel_filter <- paste0("select=between(n\\,0\\,", lim_n, ")")
      vfilter <- if (vfilter == "null") sel_filter else paste0(vfilter, ",", sel_filter)
    }
    if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    output_pattern <- file.path(run_dir, "image_%06d.bmp")

    av::av_encode_video(
      input = video,
      output = output_pattern,
      framerate = orig_fps,
      codec = "bmp",
      vfilter = vfilter,
      verbose = FALSE
    )
    if (isTRUE(verbose)) cli::cli_progress_done()
    frame_paths <- sort(list.files(run_dir, pattern = "^image_\\d+\\.bmp$", full.names = TRUE))
    num_frames <- length(frame_paths)
    if (!is.null(fps) && !is.null(v_info$video$duration) && v_info$video$duration > 0) {
      expected_max <- ceiling(v_info$video$duration * fps) + 5
      if (num_frames > expected_max) {
        keep_idx <- round(seq(1, num_frames, length.out = ceiling(v_info$video$duration * fps)))
        keep_idx <- unique(pmin(num_frames, pmax(1L, keep_idx)))
        frame_paths <- frame_paths[keep_idx]
        num_frames <- length(frame_paths)
      }
    }
    if (!is.null(max_frames) && max_frames < num_frames) {
      frame_paths <- frame_paths[seq_len(max_frames)]
      num_frames <- length(frame_paths)
    }

    if (isTRUE(verbose)) {
      if (detect) {
        cli::cli_alert_info("Running YOLO inference on {num_frames} frame(s) [{toupper(engine)} | model: {.val {basename(model_file)}}]")
      } else {
        cli::cli_alert_info("Processing {num_frames} frame(s) [inference disabled]")
      }
      if (is.null(fps) && isTRUE(v_info$video$duration > 60)) {
        cli::cli_alert_info("Tip: For faster processing of long videos, you can use e.g. {.code fps = 5} or {.code fps = 10}.")
      }
    }

    frame_records <- vector("list", num_frames)
    frame_keypoints <- vector("list", num_frames)
    t_start <- Sys.time()

    # Pre-render dimensions must be even for H.264
    out_w <- src_w - (src_w %% 2)
    out_h <- src_h - (src_h %% 2)

    # ROI coordinates for Mode B
    roi_coords <- NULL
    if (!is.null(roi) && length(roi) >= 4) {
      rx1 <- if (roi[1] <= 1.0) round(roi[1] * src_w) else round(roi[1])
      ry1 <- if (roi[2] <= 1.0) round(roi[2] * src_h) else round(roi[2])
      rx2 <- if (roi[3] <= 1.0) round(roi[3] * src_w) else round(roi[3])
      ry2 <- if (roi[4] <= 1.0) round(roi[4] * src_h) else round(roi[4])
      rx1 <- max(0L, min(src_w - 2L, as.integer(rx1)))
      ry1 <- max(0L, min(src_h - 2L, as.integer(ry1)))
      rx2 <- max(rx1 + 2L, min(src_w, as.integer(rx2)))
      ry2 <- max(ry1 + 2L, min(src_h, as.integer(ry2)))
      roi_coords <- c(rx1, ry1, rx2, ry2)
    }

    line_coords <- NULL
    if (!is.null(count_line)) {
      if (length(count_line) == 1) {
        pos_x <- if (count_line <= 1.0) round(count_line * src_w) else round(count_line)
        line_coords <- c(as.numeric(pos_x), 0.0, as.numeric(pos_x), as.numeric(src_h))
      } else if (length(count_line) >= 4) {
        lx1 <- if (count_line[1] <= 1.0) round(count_line[1] * src_w) else round(count_line[1])
        ly1 <- if (count_line[2] <= 1.0) round(count_line[2] * src_h) else round(count_line[2])
        lx2 <- if (count_line[3] <= 1.0) round(count_line[3] * src_w) else round(count_line[3])
        ly2 <- if (count_line[4] <= 1.0) round(count_line[4] * src_h) else round(count_line[4])
        line_coords <- c(as.numeric(lx1), as.numeric(ly1), as.numeric(lx2), as.numeric(ly2))
      }
    }

    if (save_video) {
      out_video_file <- normalizePath(output, winslash = "/", mustWork = FALSE)
    }

    pb <- if (isTRUE(verbose)) {
      cli::cli_progress_bar(
        if (detect) "Detecting objects in video" else "Processing video",
        total = num_frames,
        extra = list(current_fps = "0.0"),
        format = "{cli::pb_spin} Processing [{cli::pb_current}/{cli::pb_total}] | {cli::pb_percent} | {cli::pb_extra$current_fps} FPS | ETA: {cli::pb_eta}"
      )
    } else NULL

    last_df_b <- NULL
    last_raw_res <- NULL
    last_has_kpts <- FALSE
    last_has_masks <- FALSE
    last_kpt_mat_pass <- NULL
    last_track_ids_pass <- NULL
    last_flash_pass <- NULL
    last_history_x_pass <- NULL
    last_history_y_pass <- NULL
    last_total_count_val <- 0L
    last_in_zone_val <- 0L
    last_use_crop <- FALSE
    last_crop_coords <- NULL

    process_loop <- function() {
      for (i in seq_len(num_frames)) {
        f_path <- frame_paths[i]
        has_kpts <- FALSE
        has_masks <- FALSE
        kpt_mat_pass <- NULL
        track_ids_pass <- NULL
        flash_pass <- NULL
        history_x_pass <- NULL
        history_y_pass <- NULL

        use_crop <- FALSE
        crop_coords <- c(0L, 0L, src_w, src_h)
        should_infer <- (i == 1L || (i %% infer_every == 0L))

        # Decode frame once into memory for both YOLO and drawing
        bm <- load_frame_bgr_cpp(f_path)

        if (detect && should_infer) {
          use_crop <- (roi_mode == "crop" && !is.null(roi_coords))
          if (use_crop) {
            margin_x <- round(0.08 * src_w)
            margin_y <- round(0.08 * src_h)
            crop_x1 <- max(0L, as.integer(roi_coords[1] - margin_x))
            crop_y1 <- max(0L, as.integer(roi_coords[2] - margin_y))
            crop_x2 <- min(src_w, as.integer(roi_coords[3] + margin_x))
            crop_y2 <- min(src_h, as.integer(roi_coords[4] + margin_y))
            crop_coords <- c(crop_x1, crop_y1, crop_x2, crop_y2)
            bm_feed <- crop_bgr_cpp(bm, crop_x1, crop_y1, crop_x2, crop_y2)
          } else {
            crop_coords <- c(0L, 0L, src_w, src_h)
            bm_feed <- bm
          }
          feed_dims <- dim(bm_feed)
          feed_w <- feed_dims[2]
          feed_h <- feed_dims[3]

          tensor <- tryCatch(
            preprocess_yolo_cpp(bm_feed, target_size = 640L),
            error = function(e) {
              f_mat <- image_data(image_import(f_path))
              preprocess_yolo_cpp(f_mat, target_size = 640L)
            }
          )

          raw_res <- run_yolo_cpp(
            tensor_vec = tensor,
            orig_w = feed_w,
            orig_h = feed_h,
            conf_threshold = conf_threshold,
            iou_threshold = iou_threshold,
            model_path = model_file,
            lib_path = lib_path,
            num_threads = as.integer(threads),
            use_gpu = use_gpu,
            device_id = as.integer(device_id),
            need_masks = isTRUE(mask)
          )

          has_masks <- isTRUE(mask) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && any(raw_res$labels > 0L)
          num_boxes_orig <- nrow(raw_res$boxes)
          orig_idx <- seq_len(num_boxes_orig)

          # Shift box coordinates back to full frame if crop was used
          if (use_crop && nrow(raw_res$boxes) > 0) {
            raw_res$boxes[, 1] <- raw_res$boxes[, 1] + crop_coords[1]
            raw_res$boxes[, 2] <- raw_res$boxes[, 2] + crop_coords[2]
            raw_res$boxes[, 3] <- raw_res$boxes[, 3] + crop_coords[1]
            raw_res$boxes[, 4] <- raw_res$boxes[, 4] + crop_coords[2]
            if (!is.null(raw_res$keypoints) && ncol(raw_res$keypoints) >= 51) {
              raw_res$keypoints[, seq(1, 51, by = 3)] <- raw_res$keypoints[, seq(1, 51, by = 3)] + crop_coords[1]
              raw_res$keypoints[, seq(2, 51, by = 3)] <- raw_res$keypoints[, seq(2, 51, by = 3)] + crop_coords[2]
            }
          } else if (roi_mode == "filter" && is.null(count_line) && !is.null(roi_coords) && nrow(raw_res$boxes) > 0) {
            cxs <- (raw_res$boxes[, 1] + raw_res$boxes[, 3]) / 2
            cys <- (raw_res$boxes[, 2] + raw_res$boxes[, 4]) / 2
            keep_roi <- which(cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            raw_res$boxes <- raw_res$boxes[keep_roi, , drop = FALSE]
            raw_res$scores <- raw_res$scores[keep_roi]
            raw_res$class_ids <- raw_res$class_ids[keep_roi]
            orig_idx <- orig_idx[keep_roi]
            if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
              raw_res$keypoints <- raw_res$keypoints[keep_roi, , drop = FALSE]
            }
          }

          num_boxes <- nrow(raw_res$boxes)
          t_sec <- round((i - 1) / src_fps, 3)

          if (num_boxes > 0) {
            c_ids <- raw_res$class_ids
            lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))
            df_b <- data.frame(
              frame = i,
              timestamp = t_sec,
              id = seq_len(num_boxes),
              xmin = raw_res$boxes[, 1],
              ymin = raw_res$boxes[, 2],
              xmax = raw_res$boxes[, 3],
              ymax = raw_res$boxes[, 4],
              label = lbls,
              score = round(raw_res$scores, 4),
              class_id = c_ids,
              orig_idx = orig_idx,
              stringsAsFactors = FALSE
            )
            if (!is.null(classes)) {
              keep_cls <- which(df_b$label %in% classes)
              df_b <- df_b[keep_cls, , drop = FALSE]
              if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
                raw_res$keypoints <- raw_res$keypoints[keep_cls, , drop = FALSE]
              }
            }

            if ((isTRUE(track) || !is.null(count_line))) {
              if (nrow(df_b) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = df_b$xmin,
                  ymin = df_b$ymin,
                  xmax = df_b$xmax,
                  ymax = df_b$ymax,
                  labels = df_b$label,
                  scores = df_b$score,
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history
                )
                track_state <<- track_res$state
                df_b$id <- track_res$track_ids
                df_b$counted <- track_res$counted
                if (!is.null(track_res$in_zone_vec)) {
                  df_b$in_zone <- track_res$in_zone_vec
                } else if (!is.null(roi_coords)) {
                  cxs <- (df_b$xmin + df_b$xmax) / 2
                  cys <- (df_b$ymin + df_b$ymax) / 2
                  df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                }
                track_ids_pass <- track_res$track_ids
                flash_pass <- track_res$flash
                history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
                total_count_val <<- track_res$total_count
                in_zone_val <<- track_res$in_zone
              } else if (!is.null(track_state$tracks) && length(track_state$tracks) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = numeric(0),
                  ymin = numeric(0),
                  xmax = numeric(0),
                  ymax = numeric(0),
                  labels = character(0),
                  scores = numeric(0),
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history
                )
                track_state <<- track_res$state
                total_count_val <<- track_res$total_count
                in_zone_val <<- track_res$in_zone
              }
            }

            if (!("in_zone" %in% names(df_b)) && !is.null(roi_coords) && nrow(df_b) > 0) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
            }

            has_masks <- isTRUE(mask) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && any(raw_res$labels > 0L)
            if (has_masks && nrow(df_b) > 0) {
              exact_area <- tabulate(raw_res$labels, nbins = num_boxes_orig)
              conts <- extract_contours_cpp(raw_res$labels)
              meas <- poly_measures_minimal_cpp(conts)
              m_idx <- df_b$orig_idx
              df_b$area <- as.numeric(exact_area[m_idx])
              df_b$perimeter <- round(as.numeric(meas$perimeter[m_idx]), 1)
              df_b$radius_mean <- round((as.numeric(meas$maj_axis[m_idx]) + as.numeric(meas$min_axis[m_idx])) / 2, 2)
              df_b$length <- round(as.numeric(meas$length[m_idx]), 1)
              df_b$width <- round(as.numeric(meas$width[m_idx]), 1)
              df_b$circularity <- round(as.numeric(meas$circularity_norm[m_idx]), 4)
              df_b$eccentricity <- round(as.numeric(meas$eccentricity[m_idx]), 4)
              df_b$asp_ratio <- round(as.numeric(meas$length[m_idx]) / pmax(0.001, as.numeric(meas$width[m_idx])), 3)
            }

            has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && nrow(raw_res$keypoints) == nrow(df_b)
            if (has_kpts) {
              kpt_mat_pass <- raw_res$keypoints
            }
          } else {
            df_b <- data.frame(
              frame = integer(0), timestamp = numeric(0), id = integer(0),
              xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
              label = character(0), score = numeric(0), class_id = integer(0),
              stringsAsFactors = FALSE
            )
          }

          last_df_b <<- df_b
          last_raw_res <<- raw_res
          last_has_kpts <<- has_kpts
          last_has_masks <<- has_masks
          last_kpt_mat_pass <<- kpt_mat_pass
          last_track_ids_pass <<- track_ids_pass
          last_flash_pass <<- flash_pass
          last_history_x_pass <<- history_x_pass
          last_history_y_pass <<- history_y_pass
          last_total_count_val <<- total_count_val
          last_in_zone_val <<- in_zone_val
          last_use_crop <<- use_crop
          last_crop_coords <<- crop_coords
        } else if (detect && !should_infer && !is.null(last_df_b)) {
          df_b <- last_df_b
          raw_res <- last_raw_res
          has_kpts <- last_has_kpts
          has_masks <- last_has_masks
          kpt_mat_pass <- last_kpt_mat_pass
          track_ids_pass <- last_track_ids_pass
          flash_pass <- last_flash_pass
          history_x_pass <- last_history_x_pass
          history_y_pass <- last_history_y_pass
          total_count_val <<- last_total_count_val
          in_zone_val <<- last_in_zone_val
          use_crop <- last_use_crop
          crop_coords <- last_crop_coords

          # If tracking is active, propagate motion across skipped frames using track velocities
          if ((isTRUE(track) || !is.null(count_line)) && nrow(df_b) > 0 && !is.null(track_state$tracks) && length(track_state$tracks) > 0) {
            for (ti in seq_along(track_state$tracks)) {
              trk <- track_state$tracks[[ti]]
              tid <- trk$id
              match_row <- which(df_b$id == tid)
              if (length(match_row) > 0) {
                vx <- if (!is.null(trk$vx)) trk$vx else 0.0
                vy <- if (!is.null(trk$vy)) trk$vy else 0.0
                if (abs(vx) > 0.1 || abs(vy) > 0.1) {
                  df_b$xmin[match_row] <- df_b$xmin[match_row] + vx
                  df_b$xmax[match_row] <- df_b$xmax[match_row] + vx
                  df_b$ymin[match_row] <- df_b$ymin[match_row] + vy
                  df_b$ymax[match_row] <- df_b$ymax[match_row] + vy
                }
              }
            }
            last_df_b <<- df_b
          }

          if (nrow(df_b) > 0) {
            df_b$frame <- i
            df_b$timestamp <- round((i - 1) / src_fps, 3)
          }
        } else {
          df_b <- data.frame(
            frame = integer(0), timestamp = numeric(0), id = integer(0),
            xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
            label = character(0), score = numeric(0), class_id = integer(0),
            stringsAsFactors = FALSE
          )
        }

        palette_colors <- if (detect && nrow(df_b) > 0) {
          .get_class_palette(
            labels = df_b$label,
            n = nrow(df_b),
            user_col = col,
            default_col = if (!is.null(col)) col else "#00CC66",
            rainbow = rainbow,
            ids = if ("id" %in% names(df_b)) df_b$id else NULL
          )
        } else {
          character(0)
        }

        if (detect && return_data) {
          if (nrow(df_b) > 0) df_b$color <- palette_colors
          frame_records[[i]] <<- df_b
          if (has_kpts) {
            kpts_list_frame <- lapply(seq_len(nrow(df_b)), function(k) {
              vals <- raw_res$keypoints[k, ]
              data.frame(
                frame = i,
                id = df_b$id[k],
                keypoint = .coco_keypoints,
                x = vals[seq(1, 51, by = 3)],
                y = vals[seq(2, 51, by = 3)],
                conf = round(vals[seq(3, 51, by = 3)], 4),
                stringsAsFactors = FALSE
              )
            })
            frame_keypoints[[i]] <<- do.call(rbind, kpts_list_frame)
          }
        }

        # Render frame with bounding boxes and annotations
        if (save_video && !isTRUE(show)) {
          elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
          cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
          hud_txt <- if (isTRUE(show_fps)) {
            if (isTRUE(track) || !is.null(count_line)) {
              sprintf("Frame %d/%d | %.1f FPS | IN ZONE: %d | COUNT: %d", i, num_frames, cur_fps_est, in_zone_val, total_count_val)
            } else if (detect) {
              sprintf("Frame %d/%d | %.1f FPS | %d Dets", i, num_frames, cur_fps_est, nrow(df_b))
            } else {
              sprintf("Frame %d/%d | %.1f FPS", i, num_frames, cur_fps_est)
            }
          } else ""

          draw_yolo_detections_bgr_cpp(
            bm = bm,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            labels = if (detect && nrow(df_b) > 0) df_b$label else character(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            colors = palette_colors,
            lwd = as.integer(lwd),
            hud_text = hud_txt,
            font_scale = as.numeric(f_scale),
            keypoints = kpt_mat_pass,
            kpt_threshold = as.numeric(kpt_threshold),
            kpt_radius = as.integer(kpt_radius),
            draw_skeleton = isTRUE(skeleton),
            draw_boxes = isTRUE(bbox),
            track_ids = track_ids_pass,
            flash = flash_pass,
            history_x = history_x_pass,
            history_y = history_y_pass,
            roi = roi_coords,
            roi_label = roi_label,
            count_line = line_coords,
            count_line_label = if (!is.null(count_line_label)) as.character(count_line_label) else "",
            counted = if ("counted" %in% names(df_b)) df_b$counted else NULL,
            hide_outside_roi = isTRUE(hide_outside_roi),
            hide_counted = isTRUE(hide_counted),
            flash_counted = isTRUE(flash_counted),
            show_text = isTRUE(show_text),
            show_conf = isTRUE(show_conf),
            show_class = isTRUE(show_class),
            show_id = isTRUE(show_id),
            badge = badge,
            pad_scale = pad_num,
            mask_labels = if (isTRUE(mask) && has_masks) raw_res$labels else NULL,
            mask_alpha = as.numeric(alpha),
            mask_offset_x = if (use_crop) as.integer(crop_coords[1]) else 0L,
            mask_offset_y = if (use_crop) as.integer(crop_coords[2]) else 0L,
            draw_masks = isTRUE(mask),
            mask_ids = if (has_masks && nrow(df_b) > 0) as.integer(df_b$orig_idx) else NULL,
            hud_pos = as.character(hud_pos),
            hud_layout = as.character(hud_layout)
          )
          save_frame_bgr_cpp(bm, f_path, quality = 90L)
        } else if (isTRUE(show)) {
          f_mat <- image_data(image_import(f_path))
          plot(as_image(f_mat), ...)
          if (!is.null(roi_coords)) {
            graphics::rect(roi_coords[1], roi_coords[2], roi_coords[3], roi_coords[4],
                           col = grDevices::adjustcolor("#FA8072", alpha.f = 0.3),
                           border = "#FA8072", lwd = 1)
            graphics::text(roi_coords[1] + 5, roi_coords[2] - 5, roi_label, col = "#FA8072", cex = 0.7, font = 2, adj = c(0, 0))
          }
          if (!is.null(line_coords)) {
            graphics::lines(c(line_coords[1], line_coords[3]), c(line_coords[2], line_coords[4]), col = "#00E5FF", lwd = 3)
            if (!is.null(count_line_label) && nzchar(as.character(count_line_label))) {
              is_vert <- abs(line_coords[4] - line_coords[2]) >= abs(line_coords[3] - line_coords[1])
              if (is_vert) {
                graphics::text(line_coords[1] - 5, min(line_coords[2], line_coords[4]) + 15, count_line_label, col = "#00E5FF", cex = 0.7, font = 2, adj = c(1, 0))
              } else {
                graphics::text(min(line_coords[1], line_coords[3]) + 5, line_coords[2] - 5, count_line_label, col = "#00E5FF", cex = 0.7, font = 2, adj = c(0, 0))
              }
            }
          }
          if (detect && nrow(df_b) > 0) {
            # Filter bounding boxes to plot based on ROI and counted status
            keep_plot <- rep(TRUE, nrow(df_b))
            if (isTRUE(hide_outside_roi) && !is.null(roi_coords)) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              keep_plot <- keep_plot & (cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            }
            if (isTRUE(hide_counted) && "counted" %in% names(df_b)) {
              is_flashing <- if (!is.null(flash_pass) && length(flash_pass) == nrow(df_b)) (flash_pass > 0) else FALSE
              if (isTRUE(flash_counted)) {
                keep_plot <- keep_plot & (!df_b$counted | is_flashing)
              } else {
                keep_plot <- keep_plot & !df_b$counted
              }
            }

            df_plot <- df_b[keep_plot, , drop = FALSE]

            if (nrow(df_plot) > 0) {
              plot_colors <- palette_colors[keep_plot]

              if (isTRUE(mask) && has_masks && exists("conts") && length(conts) > 0) {
                for (kp_idx in which(keep_plot)) {
                  mid <- df_b$orig_idx[kp_idx]
                  if (!is.na(mid) && mid <= length(conts) && !is.null(conts[[mid]]) && nrow(conts[[mid]]) >= 3) {
                    graphics::polygon(conts[[mid]][, 1], conts[[mid]][, 2],
                                      col = grDevices::adjustcolor(palette_colors[kp_idx], alpha.f = alpha),
                                      border = NA)
                  }
                }
              }

              if (isTRUE(bbox)) {
                .plot_yolo_bboxes(
                  boxes = df_plot,
                  palette_colors = plot_colors,
                  lwd = lwd,
                  cex = text_size_num,
                  pad = pad,
                  badge = badge,
                  show_text = isTRUE(show_text),
                  show_conf = isTRUE(show_conf),
                  show_class = isTRUE(show_class),
                  show_id = isTRUE(show_id)
                )
              }
              if (isTRUE(skeleton) && has_kpts) {
                kpts_for_plot <- lapply(which(keep_plot), function(k) {
                  vals <- raw_res$keypoints[k, ]
                  data.frame(
                    id = df_b$id[k],
                    keypoint = .coco_keypoints,
                    x = vals[seq(1, 51, by = 3)],
                    y = vals[seq(2, 51, by = 3)],
                    conf = round(vals[seq(3, 51, by = 3)], 4),
                    stringsAsFactors = FALSE
                  )
                })
                .plot_yolo_keypoints(
                  keypoints_list = kpts_for_plot,
                  palette_colors = plot_colors,
                  kpt_threshold = kpt_threshold,
                  lwd = lwd,
                  kpt_radius = kpt_radius
                )
              }
            }
          }
          if (isTRUE(show_fps)) {
            elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
            cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
            hud_items <- if (isTRUE(track) || !is.null(count_line)) {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est),
                sprintf("IN ZONE: %d", in_zone_val),
                sprintf("COUNTED: %d", total_count_val))
            } else if (detect) {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est),
                sprintf("OBJECTS: %d", nrow(df_b)))
            } else {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est))
            }
            u <- graphics::par("usr")
            th <- abs(graphics::strheight("Ag", units = "user", cex = 0.65))
            pad <- 6
            is_vert <- (hud_layout == "vertical")

            if (is_vert) {
              all_lines <- c("TELEMETRY", hud_items)
              max_w <- max(abs(graphics::strwidth(all_lines, units = "user", cex = 0.65, font = 2)))
              bx_w <- max_w + 2 * pad + 10
              bx_h <- (length(all_lines) * th * 1.5) + 2 * pad
            } else {
              single_line <- paste(hud_items, collapse = "  |  ")
              bx_w <- abs(graphics::strwidth(single_line, units = "user", cex = 0.65, font = 2)) + 2 * pad
              bx_h <- th * 1.8 + 2 * pad
            }

            margin_x <- abs(u[2] - u[1]) * 0.02
            margin_y <- abs(u[4] - u[3]) * 0.02
            top_y <- max(u[3], u[4]) - margin_y
            bot_y <- min(u[3], u[4]) + margin_y
            left_x <- min(u[1], u[2]) + margin_x
            right_x <- max(u[1], u[2]) - margin_x

            if (hud_pos %in% c("top-right", "topright")) {
              x1 <- right_x - bx_w; x2 <- right_x; y2 <- top_y; y1 <- top_y - bx_h
            } else if (hud_pos %in% c("bottom-right", "bottomright")) {
              x1 <- right_x - bx_w; x2 <- right_x; y1 <- bot_y; y2 <- bot_y + bx_h
            } else if (hud_pos %in% c("bottom-left", "bottomleft")) {
              x1 <- left_x; x2 <- left_x + bx_w; y1 <- bot_y; y2 <- bot_y + bx_h
            } else if (hud_pos == "top") {
              mid_x <- (u[1] + u[2]) / 2; x1 <- mid_x - bx_w / 2; x2 <- mid_x + bx_w / 2; y2 <- top_y; y1 <- top_y - bx_h
            } else if (hud_pos == "bottom") {
              mid_x <- (u[1] + u[2]) / 2; x1 <- mid_x - bx_w / 2; x2 <- mid_x + bx_w / 2; y1 <- bot_y; y2 <- bot_y + bx_h
            } else {
              x1 <- left_x; x2 <- left_x + bx_w; y2 <- top_y; y1 <- top_y - bx_h
            }

            graphics::rect(xleft = x1, ybottom = y1, xright = x2, ytop = y2,
                           col = grDevices::adjustcolor("#0E1420", alpha.f = 0.85), border = "#233C55", lwd = 1.2)
            graphics::rect(xleft = x1, ybottom = y2 - 2, xright = x2, ytop = y2,
                           col = "#00E5FF", border = NA)
            if (is_vert) {
              graphics::points(x1 + pad + 2, y2 - pad - th / 2, pch = 16, col = "#00FF88", cex = 0.8)
              graphics::text(x = x1 + pad + 10, y = y2 - pad - th / 2, labels = "TELEMETRY", col = "#00E5FF", cex = 0.6, font = 2, adj = c(0, 0.5))
              for (k in seq_along(hud_items)) {
                item_y <- y2 - pad - th * 1.5 - (k - 0.5) * th * 1.5
                parts <- strsplit(hud_items[k], ": ")[[1]]
                lbl <- parts[1]; val <- if (length(parts) > 1) parts[2] else ""
                val_col <- if (grepl("COUNT", lbl)) "#00E5FF" else if (grepl("ZONE", lbl)) "#FFBB00" else if (grepl("SPEED", lbl)) "#00FF88" else "#FFFFFF"
                graphics::text(x = x1 + pad, y = item_y, labels = lbl, col = "#94B0C4", cex = 0.6, font = 2, adj = c(0, 0.5))
                if (nzchar(val)) {
                  graphics::text(x = x2 - pad, y = item_y, labels = val, col = val_col, cex = 0.6, font = 2, adj = c(1, 0.5))
                }
              }
            } else {
              graphics::text(x = x1 + pad, y = (y1 + y2) / 2, labels = single_line, col = "#FFFFFF", cex = 0.65, font = 2, adj = c(0, 0.5))
            }
          }
        }

        if (!is.null(pb)) {
          elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
          cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
          cli::cli_progress_update(id = pb, extra = list(current_fps = sprintf("%.1f", cur_fps_est)))
        }
      }
    }

    process_loop()

    if (save_video) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Encoding annotated video to MP4 with FFmpeg...")
      }
      av::av_encode_video(
        input = frame_paths,
        output = out_video_file,
        framerate = src_fps,
        verbose = FALSE
      )
      if (isTRUE(verbose)) cli::cli_alert_success("Video encoding complete: {.path {out_video_file}}")
    }

    all_dets <- if (length(frame_records) > 0) do.call(rbind, frame_records) else data.frame()
    all_kpts <- if (length(frame_keypoints) > 0) do.call(rbind, frame_keypoints) else NULL
    total_time <- as.numeric(Sys.time() - t_start, units = "secs")
    avg_fps <- if (total_time > 0 && num_frames > 0) (num_frames / total_time) else src_fps
  }

  if (!return_data) {
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Session complete: {if (is_camera) frame_count else num_frames} frame(s) processed at {round(avg_fps, 1)} FPS average.")
      if (save_video && !is.null(output) && file.exists(output)) {
        cli::cli_alert_success("Video saved to: {.path {output}}")
      }
    }
    return(invisible(NULL))
  }

  # ============================================================================
  # COMPREHENSIVE DATA & ECOLOGICAL SUMMARY (RICHNESS, ABUNDANCE & DIVERSITY)
  # ============================================================================
  tot_frames <- if (is_camera) frame_count else num_frames
  if (is_camera && avg_fps > 0) src_fps <- avg_fps
  duration_sec <- if (is_camera) t_elapsed else (tot_frames / src_fps)
  v_source <- if (is_camera) "Live Webcam" else normalizePath(video, winslash = "/")

  if (is.null(all_dets) || nrow(all_dets) == 0) {
    all_dets <- data.frame(
      frame = integer(0), timestamp = numeric(0), id = integer(0),
      xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
      label = character(0), score = numeric(0), class_id = integer(0),
      stringsAsFactors = FALSE
    )
    all_classes <- character(0)
    summary_df <- data.frame(
      video_source = basename(v_source),
      model = basename(model_file),
      engine = toupper(engine),
      total_frames = tot_frames,
      fps = round(avg_fps, 2),
      duration_sec = round(duration_sec, 2),
      total_detections = 0L,
      total_richness = 0L,
      unique_objects = if (isTRUE(track)) 0L else NA_integer_,
      avg_abundance_frame = 0.0,
      max_abundance_frame = 0L,
      peak_frame = 1L,
      peak_time_sec = 0.0,
      shannon_index = 0.0,
      simpson_index = 0.0,
      pielou_evenness = 0.0,
      stringsAsFactors = FALSE
    )
    by_class_df <- data.frame(
      class = character(0),
      total_detections = integer(0),
      relative_abundance = numeric(0),
      unique_objects = integer(0),
      n_frames_present = integer(0),
      occurrence_rate = numeric(0),
      mean_abundance_when_present = numeric(0),
      mean_abundance_overall = numeric(0),
      max_per_frame = integer(0),
      mean_confidence = numeric(0),
      first_seen_sec = numeric(0),
      last_seen_sec = numeric(0),
      stringsAsFactors = FALSE
    )
    temporal_df <- data.frame(
      frame = seq_len(max(1, tot_frames)),
      timestamp = round((seq_len(max(1, tot_frames)) - 1) / src_fps, 3),
      abundance = 0L,
      richness = 0L,
      shannon = 0.0,
      stringsAsFactors = FALSE
    )
    cumulative_df <- data.frame(
      frame = seq_len(max(1, tot_frames)),
      timestamp = round((seq_len(max(1, tot_frames)) - 1) / src_fps, 3),
      cumulative_detections = 0L,
      cumulative_richness = 0L,
      stringsAsFactors = FALSE
    )
    if (isTRUE(track)) cumulative_df$cumulative_unique_objects <- 0L
    diversity_metrics <- list(
      richness_S = 0L,
      total_abundance_N = 0L,
      shannon_H = 0.0,
      simpson_1_minus_D = 0.0,
      simpson_D = 0.0,
      pielou_J = 0.0,
      berger_parker_d = 0.0,
      mean_richness_per_frame = 0.0,
      mean_abundance_per_frame = 0.0
    )
  } else {
    all_classes <- sort(unique(all_dets$label))
    total_det <- nrow(all_dets)
    total_rich <- length(all_classes)

    # Cross-tabulate frame x class
    mat_counts <- table(
      factor(all_dets$frame, levels = seq_len(tot_frames)),
      factor(all_dets$label, levels = all_classes)
    )
    abundance_vec <- as.integer(rowSums(mat_counts))
    richness_vec <- as.integer(rowSums(mat_counts > 0))

    # Shannon index per frame
    props <- mat_counts / pmax(1, abundance_vec)
    log_props <- ifelse(props > 0, log(props), 0)
    shannon_vec <- -rowSums(props * log_props)

    timestamps_vec <- round((seq_len(tot_frames) - 1) / src_fps, 3)

    temporal_df <- data.frame(
      frame = seq_len(tot_frames),
      timestamp = timestamps_vec,
      abundance = abundance_vec,
      richness = richness_vec,
      shannon = round(shannon_vec, 3),
      stringsAsFactors = FALSE
    )
    temporal_df <- cbind(temporal_df, as.data.frame.matrix(mat_counts))

    # Cumulative richness and detections (Collector's Curve)
    first_frames <- tapply(all_dets$frame, all_dets$label, min)
    first_frame_counts <- tabulate(first_frames, nbins = tot_frames)
    cum_richness <- cumsum(first_frame_counts)
    cum_detections <- cumsum(abundance_vec)

    cumulative_df <- data.frame(
      frame = seq_len(tot_frames),
      timestamp = timestamps_vec,
      cumulative_detections = cum_detections,
      cumulative_richness = cum_richness,
      stringsAsFactors = FALSE
    )
    if (isTRUE(track) && "id" %in% names(all_dets) && any(all_dets$id > 0)) {
      first_seen_id <- tapply(all_dets$frame, all_dets$id, min)
      first_seen_counts <- tabulate(first_seen_id, nbins = tot_frames)
      cumulative_df$cumulative_unique_objects <- cumsum(first_seen_counts)
    }

    # Summary by class
    counts_by_cls <- if (isTRUE(track) || !is.null(count_line)) track_state$counts_by_class else NULL
    by_class_list <- lapply(all_classes, function(cls) {
      df_c <- all_dets[all_dets$label == cls, , drop = FALSE]
      n_c <- nrow(df_c)
      n_present <- sum(mat_counts[, cls] > 0)
      cnt_val <- if (!is.null(counts_by_cls) && cls %in% names(counts_by_cls)) as.integer(counts_by_cls[[cls]]) else NA_integer_
      data.frame(
        class = cls,
        total_detections = n_c,
        relative_abundance = round(n_c / total_det * 100, 2),
        unique_objects = if (isTRUE(track) && any(df_c$id > 0)) length(unique(df_c$id[df_c$id > 0])) else NA_integer_,
        n_frames_present = n_present,
        occurrence_rate = round(n_present / tot_frames * 100, 2),
        mean_abundance_when_present = round(n_c / max(1, n_present), 2),
        mean_abundance_overall = round(n_c / tot_frames, 2),
        max_per_frame = max(mat_counts[, cls]),
        mean_confidence = round(mean(df_c$score), 4),
        first_seen_sec = round(min(df_c$timestamp), 3),
        last_seen_sec = round(max(df_c$timestamp), 3),
        total_counted = cnt_val,
        stringsAsFactors = FALSE
      )
    })
    by_class_df <- do.call(rbind, by_class_list)
    by_class_df <- by_class_df[order(-by_class_df$total_detections), ]
    rownames(by_class_df) <- NULL

    if (is.null(counts_by_cls)) {
      by_class_df$total_counted <- NULL
    }

    # Overall Diversity Metrics
    p_all <- table(all_dets$label) / total_det
    shannon_total <- -sum(p_all * log(p_all))
    simpson_d <- sum(p_all^2)
    simpson_total <- 1.0 - simpson_d
    pielou_total <- if (total_rich > 1) shannon_total / log(total_rich) else if (total_rich == 1) 1.0 else 0.0

    peak_idx <- which.max(abundance_vec)

    summary_df <- data.frame(
      video_source = basename(v_source),
      model = basename(model_file),
      engine = toupper(engine),
      total_frames = tot_frames,
      fps = round(avg_fps, 2),
      duration_sec = round(duration_sec, 2),
      total_detections = total_det,
      total_richness = total_rich,
      unique_objects = if (isTRUE(track) && any(all_dets$id > 0)) length(unique(all_dets$id[all_dets$id > 0])) else NA_integer_,
      avg_abundance_frame = round(total_det / max(1, tot_frames), 2),
      max_abundance_frame = max(abundance_vec),
      peak_frame = peak_idx,
      peak_time_sec = timestamps_vec[peak_idx],
      shannon_index = round(shannon_total, 3),
      simpson_index = round(simpson_total, 3),
      pielou_evenness = round(pielou_total, 3),
      stringsAsFactors = FALSE
    )
    if (!is.null(total_count_val)) {
      summary_df$total_counted <- total_count_val
    }

    diversity_metrics <- list(
      richness_S = total_rich,
      total_abundance_N = total_det,
      shannon_H = round(shannon_total, 3),
      simpson_1_minus_D = round(simpson_total, 3),
      simpson_D = round(simpson_d, 3),
      pielou_J = round(pielou_total, 3),
      berger_parker_d = round(max(p_all), 3),
      mean_richness_per_frame = round(mean(richness_vec), 2),
      mean_abundance_per_frame = round(mean(abundance_vec), 2)
    )
  }

  # Aggregate per-unique-object summary table (tracking / shape features)
  unique_objects_df <- data.frame()
  if ((isTRUE(track) || !is.null(count_line)) && !is.null(all_dets) && nrow(all_dets) > 0 && "id" %in% names(all_dets)) {
    valid_dets <- all_dets[!is.na(all_dets$id) & all_dets$id > 0L, , drop = FALSE]
    if (nrow(valid_dets) > 0) {
      split_u <- split(valid_dets, valid_dets$id)
      has_seg_feat <- "area" %in% names(valid_dets)
      uo_list <- lapply(split_u, function(sd) {
        first_f <- min(sd$frame)
        last_f <- max(sd$frame)
        n_f <- nrow(sd)
        dur <- round((last_f - first_f + 1) / src_fps, 2)
        score_m <- round(mean(sd$score, na.rm = TRUE), 3)
        col_val <- if ("color" %in% names(sd)) sd$color[1] else "#00CC66"
        lbl_tbl <- table(sd$label)
        lbl_val <- names(lbl_tbl)[which.max(lbl_tbl)]
        row_df <- data.frame(
          id = sd$id[1],
          label = lbl_val,
          n_frames = n_f,
          first_frame = first_f,
          last_frame = last_f,
          duration_sec = dur,
          mean_score = score_m,
          stringsAsFactors = FALSE
        )
        if (has_seg_feat) {
          row_df$area <- round(mean(sd$area, na.rm = TRUE), 1)
          row_df$area_sd <- if (n_f > 1) round(stats::sd(sd$area, na.rm = TRUE), 1) else 0.0
          row_df$perimeter <- round(mean(sd$perimeter, na.rm = TRUE), 1)
          row_df$radius_mean <- round(mean(sd$radius_mean, na.rm = TRUE), 2)
          row_df$length <- round(mean(sd$length, na.rm = TRUE), 1)
          row_df$width <- round(mean(sd$width, na.rm = TRUE), 1)
          row_df$circularity <- round(mean(sd$circularity, na.rm = TRUE), 4)
          row_df$eccentricity <- round(mean(sd$eccentricity, na.rm = TRUE), 4)
          row_df$asp_ratio <- round(mean(sd$asp_ratio, na.rm = TRUE), 3)
        } else {
          row_df$bbox_width <- round(mean(sd$xmax - sd$xmin, na.rm = TRUE), 1)
          row_df$bbox_height <- round(mean(sd$ymax - sd$ymin, na.rm = TRUE), 1)
          row_df$bbox_area <- round(mean((sd$xmax - sd$xmin) * (sd$ymax - sd$ymin), na.rm = TRUE), 1)
          row_df$asp_ratio <- round(mean((sd$xmax - sd$xmin) / pmax(0.001, (sd$ymax - sd$ymin)), na.rm = TRUE), 3)
        }
        row_df$color <- col_val
        row_df
      })
      unique_objects_df <- do.call(rbind, uo_list)
      rownames(unique_objects_df) <- NULL
    }
  }

  res <- list(
    summary = summary_df,
    by_class = by_class_df,
    temporal = temporal_df,
    cumulative = cumulative_df,
    diversity = diversity_metrics,
    detections = all_dets,
    unique_objects = unique_objects_df
  )

  if (isTRUE(track) || !is.null(count_line)) {
    if (!is.null(total_count_val)) res$total_counted <- total_count_val
    if (length(track_state$crossing_events) > 0) {
      ev_df <- do.call(rbind, lapply(track_state$crossing_events, as.data.frame))
      res$crossings <- ev_df
    }
  }
  if (!is.null(all_kpts) && nrow(all_kpts) > 0) res$keypoints <- all_kpts
  if (!is.null(roi)) {
    res$roi <- roi
    res$roi_mode <- roi_mode
  }
  if (!is.null(count_line)) res$count_line <- count_line
  if (!is.null(output)) res$output <- normalizePath(output, winslash = "/", mustWork = FALSE)

  class(res) <- c("pliman_video_detect", "list")

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Inference complete: {tot_frames} frames processed at {round(avg_fps, 1)} FPS average.")
    if (!is.null(output) && file.exists(output)) {
      cli::cli_alert_success("Annotated video saved to: {.path {output}}")
    }
  }

  return(res)
}

#' @export
print.pliman_video_detect <- function(x, ...) {
  s <- x$summary
  cli::cli_h2("YOLO Video Detection & Ecological Summary")
  cli::cli_bullets(c(
    "*" = sprintf("Source: {.file %s} | Model: {.val %s} [{toupper(s$engine)}]", s$video_source, s$model),
    "*" = sprintf("Duration: %.1f s (%d frames at %.1f FPS)", s$duration_sec, s$total_frames, s$fps),
    "*" = sprintf("Total Detections: %s (Avg: %.1f objs/frame | Peak: %d objs at frame %d [%.1fs])",
                  format(s$total_detections, big.mark = ","), s$avg_abundance_frame, s$max_abundance_frame, s$peak_frame, s$peak_time_sec),
    "*" = sprintf("Richness (S): %d distinct classes observed", s$total_richness),
    "*" = sprintf("Diversity Indices: Shannon (H') = %.3f | Simpson (1-D) = %.3f | Pielou (J') = %.3f",
                  s$shannon_index, s$simpson_index, s$pielou_evenness)
  ))

  if (!is.na(s$unique_objects)) {
    cli::cli_bullets(c(
      "*" = sprintf("Persistent Tracked Objects: %d unique IDs", s$unique_objects)
    ))
  }
  if ("total_counted" %in% names(s)) {
    cli::cli_bullets(c(
      "*" = sprintf("Tripwire Counting (Line): %d objects crossed", s$total_counted)
    ))
  }

  bc <- x$by_class
  if (!is.null(bc) && nrow(bc) > 0) {
    cli::cli_h3("Class Breakdown (Abundance & Occurrence)")
    show_cols <- c("class", "total_detections", "relative_abundance", "occurrence_rate", "max_per_frame", "mean_confidence")
    if (!is.na(s$unique_objects) && "unique_objects" %in% names(bc)) {
      show_cols <- c("class", "total_detections", "relative_abundance", "unique_objects", "occurrence_rate", "max_per_frame", "mean_confidence")
    }
    if ("total_counted" %in% names(bc)) {
      show_cols <- c(show_cols, "total_counted")
    }
    print(bc[, show_cols, drop = FALSE], row.names = FALSE)
  }

  uo <- x$unique_objects
  if (!is.null(uo) && nrow(uo) > 0) {
    cli::cli_h3("Unique Tracked Objects Summary")
    if (nrow(uo) <= 10) {
      print(uo, row.names = FALSE)
    } else {
      print(head(uo, 8), row.names = FALSE)
      cli::cli_alert_info("... and {nrow(uo) - 8} more objects (access full table with $unique_objects)")
    }
  }

  if (!is.null(x$output)) {
    cli::cli_alert_info("Annotated Video: {.file {x$output}}")
  }

  cli::cli_text("{.muted List components: $summary, $by_class, $temporal, $cumulative, $diversity, $detections, $unique_objects}")
  invisible(x)
}

#' @export
summary.pliman_video_detect <- function(object, ...) {
  print(object, ...)
  invisible(list(summary = object$summary, by_class = object$by_class, diversity = object$diversity, unique_objects = object$unique_objects))
}

#' @export
as.data.frame.pliman_video_detect <- function(x, ...) {
  x$detections
}

#' @importFrom utils head tail
#' @export
head.pliman_video_detect <- function(x, n = 6L, ...) {
  head(x$detections, n, ...)
}

#' @export
tail.pliman_video_detect <- function(x, n = 6L, ...) {
  tail(x$detections, n, ...)
}

#' @export
dim.pliman_video_detect <- function(x) {
  dim(x$detections)
}

nrow.pliman_video_detect <- function(x) {
  nrow(x$detections)
}

ncol.pliman_video_detect <- function(x) {
  ncol(x$detections)
}

#' @export
`$.pliman_video_detect` <- function(x, name) {
  if (name %in% names(x)) {
    return(.subset2(x, name))
  }
  if ("detections" %in% names(x) && name %in% names(x$detections)) {
    return(.subset2(x$detections, name))
  }
  NULL
}

#' @export
`[[.pliman_video_detect` <- function(x, i, exact = TRUE) {
  if (is.character(i) && !(i %in% names(x)) && "detections" %in% names(x) && i %in% names(x$detections)) {
    return(x$detections[[i, exact = exact]])
  }
  NextMethod("[[")
}

#' @export
`[.pliman_video_detect` <- function(x, i, j, ...) {
  if (missing(j)) {
    if (is.character(i) && all(i %in% names(x))) {
      return(NextMethod("["))
    }
    return(x$detections[i, , drop = FALSE])
  }
  x$detections[i, j, ...]
}

#' @export
plot.pliman_video_detect <- function(x, type = c("temporal", "richness", "abundance", "class", "cumulative", "all"), ...) {
  type <- match.arg(type)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  temp <- x$temporal
  cum <- x$cumulative
  bcls <- x$by_class
  if (is.null(temp) || nrow(temp) == 0 || is.null(bcls) || nrow(bcls) == 0) {
    message("No detections to plot.")
    return(invisible(NULL))
  }

  classes <- bcls$class
  n_cls <- length(classes)
  colors <- .get_class_palette(classes, n_cls, rainbow = TRUE)
  names(colors) <- classes

  if (type == "temporal") {
    graphics::par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    # Panel 1: Abundance over time
    plot(temp$timestamp, temp$abundance, type = "n",
         xlab = "", ylab = "Abundance (N)",
         main = sprintf("Temporal Dynamics of Detections (Abundance & Richness) - %s", x$summary$model),
         bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$abundance, 0),
                      col = grDevices::adjustcolor("#00CC66", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$abundance, col = "#00994C", lwd = 2.5)

    # Panel 2: Richness over time
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "n",
         xlab = "Time (seconds)", ylab = "Richness (S)",
         bty = "l", las = 1, ylim = c(0, max_rich + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$richness, 0),
                      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$richness, col = "#0072B2", lwd = 2.5, type = "s")

  } else if (type == "abundance") {
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    max_c <- if (length(classes) > 0) max(temp[, classes, drop = FALSE], na.rm = TRUE) else 1
    plot(NULL, xlim = range(temp$timestamp), ylim = c(0, max_c + 1),
         xlab = "Time (seconds)", ylab = "Abundance per Class",
         main = "Class Abundance Over Time",
         bty = "l", las = 1, ...)
    graphics::grid(col = "#E0E0E0", lty = 2)
    for (cls in classes) {
      graphics::lines(temp$timestamp, temp[[cls]], col = colors[cls], lwd = 2)
    }
    graphics::legend("topright", legend = classes, col = colors, lwd = 2, bty = "n", cex = 0.85)

  } else if (type == "richness") {
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "n",
         xlab = "Time (seconds)", ylab = "Class Richness (S)",
         main = "Class Richness Over Time",
         bty = "l", las = 1, ylim = c(0, max_rich + 1), ...)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$richness, 0),
                      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$richness, col = "#0072B2", lwd = 2.5, type = "s")

  } else if (type == "cumulative") {
    graphics::par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    # Cumulative richness (species accumulation curve / collector's curve)
    max_cum_r <- max(cum$cumulative_richness, na.rm = TRUE)
    plot(cum$timestamp, cum$cumulative_richness, type = "n",
         xlab = "", ylab = "Cumulative Richness (S)",
         main = "Species Accumulation Curve (Collector's Curve)",
         bty = "l", las = 1, ylim = c(0, max_cum_r + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::lines(cum$timestamp, cum$cumulative_richness, col = "#D55E00", lwd = 2.5, type = "s")

    # Cumulative detections
    plot(cum$timestamp, cum$cumulative_detections, type = "n",
         xlab = "Time (seconds)", ylab = "Cumulative Detections",
         main = "Cumulative Detections Over Time",
         bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::lines(cum$timestamp, cum$cumulative_detections, col = "#0072B2", lwd = 2.5)

  } else if (type == "class") {
    graphics::par(mar = c(4.5, 7, 2.5, 2), mgp = c(2.5, 0.7, 0))
    ord_bcls <- bcls[order(bcls$total_detections), ]
    bp <- graphics::barplot(ord_bcls$total_detections, horiz = TRUE, names.arg = ord_bcls$class,
                            col = colors[ord_bcls$class], border = NA, las = 1,
                            xlab = "Total Detections", main = "Abundance by Class",
                            xlim = c(0, max(ord_bcls$total_detections) * 1.25), ...)
    graphics::grid(nx = NULL, ny = NA, col = "#E0E0E0", lty = 2)
    lbl_txt <- sprintf("%d (%.1f%%)", ord_bcls$total_detections, ord_bcls$relative_abundance)
    graphics::text(ord_bcls$total_detections, bp, labels = paste0("  ", lbl_txt), pos = 4, cex = 0.85, font = 2)

  } else if (type == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(3.5, 4, 2.5, 1), mgp = c(2, 0.6, 0))
    # 1. Abundance over time
    plot(temp$timestamp, temp$abundance, type = "l", col = "#00994C", lwd = 2,
         xlab = "Time (s)", ylab = "Abundance (N)", main = "Total Abundance", bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)

    # 2. Richness over time
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "s", col = "#0072B2", lwd = 2,
         xlab = "Time (s)", ylab = "Richness (S)", main = "Class Richness", bty = "l", las = 1,
         ylim = c(0, max_rich + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)

    # 3. Barplot of classes
    ord_bcls <- bcls[order(bcls$total_detections), ]
    bp <- graphics::barplot(ord_bcls$total_detections, horiz = TRUE, names.arg = ord_bcls$class,
                            col = colors[ord_bcls$class], border = NA, las = 1,
                            xlab = "Detections", main = "Abundance by Class", cex.names = 0.75)
    graphics::grid(nx = NULL, ny = NA, col = "#E0E0E0", lty = 2)

    # 4. Cumulative richness
    max_cum_r <- max(cum$cumulative_richness, na.rm = TRUE)
    plot(cum$timestamp, cum$cumulative_richness, type = "s", col = "#D55E00", lwd = 2,
         xlab = "Time (s)", ylab = "Cumulative Richness", main = "Collector's Curve", bty = "l", las = 1,
         ylim = c(0, max_cum_r + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
  }
  invisible(x)
}

#' @rdname video_detect_dl
#' @export
video_detect <- video_detect_dl

#' Universal Human Pose Estimation with YOLO26
#'
#' Performs real-time human pose estimation and 17 keypoint detection using YOLO26
#' (e.g., `yolo26n-pose`, `yolo26s-pose`, etc.) with bounding box regression and
#' anatomical skeleton visualization.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n-pose"`.
#'   Options include `"yolo26n-pose"`, `"yolo26s-pose"`, `"yolo26m-pose"`, `"yolo26l-pose"`, `"yolo26x-pose"`.
#' @param conf_threshold Minimum confidence score for candidate person detections (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param kpt_threshold Minimum confidence threshold for rendering individual keypoints and skeleton limbs (default 0.3).
#' @param col Optional color or palette for bounding boxes and keypoint markers.
#' @param lwd Line width for bounding boxes and skeleton limbs (default 2).
#' @param kpt_radius Radius/size for keypoint markers (default 4).
#' @param bbox Logical. Whether to draw bounding boxes around detected persons (default `TRUE`).
#' @param skeleton Logical. Whether to render the 17-keypoint anatomical skeleton (default `TRUE`).
#' @param cex Size factor for bounding box label text (default `1.0`).
#' @param pad Padding multiplier controlling the background badge size around label text (default `1.0`).
#' @param show_text Logical. Whether to display text labels on bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in labels (default `TRUE`).
#' @param show_class Logical. Whether to display class names in labels (default `TRUE`).
#' @param show_id Logical. Whether to display object IDs (`#1`, `#2`, ...) in labels (default `FALSE`).
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the image with pose overlays.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list with class `pliman_pose` containing:
#'   * `boxes`: Data frame containing detected person bounding boxes (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`).
#'   * `keypoints`: List of data frames (one per person) containing the 17 COCO keypoints (`keypoint`, `x`, `y`, `conf`).
#'   * `summary`: Human-readable detection summary string.
#'   * `counts`: Data frame of instance counts.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("person.png")
#'   res <- image_pose_dl(img)
#'   res$boxes
#'   res$keypoints[[1]]
#' }
image_pose_dl <- function(img,
                          model = "yolo26n-pose",
                          conf_threshold = 0.25,
                          iou_threshold = 0.45,
                          kpt_threshold = 0.3,
                          col = NULL,
                          lwd = 2,
                          cex = 1.0,
                          pad = 1.0,
                          show_text = TRUE,
                          show_conf = TRUE,
                          show_class = TRUE,
                          show_id = FALSE,
                          badge = TRUE,
                          kpt_radius = 4,
                          bbox = TRUE,
                          skeleton = TRUE,
                          threads = 0,
                          engine = c("gpu", "cpu"),
                          device_id = -1,
                          verbose = TRUE,
                          plot = TRUE,
                          dir = pliman_model_dir(),
                          ...) {
  dots <- list(...)
  if ("label_size" %in% names(dots)) cex <- dots$label_size
  if ("text_size" %in% names(dots)) cex <- dots$text_size
  if ("font_scale" %in% names(dots)) cex <- dots$font_scale
  if ("cex_scale" %in% names(dots)) cex <- dots$cex_scale
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  if (length(cex) > 1L) {
    if (missing(pad)) pad <- cex[2]
    cex <- cex[1]
  }

  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_pose_dl(x, model = model, conf_threshold = conf_threshold, iou_threshold = iou_threshold,
                    kpt_threshold = kpt_threshold, col = col, lwd = lwd, cex = cex, pad = pad,
                    show_text = show_text, show_conf = show_conf, show_class = show_class,
                    show_id = show_id, badge = badge, kpt_radius = kpt_radius,
                    bbox = bbox, skeleton = skeleton, threads = threads, engine = engine,
                    device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "yolo26n-pose"

  if (isTRUE(verbose)) {
    pb <- cli::cli_process_start(
      "Running YOLO pose estimation [{toupper(engine)}]...",
      on_exit = "failed"
    )
  }

  res <- .run_yolo(
    mat = mat,
    model = model_str,
    conf_threshold = conf_threshold,
    iou_threshold = iou_threshold,
    labels = "person",
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  if (isTRUE(verbose)) {
    sum_txt <- if (!is.null(res$summary) && nzchar(res$summary) && res$summary != "0 objects") {
      paste0("pliman detected {.bold ", res$summary, "} in the image.")
    } else {
      "No objects detected in the image."
    }
    cli::cli_process_done(id = pb, msg_done = sum_txt)
  }

  df_boxes <- res$boxes
  num_inst <- nrow(df_boxes)

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    if (num_inst > 0) {
      palette_colors <- if (!is.null(col)) {
        rep(col, length.out = num_inst)
      } else if (num_inst == 1) {
        "salmon"
      } else {
        grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
      }

      if (isTRUE(bbox)) {
        .plot_yolo_bboxes(
          boxes = df_boxes,
          palette_colors = palette_colors,
          lwd = lwd,
          cex = cex,
          pad = pad,
          show_text = show_text,
          show_conf = show_conf,
          show_class = show_class,
          show_id = show_id,
          badge = badge
        )
      }

      if (isTRUE(skeleton) && length(res$keypoints) > 0) {
        .plot_yolo_keypoints(
          keypoints_list = res$keypoints,
          palette_colors = palette_colors,
          kpt_threshold = kpt_threshold,
          lwd = lwd,
          kpt_radius = kpt_radius
        )
      }
    }
  }

  out <- list(
    boxes = df_boxes,
    keypoints = res$keypoints,
    summary = res$summary,
    counts = res$counts
  )
  class(out) <- c("pliman_pose", "list")
  invisible(out)
}

#' Real-Time and Zero-Shot Image Classification
#'
#' Classifies images using pre-trained YOLO26 classification models across
#' 1,000 ImageNet categories, or open-vocabulary text candidates using OpenAI
#' CLIP (ViT-B/32) cross-modal embeddings and cosine similarity softmax.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param candidates An optional character vector of candidate class names or natural language descriptions.
#'   When supplied, performs zero-shot open-vocabulary classification using CLIP. Defaults to `NULL`.
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n-cls"` when `candidates = NULL`,
#'   or `"clip-vit-b32"` when `candidates` is provided.
#' @param top_k Integer. Number of top class predictions to return (defaults to 5 for YOLO models,
#'   or `length(candidates)` for CLIP zero-shot classification).
#' @param labels Optional character vector of class names (defaults to 1,000 standard ImageNet classes).
#' @param format_prompt Logical. If `TRUE` (default), candidate strings that do not contain
#'   `"photo"` or `"a "` are automatically formatted as `"a photo of a {candidate}"` for optimal CLIP zero-shot accuracy.
#' @param plot Logical. Whether to display the image with top predicted class labels (or horizontal bar chart for CLIP).
#' @param col Fill color for probability bar chart when `candidates` is supplied (default `"#2b5c8f"`).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()] or [graphics::barplot()].
#' @return A data frame containing top predictions.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("sample.png")
#'   # YOLO classification
#'   image_classify_dl(img)
#'
#'   # CLIP zero-shot classification
#'   res <- image_classify_dl(img, candidates = c("healthy leaf", "rust lesion", "insect damage"))
#'   print(res)
#' }
image_classify_dl <- function(img,
                              candidates = NULL,
                              model = NULL,
                              top_k = NULL,
                              labels = NULL,
                              format_prompt = TRUE,
                              plot = FALSE,
                              col = "#2b5c8f",
                              threads = 0,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              verbose = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_classify_dl(x, candidates = candidates, model = model, top_k = top_k,
                        labels = labels, format_prompt = format_prompt,
                        plot = FALSE, col = col, threads = threads,
                        engine = engine, device_id = device_id,
                        verbose = verbose, dir = dir, ...)
    })
    return(res)
  }

  if (!is.null(candidates)) {
    if (is.null(model)) model <- "clip-vit-b32"
    if (is.null(top_k)) top_k <- length(candidates)

    if (!is.character(candidates) || length(candidates) == 0L) {
      stop("candidates must be a non-empty character vector.")
    }

    # Prompt formatting for zero-shot accuracy
    prompts <- if (isTRUE(format_prompt)) {
      vapply(candidates, function(c) {
        cl <- tolower(trimws(c))
        if (grepl("^a |^an |^photo |^image ", cl)) c else paste0("a photo of a ", cl)
      }, character(1), USE.NAMES = FALSE)
    } else {
      candidates
    }

    # 1. Image embedding
    img_emb <- image_embed_dl(
      img = img,
      model = model,
      threads = threads,
      engine = engine,
      device_id = device_id,
      verbose = verbose,
      dir = dir,
      ...
    )

    # 2. Text embeddings
    txt_embs <- text_embed_dl(
      text = prompts,
      model = model,
      threads = threads,
      engine = engine,
      device_id = device_id,
      verbose = verbose,
      dir = dir,
      ...
    )

    # 3. Scaled cosine similarity logits (OpenAI CLIP temperature scale: 100)
    logits <- as.numeric(txt_embs %*% img_emb) * 100.0

    # 4. Softmax
    max_l <- max(logits)
    exp_l <- exp(logits - max_l)
    probs <- exp_l / sum(exp_l)

    df <- data.frame(
      candidate = candidates,
      prompt = prompts,
      probability = probs,
      logit = logits,
      stringsAsFactors = FALSE
    )
    df <- df[order(df$probability, decreasing = TRUE), ]
    df$rank <- seq_len(nrow(df))
    rownames(df) <- NULL

    if (top_k < nrow(df)) {
      df <- df[1:top_k, ]
    }

    if (isTRUE(plot)) {
      plot_df <- df[order(df$probability), ]
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(4.5, max(8, max(nchar(plot_df$candidate)) * 0.6), 3, 2))
      bp <- graphics::barplot(
        plot_df$probability * 100,
        names.arg = plot_df$candidate,
        horiz = TRUE,
        las = 1,
        col = col,
        border = NA,
        xlab = "Confidence Probability (%)",
        main = "CLIP Zero-Shot Classification",
        xlim = c(0, max(100, max(plot_df$probability * 100) * 1.15))
      )
      graphics::text(
        x = plot_df$probability * 100,
        y = bp,
        labels = sprintf(" %.1f%%", plot_df$probability * 100),
        pos = 4,
        cex = 0.85,
        col = "gray20"
      )
    }

    return(df)

  } else {
    if (is.null(model)) model <- "yolo26n-cls"
    if (is.null(top_k)) top_k <- 5L

    mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
    if (is.character(img) && file.exists(img[1])) {
      img_obj <- image_import(img[1])
      mat <- image_data(img_obj)
    }
    dims <- dim(mat)
    if (is.null(dims) || length(dims) < 2) {
      stop("img must be a 2D matrix, 3D array, or an 'image' object.")
    }

    is_url <- is.character(model) && grepl("^https?://", model[1])
    model_str <- if (is_url) {
      model[1]
    } else if (is.character(model)) {
      .resolve_model_name(model[1])
    } else {
      "yolo26n-cls"
    }

    if (file.exists(model_str)) {
      model_file <- normalizePath(model_str, winslash = "/")
    } else {
      model_file <- pliman_download_model(model = model_str, dir = dir)
    }

    # Dynamically detect model input size (default 224 for YOLO-Cls, 640 for detection)
    target_size <- 224L
    tryCatch({
      m_info <- inspect_onnx_model_cpp(model_file, lib_path)
      if (!is.null(m_info$shapes) && length(m_info$shapes) > 0L) {
        s <- m_info$shapes[[1]]
        if (length(s) >= 4 && s[3] > 0) {
          target_size <- as.integer(s[3])
        }
      }
    }, error = function(e) NULL)

    tensor <- .preprocess_yolo(mat, target_size = target_size)

    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Running YOLO image classification [{toupper(engine)}]...",
        msg_done = "Image classification complete"
      )
    }

    probs <- run_yolo_cls_cpp(
      tensor_vec = tensor,
      model_path = model_file,
      lib_path = lib_path,
      num_threads = as.integer(threads),
      use_gpu = use_gpu,
      device_id = as.integer(device_id)
    )

    if (!is.null(labels)) {
      class_names <- labels
    } else {
      model_classes <- .get_onnx_classes(model_file)
      if (!is.null(model_classes) && length(model_classes) > 0L) {
        class_names <- model_classes
      } else {
        class_names <- .imagenet_classes
      }
    }
    num_classes <- length(probs)
    k <- min(as.integer(top_k), num_classes)

    top_idx <- order(probs, decreasing = TRUE)[seq_len(k)]
    top_probs <- round(probs[top_idx], 4)
    top_labels <- ifelse(top_idx <= length(class_names), class_names[top_idx], paste0("class_", top_idx - 1L))

    df_res <- data.frame(
      rank = seq_len(k),
      class = top_labels,
      probability = top_probs,
      class_id = top_idx - 1L,
      stringsAsFactors = FALSE
    )

    if (isTRUE(plot)) {
      plot(as_image(mat), ...)
      top1 <- df_res[1, ]
      graphics::title(sub = sprintf("Top 1: %s (%.1f%%)", top1$class, top1$probability * 100), col.sub = "darkgreen", font.sub = 2)
    }

    if (isTRUE(verbose)) {
      cli::cli_alert_info("Top-{k} Classification Results:")
      for (i in seq_len(nrow(df_res))) {
        cli::cli_bullets(c("*" = sprintf("%d. {.val %s}: %.2f%%", df_res$rank[i], df_res$class[i], df_res$probability[i] * 100)))
      }
    }

    invisible(df_res)
  }
}

#' Star-Convex Object Detection with StarDist
#'
#' Detects round, convex, and overlapping objects (e.g. seeds, cells, nuclei, spores)
#' by predicting radial star-convex polygons and object probability using StarDist.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"stardist"`.
#' @param prob_threshold Probability threshold for detecting object centers (default 0.5).
#' @param nms_threshold Polygon IoU threshold for Non-Maximum Suppression (default 0.3).
#' @param type Output type: `"segment"` (highlight overlay), `"mask"` (integer labels), or `"polygons"`.
#' @param col_highlight Fill color for overlay polygons (default `"salmon"`).
#' @param border Border color for polygon contours (default `"white"`).
#' @param lwd Line width for polygon contours (default 2).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the results.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A list containing:
#'   * `boxes`: Bounding boxes data frame with centers and scores.
#'   * `polygons_x`, `polygons_y`: Coordinates of radial polygon vertices for each object.
#'   * `labels`: Integer matrix with instance IDs.
#'   * `mask`: Logical foreground mask.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("grains.png")
#'   res <- image_stardist_dl(img)
#' }
image_stardist_dl <- function(img,
                              model = "stardist",
                              prob_threshold = 0.5,
                              nms_threshold = 0.3,
                              type = c("segment", "mask", "polygons"),
                              col_highlight = "salmon",
                              border = "white",
                              lwd = 2,
                              threads = 0,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  type <- match.arg(type)
  engine <- match.arg(engine)

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_stardist_dl(x, model = model, prob_threshold = prob_threshold, nms_threshold = nms_threshold,
                        type = type, col_highlight = col_highlight, border = border, lwd = lwd,
                        threads = threads, engine = engine, device_id = device_id, verbose = verbose,
                        plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "stardist"

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running StarDist polygon detection [{toupper(engine)}]...",
      msg_done = "StarDist detection complete"
    )
  }

  res <- .run_stardist(
    mat = mat,
    model = model_str,
    prob_threshold = prob_threshold,
    nms_threshold = nms_threshold,
    threads = threads,
    engine = engine,
    device_id = device_id,
    dir = dir
  )

  num_inst <- nrow(res$boxes)

  if (isTRUE(plot)) {
    if (type == "mask") {
      plot(as_image(res$labels, colormode = "Grayscale"), ...)
    } else {
      plot(as_image(mat), ...)
      if (num_inst > 0) {
        palette_colors <- if (num_inst == 1) col_highlight else grDevices::rainbow(num_inst, s = 0.85, v = 0.95)
        for (i in seq_len(num_inst)) {
          px <- res$polygons_x[[i]]
          py <- res$polygons_y[[i]]
          k_col <- palette_colors[((i - 1) %% length(palette_colors)) + 1]
          col_poly <- grDevices::adjustcolor(k_col, alpha.f = 0.4)
          graphics::polygon(x = c(px, px[1]), y = c(py, py[1]), col = col_poly, border = border, lwd = lwd)
        }
      }
    }
  }

  if (isTRUE(verbose)) {
    .print_detection_summary(res$summary, title = "StarDist Object Detection Summary")
  }

  invisible(res)
}

#' 4x Super-Resolution with Real-ESRGAN Compact
#'
#' Upscales images by 4x using generative convolutional super-resolution with
#' intelligent tiled overlapping to prevent GPU/CPU memory exhaustion.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX file. Defaults to `"realesrgan-compact"`.
#' @param scale Upscaling factor (default 4).
#' @param tile_size Processing tile size in pixels (default 256).
#' @param tile_pad Overlap margin in pixels for seamless blending (default 16).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param plot Logical. Whether to plot the super-resolved image.
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return An `Image` object with dimensions `(orig_w * 4) x (orig_h * 4) x 3`.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("low_res.png")
#'   img_hr <- image_superres_dl(img)
#' }
image_superres_dl <- function(img,
                              model = "realesrgan-compact",
                              scale = 4,
                              tile_size = 256,
                              tile_pad = 16,
                              threads = 0,
                              engine = c("gpu", "cpu"),
                              device_id = -1,
                              verbose = TRUE,
                              plot = TRUE,
                              dir = pliman_model_dir(),
                              ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_superres_dl(x, model = model, scale = scale, tile_size = tile_size,
                        tile_pad = tile_pad, threads = threads, engine = engine,
                        device_id = device_id, verbose = verbose, plot = FALSE, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]
  nch <- if (length(dims) >= 3) dims[3] else 1

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "realesrgan-compact"
  if (file.exists(model[1])) {
    model_file <- normalizePath(model[1], winslash = "/")
  } else {
    model_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running Real-ESRGAN {scale}x super-resolution [{toupper(engine)}]...",
      msg_done = "Super-resolution complete"
    )
  }

  if (is.raw(mat)) {
    val_scale <- 1 / 255.0
  } else {
    max_val <- max(mat[1:min(1000, length(mat))], na.rm = TRUE)
    val_scale <- if (max_val > 1.5) (1 / 255.0) else 1.0
  }

  if (nch >= 3) {
    R <- as.numeric(mat[, , 1]) * val_scale
    G <- as.numeric(mat[, , 2]) * val_scale
    B <- as.numeric(mat[, , 3]) * val_scale
  } else {
    R <- G <- B <- as.numeric(mat) * val_scale
  }
  tensor <- c(R, G, B)

  raw_vec <- run_super_resolution_cpp(
    tensor_vec = tensor,
    in_w = orig_w,
    in_h = orig_h,
    scale = as.integer(scale),
    tile_size = as.integer(tile_size),
    tile_pad = as.integer(tile_pad),
    model_path = model_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  out_w <- orig_w * as.integer(scale)
  out_h <- orig_h * as.integer(scale)
  plane <- out_w * out_h

  out_arr <- array(0.0, dim = c(out_w, out_h, 3L))
  out_arr[, , 1] <- matrix(raw_vec[1:plane], nrow = out_w, ncol = out_h)
  out_arr[, , 2] <- matrix(raw_vec[(plane + 1):(2 * plane)], nrow = out_w, ncol = out_h)
  out_arr[, , 3] <- matrix(raw_vec[(2 * plane + 1):(3 * plane)], nrow = out_w, ncol = out_h)

  out_img <- as_image(out_arr)

  if (isTRUE(plot)) {
    plot(out_img, ...)
  }

  invisible(out_img)
}

# ==============================================================================
# SECTION 8: OPENAI CLIP (ViT-B/32) ZERO-SHOT EMBEDDINGS & CLASSIFICATION
# ==============================================================================

# Native BPE Tokenizer for CLIP (Zero external Python dependencies)
.clip_bpe <- function(token, bpe_ranks) {
  chars <- strsplit(token, "")[[1]]
  if (length(chars) == 0L) return(character(0))
  chars[length(chars)] <- paste0(chars[length(chars)], "</w>")
  word <- chars

  get_pairs <- function(w) {
    if (length(w) < 2L) return(character(0))
    paste(w[-length(w)], w[-1], sep = " ")
  }

  pairs <- get_pairs(word)
  if (length(pairs) == 0L) return(word)

  repeat {
    ranks <- bpe_ranks[pairs]
    ranks <- ranks[!is.na(ranks)]
    if (length(ranks) == 0L) break

    min_rank <- min(ranks)
    best_pair <- names(ranks)[ranks == min_rank][1]
    pair_parts <- strsplit(best_pair, " ")[[1]]

    new_word <- character(0)
    i <- 1L
    while (i <= length(word)) {
      if (i < length(word) && word[i] == pair_parts[1] && word[i + 1L] == pair_parts[2]) {
        new_word <- c(new_word, paste0(pair_parts[1], pair_parts[2]))
        i <- i + 2L
      } else {
        new_word <- c(new_word, word[i])
        i <- i + 1L
      }
    }
    word <- new_word
    if (length(word) <= 1L) break
    pairs <- get_pairs(word)
  }
  return(word)
}

.clip_tokenize <- function(texts, vocab_path, merges_path, max_len = 77L) {
  if (!file.exists(vocab_path)) stop("CLIP vocab file not found: ", vocab_path)
  if (!file.exists(merges_path)) stop("CLIP merges file not found: ", merges_path)

  if (!exists(vocab_path, envir = .pliman_vocab_cache, inherits = FALSE)) {
    vocab <- jsonlite::fromJSON(vocab_path)
    assign(vocab_path, vocab, envir = .pliman_vocab_cache)
  } else {
    vocab <- get(vocab_path, envir = .pliman_vocab_cache)
  }

  if (!exists(merges_path, envir = .pliman_vocab_cache, inherits = FALSE)) {
    merges_lines <- readLines(merges_path, warn = FALSE)
    merges_lines <- merges_lines[!grepl("^#", merges_lines)]
    bpe_ranks <- setNames(seq_along(merges_lines), merges_lines)
    assign(merges_path, bpe_ranks, envir = .pliman_vocab_cache)
  } else {
    bpe_ranks <- get(merges_path, envir = .pliman_vocab_cache)
  }

  batch_size <- length(texts)
  ids_mat <- matrix(0L, nrow = batch_size, ncol = max_len)
  mask_mat <- matrix(0L, nrow = batch_size, ncol = max_len)

  for (b in seq_len(batch_size)) {
    text <- tolower(trimws(texts[b]))
    words <- unlist(regmatches(text, gregexpr("<\\|startoftext\\|>|<\\|endoftext\\|>|'s|'t|'re|'ve|'m|'ll|'d|[[:alpha:]]+|[[:digit:]]|[^\\s[:alnum:]]+", text, perl = TRUE)))

    tokens <- 49406L # <|startoftext|>
    for (w in words) {
      subwords <- .clip_bpe(w, bpe_ranks)
      for (sw in subwords) {
        id <- vocab[[sw]]
        if (!is.null(id)) tokens <- c(tokens, as.integer(id))
      }
    }
    tokens <- c(tokens, 49407L) # <|endoftext|>
    if (length(tokens) > max_len) {
      tokens <- tokens[1:max_len]
      tokens[max_len] <- 49407L
    }
    n_tok <- length(tokens)
    ids_mat[b, 1:n_tok] <- tokens
    mask_mat[b, 1:n_tok] <- 1L
  }

  list(input_ids = ids_mat, attention_mask = mask_mat)
}

#' Extract Image Embeddings with OpenAI CLIP ViT-B/32
#'
#' Computes L2-normalized 512-dimensional semantic embeddings for one or more images
#' using OpenAI's CLIP (ViT-B/32) vision encoder directly via C++ ONNX Runtime.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param model Model name or path to a custom ONNX vision model file (default `"clip-vit-b32"`).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param dir Directory where models are stored (default `pliman_model_dir()`).
#' @param ... Additional arguments.
#' @return A numeric vector (if single image) or numeric matrix of dimensions `N x 512`
#'   containing the L2-normalized image embeddings.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("leaf.png")
#'   emb <- image_embed_dl(img)
#'   length(emb) # 512
#' }
image_embed_dl <- function(img,
                           model = "clip-vit-b32",
                           threads = 0,
                           engine = c("gpu", "cpu"),
                           device_id = -1,
                           verbose = FALSE,
                           dir = pliman_model_dir(),
                           ...) {
  engine <- match.arg(engine)
  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_embed_dl(x, model = model, threads = threads, engine = engine,
                     device_id = device_id, verbose = verbose, dir = dir, ...)
    })
    return(do.call(rbind, res))
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  if (is.character(img) && file.exists(img[1])) {
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  }
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "clip-vit-b32"
  if (file.exists(model[1])) {
    vis_file <- normalizePath(model[1], winslash = "/")
  } else {
    vis_file <- pliman_download_model(model = model_str, dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running CLIP vision encoder [{toupper(engine)}]...",
      msg_done = "CLIP image feature extraction complete"
    )
  }

  in_size <- 224L
  tensor <- .preprocess_nchw(
    mat,
    target_size = in_size,
    mean = c(0.48145466, 0.4578275, 0.40821073),
    std = c(0.26862954, 0.26130258, 0.27577711),
    letterbox = FALSE
  )

  emb <- run_clip_vision_cpp(
    tensor_vec = tensor,
    in_w = in_size,
    in_h = in_size,
    model_path = vis_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  return(emb)
}

#' Extract Text Embeddings with OpenAI CLIP ViT-B/32
#'
#' Computes L2-normalized 512-dimensional semantic embeddings for character strings
#' using OpenAI's CLIP (ViT-B/32) text model directly via C++ ONNX Runtime.
#'
#' @param text A character vector of one or more text prompts.
#' @param model Model name or path to a custom ONNX text model file (default `"clip-vit-b32"`).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param engine Execution engine: `"cpu"` or `"gpu"` (DirectML).
#' @param device_id GPU device ID (default -1 for auto).
#' @param verbose Logical. Whether to show progress messages.
#' @param dir Directory where models are stored (default `pliman_model_dir()`).
#' @param ... Additional arguments.
#' @return A numeric matrix of dimensions `length(text) x 512` containing the L2-normalized text embeddings.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   embs <- text_embed_dl(c("healthy soybean leaf", "asian rust lesion", "pest bite"))
#'   dim(embs) # 3 x 512
#' }
text_embed_dl <- function(text,
                          model = "clip-vit-b32",
                          threads = 0,
                          engine = c("gpu", "cpu"),
                          device_id = -1,
                          verbose = FALSE,
                          dir = pliman_model_dir(),
                          ...) {
  engine <- match.arg(engine)
  if (!is.character(text) || length(text) == 0L) {
    stop("text must be a non-empty character vector.")
  }

  model_str <- if (is.character(model)) .resolve_model_name(model[1]) else "clip-vit-b32"
  dir <- pliman_model_dir(dir)

  if (file.exists(model[1])) {
    txt_file <- normalizePath(model[1], winslash = "/")
  } else {
    pliman_download_model(model = model_str, dir = dir)
    txt_file <- file.path(dir, "clip-vit-b32-text.onnx")
  }

  vocab_file <- file.path(dir, "clip_vocab.json")
  merges_file <- file.path(dir, "clip_merges.txt")

  if (!file.exists(vocab_file) || !file.exists(merges_file)) {
    pliman_download_model(model = "clip-vit-b32", dir = dir)
  }

  lib_path <- pliman_onnx_library_path()
  use_gpu <- (engine == "gpu")

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Tokenizing prompt(s) and running CLIP text encoder [{toupper(engine)}]...",
      msg_done = "CLIP text embedding complete"
    )
  }

  tok <- .clip_tokenize(text, vocab_file, merges_file, max_len = 77L)
  batch_size <- length(text)

  embs <- run_clip_text_cpp(
    input_ids = as.integer(t(tok$input_ids)),
    batch_size = as.integer(batch_size),
    seq_len = 77L,
    attention_mask = as.integer(t(tok$attention_mask)),
    model_path = txt_file,
    lib_path = lib_path,
    num_threads = as.integer(threads),
    use_gpu = use_gpu,
    device_id = as.integer(device_id)
  )

  rownames(embs) <- text
  return(embs)
}

# ==============================================================================
# SECTION 9: YOLO-WORLD ZERO-SHOT OPEN-VOCABULARY OBJECT DETECTION
# ==============================================================================

#' Open-Vocabulary Zero-Shot Object Detection with YOLO-World
#'
#' Detects arbitrary objects described in natural language without any model retraining,
#' powered by YOLO-World (v8s) and Ultralytics.
#'
#' @param img An `Image` object, a 2D/3D numeric array, a path to an image file, or a list of images.
#' @param classes A character vector of object classes or a comma-separated string (e.g. `c("dog", "frisbee")`).
#' @param conf_threshold Minimum confidence threshold (default 0.25).
#' @param iou_threshold IoU threshold for NMS (default 0.45).
#' @param model Model name or weights path (default `"yolov8s-world"`).
#' @param rainbow Logical. If `TRUE`, gives each detected object instance a distinct, vibrant color.
#'   If `FALSE` (default), objects of the same class share the exact same color.
#' @param col Optional color or palette for bounding boxes.
#' @param lwd Bounding box line width (default 2).
#' @param cex Size factor for bounding box label text (default `1.0`).
#' @param pad Padding multiplier controlling the background badge size around label text (default `1.0`).
#' @param show_text Logical. Whether to display text labels on bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in labels (default `TRUE`).
#' @param show_class Logical. Whether to display class names in labels (default `TRUE`).
#' @param show_id Logical. Whether to display object IDs (`#1`, `#2`, ...) in labels (default `FALSE`).
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param return_features Logical. If `TRUE`, extracts deep feature representations for detected
#'   bounding boxes and the global image. Attaches `"features"` (an \eqn{N \times D} matrix of per-object
#'   deep embeddings, 512-D for CLIP or 384-D for DINOv2) and `"feature_map"` (dense spatial 3-component PCA false-color `Image`)
#'   as attributes to the returned object. Defaults to `FALSE`.
#' @param feature_model Foundation model used for feature representation when `return_features = TRUE`.
#'   Options are `"clip-vit-b32"` (default, OpenAI CLIP ViT-B/32 producing 512-D semantic embeddings) or
#'   `"dinov2"` (Meta DINOv2-S ViT producing 384-D self-supervised patch embeddings).
#' @param plot_features Logical. If `TRUE` and `return_features = TRUE`, visualizes the dense spatial PCA feature map. Defaults to `FALSE`.
#' @param plot Logical. If `TRUE` (default), displays the image with annotated bounding boxes.
#' @param auto_install Logical. If `TRUE`, sets up required Python libraries if not already available.
#' @param verbose Logical. Whether to display progress messages.
#' @param dir Directory where models are stored (default `pliman_model_dir()`).
#' @param ... Additional plotting parameters.
#' @return A data frame containing detected bounding boxes (`id`, `xmin`, `ymin`, `xmax`, `ymax`, `class`, `score`).
#'   If `return_features = TRUE`, `features` (an \eqn{N \times D} matrix) and `feature_map` (dense spatial PCA `Image`)
#'   are attached as attributes.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_import("field.jpg")
#'   dets <- image_detect_world(img, classes = c("corn plant", "weed", "soil"), rainbow = TRUE)
#' }
image_detect_world <- function(img,
                               classes = c("object"),
                               conf_threshold = 0.25,
                               iou_threshold = 0.45,
                               model = "yolov8s-world",
                               rainbow = TRUE,
                               col = NULL,
                               lwd = 2,
                               cex = 1.0,
                               pad = 1.0,
                               show_text = TRUE,
                               show_conf = TRUE,
                               show_class = TRUE,
                               show_id = FALSE,
                               badge = TRUE,
                               return_features = FALSE,
                               feature_model = c("yolo", "clip-vit-b32", "dinov2"),
                               plot_features = FALSE,
                               plot = TRUE,
                               auto_install = TRUE,
                               verbose = TRUE,
                               dir = pliman_model_dir(),
                               ...) {
  dots <- list(...)
  if ("label_size" %in% names(dots)) cex <- dots$label_size
  if ("text_size" %in% names(dots)) cex <- dots$text_size
  if ("font_scale" %in% names(dots)) cex <- dots$font_scale
  if ("cex_scale" %in% names(dots)) cex <- dots$cex_scale
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("label_box" %in% names(dots)) badge <- isTRUE(dots$label_box)

  if (length(cex) > 1L) {
    if (missing(pad)) pad <- cex[2]
    cex <- cex[1]
  }

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_detect_world(x, classes = classes, conf_threshold = conf_threshold,
                         iou_threshold = iou_threshold, model = model, rainbow = rainbow,
                         col = col, lwd = lwd, cex = cex, pad = pad, show_text = show_text,
                         show_conf = show_conf, show_class = show_class, show_id = show_id,
                         badge = badge, return_features = return_features,
                         feature_model = feature_model, plot_features = plot_features,
                         plot = FALSE, auto_install = auto_install,
                         verbose = verbose, dir = dir, ...)
    })
    return(res)
  }

  mat <- if (inherits(img, c("Image", "image"))) image_data(img) else img
  temp_in <- NULL
  if (is.character(img) && file.exists(img[1])) {
    img_path <- normalizePath(img[1], winslash = "/")
    img_obj <- image_import(img[1])
    mat <- image_data(img_obj)
  } else {
    temp_in <- file.path(tempdir(), paste0("yolow_in_", as.integer(stats::runif(1, 1000, 9999)), ".png"))
    image_export(as_image(mat), temp_in)
    img_path <- normalizePath(temp_in, winslash = "/")
  }
  on.exit(if (!is.null(temp_in) && file.exists(temp_in)) unlink(temp_in), add = TRUE)

  if (length(classes) == 1L && grepl("[,;]", classes)) {
    classes <- trimws(strsplit(classes, "[,;]+")[[1]])
  }
  classes <- classes[nzchar(classes)]
  if (length(classes) == 0L) classes <- c("object")

  dir <- pliman_model_dir(dir)
  model_str <- .resolve_model_name(model[1])
  model_file <- if (file.exists(model[1])) {
    if (grepl("\\.onnx$", model[1], ignore.case = TRUE)) {
      pt_candidate <- sub("\\.onnx$", ".pt", model[1], ignore.case = TRUE)
      if (file.exists(pt_candidate)) {
        normalizePath(pt_candidate, winslash = "/")
      } else {
        cli::cli_abort("YOLO-World open-vocabulary detection requires PyTorch weights (.pt), not ONNX (.onnx).")
      }
    } else {
      normalizePath(model[1], winslash = "/")
    }
  } else {
    pliman_download_model(model = model_str, dir = dir)
  }

  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)

  run_dir <- file.path(tempdir(), "yolow_run")
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  out_json <- file.path(run_dir, paste0("yolow_out_", as.integer(stats::runif(1, 1000, 9999)), ".json"))
  py_script <- file.path(run_dir, "yolow_predict.py")

  py_code <- paste0(
    "import sys, json\n",
    "from ultralytics import YOLOWorld\n",
    "try:\n",
    "    model = YOLOWorld('", model_file, "')\n",
    "    classes = ", jsonlite::toJSON(classes), "\n",
    "    model.set_classes(classes)\n",
    "    results = model.predict('", img_path, "', conf=", conf_threshold, ", iou=", iou_threshold, ", verbose=False)\n",
    "    boxes = []\n",
    "    for r in results:\n",
    "        for box in r.boxes:\n",
    "            xyxy = box.xyxy[0].tolist()\n",
    "            cls_id = int(box.cls[0].item())\n",
    "            conf = float(box.conf[0].item())\n",
    "            name = classes[cls_id] if cls_id < len(classes) else str(cls_id)\n",
    "            boxes.append({\n",
    "                'xmin': xyxy[0], 'ymin': xyxy[1], 'xmax': xyxy[2], 'ymax': xyxy[3],\n",
    "                'score': conf, 'class': name\n",
    "            })\n",
    "    with open('", normalizePath(out_json, winslash = "/", mustWork = FALSE), "', 'w') as f:\n",
    "        json.dump(boxes, f)\n",
    "except Exception as e:\n",
    "    sys.stderr.write(f'ERROR: {e}\\n')\n",
    "    sys.exit(1)\n"
  )
  writeLines(py_code, py_script)

  if (isTRUE(verbose)) {
    pb <- cli::cli_process_start(
      "Running YOLO-World open-vocabulary object detection...",
      on_exit = "failed"
    )
  }

  status <- system2(py_exec, args = shQuote(py_script))
  if (status != 0L || !file.exists(out_json)) {
    cli::cli_abort("YOLO-World prediction failed. Check Python Ultralytics installation.")
  }

  raw_boxes <- jsonlite::fromJSON(out_json)
  unlink(out_json)

  if (length(raw_boxes) == 0L || nrow(raw_boxes) == 0L) {
    if (isTRUE(verbose)) cli::cli_process_done(id = pb, msg_done = "No objects detected above confidence threshold.")
    empty_df <- data.frame(id = integer(0), xmin = numeric(0), ymin = numeric(0),
                           xmax = numeric(0), ymax = numeric(0), class = character(0),
                           label = character(0), score = numeric(0))
    if (isTRUE(plot)) plot(as_image(mat), ...)
    empty_df <- .attach_features(empty_df, mat = mat, boxes = empty_df,
                                 return_features = return_features, feature_model = feature_model,
                                 plot_features = plot_features, verbose = verbose, dir = dir)
    return(empty_df)
  }

  df_boxes <- data.frame(
    id = seq_len(nrow(raw_boxes)),
    xmin = round(as.numeric(raw_boxes$xmin), 1),
    ymin = round(as.numeric(raw_boxes$ymin), 1),
    xmax = round(as.numeric(raw_boxes$xmax), 1),
    ymax = round(as.numeric(raw_boxes$ymax), 1),
    class = as.character(raw_boxes$class),
    label = as.character(raw_boxes$class),
    score = round(as.numeric(raw_boxes$score), 4),
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    sum_txt <- paste0("pliman detected {.bold ", .format_detection_summary(df_boxes$class), "} in the image.")
    cli::cli_process_done(id = pb, msg_done = sum_txt)
  }

  if (isTRUE(plot)) {
    plot(as_image(mat), ...)
    num_inst <- nrow(df_boxes)
    pal <- .get_class_palette(labels = df_boxes$label, rainbow = rainbow, user_col = col)
    .plot_yolo_bboxes(
      boxes = df_boxes,
      palette_colors = pal,
      lwd = lwd,
      cex = cex,
      pad = pad,
      show_text = show_text,
      show_conf = show_conf,
      show_class = show_class,
      show_id = show_id,
      badge = badge
    )
  }

  df_boxes <- .attach_features(df_boxes, mat = mat, boxes = df_boxes,
                               return_features = return_features, feature_model = feature_model,
                               plot_features = plot_features, verbose = verbose, dir = dir)
  return(df_boxes)
}

# ==============================================================================
# SECTION 10: GENERATIVE AI DIFFUSION (SD-TURBO / LCM)
# ==============================================================================

#' Ultra-Fast Generative Image Synthesis with SD-Turbo / LCM
#'
#' Generates photo-realistic synthetic images from natural language text prompts
#' in 1-2 inference steps using Stability AI's SD-Turbo (or LCM), enhanced to high
#' resolution with Real-ESRGAN super-resolution.
#'
#' @param prompt Character text prompt describing the image to generate.
#' @param negative_prompt Optional negative prompt specifying elements to avoid.
#' @param aspect Aspect ratio of the generated image: `"square"` (1:1, 512x512),
#'   `"portrait"` (9:16, 432x768), or `"landscape"` (16:9, 768x432).
#' @param n_steps Number of diffusion inference steps (default 1 for SD-Turbo, or 2-4 for LCM).
#' @param guidance_scale Classifier-Free Guidance scale (default 0.0 for SD-Turbo).
#' @param superres Logical. If `TRUE` (default), enhances the generated image resolution
#'   using Real-ESRGAN deep learning upscaler ([image_superres_dl()]).
#' @param scale Upscaling factor when `superres = TRUE` (default 4).
#' @param seed Optional integer seed for reproducible generation.
#' @param plot Logical. If `TRUE` (default), displays the generated image.
#' @param output_dir Optional directory to save generated image file.
#' @param auto_install Logical. Whether to automatically configure environment if needed.
#' @param verbose Logical. Whether to show progress messages.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A `pliman` `Image` object.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   img <- image_generate_dl("macro photo of a healthy soybean leaf with morning dew")
#' }
image_generate_dl <- function(prompt,
                              negative_prompt = NULL,
                              aspect = c("square", "portrait", "landscape"),
                              n_steps = 1,
                              guidance_scale = 0.0,
                              superres = TRUE,
                              scale = 4,
                              seed = NULL,
                              plot = TRUE,
                              output_dir = NULL,
                              auto_install = TRUE,
                              verbose = TRUE,
                              ...) {
  if (!is.character(prompt) || length(prompt) == 0L || !nzchar(prompt[1])) {
    stop("prompt must be a non-empty character string.")
  }

  if (is.character(aspect) && length(aspect) == 1L) {
    if (aspect %in% c("1:1", "1x1")) aspect <- "square"
    if (aspect %in% c("9:16", "9x16")) aspect <- "portrait"
    if (aspect %in% c("16:9", "16x9")) aspect <- "landscape"
  }
  aspect <- match.arg(aspect, c("square", "portrait", "landscape"))
  wh <- switch(aspect,
    square = c(512L, 512L),
    portrait = c(432L, 768L),
    landscape = c(768L, 432L)
  )
  width <- wh[1]
  height <- wh[2]

  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)

  run_dir <- file.path(tempdir(), "sd_turbo_run")
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  raw_file <- file.path(run_dir, paste0("sd_gen_raw_", as.integer(stats::runif(1, 1000, 9999)), ".png"))
  py_script <- file.path(run_dir, "sd_generate.py")
  seed_py <- if (!is.null(seed)) paste0("torch.Generator(device=device).manual_seed(", as.integer(seed), ")") else "None"
  neg_py <- if (!is.null(negative_prompt)) paste0("'", gsub("'", "\\\\'", negative_prompt), "'") else "None"

  py_code <- paste0(
    "import os, sys, warnings, io, contextlib\n",
    "warnings.filterwarnings('ignore')\n",
    "os.environ['PYTHONWARNINGS'] = 'ignore'\n",
    "os.environ['HF_HUB_DISABLE_SYMLINKS_WARNING'] = '1'\n",
    "os.environ['TOKENIZERS_PARALLELISM'] = 'false'\n",
    "os.environ['TRANSFORMERS_VERBOSITY'] = 'error'\n",
    "os.environ['DIFFUSERS_VERBOSITY'] = 'error'\n",
    "\n",
    "f_null = open(os.devnull, 'w')\n",
    "with contextlib.redirect_stdout(f_null), contextlib.redirect_stderr(f_null):\n",
    "    try:\n",
    "        import torch\n",
    "        from diffusers import AutoPipelineForText2Image\n",
    "        import numpy as np\n",
    "    except ImportError:\n",
    "        import subprocess\n",
    "        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'diffusers', 'transformers', 'accelerate'])\n",
    "        import torch\n",
    "        from diffusers import AutoPipelineForText2Image\n",
    "        import numpy as np\n",
    "\n",
    "    device = 'cuda' if torch.cuda.is_available() else 'cpu'\n",
    "    dtype = torch.float16 if device == 'cuda' else torch.float32\n",
    "    pipe = AutoPipelineForText2Image.from_pretrained('stabilityai/sd-turbo', torch_dtype=dtype)\n",
    "    pipe.set_progress_bar_config(disable=True)\n",
    "    pipe.to(device)\n",
    "    if device == 'cuda':\n",
    "        pipe.vae.to(dtype=torch.float32)\n",
    "        orig_decode = pipe.vae.decode\n",
    "        def safe_decode(z, *args, **kwargs):\n",
    "            return orig_decode(z.to(torch.float32), *args, **kwargs)\n",
    "        pipe.vae.decode = safe_decode\n",
    "    generator = ", seed_py, "\n",
    "    image = pipe(\n",
    "        prompt='", gsub("'", "\\\\'", prompt[1]), "',\n",
    "        negative_prompt=", neg_py, ",\n",
    "        num_inference_steps=", as.integer(n_steps), ",\n",
    "        guidance_scale=", as.numeric(guidance_scale), ",\n",
    "        width=", as.integer(width), ",\n",
    "        height=", as.integer(height), ",\n",
    "        generator=generator\n",
    "    ).images[0]\n",
    "    if np.array(image).max() == 0 and device == 'cuda':\n",
    "        pipe32 = AutoPipelineForText2Image.from_pretrained('stabilityai/sd-turbo', torch_dtype=torch.float32).to(device)\n",
    "        pipe32.set_progress_bar_config(disable=True)\n",
    "        image = pipe32(\n",
    "            prompt='", gsub("'", "\\\\'", prompt[1]), "',\n",
    "            negative_prompt=", neg_py, ",\n",
    "            num_inference_steps=", as.integer(n_steps), ",\n",
    "            guidance_scale=", as.numeric(guidance_scale), ",\n",
    "            width=", as.integer(width), ",\n",
    "            height=", as.integer(height), ",\n",
    "            generator=generator\n",
    "        ).images[0]\n",
    "    image.save('", normalizePath(raw_file, winslash = "/", mustWork = FALSE), "')\n",
    "f_null.close()\n"
  )
  writeLines(py_code, py_script)

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Generating synthetic image with SD-Turbo [{aspect} {width}x{height}, {n_steps} step(s)]...",
      msg_done = "Synthetic image generation complete"
    )
  }

  status <- system2(py_exec, args = shQuote(py_script))
  if (status != 0L || !file.exists(raw_file)) {
    cli::cli_abort("SD-Turbo image generation failed. Please check Python diffusers installation.")
  }

  gen_img <- image_import(raw_file)

  if (isTRUE(superres)) {
    gen_img <- image_superres_dl(
      gen_img,
      scale = scale,
      plot = FALSE,
      verbose = verbose
    )
  }

  if (!is.null(output_dir) && dir.exists(output_dir)) {
    out_file <- file.path(output_dir, paste0("sd_gen_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
    image_export(gen_img, out_file)
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Image saved to {out_file}")
    }
  }

  if (isTRUE(plot)) {
    plot(gen_img, ...)
  }

  return(gen_img)
}

#' Generative AI Image Inpainting with SD-Turbo
#'
#' Modifies or replaces masked regions of an image based on a natural language text prompt
#' using SD-Turbo inpainting.
#'
#' @param img An `Image` object, matrix, or file path to the input image.
#' @param mask An `Image` object, binary matrix, or file path where foreground (1/TRUE) indicates the area to inpaint.
#' @param prompt Text prompt describing what to synthesize inside the masked area.
#' @param negative_prompt Optional negative prompt.
#' @param n_steps Number of inference steps (default 2).
#' @param guidance_scale CFG scale (default 0.0).
#' @param seed Optional random seed.
#' @param plot Logical. If `TRUE` (default), plots the inpainted image.
#' @param output_dir Optional directory to save the inpainted image.
#' @param auto_install Logical. Whether to automatically configure environment if needed.
#' @param verbose Logical. Whether to show progress messages.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A `pliman` `Image` object containing the inpainted result.
#' @export
#' @examples
#' \dontrun{
#'   library(pliman)
#'   inpainted <- image_inpaint_dl(img, mask, "a fresh green sprout")
#' }
image_inpaint_dl <- function(img,
                             mask,
                             prompt,
                             negative_prompt = NULL,
                             n_steps = 2,
                             guidance_scale = 0.0,
                             seed = NULL,
                             plot = TRUE,
                             output_dir = NULL,
                             auto_install = TRUE,
                             verbose = TRUE,
                             ...) {
  if (!is.character(prompt) || length(prompt) == 0L || !nzchar(prompt[1])) {
    stop("prompt must be a non-empty character string.")
  }

  py_exec <- .ensure_python(auto_install = auto_install)
  .setup_python_yolo_env(py_exec, auto_install = auto_install)

  run_dir <- file.path(tempdir(), "sd_inpaint_run")
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  # Export img and mask to temp PNG
  in_img_file <- file.path(run_dir, "init_img.png")
  in_mask_file <- file.path(run_dir, "init_mask.png")
  out_file <- if (!is.null(output_dir) && dir.exists(output_dir)) {
    file.path(output_dir, paste0("sd_inpaint_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
  } else {
    file.path(run_dir, paste0("sd_inpaint_", as.integer(stats::runif(1, 1000, 9999)), ".png"))
  }

  img_obj <- if (is.character(img) && file.exists(img[1])) image_import(img[1]) else as_image(img)
  mask_obj <- if (is.character(mask) && file.exists(mask[1])) image_import(mask[1]) else as_image(mask)
  image_export(img_obj, in_img_file)
  image_export(mask_obj, in_mask_file)

  py_script <- file.path(run_dir, "sd_inpaint.py")
  seed_py <- if (!is.null(seed)) paste0("torch.Generator(device=device).manual_seed(", as.integer(seed), ")") else "None"
  neg_py <- if (!is.null(negative_prompt)) paste0("'", gsub("'", "\\\\'", negative_prompt), "'") else "None"

  py_code <- paste0(
    "import os, sys, warnings, io, contextlib\n",
    "warnings.filterwarnings('ignore')\n",
    "os.environ['PYTHONWARNINGS'] = 'ignore'\n",
    "os.environ['HF_HUB_DISABLE_SYMLINKS_WARNING'] = '1'\n",
    "os.environ['TOKENIZERS_PARALLELISM'] = 'false'\n",
    "os.environ['TRANSFORMERS_VERBOSITY'] = 'error'\n",
    "os.environ['DIFFUSERS_VERBOSITY'] = 'error'\n",
    "\n",
    "f_null = open(os.devnull, 'w')\n",
    "with contextlib.redirect_stdout(f_null), contextlib.redirect_stderr(f_null):\n",
    "    try:\n",
    "        import torch\n",
    "        from PIL import Image\n",
    "        from diffusers import AutoPipelineForInpainting\n",
    "        import numpy as np\n",
    "    except ImportError:\n",
    "        import subprocess\n",
    "        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'diffusers', 'transformers', 'accelerate'])\n",
    "        import torch\n",
    "        from PIL import Image\n",
    "        from diffusers import AutoPipelineForInpainting\n",
    "        import numpy as np\n",
    "\n",
    "    device = 'cuda' if torch.cuda.is_available() else 'cpu'\n",
    "    dtype = torch.float16 if device == 'cuda' else torch.float32\n",
    "    pipe = AutoPipelineForInpainting.from_pretrained('stabilityai/sd-turbo', torch_dtype=dtype)\n",
    "    pipe.set_progress_bar_config(disable=True)\n",
    "    pipe.to(device)\n",
    "    if device == 'cuda':\n",
    "        pipe.vae.to(dtype=torch.float32)\n",
    "        orig_decode = pipe.vae.decode\n",
    "        def safe_decode(z, *args, **kwargs):\n",
    "            return orig_decode(z.to(torch.float32), *args, **kwargs)\n",
    "        pipe.vae.decode = safe_decode\n",
    "    init_img = Image.open('", normalizePath(in_img_file, winslash = "/"), "').convert('RGB').resize((512, 512))\n",
    "    mask_img = Image.open('", normalizePath(in_mask_file, winslash = "/"), "').convert('L').resize((512, 512))\n",
    "    generator = ", seed_py, "\n",
    "    image = pipe(\n",
    "        prompt='", gsub("'", "\\\\'", prompt[1]), "',\n",
    "        negative_prompt=", neg_py, ",\n",
    "        image=init_img,\n",
    "        mask_image=mask_img,\n",
    "        num_inference_steps=", as.integer(n_steps), ",\n",
    "        guidance_scale=", as.numeric(guidance_scale), ",\n",
    "        generator=generator\n",
    "    ).images[0]\n",
    "    if np.array(image).max() == 0 and device == 'cuda':\n",
    "        pipe32 = AutoPipelineForInpainting.from_pretrained('stabilityai/sd-turbo', torch_dtype=torch.float32).to(device)\n",
    "        pipe32.set_progress_bar_config(disable=True)\n",
    "        image = pipe32(\n",
    "            prompt='", gsub("'", "\\\\'", prompt[1]), "',\n",
    "            negative_prompt=", neg_py, ",\n",
    "            image=init_img,\n",
    "            mask_image=mask_img,\n",
    "            num_inference_steps=", as.integer(n_steps), ",\n",
    "            guidance_scale=", as.numeric(guidance_scale), ",\n",
    "            generator=generator\n",
    "        ).images[0]\n",
    "    image.save('", normalizePath(out_file, winslash = "/", mustWork = FALSE), "')\n",
    "f_null.close()\n"
  )
  writeLines(py_code, py_script)

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Running generative inpainting with SD-Turbo [{n_steps} step(s)]...",
      msg_done = "Inpainting complete"
    )
  }

  status <- system2(py_exec, args = shQuote(py_script))
  if (status != 0L || !file.exists(out_file)) {
    cli::cli_abort("SD-Turbo inpainting failed. Please check Python diffusers installation.")
  }

  res_img <- image_import(out_file)

  if (isTRUE(plot)) {
    plot(res_img, ...)
  }

  return(res_img)
}


