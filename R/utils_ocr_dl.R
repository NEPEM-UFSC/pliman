# ==============================================================================
# SECTION 11: FIELD OPTICAL CHARACTER RECOGNITION (OCR)
# ==============================================================================

# Internal registry of OCR backbones: local file candidates (in `dir`), catalog name
# used by `pliman_download_model()` and, for CTC recognisers, the character dictionary.
.ocr_det_specs <- list(
  ppocrv5_det        = list(files = "PP-OCRv5_mobile_det.onnx", catalog = "ppocrv5_det", tag = "PP-OCRv5"),
  ppocrv5_det_server = list(files = "PP-OCRv5_server_det.onnx", catalog = "ppocrv5_det_server", tag = "PP-OCRv5-Server"),
  ppocr_det          = list(files = c("ppocr_det.onnx", "ch_PP-OCRv4_det.onnx"), catalog = "ppocr_det", tag = "PP-OCRv4"),
  ppocr_det_server   = list(files = c("ppocr_det_server.onnx", "ch_PP-OCRv4_det_server.onnx"), catalog = "ppocr_det_server", tag = "PP-OCRv4-Server"),
  craft              = list(files = "craft.onnx", catalog = "craft", tag = "CRAFT")
)

.ocr_rec_specs <- list(
  ppocrv5_rec        = list(files = "PP-OCRv5_mobile_rec.onnx", catalog = "ppocrv5_rec", dict = "ppocrv5_dict.txt", tag = "PP-OCRv5"),
  ppocrv5_rec_server = list(files = "PP-OCRv5_server_rec.onnx", catalog = "ppocrv5_rec_server", dict = "ppocrv5_dict.txt", tag = "PP-OCRv5-Server"),
  ppocrv5_rec_latin  = list(files = "latin_PP-OCRv5_mobile_rec.onnx", catalog = "ppocrv5_rec_latin", dict = "ppocrv5_latin_dict.txt", tag = "PP-OCRv5-Latin"),
  ppocr_rec          = list(files = c("ppocr_rec.onnx", "ch_PP-OCRv4_rec.onnx"), catalog = "ppocr_rec", dict = "ppocr_dict.txt", tag = "PP-OCRv4")
)

.ocr_normalize_key <- function(x, specs) {
  k <- tolower(gsub("-", "_", x))
  k <- sub("\\.onnx$", "", k)
  map <- c(
    "ppocrv5" = "ppocrv5_det", "pp_ocrv5_mobile_det" = "ppocrv5_det", "ppocrv5_det_mobile" = "ppocrv5_det",
    "pp_ocrv5_server_det" = "ppocrv5_det_server", "ppocrv5_server" = "ppocrv5_det_server",
    "pp_ocrv5_mobile_rec" = "ppocrv5_rec", "ppocrv5_rec_mobile" = "ppocrv5_rec",
    "pp_ocrv5_server_rec" = "ppocrv5_rec_server",
    "latin_pp_ocrv5_mobile_rec" = "ppocrv5_rec_latin", "ppocrv5_latin" = "ppocrv5_rec_latin",
    "ppocr" = "ppocr_det", "ch_pp_ocrv4_det" = "ppocr_det", "ppocr_server" = "ppocr_det_server",
    "ch_pp_ocrv4_det_server" = "ppocr_det_server", "ch_pp_ocrv4_rec" = "ppocr_rec",
    "craft_mlt" = "craft"
  )
  if (k %in% names(map)) k <- unname(map[k])
  if (k %in% names(specs)) k else NA_character_
}

# Resolves an OCR model to a local file: user path > exact local file > download.
# Exact matching avoids the fuzzy lookup swapping e.g. *_mobile_rec for *_mobile_det.
.ocr_model_path <- function(name, specs, dir, verbose) {
  if (file.exists(name) && !dir.exists(name)) return(normalizePath(name, winslash = "/"))
  key <- .ocr_normalize_key(name, specs)
  if (is.na(key)) {
    p <- .find_local_model(name, dir = dir)
    if (is.null(p)) p <- pliman_download_model(name, dir = dir, verbose = verbose)
    return(p)
  }
  for (cand in specs[[key]]$files) {
    p <- file.path(pliman_model_dir(dir), cand)
    if (file.exists(p) && file.info(p)$size > 1000) return(normalizePath(p, winslash = "/"))
  }
  pliman_download_model(specs[[key]]$catalog, dir = dir, verbose = verbose)
}

# Loads a CTC character dictionary: a user path, the copy in the model directory or,
# failing both, a download from `pliman_models_base_url()` (catalog entry by file name).
.ocr_load_dict <- function(dict_file, dir, verbose = TRUE) {
  path <- NULL
  if (file.exists(dict_file) && !dir.exists(dict_file)) {
    path <- dict_file
  } else {
    in_dir <- file.path(pliman_model_dir(dir), basename(dict_file))
    if (file.exists(in_dir)) {
      path <- in_dir
    } else {
      path <- tryCatch(
        pliman_download_model(basename(dict_file), dir = dir, verbose = verbose),
        error = function(e) NULL
      )
    }
  }
  if (is.null(path) || !file.exists(path)) {
    cli::cli_alert_warning("Character dictionary {.file {basename(dict_file)}} not found; recognized text will be empty.")
    return(character(0))
  }
  con <- file(path, encoding = "UTF-8")
  on.exit(close(con))
  dict <- readLines(con, warn = FALSE)
  dict <- sub("\r$", "", dict)
  # PaddleOCR appends the space character as the last class (use_space_char = TRUE).
  if (!" " %in% dict) dict <- c(dict, " ")
  dict
}

.ocr_allowed_ids <- function(charset, dict) {
  if (is.null(charset) || length(charset) == 0L || length(dict) == 0L) return(NULL)
  cs <- paste(charset, collapse = "")
  if (grepl("^\\[.*\\]$", cs)) {
    ok <- grepl(paste0("^", cs, "$"), dict, perl = TRUE)
  } else {
    ok <- dict %in% unique(strsplit(cs, "")[[1]])
  }
  if (!any(ok)) {
    cli::cli_alert_warning("No dictionary character matches {.val {cs}}; ignoring {.arg charset}.")
    return(NULL)
  }
  as.integer(which(ok))
}

# Side lengths of TL-TR-BR-BL quads (N x 8 matrix).
.ocr_quad_wh <- function(q) {
  w <- pmax(sqrt((q[, 3] - q[, 1])^2 + (q[, 4] - q[, 2])^2),
            sqrt((q[, 5] - q[, 7])^2 + (q[, 6] - q[, 8])^2))
  h <- pmax(sqrt((q[, 7] - q[, 1])^2 + (q[, 8] - q[, 2])^2),
            sqrt((q[, 5] - q[, 3])^2 + (q[, 6] - q[, 4])^2))
  list(w = pmax(1, w), h = pmax(1, h))
}

# Converts oriented component geometry (map grid) into image-space quads.
# `expand` is the per-side expansion in map pixels; `sx`, `sy` map -> image scale.
.ocr_comps_to_quads <- function(cp, expand, sx, sy) {
  u0 <- cp$u_min - 0.5 - expand
  u1 <- cp$u_max + 0.5 + expand
  v0 <- cp$v_min - 0.5 - expand
  v1 <- cp$v_max + 0.5 + expand
  vx <- -cp$uy; vy <- cp$ux
  px <- function(u, v) ((cp$cx + u * cp$ux + v * vx) + 0.5) * sx - 0.5
  py <- function(u, v) ((cp$cy + u * cp$uy + v * vy) + 0.5) * sy - 0.5
  cbind(
    x1 = px(u0, v0), y1 = py(u0, v0),
    x2 = px(u1, v0), y2 = py(u1, v0),
    x3 = px(u1, v1), y3 = py(u1, v1),
    x4 = px(u0, v1), y4 = py(u0, v1)
  )
}

# Text detection: resize to the network scale, run DBNet / CRAFT and return
# oriented quads (0-based continuous image coordinates) with detection scores.
.ocr_detect <- function(mat, pix_scale, orig_w, orig_h, is_craft, det_file, det_limit,
                        det_thresh, box_thresh, unclip_ratio, text_threshold,
                        link_threshold, low_text, min_size, lib_path, threads, use_gpu) {
  r <- det_limit / max(orig_w, orig_h)
  rw <- as.integer(max(32, round(orig_w * r / 32) * 32))
  rh <- as.integer(max(32, round(orig_h * r / 32) * 32))
  full <- matrix(c(-0.5, -0.5, orig_w - 0.5, -0.5, orig_w - 0.5, orig_h - 0.5, -0.5, orig_h - 0.5), nrow = 1)
  # PP-OCR detectors were trained on BGR (cv2) inputs; CRAFT on RGB.
  tensor <- ocr_warp_crops_cpp(
    mat, full, out_h = rh, out_w = rw, pad_w = rw,
    mean = c(0.485, 0.456, 0.406), std = c(0.229, 0.224, 0.225),
    bgr = !is_craft, scale = pix_scale
  )
  sx <- orig_w / rw
  sy <- orig_h / rh

  if (is_craft) {
    craft_res <- run_craft_det_cpp(
      tensor_vec = tensor, in_w = rw, in_h = rh, model_path = det_file,
      lib_path = lib_path, num_threads = as.integer(threads), use_gpu = use_gpu
    )
    score_text <- craft_res$score_text
    mask <- (score_text >= low_text) | (craft_res$score_link >= link_threshold)
    lab <- bwlabel_cpp(mask)
    cp <- ocr_components_cpp(lab, score_text, max(lab))
    # Official CRAFT filters (size >= 10 on the half-resolution map).
    cp <- cp[cp$n >= 40 & cp$max_score >= text_threshold, , drop = FALSE]
    if (nrow(cp) == 0L) return(NULL)
    bw <- cp$xmax - cp$xmin + 1
    bh <- cp$ymax - cp$ymin + 1
    expand <- 4 * sqrt(cp$n * pmin(bw, bh) / (2 * bw * bh))
    det_score <- cp$max_score
  } else {
    prob_map <- run_ocr_det_cpp(
      tensor_vec = tensor, in_w = rw, in_h = rh, model_path = det_file,
      lib_path = lib_path, num_threads = as.integer(threads), use_gpu = use_gpu
    )
    lab <- bwlabel_cpp(prob_map >= det_thresh)
    cp <- ocr_components_cpp(lab, prob_map, max(lab))
    cp <- cp[cp$n >= 4 & cp$mean_score >= box_thresh, , drop = FALSE]
    if (nrow(cp) == 0L) return(NULL)
    # DBNet unclip (Vatti offset of the shrunk text kernel): d = A * r / L
    L <- cp$u_max - cp$u_min + 1
    H <- cp$v_max - cp$v_min + 1
    expand <- (L * H * unclip_ratio) / (2 * (L + H))
    det_score <- cp$mean_score
  }

  q <- .ocr_comps_to_quads(cp, expand, sx, sy)
  wh <- .ocr_quad_wh(q)
  keep <- pmin(wh$w, wh$h) >= min_size
  if (!any(keep)) return(NULL)
  list(quads = q[keep, , drop = FALSE], det_score = det_score[keep],
       angle = atan2(cp$uy[keep], cp$ux[keep]) * 180 / pi,
       det_size = c(rw, rh))
}

#' Field Optical Character Recognition (OCR)
#'
#' Agronomic trials and breeding nurseries routinely track experimental entries using
#' physical plot stakes, barcode tags, seed packet labels, and ear tags in livestock.
#' In `pliman`, `image_ocr_dl()` integrates deep neural text detection (DBNet / CRAFT)
#' and text recognition (CRNN/SVTR / TrOCR) models to automatically extract text strings
#' and their spatial coordinates directly from field imagery.
#'
#' @details
#' The OCR pipeline operates in modular stages:
#' \enumerate{
#'   \item \strong{Text Detection (\code{det_model})}: Identifies oriented text regions
#'         on the image. The image is first resized so that its longest side equals
#'         \code{det_limit} (the scale the detectors were trained on), which is both faster
#'         and more accurate than running on full-resolution photographs.
#'         \itemize{
#'           \item \code{"ppocrv5_det"}: PP-OCRv5 Mobile DBNet (default, ~4.8 MB). Trained with
#'                 printed and handwritten text; fast and accurate for labels, tags and notes.
#'           \item \code{"ppocrv5_det_server"}: PP-OCRv5 Server DBNet (~88 MB). Highest detection
#'                 quality for dense or difficult handwriting, distant plot markers, or cluttered backgrounds.
#'           \item \code{"ppocr_det"}: PP-OCRv4 Mobile DBNet (~4.7 MB, legacy).
#'           \item \code{"ppocr_det_server"}: PP-OCRv4 Server DBNet (~113 MB, legacy).
#'           \item \code{"craft"}: Clova AI CRAFT (Character Region Awareness for Text Detection, ~80 MB).
#'                 Uses dual Gaussian heatmaps (character center score + inter-character affinity link).
#'                 Useful for irregular or curved scene text.
#'           \item \code{"trocr_small"}: Shortcut convenience alias. Automatically pairs the \code{"craft"}
#'                 detector with the TrOCR recognition model.
#'         }
#'         Each detected region is returned as a rotated quadrilateral (columns \code{x1..y4}),
#'         so tilted handwriting is rectified before recognition instead of being cropped together
#'         with neighbouring lines.
#'   \item \strong{Text-line orientation (\code{textline_ori})}: A tiny PP-LCNet classifier
#'         (~1 MB) detects upside-down (180 degrees) text crops and flips them before recognition.
#'   \item \strong{Text Recognition (\code{rec_model})}: Transcribes cropped image patches into character strings.
#'         Crops are rendered at 48 px height with their native aspect ratio (dynamic width, so long
#'         lines are not squeezed) and processed in batches of \code{rec_batch}.
#'         \itemize{
#'           \item \code{"ppocrv5_rec"}: PP-OCRv5 Mobile recognition (~16 MB). Single model for
#'                 printed and handwritten Chinese, English, and Latin-script text.
#'           \item \code{"ppocrv5_rec_server"}: PP-OCRv5 Server recognition (~85 MB). Best accuracy on
#'                 handwriting; slower.
#'           \item \code{"ppocrv5_rec_latin"}: PP-OCRv5 Latin Mobile recognition (~8 MB). Covers Latin-script
#'                 languages with diacritics (e.g., Portuguese \code{"\u00e7"}, \code{"\u00e3"}, \code{"\u00e9"}).
#'           \item \code{"ppocr_rec"}: PP-OCRv4 SVTR/LCNet model (~12 MB, legacy). Optimized for
#'                 printed typography and alphanumeric entry codes.
#'           \item \code{"trocr_small"}: Microsoft TrOCR Small (encoder + autoregressive
#'                 RoBERTa decoder). Transformer-based model trained on English handwritten text
#'                 lines (IAM); it tends to hallucinate words on field photographs and short codes.
#'         }
#'         When \code{rec_model} is not supplied, \code{language = "latin"} selects
#'         \code{"ppocrv5_rec_latin"} and \code{language = "english"} selects \code{"ppocrv5_rec"}.
#' }
#'
#' \strong{Restricting the character set}:
#' For plot stakes and seed packets the set of possible characters is often known (e.g. only
#' capital letters, digits, and hyphens). Supplying \code{charset} (e.g. \code{"[0-9A-Z-]"} or
#' \code{"0123456789G"}) masks every other character during CTC decoding, which removes
#' confusions such as \code{"G196"} being read as \code{"6196"} or \code{"Gi96"}.
#'
#' \strong{Reading Order Sorting}:
#' When \code{sort = TRUE} (the default), detected boxes are organized using spatial vertical IoU
#' line-clustering (> 0.35 overlap threshold) to group words belonging to the same text line,
#' followed by left-to-right sorting within each line, reproducing natural human reading order.
#'
#' \strong{Hardware}: on GPUs, ONNX Runtime may need a warm-up call per new input shape.
#' Speed gains from \code{engine = "gpu"} depend on a compatible CUDA device; on small images
#' the CPU engine can be equally fast.
#'
#' \strong{S3 Methods}:
#' The returned object of class \code{"pliman_ocr"} supports:
#' \itemize{
#'   \item \code{print(x)}: Displays a clean tabular summary with line assignments and coordinates.
#'   \item \code{plot(x)}: Renders the annotated image with colored (rotated) boxes and badges.
#'   \item \code{as.character(x)}: Returns the complete transcribed text with newline (\code{"\n"}) breaks.
#' }
#'
#' @references
#' \itemize{
#'   \item Liao, M., Wan, Z., Yao, C., Chen, K., & Bai, X. (2020).
#'   Real-time scene text detection with differentiable binarization.
#'   \emph{Proceedings of the AAAI Conference on Artificial Intelligence}, 34(07), 11474--11481.
#'   \item Du, Y., Li, C., Guo, R., Yin, X., Liu, W., Zhou, J., ... & Hu, B. (2020).
#'   PP-OCR: A practical ultra lightweight OCR system.
#'   \emph{arXiv preprint arXiv:2009.09941}.
#'   \item Cui, C., Sun, T., Lin, M., et al. (2025).
#'   PaddleOCR 3.0 Technical Report. \emph{arXiv preprint arXiv:2507.05595}.
#'   \item Baek, Y., Lee, B., Han, D., Yun, S., & Lee, H. (2019).
#'   Character region awareness for text detection.
#'   \emph{Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)}, 9365--9374.
#'   \item Li, M., Lv, T., Chen, J., Cui, L., Lu, Y., Florencio, D., Zhang, C., & Wei, F. (2023).
#'   TrOCR: Transformer-based optical character recognition with pre-trained models.
#'   \emph{AAAI Conference on Human Computation and Crowdsourcing}.
#' }
#'
#' @param img An `Image` object or character string path to an image file.
#' @param det_model Text detection backbone model. Options are:
#'   \itemize{
#'     \item \code{"ppocrv5_det"}: PP-OCRv5 Mobile DBNet (default, fast, printed and handwritten text).
#'     \item \code{"ppocrv5_det_server"}: PP-OCRv5 Server DBNet (highest accuracy).
#'     \item \code{"ppocr_det"}: PP-OCRv4 Mobile DBNet (legacy).
#'     \item \code{"ppocr_det_server"}: PP-OCRv4 Server DBNet (legacy).
#'     \item \code{"craft"}: Clova AI CRAFT (irregular/curved scene text).
#'     \item \code{"trocr_small"}: Alias running CRAFT detection with TrOCR Small recognition.
#'   }
#'   A path to a compatible DBNet ONNX file is also accepted.
#' @param rec_model Text recognition backbone model:
#'   \itemize{
#'     \item \code{"ppocrv5_rec_latin"}: PP-OCRv5 Latin recognition (default when \code{language = "latin"}).
#'     \item \code{"ppocrv5_rec"}: PP-OCRv5 Mobile recognition (default when \code{language = "english"}).
#'     \item \code{"ppocrv5_rec_server"}: PP-OCRv5 Server recognition (best for difficult handwriting).
#'     \item \code{"ppocr_rec"}: PP-OCRv4 recognition (legacy; printed text, numbers, codes).
#'     \item \code{"trocr_small"}: TrOCR Small encoder-decoder (English handwriting lines).
#'   }
#'   A path to a compatible CTC recognition ONNX file is also accepted (see \code{rec_dict}).
#' @param language Language family used to pick the default recognizer (and the TrOCR variant):
#'   \code{"latin"} (default; Latin/Portuguese text with diacritics) or \code{"english"}.
#' @param det_thresh Pixel binarization threshold for DBNet text polygon detection (default `0.30`).
#' @param box_thresh Minimum mean probability for retaining DBNet text regions (default `0.60`).
#' @param unclip_ratio Polygon expansion ratio to prevent cropping character boundaries (DBNet models,
#'   default `1.50`). Larger values (1.6--2.0) help handwriting with long ascenders/descenders.
#' @param text_threshold Character center confidence threshold for CRAFT detector (default `0.70`).
#' @param link_threshold Character affinity link threshold for CRAFT detector (default `0.40`).
#' @param low_text Lower bound text score threshold for CRAFT detector (default `0.40`).
#' @param min_size Minimum text box height in pixels to filter out noise (default `5L`).
#' @param min_score Minimum recognition confidence score to retain text boxes (default `0.0`).
#' @param sort Logical. Whether to sort detected text regions in natural reading order
#'   (top-to-bottom, line-by-line, left-to-right). Default is \code{TRUE}.
#' @param col Border color for recognized text bounding polygons (default `"forestgreen"`).
#' @param lwd Line width for rendering detected text polygon borders (default `2L`).
#' @param cex Character expansion factor for recognized text labels (default `0.80`).
#' @param badge Logical. Whether to render filled background badges behind transcribed text (default `TRUE`).
#' @param labels Label type to display on the plot: \code{"text"} (default), \code{"id"}, \code{"both"}, or \code{"none"}.
#' @param threads Number of CPU execution threads (default `0` for auto-detect).
#' @param engine Compute engine: `"gpu"` (default) or `"cpu"`.
#' @param plot Logical. Whether to render the visual OCR annotation plot (default `TRUE`).
#' @param dir Model directory (default `pliman_model_dir()`).
#' @param verbose Logical. Whether to show progress messages (default `TRUE`).
#' @param det_limit Longest image side (pixels) used for text detection. The image is resized
#'   (up or down) to this size before detection. Default `NULL` uses `960` for DBNet models and
#'   `1280` for CRAFT. Increase it (e.g. `1600`--`2000`) for very small text in large photographs.
#' @param textline_ori Logical. Whether to run the PP-LCNet text-line orientation classifier to
#'   fix upside-down crops (default `TRUE`). Ignored for TrOCR.
#' @param charset Optional allow-list of characters for CTC recognizers, either as a plain string
#'   of characters (e.g. `"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ-"`) or as a regular-expression
#'   bracket expression (e.g. `"[0-9A-Z-]"`). Default `NULL` (no restriction).
#' @param rec_dict Optional path to the character dictionary (one character per line) of a custom
#'   CTC recognition model. Default `NULL` uses the dictionary that matches `rec_model`.
#' @param rec_batch Number of text crops recognized per inference batch (default `8L`).
#' @param keep_empty Logical. Whether to keep detected regions in which no character was
#'   recognized (default `FALSE`).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return A `data.frame` of class `c("pliman_ocr", "data.frame")` containing:
#' \itemize{
#'   \item \code{id}: Sequential text region identifier.
#'   \item \code{line}: Visual text line index (when \code{sort = TRUE}).
#'   \item \code{text}: Transcribed character string.
#'   \item \code{score}: Recognition confidence score (0 to 1).
#'   \item \code{det_score}: Detection confidence score.
#'   \item \code{xmin}, \code{ymin}, \code{xmax}, \code{ymax}: Axis-aligned bounding box coordinates.
#'   \item \code{x1}, \code{y1}, ..., \code{x4}, \code{y4}: Corners (top-left, top-right, bottom-right,
#'         bottom-left) of the rotated text quadrilateral.
#'   \item \code{angle}: Text direction angle in degrees.
#' }
#' With attributes \code{"image"} (original image) and \code{"timing"} (seconds spent in each stage).
#' @export
#' @examplesIf interactive()
#'   library(pliman)
#'   img <- image_import("flaxtext")
#'   # Default: PP-OCRv5 detection + Latin recognition
#'   ocr_res <- image_ocr_dl(img)
#'   print(ocr_res)
#'
#'   # Difficult handwriting: server models
#'   ocr_hw <- image_ocr_dl(img, det_model = "ppocrv5_det_server", rec_model = "ppocrv5_rec_server")
#'
#'   # Seed packet codes: restrict the alphabet
#'   ocr_code <- image_ocr_dl(img, charset = "[0-9A-Z-]")
#'
#'   # Convert to plain text string with line breaks
#'   cat(as.character(ocr_res))
image_ocr_dl <- function(img,
                         det_model = c("ppocrv5_det", "ppocrv5_det_server", "ppocr_det", "ppocr_det_server", "craft", "trocr_small"),
                         rec_model = c("ppocrv5_rec_latin", "ppocrv5_rec", "ppocrv5_rec_server", "ppocr_rec", "trocr_small"),
                         language = c("latin", "english"),
                         det_thresh = 0.30,
                         box_thresh = 0.60,
                         unclip_ratio = 1.50,
                         text_threshold = 0.70,
                         link_threshold = 0.40,
                         low_text = 0.40,
                         min_size = 5L,
                         min_score = 0.0,
                         sort = TRUE,
                         col = "forestgreen",
                         lwd = 2L,
                         cex = 0.80,
                         badge = TRUE,
                         labels = c("text", "id", "both", "none"),
                         threads = 0,
                         engine = c("gpu", "cpu"),
                         plot = TRUE,
                         dir = pliman_model_dir(),
                         verbose = TRUE,
                         det_limit = NULL,
                         textline_ori = TRUE,
                         charset = NULL,
                         rec_dict = NULL,
                         rec_batch = 8L,
                         keep_empty = FALSE,
                         ...) {
  engine <- match.arg(engine)
  language <- if (missing(language) || is.null(language)) "latin" else match.arg(language)
  if (missing(labels) || is.null(labels) || !is.character(labels) || length(labels) == 0L) {
    labels <- "text"
  } else {
    labels <- labels[1L]
    if (!labels %in% c("text", "id", "both", "none")) labels <- "text"
  }
  user_specified_rec <- !missing(rec_model) && is.character(rec_model) && length(rec_model) == 1L && nzchar(rec_model)
  if (is.character(det_model) && length(det_model) > 1L) {
    det_model <- det_model[1L]
  }
  if (!user_specified_rec) {
    rec_model <- if (language == "latin") "ppocrv5_rec_latin" else "ppocrv5_rec"
  } else if (is.character(rec_model) && length(rec_model) > 1L) {
    rec_model <- rec_model[1L]
  }

  if (grepl("trocr", tolower(det_model))) {
    det_model <- "craft"
    if (!user_specified_rec) {
      if (isTRUE(verbose)) {
        rec_desc <- if (language == "latin") "TrOCR Latin (Portuguese/Latin handwriting)" else "TrOCR English"
        cli::cli_alert_info("TrOCR is an encoder-decoder recognition model. Using 'craft' for text detection and {rec_desc} for recognition.")
      }
      rec_model <- if (language == "latin") "trocr_lat" else "trocr"
    } else {
      rec_model <- rec_model[1L]
    }
  }

  if (is.list(img) && !inherits(img, c("Image", "image"))) {
    res <- lapply(img, function(x) {
      image_ocr_dl(
        x, det_model = det_model, rec_model = rec_model, language = language,
        det_thresh = det_thresh, box_thresh = box_thresh,
        unclip_ratio = unclip_ratio,
        text_threshold = text_threshold, link_threshold = link_threshold,
        low_text = low_text, min_size = min_size,
        min_score = min_score, sort = sort,
        col = col, lwd = lwd, cex = cex, badge = badge, labels = labels,
        threads = threads, engine = engine, plot = FALSE,
        dir = dir, verbose = verbose,
        det_limit = det_limit, textline_ori = textline_ori, charset = charset,
        rec_dict = rec_dict, rec_batch = rec_batch, keep_empty = keep_empty, ...
      )
    })
    return(res)
  }

  img_obj <- if (inherits(img, c("Image", "image"))) img else image_import(img[1])
  mat <- image_data(img_obj)
  dims <- dim(mat)
  if (is.null(dims) || length(dims) < 2) {
    stop("img must be a 2D matrix, 3D array, or an 'image' object.")
  }
  orig_w <- dims[1]
  orig_h <- dims[2]
  # Multiplier mapping stored values to the 0-255 range.
  pix_scale <- if (is.raw(mat)) 1.0 else (if (max(mat[1:min(1000, length(mat))], na.rm = TRUE) > 1.5) 1.0 else 255.0)

  is_craft <- grepl("craft", tolower(det_model))
  is_trocr <- grepl("trocr", tolower(rec_model))
  is_latin_trocr <- is_trocr && (language == "latin" || grepl("lat", tolower(rec_model)))
  if (is.null(det_limit)) det_limit <- if (is_craft) 1280L else 960L

  det_file <- .ocr_model_path(det_model, .ocr_det_specs, dir, verbose)
  det_key <- .ocr_normalize_key(det_model, .ocr_det_specs)

  if (is_trocr) {
    enc_name <- if (is_latin_trocr) "trocr_encoder_lat" else "trocr_encoder"
    dec_name <- if (is_latin_trocr) "trocr_decoder_lat" else "trocr_decoder"
    voc_name <- if (is_latin_trocr) "trocr_vocab_lat.rds" else "trocr_vocab.rds"
    voc_model <- if (is_latin_trocr) "trocr_vocab_lat" else "trocr_vocab"

    trocr_enc_file <- .find_local_model(enc_name, dir = dir)
    if (is.null(trocr_enc_file)) trocr_enc_file <- pliman_download_model(enc_name, dir = dir, verbose = verbose)
    trocr_dec_file <- .find_local_model(dec_name, dir = dir)
    if (is.null(trocr_dec_file)) trocr_dec_file <- pliman_download_model(dec_name, dir = dir, verbose = verbose)
    trocr_vocab_file <- file.path(dir, voc_name)
    if (!file.exists(trocr_vocab_file)) {
      pkg_vocab <- system.file(voc_name, package = "pliman")
      if (nzchar(pkg_vocab) && file.exists(pkg_vocab)) {
        trocr_vocab_file <- pkg_vocab
      } else {
        trocr_vocab_file <- pliman_download_model(voc_model, dir = dir, verbose = verbose)
      }
    }
    trocr_vocab <- if (file.exists(trocr_vocab_file)) readRDS(trocr_vocab_file) else character(0)
    if (!is.null(charset) && isTRUE(verbose)) {
      cli::cli_alert_warning("{.arg charset} is only supported by CTC recognizers (PP-OCR); ignored for TrOCR.")
    }
  } else {
    rec_file <- .ocr_model_path(rec_model, .ocr_rec_specs, dir, verbose)
    rec_key <- .ocr_normalize_key(rec_model, .ocr_rec_specs)
    dict_file <- if (!is.null(rec_dict)) rec_dict else if (!is.na(rec_key)) .ocr_rec_specs[[rec_key]]$dict else "ppocrv5_dict.txt"
    dict <- .ocr_load_dict(dict_file, dir, verbose)
    allowed_ids <- .ocr_allowed_ids(charset, dict)
  }

  ori_file <- NULL
  if (isTRUE(textline_ori) && !is_trocr) {
    ori_file <- tryCatch(
      .ocr_model_path("ppocr_textline_ori", list(ppocr_textline_ori = list(
        files = "PP-LCNet_x0_25_textline_ori.onnx", catalog = "ppocr_textline_ori")), dir, verbose),
      error = function(e) {
        if (isTRUE(verbose)) cli::cli_alert_warning("Text-line orientation model unavailable; skipping orientation correction.")
        NULL
      }
    )
  }

  lib_path <- pliman_onnx_lib_path(engine = engine)
  use_gpu <- (engine == "gpu")

  det_tag <- if (!is.na(det_key)) .ocr_det_specs[[det_key]]$tag else basename(det_file)
  rec_tag <- if (is_trocr) {
    if (is_latin_trocr) "TrOCR Latin" else "TrOCR English"
  } else if (!is.na(rec_key)) .ocr_rec_specs[[rec_key]]$tag else basename(rec_file)
  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = paste0("Running ", det_tag, " text detection & ", rec_tag, " recognition [{toupper(engine)}]..."),
      msg_done = "OCR inference complete"
    )
  }

  timing <- c(detection = 0, orientation = 0, recognition = 0)
  t0 <- proc.time()[["elapsed"]]
  det <- .ocr_detect(
    mat, pix_scale, orig_w, orig_h, is_craft, det_file, det_limit,
    det_thresh, box_thresh, unclip_ratio, text_threshold, link_threshold,
    low_text, min_size, lib_path, threads, use_gpu
  )
  timing[["detection"]] <- proc.time()[["elapsed"]] - t0

  empty_df <- function() {
    data.frame(
      id = integer(0), text = character(0), score = numeric(0),
      det_score = numeric(0), xmin = numeric(0), ymin = numeric(0),
      xmax = numeric(0), ymax = numeric(0),
      x1 = numeric(0), y1 = numeric(0), x2 = numeric(0), y2 = numeric(0),
      x3 = numeric(0), y3 = numeric(0), x4 = numeric(0), y4 = numeric(0),
      angle = numeric(0), stringsAsFactors = FALSE
    )
  }

  if (is.null(det)) {
    res_df <- empty_df()
    class(res_df) <- c("pliman_ocr", "data.frame")
    attr(res_df, "image") <- img_obj
    attr(res_df, "timing") <- timing
    if (isTRUE(verbose)) cli::cli_alert_info("No text detected in the image.")
    if (isTRUE(plot)) plot(res_df, col = col, lwd = lwd, cex = cex, badge = badge, ...)
    return(invisible(res_df))
  }

  quads <- det$quads
  n_b <- nrow(quads)
  angle <- det$angle

  # ---- Text-line orientation (0 / 180 degrees) -------------------------------
  if (!is.null(ori_file)) {
    t0 <- proc.time()[["elapsed"]]
    flip <- logical(n_b)
    for (s in seq(1L, n_b, by = 32L)) {
      ii <- s:min(n_b, s + 31L)
      tens <- ocr_warp_crops_cpp(
        mat, quads[ii, , drop = FALSE], out_h = 80L, out_w = rep(160L, length(ii)), pad_w = 160L,
        mean = c(0.485, 0.456, 0.406), std = c(0.229, 0.224, 0.225), bgr = TRUE, scale = pix_scale
      )
      pr <- run_ocr_cls_cpp(
        crops_vec = tens, batch_size = length(ii), crop_w = 160L, crop_h = 80L,
        model_path = ori_file, lib_path = lib_path,
        num_threads = as.integer(threads), use_gpu = use_gpu
      )
      if (ncol(pr) >= 2L) flip[ii] <- pr[, 2] > 0.90
    }
    if (any(flip)) {
      quads[flip, ] <- quads[flip, c(5, 6, 7, 8, 1, 2, 3, 4), drop = FALSE]
      angle[flip] <- ((angle[flip] + 360) %% 360) - 180
    }
    timing[["orientation"]] <- proc.time()[["elapsed"]] - t0
  }

  # ---- Text recognition ------------------------------------------------------
  t0 <- proc.time()[["elapsed"]]
  txts <- character(n_b)
  scores <- numeric(n_b)
  wh <- .ocr_quad_wh(quads)
  ratio <- wh$w / wh$h
  rec_batch <- max(1L, as.integer(rec_batch))

  if (is_trocr) {
    for (s in seq(1L, n_b, by = rec_batch)) {
      ii <- s:min(n_b, s + rec_batch - 1L)
      tens <- ocr_warp_crops_cpp(
        mat, quads[ii, , drop = FALSE], out_h = 384L, out_w = rep(384L, length(ii)), pad_w = 384L,
        mean = c(0.5, 0.5, 0.5), std = c(0.5, 0.5, 0.5), bgr = FALSE, scale = pix_scale
      )
      rec_res <- run_trocr_rec_cpp(
        crops_vec = tens,
        batch_size = length(ii),
        encoder_path = trocr_enc_file,
        decoder_path = trocr_dec_file,
        lib_path = lib_path,
        max_length = 32L,
        num_threads = as.integer(threads),
        use_gpu = use_gpu
      )
      for (j in seq_along(ii)) {
        tok_ids <- rec_res[[j]]$token_ids
        if (length(tok_ids) > 0 && length(trocr_vocab) > 0) {
          idx <- tok_ids + 1L
          pieces <- ifelse(idx >= 1L & idx <= length(trocr_vocab),
                           trocr_vocab[pmax(1L, pmin(idx, length(trocr_vocab)))], "")
          txts[ii[j]] <- trimws(gsub("\u2581", " ", paste(pieces, collapse = "")))
          scores[ii[j]] <- round(rec_res[[j]]$score, 4)
        }
      }
    }
  } else {
    imgH <- 48L
    ord <- order(ratio)
    for (s in seq(1L, n_b, by = rec_batch)) {
      ii <- ord[s:min(n_b, s + rec_batch - 1L)]
      max_ratio <- max(320 / imgH, ratio[ii])
      pad_w <- as.integer(min(3200L, ceiling(imgH * max_ratio / 32) * 32))
      cw <- as.integer(pmin(pad_w, pmax(8L, ceiling(imgH * ratio[ii]))))
      tens <- ocr_warp_crops_cpp(
        mat, quads[ii, , drop = FALSE], out_h = imgH, out_w = cw, pad_w = pad_w,
        mean = c(0.5, 0.5, 0.5), std = c(0.5, 0.5, 0.5), bgr = TRUE, scale = pix_scale
      )
      rec_res <- run_ocr_rec_cpp(
        crops_vec = tens,
        batch_size = length(ii),
        crop_w = pad_w,
        crop_h = imgH,
        model_path = rec_file,
        lib_path = lib_path,
        num_threads = as.integer(threads),
        use_gpu = use_gpu,
        allowed_ids = allowed_ids
      )
      for (j in seq_along(ii)) {
        ids <- rec_res[[j]]$char_ids
        if (length(ids) > 0 && length(dict) > 0) {
          chars <- ifelse(ids >= 1L & ids <= length(dict), dict[pmax(1L, pmin(ids, length(dict)))], "")
          txts[ii[j]] <- trimws(paste0(chars, collapse = ""))
          scores[ii[j]] <- round(rec_res[[j]]$score, 4)
        }
      }
    }
  }
  timing[["recognition"]] <- proc.time()[["elapsed"]] - t0

  # Output coordinates are 1-based pixel indices (R / pliman convention).
  qx <- quads[, c(1, 3, 5, 7), drop = FALSE] + 1
  qy <- quads[, c(2, 4, 6, 8), drop = FALSE] + 1
  res_df <- data.frame(
    id = seq_len(n_b),
    text = txts,
    score = scores,
    det_score = round(det$det_score, 4),
    xmin = pmax(1, floor(apply(qx, 1, min))),
    ymin = pmax(1, floor(apply(qy, 1, min))),
    xmax = pmin(orig_w, ceiling(apply(qx, 1, max))),
    ymax = pmin(orig_h, ceiling(apply(qy, 1, max))),
    x1 = qx[, 1], y1 = qy[, 1], x2 = qx[, 2], y2 = qy[, 2],
    x3 = qx[, 3], y3 = qy[, 3], x4 = qx[, 4], y4 = qy[, 4],
    angle = round(angle, 2),
    stringsAsFactors = FALSE
  )
  # Fragments of the probability map (isolated strokes, gaps between words) yield regions
  # with no recognized character; drop them unless requested.
  if (!isTRUE(keep_empty)) {
    res_df <- res_df[nzchar(res_df$text), , drop = FALSE]
    if (nrow(res_df) == 0L) res_df <- empty_df()
    res_df$id <- seq_len(nrow(res_df))
    rownames(res_df) <- NULL
  }

  # Sort in natural reading order (top-to-bottom, line-by-line, left-to-right)
  if (isTRUE(sort) && nrow(res_df) > 1L) {
    y_min_vec <- res_df$ymin
    x_min_vec <- res_df$xmin
    ord <- order(y_min_vec, x_min_vec)
    df_s <- res_df[ord, ]

    lines_info <- list()
    line_assign <- integer(nrow(df_s))

    for (i in seq_len(nrow(df_s))) {
      b_ymin <- df_s$ymin[i]
      b_ymax <- df_s$ymax[i]
      b_h <- b_ymax - b_ymin + 1
      b_cy <- (b_ymin + b_ymax) / 2

      best_line <- 0L
      best_ov <- 0

      if (length(lines_info) > 0) {
        for (l in seq_along(lines_info)) {
          l_ymin <- lines_info[[l]]$ymin
          l_ymax <- lines_info[[l]]$ymax
          l_h <- l_ymax - l_ymin + 1

          ov_h <- max(0, min(b_ymax, l_ymax) - max(b_ymin, l_ymin) + 1)
          ov_ratio <- ov_h / min(b_h, l_h)

          if (ov_ratio > 0.35 && ov_ratio > best_ov) {
            best_ov <- ov_ratio
            best_line <- l
          }
        }
      }

      if (best_line > 0L) {
        line_assign[i] <- best_line
        lines_info[[best_line]]$ymin <- min(lines_info[[best_line]]$ymin, b_ymin)
        lines_info[[best_line]]$ymax <- max(lines_info[[best_line]]$ymax, b_ymax)
        lines_info[[best_line]]$cy <- (lines_info[[best_line]]$ymin + lines_info[[best_line]]$ymax) / 2
      } else {
        new_l <- length(lines_info) + 1L
        line_assign[i] <- new_l
        lines_info[[new_l]] <- list(ymin = b_ymin, ymax = b_ymax, cy = b_cy)
      }
    }

    df_s$line <- line_assign
    line_cys <- vapply(lines_info, function(l) l$cy, numeric(1))
    line_order_map <- order(line_cys)
    df_s$line <- match(df_s$line, line_order_map)

    res_df <- df_s[order(df_s$line, df_s$xmin), ]
    res_df$id <- seq_len(nrow(res_df))
    rownames(res_df) <- NULL
  }

  if (min_score > 0) {
    res_df <- res_df[res_df$score >= min_score | is.na(res_df$score), , drop = FALSE]
    if (nrow(res_df) > 0L) res_df$id <- seq_len(nrow(res_df))
  }

  class(res_df) <- c("pliman_ocr", "data.frame")
  attr(res_df, "image") <- img_obj
  attr(res_df, "timing") <- timing

  if (isTRUE(verbose)) {
    n_rec <- sum(nzchar(res_df$text))
    cli::cli_alert_success("pliman OCR recognized {n_rec} text instance{?s} across {nrow(res_df)} region{?s}.")
  }

  if (isTRUE(plot)) {
    plot(res_df, col = col, lwd = lwd, cex = cex, badge = badge, labels = labels, min_score = min_score, ...)
  }

  invisible(res_df)
}

#' @export
print.pliman_ocr <- function(x, ...) {
  n <- nrow(x)
  if (n == 0L) {
    cli::cli_alert_info("No text detected in the image.")
  } else {
    n_rec <- sum(nzchar(x$text))
    cli::cli_alert_success("pliman OCR recognized {n_rec} text instance{?s} across {n} region{?s}:")
    cols <- intersect(c("id", "line", "text", "score", "xmin", "ymin", "xmax", "ymax"), names(x))
    print.data.frame(x[, cols, drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

#' @export
as.character.pliman_ocr <- function(x, collapse = "\n", ...) {
  if (nrow(x) == 0L) return("")
  if ("line" %in% names(x)) {
    lines_txt <- vapply(sort(unique(x$line)), function(l) {
      sub_t <- x$text[x$line == l & nzchar(x$text)]
      paste(sub_t, collapse = " ")
    }, character(1))
    paste(lines_txt[nzchar(lines_txt)], collapse = collapse)
  } else {
    paste(x$text[nzchar(x$text)], collapse = " ")
  }
}

#' @export
plot.pliman_ocr <- function(x,
                            col = "forestgreen",
                            lwd = 2L,
                            cex = 0.80,
                            badge = TRUE,
                            badge_col = col,
                            text_col = "white",
                            labels = c("text", "id", "both", "none"),
                            min_score = 0.0,
                            ...) {
  if (missing(labels) || is.null(labels) || !is.character(labels) || length(labels) == 0L) {
    labels <- "text"
  } else {
    labels <- labels[1L]
    if (!labels %in% c("text", "id", "both", "none")) labels <- "text"
  }

  img <- attr(x, "image")
  if (!is.null(img)) {
    plot(img, ...)
  }
  if (is.null(x) || nrow(x) == 0L) return(invisible(x))

  if (min_score > 0) {
    x <- x[x$score >= min_score | is.na(x$score), , drop = FALSE]
    if (nrow(x) == 0L) return(invisible(x))
  }

  has_quad <- all(c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4") %in% names(x))
  for (i in seq_len(nrow(x))) {
    bx1 <- x$xmin[i]; by1 <- x$ymin[i]
    bx2 <- x$xmax[i]; by2 <- x$ymax[i]
    if (has_quad && all(is.finite(c(x$x1[i], x$y1[i], x$x3[i], x$y3[i])))) {
      graphics::polygon(
        c(x$x1[i], x$x2[i], x$x3[i], x$x4[i]),
        c(x$y1[i], x$y2[i], x$y3[i], x$y4[i]),
        border = col, lwd = lwd
      )
    } else {
      graphics::rect(bx1, by1, bx2, by2, border = col, lwd = lwd)
    }

    if (labels == "none") next

    lbl <- switch(
      labels,
      "text" = x$text[i],
      "id"   = as.character(x$id[i]),
      "both" = paste0(x$id[i], ": ", x$text[i]),
      x$text[i]
    )

    if (nzchar(lbl)) {
      tw <- graphics::strwidth(lbl, cex = cex)
      th <- graphics::strheight(lbl, cex = cex)
      pad_x <- max(3, tw * 0.12)
      pad_y <- max(2, th * 0.20)

      ty1 <- if (by1 - th - 2 * pad_y >= 0) (by1 - th - 2 * pad_y) else by1
      ty2 <- ty1 + th + 2 * pad_y
      tx1 <- bx1
      tx2 <- bx1 + tw + 2 * pad_x

      if (isTRUE(badge)) {
        graphics::rect(tx1, ty1, tx2, ty2, col = badge_col, border = NA)
        graphics::text(tx1 + pad_x, ty1 + pad_y + th / 2, labels = lbl, col = text_col, cex = cex, adj = c(0, 0.5), font = 2)
      } else {
        graphics::text(tx1 + 1, ty1 + 1, labels = lbl, col = "black", cex = cex, adj = c(0, 0), font = 2)
        graphics::text(tx1, ty1, labels = lbl, col = col, cex = cex, adj = c(0, 0), font = 2)
      }
    }
  }
  invisible(x)
}
