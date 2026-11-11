# ==============================================================================
# ONNX MODEL COMPRESSION & INT8 QUANTIZATION
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Quantize ONNX Models to 8-Bit Precision (INT8)
#' @name onnx_quantize
#' @description
#' `onnx_quantize()` compresses 32-bit floating point (FP32) ONNX neural network models
#' to 8-bit integer (INT8) precision. Dynamic quantization reduces model file sizes
#' by up to 75% (e.g. from 24 MB to 6 MB) and accelerates CPU inference throughput
#' by 2–3x with negligible accuracy trade-off.
#'
#' @param model_path Character string. File path to the input FP32 `.onnx` model.
#' @param output_path Optional file path for the quantized INT8 model. If `NULL` (default),
#'   saves the model in [pliman_model_dir()] appending `"_int8.onnx"` to the original model name.
#' @param weight_type Precision type: `"QUInt8"` (unsigned 8-bit, default) or `"QInt8"` (signed 8-bit).
#' @param verbose Logical. Display compression metrics and file sizes. Defaults to `TRUE`.
#'
#' @return The normalized file path to the quantized `.onnx` model.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Quantize classification model
#' q_model <- onnx_quantize("dataset_flores_yolo_cls.onnx")
#'
#' # Run inference with 3x speedup on CPU
#' res <- image_classify_dl("rose.jpg", model = q_model)
onnx_quantize <- function(model_path,
                          output_path = NULL,
                          weight_type = c("QUInt8", "QInt8"),
                          verbose = TRUE) {
  weight_type <- match.arg(weight_type)

  # Resolve path
  actual_path <- if (file.exists(model_path)) {
    normalizePath(model_path, winslash = "/", mustWork = TRUE)
  } else {
    model_dir <- pliman_model_dir()
    cand <- file.path(model_dir, model_path)
    if (file.exists(cand)) {
      normalizePath(cand, winslash = "/", mustWork = TRUE)
    } else {
      model_path
    }
  }

  if (!file.exists(actual_path)) {
    cli::cli_abort("Model file {.path {model_path}} was not found.")
  }

  # If output_path is not specified, place in pliman_model_dir()
  if (is.null(output_path)) {
    model_dir <- pliman_model_dir()
    clean_stem <- sub("(_int8|_quantized)?\\.onnx$", "", basename(actual_path), ignore.case = TRUE)
    output_path <- file.path(model_dir, paste0(clean_stem, "_int8.onnx"))
  }
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)

  if (isTRUE(verbose)) {
    cli::cli_h2("ONNX Model Quantization (pliman)")
    cli::cli_alert_info("Input model : {.file {actual_path}}")
    cli::cli_alert_info("Output model: {.file {output_path}}")
  }

  # Clean existing destination if present
  if (file.exists(output_path)) unlink(output_path)

  # Python script content
  py_code <- sprintf("
import sys
try:
    from onnxruntime.quantization import quantize_dynamic, QuantType
    w_type = QuantType.%s
    quantize_dynamic(
        model_input=r'%s',
        model_output=r'%s',
        weight_type=w_type
    )
except Exception as e:
    sys.exit(1)
", weight_type, actual_path, output_path)

  success <- FALSE

  # Check dedicated python executable first
  py_exec <- try(.find_python_exec(), silent = TRUE)
  if (!inherits(py_exec, "try-error") && !is.null(py_exec) && file.exists(py_exec)) {
    tmp_py <- tempfile(fileext = ".py")
    writeLines(py_code, tmp_py)
    out <- try(system2(py_exec, args = c(tmp_py), stdout = FALSE, stderr = FALSE), silent = TRUE)
    unlink(tmp_py)
    if (file.exists(output_path) && file.info(output_path)$size > 0) {
      success <- TRUE
    }
  }

  # Fallback to reticulate if available and onnxruntime is present
  if (!success && requireNamespace("reticulate", quietly = TRUE)) {
    if (isTRUE(reticulate::py_module_available("onnxruntime"))) {
      res <- try(reticulate::py_run_string(py_code), silent = TRUE)
      if (!inherits(res, "try-error") && file.exists(output_path) && file.info(output_path)$size > 0) {
        success <- TRUE
      }
    }
  }

  # Fallback to system python/python3
  if (!success) {
    tmp_py <- tempfile(fileext = ".py")
    writeLines(py_code, tmp_py)
    for (cand_py in c("python", "python3")) {
      p_bin <- Sys.which(cand_py)
      if (nzchar(p_bin)) {
        suppressWarnings(system2(p_bin, args = c(shQuote(tmp_py)), stdout = FALSE, stderr = FALSE))
        if (file.exists(output_path) && file.info(output_path)$size > 0) {
          success <- TRUE
          break
        }
      }
    }
    unlink(tmp_py)
  }

  if (!success || !file.exists(output_path)) {
    cli::cli_abort(c(
      "!" = "Failed to quantize model.",
      "i" = "Please ensure {.pkg onnxruntime} is installed in Python: {.code reticulate::py_install('onnxruntime')}"
    ))
  }

  # Also ensure model is cached in pliman's model directory so image_detect_dl() finds it automatically
  model_dir <- pliman_model_dir()
  cache_target <- file.path(model_dir, basename(output_path))
  if (output_path != cache_target && file.exists(output_path)) {
    file.copy(output_path, cache_target, overwrite = TRUE)
  }

  sz_orig <- file.info(actual_path)$size / (1024 * 1024)
  sz_quant <- file.info(output_path)$size / (1024 * 1024)
  ratio <- (1 - (sz_quant / sz_orig)) * 100

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Model quantized successfully!")
    cli::cli_alert_info("Saved to    : {.file {output_path}}")
    cli::cli_text("Original size : {.val {round(sz_orig, 2)}} MB")
    cli::cli_text("Quantized size: {.val {round(sz_quant, 2)}} MB ({cli::col_green(sprintf('-%.1f%%', ratio))})")
  }

  invisible(output_path)
}

#' @title Convert Installed ONNX Models to 8-Bit Precision (INT8)
#' @name model_to_int8
#' @aliases models_to_int8
#' @description
#' `model_to_int8()` converts installed or custom ONNX neural network models to 8-bit integer
#' (INT8) precision. Dynamic quantization reduces model disk size by up to 75% and accelerates
#' CPU inference by 2–3x with virtually no loss in detection or segmentation accuracy.
#'
#' When `model = NULL` (default), `model_to_int8()` scans the pliman model directory
#' ([pliman_model_dir()]) and converts all installed FP32 `.onnx` models into `_int8.onnx` models.
#' Existing INT8 models are skipped unless `overwrite = TRUE`.
#'
#' @param model Optional character vector of model name(s) or file path(s) to quantize.
#'   If `NULL` (default), all unquantized `.onnx` models in `dir` are converted.
#' @param dir Directory where installed models are located. Defaults to [pliman_model_dir()].
#' @param weight_type Precision type: `"QUInt8"` (unsigned 8-bit, default) or `"QInt8"` (signed 8-bit).
#' @param overwrite Logical. If `TRUE`, re-quantizes models even if an `_int8.onnx` version already exists. Defaults to `FALSE`.
#' @param verbose Logical. Show progress messages and compression summary. Defaults to `TRUE`.
#'
#' @return A `data.frame` of class `pliman_quantize_summary` listing each model, original size,
#'   INT8 size, reduction percentage, and conversion status.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Convert a specific installed model
#' model_to_int8("capsulas_det.onnx")
#'
#' # Convert all installed models in pliman_model_dir()
#' summary <- model_to_int8()
#' print(summary)
model_to_int8 <- function(model = NULL,
                          dir = pliman_model_dir(),
                          weight_type = c("QUInt8", "QInt8"),
                          overwrite = FALSE,
                          verbose = TRUE) {
  weight_type <- match.arg(weight_type)
  dir <- pliman_model_dir(dir)

  # If model is NULL, discover all .onnx models in dir
  if (is.null(model)) {
    all_files <- list.files(dir, pattern = "\\.onnx$", full.names = TRUE)
    # Exclude already quantized models
    fp32_files <- all_files[!grepl("(_int8|_quantized)\\.onnx$", all_files, ignore.case = TRUE)]
    if (length(fp32_files) == 0L) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("No FP32 ONNX models found to convert in {.path {dir}}.")
      }
      return(invisible(structure(data.frame(), class = c("pliman_quantize_summary", "data.frame"))))
    }
    model_targets <- fp32_files
  } else {
    model_targets <- vapply(model, function(m) {
      if (file.exists(m)) {
        normalizePath(m, winslash = "/", mustWork = TRUE)
      } else {
        cand <- file.path(dir, m)
        if (file.exists(cand)) {
          normalizePath(cand, winslash = "/", mustWork = TRUE)
        } else if (file.exists(paste0(cand, ".onnx"))) {
          normalizePath(paste0(cand, ".onnx"), winslash = "/", mustWork = TRUE)
        } else {
          m
        }
      }
    }, character(1), USE.NAMES = FALSE)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("ONNX INT8 Batch Model Quantization (pliman)")
    cli::cli_alert_info("Found {length(model_targets)} candidate model{?s} in {.path {dir}}.")
  }

  rows <- list()

  for (i in seq_along(model_targets)) {
    m_path <- model_targets[i]
    m_name <- basename(m_path)
    clean_stem <- sub("\\.onnx$", "", m_name, ignore.case = TRUE)
    out_name <- paste0(clean_stem, "_int8.onnx")
    out_path <- file.path(dir, out_name)

    sz_orig <- if (file.exists(m_path)) file.info(m_path)$size / (1024 * 1024) else NA_real_

    # Check if already exists
    if (!isTRUE(overwrite) && file.exists(out_path) && file.info(out_path)$size > 1000) {
      sz_quant <- file.info(out_path)$size / (1024 * 1024)
      ratio <- if (!is.na(sz_orig) && sz_orig > 0) (1 - (sz_quant / sz_orig)) * 100 else NA_real_
      if (isTRUE(verbose)) {
        cli::cli_alert_info("[{i}/{length(model_targets)}] {.val {m_name}} -> already converted ({round(sz_quant, 2)} MB). Skipping.")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        int8_size_mb = round(sz_quant, 2),
        reduction_pct = round(ratio, 1),
        status = "already_exists",
        path = out_path,
        stringsAsFactors = FALSE
      )
      next
    }

    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "[{i}/{length(model_targets)}] Quantizing {.val {m_name}} to INT8 ({round(sz_orig, 1)} MB)...",
        msg_done = "[{i}/{length(model_targets)}] Quantized {.val {m_name}} successfully"
      )
    }

    res_q <- tryCatch({
      onnx_quantize(
        model_path = m_path,
        output_path = out_path,
        weight_type = weight_type,
        verbose = FALSE
      )
    }, error = function(e) {
      if (isTRUE(verbose)) {
        cli::cli_alert_danger("Failed to quantize {.val {m_name}}: {e$message}")
      }
      NULL
    })

    if (!is.null(res_q) && file.exists(out_path)) {
      sz_quant <- file.info(out_path)$size / (1024 * 1024)
      ratio <- (1 - (sz_quant / sz_orig)) * 100
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        int8_size_mb = round(sz_quant, 2),
        reduction_pct = round(ratio, 1),
        status = "converted",
        path = out_path,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        int8_size_mb = NA_real_,
        reduction_pct = NA_real_,
        status = "failed",
        path = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  class(summary_df) <- c("pliman_quantize_summary", "data.frame")

  if (isTRUE(verbose)) {
    print(summary_df)
  }

  invisible(summary_df)
}

#' @rdname model_to_int8
#' @export
models_to_int8 <- model_to_int8

#' @export
print.pliman_quantize_summary <- function(x, ...) {
  if (nrow(x) == 0L) {
    cat("Empty quantize summary.\n")
    return(invisible(x))
  }
  cli::cli_h2("ONNX Model INT8 Quantization Summary")

  tot_orig <- sum(x$orig_size_mb, na.rm = TRUE)
  tot_int8 <- sum(x$int8_size_mb, na.rm = TRUE)
  tot_ratio <- if (tot_orig > 0) (1 - (tot_int8 / tot_orig)) * 100 else 0
  n_conv <- sum(x$status %in% c("converted", "already_exists"))

  cli::cli_alert_success("{n_conv}/{nrow(x)} models available in INT8 format.")
  cli::cli_text("Total FP32 size: {.val {round(tot_orig, 1)}} MB")
  cli::cli_text("Total INT8 size: {.val {round(tot_int8, 1)}} MB ({cli::col_green(sprintf('-%.1f%% disk space', tot_ratio))})")

  print.data.frame(x[, c("model", "orig_size_mb", "int8_size_mb", "reduction_pct", "status")], row.names = FALSE)
  invisible(x)
}


# ==============================================================================
# ONNX FP16 HALF-PRECISION MODEL EXPORT & CONVERSION
# ==============================================================================

#' @title Export YOLO Model to FP16 Half-Precision ONNX
#' @name yolo_export_fp16
#' @description
#' `yolo_export_fp16()` converts a YOLO model (`.pt` PyTorch weights or `.onnx` model)
#' to 16-bit half precision (`FLOAT16`) ONNX format.
#'
#' The exported model is saved by default in [pliman_model_dir()], so it can be used
#' right away by name in any pliman inference function (e.g.
#' `image_detect_dl(img, model = "yolo26x_fp16.onnx")`).
#'
#' @section What FP16 always gives you:
#' \itemize{
#'   \item \bold{~50\% smaller files and ~50\% less VRAM} for the weights (e.g. YOLO26x:
#'     212.9 MB to 106.5 MB; YOLO11n: 10.2 MB to 5.2 MB). This lets larger models, or
#'     several models at once, fit in GPUs with limited memory (e.g. 4 GB cards).
#'   \item \bold{Practically identical predictions.} Boxes, classes and scores match the
#'     FP32 model up to small rounding differences (e.g. a score of 0.9666 in FP32 vs
#'     0.9658 in FP16).
#'   \item \bold{Portability.} The `.onnx` file runs on any machine supported by pliman
#'     and can be shared (e.g. on Hugging Face), unlike TensorRT `.engine` files, which
#'     are tied to one GPU/driver/CUDA combination (see [yolo_export_engine()]).
#' }
#'
#' @section Speed depends on the hardware:
#' FP16 is \bold{not} guaranteed to be faster. The speed-up comes from dedicated
#' half-precision units, so the result depends on the GPU:
#' \itemize{
#'   \item \bold{GPUs with Tensor Cores} (NVIDIA RTX 20xx/30xx/40xx/50xx, Quadro/RTX A-series,
#'     data-center GPUs such as T4/A100): FP16 is expected to be faster than FP32,
#'     mainly for large models (`m`, `l`, `x` variants) and high resolutions.
#'   \item \bold{GPUs without Tensor Cores} (e.g. NVIDIA GTX 10xx/16xx, many integrated
#'     GPUs): FP16 is usually \bold{equal or slower} than FP32. Example measured with
#'     YOLO26x on a GTX 1650 via DirectML (`bench::mark()`, 30 iterations): FP32 median
#'     102 ms vs FP16 median 127 ms (~25\% slower).
#'   \item \bold{CPU} (`use_gpu = FALSE`): most CPUs have no fast FP16 arithmetic, so FP16
#'     models are typically slower than FP32 on CPU. For CPU speed-ups use INT8
#'     quantization instead ([onnx_quantize()], [model_to_int8()]).
#' }
#' Small models (`n`, `s`) are often limited by pre/post-processing rather than by the
#' network itself, so gains from FP16 are small even on Tensor Core GPUs.
#'
#' \bold{Recommendation:} keep both versions and benchmark them on your own machine
#' before switching (see Examples). Use FP16 when it is faster on your GPU or when you
#' need to save VRAM/disk; otherwise keep FP32.
#'
#' @section Conversion details:
#' \itemize{
#'   \item `.pt` checkpoints are exported with Ultralytics (`format = "onnx", half = TRUE`).
#'   \item `.onnx` models are converted with \pkg{onnxconverter_common}. Models with opset
#'     < 17 are upgraded to opset 17, internal `Cast` nodes are aligned to FLOAT16,
#'     stale shape information is removed, and the original metadata (classes, stride,
#'     image size, task) is preserved so [yolo_model_details()] keeps working.
#'   \item Input and output tensors of the FP16 model are FLOAT16. pliman converts the
#'     image to half precision and the predictions back to single precision
#'     automatically, so no change is needed in your code.
#'   \item Required Python packages: `pip install onnx onnxconverter_common`
#'     (plus `ultralytics` for `.pt` files).
#' }
#'
#' @param model Character string. File path to a `.pt` PyTorch weights file or `.onnx` model,
#'   or the name of an installed model in [pliman_model_dir()].
#' @param output Optional destination file path for the exported `.onnx` model. If `NULL` (default),
#'   saves the model in [pliman_model_dir()] appending `"_fp16.onnx"` to the model name.
#' @param imgsz Integer or vector of two integers. Inference image resolution `c(height, width)`
#'   or single integer for square inputs. Defaults to `640`.
#' @param dynamic Logical. Enable dynamic input axes. Defaults to `FALSE`.
#' @param simplify Logical. Simplify ONNX graph using `onnxslim`. Defaults to `TRUE`.
#' @param opset Integer. ONNX opset version. Defaults to `17L`.
#' @param verbose Logical. Display conversion metrics and file size reduction. Defaults to `TRUE`.
#'
#' @return The normalized file path to the exported FP16 `.onnx` model.
#' @seealso [model_to_fp16()] for batch conversion, [onnx_quantize()] / [model_to_int8()]
#'   for INT8 (CPU speed-up), [yolo_export_engine()] for TensorRT engines.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Convert a YOLO ONNX model (or .pt checkpoint) to FP16.
#' #    The file is saved in pliman_model_dir() as "yolo26x_fp16.onnx".
#' fp16_model <- yolo_export_fp16("yolo26x.onnx")
#'
#' # 2. Use it by name, exactly like the FP32 model
#' img <- image_pliman("apartment.jpeg")
#' det <- image_detect_dl(img, model = "yolo26x_fp16.onnx")
#'
#' # 3. Benchmark FP16 vs FP32 on YOUR hardware before switching
#' #    (FP16 is only expected to be faster on GPUs with Tensor Cores)
#' bench::mark(
#'   fp16 = image_detect_dl(img, model = "yolo26x_fp16.onnx", plot = FALSE),
#'   fp32 = image_detect_dl(img, model = "yolo26x.onnx", plot = FALSE),
#'   iterations = 30,
#'   check = FALSE
#' )
#'
#' # 4. Compare model metadata and tensor types
#' yolo_model_details("yolo26x_fp16.onnx")
yolo_export_fp16 <- function(model,
                             output = NULL,
                             imgsz = 640,
                             dynamic = FALSE,
                             simplify = TRUE,
                             opset = 17L,
                             verbose = TRUE) {
  # Resolve model path
  actual_path <- if (file.exists(model)) {
    normalizePath(model, winslash = "/", mustWork = TRUE)
  } else {
    model_dir <- pliman_model_dir()
    cand <- file.path(model_dir, model)
    if (file.exists(cand)) {
      normalizePath(cand, winslash = "/", mustWork = TRUE)
    } else if (file.exists(paste0(cand, ".onnx"))) {
      normalizePath(paste0(cand, ".onnx"), winslash = "/", mustWork = TRUE)
    } else if (file.exists(paste0(cand, ".pt"))) {
      normalizePath(paste0(cand, ".pt"), winslash = "/", mustWork = TRUE)
    } else {
      model
    }
  }

  if (!file.exists(actual_path)) {
    cli::cli_abort("Model file {.path {model}} was not found.")
  }

  is_pt <- grepl("\\.pt$", actual_path, ignore.case = TRUE)
  is_onnx <- grepl("\\.onnx$", actual_path, ignore.case = TRUE)

  if (!is_pt && !is_onnx) {
    cli::cli_abort("Model must be a {.file .pt} or {.file .onnx} file.")
  }

  # If output is not specified, place in pliman_model_dir()
  if (is.null(output)) {
    stem <- tools::file_path_sans_ext(basename(actual_path))
    # Avoid duplicate suffixes
    clean_stem <- sub("(_fp16|_half)$", "", stem, ignore.case = TRUE)
    output <- file.path(pliman_model_dir(), paste0(clean_stem, "_fp16.onnx"))
  }
  output_path <- normalizePath(output, winslash = "/", mustWork = FALSE)

  if (isTRUE(verbose)) {
    cli::cli_h2("Exporting Model to FP16 Half-Precision (pliman)")
    cli::cli_alert_info("Input model : {.file {actual_path}}")
    cli::cli_alert_info("Output model: {.file {output_path}}")
  }

  # Clean existing destination if present
  if (file.exists(output_path)) unlink(output_path)

  py_exec <- try(.find_python_exec(), silent = TRUE)
  if (inherits(py_exec, "try-error") || is.null(py_exec) || !file.exists(py_exec)) {
    cli::cli_abort("A valid Python installation was not found.")
  }

  py_dyn_str <- if (isTRUE(dynamic)) "True" else "False"
  py_simp_str <- if (isTRUE(simplify)) "True" else "False"

  if (is_pt) {
    # PyTorch export via Ultralytics
    py_code <- sprintf("
import sys, os, shutil
from ultralytics import YOLO

try:
    m = YOLO(r'%s')
    out = m.export(format='onnx', imgsz=%d, dynamic=%s, simplify=%s, opset=%d, half=True)
    if out and os.path.exists(out):
        dest = r'%s'
        if os.path.abspath(out) != os.path.abspath(dest):
            shutil.copy2(out, dest)
            try:
                os.remove(out)
            except Exception:
                pass
        sys.exit(0)
    else:
        sys.exit(1)
except Exception as e:
    print('Export error:', e, file=sys.stderr)
    sys.exit(1)
", actual_path, as.integer(imgsz[1]), py_dyn_str, py_simp_str, as.integer(opset), output_path)
  } else {
    # Direct ONNX to FP16 conversion via onnxconverter_common + graph validator
    py_code <- sprintf("
import sys, os
try:
    import onnx
    from onnx import version_converter
    from onnxconverter_common import float16

    m_orig = onnx.load(r'%s')

    # 1. Upgrade opset to at least 17 if lower, to avoid opset-12 Resize/Concat type bugs
    current_opset = 0
    for op in m_orig.opset_import:
        if op.domain == '' or op.domain == 'ai.onnx':
            current_opset = op.version
            break

    if current_opset > 0 and current_opset < 17:
        try:
            m = version_converter.convert_version(m_orig, 17)
        except Exception:
            m = m_orig
    else:
        m = m_orig

    # 2. Convert graph to float16
    try:
        m_fp16 = float16.convert_float_to_float16(m, keep_io_types=False)
    except ValueError as ve:
        if 'already converted to float16' in str(ve).lower():
            m_fp16 = float16.convert_float_to_float16(m, keep_io_types=False, check_fp16_ready=False)
        else:
            raise ve

    # 3. Fix any internal Cast nodes that still cast to FLOAT (1)
    for n in m_fp16.graph.node:
        if n.op_type == 'Cast':
            for attr in n.attribute:
                if attr.name == 'to' and attr.i == 1:
                    attr.i = 10

    # 4. Clear stale/duplicate value_info that causes ONNXRuntime schema conflicts
    del m_fp16.graph.value_info[:]

    # 5. Restore full metadata_props from original model (classes, stride, imgsz, etc.)
    del m_fp16.metadata_props[:]
    for p in m_orig.metadata_props:
        new_p = m_fp16.metadata_props.add()
        new_p.key = p.key
        new_p.value = p.value

    onnx.save(m_fp16, r'%s')
    sys.exit(0)
except Exception as e:
    print('ONNX FP16 conversion error:', e, file=sys.stderr)
    sys.exit(1)
", actual_path, output_path)
  }

  tmp_py <- tempfile(fileext = ".py")
  writeLines(py_code, tmp_py)
  status <- system2(py_exec, args = shQuote(tmp_py), stdout = FALSE, stderr = FALSE)
  unlink(tmp_py)

  if (!file.exists(output_path) || file.info(output_path)$size == 0) {
    cli::cli_abort(c(
      "!" = "Failed to export model to FP16.",
      "i" = "Ensure required Python packages are installed: {.code pip install ultralytics onnx onnxconverter_common onnxslim}"
    ))
  }

  # Also cache in pliman model directory if not already there
  model_dir <- pliman_model_dir()
  cache_target <- file.path(model_dir, basename(output_path))
  if (output_path != cache_target && file.exists(output_path)) {
    file.copy(output_path, cache_target, overwrite = TRUE)
  }

  sz_orig <- file.info(actual_path)$size / (1024 * 1024)
  sz_fp16 <- file.info(output_path)$size / (1024 * 1024)
  ratio <- if (sz_orig > 0) (1 - (sz_fp16 / sz_orig)) * 100 else 0

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Model exported to FP16 successfully!")
    cli::cli_alert_info("Saved to    : {.file {output_path}}")
    cli::cli_text("Original size : {.val {round(sz_orig, 2)}} MB")
    cli::cli_text("FP16 size     : {.val {round(sz_fp16, 2)}} MB ({cli::col_green(sprintf('-%.1f%%', ratio))})")
    cli::cli_alert_info("Available in pliman as: {.val {basename(output_path)}}")
  }

  invisible(output_path)
}

#' @title Convert Installed Models to 16-Bit Half Precision (FP16)
#' @name model_to_fp16
#' @aliases models_to_fp16
#' @description
#' `model_to_fp16()` converts installed or custom models to 16-bit half precision (FP16).
#' This reduces model disk footprint and GPU VRAM consumption by ~50\% with practically
#' identical predictions.
#'
#' \bold{Speed depends on the hardware.} FP16 is expected to be faster only on GPUs with
#' Tensor Cores (NVIDIA RTX and data-center GPUs). On GPUs without Tensor Cores (e.g. GTX
#' 10xx/16xx) and on CPU, FP16 models are usually equal or slower than FP32. Converted
#' models are saved alongside the originals (FP32 files are never removed), so you can
#' benchmark both and choose. See [yolo_export_fp16()] for details and a benchmark example.
#'
#' When `model = NULL` (default), `model_to_fp16()` scans [pliman_model_dir()] for `.pt`
#' and unquantized `.onnx` models and converts them to `_fp16.onnx`.
#' Existing FP16 models are skipped unless `overwrite = TRUE`.
#'
#' @param model Optional character vector of model name(s) or file path(s) to convert.
#'   If `NULL` (default), all unconverted models in `dir` are processed.
#' @param dir Directory where installed models are located. Defaults to [pliman_model_dir()].
#' @param overwrite Logical. If `TRUE`, re-converts models even if an `_fp16.onnx` version already exists. Defaults to `FALSE`.
#' @param verbose Logical. Show progress messages and compression summary. Defaults to `TRUE`.
#'
#' @return A `data.frame` of class `pliman_fp16_summary` listing each model, original size,
#'   FP16 size, reduction percentage, and conversion status.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Convert all installed models to FP16
#' summary <- model_to_fp16()
#' print(summary)
model_to_fp16 <- function(model = NULL,
                          dir = pliman_model_dir(),
                          overwrite = FALSE,
                          verbose = TRUE) {
  dir <- pliman_model_dir(dir)

  if (is.null(model)) {
    all_files <- list.files(dir, pattern = "\\.(onnx|pt)$", full.names = TRUE)
    # Exclude already converted or quantized models
    fp32_files <- all_files[!grepl("(_fp16|_half|_int8|_quantized)\\.onnx$", all_files, ignore.case = TRUE)]
    if (length(fp32_files) == 0L) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("No candidate models found to convert in {.path {dir}}.")
      }
      return(invisible(structure(data.frame(), class = c("pliman_fp16_summary", "data.frame"))))
    }
    model_targets <- fp32_files
  } else {
    model_targets <- vapply(model, function(m) {
      if (file.exists(m)) {
        normalizePath(m, winslash = "/", mustWork = TRUE)
      } else {
        cand <- file.path(dir, m)
        if (file.exists(cand)) {
          normalizePath(cand, winslash = "/", mustWork = TRUE)
        } else if (file.exists(paste0(cand, ".onnx"))) {
          normalizePath(paste0(cand, ".onnx"), winslash = "/", mustWork = TRUE)
        } else if (file.exists(paste0(cand, ".pt"))) {
          normalizePath(paste0(cand, ".pt"), winslash = "/", mustWork = TRUE)
        } else {
          m
        }
      }
    }, character(1), USE.NAMES = FALSE)
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("ONNX FP16 Batch Model Conversion (pliman)")
    cli::cli_alert_info("Found {length(model_targets)} candidate model{?s} in {.path {dir}}.")
  }

  rows <- list()

  for (i in seq_along(model_targets)) {
    m_path <- model_targets[i]
    m_name <- basename(m_path)
    clean_stem <- sub("\\.(onnx|pt)$", "", m_name, ignore.case = TRUE)
    clean_stem <- sub("(_fp16|_half)$", "", clean_stem, ignore.case = TRUE)
    out_name <- paste0(clean_stem, "_fp16.onnx")
    out_path <- file.path(dir, out_name)

    sz_orig <- if (file.exists(m_path)) file.info(m_path)$size / (1024 * 1024) else NA_real_

    # Check if already exists
    if (!isTRUE(overwrite) && file.exists(out_path) && file.info(out_path)$size > 1000) {
      sz_fp16 <- file.info(out_path)$size / (1024 * 1024)
      ratio <- if (!is.na(sz_orig) && sz_orig > 0) (1 - (sz_fp16 / sz_orig)) * 100 else NA_real_
      if (isTRUE(verbose)) {
        cli::cli_alert_info("[{i}/{length(model_targets)}] {.val {m_name}} -> already converted ({round(sz_fp16, 2)} MB). Skipping.")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        fp16_size_mb = round(sz_fp16, 2),
        reduction_pct = round(ratio, 1),
        status = "already_exists",
        path = out_path,
        stringsAsFactors = FALSE
      )
      next
    }

    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "[{i}/{length(model_targets)}] Converting {.val {m_name}} to FP16 ({round(sz_orig, 1)} MB)...",
        msg_done = "[{i}/{length(model_targets)}] Converted {.val {m_name}} successfully"
      )
    }

    res_q <- tryCatch({
      yolo_export_fp16(
        model = m_path,
        output = out_path,
        verbose = FALSE
      )
    }, error = function(e) {
      if (isTRUE(verbose)) {
        cli::cli_alert_danger("Failed to convert {.val {m_name}}: {e$message}")
      }
      NULL
    })

    if (!is.null(res_q) && file.exists(out_path)) {
      sz_fp16 <- file.info(out_path)$size / (1024 * 1024)
      ratio <- (1 - (sz_fp16 / sz_orig)) * 100
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        fp16_size_mb = round(sz_fp16, 2),
        reduction_pct = round(ratio, 1),
        status = "converted",
        path = out_path,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        model = clean_stem,
        orig_size_mb = round(sz_orig, 2),
        fp16_size_mb = NA_real_,
        reduction_pct = NA_real_,
        status = "failed",
        path = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  class(summary_df) <- c("pliman_fp16_summary", "data.frame")

  if (isTRUE(verbose)) {
    print(summary_df)
  }

  invisible(summary_df)
}

#' @rdname model_to_fp16
#' @export
models_to_fp16 <- model_to_fp16

#' @export
print.pliman_fp16_summary <- function(x, ...) {
  if (nrow(x) == 0L) {
    cat("Empty FP16 summary.\n")
    return(invisible(x))
  }
  cli::cli_h2("ONNX Model FP16 Conversion Summary")

  tot_orig <- sum(x$orig_size_mb, na.rm = TRUE)
  tot_fp16 <- sum(x$fp16_size_mb, na.rm = TRUE)
  tot_ratio <- if (tot_orig > 0) (1 - (tot_fp16 / tot_orig)) * 100 else 0
  n_conv <- sum(x$status %in% c("converted", "already_exists"))

  cli::cli_alert_success("{n_conv}/{nrow(x)} models available in FP16 format.")
  cli::cli_text("Total original size: {.val {round(tot_orig, 1)}} MB")
  cli::cli_text("Total FP16 size    : {.val {round(tot_fp16, 1)}} MB ({cli::col_green(sprintf('-%.1f%% disk space / VRAM', tot_ratio))})")

  print.data.frame(x[, c("model", "orig_size_mb", "fp16_size_mb", "reduction_pct", "status")], row.names = FALSE)
  invisible(x)
}

# ==============================================================================
# NVIDIA TENSORRT ENGINE EXPORT (.engine)
# ==============================================================================

#' @title Export YOLO Model to NVIDIA TensorRT Engine (.engine)
#' @name yolo_export_engine
#' @description
#' `yolo_export_engine()` compiles a YOLO model (`.pt` or `.onnx`) into an NVIDIA TensorRT
#' optimized `.engine` binary using half-precision (FP16) or single-precision (FP32).
#'
#' \bold{Important Notes on TensorRT Engines:}
#' \itemize{
#'   \item TensorRT engines are machine-specific binaries tied to the exact GPU hardware,
#'     CUDA version, and TensorRT runtime on which they were compiled. They cannot be shared across
#'     different machines or uploaded to Hugging Face.
#'   \item Compiling a TensorRT engine requires an NVIDIA GPU, NVIDIA drivers, CUDA, and the
#'     \pkg{tensorrt} Python package (\code{pip install tensorrt}).
#'   \item For a portable, half-size model that runs on any machine without TensorRT, use
#'     [yolo_export_fp16()], which runs natively in pliman via DirectML.
#' }
#'
#' @param model Character string. File path to a `.pt` weights file or `.onnx` model,
#'   or the name of an installed model in [pliman_model_dir()].
#' @param output Optional destination path for the exported `.engine` file. If `NULL` (default),
#'   places the `.engine` in [pliman_model_dir()].
#' @param half Logical. Enable FP16 half-precision in TensorRT engine. Defaults to `TRUE`.
#' @param imgsz Integer or vector of two integers. Input image resolution. Defaults to `640`.
#' @param dynamic Logical. Enable dynamic input axes. Defaults to `FALSE`.
#' @param device Integer. GPU device ordinal (e.g. `0`). Defaults to `0`.
#' @param workspace Integer. TensorRT builder workspace memory in GB. Defaults to `4`.
#' @param verbose Logical. Show compilation progress. Defaults to `TRUE`.
#'
#' @return The file path to the compiled `.engine` file.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Compile YOLO to TensorRT engine (requires NVIDIA GPU & tensorrt package)
#' engine <- yolo_export_engine("yolo26n.pt", half = TRUE)
yolo_export_engine <- function(model,
                               output = NULL,
                               half = TRUE,
                               imgsz = 640,
                               dynamic = FALSE,
                               device = 0,
                               workspace = 4,
                               verbose = TRUE) {
  actual_path <- if (file.exists(model)) {
    normalizePath(model, winslash = "/", mustWork = TRUE)
  } else {
    model_dir <- pliman_model_dir()
    cand <- file.path(model_dir, model)
    if (file.exists(cand)) {
      normalizePath(cand, winslash = "/", mustWork = TRUE)
    } else if (file.exists(paste0(cand, ".pt"))) {
      normalizePath(paste0(cand, ".pt"), winslash = "/", mustWork = TRUE)
    } else if (file.exists(paste0(cand, ".onnx"))) {
      normalizePath(paste0(cand, ".onnx"), winslash = "/", mustWork = TRUE)
    } else {
      model
    }
  }

  if (!file.exists(actual_path)) {
    cli::cli_abort("Model file {.path {model}} was not found.")
  }

  py_exec <- try(.find_python_exec(), silent = TRUE)
  if (inherits(py_exec, "try-error") || is.null(py_exec) || !file.exists(py_exec)) {
    cli::cli_abort("A valid Python installation was not found.")
  }

  # Check if tensorrt module is available
  check_code <- "import sys; import tensorrt; sys.exit(0)"
  tmp_check <- tempfile(fileext = ".py")
  writeLines(check_code, tmp_check)
  trt_check <- suppressWarnings(system2(py_exec, args = shQuote(tmp_check), stdout = FALSE, stderr = FALSE))
  unlink(tmp_check)

  if (trt_check != 0L) {
    cli::cli_abort(c(
      "x" = "NVIDIA TensorRT is not installed in Python.",
      "i" = "TensorRT engines require an NVIDIA GPU, CUDA, and the Python {.pkg tensorrt} package:",
      " " = "  {.code pip install tensorrt}",
      ">" = "Tip: For a portable half-size model without TensorRT, use {.code yolo_export_fp16()}, which runs natively via DirectML in pliman."
    ))
  }

  if (isTRUE(verbose)) {
    cli::cli_h2("Exporting YOLO Model to NVIDIA TensorRT Engine (pliman)")
    cli::cli_alert_info("Input model : {.file {actual_path}}")
    cli::cli_alert_info("Precision   : {.val {if (isTRUE(half)) 'FP16' else 'FP32'}}")
    cli::cli_alert_info("Device      : {.val {device}}")
    cli::cli_alert_info("Workspace   : {.val {workspace}} GB")
  }

  py_dyn_str <- if (isTRUE(dynamic)) "True" else "False"
  py_half_str <- if (isTRUE(half)) "True" else "False"

  target_engine <- if (!is.null(output)) {
    normalizePath(output, winslash = "/", mustWork = FALSE)
  } else {
    file.path(pliman_model_dir(), paste0(tools::file_path_sans_ext(basename(actual_path)), ".engine"))
  }

  py_code <- sprintf("
import sys, os, shutil
from ultralytics import YOLO

try:
    m = YOLO(r'%s')
    out = m.export(
        format='engine',
        half=%s,
        dynamic=%s,
        imgsz=%d,
        device=%d,
        workspace=%d
    )
    if out and os.path.exists(out):
        dest = r'%s'
        if os.path.abspath(out) != os.path.abspath(dest):
            shutil.copy2(out, dest)
        sys.exit(0)
    else:
        sys.exit(1)
except Exception as e:
    print('TensorRT export error:', e, file=sys.stderr)
    sys.exit(1)
", actual_path, py_half_str, py_dyn_str, as.integer(imgsz[1]), as.integer(device), as.integer(workspace), target_engine)

  tmp_py <- tempfile(fileext = ".py")
  writeLines(py_code, tmp_py)
  status <- system2(py_exec, args = shQuote(tmp_py))
  unlink(tmp_py)

  if (!file.exists(target_engine) || file.info(target_engine)$size == 0) {
    cli::cli_abort("TensorRT engine compilation failed.")
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("TensorRT engine compiled successfully!")
    cli::cli_alert_info("Engine file: {.file {target_engine}}")
  }

  invisible(target_engine)
}
