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
#'   appends `"_int8.onnx"` to the original model name.
#' @param weight_type Precision type: `"QUInt8"` (unsigned 8-bit, default) or `"QInt8"` (signed 8-bit).
#' @param verbose Logical. Display compression metrics and file sizes. Defaults to `TRUE`.
#'
#' @return The normalized file path to the quantized `.onnx` model.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Quantize classification model
#' q_model <- onnx_quantize("dataset_flores_yolo_cls.onnx")
#'
#' # Run inference with 3x speedup on CPU
#' res <- image_classify_dl("rose.jpg", model = q_model)
#' }
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

  # If output_path is not specified, place in working directory (or same directory if full path passed)
  if (is.null(output_path)) {
    if (dirname(model_path) == "." || dirname(model_path) == "") {
      output_path <- file.path(getwd(), sub("\\.onnx$", "_int8.onnx", basename(actual_path)))
    } else {
      output_path <- sub("\\.onnx$", "_int8.onnx", actual_path)
    }
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
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Convert a specific installed model
#' model_to_int8("capsulas_det.onnx")
#'
#' # Convert all installed models in pliman_model_dir()
#' summary <- model_to_int8()
#' print(summary)
#' }
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
