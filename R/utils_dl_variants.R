# ==============================================================================
# PRECISION VARIANTS OF THE SAM 2.1 ENCODER AND GROUNDING DINO
# pliman: Plant Image Analysis in R
#
# Builds the alternative model files that .resolve_sam_encoder() and
# .resolve_dino_model() (utils_dl_providers.R) can select. Conversion runs in Python
# (onnx, onnxruntime, onnxconverter_common), like yolo_export_fp16() / onnx_quantize().
# ==============================================================================

.variant_py_sam_int8 <- "
import sys
from onnxruntime.quantization import quantize_dynamic, QuantType
quantize_dynamic(model_input=sys.argv[1], model_output=sys.argv[2], weight_type=QuantType.QUInt8)
"

# FP16 weights with FP32 graph inputs/outputs (keep_io_types), so the C++ runners keep
# feeding and reading float32. onnxconverter_common leaves the model's own Cast(to=FLOAT)
# nodes untouched, which then disagree with their FP16 consumers; only those original
# casts are retargeted (the converter's boundary casts must stay as they are), and the
# stale value_info is dropped.
.variant_py_fp16 <- "
import sys, onnx
import onnxruntime as ort
from onnxconverter_common import float16
src, dst = sys.argv[1], sys.argv[2]
m = onnx.load(src)
orig_casts = {n.name for n in m.graph.node if n.op_type == 'Cast'}
m16 = float16.convert_float_to_float16(m, keep_io_types=True)
for n in m16.graph.node:
    if n.op_type == 'Cast' and n.name in orig_casts:
        for a in n.attribute:
            if a.name == 'to' and a.i == onnx.TensorProto.FLOAT:
                a.i = onnx.TensorProto.FLOAT16
del m16.graph.value_info[:]
onnx.save(m16, dst)
ort.InferenceSession(dst, providers=['CPUExecutionProvider'])
"

# The public FP16 Grounding DINO export runs GridSample in half precision, which has no
# CUDA kernel in ONNX Runtime 1.20: all 48 GridSample nodes fell back to the CPU (69% of
# the time, 7.3 s per call on a GTX 1650 vs 3.2 s for the INT8 model). Wrapping them in
# FP32 casts keeps every node on the GPU (1.75 s per call).
.variant_py_gridsample_fp32 <- "
import sys, onnx
import onnxruntime as ort
from onnx import helper, TensorProto
src, dst = sys.argv[1], sys.argv[2]
m = onnx.load(src)
g = m.graph
nodes = []
n_wrapped = 0
for node in g.node:
    if node.op_type != 'GridSample':
        nodes.append(node)
        continue
    n_wrapped += 1
    tag = node.name or ('gridsample_%d' % n_wrapped)
    for k, name in enumerate(list(node.input)):
        if not name:
            continue
        cast_out = '%s__f32_%s_%d' % (name, tag, k)
        nodes.append(helper.make_node('Cast', [name], [cast_out], to=TensorProto.FLOAT, name=tag + '_in%d_f32' % k))
        node.input[k] = cast_out
    outs = list(node.output)
    for k, name in enumerate(outs):
        node.output[k] = name + '__f32'
    nodes.append(node)
    for k, name in enumerate(outs):
        nodes.append(helper.make_node('Cast', [name + '__f32'], [name], to=TensorProto.FLOAT16, name=tag + '_out%d_f16' % k))
del g.node[:]
g.node.extend(nodes)
del g.value_info[:]
onnx.save(m, dst)
ort.InferenceSession(dst, providers=['CPUExecutionProvider'])
print(n_wrapped)
"

.run_variant_python <- function(code, args, verbose = TRUE) {
  py <- try(.find_python_exec(), silent = TRUE)
  if (inherits(py, "try-error") || is.null(py) || !file.exists(py)) {
    cli::cli_abort(c(
      "x" = "A Python installation was not found.",
      "i" = "Install Python with {.code pip install onnx onnxruntime onnxconverter_common}."
    ))
  }
  tmp <- tempfile(fileext = ".py")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(code, tmp)
  out <- suppressWarnings(system2(py, c(shQuote(tmp), shQuote(args)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(c(
      "x" = "Model conversion failed.",
      "i" = "Required Python packages: {.code pip install onnx onnxruntime onnxconverter_common}",
      " " = paste(utils::tail(out, 5), collapse = "\n")
    ))
  }
  invisible(out)
}

#' Build precision variants of the SAM 2.1 encoder and Grounding DINO
#'
#' @description
#' Creates alternative versions of the large foundation models used by
#' [image_segment_dl()] (`"sam2.1"`, `"persam"`, `"grounded-sam"`) and stores them next to
#' the originals in [pliman_model_dir()]. The originals are never modified.
#'
#' * `"sam2.1-encoder-int8"` -> `sam2.1.encoder_int8.onnx`: dynamic INT8 weights (28 MB
#'   instead of 109 MB). Select it with `options(pliman.sam_encoder = "int8")`. It is
#'   used only by CPU sessions.
#' * `"sam2.1-encoder-fp16"` -> `sam2.1.encoder_fp16.onnx`: FP16 weights with FP32
#'   inputs/outputs (55 MB). Select it with `options(pliman.sam_encoder = "fp16")`.
#' * `"groundingdino-fp16"` -> `groundingdino-tiny_fp16.onnx`: FP16 Grounding DINO, built
#'   from the public FP16 export of the same model (onnx-community/grounding-dino-tiny-ONNX,
#'   about 360 MB download) with its `GridSample` operators kept in FP32. It is picked
#'   automatically with the CUDA and TensorRT providers (see [pliman_gpu_provider()]), or
#'   explicitly with `options(pliman.dino_model = "groundingdino-tiny_fp16.onnx")`.
#'
#' @section Measured effect (hardware dependent):
#' On an Intel i7 (12 threads) with a GTX 1650, where FP32 was the fastest SAM encoder
#' everywhere:
#'
#' | | FP32 | INT8 | FP16 |
#' |---|---|---|---|
#' | SAM encoder, CPU | 2.31 s | 2.31 s (mask IoU 0.995) | 5.75 s |
#' | SAM encoder, DirectML | 0.36 s | 38 s, **wrong masks** | 36.7 s |
#' | SAM encoder, CUDA | 0.39 s | 4.65 s | 1.52 s |
#' | Grounding DINO, CUDA | INT8: 3.2 s | | **1.75 s** |
#' | Grounding DINO, DirectML | INT8: 0.8 s | | 14.5 s |
#'
#' INT8 accelerates only on CPUs with VNNI/AMX instructions (recent Intel/AMD), so it is
#' never used on the GPU (under DirectML it returned wrong masks). FP16 needs GPUs with
#' native half-precision throughput (NVIDIA Tensor Cores: RTX 20xx and newer, T4,
#' A-series); without them it can be many times slower because of casting. Measure on
#' your own machine before switching (for example timing [image_segment_dl()] after
#' `onnx_clear_encoder_cache_cpp()`).
#'
#' @param variant One or more of `"sam2.1-encoder-int8"`, `"sam2.1-encoder-fp16"`,
#'   `"groundingdino-fp16"`.
#' @param dir Model directory. Defaults to [pliman_model_dir()].
#' @param overwrite Logical. Rebuild variants that already exist. Defaults to `FALSE`.
#' @param dino_fp16_url Source of the FP16 Grounding DINO export.
#' @param verbose Logical. Show progress messages.
#' @return A data.frame with the variant, its file path, size in MB and status,
#'   invisibly.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' pliman_build_model_variant("sam2.1-encoder-int8")
#' options(pliman.sam_encoder = "int8")
#' img <- image_pliman("soybean_touch.jpg")
#' res <- image_segment_dl(img, model = "sam2.1", engine = "cpu", plot = FALSE)
#' options(pliman.sam_encoder = "fp32")
pliman_build_model_variant <- function(variant = c("sam2.1-encoder-int8", "sam2.1-encoder-fp16", "groundingdino-fp16"),
                                       dir = pliman_model_dir(),
                                       overwrite = FALSE,
                                       dino_fp16_url = NULL,
                                       verbose = TRUE) {
  if (is.null(dino_fp16_url)) {
    dino_fp16_url <- paste0(
      "https://huggingface.co/onnx-community/",
      "grounding-dino-tiny-ONNX/resolve/main/onnx/model_fp16.onnx"
    )
  }
  variant <- match.arg(variant, several.ok = TRUE)
  dir <- pliman_model_dir(dir)
  rows <- list()
  for (v in variant) {
    out <- switch(v,
      "sam2.1-encoder-int8" = file.path(dir, "sam2.1.encoder_int8.onnx"),
      "sam2.1-encoder-fp16" = file.path(dir, "sam2.1.encoder_fp16.onnx"),
      "groundingdino-fp16" = file.path(dir, "groundingdino-tiny_fp16.onnx")
    )
    if (file.exists(out) && !isTRUE(overwrite)) {
      if (isTRUE(verbose)) cli::cli_alert_info("{.file {basename(out)}} already exists (use {.code overwrite = TRUE} to rebuild).")
      rows[[v]] <- data.frame(variant = v, path = out, size_mb = round(file.info(out)$size / 2^20, 1),
                              status = "already_exists", stringsAsFactors = FALSE)
      next
    }
    tmp_out <- paste0(out, ".tmp.onnx")
    on.exit(unlink(tmp_out), add = TRUE)
    if (isTRUE(verbose)) cli::cli_progress_step("Building {.file {basename(out)}}...", msg_done = "Built {.file {basename(out)}}")
    if (v %in% c("sam2.1-encoder-int8", "sam2.1-encoder-fp16")) {
      pliman_download_model("sam2.1", dir = dir)
      src <- normalizePath(file.path(dir, "sam2.1.encoder.onnx"), winslash = "/")
      code <- if (v == "sam2.1-encoder-int8") .variant_py_sam_int8 else .variant_py_fp16
      .run_variant_python(code, c(src, tmp_out), verbose = verbose)
    } else {
      dl <- file.path(tempdir(), "groundingdino_model_fp16.onnx")
      on.exit(unlink(dl), add = TRUE)
      if (!file.exists(dl) || file.info(dl)$size < 1e8) {
        old_to <- getOption("timeout")
        options(timeout = max(3600, old_to))
        on.exit(options(timeout = old_to), add = TRUE)
        ok <- .download_url_robust(dino_fp16_url, dl, min_size = 1e8)
        if (!isTRUE(ok)) cli::cli_abort("Failed to download {.url {dino_fp16_url}}.")
      }
      .run_variant_python(.variant_py_gridsample_fp32, c(normalizePath(dl, winslash = "/"), tmp_out), verbose = verbose)
    }
    if (!file.exists(tmp_out) || file.info(tmp_out)$size < 1e6) cli::cli_abort("Conversion produced no model for {.val {v}}.")
    file.rename(tmp_out, out)
    rows[[v]] <- data.frame(variant = v, path = out, size_mb = round(file.info(out)$size / 2^20, 1),
                            status = "built", stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  if (isTRUE(verbose)) print(res)
  invisible(res)
}
