#' Inspect Detailed Metadata and Specifications of YOLO Models
#'
#' @description
#' `yolo_model_details()` (and its alias `yolo_details()`) inspects a YOLO model
#' in PyTorch (`.pt`) or ONNX (`.onnx`) format and extracts extensive structural,
#' architectural, and training metadata.
#'
#' It extracts:
#' \itemize{
#'   \item \strong{Format & Provenance:} Model format (`.pt` or `.onnx`), file size,
#'     SHA-256 and MD5 cryptographic hashes, Ultralytics version, export/training date,
#'     author, description, license, and documentation URL.
#'   \item \strong{Model Specifications:} Computer vision task (`"detect"`, `"segment"`,
#'     `"pose"`, `"classify"`, or `"obb"`), architecture definition, total parameter
#'     count, feature map strides, default input image resolution (`imgsz`), batch size,
#'     and input channels.
#'   \item \strong{I/O Tensor Signatures (ONNX):} Input and output tensor names, full
#'     multidimensional shape signatures, and tensor data types.
#'   \item \strong{Classes:} Total number of classes and complete class label taxonomy.
#'   \item \strong{Training Configuration (PyTorch checkpoints):} Epochs trained,
#'     best fitness score, optimizer, learning rates, data source, and hyperparameters.
#'   \item \strong{Export Settings (ONNX):} Target opset, half-precision (FP16),
#'     dynamic axes, ONNX simplification, and non-maximum suppression (NMS) flags.
#' }
#'
#' @param model A character string specifying the path to a `.pt` or `.onnx` model file,
#'   or the name of a built-in pre-trained model (e.g. `"yolo11n"`, `"yolo26n-seg"`).
#' @param verbose Logical. If `TRUE` (default), prints a formatted diagnostic
#'   overview to the R console.
#' @param ... Currently unused.
#'
#' @return An invisible object of class `c("yolo_model_details", "list")` containing:
#' \itemize{
#'   \item `file`: Absolute path to the inspected model file.
#'   \item `file_name`: File basename.
#'   \item `format`: Model file format (`"PyTorch (.pt)"` or `"ONNX (.onnx)"`).
#'   \item `file_size`: Human-readable file size string (e.g. `"5.35 MB"`).
#'   \item `file_size_bytes`: Raw file size in bytes.
#'   \item `sha256`: SHA-256 checksum string.
#'   \item `md5`: MD5 checksum string.
#'   \item `task`: Model task (`"detect"`, `"segment"`, `"pose"`, `"classify"`, `"obb"`).
#'   \item `architecture`: Network architecture or YAML scale descriptor if available.
#'   \item `num_classes`: Total number of recognized object categories.
#'   \item `classes`: Character vector of class names.
#'   \item `parameters`: Total number of model weights/parameters (numeric).
#'   \item `parameters_formatted`: Human-readable parameter count (e.g. `"2.62 M"`).
#'   \item `stride`: Feature map stride value or vector (e.g. `c(8, 16, 32)`).
#'   \item `imgsz`: Default/training image resolution (e.g. `c(640, 640)`).
#'   \item `batch`: Default batch size.
#'   \item `channels`: Number of input color channels (typically `3` for RGB).
#'   \item `inputs`: `data.frame` detailing input tensor names, shapes, and types.
#'   \item `outputs`: `data.frame` detailing output tensor names, shapes, and types.
#'   \item `version`: Ultralytics version used during training or export.
#'   \item `date`: Creation/training timestamp.
#'   \item `author`: Model author or organization.
#'   \item `description`: Model description or training dataset reference.
#'   \item `license`: Model license string.
#'   \item `docs`: Documentation URL.
#'   \item `export_args`: List of ONNX exporter configuration settings.
#'   \item `train_args`: List of training hyperparameters and arguments (if present).
#'   \item `epoch`: Number of completed training epochs.
#'   \item `best_fitness`: Best validation fitness metric achieved during training.
#' }
#'
#' @export
#' @aliases yolo_details
#' @examples
#' \dontrun{
#'   library(pliman)
#'
#'   # Inspect a pre-trained ONNX model
#'   details <- yolo_model_details("yolo11n")
#'
#'   # Inspect a custom trained PyTorch checkpoint
#'   pt_details <- yolo_model_details("runs/detect/train/weights/best.pt")
#'
#'   # Access specific metadata fields
#'   pt_details$task
#'   pt_details$num_classes
#'   pt_details$classes
#'   pt_details$parameters
#'   pt_details$sha256
#' }
yolo_model_details <- function(model, verbose = TRUE, ...) {
  if (missing(model) || is.null(model) || !is.character(model) || length(model) == 0L || !nzchar(model[1])) {
    cli::cli_abort("Please supply a valid model path or pre-trained model name to {.fn yolo_model_details}.")
  }

  req_model <- model[1]

  # 1. Resolve model file path
  resolved_path <- NULL
  if (file.exists(req_model)) {
    resolved_path <- normalizePath(req_model, winslash = "/", mustWork = TRUE)
  } else {
    hit <- .find_local_model(req_model)
    if (!is.null(hit) && file.exists(hit)) {
      resolved_path <- normalizePath(hit, winslash = "/", mustWork = TRUE)
    } else {
      p_dir <- pliman_model_dir()
      candidates <- c(
        file.path(p_dir, req_model),
        file.path(p_dir, paste0(req_model, ".onnx")),
        file.path(p_dir, paste0(req_model, ".pt")),
        file.path("models", req_model),
        file.path("models", paste0(req_model, ".onnx")),
        file.path("models", paste0(req_model, ".pt")),
        file.path("weights", req_model),
        file.path("weights", paste0(req_model, ".onnx")),
        file.path("weights", paste0(req_model, ".pt"))
      )
      for (cand in candidates) {
        if (file.exists(cand)) {
          resolved_path <- normalizePath(cand, winslash = "/", mustWork = TRUE)
          break
        }
      }
    }
  }

  if (is.null(resolved_path) || !file.exists(resolved_path)) {
    cli::cli_abort(c(
      "x" = "Model file {.val {req_model}} was not found.",
      "i" = "Please check the path or verify models downloaded in {.path {pliman_model_dir()}}."
    ))
  }

  ext <- tolower(tools::file_ext(resolved_path))
  if (!ext %in% c("pt", "onnx")) {
    cli::cli_abort("Model format must be either {.val .pt} or {.val .onnx} (received {.val {ext}}).")
  }

  # File properties
  f_info <- file.info(resolved_path)
  f_bytes <- f_info$size
  f_size_mb <- round(f_bytes / (1024 * 1024), 2)
  f_size_fmt <- if (f_size_mb >= 1024) paste0(round(f_size_mb / 1024, 2), " GB") else paste0(f_size_mb, " MB")

  # Cryptographic checksums
  f_sha256 <- tryCatch({
    con <- file(resolved_path, "rb")
    h <- openssl::sha256(con)
    close(con)
    as.character(h)
  }, error = function(e) NA_character_)

  f_md5 <- tryCatch({
    con <- file(resolved_path, "rb")
    h <- openssl::md5(con)
    close(con)
    as.character(h)
  }, error = function(e) NA_character_)

  # Python inspection script
  py_code <- '
import sys, json, os, ast

model_path = sys.argv[1]
ext = os.path.splitext(model_path)[1].lower()

res = {
    "status": "ok",
    "task": None,
    "architecture": None,
    "classes": [],
    "num_classes": 0,
    "parameters": None,
    "stride": None,
    "imgsz": None,
    "batch": None,
    "channels": None,
    "inputs": [],
    "outputs": [],
    "date": None,
    "version": None,
    "author": None,
    "description": None,
    "license": None,
    "docs": None,
    "export_args": None,
    "train_args": None,
    "epoch": None,
    "best_fitness": None
}

if ext == ".pt":
    try:
        from ultralytics import YOLO
        y = YOLO(model_path)
        m = y.model
        res["task"] = getattr(y, "task", None)
        if hasattr(m, "yaml") and isinstance(m.yaml, dict):
            res["architecture"] = m.yaml.get("yaml_file", None) or m.yaml.get("scale", None)
            res["channels"] = m.yaml.get("ch", 3)
        names = y.names
        if isinstance(names, dict):
            res["classes"] = [names[k] for k in sorted(names.keys())]
        elif isinstance(names, list):
            res["classes"] = names
        res["num_classes"] = len(res["classes"])
        res["parameters"] = int(sum(p.numel() for p in m.parameters()))
        stride = getattr(m, "stride", None)
        if hasattr(stride, "tolist"):
            res["stride"] = stride.tolist()
        elif stride is not None:
            res["stride"] = list(stride)
        ckpt = getattr(y, "ckpt", {}) or {}
        if isinstance(ckpt, dict):
            res["epoch"] = ckpt.get("epoch", None)
            res["version"] = ckpt.get("version", None)
            res["date"] = ckpt.get("date", None)
            res["best_fitness"] = float(ckpt["best_fitness"]) if "best_fitness" in ckpt and ckpt["best_fitness"] is not None else None
            t_args = ckpt.get("train_args", {}) or ckpt.get("args", {})
            if isinstance(t_args, dict):
                res["train_args"] = t_args
                if "imgsz" in t_args: res["imgsz"] = t_args["imgsz"]
                if "batch" in t_args: res["batch"] = t_args["batch"]
    except Exception as e:
        res["status"] = "error"
        res["error"] = str(e)

elif ext == ".onnx":
    try:
        import onnxruntime as ort
        sess = ort.InferenceSession(model_path, providers=["CPUExecutionProvider"])
        meta = sess.get_modelmeta()
        custom = meta.custom_metadata_map or {}
        res["task"] = custom.get("task", None)
        res["version"] = custom.get("version", None)
        res["date"] = custom.get("date", None)
        res["author"] = custom.get("author", None)
        res["description"] = custom.get("description", meta.description or None)
        res["license"] = custom.get("license", None)
        res["docs"] = custom.get("docs", None)
        if "stride" in custom:
            try: res["stride"] = ast.literal_eval(custom["stride"])
            except: res["stride"] = int(custom["stride"]) if custom["stride"].isdigit() else custom["stride"]
        if "batch" in custom:
            try: res["batch"] = int(custom["batch"])
            except: res["batch"] = custom["batch"]
        if "channels" in custom:
            try: res["channels"] = int(custom["channels"])
            except: res["channels"] = custom["channels"]
        if "imgsz" in custom:
            try: res["imgsz"] = ast.literal_eval(custom["imgsz"])
            except: res["imgsz"] = custom["imgsz"]
        if "args" in custom:
            try: res["export_args"] = ast.literal_eval(custom["args"])
            except: pass
        if "names" in custom:
            try:
                names_raw = ast.literal_eval(custom["names"])
                if isinstance(names_raw, dict):
                    res["classes"] = [names_raw[k] for k in sorted(names_raw.keys())]
                elif isinstance(names_raw, list):
                    res["classes"] = names_raw
                res["num_classes"] = len(res["classes"])
            except Exception:
                pass
        for inp in sess.get_inputs():
            res["inputs"].append({"name": inp.name, "shape": inp.shape, "type": inp.type})
        for out in sess.get_outputs():
            res["outputs"].append({"name": out.name, "shape": out.shape, "type": out.type})

        # Infer architecture from description if not set
        if not res.get("architecture") and res.get("description"):
            import re
            m_arch = re.search(r"YOLO[0-9a-zA-Z_\\-]+", res["description"], re.IGNORECASE)
            if m_arch:
                res["architecture"] = m_arch.group(0)

        # Count parameters from ONNX initializers
        try:
            import onnx
            m_onnx = onnx.load(model_path, load_external_data=False)
            total_params = 0
            for t in m_onnx.graph.initializer:
                p_cnt = 1
                for d in t.dims:
                    p_cnt *= int(d)
                total_params += p_cnt
            if total_params > 0:
                res["parameters"] = int(total_params)
        except Exception:
            pass

        # Fallback to companion .pt if parameters still missing
        if not res.get("parameters"):
            try:
                pt_cand = os.path.splitext(model_path)[0] + ".pt"
                if not os.path.exists(pt_cand):
                    pt_cand = os.path.join(os.path.dirname(model_path), "best.pt")
                if os.path.exists(pt_cand):
                    from ultralytics import YOLO
                    y_tmp = YOLO(pt_cand)
                    res["parameters"] = int(sum(p.numel() for p in y_tmp.model.parameters()))
            except Exception:
                pass
    except Exception as e:
        res["status"] = "error"
        res["error"] = str(e)

print(json.dumps(res))
'

  parsed <- NULL
  py_res <- tryCatch({
    out <- system2("python", args = c("-c", shQuote(py_code), shQuote(resolved_path)), stdout = TRUE, stderr = FALSE)
    if (length(out) > 0) {
      jsonlite::fromJSON(paste(out, collapse = "\n"))
    } else NULL
  }, error = function(e) NULL)

  if (!is.null(py_res) && isTRUE(py_res$status == "ok")) {
    parsed <- py_res
  } else {
    # Fallback pure R binary extraction for ONNX
    if (ext == "onnx") {
      parsed <- list(
        status = "ok",
        classes = .get_onnx_classes(resolved_path)
      )
      parsed$num_classes <- length(parsed$classes)
      con <- file(resolved_path, "rb")
      seek(con, max(0, f_bytes - min(f_bytes, 524288L)))
      chk <- readBin(con, "raw", n = min(f_bytes, 524288L))
      close(con)
      txt <- rawToChar(chk[chk >= 32 & chk <= 126])

      extract_prop <- function(prop) {
        m <- regmatches(txt, regexec(paste0(prop, "[^a-zA-Z0-9_\\[{]+([a-zA-Z0-9_\\-\\.:/]+)"), txt))[[1]]
        if (length(m) >= 2) m[2] else NULL
      }
      parsed$task <- extract_prop("task")
      parsed$version <- extract_prop("version")
      parsed$date <- extract_prop("date")
      parsed$author <- extract_prop("author")
      parsed$license <- extract_prop("license")
      parsed$stride <- suppressWarnings(as.integer(extract_prop("stride")))
      parsed$channels <- suppressWarnings(as.integer(extract_prop("channels")))
    } else {
      # Fallback pure R for PT
      parsed <- list(
        status = "ok",
        task = "detect",
        classes = character(0),
        num_classes = 0L
      )
    }
  }

  n_p <- if (!is.null(parsed$parameters)) as.numeric(parsed$parameters) else NULL
  if (is.null(n_p) || is.na(n_p)) {
    pt_cand <- sub("\\.onnx$", ".pt", resolved_path, ignore.case = TRUE)
    if (!file.exists(pt_cand)) {
      pt_cand <- file.path(dirname(resolved_path), "best.pt")
    }
    if (file.exists(pt_cand)) {
      pt_res <- tryCatch({
        py_pt <- '
import sys, json
try:
    from ultralytics import YOLO
    y = YOLO(sys.argv[1])
    p = int(sum(p.numel() for p in y.model.parameters()))
    print(json.dumps({"parameters": p}))
except:
    pass
'
        out_pt <- system2("python", args = c("-c", shQuote(py_pt), shQuote(pt_cand)), stdout = TRUE, stderr = FALSE)
        if (length(out_pt) > 0) jsonlite::fromJSON(paste(out_pt, collapse = "\n"))$parameters else NULL
      }, error = function(e) NULL)
      if (!is.null(pt_res)) n_p <- as.numeric(pt_res)
    }
  }

  p_fmt <- if (!is.null(n_p)) {
    if (n_p >= 1e6) paste0(round(n_p / 1e6, 2), " M")
    else if (n_p >= 1e3) paste0(round(n_p / 1e3, 1), " K")
    else as.character(n_p)
  } else NULL

  out_list <- list(
    file = resolved_path,
    file_name = basename(resolved_path),
    format = if (ext == "pt") "PyTorch (.pt)" else "ONNX (.onnx)",
    file_size = f_size_fmt,
    file_size_bytes = f_bytes,
    sha256 = as.character(f_sha256),
    md5 = as.character(f_md5),
    task = if (!is.null(parsed$task)) parsed$task else "detect",
    architecture = parsed$architecture,
    num_classes = if (!is.null(parsed$num_classes)) parsed$num_classes else length(parsed$classes),
    classes = if (!is.null(parsed$classes)) parsed$classes else character(0),
    parameters = n_p,
    parameters_formatted = p_fmt,
    stride = parsed$stride,
    imgsz = parsed$imgsz,
    batch = parsed$batch,
    channels = parsed$channels,
    inputs = parsed$inputs,
    outputs = parsed$outputs,
    version = parsed$version,
    date = parsed$date,
    author = parsed$author,
    description = parsed$description,
    license = parsed$license,
    docs = parsed$docs,
    export_args = parsed$export_args,
    train_args = parsed$train_args,
    epoch = if (!is.null(parsed$epoch) && parsed$epoch >= 0) parsed$epoch else NULL,
    best_fitness = parsed$best_fitness
  )

  class(out_list) <- c("yolo_model_details", "list")

  if (isTRUE(verbose)) {
    print(out_list)
  }
  invisible(out_list)
}

#' @rdname yolo_model_details
#' @export
yolo_details <- yolo_model_details

#' @export
print.yolo_model_details <- function(x, ...) {
  cli::cli_h1("YOLO Model Details")

  # File information
  cli::cli_h2("File Information")
  cli::cli_ul(c(
    paste0("File: {.val ", x$file_name, "}"),
    paste0("Path: {.path ", x$file, "}"),
    paste0("Format: {.strong ", x$format, "}"),
    paste0("Size: {.val ", x$file_size, "} (", prettyNum(x$file_size_bytes, big.mark = ","), " bytes)"),
    paste0("SHA-256: {.code ", x$sha256, "}")
  ))

  # Model Specifications
  cli::cli_h2("Model Specifications")
  task_badge <- if (!is.null(x$task)) toupper(x$task) else "UNKNOWN"

  specs <- c(
    paste0("Task: {.strong ", task_badge, "}"),
    if (!is.null(x$architecture)) paste0("Architecture: {.val ", x$architecture, "}") else NULL,
    if (!is.null(x$parameters)) paste0("Parameters: {.strong ", x$parameters_formatted, "} (", prettyNum(x$parameters, big.mark = ","), " weights)") else NULL,
    if (!is.null(x$stride)) paste0("Stride: {.val ", paste(x$stride, collapse = ", "), "}") else NULL,
    if (!is.null(x$imgsz)) paste0("Input resolution (imgsz): {.val ", paste(x$imgsz, collapse = "x"), "}") else NULL,
    if (!is.null(x$channels)) paste0("Input channels: {.val ", x$channels, " (RGB)}") else NULL,
    if (!is.null(x$batch)) paste0("Batch size: {.val ", x$batch, "}") else NULL
  )
  cli::cli_ul(specs)

  # Tensor I/O Signatures (for ONNX)
  if (!is.null(x$inputs) && is.data.frame(x$inputs) && nrow(x$inputs) > 0) {
    cli::cli_h2("Tensor I/O Signatures")
    in_desc <- vapply(seq_len(nrow(x$inputs)), function(i) {
      sh <- paste0("[", paste(x$inputs$shape[[i]], collapse = ", "), "]")
      paste0("Input '", x$inputs$name[i], "': ", sh, " (", x$inputs$type[i], ")")
    }, character(1))
    out_desc <- vapply(seq_len(nrow(x$outputs)), function(i) {
      sh <- paste0("[", paste(x$outputs$shape[[i]], collapse = ", "), "]")
      paste0("Output '", x$outputs$name[i], "': ", sh, " (", x$outputs$type[i], ")")
    }, character(1))
    cli::cli_ul(c(in_desc, out_desc))
  }

  # Classes
  cli::cli_h2(paste0("Classes (", x$num_classes, ")"))
  if (x$num_classes > 0) {
    if (x$num_classes <= 20) {
      cls_str <- paste0("[", seq_along(x$classes) - 1L, "] ", x$classes)
      cli::cli_ul(cls_str)
    } else {
      preview <- paste0("[", 0:9, "] ", x$classes[1:10])
      cli::cli_ul(c(preview, paste0("... and ", x$num_classes - 10, " more classes (access via {.code $classes})")))
    }
  } else {
    cli::cli_alert_info("No embedded class names found.")
  }

  # Training Configuration (PyTorch)
  if (!is.null(x$train_args) && is.list(x$train_args) && length(x$train_args) > 0) {
    t_items <- c(
      if (!is.null(x$train_args$data)) paste0("Dataset config: {.path ", x$train_args$data, "}") else NULL,
      if (!is.null(x$train_args$epochs)) paste0("Target epochs: {.val ", x$train_args$epochs, "}") else NULL,
      if (!is.null(x$train_args$optimizer)) paste0("Optimizer: {.val ", x$train_args$optimizer, "}") else NULL,
      if (!is.null(x$train_args$lr0)) paste0("Initial learning rate (lr0): {.val ", x$train_args$lr0, "}") else NULL,
      if (!is.null(x$train_args$device)) paste0("Training device: {.val ", x$train_args$device, "}") else NULL,
      if (!is.null(x$train_args$amp)) paste0("Mixed precision (AMP): {.val ", x$train_args$amp, "}") else NULL
    )
    if (length(t_items) > 0) {
      cli::cli_h2("Training Configuration")
      cli::cli_ul(t_items)
    }
  }

  # Provenance & Training Metadata
  prov <- c(
    if (!is.null(x$author)) paste0("Author: {.val ", x$author, "}") else NULL,
    if (!is.null(x$version)) paste0("Ultralytics Version: {.val ", x$version, "}") else NULL,
    if (!is.null(x$date)) paste0("Date: {.val ", x$date, "}") else NULL,
    if (!is.null(x$epoch)) paste0("Completed Epochs: {.val ", x$epoch, "}") else NULL,
    if (!is.null(x$best_fitness)) paste0("Best Fitness: {.val ", round(x$best_fitness, 4), "}") else NULL,
    if (!is.null(x$description)) paste0("Description: {.val ", x$description, "}") else NULL,
    if (!is.null(x$license)) paste0("License: {.val ", x$license, "}") else NULL,
    if (!is.null(x$docs)) paste0("Docs: {.url ", x$docs, "}") else NULL
  )
  if (length(prov) > 0) {
    cli::cli_h2("Provenance & Training Metadata")
    cli::cli_ul(prov)
  }

  # Export arguments
  if (!is.null(x$export_args) && length(x$export_args) > 0) {
    ea <- vapply(names(x$export_args), function(nm) {
      paste0(nm, ": ", deparse(x$export_args[[nm]]))
    }, character(1))
    cli::cli_h2("Export Arguments")
    cli::cli_ul(ea)
  }

  invisible(x)
}

#' @export
summary.yolo_model_details <- function(object, ...) {
  print(object, ...)
}
