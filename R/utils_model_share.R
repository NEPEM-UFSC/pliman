# utils_model_share.R - Community Model Sharing and Hugging Face Pull Request Manager in pliman

#' Contribute and Share a Trained Deep Learning Model to the Global Community
#'
#' @description
#' `pliman_model_share()` (and its aliases `yolo_model_share()` and `pliman_share_model()`)
#' provides a seamless, one-command pipeline for researchers and practitioners to contribute
#' fine-tuned deep learning models (object detectors or segmentation models) to the official
#' NEPEM Hugging Face Model Zoo (\url{https://huggingface.co/datasets/NEPEM-UFSC/dl_models}).
#'
#' The function automatically converts PyTorch (`.pt`) models to portable ONNX
#' format, computes cryptographic SHA-256 checksums, runs visual inference on a
#' sample image, uploads binary files to Git LFS, and creates an automated
#' Community Pull Request on Hugging Face for maintainer review.
#'
#' Once approved and merged by the NEPEM team, your model becomes instantly available to
#' researchers worldwide with a single R command:
#' \preformatted{
#' pliman_download_model("your_model_name")
#' image_detect_dl("image.jpg", model = "your_model_name")
#' }
#'
#' @section Step-by-Step Workflow (Passo a Passo):
#' To successfully publish your model for global community use, follow this standard workflow:
#' \enumerate{
#'   \item \strong{Train your model:} Train your YOLO or deep learning model using [yolo_train()]
#'         or export weights from your favorite training pipeline.
#'   \item \strong{Benchmark and evaluate:} Run [yolo_model_test()] on an independent validation
#'         or test set. Passing the resulting object to `test_res` automatically embeds complete
#'         scientific metrics (mAP@50, mAP@50-95, Precision, Recall, F1, \eqn{R^2}, MAE) into the
#'         model card.
#'   \item \strong{Select a sample image:} Choose a representative, high-quality sample image
#'         (`sample_image`) showcasing typical targets under standard lighting. `pliman` will
#'         run inference on it and generate a side-by-side visual reproduction card.
#'   \item \strong{Get a free Hugging Face token:}
#'         \itemize{
#'           \item Create a free account at \url{https://huggingface.co}.
#'           \item Go to \strong{Settings -> Access Tokens} (\url{https://huggingface.co/settings/tokens}).
#'           \item Click \strong{Create New Token}, select the \strong{Write} role, and copy it.
#'           \item Pass it via the `token` argument or save it in your `.Renviron` as `HF_TOKEN="hf_..."`.
#'         }
#'   \item \strong{Share the model:} Call `pliman_model_share(...)`.
#'   \item \strong{Review the Pull Request:} Open the URL printed in the console to track discussion
#'         and maintainer review. Once approved, the model is published!
#' }
#'
#' @section What to Do (Boas Praticas / Do's):
#' \itemize{
#'   \item \strong{Use clear snake_case names:} Name your model descriptively combining the target,
#'         species, and vision task (e.g., `"flax_capsule_det"`, `"soybean_rust_seg"`, `"corn_seed_count"`).
#'   \item \strong{Always provide \code{test_res}:} Including the [yolo_model_test()] output provides
#'         vital scientific transparency so users know the model's reliability before downloading.
#'   \item \strong{Provide full attribution:} Include `author`, `institution`, and `orcid` so you and
#'         your institution receive proper academic credit and citations.
#'   \item \strong{Specify architecture:} Clearly indicate the model architecture (e.g., `"YOLO26l"`,
#'         `"YOLO11s"`, `"RT-DETR"`).
#'   \item \strong{Choose a balanced sample image:} Pick an image with clear targets and realistic
#'         background conditions (ideally 1 to 3 MB, JPEG or PNG).
#' }
#'
#' @section What NOT to Do (O que Evitar / Don'ts):
#' \itemize{
#'   \item \strong{DO NOT use a Read-only token:} Your Hugging Face token MUST have \strong{Write}
#'         permissions to open Pull Requests and stream weights. A Read token will return an HTTP 403 or 400 error.
#'   \item \strong{DO NOT use vague names:} Avoid ambiguous names like `"model"`, `"best"`, `"test"`,
#'         `"weights"`, or names containing spaces and special characters.
#'   \item \strong{DO NOT upload untested models:} Untested models with poor detection accuracy or
#'         missing metadata will be rejected during community review.
#'   \item \strong{DO NOT upload massive raw video or 50+ MB images:} Keep the `sample_image` focused
#'         and lightweight for web rendering.
#'   \item \strong{DO NOT worry about converting \code{.pt} to \code{.onnx}:} `pliman` detects PyTorch
#'         files and exports them automatically using native ONNX conversion routines.
#'   \item \strong{DO NOT install Python or Git:} The entire upload, SHA-256 calculation, and Git LFS
#'         batch protocol run natively in R and C++.
#' }
#'
#' @param model Path to the model weights (`.onnx` or `.pt` file), or a trained `yolo_train` object.
#'   If a PyTorch `.pt` file is passed, it is automatically exported to ONNX format via [yolo_export()].
#' @param sample_image Character string. Path to a representative test image (required).
#'   Both the raw original image (`{model_name}_sample.jpg`) and the annotated detection
#'   preview (`{model_name}_preview.jpg`) are uploaded to Hugging Face to ensure
#'   full visual verification and scientific reproducibility.
#' @param name Character string. A unique identifier for the model (e.g., `"flax_capsule_det"`).
#'   If `NULL` (default), derived automatically from the model filename. Spaces and special characters
#'   are automatically converted to lowercase snake_case.
#' @param architecture Character string describing the neural network architecture (e.g., `"YOLO26l"`,
#'   `"YOLO11s"`, `"YOLOv8"`, `"RT-DETR"`). If `NULL` (default), auto-detected from model metadata,
#'   filename, or train configuration, defaulting to `"YOLO"`.
#' @param classes Character vector of class names detected by the model (e.g. `c("capsula")`).
#'   If `NULL` (default), auto-extracted from `test_res`, dataset metadata, or model configuration.
#' @param task Character string. Vision task: `"detect"` (default), `"segment"`, or `"classify"`.
#'   If `NULL`, auto-detected.
#' @param author Character string. Name(s) of the model creator(s) or research team (e.g., `"Tiago Olivoto & Mateus Zunino"`).
#' @param institution Character string. Institution, laboratory, or research group (default `"NEPEM - UFSC"`).
#' @param email Optional character string. Contact email for model attribution.
#' @param orcid Optional character string. ORCID identifier of the main author.
#' @param species Optional character string. Plant, crop, or biological species (e.g., `"Linum usitatissimum"`).
#' @param organ Optional character string. Target organ, symptom, or structure (e.g., `"capsule"`, `"seed"`, `"leaf"`, `"rust lesion"`).
#' @param description Optional character string. Detailed description of the model, imaging setup, dataset size, and application domain.
#' @param test_res Optional object of class `yolo_test` (produced by [yolo_model_test()]).
#'   If provided, precision, recall, F1-score, mAP50, mAP50-95, \eqn{R^2}, and MAE are automatically
#'   extracted and embedded in the Hugging Face Model Card and metadata.
#' @param pt_file Optional character string. Direct path to PyTorch weights (`.pt`). If `NULL` (default),
#'   it is automatically searched alongside the ONNX model, in `test_res`, or in [pliman_model_dir()].
#'   If found, it is uploaded directly to the repository root alongside the `.onnx` file.
#' @param upload_assets Logical. If `TRUE` (default), uploads metadata JSON and visual preview images
#'   into the dedicated `assets/` subfolder, ensuring the repository root remains clean with only
#'   model weights (`.onnx` and `.pt`). If `FALSE`, only model weights are uploaded.
#' @param repo_id Character string. Target Hugging Face repository ID. Defaults to `"NEPEM-UFSC/dl_models"`.
#' @param repo_type Character string. Hugging Face repository type: `"dataset"` (default for `dl_models`) or `"model"`.
#' @param token Character string. Hugging Face user access token (with `Write` permissions).
#'   If `NULL` (default), checked in environment variables (`HF_TOKEN`, `HUGGINGFACE_TOKEN`) or R options (`hf_token`).
#'   If still missing, prompted interactively. Free tokens can be created at \url{https://huggingface.co/settings/tokens}.
#' @param create_pr Logical. If `TRUE` (default), creates a Community Pull Request on Hugging Face for
#'   moderation and review. If `FALSE`, attempts a direct push (for repository maintainers with write access).
#' @param auto_install Logical. Kept for backwards compatibility (uploads are executed
#'   natively in R via HTTP REST / Git LFS without requiring Python). Defaults to `TRUE`.
#' @param verbose Logical. Whether to show progress steps and formatted messages. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [image_detect_dl()] or [image_segment_dl()] when rendering `sample_image`.
#'
#' @return An S3 object of class `c("pliman_model_share", "list")` containing:
#'   * `name`: Model identifier string.
#'   * `architecture`: Neural network architecture string (e.g., `"YOLO26l"`).
#'   * `pr_url`: URL of the opened Pull Request / Discussion on Hugging Face.
#'   * `pr_number`: Pull Request number on Hugging Face.
#'   * `files`: Character vector of uploaded files (`.onnx`, `.pt`, `assets/...`).
#'   * `metadata`: Complete structured metadata list.
#'   * `repo_id`: Destination repository ID.
#'
#' @aliases yolo_model_share pliman_share_model
#' @seealso [yolo_model_test()], [yolo_export()], [pliman_download_model()], [image_detect_dl()]
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # ===========================================================================
#' # Recommended Complete Workflow
#' # ===========================================================================
#'
#' # Step 1: Benchmark model on your test dataset
#' test_res <- yolo_model_test(
#'   model = "capsule_model/weights/best.pt",
#'   data = "dataset_capsules"
#' )
#'
#' # Step 2: Share model to the NEPEM community repository
#' share_res <- pliman_model_share(
#'   model = "capsule_model/weights/best.pt",
#'   sample_image = "sample_test_image.jpg",
#'   name = "flax_capsule_det",
#'   architecture = "YOLO26l",
#'   classes = "capsula",
#'   task = "detect",
#'   author = "Tiago Olivoto & Mateus Zunino",
#'   institution = "NEPEM - UFSC",
#'   species = "Linum usitatissimum",
#'   organ = "capsule",
#'   description = "High-accuracy YOLO26l model for flax capsule detection and counting.",
#'   test_res = test_res,
#'   token = "hf_..." # Hugging Face token with Write role
#' )
#'
#' # Step 3: Inspect share results and open Pull Request URL
#' print(share_res)
#' browseURL(share_res$pr_url)

pliman_model_share <- function(model,
                               sample_image = NULL,
                               name = NULL,
                               architecture = NULL,
                               classes = NULL,
                               task = c("detect", "segment", "classify"),
                               author = NULL,
                               institution = "NEPEM - UFSC",
                               email = NULL,
                               orcid = NULL,
                               species = NULL,
                               organ = NULL,
                               description = NULL,
                               test_res = NULL,
                               pt_file = NULL,
                               upload_assets = TRUE,
                               repo_id = "NEPEM-UFSC/dl_models",
                               repo_type = c("dataset", "model"),
                               token = NULL,
                               create_pr = TRUE,
                               auto_install = TRUE,
                               verbose = TRUE,
                               ...) {
  task <- match.arg(task)
  repo_type <- match.arg(repo_type)
  t0 <- Sys.time()

  # 1. Hugging Face Access Token Resolution
  hf_token <- token
  if (is.null(hf_token) || !nzchar(hf_token)) {
    hf_token <- Sys.getenv("HF_TOKEN", unset = "")
  }
  if (!nzchar(hf_token)) {
    hf_token <- Sys.getenv("HUGGINGFACE_TOKEN", unset = "")
  }
  if (!nzchar(hf_token)) {
    hf_token <- getOption("hf_token", default = "")
  }

  if (!nzchar(hf_token)) {
    if (interactive()) {
      cli::cli_alert_info("Hugging Face user token required to create a community contribution.")
      cli::cli_text("You can create a free token (role: 'Write') at: {.url https://huggingface.co/settings/tokens}")
      hf_token <- readline("Paste your Hugging Face token (starts with hf_): ")
      hf_token <- trimws(hf_token)
    }
  }

  if (!nzchar(hf_token)) {
    cli::cli_abort(c(
      "x" = "A valid Hugging Face access token is required to contribute a model.",
      "i" = "Create a free token with {.strong Write} permissions at {.url https://huggingface.co/settings/tokens}.",
      "i" = "Provide it via {.code token = 'hf_...'}, or set {.code Sys.setenv(HF_TOKEN = 'hf_...')}."
    ))
  }

  # 2. Validate Required Arguments
  if (missing(model) || is.null(model)) {
    cli::cli_abort("Argument {.arg model} is required. Provide a model weights file (.onnx or .pt) or a {.code yolo_train} object.")
  }

  sample_image_path <- NULL
  if (isTRUE(upload_assets)) {
    if (is.null(sample_image) || !is.character(sample_image) || length(sample_image) != 1L || !nzchar(sample_image)) {
      cli::cli_abort(c(
        "x" = "Argument {.arg sample_image} is required when {.code upload_assets = TRUE}.",
        "i" = "A representative test image is needed for visual verification and side-by-side reproduction.",
        "i" = "Both the original raw image and the annotated preview will be saved under {.path assets/}.",
        "i" = "To upload strictly model weights without visual previews or metadata assets, pass {.code upload_assets = FALSE}."
      ))
    }
    sample_image_path <- normalizePath(sample_image, winslash = "/", mustWork = FALSE)
    if (!file.exists(sample_image_path)) {
      cli::cli_abort("Sample image file not found at {.path {sample_image_path}}.")
    }
  } else if (!is.null(sample_image) && nzchar(sample_image[1])) {
    sample_image_path <- normalizePath(sample_image[1], winslash = "/", mustWork = FALSE)
  }

  # 3. Resolve Model File & Weights
  model_path <- model
  pt_cand <- pt_file
  if (inherits(model, "yolo_train") || is.list(model)) {
    if (is.null(pt_cand)) {
      pt_cand <- model$pt_file %||% model$weights$best %||% model$best_pt %||% NULL
    }
    if (!is.null(model$model_file) && file.exists(model$model_file)) {
      model_path <- model$model_file
    } else if (!is.null(model$weights$best) && file.exists(model$weights$best)) {
      model_path <- model$weights$best
    } else if (!is.null(model$best_onnx) && file.exists(model$best_onnx)) {
      model_path <- model$best_onnx
    } else if (!is.null(model$model)) {
      model_path <- model$model
    }
  }

  if (!is.character(model_path) || length(model_path) != 1L) {
    cli::cli_abort("Invalid {.arg model} argument. Expected a file path or a {.code yolo_train} object.")
  }

  model_path <- normalizePath(model_path, winslash = "/", mustWork = FALSE)
  if (!file.exists(model_path)) {
    cli::cli_abort("Model file not found at {.path {model_path}}.")
  }

  # Auto-convert PyTorch (.pt) to ONNX format
  onnx_file <- NULL
  if (grepl("\\.pt$", model_path, ignore.case = TRUE)) {
    cand_onnx <- sub("\\.pt$", ".onnx", model_path, ignore.case = TRUE)
    if (file.exists(cand_onnx)) {
      onnx_file <- cand_onnx
    } else {
      local_onnx <- file.path(pliman_model_dir(), paste0(tools::file_path_sans_ext(basename(model_path)), ".onnx"))
      if (file.exists(local_onnx)) {
        onnx_file <- local_onnx
      } else {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Exporting PyTorch weights {.file {basename(model_path)}} to standalone ONNX format...")
        }
        onnx_file <- yolo_export(model = model_path, imgsz = 640, auto_install = auto_install)
      }
    }
  } else if (grepl("\\.onnx$", model_path, ignore.case = TRUE)) {
    onnx_file <- model_path
  } else {
    cli::cli_abort("Unsupported model extension. Expected {.file .onnx} or {.file .pt} file.")
  }

  if (is.null(onnx_file) || !file.exists(onnx_file)) {
    cli::cli_abort("Could not resolve valid ONNX model file for contribution.")
  }

  # Model file size
  file_mb <- round(file.size(onnx_file) / (1024 * 1024), 2)

  # Model Name Sanitization
  model_name <- name
  if (is.null(model_name) || !nzchar(model_name)) {
    model_name <- tools::file_path_sans_ext(basename(onnx_file))
  }
  # Convert to clean snake_case slug
  model_name <- tolower(gsub("[^a-zA-Z0-9_]+", "_", trimws(model_name)))
  model_name <- gsub("^_+|_+$", "", model_name)
  if (!nzchar(model_name)) model_name <- "custom_model"

  # Search for PyTorch weights (.pt) if not yet resolved
  if (is.null(pt_cand)) {
    cands_pt <- c(
      sub("\\.onnx$", ".pt", model_path, ignore.case = TRUE),
      file.path(dirname(model_path), "best.pt"),
      file.path(dirname(model_path), paste0(tools::file_path_sans_ext(basename(model_path)), ".pt")),
      file.path(dirname(dirname(model_path)), "weights", "best.pt"),
      file.path(pliman_model_dir(), paste0(tools::file_path_sans_ext(basename(model_path)), ".pt")),
      file.path(pliman_model_dir(), paste0(model_name, ".pt")),
      file.path(dirname(onnx_file), paste0(tools::file_path_sans_ext(basename(onnx_file)), ".pt")),
      file.path(dirname(onnx_file), "best.pt")
    )
    cands_pt <- cands_pt[file.exists(cands_pt) & !dir.exists(cands_pt)]
    if (length(cands_pt) > 0L) {
      pt_cand <- cands_pt[1]
    }
  }

  pt_file <- if (!is.null(pt_cand) && file.exists(pt_cand)) normalizePath(pt_cand, winslash = "/") else NULL

  # 4. Resolve Model Architecture
  dl_args <- list(...)
  if (is.null(architecture)) {
    if (!is.null(dl_args$arquiteture)) {
      architecture <- dl_args$arquiteture
    } else if (!is.null(dl_args$arquitetura)) {
      architecture <- dl_args$arquitetura
    }
  }
  dl_args$arquiteture <- NULL
  dl_args$arquitetura <- NULL

  if (!is.null(architecture) && nzchar(trimws(architecture))) {
    arch_val <- trimws(architecture)
  } else {
    cand_sources <- c(
      if (inherits(model, "yolo_train")) c(model$model, model$args$model, model$yaml$backbone),
      if (!is.null(test_res) && inherits(test_res, "yolo_test")) test_res$info$model,
      name,
      model_name,
      basename(model_path)
    )
    cand_sources <- cand_sources[!is.null(cand_sources) & nzchar(cand_sources)]
    cand_str <- paste(cand_sources, collapse = " ")

    m <- regexpr("yolo\\s*v?([0-9]+[a-z]?)", cand_str, ignore.case = TRUE)
    if (m != -1L) {
      match_val <- regmatches(cand_str, m)
      v <- sub("^yolo\\s*v?", "", match_val, ignore.case = TRUE)
      is_v <- grepl("yolov", match_val, ignore.case = TRUE)
      arch_val <- paste0("YOLO", if (is_v) "v" else "", v)
    } else if (grepl("rt-?detr", cand_str, ignore.case = TRUE)) {
      arch_val <- "RT-DETR"
    } else if (grepl("sam", cand_str, ignore.case = TRUE)) {
      arch_val <- "SAM"
    } else if (grepl("u2net", cand_str, ignore.case = TRUE)) {
      arch_val <- "U2-Net"
    } else {
      arch_val <- "YOLO"
    }
  }
  architecture <- arch_val

  # 5. Extract Classes & Task from Test Results or Arguments
  if (!is.null(test_res) && inherits(test_res, "yolo_test")) {
    if (is.null(classes) && !is.null(test_res$by_class)) {
      classes <- test_res$by_class$class
    }
    if (!is.null(test_res$info$task)) {
      task <- test_res$info$task
    }
  }

  if (is.null(classes) || length(classes) == 0L) {
    classes <- "object"
  }
  classes <- unique(as.character(classes))

  # Author & Institution defaults
  if (is.null(author) || !nzchar(author)) {
    author <- Sys.getenv("USER", Sys.getenv("USERNAME", "pliman Community Contributor"))
  }

  if (isTRUE(verbose)) {
    cli::cli_rule(
      left = cli::col_cyan("pliman: Model Share & Hugging Face Contribution"),
      right = cli::col_cyan(format(t0, "%Y-%m-%d %H:%M"))
    )
    weights_summary <- sprintf("{.strong %s.onnx} ({.val %s} MB)", model_name, file_mb)
    if (!is.null(pt_file)) {
      pt_size_mb <- round(file.size(pt_file) / (1024 * 1024), 2)
      weights_summary <- sprintf("%s + {.strong %s.pt} ({.val %s} MB)", weights_summary, model_name, pt_size_mb)
    }
    cli::cli_bullets(c(
      "*" = sprintf("Model Name: {.strong %s}", model_name),
      "*" = sprintf("Model Weights: %s", weights_summary),
      "*" = sprintf("Architecture: {.val %s} | Vision Task: {.val %s}", architecture, task),
      "*" = sprintf("Target Repo: {.url https://huggingface.co/%s/%s} [%s]", repo_type, repo_id, toupper(repo_type)),
      "*" = sprintf("Target Classes: {.val %s}", paste(classes, collapse = ", ")),
      "*" = sprintf("Author: {.val %s} (%s)", author, institution),
      "*" = sprintf("Assets Folder: %s", if (isTRUE(upload_assets)) "{.path assets/} (metadata & previews)" else "Disabled (weights only)")
    ))
  }

  # 5. Extract Metrics from test_res if provided
  metrics_info <- list()
  if (!is.null(test_res) && inherits(test_res, "yolo_test")) {
    m <- test_res$metrics
    metrics_info <- list(
      test_images = as.integer(m$images),
      total_gt_objects = as.integer(m$total_gt),
      total_predicted = as.integer(m$total_pred),
      precision = as.numeric(m$precision),
      recall = as.numeric(m$recall),
      f1_score = as.numeric(m$f1),
      mean_iou = as.numeric(m$mean_iou),
      map50 = as.numeric(m$map50),
      map75 = as.numeric(m$map75),
      map50_95 = as.numeric(m$map50_95),
      mae = as.numeric(m$mae),
      rmse = as.numeric(m$rmse),
      r2 = as.numeric(m$r2),
      slope = as.numeric(m$slope),
      intercept = as.numeric(m$intercept)
    )
  }

  # 7. Setup Packaging Directory and Prepare Files
  stage_dir <- file.path(tempdir(), paste0("pliman_share_", model_name))
  if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  dest_onnx_path <- file.path(stage_dir, paste0(model_name, ".onnx"))
  file.copy(onnx_file, dest_onnx_path, overwrite = TRUE)

  dest_pt_path <- NULL
  if (!is.null(pt_file) && file.exists(pt_file)) {
    dest_pt_path <- file.path(stage_dir, paste0(model_name, ".pt"))
    file.copy(pt_file, dest_pt_path, overwrite = TRUE)
  }

  dest_sample_path <- NULL
  dest_preview_path <- NULL

  if (isTRUE(upload_assets) && !is.null(sample_image_path) && file.exists(sample_image_path)) {
    # Copy and normalize original raw sample image
    dest_sample_path <- file.path(stage_dir, paste0(model_name, "_sample.jpg"))
    tryCatch({
      im_raw <- image_import(sample_image_path)
      image_export(im_raw, dest_sample_path)
    }, error = function(e) {
      file.copy(sample_image_path, dest_sample_path, overwrite = TRUE)
    })

    # Run inference and generate annotated preview image
    dest_preview_path <- file.path(stage_dir, paste0(model_name, "_preview.jpg"))
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Running inference on sample image for side-by-side reproduction preview...")
    }

    tryCatch({
      args_call <- c(
        list(
          img = sample_image_path,
          model = onnx_file,
          conf_threshold = 0.25,
          plot = FALSE,
          save_image = TRUE,
          dir_processed = stage_dir,
          file_name = paste0(model_name, "_preview.jpg"),
          prefix = ""
        ),
        dl_args
      )
      if (identical(task, "segment")) {
        do.call(image_segment_dl, args_call)
      } else {
        do.call(image_detect_dl, args_call)
      }
      if (file.exists(dest_preview_path) && isTRUE(verbose)) {
        cli::cli_alert_success("Annotated preview image rendered successfully!")
      }
    }, error = function(e) {
      cli::cli_abort("Failed to generate model preview on sample image: {e$message}")
    })

    if (!file.exists(dest_preview_path)) {
      cli::cli_abort("Could not generate preview image at {.path {dest_preview_path}}.")
    }
  }

  # 8. Build Metadata Structure & Model Card Description
  meta_data <- list(
    model_name = model_name,
    filename = paste0(model_name, ".onnx"),
    pt_filename = if (!is.null(dest_pt_path)) paste0(model_name, ".pt") else NULL,
    architecture = architecture,
    sample_image = if (isTRUE(upload_assets) && !is.null(dest_sample_path)) paste0("assets/", model_name, "_sample.jpg") else NULL,
    preview_image = if (isTRUE(upload_assets) && !is.null(dest_preview_path)) paste0("assets/", model_name, "_preview.jpg") else NULL,
    task = task,
    classes = classes,
    num_classes = length(classes),
    size_mb = file_mb,
    pt_size_mb = if (!is.null(dest_pt_path)) round(file.size(dest_pt_path) / (1024 * 1024), 2) else NULL,
    species = species %||% "Not specified",
    organ = organ %||% "Not specified",
    description = description %||% "Automated deep learning model trained with pliman.",
    author = author,
    institution = institution,
    email = email %||% NA_character_,
    orcid = orcid %||% NA_character_,
    date_shared = format(Sys.time(), "%Y-%m-%d"),
    framework = paste0("pliman / ONNX Runtime (", architecture, ")"),
    license = "CC-BY-4.0",
    metrics = metrics_info
  )

  # Markdown PR Description for Hugging Face
  pr_title <- paste0("[Community Model] Add ", model_name)
  size_desc <- paste0(file_mb, " MB (ONNX format)")
  if (!is.null(dest_pt_path)) {
    pt_mb <- round(file.size(dest_pt_path) / (1024 * 1024), 2)
    size_desc <- paste0(size_desc, " | ", pt_mb, " MB (PyTorch .pt format)")
  }

  pr_body <- paste0(
    "## \U0001f33f Model Contribution: `", model_name, "`\n\n",
    "This community model was submitted automatically via `pliman::pliman_model_share()`.\n\n",
    "### \U0001f4cc Model Metadata\n",
    "- **Model Name:** `", model_name, "`\n",
    "- **Architecture:** `", architecture, "`\n",
    "- **Vision Task:** `", toupper(task), "`\n",
    "- **Target Classes:** `", paste(classes, collapse = ", "), "`\n",
    "- **Species / Crop:** ", meta_data$species, "\n",
    "- **Organ / Target:** ", meta_data$organ, "\n",
    "- **Author(s):** ", author, " (", institution, ")\n",
    if (!is.na(meta_data$orcid)) paste0("- **ORCID:** [", meta_data$orcid, "](https://orcid.org/", meta_data$orcid, ")\n") else "",
    "- **Date:** ", meta_data$date_shared, "\n",
    "- **Model Size:** ", size_desc, "\n\n",
    "### \U0001f4dd Description\n",
    meta_data$description, "\n\n"
  )

  if (length(metrics_info) > 0L) {
    fmt_pct <- function(v) if (!is.null(v) && !is.na(v)) sprintf("%.2f%%", v * 100) else "\u2014"
    fmt_4d  <- function(v) if (!is.null(v) && !is.na(v)) sprintf("%.4f", v) else "\u2014"
    fmt_2d  <- function(v) if (!is.null(v) && !is.na(v)) sprintf("%.2f", v) else "\u2014"

    scale_str <- paste0(
      metrics_info$test_images %||% "?", " imgs (",
      if (!is.null(metrics_info$total_gt_objects)) format(metrics_info$total_gt_objects, big.mark = ",") else "?", " objs)"
    )

    pr_body <- paste0(
      pr_body,
      "### \U0001f4ca Test Benchmark Results\n\n",
      "| Precision | Recall | F1-Score | mAP@50 | mAP@50-95 | Mean IoU | Counting R\u00b2 | MAE | RMSE | Test Scale |\n",
      "| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n",
      "| `", fmt_pct(metrics_info$precision), "` | `",
      fmt_pct(metrics_info$recall), "` | `",
      fmt_4d(metrics_info$f1_score), "` | `",
      fmt_4d(metrics_info$map50), "` | `",
      fmt_4d(metrics_info$map50_95), "` | `",
      fmt_4d(metrics_info$mean_iou), "` | `",
      fmt_4d(metrics_info$r2), "` | `",
      fmt_2d(metrics_info$mae), "` | `",
      fmt_2d(metrics_info$rmse), "` | `",
      scale_str, "` |\n\n",
      "<sub>*Note: **Precision (P):** Correctness of detections. **Recall (R):** Proportion of targets detected. **F1-Score:** Harmonic mean of P and R. **mAP@50 / mAP@50-95:** Mean Average Precision at IoU 0.50 and 0.50\u20130.95. **Mean IoU:** Average overlap quality for True Positives. **Counting R\u00b2 / MAE / RMSE:** Correlation and error in object counting per image. **Test Scale:** Number of test images and ground-truth objects evaluated.*</sub>\n\n"
    )
  }

  model_url <- paste0("https://huggingface.co/", repo_type, "s/", repo_id, "/resolve/main/", model_name, ".onnx")
  sample_url <- paste0("https://huggingface.co/", repo_type, "s/", repo_id, "/resolve/main/assets/", model_name, "_sample.jpg")
  preview_url <- paste0("https://huggingface.co/", repo_type, "s/", repo_id, "/resolve/main/assets/", model_name, "_preview.jpg")

  infer_fn <- if (identical(task, "segment")) {
    "image_segment_dl"
  } else if (identical(task, "classify")) {
    "image_classify_dl"
  } else {
    "image_detect_dl"
  }

  if (isTRUE(upload_assets) && !is.null(dest_sample_path) && !is.null(dest_preview_path)) {
    pr_body <- paste0(
      pr_body,
      "### \U0001f5bc\ufe0f Visual Reproduction & Detection Preview\n\n",
      "| \U0001f4f7 Original Raw Image | \U0001f3af Model Detections/Segmentation |\n",
      "| :---: | :---: |\n",
      "| <a href=\"", sample_url, "\"><img src=\"", sample_url, "\" width=\"700\" alt=\"Original Raw Image\"/></a> | <a href=\"", preview_url, "\"><img src=\"", preview_url, "\" width=\"700\" alt=\"Model Detections\"/></a> |\n\n"
    )
  }

  pr_body <- paste0(
    pr_body,
    "### \U0001f52c Quick Test & Reproducibility\n\n",
    "Run an immediate verification directly in R:\n\n",
    "```r\n",
    "library(pliman)\n\n",
    "# 1. Download model weights:\n",
    "# During PR review, download directly from this PR branch:\n",
    "pliman_download_model(\"", model_url, "\", model_name = \"", model_name, "\")\n",
    "# (Once merged into main, you can simply run: pliman_download_model(\"", model_name, "\"))\n\n",
    if (isTRUE(upload_assets) && !is.null(dest_sample_path)) {
      paste0(
        "# 2. Run inference using the sample image:\n",
        "sample_url <- \"", sample_url, "\"\n",
        "res <- ", infer_fn, "(sample_url, model = \"", model_name, "\")\n"
      )
    } else {
      paste0(
        "# 2. Run inference on your test image:\n",
        "res <- ", infer_fn, "(\"my_image.jpg\", model = \"", model_name, "\")\n"
      )
    },
    "```\n\n",
    "### \U0001f680 How to use on your own images\n\n",
    "```r\n",
    "library(pliman)\n\n",
    "# 1. Download the model (if not already downloaded):\n",
    "pliman_download_model(\"", model_name, "\")\n\n",
    "# 2. Run inference on any local image:\n",
    "res <- ", infer_fn, "(\"my_image.jpg\", model = \"", model_name, "\")\n",
    "```\n",
    if (!is.null(dest_pt_path)) {
      paste0(
        "\n> \U0001f4a1 **Tip:** Both ONNX (`", model_name, ".onnx`) and PyTorch (`", model_name, ".pt`) weights are available in this repository.\n\n"
      )
    } else {
      ""
    }
  )

  # Write JSON metadata
  meta_json_path <- NULL
  if (isTRUE(upload_assets)) {
    json_str <- if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::toJSON(meta_data, pretty = TRUE, auto_unbox = TRUE)
    } else {
      .yolo_json_fallback(meta_data)
    }
    meta_json_path <- file.path(stage_dir, paste0(model_name, "_meta.json"))
    writeLines(as.character(json_str), meta_json_path)
  }

  # 9. Execute Upload and Pull Request via Native Hugging Face REST / Git LFS API (Zero Python)
  if (isTRUE(verbose)) {
    action_label <- if (isTRUE(create_pr)) "Creating Community Pull Request on Hugging Face (Native R)..." else "Uploading model to Hugging Face (Native R)..."
    cli::cli_alert_info(action_label)
  }

  upload_res <- .hf_hub_commit(
    repo_id = repo_id,
    repo_type = repo_type,
    token = hf_token,
    model_name = model_name,
    dest_onnx_path = dest_onnx_path,
    dest_pt_path = dest_pt_path,
    meta_json_path = meta_json_path,
    dest_sample_path = dest_sample_path,
    dest_preview_path = dest_preview_path,
    upload_assets = upload_assets,
    pr_title = pr_title,
    pr_body = pr_body,
    create_pr = create_pr,
    verbose = verbose
  )

  pr_url <- upload_res$pr_url
  pr_num <- upload_res$pr_num

  uploaded_files <- c(
    paste0(model_name, ".onnx"),
    if (!is.null(dest_pt_path)) paste0(model_name, ".pt"),
    if (isTRUE(upload_assets)) paste0("assets/", model_name, "_meta.json"),
    if (isTRUE(upload_assets) && !is.null(dest_sample_path)) paste0("assets/", model_name, "_sample.jpg"),
    if (isTRUE(upload_assets) && !is.null(dest_preview_path)) paste0("assets/", model_name, "_preview.jpg")
  )

  out_obj <- structure(
    list(
      name = model_name,
      architecture = architecture,
      pr_url = pr_url,
      pr_number = pr_num,
      files = uploaded_files,
      metadata = meta_data,
      repo_id = repo_id,
      repo_type = repo_type,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    class = c("pliman_model_share", "list")
  )

  if (isTRUE(verbose)) {
    elapsed_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cli::cli_alert_success(sprintf("Model {.strong %s} submitted successfully in {.val %ss}!", model_name, elapsed_sec))
    print(out_obj)
  }

  invisible(out_obj)
}

#' @rdname pliman_model_share
#' @export
yolo_model_share <- pliman_model_share

#' @rdname pliman_model_share
#' @export
pliman_share_model <- pliman_model_share


# ==============================================================================
# S3 Method: print for pliman_model_share
# ==============================================================================

#' @export
print.pliman_model_share <- function(x, ...) {
  cli::cli_h1("\U0001f33f pliman Model Contribution")
  cli::cli_bullets(c(
    "v" = sprintf("Model Name   : {.strong %s}", x$name),
    "v" = sprintf("Architecture : {.val %s}", x$architecture %||% x$metadata$architecture %||% "YOLO"),
    "v" = sprintf("Target Repo  : {.url https://huggingface.co/%ss/%s}", x$repo_type, x$repo_id),
    "v" = sprintf("Files Added  : {.file %s}", paste(x$files, collapse = ", ")),
    "v" = sprintf("Pull Request : {.url %s}", x$pr_url)
  ))

  task_fn <- if (identical(x$metadata$task, "segment")) "image_segment_dl" else if (identical(x$metadata$task, "classify")) "image_classify_dl" else "image_detect_dl"
  cli::cli_h2("Next Steps")
  cli::cli_bullets(c(
    "*" = "Your Pull Request is now awaiting maintainer review on Hugging Face.",
    "*" = "Once approved and merged by the NEPEM team, anyone will be able to run:",
    " " = cli::col_green(sprintf("  pliman_download_model(\"%s\")", x$name)),
    " " = cli::col_green(sprintf("  %s(\"image.jpg\", model = \"%s\")", task_fn, x$name)),
    "*" = sprintf("View discussion & status: {.url %s}", x$pr_url)
  ))
  invisible(x)
}

# Fallback internal JSON serializer if jsonlite is not installed
.yolo_json_fallback <- function(x) {
  format_val <- function(v) {
    if (is.null(v) || (length(v) == 1 && is.na(v))) return("null")
    if (is.logical(v)) return(if (v) "true" else "false")
    if (is.numeric(v)) return(as.character(v))
    if (is.character(v)) return(paste0("\"", gsub("\"", "\\\\\"", v), "\""))
    if (is.list(v)) {
      if (is.null(names(v))) {
        items <- vapply(v, format_val, character(1))
        return(paste0("[", paste(items, collapse = ", "), "]"))
      } else {
        items <- vapply(names(v), function(k) paste0("\"", k, "\": ", format_val(v[[k]])), character(1))
        return(paste0("{\n  ", paste(items, collapse = ",\n  "), "\n}"))
      }
    }
    paste0("\"", as.character(v), "\"")
  }
  format_val(x)
}


# ==============================================================================
# Native Hugging Face Hub Client (Zero Python / Git LFS REST API)
# ==============================================================================

.hf_calc_sha256 <- function(file_path) {
  sha256(file_path, file = TRUE)
}

.hf_calc_base64 <- function(file_path) {
  base64_encode(file_path, file = TRUE)
}

.hf_hub_commit <- function(repo_id,
                           repo_type = "dataset",
                           token,
                           model_name,
                           dest_onnx_path,
                           dest_pt_path = NULL,
                           meta_json_path = NULL,
                           dest_sample_path = NULL,
                           dest_preview_path = NULL,
                           upload_assets = TRUE,
                           pr_title,
                           pr_body,
                           create_pr = TRUE,
                           verbose = TRUE) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg curl} is required for native Hugging Face communication.")
  }

  url_prefix <- if (identical(repo_type, "dataset")) "datasets/" else if (identical(repo_type, "space")) "spaces/" else ""

  # 1. Hugging Face datasets reject committing binary files (.onnx, .pt, .jpg, .png) as raw git blobs.
  # Therefore, all binary files are uploaded via the Git LFS batch protocol.
  lfs_files <- list(
    list(rel_path = paste0(model_name, ".onnx"), abs_path = dest_onnx_path, label = "ONNX model weights")
  )

  if (!is.null(dest_pt_path) && file.exists(dest_pt_path)) {
    lfs_files <- c(lfs_files, list(
      list(rel_path = paste0(model_name, ".pt"), abs_path = dest_pt_path, label = "PyTorch .pt model weights")
    ))
  }

  if (isTRUE(upload_assets)) {
    if (!is.null(dest_sample_path) && file.exists(dest_sample_path)) {
      lfs_files <- c(lfs_files, list(
        list(rel_path = paste0("assets/", model_name, "_sample.jpg"), abs_path = dest_sample_path, label = "sample image")
      ))
    }
    if (!is.null(dest_preview_path) && file.exists(dest_preview_path)) {
      lfs_files <- c(lfs_files, list(
        list(rel_path = paste0("assets/", model_name, "_preview.jpg"), abs_path = dest_preview_path, label = "preview image")
      ))
    }
  }

  onnx_size <- as.numeric(file.size(dest_onnx_path))
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Computing SHA-256 checksum for model weights ({round(onnx_size / (1024 * 1024), 2)} MB) and images...")
  }

  for (i in seq_along(lfs_files)) {
    lfs_files[[i]]$size <- as.numeric(file.size(lfs_files[[i]]$abs_path))
    lfs_files[[i]]$sha256 <- .hf_calc_sha256(lfs_files[[i]]$abs_path)
  }

  # LFS batch instructions request
  batch_url <- sprintf("https://huggingface.co/%s%s.git/info/lfs/objects/batch", url_prefix, repo_id)
  batch_payload <- list(
    operation = "upload",
    transfers = list("basic"),
    objects = lapply(lfs_files, function(f) list(oid = f$sha256, size = f$size)),
    hash_algo = "sha256",
    ref = list(name = "main")
  )

  h_batch <- curl::new_handle()
  curl::handle_setheaders(
    h_batch,
    .list = list(
      "Accept" = "application/vnd.git-lfs+json",
      "Content-Type" = "application/vnd.git-lfs+json",
      "Authorization" = paste("Bearer", token),
      "User-Agent" = "pliman-R"
    )
  )
  curl::handle_setopt(
    h_batch,
    post = TRUE,
    postfields = if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::toJSON(batch_payload, auto_unbox = TRUE)
    } else {
      .yolo_json_fallback(batch_payload)
    }
  )

  res_batch <- tryCatch({
    curl::curl_fetch_memory(batch_url, handle = h_batch)
  }, error = function(e) {
    cli::cli_abort("Failed to connect to Hugging Face LFS batch service: {e$message}")
  })

  if (res_batch$status_code >= 400L) {
    err_txt <- rawToChar(res_batch$content)
    cli::cli_abort(c(
      "x" = "Hugging Face LFS batch error ({res_batch$status_code}).",
      "i" = "Response: {.val {err_txt}}",
      "i" = "Ensure your token has {.strong Write} permissions on repository {.val {repo_id}}."
    ))
  }

  batch_info <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(rawToChar(res_batch$content), simplifyVector = FALSE)
  } else list()

  objects_res <- if (!is.null(batch_info$objects)) batch_info$objects else list()

  # Upload any objects where upload action is requested by Hugging Face
  for (f_item in lfs_files) {
    matched <- NULL
    for (obj in objects_res) {
      if (identical(obj$oid, f_item$sha256)) {
        matched <- obj
        break
      }
    }

    actions <- if (!is.null(matched)) matched$actions else NULL
    if (!is.null(actions) && !is.null(actions$upload) && !is.null(actions$upload$href)) {
      upload_url <- actions$upload$href
      upload_headers <- actions$upload$header

      if (isTRUE(verbose)) {
        cli::cli_alert_info("Uploading {f_item$label} to Hugging Face LFS storage...")
      }

      headers_vec <- NULL
      if (!is.null(upload_headers) && length(upload_headers) > 0L) {
        headers_vec <- paste0(names(upload_headers), ": ", unlist(upload_headers))
      }

      res_upload <- tryCatch({
        curl::curl_upload(
          file = f_item$abs_path,
          url = upload_url,
          verbose = FALSE,
          httpheader = headers_vec
        )
      }, error = function(e) {
        cli::cli_abort("Failed to stream {f_item$label} to Hugging Face LFS storage: {e$message}")
      })

      if (res_upload$status_code >= 400L) {
        cli::cli_abort("LFS file upload for {f_item$label} failed with HTTP status {res_upload$status_code}.")
      }

      # Verify upload if verify action is provided
      if (!is.null(actions$verify) && !is.null(actions$verify$href)) {
        verify_url <- actions$verify$href
        h_ver <- curl::new_handle()
        curl::handle_setheaders(
          h_ver,
          .list = list(
            "Accept" = "application/vnd.git-lfs+json",
            "Content-Type" = "application/vnd.git-lfs+json",
            "Authorization" = paste("Bearer", token),
            "User-Agent" = "pliman-R"
          )
        )
        curl::handle_setopt(
          h_ver,
          post = TRUE,
          postfields = sprintf('{"oid":"%s","size":%s}', f_item$sha256, f_item$size)
        )
        res_ver <- curl::curl_fetch_memory(verify_url, handle = h_ver)
        if (res_ver$status_code >= 400L) {
          cli::cli_abort("Hugging Face LFS verification for {f_item$label} failed with HTTP status {res_ver$status_code}.")
        }
      }
    } else {
      if (isTRUE(verbose) && identical(f_item$label, "model weights")) {
        cli::cli_alert_info("Model weights already present upstream on Hugging Face (deduplicated).")
      }
    }
  }

  # 2. Prepare NDJSON Commit Payload
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Packaging metadata and files into atomic commit...")
  }

  ndjson_lines <- character()

  # Header line
  h_val <- list(summary = pr_title, description = pr_body)
  h_line <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(list(key = "header", value = h_val), auto_unbox = TRUE)
  } else {
    sprintf('{"key":"header","value":{"summary":"%s","description":"%s"}}',
            gsub("\"", "\\\\\"", pr_title), gsub("\n", "\\\\n", gsub("\"", "\\\\\"", pr_body)))
  }
  ndjson_lines <- c(ndjson_lines, as.character(h_line))

  # LFS Files (.onnx, _sample.jpg, _preview.jpg)
  for (f_item in lfs_files) {
    lfs_val <- list(
      path = f_item$rel_path,
      algo = "sha256",
      oid = f_item$sha256,
      size = f_item$size
    )
    lfs_line <- if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::toJSON(list(key = "lfsFile", value = lfs_val), auto_unbox = TRUE)
    } else {
      sprintf('{"key":"lfsFile","value":{"path":"%s","algo":"sha256","oid":"%s","size":%s}}',
              f_item$rel_path, f_item$sha256, f_item$size)
    }
    ndjson_lines <- c(ndjson_lines, as.character(lfs_line))
  }

  # Regular file helper (for text metadata JSON)
  add_regular_file <- function(rel_path, abs_path) {
    b64_content <- .hf_calc_base64(abs_path)
    r_val <- list(path = rel_path, content = b64_content, encoding = "base64")
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      as.character(jsonlite::toJSON(list(key = "file", value = r_val), auto_unbox = TRUE))
    } else {
      sprintf('{"key":"file","value":{"path":"%s","content":"%s","encoding":"base64"}}',
              rel_path, b64_content)
    }
  }

  # Add Metadata JSON in assets/
  if (isTRUE(upload_assets) && !is.null(meta_json_path) && file.exists(meta_json_path)) {
    ndjson_lines <- c(ndjson_lines, add_regular_file(paste0("assets/", model_name, "_meta.json"), meta_json_path))
  }

  ndjson_payload <- paste(ndjson_lines, collapse = "\n")

  # 3. POST Commit to Hugging Face Hub
  commit_url <- sprintf("https://huggingface.co/api/%s%s/commit/main?create_pr=%s",
                        url_prefix, repo_id, if (isTRUE(create_pr)) "1" else "0")

  h_commit <- curl::new_handle()
  curl::handle_setheaders(
    h_commit,
    .list = list(
      "Content-Type" = "application/x-ndjson",
      "Authorization" = paste("Bearer", token),
      "User-Agent" = "pliman-R"
    )
  )
  curl::handle_setopt(
    h_commit,
    post = TRUE,
    postfields = charToRaw(ndjson_payload)
  )

  res_commit <- tryCatch({
    curl::curl_fetch_memory(commit_url, handle = h_commit)
  }, error = function(e) {
    cli::cli_abort("Failed to send commit to Hugging Face: {e$message}")
  })

  if (res_commit$status_code >= 400L) {
    err_txt <- rawToChar(res_commit$content)
    cli::cli_abort(c(
      "x" = "Failed to create commit on Hugging Face (status {res_commit$status_code}).",
      "i" = "Server response: {.val {err_txt}}",
      "i" = "Ensure your Hugging Face token has {.strong Write} permissions."
    ))
  }

  commit_data <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(rawToChar(res_commit$content))
  } else list()

  pr_url <- commit_data$pullRequestUrl %||% commit_data$commitUrl
  if (is.null(pr_url) || !nzchar(pr_url)) {
    pr_url <- sprintf("https://huggingface.co/%s%s/discussions", url_prefix, repo_id)
  }

  pr_num <- NULL
  if (grepl("/discussions/([0-9]+)", pr_url)) {
    pr_num <- sub(".*/discussions/([0-9]+).*", "\\1", pr_url)
  }

  # 4. Seamless Image URL Resolution in PR Discussion
  # If a PR discussion was opened, Hugging Face renders markdown images from refs/pr/{pr_num}
  if (!is.null(pr_num)) {
    tryCatch({
      disc_url <- sprintf("https://huggingface.co/api/%s%s/discussions/%s", url_prefix, repo_id, pr_num)
      h_disc <- curl::new_handle()
      curl::handle_setheaders(
        h_disc,
        .list = list(
          "Authorization" = paste("Bearer", token),
          "User-Agent" = "pliman-R"
        )
      )
      res_disc <- curl::curl_fetch_memory(disc_url, handle = h_disc)
      if (res_disc$status_code < 400L) {
        disc_json <- jsonlite::fromJSON(rawToChar(res_disc$content))
        events <- disc_json$events
        if (!is.null(events) && nrow(events) > 0L) {
          first_comment_idx <- which(events$type == "comment")
          if (length(first_comment_idx) > 0L) {
            comment_id <- events$id[first_comment_idx[1]]
            comment_content <- events$data$latest$raw[first_comment_idx[1]]
            if (is.character(comment_content) && grepl("/resolve/main/", comment_content)) {
              fixed_content <- gsub("/resolve/main/", sprintf("/resolve/refs%%2Fpr%%2F%s/", pr_num), comment_content)
              edit_url <- sprintf("https://huggingface.co/api/%s%s/discussions/%s/comment/%s/edit",
                                  url_prefix, repo_id, pr_num, tolower(comment_id))
              h_edit <- curl::new_handle()
              curl::handle_setheaders(
                h_edit,
                .list = list(
                  "Content-Type" = "application/json",
                  "Authorization" = paste("Bearer", token),
                  "User-Agent" = "pliman-R"
                )
              )
              curl::handle_setopt(
                h_edit,
                post = TRUE,
                postfields = charToRaw(jsonlite::toJSON(list(content = fixed_content), auto_unbox = TRUE))
              )
              curl::curl_fetch_memory(edit_url, handle = h_edit)
            }
          }
        }
      }
    }, error = function(e) {
      # Non-critical, ignore preview URL patch error if network hiccup
    })
  }

  list(
    pr_url = pr_url,
    pr_num = pr_num,
    commit_url = commit_data$commitUrl,
    commit_oid = commit_data$commitOid
  )
}
