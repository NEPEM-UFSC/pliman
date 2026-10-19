# ==============================================================================
# COLMAP STRUCTURE-FROM-MOTION & 3D GAUSSIAN SPLATTING PIPELINE FOR PLIMAN
# Automated smartphone photogrammetry, camera calibration, and 3DGS training
# Author: Tiago Olivoto & pliman team
# ==============================================================================

#' @title Detect Dedicated NVIDIA GPU with CUDA Support
#' @name has_dedicated_gpu
#' @description
#' Checks whether the host operating system has an accessible dedicated NVIDIA GPU
#' with CUDA support for accelerated photogrammetry (COLMAP SIFT feature matching)
#' and differentiable 3D Gaussian Splatting optimization.
#'
#' @return Logical: `TRUE` if an accessible NVIDIA GPU is detected, `FALSE` otherwise.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' has_dedicated_gpu()
has_dedicated_gpu <- function() {
  # 1. Direct query via nvidia-smi (cross-platform: Windows and Linux)
  nvsmi <- Sys.which("nvidia-smi")
  if (nzchar(nvsmi)) {
    res <- tryCatch(
      system2(nvsmi, args = "--query-gpu=name --format=csv,noheader", stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    if (length(res) > 0 && any(nzchar(trimws(res)))) {
      return(TRUE)
    }
  }

  # 2. Alternative query on Windows via PowerShell / WMI
  sys_name <- Sys.info()["sysname"]
  if (sys_name == "Windows") {
    wmic_out <- tryCatch(
      system2("powershell",
              args = c("-NoProfile", "-Command", "Get-CimInstance Win32_VideoController | Select-Object -ExpandProperty Name"),
              stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    if (any(grepl("NVIDIA|GeForce|Quadro|RTX|Tesla|A100|H100|TITAN", wmic_out, ignore.case = TRUE))) {
      return(TRUE)
    }
  }

  # 3. Alternative query on Linux via lspci
  if (sys_name == "Linux") {
    lspci_bin <- Sys.which("lspci")
    if (nzchar(lspci_bin)) {
      pci_out <- tryCatch(system2(lspci_bin, stdout = TRUE, stderr = FALSE), error = function(e) character(0))
      if (any(grepl("NVIDIA", pci_out, ignore.case = TRUE))) {
        return(TRUE)
      }
    }
  }

  return(FALSE)
}

#' @title Locate COLMAP Executable in the System
#' @name find_colmap
#' @description
#' Searches for the COLMAP binary in the global system `PATH`, R package options
#' (`options(pliman.colmap_path = ...)`), and the persistent user data directory
#' managed by `pliman`. Configures necessary library paths in `PATH` automatically.
#'
#' @return Absolute path to the COLMAP executable or `NULL` if not found.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' colmap_bin <- find_colmap()
find_colmap <- function() {
  # Helper to configure PATH if Windows has lib/
  setup_path <- function(p) {
    p_norm <- normalizePath(p, winslash = "/", mustWork = TRUE)
    bin_dir <- dirname(p_norm)
    lib_dir <- file.path(dirname(bin_dir), "lib")
    curr_path <- Sys.getenv("PATH")
    sep <- if (Sys.info()["sysname"] == "Windows") ";" else ":"
    dirs_to_add <- character(0)
    if (dir.exists(lib_dir) && !grepl(lib_dir, curr_path, fixed = TRUE)) {
      dirs_to_add <- c(dirs_to_add, lib_dir)
    }
    if (dir.exists(bin_dir) && !grepl(bin_dir, curr_path, fixed = TRUE)) {
      dirs_to_add <- c(dirs_to_add, bin_dir)
    }
    if (length(dirs_to_add) > 0) {
      Sys.setenv(PATH = paste(c(dirs_to_add, curr_path), collapse = sep))
    }
    p_norm
  }

  # 1. User-specified option
  opt_path <- getOption("pliman.colmap_path")
  if (!is.null(opt_path) && file.exists(opt_path)) {
    return(setup_path(opt_path))
  }

  # 2. Executable in system PATH
  sys_which <- Sys.which(if (Sys.info()["sysname"] == "Windows") "colmap.exe" else "colmap")
  if (nzchar(sys_which) && file.exists(sys_which)) {
    return(setup_path(sys_which))
  }

  # 3. Persistent user data directory
  target_dir <- tools::R_user_dir("pliman", which = "data")
  colmap_dir <- file.path(target_dir, "colmap_bin")

  if (dir.exists(colmap_dir)) {
    pat <- if (Sys.info()["sysname"] == "Windows") "^colmap\\.exe$" else "^colmap$"
    candidates <- list.files(colmap_dir, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(candidates) > 0 && file.exists(candidates[1])) {
      return(setup_path(candidates[1]))
    }
  }

  return(NULL)
}

#' @title Download and Install Official COLMAP Binaries
#' @name install_colmap
#' @description
#' Detects the operating system and dedicated NVIDIA GPU (CUDA) availability.
#' Downloads official pre-compiled COLMAP binaries, unpacks them into the safe
#' persistent user data directory of `pliman`, and injects the executable and
#' library directories into the session `PATH`.
#'
#' @param force Logical. If `TRUE`, forces re-download even if the executable already exists. Defaults to `FALSE`.
#' @param version Character string. Target COLMAP GitHub release tag. Defaults to `"4.2.1"`.
#' @return Absolute path to the installed COLMAP binary.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' install_colmap()
install_colmap <- function(force = FALSE, version = "4.2.1") {
  sys_name <- Sys.info()["sysname"]
  has_cuda <- has_dedicated_gpu()

  target_dir <- tools::R_user_dir("pliman", which = "data")
  colmap_dir <- file.path(target_dir, "colmap_bin")
  dir.create(colmap_dir, recursive = TRUE, showWarnings = FALSE)

  cli::cli_alert_info("Operating system: {.val {sys_name}}")
  cli::cli_alert_info("NVIDIA GPU acceleration: {.val {ifelse(has_cuda, 'Detected (CUDA enabled)', 'Not detected (CPU mode)')}}")

  # Return existing executable if force = FALSE
  existing <- find_colmap()
  if (!is.null(existing) && !force) {
    cli::cli_alert_success("Existing COLMAP installation found at: {.path {existing}}")
    options(pliman.colmap_path = existing)
    return(invisible(existing))
  }

  base_url <- sprintf("https://github.com/colmap/colmap/releases/download/%s", version)
  bin_path <- ""

  if (sys_name == "Windows") {
    zip_flavor <- if (has_cuda) "colmap-x64-windows-cuda.zip" else "colmap-x64-windows-nocuda.zip"
    download_url <- paste0(base_url, "/", zip_flavor)
    dest_zip <- file.path(colmap_dir, zip_flavor)

    cli::cli_alert_info("Downloading official COLMAP binaries for Windows ({.val {zip_flavor}})...")
    ok <- .download_url_robust(download_url, dest_file = dest_zip, min_size = 10 * 1024 * 1024, progress = TRUE)
    if (!ok || !file.exists(dest_zip)) {
      stop("Failed to download COLMAP archive from ", download_url)
    }

    cli::cli_alert_info("Extracting archive into {.path {colmap_dir}}...")
    utils::unzip(dest_zip, exdir = colmap_dir, overwrite = TRUE)
    unlink(dest_zip)

    candidates <- list.files(colmap_dir, pattern = "^(colmap\\.exe|COLMAP\\.bat)$", recursive = TRUE, full.names = TRUE)
    if (length(candidates) > 0) {
      exe_cand <- grep("colmap\\.exe$", candidates, value = TRUE, ignore.case = TRUE)
      bin_path <- if (length(exe_cand) > 0) exe_cand[1] else candidates[1]
    }

  } else if (sys_name == "Linux") {
    sys_colmap <- Sys.which("colmap")
    if (nzchar(sys_colmap) && file.exists(sys_colmap)) {
      bin_path <- sys_colmap
    } else {
      stop(
        "On Linux, precompiled standalone binaries are not distributed for COLMAP v", version, ".\n",
        "Please install COLMAP using your system package manager or Conda:\n",
        "  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install colmap\n",
        "  - Conda:         conda install -c conda-forge colmap\n",
        "  - Arch Linux:    sudo pacman -S colmap"
      )
    }

  } else if (sys_name == "Darwin") {
    arch <- Sys.info()["machine"]
    if (arch %in% c("arm64", "aarch64")) {
      zip_flavor <- "colmap-arm64-macos.zip"
      download_url <- paste0(base_url, "/", zip_flavor)
      dest_zip <- file.path(colmap_dir, zip_flavor)

      cli::cli_alert_info("Downloading official COLMAP binaries for macOS Apple Silicon ({.val {zip_flavor}})...")
      ok <- .download_url_robust(download_url, dest_file = dest_zip, min_size = 10 * 1024 * 1024, progress = TRUE)
      if (!ok || !file.exists(dest_zip)) {
        stop("Failed to download COLMAP archive from ", download_url)
      }

      cli::cli_alert_info("Extracting archive into {.path {colmap_dir}}...")
      utils::unzip(dest_zip, exdir = colmap_dir, overwrite = TRUE)
      unlink(dest_zip)

      candidates <- list.files(colmap_dir, pattern = "^colmap$", recursive = TRUE, full.names = TRUE)
      if (length(candidates) > 0) {
        bin_path <- candidates[1]
        Sys.chmod(bin_path, mode = "0755")
      }
    } else {
      sys_colmap <- Sys.which("colmap")
      if (nzchar(sys_colmap) && file.exists(sys_colmap)) {
        bin_path <- sys_colmap
      } else {
        stop("On Intel macOS, please install COLMAP via Homebrew in Terminal:\n  brew install colmap")
      }
    }
  } else {
    stop("Unsupported operating system: ", sys_name)
  }

  if (!nzchar(bin_path) || !file.exists(bin_path)) {
    stop("Failed to locate COLMAP binary after extraction.")
  }

  bin_path <- normalizePath(bin_path, winslash = "/", mustWork = TRUE)
  bin_dir <- normalizePath(dirname(bin_path), winslash = "/", mustWork = TRUE)
  lib_dir <- file.path(dirname(bin_dir), "lib")

  # Inject directories into session PATH
  current_path <- Sys.getenv("PATH")
  sep <- if (sys_name == "Windows") ";" else ":"
  dirs_to_add <- character(0)
  if (dir.exists(lib_dir) && !grepl(lib_dir, current_path, fixed = TRUE)) {
    dirs_to_add <- c(dirs_to_add, lib_dir)
  }
  if (!grepl(bin_dir, current_path, fixed = TRUE)) {
    dirs_to_add <- c(dirs_to_add, bin_dir)
  }
  if (length(dirs_to_add) > 0) {
    Sys.setenv(PATH = paste(c(dirs_to_add, current_path), collapse = sep))
    cli::cli_alert_success("Added COLMAP directories to session PATH: {.path {paste(dirs_to_add, collapse = ', ')}}")
  }

  options(pliman.colmap_path = bin_path)

  # Sanity check
  test_check <- tryCatch(
    system2(bin_path, args = "help", stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )

  if (length(test_check) > 0) {
    cli::cli_alert_success("COLMAP binary verified and ready for execution!")
  } else {
    cli::cli_alert_warning("COLMAP installed, but sanity check returned empty. On Windows, please ensure Microsoft Visual C++ Redistributable 2015-2022 is installed.")
  }

  invisible(bin_path)
}

#' @title Execute COLMAP Structure-from-Motion (SfM) Pipeline
#' @name colmap_sfm
#' @description
#' Executes the complete photogrammetric Structure-from-Motion (SfM) pipeline using COLMAP:
#' 1. Feature extraction with SIFT, ALIKED, or LOMA (`feature_extractor`).
#' 2. Geometric two-view feature matching (`exhaustive_matcher`, `sequential_matcher`, or `spatial_matcher`).
#' 3. Sparse 3D reconstruction and camera orientation (`mapper`).
#' 4. Automated model export to plain text format and PLY sparse point cloud (`model_converter`).
#' 5. Optional dense Multi-View Stereo (MVS) depth estimation and fusion (`patch_match_stereo`, `stereo_fusion`).
#'
#' Automatically detects and adapts CLI arguments between COLMAP 4.x (v4.2.1+) and COLMAP 3.x.
#'
#' @param images Path to directory containing photos or a character vector of image file paths.
#' @param output_dir Target workspace directory for the database and reconstructions.
#'   If `NULL`, a `colmap_workspace` folder will be created inside the image directory.
#' @param camera_model Geometric camera projection model. Defaults to `"OPENCV"`. Other options:
#'   `"SIMPLE_RADIAL"`, `"RADIAL"`, `"PINHOLE"`, `"SIMPLE_PINHOLE"`, `"OPENCV_FISHEYE"`, `"FULL_OPENCV"`.
#' @param single_camera Logical. If `TRUE` (default), assumes all images were captured by the same
#'   device/camera, sharing intrinsic calibration parameters.
#' @param camera_params Optional character string specifying known prior intrinsic calibration parameters
#'   (e.g., `"fx,fy,cx,cy"`). Defaults to `NULL`.
#' @param matching Feature matching mode: `"exhaustive"` (default, matches all image pairs),
#'   `"sequential"` (recommended for continuous ordered video frames or sequential drone flights),
#'   or `"spatial"` (matches images based on GPS/spatial coordinates).
#' @param feature_type Visual feature extractor type. One of `"SIFT"` (default), `"ALIKED"`, or `"LOMA"`.
#'   (ALIKED and LOMA are supported in COLMAP 4.x+).
#' @param max_num_features Integer. Maximum number of visual features extracted per image. Defaults to `8192`.
#' @param max_image_size Integer. Maximum image dimension (in pixels) during feature extraction.
#'   `-1` (default) preserves original image resolution.
#' @param sequential_overlap Integer. Number of consecutive images to match in `"sequential"` matching mode.
#'   Defaults to `10`.
#' @param sequential_quadratic Logical. If `TRUE` (default), uses quadratic overlap for sequential matching,
#'   improving long trajectory tracking.
#' @param loop_detection Logical. If `TRUE`, attempts visual loop closure detection in sequential matching.
#'   Defaults to `FALSE`.
#' @param vocab_tree_path Optional path to a COLMAP vocabulary tree file for visual loop detection. Defaults to `NULL`.
#' @param guided_matching Logical. If `TRUE`, uses geometric guided matching to recover additional weak feature
#'   correspondences. Defaults to `FALSE`.
#' @param min_num_matches Integer. Minimum number of mutual feature matches required to verify an image pair.
#'   Defaults to `15`.
#' @param multiple_models Logical. If `TRUE` (default), allows reconstructing multiple sub-models if the scene
#'   is disconnected, automatically selecting the primary model with the largest point cloud.
#' @param extract_colors Logical. If `TRUE` (default), extracts true RGB colors from images for each 3D tie point.
#' @param ba_refine_focal_length Logical. Refine camera focal length during bundle adjustment. Defaults to `TRUE`.
#' @param ba_refine_principal_point Logical. Refine principal point during bundle adjustment. Defaults to `FALSE`.
#' @param ba_refine_extra_params Logical. Refine distortion parameters during bundle adjustment. Defaults to `TRUE`.
#' @param dense Logical. If `TRUE`, executes dense Multi-View Stereo (MVS) reconstruction
#'   via `image_undistorter`, `patch_match_stereo` (using GPU CUDA), and `stereo_fusion`,
#'   generating a dense point cloud (`fused.ply`). Defaults to `FALSE` (fast sparse SfM).
#' @param dense_geom_consistency Logical. Enforces multi-view geometric consistency during dense stereo.
#'   Defaults to `TRUE`.
#' @param dense_max_image_size Integer. Maximum image size for dense depth estimation. Defaults to `-1`.
#' @param use_gpu Logical. If `TRUE` (default), enables GPU CUDA acceleration when available.
#' @param gpu_index Character or integer. GPU device index. Defaults to `"-1"` (automatic selection).
#' @param num_threads Integer. Number of CPU threads to use. Defaults to `-1` (all available cores).
#' @param custom_args Optional character vector of extra custom CLI flags to pass to COLMAP commands.
#' @param overwrite Logical. If `TRUE`, removes any existing database or previous reconstructions in `output_dir`
#'   to start fresh. Defaults to `FALSE`.
#' @param show_colmap_output Logical. If `FALSE` (default), suppresses raw COLMAP C++ terminal output,
#'   logging it cleanly to `colmap.log` in `output_dir` while displaying clean high-level `cli` messages.
#'   If `TRUE`, streams raw COLMAP terminal output directly to the console.
#' @param colmap_path Optional path to the COLMAP executable. If `NULL`, automatically resolved via
#'   [find_colmap()] or [install_colmap()].
#' @param verbose Logical. If `TRUE` (default), displays step-by-step progress messages.
#' @return An object of class `colmap_sfm` (named list) containing paths to generated assets (`database`,
#'   `sparse_dir`, `recon_dir`, `cameras_txt`, `images_txt`, `points3d_ply`, `points3d_txt`,
#'   `dense_dir`, `fused_ply`) and reconstruction metrics (`n_cameras`, `n_points`).
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' # Fast sparse reconstruction with SIFT & OpenCV camera model
#' sfm_res <- colmap_sfm("my_photos/", camera_model = "OPENCV")
#' print(sfm_res)
#'
#' # Sequential video frame reconstruction
#' sfm_seq <- colmap_sfm("frames/", matching = "sequential", sequential_overlap = 15)
#'
#' # Full dense point cloud (PatchMatch MVS with GPU CUDA)
#' dense_res <- colmap_sfm("my_photos/", camera_model = "OPENCV", dense = TRUE)
colmap_sfm <- function(images,
                       output_dir = NULL,
                       camera_model = c("OPENCV", "SIMPLE_RADIAL", "RADIAL", "PINHOLE", "SIMPLE_PINHOLE", "OPENCV_FISHEYE", "FULL_OPENCV"),
                       single_camera = TRUE,
                       camera_params = NULL,
                       matching = c("exhaustive", "sequential", "spatial"),
                       feature_type = c("SIFT", "ALIKED", "LOMA"),
                       max_num_features = 8192,
                       max_image_size = -1,
                       sequential_overlap = 10,
                       sequential_quadratic = TRUE,
                       loop_detection = FALSE,
                       vocab_tree_path = NULL,
                       guided_matching = FALSE,
                       min_num_matches = 15,
                       multiple_models = TRUE,
                       extract_colors = TRUE,
                       ba_refine_focal_length = TRUE,
                       ba_refine_principal_point = FALSE,
                       ba_refine_extra_params = TRUE,
                       dense = FALSE,
                       dense_geom_consistency = TRUE,
                       dense_max_image_size = -1,
                       use_gpu = TRUE,
                       gpu_index = "-1",
                       num_threads = -1,
                       custom_args = NULL,
                       overwrite = FALSE,
                       show_colmap_output = FALSE,
                       colmap_path = NULL,
                       verbose = TRUE) {
  camera_model <- match.arg(camera_model)
  matching <- match.arg(matching)
  feature_type <- match.arg(feature_type)

  # 1. Resolve COLMAP executable
  if (is.null(colmap_path)) {
    colmap_path <- find_colmap()
    if (is.null(colmap_path)) {
      cli::cli_alert_info("COLMAP binary not found in system. Launching automated download and installation...")
      colmap_path <- install_colmap()
    }
  }

  if (!file.exists(colmap_path)) {
    stop("COLMAP executable not found at: ", colmap_path)
  }

  # Detect COLMAP capabilities (v4+ vs v3)
  help_ext <- tryCatch(
    system2(colmap_path, args = c("feature_extractor", "--help"), stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  v_out <- tryCatch(
    system2(colmap_path, args = "version", stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  is_v4 <- any(grepl("FeatureExtraction\\.use_gpu", help_ext)) || any(grepl("COLMAP [4-9]", v_out)) || !any(grepl("COLMAP [1-3]", v_out))

  # 2. Resolve image directory
  staging_dir <- NULL
  if (length(images) == 1 && dir.exists(images)) {
    img_dir <- normalizePath(images, winslash = "/", mustWork = TRUE)
  } else if (is.character(images) && all(file.exists(images))) {
    # Vector of individual image paths: create temporary staging folder
    staging_dir <- file.path(tempdir(), paste0("pliman_colmap_imgs_", as.integer(Sys.time())))
    dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(images, staging_dir)
    img_dir <- normalizePath(staging_dir, winslash = "/", mustWork = TRUE)
  } else {
    stop("'images' must be a valid directory path or a vector of existing image file paths.")
  }

  # 3. Resolve COLMAP workspace output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(img_dir, "colmap_workspace")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  database_path <- file.path(output_dir, "database.db")
  sparse_dir <- file.path(output_dir, "sparse")
  log_file <- file.path(output_dir, "colmap.log")

  # Handle overwrite or purge zero-byte/corrupted database from aborted previous runs
  if (isTRUE(overwrite)) {
    if (file.exists(database_path)) unlink(database_path)
    if (dir.exists(sparse_dir)) unlink(sparse_dir, recursive = TRUE)
    if (file.exists(log_file)) unlink(log_file)
  } else if (file.exists(database_path) && file.size(database_path) == 0) {
    unlink(database_path)
  }
  dir.create(sparse_dir, recursive = TRUE, showWarnings = FALSE)

  # Helper to execute COLMAP steps with clean cli steps and log redirection
  exec_colmap_step <- function(step_desc, args) {
    if (isTRUE(verbose)) {
      cli::cli_progress_step(step_desc)
    }
    out_target <- if (isTRUE(show_colmap_output)) "" else log_file
    st <- system2(
      colmap_path,
      args = args,
      stdout = out_target,
      stderr = out_target
    )
    if (st != 0) {
      if (!isTRUE(show_colmap_output) && file.exists(log_file)) {
        log_lines <- readLines(log_file, warn = FALSE)
        err_tail <- tail(log_lines[nzchar(trimws(log_lines))], 15)
        cli::cli_alert_danger("COLMAP execution error (exit code {st}):\n{paste(err_tail, collapse = '\n')}")
      }
      stop("COLMAP failed during ", step_desc, " (exit code: ", st, ").")
    }
    invisible(st)
  }

  # Detect GPU
  has_gpu <- isTRUE(use_gpu) && has_dedicated_gpu()
  gpu_flag <- if (has_gpu) "1" else "0"

  if (isTRUE(verbose)) {
    cli::cli_rule("COLMAP Structure-from-Motion (SfM)")
    cli::cli_alert_info("Input images directory: {.path {img_dir}}")
    cli::cli_alert_info("Output workspace: {.path {output_dir}}")
    cli::cli_alert_info("Camera model: {.val {camera_model}} | Shared intrinsics: {.val {single_camera}}")
    cli::cli_alert_info("Feature detector: {.val {feature_type}} | Matcher: {.val {matching}}")
    cli::cli_alert_info("Hardware acceleration: {.val {ifelse(has_gpu, 'NVIDIA GPU (CUDA)', 'CPU Mode')}}")
  }

  # --------------------------------------------------------------------------
  # Step 1: Feature Extraction
  # --------------------------------------------------------------------------
  extract_args <- c(
    "feature_extractor",
    "--database_path", database_path,
    "--image_path", img_dir,
    "--ImageReader.camera_model", camera_model,
    "--ImageReader.single_camera", if (isTRUE(single_camera)) "1" else "0"
  )

  # Filter only real image files to avoid scanning workspace database files or subdirectories
  img_files <- list.files(img_dir, pattern = "\\.(jpg|jpeg|png|tif|tiff|bmp|webp)$", ignore.case = TRUE, recursive = FALSE)
  if (length(img_files) > 0) {
    img_list_path <- file.path(output_dir, "image_list.txt")
    writeLines(img_files, img_list_path)
    extract_args <- c(extract_args, "--image_list_path", img_list_path)
  }

  if (!is.null(camera_params) && nzchar(camera_params)) {
    extract_args <- c(extract_args, "--ImageReader.camera_params", as.character(camera_params))
  }

  if (is_v4) {
    extract_args <- c(
      extract_args,
      "--FeatureExtraction.type", feature_type,
      "--FeatureExtraction.use_gpu", gpu_flag,
      "--FeatureExtraction.gpu_index", as.character(gpu_index),
      "--FeatureExtraction.max_image_size", as.character(max_image_size)
    )
    if (num_threads > 0) {
      extract_args <- c(extract_args, "--FeatureExtraction.num_threads", as.character(num_threads))
    }
  } else {
    extract_args <- c(
      extract_args,
      "--SiftExtraction.use_gpu", gpu_flag,
      "--SiftExtraction.gpu_index", as.character(gpu_index),
      "--SiftExtraction.max_image_size", as.character(max_image_size)
    )
    if (num_threads > 0) {
      extract_args <- c(extract_args, "--SiftExtraction.num_threads", as.character(num_threads))
    }
  }

  if (feature_type == "SIFT") {
    extract_args <- c(extract_args, "--SiftExtraction.max_num_features", as.character(max_num_features))
  }

  if (!is.null(custom_args)) {
    extract_args <- c(extract_args, custom_args)
  }

  n_imgs_label <- if (length(img_files) > 0) sprintf(" (%d images)", length(img_files)) else ""
  exec_colmap_step(paste0("Step 1/4: Extracting visual features (", feature_type, ") with COLMAP", n_imgs_label, "..."), extract_args)

  # --------------------------------------------------------------------------
  # Step 2: Feature Matching
  # --------------------------------------------------------------------------
  matcher_cmd <- switch(
    matching,
    "sequential" = "sequential_matcher",
    "spatial"    = "spatial_matcher",
    "exhaustive_matcher"
  )

  match_args <- c(
    matcher_cmd,
    "--database_path", database_path
  )

  if (is_v4) {
    match_args <- c(
      match_args,
      "--FeatureMatching.use_gpu", gpu_flag,
      "--FeatureMatching.gpu_index", as.character(gpu_index),
      "--FeatureMatching.guided_matching", if (isTRUE(guided_matching)) "1" else "0"
    )
    if (num_threads > 0) {
      match_args <- c(match_args, "--FeatureMatching.num_threads", as.character(num_threads))
    }
  } else {
    match_args <- c(
      match_args,
      "--SiftMatching.use_gpu", gpu_flag,
      "--SiftMatching.gpu_index", as.character(gpu_index),
      "--SiftMatching.guided_matching", if (isTRUE(guided_matching)) "1" else "0"
    )
    if (num_threads > 0) {
      match_args <- c(match_args, "--SiftMatching.num_threads", as.character(num_threads))
    }
  }

  if (matching == "sequential") {
    match_args <- c(
      match_args,
      "--SequentialMatching.overlap", as.character(sequential_overlap),
      "--SequentialMatching.quadratic_overlap", if (isTRUE(sequential_quadratic)) "1" else "0",
      "--SequentialMatching.loop_detection", if (isTRUE(loop_detection)) "1" else "0"
    )
    if (!is.null(vocab_tree_path) && file.exists(vocab_tree_path)) {
      match_args <- c(match_args, "--SequentialMatching.vocab_tree_path", vocab_tree_path)
    }
  }

  exec_colmap_step(paste0("Step 2/4: Performing geometric feature matching (", matching, " matcher)..."), match_args)

  # --------------------------------------------------------------------------
  # Step 3: Sparse Reconstruction & Mapping
  # --------------------------------------------------------------------------
  map_args <- c(
    "mapper",
    "--database_path", database_path,
    "--image_path", img_dir,
    "--output_path", sparse_dir,
    "--Mapper.min_num_matches", as.character(min_num_matches),
    "--Mapper.multiple_models", if (isTRUE(multiple_models)) "1" else "0",
    "--Mapper.extract_colors", if (isTRUE(extract_colors)) "1" else "0",
    "--Mapper.ba_refine_focal_length", if (isTRUE(ba_refine_focal_length)) "1" else "0",
    "--Mapper.ba_refine_principal_point", if (isTRUE(ba_refine_principal_point)) "1" else "0",
    "--Mapper.ba_refine_extra_params", if (isTRUE(ba_refine_extra_params)) "1" else "0"
  )
  if (num_threads > 0) {
    map_args <- c(map_args, "--Mapper.num_threads", as.character(num_threads))
  }

  exec_colmap_step("Step 3/4: Estimating camera poses and sparse tie points (Mapper)...", map_args)

  sub_dirs <- list.dirs(sparse_dir, recursive = FALSE, full.names = TRUE)
  if (length(sub_dirs) == 0) {
    stop("COLMAP could not reconstruct any registered cameras. Verify photo overlap and illumination.")
  }

  # Pick the largest reconstruction (highest registered point count / file size)
  if (length(sub_dirs) > 1) {
    sizes <- vapply(sub_dirs, function(sd) {
      f_pts <- file.path(sd, "points3D.bin")
      f_pts_txt <- file.path(sd, "points3D.txt")
      if (file.exists(f_pts)) file.size(f_pts) else if (file.exists(f_pts_txt)) file.size(f_pts_txt) else 0
    }, numeric(1))
    best_idx <- which.max(sizes)
    recon_dir <- sub_dirs[best_idx]
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Found {length(sub_dirs)} sub-reconstruction(s). Selected primary model: {.path {basename(recon_dir)}} (largest point cloud).")
    }
  } else {
    recon_dir <- sub_dirs[1]
  }

  # --------------------------------------------------------------------------
  # Step 4: Model Export to Plain Text and PLY
  # --------------------------------------------------------------------------
  # Convert binary models (cameras.bin, images.bin, points3D.bin) to TXT
  conv_txt_args <- c(
    "model_converter",
    "--input_path", recon_dir,
    "--output_path", recon_dir,
    "--output_type", "TXT"
  )
  ply_dest <- file.path(recon_dir, "points3D.ply")
  conv_ply_args <- c(
    "model_converter",
    "--input_path", recon_dir,
    "--output_path", ply_dest,
    "--output_type", "PLY"
  )

  exec_colmap_step("Step 4/4: Exporting sparse reconstruction to TXT and PLY formats...", conv_txt_args)
  out_target <- if (isTRUE(show_colmap_output)) "" else log_file
  system2(colmap_path, args = conv_ply_args, stdout = out_target, stderr = out_target)

  cameras_txt <- file.path(recon_dir, "cameras.txt")
  images_txt <- file.path(recon_dir, "images.txt")
  points3d_txt <- file.path(recon_dir, "points3D.txt")

  # Calculate registered camera and point count for user feedback
  n_cameras <- 0
  if (file.exists(images_txt)) {
    img_lines <- readLines(images_txt, warn = FALSE)
    n_cameras <- sum(!grepl("^#", img_lines)) / 2
  }
  n_points <- 0
  if (file.exists(points3d_txt)) {
    pt_lines <- readLines(points3d_txt, warn = FALSE)
    n_points <- sum(!grepl("^#", pt_lines))
  }

  dense_dir <- NULL
  fused_ply <- NULL

  # --------------------------------------------------------------------------
  # Optional Step 5: Dense Multi-View Stereo (MVS) Reconstruction
  # --------------------------------------------------------------------------
  if (isTRUE(dense)) {
    if (isTRUE(verbose)) {
      cli::cli_rule("COLMAP Dense Multi-View Stereo (MVS)")
      cli::cli_alert_info("Step 5a: Undistorting images for dense matching...")
    }

    dense_dir <- file.path(output_dir, "dense")
    dir.create(dense_dir, recursive = TRUE, showWarnings = FALSE)

    undist_args <- c(
      "image_undistorter",
      "--image_path", img_dir,
      "--input_path", recon_dir,
      "--output_path", dense_dir,
      "--output_type", "COLMAP"
    )
    if (dense_max_image_size > 0) {
      undist_args <- c(undist_args, "--max_image_size", as.character(dense_max_image_size))
    }
    if (num_threads > 0) {
      undist_args <- c(undist_args, "--num_threads", as.character(num_threads))
    }
    exec_colmap_step("Step 5a: Undistorting images for dense stereo matching...", undist_args)

    pm_args <- c(
      "patch_match_stereo",
      "--workspace_path", dense_dir,
      "--workspace_format", "COLMAP",
      "--PatchMatchStereo.geom_consistency", "true"
    )
    if (isTRUE(use_gpu)) {
      pm_args <- c(pm_args, "--PatchMatchStereo.gpu_index", "0")
    }
    exec_colmap_step("Step 5b: Executing PatchMatch dense stereo depth estimation on GPU (CUDA)...", pm_args)

    fused_ply <- file.path(dense_dir, "fused.ply")
    fuse_args <- c(
      "stereo_fusion",
      "--workspace_path", dense_dir,
      "--workspace_format", "COLMAP",
      "--input_type", "geometric",
      "--output_type", "PLY",
      "--output_path", fused_ply
    )
    if (num_threads > 0) {
      fuse_args <- c(fuse_args, "--StereoFusion.num_threads", as.character(num_threads))
    }
    exec_colmap_step(paste0("Step 5c: Performing stereo depth fusion into dense point cloud (", basename(fused_ply), ")..."), fuse_args)

    if (isTRUE(verbose)) {
      cli::cli_alert_success("Dense reconstruction complete! Dense point cloud: {.path {fused_ply}}")
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("COLMAP processing completed successfully!")
    cli::cli_alert_success("Calibrated cameras: {.val {n_cameras}} | Sparse 3D tie points: {.val {n_points}}")
    cli::cli_alert_info("Camera poses saved to: {.path {images_txt}}")
    cli::cli_alert_info("Sparse point cloud saved to: {.path {ply_dest}}")
    if (!is.null(fused_ply)) cli::cli_alert_info("Dense point cloud saved to: {.path {fused_ply}}")
  }

  # Clean up staging directory if created
  if (!is.null(staging_dir) && dir.exists(staging_dir)) {
    unlink(staging_dir, recursive = TRUE)
  }

  res <- list(
    database = database_path,
    sparse_dir = sparse_dir,
    recon_dir = recon_dir,
    cameras_txt = cameras_txt,
    images_txt = images_txt,
    points3d_ply = ply_dest,
    points3d_txt = points3d_txt,
    dense_dir = dense_dir,
    fused_ply = fused_ply,
    n_cameras = n_cameras,
    n_points = n_points,
    camera_model = camera_model,
    matching = matching,
    images = img_dir
  )
  class(res) <- "colmap_sfm"
  return(res)
}

#' @title Dense Multi-View Stereo (MVS) Point Cloud Reconstruction via COLMAP
#' @name colmap_dense
#' @description
#' Takes a calibrated Structure-from-Motion (SfM) model (from [colmap_sfm()]),
#' undistorts images, computes multi-view dense stereo depth maps on the GPU
#' via CUDA PatchMatch (`patch_match_stereo`), and fuses the depth maps into a high-density,
#' photorealistic 3D point cloud (`fused.ply`).
#'
#' This delivers dense 3D reconstructions comparable to Agisoft Metashape's "Build Dense Cloud",
#' leveraging dedicated NVIDIA CUDA GPUs in just 2 to 4 minutes.
#'
#' @param sfm An object of class `colmap_sfm` (e.g. from [colmap_sfm()]), or a character path to a COLMAP reconstruction directory.
#' @param images Optional character path to original photos directory. If `NULL`, inferred from `sfm$images`.
#' @param output_dir Destination directory for dense workspace. If `NULL`, created automatically inside `dense` subdirectory.
#' @param max_image_size Integer. Maximum image dimension for dense stereo matching. Defaults to `2000` (provides high density with fast GPU computation).
#' @param window_radius Integer. PatchMatch window radius. Defaults to `5`.
#' @param min_num_consistent Integer. Minimum number of geometrically consistent views required per fused 3D point. Defaults to `2`.
#' @param use_gpu Logical. Utilize NVIDIA CUDA GPU for PatchMatch stereo. Defaults to `TRUE`.
#' @param num_threads Integer. Number of CPU threads. Defaults to `-1` (all available cores).
#' @param verbose Logical. Print progress and execution details. Defaults to `TRUE`.
#' @param show_colmap_output Logical. Stream raw COLMAP logs directly to the R console. Defaults to `FALSE`.
#' @return An updated `colmap_sfm` object with `$dense_dir`, `$fused_ply`, `$n_dense_points`, which can be visualized directly with \code{plot()} or [gs_viewer()].
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Solve camera poses with sparse SfM
#' sfm <- colmap_sfm("photos/")
#'
#' # 2. Generate dense point cloud on GPU (in ~2-3 minutes like Metashape)
#' dense_cloud <- colmap_dense(sfm)
#'
#' # 3. Visualize dense point cloud interactively in 3D
#' plot(dense_cloud)
colmap_dense <- function(sfm,
                         images = NULL,
                         output_dir = NULL,
                         max_image_size = 2000,
                         window_radius = 5,
                         min_num_consistent = 2,
                         use_gpu = TRUE,
                         num_threads = -1,
                         verbose = TRUE,
                         show_colmap_output = FALSE) {
  colmap_path <- find_colmap()
  if (is.null(colmap_path)) {
    stop("COLMAP executable not found. Please install official GPU binaries with install_colmap().")
  }

  recon_dir <- NULL
  img_dir <- images
  is_sfm_obj <- inherits(sfm, "colmap_sfm")

  if (is_sfm_obj) {
    recon_dir <- sfm$recon_dir
    if (is.null(img_dir)) img_dir <- sfm$images
    if (is.null(output_dir)) {
      output_dir <- file.path(dirname(dirname(recon_dir)), "dense")
      if (!dir.exists(dirname(output_dir))) {
        output_dir <- file.path(dirname(recon_dir), "dense")
      }
    }
  } else if (is.character(sfm)) {
    if (dir.exists(file.path(sfm, "sparse", "0"))) {
      recon_dir <- file.path(sfm, "sparse", "0")
      if (is.null(output_dir)) output_dir <- file.path(sfm, "dense")
    } else if (dir.exists(sfm) && file.exists(file.path(sfm, "cameras.txt"))) {
      recon_dir <- sfm
      if (is.null(output_dir)) output_dir <- file.path(dirname(sfm), "dense")
    } else {
      stop("'sfm' must be a 'colmap_sfm' object or a directory containing COLMAP reconstruction files.")
    }
  } else {
    stop("'sfm' must be an object of class 'colmap_sfm' or path to a reconstruction directory.")
  }

  if (is.null(img_dir) || !dir.exists(img_dir)) {
    stop("Original images directory not found. Please specify 'images = path/to/images'.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(output_dir, "colmap_dense.log")
  exec_colmap_step <- function(step_desc, args_vec) {
    if (isTRUE(verbose)) cli::cli_alert_info(step_desc)
    out_target <- if (isTRUE(show_colmap_output)) "" else log_file
    status <- system2(colmap_path, args = args_vec, stdout = out_target, stderr = out_target)
    if (status != 0) {
      log_tail <- if (file.exists(log_file)) paste(utils::tail(readLines(log_file, warn = FALSE), 25), collapse = "\n") else ""
      stop("COLMAP execution failed with exit code ", status, " during: ", args_vec[1], "\n", log_tail)
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_rule("COLMAP Dense Multi-View Stereo (GPU CUDA Accelerated)")
    cli::cli_alert_info("Input model: {.path {recon_dir}}")
    cli::cli_alert_info("Dense workspace: {.path {output_dir}}")
  }

  # Step 1: Image Undistortion
  undist_args <- c(
    "image_undistorter",
    "--image_path", img_dir,
    "--input_path", recon_dir,
    "--output_path", output_dir,
    "--output_type", "COLMAP"
  )
  if (max_image_size > 0) {
    undist_args <- c(undist_args, "--max_image_size", as.character(max_image_size))
  }
  if (num_threads > 0) {
    undist_args <- c(undist_args, "--num_threads", as.character(num_threads))
  }
  exec_colmap_step("Step 1/3: Undistorting images for dense stereo matching...", undist_args)

  # Step 2: PatchMatch Dense Stereo on GPU (CUDA)
  pm_args <- c(
    "patch_match_stereo",
    "--workspace_path", output_dir,
    "--workspace_format", "COLMAP",
    "--PatchMatchStereo.geom_consistency", "true",
    "--PatchMatchStereo.window_radius", as.character(window_radius)
  )
  if (isTRUE(use_gpu)) {
    pm_args <- c(pm_args, "--PatchMatchStereo.gpu_index", "0")
  }
  exec_colmap_step("Step 2/3: Executing PatchMatch dense stereo depth estimation on GPU (CUDA)...", pm_args)

  # Step 3: Stereo Fusion into Dense PLY
  fused_ply <- file.path(output_dir, "fused.ply")
  fuse_args <- c(
    "stereo_fusion",
    "--workspace_path", output_dir,
    "--workspace_format", "COLMAP",
    "--input_type", "geometric",
    "--output_type", "PLY",
    "--output_path", fused_ply,
    "--StereoFusion.min_num_pixels", as.character(min_num_consistent)
  )
  if (num_threads > 0) {
    fuse_args <- c(fuse_args, "--StereoFusion.num_threads", as.character(num_threads))
  }
  exec_colmap_step(paste0("Step 3/3: Fusing depth maps into dense 3D point cloud (", basename(fused_ply), ")..."), fuse_args)

  # Read vertex count
  n_dense <- 0
  if (file.exists(fused_ply)) {
    tryCatch({
      hdr <- readLines(fused_ply, n = 30, warn = FALSE)
      elem_line <- grep("^element vertex ", hdr, value = TRUE)
      if (length(elem_line) > 0) {
        n_dense <- as.numeric(strsplit(elem_line[1], "\\s+")[[1]][3])
      }
    }, error = function(e) NULL)
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Dense reconstruction complete!")
    cli::cli_alert_success("Dense 3D points generated: {.val {format(n_dense, big.mark = ',')}}")
    cli::cli_alert_info("Dense point cloud saved to: {.path {fused_ply}}")
    cli::cli_rule()
  }

  if (is_sfm_obj) {
    sfm$dense_dir <- output_dir
    sfm$fused_ply <- fused_ply
    sfm$n_dense_points <- n_dense
    class(sfm) <- c("colmap_dense", "colmap_sfm")
    return(sfm)
  } else {
    res <- list(
      recon_dir = recon_dir,
      dense_dir = output_dir,
      fused_ply = fused_ply,
      n_dense_points = n_dense,
      images = img_dir
    )
    class(res) <- c("colmap_dense", "colmap_sfm")
    return(res)
  }
}

#' @title Reconstruct a 3D Polygonal Surface Mesh from Photos or Point Clouds
#' @name colmap_mesh
#' @description
#' Reconstructs a solid 3D watertight polygonal surface mesh (triangular faces and vertices)
#' directly from calibrated multi-view photos, an SfM reconstruction (`colmap_sfm`), or a dense
#' point cloud (`colmap_dense`).
#'
#' 100% native execution in C++ using COLMAP with direct NVIDIA GPU (CUDA) acceleration
#' for dense stereo depth estimation, with zero Python or external runtime dependencies.
#'
#' Supports:
#' * **Poisson Surface Meshing (`type = "poisson"`)**: Screened Poisson Surface Reconstruction.
#'   Creates smooth, continuous, watertight closed 3D solid models with interpolated vertex colors.
#' * **Delaunay Meshing (`type = "delaunay"`)**: Delaunay tetrahedral graph-cut meshing.
#'   Preserves sharp geometrical features and exact measured tie points. Can run on dense models
#'   or directly on sparse reconstructions.
#' * **Texture Mapping (`texture = TRUE`)**: Projects original high-resolution camera photos
#'   onto the 3D surface mesh to produce photorealistic texture atlases.
#'
#' @param sfm An object of class `colmap_sfm` (from [colmap_sfm()]), `colmap_dense` (from [colmap_dense()]),
#'   a path to a COLMAP reconstruction/dense directory, or a path to a dense `.ply` point cloud.
#' @param type Character string. Meshing algorithm: `"poisson"` (default, screened Poisson watertight mesh)
#'   or `"delaunay"` (Delaunay graph-cut surface).
#' @param images Optional character string specifying the directory with original photos. If `NULL`,
#'   automatically resolved from `sfm`.
#' @param output_file Optional character string specifying the output `.ply` mesh file path.
#'   Defaults to `"meshed-poisson.ply"` or `"meshed-delaunay.ply"` inside the dense workspace directory.
#' @param dense Logical. If `TRUE` (default for Poisson) and the dense stereo reconstruction has not
#'   yet been performed, automatically executes [colmap_dense()] on GPU first.
#' @param texture Logical. If `TRUE` (default: `FALSE`), runs COLMAP's native `mesh_texturer` to project
#'   original photo colors and generate a textured mesh with texture atlas image.
#' @param depth Integer. Octree depth for Poisson reconstruction (e.g. `10` to `14`). Higher values
#'   resolve finer geometry at the cost of more memory. Defaults to `13`.
#' @param trim Numeric. Trimming threshold for Poisson meshing to prune low-density boundary artifacts.
#'   Defaults to `10`. Set to `0` to keep the entire unclipped reconstruction.
#' @param point_weight Numeric. Weight for the screened Poisson point interpolation constraint. Defaults to `1.0`.
#' @param max_image_size Integer. Maximum dimension for image undistortion. Defaults to `2000`.
#' @param use_gpu Logical. If `TRUE` (default), leverages dedicated NVIDIA GPU (CUDA) for dense stereo matching.
#' @param num_threads Integer. Number of CPU worker threads. Defaults to `-1` (all available CPU cores).
#' @param verbose Logical. If `TRUE` (default), displays step-by-step progress alerts.
#' @param show_colmap_output Logical. If `FALSE` (default), captures raw COLMAP logs into a log file.
#' @return An object of class `colmap_mesh` (named list) containing the path to the reconstructed `.ply` mesh,
#'   vertex count, face/triangle count, and metadata. Compatible with [gs_viewer()] for interactive 3D inspection.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Reconstruct watertight 3D solid mesh from SfM
#' mesh <- colmap_mesh(sfm_completo, type = "poisson")
#'
#' # 2. Or instant Delaunay meshing on sparse tie points
#' mesh_fast <- colmap_mesh(sfm_completo, type = "delaunay")
#'
#' # 3. Visualize in 3D WebGL viewer
#' gs_viewer(mesh)
colmap_mesh <- function(sfm,
                        type = c("delaunay", "poisson"),
                        images = NULL,
                        output_file = NULL,
                        dense = NULL,
                        texture = FALSE,
                        depth = 13,
                        trim = 10,
                        point_weight = 1.0,
                        max_image_size = 2000,
                        use_gpu = TRUE,
                        num_threads = -1,
                        verbose = TRUE,
                        show_colmap_output = FALSE) {
  type <- match.arg(type)

  # For Poisson, dense normals are required; for Delaunay, sparse runs in seconds (~12s)
  if (is.null(dense)) {
    dense <- (type == "poisson")
  }

  colmap_path <- find_colmap()
  if (is.null(colmap_path)) {
    stop("COLMAP executable not found. Please install official GPU binaries with install_colmap().")
  }

  recon_dir <- NULL
  dense_dir <- NULL
  fused_ply <- NULL
  img_dir <- images
  is_sfm_obj <- inherits(sfm, "colmap_sfm")

  if (is_sfm_obj) {
    recon_dir <- sfm$recon_dir
    if (is.null(img_dir)) img_dir <- sfm$images
    dense_dir <- sfm$dense_dir
    fused_ply <- sfm$fused_ply

    if (is.null(dense_dir)) {
      cand1 <- file.path(dirname(dirname(recon_dir)), "dense")
      cand2 <- file.path(dirname(recon_dir), "dense")
      dense_dir <- if (dir.exists(cand1) || !dir.exists(cand2)) cand1 else cand2
    }
  } else if (is.character(sfm)) {
    if (dir.exists(sfm)) {
      if (dir.exists(file.path(sfm, "dense"))) {
        dense_dir <- file.path(sfm, "dense")
        recon_dir <- file.path(sfm, "sparse", "0")
      } else if (file.exists(file.path(sfm, "stereo", "fusion.cfg"))) {
        dense_dir <- sfm
        recon_dir <- file.path(dirname(sfm), "sparse", "0")
      } else if (dir.exists(file.path(sfm, "sparse", "0"))) {
        recon_dir <- file.path(sfm, "sparse", "0")
        dense_dir <- file.path(sfm, "dense")
      } else if (file.exists(file.path(sfm, "cameras.txt"))) {
        recon_dir <- sfm
        dense_dir <- file.path(dirname(sfm), "dense")
      } else {
        recon_dir <- sfm
        dense_dir <- file.path(sfm, "dense")
      }
    } else if (file.exists(sfm)) {
      ext <- tolower(tools::file_ext(sfm))
      if (ext == "ply") {
        fused_ply <- sfm
        dense_dir <- dirname(sfm)
      }
    }
  } else {
    stop("'sfm' must be an object of class 'colmap_sfm', 'colmap_dense', or path to a reconstruction directory / PLY file.")
  }

  if (is.null(dense_dir)) {
    dense_dir <- if (!is.null(recon_dir)) file.path(dirname(recon_dir), "dense") else file.path(getwd(), "dense")
  }
  dir.create(dense_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(fused_ply)) {
    cand_fused <- file.path(dense_dir, "fused.ply")
    if (file.exists(cand_fused)) fused_ply <- cand_fused
  }

  if (is.null(output_file)) {
    output_file <- file.path(dense_dir, paste0("meshed-", type, ".ply"))
  } else {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  }

  log_file <- file.path(dense_dir, "colmap_mesh.log")
  exec_colmap_step <- function(step_desc, args_vec) {
    if (isTRUE(verbose)) cli::cli_alert_info(step_desc)
    out_target <- if (isTRUE(show_colmap_output)) "" else log_file
    status <- system2(colmap_path, args = args_vec, stdout = out_target, stderr = out_target)
    if (status != 0) {
      log_tail <- if (file.exists(log_file)) paste(utils::tail(readLines(log_file, warn = FALSE), 25), collapse = "\n") else ""
      stop("COLMAP execution failed with exit code ", status, " during: ", args_vec[1], "\n", log_tail)
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_rule("COLMAP 3D Polygonal Surface Meshing (Native C++ / GPU)")
    cli::cli_alert_info("Meshing method: {.val {toupper(type)}}")
    cli::cli_alert_info("Output mesh file: {.path {output_file}}")
  }

  # --------------------------------------------------------------------------
  # Step 1: Execute Meshing
  # --------------------------------------------------------------------------
  if (type == "poisson") {
    if (is.null(fused_ply) || !file.exists(fused_ply)) {
      if (isTRUE(dense)) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Poisson meshing requires dense normals. Executing GPU dense stereo pipeline...")
        }
        sfm_dense <- colmap_dense(
          sfm = if (is_sfm_obj) sfm else recon_dir,
          images = img_dir,
          output_dir = dense_dir,
          max_image_size = max_image_size,
          use_gpu = use_gpu,
          num_threads = num_threads,
          verbose = verbose,
          show_colmap_output = show_colmap_output
        )
        fused_ply <- sfm_dense$fused_ply
      } else {
        stop("Poisson meshing requires a dense point cloud with normals ('fused.ply'). Please set 'dense = TRUE' to run PatchMatch stereo automatically on GPU, or specify type = 'delaunay'.")
      }
    }

    poisson_args <- c(
      "poisson_mesher",
      "--input_path", normalizePath(fused_ply, winslash = "/"),
      "--output_path", normalizePath(output_file, winslash = "/", mustWork = FALSE),
      "--PoissonMeshing.point_weight", as.character(point_weight),
      "--PoissonMeshing.depth", as.character(depth),
      "--PoissonMeshing.color", "1",
      "--PoissonMeshing.trim", as.character(trim)
    )
    if (num_threads > 0) {
      poisson_args <- c(poisson_args, "--PoissonMeshing.num_threads", as.character(num_threads))
    }
    exec_colmap_step("Reconstructing watertight 3D polygonal surface mesh (Poisson Mesher)...", poisson_args)
  } else {
    # Delaunay meshing
    if (isTRUE(dense)) {
      if (is.null(fused_ply) || !file.exists(fused_ply)) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Dense Delaunay meshing requires complete dense stereo. Executing GPU dense stereo pipeline...")
        }
        sfm_dense <- colmap_dense(
          sfm = if (is_sfm_obj) sfm else recon_dir,
          images = img_dir,
          output_dir = dense_dir,
          max_image_size = max_image_size,
          use_gpu = use_gpu,
          num_threads = num_threads,
          verbose = verbose,
          show_colmap_output = show_colmap_output
        )
        fused_ply <- sfm_dense$fused_ply
      }
      input_path <- dense_dir
      input_type <- "dense"
    } else {
      # Fast sparse Delaunay meshing directly on calibrated tie points (~12s)
      input_path <- recon_dir
      input_type <- "sparse"
    }

    delaunay_args <- c(
      "delaunay_mesher",
      "--input_path", normalizePath(input_path, winslash = "/"),
      "--input_type", input_type,
      "--output_path", normalizePath(output_file, winslash = "/", mustWork = FALSE)
    )
    if (num_threads > 0) {
      delaunay_args <- c(delaunay_args, "--DelaunayMeshing.num_threads", as.character(num_threads))
    }
    exec_colmap_step(paste0("Reconstructing 3D polygonal Delaunay surface mesh (", input_type, ")..."), delaunay_args)
  }

  # --------------------------------------------------------------------------
  # Step 2: Optional Texture Mapping
  # --------------------------------------------------------------------------
  textured_dir <- NULL
  if (isTRUE(texture)) {
    textured_dir <- file.path(dense_dir, "textured")
    dir.create(textured_dir, recursive = TRUE, showWarnings = FALSE)
    tex_args <- c(
      "mesh_texturer",
      "--workspace_path", normalizePath(dense_dir, winslash = "/"),
      "--input_path", normalizePath(output_file, winslash = "/"),
      "--output_path", normalizePath(textured_dir, winslash = "/", mustWork = FALSE),
      "--output_type", "BIN"
    )
    if (num_threads > 0) {
      tex_args <- c(tex_args, "--MeshTextureMapping.num_threads", as.character(num_threads))
    }
    exec_colmap_step("Baking original photo textures into 3D mesh UV texture atlas (mesh_texturer)...", tex_args)
  }

  # --------------------------------------------------------------------------
  # Read Vertex and Face Count from PLY
  # --------------------------------------------------------------------------
  n_vertices <- 0L
  n_faces <- 0L
  if (file.exists(output_file)) {
    tryCatch({
      hdr <- readLines(output_file, n = 35, warn = FALSE)
      v_line <- grep("^element vertex ", hdr, value = TRUE)
      if (length(v_line) > 0) {
        n_vertices <- as.integer(strsplit(trimws(v_line[1]), "\\s+")[[1]][3])
      }
      f_line <- grep("^element face ", hdr, value = TRUE)
      if (length(f_line) > 0) {
        n_faces <- as.integer(strsplit(trimws(f_line[1]), "\\s+")[[1]][3])
      }
    }, error = function(e) NULL)
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("3D polygonal surface mesh generated successfully!")
    cli::cli_alert_success("Vertices: {.val {format(n_vertices, big.mark = ',')}} | Triangular Faces: {.val {format(n_faces, big.mark = ',')}}")
    cli::cli_alert_info("Mesh file saved to: {.path {output_file}}")
    if (!is.null(textured_dir)) {
      cli::cli_alert_info("Textured mesh and atlas saved to: {.path {textured_dir}}")
    }
    cli::cli_rule()
  }

  res <- list(
    mesh_ply = output_file,
    type = type,
    n_vertices = n_vertices,
    n_faces = n_faces,
    dense_dir = dense_dir,
    recon_dir = recon_dir,
    fused_ply = fused_ply,
    textured_dir = textured_dir,
    images = img_dir,
    sfm = sfm
  )
  class(res) <- c("colmap_mesh", "colmap_sfm")
  res
}

#' @title Print Summary of COLMAP 3D Mesh
#' @name print.colmap_mesh
#' @param x An object of class `colmap_mesh`.
#' @param ... Additional arguments (currently ignored).
#' @export
print.colmap_mesh <- function(x, ...) {
  cli::cli_rule("COLMAP 3D Polygonal Surface Mesh")
  cli::cli_alert_info("Method: {.val {toupper(x$type)}}")
  cli::cli_alert_success("Vertices: {.val {format(x$n_vertices, big.mark = ',')}}")
  cli::cli_alert_success("Triangular Faces: {.val {format(x$n_faces, big.mark = ',')}}")
  cli::cli_alert_info("Mesh file (PLY): {.path {x$mesh_ply}}")
  if (!is.null(x$textured_dir)) {
    cli::cli_alert_info("Textured atlas directory: {.path {x$textured_dir}}")
  }
  invisible(x)
}

#' @title Plot 3D Interactive Viewer for COLMAP Mesh
#' @name plot.colmap_mesh
#' @param x An object of class `colmap_mesh`.
#' @param ... Additional arguments passed to [gs_viewer()].
#' @export
plot.colmap_mesh <- function(x, ...) {
  gs_viewer(x, ...)
}

#' @title Print Summary of COLMAP SfM Reconstruction
#' @name print.colmap_sfm
#' @param x An object of class `colmap_sfm`.
#' @param ... Additional arguments (currently ignored).
#' @export
print.colmap_sfm <- function(x, ...) {
  cli::cli_rule("COLMAP SfM Reconstruction Summary")
  cli::cli_alert_info("Reconstruction directory: {.path {x$recon_dir}}")
  cli::cli_alert_info("Camera model: {.val {x$camera_model}} | Matcher: {.val {x$matching}}")
  if (!is.null(x$n_cameras)) cli::cli_alert_success("Calibrated / registered cameras: {.val {x$n_cameras}}")
  if (!is.null(x$n_points)) cli::cli_alert_success("Sparse 3D tie points: {.val {x$n_points}}")
  if (!is.null(x$n_dense_points) && x$n_dense_points > 0) {
    cli::cli_alert_success("Dense 3D points (MVS): {.val {format(x$n_dense_points, big.mark = ',')}}")
    cli::cli_alert_info("Dense point cloud (PLY): {.path {x$fused_ply}}")
  }
  cli::cli_alert_info("Sparse point cloud (PLY): {.path {x$points3d_ply}}")
  if (!is.null(x$fused_ply) && file.exists(x$fused_ply)) {
    cli::cli_alert_success("Dense point cloud (PLY): {.path {x$fused_ply}}")
  }
  invisible(x)
}

#' @title Plot 3D Interactive Viewer for COLMAP SfM Reconstruction
#' @name plot.colmap_sfm
#' @param x An object of class `colmap_sfm`.
#' @param ... Additional arguments passed to [gs_viewer()].
#' @export
plot.colmap_sfm <- function(x, ...) {
  gs_viewer(x, ...)
}

#' @title Import Calibrated Camera Poses from COLMAP
#' @name gs_import_cameras_colmap
#' @description
#' Reads `cameras.txt` and `images.txt` exported by COLMAP, converts orientation
#' quaternions \eqn{(q_w, q_x, q_y, q_z)} to \eqn{3 \times 3} rotation matrices \eqn{R}, extracts
#' the translation vector \eqn{t}, calculates the 3D optical center \eqn{C = -R^T t}, and constructs
#' camera specifications ready for [gs_fit()].
#'
#' @param images_file Path to COLMAP `images.txt` or the reconstruction directory containing it.
#' @param cameras_file Path to COLMAP `cameras.txt`. If `NULL`, inferred from `images_file`.
#' @return A named list containing camera specifications for each calibrated image.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' cams <- gs_import_cameras_colmap("sparse/0/images.txt")
gs_import_cameras_colmap <- function(images_file, cameras_file = NULL) {
  if (dir.exists(images_file)) {
    cameras_file <- file.path(images_file, "cameras.txt")
    images_file <- file.path(images_file, "images.txt")
  } else if (is.null(cameras_file)) {
    cameras_file <- file.path(dirname(images_file), "cameras.txt")
  }

  if (!file.exists(cameras_file)) stop("cameras.txt file not found: ", cameras_file)
  if (!file.exists(images_file)) stop("images.txt file not found: ", images_file)

  # --------------------------------------------------------------------------
  # 1. Parse cameras.txt (Intrinsics: fx, fy, cx, cy, width, height)
  # --------------------------------------------------------------------------
  cam_lines <- readLines(cameras_file, warn = FALSE)
  cam_lines <- trimws(cam_lines)
  cam_lines <- cam_lines[nzchar(cam_lines) & !startsWith(cam_lines, "#")]

  sensor_map <- list()

  for (line in cam_lines) {
    tokens <- strsplit(line, "\\s+")[[1]]
    if (length(tokens) < 5) next

    cam_id <- tokens[1]
    model <- toupper(tokens[2])
    w <- as.numeric(tokens[3])
    h <- as.numeric(tokens[4])
    params <- as.numeric(tokens[5:length(tokens)])

    fx <- fy <- cx <- cy <- 0.0

    if (model == "SIMPLE_PINHOLE") {
      fx <- fy <- params[1]
      cx <- params[2]
      cy <- params[3]
    } else if (model == "PINHOLE") {
      fx <- params[1]
      fy <- params[2]
      cx <- params[3]
      cy <- params[4]
    } else if (model %in% c("SIMPLE_RADIAL", "RADIAL")) {
      fx <- fy <- params[1]
      cx <- params[2]
      cy <- params[3]
    } else if (model %in% c("OPENCV", "OPENCV_FISHEYE", "FULL_OPENCV")) {
      fx <- params[1]
      fy <- params[2]
      cx <- params[3]
      cy <- params[4]
    } else {
      # Generic fallback
      fx <- fy <- params[1]
      cx <- if (length(params) >= 2) params[2] else (w / 2.0)
      cy <- if (length(params) >= 3) params[3] else (h / 2.0)
    }

    sensor_map[[cam_id]] <- list(width = w, height = h, fx = fx, fy = fy, cx = cx, cy = cy, model = model)
  }

  # --------------------------------------------------------------------------
  # 2. Parse images.txt (Extrinsics: Quaternions QW, QX, QY, QZ, TX, TY, TZ)
  # --------------------------------------------------------------------------
  img_lines <- readLines(images_file, warn = FALSE)
  img_lines <- trimws(img_lines)
  img_lines <- img_lines[nzchar(img_lines) & !startsWith(img_lines, "#")]

  cameras_list <- list()

  # In COLMAP images.txt, each image uses 2 lines:
  # Line 1: IMAGE_ID QW QX QY QZ TX TY TZ CAMERA_ID NAME
  # Line 2: 2D feature points (POINTS2D[])
  idx <- 1
  while (idx <= length(img_lines)) {
    tokens <- strsplit(img_lines[idx], "\\s+")[[1]]
    if (length(tokens) >= 10) {
      qw <- as.numeric(tokens[2])
      qx <- as.numeric(tokens[3])
      qy <- as.numeric(tokens[4])
      qz <- as.numeric(tokens[5])
      tx <- as.numeric(tokens[6])
      ty <- as.numeric(tokens[7])
      tz <- as.numeric(tokens[8])
      cam_id <- tokens[9]
      name <- paste(tokens[10:length(tokens)], collapse = " ")

      # Normalize quaternion
      q_norm <- sqrt(qw^2 + qx^2 + qy^2 + qz^2)
      if (q_norm < 1e-9) q_norm <- 1.0
      w <- qw / q_norm
      x <- qx / q_norm
      y <- qy / q_norm
      z <- qz / q_norm

      # Convert quaternion to 3x3 World-to-Camera rotation matrix
      R <- matrix(c(
        1.0 - 2.0 * (y^2 + z^2), 2.0 * (x * y - z * w),     2.0 * (x * z + y * w),
        2.0 * (x * y + z * w),     1.0 - 2.0 * (x^2 + z^2), 2.0 * (y * z - x * w),
        2.0 * (x * z - y * w),     2.0 * (y * z + x * w),     1.0 - 2.0 * (x^2 + y^2)
      ), nrow = 3, ncol = 3, byrow = TRUE)

      t_vec <- c(tx, ty, tz)

      # Camera optical center in world space: C = -R^T * t
      pos_world <- as.numeric(-t(R) %*% t_vec)

      sens <- if (!is.null(sensor_map[[cam_id]])) sensor_map[[cam_id]] else sensor_map[[1]]

      cam_spec <- list(
        R = R,
        t = t_vec,
        fx = sens$fx,
        fy = sens$fy,
        cx = sens$cx,
        cy = sens$cy,
        width = sens$width,
        height = sens$height,
        position = pos_world,
        label = name
      )

      cameras_list[[name]] <- cam_spec
    }
    # Skip to next image entry (skipping line 2 of 2D points)
    idx <- idx + 2
  }

  cli::cli_alert_success("Imported {length(cameras_list)} calibrated camera pose(s) from COLMAP.")
  return(cameras_list)
}

#' @title Import Complete COLMAP Reconstruction (Cameras + Sparse Points)
#' @name gs_import_colmap
#' @description
#' Simultaneously loads calibrated camera poses and sparse tie points exported by
#' COLMAP, initializing a `gaussian_splat` object ready for [gs_fit()].
#'
#' @param colmap_dir Directory containing the COLMAP sparse reconstruction (e.g. `colmap_workspace/sparse/0`).
#' @param initial_scale Initial 3D splat scale. If `NULL`, estimated automatically via spatial k-NN.
#' @param initial_opacity Initial splat opacity. Defaults to `0.99`.
#' @return A named list containing `cameras` (calibrated camera list), `splat` (`gaussian_splat` object),
#'   and `colmap_dir`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#' col_data <- gs_import_colmap("my_photos/colmap_workspace/sparse/0")
gs_import_colmap <- function(colmap_dir, initial_scale = NULL, initial_opacity = 0.99) {
  fused_file_override <- NULL
  if (inherits(colmap_dir, "colmap_sfm")) {
    if (!is.null(colmap_dir$fused_ply) && file.exists(colmap_dir$fused_ply)) {
      fused_file_override <- colmap_dir$fused_ply
    }
    colmap_dir <- colmap_dir$recon_dir
  }

  if (!dir.exists(colmap_dir)) {
    sub_zero <- file.path(colmap_dir, "sparse", "0")
    if (dir.exists(sub_zero)) {
      colmap_dir <- sub_zero
    } else {
      stop("COLMAP reconstruction directory not found: ", colmap_dir)
    }
  }

  cameras <- gs_import_cameras_colmap(colmap_dir)

  # Point cloud: check for dense fused.ply first, then fallback to sparse points3D.ply
  fused_file <- if (!is.null(fused_file_override)) fused_file_override else file.path(colmap_dir, "fused.ply")
  if (!file.exists(fused_file)) {
    fused_file <- file.path(dirname(dirname(colmap_dir)), "dense", "fused.ply")
  }
  if (!file.exists(fused_file)) {
    fused_file <- file.path(dirname(colmap_dir), "dense", "fused.ply")
  }

  ply_file <- if (file.exists(fused_file)) fused_file else file.path(colmap_dir, "points3D.ply")
  txt_file <- file.path(colmap_dir, "points3D.txt")

  splat <- NULL
  if (file.exists(ply_file)) {
    cloud_type <- if (basename(ply_file) == "fused.ply") "dense MVS" else "sparse tie-point"
    cli::cli_alert_info("Importing {cloud_type} point cloud from {.path {ply_file}}...")
    splat <- gs_import_metashape(ply_file, initial_scale = initial_scale, initial_opacity = initial_opacity)
  } else if (file.exists(txt_file)) {
    cli::cli_alert_info("Parsing 3D tie points from {.path {txt_file}}...")
    lines <- readLines(txt_file, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]

    if (length(lines) > 0) {
      coords_list <- vector("list", length(lines))
      colors_list <- vector("list", length(lines))

      for (i in seq_along(lines)) {
        toks <- strsplit(lines[i], "\\s+")[[1]]
        if (length(toks) >= 7) {
          coords_list[[i]] <- as.numeric(toks[2:4])
          colors_list[[i]] <- as.numeric(toks[5:7]) / 255.0
        }
      }

      valid_idx <- which(!vapply(coords_list, is.null, logical(1)))
      if (length(valid_idx) > 0) {
        pts_mat <- do.call(rbind, coords_list[valid_idx])
        cols_mat <- do.call(rbind, colors_list[valid_idx])
        splat <- gs_create(points = pts_mat, colors = cols_mat,
                           initial_scale = initial_scale, initial_opacity = initial_opacity)
      }
    }
  }

  if (is.null(splat)) {
    cli::cli_alert_warning("No points3D.ply or points3D.txt found. Splats will be initialized from camera trajectory during gs_fit().")
  }

  return(list(
    cameras = cameras,
    splat = splat,
    colmap_dir = colmap_dir
  ))
}

#' @title End-to-End Pipeline: From Raw Photos to 3D Gaussian Splatting
#' @name gs_from_photos
#' @description
#' Executes the complete end-to-end 3D reconstruction workflow:
#' 1. Automatically detects / downloads COLMAP with GPU CUDA support.
#' 2. Executes Structure-from-Motion (SfM) to solve camera poses and sparse tie points.
#' 3. Imports poses and points directly into high-performance C++ structures.
#' 4. Optimizes 3D Gaussian Splats ([gs_fit()]) via differentiable rasterization, photometric loss, and adaptive densification.
#' 5. Optionally opens the interactive WebGL 3D viewer in Google Chrome ([gs_viewer()]).
#'
#' @param images Path to a directory containing photos (.jpg, .png) or character vector of image paths.
#' @param output_dir Destination directory for COLMAP workspace. If `NULL`, saves in `colmap_workspace` inside `images`.
#' @param camera_model COLMAP camera projection model. Defaults to `"OPENCV"`.
#' @param single_camera Logical. If `TRUE` (default), all photos share a single calibrated camera intrinsic.
#' @param matching Feature matching strategy: `"exhaustive"` (default) or `"sequential"`.
#' @param multiple_models Logical. If `FALSE` (default), forces reconstruction of a single unified 360° model.
#' @param guided_matching Logical. If `TRUE` (default), recovers weak visual correspondences.
#' @param overwrite Logical. If `TRUE`, cleans previous workspace reconstructions to start fresh. Defaults to `FALSE`.
#' @param sequential_overlap Integer. Overlap window when `matching = "sequential"`. Defaults to `15`.
#' @param max_num_features Maximum number of visual features per image. Defaults to `8192`.
#' @param iterations Total number of 3D Gaussian Splatting optimization steps. Defaults to `1000`.
#' @param downsample Image downsample factor for training acceleration (e.g. `2`, `4`, or `8`). Defaults to `4`.
#' @param densify_interval Iteration frequency for adaptive densification (split/clone/prune). Defaults to `100`.
#' @param lr_pos Learning rate for 3D Gaussian positions. Defaults to `0.00016`.
#' @param lr_scale Learning rate for 3D Gaussian scales. Defaults to `0.005`.
#' @param lr_opac Learning rate for Gaussian opacities. Defaults to `0.05`.
#' @param viewer Logical. If `TRUE` (default), automatically opens the interactive WebGL 3D viewer upon completion.
#' @param show_colmap_output Logical. If `TRUE`, streams raw COLMAP terminal output in real time to the R console. Defaults to `FALSE`.
#' @param verbose Logical. If `TRUE` (default), prints detailed step-by-step progress.
#' @param ... Additional arguments passed to [colmap_sfm()], [gs_fit()], or [gs_viewer()].
#' @return A trained `gaussian_splat` object.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Full end-to-end pipeline in a single function call!
#' splat <- gs_from_photos("frames/", matching = "exhaustive", iterations = 1000)
gs_from_photos <- function(images,
                           output_dir = NULL,
                           camera_model = "OPENCV",
                           single_camera = TRUE,
                           matching = c("exhaustive", "sequential"),
                           multiple_models = FALSE,
                           guided_matching = TRUE,
                           overwrite = FALSE,
                           sequential_overlap = 15,
                           max_num_features = 8192,
                           iterations = 1000,
                           downsample = 4,
                           densify_interval = 100,
                           lr_pos = 0.00016,
                           lr_scale = 0.005,
                           lr_opac = 0.05,
                           viewer = TRUE,
                           show_colmap_output = FALSE,
                           verbose = TRUE,
                           ...) {
  matching <- match.arg(matching)

  if (isTRUE(verbose)) {
    cli::cli_rule("Integrated 3D Gaussian Splatting Pipeline (pliman)")
  }

  dots <- list(...)
  colmap_formal_names <- names(formals(colmap_sfm))
  colmap_extra <- dots[names(dots) %in% colmap_formal_names]

  # 1. Run Structure-from-Motion with COLMAP
  sfm_args <- c(
    list(
      images = images,
      output_dir = output_dir,
      camera_model = camera_model,
      single_camera = single_camera,
      matching = matching,
      multiple_models = multiple_models,
      guided_matching = guided_matching,
      overwrite = overwrite,
      sequential_overlap = sequential_overlap,
      max_num_features = max_num_features,
      show_colmap_output = show_colmap_output,
      use_gpu = TRUE,
      verbose = verbose
    ),
    colmap_extra
  )
  sfm_res <- do.call(colmap_sfm, sfm_args)

  # 2. Import calibrated cameras and sparse points
  colmap_data <- gs_import_colmap(sfm_res$recon_dir)

  # 3. Train 3D Gaussian Splats
  fit_formal_names <- names(formals(gs_fit))
  fit_extra <- dots[names(dots) %in% fit_formal_names]

  fit_args <- c(
    list(
      images = images,
      cameras = colmap_data$cameras,
      splat = colmap_data$splat,
      iterations = iterations,
      downsample = downsample,
      densify_interval = densify_interval,
      lr_pos = lr_pos,
      lr_scale = lr_scale,
      lr_opac = lr_opac,
      verbose = verbose
    ),
    fit_extra
  )
  splat_trained <- do.call(gs_fit, fit_args)

  attr(splat_trained, "colmap_data") <- colmap_data
  attr(splat_trained, "sfm") <- sfm_res
  attr(splat_trained, "cameras") <- colmap_data$cameras
  attr(splat_trained, "images") <- if (is.character(images) && length(images) == 1 && dir.exists(images)) {
    normalizePath(images, winslash = "/")
  } else {
    images
  }
  splat_trained$cameras <- colmap_data$cameras
  splat_trained$images <- attr(splat_trained, "images")

  # 4. Open interactive WebGL viewer if requested
  if (isTRUE(viewer)) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Launching interactive WebGL 3D viewer in browser...")
    }
    gs_viewer(splat_trained, cameras = colmap_data$cameras, images = images)
  }

  return(invisible(splat_trained))
}
