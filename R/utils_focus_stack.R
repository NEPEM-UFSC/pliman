#' Multi-Focus Image Alignment, Registration and Fusion (MFIF)
#'
#' @description
#' High-performance native tools for focus stacking in digital imaging and microscopy:
#' * `align_stack()`: Sub-pixel image alignment and registration implemented natively
#'   in C++ using multi-scale Gaussian pyramid Enhanced Correlation Coefficient (ECC)
#'   maximization (affine motion or translation) to eliminate focus breathing, translation,
#'   and rotation, with optional automatic border cropping. Zero external Python/OpenCV dependencies.
#' * `fuse_multifocus()`: Multi-Focus Image Fusion (MFIF) supporting:
#'   1. Deep learning ONNX models (e.g., IFCNN) with either hierarchical tournament tree
#'      (`strategy = "pyramid"`) or linear (`strategy = "sequential"`) fusion.
#'   2. Ultra-fast native C++ Modified Laplacian Extended Depth of Field (`method = "laplacian"`),
#'      ideal for massive stacks of 50-500 images in < 1 second.
#' * `process_focus_stack()`: Orchestrates the complete end-to-end focus stacking pipeline,
#'   taking raw images, performing alignment, keyframe subsampling (`max_frames`), and fusion.
#'
#' @param images A list of `Image` objects, an array of images, or a character vector
#'   of file paths to images in the stack.
#' @param reference_index Integer specifying the index of the base reference image
#'   for alignment. If `NULL` (default), the central image of the stack
#'   (`ceiling(length(images) / 2)`) is selected as the optical reference.
#' @param method Character specifying the method:
#'   * In `align_stack()`: Alignment model (`"ecc"` for 2D Affine motion or `"translation"` for 2D shifts).
#'   * In `fuse_multifocus()`: Fusion algorithm (`"ifcnn"` for deep learning or `"laplacian"` for native C++ EDF).
#' @param align_method Alignment method for `process_focus_stack()`. One of `"ecc"` or `"translation"`.
#' @param fusion_method Multi-focus fusion method for `process_focus_stack()`. One of `"ifcnn"` (deep learning)
#'   or `"laplacian"` (native C++ Modified Laplacian Extended Depth of Field).
#' @param strategy Multi-image fusion strategy for deep learning (`"pyramid"` or `"sequential"`).
#'   Defaults to `"pyramid"`. Hierarchical tournament tree reduction pairs images at each level,
#'   bounding the number of sequential convolution passes per pixel to \eqn{\lceil \log_2(N) \rceil}
#'   and preventing contrast loss and blur over large image stacks.
#' @param blend Character specifying blending mode for Laplacian fusion: `"soft"` (default, weighted
#'   softmax power blending) or `"hard"` (maximum sharpness pixel selection).
#' @param radius Integer specifying the local window radius for Laplacian sharpness aggregation (default: 3).
#' @param power Numeric exponent for soft Laplacian blending weights (default: 6.0).
#' @param max_frames Optional integer. If provided and the image stack exceeds this number of frames,
#'   automatically subsamples the stack to `max_frames` evenly spaced focal planes.
#' @param align Logical. Should the image stack be aligned before fusion? Defaults to `TRUE`.
#' @param crop Logical. Should black/empty borders introduced by affine warping
#'   be automatically cropped out? Defaults to `TRUE`.
#' @param levels Integer specifying number of pyramid levels for coarse-to-fine ECC alignment (default: 3).
#' @param max_iter Integer specifying maximum iterations per pyramid level for ECC convergence (default: 50).
#' @param eps Numeric threshold for ECC convergence (default: 1e-4).
#' @param model Character specifying the fusion ONNX model. Defaults to `"ifcnn"`.
#'   Can also be a direct path to a custom ONNX model file.
#' @param engine Execution engine: `"cpu"` (default) or `"gpu"` (DirectML acceleration).
#' @param threads Integer number of threads for ONNX inference (0 for auto).
#' @param device_id GPU device ID (-1 for default best detected GPU).
#' @param return_aligned Logical. If `TRUE`, returns a named list with both `fused`
#'   and `aligned` images. If `FALSE` (default), returns only the fused `Image`.
#' @param verbose Logical. If `TRUE` (default), displays progress messages.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return
#' * `align_stack()`: A list of aligned `Image` objects with identical spatial dimensions.
#' * `fuse_multifocus()`: An all-in-focus `Image` object.
#' * `process_focus_stack()`: The fused `Image` object (or a list containing `fused` and `aligned` if `return_aligned = TRUE`).
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @references
#' Zhang, Y., Liu, Y., Sun, P., Yan, H., Zhao, X., & Zhang, L. (2020).
#' IFCNN: A general image fusion framework based on convolutional neural network.
#' Information Fusion, 54, 99-118.
#'
#' Evangelidis, G. D., & Psarakis, E. Z. (2008). Parametric image alignment using
#' enhanced correlation coefficient maximization. IEEE Transactions on Pattern
#' Analysis and Machine Intelligence, 30(10), 1858-1865.
#'
#' @name focus_stack
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Load or read a stack of multi-focus images
#' imgs <- list(
#'   image_read("focus_frame1.jpg"),
#'   image_read("focus_frame2.jpg"),
#'   image_read("focus_frame3.jpg")
#' )
#'
#' # Fast automated pipeline with native C++ Modified Laplacian (zero-AI, ultra-fast)
#' fused_lap <- process_focus_stack(imgs, fusion_method = "laplacian", blend = "soft")
#' plot(fused_lap)
#'
#' # Deep learning IFCNN with hierarchical tournament pyramid reduction
#' fused_dl <- process_focus_stack(imgs, fusion_method = "ifcnn", strategy = "pyramid")
#' plot(fused_dl)
#'
#' # Handling large stacks (e.g. 200 frames from a focus rack video)
#' # Automatically extract 25 key focal planes for alignment and fusion:
#' fused_rack <- process_focus_stack(imgs, max_frames = 25, fusion_method = "laplacian")
#' }
NULL

# -----------------------------------------------------------------------------
# Internal Image Stack Coercion and I/O Helpers
# -----------------------------------------------------------------------------

.coerce_image_stack <- function(images) {
  if (is.character(images)) {
    valid_paths <- vapply(images, file.exists, logical(1))
    if (!all(valid_paths)) {
      missing_files <- images[!valid_paths]
      cli::cli_abort("The following image file(s) do not exist: {.val {missing_files}}")
    }
    stack <- lapply(images, function(p) {
      img <- image_import(p)
      as_image(img, colormode = "Color")
    })
    return(stack)
  }

  if (is.list(images)) {
    if (length(images) == 0L) {
      cli::cli_abort("The input image stack is empty.")
    }
    stack <- lapply(seq_along(images), function(i) {
      im <- images[[i]]
      if (inherits(im, "Image") || inherits(im, "image")) {
        return(as_image(im, colormode = "Color"))
      }
      if (is.array(im) || is.matrix(im)) {
        return(as_image(im, colormode = "Color"))
      }
      cli::cli_abort("Element {i} in the list cannot be coerced to an Image object.")
    })
    return(stack)
  }

  if (is.array(images)) {
    dims <- dim(images)
    if (length(dims) == 4L) {
      # 4D array [W, H, C, N] or [H, W, C, N]
      n_imgs <- dims[4]
      stack <- lapply(seq_len(n_imgs), function(i) {
        as_image(images[, , , i], colormode = "Color")
      })
      return(stack)
    }
    if (length(dims) == 3L) {
      return(list(as_image(images, colormode = "Color")))
    }
  }

  cli::cli_abort("Unsupported input format for {.arg images}. Must be a list of images, character paths, or a 4D array.")
}

.image_to_nchw_vec <- function(img) {
  u <- unclass(img)
  d <- dim(u)
  if (length(d) == 2L) {
    w <- d[1]
    h <- d[2]
    c <- 1L
  } else {
    w <- d[1]
    h <- d[2]
    c <- d[3]
  }

  is_raw <- is.raw(u)

  if (c == 1L) {
    # Duplicate single channel to 3 channels (RGB)
    c_data <- if (is_raw) {
      as.numeric(u) / 255.0
    } else {
      v <- as.numeric(u)
      if (max(v, na.rm = TRUE) > 1.5) v <- v / 255.0
      pmax(0.0, pmin(1.0, v))
    }
    flat_nchw <- c(c_data, c_data, c_data)
  } else {
    extract_ch <- function(ch) {
      slice <- u[, , ch]
      if (is_raw) {
        as.numeric(slice) / 255.0
      } else {
        vals <- as.numeric(slice)
        if (max(vals, na.rm = TRUE) > 1.5) vals <- vals / 255.0
        pmax(0.0, pmin(1.0, vals))
      }
    }
    r_ch <- extract_ch(1)
    g_ch <- extract_ch(2)
    b_ch <- extract_ch(3)
    flat_nchw <- c(r_ch, g_ch, b_ch)
  }

  list(vec = flat_nchw, w = as.integer(w), h = as.integer(h))
}

.nchw_vec_to_image <- function(vec, w, h, storage = "double") {
  plane <- w * h
  arr <- array(0.0, dim = c(w, h, 3))
  arr[, , 1] <- matrix(vec[1:plane], nrow = w, ncol = h)
  arr[, , 2] <- matrix(vec[(plane + 1):(2 * plane)], nrow = w, ncol = h)
  arr[, , 3] <- matrix(vec[(2 * plane + 1):(3 * plane)], nrow = w, ncol = h)

  if (storage == "raw") {
    raw_arr <- as.raw(pmin(255, pmax(0, round(arr * 255))))
    dim(raw_arr) <- c(w, h, 3)
    return(as_image(raw_arr, colormode = "Color", storage = "raw"))
  }

  as_image(arr, colormode = "Color", storage = "double")
}

# -----------------------------------------------------------------------------
# 1. Alignment Module: align_stack
# -----------------------------------------------------------------------------

#' @rdname focus_stack
#' @export
align_stack <- function(images,
                        reference_index = NULL,
                        method = c("ecc", "translation", "orb"),
                        crop = TRUE,
                        levels = 3,
                        max_iter = 50,
                        eps = 1e-4,
                        max_frames = NULL,
                        verbose = TRUE,
                        ...) {
  method <- match.arg(method)
  if (method == "orb") {
    # Map legacy "orb" to native sub-pixel "ecc"
    method <- "ecc"
  }
  stack <- .coerce_image_stack(images)
  n <- length(stack)

  if (n <= 1L) {
    if (isTRUE(verbose)) cli::cli_alert_info("Stack contains only {n} image. Returning unmodified.")
    return(stack)
  }

  # Subsample frames if max_frames is specified
  if (!is.null(max_frames) && is.numeric(max_frames) && max_frames > 0L && n > max_frames) {
    max_f <- as.integer(max_frames[1])
    sample_idx <- sort(unique(round(seq(1, n, length.out = max_f))))
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Subsampling stack from {n} frames to {length(sample_idx)} evenly spaced focal planes.")
    }
    stack <- stack[sample_idx]
    n <- length(stack)
    if (!is.null(reference_index)) {
      reference_index <- as.integer(ceiling(n / 2))
    }
  }

  # Default reference image is the central frame
  if (is.null(reference_index)) {
    reference_index <- as.integer(ceiling(n / 2))
  } else {
    reference_index <- as.integer(reference_index[1])
    if (reference_index < 1L || reference_index > n) {
      cli::cli_abort("{.arg reference_index} ({reference_index}) must be between 1 and {n}.")
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Aligning stack of {n} images natively in C++ with {toupper(method)}...")
  }

  res <- align_stack_cpp(
    images = stack,
    ref_idx = as.integer(reference_index - 1L),
    method = method,
    crop = isTRUE(crop),
    levels = as.integer(levels),
    max_iter = as.integer(max_iter),
    eps = as.numeric(eps),
    verbose = isTRUE(verbose)
  )

  aligned_stack <- res$aligned
  statuses <- as.character(res$statuses)

  # Report warnings if any images experienced convergence or dimension issues
  for (i in seq_along(statuses)) {
    st <- statuses[i]
    if (st %in% c("fallback_identity", "dimension_mismatch")) {
      cli::cli_alert_warning("Image {i}: Alignment failed ({st}); preserving original frame.")
    }
  }

  if (isTRUE(verbose)) {
    orig_s <- dim(stack[[reference_index]])
    out_s <- dim(aligned_stack[[1]])
    cli::cli_alert_success(
      "Stack aligned natively in C++ (reference #{reference_index}). Spatial size: {orig_s[1]}x{orig_s[2]} -> {out_s[1]}x{out_s[2]}."
    )
  }

  attr(aligned_stack, "matrices") <- res$matrices
  attr(aligned_stack, "crop_box") <- res$crop_box
  attr(aligned_stack, "statuses") <- res$statuses

  return(aligned_stack)
}

#' @rdname focus_stack
#' @export
image_align_stack <- align_stack

# -----------------------------------------------------------------------------
# 2. Multi-Focus Fusion Module: fuse_multifocus
# -----------------------------------------------------------------------------

#' @rdname focus_stack
#' @export
fuse_multifocus <- function(images,
                            method = c("ifcnn", "laplacian"),
                            strategy = c("pyramid", "sequential"),
                            blend = c("soft", "hard"),
                            radius = 3,
                            power = 6.0,
                            max_frames = NULL,
                            model = "ifcnn",
                            engine = c("gpu", "cpu"),
                            threads = 0L,
                            device_id = -1L,
                            verbose = TRUE,
                            ...) {
  dots <- list(...)
  if (!missing(model) && missing(method) && model %in% c("laplacian", "ifcnn")) {
    method <- model
  }
  method <- match.arg(method)
  strategy <- match.arg(strategy)
  blend <- match.arg(blend)
  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")

  stack <- .coerce_image_stack(images)
  n <- length(stack)

  if (n == 0L) {
    cli::cli_abort("The input image stack is empty.")
  }
  if (n == 1L) {
    if (isTRUE(verbose)) cli::cli_alert_info("Single image provided. Returning as is.")
    return(stack[[1]])
  }

  # Subsample if max_frames is specified
  if (!is.null(max_frames) && is.numeric(max_frames) && max_frames > 0L && n > max_frames) {
    max_f <- as.integer(max_frames[1])
    sample_idx <- sort(unique(round(seq(1, n, length.out = max_f))))
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Subsampling stack from {n} frames to {length(sample_idx)} evenly spaced focal planes.")
    }
    stack <- stack[sample_idx]
    n <- length(stack)
  }

  # Verify consistent dimensions across stack
  first_dim <- dim(stack[[1]])
  w <- first_dim[1]
  h <- first_dim[2]
  storage_type <- if (is.raw(stack[[1]])) "raw" else "double"

  for (i in 2:n) {
    cur_dim <- dim(stack[[i]])
    if (cur_dim[1] != w || cur_dim[2] != h) {
      cli::cli_abort(
        "Image dimensions mismatch in stack: Image 1 is {w}x{h}, but Image {i} is {cur_dim[1]}x{cur_dim[2]}. Run {.fn align_stack} first."
      )
    }
  }

  # Analytical Modified Laplacian Extended Depth of Field (EDF)
  if (method == "laplacian") {
    if (isTRUE(verbose)) {
      cli::cli_progress_step(
        msg = "Fusing stack of {n} images natively in C++ via Modified Laplacian (EDF) [{blend}, r={radius}]...",
        msg_done = "Laplacian multi-focus fusion complete"
      )
    }
    out_img <- fuse_laplacian_cpp(
      images = stack,
      blend = blend,
      radius = as.integer(radius),
      power = as.numeric(power)
    )
    return(out_img)
  }

  # Deep learning ONNX (IFCNN)
  # Resolve ONNX Runtime library
  lib_path <- pliman_onnx_lib_path(engine = engine)
  if (is.null(lib_path) || !file.exists(lib_path)) {
    if (isTRUE(verbose)) cli::cli_alert_info("ONNX Runtime library ({engine}) not found. Downloading...")
    lib_path <- onnx_install(engine = engine)
  }

  # Resolve ONNX model path
  model_name <- if (is.character(model)) .resolve_model_name(model[1]) else "ifcnn"
  model_file <- if (file.exists(model[1])) {
    normalizePath(model[1], winslash = "/")
  } else {
    file_candidate <- file.path(pliman_model_dir(), paste0(model_name, ".onnx"))
    if (file.exists(file_candidate)) {
      file_candidate
    } else {
      pliman_download_model(model = model_name)
    }
  }

  if (!file.exists(model_file)) {
    cli::cli_abort("Could not find or download the ONNX model file: {.val {model_name}}.")
  }

  # Preprocess all images into flat NCHW vectors
  tensors <- lapply(seq_len(n), function(i) {
    .image_to_nchw_vec(stack[[i]])$vec
  })

  norm_model <- normalizePath(model_file, winslash = "/", mustWork = FALSE)
  norm_lib <- normalizePath(lib_path, winslash = "/", mustWork = FALSE)

  if (isTRUE(verbose)) {
    cli::cli_progress_step(
      msg = "Fusing stack of {n} images with {model_name} [{engine} | strategy: {strategy}]...",
      msg_done = "Multi-focus fusion complete"
    )
  }

  if (strategy == "pyramid") {
    # Hierarchical binary tournament tree reduction:
    # Pairs up images at each round, reducing passes from N-1 to ceil(log2(N))
    cur_tensors <- tensors
    while (length(cur_tensors) > 1L) {
      n_cur <- length(cur_tensors)
      next_tensors <- list()
      pairs_to_run <- seq(1L, n_cur, by = 2L)
      for (i in pairs_to_run) {
        if (i == n_cur) {
          # Odd leftover image carries forward to next round
          next_tensors[[length(next_tensors) + 1L]] <- cur_tensors[[i]]
        } else {
          fused_pair <- run_dual_image_inference_cpp(
            img1_vec = cur_tensors[[i]],
            img2_vec = cur_tensors[[i + 1L]],
            in_w = as.integer(w),
            in_h = as.integer(h),
            model_path = norm_model,
            lib_path = norm_lib,
            num_threads = as.integer(threads),
            use_gpu = use_gpu,
            device_id = as.integer(device_id)
          )
          next_tensors[[length(next_tensors) + 1L]] <- fused_pair
        }
      }
      cur_tensors <- next_tensors
    }
    fused_tensor <- cur_tensors[[1L]]
  } else {
    # Sequential / linear pairwise fusion:
    fused_tensor <- tensors[[1L]]
    for (k in 2:n) {
      fused_tensor <- run_dual_image_inference_cpp(
        img1_vec = fused_tensor,
        img2_vec = tensors[[k]],
        in_w = as.integer(w),
        in_h = as.integer(h),
        model_path = norm_model,
        lib_path = norm_lib,
        num_threads = as.integer(threads),
        use_gpu = use_gpu,
        device_id = as.integer(device_id)
      )
    }
  }

  # Convert back to standard pliman Image
  out_img <- .nchw_vec_to_image(fused_tensor, w = w, h = h, storage = storage_type)
  return(out_img)
}

#' @rdname focus_stack
#' @export
image_fuse_multifocus <- fuse_multifocus

# -----------------------------------------------------------------------------
# 3. High-Level Orchestrator: process_focus_stack
# -----------------------------------------------------------------------------

#' @rdname focus_stack
#' @export
process_focus_stack <- function(images,
                                align = TRUE,
                                align_method = c("ecc", "translation", "orb"),
                                reference_index = NULL,
                                crop = TRUE,
                                levels = 3,
                                max_iter = 50,
                                eps = 1e-4,
                                max_frames = NULL,
                                method = NULL,
                                fusion_method = c("ifcnn", "laplacian"),
                                strategy = c("pyramid", "sequential"),
                                blend = c("soft", "hard"),
                                radius = 3,
                                power = 6.0,
                                model = "ifcnn",
                                engine = c("gpu", "cpu"),
                                threads = 0L,
                                device_id = -1L,
                                return_aligned = FALSE,
                                verbose = TRUE,
                                ...) {
  align_method <- match.arg(align_method)
  engine <- match.arg(engine)

  dots <- list(...)
  if (!is.null(method)) {
    if (method %in% c("laplacian", "ifcnn")) {
      fusion_method <- method
    } else if (method %in% c("ecc", "translation", "orb")) {
      align_method <- method
    }
  } else if (!is.null(dots$method)) {
    m <- dots$method
    if (m %in% c("laplacian", "ifcnn")) {
      fusion_method <- m
    } else if (m %in% c("ecc", "translation", "orb")) {
      align_method <- m
    }
  }
  dots$method <- NULL

  fusion_method <- match.arg(fusion_method, c("ifcnn", "laplacian"))
  strategy <- match.arg(strategy, c("pyramid", "sequential"))
  blend <- match.arg(blend, c("soft", "hard"))

  stack <- .coerce_image_stack(images)
  n <- length(stack)

  # Step 0: Subsample if requested (before alignment to save substantial time)
  if (!is.null(max_frames) && is.numeric(max_frames) && max_frames > 0L && n > max_frames) {
    max_f <- as.integer(max_frames[1])
    sample_idx <- sort(unique(round(seq(1, n, length.out = max_f))))
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Subsampling stack from {n} frames to {length(sample_idx)} evenly spaced focal planes.")
    }
    stack <- stack[sample_idx]
    if (!is.null(reference_index)) {
      reference_index <- as.integer(ceiling(length(stack) / 2))
    }
  }

  # Step 1: Align Stack natively in C++
  if (isTRUE(align)) {
    align_args <- c(
      list(
        images = stack,
        reference_index = reference_index,
        method = align_method,
        crop = crop,
        levels = levels,
        max_iter = max_iter,
        eps = eps,
        max_frames = NULL, # Already subsampled
        verbose = verbose
      ),
      dots
    )
    aligned_images <- do.call(align_stack, align_args)
  } else {
    aligned_images <- stack
  }

  # Step 2: Multi-Focus Fusion
  fuse_args <- c(
    list(
      images = aligned_images,
      method = fusion_method,
      strategy = strategy,
      blend = blend,
      radius = radius,
      power = power,
      max_frames = NULL, # Already subsampled
      model = model,
      engine = engine,
      threads = threads,
      device_id = device_id,
      verbose = verbose
    ),
    dots
  )
  fused_image <- do.call(fuse_multifocus, fuse_args)

  # Step 3: Return Output
  if (isTRUE(return_aligned)) {
    return(list(fused = fused_image, aligned = aligned_images))
  }

  return(fused_image)
}

#' @rdname focus_stack
#' @export
image_focus_stack <- process_focus_stack
