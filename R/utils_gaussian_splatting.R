# ==============================================================================
# 3D GAUSSIAN SPLATTING FOR PLANT RECONSTRUCTION & PHENOTYPING
# pliman: Plant Image Analysis in R
# ==============================================================================

#' @title Initialize a 3D Gaussian Splatting Cloud
#' @name gs_create
#' @description
#' Creates a `gaussian_splat` object from 3D coordinates, point clouds, or bounding geometries.
#' Automatically computes initial isotropic scales based on k-nearest neighbor (k-NN) distances,
#' initializes rotation quaternions to identity, opacities to initial transmittance,
#' and Spherical Harmonics degree 0 coefficients from RGB colors.
#'
#' @param points Numeric matrix or data.frame with 3 columns (`x, y, z`).
#' @param colors Optional matrix or data.frame of RGB colors in `[0, 1]` or `[0, 255]`.
#'   Defaults to uniform green `c(0.2, 0.75, 0.25)` for plant points.
#' @param initial_scale Optional numeric scalar. If `NULL` (default), scales are automatically
#'   computed as the mean distance to the 3 nearest neighbors.
#' @param initial_opacity Numeric scalar between 0 and 1. Initial Gaussian opacity. Defaults to `0.15`.
#' @param sh_degree Integer. Maximum Spherical Harmonics degree (0, 1, or 2). Defaults to `0`.
#' @return An S3 object of class `gaussian_splat`.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Create 500 random points representing a small plant canopy
#' pts <- matrix(rnorm(1500, mean = 0, sd = 0.5), ncol = 3)
#' splat <- gs_create(pts)
#' print(splat)
#' }
gs_create <- function(points,
                      colors = NULL,
                      initial_scale = NULL,
                      initial_opacity = 0.15,
                      sh_degree = 0) {
  if (is.data.frame(points)) {
    points <- as.matrix(points[, 1:3])
  } else if (!is.matrix(points) || ncol(points) < 3) {
    stop("'points' must be a numeric matrix or data.frame with at least 3 columns (x, y, z).")
  }

  n <- nrow(points)
  if (n == 0) stop("'points' cannot be empty.")

  # Extract RGB colors
  if (is.null(colors)) {
    # Default vegetative plant green
    rgb_mat <- matrix(rep(c(0.22, 0.72, 0.28), each = n), nrow = n, ncol = 3)
  } else {
    if (is.data.frame(colors)) colors <- as.matrix(colors)
    if (ncol(colors) < 3) stop("'colors' must have 3 columns (R, G, B).")
    rgb_mat <- colors[, 1:3, drop = FALSE]
    if (max(rgb_mat, na.rm = TRUE) > 1.0) {
      rgb_mat <- rgb_mat / 255.0
    }
  }

  # Initial log-scales: use true spatial density spacing estimate or fixed scalar
  if (is.null(initial_scale)) {
    if (n >= 4) {
      n_sample <- min(2000L, n)
      idx_sample <- if (n > n_sample) sample.int(n, n_sample) else seq_len(n)
      sample_pts <- points[idx_sample, , drop = FALSE]
      knn_out <- cpp_gs_knn(sample_pts, queries = NULL, k = 1L)
      med_knn <- stats::median(knn_out$dists[, 1], na.rm = TRUE)
      if (is.finite(med_knn) && med_knn > 1e-5) {
        init_s <- min(med_knn * 1.2, 0.05)
      } else {
        init_s <- 0.01
      }
    } else {
      init_s <- 0.02
    }
  } else {
    init_s <- max(1e-5, as.numeric(initial_scale))
  }

  log_s <- log(init_s)
  scales_log <- matrix(log_s, nrow = n, ncol = 3)

  # Identity quaternions: (qw = 1, qx = 0, qy = 0, qz = 0)
  quaternions <- matrix(0.0, nrow = n, ncol = 4)
  quaternions[, 1] <- 1.0

  # Opacity logit
  init_opac <- max(0.01, min(0.99, initial_opacity))
  logit_opac <- log(init_opac / (1.0 - init_opac))
  opacities_logit <- rep(logit_opac, n)

  # Spherical Harmonics degree 0: C = 0.5 + 0.28209479 * SH0 -> SH0 = (C - 0.5) / 0.28209479
  sh_c0 <- 0.28209479177387814
  num_sh <- if (sh_degree == 0) 3L else if (sh_degree == 1) 12L else 27L
  sh_coeffs <- matrix(0.0, nrow = n, ncol = num_sh)
  sh_coeffs[, 1] <- (rgb_mat[, 1] - 0.5) / sh_c0
  sh_coeffs[, 2] <- (rgb_mat[, 2] - 0.5) / sh_c0
  sh_coeffs[, 3] <- (rgb_mat[, 3] - 0.5) / sh_c0

  semantic_ids <- rep(0L, n)
  semantic_scores <- rep(0.0, n)

  obj <- list(
    positions = points,
    scales_log = scales_log,
    quaternions = quaternions,
    opacities_logit = opacities_logit,
    sh_coeffs = sh_coeffs,
    semantic_ids = semantic_ids,
    semantic_scores = semantic_scores,
    sh_degree = as.integer(sh_degree),
    num_gaussians = n
  )
  class(obj) <- "gaussian_splat"
  obj
}

#' @title Print Summary of a 3D Gaussian Splat
#' @param x A `gaussian_splat` object.
#' @param ... Ignored.
#' @export
print.gaussian_splat <- function(x, ...) {
  cat("====================================================\n")
  cat("3D Gaussian Splat Cloud [pliman]\n")
  cat("----------------------------------------------------\n")
  cat(sprintf("Number of Gaussians:  %d\n", x$num_gaussians))
  if (isTRUE(x$has_faces) && !is.null(x$num_faces) && x$num_faces > 0) {
    cat(sprintf("3D Polygonal Mesh:    Yes (%s triangular faces)\n", format(x$num_faces, big.mark = ",")))
  }
  cat(sprintf("SH Degree:            %d (%d coefficients)\n", x$sh_degree, ncol(x$sh_coeffs)))

  # Bounding box
  bx_min <- apply(x$positions, 2, min)
  bx_max <- apply(x$positions, 2, max)
  cat(sprintf("Bounding Box X:       [%.3f, %.3f] (dx = %.3f)\n", bx_min[1], bx_max[1], bx_max[1] - bx_min[1]))
  cat(sprintf("Bounding Box Y:       [%.3f, %.3f] (dy = %.3f)\n", bx_min[2], bx_max[2], bx_max[2] - bx_min[2]))
  cat(sprintf("Bounding Box Z:       [%.3f, %.3f] (dz = %.3f)\n", bx_min[3], bx_max[3], bx_max[3] - bx_min[3]))

  # Semantic classes
  tbl <- table(x$semantic_ids)
  if (length(tbl) > 1 || names(tbl)[1] != "0") {
    cat("Semantic Classes:\n")
    for (i in seq_along(tbl)) {
      lbl <- if (names(tbl)[i] == "0") "Unclassified (0)" else if (names(tbl)[i] == "1") "Leaf (1)" else paste0("Class ", names(tbl)[i])
      cat(sprintf("  * %-18s: %d splats\n", lbl, tbl[i]))
    }
  }
  cat("====================================================\n")
  invisible(x)
}

#' @title Construct Camera Extrinsics and Intrinsics
#' @name gs_camera
#' @description
#' Creates a camera specification list with world-to-camera rotation matrix $R$,
#' translation vector $t$, focal lengths, and principal point.
#'
#' @param position Numeric vector of length 3 `c(x, y, z)` specifying the camera position in world space.
#' @param look_at Numeric vector of length 3 `c(x, y, z)` specifying the point the camera is looking at. Defaults to `c(0, 0, 0)`.
#' @param up Numeric vector of length 3 `c(x, y, z)` specifying the up vector. Defaults to `c(0, 0, 1)`.
#' @param fov Numeric scalar. Vertical field of view in degrees. Defaults to `50.0`.
#' @param width Integer. Viewport width in pixels. Defaults to `640`.
#' @param height Integer. Viewport height in pixels. Defaults to `480`.
#' @return A list with camera extrinsics and intrinsics.
#' @export
gs_camera <- function(position = c(0, -2, 1),
                      look_at = c(0, 0, 0),
                      up = c(0, 0, 1),
                      fov = 50.0,
                      width = 640,
                      height = 480) {
  pos <- as.numeric(position)
  target <- as.numeric(look_at)
  up_v <- as.numeric(up)

  # Forward vector: looking from pos to target
  fwd <- target - pos
  norm_fwd <- sqrt(sum(fwd^2))
  if (norm_fwd < 1e-8) fwd <- c(0, 1, 0) else fwd <- fwd / norm_fwd

  # Right vector: fwd x up
  right <- c(
    fwd[2] * up_v[3] - fwd[3] * up_v[2],
    fwd[3] * up_v[1] - fwd[1] * up_v[3],
    fwd[1] * up_v[2] - fwd[2] * up_v[1]
  )
  norm_r <- sqrt(sum(right^2))
  if (norm_r < 1e-8) right <- c(1, 0, 0) else right <- right / norm_r

  # True up vector: right x fwd
  cam_up <- c(
    right[2] * fwd[3] - right[3] * fwd[2],
    right[3] * fwd[1] - right[1] * fwd[3],
    right[1] * fwd[2] - right[2] * fwd[1]
  )

  # Standard camera coordinate system: X = right, Y = down (-cam_up), Z = forward
  # R matrix rows are camera axes
  R <- matrix(0.0, 3, 3)
  R[1, ] <- right
  R[2, ] <- -cam_up
  R[3, ] <- fwd

  # Translation in camera coordinates: t = -R * pos
  t_vec <- as.numeric(-R %*% pos)

  # Intrinsics
  fy <- (height / 2.0) / tan((fov * pi / 180.0) / 2.0)
  fx <- fy # Square pixels
  cx <- width / 2.0
  cy <- height / 2.0

  list(
    R = R,
    t = t_vec,
    fx = fx,
    fy = fy,
    cx = cx,
    cy = cy,
    width = as.integer(width),
    height = as.integer(height),
    position = pos,
    look_at = target
  )
}

#' @title Render 2D Views from 3D Gaussian Splats
#' @name gs_render
#' @description
#' Executes high-performance differentiable tile-based forward rasterization of a `gaussian_splat`
#' object from a specified camera pose, producing an RGB `Image` object and depth map.
#'
#' @param splat A `gaussian_splat` object.
#' @param camera A camera list created with [gs_camera()].
#' @param bg Color of background. Defaults to `"black"` (or `c(0, 0, 0)`).
#' @param return_depth Logical. If `TRUE`, attaches depth map as attribute. Defaults to `TRUE`.
#' @param plot Logical. If `TRUE`, displays the rendered image using `plot()`. Defaults to `TRUE`.
#' @return An `Image` object with attributes `depth` and `alpha`.
#' @export
gs_render <- function(splat,
                      camera = gs_camera(),
                      bg = "black",
                      return_depth = TRUE,
                      plot = TRUE) {
  if (!inherits(splat, "gaussian_splat")) {
    stop("'splat' must be an object of class 'gaussian_splat'.")
  }

  bg_rgb <- if (is.character(bg)) {
    grDevices::col2rgb(bg)[, 1] / 255.0
  } else {
    rep(as.numeric(bg), length.out = 3)
  }

  res <- cpp_gs_render(
    camera = camera,
    gaussians = splat,
    sh_degree = splat$sh_degree,
    bg_r = bg_rgb[1],
    bg_g = bg_rgb[2],
    bg_b = bg_rgb[3],
    return_depth = return_depth
  )

  img <- as_image(res$image)
  attr(img, "depth") <- res$depth
  attr(img, "alpha") <- res$alpha
  attr(img, "splats_count") <- res$splats_count
  attr(img, "cache") <- res$cache

  if (isTRUE(plot)) {
    plot(img)
  }

  invisible(img)
}

#' @title Segment Leaves and Organs in 3D using SAM
#' @name gs_segment_sam
#' @description
#' Bridges 2D foundational segmentation models (SAM 2.1 / SAM 3.1) with 3D Gaussian Splatting.
#' Projects 2D leaf or organ masks from photos onto 3D Gaussians via ray-splat correspondence,
#' updating semantic class IDs and confidence scores in the 3D cloud.
#'
#' @param splat A `gaussian_splat` object.
#' @param camera Camera list corresponding to the viewpoint of the photo.
#' @param mask 2D binary matrix or an `Image` object returned by [image_segment_dl()] where `1` represents the target leaf.
#' @param class_id Integer. Semantic class ID to assign (e.g. `1` for leaf, `2` for stem). Defaults to `1`.
#' @param threshold Double. Minimum confidence score to accept semantic label. Defaults to `0.05`.
#' @return An updated `gaussian_splat` object with assigned 3D semantic classes.
#' @export
gs_segment_sam <- function(splat,
                           camera,
                           mask,
                           class_id = 1,
                           threshold = 0.05) {
  if (!inherits(splat, "gaussian_splat")) {
    stop("'splat' must be an object of class 'gaussian_splat'.")
  }

  # Convert mask to binary IntegerMatrix
  mat_mask <- if (inherits(mask, c("Image", "image"))) {
    m <- image_data(mask)
    if (length(dim(m)) == 3) m <- m[, , 1]
    matrix(as.integer(m > 0.5), nrow = nrow(m), ncol = ncol(m))
  } else if (is.matrix(mask)) {
    matrix(as.integer(mask > 0), nrow = nrow(mask), ncol = ncol(mask))
  } else {
    stop("'mask' must be a binary matrix or 'Image' object.")
  }

  # Orient mask dimensions if passed as (height, width)
  if (!is.null(camera$width) && !is.null(camera$height)) {
    if (nrow(mat_mask) == camera$height && ncol(mat_mask) == camera$width) {
      mat_mask <- t(mat_mask)
    }
  }

  res <- cpp_gs_segment_sam(
    camera = camera,
    gaussians = splat,
    mask = mat_mask,
    class_id = as.integer(class_id),
    score_threshold = as.numeric(threshold)
  )

  splat$semantic_ids <- res$semantic_ids
  splat$semantic_scores <- res$semantic_scores

  num_labeled <- sum(splat$semantic_ids == class_id)
  cli::cli_alert_success("Assigned semantic class {class_id} to {num_labeled} 3D Gaussian splat(s).")

  splat
}

#' @title Extract 3D Plant Biometrics from Gaussian Splats
#' @name gs_measure
#' @description
#' Computes real 3D phenotypic traits from 3D Gaussian Splats, including:
#' * Total 3D surface area (e.g. true 3D leaf area index)
#' * Plant height, width, and depth
#' * Bounding box dimensions
#' * Leaf insertion angles relative to vertical
#'
#' @param splat A `gaussian_splat` object.
#' @param target_class Integer. Class ID to measure (e.g. `1` for leaf only, `-1` for all splats). Defaults to `-1`.
#' @return A data frame with 3D biometric traits.
#' @export
gs_measure <- function(splat, target_class = -1) {
  if (!inherits(splat, "gaussian_splat")) {
    stop("'splat' must be an object of class 'gaussian_splat'.")
  }

  res <- cpp_gs_measure(
    positions = splat$positions,
    scales_log = splat$scales_log,
    quaternions = splat$quaternions,
    semantic_ids = splat$semantic_ids,
    target_class = as.integer(target_class)
  )

  angles <- res$leaf_angles
  mean_angle <- if (length(angles) > 0) mean(angles) else NA_real_
  sd_angle <- if (length(angles) > 1) stats::sd(angles) else NA_real_

  df <- data.frame(
    num_splats = res$num_splats,
    surface_area = round(res$surface_area, 4),
    height = round(res$height, 4),
    width = round(res$width, 4),
    depth = round(res$depth, 4),
    mean_leaf_angle = round(mean_angle, 2),
    sd_leaf_angle = round(sd_angle, 2),
    stringsAsFactors = FALSE
  )
  attr(df, "bbox_min") <- res$bbox_min
  attr(df, "bbox_max") <- res$bbox_max
  attr(df, "all_leaf_angles") <- angles
  df
}

#' @title Export 3D Gaussian Splats to Standard PLY Format
#' @name gs_export_ply
#' @description
#' Exports a `gaussian_splat` object to the standard binary PLY format used by
#' 3D Gaussian Splatting viewers (SuperSplat, WebGL, Unity, Nerfstudio).
#'
#' @param splat A `gaussian_splat` object.
#' @param file Path to the output `.ply` file.
#' @return Invisibly returns the file path.
#' @export
gs_export_ply <- function(splat, file) {
  if (!inherits(splat, "gaussian_splat")) {
    stop("'splat' must be an object of class 'gaussian_splat'.")
  }

  con <- file(file, "wb")
  on.exit(close(con))

  n <- splat$num_gaussians

  # PLY Header
  header <- paste0(
    "ply\n",
    "format binary_little_endian 1.0\n",
    "element vertex ", n, "\n",
    "property float x\n",
    "property float y\n",
    "property float z\n",
    "property float nx\n",
    "property float ny\n",
    "property float nz\n",
    "property float f_dc_0\n",
    "property float f_dc_1\n",
    "property float f_dc_2\n",
    "property float opacity\n",
    "property float scale_0\n",
    "property float scale_1\n",
    "property float scale_2\n",
    "property float rot_0\n",
    "property float rot_1\n",
    "property float rot_2\n",
    "property float rot_3\n",
    "property int semantic_id\n",
    "end_header\n"
  )
  writeChar(header, con, eos = NULL)

  # Write binary interleaved attributes
  # Pre-pack vectors
  px <- as.numeric(splat$positions[, 1])
  py <- as.numeric(splat$positions[, 2])
  pz <- as.numeric(splat$positions[, 3])
  zeros <- rep(0.0, n)

  sh0 <- as.numeric(splat$sh_coeffs[, 1])
  sh1 <- as.numeric(splat$sh_coeffs[, 2])
  sh2 <- as.numeric(splat$sh_coeffs[, 3])

  opac <- as.numeric(splat$opacities_logit)

  s0 <- as.numeric(splat$scales_log[, 1])
  s1 <- as.numeric(splat$scales_log[, 2])
  s2 <- as.numeric(splat$scales_log[, 3])

  q0 <- as.numeric(splat$quaternions[, 1])
  q1 <- as.numeric(splat$quaternions[, 2])
  q2 <- as.numeric(splat$quaternions[, 3])
  q3 <- as.numeric(splat$quaternions[, 4])

  sem <- as.integer(splat$semantic_ids)

  for (i in seq_len(n)) {
    writeBin(as.numeric(c(px[i], py[i], pz[i], zeros[i], zeros[i], zeros[i],
                          sh0[i], sh1[i], sh2[i], opac[i],
                          s0[i], s1[i], s2[i],
                          q0[i], q1[i], q2[i], q3[i])), con, size = 4)
    writeBin(sem[i], con, size = 4)
  }

  cli::cli_alert_success("Exported {n} Gaussian splat(s) to {.path {file}}.")
  invisible(file)
}

#' @title Adaptive Densification and Pruning for 3D Gaussian Splats
#' @name gs_densify
#' @description
#' Implements the split, clone, and pruning heuristics for 3D Gaussian Splatting:
#' * **Split:** Gaussians with large view-space gradients and large scale are divided into two smaller Gaussians.
#' * **Clone:** Gaussians with large view-space gradients and small scale are duplicated to cover missing geometry.
#' * **Prune:** Gaussians with opacity below `min_opacity` or excessive scale are removed.
#'
#' @param splat A `gaussian_splat` object.
#' @param grad_2d_norm Numeric vector of view-space 2D gradient norms for each Gaussian.
#'   If `NULL`, defaults to uniform values below the threshold.
#' @param grad_threshold Double. Gradient norm threshold triggering densification. Defaults to `0.0002`.
#' @param scale_threshold Double. Maximum scale distinguishing split from clone. Defaults to `0.01`.
#' @param min_opacity Double. Minimum opacity below which splats are pruned. Defaults to `0.05`.
#' @param max_scale Double. Maximum world scale above which splats are pruned. Defaults to `0.5`.
#' @param max_gaussians Integer. Maximum allowable 3D Gaussians to prevent unbounded memory growth. Defaults to `200000`.
#' @return An updated `gaussian_splat` object with densified and pruned points.
#' @export
gs_densify <- function(splat,
                       grad_2d_norm = NULL,
                       grad_threshold = 0.0002,
                       scale_threshold = 0.01,
                       min_opacity = 0.05,
                       max_scale = 0.5,
                       max_gaussians = 200000) {
  if (!inherits(splat, "gaussian_splat")) {
    stop("'splat' must be an object of class 'gaussian_splat'.")
  }

  n <- splat$num_gaussians
  if (is.null(grad_2d_norm)) {
    grad_2d_norm <- rep(0.0, n)
  }

  res <- cpp_gs_densify(
    gaussians = splat,
    grad_2d_norm = as.numeric(grad_2d_norm),
    grad_threshold = as.numeric(grad_threshold),
    scale_threshold = as.numeric(scale_threshold),
    min_opacity = as.numeric(min_opacity),
    max_scale = as.numeric(max_scale),
    max_gaussians = as.integer(max_gaussians)
  )

  splat$positions <- res$positions
  splat$scales_log <- res$scales_log
  splat$quaternions <- res$quaternions
  splat$opacities_logit <- res$opacities_logit
  splat$sh_coeffs <- res$sh_coeffs
  splat$semantic_ids <- res$semantic_ids
  splat$semantic_scores <- res$semantic_scores
  splat$num_gaussians <- res$num_gaussians

  cli::cli_alert_info("Densification/Pruning complete: {n} -> {splat$num_gaussians} Gaussians.")
  splat
}

#' @title Import Point Clouds from Agisoft Metashape, PLY, or TXT Files
#' @name gs_import_metashape
#' @description
#' Imports a sparse or dense 3D point cloud exported by **Agisoft Metashape** (or COLMAP,
#' CloudCompare, LiDAR) and initializes a `gaussian_splat` object. Automatically reads
#' 3D coordinates ($X, Y, Z$) and RGB vertex colors from binary or ASCII PLY files,
#' as well as delimited TXT/XYZ/CSV tables.
#'
#' @param file Character string. Path to the point cloud file (`.ply`, `.txt`, `.xyz`, `.csv`, `.las`).
#' @param format Character string. Format type: `"auto"` (default, inferred from extension),
#'   `"ply"`, `"txt"`, or `"las"`.
#' @param initial_scale Optional numeric scalar. If `NULL` (default), automatically estimated
#'   via spatial k-NN sampling.
#' @param initial_opacity Numeric scalar between 0 and 1. Default is `0.99` (solid points).
#' @param subsample Optional integer. If provided, randomly downsamples the point cloud
#'   to this maximum number of points for fast prototyping or rendering.
#' @return A `gaussian_splat` object ready for [gs_viewer()], [gs_render()], [gs_segment_sam()], and [gs_measure()].
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Import a point cloud exported from Agisoft Metashape
#' splat <- gs_import_metashape("plant_metashape.ply")
#'
#' # Open interactive 3D viewer directly in Chrome
#' gs_viewer(splat)
#' }
gs_import_metashape <- function(file,
                                format = c("auto", "ply", "txt", "las"),
                                initial_scale = NULL,
                                initial_opacity = 0.99,
                                subsample = NULL) {
  format <- match.arg(format)
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))
  if (format == "auto") {
    if (ext %in% c("ply")) {
      format <- "ply"
    } else if (ext %in% c("txt", "xyz", "csv")) {
      format <- "txt"
    } else if (ext %in% c("las", "laz")) {
      format <- "las"
    } else {
      format <- "ply"
    }
  }

  cli::cli_alert_info("Reading point cloud from {.path {file}} (format: {.val {format}})...")

  faces <- NULL
  has_colors <- TRUE
  if (format == "ply") {
    res <- cpp_gs_read_ply(normalizePath(file, winslash = "/"))
    pts <- res$positions
    cols <- res$colors
    if (!is.null(res$has_colors)) has_colors <- isTRUE(res$has_colors)
    if (!is.null(res$faces) && length(res$faces) > 0) {
      faces <- res$faces
    }
  } else if (format == "txt") {
    df <- if (requireNamespace("data.table", quietly = TRUE)) {
      as.data.frame(data.table::fread(file, header = FALSE))
    } else {
      utils::read.table(file, header = FALSE, stringsAsFactors = FALSE)
    }
    if (ncol(df) < 3) stop("TXT file must have at least 3 columns (X, Y, Z).")
    pts <- as.matrix(df[, 1:3])
    if (ncol(df) >= 6) {
      cols <- as.matrix(df[, 4:6])
      if (max(cols, na.rm = TRUE) > 1.0) cols <- cols / 255.0
    } else {
      cols <- NULL
      has_colors <- FALSE
    }
  } else if (format == "las") {
    if (requireNamespace("rlas", quietly = TRUE)) {
      las_data <- rlas::read.las(file)
      pts <- cbind(las_data$X, las_data$Y, las_data$Z)
      if (!is.null(las_data$R) && !is.null(las_data$G) && !is.null(las_data$B)) {
        max_c <- max(c(las_data$R, las_data$G, las_data$B), na.rm = TRUE)
        denom <- if (max_c > 255) 65535.0 else 255.0
        cols <- cbind(las_data$R / denom, las_data$G / denom, las_data$B / denom)
      } else {
        cols <- NULL
        has_colors <- FALSE
      }
    } else {
      stop("Package 'rlas' is required to read LAS/LAZ files. Please install it or export as PLY from Metashape.")
    }
  }

  n_pts <- nrow(pts)
  if (!is.null(subsample) && is.numeric(subsample) && subsample < n_pts) {
    subsample <- as.integer(subsample)
    idx <- sample.int(n_pts, subsample)
    pts <- pts[idx, , drop = FALSE]
    if (!is.null(cols)) cols <- cols[idx, , drop = FALSE]
    cli::cli_alert_info("Subsampled {n_pts} points down to {subsample} points.")
  }

  cli::cli_alert_success("Loaded {nrow(pts)} points successfully.")
  if (!is.null(faces) && length(faces) > 0) {
    cli::cli_alert_success("Detected 3D polygonal surface mesh ({format(length(faces)/3, big.mark = ',')} triangular faces).")
  }

  splat <- gs_create(
    points = pts,
    colors = cols,
    initial_scale = initial_scale,
    initial_opacity = initial_opacity
  )

  if (!is.null(faces) && length(faces) > 0) {
    attr(splat, "faces") <- faces
    splat$faces <- faces
    splat$num_faces <- as.integer(length(faces) / 3)
    splat$has_faces <- TRUE
  } else {
    splat$faces <- integer(0)
    splat$num_faces <- 0L
    splat$has_faces <- FALSE
  }
  splat$has_colors <- has_colors

  splat
}

#' @title Import Camera Calibration and Poses from Agisoft Metashape XML
#' @name gs_import_cameras_metashape
#' @description
#' Parses camera extrinsics (rotation $R$, translation $t$) and intrinsics (focal length $f_x, f_y$,
#' principal point $c_x, c_y$, dimensions) from an Agisoft Metashape exported XML file
#' (*File -> Export -> Export Cameras -> Agisoft XML*). Allows exact 2D-to-3D projection
#' and SAM segmentation of original photos.
#'
#' @param xml_file Path to the Metashape exported camera XML file.
#' @return A named list of camera specifications (each formatted as [gs_camera()]), keyed by photo label.
#' @export
gs_import_cameras_metashape <- function(xml_file) {
  if (!file.exists(xml_file)) {
    stop("XML file not found: ", xml_file)
  }

  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required to parse Metashape camera XML files.")
  }

  doc <- xml2::read_xml(xml_file)

  # Extract sensor calibration (intrinsics)
  sensors <- xml2::xml_find_all(doc, ".//sensor")
  sensor_map <- list()

  for (s in sensors) {
    s_id <- xml2::xml_attr(s, "id")
    res_node <- xml2::xml_find_first(s, ".//resolution")
    w <- if (!is.na(res_node)) as.integer(xml2::xml_attr(res_node, "width")) else 4000L
    h <- if (!is.na(res_node)) as.integer(xml2::xml_attr(res_node, "height")) else 3000L

    f_node <- xml2::xml_find_first(s, ".//calibration/f")
    f_val <- if (!is.na(f_node)) as.numeric(xml2::xml_text(f_node)) else (w * 0.8)

    cx_node <- xml2::xml_find_first(s, ".//calibration/cx")
    cy_node <- xml2::xml_find_first(s, ".//calibration/cy")
    cx_val <- (w / 2.0) + (if (!is.na(cx_node)) as.numeric(xml2::xml_text(cx_node)) else 0.0)
    cy_val <- (h / 2.0) + (if (!is.na(cy_node)) as.numeric(xml2::xml_text(cy_node)) else 0.0)

    sensor_map[[s_id]] <- list(width = w, height = h, fx = f_val, fy = f_val, cx = cx_val, cy = cy_val)
  }

  # Extract cameras (extrinsics)
  cams_nodes <- xml2::xml_find_all(doc, ".//camera")
  cameras_list <- list()

  for (cam in cams_nodes) {
    label <- xml2::xml_attr(cam, "label")
    s_id <- xml2::xml_attr(cam, "sensor_id")
    trans_node <- xml2::xml_find_first(cam, ".//transform")

    if (is.na(trans_node)) next

    trans_vals <- as.numeric(strsplit(trimws(xml2::xml_text(trans_node)), "\\s+")[[1]])
    if (length(trans_vals) != 16) next

    T_cw <- matrix(trans_vals, nrow = 4, ncol = 4, byrow = TRUE)

    R_cw <- T_cw[1:3, 1:3]
    pos_world <- T_cw[1:3, 4]

    R_wc <- t(R_cw)
    t_wc <- as.numeric(-R_wc %*% pos_world)

    sens <- if (!is.null(sensor_map[[s_id]])) sensor_map[[s_id]] else sensor_map[[1]]

    cam_spec <- list(
      R = R_wc,
      t = t_wc,
      fx = sens$fx,
      fy = sens$fy,
      cx = sens$cx,
      cy = sens$cy,
      width = sens$width,
      height = sens$height,
      position = pos_world,
      label = label
    )

    cameras_list[[label]] <- cam_spec
  }

  cli::cli_alert_success("Imported {length(cameras_list)} calibrated camera pose(s) from {.path {xml_file}}.")
  cameras_list
}

# ==============================================================================
# GPU ACCELERATION HELPERS (CUDA, PYTORCH & CUPY NVRTC)
# ==============================================================================

.gs_detect_cuda_gpu <- function() {
  py_exec <- try(.find_python_exec(), silent = TRUE)
  if (inherits(py_exec, "try-error") || !is.character(py_exec) || !nzchar(py_exec)) {
    return(list(has_gpu = FALSE, py_exec = NULL, reason = "Python executable not found in system."))
  }

  chk_py <- paste0(
    "try:\n",
    "    import torch\n",
    "    cuda_ok = torch.cuda.is_available()\n",
    "    dev_name = torch.cuda.get_device_name(0) if cuda_ok else 'None'\n",
    "    print(f'{cuda_ok}|{dev_name}')\n",
    "except Exception as e:\n",
    "    print('False|None')\n"
  )

  out <- suppressWarnings(tryCatch(
    system2(py_exec, args = c("-c", shQuote(chk_py)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  ))

  if (length(out) > 0L && grepl("\\|", out[length(out)])) {
    parts <- strsplit(trimws(out[length(out)]), "\\|")[[1]]
    if (identical(tolower(parts[1]), "true")) {
      return(list(has_gpu = TRUE, py_exec = py_exec, name = parts[2]))
    }
  }

  list(has_gpu = FALSE, py_exec = py_exec, reason = "PyTorch CUDA is not available or no NVIDIA GPU detected.")
}

.gs_write_splat_binary <- function(splat, file) {
  con <- file(file, "wb")
  on.exit(close(con))
  n <- as.integer(splat$num_gaussians)
  writeBin(n, con, size = 4)
  writeBin(as.vector(t(as.matrix(splat$positions))), con, size = 4)
  writeBin(as.vector(t(as.matrix(splat$scales_log))), con, size = 4)
  writeBin(as.vector(t(as.matrix(splat$quaternions))), con, size = 4)
  writeBin(as.numeric(splat$opacities_logit), con, size = 4)
  writeBin(as.vector(t(as.matrix(splat$sh_coeffs))), con, size = 4)

  sem_ids <- if (!is.null(splat$semantic_ids)) as.integer(splat$semantic_ids) else integer(n)
  sem_scores <- if (!is.null(splat$semantic_scores)) as.numeric(splat$semantic_scores) else numeric(n)
  writeBin(sem_ids, con, size = 4)
  writeBin(sem_scores, con, size = 4)
}

.gs_read_splat_binary <- function(file) {
  con <- file(file, "rb")
  on.exit(close(con))
  n <- readBin(con, integer(), size = 4)
  pos <- matrix(readBin(con, numeric(), n = n * 3, size = 4), nrow = n, byrow = TRUE)
  scales_log <- matrix(readBin(con, numeric(), n = n * 3, size = 4), nrow = n, byrow = TRUE)
  quaternions <- matrix(readBin(con, numeric(), n = n * 4, size = 4), nrow = n, byrow = TRUE)
  opacities_logit <- readBin(con, numeric(), n = n, size = 4)
  sh_coeffs <- matrix(readBin(con, numeric(), n = n * 3, size = 4), nrow = n, byrow = TRUE)

  fsize <- file.info(file)$size
  bytes_left <- fsize - (4 + n * (3 + 3 + 4 + 1 + 3) * 4)
  sem_ids <- if (bytes_left >= n * 4) readBin(con, integer(), n = n, size = 4) else integer(n)
  bytes_left <- bytes_left - n * 4
  sem_scores <- if (bytes_left >= n * 4) readBin(con, numeric(), n = n, size = 4) else numeric(n)

  list(
    positions = pos,
    scales_log = scales_log,
    quaternions = quaternions,
    opacities_logit = opacities_logit,
    sh_coeffs = sh_coeffs,
    semantic_ids = sem_ids,
    semantic_scores = sem_scores,
    num_gaussians = n,
    sh_degree = 0L
  )
}

.gs_fit_gpu <- function(splat,
                        dataset,
                        py_exec,
                        iterations = 500,
                        lr_pos = 0.00016,
                        lr_scale = 0.005,
                        lr_opac = 0.05,
                        lr_sh = 0.0025,
                        densify_interval = 100,
                        densify_until = round(iterations * 0.8),
                        grad_threshold = 0.0002,
                        scale_threshold = 0.01,
                        min_opacity = 0.05,
                        max_scale = 0.5,
                        max_gaussians = 200000,
                        downsample = 4,
                        verbose = TRUE) {
  # Locate gs_trainer.py
  trainer_script <- system.file("python", "gs_trainer.py", package = "pliman")
  if (!nzchar(trainer_script) || !file.exists(trainer_script)) {
    cand <- file.path(getwd(), "inst", "python", "gs_trainer.py")
    if (file.exists(cand)) trainer_script <- cand
  }
  if (!file.exists(trainer_script)) {
    cli::cli_abort("Could not find {.file gs_trainer.py}. Ensure package pliman is properly installed.")
  }

  tmp_dir <- tempfile("gs_gpu_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  splat_in_file <- file.path(tmp_dir, "splat_in.bin")
  splat_out_file <- file.path(tmp_dir, "splat_out.bin")
  cameras_json_file <- file.path(tmp_dir, "cameras.json")

  # Write input binary splat
  .gs_write_splat_binary(splat, splat_in_file)

  # Write cameras JSON
  cam_list <- lapply(dataset, function(d) {
    cam <- d$cam
    list(
      image_path = normalizePath(d$img_source, winslash = "/"),
      label = d$label,
      fx = as.numeric(cam$fx),
      fy = as.numeric(cam$fy),
      cx = as.numeric(cam$cx),
      cy = as.numeric(cam$cy),
      width = as.integer(cam$width),
      height = as.integer(cam$height),
      R = as.matrix(cam$R),
      t = as.numeric(cam$t)
    )
  })
  jsonlite::write_json(cam_list, cameras_json_file, auto_unbox = TRUE, digits = 8)

  # Build arguments for gs_trainer.py
  py_args <- c(
    shQuote(trainer_script),
    "--splat_in", shQuote(splat_in_file),
    "--cameras_json", shQuote(cameras_json_file),
    "--splat_out", shQuote(splat_out_file),
    "--iterations", sprintf("%.0f", as.numeric(iterations)),
    "--lr_pos", sprintf("%.8f", as.numeric(lr_pos)),
    "--lr_scale", sprintf("%.8f", as.numeric(lr_scale)),
    "--lr_opac", sprintf("%.8f", as.numeric(lr_opac)),
    "--lr_sh", sprintf("%.8f", as.numeric(lr_sh)),
    "--densify_interval", sprintf("%.0f", as.numeric(densify_interval)),
    "--densify_until", sprintf("%.0f", as.numeric(densify_until)),
    "--grad_threshold", sprintf("%.8f", as.numeric(grad_threshold)),
    "--scale_threshold", sprintf("%.8f", as.numeric(scale_threshold)),
    "--min_opacity", sprintf("%.8f", as.numeric(min_opacity)),
    "--max_scale", sprintf("%.8f", as.numeric(max_scale)),
    "--max_gaussians", sprintf("%.0f", as.numeric(max_gaussians)),
    "--downsample", sprintf("%.0f", as.numeric(downsample))
  )

  # Run Python trainer and stream live progress directly to console
  status <- system2(py_exec, args = py_args, stdout = if (isTRUE(verbose)) "" else FALSE, stderr = "")

  if (status != 0 || !file.exists(splat_out_file)) {
    cli::cli_abort("GPU Gaussian Splatting training failed with exit code {status}.")
  }

  # Read trained splat
  trained_splat <- .gs_read_splat_binary(splat_out_file)
  class(trained_splat) <- "gaussian_splat"
  trained_splat
}

#' @title Optimize 3D Gaussian Splats from Multi-View Photos
#' @name gs_fit
#' @description
#' Optimizes 3D Gaussian Splatting parameters (positions, scales, opacities, and colors)
#' directly from a collection of calibrated multi-view photos (e.g. captured by smartphone)
#' and camera poses (exported from Agisoft Metashape or COLMAP).
#'
#' Evaluates the differentiable tile rasterizer, computes photometric loss against ground-truth
#' photos, backpropagates analytical gradients, applies Adam optimization, and performs
#' periodic adaptive densification (split/clone) and pruning.
#'
#' Can run with high-performance GPU hardware acceleration via CUDA/PyTorch/CuPy NVRTC
#' (up to 200x faster than CPU) or fallback to native multithreaded OpenMP C++.
#'
#' @param images Character path to directory with photos, character vector of image file paths,
#'   or a list of loaded `Image` objects.
#' @param cameras Named list of camera specifications (from [gs_import_cameras_metashape()]).
#' @param splat Optional initial `gaussian_splat` object (e.g., from [gs_import_metashape()]
#'   sparse tie points). If `NULL`, initialized automatically from camera bounding volume.
#' @param iterations Integer. Total number of optimization steps. Defaults to `500`.
#' @param lr_pos Double. Learning rate for 3D Gaussian positions. Defaults to `0.00016`.
#' @param lr_scale Double. Learning rate for 3D scales. Defaults to `0.005`.
#' @param lr_opac Double. Learning rate for opacities. Defaults to `0.05`.
#' @param lr_sh Double. Learning rate for colors (Spherical Harmonics). Defaults to `0.0025`.
#' @param densify_interval Integer. Iteration frequency for adaptive densification and pruning. Defaults to `100`.
#' @param densify_until Integer. Iteration after which densification stops. Defaults to `round(iterations * 0.8)`.
#' @param grad_threshold Double. View-space 2D gradient threshold triggering split/clone. Defaults to `0.0002`.
#' @param scale_threshold Double. Maximum scale distinguishing split from clone. Defaults to `0.01`.
#' @param min_opacity Double. Minimum opacity threshold for pruning. Defaults to `0.05`.
#' @param max_scale Double. Maximum world scale for pruning. Defaults to `0.5`.
#' @param max_gaussians Integer. Maximum number of Gaussians allowed during adaptive densification. Defaults to `200000`.
#' @param downsample Integer or numeric. Downsample factor for training images (e.g. `2`, `4`, or `8`)
#'   to accelerate optimization. Defaults to `4`.
#' @param device Character. Computing device for optimization: \code{"auto"} (default: automatically uses
#'   dedicated NVIDIA GPU via CUDA if available, otherwise CPU), \code{"gpu"} (forces GPU acceleration,
#'   aborting if CUDA is unavailable), or \code{"cpu"} (native C++ OpenMP).
#' @param verbose Logical. If `TRUE`, displays progress and training metrics. Defaults to `TRUE`.
#' @return A trained `gaussian_splat` object.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # 1. Import cameras and sparse tie points exported from Metashape
#' cams <- gs_import_cameras_metashape("cameras.xml")
#' init_splat <- gs_import_metashape("sparse_points.ply")
#'
#' # 2. Train 3D Gaussian Splats on GPU (CUDA)
#' trained_splat <- gs_fit(
#'   images = "photos/",
#'   cameras = cams,
#'   splat = init_splat,
#'   iterations = 1000,
#'   device = "gpu"
#' )
#'
#' # 3. Open interactive WebGL viewer in Google Chrome
#' gs_viewer(trained_splat)
#' }
gs_fit <- function(images,
                   cameras = NULL,
                   splat = NULL,
                   iterations = 500,
                   lr_pos = 0.00016,
                   lr_scale = 0.005,
                   lr_opac = 0.05,
                   lr_sh = 0.0025,
                   densify_interval = 100,
                   densify_until = round(iterations * 0.8),
                   grad_threshold = 0.0002,
                   scale_threshold = 0.01,
                   min_opacity = 0.05,
                   max_scale = 0.5,
                   max_gaussians = 200000,
                   downsample = 4,
                   device = c("auto", "gpu", "cpu"),
                   verbose = TRUE) {
  device <- match.arg(device)
  # 0. Intelligent parameter resolution for colmap_sfm and gaussian_splat objects
  if (inherits(images, "colmap_sfm")) {
    sfm_in <- images
    col_data <- gs_import_colmap(sfm_in$recon_dir)
    if (is.null(cameras)) cameras <- col_data$cameras
    if (is.null(splat)) splat <- col_data$splat
    images <- sfm_in$images
  } else if (inherits(images, "gaussian_splat")) {
    splat_in <- images
    if (is.null(splat)) splat <- splat_in
    if (is.null(cameras)) cameras <- splat_in$cameras %||% attr(splat_in, "cameras")
    images <- splat_in$images %||% attr(splat_in, "images")
  }

  if (inherits(splat, "colmap_sfm")) {
    sfm_splat <- splat
    col_data <- gs_import_colmap(sfm_splat$recon_dir)
    if (is.null(cameras)) cameras <- col_data$cameras
    splat <- col_data$splat
    if (missing(images) || is.null(images)) images <- sfm_splat$images
  } else if (inherits(splat, "gaussian_splat")) {
    if (is.null(cameras)) cameras <- splat$cameras %||% attr(splat, "cameras")
    if (missing(images) || is.null(images)) images <- splat$images %||% attr(splat, "images")
  }

  if (inherits(cameras, "colmap_sfm")) {
    cameras <- gs_import_cameras_colmap(cameras$recon_dir)
  }

  if (is.null(cameras) || !is.list(cameras) || length(cameras) == 0) {
    stop("'cameras' must be provided as a list of camera specifications, or extracted from a 'colmap_sfm' or 'gaussian_splat' object.")
  }

  if (missing(images) || is.null(images)) {
    stop("'images' must be a directory path, vector of file paths, or embedded in the 'colmap_sfm' / 'gaussian_splat' object.")
  }

  # 1. Resolve image paths and pair with cameras
  img_paths <- character(0)
  if (is.character(images)) {
    if (length(images) == 1 && dir.exists(images)) {
      img_paths <- list.files(images, pattern = "\\.(jpg|jpeg|png|tif|tiff)$", ignore.case = TRUE, full.names = TRUE)
      if (length(img_paths) == 0) {
        stop("No image files found in directory: ", images)
      }
    } else {
      img_paths <- images[file.exists(images)]
      if (length(img_paths) == 0) {
        stop("None of the specified image files exist.")
      }
    }
  } else if (is.list(images)) {
    img_paths <- names(images)
  } else {
    stop("'images' must be a directory path, vector of file paths, or list of Image objects.")
  }

  # Build matched pairs: camera <-> image
  dataset <- list()
  cam_names <- names(cameras)
  if (is.null(cam_names)) cam_names <- paste0("cam_", seq_along(cameras))

  for (i in seq_along(cameras)) {
    c_spec <- cameras[[i]]
    c_label <- if (!is.null(c_spec$label)) c_spec$label else cam_names[i]

    # Find matching image file
    matched_img <- NULL
    if (is.list(images) && !is.character(images)) {
      if (!is.null(images[[c_label]])) {
        matched_img <- images[[c_label]]
      } else if (i <= length(images)) {
        matched_img <- images[[i]]
      }
    } else {
      c_base <- basename(c_label)
      c_sans <- tools::file_path_sans_ext(c_base)

      for (ip in img_paths) {
        f_base <- basename(ip)
        f_sans <- tools::file_path_sans_ext(f_base)
        if (tolower(f_base) == tolower(c_base) || tolower(f_sans) == tolower(c_sans)) {
          matched_img <- ip
          break
        }
      }
      # Fallback to index if counts match
      if (is.null(matched_img) && length(img_paths) == length(cameras)) {
        matched_img <- img_paths[i]
      }
    }

    if (!is.null(matched_img)) {
      dataset[[length(dataset) + 1]] <- list(cam = c_spec, img_source = matched_img, label = c_label)
    }
  }

  if (length(dataset) == 0) {
    stop("Could not match any camera specifications with the provided images. Check file names / labels.")
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Matched {length(dataset)} multi-view camera pose(s) with ground-truth photos.")
  }

  # 2. Initialize Gaussian Splats if not provided
  if (is.null(splat)) {
    if (isTRUE(verbose)) cli::cli_alert_info("Initializing 3D Gaussians from camera trajectory...")
    cam_pos <- do.call(rbind, lapply(dataset, function(d) d$cam$position))
    center <- colMeans(cam_pos)
    dists <- apply(cam_pos, 1, function(p) sqrt(sum((p - center)^2)))
    mean_dist <- mean(dists)
    n_init <- 2500
    pts <- matrix(stats::rnorm(n_init * 3, mean = 0, sd = max(0.1, mean_dist * 0.35)), ncol = 3)
    pts[, 1] <- pts[, 1] + center[1]
    pts[, 2] <- pts[, 2] + center[2]
    pts[, 3] <- pts[, 3] + center[3]
    splat <- gs_create(pts, initial_scale = max(0.01, mean_dist * 0.03), initial_opacity = 0.25)
  } else {
    if (!inherits(splat, "gaussian_splat")) {
      stop("'splat' must be an object of class 'gaussian_splat'.")
    }
  }

  # 2.5 Resolve computing device and branch to GPU backend if applicable
  gpu_info <- .gs_detect_cuda_gpu()
  use_gpu <- FALSE
  if (device == "auto") {
    use_gpu <- isTRUE(gpu_info$has_gpu)
  } else if (device == "gpu") {
    if (!isTRUE(gpu_info$has_gpu)) {
      cli::cli_abort(paste0(
        "GPU execution requested (device = 'gpu'), but no CUDA-enabled GPU was detected. ",
        gpu_info$reason %||% ""
      ))
    }
    use_gpu <- TRUE
  }

  if (use_gpu) {
    if (isTRUE(verbose)) {
      cli::cli_rule("Optimizing 3D Gaussian Splats (GPU CUDA Acceleration)")
      cli::cli_alert_success("NVIDIA GPU detected: {.strong {gpu_info$name}} (CUDA backend) \U0001f4a5")
      cli::cli_alert_info("Initial Gaussians: {format(splat$num_gaussians, big.mark = ',')} | Total Steps: {iterations} | Views: {length(dataset)}")
    }

    splat_trained <- .gs_fit_gpu(
      splat = splat,
      dataset = dataset,
      py_exec = gpu_info$py_exec,
      iterations = iterations,
      lr_pos = lr_pos,
      lr_scale = lr_scale,
      lr_opac = lr_opac,
      lr_sh = lr_sh,
      densify_interval = densify_interval,
      densify_until = densify_until,
      grad_threshold = grad_threshold,
      scale_threshold = scale_threshold,
      min_opacity = min_opacity,
      max_scale = max_scale,
      max_gaussians = max_gaussians,
      downsample = downsample,
      verbose = verbose
    )

    splat_trained$cameras <- cameras
    splat_trained$images <- if (is.character(images)) images else NULL
    attr(splat_trained, "cameras") <- cameras
    attr(splat_trained, "images") <- splat_trained$images
    class(splat_trained) <- "gaussian_splat"
    return(splat_trained)
  }

  # 3. Preload and downsample training views
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Preparing training views (downsample factor: {downsample}x)...")
  }

  training_views <- list()
  for (idx in seq_along(dataset)) {
    d <- dataset[[idx]]
    cam <- d$cam

    img_obj <- if (inherits(d$img_source, c("Image", "image"))) {
      d$img_source
    } else {
      image_import(d$img_source)
    }

    if (downsample > 1) {
      target_rel <- max(5, round(100 / downsample))
      img_obj <- image_resize(img_obj, rel_size = target_rel)
    }

    new_w <- as.integer(dim(img_obj)[1])
    new_h <- as.integer(dim(img_obj)[2])

    cam_scaled <- cam
    cam_scaled$width <- new_w
    cam_scaled$height <- new_h
    scale_w <- new_w / cam$width
    scale_h <- new_h / cam$height
    cam_scaled$fx <- cam$fx * scale_w
    cam_scaled$fy <- cam$fy * scale_h
    cam_scaled$cx <- cam$cx * scale_w
    cam_scaled$cy <- cam$cy * scale_h

    raw_img <- image_data(img_obj)
    img_num <- as.numeric(raw_img)
    if (is.raw(raw_img) || max(img_num, na.rm = TRUE) > 1.0) {
      img_num <- img_num / 255.0
    }
    img_tensor <- img_num

    training_views[[idx]] <- list(
      cam = cam_scaled,
      target = img_tensor,
      label = d$label
    )
  }

  num_views <- length(training_views)

  # 4. Initialize Adam optimizer moments
  init_adam <- function(n_gauss, n_sh) {
    list(
      m = list(
        pos = matrix(0.0, nrow = n_gauss, ncol = 3),
        scale = matrix(0.0, nrow = n_gauss, ncol = 3),
        opac = numeric(n_gauss),
        sh = matrix(0.0, nrow = n_gauss, ncol = n_sh)
      ),
      v = list(
        pos = matrix(0.0, nrow = n_gauss, ncol = 3),
        scale = matrix(0.0, nrow = n_gauss, ncol = 3),
        opac = numeric(n_gauss),
        sh = matrix(0.0, nrow = n_gauss, ncol = n_sh)
      )
    )
  }

  num_sh <- ncol(splat$sh_coeffs)
  moments <- init_adam(splat$num_gaussians, num_sh)
  adam_m <- moments$m
  adam_v <- moments$v

  accum_grad_2d <- numeric(splat$num_gaussians)
  accum_count <- 0
  recent_losses <- numeric(0)

  if (isTRUE(verbose)) {
    cli::cli_rule("Optimizing 3D Gaussian Splats")
    cli::cli_alert_info("Initial Gaussians: {format(splat$num_gaussians, big.mark = ',')} | Total Steps: {iterations} | Views: {num_views}")
  }

  # 5. Optimization Loop
  start_time <- Sys.time()

  for (step in seq_len(iterations)) {
    v_idx <- ((step - 1) %% num_views) + 1
    view <- training_views[[v_idx]]

    step_res <- cpp_gs_train_step(
      camera = view$cam,
      gaussians = splat,
      target_img = view$target,
      adam_m = adam_m,
      adam_v = adam_v,
      iter = as.integer(step),
      lr_pos = as.numeric(lr_pos),
      lr_scale = as.numeric(lr_scale),
      lr_opac = as.numeric(lr_opac),
      lr_sh = as.numeric(lr_sh)
    )

    splat <- step_res$gaussians
    class(splat) <- "gaussian_splat"
    adam_m <- step_res$adam_m
    adam_v <- step_res$adam_v

    g_norm <- step_res$grad_2d_norm
    if (length(g_norm) == length(accum_grad_2d)) {
      accum_grad_2d <- accum_grad_2d + g_norm
      accum_count <- accum_count + 1
    }

    recent_losses <- c(recent_losses, step_res$loss)
    if (length(recent_losses) > 20) recent_losses <- utils::tail(recent_losses, 20)

    # Periodic Adaptive Densification and Pruning
    if (step %% densify_interval == 0 && step <= densify_until) {
      mean_grads <- accum_grad_2d / max(1, accum_count)
      splat <- gs_densify(
        splat = splat,
        grad_2d_norm = mean_grads,
        grad_threshold = grad_threshold,
        scale_threshold = scale_threshold,
        min_opacity = min_opacity,
        max_scale = max_scale,
        max_gaussians = max_gaussians
      )

      # Reset Adam momentum for updated Gaussian count
      moments <- init_adam(splat$num_gaussians, num_sh)
      adam_m <- moments$m
      adam_v <- moments$v
      accum_grad_2d <- numeric(splat$num_gaussians)
      accum_count <- 0

      if (isTRUE(verbose)) {
        cli::cli_alert_info("Step {step}/{iterations} | Densified to {format(splat$num_gaussians, big.mark = ',')} Gaussians (Loss: {round(mean(recent_losses), 4)})")
      }
    }

    # Dynamic ETA and live progress feedback
    log_interval <- max(1, min(25, round(iterations / 20)))
    if (isTRUE(verbose) && (step %% log_interval == 0 || step == iterations)) {
      curr_time <- Sys.time()
      elapsed_s <- as.numeric(difftime(curr_time, start_time, units = "secs"))
      it_per_sec <- step / max(0.001, elapsed_s)
      remaining_s <- (iterations - step) / max(0.001, it_per_sec)
      eta_str <- if (remaining_s >= 60) paste0(round(remaining_s / 60, 1), "m") else paste0(round(remaining_s), "s")
      pct <- round((step / iterations) * 100)
      cli::cli_alert_info("Step {step}/{iterations} [{pct}%] | {round(it_per_sec, 1)} it/s | ETA: {eta_str} | Loss: {round(mean(recent_losses), 4)} | Gaussians: {format(splat$num_gaussians, big.mark = ',')}")
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  if (isTRUE(verbose)) {
    cli::cli_alert_success("Training finished in {elapsed}s! Final Gaussians: {format(splat$num_gaussians, big.mark = ',')}")
    cli::cli_rule()
  }

  splat$cameras <- cameras
  splat$images <- if (is.character(images)) images else NULL
  attr(splat, "cameras") <- cameras
  attr(splat, "images") <- splat$images
  class(splat) <- "gaussian_splat"
  splat
}

#' @title Statistical and Radius Outlier Removal (SOR / ROR) for 3D Point Clouds
#' @name gs_filter_outliers
#' @description
#' Filters out noisy, floating tie-points and airborne outliers from point clouds
#' or Gaussian Splats based on local neighborhood density using a high-performance native
#' C++ 3D Kd-Tree (OpenMP accelerated, zero external package dependencies).
#'
#' @param splat A `gaussian_splat`, `colmap_sfm`, or numeric matrix of 3D coordinates.
#' @param k Integer. Number of nearest neighbors to evaluate (e.g. 16 to 30). Defaults to `20`.
#' @param std_ratio Numeric. Standard deviation multiplier threshold. Smaller values (e.g. 0.8 to 1.2)
#'   filter more aggressively, removing scattered floating points around the subject. Defaults to `1.0`.
#' @param radius Numeric. Optional sphere radius for density filtering. If > 0, points with fewer than
#'   `min_neighbors` within `radius` are additionally filtered out. Defaults to `0.0` (disabled).
#' @param min_neighbors Integer. Minimum neighbors required within `radius`. Defaults to `0` (disabled).
#' @param verbose Logical. Show progress and statistics. Defaults to `TRUE`.
#' @return An updated and cleaned `gaussian_splat`, `colmap_sfm`, or point matrix.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Clean floating noise from a COLMAP reconstruction
#' sfm_clean <- gs_filter_outliers(sfm_completo, k = 20, std_ratio = 1.0)
#' plot(sfm_clean)
#' }
gs_filter_outliers <- function(splat, k = 20, std_ratio = 1.0, radius = 0.0, min_neighbors = 0, verbose = TRUE) {
  is_colmap <- inherits(splat, "colmap_sfm")
  cameras_saved <- NULL
  images_saved <- NULL

  if (is_colmap) {
    images_saved <- splat$images
    col_data <- gs_import_colmap(splat$recon_dir)
    cameras_saved <- col_data$cameras
    splat <- col_data$splat
  } else if (inherits(splat, "gaussian_splat")) {
    cameras_saved <- splat$cameras %||% attr(splat, "cameras")
    images_saved <- splat$images %||% attr(splat, "images")
  }

  is_splat <- inherits(splat, "gaussian_splat")
  pts <- if (is_splat) splat$positions else as.matrix(splat)

  n <- nrow(pts)
  if (n <= k) return(splat)

  filter_res <- cpp_gs_filter_outliers(
    points = pts,
    k = as.integer(k),
    std_ratio = as.numeric(std_ratio),
    radius = as.numeric(radius),
    min_neighbors = as.integer(min_neighbors)
  )

  keep <- filter_res$keep
  n_kept <- filter_res$num_kept
  n_removed <- filter_res$num_removed

  if (isTRUE(verbose)) {
    pct <- round((n_kept / n) * 100, 1)
    cli::cli_alert_success("Statistical Outlier Filter (SOR): kept {format(n_kept, big.mark = ',')} points ({pct}%), removed {format(n_removed, big.mark = ',')} floating noise points.")
  }

  if (is_splat) {
    splat$positions <- splat$positions[keep, , drop = FALSE]
    splat$scales_log <- splat$scales_log[keep, , drop = FALSE]
    splat$quaternions <- splat$quaternions[keep, , drop = FALSE]
    splat$opacities_logit <- splat$opacities_logit[keep]
    splat$sh_coeffs <- splat$sh_coeffs[keep, , drop = FALSE]
    splat$semantic_ids <- splat$semantic_ids[keep]
    splat$semantic_scores <- splat$semantic_scores[keep]
    splat$num_gaussians <- n_kept
    attr(splat, "density_dists") <- filter_res$mean_dists[keep]

    if (!is.null(cameras_saved)) {
      attr(splat, "cameras") <- cameras_saved
      splat$cameras <- cameras_saved
    }
    if (!is.null(images_saved)) {
      attr(splat, "images") <- images_saved
      splat$images <- images_saved
    }
    return(splat)
  } else {
    res_mat <- pts[keep, , drop = FALSE]
    attr(res_mat, "density_dists") <- filter_res$mean_dists[keep]
    return(res_mat)
  }
}

#' @title Crop Point Clouds and 3D Gaussian Splats by 3D Bounding Box
#' @name gs_crop
#' @description
#' Crops a `gaussian_splat`, `colmap_sfm`, or 3D point matrix to retain only points
#' located inside a specified 3D Oriented Bounding Box (OBB) or Axis-Aligned Bounding Box (AABB).
#'
#' @param splat An object of class `gaussian_splat`, `colmap_sfm`, or a numeric matrix of 3D coordinates.
#' @param xlim Optional numeric vector of length 2 `c(min, max)` for X axis bounds.
#' @param ylim Optional numeric vector of length 2 `c(min, max)` for Y axis bounds.
#' @param zlim Optional numeric vector of length 2 `c(min, max)` for Z axis bounds.
#' @param center Optional numeric vector of length 3 `c(x, y, z)` center of the bounding box.
#' @param size Optional numeric vector of length 3 `c(dx, dy, dz)` dimensions of the bounding box.
#' @param rotation Optional numeric vector of length 3 `c(roll, pitch, yaw)` in degrees for an oriented bounding box.
#' @param verbose Logical. Print summary statistics. Defaults to `TRUE`.
#' @return An updated, cropped object of the same class as `splat`.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Crop around plant subject
#' splat_cropped <- gs_crop(sfm_completo, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), zlim = c(-0.2, 0.8))
#' plot(splat_cropped)
#' }
gs_crop <- function(splat,
                    xlim = NULL,
                    ylim = NULL,
                    zlim = NULL,
                    center = NULL,
                    size = NULL,
                    rotation = c(0, 0, 0),
                    verbose = TRUE) {
  is_colmap <- inherits(splat, "colmap_sfm")
  cameras_saved <- NULL
  images_saved <- NULL

  if (is_colmap) {
    images_saved <- splat$images
    col_data <- gs_import_colmap(splat$recon_dir)
    cameras_saved <- col_data$cameras
    splat <- col_data$splat
  } else if (inherits(splat, "gaussian_splat")) {
    cameras_saved <- splat$cameras %||% attr(splat, "cameras")
    images_saved <- splat$images %||% attr(splat, "images")
  }

  is_splat <- inherits(splat, "gaussian_splat")
  pts <- if (is_splat) splat$positions else as.matrix(splat)
  n <- nrow(pts)

  # Determine inside/outside mask
  keep <- rep(TRUE, n)

  if (!is.null(center) && !is.null(size)) {
    # Oriented Bounding Box test
    c_x <- center[1]; c_y <- center[2]; c_z <- center[3]
    s_x <- size[1]; s_y <- size[2]; s_z <- size[3]

    # Convert Euler angles (deg) to radians
    rx <- rotation[1] * pi / 180
    ry <- rotation[2] * pi / 180
    rz <- rotation[3] * pi / 180

    # Rotation matrices
    Rx <- matrix(c(1, 0, 0, 0, cos(rx), -sin(rx), 0, sin(rx), cos(rx)), 3, 3, byrow = TRUE)
    Ry <- matrix(c(cos(ry), 0, sin(ry), 0, 1, 0, -sin(ry), 0, cos(ry)), 3, 3, byrow = TRUE)
    Rz <- matrix(c(cos(rz), -sin(rz), 0, sin(rz), cos(rz), 0, 0, 0, 1), 3, 3, byrow = TRUE)
    R <- Rz %*% Ry %*% Rx
    invR <- t(R)

    # Shift to center and rotate to local box space
    diff_pts <- t(t(pts) - c(c_x, c_y, c_z))
    local_pts <- diff_pts %*% t(invR)

    keep <- (abs(local_pts[, 1]) <= s_x / 2) &
            (abs(local_pts[, 2]) <= s_y / 2) &
            (abs(local_pts[, 3]) <= s_z / 2)
  } else {
    # Axis-aligned box bounds
    if (!is.null(xlim) && length(xlim) == 2) {
      keep <- keep & (pts[, 1] >= min(xlim) & pts[, 1] <= max(xlim))
    }
    if (!is.null(ylim) && length(ylim) == 2) {
      keep <- keep & (pts[, 2] >= min(ylim) & pts[, 2] <= max(ylim))
    }
    if (!is.null(zlim) && length(zlim) == 2) {
      keep <- keep & (pts[, 3] >= min(zlim) & pts[, 3] <= max(zlim))
    }
  }

  n_kept <- sum(keep)
  n_removed <- n - n_kept

  if (isTRUE(verbose)) {
    pct <- round((n_kept / max(1, n)) * 100, 1)
    cli::cli_alert_success("3D Box Crop: kept {format(n_kept, big.mark = ',')} points ({pct}%), trimmed {format(n_removed, big.mark = ',')} outside points.")
  }

  if (is_splat) {
    splat$positions <- splat$positions[keep, , drop = FALSE]
    splat$scales_log <- splat$scales_log[keep, , drop = FALSE]
    splat$quaternions <- splat$quaternions[keep, , drop = FALSE]
    splat$opacities_logit <- splat$opacities_logit[keep]
    splat$sh_coeffs <- splat$sh_coeffs[keep, , drop = FALSE]
    splat$semantic_ids <- splat$semantic_ids[keep]
    splat$semantic_scores <- splat$semantic_scores[keep]
    splat$num_gaussians <- n_kept

    if (!is.null(cameras_saved)) {
      attr(splat, "cameras") <- cameras_saved
      splat$cameras <- cameras_saved
    }
    if (!is.null(images_saved)) {
      attr(splat, "images") <- images_saved
      splat$images <- images_saved
    }
    return(splat)
  } else {
    return(pts[keep, , drop = FALSE])
  }
}


