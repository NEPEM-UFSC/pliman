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
#'   computed per point as the root-mean-square distance to its 3 nearest neighbors
#'   (the initialization of Kerbl et al., 2023).
#' @param initial_opacity Numeric scalar between 0 and 1. Initial Gaussian opacity. Defaults to `0.15`.
#' @param sh_degree Integer. Maximum Spherical Harmonics degree (0, 1, 2, or 3). Defaults to `0`.
#' @return An S3 object of class `gaussian_splat`.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Create 500 random points representing a small plant canopy
#' pts <- matrix(rnorm(1500, mean = 0, sd = 0.5), ncol = 3)
#' splat <- gs_create(pts)
#' print(splat)
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

  # Initial log-scales. As in Kerbl et al. (2023), each Gaussian starts as an
  # isotropic blob whose radius is the root-mean-square distance to its three
  # nearest neighbours, so the initial cloud has no holes regardless of the
  # (arbitrary) metric scale of the reconstruction.
  if (is.null(initial_scale)) {
    init_s <- .gs_knn_scales(points)
  } else {
    init_s <- rep(max(1e-7, as.numeric(initial_scale)), n)
  }

  scales_log <- matrix(log(init_s), nrow = n, ncol = 3)

  # Identity quaternions: (qw = 1, qx = 0, qy = 0, qz = 0)
  quaternions <- matrix(0.0, nrow = n, ncol = 4)
  quaternions[, 1] <- 1.0

  # Opacity logit
  init_opac <- max(0.01, min(0.99, initial_opacity))
  logit_opac <- log(init_opac / (1.0 - init_opac))
  opacities_logit <- rep(logit_opac, n)

  # Spherical Harmonics degree 0: C = 0.5 + 0.28209479 * SH0 -> SH0 = (C - 0.5) / 0.28209479
  sh_c0 <- 0.28209479177387814
  sh_degree <- max(0L, min(3L, as.integer(sh_degree)))
  num_sh <- 3L * (sh_degree + 1L)^2
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

# Per-point initial scale: RMS distance to the 3 nearest neighbours (3DGS init)
.gs_knn_scales <- function(points) {
  points <- as.matrix(points[, 1:3, drop = FALSE])
  storage.mode(points) <- "double"
  n <- nrow(points)
  if (n < 4) return(rep(0.02, n))
  knn <- cpp_gs_knn(points, queries = NULL, k = 3L)
  d2 <- rowMeans(knn$dists^2)
  s <- sqrt(d2)
  pos_s <- s[is.finite(s) & s > 0]
  floor_s <- if (length(pos_s) > 0) max(1e-7, stats::quantile(pos_s, 0.01, names = FALSE)) else 1e-3
  s[!is.finite(s) | s < floor_s] <- floor_s
  s
}

# Subset every per-Gaussian field of a gaussian_splat, keeping all other metadata
.gs_subset <- function(splat, keep) {
  splat$positions <- splat$positions[keep, , drop = FALSE]
  splat$scales_log <- splat$scales_log[keep, , drop = FALSE]
  splat$quaternions <- splat$quaternions[keep, , drop = FALSE]
  splat$opacities_logit <- splat$opacities_logit[keep]
  splat$sh_coeffs <- splat$sh_coeffs[keep, , drop = FALSE]
  if (!is.null(splat$semantic_ids)) splat$semantic_ids <- splat$semantic_ids[keep]
  if (!is.null(splat$semantic_scores)) splat$semantic_scores <- splat$semantic_scores[keep]
  splat$num_gaussians <- nrow(splat$positions)
  splat
}

# Number of SH bands encoded in an SH coefficient matrix with 3 * (deg + 1)^2 columns
.gs_sh_degree_from_cols <- function(ncols) {
  k <- ncols / 3
  deg <- round(sqrt(k)) - 1
  as.integer(max(0, min(3, deg)))
}

# Scene extent following the NeRF++ normalisation used by 3DGS:
# 1.1 x the largest distance from a camera centre to the mean camera centre.
.gs_scene_extent <- function(cameras) {
  pos <- do.call(rbind, lapply(cameras, function(cm) {
    if (!is.null(cm$position)) as.numeric(cm$position) else as.numeric(-t(cm$R) %*% cm$t)
  }))
  if (is.null(pos) || nrow(pos) < 2) return(1.0)
  ctr <- colMeans(pos)
  ext <- max(sqrt(rowSums(sweep(pos, 2, ctr)^2))) * 1.1
  if (!is.finite(ext) || ext <= 0) 1.0 else ext
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
  if (isTRUE(x$trained)) {
    cat(sprintf("Optimized:            Yes (%s views%s)\n",
                if (!is.null(x$cameras)) length(x$cameras) else "?",
                if (!is.null(x$train_psnr) && is.finite(x$train_psnr)) sprintf(", train PSNR %.2f dB", x$train_psnr) else ""))
  }

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

  n <- as.integer(splat$num_gaussians)
  sh <- as.matrix(splat$sh_coeffs)
  n_bands <- ncol(sh) %/% 3                  # (deg + 1)^2 coefficients per channel
  n_rest <- n_bands - 1L

  # pliman stores SH interleaved by basis: [dc_r, dc_g, dc_b, b1_r, b1_g, b1_b, ...].
  # The 3DGS PLY convention stores f_rest channel-major: all R bands, then G, then B.
  f_rest <- NULL
  if (n_rest > 0) {
    rest_idx <- unlist(lapply(0:2, function(ch) 3L + 3L * (seq_len(n_rest) - 1L) + ch + 1L))
    f_rest <- sh[, rest_idx, drop = FALSE]
  }

  float_mat <- cbind(
    splat$positions[, 1:3, drop = FALSE],
    matrix(0, n, 3),                         # normals (unused, kept for compatibility)
    sh[, 1:3, drop = FALSE],
    f_rest,
    splat$opacities_logit,
    splat$scales_log[, 1:3, drop = FALSE],
    splat$quaternions[, 1:4, drop = FALSE]
  )
  storage.mode(float_mat) <- "double"
  float_mat[!is.finite(float_mat)] <- 0

  props <- c("x", "y", "z", "nx", "ny", "nz", "f_dc_0", "f_dc_1", "f_dc_2",
             if (n_rest > 0) paste0("f_rest_", seq_len(3 * n_rest) - 1L),
             "opacity", "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3")

  header <- paste0(
    "ply\n",
    "format binary_little_endian 1.0\n",
    "comment Exported by pliman (3D Gaussian Splatting)\n",
    "element vertex ", n, "\n",
    paste0("property float ", props, "\n", collapse = ""),
    "property int semantic_id\n",
    "end_header\n"
  )

  # Vectorised binary packing: float block + int32 semantic id per vertex
  float_raw <- writeBin(as.vector(t(float_mat)), raw(), size = 4, endian = "little")
  sem <- as.integer(splat$semantic_ids %||% integer(n))
  if (length(sem) != n) sem <- integer(n)
  sem_raw <- writeBin(sem, raw(), size = 4, endian = "little")
  body <- rbind(matrix(float_raw, nrow = 4L * ncol(float_mat)), matrix(sem_raw, nrow = 4L))

  con <- file(file, "wb")
  on.exit(close(con))
  writeChar(header, con, eos = NULL)
  writeBin(as.vector(body), con)

  cli::cli_alert_success("Exported {format(n, big.mark = ',')} Gaussian splat(s) to {.path {file}}.")
  invisible(file)
}

# Read a 3D Gaussian Splatting PLY (Kerbl et al. format, as written by
# gs_export_ply(), the original 3DGS code, Nerfstudio, gsplat, Postshot ...)
.gs_read_3dgs_ply <- function(file) {
  con <- file(file, "rb")
  on.exit(close(con))
  header <- character(0)
  repeat {
    ln <- readLines(con, n = 1, warn = FALSE)
    if (length(ln) == 0) stop("Unexpected end of file while reading PLY header.")
    ln <- sub("\r$", "", ln)
    header <- c(header, ln)
    if (ln == "end_header") break
  }
  if (!any(grepl("binary_little_endian", header))) {
    stop("Only binary little-endian 3D Gaussian Splatting PLY files are supported.")
  }
  n <- as.integer(sub("element vertex ", "", grep("^element vertex", header, value = TRUE)[1]))
  v_start <- grep("^element vertex", header)[1]
  nxt <- grep("^element ", header)
  v_end <- if (any(nxt > v_start)) min(nxt[nxt > v_start]) - 1 else length(header)
  prop_lines <- grep("^property ", header[v_start:v_end], value = TRUE)
  ptype <- vapply(strsplit(prop_lines, "\\s+"), `[`, "", 2)
  pname <- vapply(strsplit(prop_lines, "\\s+"), `[`, "", 3)
  psize <- ifelse(ptype %in% c("double", "float64"), 8L,
                  ifelse(ptype %in% c("uchar", "uint8", "char", "int8"), 1L,
                         ifelse(ptype %in% c("short", "ushort", "int16", "uint16"), 2L, 4L)))
  stride <- sum(psize)
  raw_all <- readBin(con, "raw", n = as.double(stride) * n)
  if (length(raw_all) < stride * n) stop("PLY body is truncated.")
  rm <- matrix(raw_all, nrow = stride)
  offs <- cumsum(c(0L, psize))[seq_along(psize)]
  col <- function(nm) {
    j <- match(nm, pname)
    if (is.na(j)) return(NULL)
    bytes <- as.vector(rm[(offs[j] + 1):(offs[j] + psize[j]), , drop = FALSE])
    tp <- ptype[j]
    if (tp %in% c("float", "float32")) readBin(bytes, "double", n = n, size = 4, endian = "little")
    else if (tp %in% c("double", "float64")) readBin(bytes, "double", n = n, size = 8, endian = "little")
    else if (tp %in% c("int", "int32", "uint", "uint32")) readBin(bytes, "integer", n = n, size = 4, endian = "little")
    else if (tp %in% c("uchar", "uint8")) as.integer(bytes)
    else readBin(bytes, "integer", n = n, size = psize[j], endian = "little")
  }
  pos <- cbind(col("x"), col("y"), col("z"))
  dc <- cbind(col("f_dc_0"), col("f_dc_1"), col("f_dc_2"))
  rest_names <- grep("^f_rest_", pname, value = TRUE)
  rest_names <- rest_names[order(as.integer(sub("f_rest_", "", rest_names)))]
  n_rest <- length(rest_names) %/% 3
  sh <- dc
  if (n_rest > 0) {
    rest <- do.call(cbind, lapply(rest_names, col))
    inter <- matrix(0, n, 3 * n_rest)
    for (k in seq_len(n_rest)) {
      for (ch in 0:2) inter[, 3 * (k - 1) + ch + 1] <- rest[, ch * n_rest + k]
    }
    sh <- cbind(dc, inter)
  }
  sem <- col("semantic_id")
  obj <- list(
    positions = pos,
    scales_log = cbind(col("scale_0"), col("scale_1"), col("scale_2")),
    quaternions = cbind(col("rot_0"), col("rot_1"), col("rot_2"), col("rot_3")),
    opacities_logit = col("opacity"),
    sh_coeffs = sh,
    semantic_ids = if (is.null(sem)) integer(n) else as.integer(sem),
    semantic_scores = numeric(n),
    sh_degree = .gs_sh_degree_from_cols(ncol(sh)),
    num_gaussians = n,
    faces = integer(0),
    num_faces = 0L,
    has_faces = FALSE,
    has_colors = TRUE,
    trained = TRUE
  )
  class(obj) <- "gaussian_splat"
  obj
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
#' @param verbose Logical. If `TRUE` (default), reports the change in the number of Gaussians.
#' @return An updated `gaussian_splat` object with densified and pruned points.
#' @export
gs_densify <- function(splat,
                       grad_2d_norm = NULL,
                       grad_threshold = 0.0002,
                       scale_threshold = 0.01,
                       min_opacity = 0.05,
                       max_scale = 0.5,
                       max_gaussians = 200000,
                       verbose = TRUE) {
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

  if (isTRUE(verbose)) cli::cli_alert_info("Densification/Pruning complete: {n} -> {splat$num_gaussians} Gaussians.")
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
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Import a point cloud exported from Agisoft Metashape
#' splat <- gs_import_metashape("plant_metashape.ply")
#'
#' # Open interactive 3D viewer directly in Chrome
#' gs_viewer(splat)
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

  # Trained 3D Gaussian Splatting PLY (f_dc_*, scale_*, rot_*, opacity): keep every attribute
  if (format == "ply") {
    hdr <- tryCatch({
      con_h <- file(file, "rb")
      on.exit(close(con_h), add = TRUE)
      readBin(con_h, "raw", n = 8192)
    }, error = function(e) raw(0))
    hdr_txt <- rawToChar(hdr[hdr != as.raw(0)])
    hdr_end <- regexpr("end_header", hdr_txt, fixed = TRUE, useBytes = TRUE)
    if (hdr_end > 0) hdr_txt <- substr(hdr_txt, 1, hdr_end)
    if (grepl("property float f_dc_0", hdr_txt, fixed = TRUE, useBytes = TRUE) &&
        grepl("property float rot_0", hdr_txt, fixed = TRUE, useBytes = TRUE)) {
      cli::cli_alert_info("Reading trained 3D Gaussian Splatting PLY from {.path {file}}...")
      splat <- .gs_read_3dgs_ply(file)
      if (!is.null(subsample) && is.numeric(subsample) && subsample < splat$num_gaussians) {
        idx <- sort(sample.int(splat$num_gaussians, as.integer(subsample)))
        splat <- .gs_subset(splat, idx)
      }
      cli::cli_alert_success("Loaded {format(splat$num_gaussians, big.mark = ',')} Gaussians (SH degree {splat$sh_degree}).")
      return(splat)
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
    df <- utils::read.table(file, header = FALSE, stringsAsFactors = FALSE)
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
  # v2 header: magic "GSV2", N, number of SH columns
  writeBin(c(0x32565347L, n, as.integer(ncol(as.matrix(splat$sh_coeffs)))), con, size = 4)
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
  head <- readBin(con, integer(), size = 4)
  hdr_bytes <- 4
  if (identical(head, 0x32565347L)) {
    n <- readBin(con, integer(), size = 4)
    k <- readBin(con, integer(), size = 4)
    hdr_bytes <- 12
  } else {
    n <- head
    k <- 3L
  }
  pos <- matrix(readBin(con, numeric(), n = n * 3, size = 4), nrow = n, byrow = TRUE)
  scales_log <- matrix(readBin(con, numeric(), n = n * 3, size = 4), nrow = n, byrow = TRUE)
  quaternions <- matrix(readBin(con, numeric(), n = n * 4, size = 4), nrow = n, byrow = TRUE)
  opacities_logit <- readBin(con, numeric(), n = n, size = 4)
  sh_coeffs <- matrix(readBin(con, numeric(), n = n * k, size = 4), nrow = n, byrow = TRUE)

  fsize <- file.info(file)$size
  bytes_left <- fsize - (hdr_bytes + n * (3 + 3 + 4 + 1 + k) * 4)
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
    sh_degree = .gs_sh_degree_from_cols(k)
  )
}

.gs_fit_gpu <- function(splat,
                        dataset,
                        py_exec,
                        iterations = 7000,
                        lr_pos = 0.00016,
                        lr_scale = 0.005,
                        lr_opac = 0.05,
                        lr_sh = 0.0025,
                        lr_rot = 0.001,
                        densify_interval = 100,
                        densify_from = NULL,
                        densify_until = round(iterations * 0.5),
                        opacity_reset_interval = 3000,
                        grad_threshold = 0.0002,
                        scale_threshold = 0.01,
                        min_opacity = 0.005,
                        max_scale = 0.1,
                        max_gaussians = 500000,
                        downsample = 4,
                        sh_degree = 3,
                        lambda_dssim = 0.2,
                        white_background = FALSE,
                        seed = 42,
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
  stats_json_file <- file.path(tmp_dir, "stats.json")

  # Write input binary splat
  .gs_write_splat_binary(splat, splat_in_file)

  # Write cameras JSON (R is serialised row-major by jsonlite)
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
  jsonlite::write_json(cam_list, cameras_json_file, auto_unbox = TRUE, digits = NA, matrix = "rowmajor")

  fmt_num <- function(x) format(as.numeric(x), scientific = FALSE, digits = 15, trim = TRUE)
  py_args <- c(
    shQuote(trainer_script),
    "--splat_in", shQuote(splat_in_file),
    "--cameras_json", shQuote(cameras_json_file),
    "--splat_out", shQuote(splat_out_file),
    "--stats_json", shQuote(stats_json_file),
    "--iterations", sprintf("%.0f", as.numeric(iterations)),
    "--lr_pos", fmt_num(lr_pos),
    "--lr_pos_final", fmt_num(lr_pos / 100),
    "--lr_scale", fmt_num(lr_scale),
    "--lr_opac", fmt_num(lr_opac),
    "--lr_sh", fmt_num(lr_sh),
    "--lr_rot", fmt_num(lr_rot),
    "--densify_interval", sprintf("%.0f", as.numeric(densify_interval)),
    "--densify_until", sprintf("%.0f", as.numeric(densify_until)),
    "--opacity_reset_interval", sprintf("%.0f", as.numeric(opacity_reset_interval)),
    "--grad_threshold", fmt_num(grad_threshold),
    "--percent_dense", fmt_num(scale_threshold),
    "--min_opacity", fmt_num(min_opacity),
    "--max_scale", fmt_num(max_scale),
    "--max_gaussians", sprintf("%.0f", as.numeric(max_gaussians)),
    "--downsample", fmt_num(downsample),
    "--sh_degree", sprintf("%d", as.integer(sh_degree)),
    "--lambda_dssim", fmt_num(lambda_dssim),
    "--white_background", if (isTRUE(white_background)) "1" else "0",
    "--seed", sprintf("%d", as.integer(seed))
  )
  if (!is.null(densify_from)) {
    py_args <- c(py_args, "--densify_from", sprintf("%.0f", as.numeric(densify_from)))
  }

  # Run Python trainer and stream live progress directly to console
  status <- system2(py_exec, args = py_args, stdout = if (isTRUE(verbose)) "" else FALSE, stderr = "")

  if (status != 0 || !file.exists(splat_out_file)) {
    cli::cli_abort(c(
      "GPU Gaussian Splatting training failed with exit code {status}.",
      "i" = "The trainer needs Python with {.pkg torch} (CUDA build), {.pkg cupy} and {.pkg pillow}.",
      "i" = "Use {.code device = \"cpu\"} to fall back to the native C++ optimizer."
    ))
  }

  trained_splat <- .gs_read_splat_binary(splat_out_file)
  if (file.exists(stats_json_file)) {
    st <- tryCatch(jsonlite::fromJSON(stats_json_file), error = function(e) NULL)
    if (!is.null(st)) {
      trained_splat$train_history <- as.data.frame(st$history)
      trained_splat$scene_extent <- st$extent
      if (!is.null(st$history$psnr) && length(st$history$psnr) > 0) {
        trained_splat$train_psnr <- utils::tail(st$history$psnr, 1)
      }
    }
  }
  class(trained_splat) <- "gaussian_splat"
  trained_splat
}

# Resolve undistorted PINHOLE cameras/images for a COLMAP reconstruction.
# Gaussian splatting assumes an ideal pinhole projection, so lens distortion
# estimated by OPENCV / RADIAL models must be removed from the photos first
# (COLMAP image_undistorter), exactly as in the reference 3DGS pipeline.
.gs_colmap_undistorted <- function(sfm, verbose = TRUE) {
  recon_dir <- sfm$recon_dir
  if (is.null(recon_dir) || !dir.exists(recon_dir)) return(NULL)
  cams <- tryCatch(gs_import_cameras_colmap(recon_dir), error = function(e) NULL)
  if (is.null(cams) || length(cams) == 0) return(NULL)
  models <- unique(vapply(cams, function(cm) cm$model %||% "PINHOLE", character(1)))
  if (all(models %in% c("PINHOLE", "SIMPLE_PINHOLE"))) {
    return(list(cameras = cams, images = sfm$images, undistorted = FALSE))
  }

  ws <- dirname(sfm$sparse_dir %||% dirname(recon_dir))
  und_dir <- file.path(ws, "undistorted", basename(recon_dir))
  und_sparse <- file.path(und_dir, "sparse")
  und_images <- file.path(und_dir, "images")
  stamp <- file.path(und_dir, "pliman_undistort.stamp")
  src_stamp <- as.character(file.mtime(file.path(recon_dir, "images.txt")))

  up_to_date <- file.exists(file.path(und_sparse, "cameras.txt")) &&
    file.exists(stamp) && identical(readLines(stamp, warn = FALSE)[1], src_stamp)

  if (!up_to_date) {
    colmap_path <- find_colmap()
    if (is.null(colmap_path)) colmap_path <- install_colmap()
    if (isTRUE(verbose)) {
      cli::cli_alert_info("Removing lens distortion ({.val {paste(models, collapse = ', ')}}) with COLMAP image_undistorter...")
    }
    if (dir.exists(und_dir)) unlink(und_dir, recursive = TRUE)
    dir.create(und_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(und_dir, "undistort.log")
    st <- system2(colmap_path, args = c(
      "image_undistorter",
      "--image_path", shQuote(sfm$images),
      "--input_path", shQuote(recon_dir),
      "--output_path", shQuote(und_dir),
      "--output_type", "COLMAP"
    ), stdout = log_file, stderr = log_file)
    if (st != 0 || !dir.exists(und_sparse)) {
      cli::cli_alert_warning("COLMAP image_undistorter failed (see {.path {log_file}}); training on distorted photos.")
      return(list(cameras = cams, images = sfm$images, undistorted = FALSE))
    }
    system2(colmap_path, args = c(
      "model_converter",
      "--input_path", shQuote(und_sparse),
      "--output_path", shQuote(und_sparse),
      "--output_type", "TXT"
    ), stdout = log_file, stderr = log_file)
    writeLines(src_stamp, stamp)
  } else if (isTRUE(verbose)) {
    cli::cli_alert_info("Re-using undistorted photos in {.path {und_dir}}.")
  }

  und_cams <- gs_import_cameras_colmap(und_sparse)
  list(cameras = und_cams, images = normalizePath(und_images, winslash = "/"), undistorted = TRUE)
}

#' @title Optimize 3D Gaussian Splats from Multi-View Photos
#' @name gs_fit
#' @description
#' Optimizes 3D Gaussian Splatting parameters (positions, anisotropic scales,
#' rotations, opacities, and spherical-harmonics colors) directly from a
#' collection of calibrated multi-view photos and camera poses (COLMAP or
#' Agisoft Metashape), following Kerbl et al. (2023).
#'
#' Each iteration renders one training view with a differentiable tile
#' rasterizer, computes the photometric loss
#' \eqn{(1-\lambda) L_1 + \lambda (1 - \mathrm{SSIM})} against the photo,
#' backpropagates exact analytical gradients and applies Adam. Adaptive density
#' control periodically clones small and splits large Gaussians in
#' under-reconstructed regions (high screen-space positional gradient), prunes
#' transparent or oversized ones and resets opacities to remove floaters.
#'
#' When a `colmap_sfm` object (or a splat derived from one, e.g. via
#' [gs_filter_outliers()]) is supplied and the camera model contains lens
#' distortion (e.g. `"OPENCV"`), the photos are first undistorted with COLMAP's
#' `image_undistorter`, because Gaussian splatting assumes an ideal pinhole
#' camera. Point clouds that have not been optimized yet are re-initialized as
#' in the reference implementation (opacity 0.1, isotropic scale equal to the
#' RMS distance to the 3 nearest neighbours).
#'
#' Two back ends are available: a CUDA back end (Python + PyTorch + CuPy,
#' runtime-compiled rasterization kernels) and a native multithreaded C++
#' (OpenMP) back end. The CPU back end optimizes base colors only (SH degree 0),
#' uses an \eqn{L_1} loss and is intended for small scenes. Speed-ups of the GPU
#' back end depend on the available NVIDIA hardware; low-memory GPUs may need a
#' larger `downsample` or a lower `max_gaussians`. See the section
#' *Back ends and requirements (Python)* below.
#'
#' Distance-like thresholds (`scale_threshold`, `max_scale`) and the position
#' learning rate are expressed relative to the scene extent (1.1 x the radius of
#' the camera rig), so the same defaults work for any reconstruction scale.
#'
#' @section Back ends and requirements (Python):
#' `gs_fit()` has two training back ends with different requirements and
#' results:
#'
#' | `device` | Engine | Requires Python? | Result |
#' |---|---|---|---|
#' | `"gpu"` | `inst/python/gs_trainer.py` (PyTorch + CuPy CUDA kernels) | **Yes** | Full Kerbl et al. (2023) method: SH up to degree 3 (view-dependent color), L1 + D-SSIM loss. Recommended for faithful 3D models. |
#' | `"cpu"` | Native C++ (OpenMP), bundled with pliman | No | Base colors only (SH degree 0), L1 loss, much slower. For small scenes, tests or machines without an NVIDIA GPU. |
#' | `"auto"` (default) | GPU when the requirements below are met, otherwise CPU | Only for the GPU | Falls back to the CPU back end (with a console message explaining why) if Python/PyTorch/CUDA is not available. |
#'
#' **The GPU back end runs an external Python process.** The script
#' `gs_trainer.py` ships with pliman (`system.file("python", "gs_trainer.py",
#' package = "pliman")`); R writes the images, cameras and initial splats to
#' temporary files, calls Python and reads the optimized splats back. To use
#' it you need:
#' - an NVIDIA GPU with an up-to-date driver;
#' - Python 3 available in the `PATH` (`python`, `python3` or `py`) or in a
#'   standard installation folder;
#' - PyTorch built with CUDA, e.g.
#'   `pip install torch --index-url https://download.pytorch.org/whl/cu124`;
#' - CuPy matching the CUDA major version of the driver, e.g.
#'   `pip install cupy-cuda12x`;
#' - `numpy` and `Pillow` (installed together with PyTorch in most setups).
#'
#' No CUDA Toolkit, `nvcc` or C++ compiler is needed: CuPy compiles the
#' rasterization kernels at run time (NVRTC). The `gsplat` Python package is
#' **not** required. Use `device = "gpu"` to make sure the GPU back end is used
#' (an informative error is raised if it is not available) and
#' `device = "cpu"` to run without Python at all.
#'
#' Speed depends on the hardware: as a reference, 7,000 iterations on 171
#' photos at 1000 x 560 px (`downsample = 2`) took about 15 min on a GTX 1650
#' (4 GB). GPUs with little memory may need a larger `downsample` or a lower
#' `max_gaussians`.
#'
#' @param images Character path to directory with photos, character vector of image file paths,
#'   a list of loaded `Image` objects, a `colmap_sfm` object, or a `gaussian_splat`
#'   carrying cameras and images (e.g. from [gs_filter_outliers()]).
#' @param cameras Named list of camera specifications (from [gs_import_cameras_colmap()]).
#' @param splat Optional initial `gaussian_splat` object (e.g., from [gs_import_metashape()]
#'   sparse tie points). If `NULL`, initialized automatically from camera bounding volume.
#' @param iterations Integer. Total number of optimization steps. Defaults to `7000`
#'   (the "7k" quality checkpoint of the original paper; `30000` gives the full-quality result).
#' @param lr_pos Double. Initial learning rate for 3D Gaussian positions, relative to the scene
#'   extent. Decays exponentially to `lr_pos / 100`. Defaults to `0.00016`.
#' @param lr_scale Double. Learning rate for log-scales. Defaults to `0.005`.
#' @param lr_opac Double. Learning rate for opacity logits. Defaults to `0.05`.
#' @param lr_sh Double. Learning rate for base colors (SH band 0). Higher bands use `lr_sh / 20`.
#'   Defaults to `0.0025`.
#' @param densify_interval Integer. Iteration frequency for adaptive densification and pruning. Defaults to `100`.
#' @param densify_until Integer. Iteration after which densification stops. Defaults to `round(iterations * 0.5)`.
#' @param grad_threshold Double. Average screen-space (NDC) positional gradient that triggers
#'   split/clone. Defaults to `0.0002`. Lower values produce more Gaussians and finer detail.
#' @param scale_threshold Double. Fraction of the scene extent separating clone (smaller) from
#'   split (larger) Gaussians ("percent dense"). Defaults to `0.01`.
#' @param min_opacity Double. Minimum opacity below which Gaussians are pruned. Defaults to `0.005`.
#' @param max_scale Double. Gaussians larger than `max_scale` x scene extent (or covering more
#'   than 20 pixels after the first opacity reset) are pruned. Defaults to `0.1`.
#' @param max_gaussians Integer. Maximum number of Gaussians allowed during adaptive densification.
#'   Defaults to `500000`.
#' @param downsample Numeric. Downsample factor for training images (e.g. `1`, `2`, `4`)
#'   to accelerate optimization. Defaults to `4`.
#' @param device Character. Computing device for optimization: \code{"auto"} (default: automatically uses
#'   dedicated NVIDIA GPU via CUDA if available, otherwise CPU), \code{"gpu"} (forces GPU acceleration,
#'   aborting if CUDA is unavailable), or \code{"cpu"} (native C++ OpenMP).
#'   **The GPU back end requires Python with PyTorch (CUDA) and CuPy**; the CPU
#'   back end needs no Python. See the section *Back ends and requirements (Python)*.
#' @param verbose Logical. If `TRUE`, displays progress and training metrics. Defaults to `TRUE`.
#' @param sh_degree Integer (0-3). Maximum spherical-harmonics degree for view-dependent color.
#'   One band is activated every 1000 iterations. Defaults to `3` (GPU) - the CPU back end uses `0`.
#' @param lr_rot Double. Learning rate for rotation quaternions. Defaults to `0.001`.
#' @param lambda_dssim Double. Weight of the D-SSIM term in the loss. Defaults to `0.2`.
#' @param densify_from Integer. First iteration of densification. If `NULL` (default),
#'   `min(500, iterations / 10)`.
#' @param opacity_reset_interval Integer. Opacities are reset to 0.01 every this many iterations
#'   during densification (removes floaters). `0` disables. Defaults to `3000`.
#' @param undistort Logical. If `TRUE` (default), photos of COLMAP reconstructions with lens
#'   distortion are undistorted (COLMAP `image_undistorter`) before training.
#' @param reinit Character. Re-initialization of point-cloud Gaussians before training:
#'   `"auto"` (default; only when the splat has not been optimized yet), `"always"`, or `"never"`.
#' @param white_background Logical. Composite renders over white instead of black. Defaults to `FALSE`.
#' @param seed Integer. Random seed for view ordering and densification. Defaults to `42`.
#' @param as_job Logical. If `TRUE` and running inside RStudio, the optimization runs as an
#'   RStudio background job; the trained splat is assigned to `splat_trained` in the global
#'   environment and saved to `gs_fit_result.rds` in the working directory. Defaults to `FALSE`.
#' @return A trained `gaussian_splat` object (with `$cameras`, `$images`, `$train_history`).
#' @references Kerbl, B., Kopanas, G., Leimkuehler, T., & Drettakis, G. (2023). 3D Gaussian
#'   splatting for real-time radiance field rendering. ACM Transactions on Graphics, 42(4), 1-14.
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # 1. Camera poses + sparse tie points with COLMAP
#' sfm <- colmap_sfm("photos/")
#' sfm_clean <- gs_filter_outliers(sfm, k = 5, std_ratio = 2)
#'
#' # 2. Train 3D Gaussian Splats on GPU (CUDA)
#' splat <- gs_fit(sfm_clean, iterations = 7000, downsample = 2)
#'
#' # 3. Interactive WebGL viewer with true anisotropic splat rendering
#' gs_viewer(splat)
#' gs_export_ply(splat, "plant_splat.ply")
gs_fit <- function(images,
                   cameras = NULL,
                   splat = NULL,
                   iterations = 7000,
                   lr_pos = 0.00016,
                   lr_scale = 0.005,
                   lr_opac = 0.05,
                   lr_sh = 0.0025,
                   densify_interval = 100,
                   densify_until = round(iterations * 0.5),
                   grad_threshold = 0.0002,
                   scale_threshold = 0.01,
                   min_opacity = 0.005,
                   max_scale = 0.1,
                   max_gaussians = 500000,
                   downsample = 4,
                   device = c("auto", "gpu", "cpu"),
                   verbose = TRUE,
                   sh_degree = 3,
                   lr_rot = 0.001,
                   lambda_dssim = 0.2,
                   densify_from = NULL,
                   opacity_reset_interval = 3000,
                   undistort = TRUE,
                   reinit = c("auto", "always", "never"),
                   white_background = FALSE,
                   seed = 42,
                   as_job = FALSE) {
  device <- match.arg(device)
  reinit <- match.arg(reinit)
  force(densify_until)

  # ---- Optional RStudio background job ----
  if (isTRUE(as_job)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        rstudioapi::hasFun("jobRunScript")) {
      job_dir <- file.path(tempdir(), "gs_jobs")
      dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      args_rds <- file.path(job_dir, paste0("gs_fit_args_", stamp, ".rds"))
      script <- file.path(job_dir, paste0("gs_fit_job_", stamp, ".R"))
      out_rds <- file.path(getwd(), "gs_fit_result.rds")
      job_args <- list(
        images = images, cameras = cameras, splat = splat, iterations = iterations,
        lr_pos = lr_pos, lr_scale = lr_scale, lr_opac = lr_opac, lr_sh = lr_sh,
        densify_interval = densify_interval, densify_until = densify_until,
        grad_threshold = grad_threshold, scale_threshold = scale_threshold,
        min_opacity = min_opacity, max_scale = max_scale, max_gaussians = max_gaussians,
        downsample = downsample, device = device, verbose = verbose, sh_degree = sh_degree,
        lr_rot = lr_rot, lambda_dssim = lambda_dssim, densify_from = densify_from,
        opacity_reset_interval = opacity_reset_interval, undistort = undistort,
        reinit = reinit, white_background = white_background, seed = seed, as_job = FALSE
      )
      saveRDS(list(args = job_args, wd = getwd(), out = out_rds), args_rds)
      writeLines(c(
        "# Auto-generated RStudio background job by pliman::gs_fit()",
        "library(pliman)",
        paste0("job_info <- readRDS(", deparse(args_rds), ")"),
        "setwd(job_info$wd)",
        "splat_trained <- do.call(pliman::gs_fit, job_info$args)",
        "saveRDS(splat_trained, job_info$out)",
        "cat('\\n[OK] Gaussian splatting finished. Result saved to', job_info$out, '\\n')"
      ), script)
      job_id <- rstudioapi::jobRunScript(path = script, name = "Gaussian Splatting (gs_fit)",
                                         workingDir = getwd(), importEnv = FALSE,
                                         exportEnv = "R_GlobalEnv")
      cli::cli_alert_success("Gaussian splatting dispatched to RStudio's {.strong Jobs} tab.")
      cli::cli_alert_info("When finished, the result is available as {.code splat_trained} and in {.path {out_rds}}.")
      return(invisible(job_id))
    }
    cli::cli_alert_warning("{.arg as_job = TRUE} requires RStudio; running in the current session.")
  }

  # 0. Intelligent parameter resolution for colmap_sfm and gaussian_splat objects
  sfm_src <- NULL
  if (inherits(images, "colmap_sfm")) {
    sfm_src <- images
    col_data <- gs_import_colmap(sfm_src$recon_dir)
    if (is.null(cameras)) cameras <- col_data$cameras
    if (is.null(splat)) splat <- col_data$splat
    images <- sfm_src$images
  } else if (inherits(images, "gaussian_splat")) {
    splat_in <- images
    if (is.null(splat)) splat <- splat_in
    if (is.null(cameras)) cameras <- splat_in$cameras %||% attr(splat_in, "cameras")
    sfm_src <- splat_in$sfm %||% attr(splat_in, "sfm")
    images <- splat_in$images %||% attr(splat_in, "images")
  }

  if (inherits(splat, "colmap_sfm")) {
    sfm_src <- splat
    col_data <- gs_import_colmap(sfm_src$recon_dir)
    if (is.null(cameras)) cameras <- col_data$cameras
    splat <- col_data$splat
    if (missing(images) || is.null(images)) images <- sfm_src$images
  } else if (inherits(splat, "gaussian_splat")) {
    if (is.null(cameras)) cameras <- splat$cameras %||% attr(splat, "cameras")
    if (is.null(sfm_src)) sfm_src <- splat$sfm %||% attr(splat, "sfm")
    if (missing(images) || is.null(images)) images <- splat$images %||% attr(splat, "images")
  }

  if (inherits(cameras, "colmap_sfm")) {
    if (is.null(sfm_src)) sfm_src <- cameras
    cameras <- gs_import_cameras_colmap(cameras$recon_dir)
  }

  # Lens undistortion: Gaussian splatting assumes an ideal pinhole camera
  if (isTRUE(undistort) && !is.null(sfm_src) && inherits(sfm_src, "colmap_sfm")) {
    und <- tryCatch(.gs_colmap_undistorted(sfm_src, verbose = verbose), error = function(e) {
      cli::cli_alert_warning("Undistortion skipped: {conditionMessage(e)}")
      NULL
    })
    if (!is.null(und) && isTRUE(und$undistorted)) {
      keep_lbl <- names(cameras)
      cameras <- und$cameras
      if (!is.null(keep_lbl)) cameras <- cameras[intersect(names(cameras), keep_lbl)]
      images <- und$images
      if (isTRUE(verbose)) cli::cli_alert_success("Training on undistorted PINHOLE photos ({length(cameras)} views).")
    }
  } else if (isTRUE(verbose) && !is.null(cameras) && length(cameras) > 0) {
    mdl <- unique(vapply(cameras, function(cm) cm$model %||% "PINHOLE", character(1)))
    if (!all(mdl %in% c("PINHOLE", "SIMPLE_PINHOLE"))) {
      cli::cli_alert_warning("Camera model {.val {mdl}} has lens distortion but no COLMAP workspace is available to undistort the photos; residual distortion will blur the reconstruction.")
    }
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
  img_base_lc <- tolower(basename(img_paths))
  img_sans_lc <- tolower(tools::file_path_sans_ext(basename(img_paths)))

  for (i in seq_along(cameras)) {
    c_spec <- cameras[[i]]
    c_label <- if (!is.null(c_spec$label)) c_spec$label else cam_names[i]

    matched_img <- NULL
    if (is.list(images) && !is.character(images)) {
      if (!is.null(images[[c_label]])) {
        matched_img <- images[[c_label]]
      } else if (i <= length(images)) {
        matched_img <- images[[i]]
      }
    } else {
      c_base <- tolower(basename(c_label))
      c_sans <- tolower(tools::file_path_sans_ext(c_base))
      j <- match(c_base, img_base_lc)
      if (is.na(j)) j <- match(c_sans, img_sans_lc)
      if (!is.na(j)) matched_img <- img_paths[j]
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

  extent <- .gs_scene_extent(lapply(dataset, `[[`, "cam"))

  # 2. Initialize Gaussian Splats if not provided
  if (is.null(splat)) {
    if (isTRUE(verbose)) cli::cli_alert_info("No point cloud supplied: initializing 3D Gaussians randomly inside the camera rig...")
    cam_pos <- do.call(rbind, lapply(dataset, function(d) d$cam$position))
    center <- colMeans(cam_pos)
    n_init <- 100000
    pts <- matrix(stats::runif(n_init * 3, -1, 1) * extent * 0.6, ncol = 3)
    pts <- sweep(pts, 2, center, "+")
    splat <- gs_create(pts, colors = matrix(stats::runif(n_init * 3), ncol = 3), initial_opacity = 0.1)
  } else {
    if (!inherits(splat, "gaussian_splat")) {
      stop("'splat' must be an object of class 'gaussian_splat'.")
    }
  }

  # 2.1 Re-initialize raw point clouds as in the reference implementation
  untrained <- !isTRUE(splat$trained) &&
    all(abs(splat$quaternions[, 1] - 1) < 1e-8) &&
    length(unique(round(splat$opacities_logit, 6))) == 1L
  if (reinit == "always" || (reinit == "auto" && untrained)) {
    if (isTRUE(verbose)) cli::cli_alert_info("Initializing {format(splat$num_gaussians, big.mark = ',')} Gaussians from the point cloud (opacity 0.1, k-NN scales).")
    splat$scales_log <- matrix(log(.gs_knn_scales(splat$positions)), ncol = 3, nrow = splat$num_gaussians)
    splat$quaternions <- cbind(1, matrix(0, splat$num_gaussians, 3))
    splat$opacities_logit <- rep(log(0.1 / 0.9), splat$num_gaussians)
  }

  # 2.5 Resolve computing device and branch to GPU backend if applicable
  gpu_info <- if (device == "cpu") list(has_gpu = FALSE) else .gs_detect_cuda_gpu()
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

  finalize <- function(sp) {
    sp$cameras <- cameras
    sp$images <- if (is.character(images)) images else NULL
    sp$sfm <- sfm_src
    sp$trained <- TRUE
    sp$faces <- integer(0)
    sp$num_faces <- 0L
    sp$has_faces <- FALSE
    sp$has_colors <- TRUE
    sp$scene_extent <- sp$scene_extent %||% extent
    attr(sp, "cameras") <- cameras
    attr(sp, "images") <- sp$images
    class(sp) <- "gaussian_splat"
    sp
  }

  if (use_gpu) {
    if (isTRUE(verbose)) {
      cli::cli_rule("Optimizing 3D Gaussian Splats (GPU CUDA Acceleration)")
      cli::cli_alert_success("NVIDIA GPU detected: {.strong {gpu_info$name}} (CUDA backend)")
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
      lr_rot = lr_rot,
      densify_interval = densify_interval,
      densify_from = densify_from,
      densify_until = densify_until,
      opacity_reset_interval = opacity_reset_interval,
      grad_threshold = grad_threshold,
      scale_threshold = scale_threshold,
      min_opacity = min_opacity,
      max_scale = max_scale,
      max_gaussians = max_gaussians,
      downsample = downsample,
      sh_degree = sh_degree,
      lambda_dssim = lambda_dssim,
      white_background = white_background,
      seed = seed,
      verbose = verbose
    )
    return(finalize(splat_trained))
  }

  # ===========================================================================
  # 3. Native C++ (OpenMP) back end
  # ===========================================================================
  set.seed(seed)
  if (device == "auto" && isTRUE(verbose)) {
    cli::cli_alert_warning("GPU back end not available, using the native C++ (CPU) back end. {gpu_info$reason %||% ''}")
    cli::cli_bullets(c(
      "i" = "The GPU back end needs Python with PyTorch (CUDA) and CuPy; see the section {.emph Back ends and requirements (Python)} in {.code ?gs_fit}.",
      "i" = "Use {.code device = \"cpu\"} to silence this message."
    ))
  }
  if (ncol(splat$sh_coeffs) > 3 && isTRUE(verbose)) {
    cli::cli_alert_info("CPU back end optimizes base colors only (SH degree 0).")
  }
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Preparing training views (downsample factor: {downsample}x)...")
  }

  bg_rgb <- if (isTRUE(white_background)) c(1, 1, 1) else c(0, 0, 0)
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
      new_w <- max(16L, as.integer(round(dim(img_obj)[1] / downsample)))
      new_h <- max(16L, as.integer(round(dim(img_obj)[2] / downsample)))
      img_obj <- image_resize(img_obj, width = new_w, height = new_h)
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
    if (length(dim(raw_img)) == 3 && dim(raw_img)[3] > 3) raw_img <- raw_img[, , 1:3]
    if (length(dim(raw_img)) == 2) raw_img <- array(rep(raw_img, 3), dim = c(dim(raw_img), 3))
    img_num <- as.numeric(raw_img)
    if (is.raw(raw_img) || max(img_num, na.rm = TRUE) > 1.0) {
      img_num <- img_num / 255.0
    }

    training_views[[idx]] <- list(
      cam = cam_scaled,
      target = img_num,          # W x H x 3, column-major (x fastest)
      label = d$label
    )
  }

  num_views <- length(training_views)
  if (is.null(densify_from)) densify_from <- min(500, max(1, round(iterations / 10)))
  if (is.null(opacity_reset_interval) || opacity_reset_interval <= 0) opacity_reset_interval <- .Machine$integer.max

  # Keep only the base color for the CPU optimizer
  splat$sh_coeffs <- splat$sh_coeffs[, 1:3, drop = FALSE]
  splat$sh_degree <- 0L

  init_adam <- function(n_gauss) {
    z3 <- matrix(0.0, nrow = n_gauss, ncol = 3)
    list(
      m = list(pos = z3, scale = z3, opac = numeric(n_gauss), sh = z3,
               rot = matrix(0.0, nrow = n_gauss, ncol = 4)),
      v = list(pos = z3, scale = z3, opac = numeric(n_gauss), sh = z3,
               rot = matrix(0.0, nrow = n_gauss, ncol = 4))
    )
  }

  moments <- init_adam(splat$num_gaussians)
  adam_m <- moments$m
  adam_v <- moments$v

  accum_grad_2d <- numeric(splat$num_gaussians)
  accum_denom <- numeric(splat$num_gaussians)
  recent_losses <- numeric(0)
  lr_pos_abs <- lr_pos * extent
  lr_pos_final <- lr_pos_abs / 100
  history <- list()

  if (isTRUE(verbose)) {
    cli::cli_rule("Optimizing 3D Gaussian Splats (native C++ / OpenMP)")
    cli::cli_alert_info("Initial Gaussians: {format(splat$num_gaussians, big.mark = ',')} | Total Steps: {iterations} | Views: {num_views} | Scene extent: {signif(extent, 4)}")
  }

  start_time <- Sys.time()
  view_order <- sample.int(num_views)
  vptr <- 1L
  log_interval <- max(1, min(25, round(iterations / 50)))

  for (step in seq_len(iterations)) {
    if (vptr > num_views) {
      view_order <- sample.int(num_views)
      vptr <- 1L
    }
    view <- training_views[[view_order[vptr]]]
    vptr <- vptr + 1L

    tt <- (step - 1) / max(1, iterations - 1)
    lr_now <- exp(log(lr_pos_abs) * (1 - tt) + log(lr_pos_final) * tt)

    step_res <- cpp_gs_train_step(
      camera = view$cam,
      gaussians = splat,
      target_img = view$target,
      adam_m = adam_m,
      adam_v = adam_v,
      iter = as.integer(step),
      lr_pos = as.numeric(lr_now),
      lr_scale = as.numeric(lr_scale),
      lr_opac = as.numeric(lr_opac),
      lr_sh = as.numeric(lr_sh),
      lr_rot = as.numeric(lr_rot),
      bg_r = bg_rgb[1], bg_g = bg_rgb[2], bg_b = bg_rgb[3]
    )

    splat <- step_res$gaussians
    class(splat) <- "gaussian_splat"
    adam_m <- step_res$adam_m
    adam_v <- step_res$adam_v

    if (step < densify_until) {
      g_norm <- step_res$grad_2d_norm
      vis <- step_res$visible
      if (length(g_norm) == length(accum_grad_2d)) {
        accum_grad_2d <- accum_grad_2d + g_norm
        accum_denom <- accum_denom + as.numeric(vis)
      }
    }

    recent_losses <- c(recent_losses, step_res$loss)
    if (length(recent_losses) > 20) recent_losses <- utils::tail(recent_losses, 20)

    # Periodic adaptive densification and pruning
    if (step < densify_until && step > densify_from && step %% densify_interval == 0) {
      mean_grads <- accum_grad_2d / pmax(1, accum_denom)
      n_before <- splat$num_gaussians
      splat <- gs_densify(
        splat = splat,
        grad_2d_norm = mean_grads,
        grad_threshold = grad_threshold,
        scale_threshold = scale_threshold * extent,
        min_opacity = min_opacity,
        max_scale = max_scale * extent,
        max_gaussians = max_gaussians,
        verbose = FALSE
      )

      moments <- init_adam(splat$num_gaussians)
      adam_m <- moments$m
      adam_v <- moments$v
      accum_grad_2d <- numeric(splat$num_gaussians)
      accum_denom <- numeric(splat$num_gaussians)

      if (isTRUE(verbose)) {
        cli::cli_alert_info("Step {step}/{iterations} | Densify: {format(n_before, big.mark = ',')} -> {format(splat$num_gaussians, big.mark = ',')} Gaussians (Loss: {round(mean(recent_losses), 4)})")
      }
    }
    if (step < densify_until && step %% opacity_reset_interval == 0) {
      splat$opacities_logit <- pmin(splat$opacities_logit, log(0.01 / 0.99))
      adam_m$opac[] <- 0
      adam_v$opac[] <- 0
      if (isTRUE(verbose)) cli::cli_alert_info("Step {step}/{iterations} | Opacity reset")
    }

    if (step %% log_interval == 0 || step == iterations) {
      curr_time <- Sys.time()
      elapsed_s <- as.numeric(difftime(curr_time, start_time, units = "secs"))
      it_per_sec <- step / max(0.001, elapsed_s)
      remaining_s <- (iterations - step) / max(0.001, it_per_sec)
      eta_str <- if (remaining_s >= 60) paste0(round(remaining_s / 60, 1), "m") else paste0(round(remaining_s), "s")
      pct <- round((step / iterations) * 100)
      history[[length(history) + 1]] <- data.frame(step = step, loss = mean(recent_losses), n = splat$num_gaussians)
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Step {step}/{iterations} [{pct}%] | {round(it_per_sec, 1)} it/s | ETA: {eta_str} | Loss: {round(mean(recent_losses), 4)} | Gaussians: {format(splat$num_gaussians, big.mark = ',')}")
      }
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  if (isTRUE(verbose)) {
    cli::cli_alert_success("Training finished in {elapsed}s! Final Gaussians: {format(splat$num_gaussians, big.mark = ',')}")
    cli::cli_rule()
  }

  splat$train_history <- do.call(rbind, history)
  splat$scene_extent <- extent
  finalize(splat)
}

# Voxel-grid parameters shared by gs_filter_outliers(method = "density") and the
# "Online Noise Filter" of gs_viewer(), so both produce exactly the same result.
# The cell size follows the viewer rule: 1.2% of the robust subject dimension,
# clamped to [0.008, 0.035]; the grid origin is the coordinate-wise median.
.gs_density_grid_params <- function(pos, cell_size = NULL) {
  ctr <- apply(pos, 2, stats::median, na.rm = TRUE)
  d <- sqrt(rowSums(sweep(pos, 2, ctr)^2))
  q70 <- stats::quantile(d, probs = 0.70, na.rm = TRUE, names = FALSE)
  if (!is.finite(q70) || q70 < 0.05) q70 <- 0.5
  q_min <- apply(pos, 2, stats::quantile, probs = 0.02, na.rm = TRUE)
  q_max <- apply(pos, 2, stats::quantile, probs = 0.98, na.rm = TRUE)
  global_dim <- max(pmax(0.05, q_max - q_min), 0.1)
  target_dim <- min(global_dim, max(0.2, q70 * 2.2))
  if (is.null(cell_size)) cell_size <- max(0.008, min(0.035, target_dim * 0.012))
  list(origin = as.numeric(ctr), cell_size = as.numeric(cell_size), target_dim = target_dim)
}

# Spatial density per point: number of points falling in the 3 x 3 x 3 block of
# voxels centred on the point's voxel (identical to computeSpatialDensities() in
# the viewer).
.gs_voxel_density <- function(pos, origin, cell_size) {
  ijk <- floor(sweep(pos, 2, origin) / cell_size)
  ijk[!is.finite(ijk)] <- 0
  mn <- apply(ijk, 2, min)
  ijk <- sweep(ijk, 2, mn) + 1 # >= 1, leaves room for the -1 neighbour
  dims <- apply(ijk, 2, max) + 2
  if (prod(dims) < 2^52) {
    key <- ijk[, 1] + dims[1] * (ijk[, 2] + dims[2] * ijk[, 3])
    ukeys <- unique(key)
    cnt <- tabulate(match(key, ukeys), nbins = length(ukeys))
    dsum <- numeric(length(ukeys))
    for (dx in -1:1) for (dy in -1:1) for (dz in -1:1) {
      m <- match(ukeys + dx + dims[1] * (dy + dims[2] * dz), ukeys)
      ok <- !is.na(m)
      dsum[ok] <- dsum[ok] + cnt[m[ok]]
    }
    return(dsum[match(key, ukeys)])
  }
  # Fallback for extremely sparse/huge grids: string keys
  key <- paste(ijk[, 1], ijk[, 2], ijk[, 3])
  ukeys <- unique(key)
  cnt <- tabulate(match(key, ukeys), nbins = length(ukeys))
  uijk <- ijk[match(ukeys, key), , drop = FALSE]
  dsum <- numeric(length(ukeys))
  for (dx in -1:1) for (dy in -1:1) for (dz in -1:1) {
    m <- match(paste(uijk[, 1] + dx, uijk[, 2] + dy, uijk[, 3] + dz), ukeys)
    ok <- !is.na(m)
    dsum[ok] <- dsum[ok] + cnt[m[ok]]
  }
  dsum[match(key, ukeys)]
}

#' @title Statistical and Radius Outlier Removal (SOR / ROR) for 3D Point Clouds
#' @name gs_filter_outliers
#' @description
#' Filters out noisy, floating tie-points and airborne outliers from point clouds
#' or Gaussian Splats based on local neighborhood density using a high-performance native
#' C++ 3D Kd-Tree (OpenMP accelerated, zero external package dependencies).
#'
#' Two filters are available (`method`):
#' - `"density"` (default): the same voxel-grid density filter used by the
#'   *Online Noise Filter* panel of [gs_viewer()]. Space is divided into cubic
#'   voxels (`cell_size`, by default 1.2% of the robust subject size, clamped to
#'   0.008-0.035 scene units), the density of each point is the number of points
#'   in the 3 x 3 x 3 voxel block around it, and the points whose density falls
#'   below the `density` percentile are removed. Because the cut-off is a
#'   percentile, it is not affected by a few extremely distant points, which is
#'   what makes the classical SOR filter too permissive on SfM clouds. The
#'   presets `"light"`, `"medium"` and `"strong"` correspond to the viewer
#'   buttons (15, 35 and 65%).
#' - `"sor"`: Statistical Outlier Removal (mean distance to the `k` nearest
#'   neighbours, cut-off `mean + std_ratio * sd`), optionally combined with the
#'   radius filter (`radius`, `min_neighbors`).
#' - `"both"`: a point must pass both filters.
#'
#' @param splat A `gaussian_splat`, `colmap_sfm`, or numeric matrix of 3D coordinates.
#' @param k Integer. Number of nearest neighbors to evaluate (e.g. 16 to 30). Defaults to `20`.
#'   Used by `method = "sor"` / `"both"`.
#' @param std_ratio Numeric. Standard deviation multiplier threshold. Smaller values (e.g. 0.8 to 1.2)
#'   filter more aggressively, removing scattered floating points around the subject. Defaults to `1.0`.
#'   Used by `method = "sor"` / `"both"`.
#' @param radius Numeric. Optional sphere radius for density filtering. If > 0, points with fewer than
#'   `min_neighbors` within `radius` are additionally filtered out. Defaults to `0.0` (disabled).
#' @param min_neighbors Integer. Minimum neighbors required within `radius`. Defaults to `0` (disabled).
#' @param verbose Logical. Show progress and statistics. Defaults to `TRUE`.
#' @param method Character. `"density"` (default, same filter as [gs_viewer()]),
#'   `"sor"` or `"both"`. See Description.
#' @param density Filter strength for `method = "density"`: a percentile between
#'   0 and 100 (points below it are removed) or one of the viewer presets
#'   `"light"` (15), `"medium"` (35, default) or `"strong"` (65).
#' @param cell_size Numeric. Voxel edge length (scene units) for
#'   `method = "density"`. `NULL` (default) uses the same automatic rule as
#'   [gs_viewer()].
#' @return An updated and cleaned `gaussian_splat`, `colmap_sfm`, or point matrix.
#'   The attribute `"densities"` holds the voxel density of the kept points
#'   (density filter) and `"density_dists"` the mean k-NN distance (SOR filter).
#' @export
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Clean floating noise from a COLMAP reconstruction (same as the viewer's
#' # "Medium" noise filter)
#' sfm_clean <- gs_filter_outliers(sfm_completo, density = "medium")
#'
#' # Classical statistical outlier removal
#' sfm_sor <- gs_filter_outliers(sfm_completo, k = 20, std_ratio = 1.0, method = "sor")
#' plot(sfm_clean)
gs_filter_outliers <- function(splat, k = 20, std_ratio = 1.0, radius = 0.0, min_neighbors = 0, verbose = TRUE,
                               method = c("density", "sor", "both"), density = "medium",
                               cell_size = NULL) {
  method <- match.arg(method)
  if (is.character(density)) {
    presets <- c(off = 0, light = 15, medium = 35, strong = 65)
    density <- presets[match.arg(tolower(density), names(presets))]
  }
  density <- as.numeric(density)
  if (length(density) != 1 || !is.finite(density) || density < 0 || density > 100) {
    cli::cli_abort("{.arg density} must be a percentile between 0 and 100 or one of {.val light}, {.val medium}, {.val strong}.")
  }
  if (!is.null(cell_size) && (!is.numeric(cell_size) || cell_size <= 0)) {
    cli::cli_abort("{.arg cell_size} must be a positive number or {.code NULL}.")
  }
  is_colmap <- inherits(splat, "colmap_sfm")
  cameras_saved <- NULL
  images_saved <- NULL

  sfm_saved <- NULL
  if (is_colmap) {
    images_saved <- splat$images
    sfm_saved <- splat
    col_data <- gs_import_colmap(splat$recon_dir)
    cameras_saved <- col_data$cameras
    splat <- col_data$splat
  } else if (inherits(splat, "gaussian_splat")) {
    cameras_saved <- splat$cameras %||% attr(splat, "cameras")
    images_saved <- splat$images %||% attr(splat, "images")
    sfm_saved <- splat$sfm %||% attr(splat, "sfm")
  }

  is_splat <- inherits(splat, "gaussian_splat")
  pts <- if (is_splat) splat$positions else as.matrix(splat)

  n <- nrow(pts)
  if (n <= k && method != "density") return(splat)
  if (n < 2) return(splat)

  keep <- rep(TRUE, n)
  mean_dists <- NULL
  dens <- NULL

  if (method %in% c("density", "both")) {
    gp <- .gs_density_grid_params(pts, cell_size)
    dens <- .gs_voxel_density(pts, gp$origin, gp$cell_size)
    keep_d <- rep(TRUE, n)
    if (density > 0) {
      # Same percentile indexing as the viewer: sorted[floor(frac * n)] (0-based)
      p_idx <- min(n - 1, max(0, floor(density / 100 * n)))
      cutoff_d <- sort(dens)[p_idx + 1]
      keep_d <- dens >= cutoff_d
    } else {
      cutoff_d <- 0
    }
    keep <- keep & keep_d
    if (isTRUE(verbose)) {
      n_k <- sum(keep_d)
      cli::cli_alert_success("Voxel Density Filter ({density}% percentile, cell = {signif(gp$cell_size, 3)}, cut-off = {cutoff_d} points): kept {format(n_k, big.mark = ',')} points ({round(n_k / n * 100, 1)}%), removed {format(n - n_k, big.mark = ',')} low-density points.")
    }
  }

  if (method %in% c("sor", "both") && n > k) {
    filter_res <- cpp_gs_filter_outliers(
      points = pts,
      k = as.integer(k),
      std_ratio = as.numeric(std_ratio),
      radius = as.numeric(radius),
      min_neighbors = as.integer(min_neighbors)
    )
    mean_dists <- filter_res$mean_dists
    keep <- keep & filter_res$keep
    if (isTRUE(verbose)) {
      n_k <- filter_res$num_kept
      pct <- round((n_k / n) * 100, 1)
      cli::cli_alert_success("Statistical Outlier Filter (SOR): kept {format(n_k, big.mark = ',')} points ({pct}%), removed {format(filter_res$num_removed, big.mark = ',')} floating noise points.")
    }
  }

  n_kept <- sum(keep)
  if (isTRUE(verbose) && method == "both") {
    cli::cli_alert_info("Combined filters: kept {format(n_kept, big.mark = ',')} of {format(n, big.mark = ',')} points ({round(n_kept / n * 100, 1)}%).")
  }

  if (is_splat) {
    splat <- .gs_subset(splat, keep)
    if (!is.null(sfm_saved)) {
      splat$sfm <- sfm_saved
      attr(splat, "sfm") <- sfm_saved
    }
    if (!is.null(mean_dists)) attr(splat, "density_dists") <- mean_dists[keep]
    if (!is.null(dens)) attr(splat, "densities") <- dens[keep]

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
    if (!is.null(mean_dists)) attr(res_mat, "density_dists") <- mean_dists[keep]
    if (!is.null(dens)) attr(res_mat, "densities") <- dens[keep]
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
#' @examplesIf interactive()
#' library(pliman)
#'
#' # Crop around plant subject
#' splat_cropped <- gs_crop(sfm_completo,
#'                          xlim = c(-0.5, 0.5),
#'                          ylim = c(-0.5, 0.5),
#'                          zlim = c(-0.2, 0.8))
#' plot(splat_cropped)
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

  sfm_saved <- NULL
  if (is_colmap) {
    images_saved <- splat$images
    sfm_saved <- splat
    col_data <- gs_import_colmap(splat$recon_dir)
    cameras_saved <- col_data$cameras
    splat <- col_data$splat
  } else if (inherits(splat, "gaussian_splat")) {
    cameras_saved <- splat$cameras %||% attr(splat, "cameras")
    images_saved <- splat$images %||% attr(splat, "images")
    sfm_saved <- splat$sfm %||% attr(splat, "sfm")
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
    splat <- .gs_subset(splat, keep)
    if (!is.null(sfm_saved)) {
      splat$sfm <- sfm_saved
      attr(splat, "sfm") <- sfm_saved
    }

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


