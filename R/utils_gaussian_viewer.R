# ==============================================================================
# 3D GAUSSIAN SPLATTING & POINT CLOUD INTERACTIVE VIEWER FOR R (pliman)
# High-performance Three.js / WebGL standalone desktop and browser visualization
# With Metashape-style Trackball Navigation Sphere, 3-Axis Gizmo & Surface Closure
# And COLMAP / Metashape 3D Camera Frustums with Interactive Photo Inspection
# Author: Tiago Olivoto & pliman team
# ==============================================================================

#' @title Interactive 3D Point Cloud and Gaussian Splat Viewer
#' @name gs_viewer
#' @description
#' Launches a standalone, GPU-accelerated interactive 3D WebGL viewer for both
#' **Gaussian Splats** and **Agisoft Metashape / COLMAP point clouds** with zero external
#' server dependencies.
#'
#' Features include:
#' - **Online Real-Time Noise Filter (SOR / Density Filter)**: Prunes spatial outlier noise live in the browser using a high-speed 3D spatial voxel grid.
#' - **3D Measurement Ruler (Tape Measure)**: Interactive point-to-point Euclidean 3D distance and component dimension measurement.
#' - **Camera Projection Mode Toggle**: Instant switching between Perspective and Orthographic CAD views.
#' - **RealityCapture Reconstruction Region (3D Clipping Box)**: Interactive 6-plane section clipping to isolate plant stems and foliage.
#' - **Direct Client-Side Export**: Export clean filtered point clouds (.ply binary / .txt) and calibrated camera poses (.txt COLMAP) directly to download.
#' - **Multi-Engine Shading Modes**: Photorealistic RGB, Studio Clay (geometry surface inspection), Elevation Heatmap, Point Density Map, and SAM Classes.
#' - **Zero Server Dependencies**: 100% client-side HTML5 WebGL running locally in Google Chrome / Edge without Shiny, Node.js or Python.
#' - **Direct File Import & Drag-and-Drop**: Load `.ply`, `.splat`, or COLMAP `points3D.txt` files directly into the 3D canvas anytime!
#' - **COLMAP / Metashape 3D Camera Frustums**: Pyramidal camera view cones and magenta sight lines in 3D world space.
#' - **Interactive Photo Inspection**: Click on any camera frustum in 3D to open a floating photo viewer modal with high-res thumbnails, camera metadata, and photo cycling!
#' - **"Look Through Camera" View**: Align the 3D scene camera directly with any photo pose to inspect point cloud alignment from that exact shot's perspective.
#' - **Metashape Navigation Trackball**: Translucent smooth glass 3D sphere (customizable opacity, zero wireframe clutter) with colored Pitch (X), Roll (Y), and Yaw (Z) gimbal rings.
#' - **Direct 3D Axis Dragging**: Click and drag directly on the X, Y, or Z rings in the 3D viewport to rotate the model.
#' - **3-Axis Model Alignment**: Interactive sliders, step buttons (+/-90&deg;, +/-15&deg;), and code copy to level tilted models horizontally.
#' - **Double-Click Pivot Focusing**: Double click anywhere on the model to set the orbit pivot and center the Trackball Sphere.
#' - **Dual Rendering Engines**: Switch between crisp Metashape-style point pixels and volumetric continuous Gaussian splats.
#' - **Splat Surface Closure**: Continuous alpha-blending slider to seal gaps and form solid watertight surfaces.
#' - **Subject Cluster Framing**: Automatically frames the camera tightly on the main object rather than distant room walls.
#' - **Real-Time 3D Coordinates**: Dynamic coordinate tracking in world space on mouse hover.
#'
#' @param splat Optional. An object of class `gaussian_splat` (e.g. from [gs_import_metashape()] or [gs_fit()]),
#'   `colmap_sfm`, or a character string path to a `.ply`, `.splat`, or `.txt` file. If `NULL` (default), launches
#'   a standalone viewer in import mode with a drag-and-drop interface.
#' @param cameras Optional camera poses. Can be a list of camera objects (e.g. from [gs_import_cameras_colmap()]
#'   or [gs_import_cameras_metashape()]), or a character string path to a COLMAP sparse directory
#'   (e.g. `"sparse/0"`), `images.txt`, or Metashape `cameras.xml`.
#' @param images Optional path to the directory containing the original photos (e.g. `"my_photos/"`). If provided,
#'   clicking on any 3D camera frustum will display that photo in a floating viewer modal!
#' @param show_cameras Logical. If `TRUE` (default when cameras are provided), displays 3D camera frustums.
#' @param sight_lines Logical. If `TRUE` (default), displays magenta sight rays connecting each camera to the subject center.
#' @param camera_size Numeric scalar. Initial scale factor for 3D camera frustums. Defaults to `1.0`.
#' @param file Character string. Path to the output HTML file. Defaults to a temporary file.
#' @param browser Logical. If `TRUE` (default), automatically opens the default web browser.
#' @param browser_path Optional character string specifying a custom browser executable (e.g. Google Chrome).
#' @param thumbnails Character string specifying thumbnail handling: `"on_demand"`, `"embedded"`, or `"none"` (defaults to `"on_demand"`).
#' @param mode Character string. Initial render mode: `"points"`, `"splats"`, or `"mesh"`. Defaults to `"mesh"` if a mesh is loaded, otherwise `"points"`.
#' @param color_by Character string. Initial color mode: `"rgb"`, `"clay"`, `"height"`, `"density"`, or `"semantic"`.
#' @param point_size Numeric scalar. Initial point/splat size. Defaults to `3.5`.
#' @param texture Optional character string path to a texture image (.png or .jpg) or texture atlas directory.
#' @param trackball Logical. If `TRUE` (default), displays the Metashape-style Trackball Navigation Sphere and 3-axis rotation gizmo.
#' @param sphere_opacity Numeric scalar between 0 and 1. Initial transparency of the navigation sphere. Defaults to `0.04` (very subtle, clean glass).
#' @param rotation Numeric vector of length 3 `c(x, y, z)`. Initial rotation in degrees. Defaults to `c(0, 0, 0)`.
#' @param bg Character string. Background theme: `"dark"` (default, slate dark) or `"light"` (Metashape light).
#' @param title Character string. Viewer window title.
#' @return Invisible character string containing the path to the generated HTML file.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # 1. Launch standalone viewer ready for drag & drop
#' gs_viewer()
#'
#' # 2. Or pass a model and cameras directly
#' splat <- gs_import_metashape("cuia_densa.ply")
#' gs_viewer(splat, cameras = "colmap/sparse/0", images = "photos/")
#' }
gs_viewer <- function(splat = NULL,
                      cameras = NULL,
                      images = NULL,
                      show_cameras = TRUE,
                      sight_lines = TRUE,
                      camera_size = 1.0,
                      texture = NULL,
                      thumbnails = c("on_demand", "embedded", "none"),
                      file = NULL,
                      browser = TRUE,
                      browser_path = NULL,
                      mode = c("points", "splats", "mesh"),
                      color_by = c("rgb", "clay", "height", "density", "semantic"),
                      point_size = 3.5,
                      trackball = TRUE,
                      sphere_opacity = 0.04,
                      rotation = c(0, 0, 0),
                      bg = c("dark", "light"),
                      title = "3D Point Cloud & Gaussian Splat Viewer [pliman]") {
  mode_user_specified <- !missing(mode)
  color_by_user_specified <- !missing(color_by)

  # Normalize thumbnails argument
  if (isTRUE(thumbnails)) {
    thumbnails <- "on_demand"
  } else if (isFALSE(thumbnails) || is.null(thumbnails)) {
    thumbnails <- "none"
  } else if (is.character(thumbnails)) {
    thumbnails <- match.arg(thumbnails, c("on_demand", "embedded", "none"))
  } else {
    thumbnails <- "on_demand"
  }

  has_initial_model <- !is.null(splat)
  orig_splat_textured_dir <- NULL

  if (has_initial_model) {
    if (inherits(splat, "colmap_mesh")) {
      if (is.null(cameras) && !is.null(splat$cameras)) cameras <- splat$cameras
      if (is.null(images) && !is.null(splat$images)) images <- splat$images
      orig_splat_textured_dir <- splat$textured_dir
      mesh_ply <- splat$mesh_ply
      if (file.exists(mesh_ply)) splat <- gs_import_metashape(mesh_ply)
    } else if (inherits(splat, "colmap_dense")) {
      if (is.null(cameras) && !is.null(splat$cameras)) cameras <- splat$cameras
      if (is.null(images) && !is.null(splat$images)) images <- splat$images
      fused_ply <- splat$fused_ply
      if (!is.null(fused_ply) && file.exists(fused_ply)) splat <- gs_import_metashape(fused_ply)
    } else if (inherits(splat, "colmap_sfm")) {
      if (is.null(images) && !is.null(splat$images)) images <- splat$images
      col_data <- gs_import_colmap(splat$recon_dir)
      if (is.null(cameras)) cameras <- col_data$cameras
      splat <- col_data$splat
    } else if (is.character(splat) && file.exists(splat)) {
      splat <- gs_import_metashape(splat)
    }

    if (!inherits(splat, "gaussian_splat")) {
      stop("'splat' must be an object of class 'gaussian_splat', 'colmap_sfm', 'colmap_dense', 'colmap_mesh', or a path to a valid PLY/TXT point cloud file.")
    }

    n <- splat$num_gaussians
    has_faces <- !is.null(splat$faces) && length(splat$faces) > 0
    num_faces <- if (has_faces) as.integer(length(splat$faces) / 3) else 0L
    has_colors <- if (!is.null(splat$has_colors)) isTRUE(splat$has_colors) else TRUE
    faces_str <- if (has_faces) paste(splat$faces, collapse = ",") else ""

    # Compute biometrics for HUD
    measures <- gs_measure(splat)
    bbox_min <- attr(measures, "bbox_min")
    bbox_max <- attr(measures, "bbox_max")

    # Robust centroid and bounding box:
    pos <- splat$positions
    center_x <- stats::median(pos[, 1], na.rm = TRUE)
    center_y <- stats::median(pos[, 2], na.rm = TRUE)
    center_z <- stats::median(pos[, 3], na.rm = TRUE)

    # Distance of all points from median center
    dist_from_med <- sqrt((pos[, 1] - center_x)^2 + (pos[, 2] - center_y)^2 + (pos[, 3] - center_z)^2)
    q70_dist <- stats::quantile(dist_from_med, probs = 0.70, na.rm = TRUE)
    if (!is.finite(q70_dist) || q70_dist < 0.05) q70_dist <- 0.5

    # Quantiles (2% - 98%)
    q_min <- apply(pos, 2, stats::quantile, probs = 0.02, na.rm = TRUE)
    q_max <- apply(pos, 2, stats::quantile, probs = 0.98, na.rm = TRUE)
    rob_dx <- max(0.05, q_max[1] - q_min[1])
    rob_dy <- max(0.05, q_max[2] - q_min[2])
    rob_dz <- max(0.05, q_max[3] - q_min[3])
    global_dim <- max(rob_dx, rob_dy, rob_dz, 0.1)

    # Frame camera tightly onto the central subject rather than distant room walls:
    target_dim <- min(global_dim, max(0.2, q70_dist * 2.2))

    pos_centered <- cbind(
      pos[, 1] - center_x,
      pos[, 2] - center_y,
      pos[, 3] - center_z
    )

    # Extract RGB colors from Spherical Harmonics degree 0:
    sh_c0 <- 0.28209479177387814
    r_col <- pmin(1.0, pmax(0.0, 0.5 + sh_c0 * splat$sh_coeffs[, 1]))
    g_col <- pmin(1.0, pmax(0.0, 0.5 + sh_c0 * splat$sh_coeffs[, 2]))
    b_col <- pmin(1.0, pmax(0.0, 0.5 + sh_c0 * splat$sh_coeffs[, 3]))

    # Opacities
    opacs <- 1.0 / (1.0 + exp(-splat$opacities_logit))

    # Scales (geometric mean of x, y, z scales)
    mean_scales <- exp((splat$scales_log[, 1] + splat$scales_log[, 2] + splat$scales_log[, 3]) / 3.0)

    # Semantics
    semantics <- as.integer(splat$semantic_ids)
    unique_sems <- sort(unique(semantics))

    pos_str <- paste(round(as.vector(t(pos_centered)), 5), collapse = ",")
    rgb_str <- paste(round(as.vector(t(cbind(r_col, g_col, b_col))), 4), collapse = ",")
    opac_str <- paste(round(opacs, 4), collapse = ",")
    scale_str <- paste(round(mean_scales, 5), collapse = ",")
    sem_str <- paste(semantics, collapse = ",")
  } else {
    n <- 0
    has_faces <- FALSE
    num_faces <- 0L
    has_colors <- TRUE
    faces_str <- ""
    center_x <- 0; center_y <- 0; center_z <- 0
    rob_dx <- 0; rob_dy <- 0; rob_dz <- 0
    target_dim <- 1.0
    pos_str <- ""
    rgb_str <- ""
    opac_str <- ""
    scale_str <- ""
    sem_str <- ""
    unique_sems <- 0
  }

  if (!mode_user_specified) {
    mode <- if (has_faces) "mesh" else "points"
  } else {
    mode <- match.arg(mode, c("points", "splats", "mesh"))
  }

  if (!color_by_user_specified) {
    color_by <- if (has_faces && !has_colors) "clay" else "rgb"
  } else {
    color_by <- match.arg(color_by, c("rgb", "clay", "height", "density", "semantic"))
  }

  # Texture resolution
  tex_data_uri <- ""
  if (!is.null(texture)) {
    tex_path <- NULL
    if (is.character(texture) && length(texture) == 1) {
      if (file.exists(texture) && !dir.exists(texture)) {
        tex_path <- texture
      } else if (dir.exists(texture)) {
        candidates <- list.files(texture, pattern = "\\.(png|jpg|jpeg)$", full.names = TRUE, ignore.case = TRUE)
        if (length(candidates) > 0) tex_path <- candidates[1]
      }
    }
    if (!is.null(tex_path) && file.exists(tex_path)) {
      raw_tex <- readBin(tex_path, what = "raw", n = file.info(tex_path)$size)
      mime <- if (grepl("\\.png$", tex_path, ignore.case = TRUE)) "image/png" else "image/jpeg"
      if (requireNamespace("base64enc", quietly = TRUE)) {
        tex_data_uri <- paste0("data:", mime, ";base64,", base64enc::base64encode(raw_tex))
      }
    }
  } else if (has_initial_model && !is.null(orig_splat_textured_dir)) {
    candidates <- list.files(orig_splat_textured_dir, pattern = "\\.(png|jpg|jpeg)$", full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) > 0 && file.exists(candidates[1])) {
      raw_tex <- readBin(candidates[1], what = "raw", n = file.info(candidates[1])$size)
      mime <- if (grepl("\\.png$", candidates[1], ignore.case = TRUE)) "image/png" else "image/jpeg"
      if (requireNamespace("base64enc", quietly = TRUE)) {
        tex_data_uri <- paste0("data:", mime, ";base64,", base64enc::base64encode(raw_tex))
      }
    }
  }

  bg_hex <- if (is.character(bg) && length(bg) == 1 && grepl("^#", bg)) {
    bg
  } else {
    bg_mode <- match.arg(bg, c("dark", "light"))
    if (bg_mode == "light") "#f1f5f9" else "#090d16"
  }

  # Initial rotation
  init_rot_x <- if (length(rotation) >= 1) as.numeric(rotation[1]) else 0
  init_rot_y <- if (length(rotation) >= 2) as.numeric(rotation[2]) else 0
  init_rot_z <- if (length(rotation) >= 3) as.numeric(rotation[3]) else 0
  init_sphere_opac <- max(0.0, min(0.5, as.numeric(sphere_opacity)))

  # --------------------------------------------------------------------------
  # 1. Resolve Cameras and Photos
  # --------------------------------------------------------------------------
  parsed_cameras <- NULL
  if (!is.null(cameras)) {
    if (is.character(cameras) && length(cameras) == 1) {
      if (dir.exists(cameras)) {
        if (file.exists(file.path(cameras, "cameras.txt")) || file.exists(file.path(cameras, "images.txt"))) {
          parsed_cameras <- gs_import_cameras_colmap(cameras)
        } else if (dir.exists(file.path(cameras, "sparse", "1"))) {
          parsed_cameras <- gs_import_cameras_colmap(file.path(cameras, "sparse", "1"))
        } else if (dir.exists(file.path(cameras, "sparse", "0"))) {
          parsed_cameras <- gs_import_cameras_colmap(file.path(cameras, "sparse", "0"))
        } else {
          parsed_cameras <- gs_import_cameras_colmap(cameras)
        }
      } else if (file.exists(cameras)) {
        ext <- tolower(tools::file_ext(cameras))
        if (ext == "xml") {
          parsed_cameras <- gs_import_cameras_metashape(cameras)
        } else {
          parsed_cameras <- gs_import_cameras_colmap(cameras)
        }
      }
    } else if (is.list(cameras)) {
      if (!is.null(cameras[[1]]$position) && !is.null(cameras[[1]]$R)) {
        parsed_cameras <- cameras
      } else if (!is.null(cameras$cameras)) {
        parsed_cameras <- cameras$cameras
      }
    }
  } else {
    # Check all possible locations where cameras could be stored in splat
    if (!is.null(attr(splat, "cameras"))) {
      parsed_cameras <- attr(splat, "cameras")
    } else if (!is.null(attr(splat, "colmap_data")$cameras)) {
      parsed_cameras <- attr(splat, "colmap_data")$cameras
    } else if (is.list(splat) && !is.null(splat$cameras)) {
      parsed_cameras <- splat$cameras
    } else if (!is.null(attr(splat, "sfm")$recon_dir) && dir.exists(attr(splat, "sfm")$recon_dir)) {
      tryCatch({
        cd <- gs_import_colmap(attr(splat, "sfm")$recon_dir)
        parsed_cameras <- cd$cameras
      }, error = function(e) NULL)
    }
  }

  # Resolve images folder if not explicitly passed
  if (is.null(images)) {
    if (!is.null(attr(splat, "images"))) {
      images <- attr(splat, "images")
    } else if (is.list(splat) && !is.null(splat$images)) {
      images <- splat$images
    } else if (!is.null(attr(splat, "sfm")$images)) {
      images <- attr(splat, "sfm")$images
    } else if (!is.null(attr(splat, "colmap_data")$images)) {
      images <- attr(splat, "colmap_data")$images
    }
  }

  # Resolve images folder if cameras are present
  resolved_images_dir <- NULL
  if (!is.null(parsed_cameras) && length(parsed_cameras) > 0) {
    if (!is.null(images) && is.character(images) && length(images) == 1 && dir.exists(images)) {
      resolved_images_dir <- normalizePath(images, winslash = "/")
    } else {
      # Search plausible candidate directories
      cand_dirs <- c(
        if (!is.null(images) && is.character(images)) images,
        if (!is.null(cameras) && is.character(cameras) && length(cameras) == 1) c(
          dirname(cameras),
          dirname(dirname(cameras)),
          dirname(dirname(dirname(cameras))),
          file.path(dirname(cameras), "images")
        ),
        file.path(getwd(), "frames"),
        file.path(getwd(), "images"),
        file.path(getwd(), "photos"),
        getwd()
      )
      first_label <- parsed_cameras[[1]]$label
      for (cd in cand_dirs) {
        if (dir.exists(cd)) {
          if (file.exists(file.path(cd, first_label)) || file.exists(file.path(cd, basename(first_label)))) {
            resolved_images_dir <- normalizePath(cd, winslash = "/")
            break
          }
        }
      }
    }
  }

  # Determine destination HTML file
  if (is.null(file) || !nzchar(file)) {
    if (!is.null(resolved_images_dir) && dir.exists(resolved_images_dir) && file.access(resolved_images_dir, 2) == 0) {
      file <- file.path(resolved_images_dir, "gs_viewer.html")
    } else {
      file <- tempfile(pattern = "gs_viewer_", fileext = ".html")
    }
  }

  file_dir <- normalizePath(dirname(file), winslash = "/", mustWork = FALSE)

  # Relative path helper
  get_relative_path <- function(target, base_dir) {
    norm_target <- normalizePath(target, winslash = "/", mustWork = FALSE)
    norm_base   <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
    if (norm_target == norm_base) return(".")
    if (substr(norm_target, 1, 2) != substr(norm_base, 1, 2)) return(basename(norm_target))
    target_parts <- strsplit(norm_target, "/")[[1]]
    base_parts   <- strsplit(norm_base, "/")[[1]]
    common_len <- 0
    min_len <- min(length(target_parts), length(base_parts))
    while (common_len < min_len && target_parts[common_len + 1] == base_parts[common_len + 1]) {
      common_len <- common_len + 1
    }
    ups <- length(base_parts) - common_len
    up_str <- if (ups > 0) paste(rep("..", ups), collapse = "/") else ""
    down_str <- if (common_len < length(target_parts)) paste(target_parts[(common_len + 1):length(target_parts)], collapse = "/") else ""
    if (nzchar(up_str) && nzchar(down_str)) paste0(up_str, "/", down_str)
    else if (nzchar(up_str)) up_str
    else down_str
  }

  # Process cameras into JSON
  cameras_json <- "[]"
  has_cameras <- !is.null(parsed_cameras) && length(parsed_cameras) > 0

  if (has_cameras) {
    n_cams <- length(parsed_cameras)
    cli::cli_alert_info("Prepared {n_cams} calibrated camera frustums for 3D inspection.")
    has_img_files <- thumbnails != "none" && !is.null(resolved_images_dir) && dir.exists(resolved_images_dir)

    pb <- NULL
    if (has_img_files) {
      if (thumbnails == "on_demand") {
        cli::cli_alert_success("Linked {n_cams} camera photos for instant on-demand inspection (renders upon click).")
      } else if (thumbnails == "embedded") {
        cli::cli_alert_info("Embedding photo thumbnails (Base64) from {.path {resolved_images_dir}}...")
        pb <- cli::cli_progress_bar("Encoding photo thumbnails", total = n_cams)
      }
    }

    cams_list <- vector("list", n_cams)
    for (i in seq_len(n_cams)) {
      cam <- parsed_cameras[[i]]
      lbl <- if (!is.null(cam$label)) cam$label else paste0("Camera_", i)
      c_pos <- as.numeric(cam$position)
      c_centered <- c(c_pos[1] - center_x, c_pos[2] - center_y, c_pos[3] - center_z)
      r_mat <- as.numeric(t(cam$R)) # Flattened row-major

      img_ref <- ""
      if (has_img_files) {
        img_fpath <- file.path(resolved_images_dir, lbl)
        if (!file.exists(img_fpath)) {
          img_fpath <- file.path(resolved_images_dir, basename(lbl))
        }
        if (file.exists(img_fpath)) {
          if (thumbnails == "on_demand") {
            img_ref <- get_relative_path(img_fpath, file_dir)
          } else if (thumbnails == "embedded") {
            if (!is.null(pb)) cli::cli_progress_update(id = pb)
            tryCatch({
              im <- image_import(img_fpath)
              d <- dim(im)[1:2]
              total_px <- d[1] * d[2]
              if (total_px > 50000) {
                scale_factor <- sqrt(50000 / total_px)
                tw <- max(10, as.integer(round(d[1] * scale_factor)))
                th <- max(10, as.integer(round(d[2] * scale_factor)))
                im <- image_resize(im, width = tw, height = th, filter = "nearest", plot = FALSE)
              }
              tf <- tempfile(fileext = ".jpg")
              .write_image_file(im, tf, quality = 50)
              img_ref <- paste0("data:image/jpeg;base64,", base64enc::base64encode(tf))
              unlink(tf)
            }, error = function(e) NULL)
          }
        }
      }

      cams_list[[i]] <- list(
        index = i,
        label = lbl,
        pos = round(c_centered, 5),
        R = round(r_mat, 5),
        fx = round(if (!is.null(cam$fx)) cam$fx else 1000, 2),
        fy = round(if (!is.null(cam$fy)) cam$fy else 1000, 2),
        w = as.integer(if (!is.null(cam$width)) cam$width else 1920),
        h = as.integer(if (!is.null(cam$height)) cam$height else 1080),
        img = img_ref
      )
    }

    if (requireNamespace("jsonlite", quietly = TRUE)) {
      cameras_json <- jsonlite::toJSON(cams_list, auto_unbox = TRUE)
    }
  }

  # Semantic labels mapping
  sem_labels_json <- c()
  for (s_id in unique_sems) {
    s_name <- if (s_id == 0) "Unclassified" else if (s_id == 1) "Leaf (SAM)" else if (s_id == 2) "Stem / Branch" else paste0("Class ", s_id)
    sem_labels_json <- c(sem_labels_json, sprintf('"%d": "%s"', s_id, s_name))
  }
  sem_map_js <- paste0("{", paste(sem_labels_json, collapse = ", "), "}")

  # HTML template with embedded CSS and modern Three.js WebGL application
  html_content <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>', title, '</title>
  <style>
    :root {
      --bg-dark: ', bg_hex, ';
      --card-bg: rgba(15, 23, 42, 0.88);
      --card-border: rgba(255, 255, 255, 0.12);
      --accent: #10b981;
      --accent-glow: rgba(16, 185, 129, 0.35);
      --text: #f8fafc;
      --text-muted: #94a3b8;
    }
    * {
      box-sizing: border-box;
      margin: 0;
      padding: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      user-select: none;
    }
    body, html {
      width: 100%;
      height: 100%;
      overflow: hidden;
      background-color: var(--bg-dark);
      color: var(--text);
    }
    #canvas-container {
      width: 100%;
      height: 100%;
      position: absolute;
      top: 0;
      left: 0;
    }

    /* Floating Glassmorphism HUD Panel */
    .hud-panel {
      position: absolute;
      top: 16px;
      left: 16px;
      width: 360px;
      max-height: calc(100vh - 32px);
      overflow-y: auto;
      background: var(--card-bg);
      backdrop-filter: blur(16px);
      -webkit-backdrop-filter: blur(16px);
      border: 1px solid var(--card-border);
      border-radius: 14px;
      box-shadow: 0 12px 36px rgba(0, 0, 0, 0.5), 0 0 0 1px rgba(255, 255, 255, 0.05);
      padding: 18px;
      display: flex;
      flex-direction: column;
      gap: 13px;
      z-index: 100;
      scrollbar-width: thin;
      scrollbar-color: rgba(255, 255, 255, 0.2) transparent;
    }
    .hud-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
      border-bottom: 1px solid rgba(255, 255, 255, 0.08);
      padding-bottom: 10px;
    }
    .hud-title {
      font-size: 15px;
      font-weight: 700;
      letter-spacing: -0.3px;
      color: #ffffff;
      display: flex;
      align-items: center;
      gap: 8px;
    }
    .hud-badge {
      font-size: 10px;
      font-weight: 700;
      text-transform: uppercase;
      padding: 3px 7px;
      border-radius: 6px;
      background: linear-gradient(135deg, #10b981, #059669);
      color: #ffffff;
      box-shadow: 0 2px 6px var(--accent-glow);
    }
    .fps-badge {
      font-size: 11px;
      font-family: monospace;
      color: #10b981;
      font-weight: bold;
    }

    /* Section & Controls */
    .control-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .control-label {
      font-size: 11px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      color: var(--text-muted);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .control-label span.val {
      color: #ffffff;
      font-family: monospace;
    }

    /* Buttons and segmented pills */
    .segmented-control {
      display: grid;
      grid-template-columns: repeat(2, 1fr);
      gap: 4px;
      background: rgba(0, 0, 0, 0.35);
      padding: 3px;
      border-radius: 8px;
      border: 1px solid rgba(255, 255, 255, 0.05);
    }
    .segmented-control-3 {
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 4px;
      background: rgba(0, 0, 0, 0.35);
      padding: 3px;
      border-radius: 8px;
      border: 1px solid rgba(255, 255, 255, 0.05);
    }
    .segmented-control-4 {
      display: grid;
      grid-template-columns: repeat(4, 1fr);
      gap: 3px;
      background: rgba(0, 0, 0, 0.35);
      padding: 3px;
      border-radius: 8px;
      border: 1px solid rgba(255, 255, 255, 0.05);
    }
    .segmented-control-5 {
      display: grid;
      grid-template-columns: repeat(5, 1fr);
      gap: 3px;
      background: rgba(0, 0, 0, 0.35);
      padding: 3px;
      border-radius: 8px;
      border: 1px solid rgba(255, 255, 255, 0.05);
    }
    .ruler-badge {
      background: rgba(15, 23, 42, 0.92);
      border: 1px solid #facc15;
      color: #facc15;
      padding: 4px 8px;
      border-radius: 6px;
      font-size: 11px;
      font-weight: 700;
      font-family: monospace;
      pointer-events: none;
      box-shadow: 0 4px 12px rgba(0, 0, 0, 0.6);
      white-space: nowrap;
    }
    .seg-btn {
      background: transparent;
      border: none;
      color: var(--text-muted);
      padding: 6px 4px;
      font-size: 11px;
      font-weight: 600;
      border-radius: 6px;
      cursor: pointer;
      transition: all 0.15s ease;
      text-align: center;
    }
    .seg-btn.active {
      background: #10b981;
      color: #ffffff;
      box-shadow: 0 2px 8px rgba(16, 185, 129, 0.4);
    }
    .seg-btn:hover:not(.active) {
      color: #ffffff;
      background: rgba(255, 255, 255, 0.08);
    }

    /* Axis rotation rows */
    .axis-control-row {
      display: flex;
      align-items: center;
      gap: 8px;
      margin-top: 2px;
    }
    .axis-badge {
      font-size: 10px;
      font-weight: 700;
      padding: 3px 6px;
      border-radius: 4px;
      color: #fff;
      min-width: 62px;
      text-align: center;
      letter-spacing: 0.3px;
    }
    .axis-x { background: #ef4444; }
    .axis-y { background: #22c55e; }
    .axis-z { background: #3b82f6; }
    .axis-val {
      font-size: 11px;
      font-family: monospace;
      color: #ffffff;
      min-width: 36px;
      text-align: right;
    }
    .btn-reset-sm {
      background: transparent;
      border: 1px solid rgba(255, 255, 255, 0.15);
      color: var(--text-muted);
      padding: 2px 7px;
      font-size: 10px;
      border-radius: 4px;
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .btn-reset-sm:hover {
      color: #fff;
      border-color: #10b981;
    }

    /* Range sliders */
    input[type=range] {
      -webkit-appearance: none;
      width: 100%;
      height: 5px;
      border-radius: 4px;
      background: rgba(255, 255, 255, 0.12);
      outline: none;
      cursor: pointer;
    }
    input[type=range]::-webkit-slider-thumb {
      -webkit-appearance: none;
      width: 15px;
      height: 15px;
      border-radius: 50%;
      background: #10b981;
      box-shadow: 0 0 8px #10b981;
      cursor: pointer;
      transition: transform 0.1s;
    }
    input[type=range]::-webkit-slider-thumb:hover {
      transform: scale(1.2);
    }

    /* Metrics Table HUD */
    .metrics-box {
      background: rgba(0, 0, 0, 0.35);
      border: 1px solid rgba(255, 255, 255, 0.06);
      border-radius: 8px;
      padding: 10px 12px;
      display: flex;
      flex-direction: column;
      gap: 5px;
      font-size: 11px;
    }
    .metric-row {
      display: flex;
      justify-content: space-between;
      color: var(--text-muted);
    }
    .metric-val {
      font-family: monospace;
      font-weight: 600;
      color: #ffffff;
    }

    /* Camera & Action Buttons */
    .cam-grid {
      display: grid;
      grid-template-columns: repeat(4, 1fr);
      gap: 5px;
    }
    .step-grid {
      display: grid;
      grid-template-columns: repeat(4, 1fr);
      gap: 4px;
      margin-top: 4px;
    }
    .btn-action {
      background: rgba(255, 255, 255, 0.06);
      border: 1px solid rgba(255, 255, 255, 0.1);
      color: var(--text);
      padding: 6px 3px;
      font-size: 10px;
      font-weight: 600;
      border-radius: 6px;
      cursor: pointer;
      text-align: center;
      transition: all 0.15s ease;
      white-space: nowrap;
    }
    .btn-action:hover {
      background: rgba(255, 255, 255, 0.15);
      border-color: rgba(255, 255, 255, 0.3);
      color: #ffffff;
    }

    /* Primary CTA */
    .btn-primary {
      background: linear-gradient(135deg, #10b981, #059669);
      border: none;
      color: #ffffff;
      padding: 9px 12px;
      font-size: 12px;
      font-weight: 700;
      border-radius: 8px;
      cursor: pointer;
      box-shadow: 0 4px 12px rgba(16, 185, 129, 0.35);
      transition: all 0.15s ease;
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 6px;
    }
    .btn-primary:hover {
      transform: translateY(-1px);
      box-shadow: 0 6px 16px rgba(16, 185, 129, 0.45);
    }

    /* Toast Notification */
    #toast-msg {
      position: fixed;
      top: 20px;
      right: 20px;
      background: rgba(16, 185, 129, 0.95);
      color: #ffffff;
      padding: 10px 18px;
      border-radius: 8px;
      font-size: 12px;
      font-weight: 600;
      box-shadow: 0 8px 24px rgba(0, 0, 0, 0.4);
      display: none;
      z-index: 1000;
      animation: fadeIn 0.2s ease;
    }

    /* Camera Frustum 3D Hover Tooltip */
    #cam-tooltip {
      position: absolute;
      background: rgba(15, 23, 42, 0.92);
      backdrop-filter: blur(8px);
      border: 1px solid #10b981;
      color: #ffffff;
      padding: 5px 10px;
      border-radius: 6px;
      font-size: 11px;
      font-family: monospace;
      pointer-events: none;
      display: none;
      z-index: 500;
      box-shadow: 0 4px 14px rgba(0, 0, 0, 0.5);
      transform: translate(-50%, -130%);
    }

    /* Interactive Photo Modal Lightbox */
    .modal-backdrop {
      position: fixed;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: rgba(0, 0, 0, 0.75);
      backdrop-filter: blur(10px);
      -webkit-backdrop-filter: blur(10px);
      display: none;
      justify-content: center;
      align-items: center;
      z-index: 2000;
      animation: fadeIn 0.2s ease;
    }
    .modal-card {
      background: #0f172a;
      border: 1px solid rgba(255, 255, 255, 0.15);
      border-radius: 14px;
      box-shadow: 0 24px 60px rgba(0, 0, 0, 0.8);
      max-width: 90vw;
      max-height: 90vh;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .modal-header {
      padding: 12px 18px;
      background: rgba(255, 255, 255, 0.03);
      border-bottom: 1px solid rgba(255, 255, 255, 0.08);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .modal-title-box {
      display: flex;
      align-items: center;
      gap: 10px;
    }
    .modal-cam-name {
      font-size: 14px;
      font-weight: 700;
      color: #ffffff;
      font-family: monospace;
    }
    .modal-cam-meta {
      font-size: 11px;
      color: var(--text-muted);
    }
    .modal-close-btn {
      background: rgba(255, 255, 255, 0.08);
      border: none;
      color: #ffffff;
      font-size: 14px;
      width: 28px;
      height: 28px;
      border-radius: 50%;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: center;
      transition: all 0.15s ease;
    }
    .modal-close-btn:hover {
      background: #ef4444;
    }
    .modal-body {
      padding: 12px;
      display: flex;
      justify-content: center;
      align-items: center;
      background: #020617;
      overflow: hidden;
    }
    .modal-photo {
      max-width: 82vw;
      max-height: 70vh;
      object-fit: contain;
      border-radius: 6px;
      box-shadow: 0 6px 20px rgba(0, 0, 0, 0.6);
    }
    .modal-footer {
      padding: 12px 18px;
      background: rgba(255, 255, 255, 0.03);
      border-top: 1px solid rgba(255, 255, 255, 0.08);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .btn-modal-nav {
      background: rgba(255, 255, 255, 0.08);
      border: 1px solid rgba(255, 255, 255, 0.12);
      color: #ffffff;
      padding: 7px 14px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .btn-modal-nav:hover {
      background: rgba(255, 255, 255, 0.18);
    }
    .btn-modal-action {
      background: linear-gradient(135deg, #10b981, #059669);
      border: none;
      color: #ffffff;
      padding: 7px 16px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 700;
      cursor: pointer;
      box-shadow: 0 2px 8px rgba(16, 185, 129, 0.35);
      transition: all 0.15s ease;
    }
    .btn-modal-action:hover {
      transform: translateY(-1px);
      box-shadow: 0 4px 12px rgba(16, 185, 129, 0.45);
    }

    /* Bottom Info and Real-time Coordinate HUD */
    .bottom-bar {
      position: absolute;
      bottom: 14px;
      right: 18px;
      display: flex;
      flex-direction: column;
      align-items: flex-end;
      gap: 6px;
      z-index: 100;
      pointer-events: none;
    }
    .coords-hud {
      background: var(--card-bg);
      backdrop-filter: blur(16px);
      -webkit-backdrop-filter: blur(16px);
      border: 1px solid var(--card-border);
      border-radius: 8px;
      padding: 8px 14px;
      display: flex;
      flex-direction: column;
      gap: 3px;
      box-shadow: 0 8px 24px rgba(0, 0, 0, 0.4);
      font-family: monospace;
    }
    .coords-title {
      font-size: 10px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      color: #10b981;
      display: flex;
      justify-content: space-between;
      gap: 12px;
    }
    .coords-vals {
      font-size: 13px;
      font-weight: 700;
      color: #ffffff;
      white-space: nowrap;
    }
    .bottom-hint {
      font-size: 10px;
      color: rgba(255, 255, 255, 0.7);
      background: rgba(0, 0, 0, 0.65);
      padding: 5px 12px;
      border-radius: 6px;
      font-family: monospace;
      letter-spacing: 0.2px;
    }

    /* File Drop Zone / Welcome Modal */
    .dropzone-overlay {
      position: absolute;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: rgba(9, 13, 22, 0.78);
      backdrop-filter: blur(14px);
      -webkit-backdrop-filter: blur(14px);
      display: flex;
      justify-content: center;
      align-items: center;
      z-index: 500;
      transition: all 0.25s ease;
    }
    .dropzone-card {
      background: rgba(15, 23, 42, 0.96);
      border: 2px dashed rgba(16, 185, 129, 0.55);
      border-radius: 20px;
      padding: 30px;
      text-align: center;
      max-width: 620px;
      width: 92%;
      box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.8), 0 0 35px rgba(16, 185, 129, 0.25);
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 13px;
    }
    .dropzone-card.drag-active {
      border-color: #10b981;
      background: rgba(16, 185, 129, 0.16);
      transform: scale(1.02);
    }
    .dropzone-icon {
      font-size: 42px;
      margin-bottom: 0;
      animation: floatIcon 3s ease-in-out infinite;
    }
    @keyframes floatIcon {
      0%, 100% { transform: translateY(0); }
      50% { transform: translateY(-7px); }
    }
    .dropzone-title {
      font-size: 19px;
      font-weight: 700;
      color: #f8fafc;
      letter-spacing: -0.3px;
    }
    .dropzone-desc {
      font-size: 12px;
      color: #94a3b8;
      line-height: 1.5;
    }
    .dropzone-actions {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      justify-content: center;
      margin-top: 6px;
    }
    .btn-drop-action {
      background: linear-gradient(135deg, #10b981, #059669);
      color: #ffffff;
      border: none;
      padding: 10px 18px;
      font-size: 13px;
      font-weight: 700;
      border-radius: 9px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 8px;
      box-shadow: 0 4px 14px rgba(16, 185, 129, 0.4);
      transition: all 0.2s ease;
    }
    .btn-drop-action:hover {
      transform: translateY(-2px);
      box-shadow: 0 6px 18px rgba(16, 185, 129, 0.55);
    }
    .btn-drop-secondary {
      background: rgba(255, 255, 255, 0.08);
      color: #e2e8f0;
      border: 1px solid rgba(255, 255, 255, 0.15);
      padding: 10px 16px;
      font-size: 12px;
      font-weight: 600;
      border-radius: 9px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      transition: all 0.2s ease;
    }
    .btn-drop-secondary:hover {
      background: rgba(255, 255, 255, 0.16);
      color: #ffffff;
    }
    .btn-import-header {
      background: rgba(16, 185, 129, 0.18);
      border: 1px solid rgba(16, 185, 129, 0.5);
      color: #34d399;
      font-size: 11px;
      font-weight: 700;
      padding: 4px 10px;
      border-radius: 6px;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 5px;
      transition: all 0.15s ease;
    }
    .btn-import-header:hover {
      background: #10b981;
      color: #ffffff;
      box-shadow: 0 2px 8px rgba(16, 185, 129, 0.45);
    }

    /* Distinct Asset Import Fields */
    .field-card {
      background: rgba(255, 255, 255, 0.04);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 12px;
      padding: 12px 14px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      width: 100%;
      transition: all 0.2s ease;
      text-align: left;
    }
    .field-card:hover {
      border-color: rgba(16, 185, 129, 0.4);
      background: rgba(255, 255, 255, 0.06);
    }
    .field-card.loaded {
      border-color: rgba(16, 185, 129, 0.6);
      background: rgba(16, 185, 129, 0.08);
    }
    .field-header {
      display: flex;
      align-items: center;
      gap: 12px;
      flex: 1;
      min-width: 0;
    }
    .field-icon {
      font-size: 24px;
      flex-shrink: 0;
    }
    .field-title {
      font-size: 12px;
      font-weight: 700;
      color: #f8fafc;
      margin-bottom: 2px;
    }
    .field-status {
      font-size: 11px;
      color: var(--text-muted);
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .field-status.loaded {
      color: #34d399;
      font-weight: 600;
    }
    .field-btn {
      background: rgba(255, 255, 255, 0.08);
      border: 1px solid rgba(255, 255, 255, 0.18);
      color: #f8fafc;
      padding: 7px 12px;
      border-radius: 8px;
      font-size: 11px;
      font-weight: 600;
      cursor: pointer;
      white-space: nowrap;
      transition: all 0.15s ease;
    }
    .field-btn:hover {
      background: #10b981;
      border-color: #10b981;
      color: #ffffff;
      transform: translateY(-1px);
    }
    .btn-batch-primary {
      background: linear-gradient(135deg, #10b981, #059669);
      border: none;
      color: #ffffff;
      padding: 11px 18px;
      border-radius: 10px;
      font-size: 12px;
      font-weight: 700;
      cursor: pointer;
      width: 100%;
      box-shadow: 0 4px 14px rgba(16, 185, 129, 0.4);
      transition: all 0.2s ease;
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
    }
    .btn-batch-primary:hover {
      transform: translateY(-2px);
      box-shadow: 0 6px 20px rgba(16, 185, 129, 0.6);
    }
    .btn-enter-scene {
      background: rgba(255, 255, 255, 0.12);
      border: 1px solid rgba(255, 255, 255, 0.25);
      color: #ffffff;
      padding: 8px 18px;
      border-radius: 8px;
      font-size: 12px;
      font-weight: 700;
      cursor: pointer;
      transition: all 0.2s ease;
    }
    .btn-enter-scene:hover {
      background: rgba(255, 255, 255, 0.22);
      border-color: #ffffff;
    }

    /* Scene Layers in HUD */
    .layer-list {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .layer-item {
      background: rgba(0, 0, 0, 0.28);
      border: 1px solid rgba(255, 255, 255, 0.08);
      border-radius: 8px;
      padding: 7px 10px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 8px;
    }
    .layer-info {
      display: flex;
      align-items: center;
      gap: 8px;
      flex: 1;
      min-width: 0;
    }
    .layer-icon {
      font-size: 16px;
      flex-shrink: 0;
    }
    .layer-texts {
      min-width: 0;
    }
    .layer-name {
      font-size: 11px;
      font-weight: 700;
      color: #ffffff;
      line-height: 1.2;
    }
    .layer-meta {
      font-size: 10px;
      color: var(--text-muted);
      line-height: 1.2;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .layer-actions {
      display: flex;
      align-items: center;
      gap: 5px;
    }
    .btn-layer-toggle {
      background: rgba(255, 255, 255, 0.06);
      border: 1px solid rgba(255, 255, 255, 0.12);
      color: var(--text-muted);
      border-radius: 6px;
      padding: 3px 6px;
      font-size: 12px;
      cursor: pointer;
      transition: all 0.15s ease;
      opacity: 0.45;
    }
    .btn-layer-toggle.active {
      opacity: 1;
      color: #34d399;
      background: rgba(16, 185, 129, 0.2);
      border-color: rgba(16, 185, 129, 0.4);
    }
    .btn-layer-sub {
      background: rgba(255, 255, 255, 0.06);
      border: 1px solid rgba(255, 255, 255, 0.12);
      color: #cbd5e1;
      border-radius: 6px;
      padding: 3px 6px;
      font-size: 11px;
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .btn-layer-sub:hover {
      background: rgba(255, 255, 255, 0.18);
      color: #ffffff;
    }

    /* Metashape Top Floating Toolbar */
    .metashape-toolbar {
      position: absolute;
      top: 14px;
      left: 50%;
      transform: translateX(-50%);
      background: rgba(15, 23, 42, 0.92);
      backdrop-filter: blur(16px);
      -webkit-backdrop-filter: blur(16px);
      border: 1px solid rgba(255, 255, 255, 0.14);
      border-radius: 30px;
      padding: 5px 10px;
      display: flex;
      align-items: center;
      gap: 6px;
      z-index: 200;
      box-shadow: 0 10px 30px rgba(0, 0, 0, 0.55), 0 0 0 1px rgba(255, 255, 255, 0.05);
      max-width: 96vw;
      flex-wrap: wrap;
      justify-content: center;
    }
    .tb-btn {
      background: transparent;
      border: 1px solid transparent;
      color: #94a3b8;
      padding: 6px 13px;
      font-size: 11px;
      font-weight: 700;
      border-radius: 20px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      transition: all 0.18s ease;
      white-space: nowrap;
    }
    .tb-btn:hover {
      color: #ffffff;
      background: rgba(255, 255, 255, 0.09);
    }
    .tb-btn.active {
      background: #10b981;
      border-color: #10b981;
      color: #ffffff;
      box-shadow: 0 2px 10px rgba(16, 185, 129, 0.45);
    }
    .tb-btn.active.tb-cams {
      background: #3b82f6;
      border-color: #3b82f6;
      box-shadow: 0 2px 10px rgba(59, 130, 246, 0.45);
    }
    .tb-btn.active.tb-ruler {
      background: #eab308;
      border-color: #eab308;
      color: #0f172a;
      box-shadow: 0 2px 10px rgba(234, 179, 8, 0.55);
    }
    .tb-btn.active.tb-proj {
      background: #8b5cf6;
      border-color: #8b5cf6;
      color: #ffffff;
      box-shadow: 0 2px 10px rgba(139, 92, 246, 0.55);
    }
    .tb-divider {
      width: 1px;
      height: 18px;
      background: rgba(255, 255, 255, 0.15);
      margin: 0 2px;
    }

    /* Single Stateful Toggle Buttons in Sidebar */
    .btn-toggle {
      background: rgba(255, 255, 255, 0.06);
      border: 1px solid rgba(255, 255, 255, 0.12);
      color: var(--text-muted);
      padding: 8px 12px;
      font-size: 11px;
      font-weight: 700;
      border-radius: 8px;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 7px;
      width: 100%;
      transition: all 0.18s ease;
    }
    .btn-toggle:hover {
      background: rgba(255, 255, 255, 0.12);
      color: #ffffff;
    }
    .btn-toggle.active {
      background: #10b981;
      border-color: #10b981;
      color: #ffffff;
      box-shadow: 0 2px 10px rgba(16, 185, 129, 0.35);
    }
    .btn-toggle.active.btn-cams {
      background: #3b82f6;
      border-color: #3b82f6;
      box-shadow: 0 2px 10px rgba(59, 130, 246, 0.35);
    }
  </style>

  <!-- Load Three.js, OrbitControls, and TransformControls from standard CDN -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/TransformControls.js"></script>
</head>
<body>
  <div id="canvas-container"></div>
  <div id="toast-msg">Copied!</div>
  <div id="cam-tooltip">&#128247; Camera</div>

  <!-- Metashape Viewport Top Floating Toolbar -->
  <div class="metashape-toolbar">
    <button id="tb-pts" class="tb-btn ', ifelse(mode == "points", "active", ""), '" onclick="handlePointsToolbarClick()" title="Toggle Point Cloud Visibility">&#128065;&#65039; Points</button>
    <button id="tb-mode" class="tb-btn ', ifelse(mode == "splats", "active", ""), '" onclick="handleSplatsToolbarClick()" title="Toggle Gaussian Splats Mode">&#10024; Splats</button>
    <button id="tb-mesh" class="tb-btn ', ifelse(mode == "mesh", "active", ""), '" onclick="handleMeshToolbarClick()" title="Toggle 3D Polygonal Mesh Mode">&#128314; Mesh</button>
    <div class="tb-divider"></div>
    <button id="tb-cams" class="tb-btn ', ifelse(has_cameras && show_cameras, "active tb-cams", ""), '" onclick="toggleCamerasVisibility()" title="Toggle Camera Frustums Visibility">&#128247; Cameras</button>
    <button id="tb-sight" class="tb-btn ', ifelse(has_cameras && sight_lines, "active", ""), '" onclick="toggleSightLines()" title="Toggle Sight Rays Visibility">&#128208; Sight Rays</button>
    <button id="tb-gizmo" class="tb-btn ', ifelse(trackball, "active", ""), '" onclick="toggleGizmo()" title="Toggle Navigation Sphere Gizmo">&#127760; Trackball</button>
    <button id="tb-grid" class="tb-btn active" onclick="toggleGrid()" title="Toggle Ground Grid">&#128306; Grid</button>
    <div class="tb-divider"></div>
    <button id="tb-ruler" class="tb-btn" onclick="toggleRulerTool()" title="3D Ruler: Measure euclidean distance between points">&#128207; 3D Ruler</button>
    <button id="tb-proj" class="tb-btn" onclick="toggleCameraProjection()" title="Toggle Perspective / Orthographic CAD Projection">&#128208; Ortho</button>
    <div class="tb-divider"></div>
    <button id="tb-theme" class="tb-btn" onclick="toggleBackgroundTheme()" title="Toggle Light / Dark Theme">&#9728;&#65039; Theme</button>
    <button id="tb-export" class="tb-btn" onclick="openExportModal()" title="Export Filtered Cloud (.ply/.txt) or Camera Poses" style="color: #34d399; font-weight: 700;">&#128190; Export</button>
    <button class="tb-btn" onclick="captureSnapshot()" title="Capture High-Resolution Snapshot">&#128248; Snapshot</button>
  </div>

  <!-- Hidden File Pickers -->
  <input type="file" id="input-model-file" multiple accept=".ply,.splat,.txt,.obj" style="display:none;" onchange="handleModelFileSelect(event)">
  <input type="file" id="input-mesh-file" accept=".ply,.obj" style="display:none;" onchange="handleMeshFileSelect(event)">
  <input type="file" id="input-texture-file" accept="image/*" style="display:none;" onchange="handleTextureFileSelect(event)">
  <input type="file" id="input-cams-file" multiple accept=".txt,.xml" style="display:none;" onchange="handleCamsFileSelect(event)">
  <input type="file" id="input-photos-file" multiple accept="image/*" style="display:none;" onchange="handlePhotosSelect(event)">
  <input type="file" id="input-photos-dir" webkitdirectory directory multiple style="display:none;" onchange="handlePhotosSelect(event)">
  <input type="file" id="input-batch-files" multiple accept=".ply,.splat,.txt,.obj,.xml,image/*" style="display:none;" onchange="handleBatchFilesSelect(event)">

  <!-- Welcome / Empty State Dropzone Modal (shown if no splat loaded initially) -->
  <div id="dropzone-welcome" class="dropzone-overlay" style="display: ', ifelse(has_initial_model, "none", "flex"), ';">
    <div class="dropzone-card" id="dropzone-card-welcome" style="position: relative; max-width: 620px;">
      <button class="modal-close-btn" style="position: absolute; top: 16px; right: 16px;" onclick="closeWelcomeModal()" title="Close">&#10005;</button>
      <div class="dropzone-icon">&#127793;</div>
      <div class="dropzone-title">3D Data & Layer Manager</div>
      <div class="dropzone-desc">
        Import point cloud, calibrated cameras, and photos together or individually.<br>The 3D scene will dynamically update as assets are loaded.
      </div>

      <!-- 3 Specific Asset Fields -->
      <div style="display: flex; flex-direction: column; gap: 8px; width: 100%; margin-top: 4px;">
        <!-- Field 1: Pontos 3D -->
        <div class="field-card ', ifelse(has_initial_model, "loaded", ""), '" id="welcome-model-card">
          <div class="field-header">
            <span class="field-icon">&#128230;</span>
            <div>
              <div class="field-title">Point Cloud / Splat (.ply, .splat, .txt)</div>
              <div id="welcome-model-status" class="field-status ', ifelse(has_initial_model, "loaded", "empty"), '">
                ', ifelse(has_initial_model, paste0("&#10003; ", format(n, big.mark = ","), " points loaded"), "No file selected"), '
              </div>
            </div>
          </div>
          <button class="field-btn" onclick="triggerModelImport()">&#128194; Select Points</button>
        </div>

        <!-- Field 2: C&#226;meras SFM -->
        <div class="field-card ', ifelse(has_cameras, "loaded", ""), '" id="welcome-cams-card">
          <div class="field-header">
            <span class="field-icon">&#128247;</span>
            <div>
              <div class="field-title">Camera Poses (images.txt, cameras.xml)</div>
              <div id="welcome-cams-status" class="field-status ', ifelse(has_cameras, "loaded", "empty"), '">
                ', ifelse(has_cameras, paste0("&#10003; ", length(parsed_cameras), " calibrated cameras"), "No cameras selected"), '
              </div>
            </div>
          </div>
          <button class="field-btn" onclick="triggerCamsImport()">&#128247; Select Cameras</button>
        </div>

        <!-- Field 3: Fotos -->
        <div class="field-card ', ifelse(has_cameras && has_img_files, "loaded", ""), '" id="welcome-photos-card">
          <div class="field-header">
            <span class="field-icon">&#128444;&#65039;</span>
            <div>
              <div class="field-title">Camera Photos (Directory or Images)</div>
              <div id="welcome-photos-status" class="field-status ', ifelse(has_cameras && has_img_files, "loaded", "empty"), '">
                ', ifelse(has_cameras && has_img_files, "&#10003; Embedded photos", "No photos linked"), '
              </div>
            </div>
          </div>
          <div style="display: flex; gap: 6px;">
            <button class="field-btn" onclick="triggerPhotosDirImport()" title="Link entire photo folder">&#128193; Folder</button>
            <button class="field-btn" onclick="triggerPhotosImport()" title="Select individual photo files">&#128444;&#65039; Photos</button>
          </div>
        </div>
      </div>

      <!-- Super Batch Button -->
      <button class="btn-batch-primary" onclick="triggerBatchImport()">
        &#9889; Select Multiple Files Together (Points + Cameras + Photos)
      </button>

      <div style="display: flex; align-items: center; justify-content: space-between; width: 100%; margin-top: 2px;">
        <span style="font-size: 11px; color: var(--text-muted);">
          &#9889; 100% Client-Side WebGL GPU &bull; Zero Server
        </span>
        <button class="btn-enter-scene" onclick="closeWelcomeModal()">
          &#9654; Enter 3D Scene
        </button>
      </div>
    </div>
  </div>

  <!-- Dragging Overlay (appears when user drags any file over the viewport) -->
  <div id="dropzone-drag" class="dropzone-overlay" style="display: none;">
    <div class="dropzone-card drag-active">
      <div class="dropzone-icon">&#128229;</div>
      <div class="dropzone-title">Drop your 3D Files or Cameras Here</div>
      <div class="dropzone-desc">Instant automatic asset classification and loading!</div>
    </div>
  </div>

  <!-- Photo Modal Lightbox -->
  <div id="photo-modal" class="modal-backdrop" onclick="closeModalOnBackdrop(event)">
    <div class="modal-card">
      <div class="modal-header">
        <div class="modal-title-box">
          <span style="font-size: 18px;">&#128248;</span>
          <div>
            <div id="modal-cam-name" class="modal-cam-name">Photo</div>
            <div id="modal-cam-meta" class="modal-cam-meta">Metadata</div>
          </div>
        </div>
        <button class="modal-close-btn" onclick="closePhotoModal()">&#10005;</button>
      </div>
      <div class="modal-body" style="position: relative; display: flex; align-items: center; justify-content: center; min-height: 240px; background: #020617;">
        <div id="modal-img-loading" style="display: none; position: absolute; color: #38bdf8; font-size: 13px; font-weight: 600; gap: 8px; align-items: center;">
          <span>&#9203; Loading Photo...</span>
        </div>
        <img id="modal-img" src="" alt="Camera Photo" class="modal-photo">
        <div id="modal-photo-fallback" style="display: none; flex-direction: column; align-items: center; justify-content: center; padding: 24px; gap: 10px; color: var(--text-muted); text-align: center;">
          <span style="font-size: 32px;">&#128247;</span>
          <div id="modal-fallback-msg" style="font-size: 12px; max-width: 320px; line-height: 1.5;">Photo thumbnail not loaded directly. Link the photos folder once for instant preview:</div>
          <button class="btn-action" style="background: #10b981; color: #fff; font-weight: 700; padding: 8px 16px; border-radius: 8px; cursor: pointer; margin-top: 4px;" onclick="triggerPhotosDirImport()">
            &#128193; Link Photos Folder (1 Click)
          </button>
        </div>
      </div>
      <div class="modal-footer">
        <button class="btn-modal-nav" onclick="prevPhoto()">&#9664; Previous</button>
        <button class="btn-modal-action" onclick="lookThroughCurrentCamera()">&#127909; Look Through Camera</button>
        <button class="btn-modal-nav" onclick="nextPhoto()">Next &#9654;</button>
      </div>
    </div>
  </div>

  <!-- RealityCapture Style Export Modal -->
  <div id="export-modal" class="modal-backdrop" onclick="closeExportModalOnBackdrop(event)">
    <div class="modal-card" style="max-width: 540px; width: 92%;">
      <div class="modal-header">
        <div class="modal-title-box">
          <span style="font-size: 20px;">&#128190;</span>
          <div>
            <div class="modal-cam-name">Export 3D Data & Cameras</div>
            <div class="modal-cam-meta">Direct client-side downloads of filtered assets</div>
          </div>
        </div>
        <button class="modal-close-btn" onclick="closeExportModal()" title="Close">&#10005;</button>
      </div>
      <div class="modal-body" style="background: #0f172a; padding: 18px; display: flex; flex-direction: column; gap: 12px;">
        <div style="font-size: 11px; color: var(--text-muted); line-height: 1.5;">
          Exported files will precisely contain only the active visible points (with online noise filter, bounding region, and clipping box applied).
        </div>

        <label style="display: flex; align-items: center; gap: 8px; font-size: 11px; cursor: pointer; color: #f8fafc; background: rgba(255,255,255,0.04); padding: 8px 12px; border-radius: 8px; border: 1px solid rgba(255,255,255,0.08);">
          <input type="checkbox" id="chk-export-world" checked style="accent-color: #10b981; width: 16px; height: 16px;">
          <span><strong>Preserve Real-World Coordinates</strong> (Maintains original COLMAP/SFM world scale & pose)</span>
        </label>

        <div style="display: flex; flex-direction: column; gap: 8px;">
          <!-- Export 1: Filtered PLY -->
          <div class="field-card" style="cursor: default;">
            <div class="field-header">
              <span class="field-icon">&#128230;</span>
              <div>
                <div class="field-title">Filtered Point Cloud (.PLY)</div>
                <div class="field-status">Binary Little-Endian with XYZ coordinates and RGB colors (MeshLab / CloudCompare / Gaussian Splatting)</div>
              </div>
            </div>
            <button class="field-btn" style="background: #10b981; color: #fff; font-weight: 700;" onclick="exportFilteredPLY()">Download .PLY</button>
          </div>

          <!-- Export 2: Filtered TXT -->
          <div class="field-card" style="cursor: default;">
            <div class="field-header">
              <span class="field-icon">&#128196;</span>
              <div>
                <div class="field-title">Point Cloud Table (.TXT)</div>
                <div class="field-status">Formatted table with X Y Z R G B columns for R or Python analysis</div>
              </div>
            </div>
            <button class="field-btn" onclick="exportFilteredTXT()">Download .TXT</button>
          </div>

          <!-- Export 3: Cameras TXT -->
          <div class="field-card" style="cursor: default;">
            <div class="field-header">
              <span class="field-icon">&#128247;</span>
              <div>
                <div class="field-title">Camera Poses (.TXT / COLMAP)</div>
                <div class="field-status">Standard images.txt format with normalized quaternions and camera translations</div>
              </div>
            </div>
            <button class="field-btn" onclick="exportCamerasTXT()">Download Cameras</button>
          </div>
        </div>
      </div>
      <div class="modal-footer" style="justify-content: flex-end;">
        <button class="btn-modal-nav" onclick="closeExportModal()">Close</button>
      </div>
    </div>
  </div>

  <!-- HUD Control Panel -->
  <div class="hud-panel">
    <div class="hud-header">
      <div class="hud-title">
        <span>&#127793; 3D Viewer</span>
        <span class="hud-badge">pliman</span>
      </div>
      <div style="display: flex; align-items: center; gap: 8px;">
        <button class="btn-import-header" onclick="openWelcomeModal()" title="Manage or import layers">&#128194; Import</button>
        <div id="fps-counter" class="fps-badge">60 FPS</div>
      </div>
    </div>

    <!-- Scene Layers & Asset Management -->
    <div class="control-group">
      <div class="control-label">
        <span>Scene Layers</span>
        <button class="btn-reset-sm" onclick="openWelcomeModal()" title="Open full import panel">&#128194; Import...</button>
      </div>
      <div class="layer-list">
        <!-- Layer 1: Pontos / Splat -->
        <div class="layer-item">
          <div class="layer-info">
            <span class="layer-icon">&#128230;</span>
            <div class="layer-texts">
              <div class="layer-name">Point Cloud</div>
              <div id="hud-model-status" class="layer-meta">', ifelse(has_initial_model, paste0(format(n, big.mark = ","), " pts"), "None"), '</div>
            </div>
          </div>
          <div class="layer-actions">
            <button id="btn-pts-toggle" class="btn-layer-toggle ', ifelse(has_initial_model && mode != "mesh", "active", ""), '" onclick="togglePointsVisibility()" title="Toggle Points Visibility">&#128065;&#65039;</button>
            <button class="btn-layer-sub" onclick="triggerModelImport()" title="Load / Replace Points">&#128194;</button>
          </div>
        </div>

        <!-- Layer 1b: 3D Surface Mesh -->
        <div class="layer-item">
          <div class="layer-info">
            <span class="layer-icon">&#128314;</span>
            <div class="layer-texts">
              <div class="layer-name">3D Surface Mesh</div>
              <div id="hud-mesh-status" class="layer-meta">', ifelse(has_faces, paste0(format(num_faces, big.mark = ","), " faces"), "None"), '</div>
            </div>
          </div>
          <div class="layer-actions">
            <button id="btn-mesh-toggle" class="btn-layer-toggle ', ifelse(has_faces && mode == "mesh", "active", ""), '" onclick="toggleMeshVisibility()" title="Toggle Mesh Visibility">&#128065;&#65039;</button>
            <button class="btn-layer-sub" onclick="triggerMeshImport()" title="Load / Replace Mesh">&#128194;</button>
          </div>
        </div>

        <!-- Layer 2: C&#226;meras SFM -->
        <div class="layer-item">
          <div class="layer-info">
            <span class="layer-icon">&#128247;</span>
            <div class="layer-texts">
              <div class="layer-name">SFM Cameras</div>
              <div id="hud-cams-status" class="layer-meta">', ifelse(has_cameras, paste0(length(parsed_cameras), " cams"), "None"), '</div>
            </div>
          </div>
          <div class="layer-actions">
            <button id="btn-cams-toggle" class="btn-layer-toggle ', ifelse(has_cameras && show_cameras, "active", ""), '" onclick="toggleCamerasVisibility()" title="Toggle Cameras Visibility">&#128065;&#65039;</button>
            <button class="btn-layer-sub" onclick="triggerCamsImport()" title="Load / Replace Cameras">&#128247;</button>
          </div>
        </div>

        <!-- Layer 3: Fotos -->
        <div class="layer-item">
          <div class="layer-info">
            <span class="layer-icon">&#128444;&#65039;</span>
            <div class="layer-texts">
              <div class="layer-name">Camera Photos</div>
              <div id="hud-photos-status" class="layer-meta">', ifelse(has_cameras && has_img_files, "Active photos", "None"), '</div>
            </div>
          </div>
          <div class="layer-actions">
            <button class="btn-layer-sub" onclick="triggerPhotosDirImport()" title="Link Photos Folder">&#128193;</button>
            <button class="btn-layer-sub" onclick="triggerPhotosImport()" title="Link Photo Files">&#128444;&#65039;</button>
          </div>
        </div>
      </div>
    </div>

    <!-- Color Mode -->
    <div class="control-group">
      <div class="control-label">Color & Shading Mode</div>
      <div class="segmented-control-5">
        <button id="btn-col-rgb" class="seg-btn ', ifelse(color_by == "rgb", "active", ""), '" onclick="setColorMode(\'rgb\')">RGB</button>
        <button id="btn-col-clay" class="seg-btn ', ifelse(color_by == "clay", "active", ""), '" onclick="setColorMode(\'clay\')" title="Studio Clay: gray geometry inspection">Clay</button>
        <button id="btn-col-h" class="seg-btn ', ifelse(color_by == "height", "active", ""), '" onclick="setColorMode(\'height\')" title="Elevation Heatmap (Z Gradient)">Elevation</button>
        <button id="btn-col-dens" class="seg-btn ', ifelse(color_by == "density", "active", ""), '" onclick="setColorMode(\'density\')" title="Point Density Heatmap (SOR)">Density</button>
        <button id="btn-col-sem" class="seg-btn ', ifelse(color_by == "semantic", "active", ""), '" onclick="setColorMode(\'semantic\')">SAM</button>
      </div>
    </div>

    <!-- 3D Polygonal Mesh Controls -->
    <div class="control-group" id="group-mesh-controls" style="display: ', ifelse(mode == "mesh", "flex", "none"), '; flex-direction: column; gap: 8px;">
      <div class="control-label">
        <span>&#128314; Mesh Shading & Material</span>
        <span id="lbl-mesh-mode-name" class="val">', ifelse(color_by == "clay", "Clay", "Solid"), '</span>
      </div>
      <div class="segmented-control-4">
        <button id="btn-mesh-solid" class="seg-btn ', ifelse(color_by != "clay", "active", ""), '" onclick="setMeshShadingMode(\'solid\')" title="Solid smooth shaded surface with lighting">Solid</button>
        <button id="btn-mesh-clay" class="seg-btn ', ifelse(color_by == "clay", "active", ""), '" onclick="setMeshShadingMode(\'clay\')" title="Studio Clay matte plaster surface for inspecting geometric details">Clay</button>
        <button id="btn-mesh-wire" class="seg-btn" onclick="setMeshShadingMode(\'wireframe\')" title="Wireframe triangular polygon mesh lattice">Wire</button>
        <button id="btn-mesh-solidwire" class="seg-btn" onclick="setMeshShadingMode(\'solid_wireframe\')" title="Solid surface with wireframe overlay (CAD)">Solid+Wire</button>
      </div>

      <!-- Texture & Shading Toggles -->
      <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 4px;">
        <button id="btn-mesh-texture-toggle" class="btn-action ', ifelse(tex_data_uri != "", "active", ""), '" onclick="toggleMeshTexture()" title="Toggle texture mapping on or off">&#128444;&#65039; Texture: ', ifelse(tex_data_uri != "", "ON", "OFF"), '</button>
        <button class="btn-action" onclick="triggerTextureImport()" title="Load texture image file (.png, .jpg) or texture atlas">&#128193; Load Texture</button>
      </div>

      <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 4px;">
        <button id="btn-mesh-flat-toggle" class="btn-action" onclick="toggleMeshFlatShading()" title="Toggle smooth vertex normals vs faceted flat polygons">&#10024; Smooth Shading</button>
        <button id="btn-mesh-doubleside-toggle" class="btn-action active" onclick="toggleMeshDoubleSided()" title="Toggle two-sided vs front-facing rendering">&#8596;&#65039; Double-Sided</button>
      </div>

      <!-- Wireframe Opacity Slider (when wireframe is active) -->
      <div id="mesh-wire-controls" style="display: none; flex-direction: column; gap: 4px; background: rgba(0,0,0,0.25); padding: 6px; border-radius: 6px;">
        <div style="display: flex; justify-content: space-between; font-size: 10px; color: var(--text-muted);">
          <span>Wireframe Opacity</span>
          <span id="lbl-wire-opac" class="val">45%</span>
        </div>
        <input type="range" id="slider-wire-opac" min="10" max="100" step="5" value="45" oninput="updateWireOpacity(this.value)">
      </div>
    </div>

    <!-- Point Size / Splat Fill Controls -->
    <div id="group-point-controls" style="display: ', ifelse(mode != "mesh", "flex", "none"), '; flex-direction: column; gap: 8px;">
      <div class="control-group">
        <div class="control-label">
          <span id="lbl-size-title">', ifelse(mode == "points", "Point Size", "Splat Closure / Fill"), '</span>
          <span id="lbl-size" class="val">', point_size, ifelse(mode == "points", " px", "x (Surface Fill)"), '</span>
        </div>
        <input type="range" id="slider-size" min="0.5" max="35.0" step="0.5" value="', point_size, '" oninput="updatePointSize(this.value)">
      </div>

      <!-- Point Shape Selector -->
      <div class="control-group" id="group-shape" style="display: ', ifelse(mode == "points", "flex", "none"), ';">
        <div class="control-label">Point Shape</div>
        <button id="btn-toggle-shape" class="btn-toggle" onclick="togglePointShape()">&#9898; Smooth Circle (Antialiased)</button>
      </div>
    </div>

    <!-- RealityCapture Online Noise Filter (SOR / Density) -->
    <div class="control-group">
      <div class="control-label">
        <span>&#127919; Online Noise Filter (Density / SOR)</span>
        <span id="lbl-noise-val" class="val">Off</span>
      </div>
      <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 4px; margin-bottom: 6px;">
        <button class="btn-action active" id="btn-noise-0" onclick="setNoiseFilterPreset(0)">Off</button>
        <button class="btn-action" id="btn-noise-1" onclick="setNoiseFilterPreset(15)">Light</button>
        <button class="btn-action" id="btn-noise-2" onclick="setNoiseFilterPreset(35)">Medium</button>
        <button class="btn-action" id="btn-noise-3" onclick="setNoiseFilterPreset(65)">Strong</button>
      </div>
      <input type="range" id="slider-noise" min="0" max="100" step="1" value="0" oninput="updateNoiseFilter(this.value)">
      <div style="display: flex; justify-content: space-between; font-size: 10px; color: var(--text-muted); margin-top: 2px;">
        <span>0% (Unfiltered)</span>
        <span id="lbl-noise-stats">0 noise points filtered</span>
        <span>100% (Ultra dense)</span>
      </div>
    </div>

    <!-- RealityCapture Reconstruction Region (3D Section Box) -->
    <div class="control-group">
      <div class="control-label">
        <span>&#9986;&#65039; Reconstruction Region (3D Box)</span>
        <button id="btn-toggle-clipbox" class="btn-reset-sm" onclick="toggleClipBox()">Enable Box</button>
      </div>
      <div id="clipbox-controls" style="display: none; flex-direction: column; gap: 7px; margin-top: 6px;">
        <!-- Gizmo Interactive Mode Toggle -->
        <div style="display: flex; gap: 3px; background: rgba(0,0,0,0.25); padding: 2px; border-radius: 6px; border: 1px solid rgba(255,255,255,0.08);">
          <button id="btn-gizmo-trans" class="btn-tool active" style="flex: 1; padding: 4px 2px; font-size: 10px;" onclick="setGizmoMode(\'translate\')" title="Translate / Move Box (W)">&#9995; Move</button>
          <button id="btn-gizmo-rot" class="btn-tool" style="flex: 1; padding: 4px 2px; font-size: 10px;" onclick="setGizmoMode(\'rotate\')" title="Rotate Box Angles (E)">&#128260; Rotate</button>
          <button id="btn-gizmo-scale" class="btn-tool" style="flex: 1; padding: 4px 2px; font-size: 10px;" onclick="setGizmoMode(\'scale\')" title="Resize Box Dimensions (R)">&#128208; Scale</button>
        </div>

        <!-- Quick Align Actions -->
        <div style="display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 4px;">
          <button class="btn-action" onclick="autoFitClipBox()" title="Auto-fit box around subject core">&#127793; Auto Subject</button>
          <button class="btn-action" onclick="cutClipBoxBase()" title="Lift bottom face to crop floor/table">&#9986;&#65039; Cut Floor</button>
          <button class="btn-action" onclick="levelClipBox()" title="Level box to horizon plane">&#9878;&#65039; Level</button>
        </div>

        <!-- 3D Rotation Sliders (Pitch, Roll, Yaw in Degrees) -->
        <div style="background: rgba(15,23,42,0.6); padding: 6px; border-radius: 6px; border: 1px solid rgba(255,255,255,0.06); display: flex; flex-direction: column; gap: 4px;">
          <div style="display: flex; justify-content: space-between; font-size: 10px; color: var(--text-muted); font-weight: 600;">
            <span>3D ROTATION (DEGREES)</span>
            <span id="lbl-box-rot" style="color: #38bdf8;">0&deg; / 0&deg; / 0&deg;</span>
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-x" style="min-width: 46px;">Pitch X</span>
            <input type="range" id="clip-rot-x" min="-180" max="180" value="0" step="1" oninput="updateBoxFromUI()">
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-y" style="min-width: 46px;">Roll Y</span>
            <input type="range" id="clip-rot-y" min="-180" max="180" value="0" step="1" oninput="updateBoxFromUI()">
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-z" style="min-width: 46px;">Yaw Z</span>
            <input type="range" id="clip-rot-z" min="-180" max="180" value="0" step="1" oninput="updateBoxFromUI()">
          </div>
        </div>

        <!-- Box Dimensions Sliders (Size W, L, H) -->
        <div style="background: rgba(15,23,42,0.6); padding: 6px; border-radius: 6px; border: 1px solid rgba(255,255,255,0.06); display: flex; flex-direction: column; gap: 4px;">
          <div style="display: flex; justify-content: space-between; font-size: 10px; color: var(--text-muted); font-weight: 600;">
            <span>BOX DIMENSIONS</span>
            <span id="lbl-box-dims" style="color: #4ade80;">1.0 x 1.0 x 1.0m</span>
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-x" style="min-width: 46px;">Width X</span>
            <input type="range" id="clip-size-x" min="10" max="300" value="100" step="1" oninput="updateBoxFromUI()">
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-y" style="min-width: 46px;">Length Y</span>
            <input type="range" id="clip-size-y" min="10" max="300" value="100" step="1" oninput="updateBoxFromUI()">
          </div>
          <div class="axis-control-row">
            <span class="axis-badge axis-z" style="min-width: 46px;">Height Z</span>
            <input type="range" id="clip-size-z" min="10" max="300" value="100" step="1" oninput="updateBoxFromUI()">
          </div>
        </div>

        <!-- Crop and Export Section -->
        <div style="display: flex; flex-direction: column; gap: 4px; margin-top: 2px;">
          <button id="btn-crop-cloud" class="btn-action" style="background: linear-gradient(135deg, #059669, #10b981); color: #fff; font-weight: 700; border-color: #34d399; justify-content: center; padding: 6px;" onclick="applyCropPointCloud()" title="Permanently crop and keep only points inside the 3D Box">
            &#9986;&#65039; Crop Point Cloud (Keep Inside)
          </button>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 4px;">
            <button id="btn-undo-crop" class="btn-action" style="display: none; justify-content: center;" onclick="undoCropPointCloud()" title="Revert to previous uncropped point cloud">
              &#8634; Undo Crop
            </button>
            <button id="btn-export-cropped-ply" class="btn-action" style="justify-content: center;" onclick="exportCroppedPLY()" title="Download cropped points as clean 3D PLY file">
              &#128190; Export PLY
            </button>
          </div>
        </div>
      </div>
    </div>

    <!-- 3D Ruler Measurements Panel -->
    <div class="control-group" id="group-measurements" style="display: none;">
      <div class="control-label">
        <span>&#128207; 3D Measurements (Ruler)</span>
        <button class="btn-reset-sm" onclick="clearMeasurements()">Clear All</button>
      </div>
      <div id="measurements-list" class="metrics-box" style="gap: 5px;">
        <div style="font-size: 11px; color: #facc15;" id="ruler-instruction">
          1. Click 1st point on point cloud...
        </div>
      </div>
    </div>

    <!-- Metashape Bounding Region & Noise Crop -->
    <div class="control-group">
      <div class="control-label">
        <span>&#127919; Bounding Region (Noise Crop)</span>
        <span id="lbl-crop-val" class="val">100%</span>
      </div>
      <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 4px; margin-bottom: 6px;">
        <button class="btn-action" onclick="setCropPreset(0.35)" title="Focus on central subject and crop surrounding environment walls">&#127793; Focus on Subject</button>
        <button class="btn-action" onclick="setCropPreset(1.0)" title="Display entire scene including surrounding environment">&#127760; Full Scene</button>
      </div>
      <input type="range" id="slider-crop" min="0.05" max="1.0" step="0.02" value="1.0" oninput="updateCropRadius(this.value)">
      <div style="font-size: 10px; color: var(--text-muted); margin-top: 2px;">
        Isolate the central subject and eliminate distant background noise!
      </div>
    </div>

    <!-- 3D Model Alignment / Rotation Across 3 Axes -->
    <div class="control-group">
      <div class="control-label">
        <span>Model Alignment (3 Axes)</span>
        <div style="display: flex; gap: 4px;">
          <button class="btn-reset-sm" onclick="copyRotationCode()" title="Copy rotation parameters to clipboard">Copy</button>
          <button class="btn-reset-sm" onclick="resetModelRotation()">Reset</button>
        </div>
      </div>
      <div class="segmented-control" style="margin-bottom: 8px;">
        <button id="btn-space-world" class="seg-btn active" onclick="setAlignmentSpace(\'world\')" title="Gizmo aligned with Ground Grid (World coordinates for leveling)">&#127760; World (Grid)</button>
        <button id="btn-space-local" class="seg-btn" onclick="setAlignmentSpace(\'local\')" title="Gizmo aligned with Object (Local coordinates for intrinsic tilt)">&#129482; Local (Object)</button>
      </div>
      <div class="axis-control-row">
        <span class="axis-badge axis-x">X (Pitch)</span>
        <input type="range" id="rot-x-slider" min="-180" max="180" value="', init_rot_x, '" step="1" oninput="updateModelRotation()">
        <span id="rot-x-val" class="axis-val">', init_rot_x, '&deg;</span>
      </div>
      <div class="axis-control-row">
        <span class="axis-badge axis-y">Y (Roll)</span>
        <input type="range" id="rot-y-slider" min="-180" max="180" value="', init_rot_y, '" step="1" oninput="updateModelRotation()">
        <span id="rot-y-val" class="axis-val">', init_rot_y, '&deg;</span>
      </div>
      <div class="axis-control-row">
        <span class="axis-badge axis-z">Z (Yaw)</span>
        <input type="range" id="rot-z-slider" min="-180" max="180" value="', init_rot_z, '" step="1" oninput="updateModelRotation()">
        <span id="rot-z-val" class="axis-val">', init_rot_z, '&deg;</span>
      </div>
    </div>

    <!-- Display Fine-Tuning & Scales (Sliders only - Toggles are on Top Bar) -->
    <div class="control-group">
      <div class="control-label">
        <span>Visual Adjustments & Scale</span>
      </div>
      <div class="axis-control-row">
        <span style="font-size: 10px; color: var(--text-muted); min-width: 84px;">Point Opacity</span>
        <input type="range" id="slider-opac" min="0.1" max="1.0" step="0.05" value="1.0" oninput="updateOpacityMultiplier(this.value)">
        <span id="lbl-opac" class="axis-val">1.00x</span>
      </div>
      <div class="axis-control-row" id="row-cam-scale" style="margin-top: 4px; display: ', ifelse(has_cameras, "flex", "none"), ';">
        <span style="font-size: 10px; color: var(--text-muted); min-width: 84px;">Frustum Scale</span>
        <input type="range" id="slider-cam-scale" min="0.2" max="3.0" step="0.1" value="', camera_size, '" oninput="updateCameraScale(this.value)">
        <span id="lbl-cam-scale" class="axis-val">', sprintf("%.1fx", camera_size), '</span>
      </div>
      <div class="axis-control-row" style="margin-top: 4px;">
        <span style="font-size: 10px; color: var(--text-muted); min-width: 84px;">Sphere Opacity</span>
        <input type="range" id="slider-sphere-opac" min="0.0" max="0.20" step="0.005" value="', init_sphere_opac, '" oninput="updateSphereOpacity(this.value)">
        <span id="lbl-sphere-opac" class="axis-val">', round(init_sphere_opac * 100), '%</span>
      </div>
    </div>

    <!-- Camera Presets -->
    <div class="control-group">
      <div class="control-label">Camera Angle Presets</div>
      <div class="cam-grid">
        <button class="btn-action" onclick="setCameraPreset(\'front\')">Front</button>
        <button class="btn-action" onclick="setCameraPreset(\'side\')">Side</button>
        <button class="btn-action" onclick="setCameraPreset(\'top\')">Canopy</button>
        <button class="btn-action" onclick="setCameraPreset(\'reset\')">Reset</button>
      </div>
    </div>

    <!-- Phenotypic Biometrics HUD -->
    <div class="control-group">
      <div class="control-label">Biometrics & Spatial Extent</div>
      <div class="metrics-box">
        <div class="metric-row">
          <span>Total Points:</span>
          <span id="lbl-total-pts" class="metric-val">', format(n, big.mark = ","), '</span>
        </div>
        <div class="metric-row">
          <span>Active Visible:</span>
          <span id="lbl-visible-splats" class="metric-val">', format(n, big.mark = ","), '</span>
        </div>
        <div class="metric-row">
          <span>Subject Clustered Scale:</span>
          <span id="lbl-scale-extent" class="metric-val">', sprintf("%.2f m", target_dim), '</span>
        </div>
        <div class="metric-row">
          <span>Global Bounding Box:</span>
          <span id="lbl-bbox-extent" class="metric-val">', sprintf("%.1f x %.1f x %.1f m", rob_dx, rob_dy, rob_dz), '</span>
        </div>
        <div class="metric-row">
          <span>Subject Origin (X, Y):</span>
          <span id="lbl-origin-xy" class="metric-val">', sprintf("%.2f, %.2f", center_x, center_y), '</span>
        </div>
        <div class="metric-row">
          <span>Elevation (Z):</span>
          <span id="lbl-origin-z" class="metric-val">', sprintf("%.2f m", center_z), '</span>
        </div>
      </div>
    </div>

    <!-- Snapshot / Export -->
    <button class="btn-primary" onclick="captureSnapshot()">
      &#128247; Capture High-Res Snapshot
    </button>
  </div>

  <!-- Bottom Bar: Real-time 3D Coordinate Badge + Controls Hint -->
  <div class="bottom-bar">
    <div class="coords-hud" id="coords-panel">
      <div class="coords-title">
        <span>&#128205; 3D Coordinates</span>
        <span id="coords-status" style="font-size: 9px; color: var(--text-muted); font-weight: normal;">Target Center</span>
      </div>
      <div class="coords-vals" id="coords-text">
        X: ', sprintf("%.3f", center_x), ' | Y: ', sprintf("%.3f", center_y), ' | Z: ', sprintf("%.3f", center_z), '
      </div>
    </div>
    <div class="bottom-hint">
      Left Drag: Orbit | Right Drag: Pan | Scroll: Zoom | Double Click: Set Pivot Center | Click Camera: View Photo
    </div>
  </div>

  <script>
    // Embedded Centered Geometry Data (Robust Subject Centering)
    let numGaussians = ', n, ';
    let numFaces = ', num_faces, ';
    let rawPositions = new Float32Array([', pos_str, ']);
    let rawColorsRGB = new Float32Array([', rgb_str, ']);
    let rawOpacities = new Float32Array([', opac_str, ']);
    let rawScales = new Float32Array([', scale_str, ']);
    let rawSemantics = new Int32Array([', sem_str, ']);
    let rawFaces = new Uint32Array([', faces_str, ']);
    let hasModelFaces = rawFaces.length > 0;
    let hasModelColors = ', ifelse(has_colors, "true", "false"), ';
    let initialTextureUri = "', tex_data_uri, '";
    let semanticLabels = ', sem_map_js, ';
    let rawCameras = ', cameras_json, ';
    let targetDim = ', target_dim, ';
    let centerOrigin = [', center_x, ', ', center_y, ', ', center_z, '];

    // Palettes for Semantic Classes
    const semanticColorsPalette = {
      0: [0.55, 0.60, 0.68],
      1: [0.10, 0.85, 0.35],
      2: [0.85, 0.50, 0.15],
      3: [0.95, 0.25, 0.35],
      4: [0.20, 0.60, 0.95],
      5: [0.70, 0.30, 0.85]
    };

    // Global application state
    let currentRenderMode = "', mode, '";
    let currentColorMode = "', color_by, '";
    let isMeshVisible = true;
    let arePointsVisible = true;
    let meshShadingMode = "', ifelse(color_by == "clay", "clay", "solid"), '";
    let isMeshFlatShading = false;
    let isMeshDoubleSided = true;
    let isMeshTextureActive = ', ifelse(tex_data_uri != "", "true", "false"), ';
    let meshTextureObject = null;
    let meshGroup = null;
    let meshSurfaceObject = null;
    let meshWireframeObject = null;
    let meshGeometry = null;
    let pointSizeValue = ', point_size, ';
    let pointShapeValue = 0;
    let globalOpacityMultiplier = 1.0;
    let minAlphaCutoff = 0.01;
    let activeSemanticClasses = new Set(Object.keys(semanticLabels).map(k => parseInt(k)));
    let isGizmoVisible = ', ifelse(trackball, "true", "false"), ';
    let sphereOpacityValue = ', init_sphere_opac, ';
    let sphereScaleMultiplier = 1.0;
    let areCamerasVisible = ', ifelse(show_cameras, "true", "false"), ';
    let areSightLinesVisible = ', ifelse(sight_lines, "true", "false"), ';
    let cameraScaleMultiplier = ', camera_size, ';
    let activeCameraModalIndex = 0;

    // Three.js instances
    let scene, camera, orthoCamera, activeCamera, renderer, controls;
    let isOrtho = false;
    let modelContainer, pointsMesh, geometry, material;
    let positionsAttr, colorsAttr, alphasAttr, sizesAttr;
    let gridHelper, trackballGizmo, sphereMesh;
    let camerasGroup, sightLinesGroup;
    let cameraColliders = [];
    let cameraMeshes = [];

    // RealityCapture State Variables
    let pointDensities = new Float32Array(0);
    let sortedDensities = null;
    let sortedDistances = null;
    let maxDensityVal = 1;
    let minDensityVal = 0;
    let densityPercentiles = { p10: 1, p25: 5, p50: 15, p75: 35, p90: 60 };
    let noiseFilterPercent = 0;

    let clipBoxEnabled = false;
    let clipBoxGroup = null;
    let clipBoxMesh = null;
    let clipBoxLines = null;
    let transformControls = null;
    let gizmoMode = "translate";
    let initialUncroppedState = null;

    let isRulerActive = false;
    let rulerPointA = null;
    let rulerMeasurements = [];
    let rulerGroup = null;
    let rulerHoverMarker = null;
    let rulerTempLine = null;

    // Interactive 3D Gizmo Dragging State (Physically Accurate 3D Ray-Plane Tracking)
    let alignmentSpace = "world"; // "world" or "local"
    let ringMeshes = {};
    let ringColliders = [];
    let isDraggingRing = false;
    let activeDragAxis = null;
    let dragAxisVector = new THREE.Vector3();
    let dragPlane = new THREE.Plane();
    let dragStartVector = new THREE.Vector3();
    let dragStartQuaternion = new THREE.Quaternion();
    let dragIsEdgeOn = false;
    let dragEdgeTangent2D = new THREE.Vector2();
    let lastDragMousePos = { x: 0, y: 0 };
    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();

    // Initialize WebGL Application
    function init() {
      const container = document.getElementById("canvas-container");
      scene = new THREE.Scene();
      scene.background = new THREE.Color("', bg_hex, '");

      camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.005, 500);
      // Tight camera framing right on the centered object!
      camera.position.set(0, -targetDim * 2.2, targetDim * 1.1);
      camera.up.set(0, 0, 1);
      activeCamera = camera;

      initOrthoCamera();

      renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true, powerPreference: "high-performance" });
      renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
      renderer.setSize(window.innerWidth, window.innerHeight);
      container.appendChild(renderer.domElement);

      controls = new THREE.OrbitControls(camera, renderer.domElement);
      controls.target.set(0, 0, 0);
      controls.enableDamping = true;
      controls.dampingFactor = 0.08;
      controls.update();

      // Lighting for 3D polygonal mesh surface & clay rendering
      const ambientLight = new THREE.AmbientLight(0xffffff, 0.65);
      scene.add(ambientLight);

      const keyLight = new THREE.DirectionalLight(0xffffff, 0.75);
      keyLight.position.set(5 * targetDim, 10 * targetDim, 8 * targetDim);
      scene.add(keyLight);

      const fillLight = new THREE.DirectionalLight(0x90b0e0, 0.45);
      fillLight.position.set(-6 * targetDim, -5 * targetDim, -6 * targetDim);
      scene.add(fillLight);

      const headlight = new THREE.DirectionalLight(0xffffff, 0.40);
      headlight.position.set(0, 0, 1);
      camera.add(headlight);
      scene.add(camera);

      // Model container allows interactive 3-axis rotation / leveling
      modelContainer = new THREE.Group();
      scene.add(modelContainer);

      // Measurement group attached to model container so measurements move with model
      rulerGroup = new THREE.Group();
      modelContainer.add(rulerGroup);

      // Ground grid sized to target subject
      gridHelper = new THREE.GridHelper(targetDim * 4.0, 24, 0x10b981, 0x334155);
      gridHelper.rotation.x = Math.PI / 2;
      gridHelper.position.set(0, 0, -targetDim * 0.45);
      scene.add(gridHelper);

      // Build Metashape-style Trackball Gizmo (Silky smooth transparent glass sphere with 3-Axis Rings)
      buildTrackballGizmo();

      // Build 3D Camera Frustums (COLMAP / Metashape)
      if (rawCameras && rawCameras.length > 0) {
        buildCamerasFrustums();
      }

      // Build Points & Splats Geometry if initial model is provided
      if (numGaussians > 0) {
        computeSpatialDensities();
        computeDistanceDistribution();
        buildSplatGeometry();
        if (hasModelFaces) {
          buildMeshGeometry();
        }
        if (currentRenderMode === "mesh") {
          setRenderMode("mesh");
        } else {
          setColorMode(currentColorMode);
        }
      } else if (rawCameras && rawCameras.length > 0) {
        frameCameras();
      }

      if (initialTextureUri && initialTextureUri.length > 20) {
        loadTextureFromDataUri(initialTextureUri);
      }

      updateLayersUI();

      // Initialize initial rotations if any
      updateModelRotation();

      // Setup interactive mouse handling for viewport, coordinates & gizmo dragging
      setupInteractions();

      // Setup drag and drop for 3D models and cameras
      setupDragAndDrop();

      // Keyboard navigation for modal & 3D Gizmo modes
      window.addEventListener("keydown", (e) => {
        if (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA") return;
        const modal = document.getElementById("photo-modal");
        if (modal && modal.style.display === "flex") {
          if (e.key === "Escape") closePhotoModal();
          if (e.key === "ArrowLeft") prevPhoto();
          if (e.key === "ArrowRight") nextPhoto();
          return;
        }

        // Reconstruction Region Gizmo Mode shortcuts (W = Move, E = Rotate, R = Scale)
        if (clipBoxEnabled && typeof setGizmoMode === "function") {
          if (e.key === "w" || e.key === "W") setGizmoMode("translate");
          if (e.key === "e" || e.key === "E") setGizmoMode("rotate");
          if (e.key === "r" || e.key === "R") setGizmoMode("scale");
        }
      });

      window.addEventListener("resize", onWindowResize);
      updateLayersUI();
      animate();
    }

    // Build 3D Camera Frustums & Sight Lines
    function buildCamerasFrustums() {
      if (camerasGroup) {
        modelContainer.remove(camerasGroup);
      }
      cameraColliders = [];
      cameraMeshes = [];
      camerasGroup = new THREE.Group();
      sightLinesGroup = new THREE.Group();
      camerasGroup.visible = areCamerasVisible;
      sightLinesGroup.visible = areSightLinesVisible;

      const baseFrustumLen = Math.max(0.04, targetDim * 0.12) * cameraScaleMultiplier;

      rawCameras.forEach((cam, i) => {
        const camSubGroup = new THREE.Group();
        camSubGroup.userData = { camIndex: i, camData: cam };

        const cX = cam.pos[0], cY = cam.pos[1], cZ = cam.pos[2];
        const C = new THREE.Vector3(cX, cY, cZ);

        // Rotation matrix: rows are camera Right, Up (inverted), Forward in world space
        const R = cam.R;
        const right = new THREE.Vector3(R[0], R[1], R[2]);
        const up = new THREE.Vector3(-R[3], -R[4], -R[5]); // -Y is Up
        const fwd = new THREE.Vector3(R[6], R[7], R[8]);

        // Base rectangle dimensions
        const fx = cam.fx || 1000;
        const fy = cam.fy || 1000;
        const w = baseFrustumLen * (cam.w / (2.0 * fx));
        const h = baseFrustumLen * (cam.h / (2.0 * fy));

        const baseCenter = C.clone().addScaledVector(fwd, baseFrustumLen);
        const p1 = baseCenter.clone().addScaledVector(right, -w).addScaledVector(up, h);  // Top-Left
        const p2 = baseCenter.clone().addScaledVector(right, w).addScaledVector(up, h);   // Top-Right
        const p3 = baseCenter.clone().addScaledVector(right, w).addScaledVector(up, -h);  // Bottom-Right
        const p4 = baseCenter.clone().addScaledVector(right, -w).addScaledVector(up, -h); // Bottom-Left

        // Wireframe Frustum Lines (Red #ef4444)
        const linePositions = [];
        // Apex to 4 corners
        linePositions.push(C.x, C.y, C.z, p1.x, p1.y, p1.z);
        linePositions.push(C.x, C.y, C.z, p2.x, p2.y, p2.z);
        linePositions.push(C.x, C.y, C.z, p3.x, p3.y, p3.z);
        linePositions.push(C.x, C.y, C.z, p4.x, p4.y, p4.z);
        // Base rectangle
        linePositions.push(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
        linePositions.push(p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
        linePositions.push(p3.x, p3.y, p3.z, p4.x, p4.y, p4.z);
        linePositions.push(p4.x, p4.y, p4.z, p1.x, p1.y, p1.z);
        // Top orientation indicator line
        const topMid = p1.clone().add(p2).multiplyScalar(0.5);
        const topIndicator = topMid.clone().addScaledVector(up, h * 0.35);
        linePositions.push(p1.x, p1.y, p1.z, topIndicator.x, topIndicator.y, topIndicator.z);
        linePositions.push(p2.x, p2.y, p2.z, topIndicator.x, topIndicator.y, topIndicator.z);

        const lineGeo = new THREE.BufferGeometry();
        lineGeo.setAttribute("position", new THREE.Float32BufferAttribute(linePositions, 3));
        const lineMat = new THREE.LineBasicMaterial({ color: 0xef4444, linewidth: 2 });
        const linesMesh = new THREE.LineSegments(lineGeo, lineMat);
        camSubGroup.add(linesMesh);

        // Translucent Face Mesh
        const facePositions = [];
        // 4 pyramid sides
        facePositions.push(C.x, C.y, C.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
        facePositions.push(C.x, C.y, C.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
        facePositions.push(C.x, C.y, C.z, p3.x, p3.y, p3.z, p4.x, p4.y, p4.z);
        facePositions.push(C.x, C.y, C.z, p4.x, p4.y, p4.z, p1.x, p1.y, p1.z);
        // Base plane
        facePositions.push(p1.x, p1.y, p1.z, p3.x, p3.y, p3.z, p2.x, p2.y, p2.z);
        facePositions.push(p1.x, p1.y, p1.z, p4.x, p4.y, p4.z, p3.x, p3.y, p3.z);

        const faceGeo = new THREE.BufferGeometry();
        faceGeo.setAttribute("position", new THREE.Float32BufferAttribute(facePositions, 3));
        const faceMat = new THREE.MeshBasicMaterial({
          color: 0xef4444,
          transparent: true,
          opacity: 0.16,
          side: THREE.DoubleSide,
          depthTest: false
        });
        const faceMesh = new THREE.Mesh(faceGeo, faceMat);
        camSubGroup.add(faceMesh);

        // Thick Invisible Raycast Collider for effortless click/hover detection
        const colMesh = new THREE.Mesh(faceGeo, new THREE.MeshBasicMaterial({ visible: false, depthTest: false }));
        colMesh.userData = { camIndex: i, camData: cam, lineMesh: linesMesh, faceMesh: faceMesh };
        camSubGroup.add(colMesh);
        cameraColliders.push(colMesh);

        // Sight Line connecting camera apex to subject center (Magenta #d946ef)
        const sightPos = [C.x, C.y, C.z, 0, 0, 0];
        const sightGeo = new THREE.BufferGeometry();
        sightGeo.setAttribute("position", new THREE.Float32BufferAttribute(sightPos, 3));
        const sightMat = new THREE.LineBasicMaterial({ color: 0xd946ef, transparent: true, opacity: 0.55 });
        const sightLine = new THREE.Line(sightGeo, sightMat);
        sightLinesGroup.add(sightLine);

        cameraMeshes.push({ subGroup: camSubGroup, lines: linesMesh, faces: faceMesh, sight: sightLine });
        camerasGroup.add(camSubGroup);
      });

      camerasGroup.add(sightLinesGroup);
      modelContainer.add(camerasGroup);
    }

    // Build Metashape Trackball Sphere with 3-Axis Gimbal Rings & Drag Colliders
    function buildTrackballGizmo() {
      trackballGizmo = new THREE.Group();
      trackballGizmo.visible = isGizmoVisible;

      // 1. Silky Smooth Translucent Glass Sphere (No wireframe clutter!)
      const sphereGeo = new THREE.SphereGeometry(1.0, 64, 48);
      const sphereMat = new THREE.MeshBasicMaterial({
        color: 0x38bdf8,
        transparent: true,
        opacity: sphereOpacityValue,
        depthTest: false,
        wireframe: false
      });
      sphereMesh = new THREE.Mesh(sphereGeo, sphereMat);
      trackballGizmo.add(sphereMesh);

      // Helper function to build an axis ring with visual mesh and invisible thick hit collider
      function createAxisRing(axis, colorHex, rotEuler) {
        // Sleek High-Res Visible Ring
        const ringGeo = new THREE.TorusGeometry(1.0, 0.012, 16, 96);
        const ringMat = new THREE.MeshBasicMaterial({
          color: colorHex,
          depthTest: false,
          transparent: true,
          opacity: 0.88
        });
        const rMesh = new THREE.Mesh(ringGeo, ringMat);
        rMesh.rotation.copy(rotEuler);
        trackballGizmo.add(rMesh);
        ringMeshes[axis] = rMesh;

        // Thick Invisible Collider for easy mouse grabbing
        const colGeo = new THREE.TorusGeometry(1.0, 0.09, 10, 48);
        const colMat = new THREE.MeshBasicMaterial({
          visible: false,
          depthTest: false
        });
        const colMesh = new THREE.Mesh(colGeo, colMat);
        colMesh.rotation.copy(rotEuler);
        colMesh.userData = { axis: axis, baseColor: colorHex };
        trackballGizmo.add(colMesh);
        ringColliders.push(colMesh);
      }

      // 2. Pitch Ring (X-Axis, Red #ef4444, Normal +X)
      createAxisRing("x", 0xef4444, new THREE.Euler(0, Math.PI / 2, 0));

      // 3. Roll Ring (Y-Axis, Green #22c55e, Normal +Y)
      createAxisRing("y", 0x22c55e, new THREE.Euler(-Math.PI / 2, 0, 0));

      // 4. Yaw Ring (Z-Axis, Blue #3b82f6, Normal +Z)
      createAxisRing("z", 0x3b82f6, new THREE.Euler(0, 0, 0));

      scene.add(trackballGizmo);
    }

    // Build Geometry with Shader Support for Both Crisp Pixels and Continuous Splats
    function buildSplatGeometry() {
      geometry = new THREE.BufferGeometry();

      positionsAttr = new THREE.BufferAttribute(new Float32Array(rawPositions), 3);
      colorsAttr = new THREE.BufferAttribute(new Float32Array(numGaussians * 3), 3);
      alphasAttr = new THREE.BufferAttribute(new Float32Array(numGaussians), 1);
      sizesAttr = new THREE.BufferAttribute(new Float32Array(numGaussians), 1);

      geometry.setAttribute("position", positionsAttr);
      geometry.setAttribute("customColor", colorsAttr);
      geometry.setAttribute("customAlpha", alphasAttr);
      geometry.setAttribute("customSize", sizesAttr);

      const vertexShader = `
        attribute vec3 customColor;
        attribute float customAlpha;
        attribute float customSize;
        uniform float uPointSizePx;
        uniform float uRenderMode;
        varying vec3 vColor;
        varying float vAlpha;
        void main() {
          if (customAlpha < 0.01) {
            gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
            gl_PointSize = 0.0;
            return;
          }
          vColor = customColor;
          vAlpha = customAlpha;
          vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
          gl_Position = projectionMatrix * mvPosition;

          if (uRenderMode < 0.5) {
            // Crisp Metashape Points: sharp screen pixels with depth attenuation
            float depthFactor = clamp(1.8 / -mvPosition.z, 0.35, 4.0);
            gl_PointSize = uPointSizePx * depthFactor;
          } else {
            // Gaussian Splat Mode: realistic continuous surfels scaling with camera distance
            float baseRadius = customSize;
            float worldRadius = baseRadius * (uPointSizePx * 0.25);
            float projSize = worldRadius * (650.0 / -mvPosition.z);
            gl_PointSize = clamp(projSize, 1.2, 54.0);
          }
        }
      `;

      const fragmentShader = `
        uniform float uRenderMode;
        uniform float uPointShape;
        varying vec3 vColor;
        varying float vAlpha;
        void main() {
          if (vAlpha < 0.01) discard;

          if (uRenderMode < 0.5) {
            if (uPointShape < 0.5) {
              vec2 coord = gl_PointCoord - vec2(0.5);
              float distSq = dot(coord, coord);
              if (distSq > 0.25) discard;
              float edge = smoothstep(0.25, 0.20, distSq);
              gl_FragColor = vec4(vColor, vAlpha * edge);
            } else {
              gl_FragColor = vec4(vColor, vAlpha);
            }
          } else {
            // Gaussian Splat Mode: smooth radial exponential falloff without giant circle cutoffs
            vec2 coord = gl_PointCoord - vec2(0.5);
            float r2 = dot(coord, coord) * 4.0;
            if (r2 > 1.0) discard;
            float gFactor = exp(-2.8 * r2);
            float alpha = clamp(vAlpha * gFactor, 0.0, 1.0);
            if (alpha < 0.035) discard;
            gl_FragColor = vec4(vColor, alpha);
          }
        }
      `;

      material = new THREE.ShaderMaterial({
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
          uPointSizePx: { value: pointSizeValue },
          uRenderMode: { value: currentRenderMode === "points" ? 0.0 : 1.0 },
          uPointShape: { value: pointShapeValue }
        },
        transparent: true,
        depthTest: true,
        depthWrite: currentRenderMode === "points",
        blending: THREE.NormalBlending
      });

      pointsMesh = new THREE.Points(geometry, material);
      modelContainer.add(pointsMesh);
      updateSplatBuffers();
    }

    // -------------------------------------------------------------------------
    // RealityCapture Spatial Density Precomputation (3D Spatial Hash Grid)
    // -------------------------------------------------------------------------
    function computeSpatialDensities() {
      if (numGaussians === 0) return;
      pointDensities = new Float32Array(numGaussians);

      // Fine-grained voxel cell size (3.5 cm) for sharp separation of floating outlier noise
      const cellSize = Math.max(0.008, Math.min(0.035, targetDim * 0.012));
      const invCell = 1.0 / cellSize;

      const grid = new Map();
      for (let i = 0; i < numGaussians; i++) {
        const gx = Math.floor(rawPositions[i * 3] * invCell);
        const gy = Math.floor(rawPositions[i * 3 + 1] * invCell);
        const gz = Math.floor(rawPositions[i * 3 + 2] * invCell);
        const k = `${gx},${gy},${gz}`;
        grid.set(k, (grid.get(k) || 0) + 1);
      }

      const cellDensity = new Map();
      for (let [key, count] of grid.entries()) {
        const parts = key.split(",");
        const cx = parseInt(parts[0], 10);
        const cy = parseInt(parts[1], 10);
        const cz = parseInt(parts[2], 10);
        let sum = 0;
        for (let dx = -1; dx <= 1; dx++) {
          for (let dy = -1; dy <= 1; dy++) {
            for (let dz = -1; dz <= 1; dz++) {
              const neighborKey = `${cx + dx},${cy + dy},${cz + dz}`;
              const c = grid.get(neighborKey);
              if (c) sum += c;
            }
          }
        }
        cellDensity.set(key, sum);
      }

      maxDensityVal = 1;
      minDensityVal = 999999;
      for (let i = 0; i < numGaussians; i++) {
        const gx = Math.floor(rawPositions[i * 3] * invCell);
        const gy = Math.floor(rawPositions[i * 3 + 1] * invCell);
        const gz = Math.floor(rawPositions[i * 3 + 2] * invCell);
        const k = `${gx},${gy},${gz}`;
        const d = cellDensity.get(k) || 1;
        pointDensities[i] = d;
        if (d > maxDensityVal) maxDensityVal = d;
        if (d < minDensityVal) minDensityVal = d;
      }

      // Pre-sort densities for instant, exact O(1) percentile queries
      sortedDensities = Float32Array.from(pointDensities).sort();

      densityPercentiles = {
        p10: sortedDensities[Math.floor(numGaussians * 0.10)] || 1,
        p25: sortedDensities[Math.floor(numGaussians * 0.25)] || 5,
        p50: sortedDensities[Math.floor(numGaussians * 0.50)] || 15,
        p75: sortedDensities[Math.floor(numGaussians * 0.75)] || 35,
        p90: sortedDensities[Math.floor(numGaussians * 0.90)] || 60
      };
    }

    function computeDistanceDistribution() {
      if (numGaussians === 0) return;
      const dists = new Float32Array(numGaussians);
      for (let i = 0; i < numGaussians; i++) {
        const x = rawPositions[i * 3];
        const y = rawPositions[i * 3 + 1];
        const z = rawPositions[i * 3 + 2];
        dists[i] = Math.sqrt(x * x + y * y + z * z);
      }
      sortedDistances = dists.sort();
    }

    // -------------------------------------------------------------------------
    // Online Noise Filter Controls
    // -------------------------------------------------------------------------
    function updateNoiseFilter(val) {
      noiseFilterPercent = parseInt(val, 10);
      const lbl = document.getElementById("lbl-noise-val");
      if (lbl) lbl.textContent = noiseFilterPercent === 0 ? "Off" : `${noiseFilterPercent}%`;

      [0, 15, 35, 65].forEach((preset, idx) => {
        const btn = document.getElementById(`btn-noise-${idx}`);
        if (btn) btn.classList.toggle("active", noiseFilterPercent === preset);
      });

      updateSplatBuffers();
    }

    function setNoiseFilterPreset(pct) {
      const sl = document.getElementById("slider-noise");
      if (sl) sl.value = pct;
      updateNoiseFilter(pct);
      showToast(pct === 0 ? "Noise filter disabled" : `Noise filter: ${pct}% active`);
    }

    // -------------------------------------------------------------------------
    // RealityCapture Reconstruction Region (3D Section Box)
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // RealityCapture Reconstruction Region (3D Section Box with 6-DOF & OBB)
    // -------------------------------------------------------------------------
    function toggleClipBox() {
      clipBoxEnabled = !clipBoxEnabled;
      const btn = document.getElementById("btn-toggle-clipbox");
      const ctrl = document.getElementById("clipbox-controls");
      if (btn) {
        btn.classList.toggle("active", clipBoxEnabled);
        btn.textContent = clipBoxEnabled ? "Disable Box" : "Enable Box";
      }
      if (ctrl) ctrl.style.display = clipBoxEnabled ? "flex" : "none";

      if (clipBoxEnabled && !clipBoxGroup) {
        initClipBox();
      }
      if (clipBoxGroup) clipBoxGroup.visible = clipBoxEnabled;
      if (transformControls) {
        if (clipBoxEnabled) {
          transformControls.attach(clipBoxGroup);
        } else {
          transformControls.detach();
        }
      }
      updateSplatBuffers();
      showToast(clipBoxEnabled ? "3D Reconstruction Box enabled" : "3D Reconstruction Box disabled");
    }

    function initClipBox() {
      if (clipBoxGroup) return;

      clipBoxGroup = new THREE.Group();
      clipBoxGroup.name = "ReconstructionRegion";

      // Calculate initial dimensions based on subject
      const initDim = targetDim * 1.35;
      clipBoxGroup.scale.set(initDim, initDim, initDim);
      clipBoxGroup.position.set(0, 0, 0);

      // 1. Semi-transparent volume faces (clean translucent cyan)
      const geom = new THREE.BoxGeometry(1, 1, 1);
      const faceMat = new THREE.MeshBasicMaterial({
        color: 0x38bdf8,
        transparent: true,
        opacity: 0.08,
        side: THREE.DoubleSide,
        depthWrite: false
      });
      clipBoxMesh = new THREE.Mesh(geom, faceMat);
      clipBoxGroup.add(clipBoxMesh);

      // 2. Thick cyan wireframe edges (RealityCapture style)
      const edges = new THREE.EdgesGeometry(geom);
      const lineMat = new THREE.LineBasicMaterial({
        color: 0x38bdf8,
        linewidth: 2,
        transparent: true,
        opacity: 0.88
      });
      clipBoxLines = new THREE.LineSegments(edges, lineMat);
      clipBoxGroup.add(clipBoxLines);

      modelContainer.add(clipBoxGroup);

      // 3. Initialize TransformControls if available
      if (typeof THREE.TransformControls !== "undefined") {
        transformControls = new THREE.TransformControls(camera, renderer.domElement);
        transformControls.size = 0.85;
        transformControls.setMode(gizmoMode);
        transformControls.attach(clipBoxGroup);

        transformControls.addEventListener("change", () => {
          syncUIFromBox();
          updateSplatBuffers();
        });

        transformControls.addEventListener("dragging-changed", (event) => {
          controls.enabled = !event.value; // Prevent orbit camera while dragging gizmo
        });

        scene.add(transformControls);
      }

      syncUIFromBox();
    }

    function setGizmoMode(mode) {
      gizmoMode = mode;
      const bTrans = document.getElementById("btn-gizmo-trans");
      const bRot = document.getElementById("btn-gizmo-rot");
      const bScale = document.getElementById("btn-gizmo-scale");
      if (bTrans) bTrans.classList.toggle("active", mode === "translate");
      if (bRot) bRot.classList.toggle("active", mode === "rotate");
      if (bScale) bScale.classList.toggle("active", mode === "scale");

      if (transformControls) {
        transformControls.setMode(mode);
        showToast(mode === "translate" ? "Gizmo: Move Box (W)" : (mode === "rotate" ? "Gizmo: Rotate Box (E)" : "Gizmo: Scale Box (R)"));
      }
    }

    function syncUIFromBox() {
      if (!clipBoxGroup) return;

      const rx = Math.round(clipBoxGroup.rotation.x * 180 / Math.PI);
      const ry = Math.round(clipBoxGroup.rotation.y * 180 / Math.PI);
      const rz = Math.round(clipBoxGroup.rotation.z * 180 / Math.PI);

      const rotLabel = document.getElementById("lbl-box-rot");
      if (rotLabel) rotLabel.textContent = `${rx}&deg; / ${ry}&deg; / ${rz}&deg;`;

      const rotX = document.getElementById("clip-rot-x");
      const rotY = document.getElementById("clip-rot-y");
      const rotZ = document.getElementById("clip-rot-z");
      if (rotX) rotX.value = rx;
      if (rotY) rotY.value = ry;
      if (rotZ) rotZ.value = rz;

      const sx = clipBoxGroup.scale.x;
      const sy = clipBoxGroup.scale.y;
      const sz = clipBoxGroup.scale.z;

      const dimsLabel = document.getElementById("lbl-box-dims");
      if (dimsLabel) dimsLabel.textContent = `${sx.toFixed(2)} x ${sy.toFixed(2)} x ${sz.toFixed(2)}m`;

      const sizeX = document.getElementById("clip-size-x");
      const sizeY = document.getElementById("clip-size-y");
      const sizeZ = document.getElementById("clip-size-z");
      const baseDim = Math.max(0.01, targetDim);
      if (sizeX) sizeX.value = Math.round((sx / baseDim) * 100);
      if (sizeY) sizeY.value = Math.round((sy / baseDim) * 100);
      if (sizeZ) sizeZ.value = Math.round((sz / baseDim) * 100);
    }

    function updateBoxFromUI() {
      if (!clipBoxGroup) initClipBox();

      const rx = (parseFloat(document.getElementById("clip-rot-x").value) || 0) * Math.PI / 180;
      const ry = (parseFloat(document.getElementById("clip-rot-y").value) || 0) * Math.PI / 180;
      const rz = (parseFloat(document.getElementById("clip-rot-z").value) || 0) * Math.PI / 180;
      clipBoxGroup.rotation.set(rx, ry, rz);

      const baseDim = Math.max(0.01, targetDim);
      const sx = ((parseFloat(document.getElementById("clip-size-x").value) || 100) / 100.0) * baseDim;
      const sy = ((parseFloat(document.getElementById("clip-size-y").value) || 100) / 100.0) * baseDim;
      const sz = ((parseFloat(document.getElementById("clip-size-z").value) || 100) / 100.0) * baseDim;
      clipBoxGroup.scale.set(Math.max(0.05, sx), Math.max(0.05, sy), Math.max(0.05, sz));

      const rotLabel = document.getElementById("lbl-box-rot");
      if (rotLabel) rotLabel.textContent = `${Math.round(rx*180/Math.PI)}&deg; / ${Math.round(ry*180/Math.PI)}&deg; / ${Math.round(rz*180/Math.PI)}&deg;`;

      const dimsLabel = document.getElementById("lbl-box-dims");
      if (dimsLabel) dimsLabel.textContent = `${sx.toFixed(2)} x ${sy.toFixed(2)} x ${sz.toFixed(2)}m`;

      updateSplatBuffers();
    }

    function autoFitClipBox() {
      if (!clipBoxGroup) initClipBox();
      if (!clipBoxEnabled) toggleClipBox();

      clipBoxGroup.position.set(0, 0, 0);
      clipBoxGroup.rotation.set(0, 0, 0);
      const fitDim = targetDim * 1.15;
      clipBoxGroup.scale.set(fitDim, fitDim, fitDim * 1.3);

      syncUIFromBox();
      updateSplatBuffers();
      showToast("&#127793; Reconstruction box auto-fitted to subject core");
    }

    function cutClipBoxBase() {
      if (!clipBoxGroup) initClipBox();
      if (!clipBoxEnabled) toggleClipBox();

      // Lift box center up slightly and reduce height to trim floor/table
      clipBoxGroup.position.z += clipBoxGroup.scale.z * 0.15;
      clipBoxGroup.scale.z *= 0.75;

      syncUIFromBox();
      updateSplatBuffers();
      showToast("&#9986;&#65039; Floor / table surface trimmed");
    }

    function levelClipBox() {
      if (!clipBoxGroup) return;
      // Reset pitch and roll, preserve yaw
      const currentYaw = clipBoxGroup.rotation.z;
      clipBoxGroup.rotation.set(0, 0, currentYaw);

      syncUIFromBox();
      updateSplatBuffers();
      showToast("&#9878;&#65039; Box leveled to horizontal horizon plane");
    }

    function resetClipBox() {
      if (!clipBoxGroup) initClipBox();
      clipBoxGroup.position.set(0, 0, 0);
      clipBoxGroup.rotation.set(0, 0, 0);
      const initDim = targetDim * 1.35;
      clipBoxGroup.scale.set(initDim, initDim, initDim);

      syncUIFromBox();
      updateSplatBuffers();
      showToast("&#128260; Reconstruction box reset");
    }

    function rebuildCroppedGeometry() {
      if (splatMesh) {
        modelContainer.remove(splatMesh);
        if (geometry) geometry.dispose();
      }
      buildSplatGeometry();
      updateSplatBuffers();
    }

    function applyCropPointCloud() {
      if (!clipBoxEnabled || !clipBoxGroup) {
        showToast("Please enable the 3D Box first");
        return;
      }

      clipBoxGroup.updateMatrixWorld();
      const invMat = new THREE.Matrix4().copy(clipBoxGroup.matrixWorld).invert();
      const m = invMat.elements;

      // Save initial state before first crop for Undo
      if (!initialUncroppedState) {
        initialUncroppedState = {
          positions: new Float32Array(rawPositions),
          colors: new Float32Array(rawColorsRGB),
          opacities: new Float32Array(rawOpacities),
          scales: rawScales ? new Float32Array(rawScales) : null,
          semantics: rawSemantics ? new Int32Array(rawSemantics) : null,
          densities: (pointDensities && pointDensities.length > 0) ? new Float32Array(pointDensities) : null,
          count: numGaussians
        };
      }

      const keptIndices = [];
      for (let i = 0; i < numGaussians; ++i) {
        const px = rawPositions[i * 3];
        const py = rawPositions[i * 3 + 1];
        const pz = rawPositions[i * 3 + 2];

        const lx = m[0] * px + m[4] * py + m[8] * pz + m[12];
        const ly = m[1] * px + m[5] * py + m[9] * pz + m[13];
        const lz = m[2] * px + m[6] * py + m[10] * pz + m[14];

        if (Math.abs(lx) <= 0.5 && Math.abs(ly) <= 0.5 && Math.abs(lz) <= 0.5) {
          keptIndices.push(i);
        }
      }

      const newCount = keptIndices.length;
      if (newCount === 0) {
        showToast("&#9888;&#65039; The box contains 0 points! Adjust box position to encompass your subject.");
        return;
      }

      const newPos = new Float32Array(newCount * 3);
      const newCol = new Float32Array(newCount * 3);
      const newOpac = new Float32Array(newCount);
      const newScales = rawScales ? new Float32Array(newCount) : null;
      const newSem = rawSemantics ? new Int32Array(newCount) : null;
      const newDensities = (pointDensities && pointDensities.length > 0) ? new Float32Array(newCount) : null;

      for (let j = 0; j < newCount; ++j) {
        const origIdx = keptIndices[j];
        newPos[j * 3]     = rawPositions[origIdx * 3];
        newPos[j * 3 + 1] = rawPositions[origIdx * 3 + 1];
        newPos[j * 3 + 2] = rawPositions[origIdx * 3 + 2];

        newCol[j * 3]     = rawColorsRGB[origIdx * 3];
        newCol[j * 3 + 1] = rawColorsRGB[origIdx * 3 + 1];
        newCol[j * 3 + 2] = rawColorsRGB[origIdx * 3 + 2];

        newOpac[j] = rawOpacities[origIdx];
        if (newScales) newScales[j] = rawScales[origIdx];
        if (newSem) newSem[j] = rawSemantics[origIdx];
        if (newDensities) newDensities[j] = pointDensities[origIdx];
      }

      rawPositions = newPos;
      rawColorsRGB = newCol;
      rawOpacities = newOpac;
      if (newScales) rawScales = newScales;
      if (newSem) rawSemantics = newSem;
      if (newDensities) pointDensities = newDensities;
      numGaussians = newCount;

      rebuildCroppedGeometry();

      const btnUndo = document.getElementById("btn-undo-crop");
      if (btnUndo) btnUndo.style.display = "inline-flex";

      const ptBadge = document.querySelector(".layer-item .layer-meta");
      if (ptBadge) ptBadge.textContent = `${newCount.toLocaleString()} pts (Cropped)`;

      showToast(`&#9986;&#65039; Point cloud cropped! Kept ${newCount.toLocaleString()} points inside box.`);
    }

    function undoCropPointCloud() {
      if (!initialUncroppedState) return;

      rawPositions = new Float32Array(initialUncroppedState.positions);
      rawColorsRGB = new Float32Array(initialUncroppedState.colors);
      rawOpacities = new Float32Array(initialUncroppedState.opacities);
      if (initialUncroppedState.scales) rawScales = new Float32Array(initialUncroppedState.scales);
      if (initialUncroppedState.semantics) rawSemantics = new Int32Array(initialUncroppedState.semantics);
      if (initialUncroppedState.densities) pointDensities = new Float32Array(initialUncroppedState.densities);
      numGaussians = initialUncroppedState.count;

      rebuildCroppedGeometry();

      const btnUndo = document.getElementById("btn-undo-crop");
      if (btnUndo) btnUndo.style.display = "none";

      const ptBadge = document.querySelector(".layer-item .layer-meta");
      if (ptBadge) ptBadge.textContent = `${numGaussians.toLocaleString()} pts`;

      showToast(`&#8634; Restored original point cloud (${numGaussians.toLocaleString()} points)`);
    }

    function exportCroppedPLY() {
      if (numGaussians === 0) {
        showToast("No points available to export");
        return;
      }

      showToast("&#9203; Generating clean PLY point cloud file...");

      const cx = (centerOrigin && centerOrigin.length >= 3) ? centerOrigin[0] : 0;
      const cy = (centerOrigin && centerOrigin.length >= 3) ? centerOrigin[1] : 0;
      const nl = String.fromCharCode(10);
      const header = [
        "ply",
        "format ascii 1.0",
        "comment Exported by pliman 3D Reconstruction Box",
        "element vertex " + numGaussians,
        "property float x",
        "property float y",
        "property float z",
        "property uchar red",
        "property uchar green",
        "property uchar blue",
        "end_header",
        ""
      ].join(nl);

      const chunks = [header];
      let chunk = "";
      for (let i = 0; i < numGaussians; ++i) {
        const x = (rawPositions[i * 3] + cx).toFixed(5);
        const y = (rawPositions[i * 3 + 1] + cy).toFixed(5);
        const z = (rawPositions[i * 3 + 2] + cz).toFixed(5);
        const r = Math.min(255, Math.max(0, Math.round(rawColorsRGB[i * 3] * 255)));
        const g = Math.min(255, Math.max(0, Math.round(rawColorsRGB[i * 3 + 1] * 255)));
        const b = Math.min(255, Math.max(0, Math.round(rawColorsRGB[i * 3 + 2] * 255)));

        chunk += x + " " + y + " " + z + " " + r + " " + g + " " + b + nl;
        if (chunk.length > 65536) {
          chunks.push(chunk);
          chunk = "";
        }
      }
      if (chunk.length > 0) chunks.push(chunk);

      const blob = new Blob(chunks, { type: "application/octet-stream" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `point_cloud_cropped_${numGaussians}pts.ply`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);

      showToast(`&#128190; Saved point_cloud_cropped_${numGaussians}pts.ply!`);
    }

    let cropRadiusMultiplier = 1.0;
    function setCropPreset(val) {
      cropRadiusMultiplier = parseFloat(val);
      const slider = document.getElementById("slider-crop");
      if (slider) slider.value = val;
      const lbl = document.getElementById("lbl-crop-val");
      if (lbl) lbl.textContent = Math.round(cropRadiusMultiplier * 100) + "%";
      updateSplatBuffers();
      showToast(val < 0.9 ? "&#127793; Subject focus active (distant noise cropped)" : "&#127760; Full scene restored");
    }

    function updateCropRadius(val) {
      cropRadiusMultiplier = parseFloat(val);
      const lbl = document.getElementById("lbl-crop-val");
      if (lbl) lbl.textContent = Math.round(cropRadiusMultiplier * 100) + "%";
      updateSplatBuffers();
    }

    function updateSplatBuffers() {
      const colArr = colorsAttr.array;
      const alphaArr = alphasAttr.array;
      const sizeArr = sizesAttr.array;

      const zMin = -targetDim * 0.5;
      const zRange = Math.max(0.001, targetDim);
      let visibleCount = 0;
      let filteredNoiseCount = 0;

      // Bounding Region distance threshold
      let maxCropDist = Infinity;
      if (cropRadiusMultiplier < 0.99) {
        if (sortedDistances && sortedDistances.length === numGaussians) {
          const frac = Math.min(1.0, cropRadiusMultiplier * 1.35);
          const distIdx = Math.min(numGaussians - 1, Math.max(0, Math.floor(frac * numGaussians)));
          maxCropDist = sortedDistances[distIdx];
        } else {
          maxCropDist = (targetDim * 2.5) * cropRadiusMultiplier;
        }
      }
      const cropSq = maxCropDist * maxCropDist;

      // Online noise filter cutoff threshold (Exact percentile indexing from sorted densities)
      let densityCutoff = 0;
      if (noiseFilterPercent > 0 && sortedDensities && sortedDensities.length === numGaussians) {
        const frac = noiseFilterPercent / 100.0;
        const pIdx = Math.min(numGaussians - 1, Math.max(0, Math.floor(frac * numGaussians)));
        densityCutoff = sortedDensities[pIdx];
      }

      // Pre-calculate inverse world transform matrix of the Reconstruction Box
      let invBoxMatrix = null;
      if (clipBoxEnabled && clipBoxGroup) {
        clipBoxGroup.updateMatrixWorld();
        invBoxMatrix = new THREE.Matrix4().copy(clipBoxGroup.matrixWorld).invert();
      }

      for (let i = 0; i < numGaussians; ++i) {
        const px = rawPositions[i * 3];
        const py = rawPositions[i * 3 + 1];
        const pz = rawPositions[i * 3 + 2];
        const distSq = px * px + py * py + pz * pz;

        // Metashape Bounding Region: crop distant room/wall noise
        if (cropRadiusMultiplier < 0.99 && distSq > cropSq) {
          alphaArr[i] = 0.0;
          sizeArr[i] = 0.0;
          continue;
        }

        // RealityCapture Reconstruction Region (3D Oriented Box)
        if (clipBoxEnabled && invBoxMatrix) {
          const m = invBoxMatrix.elements;
          const lx = m[0] * px + m[4] * py + m[8] * pz + m[12];
          const ly = m[1] * px + m[5] * py + m[9] * pz + m[13];
          const lz = m[2] * px + m[6] * py + m[10] * pz + m[14];

          if (Math.abs(lx) > 0.5 || Math.abs(ly) > 0.5 || Math.abs(lz) > 0.5) {
            alphaArr[i] = 0.0;
            sizeArr[i] = 0.0;
            continue;
          }
        }

        // Online Noise Filter (Spatial Density Check)
        if (densityCutoff > 0 && pointDensities[i] < densityCutoff) {
          alphaArr[i] = 0.0;
          sizeArr[i] = 0.0;
          filteredNoiseCount++;
          continue;
        }

        const sem = rawSemantics[i];
        const isClassActive = activeSemanticClasses.has(sem);
        const opac = rawOpacities[i] * globalOpacityMultiplier;

        if (!isClassActive || opac < minAlphaCutoff) {
          alphaArr[i] = 0.0;
          sizeArr[i] = 0.0;
          continue;
        }

        visibleCount++;
        alphaArr[i] = Math.min(1.0, opac);

        // Realistic Gaussian scale in meters (bounded between 0.8mm and 3.5cm)
        const sc = (rawScales && rawScales.length > i) ? rawScales[i] : 0.01;
        sizeArr[i] = Math.min(0.035, Math.max(0.0008, sc));

        if (currentColorMode === "rgb") {
          colArr[i * 3 + 0] = rawColorsRGB[i * 3 + 0];
          colArr[i * 3 + 1] = rawColorsRGB[i * 3 + 1];
          colArr[i * 3 + 2] = rawColorsRGB[i * 3 + 2];
        } else if (currentColorMode === "clay") {
          // Studio Clay shading
          const shade = 0.82 + 0.18 * Math.sin((pz - zMin) / zRange * Math.PI);
          colArr[i * 3 + 0] = 0.76 * shade;
          colArr[i * 3 + 1] = 0.75 * shade;
          colArr[i * 3 + 2] = 0.74 * shade;
        } else if (currentColorMode === "density") {
          // Point Density Heatmap
          const dVal = pointDensities[i] || 1;
          const dNorm = Math.min(1.0, Math.max(0.0, (dVal - densityPercentiles.p10) / (densityPercentiles.p90 - densityPercentiles.p10 || 1)));
          colArr[i * 3 + 0] = Math.min(1.0, Math.max(0.0, 1.6 * dNorm - 0.2));
          colArr[i * 3 + 1] = Math.min(1.0, Math.max(0.0, 1.0 - 2.0 * Math.abs(dNorm - 0.5)));
          colArr[i * 3 + 2] = Math.min(1.0, Math.max(0.0, 1.4 * (1.0 - dNorm) - 0.1));
        } else if (currentColorMode === "semantic") {
          const c = semanticColorsPalette[sem] || semanticColorsPalette[0];
          colArr[i * 3 + 0] = c[0];
          colArr[i * 3 + 1] = c[1];
          colArr[i * 3 + 2] = c[2];
        } else if (currentColorMode === "height") {
          const zNorm = Math.min(1.0, Math.max(0.0, (rawPositions[i * 3 + 2] - zMin) / zRange));
          colArr[i * 3 + 0] = Math.sin(zNorm * Math.PI * 1.5) * 0.5 + 0.5;
          colArr[i * 3 + 1] = Math.sin(zNorm * Math.PI) * 0.8 + 0.2;
          colArr[i * 3 + 2] = Math.cos(zNorm * Math.PI * 0.8) * 0.5 + 0.5;
        }
      }

      colorsAttr.needsUpdate = true;
      alphasAttr.needsUpdate = true;
      sizesAttr.needsUpdate = true;

      const visElem = document.getElementById("lbl-visible-splats");
      if (visElem) visElem.textContent = visibleCount.toLocaleString();

      const lblNoise = document.getElementById("lbl-noise-stats");
      if (lblNoise) {
        if (filteredNoiseCount > 0) {
          const pct = ((filteredNoiseCount / numGaussians) * 100).toFixed(1);
          lblNoise.textContent = `${filteredNoiseCount.toLocaleString()} noise points (${pct}%)`;
          lblNoise.style.color = "#facc15";
        } else {
          lblNoise.textContent = "0 noise points filtered";
          lblNoise.style.color = "var(--text-muted)";
        }
      }
    }

    // Setup interactive mouse handling for 3D coordinates, camera clicking, and gizmo dragging
    let hoveredCamIndex = null;
    let downPointerPos = { x: 0, y: 0 };

    function setupInteractions() {
      const coordsText = document.getElementById("coords-text");
      const coordsStatus = document.getElementById("coords-status");
      const canvasElem = renderer.domElement;
      const camTooltip = document.getElementById("cam-tooltip");
      raycaster.params.Points.threshold = Math.max(0.015, targetDim * 0.02);

      // Pointer Move
      window.addEventListener("pointermove", (event) => {
        mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

        // If dragging an axis ring directly in the 3D viewport:
        if (isDraggingRing && activeDragAxis) {
          const activeCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
          const centerPoint = trackballGizmo.position;

          if (!dragIsEdgeOn) {
            raycaster.setFromCamera(mouse, activeCam);
            const currentHit3D = new THREE.Vector3();
            if (raycaster.ray.intersectPlane(dragPlane, currentHit3D)) {
              const currentVector = currentHit3D.sub(centerPoint).normalize();
              const cross = new THREE.Vector3().crossVectors(dragStartVector, currentVector);
              const sinVal = cross.dot(dragAxisVector);
              const cosVal = Math.max(-1.0, Math.min(1.0, dragStartVector.dot(currentVector)));
              const deltaAngle = Math.atan2(sinVal, cosVal);

              const qRot = new THREE.Quaternion().setFromAxisAngle(dragAxisVector, deltaAngle);
              if (alignmentSpace === "local") {
                modelContainer.quaternion.copy(dragStartQuaternion).multiply(qRot);
              } else {
                modelContainer.quaternion.copy(qRot).multiply(dragStartQuaternion);
              }
            } else {
              dragIsEdgeOn = true;
            }
          }

          if (dragIsEdgeOn) {
            const dx = event.clientX - lastDragMousePos.x;
            const dy = event.clientY - lastDragMousePos.y;
            lastDragMousePos = { x: event.clientX, y: event.clientY };
            const proj = dx * dragEdgeTangent2D.x + dy * dragEdgeTangent2D.y;
            const incAngle = proj * 0.015;
            const qInc = new THREE.Quaternion().setFromAxisAngle(dragAxisVector, incAngle);
            if (alignmentSpace === "local") {
              modelContainer.quaternion.multiply(qInc);
            } else {
              modelContainer.quaternion.premultiply(qInc);
            }
          }

          modelContainer.rotation.setFromQuaternion(modelContainer.quaternion, "XYZ");
          syncSlidersFromModel();
          if (alignmentSpace === "local" && trackballGizmo) {
            trackballGizmo.quaternion.copy(modelContainer.quaternion);
          } else if (trackballGizmo) {
            trackballGizmo.quaternion.set(0, 0, 0, 1);
          }
          return;
        }

        if (event.target !== canvasElem) {
          if (camTooltip) camTooltip.style.display = "none";
          resetRingHighlights();
          resetCameraHighlights();
          return;
        }

        // Check if mouse is hovering over any axis ring collider
        const activeHoverCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
        if (trackballGizmo && trackballGizmo.visible) {
          raycaster.setFromCamera(mouse, activeHoverCam);
          const ringHits = raycaster.intersectObjects(ringColliders);
          if (ringHits.length > 0) {
            canvasElem.style.cursor = "grab";
            const hitAxis = ringHits[0].object.userData.axis;
            highlightRing(hitAxis);
            if (camTooltip) camTooltip.style.display = "none";
            return;
          } else {
            resetRingHighlights();
          }
        }

        // Check if mouse is hovering over any 3D camera frustum
        if (camerasGroup && camerasGroup.visible && cameraColliders.length > 0) {
          raycaster.setFromCamera(mouse, camera);
          const camHits = raycaster.intersectObjects(cameraColliders);
          if (camHits.length > 0) {
            canvasElem.style.cursor = "pointer";
            const camIdx = camHits[0].object.userData.camIndex;
            highlightCamera(camIdx);

            if (camTooltip) {
              const camData = camHits[0].object.userData.camData;
              camTooltip.style.display = "block";
              camTooltip.style.left = event.clientX + "px";
              camTooltip.style.top = event.clientY + "px";
              camTooltip.textContent = `&#128247; ${camData.label} (Click to inspect)`;
            }
            return;
          } else {
            resetCameraHighlights();
            if (camTooltip) camTooltip.style.display = "none";
          }
        }

        // Check if 3D Ruler tool is active
        const activeCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
        if (isRulerActive) {
          canvasElem.style.cursor = "crosshair";
          raycaster.setFromCamera(mouse, activeCam);
          const rHits = pointsMesh ? raycaster.intersectObject(pointsMesh) : [];
          if (rHits.length > 0) {
            const hPt = rHits[0].point;
            if (rulerHoverMarker) {
              rulerHoverMarker.position.copy(hPt);
              rulerHoverMarker.visible = true;
            }
            if (rulerPointA) {
              if (rulerTempLine) rulerGroup.remove(rulerTempLine);
              const tGeo = new THREE.BufferGeometry().setFromPoints([rulerPointA, hPt]);
              const tMat = new THREE.LineDashedMaterial({ color: 0xfacc15, dashSize: 0.02, gapSize: 0.01, depthTest: false });
              rulerTempLine = new THREE.Line(tGeo, tMat);
              rulerTempLine.computeLineDistances();
              rulerGroup.add(rulerTempLine);
            }
          }
          return;
        }

        canvasElem.style.cursor = "default";

        // Raycasting for coordinates tracking
        raycaster.setFromCamera(mouse, activeCam);
        const intersects = pointsMesh ? raycaster.intersectObject(pointsMesh) : [];

        if (intersects.length > 0) {
          const hit = intersects[0].point;
          const worldX = hit.x + centerOrigin[0];
          const worldY = hit.y + centerOrigin[1];
          const worldZ = hit.z + centerOrigin[2];
          if (coordsStatus) {
            coordsStatus.textContent = "Point Hover";
            coordsStatus.style.color = "#10b981";
          }
          if (coordsText) {
            coordsText.innerHTML = `X: <span style="color:#10b981">${worldX.toFixed(3)}</span> | Y: <span style="color:#10b981">${worldY.toFixed(3)}</span> | Z: <span style="color:#38bdf8">${worldZ.toFixed(3)}</span>`;
          }
        } else {
          const groundPlane = new THREE.Plane(new THREE.Vector3(0, 0, 1), targetDim * 0.45);
          const hitPlane = new THREE.Vector3();
          raycaster.ray.intersectPlane(groundPlane, hitPlane);
          if (hitPlane && Math.abs(hitPlane.x) < targetDim * 2.0 && Math.abs(hitPlane.y) < targetDim * 2.0) {
            const worldX = hitPlane.x + centerOrigin[0];
            const worldY = hitPlane.y + centerOrigin[1];
            const worldZ = hitPlane.z + centerOrigin[2];
            if (coordsStatus) {
              coordsStatus.textContent = "Ground Plane";
              coordsStatus.style.color = "var(--text-muted)";
            }
            if (coordsText) {
              coordsText.innerHTML = `X: ${worldX.toFixed(3)} | Y: ${worldY.toFixed(3)} | Z: ${worldZ.toFixed(3)}`;
            }
          }
        }
      });

      // Pointer Down
      window.addEventListener("pointerdown", (event) => {
        if (event.target !== canvasElem) return;
        if (event.button !== 0) return;
        downPointerPos = { x: event.clientX, y: event.clientY };

        const activeCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
        if (trackballGizmo && trackballGizmo.visible) {
          raycaster.setFromCamera(mouse, activeCam);
          const ringHits = raycaster.intersectObjects(ringColliders);
          if (ringHits.length > 0) {
            isDraggingRing = true;
            activeDragAxis = ringHits[0].object.userData.axis;
            highlightRing(activeDragAxis);
            controls.enabled = false;
            canvasElem.style.cursor = "grabbing";

            dragStartQuaternion.copy(modelContainer.quaternion);
            lastDragMousePos = { x: event.clientX, y: event.clientY };

            const centerPoint = trackballGizmo.position;

            if (alignmentSpace === "local") {
              const localA = (activeDragAxis === "x") ? new THREE.Vector3(1, 0, 0) :
                             (activeDragAxis === "y") ? new THREE.Vector3(0, 1, 0) :
                                                        new THREE.Vector3(0, 0, 1);
              dragAxisVector.copy(localA).applyQuaternion(modelContainer.quaternion).normalize();
            } else {
              if (activeDragAxis === "x") dragAxisVector.set(1, 0, 0);
              else if (activeDragAxis === "y") dragAxisVector.set(0, 1, 0);
              else dragAxisVector.set(0, 0, 1);
            }

            const dotNormal = raycaster.ray.direction.dot(dragAxisVector);
            dragIsEdgeOn = Math.abs(dotNormal) < 0.12;

            if (!dragIsEdgeOn) {
              dragPlane.setFromNormalAndCoplanarPoint(dragAxisVector, centerPoint);
              const hit3D = new THREE.Vector3();
              if (raycaster.ray.intersectPlane(dragPlane, hit3D)) {
                dragStartVector.copy(hit3D).sub(centerPoint).normalize();
              } else {
                dragIsEdgeOn = true;
              }
            }

            if (dragIsEdgeOn) {
              const cScreen = toScreenCoords(centerPoint, activeCam);
              const aPoint = centerPoint.clone().add(dragAxisVector);
              const aScreen = toScreenCoords(aPoint, activeCam);
              const axis2D = new THREE.Vector2(aScreen.x - cScreen.x, aScreen.y - cScreen.y).normalize();
              dragEdgeTangent2D.set(-axis2D.y, axis2D.x);
            }
            return;
          }
        }
      });

      // Pointer Up: 3D Ruler and Camera click detection
      window.addEventListener("pointerup", (event) => {
        if (isDraggingRing) {
          isDraggingRing = false;
          activeDragAxis = null;
          controls.enabled = true;
          canvasElem.style.cursor = "default";
          resetRingHighlights();
          return;
        }

        if (event.target !== canvasElem) return;

        const dist = Math.hypot(event.clientX - downPointerPos.x, event.clientY - downPointerPos.y);
        const activeCam = (isOrtho && orthoCamera) ? orthoCamera : camera;

        // If ruler is active, record measurement point on click
        if (isRulerActive && dist < 6) {
          mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
          mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
          raycaster.setFromCamera(mouse, activeCam);
          const rHits = pointsMesh ? raycaster.intersectObject(pointsMesh) : [];
          if (rHits.length > 0) {
            handleRulerClick(rHits[0].point);
            return;
          }
        }

        // Camera click detection
        if (dist < 4 && camerasGroup && camerasGroup.visible && cameraColliders.length > 0) {
          mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
          mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
          raycaster.setFromCamera(mouse, activeCam);
          const camHits = raycaster.intersectObjects(cameraColliders);
          if (camHits.length > 0) {
            const camIdx = camHits[0].object.userData.camIndex;
            openPhotoModal(camIdx);
          }
        }
      });

      // Double-Click anywhere to set rotation pivot center
      window.addEventListener("dblclick", (event) => {
        if (event.target !== canvasElem) return;
        mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
        const activeCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
        raycaster.setFromCamera(mouse, activeCam);
        const intersects = pointsMesh ? raycaster.intersectObject(pointsMesh) : [];
        if (intersects.length > 0) {
          const newTarget = intersects[0].point;
          controls.target.copy(newTarget);
          controls.update();
        }
      });
    }

    function highlightRing(axis) {
      resetRingHighlights();
      if (ringMeshes[axis]) {
        ringMeshes[axis].material.color.setHex(0xfacc15);
        ringMeshes[axis].material.opacity = 1.0;
      }
    }

    function resetRingHighlights() {
      if (ringMeshes.x) { ringMeshes.x.material.color.setHex(0xef4444); ringMeshes.x.material.opacity = 0.88; }
      if (ringMeshes.y) { ringMeshes.y.material.color.setHex(0x22c55e); ringMeshes.y.material.opacity = 0.88; }
      if (ringMeshes.z) { ringMeshes.z.material.color.setHex(0x3b82f6); ringMeshes.z.material.opacity = 0.88; }
    }

    function highlightCamera(idx) {
      if (hoveredCamIndex === idx) return;
      resetCameraHighlights();
      hoveredCamIndex = idx;
      const m = cameraMeshes[idx];
      if (m) {
        m.lines.material.color.setHex(0xfacc15); // Glowing yellow
        m.faces.material.color.setHex(0xfacc15);
        m.faces.material.opacity = 0.35;
        if (m.sight) m.sight.material.color.setHex(0xfacc15);
      }
    }

    function resetCameraHighlights() {
      if (hoveredCamIndex !== null && cameraMeshes[hoveredCamIndex]) {
        const m = cameraMeshes[hoveredCamIndex];
        m.lines.material.color.setHex(0xef4444);
        m.faces.material.color.setHex(0xef4444);
        m.faces.material.opacity = 0.16;
        if (m.sight) m.sight.material.color.setHex(0xd946ef);
        hoveredCamIndex = null;
      }
    }

    // Photo Modal Handling (Dynamic On-Demand Photo Rendering)
    function openPhotoModal(idx) {
      activeCameraModalIndex = parseInt(idx, 10);
      if (isNaN(activeCameraModalIndex)) activeCameraModalIndex = 0;
      const cam = rawCameras[activeCameraModalIndex];
      if (!cam) return;

      const modal = document.getElementById("photo-modal");
      const modalName = document.getElementById("modal-cam-name");
      const modalMeta = document.getElementById("modal-cam-meta");
      const modalImg = document.getElementById("modal-img");
      const modalLoading = document.getElementById("modal-img-loading");
      const modalFallback = document.getElementById("modal-photo-fallback");

      modalName.textContent = cam.label;
      modalMeta.textContent = `Camera ${activeCameraModalIndex + 1} of ${rawCameras.length} &bull; Resolution: ${cam.w}x${cam.h} &bull; fx: ${cam.fx}px`;

      let targetSrc = cam.img || "";
      if (!targetSrc && typeof localPhotoStore !== "undefined") {
        const lbl = (cam.label || "").toLowerCase();
        const base = (lbl.indexOf(String.fromCharCode(92)) !== -1 ? lbl.split(String.fromCharCode(92)) : lbl.split("/")).pop();
        if (localPhotoStore[lbl]) targetSrc = localPhotoStore[lbl];
        else if (localPhotoStore[base]) targetSrc = localPhotoStore[base];
      }

      if (targetSrc && targetSrc.length > 0) {
        if (modalLoading) modalLoading.style.display = "flex";
        if (modalFallback) modalFallback.style.display = "none";
        modalImg.style.display = "none";

        modalImg.onload = () => {
          if (modalLoading) modalLoading.style.display = "none";
          modalImg.style.display = "block";
        };
        modalImg.onerror = () => {
          if (modalLoading) modalLoading.style.display = "none";
          modalImg.style.display = "none";
          if (modalFallback) modalFallback.style.display = "flex";
        };
        modalImg.src = targetSrc;
      } else {
        if (modalLoading) modalLoading.style.display = "none";
        modalImg.style.display = "none";
        if (modalFallback) modalFallback.style.display = "flex";
      }

      modal.style.display = "flex";
      highlightCamera(activeCameraModalIndex);
    }

    function closePhotoModal() {
      const modal = document.getElementById("photo-modal");
      if (modal) modal.style.display = "none";
      resetCameraHighlights();
    }

    function closeModalOnBackdrop(event) {
      if (event.target.id === "photo-modal") {
        closePhotoModal();
      }
    }

    function prevPhoto() {
      if (!rawCameras || rawCameras.length === 0) return;
      let cur = parseInt(activeCameraModalIndex, 10);
      if (isNaN(cur)) cur = 0;
      let nextIdx = cur - 1;
      if (nextIdx < 0) nextIdx = rawCameras.length - 1;
      openPhotoModal(nextIdx);
    }

    function nextPhoto() {
      if (!rawCameras || rawCameras.length === 0) return;
      let cur = parseInt(activeCameraModalIndex, 10);
      if (isNaN(cur)) cur = 0;
      let nextIdx = cur + 1;
      if (nextIdx >= rawCameras.length) nextIdx = 0;
      openPhotoModal(nextIdx);
    }

    function lookThroughCurrentCamera() {
      const cam = rawCameras[activeCameraModalIndex];
      if (!cam) return;

      // Transform camera world pose by current model rotation
      const camPos = new THREE.Vector3(cam.pos[0], cam.pos[1], cam.pos[2]);
      const R = cam.R;
      const fwd = new THREE.Vector3(R[6], R[7], R[8]);
      const up = new THREE.Vector3(-R[3], -R[4], -R[5]);

      // Apply model rotation
      camPos.applyEuler(modelContainer.rotation);
      fwd.applyEuler(modelContainer.rotation);
      up.applyEuler(modelContainer.rotation);

      // Target ahead of the camera
      const targetPos = camPos.clone().addScaledVector(fwd, targetDim * 0.8);

      camera.position.copy(camPos);
      camera.up.copy(up);
      controls.target.copy(targetPos);
      controls.update();

      closePhotoModal();
      showToast(`Viewing from Camera: ${cam.label}`);
    }

    // Camera & Scene Toggle Controls
    function toggleCamerasVisibility() {
      areCamerasVisible = !areCamerasVisible;
      if (camerasGroup) camerasGroup.visible = areCamerasVisible;

      const tbCams = document.getElementById("tb-cams");
      if (tbCams) tbCams.classList.toggle("active", areCamerasVisible);
      const btnCams = document.getElementById("btn-toggle-cams");
      if (btnCams) btnCams.classList.toggle("active", areCamerasVisible);
      const btnCamsL = document.getElementById("btn-cams-toggle");
      if (btnCamsL) btnCamsL.classList.toggle("active", areCamerasVisible);

      showToast(areCamerasVisible ? "Cameras visible" : "Cameras hidden");
    }

    function toggleSightLines() {
      areSightLinesVisible = !areSightLinesVisible;
      if (sightLinesGroup) sightLinesGroup.visible = areSightLinesVisible;

      const tbSight = document.getElementById("tb-sight");
      if (tbSight) tbSight.classList.toggle("active", areSightLinesVisible);
      const btnSight = document.getElementById("btn-toggle-sight");
      if (btnSight) btnSight.classList.toggle("active", areSightLinesVisible);

      showToast(areSightLinesVisible ? "Sight lines visible" : "Sight lines hidden");
    }

    function toggleGizmo() {
      isGizmoVisible = !isGizmoVisible;
      if (trackballGizmo) trackballGizmo.visible = isGizmoVisible;

      const tbGizmo = document.getElementById("tb-gizmo");
      if (tbGizmo) tbGizmo.classList.toggle("active", isGizmoVisible);
      const btnGizmo = document.getElementById("btn-toggle-gizmo");
      if (btnGizmo) btnGizmo.classList.toggle("active", isGizmoVisible);

      showToast(isGizmoVisible ? "Trackball visible" : "Trackball hidden");
    }

    let isGridVisible = true;
    function toggleGrid() {
      isGridVisible = !isGridVisible;
      if (gridHelper) gridHelper.visible = isGridVisible;

      const tbGrid = document.getElementById("tb-grid");
      if (tbGrid) tbGrid.classList.toggle("active", isGridVisible);
      const btnGrid = document.getElementById("btn-toggle-grid");
      if (btnGrid) btnGrid.classList.toggle("active", isGridVisible);

      showToast(isGridVisible ? "Grid visible" : "Grid hidden");
    }

    function updateCameraScale(val) {
      cameraScaleMultiplier = parseFloat(val);
      const lbl = document.getElementById("lbl-cam-scale");
      if (lbl) lbl.textContent = cameraScaleMultiplier.toFixed(1) + "x";
      if (camerasGroup) {
        modelContainer.remove(camerasGroup);
        cameraColliders = [];
        cameraMeshes = [];
        buildCamerasFrustums();
      }
    }

    // UI Handlers for Top Toolbar & Engine Mode
    function handlePointsToolbarClick() {
      if (currentRenderMode !== "points") {
        setRenderMode("points");
        if (pointsMesh && !pointsMesh.visible) togglePointsVisibility();
      } else {
        togglePointsVisibility();
      }
    }

    function handleSplatsToolbarClick() {
      if (currentRenderMode !== "splats") {
        setRenderMode("splats");
        if (pointsMesh && !pointsMesh.visible) togglePointsVisibility();
      } else {
        togglePointsVisibility();
      }
    }

    function handleMeshToolbarClick() {
      if (!rawFaces || rawFaces.length === 0) {
        showToast("No mesh faces loaded. Import a .ply with faces.");
        triggerMeshImport();
        return;
      }
      if (currentRenderMode !== "mesh") {
        setRenderMode("mesh");
        if (meshGroup && !meshGroup.visible) toggleMeshVisibility();
      } else {
        toggleMeshVisibility();
      }
    }

    function toggleRenderMode() {
      const newMode = (currentRenderMode === "points") ? "mesh" : ((currentRenderMode === "mesh") ? "splats" : "points");
      setRenderMode(newMode);
    }

    function setRenderMode(mode) {
      currentRenderMode = mode;

      const tbMode = document.getElementById("tb-mode");
      if (tbMode) tbMode.classList.toggle("active", mode === "splats" && arePointsVisible);
      const tbPts = document.getElementById("tb-pts");
      if (tbPts) tbPts.classList.toggle("active", mode === "points" && arePointsVisible);
      const tbMesh = document.getElementById("tb-mesh");
      if (tbMesh) tbMesh.classList.toggle("active", mode === "mesh" && isMeshVisible);

      const grpMesh = document.getElementById("group-mesh-controls");
      const grpPoints = document.getElementById("group-point-controls");
      const grpShape = document.getElementById("group-shape");
      const lblT = document.getElementById("lbl-size-title");
      const sl = document.getElementById("slider-size");

      if (mode === "mesh") {
        if (!meshGroup && rawFaces && rawFaces.length > 0) {
          buildMeshGeometry();
        }
        if (meshGroup) meshGroup.visible = isMeshVisible;
        if (pointsMesh) pointsMesh.visible = false;
        if (grpMesh) grpMesh.style.display = "flex";
        if (grpPoints) grpPoints.style.display = "none";
        showToast("Mode: 3D Surface Mesh");
      } else if (mode === "points") {
        if (meshGroup) meshGroup.visible = false;
        if (pointsMesh) {
          pointsMesh.visible = arePointsVisible;
          if (material) {
            material.uniforms.uRenderMode.value = 0.0;
            material.depthWrite = true;
          }
        }
        if (grpMesh) grpMesh.style.display = "none";
        if (grpPoints) grpPoints.style.display = "flex";
        if (lblT) lblT.textContent = "Point Size";
        if (sl) { sl.min = 0.5; sl.max = 12.0; sl.step = 0.25; }
        if (grpShape) grpShape.style.display = "flex";
        if (sl) updatePointSize(sl.value);
        updateSplatBuffers();
        showToast("Mode: Points (Metashape)");
      } else if (mode === "splats") {
        if (meshGroup) meshGroup.visible = false;
        if (pointsMesh) {
          pointsMesh.visible = arePointsVisible;
          if (material) {
            material.uniforms.uRenderMode.value = 1.0;
            material.depthWrite = false;
          }
        }
        if (grpMesh) grpMesh.style.display = "none";
        if (grpPoints) grpPoints.style.display = "flex";
        if (lblT) lblT.textContent = "Splat Closure / Fill";
        if (sl) { sl.min = 0.5; sl.max = 10.0; sl.step = 0.25; }
        if (grpShape) grpShape.style.display = "none";
        if (sl) updatePointSize(sl.value);
        updateSplatBuffers();
        showToast("Mode: Gaussian Splats");
      }
      updateLayersUI();
    }

    function togglePointsVisibility() {
      if (!pointsMesh) return;
      arePointsVisible = !pointsMesh.visible;
      pointsMesh.visible = arePointsVisible;
      const btn = document.getElementById("btn-pts-toggle");
      if (btn) btn.classList.toggle("active", arePointsVisible);
      const tbPts = document.getElementById("tb-pts");
      if (tbPts && currentRenderMode === "points") tbPts.classList.toggle("active", arePointsVisible);
      const tbMode = document.getElementById("tb-mode");
      if (tbMode && currentRenderMode === "splats") tbMode.classList.toggle("active", arePointsVisible);
      showToast(arePointsVisible ? "Points: Visible" : "Points: Hidden");
    }

    function toggleMeshVisibility() {
      if (!meshGroup) {
        if (rawFaces && rawFaces.length > 0) {
          buildMeshGeometry();
        } else {
          return;
        }
      }
      isMeshVisible = !isMeshVisible;
      meshGroup.visible = isMeshVisible;
      const btn = document.getElementById("btn-mesh-toggle");
      if (btn) btn.classList.toggle("active", isMeshVisible);
      const tbMesh = document.getElementById("tb-mesh");
      if (tbMesh && currentRenderMode === "mesh") {
        tbMesh.classList.toggle("active", isMeshVisible);
      }
      showToast(isMeshVisible ? "Mesh: Visible" : "Mesh: Hidden");
    }

    function generateBoxUVs(geom) {
      if (!geom || !geom.attributes.position) return;
      if (geom.attributes.uv && geom.attributes.uv.count === geom.attributes.position.count) {
        return;
      }
      const pos = geom.attributes.position;
      const count = pos.count;
      const uvs = new Float32Array(count * 2);
      geom.computeBoundingBox();
      const bb = geom.boundingBox;
      const sizeX = Math.max(0.001, bb.max.x - bb.min.x);
      const sizeZ = Math.max(0.001, bb.max.z - bb.min.z);

      for (let i = 0; i < count; i++) {
        const x = pos.getX(i);
        const z = pos.getZ(i);
        uvs[i * 2] = (x - bb.min.x) / sizeX;
        uvs[i * 2 + 1] = (z - bb.min.z) / sizeZ;
      }
      geom.setAttribute("uv", new THREE.BufferAttribute(uvs, 2));
    }

    function createCurrentMeshSurfaceMaterial(hasValidColors) {
      if (meshShadingMode === "clay") {
        return new THREE.MeshStandardMaterial({
          color: 0xdbcfc4,
          roughness: 0.85,
          metalness: 0.05,
          flatShading: isMeshFlatShading,
          side: isMeshDoubleSided ? THREE.DoubleSide : THREE.FrontSide
        });
      }

      if (isMeshTextureActive && meshTextureObject) {
        return new THREE.MeshStandardMaterial({
          map: meshTextureObject,
          roughness: 0.60,
          metalness: 0.05,
          flatShading: isMeshFlatShading,
          side: isMeshDoubleSided ? THREE.DoubleSide : THREE.FrontSide
        });
      }

      return new THREE.MeshStandardMaterial({
        color: hasValidColors ? 0xffffff : 0xd1d5db,
        vertexColors: hasValidColors,
        roughness: 0.55,
        metalness: 0.10,
        flatShading: isMeshFlatShading,
        side: isMeshDoubleSided ? THREE.DoubleSide : THREE.FrontSide
      });
    }

    function buildMeshGeometry() {
      if (meshGroup) {
        modelContainer.remove(meshGroup);
        meshGroup = null;
        meshSurfaceObject = null;
        meshWireframeObject = null;
        meshGeometry = null;
      }

      if (!rawFaces || rawFaces.length === 0 || !rawPositions || rawPositions.length === 0) {
        return;
      }

      meshGroup = new THREE.Group();
      meshGroup.name = "meshGroup";

      meshGeometry = new THREE.BufferGeometry();
      meshGeometry.setAttribute("position", new THREE.BufferAttribute(rawPositions, 3));
      meshGeometry.setIndex(new THREE.BufferAttribute(rawFaces, 1));
      meshGeometry.computeVertexNormals();
      generateBoxUVs(meshGeometry);

      const hasValidColors = (hasModelColors && rawColorsRGB && rawColorsRGB.length === rawPositions.length);
      if (hasValidColors) {
        meshGeometry.setAttribute("color", new THREE.BufferAttribute(rawColorsRGB, 3));
      }

      const surfaceMat = createCurrentMeshSurfaceMaterial(hasValidColors);
      meshSurfaceObject = new THREE.Mesh(meshGeometry, surfaceMat);
      meshSurfaceObject.name = "meshSurface";
      meshGroup.add(meshSurfaceObject);

      const wireMat = new THREE.MeshBasicMaterial({
        color: 0x10b981,
        wireframe: true,
        transparent: true,
        opacity: 0.45
      });
      meshWireframeObject = new THREE.Mesh(meshGeometry, wireMat);
      meshWireframeObject.name = "meshWireframe";
      meshGroup.add(meshWireframeObject);

      modelContainer.add(meshGroup);
      applyMeshShadingState();

      meshGroup.visible = (currentRenderMode === "mesh" && isMeshVisible);
      if (pointsMesh && currentRenderMode === "mesh") {
        pointsMesh.visible = false;
      }
    }

    function applyMeshShadingState() {
      if (!meshSurfaceObject || !meshWireframeObject) return;

      const hasValidColors = (hasModelColors && rawColorsRGB && rawColorsRGB.length === rawPositions.length);
      meshSurfaceObject.material = createCurrentMeshSurfaceMaterial(hasValidColors);

      if (meshShadingMode === "solid") {
        meshSurfaceObject.visible = true;
        meshWireframeObject.visible = false;
      } else if (meshShadingMode === "clay") {
        meshSurfaceObject.visible = true;
        meshWireframeObject.visible = false;
      } else if (meshShadingMode === "wireframe") {
        meshSurfaceObject.visible = false;
        meshWireframeObject.visible = true;
      } else if (meshShadingMode === "solid_wireframe") {
        meshSurfaceObject.visible = true;
        meshWireframeObject.visible = true;
      }

      updateMeshUIState();
    }

    function setMeshShadingMode(mode) {
      meshShadingMode = mode;
      applyMeshShadingState();
      showToast("Mesh: " + mode.toUpperCase().replace(/_/g, " + "));
    }

    function toggleMeshFlatShading() {
      isMeshFlatShading = !isMeshFlatShading;
      applyMeshShadingState();
      const btn = document.getElementById("btn-mesh-flat-toggle");
      if (btn) {
        btn.textContent = isMeshFlatShading ? "&#128208; Faceted Flat" : "&#10024; Smooth Shading";
        btn.classList.toggle("active", isMeshFlatShading);
      }
      showToast(isMeshFlatShading ? "Shading: Faceted Flat" : "Shading: Smooth Normals");
    }

    function toggleMeshDoubleSided() {
      isMeshDoubleSided = !isMeshDoubleSided;
      applyMeshShadingState();
      const btn = document.getElementById("btn-mesh-doubleside-toggle");
      if (btn) {
        btn.textContent = isMeshDoubleSided ? "&#8596;&#65039; Double-Sided" : "&#8599;&#65039; Front-Only";
        btn.classList.toggle("active", isMeshDoubleSided);
      }
      showToast(isMeshDoubleSided ? "Faces: Two-Sided" : "Faces: Front Only");
    }

    function toggleMeshTexture() {
      if (!meshTextureObject) {
        triggerTextureImport();
        return;
      }
      isMeshTextureActive = !isMeshTextureActive;
      applyMeshShadingState();
      showToast(isMeshTextureActive ? "Texture: Enabled" : "Texture: Disabled");
    }

    function triggerTextureImport() {
      const inp = document.getElementById("input-texture-file");
      if (inp) {
        inp.value = "";
        inp.click();
      }
    }

    function handleTextureFileSelect(e) {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function(ev) {
        loadTextureFromDataUri(ev.target.result);
      };
      reader.readAsDataURL(file);
    }

    function loadTextureFromDataUri(uri) {
      new THREE.TextureLoader().load(uri, function(tex) {
        tex.wrapS = THREE.RepeatWrapping;
        tex.wrapT = THREE.RepeatWrapping;
        meshTextureObject = tex;
        isMeshTextureActive = true;
        if (meshShadingMode === "clay" || meshShadingMode === "wireframe") {
          meshShadingMode = "solid";
        }
        applyMeshShadingState();
        showToast("Texture successfully applied to 3D mesh!");
      });
    }

    function updateWireOpacity(val) {
      if (!meshWireframeObject) return;
      const op = parseInt(val, 10) / 100.0;
      meshWireframeObject.material.opacity = op;
      const lbl = document.getElementById("lbl-wire-opac");
      if (lbl) lbl.textContent = val + "%";
    }

    function triggerMeshImport() {
      const inp = document.getElementById("input-mesh-file");
      if (inp) {
        inp.value = "";
        inp.click();
      }
    }

    function handleMeshFileSelect(e) {
      const file = e.target.files[0];
      if (!file) return;
      handleModelFile(file);
    }

    function updateMeshUIState() {
      const btnSolid = document.getElementById("btn-mesh-solid");
      const btnClay = document.getElementById("btn-mesh-clay");
      const btnWire = document.getElementById("btn-mesh-wire");
      const btnSolidWire = document.getElementById("btn-mesh-solidwire");
      const wireCtrls = document.getElementById("mesh-wire-controls");
      const lblMode = document.getElementById("lbl-mesh-mode-name");
      const btnTex = document.getElementById("btn-mesh-texture-toggle");

      if (btnSolid) btnSolid.classList.toggle("active", meshShadingMode === "solid");
      if (btnClay) btnClay.classList.toggle("active", meshShadingMode === "clay");
      if (btnWire) btnWire.classList.toggle("active", meshShadingMode === "wireframe");
      if (btnSolidWire) btnSolidWire.classList.toggle("active", meshShadingMode === "solid_wireframe");
      if (wireCtrls) wireCtrls.style.display = (meshShadingMode === "wireframe" || meshShadingMode === "solid_wireframe") ? "flex" : "none";
      if (lblMode) lblMode.textContent = meshShadingMode.toUpperCase().replace(/_/g, " + ");
      if (btnTex) {
        btnTex.textContent = isMeshTextureActive && meshTextureObject ? "&#128444;&#65039; Texture: ON" : "&#128444;&#65039; Texture: OFF";
        btnTex.classList.toggle("active", isMeshTextureActive && !!meshTextureObject);
      }
    }

    function togglePointShape() {
      pointShapeValue = (pointShapeValue === 0) ? 1 : 0;
      if (material) material.uniforms.uPointShape.value = pointShapeValue;
      const btnShape = document.getElementById("btn-toggle-shape");
      if (btnShape) {
        btnShape.textContent = (pointShapeValue === 0) ? "&#9898; Smooth Circle (Antialiased)" : "&#11035; Crisp Square (Pixel)";
        btnShape.classList.toggle("active", pointShapeValue === 1);
      }
      showToast(pointShapeValue === 0 ? "Shape: Smooth Circle" : "Shape: Crisp Square");
    }

    function setPointShape(shape) {
      pointShapeValue = shape;
      if (material) material.uniforms.uPointShape.value = shape;
      const btnShape = document.getElementById("btn-toggle-shape");
      if (btnShape) {
        btnShape.textContent = (shape === 0) ? "&#9898; Smooth Circle (Antialiased)" : "&#11035; Crisp Square (Pixel)";
        btnShape.classList.toggle("active", shape === 1);
      }
    }

    function setColorMode(mode) {
      currentColorMode = mode;
      document.querySelectorAll(".segmented-control-5 .seg-btn, .segmented-control-3 .seg-btn").forEach(b => b.classList.remove("active"));
      if (mode === "rgb") {
        const b = document.getElementById("btn-col-rgb"); if (b) b.classList.add("active");
        if (meshGroup) setMeshShadingMode("solid");
      }
      if (mode === "clay") {
        const b = document.getElementById("btn-col-clay"); if (b) b.classList.add("active");
        if (meshGroup) setMeshShadingMode("clay");
      }
      if (mode === "height") { const b = document.getElementById("btn-col-h"); if (b) b.classList.add("active"); }
      if (mode === "density") { const b = document.getElementById("btn-col-dens"); if (b) b.classList.add("active"); }
      if (mode === "semantic") { const b = document.getElementById("btn-col-sem"); if (b) b.classList.add("active"); }
      updateSplatBuffers();
      showToast("Display Mode: " + mode.toUpperCase());
    }

    function updatePointSize(val) {
      pointSizeValue = parseFloat(val);
      const unit = currentRenderMode === "points" ? " px" : "x (Fill)";
      document.getElementById("lbl-size").textContent = pointSizeValue.toFixed(1) + unit;
      material.uniforms.uPointSizePx.value = pointSizeValue;
    }

    function updateOpacityMultiplier(val) {
      globalOpacityMultiplier = parseFloat(val);
      document.getElementById("lbl-opac").textContent = globalOpacityMultiplier.toFixed(2) + "x";
      updateSplatBuffers();
    }

    function updateSphereOpacity(val) {
      sphereOpacityValue = parseFloat(val);
      if (sphereMesh) {
        sphereMesh.material.opacity = sphereOpacityValue;
      }
      const pctElem = document.getElementById("lbl-sphere-opac");
      if (pctElem) pctElem.textContent = Math.round(sphereOpacityValue * 100) + "%";
    }

    function updateSphereScale(val) {
      sphereScaleMultiplier = parseFloat(val);
      const scaleElem = document.getElementById("lbl-sphere-scale");
      if (scaleElem) scaleElem.textContent = sphereScaleMultiplier.toFixed(1) + "x";
    }

    // 3-Axis Model Alignment & Leveling (World vs Local Alignment)
    function toScreenCoords(pos3D, cam) {
      const v = pos3D.clone().project(cam);
      return {
        x: (v.x * 0.5 + 0.5) * window.innerWidth,
        y: (-(v.y * 0.5) + 0.5) * window.innerHeight
      };
    }

    function setAlignmentSpace(space) {
      alignmentSpace = space;
      const bWorld = document.getElementById("btn-space-world");
      const bLocal = document.getElementById("btn-space-local");
      if (bWorld) bWorld.classList.toggle("active", space === "world");
      if (bLocal) bLocal.classList.toggle("active", space === "local");

      if (trackballGizmo) {
        if (space === "local" && modelContainer) {
          trackballGizmo.quaternion.copy(modelContainer.quaternion);
        } else {
          trackballGizmo.quaternion.set(0, 0, 0, 1);
        }
      }
      showToast(space === "world" ? "Model Alignment: World Grid (Leveling)" : "Model Alignment: Local Object (Intrinsic)");
    }

    function updateModelRotation() {
      const rxDeg = parseFloat(document.getElementById("rot-x-slider").value) || 0;
      const ryDeg = parseFloat(document.getElementById("rot-y-slider").value) || 0;
      const rzDeg = parseFloat(document.getElementById("rot-z-slider").value) || 0;

      document.getElementById("rot-x-val").textContent = Math.round(rxDeg) + "&deg;";
      document.getElementById("rot-y-val").textContent = Math.round(ryDeg) + "&deg;";
      document.getElementById("rot-z-val").textContent = Math.round(rzDeg) + "&deg;";

      const rx = rxDeg * Math.PI / 180.0;
      const ry = ryDeg * Math.PI / 180.0;
      const rz = rzDeg * Math.PI / 180.0;

      if (modelContainer) {
        modelContainer.rotation.set(rx, ry, rz, "XYZ");
      }
      if (trackballGizmo) {
        if (alignmentSpace === "local" && modelContainer) {
          trackballGizmo.quaternion.copy(modelContainer.quaternion);
        } else {
          trackballGizmo.quaternion.set(0, 0, 0, 1);
        }
      }
    }

    function syncSlidersFromModel() {
      if (!modelContainer) return;
      const rxDeg = Math.round(modelContainer.rotation.x * 180 / Math.PI);
      const ryDeg = Math.round(modelContainer.rotation.y * 180 / Math.PI);
      const rzDeg = Math.round(modelContainer.rotation.z * 180 / Math.PI);

      const sx = document.getElementById("rot-x-slider");
      const sy = document.getElementById("rot-y-slider");
      const sz = document.getElementById("rot-z-slider");
      const vx = document.getElementById("rot-x-val");
      const vy = document.getElementById("rot-y-val");
      const vz = document.getElementById("rot-z-val");

      if (sx) sx.value = rxDeg;
      if (sy) sy.value = ryDeg;
      if (sz) sz.value = rzDeg;
      if (vx) vx.textContent = rxDeg + "&deg;";
      if (vy) vy.textContent = ryDeg + "&deg;";
      if (vz) vz.textContent = rzDeg + "&deg;";
    }

    function stepModelRotation(axis, deg) {
      if (!modelContainer) return;
      const rad = deg * Math.PI / 180.0;
      let rotAx = new THREE.Vector3();
      if (alignmentSpace === "local") {
        const localA = (axis === "x") ? new THREE.Vector3(1, 0, 0) :
                       (axis === "y") ? new THREE.Vector3(0, 1, 0) :
                                        new THREE.Vector3(0, 0, 1);
        rotAx.copy(localA).applyQuaternion(modelContainer.quaternion).normalize();
      } else {
        if (axis === "x") rotAx.set(1, 0, 0);
        else if (axis === "y") rotAx.set(0, 1, 0);
        else rotAx.set(0, 0, 1);
      }
      const qDelta = new THREE.Quaternion().setFromAxisAngle(rotAx, rad);
      if (alignmentSpace === "local") {
        modelContainer.quaternion.multiply(qDelta);
      } else {
        modelContainer.quaternion.premultiply(qDelta);
      }
      modelContainer.rotation.setFromQuaternion(modelContainer.quaternion, "XYZ");
      syncSlidersFromModel();
      if (alignmentSpace === "local" && trackballGizmo) {
        trackballGizmo.quaternion.copy(modelContainer.quaternion);
      } else if (trackballGizmo) {
        trackballGizmo.quaternion.set(0, 0, 0, 1);
      }
    }

    function resetModelRotation() {
      document.getElementById("rot-x-slider").value = 0;
      document.getElementById("rot-y-slider").value = 0;
      document.getElementById("rot-z-slider").value = 0;
      updateModelRotation();
      showToast("Model rotation reset to 0&deg;");
    }

    function copyRotationCode() {
      const rx = document.getElementById("rot-x-slider").value;
      const ry = document.getElementById("rot-y-slider").value;
      const rz = document.getElementById("rot-z-slider").value;
      const code = `rotation = c(${rx}, ${ry}, ${rz})`;
      navigator.clipboard.writeText(code).then(() => {
        showToast("Copied: " + code);
      }).catch(() => {
        prompt("Copy rotation parameters:", code);
      });
    }

    function showToast(text) {
      const t = document.getElementById("toast-msg");
      if (!t) return;
      t.textContent = text;
      t.style.display = "block";
      setTimeout(() => { t.style.display = "none"; }, 2500);
    }

    let currentTheme = "dark";
    function toggleBackgroundTheme() {
      currentTheme = (currentTheme === "dark") ? "light" : "dark";
      setBackgroundTheme(currentTheme);
    }

    function setBackgroundTheme(theme) {
      currentTheme = theme;
      const btnTheme = document.getElementById("btn-toggle-theme");
      const tbTheme = document.getElementById("tb-theme");
      if (btnTheme) {
        btnTheme.classList.toggle("active", theme === "light");
        btnTheme.textContent = (theme === "light") ? "&#127769; Switch to Dark Theme" : "&#9728;&#65039; Metashape Light Theme";
      }
      if (tbTheme) tbTheme.classList.toggle("active", theme === "light");

      if (theme === "dark") {
        scene.background.set("#090d16");
        if (gridHelper) gridHelper.material.color.setHex(0x10b981);
        showToast("Dark theme enabled");
      } else {
        scene.background.set("#f1f5f9");
        if (gridHelper) gridHelper.material.color.setHex(0x94a3b8);
        showToast("Metashape light theme enabled");
      }
    }

    function setCameraPreset(preset) {
      if (preset === "front") {
        camera.position.set(0, -targetDim * 2.2, 0);
      } else if (preset === "side") {
        camera.position.set(targetDim * 2.2, 0, 0);
      } else if (preset === "top") {
        camera.position.set(0, 0, targetDim * 2.5);
      } else if (preset === "reset") {
        camera.position.set(0, -targetDim * 2.2, targetDim * 1.1);
      }
      controls.target.set(0, 0, 0);
      controls.update();
    }

    function captureSnapshot() {
      const renderCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
      renderer.render(scene, renderCam);
      const dataUrl = renderer.domElement.toDataURL("image/png");
      const a = document.createElement("a");
      a.href = dataUrl;
      a.download = "3d_model_snapshot.png";
      a.click();
    }

    function onWindowResize() {
      const aspect = window.innerWidth / window.innerHeight;
      camera.aspect = aspect;
      camera.updateProjectionMatrix();

      if (orthoCamera) {
        const dist = camera.position.distanceTo(controls.target);
        const fovRad = (camera.fov * Math.PI) / 180.0;
        const halfH = Math.tan(fovRad / 2.0) * dist;
        const halfW = halfH * aspect;
        orthoCamera.left = -halfW;
        orthoCamera.right = halfW;
        orthoCamera.top = halfH;
        orthoCamera.bottom = -halfH;
        orthoCamera.updateProjectionMatrix();
      }

      renderer.setSize(window.innerWidth, window.innerHeight);
    }

    // -------------------------------------------------------------------------
    // Camera Projection Toggle (Perspective vs Orthographic CAD)
    // -------------------------------------------------------------------------
    function initOrthoCamera() {
      const aspect = window.innerWidth / window.innerHeight;
      const d = targetDim * 1.5;
      orthoCamera = new THREE.OrthographicCamera(-d * aspect, d * aspect, d, -d, 0.001, 1000);
      orthoCamera.position.set(0, -targetDim * 2.2, targetDim * 1.1);
      orthoCamera.up.set(0, 0, 1);
    }

    function toggleCameraProjection() {
      if (!orthoCamera) initOrthoCamera();
      isOrtho = !isOrtho;
      const aspect = window.innerWidth / window.innerHeight;
      const tbProj = document.getElementById("tb-proj");

      if (isOrtho) {
        const dist = camera.position.distanceTo(controls.target);
        const fovRad = (camera.fov * Math.PI) / 180.0;
        const halfH = Math.tan(fovRad / 2.0) * dist;
        const halfW = halfH * aspect;

        orthoCamera.left = -halfW;
        orthoCamera.right = halfW;
        orthoCamera.top = halfH;
        orthoCamera.bottom = -halfH;
        orthoCamera.position.copy(camera.position);
        orthoCamera.quaternion.copy(camera.quaternion);
        orthoCamera.updateProjectionMatrix();

        controls.object = orthoCamera;
        activeCamera = orthoCamera;
        if (tbProj) {
          tbProj.classList.add("active");
          tbProj.innerHTML = "&#128208; Ortho";
        }
        showToast("Orthographic CAD projection enabled");
      } else {
        camera.position.copy(orthoCamera.position);
        camera.quaternion.copy(orthoCamera.quaternion);
        camera.updateProjectionMatrix();

        controls.object = camera;
        activeCamera = camera;
        if (tbProj) {
          tbProj.classList.remove("active");
          tbProj.innerHTML = "&#128208; Persp";
        }
        showToast("Perspective projection restored");
      }
      controls.update();
    }

    // -------------------------------------------------------------------------
    // 3D Measurement Ruler (Tape Measure)
    // -------------------------------------------------------------------------
    function toggleRulerTool() {
      isRulerActive = !isRulerActive;
      const tb = document.getElementById("tb-ruler");
      if (tb) {
        tb.classList.toggle("active", isRulerActive);
        tb.classList.toggle("tb-ruler", isRulerActive);
      }
      const groupM = document.getElementById("group-measurements");
      if (groupM) groupM.style.display = (isRulerActive || rulerMeasurements.length > 0) ? "flex" : "none";

      if (isRulerActive) {
        showToast("3D Ruler active: Click 1st point on point cloud");
        const instr = document.getElementById("ruler-instruction");
        if (instr) instr.textContent = "1. Click 1st point on point cloud...";
        if (!rulerHoverMarker) {
          const sGeo = new THREE.SphereGeometry(Math.max(0.005, targetDim * 0.015), 16, 16);
          const sMat = new THREE.MeshBasicMaterial({ color: 0xfacc15, wireframe: true, depthTest: false });
          rulerHoverMarker = new THREE.Mesh(sGeo, sMat);
          rulerGroup.add(rulerHoverMarker);
        }
        rulerHoverMarker.visible = true;
      } else {
        rulerPointA = null;
        if (rulerTempLine) { rulerGroup.remove(rulerTempLine); rulerTempLine = null; }
        if (rulerHoverMarker) rulerHoverMarker.visible = false;
        showToast("3D Ruler disabled");
      }
    }

    function handleRulerClick(point) {
      if (!rulerPointA) {
        rulerPointA = point.clone();
        const sGeo = new THREE.SphereGeometry(Math.max(0.006, targetDim * 0.016), 16, 16);
        const sMat = new THREE.MeshBasicMaterial({ color: 0x22c55e, depthTest: false });
        const markerA = new THREE.Mesh(sGeo, sMat);
        markerA.position.copy(rulerPointA);
        rulerGroup.add(markerA);

        const instr = document.getElementById("ruler-instruction");
        if (instr) instr.textContent = "2. Now click 2nd point to measure distance...";
        showToast("Point 1 placed. Click 2nd point to measure!");
      } else {
        const rulerPointB = point.clone();
        const sGeo = new THREE.SphereGeometry(Math.max(0.006, targetDim * 0.016), 16, 16);
        const sMat = new THREE.MeshBasicMaterial({ color: 0xef4444, depthTest: false });
        const markerB = new THREE.Mesh(sGeo, sMat);
        markerB.position.copy(rulerPointB);
        rulerGroup.add(markerB);

        const distM = rulerPointA.distanceTo(rulerPointB);
        const distStr = distM < 1.0 ? (distM * 100).toFixed(2) + " cm" : distM.toFixed(3) + " m";
        const dx = Math.abs(rulerPointB.x - rulerPointA.x);
        const dy = Math.abs(rulerPointB.y - rulerPointA.y);
        const dz = Math.abs(rulerPointB.z - rulerPointA.z);

        const dxStr = dx < 1.0 ? (dx * 100).toFixed(1) + "cm" : dx.toFixed(2) + "m";
        const dyStr = dy < 1.0 ? (dy * 100).toFixed(1) + "cm" : dy.toFixed(2) + "m";
        const dzStr = dz < 1.0 ? (dz * 100).toFixed(1) + "cm" : dz.toFixed(2) + "m";

        const lineGeo = new THREE.BufferGeometry().setFromPoints([rulerPointA, rulerPointB]);
        const lineMat = new THREE.LineDashedMaterial({ color: 0xfacc15, dashSize: 0.03, gapSize: 0.015, linewidth: 2, depthTest: false });
        const line = new THREE.Line(lineGeo, lineMat);
        line.computeLineDistances();
        rulerGroup.add(line);

        const midPoint = rulerPointA.clone().add(rulerPointB).multiplyScalar(0.5);
        const sprite = createMeasurementBadge(distStr, midPoint);
        rulerGroup.add(sprite);

        if (rulerTempLine) { rulerGroup.remove(rulerTempLine); rulerTempLine = null; }

        const mObj = {
          id: rulerMeasurements.length + 1,
          pA: rulerPointA,
          pB: rulerPointB,
          distStr: distStr,
          dxStr: dxStr, dyStr: dyStr, dzStr: dzStr,
          markerA: markerA,
          markerB: markerB,
          line: line,
          sprite: sprite
        };
        rulerMeasurements.push(mObj);
        updateMeasurementsListUI();

        showToast(`Distance measured: ${distStr}`);
        rulerPointA = null;
        const instr = document.getElementById("ruler-instruction");
        if (instr) instr.textContent = "1. Click 1st point for next measurement...";
      }
    }

    function createMeasurementBadge(text, position) {
      const canvas = document.createElement("canvas");
      canvas.width = 256;
      canvas.height = 70;
      const ctx = canvas.getContext("2d");

      ctx.fillStyle = "rgba(15, 23, 42, 0.92)";
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(4, 4, 248, 62, 10);
      else ctx.rect(4, 4, 248, 62);
      ctx.fill();
      ctx.lineWidth = 3;
      ctx.strokeStyle = "#facc15";
      ctx.stroke();

      ctx.fillStyle = "#facc15";
      ctx.font = "bold 26px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, 128, 35);

      const texture = new THREE.CanvasTexture(canvas);
      const spriteMat = new THREE.SpriteMaterial({ map: texture, depthTest: false });
      const sprite = new THREE.Sprite(spriteMat);
      sprite.position.copy(position);
      const badgeScale = Math.max(0.04, targetDim * 0.12);
      sprite.scale.set(badgeScale * 2.2, badgeScale * 0.6, 1);
      return sprite;
    }

    function updateMeasurementsListUI() {
      const list = document.getElementById("measurements-list");
      if (!list) return;
      if (rulerMeasurements.length === 0) {
        list.innerHTML = `<div style="font-size: 11px; color: #facc15;" id="ruler-instruction">1. Click 1st point on point cloud...</div>`;
        return;
      }
      let html = `<div style="font-size: 11px; color: #facc15;" id="ruler-instruction">1. Click 1st point for next measurement...</div>`;
      rulerMeasurements.forEach((m, idx) => {
        html += `
          <div style="display: flex; justify-content: space-between; align-items: center; background: rgba(255,255,255,0.06); padding: 5px 8px; border-radius: 6px; font-size: 11px;">
            <div>
              <span style="color: #facc15; font-weight: bold;">#${m.id}: ${m.distStr}</span>
              <span style="font-size: 10px; color: var(--text-muted); margin-left: 6px;">(&Delta;X: ${m.dxStr}, &Delta;Y: ${m.dyStr}, &Delta;Z: ${m.dzStr})</span>
            </div>
            <button class="btn-reset-sm" onclick="removeMeasurement(${idx})" style="color: #ef4444;" title="Delete measurement">&#10005;</button>
          </div>
        `;
      });
      list.innerHTML = html;
    }

    function removeMeasurement(idx) {
      const m = rulerMeasurements[idx];
      if (m) {
        if (m.markerA) rulerGroup.remove(m.markerA);
        if (m.markerB) rulerGroup.remove(m.markerB);
        if (m.line) rulerGroup.remove(m.line);
        if (m.sprite) rulerGroup.remove(m.sprite);
        rulerMeasurements.splice(idx, 1);
        updateMeasurementsListUI();
      }
    }

    function clearMeasurements() {
      rulerMeasurements.forEach(m => {
        if (m.markerA) rulerGroup.remove(m.markerA);
        if (m.markerB) rulerGroup.remove(m.markerB);
        if (m.line) rulerGroup.remove(m.line);
        if (m.sprite) rulerGroup.remove(m.sprite);
      });
      rulerMeasurements = [];
      rulerPointA = null;
      if (rulerTempLine) { rulerGroup.remove(rulerTempLine); rulerTempLine = null; }
      updateMeasurementsListUI();
      showToast("All measurements cleared.");
    }

    // -------------------------------------------------------------------------
    // RealityCapture Export Tools (Client-Side Direct Downloads)
    // -------------------------------------------------------------------------
    function openExportModal() {
      const modal = document.getElementById("export-modal");
      if (modal) modal.style.display = "flex";
    }

    function closeExportModal() {
      const modal = document.getElementById("export-modal");
      if (modal) modal.style.display = "none";
    }

    function closeExportModalOnBackdrop(e) {
      if (e.target.id === "export-modal") closeExportModal();
    }

    function exportFilteredPLY() {
      if (numGaussians === 0) {
        showToast("No model available to export.");
        return;
      }

      const alphaArr = alphasAttr.array;
      const sizeArr = sizesAttr.array;
      const colArr = colorsAttr.array;

      const visibleIndices = [];
      for (let i = 0; i < numGaussians; i++) {
        if (alphaArr[i] > 0.04 && sizeArr[i] > 0.0001) {
          visibleIndices.push(i);
        }
      }

      const count = visibleIndices.length;
      if (count === 0) {
        showToast("No visible points within active filters.");
        return;
      }

      const useWorld = document.getElementById("chk-export-world") ? document.getElementById("chk-export-world").checked : true;
      const ox = useWorld ? centerOrigin[0] : 0;
      const oy = useWorld ? centerOrigin[1] : 0;
      const oz = useWorld ? centerOrigin[2] : 0;

      const headerLines = [
        "ply",
        "format binary_little_endian 1.0",
        "comment Exported from pliman RealityCapture-grade 3D Viewer",
        "element vertex " + count,
        "property float x",
        "property float y",
        "property float z",
        "property uchar red",
        "property uchar green",
        "property uchar blue",
        "property uchar alpha",
        "end_header",
        ""
      ];
      const headerStr = headerLines.join(String.fromCharCode(10));

      const headerBytes = new TextEncoder().encode(headerStr);
      const vertexBytes = count * 16;
      const buffer = new ArrayBuffer(vertexBytes);
      const dataView = new DataView(buffer);

      let off = 0;
      for (let k = 0; k < count; k++) {
        const i = visibleIndices[k];
        dataView.setFloat32(off, rawPositions[i * 3] + ox, true);
        dataView.setFloat32(off + 4, rawPositions[i * 3 + 1] + oy, true);
        dataView.setFloat32(off + 8, rawPositions[i * 3 + 2] + oz, true);

        dataView.setUint8(off + 12, Math.min(255, Math.max(0, Math.round(colArr[i * 3] * 255))));
        dataView.setUint8(off + 13, Math.min(255, Math.max(0, Math.round(colArr[i * 3 + 1] * 255))));
        dataView.setUint8(off + 14, Math.min(255, Math.max(0, Math.round(colArr[i * 3 + 2] * 255))));
        dataView.setUint8(off + 15, Math.min(255, Math.max(0, Math.round(alphaArr[i] * 255))));
        off += 16;
      }

      const blob = new Blob([headerBytes, buffer], { type: "application/octet-stream" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `pliman_filtered_${count}_points.ply`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(url), 10000);
      showToast(`Exported: ${count.toLocaleString()} points in binary PLY!`);
      closeExportModal();
    }

    function exportFilteredTXT() {
      if (numGaussians === 0) {
        showToast("No model available to export.");
        return;
      }
      const alphaArr = alphasAttr.array;
      const sizeArr = sizesAttr.array;
      const colArr = colorsAttr.array;

      const useWorld = document.getElementById("chk-export-world") ? document.getElementById("chk-export-world").checked : true;
      const ox = useWorld ? centerOrigin[0] : 0;
      const oy = useWorld ? centerOrigin[1] : 0;
      const oz = useWorld ? centerOrigin[2] : 0;

      const lines = ["# X Y Z R G B (Exported from pliman RealityCapture 3D Viewer)"];
      let count = 0;
      for (let i = 0; i < numGaussians; i++) {
        if (alphaArr[i] > 0.04 && sizeArr[i] > 0.0001) {
          const x = (rawPositions[i * 3] + ox).toFixed(5);
          const y = (rawPositions[i * 3 + 1] + oy).toFixed(5);
          const z = (rawPositions[i * 3 + 2] + oz).toFixed(5);
          const r = Math.round(colArr[i * 3] * 255);
          const g = Math.round(colArr[i * 3 + 1] * 255);
          const b = Math.round(colArr[i * 3 + 2] * 255);
          lines.push(`${x} ${y} ${z} ${r} ${g} ${b}`);
          count++;
        }
      }

      if (count === 0) {
        showToast("No visible points within active filters.");
        return;
      }

      const blob = new Blob([lines.join(String.fromCharCode(10))], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `pliman_points_${count}.txt`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(url), 10000);
      showToast(`Exported: ${count.toLocaleString()} points in TXT!`);
      closeExportModal();
    }

    function exportCamerasTXT() {
      if (!rawCameras || rawCameras.length === 0) {
        showToast("No cameras available to export.");
        return;
      }
      const lines = [
        "# COLMAP calibrated cameras export",
        "# IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME"
      ];

      rawCameras.forEach((cam, i) => {
        const id = cam.index || (i + 1);
        const name = cam.label || `image_${id}.jpg`;
        const R = cam.R;
        const tr = R[0] + R[4] + R[8];
        let qw, qx, qy, qz;
        if (tr > 0) {
          const S = Math.sqrt(tr + 1.0) * 2;
          qw = 0.25 * S;
          qx = (R[7] - R[5]) / S;
          qy = (R[2] - R[6]) / S;
          qz = (R[3] - R[1]) / S;
        } else if ((R[0] > R[4]) && (R[0] > R[8])) {
          const S = Math.sqrt(1.0 + R[0] - R[4] - R[8]) * 2;
          qw = (R[7] - R[5]) / S;
          qx = 0.25 * S;
          qy = (R[1] + R[3]) / S;
          qz = (R[2] + R[6]) / S;
        } else if (R[4] > R[8]) {
          const S = Math.sqrt(1.0 + R[4] - R[0] - R[8]) * 2;
          qw = (R[2] - R[6]) / S;
          qx = (R[1] + R[3]) / S;
          qy = 0.25 * S;
          qz = (R[5] + R[7]) / S;
        } else {
          const S = Math.sqrt(1.0 + R[8] - R[0] - R[4]) * 2;
          qw = (R[3] - R[1]) / S;
          qx = (R[2] + R[6]) / S;
          qy = (R[5] + R[7]) / S;
          qz = 0.25 * S;
        }

        const C = cam.worldPos || [cam.pos[0] + centerOrigin[0], cam.pos[1] + centerOrigin[1], cam.pos[2] + centerOrigin[2]];
        const tx = -(R[0]*C[0] + R[1]*C[1] + R[2]*C[2]);
        const ty = -(R[3]*C[0] + R[4]*C[1] + R[5]*C[2]);
        const tz = -(R[6]*C[0] + R[7]*C[1] + R[8]*C[2]);

        lines.push(`${id} ${qw.toFixed(7)} ${qx.toFixed(7)} ${qy.toFixed(7)} ${qz.toFixed(7)} ${tx.toFixed(5)} ${ty.toFixed(5)} ${tz.toFixed(5)} 1 ${name}`);
        lines.push("");
      });

      const blob = new Blob([lines.join(String.fromCharCode(10))], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `colmap_exported_cameras_${rawCameras.length}.txt`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(url), 10000);
      showToast(`Exported: ${rawCameras.length} cameras in COLMAP format!`);
      closeExportModal();
    }

    // FPS Counter & Render Loop
    let lastTime = performance.now();
    let frameCount = 0;
    const fpsElem = document.getElementById("fps-counter");

    function animate(time) {
      requestAnimationFrame(animate);
      controls.update();

      // Trackball Gizmo smoothly tracks orbit pivot target and maintains constant elegant screen size
      if (trackballGizmo && trackballGizmo.visible) {
        trackballGizmo.position.copy(controls.target);
        const camDist = (isOrtho && orthoCamera ? orthoCamera : camera).position.distanceTo(controls.target);
        trackballGizmo.scale.setScalar(Math.max(0.02, camDist * 0.16 * sphereScaleMultiplier));
        if (alignmentSpace === "local" && modelContainer) {
          trackballGizmo.quaternion.copy(modelContainer.quaternion);
        } else {
          trackballGizmo.quaternion.set(0, 0, 0, 1);
        }
      }

      const renderCam = (isOrtho && orthoCamera) ? orthoCamera : camera;
      renderer.render(scene, renderCam);

      frameCount++;
      if (time - lastTime >= 1000) {
        const fps = Math.round((frameCount * 1000) / (time - lastTime));
        fpsElem.textContent = fps + " FPS";
        frameCount = 0;
        lastTime = time;
      }
    }

    // -------------------------------------------------------------------------
    // Client-Side Point Cloud & Gaussian Splat Parsers (PLY, SPLAT, TXT)
    // 100% Offline / Local / Zero Server Dependencies via HTML5 File API
    // -------------------------------------------------------------------------

    const localPhotoStore = {};
    let arePointsVisible = true;

    function openWelcomeModal() {
      const modal = document.getElementById("dropzone-welcome");
      if (modal) modal.style.display = "flex";
      updateLayersUI();
    }

    function closeWelcomeModal() {
      const modal = document.getElementById("dropzone-welcome");
      if (modal) modal.style.display = "none";
    }

    function togglePointsVisibility() {
      if (pointsMesh) {
        pointsMesh.visible = !pointsMesh.visible;
        arePointsVisible = pointsMesh.visible;
        const btn = document.getElementById("btn-pts-toggle");
        if (btn) btn.classList.toggle("active", pointsMesh.visible);
        const tbPts = document.getElementById("tb-pts");
        if (tbPts) tbPts.classList.toggle("active", currentRenderMode === "points" && pointsMesh.visible);
        const tbSplats = document.getElementById("tb-mode");
        if (tbSplats) tbSplats.classList.toggle("active", currentRenderMode === "splats" && pointsMesh.visible);
        showToast(pointsMesh.visible ? "3D Model visible" : "3D Model hidden");
      } else {
        showToast("No 3D model loaded yet.");
      }
    }

    function toggleCamerasVisibility() {
      if (camerasGroup) {
        camerasGroup.visible = !camerasGroup.visible;
        areCamerasVisible = camerasGroup.visible;
        const btn = document.getElementById("btn-cams-toggle");
        if (btn) btn.classList.toggle("active", camerasGroup.visible);
        const tbCams = document.getElementById("tb-cams");
        if (tbCams) tbCams.classList.toggle("active", camerasGroup.visible);
        showToast(camerasGroup.visible ? "Cameras visible" : "Cameras hidden");
      } else {
        showToast("No cameras loaded yet.");
      }
    }

    function triggerModelImport() {
      const inp = document.getElementById("input-model-file");
      if (inp) inp.click();
    }

    function triggerCamsImport() {
      const inp = document.getElementById("input-cams-file");
      if (inp) inp.click();
    }

    function triggerPhotosImport() {
      const inp = document.getElementById("input-photos-file");
      if (inp) inp.click();
    }

    function triggerPhotosDirImport() {
      const inp = document.getElementById("input-photos-dir");
      if (inp) inp.click();
    }

    function triggerBatchImport() {
      const inp = document.getElementById("input-batch-files");
      if (inp) inp.click();
    }

    function handleBatchFilesSelect(event) {
      const files = event.target.files;
      if (files && files.length > 0) {
        handleIncomingFiles(files);
      }
      event.target.value = "";
    }

    function handleModelFileSelect(event) {
      const files = event.target.files;
      if (!files || files.length === 0) return;
      if (files.length > 1) {
        handleIncomingFiles(files);
      } else {
        handleModelFile(files[0]);
      }
      event.target.value = "";
    }

    function handleCamsFileSelect(event) {
      const files = event.target.files;
      if (!files || files.length === 0) return;
      if (files.length > 1) {
        handleIncomingFiles(files);
      } else {
        handleCamerasFile(files[0]);
      }
      event.target.value = "";
    }

    function handlePhotosSelect(event) {
      const files = event.target.files;
      if (!files || files.length === 0) return;
      attachLocalPhotos(files);
      event.target.value = "";
    }

    function updateLayersUI() {
      // Model Status
      const welcomeMStatus = document.getElementById("welcome-model-status");
      const welcomeMCard = document.getElementById("welcome-model-card");
      const hudMStatus = document.getElementById("hud-model-status");
      const btnPtsToggle = document.getElementById("btn-pts-toggle");

      if (numGaussians > 0) {
        const txt = numGaussians.toLocaleString() + " pts";
        if (welcomeMStatus) {
          welcomeMStatus.textContent = "&#10003; " + txt;
          welcomeMStatus.className = "field-status loaded";
        }
        if (welcomeMCard) welcomeMCard.classList.add("loaded");
        if (hudMStatus) hudMStatus.textContent = txt;
        if (btnPtsToggle) btnPtsToggle.classList.toggle("active", arePointsVisible && currentRenderMode !== "mesh");
      } else {
        if (welcomeMStatus) {
          welcomeMStatus.textContent = "No file selected";
          welcomeMStatus.className = "field-status empty";
        }
        if (welcomeMCard) welcomeMCard.classList.remove("loaded");
        if (hudMStatus) hudMStatus.textContent = "None";
        if (btnPtsToggle) btnPtsToggle.classList.remove("active");
      }

      // Mesh Status
      const hudMeshStatus = document.getElementById("hud-mesh-status");
      const btnMeshToggle = document.getElementById("btn-mesh-toggle");
      if (rawFaces && rawFaces.length > 0) {
        const fCount = Math.floor(rawFaces.length / 3);
        const txt = fCount.toLocaleString() + " faces";
        if (hudMeshStatus) hudMeshStatus.textContent = txt;
        if (btnMeshToggle) btnMeshToggle.classList.toggle("active", isMeshVisible);
      } else {
        if (hudMeshStatus) hudMeshStatus.textContent = "None";
        if (btnMeshToggle) btnMeshToggle.classList.remove("active");
      }

      // Cameras Status
      const welcomeCStatus = document.getElementById("welcome-cams-status");
      const welcomeCCard = document.getElementById("welcome-cams-card");
      const hudCStatus = document.getElementById("hud-cams-status");
      const btnCamsToggle = document.getElementById("btn-cams-toggle");
      const rowCamScale = document.getElementById("row-cam-scale");

      if (rawCameras && rawCameras.length > 0) {
        const txt = rawCameras.length + " cameras";
        if (welcomeCStatus) {
          welcomeCStatus.textContent = "&#10003; " + txt;
          welcomeCStatus.className = "field-status loaded";
        }
        if (welcomeCCard) welcomeCCard.classList.add("loaded");
        if (hudCStatus) hudCStatus.textContent = txt;
        if (btnCamsToggle) btnCamsToggle.classList.add("active");
        if (rowCamScale) rowCamScale.style.display = "flex";
      } else {
        if (welcomeCStatus) {
          welcomeCStatus.textContent = "No cameras selected";
          welcomeCStatus.className = "field-status empty";
        }
        if (welcomeCCard) welcomeCCard.classList.remove("loaded");
        if (hudCStatus) hudCStatus.textContent = "None";
        if (btnCamsToggle) btnCamsToggle.classList.remove("active");
        if (rowCamScale) rowCamScale.style.display = "none";
      }

      // Photos Status
      const welcomePStatus = document.getElementById("welcome-photos-status");
      const welcomePCard = document.getElementById("welcome-photos-card");
      const hudPStatus = document.getElementById("hud-photos-status");

      const photoCount = Object.keys(localPhotoStore).length;
      if (photoCount > 0) {
        const txt = photoCount + " linked photos";
        if (welcomePStatus) {
          welcomePStatus.textContent = "&#10003; " + txt;
          welcomePStatus.className = "field-status loaded";
        }
        if (welcomePCard) welcomePCard.classList.add("loaded");
        if (hudPStatus) hudPStatus.textContent = txt;
      } else {
        if (welcomePStatus) {
          welcomePStatus.textContent = "No photos linked";
          welcomePStatus.className = "field-status empty";
        }
        if (welcomePCard) welcomePCard.classList.remove("loaded");
        if (hudPStatus) hudPStatus.textContent = "None";
      }
    }

    function frameCameras() {
      if (!rawCameras || rawCameras.length === 0) return;
      let cx = 0, cy = 0, cz = 0;
      rawCameras.forEach(c => {
        cx += c.pos[0]; cy += c.pos[1]; cz += c.pos[2];
      });
      cx /= rawCameras.length; cy /= rawCameras.length; cz /= rawCameras.length;
      let maxDist = 0.5;
      rawCameras.forEach(c => {
        const d = Math.sqrt(Math.pow(c.pos[0] - cx, 2) + Math.pow(c.pos[1] - cy, 2) + Math.pow(c.pos[2] - cz, 2));
        if (d > maxDist) maxDist = d;
      });
      targetDim = Math.max(0.5, maxDist * 2.2);
      camera.position.set(cx, cy - targetDim * 1.8, cz + targetDim * 1.0);
      controls.target.set(cx, cy, cz);
      controls.update();
      if (gridHelper) {
        scene.remove(gridHelper);
        gridHelper = new THREE.GridHelper(targetDim * 3.0, 24, 0x10b981, 0x334155);
        gridHelper.rotation.x = Math.PI / 2;
        gridHelper.position.set(cx, cy, cz - targetDim * 0.4);
        scene.add(gridHelper);
      }
    }

    function handleModelFile(file, onComplete) {
      const name = file.name.toLowerCase();
      showToast(`Loading ${file.name}...`);
      const reader = new FileReader();

      if (name.endsWith(".splat")) {
        reader.onload = (e) => {
          try {
            const parsed = parseSPLAT(e.target.result);
            loadModelData(parsed, file.name);
            if (onComplete) onComplete();
          } catch (err) {
            alert("Error parsing SPLAT file: " + err.message);
          }
        };
        reader.readAsArrayBuffer(file);
      } else if (name.endsWith(".ply")) {
        reader.onload = (e) => {
          try {
            const parsed = parsePLY(e.target.result);
            loadModelData(parsed, file.name);
            if (onComplete) onComplete();
          } catch (err) {
            alert("Error parsing PLY file: " + err.message);
          }
        };
        reader.readAsArrayBuffer(file);
      } else if (name.endsWith(".txt")) {
        reader.onload = (e) => {
          try {
            const parsed = parseTXT(e.target.result);
            loadModelData(parsed, file.name);
            if (onComplete) onComplete();
          } catch (err) {
            alert("Error parsing TXT file: " + err.message);
          }
        };
        reader.readAsText(file);
      } else {
        alert("Unsupported file format: " + file.name + ". Please select a .PLY, .SPLAT, or .TXT file.");
      }
    }

    function handleCamerasFile(file, onComplete) {
      const name = file.name.toLowerCase();
      const reader = new FileReader();
      reader.onload = (e) => {
        try {
          const text = e.target.result;
          let parsedCams = [];
          if (name.endsWith(".xml")) {
            parsedCams = parseMetashapeXml(text);
          } else {
            parsedCams = parseColmapImagesTxt(text);
          }
          if (parsedCams.length === 0) {
            showToast("No valid cameras found in " + file.name);
            return;
          }
          rawCameras = parsedCams;
          // Apply centering if a model is already loaded
          rawCameras.forEach(cam => {
            if (cam.worldPos) {
              cam.pos = [
                cam.worldPos[0] - centerOrigin[0],
                cam.worldPos[1] - centerOrigin[1],
                cam.worldPos[2] - centerOrigin[2]
              ];
            }
          });
          buildCamerasFrustums();
          if (numGaussians === 0) {
            frameCameras();
          }
          updateLayersUI();
          showToast(`Loaded ${rawCameras.length} camera poses from ${file.name}!`);
          if (onComplete) onComplete();
        } catch (err) {
          alert("Error parsing cameras: " + err.message);
        }
      };
      reader.readAsText(file);
    }

    function attachLocalPhotos(fileList) {
      let count = 0;
      for (let i = 0; i < fileList.length; i++) {
        const file = fileList[i];
        const low = file.name.toLowerCase();
        localPhotoStore[low] = URL.createObjectURL(file);
        count++;
      }
      let matched = 0;
      if (rawCameras && rawCameras.length > 0) {
        rawCameras.forEach(cam => {
          const lbl = (cam.label || "").toLowerCase();
          const base = (lbl.indexOf(String.fromCharCode(92)) !== -1 ? lbl.split(String.fromCharCode(92)) : lbl.split("/")).pop();
          if (localPhotoStore[lbl]) {
            cam.img = localPhotoStore[lbl];
            matched++;
          } else if (localPhotoStore[base]) {
            cam.img = localPhotoStore[base];
            matched++;
          }
        });
        showToast(`Linked ${matched} photos to camera frustums!`);
        const modal = document.getElementById("photo-modal");
        if (modal && modal.style.display === "flex") {
          openPhotoModal(activeCameraModalIndex);
        }
      } else {
        showToast(`Stored ${count} local photos in memory.`);
      }
      updateLayersUI();
    }

    // -------------------------------------------------------------------------
    // High-Performance Binary & ASCII PLY Parser
    // -------------------------------------------------------------------------
    function parsePLY(arrayBuffer) {
      const bytes = new Uint8Array(arrayBuffer);
      let headerText = "";
      let headerEndIndex = 0;
      for (let i = 0; i < Math.min(bytes.length, 65536); i++) {
        headerText += String.fromCharCode(bytes[i]);
        if (headerText.indexOf("end_header") !== -1) {
          const nextNl = headerText.indexOf(String.fromCharCode(10), headerText.indexOf("end_header"));
          if (nextNl !== -1) {
            headerEndIndex = nextNl + 1;
            break;
          }
        }
      }

      if (headerEndIndex === 0) {
        throw new Error("Invalid PLY file: end_header marker not found.");
      }

      const lines = headerText.split(String.fromCharCode(10)).map(l => l.replace(String.fromCharCode(13), ""));
      let isBinary = false;
      let isBigEndian = false;
      let vertexCount = 0;
      let faceCount = 0;
      let faceCountType = "uchar";
      let faceIndexType = "int";
      const properties = [];
      let inVertex = false;
      let inFace = false;

      for (let line of lines) {
        line = line.trim();
        if (line.startsWith("format binary_little_endian")) isBinary = true;
        else if (line.startsWith("format binary_big_endian")) { isBinary = true; isBigEndian = true; }
        else if (line.startsWith("element vertex")) {
          vertexCount = parseInt(line.split(/[ \t]+/)[2], 10);
          inVertex = true;
          inFace = false;
        } else if (line.startsWith("element face")) {
          faceCount = parseInt(line.split(/[ \t]+/)[2], 10);
          inVertex = false;
          inFace = true;
        } else if (line.startsWith("element ")) {
          inVertex = false;
          inFace = false;
        } else if (inVertex && line.startsWith("property")) {
          const parts = line.split(/[ \t]+/);
          properties.push({ type: parts[1], name: parts[2] });
        } else if (inFace && line.startsWith("property list")) {
          const parts = line.split(/[ \t]+/);
          faceCountType = parts[2];
          faceIndexType = parts[3];
        }
      }

      if (vertexCount <= 0) {
        throw new Error("No vertex elements found in PLY header.");
      }

      const positions = new Float32Array(vertexCount * 3);
      const colors = new Float32Array(vertexCount * 3);
      const opacities = new Float32Array(vertexCount);
      const scales = new Float32Array(vertexCount);
      const semantics = new Int32Array(vertexCount);

      const xIdx = properties.findIndex(p => p.name === "x");
      const yIdx = properties.findIndex(p => p.name === "y");
      const zIdx = properties.findIndex(p => p.name === "z");
      const rIdx = properties.findIndex(p => ["red", "diffuse_red", "r", "f_dc_0"].includes(p.name));
      const gIdx = properties.findIndex(p => ["green", "diffuse_green", "g", "f_dc_1"].includes(p.name));
      const bIdx = properties.findIndex(p => ["blue", "diffuse_blue", "b", "f_dc_2"].includes(p.name));
      const opacIdx = properties.findIndex(p => ["opacity", "alpha"].includes(p.name));
      const s0Idx = properties.findIndex(p => p.name === "scale_0");
      const s1Idx = properties.findIndex(p => p.name === "scale_1");
      const s2Idx = properties.findIndex(p => p.name === "scale_2");
      const semIdx = properties.findIndex(p => ["semantic_id", "class", "label"].includes(p.name));

      const isGaussianSplat = (s0Idx !== -1 && s1Idx !== -1 && s2Idx !== -1);
      const isSHColor = properties.some(p => p.name === "f_dc_0");
      const shC0 = 0.28209479177387814;
      let facesArray = new Uint32Array(0);

      if (isBinary) {
        const dataView = new DataView(arrayBuffer, headerEndIndex);
        const typeSizes = {
          "char": 1, "int8": 1, "uchar": 1, "uint8": 1,
          "short": 2, "int16": 2, "ushort": 2, "uint16": 2,
          "int": 4, "int32": 4, "uint": 4, "uint32": 4,
          "float": 4, "float32": 4, "double": 8, "float64": 8
        };
        let vertexByteSize = 0;
        const propOffsets = [];
        for (let p of properties) {
          propOffsets.push(vertexByteSize);
          vertexByteSize += (typeSizes[p.type] || 4);
        }

        function createGetter(type) {
          switch (type) {
            case "float": case "float32": return (o) => dataView.getFloat32(o, !isBigEndian);
            case "double": case "float64": return (o) => dataView.getFloat64(o, !isBigEndian);
            case "uchar": case "uint8": return (o) => dataView.getUint8(o);
            case "char": case "int8": return (o) => dataView.getInt8(o);
            case "ushort": case "uint16": return (o) => dataView.getUint16(o, !isBigEndian);
            case "short": case "int16": return (o) => dataView.getInt16(o, !isBigEndian);
            case "uint": case "uint32": return (o) => dataView.getUint32(o, !isBigEndian);
            case "int": case "int32": return (o) => dataView.getInt32(o, !isBigEndian);
            default: return (o) => dataView.getFloat32(o, !isBigEndian);
          }
        }

        const getX = xIdx >= 0 ? createGetter(properties[xIdx].type) : () => 0;
        const getY = yIdx >= 0 ? createGetter(properties[yIdx].type) : () => 0;
        const getZ = zIdx >= 0 ? createGetter(properties[zIdx].type) : () => 0;
        const getR = rIdx >= 0 ? createGetter(properties[rIdx].type) : null;
        const getG = gIdx >= 0 ? createGetter(properties[gIdx].type) : null;
        const getB = bIdx >= 0 ? createGetter(properties[bIdx].type) : null;
        const getOpac = opacIdx >= 0 ? createGetter(properties[opacIdx].type) : null;
        const getS0 = s0Idx >= 0 ? createGetter(properties[s0Idx].type) : null;
        const getS1 = s1Idx >= 0 ? createGetter(properties[s1Idx].type) : null;
        const getS2 = s2Idx >= 0 ? createGetter(properties[s2Idx].type) : null;
        const getSem = semIdx >= 0 ? createGetter(properties[semIdx].type) : null;
        const rIsUchar = rIdx >= 0 && ["uchar", "uint8"].includes(properties[rIdx].type);

        for (let i = 0; i < vertexCount; i++) {
          const vOff = i * vertexByteSize;
          positions[i * 3]     = getX(vOff + propOffsets[xIdx]);
          positions[i * 3 + 1] = getY(vOff + propOffsets[yIdx]);
          positions[i * 3 + 2] = getZ(vOff + propOffsets[zIdx]);

          if (getR && getG && getB) {
            let r = getR(vOff + propOffsets[rIdx]);
            let g = getG(vOff + propOffsets[gIdx]);
            let b = getB(vOff + propOffsets[bIdx]);
            if (isSHColor) {
              r = Math.min(1.0, Math.max(0.0, 0.5 + shC0 * r));
              g = Math.min(1.0, Math.max(0.0, 0.5 + shC0 * g));
              b = Math.min(1.0, Math.max(0.0, 0.5 + shC0 * b));
            } else if (rIsUchar || r > 1.0 || g > 1.0 || b > 1.0) {
              r /= 255.0; g /= 255.0; b /= 255.0;
            }
            colors[i * 3]     = r;
            colors[i * 3 + 1] = g;
            colors[i * 3 + 2] = b;
          } else {
            colors[i * 3] = (faceCount > 0) ? 0.82 : 0.22;
            colors[i * 3 + 1] = (faceCount > 0) ? 0.80 : 0.72;
            colors[i * 3 + 2] = (faceCount > 0) ? 0.78 : 0.28;
          }

          if (getOpac) {
            let op = getOpac(vOff + propOffsets[opacIdx]);
            if (isGaussianSplat) {
              op = 1.0 / (1.0 + Math.exp(-op));
            } else if (op > 1.0) {
              op /= 255.0;
            }
            opacities[i] = op;
          } else {
            opacities[i] = 1.0;
          }

          if (isGaussianSplat) {
            const s0 = getS0(vOff + propOffsets[s0Idx]);
            const s1 = getS1(vOff + propOffsets[s1Idx]);
            const s2 = getS2(vOff + propOffsets[s2Idx]);
            scales[i] = Math.exp((s0 + s1 + s2) / 3.0);
          } else {
            scales[i] = 0.01;
          }

          semantics[i] = getSem ? getSem(vOff + propOffsets[semIdx]) : 0;
        }

        // Read binary faces if present
        if (faceCount > 0) {
          const facesList = [];
          let fOffset = headerEndIndex + vertexCount * vertexByteSize;
          for (let f = 0; f < faceCount; f++) {
            if (fOffset >= arrayBuffer.byteLength) break;
            const nv = (faceCountType === "uchar" || faceCountType === "uint8") ? dataView.getUint8(fOffset) : dataView.getInt32(fOffset, !isBigEndian);
            fOffset += (faceCountType === "uchar" || faceCountType === "uint8") ? 1 : 4;
            const vIdx = [];
            for (let k = 0; k < nv; k++) {
              if (faceIndexType === "short" || faceIndexType === "uint16" || faceIndexType === "int16") {
                vIdx.push(dataView.getUint16(fOffset, !isBigEndian));
                fOffset += 2;
              } else {
                vIdx.push(dataView.getInt32(fOffset, !isBigEndian));
                fOffset += 4;
              }
            }
            if (nv === 3) {
              facesList.push(vIdx[0], vIdx[1], vIdx[2]);
            } else if (nv === 4) {
              facesList.push(vIdx[0], vIdx[1], vIdx[2], vIdx[0], vIdx[2], vIdx[3]);
            } else if (nv > 4) {
              for (let k = 1; k < nv - 1; k++) {
                facesList.push(vIdx[0], vIdx[k], vIdx[k + 1]);
              }
            }
          }
          facesArray = new Uint32Array(facesList);
        }
      } else {
        // ASCII format
        const textData = new TextDecoder().decode(new Uint8Array(arrayBuffer, headerEndIndex));
        const dataLines = textData.trim().split(String.fromCharCode(10)).map(l => l.replace(String.fromCharCode(13), ""));
        const limit = Math.min(vertexCount, dataLines.length);
        for (let i = 0; i < limit; i++) {
          const parts = dataLines[i].trim().split(/[ \t]+/);
          if (parts.length < 3) continue;
          positions[i * 3]     = parseFloat(parts[xIdx]);
          positions[i * 3 + 1] = parseFloat(parts[yIdx]);
          positions[i * 3 + 2] = parseFloat(parts[zIdx]);
          if (rIdx >= 0 && gIdx >= 0 && bIdx >= 0 && parts.length > Math.max(rIdx, gIdx, bIdx)) {
            let r = parseFloat(parts[rIdx]);
            let g = parseFloat(parts[gIdx]);
            let b = parseFloat(parts[bIdx]);
            if (r > 1.0 || g > 1.0 || b > 1.0) { r /= 255.0; g /= 255.0; b /= 255.0; }
            colors[i * 3] = r; colors[i * 3 + 1] = g; colors[i * 3 + 2] = b;
          } else {
            colors[i * 3] = (faceCount > 0) ? 0.82 : 0.22;
            colors[i * 3 + 1] = (faceCount > 0) ? 0.80 : 0.72;
            colors[i * 3 + 2] = (faceCount > 0) ? 0.78 : 0.28;
          }
          opacities[i] = 1.0;
          scales[i] = 0.01;
          semantics[i] = 0;
        }

        if (faceCount > 0) {
          const facesList = [];
          for (let f = limit; f < limit + faceCount && f < dataLines.length; f++) {
            const parts = dataLines[f].trim().split(/[ \t]+/);
            if (parts.length >= 4) {
              const nv = parseInt(parts[0], 10);
              const vIdx = [];
              for (let k = 1; k <= nv && k < parts.length; k++) {
                vIdx.push(parseInt(parts[k], 10));
              }
              if (nv === 3) {
                facesList.push(vIdx[0], vIdx[1], vIdx[2]);
              } else if (nv === 4) {
                facesList.push(vIdx[0], vIdx[1], vIdx[2], vIdx[0], vIdx[2], vIdx[3]);
              }
            }
          }
          facesArray = new Uint32Array(facesList);
        }
      }

      return {
        positions,
        colors,
        opacities,
        scales,
        semantics,
        faces: facesArray,
        count: vertexCount,
        faceCount: Math.floor(facesArray.length / 3),
        hasColors: (rIdx >= 0 && gIdx >= 0 && bIdx >= 0)
      };
    }

    // -------------------------------------------------------------------------
    // Binary .SPLAT Parser (32-byte layout)
    // -------------------------------------------------------------------------
    function parseSPLAT(arrayBuffer) {
      const count = Math.floor(arrayBuffer.byteLength / 32);
      if (count <= 0) throw new Error("Empty or corrupted .splat file.");
      const dataView = new DataView(arrayBuffer);
      const positions = new Float32Array(count * 3);
      const colors = new Float32Array(count * 3);
      const opacities = new Float32Array(count);
      const scales = new Float32Array(count);
      const semantics = new Int32Array(count);

      for (let i = 0; i < count; i++) {
        const off = i * 32;
        positions[i * 3]     = dataView.getFloat32(off, true);
        positions[i * 3 + 1] = dataView.getFloat32(off + 4, true);
        positions[i * 3 + 2] = dataView.getFloat32(off + 8, true);

        const sx = dataView.getFloat32(off + 12, true);
        const sy = dataView.getFloat32(off + 16, true);
        const sz = dataView.getFloat32(off + 20, true);
        scales[i] = Math.cbrt(Math.max(1e-6, sx * sy * sz));

        colors[i * 3]     = dataView.getUint8(off + 24) / 255.0;
        colors[i * 3 + 1] = dataView.getUint8(off + 25) / 255.0;
        colors[i * 3 + 2] = dataView.getUint8(off + 26) / 255.0;
        opacities[i]      = dataView.getUint8(off + 27) / 255.0;
        semantics[i]      = 0;
      }
      return { positions, colors, opacities, scales, semantics, count };
    }

    // -------------------------------------------------------------------------
    // TXT Point Cloud Parser (COLMAP points3D.txt or XYZ RGB)
    // -------------------------------------------------------------------------
    function parseTXT(text) {
      const lines = text.trim().split(String.fromCharCode(10)).map(l => l.replace(String.fromCharCode(13), ""));
      const posList = [];
      const colList = [];

      for (let line of lines) {
        line = line.trim();
        if (!line || line.startsWith("#")) continue;
        const p = line.split(/[ \t]+/);
        if (p.length >= 7 && !isNaN(parseFloat(p[1])) && !isNaN(parseFloat(p[2]))) {
          // COLMAP points3D.txt: POINT3D_ID, X, Y, Z, R, G, B, ERROR, ...
          posList.push(parseFloat(p[1]), parseFloat(p[2]), parseFloat(p[3]));
          colList.push(parseFloat(p[4]) / 255.0, parseFloat(p[5]) / 255.0, parseFloat(p[6]) / 255.0);
        } else if (p.length >= 3 && !isNaN(parseFloat(p[0]))) {
          // XYZ [R G B]
          posList.push(parseFloat(p[0]), parseFloat(p[1]), parseFloat(p[2]));
          if (p.length >= 6 && !isNaN(parseFloat(p[3]))) {
            let r = parseFloat(p[3]), g = parseFloat(p[4]), b = parseFloat(p[5]);
            if (r > 1.0 || g > 1.0 || b > 1.0) { r /= 255.0; g /= 255.0; b /= 255.0; }
            colList.push(r, g, b);
          } else {
            colList.push(0.8, 0.8, 0.8);
          }
        }
      }

      const count = Math.floor(posList.length / 3);
      if (count === 0) throw new Error("No valid 3D points found in text file.");
      const positions = new Float32Array(posList);
      const colors = new Float32Array(colList);
      const opacities = new Float32Array(count).fill(1.0);
      const scales = new Float32Array(count).fill(0.01);
      const semantics = new Int32Array(count).fill(0);
      return { positions, colors, opacities, scales, semantics, count };
    }

    // -------------------------------------------------------------------------
    // COLMAP images.txt Camera Parser
    // -------------------------------------------------------------------------
    function parseColmapImagesTxt(text) {
      const lines = text.trim().split(String.fromCharCode(10)).map(l => l.replace(String.fromCharCode(13), ""));
      const cams = [];
      for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line || line.startsWith("#")) continue;
        const parts = line.split(/[ \t]+/);
        if (parts.length >= 10) {
          const imgId = parseInt(parts[0]);
          const qw = parseFloat(parts[1]);
          const qx = parseFloat(parts[2]);
          const qy = parseFloat(parts[3]);
          const qz = parseFloat(parts[4]);
          const tx = parseFloat(parts[5]);
          const ty = parseFloat(parts[6]);
          const tz = parseFloat(parts[7]);
          const camId = parseInt(parts[8]);
          const name = parts[9];

          // Quaternion to rotation matrix (world to camera)
          const q0 = qw, q1 = qx, q2 = qy, q3 = qz;
          const r00 = 1 - 2*(q2*q2 + q3*q3), r01 = 2*(q1*q2 - q0*q3),     r02 = 2*(q1*q3 + q0*q2);
          const r10 = 2*(q1*q2 + q0*q3),     r11 = 1 - 2*(q1*q1 + q3*q3), r12 = 2*(q2*q3 - q0*q1);
          const r20 = 2*(q1*q3 - q0*q2),     r21 = 2*(q2*q3 + q0*q1),     r22 = 1 - 2*(q1*q1 + q2*q2);

          // Camera center C = -R^T * t
          const cx = -(r00*tx + r10*ty + r20*tz);
          const cy = -(r01*tx + r11*ty + r21*tz);
          const cz = -(r02*tx + r12*ty + r22*tz);

          const lowName = name.toLowerCase();
          const baseName = (lowName.indexOf(String.fromCharCode(92)) !== -1 ? lowName.split(String.fromCharCode(92)) : lowName.split("/")).pop();
          const cachedImg = localPhotoStore[lowName] || localPhotoStore[baseName] || null;

          cams.push({
            index: cams.length + 1,
            label: name,
            worldPos: [cx, cy, cz],
            pos: [cx, cy, cz],
            R: [r00, r01, r02, r10, r11, r12, r20, r21, r22],
            fx: 1000, fy: 1000, w: 1920, h: 1080,
            img: cachedImg
          });

          // Skip 2D points line
          i++;
        }
      }
      return cams;
    }

    // -------------------------------------------------------------------------
    // Agisoft Metashape cameras.xml Parser
    // -------------------------------------------------------------------------
    function parseMetashapeXml(xmlText) {
      const parser = new DOMParser();
      const xmlDoc = parser.parseFromString(xmlText, "text/xml");
      const camNodes = xmlDoc.getElementsByTagName("camera");
      const cams = [];
      for (let i = 0; i < camNodes.length; i++) {
        const camNode = camNodes[i];
        const label = camNode.getAttribute("label");
        const transformNode = camNode.getElementsByTagName("transform")[0];
        if (!transformNode) continue;
        const vals = transformNode.textContent.trim().split(/[ \t]+/).map(Number);
        if (vals.length === 16) {
          const cx = vals[3], cy = vals[7], cz = vals[11];
          const r00 = vals[0], r01 = vals[4], r02 = vals[8];
          const r10 = vals[1], r11 = vals[5], r12 = vals[9];
          const r20 = vals[2], r21 = vals[6], r22 = vals[10];

          const lowName = (label || "").toLowerCase();
          const baseName = (lowName.indexOf(String.fromCharCode(92)) !== -1 ? lowName.split(String.fromCharCode(92)) : lowName.split("/")).pop();
          const cachedImg = localPhotoStore[lowName] || localPhotoStore[baseName] || null;

          cams.push({
            index: cams.length + 1,
            label: label,
            worldPos: [cx, cy, cz],
            pos: [cx, cy, cz],
            R: [r00, r01, r02, r10, r11, r12, r20, r21, r22],
            fx: 1000, fy: 1000, w: 1920, h: 1080,
            img: cachedImg
          });
        }
      }
      return cams;
    }

    // -------------------------------------------------------------------------
    // Load Model Data into WebGL Buffers & Update HUD
    // -------------------------------------------------------------------------
    function loadModelData(parsed, fileName) {
      if (!parsed || parsed.count === 0) {
        showToast("Error: No 3D points found in file.");
        return;
      }

      // Sample points for robust centroid and quantiles
      const step = Math.max(1, Math.floor(parsed.count / 5000));
      const sampleX = [], sampleY = [], sampleZ = [];
      for (let i = 0; i < parsed.count; i += step) {
        sampleX.push(parsed.positions[i * 3]);
        sampleY.push(parsed.positions[i * 3 + 1]);
        sampleZ.push(parsed.positions[i * 3 + 2]);
      }
      sampleX.sort((a, b) => a - b);
      sampleY.sort((a, b) => a - b);
      sampleZ.sort((a, b) => a - b);

      const mid = Math.floor(sampleX.length / 2);
      const cx = sampleX[mid] || 0;
      const cy = sampleY[mid] || 0;
      const cz = sampleZ[mid] || 0;

      const q02 = Math.floor(sampleX.length * 0.02);
      const q98 = Math.floor(sampleX.length * 0.98);
      const dx = Math.max(0.05, (sampleX[q98] || 1) - (sampleX[q02] || 0));
      const dy = Math.max(0.05, (sampleY[q98] || 1) - (sampleY[q02] || 0));
      const dz = Math.max(0.05, (sampleZ[q98] || 1) - (sampleZ[q02] || 0));
      const globalDim = Math.max(dx, dy, dz, 0.1);

      const dists = [];
      for (let i = 0; i < sampleX.length; i++) {
        const d = Math.sqrt(
          Math.pow(sampleX[i] - cx, 2) +
          Math.pow(sampleY[i] - cy, 2) +
          Math.pow(sampleZ[i] - cz, 2)
        );
        dists.push(d);
      }
      dists.sort((a, b) => a - b);
      const q70Dist = Math.max(0.1, dists[Math.floor(dists.length * 0.70)] || 0.5);
      const newTargetDim = Math.min(globalDim, Math.max(0.2, q70Dist * 2.2));

      // Subtract centroid from positions
      for (let i = 0; i < parsed.count; i++) {
        parsed.positions[i * 3]     -= cx;
        parsed.positions[i * 3 + 1] -= cy;
        parsed.positions[i * 3 + 2] -= cz;
      }

      // Update globals
      numGaussians = parsed.count;
      rawPositions = parsed.positions;
      rawColorsRGB = parsed.colors;
      rawOpacities = parsed.opacities;
      rawScales = parsed.scales;
      rawSemantics = parsed.semantics;
      rawFaces = parsed.faces || new Uint32Array(0);
      numFaces = parsed.faceCount || Math.floor(rawFaces.length / 3);
      hasModelFaces = rawFaces.length > 0;
      hasModelColors = parsed.hasColors !== false;
      targetDim = newTargetDim;
      centerOrigin = [cx, cy, cz];

      // Compute spatial densities and distance distributions for real-time online noise filtering
      computeSpatialDensities();
      computeDistanceDistribution();

      // Remove old point cloud mesh
      if (pointsMesh) {
        modelContainer.remove(pointsMesh);
        if (geometry) geometry.dispose();
      }

      // Rebuild geometry and buffers
      buildSplatGeometry();
      buildMeshGeometry();
      if (hasModelFaces) {
        setRenderMode("mesh");
      } else {
        setColorMode(currentColorMode);
      }

      // Re-center cameras if loaded
      if (rawCameras && rawCameras.length > 0) {
        rawCameras.forEach(cam => {
          if (cam.worldPos) {
            cam.pos = [cam.worldPos[0] - cx, cam.worldPos[1] - cy, cam.worldPos[2] - cz];
          }
        });
        buildCamerasFrustums();
      }

      // Update Grid & Camera Framing
      if (gridHelper) scene.remove(gridHelper);
      gridHelper = new THREE.GridHelper(newTargetDim * 4.0, 24, 0x10b981, 0x334155);
      gridHelper.rotation.x = Math.PI / 2;
      gridHelper.position.set(0, 0, -newTargetDim * 0.45);
      scene.add(gridHelper);

      camera.position.set(0, -newTargetDim * 2.2, newTargetDim * 1.1);
      controls.target.set(0, 0, 0);
      controls.update();

      // Update HUD Metrics
      const elTotal = document.getElementById("lbl-total-pts");
      if (elTotal) elTotal.textContent = parsed.count.toLocaleString();
      const elVis = document.getElementById("lbl-visible-splats");
      if (elVis) elVis.textContent = parsed.count.toLocaleString();
      const elScale = document.getElementById("lbl-scale-extent");
      if (elScale) elScale.textContent = newTargetDim.toFixed(2) + " m";
      const elBbox = document.getElementById("lbl-bbox-extent");
      if (elBbox) elBbox.textContent = `${dx.toFixed(1)} x ${dy.toFixed(1)} x ${dz.toFixed(1)} m`;
      const elOxy = document.getElementById("lbl-origin-xy");
      if (elOxy) elOxy.textContent = `${cx.toFixed(2)}, ${cy.toFixed(2)}`;
      const elOz = document.getElementById("lbl-origin-z");
      if (elOz) elOz.textContent = cz.toFixed(2) + " m";
      const elCoords = document.getElementById("coords-text");
      if (elCoords) elCoords.textContent = `X: ${cx.toFixed(3)} | Y: ${cy.toFixed(3)} | Z: ${cz.toFixed(3)}`;

      // Update layers UI
      updateLayersUI();

      // Hide welcome dropzone
      const welcomeOverlay = document.getElementById("dropzone-welcome");
      if (welcomeOverlay) welcomeOverlay.style.display = "none";

      if (hasModelFaces) {
        showToast(`Loaded 3D Mesh (${parsed.count.toLocaleString()} vertices, ${numFaces.toLocaleString()} faces)!`);
      } else {
        showToast(`Loaded ${parsed.count.toLocaleString()} points (${fileName || "Model"})!`);
      }
    }

    // -------------------------------------------------------------------------
    // Drag & Drop Listener Setup
    // -------------------------------------------------------------------------
    function setupDragAndDrop() {
      const dragOverlay = document.getElementById("dropzone-drag");
      let dragCounter = 0;

      window.addEventListener("dragenter", (e) => {
        e.preventDefault();
        dragCounter++;
        if (dragOverlay) dragOverlay.style.display = "flex";
      });

      window.addEventListener("dragover", (e) => {
        e.preventDefault();
      });

      window.addEventListener("dragleave", (e) => {
        e.preventDefault();
        dragCounter--;
        if (dragCounter <= 0 && dragOverlay) {
          dragOverlay.style.display = "none";
          dragCounter = 0;
        }
      });

      window.addEventListener("drop", (e) => {
        e.preventDefault();
        dragCounter = 0;
        if (dragOverlay) dragOverlay.style.display = "none";
        const files = e.dataTransfer.files;
        if (files && files.length > 0) {
          handleIncomingFiles(files);
        }
      });
    }

    function handleIncomingFiles(files) {
      const imgFiles = [];
      let modelFile = null;
      let camFile = null;

      for (let i = 0; i < files.length; i++) {
        const f = files[i];
        const low = f.name.toLowerCase();
        if (low.endsWith(".ply") || low.endsWith(".splat")) {
          modelFile = f;
        } else if (low.endsWith(".txt")) {
          if (low.includes("camera") || low.includes("images") || low.includes("pose")) {
            camFile = f;
          } else {
            modelFile = f;
          }
        } else if (low.endsWith(".xml")) {
          camFile = f;
        } else if ([".jpg", ".jpeg", ".png", ".webp", ".bmp", ".tif", ".tiff"].some(ext => low.endsWith(ext))) {
          imgFiles.push(f);
        }
      }

      if (imgFiles.length > 0) {
        attachLocalPhotos(imgFiles);
      }
      if (modelFile && camFile) {
        handleModelFile(modelFile, () => {
          handleCamerasFile(camFile);
        });
      } else {
        if (camFile) handleCamerasFile(camFile);
        if (modelFile) handleModelFile(modelFile);
      }
    }

    window.onload = init;
  </script>
</body>
</html>
')

  writeLines(html_content, file)
  cli::cli_alert_success("Generated standalone 3D viewer at {.path {file}}.")

  # Launch browser
  if (isTRUE(browser)) {
    norm_path <- normalizePath(file, winslash = "/", mustWork = FALSE)
    file_url <- paste0("file:///", norm_path)

    if (!is.null(browser_path) && nzchar(browser_path)) {
      cli::cli_alert_info("Opening in custom browser: {.val {browser_path}}...")
      system2(browser_path, args = file_url, wait = FALSE)
    } else {
      cli::cli_alert_info("Launching 3D WebGL viewer in browser...")
      utils::browseURL(file_url)
    }
  }

  invisible(file)
}
