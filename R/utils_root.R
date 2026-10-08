#' Advanced Root System Phenotyping Architecture
#'
#' @description
#' `root_analyze()` performs comprehensive, high-throughput phenotyping of plant
#' root systems from scans, camera images, or binary masks, surpassing the analytical
#' capabilities of standalone tools like RhizoVision Explorer and WinRHIZO:
#' * **Exact Euclidean Distance Transform (EDT)**: Computes continuous local diameters
#'   for every skeletal root pixel with sub-pixel precision.
#' * **Topological Skeletonization & Pruning**: Guo-Hall thinning followed by
#'   spurious spur pruning eliminates surface roughness noise without shortening genuine lateral roots.
#' * **Root Hierarchy Extraction**: Automatically distinguishes the Primary Root Axis
#'   (order 1 taproot / seminal root) from Secondary / Lateral roots (order 2 and higher).
#' * **True Insertion / Branching Angles**: Calculates angular divergence of lateral
#'   roots relative to parent root growth vectors in degrees (0-180 deg).
#' * **3D Root Geometry**: Computes surface area (\eqn{\sum \pi \cdot d \cdot l}) and
#'   volume (\eqn{\sum \frac{\pi}{4} \cdot d^2 \cdot l}) using cylindrical integration.
#' * **Stratified Depth Profile**: Dissects the root system into horizontal depth horizons,
#'   computing root length, surface area, volume, and RhizoVision-style root crossings per layer.
#' * **Global Network Shape**: Calculates network width, depth, width-to-depth ratio,
#'   convex hull area, network solidity, and bushiness.
#' * **Custom Diameter Size Classes**: Stratifies length, area, and volume across fine,
#'   medium, coarse, and user-defined root diameter bins.
#'
#' @param img An `Image` object, a 2D binary matrix, a file path to an image, or a list of images.
#' @param index Color index used to segment roots if `img` is an RGB image (default: `"GRAY"`).
#'   Other valid options include `"R"`, `"G"`, `"B"`, `"VARI"`, `"B-G"`, etc.
#' @param threshold Segmentation threshold for binarization. Defaults to `"otsu"` (which
#'   applies Otsu thresholding scaled by `otsu_adj`, by default +20%). Can also be a numeric
#'   threshold in `[0, 1]` (or `[0, 255]`), `"adaptive"`, or `"otsu_strict"`.
#' @param otsu_adj Numeric multiplier applied to the Otsu threshold (default: 1.2, i.e., +20% adjustment).
#'   Shifts the threshold towards the background to prevent losing fine, faint lateral roots.
#'   Set to `1.0` or use `threshold = "otsu_strict"` for standard unadjusted Otsu.
#' @param fill_hull Logical. If `TRUE`, fills internal holes/cavities within the
#'   segmented binary root mask in R before passing to the C++ engine (default: `FALSE`).
#'   Useful to prevent hollow rings in thick root axes caused by glare, but should remain `FALSE`
#'   for dense root systems to avoid bridging narrow gaps between adjacent roots.
#' @param max_order Integer specifying the maximum root branching order to classify based on
#'   node hierarchy (default: 3). `1` = Primary only; `2` = Primary and Lateral; `3` = Primary,
#'   Secondary, and Tertiary+ roots.
#' @param background Background color: `"auto"` (default, senses automatically from image corners),
#'   `"white"`, or `"black"`.
#' @param pixel_size Physical size of one pixel in the chosen `unit` (default: 1.0).
#' @param dpi Optional integer scanner/camera resolution in dots per inch. If provided,
#'   overrides `pixel_size` using `pixel_size = 25.4 / dpi` (when `unit = "mm"`).
#' @param unit Unit of measurement for lengths and coordinates (default: `"mm"`).
#' @param num_classes Number of diameter size classes to compute based on quantiles (default: 4).
#'   Can also be a numeric vector of custom diameter cut points (e.g., `c(0.5, 2.0)`).
#' @param depth_slices Integer specifying the number of horizontal depth strata to compute
#'   vertical distribution profiles (default: 10).
#' @param crown_point Optional numeric vector of length 2 `c(x, y)` specifying the coordinates
#'   (in pixels) of the root crown / stem base. If `NULL` (default), the crown is automatically
#'   detected at the topmost foreground pixel of the largest connected root component.
#'   Manually providing `crown_point = c(x, y)` is especially recommended for cuttings or
#'   flat-cut stem stumps, where a wide horizontal cut surface can otherwise split into a false "Y"
#'   bifurcation at the top.
#' @param cluster_radius Numeric radius (in pixels) within which adjacent branching junctions
#'   (degree >= 3) are clustered and merged into a single consensus node (default: 3 pixels).
#'   Thick roots often generate multi-pixel clusters of junction pixels during skeletonization;
#'   increasing `cluster_radius` (e.g., 5 to 8 pixels) merges these micro-junctions into clean,
#'   single branching nodes.
#' @param prune_length Optional integer minimum length in pixels for terminal spur branches
#'   (default: `NULL`, which defaults to 10 pixels in C++). Tips with path length shorter than
#'   this are pruned away as noise. If fine lateral roots or short emerging tips are being lost,
#'   reduce this to a smaller value (e.g., `3` to `5`).
#' @param gravitropic_alpha Numeric weight (default: 3) controlling how strongly the primary
#'   root detection favours downward-growing edges over purely diameter-based selection. Higher
#'   values (e.g., 3.0--5.0) guide the taproot identification along the downward vertical axis,
#'   even when horizontal lateral roots are thicker than the primary axis. Set to `0` to use
#'   diameter alone.
#' @param closing Integer morphological closing radius (in pixels) applied to the binary mask
#'   in C++ (default: 1). Performs morphological dilation followed by erosion. Smooths bark
#'   roughness and bridges micro-fissures in segmented roots. Set to `0` when analyzing dense
#'   root systems with closely spaced or overlapping roots to prevent artificial "webbing"
#'   (connecting membranes) between adjacent roots in the diameter map and skeleton.
#' @param fill_size Integer maximum area (in pixels) of enclosed cavities to fill internally
#'   in C++ (default: 200). Converts enclosed background regions (holes) with area <= `fill_size`
#'   pixels into root foreground, preventing false loops and hollow centers in thick root axes
#'   caused by glare or bark reflectance. Operates independently of `fill_hull`. Set to `0` to
#'   completely disable internal hole filling, avoiding the fusion of narrow gaps between roots.
#' @param radius_prune Numeric multiplier applied to local parent root radius for adaptive spur
#'   pruning (default: 1.2). Terminal branches with length shorter than
#'   `max(prune_length, parent_radius * radius_prune)` are pruned away. This eliminates surface
#'   roughness, bark bumps, and root hairs along thick roots without shortening genuine lateral
#'   roots. Lower values (e.g., 0.8 to 1.0) preserve more fine lateral roots.
#' @param dissolve_cycles Logical. If `TRUE` (default), uses Kruskal's algorithm on a
#'   Maximum Spanning Forest (MSF) to dissolve all internal loops and ladder rungs, guaranteeing
#'   an acyclic tree graph. In 2D projections, overlapping or touching roots create closed loops;
#'   the algorithm removes the edge with the lowest weight (weighted by
#'   \eqn{\text{diameter}^{3.5} \cdot \text{vertical\_ratio}^{1.5} / \sqrt{\text{length}}}),
#'   preferentially preserving thick, downward-growing parent roots and cutting horizontal or
#'   looping cross-connections. Set to `FALSE` if you wish to preserve all bridging or intersecting
#'   root segments in the skeleton and total length calculations.
#' @param plot Logical. If `TRUE` (default: `FALSE`), plots the analyzed root system.
#' @param plot_type Type of visualization: `"all"` (default 4-panel dashboard), `"skeleton"`,
#'   `"hierarchy"` (or `"classification"`), `"overlay"`, `"diameter"`, `"depth"`, `"angles"`,
#'   `"classes"`, `"profile"` (or `"envelope"`), `"topology"`, `"diameter_dist"`, or `"orders"`.
#' @param verbose Logical. If `TRUE` (default), displays progress messages.
#' @param ... Additional arguments passed to [image_binary()].
#'
#' @return An object of class `root_analysis` (or `root_analysis_list` for multiple images)
#' containing:
#' * `summary`: Single-row `data.frame` with all phenotypic traits.
#' * `roots`: Detailed `data.frame` of all individual root segments (order, length, avg diameter,
#'   surface area, volume, branching angle, insertion depth).
#' * `nodes`: Coordinates, degrees, and types (`crown`, `junction`, `tip`) of all network nodes.
#' * `paths`: List of matrices containing (x, y) coordinates for each curved root segment.
#' * `diameter_classes`: Distribution of length, surface area, and volume across diameter classes.
#' * `depth_profile`: Stratified depth distribution with lengths, volumes, and root crossings.
#' * `skeleton`: 2D binary matrix of the pruned skeleton.
#' * `diameter_map`: 2D matrix of local root diameters.
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @references
#' Seethepalli, A., Guo, H., Liu, X., Ruff, T., & York, L. M. (2021).
#' RhizoVision Explorer: open-source software for generalized root phenotyping using
#' image analysis. Plant Physiology, 187(2), 739-757.
#'
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#'
#' # Analyze a root image scanned at 600 DPI
#' res <- root_analyze("root_sample.jpg", dpi = 600)
#'
#' # View global phenotypic summary
#' res$summary
#'
#' # Plot 4-panel diagnostic dashboard
#' plot(res)
#'
#' # Plot continuous diameter heatmap
#' plot(res, type = "diameter")
#' }
root_analyze <- function(img,
                         index = "GRAY",
                         threshold = "otsu",
                         otsu_adj = 1.2,
                         fill_hull = FALSE,
                         max_order = 3,
                         background = "auto",
                         pixel_size = 1.0,
                         dpi = NULL,
                         unit = "mm",
                         num_classes = 4,
                         depth_slices = 10,
                         crown_point = NULL,
                         cluster_radius = 3,
                         prune_length = NULL,
                         gravitropic_alpha = 3,
                         closing = 1,
                         fill_size = 200,
                         radius_prune = 1.2,
                         dissolve_cycles = TRUE,
                         plot = FALSE,
                         plot_type = "all",
                         verbose = TRUE,
                         ...) {

  # Calibrate pixel size from DPI if provided
  if (!is.null(dpi)) {
    if (unit == "mm") {
      pixel_size <- 25.4 / dpi
    } else if (unit == "cm") {
      pixel_size <- 2.54 / dpi
    } else if (unit == "inch") {
      pixel_size <- 1.0 / dpi
    }
  }

  # Batch processing for lists
  if (is.list(img) && !is_image(img)) {
    res_list <- lapply(seq_along(img), function(i) {
      if (verbose) {
        cli::cli_progress_step("Analyzing root image {i}/{length(img)}...")
      }
      root_analyze(
        img = img[[i]],
        index = index,
        threshold = threshold,
        otsu_adj = otsu_adj,
        fill_hull = fill_hull,
        max_order = max_order,
        background = background,
        pixel_size = pixel_size,
        dpi = NULL,
        unit = unit,
        num_classes = num_classes,
        depth_slices = depth_slices,
        crown_point = crown_point,
        cluster_radius = cluster_radius,
        prune_length = prune_length,
        gravitropic_alpha = gravitropic_alpha,
        closing = closing,
        fill_size = fill_size,
        radius_prune = radius_prune,
        dissolve_cycles = dissolve_cycles,
        plot = FALSE,
        verbose = FALSE,
        ...
      )
    })
    names(res_list) <- names(img)
    class(res_list) <- c("root_analysis_list", "list")
    return(res_list)
  }

  # Read from file if character
  if (is.character(img)) {
    if (file.exists(img)) {
      img <- image_import(img)
    } else {
      cli::cli_abort("Image file not found: {.val {img}}")
    }
  }

  # Convert to binary matrix
  bin_mat <- NULL
  orig_img <- NULL
  th_val <- NA_real_
  th_otsu_orig <- NA_real_

  if (is.matrix(img)) {
    if (is.logical(img)) {
      bin_mat <- img * 1.0
    } else if (all(img %in% c(0, 1))) {
      bin_mat <- img * 1.0
    } else {
      th_val <- if (is.numeric(threshold)) threshold else 0.5
      bin_mat <- (img > th_val) * 1.0
    }
  } else if (is_image(img)) {
    orig_img <- img
    dims <- dim(img)
    w_img <- dims[1]
    h_img <- dims[2]

    # Extract 2D normalized grayscale intensity matrix [0, 1]
    gray_mat <- matrix(0, w_img, h_img)
    if (length(dims) == 2 || (length(dims) == 3 && dims[3] == 1)) {
      arr <- image_data(img, type = "numeric")
      if (length(dim(arr)) == 3) arr <- arr[, , 1]
      gray_mat <- arr
    } else {
      arr <- image_data(img, type = "numeric")
      idx_upper <- toupper(index)
      if (idx_upper == "R") {
        gray_mat <- arr[, , 1]
      } else if (idx_upper == "G") {
        gray_mat <- arr[, , 2]
      } else if (idx_upper == "B") {
        gray_mat <- arr[, , 3]
      } else if (idx_upper %in% c("GRAY", "LUMINANCE")) {
        gray_mat <- 0.299 * arr[, , 1] + 0.587 * arr[, , 2] + 0.114 * arr[, , 3]
      } else {
        ind_res <- try(image_index(img, index = index, plot = FALSE, verbose = FALSE)[[1]], silent = TRUE)
        if (!inherits(ind_res, "try-error") && !is.null(ind_res)) {
          gray_mat <- image_data(ind_res, type = "numeric")
        } else {
          gray_mat <- 0.299 * arr[, , 1] + 0.587 * arr[, , 2] + 0.114 * arr[, , 3]
        }
      }
    }

    # Normalize gray_mat to [0, 1] if needed
    min_v <- min(gray_mat, na.rm = TRUE)
    max_v <- max(gray_mat, na.rm = TRUE)
    if (max_v > min_v && (min_v < 0 || max_v > 1)) {
      gray_mat <- (gray_mat - min_v) / (max_v - min_v)
    }

    # Determine background (light vs dark)
    if (background == "auto") {
      c_w <- min(10, w_img)
      c_h <- min(10, h_img)
      corners <- c(
        gray_mat[1:c_w, 1:c_h],
        gray_mat[(w_img - c_w + 1):w_img, 1:c_h],
        gray_mat[1:c_w, (h_img - c_h + 1):h_img],
        gray_mat[(w_img - c_w + 1):w_img, (h_img - c_h + 1):h_img]
      )
      bg_val <- stats::median(corners, na.rm = TRUE)
      is_light_bg <- (bg_val > 0.5)
    } else {
      is_light_bg <- (background == "white")
    }

    # Segment root foreground
    th_type <- if (is.character(threshold)) tolower(threshold) else "numeric"
    th_val <- 0.5
    th_otsu_orig <- NA_real_

    if (th_type == "adaptive") {
      ws <- as.integer(min(w_img, h_img) / 5)
      if (ws %% 2 == 0) ws <- ws + 1
      if (ws < 3) ws <- 3
      adapt_mat <- threshold_adaptive(gray_mat, k = 0.15, windowsize = ws)
      if (is_light_bg) {
        bin_mat <- (!adapt_mat) * 1.0
      } else {
        bin_mat <- (adapt_mat) * 1.0
      }
    } else {
      if (th_type %in% c("otsu", "otsu_adj")) {
        th_val <- help_otsu(gray_mat)
        if (th_val > 1.0 && max(gray_mat, na.rm = TRUE) <= 1.0) {
          th_val <- th_val / 255.0
        }
        th_otsu_orig <- th_val
        if (is.numeric(otsu_adj) && otsu_adj > 0) {
          if (is_light_bg) {
            # In light background, roots are darker (gray < th_val).
            # Multiplying Otsu threshold by otsu_adj (default 1.2, +20%) captures fainter root tips.
            th_val <- min(0.999, th_val * otsu_adj)
          } else {
            # In dark background, roots are brighter (gray > th_val).
            th_val <- max(0.001, th_val / otsu_adj)
          }
        }
      } else if (th_type == "otsu_strict") {
        th_val <- help_otsu(gray_mat)
        if (th_val > 1.0 && max(gray_mat, na.rm = TRUE) <= 1.0) {
          th_val <- th_val / 255.0
        }
        th_otsu_orig <- th_val
      } else if (is.numeric(threshold)) {
        th_val <- threshold
        if (th_val > 1.0 && max(gray_mat, na.rm = TRUE) <= 1.0) {
          th_val <- th_val / 255.0
        }
      } else {
        th_val <- 0.5
      }

      if (is_light_bg) {
        bin_mat <- (gray_mat < th_val) * 1.0
      } else {
        bin_mat <- (gray_mat > th_val) * 1.0
      }
    }
  } else {
    cli::cli_abort("Unsupported input type for {.fn root_analyze}.")
  }

  # Fill small internal holes/cavities to eliminate false loops in thick root axes
  if (isTRUE(fill_hull) && !is.null(bin_mat)) {
    bin_mat <- help_binary_filters_cpp(bin_mat, fill_hull = TRUE, max_size = 500) * 1.0
  }

  # Ensure binary matrix format
  storage.mode(bin_mat) <- "double"

  cluster_rad_px <- as.numeric(cluster_radius)
  if (cluster_rad_px < 0) cluster_rad_px <- 0

  prune_px <- if (!is.null(prune_length)) as.integer(prune_length) else 10L
  if (prune_px < 0) prune_px <- 0

  cx_in <- -1.0
  cy_in <- -1.0
  if (!is.null(crown_point) && length(crown_point) >= 2) {
    cx_in <- as.double(crown_point[1])
    cy_in <- as.double(crown_point[2])
  }

  if (verbose) {
    cli::cli_progress_step("Extracting root topological graph and Euclidean distance transform...")
  }

  custom_bins <- NULL
  n_cls <- 4
  if (is.numeric(num_classes)) {
    if (length(num_classes) == 1) {
      n_cls <- as.integer(max(2, num_classes))
    } else {
      custom_bins <- as.numeric(sort(num_classes))
      n_cls <- as.integer(length(custom_bins) + 1)
    }
  }

  res_cpp <- analyze_root_system_cpp(
    binary_mat = bin_mat,
    pixel_size = pixel_size,
    prune_length_px = prune_px,
    cluster_radius_px = cluster_rad_px,
    crown_x_in = cx_in,
    crown_y_in = cy_in,
    num_classes = n_cls,
    diameter_bins = custom_bins,
    n_depth_slices = as.integer(depth_slices),
    max_order = as.integer(max_order),
    gravitropic_alpha = as.double(gravitropic_alpha),
    closing_rad = as.integer(closing),
    fill_size = as.integer(fill_size),
    rad_fac = as.double(radius_prune),
    dissolve_cycles = isTRUE(dissolve_cycles)
  )

  # Package output
  out <- list(
    summary = res_cpp$summary,
    roots = res_cpp$roots,
    nodes = res_cpp$nodes,
    paths = res_cpp$paths,
    diameter_classes = {
      dc <- res_cpp$diameter_classes
      if (!is.null(dc$class) && unit != "mm") {
        dc$class <- gsub("mm", unit, dc$class)
      }
      dc
    },
    depth_profile = res_cpp$depth_profile,
    skeleton = as_image(res_cpp$skeleton, colormode = "Grayscale"),
    diameter_map = res_cpp$diameter_map,
    binary = if (!is.null(bin_mat)) as_image(bin_mat, colormode = "Grayscale") else NULL,
    original_image = orig_img,
    crown_node = res_cpp$crown_node,
    threshold = th_val,
    threshold_otsu = th_otsu_orig,
    max_order = max_order,
    pixel_size = pixel_size,
    unit = unit
  )

  class(out) <- "root_analysis"

  if (verbose) {
    cli::cli_alert_success(
      "Root phenotyping complete: {round(out$summary$total_length, 2)} {unit} total length | {out$summary$n_tips} tips | {out$summary$n_forks} forks"
    )
  }

  if (isTRUE(plot)) {
    plot(out, type = plot_type)
  }

  return(out)
}

#' Print method for root_analysis objects
#' @param x A `root_analysis` object.
#' @param ... Ignored.
#' @method print root_analysis
#' @export
print.root_analysis <- function(x, ...) {
  s <- x$summary
  u <- x$unit

  lat_bullets <- c()
  if (!is.null(s$secondary_length) && !is.null(s$tertiary_length) && s$n_tertiary > 0) {
    lat_bullets <- c(
      "*" = paste0(cli::style_bold("Secondary Roots: "), s$n_secondary, " branches (Total length: ", round(s$secondary_length, 2), " ", u, ")"),
      "*" = paste0(cli::style_bold("Tertiary+ Roots: "), s$n_tertiary, " branches (Total length: ", round(s$tertiary_length, 2), " ", u, ")")
    )
  } else {
    lat_bullets <- c(
      "*" = paste0(cli::style_bold("Lateral Roots: "), s$n_lateral, " branches (Total length: ", round(s$lateral_length, 2), " ", u, ")")
    )
  }

  cli::cli_rule(left = cli::col_cyan("Root System Phenotyping Summary"))
  cli::cli_bullets(c(
    "*" = paste0(cli::style_bold("Total Root Length: "), round(s$total_length, 2), " ", u),
    "*" = paste0(cli::style_bold("Projected Area: "), round(s$projected_area, 2), " ", u, "^2"),
    "*" = paste0(cli::style_bold("Surface Area: "), round(s$surface_area, 2), " ", u, "^2"),
    "*" = paste0(cli::style_bold("Volume: "), round(s$volume, 2), " ", u, "^3"),
    "*" = paste0(cli::style_bold("Average Diameter: "), round(s$avg_diameter, 3), " ", u),
    "*" = paste0(cli::style_bold("Primary Root Length: "), round(s$primary_length, 2), " ", u),
    lat_bullets,
    "*" = paste0(cli::style_bold("Lateral Density: "), round(s$lateral_density_per_cm, 2), " roots/cm"),
    "*" = paste0(cli::style_bold("Mean Branch Angle: "), if (!is.na(s$mean_branch_angle)) round(s$mean_branch_angle, 1) else "NA", " deg"),
    "*" = paste0(cli::style_bold("Architecture: "), "Width: ", round(s$max_width, 2), " ", u, " | Depth: ", round(s$max_depth, 2), " ", u, " | W:D Ratio: ", round(s$width_depth_ratio, 2)),
    "*" = paste0(cli::style_bold("Network Topology: "), s$n_tips, " tips | ", s$n_forks, " forks | Solidity: ", round(s$solidity, 3))
  ))
  cli::cli_rule()
  invisible(x)
}

#' Plot method for root_analysis objects
#'
#' @param x An object of class `root_analysis`.
#' @param type Type of plot:
#'   * `"all"`: 4-panel diagnostic phenotyping sheet.
#'   * `"overlay"`: Original image with translucent diameter heatmap and topological nodes overlaid.
#'   * `"diameter"`: Continuous diameter heatmap (with optional background and nodes).
#'   * `"skeleton"`: Topological skeleton with primary/lateral root paths and nodes.
#'   * `"depth"`: Stratified depth layering barplot with root crossings.
#'   * `"angles"`: Branching angles distribution histogram.
#'   * `"classes"`: Diameter classes breakdown.
#' @param type Type of plot:
#'   * `"all"` (default): 4-panel diagnostic phenotyping dashboard.
#'   * `"skeleton"`: Topological skeleton map with primary/lateral curved root paths and nodes.
#'   * `"overlay"`: Original image with translucent diameter heatmap and topological nodes overlaid.
#'   * `"diameter"`: Continuous diameter heatmap (with optional background and nodes).
#'   * `"depth"`: Stratified depth layering barplot with root crossings.
#'   * `"angles"`: Branching angles distribution histogram.
#'   * `"classes"`: Diameter classes breakdown.
#' @param which Optional alias for `type`.
#' @param background Background for skeleton, diameter, or overlay plot: `"image"` (original image
#'   if available), `"white"`, or `"black"`. Defaults to `"image"` if original image exists, otherwise `"white"`.
#' @param nodes Logical. If `TRUE`, overlays the topological nodes (Crown, Branch points, and Root tips).
#' @param alpha Numeric transparency in `[0, 1]` for the diameter heatmap overlay (default: 0.65).
#' @param col Optional vector of colors for the diameter heatmap palette.
#' @param lwd Line width for drawing root skeleton segments (default: 2.0).
#' @param legend_pos Position of the legend: `"auto"` (default, automatically places legend in
#'   the corner with fewest/no root obstacles), `"topleft"`, `"topright"`, `"bottomleft"`,
#'   `"bottomright"`, or `FALSE` to omit legends.
#' @param cex_legend Character expansion factor for the legend (default: `NULL`, which
#'   automatically uses 0.65 for multi-panel dashboards and 0.75 for individual plots).
#' @param axes Logical. If `TRUE` (default), displays coordinate axes. If `FALSE`, axes are suppressed.
#' @param grid Logical. If `TRUE` (default), displays standardized subtle background grid lines
#'   on analytical and distribution charts. If `FALSE`, grid lines are suppressed.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @method plot root_analysis
#' @export
plot.root_analysis <- function(x,
                               type = c("all", "skeleton", "hierarchy", "classification", "overlay", "diameter", "depth", "angles", "classes", "profile", "envelope", "shape", "topology", "angles_depth", "tropism", "diameter_dist", "density", "orders", "order_metrics", "binary"),
                               which = NULL,
                               background = if (!is.null(x$original_image)) "image" else "white",
                               nodes = TRUE,
                               alpha = 0.85,
                               col = NULL,
                               lwd = 2.0,
                               cex_nodes = 1.0,
                               legend_pos = "auto",
                               cex_legend = NULL,
                               axes = TRUE,
                               grid = TRUE,
                               ...) {
  if (!is.null(which)) {
    type <- which
  }
  type <- match.arg(type)
  if (type == "classification") type <- "hierarchy"
  if (type %in% c("envelope", "shape")) type <- "profile"
  if (type %in% c("angles_depth", "tropism")) type <- "topology"
  if (type %in% c("density")) type <- "diameter_dist"
  if (type %in% c("order_metrics")) type <- "orders"

  if (is.null(cex_legend)) {
    is_mf <- (type == "all") || (prod(graphics::par("mfrow")) > 1L) || (prod(graphics::par("mfcol")) > 1L)
    cex_legend <- if (is_mf) 0.65 else 0.75
  }

  if (type == "all") {
    op <- graphics::par(mfrow = c(2, 2), mar = c(4, 5, 3.5, 2))
    on.exit(graphics::par(op), add = TRUE)
    .plot_root_skeleton(x, background = background, nodes = nodes, cex_nodes = cex_nodes, lwd = lwd, legend_pos = legend_pos, cex_legend = cex_legend, main = "Root Topological Skeleton", axes = axes, grid = grid, ...)
    .plot_root_hierarchy(x, background = background, nodes = nodes, col = col, cex_nodes = cex_nodes, legend_pos = legend_pos, cex_legend = cex_legend, main = "Root Classification", axes = axes, grid = grid, ...)
    .plot_root_depth(x, main = "Depth Stratification & Crossings", grid = grid, ...)
    .plot_root_classes(x, main = "Diameter Size Classes", grid = grid, ...)
    return(invisible(x))
  }

  if (type == "overlay") {
    dots <- list(...)
    main_overlay <- if ("main" %in% names(dots)) dots$main else "Diameter Heatmap & Nodes Overlay"
    dots$main <- NULL
    do.call(.plot_root_diameter_map, c(list(x = x, background = background, nodes = nodes, alpha = alpha, col = col, cex_nodes = cex_nodes, legend_pos = legend_pos, cex_legend = cex_legend, main = main_overlay, axes = axes), dots))
  } else if (type == "skeleton") {
    .plot_root_skeleton(x, background = background, nodes = nodes, cex_nodes = cex_nodes, lwd = lwd, legend_pos = legend_pos, cex_legend = cex_legend, axes = axes, grid = grid, ...)
  } else if (type == "hierarchy") {
    .plot_root_hierarchy(x, background = background, nodes = nodes, col = col, cex_nodes = cex_nodes, legend_pos = legend_pos, cex_legend = cex_legend, axes = axes, grid = grid, ...)
  } else if (type == "diameter") {
    .plot_root_diameter_map(x, background = background, nodes = nodes, alpha = alpha, col = col, cex_nodes = cex_nodes, legend_pos = legend_pos, cex_legend = cex_legend, axes = axes, ...)
  } else if (type == "depth") {
    .plot_root_depth(x, grid = grid, ...)
  } else if (type == "angles") {
    .plot_root_angles(x, legend_pos = legend_pos, cex_legend = cex_legend, grid = grid, ...)
  } else if (type == "classes") {
    .plot_root_classes(x, grid = grid, ...)
  } else if (type == "profile") {
    .plot_root_profile(x, axes = axes, legend_pos = legend_pos, cex_legend = cex_legend, grid = grid, ...)
  } else if (type == "topology") {
    .plot_root_topology(x, axes = axes, legend_pos = legend_pos, cex_legend = cex_legend, grid = grid, ...)
  } else if (type == "diameter_dist") {
    .plot_root_diameter_dist(x, axes = axes, legend_pos = legend_pos, cex_legend = cex_legend, grid = grid, ...)
  } else if (type == "orders") {
    .plot_root_orders(x, grid = grid, ...)
  } else if (type == "binary") {
    if (!is.null(x$binary)) {
      plot(x$binary, axes = axes, main = "Root Binary Mask", ...)
    } else {
      cli::cli_warn("No binary mask stored in the root_analysis object.")
    }
  }

  invisible(x)
}

# Internal helper: Draw a fixed, non-clipping legend for any coordinate system and window size
.draw_root_legend <- function(pos = "auto",
                              labels,
                              cols = NULL,
                              types = c("point"),
                              pchs = NA,
                              lwds = NA,
                              pt_cexs = 1.0,
                              fills = NULL,
                              title = NULL,
                              cex = 0.75,
                              bg = grDevices::adjustcolor("white", alpha.f = 0.90),
                              border = "gray70",
                              margin_npc = 0.025,
                              avoid_pts = NULL,
                              ...) {
  if (isFALSE(pos) || is.null(pos)) return(invisible())

  n <- length(labels)
  if (n == 0) return(invisible())

  if (is.null(cols)) cols <- rep("black", n) else cols <- rep_len(cols, n)
  types <- rep_len(types, n)
  pchs <- rep_len(pchs, n)
  lwds <- rep_len(lwds, n)
  pt_cexs <- rep_len(pt_cexs, n)
  if (!is.null(fills)) fills <- rep_len(fills, n)

  pin <- graphics::par("pin")
  if (any(pin <= 0)) return(invisible())

  # Calculate sizes in inches to convert to NPC (Normalized Plot Coordinates: 0 to 1)
  char_w_in <- graphics::strwidth("M", units = "inches", cex = cex)
  char_h_in <- graphics::strheight("M", units = "inches", cex = cex)
  row_h_in <- char_h_in * 1.70

  max_txt_w_in <- max(graphics::strwidth(labels, units = "inches", cex = cex))
  has_title <- !is.null(title) && nzchar(title)
  if (has_title) {
    max_txt_w_in <- max(max_txt_w_in, graphics::strwidth(title, units = "inches", cex = cex * 1.05, font = 2))
  }

  sym_w_in <- char_w_in * 2.2
  pad_x_in <- char_w_in * 0.9
  pad_y_in <- char_h_in * 0.7
  title_h_in <- if (has_title) row_h_in * 1.15 else 0

  box_w_in <- pad_x_in * 2 + sym_w_in + char_w_in * 0.5 + max_txt_w_in
  box_h_in <- pad_y_in * 2 + title_h_in + n * row_h_in

  # Convert to NPC (fractions of plot area)
  box_w_npc <- min(0.96, box_w_in / pin[1])
  box_h_npc <- min(0.96, box_h_in / pin[2])
  pad_x_npc <- pad_x_in / pin[1]
  pad_y_npc <- pad_y_in / pin[2]
  sym_w_npc <- sym_w_in / pin[1]
  row_h_npc <- row_h_in / pin[2]
  title_h_npc <- title_h_in / pin[2]

  pos_str <- tolower(as.character(pos[1]))
  cands <- c("bottomright", "topright", "bottomleft", "topleft")

  if (pos_str == "auto") {
    best_cand <- "bottomright"
    if (!is.null(avoid_pts) && nrow(avoid_pts) > 0) {
      pts_x_npc <- tryCatch(graphics::grconvertX(avoid_pts[, 1], from = "user", to = "npc"), error = function(e) NULL)
      pts_y_npc <- tryCatch(graphics::grconvertY(avoid_pts[, 2], from = "user", to = "npc"), error = function(e) NULL)
      if (!is.null(pts_x_npc) && !is.null(pts_y_npc)) {
        min_overlap <- Inf
        for (cand in cands) {
          cx_left <- if (grepl("right", cand)) 1 - margin_npc - box_w_npc else margin_npc
          cx_right <- cx_left + box_w_npc
          cy_bot <- if (grepl("top", cand)) 1 - margin_npc - box_h_npc else margin_npc
          cy_top <- cy_bot + box_h_npc

          overlap_cnt <- sum(pts_x_npc >= cx_left & pts_x_npc <= cx_right &
                             pts_y_npc >= cy_bot & pts_y_npc <= cy_top, na.rm = TRUE)
          if (overlap_cnt < min_overlap) {
            min_overlap <- overlap_cnt
            best_cand <- cand
          }
        }
      }
    }
    pos_str <- best_cand
  }

  # Position box in NPC (0 to 1) so it NEVER clips outside plot area
  if (grepl("right", pos_str)) {
    x_right_npc <- 1 - margin_npc
    x_left_npc <- max(margin_npc, x_right_npc - box_w_npc)
  } else if (grepl("center", pos_str)) {
    x_left_npc <- max(margin_npc, (1 - box_w_npc) / 2)
    x_right_npc <- min(1 - margin_npc, x_left_npc + box_w_npc)
  } else { # left
    x_left_npc <- margin_npc
    x_right_npc <- min(1 - margin_npc, x_left_npc + box_w_npc)
  }

  if (grepl("top", pos_str)) {
    y_top_npc <- 1 - margin_npc
    y_bot_npc <- max(margin_npc, y_top_npc - box_h_npc)
  } else if (grepl("center", pos_str)) {
    y_bot_npc <- max(margin_npc, (1 - box_h_npc) / 2)
    y_top_npc <- min(1 - margin_npc, y_bot_npc + box_h_npc)
  } else { # bottom
    y_bot_npc <- margin_npc
    y_top_npc <- min(1 - margin_npc, y_bot_npc + box_h_npc)
  }

  # Convert box boundaries from NPC to active user coordinates
  x0 <- graphics::grconvertX(x_left_npc, from = "npc", to = "user")
  x1 <- graphics::grconvertX(x_right_npc, from = "npc", to = "user")
  y0 <- graphics::grconvertY(y_bot_npc, from = "npc", to = "user")
  y1 <- graphics::grconvertY(y_top_npc, from = "npc", to = "user")

  # Draw background box
  graphics::rect(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1),
                 col = bg, border = border, lwd = 0.8)

  cur_y_npc <- y_top_npc - pad_y_npc
  if (has_title) {
    t_y <- graphics::grconvertY(cur_y_npc - title_h_npc * 0.4, from = "npc", to = "user")
    t_x <- graphics::grconvertX(x_left_npc + pad_x_npc, from = "npc", to = "user")
    graphics::text(t_x, t_y, title, font = 2, cex = cex * 1.05, adj = c(0, 0.5))
    cur_y_npc <- cur_y_npc - title_h_npc
  }

  sym_cx_npc <- x_left_npc + pad_x_npc + sym_w_npc / 2
  txt_x_npc <- x_left_npc + pad_x_npc + sym_w_npc + char_w_in * 0.5 / pin[1]

  sym_cx <- graphics::grconvertX(sym_cx_npc, from = "npc", to = "user")
  txt_x <- graphics::grconvertX(txt_x_npc, from = "npc", to = "user")

  for (i in seq_len(n)) {
    line_y_npc <- cur_y_npc - (i - 0.5) * row_h_npc
    y_i <- graphics::grconvertY(line_y_npc, from = "npc", to = "user")
    t <- types[i]

    if (t == "line") {
      lx0 <- graphics::grconvertX(sym_cx_npc - sym_w_npc * 0.42, from = "npc", to = "user")
      lx1 <- graphics::grconvertX(sym_cx_npc + sym_w_npc * 0.42, from = "npc", to = "user")
      l_w <- if (!is.na(lwds[i])) lwds[i] else 2.0
      graphics::lines(c(lx0, lx1), c(y_i, y_i), col = cols[i], lwd = l_w)
    } else if (t == "point") {
      p_lwd <- if (!is.na(pchs[i]) && pchs[i] %in% c(8, 21, 22, 23, 24, 25)) 1.5 else 1.0
      f_col <- if (!is.null(fills) && !is.na(fills[i])) fills[i] else cols[i]
      graphics::points(sym_cx, y_i, pch = pchs[i], col = cols[i], bg = f_col,
                       cex = pt_cexs[i], lwd = p_lwd)
    } else if (t %in% c("polygon", "fill")) {
      bx0 <- graphics::grconvertX(sym_cx_npc - sym_w_npc * 0.40, from = "npc", to = "user")
      bx1 <- graphics::grconvertX(sym_cx_npc + sym_w_npc * 0.40, from = "npc", to = "user")
      by0 <- graphics::grconvertY(line_y_npc - row_h_npc * 0.30, from = "npc", to = "user")
      by1 <- graphics::grconvertY(line_y_npc + row_h_npc * 0.30, from = "npc", to = "user")
      f_col <- if (!is.null(fills) && !is.na(fills[i])) fills[i] else cols[i]
      graphics::rect(min(bx0, bx1), min(by0, by1), max(bx0, bx1), max(by0, by1),
                     col = f_col, border = cols[i], lwd = 1.0)
    } else if (t == "both") {
      lx0 <- graphics::grconvertX(sym_cx_npc - sym_w_npc * 0.42, from = "npc", to = "user")
      lx1 <- graphics::grconvertX(sym_cx_npc + sym_w_npc * 0.42, from = "npc", to = "user")
      l_w <- if (!is.na(lwds[i])) lwds[i] else 2.0
      graphics::lines(c(lx0, lx1), c(y_i, y_i), col = cols[i], lwd = l_w)
      f_col <- if (!is.null(fills) && !is.na(fills[i])) fills[i] else cols[i]
      graphics::points(sym_cx, y_i, pch = pchs[i], col = cols[i], bg = f_col,
                       cex = pt_cexs[i], lwd = 1.5)
    }

    graphics::text(txt_x, y_i, labels[i], adj = c(0, 0.5), cex = cex, col = "black")
  }
}

# Internal helper: Draw a continuous gradient colorbar legend with optional node symbols
.draw_root_colorbar <- function(pos = "auto",
                                range = c(0, 1.0),
                                cols = NULL,
                                unit = "mm",
                                title = "Diameter",
                                cex = 0.75,
                                node_items = NULL,
                                bg = grDevices::adjustcolor("white", alpha.f = 0.92),
                                border = "gray75",
                                margin_npc = 0.025,
                                avoid_pts = NULL) {
  if (isFALSE(pos) || is.null(pos)) return(invisible())

  if (is.null(cols)) {
    cols <- grDevices::colorRampPalette(c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"))(100)
  }

  if (length(range) < 2 || range[2] <= range[1]) {
    range <- c(0, max(1.0, range[1]))
  }

  pin <- graphics::par("pin")
  if (any(pin <= 0)) return(invisible())

  char_w_in <- graphics::strwidth("M", units = "inches", cex = cex)
  char_h_in <- graphics::strheight("M", units = "inches", cex = cex)

  # Colorbar dimensions in inches
  bar_w_in <- char_w_in * 1.5
  bar_h_in <- char_h_in * 6.5
  tick_len_in <- char_w_in * 0.45

  # 3 clean reference ticks: min, mid, max
  max_d <- range[2]
  min_d <- range[1]
  mid_d <- (min_d + max_d) / 2
  ticks <- c(min_d, round(mid_d, 2), round(max_d, 2))
  tick_labels <- paste0(ticks, if (!is.null(unit) && nzchar(unit)) paste0(" ", unit) else "")

  max_tick_w_in <- max(graphics::strwidth(tick_labels, units = "inches", cex = cex * 0.85))

  has_title <- !is.null(title) && nzchar(title)
  title_txt <- title
  title_w_in <- if (has_title) graphics::strwidth(title_txt, units = "inches", cex = cex * 1.05, font = 2) else 0
  title_h_in <- if (has_title) char_h_in * 1.8 else 0

  # Node items dimensions if present
  n_nodes <- if (!is.null(node_items)) nrow(node_items) else 0
  node_row_h_in <- char_h_in * 1.6
  node_section_h_in <- if (n_nodes > 0) (n_nodes * node_row_h_in + char_h_in * 0.8) else 0
  max_node_txt_w_in <- if (n_nodes > 0) max(graphics::strwidth(node_items$label, units = "inches", cex = cex)) else 0
  node_sym_w_in <- if (n_nodes > 0) char_w_in * 2.0 else 0
  node_w_in <- if (n_nodes > 0) (node_sym_w_in + max_node_txt_w_in) else 0

  pad_x_in <- char_w_in * 1.0
  pad_y_in <- char_h_in * 0.8

  content_w_in <- max(title_w_in, bar_w_in + tick_len_in + char_w_in * 0.5 + max_tick_w_in, node_w_in)
  box_w_in <- pad_x_in * 2 + content_w_in
  box_h_in <- pad_y_in * 2 + title_h_in + bar_h_in + node_section_h_in

  # NPC conversion
  box_w_npc <- min(0.96, box_w_in / pin[1])
  box_h_npc <- min(0.96, box_h_in / pin[2])

  pos_str <- tolower(as.character(pos[1]))
  cands <- c("bottomright", "topright", "bottomleft", "topleft")
  if (pos_str == "auto") {
    best_cand <- "bottomright"
    if (!is.null(avoid_pts) && nrow(avoid_pts) > 0) {
      pts_x_npc <- tryCatch(graphics::grconvertX(avoid_pts[, 1], from = "user", to = "npc"), error = function(e) NULL)
      pts_y_npc <- tryCatch(graphics::grconvertY(avoid_pts[, 2], from = "user", to = "npc"), error = function(e) NULL)
      if (!is.null(pts_x_npc) && !is.null(pts_y_npc)) {
        min_overlap <- Inf
        for (cand in cands) {
          cx_left <- if (grepl("right", cand)) 1 - margin_npc - box_w_npc else margin_npc
          cx_right <- cx_left + box_w_npc
          cy_bot <- if (grepl("top", cand)) 1 - margin_npc - box_h_npc else margin_npc
          cy_top <- cy_bot + box_h_npc
          overlap_cnt <- sum(pts_x_npc >= cx_left & pts_x_npc <= cx_right &
                             pts_y_npc >= cy_bot & pts_y_npc <= cy_top, na.rm = TRUE)
          if (overlap_cnt < min_overlap) {
            min_overlap <- overlap_cnt
            best_cand <- cand
          }
        }
      }
    }
    pos_str <- best_cand
  }

  if (grepl("right", pos_str)) {
    x_right_npc <- 1 - margin_npc
    x_left_npc <- max(margin_npc, x_right_npc - box_w_npc)
  } else if (grepl("center", pos_str)) {
    x_left_npc <- max(margin_npc, (1 - box_w_npc) / 2)
    x_right_npc <- min(1 - margin_npc, x_left_npc + box_w_npc)
  } else {
    x_left_npc <- margin_npc
    x_right_npc <- min(1 - margin_npc, x_left_npc + box_w_npc)
  }

  if (grepl("top", pos_str)) {
    y_top_npc <- 1 - margin_npc
    y_bot_npc <- max(margin_npc, y_top_npc - box_h_npc)
  } else if (grepl("center", pos_str)) {
    y_bot_npc <- max(margin_npc, (1 - box_h_npc) / 2)
    y_top_npc <- min(1 - margin_npc, y_bot_npc + box_h_npc)
  } else {
    y_bot_npc <- margin_npc
    y_top_npc <- min(1 - margin_npc, y_bot_npc + box_h_npc)
  }

  # Draw background box in user coordinates
  x0 <- graphics::grconvertX(x_left_npc, from = "npc", to = "user")
  x1 <- graphics::grconvertX(x_right_npc, from = "npc", to = "user")
  y0 <- graphics::grconvertY(y_bot_npc, from = "npc", to = "user")
  y1 <- graphics::grconvertY(y_top_npc, from = "npc", to = "user")
  graphics::rect(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1),
                 col = bg, border = border, lwd = 1.0)

  # Title
  cur_y_in <- (y_top_npc * pin[2]) - pad_y_in
  cur_x_in <- (x_left_npc * pin[1]) + pad_x_in

  if (has_title) {
    t_y <- graphics::grconvertY((cur_y_in - char_h_in * 0.5) / pin[2], from = "npc", to = "user")
    t_x <- graphics::grconvertX(cur_x_in / pin[1], from = "npc", to = "user")
    graphics::text(t_x, t_y, title_txt, font = 2, cex = cex * 1.05, adj = c(0, 0.5), col = "gray15")
    cur_y_in <- cur_y_in - title_h_in
  }

  # Colorbar rectangle in inches
  bar_top_in <- cur_y_in
  bar_bot_in <- cur_y_in - bar_h_in
  bar_left_in <- cur_x_in
  bar_right_in <- cur_x_in + bar_w_in

  bar_x0 <- graphics::grconvertX(bar_left_in / pin[1], from = "npc", to = "user")
  bar_x1 <- graphics::grconvertX(bar_right_in / pin[1], from = "npc", to = "user")
  bar_y_top <- graphics::grconvertY(bar_top_in / pin[2], from = "npc", to = "user")
  bar_y_bot <- graphics::grconvertY(bar_bot_in / pin[2], from = "npc", to = "user")

  # Vertical ramp matrix: row 1 is top (max value, yellow), last row is bottom (0, purple)
  ramp_mat <- matrix(rev(cols), ncol = 1)
  graphics::rasterImage(ramp_mat, xleft = min(bar_x0, bar_x1), ybottom = bar_y_bot,
                         xright = max(bar_x0, bar_x1), ytop = bar_y_top, interpolate = TRUE)
  graphics::rect(min(bar_x0, bar_x1), min(bar_y_bot, bar_y_top),
                 max(bar_x0, bar_x1), max(bar_y_bot, bar_y_top),
                 border = "gray40", lwd = 0.8)

  # Ticks and tick labels
  tick_x_start <- max(bar_x0, bar_x1)
  tick_x_end_npc <- (bar_right_in + tick_len_in) / pin[1]
  tick_x_end <- graphics::grconvertX(tick_x_end_npc, from = "npc", to = "user")
  txt_x_npc <- (bar_right_in + tick_len_in + char_w_in * 0.4) / pin[1]
  txt_x <- graphics::grconvertX(txt_x_npc, from = "npc", to = "user")

  for (k in seq_along(ticks)) {
    frac <- (ticks[k] - range[1]) / (range[2] - range[1])
    frac <- max(0, min(1, frac))
    tick_y_in <- bar_bot_in + frac * bar_h_in
    tick_y <- graphics::grconvertY(tick_y_in / pin[2], from = "npc", to = "user")

    # Tick mark
    graphics::lines(c(tick_x_start, tick_x_end), c(tick_y, tick_y), col = "gray30", lwd = 0.8)
    # Tick label
    graphics::text(txt_x, tick_y, tick_labels[k], adj = c(0, 0.5), cex = cex * 0.85, col = "gray20")
  }

  cur_y_in <- bar_bot_in - char_h_in * 0.8

  # If node items are provided, draw separator and nodes
  if (n_nodes > 0) {
    sep_x0 <- graphics::grconvertX(cur_x_in / pin[1], from = "npc", to = "user")
    sep_x1 <- graphics::grconvertX((x_right_npc * pin[1] - pad_x_in) / pin[1], from = "npc", to = "user")
    sep_y <- graphics::grconvertY(cur_y_in / pin[2], from = "npc", to = "user")
    graphics::lines(c(sep_x0, sep_x1), c(sep_y, sep_y), col = "gray80", lwd = 0.8)

    cur_y_in <- cur_y_in - char_h_in * 0.4
    sym_x_npc <- (cur_x_in + char_w_in * 0.8) / pin[1]
    sym_x <- graphics::grconvertX(sym_x_npc, from = "npc", to = "user")
    node_txt_x_npc <- (cur_x_in + char_w_in * 2.0) / pin[1]
    node_txt_x <- graphics::grconvertX(node_txt_x_npc, from = "npc", to = "user")

    for (i in seq_len(n_nodes)) {
      row_y_in <- cur_y_in - (i - 0.5) * node_row_h_in
      r_y <- graphics::grconvertY(row_y_in / pin[2], from = "npc", to = "user")

      p_lwd <- if (node_items$pch[i] %in% c(8, 21, 22, 23, 24, 25)) 1.8 else 1.0
      graphics::points(sym_x, r_y, pch = node_items$pch[i], col = node_items$col[i],
                       bg = node_items$bg[i], cex = node_items$cex[i], lwd = p_lwd)
      graphics::text(node_txt_x, r_y, node_items$label[i], adj = c(0, 0.5), cex = cex * 0.9, col = "gray15")
    }
  }

  invisible()
}

# Internal helper: Plot topological skeleton and network nodes (with adjustable line width lwd)
.plot_root_skeleton <- function(x,
                                background = "image",
                                nodes = TRUE,
                                cex_nodes = 1.0,
                                lwd = 2.0,
                                col_skeleton = NULL,
                                legend_pos = "auto",
                                cex_legend = 0.75,
                                main = "Root Topological Skeleton",
                                axes = TRUE,
                                grid = FALSE,
                                ...) {
  skel <- as.array(x$skeleton)
  w <- nrow(skel)
  h <- ncol(skel)

  has_bg_img <- (background == "image" && !is.null(x$original_image))

  if (has_bg_img) {
    plot(x$original_image, axes = axes, main = main, ...)
  } else {
    bg_col <- if (background == "black") "black" else "white"
    plot(1, 1, type = "n", xlim = c(1, w), ylim = c(h, 1),
         xlab = if (isTRUE(axes)) paste0("X (", x$unit, ")") else "",
         ylab = if (isTRUE(axes)) paste0("Depth (", x$unit, ")") else "",
         axes = axes,
         main = main, asp = 1, ...)
    graphics::rect(0, 0, w + 1, h + 1, col = bg_col, border = NA)
    if (isTRUE(grid) && isTRUE(axes)) {
      grid_col <- if (background == "black") "gray30" else "gray90"
      graphics::grid(col = grid_col, lty = 1)
    }
  }

  skel_col <- if (!is.null(col_skeleton)) {
    col_skeleton
  } else if (has_bg_img) {
    grDevices::adjustcolor("#00FFFF", alpha.f = 0.90)
  } else {
    "#333333"
  }

  # Draw skeleton paths with controllable line width (lwd)
  if (!is.null(x$paths) && length(x$paths) > 0) {
    for (p in x$paths) {
      if (!is.null(p) && nrow(p) > 1) {
        graphics::lines(p[, 1], p[, 2], col = skel_col, lwd = lwd)
      } else if (!is.null(p) && nrow(p) == 1) {
        graphics::points(p[1, 1], p[1, 2], col = skel_col, pch = 16, cex = lwd * 0.4)
      }
    }
  } else {
    # Fallback to raster image if paths are unavailable
    fg <- which(skel > 0, arr.ind = TRUE)
    if (nrow(fg) > 0) {
      col_skel <- matrix("#00000000", nrow = h, ncol = w)
      col_skel[cbind(fg[, 2], fg[, 1])] <- skel_col
      graphics::rasterImage(col_skel, 0.5, h + 0.5, w + 0.5, 0.5, interpolate = FALSE)
    }
  }

  # Overlay topological nodes: junctions (orange), tips (green), crown (star)
  nds <- x$nodes
  if (isTRUE(nodes) && !is.null(nds) && nrow(nds) > 0) {
    juncs <- nds[nds$type == "junction", ]
    tips <- nds[nds$type == "tip", ]
    crowns <- nds[nds$type == "crown", ]

    if (nrow(juncs) > 0) {
      graphics::points(juncs$x, juncs$y, col = "#D55E00", pch = 19, cex = cex_nodes)
    }
    if (nrow(tips) > 0) {
      graphics::points(tips$x, tips$y, col = "#009E73", pch = 15, cex = cex_nodes)
    }
    if (nrow(crowns) > 0) {
      graphics::points(crowns$x, crowns$y, col = "#CC79A7", pch = 8, cex = cex_nodes * 2.2, lwd = 3)
    }
  }

  # Legend for skeleton plot:
  if (!isFALSE(legend_pos) && !is.null(legend_pos)) {
    avoid <- if (!is.null(nds) && nrow(nds) > 0) nds[, c("x", "y")] else NULL
    if (isTRUE(nodes)) {
      .draw_root_legend(
        pos = legend_pos,
        labels = c("Skeleton Map", "Crown", "Branch Point", "Root Tip"),
        cols = c(skel_col, "#CC79A7", "#D55E00", "#009E73"),
        types = c("line", "point", "point", "point"),
        pchs = c(NA, 8, 19, 15),
        lwds = c(lwd, NA, NA, NA),
        pt_cexs = c(1.0, 1.8 * cex_nodes, 1.2 * cex_nodes, 1.2 * cex_nodes),
        cex = cex_legend,
        img_w = w,
        img_h = h,
        avoid_pts = avoid
      )
    } else {
      .draw_root_legend(
        pos = legend_pos,
        labels = c("Skeleton Map"),
        cols = c(skel_col),
        types = c("line"),
        pchs = c(NA),
        lwds = c(lwd),
        pt_cexs = c(1.0),
        cex = cex_legend,
        img_w = w,
        img_h = h,
        avoid_pts = avoid
      )
    }
  }
}

# Internal helper: Plot root classification (Primary vs Secondary vs Tertiary+)
.plot_root_hierarchy <- function(x, background = "image", nodes = TRUE, col = NULL, cex_nodes = 1.0, legend_pos = "auto", cex_legend = 0.75, main = "Root Classification", axes = TRUE, grid = FALSE, ...) {
  skel <- as.array(x$skeleton)
  w <- nrow(skel)
  h <- ncol(skel)

  has_bg_img <- (background == "image" && !is.null(x$original_image))

  if (has_bg_img) {
    plot(x$original_image, axes = axes, main = main, ...)
  } else {
    bg_col <- if (background == "black") "black" else "white"
    plot(1, 1, type = "n", xlim = c(1, w), ylim = c(h, 1),
         xlab = if (isTRUE(axes)) paste0("X (", x$unit, ")") else "",
         ylab = if (isTRUE(axes)) paste0("Depth (", x$unit, ")") else "",
         axes = axes,
         main = main, asp = 1, ...)
    graphics::rect(0, 0, w + 1, h + 1, col = bg_col, border = NA)
    if (isTRUE(grid) && isTRUE(axes)) {
      grid_col <- if (background == "black") "gray30" else "gray90"
      graphics::grid(col = grid_col, lty = 1)
    }
  }

  orders <- if (!is.null(x$roots$order)) x$roots$order else rep(2, length(x$paths))
  has_tertiary <- any(orders >= 3)

  col_prim <- "#0072B2"
  col_sec <- "#E69F00"
  col_tert <- "#CC79A7"

  if (!is.null(col)) {
    if (length(col) >= 3) {
      col_prim <- col[1]
      col_sec <- col[2]
      col_tert <- col[3]
    } else if (length(col) == 2) {
      col_prim <- col[1]
      col_sec <- col[2]
    }
  }

  # Draw curved root paths by order: Tertiary+ first, Secondary, Primary on top
  if (!is.null(x$paths) && length(x$paths) > 0) {
    # 1. Tertiary+ roots first (order >= 3)
    for (i in which(orders >= 3)) {
      p <- x$paths[[i]]
      if (!is.null(p) && nrow(p) > 1) {
        graphics::lines(p[, 1], p[, 2], col = col_tert, lwd = 1.8)
      }
    }
    # 2. Secondary roots (order == 2)
    for (i in which(orders == 2)) {
      p <- x$paths[[i]]
      if (!is.null(p) && nrow(p) > 1) {
        graphics::lines(p[, 1], p[, 2], col = col_sec, lwd = 2.6)
      }
    }
    # 3. Primary root on top (order == 1)
    for (i in which(orders == 1)) {
      p <- x$paths[[i]]
      if (!is.null(p) && nrow(p) > 1) {
        graphics::lines(p[, 1], p[, 2], col = col_prim, lwd = 3.8)
      }
    }
  } else if (!is.null(x$roots) && nrow(x$roots) > 0) {
    roots <- x$roots
    nds <- x$nodes
    for (i in which(roots$order >= 3)) {
      r <- roots[i, ]
      n1 <- nds[nds$node_id == r$node_start, ]
      n2 <- nds[nds$node_id == r$node_end, ]
      if (nrow(n1) > 0 && nrow(n2) > 0) {
        graphics::lines(c(n1$x, n2$x), c(n1$y, n2$y), col = col_tert, lwd = 1.8)
      }
    }
    for (i in which(roots$order == 2)) {
      r <- roots[i, ]
      n1 <- nds[nds$node_id == r$node_start, ]
      n2 <- nds[nds$node_id == r$node_end, ]
      if (nrow(n1) > 0 && nrow(n2) > 0) {
        graphics::lines(c(n1$x, n2$x), c(n1$y, n2$y), col = col_sec, lwd = 2.6)
      }
    }
    for (i in which(roots$order == 1)) {
      r <- roots[i, ]
      n1 <- nds[nds$node_id == r$node_start, ]
      n2 <- nds[nds$node_id == r$node_end, ]
      if (nrow(n1) > 0 && nrow(n2) > 0) {
        graphics::lines(c(n1$x, n2$x), c(n1$y, n2$y), col = col_prim, lwd = 3.8)
      }
    }
  }

  # Overlay Crown point (and nodes if requested)
  nds <- x$nodes
  if (!is.null(nds) && nrow(nds) > 0) {
    crowns <- nds[nds$type == "crown", ]
    if (nrow(crowns) > 0) {
      graphics::points(crowns$x, crowns$y, col = "#CC79A7", pch = 8, cex = cex_nodes * 2.2, lwd = 3)
    }
    if (isTRUE(nodes)) {
      juncs <- nds[nds$type == "junction", ]
      tips <- nds[nds$type == "tip", ]
      if (nrow(juncs) > 0) {
        graphics::points(juncs$x, juncs$y, col = "#D55E00", pch = 19, cex = cex_nodes * 0.8)
      }
      if (nrow(tips) > 0) {
        graphics::points(tips$x, tips$y, col = "#009E73", pch = 15, cex = cex_nodes * 0.8)
      }
    }
  }

  # Legend for root classification
  if (!isFALSE(legend_pos) && !is.null(legend_pos)) {
    avoid <- if (!is.null(nds) && nrow(nds) > 0) nds[, c("x", "y")] else NULL
    if (has_tertiary) {
      if (isTRUE(nodes)) {
        .draw_root_legend(
          pos = legend_pos,
          labels = c("Primary Root", "Secondary Root", "Tertiary+ Root", "Crown", "Branch Point", "Root Tip"),
          cols = c(col_prim, col_sec, col_tert, "#CC79A7", "#D55E00", "#009E73"),
          types = c("line", "line", "line", "point", "point", "point"),
          pchs = c(NA, NA, NA, 8, 19, 15),
          lwds = c(3.8, 2.6, 1.8, NA, NA, NA),
          pt_cexs = c(1.0, 1.0, 1.0, 1.8 * cex_nodes, 1.0 * cex_nodes, 1.0 * cex_nodes),
          cex = cex_legend,
          img_w = w,
          img_h = h,
          avoid_pts = avoid
        )
      } else {
        .draw_root_legend(
          pos = legend_pos,
          labels = c("Primary Root", "Secondary Root", "Tertiary+ Root", "Crown"),
          cols = c(col_prim, col_sec, col_tert, "#CC79A7"),
          types = c("line", "line", "line", "point"),
          pchs = c(NA, NA, NA, 8),
          lwds = c(3.8, 2.6, 1.8, NA),
          pt_cexs = c(1.0, 1.0, 1.0, 1.8 * cex_nodes),
          cex = cex_legend,
          img_w = w,
          img_h = h,
          avoid_pts = avoid
        )
      }
    } else {
      if (isTRUE(nodes)) {
        .draw_root_legend(
          pos = legend_pos,
          labels = c("Primary Root", "Lateral Root", "Crown", "Branch Point", "Root Tip"),
          cols = c(col_prim, col_sec, "#CC79A7", "#D55E00", "#009E73"),
          types = c("line", "line", "point", "point", "point"),
          pchs = c(NA, NA, 8, 19, 15),
          lwds = c(3.8, 2.6, NA, NA, NA),
          pt_cexs = c(1.0, 1.0, 1.8 * cex_nodes, 1.0 * cex_nodes, 1.0 * cex_nodes),
          cex = cex_legend,
          img_w = w,
          img_h = h,
          avoid_pts = avoid
        )
      } else {
        .draw_root_legend(
          pos = legend_pos,
          labels = c("Primary Root", "Lateral Root", "Crown"),
          cols = c(col_prim, col_sec, "#CC79A7"),
          types = c("line", "line", "point"),
          pchs = c(NA, NA, 8),
          lwds = c(3.8, 2.6, NA),
          pt_cexs = c(1.0, 1.0, 1.8 * cex_nodes),
          cex = cex_legend,
          img_w = w,
          img_h = h,
          avoid_pts = avoid
        )
      }
    }
  }
}

# Internal helper: Plot diameter heatmap
.plot_root_diameter_map <- function(x, background = "image", nodes = TRUE, alpha = 0.65, col = NULL, cex_nodes = 1.0, legend_pos = "auto", cex_legend = 0.75, main = "Root Diameter Distribution", axes = TRUE, ...) {
  dmat <- x$diameter_map
  w <- nrow(dmat)
  h <- ncol(dmat)

  max_d <- max(dmat, na.rm = TRUE)
  if (max_d <= 0) max_d <- 1.0

  if (is.null(col)) {
    cols <- grDevices::colorRampPalette(c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"))(100)
  } else {
    cols <- grDevices::colorRampPalette(col)(100)
  }

  has_bg_img <- (background == "image" && !is.null(x$original_image))

  if (has_bg_img) {
    plot(x$original_image, axes = axes, main = main, ...)
    cols_draw <- grDevices::adjustcolor(cols, alpha.f = alpha)
  } else {
    bg_col <- if (background == "black") "black" else "white"
    plot(1, 1, type = "n", xlim = c(1, w), ylim = c(h, 1),
         xlab = if (isTRUE(axes)) paste0("X (", x$unit, ")") else "",
         ylab = if (isTRUE(axes)) paste0("Depth (", x$unit, ")") else "",
         axes = axes,
         main = main, asp = 1, ...)
    graphics::rect(0, 0, w + 1, h + 1, col = bg_col, border = NA)
    cols_draw <- cols
  }

  # Draw diameter heatmap using fast rasterImage with transparency
  col_mat <- matrix("#00000000", nrow = h, ncol = w)
  fg_idx <- which(dmat > 0, arr.ind = TRUE)
  if (nrow(fg_idx) > 0) {
    diams <- dmat[fg_idx]
    c_idx <- pmin(100, pmax(1, ceiling((diams / max_d) * 100)))
    col_mat[cbind(fg_idx[, 2], fg_idx[, 1])] <- cols_draw[c_idx]
    graphics::rasterImage(col_mat, 0.5, h + 0.5, w + 0.5, 0.5, interpolate = FALSE)
  }

  # Overlay nodes if requested
  nds <- x$nodes
  if (isTRUE(nodes)) {
    if (!is.null(nds) && nrow(nds) > 0) {
      juncs <- nds[nds$type == "junction", ]
      tips <- nds[nds$type == "tip", ]
      crowns <- nds[nds$type == "crown", ]

      if (nrow(juncs) > 0) {
        graphics::points(juncs$x, juncs$y, col = "#D55E00", pch = 19, cex = cex_nodes)
      }
      if (nrow(tips) > 0) {
        graphics::points(tips$x, tips$y, col = "#009E73", pch = 15, cex = cex_nodes)
      }
      if (nrow(crowns) > 0) {
        graphics::points(crowns$x, crowns$y, col = "#CC79A7", pch = 8, cex = cex_nodes * 2.2, lwd = 3)
      }
    }
  }

  # Continuous gradient colorbar legend with optional nodes
  if (!isFALSE(legend_pos) && !is.null(legend_pos)) {
    avoid <- if (!is.null(nds) && nrow(nds) > 0) nds[, c("x", "y")] else NULL
    node_items <- NULL
    if (isTRUE(nodes) && !is.null(nds) && nrow(nds) > 0) {
      node_items <- data.frame(
        label = c("Crown", "Branch Point", "Root Tip"),
        col = c("#CC79A7", "#D55E00", "#009E73"),
        bg = c(NA, NA, NA),
        pch = c(8, 19, 15),
        cex = c(1.8 * cex_nodes, 1.2 * cex_nodes, 1.2 * cex_nodes),
        stringsAsFactors = FALSE
      )
    }

    .draw_root_colorbar(
      pos = legend_pos,
      range = c(0, max_d),
      cols = cols,
      unit = x$unit,
      title = "Diameter",
      cex = cex_legend,
      node_items = node_items,
      avoid_pts = avoid
    )
  }
}

# Internal helper: Plot depth profile with dynamic margin calculation
.plot_root_depth <- function(x,
                             var = c("length", "volume", "area"),
                             crossings = TRUE,
                             main = NULL,
                             col = "#56B4E9",
                             border = "#0072B2",
                             col_crossings = "#D55E00",
                             grid = TRUE,
                             ...) {
  var <- match.arg(var)
  dp <- x$depth_profile
  u <- x$unit

  y_vals <- dp[[var]]
  var_label <- switch(var,
    length = paste0("Root Length (", u, ")"),
    volume = paste0("Root Volume (", u, "\u00B3)"),
    area   = paste0("Root Surface Area (", u, "\u00B2)")
  )
  if (is.null(main)) {
    main <- switch(var,
      length = "Root Depth Stratification",
      volume = "Root Volume Stratification",
      area   = "Root Surface Area Stratification"
    )
  }

  labels <- paste0(round(dp$depth_start, 1), " - ", round(dp$depth_end, 1))

  # Dynamic margin computation to prevent clipping of Y axis labels
  cur_mar <- graphics::par("mar")
  max_ch <- max(nchar(labels), default = 8)
  needed_left <- max(cur_mar[2], max_ch * 0.50)
  has_crossings <- isTRUE(crossings) && ("crossings" %in% names(dp)) && max(dp$crossings, na.rm = TRUE) > 0
  needed_top  <- if (has_crossings) max(cur_mar[3], 6) else max(cur_mar[3], 4)
  needed_bot  <- max(cur_mar[1], 4.2)
  needed_right <- max(cur_mar[4], 2.2)

  op <- graphics::par(mar = c(needed_bot, needed_left, needed_top, needed_right))
  on.exit(graphics::par(op), add = TRUE)

  max_val <- if (length(y_vals) > 0 && max(y_vals, na.rm = TRUE) > 0) max(y_vals, na.rm = TRUE) else 1

  bp <- graphics::barplot(
    rev(y_vals),
    horiz = TRUE,
    names.arg = rev(labels),
    las = 1,
    col = if (isTRUE(grid)) NA else col,
    border = if (isTRUE(grid)) NA else border,
    xlab = var_label,
    ylab = "",
    main = "",
    xlim = c(0, max_val * 1.1),
    cex.names = 0.78,
    ...
  )

  if (isTRUE(grid)) {
    graphics::grid(nx = NULL, ny = NA, col = "gray90", lty = 1)
    graphics::barplot(
      rev(y_vals),
      horiz = TRUE,
      col = col,
      border = border,
      xlim = c(0, max_val * 1.1),
      add = TRUE,
      axes = FALSE
    )
  }

  # Y axis title positioned cleanly outside rotated interval labels
  graphics::title(ylab = paste0("Depth Interval (", u, ")"), line = needed_left - 1.3, cex.lab = 0.95)

  if (has_crossings) {
    graphics::title(main = main, line = 3.6, cex.main = 1.05)
    max_cr <- max(dp$crossings, na.rm = TRUE)
    scaled_cr <- (rev(dp$crossings) / max_cr) * max_val * 0.95
    graphics::lines(scaled_cr, bp, col = col_crossings, lwd = 2.5, type = "o", pch = 19)
    graphics::axis(3, at = seq(0, max_val * 0.95, length.out = 4),
                   labels = round(seq(0, max_cr, length.out = 4)),
                   col = col_crossings, col.axis = col_crossings, line = 0.2)
    graphics::mtext("Root Crossings", side = 3, line = 1.9, col = col_crossings, cex = 0.85, font = 2)
  } else {
    graphics::title(main = main, line = 1.5, cex.main = 1.05)
  }
}

# Internal helper: Plot diameter classes
.plot_root_classes <- function(x, main = "Root Diameter Classes", grid = TRUE, ...) {
  dc <- x$diameter_classes
  u <- x$unit

  cur_mar <- graphics::par("mar")
  needed_bot <- max(cur_mar[1], 4.2)
  needed_left <- max(cur_mar[2], 4.5)
  needed_top <- max(cur_mar[3], 3.2)
  op <- graphics::par(mar = c(needed_bot, needed_left, needed_top, cur_mar[4]))
  on.exit(graphics::par(op), add = TRUE)

  n_bars <- nrow(dc)
  bar_cols <- if (n_bars == 4) {
    c("#56B4E9", "#009E73", "#E69F00", "#D55E00")
  } else {
    grDevices::colorRampPalette(c("#56B4E9", "#009E73", "#E69F00", "#D55E00"))(n_bars)
  }

  max_len <- if (length(dc$length) > 0 && max(dc$length, na.rm = TRUE) > 0) max(dc$length, na.rm = TRUE) else 1

  bp <- graphics::barplot(
    dc$length,
    names.arg = dc$class,
    col = if (isTRUE(grid)) NA else bar_cols,
    border = if (isTRUE(grid)) NA else "black",
    ylab = paste0("Total Length (", u, ")"),
    main = main,
    las = 1,
    ylim = c(0, max_len * 1.15),
    cex.names = 0.85,
    ...
  )

  if (isTRUE(grid)) {
    graphics::grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    graphics::barplot(
      dc$length,
      col = bar_cols,
      ylim = c(0, max_len * 1.15),
      add = TRUE,
      axes = FALSE
    )
  }

  graphics::text(bp, dc$length + max_len * 0.04,
                 labels = paste0(round(dc$length, 1), " ", u),
                 cex = 0.8)
}

# Internal helper: Plot branching angles distribution
.plot_root_angles <- function(x, main = "Lateral Root Branching Angles", legend_pos = "auto", cex_legend = 0.75, grid = TRUE, ...) {
  angles <- stats::na.omit(x$roots$branching_angle)
  if (length(angles) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = main)
    graphics::text(1, 1, "No lateral branching angles recorded", cex = 1.1)
    return()
  }

  cur_mar <- graphics::par("mar")
  op <- graphics::par(mar = c(max(cur_mar[1], 4.2), max(cur_mar[2], 4.5), max(cur_mar[3], 3.2), cur_mar[4]))
  on.exit(graphics::par(op), add = TRUE)

  h_ang <- graphics::hist(angles, breaks = seq(0, 180, by = 15), plot = FALSE)
  plot(h_ang, col = NA, border = NA,
       xlab = "Insertion Angle (degrees)",
       ylab = "Frequency",
       main = main,
       xlim = c(0, 180),
       las = 1,
       ...)

  if (isTRUE(grid)) {
    graphics::grid(col = "gray90", lty = 1)
  }

  plot(h_ang, col = "#009E73", border = "white", add = TRUE)
  graphics::abline(v = mean(angles), col = "#D55E00", lwd = 2.5, lty = 2)

  leg_pos <- if (legend_pos == "auto") "topright" else legend_pos
  if (!isFALSE(leg_pos) && !is.null(leg_pos)) {
    .draw_root_legend(
      pos = leg_pos,
      labels = paste0("Mean: ", round(mean(angles), 1), "\u00B0"),
      cols = "#D55E00",
      types = "line",
      lwds = 2.5,
      cex = cex_legend
    )
  }
}

# Internal helper: Plot spatial architectural profile and convex hull envelope
.plot_root_profile <- function(x,
                               main = "Root System Spatial Profile & Envelope",
                               col_envelope = "#56B4E933",
                               col_skeleton = "#0072B2",
                               axes = TRUE,
                               legend_pos = "auto",
                               cex_legend = 0.75,
                               grid = TRUE,
                               ...) {
  u <- x$unit
  paths <- x$paths
  nodes <- x$nodes
  s <- x$summary
  ps <- x$pixel_size

  all_pts <- do.call(rbind, paths)
  if (is.null(all_pts) || nrow(all_pts) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = main)
    graphics::text(1, 1, "No root skeleton coordinates available", cex = 1.1)
    return()
  }

  xs <- all_pts[, 1] * ps
  ys <- all_pts[, 2] * ps

  crown_node <- nodes[nodes$type %in% c("crown", 3), ]
  crown_x <- if (nrow(crown_node) > 0) crown_node$x[1] * ps else mean(xs)
  crown_y <- if (nrow(crown_node) > 0) crown_node$y[1] * ps else min(ys)

  xlim <- range(xs) + c(-1, 1) * diff(range(xs)) * 0.08
  ylim <- c(max(ys) * 1.05, max(0, min(ys) * 0.95))

  cur_mar <- graphics::par("mar")
  op <- graphics::par(mar = c(max(cur_mar[1], 4.5), max(cur_mar[2], 4.8), max(cur_mar[3], 3.5), max(cur_mar[4], 2.0)))
  on.exit(graphics::par(op), add = TRUE)

  plot(xs, ys, type = "n", xlim = xlim, ylim = ylim,
       xlab = paste0("Horizontal Span (", u, ")"),
       ylab = paste0("Depth below Crown (", u, ")"),
       main = main, las = 1, axes = axes, ...)
  if (isTRUE(grid) && isTRUE(axes)) graphics::grid(col = "gray90", lty = 1)

  ch_idx <- grDevices::chull(xs, ys)
  graphics::polygon(xs[ch_idx], ys[ch_idx], col = col_envelope, border = "#56B4E9", lwd = 1.5, lty = 2)

  for (p in paths) {
    graphics::lines(p[, 1] * ps, p[, 2] * ps, col = col_skeleton, lwd = 2)
  }

  cx <- if (!is.null(s$centroid_x) && !is.na(s$centroid_x)) s$centroid_x else mean(xs)
  cy <- if (!is.null(s$centroid_y) && !is.na(s$centroid_y)) s$centroid_y else mean(ys)
  if (!is.null(cx) && !is.null(cy)) {
    graphics::points(cx, cy, pch = 23, bg = "#E69F00", col = "black", cex = 1.8, lwd = 1.5)
  }

  graphics::points(crown_x, crown_y, pch = 21, bg = "#D55E00", col = "black", cex = 1.8, lwd = 1.5)

  if (!isFALSE(legend_pos) && !is.null(legend_pos)) {
    avoid <- if (!is.null(cx)) rbind(c(crown_x, crown_y), c(cx, cy)) else matrix(c(crown_x, crown_y), nrow = 1)
    .draw_root_legend(
      pos = legend_pos,
      labels = c("Crown Node", "Centroid", "Root Skeleton", "Convex Envelope"),
      cols = c("black", "black", col_skeleton, "#56B4E9"),
      fills = c("#D55E00", "#E69F00", NA, col_envelope),
      types = c("point", "point", "line", "polygon"),
      pchs = c(21, 23, NA, NA),
      lwds = c(NA, NA, 2.0, 1.5),
      pt_cexs = c(1.3, 1.3, NA, NA),
      cex = cex_legend,
      avoid_pts = avoid
    )
  }

  metrics_txt <- paste0(
    "Max Width: ", round(s$max_width, 1), " ", u, "\n",
    "Max Depth: ", round(s$max_depth, 1), " ", u, "\n",
    "W:D Ratio: ", round(s$width_depth_ratio, 2), "\n",
    "Solidity: ", round(s$solidity, 3)
  )
  graphics::text(xlim[1] + diff(xlim)*0.02, ylim[2] + diff(ylim)*0.02,
                 labels = metrics_txt, adj = c(0, 1), cex = 0.78, font = 3, col = "gray30")
}

# Internal helper: Plot gravitropic branching angle vs insertion depth
.plot_root_topology <- function(x,
                                main = "Lateral Branching Angles vs Insertion Depth",
                                axes = TRUE,
                                legend_pos = "auto",
                                cex_legend = 0.75,
                                grid = TRUE,
                                ...) {
  r <- x$roots
  u <- x$unit
  s <- x$summary
  sub_r <- r[!is.na(r$branching_angle) & !is.na(r$insertion_depth), ]
  if (nrow(sub_r) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = main)
    graphics::text(1, 1, "No lateral branching angles and insertion depths recorded", cex = 1.1)
    return()
  }

  cur_mar <- graphics::par("mar")
  op <- graphics::par(mar = c(max(cur_mar[1], 4.5), max(cur_mar[2], 4.8), max(cur_mar[3], 3.5), max(cur_mar[4], 2.0)))
  on.exit(graphics::par(op), add = TRUE)

  cols <- ifelse(sub_r$order == 2, "#009E73", "#E69F00")
  cexs <- 1.0 + (sub_r$length / max(sub_r$length, na.rm = TRUE)) * 1.5

  max_d <- if (!is.null(s$max_depth) && !is.na(s$max_depth)) s$max_depth else max(sub_r$insertion_depth)
  ylim <- c(max_d * 1.05, 0)

  plot(sub_r$branching_angle, sub_r$insertion_depth,
       type = "n",
       xlim = c(0, 180), ylim = ylim,
       xlab = "Branching Angle (degrees)",
       ylab = paste0("Insertion Depth (", u, ")"),
       main = main, las = 1, axes = axes, ...)
  if (isTRUE(grid) && isTRUE(axes)) {
    graphics::grid(col = "gray90", lty = 1)
  }
  if (isTRUE(axes)) {
    graphics::abline(v = 90, lty = 3, col = "gray50")
  }

  graphics::points(sub_r$branching_angle, sub_r$insertion_depth,
                   pch = 21, bg = cols, col = "black", cex = cexs)

  if (nrow(sub_r) >= 4) {
    fit <- tryCatch(stats::loess(insertion_depth ~ branching_angle, data = sub_r, span = 0.9), error = function(e) NULL)
    if (!is.null(fit)) {
      grid_x <- seq(min(sub_r$branching_angle), max(sub_r$branching_angle), length.out = 100)
      pred_y <- stats::predict(fit, newdata = data.frame(branching_angle = grid_x))
      graphics::lines(grid_x, pred_y, col = "#D55E00", lwd = 2.5)
    }
  }

  if (!isFALSE(legend_pos) && !is.null(legend_pos)) {
    avoid_pts <- cbind(sub_r$branching_angle, sub_r$insertion_depth)
    .draw_root_legend(
      pos = legend_pos,
      labels = c("Secondary Root (2nd)", "Tertiary+ Root (3rd+)", "Orthogonal (90\u00B0)"),
      cols = c("black", "black", "gray50"),
      fills = c("#009E73", "#E69F00", NA),
      types = c("point", "point", "line"),
      pchs = c(21, 21, NA),
      lwds = c(NA, NA, 2.0),
      pt_cexs = c(1.3, 1.3, NA),
      cex = cex_legend,
      avoid_pts = avoid_pts
    )
  }
}

# Internal helper: Plot continuous diameter distribution and cumulative length
.plot_root_diameter_dist <- function(x,
                                     main = "Continuous Root Diameter Distribution",
                                     axes = TRUE,
                                     legend_pos = "auto",
                                     cex_legend = 0.75,
                                     grid = TRUE,
                                     ...) {
  r <- x$roots
  u <- x$unit
  diams <- r$avg_diameter[!is.na(r$avg_diameter) & r$avg_diameter > 0]
  lengths <- r$length[!is.na(r$avg_diameter) & r$avg_diameter > 0]

  if (length(diams) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = main)
    graphics::text(1, 1, "No root diameter measurements available", cex = 1.1)
    return()
  }

  cur_mar <- graphics::par("mar")
  op <- graphics::par(mar = c(max(cur_mar[1], 4.5), max(cur_mar[2], 4.8), max(cur_mar[3], 3.5), max(cur_mar[4], 4.8)))
  on.exit(graphics::par(op), add = TRUE)

  h <- graphics::hist(diams, breaks = 15, plot = FALSE)
  max_dens <- max(h$density)
  if (length(diams) >= 3) {
    dens <- stats::density(diams)
    max_dens <- max(max_dens, max(dens$y))
  }
  ylim_max <- max_dens * 1.30

  # 1. Initialize canvas without drawing bars
  plot(h, freq = FALSE, col = NA, border = NA,
       xlab = paste0("Root Diameter (", u, ")"),
       ylab = "Probability Density",
       main = main, las = 1, ylim = c(0, ylim_max), axes = axes, ...)

  # 2. Draw background grid
  if (isTRUE(grid) && isTRUE(axes)) {
    graphics::grid(col = "gray90", lty = 1)
  }

  # 3. Draw histogram bars over grid
  plot(h, freq = FALSE, col = "#56B4E966", border = "#0072B2", add = TRUE)

  # 4. Density curve over bars
  if (length(diams) >= 3) {
    graphics::lines(dens, col = "#0072B2", lwd = 2.5)
  }

  mean_d <- x$summary$avg_diameter
  med_d <- x$summary$median_diameter
  graphics::abline(v = mean_d, col = "#D55E00", lwd = 2, lty = 2)
  graphics::abline(v = med_d, col = "#009E73", lwd = 2, lty = 3)

  ord <- order(diams)
  cum_len <- cumsum(lengths[ord]) / sum(lengths) * 100
  cum_y <- (cum_len / 100) * (ylim_max * 0.95)
  graphics::lines(diams[ord], cum_y, col = "#CC79A7", lwd = 2)
  graphics::axis(4, at = seq(0, ylim_max * 0.95, length.out = 5),
                 labels = paste0(seq(0, 100, by = 25), "%"),
                 col = "#CC79A7", col.axis = "#CC79A7", las = 1)
  graphics::mtext("Cumulative Length (%)", side = 4, line = 2.8, col = "#CC79A7", cex = 0.85)

  # Default to topright as requested, unless explicitly specified
  leg_pos <- if (legend_pos == "auto") "topright" else legend_pos
  if (!isFALSE(leg_pos) && !is.null(leg_pos)) {
    .draw_root_legend(
      pos = leg_pos,
      labels = c(paste0("Mean: ", round(mean_d, 2), " ", u),
                 paste0("Median: ", round(med_d, 2), " ", u),
                 "KDE Density",
                 "Cumulative %"),
      cols = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"),
      types = c("line", "line", "line", "line"),
      lwds = c(2.0, 2.0, 2.5, 2.0),
      cex = cex_legend
    )
  }
}

# Internal helper: Plot morphometric breakdown across root branching orders
.plot_root_orders <- function(x,
                              main = "Root Orders Morphometric Breakdown",
                              grid = TRUE,
                              ...) {
  r <- x$roots
  u <- x$unit
  if (nrow(r) == 0) return()

  op <- graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.8, 3.5, 1.5))
  on.exit(graphics::par(op), add = TRUE)

  ord_split <- split(r, r$order)
  orders <- names(ord_split)
  ord_names <- sapply(orders, function(o) {
    if (o == "1") "Primary"
    else if (o == "2") "Secondary"
    else paste0("Tertiary (", o, ")")
  })
  cols <- c("#D55E00", "#009E73", "#E69F00", "#56B4E9")[seq_along(orders)]

  tot_lens <- sapply(ord_split, function(df) sum(df$length, na.rm = TRUE))
  max_l <- max(tot_lens, na.rm = TRUE) * 1.15
  bp1 <- graphics::barplot(tot_lens, names.arg = ord_names,
                          col = if (isTRUE(grid)) NA else cols,
                          border = if (isTRUE(grid)) NA else "black",
                          ylim = c(0, max_l),
                          ylab = paste0("Total Length (", u, ")"),
                          main = "Length by Order", las = 1, cex.names = 0.85, ...)
  if (isTRUE(grid)) {
    graphics::grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    graphics::barplot(tot_lens, col = cols, ylim = c(0, max_l), add = TRUE, axes = FALSE)
  }
  graphics::text(bp1, tot_lens + max(tot_lens)*0.04, labels = paste0(round(tot_lens, 1)), cex = 0.8)

  avg_diams <- sapply(ord_split, function(df) stats::weighted.mean(df$avg_diameter, df$length, na.rm = TRUE))
  max_d <- max(avg_diams, na.rm = TRUE) * 1.15
  bp2 <- graphics::barplot(avg_diams, names.arg = ord_names,
                          col = if (isTRUE(grid)) NA else cols,
                          border = if (isTRUE(grid)) NA else "black",
                          ylim = c(0, max_d),
                          ylab = paste0("Mean Diameter (", u, ")"),
                          main = "Diameter by Order", las = 1, cex.names = 0.85, ...)
  if (isTRUE(grid)) {
    graphics::grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    graphics::barplot(avg_diams, col = cols, ylim = c(0, max_d), add = TRUE, axes = FALSE)
  }
  graphics::text(bp2, avg_diams + max(avg_diams)*0.04, labels = paste0(round(avg_diams, 2)), cex = 0.8)

  counts <- sapply(ord_split, nrow)
  max_c <- max(counts, na.rm = TRUE) * 1.15
  bp3 <- graphics::barplot(counts, names.arg = ord_names,
                          col = if (isTRUE(grid)) NA else cols,
                          border = if (isTRUE(grid)) NA else "black",
                          ylim = c(0, max_c),
                          ylab = "Branch Count (n)",
                          main = "Branch Count by Order", las = 1, cex.names = 0.85, ...)
  if (isTRUE(grid)) {
    graphics::grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    graphics::barplot(counts, col = cols, ylim = c(0, max_c), add = TRUE, axes = FALSE)
  }
  graphics::text(bp3, counts + max(counts)*0.04, labels = counts, cex = 0.8)
}
