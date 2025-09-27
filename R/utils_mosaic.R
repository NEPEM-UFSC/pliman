validate_and_replicate <- function(argument, created_shapes, verbose = TRUE) {
  if ((length(argument) != length(created_shapes)) & verbose) {
    cli::cli_warn(c(
      "!" = "`{.arg {deparse(substitute(argument))}}` must have length 1 or {.val {length(created_shapes)}} (the number of drawn polygons)."
    ))
  }

  if (length(argument) == 1 & length(created_shapes) != 1) {
    argument <- rep(argument, length(created_shapes))
  }
  return(argument)
}
validate_and_replicate2 <- function(argument, created_shapes, verbose = TRUE) {
  if ((!is.null(argument) && (length(argument) != nrow(created_shapes))) && verbose) {
    cli::cli_warn(c(
      "!" = "`{.arg {deparse(substitute(argument))}}` must have length 1 or {.val {nrow(created_shapes)}} (the number of drawn polygons)."
    ))
  }

  if (length(argument) == 1 & nrow(created_shapes) != 1) {
    argument <- rep(argument, nrow(created_shapes))
  }
  return(argument)
}
sf_to_polygon <- function(shps) {
  if(inherits(shps, "list")){
    shps <- do.call(rbind, shps)
  }
  classes <- sapply(lapply(sf::st_geometry(shps$geometry), class), function(x){x[2]})
  shps[classes %in% c("POINT", "LINESTRING"), ] <-
    shps[classes %in% c("POINT", "LINESTRING"), ] |>
    sf::st_buffer(0.0000001) |>
    sf::st_cast("POLYGON") |>
    sf::st_simplify(preserveTopology = TRUE)
  return(shps)
}

find_aggrfact <- function(mosaic, max_pixels = 1000000){
  compute_downsample <- function(nr, nc, n) {
    if (n == 0) {
      invisible(nr * nc)
    } else if (n == 1) {
      invisible(ceiling(nr/2) * ceiling(nc/2))
    } else if (n > 1) {
      invisible(ceiling(nr/(n+1)) * ceiling(nc/(n+1)))
    } else {
      cli::cli_abort("{.arg n} must be a non-negative integer.")

    }
  }
  nr <- nrow(mosaic)
  nc <- ncol(mosaic)
  npixel <- nr * nc
  possible_downsamples <- 0:20
  possible_npix <- sapply(possible_downsamples, function(x){
    compute_downsample(nr, nc, x)
  })
  downsample <- which.min(abs(possible_npix - max_pixels))
  downsample <- ifelse(downsample == 1, 0, downsample)
  return(downsample)
}
compute_measures_mosaic <- function(contour){
  lw <- help_lw(contour)
  cdist <- help_centdist(contour)
  data.frame(area = help_area(contour),
             perimeter = sum(help_distpts(contour)),
             length = lw[[1]],
             width = lw[[2]],
             diam_min = min(cdist) * 2,
             diam_mean = mean(cdist) * 2,
             diam_max = max(cdist) * 2)

}
compute_dists <- function(subset_coords, direction = c("horizontal", "vertical")){
  optdirec <- c("horizontal", "vertical")
  optdirec <- pmatch(direction[[1]], optdirec)
  n <- nrow(subset_coords)
  subset_coords <- subset_coords |> dplyr::select(x, y) |> as.data.frame()
  nearest <- order(subset_coords[, optdirec])
  subset_distances <- numeric(n - 1)
  for (j in 1:(n - 1)) {
    x1 <- subset_coords[nearest[j], 1]
    y1 <- subset_coords[nearest[j], 2]
    x2 <- subset_coords[nearest[j+1], 1]
    y2 <- subset_coords[nearest[j+1], 2]
    distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    subset_distances[j] <- distance
  }
  subset_distances
}

linear_iterpolation <- function(mosaic, points, method = "loess"){
  if(inherits(points, "list")){
    points <- do.call(rbind, points)
  }
  xy <- sf::st_coordinates(points)[, 1:2]
  vals <- terra::values(mosaic)[terra::cellFromXY(mosaic, xy), ]
  vals <- data.frame(cbind(xy, vals))
  names(vals) <- c("x", "y", "z")
  newdata <- as.data.frame(terra::xyFromCell(mosaic, 1:terra::ncell(mosaic)))
  new_ras <-
    terra::rast(
      lapply(3:ncol(vals), function(i){
        if(method == "loess"){
          mod <- loess(vals[, i] ~ x + y, data = vals)
        } else{
          mod <- lm(vals[, i] ~ x + y, data = vals)
        }
        terra::rast(matrix(predict(mod, newdata = newdata),
                           nrow = nrow(mosaic),
                           ncol = ncol(mosaic),
                           byrow = TRUE))
      })
    )
  terra::crs(new_ras) <- terra::crs(mosaic)
  terra::ext(new_ras) <- terra::ext(mosaic)
  terra::resample(new_ras, mosaic)
}
idw_interpolation <- function(mosaic, points){
  downsample <- find_aggrfact(mosaic, max_pixels = 200000)
  if(downsample > 0){
    magg <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
  } else{
    magg <- mosaic
  }
  if(inherits(points, "list")){
    points <- do.call(rbind, points)
  }
  xy <- sf::st_coordinates(points)[, 1:2]
  vals <- terra::values(magg)[terra::cellFromXY(magg, xy), ]
  vals <- data.frame(cbind(xy, vals))

  xy_grid <- terra::xyFromCell(magg, 1:terra::ncell(magg))
  newx <- seq(min(xy_grid[,1]), max(xy_grid[,1]), length.out = 1000)
  newy <- seq(min(xy_grid[,2]), max(xy_grid[,2]), length.out = 1000)

  new_ras <-
    terra::rast(
      lapply(3:ncol(vals), function(i){
        interp <- idw_interpolation_cpp(vals[, 1], vals[, 2], vals[, i], xy_grid[, 1], xy_grid[, 2])
        ra3 <-
          terra::rast(matrix(interp,
                             nrow = nrow(magg),
                             ncol = ncol(magg),
                             byrow = TRUE))
      })
    )

  terra::crs(new_ras) <- terra::crs(mosaic)
  terra::ext(new_ras) <- terra::ext(mosaic)
  terra::resample(new_ras, mosaic)
}

# Helper function to check and align DSM with mosaic
align_dsm <- function(dsm, mosaic) {
  # Check if extent, resolution, and CRS match
  if (!terra::ext(dsm) == terra::ext(mosaic) ||
      !all(terra::res(dsm) == terra::res(mosaic)) ||
      !terra::crs(dsm) == terra::crs(mosaic)) {
    cli::cli_alert_info("Adjusting DSM to match mosaic extent, resolution, and CRS.")
    # Align DSM to match the mosaic properties
    dsm <- terra::resample(dsm, mosaic, method = "lanczos")
  }
  return(dsm)
}

#' Mosaic interpolation
#'
#' Performs the interpolation of points from a raster object.
#'
#' @param mosaic An `SpatRaster` object
#' @param points An `sf` object with the points for x and y coordinates, usually
#'   obtained with [shapefile_build()]. Alternatively, an external shapefile
#'   imported with [shapefile_input()] containing the x and y coordinates can be
#'   used. The function will handle most used shapefile formats (eg.,
#'   .shp, .rds) and convert the imported shapefile to an sf object.
#' @param method One of "bilinear" (default), "loess" (local regression) or
#'   "idw" (Inverse Distance Weighting).
#' @importFrom stats loess
#'
#' @return An `SpatRaster` object with the same extent and crs from `mosaic`
#' @export
#'
mosaic_interpolate <- function(mosaic, points, method = c("bilinear", "loess", "idw")){
  if(terra::crs(points) != terra::crs(mosaic)){
    terra::crs(points) <- terra::crs(mosaic)
  }
  if(!method[[1]] %in% c("bilinear", "idw", "loess")){
    cli::cli_abort(c(
      "x" = "`method` must be one of {.val bilinear}, {.val loess}, or {.val idw}."
    ))

  }
  if(method[[1]]  %in%  c("bilinear", "loess")){
    linear_iterpolation(mosaic, points, method = method[[1]])
  } else{
    idw_interpolation(mosaic, points)
  }
}


#' Analyze a mosaic of remote sensing data
#'
#' This function analyzes a mosaic of remote sensing data (UVAs or satellite
#' imagery), extracting information from specified regions of interest (ROIs)
#' defined in a shapefile or interactively drawn on the mosaic. It allows
#' counting and measuring individuals (eg., plants), computing canopy coverage,
#' and statistical summaries (eg., mean, coefficient of variation) for
#' vegetation indices (eg, NDVI) at a block, plot, individual levels or even
#' extract the raw results at pixel level.
#'
#' @details
#' Since multiple blocks can be analyzed, the length of arguments `grid`,
#' `nrow`, `ncol`, `buffer_edge`, , `buffer_col`, `buffer_row`, `segment_plot`,
#' `segment_i, ndividuals`, `includ_if`, `threshold`, `segment_index`, `invert`,
#' `filter`, `threshold`, `lower_size`, `upper_size`, `watershed`, and
#' `lower_noise`, can be either an scalar (the same argument applied to all the
#' drawn blocks), or a vector with the same length as the number of drawn. In
#' the last, each block can be analyzed with different arguments.
#'
#' When `segment_individuals = TRUE` is enabled, individuals are included within
#' each plot based on the `include_if` argument. The default value
#' (`'centroid'`) includes an object in a given plot if the centroid of that
#' object is within the plot. This makes the inclusion mutually exclusive (i.e.,
#' an individual is included in only one plot). If `'covered'` is selected,
#' objects are included only if their entire area is covered by the plot. On the
#' other hand, selecting `overlap` is the complement of `covered`; in other
#' words, objects that overlap the plot boundary are included. Finally, when
#' `intersect` is chosen, objects that intersect the plot boundary are included.
#' This makes the inclusion ambiguous (i.e., an object can be included in more
#' than one plot).

#'
#' @inheritParams mosaic_view
#' @inheritParams analyze_objects
#' @inheritParams image_binary
#' @inheritParams plot_id
#' @param r,g,b,re,nir,swir,tir The red, green, blue, red-edge,  near-infrared,
#'   shortwave Infrared, and thermal infrared bands of the image, respectively.
#'   By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
#'   multispectral image is provided up to seven bands can be used to compute
#'   built-in indexes. There are no limitation of band numbers if the index is
#'   computed using the band name.
#' @param crop_to_shape_ext Crop the mosaic to the extension of shapefile?
#'   Defaults to `TRUE`. This allows for a faster index computation when the
#'   region of the built shapefile is much smaller than the entire mosaic
#'   extension.
#' @param grid Logical, indicating whether to use a grid for segmentation
#'   (default: TRUE).
#' @param nrow Number of rows for the grid (default: 1).
#' @param ncol Number of columns for the grid (default: 1).
#' @param plot_width,plot_height The width and height of the plot shape (in the
#'   mosaic unit). It is mutually exclusiv with `buffer_col` and `buffer_row`.
#' @param indexes An optional `SpatRaster` object with the image indexes,
#'   computed with [mosaic_index()].
#' @param shapefile An optional shapefile containing regions of interest (ROIs)
#'   for analysis.
#' @param basemap An optional basemap generated with [mosaic_view()].
#' @param build_shapefile Logical, indicating whether to interactively draw ROIs
#'   if the shapefile is `NULL` (default: TRUE).
#' @param check_shapefile Logical, indicating whether to validate the shapefile
#'   with an interactive map view (default: TRUE). This enables live editing of
#'   the drawn shapefile by deleting or changing the drawn grids.
#' @param buffer_edge Width of the buffer around the shapefile (default: 5).
#' @param buffer_col,buffer_row Buffering factor for the columns and rows,
#'   respectively, of each individual plot's side. A value between 0 and 0.5
#'   where 0 means no buffering and 0.5 means complete buffering (default: 0). A
#'   value of 0.25 will buffer the plot by 25% on each side.
#' @param segment_plot Logical, indicating whether to segment plots (default:
#'   FALSE). If `TRUE`, the `segment_index` will be computed, and pixels with
#'   values below the `threshold` will be selected.
#' @param segment_individuals Logical, indicating whether to segment individuals
#'   within plots (default: FALSE). If `TRUE`, the `segment_index` will be
#'   computed, and pixels with values below the `threshold` will be selected, and
#'   a watershed-based segmentation will be performed.
#' @param segment_pick When `segment_plot` or `segment_individuals` are `TRUE`,
#'   `segment_pick` allows segmenting background (eg., soil) and foreground
#'   (eg., plants) interactively by picking samples from background and
#'   foreground using [mosaic_segment_pick()]
#' @param mask An optional mask (SpatRaster) to mask the mosaic.
#' @param dsm A SpatRaster object representing the digital surface model. Must
#'   be a single-layer raster. If a DSM is informed, a mask will be derived from
#'   it using [mosaic_chm_mask()].
#' @param dsm_lower A numeric value specifying the lower height threshold. All
#'   heights greater than this value are retained.
#' @param dsm_upper An optional numeric value specifying the upper height
#'   threshold. If provided, only heights between lower and upper are retained.
#' @param dsm_window_size An integer (meters) specifying the window size (rows
#'   and columns, respectively) for creating a DTM using a moving window.
#'   Default is c(5, 5).
#' @param simplify Removes vertices in polygons to form simpler shapes. The
#'   function implementation uses the Douglas-Peucker algorithm using
#'   [sf::st_simplify()] for simplification.
#' @param map_individuals If `TRUE`, the distance between objects within plots
#'   is computed. The distance can be mapped either in the horizontal or vertical
#'   direction. The distances, coefficient of variation (CV), and mean of
#'   distances are then returned.
#' @param map_direction The direction for mapping individuals within plots.
#'   Should be one of `"horizontal"` or `"vertical"` (default).
#' @param watershed If `TRUE` (default), performs watershed-based object
#'   detection. This will detect objects even when they are touching one another.
#'   If FALSE, all pixels for each connected set of foreground pixels are set to
#'   a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. A higher value smooths out small objects.
#' @param include_if Character vector specifying the type of intersection.
#'   Defaults to "centroid" (individuals in which the centroid is included within
#'   the drawn plot will be included in that plot). Other possible values include
#'   `"covered"`, `"overlap"`, and `"intersect"`. See Details for a detailed
#'   explanation of these intersecting controls.
#' @param plot_index The index(es) to be computed for the drawn plots. Either a
#'   single vegetation index (e.g., `"GLAI"`), a vector of indexes (e.g.,
#'   `c("GLAI", "NGRDI", "HUE")`), or a custom index based on the available
#'   bands (e.g., `"(R-B)/(R+B)"`). See [pliman_indexes()] and [image_index()]
#'   for more details.
#' @param segment_index The index used for segmentation. The same rule as
#'   `plot_index`. Defaults to `NULL`
#' @param threshold By default (threshold = "Otsu"), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is provided, this value will be used as a threshold.
#' @param summarize_fun The function to compute summaries for the pixel values.
#'   Defaults to "mean," i.e., the mean value of the pixels (either at a plot- or
#'   individual-level) is returned.
#' @param summarize_quantiles quantiles to be computed when 'quantile' is on `summarize_fun`.
#' @param attribute The attribute to be shown at the plot when `plot` is `TRUE`. Defaults to the first `summary_fun` and first `segment_index`.
#' @param invert Logical, indicating whether to invert the mask. Defaults to
#'   `FALSE`, i.e., pixels with intensity greater than the threshold values are
#'   selected.
#' @param color_regions The color palette for regions (default:
#'   rev(grDevices::terrain.colors(50))).
#' @param plot Logical, indicating whether to generate plots (default: TRUE).
#' @param verbose Logical, indicating whether to display verbose output
#'   (default: TRUE).
#'
#' @return A list containing the following objects:
#' * `result_plot`: The results at a plot level.
#' *  `result_plot_summ`: The summary of results at a plot level. When
#'  `segment_individuals = TRUE`, the number of individuals, canopy coverage,
#'  and mean values of some shape statistics such as perimeter, length, width,
#'  and diameter are computed.
#' * `result_individ`: The results at an individual level.
#' * `map_plot`: An object of class `mapview` showing the plot-level results.
#' * `map_individual`: An object of class `mapview` showing the individual-level
#'   results.
#' * `shapefile`: The generated shapefile, with the drawn grids/blocks.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' url <- "https://github.com/TiagoOlivoto/images/raw/master/pliman/rice_field/rice_ex.tif"
#' mosaic <- mosaic_input(url)
#' # Draw a polygon (top left, top right, bottom right, bottom left, top left)
#' # include 8 rice lines and one column
#'res <-
#'  mosaic_analyze(mosaic,
#'                 r = 1, g = 2, b = 3,
#'                 segment_individuals = TRUE,     # segment the individuals
#'                 segment_index = "(G-B)/(G+B-R)",# index for segmentation
#'                 filter = 4,
#'                 nrow = 8,
#'                 map_individuals = TRUE)
#'# map with individual results
#'res$map_indiv
#' }

mosaic_analyze <- function(mosaic,
                           r = 3,
                           g = 2,
                           b = 1,
                           re = NA,
                           nir = NA,
                           swir = NA,
                           tir = NA,
                           crop_to_shape_ext = TRUE,
                           grid = TRUE,
                           nrow = 1,
                           ncol = 1,
                           plot_width = NULL,
                           plot_height = NULL,
                           layout = "lrtb",
                           indexes = NULL,
                           shapefile = NULL,
                           basemap = NULL,
                           build_shapefile = TRUE,
                           check_shapefile = TRUE,
                           buffer_edge = 1,
                           buffer_col = 0,
                           buffer_row = 0,
                           segment_plot = FALSE,
                           segment_individuals = FALSE,
                           segment_pick = FALSE,
                           mask = NULL,
                           dsm = NULL,
                           dsm_lower = 0.2,
                           dsm_upper = NULL,
                           dsm_window_size = c(5, 5),
                           simplify = FALSE,
                           map_individuals = FALSE,
                           map_direction = c("horizontal", "vertical"),
                           watershed = TRUE,
                           tolerance = 1,
                           extension = 1,
                           include_if = "centroid",
                           plot_index = "GLI",
                           segment_index = NULL,
                           threshold = "Otsu",
                           opening = FALSE,
                           closing = FALSE,
                           filter = FALSE,
                           erode = FALSE,
                           dilate = FALSE,
                           lower_noise = 0.15,
                           lower_size = NULL,
                           upper_size = NULL,
                           topn_lower = NULL,
                           topn_upper = NULL,
                           summarize_fun = "mean",
                           summarize_quantiles = NULL,
                           attribute = NULL,
                           invert = FALSE,
                           color_regions = rev(grDevices::terrain.colors(50)),
                           alpha = 1,
                           max_pixels = 2e6,
                           downsample = NULL,
                           quantiles = c(0, 1),
                           plot = TRUE,
                           verbose = TRUE){
  if(verbose){
    cli::cli_rule(
      left = cli::col_blue("Analyzing the mosaic"),
      right = cli::col_blue("Started on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
    )
  }
  if(!is.null(dsm)){
    dsm <- align_dsm(dsm, mosaic)
    if(verbose){
      msg <- "Creating the mask based on the digital surface model..."
      cli::cli_progress_step(
        msg = msg,
        msg_done = sub("\\.\\.\\.$", "", msg),
        msg_failed = "Mask creation failed"
      )
    }
    mask <- mosaic_chm_mask(dsm, lower = dsm_lower, upper = dsm_upper, window_size = dsm_window_size)
  }
  includeopt <- c("intersect", "covered", "overlap", "centroid")
  includeopt <- includeopt[sapply(include_if, function(x){pmatch(x, includeopt)})]
  if(is.null(plot_index) & !is.null(segment_index)){
    plot_index <- segment_index
  }
  if(!is.null(indexes)){
    if(!inherits(indexes, "SpatRaster")){
      cli::cli_abort("Object {.arg indexes} must be an object of class `SpatRaster`.")
    } else{
      plot_index <- names(indexes)
    }
  }
  if(!is.null(plot_index) & is.null(segment_index)){
    segment_index <- plot_index[[1]]
  }
  if(any(segment_individuals) | any(segment_plot) & !is.null(plot_index) & !segment_index %in% plot_index){
    plot_index <- unique(append(plot_index, segment_index))
  }
  if(is.null(attribute)){
    attribute <- paste(summarize_fun[[1]], segment_index[[1]], sep = ".")
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:4326")
  }
  nlyrs <- terra::nlyr(mosaic)
  if(is.null(basemap)){
    if(verbose){
      msg <- "Building the basemap..."
      cli::cli_progress_step(
        msg = msg,
        msg_done = sub("\\.\\.\\.$", "", msg),
        msg_failed = "Basemap creation failed"
      )
    }
    basemap <-
      suppressWarnings(
        mosaic_view(mosaic,
                    r = r,
                    g = g,
                    b = b,
                    re = re,
                    nir = nir,
                    swir = swir,
                    tir = tir,
                    max_pixels = max_pixels,
                    verbose = verbose,
                    downsample = downsample,
                    quantiles = quantiles,
                    edit = FALSE)
      )
  }
  if(is.null(shapefile)){
    if(verbose){
      msg <- "Building the shapefiles..."
      cli::cli_progress_step(
        msg = msg,
        msg_done = sub("\\.\\.\\.$", "", msg),
        msg_failed = "Shapefile creation failed"
      )
    }
    created_shapes <-
      suppressWarnings(
        shapefile_build(mosaic,
                        basemap = basemap,
                        grid = grid,
                        nrow = nrow,
                        ncol = ncol,
                        plot_width = plot_width,
                        plot_height = plot_height,
                        layout = layout,
                        build_shapefile = build_shapefile,
                        check_shapefile = check_shapefile,
                        sf_to_polygon =  TRUE,
                        buffer_edge = buffer_edge,
                        buffer_col = buffer_col,
                        buffer_row = buffer_row,
                        max_pixels = max_pixels,
                        verbose = FALSE,
                        downsample = downsample,
                        quantiles = quantiles)
      )
    # crop to the analyzed area
    if(crop_to_shape_ext){
      if(verbose){
        msg <- "Cropping the mosaic to the shapefile extent..."
        cli::cli_progress_step(
          msg = msg,
          msg_done = sub("\\.\\.\\.$", "", msg),
          msg_failed = "Mosaic cropping failed"
        )
      }
      ress <- terra::res(mosaic)
      if(sum(ress) != 2){
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      } else{
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      }
      mosaiccr <- terra::crop(mosaic, poly_ext)
    } else{
      mosaiccr <- mosaic
    }
  } else{
    if(inherits(shapefile, "list")){
      created_shapes <- lapply(shapefile, function(x){
        x
      })
    } else{
      if(inherits(shapefile, "SpatVector")){
        created_shapes <- sf::st_as_sf(shapefile) |> sf_to_polygon()
      }
      if(!"block" %in% colnames(shapefile)){
        cli::cli_abort("{.arg block} and {.arg plot_id} must be in the column names of shapefile")
      }
      created_shapes <- split(shapefile, shapefile$block)
    }
    if(crop_to_shape_ext){
      if(verbose){
        msg <- "Cropping the mosaic to the shapefile extent..."
        cli::cli_progress_step(
          msg = msg,
          msg_done = sub("\\.\\.\\.$", "", msg),
          msg_failed = "Mosaic cropping failed"
        )
      }
      ress <- terra::res(mosaic)
      if(sum(ress) != 2){
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      } else{
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      }

      mosaiccr <- terra::crop(mosaic, poly_ext)

    } else{
      mosaiccr <- mosaic
    }
  }
  segment_plot <- validate_and_replicate(segment_plot, created_shapes, verbose = verbose)
  segment_individuals <- validate_and_replicate(segment_individuals, created_shapes, verbose = verbose)
  threshold <- validate_and_replicate(threshold, created_shapes, verbose = verbose)
  watershed <- validate_and_replicate(watershed, created_shapes, verbose = verbose)
  segment_index <- validate_and_replicate(segment_index, created_shapes, verbose = verbose)
  invert <- validate_and_replicate(invert, created_shapes, verbose = verbose)
  includeopt <- validate_and_replicate(includeopt, created_shapes, verbose = verbose)
  opening <- validate_and_replicate(opening, created_shapes, verbose = verbose)
  closing <- validate_and_replicate(closing, created_shapes, verbose = verbose)
  filter <- validate_and_replicate(filter, created_shapes, verbose = verbose)
  erode <- validate_and_replicate(erode, created_shapes, verbose = verbose)
  dilate <- validate_and_replicate(dilate, created_shapes, verbose = verbose)
  grid <- validate_and_replicate(grid, created_shapes, verbose = verbose)
  lower_noise <- validate_and_replicate(lower_noise, created_shapes, verbose = verbose)

  if(!is.null(lower_size)){
    lower_size <- validate_and_replicate(lower_size, created_shapes, verbose = verbose)
  }
  if(!is.null(upper_size)){
    upper_size <- validate_and_replicate(upper_size, created_shapes, verbose = verbose)
  }
  if(!is.null(topn_lower)){
    topn_lower <- validate_and_replicate(topn_lower, created_shapes, verbose = verbose)
  }
  if(!is.null(topn_upper)){
    topn_upper <- validate_and_replicate(topn_upper, created_shapes, verbose = verbose)
  }
  #
  if(is.null(indexes)){
    if(verbose){
      msg <- "Computing vegetation indexes..."
      cli::cli_progress_step(
        msg = msg,
        msg_done = sub("\\.\\.\\.$", "", msg),
        msg_failed = "Index computation failed"
      )
    }
    if(nlyrs > 1 | !all(plot_index %in% names(mosaiccr))){
      mind <- terra::rast(
        Map(c,
            lapply(seq_along(plot_index), function(i){
              mosaic_index(mosaiccr,
                           index = plot_index[[i]],
                           r = r,
                           g = g,
                           b = b,
                           re = re,
                           nir = nir,
                           swir = swir,
                           tir = tir,
                           plot = FALSE)
            })
        )
      )
    } else{
      plot_index <- names(mosaiccr)
      mind <- mosaiccr
    }
  } else{
    mind <- indexes
    if(!all(segment_index %in% names(mind))){
      cli::cli_warn("{.val segment_index} must be present in `indexes`")
    }
  }

  results <- list()
  result_indiv <- list()
  extents <- terra::ext(mosaiccr)
  usepickmask <- segment_pick & (segment_individuals[[1]] | segment_plot[[1]])
  if(usepickmask){
    if(build_shapefile & is.null(shapefile)){
      mapview::mapview() |> mapedit::editMap()
    }
    mask <- suppressWarnings(
      mosaic_segment_pick(mosaic,
                          basemap = basemap,
                          r = r,
                          g = g,
                          b = b,
                          max_pixels = max_pixels,
                          return = "mask")
    )
  }
  ihaveamask <- !is.null(mask) & (segment_individuals[[1]] | segment_plot[[1]])
  if(ihaveamask){
    mask <- mask
  }
  if(verbose){
    cli::cli_progress_bar(
      format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
      total = length(created_shapes),
      clear = FALSE
    )
  }
  for(j in seq_along(created_shapes)){
    if(verbose){
      cli::cli_h2("Analyzing block {.val {j}}")
      cli::cli_progress_update()
      if(segment_plot[j] & segment_individuals[j]){
        cli::cli_abort("Only {.arg segment_plot} OR {.arg segment_individuals} can be used")
      }
    }
    if(inherits(created_shapes[[j]]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(created_shapes[[j]]$geometry[[1]])) == 5 & grid[[j]]){
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      if(crop_to_shape_ext){
        ress <- terra::res(mosaic)
        if(sum(ress) != 2){
          plot_grid <-
            plot_grid |>
            sf::st_transform(crs = sf::st_crs(terra::crs(mosaic)))
        }
        ext_anal <-
          plot_grid |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()

        mind_temp <- terra::crop(mind, terra::ext(ext_anal))
        if(!is.null(mask)){
          mask <- terra::crop(mask, terra::ext(ext_anal))
        }
      } else{
        mind_temp <- mind
      }
      extents <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(verbose){
          msg <- "Masking vegetation from ground..."
          cli::cli_progress_step(
            msg        = msg,
            msg_done   = "Vegetation masking completed",
            msg_failed = "Failed to mask vegetation from ground"
          )
        }
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            cli::cli_abort("{.arg segment_index} must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] > thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] < thresh
          }
        }
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
        # compute plot coverage

        tmp <- exactextractr::exact_extract(mind_temp,
                                            plot_grid,
                                            coverage_area = TRUE,
                                            force_df = TRUE,
                                            progress = FALSE)
        covered_area <-
          purrr::map_dfr(tmp, function(x){
            data.frame(covered_area = sum(na.omit(x)[, "coverage_area"]),
                       plot_area = sum(x[, "coverage_area"]))
          }) |>
          dplyr::mutate(coverage = covered_area / plot_area)


        plot_grid <- dplyr::bind_cols(plot_grid, covered_area)
        if(simplify){
          plot_grid <- plot_grid |> sf::st_simplify(preserveTopology = TRUE)
        }
        rm(tmp)

      }

      # check if segmentation is performed (analyze individuals)
      if(segment_individuals[j]){
        if(verbose){
          msg <- "Segmenting individuals within plots..."
          cli::cli_progress_step(
            msg        = msg,
            msg_done = sub("\\.\\.\\.$", "", msg),
            msg_failed = "Segmentation failed"
          )
        }
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }

        } else{
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        dmask <- EBImage::Image(matrix(mask, ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(is.numeric(erode[j]) & erode[j] > 0){
          dmask <- image_erode(dmask, size = erode[j])
        }
        if(is.numeric(dilate[j]) & dilate[j] > 0){
          dmask <- image_dilate(dmask, size = dilate[j])
        }
        if(is.numeric(opening[j]) & opening[j] > 0){
          dmask <- image_opening(dmask, size = opening[j])
        }
        if(is.numeric(closing[j]) & closing[j] > 0){
          dmask <- image_closing(dmask, size = closing[j])
        }
        if(watershed[j]){
          dmask <- EBImage::watershed(EBImage::distmap(dmask), tolerance = tolerance, ext = extension)
        } else{
          dmask <- EBImage::bwlabel(dmask)
        }
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        conts <- EBImage::ocontour(matrix(dmask, ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        conts <- conts[sapply(conts, nrow) > 2]
        sf_df <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extents[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extents[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        )
        if(simplify){
          sf_df <- sf_df |> sf::st_simplify(preserveTopology = TRUE)
        }
        centroids <- suppressWarnings(sf::st_centroid(sf_df))
        intersects <-
          switch (includeopt[j],
                  "intersect" = sf::st_intersects(sf_df, plot_grid),
                  "centroid" =  sf::st_within(centroids, plot_grid),
                  "covered" = sf::st_covered_by(sf_df, plot_grid),
                  "overlap" = sf::st_overlaps(sf_df, plot_grid),

          )
        plot_gridtmp <-
          plot_grid |>
          dplyr::mutate(plot_id_seq = paste0("P", leading_zeros(1:nrow(plot_grid), 4)))
        plot_id <- data.frame(plot_id_seq = paste0(intersects))
        valid_rows <- plot_id$plot_id_seq != "integer(0)"
        sf_df <- sf_df[valid_rows, ]
        plot_id <- paste0("P", leading_zeros(as.numeric(plot_id[valid_rows, ]), n = 4))

        gridindiv <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  })) |>
          dplyr::mutate(plot_id_seq = plot_id,
                        individual = paste0(1:nrow(sf_df)),
                        geometry = sf_df$geometry) |>
          dplyr::left_join(plot_gridtmp |> sf::st_drop_geometry(), by = dplyr::join_by(plot_id_seq)) |>
          dplyr::select(-plot_id_seq) |>
          dplyr::relocate(block, plot_id, individual, .before = 1) |>
          sf::st_sf()

        # control noise removing
        if(!is.null(lower_size[j]) & !is.null(topn_lower[j]) | !is.null(upper_size[j]) & !is.null(topn_upper[j])){
          cli::cli_abort("Only one of {.arg lower_* or} {.arg topn_*} can be used.")
        }
        ifelse(!is.null(lower_size[j]),
               gridindiv <- gridindiv[gridindiv$area > lower_size[j], ],
               gridindiv <- gridindiv[gridindiv$area > mean(gridindiv$area) * lower_noise[j], ])

        if(!is.null(upper_size[j])){
          gridindiv <- gridindiv[gridindiv$area < upper_size[j], ]
        }
        if(!is.null(topn_lower[j])){
          gridindiv <- gridindiv[order(gridindiv$area),][1:topn_lower[j],]
        }
        if(!is.null(topn_upper[j])){
          gridindiv <- gridindiv[order(gridindiv$area, decreasing = TRUE),][1:topn_upper[j],]
        }
        if(verbose){
          msg <- "Extracting features from segmented individuals..."
          cli::cli_progress_step(
            msg        = msg,
            msg_done = sub("\\.\\.\\.$", "", msg),
            msg_failed = "Feature extraction failed"
          )
        }

        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = sf::st_sf(gridindiv),
                                       fun = summarize_fun,
                                       quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          if(is.null(summarize_fun)){
            valindiv <- dplyr::bind_rows(valindiv, .id = "individual")
            if("coverage_fraction" %in% colnames(valindiv)){
              valindiv$coverage_fraction <- NULL
            }
            if("value" %in% colnames(valindiv)){
              colnames(valindiv)[2] <- plot_index
            }
            valindiv <- valindiv |> dplyr::nest_by(individual)
          } else{
            valindiv <-
              do.call(rbind, lapply(1:length(valindiv), function(i){
                tmp <- transform(valindiv[[i]],
                                 individual = paste0(i))
                tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

              }
              ))

            if(length(plot_index) == 1){
              colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
            } else{
              # colnames(valindiv) <- c("block", "plot_id", plot_index)
              colnames(vals) <- paste0(colnames(vals), ".", plot_index)
            }
          }
        } else{
          if(length(plot_index) == 1){
            colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
          }
        }
        if(!is.null(summarize_fun)){
          valindiv <-
            dplyr::bind_cols(gridindiv, valindiv) |>
            dplyr::mutate(individual = paste0(1:nrow(gridindiv)), .before = area) |>
            sf::st_sf()
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        } else{
          valindiv <- dplyr::bind_cols(dplyr::left_join(gridindiv, valindiv, by = dplyr::join_by(individual))) |> sf::st_sf()
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        }
      } else{
        dmask <- NULL
        result_indiv[[j]] <- NULL
      }

      # extract the values for the individual plots
      # check if a mask is used and no segmentation
      if(!is.null(mask) & (!segment_individuals[[1]] & !segment_plot[[1]])){
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
      }
      if(verbose){
        msg <- "Extracting plot-level features..."
        cli::cli_progress_step(
          msg        = msg,
          msg_done = sub("\\.\\.\\.$", "", msg),
          msg_failed = "Feature extraction of plots failed"
        )
      }
      vals <-
        exactextractr::exact_extract(x = mind_temp,
                                     y = plot_grid,
                                     fun = summarize_fun,
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

    } else{
      ####### ANY TYPE OF POLYGON ########
      # check if segmentation is performed
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      if(crop_to_shape_ext){
        ext_anal <-
          plot_grid |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
        mind_temp <- terra::crop(mind, terra::ext(ext_anal))
        if(!is.null(mask)){
          mask <- terra::crop(mask, terra::ext(ext_anal))
        }
      } else{
        mind_temp <- mind
      }
      extents <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            cli::cli_abort("{.arg segment_index} must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] > thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] < thresh
          }
        }
        # compute plot coverage
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
        tmp <- exactextractr::exact_extract(mind_temp,
                                            plot_grid,
                                            coverage_area = TRUE,
                                            force_df = TRUE,
                                            progress = FALSE)


        covered_area <-
          purrr::map_dfr(tmp, function(x){
            data.frame(covered_area = sum(na.omit(x)[, "coverage_area"]),
                       plot_area = sum(x[, "coverage_area"]))
          }) |>
          dplyr::mutate(coverage = covered_area / plot_area)
        plot_grid <- dplyr::bind_cols(plot_grid, covered_area)
        if(simplify){
          plot_grid <- plot_grid |> sf::st_simplify(preserveTopology = TRUE)
        }
        rm(tmp)

      }

      if(segment_individuals[j]){
        if(verbose){
          msg <- "Segmenting individuals..."
          cli::cli_progress_step(
            msg        = msg,
            msg_done = sub("\\.\\.\\.$", "", msg),
            msg_failed = sub("\\.\\.\\.$", "", msg)
          )
        }
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        extents <- terra::ext(mind_temp)
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(is.numeric(erode[j]) & erode[j] > 0){
          dmask <- image_erode(dmask, size = erode[j])
        }
        if(is.numeric(dilate[j]) & dilate[j] > 0){
          dmask <- image_dilate(dmask, size = dilate[j])
        }
        if(is.numeric(opening[j]) & opening[j] > 0){
          dmask <- image_opening(dmask, size = opening[j])
        }
        if(is.numeric(closing[j]) & closing[j] > 0){
          dmask <- image_closing(dmask, size = closing[j])
        }
        if(watershed[j]){
          dmask <- EBImage::watershed(EBImage::distmap(dmask), tolerance = tolerance, ext = extension)
        } else{
          dmask <- EBImage::bwlabel(dmask)
        }
        conts <- EBImage::ocontour(dmask)
        conts <- conts[sapply(conts, nrow) > 2]
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        sf_df <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extents[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extents[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        )
        if(simplify){
          sf_df <- sf_df |> sf::st_simplify(preserveTopology = TRUE)
        }
        centroids <- suppressWarnings(sf::st_centroid(sf_df))
        intersect_individ <-
          switch (includeopt[j],
                  "intersect" = sf::st_intersects(sf_df, plot_grid, sparse = FALSE)[,1],
                  "centroid" = sf::st_within(centroids, plot_grid, sparse = FALSE)[,1],
                  "covered" = sf::st_covered_by(sf_df, plot_grid, sparse = FALSE)[,1],
                  "overlap" = sf::st_overlaps(sf_df, plot_grid, sparse = FALSE)[,1],

          )
        sf_df <- sf_df[intersect_individ, ]
        addmeasures <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  }))
        gridindiv <- cbind(sf_df, addmeasures)

        # control noise removing
        if(!is.null(lower_size[j]) & !is.null(topn_lower[j]) | !is.null(upper_size[j]) & !is.null(topn_upper[j])){
          cli::cli_abort("Only one of {.arg lower_*} or {.arg topn_*} can be used.")
        }
        ifelse(!is.null(lower_size[j]),
               gridindiv <- gridindiv[gridindiv$area > lower_size[j], ],
               gridindiv <- gridindiv[gridindiv$area > mean(gridindiv$area) * lower_noise[j], ])
        if(!is.null(upper_size[j])){
          gridindiv <- gridindiv[gridindiv$area < upper_size[j], ]
        }
        if(!is.null(topn_lower[j])){
          gridindiv <- gridindiv[order(gridindiv$area),][1:topn_lower[j],]
        }
        if(!is.null(topn_upper[j])){
          gridindiv <- gridindiv[order(gridindiv$area, decreasing = TRUE),][1:topn_upper[j],]
        }
        if(verbose){
          msg <- "Extracting plant-level features..."
          cli::cli_progress_step(
            msg        = msg,
            msg_done = sub("\\.\\.\\.$", "", msg),
            msg_failed = sub("\\.\\.\\.$", "", msg)
          )
        }
        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = gridindiv,
                                       fun = summarize_fun,
                                       # quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          if(is.null(summarize_fun)){
            valindiv <- dplyr::bind_rows(valindiv, .id = "individual")
            if("coverage_fraction" %in% colnames(valindiv)){
              valindiv$coverage_fraction <- NULL
            }
            if("value" %in% colnames(valindiv)){
              colnames(valindiv)[2] <- plot_index
            }
            valindiv <- valindiv |> dplyr::nest_by(individual) |> dplyr::ungroup()
          } else{
            valindiv <-
              do.call(rbind, lapply(1:length(valindiv), function(i){
                tmp <- transform(valindiv[[i]],
                                 individual = paste0(i),
                                 block = paste0("B", leading_zeros(j, n = 2)))
                tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

              }
              ))

            if(length(plot_index) == 1){
              colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
            } else{
              # colnames(valindiv) <- c("block", "plot_id", plot_index)
              colnames(vals) <- paste0(colnames(vals), ".", plot_index)
            }
          }
        } else{
          if(length(plot_index) == 1){
            colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
          }
        }

        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), plot_id = "P0001", gridindiv, valindiv, check.names = FALSE)
          result_indiv[[j]] <- valindiv
        } else{
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), plot_id = "P0001", dplyr::left_join(gridindiv, valindiv, by = dplyr::join_by(individual)), check.names = FALSE)
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        }
      } else{
        result_indiv[[j]] <- NULL
      }

      # extract the values for the individual plots
      # check if a mask is used and no segmentation
      if(!is.null(mask) & (!segment_individuals[[1]] & !segment_plot[[1]])){
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
      }
      if(verbose){
        msg <- "Extracting plot-level features..."
        cli::cli_progress_step(
          msg        = msg,
          msg_done = sub("\\.\\.\\.$", "", msg),
          msg_failed = sub("\\.\\.\\.$", "", msg)
        )
      }
      vals <-
        exactextractr::exact_extract(x = mind_temp,
                                     y = plot_grid,
                                     fun = summarize_fun,
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))
    }
    # bind the results
    if(verbose){
      msg <- "Binding the extracted features..."
      cli::cli_progress_step(
        msg        = msg,
        msg_done = sub("\\.\\.\\.$", "", msg),
        msg_failed = sub("\\.\\.\\.$", "", msg)
      )
    }
    if(inherits(vals, "list")){
      names(vals) <- paste0(plot_grid$block, "_", plot_grid$plot_id)
      vals <- dplyr::bind_rows(vals, .id = "plot") |> pliman::separate_col(plot, c("block", "plot_id"))
      if("coverage_fraction" %in% colnames(vals)){
        vals$coverage_fraction <- NULL
      }
      if(length(plot_index) == 1){
        if(ncol(vals) == 3){
          colnames(vals)[3] <- plot_index
        }
      }
      vals <-
        vals |>
        dplyr::nest_by(block, plot_id, row, column) |>
        dplyr::ungroup() |>
        dplyr::left_join(plot_grid, by = dplyr::join_by(block, plot_id, row, column))
    } else{
      if(length(plot_index) == 1){
        if(ncol(vals) == 1){
          colnames(vals) <- paste0(colnames(vals), ".", plot_index)
        } else{
          colnames(vals) <- paste0(colnames(vals), ".", plot_index)
        }
      }
      vals <- dplyr::bind_cols(plot_grid, vals)
    }
    results[[j]] <- vals
  }


  # bind the results  ## at a level plot
  if(verbose){
    msg <- "Summarizing the results..."
    cli::cli_progress_step(
      msg        = msg,
      msg_done = sub("\\.\\.\\.$", "", msg),
      msg_failed = sub("\\.\\.\\.$", "", msg)
    )
  }

  results <- dplyr::bind_rows(results) |> sf::st_sf()
  if(any(segment_individuals)){
    result_indiv <- do.call(rbind, result_indiv)
    blockid <- unique(result_indiv$block)

    summres <-
      lapply(1:length(blockid), function(i){
        # Retain geometry column separately
        nonndata <-
          result_indiv |>
          dplyr::filter(block == blockid[i]) |>
          sf::st_drop_geometry() |>
          dplyr::select(plot_id, !where(is.numeric), - block) |>
          dplyr::group_by(plot_id) |>
          dplyr::slice(1)

        result_indiv |>
          dplyr::filter(block == blockid[i]) |>
          as.data.frame() |>
          dplyr::group_by(plot_id) |>
          dplyr::summarise(
            area_sum = sum(area, na.rm = TRUE),
            n = length(area),
            dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE))
          ) |>
          dplyr::mutate(block = blockid[i], .before = 1) |>
          dplyr::ungroup() |>
          dplyr::relocate(n, .after = plot_id) |>   # Join the geometry column back
          dplyr::left_join(nonndata, by = "plot_id")
      })
    names(summres) <- blockid
    # compute plot area
    plot_area <-
      results |>
      dplyr::mutate(plot_area = sf::st_area(geometry)) |>
      as.data.frame(check.names = FALSE) |>
      dplyr::select(block, plot_id, row, column, plot_area)

    # compute coverage area
    result_plot_summ <- do.call(rbind, lapply(summres, function(x){x}))
    if(!"row" %in% colnames(result_plot_summ)){
      result_plot_summ <-
        result_plot_summ |>
        dplyr::mutate(row = 1,
                      column = 1,
                      .after = plot_id)
    }

    result_plot_summ <-
      result_plot_summ |>
      dplyr::left_join(plot_area, by = dplyr::join_by(block, plot_id, row, column)) |>
      dplyr::mutate(coverage = as.numeric(area_sum / plot_area), .after = area) |>
      dplyr::left_join(results |>   dplyr::select(block, plot_id, row, column, geometry),
                       by = dplyr::join_by(block, plot_id, row, column)) |>
      sf::st_as_sf()

    centroid <-
      suppressWarnings(sf::st_centroid(sf::st_sf(result_indiv)) |>
                         sf::st_coordinates() |>
                         as.data.frame() |>
                         setNames(c("x", "y")))
    result_indiv <-
      dplyr::bind_cols(result_indiv, centroid) |>
      dplyr::relocate(x, y, .after = individual)

    if(map_individuals){
      if(verbose){
        msg <- "Mapping individuals within plots..."
        cli::cli_progress_step(
          msg        = msg,
          msg_done = sub("\\.\\.\\.$", "", msg),
          msg_failed = sub("\\.\\.\\.$", "", msg)
        )
      }
      dists <-
        result_indiv |>
        sf::st_drop_geometry() |>
        dplyr::select(block, plot_id, row, column, x, y) |>
        dplyr::group_by(block, plot_id, row, column)

      splits <- dplyr::group_split(dists)
      names(splits) <- dplyr::group_keys(dists) |> dplyr::mutate(key = paste0(block, "_", plot_id)) |> dplyr::pull()
      dists <- lapply(splits, compute_dists)

      cvs <- sapply(dists, function(x){
        (sd(x) / mean(x)) * 100
      })
      means <- sapply(dists, mean)

      result_plot_summ <-
        result_plot_summ |>
        dplyr::mutate(mean_distance = means,
                      cv = cvs,
                      .before = n)
      result_individ_map <- list(distances = dists,
                                 means = means,
                                 cvs = cvs)
    } else{
      result_individ_map <- NULL
    }

  } else{
    result_plot_summ <- NULL
    result_indiv <- NULL
    result_individ_map <- NULL

  }

  if(isTRUE(plot)){
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- mosaic_aggregate(mosaiccr, pct = round(100 / downsample))
    }
    if(any(segment_individuals)){
      dfplot <- result_plot_summ
      if(!attribute %in% colnames(dfplot)){
        attribute <- "area"
      }
    } else{
      dfplot <- results
      if(!attribute %in% colnames(dfplot)){
        attribute <- NULL
      }
    }
    if(is.null(summarize_fun)){
      check_and_install_package("tidyr")
      dfplot <-
        dfplot |>
        sf::st_drop_geometry() |>
        tidyr::unnest(cols = data) |>
        dplyr::group_by(block, plot_id, row, column) |>
        dplyr::summarise(dplyr::across(where(is.numeric), \(x){mean(x, na.rm = TRUE)}), .groups = "drop") |>
        dplyr::left_join(dfplot |> dplyr::select(block, plot_id, row, column, geometry),
                         by = dplyr::join_by(block, plot_id, row, column)) |>
        sf::st_sf()

    }
    map <-
      basemap +
      suppressWarnings(
        mapview::mapview(dfplot,
                         zcol = attribute,
                         layer.name = attribute,
                         col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                         alpha.regions = 0.75,
                         na.color = "#00000000",
                         maxBytes = 64 * 1024 * 1024,
                         verbose = FALSE)
      )

    if(any(segment_individuals)){
      attribute <- ifelse(!attribute %in% colnames(result_indiv), "area", attribute)
      mapindivid <-
        basemap +
        suppressWarnings(
          mapview::mapview(result_plot_summ,
                           legend = FALSE,
                           alpha.regions = 0.4,
                           zcol = "block",
                           map.types = "OpenStreetMap") +
            mapview::mapview(result_indiv,
                             zcol = attribute,
                             layer.name = attribute,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
    } else{
      mapindivid <- NULL
    }
  } else{
    map <- NULL
    mapindivid <- NULL
  }
  if(verbose){
    cli::cli_rule(
      left = cli::col_blue("Mosaic successfully analyzed"),
      right = cli::col_blue("Finished on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
    )
  }
  return(list(result_plot = results,
              result_plot_summ = result_plot_summ,
              result_indiv = result_indiv,
              result_individ_map = result_individ_map,
              map_plot = map,
              map_indiv = mapindivid,
              shapefile = created_shapes))
}

#' Analyze mosaics iteratively
#'
#' High-resolution mosaics can take a significant amount of time to analyze,
#' especially when `segment_individuals = TRUE` is used in mosaic_analyze().
#' This is because the function needs to create in-memory arrays to segment
#' individual using the watershed algorithm. This process utilizes a for-loop
#' approach, iteratively analyzing each shape within the mosaic one at a time.
#' To speed up processing, the function crops the original mosaic to the extent
#' of the current shape before analyzing it. This reduces the resolution for
#' that specific analysis, sacrificing some detail for faster processing.
#'
#' @inheritParams mosaic_analyze
#' @inheritParams analyze_objects
#' @param ... Further arguments passed on to [mosaic_analyze()]
#'
#' @return A list containing the following objects:
#' * `result_plot`: The results at a plot level.
#' *  `result_plot_summ`: The summary of results at a plot level. When
#'  `segment_individuals = TRUE`, the number of individuals, canopy coverage,
#'  and mean values of some shape statistics such as perimeter, length, width,
#'  and diameter are computed.
#' * `result_individ`: The results at an individual level.
#' * `map_plot`: An object of class `mapview` showing the plot-level results.
#' * `map_individual`: An object of class `mapview` showing the individual-level
#'   results.
#' @export
#'


mosaic_analyze_iter <- function(mosaic,
                                shapefile,
                                basemap = NULL,
                                r = 3,
                                g = 2,
                                b = 1,
                                re = NA,
                                nir = NA,
                                swir = NA,
                                tir = NA,
                                plot = TRUE,
                                verbose = TRUE,
                                max_pixels = 3e6,
                                attribute = NULL,
                                summarize_fun = "mean",
                                segment_plot = FALSE,
                                segment_individuals = FALSE,
                                segment_index = "VARI",
                                plot_index =  "VARI",
                                color_regions = rev(grDevices::terrain.colors(50)),
                                alpha = 0.75,
                                quantiles = c(0, 1),
                                parallel = FALSE,
                                workers = NULL,
                                ...){
  pind <- unique(c(plot_index, segment_index))
  if(is.null(attribute)){
    attribute <- paste(summarize_fun, pind[[1]], sep = ".")
  }
  shapefile <- shapefile_input(shapefile, info = FALSE)

  if(terra::inMemory(mosaic)){
    tf <- paste0(tempfile(), ".tif")
    on.exit(file.remove(tf))
    terra::writeRaster(mosaic, filename = tf)
    tempf <- tf
  } else{
    tempf <- terra::sources(mosaic)
  }
  rast_file <- terra::rast(tempf)
  if (verbose) {
    cli::cli_progress_step(
      msg = "Clipping mosaic plots to temporary files...",
      msg_done = "Clipping done.",
      msg_failed = "Oops, something went wrong."
    )
  }

  tempfiles <- mosaic_clip(mosaic, shapefile, out_dir = tempdir(), verbose = FALSE, overwrite = TRUE)
  on.exit(unlink(tempfiles), add = TRUE)

  if(parallel){
    nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.4), workers)
    mirai::daemons(nworkers)
    on.exit(mirai::daemons(0))

    if (verbose) {
      cli::cli_rule(
        left = cli::col_blue("Parallel processing using {nworkers} cores"),
        right = cli::col_blue("Started on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
      )
      cli::cli_progress_step(
        msg        = "Initializing {.strong {nworkers}} Mirai workers for parallel processing...",
        msg_done   = "Parallel environment configured",
        msg_failed = "{.cross} Failed to initialize workers: {.emph {err$message}}"
      )
    }
    shapes <- split(shapefile, seq_len(nrow(shapefile)))
    worker_fun <- function(path, shp) {
      terra::rast(path) |>
        pliman::mosaic_analyze(shapefile = shp,
                               basemap               = NULL,
                               r                     = r, g = g, b = b,
                               re = re, nir = nir, swir = swir, tir = tir,
                               segment_individuals   = segment_individuals,
                               segment_index         = segment_index,
                               segment_plot          = segment_plot,
                               plot_index            = pind,
                               build_shapefile       = FALSE,
                               plot                  = FALSE,
                               grid                  = FALSE,
                               verbose               = FALSE,
                               crop_to_shape_ext     = FALSE,
                               ...)
    }

    bind <-
      mirai::mirai_map(
        seq_along(tempfiles),
        function(i){worker_fun(tempfiles[i], shapes[[i]])}
      )[.progress]
    # return(bind)

  } else{
    bind <- list()
    if(verbose){
      cli::cli_rule(
        left = cli::col_blue("Analyzing the mosaic by plot using a sequential approach"),
        right = cli::col_blue("Started on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        total  = nrow(shapefile),
        clear  = TRUE
      )
    }
    for (i in seq_along(tempfiles)) {
      if(verbose){
        cli::cli_progress_update()
      }
      bind[[paste0("P", leading_zeros(i, 4))]] <-
        mosaic_analyze(terra::rast(tempfiles[[i]]),
                       basemap = basemap,
                       r = r, g = g, b = b, re = re, nir = nir, swir = swir, tir = tir,
                       shapefile = shapefile[i, ],
                       segment_individuals = segment_individuals,
                       segment_index = segment_index,
                       segment_plot = segment_plot,
                       plot_index = pind,
                       build_shapefile = FALSE,
                       plot = FALSE,
                       grid = FALSE,
                       verbose = FALSE,
                       crop_to_shape_ext = FALSE,
                       ...)
    }
  }
  if(is.null(bind[[1]]$result_individ_map)){
    result_individ_map <- NULL
  }
  if(is.null(bind[[1]]$result_indiv)){
    result_indiv <- result_plot_summ <- NULL
  } else{
    result_indiv <- dplyr::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_indiv
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      dplyr::relocate(plot_id, .after = block) |>
      sf::st_as_sf()

    result_plot_summ <- dplyr::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_plot_summ
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      dplyr::relocate(plot_id, .after = block) |>
      sf::st_as_sf()
  }

  result_plot <- dplyr::bind_rows(
    lapply(bind, function(x){
      tmp <- x$result_plot
      tmp$plot_id <- NULL
      tmp
    }),
    .id = "plot_id"
  ) |>
    dplyr::relocate(plot_id, .after = block) |>
    sf::st_as_sf()



  if(isTRUE(plot)){
    if(is.null(basemap)){
      if(terra::nlyr(mosaic) < 3){
        basemap <-
          suppressWarnings(
            mapview::mapview(mosaic,
                             maxpixels = 5e6,
                             legend = FALSE,
                             map.types = "CartoDB.Positron",
                             alpha.regions = 1,
                             na.color = "transparent",
                             verbose = FALSE)
          )
      } else{
        basemap <-
          suppressWarnings(
            mapview::viewRGB(
              as(mosaic, "Raster"),
              layer.name = "base",
              r = r,
              g = g,
              b = b,
              na.color = "#00000000",
              maxpixels = 5e6,
              quantiles = quantiles
            )
          )
      }
    }
    if(!is.null(result_indiv)){
      dfplot <- result_plot_summ
    } else{
      dfplot <- result_plot
    }
    # plot level
    map <-
      basemap +
      suppressWarnings(
        mapview::mapview(dfplot,
                         zcol = ifelse(!attribute %in% colnames(dfplot), "area", attribute),
                         layer.name = ifelse(!attribute %in% colnames(dfplot), "area", attribute),
                         col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                         alpha.regions = 0.75,
                         na.color = "#00000000",
                         maxBytes = 64 * 1024 * 1024,
                         verbose = FALSE)
      )
    # individual plot
    if(!is.null(result_indiv)){
      mapindivid <-
        basemap +
        suppressWarnings(
          mapview::mapview(result_plot_summ,
                           alpha.regions = 0.4,
                           zcol = attribute,
                           col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                           map.types = "OpenStreetMap") +
            mapview::mapview(result_indiv,
                             zcol = ifelse(!attribute %in% colnames(result_indiv), "area", attribute),
                             layer.name = ifelse(!attribute %in% colnames(result_indiv), "area", attribute),
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
    } else{
      mapindivid <- NULL
    }
  } else{
    map <- NULL
    mapindivid <- NULL
  }
  if (verbose) {
    cli::cli_rule(
      left = cli::col_blue("Features successfully extracted"),
      right = cli::col_blue("Finished on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
    )
  }

  return(list(result_plot = result_plot,
              result_plot_summ = result_plot_summ,
              result_indiv = result_indiv,
              result_individ_map = result_individ_map,
              map_plot = map,
              map_indiv = mapindivid))
}



#' Mosaic View
#'
#' @details
#' The function can generate either an interactive map using the 'mapview'
#' package or a static plot using the 'base' package, depending on the `viewer`
#' and `show` parameters. If show = "index" is used, the function first computes
#' an image index that can be either an RGB-based index or a multispectral
#' index, if a multispectral mosaic is provided.
#'
#' @param mosaic A mosaic of class `SpatRaster`, generally imported with
#'   [mosaic_input()].
#' @inheritParams image_view
#' @inheritParams image_align
#' @inheritParams mosaic_index
#' @param edit If `TRUE` enable editing options using [mapedit::editMap()].
#' @param title A title for the generated map or plot (default: "").
#' @param shapefile An optional shapefile of class `sf` to be plotted over the
#'   mosaic. It can be, for example, a plot-level result returned by
#'   [mosaic_analyze()].
#' @param attribute The attribute name(s) or column number(s) in shapefile table
#'   of the column(s) to be rendered.
#' @param max_pixels Maximum number of pixels to render in the map or plot
#'   (default: 500000).
#' @param downsample Downsampling factor to reduce the number of pixels
#'   (default: NULL). In this case, if the number of pixels in the image (width
#'   x height) is greater than `max_pixels` a downsampling factor will be
#'   automatically chosen so that the number of plotted pixels approximates the
#'   `max_pixels`.
#' @param downsample_fun The resampling function. Defaults to nearest. See further details in [mosaic_aggregate()].
#' @param alpha opacity of the fill color of the raster layer(s).
#' @param quantiles the upper and lower quantiles used for color stretching.
#' @param axes logical. Draw axes? Defaults to `FALSE`.
#' @param ... Additional arguments passed on to [terra::plot()] when `viewer =
#'   "base"`.
#' @return An sf object, the same object returned by [mapedit::editMap()].
#'
#' @importFrom terra rast crs nlyr terraOptions
#' @importFrom methods as
#' @importFrom sf st_crs st_transform st_make_grid st_intersection st_make_valid
#' @importFrom dplyr summarise across mutate arrange left_join bind_cols
#'   bind_rows contains ends_with everything between where select filter
#'   relocate rename
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' # Load a raster showing the elevation of Luxembourg
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # Generate an interactive map using 'mapview'
#' mosaic_view(mosaic)
#'
#' # Generate a static plot using 'base'
#' mosaic_view(mosaic, viewer = "base")
#' }
#'
#'
#' @export
mosaic_view <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        edit = FALSE,
                        title = "",
                        shapefile = NULL,
                        attribute = NULL,
                        viewer = c("mapview", "base"),
                        show = c("rgb", "index"),
                        index = "B",
                        max_pixels = 1000000,
                        downsample = NULL,
                        downsample_fun = "nearest",
                        alpha = 1,
                        quantiles = c(0, 1),
                        color_regions = custom_palette(c("red", "yellow", "forestgreen")),
                        axes = FALSE,
                        ...){
  terra::terraOptions(progress = 0)
  on.exit(terra::terraOptions(progress = 1))
  # check_mapview()
  mapview::mapviewOptions(layers.control.pos = "topright", raster.size = 64 * 1024 * 1024)
  on.exit(mapview::mapviewOptions(default = TRUE))
  viewopt <- c("rgb", "index")
  viewopt <- viewopt[pmatch(show[[1]], viewopt)]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]

  if(inherits(mosaic, "Image")){
    mosaic <- terra::rast(EBImage::transpose(mosaic)@.Data)
  }
  if(viewopt == "rgb" & vieweropt == "base" & terra::nlyr(mosaic) > 1){
    cli::cli_warn("{.arg viewer = 'base'} can only be used with {.arg show = 'index'}. Defaulting to {.arg viewer = 'mapview'}")
    vieweropt <- "mapview"
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("+proj=utm +zone=05 +datum=WGS84 +units=m")
  }
  dimsto <- dim(mosaic)
  nr <- dimsto[1]
  nc <- dimsto[2]

  if (max_pixels > 2000000) {
    cli::cli_inform(c("i" = "The number of pixels is {.strong very high}, which might slow the rendering process."))
  }

  dwspf <- find_aggrfact(mosaic, max_pixels = max_pixels)

  if (dwspf > 0 && is.null(downsample)) {
    cli::cli_inform(c("i" = "Using {.code downsample = {dwspf}} to match the {.field max_pixels} constraint."))
    mosaic <- mosaic_aggregate(mosaic, pct = round(100 / dwspf), fun = downsample_fun)
  }

  if(viewopt == "index" & terra::nlyr(mosaic) > 2){
    mosaic <- mosaic_index(mosaic, index = index, plot = FALSE)
  }
  if(viewopt == "rgb"){
    if(terra::nlyr(mosaic) > 2){
      if(is.null(shapefile)){
        map <-
          mapview::viewRGB(as(mosaic[[c(r, g, b)]], "Raster"),
                           na.color = "#00000000",
                           layer.name = "base",
                           r = 1,
                           g = 2,
                           b = 3,
                           maxpixels = 60000000,
                           quantiles = quantiles)
        if(edit){
          map <-
            mapedit::editMap(map,
                             editor = "leafpm",
                             title = title)
        }
        map
      } else{
        mapview::viewRGB(as(mosaic[[c(r, g, b)]], "Raster"),
                         na.color = "#00000000",
                         layer.name = "base",
                         r = 1,
                         g = 2,
                         b =3,
                         maxpixels = 60000000,
                         quantiles = quantiles) +
          suppressWarnings(
            mapview::mapview(shapefile,
                             zcol = attribute,
                             layer.name = attribute,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
      }

    } else{
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        if(is.null(shapefile)){
          index <- gsub("[/\\\\]", "_", index, perl = TRUE)
          map <-
            mapview::mapview(mosaic,
                             map.types = mapview::mapviewGetOption("basemaps"),
                             layer.name = index,
                             maxpixels =  max_pixels,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE)

          if(edit){
            map <-
              map |>
              mapedit::editMap(editor = "leafpm",
                               title = title)
          }
          map
        } else{
          mapview::mapview(as(mosaic, "Raster"),
                           na.color = "#00000000",
                           layer.name = "",
                           maxpixels = 60000000) +
            suppressWarnings(
              mapview::mapview(shapefile,
                               zcol = attribute,
                               layer.name = attribute,
                               col.regions = color_regions,
                               alpha.regions = alpha,
                               na.color = "#00000000",
                               maxBytes = 64 * 1024 * 1024,
                               verbose = FALSE))

        }
      }
    }
  } else{
    if(terra::nlyr(mosaic) > 2){
      index <- gsub("[/\\\\]", "_", index, perl = TRUE)
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(mosaic,
                           layer.name = index,
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
        if(edit){
          map <-
            map |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        }
        map
      }
    } else{
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(mosaic,
                           layer.name = names(mosaic),
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
        if(edit){
          map <-
            map |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        }
        map
      }
    }
  }
}

#' Create and Export mosaics
#' @details
#' * `mosaic_input()` is a simply wrapper around [terra::rast()]. It creates a
#' `SpatRaster` object from scratch, from a filename, or from another object.
#' * `mosaic_export()` is a simply wrapper around [terra::writeRaster()]. It write
#' a `SpatRaster` object to a file.
#'
#' @name mosaic_input
#' @param mosaic
#'  * For `mosaic_input()`, a file path to the raster to imported, a matrix,
#'    array or a list of `SpatRaster` objects.
#'  * For `mosaic_export()`, an `SpatRaster` object.
#' @param mosaic_pattern A pattern name to import multiple mosaics into a list.
#' @param info Print the mosaic informations (eg., CRS, extent). Defaults to `TRUE`
#' @param check_16bits Checks if mosaic has maximum value in the 16-bits format
#'   (65535), and replaces it by NA. Defaults to `FALSE`.
#' @param check_datatype Logical. If \code{TRUE}, checks and suggests the
#'   appropriate data type based on the raster values.
#' @param filename character. The Output filename.
#' @param datatype The datatype. By default, the function will try to guess the
#'   data type that saves more memory usage and file size. See
#'   [terra::writeRaster()] and [terra::datatype()] for more details.
#' @param overwrite logical. If `TRUE`, filename is overwritten.
#' @param ... Additional arguments passed to [terra::rast()] (`mosaic_input()`)
#'   or  [terra::writeRaster()] (`mosaic_output()`)
#'
#' @return
#' * `mosaic_input()` returns an `SpatRaster` object.
#' * `mosaic_export()` do not return an object.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'
#' # create an SpatRaster object based on a matrix
#' x <- system.file("ex/logo.tif", package="terra")
#' rast <- mosaic_input(x)
#' mosaic_plot(rast)
#'
#' # create a temporary filename for the example
#' f <- file.path(tempdir(), "test.tif")
#' mosaic_export(rast, f, overwrite=TRUE)
#' list.files(tempdir())
#' }
#'
mosaic_input <- function(mosaic,
                         mosaic_pattern = NULL,
                         info = TRUE,
                         check_16bits = FALSE,
                         check_datatype = FALSE,
                         ...){
  if(file_extension(mosaic) %in% c("jpg", "jpeg", "png")){
    flip <- TRUE
  } else{
    flip <- FALSE
  }
  if(!is.null(mosaic_pattern)){
    if(mosaic_pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      mosaic_pattern <- "^[0-9].*$"
    }
    path <- getwd()
    imgs <- list.files(pattern = mosaic_pattern, path)
    if(length(grep(mosaic_pattern, imgs)) == 0){
      cli::cli_abort(c(
        "!" = "The specified {.arg mosaic_pattern} was not found.",
        "x" = "Pattern {.val {mosaic_pattern}} not found in directory {.path {dir}}."
      ))

    }
    list_img <-
      lapply(imgs, function(x){
        mosaic_input(x, info = FALSE)
      })
    names(list_img) <- imgs
    invisible(list_img)
  } else{
    mosaic <- suppressWarnings(terra::rast(mosaic, ...))
    if(terra::crs(mosaic) == ""){
      cli::cli_alert_info("Missing Coordinate Reference System. Setting to EPSG:3857")
      terra::crs(mosaic) <- terra::crs("EPSG:3857")
    }
    if(flip){
      mosaic <- mosaic_rotate(mosaic, 180, "anticlockwise")
    }
    if (terra::is.lonlat(mosaic)) {
      eps <- mosaic_epsg(mosaic)
      cli::cli_warn(c(
        "!" = "The current raster is in a {.emph lat/lon} coordinate system, which may lead to processing errors in {.fn mosaic_analyze()}.",
        "i" = "It is highly recommended to reproject the raster using {.fn mosaic_project()} with {.val {eps}}."
      ))
    }

    if(check_16bits | check_datatype){
      cels <- sample(1:terra::ncell(mosaic), 2000, replace = TRUE)
      a <- na.omit(unlist(terra::extract(mosaic, cels)))
      a <- a[!is.infinite(a)]
      if(inherits(a, "numeric")){
        if(check_datatype){
          if(length(a[a - floor(a) != 0]) == 0){
            minv <- min(a)
            maxv <- max(a)
            if(all(minv >= 0) & all(maxv <= 255)){
              datatype <- "INT1U"
            } else{
              datatype <- "INT2U"
            }
          } else{
            datatype <- "FLT4S"
          }
          dtterra <- terra::datatype(mosaic)[[1]]
          if (datatype != dtterra) {
            cli::cli_warn(c(
              "!" = "The detected datatype is {.val {dtterra}}, but the suggested datatype is {.val {datatype}}.",
              "i" = "Consider using {.fn mosaic_export()} to convert the datatype.",
              "i" = "Using the suggested datatype may reduce file size and improve memory efficiency during index computation."
            ))
          }

        }
        if(check_16bits){
          if(max(suppressWarnings(terra::minmax(mosaic)), na.rm = TRUE) == 65535){
            mosaic[mosaic == 65535] <- NA
          }
        }
      }
    }
    if(info){
      print(mosaic)
    }
    return(mosaic)
  }
}

#' @export
#' @name mosaic_input
mosaic_export <- function(mosaic,
                          filename,
                          datatype = NULL,
                          overwrite = FALSE,
                          ...){
  cels <- sample(1:terra::ncell(mosaic), 2000)
  a <- na.omit(unlist(terra::extract(mosaic, cels)))
  a <- a[!is.infinite(a)]

  if(is.null(datatype)){
    if(length(a[a - floor(a) != 0]) == 0){
      minv <- min(a)
      maxv <- max(a)
      if(all(minv >= 0) & all(maxv <= 255)){
        datatype <- "INT1U"
      } else{
        datatype <- "INT2U"
      }
    } else{
      datatype <- "FLT4S"
    }
  }
  cli::cli_inform("Exporting the mosaic using {.arg datatype} = {.val {datatype}}.")

  terra::writeRaster(mosaic,
                     filename = filename,
                     overwrite = overwrite,
                     datatype = datatype,
                     gdal=c("COMPRESS=DEFLATE", "BIGTIFF=IF_NEEDED"),
                     ...)
}


#' A wrapper around terra::resample()
#'
#' Transfers values between SpatRaster objects that do not align (have a
#' different origin and/or resolution). See [terra::resample()] for more details
#'
#' @param mosaic SpatRaster to be resampled
#' @param y SpatRaster with the geometry that x should be resampled to
#' @param ... Further arguments passed on to [terra::resample()].
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' library(terra)
#' r <- rast(nrows=3, ncols=3, xmin=0, xmax=10, ymin=0, ymax=10)
#' values(r) <- 1:ncell(r)
#' s <- rast(nrows=25, ncols=30, xmin=1, xmax=11, ymin=-1, ymax=11)
#' x <- mosaic_resample(r, s, method="bilinear")
#' opar <- par(no.readonly =TRUE)
#' par(mfrow=c(1,2))
#' plot(r)
#' plot(x)
#' par(opar)
#' }
mosaic_resample <- function(mosaic, y, ...){
  terra::resample(mosaic, y, ...)
}


#' SpatRaster aggregation
#'
#' Aggregate a SpatRaster to create a new SpatRaster with a lower resolution
#' (larger cells), using the GDAL's gdal_translate utility
#' https://gdal.org/programs/gdal_translate.html
#'
#' @param mosaic SpatRaster
#' @param pct The size as a fraction (percentage) of the input image size.
#'   Either a scalar (eg., 50), or a length-two numeric vector. In the last,
#'   different percentage reduction/expansion can be used for columns, and rows,
#'   respectively.
#' @param fun The resampling function. Defaults to `nearest`, which applies the
#'   nearest neighbor (simple sampling) resampler. Other accepted values are:
#'   'average', 'rms', 'bilinear', 'cubic', 'cubicspline', 'lanczos', and
#'   'mode'. See Details for a detailed explanation.
#' @param in_memory Wheter to return an 'in-memory' `SpatRaster`. If `FALSE`,
#'   the aggregated raster will be returned as an 'in-disk' object.
#' @return SpatRaster
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' library(terra)
#' r <- rast()
#' values(r) <- 1:ncell(r)
#' r2 <- mosaic_aggregate(r, pct = 10)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2))
#' mosaic_plot(r)
#' mosaic_plot(r2)
#' par(opar)
#' }
mosaic_aggregate <- function(mosaic,
                             pct = 50,
                             fun = "nearest",
                             in_memory = TRUE){
  outsize <- compute_outsize(pct)
  td <- tempdir()
  if(terra::inMemory(mosaic)[[1]]){
    in_raster <- file.path(td, "tmp_aggregate.tif")
    terra::writeRaster(mosaic, in_raster, overwrite = TRUE)
    on.exit({
      file.remove(in_raster)
      if(in_memory){
        file.remove(out_raster)
      }
    })
  } else{
    in_raster <- terra::sources(mosaic)[[1]]
    on.exit(
      if(in_memory){
        file.remove(out_raster)
      }
    )
  }
  out_raster <- file.path(td, "tmp_aggregate_small.tif")
  sf::gdal_utils(
    util = "translate",
    source = in_raster,
    destination = out_raster,
    options = strsplit(paste("-r", fun, "-outsize", outsize[1], outsize[2]), split = "\\s")[[1]]
  )
  if(in_memory){
    terra::rast(out_raster) |> terra::wrap() |> terra::unwrap()
  } else{
    terra::rast(out_raster)
  }
}


#' A wrapper around terra::plot()
#'
#' Plot the values of a SpatRaster
#'
#' @param mosaic SpatRaster
#' @param col character vector to specify the colors to use. Defaults to
#'   `custom_palette(c("red", "yellow", "forestgreen"))`.
#' @param ... Further arguments passed on to [terra::plot()].
#' @param smooth logical. If TRUE (default) the cell values are smoothed (only
#'   if a continuous legend is used).
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' r <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_plot(r)
#' }
mosaic_plot <- function(mosaic,
                        col = custom_palette(c("red", "yellow", "forestgreen"), n = 200),
                        smooth = TRUE,
                        ...){
  if(!inherits(mosaic, "SpatRaster")){
    cli::cli_abort("{.arg mosaic} must be an object of class {.cls SpatRaster}.")

  }
  terra::plot(mosaic,
              col = col,
              smooth = smooth,
              ...)
}

#' A wrapper around terra::hist()
#'
#' Create a histogram of the values of a `SpatRaster`.
#'
#' @param mosaic SpatRaster
#' @param layer positive integer or character to indicate layer numbers (or
#'   names). If missing, all layers are used
#' @param ... Further arguments passed on to [terra::hist()].
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' r <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_hist(r)
#' }
mosaic_hist <- function(mosaic, layer, ...){
  if(!inherits(mosaic, "SpatRaster")){
    cli::cli_abort("{.arg mosaic} must be an object of class {.cls SpatRaster}.")
  }
  terra::hist(mosaic, layer, ...)
}

#' A wrapper around terra::plotRGB()
#'
#' Plot the RGB of a SpatRaster
#'
#' @param mosaic SpatRaster
#' @param ... Further arguments passed on to [terra::plotRGB()].
#'
#' @return A `NULL` object
#' @export
#'
mosaic_plot_rgb <- function(mosaic, ...){
  if(!inherits(mosaic, "SpatRaster")){
    cli::cli_abort("{.arg mosaic} must be an object of class {.cls SpatRaster}.")
  }
  terra::plotRGB(mosaic, ...)
}


#' Crop or Mask a Mosaic Raster
#'
#' Crop or mask a `SpatRaster` object (`mosaic`) based on user input from an
#' interactive map or by using a provided shapefile or another raster.
#'
#' @description
#' This function allows cropping of a raster mosaic interactively or programmatically:
#'
#' - **Interactive Mode**: If neither `shapefile` nor `mosaic2` is provided, an interactive map
#'   is shown via [mosaic_view()], allowing users to draw a rectangle to define the cropping area.
#' - **Shapefile Mode**: If a `SpatVector` is provided in `shapefile`, cropping or masking is performed
#'   based on its extent or exact shape, optionally with a buffer.
#' - **Raster Mode**: If `mosaic2` is provided, `mosaic` will be cropped to match the extent of `mosaic2`.
#'
#' For disk-based mosaics, cropping with shapefiles uses GDAL (`sf::gdal_utils()`) to improve efficiency.
#'
#' @param mosaic A `SpatRaster` object to be cropped.
#' @param r,g,b,re,nir Integer indices representing the red, green, blue, red-edge, and near-infrared
#'   bands of the input mosaic. Default assumes BGR format (b = 1, g = 2, r = 3).
#' @param shapefile An optional `SpatVector` (or `sf` object) to use as cropping/masking geometry.
#'   Can be created interactively with [shapefile_input()].
#' @param mosaic2 A second `SpatRaster` whose extent will be used to crop `mosaic`.
#' @param buffer A numeric value indicating a buffer (in CRS units) to apply around the shapefile geometry.
#' @param in_memory Logical. If `TRUE`, raster processing will occur entirely in memory using `terra`.
#'   If `FALSE` (default), disk-based processing with GDAL will be used when appropriate.
#' @param show A character value indicating what to display in the interactive viewer. Either `"rgb"` or `"index"`.
#' @param index The index to show if `show = "index"`. Default is `"R"`.
#' @param max_pixels Maximum number of pixels to render in the interactive viewer.
#' @param downsample Optional downsampling factor for display purposes.
#' @param type Either `"crop"` (default) or `"mask"`:
#'   - `"crop"` crops the mosaic to the bounding box of the shapefile.
#'   - `"mask"` sets pixels outside the shapefile geometry to `NA` (recommended when using exact shapes).
#' @param ... Additional arguments passed to [mosaic_view()].
#'
#' @return A cropped or masked `SpatRaster` object.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#'   library(pliman)
#'   # Load a sample raster
#'   mosaic <- mosaic_input(system.file("ex/elev.tif", package = "terra"))
#'
#'   # Interactive cropping with drawn rectangle
#'   cropped <- mosaic_crop(mosaic)
#'
#'   # View result
#'   mosaic_view(cropped)
#' }

mosaic_crop <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        re = 4,
                        nir = 5,
                        shapefile = NULL,
                        in_memory = FALSE,
                        mosaic2 = NULL,
                        buffer = 0,
                        show = c("rgb", "index"),
                        index = "R",
                        max_pixels = 500000,
                        downsample = NULL,
                        type = c("crop", "mask"),
                        ...){
  if(is.null(shapefile) & is.null(mosaic2)){
    showopt <- c("rgb", "index")
    showopt <- showopt[pmatch(show[[1]], showopt)]
    controls <- mosaic_view(mosaic,
                            show = showopt,
                            index = index,
                            r = r,
                            g = g,
                            b = b,
                            nir = nir,
                            re = re,
                            max_pixels = max_pixels,
                            downsample = downsample,
                            edit = TRUE,
                            title = "Use the 'Draw rectangle' tool to select the cropping area.",
                            ...)

    if(!is.na(sf::st_crs(mosaic))){
      grids <-
        sf::st_make_grid(controls$finished, n = c(1, 1)) |>
        sf::st_transform(sf::st_crs(mosaic))
    } else{
      terra::crs(mosaic) <- terra::crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
      grids <-
        sf::st_make_grid(controls$finished, n = c(1, 1)) |>
        sf::st_transform(sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m"))
    }
    if(!in_memory){
      cropped <- terra::crop(mosaic, grids)
    } else{

    }
  } else{
    if(!is.null(shapefile)){
      crop_gdal <- function(mosaic, shp, exact, buffer = 0) {
        bb   <- shp |> sf::st_buffer(buffer) |> sf::st_bbox()
        opts <- c(
          "-overwrite",
          "-te",
          bb[["xmin"]], bb[["ymin"]],
          bb[["xmax"]], bb[["ymax"]]
        )
        if (exact) {
          temp_path <- tempfile(fileext = ".geojson")
          shapefile |>
            sf::st_union() |>
            sf::st_buffer(buffer) |>
            sf::st_convex_hull() |>
            sf::st_write(temp_path, driver = "GeoJSON", quiet = TRUE)
          cutline_path <- temp_path
          on.exit(unlink(temp_path), add = TRUE)
          opts <- c(
            "-cutline", cutline_path,
            "-crop_to_cutline",
            "-dstnodata", "nan",    # or a numeric nodata value like "255" or "-9999"
            "-overwrite"
          )
        }
        if(terra::inMemory(mosaic)){
          tf <- paste0(tempfile(), ".tif")
          on.exit(file.remove(tf))
          terra::writeRaster(mosaic, filename = tf)
          tempf <- tf
        } else{
          tempf <- terra::sources(mosaic)
        }
        out <- tempfile(fileext = ".tif")
        suppressWarnings(
          sf::gdal_utils(
            util        = "warp",
            source      = tempf,
            destination = out,
            options     = opts
          )
        )
        return(terra::rast(out))
      }

      if(type[[1]] == "crop"){
        if(in_memory){
          shp <- shapefile |> terra::vect() |> terra::buffer(buffer)
          cropped <- terra::crop(mosaic, shp)
        } else{
          cropped <- crop_gdal(mosaic, shapefile, exact = FALSE, buffer)
        }

      } else{
        if(in_memory){
          shp <- shapefile |> terra::vect() |> terra::buffer(buffer)
          cropped <- terra::mask(mosaic, shp)
        } else{
          cropped <- crop_gdal(mosaic, shapefile, exact = TRUE, buffer)
        }

      }

    }
    if(!is.null(mosaic2)){
      cropped <- terra::crop(mosaic, mosaic2)
    }
  }
  invisible(cropped)

}


#' Mosaic Index
#'
#' Compute or extract an index layer from a multi-band mosaic raster.
#' @inheritParams mosaic_view
#' @inheritParams image_index
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to a binary image. Use [pliman_indexes_rgb()]
#'   and [pliman_indexes_me()] to see the available RGB and multispectral
#'   indexes, respectively. Users can also calculate their own index using  `R,
#'   G, B, RE, NIR, SWIR, and TIR` bands (eg., `index = "R+B/G"`) or using the
#'   names of the mosaic's layers (ex., "(band_1 + band_2) / 2").
#' @param r,g,b,re,nir,swir,tir The red, green, blue, red-edge,  near-infrared,
#'   shortwave Infrared, and thermal infrared bands of the image, respectively.
#'   By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
#'   multispectral image is provided up to seven bands can be used to compute
#'   built-in indexes. There are no limitation of band numbers if the index is
#'   computed using the band name.
#' @param plot Plot the computed index? Defaults to `TRUE`.
#' @param in_memory Logical, indicating whether the indexes should be computed
#'   in memory. Defaults to `TRUE`. In most cases, this is 2-3 times faster, but
#'   errors can occur if `mosaic` is a large `SpatRaster`. If `FALSE`, raster
#'   algebra operations are performed on temporary files.
#' @param output Character(1), either \code{"memory"} or \code{"disk"}.
#'   If \code{"memory"}, the function returns a \code{terra::SpatRaster} object assembled in memory.
#'   If \code{"disk"}, each index layer is written out to a temporary GeoTIFF and the function returns
#'   a \code{terra::SpatRaster} object that points to those rasters. Default is \code{"memory"}.
#' @param workers numeric. The number of workers you want to use for parallel
#'   processing when computing multiple indexes.
#' @param verbose Whether to display progress messages.
#' @return An index layer extracted/computed from the mosaic raster.
#'
#' @details This function computes or extracts an index layer from the input
#'   mosaic raster based on the specified index name. If the index is not found
#'   in the package's predefined index list (see [image_index()] for more
#'   details), it attempts to compute the index using the specified band
#'   indices. The resulting index layer is returned as an `SpatRaster` object.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' names(mosaic)
#' elev2 <- mosaic_index(mosaic, "elevation * 5", plot = FALSE)
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#'
#' mosaic_plot(mosaic)
#' mosaic_plot(elev2)
#'
#' # return the original parameters
#' par(oldpar)
#' }
#'

mosaic_index <- function(mosaic,
                         index = "NGRDI",
                         r = 3,
                         g = 2,
                         b = 1,
                         re = NA,
                         nir = NA,
                         swir = NA,
                         tir = NA,
                         plot = TRUE,
                         in_memory = TRUE,
                         output = c("memory", "disk"),
                         workers = 1,
                         verbose = TRUE){

  indices <- c(r = r, g = g, b = b, re = re, nir = nir, swir = swir, tir = tir)
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
  indexname <- index
  for(i in seq_along(index)){
    if(index[i] %in% ind$Index && ind[which(index[i] == ind$Index), 3] == "HYP"){
      index[i] <- indexband_to_formula(names(mosaic), ind[which(index[i] == ind$Index), 2])
    }
  }
  # GENERAL BAND-AVAILABILITY CHECK
  # for each index expression, find which bands it uses,
  # Compute which bands are unavailable for each index
  required <- lapply(index, function(expr) {
    # pick formula: either named equation or literal expr
    formula <- if (expr %in% ind$Index) {
      as.character(ind$Equation[ind$Index == expr])
    } else {
      expr
    }
    # get positions of layers referenced in the formula
    lyrs <- unavailable_layers(indices, formula)
    names(indices)[lyrs]
  })

  # Build a list of missing-band info per index
  missing_info <- lapply(seq_along(required), function(i) {
    mb <- required[[i]][is.na(indices[ required[[i]] ])]
    if (length(mb) > 0) {
      list(name = indexname[i], bands = mb)
    } else {
      NULL
    }
  })
  missing_info <- Filter(Negate(is.null), missing_info)

  # Assemble into a data.frame
  df_missing <- if (length(missing_info) > 0) {
    data.frame(
      index         = vapply(missing_info, `[[`, "", "name"),
      missing_bands = vapply(missing_info, function(x) paste(x$bands, collapse = ", "), ""),
      stringsAsFactors = FALSE
    )
  } else {
    # no missing bands
    data.frame(index = character(0), missing_bands = character(0))
  }

  if (nrow(df_missing) > 0) {
    # build each line with an "i" bullet and the index in blue
    lines <- vapply(seq_len(nrow(df_missing)), function(i) {
      paste0(cli::col_blue(df_missing$index[i]), ": ", df_missing$missing_bands[i])
    }, character(1))

    msgs <- c(
      "The following indices were skipped due to missing bands:",
      lines
    )
    names(msgs) <- c("x", rep("i", length(lines)))

    cli::cli_warn(msgs)
  }
  # Skip any indices with missing bands
  keep <- !indexname %in% df_missing$index
  index     <- index[keep]
  indexname <- indexname[keep]
  if(length(indexname) == 0){
    cli::cli_abort("None of the requested indices can be computed because all required bands are missing.")
  }

  if(length(index) == 1){
    if(inherits(mosaic, "Image")){
      ras <- t(terra::rast(mosaic@.Data))
    } else{
      ras <- mosaic
    }
    checkind <- index %in% ind$Index
    if (!all(checkind)) {
      missing <- index[!checkind]
      cli::cli_warn(
        "Index{?es} {.val {paste(missing, collapse = \", \")}} not available. Trying to compute your own index{?es}."
      )
    }
    pattern <- "\\b\\w+\\b"
    reserved <- c("exp", "abs", "min", "max", "median", "sum", "sqrt", "cos", "sin", "tan", "log", "log10")
    layersused <- setdiff(unlist(regmatches(index, gregexpr(pattern, index, perl = TRUE))), reserved)
    onlychar <- suppressWarnings(is.na(as.numeric(layersused)))
    layers_used <- layersused[onlychar]
    if(!any(index  %in% ind$Index) & !all(layers_used  %in% c("R", "G", "B", "RE", "NIR", "SWIR", "TIR"))){
      # Extract individual layers based on the expression
      layers_used <- layers_used[is.na(suppressWarnings(as.numeric(layers_used)))]
      layers <-
        lapply(layers_used, function(x){
          mosaic[[x]]
        })
      names(layers) <- layers_used
      mosaic_gray <- eval(parse(text = index), envir = layers)
    } else{
      if(in_memory){
        if(index %in% ind$Index){
          formula <- as.character(ind$Equation[as.character(ind$Index)==index])
          lyrs <- used_layers(indices, formula)
          names(lyrs) <- names(indices)[indices %in% lyrs]
          mosaic_gray <- terra::lapp(mosaic[[lyrs]], parse_formula(formula, lyrs))
        } else{
          lyrs <- used_layers(indices, index)
          names(lyrs) <- names(indices)[indices %in% lyrs]
          mosaic_gray <- terra::lapp(mosaic[[lyrs]], parse_formula(index, lyrs))
        }
      } else{
        R <- try(ras[[indices[["r"]]]], TRUE)
        G <- try(ras[[indices[["g"]]]], TRUE)
        B <- try(ras[[indices[["b"]]]], TRUE)
        RE <- try(ras[[indices[["re"]]]], TRUE)
        NIR <- try(ras[[indices[["nir"]]]], TRUE)
        SWIR <- try(ras[[indices[["swir"]]]], TRUE)
        TIR <- try(ras[[indices[["tir"]]]], TRUE)
        if(index %in% ind$Index){
          mosaic_gray <- eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index])))
        } else{
          mosaic_gray <- eval(parse(text = index))
        }
      }
    }
    names(mosaic_gray) <- indexname
    if(!is.na(terra::crs(mosaic))){
      if(terra::crs(mosaic_gray) != terra::crs(mosaic)){
        suppressWarnings(terra::crs(mosaic_gray) <- terra::crs(mosaic))
      }
    } else{
      suppressWarnings(terra::crs(mosaic_gray) <- "+proj=utm +zone=32 +datum=WGS84 +units=m")
    }
    if(output[1] == "disk"){
      tfil <- paste0(tempfile(), ".tif")
      terra::writeRaster(mosaic_gray, tfil, overwrite = TRUE)
      mosaic_gray <-terra::rast(tfil)
    }
  } else{
    if(workers > 1){
      # CLI header + progress step
      if(output[1] == "memory"){
        cli::cli_alert_info("Argument {.arg output = {output[1]}} can only be used in a sequential strategy. Parallel processing will return the vegetation indexes in disk by default.")
      }
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Creating rasters for {.val {length(unique(index))}} indices"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
        cli::cli_progress_step(
          msg        = "Processing {.val {length(unique(index))}} indices in parallel...",
          msg_done   = "Vegetation indices successfully computed.",
          msg_failed = "Raster creation failed."
        )
      }

      # start mirai daemons
      mirai::daemons(workers)
      on.exit(mirai::daemons(0), add = TRUE)

      if(terra::inMemory(mosaic)){
        tf <- paste0(tempfile(), ".tif")
        on.exit(file.remove(tf))
        terra::writeRaster(mosaic, filename = tf)
        tempf <- tf
      } else{
        tempf <- terra::sources(mosaic)
      }

      # run mosaic_index in parallel and capture temp file paths
      tdir <- tempdir()
      tif_list <- mirai::mirai_map(
        .x = unique(index),
        .f = function(idx) {
          tfil <- paste0(tempfile(), ".tif")
          pliman::mosaic_index(
            terra::rast(tempf),
            index     = idx,
            r         = r,
            g         = g,
            b         = b,
            re        = re,
            nir       = nir,
            swir      = swir,
            tir       = tir,
            in_memory = in_memory,
            plot      = FALSE
          ) |>
            terra::writeRaster(tfil, overwrite = TRUE)
          tfil
        }
      )[.progress]

      # assemble the final mosaic
      mosaic_gray <-terra::rast(unlist(tif_list))

    } else{
      if (verbose) {
        cli::cli_rule(
          left  = cli::col_blue("Computing rasters for {.val {length(unique(index))}} indices"),
          right = cli::col_blue("Started at {.val {format(Sys.time(), '%H:%M:%OS0')}}")
        )
        cli::cli_progress_bar(
          format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
          total  = length(unique(index)),
          clear  = TRUE
        )
      }

      if(output[1] == "disk"){
        idxs    <- unique(index)
        tif_list <- vector("character", length(idxs))
        for (i in seq_along(idxs)) {
          idx  <- idxs[i]
          tfil <- tempfile(fileext = ".tif")

          # compute & write each index directly to disk
          r_i <- pliman::mosaic_index(
            mosaic,
            index     = idx,
            r         = r,
            g         = g,
            b         = b,
            re        = re,
            nir       = nir,
            swir      = swir,
            tir       = tir,
            in_memory = in_memory,
            plot      = FALSE
          )
          terra::writeRaster(r_i, filename = tfil, overwrite = TRUE)
          rm(r_i)       # drop from memory immediately
          tif_list[i] <- tfil

          if (verbose) {
            cli::cli_progress_update()
          }
        }
        # assemble your final multi-layer raster from the on-disk files
        mosaic_gray <- terra::rast(tif_list)
        names(mosaic_gray) <- unique(indexname)
      } else{

        # compute each index raster with progress updates
        res_list <- vector("list", length(unique(index)))
        for (i in seq_along(unique(index))) {
          res_list[[i]] <- mosaic_index(
            mosaic,
            index     = unique(index)[[i]],
            r         = r,
            g         = g,
            b         = b,
            re        = re,
            nir       = nir,
            swir      = swir,
            tir       = tir,
            in_memory = in_memory,
            plot      = FALSE
          )
          if (verbose) {
            cli::cli_progress_update()
          }
        }

        if (verbose) {
          cli::cli_progress_done()
        }

        # assemble into a SpatRaster
        mosaic_gray <- terra::rast(Map(c, res_list))
        names(mosaic_gray) <- unique(indexname)

      }

      if (verbose) {
        cli::cli_progress_done()
        cli::cli_rule(
          left  = cli::col_blue("{.val {length(unique(indexname))}} vegetation indices computed"),
          right = cli::col_blue("Ended at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
        )
      }


    }
  }
  if(plot){
    mosaic_plot(mosaic_gray)
  }
  invisible(mosaic_gray)
}

#' Mosaic Index with GDAL
#'
#' Compute or extract an index layer from a multi-band mosaic raster using
#' gdal_calc.py (https://gdal.org/programs/gdal_calc.html). This requires a
#' Python and GDAL installation.
#' @inheritParams mosaic_index
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
#' @param python The PATH for python.exe
#' @param gdal The PATH for gdal_calc.py
#' @return An index layer extracted/computed from the mosaic raster.
#' @export
#' @examples
#' if(interactive() & (Sys.which('python.exe') != '' ) & (Sys.which('gdal_calc.py') != '' )){
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' names(mosaic) <- "R"
#' elev2 <- mosaic_index2(mosaic, "R * 5", plot = FALSE)
#' oldpar <- par(no.readonly=TRUE)
#' mosaic_plot(mosaic)
#' mosaic_plot(elev2)
#' par(mfrow=c(1,2))
#' }

mosaic_index2 <- function(mosaic,
                          index = "B",
                          r = 3,
                          g = 2,
                          b = 1,
                          re = 4,
                          nir = 5,
                          plot = TRUE,
                          python = Sys.which('python.exe'),
                          gdal = Sys.which('gdal_calc.py')) {
  if(python == ''){
    cli::cli_abort(c(
      "!" = "Python executable {.file python.exe} not found on the system.",
      "i" = "Please install Python and ensure it is added to the system {.envvar PATH} variable."
    ))
  }
  if(gdal==''){
    cli::cli_abort(c(
      "!" = "{.file gdal_calc.py} not found on the system.",
      "i" = "Make sure GDAL is installed and available in the system {.envvar PATH}.",
      ">" = "You can install GDAL via Miniconda: {.code conda install -c conda-forge gdal}",
      "i" = "Download Miniconda at: {.url https://docs.conda.io/projects/miniconda/en/latest/}",
      "i" = "Also ensure that the Python Scripts directory (e.g., {.path C:/Users/yourname/AppData/Local/Programs/Python/PythonXX/Scripts}) is in the {.envvar PATH}."
    ))

  }
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
  checkind <- index %in% ind$Index
  if (!checkind) {
    cli::cli_inform("Index {.val {paste0(index[!checkind], collapse = ', ')}} is not available. Trying to compute your own index.")

  } else{
    index <- as.character(ind$Equation[as.character(ind$Index)==index])
  }
  if(terra::inMemory(mosaic)){
    tf <- tempfile(fileext = ".tif")
    mosaic_export(mosaic, tf)
    on.exit(file.remove(tf))
    infile <- tf
  } else{
    infile <- terra::sources(mosaic)
  }
  mosaicbands <- c("R", "G", "B", "E", "I")
  outfile <- tempfile(fileext = ".tif")
  nbands <- terra::nlyr(mosaic)
  inputs <- paste0('-',
                   mosaicbands[seq_len(nbands)], ' ', infile, ' --',
                   mosaicbands[seq_len(nbands)], '_band ', seq_len(nbands), collapse=' ')
  system2(python,
          args=c(gdal,
                 inputs,
                 sprintf("--outfile=%s", outfile),
                 sprintf('--calc="%s"', index),
                 '--co="COMPRESS=DEFLATE"',
                 '--co="BIGTIFF=IF_NEEDED"',
                 '--overwrite'),
          stdout=FALSE
  )
  res <- terra::rast(outfile)
  if(plot){
    terra::plot(res)
  }
  return(res)
}

#' Segment a mosaic
#'
#' Segment a `SpatRaster` using a computed image index. By default, values
#' greater than `threshold` are kept in the mask.
#'
#' @inheritParams mosaic_index
#' @inheritParams mosaic_analyze
#' @param return The output of the function. Either 'mosaic' (the segmented
#'   mosaic), or 'mask' (the binary mask).
#'
#' @return The segmented mosaic (`SpatRaster` object)
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' seg <-
#' mosaic_segment(mosaic,
#'                index = "elevation",
#'                threshold = 350)
#' mosaic_plot(seg)
#' }
mosaic_segment <- function(mosaic,
                           index = "R",
                           r = 3,
                           g = 2,
                           b = 1,
                           re = NA,
                           nir = NA,
                           swir = NA,
                           tir = NA,
                           threshold = "Otsu",
                           invert = FALSE,
                           return = c("mosaic", "mask")){
  if(!return[[1]] %in% c("mosaic", "mask")){
    cli::cli_abort(c(
      "!" = "{.arg return} must be one of:",
      "x" = "{.val mosaic}, {.val mask}"
    ))

  }
  ind <- mosaic_index(mosaic,
                      index = index,
                      r = r,
                      g = g,
                      b = b,
                      re = re,
                      nir = nir,
                      swir = swir,
                      tir = tir,
                      plot = FALSE)
  thresh <- ifelse(threshold == "Otsu", otsu(na.omit(terra::values(ind)[, index])), threshold)
  if(invert){
    mask <- ind[[index]] > thresh
  } else{
    mask <- ind[[index]] < thresh
  }
  if(return[[1]] == 'mosaic'){
    terra::mask(mosaic, mask, maskvalue = TRUE)
  } else{
    mask
  }
}


#' Segments a mosaic interactively
#'
#' The function segments a mosaic using an interative process where the user
#' picks samples from background (eg., soil) and foreground (eg., plants).
#'
#' @inheritParams mosaic_index
#' @inheritParams mosaic_view
#' @param basemap An optional `mapview` object.
#' @param return The output of the function. Either 'mosaic' (the segmented
#'   mosaic), or 'mask' (the binary mask).
#'
#' @return An `SpatRaster` object with the segmented `mosaic` (if `return =
#'   'mosaic'`) or a mask (if `return = 'mask'`).
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#'  mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'  seg <- mosaic_segment_pick(mosaic)
#'  mosaic_plot(seg)
#' }
mosaic_segment_pick <- function(mosaic,
                                basemap = NULL,
                                g = 2,
                                r = 3,
                                b = 1,
                                max_pixels = 2e6,
                                downsample = NULL,
                                quantiles = c(0, 1),
                                return = c("mosaic", "mask")){
  if(!return[[1]] %in% c("mosaic", "mask")){
    cli::cli_abort(c(
      "!" = "{.arg return} must be one of:",
      "v" = "{.val mosaic} or {.val mask}"
    ))
  }
  downsample <- ifelse(is.null(downsample), find_aggrfact(mosaic, max_pixels = max_pixels), downsample)
  if(downsample > 0){
    mosaic <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
  }
  if(is.null(basemap)){
    basemap <-
      mosaic_view(mosaic,
                  r = r,
                  g = g,
                  b = b,
                  max_pixels = max_pixels,
                  downsample = downsample,
                  quantiles = quantiles)
  }
  soil <- mapedit::editMap(basemap,
                           title = "Use the 'Draw Rectangle' tool to pick up background fractions",
                           editor = "leafpm")$finished
  soil <- soil |> sf::st_transform(sf::st_crs(mosaic))
  soil_sample <-
    exactextractr::exact_extract(mosaic, soil, progress = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::select(-coverage_fraction) |>
    dplyr::mutate(class = 0)

  mapview::mapview() |> mapedit::editMap()

  plant <- mapedit::editMap(basemap,
                            title = "Use the 'Draw Rectangle' tool to pick up foreground fractions",
                            editor = "leafpm")$finished
  plant <- plant |> sf::st_transform(sf::st_crs(mosaic))
  plant_sample <-
    exactextractr::exact_extract(mosaic, plant, progress = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::select(-coverage_fraction) |>
    dplyr::mutate(class = 1)
  df_train <- dplyr::bind_rows(plant_sample, soil_sample)
  if(ncol(df_train) == 2){
    names(df_train)[[1]] <- names(mosaic)
  }
  mod <- suppressWarnings(
    glm(class ~.,
        data = df_train,
        family = binomial("logit"))
  )
  mask <- terra::predict(mosaic, mod, type = "response")
  mask[mask < 0.5] <- 0
  mask[mask > 0.5] <- 1
  if(return[[1]] == 'mosaic'){
    terra::mask(mosaic, mask, maskvalue = TRUE)
  } else{
    mask
  }
}


#' Mosaic to pliman
#'
#' Convert an `SpatRaster` object to a `Image` object with optional scaling.
#' @inheritParams mosaic_view
#' @inheritParams mosaic_index
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
#' @param rescale Rescale the final values? If `TRUE` the final values are
#'   rescaled so that the maximum value is 1.
#' @param coef An addition coefficient applied to the resulting object. This is
#'   useful to adjust the brightness of the final image. Defaults to 0.
#'
#' @return An `Image` object with the same number of layers as `mosaic`.
#'
#' @details This function converts `SpatRaster` into an `Image` object, which
#'   can be used for image analysis in `pliman`. Note that if a large
#'   `SpatRaster` is loaded, the resulting object may increase considerably the
#'   memory usage.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' # Convert a mosaic raster to an Image object
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' pliman_image <- mosaic_to_pliman(mosaic)
#' plot(pliman_image)
#' }
#'
mosaic_to_pliman <- function(mosaic,
                             r = 3,
                             g = 2,
                             b = 1,
                             re = 4,
                             nir = 5,
                             rescale =  TRUE,
                             coef = 0){
  if(class(mosaic) %in% c("RasterStack","RasterLayer","RasterBrick")){
    mosaic <- terra::rast(mosaic)
  }
  nlr <- terra::nlyr(mosaic)
  if(nlr == 5){
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))[,, c(r, g, b, re, nir)]
  } else if(nlr == 3){
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))[,, c(r, g, b)]
  } else{
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))
  }
  if(isTRUE(rescale)){
    mosaic <- mosaic / max(mosaic, na.rm = TRUE)
  }
  if(nlr == 3){
    EBImage::colorMode(mosaic) <- "color"
  }
  return(mosaic + coef)
}

#' Mosaic to RGB
#'
#' Convert an `SpatRaster` to a three-band RGB image of class `Image`.
#'
#' @inheritParams mosaic_to_pliman
#' @param plot Logical, whether to display the resulting RGB image (default:
#'   TRUE).
#' @param r,g,b The red, green, blue bands.
#' @return A three-band RGB image represented as a pliman (EBImage) object.
#'
#' @details This function converts `SpatRaster` that contains the RGB bands into
#'   a three-band RGB image using pliman (EBImage). It allows you to specify the
#'   band indices for the red, green, and blue channels, as well as apply a
#'   scaling coefficient to the final image. By default, the resulting RGB image
#'   is displayed, but this behavior can be controlled using the `plot`
#'   parameter.
#'
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#'
#' library(pliman)
#' # Convert a mosaic raster to an RGB image and display it
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # Convert a mosaic raster to an RGB image without displaying it
#' rgb_image <- mosaic_to_rgb(c(mosaic * 2, mosaic - 0.3, mosaic * 0.8))
#' plot(rgb_image)
#' }
#'
#'
mosaic_to_rgb <- function(mosaic,
                          r = 3,
                          g = 2,
                          b = 1,
                          coef = 0,
                          plot = TRUE){
  ebim <- mosaic_to_pliman(mosaic,
                           r = r,
                           g = g,
                           b = b,
                           coef = coef)[,,c(r, g, b)]
  EBImage::colorMode(ebim) <- "color"
  invisible(ebim)
}


#' Prepare a mosaic
#'
#' Prepare an `SpatRaster` object to be analyzed in pliman. This includes
#' cropping the original mosaic, aligning it, and cropping the aligned object.
#' The resulting object is an object of class `Image` that can be further
#' analyzed.
#' @inheritParams mosaic_view
#' @inheritParams mosaic_index
#' @inheritParams mosaic_to_pliman
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
#' @param crop_mosaic Logical, whether to crop the mosaic interactively before
#'   aligning it (default: FALSE).
#' @param align Logical, whether to align the mosaic interactively (default:
#'   TRUE).
#' @param crop_aligned Logical, whether to crop the aligned mosaic interactively
#'   (default: TRUE).
#'
#' @return A prepared object of class `Image`.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_prepare(mosaic)
#' }
#'
mosaic_prepare <- function(mosaic,
                           r = 3,
                           g = 2,
                           b = 1,
                           re = 4,
                           nir = 5,
                           crop_mosaic = TRUE,
                           align = TRUE,
                           crop_aligned = TRUE,
                           rescale =  TRUE,
                           coef = 0,
                           viewer = "mapview",
                           max_pixels = 500000,
                           show = "rgb",
                           index = "R"){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(isTRUE(crop_mosaic)){
    cropped <- mosaic_crop(mosaic,
                           show = show,
                           index = index,
                           max_pixels = max_pixels,
                           r = r,
                           g = g,
                           b = b,
                           nir = nir,
                           re = re)
    ebimg <- mosaic_to_pliman(cropped,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              rescale = rescale,
                              coef = coef)
    if(vieweropt != "base"){
      image_view(ebimg[1:5, 1:5,], edit = TRUE)
    }
  } else{
    ebimg <- mosaic_to_pliman(mosaic,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              rescale = rescale,
                              coef = coef)
  }
  if(isTRUE(align)){
    aligned <- image_align(ebimg, viewer = vieweropt)
    if(vieweropt != "base"){
      image_view(aligned[1:5, 1:5,], edit = TRUE)
    }
  } else{
    aligned <- ebimg
  }

  if(isTRUE(crop_aligned)){
    cropped <- image_crop(aligned, viewer = vieweropt)
  } else{
    cropped <- aligned
  }
  if(dim(cropped)[3] == 3){
    EBImage::colorMode(cropped) <- "color"
  }
  invisible(cropped)
}




#' Drawing Lines or Polygons with Raster Information
#'
#' @details
#' The `mosaic_draw` function enables you to create mosaic drawings from
#' remote sensing data and compute vegetation indices.
#'
#'  * If a line is drawn using the "Draw Polyline" tool, the profile of `index` is
#'  displayed on the y-axis along the line's distance, represented in meter
#'  units. It is important to ensure that the Coordinate Reference System (CRS)
#'  of `mosaic` has latitude/longitude units for accurate distance
#'  representation.
#'
#'  * If a rectangle or polygon is drawn using the "Draw Rectangle" or "Draw Polygon"
#'  tools, the `index` values are calculated for each object. By default, the
#'  raw data is returned. You can set the `summarize_fun` to compute a summary
#'  statistic for each object.
#'
#' @inheritParams mosaic_view
#' @inheritParams mosaic_index
#' @inheritParams analyze_objects
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
#' @param color_regions The color palette for displaying index values. Defaults
#'   to `rev(grDevices::terrain.colors(50))`.
#' @param threshold By default (threshold = "Otsu"), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is informed, this value will be used as a threshold.
#' @param invert Inverts the mask if desired. Defaults to `FALSE`.
#' @param segment Should the raster object be segmented? If set to `TRUE`,
#'   pixels within each polygon/rectangle will be segmented based on the
#'   `threshold` argument.
#' @param summarize_fun An optional function or character vector. When
#'   `summarize_fun = "mean"`, the mean values of `index` are calculated within
#'   each object. For more details on available functions, refer to
#'  [exactextractr::exact_extract()].
#' @param buffer Adds a buffer around the geometries of the SpatVector created.
#'   Note that the distance unit of `buffer` will vary according to the CRS of
#'   `mosaic`.
#' @param plot Plots the draw line/rectangle? Defaults to `TRUE`.
#' @param plot_layout The de plot layout. Defaults to `plot_layout = c(1, 2, 3,
#'   3)`. Ie., the first row has two plots, and the second row has one plot.
#' @importFrom terra crop vect extract
#' @importFrom exactextractr exact_extract
#' @importFrom graphics layout
#' @importFrom stats smooth
#' @return An invisible list containing the mosaic, draw_data, distance,
#'   distance_profile, geometry, and map.
#' @export
#' @examples
#'
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' # Load a raster showing the elevation of Luxembourg
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # draw a polyline to see the elevation profile along the line
#' mosaic_draw(mosaic, buffer = 1500)
#' }
mosaic_draw <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        re = 4,
                        nir = 5,
                        index = "NGRDI",
                        show = "rgb",
                        segment = FALSE,
                        viewer = c("mapview", "base"),
                        threshold = "Otsu",
                        invert = FALSE,
                        summarize_fun = NULL,
                        buffer = 2,
                        color_regions = rev(grDevices::terrain.colors(50)),
                        alpha = 1,
                        max_pixels = 1000000,
                        downsample = NULL,
                        quantiles = c(0, 1),
                        plot = TRUE,
                        plot_layout = c(1, 2, 3, 3)){
  if(is.null(terra::crs(mosaic))){
    terra::crs(mosaic) <- "+proj=utm +zone=32 +datum=WGS84 +units=m"
  }
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]

  points <- mosaic_view(mosaic,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        max_pixels = max_pixels,
                        downsample = downsample,
                        alpha = alpha,
                        quantiles = quantiles,
                        index = index[[1]],
                        show = show,
                        edit = TRUE)

  nlyrs <- terra::nlyr(mosaic)
  polygons <- points$finished$geometry
  polygons_spv <- sf::st_transform(polygons, crs = sf::st_crs(mosaic))
  polygons_ext <- terra::vect(polygons_spv)
  ext <- terra::buffer(polygons_ext, buffer) |> terra::ext()
  mosaiccr <- terra::crop(mosaic, ext)

  # Compute the image indexes
  if(nlyrs > 1){
    mind <- terra::rast(
      Map(c,
          lapply(seq_along(index), function(i){
            mosaic_index(mosaiccr,
                         index = index[[i]],
                         r = r,
                         g = g,
                         b = b,
                         re = re,
                         nir = nir,
                         plot = FALSE)
          })
      )
    )
  } else{
    index <- names(mosaiccr)
    mind <- mosaiccr
  }

  if(inherits(polygons, "sfc_LINESTRING")){
    vals <-
      terra::extractAlong(x = mind,
                          y = polygons_ext,
                          ID = FALSE)
    coords <- as.matrix(polygons_spv[[1]])
    n <- nrow(coords)
    distances <- NULL
    for (j in 1:(n - 1)) {
      x1 <- coords[j, 1]
      y1 <- coords[j, 2]
      x2 <- coords[j + 1, 1]
      y2 <- coords[j + 1, 2]
      distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
      distances[j] <- distance
    }
    # distances
    dists <- cumsum(distances)
    dist <- max(dists)

    if(plot){
      if(nlyrs > 2){
        layout(
          matrix(plot_layout, nrow = 2, byrow = TRUE),
          heights = c(3, 3)
        )
        on.exit(layout(1))
        scl <- max(terra::minmax(mosaiccr))

        terra::plotRGB(mosaiccr,
                       r = r,
                       g = g,
                       b = b,
                       scale = ifelse(scl < 255, 255, scl),
                       colNA = "#00000000",
                       mar = c(2, 2, 2, 2),
                       axes = FALSE)
        lines(coords,
              col = "red",
              lwd = 3)

        terra::plot(mind[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)
        lines(coords,
              col = "red",
              lwd = 3)

        plot(x = seq(0, dist, length.out = nrow(vals)),
             y = smooth(vals[, index[[1]]]),
             xlab = "Distance",
             ylab = index[[1]],
             # mar = c(0, 0, 0, 0),
             type = "l",
             xlim = c(0, dist),
             col = "red")
      } else{
        layout(
          matrix(plot_layout, nrow = 2, byrow = TRUE),
          heights = c(3, 3)
        )
        on.exit(layout(1))
        ext2 <- terra::buffer(polygons_ext, buffer * 3) |> terra::ext()
        mosaiccr2 <- terra::crop(mosaic[[1]], ext2)
        terra::plot(mosaiccr2[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)


        terra::plot(mind[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)
        lines(coords,
              col = "red",
              lwd = 3)

        plot(x = seq(0, dist, length.out = nrow(vals)),
             y = smooth(vals[, index[[1]]]),
             xlab = "Distance",
             ylab = index[[1]],
             type = "l",
             xlim = c(0, dist),
             col = "red")

      }
    }
    map <- NULL
  } else{
    if(segment){
      if(invert){
        mask <- mind[[1]] < otsu(na.omit(terra::values(mind)[, index[1]]))
      } else{
        mask <- mind[[1]] < otsu(na.omit(terra::values(mind)[, index[1]]))
      }
      mind <- terra::mask(mind, mask, maskvalues = TRUE)
    }
    mind <- terra::mask(mind, polygons_ext)
    vals <-
      suppressWarnings(
        exactextractr::exact_extract(x = mind,
                                     y = polygons,
                                     fun = summarize_fun,
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))
      )
    if(inherits(vals, "list")){
      vals <-
        do.call(rbind, lapply(1:length(vals), function(i){
          tmp <- transform(vals[[i]], id = i)
          tmp[, c(ncol(tmp), 1:(ncol(tmp) - 2))]
        }
        ))
    } else{
      vals <- transform(vals, id = paste0(1:nrow(vals)))
      vals <- vals[, c(ncol(vals), 1:(ncol(vals) - 2))]
    }
    colnames(vals) <- c("id", index)
    vals <- na.omit(vals)

    if(nlyrs > 2){
      if(vieweropt == "base"){
        terra::plot(mind, axes = FALSE)
        map <- NULL
      } else{

        map <-
          mapview::viewRGB(as(mosaiccr, "Raster"),
                           na.color = "#00000000",
                           layer.name = "base",
                           r = r,
                           g = g,
                           b = b,
                           maxpixels = 10000000,
                           quantiles = quantiles)

        for (i in 1:length(index)) {
          map <-
            mapview::mapview(mind[[i]],
                             hide = TRUE,
                             map = map,
                             layer.name = index[[i]],
                             map.types = mapview::mapviewGetOption("basemaps"),
                             maxpixels =  max_pixels,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE)
        }
        map

      }
    } else{
      if(vieweropt == "base"){
        terra::plot(mind, axes = FALSE)
        map <- NULL
      } else{
        map <-
          mapview::mapview(mind,
                           layer.name = names(mind),
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
      }
    }

    dists <- NULL
    dist <- NULL
  }

  invisible(list(
    mosaic = mosaic,
    draw_data = vals,
    distance = dist,
    distance_profile = dists,
    geometry = points,
    map = map

  ))
}


#' Convert Sentinel data to GeoTIFF format
#'
#' This function converts Sentinel satellite data files to GeoTIFF format.
#'
#' @param layers (character) Vector of file paths to Sentinel data files. If
#'   NULL, the function searches for files in the specified path with names
#'   containing "B".
#' @param path (character) Directory path where Sentinel data files are located.
#'   Default is the current directory.
#' @param destination (character) File path for the output GeoTIFF file.
#' @param spat_res (numeric) Spatial resolution of the output GeoTIFF file.
#'   Default is 10 meters.
#'
#' @details The function converts Sentinel satellite data files to GeoTIFF
#'   format using GDAL utilities. It builds a virtual raster file (VRT) from the
#'   input files and then translates it to GeoTIFF format. Compression is
#'   applied to the output GeoTIFF file using DEFLATE method.
#'
#'
#'
#' @export
sentinel_to_tif <- function(layers = NULL,
                            path = ".",
                            destination,
                            spat_res = 10){

  if(is.null(layers)){
    files <- list.files(path = path, pattern = "B")
  } else{
    files <- layers
  }
  tf <- tempfile(fileext = ".vrt")
  sf::gdal_utils(
    util = "buildvrt",
    source  = files,
    destination = tf,
    options = strsplit(paste("-separate -tr", spat_res, spat_res), split = "\\s")[[1]]
  )

  sf::gdal_utils(
    util = "translate",
    source  = tf,
    destination = paste0(path, "/", destination),
    options = strsplit(paste("-co COMPRESS=DEFLATE", "-of GTiff"), split = "\\s")[[1]]
  )
}

# helper to run code and capture all cat() into a temp file
capture_to_tempfile <- function(expr) {
  tmp <- tempfile(pattern = "fields_log_", fileext = ".txt")
  con <- file(tmp, open = "wt")
  sink(con)            # divert stdout (where cat() goes)
  on.exit({            # ensure we restore and close
    sink()
    close(con)
  })
  result <- eval(expr) # run your code
  invisible(list(fit = result, logfile = tmp))
}

a <- capture_to_tempfile({
  cat("a")
  3
})
#' Calculate Canopy Height Model and Volume
#'
#' This function calculates the canopy height model (CHM) and the volume for a
#' given digital surface model (DSM) raster layer. Optionally, a digital terrain
#' model (DTM) can be provided or interpolated using a set of points or a moving
#' window.
#'
#' @param dsm A `SpatRaster` object representing the digital surface model. Must
#'   be a single-layer raster.
#' @param dtm (optional) A `SpatRaster` object representing the digital terrain
#'   model. Must be a single-layer raster. If not provided, it can be
#'   interpolated from points or created using a moving window.
#' @param points (optional) An `sf` object representing sample points for DTM
#'   interpolation. If provided, `dtm` will be interpolated using these points.
#' @param interpolation (optional) A character string specifying the
#'   interpolation method to use when `points` are provided. Options are
#'   "Kriging" (default) or "Tps" (Thin Plate Spline).
#' @param window_size An integer  (meters) specifying the window size (rows and
#'   columns, respectively) for creating a DTM using a moving window. Default is
#'   c(10, 10).
#' @param ground_quantile Numeric value between `0` and `1` indicating the
#'   quantile threshold for ground point selection in the CHM computation. Lower
#'   values (e.g., `0`) retain the lowest ground points, while higher values
#'   (e.g., `1`) consider higher ground elevations. Default is `0`, which uses
#'   the lowest points within each window.
#' @param mask (optional) A `SpatRaster` object used to mask the CHM and volume
#'   results. Default is NULL.
#' @param mask_soil Is `mask` representing a soil mask (eg., removing plants)? Default is TRUE.
#' @param verbose Return the progress messages. Default is TRUE.
#'
#' @return A `SpatRaster` object with three layers: `dtm` (digital terrain
#'   model), `height` (canopy height model), and `volume`.
#'
#' @details
#' The function first checks if the input `dsm` is a valid single-layer
#' `SpatRaster` object. If `dtm` is not provided, The function generates a
#' Digital Terrain Model (DTM) from a Digital Surface Model (DSM) by
#' downsampling and smoothing the input raster data. It iterates over the DSM
#' matrix in windows of specified size, finds the minimum value within each
#' window, and assigns these values to a downsampled matrix. After downsampling,
#' the function applies a mean filter to smooth the matrix, enhancing the visual
#' and analytical quality of the DTM. Afterwards, DTM is resampled with the
#' original DSM.
#'
#' If both `dsm` and `dtm` are provided, the function ensures they have the same
#' extent and number of cells, resampling `dtm` if necessary. The CHM is then
#' calculated as the difference between `dsm` and `dtm`, and the volume is
#' calculated by multiplying the CHM by the pixel size. The results are
#' optionally masked using the provided `mask`.
#'
#' @export

mosaic_chm <- function(dsm,
                       dtm = NULL,
                       points = NULL,
                       interpolation = c("Tps", "Kriging"),
                       window_size = c(5, 5),
                       ground_quantile = 0,
                       mask = NULL,
                       mask_soil = TRUE,
                       verbose = TRUE){
  if(!interpolation[[1]] %in% c("Tps", "Kriging")){
    cli::cli_abort("{.arg interpolation} must be one of {.val Tps} or {.val Kriging}.")
  }
  # Check if fields is installed
  check_and_install_package("fields")
  sampp <- NULL
  ch1 <- !inherits(dsm,"SpatRaster") || !terra::nlyr(dsm) == 1 || terra::is.bool(dsm) || is.list(dsm)
  if(ch1){
    cli::cli_abort("{.arg dsm} must be single-layer {.code SpatRaster} objects")
  }
  # mensagens CLI
  if (verbose) {
    cli::cli_rule(
      left  = cli::col_blue("Canopy Height-Model generation"),
      right = cli::col_blue("{.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
    )
  }
  # interpolate dtm using sample of points
  if(is.null(dtm) & !is.null(points)){
    # sampling points
    points <- points |> sf::st_transform(sf::st_crs(dsm))


    if(verbose){
      cli::cli_progress_step(
        msg        = "Extracting values...",
        msg_done   = "Extracting values",
        msg_failed = "Extracting values failed"
      )
    }
    vals <- terra::extract(dsm, terra::vect(points), xy = TRUE)
    xy <- cbind(vals$x,vals$y)
    z <- vals[, 2]

    if(verbose){
      cli::cli_progress_step(
        msg        = "Interpolating the raster...",
        msg_done   = "Interpolating the raster",
        msg_failed = "Interpolation failed"
      )
    }
    if(interpolation[[1]] == "Kriging"){
      fit <- fields::Krig(xy, z, aRange=20, give.warnings = FALSE)
    }
    if(interpolation[[1]] == "Tps"){
      fit <- fields::Tps(xy, z, give.warnings = FALSE)
    }

    sampp <- NULL
    # low resolution to interpolate
    aggr <- find_aggrfact(dsm, 4e5)
    if(aggr > 0){
      mosaicintp <- mosaic_aggregate(dsm, round(100 / aggr))
    } else{
      mosaicintp <- dsm
    }
    dtm <- terra::interpolate(terra::rast(mosaicintp), fit)
    gc()
    terra::crs(dtm) <- terra::crs(dsm)
    if(verbose){
      cli::cli_progress_step(
        msg        = "Resampling and masking the interpolated raster...",
        msg_done   = "Resampling and masking the interpolated raster",
        msg_failed = "Resampling failed"
      )
    }
    dtm <- terra::resample(dtm, dsm)
    gc()
    dtm <- terra::mask(dtm, dsm)
  }
  # create a dtm using a moving window that extract the minimum values
  # from dsm
  if(is.null(dtm) & is.null(points)){
    resolu <- terra::res(dsm)
    extens <- terra::ext(dsm)
    wide <- extens[2] - extens[1]
    heig <- extens[4] - extens[3]
    nr <- ceiling( heig / window_size[[1]])
    nc <- ceiling( wide / window_size[[2]])
    if(verbose){
      cli::cli_progress_step(
        msg        = "Extracting ground points for each moving window...",
        msg_done   = "Extracting ground points for each moving window",
        msg_failed = "Extraction failed"
      )
    }
    shp <- shapefile_build(dsm,
                           nrow = nr,
                           ncol = nc,
                           build_shapefile = FALSE,
                           verbose = FALSE)
    vals <- exactextractr::exact_extract(dsm,
                                         shp[[1]],
                                         fun = "quantile",
                                         quantiles = ground_quantile,
                                         progress = FALSE)
    gc()
    cent <- suppressWarnings(sf::st_centroid(shp[[1]]))
    sampp <-
      cent |>
      dplyr::mutate(dtm = vals) |>
      dplyr::filter(!is.na(dtm))

    xy <- sf::st_coordinates(sampp)
    z <- sampp$dtm

    if(verbose){
      cli::cli_progress_step(
        msg        = "Interpolating ground points...",
        msg_done   = "Interpolating ground points",
        msg_failed = "Interpolation failed"
      )
    }
    if(interpolation[[1]] == "Kriging"){
      fit <- fields::Krig(xy, z, aRange=20, give.warnings = FALSE)
    }
    if(interpolation[[1]] == "Tps"){
      fit <- fields::Tps(xy, z, give.warnings = FALSE)
    }

    # low resolution to interpolate
    aggr <- find_aggrfact(dsm, 4e5)
    if(aggr > 0){
      mosaicintp <- mosaic_aggregate(dsm, round(100 / aggr))
    } else{
      mosaicintp <- dsm
    }
    dtm <- terra::interpolate(terra::rast(mosaicintp), fit)
    terra::crs(dtm) <- terra::crs(dsm)
    if(verbose){
      cli::cli_progress_step(
        msg        = "Resampling and masking the interpolated raster...",
        msg_done   = "Resampling and masking the interpolated raster",
        msg_failed = "Resampling failed"
      )
    }
    dtm <- terra::resample(dtm, dsm)
    gc()
    dtm <- terra::mask(dtm, dsm)
    gc()
  }
  # now, create a chm
  if(!is.null(dtm)){
    ch2 <- !inherits(dtm,"SpatRaster") || !terra::nlyr(dtm) == 1 || terra::is.bool(dtm) || is.list(dtm)
    if(ch2){
      cli::cli_abort("{.arg dtm} must be single-layer SpatRaster objects")
    }
    if((terra::ext(dsm) != terra::ext(dtm)) || (terra::ncell(dtm) != terra::ncell(dsm))){
      if(verbose){
        cli::cli_progress_step(
          msg        = "Putting dtm and dsm in the same resolution...",
          msg_done   = "Putting dtm and dsm in the same resolution",
          msg_failed = "Failed"
        )
      }
      dtm <- terra::resample(dtm, dsm)
    }
    if(verbose){
      cli::cli_progress_step(
        msg        = "Building the canopy height model...",
        msg_done   = "Building the canopy height model",
        msg_failed = "Failed"
      )
    }
    chm <- dsm - dtm
    gc()
    if(!is.null(mask)){
      if((terra::ext(mask) != terra::ext(dsm)) || (terra::ncell(mask) != terra::ncell(dsm))){
        mask <- terra::resample(mask, dsm)
      }
      chm <- terra::mask(chm,  mask, maskvalues = mask_soil)
    }
    chm <- c(dtm, chm)
    names(chm) <- c("dtm", "height")
  }
  return(list(chm = chm,
              sampling_points = sampp,
              mask = ifelse(is.null(mask), FALSE, TRUE),
              res =  terra::res(dsm)))
}

#' Extracts height metrics and plot quality from a Canopy Height Model (CHM)
#'
#' This function extracts height-related summary statistics from a CHM using a
#' given shapefile.
#'
#' @param chm An object computed with [mosaic_chm()].
#' @param shapefile An `sf` object containing the polygons over which height
#'   metrics are extracted.
#' @param chm_threshold A numeric value representing the height threshold for
#'   calculating coverage. If `NULL`, coverage is not computed.
#' @return An `sf` object containing height summary statistics for each plot,
#'   including:
#' * `min`: Minimum height value.
#' * `q05`: 5th percentile height value.
#' * `q50`: Median height value.
#' * `q95`: 95th percentile height value.
#' * `max`: Maximum height value.
#' * `mean`: Mean height value.
#' * `volume`: Total sum of heights multiplied by CHM resolution.
#' * `coverage`: If a mask is used in [mosaic_chm()] or `chm_threshold` is
#' informed, returns the proportion of pixels covered within the plot.
#' Otherwise, returns 1.
#'
#' @export

mosaic_chm_extract <- function(chm, shapefile, chm_threshold = NULL) {
  custom_summary <- function(values, coverage_fractions, ...) {
    valids <- na.omit(values)
    sumvalids <- sum(valids)
    quantiles <- quantile(valids, c(0, 0.05, 0.5, 0.95, 1))
    mean_val <- sumvalids / length(valids)
    volume <- sumvalids * prod(chm[["res"]])
    cv <- mean(valids) / sd(valids)
    entropy <- entropy(valids)
    if (!is.null(chm_threshold)) {
      coverage <- sum(valids > chm_threshold) / length(valids)
      data.frame(
        min = quantiles[[1]],
        q05 = quantiles[[2]],
        q50 = quantiles[[3]],
        q95 = quantiles[[4]],
        max = quantiles[[5]],
        mean = mean_val,
        cv = cv,
        entropy = entropy,
        volume = volume,
        coverage = coverage
      )
    } else{
      data.frame(
        min = quantiles[[1]],
        q05 = quantiles[[2]],
        q50 = quantiles[[3]],
        q95 = quantiles[[4]],
        max = quantiles[[5]],
        mean = mean_val,
        cv = cv,
        entropy = entropy,
        volume = volume
      )
    }

  }
  height <- exactextractr::exact_extract(chm$chm[[2]],
                                         shapefile,
                                         fun = custom_summary,
                                         force_df = TRUE,
                                         progress = FALSE)
  if (chm$mask) {
    area2 <- exactextractr::exact_extract(chm$chm[[2]],
                                          shapefile,
                                          coverage_area = TRUE,
                                          force_df = TRUE,
                                          progress = FALSE)
    covered_area <- purrr::map_dfr(area2, function(x) {
      data.frame(covered_area = sum(na.omit(x)[, "coverage_area"]),
                 plot_area = sum(x[, "coverage_area"]))
    }) |>
      dplyr::mutate(coverage = covered_area / plot_area)
  } else {
    area <- as.numeric(sf::st_area(shapefile))
    if (is.null(chm_threshold)) {
      covered_area <- data.frame(plot_area = area)
    } else {
      covered_area <- data.frame(covered_area = area * height$coverage,
                                 plot_area = area)
    }
  }
  shapefile <-
    shapefile |>
    dplyr::select(-suppressWarnings(dplyr::any_of(c("x", "y"))))
  centroids <- suppressWarnings(sf::st_centroid(shapefile)) |> sf::st_coordinates()
  colnames(centroids) <- c("x", "y")
  dftmp <-
    dplyr::bind_cols(height, covered_area, centroids, shapefile) |>
    sf::st_as_sf() |>
    dplyr::relocate(unique_id, block, plot_id, row, column, x, y, .before = 1)

  return(dftmp)
}




#' Apply a height mask to CHM data
#'
#' This function applies a height-based mask to a Canopy Height Model (CHM),
#' focusing on areas with heights above a specified `lower` threshold and,
#' optionally, below an `upper` threshold.
#'
#' The `mosaic_chm` function, used internally, generates the DTM from the DSM by
#' downsampling and smoothing raster data, applying a moving window to extract
#' minimum values and then interpolating the results. The CHM is computed as the
#' height difference between the DSM and DTM. This function calculates and
#' applies a mask based on height thresholds.
#'
#' @inheritParams mosaic_chm
#' @param lower A numeric value specifying the lower height threshold. All
#'   heights greater than this value are retained.
#' @param upper An optional numeric value specifying the upper height threshold.
#'   If provided, only heights between `lower` and `upper` are retained.
#'
#' @return An `SpatRaster` object representing the masked CHM.
#' @export
#'
mosaic_chm_mask <- function(dsm,
                            lower,
                            upper = NULL,
                            window_size = c(5, 5),
                            interpolation = "Tps"){
  chm <- mosaic_chm(dsm,
                    window_size = window_size,
                    verbose = FALSE)
  if(is.null(upper)){
    chm$chm$height > lower
  } else {
    chm$chm$height > lower & chm$chm$height < upper
  }
}


#' Determine EPSG Code for a Mosaic
#'
#' This function calculates the EPSG code for a given mosaic based on its
#' geographic extent.
#'
#' @param mosaic A raster object representing the mosaic for which the EPSG code
#'   is to be determined.
#'
#' @return A character string representing the EPSG code corresponding to the
#'   UTM zone and hemisphere of the mosaic's centroid. If the mosaic is not in
#'   the lon/lat coordinate system, a warning is issued.
#'
#' @details The function calculates the centroid of the mosaic's extent,
#'   determines the UTM zone based on the centroid's longitude, and identifies
#'   the hemisphere based on the centroid's latitude. The EPSG code is then
#'   constructed accordingly.
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' library(terra)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#'
#' # Get the EPSG code for the mosaic
#' mosaic_epsg(mosaic)
#' }
#'
#' @export
mosaic_epsg <- function(mosaic) {
  if(terra::is.lonlat(mosaic)){
    extens <- terra::ext(mosaic)
    latitude <- mean(c(terra::ymin(mosaic), terra::ymax(mosaic)))
    longitude <- mean(c(terra::xmin(mosaic), terra::xmax(mosaic)))
    utm_zone <- floor((longitude + 180) / 6) + 1
    hemisphere <- ifelse(latitude >= 0, "N", "S")
    epsg_code <- if (hemisphere == "N") {
      32600 + utm_zone
    } else {
      32700 + utm_zone
    }
    return(paste0("EPSG:", epsg_code))
  } else{
    cli::cli_alert_danger("`mosaic` is not in the lon/lat coordinate system.")
  }
}

#' Project a Mosaic to a New Coordinate Reference System (CRS)
#'
#' This function projects a given mosaic to a specified CRS.
#'
#' @param mosaic A raster object representing the mosaic to be projected.
#' @param y The target CRS to which the mosaic should be projected. This can be
#'   specified in various formats accepted by the [terra::project()] function.
#' @param ... Additional arguments passed to the [terra::project()] function.
#'
#' @return A raster object representing the projected mosaic.
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(terra)
#' library(pliman)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#' mosaic
#' # Define target CRS (EPSG code for WGS 84 / UTM zone 33N)
#' target_crs <- "EPSG:32633"
#'
#' # Project the mosaic
#' projected_mosaic <- mosaic_project(mosaic, "EPSG:32633")
#' projected_mosaic
#' }
#'
#' @export
mosaic_project <- function(mosaic, y, ...){
  return(terra::project(mosaic, y, ...))
}

#' Project a Mosaic from Lon/Lat to EPSG-based CRS
#'
#' This function projects a given mosaic from the lon/lat coordinate system to
#' an EPSG-based CRS determined by the mosaic's extent.
#'
#' @param mosaic A raster object representing the mosaic to be projected. The
#'   mosaic must be in the lon/lat coordinate system.
#'
#' @return A raster object representing the projected mosaic. If the mosaic is
#'   not in the lon/lat coordinate system, a warning is issued.
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(terra)
#' library(pliman)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#'
#' # Project the mosaic to the appropriate UTM zone
#' mosaic_lonlat2epsg(mosaic)
#' }
#'
#' @export
mosaic_lonlat2epsg <- function(mosaic){
  if(terra::is.lonlat(mosaic)){
    epsg <- mosaic_epsg(mosaic)
    return(terra::project(mosaic, epsg))
  } else{
    cli::cli_alert_danger("`mosaic` is not in the lon/lat coordinate system.")
  }
}

#' Extract Values from a Raster Mosaic Using a Shapefile
#'
#' This function extracts values from a raster mosaic based on the regions
#' defined in a shapefile using [exactextractr::exact_extract()].
#'
#' @param mosaic A `SpatRaster` object representing the raster mosaic from which
#'   values will be extracted.
#' @param shapefile A shapefile, which can be a `SpatVector` or an `sf` object,
#'   defining the regions of interest for extraction.
#' @param fun A character string specifying the summary function to be used for
#'   extraction. Default is `"median"`.
#' @param ... Additional arguments to be passed to [exactextractr::exact_extract()].
#' @return A data frame containing the extracted values for each region defined in the shapefile.
#' @export
#'
mosaic_extract <- function(mosaic,
                           shapefile,
                           fun = "median",
                           ...){
  if(inherits(shapefile, "SpatVector")){
    shapefile <- sf::st_as_sf(shapefile)
  }
  results <-
    exactextractr::exact_extract(mosaic,
                                 shapefile,
                                 fun = fun,
                                 force_df = TRUE,
                                 ...)
  sf::st_as_sf(dplyr::bind_cols(results, shapefile)) |> dplyr::relocate(unique_id:column, .before = 1)
}

#' Vectorize a `SpatRaster` mask to an `sf` object
#'
#' Converts a raster mask into a vectorized `sf` object, with various options
#' for morphological operations and filtering.
#'
#' @inheritParams mosaic_analyze
#' @param aggregate The size as a fraction (percentage) of the input image size.
#'   Either a scalar (eg., 50), or a length-two numeric vector. In the last,
#'   different percentage reduction/expansion can be used for columns, and rows,
#'   respectively.
#' @param fill_hull Fill holes in the binary image? Defaults to `FALSE`.
#' @return An `sf` object containing vectorized features from the raster mask,
#'   with added area measurements.
#' @param smooth Smoothes the contours using a moving average filter. Default is
#'   `FALSE`.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' mask <- image_pliman("mask.tif")
#' shp <- mosaic_vectorize(mask, watershed = FALSE)
#' mosaic_plot(mask)
#' shapefile_plot(shp, add = TRUE, lwd = 3)
#' }
#'
#'
mosaic_vectorize <- function(mask,
                             aggregate = NULL,
                             watershed = TRUE,
                             tolerance = 1,
                             extension = 1,
                             opening = FALSE,
                             closing = FALSE,
                             filter = FALSE,
                             erode = FALSE,
                             dilate = FALSE,
                             fill_hull = FALSE,
                             lower_size = NULL,
                             upper_size = NULL,
                             topn_lower = NULL,
                             topn_upper = NULL,
                             smooth = FALSE){
  if(!is.null(aggregate)){
    mask <- mosaic_aggregate(mask, aggregate)
  }
  dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mask), nrow = ncol(mask)))
  extents <- terra::ext(mask)
  dmask[is.na(dmask) == TRUE] <- 0
  if(!isFALSE(fill_hull)){
    dmask <- EBImage::fillHull(dmask)
  }
  if(!isFALSE(filter) & filter > 1){
    dmask <- EBImage::medianFilter(dmask, filter)
  }
  if(is.numeric(erode) & erode > 0){
    dmask <- image_erode(dmask, size = erode)
  }
  if(is.numeric(dilate) & dilate > 0){
    dmask <- image_dilate(dmask, size = dilate)
  }
  if(is.numeric(opening) & opening > 0){
    dmask <- image_opening(dmask, size = opening)
  }
  if(is.numeric(closing) & closing > 0){
    dmask <- image_closing(dmask, size = closing)
  }
  if(watershed){
    dmask <- EBImage::watershed(EBImage::distmap(dmask), tolerance = tolerance, ext = extension)
  } else{
    dmask <- EBImage::bwlabel(dmask)
  }
  conts <- EBImage::ocontour(dmask)
  conts <- conts[sapply(conts, nrow) > 2]
  if(is.numeric(smooth) & smooth > 0){
    conts <- smoothContours(conts, smooth)
  }
  resx <- terra::res(mask)[1]
  resy <- terra::res(mask)[1]
  sf_df <- sf::st_sf(
    geometry = lapply(conts, function(x) {
      tmp <- x + 1
      tmp[, 2] <-  extents[3] + (nrow(mask) - tmp[, 2]) * resy
      tmp[, 1] <- extents[1] + tmp[, 1] * resy
      geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
    }),
    data = data.frame(unique_id = paste0(1:length(conts))),
    crs = terra::crs(mask)
  )
  addmeasures <-
    do.call(rbind,
            lapply(1:nrow(sf_df), function(i){
              compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
            }))
  gridindiv <- cbind(sf_df, addmeasures)

  if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
    cli::cli_abort("{.msg Only one of {.arg lower_*} or {.arg topn_*} can be used.}")

  }
  if(!is.null(lower_size)){
    gridindiv <- gridindiv[gridindiv$area > lower_size, ]
  }
  if(!is.null(upper_size)){
    gridindiv <- gridindiv[gridindiv$area < upper_size, ]
  }
  if(!is.null(topn_lower)){
    gridindiv <- gridindiv[order(gridindiv$area),][1:topn_lower,]
  }
  if(!is.null(topn_upper)){
    gridindiv <- gridindiv[order(gridindiv$area, decreasing = TRUE),][1:topn_upper,]
  }
  return(gridindiv |> check_cols_shp())
}

#' Rotate a mosaic image by specified angles
#'
#' This function rotates a mosaic image by 90, 180, or 270 degrees.
#'
#' @param mosaic A `SpatRaster` object representing the mosaic image.
#' @param angle An integer specifying the rotation angle. Must be one of 90,
#'   180, or 270.
#' @param direction A string specifying the rotation direction. Must be either
#'   "clockwise" or "anticlockwise".
#' @return A `SpatRaster` object with the rotated mosaic image.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' # Convert a mosaic raster to an Image object
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' r90 <- mosaic_rotate(mosaic, 90)
#' r180 <- mosaic_rotate(mosaic, 180)
#' r270 <- mosaic_rotate(mosaic, 270)
#' # Plot all rotations side by side
#' par(mfrow = c(2, 2))
#' mosaic_plot(mosaic, main = "Original")
#' mosaic_plot(r90, main = "90 Degrees")
#' mosaic_plot(r180, main = "180 Degrees")
#' mosaic_plot(r270, main = "270 Degrees")
#' par(mfrow = c(1, 1))
#'
#' }
mosaic_rotate <- function(mosaic, angle, direction = "clockwise") {
  # Ensure angle is one of the allowed values
  if (!angle %in% c(90, 180, 270)) {
    cli::cli_abort("{.arg angle} must be one of {.val 90}, {.val 180}, or {.val 270}.")
  }

  # Ensure direction is valid
  if (!direction %in% c("clockwise", "anticlockwise")) {
    cli::cli_abort("{.arg direction} must be either {.val clockwise} or {.val anticlockwise}.")
  }

  # Apply the appropriate rotation based on direction
  rotated_mosaic <- if (direction == "clockwise") {
    switch(
      as.character(angle),
      "90"  = terra::t(terra::flip(mosaic)),
      "180" = terra::flip(terra::flip(mosaic), "horizontal"),
      "270" = terra::t(terra::flip(mosaic, "horizontal"))
    )
  } else {
    switch(
      as.character(angle),
      "90"  = terra::flip(terra::t(mosaic)),
      "180" = terra::flip(mosaic, "vertical"),
      "270" = terra::t(terra::flip(mosaic))
    )
  }

  # Return the rotated mosaic for further use
  return(rotated_mosaic)
}

#' Classify a Mosaic Based on Index Breaks
#'
#' This function classifies a given raster mosaic based on user-defined breaks.
#' It provides an option to calculate the frequency and area of each class, as
#' well as plot the classified mosaic.
#'
#' @param mosaic A `SpatRaster` object representing the mosaic to be classified.
#' @param breaks A numeric vector specifying the breakpoints for classification.
#' @param frequency Logical. If `TRUE`, computes the class frequency and area (in hectares).
#' @param plot Logical. If `TRUE`, plots the classified mosaic.
#'
#' @return A list with two elements:
#'   - `classified`: A `SpatRaster` object containing the classified mosaic.
#'   - `class_freq`: A data frame containing class frequencies, areas (ha), and percentages (if `frequency = TRUE`).
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' library(terra)
#'
#' # Create an example raster
#' r <- terra::rast(matrix(runif(100, min = 0, max = 1), nrow=10, ncol=10))
#'
#' # Classify the raster
#' result <- mosaic_classify(r, breaks = c(0.3, 0.6))
#'
#' # View results
#' result$classified
#' result$class_freq
#' }
#'
mosaic_classify <- function(mosaic, breaks, frequency = TRUE, plot = TRUE) {
  breaks <- c(-Inf, breaks, Inf)
  classified <- terra::classify(mosaic, breaks, brackets = TRUE, include.lowest = FALSE)

  if (frequency) {
    class_freq <- terra::freq(classified)
    if (terra::crs(mosaic) == "" | terra::is.lonlat(mosaic)) {
      cli::cli_alert_danger("Unknown CRS or lat/lon used. The area will be calculated in square units of raster")
      class_freq$area_ha <- class_freq$count * prod(terra::res(mosaic))
    } else {
      class_freq$area_ha <- (class_freq$count * prod(terra::res(mosaic))) / 10000
    }
    total_pixels <- sum(class_freq$count)
    class_freq$percentage <- (class_freq$count / total_pixels) * 100
  } else {
    class_freq <- NULL
  }

  if (plot) {
    terra::plot(
      classified,
      col = pliman::custom_palette(c("darkred", "yellow", "darkgreen"), n = length(breaks)),
      maxcell = 1e6
    )
  }
  return(list(classified = classified, class_freq = class_freq))
}


#' Clip a Raster Mosaic by Polygons
#'
#' @title Clip Raster Mosaic by Polygons
#'
#' @description
#' Quickly partition a large raster mosaic into individual tiles using a
#' polygon layer.  Each tile is clipped by either the polygon's bounding box
#' or (optionally) the exact feature geometry, and written to disk as a separate
#' GeoTIFF named by the feature's `unique_id`.
#'
#' @details
#' This function wraps GDAL's `warp` utility for efficient raster clipping.
#' When `parallel = TRUE`, it will spawn multiple workers via `mirai` and
#' process tiles in batches.  Use `exact = TRUE` to clip to the true polygon
#' shape (at some extra cost), or leave `exact = FALSE` for a faster
#' bounding-box crop.
#'
#' @param mosaic
#'   A [terra::SpatRaster] object or a file path pointing to a raster.
#'   In-memory rasters are first written to a temporary GeoTIFF.
#' @param shapefile
#'   An [sf::sf], [terra::SpatVector], or path to a vector file.  Must
#'   contain a column named `unique_id` for naming each output tile.
#' @param unique_id A column present in `shapefile` that uniquely identifies the
#'   plots to be clipped.
#' @param out_dir
#'   Directory where clipped rasters will be saved.  Defaults to the current
#'   working directory.  Created recursively if it does not exist.
#' @param overwrite
#'   Logical; if `TRUE` (the default), existing files in `out_dir` with the
#'   same name will be overwritten.
#' @param verbose
#'   Logical; if `TRUE` (default), progress bars and status messages will be
#'   shown.
#' @param exact
#'   Logical; if `FALSE` (default), tiles are cropped by each feature's
#'   bounding box.  If `TRUE`, the function extracts each polygon as a cutline
#'   for an exact crop (slower, but shape-accurate).
#' @param parallel
#'   Logical; if `TRUE` (default), processing is parallelized using `mirai`.
#'   Set to `FALSE` for purely sequential execution.
#' @param workers
#'   Integer; number of parallel daemons to launch when `parallel = TRUE`.
#'   Defaults to 70% of available cores.
#' @export
#' @return
#'   Invisibly returns a character vector of file paths to all clipped GeoTIFFs.
mosaic_clip <- function(mosaic,
                        shapefile,
                        unique_id =  "unique_id",
                        out_dir = NULL,
                        overwrite = TRUE,
                        verbose = TRUE,
                        exact = FALSE,
                        parallel = FALSE,
                        workers = NULL) {
  odir <- if (is.null(out_dir)) "./" else {
    if (grepl("[/\\\\]", out_dir)) out_dir else file.path(".", out_dir)
  }
  if (!dir.exists(odir)) {
    cli::cli_alert_info("Creating output directory {.path {odir}}")
    dir.create(odir, recursive = TRUE)
  }
  shp <- pliman::shapefile_input(shapefile, info = FALSE)
  if (!unique_id %in% names(shp)) {
    cli::cli_abort(
      c(
        "x" = "{.val {unique_id}} column not found in {.arg shapefile}.",
        "i" = "Available columns are: {.val {names(shapefile)}}."
      )
    )

  }
  if (terra::inMemory(mosaic)) {
    tf <- paste0(tempfile(), ".tif")
    on.exit(file.remove(tf), add = TRUE)
    terra::writeRaster(mosaic, filename = tf, overwrite = TRUE)
    tempf <- tf
  } else {
    tempf <- terra::sources(mosaic)
  }

  n   <- nrow(shp)
  ids <- as.character(shp$unique_id)
  # check if ids are unique
  if (length(unique(ids)) != length(ids)) {
    cli::cli_abort("The {.val {unique_id}} column in {.arg shapefile} must contain unique values.")
  }
  out_files <- character(n)
  out <- file.path(odir, paste0(ids, ".tif"))
  if(any(basename(out) %in% list.files(path = odir)) & isFALSE(overwrite)){
    cli::cli_abort("Some raster tiles are already in the destination directory. Use {.arg overwrite = TRUE} to overwrite those files.")
  }

  clip_fun <- function(i, tempf, shp, ids, odir) {
    bb   <- sf::st_bbox(shp[i, ])
    if(overwrite){
      opts <- c(
        "-overwrite",
        "-te",
        bb[["xmin"]], bb[["ymin"]],
        bb[["xmax"]], bb[["ymax"]]
      )
    } else{
      opts <- c(
        "-te",
        bb[["xmin"]], bb[["ymin"]],
        bb[["xmax"]], bb[["ymax"]]
      )
    }
    if (exact) {
      # ensure shapefile is a file path
      cutline_path <- if (is.character(shapefile) && file.exists(shapefile)) {
        shapefile
      } else {
        temp_path <- tempfile(fileext = ".geojson")
        sf::st_write(shapefile, temp_path, driver = "GeoJSON", quiet = TRUE)
        on.exit(unlink(temp_path), add = TRUE)
      }
      opts <- c(
        "-cutline", cutline_path,
        "-cwhere", sprintf("unique_id = '%s'", ids[i]),
        "-crop_to_cutline",
        "-dstalpha",
        if(overwrite){"-overwrite"}
      )
    }

    out <- file.path(odir, paste0(ids[i], ".tif"))
    suppressWarnings(
      sf::gdal_utils(
        util        = "warp",
        source      = tempf,
        destination = out,
        options     = opts
      )
    )
    invisible(out)
  }

  if(parallel){
    nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.7), workers)
    mirai::daemons(nworkers)
    on.exit(mirai::daemons(0))

    cli::cli_rule(
      left  = cli::col_blue("Clipping mosaic into tiles using  {.val {n}} polygons"),
      right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
    )
    cli::cli_progress_step(
      msg        = "Configuring the parallel backend using {.val {nworkers}} clusters",
      msg_done   = "Parallel process finished",
      msg_failed = "Raster clipping failed."
    )


    # break features into 3 chunks
    chunks <- split(seq_len(n), ceiling(seq_len(n)/ceiling(n/nworkers)))

    # return(chunks)
    out_files <- mirai::mirai_map(
      chunks,
      function(idxs, tempf, shp, ids, odir, overwrite, exact) {
        out_paths <- character(length(idxs))

        for (k in seq_along(idxs)) {
          i <- idxs[k]

          # decide between bounding-box vs cutline
          if (exact) {
            cutline_path <- tempfile(fileext = ".geojson")
            sf::st_write(shp[i, ], cutline_path,
                         driver = "GeoJSON", quiet = TRUE)

            opts <- c(
              "-cutline", cutline_path,
              "-cwhere", sprintf("unique_id = '%s'", ids[i]),
              "-crop_to_cutline",
              "-dstalpha",
              if (overwrite) "-overwrite"
            )
          } else {
            bb <- sf::st_bbox(shp[i, ])
            opts <- c(
              "-te",
              bb[["xmin"]], bb[["ymin"]],
              bb[["xmax"]], bb[["ymax"]],
              if (overwrite) "-overwrite"
            )
          }
          out_paths[k] <- file.path(odir, paste0(ids[i], ".tif"))
          sf::gdal_utils(
            util        = "warp",
            source      = tempf,
            destination = out_paths[k],
            options     = opts
          )
        }
        out_paths
      },
      .args = list(
        tempf     = tempf,
        shp       = shp,
        ids       = ids,
        odir      = odir,
        overwrite = overwrite,
        exact     = exact
      )
    )[.progress]
    out_files <- unlist(out_files)
  } else{
    if(verbose){
      cli::cli_rule(
        left  = cli::col_blue("Clipping mosaic into tiles using  {.val {n}} polygons"),
        right = cli::col_blue("Started at {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
      )
      cli::cli_progress_bar(
        format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        total = n,
        clear = TRUE
      )
    }

    for (i in seq_len(n)) {
      if(verbose){
        cli::cli_progress_update()
      }
      out_files[i] <- clip_fun(i, tempf, shp, ids, odir)
    }
  }
  if(verbose){
    cli::cli_progress_done()
    cli::cli_rule(
      left  = cli::col_blue("Raster clipping finished"),
      right = cli::col_blue("Finished on {.val {format(Sys.time(), '%Y-%m-%d | %H:%M:%OS0')}}")
    )
  }
  invisible(out_files)
}

