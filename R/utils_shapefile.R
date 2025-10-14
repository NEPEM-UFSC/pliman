point_to_polygon <- function(sf_object, n_sides = 500) {
  # Extract CRS of the input sf object
  crsobj <- sf::st_crs(sf_object)
  # Create a new geometry list
  new_geometries <- lapply(seq_len(nrow(sf_object)), function(i) {
    geom_type <- sf::st_geometry_type(sf_object[i, ])
    if (geom_type == "POINT") {
      # Get the point coordinates
      point <- sf::st_coordinates(sf_object[i, ])
      radius <- sf_object[["radius"]][i]

      if (is.na(radius)) {
        cli::cli_abort("Radius is missing for a POINT geometry!")
      }
      angles <- seq(0, 2 * pi, length.out = n_sides + 1)
      circle_coords <- cbind(
        point[1] + radius * cos(angles),  # X coordinates
        point[2] + radius * sin(angles)   # Y coordinates
      )
      sf::st_polygon(list(circle_coords))
    } else {
      sf::st_geometry(sf_object[i, ])
    }
  })
  # Function to ensure all geometries in a list are valid sfg objects
  validate_geometries <- function(geometry_list) {
    lapply(geometry_list, function(geom) {
      if (inherits(geom, "sfg")) {
        return(geom)  # Valid sfg object, return as is
      } else if (inherits(geom, "sfc")) {
        return(geom[[1]])  # Unnest if it's an sfc object
      } else if (is.list(geom) && inherits(geom[[1]], "sfg")) {
        return(geom[[1]])  # Handle nested lists containing sfg objects
      } else {
        cli::cli_abort("Invalid geometry found in the list")
      }
    })
  }
  sf_object <-
    sf::st_set_geometry(sf_object, sf::st_sfc(validate_geometries(new_geometries))) |>
    sf::st_set_crs(crsobj)
  return(sf_object)
}
add_width_height <- function(grid, width, height, mosaic, points_align) {
  gridl <-lapply(sf::st_geometry(grid), sf::st_coordinates)
  gridadj <- add_width_height_cpp(gridl, height, width, points_align)
  return(sf::st_sf(geometry = sf::st_as_sfc(gridadj, crs = sf::st_crs(grid))))
}
create_buffer <- function(coords, buffer_col, buffer_row) {
  # Calculate the new x-min, x-max, y-min, and y-max after adjustment
  coords <- sf::st_coordinates(coords)
  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  y_min <- min(coords[, 2])
  y_max <- max(coords[, 2])
  new_x_min <- x_min - buffer_col * (x_max - x_min)
  new_x_max <- x_max + buffer_col * (x_max - x_min)
  new_y_min <- y_min - buffer_row * (y_max - y_min)
  new_y_max <- y_max + buffer_row * (y_max - y_min)

  # Calculate the scaling factors for x and y
  x_scale_factor <- (new_x_max - new_x_min) / (x_max - x_min)
  y_scale_factor <- (new_y_max - new_y_min) / (y_max - y_min)

  # Apply the scaling to the coordinates
  resized_coords <- coords
  resized_coords[, 1] <- (resized_coords[, 1] - x_min) * x_scale_factor + new_x_min
  resized_coords[, 2] <- (resized_coords[, 2] - y_min) * y_scale_factor + new_y_min
  sf::st_polygon(list(resized_coords[, 1:2]))
}

make_grid <- function(points, nrow, ncol, mosaic, buffer_col = 0, buffer_row = 0, plot_width = NULL, plot_height = NULL) {
  points_align <-
    sf::st_transform(points, sf::st_crs(mosaic)) |>
    sf::st_coordinates()

  grids <-
    sf::st_make_grid(points, n = c(nrow, ncol)) |>
    sf::st_transform(sf::st_crs(mosaic))

  sxy <-
    points |>
    sf::st_make_grid(n = c(1, 1)) |>
    sf::st_cast("POINT") |>
    rev() |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_coordinates()
  txy <-
    points |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_coordinates()
  txy <- txy[1:4, 1:2]

  cvm <- lm(txy ~ sxy[1:4, ])
  parms <- cvm$coefficients[2:3, ]
  intercept <- cvm$coefficients[1, ]
  geometry <- sf::st_sf(geometry = grids * parms + intercept, crs = sf::st_crs(mosaic))

  if (buffer_row != 0 | buffer_col != 0) {
    geometry <- lapply(geometry, function(g) {
      create_buffer(g, buffer_col, buffer_row)
    }) |>
      sf::st_sfc(crs = sf::st_crs(mosaic))
  }
  if (!is.null(plot_width) & !is.null(plot_height)) {
    geometry <- add_width_height(grid = geometry, width = plot_width, height = plot_height, points_align = points_align[2:3, 1:2])
  }
  return(geometry)
}

#' Generate plot IDs with different layouts
#'
#' Based on a shapefile, number of columns and rows, generate plot IDs with
#' different layouts.
#'
#' @param shapefile An object computed with [shapefile_build()]
#' @param nrow The number of columns
#' @param ncol The number of rows
#' @param layout Character: one of
#'  * `'tblr'` for top/bottom left/right orientation
#'  * `'tbrl'` for top/bottom right/left orientation
#'  * `'btlr'` for bottom/top left/right orientation
#'  * `'btrl'` for bottom/top right/left orientation
#'  * `'lrtb'` for left/right top/bottom orientation
#'  * `'lrbt'` for left/right bottom/top orientation
#'  * `'rltb'` for right/left top/bottom orientation
#'  * `'rlbt'` for right/left bottom/top orientation
#' @param plot_prefix The plot_id prefix. Defaults to `'P'`.
#' @param serpentine Create a serpentine-based layout? Defaults to `FALSE`.
#' @return A list of plot IDs with specified layout and updated rows/columns.
#' @export
#'
plot_id <- function(shapefile,
                    nrow,
                    ncol,
                    layout = c("tblr", "tbrl", "btlr", "btrl", "lrtb", "lrbt", "rltb", "rlbt"),
                    plot_prefix = "P",
                    serpentine = FALSE) {
  # Ensure the specified layout is valid
  allowed <- c("tblr", "tbrl", "btlr", "btrl", "lrtb", "lrbt", "rltb", "rlbt")
  layout <- layout[[1]]
  if (!layout %in% allowed) {
    cli::cli_abort(c(
      "!" = "{.arg layout} must be one of the following:",
      "v" = "{.val {paste(allowed, collapse = ', ')}}"
    ))
  }

  # Ensure that the number of rows in the shapefile matches expected dimensions
  expected_rows <- nrow * ncol
  if (nrow(shapefile) != expected_rows) {
    cli::cli_abort("Expected {.val {expected_rows}} rows, but {.arg shapefile} has {.val {nrow(shapefile)}} rows.")
  }


  # Helper function for generating plot names
  leading_zeros <- function(x, n) {
    sprintf(paste0("%0", n, "d"), x)
  }

  plots_tblr <- paste0(plot_prefix, leading_zeros(1:nrow(shapefile), 4))
  rows <- shapefile$row
  cols <- shapefile$column

  # Define layout functions
  make_tblr <- function() {
    return(list(plots = plots_tblr, row = rows, col = cols))
  }

  make_tbrl <- function() {
    plots_tblr_rev <- rev(plots_tblr)
    plots_tbrl <- NULL
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_tbrl <- c(plots_tbrl, rev(plots_tblr_rev[start:end]))
    }
    return(list(plots = plots_tbrl, row = rows, col = rev(cols)))


  }

  make_btrl<- function() {
    plots_rev <- rev(plots_tblr)
    plots_btlr <- NULL
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_btlr <- c(plots_btlr, plots_rev[start:end])
    }
    return(list(plots = plots_btlr, row = rev(rows), col = rev(cols)))

  }

  make_btlr <- function() {
    plots_btlr_rev <- rev(make_btrl()$plots)
    plots_btrl <- NULL
    for (i in seq_len(ncol)) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_btrl <- c(plots_btrl, rev(plots_btlr_rev[start:end]))
    }
    return(list(plots = plots_btrl, row = rev(rows), col = cols))

  }

  make_lrtb <- function() {
    plots_lrtb <- NULL
    for (i in 1:ncol) {
      plots_lrtb <- c(plots_lrtb, plots_tblr[seq(i, length(plots_tblr), by = ncol)])
    }
    return(list(plots = plots_lrtb, row = rows, col = cols))
  }

  make_lrbt <- function() {
    plots_lrbt <- NULL
    plots_lrtb <- make_lrtb()$plots
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_lrbt <- c(plots_lrbt, rev(plots_lrtb[start:end]))
    }
    return(list(plots = plots_lrbt, row = rev(rows), col = cols))
  }

  make_rltb <- function() {
    plots_rltb <- NULL
    for (i in 1:ncol) {
      # Columns from right to left
      plots_rltb <- c(plots_rltb, plots_tblr[seq(ncol - i + 1, length(plots_tblr), by = ncol)])
    }
    return(list(plots = plots_rltb, row = rows, col = rev(cols)))
  }

  make_rlbt <- function() {
    plots_rltb <- make_rltb()$plots
    plots_rlbt <- NULL
    for (i in seq_len(ncol)) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_rlbt <- c(plots_rlbt, rev(plots_rltb[start:end]))
    }
    return(list(plots = plots_rlbt, row = rev(rows), col = rev(cols)))
  }

  # Return the appropriate layout
  plots <-  switch(layout,
                   "tblr" = make_tblr()$plots,
                   "tbrl" = make_tbrl()$plots,
                   "btlr" = make_btlr()$plots,
                   "btrl" = make_btrl()$plots,
                   "lrtb" = make_lrtb()$plots,
                   "lrbt" = make_lrbt()$plots,
                   "rltb" = make_rltb()$plots,
                   "rlbt" = make_rlbt()$plots)
  newrows <-  switch(layout,
                     "tblr" = make_tblr()$row,
                     "tbrl" = make_tbrl()$row,
                     "btlr" = make_btlr()$row,
                     "btrl" = make_btrl()$row,
                     "lrtb" = make_lrtb()$row,
                     "lrbt" = make_lrbt()$row,
                     "rltb" = make_rltb()$row,
                     "rlbt" = make_rlbt()$row)
  newcols <-  switch(layout,
                     "tblr" = make_tblr()$col,
                     "tbrl" = make_tbrl()$col,
                     "btlr" = make_btlr()$col,
                     "btrl" = make_btrl()$col,
                     "lrtb" = make_lrtb()$col,
                     "lrbt" = make_lrbt()$col,
                     "rltb" = make_rltb()$col,
                     "rlbt" = make_rlbt()$col)
  mat <- matrix(plots, ncol = ncol, nrow = nrow)
  if(serpentine){
    # column serpentine
    if(layout %in% c("tblr", "btlr")){
      mat2 <- mat
      for (j in 1:ncol(mat)) {
        if(j %% 2 == 0){
          mat2[, j] <- rev(mat[, j])
        } else{
          mat2[, j]
        }
      }
    }
    if(layout %in% c( "tbrl", "btrl")){
      mat2 <- mat
      cols <- ncol(mat):1
      for (j in cols[seq(2, ncol, by = 2)]) {
        mat2[, j] <- rev(mat[, j])
      }
    }
    # row serpentine
    if(layout %in% c("lrtb", "rltb")){
      mat2 <- mat
      for (j in 1:nrow(mat)) {
        if(j %% 2 == 0){
          mat2[j, ] <- rev(mat[j, ])
        } else{
          mat2[j, ]
        }
      }
    }
    if(layout %in% c("lrbt", "rlbt")){
      mat2 <- mat
      rows <- nrow(mat):1
      for (j in rows[seq(2, nrow, by = 2)]) {
        mat2[j, ] <- rev(mat[j, ])
      }
    }
  } else{
    mat2 <- mat
  }
  return(list(plots = as.vector(mat2), rows = newrows, cols = newcols))
}

#' Build a shapefile from a mosaic raster
#'
#' This function takes a mosaic raster to create a shapefile containing polygons
#' for the specified regions. Users can drawn Areas of Interest (AOIs) that can
#' be either a polygon with n sides, or a grid, defined by `nrow`, and `ncol`
#' arguments.
#' @details
#' Since multiple blocks can be created, the length of arguments `grid`, `nrow`,
#' `ncol`, `buffer_edge`, `buffer_col`, and `buffer_row` can be either an scalar
#' (the same argument applied to all the drawn blocks), or a vector with the
#' same length as the number of drawn blocks. In the last, shapefiles in each
#' block can be created with different dimensions.
#' @param sf_to_polygon Convert sf geometry like POINTS and LINES to POLYGONS?
#'   Defaults to `FALSE`. Using `TRUE` allows using POINTS to extract values
#'   from a raster using `exactextractr::exact_extract()`.
#' @param mosaic A `SpatRaster` object, typically imported using
#'   [mosaic_input()].  If not provided, a latitude/longitude basemap will be
#'   generated in the "EPSG:4326" coordinate reference system.
#' @param basemap An optional `mapview` object.
#' @param controlpoints An `sf` object created with [mapedit::editMap()],
#'   containing the polygon that defines the region of interest to be analyzed.
#' @param nsides The number of sides if the geometry is generated with `Draw
#'   Circle` tool.
#' @inheritParams mosaic_analyze
#' @inheritParams mosaic_index
#' @inheritParams mosaic_view
#' @inheritParams utils_shapefile
#' @inheritParams plot_id
#' @return A list with the built shapefile. Each element is an `sf` object with
#'   the coordinates of the drawn polygons.
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' shps <-
#'       shapefile_build(mosaic,
#'                       nrow = 6,
#'                       ncol = 3,
#'                       buffer_row = -0.05,
#'                       buffer_col = -0.25,
#'                       check_shapefile = FALSE,
#'                       build_shapefile = FALSE) ## Use TRUE to interactively build the plots
#' mosaic_plot(mosaic)
#' shapefile_plot(shps[[1]], add = TRUE)
#' }
#'

shapefile_build <- function(mosaic,
                            basemap = NULL,
                            controlpoints = NULL,
                            r = 3,
                            g = 2,
                            b = 1,
                            crop_to_shape_ext = FALSE,
                            grid = TRUE,
                            nrow = 1,
                            ncol = 1,
                            nsides = 200,
                            plot_width = NULL,
                            plot_height = NULL,
                            layout = "lrtb",
                            serpentine = TRUE,
                            build_shapefile = TRUE,
                            check_shapefile = FALSE,
                            sf_to_polygon = FALSE,
                            buffer_edge = 1,
                            buffer_col = 0,
                            buffer_row = 0,
                            as_sf = TRUE,
                            verbose = TRUE,
                            max_pixels = 1000000,
                            downsample = NULL,
                            quantiles =  c(0, 1)){
  nomosaic <- missing(mosaic)
  check_mapview()
  if(nomosaic){
    mosaic <- terra::rast(nrows=180, ncols=360, nlyrs=3, crs = "EPSG:4326")
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:4326")
  }
  ress <- terra::res(mosaic)
  if(build_shapefile){
    if(verbose){
      cli::cli_progress_step("Building the mosaic",
                             msg_done = "Mosaic built",
                             msg_failed = "Failed to build mosaic")
    }
    if(is.null(basemap)){
      basemap <- mosaic_view(mosaic,
                             r = r,
                             g = g,
                             b = b,
                             max_pixels = max_pixels,
                             downsample = downsample,
                             quantiles = quantiles)
      if(nomosaic){
        basemap@map <- leaflet::setView(basemap@map, lng = -51, lat = -14, zoom = 4)
      }
    }
    if(is.null(controlpoints)){
      points <- mapedit::editMap(basemap,
                                 editor = "leafpm",
                                 editorOptions = list(toolbarOptions = list(
                                   drawMarker = TRUE,
                                   drawPolygon = TRUE,
                                   drawPolyline = TRUE,
                                   drawCircle = TRUE,
                                   drawRectangle = TRUE,
                                   editMode = TRUE,
                                   cutPolygon = TRUE,
                                   removalMode = TRUE,
                                   position = "topleft"
                                 ))
      )
      cpoints <-
        points$finished |>
        sf::st_transform(sf::st_crs(mosaic)) |>
        point_to_polygon(n_sides = nsides)
    } else{
      cpoints <- controlpoints
    }
    if(sf_to_polygon){
      cpoints <- cpoints |> sf_to_polygon()
    }
  } else{
    extm <- terra::ext(mosaic)
    xmin <- extm[1]
    xmax <- extm[2]
    ymin <- extm[3]
    ymax <- extm[4]
    coords <- matrix(c(xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin, xmin, ymax), ncol = 2, byrow = TRUE)
    # Create a Polygon object
    polygon <- sf::st_polygon(list(coords))
    # Create an sf object with a data frame that includes the 'geometry' column
    if(sum(ress) == 2){
      crop_to_shape_ext <- FALSE

    }
    cpoints <- sf::st_sf(data.frame(id = 1),
                         geometry = sf::st_sfc(polygon),
                         crs = sf::st_crs(mosaic))
  }

  # crop to the analyzed area
  if(crop_to_shape_ext){
    if(verbose){
      cli::cli_progress_step("Cropping the mosaic",
                             msg_done = "Mosaic cropped",
                             msg_failed = "Failed to crop mosaic")
    }
    if(sum(ress) != 2){
      cpoints <- cpoints |> sf::st_transform(crs = sf::st_crs(terra::crs(mosaic)))
    }
    poly_ext <-
      cpoints |>
      terra::vect() |>
      terra::buffer(buffer_edge) |>
      terra::ext()
    mosaiccr <- terra::crop(mosaic, poly_ext)

  } else{
    mosaiccr <- mosaic
  }
  # check the parameters
  nrow <- validate_and_replicate2(nrow, cpoints, verbose = verbose)
  ncol <- validate_and_replicate2(ncol, cpoints, verbose = verbose)
  layout <- validate_and_replicate2(layout, cpoints, verbose = verbose)
  buffer_col <- validate_and_replicate2(buffer_col, cpoints, verbose = verbose)
  buffer_row <- validate_and_replicate2(buffer_row, cpoints, verbose = verbose)
  plot_width <- validate_and_replicate2(plot_width, cpoints, verbose = verbose)
  plot_height <- validate_and_replicate2(plot_height, cpoints, verbose = verbose)
  serpentine <- validate_and_replicate2(serpentine, cpoints, verbose = verbose)
  grid <- validate_and_replicate2(grid, cpoints, verbose = verbose)

  # check the created shapes?
  if(verbose){
    cli::cli_progress_step("Creating the shapes",
                           msg_done = "Shapes created",
                           msg_failed = "Failed to create shapes")
  }
  created_shapes <- list()
  for(k in 1:nrow(cpoints)){
    if(inherits(cpoints[k, ]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(cpoints[k, ])) == 5 & grid[[k]]){
      pg <-
        make_grid(cpoints[k, ],
                  nrow = nrow[k],
                  ncol = ncol[k],
                  mosaic = mosaic,
                  buffer_col = buffer_col[k],
                  buffer_row = buffer_row[k],
                  plot_width = plot_width[k],
                  plot_height = plot_height[k]) |>
        dplyr::mutate(row = rep(1:nrow[k], ncol[k]),
                      column = rep(1:ncol[k], each = nrow[k]),
                      .before = geometry)
      updateids <- plot_id(pg, nrow = nrow[k], ncol = ncol[k], layout = layout[k], serpentine = serpentine[k])
      pg <-
        pg |>
        dplyr::mutate(unique_id = dplyr::row_number(),
                      block = paste0("B", leading_zeros(1, 2)),
                      plot_id = updateids$plots,
                      row = updateids$rows,
                      column = updateids$cols,
                      .before = 1)
    } else{
      pg <-
        cpoints[k, ] |>
        sf::st_transform(sf::st_crs(mosaic)) |>
        dplyr::select(geometry) |>
        dplyr::mutate(unique_id = dplyr::row_number(),
                      block = paste0("B", leading_zeros(k, 2)),
                      plot_id = "P0001",
                      row = 1,
                      column = 1,
                      .before = 1)
    }

    created_shapes[[k]] <- pg

  }
  if(check_shapefile){
    if(verbose){
      cli::cli_progress_step("Checking the built shapefile",
                             msg_done = "Shapefile checked",
                             msg_failed = "Failed to check shapefile")
    }
    lengths <- sapply(created_shapes, nrow)
    pg_edit <-
      do.call(rbind, lapply(seq_along(created_shapes), function(i){
        created_shapes[[i]] |>
          dplyr::mutate(`_leaflet_id` = 1:nrow(created_shapes[[i]]),
                        feature_type = "polygon") |>
          dplyr::relocate(geometry, .after = 4) |>
          sf::st_transform(crs = 4326)
      }))
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- mosaic_aggregate(mosaiccr, pct = round(100 / downsample))
    }
    if(build_shapefile){
      mapview::mapview() |> mapedit::editMap()
    }
    edited <-
      mapedit::editFeatures(pg_edit, basemap) |>
      dplyr::select(geometry, block, plot_id, row, column) |>
      dplyr::mutate(unique_id = dplyr::row_number(), .before = 1) |>
      sf::st_transform(sf::st_crs(mosaic))
    sfeat <- sf::st_as_sf(edited)
    sf::st_geometry(sfeat) <- "geometry"
    created_shapes <- split(sfeat, edited$block)
  }
  if(verbose){
    cli::cli_progress_step("Finishing the shapefile",
                           msg_done = "Shapefile built",
                           msg_failed = "Failed to export shapefile")
  }
  if(!as_sf){
    return(shapefile_input(created_shapes, info = FALSE, as_sf = FALSE))
  } else{
    return(created_shapes)
  }
}


#' A wrapper around terra::plot()
#'
#' Plot the values of a SpatVector
#'
#' @param shapefile An SpatVector of sf object.
#' @param ... Further arguments passed on to [terra::plot()].
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' r <- shapefile_input(system.file("ex/lux.shp", package="terra"))
#' shapefile_plot(r)
#' }
shapefile_plot <- function(shapefile, ...){
  if(!inherits(shapefile, "SpatVector") & !inherits(shapefile, "sf") ){
    cli::cli_abort(c(
      "!" = "{.arg shapefile} must be an object of class {.cls SpatVector} or {.cls sf}."
    ))
  }
  if(inherits(shapefile, "sf")){
    shapefile <- terra::vect(shapefile)
  }
  terra::plot(shapefile, ...)
}


#' Import/export shapefiles.
#' @description
#'
#' * `shapefile_input()` creates or imports a shapefile and optionally converts
#'  it to an `sf` object. It can also cast `POLYGON` or `MULTIPOLYGON` geometries
#'  to `MULTILINESTRING` if required.
#' * `shapefile_export()` exports an object (`sf` or `SpatVector`) to a file.
#' * `shapefile_view()` is a simple wrapper around `mapview()` to plot a shapefile.
#'
#' @name utils_shapefile
#'
#' @param shapefile For `shapefile_input()`, character (filename), or an object that can be
#'   coerced to a SpatVector, such as an `sf` (simple features) object. See
#'   [terra::vect()] for more details.
#'
#'   For `shapefile_export()`, a `SpatVector` or an `sf` object to be exported as
#'   a shapefile.
#'
#' @param info Logical value indicating whether to print information about the
#'   imported shapefile (default is `TRUE`).
#' @param as_sf Logical value indicating whether to convert the imported
#'   shapefile to an `sf` object (default is `TRUE`).
#' @param multilinestring Logical value indicating whether to cast polygon geometries
#'   to `MULTILINESTRING` geometries (default is `FALSE`).
#' @param type A character string specifying whether to visualize the shapefile
#'   as `"shape"` or as `"centroid"`. Partial matching is allowed. If set to
#'   `"centroid"`, the function will convert the shapefile's geometry to
#'   centroids before displaying. Defaults to `"shape"`.
#' @param filename The path to the output shapefile.
#' @param attribute The attribute to be shown in the color key. It must be a
#'   variable present in `shapefile`.
#' @param color_regions The color palette to represent `attribute`.
#' @param ... Additional arguments to be passed to [terra::vect()]
#'   (`shapefile_input()`), [terra::writeVector()] (`shapefile_export()`) or
#'   [mapview::mapview()] (`shapefile_view()`).
#'
#' @return
#'  * `shapefile_input()` returns an object of class `sf` (default) representing
#'  the imported shapefile.
#'
#'  * `shapefile_export()` returns a `NULL` object.
#'
#'  * `shapefile_view()` returns an object of class `mapview`.
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' shp <- system.file("ex/lux.shp", package="terra")
#' shp_file <- shapefile_input(shp, as_sf = FALSE)
#' shapefile_view(shp_file)
#' }
#'
#' @export
shapefile_input <- function(shapefile,
                            info = TRUE,
                            as_sf = TRUE,
                            multilinestring = FALSE,
                            ...) {
  check_mapview()
  # Check if shapefile is a URL and download it
  if (is.character(shapefile) && grepl("^http", shapefile)) {
    check_and_install_package("curl")
    temp_shapefile <- tempfile(fileext = ".rds")
    curl::curl_download(shapefile, temp_shapefile)
    shapefile <- temp_shapefile
  }

  create_shp <- function(shapefile, info, as_sf, ...){
    shp <- terra::vect(shapefile, ...)
    if (terra::crs(shp) == "") {
      cli::cli_abort("Missing Coordinate Reference System. Setting to EPSG:3857")
      terra::crs(shp) <- terra::crs("EPSG:3857")
    }
    if (as_sf) {
      shp <- sf::st_as_sf(shp)
      if(multilinestring){
        shp <- sf::st_cast(shp, "MULTILINESTRING")
      }
    }
    if (info) {
      print(shp)
    }
    return(shp)
  }

  if(inherits(shapefile, "list")){
    shapes <- do.call(rbind, lapply(shapefile, function(x){x}))
    create_shp(shapes, info, as_sf, ...) |> add_missing_columns()
  } else{
    create_shp(shapefile, info, as_sf, ...) |> add_missing_columns()
  }
}
#' @name utils_shapefile
#' @export
shapefile_export <- function(shapefile, filename, ...) {
  if (inherits(shapefile, "list")) {
    shapefile <- shapefile_input(shapefile, info = FALSE)
  }
  if (!inherits(shapefile, "SpatVector")) {
    shapefile <- try(terra::vect(shapefile))
  }
  if (inherits(shapefile, "try-error")) {
    cli::cli_abort(c(
      "!" = "{.arg shapefile} must be an object of class {.cls SpatVector} or {.cls sf}."
    ))
  }
  terra::writeVector(shapefile, filename, ...)
}

#' @name utils_shapefile
#' @export
shapefile_view <- function(shapefile,
                           attribute = NULL,
                           type = c("shape", "centroid"),
                           color_regions = custom_palette(c("red", "yellow", "forestgreen")),
                           ...){
  type <- match.arg(type, choices = c("shape", "centroid"))
  # Example usage of the checked 'type'
  if (type == "centroid") {
    shapefile <- suppressWarnings(sf::st_centroid(shapefile))
  }
  if(!is.null(attribute) && attribute == "plot_id"){
    if(inherits(shapefile, "list")){
      shapefile <-
        lapply(shapefile, function(x){
          x |> dplyr::mutate(plot_id = as.numeric(gsub("\\D", "", plot_id)))
        })
    } else{
      shapefile <- shapefile |> dplyr::mutate(plot_id = as.numeric(gsub("\\D", "", plot_id)))
    }
  }
  suppressWarnings(
    mapview::mapview(shapefile,
                     zcol = attribute,
                     col.regions = color_regions,
                     ...)
  )
}

#' Edit Features in a Shapefile
#'
#' This function allows you to interactively edit features in a shapefile using
#' the mapedit package.
#'
#' @param shapefile A shapefile (`sf` object) that can be created with
#'   [shapefile_input()].
#' @param mosaic Optionally, a mosaic (SpatRaster) to be displayed as a
#'   background.
#' @param basemap An optional `mapview` object.
#' @param r Red band index for RGB display (default is 3).
#' @param g Green band index for RGB display (default is 2).
#' @param b Blue band index for RGB display (default is 1).
#' @param max_pixels Maximum number of pixels for down-sampling the mosaic
#'   (default is 3e6).
#' @export
#' @return A modified shapefile with user-edited features.
#'
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#' shp <- shapefile_input(system.file("ex/lux.shp", package="terra"))
#' edited <- shapefile_edit(shp)
#' }
shapefile_edit <- function(shapefile,
                           mosaic = NULL,
                           basemap = NULL,
                           r = 3,
                           g = 2,
                           b = 1,
                           max_pixels = 3e6){
  shapefile <- shapefile_input(shapefile, info = FALSE)
  if(!is.null(mosaic)){
    if(is.null(basemap)){
      downsample <- find_aggrfact(mosaic, max_pixels = max_pixels)
      if(downsample > 0){
        mosaic <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
      }
      nlyrs <- terra::nlyr(mosaic)
      if(nlyrs > 2){
        map <-
          mapview::viewRGB(
            x = as(mosaic, "Raster"),
            layer.name = "base",
            r = r,
            g = g,
            b = b,
            na.color = "#00000000",
            maxpixels = 5e6
          )
      } else{
        map <-
          mapview::mapview() %>%
          leafem::addGeoRaster(x = as(mosaic[[1]], "Raster"),
                               colorOptions = leafem::colorOptions(palette = custom_palette(),
                                                                   na.color = "transparent"))
      }
    } else{
      map <- basemap
    }
    edited <- mapedit::editFeatures(shapefile |> sf::st_transform(crs = 4326), map)
  } else{
    edited <- mapedit::editFeatures(shapefile)
  }
  return(edited)
}

#' Extract geometric measures from a shapefile object
#'
#' `shapefile_measures()` calculates key geometric measures such as the number
#' of points, area, perimeter, width, height, and centroid coordinates for a
#' given shapefile (polygon) object.
#'
#' @param shapefile An `sf` object representing the shapefile. It should contain
#'   polygonal geometries for which the measures will be calculated.
#' @param n An integer specifying the number of polygons to process. If `NULL`,
#'   all polygons are considered.
#' @return A modified `sf` object with added columns for:
#' - `xcoord`: The x-coordinate of the centroid.
#' - `ycoord`: The y-coordinate of the centroid.
#' - `area`: The area of the polygon (in square units).
#' - `perimeter`: The perimeter of the polygon (in linear units).
#' - `width`: The calculated width based on sequential distances between points.
#'  The result will only be accurate if the polygon is rectangular.
#' - `height`: The calculated height based on sequential distances between points.
#'  The result will only be accurate if the polygon is rectangular.
#'
#' @details
#' This function processes a single or multi-polygon `sf` object and computes
#' geometric properties. It calculates distances between points, extracts the
#' centroid coordinates, and computes the area and perimeter of the polygons.
#' The width and height are derived from sequential distances between points.
#'
#' @export
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'
#' path_shp <- paste0(image_pliman(), "/soy_shape.rds")
#' shp <- shapefile_input(path_shp)
#' shapefile_measures(shp)
#' }
#'
#'

shapefile_measures <- function(shapefile, n = NULL) {
  if (inherits(shapefile, "list")) {
    shapefile <- shapefile_input(shapefile, info = FALSE)
  }
  if(!is.null(n)){
    shapefile <- shapefile |> dplyr::slice(1:n)
  }
  results <- lapply(1:nrow(shapefile), function(i) {
    # Extract points from the current geometry
    geom <- shapefile[i, ]$geometry
    if(nrow(geom[[1]][[1]]) == 5){
      points <- sf::st_cast(geom, "POINT")
      # Calculate pairwise distances between points
      dists <- suppressWarnings(as.matrix(sf::st_distance(points)))
      # Extract width and height
      width <- as.numeric(round(dists[[4]], 3))
      height <- as.numeric(round(dists[[2]], 3))
      c(width, height)
    } else{
      c(NA, NA)
    }
  })
  wh <- do.call(rbind, results)

  # Calculate the centroid and add measurements
  coords <- suppressWarnings(sf::st_centroid(shapefile))|> sf::st_coordinates()
  measures <-
    shapefile |>
    dplyr::mutate(
      xcoord = coords[, 1],
      ycoord = coords[, 2],
      area = as.numeric(sf::st_area(shapefile)),
      perimeter = rcpp_st_perimeter(as.list(sf::st_geometry(shapefile))),
      width = wh[, 1],
      height = wh[, 2],
      .before = geometry
    )

  return(measures)
}


#' Interpolate values at specific points based on coordinates and a target variable
#'
#' This function interpolates values at specified points using x, y coordinates and a target variable
#' from a shapefile. It supports "Kriging" and "Tps" interpolation methods.
#'
#' @param shapefile An sf object containing the x, y, and target variable (z)
#'   columns. It is highly recommended to use `shapefile_measures()` to obtain
#'   this data.
#' @param z A string specifying the name of the column in the shapefile that
#'   contains the target variable to be interpolated.
#' @param x A string specifying the name of the column containing x-coordinates.
#'   Default is 'x'.
#' @param y A string specifying the name of the column containing y-coordinates.
#'   Default is 'y'.
#' @param interpolation A character vector specifying the interpolation method.
#'   Options are "Kriging" or "Tps".
#' @param verbose Logical; if TRUE, progress messages will be displayed.
#'
#' @return A vector of interpolated values at the specified points.
#' @export
shapefile_interpolate <- function(shapefile,
                                  z,
                                  x = 'x',
                                  y = 'y',
                                  interpolation = c("Kriging", "Tps"),
                                  verbose = FALSE) {

  check_and_install_package("fields")

  # Validate shapefile input
  if (!all(c(x, y, z) %in% names(shapefile))) {
    cli::cli_abort(c(
      "!" = "{.arg shapefile} must contain the columns {.val x}, {.val y}, and {.val z}."
    ))
  }

  # Extract coordinates and values from shapefile object
  xy <- cbind(shapefile[[x]], shapefile[[y]])
  values <- shapefile[[z]]  # The target variable to be interpolated

  if (verbose){
    cli::cli_alert_info("Interpolating the points using ", interpolation[[1]], "...")
  }

  # Perform interpolation
  if (interpolation[[1]] == "Kriging") {
    fit <- suppressMessages(suppressWarnings(fields::Krig(xy, values, aRange = 20)))
  } else if (interpolation[[1]] == "Tps") {
    fit <- suppressMessages(suppressWarnings(fields::Tps(xy, values)))
  } else {
    cli::cli_abort("Invalid interpolation method. Choose {.val Kriging} or {.val Tps}.")
  }
  return(fit)
}

#' Generate a spatial surface plot based on interpolated values
#'
#' This function creates a surface plot from an interpolated spatial model, with options to customize
#' plot appearance, grid resolution, and color palette.
#'
#' @param model An interpolated spatial object (e.g., from `shapefile_interpolate()`) containing the data for plotting.
#' @param curve Logical; if TRUE, a contour plot is generated (`type = "C"`), otherwise an image plot (`type = "I"`). Default is TRUE.
#' @param nx Integer; the number of grid cells in the x-direction. Default is 300.
#' @param ny Integer; the number of grid cells in the y-direction. Default is 300.
#' @param xlab Character; label for the x-axis. Default is "Longitude (UTM)".
#' @param ylab Character; label for the y-axis. Default is "Latitude (UTM)".
#' @param col A color palette function for the surface plot. Default is a custom palette from dark red to yellow to forest green.
#' @param ... Additional parameters to pass to `fields::surface`.
#'
#' @return A surface plot showing spatially interpolated data.
#' @export

shapefile_surface <- function(model,
                              curve = TRUE,
                              nx = 300,
                              ny = 300,
                              xlab = "Longitude (UTM)",
                              ylab = "Latitude (UTM)",
                              col = custom_palette(c("darkred", "yellow", "forestgreen"), n = 100),
                              ...) {

  check_and_install_package("fields")

  # Generate the surface plot
  fields::surface(model,
                  type = ifelse(curve, "C", "I"), # "C" for contour, "I" for image plot
                  nx = nx,
                  ny = ny,
                  xlab = xlab,
                  ylab = ylab,
                  col = col,
                  ...)
}

check_cols_shp <- function(shpimp){
  if(!"unique_id" %in% colnames(shpimp)){
    shpimp <- shpimp |> dplyr::mutate(unique_id = dplyr::row_number())
  }
  if(!"block" %in% colnames(shpimp)){
    shpimp <- shpimp |> dplyr::mutate(block = "B01")
  }
  if(!"plot_id" %in% colnames(shpimp)){
    shpimp <- shpimp |> dplyr::mutate(plot_id = paste0("P", leading_zeros(1:nrow(shpimp), 3)))
  }
  if(!"row" %in% colnames(shpimp)){
    shpimp <- shpimp |> dplyr::mutate(row = 1)
  }
  if(!"column" %in% colnames(shpimp)){
    shpimp <- shpimp |> dplyr::mutate(column = 1)
  }
  shpimp |> dplyr::relocate(geometry, .after = dplyr::last_col())
}

#' Spatial Operations on Shapefiles
#'
#' These functions perform various spatial operations on two shapefiles, including determining which geometries fall within, outside, touch, cross, overlap, or intersect another geometry. They also include functions for geometric operations such as intersection, difference, and union.
#'
#' @param shp1 An `sf` object representing the first shapefile.
#' @param shp2 An `sf` object representing the second shapefile.
#'
#' @details All functions ensure that the coordinate reference systems (CRS) of both shapefiles are the same before performing operations. If the CRSs are different, `shp2` will be transformed to match the CRS of `shp1`.
#' - `shapefile_within()`: Filters features in `shp1` that are fully within `shp2`.
#' - `shapefile_outside()`: Filters features in `shp1` that are outside or do not overlap `shp2`.
#' - `shapefile_overlaps()`: Filters features in `shp1` that overlap with `shp2`.
#' - `shapefile_touches()`: Filters features in `shp1` that touch the boundary of `shp2`.
#' - `shapefile_crosses()`: Filters features in `shp1` that cross through `shp2`.
#' - `shapefile_intersection()`: Computes the geometric intersection of `shp1` and `shp2`.
#' - `shapefile_difference()`: Computes the geometric difference of `shp1` minus `shp2`.
#' - `shapefile_union()`: Computes the geometric union of `shp1` and `shp2`.
#'
#' @return A filtered `sf` object or the result of the geometric operation.
#' @examples
#' if (interactive() && requireNamespace("EBImage")) {
#' library(pliman)
#'
#' shp1 <- shapefile_input(paste0(image_pliman(), "/shp1.rds"))
#' shp2 <- shapefile_input(paste0(image_pliman(), "/shp2.rds"))
#' shapefile_view(shp1) + shapefile_view(shp1)
#'
#' # Apply operations
#' shapefile_within(shp1, shp2)
#' shapefile_outside(shp1, shp2)
#' shapefile_overlaps(shp1, shp2)
#' shapefile_touches(shp1, shp2)
#' shapefile_crosses(shp1, shp2)
#' shapefile_intersection(shp1, shp2)
#' shapefile_difference(shp1, shp2)
#' shapefile_union(shp1, shp2)
#' }
#' @name shapefile_operations
#' @export
#'
shapefile_within <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  bin <- sf::st_within(shp1, shp2) |> as.logical()
  bin[is.na(bin)] <- FALSE
  shp1 |> dplyr::filter(bin)
}
#' @name shapefile_operations
#' @export
shapefile_outside <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  bin <- sf::st_within(shp1, shp2) |> as.logical()
  over <- sf::st_overlaps(shp1, shp2) |> as.logical()
  bin[is.na(bin)] <- FALSE
  over[is.na(over)] <- FALSE
  filt <- as.logical(bin + over)
  shp1 |> dplyr::filter(!filt)
}
#' @name shapefile_operations
#' @export
shapefile_overlaps <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  bin <- sf::st_overlaps(shp1, shp2) |> as.logical()
  bin[is.na(bin)] <- FALSE
  shp1 |> dplyr::filter(bin)
}
#' @name shapefile_operations
#' @export
shapefile_touches <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  bin <- sf::st_touches(shp1, shp2) |> as.logical()
  bin[is.na(bin)] <- FALSE
  shp1 |> dplyr::filter(bin)
}
#' @name shapefile_operations
#' @export
shapefile_crosses <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  bin <- sf::st_crosses(shp1, shp2) |> as.logical()
  bin[is.na(bin)] <- FALSE
  shp1 |> dplyr::filter(bin)
}
#' @name shapefile_operations
#' @export
shapefile_intersection <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  suppressWarnings(sf::st_intersection(shp1, shp2))
}
#' @name shapefile_operations
#' @export
shapefile_difference <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  suppressWarnings(sf::st_difference(shp1, shp2))
}
#' @name shapefile_operations
#' @export
shapefile_union <- function(shp1, shp2){
  if(sf::st_crs(shp1) != sf::st_crs(shp2)){
    cli::cli_alert_warning("CRS are different. Matching CRS of {.arg shp2} to {.arg shp1}.")
    shp2 <- shp2 |> sf::st_transform(sf::st_crs(shp1))
  }
  suppressWarnings(sf::st_union(shp1, shp2))
}

#' Extract mid‐lines from half‐plots
#'
#' For each polygon in an `sf` object, computes the line segment joining
#' the midpoints of the longer pair of opposite edges (the “half‐plot line”).
#'
#' @param shapefile An `sf` object of polygons. Each geometry must be closed
#'   (first and last coordinate coincide) so that `st_coordinates(...)`
#'   yields a repeating start point.
#' @return A `SpatVector` (from the **terra** package) of line geometries
#'   representing the half‐plot midlines.
#' @export
#'
#' @examples
#'
#' if(interactive()){
#' library(pliman)
#' shp <- shapefile_input( paste0(image_pliman(), "/soy_shape.rds"))
#' mosaic <- mosaic_input( paste0(image_pliman(), "/soy_dsm.tif"))
#' mosaic_plot(mosaic)
#' half <- line_on_halfplot(shp)
#' shapefile_plot(half, add = TRUE, col = "blue")
#'
#' }
line_on_halfplot <- function(shapefile) {
  coords <- sf::st_coordinates(shapefile)
  corners_list <- split(coords[, c("X","Y")], coords[, "L2"])
  wkt_lines <- corners_to_wkt(corners_list)
  terra::vect(wkt_lines)
}
