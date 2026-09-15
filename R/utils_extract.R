#' Get Available System Physical Memory (RAM)
#'
#' @description
#' `get_free_ram()` retrieves the current available physical RAM of the system in real-time
#' using native OS system calls in C++ (`GlobalMemoryStatusEx` on Windows, `sysconf` on Linux,
#' and Mach Kernel `host_statistics` on macOS). It executes in microseconds with zero memory allocation.
#'
#' @param unit Character string specifying memory unit: `"GB"` (default), `"MB"`, `"KB"`, or `"B"` (bytes).
#'
#' @return A numeric value representing available system physical memory in the specified unit.
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#' # Get free RAM in GB
#' get_free_ram()
#'
#' # Get free RAM in MB
#' get_free_ram("MB")
#' }
get_free_ram <- function(unit = c("GB", "MB", "KB", "B")) {
  unit <- match.arg(unit)
  ram_bytes <- tryCatch(get_free_ram_cpp(), error = function(e) -1.0)
  if (is.null(ram_bytes) || !is.numeric(ram_bytes) || ram_bytes <= 0) {
    ram_bytes <- 8 * 1024^3
  }

  res <- switch(
    unit,
    "GB" = ram_bytes / (1024^3),
    "MB" = ram_bytes / (1024^2),
    "KB" = ram_bytes / 1024,
    "B"  = ram_bytes
  )
  return(res)
}

#' Extract Values from a Raster Mosaic Using a Shapefile
#'
#' @description
#' `mosaic_extract()` (and its alias `extract_cpp()`) extracts pixel values or summary statistics from a
#' `terra::SpatRaster` object using an `sf` or `SpatVector` shapefile as spatial reference.
#' It runs on a highly optimized C++ engine with OpenMP multi-threading and exact coverage calculations,
#' completely independent of external GDAL libraries or exactextractr.
#' For rasters on disk, it dynamically monitors available system RAM via native C++ API (`get_free_ram()`) and automatically selects between
#' full in-memory processing and memory-safe spatial block processing (when raster size exceeds 80% of available RAM).
#'
#' @param mosaic A `SpatRaster` object (in memory or on disk).
#' @param shapefile An `sf` data frame, `sfc` geometry list, or `SpatVector` object defining regions.
#' @param fun Optional summary function(s) to compute for each polygon: `NULL` (default, returns raw pixel values list),
#'   `"mean"`, `"median"`, `"sum"`, `"min"`, `"max"`, `"sd"`, `"count"`, `"quantiles"` (or `"quantile"`), or a vector of functions such as `c("mean", "quantiles")`.
#' @param return Character string specifying the return object type when statistics are computed:
#'   `"df"` (default, returns summary data frame) or `"shapefile"` (returns original shapefile object with extracted columns bound to it).
#' @param exact Logical. If `TRUE`, computes exact pixel coverage area fractions. Defaults to `FALSE`.
#' @param coverage_area Logical. Include polygon coverage area metrics (`covered_area`, `plot_area`, `coverage`) when computing summary statistics, or include raw `coverage_area` column when `fun = NULL`? Defaults to `FALSE`.
#' @param subdiv Integer subdivision grid size per dimension when `exact = TRUE`. Defaults to `5`.
#' @param summarize_quantiles Numeric vector of quantile probabilities (values between 0 and 1) to compute when `"quantiles"` (or `"quantile"`) is included in `fun`. Defaults to `c(0.05, 0.975)`.
#' @param verbose Logical. Displays progress steps and status messages? Defaults to `FALSE`.
#' @param max_cells_in_memory Optional override for maximum cell count. Defaults to `NULL` (dynamically computed based on 80% of available system RAM).
#' @param ... Additional arguments.
#'
#' @return
#' * If `fun = NULL`, a `list` of `data.frame`s (one for each polygon) containing columns `cell`, `row`, `col`, and raster layer values (plus `coverage_area` if `coverage_area = TRUE`).
#' * If `fun` is specified and `return = "shapefile"`, returns the **original shapefile** (`sf` or `SpatVector`) with the extracted summary columns attached (including `covered_area`, `plot_area`, `coverage` if `coverage_area = TRUE`).
#' * If `fun` is specified and `return = "df"` (default), returns a `data.frame` of summary statistics (including `covered_area`, `plot_area`, `coverage` if `coverage_area = TRUE`).
#'
#' @export
#' @examples
#' \dontrun{
#' library(pliman)
#' shp <- shapefile_input( paste0(image_pliman(), "/soy_shape.rds"))
#' mosaic <- mosaic_input( paste0(image_pliman(), "/soy_ortho.tif"))
#'
#' # Extract summary statistics attached to shapefile
#' mosaic_extract(mosaic, shp, fun = "mean", return = "shapefile")
#'
#' # Extract mean and quantiles with plot area coverage metrics
#' mosaic_extract(mosaic, shp, fun = c("mean", "quantiles"), coverage_area = TRUE)
#' }
mosaic_extract <- function(mosaic,
                           shapefile,
                           fun = NULL,
                           return = "df",
                           exact = FALSE,
                           coverage_area = FALSE,
                           subdiv = 5L,
                           summarize_quantiles = c(0.05, 0.975),
                           verbose = FALSE,
                           max_cells_in_memory = NULL,
                           ...) {
  t_start <- Sys.time()
  orig_shapefile <- shapefile
  return <- match.arg(return, c("shapefile", "df"))

  if (!inherits(mosaic, "SpatRaster")) {
    cli::cli_abort("Argument {.arg mosaic} must be an object of class {.cls SpatRaster}.")
  }

  if (inherits(shapefile, "SpatVector")) {
    sf_poly <- sf::st_as_sf(shapefile)
  } else if (inherits(shapefile, "sf") || inherits(shapefile, "sfc")) {
    sf_poly <- shapefile
  } else {
    cli::cli_abort("Argument {.arg shapefile} must be an {.cls sf} or {.cls SpatVector} object.")
  }

  if (is.function(fun)) {
    if (verbose) {
      cli::cli_progress_step("Extracting pixel values for custom R function...")
    }
    raw_res <- mosaic_extract(
      mosaic = mosaic,
      shapefile = orig_shapefile,
      fun = NULL,
      exact = exact,
      coverage_area = coverage_area,
      subdiv = subdiv,
      verbose = verbose,
      max_cells_in_memory = max_cells_in_memory
    )

    if (verbose) {
      cli::cli_progress_step("Applying custom R summary function across features...")
    }

    first_arg <- if (length(formals(fun)) > 0L) names(formals(fun))[1L] else "x"

    res_list <- lapply(seq_along(raw_res), function(i) {
      df_i <- raw_res[[i]]
      if (nrow(df_i) == 0L) {
        return(data.frame())
      }

      val_cols <- setdiff(names(df_i), c("cell", "coverage_area", "row", "col"))

      if (first_arg %in% c("df", "data", "data_frame", "x_df")) {
        res_i <- fun(df_i, ...)
      } else {
        input_vals <- if (length(val_cols) == 1L) df_i[[val_cols[1L]]] else as.matrix(df_i[, val_cols, drop = FALSE])
        res_i <- fun(input_vals, ...)
      }

      if (is.numeric(res_i) && !is.data.frame(res_i)) {
        res_i <- as.data.frame(as.list(res_i))
      }
      res_i
    })

    df_res <- dplyr::bind_rows(res_list)

    res_out <- if (return == "shapefile") {
      if (inherits(orig_shapefile, "SpatVector")) {
        res_sf <- dplyr::bind_cols(sf::st_as_sf(orig_shapefile), df_res)
        terra::vect(res_sf)
      } else if (inherits(orig_shapefile, "sf")) {
        dplyr::bind_cols(orig_shapefile, df_res)
      } else {
        dplyr::bind_cols(sf_poly, df_res)
      }
    } else {
      df_res
    }

    if (verbose) {
      cli::cli_progress_done()
      help_print_time(t_start)
    }

    return(res_out)
  }

  valid_funs <- c("mean", "median", "sum", "min", "max", "sd", "count", "quantile", "quantiles", "none")

  if (is.null(fun)) {
    fun_vec <- "none"
  } else {
    fun_raw <- tolower(as.character(fun))

    invalid_funs <- setdiff(fun_raw, valid_funs)
    if (length(invalid_funs) > 0L) {
      cli::cli_abort(c(
        "x" = "Invalid summary statistic{?s} specified in {.arg fun}: {.val {invalid_funs}}.",
        "i" = "Supported summary statistics are: {.val {valid_funs[valid_funs != 'none']}}."
      ))
    }
    fun_vec <- fun_raw
    if (!all(fun_vec == "none")) {
      fun_vec <- fun_vec[fun_vec != "none"]
    }
  }

  if (!is.numeric(summarize_quantiles) || any(summarize_quantiles < 0 | summarize_quantiles > 1)) {
    cli::cli_abort("Argument {.arg summarize_quantiles} must be a numeric vector with values between 0 and 1.")
  }

  if (verbose) {
    cli::cli_progress_step("Parsing polygon geometries...")
  }
  parsed_geoms <- help_parse_sf(sf_poly)

  e <- terra::ext(mosaic)
  r_xmin <- e[1]
  r_xmax <- e[2]
  r_ymin <- e[3]
  r_ymax <- e[4]
  r_bbox <- c(r_xmin, r_xmax, r_ymin, r_ymax)

  n_row <- terra::nrow(mosaic)
  n_col <- terra::ncol(mosaic)
  n_lyr <- terra::nlyr(mosaic)
  n_cell <- as.numeric(n_row) * as.numeric(n_col) * as.numeric(n_lyr)

  lyr_names <- names(mosaic)
  if (is.null(lyr_names) || any(nchar(lyr_names) == 0)) {
    lyr_names <- paste0("lyr.", seq_len(n_lyr))
  }

  in_memory_already <- all(terra::inMemory(mosaic))

  if (in_memory_already) {
    in_mem_mode <- TRUE
  } else {
    # Dynamic RAM Monitoring: Check available system RAM & 80% safe threshold only for rasters on disk
    required_ram_bytes <- n_cell * 8
    required_ram_gb <- required_ram_bytes / (1024^3)

    avail_ram_bytes <- get_free_ram(unit = "B")
    avail_ram_gb <- avail_ram_bytes / (1024^3)
    safe_ram_threshold_bytes <- avail_ram_bytes * 0.80
    safe_ram_threshold_gb <- safe_ram_threshold_bytes / (1024^3)

    if (verbose) {
      cli::cli_alert_info("Mosaic size: {.val {n_row}} x {.val {n_col}} x {.val {n_lyr}} ({.val {n_cell}} values, ~{.val {round(required_ram_gb, 3)}} GB RAM)")
      cli::cli_alert_info("Available RAM: {.val {round(avail_ram_gb, 2)}} GB (80% safe threshold: {.val {round(safe_ram_threshold_gb, 2)}} GB)")
    }

    if (is.null(max_cells_in_memory)) {
      in_mem_mode <- (required_ram_bytes <= safe_ram_threshold_bytes)
    } else {
      in_mem_mode <- (n_cell <= max_cells_in_memory)
    }

    if (!in_mem_mode && verbose) {
      cli::cli_alert_info(
        "Raster size ({.val {n_cell}} values, ~{.val {round(required_ram_gb, 2)}} GB RAM) exceeds 80% of available system memory ({.val {round(safe_ram_threshold_gb, 2)}} GB of {.val {round(avail_ram_gb, 2)}} GB free)."
      )
      cli::cli_alert_info(
        "Processing raster in spatial blocks to optimize memory usage..."
      )
    }
  }

  # Helper for naming summary matrix columns
  format_summary_colnames <- function(fun_vec, lyr_names, summarize_quantiles, coverage_area = FALSE) {
    col_names <- c()
    n_lyr <- length(lyr_names)

    for (f in fun_vec) {
      if (f %in% c("quantile", "quantiles")) {
        q_str <- as.character(summarize_quantiles)
        for (q_s in q_str) {
          if (n_lyr == 1L) {
            col_names <- c(col_names, paste0("q", q_s))
          } else {
            col_names <- c(col_names, paste0(lyr_names, "_q", q_s))
          }
        }
      } else {
        if (n_lyr == 1L) {
          col_names <- c(col_names, f)
        } else {
          if (length(fun_vec) == 1L) {
            col_names <- c(col_names, lyr_names)
          } else {
            col_names <- c(col_names, paste0(lyr_names, "_", f))
          }
        }
      }
    }
    if (isTRUE(coverage_area)) {
      col_names <- c(col_names, "covered_area", "plot_area", "coverage")
    }
    return(col_names)
  }

  # MODO 1: Processing fully in memory
  if (in_mem_mode) {
    if (verbose) {
      msg_mem <- if (in_memory_already) "Fetching in-memory raster matrix..." else "Reading raster values into memory matrix..."
      cli::cli_progress_step(msg_mem)
    }
    val_mat <- terra::values(mosaic, mat = TRUE)

    if (verbose) {
      cli::cli_progress_step("Running single-pass extraction engine...")
    }
    res_cpp <- cpp_extract_raster(
      values = val_mat,
      n_row = as.integer(n_row),
      n_col = as.integer(n_col),
      n_lyr = as.integer(n_lyr),
      bbox_raster = r_bbox,
      geoms_r = parsed_geoms,
      fun = fun_vec,
      exact = isTRUE(exact),
      return_coverage_area = isTRUE(coverage_area),
      subdiv = as.integer(subdiv),
      summarize_quantiles = as.numeric(summarize_quantiles)
    )

    if (verbose) {
      cli::cli_progress_step("Formatting extracted results...")
    }

    if (fun_vec[1L] == "none") {
      # If custom layer names exist on mosaic, rename layer columns in res_cpp directly
      if (!is.null(lyr_names) && !identical(lyr_names, paste0("lyr.", seq_len(n_lyr)))) {
        cols <- if (isTRUE(coverage_area)) {
          c("cell", "coverage_area", "row", "col", lyr_names)
        } else {
          c("cell", "row", "col", lyr_names)
        }
        res <- lapply(res_cpp, function(df) {
          names(df) <- cols
          df
        })
      } else {
        res <- res_cpp
      }

      if (verbose) {
        cli::cli_progress_done()
        help_print_time(t_start)
      }
      return(res)
    } else {
      df_res <- as.data.frame(res_cpp)
      col_names <- format_summary_colnames(fun_vec, lyr_names, summarize_quantiles, coverage_area = coverage_area)
      colnames(df_res) <- col_names

      res_out <- if (return == "shapefile") {
        if (inherits(orig_shapefile, "SpatVector")) {
          res_sf <- dplyr::bind_cols(sf::st_as_sf(orig_shapefile), df_res)
          terra::vect(res_sf)
        } else if (inherits(orig_shapefile, "sf")) {
          dplyr::bind_cols(orig_shapefile, df_res)
        } else {
          dplyr::bind_cols(sf_poly, df_res)
        }
      } else {
        df_res
      }

      if (verbose) {
        cli::cli_progress_done()
        help_print_time(t_start)
      }

      return(res_out)
    }
  }

  # MODO 2: Processing in spatial blocks to protect RAM
  if (verbose) {
    cli::cli_progress_step("Clustering polygons into spatial blocks...")
  }
  block_clusters <- help_block_cluster(sf_poly, max_per_block = 500L)
  n_blocks <- length(block_clusters)
  n_feat <- length(parsed_geoms)
  res_vec <- terra::res(mosaic)

  if (verbose) {
    cli::cli_progress_step("Processing {.val {n_blocks}} spatial block{?s} for raster with {.val {n_cell}} cells...")
    pb <- cli::cli_progress_bar(
      name = "Extracting spatial blocks",
      total = n_blocks
    )
  }

  if (fun_vec[1L] != "none") {
    col_names <- format_summary_colnames(fun_vec, lyr_names, summarize_quantiles, coverage_area = coverage_area)
    res_mat <- matrix(NA_real_, nrow = n_feat, ncol = length(col_names))

    for (k in seq_along(block_clusters)) {
      idx_k <- block_clusters[[k]]
      sf_k <- sf_poly[idx_k, ]
      bbox_k <- sf::st_bbox(sf_k)

      crop_box_k <- terra::ext(
        max(r_xmin, bbox_k$xmin - 2 * res_vec[1]),
        min(r_xmax, bbox_k$xmax + 2 * res_vec[1]),
        max(r_ymin, bbox_k$ymin - 2 * res_vec[2]),
        min(r_ymax, bbox_k$ymax + 2 * res_vec[2])
      )

      sub_raster_k <- terra::crop(mosaic, crop_box_k, snap = "out")
      e_k <- terra::ext(sub_raster_k)
      r_bbox_k <- c(e_k[1], e_k[2], e_k[3], e_k[4])

      val_mat_k <- terra::values(sub_raster_k, mat = TRUE)

      res_cpp_k <- cpp_extract_raster(
        values = val_mat_k,
        n_row = as.integer(terra::nrow(sub_raster_k)),
        n_col = as.integer(terra::ncol(sub_raster_k)),
        n_lyr = as.integer(n_lyr),
        bbox_raster = r_bbox_k,
        geoms_r = parsed_geoms[idx_k],
        fun = fun_vec,
        exact = isTRUE(exact),
        return_coverage_area = isTRUE(coverage_area),
        subdiv = as.integer(subdiv),
        summarize_quantiles = as.numeric(summarize_quantiles)
      )

      res_mat[idx_k, ] <- res_cpp_k

      if (verbose) {
        cli::cli_progress_update(id = pb)
      }
    }

    if (verbose) {
      cli::cli_progress_done(id = pb)
    }

    if (verbose) {
      cli::cli_progress_step("Formatting extracted results...")
    }

    df_res <- as.data.frame(res_mat)
    colnames(df_res) <- col_names

    res_out <- if (return == "shapefile") {
      if (inherits(orig_shapefile, "SpatVector")) {
        res_sf <- dplyr::bind_cols(sf::st_as_sf(orig_shapefile), df_res)
        terra::vect(res_sf)
      } else if (inherits(orig_shapefile, "sf")) {
        dplyr::bind_cols(orig_shapefile, df_res)
      } else {
        dplyr::bind_cols(sf_poly, df_res)
      }
    } else {
      df_res
    }

    if (verbose) {
      cli::cli_progress_done()
      help_print_time(t_start)
    }

    return(res_out)
  } else {
    res_list <- vector("list", n_feat)

    for (k in seq_along(block_clusters)) {
      idx_k <- block_clusters[[k]]
      sf_k <- sf_poly[idx_k, ]
      bbox_k <- sf::st_bbox(sf_k)

      crop_box_k <- terra::ext(
        max(r_xmin, bbox_k$xmin - 2 * res_vec[1]),
        min(r_xmax, bbox_k$xmax + 2 * res_vec[1]),
        max(r_ymin, bbox_k$ymin - 2 * res_vec[2]),
        min(r_ymax, bbox_k$ymax + 2 * res_vec[2])
      )

      sub_raster_k <- terra::crop(mosaic, crop_box_k, snap = "out")
      e_k <- terra::ext(sub_raster_k)
      r_bbox_k <- c(e_k[1], e_k[2], e_k[3], e_k[4])

      val_mat_k <- terra::values(sub_raster_k, mat = TRUE)

      res_cpp_k <- cpp_extract_raster(
        values = val_mat_k,
        n_row = as.integer(terra::nrow(sub_raster_k)),
        n_col = as.integer(terra::ncol(sub_raster_k)),
        n_lyr = as.integer(n_lyr),
        bbox_raster = r_bbox_k,
        geoms_r = parsed_geoms[idx_k],
        fun = fun_vec,
        exact = isTRUE(exact),
        return_coverage_area = isTRUE(coverage_area),
        subdiv = as.integer(subdiv),
        summarize_quantiles = as.numeric(summarize_quantiles)
      )

      if (!is.null(lyr_names) && !identical(lyr_names, paste0("lyr.", seq_len(n_lyr)))) {
        cols <- if (isTRUE(coverage_area)) {
          c("cell", "coverage_area", "row", "col", lyr_names)
        } else {
          c("cell", "row", "col", lyr_names)
        }
        for (j in seq_along(idx_k)) {
          df <- res_cpp_k[[j]]
          names(df) <- cols
          res_list[[idx_k[j]]] <- df
        }
      } else {
        for (j in seq_along(idx_k)) {
          res_list[[idx_k[j]]] <- res_cpp_k[[j]]
        }
      }

      if (verbose) {
        cli::cli_progress_update(id = pb)
      }
    }

    if (verbose) {
      cli::cli_progress_done(id = pb)
    }

    if (verbose) {
      cli::cli_progress_step("Formatting extracted results...")
    }

    if (verbose) {
      cli::cli_progress_done()
      help_print_time(t_start)
    }

    return(res_list)
  }
}

#' @rdname mosaic_extract
#' @export
extract_cpp <- mosaic_extract

# Helper to parse sf geometries into fast C++ geometry representation
help_parse_sf <- function(sf_poly) {
  geoms <- sf::st_geometry(sf_poly)
  n_feat <- length(geoms)
  res <- vector("list", n_feat)

  for (i in seq_len(n_feat)) {
    geom_i <- geoms[[i]]
    g_type <- as.character(sf::st_geometry_type(geom_i))

    poly_list <- if (g_type == "POLYGON") {
      list(geom_i)
    } else if (g_type == "MULTIPOLYGON") {
      as.list(geom_i)
    } else {
      list()
    }

    parts_res <- vector("list", length(poly_list))
    for (p in seq_along(poly_list)) {
      rings <- poly_list[[p]]
      outer_m <- as.matrix(rings[[1]])
      n_rings <- length(rings)
      holes_l <- vector("list", max(0, n_rings - 1))
      if (n_rings > 1) {
        for (h in 2:n_rings) {
          holes_l[[h - 1]] <- as.matrix(rings[[h]])
        }
      }
      parts_res[[p]] <- list(outer = outer_m, holes = holes_l)
    }
    res[[i]] <- list(parts = parts_res)
  }
  return(res)
}

# Helper to cluster polygons into spatial blocks for memory-safe extraction
help_block_cluster <- function(sf_poly, max_per_block = 500L) {
  n_feat <- nrow(sf_poly)
  if (n_feat <= max_per_block) {
    return(list(seq_len(n_feat)))
  }
  split(seq_len(n_feat), ceiling(seq_len(n_feat) / max_per_block))
}

# Internal helper to print formatted total processing time via CLI
help_print_time <- function(t_start) {
  t_diff <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  t_str <- if (t_diff < 1) paste0(round(t_diff * 1000, 0), "ms") else paste0(round(t_diff, 2), "s")
  cli::cli_progress_step("Extraction completed in {.val {t_str}}.")
}

# Internal helper to retrieve available system memory in bytes cross-platform via C++ (<0.01 ms)
help_get_free_ram_bytes <- function() {
  ram <- tryCatch(get_free_ram_cpp(), error = function(e) -1.0)
  if (!is.null(ram) && is.numeric(ram) && ram > 0) {
    return(ram)
  }

  # Fallback to 8 GB if RAM cannot be queried
  return(8 * 1024^3)
}
