#' @title Sample images
#' @description Sample images installed with the \pkg{pliman} package
#' @format `*.jpg` format
#' * `colorcheck.jpg` An image with RGB colorcheck, available at https://github.com/HarryCWright/PlantSizeClr
#' * `flax_leaves.jpg` Flax leaves in a white background
#' * `flax_grains.jpg` Flax grains with background light.
#' *  `la_back.jpg` A cyan palette representing the background of images
#' la_pattern, la_leaves, and soybean_touch.
#' * `la_leaf.jpg` A sample of the leaves in `la_leaves`
#' * `la_leaves.jpg` Tree leaves with a sample of known area.
#' * `mult_leaves.jpg` Three soybean leaflets with soybean rust symptoms.
#' * `objects_300dpi.jpg` An image with 300 dpi resolution.
#' * `potato_leaves.jpg` Three potato leaves, which were gathered from Gupta et
#' al. (2020).
#' * `sev_leaf.jpg` A soybean leaf with a blue background.
#' * `sev_leaf_nb.jpg` A soybean leaf without background.
#' * `sev_back.jpg` A blue palette representing the background of `sev_leaf`.
#' * `sev_healthy.jpg` Healthy area of `sev_leaf`.
#' * `sev_sympt.jpg` The symptomatic area `sev_leaf`.
#' * `shadow.jpg` A shaded leaf, useful to test adaptive thresholding
#' * `soy_green.jpg` Soybean grains with a white background.
#' * `soybean_grain.jpg` A sample palette of the grains in `soy_green`.
#' * `soybean_touch.jpg` Soybean grains with a cyan background touching one each
#' other.
#' * `field_mosaic.jpg` An UVA image from a soybean field.
#' @format `*.tif` format
#'
#' The following `.tif` files are provided as sample data, representing a slice
#' from a large orthomosaic with soybean plots in the vegetative stage. These
#' files were kindly provided by Arthur Bernardeli.
#' - **`ortho.tif`**: An orthomosaic with soybean plots (5 rows and 3 columns).
#' - **`dsm.tif`**: A digital surface model (DSM) for the soybean plots.
#' - **`dtm.tif`**: A digital terrain model (DTM) for the area.
#' - **`mask.tif`**: A mask that represents the soybean plants.
#' @md
#' @source Personal data, Gupta et al. (2020).
#' @references Gupta, S., Rosenthal, D. M., Stinchcombe, J. R., & Baucom, R. S.
#'   (2020). The remarkable morphological diversity of leaf shape in sweet
#'   potato (Ipomoea batatas): the influence of genetics, environment, and G×E.
#'   New Phytologist, 225(5), 2183–2195. \doi{10.1111/NPH.16286}
#' @name pliman_images
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords images
NULL


.onAttach <- function(libname, pkgname) {
  vers <- as.character(utils::packageVersion("pliman"))
  show_startup <- interactive() || !is.null(getOption("knitr.in.progress"))
  if(show_startup){
    fmt <- cli::ansi_columns(
      c(
        cli::format_inline("{.strong Developed collaboratively by} NEPEM {.url https://nepemufsc.com}"),
        cli::format_inline("{.strong Group lead:} Prof. Tiago Olivoto"),
        cli::format_inline("For citation, type {.code citation('pliman')}"),
        cli::format_inline("We welcome your feedback and suggestions!")
      ),
      width = 60,
      fill = "rows",
      align = "left"
    )
    packageStartupMessage(
      cli::boxx(fmt, ,
                header = cli::format_inline("{.cyan Welcome to pliman version {.val {vers}}!}"),
                footer = cli::format_inline("{.emph Simplifying high-throughput plant phenotyping in R}"),
                border_style = "round")
    )
  }
  check_ebi()
  check_mapview()
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Contorno", "display", "CODE", "dir_original" ,"dir_processada",
      "Spectrum", "value", "area", "id", ".", "object", "s.radius.max",
      "s.radius.min", "y", "s.area", "s.perimeter", "symptomatic", "m.eccentricity",
      "m.majoraxis", "s.radius.mean", "n_greater", "n_less", "setNames", "s.radius.sd",
      "perimeter", "radius_max", "radius_mean", "radius_min", "radius_sd", "X1",
      "X2", "i", "img", "plotn", "x", "leaf", "Band", "block", "mosaic",
      "geometry", "n", "area_sum", "individual", "compute_downsample", "plot_id",
      "re", "nir", "coverage_fraction", "sigma", "summarize_quantiles", "prop", "plot_id_seq",
      "B1", "B2", "B3", "cluster", "h", "s", "column", "data", "plot_area", "unique_id",
      "diam_max", "uuids", ".progress", "show"))
}


#' @title Contour outlines from five leaves
#' @description A list of contour outlines from five leaves. It may be used as
#'   example in some functions such as [efourier()]
#' @format A list with five objects
#' * `leaf_1`
#' * `leaf_2`
#' * `leaf_3`
#' * `leaf_4`
#' * `leaf_5`
#'
#' Each object is a `data.frame` with the coordinates for the outline perimeter
#' @md
#' @name contours
#' @source Personal data. The images were obtained in the Flavia data set
#'   downlodable at <https://flavia.sourceforge.net/>
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' Global option for controlling the viewer in pliman package
#'
#' Users can set the value of this option using `options("pliman_viewer", value)`.
#' The default value is "base". Use "mapview" to allow image to be
#' plotted/edited using the R packages mapview and mapedit
#'
#' @name pliman_viewer
#' @aliases pliman_viewer
#' @rdname pliman_viewer
NULL


#' Get the value of the pliman_viewer option
#'
#' Retrieves the current value of the pliman_viewer option used in the package.
#'
#' @return The current value of the pliman_viewer option.
#' @export
get_pliman_viewer <- function() {
  option_value <- options("pliman_viewer")[[1]]
  if (!is.null(option_value)) {
    option_value
  } else {
    "base"  # Default value
  }
}

#' Set the value of the pliman_viewer option
#'
#' Sets the value of the pliman_viewer option used in the package.
#'
#' @param value The value to be set for the pliman_viewer option.
#' @export
set_pliman_viewer <- function(value) {
  if(!value %in% c("base", "mapview")){
    cli::cli_abort("{.arg value} must be either {.val base} or {.val mapview}")
  }
  options(pliman_viewer = value)
}

### Utilities for defining viewer options
.onLoad <- function(libname, pkgname) {
  # Define as opções apenas se ainda não estiverem definidas
  op <- options()
  op.pliman <- list(
    pliman_viewer = "base",
    pliman_quiet = FALSE
  )
  toset <- !(names(op.pliman) %in% names(op))
  if (any(toset)) options(op.pliman[toset])
}


