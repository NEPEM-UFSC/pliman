% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_mosaic.R
\name{mosaic_analyze}
\alias{mosaic_analyze}
\title{Analyze a mosaic of remote sensing data}
\usage{
mosaic_analyze(
  mosaic,
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
  max_pixels = 2e+06,
  downsample = NULL,
  quantiles = c(0, 1),
  plot = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{mosaic}{A mosaic of class \code{SpatRaster}, generally imported with
\code{\link[=mosaic_input]{mosaic_input()}}.}

\item{r, g, b, re, nir, swir, tir}{The red, green, blue, red-edge,  near-infrared,
shortwave Infrared, and thermal infrared bands of the image, respectively.
By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
multispectral image is provided up to seven bands can be used to compute
built-in indexes. There are no limitation of band numbers if the index is
computed using the band name.}

\item{crop_to_shape_ext}{Crop the mosaic to the extension of shapefile?
Defaults to \code{TRUE}. This allows for a faster index computation when the
region of the built shapefile is much smaller than the entire mosaic
extension.}

\item{grid}{Logical, indicating whether to use a grid for segmentation
(default: TRUE).}

\item{nrow}{Number of rows for the grid (default: 1).}

\item{ncol}{Number of columns for the grid (default: 1).}

\item{plot_width, plot_height}{The width and height of the plot shape (in the
mosaic unit). It is mutually exclusiv with \code{buffer_col} and \code{buffer_row}.}

\item{layout}{Character: one of
\itemize{
\item \code{'tblr'} for top/bottom left/right orientation
\item \code{'tbrl'} for top/bottom right/left orientation
\item \code{'btlr'} for bottom/top left/right orientation
\item \code{'btrl'} for bottom/top right/left orientation
\item \code{'lrtb'} for left/right top/bottom orientation
\item \code{'lrbt'} for left/right bottom/top orientation
\item \code{'rltb'} for right/left top/bottom orientation
\item \code{'rlbt'} for right/left bottom/top orientation
}}

\item{indexes}{An optional \code{SpatRaster} object with the image indexes,
computed with \code{\link[=mosaic_index]{mosaic_index()}}.}

\item{shapefile}{An optional shapefile containing regions of interest (ROIs)
for analysis.}

\item{basemap}{An optional basemap generated with \code{\link[=mosaic_view]{mosaic_view()}}.}

\item{build_shapefile}{Logical, indicating whether to interactively draw ROIs
if the shapefile is \code{NULL} (default: TRUE).}

\item{check_shapefile}{Logical, indicating whether to validate the shapefile
with an interactive map view (default: TRUE). This enables live editing of
the drawn shapefile by deleting or changing the drawn grids.}

\item{buffer_edge}{Width of the buffer around the shapefile (default: 5).}

\item{buffer_col, buffer_row}{Buffering factor for the columns and rows,
respectively, of each individual plot's side. A value between 0 and 0.5
where 0 means no buffering and 0.5 means complete buffering (default: 0). A
value of 0.25 will buffer the plot by 25\% on each side.}

\item{segment_plot}{Logical, indicating whether to segment plots (default:
FALSE). If \code{TRUE}, the \code{segment_index} will be computed, and pixels with
values below the \code{threshold} will be selected.}

\item{segment_individuals}{Logical, indicating whether to segment individuals
within plots (default: FALSE). If \code{TRUE}, the \code{segment_index} will be
computed, and pixels with values below the \code{threshold} will be selected, and
a watershed-based segmentation will be performed.}

\item{segment_pick}{When \code{segment_plot} or \code{segment_individuals} are \code{TRUE},
\code{segment_pick} allows segmenting background (eg., soil) and foreground
(eg., plants) interactively by picking samples from background and
foreground using \code{\link[=mosaic_segment_pick]{mosaic_segment_pick()}}}

\item{mask}{An optional mask (SpatRaster) to mask the mosaic.}

\item{dsm}{A SpatRaster object representing the digital surface model. Must
be a single-layer raster. If a DSM is informed, a mask will be derived from
it using \code{\link[=mosaic_chm_mask]{mosaic_chm_mask()}}.}

\item{dsm_lower}{A numeric value specifying the lower height threshold. All
heights greater than this value are retained.}

\item{dsm_upper}{An optional numeric value specifying the upper height
threshold. If provided, only heights between lower and upper are retained.}

\item{dsm_window_size}{An integer (meters) specifying the window size (rows
and columns, respectively) for creating a DTM using a moving window.
Default is c(5, 5).}

\item{simplify}{Removes vertices in polygons to form simpler shapes. The
function implementation uses the Douglas–Peucker algorithm using
\code{\link[sf:geos_unary]{sf::st_simplify()}} for simplification.}

\item{map_individuals}{If \code{TRUE}, the distance between objects within plots
is computed. The distance can be mapped either in the horizontal or vertical
direction. The distances, coefficient of variation (CV), and mean of
distances are then returned.}

\item{map_direction}{The direction for mapping individuals within plots.
Should be one of \code{"horizontal"} or \code{"vertical"} (default).}

\item{watershed}{If \code{TRUE} (default), performs watershed-based object
detection. This will detect objects even when they are touching one another.
If FALSE, all pixels for each connected set of foreground pixels are set to
a unique object. This is faster but is not able to segment touching
objects.}

\item{tolerance}{The minimum height of the object in the units of image
intensity between its highest point (seed) and the point where it contacts
another object (checked for every contact pixel). If the height is smaller
than the tolerance, the object will be combined with one of its neighbors,
which is the highest.}

\item{extension}{Radius of the neighborhood in pixels for the detection of
neighboring objects. A higher value smooths out small objects.}

\item{include_if}{Character vector specifying the type of intersection.
Defaults to "centroid" (individuals in which the centroid is included within
the drawn plot will be included in that plot). Other possible values include
\code{"covered"}, \code{"overlap"}, and \code{"intersect"}. See Details for a detailed
explanation of these intersecting controls.}

\item{plot_index}{The index(es) to be computed for the drawn plots. Either a
single vegetation index (e.g., \code{"GLAI"}), a vector of indexes (e.g.,
\code{c("GLAI", "NGRDI", "HUE")}), or a custom index based on the available
bands (e.g., \code{"(R-B)/(R+B)"}). See \code{\link[=pliman_indexes]{pliman_indexes()}} and \code{\link[=image_index]{image_index()}}
for more details.}

\item{segment_index}{The index used for segmentation. The same rule as
\code{plot_index}. Defaults to \code{NULL}}

\item{threshold}{By default (threshold = "Otsu"), a threshold value based on
Otsu's method is used to reduce the grayscale image to a binary image. If a
numeric value is provided, this value will be used as a threshold.}

\item{opening, closing, filter, erode, dilate}{\strong{Morphological operations (brush size)}
\itemize{
\item \code{dilate} puts the mask over every background pixel, and sets it to
foreground if any of the pixels covered by the mask is from the foreground.
\item \code{erode} puts the mask over every foreground pixel, and sets it to
background if any of the pixels covered by the mask is from the background.
\item \code{opening} performs an erosion followed by a dilation. This helps to
remove small objects while preserving the shape and size of larger objects.
\item \code{closing} performs a dilatation followed by an erosion. This helps to
fill small holes while preserving the shape and size of larger objects.
\item \code{filter} performs median filtering in the binary image. Provide a positive
integer > 1 to indicate the size of the median filtering. Higher values are
more efficient to remove noise in the background but can dramatically impact
the perimeter of objects, mainly for irregular perimeters such as leaves
with serrated edges.
}}

\item{lower_noise}{To prevent noise from affecting the image analysis, objects
with lesser than 10\% of the mean area of all objects are removed
(\code{lower_noise = 0.1}). Increasing this value will remove larger noises (such
as dust points), but can remove desired objects too. To define an explicit
lower or upper size, use the \code{lower_size} and \code{upper_size} arguments.}

\item{lower_size, upper_size}{Lower and upper limits for size for the image
analysis. Plant images often contain dirt and dust.  Upper limit is set to
\code{NULL}, i.e., no upper limit used. One can set a known area or use
\code{lower_size = 0} to select all objects (not advised). Objects that matches
the size of a given range of sizes can be selected by setting up the two
arguments. For example, if \code{lower_size = 120} and \code{upper_size = 140},
objects with size greater than or equal 120 and less than or equal 140 will
be considered.}

\item{topn_lower, topn_upper}{Select the top \code{n} objects based on its area.
\code{topn_lower} selects the \code{n} elements with the smallest area whereas
\code{topn_upper} selects the \code{n} objects with the largest area.}

\item{summarize_fun}{The function to compute summaries for the pixel values.
Defaults to "mean," i.e., the mean value of the pixels (either at a plot- or
individual-level) is returned.}

\item{summarize_quantiles}{quantiles to be computed when 'quantile' is on \code{summarize_fun}.}

\item{attribute}{The attribute to be shown at the plot when \code{plot} is \code{TRUE}. Defaults to the first \code{summary_fun} and first \code{segment_index}.}

\item{invert}{Logical, indicating whether to invert the mask. Defaults to
\code{FALSE}, i.e., pixels with intensity greater than the threshold values are
selected.}

\item{color_regions}{The color palette for regions (default:
rev(grDevices::terrain.colors(50))).}

\item{alpha}{opacity of the fill color of the raster layer(s).}

\item{max_pixels}{Maximum number of pixels to render in the map or plot
(default: 500000).}

\item{downsample}{Downsampling factor to reduce the number of pixels
(default: NULL). In this case, if the number of pixels in the image (width
x height) is greater than \code{max_pixels} a downsampling factor will be
automatically chosen so that the number of plotted pixels approximates the
\code{max_pixels}.}

\item{quantiles}{the upper and lower quantiles used for color stretching.}

\item{plot}{Logical, indicating whether to generate plots (default: TRUE).}

\item{verbose}{Logical, indicating whether to display verbose output
(default: TRUE).}
}
\value{
A list containing the following objects:
\itemize{
\item \code{result_plot}: The results at a plot level.
\item \code{result_plot_summ}: The summary of results at a plot level. When
\code{segment_individuals = TRUE}, the number of individuals, canopy coverage,
and mean values of some shape statistics such as perimeter, length, width,
and diameter are computed.
\item \code{result_individ}: The results at an individual level.
\item \code{map_plot}: An object of class \code{mapview} showing the plot-level results.
\item \code{map_individual}: An object of class \code{mapview} showing the individual-level
results.
\item \code{shapefile}: The generated shapefile, with the drawn grids/blocks.
}
}
\description{
This function analyzes a mosaic of remote sensing data (UVAs or satellite
imagery), extracting information from specified regions of interest (ROIs)
defined in a shapefile or interactively drawn on the mosaic. It allows
counting and measuring individuals (eg., plants), computing canopy coverage,
and statistical summaries (eg., mean, coefficient of variation) for
vegetation indices (eg, NDVI) at a block, plot, individual levels or even
extract the raw results at pixel level.
}
\details{
Since multiple blocks can be analyzed, the length of arguments \code{grid},
\code{nrow}, \code{ncol}, \code{buffer_edge}, , \code{buffer_col}, \code{buffer_row}, \code{segment_plot},
\verb{segment_i, ndividuals}, \code{includ_if}, \code{threshold}, \code{segment_index}, \code{invert},
\code{filter}, \code{threshold}, \code{lower_size}, \code{upper_size}, \code{watershed}, and
\code{lower_noise}, can be either an scalar (the same argument applied to all the
drawn blocks), or a vector with the same length as the number of drawn. In
the last, each block can be analyzed with different arguments.

When \code{segment_individuals = TRUE} is enabled, individuals are included within
each plot based on the \code{include_if} argument. The default value
(\code{'centroid'}) includes an object in a given plot if the centroid of that
object is within the plot. This makes the inclusion mutually exclusive (i.e.,
an individual is included in only one plot). If \code{'covered'} is selected,
objects are included only if their entire area is covered by the plot. On the
other hand, selecting \code{overlap} is the complement of \code{covered}; in other
words, objects that overlap the plot boundary are included. Finally, when
\code{intersect} is chosen, objects that intersect the plot boundary are included.
This makes the inclusion ambiguous (i.e., an object can be included in more
than one plot).
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
library(pliman)
url <- "https://github.com/TiagoOlivoto/images/raw/master/pliman/rice_field/rice_ex.tif"
mosaic <- mosaic_input(url)
# Draw a polygon (top left, top right, bottom right, bottom left, top left)
# include 8 rice lines and one column
res <-
 mosaic_analyze(mosaic,
                r = 1, g = 2, b = 3,
                segment_individuals = TRUE,     # segment the individuals
                segment_index = "(G-B)/(G+B-R)",# index for segmentation
                filter = 4,
                nrow = 8,
                map_individuals = TRUE)
# map with individual results
res$map_indiv
}
}
