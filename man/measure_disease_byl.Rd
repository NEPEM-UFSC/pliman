% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measure_disease_byl.R
\name{measure_disease_byl}
\alias{measure_disease_byl}
\title{Performs plant disease measurements by leaf}
\usage{
measure_disease_byl(
  img,
  index = "B",
  index_lb = "B",
  index_dh = "NGRDI",
  lower_size = NULL,
  watershed = TRUE,
  invert = FALSE,
  fill_hull = FALSE,
  opening = c(10, 0),
  closing = c(0, 0),
  filter = c(0, 0),
  erode = c(0, 0),
  dilate = c(0, 0),
  threshold = "Otsu",
  extension = NULL,
  tolerance = NULL,
  object_size = "large",
  img_healthy = NULL,
  img_symptoms = NULL,
  plot = TRUE,
  save_image = FALSE,
  dir_original = NULL,
  dir_processed = NULL,
  pattern = NULL,
  parallel = FALSE,
  workers = NULL,
  show_features = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{img}{The image to be analyzed.}

\item{index}{A character value specifying the target mode for
conversion to binary to segment the leaves from background. Defaults to "B"
(blue). See \code{\link[=image_index]{image_index()}} for more details. Personalized indexes can be
informed as, e.g., \verb{index = "R*G/B}.}

\item{index_lb}{The index used to segment the foreground (e.g., leaf) from
the background. If not declared, the entire image area (pixels) will be
considered in the computation of the severity.}

\item{index_dh}{The index used to segment diseased from healthy tissues when
\code{img_healthy} and \code{img_symptoms} are not declared. Defaults to \code{"GLI"}. See
\code{\link[=image_index]{image_index()}} for more details.}

\item{lower_size}{To prevent dust from affecting object segmentation, objects
with lesser than \verb{10\%} of the mean of all objects are removed. . One can
set a known area or use \code{lower_limit = 0} to select all objects (not
advised).}

\item{watershed}{If \code{TRUE} (default) performs watershed-based object
detection. This will detect objects even when they are touching one other.
If \code{FALSE}, all pixels for each connected set of foreground pixels are set
to a unique object. This is faster but is not able to segment touching
objects.}

\item{invert}{Inverts the binary image if desired. This is useful to process
images with a black background. Defaults to \code{FALSE}. If \code{reference = TRUE}
is use, \code{invert} can be declared as a logical vector of length 2 (eg.,
\verb{invert = c(FALSE, TRUE}). In this case, the segmentation of objects and
reference from the foreground using \code{back_fore_index} is performed using the
default (not inverted), and the segmentation of objects from the reference
is performed by inverting the selection (selecting pixels higher than the
threshold).}

\item{fill_hull}{Fill holes in the binary image? Defaults to \code{FALSE}. This is
useful to fill holes in objects that have portions with a color similar to
the background. IMPORTANT: Objects touching each other can be combined into
one single object, which may underestimate the number of objects in an
image.}

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

\item{threshold}{The theshold method to be used.
\itemize{
\item By default (\code{threshold = "Otsu"}), a threshold value based
on Otsu's method is used to reduce the grayscale image to a binary image. If
a numeric value is informed, this value will be used as a threshold.
\item If \code{threshold = "adaptive"}, adaptive thresholding (Shafait et al. 2008)
is used, and will depend on the \code{k} and \code{windowsize} arguments.
\item If any non-numeric value different than \code{"Otsu"} and \code{"adaptive"} is used,
an iterative section will allow you to choose the threshold based on a
raster plot showing pixel intensity of the index.
}}

\item{extension}{Radius of the neighborhood in pixels for the detection of
neighboring objects. Higher value smooths out small objects.}

\item{tolerance}{The minimum height of the object in the units of image
intensity between its highest point (seed) and the point where it contacts
another object (checked for every contact pixel). If the height is smaller
than the tolerance, the object will be combined with one of its neighbors,
which is the highest.}

\item{object_size}{The size of the object. Used to automatically set up
\code{tolerance} and \code{extension} parameters. One of the following. \code{"small"}
(e.g, wheat grains), \code{"medium"} (e.g, soybean grains), \code{"large"}(e.g, peanut
grains), and \code{"elarge"} (e.g, soybean pods)`.}

\item{img_healthy}{A color palette of healthy tissues.}

\item{img_symptoms}{A color palette of lesioned tissues.}

\item{plot}{Show image after processing?}

\item{save_image}{Save the image after processing? The image is saved in the
current working directory named as \verb{proc_*} where \code{*} is the image name
given in \code{img}.}

\item{dir_original, dir_processed}{The directory containing the original and
processed images. Defaults to \code{NULL}. In this case, the function will
search for the image \code{img} in the current working directory. After
processing, when \code{save_image = TRUE}, the processed image will be also
saved in such a directory. It can be either a full path, e.g.,
\code{"C:/Desktop/imgs"}, or a subfolder within the current working directory,
e.g., \code{"/imgs"}.}

\item{pattern}{A pattern of file name used to identify images to be
processed. For example, if \code{pattern = "im"} all images that the name
matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
Providing any number as pattern (e.g., \code{pattern = "1"}) will select
images that are named as 1.-, 2.-, and so on.}

\item{parallel}{Processes the images asynchronously (in parallel) in separate
R sessions running in the background on the same machine. It may speed up
the processing time, especially when \code{pattern} is used is informed. The
number of sections is set up to 30\% of available cores.}

\item{workers}{A positive numeric scalar or a function specifying the maximum
number of parallel processes that can be active at the same time.}

\item{show_features}{If \code{TRUE} returnS the lesion features such as number,
area, perimeter, and radius. Defaults to \code{FALSE}.}

\item{verbose}{If \code{TRUE} (default) a summary is shown in the console.}

\item{...}{Additional arguments passed on to \code{\link[=measure_disease]{measure_disease()}}.}
}
\value{
\itemize{
\item A list with the following objects:
\itemize{
\item \code{severity} A data frame with the percentage of healthy and symptomatic
areas for each leaf in the image(s).
\item \code{shape},\code{statistics} If \code{show_features = TRUE} is used, returns the shape
(area, perimeter, etc.) for each lesion and a summary statistic of the
results.
}
}
}
\description{
Computes the percentage of symptomatic leaf area using color palettes or RGB
indexes \code{by} each \code{l}eaf of an image. This allows, for example, processing
replicates of the same treatment  and obtaining the results for each
replication with a single image. To do that, leaf samples are first splitten
with \code{\link[=object_split]{object_split()}} and then, \code{\link[=measure_disease]{measure_disease()}} is applied to the list of
leaves.
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
library(pliman)
img <- image_pliman("mult_leaves.jpg", plot = TRUE)
sev <-
 measure_disease_byl(img = img,
                     index_lb = "B",
                     index_dh = "NGRDI",
                     workers = 2)
sev$severity
}


}
