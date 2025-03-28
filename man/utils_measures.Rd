% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_measures.R
\name{utils_measures}
\alias{utils_measures}
\alias{get_measures}
\alias{plot_measures}
\title{Utilities for object measures}
\usage{
get_measures(
  object,
  measure = NULL,
  id = NULL,
  dpi = NULL,
  sep = "\\\\_|-",
  verbose = TRUE,
  digits = 5
)

plot_measures(
  object,
  measure = "id",
  id = NULL,
  hjust = NULL,
  vjust = NULL,
  digits = 2,
  size = 0.9,
  col = "white",
  ...
)
}
\arguments{
\item{object}{An object computed with \code{\link[=analyze_objects]{analyze_objects()}}.}

\item{measure}{For \code{plot_measures()}, a character string; for
\code{get_measures()}, a two-sided formula, e.g., \code{measure = area ~ 100}
indicating the known value of object \code{id}. The right-hand side is the known
value and the left-hand side can be one of the following.
\itemize{
\item \code{area} The known area of the object.
\item \code{perimeter} The known perimeter of the object.
\item \code{radius_mean} The known radius of the object.
\item \code{radius_min} The known minimum radius of the object. If the object is a
square, then the \code{radius_min} of such object will be \code{L/2} where \code{L} is the
length of the square side.
\item \code{radius_max} The known maximum radius of the object. If the object is a
square, then the \code{radius_max} of such object according to the Pythagorean
theorem will be \verb{L x sqrt(2) / 2} where \code{L} is the length of the square side.
}}

\item{id}{An object in the image to indicate a known value.}

\item{dpi}{A known resolution of the image in DPI (dots per inch).}

\item{sep}{Regular expression to manage file names. The function combines in
the \code{merge} object the object measures (sum of area and mean of all the
other measures) of all images that share the same filename prefix, defined
as the part of the filename preceding the first hyphen (-) or underscore
(_) (no hyphen or underscore is required). For example, the measures of
images named \verb{L1-1.jpeg}, \verb{L1-2.jpeg}, and \verb{L1-3.jpeg} would be combined
into a single image information (L1). This feature allows the user to treat
multiple images as belonging to a single sample, if desired. Defaults to
\code{sep = "\\\\_|-"}.}

\item{verbose}{If \code{FALSE}, runs the code silently.}

\item{digits}{The number of significant figures. Defaults to \code{2.}}

\item{hjust, vjust}{A numeric value to adjust the labels horizontally and
vertically. Positive values will move labels to right (hjust) and top
(vjust). Negative values will move the labels to left and bottom,
respectively.}

\item{size}{The size of the text. Defaults to \code{0.9}.}

\item{col}{The color of the text. Defaults to \code{"white"}.}

\item{...}{Further arguments passed on to \code{\link[graphics:text]{graphics::text()}}.}
}
\value{
\itemize{
\item For \code{get_measures()}, if \code{measure} is informed, the pixel values will be
corrected by the value of the known object, given in the unit of the
right-hand side of \code{meae}. If \code{dpi} is informed, then all the measures
will be adjusted to the knosurwn \code{dpi}.
}
\itemize{
\item If applied to an object of class \code{anal_obj}, returns a data frame with the
object \code{id} and the (corrected) measures.
\itemize{
\item If applied to an object of class \code{anal_obj_ls}, returns a list of class
\code{measures_ls}, with two objects: (i) \code{results}, a data frame containing
the identification of each image (img) and object within each image (id);
and (ii) \code{summary} a data frame containing the values for each image. If
more than one object is detected in a given image, the number of objects
(\code{n}), total area (\code{area_sum}), mean area (\code{area_mean}) and the standard
deviation of the area (\code{area_sd}) will be computed. For the other measures
(perimeter and radius), the mean values are presented.
}
}
\itemize{
\item \code{plot_measures()} returns a \code{NULL} object, drawing the text according to
the x and y coordinates of the objects in \code{object}.
}
}
\description{
\itemize{
\item \code{get_measures()} computes object measures (area, perimeter, radius) by using
either a known resolution (dpi) or an object with known measurements.
\item \code{plot_measures()} draws the object measures given in an object to the
current plot. The object identification (\code{"id"}) is drawn by default.
}
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
library(pliman)
img <- image_pliman("objects_300dpi.jpg")
plot(img)
# Image with four objects with a known resolution of 300 dpi
# Higher square: 10 x 10 cm
# Lower square: 5 x 5 cm
# Rectangle: 4 x 2 cm
# Circle: 3 cm in diameter

# Count the objects using the blue band to segment the image
results <-
   analyze_objects(img,
                 index = "B",
                 lower_noise = 0.1)
plot_measures(results, measure = "id")

# Get object measures by declaring the known resolution in dots per inch
(measures <- get_measures(results, dpi = 300))

# Calculated diagonal of the object 1
# 10 * sqrt(2) = 14.14

# Observed diagonal of the object 1
measures[1, "radius_max"] * 2


# Get object measures by declaring the known area of object 1
get_measures(results,
             id = 1,
             area ~ 100)
}


}
\author{
Tiago Olivoto \email{tiagoolivoto@gmail.com}
}
