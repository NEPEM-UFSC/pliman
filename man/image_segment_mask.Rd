% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_mask.R
\name{image_segment_mask}
\alias{image_segment_mask}
\title{Segment an \code{Image} object using a brush mask}
\usage{
image_segment_mask(
  img,
  size,
  shape = "disc",
  rel_pos_x = 0.5,
  rel_pos_y = 0.5,
  type = c("binary", "shadow"),
  col_background = "white",
  plot = TRUE,
  ...
)
}
\arguments{
\item{img}{A \code{Image} object}

\item{size}{A numeric containing the size of the brush in pixels. This should
be an odd number; even numbers are rounded to the next odd one.}

\item{shape}{A character vector indicating the shape of the brush. Can be
\code{"box"}, \code{"disc"}, \code{"diamond"}, \code{"Gaussian"} or \code{"line"} Defaults to
\code{"disc"}.}

\item{rel_pos_x, rel_pos_y}{A relative position to include the brush in the
image. Defaults to 0.5. This means that the brush will be centered in the
original image. Smaller values move the brush toward the left and top,
respectively.}

\item{type}{Defines the type of the mask. By default, a binary mask is
applied. This results in white pixels in the original image that matches
the 0s pixels in the brush. If \code{type = "shadow"} is used, a shadow mask is produced}

\item{col_background}{Background color after image segmentation. Defaults to
\code{"white"}.}

\item{plot}{Plots the generated mask? Defaults to \code{TRUE}.}

\item{...}{Further arguments passed on to \code{\link[EBImage:morphology]{EBImage::makeBrush()}}.}
}
\value{
A color \code{Image} object
}
\description{
It combines \code{\link[=make_mask]{make_mask()}} and \code{\link[=make_brush]{make_brush()}} to segment an \code{Image} object
using a brush of desired size, shape, and position.
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
img <- image_pliman("soybean_touch.jpg")
plot(img)
image_segment_mask(img, size = 601)
image_segment_mask(img,
                   size = 401,
                   shape = "diamond",
                   rel_pos_x = 0,
                   rel_pos_y = 0,
                   type = "shadow")
}
}
