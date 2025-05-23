% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_mosaic.R
\name{mosaic_index}
\alias{mosaic_index}
\title{Mosaic Index}
\usage{
mosaic_index(
  mosaic,
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
  workers = 1
)
}
\arguments{
\item{mosaic}{A mosaic of class \code{SpatRaster}, generally imported with
\code{\link[=mosaic_input]{mosaic_input()}}.}

\item{index}{A character value (or a vector of characters) specifying the
target mode for conversion to a binary image. Use \code{\link[=pliman_indexes_rgb]{pliman_indexes_rgb()}}
and \code{\link[=pliman_indexes_me]{pliman_indexes_me()}} to see the available RGB and multispectral
indexes, respectively. Users can also calculate their own index using  \verb{R, G, B, RE, NIR, SWIR, and TIR} bands (eg., \code{index = "R+B/G"}) or using the
names of the mosaic's layers (ex., "(band_1 + band_2) / 2").}

\item{r, g, b, re, nir, swir, tir}{The red, green, blue, red-edge,  near-infrared,
shortwave Infrared, and thermal infrared bands of the image, respectively.
By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
multispectral image is provided up to seven bands can be used to compute
built-in indexes. There are no limitation of band numbers if the index is
computed using the band name.}

\item{plot}{Plot the computed index? Defaults to \code{TRUE}.}

\item{in_memory}{Logical, indicating whether the indexes should be computed
in memory. Defaults to \code{TRUE}. In most cases, this is 2-3 times faster, but
errors can occur if \code{mosaic} is a large \code{SpatRaster}. If \code{FALSE}, raster
algebra operations are performed on temporary files.}

\item{workers}{numeric. The number of workers you want to use for parallel
processing when computing multiple indexes.}
}
\value{
An index layer extracted/computed from the mosaic raster.
}
\description{
Compute or extract an index layer from a multi-band mosaic raster.
}
\details{
This function computes or extracts an index layer from the input
mosaic raster based on the specified index name. If the index is not found
in the package's predefined index list (see \code{\link[=image_index]{image_index()}} for more
details), it attempts to compute the index using the specified band
indices. The resulting index layer is returned as an \code{SpatRaster} object.
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
library(pliman)
mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
names(mosaic)
elev2 <- mosaic_index(mosaic, "elevation * 5", plot = FALSE)
oldpar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))

mosaic_plot(mosaic)
mosaic_plot(elev2)

# return the original parameters
par(oldpar)
}

}
