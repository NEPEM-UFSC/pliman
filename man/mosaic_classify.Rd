% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_mosaic.R
\name{mosaic_classify}
\alias{mosaic_classify}
\title{Classify a Mosaic Based on Index Breaks}
\usage{
mosaic_classify(mosaic, breaks, frequency = TRUE, plot = TRUE)
}
\arguments{
\item{mosaic}{A \code{SpatRaster} object representing the mosaic to be classified.}

\item{breaks}{A numeric vector specifying the breakpoints for classification.}

\item{frequency}{Logical. If \code{TRUE}, computes the class frequency and area (in hectares).}

\item{plot}{Logical. If \code{TRUE}, plots the classified mosaic.}
}
\value{
A list with two elements:
\itemize{
\item \code{classified}: A \code{SpatRaster} object containing the classified mosaic.
\item \code{class_freq}: A data frame containing class frequencies, areas (ha), and percentages (if \code{frequency = TRUE}).
}
}
\description{
This function classifies a given raster mosaic based on user-defined breaks.
It provides an option to calculate the frequency and area of each class, as
well as plot the classified mosaic.
}
\examples{
if(interactive()){
library(pliman)
library(terra)

# Create an example raster
r <- terra::rast(matrix(runif(100, min = 0, max = 1), nrow=10, ncol=10))

# Classify the raster
result <- mosaic_classify(r, breaks = c(0.3, 0.6))

# View results
result$classified
result$class_freq
}

}
