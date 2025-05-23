% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_shapefile.R
\name{shapefile_interpolate}
\alias{shapefile_interpolate}
\title{Interpolate values at specific points based on coordinates and a target variable}
\usage{
shapefile_interpolate(
  shapefile,
  z,
  x = "x",
  y = "y",
  interpolation = c("Kriging", "Tps"),
  verbose = FALSE
)
}
\arguments{
\item{shapefile}{An sf object containing the x, y, and target variable (z)
columns. It is highly recommended to use \code{shapefile_measures()} to obtain
this data.}

\item{z}{A string specifying the name of the column in the shapefile that
contains the target variable to be interpolated.}

\item{x}{A string specifying the name of the column containing x-coordinates.
Default is 'x'.}

\item{y}{A string specifying the name of the column containing y-coordinates.
Default is 'y'.}

\item{interpolation}{A character vector specifying the interpolation method.
Options are "Kriging" or "Tps".}

\item{verbose}{Logical; if TRUE, progress messages will be displayed.}
}
\value{
A vector of interpolated values at the specified points.
}
\description{
This function interpolates values at specified points using x, y coordinates and a target variable
from a shapefile. It supports "Kriging" and "Tps" interpolation methods.
}
