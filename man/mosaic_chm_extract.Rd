% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_mosaic.R
\name{mosaic_chm_extract}
\alias{mosaic_chm_extract}
\title{Extracts height metrics and plot quality from a Canopy Height Model (CHM)}
\usage{
mosaic_chm_extract(
  chm,
  shapefile,
  chm_threshold = NULL,
  plot_quality = c("absolute", "relative")
)
}
\arguments{
\item{chm}{An object computed with \code{\link[=mosaic_chm]{mosaic_chm()}}.}

\item{shapefile}{An \code{sf} object containing the polygons over which height
metrics are extracted.}

\item{chm_threshold}{A numeric value representing the height threshold for
calculating coverage. If \code{NULL}, coverage is not computed.}

\item{plot_quality}{A character string specifying whether plot quality should
be returned as \code{"absolute"} (raw Euclidean distance) or \code{"relative"}
(normalized between 0 and 1, where 1 indicates the best plot quality and 0
the worst). Defaults to \code{"absolute"}.}
}
\value{
An \code{sf} object containing height summary statistics for each plot, including:
\item{min}{Minimum height value.}
\item{q05}{5th percentile height value.}
\item{q50}{Median height value.}
\item{q95}{95th percentile height value.}
\item{max}{Maximum height value.}
\item{mean}{Mean height value.}
\item{cv}{Coefficient of variation (standard deviation divided by the mean).}
\item{entropy}{Shannon entropy of height values, representing height distribution complexity.}
\item{volume}{Total sum of heights multiplied by CHM resolution.}
\item{coverage}{Proportion of pixels exceeding \code{chm_threshold}. Only computed if \code{chm_threshold} is provided.}
\item{plot_quality}{Plot quality index. If \code{"relative"}, values are normalized between 0 and 1.}
}
\description{
This function extracts height-related summary statistics from a CHM using a given shapefile
and computes a plot quality metric based on coefficient of variation (CV), entropy, and coverage.
The plot quality can be returned in absolute or relative terms.
}
