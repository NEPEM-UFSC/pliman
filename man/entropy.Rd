% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{entropy}
\alias{entropy}
\title{Compute Shannon Entropy}
\usage{
entropy(x)
}
\arguments{
\item{x}{A numeric vector containing the values for which entropy will be computed.}
}
\value{
A numeric value representing the entropy.
}
\description{
This function calculates the Shannon entropy of a numeric vector.
}
\examples{
library(pliman)
x <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
entropy(x)
}
