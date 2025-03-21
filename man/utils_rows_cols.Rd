% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{utils_rows_cols}
\alias{utils_rows_cols}
\alias{column_to_rownames}
\alias{rownames_to_column}
\alias{remove_rownames}
\alias{round_cols}
\title{Utilities for handling with rows and columns}
\usage{
column_to_rownames(.data, var = "rowname")

rownames_to_column(.data, var = "rowname")

remove_rownames(.data)

round_cols(.data, digits = 2)
}
\arguments{
\item{.data}{A data frame}

\item{var}{Name of column to use for rownames.}

\item{digits}{The number of significant figures. Defaults to \code{2.}}
}
\description{
\itemize{
\item \code{columns_to_rownames()}: Move a column of \code{.data} to its row
names.
\item \code{rownames_to_column()}: Move the row names of \code{.data} to a new
column.
\item \code{remove_rownames()}: Remove the row names of \code{.data}.
\item \code{round_cols()} Rounds the values of all numeric variables to the specified
number of decimal places (default 2).
}
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
library(pliman)
iris2 <- iris |> rownames_to_column()
head(iris2)
iris2$rowname <- paste0("r", iris2$rowname)
iris2 |> column_to_rownames("rowname") |> head()
}

}
\author{
Tiago Olivoto \email{tiagoolivoto@gmail.com}
}
