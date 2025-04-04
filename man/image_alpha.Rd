% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_imagem.R
\name{image_alpha}
\alias{image_alpha}
\title{Add Alpha Layer to an RGB Image}
\usage{
image_alpha(img, mask)
}
\arguments{
\item{img}{An RGB image of class \code{Image} from the EBImage package. The image must be in RGB format (color mode 2).}

\item{mask}{A numeric value or matrix/array specifying the alpha layer:
* If \code{mask} is a single numeric value, it sets a uniform transparency level (0 for fully transparent, 1 for fully opaque).
* If \code{mask} is a matrix or array, it must have the same dimensions as the image channels, allowing for varying transparency.}
}
\value{
An \code{Image} object with an added alpha layer, maintaining the RGBA format.
}
\description{
This function adds an alpha (transparency) layer to an RGB image using the EBImage package.
The alpha layer can be specified as a single numeric value for uniform transparency
or as a matrix/array matching the dimensions of the image for varying transparency.
}
\examples{
if (interactive() && requireNamespace("EBImage")) {
# Load the EBImage package
library(pliman)

# Load a sample RGB image
img <- image_pliman("soybean_touch.jpg")

# 50\% transparency
image_alpha(img, 0.5) |> plot()

# transparent background
mask <- image_binary(img, "NB")[[1]]
img_tb <- image_alpha(img, mask)
plot(img_tb)

}

}
