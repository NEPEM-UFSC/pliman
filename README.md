
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pliman <img src="man/figures/logo.png" align="right" height="140/"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version-ago/pliman)](https://CRAN.R-project.org/package=pliman)
[![Lifecycle:
stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-week/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-day/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![DOI](https://zenodo.org/badge/DOI/10.32614/CRAN.package.pliman.svg)](https://doi.org/10.32614/CRAN.package.pliman)

<!-- badges: end -->

The pliman package offers tools for both single and batch image
manipulation and analysis ([Olivoto,
2022](https://onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13803)),
including the quantification of leaf area, disease severity assessment
([Olivoto et al.,
2022](https://link.springer.com/article/10.1007/s40858-021-00487-5)),
object counting, extraction of image indexes, shape measurement, object
landmark identification, and Elliptical Fourier Analysis of object
outlines, as detailed by Claude
([2008](https://link.springer.com/book/10.1007/978-0-387-77789-4)). The
package also provides a comprehensive pipeline for generating shapefiles
with complex layouts and supports high-throughput phenotyping of RGB,
multispectral, and hyperspectral orthomosaics. This functionality
facilitates field phenotyping using UAV- or satellite-based imagery.

`pliman` also provides useful functions for image transformation,
binarization, segmentation, and resolution. Please visit the
[Examples](https://nepem-ufsc.github.io/pliman/reference/index.html)
page on the `pliman` website for detailed documentation of each
function.

# Installation

Install the latest stable version of `pliman` from
[CRAN](https://CRAN.R-project.org/package=pliman) with:

``` r
install.packages("pliman")
```

The development version of `pliman` can be installed from
[GitHub](https://github.com/nepem-ufsc/pliman) using the
[pak](https://github.com/r-lib/pak) package:

``` r
pak::pkg_install("nepem-ufsc/pliman")
```

*Note*: If you are a Windows user, you should also first download and
install the latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

# Analyze objects

The function `analyze_objects()` can be used to analyze objects such as
leaves, grains, pods, and pollen in an image. By default, all measures
are returned in pixel units. Users can [adjust the object
measures](https://nepem-ufsc.github.io/pliman/articles/analyze_objects.html#adjusting-object-measures)
with `get_measures()` provided that the image resolution (Dots Per Inch)
is known. Another option is to use a reference object in the image. In
this last case, the argument `reference` must be set to `TRUE`. There
are two options to identify the reference object:

1.  By its color, using the arguments `back_fore_index` and
    `fore_ref_index`  
2.  By its size, using the arguments `reference_larger` or
    `reference_smaller`

In both cases, the `reference_area` must be declared. Let’s see how to
analyze an image with flax grains containing a reference object
(rectangle with 2x3 cm). Here, we’ll identify the reference object by
its size; so, the final results in this case will be in metric units
(cm).

``` r
library(pliman)
img <- image_pliman("flax_grains.jpg")
flax <- 
  analyze_objects(img,
                  index = "GRAY",
                  reference = TRUE,
                  reference_larger = TRUE,
                  reference_area = 6,
                  marker = "point",
                  marker_size = 0.5,
                  marker_col = "red", # default is white
                  show_contour = FALSE) # default is TRUE
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
# summary statistics
flax$statistics
#        stat        value
# 1         n 2.680000e+02
# 2  min_area 3.606989e-02
# 3 mean_area 6.250403e-02
# 4  max_area 1.262446e-01
# 5   sd_area 8.047152e-03
# 6  sum_area 1.675108e+01
# 7  coverage 5.388462e-02
```

# Disease severity

## Using image indexes

To compute the percentage of symptomatic leaf area you can use the
`measure_disease()` function you can use an image index to segment the
entire leaf from the background and then separate the diseased tissue
from the healthy tissue. Alternatively, you can provide color palette
samples to the `measure_disease()` function. In this approach, the
function fits a general linear model (binomial family) to the RGB values
of the image. It then uses the color palette samples to segment the
lesions from the healthy leaf.

In the following example, we compute the symptomatic area of a soybean
leaf. The proportion of healthy and symptomatic areas is given as a
proportion of the total leaf area after segmenting the leaf from the
background (blue).

``` r
img <- image_pliman("sev_leaf.jpg")
# Computes the symptomatic area
sev <- 
  measure_disease(img = img,
                  index_lb = "B", # to remove the background
                  index_dh = "NGRDI", # to isolate the diseased area
                  threshold = c("Otsu", 0), # You can also use the Otsu algorithm in both indexes (default)
                  plot = TRUE)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
sev$severity
#    healthy symptomatic
# 1 92.63213    7.367868
```

## Interactive disease measurements

An alternative approach to measuring disease percentage is available
through the `measure_disease_iter()` function. This function offers an
interactive interface that empowers users to manually select sample
colors directly from the image. By doing so, it provides a highly
customizable analysis method.

One advantage of using `measure_disease_iter()` is the ability to
utilize the “mapview” viewer, which enhances the analysis process by
offering zoom-in options. This feature allows users to closely examine
specific areas of the image, enabling detailed inspection and accurate
disease measurement.

``` r
img <- image_pliman("sev_leaf.jpg", plot = TRUE)

# works only in an interactive section
measure_disease_iter(img, viewer = "mapview")
```

# Citation

``` r
citation("pliman")
Please, support this project by citing it in your publications!

  Olivoto T (2022). "Lights, camera, pliman! An R package for plant
  image analysis." _Methods in Ecology and Evolution_, *13*(4),
  789-798. doi:10.1111/2041-210X.13803
  <https://doi.org/10.1111/2041-210X.13803>.

Uma entrada BibTeX para usuários(as) de LaTeX é

  @Article{,
    title = {Lights, camera, pliman! An R package for plant image analysis},
    author = {Tiago Olivoto},
    year = {2022},
    journal = {Methods in Ecology and Evolution},
    volume = {13},
    number = {4},
    pages = {789-798},
    doi = {10.1111/2041-210X.13803},
  }
```

# Getting help

If you come across any clear bugs while using the package, please
consider filing a minimal reproducible example on
[github](https://github.com/TiagoOlivoto/pliman/issues). This will help
the developers address the issue promptly.

Suggestions and criticisms aimed at improving the quality and usability
of the package are highly encouraged. Your feedback is valuable in
making {pliman} even better!

# Code of Conduct

Please note that the pliman project is released with a [Contributor Code
of Conduct](https://nepem-ufsc.github.io/pliman/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
