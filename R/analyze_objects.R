#'Analyzes objects in an image
#'
#'@description
#' * [analyze_objects()] provides tools for counting and extracting object
#'features (e.g., area, perimeter, radius, pixel intensity) in an image. See
#'more at the **Details** section.
#' * [analyze_objects_iter()] provides an iterative section to measure object
#'features using an object with a known area.
#' * [plot.anal_obj()] produces a histogram for the R, G, and B values when
#'argument `object_index` is used in the function [analyze_objects()].
#'
#'@details A binary image is first generated to segment the foreground and
#'  background. The argument `index` is useful to choose a proper index to
#'  segment the image (see [image_binary()] for more details). It is also
#'  possible to provide color palettes for background and foreground (arguments
#'  `background` and `foreground`, respectively). When this is used, a general
#'  linear model (binomial family) fitted to the RGB values to segment fore- and
#'  background.
#'
#'  Then, the number of objects in the  foreground is counted. By setting up
#'  arguments such as `lower_size` and `upper_size`, it is possible to set a
#'  threshold for lower and upper sizes of the objects, respectively. The
#'  argument `object_size` can be used to set up pre-defined values of
#'  `tolerance` and `extension` depending on the image resolution. This will
#'  influence the watershed-based object segmentation. Users can also tune up
#'  `tolerance` and `extension` explicitly for a better precision of watershed
#'  segmentation.
#'
#'  If `watershed = FALSE` is used, all pixels for each connected set of
#'  foreground pixels in `img` are set to a unique object. This is faster,
#'  especially for a large number of objects, but it is not able to segment
#'  touching objects.
#'
#'  There are some ways to correct the measures based on a reference object. If
#'  a reference object with a known area (`reference_area`) is used in the image
#'  and `reference = TRUE` is used, the measures of the objects will be
#'  corrected, considering the unit of measure informed in `reference_area`.
#'  There are two main ways to work with reference objects.
#'    * The first, is to provide a reference object that has a contrasting color with
#'  both the background and object of interest. In this case, the arguments
#'  `back_fore_index` and `fore_ref_index` can be used to define an index to
#'  first segment the reference object and objects to be measured from the
#'  background, then the reference object from objects to be measured.
#'
#'
#'    * The second one is to use a reference object that has a similar color to the
#'  objects to be measured, but has a contrasting size. For example, if we are
#'  counting small brown grains, we can use a brown reference template that has
#'  an area larger (says 3 times the area of the grains) and then uses
#'  `reference_larger = TRUE`. With this, the larger object in the image will be
#'  used as the reference object. This is particularly useful when images are
#'  captured with background light, such as the example  2. Some types: (i) It
#'  is suggested that the reference object is not too much larger than the
#'  objects of interest (mainly when the `watershed = TRUE`). In some cases, the
#'  reference object can be broken into several pieces due to the watershed
#'  algorithm. (ii) Since the reference object will increase the mean area of
#'  the object, the argument `lower_noise` can be increased. By default
#'  (`lower_noise = 0.1`) objects with lesser than 10% of the mean area of all
#'  objects are removed. Since the mean area will be increased, increasing
#'  `lower_noise` will remove dust and noises more reliably. The argument
#'  `reference_smaller` can be used in the same way
#'
#'  By using `pattern`, it is possible to process several images with common
#'  pattern names that are stored in the current working directory or in the
#'  subdirectory informed in `dir_original`. To speed up the computation time,
#'  one can set `parallel = TRUE`.
#'
#'  [analyze_objects_iter()] can be used to process several images using an
#'  object with a known area as a template. In this case, all the images in the
#'  current working directory that match the `pattern` will be processed. For
#'  each image, the function will compute the features for the objects and show
#'  the identification (id) of each object. The user only needs to inform which
#'  is the id of the known object. Then, given the `known_area`, all the
#'  measures will be adjusted. In the end, a data.frame with the adjusted
#'  measures will be returned. This is useful when the images are taken at
#'  different heights. In such cases, the image resolution cannot be conserved.
#'  Consequently, the measures cannot be adjusted using the argument `dpi` from
#'  [get_measures()], since each image will have a different resolution. NOTE:
#'  This will only work in an interactive section.
#'
#' * Additional measures: By default, some measures are not computed, mainly due to
#'  computational efficiency when the user only needs simple measures such as
#'  area, length, and width.
#'
#'   - If `haralick = TRUE`, The function computes 13 Haralick texture features for
#'  each object based on a gray-level co-occurrence matrix (Haralick et al.
#'  1979). Haralick features depend on the configuration of the parameters
#'  `har_nbins` and `har_scales`. `har_nbins` controls the number of bins used
#'  to compute the Haralick matrix. A smaller `har_nbins` can give more accurate
#'  estimates of the correlation because the number of events per bin is higher.
#'  While a higher value will give more sensitivity. `har_scales` controls the
#'  number of scales used to compute the Haralick features. Since Haralick
#'  features compute the correlation of intensities of neighboring pixels it is
#'  possible to identify textures with different scales, e.g., a texture that is
#'  repeated every two pixels or 10 pixels. By default, the Haralick features
#'  are computed with the R band. To chance this default, use the argument
#'  `har_band`. For example, `har_band = 2` will compute the features with the
#'  green band. Additionaly, har_band = "GRAY" can be used. In this case, a
#'  grayscale (0.299 * R + 0.587 * G + 0.114 * B) is used.
#'
#'   - If `efourier = TRUE` is used, an Elliptical Fourier Analysis (Kuhl and
#'  Giardina, 1982) is computed for each object contour using [efourier()].
#'
#'    - If `veins = TRUE`, vein features are computed for each object using
#'  [object_veins()]. Veins are detected using Difference of Gaussians (DoG)
#'  bandpass filtering, dynamic resolution-independent boundary erosion to eliminate
#'  outer leaf contour artifacts, and per-object Otsu adaptive thresholding.
#'
#'     - If `ab_angles = TRUE` the apex and base angles of each object are
#'  computed with [poly_apex_base_angle()]. By default, the function computes
#'  the angle from the first pixel of the apex of the object to the two pixels
#'  that slice the object at the 25th percentile of the object height (apex
#'  angle). The base angle is computed in the same way but from the first base
#'  pixel.
#'
#'    - If `width_at = TRUE`, the width at the  5th, 25th, 50th, 75th, and 95th
#' percentiles of the object height are computed by default. These quantiles can
#' be adjusted with the `width_at_percentiles` argument.
#'
#' @inheritParams image_binary
#' @inheritParams image_index
#'
#'@param img The image to be analyzed.
#'@param foreground,background,reference_img A color palette for the foregrond,
#'  background, and reference object, respectively (optional). If a chacarceter
#'  is used (eg., `foreground = "fore"`), the function will search in the
#'  current working directory a valid image named "fore".
#' @param opening,closing,filter,erode,dilate **Morphological operations (brush
#'   size)**, applied in hierarchical order.
#'  * `opening` performs an erosion followed by a dilation. This helps to
#'   remove small objects while preserving the shape and size of larger objects.
#'  * `closing` performs a dilatation followed by an erosion. This helps to
#'   fill small holes while preserving the shape and size of larger objects.
#'  * `filter` performs median filtering in the binary image. Provide a positive
#'  integer > 1 to indicate the size of the median filtering. Higher values are
#'  more efficient to remove noise in the background but can dramatically impact
#'  the perimeter of objects, mainly for irregular perimeters such as leaves
#'  with serrated edges.
#'  * `erode` puts the mask over every foreground pixel, and sets it to
#'  background if any of the pixels covered by the mask is from the background.
#'  * `dilate` puts the mask over every background pixel, and sets it to
#'  foreground if any of the pixels covered by the mask is from the foreground.
#' @param pick_palettes  Logical argument indicating wheater the user needs to
#'   pick up the color palettes for foreground and background for the image. If
#'   `TRUE` [pick_palette()] will be called internally so that the user can sample
#'   color points representing foreground and background.
#' @param segment_objects Segment objects in the image? Defaults to `TRUE`. In
#'   this case, objects are segmented using the index defined in the `index`
#'   argument, and each object is analyzed individually. If `segment_objects =
#'   FALSE` is used, the objects are not segmented and the entire image is
#'   analyzed. This is useful, for example, when analyzing an image without
#'   background, where an `object_index` could be computed for the entire image,
#'   like the index of a crop canopy.
#' @param viewer The viewer option. This option controls the type of viewer to
#'   use for interactive plotting (eg., when `pick_palettes = TRUE`).  If not
#'   provided, the value is retrieved using [get_pliman_viewer()].
#'@param reference Logical to indicate if a reference object is present in the
#'  image. This is useful to adjust measures when images are not obtained with
#'  standard resolution (e.g., field images). See more in the details section.
#'@param reference_area The known area of the reference objects. The measures of
#'  all the objects in the image will be corrected using the same unit of the
#'  area informed here.
#'@param back_fore_index A character value to indicate the index to segment the
#'  foreground (objects and reference) from the background. Defaults to
#'  `"R/(G/B)"`. This index is optimized to segment white backgrounds from green
#'  leaves and a blue reference object.
#'@param fore_ref_index A character value to indicate the index to segment
#'  objects and the reference object. It can be either an available index in
#'  `pliman` (see [pliman_indexes()] or an own index computed with the R, G, and
#'  B bands. Defaults to `"B-R"`. This index is optimized to segment green
#'  leaves from a blue reference object after a white background has been
#'  removed.
#'@param reference_larger,reference_smaller Logical argument indicating when the
#'  larger/smaller object in the image must be used as the reference object.
#'  This only is valid when `reference` is set to `TRUE` and `reference_area`
#'  indicates the area of the reference object. IMPORTANT. When
#'  `reference_smaller` is used, objects with an area smaller than 1% of the
#'  mean of all the objects are ignored. This is used to remove possible noise
#'  in the image such as dust. So, be sure the reference object has an area that
#'  will be not removed by that cutpoint.
#'@param pattern A pattern of file name used to identify images to be imported.
#'  For example, if `pattern = "im"` all images in the current working directory
#'  that the name matches the pattern (e.g., img1.-, image1.-, im2.-) will be
#'  imported as a list. Providing any number as pattern (e.g., `pattern = "1"`)
#'  will select images that are named as 1.-, 2.-, and so on. An error will be
#'  returned if the pattern matches any file that is not supported (e.g.,
#'  img1.pdf).
#'@param parallel If `TRUE` processes the images asynchronously (in parallel) in
#'  separate R sessions running in the background on the same machine. It may
#'  speed up the processing time, especially when `pattern` is used is informed.
#'  When `object_index` is informed, multiple sections will be used to extract
#'  the RGB values for each object in the image. This may significantly speed up
#'  processing time when an image has lots of objects (say >1000).
#'@param workers A positive numeric scalar or a function specifying the number
#'  of parallel processes that can be active at the same time. By default, the
#'  number of sections is set up to 30% of available cores.
#'@param watershed If `TRUE` (default) performs watershed-based object
#'  detection. This will detect objects even when they are touching one other.
#'  If `FALSE`, all pixels for each connected set of foreground pixels are set
#'  to a unique object. This is faster but is not able to segment touching
#'  objects.
#'@param veins Logical argument indicating whether vein features are computed per object.
#'  Defaults to `FALSE`. If `TRUE`, detects veins using Difference of Gaussians (DoG)
#'  bandpass filtering and dynamic resolution-independent boundary erosion.
#'@param sigma_veins Fine Gaussian scale parameter for vein detection (default `1`).
#'  Coarse scale is automatically set to `sigma_veins * 5`.
#'@param veins_thinning Logical value indicating whether Guo-Hall skeleton thinning should
#'  be applied to detected veins (default `FALSE`). If `TRUE`, measures vein length
#'  skeleton density instead of vein area fraction.
#'@param veins_rel_erode Relative boundary erosion fraction of `min(width, height)` (default `0.01` = 1%).
#'  Trims leaf boundary pixels to eliminate outer edge artifacts independently of image resolution.
#' @param ab_angles  Logical argument indicating whether apex and base angles
#'   should be computed. Defaults to `FALSE`. If `TRUE`, `poly_apex_base_angle()`
#'   are called and the base and apex angles are computed considering the 25th
#'   and 75th percentiles of the object height. These percentiles can be changed
#'   with the argument `ab_angles_percentiles`.
#' @param ab_angles_percentiles The percentiles indicating the heights of the
#'   object for which the angle should be computed (from the apex and the
#'   bottom). Defaults to c(0.25, 0.75), which means considering the 25th and
#'   75th percentiles of the object height.
#' @param width_at Logical. If `TRUE`, the widths of the object at a given set
#'   of quantiles of the height are computed.
#' @param width_at_percentiles  A vector of heights along the vertical axis of
#'   the object at which the width will be computed. The default value is
#'   c(0.05, 0.25, 0.5, 0.75, 0.95), which means the function will return the
#'   width at the 5th, 25th, 50th, 75th, and 95th percentiles of the object's
#'   height.
#'@param haralick Logical value indicating whether Haralick features are
#'  computed. Defaults to `FALSE`.
#'@param har_nbins An integer indicating the number of bins using to compute the
#'  Haralick matrix. Defaults to 32. See Details
#'@param har_scales A integer vector indicating the number of scales to use to
#'  compute the Haralick features. See Details.
#'@param har_band The band to compute the Haralick features (1 = R, 2 = G, 3 =
#'  B). Defaults to 1. Other allowed value is `har_band = "GRAY"`.
#'@param smooth whether the object contours should be smoothed with
#'  [poly_smooth()]. Defaults to `FALSE`. To smooth use a numeric value
#'  indicating the number of interactions used to smooth the contours.
#' @param pcv Computes the Perimeter Complexity Value? Defaults to `FALSE`.
#' @param pcv_niter An integer specifying the number of smoothing iterations for
#'   computing the  Perimeter Complexity Value. Defaults to 100.
#'@param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'  numeric value of range 0-100 (proportion of the size of the original image).
#'@param trim Number of pixels removed from edges in the analysis. The edges of
#'  images are often shaded, which can affect image analysis. The edges of
#'  images can be removed by specifying the number of pixels. Defaults to
#'  `FALSE` (no trimmed edges).
#'@param fill_hull Fill holes in the binary image? Defaults to `FALSE`. This is
#'  useful to fill holes in objects that have portions with a color similar to
#'  the background. IMPORTANT: Objects touching each other can be combined into
#'  one single object, which may underestimate the number of objects in an
#'  image.
#'@param invert Inverts the binary image if desired. This is useful to process
#'  images with a black background. Defaults to `FALSE`. If `reference = TRUE`
#'  is use, `invert` can be declared as a logical vector of length 2 (eg.,
#'  `invert = c(FALSE, TRUE`). In this case, the segmentation of objects and
#'  reference from the foreground using `back_fore_index` is performed using the
#'  default (not inverted), and the segmentation of objects from the reference
#'  is performed by inverting the selection (selecting pixels higher than the
#'  threshold).
#'@param object_size The size of the object. Used to automatically set up
#'  `tolerance` and `extension` parameters. One of the following. `"small"`
#'  (e.g, wheat grains), `"medium"` (e.g, soybean grains), `"large"`(e.g, peanut
#'  grains), and `"elarge"` (e.g, soybean pods)`.
#'@param index A character value specifying the target mode for conversion to
#'  binary image when `foreground` and `background` are not declared. Defaults
#'  to `"NB"` (normalized blue). See [image_index()] for more details. User can
#'  also calculate your own index using the bands names, e.g. `index = "R+B/G"`
#'@param object_index Defaults to `FALSE`. If an index is informed, the average
#'  value for each object is returned. It can be the R, G, and B values or any
#'  operation involving them, e.g., `object_index = "R/B"`. In this case, it
#'  will return for each object in the image, the average value of the R/B
#'  ratio. Use [pliman_indexes_eq()] to see the equations of available indexes.
#' @param pixel_level_index Return the indexes computed in `object_index` in the
#'   pixel level? Defaults to `FALSE` to avoid returning large data.frames.
#' @param return_mask Returns the mask for the analyzed image? Defaults to `FALSE`.
#' @param return_exact Logical. If `TRUE`, exact morphological reconstruction via watershed-partitioned seed matching is performed when `opening > 0`. This preserves 100% of exact original shapes and boundaries of surviving objects while removing attached noise.
#' @param rel_size Relative size threshold as a fraction of the mean object area (e.g., `0.1` for 10% of mean object area). Setting `rel_size = 0.1` is equivalent to setting `lower_noise = 0.1`.
#' @param efourier Logical argument indicating if Elliptical Fourier should be
#'  computed for each object. This will call [efourier()] internally. It
#'  `efourier = TRUE` is used, both standard and normalized Fourier coefficients
#'  are returned.
#'@param nharm An integer indicating the number of harmonics to use. Defaults to
#'  10. For more details see [efourier()].
#'@param tolerance The minimum height of the object in the units of image
#'  intensity between its highest point (seed) and the point where it contacts
#'  another object (checked for every contact pixel). If the height is smaller
#'  than the tolerance, the object will be combined with one of its neighbors,
#'  which is the highest. Values < 1.0 are interpreted as relative tolerance
#'  (percentage of the peak height). Values >= 1.0 are interpreted as absolute
#'  tolerance in pixels (linear Euclidean distance).
#'@param extension Radius of the neighborhood in pixels for the detection of
#'  neighboring objects. Higher value smooths out small objects.
#'@param lower_noise To prevent noise from affecting the image analysis, objects
#'  with lesser than 10% of the mean area of all objects are removed
#'  (`lower_noise = 0.1`). Increasing this value will remove larger noises (such
#'  as dust points), but can remove desired objects too. To define an explicit
#'  lower or upper size, use the `lower_size` and `upper_size` arguments.
#'@param lower_size,upper_size Lower and upper limits for size for the image
#'  analysis. Plant images often contain dirt and dust.  Upper limit is set to
#'  `NULL`, i.e., no upper limit used. One can set a known area or use
#'  `lower_size = 0` to select all objects (not advised). Objects that matches
#'  the size of a given range of sizes can be selected by setting up the two
#'  arguments. For example, if `lower_size = 120` and `upper_size = 140`,
#'  objects with size greater than or equal 120 and less than or equal 140 will
#'  be considered.
#'@param topn_lower,topn_upper Select the top `n` objects based on its area.
#'  `topn_lower` selects the `n` elements with the smallest area whereas
#'  `topn_upper` selects the `n` objects with the largest area.
#'@param lower_eccent,upper_eccent,lower_circ,upper_circ Lower and upper limit
#'  for object eccentricity/circularity for the image analysis. Users may use
#'  these arguments to remove objects such as square papers for scale (low
#'  eccentricity) or cut petioles (high eccentricity) from the images. Defaults
#'  to `NULL` (i.e., no lower and upper limits).
#'@param randomize Randomize the lines before training the model?
#'@param nrows The number of lines to be used in training step. Defaults to
#'  2000.
#'@param plot Show image after processing?
#'@param show_original Show the count objects in the original image?
#'@param show_chull Show the convex hull around the objects? Defaults to
#'  `FALSE`.
#'@param show_contour Show a contour line around the objects? Defaults to
#'  `TRUE`.
#'@param show_bbox Show the bounding box around the objects? Defaults to `FALSE`.
#'@param contour_col,contour_size The color and size for the contour line around
#'  objects. Defaults to `contour_col = "red"` and `contour_size = 1`.
#'@param show_lw If `TRUE`, plots the length and width lines on each object
#'  calling [plot_lw()].
#'@param show_background Show the background? Defaults to `TRUE`. A white
#'  background is shown by default when `show_original = FALSE`.
#'@param show_segmentation Shows the object segmentation colored with random
#'  permutations. Defaults to `FALSE`.
#'@param col_foreground,col_background Foreground and background color after
#'  image processing. Defaults to `NULL`, in which `"black"`, and `"white"` are
#'  used, respectively.
#'@param marker,marker_col,marker_size The type, color and size of the object
#'  marker. Defaults to 'point', which plots a red point in each object. Use
#'  `marker = FALSE` to omit object marker.
#'@param save_image Save the image after processing? The image is saved in the
#'  current working directory named as `proc_*` where `*` is the image name
#'  given in `img`.
#'@param max_pixels Maximum number of pixels to render when plotting or saving
#'  the processed image. Defaults to `2e6`. Set to `NULL` for full rendering.
#'@param prefix The prefix to be included in the processed images. Defaults to
#'  `"proc_"`.
#'@param dir_original,dir_processed The directory containing the original and
#'  processed images. Defaults to `NULL`. In this case, the function will search
#'  for the image `img` in the current working directory. After processing, when
#'  `save_image = TRUE`, the processed image will be also saved in such a
#'  directory. It can be either a full path, e.g., `"C:/Desktop/imgs"`, or a
#'  subfolder within the current working directory, e.g., `"/imgs"`.
#'@param verbose If `TRUE` (default) a summary is shown in the console.
#'@param known_area The known area of the template object.
#'@param ... Depends on the function:
#' * For [analyze_objects_iter()], further arguments passed on to
#'  [analyze_objects()].
#'@return `analyze_objects()` returns a list with the following objects:
#'  * `results` A data frame with the following variables for each object in the
#'  image:
#'     - `id`:  object identification.
#'
#'     - `x`,`y`:  x and y coordinates for the center of mass of the object.
#'     - `area`:  area of the object (in pixels).
#'
#'     - `area_ch`:  the area of the convex hull around object (in pixels).
#'     - `perimeter`: perimeter (in pixels).
#'
#'     - `radius_min`, `radius_mean`, and `radius_max`: The minimum, mean, and
#'  maximum radius (in pixels), respectively.
#'
#'     - `radius_sd`: standard deviation of the mean radius (in pixels).
#'
#'     - `diam_min`, `diam_mean`, and `diam_max`: The minimum, mean, and
#'  maximum diameter (in pixels), respectively.
#'
#'     - `major_axis`, `minor_axis`: elliptical fit for major and minor axes (in
#'  pixels).
#'
#'     - `caliper`: The longest distance between any two points on the margin
#'  of the object. See [poly_caliper()] for more details
#'
#'     - `length`, `width` The length and width of objects (in pixels). These
#'  measures are obtained as the range of x and y coordinates after aligning
#'  each object with [poly_align()].
#'
#'     - `radius_ratio`: radius ratio given by `radius_max / radius_min`.
#'
#'     - `theta`: object angle (in radians).
#'
#'     - `eccentricity`: elliptical eccentricity computed using the
#'  ratio of the eigen values (inertia axes of coordinates).
#'
#'     - `form_factor` (Wu et al., 2007):  the difference between a leaf and a
#'  circle. It is defined as `4*pi*A/P`, where A is the area and P is the
#'  perimeter of the object.
#'
#'     - `narrow_factor` (Wu et al., 2007): Narrow factor (`caliper / length`).
#'
#'     - `asp_ratio` (Wu et al., 2007): Aspect ratio (`length / width`).
#'
#'     - `rectangularity` (Wu et al., 2007): The similarity between a leaf and
#'  a rectangle (`length * width/ area`).
#'
#'     - `pd_ratio` (Wu et al., 2007): Ratio of perimeter to diameter
#'  (`perimeter / caliper`)
#'
#'     - `plw_ratio` (Wu et al., 2007): Perimeter ratio of length and width
#'  (`perimeter / (length + width)`)

#'     - `solidity`: object solidity given by `area / area_ch`.
#'
#'     - `convexity`: The convexity of the object computed using the ratio
#'    between the perimeter of the convex hull and the perimeter of the polygon.
#'
#'     - `elongation`: The elongation of the object computed as `1 - width /
#'    length`.
#'
#'     - `circularity`: The object circularity given by `perimeter ^ 2 / area`.
#'
#'     - `circularity_haralick`: The Haralick's circularity (CH), computed as
#'     `CH =  m/sd`, where `m` and `sd` are the mean and standard deviations
#'     from each pixels of the perimeter to the centroid of the object.
#'
#'     - `circularity_norm`: The normalized circularity (Cn), to be unity for a
#'     circle. This measure is computed as `Cn = perimeter ^ 2 / 4*pi*area` and
#'     is invariant under translation, rotation, scaling transformations, and
#'     dimensionless.
#'
#'     - `asm`: The angular second-moment feature.
#'
#'     - `con`: The contrast feature
#'
#'     - `cor`: Correlation measures the linear dependency of gray levels of
#'     neighboring pixels.
#'
#'     - `var`: The variance of gray levels pixels.
#'
#'     - `idm`: The Inverse Difference Moment (IDM), i.e., the local
#'     homogeneity.
#'
#'     - `sav`: The Sum Average.
#'
#'     - `sva`: The Sum Variance.
#'
#'     - `sen`: Sum Entropy.
#'
#'     - `dva`: Difference Variance.
#'
#'     - `den`: Difference Entropy
#'
#'     - `f12`: Difference Variance.
#'
#'     - `f13`: The angular second-moment feature.
#'
#'  * `statistics`: A data frame with the summary statistics for the area of the
#'  objects.
#'  * `count`: If `pattern` is used, shows the number of objects in each image.
#'  * `obj_rgb`: If `object_index` is used, returns the R, G, and B values
#'  for each pixel of each object.
#'  * `object_index`: If `object_index` is used, returns the index computed for
#'  each object.
#'
#'  * Elliptical Fourier Analysis: If `efourier = TRUE` is used, the following
#'     objects are returned.
#'     - `efourier`: The Fourier coefficients.  For more details see
#'        [efourier()].
#'     - `efourier_norm`: The normalized Fourier coefficients. For more details
#'        see [efourier_norm()].
#'     - `efourier_error`: The error between original data and  reconstructed
#'        outline. For more details see [efourier_error()].
#'     - `efourier_power`: The spectrum of harmonic Fourier power.
#'        For more details see [efourier_power()].
#'
#'  * `veins`: If `veins = TRUE` is used, returns a data.frame containing the `id` and
#'  `prop_veins` (proportion of vein pixels per object area) for each object in the image.
#'
#'  * `analyze_objects_iter()` returns a data.frame containing the features
#'  described in the `results` object of [analyze_objects()].
#'
#'  * `plot.anal_obj()` returns a `trellis` object containing the distribution
#'  of the pixels, optionally  for each object when `facet = TRUE` is used.
#'
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#'
#' Gupta, S., Rosenthal, D. M., Stinchcombe, J. R., & Baucom, R. S. (2020). The
#' remarkable morphological diversity of leaf shape in sweet potato (Ipomoea
#' batatas): the influence of genetics, environment, and G×E. New Phytologist,
#' 225(5), 2183–2195. \doi{10.1111/NPH.16286}
#'
#' Haralick, R.M., K. Shanmugam, and I. Dinstein. 1973. Textural Features for Image
#' Classification. IEEE Transactions on Systems, Man, and Cybernetics SMC-3(6): 610–621.
#' \doi{10.1109/TSMC.1973.4309314}
#'
#' Kuhl, F. P., and Giardina, C. R. (1982). Elliptic Fourier features of a
#' closed contour. Computer Graphics and Image Processing 18, 236-258. doi:
#' \doi{10.1016/0146-664X(82)90034-X}
#'
#' Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area of a Polygon
#' and the Vector Cross Product. The Mathematics Teacher, 110(8), 631–636.
#' \doi{10.5951/mathteacher.110.8.0631}
#'
#' Montero, R. S., Bribiesca, E., Santiago, R., & Bribiesca, E. (2009). State
#' of the Art of Compactness and Circularity Measures. International
#' Mathematical Forum, 4(27), 1305–1335.
#'
#' Chen, C.H., and P.S.P. Wang. 2005. Handbook of Pattern Recognition and
#' Computer Vision. 3rd ed. World Scientific.
#'
#' Wu, S. G., Bao, F. S., Xu, E. Y., Wang, Y.-X., Chang, Y.-F., and Xiang, Q.-L.
#' (2007). A Leaf Recognition Algorithm for Plant Classification Using
#' Probabilistic Neural Network. in 2007 IEEE International Symposium on Signal
#' Processing and Information Technology, 11–16.
#' \doi{10.1109/ISSPIT.2007.4458016}
#'
#' @export
#' @name analyze_objects
#' @importFrom  utils install.packages
#' @importFrom grDevices col2rgb dev.off jpeg png rgb hcl.colors
#' @importFrom graphics lines par points rect text hist
#' @importFrom stats aggregate binomial glm kmeans predict sd runif dist var density
#' @importFrom utils menu
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' if (interactive()) {
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' obj <- analyze_objects(img)
#' obj$statistics
#'
#' ########################### Example 1 #########################
#' # Enumerate the objects in the original image
#' # Return the top-5 grains with the largest area
#' top <-
#'  analyze_objects(img,
#'                  marker = "id",
#'                  topn_upper = 5)
#' top$results
#'
#'
#' #' ########################### Example 1 #########################
#' # Correct the measures based on the area of the largest ob
#' ject
#' # note that since the reference object
#'
#' img <- image_pliman("flax_grains.jpg")
#' res <-
#'   analyze_objects(img,
#'                   index = "GRAY",
#'                   marker = "point",
#'                   show_contour = FALSE,
#'                   reference = TRUE,
#'                   reference_area = 6,
#'                   reference_larger = TRUE,
#'                   lower_noise = 0.3)
#' }
#'
analyze_objects <- function(img,
                            foreground = NULL,
                            background = NULL,
                            pick_palettes = FALSE,
                            segment_objects = TRUE,
                            viewer = get_pliman_viewer(),
                            reference = FALSE,
                            reference_area = NULL,
                            back_fore_index = "R/(G/B)",
                            fore_ref_index = "B-R",
                            reference_img = NULL,
                            reference_larger = FALSE,
                            reference_smaller = FALSE,
                            pattern = NULL,
                            parallel = FALSE,
                            workers = NULL,
                            watershed = TRUE,
                            veins = FALSE,
                            sigma_veins = 1,
                            veins_thinning = FALSE,
                            veins_rel_erode = 0.005,
                            ab_angles = FALSE,
                            ab_angles_percentiles = c(0.25, 0.75),
                            width_at = FALSE,
                            width_at_percentiles = c(0.05, 0.25, 0.50, 0.75, 0.95),
                            haralick = FALSE,
                            har_nbins = 32,
                            har_scales = 1,
                            har_band = 1,
                            smooth = FALSE,
                            pcv = FALSE,
                            pcv_niter = 100,
                            resize = FALSE,
                            trim = FALSE,
                            fill_hull = FALSE,
                            erode = FALSE,
                            dilate = FALSE,
                            opening = FALSE,
                            closing = FALSE,
                            filter = FALSE,
                            return_exact = FALSE,
                            filter_order = c("erode", "dilate", "opening", "closing", "filter", "fill_hull"),
                            invert = FALSE,
                            object_size = "medium",
                            area_mode = c("contour", "pixel"),
                            index = "NB",
                            r = 1,
                            g = 2,
                            b = 3,
                            re = 4,
                            nir = 5,
                            object_index = NULL,
                            pixel_level_index = FALSE,
                            return_mask = FALSE,
                            efourier = FALSE,
                            nharm = 10,
                            threshold = "Otsu",
                            k = 0.1,
                            windowsize = NULL,
                            tolerance = NULL,
                            extension = NULL,
                            lower_noise = 0.10,
                            rel_size = NULL,
                            lower_size = NULL,
                            upper_size = NULL,
                            topn_lower = NULL,
                            topn_upper = NULL,
                            lower_eccent = NULL,
                            upper_eccent = NULL,
                            lower_circ = NULL,
                            upper_circ = NULL,
                            randomize = TRUE,
                            nrows = 1000,
                            plot = TRUE,
                            show_original = TRUE,
                            show_chull = FALSE,
                            show_contour = FALSE,
                            show_bbox = FALSE,
                            contour_col = "red",
                            contour_size = 1,
                            show_lw = FALSE,
                            show_background = TRUE,
                            show_segmentation = FALSE,
                            col_foreground = NULL,
                            col_background = NULL,
                            marker = "point",
                            marker_col = "red",
                            marker_size = NULL,
                            save_image = FALSE,
                            max_pixels = 2e6,
                            prefix = "proc_",
                            dir_original = NULL,
                            dir_processed = NULL,
                            verbose = TRUE){

  check_filter_order(filter_order, verbose, erode, dilate, opening, closing, filter, fill_hull)
  if (!is.null(rel_size)) {
    lower_noise <- rel_size
  }
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)
  if (!object_size %in% c("small", "medium", "large", "elarge")) {
    cli::cli_abort("Argument {.arg object_size} must be one of {.val small}, {.val medium}, {.val large}, or {.val elarge}.")
  }

  if (!missing(img) && !missing(pattern)) {
    cli::cli_abort("Only one of {.arg img} or {.arg pattern} can be used.")
  }

  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  help_count <-
    function(img, foreground, background, pick_palettes, resize, fill_hull, threshold, erode, dilate, opening, closing, filter, filter_order,
             tolerance, extension, randomize, nrows, plot, show_original,
             show_background, marker, marker_col, marker_size, save_image,
             prefix, dir_original, dir_processed, verbose, col_background,
             col_foreground, lower_noise, ab_angles, ab_angles_percentiles, width_at, width_at_percentiles, return_mask, pcv,
             object_index, max_pixels = 4e6, return_exact = FALSE){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      if(trim != FALSE){
        if(!is.numeric(trim)){
          cli::cli_abort("Argument {.arg trim} must be numeric.")
        }
        img <- image_trim(img, trim)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          cli::cli_abort("Argument {.arg resize} must be numeric.")
        }
        img <- image_resize(img, resize)
      }
      # when reference is not used
      if(isFALSE(reference)){
        if(isTRUE(pick_palettes)){
          viewopt <- c("base", "mapview")
          viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]

          if(interactive()){
            if(viewopt == "base"){
              plot(img, max_pixels = max_pixels)
            }
            if(viewopt == "base"){
              cli::cli_alert_info("Use the first mouse button to pick up {.strong BACKGROUND} colors. Press ESC to exit.")
            }
            background <- pick_palette(img,
                                       r = 5,
                                       verbose = FALSE,
                                       palette  = FALSE,
                                       plot = FALSE,
                                       col = "blue",
                                       external_device = FALSE,
                                       title = "Use the first mouse button to pick up BACKGROUND colors. Click 'Done' to finish",
                                       viewer = viewer)
            if(viewopt != "base"){
              image_view(img[1:10, 1:10,], edit = TRUE)
            }
            if(viewopt == "base"){
              cli::cli_alert_info("Use the first mouse button to pick up {.strong FOREGROUND} colors. Press ESC to exit.")
            }
            foreground <- pick_palette(img,
                                       r = 5,
                                       verbose = FALSE,
                                       palette  = FALSE,
                                       plot = FALSE,
                                       col = "salmon",
                                       external_device = FALSE,
                                       title = "Use the first mouse button to pick up FOREGROUND colors. Click 'Done' to finish",
                                       viewer = viewer)
          }
        }
        if(!is.null(foreground) && !is.null(background)){
          if(is.character(foreground)){
            all_files <- sapply(list.files(getwd()), file_name)
            imag <- list.files(getwd(), pattern = foreground)
            check_names_dir(foreground, all_files, getwd())
            name <- file_name(imag)
            extens <- file_extension(imag)
            foreground <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
          }
          if(is.character(background)){
            all_files <- sapply(list.files(getwd()), file_name)
            imag <- list.files(getwd(), pattern = background)
            check_names_dir(background, all_files, getwd())
            name <- file_name(imag)
            extens <- file_extension(imag)
            background <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
          }
          img_num <- image_data(img, type = "numeric")
          fore_num <- image_data(foreground, type = "numeric")
          back_num <- image_data(background, type = "numeric")

          foreground <-
            data.frame(CODE = "foreground",
                       R = c(fore_num[,,1]),
                       G = c(fore_num[,,2]),
                       B = c(fore_num[,,3]))
          background <-
            data.frame(CODE = "background",
                       R = c(back_num[,,1]),
                       G = c(back_num[,,2]),
                       B = c(back_num[,,3]))
          back_fore <-
            transform(rbind(foreground[sample(1:nrow(foreground)),][1:nrows,],
                            background[sample(1:nrow(background)),][1:nrows,]),
                      Y = ifelse(CODE == "background", 0, 1))

          formula <- as.formula(paste("Y ~ ", back_fore_index))

          modelo1 <- suppressWarnings(glm(formula,
                                          family = binomial("logit"),
                                          data = back_fore))
          foreground_background <- predict_binary_glm(modelo1, img)
          for (op in filter_order) {
            if (op == "erode" && is.numeric(erode) && erode > 0) {
              foreground_background <- image_erode(foreground_background, size = erode, verbose = FALSE)
            } else if (op == "dilate" && is.numeric(dilate) && dilate > 0) {
              foreground_background <- image_dilate(foreground_background, size = dilate, verbose = FALSE)
            } else if (op == "opening" && is.numeric(opening) && opening > 0) {
              foreground_background <- image_opening(foreground_background, size = opening, verbose = FALSE)
            } else if (op == "closing" && is.numeric(closing) && closing > 0) {
              foreground_background <- image_closing(foreground_background, size = closing, verbose = FALSE)
            } else if (op == "filter" && !isFALSE(filter) && filter > 1) {
              foreground_background <- image_filter(foreground_background, size = filter, verbose = FALSE)
            } else if (op == "fill_hull" && isTRUE(fill_hull)) {
              foreground_background <- image_fill_hull(foreground_background, verbose = FALSE)
            }
          }


          ID <- c(foreground_background == 1)
          ID2 <- c(foreground_background == 0)
          if(isTRUE(watershed)){
            nmask <- image_watershed(foreground_background, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
          } else{
            nmask <- image_bwlabel(foreground_background)
          }

        } else{
          if(isTRUE(segment_objects)){
            img2 <- help_binary(img,
                                index = index,
                                r = r,
                                g = g,
                                b = b,
                                re = re,
                                nir = nir,
                                invert = invert,
                                fill_hull = fill_hull,
                                threshold = threshold,
                                k = k,
                                windowsize = windowsize,
                                erode = erode,
                                dilate = dilate,
                                opening = opening,
                                closing = closing,
                                filter = filter,
                                resize = FALSE,
                                filter_order = filter_order,
                                return_exact = return_exact)
            if(isTRUE(watershed)){
              nmask <- image_watershed(img2, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
            } else{
              nmask <- image_bwlabel(img2)
            }
          } else{
            img2 <- img[,,1]
            img2[image_data(img2) == 0 | image_data(img2) != 0] <- TRUE
            nmask <- image_bwlabel(img2)
          }

          ID <- which(img2 != 0)
          ID2 <- which(img2 == 0)
        }
        if(isTRUE(fill_hull)){
          nmask <- image_fill_hull(nmask)
        }
        shape <- compute_measures(mask = nmask,
                                  img = img,
                                  haralick = haralick,
                                  har_nbins = har_nbins,
                                  har_scales = har_scales,
                                  har_band = har_band,
                                  smooth = smooth)
        object_contour <- shape$cont
        # return(list(object_contour = object_contour, shape = shape))
        ch <- shape$ch
        shape <- shape$shape

      } else{
        # when reference is used
        if(is.null(reference_area)){
          cli::cli_abort("A known {.strong area} must be declared when a {.strong template} is used.")
        }
        if(isFALSE(reference_larger) & isFALSE(reference_smaller)){
          # segment back and fore
          if(!isFALSE(invert)){
            invert1 <- ifelse(length(invert) == 1, invert, invert[1])
          } else{
            invert1 <- FALSE
          }
          if(!is.null(background) & !is.null(foreground) & !is.null(reference_img)){
            if(is.character(foreground)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = foreground)
              check_names_dir(foreground, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              foreground <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            if(is.character(background)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = background)
              check_names_dir(background, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              background <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            if(is.character(reference_img)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = reference_img)
              check_names_dir(reference_img, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              reference_img <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            img_num <- image_data(img, type = "numeric")
            fore_num <- image_data(foreground, type = "numeric")
            ref_num <- image_data(reference_img, type = "numeric")
            back_num <- image_data(background, type = "numeric")

            original <-
              data.frame(CODE = "img",
                         R = c(img_num[,,1]),
                         G = c(img_num[,,2]),
                         B = c(img_num[,,3]))
            fore <-
              data.frame(CODE = "foreground",
                         R = c(fore_num[,,1]),
                         G = c(fore_num[,,2]),
                         B = c(fore_num[,,3]))
            ref <-
              data.frame(CODE = "reference",
                         R = c(ref_num[,,1]),
                         G = c(ref_num[,,2]),
                         B = c(ref_num[,,3]))
            back <-
              data.frame(CODE = "background",
                         R = c(back_num[,,1]),
                         G = c(back_num[,,2]),
                         B = c(back_num[,,3]))
            back_fore <-
              transform(rbind(fore[sample(1:nrow(fore)),][1:1000,],
                              ref[sample(1:nrow(ref)),][1:1000,],
                              back[sample(1:nrow(back)),][1:1000,]),
                        Y = ifelse(CODE == "background", 0, 1))

            formula <- as.formula(paste("Y ~ ", "R+G+B"))

            modelo1 <- suppressWarnings(glm(formula,
                                            family = binomial("logit"),
                                            data = back_fore))
            img_bf <- as_image(predict_binary_glm(modelo1, img))
            for (op in filter_order) {
              if (op == "erode" && is.numeric(erode) && erode > 0) {
                img_bf <- image_erode(img_bf, size = erode, verbose = FALSE)
              } else if (op == "dilate" && is.numeric(dilate) && dilate > 0) {
                img_bf <- image_dilate(img_bf, size = dilate, verbose = FALSE)
              } else if (op == "opening" && is.numeric(opening) && opening > 0) {
                img_bf <- image_opening(img_bf, size = opening, verbose = FALSE)
              } else if (op == "closing" && is.numeric(closing) && closing > 0) {
                img_bf <- image_closing(img_bf, size = closing, verbose = FALSE)
              } else if (op == "filter" && !isFALSE(filter) && filter > 1) {
                img_bf <- image_filter(img_bf, filter, verbose = FALSE)
              } else if (op == "fill_hull" && isTRUE(fill_hull)) {
                img_bf <- image_fill_hull(img_bf, verbose = FALSE)
              }
            }
          } else{
            img_bf <-
              help_binary(img,
                          threshold = threshold,
                          index = back_fore_index,
                          erode = erode,
                          dilate = dilate,
                          opening = opening,
                          closing = closing,
                          filter = filter,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          k = k,
                          windowsize = windowsize,
                          invert = invert1,
                          fill_hull = fill_hull,
                          filter_order = filter_order)
          }

          img3 <- img
          img3[,,1][which(img_bf != 1)] <- 2
          img3[,,2][which(img_bf != 1)] <- 2
          img3[,,3][which(img_bf != 1)] <- 2
          ID <-  which(img_bf == 1) # IDs for foreground
          ID2 <- which(img_bf == 0) # IDs for background
          # segment fore and ref
          if(!isFALSE(invert)){
            invert2 <- ifelse(length(invert) == 1, invert, invert[2])
          } else{
            invert2 <- FALSE
          }
          if(!is.null(background) & !is.null(foreground) & !is.null(reference_img)){
            back_fore <-
              transform(rbind(fore[sample(1:nrow(fore)),][1:1000,],
                              ref[sample(1:nrow(ref)),][1:1000,],
                              back[sample(1:nrow(back)),][1:1000,]),
                        Y = ifelse(CODE == "reference", 0, 1))

            formula <- as.formula(paste("Y ~ ", "R+G+B"))

            modelo1 <- suppressWarnings(glm(formula,
                                            family = binomial("logit"),
                                            data = back_fore))
            img4 <- as_image(predict_binary_glm(modelo1, img))
            for (op in filter_order) {
              if (op == "erode" && is.numeric(erode) && erode > 0) {
                img4 <- image_erode(img4, size = erode, verbose = FALSE)
              } else if (op == "dilate" && is.numeric(dilate) && dilate > 0) {
                img4 <- image_dilate(img4, size = dilate, verbose = FALSE)
              } else if (op == "opening" && is.numeric(opening) && opening > 0) {
                img4 <- image_opening(img4, size = opening, verbose = FALSE)
              } else if (op == "closing" && is.numeric(closing) && closing > 0) {
                img4 <- image_closing(img4, size = closing, verbose = FALSE)
              } else if (op == "filter" && !isFALSE(filter) && filter > 1) {
                img4 <- image_filter(img4, filter, verbose = FALSE)
              } else if (op == "fill_hull" && isTRUE(fill_hull)) {
                img4 <- image_fill_hull(img4, verbose = FALSE)
              }
            }
          } else{
            img4 <-
              help_binary(img3,
                          threshold = threshold,
                          index = fore_ref_index,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          erode = erode,
                          dilate = dilate,
                          opening = opening,
                          closing = closing,
                          filter = filter,
                          k = k,
                          windowsize = windowsize,
                          invert = invert2,
                          filter_order = filter_order)
          }

          mask <- img_bf
          pix_ref <- which(img4 != 1)
          img[,,1][pix_ref] <- 255
          img[,,2][pix_ref] <- 0
          img[,,3][pix_ref] <- 0
          npix_ref <- length(pix_ref)
          mask[pix_ref] <- 0
          if(is.numeric(filter) & filter > 1){
            mask <- image_filter(mask, size = filter)
          }
          if(isTRUE(watershed)){
            nmask <- image_watershed(mask, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
          } else{
            nmask <- image_bwlabel(mask)
          }
          shape <- compute_measures(mask = nmask,
                                    img = img,
                                    haralick = haralick,
                                    har_nbins = har_nbins,
                                    har_scales = har_scales,
                                    har_band = har_band,
                                    smooth = smooth)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape
          if(isTRUE(show_lw)){
            shape_ori <- shape
          }
          # correct measures based on the area of the reference object
          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side^2
          shape$area_ch <- shape$area_ch * px_side^2
          shape[6:18] <- apply(shape[6:18], 2, function(x){
            x * px_side
          })
        } else{
          # correct the measures based on larger or smaller objects

          if((!is.null(foreground) && !is.null(background)) | isTRUE(pick_palettes)){

            if(isTRUE(pick_palettes)){
              viewopt <- c("base", "mapview")
              viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]

              if(interactive()){
                if(viewopt == "base"){
                  plot(img, max_pixels = max_pixels)
                }
                if(viewopt == "base"){
                  cli::cli_alert_info("Use the first mouse button to pick up {.strong BACKGROUND} colors. Press ESC to exit.")
                }
                background <- pick_palette(img,
                                           r = 5,
                                           verbose = FALSE,
                                           palette  = FALSE,
                                           plot = FALSE,
                                           col = "blue",
                                           external_device = FALSE,
                                           title = "Use the first mouse button to pick up BACKGROUND colors. Click 'Done' to finish",
                                           viewer = viewer)
                if(viewopt != "base"){
                  image_view(img[1:10, 1:10,], edit = TRUE)
                }
                if(viewopt == "base"){
                  cli::cli_alert_info("Use the first mouse button to pick up {.strong FOREGROUND} colors. Press ESC to exit.")
                }
                foreground <- pick_palette(img,
                                           r = 5,
                                           verbose = FALSE,
                                           palette  = FALSE,
                                           plot = FALSE,
                                           col = "salmon",
                                           external_device = FALSE,
                                           title = "Use the first mouse button to pick up FOREGROUND colors. Click 'Done' to finish",
                                           viewer = viewer)
              }
            }

            if(is.character(foreground)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = foreground)
              check_names_dir(foreground, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              foreground <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            if(is.character(background)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = background)
              check_names_dir(background, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              background <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            foreground <-
              data.frame(CODE = "foreground",
                         R = c(as.numeric(image_data(foreground)[,,1])),
                         G = c(as.numeric(image_data(foreground)[,,2])),
                         B = c(as.numeric(image_data(foreground)[,,3])))
            background <-
              data.frame(CODE = "background",
                         R = c(as.numeric(image_data(background)[,,1])),
                         G = c(as.numeric(image_data(background)[,,2])),
                         B = c(as.numeric(image_data(background)[,,3])))
            back_fore <-
              transform(rbind(foreground[sample(1:nrow(foreground)),][1:nrows,],
                              background[sample(1:nrow(background)),][1:nrows,]),
                        Y = ifelse(CODE == "background", 0, 1))

            formula <- as.formula(paste("Y ~ ", back_fore_index))

            modelo1 <- suppressWarnings(glm(formula,
                                            family = binomial("logit"),
                                            data = back_fore))
            foreground_background <- predict_binary_glm(modelo1, img)
            if(is.numeric(filter) & filter > 1){
              foreground_background <- image_filter(foreground_background, size = filter)
            }
            ID <- c(foreground_background == 1)
            ID2 <- c(foreground_background == 0)
            if(isTRUE(watershed)){
              nmask <- image_watershed(foreground_background, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
            } else{
              nmask <- image_bwlabel(foreground_background)
            }

          } else{


            mask <-
              help_binary(img,
                          threshold = threshold,
                          index = index,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          k = k,
                          windowsize = windowsize,
                          erode = erode,
                          dilate = dilate,
                          opening = opening,
                          closing = closing,
                          filter = filter,
                          invert = invert,
                          fill_hull = fill_hull,
                          filter_order = filter_order)
            ID <-   which(mask != 0) # IDs for foreground
            ID2 <- which(mask == 0) # IDs for background
            if(isTRUE(watershed)){
              nmask <- image_watershed(mask, tolerance = tolerance, ext = ifelse(is.null(extension), 1, extension))
            } else{
              nmask <- image_bwlabel(mask)
            }

          }


          shape <- compute_measures(mask = nmask,
                                    img = img,
                                    haralick = haralick,
                                    har_nbins = har_nbins,
                                    har_scales = har_scales,
                                    har_band = har_band,
                                    smooth = smooth)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape

          if (area_mode[[1]] == "pixel") {
            pix_areas <- get_area_mask(image_data(nmask))
            shape$area <- pix_areas[shape$id]
            shape$coverage <- shape$area / length(nmask)
            shape$form_factor <- 4 * pi * shape$area / (shape$perimeter ^ 2)
            shape$rectangularity <- (shape$length * shape$width) / shape$area
            shape$solidity <- shape$area / shape$area_ch
            shape$circularity <- (shape$perimeter ^ 2) / shape$area
            shape$circularity_norm <- (shape$area * 4 * pi) / (shape$perimeter ^ 2)
          }


          r_val <- if (is.raw(img)) 255 else 1
          if(isTRUE(reference_larger)){
            lineid <- which.max(shape$area)
            id_ref <- shape[lineid, "id"]
            pix_ref <- which(image_data(nmask) == id_ref)
            img[,,1][pix_ref] <- r_val
            img[,,2][pix_ref] <- 0
            img[,,3][pix_ref] <- 0
            npix_ref <- shape[lineid, 4]
            shape <- shape[-lineid,]
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
          } else{
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
            lineid <- which.min(shape$area)
            id_ref <- shape[lineid, "id"]
            pix_ref <- which(image_data(nmask) == id_ref)
            img[,,1][pix_ref] <- r_val
            img[,,2][pix_ref] <- 0
            img[,,3][pix_ref] <- 0
            npix_ref <- shape[lineid, 4]
            shape <- shape[-lineid,]
          }
          if(isTRUE(show_lw)){
            shape_ori <- shape
          }

          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side ^ 2
          shape$area_ch <- shape$area_ch * px_side ^ 2
          shape[6:18] <- apply(shape[6:18], 2, function(x){
            x * px_side
          })
        }
      }


      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        cli::cli_abort("x" = "Only one of {.arg lower_*} or {.arg topn_*} can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$area > lower_size, ],
             shape <- shape[shape$area > mean(shape$area) * lower_noise, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$area, decreasing = TRUE),][1:topn_upper,]
      }
      if(!is.null(lower_eccent)){
        shape <- shape[shape$eccentricity > lower_eccent, ]
      }
      if(!is.null(upper_eccent)){
        shape <- shape[shape$eccentricity < upper_eccent, ]
      }
      if(!is.null(lower_circ)){
        shape <- shape[shape$circularity > lower_circ, ]
      }
      if(!is.null(upper_circ)){
        shape <- shape[shape$circularity < upper_circ, ]
      }
      object_contour <- object_contour[as.character(shape$id)]
      ch <- ch[as.character(shape$id)]
      nmask[!nmask %in% shape$id] <- 0L

      # check if fourier is computed
      if(isTRUE(efourier)){
        efr <-
          efourier(object_contour,
                   nharm = nharm)
        efer <- efourier_error(efr, plot = FALSE)$stats
        efpowwer <- efourier_power(efr, plot = FALSE)
        efpow <- efpowwer$cum_power
        min_harm <- efpowwer$min_harm
        efrn <- efourier_norm(efr)
        efr <- efourier_coefs(efr)
        names(efr)[1] <- "id"
        efrn <- efourier_coefs(efrn)
        names(efrn)[1] <- "id"
      } else{
        efr <- NULL
        efrn <- NULL
        efer <- NULL
        efpow <- NULL
        min_harm <- NULL
      }

      # check if angles should be computed
      if(isTRUE(ab_angles)){
        angles <- poly_apex_base_angle(object_contour, ab_angles_percentiles)
      } else{
        angles <- NULL
      }

      # check if width at should be computed
      if(isTRUE(width_at)){
        widths <-
          do.call(rbind, lapply(object_contour, function(x){
            x |> poly_align(plot = FALSE) |> poly_width_at(width_at_percentiles)
          })) |>
          as.data.frame() |>
          rownames_to_column("id")
        names(widths) <- c("id", paste0("width", width_at_percentiles))
      } else{
        widths <- NULL
      }

      # check if veins is computed
      if(isTRUE(veins)){
        # detect_veins_cpp uses Difference of Gaussians (DoG) bandpass filter
        # + dynamic resolution-independent boundary erosion + per-object Otsu threshold
        # + optional Guo-Hall thinning.
        vein_result <- detect_veins_cpp(
          R_sexp      = img[,,1],
          G_sexp      = img[,,2],
          B_sexp      = img[,,3],
          labels_sexp = as.matrix(image_data(nmask)),
          sigma1     = sigma_veins,
          sigma2     = sigma_veins * 5,
          threshold  = -1,             # -1 = adaptive per-object Otsu
          channel    = 0,              # 0 = green channel (best for leaves)
          erode_size = -1,             # -1 = use rel_erode
          rel_erode  = veins_rel_erode, # 1% of min(W, H)
          thinning   = veins_thinning,  # TRUE for skeleton density, FALSE for area fraction
          return_map = FALSE
        )
        # Build a data.frame indexed by object ID (skip label 0 = background)
        vein_props <- vein_result$proportion
        valid_labs <- which(!is.na(vein_props) & seq_along(vein_props) > 0)
        prop_veins <- data.frame(
          id         = valid_labs,
          prop_veins = vein_props[valid_labs]
        )
        # Keep only objects that appear in shape (same filter as other features)
        prop_veins <- prop_veins[prop_veins$id %in% shape$id, ]
      } else{
        prop_veins <- NULL
      }

      # check if perimeter complexity value is computed
      if(isTRUE(pcv)){
        pcv <- poly_pcv(object_contour, niter = pcv_niter)
      } else{
        pcv <- NULL
      }

      if(!is.null(object_index)){
        object_index_used <- object_index[1]
        if(!is.character(object_index)){
          cli::cli_abort("{.arg object_index} must be a character.")
        }
        ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        if(any(object_index %in% ind$Index)){
          ind_formula <- ind[match(object_index, ind$Index), 2]
        } else{
          ind_formula <- object_index
        }
        ind_name  <- object_index
        data_mask <- as.numeric(image_data(nmask))

        # ── Optimised index pipeline ─────────────────────────────────────────
        # Instead of:  object_rgb → rgb_to_hsb → by() → aggregate()
        # We:
        #   1. Set up full-image evaluation environment (no pixel extraction)
        #   2. Evaluate each formula on the full W×H matrices (vectorised, fast)
        #   3. Aggregate per-object means in C++ (two-pass, zero allocations)
        #
        # This avoids: data.frame creation, by() object splitting, do.call(rbind),
        # aggregate(), and multiple passes over the pixel data.

        R_ch <- img[,,1]; G_ch <- img[,,2]; B_ch <- img[,,3]

        # Build HSB arrays once for the full image (C++ vectorised)
        hsb_arrays <- rgb_to_hsb_cpp(R_ch, G_ch, B_ch)

        # Evaluation environment: named arrays matching formula variable names.
        # rgb_to_hsb() returns lowercase column names: h, s, b (brightness).
        # The original code used these via envir = obj_rgb (a data.frame with
        # columns id, R, G, B, h, s, b).  We must expose the same names here.
        # Note: B = Blue channel (uppercase), b = brightness from HSB (lowercase).
        eval_env <- list(
          R = R_ch, G = G_ch, B = B_ch,     # uppercase: RGB channels [0,1]
          h = hsb_arrays$H,                  # lowercase h, s, b: as in rgb_to_hsb()
          s = hsb_arrays$S,
          b = hsb_arrays$B,
          H = hsb_arrays$H,                  # uppercase aliases (defensive)
          S = hsb_arrays$S
        )

        # Parse formulas once (avoid repeated parse() inside a loop)
        parsed <- lapply(ind_formula, function(f) parse(text = f)[[1]])

        # Evaluate every formula on the full image (one vectorised call each)
        idx_mat <- do.call(cbind, lapply(parsed, function(expr) {
          as.vector(eval(expr, envir = list2env(eval_env, parent = baseenv())))
        }))

        # valid object IDs from shape (already filtered by the analysis)
        valid_ids <- as.integer(shape$id)

        # C++ two-pass aggregation: one sequential scan per formula column
        means_mat <- compute_index_means_cpp(idx_mat, as.vector(data_mask), valid_ids)
        colnames(means_mat) <- ind_name

        indexes <- data.frame(id = valid_ids, means_mat, check.names = FALSE)

        # Pixel-level output only when explicitly requested (large allocation)
        if(isTRUE(pixel_level_index)){
          obj_rgb <- object_rgb(img, data_mask)
          obj_rgb <- subset(obj_rgb, id %in% shape$id)
          obj_rgb <- cbind(obj_rgb, rgb_to_hsb(obj_rgb[, 2:4]))
          # Add per-pixel index values
          pixel_idx <- do.call(cbind, lapply(parsed, function(expr) {
            as.vector(eval(expr, envir = list2env(eval_env, parent = baseenv())))
          }))
          # Subset to pixels belonging to valid objects
          mask_flat   <- as.vector(data_mask)
          pixel_keep  <- mask_flat %in% valid_ids
          obj_rgb <- cbind(obj_rgb,
                           as.data.frame(
                             pixel_idx[pixel_keep, , drop = FALSE],
                             col.names = ind_name))
        } else {
          obj_rgb <- NULL
        }
      } else{
        obj_rgb <- NULL
        indexes <- NULL
        object_index_used <- NULL
      }
      if(isTRUE(return_mask)){
        mask <- nmask
      } else{
        mask <- NULL
      }
      stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area",
                                   "sd_area", "sum_area", "coverage"),
                          value = c(length(shape$area),
                                    min(shape$area),
                                    mean(shape$area),
                                    max(shape$area),
                                    sd(shape$area),
                                    sum(shape$area),
                                    sum(shape$coverage)))
      results <- list(results = shape,
                      statistics = stats,
                      object_rgb = obj_rgb,
                      object_index = indexes,
                      efourier = efr,
                      efourier_norm = efrn,
                      efourier_error = efer,
                      efourier_power = efpow,
                      efourier_minharm = min_harm,
                      veins = prop_veins,
                      angles = angles,
                      width_at = widths,
                      mask = mask,
                      pcv = pcv,
                      contours = object_contour,
                      parms = list(
                        index = index,
                        object_index = object_index_used,
                        reference = isTRUE(reference),
                        reference_area = reference_area,
                        reference_larger = reference_larger,
                        reference_smaller = reference_smaller,
                        npix_ref = if (exists("npix_ref")) npix_ref else NULL,
                        px_side = if (exists("px_side")) px_side else NULL
                      ))
      class(results) <- "anal_obj"
      if(!is.null(object_index)){
        shape_markers <- cbind(shape, indexes)
      } else{
        shape_markers <- shape
      }
      if(plot == TRUE | save_image == TRUE){
        backg <- !is.null(col_background)
        # color for background
        if (is.null(col_background)){
          col_background <- col2rgb("white")
        } else{
          ifelse(is.character(col_background),
                 col_background <- col2rgb(col_background),
                 col_background <- col_background)
        }
        # color for lesions
        if (is.null(col_foreground)){
          col_foreground <- col2rgb("gray")
        } else{
          ifelse(is.character(col_foreground),
                 col_foreground <- col2rgb(col_foreground),
                 col_foreground <- col_foreground)
        }

        if(show_original == TRUE & show_segmentation == FALSE){
          im2 <- img[,,1:3]
          if(backg){
            im3 <- image_color_labels(nmask)
            im2[,,1][which(im3[,,1]==0)] <- col_background[1]
            im2[,,2][which(im3[,,2]==0)] <- col_background[2]
            im2[,,3][which(im3[,,3]==0)] <- col_background[3]
          }
        }
        if(show_original == TRUE & show_segmentation == TRUE){
          im2 <- image_color_labels(nmask)
          if(backg){
            im2[,,1][which(im2[,,1]==0)] <- col_background[1]
            im2[,,2][which(im2[,,2]==0)] <- col_background[2]
            im2[,,3][which(im2[,,3]==0)] <- col_background[3]
          } else{
            im2[,,1][which(im2[,,1]==0)] <- as.numeric(img[,,1])[which(im2[,,1]==0)]
            im2[,,2][which(im2[,,2]==0)] <- as.numeric(img[,,2])[which(im2[,,2]==0)]
            im2[,,3][which(im2[,,3]==0)] <- as.numeric(img[,,3])[which(im2[,,3]==0)]
          }
        }
        if(show_original == FALSE){
          if(show_segmentation == TRUE){
            im2 <- image_color_labels(nmask)
            im2[,,1][which(im2[,,1]==0)] <- col_background[1]
            im2[,,2][which(im2[,,2]==0)] <- col_background[2]
            im2[,,3][which(im2[,,3]==0)] <- col_background[3]
          } else{
            im2 <- img[,,1:3]
            im2[,,1][ID] <- col_foreground[1]
            im2[,,2][ID] <- col_foreground[2]
            im2[,,3][ID] <- col_foreground[3]
            im2[,,1][ID2] <- col_background[1]
            im2[,,2][ID2] <- col_background[2]
            im2[,,3][ID2] <- col_background[3]
          }
        }
        show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
        marker <- ifelse(is.null(marker), "id", marker)
        if (!isFALSE(show_mark) && marker != "point" && !marker %in% colnames(shape_markers)) {
          cli::cli_warn(
            "Accepted {.arg marker} values are: {.val {paste(colnames(shape_markers), collapse = \", \")}}. Drawing the object id instead."
          )
          marker <- "id"
        }

        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
        # correct the contour

        if(plot == TRUE){
          if(marker != "point"){
            plot(im2, max_pixels = max_pixels)
            if(isTRUE(show_contour)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_bbox){
              plot_bbox(object_contour, col = contour_col)
            }
            if(show_mark){
              text(shape_markers[, 2],
                   shape_markers[, 3],
                   round(shape_markers[, marker], 3),
                   col = marker_col,
                   cex = marker_size)
            }
            if(isTRUE(show_chull)){
              plot_contour(ch |> poly_close(), col = "black")
            }
          } else{
            plot(im2, max_pixels = max_pixels)
            if(isTRUE(show_contour)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_bbox){
              plot_bbox(object_contour, col = contour_col)
            }
            if(show_mark){
              points(shape_markers[, 2],
                     shape_markers[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          # plot length and width
          if(isTRUE(show_lw)){
            if(isTRUE(reference)){
              plot_lw(shape_ori)
            } else{
              plot_lw(results)
            }
          }
        }

        if(save_image == TRUE){
          if(dir.exists(diretorio_processada) == FALSE){
            dir.create(diretorio_processada, recursive = TRUE)
          }
          img_d <- dim(image_data(im2))
          tot_pix <- img_d[1] * img_d[2]
          if (!is.null(max_pixels) && tot_pix > max_pixels) {
            scale_factor <- sqrt(max_pixels / tot_pix)
            png_w <- round(img_d[1] * scale_factor)
            png_h <- round(img_d[2] * scale_factor)
          } else {
            png_w <- img_d[1]
            png_h <- img_d[2]
          }
          png(paste0(diretorio_processada, "/",
                     prefix,
                     name_ori, ".",
                     extens_ori),
              width = png_w,
              height = png_h)
          if(marker != "point"){
            plot(im2, max_pixels = max_pixels)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_bbox){
              plot_bbox(object_contour, col = contour_col)
            }
            if(show_mark){
              text(shape_markers[, 2],
                   shape_markers[, 3],
                   round(shape_markers[, marker], 3),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2, max_pixels = max_pixels)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_bbox){
              plot_bbox(object_contour, col = contour_col)
            }
            if(show_mark){
              points(shape_markers[, 2],
                     shape_markers[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          if(isTRUE(show_lw)){
            if(isTRUE(reference)){
              plot_lw(shape_ori)
            } else{
              plot_lw(results)
            }
          }
          # Discrete footer with summary stats
          n_obj   <- nrow(shape)
          a_mean  <- round(mean(shape$area, na.rm = TRUE), 1)
          a_min   <- round(min(shape$area, na.rm = TRUE), 1)
          a_max   <- round(max(shape$area, na.rm = TRUE), 1)
          footer_txt <- sprintf("N: %d  |  Area (px\u00b2)  mean: %s  |  min: %s  |  max: %s",
                                n_obj, a_mean, a_min, a_max)
          usr <- par("usr")
          text(usr[1] + diff(usr[1:2]) * 0.01,
               usr[3] + diff(usr[3:4]) * 0.01,
               footer_txt,
               adj = c(0, 0), cex = 1, col = "#555555", font = 1, xpd = NA)
          dev.off()
        }
      }
      invisible(results)
    }

  if(missing(pattern)){
    if(verbose){
      cli::cli_progress_step(
        msg = "Processing a single image. Please, wait.",
        msg_done = "Image {.emph Successfully} analyzed!",
        msg_failed = "Oops, something went wrong."
      )
    }
    help_count(img, foreground, background, pick_palettes, resize, fill_hull, threshold, erode, dilate, opening, closing, filter, filter_order,
               tolerance, extension, randomize, nrows, plot, show_original,
               show_background, marker, marker_col, marker_size, save_image,
               prefix, dir_original, dir_processed, verbose, col_background,
               col_foreground, lower_noise, ab_angles, ab_angles_percentiles, width_at, width_at_percentiles, return_mask, pcv,
               object_index, max_pixels, return_exact)
  } else{
    if(pattern %in% as.character(0:9)){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    imgpath <- file.path(getwd(), sub('./', '', diretorio_original)) |> trunc_path(max_chars = 50)

    if (length(grep(pattern, names_plant)) == 0) {
      cli::cli_abort("Pattern {.val {pattern}} not found in directory {.path {imgpath}}.")
    }

    allowed_ext <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF")

    if (!all(extensions %in% allowed_ext)) {
      cli::cli_abort("Allowed extensions are {.val {allowed_ext}}.")
    }

    old_opt <- options(cli.progress_bar_style = "bar")
    on.exit(options(old_opt), add = TRUE)


    if (parallel == TRUE) {
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 0.3), workers)
      mirai::daemons(nworkers)
      on.exit(mirai::daemons(0), add = TRUE)
      pkg_root <- tryCatch(rprojroot::find_package_root_file(), error = function(e) getwd())
      is_dev_mode <- file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE)

      if (is_dev_mode) {
        mirai::everywhere(
          {
            .libPaths(lp)
            pkgload::load_all(p_dir, quiet = TRUE, helpers = FALSE)
          },
          lp = .libPaths(),
          p_dir = pkg_root
        )
      } else {
        mirai::everywhere(
          {
            .libPaths(lp)
            library(pliman)
          },
          lp = .libPaths()
        )
      }
      if (verbose) {
        cli::cli_rule(
          left = cli::col_blue("Parallel processing using {nworkers} cores"),
          right = cli::col_blue("Started on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
        )
      }

      # Cria uma função wrapper com todos os argumentos
      # process_image <- function(img) {
      #   help_count(
      #     img = img,
      #     foreground, background, pick_palettes, resize, fill_hull, threshold,
      #     erode, dilate, opening, closing, filter, tolerance, extension,
      #     randomize, nrows, plot, show_original, show_background, marker,
      #     marker_col, marker_size, save_image, prefix, dir_original,
      #     dir_processed, verbose, col_background, col_foreground,
      #     lower_noise, ab_angles, ab_angles_percentiles,
      #     width_at, width_at_percentiles, return_mask, pcv, object_index
      #   )
      # }
      results <- mirai::mirai_map(
        .x = names_plant,
        .f = help_count,
        .args = list(foreground, background, pick_palettes, resize, fill_hull,
                     threshold, erode, dilate, opening, closing, filter, filter_order,
                     tolerance, extension, randomize, nrows, plot, show_original,
                     show_background, marker, marker_col, marker_size, save_image,
                     prefix, dir_original, dir_processed, verbose, col_background,
                     col_foreground, lower_noise, ab_angles, ab_angles_percentiles,
                     width_at, width_at_percentiles, return_mask, pcv,
                     object_index, max_pixels, return_exact)
      )[.progress]
      cli::cli_progress_step(
        msg = "Processing {.val {length(names_plant)}} images found on {.path {imgpath}}. Please, wait.",
        msg_done = "Batch processing finished",
        msg_failed = "Oops, something went wrong."
      )
    } else {
      if (verbose) {
        cli::cli_rule(
          left = cli::col_blue("Analyzing {length(names_plant)} images"),
          right = cli::col_blue("Started at {format(Sys.time(), '%H:%M:%S')}")
        )

        cli::cli_progress_bar(
          format = "{cli::pb_spin} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
          total = length(names_plant),
          clear = FALSE
        )
      }
      results <- vector("list", length(names_plant))
      for (i in seq_along(names_plant)) {
        if (verbose) cli::cli_progress_update()
        results[[i]] <-
          help_count(
            img = names_plant[i],
            foreground, background, pick_palettes, resize, fill_hull, threshold,
            erode, dilate, opening, closing, filter, filter_order,
            tolerance, extension, randomize, nrows, plot, show_original,
            show_background, marker, marker_col, marker_size, save_image,
            prefix, dir_original, dir_processed, verbose, col_background,
            col_foreground, lower_noise, ab_angles, ab_angles_percentiles,
            width_at, width_at_percentiles, return_mask, pcv,
            object_index, max_pixels, return_exact
          )
      }
      if (verbose) {
        cli::cli_progress_done()
      }
    }

    ## bind the results

    names(results) <- names_plant
    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["statistics"]],
                          img =  names(results[i]))[,c(3, 1, 2)]
              })
      )

    if(!is.null(object_index)){
      if(!is.null(results[[1]][["object_rgb"]])){
        obj_rgb <-
          do.call(rbind,
                  lapply(seq_along(results), function(i){
                    transform(results[[i]][["object_rgb"]],
                              img =  names(results[i]))
                  })
          )
        obj_rgb <- obj_rgb[, c(ncol(obj_rgb), 1:ncol(obj_rgb) - 1)]
      } else{
        obj_rgb <- NULL
      }
      object_index <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["object_index"]],
                            img =  names(results[i]))
                })
        )
      object_index <- object_index[, c(ncol(object_index), 1:ncol(object_index) - 1)]
    } else{
      obj_rgb <- NULL
      object_index <- NULL
    }

    if(!isFALSE(efourier)){
      efourier <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier"]],
                            img =  names(results[i]))
                })
        )
      efourier <- efourier[, c(ncol(efourier), 1:ncol(efourier)-1)]
      names(efourier)[2] <- "id"

      efourier_norm <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_norm"]],
                            img =  names(results[i]))
                })
        )
      efourier_norm <- efourier_norm[, c(ncol(efourier_norm), 1:ncol(efourier_norm)-1)]
      names(efourier_norm)[2] <- "id"


      efourier_error <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_error"]],
                            img =  names(results[i]))
                })
        )
      efourier_error <- efourier_error[, c(ncol(efourier_error), 1:ncol(efourier_error)-1)]
      names(efourier_error)[2] <- "id"

      efourier_power <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_power"]],
                            img =  names(results[i]))
                })
        )
      efourier_power <- efourier_power[, c(ncol(efourier_power), 1:ncol(efourier_power)-1)]
      names(efourier_power)[2] <- "id"

      efourier_minharm <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_minharm"]],
                            img =  names(results[i]))
                })
        )
      efourier_minharm <- efourier_minharm[, c(ncol(efourier_minharm), 1:ncol(efourier_minharm)-1)]
      names(efourier_minharm)[2] <- "id"

    } else{
      efourier <- NULL
      efourier_norm <- NULL
      efourier_error <- NULL
      efourier_power <- NULL
      efourier_minharm <- NULL
    }



    if(isTRUE(veins)){
      veins <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["veins"]],
                            img =  names(results[i]))
                })
        )

      veins <- veins[, c(ncol(veins), 1:ncol(veins)-1)]
    } else{
      veins <- NULL
    }



    if(isTRUE(ab_angles)){
      angles <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["angles"]],
                            img =  names(results[i]))
                })
        )

      angles <- angles[, c(ncol(angles), 1:ncol(angles)-1)]
    } else{
      angles <- NULL
    }

    if(isTRUE(width_at)){
      width_at <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["width_at"]],
                            img =  names(results[i]))
                })
        )

      width_at <- width_at[, c(ncol(width_at), 1:ncol(width_at)-1)]
    } else{
      width_at <- NULL
    }

    if(isTRUE(pcv)){
      pcv <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  data.frame(pcv = results[[i]][["pcv"]]) |>
                    transform(img =  names(results[i]))
                })
        )

      pcv <- pcv[, c("img", "pcv")]
    } else{
      pcv <- NULL
    }

    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["results"]],
                          img =  names(results[i]))
              })
      )

    if("img" %in% colnames(results)){
      results <- results[, c(ncol(results), 1:ncol(results) - 1)]
    }
    nimages <- length(unique(stats$id))
    n_img <-
      results |>
      dplyr::group_by(img) |>
      dplyr::summarise(
        n = dplyr::n(),
        area_mean = mean(area, na.rm = TRUE),
        area_min = min(area, na.rm = TRUE),
        area_max = max(area, na.rm = TRUE),
        area_sum = sum(area, na.rm = TRUE),
        area_sd = sd(area, na.rm = TRUE)
      )


    if(verbose == TRUE){
      average_n <- mean(n_img$n)
      min_n <- min(n_img$n)
      max_n <- max(n_img$n)
      average_area <- mean(n_img$area_mean)
      min_area <- min(n_img$area_max)
      max_area <- max(n_img$area_min)

      # Global statistics
      glob_stat <- cli::ansi_columns(
        paste(
          c(
            "Total objects:",
            "Total area:",
            "Overall mean area:",
            "Overall SD:",
            "Min area:",
            "Max area:"
          ),
          c(
            sum(n_img$n),
            round(sum(n_img$area_sum, na.rm = TRUE), 2),
            round(mean(results$area, na.rm = TRUE), 2),
            round(sd(results$area, na.rm = TRUE), 2),
            round(min(results$area, na.rm = TRUE), 2),
            round(max(results$area, na.rm = TRUE), 2)
          )
        ),
        width = 60,
        fill = "rows",
        align = "left",
        sep = "",
        max_cols = 2
      )
      cli::boxx(glob_stat, header = "Global statistics ")  |> cat(sep = "\n")

      cross_imgstat <-
        cli::ansi_columns(
          paste(
            c(
              "Avg objects:",
              "Avg sum area:",
              "Min objects:",
              "Max objects:",
              "Avg area:",
              "Avg SD area:",
              "Min mean area:",
              "Max mean area:"
            ),
            c(
              round(mean(n_img$n), 2),
              round(mean(n_img$area_sum, na.rm = TRUE), 2),
              min(n_img$n),
              max(n_img$n),
              round(mean(n_img$area_mean, na.rm = TRUE), 2),
              round(mean(n_img$area_sd, na.rm = TRUE), 2),
              round(min(n_img$area_mean, na.rm = TRUE), 2),
              round(max(n_img$area_mean, na.rm = TRUE), 2)
            )
          ),
          width = 60,
          fill = "rows",
          align = "left",
          sep = "",
          max_cols = 2
        )

      cli::boxx(cross_imgstat,
                header = "Across-image statistics (per-image averages)",
                footer = paste0("Based on ", nimages, " images")) |>
        cat(sep = "\n")

      cli::cli_rule(
        left = cli::col_blue("Processing successfully finished"),
        right = cli::col_blue("on {format(Sys.time(), format = '%Y-%m-%d | %H:%M:%OS0')}")
      )
    }

    invisible(
      structure(
        list(statistics = n_img,
             count = stats[stats$stat == "n", c(1, 3)],
             results = results,
             obj_rgb = obj_rgb,
             object_index = object_index,
             efourier = efourier,
             efourier_norm = efourier_norm,
             efourier_error = efourier_error,
             efourier_minharm = efourier_minharm,
             veins = veins,
             angles = angles,
             width_at = width_at,
             pcv = pcv),
        class = "anal_obj_ls"
      )
    )
  }
}


#' @name analyze_objects
#' @param x An object of class `anal_obj`.
#' @param which Which to plot. Either 'measure' (object measures) or 'index'
#'   (object index). Defaults to `"measure"`.
#' @param measure The measure to plot. Defaults to `"area"`.
#' @param type The type of plot. Either `"hist"` or `"density"`. Partial matches
#'   are recognized.
#' @method plot anal_obj
#' @export
#'
#' @examples
#' if (interactive()) {
#' library(pliman)
#'
#' img <- image_pliman("soy_green.jpg")
#' # Segment the foreground (grains) using the normalized blue index (NB, default)
#' # Shows the average value of the blue index in each object
#'
#' rgb <-
#'    analyze_objects(img,
#'                    marker = "id",
#'                    object_index = "B",
#'                    pixel_level_index = TRUE)
#' # density of area
#' plot(rgb)
#'
#' # histogram of perimeter
#' plot(rgb, measure = "perimeter", type = "histogram") # or 'hist'
#'
#' # density of the blue (B) index
#' plot(rgb, which = "index")
#' }
plot.anal_obj <- function(x,
                          which = "measure",
                          measure = "area",
                          type = c("density", "histogram"),
                          ...){
  if(!which %in% c("measure", "index")){
    cli::cli_abort("{.arg which} must be one of {.val measure} or {.val index}.")
  }
  if(which == "measure"){
    nam <- colnames(x$results)
    if(!measure %in% nam){
      cli::cli_abort(c(
        "x" = "Measure {.val {measure}} not available in {.arg x}.",
        "i" = "Try one of {.val {paste(nam, collapse = \", \")}}."
      ))

    }
    temp <- x$results[[measure]]
    types <- c("density", "histogram")
    matches <- grepl(type[1], types)
    type <- types[matches]
    if(type == "histogram"){
      hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
    } else{
      density_data <- density(temp)  # Calculate the density for the column
      plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
      points(x = temp, y = rep(0, length(temp)), col = "red")
    }
  } else{
    rgb <- x$object_rgb
    if(is.null(rgb)){
      cli::cli_abort(c(
        "x" = "RGB values not found. Use {.arg object_index} in {.fn analyze_objects}().",
        "i" = "Have you accidentally missed the {.arg pixel_level_index} = TRUE argument?"
      ))
    }
    plot(density(rgb$R),
         main = NA,
         col = "red",
         lwd = 2,
         xlim = c(min(rgb$R, rgb$G, rgb$B), max(rgb$R, rgb$G, rgb$B)),
         ylim = c(0, max(density(rgb$R)$y, density(rgb$G)$y, density(rgb$B)$y)),
         xlab = "Pixel value",
         ylab = "Density")

    # Add the density curves for G and B
    lines(density(rgb$G), col = "green", lwd = 2)
    lines(density(rgb$B), col = "blue", lwd = 2)
    # Add a legend
    legend("topright", legend = c("R", "G", "B"), col = c("red", "green", "blue"), lty = 1,
           lwd = 2)

  }
}
#' @name analyze_objects
#' @export
plot.anal_obj_ls <- function(x,
                             which = "measure",
                             measure = "area",
                             type = c("density", "histogram"),
                             ...){
  if(!which %in% c("measure", "index")){
    cli::cli_abort("{.arg which} must be one of {.val measure} or {.val index}.")
  }
  if(which == "measure"){
    nam <- colnames(x$results)
    if (!measure %in% nam) {
      cli::cli_abort(c(
        "x" = "Measure {.val {measure}} not available in {.arg x}.",
        "i" = "Try one of {.val {paste(nam, collapse = \", \")}}."
      ))
    }

    temp <- x$results[[measure]]
    types <- c("density", "histogram")
    matches <- grepl(type[1], types)
    type <- types[matches]
    if(type == "histogram"){
      hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
    } else{
      density_data <- density(temp)  # Calculate the density for the column
      plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
      points(x = temp, y = rep(0, length(temp)), col = "red")
    }
  } else{
    rgb <- x$object_rgb
    if (is.null(rgb)) {
      cli::cli_abort(c(
        "x" = "RGB values not found. Use {.arg object_index} in {.fn analyze_objects}().",
        "i" = "Have you accidentally missed the {.arg pixel_level_index} = TRUE argument?"
      ))
    }

    plot(density(rgb$R),
         main = NA,
         col = "red",
         lwd = 2,
         xlim = c(min(rgb$R, rgb$G, rgb$B), max(rgb$R, rgb$G, rgb$B)),
         ylim = c(0, max(density(rgb$R)$y, density(rgb$G)$y, density(rgb$B)$y)),
         xlab = "Pixel value",
         ylab = "Density")

    # Add the density curves for G and B
    lines(density(rgb$G), col = "green", lwd = 2)
    lines(density(rgb$B), col = "blue", lwd = 2)
    # Add a legend
    legend("topright", legend = c("R", "G", "B"), col = c("red", "green", "blue"), lty = 1,
           lwd = 2)

  }
}



#' @export
#' @name analyze_objects
analyze_objects_iter <- function(pattern,
                                 known_area,
                                 verbose = TRUE,
                                 ...){
  if (interactive()) {
    imgs <- list.files(pattern = pattern)
    measures <- list()
    for (i in 1:length(imgs)) {
      tmp <-
        analyze_objects(img = file_name(imgs[[i]]),
                        marker = "id",
                        ...)
      object <- as.numeric(readline("Known object: "))
      tmp <- get_measures(tmp,
                          id = object,
                          measure = area ~ known_area,
                          verbose = verbose)
      tmp$img <- file_name(imgs[[i]])
      tmp <- tmp[, c(ncol(tmp), 1:ncol(tmp) - 1)]
      measures[[i]] <- tmp
    }
    do.call(rbind, lapply(measures, function(x){x}))
  }
}


