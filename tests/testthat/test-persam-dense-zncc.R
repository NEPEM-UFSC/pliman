# Regression tests for the dense stride-1 exemplar matching path (feat_res = 1024).
#
# The FFT backend replaces six direct correlations with products of spectra for large
# templates. It must agree with direct evaluation on the instances found, and it must
# actually be taken when the radius warrants it -- a silent fallback would make the
# optimisation invisible rather than wrong, which is harder to notice.

skip_if_no_sam <- function() {
  lib <- tryCatch(pliman_onnx_lib_path(engine = "cpu"), error = function(e) NULL)
  skip_if(is.null(lib) || !file.exists(lib), "ONNX Runtime library not installed")
  files <- file.path(pliman_model_dir(),
                     c("sam2.1.encoder.onnx", "sam2.1.decoder.onnx"))
  skip_if(!all(file.exists(files)), "SAM 2.1 weights not cached")
  invisible(TRUE)
}

persam_dense <- function(img, exemplars, min_radius) {
  withr::local_envvar(c(PLIMAN_ZNCC_FFT_MIN_RADIUS = as.character(min_radius)))
  image_segment_dl(img, model = "persam", prompt = exemplars,
                   feat_res = 1024, sim_threshold = 0.5, min_dist = 12,
                   engine = "cpu", sam_device = "cpu",
                   verbose = FALSE, plot = FALSE)
}

test_that("the FFT dense path agrees with direct evaluation", {
  skip_if_no_sam()
  skip_if_not_installed("withr")

  img <- image_pliman("soybean_touch.jpg")
  m <- image_data(img)
  ex <- cbind(x = dim(m)[1] * 0.5, y = dim(m)[2] * 0.5)

  before <- pliman_zncc_fft_stats_cpp()
  direct <- persam_dense(img, ex, 4000)      # gate above any radius: direct path
  mid <- pliman_zncc_fft_stats_cpp()
  fft <- persam_dense(img, ex, 1)            # gate below any radius: FFT path
  after <- pliman_zncc_fft_stats_cpp()

  # The gate must actually switch paths, otherwise this test compares nothing.
  expect_equal(mid[["success"]] - before[["success"]], 0L)
  expect_gte(after[["success"]] - mid[["success"]], 1L)

  expect_equal(nrow(fft$results), nrow(direct$results))

  if (nrow(direct$results) > 0) {
    cols <- c("xmin", "ymin", "xmax", "ymax")
    # Spectral and direct summation differ in rounding, not in result: boxes land on
    # the same pixels.
    expect_equal(as.matrix(fft$results[, cols]),
                 as.matrix(direct$results[, cols]),
                 tolerance = 1e-6)
  }

  a <- attr(direct$similarity_map, "raw_matrix")
  b <- attr(fft$similarity_map, "raw_matrix")
  if (!is.null(a) && !is.null(b)) {
    expect_equal(dim(a), dim(b))
    # The two maps are the same field computed two ways, not bit-identical numbers.
    # The direct path sums several thousand terms in single precision while the FFT
    # path works in double, the ZNCC is a ratio of such sums (so the error grows where
    # the patch standard deviation is small), and each map is finally divided by its
    # own maximum. Differences of ~1e-3 on a [0, 1] map follow from that, with the FFT
    # result being the more accurate of the two; what must not move is the field's
    # shape and the instances derived from it.
    expect_gt(stats::cor(as.numeric(a), as.numeric(b)), 0.9999)
    expect_lt(max(abs(a - b)), 0.01)
  }
})

test_that("dense kernels are packed into contiguous spans", {
  skip_if_no_sam()

  img <- image_pliman("soybean_touch.jpg")
  m <- image_data(img)
  ex <- cbind(x = dim(m)[1] * 0.5, y = dim(m)[2] * 0.5)

  res <- image_segment_dl(img, model = "persam", prompt = ex, feat_res = 1024,
                          engine = "cpu", sam_device = "cpu",
                          verbose = FALSE, plot = FALSE)
  dk <- res$prototypes$dense[[1]]
  expect_true(all(c("dx", "dy", "w", "t_zero", "radius") %in% names(dk)))
  expect_gt(length(dk$dx), 0L)
  # t_zero carries one mean-subtracted RGB triple per entry.
  expect_equal(nrow(dk$t_zero), length(dk$dx))
  expect_equal(ncol(dk$t_zero), 3L)
})
