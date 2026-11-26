# Regression tests for the mapping between instances, their contours and their traits.
#
# `extract_contours_cpp()` returns one slot per label id, so slot k belongs to instance
# k and is NULL when that instance kept no pixels. An instance can end up empty even
# though it was detected: overlapping pixels keep the best-scoring label, so a weaker
# instance fully covered by a stronger one disappears from the label matrix, and
# morphological cleaning can remove a small mask entirely.
#
# Two things used to go wrong as a result, and both are silent or fatal rather than
# obvious:
#   * the list came back shorter than the number of boxes, and indexing it per box
#     overran ("subscript out of bounds");
#   * dropping the empty slots before measuring renumbered every later instance, so
#     per-instance traits were joined onto the wrong objects.

make_labels <- function(w = 200, h = 200, empty_id = 3L) {
  lbl <- matrix(0L, nrow = w, ncol = h)
  boxes <- list(c(10, 10, 40, 40), c(60, 10, 90, 40), c(110, 10, 140, 40),
                c(10, 60, 40, 90), c(60, 60, 90, 90))
  for (k in seq_along(boxes)) {
    if (identical(as.integer(k), as.integer(empty_id))) next
    b <- boxes[[k]]
    lbl[b[1]:b[3], b[2]:b[4]] <- k
  }
  lbl
}

test_that("contour slots stay aligned with instance ids when one instance is empty", {
  lbl <- make_labels(empty_id = 3L)
  conts <- extract_contours_cpp(lbl)

  # Slot 3 is empty but slots 4 and 5 still exist, so the list keeps its full length.
  expect_equal(length(conts), 5L)
  expect_true(is.null(conts[[3]]))
  expect_true(is.matrix(conts[[4]]))

  # Measuring the full list recovers the true ids; measuring a compacted list does not.
  feats <- poly_measures(conts)
  expect_equal(as.integer(feats$id), c(1L, 2L, 4L, 5L))

  compacted <- poly_measures(conts[!vapply(conts, is.null, logical(1))])
  expect_false(identical(as.integer(compacted$id), c(1L, 2L, 4L, 5L)))

  # Joining traits by id must leave the empty instance without traits rather than
  # borrowing a neighbour's.
  area_by_id <- feats$area[match(1:5, feats$id)]
  expect_true(is.na(area_by_id[3]))
  expect_false(any(is.na(area_by_id[c(1, 2, 4, 5)])))
})

test_that("trailing empty instances do not shorten the contour list past the ids", {
  # The last instance is the empty one, so max(label) < number of instances and the
  # returned list is shorter than the box table.
  lbl <- make_labels(empty_id = 5L)
  conts <- extract_contours_cpp(lbl)
  expect_equal(length(conts), 4L)

  # This is the shape that used to crash: padding to the box count keeps per-box
  # indexing in range while preserving the alignment of the slots that do exist.
  num_inst <- 5L
  length(conts) <- num_inst
  expect_equal(length(conts), num_inst)
  expect_true(is.null(conts[[5]]))
  expect_silent(for (k in seq_len(num_inst)) invisible(conts[[k]]))
})

test_that("instance segmentation survives a detection whose mask is fully overlapped", {
  lib <- tryCatch(pliman_onnx_lib_path(engine = "cpu"), error = function(e) NULL)
  skip_if(is.null(lib) || !file.exists(lib), "ONNX Runtime library not installed")
  skip_if(!file.exists(file.path(pliman_model_dir(), "yolo26n-seg.onnx")),
          "yolo26n-seg weights not cached")

  img <- image_pliman("soybean_touch.jpg")

  # A low confidence threshold admits weak, heavily overlapping detections, which is
  # what produces instances with no surviving pixels.
  res <- image_segment_dl(img, model = "yolo26n-seg", type = "highlight",
                          conf_threshold = 0.05, iou_threshold = 0.9,
                          engine = "cpu", verbose = FALSE, plot = FALSE)

  expect_s3_class(res, "image_segment_dl")
  expect_true(is.data.frame(res$results))
  # Every box keeps its own row, whether or not a contour was recoverable for it.
  expect_equal(nrow(res$results), nrow(res$boxes))
  expect_true(all(c("area", "perimeter") %in% names(res$results)))
})
