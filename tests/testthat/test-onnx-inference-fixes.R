# Regression tests for the deep-learning inference engine.
#
# These cover defects that produced *wrong results* rather than merely slow ones, so
# they are the ones most worth guarding. Each test needs the ONNX Runtime library and
# the relevant model weights; all of them skip cleanly when those are absent, so the
# file is safe to run on CRAN and on machines that have never called onnx_setup().

skip_if_no_onnx <- function() {
  lib <- tryCatch(pliman_onnx_lib_path(engine = "cpu"), error = function(e) NULL)
  skip_if(is.null(lib) || !file.exists(lib), "ONNX Runtime library not installed")
  lib
}

skip_if_no_model <- function(...) {
  files <- file.path(pliman_model_dir(), c(...))
  skip_if(!all(file.exists(files)),
          paste("model weights not cached:", paste(basename(files), collapse = ", ")))
  normalizePath(files, winslash = "/", mustWork = FALSE)
}

sam_tensor <- function(mat) {
  as.numeric(pliman:::.preprocess_nchw(
    mat, target_size = 1024L,
    mean = c(0.485, 0.456, 0.406),
    std  = c(0.229, 0.224, 0.225)
  ))
}

test_that("a per-call threads argument does not change the global OpenMP default", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("sam2.1.encoder.onnx", "sam2.1.decoder.onnx")

  mat <- image_data(image_pliman("soybean_touch.jpg"))
  tv <- sam_tensor(mat)
  w <- dim(mat)[1]
  h <- dim(mat)[2]

  before <- pliman_omp_max_threads_cpp()
  invisible(run_sam2_persam_cpp(
    tensor_vec = tv,
    exemplar_x = w / 2, exemplar_y = h / 2,
    orig_w = as.double(w), orig_h = as.double(h),
    encoder_path = paths[1], decoder_path = paths[2], lib_path = lib,
    max_objects = 5L, feat_res = 64L, num_threads = 2L,
    use_gpu = FALSE, device_id = -1L
  ))
  expect_identical(pliman_omp_max_threads_cpp(), before)
})

test_that("the SAM encoder cache distinguishes tensors that differ slightly", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("sam2.1.encoder.onnx", "sam2.1.decoder.onnx")

  mat <- image_data(image_pliman("soybean_touch.jpg"))
  w <- dim(mat)[1]
  h <- dim(mat)[2]
  ten_a <- sam_tensor(mat)

  # Differs from ten_a only inside a small patch. A hash over sampled elements, as an
  # earlier implementation used, would very likely miss this and serve ten_a's
  # embeddings for ten_b: two different frames of a fixed-camera scene colliding.
  ten_b <- ten_a
  idx <- seq.int(700000, 700400)
  ten_b[idx] <- ten_b[idx] + 1.5

  persam <- function(tv) {
    run_sam2_persam_cpp(
      tensor_vec = tv,
      exemplar_x = w / 2, exemplar_y = h / 2,
      orig_w = as.double(w), orig_h = as.double(h),
      encoder_path = paths[1], decoder_path = paths[2], lib_path = lib,
      max_objects = 10L, feat_res = 64L, num_threads = 0L,
      use_gpu = FALSE, device_id = -1L
    )$similarity_map
  }

  onnx_clear_encoder_cache_cpp()
  map_a <- persam(ten_a)
  map_a_again <- persam(ten_a)
  map_b <- persam(ten_b)

  # Identical input must give an identical map (the cache is a cache, not a guess)...
  expect_equal(map_a, map_a_again)
  # ...and a changed input must not be served from it.
  expect_false(isTRUE(all.equal(map_a, map_b)))
})

test_that("PerSAM returns dense kernels that survive a prototype round trip", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("sam2.1.encoder.onnx", "sam2.1.decoder.onnx")

  mat <- image_data(image_pliman("soybean_touch.jpg"))
  w <- dim(mat)[1]
  h <- dim(mat)[2]
  tv <- sam_tensor(mat)

  first <- run_sam2_persam_cpp(
    tensor_vec = tv,
    exemplar_x = w / 2, exemplar_y = h / 2,
    orig_w = as.double(w), orig_h = as.double(h),
    encoder_path = paths[1], decoder_path = paths[2], lib_path = lib,
    max_objects = 10L, feat_res = 256L, num_threads = 0L,
    use_gpu = FALSE, device_id = -1L
  )

  proto <- first$prototypes
  expect_true(!is.null(proto$dense))
  expect_gte(length(proto$dense), 1L)
  expect_true(all(c("dx", "dy", "w", "t_zero", "radius", "sum_w", "tmpl_std") %in%
                    names(proto$dense[[1]])))

  # With carried prototypes and no clicks, feat_res = 1024 used to scan a million
  # pixels with no templates and return an all-zero map, hence zero instances.
  carried <- run_sam2_persam_cpp(
    tensor_vec = tv,
    exemplar_x = numeric(0), exemplar_y = numeric(0),
    orig_w = as.double(w), orig_h = as.double(h),
    encoder_path = paths[1], decoder_path = paths[2], lib_path = lib,
    sim_threshold = 0.5, max_objects = 50L, feat_res = 1024L, num_threads = 0L,
    use_gpu = FALSE, device_id = -1L,
    precomputed_prototypes = proto
  )
  expect_gt(length(carried$scores), 0L)
})

test_that("a SAM box prompt segments inside the box", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("sam2.1.encoder.onnx", "sam2.1.decoder.onnx")

  mat <- image_data(image_pliman("soybean_touch.jpg"))
  tv <- sam_tensor(mat)

  box <- c(0.30, 0.30, 0.70, 0.70)
  m <- run_sam2_inference_cpp(
    tensor_vec = tv,
    points_x = c(box[1], box[3]) * 1024,
    points_y = c(box[2], box[4]) * 1024,
    point_labels = c(2L, 3L),
    encoder_path = paths[1], decoder_path = paths[2], lib_path = lib,
    num_threads = 0L, use_gpu = FALSE, device_id = -1L
  )

  inside <- m[round(box[1] * 256):round(box[3] * 256),
              round(box[2] * 256):round(box[4] * 256)]
  # The corner points define one box prompt. Sent as two independent single-point
  # prompts with labels 2 and 3, as an earlier version did, the mask bore no relation
  # to the requested region.
  expect_gt(mean(inside > 0.5), mean(m > 0.5))
})

test_that("Grounding DINO honours text_threshold", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("groundingdino-tiny.onnx", "vocab.txt")

  mat <- image_data(image_pliman("soybean_touch.jpg"))
  pv <- as.numeric(pliman:::.preprocess_nchw(
    mat, target_size = 800L,
    mean = c(0.485, 0.456, 0.406), std = c(0.229, 0.224, 0.225)
  ))
  tok <- pliman:::.bert_tokenize("seed . leaf . soil", paths[2])

  n_boxes <- function(tt) {
    length(run_grounding_dino_cpp(
      pixel_values = pv,
      input_ids = as.integer(tok$input_ids),
      token_type_ids = as.integer(tok$token_type_ids),
      attention_mask = as.integer(tok$attention_mask),
      model_path = paths[1], lib_path = lib,
      box_threshold = 0.15, text_threshold = tt, iou_threshold = 0.5,
      num_threads = 0L, use_gpu = FALSE, device_id = -1L
    )$scores)
  }

  # text_threshold was accepted and never used, so it had no effect whatsoever.
  expect_lte(n_boxes(0.95), n_boxes(0.05))
})

test_that("one model can be evicted without flushing the session cache", {
  lib <- skip_if_no_onnx()
  paths <- skip_if_no_model("yolo26n.onnx", "u2netp.onnx")

  img <- image_pliman("soybean_touch.jpg")
  invisible(image_detect_dl(img, model = "yolo26n", engine = "cpu",
                            verbose = FALSE, plot = FALSE))
  invisible(image_segment_dl(img, model = "u2netp", type = "mask", engine = "cpu",
                             verbose = FALSE, plot = FALSE))

  before <- nrow(pliman_session_info())
  skip_if(before < 2, "fewer than two sessions cached")

  evicted <- onnx_evict_session_cpp(paths[1])
  after <- nrow(pliman_session_info())

  expect_gte(evicted, 1L)
  expect_equal(after, before - evicted)
  expect_gt(after, 0L)   # the other model is still warm
})
