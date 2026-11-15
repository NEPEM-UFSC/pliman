# Regression tests for the execution-provider layer, the optimized-graph cache, the
# YOLO inference prefetch / tile batching paths and the native batch renderer
# (OPTIMIZATION_PLAN phases 3-5). The speed paths must never change results, so most
# tests compare a fast path against the plain one. Everything skips when ONNX Runtime or
# the weights are not available.

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

# 640 x 360 interleaved RGB frame (non-square on purpose: W x H convention)
test_frame <- function() {
  img <- image_resize(image_pliman("soybean_touch.jpg"), width = 640)
  img <- image_crop(img, width = 1:640, height = 1:360)
  image_to_rgb_buffer_cpp(image_data(img))
}

test_that("pliman_gpu_provider() validates, sets and restores the provider", {
  skip_if_not_installed("withr")
  withr::local_envvar(c(PLIMAN_GPU_PROVIDER = NA, PLIMAN_TRT_FP16 = NA))
  expect_equal(pliman_gpu_provider(), "auto")
  prev <- pliman_gpu_provider("directml")
  expect_equal(prev, "auto")
  expect_equal(pliman_gpu_provider(), "dml")
  pliman_gpu_provider(trt_fp16 = TRUE)
  expect_equal(Sys.getenv("PLIMAN_TRT_FP16"), "1")
  expect_error(pliman_gpu_provider("not-a-provider"))
  pliman_gpu_provider("auto")
  expect_equal(pliman_gpu_provider(), "auto")
})

test_that("CPU sessions report the cpu provider and the runtime lists it", {
  lib <- skip_if_no_onnx()
  m <- skip_if_no_model("yolo26n.onnx")
  expect_true("CPUExecutionProvider" %in% onnx_available_providers_cpp(lib))
  invisible(run_yolo_cpp(test_frame(), 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE))
  info <- pliman_session_info()
  row <- info[info$model == "yolo26n.onnx" & !info$gpu, , drop = FALSE]
  expect_gte(nrow(row), 1L)
  expect_true(all(row$provider == "cpu"))
  expect_true(all(row$create_ms >= 0))
})

test_that("a GPU session that fails to build warns, falls back and is retried", {
  skip_if_not_installed("withr")
  lib <- tryCatch(pliman_onnx_lib_path(engine = "gpu"), error = function(e) NULL)
  skip_if(is.null(lib) || !file.exists(lib), "GPU ONNX Runtime library not installed")
  skip_if(!any(c("DmlExecutionProvider", "CUDAExecutionProvider") %in% onnx_available_providers_cpp(lib)),
          "no GPU execution provider")
  m <- skip_if_no_model("yolo26n.onnx")
  fr <- test_frame()
  pliman_clear_sessions()
  withr::local_envvar(c(PLIMAN_TEST_FAIL_GPU_SESSION = "1", PLIMAN_GPU_RETRY_S = "0"))
  expect_warning(
    r_cpu <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = TRUE),
    "running on the CPU"
  )
  info <- pliman_session_info()
  expect_false(info$gpu[1])
  expect_true(info$fell_back_to_cpu[1])
  expect_match(info$fallback_reason[1], "simulated GPU session failure")

  # once the GPU is available again the next call rebuilds the session on it
  Sys.unsetenv("PLIMAN_TEST_FAIL_GPU_SESSION")
  r_gpu <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = TRUE)
  info <- pliman_session_info()
  expect_true(info$gpu[1])
  expect_equal(info$fallback_reason[1], "")
  expect_equal(nrow(r_gpu$boxes), nrow(r_cpu$boxes))
  pliman_clear_sessions()
})

test_that("the optimized-graph cache is reused and does not change detections", {
  skip_if_not_installed("withr")
  lib <- skip_if_no_onnx()
  m <- skip_if_no_model("yolo26n.onnx")
  cache <- withr::local_tempdir()
  withr::local_envvar(c(PLIMAN_ORT_OPT_CACHE = NA))
  fr <- test_frame()

  pliman_clear_sessions()
  ref <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)

  pliman_onnx_graph_cache(TRUE, dir = cache)
  pliman_clear_sessions()
  first <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)
  expect_equal(pliman_session_info()$graph_cache[1], "saved")
  expect_length(list.files(cache, "\\.opt\\.onnx$"), 1L)

  pliman_clear_sessions()
  hit <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)
  expect_equal(pliman_session_info()$graph_cache[1], "hit")
  pliman_onnx_graph_cache(FALSE)
  pliman_clear_sessions()

  expect_equal(first$boxes, ref$boxes, tolerance = 1e-5)
  expect_equal(hit$boxes, ref$boxes, tolerance = 1e-5)
  expect_equal(hit$class_ids, ref$class_ids)
})

test_that("a prefetched frame gives exactly the synchronous result", {
  lib <- skip_if_no_onnx()
  m <- skip_if_no_model("yolo26n.onnx")
  fr <- test_frame()
  ref <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)

  s0 <- yolo_prefetch_stats_cpp()
  expect_true(yolo_prefetch_cpp(fr, model_path = m, lib_path = lib, use_gpu = FALSE))
  pf <- run_yolo_cpp(fr, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)
  s1 <- yolo_prefetch_stats_cpp()
  expect_equal(s1[["used"]] - s0[["used"]], 1L)
  expect_identical(pf$boxes, ref$boxes)
  expect_identical(pf$scores, ref$scores)
  expect_identical(pf$class_ids, ref$class_ids)

  # a different frame must not consume the prefetched outputs
  other <- fr
  other[1:3000] <- as.raw(0)
  expect_true(yolo_prefetch_cpp(fr, model_path = m, lib_path = lib, use_gpu = FALSE))
  r_other <- run_yolo_cpp(other, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)
  s2 <- yolo_prefetch_stats_cpp()
  expect_equal(s2[["used"]] - s1[["used"]], 0L)
  expect_equal(s2[["discarded"]] - s1[["discarded"]], 1L)
  ref_other <- run_yolo_cpp(other, 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE)
  expect_identical(r_other$boxes, ref_other$boxes)
  yolo_prefetch_cancel_cpp()
})

test_that("batched tiles match per-tile inference with a batch-dynamic model", {
  skip_if_not_installed("withr")
  lib <- skip_if_no_onnx()
  dyn <- Sys.getenv("PLIMAN_TEST_DYNAMIC_YOLO", "")
  skip_if(!nzchar(dyn) || !file.exists(dyn), "set PLIMAN_TEST_DYNAMIC_YOLO to a batch-dynamic YOLO export")
  img <- image_resize(image_pliman("soybean_touch.jpg"), width = 1600)
  bm <- image_to_rgb_buffer_cpp(image_data(img))
  W <- dim(bm)[2]; H <- dim(bm)[3]
  gx <- rep(c(0L, 480L, W - 640L), times = 2); gy <- rep(c(0L, H - 640L), each = 3)
  run <- function(b) {
    withr::local_envvar(c(PLIMAN_YOLO_TILE_BATCH = as.character(b)))
    run_yolo_tiles_cpp(bm, gx, gy, 640L, 640L, 0.25, 0.45, dyn, lib, 0L, FALSE, -1L, 640L, FALSE)
  }
  r1 <- run(1); r4 <- run(4)
  expect_equal(nrow(r4$boxes), nrow(r1$boxes))
  o1 <- order(r1$tile_index, r1$boxes[, 1]); o4 <- order(r4$tile_index, r4$boxes[, 1])
  expect_equal(r4$boxes[o4, , drop = FALSE], r1$boxes[o1, , drop = FALSE], tolerance = 1e-3)
  expect_equal(r4$class_ids[o4], r1$class_ids[o1])
})

test_that("image_to_rgb_buffer_cpp() is the inverse of rgb_buffer_to_image_cpp()", {
  img <- image_pliman("soybean_touch.jpg")
  d <- dim(img)
  bm <- image_to_rgb_buffer_cpp(image_data(img))
  expect_equal(dim(bm), c(3L, d[1], d[2]))
  back <- rgb_buffer_to_image_cpp(bm, d[1], d[2])
  v <- as.numeric(image_data(img))
  if (max(v) <= 1.5) v <- v * 255          # same scale rule as the converter
  expect_equal(as.integer(back), as.integer(round(v)))
  # grayscale is replicated to three channels
  g <- matrix(c(0, 0.5, 1, 0.25), 2, 2)
  gb <- image_to_rgb_buffer_cpp(g)
  expect_equal(as.integer(gb[1:3]), c(0L, 0L, 0L))
  expect_equal(as.integer(gb[4:6]), rep(128L, 3))
})

test_that("the native batch renderer writes JPEG and PNG and declines other formats", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  img <- image_pliman("soybean_touch.jpg")
  res <- data.frame(xmin = 10, ymin = 10, xmax = 60, ymax = 40, label = "grain", score = 0.9)
  for (ext in c("jpg", "png")) {
    f <- file.path(dir, paste0("out.", ext))
    expect_true(pliman:::.save_detections_native(img, res, f, list()))
    expect_gt(file.size(f), 1000)
    back <- image_import(f)
    expect_equal(dim(back)[1:2], dim(img)[1:2])
  }
  expect_false(pliman:::.save_detections_native(img, res, file.path(dir, "out.tif"), list()))
  withr::local_options(pliman.native_save = FALSE)
  expect_false(pliman:::.save_detections_native(img, res, file.path(dir, "out2.jpg"), list()))
})

test_that("encoder variant selection never puts INT8 on the GPU and falls back to FP32", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("sam2.1.encoder.onnx", "sam2.1.encoder_int8.onnx")))
  base <- file.path(dir, "sam2.1.encoder.onnx")
  withr::local_options(pliman.sam_encoder = "int8")
  expect_match(pliman:::.resolve_sam_encoder(dir, use_gpu = FALSE), "encoder_int8\\.onnx$")
  expect_equal(suppressMessages(pliman:::.resolve_sam_encoder(dir, use_gpu = TRUE)), base)
  withr::local_options(pliman.sam_encoder = "fp16")       # file missing -> FP32
  expect_equal(suppressMessages(pliman:::.resolve_sam_encoder(dir, use_gpu = TRUE)), base)
  withr::local_options(pliman.sam_encoder = NULL)
  expect_equal(pliman:::.resolve_sam_encoder(dir, use_gpu = FALSE), base)
})

test_that("the profiler writes a trace that onnx_profile_summary() can read", {
  skip_if_not_installed("withr")
  skip_if_not_installed("jsonlite")
  lib <- skip_if_no_onnx()
  m <- skip_if_no_model("yolo26n.onnx")
  dir <- withr::local_tempdir()
  withr::local_envvar(c(PLIMAN_ORT_PROFILE = NA))
  pliman_clear_sessions()
  pliman_onnx_profiling(TRUE, dir = dir)
  invisible(run_yolo_cpp(test_frame(), 640, 360, model_path = m, lib_path = lib, use_gpu = FALSE))
  pliman_clear_sessions()                      # releasing the session writes the trace
  pliman_onnx_profiling(FALSE)
  f <- list.files(dir, "\\.json$", full.names = TRUE)
  expect_length(f, 1L)
  s <- onnx_profile_summary(f)
  expect_true(all(c("provider", "op_type", "nodes", "total_ms", "share") %in% names(s)))
  expect_true("CPU" %in% attr(s, "providers")$provider)
  expect_equal(sum(s$share), 1, tolerance = 1e-3)
})
