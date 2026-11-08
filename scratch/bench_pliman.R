devtools::load_all("d:/Desktop/pliman")

onnx_path <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"
lib_path <- pliman:::pliman_onnx_library_path()

cat("=================================================================\n")
cat("BENCHMARK PLIMAN ONNX RUNTIME (PURAMENTE INFERÊNCIA)\n")
cat("=================================================================\n")

# Prepare dummy tensor letterboxed 640x640 (3 x 640 x 640)
n_floats <- 3 * 640 * 640
tensor <- rep(114.0 / 255.0, n_floats)
attr(tensor, "target_size") <- 640L
attr(tensor, "orig_w") <- 640L
attr(tensor, "orig_h") <- 640L
attr(tensor, "gain") <- 1.0
attr(tensor, "pad_x") <- 0.0
attr(tensor, "pad_y") <- 0.0

n_iter <- 100

# 1. PLIMAN GPU (DirectML)
cat("Testando pliman run_yolo_cpp (GPU DirectML)...\n")
# Warmup
for (i in 1:15) {
  dummy <- pliman:::run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = 640,
    orig_h = 640,
    conf_threshold = 0.5,
    iou_threshold = 0.5,
    model_path = onnx_path,
    lib_path = lib_path,
    use_gpu = TRUE,
    need_masks = FALSE
  )
}

t0 <- Sys.time()
for (i in 1:n_iter) {
  res_gpu <- pliman:::run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = 640,
    orig_h = 640,
    conf_threshold = 0.5,
    iou_threshold = 0.5,
    model_path = onnx_path,
    lib_path = lib_path,
    use_gpu = TRUE,
    need_masks = FALSE
  )
}
t_gpu_total <- as.numeric(Sys.time() - t0, units = "secs")
t_gpu_ms <- (t_gpu_total / n_iter) * 1000
fps_gpu <- 1000.0 / t_gpu_ms

cat(sprintf("1) pliman ONNXRuntime (best.onnx | GPU DirectML):\n"))
cat(sprintf("   - run_yolo_cpp() puro       : %.2f ms/frame  (%.1f FPS)\n", t_gpu_ms, fps_gpu))

# 2. PLIMAN CPU
cat("\nTestando pliman run_yolo_cpp (CPU)...\n")
for (i in 1:5) {
  dummy <- pliman:::run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = 640,
    orig_h = 640,
    conf_threshold = 0.5,
    iou_threshold = 0.5,
    model_path = onnx_path,
    lib_path = lib_path,
    use_gpu = FALSE,
    num_threads = 4L,
    need_masks = FALSE
  )
}

n_iter_cpu <- 30
t0 <- Sys.time()
for (i in 1:n_iter_cpu) {
  res_cpu <- pliman:::run_yolo_cpp(
    tensor_vec = tensor,
    orig_w = 640,
    orig_h = 640,
    conf_threshold = 0.5,
    iou_threshold = 0.5,
    model_path = onnx_path,
    lib_path = lib_path,
    use_gpu = FALSE,
    num_threads = 4L,
    need_masks = FALSE
  )
}
t_cpu_total <- as.numeric(Sys.time() - t0, units = "secs")
t_cpu_ms <- (t_cpu_total / n_iter_cpu) * 1000
fps_cpu <- 1000.0 / t_cpu_ms

cat(sprintf("2) pliman ONNXRuntime (best.onnx | CPU 4 threads):\n"))
cat(sprintf("   - run_yolo_cpp() CPU        : %.2f ms/frame  (%.1f FPS)\n", t_cpu_ms, fps_cpu))
cat("=================================================================\n")
