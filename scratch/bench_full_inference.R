devtools::load_all("d:/Desktop/pliman")

onnx_path <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"
lib_path <- pliman:::pliman_onnx_library_path()

# Create dummy image 1920x1080 (3 x 1920 x 1080 raw bitmap)
w <- 1920
h <- 1080
bm <- raw(3 * w * h)
dim(bm) <- c(3, w, h)

# Warmup
for (i in 1:10) {
  t_vec <- pliman:::preprocess_yolo_cpp(bm, target_size = 640L)
  dummy <- pliman:::run_yolo_cpp(t_vec, w, h, 0.5, 0.5, onnx_path, lib_path, use_gpu = TRUE, need_masks = FALSE)
}

n_iter <- 50
t0 <- Sys.time()
t_pre_total <- 0
t_inf_total <- 0

for (i in 1:n_iter) {
  t_a <- Sys.time()
  t_vec <- pliman:::preprocess_yolo_cpp(bm, target_size = 640L)
  t_b <- Sys.time()
  dummy <- pliman:::run_yolo_cpp(t_vec, w, h, 0.5, 0.5, onnx_path, lib_path, use_gpu = TRUE, need_masks = FALSE)
  t_c <- Sys.time()
  
  t_pre_total <- t_pre_total + as.numeric(t_b - t_a, units = "secs")
  t_inf_total <- t_inf_total + as.numeric(t_c - t_b, units = "secs")
}

ms_pre <- (t_pre_total / n_iter) * 1000
ms_inf <- (t_inf_total / n_iter) * 1000
ms_total <- ms_pre + ms_inf
fps_total <- 1000.0 / ms_total

cat(sprintf("pliman com Pré-processamento de imagem 1080p:\n"))
cat(sprintf("  - Preprocess 1080p -> 640x640: %.2f ms\n", ms_pre))
cat(sprintf("  - Inferência GPU ONNX:        %.2f ms\n", ms_inf))
cat(sprintf("  - Total Preprocess + Infer:   %.2f ms  (%.1f FPS)\n", ms_total, fps_total))
