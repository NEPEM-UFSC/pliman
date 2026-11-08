devtools::load_all(".")
setwd("D:/Downloads/pliman4")
res_model <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"

cat("Testing inference times:\n")
f_path <- "D:/Downloads/pliman4/forest_guard_rf/train/images/ds1_001_jpg.rf.859664dae638dbf614457f920f023f00.jpg"
if (!file.exists(f_path)) {
  # find an image
  imgs <- list.files("D:/Downloads/pliman4", pattern = "\\.jpg$", recursive = TRUE, full.names = TRUE)
  f_path <- imgs[1]
}
cat("Image:", f_path, "\n")
bm <- load_frame_bgr_cpp(f_path)
dims <- dim(bm)
w <- dims[2]
h <- dims[3]

# 1. 640x640 with GPU
t1 <- preprocess_yolo_cpp(bm, target_size = 640)
t0 <- Sys.time()
for(k in 1:5) {
  r1 <- run_yolo_cpp(t1, w, h, 0.25, 0.45, res_model, pliman_onnx_library_path(), 0L, use_gpu = TRUE, device_id = -1, need_masks = FALSE)
}
time_gpu_640 <- as.numeric(Sys.time() - t0)/5 * 1000
cat("GPU 640ms:", time_gpu_640, "ms (", 1000/time_gpu_640, "FPS)\n")

# 2. 1024x1024 with GPU
t2 <- preprocess_yolo_cpp(bm, target_size = 1024)
t0 <- Sys.time()
for(k in 1:5) {
  r2 <- run_yolo_cpp(t2, w, h, 0.25, 0.45, res_model, pliman_onnx_library_path(), 0L, use_gpu = TRUE, device_id = -1, need_masks = FALSE)
}
time_gpu_1024 <- as.numeric(Sys.time() - t0)/5 * 1000
cat("GPU 1024ms:", time_gpu_1024, "ms (", 1000/time_gpu_1024, "FPS)\n")

# 3. What if need_masks = TRUE?
t0 <- Sys.time()
for(k in 1:5) {
  r3 <- run_yolo_cpp(t1, w, h, 0.25, 0.45, res_model, pliman_onnx_library_path(), 0L, use_gpu = TRUE, device_id = -1, need_masks = TRUE)
}
time_gpu_mask <- as.numeric(Sys.time() - t0)/5 * 1000
cat("GPU need_masks=TRUE:", time_gpu_mask, "ms\n")

# 4. What about CPU?
t0 <- Sys.time()
for(k in 1:5) {
  r4 <- run_yolo_cpp(t1, w, h, 0.25, 0.45, res_model, pliman_onnx_library_path(), 0L, use_gpu = FALSE, device_id = -1, need_masks = FALSE)
}
time_cpu_640 <- as.numeric(Sys.time() - t0)/5 * 1000
cat("CPU 640ms:", time_cpu_640, "ms\n")
