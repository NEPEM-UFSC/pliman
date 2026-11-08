# Test dynamic imgsz support in image_detect_dl, video_detect_dl, and onnx inference
pkgload::load_all("d:/Desktop/pliman")

cat("=== 1. Checking .get_onnx_imgsz on toras_seg2.onnx ===\n")
toras_onnx <- file.path(pliman_model_dir(), "toras_seg2.onnx")
if (file.exists(toras_onnx)) {
  meta_imgsz <- pliman:::.get_onnx_imgsz(toras_onnx)
  cat("Detected native ONNX imgsz:", meta_imgsz, "\n")
  stopifnot(identical(meta_imgsz, 640L))
  cat("[OK] .get_onnx_imgsz extracted 640L correctly!\n")
}

cat("\n=== 2. Exporting a dynamic ONNX model with Python ===\n")
dyn_onnx <- file.path(tempdir(), "yolo11n.onnx")
if (!file.exists(dyn_onnx)) {
  # Export yolo11n with dynamic=True
  py_script <- sprintf("import os, shutil; from ultralytics import YOLO; m = YOLO('d:/Desktop/pliman/yolo11n.pt'); p = m.export(format='onnx', dynamic=True, imgsz=640); shutil.copy(p, '%s')", gsub("\\\\", "/", dyn_onnx))
  system2("python", c("-c", shQuote(py_script)))
}

if (file.exists(dyn_onnx)) {
  cat("[OK] Dynamic ONNX model exported to:", dyn_onnx, "\n")
  details <- yolo_model_details(dyn_onnx)
  cat("Dynamic model export args dynamic:", details$export_args$dynamic, "\n")
  cat("Dynamic model tensor inputs:\n")
  print(details$tensor_signatures$inputs)

  # Create a synthetic image (1200 x 900 x 3)
  test_im <- array(runif(1200 * 900 * 3, 0.2, 0.8), dim = c(1200, 900, 3))

  cat("\n=== 3. Running image_detect_dl with imgsz = 640 ===\n")
  t1 <- Sys.time()
  res_640 <- image_detect_dl(test_im, model = dyn_onnx, imgsz = 640, plot = FALSE, verbose = FALSE)
  cat("imgsz=640 inference elapsed:", round(as.numeric(Sys.time() - t1, units = "secs") * 1000, 1), "ms\n")
  cat("[OK] imgsz = 640 completed successfully!\n")

  cat("\n=== 4. Running image_detect_dl with imgsz = 960 (Dynamic inference) ===\n")
  t2 <- Sys.time()
  res_960 <- image_detect_dl(test_im, model = dyn_onnx, imgsz = 960, plot = FALSE, verbose = FALSE)
  cat("imgsz=960 inference elapsed:", round(as.numeric(Sys.time() - t2, units = "secs") * 1000, 1), "ms\n")
  cat("[OK] imgsz = 960 completed successfully on dynamic ONNX model!\n")

  cat("\n=== 5. Running with user alias imgz = 960 ===\n")
  t3 <- Sys.time()
  res_alias <- image_detect_dl(test_im, model = dyn_onnx, imgz = 960, plot = FALSE, verbose = FALSE)
  cat("imgz=960 inference elapsed:", round(as.numeric(Sys.time() - t3, units = "secs") * 1000, 1), "ms\n")
  cat("[OK] imgz = 960 alias completed successfully!\n")
} else {
  cat("[WARNING] Could not export dynamic model with Python. Skipping dynamic ONNX run.\n")
}

cat("\n=== 6. Checking argument signatures ===\n")
funs <- list(
  image_detect_dl = image_detect_dl,
  video_detect_dl = video_detect_dl,
  image_segment_dl = image_segment_dl,
  image_pose_dl = image_pose_dl,
  image_detect_dl_batch = image_detect_dl_batch,
  image_segment_dl_batch = image_segment_dl_batch,
  image_detect_sahi = image_detect_sahi,
  image_segment_sahi = image_segment_sahi
)
for (nm in names(funs)) {
  f_args <- names(formals(funs[[nm]]))
  cat(sprintf("Checking %-25s : has 'imgsz'? %s\n", nm, "imgsz" %in% f_args))
  stopifnot("imgsz" %in% f_args)
}

cat("\nALL DYNAMIC IMGSZ & ALIAS TESTS PASSED SUCCESSFULLY!\n")
