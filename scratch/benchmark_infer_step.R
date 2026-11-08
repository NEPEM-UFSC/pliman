devtools::load_all("d:/Desktop/pliman")
res_model <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"

cat("=== BENCHMARK: infer_every = 1L (cada frame) ===\n")
t1 <- Sys.time()
res1 <- video_detect_dl(
  conf_threshold = 0.5,
  video = "D:/Downloads/pliman4/toras2.mp4",
  model = res_model,
  roi = c(0.5, 1),
  roi_orientation = "v",
  count_line = 0.5,
  count_line_orientation = "v",
  track = TRUE,
  trail = FALSE,
  infer_every = 1L,
  output = "D:/Downloads/pliman4/toras2_infer1.mp4",
  show_class = TRUE,
  telemetry = "detailed",
  hud_pos = "top-left",
  min_frames = 20L,
  show_text = FALSE,
  as_job = FALSE,
  max_frames = 60,
  verbose = FALSE
)
el1 <- as.numeric(Sys.time() - t1, units = "secs")
fps1 <- 60 / el1
cat(sprintf("infer_every = 1L: %.2f s (%.1f FPS) | Counted: %s\n", el1, fps1, res1$total_count))

cat("\n=== BENCHMARK: infer_every = 2L (frame alternado com propagação de velocidade) ===\n")
t2 <- Sys.time()
res2 <- video_detect_dl(
  conf_threshold = 0.5,
  video = "D:/Downloads/pliman4/toras2.mp4",
  model = res_model,
  roi = c(0.5, 1),
  roi_orientation = "v",
  count_line = 0.5,
  count_line_orientation = "v",
  track = TRUE,
  trail = FALSE,
  infer_every = 2L,
  output = "D:/Downloads/pliman4/toras2_infer2.mp4",
  show_class = TRUE,
  telemetry = "detailed",
  hud_pos = "top-left",
  min_frames = 20L,
  show_text = FALSE,
  as_job = FALSE,
  max_frames = 60,
  verbose = FALSE
)
el2 <- as.numeric(Sys.time() - t2, units = "secs")
fps2 <- 60 / el2
cat(sprintf("infer_every = 2L: %.2f s (%.1f FPS) | Counted: %s\n", el2, fps2, res2$total_count))
