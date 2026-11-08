devtools::load_all("d:/Desktop/pliman")
res_model <- "D:/Downloads/pliman4/wood_log_model/weights/best.onnx"
res <- video_detect_dl(
  conf_threshold = 0.5,
  video = "D:/Downloads/pliman4/toras2.mp4",
  model = res_model,
  roi = c(0.5, 1),
  roi_orientation = "v",
  roi_label = "WOOD LOG DETECTOR",
  count_line = 0.5,
  count_line_orientation = "v",
  track = TRUE,
  trail = FALSE,
  output = "D:/Downloads/pliman4/toras2_track_test.mp4",
  text_size_objects = 2,
  text_size_telemetry = 2,
  text_size_roi = 2,
  show_class = TRUE,
  telemetry = "detailed",
  hud_pos = "top-left",
  min_frames = 20L,
  show_text = FALSE,
  as_job = FALSE,
  max_frames = 60
)
cat("SUCCESS! Counted:", res$total_count, "\n")
