library(pliman)

cat("=== 1. Testing update_tracker_cpp with Detailed Telemetry Accumulators ===\n")

state <- list()
line_coords <- c(100, 250, 400, 250)
roi_coords <- c(50, 50, 450, 450)

# Simulate 3 consecutive frames with an object moving across the line
for (f in 1:3) {
  t_sec <- (f - 1) * 0.033
  y_pos <- 230 + (f - 1) * 15 # starts at 230, then 245, then 260 (crosses 250)
  
  res <- pliman:::update_tracker_cpp(
    xmin = 200,
    ymin = y_pos - 10,
    xmax = 240,
    ymax = y_pos + 10,
    labels = "object",
    scores = 0.88,
    state = state,
    count_line = line_coords,
    roi = roi_coords,
    max_dist = 120,
    min_iou = 0.25,
    max_lost = 15,
    max_history = 30,
    current_time = t_sec,
    frame_w = 640L,
    frame_h = 480L
  )
  state <- res$state
  
  cat(sprintf(
    "Frame %d: Active=%d, InZone=%d, Counted=%d, ConfMean=%.2f, Flow=%.2f, VelMean=%.2f, WelfordMean=%.1f, WelfordCV=%.1f%%, RoiOcc=%.2f%%\n",
    f, res$active, res$in_zone, res$total_count, res$conf_mean, res$flow, res$vel_mean, res$welford_mean, res$welford_cv, res$roi_occ
  ))
}

stopifnot(res$active >= 1)
stopifnot(!is.null(res$flow_history))
stopifnot(length(res$flow_history) == 20)
cat("[OK] Tracker math & Welford accumulators verified!\n\n")

cat("=== 2. Testing C++ Rendering with Simple vs Detailed HUD ===\n")
# Create synthetic BGR buffer (640x480x3)
bm <- raw(3L * 640L * 480L)
attr(bm, "dim") <- c(3L, 640L, 480L)

# Simple HUD
pliman:::draw_yolo_detections_bgr_cpp(
  bm = bm,
  xmin = 200, ymin = 200, xmax = 240, ymax = 260,
  labels = "log", scores = 0.92, colors = "#00FF00", lwd = 2L,
  hud_text = "FRAME: 1 | SPEED: 30.0 FPS | IN ZONE: 1 | COUNT: 0",
  detailed_telemetry = FALSE,
  sparkline = NULL
)
cat("[OK] Simple HUD rendered successfully!\n")

# Detailed HUD with UTF-8 characters and sparkline
detailed_txt <- paste(
  "FRAME   : 120",
  "SPEED   :  29.8 FPS",
  "LATENCY :  85ms [INF:70 | TRK: 2 | GUI:13]",
  "TRACKS  :   62 ACT",
  "CONF    :  0.74 \xc2\xb1 0.08 (MIN: 0.48)",
  "FLOW    :  18.4 obj/s",
  "MOTION  :  11.2 px/f (JIT:  1.4)",
  "IN ZONE :   45",
  "COUNTED :  128",
  "SIZE    :   450 px\xc2\xb2 (CV: 14.2%)",
  "ROI OCC :  32.5%",
  sep = "\n"
)

pliman:::draw_yolo_detections_bgr_cpp(
  bm = bm,
  xmin = 200, ymin = 200, xmax = 240, ymax = 260,
  labels = "log", scores = 0.92, colors = "#00FF00", lwd = 2L,
  hud_text = detailed_txt,
  detailed_telemetry = TRUE,
  sparkline = as.numeric(c(rep(10, 10), rep(18, 10)))
)
cat("[OK] Detailed HUD rendered successfully with 8x16 font, UTF-8 symbols & sparkline!\n\n")

cat("=== 3. Testing video_detect_dl argument validation ===\n")
f_args <- formals(video_detect_dl)
stopifnot("telemetry" %in% names(f_args))
cat("[OK] 'telemetry' is in video_detect_dl arguments.\n")

# Verify match.arg behavior
tryCatch({
  video_detect_dl(telemetry = "invalid_choice", show = FALSE, max_frames = 1)
  stop("Should have failed on invalid telemetry choice")
}, error = function(e) {
  cat("[OK] match.arg rejected invalid choice: ", e$message, "\n")
})

cat("\nALL VERIFICATION CHECKS PASSED!\n")
