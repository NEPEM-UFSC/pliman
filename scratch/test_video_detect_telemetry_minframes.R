pkgload::load_all("d:/Desktop/pliman")

cat("=== 1. Checking video_detect_dl formals ===\n")
f_args <- formals(video_detect_dl)
cat("Arguments in video_detect_dl:\n")
stopifnot("telemetry" %in% names(f_args))
stopifnot("as_job" %in% names(f_args))
stopifnot("min_frames" %in% names(f_args))
stopifnot("show_box" %in% names(f_args))
stopifnot("show_point" %in% names(f_args))
cat("telemetry choices:", paste(eval(f_args$telemetry), collapse = ", "), "\n")
cat("as_job default:", deparse(f_args$as_job), "\n")
cat("min_frames default:", deparse(f_args$min_frames), "\n")
cat("[OK] Formals verified!\n\n")

cat("=== 2. Testing min_frames suppression of transient noise in tracking ===\n")
# If min_hits = 5, an object with 1, 2, 3, 4 hits must NOT appear in track_ids or confirmed tracks.
# On hit 5, it must be confirmed and appear!
state <- list()
for (f in 1:6) {
  res <- pliman:::update_tracker_cpp(
    xmin = 100, ymin = 100, xmax = 200, ymax = 200,
    labels = "log", scores = 0.9,
    state = state,
    min_hits = 5L,
    current_time = f * 0.033,
    frame_w = 640L, frame_h = 480L
  )
  state <- res$state
  cat(sprintf("Frame %d: detected 1 box -> confirmed output tracks: %d (hits=%d)\n",
              f, length(res$track_ids), state$tracks[[1]]$hits))
  if (f < 5) {
    stopifnot(length(res$track_ids) == 0) # Suppressed!
  } else {
    stopifnot(length(res$track_ids) == 1) # Confirmed and visible!
  }
}
cat("[OK] min_frames (min_hits) suppression works perfectly! Objects only render once min_frames is reached.\n\n")

cat("=== 3. Testing devtools::document() ===\n")
devtools::document("d:/Desktop/pliman")
cat("[OK] devtools::document() updated man pages successfully!\n")
