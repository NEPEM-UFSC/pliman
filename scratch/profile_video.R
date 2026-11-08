
pkgload::load_all("d:/Desktop/pliman")
vid <- "D:/Downloads/pliman4/toras2.mp4"
cat("File exists:", file.exists(vid), "\n")
cat("has_native_video_reader_cpp:", pliman:::has_native_video_reader_cpp(), "\n")
cat("has_native_video_writer_cpp:", pliman:::has_native_video_writer_cpp(), "\n")
vr <- tryCatch(pliman:::open_native_video_reader_cpp(normalizePath(vid, winslash="/")), error=function(e) e$message)
if (inherits(vr, "externalptr")) {
  info <- pliman:::get_video_reader_info_cpp(vr)
  cat(sprintf("Native video reader info: %dx%d, fps=%.2f, total_frames=%d, duration=%.2fs\n",
              info$width, info$height, info$fps, info$total_frames, info$duration))
  pliman:::close_native_video_reader_cpp(vr)
} else {
  cat("Native video reader error/fallback:", vr, "\n")
}
