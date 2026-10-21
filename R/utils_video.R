#' Split a video into image frames
#'
#' @description
#' `video_split()` extracts individual frames from a video file and saves them
#' directly to disk as image files (`.jpg` or `.png`). It provides native
#' hardware-accelerated decoding on Windows via Windows Media Foundation and C++,
#' saving images directly to disk without intermediate temporary files or
#' external command-line dependencies.
#'
#' On Linux and macOS, or when explicitly requested via `backend = "av"`, it
#' seamlessly delegates decoding to the `av` package.
#'
#' You can control the extraction rate by either:
#' * Specifying `fps` (e.g., `fps = 1` for 1 frame per second).
#' * Specifying `n` for a fixed number of uniformly spaced frames sampled across
#'   the entire video duration.
#' * Leaving both `NULL` to extract every frame of the video.
#'
#' `video_to_frames()` is an alias for `video_split()`.
#'
#' @aliases video_to_frames
#' @param video Character. Path to the video file (e.g., `"video.mp4"`, `"clip.avi"`).
#' @param dir Character. Directory where the extracted frames will be saved.
#'   Defaults to `"frames"`. If the directory does not exist, it will be created.
#' @param fps Numeric. Desired frame rate (frames per second) for extraction.
#'   For example, `fps = 1` extracts one frame per second, `fps = 0.5` extracts
#'   one frame every two seconds. Mutually exclusive with `n`.
#' @param n Integer. Desired total number of frames to sample uniformly across
#'   the entire video duration. For example, `n = 50` will evenly sample 50
#'   frames from the start to the end of the video. Mutually exclusive with `fps`.
#' @param format Character. Output image format: `"jpg"` (or `"jpeg"`) or `"png"`.
#'   Defaults to `"jpg"`.
#' @param quality Integer (1-100). JPEG compression quality. Defaults to `95`.
#'   Only used when `format = "jpg"`.
#' @param prefix Character. Prefix for generated frame filenames. Defaults to
#'   `"frame_"`.
#' @param digits Integer. Number of digits with zero-padding in frame filenames
#'   (e.g., `digits = 6` produces `"frame_000001.jpg"`). Defaults to `6`.
#' @param backend Character. Decoding engine to use: `"auto"` (default),
#'   `"native"` (Windows Media Foundation C++), or `"av"` (using package `av`).
#' @param verbose Logical. If `TRUE` (default), displays extraction progress and
#'   summary.
#'
#' @return A data frame of class `c("pliman_video_frames", "data.frame")` containing:
#' * `frame`: Frame sequential index.
#' * `file`: Full path to the saved image file.
#' * `time`: Timestamp of the frame in seconds.
#' * `width`: Width of the image in pixels.
#' * `height`: Height of the image in pixels.
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examplesIf interactive()
#' # Extract 1 frame per second
#' video_split("my_video.mp4", dir = "frames_fps", fps = 1)
#'
#' # Sample exactly 20 frames uniformly across the video
#' video_split("my_video.mp4", dir = "frames_n", n = 20)
#'
#' # Extract all frames as PNG
#' video_to_frames("my_video.mp4", dir = "all_frames", format = "png")
video_split <- function(video,
                        dir = "frames",
                        fps = NULL,
                        n = NULL,
                        format = c("jpg", "png"),
                        quality = 95,
                        prefix = "frame_",
                        digits = 6,
                        backend = c("auto", "native", "av"),
                        verbose = TRUE) {
  format <- match.arg(format)
  backend <- match.arg(backend)

  if (!file.exists(video)) {
    cli::cli_abort("Video file {.path {video}} does not exist.")
  }

  if (!is.null(fps) && !is.null(n)) {
    cli::cli_abort("Please specify either {.arg fps} or {.arg n}, not both.")
  }

  if (!is.null(fps)) {
    if (!is.numeric(fps) || fps <= 0) {
      cli::cli_abort("{.arg fps} must be a positive number.")
    }
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || n < 1) {
      cli::cli_abort("{.arg n} must be an integer >= 1.")
    }
    n <- as.integer(round(n))
  }

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir_full <- normalizePath(dir, winslash = "/", mustWork = FALSE)
  video_full <- normalizePath(video, winslash = "/", mustWork = TRUE)

  is_windows <- tolower(Sys.info()[["sysname"]]) == "windows"
  can_native <- is_windows && isTRUE(has_native_video_cpp())
  has_av <- requireNamespace("av", quietly = TRUE)

  use_native <- switch(backend,
    "native" = {
      if (!can_native) {
        cli::cli_abort("Native video backend is only supported on Windows with Media Foundation.")
      }
      TRUE
    },
    "av" = {
      if (!has_av) {
        cli::cli_abort(c(
          "Package {.pkg av} is required for this backend.",
          "i" = "Install it with {.run install.packages('av')}."
        ))
      }
      FALSE
    },
    "auto" = {
      # Prefer av if available for maximum hardware SIMD speed (10x faster), otherwise fallback to native
      if (has_av) FALSE else can_native
    }
  )

  t0 <- proc.time()

  mode_txt <- if (!is.null(fps)) {
    paste0(fps, " fps")
  } else if (!is.null(n)) {
    paste0(n, " uniform frames")
  } else {
    "all frames"
  }

  dim_txt <- ""
  dur_val <- NULL
  if (use_native) {
    info <- tryCatch(get_native_video_info_cpp(video_full), error = function(e) NULL)
    if (!is.null(info)) {
      dur_val <- info$duration
      dim_txt <- paste0(" [", info$width, "x", info$height, ", ~", round(info$duration, 1), "s]")
    }
  } else {
    v_info_pre <- tryCatch(av::av_video_info(video_full), error = function(e) NULL)
    if (!is.null(v_info_pre)) {
      dur_val <- if (!is.null(v_info_pre$duration)) v_info_pre$duration else if (!is.null(v_info_pre$video$duration)) v_info_pre$video$duration else NULL
      w_val <- if (!is.null(v_info_pre$video$width) && length(v_info_pre$video$width) > 0) v_info_pre$video$width[1] else 0L
      h_val <- if (!is.null(v_info_pre$video$height) && length(v_info_pre$video$height) > 0) v_info_pre$video$height[1] else 0L
      dur_txt <- if (!is.null(dur_val) && is.numeric(dur_val)) paste0(", ~", round(dur_val, 1), "s") else ""
      if (w_val > 0 && h_val > 0) dim_txt <- paste0(" [", w_val, "x", h_val, dur_txt, "]")
    }
  }

  total_expected <- if (!is.null(n)) {
    as.integer(n)
  } else if (!is.null(fps)) {
    if (!is.null(dur_val) && is.numeric(dur_val) && dur_val > 0) {
      as.integer(max(1, round(dur_val * fps)))
    } else {
      NA
    }
  } else {
    if (!use_native && !is.null(v_info_pre$video$frames)) {
      as.integer(v_info_pre$video$frames[1])
    } else if (use_native && !is.null(info$estimated_frames)) {
      as.integer(info$estimated_frames)
    } else {
      NA
    }
  }

  pb <- NULL
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Extracting video frames ({mode_txt}){dim_txt}...")
    if (!use_native) {
      fmt_str <- if (!is.na(total_expected)) {
        "{cli::pb_spin} Extracting frames [{cli::pb_bar}] {cli::pb_current}/{cli::pb_total} ({cli::pb_percent}) [ETA: {cli::pb_eta}]"
      } else {
        "{cli::pb_spin} Extracting frames [{cli::pb_bar}] {cli::pb_current} frames [Elapsed: {cli::pb_elapsed}]"
      }
      pb <- cli::cli_progress_bar(
        total = total_expected,
        format = fmt_str,
        clear = FALSE
      )
    }
  }

  if (use_native) {
    target_fps <- if (is.null(fps)) 0.0 else as.numeric(fps)
    target_n <- if (is.null(n)) 0L else as.integer(n)

    res <- split_video_native_cpp(
      video_path = video_full,
      output_dir = dir_full,
      prefix = prefix,
      digits = as.integer(digits),
      format = format,
      quality = as.integer(quality),
      target_fps = target_fps,
      target_n = target_n,
      verbose = isTRUE(verbose)
    )
  } else {
    v_info <- av::av_video_info(video_full)
    dur <- if (!is.null(v_info$duration)) v_info$duration else if (!is.null(v_info$video$duration)) v_info$video$duration else NULL
    framerate <- if (!is.null(v_info$video$framerate) && length(v_info$video$framerate) > 0) v_info$video$framerate[1] else 30
    w <- if (!is.null(v_info$video$width) && length(v_info$video$width) > 0) v_info$video$width[1] else 0L
    h <- if (!is.null(v_info$video$height) && length(v_info$video$height) > 0) v_info$video$height[1] else 0L

    av_fps <- fps
    if (!is.null(n)) {
      if (!is.null(dur) && is.numeric(dur) && dur > 0) {
        av_fps <- if (n == 1) (1 / dur) else ((n / dur) * 1.02)
      }
    }

    # Build filter
    vfilter <- if (!is.null(av_fps) && av_fps > 0) {
      paste0("fps=fps=", av_fps)
    } else {
      "null"
    }

    av_ext <- if (format == "jpg") "jpg" else "png"
    av_codec <- if (format == "jpg") "mjpeg" else "png"
    out_pattern <- file.path(dir_full, paste0(prefix, "%0", digits, "d.", av_ext))
    pat <- paste0("^", prefix, "\\d{", digits, ",}\\.", av_ext, "$")

    old_log_level <- av::av_log_level(-8)
    on.exit(av::av_log_level(old_log_level), add = TRUE)

    av::av_encode_video(
      input = video_full,
      output = out_pattern,
      framerate = framerate,
      codec = av_codec,
      vfilter = vfilter,
      verbose = FALSE
    )
    av::av_log_level(old_log_level)

    av_files <- list.files(dir_full, pattern = pat, full.names = TRUE)

    if (!is.null(n) && length(av_files) > n) {
      excess <- av_files[(n + 1):length(av_files)]
      unlink(excess)
      av_files <- av_files[1:n]
    }

    if (!is.null(pb)) {
      tryCatch(cli::cli_progress_update(id = pb, set = length(av_files)), error = function(e) NULL)
      tryCatch(cli::cli_progress_done(id = pb), error = function(e) NULL)
    }

    res <- data.frame(
      frame = seq_along(av_files),
      file = normalizePath(av_files, winslash = "/", mustWork = FALSE),
      time = if (!is.null(n) && !is.null(dur) && dur > 0) {
        if (n == 1) 0 else round(seq(0, dur, length.out = length(av_files)), 3)
      } else if (!is.null(av_fps) && av_fps > 0) {
        round((seq_along(av_files) - 1) / av_fps, 3)
      } else if (!is.null(framerate) && framerate > 0) {
        round((seq_along(av_files) - 1) / framerate, 3)
      } else {
        NA_real_
      },
      width = rep(as.integer(w), length(av_files)),
      height = rep(as.integer(h), length(av_files)),
      stringsAsFactors = FALSE
    )
  }

  elapsed <- round(as.numeric((proc.time() - t0)[3]), 2)

  if (isTRUE(verbose)) {
    cli::cli_alert_success(
      "Successfully extracted {nrow(res)} frame{?s} to {.path {dir}} in {elapsed}s."
    )
  }

  attr(res, "video") <- video_full
  attr(res, "dir") <- dir_full
  attr(res, "format") <- format
  attr(res, "fps") <- fps
  attr(res, "n") <- n
  attr(res, "elapsed") <- elapsed
  class(res) <- c("pliman_video_frames", "data.frame")

  invisible(res)
}

#' @export
video_to_frames <- video_split

#' Print method for pliman_video_frames
#' @param x An object of class `pliman_video_frames`.
#' @param ... Additional arguments (not used).
#' @export
print.pliman_video_frames <- function(x, ...) {
  n_frames <- nrow(x)
  dir_path <- attr(x, "dir")
  vid_path <- attr(x, "video")
  fmt <- attr(x, "format")
  elapsed <- attr(x, "elapsed")

  cli::cli_rule(left = "Pliman Video Frames")
  cli::cli_bullets(c(
    "*" = paste0("Video: ", basename(vid_path)),
    "*" = paste0("Output dir: ", dir_path),
    "*" = paste0("Frames extracted: ", n_frames, " (", toupper(fmt), ")"),
    "*" = if (n_frames > 0) paste0("Resolution: ", x$width[1], " x ", x$height[1]) else NULL,
    "*" = if (!is.null(elapsed)) paste0("Elapsed time: ", elapsed, "s") else NULL
  ))
  cli::cli_rule()
  print.data.frame(head(x, 10))
  if (n_frames > 10) {
    cli::cli_alert_info("... with {n_frames - 10} more rows. Use as.data.frame() to view all.")
  }
  invisible(x)
}

#' Extract video frames directly into Image objects
#'
#' @description
#' `frame_to_img()` extracts one or more frames from a saved video file directly
#' into memory as an object of class `image` (or a named list of `image` objects
#' if multiple frames are requested).
#'
#' On Windows, decoding runs with zero disk overhead via high-speed native
#' Media Foundation C++ integration. On other operating systems (or when native
#' reader is unavailable), it seamlessly falls back to `av`.
#'
#' @param video Character. Path to the video file (e.g., `"video.mp4"`).
#' @param frames Frame index or indices to extract (1-indexed). Can be:
#'   * A single integer (e.g., `1` for the first frame, `50` for frame 50).
#'   * A vector of integers (e.g., `c(1, 15, 30)` or `1:5`).
#'   * `"first"` for the first frame (`1L`).
#'   * `"last"` for the last frame of the video.
#'   * `"all"` for all frames in the video.
#'   Defaults to `1L`. Mutually exclusive with `time`.
#' @param time Numeric. Video timestamp(s) in seconds from which to extract frames
#'   (e.g., `time = 1.5` or `c(0.5, 2.0)`). Mutually exclusive with `frames`.
#' @param show Logical. If `TRUE`, plots the extracted image(s). Default is `FALSE`.
#' @param output_dir Optional character. Directory where extracted images should
#'   be saved to disk. If `NULL` (default), images are kept in memory only.
#' @param verbose Logical. If `TRUE`, prints informative messages. Default is `FALSE`.
#' @param ... Additional arguments passed to [plot.image()] or [image_combine()]
#'   when `show = TRUE`.
#'
#' @return An object of class `image` (if a single frame is requested) or a named
#'   list of `image` objects (if multiple frames are requested, named
#'   `"frame_1"`, `"frame_2"`, etc.).
#'
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examplesIf interactive()
#' # Extract the first frame of a video
#' img <- frame_to_img("my_video.mp4")
#' plot(img)
#'
#' # Extract frame at 2.5 seconds
#' img <- frame_to_img("my_video.mp4", time = 2.5)
#'
#' # Extract multiple frames as a list of images
#' imgs <- frame_to_img("my_video.mp4", frames = c(10, 30, 50))
#' image_combine(imgs)
frame_to_img <- function(video,
                         frames = 1L,
                         time = NULL,
                         show = FALSE,
                         output_dir = NULL,
                         verbose = FALSE,
                         ...) {
  if (!is.character(video) || length(video) != 1L || !file.exists(video)) {
    cli::cli_abort("Video file {.path {video}} does not exist.")
  }

  video_full <- normalizePath(video, winslash = "/", mustWork = TRUE)
  is_windows <- tolower(Sys.info()[["sysname"]]) == "windows"
  can_native <- is_windows && isTRUE(has_native_video_cpp())

  reader <- if (can_native) tryCatch(open_native_video_reader_cpp(video_full), error = function(e) NULL) else NULL
  use_native <- !is.null(reader)

  # Determine total frames and fps
  fps_val <- 30.0
  tot_val <- 1000000L

  if (use_native) {
    r_info <- get_video_reader_info_cpp(reader)
    if (!is.null(r_info$fps) && r_info$fps > 0) fps_val <- as.numeric(r_info$fps)
    if (!is.null(r_info$total_frames) && r_info$total_frames > 0) tot_val <- as.integer(r_info$total_frames)
    w_val <- as.integer(r_info$width)
    h_val <- as.integer(r_info$height)
  } else {
    v_info <- tryCatch(suppressWarnings(av::av_video_info(video_full)), error = function(e) NULL)
    if (!is.null(v_info)) {
      if (!is.null(v_info$video$framerate) && v_info$video$framerate[1] > 0) fps_val <- as.numeric(v_info$video$framerate[1])
      if (!is.null(v_info$video$frames) && v_info$video$frames[1] > 0) tot_val <- as.integer(v_info$video$frames[1])
      w_val <- if (!is.null(v_info$video$width)) as.integer(v_info$video$width[1]) else 0L
      h_val <- if (!is.null(v_info$video$height)) as.integer(v_info$video$height[1]) else 0L
    }
  }

  # Parse requested frames
  if (!is.null(time)) {
    if (!is.numeric(time) || any(time < 0)) {
      if (use_native) close_native_video_reader_cpp(reader)
      cli::cli_abort("{.arg time} must be a positive numeric value in seconds.")
    }
    raw_frames <- as.integer(pmax(1L, round(as.numeric(time) * fps_val) + 1L))
  } else if (is.character(frames) && length(frames) == 1L) {
    f_mode <- tolower(frames)
    raw_frames <- switch(f_mode,
      "first" = 1L,
      "last" = if (tot_val < 1000000L) tot_val else 1L,
      "all" = if (tot_val < 1000000L) seq_len(tot_val) else seq_len(100L),
      {
        val <- suppressWarnings(as.integer(frames))
        if (is.na(val)) 1L else val
      }
    )
  } else {
    raw_frames <- as.integer(frames)
  }

  is_single <- (length(raw_frames) == 1L)
  target_unique <- sort(unique(raw_frames[raw_frames >= 1L]))

  if (length(target_unique) == 0L) {
    if (use_native) close_native_video_reader_cpp(reader)
    cli::cli_abort("No valid frame indices requested (indices must be >= 1).")
  }

  imgs_map <- list()

  if (use_native) {
    for (f_idx in target_unique) {
      img <- extract_native_frame_to_image_cpp(reader, as.integer(f_idx))
      if (!is.null(img)) {
        imgs_map[[as.character(f_idx)]] <- img
      }
    }
    close_native_video_reader_cpp(reader)
  } else {
    # Cross-platform fallback via av
    tmp_dir <- file.path(tempdir(), paste0("pliman_frame_to_img_", as.integer(stats::runif(1, 100000, 999999))))
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    suppressWarnings(av::av_video_images(video_full, destdir = tmp_dir, format = "jpg"))
    img_files <- list.files(tmp_dir, pattern = "\\.jpe?g$", full.names = TRUE)

    for (f_idx in target_unique) {
      if (f_idx <= length(img_files)) {
        imgs_map[[as.character(f_idx)]] <- image_import(img_files[f_idx])
      }
    }
  }

  if (length(imgs_map) == 0L) {
    cli::cli_abort("Failed to extract any frames from {.path {video}}.")
  }

  # Preserve requested ordering
  result_list <- list()
  for (rf in raw_frames) {
    key <- as.character(rf)
    if (!is.null(imgs_map[[key]])) {
      result_list[[paste0("frame_", rf)]] <- imgs_map[[key]]
    }
  }

  # Optional export to disk
  if (!is.null(output_dir) && is.character(output_dir) && nzchar(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(result_list)) {
      out_f <- file.path(output_dir, paste0(nm, ".jpg"))
      image_export(result_list[[nm]], out_f, verbose = FALSE)
    }
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Saved {length(result_list)} frame{?s} to {.path {output_dir}}.")
    }
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Extracted {length(result_list)} frame{?s} from {.file {basename(video)}}.")
  }

  # Return single image or list of images
  if (is_single && length(result_list) == 1L) {
    out_obj <- result_list[[1L]]
    if (isTRUE(show)) {
      plot(out_obj, ...)
    }
    return(out_obj)
  }

  if (isTRUE(show)) {
    tryCatch({
      image_combine(result_list, ...)
    }, error = function(e) {
      plot(result_list[[1L]], ...)
    })
  }

  return(result_list)
}
