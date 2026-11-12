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
      if (can_native) {
        TRUE
      } else if (has_av) {
        FALSE
      } else {
        cli::cli_abort(c(
          "No video decoding backend available.",
          "i" = "Install {.pkg av} with {.run install.packages('av')}."
        ))
      }
    }
  )

  t0 <- proc.time()

  av_ext <- if (format == "jpg") "jpg" else "png"
  pat <- paste0("^", prefix, "\\d{", digits, ",}\\.", av_ext, "$")
  old_matching <- list.files(dir_full, pattern = pat, full.names = TRUE)
  if (length(old_matching) > 0L) {
    unlink(old_matching)
  }

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
    } else if (backend == "auto" && has_av) {
      use_native <- FALSE
    }
  }
  if (!use_native) {
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
      NA_integer_
    }
  } else {
    if (!use_native && !is.null(v_info_pre$video$frames)) {
      as.integer(v_info_pre$video$frames[1])
    } else if (use_native && !is.null(info$estimated_frames)) {
      as.integer(info$estimated_frames)
    } else {
      NA_integer_
    }
  }

  pb <- NULL
  if (isTRUE(verbose)) {
    cli::cli_alert_info("Extracting video frames ({mode_txt}){dim_txt}...")
    fmt_str <- if (!is.na(total_expected)) {
      "{cli::pb_spin} Extracting frames [{cli::pb_bar}] {cli::pb_current}/{cli::pb_total} ({cli::pb_percent}) [ETA: {cli::pb_eta}]"
    } else {
      "{cli::pb_spin} Extracting frames [{cli::pb_bar}] {cli::pb_current} frames [Elapsed: {cli::pb_elapsed}]"
    }
    pb <- cli::cli_progress_bar(
      total = total_expected,
      format = fmt_str,
      clear = FALSE,
      auto_terminate = FALSE
    )
  }

  res <- NULL
  if (use_native) {
    target_fps <- if (is.null(fps)) 0.0 else as.numeric(fps)
    target_n <- if (is.null(n)) 0L else as.integer(n)
    cur_frame <- 0L
    cb_fn <- if (!is.null(pb)) {
      function() {
        cur_frame <<- cur_frame + 1L
        tot_now <- if (!is.na(total_expected)) max(total_expected, cur_frame) else NA_integer_
        tryCatch(
          cli::cli_progress_update(id = pb, set = cur_frame, total = tot_now),
          error = function(e) NULL
        )
      }
    } else {
      NULL
    }

    res <- tryCatch(
      split_video_native_cpp(
        video_path = video_full,
        output_dir = dir_full,
        prefix = prefix,
        digits = as.integer(digits),
        format = format,
        quality = as.integer(quality),
        target_fps = target_fps,
        target_n = target_n,
        verbose = FALSE,
        callback = cb_fn
      ),
      error = function(e) {
        if (backend == "auto" && has_av) {
          NULL
        } else {
          stop(e)
        }
      }
    )
    if (is.null(res)) {
      use_native <- FALSE
    }
  }

  if (!use_native) {
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

    old_log_level <- av::av_log_level(-8)
    on.exit(av::av_log_level(old_log_level), add = TRUE)

    if (format == "jpg") {
      # Decode via lossless uncompressed BMP first to avoid av's hardcoded 200kbps MJPEG bitrate artifacting,
      # then encode to high-quality JPEG in parallel C++ respecting `quality`.
      tmp_av_dir <- file.path(tempdir(), paste0("pliman_av_split_", as.integer(stats::runif(1, 100000, 999999))))
      dir.create(tmp_av_dir, recursive = TRUE, showWarnings = FALSE)
      on.exit(unlink(tmp_av_dir, recursive = TRUE, force = TRUE), add = TRUE)

      bmp_pattern <- file.path(tmp_av_dir, paste0(prefix, "%0", digits, "d.bmp"))
      bmp_pat <- paste0("^", prefix, "\\d{", digits, ",}\\.bmp$")

      av::av_encode_video(
        input = video_full,
        output = bmp_pattern,
        framerate = framerate,
        codec = "bmp",
        vfilter = vfilter,
        verbose = FALSE
      )
      av::av_log_level(old_log_level)

      bmp_files <- sort(list.files(tmp_av_dir, pattern = bmp_pat, full.names = TRUE))
      if (!is.null(n) && length(bmp_files) > n) {
        bmp_files <- bmp_files[seq_len(n)]
      }

      if (length(bmp_files) > 0L) {
        jpg_names <- sprintf(paste0("%s%0", digits, "d.jpg"), prefix, seq_along(bmp_files))
        av_files <- file.path(dir_full, jpg_names)
        ok_vec <- convert_frames_to_jpg_cpp(
          src_files = bmp_files,
          dst_files = av_files,
          quality = as.integer(quality),
          remove_src = TRUE
        )
        av_files <- av_files[ok_vec]
      } else {
        av_files <- character(0)
      }
      unlink(tmp_av_dir, recursive = TRUE, force = TRUE)
    } else {
      out_pattern <- file.path(dir_full, paste0(prefix, "%0", digits, "d.png"))
      av::av_encode_video(
        input = video_full,
        output = out_pattern,
        framerate = framerate,
        codec = "png",
        vfilter = vfilter,
        verbose = FALSE
      )
      av::av_log_level(old_log_level)

      av_files <- sort(list.files(dir_full, pattern = pat, full.names = TRUE))
      if (!is.null(n) && length(av_files) > n) {
        excess <- av_files[(n + 1):length(av_files)]
        unlink(excess)
        av_files <- av_files[seq_len(n)]
      }
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

  if (!is.null(pb)) {
    n_done <- nrow(res)
    tryCatch(cli::cli_progress_update(id = pb, set = n_done, total = max(1L, n_done), force = TRUE), error = function(e) NULL)
    tryCatch(cli::cli_progress_done(id = pb), error = function(e) NULL)
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

    suppressWarnings(av::av_video_images(video_full, destdir = tmp_dir, format = "png"))
    img_files <- list.files(tmp_dir, pattern = "\\.png$", full.names = TRUE)

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

.resolve_count_line <- function(count_line, width, height, orientation = "vertical") {
  if (is.null(count_line)) return(NULL)
  count_line <- as.numeric(count_line)
  if (length(count_line) == 1L) {
    orientation <- match.arg(tolower(orientation), c("vertical", "horizontal"))
    if (orientation == "vertical") {
      px <- if (count_line <= 1.0) round(count_line * width) else round(count_line)
      c(as.numeric(px), 0.0, as.numeric(px), as.numeric(height))
    } else {
      py <- if (count_line <= 1.0) round(count_line * height) else round(count_line)
      c(0.0, as.numeric(py), as.numeric(width), as.numeric(py))
    }
  } else if (length(count_line) >= 4L) {
    lx1 <- if (count_line[1] <= 1.0) round(count_line[1] * width) else round(count_line[1])
    ly1 <- if (count_line[2] <= 1.0) round(count_line[2] * height) else round(count_line[2])
    lx2 <- if (count_line[3] <= 1.0) round(count_line[3] * width) else round(count_line[3])
    ly2 <- if (count_line[4] <= 1.0) round(count_line[4] * height) else round(count_line[4])
    c(as.numeric(lx1), as.numeric(ly1), as.numeric(lx2), as.numeric(ly2))
  } else {
    cli::cli_abort("{.arg count_line} must be a single numeric value (0 to 1) or a vector of 4 coordinates c(x1, y1, x2, y2).")
  }
}

.resolve_roi <- function(roi, width, height, orientation = "vertical") {
  if (is.null(roi)) return(NULL)
  roi <- as.numeric(roi)
  if (length(roi) == 2L) {
    orientation <- match.arg(tolower(orientation), c("vertical", "horizontal"))
    r_min <- min(roi)
    r_max <- max(roi)
    if (orientation == "vertical") {
      rx1 <- if (r_min <= 1.0) round(r_min * width) else round(r_min)
      ry1 <- 0L
      rx2 <- if (r_max <= 1.0) round(r_max * width) else round(r_max)
      ry2 <- as.integer(height)
    } else {
      rx1 <- 0L
      ry1 <- if (r_min <= 1.0) round(r_min * height) else round(r_min)
      rx2 <- as.integer(width)
      ry2 <- if (r_max <= 1.0) round(r_max * height) else round(r_max)
    }
  } else if (length(roi) >= 4L) {
    rx1 <- if (roi[1] <= 1.0) round(roi[1] * width) else round(roi[1])
    ry1 <- if (roi[2] <= 1.0) round(roi[2] * height) else round(roi[2])
    rx2 <- if (roi[3] <= 1.0) round(roi[3] * width) else round(roi[3])
    ry2 <- if (roi[4] <= 1.0) round(roi[4] * height) else round(roi[4])
  } else {
    cli::cli_abort("{.arg roi} must be a numeric vector of length 2 c(min, max) or length 4 c(xmin, ymin, xmax, ymax).")
  }
  rx1 <- max(0L, min(width - 2L, as.integer(rx1)))
  ry1 <- max(0L, min(height - 2L, as.integer(ry1)))
  rx2 <- max(rx1 + 2L, min(width, as.integer(rx2)))
  ry2 <- max(ry1 + 2L, min(height, as.integer(ry2)))
  c(rx1, ry1, rx2, ry2)
}

# Factory to create a moving average smoother for telemetry metrics
# - Internally updated on EVERY frame using a rolling window of length `smoot_stat`
# - Displayed values are latched and only updated every `smoot_stat` frames (and on frame 1)
.make_stats_smoother <- function(smoot_stat = 30L) {
  win <- as.integer(smoot_stat[1])
  if (is.na(win) || win < 1L) win <- 1L

  # Buffer for rolling window of length win
  history <- list()
  displayed_vals <- list()
  last_valid <- list()
  frame_idx <- 0L

  # Metrics that represent object properties and should retain their last valid value
  # when no objects are present in the zone/frame, instead of plunging to 0
  persist_keys <- c("flow", "vel_mean", "vel_jit", "size_mean", "size_cv", "conf_mean", "conf_sd", "conf_min")

  function(metrics_list) {
    if (!is.list(metrics_list) || length(metrics_list) == 0L) {
      return(metrics_list)
    }
    if (win <= 1L) {
      return(metrics_list)
    }

    frame_idx <<- frame_idx + 1L

    # 1. Update rolling history for every metric on every frame
    for (nm in names(metrics_list)) {
      v <- metrics_list[[nm]]
      if (is.numeric(v) && length(v) == 1L && !is.na(v)) {
        cur <- c(history[[nm]], v)
        if (length(cur) > win) {
          cur <- cur[(length(cur) - win + 1L):length(cur)]
        }
        history[[nm]] <<- cur
        if (nm %in% persist_keys && v > 1e-4) {
          last_valid[[nm]] <<- v
        }
      }
    }

    # 2. On frame 1, populate displayed_vals immediately so HUD is not blank/zero
    if (length(displayed_vals) == 0L) {
      for (nm in names(metrics_list)) {
        displayed_vals[[nm]] <<- metrics_list[[nm]]
      }
    }

    # 3. Only update displayed values every `win` frames
    if (frame_idx == 1L || (frame_idx %% win == 0L)) {
      for (nm in names(history)) {
        vals <- history[[nm]]
        if (length(vals) > 0L) {
          if (nm == "conf_min") {
            pos_vals <- vals[vals > 1e-4]
            if (length(pos_vals) > 0L) {
              m <- min(pos_vals, na.rm = TRUE)
              last_valid[[nm]] <<- m
              displayed_vals[[nm]] <<- m
            } else if (!is.null(last_valid[[nm]])) {
              displayed_vals[[nm]] <<- last_valid[[nm]]
            } else {
              displayed_vals[[nm]] <<- 0.0
            }
          } else if (nm %in% persist_keys) {
            pos_vals <- vals[vals > 1e-4]
            if (length(pos_vals) > 0L) {
              m <- mean(pos_vals, na.rm = TRUE)
              last_valid[[nm]] <<- m
              displayed_vals[[nm]] <<- m
            } else if (!is.null(last_valid[[nm]])) {
              displayed_vals[[nm]] <<- last_valid[[nm]]
            } else {
              displayed_vals[[nm]] <<- 0.0
            }
          } else {
            displayed_vals[[nm]] <<- mean(vals, na.rm = TRUE)
          }
        }
      }
    }

    # 4. Return the latched displayed values
    res <- metrics_list
    for (nm in names(displayed_vals)) {
      res[[nm]] <- displayed_vals[[nm]]
    }
    res
  }
}

#' Real-Time Video Object Detection with YOLO (ONNX Runtime)
#'
#' Performs real-time object detection and tracking on video files (.mp4, .avi,
#' .mov) or live camera feeds (webcam) using YOLO models via native C++ ONNX
#' Runtime. When `return_data = TRUE`, it computes an extensive ecological and
#' temporal community profile (richness, abundance, Shannon diversity, Simpson
#' index, Pielou evenness, species accumulation curves, and class persistence
#' rates).
#'
#' @details
#' \subsection{Ecological and Analytical Framework}{
#' In automated phenotyping, precision agriculture, ecological monitoring, and industrial inspection,
#' a video stream represents a continuous temporal transect composed of \eqn{F} sequential frames
#' recorded at a frame rate of \eqn{\text{FPS}}{FPS} frames per second. In each frame \eqn{t \in \{1, \dots, F\}},
#' deep neural network inference detects candidate bounding boxes belonging to \eqn{S_t} recognized classes.
#'
#' By aggregating detections over space and time, `video_detect_dl` translates raw computer vision
#' detections into standardized ecological metrics, quantifying instantaneous community state,
#' temporal persistence, class dominance, and cumulative discovery rates.
#' }
#'
#' \subsection{Abundance Metrics and Density Dynamics}{
#' Abundance characterizes the numerical density of detected objects:
#' \itemize{
#'   \item \strong{Instantaneous Abundance (\eqn{N_t})}: The total count of object instances detected in frame \eqn{t}:
#'     \deqn{N_t = \sum_{c=1}^S n_{t, c}}{N_t = \sum_{c=1}^S n_{t, c}}
#'     where \eqn{n_{t, c}} is the count of detected objects belonging to class \eqn{c} in frame \eqn{t}, and \eqn{S} is the total number of candidate classes.
#'   \item \strong{Total Video Abundance (\eqn{N})}: The cumulative sum of all detections identified across the processed footage:
#'     \deqn{N = \sum_{t=1}^F N_t = \sum_{c=1}^S n_c}{N = \sum_{t=1}^F N_t = \sum_{c=1}^S n_c}
#'     where \eqn{n_c = \sum_{t=1}^F n_{t, c}} is the total detection count for class \eqn{c}.
#'   \item \strong{Average Abundance per Frame (\eqn{\bar{N}})}: Quantifies the baseline object density per frame:
#'     \deqn{\bar{N} = \frac{N}{F}}{mean_N = N / F}
#'   \item \strong{Peak Abundance (\eqn{N_{\text{peak}}})}: The maximum instantaneous count observed in any single frame:
#'     \deqn{N_{\text{peak}} = \max_{1 \le t \le F} N_t}{N_peak = max(N_t)}
#'     The exact frame index \eqn{t_{\text{peak}} = \operatorname{argmax}_t N_t} and corresponding timestamp (in seconds)
#'     \eqn{T_{\text{peak}} = (t_{\text{peak}} - 1) / \text{FPS}} identify the moment of maximum congestion or cluster density.
#' }
#' }
#'
#' \subsection{Richness and Species Accumulation Dynamics}{
#' Richness quantifies class variety and categorical complexity over time:
#' \itemize{
#'   \item \strong{Instantaneous Richness (\eqn{S_t})}: The number of distinct classes co-occurring simultaneously in frame \eqn{t}:
#'     \deqn{S_t = \sum_{c=1}^S \mathbb{I}(n_{t, c} > 0)}{S_t = \sum_{c=1}^S I(n_{t, c} > 0)}
#'     where \eqn{\mathbb{I}(\cdot)} is the indicator function equal to 1 if class \eqn{c} has at least one detection in frame \eqn{t}, and 0 otherwise.
#'   \item \strong{Total Video Richness (\eqn{S})}: The total count of unique classes recognized across the entire video:
#'     \deqn{S = |\{c : n_c > 0\}|}{S = count of unique classes detected}
#'   \item \strong{Cumulative Richness / Species Accumulation Curve (Collector's Curve)}:
#'     \deqn{S_{\text{cum}}(t) = |\{c : \exists t' \le t \text{ such that } n_{t', c} > 0\}|}{S_cum(t) = count of distinct classes discovered up to frame t}
#'     The Collector's Curve plots the progressive discovery of novel classes over time. A plateau indicates that sampling effort
#'     (video duration) is sufficient to capture the full assemblage of classes present in the environment.
#'   \item \strong{Cumulative Detections (\eqn{N_{\text{cum}}(t)})}:
#'     \deqn{N_{\text{cum}}(t) = \sum_{t'=1}^t N_{t'}}{N_cum(t) = \sum_{t'=1}^t N_{t'}}
#'     Tracking the cumulative volume of detections provides insights into detection throughput and arrival rates.
#' }
#' }
#'
#' \subsection{Per-Class Ecological Metrics and Constancy}{
#' For each class \eqn{c}, `video_detect_dl` extracts detailed persistence and distribution metrics:
#' \itemize{
#'   \item \strong{Total Detections (\eqn{n_c})}: Total bounding box instances identified for class \eqn{c}.
#'   \item \strong{Relative Abundance (\eqn{p_c})}: The percentage contribution of class \eqn{c} to total video detections:
#'     \deqn{p_c = \frac{n_c}{N} \times 100\%}{p_c = (n_c / N) * 100\%}
#'   \item \strong{Frames Present (\eqn{F_c})}: The number of frames in which class \eqn{c} was detected at least once:
#'     \deqn{F_c = \sum_{t=1}^F \mathbb{I}(n_{t, c} > 0)}{F_c = \sum_{t=1}^F I(n_{t, c} > 0)}
#'   \item \strong{Occurrence Rate / Constancy (\eqn{O_c})}: The temporal persistence of class \eqn{c} across the video footage:
#'     \deqn{O_c = \frac{F_c}{F} \times 100\%}{O_c = (F_c / F) * 100\%}
#'     In ecological literature, species with \eqn{O_c \ge 50\%} are classified as constant, \eqn{25\% \le O_c < 50\%} as accessory, and \eqn{O_c < 25\%} as accidental.
#'   \item \strong{Mean Abundance When Present}: The expected number of individuals per frame given that the class is present:
#'     \deqn{\bar{n}_{c, \text{present}} = \frac{n_c}{\max(1, F_c)}}{mean_n_present = n_c / F_c}
#'   \item \strong{Mean Abundance Overall}: The unconditional expected number of individuals per frame across all processed frames:
#'     \deqn{\bar{n}_{c, \text{overall}} = \frac{n_c}{F}}{mean_n_overall = n_c / F}
#'   \item \strong{Maximum per Frame (\eqn{\max_t n_{t, c}})}: Peak instantaneous count of class \eqn{c}.
#'   \item \strong{Mean Confidence (\eqn{\bar{s}_c})}: The arithmetic mean of detection confidence scores for class \eqn{c}:
#'     \deqn{\bar{s}_c = \frac{1}{n_c} \sum_{i=1}^{n_c} \text{score}_i}{mean_conf = (1 / n_c) * \sum score_i}
#'   \item \strong{Temporal Span}: First observed timestamp (\eqn{t_{\text{first}, c}}) and last observed timestamp (\eqn{t_{\text{last}, c}}) in seconds.
#' }
#' }
#'
#' \subsection{Diversity, Dominance, and Evenness Indices}{
#' Summary diversity indices synthesize class variety and equitability into scalar coefficients:
#' \itemize{
#'   \item \strong{Shannon-Wiener Diversity Index (\eqn{H'})}:
#'     \deqn{H' = -\sum_{c=1}^S p_c \ln(p_c)}{H' = - \sum_{c=1}^S p_c * ln(p_c)}
#'     where \eqn{p_c = n_c / N}. \eqn{H'} quantifies the uncertainty in predicting the class identity of a randomly sampled detection.
#'     Higher values indicate higher diversity and a more equitable distribution across classes.
#'     In addition to the global video-level index, \eqn{H'_t} is computed dynamically for each frame in the `$temporal` data frame.
#'   \item \strong{Simpson's Diversity Index (\eqn{1 - D}) and Dominance (\eqn{D})}:
#'     \deqn{D = \sum_{c=1}^S p_c^2, \qquad 1 - D = 1 - \sum_{c=1}^S p_c^2}{D = \sum_{c=1}^S p_c^2, 1 - D = 1 - \sum_{c=1}^S p_c^2}
#'     \eqn{D} represents the probability that two randomly chosen detections belong to the same class (dominance),
#'     while \eqn{1 - D} (Gini-Simpson index) represents the probability that they belong to different classes.
#'     Values range from 0 (monoculture / single class dominating entirely) to \eqn{1 - 1/S} (maximal diversity).
#'   \item \strong{Pielou's Evenness Index (\eqn{J'})}:
#'     \deqn{J' = \frac{H'}{\ln(S)}}{J' = H' / ln(S)}
#'     for \eqn{S > 1} (with \eqn{J' = 1.0} if \eqn{S = 1}, and \eqn{J' = 0.0} if \eqn{S = 0}).
#'     \eqn{J'} normalizes Shannon diversity against its theoretical maximum (\eqn{H'_{\text{max}} = \ln(S)}),
#'     providing an equitability index bounded between 0 and 1. A value of 1.0 indicates that all classes are detected with equal frequency.
#'   \item \strong{Berger-Parker Dominance Index (\eqn{d})}:
#'     \deqn{d = \max_{1 \le c \le S} (p_c)}{d = max(p_c)}
#'     Quantifies the proportional dominance of the most abundant single class in the assemblage.
#' }
#' }
#'
#' \subsection{Multi-Object Tracking and Unique Entity Counting}{
#' When `track = TRUE`, a spatial-temporal Kalman filter tracking algorithm (SORT/ByteTrack)
#' associates detections across consecutive frames using Intersection over Union (IoU) cost matrices:
#' \itemize{
#'   \item \strong{Unique Objects (\eqn{U})}: Persistent IDs (\verb{#1}, \verb{#2}, ...) track physical entities over time,
#'     preventing double-counting when calculating true entity abundance:
#'     \deqn{U = |\{ \text{distinct persistent IDs} > 0 \}|}{U = count of distinct persistent IDs}
#'     For each class \eqn{c}, \eqn{U_c} records the unique tracked individuals assigned to that class.
#'   \item \strong{Cumulative Unique Objects Curve (\eqn{U_{\text{cum}}(t)})}:
#'     Tracks the progressive entry of novel unique physical entities into the monitored frame.
#'   \item \strong{Virtual Tripwire Counting (\code{count_line})}:
#'     When a virtual line is provided as \code{c(x1, y1, x2, y2)} in normalized coordinates \eqn{[0, 1]},
#'     the trajectory segments between consecutive centroids \eqn{\mathbf{p}_{t-1}} and \eqn{\mathbf{p}_t} of tracked objects
#'     are evaluated for intersection with the virtual line segment.
#'     Crossing events are registered precisely when the signed 2D cross-product changes orientation,
#'     cataloguing the event frame, timestamp, persistent ID, and class label in `$crossings`.
#' }
#' }
#'
#' \subsection{Visualization and S3 Methods}{
#' The returned object has class `c("pliman_video_detect", "list")` and provides built-in S3 methods:
#' \itemize{
#'   \item `plot(x, type = ...)`: Renders clean, publication-ready graphics. Supported `type` options:
#'     \itemize{
#'       \item `"temporal"` (default): Two-panel stacked time series displaying total abundance \eqn{N(t)} and richness \eqn{S(t)}.
#'       \item `"abundance"`: Multi-class abundance trajectories over time.
#'       \item `"richness"`: Step-function curve of class richness over time.
#'       \item `"cumulative"`: Two-panel plot of the Species Accumulation Curve (Collector's Curve) and cumulative detection volume.
#'       \item `"class"`: Horizontal barplot displaying total detections and relative abundance percentages for each class.
#'       \item `"all"`: 2x2 multi-panel overview combining abundance, richness, class distribution, and the Collector's Curve.
#'     }
#'   \item `summary(x)`: Formats and displays an overview table with video metadata, ecological metrics, and class distribution.
#'   \item `print(x)`: Displays an informative console summary with key video metrics and class rankings.
#'   \item `as.data.frame(x)`, `head(x)`, `tail(x)`, `dim(x)`, `nrow(x)`, `ncol(x)`, `x$column`, `x[i, j]`:
#'     Provide transparent backward compatibility with data frame operations on raw detections.
#' }
#' }
#'
#' @param video Input video source. Either a camera index (e.g. `0`, `1`), a camera name
#'   or substring (e.g. `"EMEET"`, `"HD User Facing"`), `"camera"` (default camera `0`),
#'   or a path to a video file (e.g., `"field_drone.mp4"`).
#' @param model Model name or path to a custom ONNX file. Defaults to `"yolo26n"`.
#'   Pre-trained options include `"yolo26n"`, `"yolo26s"`, `"yolov8n"`, `"yolov8s"`,
#'   `"yolov11n"`, etc., or a full path to a custom user-trained `.onnx` model.
#' @param conf_threshold Minimum confidence threshold for candidate detections (default 0.25).
#' @param iou_threshold IoU threshold for Non-Maximum Suppression (default 0.45).
#' @param labels Optional character vector of class names. If `NULL` (default), labels are
#'   automatically detected from model metadata or fallback to standard classes.
#' @param classes Optional character vector of specific class names to filter and retain.
#'   If specified, only detections matching these classes will be kept.
#' @param detect Logical. Whether to run object detection model inference. If `FALSE`,
#'   the camera/video stream runs in pass-through mode without neural network inference, allowing
#'   fast real-time streaming and recording. Defaults to `TRUE`.
#' @param output Optional output file path (e.g. `"annotated_video.mp4"`) to save the processed video.
#' @param save_video Logical. Whether to save the processed video. Automatically `TRUE` if `output` is specified.
#' @param return_data Logical. Whether to compute and return the detection dataset and ecological summary statistics (default `TRUE`). If `FALSE`, returns `invisible(NULL)`.
#' @param track Logical. Whether to enable multi-object tracking and assign persistent IDs (`#1`, `#2`, etc.) across frames. Defaults to `FALSE`.
#' @param min_frames Minimum number of frames an object must persist across the video to be retained in tracking statistics. Defaults to `10L`, discarding fleeting, transient noise artifacts (e.g. spurious false positives appearing in 10 or fewer frames). Set to `0` or `NULL` to retain all detections.
#' @param interpolate_box Logical. Whether to temporally interpolate bounding box positions for tracked objects across briefly lost or occluded frames using velocity motion vectors (default `FALSE`).
#' @param interpolate_mask Logical. Whether to interpolate instance segmentation masks alongside interpolated bounding boxes (default `FALSE`).
#' @param trail Logical. Whether to draw historical trajectory motion trail lines behind tracked objects when `track = TRUE`. Defaults to `FALSE`.
#' @param cex Size factor for label text (default `1.0`). Use smaller values (e.g., `cex = 0.2` or `0.5`)
#'   to reduce text size for dense scenes or many objects.
#' @param pad Padding multiplier controlling the area/background badge size around label text (default `1.0`).
#'   Can be a single value or a length-2 vector `c(pad_x, pad_y)`. Use `pad = 0` or smaller values for a compact area.
#' @param badge Logical. Whether to draw a solid background badge behind label text (default `TRUE`).
#'   If `FALSE`, only the text is drawn directly.
#' @param text_size Numeric font scale for bounding box text labels and/or telemetry HUD (alias for `cex`, default `1.0`).
#'   Can be a single numeric value (applied to both objects and telemetry HUD) or a numeric vector of length 2
#'   `c(objects, telemetry)` to scale them independently (e.g. `text_size = c(1, 3)`).
#' @param text_size_telemetry Numeric font scale specifically for the telemetry HUD overlay (e.g., `text_size_telemetry = 3`).
#'   If provided, overrides the telemetry component of `text_size`. Set to `0` to hide telemetry text.
#' @param text_size_objects Numeric font scale specifically for detected object bounding box labels and badges (e.g., `text_size_objects = 1`).
#'   If provided, overrides the object component of `text_size`.
#' @param text_size_roi Numeric font scale specifically for Region of Interest (ROI).
#' @param show_text Logical. Whether to display text badges above bounding boxes (default `TRUE`).
#' @param show_conf Logical. Whether to display confidence scores in badges (default `TRUE`).
#' @param show_class Logical. Whether to display class label names in badges (default `TRUE`).
#' @param show_id Logical. Whether to display persistent tracking IDs (`#1`, `#2`) in badges when `track = TRUE` (default `FALSE`).
#' @param count_line Optional virtual counting tripwire line. Can be specified as a single numeric
#'   value between `0` and `1` (indicating relative position along the axis, e.g. `0.5` for center)
#'   combined with `count_line_orientation`, or legacy 4-element coordinates `c(x1, y1, x2, y2)`.
#'   Objects crossing this virtual line are counted.
#' @param count_line_orientation Orientation of the counting line when `count_line` is a single scalar:
#'   either `"vertical"` (default) or `"horizontal"`. Can also be passed as `count_line_dir`.
#' @param count_line_label Optional label displayed along the counting line (default `NULL`, no text).
#' @param count_direction Character string specifying which crossing direction(s) trigger counting
#'   when `count_line` is set: `"both"` (default), `"forward"`, or `"backward"`.
#' @param roi Optional Region of Interest coordinates. Can be specified as a numeric vector of
#'   length 2 `c(min, max)` in relative coordinates (`0` to `1`) combined with `roi_orientation`
#'   (e.g., `roi = c(0.2, 0.8), roi_orientation = "vertical"`), or legacy 4-element coordinates `c(xmin, ymin, xmax, ymax)`.
#' @param roi_orientation Orientation of the ROI band when `roi` has length 2:
#'   either `"vertical"` (default) or `"horizontal"`. Can also be passed as `roi_dir`.
#' @param roi_mode Either `"crop"` (infer only inside ROI) or `"filter"` (detect everywhere, retain only detections inside ROI).
#' @param roi_label Label displayed on the ROI box (default `"ROI / MONITORED ZONE"`).
#' @param hide_outside_roi Logical. Whether to remove / suppress bounding boxes and visual annotations for objects outside the ROI zone. Defaults to `TRUE` when `roi` is provided.
#' @param hide_counted Logical. Whether to remove / suppress bounding boxes for objects that have already crossed the virtual counting line. Defaults to `TRUE` when counting / tracking is active.
#' @param flash_counted Logical. Whether to show a brief neon green confirmation flash (~12 frames) when an object crosses the counting line before removing its bounding box. Defaults to `TRUE`.
#' @param fullscreen Logical. Whether to display the video feed in full screen mode (default `FALSE`).
#' @param window_size Optional numeric vector `c(width, height)` for the display window.
#' @param window_scale Scaling factor for the display window (default 1.0).
#' @param resolution Optional camera capture resolution as a numeric vector `c(width, height)`,
#'   e.g. `c(1920, 1080)` or `c(1280, 720)`. Default is `NULL` (uses native sensor resolution).
#' @param crop Optional frame cropping. Set to `"1:1"`, `"square"`, or `TRUE` for an automatic 1:1 centered square crop (e.g. 1080x1080 or 720x720), or `c(width, height)` for a custom centered rectangular crop. Default is `NULL`.
#' @param backend Camera and display backend: `"native"` (default, uses native Windows Media Foundation / Win32 GDI with zero external dependencies).
#' @param skeleton Logical. Whether to render the 17-keypoint anatomical skeleton when running pose estimation models (e.g. `"yolo26n-pose"`). Defaults to `TRUE`.
#' @param kpt_threshold Minimum confidence threshold for rendering individual keypoints and skeleton limbs (default 0.3).
#' @param kpt_radius Integer radius of keypoint joint markers in pixels (default 4).
#' @param bbox,show_box Logical. Whether to render bounding boxes around detected objects or people (default `TRUE`).
#' @param show_point Logical. Whether to render centroid markers/points for detected objects (default `FALSE`).
#' @param point_radius Integer radius of centroid markers in pixels (default 4).
#' @param point_col Optional color for centroid markers. If `NULL` (default), matches the object's class/instance palette.
#' @param mask Logical. Whether to render semi-transparent instance segmentation masks when using segmentation models (default `TRUE`).
#' @param alpha Opacity/transparency of instance segmentation masks (default `0.4`, range `0` to `1`).
#' @param saturation Numeric color saturation multiplier applied in-place via native C++ OpenMP fixed-point arithmetic (default `1.0`). Values `> 1.0` (e.g., `1.3` or `1.5`) boost vividness on flat/washed-out `.raw` or log-profile footage with zero FPS penalty.
#' @param contrast Numeric contrast multiplier around midpoint luma applied in-place via native C++ OpenMP fixed-point arithmetic (default `1.0`). Values `> 1.0` (e.g., `1.2` or `1.4`) enhance dynamic range on flat footage with zero FPS penalty.
#' @param heatmap Logical. Whether to render a cumulative spatial density heatmap overlay for detected objects. Defaults to `FALSE`.
#' @param heatmap_alpha Numeric opacity in `[0, 1]` for the spatial density heatmap overlay. Defaults to `0.40`.
#' @param heatmap_decay Numeric decay factor in `(0, 1]` controlling how quickly older detections fade from the heatmap per frame. Defaults to `0.98`.
#' @param sahi Logical. Whether to enable SAHI (Slicing Aided Hyper Inference) tiled detection on each video frame for small/distant objects. Defaults to `FALSE`.
#' @param slice_height Optional integer height of SAHI sliding window slices (in pixels). Defaults to `NULL` (auto-sized).
#' @param slice_width Optional integer width of SAHI sliding window slices (in pixels). Defaults to `NULL` (auto-sized).
#' @param overlap_ratio Numeric overlap fraction in `[0, 0.5)` between adjacent SAHI slices. Defaults to `0.2`.
#' @param save_crops Logical. Whether to save cropped images of detected objects to disk. Defaults to `FALSE`.
#' @param crops_dir Optional character directory path where cropped object images are saved when `save_crops = TRUE`. Defaults to `NULL`.
#' @param start_time Optional numeric start timestamp (in seconds) from which to begin processing a video file. Defaults to `NULL`.
#' @param end_time Optional numeric end timestamp (in seconds) at which to stop processing a video file. Defaults to `NULL`.
#' @param show Logical. Whether to display the video feed with bounding boxes in real time.
#'   Default is `interactive()`.
#' @param show_fps Logical. Whether to overlay a real-time FPS and detection counter badge (default `TRUE`).
#' @param telemetry Telemetry HUD display mode: `"detailed"` (default, glassmorphic diagnostic HUD card with real-time FPS, latency breakdown (inference, tracking, GUI), object counts, confidence statistics, flow rate, motion velocity, ROI occupancy, and 20-frame sparkline trend), `"simple"` (compact FPS, in-zone, and count badges), or `"none"` (disables HUD overlay).
#' @param smoot_stat Window size (in frames) for moving average of telemetry statistics (default `20`).
#'   Internally, the moving average of metrics (SPEED, LATENCY, INF, TRK, GUI, CONF, FLOW, MOTION, SIZE, ROI OCC)
#'   is updated on every single frame, but displayed values on the telemetry HUD are refreshed every `smoot_stat` frames.
#'   Frame number, in-zone count, counted total, and active tracks are displayed instantaneously on every frame.
#'   Supported aliases: `smooth_stat`, `smooth_stats`, `smoot_stats`.
#' @param min_frames Minimum number of consecutive/persistent frames an object must be tracked before rendering its bounding box or point (default `10L`). Suppresses transient false positives and noise from flickering on screen.
#' @param as_job Logical. If `TRUE` and running inside RStudio, dispatches the video processing task to an asynchronous RStudio Background Job, keeping the R console free. Default is `FALSE`.
#' @param hud_pos Screen position of the telemetry HUD card: either a preset string
#'   (`"top-right"` (default, avoids top-left ROI zone), `"top-left"`, `"bottom-left"`,
#'   `"bottom-right"`, `"top"`, or `"bottom"`), or a numeric vector of length 2 `c(x, y)`
#'   specifying relative coordinates in range 0 to 1 across frame width and height (e.g., `c(0.2, 0.5)`).
#' @param hud_layout Layout style of the telemetry HUD: `"vertical"` (default modern multi-row glassmorphic card)
#'   or `"horizontal"` (sleek capsule pill bar).
#' @param infer_every Integer stride for neural network inference (default `1L`, running inference on every frame).
#'   Setting e.g. `infer_every = 2L` (or using alias `stride = 2` / `skip_frames = 1`) runs the heavy neural network
#'   on odd frames and smoothly propagates tracks and detections on intermediate frames, cutting compute time in half
#'   and maximizing real-time video display smoothness up to the camera's physical frame rate.
#' @param max_frames Optional maximum number of frames to process. Default is `NULL` (all frames).
#' @param fps Target framerate for video output. If `NULL` (default), matches the source video or camera.
#' @param col Optional color or palette for bounding boxes. If `NULL`, distinct colors are assigned.
#' @param lwd Line width for bounding boxes (default 2).
#' @param rainbow Logical. If `TRUE`, assigns distinct colors per detected instance.
#' @param engine Execution engine: `"gpu"` (DirectML, default) or `"cpu"`.
#' @param device_id GPU device ID (default -1 for auto).
#' @param threads Number of CPU threads (default 0 for auto-tuning).
#' @param dir Directory where pre-trained models are stored (default `pliman_model_dir()`).
#' @param imgsz Target input image resolution for YOLO inference (e.g., `640`, `960`, or `1024`).
#'   If `NULL` (default), inferred automatically from the model or defaults to `640`.
#' @param verbose Logical. Whether to show progress messages.
#' @param ... Additional arguments passed to graphics plotting functions or aliases
#'   (`camera`, `cam`, `device`, `infer`, `inference`, `save`, `return_df`, `return_detections`).
#' @return If `return_data = TRUE`, an object of class `c("pliman_video_detect", "list")` containing:
#' \itemize{
#'   \item `summary`: A `data.frame` of video-level metadata and ecological diversity metrics:
#'     `video_source`, `model`, `engine`, `total_frames`, `fps`, `duration_sec`,
#'     `total_detections` (\eqn{N}), `total_richness` (\eqn{S}), `unique_objects` (\eqn{U}),
#'     `avg_abundance_frame` (\eqn{\bar{N}}), `max_abundance_frame` (\eqn{N_{\text{peak}}}),
#'     `peak_frame` (\eqn{t_{\text{peak}}}), `peak_time_sec` (\eqn{T_{\text{peak}}}),
#'     `shannon_index` (\eqn{H'}), `simpson_index` (\eqn{1-D}), `pielou_evenness` (\eqn{J'}),
#'     and optionally `total_counted` (when `count_line` is used).
#'   \item `by_class`: A `data.frame` sorted by total detections detailing per-class metrics:
#'     `class`, `total_detections` (\eqn{n_c}), `relative_abundance` (\eqn{p_c}, \%),
#'     `unique_objects` (\eqn{U_c}), `n_frames_present` (\eqn{F_c}), `occurrence_rate` (\eqn{O_c}, \%),
#'     `mean_abundance_when_present` (\eqn{\bar{n}_{c, \text{present}}}), `mean_abundance_overall` (\eqn{\bar{n}_{c, \text{overall}}}),
#'     `max_per_frame`, `mean_confidence` (\eqn{\bar{s}_c}), `first_seen_sec` (\eqn{t_{\text{first}, c}}),
#'     `last_seen_sec` (\eqn{t_{\text{last}, c}}), and optionally `total_counted`.
#'   \item `unique_objects`: A `data.frame` aggregating every unique tracked object across all observed frames:
#'     `id`, `label`, `n_frames`, `first_frame`, `last_frame`, `duration_sec`, `mean_score`,
#'     and (when using instance segmentation models) detailed morphological shape features:
#'     `area` (in pixels), `area_sd`, `perimeter`, `radius_mean`, `length`, `width`, `circularity`, `eccentricity`, and `asp_ratio`.
#'   \item `temporal`: A `data.frame` of frame-by-frame time series at resolution \eqn{t \in \{1, \dots, F\}}:
#'     `frame`, `timestamp` (in seconds), instantaneous `abundance` (\eqn{N_t}), instantaneous `richness` (\eqn{S_t}),
#'     instantaneous `shannon` diversity (\eqn{H'_t}), and individual count columns for every detected class.
#'   \item `cumulative`: A `data.frame` of cumulative progression:
#'     `frame`, `timestamp`, `cumulative_detections` (\eqn{N_{\text{cum}}(t)}),
#'     `cumulative_richness` (Species Accumulation Curve / Collector's Curve, \eqn{S_{\text{cum}}(t)}),
#'     and optionally `cumulative_unique_objects` (\eqn{U_{\text{cum}}(t)}).
#'   \item `diversity`: A `list` containing global ecological diversity indices:
#'     `richness_S` (\eqn{S}), `total_abundance_N` (\eqn{N}), `shannon_H` (\eqn{H'}),
#'     `simpson_1_minus_D` (\eqn{1 - D}), `simpson_D` (\eqn{D}), `pielou_J` (\eqn{J'}),
#'     `berger_parker_d` (\eqn{d}), `mean_richness_per_frame`, and `mean_abundance_per_frame`.
#'   \item `detections`: A `data.frame` of all raw candidate detections across all frames:
#'     `frame`, `timestamp`, `id`, `xmin`, `ymin`, `xmax`, `ymax`, `label`, `score`, and `class_id`
#'     (plus shape features when running segmentation models).
#'   \item `crossings`: (Optional) A `data.frame` of virtual tripwire crossing events when `count_line` is specified.
#'   \item `keypoints`: (Optional) A `data.frame` of 17 keypoint anatomical landmarks when running pose estimation models.
#'   \item `roi`, `roi_mode`, `count_line`, `output`: Additional metadata reflecting execution configuration.
#' }
#' If `return_data = FALSE`, returns `invisible(NULL)`.
#' @export
#' @examplesIf interactive()
#'   library(pliman)
#'
#'   # 1. Video detection and ecological diversity analysis
#'   dets <- video_detect_dl(
#'     video = "field_drone.mp4",
#'     model = "yolo26n",
#'     conf_threshold = 0.5,
#'     return_data = TRUE,
#'     output = "annotated_field.mp4"
#'   )
#'
#'   # Summary table of ecological metrics
#'   summary(dets)
#'
#'   # Plot community temporal dynamics (abundance & richness)
#'   plot(dets, type = "temporal")
#'
#'   # Plot Species Accumulation Curve (Collector's Curve)
#'   plot(dets, type = "cumulative")
#'
#'   # Plot multi-class abundance trajectories
#'   plot(dets, type = "abundance")
#'
#'   # Multi-panel 2x2 diagnostic overview
#'   plot(dets, type = "all")
#'
#'   # Access per-class breakdown
#'   dets$by_class
#'
#'   # 2. Multi-object tracking with virtual tripwire line counting
#'   tracked_dets <- video_detect_dl(
#'     video = "conveyor_belt.mp4",
#'     model = "seed_weights.onnx",
#'     track = TRUE,
#'     count_line = c(0.5, 0, 0.5, 1), # vertical counting tripwire
#'     return_data = TRUE
#'   )
#' @export

video_detect_dl <- function(video = 0,
                            model = "yolo26n",
                            conf_threshold = 0.25,
                            iou_threshold = 0.45,
                            labels = NULL,
                            classes = NULL,
                            detect = TRUE,
                            output = NULL,
                            save_video = !is.null(output),
                            return_data = TRUE,
                            track = FALSE,
                            min_frames = 20L,
                            interpolate_box = FALSE,
                            interpolate_mask = FALSE,
                            trail = FALSE,
                            cex = 1.0,
                            pad = 1.0,
                            badge = TRUE,
                            show_text = TRUE,
                            show_conf = TRUE,
                            show_class = TRUE,
                            show_id = FALSE,
                            text_size = 1.0,
                            text_size_telemetry = NULL,
                            text_size_objects = NULL,
                            text_size_roi = NULL,
                            count_line = NULL,
                            count_line_orientation = c("vertical", "horizontal"),
                            count_line_label = NULL,
                            count_direction = c("both", "forward", "backward"),
                            roi = NULL,
                            roi_orientation = c("vertical", "horizontal"),
                            roi_mode = c("crop", "filter"),
                            roi_label = "ROI / MONITORED ZONE",
                            hide_outside_roi = TRUE,
                            hide_counted = TRUE,
                            flash_counted = TRUE,
                            fullscreen = FALSE,
                            window_size = NULL,
                            window_scale = 1.0,
                            resolution = NULL,
                            crop = NULL,
                            backend = c("native", "auto"),
                            skeleton = TRUE,
                            kpt_threshold = 0.3,
                            kpt_radius = 4,
                            bbox = TRUE,
                            show_box = NULL,
                            show_point = FALSE,
                            point_radius = 4,
                            point_col = NULL,
                            mask = TRUE,
                            alpha = 0.4,
                            saturation = 1,
                            contrast = 1,
                            heatmap = FALSE,
                            heatmap_alpha = 0.40,
                            heatmap_decay = 0.98,
                            sahi = FALSE,
                            slice_height = NULL,
                            slice_width = NULL,
                            overlap_ratio = 0.2,
                            save_crops = FALSE,
                            crops_dir = NULL,
                            start_time = NULL,
                            end_time = NULL,
                            show = NULL,
                            show_fps = TRUE,
                            telemetry = c("detailed", "simple", "none"),
                            smoot_stat = 20,
                            hud_pos = c("top-right", "top-left", "bottom-left", "bottom-right", "top", "bottom"),
                            hud_layout = c("vertical", "horizontal"),
                            infer_every = 1L,
                            max_frames = NULL,
                            fps = NULL,
                            col = NULL,
                            lwd = 2,
                            rainbow = TRUE,
                            engine = c("gpu", "cpu"),
                            device_id = -1,
                            threads = 0,
                            dir = pliman_model_dir(),
                            imgsz = NULL,
                            as_job = FALSE,
                            verbose = TRUE,
                            ...) {
  if (isTRUE(as_job)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() && rstudioapi::hasFun("jobRunScript")) {
      call_args <- match.call(expand.dots = TRUE)
      env <- parent.frame()
      eval_args <- as.list(call_args)[-1]
      eval_args <- lapply(eval_args, function(x) eval(x, envir = env))
      eval_args$as_job <- FALSE

      job_dir <- file.path(tempdir(), paste0("pliman_vid_job_", as.integer(stats::runif(1, 10000, 99999))))
      dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
      job_args_rds <- file.path(job_dir, "job_args.rds")
      job_script_path <- file.path(job_dir, "run_job.R")

      job_payload <- list(
        getwd = getwd(),
        args = eval_args
      )
      saveRDS(job_payload, job_args_rds)

      v_label <- if (is.character(video)) basename(video[1]) else paste0("Camera ", video[1])
      m_label <- if (is.character(model)) basename(model[1]) else "YOLO"
      j_name <- paste0("video_detect_dl: ", v_label)

      job_code <- paste0(
        "# Auto-generated RStudio Background Job by pliman::video_detect_dl()\n",
        "library(pliman)\n",
        "cat('================================================================\\n')\n",
        "cat('   pliman video_detect_dl() - RStudio Background Job\\n')\n",
        "cat('================================================================\\n\\n')\n",
        "job_info <- readRDS(", base::deparse(job_args_rds), ")\n",
        "setwd(job_info$getwd)\n",
        "res <- do.call(pliman::video_detect_dl, job_info$args)\n",
        "cat('\\n[OK] Video detection background job finished successfully!\\n')\n"
      )
      writeLines(job_code, job_script_path)

      job_id <- rstudioapi::jobRunScript(
        path = job_script_path,
        name = j_name,
        workingDir = getwd(),
        importEnv = FALSE
      )

      cli::cli_h2("Video Detection (pliman)")
      cli::cli_alert_success("Video processing task has been dispatched to RStudio's {.strong Jobs} tab!")
      cli::cli_ul(c(
        "Job name: {.val {j_name}}",
        "Video: {.val {v_label}}",
        "Model: {.val {m_label}}",
        "Working directory: {.path {getwd()}}"
      ))
      cli::cli_alert_info("Your R console is free! You can track live progress in the {.strong Jobs} tab.")
      return(invisible(job_id))
    } else {
      cli::cli_alert_warning("RStudio Job runner not available. Running synchronously in current session...")
    }
  }

  positions <- c("vertical", "horizontal")
  count_line_orientation <- positions[pmatch(count_line_orientation, positions)]
  roi_orientation <- positions[pmatch(roi_orientation, positions)]
  hud_layout <- positions[pmatch(hud_layout, positions)]

  if (!missing(cex) && missing(text_size)) {
    text_size <- cex
  } else if (!missing(text_size) && missing(cex)) {
    cex <- text_size
  }
  dots <- list(...)
  if (is.null(imgsz)) {
    if ("imgz" %in% names(dots)) imgsz <- dots$imgz
    if ("img_size" %in% names(dots)) imgsz <- dots$img_size
    if ("image_size" %in% names(dots)) imgsz <- dots$image_size
  }
  if ("camera" %in% names(dots)) video <- dots$camera
  if ("cam" %in% names(dots)) video <- dots$cam
  if ("device" %in% names(dots)) video <- dots$device
  if ("square" %in% names(dots) && isTRUE(dots$square)) crop <- "1:1"
  if ("aspect_ratio" %in% names(dots) && dots$aspect_ratio %in% c("1:1", "square")) crop <- "1:1"
  if ("trail" %in% names(dots)) trail <- dots$trail
  if ("draw_trail" %in% names(dots)) trail <- dots$draw_trail
  if ("remove_counted" %in% names(dots)) hide_counted <- isTRUE(dots$remove_counted)
  if ("bbox_counted" %in% names(dots)) hide_counted <- !isTRUE(dots$bbox_counted)
  if ("show_counted" %in% names(dots)) hide_counted <- !isTRUE(dots$show_counted)
  if ("bbox_outside_roi" %in% names(dots)) hide_outside_roi <- !isTRUE(dots$bbox_outside_roi)
  if ("show_outside_roi" %in% names(dots)) hide_outside_roi <- !isTRUE(dots$show_outside_roi)
  if ("remove_outside_roi" %in% names(dots)) hide_outside_roi <- isTRUE(dots$remove_outside_roi)
  if ("track_trail" %in% names(dots)) trail <- dots$track_trail
  if ("cex" %in% names(dots)) { text_size <- dots$cex; cex <- dots$cex }
  if ("text_size" %in% names(dots)) { text_size <- dots$text_size; cex <- dots$text_size }
  if ("font_scale" %in% names(dots)) { text_size <- dots$font_scale; cex <- dots$font_scale }
  if ("label_size" %in% names(dots)) { text_size <- dots$label_size; cex <- dots$label_size }
  if ("cex_scale" %in% names(dots)) { text_size <- dots$cex_scale; cex <- dots$cex_scale }
  for (nm in c("text_size_objects", "text_size_obj", "text_size_boxes", "objects_text_size", "object_text_size", "obj_size", "cex_objects", "cex_obj")) {
    if (nm %in% names(dots)) {
      text_size_objects <- dots[[nm]]
      break
    }
  }
  for (nm in c("text_size_telemetry", "text_size_hud", "telemetry_text_size", "hud_text_size", "telemetry_size", "hud_size", "telemetry_scale", "hud_scale", "cex_telemetry", "cex_hud")) {
    if (nm %in% names(dots)) {
      text_size_telemetry <- dots[[nm]]
      break
    }
  }
  if ("pad" %in% names(dots)) pad <- dots$pad
  if ("label_pad" %in% names(dots)) pad <- dots$label_pad
  if ("badge_pad" %in% names(dots)) pad <- dots$badge_pad
  if ("pad_scale" %in% names(dots)) pad <- dots$pad_scale
  if ("badge" %in% names(dots)) badge <- dots$badge
  if ("show_labels" %in% names(dots)) show_class <- isTRUE(dots$show_labels)
  if ("labels" %in% names(dots) && is.logical(dots$labels)) show_class <- isTRUE(dots$labels)
  if ("show_scores" %in% names(dots)) show_conf <- isTRUE(dots$show_scores)
  if ("scores" %in% names(dots) && is.logical(dots$scores)) show_conf <- isTRUE(dots$scores)
  if ("conf" %in% names(dots) && is.logical(dots$conf)) show_conf <- isTRUE(dots$conf)
  if ("show_text" %in% names(dots)) show_text <- isTRUE(dots$show_text)
  if ("show_id" %in% names(dots)) show_id <- isTRUE(dots$show_id)
  if ("mask" %in% names(dots)) mask <- dots$mask
  if ("masks" %in% names(dots)) mask <- dots$masks
  if ("draw_mask" %in% names(dots)) mask <- dots$draw_mask
  if ("draw_masks" %in% names(dots)) mask <- dots$draw_masks
  if ("alpha" %in% names(dots)) alpha <- dots$alpha
  if ("mask_alpha" %in% names(dots)) alpha <- dots$mask_alpha
  if ("infer" %in% names(dots)) detect <- dots$infer
  if ("inference" %in% names(dots)) detect <- dots$inference
  if ("save" %in% names(dots)) save_video <- dots$save
  if ("return_df" %in% names(dots)) return_data <- dots$return_df
  if ("return_detections" %in% names(dots)) return_data <- dots$return_detections

  detect <- isTRUE(detect)
  save_video <- isTRUE(save_video)
  return_data <- isTRUE(return_data)
  trail <- isTRUE(trail)
  mask <- isTRUE(mask)
  alpha <- as.numeric(alpha)
  if (is.na(alpha) || alpha < 0) alpha <- 0.4
  if (alpha > 1.0) alpha <- 1.0
  src_fps <- if (!is.null(fps) && is.numeric(fps) && fps > 0) as.numeric(fps) else 30.0

  if ("telemetry" %in% names(dots)) telemetry <- dots$telemetry
  if ("hud" %in% names(dots)) telemetry <- dots$hud
  if (is.logical(telemetry)) {
    telemetry <- if (isTRUE(telemetry)) "detailed" else "none"
  } else if (is.character(telemetry) && length(telemetry) > 0) {
    telemetry <- match.arg(tolower(telemetry[1]), c("detailed", "simple", "none"))
  } else {
    telemetry <- "detailed"
  }
  if (!isTRUE(show_fps)) {
    telemetry <- "none"
  }
  is_detailed <- (telemetry == "detailed")

  for (nm in c("smoot_stat", "smooth_stat", "smooth_stats", "smoot_stats", "stats_window", "smooth_window", "telemetry_window")) {
    if (nm %in% names(dots)) {
      smoot_stat <- dots[[nm]]
      break
    }
  }
  if (is.logical(smoot_stat)) {
    smoot_stat <- if (isTRUE(smoot_stat)) 30L else 1L
  }
  smoot_stat <- suppressWarnings(as.integer(smoot_stat[1]))
  if (is.na(smoot_stat) || smoot_stat < 1L) smoot_stat <- 1L

  if (is.null(show_box)) {
    if ("show_boxes" %in% names(dots)) {
      show_box <- isTRUE(dots$show_boxes)
    } else {
      show_box <- isTRUE(bbox)
    }
  } else {
    show_box <- isTRUE(show_box)
  }
  if ("box" %in% names(dots)) show_box <- isTRUE(dots$box)
  if ("boxes" %in% names(dots)) show_box <- isTRUE(dots$boxes)
  if ("bbox" %in% names(dots)) show_box <- isTRUE(dots$bbox)
  if ("show_box" %in% names(dots)) show_box <- isTRUE(dots$show_box)
  bbox <- isTRUE(show_box)

  if ("show_point" %in% names(dots)) show_point <- isTRUE(dots$show_point)
  if ("show_points" %in% names(dots)) show_point <- isTRUE(dots$show_points)
  if ("point" %in% names(dots)) show_point <- isTRUE(dots$point)
  if ("points" %in% names(dots)) show_point <- isTRUE(dots$points)
  show_point <- isTRUE(show_point)
  if ("point_radius" %in% names(dots)) point_radius <- dots$point_radius
  if ("pt_radius" %in% names(dots)) point_radius <- dots$pt_radius
  point_radius <- as.integer(point_radius)
  if ("point_col" %in% names(dots)) point_col <- dots$point_col
  if ("point_color" %in% names(dots)) point_col <- dots$point_color

  if ("min_frames" %in% names(dots)) min_frames <- dots$min_frames
  if ("min_hits" %in% names(dots)) min_frames <- dots$min_hits
  min_hits_val <- if (!is.null(min_frames) && as.integer(min_frames) > 0L) as.integer(min_frames) else 1L
  if ((!missing(min_frames) || "min_frames" %in% names(dots) || "min_hits" %in% names(dots)) && min_hits_val > 1L) {
    track <- TRUE
  }

  if ("hud_pos" %in% names(dots)) hud_pos <- dots$hud_pos
  if ("pos" %in% names(dots)) hud_pos <- dots$pos
  if ("position" %in% names(dots)) hud_pos <- dots$position
  if ("hud_position" %in% names(dots)) hud_pos <- dots$hud_position
  if ("hud_layout" %in% names(dots)) hud_layout <- dots$hud_layout
  if ("layout" %in% names(dots)) hud_layout <- dots$layout
  if ("hud_style" %in% names(dots)) hud_layout <- dots$hud_style
  if ("style" %in% names(dots)) hud_layout <- dots$style

  valid_pos <- c("top-right", "top-left", "bottom-left", "bottom-right", "top", "bottom")
  hud_pos <- if (is.numeric(hud_pos) && length(hud_pos) >= 2) {
    paste0(round(as.numeric(hud_pos[1:2]), 4), collapse = ",")
  } else if (is.character(hud_pos) && length(hud_pos) > 0) {
    p_first <- tolower(hud_pos[1])
    if (grepl("^[0-9.]+[ ,_][0-9.]+$|^[0-9.]+,[0-9.]+$", p_first)) {
      p_first
    } else {
      match.arg(p_first, valid_pos)
    }
  } else {
    "top-right"
  }

  valid_layout <- c("vertical", "horizontal")
  hud_layout <- if (is.character(hud_layout) && length(hud_layout) > 0) {
    l_match <- match.arg(tolower(hud_layout[1]), valid_layout)
    l_match
  } else "vertical"

  if ("infer_every" %in% names(dots)) infer_every <- dots$infer_every
  if ("infer_each" %in% names(dots)) infer_every <- dots$infer_each
  if ("infer_step" %in% names(dots)) infer_every <- dots$infer_step
  if ("step" %in% names(dots)) infer_every <- dots$step
  if ("stride" %in% names(dots)) infer_every <- dots$stride
  if ("skip_frames" %in% names(dots)) infer_every <- as.integer(dots$skip_frames) + 1L
  infer_every <- max(1L, as.integer(infer_every))

  default_obj_size <- 1.0
  default_tel_size <- 1.0
  default_roi_size <- 1.0

  if (is.numeric(text_size) && length(text_size) >= 2) {
    nms <- names(text_size)
    if (!is.null(nms)) {
      obj_idx <- grep("obj|box|label", nms, ignore.case = TRUE)
      tel_idx <- grep("telemetry|hud", nms, ignore.case = TRUE)
      if (length(obj_idx) > 0) default_obj_size <- as.numeric(text_size[obj_idx[1]])
      if (length(tel_idx) > 0) default_tel_size <- as.numeric(text_size[tel_idx[1]])
      if (length(obj_idx) == 0 && length(tel_idx) == 0) {
        default_obj_size <- as.numeric(text_size[1])
        default_tel_size <- as.numeric(text_size[2])
      }
    } else {
      default_obj_size <- as.numeric(text_size[1])
      default_tel_size <- as.numeric(text_size[2])
    }
  } else {
    val <- suppressWarnings(as.numeric(text_size[1]))
    if (!is.na(val)) {
      default_obj_size <- val
      default_tel_size <- val
    }
  }

  obj_sz <- if (!is.null(text_size_objects)) suppressWarnings(as.numeric(text_size_objects[1])) else default_obj_size
  if (is.na(obj_sz)) obj_sz <- 1.0

  tel_sz <- if (!is.null(text_size_telemetry)) suppressWarnings(as.numeric(text_size_telemetry[1])) else default_tel_size
  if (is.na(tel_sz)) tel_sz <- 1.0

  roi_sz <- if (!is.null(text_size_roi)) suppressWarnings(as.numeric(text_size_roi[1])) else default_roi_size
  if (is.na(roi_sz)) roi_sz <- 1.0

  if (isFALSE(show_text) || obj_sz <= 0.05) {
    show_text <- FALSE
    f_scale_obj <- 0.0
  } else {
    f_scale_obj <- obj_sz
  }

  f_scale_hud <- if (tel_sz <= 0.05) 0.0 else tel_sz
  text_size_num <- f_scale_obj
  pad_num <- if (is.numeric(pad) && length(pad) > 0) pad[1] else 1.0
  pad_num <- max(0.0, as.numeric(pad_num))
  badge <- isTRUE(badge)

  # For small object text (<= 0.3), hide redundant class name and confidence score by default
  # unless the user explicitly requested them in arguments / dots
  if (!"show_conf" %in% names(dots) && !"show_scores" %in% names(dots) && !"scores" %in% names(dots) && f_scale_obj <= 0.3) {
    show_conf <- FALSE
  }
  if (!"show_class" %in% names(dots) && !"show_labels" %in% names(dots) && f_scale_obj <= 0.3 && (isTRUE(track) || !is.null(count_line))) {
    show_class <- FALSE
  }
  track_max_dist <- if ("max_dist" %in% names(dots)) as.numeric(dots$max_dist) else 150.0
  track_min_iou <- if ("min_iou" %in% names(dots)) as.numeric(dots$min_iou) else 0.25
  track_max_lost <- if ("max_lost" %in% names(dots)) as.integer(dots$max_lost) else 15L
  track_max_history <- if ("max_history" %in% names(dots)) as.integer(dots$max_history) else 30L
  if ("interpolate_box" %in% names(dots)) interpolate_box <- dots$interpolate_box
  if ("interp_box" %in% names(dots)) interpolate_box <- dots$interp_box
  if ("interpolate_mask" %in% names(dots)) interpolate_mask <- dots$interpolate_mask
  if ("interp_mask" %in% names(dots)) interpolate_mask <- dots$interp_mask
  interpolate_box <- isTRUE(interpolate_box)
  interpolate_mask <- isTRUE(interpolate_mask)
  if ("count_line_dir" %in% names(dots)) {
    count_line_orientation <- dots$count_line_dir
  }
  count_line_orientation <- match.arg(count_line_orientation, c("vertical", "horizontal"))
  count_direction <- match.arg(tolower(count_direction[1]), c("both", "forward", "backward"))
  if ("roi_dir" %in% names(dots)) {
    roi_orientation <- dots$roi_dir
  }
  roi_orientation <- match.arg(roi_orientation, c("vertical", "horizontal"))
  roi_mode <- match.arg(roi_mode)
  backend <- match.arg(backend)
  if (backend == "auto") backend <- "native"
  if ("sat" %in% names(dots) && missing(saturation)) saturation <- dots$sat
  if ("cont" %in% names(dots) && missing(contrast)) contrast <- dots$cont
  saturation_val <- max(0.0, as.numeric(saturation[1]))
  contrast_val <- max(0.0, as.numeric(contrast[1]))
  if (is.na(saturation_val)) saturation_val <- 1.0
  if (is.na(contrast_val)) contrast_val <- 1.0
  need_color_adj <- (abs(saturation_val - 1.0) > 0.005 || abs(contrast_val - 1.0) > 0.005)
  heatmap <- isTRUE(heatmap)
  heatmap_alpha <- max(0.05, min(1.0, as.numeric(heatmap_alpha[1])))
  heatmap_decay <- max(0.0, min(1.0, as.numeric(heatmap_decay[1])))
  sahi <- isTRUE(sahi)
  save_crops <- isTRUE(save_crops) || !is.null(crops_dir)
  if (save_crops) {
    if (is.null(crops_dir) || !nzchar(as.character(crops_dir[1]))) {
      crops_dir <- file.path(getwd(), "crops")
    }
    crops_dir <- normalizePath(crops_dir[1], winslash = "/", mustWork = FALSE)
    if (!dir.exists(crops_dir)) dir.create(crops_dir, recursive = TRUE, showWarnings = FALSE)
  }
  best_crop_scores <- new.env(parent = emptyenv())
  start_sec <- if (!is.null(start_time) && is.numeric(start_time) && start_time[1] > 0) as.numeric(start_time[1]) else 0.0
  end_sec <- if (!is.null(end_time) && is.numeric(end_time) && end_time[1] > 0) as.numeric(end_time[1]) else 0.0

  .populate_seg_stats_df <- function(df_b, raw_res, num_boxes_orig) {
    if (nrow(df_b) == 0L) return(df_b)
    if (!is.null(raw_res$seg_stats) && is.matrix(raw_res$seg_stats) && ncol(raw_res$seg_stats) >= 8L) {
      st <- raw_res$seg_stats
      m_idx <- df_b$orig_idx
      valid_m <- which(!is.na(m_idx) & m_idx > 0L & m_idx <= nrow(st))
      df_b$area <- NA_real_
      df_b$perimeter <- NA_real_
      df_b$radius_mean <- NA_real_
      df_b$length <- NA_real_
      df_b$width <- NA_real_
      df_b$circularity <- NA_real_
      df_b$eccentricity <- NA_real_
      df_b$asp_ratio <- NA_real_
      if (length(valid_m) > 0L) {
        idx_v <- m_idx[valid_m]
        df_b$area[valid_m] <- round(as.numeric(st[idx_v, 1]), 1)
        df_b$perimeter[valid_m] <- round(as.numeric(st[idx_v, 2]), 1)
        df_b$radius_mean[valid_m] <- round(as.numeric(st[idx_v, 3]), 2)
        df_b$length[valid_m] <- round(as.numeric(st[idx_v, 4]), 1)
        df_b$width[valid_m] <- round(as.numeric(st[idx_v, 5]), 1)
        df_b$circularity[valid_m] <- round(as.numeric(st[idx_v, 6]), 4)
        df_b$eccentricity[valid_m] <- round(as.numeric(st[idx_v, 7]), 4)
        df_b$asp_ratio[valid_m] <- round(as.numeric(st[idx_v, 8]), 3)
      }
    }
    df_b
  }

  .save_object_crops_frame <- function(bm, df_b, frame_idx) {
    if (!save_crops || nrow(df_b) == 0L) return(invisible(NULL))
    b_dims <- dim(bm)
    fw <- b_dims[2]
    fh <- b_dims[3]
    for (r in seq_len(nrow(df_b))) {
      oid <- df_b$id[r]
      sc <- df_b$score[r]
      lbl <- gsub("[^a-zA-Z0-9_-]", "_", as.character(df_b$label[r]))
      if (isTRUE(track) && !is.na(oid) && oid > 0L) {
        key <- as.character(oid)
        prev_sc <- if (exists(key, envir = best_crop_scores, inherits = FALSE)) get(key, envir = best_crop_scores) else -1.0
        if (sc <= prev_sc) next
        assign(key, sc, envir = best_crop_scores)
        fn <- file.path(crops_dir, sprintf("id_%04d_%s.jpg", oid, lbl))
      } else {
        fn <- file.path(crops_dir, sprintf("frame_%06d_obj_%03d_%s.jpg", frame_idx, r, lbl))
      }
      cx1 <- max(0L, min(fw - 2L, as.integer(floor(df_b$xmin[r]))))
      cy1 <- max(0L, min(fh - 2L, as.integer(floor(df_b$ymin[r]))))
      cx2 <- max(cx1 + 1L, min(fw, as.integer(ceiling(df_b$xmax[r]))))
      cy2 <- max(cy1 + 1L, min(fh, as.integer(ceiling(df_b$ymax[r]))))
      crop_bm <- crop_bgr_cpp(bm, cx1, cy1, cx2, cy2)
      save_frame_bgr_cpp(crop_bm, fn, quality = 92L)
    }
  }

  # Safely wipe any dangling ghost progress bars from previously interrupted sessions
  try({
    cli_ns <- asNamespace("cli")
    if (exists("clienv", envir = cli_ns, inherits = FALSE)) {
      assign("progress", list(), envir = cli_ns$clienv)
    }
  }, silent = TRUE)

  on.exit({
    try(cli::cli_progress_cleanup(), silent = TRUE)
    try({
      cli_ns <- asNamespace("cli")
      if (exists("clienv", envir = cli_ns, inherits = FALSE)) {
        assign("progress", list(), envir = cli_ns$clienv)
      }
    }, silent = TRUE)
  }, add = TRUE)

  if (save_video && is.null(output)) {
    output <- if (detect) "annotated_video.mp4" else "recorded_video.mp4"
  }
  if (!save_video) {
    output <- NULL
  }

  out_video_file <- if (save_video && !is.null(output)) {
    normalizePath(output, winslash = "/", mustWork = FALSE)
  } else {
    ""
  }

  engine <- match.arg(engine)
  use_gpu <- (engine == "gpu")
  lib_path <- pliman_onnx_library_path()

  if (detect) {
    # Model resolution: custom ONNX file or pre-trained download
    local_hit <- .find_local_model(model[1], dir = dir)
    if (!is.null(local_hit)) {
      model_file <- local_hit
    } else {
      model_file <- pliman_download_model(model = model[1], dir = dir)
    }

    # Label resolution: custom vector -> ONNX embedded metadata -> COCO fallback
    if (!is.null(labels)) {
      class_names <- labels
    } else {
      model_classes <- .get_onnx_classes(model_file)
      if (!is.null(model_classes) && length(model_classes) > 0L) {
        class_names <- model_classes
      } else {
        class_names <- .coco_classes
      }
    }
    # Inference resolution (imgsz)
    model_shape_info <- tryCatch(
      get_onnx_model_io_shape_cpp(model_file, lib_path, as.integer(threads), use_gpu, as.integer(device_id)),
      error = function(e) list(dynamic = TRUE, expected_w = 640L, expected_h = 640L)
    )
    target_size <- if (!is.null(imgsz)) {
      as.integer(imgsz[1])
    } else {
      model_imgsz <- .get_onnx_imgsz(model_file)
      if (!is.null(model_imgsz) && length(model_imgsz) > 0L) as.integer(model_imgsz[1]) else model_shape_info$expected_w
    }
    if (is.na(target_size) || target_size <= 0L) target_size <- 640L

    if (!isTRUE(model_shape_info$dynamic) && target_size != model_shape_info$expected_w) {
      cand_pt <- file.path(dirname(model_file), "best.pt")
      if (!file.exists(cand_pt)) cand_pt <- sub("\\.onnx$", ".pt", model_file)
      reexported <- FALSE
      if (file.exists(cand_pt)) {
        cli::cli_alert_info("Model {.file {basename(model_file)}} has fixed static input shape ({model_shape_info$expected_w}x{model_shape_info$expected_h}). Auto-upgrading to dynamic ONNX via {.file {basename(cand_pt)}}...")
        tryCatch({
          new_ox <- yolo_export(cand_pt, output_dir = dirname(model_file), output_name = basename(model_file), dynamic = TRUE, auto_install = FALSE)
          if (file.exists(new_ox)) {
            model_file <- new_ox
            reexported <- TRUE
            model_shape_info$dynamic <- TRUE
            cli::cli_alert_success("Model upgraded to dynamic ONNX: now supports any imgsz (e.g., {target_size})!")
          }
        }, error = function(e) NULL)
      }
      if (!reexported) {
        cli::cli_alert_warning(c(
          "!" = "Model {.file {basename(model_file)}} has fixed static input resolution of {model_shape_info$expected_w}x{model_shape_info$expected_h}.",
          "i" = "Requested {.code imgsz = {target_size}} cannot be applied to a static model. Automatically using native resolution ({model_shape_info$expected_w}).",
          ">" = "To enable dynamic multi-resolution inference, export with: {.code yolo_export(model, dynamic = TRUE)}"
        ))
        target_size <- model_shape_info$expected_w
      }
    }
  } else {
    model_file <- "None"
    class_names <- character(0)
    target_size <- 640L
  }

  .filter_nested_mat <- function(boxes, scores, class_ids, stats = NULL,
                                 containment_thresh = 0.85, max_area_ratio = 0.70, conf_floor = 0.50) {
    n <- nrow(boxes)
    if (is.null(n) || n <= 1L) return(list(boxes = boxes, scores = scores, class_ids = class_ids, seg_stats = stats))
    x1 <- boxes[, 1]; y1 <- boxes[, 2]; x2 <- boxes[, 3]; y2 <- boxes[, 4]
    areas <- pmax(0.0, x2 - x1) * pmax(0.0, y2 - y1)
    valid_child <- which(scores >= conf_floor & areas > 0.0)
    if (length(valid_child) == 0L) return(list(boxes = boxes, scores = scores, class_ids = class_ids, seg_stats = stats))
    keep <- rep(TRUE, n)
    for (i in seq_len(n)) {
      ai <- areas[i]
      if (ai <= 0.0) next
      bx1 <- x1[i]; by1 <- y1[i]; bx2 <- x2[i]; by2 <- y2[i]
      ci <- valid_child[valid_child != i & areas[valid_child] < ai * max_area_ratio]
      if (length(ci) == 0L) next
      ci <- ci[x2[ci] > bx1 & x1[ci] < bx2 & y2[ci] > by1 & y1[ci] < by2]
      if (length(ci) == 0L) next
      iw <- pmin(bx2, x2[ci]) - pmax(bx1, x1[ci])
      ih <- pmin(by2, y2[ci]) - pmax(by1, y1[ci])
      inter <- pmax(0.0, iw) * pmax(0.0, ih)
      if (any(inter / areas[ci] >= containment_thresh)) keep[i] <- FALSE
    }
    if (all(keep)) return(list(boxes = boxes, scores = scores, class_ids = class_ids, seg_stats = stats))
    list(
      boxes = boxes[keep, , drop = FALSE],
      scores = scores[keep],
      class_ids = class_ids[keep],
      seg_stats = if (!is.null(stats) && nrow(stats) == n) stats[keep, , drop = FALSE] else stats
    )
  }

  .sahi_merge_video <- function(boxes, scores, class_ids, is_full, stats = NULL,
                                match_threshold = 0.50) {
    n <- nrow(boxes)
    if (is.null(n) || n <= 1L) return(list(boxes = boxes, scores = scores, class_ids = class_ids, seg_stats = stats))
    has_stats <- !is.null(stats) && nrow(stats) == n
    cls_groups <- split(seq_len(n), class_ids)

    out_b_list <- list()
    out_s_list <- list()
    out_c_list <- list()
    out_st_list <- list()

    for (g_idx in cls_groups) {
      x1 <- boxes[g_idx, 1]; y1 <- boxes[g_idx, 2]
      x2 <- boxes[g_idx, 3]; y2 <- boxes[g_idx, 4]
      cf <- scores[g_idx]
      cid <- class_ids[g_idx[1L]]
      isf <- is_full[g_idx]
      st_g <- if (has_stats) stats[g_idx, , drop = FALSE] else NULL

      # Step 1: Subsume slice fragments into full-image detections
      if (any(isf)) {
        full_idx <- which(isf)
        slice_idx <- which(!isf)
        if (length(slice_idx) > 0L) {
          sx1 <- x1[slice_idx]; sy1 <- y1[slice_idx]
          sx2 <- x2[slice_idx]; sy2 <- y2[slice_idx]
          scf <- cf[slice_idx]
          s_area <- (sx2 - sx1 + 1.0) * (sy2 - sy1 + 1.0)
          absorbed <- logical(length(slice_idx))
          for (fi in full_idx) {
            rem_s <- which(!absorbed)
            if (length(rem_s) == 0L) break
            fx1 <- x1[fi]; fy1 <- y1[fi]; fx2 <- x2[fi]; fy2 <- y2[fi]
            rem_s <- rem_s[sx2[rem_s] >= fx1 & sx1[rem_s] <= fx2 & sy2[rem_s] >= fy1 & sy1[rem_s] <= fy2]
            if (length(rem_s) == 0L) next
            f_area <- (fx2 - fx1 + 1.0) * (fy2 - fy1 + 1.0)
            iw <- pmin(fx2, sx2[rem_s]) - pmax(fx1, sx1[rem_s]) + 1.0
            ih <- pmin(fy2, sy2[rem_s]) - pmax(fy1, sy1[rem_s]) + 1.0
            inter <- iw * ih
            sa <- s_area[rem_s]
            ios <- inter / pmin(f_area, sa)
            iou <- inter / (f_area + sa - inter)
            hit <- rem_s[ios >= 0.40 | iou >= 0.35]
            if (length(hit) > 0L) {
              max_s <- max(scf[hit])
              if (max_s > cf[fi]) cf[fi] <- max_s
              absorbed[hit] <- TRUE
            }
          }
          if (any(absorbed)) {
            keep_r <- rep(TRUE, length(x1))
            keep_r[slice_idx[absorbed]] <- FALSE
            x1 <- x1[keep_r]; y1 <- y1[keep_r]; x2 <- x2[keep_r]; y2 <- y2[keep_r]
            cf <- cf[keep_r]
            if (has_stats) st_g <- st_g[keep_r, , drop = FALSE]
          }
        }
      }

      n_rem <- length(x1)
      if (n_rem == 0L) next

      # Step 2: Vectorized Greedy NMM with IOS
      ord <- order(cf, decreasing = TRUE)
      x1 <- x1[ord]; y1 <- y1[ord]; x2 <- x2[ord]; y2 <- y2[ord]; cf <- cf[ord]
      if (has_stats) st_g <- st_g[ord, , drop = FALSE]
      areas <- (x2 - x1 + 1.0) * (y2 - y1 + 1.0)

      alive <- rep(TRUE, n_rem)
      keep_seeds <- integer(n_rem)
      ox1 <- numeric(n_rem); oy1 <- numeric(n_rem)
      ox2 <- numeric(n_rem); oy2 <- numeric(n_rem)
      ocf <- numeric(n_rem)
      n_out <- 0L

      for (seed in seq_len(n_rem)) {
        if (!alive[seed]) next
        alive[seed] <- FALSE
        n_out <- n_out + 1L
        keep_seeds[n_out] <- seed
        bx1 <- x1[seed]; by1 <- y1[seed]; bx2 <- x2[seed]; by2 <- y2[seed]
        b_area <- areas[seed]

        cand <- which(alive)
        if (length(cand) > 0L) {
          cand <- cand[x2[cand] >= bx1 & x1[cand] <= bx2]
          if (length(cand) > 0L) cand <- cand[y2[cand] >= by1 & y1[cand] <= by2]
        }
        if (length(cand) > 0L) {
          iw <- pmin(bx2, x2[cand]) - pmax(bx1, x1[cand]) + 1.0
          ih <- pmin(by2, y2[cand]) - pmax(by1, y1[cand]) + 1.0
          inter <- iw * ih
          ca <- areas[cand]
          iou <- inter / (b_area + ca - inter)
          ios <- inter / pmin(b_area, ca)
          matched <- cand[ios >= match_threshold & iou >= 0.30]
        } else {
          matched <- integer(0)
        }
        if (length(matched) > 0L) {
          alive[matched] <- FALSE
          m_idx <- c(seed, matched)
          ox1[n_out] <- min(x1[m_idx]); oy1[n_out] <- min(y1[m_idx])
          ox2[n_out] <- max(x2[m_idx]); oy2[n_out] <- max(y2[m_idx])
          ocf[n_out] <- max(cf[m_idx])
        } else {
          ox1[n_out] <- bx1; oy1[n_out] <- by1
          ox2[n_out] <- bx2; oy2[n_out] <- by2
          ocf[n_out] <- cf[seed]
        }
      }
      if (n_out > 0L) {
        s_idx <- seq_len(n_out)
        out_b_list[[length(out_b_list) + 1L]] <- cbind(ox1[s_idx], oy1[s_idx], ox2[s_idx], oy2[s_idx])
        out_s_list[[length(out_s_list) + 1L]] <- ocf[s_idx]
        out_c_list[[length(out_c_list) + 1L]] <- rep(cid, n_out)
        if (has_stats) out_st_list[[length(out_st_list) + 1L]] <- st_g[keep_seeds[s_idx], , drop = FALSE]
      }
    }

    list(
      boxes = if (length(out_b_list) == 1L) out_b_list[[1L]] else do.call(rbind, out_b_list),
      scores = unlist(out_s_list),
      class_ids = unlist(out_c_list),
      seg_stats = if (has_stats && length(out_st_list) > 0L) do.call(rbind, out_st_list) else stats
    )
  }

  .run_video_yolo <- function(bm_feed, feed_w, feed_h, need_masks_flag) {
    if (!sahi) {
      res_single <- run_yolo_cpp(
        tensor_vec = bm_feed,
        orig_w = feed_w,
        orig_h = feed_h,
        conf_threshold = conf_threshold,
        iou_threshold = iou_threshold,
        model_path = model_file,
        lib_path = lib_path,
        num_threads = as.integer(threads),
        use_gpu = use_gpu,
        device_id = as.integer(device_id),
        need_masks = need_masks_flag,
        target_size = as.integer(target_size)
      )
      if (!need_masks_flag && nrow(res_single$boxes) > 1L) {
        fn <- .filter_nested_mat(res_single$boxes, res_single$scores, res_single$class_ids, res_single$seg_stats)
        res_single$boxes <- fn$boxes
        res_single$scores <- fn$scores
        res_single$class_ids <- fn$class_ids
        res_single$seg_stats <- fn$seg_stats
      }
      return(res_single)
    }
    sh <- if (!is.null(slice_height) && slice_height > 0) {
      as.integer(slice_height[1])
    } else {
      min(640L, max(256L, as.integer(round(feed_h * 0.60))))
    }
    sw <- if (!is.null(slice_width) && slice_width > 0) {
      as.integer(slice_width[1])
    } else {
      min(640L, max(256L, as.integer(round(feed_w * 0.58))))
    }
    sh <- min(sh, feed_h)
    sw <- min(sw, feed_w)
    tile_target_size <- if (isTRUE(model_shape_info$dynamic)) {
      min(as.integer(target_size), max(640L, ((max(sw, sh) + 31L) %/% 32L) * 32L))
    } else {
      as.integer(target_size)
    }
    ov <- max(0.0, min(0.8, as.numeric(overlap_ratio[1])))
    step_x <- max(1L, as.integer(round(sw * (1.0 - ov))))
    step_y <- max(1L, as.integer(round(sh * (1.0 - ov))))
    xs <- unique(c(seq(0L, max(0L, feed_w - sw), by = step_x), max(0L, feed_w - sw)))
    ys <- unique(c(seq(0L, max(0L, feed_h - sh), by = step_y), max(0L, feed_h - sh)))

    full_res <- run_yolo_cpp(
      tensor_vec = bm_feed,
      orig_w = feed_w,
      orig_h = feed_h,
      conf_threshold = conf_threshold,
      iou_threshold = iou_threshold,
      model_path = model_file,
      lib_path = lib_path,
      num_threads = as.integer(threads),
      use_gpu = use_gpu,
      device_id = as.integer(device_id),
      need_masks = need_masks_flag,
      target_size = as.integer(target_size)
    )
    if (nrow(full_res$boxes) > 1L) {
      fn <- .filter_nested_mat(full_res$boxes, full_res$scores, full_res$class_ids, full_res$seg_stats)
      full_res$boxes <- fn$boxes
      full_res$scores <- fn$scores
      full_res$class_ids <- fn$class_ids
      full_res$seg_stats <- fn$seg_stats
    }
    boxes_list <- list(full_res$boxes)
    scores_list <- list(full_res$scores)
    cls_list <- list(full_res$class_ids)
    isf_list <- list(rep(TRUE, nrow(full_res$boxes)))
    stats_list <- if (!is.null(full_res$seg_stats)) list(full_res$seg_stats) else NULL

    for (ty in ys) {
      for (tx in xs) {
        tile_bm <- crop_bgr_cpp(bm_feed, tx, ty, tx + sw, ty + sh)
        tres <- run_yolo_cpp(
          tensor_vec = tile_bm,
          orig_w = sw,
          orig_h = sh,
          conf_threshold = conf_threshold,
          iou_threshold = iou_threshold,
          model_path = model_file,
          lib_path = lib_path,
          num_threads = as.integer(threads),
          use_gpu = use_gpu,
          device_id = as.integer(device_id),
          need_masks = FALSE,
          target_size = tile_target_size
        )
        if (nrow(tres$boxes) > 0L) {
          if (nrow(tres$boxes) > 1L) {
            fn_t <- .filter_nested_mat(tres$boxes, tres$scores, tres$class_ids, tres$seg_stats)
            tres$boxes <- fn_t$boxes
            tres$scores <- fn_t$scores
            tres$class_ids <- fn_t$class_ids
            tres$seg_stats <- fn_t$seg_stats
          }
          tb <- tres$boxes
          tb[, 1] <- tb[, 1] + tx
          tb[, 2] <- tb[, 2] + ty
          tb[, 3] <- tb[, 3] + tx
          tb[, 4] <- tb[, 4] + ty
          boxes_list[[length(boxes_list) + 1L]] <- tb
          scores_list[[length(scores_list) + 1L]] <- tres$scores
          cls_list[[length(cls_list) + 1L]] <- tres$class_ids
          isf_list[[length(isf_list) + 1L]] <- rep(FALSE, nrow(tb))
          if (!is.null(stats_list) && !is.null(tres$seg_stats)) {
            stats_list[[length(stats_list) + 1L]] <- tres$seg_stats
          }
        }
      }
    }
    all_b <- do.call(rbind, boxes_list)
    all_s <- unlist(scores_list)
    all_c <- unlist(cls_list)
    all_isf <- unlist(isf_list)
    all_st <- if (!is.null(stats_list)) do.call(rbind, stats_list) else NULL
    if (nrow(all_b) > 0L) {
      m_res <- .sahi_merge_video(all_b, all_s, all_c, all_isf, all_st, match_threshold = 0.50)
      full_res$boxes <- m_res$boxes
      full_res$scores <- m_res$scores
      full_res$class_ids <- m_res$class_ids
      if (!is.null(stats_list)) {
        full_res$seg_stats <- m_res$seg_stats
      }
    }
    full_res
  }

  has_video_ext <- is.character(video) && tolower(tools::file_ext(video[1])) %in% c("mp4", "avi", "mov", "mkv", "wmv", "flv", "webm", "m4v", "raw", "braw", "mts", "m2ts", "ts", "mpg", "mpeg", "3gp")
  is_camera <- is.numeric(video) || (is.character(video) && !has_video_ext && (
    tolower(video[1]) %in% c("0", "1", "2", "3", "4", "5", "camera", "live", "webcam") ||
    !file.exists(video[1])
  ))
  if (is.null(show)) {
    show <- if (is_camera) interactive() else FALSE
  }
  if (!is_camera && isTRUE(show) && isTRUE(verbose)) {
    cli::cli_alert_info("Interactive frame display enabled (rendering plots to R graphics device reduces throughput).")
  }

  # Initialize object tracking state
  track_state <- list(
    next_id = 1L,
    total_count = 0L,
    tracks = list(),
    counts_by_class = list(),
    crossing_events = list()
  )
  total_count_val <- 0L
  in_zone_val <- 0L

  # ============================================================================
  # MODE A: LIVE WEBCAM STREAMING (NATIVE C++ MEDIA FOUNDATION)
  # ============================================================================
  if (is_camera) {
    if (.Platform$OS.type != "windows" || !has_native_camera_cpp()) {
      cli::cli_abort(c(
        "x" = "Live camera streaming is currently supported natively on Windows (via Windows Media Foundation).",
        "i" = "To analyze video on Linux/macOS, please pass a video file path (e.g., {.path 'sample.mp4'})."
      ))
    }

    if (isTRUE(verbose)) {
      if (detect) {
        cli::cli_alert_info("Starting real-time camera inference with model {.val {basename(model_file)}} [{toupper(engine)} | native C++ engine]")
      } else {
        cli::cli_alert_info("Starting real-time camera stream [pass-through mode, inference disabled | native C++ engine]")
      }
      if (isTRUE(fullscreen)) {
        cli::cli_alert_info("Display mode: Fullscreen enabled.")
      }
      if (!is.null(roi)) {
        cli::cli_alert_info("Region of Interest (ROI) active: mode = {.val {roi_mode}}.")
      }
      if (isTRUE(track) || !is.null(count_line)) {
        cli::cli_alert_info("Real-time object tracking and counting tripwire enabled.")
      }
      cli::cli_alert_info("Press {.kbd ESC} or close the video window to stop.")
    }

    if (save_video) {
      out_video_file <- normalizePath(output, winslash = "/", mustWork = FALSE)
      cam_dir <- file.path(tempdir(), paste0("cam_yolo_", as.integer(stats::runif(1, 1000, 9999))))
      if (!dir.exists(cam_dir)) dir.create(cam_dir, recursive = TRUE, showWarnings = FALSE)
      on.exit(unlink(cam_dir, recursive = TRUE, force = TRUE), add = TRUE)
      recorded_frames <- list()
    }

    frame_records <- list()
    frame_keypoints <- list()
    frame_count <- 0L
    t_start <- Sys.time()
    last_t <- Sys.time()
    fps_smooth <- 0.0

    cam_id <- 0L
    cam_name <- "Default Camera"
      avail_cams <- tryCatch(list_native_cameras_cpp(), error = function(e) list(id = integer(0), name = character(0)))
      if (length(avail_cams$id) > 0) {
        if (is.numeric(video)) {
          idx_req <- as.integer(video[1])
          match_idx <- which(avail_cams$id == idx_req)
          if (length(match_idx) > 0) {
            cam_id <- avail_cams$id[match_idx[1]]
            cam_name <- avail_cams$name[match_idx[1]]
          } else {
            cli::cli_alert_warning("Camera index {idx_req} not found. Available: {paste0('[', avail_cams$id, '] ', avail_cams$name, collapse = ', ')}. Defaulting to [{avail_cams$id[1]}].")
            cam_id <- avail_cams$id[1]
            cam_name <- avail_cams$name[1]
          }
        } else if (is.character(video)) {
          v_num <- suppressWarnings(as.integer(video[1]))
          if (!is.na(v_num) && v_num %in% avail_cams$id) {
            cam_id <- v_num
            cam_name <- avail_cams$name[which(avail_cams$id == v_num)[1]]
          } else {
            matched <- grep(video[1], avail_cams$name, ignore.case = TRUE)
            if (length(matched) > 0) {
              cam_id <- avail_cams$id[matched[1]]
              cam_name <- avail_cams$name[matched[1]]
            } else {
              cli::cli_alert_warning("No camera matching '{video[1]}' found. Available: {paste0('[', avail_cams$id, '] ', avail_cams$name, collapse = ', ')}. Defaulting to [{avail_cams$id[1]}].")
              cam_id <- avail_cams$id[1]
              cam_name <- avail_cams$name[1]
            }
          }
        }
      }

      target_w <- 0L
      target_h <- 0L
      if (!is.null(resolution) && length(resolution) >= 2) {
        target_w <- as.integer(resolution[1])
        target_h <- as.integer(resolution[2])
      }

      cam <- tryCatch({
        open_native_camera_cpp(cam_id, target_w, target_h)
      }, error = function(e) {
        cli::cli_abort(c(
          "x" = "Failed to open native camera {cam_id} ({cam_name}): {e$message}",
          "i" = "Please check if your webcam is connected and not in use by another application."
        ))
      })
      on.exit(close_native_camera_cpp(cam), add = TRUE)

      cam_dims <- get_camera_dims_cpp(cam)
      cam_w <- cam_dims$width
      cam_h <- cam_dims$height

      # Detect if square aspect ratio or custom crop requested
      is_square <- isTRUE(crop) ||
                   (is.character(crop) && tolower(crop[1]) %in% c("1:1", "square")) ||
                   (!is.null(resolution) && length(resolution) == 2 && resolution[1] == resolution[2]) ||
                   (!is.null(resolution) && length(resolution) == 1 && is.numeric(resolution))

      crop_x <- 0L
      crop_y <- 0L
      crop_w <- 0L
      crop_h <- 0L
      orig_w <- cam_w
      orig_h <- cam_h

      if (is_square) {
        req_size <- if (!is.null(resolution)) as.integer(resolution[1]) else min(cam_w, cam_h)
        sq_size <- min(req_size, cam_w, cam_h)
        crop_w <- as.integer(sq_size)
        crop_h <- as.integer(sq_size)
        crop_x <- as.integer(round((cam_w - sq_size) / 2))
        crop_y <- as.integer(round((cam_h - sq_size) / 2))
        orig_w <- crop_w
        orig_h <- crop_h
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Square frame active: 1:1 aspect ratio ({orig_w}x{orig_h}) centered from {cam_w}x{cam_h} camera feed.")
        }
      } else if (is.numeric(crop) && length(crop) >= 2) {
        req_w <- min(as.integer(crop[1]), cam_w)
        req_h <- min(as.integer(crop[2]), cam_h)
        crop_w <- req_w
        crop_h <- req_h
        crop_x <- as.integer(round((cam_w - req_w) / 2))
        crop_y <- as.integer(round((cam_h - req_h) / 2))
        orig_w <- crop_w
        orig_h <- crop_h
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Custom cropped frame: {orig_w}x{orig_h} centered from {cam_w}x{cam_h} camera feed.")
        }
      } else {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Active camera: [{cam_id}] {.strong {cam_name}} ({orig_w}x{orig_h} | {round(orig_w/orig_h, 2)}:1 FOV)")
        }
      }

      bm <- raw(3L * orig_w * orig_h)
      attr(bm, "dim") <- c(3L, orig_w, orig_h)
      heatmap_accum <- if (heatmap) matrix(0.0, nrow = orig_w, ncol = orig_h) else NULL
      need_masks_cam <- isTRUE(mask) && (isTRUE(save_video) || isTRUE(show))

      cam_writer <- NULL
      if (save_video && has_native_video_writer_cpp()) {
        rec_fps_init <- if (!is.null(fps) && fps > 0) as.numeric(fps) else 24.0
        cam_writer <- tryCatch(
          open_native_video_writer_cpp(
            video_path = out_video_file,
            width = orig_w,
            height = orig_h,
            fps = rec_fps_init
          ),
          error = function(e) NULL
        )
        if (!is.null(cam_writer)) {
          on.exit(try(close_native_video_writer_cpp(cam_writer), silent = TRUE), add = TRUE)
        }
      }

      win <- NULL
      if (isTRUE(show)) {
        w_w <- if (!is.null(window_size)) window_size[1] else round(orig_w * window_scale)
        w_h <- if (!is.null(window_size)) window_size[2] else round(orig_h * window_scale)
        win_title <- if (detect) {
          sprintf("pliman - YOLO [%s | %s] (ESC to exit, Space to pause, F for fullscreen)", toupper(engine), basename(model_file))
        } else {
          sprintf("pliman - Camera Stream [%s] (ESC to exit, Space to pause, F for fullscreen)", cam_name)
        }
        win <- create_native_window_cpp(title = win_title, width = as.integer(w_w), height = as.integer(w_h), fullscreen = isTRUE(fullscreen))
        on.exit(close_native_window_cpp(win), add = TRUE)
      }

      last_df_b <- NULL
      last_raw_res <- NULL
      last_has_kpts <- FALSE
      last_has_masks <- FALSE
      last_kpt_mat_pass <- NULL
      last_track_ids_pass <- NULL
      last_flash_pass <- NULL
      last_history_x_pass <- NULL
      last_history_y_pass <- NULL
      last_total_count_val <- 0L
      last_in_zone_val <- 0L
      last_use_crop <- FALSE
      last_track_res <- NULL
      t_inf_ms <- 0.0
      t_trk_ms <- 0.0
      last_gui_ms <- 14.0
      stats_smoother <- .make_stats_smoother(smoot_stat)

      while (TRUE) {
        if (!is.null(max_frames) && frame_count >= max_frames) {
          break
        }

        has_frame <- grab_native_frame_cpp(cam, bm, crop_x1 = crop_x, crop_y1 = crop_y, crop_w = crop_w, crop_h = crop_h)
        if (!has_frame) {
          break
        }

        t_cur <- Sys.time()
        dt <- as.numeric(t_cur - last_t, units = "secs")
        last_t <- t_cur
        cur_fps <- if (dt > 0) (1.0 / dt) else 30.0
        fps_smooth <- if (fps_smooth == 0) cur_fps else (0.85 * fps_smooth + 0.15 * cur_fps)

        frame_count <- frame_count + 1L
        t_elapsed <- as.numeric(t_cur - t_start, units = "secs")

        # ROI coordinates calculation
        roi_coords <- .resolve_roi(roi, orig_w, orig_h, roi_orientation)

        # Virtual count line coordinates calculation
        line_coords <- .resolve_count_line(count_line, orig_w, orig_h, count_line_orientation)

        df_b <- data.frame(
          frame = integer(0), timestamp = numeric(0), id = integer(0),
          xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
          label = character(0), score = numeric(0), class_id = integer(0),
          stringsAsFactors = FALSE
        )
        has_kpts <- FALSE
        has_masks <- FALSE
        kpt_mat_pass <- NULL
        track_ids_pass <- NULL
        flash_pass <- NULL
        history_x_pass <- NULL
        history_y_pass <- NULL

        use_crop <- FALSE
        should_infer <- (frame_count == 1L || (frame_count %% infer_every == 0L))

        if (detect && should_infer) {
          use_crop <- (roi_mode == "crop" && !is.null(roi_coords))
          if (use_crop) {
            margin_x <- round(0.08 * orig_w)
            margin_y <- round(0.08 * orig_h)
            crop_x1 <- max(0L, as.integer(roi_coords[1] - margin_x))
            crop_y1 <- max(0L, as.integer(roi_coords[2] - margin_y))
            crop_x2 <- min(orig_w, as.integer(roi_coords[3] + margin_x))
            crop_y2 <- min(orig_h, as.integer(roi_coords[4] + margin_y))
            crop_coords <- c(crop_x1, crop_y1, crop_x2, crop_y2)
            bm_feed <- crop_bgr_cpp(bm, crop_x1, crop_y1, crop_x2, crop_y2)
          } else {
            crop_coords <- c(0L, 0L, orig_w, orig_h)
            bm_feed <- bm
          }
          feed_dims <- dim(bm_feed)
          feed_w <- feed_dims[2]
          feed_h <- feed_dims[3]

          t_infer_start <- Sys.time()
          raw_res <- .run_video_yolo(bm_feed, feed_w, feed_h, need_masks_cam)
          t_infer_end <- Sys.time()
          t_inf_ms <- as.numeric(t_infer_end - t_infer_start, units = "secs") * 1000

          has_masks <- need_masks_cam && isTRUE(raw_res$is_seg) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && nrow(raw_res$boxes) > 0L
          num_boxes_orig <- nrow(raw_res$boxes)
          orig_idx <- seq_len(num_boxes_orig)

          # Shift box coordinates back to full frame if crop was used
          if (use_crop && nrow(raw_res$boxes) > 0) {
            raw_res$boxes[, 1] <- raw_res$boxes[, 1] + crop_coords[1]
            raw_res$boxes[, 2] <- raw_res$boxes[, 2] + crop_coords[2]
            raw_res$boxes[, 3] <- raw_res$boxes[, 3] + crop_coords[1]
            raw_res$boxes[, 4] <- raw_res$boxes[, 4] + crop_coords[2]
            if (!is.null(raw_res$keypoints) && ncol(raw_res$keypoints) >= 51) {
              raw_res$keypoints[, seq(1, 51, by = 3)] <- raw_res$keypoints[, seq(1, 51, by = 3)] + crop_coords[1]
              raw_res$keypoints[, seq(2, 51, by = 3)] <- raw_res$keypoints[, seq(2, 51, by = 3)] + crop_coords[2]
            }
          } else if (roi_mode == "filter" && is.null(count_line) && !is.null(roi_coords) && nrow(raw_res$boxes) > 0) {
            cxs <- (raw_res$boxes[, 1] + raw_res$boxes[, 3]) / 2
            cys <- (raw_res$boxes[, 2] + raw_res$boxes[, 4]) / 2
            keep_roi <- which(cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            raw_res$boxes <- raw_res$boxes[keep_roi, , drop = FALSE]
            raw_res$scores <- raw_res$scores[keep_roi]
            raw_res$class_ids <- raw_res$class_ids[keep_roi]
            orig_idx <- orig_idx[keep_roi]
            if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
              raw_res$keypoints <- raw_res$keypoints[keep_roi, , drop = FALSE]
            }
          }

          num_boxes <- nrow(raw_res$boxes)
          if (num_boxes > 0) {
            c_ids <- raw_res$class_ids
            lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))
            df_b <- data.frame(
              frame = frame_count,
              timestamp = round(t_elapsed, 3),
              id = seq_len(num_boxes),
              xmin = raw_res$boxes[, 1],
              ymin = raw_res$boxes[, 2],
              xmax = raw_res$boxes[, 3],
              ymax = raw_res$boxes[, 4],
              label = lbls,
              score = round(raw_res$scores, 4),
              class_id = c_ids,
              orig_idx = orig_idx,
              stringsAsFactors = FALSE
            )
            if (!is.null(classes)) {
              keep_cls <- which(df_b$label %in% classes)
              df_b <- df_b[keep_cls, , drop = FALSE]
              if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
                raw_res$keypoints <- raw_res$keypoints[keep_cls, , drop = FALSE]
              }
            }

            if ((isTRUE(track) || !is.null(count_line))) {
              t_trk_start <- Sys.time()
              if (nrow(df_b) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = df_b$xmin,
                  ymin = df_b$ymin,
                  xmax = df_b$xmax,
                  ymax = df_b$ymax,
                  labels = df_b$label,
                  scores = df_b$score,
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history,
                  current_time = round(t_elapsed, 3),
                  frame_w = orig_w,
                  frame_h = orig_h,
                  class_ids = df_b$class_id,
                  mask_ids = if ("orig_idx" %in% names(df_b)) df_b$orig_idx else seq_len(nrow(df_b)),
                  interpolate_box = interpolate_box,
                  interpolate_mask = interpolate_mask,
                  min_hits = min_hits_val
                )
                track_state <- track_res$state
                total_count_val <- track_res$total_count
                in_zone_val <- track_res$in_zone
                n_trk <- length(track_res$track_ids)
                if (n_trk > 0) {
                  df_b <- data.frame(
                    frame = rep(frame_count, n_trk),
                    timestamp = rep(round(t_elapsed, 3), n_trk),
                    id = as.integer(track_res$track_ids),
                    xmin = as.numeric(track_res$xmin),
                    ymin = as.numeric(track_res$ymin),
                    xmax = as.numeric(track_res$xmax),
                    ymax = as.numeric(track_res$ymax),
                    label = as.character(track_res$labels),
                    score = as.numeric(track_res$scores),
                    class_id = as.integer(track_res$class_ids),
                    counted = as.logical(track_res$counted),
                    orig_idx = as.integer(track_res$mask_ids),
                    stringsAsFactors = FALSE
                  )
                  if (!is.null(track_res$in_zone_vec)) {
                    df_b$in_zone <- as.logical(track_res$in_zone_vec)
                  } else if (!is.null(roi_coords)) {
                    cxs <- (df_b$xmin + df_b$xmax) / 2
                    cys <- (df_b$ymin + df_b$ymax) / 2
                    df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                  }
                } else {
                  df_b <- df_b[0, , drop = FALSE]
                }
                track_ids_pass <- track_res$track_ids
                flash_pass <- track_res$flash
                history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
              } else if (!is.null(track_state$tracks) && length(track_state$tracks) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = numeric(0),
                  ymin = numeric(0),
                  xmax = numeric(0),
                  ymax = numeric(0),
                  labels = character(0),
                  scores = numeric(0),
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history,
                  current_time = round(t_elapsed, 3),
                  frame_w = orig_w,
                  frame_h = orig_h,
                  interpolate_box = interpolate_box,
                  interpolate_mask = interpolate_mask,
                  min_hits = min_hits_val
                )
                track_state <- track_res$state
                total_count_val <- track_res$total_count
                in_zone_val <- track_res$in_zone
                n_trk <- length(track_res$track_ids)
                if (n_trk > 0) {
                  df_b <- data.frame(
                    frame = rep(frame_count, n_trk),
                    timestamp = rep(round(t_elapsed, 3), n_trk),
                    id = as.integer(track_res$track_ids),
                    xmin = as.numeric(track_res$xmin),
                    ymin = as.numeric(track_res$ymin),
                    xmax = as.numeric(track_res$xmax),
                    ymax = as.numeric(track_res$ymax),
                    label = as.character(track_res$labels),
                    score = as.numeric(track_res$scores),
                    class_id = as.integer(track_res$class_ids),
                    counted = as.logical(track_res$counted),
                    orig_idx = as.integer(track_res$mask_ids),
                    stringsAsFactors = FALSE
                  )
                  if (!is.null(track_res$in_zone_vec)) {
                    df_b$in_zone <- as.logical(track_res$in_zone_vec)
                  } else if (!is.null(roi_coords)) {
                    cxs <- (df_b$xmin + df_b$xmax) / 2
                    cys <- (df_b$ymin + df_b$ymax) / 2
                    df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                  }
                  track_ids_pass <- track_res$track_ids
                  flash_pass <- track_res$flash
                  history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                  history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
                }
              }
              t_trk_end <- Sys.time()
              t_trk_ms <- as.numeric(t_trk_end - t_trk_start, units = "secs") * 1000
              last_track_res <- track_res
            }

            if (!("in_zone" %in% names(df_b)) && !is.null(roi_coords) && nrow(df_b) > 0) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
            }

            if (isTRUE(mask) && isTRUE(raw_res$is_seg) && nrow(df_b) > 0) {
              df_b <- .populate_seg_stats_df(df_b, raw_res, num_boxes_orig)
            }

            has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && nrow(raw_res$keypoints) == nrow(df_b)
            if (has_kpts) {
              kpt_mat_pass <- raw_res$keypoints
            }
          }

          if (save_crops && nrow(df_b) > 0L) {
            .save_object_crops_frame(bm, df_b, frame_count)
          }

          last_df_b <- df_b
          last_raw_res <- raw_res
          last_has_kpts <- has_kpts
          last_has_masks <- has_masks
          last_kpt_mat_pass <- kpt_mat_pass
          last_track_ids_pass <- track_ids_pass
          last_flash_pass <- flash_pass
          last_history_x_pass <- history_x_pass
          last_history_y_pass <- history_y_pass
          last_total_count_val <- total_count_val
          last_in_zone_val <- in_zone_val
          last_use_crop <- use_crop
        } else if (detect && !should_infer && !is.null(last_df_b)) {
          df_b <- last_df_b
          raw_res <- last_raw_res
          has_kpts <- last_has_kpts
          has_masks <- last_has_masks
          kpt_mat_pass <- last_kpt_mat_pass
          track_ids_pass <- last_track_ids_pass
          flash_pass <- last_flash_pass
          history_x_pass <- last_history_x_pass
          history_y_pass <- last_history_y_pass
          total_count_val <- last_total_count_val
          in_zone_val <- last_in_zone_val
          use_crop <- last_use_crop
          if (nrow(df_b) > 0) {
            df_b$frame <- frame_count
            df_b$timestamp <- round(t_elapsed, 3)
          }
        }

        palette_colors <- if (detect && nrow(df_b) > 0) {
          .get_class_palette(
            labels = df_b$label,
            n = nrow(df_b),
            user_col = col,
            default_col = if (!is.null(col)) col else "#00CC66",
            rainbow = rainbow,
            ids = if ("id" %in% names(df_b)) df_b$id else NULL
          )
        } else {
          character(0)
        }

        if (detect && return_data) {
          if (nrow(df_b) > 0) {
            df_rec <- df_b
            df_rec$orig_idx <- NULL
            df_rec$color <- palette_colors
            frame_records[[frame_count]] <- df_rec
          }
          if (has_kpts) {
            kpts_list_frame <- lapply(seq_len(nrow(raw_res$keypoints)), function(k) {
              vals <- raw_res$keypoints[k, ]
              data.frame(
                frame = frame_count,
                id = df_b$id[k],
                keypoint = .coco_keypoints,
                x = vals[seq(1, 51, by = 3)],
                y = vals[seq(2, 51, by = 3)],
                conf = round(vals[seq(3, 51, by = 3)], 4),
                stringsAsFactors = FALSE
              )
            })
            frame_keypoints[[frame_count]] <- do.call(rbind, kpts_list_frame)
          }
        }

        if (need_color_adj && (isTRUE(show) || isTRUE(save_video))) {
          adjust_frame_color_cpp(bm, saturation = saturation_val, contrast = contrast_val)
        }

        if (heatmap && !is.null(heatmap_accum)) {
          update_video_heatmap_cpp(
            heatmap = heatmap_accum,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            decay = heatmap_decay,
            draw_on_bm = isTRUE(show) || isTRUE(save_video),
            bm = bm,
            alpha = heatmap_alpha
          )
        }

        t_lat_ms <- max(1.0, t_inf_ms + t_trk_ms + last_gui_ms)
        if (telemetry != "none" && f_scale_hud > 0.05) {
          trk_info <- if (!is.null(last_track_res)) last_track_res else list()
          trk_active <- if (!is.null(trk_info$active)) trk_info$active else in_zone_val
          conf_mean <- if (!is.null(trk_info$conf_mean)) trk_info$conf_mean else if (nrow(df_b) > 0) mean(df_b$score) else 0.0
          conf_sd <- if (!is.null(trk_info$conf_sd)) trk_info$conf_sd else 0.0
          conf_min <- if (!is.null(trk_info$conf_min)) trk_info$conf_min else if (nrow(df_b) > 0) min(df_b$score) else 0.0
          flow_val <- if (!is.null(trk_info$flow)) trk_info$flow else 0.0
          vel_mean <- if (!is.null(trk_info$vel_mean)) trk_info$vel_mean else 0.0
          vel_jit <- if (!is.null(trk_info$vel_jit)) trk_info$vel_jit else 0.0
          welford_mean <- if (!is.null(trk_info$welford_mean)) trk_info$welford_mean else 0.0
          welford_cv <- if (!is.null(trk_info$welford_cv)) trk_info$welford_cv else 0.0
          roi_occ <- if (!is.null(trk_info$roi_occ)) trk_info$roi_occ else 0.0

          raw_stats <- list(
            fps = cur_fps,
            lat = t_lat_ms,
            inf = t_inf_ms,
            trk = t_trk_ms,
            gui = last_gui_ms,
            conf_mean = conf_mean,
            conf_sd = conf_sd,
            conf_min = conf_min,
            flow = flow_val,
            vel_mean = vel_mean,
            vel_jit = vel_jit,
            size_mean = welford_mean,
            size_cv = welford_cv,
            roi_occ = roi_occ
          )
          sm_stats <- stats_smoother(raw_stats)

          if (is_detailed) {
            hud_msg <- paste(
              sprintf("FRAME: %d", frame_count),
              sprintf("SPEED: %.1f FPS", sm_stats$fps),
              sprintf("LATENCY: %dms", round(sm_stats$lat)),
              sprintf("  INF: %dms", round(sm_stats$inf)),
              sprintf("  TRK: %dms", round(sm_stats$trk)),
              sprintf("  GUI: %dms", round(sm_stats$gui)),
              sprintf("TRACKS: %d ACT", trk_active),
              sprintf("CONF: %.2f \u00b1 %.2f (MIN: %.2f)", sm_stats$conf_mean, sm_stats$conf_sd, sm_stats$conf_min),
              sprintf("FLOW: %.1f obj/s", sm_stats$flow),
              sprintf("MOTION: %.1f px/f (JIT: %.1f)", sm_stats$vel_mean, sm_stats$vel_jit),
              sprintf("IN ZONE: %d", in_zone_val),
              sprintf("COUNTED: %d", total_count_val),
              sprintf("SIZE: %.0f px\u00b2 (CV: %.1f%%)", sm_stats$size_mean, sm_stats$size_cv),
              sprintf("ROI OCC: %.1f%%", sm_stats$roi_occ),
              sep = "\n"
            )
            sparkline_pass <- if (!is.null(trk_info$flow_history)) as.numeric(trk_info$flow_history) else NULL
          } else {
            hud_msg <- if (isTRUE(track) || !is.null(count_line)) {
              sprintf("%.1f FPS | IN ZONE: %d | COUNT: %d", sm_stats$fps, in_zone_val, total_count_val)
            } else if (detect) {
              sprintf("%.1f FPS | %d Dets", sm_stats$fps, nrow(df_b))
            } else {
              sprintf("%.1f FPS", sm_stats$fps)
            }
            sparkline_pass <- NULL
          }
        } else {
          hud_msg <- ""
          sparkline_pass <- NULL
        }

        if ((detect && nrow(df_b) > 0) || nzchar(hud_msg) || !is.null(roi_coords) || !is.null(line_coords)) {
          t_gui_start <- Sys.time()
          draw_yolo_detections_bgr_cpp(
            bm = bm,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            labels = if (detect && nrow(df_b) > 0) df_b$label else character(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            colors = palette_colors,
            lwd = as.integer(lwd),
            hud_text = hud_msg,
            font_scale = as.numeric(f_scale_obj),
            keypoints = kpt_mat_pass,
            kpt_threshold = as.numeric(kpt_threshold),
            kpt_radius = as.integer(kpt_radius),
            draw_skeleton = isTRUE(skeleton),
            draw_boxes = isTRUE(show_box),
            track_ids = track_ids_pass,
            flash = flash_pass,
            history_x = history_x_pass,
            history_y = history_y_pass,
            roi = roi_coords,
            roi_label = roi_label,
            count_line = line_coords,
            count_line_label = if (!is.null(count_line_label)) as.character(count_line_label) else "",
            counted = if ("counted" %in% names(df_b)) df_b$counted else NULL,
            hide_outside_roi = isTRUE(hide_outside_roi),
            hide_counted = isTRUE(hide_counted),
            flash_counted = isTRUE(flash_counted),
            show_text = isTRUE(show_text),
            show_conf = isTRUE(show_conf),
            show_class = isTRUE(show_class),
            show_id = isTRUE(show_id),
            badge = badge,
            pad_scale = pad_num,
            mask_labels = if (isTRUE(mask) && has_masks) raw_res$labels else NULL,
            mask_alpha = as.numeric(alpha),
            mask_offset_x = if (use_crop) as.integer(crop_coords[1]) else 0L,
            mask_offset_y = if (use_crop) as.integer(crop_coords[2]) else 0L,
            draw_masks = isTRUE(mask),
            mask_ids = if (has_masks && nrow(df_b) > 0) as.integer(df_b$orig_idx) else NULL,
            hud_pos = as.character(hud_pos),
            hud_layout = as.character(hud_layout),
            font_scale_hud = as.numeric(f_scale_hud),
            font_scale_roi = as.numeric(roi_sz),
            detailed_telemetry = is_detailed,
            sparkline = sparkline_pass,
            show_point = isTRUE(show_point),
            point_radius = as.integer(point_radius),
            point_col = if (!is.null(point_col)) as.character(point_col) else ""
          )
          t_gui_end <- Sys.time()
          last_gui_ms <- max(1.0, as.numeric(t_gui_end - t_gui_start, units = "secs") * 1000)
        }

        if (save_video) {
          if (!is.null(cam_writer)) {
            write_native_video_frame_cpp(cam_writer, bm)
          } else {
            f_rec_path <- file.path(cam_dir, sprintf("rec_%06d.jpg", frame_count))
            save_frame_bgr_cpp(bm, f_rec_path, quality = 90L)
            recorded_frames[[frame_count]] <- f_rec_path
          }
        }

        if (isTRUE(show)) {
          show_native_frame_cpp(win, bm, preserve_aspect = TRUE)
          status <- poll_native_window_cpp(win)
          if (status == 1L) {
            break
          }
          while (status == 2L) {
            Sys.sleep(0.05)
            status <- poll_native_window_cpp(win)
            if (status == 1L) break
          }
          if (status == 1L) break
        }
      }

    valid_records <- frame_records[!sapply(frame_records, is.null)]
    valid_kpts <- frame_keypoints[!sapply(frame_keypoints, is.null)]
    all_dets <- if (length(valid_records) > 0) do.call(rbind, valid_records) else data.frame()
    all_kpts <- if (length(valid_kpts) > 0) do.call(rbind, valid_kpts) else NULL
    total_time <- as.numeric(Sys.time() - t_start, units = "secs")
    avg_fps <- if (total_time > 0 && frame_count > 0) (frame_count / total_time) else fps_smooth

    if (save_video) {
      if (!is.null(cam_writer)) {
        try(close_native_video_writer_cpp(cam_writer), silent = TRUE)
        if (isTRUE(verbose)) cli::cli_alert_success("Camera video saved: {.path {out_video_file}}")
      } else if (length(recorded_frames) > 0) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Encoding recorded camera session to {.path {out_video_file}}...")
        }
        rec_fps <- if (!is.null(fps)) fps else (if (avg_fps > 0) round(avg_fps, 1) else 24.0)
        av::av_encode_video(
          input = unlist(recorded_frames),
          output = out_video_file,
          framerate = rec_fps,
          verbose = FALSE
        )
        if (isTRUE(verbose)) cli::cli_alert_success("Camera video saved: {.path {out_video_file}}")
      }
    }

  # ============================================================================
  # MODE B: VIDEO FILE PROCESSING (AV / FFMPEG)
  # ============================================================================
  } else {
    if (!file.exists(video)) {
      cli::cli_abort("Video file not found: {.path {video}}")
    }
    ext <- tolower(tools::file_ext(video))
    if (ext %in% c("jpg", "jpeg", "png", "bmp", "webp", "tif", "tiff")) {
      cli::cli_abort(c(
        "x" = "File {.path {video}} is a static image, not a video file.",
        "i" = "For static image detection or pose estimation, please use {.code image_detect_dl()} or {.code image_pose_dl()}."
      ))
    }
    if (!requireNamespace("av", quietly = TRUE)) {
      cli::cli_abort(c(
        "!" = "Package {.pkg av} is required for video file inference.",
        "i" = "Please install it with: {.code install.packages('av')}"
      ))
    }

    out_video_file <- if (save_video) normalizePath(output, winslash = "/", mustWork = FALSE) else NULL

    use_native_stream <- (.Platform$OS.type == "windows" && has_native_video_reader_cpp())
    v_reader <- if (use_native_stream) {
      tryCatch(open_native_video_reader_cpp(normalizePath(video, winslash = "/", mustWork = TRUE)),
               error = function(e) NULL)
    } else NULL

    use_native_writer <- FALSE
    v_writer <- NULL
    v_info <- NULL
    v_info_native <- NULL

    if (!is.null(v_reader)) {
      use_native_stream <- TRUE
      v_info_native <- get_video_reader_info_cpp(v_reader)
      src_w <- v_info_native$width
      src_h <- v_info_native$height
      orig_fps <- v_info_native$fps
      if (is.null(orig_fps) || is.na(orig_fps) || orig_fps <= 0) orig_fps <- 30.0
      src_fps <- if (!is.null(fps)) fps else orig_fps
      if (is.null(src_fps) || is.na(src_fps) || src_fps <= 0) src_fps <- 30.0

      eff_dur <- v_info_native$duration
      if (is.null(eff_dur) || is.na(eff_dur) || eff_dur <= 0) {
        eff_dur <- v_info_native$total_frames / orig_fps
      }
      if (end_sec > 0 && end_sec < eff_dur) eff_dur <- end_sec
      if (start_sec > 0 && start_sec < eff_dur) eff_dur <- eff_dur - start_sec
      est_total <- if (!is.null(fps) && fps > 0 && fps < orig_fps * 0.98) {
        max(1L, as.integer(ceiling(eff_dur * fps)))
      } else if (start_sec > 0 || end_sec > 0) {
        max(1L, as.integer(ceiling(eff_dur * orig_fps)))
      } else {
        as.integer(v_info_native$total_frames)
      }
      num_frames <- if (!is.null(max_frames) && max_frames < est_total) as.integer(max_frames) else est_total

      start_native_video_prefetch_cpp(
        v_reader,
        target_fps = if (!is.null(fps) && fps > 0) as.numeric(fps) else 0.0,
        start_time = start_sec,
        end_time = end_sec
      )
      on.exit(try(close_native_video_reader_cpp(v_reader), silent = TRUE), add = TRUE)

      if (save_video && has_native_video_writer_cpp()) {
        v_writer <- tryCatch(
          open_native_video_writer_cpp(
            video_path = out_video_file,
            width = src_w,
            height = src_h,
            fps = src_fps
          ),
          error = function(e) NULL
        )
        if (!is.null(v_writer)) {
          use_native_writer <- TRUE
          on.exit(try(close_native_video_writer_cpp(v_writer), silent = TRUE), add = TRUE)
        }
      }

      run_dir <- if (save_video && !use_native_writer) file.path(tempdir(), paste0("vid_yolo_", as.integer(stats::runif(1, 1000, 9999)))) else NULL
      if (!is.null(run_dir)) {
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
        on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
      }
      frame_paths <- character(0)
    } else {
      use_native_stream <- FALSE
      v_info <- suppressWarnings(av::av_video_info(video))
      src_w <- v_info$video$width
      src_h <- v_info$video$height
      orig_fps <- v_info$video$framerate
      if (is.null(orig_fps) || is.na(orig_fps) || orig_fps <= 0) orig_fps <- 30.0
      src_fps <- if (!is.null(fps)) fps else orig_fps
      if (is.null(src_fps) || is.na(src_fps) || src_fps <= 0) src_fps <- 30.0
      # Destination directory for decoded frames
      run_dir <- file.path(tempdir(), paste0("vid_yolo_", as.integer(stats::runif(1, 1000, 9999))))
      on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)

      if (isTRUE(verbose)) {
        cli::cli_progress_step("Decoding video frames with FFmpeg...", msg_done = "Video frames decoded")
      }
      filter_fps <- if (length(fps) && !is.null(fps)) paste0("fps=fps=", fps) else NULL
      vfilter <- if (length(filter_fps)) filter_fps else "null"
      if (!is.null(max_frames) && is.numeric(max_frames) && max_frames > 0) {
        lim_n <- as.integer(max_frames) + 2L
        sel_filter <- paste0("select=between(n\\,0\\,", lim_n, ")")
        vfilter <- if (vfilter == "null") sel_filter else paste0(vfilter, ",", sel_filter)
      }
      if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      output_pattern <- file.path(run_dir, "image_%06d.bmp")

      av::av_encode_video(
        input = video,
        output = output_pattern,
        framerate = orig_fps,
        codec = "bmp",
        vfilter = vfilter,
        verbose = FALSE
      )
      if (isTRUE(verbose)) cli::cli_progress_done()
      frame_paths <- sort(list.files(run_dir, pattern = "^image_\\d+\\.bmp$", full.names = TRUE))
      num_frames <- length(frame_paths)
      if (!is.null(fps) && !is.null(v_info$video$duration) && v_info$video$duration > 0) {
        expected_max <- ceiling(v_info$video$duration * fps) + 5
        if (num_frames > expected_max) {
          keep_idx <- round(seq(1, num_frames, length.out = ceiling(v_info$video$duration * fps)))
          keep_idx <- unique(pmin(num_frames, pmax(1L, keep_idx)))
          frame_paths <- frame_paths[keep_idx]
          num_frames <- length(frame_paths)
        }
      }
      if (!is.null(max_frames) && max_frames < num_frames) {
        frame_paths <- frame_paths[seq_len(max_frames)]
        num_frames <- length(frame_paths)
      }
    }

    if (isTRUE(verbose)) {
      if (detect) {
        cli::cli_alert_info("Running YOLO inference on {num_frames} frame(s) [{toupper(engine)} | model: {.val {basename(model_file)}}]")
      } else {
        cli::cli_alert_info("Processing {num_frames} frame(s) [inference disabled]")
      }
      vid_duration <- if (!is.null(v_info_native)) v_info_native$duration else if (!is.null(v_info)) v_info$video$duration else 0
      if (is.null(fps) && isTRUE(vid_duration > 60)) {
        cli::cli_alert_info("Tip: For faster processing of long videos, you can use e.g. {.code fps = 5} or {.code fps = 10}.")
      }
    }

    frame_records <- vector("list", num_frames)
    frame_keypoints <- vector("list", num_frames)
    t_start <- Sys.time()

    # Pre-render dimensions must be even for H.264
    out_w <- src_w - (src_w %% 2)
    out_h <- src_h - (src_h %% 2)
    heatmap_accum <- if (heatmap) matrix(0.0, nrow = src_w, ncol = src_h) else NULL
    need_masks_vid <- isTRUE(mask) && (isTRUE(save_video) || isTRUE(show))
    use_native_win_b <- isTRUE(show) && (.Platform$OS.type == "windows") && (backend == "native")
    win_b <- NULL
    if (use_native_win_b) {
      w_w <- if (!is.null(window_size)) window_size[1] else round(src_w * window_scale)
      w_h <- if (!is.null(window_size)) window_size[2] else round(src_h * window_scale)
      win_title <- sprintf("pliman - YOLO Video [%s | %s] (ESC to exit, Space to pause, F for fullscreen)", toupper(engine), basename(model_file))
      win_b <- tryCatch(
        create_native_window_cpp(title = win_title, width = as.integer(w_w), height = as.integer(w_h), fullscreen = isTRUE(fullscreen)),
        error = function(e) NULL
      )
      if (!is.null(win_b)) {
        on.exit(try(close_native_window_cpp(win_b), silent = TRUE), add = TRUE)
      } else {
        use_native_win_b <- FALSE
      }
    }

    # ROI coordinates for Mode B
    roi_coords <- .resolve_roi(roi, src_w, src_h, roi_orientation)

    line_coords <- .resolve_count_line(count_line, src_w, src_h, count_line_orientation)

    if (save_video) {
      out_video_file <- normalizePath(output, winslash = "/", mustWork = FALSE)
    }

    pb <- if (isTRUE(verbose)) {
      cli::cli_progress_bar(
        if (detect) "Detecting objects in video" else "Processing video",
        total = num_frames,
        extra = list(current_fps = "0.0"),
        format = "{cli::pb_spin} Processing [{cli::pb_current}/{cli::pb_total}] | {cli::pb_percent} | {cli::pb_extra$current_fps} FPS | ETA: {cli::pb_eta}"
      )
    } else NULL
    pb_step <- max(1L, min(5L, as.integer(ceiling(num_frames / 100))))

    last_df_b <- NULL
    last_raw_res <- NULL
    last_has_kpts <- FALSE
    last_has_masks <- FALSE
    last_kpt_mat_pass <- NULL
    last_track_ids_pass <- NULL
    last_flash_pass <- NULL
    last_history_x_pass <- NULL
    last_history_y_pass <- NULL
    last_total_count_val <- 0L
    last_in_zone_val <- 0L
    last_use_crop <- FALSE
    last_crop_coords <- NULL
    last_track_res <- NULL
    t_inf_ms <- 0.0
    t_trk_ms <- 0.0
    last_gui_ms <- 14.0
    last_frame_t <- Sys.time()
    stats_smoother <- .make_stats_smoother(smoot_stat)

    recorded_frames <- vector("list", num_frames)
    frames_processed <- 0L
    bm_stream <- if (use_native_stream) {
      b_raw <- raw(3L * src_w * src_h)
      attr(b_raw, "dim") <- c(3L, src_w, src_h)
      b_raw
    } else NULL

    process_loop <- function() {
      for (i in seq_len(num_frames)) {
        if (use_native_stream) {
          has_frame <- read_native_video_frame_cpp(v_reader, bm_stream)
          if (!has_frame) break
          bm <- bm_stream
          f_path <- ""
        } else {
          f_path <- frame_paths[i]
          bm <- load_frame_bgr_cpp(f_path)
        }
        has_kpts <- FALSE
        has_masks <- FALSE
        kpt_mat_pass <- NULL
        track_ids_pass <- NULL
        flash_pass <- NULL
        history_x_pass <- NULL
        history_y_pass <- NULL

        use_crop <- FALSE
        crop_coords <- c(0L, 0L, src_w, src_h)
        should_infer <- (i == 1L || (i %% infer_every == 0L))
        t_sec <- if (use_native_stream) round(get_video_reader_timestamp_cpp(v_reader), 3) else round((i - 1) / src_fps, 3)

        if (detect && should_infer) {
          use_crop <- (roi_mode == "crop" && !is.null(roi_coords))
          if (use_crop) {
            margin_x <- round(0.08 * src_w)
            margin_y <- round(0.08 * src_h)
            crop_x1 <- max(0L, as.integer(roi_coords[1] - margin_x))
            crop_y1 <- max(0L, as.integer(roi_coords[2] - margin_y))
            crop_x2 <- min(src_w, as.integer(roi_coords[3] + margin_x))
            crop_y2 <- min(src_h, as.integer(roi_coords[4] + margin_y))
            crop_coords <- c(crop_x1, crop_y1, crop_x2, crop_y2)
            bm_feed <- crop_bgr_cpp(bm, crop_x1, crop_y1, crop_x2, crop_y2)
          } else {
            crop_coords <- c(0L, 0L, src_w, src_h)
            bm_feed <- bm
          }
          feed_dims <- dim(bm_feed)
          feed_w <- feed_dims[2]
          feed_h <- feed_dims[3]

          t_inf_start <- Sys.time()
          raw_res <- .run_video_yolo(bm_feed, feed_w, feed_h, need_masks_vid)
          t_inf_end <- Sys.time()
          t_inf_ms <<- as.numeric(t_inf_end - t_inf_start, units = "secs") * 1000

          has_masks <- need_masks_vid && isTRUE(raw_res$is_seg) && !is.null(raw_res$labels) && is.matrix(raw_res$labels) && nrow(raw_res$boxes) > 0L
          num_boxes_orig <- nrow(raw_res$boxes)
          orig_idx <- seq_len(num_boxes_orig)

          # Shift box coordinates back to full frame if crop was used
          if (use_crop && nrow(raw_res$boxes) > 0) {
            raw_res$boxes[, 1] <- raw_res$boxes[, 1] + crop_coords[1]
            raw_res$boxes[, 2] <- raw_res$boxes[, 2] + crop_coords[2]
            raw_res$boxes[, 3] <- raw_res$boxes[, 3] + crop_coords[1]
            raw_res$boxes[, 4] <- raw_res$boxes[, 4] + crop_coords[2]
            if (!is.null(raw_res$keypoints) && ncol(raw_res$keypoints) >= 51) {
              raw_res$keypoints[, seq(1, 51, by = 3)] <- raw_res$keypoints[, seq(1, 51, by = 3)] + crop_coords[1]
              raw_res$keypoints[, seq(2, 51, by = 3)] <- raw_res$keypoints[, seq(2, 51, by = 3)] + crop_coords[2]
            }
          } else if (roi_mode == "filter" && is.null(count_line) && !is.null(roi_coords) && nrow(raw_res$boxes) > 0) {
            cxs <- (raw_res$boxes[, 1] + raw_res$boxes[, 3]) / 2
            cys <- (raw_res$boxes[, 2] + raw_res$boxes[, 4]) / 2
            keep_roi <- which(cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            raw_res$boxes <- raw_res$boxes[keep_roi, , drop = FALSE]
            raw_res$scores <- raw_res$scores[keep_roi]
            raw_res$class_ids <- raw_res$class_ids[keep_roi]
            orig_idx <- orig_idx[keep_roi]
            if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
              raw_res$keypoints <- raw_res$keypoints[keep_roi, , drop = FALSE]
            }
          }

          num_boxes <- nrow(raw_res$boxes)

          if (num_boxes > 0) {
            c_ids <- raw_res$class_ids
            lbls <- ifelse(c_ids >= 0 & c_ids < length(class_names), class_names[c_ids + 1], as.character(c_ids))
            df_b <- data.frame(
              frame = i,
              timestamp = t_sec,
              id = seq_len(num_boxes),
              xmin = raw_res$boxes[, 1],
              ymin = raw_res$boxes[, 2],
              xmax = raw_res$boxes[, 3],
              ymax = raw_res$boxes[, 4],
              label = lbls,
              score = round(raw_res$scores, 4),
              class_id = c_ids,
              orig_idx = orig_idx,
              stringsAsFactors = FALSE
            )
            if (!is.null(classes)) {
              keep_cls <- which(df_b$label %in% classes)
              df_b <- df_b[keep_cls, , drop = FALSE]
              if (!is.null(raw_res$keypoints) && nrow(raw_res$keypoints) > 0) {
                raw_res$keypoints <- raw_res$keypoints[keep_cls, , drop = FALSE]
              }
            }

            if ((isTRUE(track) || !is.null(count_line))) {
              t_trk_start <- Sys.time()
              if (nrow(df_b) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = df_b$xmin,
                  ymin = df_b$ymin,
                  xmax = df_b$xmax,
                  ymax = df_b$ymax,
                  labels = df_b$label,
                  scores = df_b$score,
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history,
                  current_time = t_sec,
                  frame_w = src_w,
                  frame_h = src_h,
                  class_ids = df_b$class_id,
                  mask_ids = if ("orig_idx" %in% names(df_b)) df_b$orig_idx else seq_len(nrow(df_b)),
                  interpolate_box = interpolate_box,
                  interpolate_mask = interpolate_mask,
                  min_hits = min_hits_val
                )
                track_state <<- track_res$state
                total_count_val <<- track_res$total_count
                in_zone_val <<- track_res$in_zone
                n_trk <- length(track_res$track_ids)
                if (n_trk > 0) {
                  df_b <- data.frame(
                    frame = rep(i, n_trk),
                    timestamp = rep(t_sec, n_trk),
                    id = as.integer(track_res$track_ids),
                    xmin = as.numeric(track_res$xmin),
                    ymin = as.numeric(track_res$ymin),
                    xmax = as.numeric(track_res$xmax),
                    ymax = as.numeric(track_res$ymax),
                    label = as.character(track_res$labels),
                    score = as.numeric(track_res$scores),
                    class_id = as.integer(track_res$class_ids),
                    counted = as.logical(track_res$counted),
                    orig_idx = as.integer(track_res$mask_ids),
                    stringsAsFactors = FALSE
                  )
                  if (!is.null(track_res$in_zone_vec)) {
                    df_b$in_zone <- as.logical(track_res$in_zone_vec)
                  } else if (!is.null(roi_coords)) {
                    cxs <- (df_b$xmin + df_b$xmax) / 2
                    cys <- (df_b$ymin + df_b$ymax) / 2
                    df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                  }
                } else {
                  df_b <- df_b[0, , drop = FALSE]
                }
                track_ids_pass <- track_res$track_ids
                flash_pass <- track_res$flash
                history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
              } else if (!is.null(track_state$tracks) && length(track_state$tracks) > 0) {
                track_res <- update_tracker_cpp(
                  xmin = numeric(0),
                  ymin = numeric(0),
                  xmax = numeric(0),
                  ymax = numeric(0),
                  labels = character(0),
                  scores = numeric(0),
                  state = track_state,
                  count_line = line_coords,
                  roi = roi_coords,
                  max_dist = track_max_dist,
                  min_iou = track_min_iou,
                  max_lost = track_max_lost,
                  max_history = track_max_history,
                  current_time = t_sec,
                  frame_w = src_w,
                  frame_h = src_h,
                  interpolate_box = interpolate_box,
                  interpolate_mask = interpolate_mask,
                  min_hits = min_hits_val
                )
                track_state <<- track_res$state
                total_count_val <<- track_res$total_count
                in_zone_val <<- track_res$in_zone
                n_trk <- length(track_res$track_ids)
                if (n_trk > 0) {
                  df_b <- data.frame(
                    frame = rep(i, n_trk),
                    timestamp = rep(t_sec, n_trk),
                    id = as.integer(track_res$track_ids),
                    xmin = as.numeric(track_res$xmin),
                    ymin = as.numeric(track_res$ymin),
                    xmax = as.numeric(track_res$xmax),
                    ymax = as.numeric(track_res$ymax),
                    label = as.character(track_res$labels),
                    score = as.numeric(track_res$scores),
                    class_id = as.integer(track_res$class_ids),
                    counted = as.logical(track_res$counted),
                    orig_idx = as.integer(track_res$mask_ids),
                    stringsAsFactors = FALSE
                  )
                  if (!is.null(track_res$in_zone_vec)) {
                    df_b$in_zone <- as.logical(track_res$in_zone_vec)
                  } else if (!is.null(roi_coords)) {
                    cxs <- (df_b$xmin + df_b$xmax) / 2
                    cys <- (df_b$ymin + df_b$ymax) / 2
                    df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
                  }
                  track_ids_pass <- track_res$track_ids
                  flash_pass <- track_res$flash
                  history_x_pass <- if (isTRUE(trail)) track_res$history_x else NULL
                  history_y_pass <- if (isTRUE(trail)) track_res$history_y else NULL
                }
              }
              t_trk_end <- Sys.time()
              t_trk_ms <<- as.numeric(t_trk_end - t_trk_start, units = "secs") * 1000
              last_track_res <<- track_res
            }

            if (!("in_zone" %in% names(df_b)) && !is.null(roi_coords) && nrow(df_b) > 0) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              df_b$in_zone <- cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4]
            }

            if (isTRUE(mask) && isTRUE(raw_res$is_seg) && nrow(df_b) > 0) {
              df_b <- .populate_seg_stats_df(df_b, raw_res, num_boxes_orig)
            }

            has_kpts <- !is.null(raw_res$keypoints) && is.matrix(raw_res$keypoints) && ncol(raw_res$keypoints) == 51 && nrow(raw_res$keypoints) == nrow(df_b)
            if (has_kpts) {
              kpt_mat_pass <- raw_res$keypoints
            }
          } else {
            df_b <- data.frame(
              frame = integer(0), timestamp = numeric(0), id = integer(0),
              xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
              label = character(0), score = numeric(0), class_id = integer(0),
              stringsAsFactors = FALSE
            )
          }

          if (save_crops && nrow(df_b) > 0L) {
            .save_object_crops_frame(bm, df_b, i)
          }

          last_df_b <<- df_b
          last_raw_res <<- raw_res
          last_has_kpts <<- has_kpts
          last_has_masks <<- has_masks
          last_kpt_mat_pass <<- kpt_mat_pass
          last_track_ids_pass <<- track_ids_pass
          last_flash_pass <<- flash_pass
          last_history_x_pass <<- history_x_pass
          last_history_y_pass <<- history_y_pass
          last_total_count_val <<- total_count_val
          last_in_zone_val <<- in_zone_val
          last_use_crop <<- use_crop
          last_crop_coords <<- crop_coords
        } else if (detect && !should_infer && !is.null(last_df_b)) {
          df_b <- last_df_b
          raw_res <- last_raw_res
          has_kpts <- last_has_kpts
          has_masks <- last_has_masks
          kpt_mat_pass <- last_kpt_mat_pass
          track_ids_pass <- last_track_ids_pass
          flash_pass <- last_flash_pass
          history_x_pass <- last_history_x_pass
          history_y_pass <- last_history_y_pass
          total_count_val <<- last_total_count_val
          in_zone_val <<- last_in_zone_val
          use_crop <- last_use_crop
          crop_coords <- last_crop_coords

          # If tracking is active, propagate motion across skipped frames using track velocities
          if ((isTRUE(track) || !is.null(count_line)) && nrow(df_b) > 0 && !is.null(track_state$tracks) && length(track_state$tracks) > 0) {
            for (ti in seq_along(track_state$tracks)) {
              trk <- track_state$tracks[[ti]]
              tid <- trk$id
              match_row <- which(df_b$id == tid)
              if (length(match_row) > 0) {
                vx <- if (!is.null(trk$vx)) trk$vx else 0.0
                vy <- if (!is.null(trk$vy)) trk$vy else 0.0
                if (abs(vx) > 0.1 || abs(vy) > 0.1) {
                  df_b$xmin[match_row] <- df_b$xmin[match_row] + vx
                  df_b$xmax[match_row] <- df_b$xmax[match_row] + vx
                  df_b$ymin[match_row] <- df_b$ymin[match_row] + vy
                  df_b$ymax[match_row] <- df_b$ymax[match_row] + vy
                }
              }
            }
            last_df_b <<- df_b
          }

          if (nrow(df_b) > 0) {
            df_b$frame <- i
            df_b$timestamp <- t_sec
          }
        } else {
          df_b <- data.frame(
            frame = integer(0), timestamp = numeric(0), id = integer(0),
            xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
            label = character(0), score = numeric(0), class_id = integer(0),
            stringsAsFactors = FALSE
          )
        }

        palette_colors <- if (detect && nrow(df_b) > 0) {
          .get_class_palette(
            labels = df_b$label,
            n = nrow(df_b),
            user_col = col,
            default_col = if (!is.null(col)) col else "#00CC66",
            rainbow = rainbow,
            ids = if ("id" %in% names(df_b)) df_b$id else NULL
          )
        } else {
          character(0)
        }

        if (detect && return_data) {
          if (nrow(df_b) > 0) {
            df_rec <- df_b
            df_rec$orig_idx <- NULL
            df_rec$color <- palette_colors
            frame_records[[i]] <<- df_rec
          }
          if (has_kpts) {
            kpts_list_frame <- lapply(seq_len(nrow(df_b)), function(k) {
              vals <- raw_res$keypoints[k, ]
              data.frame(
                frame = i,
                id = df_b$id[k],
                keypoint = .coco_keypoints,
                x = vals[seq(1, 51, by = 3)],
                y = vals[seq(2, 51, by = 3)],
                conf = round(vals[seq(3, 51, by = 3)], 4),
                stringsAsFactors = FALSE
              )
            })
            frame_keypoints[[i]] <<- do.call(rbind, kpts_list_frame)
          }
        }

        if (need_color_adj && (isTRUE(show) || isTRUE(save_video))) {
          adjust_frame_color_cpp(bm, saturation = saturation_val, contrast = contrast_val)
        }

        if (heatmap && !is.null(heatmap_accum)) {
          update_video_heatmap_cpp(
            heatmap = heatmap_accum,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            decay = heatmap_decay,
            draw_on_bm = isTRUE(save_video) || use_native_win_b,
            bm = bm,
            alpha = heatmap_alpha
          )
        }

        # Render frame with bounding boxes and annotations
        if (save_video || use_native_win_b) {
          elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
          cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
          t_now <- Sys.time()
          dt_frame <- as.numeric(t_now - last_frame_t, units = "secs")
          last_frame_t <<- t_now
          cur_fps_inst <- if (dt_frame > 0 && dt_frame < 10) (1.0 / dt_frame) else cur_fps_est
          t_lat_ms <- max(1.0, t_inf_ms + t_trk_ms + last_gui_ms)

          if (telemetry != "none" && f_scale_hud > 0.05) {
            trk_info <- if (!is.null(last_track_res)) last_track_res else list()
            trk_active <- if (!is.null(trk_info$active)) trk_info$active else in_zone_val
            conf_mean <- if (!is.null(trk_info$conf_mean)) trk_info$conf_mean else if (nrow(df_b) > 0) mean(df_b$score) else 0.0
            conf_sd <- if (!is.null(trk_info$conf_sd)) trk_info$conf_sd else 0.0
            conf_min <- if (!is.null(trk_info$conf_min)) trk_info$conf_min else if (nrow(df_b) > 0) min(df_b$score) else 0.0
            flow_val <- if (!is.null(trk_info$flow)) trk_info$flow else 0.0
            vel_mean <- if (!is.null(trk_info$vel_mean)) trk_info$vel_mean else 0.0
            vel_jit <- if (!is.null(trk_info$vel_jit)) trk_info$vel_jit else 0.0
            welford_mean <- if (!is.null(trk_info$welford_mean)) trk_info$welford_mean else 0.0
            welford_cv <- if (!is.null(trk_info$welford_cv)) trk_info$welford_cv else 0.0
            roi_occ <- if (!is.null(trk_info$roi_occ)) trk_info$roi_occ else 0.0

            raw_stats <- list(
              fps = cur_fps_inst,
              lat = t_lat_ms,
              inf = t_inf_ms,
              trk = t_trk_ms,
              gui = last_gui_ms,
              conf_mean = conf_mean,
              conf_sd = conf_sd,
              conf_min = conf_min,
              flow = flow_val,
              vel_mean = vel_mean,
              vel_jit = vel_jit,
              size_mean = welford_mean,
              size_cv = welford_cv,
              roi_occ = roi_occ
            )
            sm_stats <- stats_smoother(raw_stats)

            if (is_detailed) {
              hud_txt <- paste(
                sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", sm_stats$fps),
                sprintf("LATENCY: %dms", round(sm_stats$lat)),
                sprintf("  INF: %dms", round(sm_stats$inf)),
                sprintf("  TRK: %dms", round(sm_stats$trk)),
                sprintf("  GUI: %dms", round(sm_stats$gui)),
                sprintf("TRACKS: %d ACT", trk_active),
                sprintf("CONF: %.2f \u00b1 %.2f (MIN: %.2f)", sm_stats$conf_mean, sm_stats$conf_sd, sm_stats$conf_min),
                sprintf("FLOW: %.1f obj/s", sm_stats$flow),
                sprintf("MOTION: %.1f px/f (JIT: %.1f)", sm_stats$vel_mean, sm_stats$vel_jit),
                sprintf("IN ZONE: %d", in_zone_val),
                sprintf("COUNTED: %d", total_count_val),
                sprintf("SIZE: %.0f px\u00b2 (CV: %.1f%%)", sm_stats$size_mean, sm_stats$size_cv),
                sprintf("ROI OCC: %.1f%%", sm_stats$roi_occ),
                sep = "\n"
              )
              sparkline_pass <- if (!is.null(trk_info$flow_history)) as.numeric(trk_info$flow_history) else NULL
            } else {
              hud_txt <- if (isTRUE(track) || !is.null(count_line)) {
                sprintf("Frame %d/%d | %.1f FPS | IN ZONE: %d | COUNT: %d", i, num_frames, sm_stats$fps, in_zone_val, total_count_val)
              } else if (detect) {
                sprintf("Frame %d/%d | %.1f FPS | %d Dets", i, num_frames, sm_stats$fps, nrow(df_b))
              } else {
                sprintf("Frame %d/%d | %.1f FPS", i, num_frames, sm_stats$fps)
              }
              sparkline_pass <- NULL
            }
          } else {
            hud_txt <- ""
            sparkline_pass <- NULL
          }

          t_gui_start <- Sys.time()
          draw_yolo_detections_bgr_cpp(
            bm = bm,
            xmin = if (detect && nrow(df_b) > 0) df_b$xmin else numeric(0),
            ymin = if (detect && nrow(df_b) > 0) df_b$ymin else numeric(0),
            xmax = if (detect && nrow(df_b) > 0) df_b$xmax else numeric(0),
            ymax = if (detect && nrow(df_b) > 0) df_b$ymax else numeric(0),
            labels = if (detect && nrow(df_b) > 0) df_b$label else character(0),
            scores = if (detect && nrow(df_b) > 0) df_b$score else numeric(0),
            colors = palette_colors,
            lwd = as.integer(lwd),
            hud_text = hud_txt,
            font_scale = as.numeric(f_scale_obj),
            keypoints = kpt_mat_pass,
            kpt_threshold = as.numeric(kpt_threshold),
            kpt_radius = as.integer(kpt_radius),
            draw_skeleton = isTRUE(skeleton),
            draw_boxes = isTRUE(show_box),
            track_ids = track_ids_pass,
            flash = flash_pass,
            history_x = history_x_pass,
            history_y = history_y_pass,
            roi = roi_coords,
            roi_label = roi_label,
            count_line = line_coords,
            count_line_label = if (!is.null(count_line_label)) as.character(count_line_label) else "",
            counted = if ("counted" %in% names(df_b)) df_b$counted else NULL,
            hide_outside_roi = isTRUE(hide_outside_roi),
            hide_counted = isTRUE(hide_counted),
            flash_counted = isTRUE(flash_counted),
            show_text = isTRUE(show_text),
            show_conf = isTRUE(show_conf),
            show_class = isTRUE(show_class),
            show_id = isTRUE(show_id),
            badge = badge,
            pad_scale = pad_num,
            mask_labels = if (isTRUE(mask) && has_masks) raw_res$labels else NULL,
            mask_alpha = as.numeric(alpha),
            mask_offset_x = if (use_crop) as.integer(crop_coords[1]) else 0L,
            mask_offset_y = if (use_crop) as.integer(crop_coords[2]) else 0L,
            draw_masks = isTRUE(mask),
            mask_ids = if (has_masks && nrow(df_b) > 0) as.integer(df_b$orig_idx) else NULL,
            hud_pos = as.character(hud_pos),
            hud_layout = as.character(hud_layout),
            font_scale_hud = as.numeric(f_scale_hud),
            font_scale_roi = as.numeric(roi_sz),
            detailed_telemetry = is_detailed,
            sparkline = sparkline_pass,
            show_point = isTRUE(show_point),
            point_radius = as.integer(point_radius),
            point_col = if (!is.null(point_col)) as.character(point_col) else ""
          )
          if (save_video) {
            if (use_native_writer) {
              write_native_video_frame_cpp(v_writer, bm)
            } else if (use_native_stream) {
              f_rec_path <- file.path(run_dir, sprintf("rec_%06d.bmp", i))
              save_frame_bgr_cpp(bm, f_rec_path)
              recorded_frames[[i]] <<- f_rec_path
            } else {
              save_frame_bgr_cpp(bm, f_path, quality = 90L)
            }
          }
          if (use_native_win_b && !is.null(win_b)) {
            show_native_frame_cpp(win_b, bm, preserve_aspect = TRUE)
            status <- poll_native_window_cpp(win_b)
            if (status == 1L) break
            while (status == 2L) {
              Sys.sleep(0.05)
              status <- poll_native_window_cpp(win_b)
              if (status == 1L) break
            }
            if (status == 1L) break
          }
          t_gui_end <- Sys.time()
          last_gui_ms <<- max(1.0, as.numeric(t_gui_end - t_gui_start, units = "secs") * 1000)
        } else if (isTRUE(show) && !use_native_win_b) {
          f_mat <- if (use_native_stream) image_data(as_image(bm)) else image_data(image_import(f_path))
          plot(as_image(f_mat), ...)
          if (!is.null(roi_coords)) {
            graphics::rect(roi_coords[1], roi_coords[2], roi_coords[3], roi_coords[4],
                           col = grDevices::adjustcolor("#FA8072", alpha.f = 0.3),
                           border = "#FA8072", lwd = 1)
            graphics::text(roi_coords[1] + 5, roi_coords[2] - 5, roi_label, col = "#FA8072", cex = 2, font = 2, adj = c(0, 0))
          }
          if (!is.null(line_coords)) {
            graphics::lines(c(line_coords[1], line_coords[3]), c(line_coords[2], line_coords[4]), col = "#00E5FF", lwd = 3)
            if (!is.null(count_line_label) && nzchar(as.character(count_line_label))) {
              is_vert <- abs(line_coords[4] - line_coords[2]) >= abs(line_coords[3] - line_coords[1])
              if (is_vert) {
                graphics::text(line_coords[1] - 5, min(line_coords[2], line_coords[4]) + 15, count_line_label, col = "#00E5FF", cex = 0.7, font = 2, adj = c(1, 0))
              } else {
                graphics::text(min(line_coords[1], line_coords[3]) + 5, line_coords[2] - 5, count_line_label, col = "#00E5FF", cex = 0.7, font = 2, adj = c(0, 0))
              }
            }
          }
          if (detect && nrow(df_b) > 0) {
            # Filter bounding boxes to plot based on ROI and counted status
            keep_plot <- rep(TRUE, nrow(df_b))
            if (isTRUE(hide_outside_roi) && !is.null(roi_coords)) {
              cxs <- (df_b$xmin + df_b$xmax) / 2
              cys <- (df_b$ymin + df_b$ymax) / 2
              keep_plot <- keep_plot & (cxs >= roi_coords[1] & cxs <= roi_coords[3] & cys >= roi_coords[2] & cys <= roi_coords[4])
            }
            if (isTRUE(hide_counted) && "counted" %in% names(df_b)) {
              is_flashing <- if (!is.null(flash_pass) && length(flash_pass) == nrow(df_b)) (flash_pass > 0) else FALSE
              if (isTRUE(flash_counted)) {
                keep_plot <- keep_plot & (!df_b$counted | is_flashing)
              } else {
                keep_plot <- keep_plot & !df_b$counted
              }
            }

            df_plot <- df_b[keep_plot, , drop = FALSE]

            if (nrow(df_plot) > 0) {
              plot_colors <- palette_colors[keep_plot]

              conts_frame <- if (!is.null(raw_res$contours)) raw_res$contours else list()
              if (isTRUE(mask) && has_masks && length(conts_frame) > 0L) {
                for (kp_idx in which(keep_plot)) {
                  mid <- df_b$orig_idx[kp_idx]
                  if (!is.na(mid) && mid <= length(conts_frame) && !is.null(conts_frame[[mid]]) && nrow(conts_frame[[mid]]) >= 3) {
                    graphics::polygon(conts_frame[[mid]][, 1], conts_frame[[mid]][, 2],
                                      col = grDevices::adjustcolor(palette_colors[kp_idx], alpha.f = alpha),
                                      border = NA)
                  }
                }
              }

              if (isTRUE(bbox)) {
                .plot_yolo_bboxes(
                  boxes = df_plot,
                  palette_colors = plot_colors,
                  lwd = lwd,
                  cex = f_scale_obj,
                  pad = pad,
                  badge = badge,
                  show_text = isTRUE(show_text),
                  show_conf = isTRUE(show_conf),
                  show_class = isTRUE(show_class),
                  show_id = isTRUE(show_id)
                )
              }
              if (isTRUE(skeleton) && has_kpts) {
                kpts_for_plot <- lapply(which(keep_plot), function(k) {
                  vals <- raw_res$keypoints[k, ]
                  data.frame(
                    id = df_b$id[k],
                    keypoint = .coco_keypoints,
                    x = vals[seq(1, 51, by = 3)],
                    y = vals[seq(2, 51, by = 3)],
                    conf = round(vals[seq(3, 51, by = 3)], 4),
                    stringsAsFactors = FALSE
                  )
                })
                .plot_yolo_keypoints(
                  keypoints_list = kpts_for_plot,
                  palette_colors = plot_colors,
                  kpt_threshold = kpt_threshold,
                  lwd = lwd,
                  kpt_radius = kpt_radius
                )
              }
            }
          }
          if (isTRUE(show_fps) && f_scale_hud > 0.05) {
            elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
            cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
            hud_items <- if (isTRUE(track) || !is.null(count_line)) {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est),
                sprintf("IN ZONE: %d", in_zone_val),
                sprintf("COUNTED: %d", total_count_val))
            } else if (detect) {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est),
                sprintf("OBJECTS: %d", nrow(df_b)))
            } else {
              c(sprintf("FRAME: %d/%d", i, num_frames),
                sprintf("SPEED: %.1f FPS", cur_fps_est))
            }
            u <- graphics::par("usr")
            th <- abs(graphics::strheight("Ag", units = "user", cex = 0.65 * f_scale_hud))
            pad <- 6 * f_scale_hud
            is_vert <- (hud_layout == "vertical")

            if (is_vert) {
              all_lines <- c("TELEMETRY", hud_items)
              max_w <- max(abs(graphics::strwidth(all_lines, units = "user", cex = 0.65 * f_scale_hud, font = 2)))
              bx_w <- max_w + 2 * pad + 10 * f_scale_hud
              bx_h <- (length(all_lines) * th * 1.5) + 2 * pad
            } else {
              single_line <- paste(hud_items, collapse = "  |  ")
              bx_w <- abs(graphics::strwidth(single_line, units = "user", cex = 0.65 * f_scale_hud, font = 2)) + 2 * pad
              bx_h <- th * 1.8 + 2 * pad
            }

            margin_x <- abs(u[2] - u[1]) * 0.02
            margin_y <- abs(u[4] - u[3]) * 0.02
            top_y <- max(u[3], u[4]) - margin_y
            bot_y <- min(u[3], u[4]) + margin_y
            left_x <- min(u[1], u[2]) + margin_x
            right_x <- max(u[1], u[2]) - margin_x

            if (grepl("^[0-9.]+[ ,_][0-9.]+$|^[0-9.]+,[0-9.]+$", hud_pos)) {
              coords <- as.numeric(strsplit(hud_pos, "[, _]+")[[1]])
              cx <- if (coords[1] <= 1.0) left_x + coords[1] * (right_x - left_x) else left_x + coords[1]
              cy <- if (coords[2] <= 1.0) top_y - coords[2] * (top_y - bot_y) else top_y - coords[2]
              x1 <- cx; x2 <- cx + bx_w; y2 <- cy; y1 <- cy - bx_h
              if (x2 > right_x) { x2 <- right_x; x1 <- x2 - bx_w }
              if (x1 < left_x) { x1 <- left_x; x2 <- x1 + bx_w }
              if (y1 < bot_y) { y1 <- bot_y; y2 <- y1 + bx_h }
              if (y2 > top_y) { y2 <- top_y; y1 <- y2 - bx_h }
            } else if (hud_pos %in% c("top-right", "topright")) {
              x1 <- right_x - bx_w; x2 <- right_x; y2 <- top_y; y1 <- top_y - bx_h
            } else if (hud_pos %in% c("bottom-right", "bottomright")) {
              x1 <- right_x - bx_w; x2 <- right_x; y1 <- bot_y; y2 <- bot_y + bx_h
            } else if (hud_pos %in% c("bottom-left", "bottomleft")) {
              x1 <- left_x; x2 <- left_x + bx_w; y1 <- bot_y; y2 <- bot_y + bx_h
            } else if (hud_pos == "top") {
              mid_x <- (u[1] + u[2]) / 2; x1 <- mid_x - bx_w / 2; x2 <- mid_x + bx_w / 2; y2 <- top_y; y1 <- top_y - bx_h
            } else if (hud_pos == "bottom") {
              mid_x <- (u[1] + u[2]) / 2; x1 <- mid_x - bx_w / 2; x2 <- mid_x + bx_w / 2; y1 <- bot_y; y2 <- bot_y + bx_h
            } else {
              x1 <- left_x; x2 <- left_x + bx_w; y2 <- top_y; y1 <- top_y - bx_h
            }

            graphics::rect(xleft = x1, ybottom = y1, xright = x2, ytop = y2,
                           col = grDevices::adjustcolor("#0E1420", alpha.f = 0.85), border = "#233C55", lwd = 1.2)
            graphics::rect(xleft = x1, ybottom = y2 - 2, xright = x2, ytop = y2,
                           col = "#00E5FF", border = NA)
            if (is_vert) {
              graphics::points(x1 + pad + 2, y2 - pad - th / 2, pch = 16, col = "#00FF88", cex = 0.8 * f_scale_hud)
              graphics::text(x = x1 + pad + 10 * f_scale_hud, y = y2 - pad - th / 2, labels = "TELEMETRY", col = "#00E5FF", cex = 0.6 * f_scale_hud, font = 2, adj = c(0, 0.5))
              for (k in seq_along(hud_items)) {
                item_y <- y2 - pad - th * 1.5 - (k - 0.5) * th * 1.5
                parts <- strsplit(hud_items[k], ": ")[[1]]
                lbl <- parts[1]; val <- if (length(parts) > 1) parts[2] else ""
                val_col <- if (grepl("COUNT", lbl)) "#00E5FF" else if (grepl("ZONE", lbl)) "#FFBB00" else if (grepl("SPEED", lbl)) "#00FF88" else "#FFFFFF"
                graphics::text(x = x1 + pad, y = item_y, labels = lbl, col = "#94B0C4", cex = 0.6 * f_scale_hud, font = 2, adj = c(0, 0.5))
                if (nzchar(val)) {
                  graphics::text(x = x2 - pad, y = item_y, labels = val, col = val_col, cex = 0.6 * f_scale_hud, font = 2, adj = c(1, 0.5))
                }
              }
            } else {
              graphics::text(x = x1 + pad, y = (y1 + y2) / 2, labels = single_line, col = "#FFFFFF", cex = 0.65 * f_scale_hud, font = 2, adj = c(0, 0.5))
            }
          }
        }

        frames_processed <<- i

        if (!is.null(pb) && (i %% pb_step == 0L || i == num_frames)) {
          elapsed_now <- as.numeric(Sys.time() - t_start, units = "secs")
          cur_fps_est <- if (elapsed_now > 0) (i / elapsed_now) else src_fps
          cli::cli_progress_update(id = pb, set = i, extra = list(current_fps = sprintf("%.1f", cur_fps_est)))
        }
      }
    }

    process_loop()
    if (!is.null(pb)) try(cli::cli_progress_done(id = pb), silent = TRUE)
    t_loop_end <- Sys.time()
    loop_time <- max(0.001, as.numeric(t_loop_end - t_start, units = "secs"))
    actual_n <- frames_processed
    avg_fps <- if (loop_time > 0 && actual_n > 0) (actual_n / loop_time) else src_fps

    if (use_native_stream) {
      try(close_native_video_reader_cpp(v_reader), silent = TRUE)
    }

    if (use_native_writer) {
      try(close_native_video_writer_cpp(v_writer), silent = TRUE)
      if (isTRUE(verbose)) cli::cli_alert_success("Video encoding complete: {.path {out_video_file}}")
    } else if (save_video) {
      if (isTRUE(verbose)) {
        cli::cli_alert_info("Encoding annotated video to MP4 with FFmpeg...")
      }
      frames_to_encode <- if (use_native_stream) {
        valid_rec <- recorded_frames[!sapply(recorded_frames, is.null)]
        unlist(valid_rec)
      } else {
        frame_paths
      }
      av::av_encode_video(
        input = frames_to_encode,
        output = out_video_file,
        framerate = src_fps,
        verbose = FALSE
      )
      if (isTRUE(verbose)) cli::cli_alert_success("Video encoding complete: {.path {out_video_file}}")
    }

    valid_records <- frame_records[!sapply(frame_records, is.null)]
    valid_kpts <- frame_keypoints[!sapply(frame_keypoints, is.null)]
    all_dets <- if (length(valid_records) > 0) do.call(rbind, valid_records) else data.frame()
    all_kpts <- if (length(valid_kpts) > 0) do.call(rbind, valid_kpts) else NULL
    total_time <- as.numeric(Sys.time() - t_start, units = "secs")
  }

  if (!return_data) {
    if (isTRUE(verbose)) {
      cli::cli_alert_success("Session complete: {if (is_camera) frame_count else actual_n} frame(s) processed at {round(avg_fps, 1)} FPS average.")
      if (save_video && !is.null(output) && file.exists(output)) {
        cli::cli_alert_success("Video saved to: {.path {output}}")
      }
    }
    return(invisible(NULL))
  }

  # ============================================================================
  # COMPREHENSIVE DATA & ECOLOGICAL SUMMARY (RICHNESS, ABUNDANCE & DIVERSITY)
  # ============================================================================
  tot_frames <- if (is_camera) frame_count else actual_n
  if (is_camera && avg_fps > 0) src_fps <- avg_fps
  duration_sec <- if (is_camera) t_elapsed else (tot_frames / src_fps)
  v_source <- if (is_camera) "Live Webcam" else normalizePath(video, winslash = "/")

  # Filter transient noise detections if tracking is enabled (removes objects in <= min_frames)
  valid_track_ids <- NULL
  if (isTRUE(track) && !is.null(all_dets) && nrow(all_dets) > 0 && "id" %in% names(all_dets) &&
      !is.null(min_frames) && as.integer(min_frames) > 0L) {
    min_f_val <- as.integer(min_frames)
    id_tab <- table(all_dets$id[all_dets$id > 0L])
    valid_track_ids <- as.integer(names(id_tab)[id_tab > min_f_val])
    all_dets <- all_dets[all_dets$id %in% valid_track_ids, , drop = FALSE]
  }

  if (is.null(all_dets) || nrow(all_dets) == 0) {
    all_dets <- data.frame(
      frame = integer(0), timestamp = numeric(0), id = integer(0),
      xmin = numeric(0), ymin = numeric(0), xmax = numeric(0), ymax = numeric(0),
      label = character(0), score = numeric(0), class_id = integer(0),
      stringsAsFactors = FALSE
    )
    all_classes <- character(0)
    summary_df <- data.frame(
      video_source = basename(v_source),
      model = basename(model_file),
      engine = toupper(engine),
      total_frames = tot_frames,
      fps = round(avg_fps, 2),
      duration_sec = round(duration_sec, 2),
      total_detections = 0L,
      total_richness = 0L,
      unique_objects = if (isTRUE(track)) 0L else NA_integer_,
      avg_abundance_frame = 0.0,
      max_abundance_frame = 0L,
      peak_frame = 1L,
      peak_time_sec = 0.0,
      shannon_index = 0.0,
      simpson_index = 0.0,
      pielou_evenness = 0.0,
      stringsAsFactors = FALSE
    )
    by_class_df <- data.frame(
      class = character(0),
      total_detections = integer(0),
      relative_abundance = numeric(0),
      unique_objects = integer(0),
      n_frames_present = integer(0),
      occurrence_rate = numeric(0),
      mean_abundance_when_present = numeric(0),
      mean_abundance_overall = numeric(0),
      max_per_frame = integer(0),
      mean_confidence = numeric(0),
      first_seen_sec = numeric(0),
      last_seen_sec = numeric(0),
      stringsAsFactors = FALSE
    )
    temporal_df <- data.frame(
      frame = seq_len(max(1, tot_frames)),
      timestamp = round((seq_len(max(1, tot_frames)) - 1) / src_fps, 3),
      abundance = 0L,
      richness = 0L,
      shannon = 0.0,
      stringsAsFactors = FALSE
    )
    cumulative_df <- data.frame(
      frame = seq_len(max(1, tot_frames)),
      timestamp = round((seq_len(max(1, tot_frames)) - 1) / src_fps, 3),
      cumulative_detections = 0L,
      cumulative_richness = 0L,
      stringsAsFactors = FALSE
    )
    if (isTRUE(track)) cumulative_df$cumulative_unique_objects <- 0L
    diversity_metrics <- list(
      richness_S = 0L,
      total_abundance_N = 0L,
      shannon_H = 0.0,
      simpson_1_minus_D = 0.0,
      simpson_D = 0.0,
      pielou_J = 0.0,
      berger_parker_d = 0.0,
      mean_richness_per_frame = 0.0,
      mean_abundance_per_frame = 0.0
    )
  } else {
    all_classes <- sort(unique(all_dets$label))
    total_det <- nrow(all_dets)
    total_rich <- length(all_classes)

    # Cross-tabulate frame x class
    mat_counts <- table(
      factor(all_dets$frame, levels = seq_len(tot_frames)),
      factor(all_dets$label, levels = all_classes)
    )
    abundance_vec <- as.integer(rowSums(mat_counts))
    richness_vec <- as.integer(rowSums(mat_counts > 0))

    # Shannon index per frame
    props <- mat_counts / pmax(1, abundance_vec)
    log_props <- ifelse(props > 0, log(props), 0)
    shannon_vec <- -rowSums(props * log_props)

    timestamps_vec <- round((seq_len(tot_frames) - 1) / src_fps, 3)

    temporal_df <- data.frame(
      frame = seq_len(tot_frames),
      timestamp = timestamps_vec,
      abundance = abundance_vec,
      richness = richness_vec,
      shannon = round(shannon_vec, 3),
      stringsAsFactors = FALSE
    )
    temporal_df <- cbind(temporal_df, as.data.frame.matrix(mat_counts))

    # Cumulative richness and detections (Collector's Curve)
    first_frames <- tapply(all_dets$frame, all_dets$label, min)
    first_frame_counts <- tabulate(first_frames, nbins = tot_frames)
    cum_richness <- cumsum(first_frame_counts)
    cum_detections <- cumsum(abundance_vec)

    cumulative_df <- data.frame(
      frame = seq_len(tot_frames),
      timestamp = timestamps_vec,
      cumulative_detections = cum_detections,
      cumulative_richness = cum_richness,
      stringsAsFactors = FALSE
    )
    if (isTRUE(track) && "id" %in% names(all_dets) && any(all_dets$id > 0)) {
      first_seen_id <- tapply(all_dets$frame, all_dets$id, min)
      first_seen_counts <- tabulate(first_seen_id, nbins = tot_frames)
      cumulative_df$cumulative_unique_objects <- cumsum(first_seen_counts)
    }

    # Summary by class
    counts_by_cls <- if (isTRUE(track) || !is.null(count_line)) track_state$counts_by_class else NULL
    by_class_list <- lapply(all_classes, function(cls) {
      df_c <- all_dets[all_dets$label == cls, , drop = FALSE]
      n_c <- nrow(df_c)
      n_present <- sum(mat_counts[, cls] > 0)
      cnt_val <- if (!is.null(counts_by_cls) && cls %in% names(counts_by_cls)) as.integer(counts_by_cls[[cls]]) else NA_integer_
      data.frame(
        class = cls,
        total_detections = n_c,
        relative_abundance = round(n_c / total_det * 100, 2),
        unique_objects = if (isTRUE(track) && any(df_c$id > 0)) length(unique(df_c$id[df_c$id > 0])) else NA_integer_,
        n_frames_present = n_present,
        occurrence_rate = round(n_present / tot_frames * 100, 2),
        mean_abundance_when_present = round(n_c / max(1, n_present), 2),
        mean_abundance_overall = round(n_c / tot_frames, 2),
        max_per_frame = max(mat_counts[, cls]),
        mean_confidence = round(mean(df_c$score), 4),
        first_seen_sec = round(min(df_c$timestamp), 3),
        last_seen_sec = round(max(df_c$timestamp), 3),
        total_counted = cnt_val,
        stringsAsFactors = FALSE
      )
    })
    by_class_df <- do.call(rbind, by_class_list)
    by_class_df <- by_class_df[order(-by_class_df$total_detections), ]
    rownames(by_class_df) <- NULL

    if (is.null(counts_by_cls)) {
      by_class_df$total_counted <- NULL
    }

    # Overall Diversity Metrics
    p_all <- table(all_dets$label) / total_det
    shannon_total <- -sum(p_all * log(p_all))
    simpson_d <- sum(p_all^2)
    simpson_total <- 1.0 - simpson_d
    pielou_total <- if (total_rich > 1) shannon_total / log(total_rich) else if (total_rich == 1) 1.0 else 0.0

    peak_idx <- which.max(abundance_vec)

    summary_df <- data.frame(
      video_source = basename(v_source),
      model = basename(model_file),
      engine = toupper(engine),
      total_frames = tot_frames,
      fps = round(avg_fps, 2),
      duration_sec = round(duration_sec, 2),
      total_detections = total_det,
      total_richness = total_rich,
      unique_objects = if (isTRUE(track) && any(all_dets$id > 0)) length(unique(all_dets$id[all_dets$id > 0])) else NA_integer_,
      avg_abundance_frame = round(total_det / max(1, tot_frames), 2),
      max_abundance_frame = max(abundance_vec),
      peak_frame = peak_idx,
      peak_time_sec = timestamps_vec[peak_idx],
      shannon_index = round(shannon_total, 3),
      simpson_index = round(simpson_total, 3),
      pielou_evenness = round(pielou_total, 3),
      stringsAsFactors = FALSE
    )
    if (!is.null(total_count_val)) {
      summary_df$total_counted <- total_count_val
    }

    diversity_metrics <- list(
      richness_S = total_rich,
      total_abundance_N = total_det,
      shannon_H = round(shannon_total, 3),
      simpson_1_minus_D = round(simpson_total, 3),
      simpson_D = round(simpson_d, 3),
      pielou_J = round(pielou_total, 3),
      berger_parker_d = round(max(p_all), 3),
      mean_richness_per_frame = round(mean(richness_vec), 2),
      mean_abundance_per_frame = round(mean(abundance_vec), 2)
    )
  }

  # Aggregate per-unique-object summary table (tracking / shape / kinematic features)
  unique_objects_df <- data.frame()
  if ((isTRUE(track) || !is.null(count_line)) && !is.null(all_dets) && nrow(all_dets) > 0 && "id" %in% names(all_dets)) {
    valid_dets <- all_dets[!is.na(all_dets$id) & all_dets$id > 0L, , drop = FALSE]
    if (nrow(valid_dets) > 0) {
      split_u <- split(valid_dets, valid_dets$id)
      has_seg_feat <- "area" %in% names(valid_dets)
      uo_list <- lapply(split_u, function(sd) {
        sd <- sd[order(sd$frame), , drop = FALSE]
        first_f <- sd$frame[1]
        last_f <- sd$frame[nrow(sd)]
        n_f <- nrow(sd)
        dur <- round((last_f - first_f + 1) / src_fps, 2)
        score_m <- round(mean(sd$score, na.rm = TRUE), 3)
        col_val <- if ("color" %in% names(sd)) sd$color[1] else "#00CC66"
        lbl_tbl <- table(sd$label)
        lbl_val <- names(lbl_tbl)[which.max(lbl_tbl)]
        cx_vec <- 0.5 * (sd$xmin + sd$xmax)
        cy_vec <- 0.5 * (sd$ymin + sd$ymax)
        start_x <- round(cx_vec[1], 1)
        start_y <- round(cy_vec[1], 1)
        end_x <- round(cx_vec[n_f], 1)
        end_y <- round(cy_vec[n_f], 1)
        disp_px <- round(sqrt((cx_vec[n_f] - cx_vec[1])^2 + (cy_vec[n_f] - cy_vec[1])^2), 1)
        dist_px <- if (n_f > 1L) {
          round(sum(sqrt(diff(cx_vec)^2 + diff(cy_vec)^2)), 1)
        } else 0.0
        straight_val <- if (dist_px > 0) round(min(1.0, disp_px / dist_px), 3) else 0.0
        dt_span <- max(0.001, (last_f - first_f) / src_fps)
        speed_px_s <- if (n_f > 1L) round(dist_px / dt_span, 1) else 0.0

        row_df <- data.frame(
          id = sd$id[1],
          label = lbl_val,
          n_frames = n_f,
          first_frame = first_f,
          last_frame = last_f,
          duration_sec = dur,
          mean_score = score_m,
          start_x = start_x,
          start_y = start_y,
          end_x = end_x,
          end_y = end_y,
          distance_px = dist_px,
          displacement_px = disp_px,
          straightness = straight_val,
          mean_speed_px_s = speed_px_s,
          stringsAsFactors = FALSE
        )
        if (has_seg_feat) {
          row_df$area <- round(mean(sd$area, na.rm = TRUE), 1)
          row_df$area_sd <- if (n_f > 1) round(stats::sd(sd$area, na.rm = TRUE), 1) else 0.0
          row_df$perimeter <- round(mean(sd$perimeter, na.rm = TRUE), 1)
          row_df$radius_mean <- round(mean(sd$radius_mean, na.rm = TRUE), 2)
          row_df$length <- round(mean(sd$length, na.rm = TRUE), 1)
          row_df$width <- round(mean(sd$width, na.rm = TRUE), 1)
          row_df$circularity <- round(mean(sd$circularity, na.rm = TRUE), 4)
          row_df$eccentricity <- round(mean(sd$eccentricity, na.rm = TRUE), 4)
          row_df$asp_ratio <- round(mean(sd$asp_ratio, na.rm = TRUE), 3)
        } else {
          row_df$bbox_width <- round(mean(sd$xmax - sd$xmin, na.rm = TRUE), 1)
          row_df$bbox_height <- round(mean(sd$ymax - sd$ymin, na.rm = TRUE), 1)
          row_df$bbox_area <- round(mean((sd$xmax - sd$xmin) * (sd$ymax - sd$ymin), na.rm = TRUE), 1)
          row_df$asp_ratio <- round(mean((sd$xmax - sd$xmin) / pmax(0.001, (sd$ymax - sd$ymin)), na.rm = TRUE), 3)
        }
        row_df$color <- col_val
        row_df
      })
      unique_objects_df <- do.call(rbind, uo_list)
      rownames(unique_objects_df) <- NULL
    }
  }

  frame_w_out <- if (is_camera) orig_w else src_w
  frame_h_out <- if (is_camera) orig_h else src_h

  res <- list(
    summary = summary_df,
    by_class = by_class_df,
    temporal = temporal_df,
    cumulative = cumulative_df,
    diversity = diversity_metrics,
    detections = all_dets,
    unique_objects = unique_objects_df,
    frame_dims = c(width = as.integer(frame_w_out), height = as.integer(frame_h_out))
  )

  if (!is.null(heatmap_accum)) {
    res$heatmap <- heatmap_accum
  }
  if (isTRUE(save_crops) && !is.null(crops_dir)) {
    res$crops_dir <- normalizePath(crops_dir, winslash = "/", mustWork = FALSE)
  }

  if (isTRUE(track) || !is.null(count_line)) {
    if (length(track_state$crossing_events) > 0) {
      ev_df <- do.call(rbind, lapply(track_state$crossing_events, as.data.frame))
      if (!is.null(valid_track_ids)) {
        ev_df <- ev_df[ev_df$id %in% valid_track_ids, , drop = FALSE]
      }
      if ("direction" %in% names(ev_df)) {
        res$forward_count <- sum(ev_df$direction == "forward", na.rm = TRUE)
        res$backward_count <- sum(ev_df$direction == "backward", na.rm = TRUE)
        if (count_direction != "both") {
          ev_df <- ev_df[ev_df$direction == count_direction, , drop = FALSE]
        }
      }
      res$crossings <- ev_df
      if (!is.null(total_count_val)) {
        res$total_counted <- nrow(ev_df)
        res$summary$total_counted <- nrow(ev_df)
      }
    } else if (!is.null(total_count_val)) {
      res$total_counted <- total_count_val
    }
  }
  if (isTRUE(track)) res$min_frames <- min_frames
  if (!is.null(all_kpts) && nrow(all_kpts) > 0) res$keypoints <- all_kpts
  if (!is.null(roi)) {
    res$roi <- roi
    res$roi_coords <- roi_coords
    res$roi_orientation <- roi_orientation
    res$roi_mode <- roi_mode
  }
  if (!is.null(count_line)) {
    res$count_line <- count_line
    res$count_line_orientation <- count_line_orientation
    res$count_direction <- count_direction
  }
  if (!is.null(output)) res$output <- normalizePath(output, winslash = "/", mustWork = FALSE)

  class(res) <- c("pliman_video_detect", "list")

  if (isTRUE(verbose)) {
    cli::cli_alert_success("Inference complete: {tot_frames} frames processed at {round(avg_fps, 1)} FPS average.")
    if (!is.null(output) && file.exists(output)) {
      cli::cli_alert_success("Annotated video saved to: {.path {output}}")
    }
  }

  return(res)
}

#' @export
print.pliman_video_detect <- function(x, ...) {
  s <- x$summary
  cli::cli_h2("YOLO Video Detection & Ecological Summary")
  cli::cli_bullets(c(
    "*" = sprintf("Source: {.file %s} | Model: {.val %s} [{toupper(s$engine)}]", s$video_source, s$model),
    "*" = sprintf("Duration: %.1f s (%d frames at %.1f FPS)", s$duration_sec, s$total_frames, s$fps),
    "*" = sprintf("Total Detections: %s (Avg: %.1f objs/frame | Peak: %d objs at frame %d [%.1fs])",
                  format(s$total_detections, big.mark = ","), s$avg_abundance_frame, s$max_abundance_frame, s$peak_frame, s$peak_time_sec),
    "*" = sprintf("Richness (S): %d distinct classes observed", s$total_richness),
    "*" = sprintf("Diversity Indices: Shannon (H') = %.3f | Simpson (1-D) = %.3f | Pielou (J') = %.3f",
                  s$shannon_index, s$simpson_index, s$pielou_evenness)
  ))

  if (!is.na(s$unique_objects)) {
    cli::cli_bullets(c(
      "*" = sprintf("Persistent Tracked Objects: %d unique IDs", s$unique_objects)
    ))
  }
  if ("total_counted" %in% names(s)) {
    cli::cli_bullets(c(
      "*" = sprintf("Tripwire Counting (Line): %d objects crossed", s$total_counted)
    ))
  }

  bc <- x$by_class
  if (!is.null(bc) && nrow(bc) > 0) {
    cli::cli_h3("Class Breakdown (Abundance & Occurrence)")
    show_cols <- c("class", "total_detections", "relative_abundance", "occurrence_rate", "max_per_frame", "mean_confidence")
    if (!is.na(s$unique_objects) && "unique_objects" %in% names(bc)) {
      show_cols <- c("class", "total_detections", "relative_abundance", "unique_objects", "occurrence_rate", "max_per_frame", "mean_confidence")
    }
    if ("total_counted" %in% names(bc)) {
      show_cols <- c(show_cols, "total_counted")
    }
    print(bc[, show_cols, drop = FALSE], row.names = FALSE)
  }

  uo <- x$unique_objects
  if (!is.null(uo) && nrow(uo) > 0) {
    cli::cli_h3("Unique Tracked Objects Summary")
    if (nrow(uo) <= 10) {
      print(uo, row.names = FALSE)
    } else {
      print(head(uo, 8), row.names = FALSE)
      cli::cli_alert_info("... and {nrow(uo) - 8} more objects (access full table with $unique_objects)")
    }
  }

  if (!is.null(x$output)) {
    cli::cli_alert_info("Annotated Video: {.file {x$output}}")
  }

  cli::cli_text("{.muted List components: $summary, $by_class, $temporal, $cumulative, $diversity, $detections, $unique_objects}")
  invisible(x)
}

#' @export
summary.pliman_video_detect <- function(object, ...) {
  print(object, ...)
  invisible(list(summary = object$summary, by_class = object$by_class, diversity = object$diversity, unique_objects = object$unique_objects))
}

#' @export
as.data.frame.pliman_video_detect <- function(x, ...) {
  x$detections
}

#' @importFrom utils head tail
#' @export
head.pliman_video_detect <- function(x, n = 6L, ...) {
  head(x$detections, n, ...)
}

#' @export
tail.pliman_video_detect <- function(x, n = 6L, ...) {
  tail(x$detections, n, ...)
}

#' @export
dim.pliman_video_detect <- function(x) {
  dim(x$detections)
}

nrow.pliman_video_detect <- function(x) {
  nrow(x$detections)
}

ncol.pliman_video_detect <- function(x) {
  ncol(x$detections)
}

#' @export
`$.pliman_video_detect` <- function(x, name) {
  if (name %in% names(x)) {
    return(.subset2(x, name))
  }
  if ("detections" %in% names(x) && name %in% names(x$detections)) {
    return(.subset2(x$detections, name))
  }
  NULL
}

#' @export
`[[.pliman_video_detect` <- function(x, i, exact = TRUE) {
  if (is.character(i) && !(i %in% names(x)) && "detections" %in% names(x) && i %in% names(x$detections)) {
    return(x$detections[[i, exact = exact]])
  }
  NextMethod("[[")
}

#' @export
`[.pliman_video_detect` <- function(x, i, j, ...) {
  if (missing(j)) {
    if (is.character(i) && all(i %in% names(x))) {
      return(NextMethod("["))
    }
    return(x$detections[i, , drop = FALSE])
  }
  x$detections[i, j, ...]
}

#' @export
plot.pliman_video_detect <- function(x, type = c("temporal", "richness", "abundance", "class", "cumulative", "trajectories", "heatmap", "all"), ...) {
  type <- match.arg(type)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  temp <- x$temporal
  cum <- x$cumulative
  bcls <- x$by_class
  if (is.null(temp) || nrow(temp) == 0 || is.null(bcls) || nrow(bcls) == 0) {
    message("No detections to plot.")
    return(invisible(NULL))
  }

  classes <- bcls$class
  n_cls <- length(classes)
  colors <- .get_class_palette(classes, n_cls, rainbow = TRUE)
  names(colors) <- classes

  if (type == "temporal") {
    graphics::par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    # Panel 1: Abundance over time
    plot(temp$timestamp, temp$abundance, type = "n",
         xlab = "", ylab = "Abundance (N)",
         main = sprintf("Temporal Dynamics of Detections (Abundance & Richness) - %s", x$summary$model),
         bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$abundance, 0),
                      col = grDevices::adjustcolor("#00CC66", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$abundance, col = "#00994C", lwd = 2.5)

    # Panel 2: Richness over time
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "n",
         xlab = "Time (seconds)", ylab = "Richness (S)",
         bty = "l", las = 1, ylim = c(0, max_rich + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$richness, 0),
                      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$richness, col = "#0072B2", lwd = 2.5, type = "s")

  } else if (type == "abundance") {
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    max_c <- if (length(classes) > 0) max(temp[, classes, drop = FALSE], na.rm = TRUE) else 1
    plot(NULL, xlim = range(temp$timestamp), ylim = c(0, max_c + 1),
         xlab = "Time (seconds)", ylab = "Abundance per Class",
         main = "Class Abundance Over Time",
         bty = "l", las = 1, ...)
    graphics::grid(col = "#E0E0E0", lty = 2)
    for (cls in classes) {
      graphics::lines(temp$timestamp, temp[[cls]], col = colors[cls], lwd = 2)
    }
    graphics::legend("topright", legend = classes, col = colors, lwd = 2, bty = "n", cex = 0.85)

  } else if (type == "richness") {
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "n",
         xlab = "Time (seconds)", ylab = "Class Richness (S)",
         main = "Class Richness Over Time",
         bty = "l", las = 1, ylim = c(0, max_rich + 1), ...)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::polygon(c(temp$timestamp[1], temp$timestamp, temp$timestamp[nrow(temp)]),
                      c(0, temp$richness, 0),
                      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.25), border = NA)
    graphics::lines(temp$timestamp, temp$richness, col = "#0072B2", lwd = 2.5, type = "s")

  } else if (type == "cumulative") {
    graphics::par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    max_cum_r <- max(cum$cumulative_richness, na.rm = TRUE)
    plot(cum$timestamp, cum$cumulative_richness, type = "n",
         xlab = "", ylab = "Cumulative Richness (S)",
         main = "Species Accumulation Curve (Collector's Curve)",
         bty = "l", las = 1, ylim = c(0, max_cum_r + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::lines(cum$timestamp, cum$cumulative_richness, col = "#D55E00", lwd = 2.5, type = "s")

    plot(cum$timestamp, cum$cumulative_detections, type = "n",
         xlab = "Time (seconds)", ylab = "Cumulative Detections",
         main = "Cumulative Detections Over Time",
         bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)
    graphics::lines(cum$timestamp, cum$cumulative_detections, col = "#0072B2", lwd = 2.5)

  } else if (type == "class") {
    graphics::par(mar = c(4.5, 7, 2.5, 2), mgp = c(2.5, 0.7, 0))
    ord_bcls <- bcls[order(bcls$total_detections), ]
    bp <- graphics::barplot(ord_bcls$total_detections, horiz = TRUE, names.arg = ord_bcls$class,
                            col = colors[ord_bcls$class], border = NA, las = 1,
                            xlab = "Total Detections", main = "Abundance by Class",
                            xlim = c(0, max(ord_bcls$total_detections) * 1.25), ...)
    graphics::grid(nx = NULL, ny = NA, col = "#E0E0E0", lty = 2)
    lbl_txt <- sprintf("%d (%.1f%%)", ord_bcls$total_detections, ord_bcls$relative_abundance)
    graphics::text(ord_bcls$total_detections, bp, labels = paste0("  ", lbl_txt), pos = 4, cex = 0.85, font = 2)

  } else if (type == "trajectories") {
    dets <- x$detections
    if (is.null(dets) || nrow(dets) == 0 || !("id" %in% names(dets)) || all(dets$id <= 0L)) {
      message("No tracked object trajectories to plot. Run with `track = TRUE`.")
      return(invisible(NULL))
    }
    f_dims <- if (!is.null(x$frame_dims)) x$frame_dims else c(width = max(dets$xmax, na.rm = TRUE), height = max(dets$ymax, na.rm = TRUE))
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    plot(NULL, xlim = c(0, f_dims[1]), ylim = c(f_dims[2], 0),
         xlab = "X Coordinate (px)", ylab = "Y Coordinate (px)",
         main = "2D Spatial Trajectories of Tracked Objects",
         asp = 1, bty = "o", las = 1, ...)
    graphics::grid(col = "#E0E0E0", lty = 2)
    valid_d <- dets[dets$id > 0L, , drop = FALSE]
    sp_ids <- split(valid_d, valid_d$id)
    for (sd in sp_ids) {
      sd <- sd[order(sd$frame), , drop = FALSE]
      cx <- 0.5 * (sd$xmin + sd$xmax)
      cy <- 0.5 * (sd$ymin + sd$ymax)
      col_tr <- if ("color" %in% names(sd) && !is.na(sd$color[1])) sd$color[1] else colors[sd$label[1]]
      if (length(cx) > 1L) {
        graphics::lines(cx, cy, col = col_tr, lwd = 2)
        graphics::arrows(cx[length(cx) - 1L], cy[length(cy) - 1L], cx[length(cx)], cy[length(cy)],
                         length = 0.08, col = col_tr, lwd = 2)
      }
      graphics::points(cx[1], cy[1], pch = 16, col = col_tr, cex = 0.9)
    }

  } else if (type == "heatmap") {
    dets <- x$detections
    f_dims <- if (!is.null(x$frame_dims)) x$frame_dims else c(width = max(dets$xmax, na.rm = TRUE), height = max(dets$ymax, na.rm = TRUE))
    hm_mat <- x$heatmap
    if (is.null(hm_mat)) {
      gx <- 96L
      gy <- max(24L, as.integer(round(96L * f_dims[2] / max(1, f_dims[1]))))
      cxs <- 0.5 * (dets$xmin + dets$xmax)
      cys <- 0.5 * (dets$ymin + dets$ymax)
      ix <- pmin(gx, pmax(1L, as.integer(ceiling(cxs / f_dims[1] * gx))))
      iy <- pmin(gy, pmax(1L, as.integer(ceiling(cys / f_dims[2] * gy))))
      hm_mat <- matrix(0, nrow = gx, ncol = gy)
      for (k in seq_along(ix)) {
        hm_mat[ix[k], iy[k]] <- hm_mat[ix[k], iy[k]] + 1
      }
    }
    graphics::par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.2, 0.7, 0))
    pal_hm <- grDevices::hcl.colors(64, "Inferno")
    xs <- seq(0, f_dims[1], length.out = nrow(hm_mat))
    ys <- seq(0, f_dims[2], length.out = ncol(hm_mat))
    graphics::image(xs, ys, hm_mat[, rev(seq_len(ncol(hm_mat))), drop = FALSE],
                    col = pal_hm, xlab = "X Coordinate (px)", ylab = "Y Coordinate (px (inverted))",
                    main = "2D Spatial Presence Density Heatmap", asp = f_dims[2] / f_dims[1], las = 1, ...)

  } else if (type == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(3.5, 4, 2.5, 1), mgp = c(2, 0.6, 0))
    plot(temp$timestamp, temp$abundance, type = "l", col = "#00994C", lwd = 2,
         xlab = "Time (s)", ylab = "Abundance (N)", main = "Total Abundance", bty = "l", las = 1)
    graphics::grid(col = "#E0E0E0", lty = 2)
    max_rich <- max(temp$richness, na.rm = TRUE)
    plot(temp$timestamp, temp$richness, type = "s", col = "#0072B2", lwd = 2,
         xlab = "Time (s)", ylab = "Richness (S)", main = "Class Richness", bty = "l", las = 1,
         ylim = c(0, max_rich + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
    ord_bcls <- bcls[order(bcls$total_detections), ]
    bp <- graphics::barplot(ord_bcls$total_detections, horiz = TRUE, names.arg = ord_bcls$class,
                            col = colors[ord_bcls$class], border = NA, las = 1,
                            xlab = "Detections", main = "Abundance by Class", cex.names = 0.75)
    graphics::grid(nx = NULL, ny = NA, col = "#E0E0E0", lty = 2)
    max_cum_r <- max(cum$cumulative_richness, na.rm = TRUE)
    plot(cum$timestamp, cum$cumulative_richness, type = "s", col = "#D55E00", lwd = 2,
         xlab = "Time (s)", ylab = "Cumulative Richness", main = "Collector's Curve", bty = "l", las = 1,
         ylim = c(0, max_cum_r + 1))
    graphics::grid(col = "#E0E0E0", lty = 2)
  }
  invisible(x)
}

#' @rdname video_detect_dl
#' @export
video_detect <- video_detect_dl
