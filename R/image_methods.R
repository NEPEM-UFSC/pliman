#' Plot method for `image` objects
#'
#' Render an `image` object using native base R graphics, with support for title,
#' custom bounding limits (`xlim`, `ylim`), and axes.
#'
#' @param x An `image` object.
#' @param y Unused.
#' @param colormode Character specifying colormode override (`"Color"` or `"Grayscale"`).
#' @param colormode Unused; reserved for compatibility.
#' @param main Title for the plot.
#' @param axes Logical; should axes be drawn? Default `FALSE`.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param xlim X-axis limits.
#' @param ylim Y-axis limits.
#' @param asp Aspect ratio. Default `1`.
#' @param interpolate Logical; should raster image be interpolated? Default `TRUE`.
#' @param max_pixels Maximum number of pixels to render on screen plot. Default
#'   `2e6` (2 MP) for fast rendering. Set to `NULL` for full resolution.
#' @param add Logical; add to existing plot window? Default `FALSE`.
#' @param ... Additional graphical arguments passed to `graphics::rasterImage()`.
#'
#' @return Invisibly returns the plotted `image` object.
#' @export
#'
#' @examples
#' \dontrun{
#' img <- image(array(as.raw(sample(0:255, 30000)), dim = c(100, 100, 3)))
#' plot(img, main = "Sample Image")
#' plot(img, xlim = c(20, 80), ylim = c(20, 80), main = "Zoomed Region")
#' }
plot.image <- function(x,
                       y = NULL,
                       colormode = NULL,
                       main = NULL,
                       axes = FALSE,
                       xlab = "",
                       ylab = "",
                       xlim = NULL,
                       ylim = NULL,
                       asp = 1,
                       interpolate = TRUE,
                       max_pixels = 2e6,
                       add = FALSE,
                       ...) {
  if (!is_image(x)) {
    cli::cli_abort("Object passed to {.fn plot.image} must be of class {.code image}.")
  }

  if (storage.mode(x) == "integer") {
    x <- image_color_labels(x)
  }

  dims <- dim(x)
  w <- dims[1]
  h <- dims[2]

  if (is.null(xlim)) xlim <- c(1, w)
  if (is.null(ylim)) ylim <- c(1, h)

  # Validate and clamp xlim & ylim
  xlim <- c(max(1, min(xlim)), min(w, max(xlim)))
  ylim <- c(max(1, min(ylim)), min(h, max(ylim)))

  # Crop if xlim or ylim are subsets of the image
  sub_img <- if (xlim[1] > 1 || xlim[2] < w || ylim[1] > 1 || ylim[2] < h) {
    if (length(dims) == 3) {
      x[xlim[1]:xlim[2], ylim[1]:ylim[2], , drop = FALSE]
    } else {
      x[xlim[1]:xlim[2], ylim[1]:ylim[2], drop = FALSE]
    }
  } else {
    x
  }

  sub_dims <- dim(sub_img)
  dw <- sub_dims[1]
  dh <- sub_dims[2]
  tot_pix <- dw * dh

  disp_img <- if (!is.null(max_pixels) && tot_pix > max_pixels) {
    step <- ceiling(sqrt(tot_pix / max_pixels))
    if (length(sub_dims) == 3) {
      sub_img[seq(1, dw, by = step), seq(1, dh, by = step), , drop = FALSE]
    } else {
      sub_img[seq(1, dw, by = step), seq(1, dh, by = step), drop = FALSE]
    }
  } else {
    sub_img
  }

  u_disp <- unclass(disp_img)
  rst <- cpp_as_native_raster(u_disp, dim(disp_img))

  if (!isTRUE(add)) {
    is_mf <- (prod(graphics::par("mfrow")) > 1L) || (prod(graphics::par("mfcol")) > 1L)
    if (!is_mf) {
      if (!isTRUE(axes) && is.null(main)) {
        op <- graphics::par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
        on.exit(graphics::par(op), add = TRUE)
      } else {
        top_mar <- if (is.null(main)) 2 else 4
        op <- graphics::par(mar = c(3, 3, top_mar, 1))
        on.exit(graphics::par(op), add = TRUE)
      }
    }

    graphics::plot.new()
    graphics::plot.window(
      xlim = c(xlim[1] - 0.5, xlim[2] + 0.5),
      ylim = c(ylim[2] + 0.5, ylim[1] - 0.5),
      asp = asp,
      xaxs = "i",
      yaxs = "i"
    )
  }

  graphics::rasterImage(
    rst,
    xleft = xlim[1] - 0.5,
    ybottom = ylim[2] + 0.5,
    xright = xlim[2] + 0.5,
    ytop = ylim[1] - 0.5,
    interpolate = interpolate,
    ...
  )

  if (isTRUE(axes)) {
    graphics::axis(1)
    graphics::axis(2)
    graphics::box(bty = "l")
  }

  if (!is.null(main)) {
    graphics::title(main = main)
  }

  invisible(x)
}

#' Generic function for interactive viewing of images and spatial objects
#'
#' @param x An image or spatial object to view.
#' @param ... Additional arguments passed to specific view methods.
#'
#' @export
view <- function(x, ...) {
  UseMethod("view")
}

#' @export
view.default <- function(x, ...) {
  if (exists("image_view", mode = "function")) {
    image_view(x, ...)
  } else {
    cli::cli_abort("No viewer method available for class {.val {class(x)[1]}}.")
  }
}

#' Interactive viewer method for `image` objects with Zoom In and Zoom Out
#'
#' Provides an interactive viewer for `image` objects in R. When using base viewer mode,
#' users can interactively zoom in (selecting a bounding rectangle), zoom out, reset,
#' or finish viewing.
#'
#' @param x An `image` object.
#' @param viewer Character specifying the viewer back-end: `"base"` (interactive console plot)
#'   or `"mapview"` (interactive Leaflet map). Default is retrieved via `get_pliman_viewer()`.
#' @param ... Additional arguments passed to [image_view()] or [plot.image()].
#'
#' @return Invisibly returns the current zoomed sub-image or full image object.
#' @export
#'
#' @examples
#' \dontrun{
#' img <- image_import(image_pliman("sev_leaf.jpg"))
#' view(img)
#' }
view.image <- function(x, viewer = get_pliman_viewer(), ...) {
  if (missing(viewer) && requireNamespace("mapview", quietly = TRUE)) {
    viewer <- "mapview"
  } else {
    viewer_opt <- c("base", "mapview")
    viewer <- viewer_opt[pmatch(viewer[[1]], viewer_opt)]
    if (is.na(viewer)) viewer <- "base"
  }

  if (viewer == "mapview" && requireNamespace("mapview", quietly = TRUE)) {
    return(image_view(x, ...))
  }

  # Interactive Base Viewer with Zoom In & Zoom Out
  dims <- dim(x)
  w <- dims[1]
  h <- dims[2]

  curr_xlim <- c(1, w)
  curr_ylim <- c(1, h)

  if (!interactive()) {
    plot.image(x, xlim = curr_xlim, ylim = curr_ylim, axes = TRUE, main = "Image View", ...)
    return(invisible(x))
  }

  cli::cli_inform(c(
    "i" = "Entering interactive Image Viewer.",
    "*" = "Controls: [1] Zoom In (Click 2 points), [2] Zoom Out, [3] Reset Zoom, [0/Q] Exit."
  ))

  repeat {
    plot.image(
      x,
      xlim = curr_xlim,
      ylim = curr_ylim,
      axes = TRUE,
      main = sprintf("View [%d:%d, %d:%d] | W:%d H:%d",
                     round(curr_xlim[1]), round(curr_xlim[2]),
                     round(curr_ylim[1]), round(curr_ylim[2]), w, h),
      ...
    )

    cat("\nInteractive Viewer Actions:\n")
    cat("  1: Zoom In (Click 2 corners on image)\n")
    cat("  2: Zoom Out (2x)\n")
    cat("  3: Reset Zoom (Full Image)\n")
    cat("  0: Done / Exit\n")
    ans <- readline(prompt = "Select option [0-3]: ")
    ans <- trimws(tolower(ans))

    if (ans %in% c("0", "q", "exit", "done", "")) {
      cli::cli_inform("Exited viewer.")
      break
    } else if (ans == "1") {
      cat("Click TWO points on the graphics window to define Zoom area...\n")
      pts <- graphics::locator(n = 2, type = "p", col = "red", pch = 3)
      if (!is.null(pts) && length(pts$x) == 2) {
        new_x <- sort(pts$x)
        new_y <- sort(pts$y)
        curr_xlim <- c(max(1, round(new_x[1])), min(w, round(new_x[2])))
        curr_ylim <- c(max(1, round(new_y[1])), min(h, round(new_y[2])))
        if (curr_xlim[2] <= curr_xlim[1]) curr_xlim[2] <- min(w, curr_xlim[1] + 10)
        if (curr_ylim[2] <= curr_ylim[1]) curr_ylim[2] <- min(h, curr_ylim[1] + 10)
      }
    } else if (ans == "2") {
      center_x <- mean(curr_xlim)
      center_y <- mean(curr_ylim)
      span_x <- (curr_xlim[2] - curr_xlim[1]) * 2
      span_y <- (curr_ylim[2] - curr_ylim[1]) * 2
      curr_xlim <- c(max(1, round(center_x - span_x / 2)), min(w, round(center_x + span_x / 2)))
      curr_ylim <- c(max(1, round(center_y - span_y / 2)), min(h, round(center_y + span_y / 2)))
    } else if (ans == "3") {
      curr_xlim <- c(1, w)
      curr_ylim <- c(1, h)
    }
  }

  res_sub <- x[curr_xlim[1]:curr_xlim[2], curr_ylim[1]:curr_ylim[2], , drop = FALSE]
  invisible(res_sub)
}
