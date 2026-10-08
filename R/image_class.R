#' The `image` class constructor and utilities
#'
#' @description
#' `image` creates or coerces objects into a memory-optimized `image` S3 object
#' in `pliman`. Storing pixel data efficiently (using 8-bit `raw` vectors by default)
#' provides up to 8x memory reduction compared to standard float64 double arrays.
#'
#' @param data An array, matrix, raw vector, or existing image object.
#' @param colormode Character specifying color mode: `"Color"` (default for 3D arrays)
#'   or `"Grayscale"` (default for 2D matrices).
#' @param colspace Character specifying color space (default `"RGB"`).
#' @param gamma Numeric gamma value (default `1.0`).
#' @param storage Character specifying storage type: `"raw"` (uint8 0..255, memory optimized)
#'   or `"double"` (0..1 floating point).
#' @param dim Optional dimensions vector to reshape data.
#' @param x An object to check if it inherits from class `image`.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class `c("image", "array")`.
#' @export
#'
#' @examples
#' \dontrun{
#' img <- image(array(as.raw(sample(0:255, 300, replace = TRUE)), dim = c(10, 10, 3)))
#' print(img)
#' summary(img)
#' is_image(img)
#' }
image <- function(data = raw(),
                  colormode = NULL,
                  colspace = "RGB",
                  gamma = 1.0,
                  storage = c("raw", "double", "integer", "logical"),
                  ...) {
  storage <- match.arg(storage)
  as_image(data, colormode = colormode, colspace = colspace, gamma = gamma, storage = storage, ...)
}

#' @rdname image
#' @export
as_image <- function(data,
                     colormode = NULL,
                     colspace = "RGB",
                     gamma = 1.0,
                     storage = c("raw", "double", "integer", "logical"),
                     dim = NULL,
                     ...) {
  dots <- list(...)
  if (is.null(dim) && !is.null(dots$dim)) dim <- dots$dim

  if (is.null(dim(data)) && !is.null(dim)) {
    dim(data) <- dim
  }

  if (!is.null(colormode)) {
    cm_check <- tolower(colormode)
    if (cm_check == "color" && (length(dim(data)) == 2 || (length(dim(data)) == 3 && dim(data)[3] == 1))) {
      d_cur <- dim(data)
      mat_src <- if (length(d_cur) == 3) data[, , 1] else data
      data <- array(rep(mat_src, 3), dim = c(d_cur[1], d_cur[2], 3))
    }
  }

  if (inherits(data, "image")) {
    target_storage <- match.arg(storage)
    current_storage <- if (is.raw(data)) "raw" else if (is.double(data)) "double" else if (is.integer(data)) "integer" else if (is.logical(data)) "logical" else "unknown"
    if (target_storage == current_storage) {
      if (!is.null(colormode)) attr(data, "colormode") <- colormode
      if (!is.null(colspace)) attr(data, "colspace") <- colspace
      if (!is.null(gamma)) attr(data, "gamma") <- gamma
      return(data)
    }
  }
  if (isS4(data) && .hasSlot(data, ".Data")) {
    cm <- if (is.null(colormode)) {
      if (.hasSlot(data, "colormode")) {
        cm_val <- slot(data, "colormode")
        if (is.numeric(cm_val)) ifelse(cm_val == 2, "Color", "Grayscale") else as.character(cm_val)
      } else {
        "Color"
      }
    } else colormode
    data <- slot(data, ".Data")
    colormode <- cm
  }

  dims <- dim(data)
  if (is.null(dims)) {
    cli::cli_abort("Data provided to {.fn as_image} must have dimensions (array or matrix).")
  }

  if (length(dims) < 2 || length(dims) > 3) {
    cli::cli_abort("Image dimensions must be 2D (grayscale) or 3D (multi-channel).")
  }

  if (is.null(colormode)) {
    colormode <- if (length(dims) == 3 && dims[3] >= 3) "Color" else "Grayscale"
  } else {
    colormode <- match.arg(colormode, c("Color", "Grayscale", "color", "grayscale"))
    if (tolower(colormode) == "color") colormode <- "Color"
    if (tolower(colormode) == "grayscale") colormode <- "Grayscale"
  }

  storage <- match.arg(storage)

  if (storage == "raw" && !is.raw(data)) {
    if (is.numeric(data)) {
      if (anyNA(data)) {
        storage <- "double"
      } else {
        mx <- max(data, na.rm = TRUE)
        if (mx <= 1.0 && mx >= 0) {
          data <- as.raw(round(data * 255))
        } else {
          data <- as.raw(pmin(pmax(round(data), 0), 255))
        }
        dim(data) <- dims
      }
    } else if (is.logical(data)) {
      data <- as.raw(data * 255)
      dim(data) <- dims
    }
  } else if (storage == "double" && is.raw(data)) {
    data <- as.numeric(data) / 255
    dim(data) <- dims
  } else if (storage == "integer" && !is.integer(data)) {
    data <- as.integer(data)
    dim(data) <- dims
  } else if (storage == "logical" && !is.logical(data)) {
    data <- data != 0
    dim(data) <- dims
  }

  class(data) <- c("image", "array")
  attr(data, "colormode") <- colormode
  attr(data, "colspace") <- colspace
  attr(data, "gamma") <- gamma
  return(data)
}

#' @rdname image
#' @export
is_image <- function(x) {
  inherits(x, "image")
}

#' Print method for `image` objects
#'
#' @param x An `image` object.
#' @param ... Unused.
#' @export
print.image <- function(x, ...) {
  dims <- dim(x)
  w <- dims[1]
  h <- dims[2]
  ch <- if (length(dims) == 3) dims[3] else 1
  cm <- attr(x, "colormode") %||% "Color"
  st <- typeof(x)
  sz <- format(utils::object.size(x), units = "auto")

  cli::cli_verbatim(cli::col_green(cli::style_bold("Image Object (pliman)")))
  cli::cli_verbatim(paste0(cli::col_cyan("  Dimensions: "), w, " x ", h, " (width x height)"))
  cli::cli_verbatim(paste0(cli::col_cyan("  Channels  : "), ch))
  cli::cli_verbatim(paste0(cli::col_cyan("  Colormode : "), cm))
  cli::cli_verbatim(paste0(cli::col_cyan("  Storage   : "), st, " (", sz, ")"))

  if (is.raw(x)) {
    rng <- range(as.integer(x), na.rm = TRUE)
    cli::cli_verbatim(paste0(cli::col_cyan("  Range     : "), "[", rng[1], ", ", rng[2], "] (uint8)"))
  } else if (is.logical(x)) {
    cli::cli_verbatim(paste0(cli::col_cyan("  Range     : "), "[FALSE, TRUE] (binary logical)"))
  } else if (is.numeric(x)) {
    rng <- range(x, na.rm = TRUE)
    cli::cli_verbatim(paste0(cli::col_cyan("  Range     : "), "[", round(rng[1], 4), ", ", round(rng[2], 4), "] (double)"))
  }
  invisible(x)
}

#' Subsetting method for `image` objects
#'
#' @param x An `image` object.
#' @param i Row/X index or missing.
#' @param j Column/Y index or missing.
#' @param k Channel index or missing.
#' @param ... Additional arguments.
#' @param drop Logical; whether to drop dimensions.
#' @export
`[.image` <- function(x, i, j, k, ..., drop = FALSE) {
  cm <- attr(x, "colormode")
  cs <- attr(x, "colspace")
  gm <- attr(x, "gamma")

  res <- NextMethod("[", drop = drop)

  if (!drop && !is.null(dim(res)) && length(dim(res)) >= 2) {
    class(res) <- c("image", "array")
    attr(res, "colormode") <- cm
    attr(res, "colspace") <- cs
    attr(res, "gamma") <- gm
  } else if (drop && !is.null(dim(res)) && length(dim(res)) >= 2) {
    class(res) <- c("image", "array")
    attr(res, "colormode") <- if (length(dim(res)) == 2) "Grayscale" else cm
    attr(res, "colspace") <- cs
    attr(res, "gamma") <- gm
  }
  res
}

#' @export
`[<-.image` <- function(x, ..., value) {
  cm <- attr(x, "colormode")
  cs <- attr(x, "colspace")
  gm <- attr(x, "gamma")

  u     <- unclass(x)
  val_u <- unclass(value)

  if (is.raw(u)) {
    if (is.numeric(val_u)) {
      mx <- suppressWarnings(max(val_u, na.rm = TRUE))
      if (is.finite(mx) && mx <= 1.0 && mx >= 0 && any(val_u > 0 & val_u < 1)) {
        val_u <- as.raw(round(val_u * 255))
      } else {
        val_u <- as.raw(pmin(pmax(round(val_u), 0), 255))
      }
    } else if (is.logical(val_u)) {
      val_u <- as.raw(as.integer(val_u))
    }
  }

  cl         <- match.call(expand.dots = TRUE)
  cl[[1]]    <- quote(`[<-`)
  cl[[2]]    <- quote(u)                     # first positional arg = object

  # Replace the value argument (last named element) with val_u symbol
  val_idx <- which(names(cl) == "value")
  if (length(val_idx) > 0L) {
    cl[[val_idx]] <- quote(val_u)
  }

  ee       <- new.env(parent = parent.frame())
  ee$u     <- u
  ee$val_u <- val_u

  u <- eval(cl, envir = ee)

  class(u)              <- c("image", "array")
  attr(u, "colormode") <- cm
  attr(u, "colspace")  <- cs
  attr(u, "gamma")     <- gm
  u
}

#' @export
as.array.image <- function(x, ...) {
  class(x) <- setdiff(class(x), "image")
  x
}

#' @export
as.matrix.image <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 3 && d[3] == 1) {
    res <- x[, , 1, drop = TRUE]
    return(res)
  }
  if (length(d) == 2) {
    class(x) <- setdiff(class(x), "image")
    return(x)
  }
  cli::cli_abort("Cannot coerce multi-channel 3D image to 2D matrix directly without indexing or channel selection.")
}

#' Convert `image` to standard R raster matrix
#'
#' @param x An `image` object.
#' @param ... Unused.
#' @export
as.raster.image <- function(x, ...) {
  d <- dim(x)
  u <- unclass(x)

  arr <- if (is.raw(u)) {
    as.numeric(u) / 255
  } else if (is.logical(u)) {
    as.numeric(u)
  } else {
    mx <- max(u, na.rm = TRUE)
    if (mx > 1.0) as.numeric(u) / mx else as.numeric(u)
  }

  dim(arr) <- d

  if (length(d) == 3 && d[3] >= 3) {
    arr_t <- aperm(arr[, , 1:min(d[3], 4), drop = FALSE], c(2, 1, 3))
  } else {
    mat2d <- if (length(d) == 3) arr[, , 1] else arr
    arr_t <- t(mat2d)
  }

  grDevices::as.raster(arr_t)
}

#' Group generic Summary method for `image` objects (min, max, range, sum, prod, any, all)
#'
#' @param ... `image` objects or numbers.
#' @param na.rm Logical; should NAs be removed?
#' @export
Summary.image <- function(..., na.rm = FALSE) {
  args <- list(...)
  vals <- unlist(lapply(args, function(arg) {
    u <- unclass(arg)
    if (is.raw(u)) {
      as.integer(u)
    } else if (is.logical(u)) {
      u
    } else if (.Generic %in% c("any", "all")) {
      as.logical(u)
    } else {
      as.numeric(u)
    }
  }), use.names = FALSE)

  op_func <- get(.Generic, envir = baseenv())
  op_func(vals, na.rm = na.rm)
}

#' Mean method for `image` objects
#'
#' @param x An `image` object.
#' @param ... Additional arguments.
#' @param na.rm Logical; should NAs be removed?
#' @export
mean.image <- function(x, ..., na.rm = FALSE) {
  u <- unclass(x)
  vals <- if (is.raw(u)) as.integer(u) else as.numeric(u)
  mean(vals, na.rm = na.rm)
}

#' Coerce `image` pixel values to numeric vector
#'
#' @param x An `image` object.
#' @param raw_to_scale Logical; if `TRUE` and image is `raw`, scales 0..255 to 0..1. Default `TRUE`.
#' @param ... Unused.
#' @export
as.numeric.image <- function(x, raw_to_scale = TRUE, ...) {
  u <- unclass(x)
  if (is.raw(u)) {
    if (isTRUE(raw_to_scale)) as.numeric(u) / 255 else as.numeric(u)
  } else {
    as.numeric(u)
  }
}

#' @export
as.double.image <- function(x, ...) {
  as.numeric.image(x, ...)
}

#' @export
as.integer.image <- function(x, ...) {
  u <- unclass(x)
  if (is.raw(u)) as.integer(u) else as.integer(u * 255)
}

#' Extract raw or numeric image data array/matrix
#'
#' @param x An `image` or `image` object.
#' @param ... Additional arguments passed to methods.
#' @export
image_data <- function(x, ...) {
  UseMethod("image_data")
}

#' @export
image_data.default <- function(x, ...) {
  if (isS4(x) && .hasSlot(x, ".Data")) {
    return(slot(x, ".Data"))
  }
  unclass(x)
}

#' @rdname image_data
#' @param type Storage format to extract:
#'   * `"same"` — unmodified array (raw bytes or double as-is; default, fastest).
#'   * `"numeric"` — integer `0..255` regardless of internal storage
#'     (useful for arithmetic like `image_data(img, "numeric") / 2` or index calculation).
#'   * `"normalized"` — double `0..1` (raw images are divided by 255, double images returned as-is).
#'   * `"raw"` — raw bytes `0..255` (double images are scaled and coerced to raw).
#' @export
image_data.image <- function(x, type = c("same", "numeric", "normalized", "raw"), ...) {
  type <- match.arg(type)
  u <- unclass(x)
  d <- dim(u)

  res <- if (type == "numeric") {
    if (is.raw(u)) as.integer(u) else as.integer(round(u * 255))
  } else if (type == "normalized") {
    if (is.raw(u)) as.numeric(u) / 255 else as.numeric(u)
  } else if (type == "raw") {
    if (is.raw(u)) u else as.raw(pmin(pmax(round(as.numeric(u) * 255), 0), 255))
  } else {
    u
  }

  dim(res) <- d
  res
}

#' Convert `image` between raw (0..255 uint8) and numeric (0..1 double)
#'
#' @param x An `image` object.
#' @return An `image` object converted to numeric storage.
#' @export
to_numeric <- function(x) {
  as_image(x, storage = "double")
}

#' @rdname to_numeric
#' @export
to_raw <- function(x) {
  as_image(x, storage = "raw")
}

#' Group generic Ops method for `image` objects
#'
#' @param e1 Left operand.
#' @param e2 Right operand.
#' @export
Ops.image <- function(e1, e2) {
  d1 <- if (is.array(e1)) dim(e1) else NULL
  d2 <- if (!missing(e2) && is.array(e2)) dim(e2) else NULL
  target_dim <- if (!is.null(d1)) d1 else d2

  cm1 <- attr(e1, "colormode")
  cm2 <- if (!missing(e2)) attr(e2, "colormode") else NULL
  cm <- if (!is.null(cm1)) cm1 else if (!is.null(cm2)) cm2 else "Color"

  cs1 <- attr(e1, "colspace")
  cs2 <- if (!missing(e2)) attr(e2, "colspace") else NULL
  cs <- if (!is.null(cs1)) cs1 else if (!is.null(cs2)) cs2 else "RGB"

  gm1 <- attr(e1, "gamma")
  gm2 <- if (!missing(e2)) attr(e2, "gamma") else NULL
  gm <- if (!is.null(gm1)) gm1 else if (!is.null(gm2)) gm2 else 1.0

  u1 <- unclass(e1)
  op_func <- get(.Generic, envir = baseenv())

  # Fast-path for relational comparisons (>, <, >=, <=, ==, !=) via C++ kernel
  if (.Generic %in% c("<", "<=", ">", ">=", "==", "!=") && !missing(e2) && is.numeric(e2) && length(e2) == 1) {
    op_code <- switch(.Generic, "<" = 1L, "<=" = 2L, ">" = 3L, ">=" = 4L, "==" = 5L, "!=" = 6L)
    res <- cpp_binary_threshold(e1, as.double(e2), op_code, FALSE)
    if (!is.null(target_dim)) dim(res) <- target_dim
    res <- as_image(res, colormode = cm, colspace = cs, gamma = gm, storage = "logical")
    return(res)
  }

  v1 <- if (is.raw(u1)) as.numeric(u1) / 255 else u1

  if (missing(e2)) {
    v2 <- NULL
  } else {
    u2 <- unclass(e2)
    v2 <- if (is.raw(u2)) as.numeric(u2) / 255 else u2
  }

  res <- if (is.null(v2)) op_func(v1) else op_func(v1, v2)

  if (!is.null(target_dim) && length(res) == prod(target_dim)) {
    dim(res) <- target_dim
    if (.Generic %in% c("<", "<=", ">", ">=", "==", "!=")) {
      res <- as_image(res, colormode = cm, colspace = cs, gamma = gm, storage = "logical")
    } else {
      class(res) <- c("image", "array")
      attr(res, "colormode") <- cm
      attr(res, "colspace") <- cs
      attr(res, "gamma") <- gm
    }
  }
  res
}

#' Compute and Plot Histogram for an `image` Object
#'
#' Computes an ultra-fast C++ pixel intensity histogram for `image` objects
#' of any storage mode (`raw`, `logical`, `double`, `integer`).
#'
#' @param x An `image` object.
#' @param nbins Integer specifying the number of bins (default 256 for uint8/normalized data).
#' @param plot Logical; if `TRUE` (default), plots the histogram.
#' @param col Color(s) used to fill the histogram bars.
#' @param main Main title for the plot.
#' @param xlab X-axis label.
#' @param ... Additional graphical arguments passed to `plot()`.
#'
#' @return A `histogram` S3 object (or list of `histogram` objects for multi-channel images) invisibly.
#' @export
hist.image <- function(x, nbins = 256, plot = TRUE, col = NULL, main = NULL, xlab = "Pixel Value", ...) {
  if (!is_image(x)) {
    cli::cli_abort("Input must be an {.cls image} object.")
  }

  u <- unclass(x)
  d <- dim(u)
  nch <- if (length(d) == 3 && d[3] >= 3) d[3] else 1L
  st <- typeof(u)

  counts_list <- image_histogram_cpp(x, as.integer(nbins))

  make_hist_obj <- function(counts, channel_name = NULL, is_raw = FALSE, is_lgl = FALSE) {
    nb <- length(counts)
    if (is_raw) {
      breaks <- (0:256) - 0.5
      mids <- 0:255
    } else if (is_lgl) {
      breaks <- c(-0.5, 0.5, 1.5)
      mids <- c(0, 1)
    } else {
      breaks <- seq(0, 1, length.out = nb + 1)
      mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
    }

    structure(
      list(
        breaks   = breaks,
        counts   = counts,
        density  = counts / (sum(counts) * diff(breaks[1:2])),
        mids     = mids,
        xname    = channel_name %||% "Pixel values",
        equidist = TRUE
      ),
      class = "histogram"
    )
  }

  is_raw_storage <- is.raw(u)
  is_lgl_storage <- is.logical(u)

  hists <- vector("list", nch)
  ch_names <- if (nch >= 3) c("Red", "Green", "Blue", paste0("Channel_", 4:nch)) else "Gray"

  for (ch in seq_len(nch)) {
    hists[[ch]] <- make_hist_obj(counts_list[[ch]], ch_names[ch], is_raw_storage, is_lgl_storage)
  }

  if (nch == 1) hists <- hists[[1]]

  if (isTRUE(plot)) {
    if (nch == 1) {
      bar_col <- col %||% if (is_lgl_storage) "lightgray" else "gray30"
      bar_main <- main %||% paste0("Histogram of Pixel Values (", st, ")")
      plot(hists, col = bar_col, main = bar_main, xlab = xlab, ...)
    } else {
      op <- graphics::par(mfrow = c(1, min(nch, 3)))
      on.exit(graphics::par(op))
      default_cols <- c("red", "green", "blue")
      for (ch in seq_len(min(nch, 3))) {
        bar_col <- if (!is.null(col)) col[min(ch, length(col))] else default_cols[ch]
        plot(hists[[ch]], col = bar_col, main = paste0("Channel: ", ch_names[ch]), xlab = xlab, ...)
      }
    }
  }

  invisible(hists)
}

#' Transpose an Image
#'
#' Transposes the spatial dimensions (width and height) of an image while preserving
#' multi-channel spatial structure, color modes, and metadata attributes.
#'
#' @param img,x An `image` object or 2D/3D array.
#' @return A transposed `image` object with swapped row/column dimensions.
#' @export
image_transpose <- function(img) {
  if (is.null(dim(img))) {
    cli::cli_abort("Data provided to {.fn image_transpose} must have dimensions (array or matrix).")
  }
  cpp_image_transpose(img)
}

#' @rdname image_transpose
#' @export
t.image <- function(x) {
  cpp_image_transpose(x)
}
