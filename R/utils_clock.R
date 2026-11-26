# ==============================================================================
# HIGH-RESOLUTION CLOCK AND BENCHMARKING
# pliman: Plant Image Analysis in R
# ==============================================================================

#' High-resolution clock for timing code
#'
#' @description
#' `pliman_clock()` reads a high-resolution monotonic clock implemented in C++
#' (`std::chrono::steady_clock`). It can be used in two ways:
#'
#' * Without `expr`, it returns the current clock reading. The difference between two
#'   readings is the elapsed time.
#' * With `expr`, it evaluates the expression `times` times and returns the time taken
#'   by each evaluation.
#'
#' To compare several expressions, including their memory use, see [benchmark()].
#'
#' @details
#' [base::proc.time()] and [base::Sys.time()] are enough for long tasks, but on Windows
#' [base::proc.time()] only advances in steps of about 10-16 ms, which is too coarse to time
#' operations that take a few milliseconds, such as a single deep learning inference.
#' `std::chrono::steady_clock` uses `QueryPerformanceCounter` on Windows and
#' `clock_gettime(CLOCK_MONOTONIC)` on Linux and macOS. Its resolution is usually
#' 100 ns or better (use `pliman_clock(info = TRUE)` to measure it on your machine),
#' and it never jumps backwards or forwards when the system time changes.
#'
#' A clock reading has no meaning by itself (its origin is arbitrary, often the time
#' the computer was started); only differences between readings are useful.
#'
#' Timings include the small overhead of calling the clock from R (typically below
#' 1 microsecond) and, with `expr`, of evaluating the expression in R.
#'
#' @param expr Optional R expression to time. If missing, the current clock reading
#'   is returned.
#' @param times Integer. Number of times `expr` is evaluated. Defaults to `1`.
#' @param unit Time unit of the result: `"ms"` (milliseconds, default), `"us"`
#'   (microseconds), `"ns"` (nanoseconds) or `"s"` (seconds).
#' @param info Logical. If `TRUE`, returns information about the clock (name, whether it
#'   is monotonic, tick length and measured resolution in nanoseconds) instead of a
#'   reading. Defaults to `FALSE`.
#' @param envir Environment in which `expr` is evaluated. Defaults to the calling
#'   environment.
#'
#' @return
#' * Without `expr`: a single number, the current clock reading in `unit`.
#' * With `expr`: a numeric vector of length `times` with the elapsed time of each
#'   evaluation in `unit`. The value returned by the last evaluation is stored in the
#'   `"value"` attribute.
#' * With `info = TRUE`: a list with elements `clock`, `is_steady`, `tick_ns` and
#'   `resolution_ns`.
#' @export
#' @examples
#' library(pliman)
#'
#' # Elapsed time between two readings
#' t0 <- pliman_clock()
#' x <- sum(sqrt(seq_len(1e5)))
#' pliman_clock() - t0  # milliseconds
#'
#' # Time an expression 10 times, in microseconds
#' tm <- pliman_clock(sum(sqrt(seq_len(1e5))), times = 10, unit = "us")
#' summary(as.numeric(tm))
#' attr(tm, "value")
#'
#' # Clock details and measured resolution
#' pliman_clock(info = TRUE)
#'
#' @examplesIf interactive()
#' # Time a single YOLO inference
#' img <- image_pliman("farm.jpg")
#' pliman_clock(image_detect_dl(img, plot = FALSE, verbose = FALSE), times = 5)
pliman_clock <- function(expr,
                         times = 1,
                         unit = c("ms", "us", "ns", "s"),
                         info = FALSE,
                         envir = parent.frame()) {
  unit <- match.arg(unit)
  if (isTRUE(info)) {
    return(pliman_clock_info_cpp())
  }
  scale <- .clock_scale(unit)
  if (missing(expr)) {
    return(pliman_clock_ms_cpp() * scale)
  }
  times <- .check_times(times)
  ex <- substitute(expr)
  elapsed <- numeric(times)
  value <- NULL
  for (i in seq_len(times)) {
    t0 <- pliman_clock_ms_cpp()
    value <- eval(ex, envir = envir)
    elapsed[i] <- pliman_clock_ms_cpp() - t0
  }
  out <- elapsed * scale
  attr(out, "value") <- value
  out
}


#' Benchmark R expressions: time and memory
#'
#' @description
#' `benchmark()` evaluates one or more R expressions several times and returns a data
#' frame comparing them: execution time (minimum, median, mean, maximum, standard
#' deviation, iterations per second and speed relative to the fastest expression) and
#' memory use (memory allocated by R and memory used by the whole R process).
#'
#' Time is measured with the high-resolution clock of [pliman_clock()].
#'
#' @details
#' For each expression, in the order given:
#' 1. If `gc = TRUE`, [base::gc()] is run so that garbage left by previous expressions does
#'    not affect the measurements.
#' 2. The expression is run `warmup` times without timing. This excludes one-time costs
#'    such as loading a deep learning model or compiling it for the GPU.
#' 3. If `memory = TRUE`, the expression is run once more under R's memory profiler
#'    ([utils::Rprofmem()]) to measure the memory it allocates. This run is not timed,
#'    because the profiler slows the code down.
#' 4. The expression is run `times` times, timing each run.
#'
#' Two kinds of memory are reported:
#' * `mem_alloc_mb`: memory allocated by R during one evaluation (sum of all R vector
#'   allocations, even if they are released afterwards). It requires R to be built
#'   with memory profiling (`capabilities("profmem")`), which is the case for the CRAN
#'   builds of R. Small allocations that R serves from its own memory pool are not
#'   counted. It is `NA` when memory profiling is not available or `memory = FALSE`.
#' * `rss_mb`, `rss_delta_mb` and `peak_mb`: memory used by the whole R process
#'   (resident set size), read from the operating system. This includes memory
#'   allocated outside R, by C/C++ code and libraries such as ONNX Runtime, which R
#'   does not see. `rss_delta_mb` is the change from just before the expression was
#'   first run (including warm-up) to the end of the timed runs; a large value usually
#'   means that the expression keeps something in memory, such as a cached model.
#'   `peak_mb` is the highest memory used by the process up to the end of the
#'   expression's runs; it never decreases, so it only reflects an expression if that
#'   expression set a new peak. GPU memory is not included.
#'
#' Timings depend on everything else running on the computer. Use enough `times` and
#' compare medians, which are less affected by occasional slow runs than means.
#'
#' @param ... Expressions to benchmark. Name them to set the labels shown in the
#'   result (e.g. `benchmark(cpu = f(x, engine = "cpu"), gpu = f(x, engine = "gpu"))`);
#'   unnamed expressions are labelled with their code.
#' @param times Integer. Number of timed runs of each expression. Defaults to `10`.
#' @param warmup Integer. Number of untimed runs of each expression before timing.
#'   Defaults to `1`.
#' @param unit Time unit of the result: `"ms"` (milliseconds, default), `"us"`
#'   (microseconds), `"ns"` (nanoseconds) or `"s"` (seconds).
#' @param memory Logical. If `TRUE` (default), measures the memory allocated by R with
#'   one extra, untimed run of each expression.
#' @param gc Logical. If `TRUE` (default), runs [base::gc()] before each expression.
#' @param check Logical. If `TRUE`, checks with [base::all.equal()] whether every expression
#'   returns the same result as the first one, and adds the column `same_result`.
#'   Defaults to `FALSE`.
#' @param envir Environment in which the expressions are evaluated. Defaults to the
#'   calling environment.
#'
#' @return A data frame of class `pliman_bench`, with one row per expression and the
#'   columns:
#' * `expression`: the label of the expression.
#' * `n`: number of timed runs.
#' * `min`, `median`, `mean`, `max`, `sd`: execution time, in `unit`.
#' * `itr_sec`: runs per second, based on the median time.
#' * `relative`: median time divided by the median time of the fastest expression
#'   (`1` for the fastest).
#' * `mem_alloc_mb`: memory allocated by R in one run, in MB.
#' * `rss_mb`: memory used by the R process after the runs, in MB.
#' * `rss_delta_mb`: change in the memory used by the R process, in MB.
#' * `peak_mb`: highest memory used by the R process so far, in MB.
#' * `same_result`: only when `check = TRUE`; whether the result is equal to the result
#'   of the first expression.
#'
#' The time of every run is stored in the `"timings"` attribute (a list with one
#' numeric vector per expression), and the unit in the `"unit"` attribute.
#' @seealso [plot.pliman_bench()] to plot the time distributions, [pliman_clock()], [yolo_benchmark()]
#' @export
#' @examples
#' library(pliman)
#' x <- runif(1e5)
#' benchmark(
#'   sqrt = sqrt(x),
#'   power = x^0.5,
#'   exp_log = exp(0.5 * log(x)),
#'   times = 20,
#'   check = TRUE
#' )
#'
#' @examplesIf interactive()
#' # Compare models and engines on the same image
#' img <- image_pliman("farm.jpg")
#' res <- benchmark(
#'   yolo26n_cpu = image_detect_dl(img, model = "yolo26n", engine = "cpu", plot = FALSE, verbose = FALSE),
#'   yolo26n_gpu = image_detect_dl(img, model = "yolo26n", engine = "gpu", plot = FALSE, verbose = FALSE),
#'   yolo26s_gpu = image_detect_dl(img, model = "yolo26s", engine = "gpu", plot = FALSE, verbose = FALSE),
#'   times = 20,
#'   warmup = 3
#' )
#' res
#'
#' # Time of every run of the first expression
#' attr(res, "timings")[[1]]
benchmark <- function(...,
                      times = 10,
                      warmup = 1,
                      unit = c("ms", "us", "ns", "s"),
                      memory = TRUE,
                      gc = TRUE,
                      check = FALSE,
                      envir = parent.frame()) {
  unit <- match.arg(unit)
  scale <- .clock_scale(unit)
  times <- .check_times(times)
  warmup <- as.integer(warmup)
  if (length(warmup) != 1L || is.na(warmup) || warmup < 0L) {
    cli::cli_abort("{.arg warmup} must be a non-negative integer.")
  }

  exprs <- as.list(substitute(list(...)))[-1L]
  if (length(exprs) == 0L) {
    cli::cli_abort("Provide at least one expression to benchmark.")
  }
  labels <- names(exprs)
  if (is.null(labels)) labels <- rep("", length(exprs))
  for (i in seq_along(exprs)) {
    if (!nzchar(labels[i])) {
      lab <- paste(deparse(exprs[[i]], width.cutoff = 500L), collapse = " ")
      labels[i] <- if (nchar(lab) > 50L) paste0(substr(lab, 1L, 47L), "...") else lab
    }
  }
  labels <- make.unique(labels)

  can_profmem <- isTRUE(memory) && isTRUE(capabilities("profmem"))
  if (isTRUE(memory) && !can_profmem) {
    cli::cli_alert_warning("This R build has no memory profiling ({.code capabilities(\"profmem\")} is {.val FALSE}); {.field mem_alloc_mb} will be {.val NA}.")
  }

  mb <- 1024^2
  n_expr <- length(exprs)
  timings <- vector("list", n_expr)
  names(timings) <- labels
  res <- data.frame(
    expression = labels,
    n = times,
    min = NA_real_, median = NA_real_, mean = NA_real_, max = NA_real_, sd = NA_real_,
    itr_sec = NA_real_, relative = NA_real_,
    mem_alloc_mb = NA_real_, rss_mb = NA_real_, rss_delta_mb = NA_real_, peak_mb = NA_real_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(check)) res$same_result <- NA
  first_value <- NULL

  for (j in seq_len(n_expr)) {
    ex <- exprs[[j]]
    if (isTRUE(gc)) invisible(base::gc(verbose = FALSE))
    mem0 <- pliman_process_memory_cpp()

    for (i in seq_len(warmup)) eval(ex, envir = envir)

    if (can_profmem) {
      res$mem_alloc_mb[j] <- .profmem_bytes(ex, envir) / mb
    }

    t <- numeric(times)
    value <- NULL
    for (i in seq_len(times)) {
      t0 <- pliman_clock_ms_cpp()
      value <- eval(ex, envir = envir)
      t[i] <- pliman_clock_ms_cpp() - t0
    }
    mem1 <- pliman_process_memory_cpp()

    t <- t * scale
    timings[[j]] <- t
    res$min[j] <- min(t)
    res$median[j] <- stats::median(t)
    res$mean[j] <- mean(t)
    res$max[j] <- max(t)
    res$sd[j] <- if (times > 1L) stats::sd(t) else NA_real_
    res$itr_sec[j] <- if (res$median[j] > 0) scale / (res$median[j] / 1000) else NA_real_
    res$rss_mb[j] <- mem1[["rss"]] / mb
    res$rss_delta_mb[j] <- (mem1[["rss"]] - mem0[["rss"]]) / mb
    res$peak_mb[j] <- mem1[["peak"]] / mb

    if (isTRUE(check)) {
      if (j == 1L) {
        first_value <- value
        res$same_result[j] <- TRUE
      } else {
        res$same_result[j] <- isTRUE(all.equal(first_value, value))
      }
    }
  }

  fastest <- min(res$median, na.rm = TRUE)
  res$relative <- if (is.finite(fastest) && fastest > 0) res$median / fastest else NA_real_

  attr(res, "timings") <- timings
  attr(res, "unit") <- unit
  attr(res, "warmup") <- warmup
  class(res) <- c("pliman_bench", "data.frame")
  res
}

#' Plot the time distribution of a benchmark
#'
#' Plots the time of every run of each expression benchmarked with [benchmark()]: a
#' horizontal boxplot per expression with the individual runs drawn on top. Expressions
#' are ordered from the fastest (top) to the slowest (bottom) median time, and the
#' median of each one is written at the right of the plot.
#'
#' @param x An object of class `pliman_bench`, returned by [benchmark()].
#' @param log Logical or `"auto"` (default). If `TRUE`, the time axis is on a log scale,
#'   which helps when the expressions differ a lot in speed. With `"auto"`, the log scale
#'   is used when the slowest run is more than 20 times slower than the fastest one.
#' @param points Logical. If `TRUE` (default), draws every run as a point on top of the
#'   boxplot.
#' @param col Fill colors of the boxes, one per expression, recycled if needed. If `NULL`
#'   (default), a qualitative palette is used.
#' @param show_median Logical. If `TRUE` (default), writes the median time of each
#'   expression at the right of the plot.
#' @param main,xlab Title and time-axis label. By default, `xlab` shows the time unit.
#' @param ... Further arguments passed to [graphics::boxplot()].
#' @return `x`, invisibly.
#' @method plot pliman_bench
#' @export
#' @examples
#' library(pliman)
#' x <- runif(1e5)
#' res <- benchmark(
#'   sqrt = sqrt(x),
#'   power = x^0.5,
#'   exp_log = exp(0.5 * log(x)),
#'   times = 30
#' )
#' plot(res)
#' plot(res, log = FALSE, points = FALSE)
plot.pliman_bench <- function(x,
                              log = "auto",
                              points = TRUE,
                              col = NULL,
                              show_median = TRUE,
                              main = "Benchmark",
                              xlab = NULL,
                              ...) {
  timings <- attr(x, "timings")
  if (is.null(timings) || length(timings) == 0L) {
    cli::cli_abort("{.arg x} has no timings. Use an object returned by {.fn benchmark}.")
  }
  unit <- attr(x, "unit")
  if (is.null(unit)) unit <- "ms"
  unit_lab <- switch(unit, ms = "milliseconds", us = "microseconds", ns = "nanoseconds", s = "seconds", unit)
  if (is.null(xlab)) xlab <- paste0("Time (", unit_lab, ")")

  # fastest expression on top: boxplot() draws the first group at the bottom
  meds <- vapply(timings, stats::median, numeric(1))
  ord <- rev(order(meds))
  timings <- timings[ord]
  meds <- meds[ord]
  n <- length(timings)

  if (is.null(col)) {
    col <- grDevices::hcl.colors(max(n, 3L), palette = "Set 2")[seq_len(n)]
    col <- rev(col)
  } else {
    col <- rep_len(col, n)[ord]
  }

  all_t <- unlist(timings, use.names = FALSE)
  all_t <- all_t[is.finite(all_t)]
  use_log <- if (identical(log, "auto")) {
    length(all_t) > 0 && min(all_t) > 0 && max(all_t) / min(all_t) > 20
  } else {
    isTRUE(log)
  }
  if (use_log && any(all_t <= 0)) {
    cli::cli_alert_warning("Some times are zero or negative; using a linear scale.")
    use_log <- FALSE
  }

  # left margin wide enough for the labels, right margin for the medians
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  lab_w <- max(graphics::strwidth(names(timings), units = "inches", cex = 0.9))
  med_txt <- paste(vapply(meds, function(m) format(signif(m, 3), big.mark = ","), character(1)), unit)
  med_w <- if (isTRUE(show_median)) max(graphics::strwidth(med_txt, units = "inches", cex = 0.8)) else 0
  mai <- graphics::par("mai")
  mai[2] <- max(mai[2], lab_w + 0.3)
  mai[4] <- max(0.2, med_w + 0.25)
  graphics::par(mai = mai)

  # time axis covering every run (with outline = FALSE, boxplot() would ignore outliers
  # that are still drawn as points). For horizontal boxplots it is set through `ylim`.
  # A user-supplied `ylim` (in `...`) takes precedence.
  dots <- list(...)
  time_lim <- if (!is.null(dots$ylim)) dots$ylim else range(all_t)
  dots$ylim <- NULL
  draw_box <- function(add) {
    args <- list(
      timings,
      horizontal = TRUE,
      log = if (use_log) "x" else "",
      col = grDevices::adjustcolor(col, alpha.f = 0.6),
      border = "gray25",
      outline = !isTRUE(points),
      las = 1,
      cex.axis = 0.9,
      xlab = if (add) "" else xlab,
      main = if (add) "" else main,
      add = add,
      axes = !add,
      ylim = time_lim
    )
    do.call(graphics::boxplot, utils::modifyList(args, dots))
  }
  # draw once to set up the axes, add the grid, then redraw the boxes over the grid
  draw_box(add = FALSE)
  graphics::grid(nx = NULL, ny = NA, col = "gray88", lty = 1)
  draw_box(add = TRUE)

  if (isTRUE(points)) {
    for (i in seq_len(n)) {
      ti <- timings[[i]]
      yj <- i + stats::runif(length(ti), -0.18, 0.18)
      graphics::points(ti, yj, pch = 16, cex = 0.7,
                       col = grDevices::adjustcolor(col[i], red.f = 0.6, green.f = 0.6, blue.f = 0.6, alpha.f = 0.7))
    }
  }

  if (isTRUE(show_median)) {
    usr <- graphics::par("usr")
    xr <- if (use_log) 10^usr[2] else usr[2]
    graphics::text(xr, seq_len(n), labels = med_txt, pos = 4, xpd = TRUE, cex = 0.8, col = "gray20")
  }
  invisible(x)
}

#' @export
print.pliman_bench <- function(x, digits = 3, ...) {
  unit <- attr(x, "unit")
  cli::cli_rule(left = "Benchmark ({nrow(x)} expression{?s}, {x$n[1]} run{?s} each)")
  cli::cli_alert_info("Times in {.field {unit}} | memory in {.field MB} | warm-up runs: {attr(x, 'warmup')}")
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(df, "timings") <- NULL
  attr(df, "unit") <- NULL
  attr(df, "warmup") <- NULL
  num <- vapply(df, is.numeric, logical(1)) & names(df) != "n"
  df[num] <- lapply(df[num], signif, digits = digits)
  print(df, row.names = FALSE)
  invisible(x)
}

# --- internal helpers ---------------------------------------------------------

.clock_scale <- function(unit) {
  switch(unit, ms = 1, us = 1e3, ns = 1e6, s = 1e-3)
}

.check_times <- function(times) {
  times <- as.integer(times)
  if (length(times) != 1L || is.na(times) || times < 1L) {
    cli::cli_abort("{.arg times} must be a positive integer.")
  }
  times
}

# Bytes allocated by R while evaluating `ex` once, from utils::Rprofmem().
# Each line of the profile starts with the size of a vector allocation in bytes;
# "new page" lines (small vectors served from R's own pool) are not counted.
.profmem_bytes <- function(ex, envir) {
  tf <- tempfile("pliman_profmem_", fileext = ".out")
  on.exit({
    utils::Rprofmem(NULL)
    unlink(tf)
  }, add = TRUE)
  utils::Rprofmem(tf, threshold = 0)
  eval(ex, envir = envir)
  utils::Rprofmem(NULL)
  if (!file.exists(tf)) return(NA_real_)
  lines <- readLines(tf, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub("^([0-9]+) ?:.*$", "\\1", lines)))
  sum(bytes, na.rm = TRUE)
}
