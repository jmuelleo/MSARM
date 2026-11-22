
#' @title MSARM.plot
#' @description MSARM.plot creates forecast plots for MSARM.fit results
#' @param res_MSARM.predict An object saving the results from MSARM.predict
#' @param conf Set to TRUE if MSARM.predict was implemented with the option boot = TRUE, then bootstrap confidence intervals are plotted
#' @param start Start time. Defaults to the time series start when available.
#' @param freq Series frequency. Defaults to the input frequency when available.
#'
#' @return Invisibly returns NULL. Draws the forecast and optional prediction bands.
#' @export
#'
MSARM.plot <- function(res_MSARM.predict, conf = FALSE,
                       start = c(1, 1), freq = 1) {
  observed <- res_MSARM.predict$Y_T
  if (missing(start) && stats::is.ts(observed)) {
    start <- stats::start(observed)
  }
  if (missing(freq) && stats::is.ts(observed)) {
    freq <- stats::frequency(observed)
  }
  n <- length(observed)
  future <- res_MSARM.predict$fcast
  intervals <- res_MSARM.predict$fcast_boot
  if (conf && is.null(intervals)) {
    stop("Confidence bands require MSARM.predict(..., boot = TRUE).",
         call. = FALSE)
  }
  times <- as.numeric(stats::time(stats::ts(c(observed, future),
                                            start = start, frequency = freq)))
  ahead <- times[n + seq_along(future)]
  draw <- function(first) {
    shown <- seq.int(first, n)
    limits <- range(observed[shown], future, if (conf) intervals)
    graphics::plot(times[shown], observed[shown], type = "l", xlab = "Time",
                   ylab = "Value", xlim = range(times[first], ahead),
                   ylim = limits)
    center <- future
    if (conf) {
      count <- (ncol(intervals) - 1L) / 2L
      shades <- grDevices::gray(seq(0.85, 0.5, length.out = count))
      for (i in seq_len(count)) {
        graphics::polygon(c(ahead, rev(ahead)),
                          c(intervals[, i + 1L],
                            rev(intervals[, ncol(intervals) - i + 1L])),
                          col = shades[i], border = NA)
      }
      center <- intervals[, 1L]
    }
    graphics::lines(c(times[n], ahead), c(observed[n], center),
                    col = "blue", lwd = 2)
  }
  draw(1L)
  if (conf) {
    draw(max(1L, n - 29L))
  }
  invisible(NULL)
}
