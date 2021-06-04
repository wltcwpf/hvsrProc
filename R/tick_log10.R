#' A function to add log10 ticks
#'
#' This function generates log10 ticks for the figure in x-axis / y-axis
#' @param x The array of x values
#' @param y The array of y values
#' @param log10_scale A string takes "x", "y", or "xy". It indicates the axix that will be plotted in log10 scale
#' @param tick_type A string takes "exp" (the ticks are expressed by 10^) or "lin" (the ticks are expressed by 100..)
#' @param minor_tick Binary, indicates if minor ticks are plotted
#' @importFrom graphics axis par
#' @export
log10_ticks <- function(x, y, log10_scale = 'x', tick_type = 'exp', minor_tick = TRUE) {

  if (grepl('x', log10_scale)) {
    x <- x[x > 0]
    x_range <- range(x)
    x_range_log10 <- c(ceiling(log10(x_range[1])), floor(log10(x_range[2])))
    # minor ticks
    if (minor_tick) {
      x_range_minor_log10 <- c(floor(log10(x_range[1])), ceiling(log10(x_range[2])))
      atx_minor <- outer(1:9, 10^(x_range_minor_log10[1]:x_range_minor_log10[2]))
      axis(1, at = atx_minor, labels = FALSE, tcl = par("tcl") * 0.5)
    }
    # major ticks
    if (tick_type == 'exp') {
      atx <- 10^(x_range_log10[1]:x_range_log10[2])
      labels <- sapply(x_range_log10[1]:x_range_log10[2], function(i) as.expression(bquote(10^ .(i))))
      axis(1, at = atx, labels = labels)
    }else if (tick_type == 'lin') {
      atx <- 10^(x_range_log10[1]:x_range_log10[2])
      axis(1, at = atx, labels = atx)
    }
  }

  if (grepl('y', log10_scale)) {
    y <- y[y > 0]
    y_range <- range(y)
    y_range_log10 <- c(ceiling(log10(y_range[1])), floor(log10(y_range[2])))
    # minor ticks
    if (minor_tick) {
      y_range_minor_log10 <- c(floor(log10(y_range[1])), ceiling(log10(y_range[2])))
      aty_minor <- outer(1:9, 10^(y_range_minor_log10[1]:y_range_minor_log10[2]))
      axis(2, at = aty_minor, labels = FALSE, tcl = par("tcl") * 0.5)
    }
    # major ticks
    if (tick_type == 'exp') {
      aty <- 10^(y_range_log10[1]:y_range_log10[2])
      labels <- sapply(y_range_log10[1]:y_range_log10[2], function(i) as.expression(bquote(10^ .(i))))
      axis(2, at = aty, labels = labels)
    }else if (tick_type == 'lin') {
      aty <- 10^(y_range_log10[1]:y_range_log10[2])
      axis(2, at = aty, labels = aty)
    }
  }
}
