#' A function to plot a figure in log10 scale
#'
#' This function generates a figure in log10 scale
#' @param x The array of x values
#' @param y The array of y values
#' @param log10_scale A string takes "x", "y", or "xy". It indicates the axix that will be plotted in log10 scale
#' @param tick_type A string takes "exp" (the ticks are expressed by 10^) or "lin" (the ticks are expressed by 100..)
#' @param minor_tick Binary, indicates if minor ticks are plotted
#' @param xlab A title for the x axis
#' @param ylab A title for the y axis
#' @param type What type of plot should be drawn. The setting is the same as plot()
#' @param main An overall title for the plot
#' @param xlim The limit range of x axis
#' @param ylim The limit range of y axis
#' @importFrom graphics axis par
#' @export
plot_log10 <- function(x, y, log10_scale = 'x', tick_type = 'lin', minor_tick = TRUE,
                       xlab = NULL, ylab = NULL, type = NULL, main = NULL, xlim = NULL,
                       ylim = NULL) {

  if (log10_scale == 'x') {

    plot(x, y, type = type, xaxt = 'n', xlab = xlab, ylab = ylab, main = main, log = 'x',
         xlim = xlim, ylim = ylim)

    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type, minor_tick = minor_tick)
  }

  if (log10_scale == 'y') {

    plot(x, y, type = type, yaxt = 'n', xlab = xlab, ylab = ylab, main = main, log = 'y',
         xlim = xlim, ylim = ylim)

    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type, minor_tick = minor_tick)
  }

  if (log10_scale == 'xy') {

    plot(x, y, type = type, yaxt = 'n', xaxt = 'n', xlab = xlab, ylab = ylab, main = main,
         log = 'xy', xlim = xlim, ylim = ylim)

    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type, minor_tick = minor_tick)
  }
}
