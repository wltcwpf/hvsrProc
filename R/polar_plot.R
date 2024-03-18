#' A function to generate polar plot
#'
#' This function generates a polar plot figure (heatmap)
#' @param x The array of x values, e.g., frequency array
#' @param y The array of y values, e.g., azimuth angles
#' @param z The numeric matrix containing the values to be plotted. It must have the dimension of len(x) by len(y).
#' @param log10_scale A string takes "x", "y", or "xy". It indicates the axix that will be plotted in log10 scale
#' @param tick_type A string takes "exp" (the ticks are expressed by 10^) or "lin" (the ticks are expressed by 100..)
#' @param minor_tick Binary, indicates if minor ticks are plotted
#' @param xlab A title for the x axis
#' @param ylab A title for the y axis
#' @importFrom graphics axis par image rect text
#' @importFrom grDevices terrain.colors
#' @importFrom stats quantile
#' @export
polar_plot <- function(x, y, z, log10_scale = "x", tick_type = "lin",
                       minor_tick = TRUE,
                       xlab = 'Freq (Hz)', ylab = 'Azimuth (deg)') {
  color_s <- terrain.colors(100)
  if (log10_scale == 'x') {
    image(x = x, y = y, z = z, xaxt = 'n', log = 'x', col = color_s,
          xlab = xlab, ylab = ylab)
    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type,
                minor_tick = minor_tick)
    xlim <- 10^c(4*(range(log10(x))[2] - range(log10(x))[1]) / 5 + range(log10(x))[1],
                 4*(range(log10(x))[2] - range(log10(x))[1]) / 5 +
                   (range(log10(x))[2] - range(log10(x))[1]) / 20 + range(log10(x))[1])
    ylim <- c(quantile(y, 0.25), quantile(y, 0.75))
  } else if (log10_scale == 'y') {
    image(x = x, y = y, z = z, yaxt = 'n', log = 'y', col = color_s,
          xlab = xlab, ylab = ylab)
    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type,
                minor_tick = minor_tick)
    xlim <- c(quantile(x, 0.75), quantile(x, 0.85))
    ylim <- c(2*(range(log10(y))[2] - range(log10(y))[1]) / 5,
              4*(range(log10(y))[2] - range(log10(y))[1]) / 5)
  } else if (log10_scale == 'xy') {
    image(x = x, y = y, z = z, xaxt = 'n', yaxt = 'n', log = 'xy', col = color_s,
          xlab = xlab, ylab = ylab)
    log10_ticks(x, y, log10_scale = log10_scale, tick_type = tick_type,
                minor_tick = minor_tick)
    xlim <- c(4*(range(log10(x))[2] - range(log10(x))[1]) / 5,
              4*(range(log10(x))[2] - range(log10(x))[1]) / 5 +
                (range(log10(x))[2] - range(log10(x))[1]) / 10)
    ylim <- c(2*(range(log10(y))[2] - range(log10(y))[1]) / 5,
              4*(range(log10(y))[2] - range(log10(y))[1]) / 5)
  } else {
    image(x = x, y = y, z = z, xlab = xlab, ylab = ylab)
    xlim <- c(quantile(x, 0.75), quantile(x, 0.85))
    ylim <- c(quantile(y, 0.25), quantile(y, 0.75))
  }

  # add color bar
  len <- length(color_s)
  y_seq <- seq(ylim[1], ylim[2], length.out = len + 1)
  rect(rep(xlim[1], len), y_seq[1:len],
       rep(xlim[2], len), y_seq[-1], col = color_s,
       border = color_s)
  rect(xlim[1], ylim[1], xlim[2], ylim[2],
       border = "black")
  text(x = rep(xlim[2], 2), y = c(ylim[1], ylim[2]),
       labels = round(range(z), digits = 2), pos = 4)
}
