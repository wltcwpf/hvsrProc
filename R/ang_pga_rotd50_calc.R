#' Rotation angle of PGA RotD50 Calculation
#'
#' This function calculates the rotation angle (relative to the azimuth of h1) of PGA RotD50
#' @param h1 A vector of time series for the 1st horizontal component
#' @param h2 A vector of time series for the 2nd horizontal component
#' @return The rotation angle for PGA RotD50
#' @export
ang_pga_rotd50_calc <- function(h1, h2) {
  pga_s <- rep(NA, 180)
  for (i in 0:179) {
    ts_tmp <- h1 * cos(i / 180 * pi) + h2 * sin(i / 180 * pi)
    pga_s[i + 1] <- max(abs(ts_tmp))
  }
  which.min(abs(pga_s - stats::median(pga_s))) - 1
}
