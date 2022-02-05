#' A function for Fourier Amplitude Spectrum (FAS) calculation
#'
#' This function computes FAS (with and without 0 padding)
#' @param ts An array of time series
#' @param dt The time step
#' @param max_order An integer, the maximum factor after padding 0 at the end;
#' if max_order = 0, then no zero padding;
#' if max_order = 3, the length of ts becomes at least 4 times.
#' @return A list is returned with FFT of ts, FAS, the associated frequency of FAS, the associated phase
#' @keywords FAS
#' @importFrom stats fft
#' @export
fas_cal <- function(ts, dt, max_order = 0){
  npts <- length(ts)
  nfft <- ifelse(max_order > 0, 2^(floor(logb(npts, 2)) + max_order), npts)
  nNyq <- nfft / 2 + 1
  df <- 1 / (nfft * dt)
  ts <- c(ts, rep(0, nfft - npts))
  fft_ts <- stats::fft(ts)
  amp <- abs(fft_ts[1:nNyq]) * dt
  phase <- Arg(ts)[1:nNyq]
  freq <- ((1:nNyq) - 1) * df
  output <- list()
  output$fft_ts <- fft_ts
  output$amp <- amp
  output$freq <- freq
  output$phase <- phase
  output
}
