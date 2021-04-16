#' Low-pass Butterworth filter
#'
#' This function applies low-pass Butterworth filter to process time series (for both causal and acausal)
#' @param ts An array of time series
#' @param dt The time step
#' @param fc The corner frequency for low-pass/high-pass
#' @param nPole The pole parameter for the order, where the order is 2 * nPole. Positive for low-pass; Negative for high-pass
#' Generally, nPole = 4 for low-pass and nPole = -5 for high-pass
#' @param is_causal A binary to indicator if causal is applied. Apply causal if is_causal = TRUE and apply acausal if is_causal = FALSE
#' @param order_zero_padding The order needs to be added for zeroes padding at the end of recordings to
#' increase the number of data points to a power of 2 (Kishida T., Darragh R.B., Chiou B.S.J., Bozorgnia Y., Mazzoni S.,
#' Contreras V., Boroschek R., Rojas F., Stewart J.P., 2020.
#' Chapter 3: Ground Motions and Intensity Measures, in Data Resources for NGA-Subduction Project,
#' PEER Report 2020/02, J.P. Stewart (editor), Pacific Earthquake Engineering Research Center, UC Berkeley (headquarters).)
#' @return The filtered time series
#' @export
bw_pass <- function(ts, dt, fc, nPole = 4, is_causal = FALSE, order_zero_padding = 4) {

  npts <- length(ts)
  nfft <- 2^(floor(logb(npts, 2)) + order_zero_padding)
  nNyq <- nfft / 2 + 1
  df <- 1.0 / (nfft * dt)
  freq <- ((1:nNyq) - 1) * df
  # zero padding at the end of recordings
  ts <- c(ts, rep(0, nfft - npts))
  ts_fft <- stats::fft(ts)

  # transfer response
  Tmp_resp <- bw_resp(freq[1:nNyq], fc, nPole, is_causal)
  ts_fft[1:nNyq] <- ts_fft[1:nNyq] * Tmp_resp

  # Calculate inverse Fourier spectrum: keep the original time length
  ts_fft[1] <- complex(real = Re(ts_fft[1]), imaginary = 0)
  ts_fft[nNyq] <- complex(real = Re(ts_fft[nNyq]), imaginary = 0)
  ts_fft[nfft + 2 - (2:(nfft / 2))] <- Conj(ts_fft[2:(nfft / 2)])
  ts_flt <- Re(stats::fft(ts_fft, inverse = T))[1:npts] / nfft

  # filtered FAS
  ts_flt_amp <- abs(ts_fft[1:nNyq]) * dt

  res <- list()
  res$flt_ts <- ts_flt
  res$flt_amp <- ts_flt_amp
  res$flt_resp <- Tmp_resp
  res$flt_fft <- ts_fft
  return(res)
}
