#' Pre-processing function
#'
#' This function pre-process time series (noise). The pre-processing includes mean removal, Taper, and Butterworth filter
#' @param ts An array of time series
#' @param dt The time step
#' @param dc_flag The flag indicates if mean removal is applied
#' @param taper_flag The flag indicates if Taper is applied
#' @param t_front A number, the percentage of taperring for the beginning of the time series
#' @param t_end A number, the percentage of taperring for the end of the time series
#' @param filter_flag The flag indicates if filter is applied
#' @param fc The corner frequency for low-pass/high-pass
#' @param nPole The parameter for the order, where the order is 2 * nPole. Positive for low-pass; Negative for high-pass
#' Generally, nPole = 4 for low-pass and nPole = -5 for high-pass
#' @param is_causal A binary to indicator if causal is applied. Apply causal if is_causal = TRUE and apply acausal if is_causal = FALSE
#' @param order_zero_padding The order needs to be added for zeroes padding at the end of recordings to
#' increase the number of data points to a power of 2 (Kishida T., Darragh R.B., Chiou B.S.J., Bozorgnia Y., Mazzoni S.,
#' Contreras V., Boroschek R., Rojas F., Stewart J.P., 2020.
#' Chapter 3: Ground Motions and Intensity Measures, in Data Resources for NGA-Subduction Project,
#' PEER Report 2020/02, J.P. Stewart (editor), Pacific Earthquake Engineering Research Center, UC Berkeley (headquarters).)
#' @return The filtered time series
#' @export
pre_proc <- function(ts, dt, dc_flag = TRUE, taper_flag = TRUE, t_front, t_end, filter_flag = TRUE,
                     fc, nPole, is_causal, order_zero_padding) {

  # mean removal
  if (dc_flag)
    ts <- mean_removal(ts)

  # Taper
  if (taper_flag)
    ts <- taper(ts, t_front = t_front, t_end = t_end)

  # Filter
  if (filter_flag) {
    res <- bw_pass(ts = ts, dt = dt, fc = fc, nPole = nPole, is_causal = is_causal,
                   order_zero_padding = order_zero_padding)
    ts <- res$flt_ts
  }

  # return
  ts
}
