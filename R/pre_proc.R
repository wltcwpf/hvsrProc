#' Pre-processing function
#'
#' This function pre-process time series (noise). The pre-processing includes mean removal, Taper, and Butterworth filter
#' @param ts An array of time series
#' @param dt The time step
#' @param detrend The indicator specifies which detrend method is used. 0: no detrend; 1: mean removal; 2: linear trend removal.
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
#' @importFrom stats lm
#' @importFrom stats weighted.mean
#' @export
pre_proc <- function(ts, dt, detrend = 1, taper_flag = TRUE, t_front, t_end, filter_flag = TRUE,
                     fc, nPole, is_causal, order_zero_padding) {

  # detrend
  if (detrend == 0) {
    ts <- ts
  } else if (detrend == 1) {
    # ts <- mean_removal(ts)

    win <- taper(rep(1, length(ts)), t_front = t_front, t_end = t_end)
    ts_avg <- weighted.mean(ts, win)
    ts <- (ts - ts_avg)

  } else if (detrend == 2) {
    ts <- lm(ts ~ seq(1, length(ts)))$residuals
  } else if (detrend == 6) {
    x <- seq(1, length(ts))
    ts <- lm(ts ~ I(x^5) + I(x^4) + I(x^3) + I(x^2) + x)$residuals
  }

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
