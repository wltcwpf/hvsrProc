#' HVSR calculation for the given window
#'
#' This function calculates HVSR and polar HVSR for the given window
#' @param i_win The index of running window
#' @param h1_wins A list of 1st horizontal component time series. Each element is for one window
#' @param h2_wins A list of 2nd horizontal component time series. Each element is for one window
#' @param v_wins A list of vertical horizontal component time series. Each element is for one window
#' @param dt The time step
#' @param horizontal_comb The parameter specifies the combination of two horizontal components. ps_RotD50: rotated combination at the angle where PGA is median; geometric_mean: geometric mean (sqrt(h1(f) * h2(f))); squared_average: squared average (sqrt((h1(f)^2 + h2(f)^2)/2))
#' @param freq_hv_mean The target frequencys of HVSR
#' @param polar_curves_flag The flag indicates if polar curves are calculated
#' @param freq_polar The target frequencys of polar HVSR
#' @param deg_increment The degree increment for HVSR polar curves
#' @return The HVSR for the given window
#' @importFrom stats approx fft
#' @export
hvsr_win_calc <- function(i_win, h1_wins, h2_wins, v_wins, dt, horizontal_comb = 'geometric_mean', freq_hv_mean,
                          polar_curves_flag = TRUE, freq_polar, deg_increment = 10) {

  h1_sub <- h1_wins[[ i_win ]]
  h2_sub <- h2_wins[[ i_win ]]
  v_sub <- v_wins[[ i_win ]]
  res <- list()

  # combined horizontal curve:
  if (horizontal_comb == 'ps_RotD50') {
    # find rotation angle for PGA RotD50 and calculate RotD50 HVSR
    angle_idx <- ang_pga_rotd50_calc(h1 = h1_sub, h2 = h2_sub)
    h_combin <- h1_sub * cos(pi * angle_idx / 180) + h2_sub * sin(pi * angle_idx / 180)
    fas_h <- fas_cal(ts = h_combin, dt = dt)
    freq <- fas_h$freq
    fas_h <- fas_h$amp
  } else if (horizontal_comb == 'squared_average') {
    fas_h1 <- fas_cal(ts = h1_sub, dt = dt)
    freq <- fas_h1$freq
    fas_h2 <- fas_cal(ts = h2_sub, dt = dt)
    fas_h <- sqrt((fas_h1$amp^2 + fas_h2$amp^2)/2)
  } else {
    # geometric mean
    fas_h1 <- fas_cal(ts = h1_sub, dt = dt)
    freq <- fas_h1$freq
    fas_h2 <- fas_cal(ts = h2_sub, dt = dt)
    fas_h <- sqrt(fas_h1$amp * fas_h2$amp)
  }

  h_smooth <- ko_smooth(freq = freq, amp = fas_h)
  fas_v <- fas_cal(ts = v_sub, dt = dt)
  v_smooth <- ko_smooth(freq = fas_v$freq, amp = fas_v$amp)
  hv_ratio <- h_smooth / v_smooth

  fas_h1 <- fas_cal(ts = h1_sub, dt = dt)
  fas_h2 <- fas_cal(ts = h2_sub, dt = dt)
  h1_smooth <- ko_smooth(freq = freq, amp = fas_h1)
  h2_smooth <- ko_smooth(freq = freq, amp = fas_h2)

  # interpolation
  hv_ratio <- approx(freq, hv_ratio, freq_hv_mean, rule = 2)$y
  res$h1_smooth <- h1_smooth
  res$h2_smooth <- h2_smooth
  res$v_smooth <- v_smooth
  res$hv_ratio <- hv_ratio


  # polar curve:
  if (polar_curves_flag) {
    polar_degs <- seq(0, 179, by = deg_increment)
    polar_hv_ratio <- matrix(data = NA, nrow = length(freq_polar), ncol = length(polar_degs))
    h1_fft <- fft(h1_sub)
    h2_fft <- fft(h2_sub)
    for (i in 1:length(polar_degs)) {
      if (i == 1) {
        fas_v <- fas_cal(ts = v_sub, dt = dt)
        freq <- fas_v$freq
        v_smooth <- ko_smooth(freq = fas_v$freq, amp = fas_v$amp)
      }
      angle_idx <- polar_degs[i]
      h_fft <- h1_fft * cos(pi * angle_idx / 180) + h2_fft * sin(pi * angle_idx / 180)
      h_smooth <- ko_smooth(freq = freq, amp = abs(h_fft) * dt)
      hv_ratio <- h_smooth / v_smooth
      polar_hv_ratio[, i] <- approx(freq, hv_ratio, freq_polar, rule = 2)$y
    }
    res$polar_hv_ratio <- polar_hv_ratio
  }
  return(res)
}
