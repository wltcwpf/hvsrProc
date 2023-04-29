#' Frequency-domain plot and window selection
#'
#' This function plots each window in frequency domain and let users to manually select windows for rejection
#' @param hvsr_list A list. Each element is a window HVSR
#' @param freq_hv_mean The frequency
#' @param freq_min The minimum frequency of interest
#' @param freq_max The maximum frequency of interest
#' @param pre_filter_flag The flag indicate if pre-filter is applied
#' @param pre_filter_hpass_fc High-pass corner frequency applied in pre-filter
#' @param pre_filter_lpass_fc Low-pass corner frequency applied in pre-filter
#' @param filter_flag The flag indicate if filter is applied
#' @param hpass_fc High-pass corner frequency applied in filter
#' @param lpass_fc Low-pass corner frequency applied in filter
#' @param distribution The distribution assumption on the HVSR amplitudes. It can take "normal" or "log_normal".
#' @return The lefted indices of windows after frequency-domain rejection
#' @importFrom graphics plot abline lines legend locator text
#' @importFrom grDevices dev.off adjustcolor
#' @importFrom stats sd median mad IQR
#' @export
fd_plt_select <- function(hvsr_list, robust_est = FALSE, freq_hv_mean, freq_min, freq_max,
                          pre_filter_flag, pre_filter_hpass_fc, pre_filter_lpass_fc,
                          filter_flag, hpass_fc, lpass_fc, distribution) {

  visual_freq_min <- max(c(freq_min, min(freq_hv_mean))) / 2
  x_range <- c(visual_freq_min, min(c(freq_max, max(freq_hv_mean))))
  hvsr_mat <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(hvsr_list))
  h1_mat <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(hvsr_list))
  h2_mat <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(hvsr_list))
  v_mat <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(hvsr_list))
  for (i in 1:length(hvsr_list)) {
    hvsr_mat[, i] <- hvsr_list[[ i ]]$hv_ratio
    h1_mat[, i] <- hvsr_list[[ i ]]$h1_smooth
    h2_mat[, i] <- hvsr_list[[ i ]]$h2_smooth
    v_mat[, i] <- hvsr_list[[ i ]]$v_smooth
  }
  y_range <- range(hvsr_mat[freq_hv_mean >= visual_freq_min, ])

  ind_y_range <- range(c(h1_mat[freq_hv_mean >= visual_freq_min, ],
                         h2_mat[freq_hv_mean >= visual_freq_min, ],
                         v_mat[freq_hv_mean >= visual_freq_min, ]))
  par(mfrow = c(1, 2))
  plot(x_range, c(ind_y_range[1], ind_y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
       ylab = 'FAS', log = 'xy', xaxt = 'n', yaxt = 'n')
  log10_ticks(x = x_range, y = c(ind_y_range[1], ind_y_range[2]*1.5), log10_scale = 'xy', tick_type = 'lin')
  abline(v = x_range[2], lwd = 2)
  if (pre_filter_flag) {
    if (!is.na(pre_filter_hpass_fc))
      abline(v = pre_filter_hpass_fc, lwd = 2, col = 'purple', lty = 2)
    if (!is.na(pre_filter_lpass_fc))
      abline(v = pre_filter_lpass_fc, lwd = 2, col = 'purple', lty = 2)
  }
  if (filter_flag) {
    if (!is.na(hpass_fc))
      abline(v = hpass_fc, lwd = 2, col = 'purple')
    if (!is.na(lpass_fc))
      abline(v = lpass_fc, lwd = 2, col = 'purple')
  }
  for(i_plot in 1:ncol(hvsr_mat)){
    lines(freq_hv_mean, h1_mat[,i_plot], col = adjustcolor('red', alpha.f = 0.2), lwd = 0.5)
    lines(freq_hv_mean, h2_mat[,i_plot], col = adjustcolor('green', alpha.f = 0.2), lwd = 0.5)
    lines(freq_hv_mean, v_mat[,i_plot], col = adjustcolor('blue', alpha.f = 0.2), lwd = 0.5)
  }
  if (robust_est) {
    lines(freq_hv_mean, apply(h1_mat, 1, median), col = 'red', lwd = 2)
    lines(freq_hv_mean, apply(h2_mat, 1, median), col = 'green', lwd = 2)
    lines(freq_hv_mean, apply(v_mat, 1, median), col = 'blue', lwd = 2)
  } else {
    lines(freq_hv_mean, apply(h1_mat, 1, function(x) exp(mean(log(x)))), col = 'red', lwd = 2)
    lines(freq_hv_mean, apply(h2_mat, 1, function(x) exp(mean(log(x)))), col = 'green', lwd = 2)
    lines(freq_hv_mean, apply(v_mat, 1, function(x) exp(mean(log(x)))), col = 'blue', lwd = 2)
  }
  legend('topright', legend = c('H1', 'H2', 'V'), col = c('red', 'green', 'blue'), lty = 1,
         lwd = 2)


  plot(x_range, c(y_range[1], y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
       ylab = 'HVSR', log = 'x', xaxt = 'n')
  log10_ticks(x = x_range, y = c(y_range[1], y_range[2]*1.5), log10_scale = 'x', tick_type = 'lin')
  abline(v = x_range[2], lwd = 2)
  if (pre_filter_flag) {
    if (!is.na(pre_filter_hpass_fc))
      abline(v = pre_filter_hpass_fc, lwd = 2, col = 'purple', lty = 2)
    if (!is.na(pre_filter_lpass_fc))
      abline(v = pre_filter_lpass_fc, lwd = 2, col = 'purple', lty = 2)
  }
  if (filter_flag) {
    if (!is.na(hpass_fc))
      abline(v = hpass_fc, lwd = 2, col = 'purple')
    if (!is.na(lpass_fc))
      abline(v = lpass_fc, lwd = 2, col = 'purple')
  }
  for(i_plot in 1:ncol(hvsr_mat)){
    lines(freq_hv_mean, hvsr_mat[,i_plot])
  }

  if (robust_est){
    hvsr_mean <- apply(hvsr_mat, 1, median)
    hvsr_sd <- apply(hvsr_mat, 1, mad) # /sqrt(ncol(hvsr_mat))
    hvsr_sd1 <- apply(hvsr_mat, 1, IQR)
    lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
    lines(freq_hv_mean, hvsr_mean - hvsr_sd, lwd = 2, col = 'blue', lty = 2)
    lines(freq_hv_mean, hvsr_mean + hvsr_sd, lwd = 2, col = 'blue', lty = 2)
    legend('top', 'Please select HVSR curves for removal.', bg = 'white')
  } else {
    if (distribution == 'normal') {
      hvsr_mean <- apply(hvsr_mat, 1, mean)
      hvsr_sd <- apply(hvsr_mat, 1, sd) # /sqrt(ncol(hvsr_mat))
      lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
      lines(freq_hv_mean, hvsr_mean - hvsr_sd, lwd = 2, col = 'blue', lty = 2)
      lines(freq_hv_mean, hvsr_mean + hvsr_sd, lwd = 2, col = 'blue', lty = 2)
      legend('top', 'Please select HVSR curves for removal.', bg = 'white')
    } else if (distribution == 'log_normal') {
      hvsr_mat[hvsr_mat <= 0] <- 10e-5  # avoid log(0 or negative)
      hvsr_mean <- exp(apply(log(hvsr_mat), 1, mean))
      hvsr_sd <- apply(log(hvsr_mat), 1, sd) # /sqrt(ncol(hvsr_mat))
      lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
      lines(freq_hv_mean, hvsr_mean / exp(hvsr_sd), lwd = 2, col = 'blue', lty = 2)
      lines(freq_hv_mean, hvsr_mean * exp(hvsr_sd), lwd = 2, col = 'blue', lty = 2)
      legend('top', 'Please select HVSR curves for removal.', bg = 'white')
    }
  }

  p <- locator(1)
  idx_select <- 1:ncol(hvsr_mat)
  idx_remove <- c()
  while(p$x <= x_range[2]){
    while(p$x <= x_range[2]){
      short_dist <- sapply(idx_select, function(i){
        min((log(freq_hv_mean) - log(p$x))^2 + (hvsr_mat[, i] - p$y)^2, na.rm = TRUE)
      })
      idx_temp <- which.min(short_dist)
      lines(freq_hv_mean, hvsr_mat[, idx_select[idx_temp]], col = 'gray')
      idx_remove <- c(idx_remove, idx_select[idx_temp])
      idx_select <- idx_select[-idx_temp]
      p <- locator(1)
    }
    y_range <- range(hvsr_mat[freq_hv_mean >= visual_freq_min, idx_select])

    par(mfrow = c(1, 2))
    plot(x_range, c(ind_y_range[1], ind_y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
         ylab = 'FAS', log = 'xy', xaxt = 'n', yaxt = 'n')
    log10_ticks(x = x_range, y = c(ind_y_range[1], ind_y_range[2]*1.5), log10_scale = 'xy', tick_type = 'lin')
    abline(v = x_range[2], lwd = 2)
    if (pre_filter_flag) {
      if (!is.na(pre_filter_hpass_fc))
        abline(v = pre_filter_hpass_fc, lwd = 2, col = 'purple', lty = 2)
      if (!is.na(pre_filter_lpass_fc))
        abline(v = pre_filter_lpass_fc, lwd = 2, col = 'purple', lty = 2)
    }
    if (filter_flag) {
      if (!is.na(hpass_fc))
        abline(v = hpass_fc, lwd = 2, col = 'purple')
      if (!is.na(lpass_fc))
        abline(v = lpass_fc, lwd = 2, col = 'purple')
    }
    for(i_plot in idx_select){
      lines(freq_hv_mean, h1_mat[,i_plot], col = adjustcolor('red', alpha.f = 0.2), lwd = 0.5)
      lines(freq_hv_mean, h2_mat[,i_plot], col = adjustcolor('green', alpha.f = 0.2), lwd = 0.5)
      lines(freq_hv_mean, v_mat[,i_plot], col = adjustcolor('blue', alpha.f = 0.2), lwd = 0.5)
    }

    if (robust_est) {
      lines(freq_hv_mean, apply(h1_mat[, idx_select], 1, median), col = 'red', lwd = 2)
      lines(freq_hv_mean, apply(h2_mat[, idx_select], 1, median), col = 'green', lwd = 2)
      lines(freq_hv_mean, apply(v_mat[, idx_select], 1, median), col = 'blue', lwd = 2)
    } else {
      lines(freq_hv_mean, apply(h1_mat[, idx_select], 1, function(x) exp(mean(log(x)))), col = 'red', lwd = 2)
      lines(freq_hv_mean, apply(h2_mat[, idx_select], 1, function(x) exp(mean(log(x)))), col = 'green', lwd = 2)
      lines(freq_hv_mean, apply(v_mat[, idx_select], 1, function(x) exp(mean(log(x)))), col = 'blue', lwd = 2)
    }
    legend('topright', legend = c('H1', 'H2', 'V'), col = c('red', 'green', 'blue'), lty = 1,
           lwd = 2)

    plot(x_range, c(y_range[1], y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
         ylab = 'HVSR', log = 'x', xaxt = 'n')
    log10_ticks(x = x_range, y = c(y_range[1], y_range[2]*1.5), log10_scale = 'x', tick_type = 'lin')
    abline(v = x_range[2], lwd = 2)
    if (pre_filter_flag) {
      if (!is.na(pre_filter_hpass_fc))
        abline(v = pre_filter_hpass_fc, lwd = 2, col = 'purple', lty = 2)
      if (!is.na(pre_filter_lpass_fc))
        abline(v = pre_filter_lpass_fc, lwd = 2, col = 'purple', lty = 2)
    }
    if (filter_flag) {
      if (!is.na(hpass_fc))
        abline(v = hpass_fc, lwd = 2, col = 'purple')
      if (!is.na(lpass_fc))
        abline(v = lpass_fc, lwd = 2, col = 'purple')
    }
    for(i_plot in idx_select){
      lines(freq_hv_mean, hvsr_mat[,i_plot])
    }

    if (robust_est) {
      hvsr_mean <- apply(hvsr_mat[, idx_select], 1, median)
      hvsr_sd <- apply(hvsr_mat[, idx_select], 1, mad) # /sqrt(ncol(hvsr_mat))
      hvsr_sd1 <- apply(hvsr_mat, 1, IQR)
      lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
      lines(freq_hv_mean, hvsr_mean - hvsr_sd, lwd = 2, col = 'blue', lty = 2)
      lines(freq_hv_mean, hvsr_mean + hvsr_sd, lwd = 2, col = 'blue', lty = 2)
      legend('top', 'Please select HVSR curves for removal.', bg = 'white')
    } else {
      if (distribution == 'normal') {
        hvsr_mean <- apply(hvsr_mat[, idx_select], 1, mean)
        hvsr_sd <- apply(hvsr_mat[, idx_select], 1, sd) # /sqrt(ncol(hvsr_mat))
        lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
        lines(freq_hv_mean, hvsr_mean - hvsr_sd, lwd = 2, col = 'blue', lty = 2)
        lines(freq_hv_mean, hvsr_mean + hvsr_sd, lwd = 2, col = 'blue', lty = 2)
        legend('top', 'Please select HVSR curves for removal.', bg = 'white')
      } else if (distribution == 'log_normal') {
        hvsr_mean <- exp(apply(log(hvsr_mat[, idx_select]), 1, mean))
        hvsr_sd <- apply(log(hvsr_mat[, idx_select]), 1, sd) # /sqrt(ncol(hvsr_mat))
        lines(freq_hv_mean, hvsr_mean, lwd = 3, col = 'red')
        lines(freq_hv_mean, hvsr_mean / exp(hvsr_sd), lwd = 2, col = 'blue', lty = 2)
        lines(freq_hv_mean, hvsr_mean * exp(hvsr_sd), lwd = 2, col = 'blue', lty = 2)
        legend('top', 'Please select HVSR curves for removal.', bg = 'white')
      }
    }
    p <- locator(1)
  }
  legend('center', 'HVSR curves selection is DONE !', box.col = "lightblue", bg = "lightblue")
  res <- list()
  res$idx_select <- idx_select
  res$hvsr_mean <- hvsr_mean
  res$hvsr_sd <- hvsr_sd
  if (robust_est) res$hvsr_sd1 <- hvsr_sd1
  return(res)
}
