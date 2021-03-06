#' Frequency-domain plot and window selection
#'
#' This function plots each window in frequency domain and let users to manually select windows for rejection
#' @param hvsr_list A list. Each element is a window HVSR
#' @param freq_hv_mean The frequency
#' @param freq_min The minimum frequency of interest
#' @param freq_max The maximum frequency of interest
#' @param hpass_fc High-pass corner frequency applied in filter
#' @param lpass_fc Low-pass corner frequency applied in filter
#' @param distribution The distribution assumption on the HVSR amplitudes. It can take "normal" or "log_normal".
#' @return The lefted indices of windows after frequency-domain rejection
#' @importFrom graphics plot abline lines legend locator text
#' @importFrom grDevices dev.off
#' @importFrom stats sd
#' @export
fd_plt_select <- function(hvsr_list, freq_hv_mean, freq_min, freq_max, hpass_fc, lpass_fc, distribution) {

  visual_freq_min <- max(c(freq_min, min(freq_hv_mean))) / 2
  x_range <- c(visual_freq_min, min(c(freq_max, max(freq_hv_mean))))
  hvsr_mat <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(hvsr_list))
  for (i in 1:length(hvsr_list)) {
    hvsr_mat[, i] <- hvsr_list[[ i ]]$rotd50_hv_ratio
  }
  y_range <- range(hvsr_mat[freq_hv_mean >= visual_freq_min, ])
  plot(x_range, c(y_range[1], y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
       ylab = 'HVSR', log = 'x', xaxt = 'n')
  log10_ticks(x = x_range, y = c(y_range[1], y_range[2]*1.5), log10_scale = 'x', tick_type = 'lin')
  abline(v = x_range[2], lwd = 2)
  if (!is.na(hpass_fc))
    abline(v = hpass_fc, lwd = 2, col = 'purple')
  if (!is.na(lpass_fc))
    abline(v = lpass_fc, lwd = 2, col = 'purple')
  for(i_plot in 1:ncol(hvsr_mat)){
    lines(freq_hv_mean, hvsr_mat[,i_plot])
  }
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
    plot(x_range, c(y_range[1], y_range[2]*1.5), type = 'n', xlab = 'Freq. (Hz)',
         ylab = 'HVSR', log = 'x', xaxt = 'n')
    log10_ticks(x = x_range, y = c(y_range[1], y_range[2]*1.5), log10_scale = 'x', tick_type = 'lin')
    abline(v = x_range[2], lwd = 2)
    if (!is.na(hpass_fc))
      abline(v = hpass_fc, lwd = 2, col = 'purple')
    if (!is.na(lpass_fc))
      abline(v = lpass_fc, lwd = 2, col = 'purple')
    for(i_plot in idx_select){
      lines(freq_hv_mean, hvsr_mat[,i_plot])
    }
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
    p <- locator(1)
  }
  legend('center', 'HVSR curves selection is DONE !', box.col = "lightblue", bg = "lightblue")
  res <- list()
  res$idx_select <- idx_select
  res$hvsr_mean <- hvsr_mean
  res$hvsr_sd <- hvsr_sd
  return(res)
}
