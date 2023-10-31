#' Time-domain plot and window selection
#'
#' This function plots each window in time domain and let users to manually select windows for rejection
#' @param h1_wins A list. Each element is a window time series for 1st horizontal component
#' @param h2_wins A list. Each element is a window time series for 2nd horizontal component
#' @param v_wins A list. Each element is a window time series for vertical component
#' @param dt The time step
#' @param sta_lta_flag The flag indicates if STA/LTA is calculated
#' @param h1_stalta A list. Each element is a STA/LTA for a window of 1st horizontal component
#' @param h2_stalta A list. Each element is a STA/LTA for a window of 2nd horizontal component
#' @param v_stalta A list. Each element is a STA/LTA for a window of vertical component
#' @param visual_deci_flag The flag indicates if decimation is applied.
#' Note it will be very slow if the input data is large but without decimation
#' @param visual_deci A decimation factor to speed up plotting. The larger it is, the faster the plot is generated.
#' If visual_deci is NA, then the code automatically determines how much to decimate.
#' @return The lefted indices of windows after time-domain rejection
#' @importFrom graphics plot abline lines legend locator text
#' @importFrom grDevices dev.off
#' @export
td_plt_select <- function(h1_wins, h2_wins, v_wins, dt, sta_lta_flag, h1_stalta, h2_stalta, v_stalta, visual_deci_flag = TRUE,
                          visual_deci = NA) {

  num_wins <- length(h1_wins)
  cols <- rep(c('pink', 'lightgreen'), num_wins)
  len_wins <- unlist(lapply(h1_wins, length))
  tt_len <- sum(len_wins)
  if (tt_len <= 108000) {
    auto_visual_deci <- 1
  } else {
    auto_visual_deci <- floor(tt_len / 108000)  # 216000 is the max number of data points will be plotted
  }
  t_begin <- (c(0, cumsum(len_wins)[-num_wins]) + 1) * dt
  t_end <- cumsum(len_wins) * dt

  # shift h2_wins and v_wins
  h1_range <- range(unlist(lapply(h1_wins, range)))
  h2_range <- range(unlist(lapply(h2_wins, range)))
  shift_1 <- h2_range[2] - h1_range[1]
  h2_wins_sh <- lapply(h2_wins, function(x, shift_1) { x <- x - shift_1 }, shift_1 = shift_1)
  h2_range <- range(unlist(lapply(h2_wins_sh, range)))
  v_range <- range(unlist(lapply(v_wins, range)))
  shift_2 <- v_range[2] - h2_range[1]
  v_wins_sh <- lapply(v_wins, function(x, shift_2) { x <- x - shift_2 }, shift_2 = shift_2)
  x_range <- range(c(0, dt * tt_len))
  y_range <- range(c(range(unlist(lapply(h1_wins, range))), range(unlist(lapply(h2_wins_sh, range))),
                     range(unlist(lapply(v_wins_sh, range)))))

  # generate the plot
  par(mfrow = c(1,1))
  plot(x_range, y_range, type = 'n', xlab = 'Time (s)', ylab = 'Acc.', yaxt = 'n')
  abline(v = x_range[2], lwd = 2)
  for (i_plot in 1:num_wins) {
    t_seq <- seq(1, length(h1_wins[[ i_plot ]])) * dt + (i_plot - 1) * length(h1_wins[[ i_plot ]]) * dt
    ifelse (visual_deci_flag,
      ifelse(!is.na(visual_deci) & visual_deci > 0,
             idx <- seq(1, length(t_seq), by = visual_deci),
             idx <- seq(1, length(t_seq), by = auto_visual_deci)),
      idx <- seq(1, length(t_seq)))
    tmp_ts <- h1_wins[[ i_plot ]]
    lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])
    tmp_ts <- h2_wins_sh[[ i_plot ]]
    lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])
    tmp_ts <- v_wins_sh[[ i_plot ]]
    lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])

    # add STA/LTA
    if (sta_lta_flag) {
      range_h1 <- max(h1_wins[[ i_plot ]]) - min(h1_wins[[ i_plot ]])
      text(mean(t_seq[idx]), range(h1_wins[[ i_plot ]])[1] + range_h1*0.2, round(min(h1_stalta[[ i_plot ]]), digits = 1))
      text(mean(t_seq[idx]), range(h1_wins[[ i_plot ]])[2] - range_h1*0.2, round(max(h1_stalta[[ i_plot ]]), digits = 1))

      range_h2 <- max(h2_wins_sh[[ i_plot ]]) - min(h2_wins_sh[[ i_plot ]])
      text(mean(t_seq[idx]), range(h2_wins_sh[[ i_plot ]])[1] + range_h2*0.2, round(min(h2_stalta[[ i_plot ]]), digits = 1))
      text(mean(t_seq[idx]), range(h2_wins_sh[[ i_plot ]])[2] - range_h2*0.2, round(max(h2_stalta[[ i_plot ]]), digits = 1))

      range_v <- max(v_wins_sh[[ i_plot ]]) - min(v_wins_sh[[ i_plot ]])
      text(mean(t_seq[idx]), range(v_wins_sh[[ i_plot ]])[1] + range_v*0.2, round(min(v_stalta[[ i_plot ]]), digits = 1))
      text(mean(t_seq[idx]), range(v_wins_sh[[ i_plot ]])[2] - range_v*0.2, round(max(v_stalta[[ i_plot ]]), digits = 1))

      # text(mean(t_seq[idx]), 0, round(mean(h1_stalta[[ i_plot ]]), digits = 1))
      # text(mean(t_seq[idx]), -shift_1, round(mean(h2_stalta[[ i_plot ]]), digits = 1))
      # text(mean(t_seq[idx]), -shift_2, round(mean(v_stalta[[ i_plot ]]), digits = 1))
    }
  }
  legend('top', 'Please select time series for removal.', bg = 'white')

  # start selection process
  p <- locator(1)
  idx_select <- 1:num_wins
  idx_remove <- c()
  while (p$x <= x_range[2]) {
    while (p$x <= x_range[2]) {
      idx_temp <- which(p$x < t_end & p$x > t_begin)
      if(idx_temp %in% idx_select){
        idx_remove <- c(idx_remove, idx_temp)
        idx_select <- idx_select[-which(idx_temp == idx_select)]
      }
      t_seq <- seq(1, length(h1_wins[[ idx_temp ]])) * dt + (idx_temp - 1) * length(h1_wins[[ idx_temp ]]) * dt
      ifelse (visual_deci_flag,
              ifelse(!is.na(visual_deci) & visual_deci > 0,
                     idx <- seq(1, length(t_seq), by = visual_deci),
                     idx <- seq(1, length(t_seq), by = auto_visual_deci)),
              idx <- seq(1, length(t_seq)))
      tmp_ts <- h1_wins[[ idx_temp ]]
      lines(t_seq[idx], tmp_ts[idx], col = 'gray')
      tmp_ts <- h2_wins_sh[[ idx_temp ]]
      lines(t_seq[idx], tmp_ts[idx], col = 'gray')
      tmp_ts <- v_wins_sh[[ idx_temp ]]
      lines(t_seq[idx], tmp_ts[idx], col = 'gray')
      p <- locator(1)
    }

    iidx_select <- c(sapply(idx_select, function(x) 2 * (x - 1) + seq(1, 2)))
    h1_range <- range(unlist(lapply(h1_wins, range))[iidx_select])
    h2_range <- range(unlist(lapply(h2_wins, range))[iidx_select])
    v_range <- range(unlist(lapply(v_wins, range))[iidx_select])
    shift_1 <- h2_range[2] - h1_range[1]
    h2_wins_sh <- lapply(h2_wins, function(x, shift_1) { x <- x - shift_1 }, shift_1 = shift_1)
    h2_range_sh <- range(unlist(lapply(h2_wins_sh, range))[iidx_select])
    shift_2 <- v_range[2] - h2_range_sh[1]
    v_wins_sh <- lapply(v_wins, function(x, shift_2) { x <- x - shift_2 }, shift_2 = shift_2)
    x_range <- range(c(0, dt * tt_len))
    y_range <- range(c(range(unlist(lapply(h1_wins, range))[iidx_select]),
                       range(unlist(lapply(h2_wins_sh, range))[iidx_select]),
                       range(unlist(lapply(v_wins_sh, range))[iidx_select])))

    plot(x_range, y_range, type = 'n', xlab = 'Time (s)', ylab = 'Acc.', yaxt = 'n')
    abline(v = x_range[2], lwd = 2)
    legend('top', 'Please select time series for removal.', bg = 'white')
    for(i_plot in idx_select){
      t_seq <- seq(1, length(h1_wins[[ i_plot ]])) * dt + (i_plot - 1) * length(h1_wins[[ i_plot ]]) * dt
      ifelse (visual_deci_flag,
              ifelse(!is.na(visual_deci) & visual_deci > 0,
                     idx <- seq(1, length(t_seq), by = visual_deci),
                     idx <- seq(1, length(t_seq), by = auto_visual_deci)),
              idx <- seq(1, length(t_seq)))
      tmp_ts <- h1_wins[[ i_plot ]]
      lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])
      tmp_ts <- h2_wins_sh[[ i_plot ]]
      lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])
      tmp_ts <- v_wins_sh[[ i_plot ]]
      lines(t_seq[idx], tmp_ts[idx], col = cols[i_plot])

      # add STA/LTA
      if (sta_lta_flag) {
        range_h1 <- max(h1_wins[[ i_plot ]]) - min(h1_wins[[ i_plot ]])
        text(mean(t_seq[idx]), range(h1_wins[[ i_plot ]])[1] + range_h1*0.2, round(min(h1_stalta[[ i_plot ]]), digits = 1))
        text(mean(t_seq[idx]), range(h1_wins[[ i_plot ]])[2] - range_h1*0.2, round(max(h1_stalta[[ i_plot ]]), digits = 1))

        range_h2 <- max(h2_wins_sh[[ i_plot ]]) - min(h2_wins_sh[[ i_plot ]])
        text(mean(t_seq[idx]), range(h2_wins_sh[[ i_plot ]])[1] + range_h2*0.2, round(min(h2_stalta[[ i_plot ]]), digits = 1))
        text(mean(t_seq[idx]), range(h2_wins_sh[[ i_plot ]])[2] - range_h2*0.2, round(max(h2_stalta[[ i_plot ]]), digits = 1))

        range_v <- max(v_wins_sh[[ i_plot ]]) - min(v_wins_sh[[ i_plot ]])
        text(mean(t_seq[idx]), range(v_wins_sh[[ i_plot ]])[1] + range_v*0.2, round(min(v_stalta[[ i_plot ]]), digits = 1))
        text(mean(t_seq[idx]), range(v_wins_sh[[ i_plot ]])[2] - range_v*0.2, round(max(v_stalta[[ i_plot ]]), digits = 1))

        # text(mean(t_seq[idx]), 0, round(mean(h1_stalta[[ i_plot ]]), digits = 1))
        # text(mean(t_seq[idx]), -shift_1, round(mean(h2_stalta[[ i_plot ]]), digits = 1))
        # text(mean(t_seq[idx]), -shift_2, round(mean(v_stalta[[ i_plot ]]), digits = 1))
      }
    }
    p <- locator(1)
  }
  legend('center', 'Time series selection is DONE !', box.col = "lightblue", bg = "lightblue")
  idx_select
}
