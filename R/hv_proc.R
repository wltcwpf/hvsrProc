#' HVSR-processing main function
#'
#' This is the main function for HVSR processing (ambient-noise or corrected strong motions)
#' @param is_noise A binary indicating if the input data is ambient noise. is_noise = TRUE if input is noise;
#' otherwise, if is_noise = FALSE, the input is corrected strong motion. Note when is_noise = TRUE, h1,
#' h2, v, and, dt are needed; when is_noise = FALSE, eqk_filepath is needed.
#' @param h1 The first component of horizontal time series
#' @param h2 The second component of horizontal time series
#' @param v The vertical component of time series
#' @param dt The time step
#' @param eqk_filepath The file path where it saves the file of corrected strong motions data in standard format
#' @param output_dir The directory of output results
#' @param output_pf_flnm The prefix name of output file names
#' @param distribution The distribution assumption on the HVSR amplitudes. It can take "normal" or "log_normal".
#' We recommend normal for noise data and log_normal for corrected strong motions
#' @param filter_flag The flag indicates if filter is applied to noise. The corrected strong motions do not apply
#' @param is_causal A binary to indicator if causal is applied to noise. Apply causal if is_causal = TRUE and apply acausal if is_causal = FALSE
#' @param hpass_fc The corner frequency for high-pass filter
#' @param lpass_fc The corner frequency for low-pass filter
#' @param nPole_hp The pole parameter for high-pass. Default is 5
#' @param nPole_lp The pole parameter for low-pass. Default is 4
#' @param order_zero_padding The order needs to be added for zeroes padding at the end of recordings to
#' increase the number of data points to a power of 2. Default is 2
#' @param detrend The indicator specifies which detrend method is used. 0: no detrend; 1: mean removal; 2: linear trend removal; 6: fifth order polynomial detrend.
#' @param taper_flag The flag indicates if Taper is applied to noise. The corrected strong motions do not apply
#' @param t_front The percentage of taperring for the beginning of the time series.
#' @param t_end The percentage of taperring for the end of the time series.
#' @param horizontal_comb The parameter specifies the combination of two horizontal components. ps_RotD50: rotated combination at the angle where PGA is median; geometric_mean: geometric mean (sqrt(h1(f) * h2(f))); squared_average: squared average (sqrt((h1(f)^2 + h2(f)^2)/2))
#' @param ko_smooth_flag The flag indicates if KO smoothing is applied.
#' @param ko_smooth_b The coefficient of bandwidth. Default is 20. A smaller value will lead to more smoothing
#' @param parzen_flag The flag indicates if Parzen smoothing is applied. Default is FALSE since KO smoothing is applied by default
#' @param parzen_bwidth The bandwidth. A larger value will lead to more smoothing
#' @param win_width The window length for noise data. The unit is second. The win_width equals to one earthquake duration if
#' the input is corrected strong motion
#' @param overlapping Overlapping between two windows, unit is second. Default is 0 without any overlapping
#' @param sta_lta_flag The flag indicates if STA/LTA is calculated
#' @param sta_lta_moving_term_len The length of moving step for STA/LTA, in second
#' @param short_term_len The short term length for STA, in second
#' @param long_term_len The long term length for LTA, in second
#' @param polar_curves_flag The flag indicates if polar curves are calculated
#' @param deg_increment The degree increment for HVSR polar curves
#' @param resample_lin2log The flag indicates if resampling frequency from equal-spaced linear scale to equal-spaced log scale
#' @param deci_mean_factor An integer for decimation factor for HVSR mean curve
#' @param deci_polar_factor An integer for decimation factor for HVSR polar curves
#' @param output_freq_min The minimum output frequency of HVSR
#' @param output_freq_max The maximum output frequency of HVSR
#' @param output_selected_ts The flag indicates if output selected time seriese windows
#' @param output_removed_ts The flag indicates if output rejected time seriese windows
#' @param output_selected_hvsr The flag indicates if output selected HVSR curves
#' @param output_removed_hvsr The flag indicates if output rejected HVSR curves
#' @param output_mean_curve The flag indicates if output HVSR mean curve
#' @param output_polar_curves The flag indicates if output HVSR polar curves
#' @param output_metadata The flag indicates if output HVSR processing meta data
#' @return The results are saved and written in files in the specified output_dir
#' @importFrom utils read.csv write.csv
#' @export
hv_proc <- function(is_noise = TRUE, h1, h2, v, dt, eqk_filepath, output_dir, output_pf_flnm = 'Test_',
                    distribution = 'normal', filter_flag = TRUE, is_causal = FALSE,
                    hpass_fc = 0.1, lpass_fc = NA, nPole_hp = 5, nPole_lp = 4, order_zero_padding = 2, detrend = 1,
                    taper_flag = TRUE, t_front = 5, t_end = 5, horizontal_comb = 'geometric_mean', ko_smooth_flag = TRUE, ko_smooth_b = 20,
                    parzen_flag = FALSE, parzen_bwidth = 1.5, win_width = 150, overlapping = 0,
                    sta_lta_flag = TRUE, sta_lta_moving_term_len = 1, short_term_len = 1, long_term_len = 30,
                    polar_curves_flag = TRUE, deg_increment = 10, resample_lin2log = TRUE, deci_mean_factor = 10,
                    deci_polar_factor = 10, output_freq_min = 0.01, output_freq_max = 50,
                    output_selected_ts = FALSE, output_removed_ts = FALSE,
                    output_selected_hvsr = TRUE, output_removed_hvsr = FALSE,
                    output_mean_curve = TRUE, output_polar_curves = TRUE,
                    output_metadata = TRUE) {

  dir.create(output_dir, showWarnings = FALSE)

  if (!(detrend %in% c(0, 1, 2, 6)))
    stop('Please specify a correct detrend code!')

  # pre-process noise data
  if (is_noise) {
    npts_win <- win_width / dt
    npts_over <- overlapping / dt
    win_moving <- npts_win - npts_over  # the length to move for the next window
    num_wins <- floor((length(h1) - npts_win) / win_moving + 1)
    if (sta_lta_flag) {
      short_term <- floor(short_term_len / dt)
      long_term <- floor(long_term_len / dt)
      sta_lta_moving <- floor(sta_lta_moving_term_len / dt)
    }

    # split data into num_wins windows
    h1_wins <- rep(list(), length.out = num_wins)
    h2_wins <- rep(list(), length.out = num_wins)
    v_wins <- rep(list(), length.out = num_wins)
    for (i in 1:num_wins) {
      h1_wins[[ i ]] <- h1[((i - 1) * win_moving + 1):(i * win_moving)]
      h2_wins[[ i ]] <- h2[((i - 1) * win_moving + 1):(i * win_moving)]
      v_wins[[ i ]] <- v[((i - 1) * win_moving + 1):(i * win_moving)]
    }

    # pre-processing
    if (!is.na(hpass_fc) & filter_flag) {
      h1_wins <- lapply(h1_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                         t_end = t_end, filter_flag = filter_flag, fc = hpass_fc, nPole = -nPole_hp,
                         is_causal = is_causal, order_zero_padding = order_zero_padding)
      h2_wins <- lapply(h2_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                         t_end = t_end, filter_flag = filter_flag, fc = hpass_fc, nPole = -nPole_hp,
                         is_causal = is_causal, order_zero_padding = order_zero_padding)
      v_wins <- lapply(v_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                        t_end = t_end, filter_flag = filter_flag, fc = hpass_fc, nPole = -nPole_hp,
                        is_causal = is_causal, order_zero_padding = order_zero_padding)
    }

    if (!is.na(lpass_fc) & filter_flag) {
      h1_wins <- lapply(h1_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                        t_end = t_end, filter_flag = filter_flag, fc = lpass_fc, nPole = nPole_lp,
                        is_causal = is_causal, order_zero_padding = order_zero_padding)
      h2_wins <- lapply(h2_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                        t_end = t_end, filter_flag = filter_flag, fc = lpass_fc, nPole = nPole_lp,
                        is_causal = is_causal, order_zero_padding = order_zero_padding)
      v_wins <- lapply(v_wins, pre_proc, dt = dt, detrend = detrend, taper_flag = taper_flag, t_front = t_front,
                        t_end = t_end, filter_flag = filter_flag, fc = lpass_fc, nPole = nPole_lp,
                        is_causal = is_causal, order_zero_padding = order_zero_padding)
    }
    print("Pre-processing noise data is DONE!")
  }

  # assemble corrected earthquake strong motions
  if (!is_noise) {
    eqk_data <- read.csv(eqk_filepath)
    num_wins <- dim(eqk_data)[2] / 4
    h1_wins <- rep(list(), length.out = num_wins)
    h2_wins <- rep(list(), length.out = num_wins)
    v_wins <- rep(list(), length.out = num_wins)
    for (i in 1:num_wins) {
      h1_wins[[ i ]] <- subset(eqk_data[, 2 + (i - 1) * 4], !is.na(eqk_data[, 2 + (i - 1) * 4]))
      h2_wins[[ i ]] <- subset(eqk_data[, 3 + (i - 1) * 4], !is.na(eqk_data[, 3 + (i - 1) * 4]))
      v_wins[[ i ]] <- subset(eqk_data[, i * 4], !is.na(eqk_data[, i * 4]))
    }
    print("Assembling earthquake strong motion is DONE!")
  }

  # time-domain plot and selection
  if (sta_lta_flag) {
    h1_stalta <- lapply(h1_wins, sta_lta_calc, short_term = short_term, long_term = long_term, moving_term = sta_lta_moving)
    h2_stalta <- lapply(h2_wins, sta_lta_calc, short_term = short_term, long_term = long_term, moving_term = sta_lta_moving)
    v_stalta <- lapply(v_wins, sta_lta_calc, short_term = short_term, long_term = long_term, moving_term = sta_lta_moving)
  }
  idx_select <- td_plt_select(h1_wins = h1_wins, h2_wins = h2_wins, v_wins = v_wins, dt = dt, sta_lta_flag = sta_lta_flag,
                              h1_stalta = h1_stalta, h2_stalta = h2_stalta, v_stalta = v_stalta, visual_deci_flag = TRUE,
                              visual_deci = NA)
  print("Time-domain selection is DONE!")

  # time-domain data output
  if (length(idx_select) == 0) {
    stop('No window is selected, please try different data!')
  } else {
    outputflname_ts_sel <- paste0(output_pf_flnm, 'ts_sel.csv')
    max_win_len <- max(unlist(lapply(h1_wins, length)))
    min_win_len <- min(unlist(lapply(h1_wins, length)))
    if (output_selected_ts) {
      ts_sel_out <- matrix(data = NA, nrow = max_win_len, ncol = 4*length(idx_select))
      for (i in 1:length(idx_select)) {
        ts_sel_out[, 4 * (i - 1) + 1] <- seq(1, length(h1_wins[[ idx_select[i] ]])) * dt
        ts_sel_out[, 4 * (i - 1) + 2] <- h1_wins[[ idx_select[i] ]]
        ts_sel_out[, 4 * (i - 1) + 3] <- h2_wins[[ idx_select[i] ]]
        ts_sel_out[, 4 * (i - 1) + 4] <- v_wins[[ idx_select[i] ]]
      }
      colnames(ts_sel_out) <- rep(NA, ncol(ts_sel_out))
      colnames(ts_sel_out)[seq(1, ncol(ts_sel_out), by = 4)] <- 'Time_s'
      colnames(ts_sel_out)[seq(2, ncol(ts_sel_out), by = 4)] <- 'H1_acc'
      colnames(ts_sel_out)[seq(3, ncol(ts_sel_out), by = 4)] <- 'H2_acc'
      colnames(ts_sel_out)[seq(4, ncol(ts_sel_out), by = 4)] <- 'V_acc'
      write.csv(ts_sel_out, paste(output_dir, outputflname_ts_sel, sep = '/'), row.names = FALSE)
    }
    if (output_removed_ts) {
      idx_remove <- seq(1, num_wins)[-idx_select]
      if (length(idx_remove) > 0) {
        outputflname_ts_unsel <- paste0(output_pf_flnm, 'ts_unsel.csv')
        max_win_len <- max(unlist(lapply(h1_wins, length)))
        min_win_len <- min(unlist(lapply(h1_wins, length)))
        ts_unsel_out <- matrix(data = NA, nrow = max_win_len, ncol = 4*length(idx_remove))
        for (i in 1:length(idx_remove)) {
          ts_unsel_out[, 4 * (i - 1) + 1] <- seq(1, length(h1_wins[[ idx_remove[i] ]])) * dt
          ts_unsel_out[, 4 * (i - 1) + 2] <- h1_wins[[ idx_remove[i] ]]
          ts_unsel_out[, 4 * (i - 1) + 3] <- h2_wins[[ idx_remove[i] ]]
          ts_unsel_out[, 4 * (i - 1) + 4] <- v_wins[[ idx_remove[i] ]]
        }
        colnames(ts_unsel_out) <- rep(NA, ncol(ts_unsel_out))
        colnames(ts_unsel_out)[seq(1, ncol(ts_unsel_out), by = 4)] <- 'Time_s'
        colnames(ts_unsel_out)[seq(2, ncol(ts_unsel_out), by = 4)] <- 'H1_acc'
        colnames(ts_unsel_out)[seq(3, ncol(ts_unsel_out), by = 4)] <- 'H2_acc'
        colnames(ts_unsel_out)[seq(4, ncol(ts_unsel_out), by = 4)] <- 'V_acc'
        write.csv(ts_unsel_out, paste(output_dir, outputflname_ts_unsel, sep = '/'), row.names = FALSE)
      }
    }

    # frequency-domain plot and selection
    print("Preparing for frequency-domain, please wait...")
    nNyq_min <- min_win_len / 2 + 1
    df <- 1 / (min_win_len * dt)
    freq <- ((1:nNyq_min) - 1) * df
    if (freq[1] == 0)
      freq <- freq[-1]
    if (resample_lin2log)
      freq <- logspace(x1 = min(freq), x2 = max(freq), n = length(freq))
    ifelse (deci_mean_factor > 0,
            freq_hv_mean <- freq[seq(1, length(freq), by = floor(deci_mean_factor))],
            freq_hv_mean <- freq)
    hvsr_list <- lapply(idx_select, hvsr_win_calc, h1_wins = h1_wins, h2_wins = h2_wins,
                        v_wins = v_wins, dt = dt, horizontal_comb = horizontal_comb, freq_hv_mean = freq_hv_mean,
                        polar_curves_flag = FALSE)
    fd_select <- fd_plt_select(hvsr_list = hvsr_list, freq_hv_mean = freq_hv_mean, freq_min = output_freq_min,
                               freq_max = output_freq_max, hpass_fc = hpass_fc, lpass_fc = lpass_fc,
                               distribution = distribution)
    iidx_select <- fd_select$idx_select
    print("Frequency-domain selection is DONE!")

    # frequency-domain data output
    if (length(iidx_select) == 0) {
      stop('No window is selected, please try different data!')
    } else {
      if (output_selected_hvsr) {
        hvsr_sel_out <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(iidx_select) + 1)
        hvsr_sel_out[, 1] <- freq_hv_mean
        for (i in 1:length(iidx_select))
          hvsr_sel_out[, i + 1] <- hvsr_list[[ iidx_select[i] ]]$hv_ratio
        colnames(hvsr_sel_out) <- rep(NA, ncol(hvsr_sel_out))
        colnames(hvsr_sel_out)[1] <- 'Freq_Hz'
        colnames(hvsr_sel_out)[-1] <- paste0(rep('HVSR_', length(iidx_select)), seq(1, length(iidx_select)))
        outputflname_hvsr_sel <- paste0(output_pf_flnm, 'hvsr_sel.csv')
        if (!is.na(output_freq_min) & output_freq_min > min(hvsr_sel_out[, 1])) {
          idx_min <- which(hvsr_sel_out[, 1] >= output_freq_min)
          hvsr_sel_out <- hvsr_sel_out[idx_min, ]
        }
        if (!is.na(output_freq_max) & output_freq_max < max(hvsr_sel_out[, 1])) {
          idx_max <- which(hvsr_sel_out[, 1] <= output_freq_max)
          hvsr_sel_out <- hvsr_sel_out[idx_max, ]
        }
        write.csv(hvsr_sel_out, paste(output_dir, outputflname_hvsr_sel, sep = '/'), row.names = FALSE)
      }

      if (output_mean_curve) {
        outputflname_hvsr_mean <- paste0(output_pf_flnm, 'hvsr_mean.csv')
        hvsr_mean_out <- cbind(freq_hv_mean, fd_select$hvsr_mean, fd_select$hvsr_sd)
        colnames(hvsr_mean_out) <- c('freq_Hz', 'HVSR mean', 'HVSR sd')
        if (!is.na(output_freq_min) & output_freq_min > min(hvsr_mean_out[, 1])) {
          idx_min <- which(hvsr_mean_out[, 1] >= output_freq_min)
          hvsr_mean_out <- hvsr_mean_out[idx_min, ]
        }
        if (!is.na(output_freq_max) & output_freq_max < max(hvsr_mean_out[, 1])) {
          idx_max <- which(hvsr_mean_out[, 1] <= output_freq_max)
          hvsr_mean_out <- hvsr_mean_out[idx_max, ]
        }
        write.csv(hvsr_mean_out, paste(output_dir, outputflname_hvsr_mean, sep = '/'), row.names = FALSE)
      }

      if (output_removed_hvsr) {
        if (length(idx_select) > length(iidx_select)) {
          idx_remove <- seq(1, length(idx_select))[-iidx_select]
          outputflname_hvsr_unsel <- paste0(output_pf_flnm, 'hvsr_unsel.csv')
          hvsr_unsel_out <- matrix(data = NA, nrow = length(freq_hv_mean), ncol = length(idx_remove) + 1)
          hvsr_unsel_out[, 1] <- freq_hv_mean
          for (i in 1:length(idx_remove))
            hvsr_unsel_out[, i + 1] <- hvsr_list[[ idx_remove[i] ]]$hv_ratio
          colnames(hvsr_unsel_out) <- rep(NA, ncol(hvsr_unsel_out))
          colnames(hvsr_unsel_out)[1] <- 'Freq_Hz'
          colnames(hvsr_unsel_out)[-1] <- paste0(rep('HVSR_', length(idx_remove)), seq(1, length(idx_remove)))
          if (!is.na(output_freq_min) & output_freq_min > min(hvsr_unsel_out[, 1])) {
            idx_min <- which(hvsr_unsel_out[, 1] >= output_freq_min)
            hvsr_unsel_out <- hvsr_unsel_out[idx_min, ]
          }
          if (!is.na(output_freq_max) & output_freq_max < max(hvsr_unsel_out[, 1])) {
            idx_max <- which(hvsr_unsel_out[, 1] <= output_freq_max)
            hvsr_unsel_out <- hvsr_unsel_out[idx_max, ]
          }
          write.csv(hvsr_unsel_out, paste(output_dir, outputflname_hvsr_unsel, sep = '/'), row.names = FALSE)
        }
      }

      # generate polar curves and output
      if (output_polar_curves) {
        print("Calculating and generating polar curve data, please wait......")
        outputflname_hvsr_polar <- paste0(output_pf_flnm, 'hvsr_polar.csv')
        ifelse (deci_polar_factor > 0,
                freq_polar <- freq[seq(1, length(freq), by = floor(deci_polar_factor))],
                freq_polar <- freq)
        hvsr_list <- lapply(idx_select[iidx_select], hvsr_win_calc, h1_wins = h1_wins, h2_wins = h2_wins,
                            v_wins = v_wins, dt = dt, horizontal_comb = horizontal_comb, polar_curves_flag = TRUE,
                            freq_polar = freq_polar, deg_increment = deg_increment)
        polar_degs <- seq(0, 179, by = deg_increment)
        polar_hvsr_mat <- matrix(data = NA, nrow = length(freq_polar), ncol = length(polar_degs) * 3)
        colnames(polar_hvsr_mat) <- seq(1, ncol(polar_hvsr_mat))
        tmp_hvsr_mat <- matrix(data = NA, nrow = length(freq_polar), ncol = length(iidx_select))
        for (i in 1:length(polar_degs)) {
          for (j in 1:length(iidx_select))
            tmp_hvsr_mat[, j] <- hvsr_list[[ j ]]$polar_hv_ratio[, i]
          if (distribution == 'normal') {
            tmp_mean <- apply(tmp_hvsr_mat, 1, mean)
            tmp_sd <- apply(tmp_hvsr_mat, 1, sd)
          } else if (distribution == 'log_normal') {
            tmp_hvsr_mat[tmp_hvsr_mat <= 0] <- 10e-5  # avoid log(0 or negative)
            hvsr_mean <- exp(apply(log(tmp_hvsr_mat), 1, mean))
            hvsr_sd <- apply(log(tmp_hvsr_mat), 1, sd)
          }
          polar_hvsr_mat[, 3 * (i - 1) + 1] <- freq_polar
          polar_hvsr_mat[, 3 * (i - 1) + 2] <- tmp_mean
          polar_hvsr_mat[, 3 * i] <- tmp_sd
          colnames(polar_hvsr_mat)[3 * (i - 1) + seq(1, 3)] <-
            c(paste0('Freq_', polar_degs[i]), paste0('HVSR_', polar_degs[i]), paste0('HVSR_sd_', polar_degs[i]))
        }
        if (!is.na(output_freq_min) & output_freq_min > min(polar_hvsr_mat[, 1])) {
          idx_min <- which(polar_hvsr_mat[, 1] >= output_freq_min)
          polar_hvsr_mat <- polar_hvsr_mat[idx_min, ]
        }
        if (!is.na(output_freq_max) & output_freq_max < max(polar_hvsr_mat[, 1])) {
          idx_max <- which(polar_hvsr_mat[, 1] <= output_freq_max)
          polar_hvsr_mat <- polar_hvsr_mat[idx_max, ]
        }
        write.csv(polar_hvsr_mat, paste(output_dir, outputflname_hvsr_polar, sep = '/'), row.names = FALSE)
      }

      # output Metadata
      if (output_metadata) {
        outputflname_meta <- paste0(output_pf_flnm, 'metadata.csv')
        Meta_names <- c('sample freq (Hz)', 'record duration (min)', 'detrend',
                        'time window (sec)', 'window overlap (sec)', 'taper type', 'front taper width (percentage)',
                        'end taper width (percentage)', 'horizontal combination', 'number of windows (total)',
                        'number of windows (selected)', 'high pass filter',
                        'high pass filter corner frequency (Hz)', 'high pass filter type',
                        'smoothing type', 'smoothing constant', 'data type', 'distribution')
        Meta_output <- data.frame(matrix(data = NA, nrow = 1, ncol = length(Meta_names))) # sample freq (Hz), record duration (min), mean removal,
        colnames(Meta_output) <- Meta_names
        Meta_output[1] <- 1/dt   # Hz
        Meta_output[2] <- length(h1) * dt /60   # total duration
        Meta_output[3] <- ifelse(detrend == 0, 'no detrend', ifelse(detrend == 1, 'mean removal', ifelse(detrend == 2, 'linear detrend', 'fifth order polynomial detrend')))  # detrend
        Meta_output[5] <- overlapping
        Meta_output[6] <- 'Tukey'
        Meta_output[7] <- t_front
        Meta_output[8] <- t_end
        Meta_output[9] <- horizontal_comb
        Meta_output[10] <- num_wins
        Meta_output[11] <- length(iidx_select)
        Meta_output[12] <- ifelse(filter_flag & !is.na(hpass_fc), 0, 1)
        Meta_output[13] <- hpass_fc
        Meta_output[14] <- 'Butterworth'
        Meta_output[15] <- 'KonnoOhmachi'
        Meta_output[16] <- ko_smooth_b
        Meta_output[17] <- ifelse(is_noise, 0, 1) # 0: noise; 1: earthquake strong motions
        Meta_output[18] <- distribution
        write.csv(Meta_output, paste(output_dir, outputflname_meta, sep = '/'), row.names = FALSE)
      }
    }
  }
  print("Everything is DONE, check out the results in output folder!")
}
