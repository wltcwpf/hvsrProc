#' A function for HVSR peak detection and fitting automatically
#'
#' This function automatically detect if a clear peak exits. If there is a clear peak,
#' a Gaussian pulse function is used to fit the peak. The fitted parameters are returned.
#' @param x_freq The array of frequnecy (in Hz)
#' @param hvsr_mean The array of HVSR mean curve
#' @param hvsr_sd The array of HVSR amplitude standard deviation
#' @param ratio_thres The ratio threshold to determine if a peak is "clear" enough
#' @param amp_thres The minimum required peak amplitude for a clear peak
#' @param k_cv_thres The threshold of coefficient of variation for a clear peak
#' @param max_freq The maximum frequency to be considered as a peak
#' @param min_freq The minimum frequency to be considered as a peak
#' @param cp The complexity parameter for regression tree. A smaller number leads to a finer fitting.
#' @param x_jump The step width to determine a stable plateau
#' @param peak_width_thres_L The threshold of peak width (in natural log of frequency) for a clear peak
#' @importFrom rpart rpart
#' @importFrom stats approx optim predict
#' @export
peak_fit_auto <- function( x_freq, hvsr_mean, hvsr_sd, ratio_thres = 0.7, amp_thres = 1.5,
                           k_cv_thres = 1, max_freq = 15, min_freq = 0.1, cp = 0.005,
                           x_jump = 0.35, peak_width_thres_L = 2.3 ){

  # tree regression to get the step functions
  res <- tree_reg( x_freq = x_freq, y_hvsr = hvsr_mean, y_hvsrsd = hvsr_sd, ratio_thres = ratio_thres,
                   amp_thres = amp_thres, k_cv_thres = k_cv_thres, cp = cp, x_jump = x_jump,
                   peak_width_thres_L = peak_width_thres_L )

  if( length( res[[ 1 ]] ) > 0 ){
    # there exists at least one significant peak

    # we take the first (with lowest peak freq.) peak
    res[[ 1 ]] <- res[[ 1 ]][[ 1 ]]

    # Gaussian fitting function
    # ----------------------------------------------------------------------------------------------------------
    res_Gauss <- res[[ 1 ]]

    c0 <- min( res[[ 4 ]][ res_Gauss$idx_left:res_Gauss$idx_right ] )

    # resample in linearly-natural log scale
    freq_log <- seq( min( log( x_freq ) ), max( log( x_freq ) ), length.out = length( x_freq ) )

    hvsr_mean <- approx( x = log( x_freq ), y = hvsr_mean, xout = freq_log )$y

    hvsr_sd <- approx( x = log( x_freq ), y = hvsr_sd, xout = freq_log )$y

    f1 <- res[[ 3 ]][ res_Gauss$idx_left ]

    f2 <- res[[ 3 ]][ res_Gauss$idx_right ]

    f1_pk <- res[[ 3 ]][ res_Gauss$peak_left ]

    f2_pk <- res[[ 3 ]][ res_Gauss$peak_right ]

    f_middle <- res_Gauss$peak_freq

    # low bound of parameters: c0, c1, mu, sigma
    lb <- c( c0, 0.1, exp(f1), 0.01 )

    # upper bound of parameters
    ub <- c( 10, 10, exp(f2), 0.5 )

    # find the indices of selected range
    idx1 <- which.min( abs( freq_log - f1 ) )

    idx2 <- which.min( abs( freq_log - f2 ) )

    # find the indices of peak range
    idx1_pk <- which.min( abs( freq_log - f1_pk ) )

    idx2_pk <- which.min( abs( freq_log - f2_pk ) )

    # get freq within the selected range
    freq_temp <- exp(freq_log[ idx1:idx2 ])

    # get hvsr within the selected range
    hv_temp <- hvsr_mean[ idx1:idx2 ]

    hv_sd_temp <- mean( hvsr_sd[ idx1_pk:idx2_pk ] )

    Res_ratio <- max( c( res_Gauss$left_ratio, res_Gauss$right_ratio ) )

    Res_stepAmp <- res_Gauss$peak_amp

    Res_cv <- hv_sd_temp / Res_stepAmp

    par_temp <- c( 0, 1, f_middle, 0.1 )

    par_temp <- optim( par = par_temp, fn = wavelet_res, x = freq_temp, y = hv_temp, method = 'L-BFGS-B',
                       lower = lb, upper = ub )

    newGauss_pars <- par_temp$par

    Res_mu <- newGauss_pars[ 3 ]

    if( Res_mu > max_freq | Res_mu < min_freq ){

      Res_c0 <- NA

      Res_c1 <- NA

      Res_sigma <- NA

      Res_Amp <- NA

      Peak_presence <- 'N'

    }else{

      Amp <- max( gaussian_gen_fit( newGauss_pars, exp(freq_log[ idx1:idx2 ] )) )

      Res_c0 <- newGauss_pars[ 1 ]

      Res_c1 <- newGauss_pars[ 2 ]

      Res_sigma <- newGauss_pars[ 4 ]

      Res_Amp <- Amp

      Peak_presence <- 'Y'
    }

  }else{

    Res_c0 <- NA

    Res_c1 <- NA

    Res_mu <- NA

    Res_sigma <- NA

    Res_Amp <- NA

    Res_ratio <- NA

    Res_stepAmp <- NA

    Res_cv <- NA

    Peak_presence <- 'N'
  }

  res <- data.frame( row.names = 1:length( Res_c0 ) )

  res$Peak_presence <- Peak_presence

  res$c0 <- Res_c0

  res$c1 <- Res_c1

  res$fp <- Res_mu

  res$w <- Res_sigma

  return( res )

}


#' A function for HVSR peak detection and fitting manually
#'
#' This function takes the user defined frequency range (lower bound and upper bound) and
#' use the Gaussian pulse function to fit the peak. The fitted parameters are returned.
#' @param x_freq The array of frequnecy (in Hz)
#' @param hvsr_mean The array of HVSR mean curve
#' @param fit_range_lb The lower frequency bound for a HVSR peak (in Hz)
#' @param fit_range_ub The upper frequency bound for a HVSR peak (in Hz)
#' @importFrom stats approx optim predict
#' @export
peak_fit_manual <- function( x_freq, hvsr_mean, fit_range_lb, fit_range_ub ){

  # low bound of parameters: c0, c1, mu, sigma
  lb <- c( 0.1, 0.1, exp(fit_range_lb), 0.01 )

  # upper bound of parameters
  ub <- c( 10, 10, exp(fit_range_ub), 0.5 )

  par_temp <- c( 0, 1, mean(c(fit_range_lb, fit_range_ub)), 0.1 )

  idx <- which(x_freq >= fit_range_lb & x_freq <= fit_range_ub)

  freq_temp <- x_freq[idx]

  hv_temp <- hvsr_mean[idx]

  par_temp <- optim( par = par_temp, fn = wavelet_res, x = freq_temp, y = hv_temp, method = 'L-BFGS-B',
                     lower = lb, upper = ub )

  res <- data.frame( row.names = 1 )

  res$c0 <- par_temp$par[1]

  res$c1 <- par_temp$par[2]

  res$fp <- par_temp$par[3]

  res$w <- par_temp$par[4]

  return( res )

}


# helper functions
# ------------------------------------------------------------------------------------------------------------
tree_reg <- function( x_freq, y_hvsr, y_hvsrsd = NULL, log.freq = FALSE, make.plot = FALSE, ratio_thres = 0.7,
                      amp_thres = 1.5, k_cv_thres = 0.5, cp = 0.005, x_jump = 0.35,
                      peak_width_thres_L = 2.3 ){
  # input: x_freq - freq, either in log or arithmetic scale
  #        y_hvsr - hvsr mean curve amplitudes
  #        y_hvsrsd - hvsr sd, could be NULL
  #        log.freq - indicator if x_freq is in log equal-spaced or arithmetic scale equal-spaced.
  #                   TRUE - log; FALSE - arithmetic
  #        make.plot - indicate if make a plot
  #        ratio_thres - the threshold of left_ratio and right_ratio to be considered as a peak
  #        amp_thres - the threshold of amplitude to be considered as a peak
  #        k_cv_thres - the threshold of coef. of variation (only work when y_hvsrsd has values) ->
  #                     lower bound of peak to be considered as a peak
  #        cp - the penalty parameter
  #        x_jump - the step to jump to next step for the left_min and right_min
  # output: a list contains:
  #        1 - the parameters for step functions
  #        2 - tree predictions
  #        3 - frequency in natural log scale
  #        4 - the corresponding hvsr at freq given in list 3

  if( !log.freq ){

    x_freqlog <- seq( min( log( x_freq ) ), max( log( x_freq ) ), length.out = length( x_freq ) )

    y_hvsrlog <- approx( x = log( x_freq ), y = y_hvsr, xout = x_freqlog )$y

    if( !(length( y_hvsrsd ) == 0) )
      y_hvsrsdlog <- approx( x = log( x_freq ), y = y_hvsrsd, xout = x_freqlog )$y

  }else{

    x_freqlog <- x_freq

    y_hvsrlog <- y_hvsr

    y_hvsrsdlog <- y_hvsrsd
  }

  df <- data.frame( x = x_freqlog, y = y_hvsrlog )

  tree <- rpart( y ~ x, data = df, control = list( cp = cp ) )

  tree_pred <- plot_tree( tree = tree, x = x_freqlog, y = y_hvsrlog, make.plot = make.plot )

  tree_pred_same <- plot_tree( tree = tree, x = x_freqlog, y = y_hvsrlog, dx = NA, make.plot = make.plot )

  y_pred <- tree_pred$y_pred

  x_inv <- sort( tree$splits[ , 4 ] )

  x_inv <- c( min( x_freqlog ), x_inv, max( x_freqlog ) )

  x_range <- cut( x_freqlog, breaks = x_inv )

  y_inv <- unique( y_pred )

  df <- data.frame( x_range = levels( x_range ), x_width = x_inv[ -1 ] - x_inv[ -length( x_inv ) ],
                    y_val = y_inv, x_middle = ( x_inv[ -1 ] + x_inv[ -length( x_inv ) ] ) / 2 )



  # loop each interval to check the potential significant peak
  num_pk <- 0

  res <- list()

  for( i in 1:nrow( df ) ){

    Tmp <- range_ratio( df = df, idx = i, x_jump = x_jump )

    Res_ratio <- mean( c( Tmp$left_ratio, Tmp$right_ratio ) )

    Res_ratio_max <- max( c( Tmp$left_ratio, Tmp$right_ratio ) )

    Res_stepAmp <- Tmp$peak_amp

    fit_range <- range( which( x_range %in% df$x_range[ Tmp$k_left:Tmp$k_right ] ) )

    peak_range <- range( which( x_range %in% df$x_range[ Tmp$idx ] ) )

    Tmp$idx_left <- fit_range[ 1 ]

    Tmp$idx_right <- fit_range[ 2 ]

    Tmp$peak_left <- peak_range[ 1 ]

    Tmp$peak_right <- peak_range[ 2 ]

    peak_width_L <- x_freqlog[ max( c( min( which( tree_pred_same$y_pred == Tmp$right_min ) ),
                                       Tmp$peak_right ) ) ] -
      x_freqlog[ min( c( max( which( tree_pred_same$y_pred == Tmp$left_min ) ),
                         Tmp$peak_left ) ) ]

    # mean of stds in the peak step range
    hv_sd_temp <- mean( y_hvsrsdlog[ Tmp$peak_left:Tmp$peak_right ] )

    # cv based on steps
    Res_cv <- hv_sd_temp / Res_stepAmp

    # the nearby low plateau
    Res_plateau <- max( c( Tmp$left_min, Tmp$right_min ) )

    # the lower bound of the peak step
    Res_lb <- Res_stepAmp * ( 1 - k_cv_thres * Res_cv )

    if( Res_ratio < ratio_thres & Res_stepAmp > amp_thres & Res_lb > Res_plateau &
        Res_ratio_max < 1 & peak_width_L < peak_width_thres_L ){
      # there is a significant peak
      num_pk <- num_pk + 1

      res[[ num_pk ]] <- Tmp
    }
  }

  return( list( res, tree_pred, x_freqlog, y_hvsrlog ) )
}


plot_tree <- function( tree, x, y, dx = 0.001, make.plot ){

  ifelse( !is.na( dx ), s <- seq( min( x ), max( x ), by = dx ),
          s <- x )

  y_pred <- predict( tree, data.frame( x = s ) )

  if( make.plot ){

    plot( x, y, xlab = 'Log of freq.', ylab = 'HVSR' )

    lines( s, y_pred, col = 'red', lwd = 2 )

  }

  res <- data.frame( s = s, y_pred = y_pred )

  return( res )
}

range_ratio <- function( df, idx, x_jump ){

  k_left <- idx

  k_right <- idx

  if( idx > 1 ){

    if( df$y_val[ idx - 1 ] < df$y_val[ idx ] ){

      k_left <- idx - 1

      while( ( k_left > 1 ) & ( df$x_width[ k_left ] < x_jump ) )
        ifelse( df$y_val[ k_left - 1 ] < df$y_val[ k_left ],
                k_left <- k_left - 1 ,
                break )
    }
  }

  if( idx < nrow( df ) ){

    if( df$y_val[ idx + 1 ] < df$y_val[ idx ] ){

      k_right <- idx + 1

      while( ( k_right < nrow( df ) ) & ( df$x_width[ k_right ] < x_jump ) )
        ifelse( df$y_val[ k_right + 1 ] < df$y_val[ k_right ],
                k_right <- k_right + 1,
                break )
    }
  }

  left_ratio <- df$y_val[ k_left ] / df$y_val[ idx ]

  right_ratio <- df$y_val[ k_right ] / df$y_val[ idx ]

  res <- data.frame( idx = idx, k_left = k_left, k_right = k_right, left_min = df$y_val[ k_left ],
                     right_min = df$y_val[ k_right ], peak_amp =df$y_val[ idx ],
                     peak_freq = exp( df$x_middle[ idx ] ),
                     left_ratio = left_ratio, right_ratio = right_ratio )

  return( res )
}

wavelet_res <- function( par, x, y ){

  ## Note: x - freq; y - arithmetic, linear
  y_hat <- gaussian_gen_fit( par, x )

  resd <- sum(( y - y_hat )^2)

  return(resd)
}


#' Generalized Gaussian pulse function
#'
#' This is a generalized Gaussian pulse function used to fit HVSR peaks
#' @param par The array of fitted Gaussian pulse function parameters, c0, c1, fp, and w
#' @param x_freq An array of frequency where we predict HVSR
#' @export
gaussian_gen_fit <- function( par, x_freq ){

  c0 <- par[ 1 ]  # constant term

  c1 <- par[ 2 ]  # multipler term

  mu <- par[ 3 ]  # peak x-coord

  sigma <- par[ 4 ]  # width

  y_hat <- c0 + c1 * exp(-(log(x_freq) - log(mu))^2 / 2 / sigma^2)

  return(y_hat)
}
