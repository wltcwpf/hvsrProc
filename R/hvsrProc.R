#' HVSR Processing Package
#'
#' This package is developed for processing and calculating
#' HVSR from noise-based and earthquake-based time series data and
#' parameterizing HVSR peaks (if exist) using a generalized Gaussian pulse function.
#'
#' @docType package
#'
#' @name hvsrProc-package
#' @aliases hvsrProc hvsrProc-package
#'
#' @note There are 20 hvsrProc functions:
#' \code{ang_pga_rotd50_calc, bw_pass, bw_resp, fas_cal, fd_plt_select,
#' hv_proc, hvsr_win_calc, ko_smooth, log10_ticks, logspace,
#' mean_removal, parzen_smooth, plot_log10, pre_proc,
#' sta_lta_calc, taper, td_plt_select, peak_fit_auto, peak_fit_manual,
#' gaussian_gen_fit}
#'
#' Of which \code{hv_proc} is main function for HVSR processing.
#'
#' \code{peak_fit_auto} is the function implements regression tree method to automatically
#' detect if a clear peak exists and parameterize the clear peak by a generalized Gaussian
#' pulse function.
#'
#' \code{peak_fit_manual} is the function takes user defined peak range to fit the peak
#' by a generalized Gaussian pulse function.
#'
#' @examples
#' ## Not run:
#' ##   See examples in the help files for all functions.
#' ## End(Not run)
#'
#' @author Pengfei Wang
NULL
#> NULL
