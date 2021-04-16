#' Generate linear log space
#'
#' This function generates n points with equal space in log scale
#' @param x1 The start point
#' @param x2 The end point
#' @param log_type The type of log scale. "log10" is for log10(), "ln" is for natural log()
#' @param n The number of points to be generated
#' @return An array of n points with equal space in log scale
#' @export
logspace <- function (x1, x2, log_type = "log10", n = 50) {
  n <- floor(n)
  if (log_type == 'log10') {
    return(10^seq(log10(x1), log10(x2), length.out = n))
  }
  if (log_type == 'ln') {
    return(exp(seq(log(x1), log(x2), length.out = n)))
  }
}
