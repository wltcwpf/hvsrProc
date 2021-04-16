#include <Rmath.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' STA to LTA calculation
//'
//' Calculate the ratio of Short Term Average (STA) to Long Term Average (LTA)
//' @useDynLib hvsrProc
//' @importFrom Rcpp sourceCpp
//' @param ts An array of time series
//' @param short_term An integer, short term length for STA (number of points)
//' @param long_term An integer, long term length for LTA (number of points)
//' @param moving_term An integer, moving step of the preceding point (number of points)
//' @return A vector, the STA/LTA ratios
//' @export
// [[Rcpp::export]]
arma::vec sta_lta_calc(arma::vec ts, int short_term, int long_term, int moving_term) {

  int num_moving = floor((ts.size() - short_term) / moving_term) + 1;
  arma::vec sta_lta_ratio(num_moving);
  arma::vec short_temp, long_temp;
  int idx1;
  for (int i = 0; i < num_moving; i++) {
    short_temp = ts.subvec(i * moving_term, i * moving_term + short_term - 1);
    idx1 = (i * moving_term + short_term - 1) - (long_term - 1);
    if (idx1 < 0) {
      long_temp = ts;
    } else {
      long_temp = ts.subvec(idx1, i * moving_term + short_term - 1);
    }
    sta_lta_ratio[i] = mean(abs(short_temp)) / mean(abs(long_temp));
  }

  return sta_lta_ratio;
}
