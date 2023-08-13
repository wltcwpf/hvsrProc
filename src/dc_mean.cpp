#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//' Mean removal
//'
//' Remove the mean
//' @useDynLib hvsrProc
//' @importFrom Rcpp sourceCpp
//' @param ts An array of time series
//' @param nDC a flag number, if negative, DC is based on the whole windowed record;
//' if positive, DC is based on the first nDC points.
//' @return The DC time series.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mean_removal(Rcpp::NumericVector ts, int nDC = -1) {

  double mean = 0;
  int ts_len = ts.size(), mean_len = ts.size();
  NumericVector ts_new(ts_len);
  if (nDC > 0) {
    mean_len = nDC;
  }
  for (int i = 0; i < mean_len; i++) {
    mean += ts[i];
  }
  mean = mean / mean_len;

  for (int i = 0; i < ts_len; i++) {
    ts_new[i] = ts[i] - mean;
  }

  return ts_new;
}
