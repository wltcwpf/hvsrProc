#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//' Taper function
//'
//' Apply Taper to the time series
//' @useDynLib hvsrProc
//' @importFrom Rcpp sourceCpp
//' @param ts An array of time series
//' @param t_front A number, the percentage of taperring for the beginning of the time series
//' @param t_end A number, the percentage of taperring for the end of the time series
//' @return The Taperred time series.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector taper(Rcpp::NumericVector ts, double t_front = 5, double t_end = 5) {

  int ts_len = ts.size();
  NumericVector ts_new(ts_len);
  for (int i = 0; i < ts_len; i++) {
    ts_new[i] = ts[i];
  }
  if (t_front > 0) {
    int num_front = floor(ts_len * t_front / 100);
    for (int i = 0; i < num_front; i++){
      ts_new[i] = ts[i] * (1 + cos(M_PI * i / num_front + M_PI)) / 2;
    }
  }
  if (t_end > 0) {
    int num_end = floor(ts_len * t_end / 100);
    for (int i = 0; i < num_end; i++){
      ts_new[ts_len - 1 - i] = ts[ts_len - 1 - i] * (1 + cos(M_PI * i / num_end + M_PI)) / 2;
    }
  }

  return ts_new;
}
