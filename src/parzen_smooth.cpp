#include <Rmath.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Parzen smoothing
//'
//' Apply Parzen smoothing to the Fourier Amplitude Spectra (FAS). The key formula is referred from Eq (5) on the page: http://www2.kobe-u.ac.jp/~nagaotak/poster2017wcee.pdf
//' @useDynLib hvsrProc
//' @importFrom Rcpp sourceCpp
//' @param freq An array of the frequency
//' @param amp An array of the corresponding FAS
//' @param b The bandwidth. A larger value will lead to more smoothing
//' @return The smoothed FAS
//' @export
// [[Rcpp::export]]
arma::vec parzen_smooth(arma::vec freq, arma::vec amp, float b=1.5) {

  int len_f = freq.size();
  double u = 151 * b / 280;
  double fc;
  arma::vec smoothed(len_f);
  arma::vec temp(len_f);
  arma::vec wb(len_f);
  arma::vec num(len_f);
  for(int i = 0; i < len_f; i++){
    fc = freq(i);
    temp = M_PI * u * (freq - fc) / 2;
    wb = pow(sin(temp) / temp, 4) * 3 / 4 * u;
    wb(i) = 1;
    num = wb % amp;
    smoothed(i) = sum(num) / sum(wb);
  }

  return smoothed;
}
