// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mean_removal
Rcpp::NumericVector mean_removal(Rcpp::NumericVector ts, int nDC);
RcppExport SEXP _hvsrProc_mean_removal(SEXP tsSEXP, SEXP nDCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type nDC(nDCSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_removal(ts, nDC));
    return rcpp_result_gen;
END_RCPP
}
// ko_smooth
arma::vec ko_smooth(arma::vec freq, arma::vec amp, float b, float rate);
RcppExport SEXP _hvsrProc_ko_smooth(SEXP freqSEXP, SEXP ampSEXP, SEXP bSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type amp(ampSEXP);
    Rcpp::traits::input_parameter< float >::type b(bSEXP);
    Rcpp::traits::input_parameter< float >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(ko_smooth(freq, amp, b, rate));
    return rcpp_result_gen;
END_RCPP
}
// parzen_smooth
arma::vec parzen_smooth(arma::vec freq, arma::vec amp, float b);
RcppExport SEXP _hvsrProc_parzen_smooth(SEXP freqSEXP, SEXP ampSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type amp(ampSEXP);
    Rcpp::traits::input_parameter< float >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(parzen_smooth(freq, amp, b));
    return rcpp_result_gen;
END_RCPP
}
// sta_lta_calc
arma::vec sta_lta_calc(arma::vec ts, int short_term, int long_term, int moving_term);
RcppExport SEXP _hvsrProc_sta_lta_calc(SEXP tsSEXP, SEXP short_termSEXP, SEXP long_termSEXP, SEXP moving_termSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type short_term(short_termSEXP);
    Rcpp::traits::input_parameter< int >::type long_term(long_termSEXP);
    Rcpp::traits::input_parameter< int >::type moving_term(moving_termSEXP);
    rcpp_result_gen = Rcpp::wrap(sta_lta_calc(ts, short_term, long_term, moving_term));
    return rcpp_result_gen;
END_RCPP
}
// taper
Rcpp::NumericVector taper(Rcpp::NumericVector ts, double t_front, double t_end);
RcppExport SEXP _hvsrProc_taper(SEXP tsSEXP, SEXP t_frontSEXP, SEXP t_endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< double >::type t_front(t_frontSEXP);
    Rcpp::traits::input_parameter< double >::type t_end(t_endSEXP);
    rcpp_result_gen = Rcpp::wrap(taper(ts, t_front, t_end));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hvsrProc_mean_removal", (DL_FUNC) &_hvsrProc_mean_removal, 2},
    {"_hvsrProc_ko_smooth", (DL_FUNC) &_hvsrProc_ko_smooth, 4},
    {"_hvsrProc_parzen_smooth", (DL_FUNC) &_hvsrProc_parzen_smooth, 3},
    {"_hvsrProc_sta_lta_calc", (DL_FUNC) &_hvsrProc_sta_lta_calc, 4},
    {"_hvsrProc_taper", (DL_FUNC) &_hvsrProc_taper, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_hvsrProc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
