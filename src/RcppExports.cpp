// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rmultinom
arma::imat rmultinom(int n_boot, int num, int length, arma::vec weight);
RcppExport SEXP _safc_rmultinom(SEXP n_bootSEXP, SEXP numSEXP, SEXP lengthSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinom(n_boot, num, length, weight));
    return rcpp_result_gen;
END_RCPP
}
// get_demboot
arma::imat get_demboot(arma::vec agg_index, arma::vec num, arma::vec weight, int n_boot, int n_dest, int n_region);
RcppExport SEXP _safc_get_demboot(SEXP agg_indexSEXP, SEXP numSEXP, SEXP weightSEXP, SEXP n_bootSEXP, SEXP n_destSEXP, SEXP n_regionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type agg_index(agg_indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type num(numSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type n_boot(n_bootSEXP);
    Rcpp::traits::input_parameter< int >::type n_dest(n_destSEXP);
    Rcpp::traits::input_parameter< int >::type n_region(n_regionSEXP);
    rcpp_result_gen = Rcpp::wrap(get_demboot(agg_index, num, weight, n_boot, n_dest, n_region));
    return rcpp_result_gen;
END_RCPP
}
// dist
arma::vec dist(arma::mat origin, arma::mat destinations);
RcppExport SEXP _safc_dist(SEXP originSEXP, SEXP destinationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type origin(originSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type destinations(destinationsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist(origin, destinations));
    return rcpp_result_gen;
END_RCPP
}
// get_distmat
arma::field<arma::uvec> get_distmat(arma::mat origins, arma::mat destinations, double max_dist, int n_origin);
RcppExport SEXP _safc_get_distmat(SEXP originsSEXP, SEXP destinationsSEXP, SEXP max_distSEXP, SEXP n_originSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type origins(originsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type destinations(destinationsSEXP);
    Rcpp::traits::input_parameter< double >::type max_dist(max_distSEXP);
    Rcpp::traits::input_parameter< int >::type n_origin(n_originSEXP);
    rcpp_result_gen = Rcpp::wrap(get_distmat(origins, destinations, max_dist, n_origin));
    return rcpp_result_gen;
END_RCPP
}
// calc_avail
arma::vec calc_avail(arma::vec supply, arma::mat demand, arma::field<arma::uvec> distmat, arma::vec agg_index, int n_region, int n_origin);
RcppExport SEXP _safc_calc_avail(SEXP supplySEXP, SEXP demandSEXP, SEXP distmatSEXP, SEXP agg_indexSEXP, SEXP n_regionSEXP, SEXP n_originSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::uvec> >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type agg_index(agg_indexSEXP);
    Rcpp::traits::input_parameter< int >::type n_region(n_regionSEXP);
    Rcpp::traits::input_parameter< int >::type n_origin(n_originSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_avail(supply, demand, distmat, agg_index, n_region, n_origin));
    return rcpp_result_gen;
END_RCPP
}
// get_availability
arma::vec get_availability(arma::vec supply, arma::mat demand, arma::field<arma::uvec> distmat, arma::vec agg_index, int n_region, int n_origin);
RcppExport SEXP _safc_get_availability(SEXP supplySEXP, SEXP demandSEXP, SEXP distmatSEXP, SEXP agg_indexSEXP, SEXP n_regionSEXP, SEXP n_originSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::uvec> >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type agg_index(agg_indexSEXP);
    Rcpp::traits::input_parameter< int >::type n_region(n_regionSEXP);
    Rcpp::traits::input_parameter< int >::type n_origin(n_originSEXP);
    rcpp_result_gen = Rcpp::wrap(get_availability(supply, demand, distmat, agg_index, n_region, n_origin));
    return rcpp_result_gen;
END_RCPP
}
// progress
void progress(int s, int n_iter);
RcppExport SEXP _safc_progress(SEXP sSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    progress(s, n_iter);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_safc_rmultinom", (DL_FUNC) &_safc_rmultinom, 4},
    {"_safc_get_demboot", (DL_FUNC) &_safc_get_demboot, 6},
    {"_safc_dist", (DL_FUNC) &_safc_dist, 2},
    {"_safc_get_distmat", (DL_FUNC) &_safc_get_distmat, 4},
    {"_safc_calc_avail", (DL_FUNC) &_safc_calc_avail, 6},
    {"_safc_get_availability", (DL_FUNC) &_safc_get_availability, 6},
    {"_safc_progress", (DL_FUNC) &_safc_progress, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_safc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
