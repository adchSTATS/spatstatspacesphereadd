// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// engine_K
NumericMatrix engine_K(NumericVector r, NumericVector s, NumericMatrix dists_3d, NumericMatrix dists_sph, NumericMatrix Dmat);
RcppExport SEXP _spatstatspacesphereadd_engine_K(SEXP rSEXP, SEXP sSEXP, SEXP dists_3dSEXP, SEXP dists_sphSEXP, SEXP DmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists_3d(dists_3dSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists_sph(dists_sphSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dmat(DmatSEXP);
    rcpp_result_gen = Rcpp::wrap(engine_K(r, s, dists_3d, dists_sph, Dmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatstatspacesphereadd_engine_K", (DL_FUNC) &_spatstatspacesphereadd_engine_K, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatstatspacesphereadd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}