// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// engine_k3s
NumericMatrix engine_k3s(NumericVector r, NumericVector s, NumericMatrix dists_3d, NumericMatrix dists_sph, NumericMatrix Dmat);
RcppExport SEXP _spatstatspacesphereadd_engine_k3s(SEXP rSEXP, SEXP sSEXP, SEXP dists_3dSEXP, SEXP dists_sphSEXP, SEXP DmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists_3d(dists_3dSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists_sph(dists_sphSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dmat(DmatSEXP);
    rcpp_result_gen = Rcpp::wrap(engine_k3s(r, s, dists_3d, dists_sph, Dmat));
    return rcpp_result_gen;
END_RCPP
}
// engine_k3s2
NumericMatrix engine_k3s2(NumericVector r, NumericVector s, NumericVector x_vec, NumericVector y_vec, NumericVector z_vec, NumericMatrix dists_sph, NumericMatrix Dmat);
RcppExport SEXP _spatstatspacesphereadd_engine_k3s2(SEXP rSEXP, SEXP sSEXP, SEXP x_vecSEXP, SEXP y_vecSEXP, SEXP z_vecSEXP, SEXP dists_sphSEXP, SEXP DmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_vec(x_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_vec(y_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_vec(z_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists_sph(dists_sphSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dmat(DmatSEXP);
    rcpp_result_gen = Rcpp::wrap(engine_k3s2(r, s, x_vec, y_vec, z_vec, dists_sph, Dmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatstatspacesphereadd_engine_k3s", (DL_FUNC) &_spatstatspacesphereadd_engine_k3s, 5},
    {"_spatstatspacesphereadd_engine_k3s2", (DL_FUNC) &_spatstatspacesphereadd_engine_k3s2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatstatspacesphereadd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
