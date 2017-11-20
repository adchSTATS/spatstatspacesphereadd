#include <Rcpp.h>
using namespace Rcpp;

//' Engine for computing K-function
//'
//' Engine for calculating K-function of point processes on the product space of R^3 and S^2.
//' @param r Vector of values for the argument r at which K(r, s) should be evaluated.
//' @param s Vector of values for the argument s at which K(r, s) should be evaluated.
//' @param dists_3d Matrix of distances between points in R^3.
//' @param dists_sph Matrix of distances between points on S^2.
//' @param Dmat Matrix of all possible values to sum over to.
//' The values of the K-function before multiplying by the indicator function.
//' Constants may be multiplied on the resulting matrix.
// [[Rcpp::export]]
NumericMatrix engine_K(NumericVector r, NumericVector s, NumericMatrix dists_3d, NumericMatrix dists_sph, NumericMatrix Dmat) {
  int nrow_out = r.size();
  int ncol_out = s.size();
  int np = dists_3d.nrow();
  double dist_3d;
  double dist_sph;
  NumericMatrix mat_out(nrow_out, ncol_out);
  /* For each combination of r and s calculate the K-function */
  for(int ridx = 0; ridx < nrow_out; ridx++) {
    for(int sidx = 0; sidx < ncol_out; sidx++) {
      double val = 0;
      /* Extract distances */
      for(int i = 0; i < (np-1); i++) {
        for(int j = i+1; j < np; j++) {
          dist_3d = dists_3d(i, j);
          dist_sph = dists_sph(i, j);
          /* If condition is satisfied add to sum */
          if((dist_3d <= r[ridx]) & (dist_sph <= s[sidx])) {
            val += Dmat(i, j);
          }
        }
      }
      /* 2 because they are ordered pairs */
      mat_out(ridx, sidx) = 2*val;
    }
  }
  return mat_out;
}
