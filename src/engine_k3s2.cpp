#include <Rcpp.h>
using namespace Rcpp;

//' Engine for computing K-function
//'
//' Engine for calculating K-function of point processes on the product space of R^3 and S^2.
//' @param r Vector of values for the argument r at which K(r, s) should be evaluated.
//' @param s Vector of values for the argument s at which K(r, s) should be evaluated.
//' @param x_vec x coordinates
//' @param y_vec y coordinates
//' @param z_vec z coordinates
//' @param dists_sph Matrix of distances between points on S^2.
//' @param Dmat Matrix of all possible values to sum over to.
//' The values of the K-function before multiplying by the indicator function.
//' Constants may be multiplied on the resulting matrix.
//' @useDynLib spatstatspacesphereadd, .registration = TRUE
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix engine_k3s2(NumericVector r, NumericVector s, 
                          NumericVector x_vec, NumericVector y_vec, NumericVector z_vec, 
                          NumericMatrix dists_sph, 
                          NumericMatrix Dmat) {
  int nrow_out = r.size();
  int ncol_out = s.size();
  int np = x_vec.size();
  double xi, yi, zi;
  double dx, dy, dz;
  double dist_3d;
  double dist_sph;
  NumericMatrix mat_out(nrow_out, ncol_out);
  /* For each combination of r and s calculate the K-function */
  for(int ridx = 0; ridx < nrow_out; ridx++) {
    for(int sidx = 0; sidx < ncol_out; sidx++) {
      double val = 0;
      /* Extract distances */
      for(int i = 0; i < (np-1); i++) {
        xi = x_vec[i];
        yi = y_vec[i];
        zi = z_vec[i];
        for(int j = i+1; j < np; j++) {
          dx = xi - x_vec[j];
          dy = yi - y_vec[j];
          dz = zi - z_vec[j];
          dist_3d = sqrt(dx * dx + dy * dy + dz * dz);
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
