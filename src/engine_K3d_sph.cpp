#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' Engine for computing K-function
//'
//' Engine for calculating K-function of point processes on the product space of R^3 and S^2.
//' @param r Vector of values for the argument r at which K(r, s) should be evaluated.
//' @param s Vector of values for the argument s at which K(r, s) should be evaluated.
//' @param x_vec Vector of x coordinates.
//' @param y_vec Vector of y coordinates.
//' @param z_vec Vector of z coordinates.
//' @param dists_sph Matrix of distances between points on S^2.
//' @param Dmat Matrix of all possible values to sum over to.
//' The values of the K-function before multiplying by the indicator function.
//' Constants may be multiplied on the resulting matrix.
//' @useDynLib spatstatspacesphereadd, .registration = TRUE
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix engine_K3d_sph(NumericVector r, NumericVector s, 
                             NumericVector x_vec, NumericVector y_vec, NumericVector z_vec, 
                             NumericMatrix dists_sph,
                             NumericMatrix Dmat) {
  int nrow_out = r.size();
  int ncol_out = s.size();
  int np = x_vec.size();
  double xi, yi, zi, dx, dy, dz;
  double dist_3d, dist_sph;
  double val;
  NumericMatrix mat_out(nrow_out, ncol_out);
  /* Extract distances */
  for(int i = 0; i < (np-1); i++) {
    /* Distance in 3d */
    xi = x_vec[i];
    yi = y_vec[i];
    zi = z_vec[i];
    for(int j = i+1; j < np; j++) {
      /* Distance in 3d */
      dx = xi - x_vec[j];
      dy = yi - y_vec[j];
      dz = zi - z_vec[j];
      dist_3d = sqrt(dx * dx + dy * dy + dz * dz);
      /* Distance on sphere */
      dist_sph = dists_sph(i, j);
      /* Extract the relevant term */
      val = 2 * Dmat(i, j);
      /* Add proper value to those points satisfying the condition */
      for(int ridx = 0; ridx < nrow_out; ridx++) {
        if(dist_3d <= r[ridx]) {
          for(int sidx = 0; sidx < ncol_out; sidx++) {
            if(dist_sph <= s[sidx]) {
              mat_out(ridx, sidx) += val;
            }
          }
        }
      }
    }
  }
  return(mat_out);
}
