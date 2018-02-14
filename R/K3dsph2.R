#' K-function
#'
#' Estimates the space-sphere K-function from a point pattern with points in \code{R^3 x S^2}
#' in a window of arbitry shape in space and on the entire sphere.
#' @param X Point pattern of class \code{\link{pp3}}.
#' @param Y Point pattern of class \code{\link{ppc}} or \code{\link{pps}}.
#' @param r Optional. Vector of values for the argument r at which K(r, s) should be evaluated.
#' @param s Optional. Vector of values for the argument s at which K(r, s) should be evaluated.
#' @param rmax Optional. Maximum desired value of the argument r.
#' @param smax Optional. Maximum desired value of the argument s.
#' @param nrval Optional. Number of values of r for which K(r, s) will be estimated.
#' @param nsval Optional. Number of values of s for which K(r, s) will be estimated.
#' @param intenssX Optional. Values of the estimated intensity function in 3d space.
#' Either a \code{NULL}, vector, matrix, or function.
#' If \code{NULL} the function will return the space-sphere K-function under the assumption of homogeneity in 3d space.
#' If a vector it should contain the intensities at each of the observed points.
#' If a matrix it should contain the product of the intensities of every pair of observed points.
#' If a function it should return a vector containing the intensities at each of the observed points.
#' @param intenssY Optional. Values of the estimated intensity function on the unitsphere.
#' Either a \code{NULL}, vector, matrix, or function.
#' If \code{NULL} the function will return the space-sphere K-function under the assumption of homogeneity on the unitsphere.
#' If a vector it should contain the intensities at each of the observed points.
#' If a matrix it should contain the product of the intensities of every pair of observed points.
#' If a function it should return a vector containing the intensities at each of the observed points.
#' @param parmsX List of additional arguments passed to \code{intenssX} if it is given as a function.
#' @param parmsY List of additional arguments passed to \code{intenssY} if it is given as a function.
#' @details If inhomogeneous it is assumed that the intensity of the space-sphere points is the product of the intensity in space and intensity on the sphere.
#' @return A list containing \code{r}, \code{s}, \code{theo}, and \code{K3dsph}.
#' \code{theo} is the theoretical space-sphere K-function under Stationary Poisson.
#' \code{theo} and \code{K3dsph} are matrices with \code{length(r)} rows and \code{length(s)} columns.
#' @import spatstat spherstat spatstatsphadd spatstat3dadd spatstatciradd
#' @export
K3dsph2 <- function(X, Y,
                    r = NULL, s = NULL, rmax = NULL, smax = NULL, nrval = 128, nsval = nrval,
                    intenssX = NULL, intenssY = NULL, parmsX, parmsY) {
  stopifnot(inherits(X, "pp3"))
  stopifnot(inherits(Y, "pps") || inherits(Y, "ppc"))
  stopifnot(npoints(X) == npoints(Y))
  
  if (is.null(r)) {
    if (is.null(rmax)) {
      rmax <- diameter(X$domain) / 2
    }
    r_vec <- seq(from = 0, to = rmax, length.out = nrval)
  } else {
    stopifnot(is.vector(r))
    r_vec <- r
  }
  if (is.null(s)) {
    if (is.null(smax)) {
      smax <- pi / 2
    }
    s_vec <- seq(from = 0, to = smax, length.out = nsval)
  } else {
    stopifnot(is.vector(s))
    s_vec <- s
  }
  
  if (inherits(Y, "ppc")) {
    out <- list(r = r_vec,
                s = s_vec,
                theo = outer(r_vec, s_vec, function(r, s) {
                  r^3 * pi^(5/2) / (gamma(5/2) * gamma(1)) * (2 * s / pi)
                }))
    
    angs <- Y$data$angs
    dists_sph <- abs(abs(outer(X = angs, Y = angs, FUN = "-")) %% (-pi))
    
    win_area_sph <- 2 * pi
  } else if (inherits(Y, "pps")) {
    out <- list(r = r_vec,
                s = s_vec,
                theo = outer(r_vec, s_vec, function(r, s) r^3 * pi^3 / (gamma(5/2) * gamma(3/2)) * (1 - cos(s))))
    
    dists_sph <- pairdistsph(pps2sp2(Y))
    
    win_area_sph <- 4 * pi
  }
  
  edge_factors_3d <- edge.Trans.pp3(X)
  np <- npoints(X)
  
  intenssX_mat <- switch(class(intenssX),
                         "NULL" = {
                           matrix(rep(np * (np - 1) / volume(X$domain)^2, np^2), ncol = np)
                         },
                         numeric = {
                           stopifnot(length(intenssX) == np)
                           tcrossprod(intenssX)
                         },
                         matrix = {
                           stopifnot(is.numeric(intenssX))
                           stopifnot(all(c(np, np) == dim(intenssX)))
                           intenssX
                         },
                         "function" = {
                           tcrossprod(apply(data.frame(X$data$x, X$data$y, X$data$z),
                                            MARGIN = 1,
                                            FUN = function(x) do.call(intenssX, c(list(x), parmsX))))
                         },
                         {
                           stop("intenssX should be either NULL, a vector, a matrix, or a function.")
                         })
  
  intenssY_mat <- switch(class(intenssY), 
                         "NULL" = {
                           matrix(rep(np * (np - 1) / win_area_sph^2, np^2), ncol = np)
                         },
                         numeric = {
                           stopifnot(length(intenssY) == np)
                           tcrossprod(intenssY)
                         },
                         matrix = {
                           stopifnot(is.numeric(intenssY))
                           stopifnot(all(c(np, np) == dim(intenssY)))
                           intenssY
                         },
                         "function" = {
                           tcrossprod(apply(data.frame(Y$data$long, Y$data$lat),
                                            MARGIN = 1,
                                            FUN = function(x) do.call(intenssY, c(list(x), parmsY))))
                         },
                         {
                           stop("intenssY should be either NULL, a vector, a matrix, or a function.")
                         })
  
  tmp_mat <- edge_factors_3d / (intenssX_mat * intenssY_mat / np^2)
  K <- engine_K3d_sph(r = r_vec, s = s_vec, 
                      x_vec = X$data$x, y_vec = X$data$y, z_vec = X$data$z, 
                      dists_sph = dists_sph,
                      Dmat = tmp_mat) 
  out$K3dsph <- K / win_area_sph
  
  return(out)
}
