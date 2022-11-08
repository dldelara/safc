#include <RcppArmadillo.h>
#include "progress.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec calc_avail(
  arma::vec supply,
  arma::mat demand, 
  arma::field<arma::uvec> distmat,
  arma::vec agg_index,
  int n_region,
  int n_origin
) {
  vec avail(n_region);
  for (int i = 0; i < n_origin; i++) {
    uvec dist_i  = distmat[i];
    int n_dist_i = dist_i.n_elem;
    if (n_dist_i > 0) {
      mat catchment(n_dist_i, 2);
      catchment.col(0) = agg_index.rows(dist_i);
      catchment.col(1) = demand.rows(dist_i);
      mat cumdemand = cumsum(catchment.col(1));
      if (cumdemand[n_dist_i - 1] > supply[i]) {
        uvec trunc_max = find(cumdemand > supply[i]);
        catchment = catchment.rows(0, trunc_max[0]);
      }
      for (unsigned int i = 0; i < catchment.n_rows; i++) {
        avail[catchment(i, 0)] += catchment(i, 1);
      }
    }
  }
  return avail;
}

//[[Rcpp::export]]
arma::vec get_availability(
  arma::vec supply, 
  arma::mat demand, 
  arma::field<arma::uvec> distmat,
  arma::vec agg_index, 
  int n_region, 
  int n_origin
) {
  int n_boot = demand.n_cols;
  if (n_boot > 1) {
    Rcout << "Generating median availability estimates...\n";
    mat avail(n_region, n_boot);
    for (int i = 0; i < n_boot; i++) {
      avail.col(i) = calc_avail(supply, demand.col(i), distmat, agg_index, n_region, n_origin);
      progress(i, n_boot);
    }
    return median(avail, 1);
  } else {
    return calc_avail(supply, demand, distmat, agg_index, n_region, n_origin);
  }
}
