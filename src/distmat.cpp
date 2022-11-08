#include <RcppArmadillo.h>
#include "progress.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec dist(arma::mat origin, arma::mat destinations) {
  int d_length = destinations.n_rows;
  vec distance(d_length, 1, fill::zeros);
  for (int i = 0; i < d_length; i++) {
    distance[i] = sqrt(pow(origin[0] - destinations.col(0)[i], 2) + pow(origin[1] - destinations.col(1)[i], 2));
  }
  return distance;
}

//[[Rcpp::export]]
arma::field<arma::uvec> get_distmat(arma::mat origins, arma::mat destinations, double max_dist, int n_origin) {
  Rcout << "Generating distance matrix...\n";
  field<uvec> distmat(n_origin);
  for (int i = 0; i < n_origin; i++) {
    uvec ind   = find(abs(origins(i, 0) - destinations.col(0)) <= max_dist);
    ind        = ind.rows(find(abs(origins(i, 1) - destinations(ind, (uvec){1})) <= max_dist));
    vec dist_i = dist(origins.row(i), destinations.rows(ind));
    uvec trunc = find(dist_i <= max_dist);
    ind        = ind.rows(trunc);
    dist_i     = dist_i.rows(trunc);
    ind        = ind.rows(sort_index(dist_i));
    distmat[i] = ind;
    progress(i, n_origin);
  }
  return distmat;
}