#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "progress.h"

//[[Rcpp::export]]
arma::imat rmultinom(int n_boot, int num, int length, arma::vec weight) {
  weight /= sum(weight);
  int lw = weight.n_elem;
  ivec output(lw);
  imat out_full(lw, n_boot);
  for (int i = 0; i < n_boot; i++) {
    R::rmultinom(num, weight.begin(), length, output.begin());
    out_full.col(i) = output;
  }
  return out_full;
}

//[[Rcpp::export]]
arma::imat get_demboot(
  arma::vec agg_index, 
  arma::vec num, 
  arma::vec weight, 
  int n_boot, 
  int n_dest, 
  int n_region
) {
  Rcout << "Generating de-aggregated demand...\n";
  imat demand_boot(n_dest, n_boot);
  for (int i = 0; i < n_region; i++) {
    uvec ind = find(agg_index == i);
    demand_boot.rows(ind) = rmultinom(n_boot, num[i], ind.n_elem, weight.rows(ind));
    if (i < n_region - 1) progress(i, n_region); else progress(49, 50);
  }
  return demand_boot;
}