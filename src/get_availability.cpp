//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void progress(int s, int n_iter) {
  if (s == 0) {
    String progress(" Progress: |..................................................|");
    Rcout << progress.get_cstring() << "\r";
  }
  if ((s + 1) % (n_iter / 50) == 0) {
    String progress(" Progress: |..................................................|");
    for (int i = 0; i < round((s + 1) * 50.0 / n_iter); i++) progress.replace_first(".", "*");
    Rcout << progress.get_cstring() << "\r";
  }
  if (s == n_iter - 1) Rcout << "\n";
}

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

