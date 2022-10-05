//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::vec dist(arma::mat origin, arma::mat destinations) {
	int d_length = destinations.n_rows;
	vec distance(d_length, 1, fill::zeros);
	for (int i = 0; i < d_length; i++) {
		distance[i] = sqrt(pow(origin[0] - destinations.col(0)[i], 2) + pow(origin[1] - destinations.col(1)[i], 2));
	}
	return distance;
}
//[[Rcpp::export(.catchment_areas)]]
arma::mat catchment_areas(arma::mat origins, arma::mat destinations, arma::vec supply, arma::vec demand, arma::vec fips) {
	int n_origins = origins.n_rows;
	mat catchment;
	uvec trunc0(destinations.n_rows, fill::zeros);
	for (unsigned int i = 1; i < trunc0.n_elem; i++) trunc0[i] = i;
	for (int i = 0; i < n_origins; i++) {
		mat origin   = origins.row(i);
		uvec trunc   = trunc0.rows(find(abs(origin[0] - destinations.col(0)) < 1));
		trunc        = trunc.rows(find(abs(origin[1] - destinations(trunc, (uvec){1})) < 1));
		mat distance = dist(origin, destinations.rows(trunc));
		if (trunc.n_elem > 0) {
			mat catchment_temp(trunc.n_elem, 4, fill::zeros);
			catchment_temp.col(0) = fips.rows(trunc);
			catchment_temp.col(1) = mat(trunc.n_elem, 1, fill::zeros);
			catchment_temp.col(1) += i;
			catchment_temp.col(2) = demand.rows(trunc);
			catchment_temp = catchment_temp.rows(sort_index(distance));
			catchment_temp.col(3) = cumsum(catchment_temp.col(2));
			if (max(catchment_temp.col(3)) > supply[i]) {
				uvec trunc_max = find(catchment_temp.col(3) > supply[i]);
				catchment_temp = catchment_temp.rows(0, trunc_max[0]);
			}
			if (i == 0) catchment = catchment_temp; 
			else catchment = join_vert(catchment, catchment_temp);
		}
	}
	return catchment;
}
