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