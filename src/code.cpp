#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
IntegerVector drop_dense_internal(NumericVector sorted_sample, NumericVector sorted_theoretical) {
  NumericVector x = clone(sorted_sample);
  NumericVector y = clone(sorted_theoretical);
  // Do the computation for comparing the speed gain
  IntegerVector keep_inds = {1,2,3};
  return keep_inds;
}


