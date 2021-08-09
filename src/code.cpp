#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
DataFrame drop_dense_internal(NumericVector sorted_sample,
                              NumericVector sorted_theoretical,
                              int N_hard) {
  NumericVector x = clone(sorted_sample);
  NumericVector y = clone(sorted_theoretical);
  size_t len = x.length();
  double minx = x[len-1];
  double maxx = x[0];
  double miny = y[len-1];
  double maxy = y[0];

  double x_width = maxx - minx;
  double y_width = maxy - miny;
  double distThreshold = std::min(x_width, y_width)/((double)(N_hard));

  // Start with the first and last
  std::vector<double> final_x = {x[0],x[len-1]};
  std::vector<double> final_y = {y[0],y[len-1]};

  double d = (x[0] - x[1]) + (y[0] - y[1]);
  for(size_t i = 1; i < len-1; i++){
    double d_x = (x[i] - x[i+1]) + (y[i] - y[i+1]);
    if(d < distThreshold){
      d = d + d_x;
    }else{
      final_x.push_back(x[i]);
      final_y.push_back(y[i]);
      d = d_x;
    }
  }
  // Do the computation for comparing the speed gain
  DataFrame df = DataFrame::create( Named("sorted_pruned_sample") = NumericVector(final_x.begin(), final_x.end()),         // simple assign
                                    Named("sorted_pruned_theoretical") = NumericVector(final_y.begin(), final_y.end())); // using clone()
  return df;
}


