#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
using namespace Rcpp;

struct greater
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};

//' @export
// [[Rcpp::export]]
DataFrame drop_dense_internal(NumericVector sorted_sample,
                              NumericVector sorted_theoretical,
                              int N_hard) {
  NumericVector x = clone(sorted_sample);
  NumericVector y = clone(sorted_theoretical);
  size_t len = y.length();
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
  DataFrame df = DataFrame::create(
    Named("sorted_pruned_sample") = NumericVector(final_x.begin(), final_x.end()),
    Named("sorted_pruned_theoretical") = NumericVector(final_y.begin(), final_y.end())
  );
  return df;
}

//' @noRd
// [[Rcpp::export]]
DataFrame drop_dense_qq(NumericVector sample) {
  int N_hard = 10000;
  // Step 1:
  // Remove inf, NA, NaN, <0 and >1 and add to std::vector
  std::vector<double> x(sample.size());
  size_t k=0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (!NumericVector::is_na(sample[i]) &
        !traits::is_nan<REALSXP>(sample[i]) &
        !traits::is_infinite<REALSXP>(sample[i]) &
        (sample[i] >= 0.0) & (sample[i] <= 1.0)) {
      x[k] = -log10(sample[i]);
      k++;
    }
  }
  x.resize(k);

  // Step 2:
  // Sort reverse order and -log10.
  std::sort(x.begin(), x.end(),greater());

  // Step 3:
  // Create y
  std::vector<double> y(k);
  double max_n = (double)(k);
  for(size_t i = 0; i < k; ++i){
    y[i] = -log10(((double)(i)+0.5)/(max_n));
  }

  size_t len = k;
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
  DataFrame df = DataFrame::create(
    Named("sorted_pruned_sample") = NumericVector(final_x.begin(), final_x.end()),
    Named("sorted_pruned_theoretical") = NumericVector(final_y.begin(), final_y.end())
  );
  return df;
}
