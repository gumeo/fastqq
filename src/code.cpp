#include <Rcpp.h>
using namespace Rcpp;

struct greater
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};

//' @noRd
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
  std::vector<double> final_x;
  final_x.push_back(x[0]);
  final_x.push_back(x[len-1]);
  std::vector<double> final_y;
  final_y.push_back(y[0]);
  final_y.push_back(y[len-1]);

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
    Named("x") = NumericVector(final_x.begin(), final_x.end()),
    Named("y") = NumericVector(final_y.begin(), final_y.end())
  );
  return df;
}

//' @noRd
// [[Rcpp::export]]
DataFrame drop_dense_qq(NumericVector sample, int N_hard) {
  // Step 1:
  // Remove inf, NA, NaN, <0 and >1 and add to std::vector
  // This is to be consistent with the qqman::qq functionality, this
  // is considerably faster in cpp.
  std::vector<double> x(sample.size());
  size_t k=0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (!NumericVector::is_na(sample[i]) &&
        !traits::is_nan<REALSXP>(sample[i]) &&
        !traits::is_infinite<REALSXP>(sample[i]) &&
        (sample[i] >= 0.0) && (sample[i] <= 1.0)) {
      x[k] = -log10(sample[i]);
      k++;
    }
  }
  x.resize(k);

  // Step 2:
  // Sort reverse order
  std::sort(x.begin(), x.end(),greater());

  // Step 3:
  // Create y, -log10 transformed theoretical uniform quantiles.
  std::vector<double> y(k);
  double max_n = (double)(k);
  for(size_t i = 0; i < k; ++i){
    y[i] = -log10(((double)(i)+0.5)/(max_n));
  }

  // Step 4:
  // Prune.
  size_t len = k;
  double minx = x[len-1];
  double maxx = x[0];
  double miny = y[len-1];
  double maxy = y[0];

  double x_width = maxx - minx;
  double y_width = maxy - miny;
  double distThreshold = std::min(x_width, y_width)/((double)(N_hard));

  // Start with the first and last
  std::vector<double> final_x;
  final_x.push_back(x[0]);
  final_x.push_back(x[len-1]);
  std::vector<double> final_y;
  final_y.push_back(y[0]);
  final_y.push_back(y[len-1]);

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

  DataFrame df = DataFrame::create(
    Named("sorted_pruned_sample") = NumericVector(final_x.begin(), final_x.end()),
    Named("sorted_pruned_theoretical") = NumericVector(final_y.begin(), final_y.end())
  );
  return df;
}
