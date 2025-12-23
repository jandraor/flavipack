#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void truncate_at_zero(NumericMatrix &x) {
  const int n = x.nrow();
  const int p = x.ncol();

  for (int j = 0; j < p; ++j) {
    double* col = &x(0, j);
    for (int i = 0; i < n; ++i) {
      if (col[i] < 0.0) col[i] = 0.0;
    }
  }
}
