#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix simulate_DENV_long_decay_titres_cpp(
    List inf_times_list,
    NumericVector decay_rate_vec,
    double log_first_peak,
    double phi,
    double beta,
    int final_age)
{
  int n_people = inf_times_list.size();

  NumericMatrix titres(n_people, final_age);
  std::fill(titres.begin(), titres.end(), 5.0);

  for (int i = 0; i < n_people; ++i) {

    IntegerVector inf_times = inf_times_list[i];
    int n_inf = inf_times.size();
    if (n_inf == 0) continue;

    for (int inf_idx = 0; inf_idx < n_inf; ++inf_idx) {

      int inf_age = inf_times[inf_idx];
      int end_age = (inf_idx < n_inf - 1) ?
        inf_times[inf_idx + 1] - 1: final_age;

      double A0;
      if (inf_idx == 0)
      {
        A0 = 10.0 * std::pow(2.0, log_first_peak - 1.0);
      } else {
        double log_A0 =
          phi - (phi - log_first_peak) * std::exp(-beta * inf_idx);
        A0 = 10.0 * std::pow(2.0, log_A0 - 1.0);
      }

      double decay_rate = decay_rate_vec[std::min(inf_idx, 3)];
      int len = end_age - inf_age + 1;

      for (int t = 0; t < len; ++t)
      {
        double titre = A0 * std::exp(-decay_rate * t);
        titres(i, inf_age - 1 + t) = titre;
      }
    }
  }

  return titres;
}
