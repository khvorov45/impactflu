#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_ideal_cpp(const int& init_pop_size,
                        const IntegerVector& vaccinations,
                        const IntegerVector& cases_novac,
                        const NumericVector& ve,
                        const int& lag) {
  int nt = vaccinations.size();
  // Day counter
  IntegerVector timepoint(nt);
  for (int i = 0; i < nt; i++) timepoint[i] = i + 1;

  DataFrame ideal_pop = DataFrame::create(
    _["timepoint"] = timepoint
  );
  return ideal_pop;
}
