#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame method1_cpp(const int init_pop_size,
                      const IntegerVector& vaccinations,
                      const IntegerVector& cases,
                      const NumericVector& ve) {
  int nt = cases.size();
  IntegerVector pops(nt), popn(nt), cases_novac(nt), avert(nt);
  NumericVector pvac(nt), pflu(nt), vc_lag(nt);
  for (int i = 0; i < nt; i++)
    pvac[i] = vaccinations[i] / double(init_pop_size);
  vc_lag[0] = pvac[0] / 2;
  pops[0] = init_pop_size * (1 - vc_lag[0] * ve[0]);
  pflu[0] = double(cases[0]) / pops[0];
  popn[0] = init_pop_size - cases[0];
  cases_novac[0] = pflu[0] * popn[0];
  avert[0] = cases_novac[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    vc_lag[i] = (pvac[i] + pvac[i - 1]) / 2;
    pops[i] = (pops[i - 1] - cases[i - 1]) * (1 - vc_lag[i] * ve[i]);
    pflu[i] = double(cases[i]) / pops[i];
    popn[i] = popn[i - 1] - cases_novac[i - 1];
    cases_novac[i] = pflu[i] * popn[i];
    avert[i] = cases_novac[i] - cases[i];
  }
  DataFrame method1 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["vc_lag"] = vc_lag,
    _["pops"] = pops,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["cases_novac"] = cases_novac,
    _["avert"] = avert
  );
  return method1;
}
