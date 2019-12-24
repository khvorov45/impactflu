#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_ideal_cpp(const int& init_pop_size,
                        const IntegerVector& vaccinations,
                        const IntegerVector& cases_novac,
                        const NumericVector& ve,
                        const int& lag) {
  int nt = vaccinations.size();
  IntegerVector timepoint(nt);
  for (int i = 0; i < nt; i++) timepoint[i] = i + 1;

  IntegerVector cases(nt), popn(nt), avert(nt), A_to_E(nt), A(nt), b(nt),
    bprev(lag), B(nt), C(nt), D(nt), E(nt), F(nt);
  NumericVector pflu(nt), pvac(nt);

  pflu[0] = cases_novac[0] / double(init_pop_size);
  cases[0] = cases_novac[0];
  popn[0] = init_pop_size - cases_novac[0];
  pvac[0] = vaccinations[0] / double(init_pop_size);
  b[0] = R::rbinom(init_pop_size, pvac[0]);
  A_to_E[0] = R::rbinom(init_pop_size, pflu[0]);
  A[0] = init_pop_size - A_to_E[0] - b[0];
  B[0] = b[0];
  E[0] = A_to_E[0];

  DataFrame ideal_pop = DataFrame::create(
    _["timepoint"] = timepoint,
    _["vaccinations"] = vaccinations,
    _["cases_novac"] = cases_novac,
    _["ve"] = ve,
    _["pflu"] = pflu,
    _["cases"] = cases,
    _["popn"] = popn,
    _["pvac"] = pvac,
    _["b"] = b,
    _["A_to_E"] = A_to_E,
    _["A"] = A,
    _["B"] = B,
    _["E"] = E
  );
  return ideal_pop;
}
