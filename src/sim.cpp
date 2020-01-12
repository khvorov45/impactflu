#include <Rcpp.h>
using namespace Rcpp;

int my_rbinom(int n, double p, bool deterministic) {
  if (deterministic) return int(R::fround(n * p, 0));
  return int(R::rbinom(n, p));
}

// [[Rcpp::export]]
DataFrame sim_reference_cpp(const int init_pop_size,
                            const IntegerVector& vaccinations,
                            const IntegerVector& infections_novac,
                            const NumericVector& ve,
                            const int lag,
                            const bool deterministic) {
  int nt = vaccinations.size();
  IntegerVector timepoint(nt);
  for (int i = 0; i < nt; i++) timepoint[i] = i + 1;

  IntegerVector popn(nt), b(nt), b_og(nt), A(nt),
    C(nt), D(nt), E(nt), F(nt), cases(nt), avert(nt);
  NumericVector pflu(nt), pvac(nt);

  pflu[0] = infections_novac[0] / double(init_pop_size);
  popn[0] = init_pop_size - infections_novac[0];
  pvac[0] = vaccinations[0] / double(init_pop_size);
  b[0] = my_rbinom(init_pop_size, pvac[0], deterministic);
  b_og[0] = b[0];
  int A_to_E = my_rbinom(init_pop_size, pflu[0], deterministic);
  A[0] = init_pop_size - A_to_E - b[0];
  if (lag == 0) {
    C[0] = my_rbinom(b[0], 1 - ve[0], deterministic);
    D[0] = b[0] - C[0];
  }
  E[0] = A_to_E;
  cases[0] = A_to_E;
  avert[0] = infections_novac[0] - cases[0];

  for (int i = 1; i < nt; i++) {
    pflu[i] = infections_novac[i] / double(popn[i - 1]);
    popn[i] = popn[i - 1] - infections_novac[i];
    pvac[i] = double(vaccinations[i]) / (A[i - 1] + E[i - 1]);
    b[i] = my_rbinom(A[i - 1], pvac[i], deterministic);
    b_og[i] = b[i];
    A_to_E = my_rbinom(A[i - 1], pflu[i], deterministic);
    cases[i] = A_to_E;
    A[i] = A[i - 1] - A_to_E - b[i];
    F[i] = F[i - 1];
    for (int j = 1; (j <= lag) && (i - j >= 0); j++) {
      int bimj_to_F = my_rbinom(b[i - j], pflu[i], deterministic);
      b[i - j] -= bimj_to_F;
      F[i] += bimj_to_F;
      cases[i] += bimj_to_F;
    }
    int blag_to_C, blag_to_D;
    if (i - lag >= 0) {
      blag_to_C = my_rbinom(b[i - lag], 1 - ve[i], deterministic);
      blag_to_D = b[i - lag] - blag_to_C;
    } else blag_to_C = blag_to_D = 0;
    int C_to_F = my_rbinom(C[i - 1], pflu[i], deterministic);
    C[i] = C[i - 1] - C_to_F + blag_to_C;
    D[i] = D[i - 1] + blag_to_D;
    int E_to_F = vaccinations[i] - b[i];
    E[i] = E[i - 1] + A_to_E - E_to_F;
    F[i] += E_to_F + C_to_F;
    cases[i] += C_to_F;
    avert[i] = infections_novac[i] - cases[i];
  }

  DataFrame ideal_pop = DataFrame::create(
    _["timepoint"] = timepoint,
    _["vaccinations"] = vaccinations,
    _["infections_novac"] = infections_novac,
    _["ve"] = ve,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["pvac"] = pvac,
    _["b"] = b,
    _["b_og"] = b_og,
    _["A"] = A,
    _["C"] = C,
    _["D"] = D,
    _["E"] = E,
    _["F"] = F,
    _["cases"] = cases,
    _["avert"] = avert
  );
  return ideal_pop;
}
