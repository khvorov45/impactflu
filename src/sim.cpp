#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_reference_cpp(const int init_pop_size,
                            const IntegerVector& vaccinations,
                            const IntegerVector& infections_novac,
                            const NumericVector& ve,
                            const int lag) {
  int nt = vaccinations.size();
  IntegerVector timepoint(nt);
  for (int i = 0; i < nt; i++) timepoint[i] = i + 1;

  IntegerVector popn(nt), A(nt), b(nt), b_og(nt),
    C(nt), D(nt), E(nt), F(nt), infections(nt), avert(nt);
  NumericVector pflu(nt), pvac(nt);

  pflu[0] = infections_novac[0] / double(init_pop_size);
  int A_to_E = R::fround(init_pop_size * pflu[0], 0);
  infections[0] = A_to_E;
  popn[0] = init_pop_size - infections_novac[0];
  avert[0] = infections_novac[0] - infections[0];
  pvac[0] = vaccinations[0] / double(init_pop_size);
  b[0] = R::fround(init_pop_size * pvac[0], 0);
  b_og[0] = b[0];
  A[0] = init_pop_size - A_to_E - b[0];
  if (lag == 0) {
    C[0] = b[0] - R::fround(b[0] * ve[0], 0);
    D[0] = b[0] - C[0];
  }
  E[0] = A_to_E;

  for (int i = 1; i < nt; i++) {
    pflu[i] = infections_novac[i] / double(popn[i - 1]);
    A_to_E = R::fround(A[i - 1] * pflu[i], 0);
    infections[i] = A_to_E;
    F[i] = F[i - 1];
    for (int j = 1; (j <= lag) && (i - j >= 0); j++) {
      int bimj_to_F = R::fround(b[i - j] * pflu[i], 0);
      b[i - j] -= bimj_to_F;
      F[i] += bimj_to_F;
      infections[i] += bimj_to_F;
    }
    int C_to_F = R::fround(C[i - 1] * pflu[i], 0);
    infections[i] += C_to_F;
    popn[i] = popn[i - 1] - infections_novac[i];
    avert[i] = infections_novac[i] - infections[i];
    pvac[i] = double(vaccinations[i]) / (A[i - 1] + E[i - 1]);
    b[i] = R::fround(A[i - 1] * pvac[i], 0);
    b_og[i] = b[i];
    A[i] = A[i - 1] - A_to_E - b[i];
    int blag_to_C, blag_to_D;
    if (i - lag >= 0) {
      blag_to_C = b[i - lag] - R::fround(b[i - lag] * ve[i], 0);
      blag_to_D = b[i - lag] - blag_to_C;
    } else blag_to_C = blag_to_D = 0;
    C[i] = C[i - 1] - C_to_F + blag_to_C;
    D[i] = D[i - 1] + blag_to_D;
    int E_to_F = R::fround(E[i - 1] * pvac[i], 0);
    E[i] = E[i - 1] + A_to_E - E_to_F;
    F[i] += C_to_F + E_to_F;
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
    _["infections"] = infections,
    _["avert"] = avert
  );
  return ideal_pop;
}
