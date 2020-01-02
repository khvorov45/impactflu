// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// method1_cpp
DataFrame method1_cpp(const int init_pop_size, const IntegerVector& vaccinations, const IntegerVector& cases, const NumericVector& ve);
RcppExport SEXP _impactflu_method1_cpp(SEXP init_pop_sizeSEXP, SEXP vaccinationsSEXP, SEXP casesSEXP, SEXP veSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type init_pop_size(init_pop_sizeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type vaccinations(vaccinationsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cases(casesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type ve(veSEXP);
    rcpp_result_gen = Rcpp::wrap(method1_cpp(init_pop_size, vaccinations, cases, ve));
    return rcpp_result_gen;
END_RCPP
}
// sim_ideal_cpp
DataFrame sim_ideal_cpp(const int init_pop_size, const IntegerVector& vaccinations, const IntegerVector& cases_novac, const NumericVector& ve, const int lag, bool deterministic);
RcppExport SEXP _impactflu_sim_ideal_cpp(SEXP init_pop_sizeSEXP, SEXP vaccinationsSEXP, SEXP cases_novacSEXP, SEXP veSEXP, SEXP lagSEXP, SEXP deterministicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type init_pop_size(init_pop_sizeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type vaccinations(vaccinationsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cases_novac(cases_novacSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type ve(veSEXP);
    Rcpp::traits::input_parameter< const int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< bool >::type deterministic(deterministicSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ideal_cpp(init_pop_size, vaccinations, cases_novac, ve, lag, deterministic));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_impactflu_method1_cpp", (DL_FUNC) &_impactflu_method1_cpp, 4},
    {"_impactflu_sim_ideal_cpp", (DL_FUNC) &_impactflu_sim_ideal_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_impactflu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
