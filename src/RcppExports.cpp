// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calcCrossMean
arma::mat calcCrossMean(arma::mat& geno, arma::vec& a, arma::vec& d, arma::uword ploidy);
RcppExport SEXP _gpcp_calcCrossMean(SEXP genoSEXP, SEXP aSEXP, SEXP dSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(calcCrossMean(geno, a, d, ploidy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpcp_calcCrossMean", (DL_FUNC) &_gpcp_calcCrossMean, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpcp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
