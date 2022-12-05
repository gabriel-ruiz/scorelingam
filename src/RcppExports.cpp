// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scorelingam_
arma::uvec scorelingam_(const arma::mat& Xmat, const std::vector <std::vector<int>>& mb, const int& numUpdates, const std::string& family, const double& df);
RcppExport SEXP _scorelingam_scorelingam_(SEXP XmatSEXP, SEXP mbSEXP, SEXP numUpdatesSEXP, SEXP familySEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< const std::vector <std::vector<int>>& >::type mb(mbSEXP);
    Rcpp::traits::input_parameter< const int& >::type numUpdates(numUpdatesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const double& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(scorelingam_(Xmat, mb, numUpdates, family, df));
    return rcpp_result_gen;
END_RCPP
}
// getWeights_
arma::mat getWeights_(const arma::mat& X, const std::vector <std::vector<unsigned int>>& pa);
RcppExport SEXP _scorelingam_getWeights_(SEXP XSEXP, SEXP paSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const std::vector <std::vector<unsigned int>>& >::type pa(paSEXP);
    rcpp_result_gen = Rcpp::wrap(getWeights_(X, pa));
    return rcpp_result_gen;
END_RCPP
}
// corrMat_
arma::mat corrMat_(const arma::mat& Xmat);
RcppExport SEXP _scorelingam_corrMat_(SEXP XmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    rcpp_result_gen = Rcpp::wrap(corrMat_(Xmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scorelingam_scorelingam_", (DL_FUNC) &_scorelingam_scorelingam_, 5},
    {"_scorelingam_getWeights_", (DL_FUNC) &_scorelingam_getWeights_, 2},
    {"_scorelingam_corrMat_", (DL_FUNC) &_scorelingam_corrMat_, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_scorelingam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}