// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sortllr
arma::uvec sortllr(const arma::mat& Xmat, const std::vector <std::vector<int>>& mb, const int& numUpdates, const std::string& family, const double& df);
RcppExport SEXP _scorelingam_sortllr(SEXP XmatSEXP, SEXP mbSEXP, SEXP numUpdatesSEXP, SEXP familySEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< const std::vector <std::vector<int>>& >::type mb(mbSEXP);
    Rcpp::traits::input_parameter< const int& >::type numUpdates(numUpdatesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family(familySEXP);
    Rcpp::traits::input_parameter< const double& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(sortllr(Xmat, mb, numUpdates, family, df));
    return rcpp_result_gen;
END_RCPP
}
// getWeights
arma::mat getWeights(const arma::mat& X, const std::vector <std::vector<unsigned int>>& pa);
RcppExport SEXP _scorelingam_getWeights(SEXP XSEXP, SEXP paSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const std::vector <std::vector<unsigned int>>& >::type pa(paSEXP);
    rcpp_result_gen = Rcpp::wrap(getWeights(X, pa));
    return rcpp_result_gen;
END_RCPP
}
// corrMat
arma::mat corrMat(const arma::mat& Xmat);
RcppExport SEXP _scorelingam_corrMat(SEXP XmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xmat(XmatSEXP);
    rcpp_result_gen = Rcpp::wrap(corrMat(Xmat));
    return rcpp_result_gen;
END_RCPP
}
// getResids
arma::mat getResids(const arma::mat& X, const arma::mat& B);
RcppExport SEXP _scorelingam_getResids(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(getResids(X, B));
    return rcpp_result_gen;
END_RCPP
}
// getRsq
arma::vec getRsq(const arma::mat& X, const arma::mat& B);
RcppExport SEXP _scorelingam_getRsq(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(getRsq(X, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scorelingam_sortllr", (DL_FUNC) &_scorelingam_sortllr, 5},
    {"_scorelingam_getWeights", (DL_FUNC) &_scorelingam_getWeights, 2},
    {"_scorelingam_corrMat", (DL_FUNC) &_scorelingam_corrMat, 1},
    {"_scorelingam_getResids", (DL_FUNC) &_scorelingam_getResids, 2},
    {"_scorelingam_getRsq", (DL_FUNC) &_scorelingam_getRsq, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_scorelingam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
