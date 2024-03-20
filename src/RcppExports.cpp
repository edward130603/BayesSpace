// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// iterate
List iterate(const arma::mat& Y, const List& df_j, int nrep, int thin, int n, int d, double gamma, int q, const arma::uvec& init, const NumericVector& mu0, const arma::mat& lambda0, double alpha, double beta);
RcppExport SEXP _BayesSpace_iterate(SEXP YSEXP, SEXP df_jSEXP, SEXP nrepSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP dSEXP, SEXP gammaSEXP, SEXP qSEXP, SEXP initSEXP, SEXP mu0SEXP, SEXP lambda0SEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const List& >::type df_j(df_jSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate(Y, df_j, nrep, thin, n, d, gamma, q, init, mu0, lambda0, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// iterate_vvv
List iterate_vvv(const arma::mat& Y, const List& df_j, int nrep, int thin, int n, int d, double gamma, int q, const arma::uvec& init, const NumericVector& mu0, const arma::mat& lambda0, double alpha, double beta);
RcppExport SEXP _BayesSpace_iterate_vvv(SEXP YSEXP, SEXP df_jSEXP, SEXP nrepSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP dSEXP, SEXP gammaSEXP, SEXP qSEXP, SEXP initSEXP, SEXP mu0SEXP, SEXP lambda0SEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const List& >::type df_j(df_jSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate_vvv(Y, df_j, nrep, thin, n, d, gamma, q, init, mu0, lambda0, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// iterate_t
List iterate_t(const arma::mat& Y, const List& df_j, int nrep, int thin, int n, int d, double gamma, int q, const arma::uvec& init, const NumericVector& mu0, const arma::mat& lambda0, double alpha, double beta);
RcppExport SEXP _BayesSpace_iterate_t(SEXP YSEXP, SEXP df_jSEXP, SEXP nrepSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP dSEXP, SEXP gammaSEXP, SEXP qSEXP, SEXP initSEXP, SEXP mu0SEXP, SEXP lambda0SEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const List& >::type df_j(df_jSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate_t(Y, df_j, nrep, thin, n, d, gamma, q, init, mu0, lambda0, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// iterate_t_vvv
List iterate_t_vvv(const arma::mat& Y, const List& df_j, int nrep, int thin, int n, int d, double gamma, int q, const arma::uvec& init, const NumericVector& mu0, const arma::mat& lambda0, double alpha, double beta);
RcppExport SEXP _BayesSpace_iterate_t_vvv(SEXP YSEXP, SEXP df_jSEXP, SEXP nrepSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP dSEXP, SEXP gammaSEXP, SEXP qSEXP, SEXP initSEXP, SEXP mu0SEXP, SEXP lambda0SEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const List& >::type df_j(df_jSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate_t_vvv(Y, df_j, nrep, thin, n, d, gamma, q, init, mu0, lambda0, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// iterate_deconv
List iterate_deconv(const arma::mat& subspot_positions, const double dist, const CharacterVector& spot_neighbors, arma::mat& Y, bool tdist, int nrep, int thin, int n, int n0, int d, double gamma, int q, const arma::uvec& init, int subspots, bool verbose, double jitter_scale, int adapt_before, double c, const NumericVector& mu0, const arma::mat& lambda0, double alpha, double beta, int thread_num);
RcppExport SEXP _BayesSpace_iterate_deconv(SEXP subspot_positionsSEXP, SEXP distSEXP, SEXP spot_neighborsSEXP, SEXP YSEXP, SEXP tdistSEXP, SEXP nrepSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP n0SEXP, SEXP dSEXP, SEXP gammaSEXP, SEXP qSEXP, SEXP initSEXP, SEXP subspotsSEXP, SEXP verboseSEXP, SEXP jitter_scaleSEXP, SEXP adapt_beforeSEXP, SEXP cSEXP, SEXP mu0SEXP, SEXP lambda0SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP thread_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type subspot_positions(subspot_positionsSEXP);
    Rcpp::traits::input_parameter< const double >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type spot_neighbors(spot_neighborsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type tdist(tdistSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type subspots(subspotsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type jitter_scale(jitter_scaleSEXP);
    Rcpp::traits::input_parameter< int >::type adapt_before(adapt_beforeSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type thread_num(thread_numSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate_deconv(subspot_positions, dist, spot_neighbors, Y, tdist, nrep, thin, n, n0, d, gamma, q, init, subspots, verbose, jitter_scale, adapt_before, c, mu0, lambda0, alpha, beta, thread_num));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesSpace_iterate", (DL_FUNC) &_BayesSpace_iterate, 13},
    {"_BayesSpace_iterate_vvv", (DL_FUNC) &_BayesSpace_iterate_vvv, 13},
    {"_BayesSpace_iterate_t", (DL_FUNC) &_BayesSpace_iterate_t, 13},
    {"_BayesSpace_iterate_t_vvv", (DL_FUNC) &_BayesSpace_iterate_t_vvv, 13},
    {"_BayesSpace_iterate_deconv", (DL_FUNC) &_BayesSpace_iterate_deconv, 23},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesSpace(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
