// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void
inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat) {
  arma::uword const n = trimat.n_cols;

  for (unsigned j = n; j-- > 0;) {
    double tmp(0.);
    for (unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// Borrowed with appreciation from Nino Hardt, Dicko Ahmadou, Ben Christofferson
// https://gallery.rcpp.org/articles/dmvnorm_arma/
arma::vec
dmvnrm_arma_fast(
    arma::mat const &x, arma::rowvec const &mean, arma::mat const &sigma,
    bool const logd = false
) {
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti    = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum    = arma::sum(log(rooti.diag())),
               constants   = -(double) xdim / 2.0 * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
List
iterate(
    arma::mat Y, List df_j, int nrep, int n, int d, double gamma, int q,
    arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha,
    double beta
) {

  // Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q * d, fill::zeros);
  List df_sim_lambda(nrep);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0)  = init.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    mat lambda_prev = df_sim_lambda[i - 1];
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      int n_i       = index_1k.n_elem;
      NumericVector Ysums;
      mat Yrows  = Y.rows(index_1k);
      Ysums      = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_prev) *
                   (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums));
      mat var_i       = inv(lambda0 + n_i * lambda_prev);
      mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(df_sim_z(i - 1, j) - 1);
    }
    mat sumofsq = (Y - mu_i_long).t() * (Y - mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv         = diagmat(beta_d);
    mat lambda_i     = rwish(n + alpha, inv(Vinv + sumofsq));
    df_sim_lambda[i] = lambda_i;
    mat sigma_i      = inv(lambda_i);

    // Update z
    df_sim_z.row(i)    = df_sim_z.row(i - 1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev         = df_sim_z(i, j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      uvec i_vector(1);
      i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev = gamma / j_vector.size() * 2 *
                       accu((df_sim_z(i_vector, j_vector) == z_j_prev)) +
                   dmvnrm_arma_fast(
                       Y.row(j), mu_i.row(z_j_prev - 1), sigma_i, true
                   )[0];
        h_z_new =
            gamma / j_vector.size() * 2 *
                accu((df_sim_z(i_vector, j_vector) == z_j_new)) +
            dmvnrm_arma_fast(Y.row(j), mu_i.row(z_j_new - 1), sigma_i, true)[0];
      } else {
        h_z_prev = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), sigma_i, true
        )[0];
        h_z_new =
            dmvnrm_arma_fast(Y.row(j), mu_i.row(z_j_new - 1), sigma_i, true)[0];
      }
      double prob_j = exp(h_z_new - h_z_prev);
      if (prob_j > 1) {
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs   = {1 - prob_j, prob_j};
      df_sim_z(i, j)        = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["plogLik"] = plogLik
  );
  return (out);
}

// [[Rcpp::export]]
List
iterate_vvv(
    arma::mat Y, List df_j, int nrep, int n, int d, double gamma, int q,
    arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha,
    double beta
) {

  // Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q * d, fill::zeros);
  List df_sim_lambda(nrep);
  List lambda_list(q);
  List sigma_list(q);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  for (int k = 1; k <= q; k++) {
    lambda_list[k - 1] = lambda0;
    sigma_list[k - 1]  = inv(lambda0);
  }
  df_sim_lambda[0] = lambda_list;
  df_sim_z.row(0)  = init.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    List lambda_prev = df_sim_lambda[i - 1];
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      int n_i       = index_1k.n_elem;
      mat Yrows     = Y.rows(index_1k);
      Ysums         = sum(Yrows, 0);
      mat lambda_k  = lambda_prev[k - 1];
      vec mean_i    = inv(lambda0 + n_i * lambda_k) *
                   (lambda0 * mu0vec + lambda_k * as<colvec>(Ysums));
      mat var_i       = inv(lambda0 + n_i * lambda_k);
      mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(df_sim_z(i - 1, j) - 1);
    }
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      int n_i       = index_1k.n_elem;
      mat sumofsq   = (Y.rows(index_1k) - mu_i_long.rows(index_1k)).t() *
                    (Y.rows(index_1k) - mu_i_long.rows(index_1k));
      mat lambda_i       = rwish(n_i + alpha, inv(Vinv + sumofsq));
      mat sigma_i        = inv(lambda_i);
      lambda_list[k - 1] = lambda_i;
      sigma_list[k - 1]  = sigma_i;
    }
    df_sim_lambda[i] = lambda_list;

    // Update z
    df_sim_z.row(i)    = df_sim_z.row(i - 1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev         = df_sim_z(i, j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      uvec i_vector(1);
      i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev =
            gamma / j_vector.size() * 2 *
                accu((df_sim_z(i_vector, j_vector) == z_j_prev)) +
            dmvnrm_arma_fast(
                Y.row(j), mu_i.row(z_j_prev - 1), sigma_list[z_j_prev - 1], true
            )[0];
        h_z_new =
            gamma / j_vector.size() * 2 *
                accu((df_sim_z(i_vector, j_vector) == z_j_new)) +
            dmvnrm_arma_fast(
                Y.row(j), mu_i.row(z_j_new - 1), sigma_list[z_j_new - 1], true
            )[0];
      } else {
        h_z_prev = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), sigma_list[z_j_prev - 1], true
        )[0];
        h_z_new = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_new - 1), sigma_list[z_j_new - 1], true
        )[0];
      }
      double prob_j = exp(h_z_new - h_z_prev);
      if (prob_j > 1) {
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs   = {1 - prob_j, prob_j};
      df_sim_z(i, j)        = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["plogLik"] = plogLik
  );
  return (out);
}

// [[Rcpp::export]]
List
iterate_t(
    arma::mat Y, List df_j, int nrep, int n, int d, double gamma, int q,
    arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha,
    double beta
) {

  // Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q * d, fill::zeros);
  List df_sim_lambda(nrep);
  mat df_sim_w(nrep, n);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0)  = init.t();
  vec w            = ones<vec>(n);
  df_sim_w.row(0)  = w.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    mat lambda_prev = df_sim_lambda[i - 1];
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      double n_i    = sum(w(index_1k));
      mat Yrows     = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums      = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_prev) *
                   (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums)
                   );                                 // posterior mean
      mat var_i = inv(lambda0 + n_i * lambda_prev);   // posterior variance
      mu_i.row(k - 1) =
          rmvnorm(1, mean_i, var_i);   // sample from posterior for mu
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(df_sim_z(i - 1, j) - 1);
    }
    mat sumofsq = (Y - mu_i_long).t() * diagmat(w) * (Y - mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv         = diagmat(beta_d);
    mat lambda_i     = rwish(n + alpha, inv(Vinv + sumofsq));
    df_sim_lambda[i] = lambda_i;
    mat sigma_i      = inv(lambda_i);

    // Update z and w
    double w_alpha = (d + 4) / 2;   // shape parameter
    double w_beta;
    df_sim_z.row(i)    = df_sim_z.row(i - 1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      w_beta = as_scalar(
          2 / ((Y.row(j) - mu_i_long.row(j)) * lambda_i *
                   (Y.row(j) - mu_i_long.row(j)).t() +
               4)
      );                                   // scale parameter
      w[j] = R::rgamma(w_alpha, w_beta);   // sample from posterior for w

      int z_j_prev         = df_sim_z(i, j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      uvec i_vector(1);
      i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev = gamma / j_vector.size() * 2 *
                       accu((df_sim_z(i_vector, j_vector) == z_j_prev)) +
                   dmvnrm_arma_fast(
                       Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
                   )[0];
        h_z_new = gamma / j_vector.size() * 2 *
                      accu((df_sim_z(i_vector, j_vector) == z_j_new)) +
                  dmvnrm_arma_fast(
                      Y.row(j), mu_i.row(z_j_new - 1), sigma_i / w[j], true
                  )[0];
      } else {
        h_z_prev = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
        )[0];
        h_z_new = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_new - 1), sigma_i / w[j], true
        )[0];
      }
      double prob_j = exp(h_z_new - h_z_prev);
      if (prob_j > 1) {
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs   = {1 - prob_j, prob_j};
      df_sim_z(i, j)        = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    df_sim_w.row(i) = w.t();
    plogLik[i]      = sum(plogLikj);
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = df_sim_w, _["plogLik"] = plogLik
  );
  return (out);
}

// [[Rcpp::export]]
List
iterate_t_vvv(
    arma::mat Y, List df_j, int nrep, int n, int d, double gamma, int q,
    arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha,
    double beta
) {

  // Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q * d, fill::zeros);
  List df_sim_lambda(nrep);
  List lambda_list(q);
  List sigma_list(q);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  for (int k = 1; k <= q; k++) {
    lambda_list[k - 1] = lambda0;
    sigma_list[k - 1]  = inv(lambda0);
  }
  df_sim_lambda[0] = lambda_list;
  df_sim_z.row(0)  = init.t();
  vec w            = ones<vec>(n);

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    List lambda_prev = df_sim_lambda[i - 1];
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      double n_i    = sum(w(index_1k));
      mat Yrows     = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums        = sum(Yrows, 0);
      mat lambda_k = lambda_prev[k - 1];
      vec mean_i =
          inv(lambda0 + n_i * lambda_k) *
          (lambda0 * mu0vec + lambda_k * as<colvec>(Ysums));   // posterior mean
      mat var_i = inv(lambda0 + n_i * lambda_k);   // posterior variance
      mu_i.row(k - 1) =
          rmvnorm(1, mean_i, var_i);   // sample from posterior for mu
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(df_sim_z(i - 1, j) - 1);
    }
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(df_sim_z.row(i - 1) == k);
      int n_i       = index_1k.n_elem;
      mat sumofsq   = (Y.rows(index_1k) - mu_i_long.rows(index_1k)).t() *
                    diagmat(w(index_1k)) *
                    (Y.rows(index_1k) - mu_i_long.rows(index_1k));
      mat lambda_i       = rwish(n_i + alpha, inv(Vinv + sumofsq));
      mat sigma_i        = inv(lambda_i);
      lambda_list[k - 1] = lambda_i;
      sigma_list[k - 1]  = sigma_i;
    }
    df_sim_lambda[i] = lambda_list;

    // Update z and w
    double w_alpha = (d + 4) / 2;   // shape parameter
    double w_beta;
    df_sim_z.row(i)    = df_sim_z.row(i - 1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev = df_sim_z(i, j);
      mat lambda_i = lambda_list[z_j_prev - 1];
      mat sigma_i  = sigma_list[z_j_prev - 1];
      w_beta       = as_scalar(
          2 / ((Y.row(j) - mu_i_long.row(j)) * lambda_i *
                   (Y.row(j) - mu_i_long.row(j)).t() +
               4)
      );                             // scale parameter
      w[j] = R::rgamma(w_alpha, w_beta);   // sample from posterior for w

      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      mat sigma_i_new      = sigma_list[z_j_new - 1];
      uvec j_vector        = df_j[j];
      uvec i_vector(1);
      i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev = gamma / j_vector.size() * 2 *
                       accu((df_sim_z(i_vector, j_vector) == z_j_prev)) +
                   dmvnrm_arma_fast(
                       Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
                   )[0];
        h_z_new = gamma / j_vector.size() * 2 *
                      accu((df_sim_z(i_vector, j_vector) == z_j_new)) +
                  dmvnrm_arma_fast(
                      Y.row(j), mu_i.row(z_j_new - 1), sigma_i_new / w[j], true
                  )[0];
      } else {
        h_z_prev = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
        )[0];
        h_z_new = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_new - 1), sigma_i_new / w[j], true
        )[0];
      }
      double prob_j = exp(h_z_new - h_z_prev);
      if (prob_j > 1) {
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs   = {1 - prob_j, prob_j};
      df_sim_z(i, j)        = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = w, _["plogLik"] = plogLik
  );
  return (out);
}

// [[Rcpp::export]]
/**
 * @brief
 *
 * @param Y
 * @param df_j
 * @param tdist
 * @param nrep
 * @param n
 * @param n0
 * @param d
 * @param gamma
 * @param q
 * @param init
 * @param subspots
 * @param verbose
 * @param jitter_scale
 * @param c
 * @param mu0
 * @param lambda0
 * @param alpha
 * @param beta
 * @return List
 */
List
iterate_deconv(
    arma::mat Y, List df_j, bool tdist, int nrep, int n, int n0, int d,
    double gamma, int q, arma::vec init, int subspots, bool verbose,
    double jitter_scale, double c, NumericVector mu0, arma::mat lambda0,
    double alpha, double beta
) {

  // Initalize matrices storing iterations
  mat Y0 = Y.rows(0, n0 - 1);
  mat df_sim_z(nrep / 100 + 1, n, fill::zeros);
  mat df_sim_mu(nrep, q * d, fill::zeros);
  List df_sim_lambda(nrep / 100 + 1);
  List df_sim_Y(nrep / 100 + 1);
  mat df_sim_w(nrep / 100 + 1, n);
  NumericVector Ychange(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  mat lambda_i     = lambda0;
  df_sim_z.row(0)  = init.t();
  vec z            = init;
  df_sim_Y[0]      = Y;
  vec w            = ones<vec>(n);
  df_sim_w.row(0)  = w.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  mat mu_i(q, d);
  mat mu_i_long(n, d);
  uvec j0_vector;
  if (subspots == 6) {
    j0_vector = {0, 1, 2, 3, 4, 5};
  } else {
    j0_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  }
  mat Y_j_prev(subspots, d);
  mat Y_j_new(subspots, d);
  mat error(subspots, d);
  vec zero_vec  = zeros<vec>(d);
  vec one_vec   = ones<vec>(d);
  mat error_var = diagmat(one_vec) / d * jitter_scale;
  Progress p(nrep - 1, verbose);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (0.24s/iter based on 6k subspots)
    if (i % 4 == 0)
      Rcpp::checkUserInterrupt();

    p.increment();
    if (i % 10 == 0) {
      if (Progress::check_abort())
        return (List::create(
            _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
            _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange
        ));
    }

    // Update mu
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      double n_i    = sum(w(index_1k));
      mat Yrows     = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums      = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_i) *
                   (lambda0 * mu0vec + lambda_i * as<colvec>(Ysums));
      mat var_i       = inv(lambda0 + n_i * lambda_i);
      mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    // Update lambda
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(z[j] - 1);
    }
    mat sumofsq = (Y - mu_i_long).t() * diagmat(w) * (Y - mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv    = diagmat(beta_d);
    lambda_i    = rwish(n + alpha, inv(Vinv + sumofsq));
    mat sigma_i = inv(lambda_i);

    // Update Y
    int updateCounter = 0;
    for (int j0 = 0; j0 < n0; j0++) {
      Y_j_prev          = Y.rows(j0_vector * n0 + j0);
      error             = rmvnorm(subspots, zero_vec, error_var);
      rowvec error_mean = sum(error, 0) / subspots;
      for (int r = 0; r < subspots; r++) {
        error.row(r) = error.row(r) - error_mean;
      }
      Y_j_new    = Y_j_prev + error;
      mat mu_i_j = mu_i_long.rows(j0_vector * n0 + j0);
      vec p_prev = {0.0};
      vec p_new  = {0.0};
      for (int r = 0; r < subspots; r++) {
        p_prev +=
            dmvnrm_arma_fast(
                Y_j_prev.row(r), mu_i_j.row(r), sigma_i / w[j0 + n0 * r], true
            ) -
            c * (accu(pow(Y_j_prev.row(r) - Y0.row(j0), 2)));
        p_new +=
            dmvnrm_arma_fast(
                Y_j_new.row(r), mu_i_j.row(r), sigma_i / w[j0 + n0 * r], true
            ) -
            c * (accu(pow(Y_j_new.row(r) - Y0.row(j0), 2)));
      }
      double probY_j = as_scalar(exp(p_new - p_prev));
      if (probY_j > 1) {
        probY_j = 1;
      }
      IntegerVector Ysample = {0, 1};
      NumericVector probsY  = {1 - probY_j, probY_j};
      int yesUpdate         = sample(Ysample, 1, true, probsY)[0];
      if (yesUpdate == 1) {
        Y.rows(j0_vector * n0 + j0) = Y_j_new;
        updateCounter++;
      }
    }
    Ychange[i] = updateCounter * 1.0 / n0;

    // Update w and z
    double w_alpha = (d + 4) / 2;   // shape parameter
    double w_beta;
    IntegerVector qvec = seq_len(q);
    for (int j = 0; j < n; j++) {
      if (tdist) {
        w_beta = as_scalar(
            2 / ((Y.row(j) - mu_i_long.row(j)) * lambda_i *
                     (Y.row(j) - mu_i_long.row(j)).t() +
                 4)
        );                                   // scale parameter
        w[j] = R::rgamma(w_alpha, w_beta);   // sample from posterior for w
      }
      int z_j_prev         = z[j];
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev =
            gamma / j_vector.size() * 2 * accu((z(j_vector) == z_j_prev)) +
            dmvnrm_arma_fast(
                Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
            )[0];
        h_z_new = gamma / j_vector.size() * 2 * accu((z(j_vector) == z_j_new)) +
                  dmvnrm_arma_fast(
                      Y.row(j), mu_i.row(z_j_new - 1), sigma_i / w[j], true
                  )[0];
      } else {
        h_z_prev = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), sigma_i / w[j], true
        )[0];
        h_z_new = dmvnrm_arma_fast(
            Y.row(j), mu_i.row(z_j_new - 1), sigma_i / w[j], true
        )[0];
      }
      double prob_j = exp(h_z_new - h_z_prev);
      if (prob_j > 1) {
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs   = {1 - prob_j, prob_j};
      z[j]                  = sample(zsample, 1, true, probs)[0];
    }
    if ((i + 1) % 100 == 0) {
      df_sim_lambda[(i + 1) / 100] = lambda_i;
      df_sim_Y[(i + 1) / 100]      = Y;
      df_sim_w.row((i + 1) / 100)  = w.t();
      df_sim_z.row((i + 1) / 100)  = z.t();
    }
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange
  );
  return (out);
}
