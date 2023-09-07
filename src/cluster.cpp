// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include "double_states_vector.h"
#include "neighbor.h"
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <chrono>
#include <csignal>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]

static double const log2pi              = std::log(2.0 * M_PI);
static volatile sig_atomic_t early_stop = 0;

template <typename T>
void
print_thread_hits(const std::vector<T> &arr) {
  if (arr.size() > 0) {
    for (size_t i = 0; i < arr.size(); i++)
      std::cout << "[DEBUG] Thread " << i << " is hit " << arr[i]
                << " times.\n";
    std::cout << std::endl;
  }
}

static void
sig_handler(int _) {
  (void) _;
  std::cerr << "\nStopping..." << std::endl;

  early_stop = 1;
}

/* C++ version of the dtrmv BLAS function */
void
inplace_tri_mat_mult_t(arma::rowvec &x, arma::mat const &trimat) {
  arma::uword const n = trimat.n_cols;

  for (unsigned i = 0; i < n; i++) {
    double tmp(0.);
    for (unsigned j = i; j < n; j++)
      tmp += trimat.at(i, j) * x[j];
    x[i] = tmp;
  }
}

// Borrowed with appreciation from Nino Hardt, Dicko Ahmadou, Ben Christofferson
// https://gallery.rcpp.org/articles/dmvnorm_arma/
// Use covariance matrix
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
    inplace_tri_mat_mult_t(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

// Use precision matrix
arma::vec
dmvnrm_prec_arma_fast(
    arma::mat const &x, arma::rowvec const &mean, arma::mat const &lambda,
    bool const logd = false
) {
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti    = arma::trimatu(arma::chol(lambda));
  double const rootisum    = arma::sum(log(rooti.diag())),
               constants   = -(double) xdim / 2.0 * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult_t(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

void
convert_neighbors(const List &ori, std::vector<Neighbor> &des) {
  for (auto i = 0; i != ori.size(); i++) {
    des.emplace_back(Neighbor(NumericVector(ori[i])));
  }
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

/**
 * @brief
 *
 * @param Y the initialized principal components (num_subspots * num_pcs)
 * @param df_j the indicies of neighbors of each spot (indices adjusted to
 * 0-based)
 * @param tdist whether to use multivariate t distribution or not
 * @param nrep the number of MCMC iterations
 * @param n the number of subspots (after deconvolution)
 * @param n0 the number of spots (before deconvolution)
 * @param d the number of PCs
 * @param gamma smoothing parameter
 * @param q the number of clusters
 * @param init the initialized clustering of subspots
 * @param subspots the number of subspots of each spot
 * @param verbose whether to print more information
 * @param jitter_scale the amount of jittering (the variance) for the proposal
 * distribution
 * @param c the amount of jittering (the variance) for the prior
 * distribution
 * @param mu0 the mean hyperparameetr of mu
 * @param lambda0 the precision hyperparameter of mu
 * @param alpha one of the hyperparamters of lambda
 * @param beta one of the hyperparamters of lambda
 * @param thread_num the number of threads to be used
 * @param verbose information for debugging
 * @return List MCMC samples of latent variables in a list
 */
// [[Rcpp::export]]
List
iterate_deconv(
    arma::mat Y, List df_j, bool tdist, int nrep, int n, int n0, int d,
    double gamma, int q, arma::uvec init, int subspots, bool verbose,
    double jitter_scale, double c, NumericVector mu0, arma::mat lambda0,
    double alpha, double beta, int thread_num = 1
) {
  std::vector<int> thread_hits;

#ifdef _OPENMP
  omp_set_max_active_levels(2);
  omp_set_num_threads(thread_num);

  for (int i = 0; i < thread_num; i++)
    thread_hits.emplace_back(0);

  if (verbose) {
    std::cout << "[DEBUG] The number of threads is " << thread_num << std::endl;
  }
#endif

  // Initalize matrices storing iterations
  const mat Y0        = Y.rows(0, n0 - 1);   // The input PCs on spot level.
  mat Y_new           = mat(Y.n_rows,
                            Y.n_cols);   // The proposed PCs on subspot level.
  uvec z_new          = uvec(n);         // The proposed zs on subspot level.
  vec acceptance_prob = vec(n);   // The probability of accepting the proposals
                                  // on subspot level.
  DoubleStatesVector<double> log_likelihoods(n
  );   // The log-likelihoods on subspot level.
  umat df_sim_z(nrep / 100 + 1, n, fill::zeros);
  mat df_sim_mu(nrep, d * q, fill::zeros);
  List df_sim_lambda(nrep / 100 + 1);
  List df_sim_Y(nrep / 100 + 1);
  mat df_sim_w(nrep / 100 + 1, n);
  NumericVector Ychange(nrep, NA_REAL);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  df_sim_mu.row(0) = rowvec(rep(mu0, q));
  df_sim_lambda[0] = lambda0;
  const mat Vinv   = diagmat(vec(d, fill::value(beta)));
  mat lambda_i     = lambda0;
  df_sim_z.row(0)  = init.t();
  uvec z           = init;
  df_sim_Y[0]      = Y;
  vec w            = ones<vec>(n);
  df_sim_w.row(0)  = w.t();
  std::vector<Neighbor> neighbors;
  convert_neighbors(df_j, neighbors);

  // Iterate
  const colvec mu0vec = as<colvec>(mu0);
  mat mu_i(q, d);
  mat mu_i_long(n, d);
  uvec j0_vector;
  if (subspots == 6) {
    j0_vector = {0, 1, 2, 3, 4, 5};
  } else {
    j0_vector = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  }
  mat error(n, d);
  const double w_alpha     = (d + 4) / 2;   // shape parameter
  const IntegerVector qvec = seq_len(q);
  const vec zero_vec       = zeros<vec>(d);
  const vec one_vec        = ones<vec>(d);
  const mat error_var =
      diagmat(one_vec) / d *
      jitter_scale;   // d here does not seem to have an effect

  // Progress bar
  indicators::show_console_cursor(false);
  indicators::ProgressBar pb{
      indicators::option::MaxProgress{nrep - 1},
      indicators::option::BarWidth{50},
      indicators::option::Start{" ["},
      indicators::option::Fill{"█"},
      indicators::option::Lead{"█"},
      indicators::option::Remainder{"-"},
      indicators::option::End{"]"},
      indicators::option::PrefixText{"Enhancing"},
      indicators::option::ForegroundColor{indicators::Color::blue},
      indicators::option::ShowElapsedTime{true},
      indicators::option::ShowRemainingTime{true},
      indicators::option::FontStyles{
          std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
      }
  };

  // Keyboard interruption
  signal(SIGABRT, sig_handler);

#pragma omp parallel shared(                                                   \
        early_stop, Y, Y_new, error, acceptance_prob, j0_vector, subspots,     \
            zero_vec, mu_i_long, Y0, lambda_i, n0, c, thread_hits, n, z,       \
            z_new, neighbors, gamma, mu_i, w, log_likelihoods                  \
)
  {
    for (int i = 1; i < nrep; i++) {
#pragma omp barrier
#pragma omp single
      {
        pb.tick();

        // Update mu
        NumericVector Ysums;
        for (int k = 1; k <= q; k++) {
          const uvec index_1k = find(z == k);
          mat Yrows           = Y.rows(index_1k);
          Yrows.each_col() %= w(index_1k);
          Ysums           = sum(Yrows, 0);
          const mat var_i = inv(lambda0 + sum(w(index_1k)) * lambda_i);
          const vec mean_i =
              var_i * (lambda0 * mu0vec + lambda_i * as<colvec>(Ysums));
          mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
        }
        df_sim_mu.row(i) = vectorise(mu_i, 1);

        // Update lambda
        for (int j = 0; j < n; j++) {
          mu_i_long.row(j) = mu_i.row(z(j) - 1);
        }
        lambda_i = rwish(
            n + alpha,
            inv(Vinv + (Y - mu_i_long).t() * diagmat(w) * (Y - mu_i_long))
        );

        // Propose new values for Y.
        error = rmvnorm(
            n, zero_vec, error_var
        );   // Generate random numbers before entering multithreading.
      }

      // Multithreading to compute the MCMC kernel of Y.
#pragma omp for
      for (int j0 = 0; j0 < n0; j0++) {
#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        const mat Y_j_prev = Y.rows(j0_vector * n0 + j0);
        mat error_j        = error.rows(j0_vector * n0 + j0);

        // Make sure that the sum of the error terms is zero.
        const rowvec error_mean = sum(error_j, 0) / subspots;
        for (int r = 0; r < subspots; r++) {
          error_j.row(r) = error_j.row(r) - error_mean;
        }

        const mat Y_j_new = Y_j_prev + error_j;
        const mat mu_i_j  = mu_i_long.rows(j0_vector * n0 + j0);
        vec p_prev        = {0.0};
        vec p_new         = {0.0};
        for (int r = 0; r < subspots; r++) {
          p_prev += dmvnrm_prec_arma_fast(
                        Y_j_prev.row(r), mu_i_j.row(r),
                        lambda_i * w(j0 + n0 * r), true
                    ) -
                    c * (accu(pow(Y_j_prev.row(r) - Y0.row(j0), 2)));
          p_new +=
              dmvnrm_prec_arma_fast(
                  Y_j_new.row(r), mu_i_j.row(r), lambda_i * w(j0 + n0 * r), true
              ) -
              c * (accu(pow(Y_j_new.row(r) - Y0.row(j0), 2)));
        }
        double probY_j = as_scalar(exp(p_new - p_prev));
        if (probY_j > 1) {
          probY_j = 1;
        }
#pragma omp critical(iter_y)
        {
          acceptance_prob(j0)             = probY_j;
          Y_new.rows(j0_vector * n0 + j0) = Y_j_new;
        }
      }

#pragma omp single
      {
        int updateCounter = 0;

        // Accept or reject proposals of Y; update w; propose new values for
        // z.
        for (int j0 = 0; j0 < n0; j0++) {
          const IntegerVector Ysample = {0, 1};
          const NumericVector probsY  = {
              1 - acceptance_prob(j0), acceptance_prob(j0)
          };
          const int yesUpdate = sample(Ysample, 1, true, probsY)[0];
          if (yesUpdate == 1) {
            Y.rows(j0_vector * n0 + j0) = Y_new.rows(j0_vector * n0 + j0);
            updateCounter++;
          }

          for (int r = 0; r < subspots; r++) {
            // Update w.
            if (tdist) {
              const double w_beta = as_scalar(
                  2 /
                  ((Y.row(r * n0 + j0) - mu_i_long.row(r * n0 + j0)) *
                       lambda_i *
                       (Y.row(r * n0 + j0) - mu_i_long.row(r * n0 + j0)).t() +
                   4)
              );   // scale parameter
              w(r * n0 + j0) =
                  R::rgamma(w_alpha, w_beta);   // sample from posterior for w
            }

            // Propose new values for z.
            const IntegerVector qlessk = qvec[qvec != z(r * n0 + j0)];
            z_new(r * n0 + j0)         = sample(qlessk, 1)[0];
          }
        }
        Ychange[i] = updateCounter * 1.0 / n0;
      }

      // Update z
#pragma omp for
      for (int j = 0; j < n; j++) {
#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        const int z_j_prev      = z(j);
        const int z_j_new       = z_new(j);
        const Neighbor j_vector = neighbors[j];

        // log likelihood; prior
        vec h_z_prev(2, arma::fill::zeros), h_z_new(2, arma::fill::zeros);

        h_z_prev(0) = dmvnrm_prec_arma_fast(
            Y.row(j), mu_i.row(z_j_prev - 1), lambda_i * w(j), true
        )[0];
        h_z_new(0) = dmvnrm_prec_arma_fast(
            Y.row(j), mu_i.row(z_j_new - 1), lambda_i * w(j), true
        )[0];

        if (j_vector.get_size() != 0) {
          h_z_prev(1) = gamma / j_vector.get_size() * 2 *
                        accu((z(j_vector.get_neighbors()) == z_j_prev));
          h_z_new(1) = gamma / j_vector.get_size() * 2 *
                       accu((z(j_vector.get_neighbors()) == z_j_new));
        }
        double prob_j = exp(accu(h_z_new) - accu(h_z_prev));
        if (prob_j > 1) {
          prob_j = 1;
        }

#pragma omp critical(iter_z)
        {
          acceptance_prob(j)     = prob_j;
          log_likelihoods.row(j) = {h_z_prev(0), h_z_new(0)};
        }
      }

#pragma omp single
      {
        // Accept or reject proposals of z.
        for (int j = 0; j < n; j++) {
          const IntegerVector zsample = {0, 1};
          const NumericVector probs   = {
              1 - acceptance_prob(j), acceptance_prob(j)
          };
          const uword yesUpdate = sample(zsample, 1, true, probs)[0];
          log_likelihoods.set_col_idx(j, yesUpdate);
          if (yesUpdate == 1) {
            z(j) = z_new(j);
          }
        }
        plogLik[i] = accu(log_likelihoods.get_current_values());

        // Save samples for every 100 iterations.
        if ((i + 1) % 100 == 0) {
          df_sim_lambda[(i + 1) / 100] = lambda_i;
          df_sim_Y[(i + 1) / 100]      = Y;
          df_sim_w.row((i + 1) / 100)  = w.t();
          df_sim_z.row((i + 1) / 100)  = z.t();
        }
      }

      if ((i + 1) % 100 == 0 && early_stop > 0)
        i = nrep;
    }
  }

  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange,
      _["plogLik"] = plogLik
  );

  indicators::show_console_cursor(true);

#ifdef _OPENMP
  if (verbose) {
    print_thread_hits(thread_hits);
  }
#endif

  return (out);
}
