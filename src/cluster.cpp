// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include "double_states_vector.h"
#include "neighbor.h"
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <algorithm>
#include <chrono>
#include <cmath>
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

mat
adaptive_mcmc(
    double accpetance_rate, double target_acceptance_rate, size_t curr_iter,
    const mat &identity_mtx, const mat &adaptive_mtx, const rowvec &samples
) {
  const double step_size =
      std::min(1.0, samples.n_elem * std::pow(curr_iter, -2.0 / 3));

  return chol(
      adaptive_mtx *
          (identity_mtx + step_size *
                              (accpetance_rate - target_acceptance_rate) *
                              samples.t() * samples / accu(samples % samples)) *
          adaptive_mtx.t(),
      "lower"
  );
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
convert_neighbors(const CharacterVector &ori, std::vector<Neighbor> &des) {
  for (auto i = 0; i != ori.size(); i++) {
    if (CharacterVector::is_na(ori(i))) {
      des.emplace_back(Neighbor());
    } else {
      des.emplace_back(
          Neighbor(static_cast<std::string>(ori(i)), std::string(","), false)
      );
    }
  }
}

umat
convert_neighbors2mtx(std::vector<Neighbor> neighbors) {
  umat ret(neighbors.size(), 4, fill::zeros);

  for (size_t i = 0; i < neighbors.size(); i++) {
    if (neighbors[i].get_size() > 0) {
      ret(i, span(0, neighbors[i].get_size() - 1)) =
          neighbors[i].get_neighbors().t() + 1;
    }
  }

  return ret;
}

void
find_subspot_neighbors(
    const std::vector<Neighbor> &spot_neighbors, const mat &subspot_positions,
    const double dist, std::vector<Neighbor> &subspot_neighbors
) {
  const uword num_subspots = subspot_positions.n_rows / spot_neighbors.size();

  if (std::abs(
          static_cast<double>(subspot_positions.n_rows) /
              spot_neighbors.size() -
          num_subspots
      ) > 1e-6) {
    throw std::runtime_error("Invalid arguments!");
  }

  const uvec subspots = linspace<uvec>(0, num_subspots - 1, num_subspots);

  for (uword subspot_idx = 0; subspot_idx < subspot_positions.n_rows;
       subspot_idx++) {
    // Get spot index.
    const uword spot_idx = subspot_idx % spot_neighbors.size();

    // Get candidate spot neighbors.
    urowvec __spot_neighbors(spot_neighbors[spot_idx].get_neighbors().t());
    __spot_neighbors.resize(__spot_neighbors.n_elem + 1);
    __spot_neighbors(__spot_neighbors.n_elem - 1) = spot_idx;

    // Get candidate subspot neighbors.
    const uvec candidate_neighbors =
        ((umat(subspots) *
          umat(ones<urowvec>(__spot_neighbors.n_elem) * spot_neighbors.size()))
             .eval()
             .each_row() +
         __spot_neighbors)
            .as_col();

    // Compute Manhattan distance.
    const mat man_dist =
        sum(abs(subspot_positions.rows(candidate_neighbors).eval().each_row() -
                subspot_positions.rows(uvec(std::vector<arma::uword>{subspot_idx
                }))),
            1);

    // Get neighbors.
    const uvec neighbor_idx = (man_dist < dist) && (man_dist > 1e-6);

    if (accu(neighbor_idx) < 2) {
      throw std::runtime_error("Error in finding neighbors of subspots!");
    }

    subspot_neighbors.emplace_back(
        Neighbor(candidate_neighbors(find(neighbor_idx == 1)))
    );
  }
}

// [[Rcpp::export]]
List
iterate(
    const arma::mat &Y, const List &df_j, int nrep, int thin, int n, int d,
    double gamma, int q, const arma::uvec &init, const NumericVector &mu0,
    const arma::mat &lambda0, double alpha, double beta
) {

  // Initalize matrices storing iterations
  umat df_sim_z(nrep / thin + 1, n, fill::zeros);
  mat df_sim_mu(nrep / thin + 1, q * d, fill::zeros);
  List df_sim_lambda(nrep / thin + 1);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  mat lambda_i     = lambda0;
  df_sim_lambda[0] = lambda0;
  uvec z           = init;
  df_sim_z.row(0)  = init.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      int n_i       = index_1k.n_elem;
      NumericVector Ysums;
      mat Yrows  = Y.rows(index_1k);
      Ysums      = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_i) *
                   (lambda0 * mu0vec + lambda_i * as<colvec>(Ysums));
      mat var_i       = inv(lambda0 + n_i * lambda_i);
      mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
    }

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(z(j) - 1);
    }
    mat sumofsq = (Y - mu_i_long).t() * (Y - mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv          = diagmat(beta_d);
    lambda_i          = rwish(n + alpha, inv(Vinv + sumofsq));
    const mat sigma_i = inv(lambda_i);

    // Update z
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev         = z(j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev =
            gamma / j_vector.size() * 2 * accu((z(j_vector) == z_j_prev)) +
            dmvnrm_arma_fast(
                Y.row(j), mu_i.row(z_j_prev - 1), sigma_i, true
            )[0];
        h_z_new =
            gamma / j_vector.size() * 2 * accu((z(j_vector) == z_j_new)) +
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
      z(j)                  = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);

    // Save samples for every a few iterations.
    if ((i + 1) % thin == 0) {
      df_sim_mu.row((i + 1) / thin) = vectorise(mu_i, 1);
      df_sim_lambda[(i + 1) / thin] = lambda_i;
      df_sim_z.row((i + 1) / thin)  = z.t();
    }
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
    const arma::mat &Y, const List &df_j, int nrep, int thin, int n, int d,
    double gamma, int q, const arma::uvec &init, const NumericVector &mu0,
    const arma::mat &lambda0, double alpha, double beta
) {

  // Initalize matrices storing iterations
  umat df_sim_z(nrep / thin + 1, n, fill::zeros);
  mat df_sim_mu(nrep / thin + 1, q * d, fill::zeros);
  List df_sim_lambda(nrep / thin + 1);
  List lambda_list(q);
  lambda_list[0] = lambda0;
  List sigma_list(q);
  sigma_list[0] = inv(lambda0);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  for (int k = 1; k < q; k++) {
    lambda_list[k] = lambda_list[0];
    sigma_list[k]  = sigma_list[0];
  }
  df_sim_lambda[0] = lambda_list;
  uvec z           = init;
  df_sim_z.row(0)  = init.t();

  // Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++) {
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Update mu
    mat mu_i(q, d);
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      int n_i       = index_1k.n_elem;
      mat Yrows     = Y.rows(index_1k);
      Ysums         = sum(Yrows, 0);
      mat lambda_k  = lambda_list[k - 1];
      vec mean_i    = inv(lambda0 + n_i * lambda_k) *
                   (lambda0 * mu0vec + lambda_k * as<colvec>(Ysums));
      mat var_i       = inv(lambda0 + n_i * lambda_k);
      mu_i.row(k - 1) = rmvnorm(1, mean_i, var_i);
    }

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(z(j) - 1);
    }
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      int n_i       = index_1k.n_elem;
      mat sumofsq   = (Y.rows(index_1k) - mu_i_long.rows(index_1k)).t() *
                    (Y.rows(index_1k) - mu_i_long.rows(index_1k));
      mat lambda_i       = rwish(n_i + alpha, inv(Vinv + sumofsq));
      lambda_list[k - 1] = lambda_i;
      sigma_list[k - 1]  = inv(lambda_i);
    }

    // Update z
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev         = z(j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new          = sample(qlessk, 1)[0];
      uvec j_vector        = df_j[j];
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0) {
        h_z_prev =
            gamma / j_vector.size() * 2 * accu(((j_vector) == z_j_prev)) +
            dmvnrm_arma_fast(
                Y.row(j), mu_i.row(z_j_prev - 1), sigma_list[z_j_prev - 1], true
            )[0];
        h_z_new =
            gamma / j_vector.size() * 2 * accu((z(j_vector) == z_j_new)) +
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
      z(j)                  = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);

    // Save samples for every a few iterations.
    if ((i + 1) % thin == 0) {
      df_sim_mu.row((i + 1) / thin) = vectorise(mu_i, 1);
      df_sim_lambda[(i + 1) / thin] = lambda_list;
      df_sim_z.row((i + 1) / thin)  = z.t();
    }
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
    const arma::mat &Y, const List &df_j, int nrep, int thin, int n, int d,
    double gamma, int q, const arma::uvec &init, const NumericVector &mu0,
    const arma::mat &lambda0, double alpha, double beta
) {

  // Initalize matrices storing iterations
  umat df_sim_z(nrep / thin + 1, n, fill::zeros);
  mat df_sim_mu(nrep / thin + 1, q * d, fill::zeros);
  List df_sim_lambda(nrep / thin + 1);
  mat df_sim_w(nrep / thin + 1, n);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  mat lambda_i     = lambda0;
  df_sim_lambda[0] = lambda0;
  uvec z           = init;
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
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      double n_i    = sum(w(index_1k));
      mat Yrows     = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums = sum(Yrows, 0);
      vec mean_i =
          inv(lambda0 + n_i * lambda_i) *
          (lambda0 * mu0vec + lambda_i * as<colvec>(Ysums));   // posterior mean
      mat var_i = inv(lambda0 + n_i * lambda_i);   // posterior variance
      mu_i.row(k - 1) =
          rmvnorm(1, mean_i, var_i);   // sample from posterior for mu
    }

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(z(j) - 1);
    }
    mat sumofsq = (Y - mu_i_long).t() * diagmat(w) * (Y - mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv          = diagmat(beta_d);
    lambda_i          = rwish(n + alpha, inv(Vinv + sumofsq));
    const mat sigma_i = inv(lambda_i);

    // Update z and w
    double w_alpha = (d + 4) / 2;   // shape parameter
    double w_beta;
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      w_beta = as_scalar(
          2 / ((Y.row(j) - mu_i_long.row(j)) * lambda_i *
                   (Y.row(j) - mu_i_long.row(j)).t() +
               4)
      );                                   // scale parameter
      w[j] = R::rgamma(w_alpha, w_beta);   // sample from posterior for w

      int z_j_prev         = z(j);
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
      z(j)                  = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);

    // Save samples for every a few iterations.
    if ((i + 1) % thin == 0) {
      df_sim_mu.row((i + 1) / thin) = vectorise(mu_i, 1);
      df_sim_lambda[(i + 1) / thin] = lambda_i;
      df_sim_w.row((i + 1) / thin)  = w.t();
      df_sim_z.row((i + 1) / thin)  = z.t();
    }
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
    const arma::mat &Y, const List &df_j, int nrep, int thin, int n, int d,
    double gamma, int q, const arma::uvec &init, const NumericVector &mu0,
    const arma::mat &lambda0, double alpha, double beta
) {

  // Initalize matrices storing iterations
  umat df_sim_z(nrep / thin + 1, n, fill::zeros);
  mat df_sim_mu(nrep / thin + 1, q * d, fill::zeros);
  List df_sim_lambda(nrep / thin + 1);
  mat df_sim_w(nrep / thin + 1, n);
  List lambda_list(q);
  lambda_list[0] = lambda0;
  List sigma_list(q);
  sigma_list[0] = inv(lambda0);
  NumericVector plogLik(nrep, NA_REAL);

  // Initialize parameters
  rowvec initmu    = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  for (int k = 1; k < q; k++) {
    lambda_list[k] = lambda_list[0];
    sigma_list[k]  = sigma_list[0];
  }
  df_sim_lambda[0] = lambda_list;
  uvec z           = init;
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
    NumericVector Ysums;
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      double n_i    = sum(w(index_1k));
      mat Yrows     = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums        = sum(Yrows, 0);
      mat lambda_k = lambda_list[k - 1];
      vec mean_i =
          inv(lambda0 + n_i * lambda_k) *
          (lambda0 * mu0vec + lambda_k * as<colvec>(Ysums));   // posterior mean
      mat var_i = inv(lambda0 + n_i * lambda_k);   // posterior variance
      mu_i.row(k - 1) =
          rmvnorm(1, mean_i, var_i);   // sample from posterior for mu
    }

    // Update lambda
    mat mu_i_long(n, d);
    for (int j = 0; j < n; j++) {
      mu_i_long.row(j) = mu_i.row(z(j) - 1);
    }
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    for (int k = 1; k <= q; k++) {
      uvec index_1k = find(z == k);
      int n_i       = index_1k.n_elem;
      mat sumofsq   = (Y.rows(index_1k) - mu_i_long.rows(index_1k)).t() *
                    diagmat(w(index_1k)) *
                    (Y.rows(index_1k) - mu_i_long.rows(index_1k));
      mat lambda_i       = rwish(n_i + alpha, inv(Vinv + sumofsq));
      lambda_list[k - 1] = lambda_i;
      sigma_list[k - 1]  = inv(lambda_i);
    }

    // Update z and w
    double w_alpha = (d + 4) / 2;   // shape parameter
    double w_beta;
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++) {
      int z_j_prev = z(j);
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
      z(j)                  = sample(zsample, 1, true, probs)[0];
      plogLikj[j]           = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);

    // Save samples for every a few iterations.
    if ((i + 1) % thin == 0) {
      df_sim_mu.row((i + 1) / thin) = vectorise(mu_i, 1);
      df_sim_lambda[(i + 1) / thin] = lambda_list;
      df_sim_w.row((i + 1) / thin)  = w.t();
      df_sim_z.row((i + 1) / thin)  = z.t();
    }
  }
  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = df_sim_w, _["plogLik"] = plogLik
  );
  return (out);
}

/**
 * @brief
 *
 * @param subspot_positions: coordinates of subspots
 * @param dist: the maximum distance used to identify neighbors of each subspot
 * @param spot_neighbors: indices of neighbors of each spot (1-based)
 * @param Y: the initialized principal components (num_subspots * num_pcs)
 * @param tdist: whether to use multivariate t distribution or not
 * @param nrep: the number of MCMC iterations
 * @param n: the number of subspots (after deconvolution)
 * @param n0: the number of spots (before deconvolution)
 * @param d: the number of PCs
 * @param gamma: smoothing parameter
 * @param q: the number of clusters
 * @param init: the initialized clustering of subspots
 * @param subspots: the number of subspots of each spot
 * @param verbose: whether to print more information
 * @param jitter_scale: the amount of jittering (the variance) for the proposal
 * distribution
 * @param c: the amount of jittering (the variance) for the prior
 * distribution
 * @param mu0: the mean hyperparameetr of mu
 * @param lambda0: the precision hyperparameter of mu
 * @param alpha: one of the hyperparamters of lambda
 * @param beta: one of the hyperparamters of lambda
 * @param thread_num: the number of threads to be used
 * @param verbose: information for debugging
 *
 * @return List MCMC samples of latent variables in a list
 */
// [[Rcpp::export]]
List
iterate_deconv(
    const arma::mat &subspot_positions, const double dist,
    const CharacterVector &spot_neighbors, arma::mat &Y, bool tdist, int nrep,
    int thin, int n, int n0, int d, double gamma, int q, const arma::uvec &init,
    int subspots, bool verbose, double jitter_scale, int adapt_before, double c,
    const NumericVector &mu0, const arma::mat &lambda0, double alpha,
    double beta, int thread_num = 1
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

  // To identify neighbors of subspots.
  if (verbose) {
    std::cout << "[DEBUG] Identifying neighbors of subspots..." << std::endl;
  }
  std::vector<Neighbor> __spot_neighbors, __subspot_neighbors;
  convert_neighbors(spot_neighbors, __spot_neighbors);
  find_subspot_neighbors(
      __spot_neighbors, subspot_positions, dist, __subspot_neighbors
  );

  // Initalize matrices storing iterations
  const mat Y0        = Y.rows(0, n0 - 1);   // The input PCs on spot level.
  mat Y_new           = mat(Y.n_rows,
                            Y.n_cols);   // The proposed PCs on subspot level.
  vec acceptance_prob = vec(n0);   // The probability of accepting the proposals
                                   // on subspot level.
  DoubleStatesVector<double> log_likelihoods(n
  );   // The log-likelihoods on subspot level.
  umat df_sim_z(nrep / thin + 1, n, fill::zeros);
  mat df_sim_mu(nrep / thin + 1, d * q, fill::zeros);
  List df_sim_lambda(nrep / thin + 1);
  List df_sim_Y(nrep / thin + 1);
  mat df_sim_w(nrep / thin + 1, n);
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

  // Iterate
  const colvec mu0vec = as<colvec>(mu0);
  mat mu_i(q, d);
  mat mu_i_long(n, d);
  const uvec j0_vector = linspace<uvec>(0, subspots - 1, subspots);
  mat error(n, d);
  const double w_alpha     = (d + 4) / 2;   // shape parameter
  const IntegerVector qvec = seq_len(q);
  const vec zero_vec       = zeros<vec>(d);
  const vec one_vec        = ones<vec>(d);
  const mat error_var      = diagmat(one_vec);
  jitter_scale /= d;

  // For adaptive MCMC
  std::vector<size_t> num_accepts(n0, 0);
  std::vector<size_t> num_rejects(n0, 0);
  std::vector<mat> adaptive_mtx(n);
  if (jitter_scale == 0.0) {
    if (verbose) {
      std::cout << "[DEBUG] Turning on adaptive MCMC ";

      if (adapt_before == 0) {
        std::cout << "throughout the entire chain." << std::endl;
      } else {
        std::cout << "only in the first " << adapt_before << " iterations."
                  << std::endl;
      }
    }

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
      adaptive_mtx[i] = diagmat(one_vec);
    }
  }

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
  signal(SIGTERM, sig_handler);

  // Time measurement
  // double t_start, t_end;

#pragma omp parallel shared(                                                   \
        early_stop, Y, Y_new, error, error_var, num_accepts, num_rejects,      \
            adaptive_mtx, acceptance_prob, j0_vector, subspots, zero_vec,      \
            mu_i_long, Y0, lambda_i, n0, c, thread_hits, n, z,                 \
            __subspot_neighbors, gamma, mu_i, w, log_likelihoods               \
)
  {
    for (int i = 1; i < nrep; i++) {
#pragma omp single
      {

        // #ifdef _OPENMP
        // if (i == 1)
        // t_start = omp_get_wtime();
        // #endif

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
        // df_sim_mu.row(i) = vectorise(mu_i, 1);

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
            n, zero_vec,
            jitter_scale == 0.0 ? error_var : error_var * jitter_scale
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

        if (jitter_scale == 0.0) {
          for (int r = 0; r < subspots; r++) {
            error_j.row(r) = trans(
                adaptive_mtx[r * n0 + j0] *
                resize(error_j.row(r), error_j.n_cols, 1)
            );
          }
        }

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

        // Accept or reject proposals of Y; update w and z.
        for (int j0 = 0; j0 < n0; j0++) {
          const IntegerVector Ysample = {0, 1};
          const NumericVector probsY  = {
              1 - acceptance_prob(j0), acceptance_prob(j0)
          };

          uword yesUpdate = sample(Ysample, 1, true, probsY)[0];
          if (yesUpdate == 1) {
            Y.rows(j0_vector * n0 + j0) = Y_new.rows(j0_vector * n0 + j0);

            num_accepts[j0]++;
            updateCounter++;
          } else
            num_rejects[j0]++;

          for (int r = 0; r < subspots; r++) {

            // Adaptive MCMC.
            if (jitter_scale == 0.0 && i > 10 &&
                (adapt_before == 0 || i <= adapt_before)) {
              adaptive_mtx[r * n0 + j0] = adaptive_mcmc(
                  static_cast<double>(num_accepts[j0]) /
                      (num_accepts[j0] + num_rejects[j0]),
                  0.234, i, error_var, adaptive_mtx[r * n0 + j0],
                  error.row(r * n0 + j0)
              );
            }

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

            // Update z.
            const IntegerVector qlessk = qvec[qvec != z(r * n0 + j0)];
            const int z_j_new          = sample(qlessk, 1)[0];

            const int z_j_prev      = z(r * n0 + j0);
            const Neighbor j_vector = __subspot_neighbors[r * n0 + j0];

            // log likelihood; prior
            vec h_z_prev(2, arma::fill::zeros), h_z_new(2, arma::fill::zeros);

            h_z_prev(0) = dmvnrm_prec_arma_fast(
                Y.row(r * n0 + j0), mu_i.row(z_j_prev - 1),
                lambda_i * w(r * n0 + j0), true
            )[0];
            h_z_new(0) = dmvnrm_prec_arma_fast(
                Y.row(r * n0 + j0), mu_i.row(z_j_new - 1),
                lambda_i * w(r * n0 + j0), true
            )[0];

            if (j_vector.get_size() != 0) {
              h_z_prev(1) = gamma / j_vector.get_size() * 2 *
                            accu((z(j_vector.get_neighbors()) == z_j_prev));
              h_z_new(1) = gamma / j_vector.get_size() * 2 *
                           accu((z(j_vector.get_neighbors()) == z_j_new));
            }
            log_likelihoods.row(r * n0 + j0) = {h_z_prev(0), h_z_new(0)};

            const double prob_j = exp(accu(h_z_new) - accu(h_z_prev));

            if (prob_j >= 1) {
              yesUpdate = 1;
            } else {
              const IntegerVector zsample = {0, 1};
              const NumericVector probs   = {1 - prob_j, prob_j};

              yesUpdate = sample(zsample, 1, true, probs)[0];
            }

            log_likelihoods.set_col_idx(r * n0 + j0, yesUpdate);
            if (yesUpdate == 1) {
              z(r * n0 + j0) = z_j_new;
            }
          }
        }

        Ychange[i] = updateCounter * 1.0 / n0;
        plogLik[i] = accu(log_likelihoods.get_current_values());
      }

      // Adaptive MCMC.
      // #pragma omp for
      //       for (int j = 0; j < n; j++) {
      // #ifdef _OPENMP
      // #pragma omp atomic update
      //         thread_hits[omp_get_thread_num()]++;
      // #endif

      //         if (jitter_scale == 0.0 && i > 10 &&
      //             (adapt_before == 0 || i <= adapt_before)) {
      //           adaptive_mtx[j] = adaptive_mcmc(
      //               static_cast<double>(num_accepts[j % n0]) /
      //                   (num_accepts[j % n0] + num_rejects[j % n0]),
      //               0.234, i, error_var, adaptive_mtx[j], error.row(j)
      //           );
      //         }
      //       }

#pragma omp single
      {
        // Save samples for every a few iterations.
        if ((i + 1) % thin == 0) {
          df_sim_mu.row((i + 1) / thin) = vectorise(mu_i, 1);
          df_sim_lambda[(i + 1) / thin] = lambda_i;
          df_sim_Y[(i + 1) / thin]      = Y;
          df_sim_w.row((i + 1) / thin)  = w.t();
          df_sim_z.row((i + 1) / thin)  = z.t();
        }

        // #ifdef _OPENMP
        //         if (i == 1) {
        //           t_end = omp_get_wtime();

        //           if (verbose) {
        //             const double t_per_iter = t_end - t_start;

        //             std::cout << "[DEBUG] " << std::setprecision(2) <<
        //             t_per_iter
        //                       << "s per iteration (expecting "
        //                       << t_per_iter * nrep / 3600 << " hours in
        //                       total)."
        //                       << std::endl;
        //           }
        //         }
        // #endif
      }

      if ((i + 1) % thin == 0 && early_stop > 0)
        i = nrep;

#pragma omp barrier
    }
  }

  // #ifdef _OPENMP
  //   t_end = omp_get_wtime();

  //   if (verbose) {
  //     const double t_all = t_end - t_start;

  //     std::cout << "[DEBUG] Finished in " << std::setprecision(2) << t_all /
  //     3600
  //               << " hours." << std::endl;
  //   }
  // #endif

  List out = List::create(
      _["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda,
      _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange,
      _["plogLik"] = plogLik,
      _["df_j"]    = convert_neighbors2mtx(__subspot_neighbors)
  );

  indicators::show_console_cursor(true);

#ifdef _OPENMP
  if (verbose) {
    print_thread_hits(thread_hits);
  }
#endif

  return (out);
}
