// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;

arma::mat update_mu(const arma::mat& Y, int q, const arma::rowvec& z_prev, 
                    const arma::vec& w, const arma::mat& lambda0, 
                    const arma::mat& lambda_prev, const arma::colvec& mu0vec) {

  mat mu_i(q, Y.n_cols);
  NumericVector Ysums;

  for (int k = 1; k <= q; k++){
    uvec index_1k = find(z_prev == k);
    int n_i = sum(w(index_1k)); 
    mat Yrows = Y.rows(index_1k);
    Yrows.each_col() %= w(index_1k);
    Ysums = sum(Yrows, 0);
    vec mean_i = inv(lambda0 + n_i * lambda_prev) * (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums));  //posterior mean
    mat var_i = inv(lambda0 + n_i * lambda_prev);  //posterior variance
    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);  //sample from posterior for mu
  }

  return mu_i;
}

arma::mat update_lambda(const arma::mat& resid, const arma::rowvec& z_prev,
                        const arma::vec& w, double alpha, double beta) {

  int d = resid.n_cols;
  int n = resid.n_rows;

  mat sumofsq = resid.t() * diagmat(w) * resid;
  vec beta_d(d); 
  beta_d.fill(beta);
  mat Vinv = diagmat(beta_d);
  mat lambda_i = rwish(n + alpha, inv(Vinv + sumofsq));

  return lambda_i;

}

// Per spot update
// Can't extract from loop b/c will change output relative to original
// TODO: replace magic 4 with defined constant
// TODO rm int j and pass row directly
double update_w_j(const arma::mat& resid, const arma::mat& lambda_i, int j) {
  int d = resid.n_cols;

  // shape parameter
  double w_alpha = (d + 4) / 2;  

  // scale parameter
  double w_beta = as_scalar(2/(resid.row(j) * lambda_i * resid.row(j).t() + 4));

  // sample from posterior for w
  double w_j = R::rgamma(w_alpha, w_beta);  

  return w_j;
}

double compute_h_z() {

}

// Spatial smoothing prior
double compute_pi_z() {
  
}

int sample_z_new(const arma::vec& qvec, int z_j_prev) {
  vec qlessk = qvec(find(qvec != z_j_prev));
  int z_j_new = Rcpp::RcppArmadillo::sample(qlessk, 1, false)[0];
  
  return z_j_new;
}

// int update_z_j(const arma::mat& Y, const arma::mat& mu_i, List df_j,
//                const arma::mat& df_sim_z, const arma::mat& sigma_i,
//                const arma::vec& w, const arma::vec& plogLik,
//                int q, int i, double gamma) {
// 
//   int n = Y.n_rows;
//   NumericVector plogLikj(n, NA_REAL);
// 
//   for (int j = 0; j < n; j++){
//     uvec j_vector = df_j[j];
//     uvec i_vector(1); i_vector.fill(i);
//     double h_z_prev;
//     double h_z_new;
// 
//     // Has neighbors
//     if (j_vector.size() != 0){
//       h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i/w[j], true)[0];
//       h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i/w[j], true)[0];
//     } else {
//       h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i/w[j], true)[0];
//       h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i/w[j], true)[0];
//     }
// 
//     double prob_j = exp(h_z_new-h_z_prev);
//     if (prob_j > 1){
//       prob_j = 1;
//     }
// 
//     IntegerVector zsample = {z_j_prev, z_j_new};
//     NumericVector probs = {1 - prob_j, prob_j};
//     df_sim_z(i,j) = sample(zsample, 1, true, probs)[0];
//     plogLikj[j] = h_z_prev;
//   }
// 
//   plogLik[i] = sum(plogLikj);
// }

// [[Rcpp::export]]
List iterate_t_refactor(arma::mat Y, List df_j, int nrep, int n, int d, double gamma, int q, arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha, double beta){
  
  //Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep);
  mat df_sim_w(nrep, n);
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  vec w = ones<vec>(n);
  df_sim_w.row(0) = w.t();
  
  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++){
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();
    
    //Update mu
    mat mu_i = update_mu(Y, q, df_sim_z.row(i-1), w, lambda0, df_sim_lambda[i-1], mu0vec);
    df_sim_mu.row(i) = vectorise(mu_i, 1);
    
    // Compute residuals
    mat mu_i_long(n,d);
    for (int j = 0; j < n; j++){
      mu_i_long.row(j) = mu_i.row(df_sim_z(i-1, j) - 1);
    }
    mat resid = Y - mu_i_long;

    //Update lambda
    mat lambda_i = update_lambda(resid, df_sim_z.row(i-1), w, alpha, beta);
    df_sim_lambda[i] = lambda_i;
    mat sigma_i = inv(lambda_i);
    
    //Update z and w
    df_sim_z.row(i) = df_sim_z.row(i-1);
    const vec qvec = linspace(1, q, q);
    NumericVector plogLikj(n, NA_REAL);

    for (int j = 0; j < n; j++){
      w[j] = update_w_j(resid, lambda_i, j);
      
      int z_j_prev = df_sim_z(i,j);
      int z_j_new = sample_z_new(qvec, z_j_prev);
      
      uvec j_vector = df_j[j];
      uvec i_vector(1); i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i/w[j], true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i/w[j], true)[0];
      } else {
        h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i/w[j], true)[0];
        h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i/w[j], true)[0];
      }
      double prob_j = exp(h_z_new-h_z_prev);
      if (prob_j > 1){
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs = {1-prob_j, prob_j};
      df_sim_z(i,j) = sample(zsample, 1, true, probs)[0];
      plogLikj[j] = h_z_prev;
    }
    df_sim_w.row(i) = w.t();
    plogLik[i] = sum(plogLikj);
  }

  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["weights"] = df_sim_w, _["plogLik"] = plogLik);

  return(out);
}
