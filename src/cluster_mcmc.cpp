// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <testthat.h>
using namespace Rcpp;

struct ModelState {
  arma::mat mu;
  arma::mat lambda;
  arma::vec z;
};

// arma::mat update_mu(ModelState curr, ModelState prev) {
arma::mat update_mu(arma::mat Y, arma::rowvec z_prev, arma::mat lambda_prev, 
                    arma::mat lambda0, arma::vec mu0, int n_clusters, int n_dims) {
  
  int n_i;
  arma::uvec index_1k;
  arma::vec mean_i;
  arma::mat var_i;
  NumericVector Ysums;
    
  arma::mat mu_i(n_clusters, n_dims, arma::fill::zeros);
    
  for (int k = 1; k <= n_clusters; k++){
    index_1k = arma::find(z_prev == k);
    n_i = index_1k.n_elem;
    Ysums = sum(Y.rows(index_1k), 0);

    mean_i = inv(lambda0 + n_i * lambda_prev) * (lambda0 * mu0 + lambda_prev * as<arma::colvec>(Ysums));
    var_i = inv(lambda0 + n_i * lambda_prev);

    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
  }
    
  return mu_i;
}

arma::mat update_lambda(arma::mat Y, arma::mat mu_i, arma::rowvec z_prev, 
                        int n_spots, int n_dims, double alpha, double beta) {
  
  arma::mat mu_i_long(n_spots, n_dims, arma::fill::zeros);
  
  for (int j = 0; j < n_spots; j++){
    mu_i_long.row(j) = mu_i.row(z_prev[j] - 1);
  }
  
  arma::mat sumofsq = (Y - mu_i_long).t() * (Y - mu_i_long);
  
  arma::vec beta_d(n_dims); 
  beta_d.fill(beta);
  
  arma::mat lambda_i = rwish(n_spots + alpha, inv(arma::diagmat(beta_d) + sumofsq));

  return lambda_i;
}

arma::mat update_z(arma::mat Y, arma::rowvec z_prev, List df_j, IntegerVector qvec, 
                      arma::mat mu_i, arma::mat sigma_i,
                      int n_spots, double gamma) {
  
  arma::rowvec z = z_prev;
  arma::rowvec plogLikj(n_spots, arma::fill::zeros);
  double h_z_prev;
  double h_z_new;
  
  for (int j = 0; j < n_spots; j++){
    int z_j_prev = z(j);
    IntegerVector qlessk = qvec[qvec != z_j_prev];
    int z_j_new = sample(qlessk, 1)[0];
    
    arma::uvec j_vector = df_j[j];
    
    // has neighbors
    if (j_vector.size() != 0){
      h_z_prev = gamma/j_vector.size() * 2*arma::accu((z(j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
      h_z_new  = gamma/j_vector.size() * 2*arma::accu((z(j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
    } else {
      h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
      h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
    }
    
    double prob_j = exp(h_z_new - h_z_prev);
    if (prob_j > 1){
      prob_j = 1;
    }
    
    IntegerVector zsample = {z_j_prev, z_j_new};
    NumericVector probs = {1-prob_j, prob_j};
    
    z(j) = sample(zsample, 1, true, probs)[0];
    plogLikj(j) = h_z_prev;
  }
  
  arma::mat result = join_cols(z, plogLikj);
  return result;
}

// [[Rcpp::export]]
List cluster_mcmc(arma::mat Y, List df_j, int nrep, int n_spots, int n_dims, double gamma, 
                  int n_clusters, arma::vec init, NumericVector mu0, arma::mat lambda0, 
                  double alpha, double beta){
  
  //Initalize matrices storing iterations
  arma::mat df_sim_z(nrep, n_spots, arma::fill::zeros);
  arma::mat df_sim_mu(nrep, n_clusters * n_dims, arma::fill::zeros);
  List df_sim_lambda(nrep);
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  arma::rowvec initmu = rep(mu0, n_clusters);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  
  //Iterate
  arma::vec mu0vec = as<arma::vec>(mu0);
  arma::mat mu_i(n_clusters, n_dims, arma::fill::zeros);  // q x d
  arma::mat lambda_i;  // d x d
  arma::mat sigma_i;
  
  for (int i = 1; i < nrep; i++){
    //Update mu
    mu_i = update_mu(Y, df_sim_z.row(i - 1), df_sim_lambda[i - 1], lambda0, mu0vec, n_clusters, n_dims);
    df_sim_mu.row(i) = mu_i.as_row();
    
    //Update lambda
    lambda_i = update_lambda(Y, mu_i, df_sim_z.row(i - 1), n_spots, n_dims, alpha, beta);
    df_sim_lambda[i] = lambda_i;
    sigma_i = inv(lambda_i);
    
    //Update z
    IntegerVector qvec = seq_len(n_clusters);
    arma::mat result = update_z(Y, df_sim_z.row(i - 1), df_j, qvec, mu_i, sigma_i, n_spots, gamma);
    df_sim_z.row(i) = result.row(0);
    plogLik[i] = arma::sum(result.row(1));
  }
  
  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["plogLik"] = plogLik);
  
  return(out);
}

context("Sample unit tests v3") {
  
  arma::mat lambda0 = arma::mat(5, 5, arma::fill::eye);
  
  test_that("Covariance is symmetric") {
    expect_true(lambda0.is_symmetric());
  }
}