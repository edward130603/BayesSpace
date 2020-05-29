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
                    arma::mat lambda0, arma::vec mu0, int q, int d) {
  
  int n_i;
  arma::uvec index_1k;
  arma::vec mean_i;
  arma::mat var_i;
  NumericVector Ysums;
    
  arma::mat mu_i(q, d, arma::fill::zeros);
    
  for (int k = 1; k <= q; k++){
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
                        int n, int d, double alpha, double beta) {
  
  arma::mat mu_i_long(n, d, arma::fill::zeros);
  
  for (int j = 0; j < n; j++){
    mu_i_long.row(j) = mu_i.row(z_prev[j] - 1);
  }
  
  arma::mat sumofsq = (Y - mu_i_long).t() * (Y - mu_i_long);
  
  arma::vec beta_d(d); 
  beta_d.fill(beta);
  
  arma::mat lambda_i = rwish(n + alpha, inv(arma::diagmat(beta_d) + sumofsq));

  return lambda_i;
}


// [[Rcpp::export]]
List cluster_mcmc(arma::mat Y, List df_j, int nrep, int n, int d, double gamma, 
                  int q, arma::vec init, NumericVector mu0, arma::mat lambda0, 
                  double alpha, double beta){
  
  //Initalize matrices storing iterations
  arma::mat df_sim_z(nrep, n, arma::fill::zeros);
  arma::mat df_sim_mu(nrep, q*d, arma::fill::zeros);
  List df_sim_lambda(nrep);
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  arma::rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  
  //Iterate
  arma::vec mu0vec = as<arma::vec>(mu0);
  arma::mat mu_i(q , d, arma::fill::zeros);
  arma::mat lambda_i(n, n, arma::fill::zeros);
  arma::mat sigma_i(n, n, arma::fill::zeros);
  for (int i = 1; i < nrep; i++){
    //Update mu
    mu_i = update_mu(Y, df_sim_z.row(i - 1), df_sim_lambda[i - 1], lambda0, mu0vec, q, d);
    df_sim_mu.row(i) = mu_i.as_row();
    
    //Update lambda
    lambda_i = update_lambda(Y, mu_i, df_sim_z.row(i - 1), n, d, alpha, beta);
    df_sim_lambda[i] = lambda_i;
    sigma_i = inv(lambda_i);
    
    //Update z
    df_sim_z.row(i) = df_sim_z.row(i-1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++){
      int z_j_prev = df_sim_z(i,j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      arma::uvec j_vector = df_j[j];
      arma::uvec i_vector(1); i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
      } else {
        h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
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
    plogLik[i] = sum(plogLikj);
  }
  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["plogLik"] = plogLik);
  
  return(out);
}