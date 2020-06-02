#include "PottsModel.h"

#include <RcppDist.h>
#include "ClusterParams.h"

using namespace Rcpp;

// Compute next mean vector
arma::mat updateMu(ClusterParams params, PottsModel prev) {
  arma::mat mu_i(params.q, params.n_dims(), arma::fill::zeros);
  
  for (int k = 1; k <= params.q; k++){
    arma::uvec cluster_members = arma::find(prev.z == k);
    
    // Cov  = (Λ_0 + n_k * Λ)^{-1}
    arma::mat var_i = arma::inv(params.lambda0 + cluster_members.n_elem * prev.lambda);
    
    // Mean = (Λ_0 + n_k * Λ)^{-1} (Λ_0 * μ_0 + Λ * Σ(y_i))
    arma::vec Ysums = arma::sum(params.Y.rows(cluster_members), 0).t();
    arma::vec mean_i = var_i * (params.lambda0 * params.mu0 + prev.lambda * Ysums);
    
    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
  }
  
  return mu_i;
}

// Compute next covariance matrix
arma::mat updateLambda(ClusterParams params, PottsModel curr, PottsModel prev) {
  
  arma::mat mu_i_long(params.n_spots(), params.n_dims(), arma::fill::zeros);
  
  for (int j = 0; j < params.n_spots(); j++){
    mu_i_long.row(j) = curr.mu.row(prev.z[j] - 1);
  }
  
  arma::mat sumofsq = (params.Y - mu_i_long).t() * (params.Y - mu_i_long);
  
  arma::vec beta_d(params.n_dims()); 
  beta_d.fill(params.beta);
  
  arma::mat lambda_i = rwish(params.n_spots() + params.alpha, inv(arma::diagmat(beta_d) + sumofsq));
  
  return lambda_i;
}


// TODO store neighbors in params
double computeEnergy(arma::rowvec z, arma::uvec neighbors, 
                     int label, ClusterParams params) {
  
  double H;
  H = (params.gamma / neighbors.n_elem) * 2 * arma::accu(z(neighbors) == label);
  
  return H;
}

double computeLogLikelihood(arma::rowvec x, int label, PottsModel state) {
  double logP_YGivenZ;
  
  // Mean vector of given cluster
  arma::vec mean = arma::vectorise(state.mu.row(label - 1));
  
  logP_YGivenZ = dmvnorm(x, mean, state.sigma(), true)[0];
  return logP_YGivenZ;
}

// Compute next set of cluster assignments and corresponding log-likelihoods
// TODO: extract energy/likelihood computation
arma::mat updateZ(ClusterParams params, PottsModel curr, PottsModel prev, List df_j) {
  
  arma::rowvec z = prev.z;
  arma::rowvec plogLikj(params.n_spots(), arma::fill::zeros);
  IntegerVector qvec = seq_len(params.q);
  double h_z_prev;
  double h_z_new;
  
  for (int j = 0; j < params.n_spots(); j++){
    int z_j_prev = z(j);
    IntegerVector qlessk = qvec[qvec != z_j_prev];
    int z_j_new = sample(qlessk, 1)[0];
    
    // Locations of neighbors
    arma::uvec neighbors = df_j[j];
    
    // has neighbors
    if (neighbors.n_elem > 0) {
      // TODO: move this into PottsModel::computeLogPosterior() and PottsModel::computeAcceptance()
      h_z_prev = computeEnergy(z, neighbors, z_j_prev, params) + computeLogLikelihood(params.Y.row(j), z_j_prev, curr);
      h_z_new  = computeEnergy(z, neighbors, z_j_new, params) + computeLogLikelihood(params.Y.row(j), z_j_new, curr);
    } else {
      h_z_prev = computeLogLikelihood(params.Y.row(j), z_j_prev, curr);
      h_z_new  = computeLogLikelihood(params.Y.row(j), z_j_new, curr);
    }
    
    double prob_j = std::min(exp(h_z_new - h_z_prev), 1.0);
    
    IntegerVector zsample = {z_j_prev, z_j_new};
    NumericVector probs = {1-prob_j, prob_j};
    
    z(j) = sample(zsample, 1, true, probs)[0];
    
    // This is logPosterior
    plogLikj(j) = h_z_prev;
  }
  
  arma::mat result = join_cols(z, plogLikj);
  
  return result;
}