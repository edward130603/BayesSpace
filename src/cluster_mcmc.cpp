// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <testthat.h>
using namespace Rcpp;

// Hyperparameter defaults
const double kDefaultAlpha = 1;
const double kDefaultBeta = 0.01;
const double kDefaultGamma = 2;
const double kDefaultLambda = 0.01;

class PottsState {
public:
  arma::mat mu;
  arma::mat lambda;
  arma::rowvec z;
  // TODO arma::vec pLogLik
  
  PottsState(): mu(arma::zeros(0, 0)), lambda(arma::zeros(0, 0)), z(arma::zeros<arma::rowvec>(0)) {}
  PottsState(arma::mat mu, arma::mat lambda, arma::rowvec z): mu(mu), lambda(lambda), z(z) {}
  
  arma::mat sigma() const { return inv(lambda); }
};

// TODO: add lambda0, mu0?
class ClusterParams {
public:
  // Keep public for ease of access (TODO: set up accessors later)
  arma::mat Y;         // expression/feature matrix (n x d)
  int q;               // number of clusters
  double alpha, beta;  // Wishart hyperparameters
  double gamma;        // Smoothing hyperparameter
  
  // Initial cluster means and covariance
  arma::vec mu0;
  arma::mat lambda0;
  
  ClusterParams () = delete;
  
  ClusterParams (arma::mat features, int n_clusters) {
    Y = features;
    q = n_clusters;
    alpha = kDefaultAlpha;
    beta = kDefaultBeta;
    gamma = kDefaultGamma;
    
    setDefaultPriorParameters();
  } 
  
  ClusterParams (arma::mat features, int n_clusters,
                 double a, double b, double g) {
    Y = features;
    q = n_clusters;
    alpha = a;
    beta = b;
    gamma = g;
    
    setDefaultPriorParameters();
  }
  
  // TODO: revisit this - 
  // maybe make  mu0 and lambda0 optional arguments to cluster_mcmc()?
  void setDefaultPriorParameters () {
    // Default mu0 are means across all spots.
    // Y is (spots x features), so mean() returns a rowvec that we need to transpose
    mu0 = mean(Y, 0).t();
    lambda0 = arma::eye(arma::size(n_dims(), n_dims())) * kDefaultLambda;
    
    assert(mu0.n_elem == n_dims());
  }
  
  // No need to track n and d directly
  int n_spots() const { return Y.n_rows; }
  int n_dims() const { return Y.n_cols; }
};

// Compute next mean vector
arma::mat update_mu(ClusterParams params, PottsState prev) {
  
  int n_i;
  arma::uvec index_1k;
  arma::vec mean_i;
  arma::mat var_i;
  NumericVector Ysums;
    
  arma::mat mu_i(params.q, params.n_dims(), arma::fill::zeros);
    
  for (int k = 1; k <= params.q; k++){
    index_1k = arma::find(prev.z == k);
    n_i = index_1k.n_elem;
    Ysums = sum(params.Y.rows(index_1k), 0);

    mean_i = arma::inv(params.lambda0 + n_i * prev.lambda) * (params.lambda0 * params.mu0 + prev.lambda * as<arma::colvec>(Ysums));
    var_i = arma::inv(params.lambda0 + n_i * prev.lambda);

    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
  }
    
  return mu_i;
}

// Compute next covariance matrix
arma::mat update_lambda(ClusterParams params, PottsState curr, PottsState prev) {
  
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

double computeLogLikelihood(arma::rowvec x, int label, PottsState state) {
  double logP_YGivenZ;
  
  // Mean vector of given cluster
  arma::vec mean = arma::vectorise(state.mu.row(label - 1));
  
  logP_YGivenZ = dmvnorm(x, mean, state.sigma(), true)[0];
  return logP_YGivenZ;
}

// Compute next set of cluster assignments and corresponding log-likelihoods
// TODO: extract energy/likelihood computation
arma::mat update_z(ClusterParams params, PottsState curr, PottsState prev, List df_j) {
  
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
      h_z_prev = computeEnergy(z, neighbors, z_j_prev, params) + computeLogLikelihood(params.Y.row(j), z_j_prev, curr);
      h_z_new  = computeEnergy(z, neighbors, z_j_new, params) + computeLogLikelihood(params.Y.row(j), z_j_new, curr);
    } else {
      h_z_prev = computeLogLikelihood(params.Y.row(j), z_j_prev, curr);
      h_z_new  = computeLogLikelihood(params.Y.row(j), z_j_new, curr);
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

// Cluster input matrix Y
// [[Rcpp::export]]
List cluster_mcmc(arma::mat Y, List df_j, int nrep, int n_spots, int n_dims, double gamma, 
                  int n_clusters, arma::vec init, NumericVector mu0, arma::mat lambda0, 
                  double alpha, double beta){
  
  ClusterParams params(Y, n_clusters, alpha, beta, gamma);
  std::vector<PottsState> chain;
  
  //Initalize matrices storing iterations
  arma::mat df_sim_z(nrep, n_spots, arma::fill::zeros);
  arma::mat df_sim_mu(nrep, n_clusters * n_dims, arma::fill::zeros);
  List df_sim_lambda(nrep);
  
  // TODO: add this attribute and calculation to PottsState class
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  arma::mat mu_init = repmat(params.mu0, n_clusters, 1);
  PottsState curr(mu_init, lambda0, init.t());
  chain.push_back(curr);
  
  for (int i = 1; i < nrep; i++){
    curr = PottsState();
    
    //Update mu
    curr.mu = update_mu(params, chain.back());
    
    //Update lambda
    curr.lambda = update_lambda(params, curr, chain.back());
    
    //Update z
    arma::mat result = update_z(params, curr, chain.back(), df_j);
    curr.z = result.row(0);
    plogLik[i] = arma::sum(result.row(1));
    
    // Original output format
    // TODO: write accumulator to genarate df_sim_* matrices from chain
    df_sim_mu.row(i) = curr.mu.as_row();
    df_sim_lambda[i] = curr.lambda;
    df_sim_z.row(i) = curr.z; 
    
    chain.push_back(curr);
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

context("Inputs class") {
  arma::mat Y = arma::randu<arma::mat>(5, 7);
  int q = 3;
  
  ClusterParams params(Y, q);
  
  test_that("rows match input") {
    expect_true(params.n_spots() == Y.n_rows);
  }
  test_that("columns match input") {
    expect_true(params.n_dims() == Y.n_cols); 
  }
  test_that("clusters match input") {
    expect_true(params.q == q);
  }
  test_that("default alpha matches") {
    expect_true(params.alpha == kDefaultAlpha);
  }
}