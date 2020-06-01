// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <testthat.h>
using namespace Rcpp;

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

// TODO: add alpha, beta, gamma
// TODO: add lambda0, mu0
class ClusterInputs {
public:
  // Constructors
  ClusterInputs (): Y(arma::zeros(0, 0)), q(0) {}
  ClusterInputs (arma::mat Y, int q): Y(Y), q(q) {}
  
  // Getters
  int n_spots() const { return Y.n_rows; }
  int n_dims() const { return Y.n_cols; }
  int n_clusters() const { return q; }
  arma::mat features() const { return Y; }
  
  // Setters
  void set_Y(arma::mat features) { Y = features; }
  void set_q(int n_clusters) { q = n_clusters; }
  
private:
  arma::mat Y;  // n x d Expression/feature matrix
  int q;        // Number of clusters
};

// arma::mat update_mu(ModelState curr, ModelState prev) {
arma::mat update_mu(ClusterInputs inputs, PottsState prev,
                    arma::mat lambda0, arma::vec mu0) {
  
  int n_i;
  arma::uvec index_1k;
  arma::vec mean_i;
  arma::mat var_i;
  NumericVector Ysums;
    
  arma::mat mu_i(inputs.n_clusters(), inputs.n_dims(), arma::fill::zeros);
    
  for (int k = 1; k <= inputs.n_clusters(); k++){
    index_1k = arma::find(prev.z == k);
    n_i = index_1k.n_elem;
    Ysums = sum(inputs.features().rows(index_1k), 0);

    mean_i = inv(lambda0 + n_i * prev.lambda) * (lambda0 * mu0 + prev.lambda * as<arma::colvec>(Ysums));
    var_i = inv(lambda0 + n_i * prev.lambda);

    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
  }
    
  return mu_i;
}

arma::mat update_lambda(ClusterInputs inputs, PottsState curr, PottsState prev, 
                        double alpha, double beta) {
  
  arma::mat mu_i_long(inputs.n_spots(), inputs.n_dims(), arma::fill::zeros);
  
  for (int j = 0; j < inputs.n_spots(); j++){
    mu_i_long.row(j) = curr.mu.row(prev.z[j] - 1);
  }
  
  arma::mat sumofsq = (inputs.features() - mu_i_long).t() * (inputs.features() - mu_i_long);
  
  arma::vec beta_d(inputs.n_dims()); 
  beta_d.fill(beta);
  
  arma::mat lambda_i = rwish(inputs.n_spots() + alpha, inv(arma::diagmat(beta_d) + sumofsq));

  return lambda_i;
}

arma::mat update_z(ClusterInputs inputs, PottsState curr, PottsState prev, 
                   List df_j, double gamma) {
  
  arma::rowvec z = prev.z;
  arma::rowvec plogLikj(inputs.n_spots(), arma::fill::zeros);
  IntegerVector qvec = seq_len(inputs.n_clusters());
  double h_z_prev;
  double h_z_new;
  
  for (int j = 0; j < inputs.n_spots(); j++){
    int z_j_prev = z(j);
    IntegerVector qlessk = qvec[qvec != z_j_prev];
    int z_j_new = sample(qlessk, 1)[0];
    
    arma::uvec j_vector = df_j[j];
    
    // has neighbors
    if (j_vector.size() != 0){
      h_z_prev = gamma/j_vector.size() * 2*arma::accu((z(j_vector) == z_j_prev)) + dmvnorm(inputs.features().row(j), vectorise(curr.mu.row(z_j_prev-1)), curr.sigma(), true)[0];
      h_z_new  = gamma/j_vector.size() * 2*arma::accu((z(j_vector) == z_j_new )) + dmvnorm(inputs.features().row(j), vectorise(curr.mu.row(z_j_new -1)), curr.sigma(), true)[0];
    } else {
      h_z_prev = dmvnorm(inputs.features().row(j), vectorise(curr.mu.row(z_j_prev-1)), curr.sigma(), true)[0];
      h_z_new  = dmvnorm(inputs.features().row(j), vectorise(curr.mu.row(z_j_new -1)), curr.sigma(), true)[0];
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
  
  ClusterInputs inputs(Y, n_clusters);
  std::vector<PottsState> chain;
  
  //Initalize matrices storing iterations
  arma::mat df_sim_z(nrep, n_spots, arma::fill::zeros);
  arma::mat df_sim_mu(nrep, n_clusters * n_dims, arma::fill::zeros);
  List df_sim_lambda(nrep);
  
  // TODO: add this attribute and calculation to PottsState class
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  arma::vec mu0vec = as<arma::vec>(mu0);
  arma::mat mu_init = repmat(mu0vec, n_clusters, 1);
  PottsState curr(mu_init, lambda0, init.t());
  chain.push_back(curr);
  
  
  for (int i = 1; i < nrep; i++){
    curr = PottsState();
    
    //Update mu
    curr.mu = update_mu(inputs, chain.back(), lambda0, mu0vec);
    
    //Update lambda
    curr.lambda = update_lambda(inputs, curr, chain.back(), alpha, beta);
    
    //Update z
    arma::mat result = update_z(inputs, curr, chain.back(), df_j, gamma);
    curr.z = result.row(0);
    plogLik[i] = arma::sum(result.row(1));
    
    // Original output format
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
  
  ClusterInputs inputs(Y, q);
  
  test_that("rows match input") {
    expect_true(inputs.n_spots() == Y.n_rows);
  }
  test_that("columns match input") {
    expect_true(inputs.n_dims() == Y.n_cols); 
  }
  test_that("clusters match input") {
    expect_true(inputs.n_clusters() == q);
  }
}