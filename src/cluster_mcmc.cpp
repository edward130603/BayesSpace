// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <testthat.h>
#include "ClusterParams.h"
#include "PottsModel.h"
using namespace Rcpp;


// [[Rcpp::export]]
List cluster_mcmc(arma::mat Y, List df_j, int nrep, int n_spots, int n_dims, double gamma, 
                  int n_clusters, arma::vec init, NumericVector mu0, arma::mat lambda0, 
                  double alpha, double beta){
  
  ClusterParams params(Y, n_clusters, alpha, beta, gamma);
  std::vector<PottsModel> chain;
  
  //Initalize matrices storing iterations
  arma::mat df_sim_z(nrep, n_spots, arma::fill::zeros);
  arma::mat df_sim_mu(nrep, n_clusters * n_dims, arma::fill::zeros);
  List df_sim_lambda(nrep);
  
  // TODO: add this attribute and calculation to PottsModel class
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  arma::mat mu_init = repmat(params.mu0, n_clusters, 1);
  PottsModel curr(mu_init, lambda0, init.t());
  chain.push_back(curr);
  
  for (int i = 1; i < nrep; i++){
    curr = PottsModel();
    
    //Update mu
    curr.mu = curr.updateMu(params, chain.back());
    
    //Update lambda
    curr.lambda = curr.updateLambda(params, chain.back());
    
    //Update z
    arma::mat result = curr.updateZ(params, chain.back(), df_j);
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