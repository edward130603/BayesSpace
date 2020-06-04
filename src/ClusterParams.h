#ifndef CLUSTER_PARAMS_H
#define CLUSTER_PARAMS_H

#include <RcppArmadillo.h>
#include <testthat.h>

// Hyperparameter defaults
const double kDefaultAlpha = 1;
const double kDefaultBeta = 0.01;
const double kDefaultGamma = 2;
const double kDefaultLambda = 0.01;

// TODO: add lambda0, mu0?
class ClusterParams {
public:
  // Keep public for ease of access (TODO: set up accessors later)
  arma::mat Y;         // expression/feature matrix (n x d)
  int q;               // number of clusters
  double alpha, beta;  // Wishart hyperparameters
  double gamma;        // Smoothing hyperparameter
  // TODO: add neighbors matrix/list
  
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

#endif