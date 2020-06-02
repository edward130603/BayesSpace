/* Copyright Fred Hutch
 * License <TODO>
 * 
 * Parameters associated with a given Potts model and functions to compute
 * their likelihoods
 */

#ifndef POTTS_MODEL_H
#define POTTS_MODEL_H

#include <RcppArmadillo.h>
#include "ClusterParams.h"
using namespace Rcpp;

// TODO: store reference to params in model so it doesn't need to be passed to each update
class PottsModel {
public:
  arma::mat mu;
  arma::mat lambda;
  arma::rowvec z;
  
  PottsModel(): mu(arma::zeros(0, 0)), lambda(arma::zeros(0, 0)), z(arma::zeros<arma::rowvec>(0)) {}
  PottsModel(arma::mat mu, arma::mat lambda, arma::rowvec z): mu(mu), lambda(lambda), z(z) {}
  
  arma::mat sigma() const { return inv(lambda); }
};

// TODO: convert updates to member functions
arma::mat updateMu(ClusterParams params, PottsModel prev);
arma::mat updateLambda(ClusterParams params, PottsModel curr, PottsModel prev);
arma::mat updateZ(ClusterParams params, PottsModel curr, PottsModel prev, List df_j); // TODO: get rid of List

// TODO: add likelihood/energy/poster computations

#endif