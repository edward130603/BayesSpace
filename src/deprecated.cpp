// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List iterate2(arma::mat Y, List df_j, int nrep, int n, int n0, int d, double gamma, int q, arma::vec init, NumericVector mu0, arma::mat lambda0, double alpha, double beta){
  
  //Initalize matrices storing iterations
  mat Y0 = Y.rows(0, n0-1);
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep);
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  mat mu_i(q,d);
  for (int k = 1; k <= q; k++){
    uvec index_1k = find(df_sim_z.row(0) == k);
    int n_i = index_1k.n_elem;
    mat Yrows = Y.rows(index_1k);
    rowvec mean_i = sum(Yrows, 0)/n_i;
    mu_i.row(k-1) = mean_i;
  }
  for (int i = 1; i < nrep; i++){
    
    //Update mu
    df_sim_mu.row(i) = vectorise(mu_i, 1);
    
    //Update lambda
    mat lambda_prev = df_sim_lambda[i-1];
    mat mu_i_long(n0,d);
    uvec i_vector(1); i_vector.fill(i);
    uvec j0_vector = {0,1,2,3,4,5,6,7,8};
    for (int j = 0; j < n0; j++){
      uvec z_js = conv_to<uvec>::from(df_sim_z(i_vector-1, j0_vector*n0+j) -1);
      mat mu_j = mu_i.rows(z_js);
      mu_i_long.row(j) = sum(mu_j)/9;
    }
    mat sumofsq = (Y0-mu_i_long).t() * (Y0-mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    mat lambda_i = rwish(n0 + alpha, inv(Vinv + sumofsq));
    df_sim_lambda[i] = lambda_i;
    mat sigma_i = inv(lambda_i);
    
    //Update z
    df_sim_z.row(i) = df_sim_z.row(i-1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++){
      int z_j_prev = df_sim_z(i,j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      uvec j_vector = df_j[j];
      // uvec i_vector(1); i_vector.fill(i);
      
      uvec z_js = conv_to<uvec>::from(df_sim_z(i_vector, j0_vector*n0+j%n0));
      mat mu_js_prev = mu_i.rows(z_js-1);
      vec mu_j_prev = vectorise(sum(mu_js_prev)/9);
      // Rcout << "j = "<< j <<std::endl <<"old: "  << mu_j_prev  << "z_j_prev: " << z_j_prev << std::endl;
      z_js[j/n0] = z_j_new;
      mat mu_js_new = mu_i.rows(z_js-1);
      vec mu_j_new = vectorise(sum(mu_js_new)/9);
      // Rcout << "new: " << mu_j_new << "z_j_new: " << z_j_new << std::endl;
      
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
      } else {
        h_z_prev = dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
        h_z_new  = dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
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
