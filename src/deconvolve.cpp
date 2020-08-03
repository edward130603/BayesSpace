// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List iterate_deconv(arma::mat Y, List df_j, bool tdist, int nrep, int n, int n0, 
  int d, double gamma, int q, arma::vec init, int subspots, bool verbose, 
  double jitter_scale, double c, NumericVector mu0, arma::mat lambda0, 
  double alpha, double beta) {
  
  //Initalize matrices storing iterations
  mat Y0 = Y.rows(0, n0-1);
  mat df_sim_z(nrep/100 + 1, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep/100 + 1);
  List df_sim_Y(nrep/100 + 1);
  mat df_sim_w(nrep/100 + 1, n);
  NumericVector Ychange(nrep, NA_REAL);
  
  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  mat lambda_i = lambda0;
  df_sim_z.row(0) = init.t();
  vec z = init;
  df_sim_Y[0] = Y;
  vec w = ones<vec>(n);
  df_sim_w.row(0) = w.t();
  
  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  mat mu_i(q,d);
  mat mu_i_long(n,d);
  uvec j0_vector;
  if (subspots == 6){
    j0_vector = {0,1,2,3,4,5};
  } else {
    j0_vector = {0,1,2,3,4,5,6,7,8};
  }
  mat Y_j_prev(subspots,d);
  mat Y_j_new(subspots,d);
  mat error(subspots,d);
  vec zero_vec = zeros<vec>(d);
  vec one_vec = ones<vec>(d);
  mat error_var = diagmat(one_vec)/d*jitter_scale; 
  Progress p(nrep - 1, verbose);
  for (int i = 1; i < nrep; i++){
    p.increment();
    if (i % 10 == 0){
      if (Progress::check_abort())
        return(List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange));
    }
    
    //Update mu
    NumericVector Ysums;
    for (int k = 1; k <= q; k++){
      uvec index_1k = find(z == k);
      int n_i = sum(w(index_1k)); 
      mat Yrows = Y.rows(index_1k);
      Yrows.each_col() %= w(index_1k);
      Ysums = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_i) * (lambda0 * mu0vec + lambda_i * as<colvec>(Ysums));
      mat var_i = inv(lambda0 + n_i * lambda_i);
      mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);
    
    //Update lambda
    for (int j = 0; j < n; j++){
      mu_i_long.row(j) = mu_i.row(z[j]-1);
    }
    mat sumofsq = (Y-mu_i_long).t() * diagmat(w) *  (Y-mu_i_long);
    vec beta_d(d);
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    lambda_i = rwish(n + alpha, inv(Vinv + sumofsq));
    mat sigma_i = inv(lambda_i);
    
    //Update Y
    int updateCounter = 0;
    for (int j0 = 0; j0 < n0; j0++){
      Y_j_prev = Y.rows(j0_vector*n0+j0);
      error = rmvnorm(subspots, zero_vec, error_var); 
      rowvec error_mean = sum(error, 0)/subspots;
      for (int r = 0; r < subspots; r++){
        error.row(r) = error.row(r) - error_mean;
      }
      Y_j_new = Y_j_prev + error;
      mat mu_i_j = mu_i_long.rows(j0_vector*n0+j0);
      vec p_prev = {0.0};
      vec p_new = {0.0};
      for (int r = 0; r< subspots; r++){
        p_prev += dmvnorm(Y_j_prev.row(r), vectorise(mu_i_j.row(r)), sigma_i/w[j0 + n0*r], true) - c*(accu(pow(Y_j_prev.row(r)-Y0.row(j0),2)));
        p_new  += dmvnorm(Y_j_new.row(r),  vectorise(mu_i_j.row(r)), sigma_i/w[j0 + n0*r], true) - c*(accu(pow(Y_j_new.row(r) -Y0.row(j0),2)));
      }
      double probY_j = as_scalar(exp(p_new - p_prev));
      if (probY_j > 1){
        probY_j = 1;
      }
      IntegerVector Ysample = {0, 1};
      NumericVector probsY = {1-probY_j, probY_j};
      int yesUpdate = sample(Ysample, 1, true, probsY)[0];
      if (yesUpdate == 1){
        Y.rows(j0_vector*n0 + j0) = Y_j_new;
        updateCounter++;
      }
    }
    Ychange[i] = updateCounter * 1.0 / n0;
    
    //Update w and z
    double w_alpha = (d+4)/2; //shape parameter
    double w_beta;
    IntegerVector qvec = seq_len(q);
    for (int j = 0; j < n; j++){
      if (tdist){
        w_beta = as_scalar(2/((Y.row(j)-mu_i_long.row(j))* lambda_i * (Y.row(j)-mu_i_long.row(j)).t() + 4)); //scale parameter
        w[j] = R::rgamma(w_alpha, w_beta); //sample from posterior for w
      }
      int z_j_prev = z[j];
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      uvec j_vector = df_j[j];
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((z(j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i/w[j], true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((z(j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i/w[j], true)[0];
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
      z[j] = sample(zsample, 1, true, probs)[0];
    }
    if ((i+1) % 100 == 0){
      df_sim_lambda[(i+1)/100] = lambda_i;
      df_sim_Y[(i+1)/100] = Y;
      df_sim_w.row((i+1)/100) = w.t();
      df_sim_z.row((i+1)/100) = z.t();
    }
  }
  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["weights"] = df_sim_w, _["Y"] = df_sim_Y, _["Ychange"] = Ychange);
  return(out);
}
