// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// Borrowed with appreciation from Nino Hardt, Dicko Ahmadou, Ben Christofferson
// https://gallery.rcpp.org/articles/dmvnorm_arma/
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd=false) { 
    using arma::uword;
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}

arma::mat update_mu(const arma::mat& Y, int q, const arma::rowvec& z_prev, 
                    const arma::vec& w, const arma::mat& lambda0, 
                    const arma::mat& lambda_prev, const arma::colvec& mu0vec) {

  mat mu_i(q, Y.n_cols);
  NumericVector Ysums;

  for (int k = 1; k <= q; k++){
    uvec index_1k = find(z_prev == k);
    int n_i = sum(w(index_1k)); 
    mat Yrows = Y.rows(index_1k);
    Yrows.each_col() %= w(index_1k);
    Ysums = sum(Yrows, 0);
    vec mean_i = inv(lambda0 + n_i * lambda_prev) * (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums));  //posterior mean
    mat var_i = inv(lambda0 + n_i * lambda_prev);  //posterior variance
    mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);  //sample from posterior for mu
  }

  return mu_i;
}

arma::mat update_lambda(const arma::mat& resid, const arma::rowvec& z_prev,
                        const arma::vec& w, double alpha, double beta) {

  int d = resid.n_cols;
  int n = resid.n_rows;

  mat sumofsq = resid.t() * diagmat(w) * resid;
  vec beta_d(d); 
  beta_d.fill(beta);
  mat Vinv = diagmat(beta_d);
  mat lambda_i = rwish(n + alpha, inv(Vinv + sumofsq));

  return lambda_i;

}

// Per spot update
// Can't extract from loop b/c will change output relative to original
// TODO: replace magic 4 with defined constant
// TODO rm int j and pass row directly
double update_w_j(const arma::mat& resid, const arma::mat& lambda_i, int j) {
  int d = resid.n_cols;

  // shape parameter
  double w_alpha = (d + 4) / 2;  

  // scale parameter
  double w_beta = as_scalar(2/(resid.row(j) * lambda_i * resid.row(j).t() + 4));

  // sample from posterior for w
  double w_j = R::rgamma(w_alpha, w_beta);  

  return w_j;
}

// TODO pass df_sim_z.row(i) instead of i_vector
double compute_h_z(const arma::mat& df_sim_z, const arma::uvec& i_vector,
                   const arma::uvec& j_vector, int z_j,
                   const arma::rowvec& Y_j, const arma::rowvec& mu_k, 
                   double w_j, const arma::mat& sigma_i, double gamma) {
  
  // Spatial smoothing prior
  double pi_z = 0;
  if (j_vector.size() != 0) {
    pi_z = gamma / j_vector.size() * 2 * accu((df_sim_z(i_vector, j_vector) == z_j));
  }
  
  double h_z = pi_z + dmvnrm_arma_fast(Y_j, mu_k, sigma_i/w_j, true)[0];
  
  return h_z;
}

// [[Rcpp::export]]
List iterate_t_refactor(const arma::mat& Y, List df_j, int nrep, int n, int d,
                        double gamma, int q, arma::vec init, NumericVector mu0,
                        arma::mat lambda0, double alpha, double beta, 
                        std::string model){
  
  //Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep);
  mat df_sim_w(nrep, n);
  NumericVector plogLik(nrep, NA_REAL);
  
  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  vec w = ones<vec>(n);
  df_sim_w.row(0) = w.t();
  
  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++){
    // Check for interrupt every ~1-2 seconds (timing based on 5k spots)
    if (i % 10 == 0)
      Rcpp::checkUserInterrupt();
    
    //Update mu
    mat mu_i = update_mu(Y, q, df_sim_z.row(i-1), w, lambda0, df_sim_lambda[i-1], mu0vec);
    df_sim_mu.row(i) = vectorise(mu_i, 1);
    
    // Compute residuals
    mat mu_i_long(n,d);
    for (int j = 0; j < n; j++){
      mu_i_long.row(j) = mu_i.row(df_sim_z(i-1, j) - 1);
    }
    mat resid = Y - mu_i_long;

    //Update lambda
    mat lambda_i = update_lambda(resid, df_sim_z.row(i-1), w, alpha, beta);
    df_sim_lambda[i] = lambda_i;
    mat sigma_i = inv(lambda_i);
    
    // Update z and w
    df_sim_z.row(i) = df_sim_z.row(i-1);
    const vec qvec = linspace(1, q - 1, q - 1);
    NumericVector plogLikj(n, NA_REAL);

    for (int j = 0; j < n; j++){
      // Update spot weight (t-distributed error model)
      // Under normal error model, weights remain one
      if (model == "t") {
          w[j] = update_w_j(resid, lambda_i, j);
      }
      
      // Sample new cluster assignment (one-indexed)
      int z_j_prev = df_sim_z(i,j);
      int z_j_new = Rcpp::RcppArmadillo::sample(qvec, 1, false)[0];
      if (z_j_new >= z_j_prev) {
        z_j_new++;
      }
      
      uvec j_vector = df_j[j];
      uvec i_vector(1); i_vector.fill(i);
      double h_z_prev = compute_h_z(df_sim_z, i_vector, j_vector, z_j_prev, Y.row(j), mu_i.row(z_j_prev - 1), w[j], sigma_i, gamma);
      double h_z_new = compute_h_z(df_sim_z, i_vector, j_vector, z_j_new, Y.row(j), mu_i.row(z_j_new - 1), w[j], sigma_i, gamma);
      double prob_j = std::min(std::exp(h_z_new - h_z_prev), 1.0);

      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs = {1-prob_j, prob_j};
      df_sim_z(i,j) = Rcpp::RcppArmadillo::sample(zsample, 1, true, probs)[0];
      plogLikj[j] = h_z_prev;
    }
    
    df_sim_w.row(i) = w.t();
    plogLik[i] = sum(plogLikj);
  }

  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["weights"] = df_sim_w, _["plogLik"] = plogLik);

  return(out);
}



