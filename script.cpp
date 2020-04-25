#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
List iterate(mat Y, List df_j, int nrep, int n, int d, double gamma, int q, vec init, NumericVector mu0, mat lambda0, double alpha, double beta){

  //Initalize matrices storing iterations
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep);
  NumericVector plogLik(nrep, NA_REAL);
  // vec gammavec = linspace(0, gamma, nrep);
  // Rcout << gammavec;

  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();

  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++){
    // gamma = gammavec[i];
    
    //Update mu
    mat mu_i(q,d);
    mat lambda_prev = df_sim_lambda[i-1];
    for (int k = 1; k <= q; k++){
      uvec index_1k = find(df_sim_z.row(i-1) == k);
      int n_i = index_1k.n_elem;
      NumericVector Ysums;
      for (int di = 0; di < d; di++){
        mat Yrows = Y.rows(index_1k);
        Ysums.push_back(sum(Yrows.col(di)));
      }
      vec mean_i = inv(lambda0 + n_i * lambda_prev) * (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums));
      mat var_i = inv(lambda0 + n_i * lambda_prev);
      mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);

    //Update lambda
    mat mu_i_long(n,d);
    for (int j = 0; j < n; j++){
      mu_i_long.row(j) = mu_i.row(df_sim_z(i-1, j)-1);
    }
    mat sumofsq = (Y-mu_i_long).t() * (Y-mu_i_long);
    vec beta_d(d); 
    beta_d.fill(beta);
    mat Vinv = diagmat(beta_d);
    mat lambda_i = rwish(n + alpha, inv(Vinv + sumofsq));
    df_sim_lambda[i] = lambda_i;
    mat sigma_i = inv(lambda_i);

    //Update z
    df_sim_z.row(i) = df_sim_z.row(i-1);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++){
      int z_j_prev = df_sim_z(i,j);
      IntegerVector qvec = seq_len(q);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      uvec j_vector = df_j[j];
      uvec i_vector(1); i_vector.fill(i);
      double h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
      // Rcout << "prior prev: " << gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) << std:: endl;
      // Rcout << "lik prev: " << dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0] << std:: endl;
      double h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
      double prob_j = exp(h_z_new-h_z_prev);
      // Rcout << "prior new: " << gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new)) << std:: endl;
      // Rcout << "lik new: " << dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new-1)), sigma_i, true)[0] << std:: endl;
      if (prob_j > 1){
        prob_j = 1;
      }
      IntegerVector zsample = {z_j_prev, z_j_new};
      NumericVector probs = {1-prob_j, prob_j};
      // NumericVector probs = {0, 1};
      df_sim_z(i,j) = sample(zsample, 1, true, probs)[0];
      // Rcout << "probs: " << probs << " for " << zsample << std::endl;
      // Rcout << sample(zsample, 1, true, probs)[0] << std::endl;
      plogLikj[j] = h_z_prev;
    }
    plogLik[i] = sum(plogLikj);
  }
  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["plogLik"] = plogLik);
  return(out);
}

/*** R

# df_sim_z[i,] = df_sim_z[i-1, ]
# plogLikj = rep(NA, n)
# denom = rep(NA, n)
# for (j in 1:n){
#   z_j_prev = df_sim_z[i,j]
#   qlessk = setdiff(1:q, z_j_prev)
#   z_j_new = sample(qlessk, 1)
#   j_vector = df_j[[j]]
#   h_z_prev = gamma/length(j_vector)* 2*sum(((z_j_prev == df_sim_z[i, j_vector])-0.5)) + mvnfast::dmvn(Y[j,], mu = mu_i[z_j_prev,], sigma = sigma_i, log = T)
#   h_z_new = gamma/length(j_vector) * 2*sum(((z_j_new  == df_sim_z[i, j_vector])-0.5)) + mvnfast::dmvn(Y[j,], mu = mu_i[z_j_new, ], sigma = sigma_i, log = T)

#   prob_j = min(exp(h_z_new - h_z_prev),1)
#   df_sim_z[i, j] = sample(x = c(z_j_prev, z_j_new), size = 1, prob = c(1-prob_j, prob_j))
#   plogLikj[j] = h_z_prev
# }
# set.seed(100)
# iterate(as.matrix(Y), df_j, nrep = 2, n = 3639, d = 10, gamma, q, init, mu0, lambda0, alpha, beta)
*/
