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
  
  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();

  //Iterate
  colvec mu0vec = as<colvec>(mu0);
  for (int i = 1; i < nrep; i++){

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
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++){
      int z_j_prev = df_sim_z(i,j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      uvec j_vector = df_j[j];
      uvec i_vector(1); i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
      } else {
        h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
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

// [[Rcpp::export]]
List iterate_deconv(mat Y, List df_j, int nrep, int n, int n0, int d, double gamma, int q, vec init, NumericVector mu0, mat lambda0, double alpha, double beta){

  //Initalize matrices storing iterations
  mat Y0 = Y.rows(0, n0-1);
  mat df_sim_z(nrep, n, fill::zeros);
  mat df_sim_mu(nrep, q*d, fill::zeros);
  List df_sim_lambda(nrep);
  List df_sim_Y(nrep);
  // NumericVector plogLik(nrep, NA_REAL);

  //Initialize parameters
  rowvec initmu = rep(mu0, q);
  df_sim_mu.row(0) = initmu;
  df_sim_lambda[0] = lambda0;
  df_sim_z.row(0) = init.t();
  //Iterate
  mat mu_i(q,d);
  mat mu_i_long(n,d);
  for (int i = 1; i < nrep; i++){

    //Update mu
    mat lambda_prev = df_sim_lambda[i-1];
    NumericVector Ysums;
    for (int k = 1; k <= q; k++){
      uvec index_1k = find(df_sim_z.row(i-1) == k);
      int n_i = index_1k.n_elem;
      mat Yrows = Y.rows(index_1k);
      Ysums = sum(Yrows, 0);
      vec mean_i = inv(lambda0 + n_i * lambda_prev) * (lambda0 * mu0vec + lambda_prev * as<colvec>(Ysums));
      mat var_i = inv(lambda0 + n_i * lambda_prev);
      mu_i.row(k-1) = rmvnorm(1, mean_i, var_i);
    }
    df_sim_mu.row(i) = vectorise(mu_i, 1);
    
    //Update lambda
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
    
    //Update Y
    for (int j0 = 0; j0 < n0; j0++){
      
    }
    
    // for (j in 1:n0){
    //   Y_j_prev = df_sim_Y[[i-1]][four_map[[j]],]
    //   error0 = mvnfast::rmvn(n = 9, rep(0, d), sigma = diag(d)/200)
    //   error =  t(t(error0) - colMeans(error0))
    //   Y_j_new = Y_j_prev + error
    //   mu_i_four = mu_i[df_sim_z[i-1,j + 0:8 * n0],]
    //   p_prev = sum(sapply(1:9, function(x){mvnfast::dmvn(Y_j_prev[x,], mu_i_four[x,], sigma_i, log = T)}))
    //   p_new  = sum(sapply(1:9, function(x){mvnfast::dmvn(Y_j_new[x,] , mu_i_four[x,], sigma_i, log = T)}))
    //   probY_j = min(exp(p_new - p_prev) * exp(-0.1*(sum(diag(crossprod(df_sim_Y[[1]][rep(j,9),] - Y_j_new))) -
    //     sum(diag(crossprod(df_sim_Y[[1]][rep(j,9),] - Y_j_prev))))), 1)
    //   test[j] = sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j)) #(testing purpose only!!)
    //   if (sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j))){
    //     df_sim_Y[[i]][four_map[[j]],] = Y_j_new
    //   } else {
    //     df_sim_Y[[i]][four_map[[j]],] = Y_j_prev
    //   }
    // }
    
    //Update z
    df_sim_z.row(i) = df_sim_z.row(i-1);
    IntegerVector qvec = seq_len(q);
    NumericVector plogLikj(n, NA_REAL);
    for (int j = 0; j < n; j++){
      int z_j_prev = df_sim_z(i,j);
      IntegerVector qlessk = qvec[qvec != z_j_prev];
      int z_j_new = sample(qlessk, 1)[0];
      uvec j_vector = df_j[j];
      uvec i_vector(1); i_vector.fill(i);
      double h_z_prev;
      double h_z_new;
      if (j_vector.size() != 0){
        h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
      } else {
        h_z_prev = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_prev-1)), sigma_i, true)[0];
        h_z_new  = dmvnorm(Y.row(j), vectorise(mu_i.row(z_j_new -1)), sigma_i, true)[0];
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
    
    // df_sim_z.row(i) = df_sim_z.row(i-1);
    // IntegerVector qvec = seq_len(q);
    // NumericVector plogLikj(n, NA_REAL);
    // for (int j = 0; j < n; j++){
    //   int z_j_prev = df_sim_z(i,j);
    //   IntegerVector qlessk = qvec[qvec != z_j_prev];
    //   int z_j_new = sample(qlessk, 1)[0];
    //   uvec j_vector = df_j[j];
    //   // uvec i_vector(1); i_vector.fill(i);
    // 
    //   uvec z_js = conv_to<uvec>::from(df_sim_z(i_vector, j0_vector*n0+j%n0));
    //   mat mu_js_prev = mu_i.rows(z_js-1);
    //   vec mu_j_prev = vectorise(sum(mu_js_prev)/9);
    //   // Rcout << "j = "<< j <<std::endl <<"old: "  << mu_j_prev  << "z_j_prev: " << z_j_prev << std::endl;
    //   z_js[j/n0] = z_j_new;
    //   mat mu_js_new = mu_i.rows(z_js-1);
    //   vec mu_j_new = vectorise(sum(mu_js_new)/9);
    //   // Rcout << "new: " << mu_j_new << "z_j_new: " << z_j_new << std::endl;
    // 
    //   double h_z_prev;
    //   double h_z_new;
    //   if (j_vector.size() != 0){
    //     h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
    //     h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
    //   } else {
    //     h_z_prev = dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
    //     h_z_new  = dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
    //   }
    //   double prob_j = exp(h_z_new-h_z_prev);
    //   if (prob_j > 1){
    //     prob_j = 1;
    //   }
    //   IntegerVector zsample = {z_j_prev, z_j_new};
    //   NumericVector probs = {1-prob_j, prob_j};
    //   df_sim_z(i,j) = sample(zsample, 1, true, probs)[0];
    //   plogLikj[j] = h_z_prev;
    // }
  }
  List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda);
  return(out);
}



// // [[Rcpp::export]]
// List iterate2(mat Y, List df_j, int nrep, int n, int n0, int d, double gamma, int q, vec init, NumericVector mu0, mat lambda0, double alpha, double beta){
//   
//   //Initalize matrices storing iterations
//   mat Y0 = Y.rows(0, n0-1);
//   mat df_sim_z(nrep, n, fill::zeros);
//   mat df_sim_mu(nrep, q*d, fill::zeros);
//   List df_sim_lambda(nrep);
//   NumericVector plogLik(nrep, NA_REAL);
//   
//   //Initialize parameters
//   rowvec initmu = rep(mu0, q);
//   df_sim_mu.row(0) = initmu;
//   df_sim_lambda[0] = lambda0;
//   df_sim_z.row(0) = init.t();
//   //Iterate
//   colvec mu0vec = as<colvec>(mu0);
//   mat mu_i(q,d);
//   for (int k = 1; k <= q; k++){
//     uvec index_1k = find(df_sim_z.row(0) == k);
//     int n_i = index_1k.n_elem;
//     mat Yrows = Y.rows(index_1k);
//     rowvec mean_i = sum(Yrows, 0)/n_i;
//     mu_i.row(k-1) = mean_i;
//   }
//   for (int i = 1; i < nrep; i++){
//     
//     //Update mu
//     df_sim_mu.row(i) = vectorise(mu_i, 1);
// 
//     //Update lambda
//     mat lambda_prev = df_sim_lambda[i-1];
//     mat mu_i_long(n0,d);
//     uvec i_vector(1); i_vector.fill(i);
//     uvec j0_vector = {0,1,2,3,4,5,6,7,8};
//     for (int j = 0; j < n0; j++){
//       uvec z_js = conv_to<uvec>::from(df_sim_z(i_vector-1, j0_vector*n0+j) -1);
//       mat mu_j = mu_i.rows(z_js);
//       mu_i_long.row(j) = sum(mu_j)/9;
//     }
//     mat sumofsq = (Y0-mu_i_long).t() * (Y0-mu_i_long);
//     vec beta_d(d);
//     beta_d.fill(beta);
//     mat Vinv = diagmat(beta_d);
//     mat lambda_i = rwish(n0 + alpha, inv(Vinv + sumofsq));
//     df_sim_lambda[i] = lambda_i;
//     mat sigma_i = inv(lambda_i);
// 
//     //Update z
//     df_sim_z.row(i) = df_sim_z.row(i-1);
//     IntegerVector qvec = seq_len(q);
//     NumericVector plogLikj(n, NA_REAL);
//     for (int j = 0; j < n; j++){
//       int z_j_prev = df_sim_z(i,j);
//       IntegerVector qlessk = qvec[qvec != z_j_prev];
//       int z_j_new = sample(qlessk, 1)[0];
//       uvec j_vector = df_j[j];
//       // uvec i_vector(1); i_vector.fill(i);
//       
//       uvec z_js = conv_to<uvec>::from(df_sim_z(i_vector, j0_vector*n0+j%n0));
//       mat mu_js_prev = mu_i.rows(z_js-1);
//       vec mu_j_prev = vectorise(sum(mu_js_prev)/9);
//       // Rcout << "j = "<< j <<std::endl <<"old: "  << mu_j_prev  << "z_j_prev: " << z_j_prev << std::endl;
//       z_js[j/n0] = z_j_new;
//       mat mu_js_new = mu_i.rows(z_js-1);
//       vec mu_j_new = vectorise(sum(mu_js_new)/9);
//       // Rcout << "new: " << mu_j_new << "z_j_new: " << z_j_new << std::endl;
// 
//       double h_z_prev;
//       double h_z_new;
//       if (j_vector.size() != 0){
//         h_z_prev = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_prev)) + dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
//         h_z_new  = gamma/j_vector.size() * 2*accu((df_sim_z(i_vector, j_vector) == z_j_new )) + dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
//       } else {
//         h_z_prev = dmvnorm(Y.row(j), mu_j_prev, sigma_i, true)[0];
//         h_z_new  = dmvnorm(Y.row(j), mu_j_new , sigma_i, true)[0];
//       }
//       double prob_j = exp(h_z_new-h_z_prev);
//       if (prob_j > 1){
//         prob_j = 1;
//       }
//       IntegerVector zsample = {z_j_prev, z_j_new};
//       NumericVector probs = {1-prob_j, prob_j};
//       df_sim_z(i,j) = sample(zsample, 1, true, probs)[0];
//       plogLikj[j] = h_z_prev;
//     }
//     plogLik[i] = sum(plogLikj);
//   }
//   List out = List::create(_["z"] = df_sim_z, _["mu"] = df_sim_mu, _["lambda"] = df_sim_lambda, _["plogLik"] = plogLik);
//   return(out);
// }

/*** R
#iterate2(Y = Y, df_j = df_j, nrep = nrep, n = n, n0 = n0, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
*/
