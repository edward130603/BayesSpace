/* Fast dmnvorm implementation
 *
 * Borrowed with appreciation from Nino Hardt, Dicko Ahmadou, Ben Christofferson
 * 2020-08-04
 * https://gallery.rcpp.org/articles/dmvnorm_arma/
 *
 */

#ifndef DMVNRM_H
#define DMVNRM_H

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
  
arma::vec dmvnrm_arma(arma::mat const &x,  
                      arma::rowvec const &mean,  
                      arma::mat const &sigma, 
                      bool const logd = false);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false);

#endif
