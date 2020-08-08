/* 
 * Fast dmnvorm implementation
 *
 * Borrowed with appreciation from Nino Hardt, Dicko Ahmadou, Ben Christofferson
 * 2020-08-04
 * https://gallery.rcpp.org/articles/dmvnorm_arma/
 *
 */

#ifndef DMVNORM_H
#define DMVNORM_H

#include <RcppArmadillo.h>

arma::vec Mahalanobis(arma::mat const &x, 
                      arma::vec const &center, 
                      arma::mat const &cov);

arma::vec dmvnorm_arma(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false);

#endif
