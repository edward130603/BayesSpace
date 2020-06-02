#include <testthat.h>
#include "ClusterParams.h"

context("Test Params construction") {
  arma::mat Y = arma::randu<arma::mat>(5, 7);
  int q = 3;
  
  ClusterParams params(Y, q);
  
  test_that("rows match input") {
    expect_true(params.n_spots() == Y.n_rows);
  }
  test_that("columns match input") {
    expect_true(params.n_dims() == Y.n_cols); 
  }
  test_that("clusters match input") {
    expect_true(params.q == q);
  }
  test_that("default alpha matches") {
    expect_true(params.alpha == kDefaultAlpha);
  }
}