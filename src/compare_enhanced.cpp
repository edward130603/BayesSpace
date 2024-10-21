// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T> class RedMat {
private:
  arma::Mat<T> memb;

public:
  RedMat()                          = delete;
  RedMat(const RedMat &)            = default;
  RedMat(RedMat &&)                 = delete;
  RedMat &operator=(const RedMat &) = delete;
  RedMat &operator=(RedMat &&)      = delete;

  explicit RedMat(arma::uword n_row, arma::uword n_col) {
    this->memb = arma::Mat<T>(n_row, n_col);
  };

  explicit RedMat(arma::Mat<T> &&v) { this->memb = arma::Mat<T>(std::move(v)); };

  RedMat<T> &operator+=(const RedMat<T> &v) {
    this->memb = arma::join_cols(this->memb, v.memb);

    return *this;
  };

  RedMat<T> &operator+=(RedMat<T> &&v) {
    this->memb = arma::join_cols(this->memb, v.memb);

    return *this;
  };

  arma::Mat<T> get_mat() const { return this->memb; };
};

#pragma omp declare reduction(comb_mat : RedMat<double> : omp_out += omp_in)   \
    initializer(omp_priv = omp_orig)

// [[Rcpp::export]]
arma::mat
map_subspot2ref(
    const arma::mat &subspot_coords, const arma::mat &ref_coords,
    int thread_num = 1
) {
  RedMat<double> ret(0, 3);

#ifdef _OPENMP
  double t1 = 0, t2 = 0;

  omp_set_max_active_levels(2);
  omp_set_num_threads(thread_num);

  t1 = omp_get_wtime();
#endif

#pragma omp parallel for reduction(comb_mat : ret)
  for (arma::uword i = 0; i < subspot_coords.n_rows; i++) {

    const arma::vec dists = arma::sum(
        (subspot_coords.row(i) - ref_coords.each_row()).transform([](double v) {
          return std::pow(v, 2);
        }),
        1
    );

    const arma::uvec min_indices = arma::find(dists == arma::min(dists));

    arma::mat __ret(min_indices.n_elem, 3);
    __ret.col(0) = arma::ones<arma::mat>(min_indices.n_elem, 1) * i;
    __ret.col(1) = arma::ones<arma::mat>(min_indices.n_elem, 1) % min_indices;
    __ret.col(2) =
        arma::ones<arma::mat>(min_indices.n_elem, 1) * arma::min(dists);

    ret += RedMat<double>(std::move(__ret));
  }

#ifdef _OPENMP
  t2 = omp_get_wtime();

  Rcpp::Rcout << "[DEBUG] Mapping takes " << t2 - t1 << " seconds." << std::endl;
#endif

  return ret.get_mat();
}

// [[Rcpp::export]]
arma::mat
compute_corr(
    const arma::mat &m1,
    const arma::mat &m2,
    int thread_num = 1
) {
  if (m1.n_rows != m2.n_rows || m1.n_cols != m2.n_cols) {
    Rcpp::stop("Input matrices should be of the same shape.");
  }

  RedMat<double> ret(0, 2);

#ifdef _OPENMP
  double t1, t2;

  omp_set_max_active_levels(2);
  omp_set_num_threads(thread_num);

  t1 = omp_get_wtime();
#endif

#pragma omp parallel for reduction(comb_mat : ret)
  for (arma::uword i = 0; i < m1.n_rows; i++) {
    arma::mat __ret(1, 2);
    __ret(0, 0) = i;
    __ret(0, 1) = arma::cor(m1.row(i), m2.row(i)).eval()(0, 0);

    ret += RedMat<double>(std::move(__ret));
  }

#ifdef _OPENMP
  t2 = omp_get_wtime();

  Rcpp::Rcout << "[DEBUG] Computing correlation takes " << t2 - t1 << " seconds." << std::endl;
#endif

  return ret.get_mat();
}
