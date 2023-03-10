#ifndef BAYESSPACE_DOUBLE_STATES_VECTOR_H
#define BAYESSPACE_DOUBLE_STATES_VECTOR_H

#include <RcppArmadillo.h>

template <class T> class DoubleStatesVector {
private:
  arma::uvec __row_idx;
  arma::uvec __col_idx;

  // Two columns for the old (column 0) and the new (column 1) value
  arma::Mat<T> values;

public:
  explicit DoubleStatesVector(int);
  DoubleStatesVector(const DoubleStatesVector &);
  DoubleStatesVector(DoubleStatesVector &&);
  ~DoubleStatesVector();

  T &operator()(const arma::uword, const arma::uword);

  arma::subview_row<T> row(const arma::uword);
  const arma::subview_row<T> row(const arma::uword) const;
  arma::subview_col<T> col(const arma::uword);
  const arma::subview_col<T> col(const arma::uword) const;
  arma::subview<T> rows(const arma::uvec &);
  const arma::subview<T> rows(const arma::uvec &) const;
  arma::subview<T> cols(const arma::uvec &);
  const arma::subview<T> cols(const arma::uvec &) const;

  void set_col_idx(const arma::uword, const arma::uword);
  const arma::Col<T> get_current_values() const;
};

template <class T>
DoubleStatesVector<T>::DoubleStatesVector(int n_row)
    : __col_idx(n_row, arma::fill::zeros), values(n_row, 2) {
  std::vector<u_int32_t> tmp(n_row);
  std::iota(tmp.begin(), tmp.end(), 0);
  __row_idx = arma::uvec(tmp);
}

template <class T>
DoubleStatesVector<T>::DoubleStatesVector(const DoubleStatesVector &v)
    : __row_idx(v.__row_idx), __col_idx(v.__col_idx), values(v.values) {}

template <class T>
DoubleStatesVector<T>::DoubleStatesVector(DoubleStatesVector &&v)
    : __row_idx(std::move(v.__row_idx)), __col_idx(std::move(v.__col_idx)),
      values(std::move(v.values)) {}

template <class T> DoubleStatesVector<T>::~DoubleStatesVector() {
  __row_idx.reset();
  __col_idx.reset();
  values.reset();
}

template <class T>
T &
DoubleStatesVector<T>::operator()(
    const arma::uword row_idx, const arma::uword col_idx
) {
  if (col_idx > 1 || row_idx > __row_idx.n_elem) {
    throw "Invalid matrix access.";
  }

  return values(row_idx, col_idx);
}

template <class T>
arma::subview_row<T>
DoubleStatesVector<T>::row(const arma::uword idx) {
  return values.row(idx);
}

template <class T>
const arma::subview_row<T>
DoubleStatesVector<T>::row(const arma::uword idx) const {
  return values.row(idx);
}

template <class T>
arma::subview_col<T>
DoubleStatesVector<T>::col(const arma::uword idx) {
  return values.col(idx);
}

template <class T>
const arma::subview_col<T>
DoubleStatesVector<T>::col(const arma::uword idx) const {
  return values.row(idx);
}

template <class T>
arma::subview<T>
DoubleStatesVector<T>::rows(const arma::uvec &idx) {
  return values.rows(idx);
}

template <class T>
const arma::subview<T>
DoubleStatesVector<T>::rows(const arma::uvec &idx) const {
  return values.rows(idx);
}

template <class T>
arma::subview<T>
DoubleStatesVector<T>::cols(const arma::uvec &idx) {
  return values.cols(idx);
}

template <class T>
const arma::subview<T>
DoubleStatesVector<T>::cols(const arma::uvec &idx) const {
  return values.cols(idx);
}

template <class T>
void
DoubleStatesVector<T>::set_col_idx(
    const arma::uword idx, const arma::uword val
) {
  if (val > 1 || idx > __row_idx.n_elem) {
    throw "Invalid vector access.";
  }

  __col_idx(idx) = val;
}

template <class T>
const arma::Col<T>
DoubleStatesVector<T>::get_current_values() const {
  return values(__col_idx * __row_idx.n_elem + __row_idx);
}

#endif   // BAYESSPACE_DOUBLE_STATES_VECTOR_H
