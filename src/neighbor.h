#ifndef BAYESSPACE_NEIGHBOR_H
#define BAYESSPACE_NEIGHBOR_H

#include <RcppArmadillo.h>
#include <string>

#include "utils.h"

class Neighbor {
private:
  arma::uvec neighbors;

public:
  Neighbor();
  explicit Neighbor(const arma::uvec &);
  explicit Neighbor(arma::uvec &&);
  explicit Neighbor(const arma::uvec &&);
  explicit Neighbor(const Rcpp::NumericVector &, const bool);
  explicit Neighbor(const std::string &, const std::string &, const bool);
  Neighbor(const Neighbor &);
  Neighbor(Neighbor &&);
  Neighbor(const Neighbor &&);
  ~Neighbor();

  void set_neighbors(const arma::uvec &);
  void set_neighbors(arma::uvec &&);
  void set_neighbors(const arma::uvec &&);
  void set_neighbors(const Rcpp::NumericVector &, const bool);
  void set_neighbors(const std::string &, const std::string &, const bool);

  size_t get_size() const;
  const arma::uvec &get_neighbors() const;
};

Neighbor::Neighbor() {
  this->neighbors.reset();
}

Neighbor::Neighbor(const arma::uvec &v) {
  if (v.size() > 0) {
    this->neighbors = v;
  } else {
    this->neighbors.reset();
  }
}

Neighbor::Neighbor(arma::uvec &&v) : neighbors(std::move(v)) {}

Neighbor::Neighbor(const arma::uvec &&v) : neighbors(std::move(v)) {}

Neighbor::Neighbor(const Rcpp::NumericVector &v, const bool is_zero_based) {
  if (v.size() > 0) {
    this->neighbors = arma::uvec(v.size());

    for (auto i = 0; i != v.size(); i++) {
      this->neighbors(i) = v[i];
    }

    if (!is_zero_based) {
      this->neighbors = -this->neighbors - 1;
    }
  } else {
    this->neighbors.reset();
  }
}

Neighbor::Neighbor(
    const std::string &v, const std::string &separator, const bool is_zero_based
) {
  const std::vector<arma::uword> vals = str_split<arma::uword>(v, separator);

  if (vals.size() > 0) {
    this->neighbors = arma::uvec(std::move(vals));

    if (!is_zero_based)
      this->neighbors = this->neighbors - 1;
  } else {
    this->neighbors.reset();
  }
}

Neighbor::Neighbor(const Neighbor &v) { this->neighbors = v.neighbors; }

Neighbor::Neighbor(Neighbor &&v) : neighbors(std::move(v.neighbors)) {}

Neighbor::Neighbor(const Neighbor &&v) : neighbors(std::move(v.neighbors)) {}

Neighbor::~Neighbor() { this->neighbors.reset(); }

void
Neighbor::set_neighbors(const arma::uvec &v) {
  if (v.size() > 0) {
    this->neighbors = v;
  } else {
    this->neighbors.reset();
  }
}

void
Neighbor::set_neighbors(arma::uvec &&v) {
  this->neighbors = arma::uvec(std::move(v));
}

void
Neighbor::set_neighbors(const arma::uvec &&v) {
  this->neighbors = arma::uvec(std::move(v));
}

void
Neighbor::set_neighbors(const Rcpp::NumericVector &v, const bool is_zero_based) {
  if (v.size() > 0) {
    this->neighbors = arma::uvec(v.size());

    for (auto i = 0; i != v.size(); i++) {
      this->neighbors(i) = v[i];
    }

    if (!is_zero_based)
      this->neighbors = this->neighbors - 1;
  } else {
    this->neighbors.reset();
  }
}

void
Neighbor::set_neighbors(
    const std::string &v, const std::string &separator, const bool is_zero_based
) {
  const std::vector<arma::uword> vals = str_split<arma::uword>(v, separator);

  if (vals.size() > 0) {
    this->neighbors = arma::uvec(std::move(vals));

    if (!is_zero_based)
      this->neighbors = this->neighbors - 1;
  } else {
    this->neighbors.reset();
  }
}

size_t
Neighbor::get_size() const {
  return this->neighbors.n_elem;
}

const arma::uvec &
Neighbor::get_neighbors() const {
  return this->neighbors;
}

#endif   // BAYESSPACE_NEIGHBOR_H
