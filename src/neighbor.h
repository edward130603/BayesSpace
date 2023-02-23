#ifndef BAYESSPACE_NEIGHBOR_H
#define BAYESSPACE_NEIGHBOR_H

#include <RcppArmadillo.h>

class Neighbor {
private:
  arma::uvec neighbors;

public:
  explicit Neighbor(const arma::uvec &);
  explicit Neighbor(arma::uvec &&);
  explicit Neighbor(const Rcpp::NumericVector &);
  Neighbor(const Neighbor &);
  Neighbor(Neighbor &&);
  ~Neighbor();

  void set_neighbors(const arma::uvec &);
  void set_neighbors(arma::uvec &&);
  void set_neighbors(const Rcpp::NumericVector &);

  size_t get_size() const;
  const arma::uvec &get_neighbors() const;
};

Neighbor::Neighbor(const arma::uvec &v) {
  if (v.size() > 0) {
    this->neighbors = v;
  } else {
    this->neighbors.reset();
  }
}

Neighbor::Neighbor(arma::uvec &&v) : neighbors(std::move(v)) {}

Neighbor::Neighbor(const Rcpp::NumericVector &v) {
  if (v.size() > 0) {
    this->neighbors = arma::uvec(v.size());

    for (auto i = 0; i != v.size(); i++) {
      this->neighbors(i) = v[i];
    }
  } else {
    this->neighbors.reset();
  }
}

Neighbor::Neighbor(const Neighbor &v) { this->neighbors = v.neighbors; }

Neighbor::Neighbor(Neighbor &&v) : neighbors(std::move(v.neighbors)) {}

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
  this->neighbors(std::move(v));
}

void
Neighbor::set_neighbors(const Rcpp::NumericVector &v) {
  if (v.size() > 0) {
    this->neighbors = arma::uvec(v.size());

    for (auto i = 0; i != v.size(); i++) {
      this->neighbors(i) = v[i];
    }
  } else {
    this->neighbors.reset();
  }
}

size_t
Neighbor::get_size() const {
  return this->neighbors.size();
}

const arma::uvec &
Neighbor::get_neighbors() const {
  return this->neighbors;
}

#endif   // BAYESSPACE_NEIGHBOR_H
