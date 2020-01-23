#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp;


// **********************************************************//
//     	           Likelihood function            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull3 (arma::vec Y,
                    arma::vec Y0,
                    arma::vec eXB,
                    arma::vec delta,
                    arma::vec C,
                    double rho) {
  arma::vec dexp0 = exp(-pow(eXB % Y0, rho));
  arma::vec dexp1 = exp(-pow(eXB % Y, rho));
  arma::vec dexp2 = pow(eXB % Y, rho - 1);
  arma::vec llik1 = delta + (1 - delta) % dexp1;
  arma::vec ldelta = log(1 - delta);
  arma::vec leXB = log(eXB);
  arma::uvec ids5 = find(dexp0 == 0);
  dexp0.elem(ids5).fill(exp(-740));
  arma::uvec ids0 = find(llik1 == 0);
  llik1.elem(ids0).fill(exp(-740));
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp2 == 0);
  dexp2.elem(ids2).fill(exp(-740));
  arma::uvec ids3 = find(ldelta == -arma::datum::inf);
  ldelta.elem(ids3).fill(-740);
  arma::uvec ids4 = find(leXB == -arma::datum::inf);
  leXB.elem(ids4).fill(-740);
  arma::vec llik = C % (ldelta + log(rho) + leXB + log(dexp2) - pow(eXB % Y, rho)) + (1 - C) % log(llik1) -log(dexp0);
  return sum(llik);
}
