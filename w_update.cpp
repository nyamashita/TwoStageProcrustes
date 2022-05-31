#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double posit(double val){
  double ret;
  if (val > 0){
    ret = val;
  }else{
    ret = 0;
  }
  return ret;
}

// [[Rcpp::export]]
mat W_update(mat X, mat Y, double lambda){
  mat W_return(X.n_cols, Y.n_cols, fill::zeros);
  mat W_ols = inv(trans(X)*X)*trans(X)*Y;
  for (uword i=0; i <= W_ols.n_rows-1; ++i){
    for (uword j=0; j <= W_ols.n_cols-1; ++j){
      W_return(i,j) = sign(W_ols(i,j)) * posit(abs(W_ols(i,j)) - lambda/(2*accu(pow(X.col(i), 2))));
    }
  }
  return W_return;
}
