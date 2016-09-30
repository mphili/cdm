#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP calc_lxa(SEXP Xs, SEXP Ps) {
BEGIN_RCPP

  NumericMatrix Xr(Xs);
  NumericMatrix Pr(Ps);

  int n = Xr.nrow(), j = Xr.ncol(), l = Pr.ncol();
  
  mat x(Xr.begin(), n, j, false);
  mat pa(Pr.begin(), j, l, false);
  
  mat post = exp(x * log(pa) + (1-x) * log(1-pa));
  return wrap(post);
  
END_RCPP
}

RcppExport SEXP calc_lxa2(SEXP data, SEXP Palpha) {
BEGIN_RCPP

  NumericMatrix x(data);
  NumericMatrix pa(Palpha);

  unsigned int I = x.nrow();
  unsigned int J = x.ncol();
  unsigned int L = pa.ncol();

  NumericMatrix post(I, L);
  std::fill(post.begin(), post.end(), 1.0);
     
  for (int i = 0; i < I; i++) {
    for (int l = 0; l < L; l++) {
      for (int j = 0; j < J; j++) {
        post(i,l) *= ((x(i,j) == 1) ? pa(j,l) : (1-pa(j,l)));
      }
    }
  }
     
  return post;

END_RCPP
}
