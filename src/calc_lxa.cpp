#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP calc_lxa(SEXP data, SEXP Palpha) {
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
