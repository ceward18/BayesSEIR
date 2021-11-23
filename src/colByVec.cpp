// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @useDynLib BayesSEIR, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]] 
arma::vec colByVec(arma::mat x, arma::vec y) {
  
  int I = x.n_rows;
  int J = x.n_cols;
  
  
  arma::vec out = arma::zeros(I);
  
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      
      out(i) += x(i, j) * y(j);
      
      
    }
  }
  
  return out;
}