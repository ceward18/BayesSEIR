#include <Rcpp.h>

//' @useDynLib BayesSEIR, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::NumericVector getS(int S0, Rcpp::NumericVector Estar) {
  
  int N = Estar.length();
  
  Rcpp::NumericVector out(N);
  int cumsum = 0;
  
  for (int i = 0; i < N; i++) {
    
    out[i] = S0 - cumsum;
    cumsum = cumsum + Estar[i];
    
  }
  
  return out;
}

//' @useDynLib BayesSEIR, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]] 
Rcpp::NumericVector getE(int E0, Rcpp::NumericVector Estar, Rcpp::NumericVector Istar) {
  
  int N = Estar.length();
  
  Rcpp::NumericVector out(N);
  int cumsum = 0;
  
  for (int i = 0; i < N; i++) {
    
    out[i] = E0 - cumsum;
    cumsum = cumsum - Estar[i] + Istar[i];
    
  }
  
  return out;
}

//' @useDynLib BayesSEIR, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]] 
Rcpp::NumericVector getI(int I0, Rcpp::NumericVector Istar, Rcpp::NumericVector Rstar) {
  
  int N = Istar.length();
  
  Rcpp::NumericVector out(N);
  int cumsum = 0;
  
  for (int i = 0; i < N; i++) {
    
    out[i] = I0 - cumsum;
    cumsum = cumsum - Istar[i] + Rstar[i];
    
  }
  
  return out;
}
