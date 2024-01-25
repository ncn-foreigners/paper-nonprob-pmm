#include <Rcpp.h>
#include <random>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector syspss_cpp(NumericVector x, int n) {
  int N = x.size();
  IntegerVector U = Rcpp::seq(1,N);
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(U.begin(), U.end(), g);
  
  NumericVector xx(N);
  for(int i=0; i<N; i++) {
    xx[i] = x[U[i]-1];
  }
  NumericVector xx_cs = cumsum(x);
  NumericVector z = n*xx_cs/sum(x);
  
  double r = runif(1)[0];
  NumericVector s(n);
  int k = 0;
  for(int i=0; i<N; i++) {
    if(z[i] >= r) {
      s[k] = U[i];
      r = r + 1;
      k = k + 1;
    }
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix syspss_pij_cpp(NumericVector x, NumericVector s, Int32 d) {
  int n = s.size();
  NumericMatrix p(n, n);
  std::fill( p.begin(), p.end(), 0);
  for (int k = 0; k < d; ++k) {
    NumericVector ss = syspss_cpp(x, n);
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (min(abs(ss - s[i])) + min(abs(ss - s[j])) == 0) {
          p(i,j) += 1;
        }
      }
    }
  }
  return p;
}

/*** R
set.seed(123)
x <- rnorm(2000)
z <- x + max(x)
syspss_cpp(z, 10)
*/


