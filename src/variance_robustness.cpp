#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
unsigned nChoosek( unsigned n, unsigned k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}


//' @title Compute the Conditional Variance of the Robustness
//'
//' @description Compute Var_A(E_S[AUC | A])
//' @inheritParams auc_robustness_lbm
//' @export
// [[Rcpp::export]]
double var_auc_unif_lbm_cpp(NumericMatrix con, NumericVector pi,
                            NumericVector rho, int nr, int nc) {
  int K = pi.length();
  int Q = rho.length();
  double res = 0;
  double aucsq = 0;
  double tmp = 0;
  NumericVector etaq(Q);
  NumericMatrix etaqq(Q, Q);
  for(int q = 0;  q < Q; q++) {
    for(int k = 0;  k < K; k++) {
      etaq(q) += pi(k)*(1-con(k,q));
    }
  }
  for(int q = 0;  q < Q; q++) {
    for(int qp = 0; qp < Q; qp++) {
      for(int k = 0;  k < K; k++) {
        etaqq(q, qp) += pi(k)*(1-con(k,q))*(1-con(k, qp));
      }
    }
  }
  for (int m = 0; m <= nr; m++) {
      for (int q = 0; q < Q; q++) {
        aucsq += pow(etaq(q), m) * rho(q);
      }
  }
  aucsq = aucsq*aucsq/pow(nr, 2);
  for (int m = 0; m <= nr; m++) {
    for (int mp = 0; mp <= nr; mp++) {
      for (int l = std::max(m, mp); l <= std::min(nr, m + mp); l++ ) {
        for(int q = 0;  q < Q; q++) {
          for(int qp = 0; qp < Q; qp++) {
            tmp  += (1-1/double(nc)) *
                  rho(q) * rho(qp) * double(pow(etaq(q), l-mp)) *
                  double(pow(etaq(qp),l-m)) *
                  double(pow(etaqq(q, qp), m+mp-l));
            // Rcout << pow(etaqq(q, qp),m+mp-l) << "\t";
          }
          tmp += (1/double(nc)) * double(pow(etaq(q),l)) * rho(q);
        }
        res = res + tmp * nChoosek(m, l-mp) * nChoosek(nr-m, l-m) /
          double(nChoosek(nr, mp));
        tmp = 0;
      }
    }
  }
  res = res/double(nr*nr) - aucsq;
  return res;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
pi <- c(.3,.7)
rho <- c(.3,.7)
con <- matrix(c(.5, .3, .3, .1), 2, 2)
nr <- 50
nc <- 5
var_auc_unif_lbm_cpp(con = con, pi =  pi,
                     rho =  rho, n =  nr, nc =  nc)
*/
