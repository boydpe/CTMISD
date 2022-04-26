#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_theta(NumericVector min_int, 
                        NumericVector max_int,
                        NumericMatrix diff_mat, 
                        NumericMatrix p0,
                        bool mark_check) {
  
  int n = pow(p0.size(), 0.5);
  int m = min_int.size();
  NumericVector theta(m);
  
  if(mark_check == 0) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < j; k++) {
        if (min_int[i] < diff_mat(j,k) && 
            diff_mat(j,k) <= max_int[i]) {
          theta[i] = theta[i] + p0(j,k);
        }
      }
    }
  }
  return theta;
  } else {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
          if (min_int[i] < diff_mat[j] && 
              diff_mat[j] <= max_int[i]) {
            theta[i] = theta[i] + p0(j,k);
          }
        }
      }
    }
    return theta;
  }
}


