#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_g(NumericMatrix p0, 
                       NumericMatrix  time_mat, 
                       NumericMatrix time_ints,
                       int n_ints) {
  
  // numerator
  int n_l = n_ints;
  int n_i = pow(p0.size(), 0.5);
  NumericVector num_g(n_l);
  NumericVector den_g(n_l);
  double diag_sum = 0;
  NumericMatrix g_vals(n_l, 2);
  
  for (int l = 0; l < n_l; l++) {
    for (int i = 0; i < n_i; i++){
      for (int j = i; j < n_i; j++){
        if(time_ints(l,0) < fabs(time_mat(j,i)) && 
           fabs(time_mat(j,i)) <= time_ints(l,1)){
          num_g[l] += p0(j,i);
        } else{
          num_g[l] += 0;
        }
      }
    }
  }
  
  // denominator
  
  for (int i = 0; i < n_l; i++){
    den_g[i] = (time_ints(i,1) - time_ints(i,0));
  }
  
  for (int i = 0; i < n_i; i++) {
    diag_sum += p0(i,i);
  }
  
  den_g = den_g * (sum(p0) - diag_sum);
  // for (int i = 0; i < n_l; i++) {
  //   time_ints(i,2) = num_g[i] / den_g[i];
  // }
  for (int i = 0; i < n_l; i++) {
    g_vals(i,0) = num_g[i];
    g_vals(i,1) = den_g[i];
  }
  return g_vals;
}

