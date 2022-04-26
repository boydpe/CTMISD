#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_h(NumericMatrix p0, 
                       NumericMatrix dist_mat,
                       NumericMatrix dist_ints,
                       int n_ints) {
  
  // numerator
  int n_l = n_ints;
  int n_i = pow(p0.size(), 0.5);
  NumericVector num_h(n_l);
  NumericVector den_h(n_l);
  double diag_sum = 0;
  
  if( sum(dist_mat) == 0){
    for(int i = 0; i < n_l; i++) {
      dist_ints(i,2) = 1;
    }
  }
  
  else {
    
    for (int l = 0; l < n_l; l++) {
      for (int i = 0; i < n_i; i++){
        for (int j = i; j < n_i; j++){
          if(dist_ints(l,0) < dist_mat(j,i) && 
             dist_mat(j,i) <= dist_ints(l,1)){
            num_h[l] += p0(j,i);
          } else{
            num_h[l] += 0;
          }
        }
      }
    }
    
    // denominator
    
    for (int i = 0; i < n_l; i++){
      den_h[i] = (dist_ints(i,1) - dist_ints(i,0));
    }
    
    for (int i = 0; i < n_i; i++) {
      diag_sum += p0(i,i);
    }
    den_h = den_h * (n_i - diag_sum);
    for (int i = 0; i < n_l; i++) {
      dist_ints(i,2) = num_h[i] / den_h[i];
    }
  }
  return dist_ints;
}
