#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix update_p(NumericMatrix p0, NumericMatrix dist_mat,
                       NumericVector br, NumericVector predg, NumericVector predh,
                       NumericVector predk, int sp){
  
  int n_i = pow(p0.size(), 0.5);
  double num_p;
  double den_p ;
  NumericMatrix p(n_i, n_i);
  double z;
  // const double pi = 3.14159265358979323846;
  double pi = M_PI;
  double sum_mat;
  
  // for (int i = 1; i < n_i; i++) {
  //   for ( int j = 0; j < i; j++){
  //     sum_mat += dist_mat(i,j);
  //   }
  // }
  
  if ( sp == 1){
    for (int i = 0; i < n_i; i++){
      for (int j = 0; j < n_i; j++){
        if  ( j < i){
          num_p = (predg[i*(i+1)/2 + j]*predh[i*(i+1)/2 + j]*predk[j]);// /
            // (2*pi*dist_mat(i,j)));
          
          for (int  l = 0; l < i; l++){
            z += predg[i*(i+1)/2 + l]*predh[i*(i+1)/2 + l]*predk[l];// /
             //  (2*pi*dist_mat(i,l));
          }
          den_p = br[i] + z; 
          p(i,j) = num_p / den_p;
          z = 0;
          
        } else if (i == j) { 
          num_p = br[i];
          
          for (int  l = 0; l < i; l++){
            z += predg[i*(i+1)/2 + l]*predh[i*(i+1)/2 + l]*predk[l];// /
             // (2*pi*dist_mat(i,l));
          }
          den_p = br[i] + z;
          p(i,j) = num_p / den_p;
          z = 0;
        } else {
          p(i,j) = 0;
        }
      }
    }
    p(0,0) = 1;
    return p;
  } else {
    
    for (int i = 0; i < n_i; i++){
      for (int j = 0; j < n_i; j++){
        if  ( j < i){
          num_p = predg[i*(i+1)/2 + j]*predh[i*(i+1)/2 + j]*predk[j];
          
          for (int  l = 0; l < i; l++){
            z += predg[i*(i+1)/2 + l]*predh[i*(i+1)/2 + l]*predk[l];
          }
          den_p = br[i] +z;
          p(i,j) = num_p / den_p;
          z = 0;
          
        }else if (i == j) { 
          num_p = br[i];
          
          for (int  l = 0; l < i; l++){
            z += predg[i*(i+1)/2 + l]*predh[i*(i+1)/2 + l]*predk[l];
          }
          den_p = br[i] + z;
          p(i,j) = num_p / den_p;
          z = 0;
        }
        
        else {
          p(i,j) = 0;
        }
      }
    }
    p(0,0) = 1;
    return p;
  }
}
