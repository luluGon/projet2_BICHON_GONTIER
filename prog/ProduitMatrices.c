#include "ProduitMatrice.h"

double prod_mat(double A[N1+1,N1+1],double B[N1+1]){
  int i,j;
  double C[N1+1];
  for (i=0;i<=N1;i++){
    C[i]=0.;
    for (j=0;j<=N1;j++){
      C[i]+=A[i,j]*B[j];
    }
  }
  return C
}
