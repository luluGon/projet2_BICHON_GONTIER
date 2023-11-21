#include "OpérationsMatrices.h"

//On écrit un programme permettant d'effectuer un produit matriciel
void prod_mat(double A[N1+1][N1+1],double B[N1+1]){
  int i,j;
  double C[N1+1];
  for (i=0;i<=N1;i++){
    C[i]=0.;
    for (j=0;j<=N1;j++){
      C[i]+=A[i][j]*B[j];
    }
    B[i]=C[i];
  }
}

//On écrit un programme permettant d'effectuer une multiplication matrice/scalaire
void mult_mat_scal(double lambda,double A[N1+1][N1+1]){
  int i,j;
  for (i=0;i<=N1;i++){
    for (j=0;j<=N1;j++){
      A[i][j]*=lambda;
    }
  }
}

//On écrit un programme permettant d'effectuer une soustraction de deux matrices
void sous_matrice(double I[N1+1][N1+1],double A[N1+1][N1+1]){
  int i,j;
  for (i=0;i<=N1;i++){
    for (j=0;j<=N1;j++){
      A[i][j]=I[i][j]-A[i][j];
    }
  }
}
