#include "OperationsMatrices.h"

//Dans ce fichier, on écrit les programmes de calcul pour les opérations sur les matrices

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
void sous_matrice(double In[N1+1][N1+1],double A[N1+1][N1+1]){
  int i,j;
  for (i=0;i<=N1;i++){
    for (j=0;j<=N1;j++){
      A[i][j]=In[i][j]-A[i][j];
    }
  }
}
