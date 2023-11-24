#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"
#include "InversionMatrice.h"

int main(){
  

  int n=100;
  double y,z;
  y=plot_err_alpha_adov2(0.01,10.,0.01,n);
  printf("alpha_min pour adomainv2= %lf\n", y);
  
  z=plot_err_alpha_ker_sep(0.000001,0.001,0.000001);
  printf("alpha_min pour ker_sep= %lf\n", z);
  
  double A[2][2];
  
  A[0][0]=-1.;
  A[0][1]=0.;
  A[1][0]=-11.;
  A[1][1]=10.;
  inv_mat(A);
  
  printf("A= %lf\n",A[0][0]);
  printf("A= %lf\n",A[0][1]);
  printf("A= %lf\n",A[1][0]);
  printf("A= %lf\n",A[1][1]);
  return 0;
  }
  

