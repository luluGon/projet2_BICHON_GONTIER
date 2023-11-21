#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"

int main(){
  
  int n=5;
  double alpha_opti;
  
  double alpha_deb=0.1;
  double alpha_fin=100.0;
  double pas =0.05;
  double f(double x){
   return 0;}
   
  alpha_opti=plot_err_alpha_approx_succesive(f,alpha_deb, alpha_fin,pas,n);
  
  printf("%lf\n",alpha_opti);
  return 0;
  }
  

