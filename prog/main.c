#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"

int main(){
  
  int n=11;
  double t1,t2;
  double alpha_opt_adov1=0.272000;
  double alpha=1.0;
  double alpha_opt_app_succ=1.632200;

  int nx=100;
  
  double f(double x){
  return 0.0;}
  
  //alpha_opt_app_succ=plot_err_alpha_approx_succesive(f,1.0,3.0,0.0001,n);
  //printf("%lf\n",alpha_opt_app_succ);
  //t1=adomain_v1(0.0, alpha,n);
  //printf("%lf\n%lf\n", h(0.0)/alpha,t1);

  //plot_f3_app_succ(f, nx, n,alpha_opt_app_succ);
  precalcul_KUn(10);
  
  return 0;
  }
  

