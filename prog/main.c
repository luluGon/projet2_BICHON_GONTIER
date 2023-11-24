#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"

int main(){
  int n=2;
  double t1,t2;
  double alpha_opt_adov1=0.272000;
  double alpha=1.0;
  double alpha_opt_app_succ=1.632200;
  double x=0.3;
  int nx=100;
  
  double f(double x){
  return 0.0;}
  
  //alpha_opt_app_succ=plot_err_alpha_approx_succesive(f,0.1,1,0.1,n);
  alpha_opt_adov1=plot_err_alpha_adov1(0.1,10,0.01,n);
  //printf("%lf\n",alpha_opt_app_succ);
  //t1=adomain_v1(x, alpha_opt_adov1,n);
  //t1=adomain_v1_2(x, Un, alpha_opt_adov1,n);
  //printf("%lf\n%lf\n", Uex(x,H),t1);
  //double tab[10];
  //tab_sin(tab);
  //printf("%lf\n", GL2(tab, 0.0,L, 1));

  //plot_f3_app_succ(f, nx, n,alpha_opt_app_succ);
  //precalcul_KUn(10);
  
  return 0;
  }
  

