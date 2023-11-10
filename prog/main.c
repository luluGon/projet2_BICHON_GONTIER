#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"


int main(){
  
  double alpha=5;
  int n=100;
  double x=0.2;
  double t=Uex(x,1);
  double y=adomain_v2(h,x,alpha,n);
  printf("test =%lf\n",y);
  printf("exact =%lf\n",t);
  return 0;
  }
  

