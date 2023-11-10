#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"


int main(){
  
  double alpha;
  int n=1;
  double x=0.5;
  
  double y=adomain_v2(h,x,alpha,n);
  
  printf("%lf\n",y);
  
  return 0;
  }
  

