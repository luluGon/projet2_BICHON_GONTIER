#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"

int main(){
  
  int n=4;
  double y;
  double t;
  t=Uex(0.2, H);
  y=adomain_v1(h,0.2,10.0,n);
  printf("%lf\n%lf\n",y,t);
  return 0;
  }
  

