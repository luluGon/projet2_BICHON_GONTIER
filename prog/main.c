#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"

int main(){
  
  int n=100;
  double y;
  y=plot_err_alpha_adov1(0.1,10.0,0.1,n);
  printf("alpha_min = %lf\n", y);
  return 0;
  }
  

