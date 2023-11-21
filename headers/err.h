#ifndef err
#define err
  #include "constante.h"
  #include "methodes.h"
  #include "Gauss_Legendre.h"
  #include "fonctions.h"
  
  
  double L2_err_f3_adov1(double alpha,int n);
  
  double L2_err_f3_approx_succesive(double (*f)(double),double alpha,int n);
#endif
