//Header du fichier calculant les erreurs L2 sur [0,L]x{H} pour chaque méthode

#ifndef err
#define err
  #include "constante.h"
  #include "methodes.h"
  #include "Gauss_Legendre.h"
  #include "fonctions.h"
  
  
  double L2_err_f3_adov1(double alpha,int n);
  
  double L2_err_f3_adov2(double alpha, int n);
  
  double L2_err_f3_approx_succesive(double (*f)(double),double alpha,int n);
  
  double L2_err_ker_sep_f3(double alpha);
#endif
