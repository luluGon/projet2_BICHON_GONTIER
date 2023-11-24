#ifndef plot
#define plot

  #include "err.h"
  #include "constante.h"
  
  double plot_err_alpha_adov1(double alpha_deb,double alpha_fin, double pas,int n);
  
  double plot_err_alpha_approx_succesive(double (*f)(double),double alpha_deb,double alpha_fin, double pas,int n);
  
  int plot_f3_adov1( int nx,int n,double alpha);
  
  int plot_f3ex( int nx);
  
  int plot_f3_app_succ(double (*f)(double), int nx,int n,double alpha);
#endif

