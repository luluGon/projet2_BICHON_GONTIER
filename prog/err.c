#include "err.h"

double L2_err_f3_adov1(double alpha, int n){

    double res1;

    double diff(double x){
      double res;
      res=adomain_v1(x,alpha,n)- Uex(x,H);
      return res*res;
      }
    res1=integrale_GL(diff, 0.0,L,0);

    return sqrt(res1);
}

double L2_err_f3_approx_succesive(double (*f)(double),double alpha,int n){

  double res1;

    double diff(double x){
      double res;
      res=approximation_succesive(f,x,alpha,n)- Uex(x,H);
      return res*res;
      }
    res1=integrale_GL(diff, 0.0,L,0);

    return sqrt(res1);
}
