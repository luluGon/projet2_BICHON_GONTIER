#include "err.h"

double L2_err_f3(double alpha, int n){

    double res1;

    double diff(double x){
      double res;
      res=adomain_v2(h,x,alpha,n)- Uex(x,H);
      return res*res;
      }
    res1=integrale_GL(diff, 0.0,L,0);

    return sqrt(res1);
}
