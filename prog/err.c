#include "err.h"


//Dans ce fichier on calcule l'erreur L2 sur l'intervalle [0,L]x{H}, associée à un paramètre alpha pour chaque méthode
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

double L2_err_f3_adov2(double alpha, int n){

    double res1;

    double diff(double x){
      double res;
      res=adomain_v2(x,alpha,n)- Uex(x,H);
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

double L2_err_ker_sep_f3(double alpha){

    double res1;

    double diff(double x){
      double res;
      res=ker_sep(x,alpha)- Uex(x,H);
      return res*res;
      }
    res1=integrale_GL(diff, 0.0,L,0);

    return sqrt(res1);
}
