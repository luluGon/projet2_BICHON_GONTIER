#ifndef newton
#define newton

#include <stdio.h>
#include <math.h>
#include "fonctions.h"
#include "T_alpha.h"

double newton_T_alpha_ado1(double x,double y_init, double err, int ite_max,double alpha,int n) ;

double newton_T_alpha_ado2(double x,double y_init, double err, int ite_max,double alpha,int n) ;

double newton_T_alpha_AS(double x,double y_init, double err, int ite_max,double (*f)(double),double alpha,int n);

double newton_T_alpha_KS(double x,double y_init, double err, int ite_max,double alpha) ;


#endif
