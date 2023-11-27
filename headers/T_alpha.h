#ifndef T_alpha
#define T_alpha
#include "coefficients.h"
#include <math.h>

//ADOMAIN V1

//Calcul de T_alpha(x,y) pour f3 par adomain v1
double Talpha_ado1(double x, double y,double alpha,int n);

//T'(x,y) pour adomain v1
double T_prime_ado1(double x, double y, double alpha, int n);


//ADOMAIN V2

//Calcul de T_alpha(x,y) pour f3 par adomain v2
double Talpha_ado2(double x, double y,double alpha,int n);

//T'(x,y) pour adomain v2
double T_prime_ado2(double x, double y, double alpha, int n);


//APPROXIMATION SUCCESIVES

//Calcul de T_alpha(x,y) pour f3 par approximation succesives
double Talpha_AS(double x, double y,double (*f)(double),double alpha,int n);

//T'(x,y) pour approxiamtion succesives
double T_prime_AS(double x, double y,double (*f)(double), double alpha, int n);


//NOYAUX SÉPARABLES

//Calcul de T_alpha(x,y) pour f3 par anoyaux séprables
double Talpha_KS(double x, double y,double alpha);

//T'(x,y) pour noyaux séparables
double T_prime_KS(double x, double y, double alpha);


#endif
