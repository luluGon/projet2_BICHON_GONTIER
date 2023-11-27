//Dans ce headers, on retrouve toutes nos fonctions pour calculer les coefficients,
//selon la méthode employé pour le calcul de f3
#ifndef coefficients
#define coefficients

#include "fonctions.h"
#include "methodes.h"
#include <math.h>

//GLOBAL
/*fonction calculant le coefficient A0*/
double A0 ();


//ADOMAIN V1

//calcule le coefficient B0 avec f3 calculé par la première version de la méthode de Adomain
double B0_ado1 ( double alpha,int n);

//de même pour Am
double Am_adov1 ( int m, double alpha,int n);

//de même pour Bm
double Bm_ado1(int m,double alpha,int n);


//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par adomain v1
void gen_tab_ado1(double A_m[],double B_m[],double alpha,int n);


//ADOMAIN V2

//calcule B0 avec f3 par la seconde version de Adomain
double B0_ado2 ( double alpha,int n);

//de même pour Am
double Am_adov2 ( int m, double alpha,int n);

//de même pour Bm
double Bm_ado2(int m,double alpha,int n);

//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par adomain v2
void gen_tab_ado2(double A_m[],double B_m[],double alpha,int n);


//APPROXIMATION SUCCESIVESES

//calcule B0 avec f3 par approximations succesives
double B0_AS (double (*f)(double), double alpha,int n);

//de même pour Am
double Am_AS ( double (*f)(double),int m, double alpha,int n);

//de même pour Bm
double Bm_AS(double (*f)(double),int m,double alpha,int n);

//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par approximations succesives 
void gen_tab_AS(double A_m[],double B_m[],double (*f)(double),double alpha,int n);


//NOYAUX SÉPARABLES

//calcul de B0 avec f3 par la méthode de noyaux séparables
double B0_KS ( double alpha);

//de même pour Am
double Am_KS ( int m, double alpha);

//de même pour Bm
double Bm_KS(int m,double alpha);

//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par noyaux séparables
void gen_tab_KS(double A_m[],double B_m[],double alpha);


#endif
