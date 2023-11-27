#include "coefficients.h"
#include <math.h>

/*Les fonctions suivantes calculent les coefficient dans le cas général, sans analyse qui permettrait de d'économiser des calculs */


//GLOBAL


double A0 (){
  
  /*On initie la variable locale de la fonction*/
  double res;
  
  
  /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à f0 qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
  res= integrale_GL(f0,0.0,L,0);
  res=res/L;
  
  return res;
  }




//ADOMAIN V1

//calcule le coefficient B0 avec f3 calculé par la première version de la méthode de Adomain
double B0_ado1 ( double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
 
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à f3-f0 qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

  double f3_bis(double x1){
      return adomain_v1(x1, alpha,n);
    } 

/*On décompose l'intégrale par linéarité pour appliquer GL sur chacune des fonctions */
  inte1= integrale_GL(f3_bis,0.0,L,0);
  inte2=integrale_GL(f0,0.0,L,0);
  
  res=inte1-inte2;
  res=res/(H*L);
  
  return res;
  }

//de même pour Am
double Am_ado1 ( int m, double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double scal;
  double inte1;
  double inte2;
  
  scal=m*pi/L;
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (f3-exp(-m*pi*H/L)*f0)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

    
   double f3_bis(double x1){
      return adomain_v1(x1, alpha,n);
    } 

 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(-scal*H)*inte2;
  
  res=inte1-inte2;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }

//de même pour Bm
double Bm_ado1(int m,double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
  double scal;
  
  scal=m*pi/L;
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (exp(m*pi*H/L)*f0-f3)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

    
   double f3_bis(double x1){
      return adomain_v1(x1, alpha,n);
    } 
 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(scal*H)*inte2;
  
  res=inte2-inte1;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }

//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par adomain v1
void gen_tab_ado1(double A_m[],double B_m[],double alpha,int n){
  
  int i;
  
  for (i=0; i<N2; i++){
    A_m[i]=Am_ado1(i+1,alpha,n);
    
    B_m[i]=Bm_ado1(i+1,alpha,n);
    }
  }



//ADOMAIN V2

//calcule B0 avec f3 par la seconde version de Adomain
double B0_ado2 ( double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
 
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à f3-f0 qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

  double f3_bis(double x1){
      return adomain_v2(x1, alpha,n);
    } 

/*On décompose l'intégrale par linéarité pour appliquer GL sur chacune des fonctions */
  inte1= integrale_GL(f3_bis,0.0,L,0);
  inte2=integrale_GL(f0,0.0,L,0);
  
  res=inte1-inte2;
  res=res/(H*L);
  
  return res;
  }

//de même pour Am
double Am_ado2 ( int m, double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double scal;
  double inte1;
  double inte2;
  
  scal=m*pi/L;
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (f3-exp(-m*pi*H/L)*f0)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

    
   double f3_bis(double x1){
      return adomain_v2(x1, alpha,n);
    } 

 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(-scal*H)*inte2;
  
  res=inte1-inte2;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }

//de même pour Bm
double Bm_ado2(int m,double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
  double scal;
  
  scal=m*pi/L;
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (exp(m*pi*H/L)*f0-f3)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

    
   double f3_bis(double x1){
      return adomain_v2(x1, alpha,n);
    } 
 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(scal*H)*inte2;
  
  res=inte2-inte1;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }

//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par adomain v2
void gen_tab_ado2(double A_m[],double B_m[],double alpha,int n){
  
  int i;
  
  for (i=0; i<N2; i++){
    A_m[i]=Am_ado2(i+1,alpha,n);
    
    B_m[i]=Bm_ado2(i+1,alpha,n);
    }
  }  
  
  
 
//APPROXIMATIONS SUCCESIVES

//calcule B0 avec f3 par approximations succesives
double B0_AS (double (*f)(double), double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
 
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à f3-f0 qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

  double f3_bis(double x1){
      return approximation_succesive(f,x1, alpha,n);
    } 

/*On décompose l'intégrale par linéarité pour appliquer GL sur chacune des fonctions */
  inte1= integrale_GL(f3_bis,0.0,L,0);
  inte2=integrale_GL(f0,0.0,L,0);
  
  res=inte1-inte2;
  res=res/(H*L);
  
  return res;
  }

//de même pour Am
double Am_AS ( double (*f)(double),int m, double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double scal;
  double inte1;
  double inte2;
  
  scal=m*pi/L;
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (f3-exp(-m*pi*H/L)*f0)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 
 
    
   double f3_bis(double x1){
      return approximation_succesive(f,x1, alpha,n);
    } 

 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(-scal*H)*inte2;
  
  res=inte1-inte2;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }


//de même pour Bm
double Bm_AS(double (*f)(double),int m,double alpha,int n){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
  double scal;
  
  scal=m*pi/L;
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (exp(m*pi*H/L)*f0-f3)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

   double f3_bis(double x1){
      return approximation_succesive(f,x1, alpha,n);
    } 
 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(scal*H)*inte2;
  
  res=inte2-inte1;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }
  
//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par approximations succesives 
void gen_tab_AS(double A_m[],double B_m[],double (*f)(double),double alpha,int n){
  
  int i;
  
  for (i=0; i<N2; i++){
    A_m[i]=Am_AS(f,i+1,alpha,n);
    
    B_m[i]=Bm_AS(f,i+1,alpha,n);
    }
  }   
  
  
//NOYAUX SÉPARABLES

//calcul de B0 avec f3 par la méthode de noyaux séparables
double B0_KS ( double alpha){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
 
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à f3-f0 qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

  double f3_bis(double x1){
      return ker_sep(x1, alpha);
    } 

/*On décompose l'intégrale par linéarité pour appliquer GL sur chacune des fonctions */
  inte1= integrale_GL(f3_bis,0.0,L,0);
  inte2=integrale_GL(f0,0.0,L,0);
  
  res=inte1-inte2;
  res=res/(H*L);
  
  return res;
  }

//de même pour Am
double Am_KS ( int m, double alpha){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double scal;
  double inte1;
  double inte2;
  
  scal=m*pi/L;
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (f3-exp(-m*pi*H/L)*f0)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 
    
   double f3_bis(double x1){
      return ker_sep(x1, alpha);
    }  

 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(-scal*H)*inte2;
  
  res=inte1-inte2;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }


//de même pour Bm
double Bm_KS(int m,double alpha){
  
  /*On initie la variable locale de la fonction*/
  double res;
  double inte1;
  double inte2;
  double scal;
  
  scal=m*pi/L;
    
 /* On calcule l'intégrale grâce à la fonction intégrale_GL appliqué à (exp(m*pi*H/L)*f0-f3)cos(m*pi*x/L) qui utilise la quadrature de Gauss-Legendre, cette fonction se trouve dans sous-programme Gauss_Legendre*/
 

    double f3_bis(double x1){
      return ker_sep(x1, alpha);
    }  
 /*On va utiliser la linearité de l'integrale pour calculer l'integrale*/
  
  inte1= integrale_GL(f3_bis,0.0,L,m);
  inte2=integrale_GL(f0,0.0,L,m);
  inte2=exp(scal*H)*inte2;
  
  res=inte2-inte1;
  
  res=res/(L*sinh(scal*H));
  
  return res;
  }
  
//initialise deux tableaux prenant les coefficients Am et Bm en valeur, pour f3 par noyaux séparables
void gen_tab_KS(double A_m[],double B_m[],double alpha){
  
  int i;
  
  for (i=0; i<N2; i++){
    A_m[i]=Am_KS(i+1,alpha);
    
    B_m[i]=Bm_KS(i+1,alpha);
    }
  }  





