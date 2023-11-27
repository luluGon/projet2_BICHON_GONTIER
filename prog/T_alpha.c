#include "T_alpha.h"

//ADOMAIN V1

//Calcul de T_alpha(x,y) pour f3 par adomain v1
double Talpha_ado1(double x, double y,double alpha,int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_ado1(A_m, B_m,  alpha, n);
    
  res=A0()+B0_ado1(alpha,n)*y;
  
  for (i=1; i<=N2;i++){
  
    scal=i*pi/L;
    res=res + (A_m[i-1]*exp(scal*y)+B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
    
  return res;
  
  }
  
//T'(x,y) pour adomain v1
double T_prime_ado1(double x, double y, double alpha, int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_ado1(A_m, B_m,alpha, n);
  
  res= B0_ado1(alpha,n);
  
  for (i=1;i<=N2;i++){
  
    scal=i*pi/L;
    res=res +(scal*A_m[i-1]*exp(scal*y) -scal*B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
  return res;
  }


//ADOMAIN V2

//Calcul de T_alpha(x,y) pour f3 par adomain v2
double Talpha_ado2(double x, double y,double alpha,int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_ado2(A_m, B_m,alpha, n);
    
  res=A0()+B0_ado2(alpha,n)*y;
  
  for (i=1; i<=N2;i++){
  
    scal=i*pi/L;
    res=res + (A_m[i-1]*exp(scal*y)+B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
    
  return res;
  
  }
 
 
//T'(x,y) pour adomain v2
double T_prime_ado2(double x, double y, double alpha, int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_ado2(A_m, B_m,alpha, n);
  
  res= B0_ado2(alpha,n);
  
  for (i=1;i<=N2;i++){
  
    scal=i*pi/L;
    res=res +(scal*A_m[i-1]*exp(scal*y) -scal*B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
  return res;
  }

//APPROXIMATION SUCCESIVES

//Calcul de T_alpha(x,y) pour f3 par approximation succesives
double Talpha_AS(double x, double y,double (*f)(double),double alpha,int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_AS(A_m, B_m,f,alpha, n);
    
  res=A0()+B0_AS(f,alpha,n)*y;
  
  for (i=1; i<=N2;i++){
  
    scal=i*pi/L;
    res=res + (A_m[i-1]*exp(scal*y)+B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
    
  return res;
  
  }
 
 
//T'(x,y) pour approxiamtion succesives
double T_prime_AS(double x, double y,double (*f)(double), double alpha, int n){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_AS(A_m, B_m,f,alpha, n);
  
  res= B0_AS(f,alpha,n);
  
  for (i=1;i<=N2;i++){
  
    scal=i*pi/L;
    res=res +(scal*A_m[i-1]*exp(scal*y) -scal*B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
  return res;
  }
  

//NOYAUX SÉPARABLES

//Calcul de T_alpha(x,y) pour f3 par anoyaux séprables
double Talpha_KS(double x, double y,double alpha){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_KS(A_m, B_m,alpha);
    
  res=A0()+B0_KS(alpha)*y;
  
  for (i=1; i<=N2;i++){
  
    scal=i*pi/L;
    res=res + (A_m[i-1]*exp(scal*y)+B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
    
  return res;
  
  }
 
 
//T'(x,y) pour noyaux séparables
double T_prime_KS(double x, double y, double alpha){

  double res;
  double scal;
  int i;
  double A_m[N2];
  double B_m[N2];
  
  gen_tab_KS(A_m, B_m,alpha);
  
  res= B0_KS(alpha);
  
  for (i=1;i<=N2;i++){
  
    scal=i*pi/L;
    res=res +(scal*A_m[i-1]*exp(scal*y) -scal*B_m[i-1]*exp(-scal*y))*cos(scal*x);
    }
  return res;
  }
