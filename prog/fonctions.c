#include "fonctions.h"

double Uex(double x, double y){

  double scal;
  double res;
  
  scal= k*pi/L;
  
  switch (num)
    {
    case 0: 
      res=cosh(scal*y)*cos(scal*x);
      break;
    case 1:
      res=cosh(scal*y)*cos(scal*x)*x/k;
      break;
    default: printf(" Pas d'autres fonctions\n");
    }
    return res;
    }
    
  
double f0(double x){
  
    double res;
    switch (num)
    {
    case 0: 
      res=cos(k*pi*x/L);
      break;
    case 1:
      res=cos(k*pi*x/L)*x/k;
      break;
    default: printf(" Pas d'autres fonctions\n");
    }
    return res;
    }
    
double q0(double x){

  double res;
  switch (num)
  {
  case 0:
    /*q0(x)=-pi*sinh(pi*0)*cos(pi*x)*/
    /*Donc q0(x)=0 pour tout x */
    res=0;
    break;
  case 1:
    res=0;
    break;
    
  default : printf(" Pas d'autres fonctions");
  }
  return res;
  }


double kbis(double x, double t){
  double somme;
  double scal;
  double res;
  int i;
  
  somme=0.0;
  res=1.0/(H*L);
  
  for (i=1; i<=N1; i++)
  {
    scal=i*pi/L;
    
    i*cos(scal*x)*cos(scal*t)/sinh(scal*H);
  }
  somme=2.0*pi*somme/(L*L);
  
  res=res +somme;
  
  return res;
  }

double h(double x){
  /*On décompose le calcul de h en sous éléments, int1 la première intégrale, int2 la seconde, somme la somme*/
  /*De plus on pose res pour le résultat global, et scal pour les scalaires récurrent des boucles */
  
 
  double res;
  double int1;
  double int2;
  double somme;
  double scal;
  int i;
  
  somme=0.0;
  res=0.0;
    
  int1=integrale_GL(f0, 0.0, L,0);
  int1=int1/(H*L);
  
  for (i=1;i<=N1;i++){
    
    scal =pi*i/L;
    int2=integrale_GL(f0,0.0,L,i);
    
    int2=(f0(x))*cosh(scal*H)*int2;
    
    somme=somme +int2*i*cos(scal*x)/sinh(scal*H);
  }
  
  somme=somme*2.0*pi/(L*L);
  
  res = -q0(x) + int1 +somme;
  
  return res;
}

