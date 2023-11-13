#include "Gauss_Legendre.h"


double w[10]={
0.2955242247147529,
0.2955242247147529,
0.2692667193099963,
0.2692667193099963,
0.2190863625159820,
0.2190863625159820,
0.1494513491505806,
0.1494513491505806,
0.0666713443086881,
0.0666713443086881};

double X[10]={
-0.1488743389816312,
0.1488743389816312,
-0.4333953941292472,
0.4333953941292472,
-0.6794095682990244,
0.6794095682990244,
-0.8650633666889845,
0.8650633666889845,
-0.9739065285171717,
0.9739065285171717};

double integrale_GL( double (*f)(double), double a, double b,int m){
  
  /*On definie les constante interne au programme, b_moins_a pour (b-a)/2 et b_plus_a pour (b+a)/2 qui sont issue de notre changement de variable pour appliquer Gauss_Legendre sur tout intervalle [a;b] et pas sur [-1;1]*/
  
  double b_moins_a;
  double b_plus_a;
  double scal;
  double Xbis[N];
  
  int i;
  double res;

  res=0.0;
  scal=m*pi/L;
  b_moins_a=0.5*(b-a);
  b_plus_a=0.5*(b+a);
  
  /*Comme indiqué dans le header, on fait deux cas, un pour Gauss_Legendre classique, m=0*/
  if (m==0){
    for (i=0;i<N;i++){
    
      Xbis[i]=b_moins_a*X[i] + b_plus_a;
      res=res + w[i]*((*f)(Xbis[i]));
      }
      }
      
  /* Un autre pour multiplier la fonction f(x) par cos(m*pi*x/L) pour tout m =!0 */
  
  /*Note: On aurait pu ne pas faire de distinction de cas car cos(0)=1, mais j'ai préféré faire une distinction de cas par soucis d'optimisation */
  
  else {
    for (i=0;i<N;i++){
    
      Xbis[i]=b_moins_a*X[i] + b_plus_a;
      res=res + w[i]*((*f)(Xbis[i]))*cos(scal*Xbis[i]);
      }
      } 
  res=res*(b_moins_a);
  
  return res;
}
  
  
double GL2( double tab[N], double a, double b,int m){
  
  /*On definie les constante interne au programme, b_moins_a pour (b-a)/2 et b_plus_a pour (b+a)/2 qui sont issue de notre changement de variable pour appliquer Gauss_Legendre sur tout intervalle [a;b] et pas sur [-1;1]*/
  
  double b_moins_a;
  double b_plus_a;
  double scal;
  double Xbis[N];
  int i;
  double res;

  res=0.0;
  scal=m*pi/L;
  b_moins_a=0.5*(b-a);
  b_plus_a=0.5*(b+a);
  
  if (m==0){
  for (i=0;i<N;i++){
    
      Xbis[i]=b_moins_a*X[i] + b_plus_a;
      
      res=res + w[i]*(tab[i]);
      }
      }
      
  /* Un autre pour multiplier la fonction f(x) par cos(m*pi*x/L) pour tout m =!0 */
  
  /*Note: On aurait pu ne pas faire de distinction de cas car cos(0)=1, mais j'ai préféré faire une distinction de cas par soucis d'optimisation */
  
  else {
    for (i=0;i<N;i++){
    
      Xbis[i]=b_moins_a*X[i] + b_plus_a;
      res=res + w[i]*(tab[i])*cos(scal*X[i]);
      }
      } 
  res=res*(b_moins_a);
  
  return res;
}
