#include "methodes.h"

double adomain(double (*f)(double), double x){
  
  double somme;
  double inte;
  int i;
  double scal;
  
  somme=f(x);
  double u0(double x1){
    return f(x1);
    }
  
  for (i=1;i<=n;i++)
  {
    double fct(double x2){
      double res;
      int j;
      double somme2=0.0;
      
      double inte(m){
        double somme1=0.0;
        
        somme1=integrale_GL(u0, 0.0,L,m);
       }
      
      for (j=1;j<n=;j++){
        scal=j*pi/L;
        somme2=somme2+inte(j)*j*cos(scal*x2)/sinh(scal*H);
        somme2=somme2*2.0*pi/(L*L);
        }
      res=integrale_GL(u0,0.0,L,0);
      res=res/(H*L);
      
      res= res+somme2;
      return res;
      }
    
    somme=somme +fct(x);
    
    double u0(x){
      return fct(x);
      }
      
    }
    return somme;  
    }
   
      

