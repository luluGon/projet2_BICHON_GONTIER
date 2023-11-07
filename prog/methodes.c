#include "methodes.h"

double adomain(double (*f)(double), double x,double alpha){
  
  //On crée nos variables
  double somme;
  double inte;
  int i;
  double scal;
  
  //On initialise pour la première itération
  somme=f(x)/alpha;
  
  //On crée notre premiere fonction u0
  double u0(double x1){
    return f(x1);
    }
  
  //On fait une boucle pour définir les Un+1 et faire la somme
  for (i=1;i<=n;i++)
  {
  
    //On crée une fonction qui dépend du Un qui renverra son intégrale sur  diviser par alpha
    double fct(double x2){
    
      //On crée nos variables locales à la boucle
      double res;
      double intej;
      int j;
      double somme2=0.0;
      
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      
      //On calculel'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (j=1;j<n=;j++){
        
        scal=j*pi/L;
        intej=integrale_GL(u0, 0.0,L,m)
        
        
        somme2=somme2+intej*j*cos(scal*x2)/sinh(scal*H);
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=integrale_GL(u0,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      return res/alpha;
      }
    
    
    
    somme=somme +fct(x);
    
    //On redéfinie notre fonction Un+1
    double u0(x){
      return fct(x);
      }
      
    }
    return somme;  
    }
   
   
   
      

