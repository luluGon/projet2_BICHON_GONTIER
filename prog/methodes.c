#include "methodes.h"


/*ne marche pas */
/*double adomain(double (*f)(double), double x,double alpha,int n){
  
  //On crée nos variables
  double somme;
  double inte;
  int i;
  double scal;
  double res;
  double intej;
  int m;
  double somme2=0.0;
  
  
  
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
      
      
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (m=1;m<=N1;m++){
        
        scal=m*pi/L;
        intej=integrale_GL(u0, 0.0,L,m)
        
        
        somme2=somme2+intej*m*cos(scal*x2)/sinh(scal*H);
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=integrale_GL(u0,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      return -res/alpha;
      }
    
    
    
    somme=somme +fct(x);
    
    //On redéfinie notre fonction Un+1
    double u0(x){
      return fct(x);
      }
      
    }
    return somme;  
    }
*/

double adomain_v2(double (*f)(double), double x,double alpha, int n){
  
  //On crée nos variables
  double somme;
  double inte;
  int i;
  double scal;
  double res2;
  double res;
  double intej;
  int m;
  double somme2=0.0;
  
  
  //On initialise pour la première itération
  somme=-f(x)/alpha;
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (int j=0;j<N;j++){
  
    //on oublie pas le changement de variable
    u[j]=(X[j]+1)*L/2.0;
    ubis[j]=(X[j]+1)*L/2.0;
    }
  
  //On va maintenant boucler sur le n dans Un
  for (i=1;i<=n;i++)
  {
    res2=0.0;
    for (int m1=1;m1<=N1;m1++){
        scal=m1*pi/L;
        intej=GL2(ubis,0.0,L,m1);
        
        
        somme2=somme2+intej*m*cos(scal*x)/sinh(scal*H);
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      res2=-res/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme=somme+res2;
    printf("somme=%lf\n",somme);
    
    
    for (int K=0;K<N;K++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme2=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (m=1;m<=N1;m++){
        scal=m*pi/L;
        intej=GL2(ubis,0.0,L,m);
        
        
        somme2=somme2+intej*m*cos(scal*((1+X[K])*L/2.0)/sinh(scal*H));
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      u[K]=-res/alpha;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(int l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
    return somme;  
    }
         
double adomain_v3( double x,double alpha, int n){
  
  //On crée nos variables
  double somme;
  double inte;
  int i;
  double scal;
  double res2;
  double res;
  double intej;
  int m;
  double somme2=0.0;
  
  
  //On initialise pour la première itération
  somme=-h(x)/alpha;
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (int j=0;j<N;j++){
  
    //on oublie pas le changement de variable
    u[j]=(X[j]+1)*L/2.0;
    ubis[j]=(X[j]+1)*L/2.0;
    }
  
  //On va maintenant boucler sur le n dans Un
  for (i=1;i<=n;i++)
  {
    res2=0.0;
    for (int m1=1;m1<=N1;m1++){
        scal=m1*pi/L;
        intej=GL2(ubis,0.0,L,m1);
        
        
        somme2=somme2+intej*m*cos(scal*x)/sinh(scal*H);
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      res2=-res/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme=somme+res2;
    printf("somme=%lf\n",somme);
    
    
    for (int K=0;K<N;K++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (m=1;m<=N1;m++){
        scal=m*pi/L;
        intej=GL2(ubis,0.0,L,m);
        
        
        somme2=somme2+intej*m*cos(scal*((1+X[K])*L/2.0)/sinh(scal*H));
        somme2=somme2*2.0*pi/(L*L);
        }
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme2;
      u[K]=-res/alpha;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(int l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
    return somme;  
    }
