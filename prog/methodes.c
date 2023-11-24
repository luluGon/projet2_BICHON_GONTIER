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
  somme=f(x)/alpha;
  
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
        
        
        somme2=somme2+intej*m1*cos(scal*x)/sinh(scal*H);
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


double ker_sep(double (*h)(double), double x,double alpha){   
  double A[N1+1][N1+1];
  double B[N1+1];
  //double c[N1+1];
  double In[N1+1][N1+1];
  double lambda=-1./alpha;
  int i,m;
  //On construit f
  double f(double x){
    return (1./alpha)*h(x);     
  }
  //On initialise A et B
  for (i=1;i<=N1+1;i++){
    for (m=1;m<=N1+1;m++){
      //On définit la fonction à placer dans l'intégrale pour calculer A
      double f_A(double x){           
        return beta_i(x,m)*alpha_i(x,i);
      }
      A[m-1][i-1]=integrale_GL(f_A,0.,L,0);      
        //On choisit m=0 car cela implique une multiplication de f_A par cos(0)=1
      //On initialise ensuite I
      if (m==i){
        In[i-1][i-1]=1.;
      }
      else{
        In[m-1][i-1]=0.;
      }
    }
    //On définit la fonction à placer dans l'intégrale pour calculer B
    double f_B(double x){
      return beta_i(x,i)*f(x);
    }
    B[i-1]=integrale_GL(f_B,0.,L,0);   
    //On choisit m=0 car cela implique une multiplication de f_B par cos(0)=1
  }
  //On fait intervenir un programme de multiplication de matrice par un scalaire qui renvoie le résultat dans A
  //mult_mat_scal(lambda,A);
  //On fait intervenir un programme de soustraction de matrices qui renvoie le résultat dans A
  //sous_matrice(In,A);
  //On fait intervenir notre programme d'inversion de matrice
  // Déclarer les matrices et les vecteurs
  int IPIV[N1+1] = {0, 0};
  lapack_int info = 0;

  // Déclarer des variables pour stocker des informations supplémentaires
  lapack_int n = N1+1; // Dimension de la matrice A
  lapack_int nrhs = N1+1; // Nombre de colonnes de la matrice B

  // Effectuer la factorisation LU et résoudre le système
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, IPIV, B, n);
  //On définit la fonction u, solution du problème, à l'aide de ce que l'on a calculé précédemment
  //printf("c= %lf\n",c);
  double u(double x){
    double som=0.;
    int m;
    for (m=1;m<=N1+1;m++){
      som+=B[m-1]*alpha_i(x,m);
    }    
    som=f(x)+lambda*som;
    return som;
  }
  return u(x);
}
