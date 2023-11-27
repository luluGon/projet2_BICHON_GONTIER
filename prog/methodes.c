#include "methodes.h"

//ADOMAIN V1
double adomain_v1( double x,double alpha, int n){
  
  //On définie nos variables
  
  //somme_Un la variable qui représentera la somme des Un
  double somme_Un;
  //scal définie pour m*pi/L qui est très utilisé
  double scal;
  //res sera la variable permettant de calculer 
  double res;
  //inte_Uncos sera la variable pour calculer les integrales de Un(t)*cos(m*pi*t/L)
  double inte_Uncos;
  //somme_K donnera 
  double somme_K=0.0;
  
  double lambda= -1.0/alpha;
  
  int i,m,K1,l,j;
  
  
  //On initialise pour la première itération
  somme_Un=-lambda*h(x);
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=-h((X[i]+1.0)*L*0.5)*lambda;
    
    ubis[i]=-h((X[i]+1.0)*L*0.5)*lambda;
    }
  
  //On va maintenant boucler sur le n pour calculer les Un 
  for (j=1;j<=n;j++)
  {
    res=0.0;
    
    //On va venir calculer la valeur de Ui+1(x) et l'ajouter à la somme
    for (m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*x)/sinh(scal*H);
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      res=lambda*res;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme_Un=somme_Un+res;
    
    
    //On calcule les points pour utiliser Gauss Legendre à la prochaine itération
    for (K1=0;K1<N;K1++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme_K=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (int m1=1;m1<=N1;m1++){
        scal=m1*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m1);
        
        somme_K=somme_K+inte_Uncos*m1*cos(scal*((1.0+X[K1])*L*0.5)/sinh(scal*H));
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      u[K1]=lambda*res;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
    return somme_Un;  
    }
      
//ADOMAIN V2

//Calcul de K(s0,s1,..,sn) décrit dans l'énoncé, sur les points de Gauss Legendre, multiplier par le produit des poids w[s[i]]
//s={s0,...sn}, où s0 a sn sont des indices de 0 a 10, indiquant sur quel point de Gauss Legendre on se place
double Kn(double x, int s[], int n)
{   
    double prod=h(0.5*L*(X[s[n]]+1)); //on initialise à h(sn)
    int j;
      for (j=n;j>=2;j--)
      {
        prod=prod*K(0.5*L*(X[s[j-1]]+1),0.5*L*(X[s[j]]+1))*w[s[j]]; //on applique directement le changement de variable pour être de 0 à L
      }
    prod=prod*K(x,0.5*L*(X[s[1]]+1))*w[s[1]]; //on vient mettre le terme dépendant de x de côté
    
    return prod;
    
}


//On calcule les Un(x) sans le coefficient selon alpha devant
double U_n(double x,int n){
  int s[n];
  double somme=0.0;
  double c=L/2.; //correspond à la constante multiplicative du changement de variable induit par Gauss Legendre
  for (int it=0;it<=n;it++){
    s[it]=0;
    }
  while (s[0]==0){
    
    somme+=Kn(x,s,n);
    s[n] ++;
    for (int j=n;j>=1;j--)
    {
      if (s[j]==10){
        s[j]=0;
        s[j-1] ++ ;
        }
    }
  }
  return pow(c,n)*somme; 
  }
  
double adomain_v2(double x, double alpha, int n){
  double somme;
  double lambda= -1.0/alpha;
  double res;
  for (int i=1; i<=n;i++){
    somme+=pow(lambda,n)*U_n(x,i);
    }
  res=lambda*(somme+h(x)); 
  return res; 
  }


//APPROXIMATIONS SUCCESIVES
double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n){
  
  //On définie nos variables
  
  
  //scal définie pour m*pi/L qui est très utilisé
  double scal;
  //res sera la variable permettant de calculer 
  double res;
  //inte_Uncos sera la variable pour calculer les integrales de Un(t)*cos(m*pi*t/L)
  double inte_Uncos;
  //somme_K donnera la somme des intégrale du terme de droite de k fois Un
  double somme_K=0.0;
  
  int i,m,K1,l;
  
  
  
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=fonction_initial((X[i]+1.0)*L*0.5);
    ubis[i]=fonction_initial((X[i]+1.0)*L*0.5);
    }
  
  //On va maintenant boucler sur le n pour calculer les Un Un
  for (i=1;i<n;i++)
  {
    
    //On calcule les points pour utiliser Gauss Legendre à la prochaine itération
    for (K1=0;K1<N;K1++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme_K=0.0;
      
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (int m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*((1.0+X[K1])*L*0.5)/sinh(scal*H));
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      u[K1]=(h((1.0+X[K1])*L*0.5)-res)/alpha;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
   //On va venir calculer la valeur de Un(x) et l'ajouter à la somme
    for (m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*x)/sinh(scal*H);
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      res=(h(x)-res)/alpha;
    return res;  
    }
  
  
//NOYAUX SÉPARÉS
double ker_sep(double x,double alpha){   
  // Déclarer des variables pour stocker des informations supplémentaires
  lapack_int nnew = N1+1; // Dimension de la matrice A
  lapack_int nrhs = 1; // Nombre de colonnes de la matrice B
  // Déclarer les matrices et les vecteurs
  int IPIV[N1+1] = {0, 0};
  lapack_int info = 0;
  double A[nnew][nnew];
  double B[nnew];
  double In[nnew][nnew];
  double lambda=-1./alpha;
  double Abis[(nnew)*(nnew)];
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
  mult_mat_scal(lambda,A);
  //On fait intervenir un programme de soustraction de matrices qui renvoie le résultat dans A
  sous_matrice(In,A);
 
 //o transforme la matrice en un tableau pour utiliser dgesv
 for (int l=0;l<(N1+1);l++){
    for (int c=0;c<(N1+1);c++){
      Abis[(N1+1)*l +c]=A[l][c];
      }
    }
  
  //On fait intervenir notre programme d'inversion de matrice
  // Effectuer la factorisation LU et résoudre le système
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, nnew, nrhs, Abis, nnew, IPIV, B, nnew);
  //On définit la fonction u, solution du problème, à l'aide de ce que l'on a calculé précédemment
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
