#include "methodes.h"


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
  
  int i,m,K,l;
  
  
  //On initialise pour la première itération
  somme_Un=h(x)/alpha;
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=h((X[i]+1.0)*L/2.0);
    ubis[i]=h((X[i]+1.0)*L/2.0);
    }
  
  //On va maintenant boucler sur le n pour calculer les Un Un
  for (i=1;i<=n;i++)
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
      res=-res/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme_Un=somme_Un+res;
    
    
    //On calcule les points pour utiliser Gauss Legendre à la prochaine itération
    for (K=0;K<N;K++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme_K=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (int m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*((1.0+X[K])*L/2.0)/sinh(scal*H));
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      u[K]=-res/alpha;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
    return somme_Un;  
    }
         
/*double adomain_v2(double x, double alpha, int n){
  double 
  for (int i; i<n;i++){
    */
    
double adomain_v3(double x, double alpha, int n){
  double somme_K;
  double Un[10];
  int compteur=0;
  double res;
  double inte_Uncos;
  
  char *c="dat/Un.dat";
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=h((X[i]+1.0)*L/2.0);
    ubis[i]=h((X[i]+1.0)*L/2.0);
    }
  
  FILE *fichier = fopen(c, "rw");
  while (fscanf(fichier, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Un[0],&Un[1],&Un[2],&Un[3],&Un[4],&Un[5],&Un[6],&Un[7],&Un[8],&Un[9]) != EOF && compteur<n){
    somme =0.0;
    compteur ++;
    for (m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(Un,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*x)/sinh(scal*H);
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(Un,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      res=-res/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme_Un=somme_Un+res;
    }
  fclose(fichier);
  if (n>compteur){
    printf("Pas encore assez d'itérations pré-enregistré,\n ici %d itération de Un enregistré.\nVoulez vous aller plus loin ? (Yes=1/Non=0)");
    scanf("%c\n",
      
      
      
double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n){
  
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
  
  int i,m,K,l;
  
  
  //On initialise pour la première itération
  somme_Un=h(x)/alpha;
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=fonction_initial((X[i]+1.0)*L/2.0);
    ubis[i]=fonction_initial((X[i]+1.0)*L/2.0);
    }
  
  //On va maintenant boucler sur le n pour calculer les Un Un
  for (i=1;i<=n;i++)
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
      res=(h(x)-res)/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme_Un=somme_Un+res;
    
    
    //On calcule les points pour utiliser Gauss Legendre à la prochaine itération
    for (K=0;K<N;K++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme_K=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (int m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*((1.0+X[K])*L/2.0)/sinh(scal*H));
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(ubis,0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par alpha
      res= res+somme_K;
      u[K]=(h((1.0+X[K])*L/2.0)-res)/alpha;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(l=0;l<N;l++){
      ubis[l]=u[l];
      }
  }
    return somme_Un;  
    }
  
    
