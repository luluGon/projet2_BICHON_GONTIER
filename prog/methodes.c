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
  
  int i,m,K,l,j;
  
  
  //On initialise pour la première itération
  somme_Un=h(x)/alpha;
  
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=h((X[i]+1.0)*L/2.0)/alpha;
    
    ubis[i]=h((X[i]+1.0)*L/2.0)/alpha;
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
      res=-res/alpha;
      
      
    
    //Ensuite on le somme à la somme des Un(x) 
    somme_Un=somme_Un+res;
    
    
    //On calcule les points pour utiliser Gauss Legendre à la prochaine itération
    for (K=0;K<N;K++){
      //On va décomposer l'intégrale de Un(x)*k(x,t) en deux morceaux
      
      res=0.0;
      somme_K=0.0;
      //On calcule l'intégrale de Un(x)*(le terme de droite de k(x,t) ) de 0 à L.
      for (int m1=1;m1<=N1;m1++){
        scal=m1*pi/L;
        inte_Uncos=GL2(ubis,0.0,L,m1);
        
        somme_K=somme_K+inte_Uncos*m1*cos(scal*((1.0+X[K])*L/2.0)/sinh(scal*H));
        
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
void precalcul_KUn(int n){
  char *c="dat/Un.dat";
  //On crée un tableau pour les Un pris en les points de Gauss Legendre qui vont nous servir pour calculer les Un+1
  double u[N];
  //On crée aussi un tableau tampon qui sera égal à celui de Un, que l'on utilisera pour les calcul
  double ubis[N];
  double  inte_Uncos;
  double somme_K;
  double res;
  int i, K,m,l;
  double scal;
  
  //On initie les tableaux
  for (i=0;i<N;i++){
  
    //on oublie pas le changement de variable
    u[i]=h((X[i]+1.0)*L/2.0);
    ubis[i]=h((X[i]+1.0)*L/2.0);
    }
  
  FILE *fichier = fopen(c, "w");
  fprintf(fichier, "{%lf,%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf}\n", u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],u[9]);
  for (i=1;i<=n;i++)
  {
    res=0.0;
    
    
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
      
      //On somme le tout et on applique le signe moins
      res= -(res+somme_K);
      u[K]=res;
      }
      
    
    //on rédéfinie notre liste intermédiare ubis
    for(l=0;l<N;l++){
      ubis[l]=u[l];
      }
    fprintf(fichier, ",{%lf,%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf}\n", u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],u[9]);}
    fclose(fichier);
    }
  

//l'idée derrière cette version de adomain est qu'on peut factoriser les Un par 1/alpha puissance n
// Donc on peut créer un tableau prenant les N points de Gauss Legendre au préalable qui ne dépendra pas de alpha
//Une fois ce tableau fait, on aura une méthode de Gauss Legendre plus rapide
double adomain_v1_2(double x, double Un[][10],double alpha, int n){
  double somme_K;
  double somme_Un;
  double res;
  double inte_Uncos;
  double scal;
  
  
    somme_Un=h(x)/alpha;

    
  for (int i=0;i<n;i++){
    for (int m=1;m<=N1;m++){
        scal=m*pi/L;
        inte_Uncos=GL2(Un[i],0.0,L,m);
        
        
        somme_K=somme_K+inte_Uncos*m*cos(scal*x)/sinh(scal*H);
        
        }
      somme_K=somme_K*2.0*pi/(L*L);
      
      //On calcule l'intégrale de Un(x)*(terme de gauche de k(x,t) ) de 0 à L
      res=GL2(Un[i],0.0,L,0);
      res=res/(H*L);
      
      //On somme le tout et on divise par (alpha)** i
      //res correspondra à la valeur de Ui+1(x)
      res= (res+somme_K)/pow(alpha,(double) i);
     //Ensuite on le somme à la somme des Un(x) 
      somme_Un=somme_Un+res;}
  //On redivise par alpha car tout est factorisable par 1/alpha
  somme_Un=somme_Un/alpha;
  return somme_Un;
  }
 
      
      
      
double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n){
  
  //On définie nos variables
  
  
  //scal définie pour m*pi/L qui est très utilisé
  double scal;
  //res sera la variable permettant de calculer 
  double res;
  //inte_Uncos sera la variable pour calculer les integrales de Un(t)*cos(m*pi*t/L)
  double inte_Uncos;
  //somme_K donnera 
  double somme_K=0.0;
  
  int i,m,K,l;
  
  
  
  
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
  for (i=1;i<n;i++)
  {
    res=0.0;
    
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
  
    
