#include "plot.h"


double plot_err_alpha_adov1(double alpha_deb,double alpha_fin,double pas,int n){

  double alpha;
  double res;
  double min;
  double alpha_min;
  
  
  //on commence à alpha=1e-15 pour rester en double
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_adov1.dat";
  
  //ouverture du fichier err_alpha.dat dans le dossier /dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1.0;
      }
  else {
  
  //on va jusque alpha =10 par pas de 0.0001, cela nous permet une première approximation
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos donnée initiales, celle ci est equivalente à celle sur T
    res=L2_err_f3_adov1(alpha,n);
    
    //On crée une condition pour récupérer le alpha qui a l'erreur minimal
    if (res<min){
      min=res;
      alpha_min=alpha;
      }
      
    //On écrit en colonne dans le fichier le alpha et son erreur correspondante
    fprintf(fichier, "%lf %lf\n", alpha, res);
    
    //on incremente alpha du pas
    alpha=alpha+pas;
    
    }
  //fermuture du fichier
  fclose(fichier);
  }
  //On renvoie le alpha tel que l'erreur est la plus petite
  return alpha_min;
  }

//Meme chose mais avec methode d'approximation sucessive
double plot_err_alpha_approx_succesive(double (*f)(double),double alpha_deb,double alpha_fin,double pas,int n){

  double alpha;
  double res;
  double min;
  double alpha_min;
  
  
  //on commence à alpha=1e-15 pour rester en double
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_approx_succesive.dat";
  
  //ouverture du fichier err_alpha.dat dans le dossier /dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1.0;
      }
  else {
  
  //on va jusque alpha =10 par pas de 0.0001, cela nous permet une première approximation
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos donnée initiales, celle ci est equivalente à celle sur T
    res=L2_err_f3_approx_succesive(f,alpha,n);
    
    //On crée une condition pour récupérer le alpha qui a l'erreur minimal
    if (res<min){
      min=res;
      alpha_min=alpha;
      }
      
    //On écrit en colonne dans le fichier le alpha et son erreur correspondante
    fprintf(fichier, "%lf %lf\n", alpha, res);
    
    //on incremente alpha du pas
    alpha=alpha+pas;
    
    }
  //fermuture du fichier
  fclose(fichier);
  }
  //On renvoie le alpha tel que l'erreur est la plus petite
  return alpha_min;
  }
  
  
int plot_f3_adov1( int nx,int n,double alpha){
  
  double x;
  double f3_adov1;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_adov1.dat";
  
  //ouverture d'un fichier du nom T(x,H).dat dans le repertoire /Projet1_Gontier/dat"
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
  else {
  for (i = 0; i < nx; i++) {
        //discretisation selon x
        x=i*pas;
        f3_adov1=adomain_v1(x, alpha,n);
      
        //ecriture des x T(x,H) en colonne dans le fichier
        fprintf(fichier, "%lf %lf\n", x,f3_adov1);
        
        }
  //fermeture du fichier
  fclose(fichier);
  
  }
  return 0;
    
  }
  
int plot_f3ex( int nx){
  
  double x;
  double f3ex;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3ex.dat";
  
  //ouverture d'un fichier du nom T(x,H).dat dans le repertoire /Projet1_Gontier/dat"
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
  else {
  for (i = 0; i < nx; i++) {
        //discretisation selon x
        x=i*pas;
        f3ex=Uex(x,H);
      
        //ecriture des x T(x,H) en colonne dans le fichier
        fprintf(fichier, "%lf %lf\n", x,f3ex);
        
        }
  //fermeture du fichier
  fclose(fichier);
  
  }
  return 0;
    
  }

int plot_f3_app_succ(double (*f)(double), int nx,int n,double alpha){
  
  double x;
  double f3_app_succ;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_app_succ.dat";
  
  //ouverture d'un fichier du nom T(x,H).dat dans le repertoire /Projet1_Gontier/dat"
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
  else {
  for (i = 0; i < nx; i++) {
        //discretisation selon x
        x=i*pas;
        f3_app_succ=approximation_succesive(f,x, alpha,n);
      
        //ecriture des x T(x,H) en colonne dans le fichier
        fprintf(fichier, "%lf %lf\n", x,f3_app_succ);
        
        }
  //fermeture du fichier
  fclose(fichier);
  
  }
  return 0;
    
  }
