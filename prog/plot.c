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
  
  char *c="dat/err_alpha.dat";
  
  //ouverture du fichier err_alpha.dat dans le dossier Projet1_Gontier/dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1.0;
      }
  else {
  
  //on va jusque alpha =10 par pas de 0.0001, cela nous permet une première approximation
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos donnée initiales, celle ci est equivalente à celle sur T
    res=L2_err_f3(alpha,n);
    
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
 
