#include "plot.h"

//Algorithme pour tracé l'erreur en utilisant la méthode d'Adomain

double plot_err_alpha_adov1(double alpha_deb,double alpha_fin,double pas,int n){

  double alpha;
  double res;
  double min;
  double alpha_min;
  
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_adov1.dat";
  
  //ouverture du fichier err_alpha_adov1.dat dans le dossier dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return -1.0;
      }
  else {
  
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos données initiales, celle-ci est équivalente à celle sur T
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

//Algorithme pour tracé l'erreur en utilisant la méthode d'Adomain par la méthode proposée dans l'énoncé du DM

double plot_err_alpha_adov2(double alpha_deb,double alpha_fin,double pas,int n){

  double alpha;
  double res;
  double min;
  double alpha_min;
  
  
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_adov2.dat";
  
  //ouverture du fichier err_alpha.dat dans le dossier dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return -1.0;
      }
  else {
  
 
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos donnée initiales, celle ci est equivalente à celle sur T
    res=L2_err_f3_adov2(alpha,n);
    
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
  
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_approx_succesive.dat";
  
  //ouverture du fichier err_alpha_approx_succesive.dat dans le dossier dat 
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return -1.0;
      }
  else {
  
  
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos données initiales, celle-ci est équivalente à celle sur T
    res=L2_err_f3_approx_succesive(f,alpha,n);
    
    //On crée une condition pour récupérer le alpha qui a l'erreur minimal
    if (res<min){
      min=res;
      alpha_min=alpha;
      }
      
    //On écrit en colonne dans le fichier le alpha et son erreur correspondante
    fprintf(fichier, "%lf %lf\n", alpha, res);
    
    //on incrémente alpha du pas
    alpha=alpha+pas;
    
    }
  //fermeture du fichier
  fclose(fichier);
  }
  //On renvoie le alpha tel que l'erreur est la plus petite
  return alpha_min;
}

//Algorithme pour tracé l'erreur en utilisant la méthode des noyaux séparables
  
double plot_err_alpha_ker_sep(double alpha_deb,double alpha_fin,double pas){

  double alpha;
  double res;
  double min;
  double alpha_min;
  
  alpha=alpha_deb;
  
  min=100.0;
  alpha_min=0.0;
  
  char *c="dat/err_alpha_ker_sep.dat";
  
  //ouverture du fichier err_alpha_ker_sep.dat dans le dossier dat 
  FILE *fichier2 = fopen(c, "w");
  
  if (fichier2 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return -1.0;
      }
  else {
  
  while (alpha<=alpha_fin){
  
    //on fait l'erreur sur les f3 car, avec nos donnée initiales, celle ci est equivalente à celle sur T
    res=L2_err_ker_sep_f3(alpha);
    
    //On crée une condition pour récupérer le alpha qui a l'erreur minimal
    if (res<min){
      min=res;
      alpha_min=alpha;
      }
      
    //On écrit en colonne dans le fichier le alpha et son erreur correspondante
    fprintf(fichier2, "%lf %lf\n", alpha, res);
    
    //on incremente alpha du pas
    alpha=alpha+pas;
    
    }
  //fermuture du fichier
  fclose(fichier2);
  }
  //On renvoie le alpha tel que l'erreur est la plus petite
  return alpha_min;
  }

//Algorithme pour tracé l'approximation de f3^alpha en utilisant la méthode d'Adomain

int plot_f3_adov1( int nx,int n,double alpha){
  
  double x;
  double f3_adov1;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_adov1.dat";
  
  //ouverture d'un fichier du nom f3_adov1.dat dans le repertoire /dat"
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

//Algorithme pour tracé l'approximation de f3^alpha en utilisant la méthode d'Adomain décrite dans le sujet du DM

int plot_f3_adov2( int nx,int n,double alpha){
  
  double x;
  double f3_adov2;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_adov2.dat";
  
  //ouverture d'un fichier du nom T(x,H).dat dans le repertoire dat"
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
  else {
  for (i = 0; i < nx; i++) {
        //discretisation selon x
        x=i*pas;
        f3_adov2=adomain_v2(x, alpha,n);
      
        //ecriture des x T(x,H) en colonne dans le fichier
        fprintf(fichier, "%lf %lf\n", x,f3_adov2);
        
        }
  //fermeture du fichier
  fclose(fichier);
  
  }
  return 0;
    
  }

//Algorithme pour tracé f3
int plot_f3ex( int nx){
  
  double x;
  double f3ex;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3ex.dat";
  
  //ouverture d'un fichier du nom f3ex.dat dans le repertoire dat"
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

//Algorithme pour tracé l'approximation de f3^alpha en utilisant la méthode des approximations successives
int plot_f3_app_succ(double (*f)(double), int nx,int n,double alpha){
  
  double x;
  double f3_app_succ;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_app_succ.dat";
  
  //ouverture d'un fichier du nom f3_app_succ.dat dans le repertoire dat"
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
  
//Algorithme pour tracé l'approximation de f3^alpha en utilisant la méthode des noyaux séparables

int plot_f3_ker_sep(int nx,double alpha){
  
  double x;
  double f3_ker_sep;
  //on caclule le pas pour la discretisation selon x
  double pas=L/((double) (nx-1));
  int i;
  
  char *c="dat/f3_ker_sep.dat";
  
  //ouverture d'un fichier du nom f3_ker_sep.dat dans le repertoire dat"
  FILE *fichier = fopen(c, "w");
  
  if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
  else {
  for (i = 0; i < nx; i++) {
        //discretisation selon x
        x=i*pas;
        f3_ker_sep=ker_sep(x, alpha);
      
        //ecriture des x T(x,H) en colonne dans le fichier
        fprintf(fichier, "%lf %lf\n", x,f3_ker_sep);
        
        }
  //fermeture du fichier
  fclose(fichier);
  
  }
  return 0;
    
}


/*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
//Pour adomain v1
int plot_Tout_ado1(int nx,int ny,double alpha,int n){
  int i;
  int j;
  double pas_x;
  double pas_y;
  double x;
  double y;
  double Tex;
  double Tcal;
  double Tex_Tcal;
  
  //On crée le pas pour la discretisation en x,y
  pas_x=L/((double) (nx-1));
  pas_y=H/((double) (ny-1));
  
  
  //On va travailler dans 3 fichier en même temps pour eviter de faire un programme trop long
  char *c1="dat/Tex(x,y).dat";
  char *c2="dat/Tcal(x,y)_ado1.dat";
  char *c3="dat/Tex_Tcal_ado1.dat";
  
  //On ouvre Tex(x,y).dat, Tcal(x,y) et Tex_Tcal.dat dans le répertoire Projet1_Gontier/dat 
  FILE *fichier1 = fopen(c1, "w");
  FILE *fichier2 = fopen(c2, "w");
  FILE *fichier3 = fopen(c3, "w");
  
  //si il y'a un problème à l'ouverture d'un des fichiers, on renvoie une erreur
  if (fichier1 == NULL || fichier2 == NULL || fichier3 == NULL ) {
        fprintf(stderr, "Erreur lors de l'ouverture d'un fichier.\n");
        return 1;
      }
      
  else {
  //boucle sur x
  for (i=0;i<nx;i++){
    //discretisation de x
    x=i*pas_x;
    
    //boucle sur y
    for(j=0;j<ny;j++){
      //discretisation de y
      y=j*pas_y;
      
      //on calcule le Tex correspondant et on le met dans son .dat
      Tex=Uex(x,y);
      fprintf(fichier1, "%lf %lf %lf\n",x,y,Tex);
      
      //on calcule le Tcal correspondant et on le met dans son .dat
      Tcal=Talpha_ado1(x,y,alpha,n);
      fprintf(fichier2, "%lf %lf %lf\n",x,y,Tcal);
      
      //on calcule le Tex-Tcal correspondant et on le met dans son .dat
      Tex_Tcal=Tex-Tcal;
      fprintf(fichier3, "%lf %lf %lf\n",x,y,Tex_Tcal);
      }
    //on fait un ligne vide entre chaque itération de x pour que notre gnuplot puisse le déchiffrer pour faire une heatmap
    fprintf(fichier1, "\n");
    fprintf(fichier2, "\n");
    fprintf(fichier3, "\n");
    }
  } 
  //on ferme les 3 fichiers
  fclose(fichier1);
  fclose(fichier2);
  fclose(fichier3);
  return 0;
  }
  
/*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
//Pour adomain v2
int plot_Tout_ado2(int nx,int ny,double alpha,int n){
  int i;
  int j;
  double pas_x;
  double pas_y;
  double x;
  double y;
  double Tex;
  double Tcal;
  double Tex_Tcal;
  
  //On crée le pas pour la discretisation en x,y
  pas_x=L/((double) (nx-1));
  pas_y=H/((double) (ny-1));
  
  
  //On va travailler dans 3 fichier en même temps pour eviter de faire un programme trop long
  char *c1="dat/Tex(x,y).dat";
  char *c2="dat/Tcal(x,y)_ado2.dat";
  char *c3="dat/Tex_Tcal_ado2.dat";
  
  //On ouvre Tex(x,y).dat, Tcal(x,y) et Tex_Tcal.dat dans le répertoire Projet1_Gontier/dat 
  FILE *fichier1 = fopen(c1, "w");
  FILE *fichier2 = fopen(c2, "w");
  FILE *fichier3 = fopen(c3, "w");
  
  //si il y'a un problème à l'ouverture d'un des fichiers, on renvoie une erreur
  if (fichier1 == NULL || fichier2 == NULL || fichier3 == NULL ) {
        fprintf(stderr, "Erreur lors de l'ouverture d'un fichier.\n");
        return 1;
      }
      
  else {
  //boucle sur x
  for (i=0;i<nx;i++){
    //discretisation de x
    x=i*pas_x;
    
    //boucle sur y
    for(j=0;j<ny;j++){
      //discretisation de y
      y=j*pas_y;
      
      //on calcule le Tex correspondant et on le met dans son .dat
      Tex=Uex(x,y);
      fprintf(fichier1, "%lf %lf %lf\n",x,y,Tex);
      
      //on calcule le Tcal correspondant et on le met dans son .dat
      Tcal=Talpha_ado2(x,y,alpha,n);
      fprintf(fichier2, "%lf %lf %lf\n",x,y,Tcal);
      
      //on calcule le Tex-Tcal correspondant et on le met dans son .dat
      Tex_Tcal=Tex-Tcal;
      fprintf(fichier3, "%lf %lf %lf\n",x,y,Tex_Tcal);
      }
    //on fait un ligne vide entre chaque itération de x pour que notre gnuplot puisse le déchiffrer pour faire une heatmap
    fprintf(fichier1, "\n");
    fprintf(fichier2, "\n");
    fprintf(fichier3, "\n");
    }
  } 
  //on ferme les 3 fichiers
  fclose(fichier1);
  fclose(fichier2);
  fclose(fichier3);
  return 0;
  }
  
/*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
//Pour approximation succesives
int plot_Tout_AS(int nx,int ny,double (*f)(double),double alpha,int n){
  int i;
  int j;
  double pas_x;
  double pas_y;
  double x;
  double y;
  double Tex;
  double Tcal;
  double Tex_Tcal;
  
  //On crée le pas pour la discretisation en x,y
  pas_x=L/((double) (nx-1));
  pas_y=H/((double) (ny-1));
  
  
  //On va travailler dans 3 fichier en même temps pour eviter de faire un programme trop long
  char *c1="dat/Tex(x,y).dat";
  char *c2="dat/Tcal(x,y)_AS.dat";
  char *c3="dat/Tex_Tcal_AS.dat";
  
  //On ouvre Tex(x,y).dat, Tcal(x,y) et Tex_Tcal.dat dans le répertoire Projet1_Gontier/dat 
  FILE *fichier1 = fopen(c1, "w");
  FILE *fichier2 = fopen(c2, "w");
  FILE *fichier3 = fopen(c3, "w");
  
  //si il y'a un problème à l'ouverture d'un des fichiers, on renvoie une erreur
  if (fichier1 == NULL || fichier2 == NULL || fichier3 == NULL ) {
        fprintf(stderr, "Erreur lors de l'ouverture d'un fichier.\n");
        return 1;
      }
      
  else {
  //boucle sur x
  for (i=0;i<nx;i++){
    //discretisation de x
    x=i*pas_x;
    
    //boucle sur y
    for(j=0;j<ny;j++){
      //discretisation de y
      y=j*pas_y;
      
      //on calcule le Tex correspondant et on le met dans son .dat
      Tex=Uex(x,y);
      fprintf(fichier1, "%lf %lf %lf\n",x,y,Tex);
      
      //on calcule le Tcal correspondant et on le met dans son .dat
      Tcal=Talpha_AS(x,y,f,alpha,n);
      fprintf(fichier2, "%lf %lf %lf\n",x,y,Tcal);
      
      //on calcule le Tex-Tcal correspondant et on le met dans son .dat
      Tex_Tcal=Tex-Tcal;
      fprintf(fichier3, "%lf %lf %lf\n",x,y,Tex_Tcal);
      }
    //on fait un ligne vide entre chaque itération de x pour que notre gnuplot puisse le déchiffrer pour faire une heatmap
    fprintf(fichier1, "\n");
    fprintf(fichier2, "\n");
    fprintf(fichier3, "\n");
    }
  } 
  //on ferme les 3 fichiers
  fclose(fichier1);
  fclose(fichier2);
  fclose(fichier3);
  return 0;
  }
  
/*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
//Pour noyaux séparables
int plot_Tout_KS(int nx,int ny,double alpha){
  int i;
  int j;
  double pas_x;
  double pas_y;
  double x;
  double y;
  double Tex;
  double Tcal;
  double Tex_Tcal;
  
  //On crée le pas pour la discretisation en x,y
  pas_x=L/((double) (nx-1));
  pas_y=H/((double) (ny-1));
  
  
  //On va travailler dans 3 fichier en même temps pour eviter de faire un programme trop long
  char *c1="dat/Tex(x,y).dat";
  char *c2="dat/Tcal(x,y)_KS.dat";
  char *c3="dat/Tex_Tcal_KS.dat";
  
  //On ouvre Tex(x,y).dat, Tcal(x,y) et Tex_Tcal.dat dans le répertoire Projet1_Gontier/dat 
  FILE *fichier1 = fopen(c1, "w");
  FILE *fichier2 = fopen(c2, "w");
  FILE *fichier3 = fopen(c3, "w");
  
  //si il y'a un problème à l'ouverture d'un des fichiers, on renvoie une erreur
  if (fichier1 == NULL || fichier2 == NULL || fichier3 == NULL ) {
        fprintf(stderr, "Erreur lors de l'ouverture d'un fichier.\n");
        return 1;
      }
      
  else {
  //boucle sur x
  for (i=0;i<nx;i++){
    //discretisation de x
    x=i*pas_x;
    
    //boucle sur y
    for(j=0;j<ny;j++){
      //discretisation de y
      y=j*pas_y;
      
      //on calcule le Tex correspondant et on le met dans son .dat
      Tex=Uex(x,y);
      fprintf(fichier1, "%lf %lf %lf\n",x,y,Tex);
      
      //on calcule le Tcal correspondant et on le met dans son .dat
      Tcal=Talpha_KS(x,y,alpha);
      fprintf(fichier2, "%lf %lf %lf\n",x,y,Tcal);
      
      //on calcule le Tex-Tcal correspondant et on le met dans son .dat
      Tex_Tcal=Tex-Tcal;
      fprintf(fichier3, "%lf %lf %lf\n",x,y,Tex_Tcal);
      }
    //on fait un ligne vide entre chaque itération de x pour que notre gnuplot puisse le déchiffrer pour faire une heatmap
    fprintf(fichier1, "\n");
    fprintf(fichier2, "\n");
    fprintf(fichier3, "\n");
    }
  } 
  //on ferme les 3 fichiers
  fclose(fichier1);
  fclose(fichier2);
  fclose(fichier3);
  return 0;
  }
  
  
//fonction renvoyant le fichier .dat contenant les x,y de la fontiere approchée Gamma_app, les y sont calculé avec Newton
int plot_frontiere_app_KS(int nx ,double alpha, char *c){
    
    int i;
    //On crée le pas pour discrétiser les x
    double pas=L/((double) (nx-1));
    double x;
    double y;
    
  
    //ouverture des fichier du nom Gamma_cal.dat  dans le repertoire /Projet1_Gontier/dat"
    FILE *fichier = fopen(c, "w");
    
    
    if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
    else {
  
    for (i=0;i<nx;i++){
    
      x=i*pas;
      y=newton_T_alpha_KS(x, 0.2, 1e-6,  100 , alpha);
      
      
      //ecriture des x y qui appartiennent à Gamma_app en colonne dans le fichier
      fprintf(fichier, "%lf %lf\n", x,y);
      
      }
    }
    
    //On ferme le fichier 
    fclose(fichier);
    return 0;
    }
    
//fonction renvoyant le fichier .dat contenant les x,y de la fontiere exacte
int plot_frontiere_ex(int nx){
    
    int i;
    //On crée le pas pour discrétiser les x
    double pas=L/((double) (nx-1));
    double x;
    double y;
    
    char *c1="dat/Gamma_ex.dat";
    
  
    //ouverture des fichier du nom Gamma_cal.dat  dans le repertoire /Projet1_Gontier/dat"
    FILE *fichier = fopen(c1, "w");
    
    
    if (fichier == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier.\n");
        return 1;
      }
      
    else {
  
    for (i=0;i<nx;i++){
    
      x=i*pas;
      y=0.7;
      
      //ecriture des x y qui appartiennent à Gamma_app en colonne dans le fichier
      fprintf(fichier, "%lf %lf\n", x,y);
      
      }
    }
    
    //On ferme le fichier 
    fclose(fichier);
    return 0;
    }
