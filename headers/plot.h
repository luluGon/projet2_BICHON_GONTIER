//Header des fonctions que l'on utilise pour les tracés

#ifndef plot
#define plot

  #include "err.h"
  #include "constante.h"
  #include "fonctions.h"
  #include "T_alpha.h"
  #include "newton.h"
  
  //fonction pour écrire les valeurs de l'erreur par la méthode d'Adomain
  double plot_err_alpha_adov1(double alpha_deb,double alpha_fin, double pas,int n);

  //fonction pour écrire les valeurs de l'erreur par la méthode d'Adomain donnée dans l'énoncé du DM
  double plot_err_alpha_adov2(double alpha_deb,double alpha_fin,double pas,int n);
  
  //fonction pour écrire les valeurs de l'erreur par la méthode des approximations successives
  double plot_err_alpha_approx_succesive(double (*f)(double),double alpha_deb,double alpha_fin, double pas,int n);
  
  //fonction pour écrire les valeurs de l'erreur par la méthode des noyaux séparables
  double plot_err_alpha_ker_sep(double alpha_deb,double alpha_fin,double pas);
  
  //fonction pour écrire les valeurs de f3 approchée par la méthode d'Adomain
  int plot_f3_adov1( int nx,int n,double alpha);
  
  //fonction pour écrire les valeurs de f3 approchée par la méthode d'Adomain donnée dans l'énoncé du DM
  int plot_f3_adov2( int nx,int n,double alpha);

  //fonction pour écrire les valeurs de f3 exacte
  int plot_f3ex( int nx);

  //fonction pour écrire les valeurs de f3 approchée par la méthode des approximations successives
  int plot_f3_app_succ(double (*f)(double), int nx,int n,double alpha);
  
  //fonction pour écrire les valeurs de f3 approchée par la méthode des noyaux séparables
  int plot_f3_ker_sep(int nx,double alpha);
  
  /*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
  //Pour adomain v1
  int plot_Tout_ado1(int nx,int ny,double alpha,int n);
  
  /*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
  //Pour adomain v2
  int plot_Tout_ado2(int nx,int ny,double alpha,int n);
  
  /*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
  //Pour approximation succesives
  int plot_Tout_AS(int nx,int ny,double (*f)(double),double alpha,int n);  
  
  /*fonction renvoyant le fichier .dat associé au calcul de Tex,Tcal et Tex-Tcal avec comme entrée le nombre de points de discretisation selon l'axe x et celui de y*/
  //Pour noyaux séparables
  int plot_Tout_KS(int nx,int ny,double alpha);
  
  //fonction renvoyant le fichier .dat contenant les x,y de la fontiere approchée Gamma_app, les y sont calculé avec Newton
  int plot_frontiere_app_KS(int nx ,double alpha, char *c);
  
  //fonction renvoyant le fichier .dat contenant les x,y de la fontiere exacte
  int plot_frontiere_ex(int nx);
  
#endif

