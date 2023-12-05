#include <stdio.h>
#include "methodes.h"
#include "constante.h"
#include "Gauss_Legendre.h"
#include "fonctions.h"
#include "err.h"
#include "plot.h"
#include "OperationsMatrices.h"
#include "coefficients.h"
#include "T_alpha.h"



int nx=100; //nombres d'intervalle sur l'axe des abcisses pour le tracé de f3 et de T(x,y)
int ny=100; //nombres d'intervalle sur l'axe des ordonnées pour le tracé de f3 et de T(x,y) 


//DECLARATION DE VARIABLES POUR F3 ADOMAIN V1
double alpha_opt_adov1;
double alpha_debut_adov1=0.001;
double alpha_fin_adov1=10.;
double pas_alpha_adov1=0.001;
int n_adov1=5;
double ado1_err_min;

//DECLARATION DE VARIABLES POUR F3 ADOMAIN V2
double alpha_opt_adov2;
double alpha_debut_adov2=0.01;
double alpha_fin_adov2=10.;
double pas_alpha_adov2=0.01;
int n_adov2=3;
double ado2_err_min;
  
//DECLARATION DE VARIABLES POUR F3 APPROXIMATIONS SUCCESIVES
double alpha_opt_AS;
double alpha_debut_AS=0.001;
double alpha_fin_AS=10.;
double pas_alpha_AS=0.001;
int n_AS=5;
double AS_err_min;
double f(double x){return 0;}

//DECLARATION DE VARIABLES POUR F3 NOYAUX SEPARÉS
double alpha_opt_ker_sep;
double alpha_debut_ker_sep=0.000001;
double alpha_fin_ker_sep=0.00001;
double pas_alpha_ker_sep=0.00000001;
double ker_sep_err_min;

int main(){
  
  //TRACÉ DE LA SOLUTION EXACTE
  plot_f3ex(nx);
  
  
  /*
  //COMPARAISON TEMPS ADOMAIN V1 ET ADOMAIN V2
  // on commente la ligne que l'on veut tester, et ensuite on compile puis fait time ./exe
  //adomain_v1(0.2, 1., 7);
  //adomain_v2(0.2, 1., 7);
  */
  
  
  //ALPHA OPTIMAL POUR ADOMAIN V1 ET TRACÉ DE F3
    
  //on calcule le alpha optimal et son erreur associé
  alpha_opt_adov1=plot_err_alpha_adov1(alpha_debut_adov1,alpha_fin_adov1,pas_alpha_adov1,n_adov1);
  ado1_err_min = L2_err_f3_adov1(alpha_opt_adov1 , n_adov1);
  
  //on en profite pour tracer la courbe de f3 pour adomain v1
  plot_f3_adov1( nx,n_adov1, alpha_opt_adov1);
  
  printf("Dans la première version de la méthode de Adomain.\n");
  printf("Pour n = %d, avec un pas selon alpha de : %lf\n",n_adov1,pas_alpha_adov1);
  printf("L'erreur en norme L2 de f3-Uex(.,H) est minimal pour :\n");
  printf("alpha = %lf \net vaut : %lf\n", alpha_opt_adov1,ado1_err_min );
  
  
  
  //ALPHA OPTIMAL POUR ADOMAIN V2 ET TRACÉ DE F3
  
  
  //on calcule le alpha optimal et son erreur associé
  alpha_opt_adov2=plot_err_alpha_adov2(alpha_debut_adov2,alpha_fin_adov2,pas_alpha_adov2,n_adov2);
  ado2_err_min = L2_err_f3_adov2(alpha_opt_adov2 , n_adov2);
  
  
  //on en profite pour tracer la courbe de f3 pour approximations succesives
  plot_f3_adov2(nx,n_adov2,alpha_opt_adov2);
  
  printf("Dans la seconde version de la méthode de Adomain.\n");
  printf("Pour n = %d, avec un pas selon alpha de : %lf\n",n_adov2,pas_alpha_adov2);
  printf("L'erreur en norme L2 de f3-Uex(.,H) est minimal pour :\n");
  printf("alpha = %lf \net vaut : %lf\n", alpha_opt_adov2,ado2_err_min );
  
  
  
  
  //ALPHA OPTIMAL POUR APPROXIMATIONS SUCCESIVES ET TRACÉ DE F3
  
  
  //on calcule le alpha optimal et son erreur associé
  alpha_opt_AS=plot_err_alpha_approx_succesive(f,alpha_debut_AS,alpha_fin_AS,pas_alpha_AS,n_AS);
  AS_err_min = L2_err_f3_approx_succesive(f,alpha_opt_AS , n_AS);
  
  
  //on en profite pour tracer la courbe de f3 pour approximations succesives
  plot_f3_app_succ(f, nx,n_AS,alpha_opt_AS);
  
  
  printf("Dans la méthode d'approximation succesive.\n");
  printf("Avec comme fonction de départ, la fonction nulle\n");//à modifier si changement de fonction de départ
  printf("Pour n = %d, avec un pas selon alpha de : %lf\n",n_AS,pas_alpha_AS);
  printf("L'erreur en norme L2 de f3-Uex(.,H) est minimal pour :\n");
  printf("alpha = %lf \net vaut : %lf\n", alpha_opt_AS,AS_err_min );
  
  
  
  //ALPHA OPTIMAL POUR NOYAUX SÉPARÉS ET TRACÉ DE F3
  
  //on calcule le alpha optimal et son erreur associé
  alpha_opt_ker_sep=plot_err_alpha_ker_sep( alpha_debut_ker_sep , alpha_fin_ker_sep , pas_alpha_ker_sep );
  ker_sep_err_min=L2_err_ker_sep_f3(alpha_opt_ker_sep);
  
  //On en profite pour tracer la courbe de f3 pour les noyaux séparés
  plot_f3_ker_sep(nx,alpha_opt_ker_sep);
  
  printf("Dans la méthode des noyaux séparés.\n");
  printf("Pour un pas selon alpha de : %.8lf\n",pas_alpha_ker_sep);
  printf("L'erreur en norme L2 de f3-Uex(.,H) est minimal pour :\n");
  printf("alpha = %.8lf \net vaut : %.8lf\n", alpha_opt_ker_sep,ker_sep_err_min);
  
  
  //
  //TRACÉ DES T(X,Y)
  //
  
  
  //TRACÉ DE T(X,Y) PAR LA MÉTHODE D'ADOMAIN 1, AVEC SON ALPHA OPTIMAL PRÉCALCULÉ
  alpha_opt_adov1=7.449000;
  plot_Tout_ado1(nx,ny,alpha_opt_adov1,n_adov1);
  
  
  
  //TRACÉ DE T(X,Y) PAR LA MÉTHODE D'ADOMAIN 2, AVEC SON ALPHA OPTIMAL PRÉCALCULÉ
  alpha_opt_adov2=0.27;
  plot_Tout_ado2(nx,ny,alpha_opt_adov2,n_adov2);
  
  
  
  //TRACÉ DE T(X,Y) PAR LA MÉTHODE D'APPROXIMATION SUCCESIVES, AVEC SON ALPHA OPTIMAL PRÉCALCULÉ
  alpha_opt_AS=0.516000;
  plot_Tout_AS(nx,ny,f,alpha_opt_AS,n_AS);
  
  
  
  //TRACÉ DE T(X,Y) PAR LA MÉTHODE DE NOYAUX SÉPARABLES, AVEC SON ALPHA OPTIMAL PRÉCALCULÉ
  alpha_opt_ker_sep=0.00000104;
  plot_Tout_KS(nx,ny,alpha_opt_ker_sep);
  
  
  
  //
  //TRACER DE LA FRONTIERE APPROCHÉE
  //
  alpha_opt_ker_sep=0.00000104;
  plot_frontiere_app_KS(nx ,alpha_opt_ker_sep, "dat/frontiere_alpha_opt.dat");
  
  
  
  //TRACER DE LA FRONTIERE EXACTE
  plot_frontiere_ex(nx);
  
  
  return 0;
}

