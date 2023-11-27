/*Header servant a définir les constantes */

#ifndef constante
#define constante

  #include <lapacke.h>
  #include <stdio.h>
  #define pi 3.141592653589793 
  #define H 1.0 
  #define L 1.0
  #define N 10 //Le nombre de poids et points pour Gauss_Legendre //
  #define N1 1 //N1 est le nombre d'itérations associée à la somme dans k(x,t) et dans h(x)
               //On note que lorsque N1>1, les termes des sommes sont égaux à 0, donc on choisit N1=1
  
  #define N2 1 //nombre d'itération pour le calcul de T(x,y)
  #define k 1  //k est un paramètre associé à la fonction Uex
  #define num 0 //num est un paramètre pour choisir la fonction Uex à approcher
  

  
  
#endif
