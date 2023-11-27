//Dans ce Headers, on définie les fonctions qui nous seront utiles aux projets
#ifndef fonctions
#define fonctions
    
    #include "Gauss_Legendre.h"
    #include "constante.h"
    #include <stdio.h>
    
    //fonction prenant deux double x et y en entrée et renvoie la solution exacte en x,y
    double Uex(double x, double y);

    //fonction prenant en entrée un double x et renvoyant f0(x)
    double f0(double x);

    //fonction prenant en entrée un double x et renvoyant q0(x)
    double q0(double x);

    //fonction prenant en entrée un double x et un entier i et renvoyant alpha_i(x), la fonction implicite de k(x,t)
    double alpha_i(double x,int i);
    
    //fonction prenant en entrée un double t et un entier i et renvoyant beta_i(x), la fonction implicite de k(x,t)
    double beta_i(double t,int i);
    
    //fonction pour k(x,t)
    double K(double x, double t);

    //fonction prenant en entrée un double x et renvoie h(x)
    double h(double x);


#endif
