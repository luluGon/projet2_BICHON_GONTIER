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

    //fonction intermédiaire qui prend entrée deux double x et t et renvoie k(x,t)
    double kbis(double x, double t);

    //fonction prenant en entrée un double x et renvoie h(x)
    double h(double x);


#endif
