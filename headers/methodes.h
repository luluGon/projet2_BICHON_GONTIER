//Dans ce headers, on définie les fontions implémentant les différentes méthodes 
//pour calculer les intgérales de Fredholm

#ifndef methodes
#define methodes

    #include "fonctions.h"
    #include "Gauss_Legendre.h"
    #include <math.h>


    //Cette fonction prend comme argument d'entrée, f une fonction de double, et x un double
    //Elle renvoie la valeur de f_3(x) approchée par la méthode de adomain
    double adomain(double (*f)(double), double x);

#endif

