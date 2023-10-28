#ifndef Gauss_Legendre
#define Gauss_Legendre


#include "constante.h"
#include <math.h>


/* fonction approximant l'intégrale de la fonction f sur [a;b] avec la quadrature de Gauss-Legendre*/

/*Elle prend en entrée dans l'ordre f une fonction d'un double, a la borne inferieur, b la borne superieur, et m un entier qui si m =! 0,multiplie la fonction d'entrée f par la fonction qui a x associe cos(m*pi*x/L) */

/*Cette dernière option est là car dans ce problème cela arrive souvent d'avoir un facteur cos(m*pi*x/L) dans l'integrale, mettre 0 pour avoir un Gauss Legendre classique */

double integrale_GL(double (*f)(double), double a, double b,int m);

#endif
