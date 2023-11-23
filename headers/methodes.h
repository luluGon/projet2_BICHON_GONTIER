//Dans ce headers, on définie les fontions implémentant les différentes méthodes 
//pour calculer les intgérales de Fredholm

#ifndef methodes
#define methodes

    #include "constante.h"
    #include "fonctions.h"
    #include "Gauss_Legendre.h"
    #include <math.h>
    
    extern double X[10];
    extern double w[10];

    //Cette fonction prend comme argument d'entrée, f une fonction de double, et x un double
    //Elle renvoie la valeur de f_3(x) approchée par la méthode de adomain
    /*double adomain(double (*f)(double), double x);*/
    
    double adomain_v1( double x,double alpha, int n);
    
    void precalcul_KUn(int n);
    double adomain_v1_2(double x, double Un[][10],double alpha, int n);
    double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n);
    
    
#endif

