//Dans ce headers, on définie les fontions implémentant les différentes méthodes 
//pour calculer les intgérales de Fredholm

#ifndef methodes
#define methodes

    #include "constante.h"
    #include "fonctions.h"
    #include "Gauss_Legendre.h"
    #include "OperationsMatrices.h"
    #include <math.h>
    #include <stdio.h>
    #include <lapacke.h>
    
    //Les vecteurs utilisés dans Gauss_Legendre.c
    extern double X[10];
    extern double w[10];

    //Cette fonction prend comme argument d'entrée x et alpha des doubles et un entier n
    //Elle renvoie la valeur de f_3(x) approchée par la méthode de adomain
    double adomain_v1( double x,double alpha, int n);
    
    //cette fonction prend en entrée x, in vecteur s d'entier allant de 0 à 9, et un entier n
    //fonction intermediaire pour la deuxieme version de Adomain
    //elle renvoie K(s0,s1,...,sn) définie dans le rapport aux points de Gauss Legendre definie par le vecteur s, le tout multiplie par le produit des poids w[s[i]]
    double Kn(double x, int s[], int n);
    
    //cette fonction prend en entrée un réel x et un entier n 
    //cette fontion renvoie la valeur de Un(x) sans le terme dépendant de alpha pour chaque n
    double U_n(double x,int n);
    
    //Cette fonction prend comme argument d'entrée un réel x, un réel alpha et un entier n.
    //Elle renvoie la valeur de f3(x) calculé par la seconde version de Adomain
    double adomain_v2(double x, double alpha, int n);
    
    //Cette fonction prend comme argument d'entrée, f une fonction de double, x et alpha des doubles et un entier n
    //Elle renvoie la valeur de f_3(x) approchée par la méthode des approximations successives
    double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n);
    
    //Cette fonction prend comme argument d'entrée x et alpha des doubles
    //Elle renvoie la valeur de f_3(x) approchée par la méthode des noyaux séparables
    double ker_sep(double x, double alpha);
#endif
