//Dans ce Headers, on définie les fonctions utilisées pour l'inversion de matrice dans la méthode des noyaux séparables
#ifndef OperationsMatrices
#define OperationsMatrices

    #include "constante.h"
    
    //Fonction permettant d'effectuer le produit d'une matrice par un scalaire
    void mult_mat_scal(double lambda,double A[N1+1][N1+1]);
    
    //Fonction permettant d'effectuer la soustraction de deux matrices
    void sous_matrice(double In[N1+1][N1+1],double A[N1+1][N1+1]);


#endif
