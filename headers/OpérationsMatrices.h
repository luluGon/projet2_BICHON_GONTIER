//Dans ce Headers, on définie les fonctions utilisées pour les opérations simple sur des matrices de taille N1+1xN1+1, dans la méthode des noyaux séparables
#ifndef fonctions
#define fonctions

    #include "constante.h"
    
    
    void prod_mat(double A[N1+1][N1+1],double B[N1+1]);
    void mult_mat_scal(double lambda,double A[N1+1][N1+1]);
    void sous_matrice(double I[N1+1][N1+1],double A[N1+1][N1+1]);
    

#endif
