//Dans ce Headers, on définie les fonctions utilisées pour l'inversion de matrice dans la méthode des noyaux séparables
#ifndef fonctions
#define fonctions

    #include "constante.h"
    
    double determinant(double [][N1+1], double);

    void cofactor(double [][N1+1], double);

    void transpose(double [][N1+1], double [][N1+1], double);

    void inv_mat(double a[N1+1][N1+1]);
    


#endif
