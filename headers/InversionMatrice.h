//Dans ce Headers, on définie les fonctions utilisées pour l'inversion de matrice dans la méthode des noyaux séparables
#ifndef fonctions
#define fonctions

    #include constante.h
    
    float determinant(float [][N1+1], float);

    void cofactor(float [][N1+1], float);

    void transpose(float [][N1+1], float [][N1+1], float);

    double[N1+1] inv_mat(double a[N1+1,N1+1]);
    


#endif
