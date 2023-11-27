#include "newton.h"

double newton_T_alpha_ado1(double x,double y_init, double err, int ite_max,double alpha,int n) {

    int i = 0;
    double y=y_init;
    
    //on s'impose un nombre d'itération maximum pour eviter les boucle infinies 
    while (i < ite_max) {
        double fy = Talpha_ado1(x,y,alpha,n)-Uex(x,x);
        double fpy = T_prime_ado1(x,y,alpha,n);
        
        if (fabs(fy) < err) {
            
            return y;
        }
        
        y = y - fy / fpy;
        i++;
    }
    
    return y;  
}


double newton_T_alpha_ado2(double x,double y_init, double err, int ite_max,double alpha,int n) {

    int i = 0;
    double y=y_init;
    
    //on s'impose un nombre d'itération maximum pour eviter les boucle infinies 
    while (i < ite_max) {
        double fy = Talpha_ado2(x,y,alpha,n)-Uex(x,x);
        double fpy = T_prime_ado2(x,y,alpha,n);
        
        if (fabs(fy) < err) {
            
            return y;
        }
        
        y = y - fy / fpy;
        i++;
    }
    
    return y;  
}

double newton_T_alpha_AS(double x,double y_init, double err, int ite_max,double (*f)(double),double alpha,int n) {

    int i = 0;
    double y=y_init;
    
    //on s'impose un nombre d'itération maximum pour eviter les boucle infinies 
    while (i < ite_max) {
        double fy = Talpha_AS(x,y,f,alpha,n)-Uex(x,x);
        double fpy = T_prime_AS(x,y,f,alpha,n);
        
        if (fabs(fy) < err) {
            
            return y;
        }
        
        y = y - fy / fpy;
        i++;
    }
    
    return y;  
}



double newton_T_alpha_KS(double x,double y_init, double err, int ite_max,double alpha) {

    int i = 0;
    double y=y_init;
    
    //on s'impose un nombre d'itération maximum pour eviter les boucle infinies 
    while (i < ite_max) {
        double fy = Talpha_KS(x,y,alpha)-Uex(x,0.7);
        double fpy = T_prime_KS(x,y,alpha);
        
        if (fabs(fy) < err) {
            
            return y;
        }
        
        y = y - fy / fpy;
        i++;
    }
    
    return y;  
}
