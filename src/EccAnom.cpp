#include "../include/EccAnom.hpp"

double EccAnom(double M, double e){
    int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, 2*M_PI);
    if (M < 0) {
        M += 2*M_PI;
    }

    double E = 0;
    if (e<0.8){
        E = M; 
    }else{
        E = M_PI;
    }

    double f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );

    // Iteration
    while (fabs(f) > 1e2 * DBL_EPSILON){   
        f = E - e*sin(E) - M;
        E = E - f / ( 1.0 - e*cos(E) );
        i = i+1;
        if (i==maxit){
            cout << " convergence problems in EccAnom" << endl;
            exit(EXIT_FAILURE);
        }
    }

    return E;
}