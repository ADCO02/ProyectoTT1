/**
 * @file Cheb3D.cpp
 * @brief Implementation of the Cheb3D function.
 */

#include "../include/Cheb3D.hpp"

Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz){
    // Check validity
    if ( (t<Ta) || (Tb<t) ){
        cout << "ERROR: Time out of range in Cheb3D::Value\n";
		exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);  

    Matrix& f1 = zeros(3);
    Matrix& f2 = zeros(3);
    Matrix& old_f1 = zeros(3);

    Matrix& vaux = zeros(3);
    for (int i = N; i >= 2; i--) {
        old_f1 = f1;
        vaux(1)=Cx(i); vaux(2)=Cy(i); vaux(3)=Cz(i);
        f1 = ((f1*(2*tau))-f2)+vaux;
        f2 = old_f1;
    }


    vaux(1)=Cx(1); vaux(2)=Cy(1); vaux(3)=Cz(1);
    Matrix& ChebApp = (((f1*tau)-f2)+vaux);
    return ChebApp;
}