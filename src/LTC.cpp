/**
 * @file LTC.cpp
 * @brief Implementation of the LTC function.
 */

#include "../include/LTC.hpp"

Matrix& LTC(double lon, double lat){
    Matrix& M = (R_y(lat*(-1.0)))*(R_z(lon));

    for (int j=1; j<=3; j++){
        double Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;
    }

    return M;
}