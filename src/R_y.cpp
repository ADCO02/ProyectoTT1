/**
 * @file R_y.cpp
 * @brief Implementation of the R_y function.
 */

#include "../include/R_x.hpp"

Matrix& R_y(double angle){
    double C = cos(angle);
    double S = sin(angle);
    Matrix& rotmat = zeros(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}