/**
 * @file NutMatrix.cpp
 * @brief Implementation of the NutMatrix function.
 */

#include "../include/NutMatrix.hpp"

Matrix& NutMatrix(double Mjd_TT){
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    auto [dpsi, deps] = NutAngles (Mjd_TT);

    // Transformation from mean to true equator and equinox
    Matrix& NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

    return NutMat;
}