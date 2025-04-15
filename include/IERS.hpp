#ifndef _IERS_
#define _IERS_

#include "matrix.hpp"
#include "SAT_Const.hpp"

void IERS(Matrix &eop, double Mjd_UTC, char interp, double* x_pole, double* y_pole, 
    double* UT1_UTC, double* LOD, double* dpsi, double* deps, double* dx_pole, 
    double* dy_pole, double* TAI_UTC );

inline void IERS(Matrix &eop, double Mjd_UTC,
    double* x_pole, double* y_pole,
    double* UT1_UTC, double* LOD, double* dpsi,
    double* deps, double* dx_pole, double* dy_pole, double* TAI_UTC)
{
IERS(eop, Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi,
   deps, dx_pole, dy_pole, TAI_UTC);
}

#endif