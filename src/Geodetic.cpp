#include "../include/Geodetic.hpp"

tuple<double, double, double> Geodetic(Matrix& r){
    double R_equ = R_Earth;
    double f     = f_Earth;

    double epsRequ = eps*R_equ;        // Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = r(1);                   // Cartesian coordinates
    double Y = r(2);
    double Z = r(3);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis

    // Check validity of input data
    if (norm(r)==0.0){
        cout << "ERROR: invalid input in Geodetic constructor/n";
        double lon = 0.0;
        double lat = 0.0;
        double h   = -R_Earth;
        return tie(lon, lat, h);
    }

    // Iteration 
    double dZ = e2*Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;

    while(1){
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( abs(dZ-dZ_new) < epsRequ ){
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    double lon = atan2 ( Y, X );
    double lat = atan2 ( ZdZ, sqrt(rho2) );
    double h   = Nh - N;
    return tie(lon,lat,h);

}