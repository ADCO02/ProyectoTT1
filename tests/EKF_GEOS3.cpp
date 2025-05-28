#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Position.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\Accel.hpp"
#include "..\include\R_z.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\LTC.hpp"
#include "..\include\anglesg.hpp"
#include <iostream>
#include <tuple>
#include <string.h>

using namespace std;

int main() {
    //--------------------------------------------------------------------------
    //
    //  Initial Orbit Determination using Gauss and Extended Kalman Filter methods
    //
    // References:
    //   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
    //   Applications", Springer Verlag, Heidelberg, 2000.
    //   
    //   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
    //   4th Edition, 2013.
    //
    //   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
    //
    // Last modified:   2020/03/16   Meysam Mahooti
    //--------------------------------------------------------------------------
    
    //en matlab usar tik tok despues del clean para medir el tiempo

    DE430Coeff(2285, 1020);
    
    GGM03S(181);

    // Model parameters
    // AuxParam = struct ('Mjd_UTC',0,'n',0,'m',0);

    // read Earth orientation parameters
    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s 
    //  ----------------------------------------------------------------------------------------------------
    eop19620101(21413);

    int nobs = 46;
    Matrix& obs = zeros(nobs,4);

    // read observations
    FILE *fid = fopen("../data/GEOS3.txt","r");
    if(fid == NULL){
        printf("Fail open GEOS3.txt file\n");
        exit(EXIT_FAILURE);
    }

    char tline[256];
    char sub[10];
    double YY,M,D,h,m,ss,az,el,Dist;
    for (int i=1;i<=nobs;i++){
        fgets(tline, sizeof(tline), fid);
        if (strlen(tline) < 4){
            break;
        }
        strncpy(sub, &tline[0], 4);
        sub[4] = '\0';
        YY = atof(sub);
        strncpy(sub, &tline[5], 2);
        sub[2] = '\0';
        M = atof(sub);
        strncpy(sub, &tline[8], 2);
        sub[2] = '\0';
        D = atof(sub);
        strncpy(sub, &tline[12], 2);
        sub[2] = '\0';
        h = atof(sub);
        strncpy(sub, &tline[15], 2);
        sub[2] = '\0';
        m = atof(sub);
        strncpy(sub, &tline[18], 6);
        sub[6] = '\0';
        ss = atof(sub);
        strncpy(sub, &tline[25], 8);
        sub[8] = '\0';
        az = atof(sub);
        strncpy(sub, &tline[35], 7);
        sub[7] = '\0';
        el = atof(sub);
        strncpy(sub, &tline[44], 10);
        sub[10] = '\0';
        Dist = atof(sub);

        obs(i,1) = Mjday(YY,M,D,h,m,ss);
        obs(i,2) = Rad*az;
        obs(i,3) = Rad*el;
        obs(i,4) = 1e3*Dist;
    }

    fclose(fid);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix& Rs = Position(lon, lat, alt).transpose();

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    Matrix& r2 = zeros(1);
    Matrix& v2 = zeros(1);
    tie(r2,v2) = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),
                    Mjd1,Mjd2,Mjd3,Rs.transpose(),Rs.transpose(),Rs.transpose());
    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Matrix& Y0_apr = r2.union_vector(v2).transpose();

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    int n_eqn  = 6;

    Matrix& Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr.transpose());

    Matrix P = zeros(6,6);
    
    for (int i=1; i<=3; i++){
        P(i,i)=1e8;
    }
    for (int i=4; i<=6; i++){
        P(i,i)=1e3;
    }

    Matrix& LT = LTC(lon,lat);

    Matrix& yPhi = zeros(42);
    Matrix& Phi  = zeros(6,6);

    // Measurement loop
    double t = 0.0;

    //Declaración de variables
    double t_old;

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,
        UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_TT,Mjd_UT1,theta,Azim,Elev;

    Matrix& Y_old = zeros(1);
    Matrix& U = zeros(1);
    Matrix& r = zeros(1);
    Matrix& s = zeros(1);
    Matrix& dAds = zeros(1);
    Matrix& dEds = zeros(1);
    Matrix& dAdY = zeros(1);
    Matrix& K = zeros(1);
    Matrix& dEdY = zeros(1);
    Matrix& dDds = zeros(1);
    Matrix& dDdY = zeros(1);
    
    for (int i=1; i<=nobs; i++){
        // Previous step
        t_old = t;
        Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = obs(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
        
        tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
        tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
            
        for (int ii=1; ii<=6; ii++){
            yPhi(ii) = Y_old(ii);
            for (int j=1; j<=6; j++){
                if (ii==j){
                    yPhi(6*j+ii) = 1.0; 
                }else{
                    yPhi(6*j+ii) = 0.0;
                }
            }
        }
        
        yPhi = DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
        
        // Extract state transition matrices
        for (int j=1; j<=6; j++){
            Phi.assign_column(j,yPhi.extract_vector(6*j+1,6*j+6));
        }

        Y = DEInteg (Accel,0.0,t-t_old,1e-13,1e-6,6,Y_old);
        
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        U = R_z(theta);
        r = Y.extract_vector(1,3).transpose();

        s = LT*(U*r-Rs);                          // Topocentric position [m]
        
        // Time update
        P = TimeUpdate(P, Phi);
            
        // Azimuth and partials
        tie(Azim, Elev, dAds, dEds) = AzElPa(s.transpose());     // Azimuth, Elevation
        dAdY = (dAds*LT*U).union_vector(zeros(3));
        
        // Measurement update
        tie(K, Y, P) = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
        
        // Elevation and partials
        r = Y.extract_vector(1,3).transpose();
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        tie(Azim, Elev, dAds, dEds) = AzElPa(s.transpose());     // Azimuth, Elevation
        dEdY = (dEds*LT*U).union_vector(zeros(3));
        
        // Measurement update
        tie(K, Y, P) = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
        
        // Range and partials
        r = Y.extract_vector(1,3).transpose();
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        Dist = norm(s.transpose()); dDds = (s/Dist).transpose();         // Range
        dDdY=(dDds*LT*U).union_vector(zeros(3));
        
        // Measurement update
        tie(K, Y, P) = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
    }

    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(obs(46,1),'l');
    tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix& Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

    Matrix& Y_true = zeros(6,1);
    Y_true(1) = 5753.173e3;
    Y_true(2) = 2673.361e3;
    Y_true(3) = 3440.304e3;
    Y_true(4) = 4.324207e3;
    Y_true(5) = -1.924299e3;
    Y_true(6) = -5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n",Y0(1)-Y_true(1));
    printf("dY%10.1f [m]\n",Y0(2)-Y_true(2));
    printf("dZ%10.1f [m]\n",Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n",Y0(4)-Y_true(4));
    printf("dVy%8.1f [m/s]\n",Y0(5)-Y_true(5));
    printf("dVz%8.1f [m/s]\n",Y0(6)-Y_true(6));

    return 0;
}