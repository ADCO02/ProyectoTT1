#include "..\include\anglesg.hpp"

tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3,
    double Mjd1, double Mjd2, double Mjd3, Matrix& Rs1, Matrix& Rs2, Matrix& Rs3){

    Rs1 = Rs1.transpose();
    Rs2 = Rs2.transpose();
    Rs3 = Rs3.transpose();

    Matrix& L1 = zeros(3,1);
    L1(1) = cos(el1)*sin(az1); L1(2) = cos(el1)*cos(az1); L1(3) = sin(el1);
    Matrix& L2 = zeros(3,1);
    L2(1) = cos(el2)*sin(az2); L2(2) = cos(el2)*cos(az2); L2(3) = sin(el2);
    Matrix& L3 = zeros(3,1);
    L3(1) = cos(el3)*sin(az3); L3(2) = cos(el3)*cos(az3); L3(3) = sin(el3);

    double lon1, lat1, h1;
    tie(lon1, lat1, h1) = Geodetic(Rs1.transpose());
    double lon2, lat2, h2;
    tie(lon2, lat2, h2) = Geodetic(Rs2.transpose());
    double lon3, lat3, h3;
    tie(lon3, lat3, h3) = Geodetic(Rs3.transpose());

    Matrix& M1 = LTC(lon1, lat1);
    Matrix& M2 = LTC(lon2, lat2);
    Matrix& M3 = LTC(lon3, lat3);

    // body-fixed system
    Matrix& Lb1 = M1.transpose()*L1;
    Matrix& Lb2 = M1.transpose()*L2;
    Matrix& Lb3 = M1.transpose()*L3;

    // mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix& P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix& N = NutMatrix(Mjd_TT);
    Matrix& T = N * P;
    Matrix& E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
    Matrix& Lm1 = E.transpose()*Lb1;
    Rs1 = E.transpose()*(Rs1);

    Mjd_UTC = Mjd2;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
    tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix& Lm2 = E.transpose()*Lb2;
    Rs2 = E.transpose()*Rs2;

    Mjd_UTC = Mjd3;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
    tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix& Lm3 = E.transpose()*Lb3;
    Rs3 = E.transpose()*Rs3;

    // geocentric inertial position
    double tau1 = (Mjd1-Mjd2)*86400;
    double tau3 = (Mjd3-Mjd2)*86400;

    double a1 = tau3/(tau3-tau1);
    double a3 =-tau1/(tau3-tau1);

    double b1 = tau3/(6*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau3,2));
    double b3 =-tau1/(6*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau1,2));

    Lm1 = Lm1.transpose();
    Lm2 = Lm2.transpose();
    Lm3 = Lm3.transpose();

    Rs1 = Rs1.transpose();
    Rs2 = Rs2.transpose();
    Rs3 = Rs3.transpose();

    Matrix& aux_D1 = zeros(3,3);
    aux_D1.assign_column(1,Lm1);
    aux_D1.assign_column(2,Lm2);
    aux_D1.assign_column(3,Lm3);
    Matrix& aux_D2 = zeros(3,3);
    aux_D2.assign_column(1,Rs1);
    aux_D2.assign_column(2,Rs2);
    aux_D2.assign_column(3,Rs2);
    Matrix& D = inv(aux_D1)*aux_D2;

    double d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    double d2s = D(2,1)*b1+D(2,3)*b3;

    double Ccye = 2*dot(Lm2,Rs2);

    double poly[10], zeror[10], zeroi[10];

    poly[0]=  1.0;  // R2^8... polynomial
    poly[1]=  0.0;
    poly[2]=  -( pow(d1s,2) + d1s*Ccye + pow((norm(Rs2)),2));
    poly[3]=  0.0;
    poly[4]=  0.0;
    poly[5]=  -GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    poly[6]=  0.0;
    poly[7]=  0.0;
    poly[8]=  -pow(GM_Earth,2)*pow(d2s,2);
    int nr = real_poly_roots( poly, 8, zeror, zeroi );

    double bigr2= -99999990.0;

    for (int j=1;j<=nr;j++){
        if (( zeror[j-1] > bigr2 ) && ( fabs(zeroi[j-1]) < 10e-12) ){
            bigr2= zeror[j-1];
        }  
    }

    double u = GM_Earth/(pow(bigr2,3));

    double C1 = a1+b1*u;
    double C2 = -1;
    double C3 = a3+b3*u;

    Matrix& aux_temp = zeros(3);
    aux_temp(1) = C1; aux_temp(2) = C2; aux_temp(3) = C3;
    Matrix& temp = D*(aux_temp.transpose())*(-1);
    double rho1 = temp(1)/(a1+b1*u);
    double rho2 = -temp(2);
    double rho3 = temp(3)/(a3+b3*u);

    double rhoold1 = rho1;
    double rhoold2 = rho2;
    double rhoold3 = rho3;

    rho2 = 99999999.9;
    int ll   = 0;

    Matrix& v2 = zeros(3);
    while ((abs(rhoold2-rho2) > 1e-12) && (ll <= 0 )){
        ll = ll + 1;
        rho2 = rhoold2;
        
        Matrix& r1 = Rs1+(Lm1*rho1);
        Matrix& r2 = Rs2+(Lm2*rho2);
        Matrix& r3 = Rs3+(Lm3*rho3);
        
        double magr1 = norm(r1);
        double magr2 = norm(r2);
        double magr3 = norm(r3);
        double theta, theta1, copa;
        string error;
        tie(v2, theta,theta1,copa,error) = gibbs(r1,r2,r3);        
        if ( (error != "          ok") & (copa < M_PI/180) ){      
            tie(v2,theta,theta1,copa,error) = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);
        }
        
        double p, a, e, i, Omega, omega, M;
        tie(p, a, e, i, Omega, omega, M) = elements (r2.union_vector(v2));

        double f1, g1, f3, g3;
        
        if ( ll <= 8 ){
            u = GM_Earth/pow(magr2,3);
            double rdot= dot(r2,v2)/magr2;
            double udot= (-3*GM_Earth*rdot)/(pow(magr2,4));            
            double tausqr= tau1*tau1;
            f1=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau1
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau1;
            g1= tau1 - (1/6)*u*tau1*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau1 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau3*tau3;
            f3=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau3 
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau3;
            g3= tau3 - (1/6)*u*tau3*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau3 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
        }else{
            
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );
            
            f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
            f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
        }
        
        C1 = g3/(f1*g3-f3*g1);
        C2 = -1;
        C3 =-g1/(f1*g3-f3*g1);
        
        double H1 = GM_Earth*tau3/12;
        double H3 =-GM_Earth*tau1/12;
        double H2 = H1-H3;
        
        double G1 = -tau3/(tau1*(tau3-tau1));
        double G3 = -tau1/(tau3*(tau3-tau1));
        double G2 = G1-G3;
        
        double D1 = G1+H1/pow(magr1,3);
        double D2 = G2+H2/pow(magr2,3);
        double D3 = G3+H3/pow(magr3,3);
        
        Matrix& aux_temp1 = zeros(3);
        aux_temp1(1) = D1; aux_temp1(2) = D2; aux_temp1(3) = D3;
        Matrix& aux_temp2 = zeros(3);
        aux_temp2(1) = C1; aux_temp2(2) = C2; aux_temp2(3) = C3;
        temp = aux_temp1*aux_temp2.transpose()*(-1);
        rhoold1 = temp(1)/(a1+b1*u);
        rhoold2 = -temp(1);
        rhoold3 = temp(1)/(a3+b3*u);
        
        r1 = Rs1+Lm1*rhoold1;
        r2 = Rs2+Lm2*rhoold2;
        r3 = Rs3+Lm3*rhoold3;
        
    }

    Matrix& r1 = Rs1+Lm1*rho1;
    Matrix& r2 = Rs2+Lm2*rho2;
    Matrix& r3 = Rs3+Lm3*rho3;

    return tie(r2,v2);
}