#include "..\include\IERS.hpp"

tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix &eop, double Mjd_UTC, char interp){
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    if (interp =='l'){
        // linear interpolation
        double mjd = (floor(Mjd_UTC));
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }
        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
		    exit(EXIT_FAILURE);
        }
        Matrix preeop = eop.extract_column(i);
        Matrix nexteop = eop.extract_column(i+1);
        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
        y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
        UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
        LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
        dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
        deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
        dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
        dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
        TAI_UTC = preeop(13);
        
        x_pole  = x_pole/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = y_pole/Const::Arcs;  // Pole coordinate [rad]
        dpsi    = dpsi/Const::Arcs;
        deps    = deps/Const::Arcs;
        dx_pole = dx_pole/Const::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole/Const::Arcs; // Pole coordinate [rad]
    }else if (interp =='n'){    
        double mjd = (floor(Mjd_UTC));
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }
        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
		    exit(EXIT_FAILURE);
        }
        eop = eop.extract_column(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = eop(5)/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = eop(6)/Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = eop(7);             // UT1-UTC time difference [s]
        LOD     = eop(8);             // Length of day [s]
        dpsi    = eop(9)/Const::Arcs;
        deps    = eop(10)/Const::Arcs;
        dx_pole = eop(11)/Const::Arcs; // Pole coordinate [rad]
        dy_pole = eop(12)/Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = eop(13);            // TAI-UTC time difference [s]
    }

    return make_tuple(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
}