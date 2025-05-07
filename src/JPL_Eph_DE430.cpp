#include "..\include\JPL_Eph_DE430.hpp"

tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB){
    double JD = Mjd_TDB + 2400000.5;
    double i = -1;
    bool found = false;
    int i_aux = 1;
    while(!found && i_aux<=PC.n_row){
        if(PC(i_aux,1)<=JD && JD<=PC(i_aux,2)){
            i=i_aux;
            found=true;
        }
        i_aux++;
    }
    Matrix& PCtemp = PC.extract_row(i);

    double t1 = PCtemp(1)-2400000.5; // MJD at start of interval

    double dt = Mjd_TDB - t1;

    int temp_aux = (270-231)/13 + 1;
    Matrix& temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*13 + 231;
    }
    Matrix& Cx_Earth = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Earth = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Earth = PCtemp.extract_vector(temp(3),temp(4)-1);
    temp = temp+39;
    Matrix& Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz = PCtemp.extract_vector(temp(3),temp(4)-1);
    Cx_Earth = Cx_Earth.union_vector(Cx);
    Cy_Earth = Cy_Earth.union_vector(Cy);
    Cz_Earth = Cz_Earth.union_vector(Cz);  
    int j;
    double Mjd0;  
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix& r_Earth = (Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth.extract_vector(13*j+1,13*j+13),
                        Cy_Earth.extract_vector(13*j+1,13*j+13), Cz_Earth.extract_vector(13*j+1,13*j+13))
                        .transpose())*1e3;

    temp_aux = (480-441)/13 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*13 + 441;
    }
    Matrix& Cx_Moon = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Moon = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Moon = PCtemp.extract_vector(temp(3),temp(4)-1);
    for (int i=1; i<=7; i++){
        temp = temp+39;
        Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3),temp(4)-1);   
        Cx_Moon = Cx_Moon.union_vector(Cx);
        Cy_Moon = Cy_Moon.union_vector(Cy);
        Cz_Moon = Cz_Moon.union_vector(Cz);    
    }
    if (0<=dt && dt<=4){
        j=0;
        Mjd0 = t1;
    }else if(4<dt && dt<=8){
        j=1;
        Mjd0 = t1+4*j;
    }else if(8<dt && dt<=12){
        j=2;
        Mjd0 = t1+4*j;
    }else if(12<dt && dt<=16){
        j=3;
        Mjd0 = t1+4*j;
    }else if(16<dt && dt<=20){
        j=4;
        Mjd0 = t1+4*j;
    }else if(20<dt && dt<=24){
        j=5;
        Mjd0 = t1+4*j;
    }else if(24<dt && dt<=28){
        j=6;
        Mjd0 = t1+4*j;
    }else if(28<dt && dt<=32){
        j=7;
        Mjd0 = t1+4*j;
    }
    Matrix& r_Moon = (Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon.extract_vector(13*j+1,13*j+13),
                        Cy_Moon.extract_vector(13*j+1,13*j+13), Cz_Moon.extract_vector(13*j+1,13*j+13))
                        .transpose())*1e3;

    temp_aux = (786-753)/11 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*11 + 753;
    }
    Matrix& Cx_Sun = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Sun = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Sun = PCtemp.extract_vector(temp(3),temp(4)-1);
    temp = temp+33;
    Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
    Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
    Cz = PCtemp.extract_vector(temp(3),temp(4)-1);   
    Cx_Sun = Cx_Sun.union_vector(Cx);
    Cy_Sun = Cy_Sun.union_vector(Cy);
    Cz_Sun = Cz_Sun.union_vector(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix& r_Sun = (Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun.extract_vector(11*j+1,11*j+11),
                    Cy_Sun.extract_vector(11*j+1,11*j+11), Cz_Sun.extract_vector(11*j+1,11*j+11))
                    .transpose())*1e3;

    temp_aux = (45-3)/14 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*14 + 3;
    }
    Matrix& Cx_Mercury = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Mercury = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Mercury = PCtemp.extract_vector(temp(3),temp(4)-1);
    for (int i=1; i<=3; i++){
        temp = temp+42;
        Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3),temp(4)-1);
        Cx_Mercury = Cx_Mercury.union_vector(Cx);
        Cy_Mercury = Cy_Mercury.union_vector(Cy);
        Cz_Mercury = Cz_Mercury.union_vector(Cz);    
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix& r_Mercury = (Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury.extract_vector(14*j+1,14*j+14),
                        Cy_Mercury.extract_vector(14*j+1,14*j+14), Cz_Mercury.extract_vector(14*j+1,14*j+14))
                        .transpose())*1e3;

    temp_aux = (201-171)/10 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*10 + 171;
    }
    Matrix& Cx_Venus = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Venus = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Venus = PCtemp.extract_vector(temp(3),temp(4)-1);
    temp = temp+30;
    Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
    Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
    Cz = PCtemp.extract_vector(temp(3),temp(4)-1);
    Cx_Venus = Cx_Venus.union_vector(Cx);
    Cy_Venus = Cy_Venus.union_vector(Cy);
    Cz_Venus = Cz_Venus.union_vector(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix& r_Venus = (Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus.extract_vector(10*j+1,10*j+10),
                        Cy_Venus.extract_vector(10*j+1,10*j+10), Cz_Venus.extract_vector(10*j+1,10*j+10))
                        .transpose())*1e3;

    temp_aux = (342-309)/11 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*11 + 309;
    }
    Matrix& Cx_Mars = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Mars = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Mars = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Mars = (Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars.extract_vector(11*j+1,11*j+11),
                        Cy_Mars.extract_vector(11*j+1,11*j+11), Cz_Mars.extract_vector(11*j+1,11*j+11))
                        .transpose())*1e3;

    temp_aux = (366-342)/8 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*8 + 342;
    }
    Matrix& Cx_Jupiter = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Jupiter = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Jupiter = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Jupiter = (Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter.extract_vector(8*j+1,8*j+8),
                        Cy_Jupiter.extract_vector(8*j+1,8*j+8), Cz_Jupiter.extract_vector(8*j+1,8*j+8))
                        .transpose())*1e3;

    temp_aux = (387-366)/7 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*7 + 366;
    }
    Matrix& Cx_Saturn = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Saturn = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Saturn = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Saturn = (Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn.extract_vector(7*j+1,7*j+7),
                        Cy_Saturn.extract_vector(7*j+1,7*j+7), Cz_Saturn.extract_vector(7*j+1,7*j+7))
                        .transpose())*1e3;

    temp_aux = (405-387)/6 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*6 + 387;
    }
    Matrix& Cx_Uranus = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Uranus = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Uranus = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Uranus = (Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus.extract_vector(6*j+1,6*j+6),
                        Cy_Uranus.extract_vector(6*j+1,6*j+6), Cz_Uranus.extract_vector(6*j+1,6*j+6))
                        .transpose())*1e3;

    temp_aux = (423-405)/6 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*6 + 405;
    }
    Matrix& Cx_Neptune = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Neptune = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Neptune = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Neptune = (Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune.extract_vector(6*j+1,6*j+6),
                        Cy_Neptune.extract_vector(6*j+1,6*j+6), Cz_Neptune.extract_vector(6*j+1,6*j+6))
                        .transpose())*1e3;

    temp_aux = (441-423)/6 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*6 + 423;
    }
    Matrix& Cx_Pluto = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Pluto = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Pluto = PCtemp.extract_vector(temp(3),temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix& r_Pluto = (Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto.extract_vector(6*j+1,6*j+6),
                        Cy_Pluto.extract_vector(6*j+1,6*j+6), Cz_Pluto.extract_vector(6*j+1,6*j+6))
                        .transpose())*1e3;

    temp_aux = (839-819)/10 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*10 + 819;
    }
    Matrix& Cx_Nutations = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Nutations = PCtemp.extract_vector(temp(2),temp(3)-1);
    for (int i=1; i<=3; i++){
        temp = temp+20;
        Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
        Cx_Nutations = Cx_Nutations.union_vector(Cx);
        Cy_Nutations = Cy_Nutations.union_vector(Cy);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix& Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations.extract_vector(10*j+1,10*j+10),
                    Cy_Nutations.extract_vector(10*j+1,10*j+10),zeros(10,1)).transpose();

    temp_aux = (929-899)/10 + 1;
    temp = zeros(temp_aux);
    for (int i = 1; i <= temp_aux; i ++) {
        temp(i)=(i-1)*10 + 899;
    }
    Matrix& Cx_Librations = PCtemp.extract_vector(temp(1),temp(2)-1);
    Matrix& Cy_Librations = PCtemp.extract_vector(temp(2),temp(3)-1);
    Matrix& Cz_Librations = PCtemp.extract_vector(temp(3),temp(4)-1);
    for (int i=1; i<=3; i++){
        temp = temp+30;
        Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3),temp(4)-1);
        Cx_Librations = Cx_Librations.union_vector(Cx);
        Cy_Librations = Cy_Librations.union_vector(Cy);
        Cz_Librations = Cz_Librations.union_vector(Cz);    
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix& Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations.extract_vector(10*j+1,10*j+10),
                        Cy_Librations.extract_vector(10*j+1,10*j+10), Cz_Librations.extract_vector(10*j+1,10*j+10))
                        .transpose();
    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1.0/(1+EMRAT);
    r_Earth = r_Earth-(r_Moon*EMRAT1);
    r_Mercury = r_Mercury-r_Earth;
    r_Venus = r_Venus-r_Earth;
    r_Mars = r_Mars-r_Earth;
    r_Jupiter = r_Jupiter-r_Earth;
    r_Saturn = r_Saturn-r_Earth;
    r_Uranus = r_Uranus-r_Earth;
    r_Neptune = r_Neptune-r_Earth;
    r_Pluto = r_Pluto-r_Earth;
    r_Sun = r_Sun-r_Earth;

    return tie(r_Mercury.transpose(),r_Venus.transpose(),r_Earth.transpose(),r_Mars.transpose(),
        r_Jupiter.transpose(),r_Saturn.transpose(),r_Uranus.transpose(),r_Neptune.transpose(),
        r_Pluto.transpose(),r_Moon.transpose(),r_Sun.transpose());
}