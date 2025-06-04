#include "../include/gast.hpp"

double gast(double Mjd_UT1){
    double fmod_gstime = fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI );
    double gstime = (fmod_gstime < 0) ? fmod_gstime + 2 * M_PI : fmod_gstime;
    return gstime;
}