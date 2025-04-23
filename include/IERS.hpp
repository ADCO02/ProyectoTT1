#ifndef _IERS_
#define _IERS_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <tuple>

tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix &eop, double Mjd_UTC, char interp = 'n');

#endif