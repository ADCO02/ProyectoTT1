#ifndef _IERS_
#define _IERS_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "global.hpp"
#include <cmath>
#include <tuple>

tuple<double,double,double,double,double,double,double,double,double> IERS(double Mjd_UTC, char interp = 'n');

#endif