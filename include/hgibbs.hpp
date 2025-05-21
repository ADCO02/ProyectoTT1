#ifndef _HGIBBS_
#define _HGIBBS_

#include "matrix.hpp"
#include "angl.hpp"
#include "SAT_Const.hpp"
#include "unit.hpp"
#include <tuple>

tuple<Matrix&, double, double, double, string> hgibbs(Matrix& r1, Matrix& r2, Matrix& r3,
    double Mjd1, double Mjd2, double Mjd3);

#endif