#ifndef _GIBBS_
#define _GIBBS_

#include "matrix.hpp"
#include "angl.hpp"
#include "SAT_Const.hpp"
#include "unit.hpp"
#include <tuple>

tuple<Matrix&, double, double, double, string> gibbs(Matrix& r1, Matrix& r2, Matrix& r3);

#endif