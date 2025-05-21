#ifndef _ANGLESG_
#define _ANGLESG_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "LTC.hpp"
#include "IERS.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "Geodetic.hpp"
#include "gibbs.hpp"
#include "hgibbs.hpp"
#include "elements.hpp"
#include "rpoly.hpp"
#include <tuple>

tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3,
    double Mjd1, double Mjd2, double Mjd3, Matrix& Rs1, Matrix& Rs2, Matrix& Rs3);

#endif