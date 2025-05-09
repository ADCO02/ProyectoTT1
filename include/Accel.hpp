#ifndef _ACCEL_
#define _ACCEL_

#include "matrix.hpp"
#include "global.hpp"
#include "IERS.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "SAT_Const.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "Mjday_TDB.hpp"
#include "JPL_Eph_DE430.hpp"
#include "AccelHarmonic.hpp"
#include "AccelPointMass.hpp"

Matrix& Accel(double x, Matrix& Y);

#endif