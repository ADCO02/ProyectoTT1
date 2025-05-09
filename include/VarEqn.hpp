#ifndef _VAREQN_
#define _VAREQN_

#include "matrix.hpp"
#include "IERS.hpp"
#include "global.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "SAT_Const.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "AccelHarmonic.hpp"
#include "G_AccelHarmonic.hpp"

Matrix& VarEqn(double x, Matrix& yPhi);

#endif