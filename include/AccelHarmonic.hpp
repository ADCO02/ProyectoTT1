#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "matrix.hpp"
#include "Legendre.hpp"
#include "global.hpp"
#include <cmath>

Matrix AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max);

#endif