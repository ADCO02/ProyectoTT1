#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "matrix.hpp"
#include <cmath>

Matrix& TimeUpdate(Matrix P, Matrix Phi, double Qdt = 0.0);

#endif