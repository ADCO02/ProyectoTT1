#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "matrix.hpp"
#include <cmath>

Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt);
Matrix TimeUpdate(Matrix P, Matrix Phi);

#endif