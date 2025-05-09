#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include <tuple>
#include "matrix.hpp"

using namespace std;

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n);

#endif