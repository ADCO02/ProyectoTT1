#ifndef _AZELPA_
#define _AZELPA_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <cmath>
#include<tuple>

/**
 * @brief Computes azimuth, elevation and their partial derivatives from local tangent coordinates.
 * 
 * @param s Topocentric local tangent coordinates vector (East-North-Zenith frame).
 * 
 * @return A tuple containing:
 *         - Azimuth angle [rad]
 *         - Elevation angle [rad]
 *         - Partial derivatives of azimuth w.r.t s (Matrix&)
 *         - Partial derivatives of elevation w.r.t s (Matrix&)
 */
tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s);

#endif