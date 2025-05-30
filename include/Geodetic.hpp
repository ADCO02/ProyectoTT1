/**
 * @file Geodetic.hpp
 * 
 * @brief Converts Cartesian coordinates to geodetic coordinates.
 * 
 * Computes geodetic longitude, latitude, and altitude from a position vector.
 */

#ifndef _GEODETIC_
#define _GEODETIC_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <tuple>

/**
 * @brief Computes geodetic coordinates from position vector.
 * 
 * @param r Position vector in Cartesian coordinates [m].
 * 
 * @return tuple containing:
 *   - longitude [rad]
 *   - latitude [rad]
 *   - altitude [m]
 */
tuple<double, double, double> Geodetic(Matrix& r);

#endif