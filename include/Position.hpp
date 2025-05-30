/**
 * @file Position.hpp
 * 
 * @brief Computes the position vector from geodetic coordinates.
 *
 * Converts geodetic coordinates (longitude, latitude, altitude) to a 
 * position vector in an Earth-centered Earth-fixed (ECEF) Cartesian frame.
 */

#ifndef _POSITION_
#define _POSITION_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <cmath>

/**
 * @brief Computes the position vector from geodetic coordinates.
 *
 * @param lon Geodetic longitude in radians.
 * @param lat Geodetic latitude in radians.
 * @param h   Altitude above the reference ellipsoid in meters.
 * 
 * @return Reference to a Matrix representing the 3D position vector in meters.
 */
Matrix& Position(double lon, double lat, double h);

#endif