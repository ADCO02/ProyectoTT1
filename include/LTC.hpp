/**
 * @file LTC.hpp
 * 
 * @brief Computes the rotation matrix from ECEF to local tangent coordinates.
 *
 * Defines the transformation from the Earth-centered, Earth-fixed (ECEF) 
 * reference frame aligned with the Greenwich meridian, to a local East-North-Up 
 * coordinate system centered at a given geodetic longitude and latitude. 
 * This rotation is essential for converting global positions into local horizon-based
 * observations, typically used in ground-based tracking and geolocation.
 */

#ifndef _LTC_
#define _LTC_

#include "matrix.hpp"
#include "R_x.hpp"
#include "R_y.hpp"
#include "R_z.hpp"

/**
 * @brief Computes the rotation matrix from ECEF to local tangent coordinates.
 *
 * @param lon Geodetic East longitude in radians.
 * @param lat Geodetic latitude in radians.
 * 
 * @return Reference to a 3×3 rotation matrix that transforms vectors from 
 *         the Earth equator and Greenwich meridian system to the local 
 *         East-North-Zenith (tangent) coordinate frame.
 */
Matrix& LTC(double lon, double lat);

#endif