/**
 * @file GHAMatrix.hpp
 * 
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system.
 *
 * Generates the Greenwich Hour Angle (GHA) rotation matrix based on the 
 * Greenwich Apparent Sidereal Time at the given UT1 date.
 */

#ifndef _GHAMATRIX_
#define _GHAMATRIX_

#include "matrix.hpp"
#include "R_z.hpp"
#include "gast.hpp"

/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system.
 *
 * @param Mjd_UT1  Modified Julian Date in UT1 time scale.
 * 
 * @return Reference to the 3×3 Greenwich Hour Angle rotation matrix.
 */
Matrix& GHAMatrix(double Mjd_UT1);

#endif