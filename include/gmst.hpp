/**
 * @file gmst.hpp
 * 
 * @brief Computes Greenwich Mean Sidereal Time (GMST) from UT1.
 *
 * Returns the Greenwich Mean Sidereal Time (GMST) in radians, based on the 
 * input Modified Julian Date in UT1 time scale. GMST is essential for transforming 
 * between inertial and Earth-fixed coordinate systems, as it reflects the rotation 
 * of the Earth relative to the fixed stars.
 */

#ifndef _GMST_
#define _GMST_

#include "SAT_Const.hpp"
#include "Frac.hpp"

/**
 * @brief Computes Greenwich Mean Sidereal Time (GMST) from UT1.
 *
 * @param Mjd_UT1 Modified Julian Date (UT1 time scale).
 * 
 * @return GMST in radians within the range [0, 2π).
 */
double gmst(double Mjd_UT1);

#endif