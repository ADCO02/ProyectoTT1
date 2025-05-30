/**
 * @file gast.hpp
 * 
 * @brief Computes Greenwich Apparent Sidereal Time (GAST) from UT1.
 *
 * Returns the Greenwich Apparent Sidereal Time (GAST) in radians, based on the 
 * input Modified Julian Date in the UT1 time scale. GAST accounts for both the 
 * Earth's rotation (GMST) and the equation of the equinoxes, reflecting the 
 * apparent position of the vernal equinox due to nutation.
 */

#ifndef _GAST_
#define _GAST_

#include "SAT_Const.hpp"
#include "EqnEquinox.hpp"
#include "gmst.hpp"

/**
 * @brief Computes Greenwich Apparent Sidereal Time (GAST) from UT1.
 *
 * @param Mjd_UT1 Modified Julian Date (UT1 time scale).
 * 
 * @return GAST in radians within the range [0, 2π).
 */
double gast(double Mjd_UT1);

#endif