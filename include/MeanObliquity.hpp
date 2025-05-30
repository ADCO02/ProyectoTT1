/**
 * @file MeanObliquity.hpp
 * 
 * @brief Computes the mean obliquity of the ecliptic.
 *
 * Calculates the mean obliquity based on the Modified Julian Date in Terrestrial Time (TT),
 * using a standard IAU polynomial approximation.
 */

#ifndef _MEANOBLIQUITY_
#define _MEANOBLIQUITY_

#include "SAT_Const.hpp"

/**
 * @brief Computes the mean obliquity of the ecliptic.
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * 
 * @return Mean obliquity of the ecliptic in radians.
 */
double MeanObliquity(double Mjd_TT);

#endif