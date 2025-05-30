/**
 * @file NutAngles.hpp
 * 
 * @brief Computes the nutation angles in longitude and obliquity.
 *
 * This function calculates the nutation angles (Δψ and Δε) in radians
 * for a given Terrestrial Time (TT) expressed in Modified Julian Date (MJD).
 */

#ifndef _NUTANGLES_
#define _NUTANGLES_

#include "SAT_Const.hpp"
#include <cmath>
#include <tuple>

using namespace std;

/**
 * @brief Computes the nutation angles in longitude and obliquity.
 * 
 * @param Mjd_TT Terrestrial Time expressed in Modified Julian Date.
 * 
 * @return A std::tuple containing:
 *         - nutation in longitude (Δψ) in radians,
 *         - nutation in obliquity (Δε) in radians.
 */
tuple<double,double> NutAngles(double Mjd_TT);

#endif