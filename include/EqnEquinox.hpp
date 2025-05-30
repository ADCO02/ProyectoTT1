/**
 * @file EqnEquinox.hpp
 * 
 * @brief Computes the equation of the equinoxes from nutation and obliquity values.
 * 
 * This function calculates the equation of the equinoxes, defined as the product of 
 * the nutation in longitude and the cosine of the mean obliquity. It represents the 
 * difference between apparent and mean sidereal time, or equivalently, the right ascension 
 * of the mean equinox referred to the true equator and equinox.
 */

#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_

#include "NutAngles.hpp"
#include "MeanObliquity.hpp"

/**
 * @brief Computes the equation of the equinoxes from nutation and obliquity values.
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * 
 * @return The equation of the equinoxes [rad].
 */
double EqnEquinox(double Mjd_TT);

#endif