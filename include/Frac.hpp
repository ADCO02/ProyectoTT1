/**
 * @file Frac.hpp
 * 
 * @brief Returns the fractional part of a number.
 *
 * Computes the fractional part of a real number, defined as \f$ y = x - \lfloor x \rfloor \f$.
 */

#ifndef _FRAC_
#define _FRAC_

#include <cmath>

/**
 * @brief Returns the fractional part of a number.
 *
 * @param x Input real number.
 * 
 * @return Fractional part of the input.
 */
double Frac(double x);

#endif