/**
 * @file sign_.hpp
 * 
 * @brief Returns the absolute value of a with the sign of b.
 * 
 * Provides a utility function that returns the magnitude of one value 
 * combined with the sign of another.
 */

#ifndef _SIGN__
#define _SIGN__

#include <cmath>

/**
 * @brief Returns the absolute value of a with the sign of b.
 *
 * @param a The value whose magnitude is used.
 * @param b The value whose sign is used.
 * @return The value of a with the sign of b.
 */
double sign_(double a, double b);

#endif