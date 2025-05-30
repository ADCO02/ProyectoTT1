/**
 * @file unit.hpp
 * 
 * @brief Computes the unit vector of a given vector.
 * 
 * This function calculates the unit vector corresponding to the input vector.
 * If the input vector is a zero vector (or its magnitude is very small), the
 * output is set to the zero vector.
 */

#ifndef _UNIT_
#define _UNIT_

#include "matrix.hpp"

/**
 * @brief Computes the unit vector of a given vector.
 * 
 * @param vec Input vector.
 * 
 * @return Reference to the unit vector.
 */
Matrix& unit(Matrix& vec);

#endif