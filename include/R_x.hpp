/**
 * @file R_x.hpp
 * 
 * @brief Computes the rotation matrix for a rotation about the x-axis.
 *
 * Generates a 3×3 rotation matrix representing a right-handed rotation 
 * around the x-axis by the specified angle.
 */

#ifndef _R_X_
#define _R_X_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Computes the rotation matrix for a rotation about the x-axis.
 *
 * @param angle Rotation angle in radians.
 * 
 * @return Reference to a 3×3 rotation matrix.
 */
Matrix& R_x(double angle);

#endif