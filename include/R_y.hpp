/**
 * @file R_y.hpp
 * 
 * @brief Computes the rotation matrix for a rotation about the y-axis.
 *
 * Generates a 3×3 rotation matrix representing a right-handed rotation 
 * around the y-axis by the specified angle.
 */

#ifndef _R_Y_
#define _R_Y_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Computes the rotation matrix for a rotation about the y-axis.
 *
 * @param angle Rotation angle in radians.
 * 
 * @return Reference to a 3×3 rotation matrix.
 */
Matrix& R_y(double angle);

#endif