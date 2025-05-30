/**
 * @file R_z.hpp
 * 
 * @brief Computes the rotation matrix for a rotation about the z-axis.
 *
 * Generates a 3×3 rotation matrix representing a right-handed rotation 
 * around the z-axis by the specified angle.
 */

#ifndef _R_Z_
#define _R_Z_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Computes the rotation matrix for a rotation about the z-axis.
 *
 * @param angle Rotation angle in radians.
 * 
 * @return Reference to a 3×3 rotation matrix.
 */
Matrix& R_z(double angle);

#endif