/**
 * @file angl.hpp
 * 
 * @brief Computes the angle between two vectors.
 */

#ifndef _ANGL_
#define _ANGL_

#include "matrix.hpp"

/**
 * @brief Calculates the angle between two vectors in radians.
 * 
 * @param vec1 First vector.
 * @param vec2 Second vector.
 * @return Angle between vec1 and vec2 in radians, from 0 to pi.
 *         Returns a large undefined value if any vector is near zero length.
 */
double angl(Matrix& vec1, Matrix& vec2);

#endif