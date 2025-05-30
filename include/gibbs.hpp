/**
 * @file gibbs.hpp
 * 
 * @brief Performs the Gibbs method of orbit determination to compute velocity.
 * 
 * This function implements the Gibbs method to determine the velocity vector 
 * at the middle of three given position vectors in space. It uses vector 
 * operations to check coplanarity and consistency of the orbit and returns 
 * the computed velocity along with angles between vectors and an error status.
 */

#ifndef _GIBBS_
#define _GIBBS_

#include "matrix.hpp"
#include "angl.hpp"
#include "SAT_Const.hpp"
#include "unit.hpp"
#include <tuple>

/**
 * @brief Performs the Gibbs method of orbit determination to compute velocity.
 *
 * @param r1 Position vector #1 (meters).
 * @param r2 Position vector #2 (meters).
 * @param r3 Position vector #3 (meters).
 * 
 * @return Tuple containing:
 *  - Velocity vector at r2 (3×1 Matrix) [m/s]
 *  - Angle between r1 and r2 vectors (radians)
 *  - Angle between r2 and r3 vectors (radians)
 *  - Coplanarity angle (radians)
 *  - Error status string ("ok", "not coplanar", "impossible")
 */
tuple<Matrix&, double, double, double, string> gibbs(Matrix& r1, Matrix& r2, Matrix& r3);

#endif