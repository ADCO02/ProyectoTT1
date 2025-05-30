/**
 * @file hgibbs.hpp
 * 
 * @brief Implements the Herrick-Gibbs approximation for orbit determination.
 * 
 * This function computes the velocity vector at the middle position vector
 * using the Herrick-Gibbs method, based on three position vectors and their
 * corresponding Julian dates. It checks for coplanarity and angle tolerances
 * to ensure the validity of the orbit determination.
 */

#ifndef _HGIBBS_
#define _HGIBBS_

#include "matrix.hpp"
#include "angl.hpp"
#include "SAT_Const.hpp"
#include "unit.hpp"
#include <tuple>

/**
 * @brief Implements the Herrick-Gibbs approximation for orbit determination.
 *
 * @param r1 First position vector (meters).
 * @param r2 Second position vector (meters).
 * @param r3 Third position vector (meters).
 * @param Mjd1 Julian date corresponding to r1 (days).
 * @param Mjd2 Julian date corresponding to r2 (days).
 * @param Mjd3 Julian date corresponding to r3 (days).
 * 
 * @return Tuple containing:
 *  - Velocity vector at r2 (3×1 Matrix) [m/s]
 *  - Angle between r1 and r2 vectors (radians)
 *  - Angle between r2 and r3 vectors (radians)
 *  - Coplanarity angle (radians)
 *  - Error status string ("ok", "not coplanar", "angl > 1ø")
 */
tuple<Matrix&, double, double, double, string> hgibbs(Matrix& r1, Matrix& r2, Matrix& r3,
    double Mjd1, double Mjd2, double Mjd3);

#endif