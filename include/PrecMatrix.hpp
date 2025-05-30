/**
 * @file PrecMatrix.hpp
 * 
 * @brief Computes the precession matrix between two epochs in Terrestrial Time (TT).
 *
 * Calculates the 3×3 rotation matrix that accounts for the precession of the Earth's 
 * rotation axis between two epochs, expressed as Modified Julian Dates in TT. 
 * Precession reflects the long-term, gradual change in the orientation of the Earth's 
 * rotational axis due to gravitational torques, and is crucial for transforming 
 * equatorial coordinates between different times.
 */

#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "matrix.hpp"
#include "R_x.hpp"
#include "R_y.hpp"
#include "R_z.hpp"
#include "SAT_Const.hpp"

/**
 * @brief Computes the precession matrix between two epochs in Terrestrial Time (TT).
 *
 * @param Mjd_1 Initial epoch as Modified Julian Date (TT).
 * @param Mjd_2 Target epoch to precess to, as Modified Julian Date (TT).
 * 
 * @return Reference to the 3×3 precession matrix.
 */
Matrix& PrecMatrix(double Mjd_1, double Mjd_2);

#endif