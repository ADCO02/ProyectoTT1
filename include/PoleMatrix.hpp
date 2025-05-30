/**
 * @file PoleMatrix.hpp
 * 
 * @brief Computes the pole matrix for transforming pseudo Earth-fixed to true Earth-fixed coordinates.
 *
 * Generates the rotation matrix that corrects for polar motion, converting coordinates from 
 * the pseudo Earth-fixed frame (which assumes a fixed pole) to the actual Earth-fixed frame 
 * accounting for the small variations in Earth's rotation axis due to polar motion.
 * This correction is essential in precise geodetic and astronomical calculations.
 */

#ifndef _POLEMATRIX_
#define _POLEMATRIX_

#include "matrix.hpp"
#include "R_x.hpp"
#include "R_y.hpp"

/**
 * @brief Computes the pole matrix for transforming pseudo Earth-fixed to true Earth-fixed coordinates.
 *
 * @param xp Polar motion coordinate in radians (rotation about the y-axis).
 * @param yp Polar motion coordinate in radians (rotation about the x-axis).
 * 
 * @return Reference to the 3×3 pole matrix.
 */
Matrix& PoleMatrix(double xp, double yp);

#endif