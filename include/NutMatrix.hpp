/**
 * @file NutMatrix.hpp
 * 
 * @brief Computes the nutation matrix from mean to true equator and equinox.
 *
 * Constructs the transformation matrix that accounts for the effects of nutation, 
 * converting coordinates from the mean equator and equinox of date to the true 
 * equator and equinox. This transformation is essential for high-precision 
 * astronomical calculations, as it includes periodic variations caused by the 
 * gravitational influence of the Moon and Sun on Earth's equatorial bulge.
 */

#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "matrix.hpp"
#include "NutAngles.hpp"
#include "MeanObliquity.hpp"
#include "R_x.hpp"
#include "R_z.hpp"

/**
 * @brief Computes the nutation matrix from mean to true equator and equinox.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * 
 * @return Reference to the 3×3 nutation matrix.
 */
Matrix& NutMatrix(double Mjd_TT);

#endif