/**
 * @file JPL_Eph_DE430.hpp
 * 
 * @brief Provides planetary ephemerides from the JPL DE430 dataset.
 *
 * Computes planetary positions and velocities using the JPL DE430 ephemeris,
 * returning Chebyshev-polynomial-interpolated vectors for solar system bodies 
 * at a given TDB-modified Julian date. This function interfaces with the underlying
 * interpolation and data-loading mechanisms defined in Cheb3D.hpp and global.hpp.
 */

#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include <tuple>
#include "global.hpp"
#include "matrix.hpp"
#include "Cheb3D.hpp"

/**
 * @brief Provides planetary ephemerides from the JPL DE430 dataset.
 *
 * @param Mjd_TDB Modified Julian Date in Barycentric Dynamical Time (TDB).
 * 
 * @return A tuple of references to 11 `Matrix` objects, each representing
 *         the position and/or velocity vectors of specific solar system bodies,
 *         interpolated using Chebyshev polynomials from the DE430 dataset.
 *
 *         The matrices typically correspond to:
 *         Sun, Mercury, Venus, Earth-Moon barycenter, Mars, Jupiter, Saturn,
 *         Uranus, Neptune, Pluto, and Moon or Earth (depending on implementation).
 */
tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif