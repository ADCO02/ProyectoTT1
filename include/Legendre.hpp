/**
 * @file Legendre.hpp
 * 
 * @brief Computes associated Legendre polynomials and their derivatives.
 *
 * Given degree n, order m, and angle fi (in radians), this function computes
 * the normalized associated Legendre functions \( P_{nm}(\phi) \) and their
 * derivatives with respect to \(\phi\).
 */

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "matrix.hpp"
#include <cmath>
#include<tuple>

/**
 * @brief Computes associated Legendre polynomials and their derivatives.
 *
 * @param n Degree of the Legendre polynomial (non-negative integer).
 * @param m Order of the Legendre polynomial (0 ≤ m ≤ n).
 * @param fi Angle in radians (colatitude or similar).
 *
 * @return std::tuple containing references to two matrices (n+1)×(m+1):
 *         - pnm: matrix of normalized associated Legendre polynomials \( P_{nm}(\phi) \)
 *         - dpnm: matrix of their derivatives with respect to \(\phi\).
 */
tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi);

#endif