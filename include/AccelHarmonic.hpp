/**
 * @file AccelHarmonic.hpp
 * 
 * @brief Computes the acceleration due to the harmonic gravity field of the central body.
 * 
 * This function calculates the gravitational acceleration exerted on a satellite 
 * by the Earth's harmonic gravity field up to specified degree and order.
 * It transforms the satellite position into the body-fixed frame, evaluates 
 * the spherical harmonic terms using associated Legendre polynomials, and returns 
 * the acceleration vector in the inertial frame.
 */

#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "matrix.hpp"
#include "Legendre.hpp"
#include "global.hpp"
#include <cmath>

/**
 * @brief Computes the acceleration due to the harmonic gravity field of the central body.
 * 
 * @param r Satellite position vector in the inertial reference frame.
 * @param E Transformation matrix from inertial to body-fixed frame.
 * @param n_max Maximum degree of the spherical harmonic expansion.
 * @param m_max Maximum order of the spherical harmonic expansion (m_max <= n_max).
 * 
 * @return Reference to a Matrix representing the acceleration vector in the inertial frame.
 */
Matrix& AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif