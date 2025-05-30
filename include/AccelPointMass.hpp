/**
 * @file AccelPointMass.hpp
 * @brief Computes the perturbational acceleration due to a point mass.
 *
 * Calculates the acceleration experienced by a satellite 
 * due to the gravitational attraction of a point mass.
 */

#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Computes the perturbational acceleration due to a point mass.
 * 
 * @param r Satellite position vector (relative to an inertial frame).
 * @param s Point mass position vector (relative to the same frame as r).
 * @param GM Gravitational parameter (G * mass) of the point mass.
 * 
 * @return Reference to a Matrix representing the resulting acceleration vector.
 */
Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM);

#endif