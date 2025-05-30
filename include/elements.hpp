/**
 * @file elements.hpp
 * 
 * @brief Computes the osculating Keplerian elements from a satellite state vector.
 * 
 * This module calculates the classical Keplerian orbital elements — such as 
 * semilatus rectum, semi-major axis, eccentricity, inclination, longitude of 
 * ascending node, argument of perihelion, and mean anomaly — from the satellite's 
 * position and velocity vectors. It assumes elliptical orbits and is not designed 
 * for circular or equatorial orbits.
 */

#ifndef _ELEMENTS_
#define _ELEMENTS_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <math.h>
#include <tuple>

/**
 * @brief Computes the osculating Keplerian elements from a satellite state vector.
 *
 * @param y State vector containing position (x, y, z) and velocity (vx, vy, vz).
 * 
 * @return Tuple containing:
 *  - semilatus rectum p [m]
 *  - semi-major axis a [m]
 *  - eccentricity e [-]
 *  - inclination i [rad]
 *  - longitude of ascending node Omega [rad]
 *  - argument of perihelion omega [rad]
 *  - mean anomaly M [rad]
 */
tuple<double, double, double, double, double, double, double> elements(Matrix& y);

#endif