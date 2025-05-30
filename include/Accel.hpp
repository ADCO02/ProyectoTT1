/**
 * @file Accel.hpp
 * 
 * @brief Computes the acceleration of an Earth orbiting satellite.
 *
 * Calculates the satellite acceleration considering:
 * - Earth's harmonic gravity field,
 * - Gravitational perturbations from the Sun and Moon,
 * - Planetary perturbations from major planets.
 */

#ifndef _ACCEL_
#define _ACCEL_

#include "matrix.hpp"
#include "global.hpp"
#include "IERS.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "SAT_Const.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "Mjday_TDB.hpp"
#include "JPL_Eph_DE430.hpp"
#include "AccelHarmonic.hpp"
#include "AccelPointMass.hpp"

/**
 * @brief Computes the acceleration of an Earth orbiting satellite.
 *
 * @param x  Time offset in seconds since epoch (UTC).
 * @param Y  Satellite state vector (position and velocity) in ICRF/EME2000 system.
 * 
 * @return Reference to the satellite acceleration vector in ICRF/EME2000 system.
 */
Matrix& Accel(double x, Matrix& Y);

#endif