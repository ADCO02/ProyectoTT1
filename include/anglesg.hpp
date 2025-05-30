/**
 * @file anglesg.hpp
 * 
 * @brief Orbit determination from three optical sightings using the Gauss method.
 * 
 * This module implements the solution of the classical orbit determination problem 
 * from three optical observations (azimuth and elevation angles) at three different times.
 * It returns the position and velocity vectors in the inertial reference frame at the second observation time.
 * 
 * The implementation is based on the classical angles-only Gaussian method, enhanced 
 * with Earth rotation and observational site position corrections.
 * 
 * The function depends on several auxiliary modules for Earth rotation parameters,
 * coordinate transformations, and orbit mechanics.
 */

#ifndef _ANGLESG_
#define _ANGLESG_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "LTC.hpp"
#include "IERS.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "Geodetic.hpp"
#include "gibbs.hpp"
#include "hgibbs.hpp"
#include "elements.hpp"
#include "rpoly.hpp"
#include <tuple>

/**
 * @brief Solves the orbit determination problem using three optical sightings.
 * 
 * @param az1 Azimuth angle at time t1 [rad]
 * @param az2 Azimuth angle at time t2 [rad]
 * @param az3 Azimuth angle at time t3 [rad]
 * @param el1 Elevation angle at time t1 [rad]
 * @param el2 Elevation angle at time t2 [rad]
 * @param el3 Elevation angle at time t3 [rad]
 * @param Mjd1 Modified Julian Date of the first observation
 * @param Mjd2 Modified Julian Date of the second observation
 * @param Mjd3 Modified Julian Date of the third observation
 * @param Rs1 Position vector of site 1 in IJk frame [m]
 * @param Rs2 Position vector of site 2 in IJk frame [m]
 * @param Rs3 Position vector of site 3 in IJk frame [m]
 * 
 * @return A tuple containing:
 *  - r: Position vector at time t2 in IJk frame [m]
 *  - v: Velocity vector at time t2 in IJk frame [m/s]
 */
tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3,
    double Mjd1, double Mjd2, double Mjd3, Matrix& Rs1, Matrix& Rs2, Matrix& Rs3);

#endif