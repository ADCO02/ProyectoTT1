/**
 * @file EccAnom.hpp
 * 
 * @brief Computes the eccentric anomaly for elliptic orbits.
 *
 * Solves Kepler's equation using iterative Newton-Raphson method 
 * to find the eccentric anomaly corresponding to a given mean anomaly.
 */

#ifndef _ECCANOM_
#define _ECCANOM_

#include "SAT_Const.hpp"
#include <cfloat>
#include <cmath>
#include <iostream>

using namespace std;

/** 
 * @brief Computes the eccentric anomaly for elliptic orbits.
 *
 * @param M Mean anomaly in radians.
 * @param e Eccentricity of the orbit (range: [0, 1]).
 * 
 * @return Eccentric anomaly in radians.
 */
double EccAnom(double M, double e);

#endif