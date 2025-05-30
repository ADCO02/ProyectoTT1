/**
 * @file VarEqn.hpp
 * 
 * @brief Computes the variational equations for satellite dynamics.
 *
 * Calculates the time derivative of the combined state vector and 
 * the state transition matrix.
 */

#ifndef _VAREQN_
#define _VAREQN_

#include "matrix.hpp"
#include "IERS.hpp"
#include "global.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "SAT_Const.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "AccelHarmonic.hpp"
#include "G_AccelHarmonic.hpp"

/**
 * @brief Computes the derivative of the combined state vector and state transition matrix.
 *
 * @param x      Time since epoch in seconds.
 * @param yPhi   Vector containing the state vector (position & velocity) and the state transition matrix.
 * 
 * @return Reference to the derivative vector of yPhi.
 */
Matrix& VarEqn(double x, Matrix& yPhi);

#endif