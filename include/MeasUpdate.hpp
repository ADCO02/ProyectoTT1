/**
 * @file MeasUpdate.hpp
 * 
 * @brief Performs a Kalman filter measurement update step.
 *
 * Implements the measurement update step of a Kalman filter. Given the 
 * prior state and covariance, along with a new measurement and its expected 
 * value and sensitivity, it computes the Kalman gain, updates the state vector 
 * and covariance matrix accordingly.
 */

#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include <tuple>
#include "matrix.hpp"

using namespace std;

/**
 * @brief Performs a Kalman filter measurement update step.
 *
 * @param x  State vector (will be updated).
 * @param z  Measurement scalar.
 * @param g  Expected measurement scalar (based on predicted state).
 * @param s  Measurement noise standard deviation.
 * @param G  Measurement sensitivity matrix.
 * @param P  State covariance matrix (will be updated).
 * @param n  Dimension of the state vector.
 * 
 * @return Tuple containing references to:
 *         - Kalman gain matrix K,
 *         - Updated state vector x,
 *         - Updated covariance matrix P.
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n);

#endif