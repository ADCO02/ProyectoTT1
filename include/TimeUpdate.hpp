/**
 * @file TimeUpdate.hpp
 * 
 * @brief Updates the covariance matrix over a time step using the state transition matrix.
 * 
 * This function performs the time update (also known as the prediction step) 
 * of a covariance matrix P based on the state transition matrix Phi and optionally adds 
 * process noise scaled by Qdt.
 */

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Updates the covariance matrix over a time step using the state transition matrix.
 * 
 * @param P The covariance matrix prior to the update.
 * @param Phi The state transition matrix for the time interval.
 * @param Qdt Optional scalar value representing process noise contribution over the time interval (default is 0.0).
 * 
 * @return Reference to the updated covariance matrix.
 */
Matrix& TimeUpdate(Matrix& P, Matrix& Phi, double Qdt = 0.0);

#endif