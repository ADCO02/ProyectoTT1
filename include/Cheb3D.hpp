#ifndef _CHEB3D_
#define _CHEB3D_

#include "matrix.hpp"

/**
 * @brief Evaluates a Chebyshev approximation for a 3D vector at a given time.
 *
 * This function uses the Clenshaw algorithm to compute the value of a 
 * 3D vector function approximated by Chebyshev polynomials.
 *
 * @param t  Evaluation time (must be in the interval [Ta, Tb]).
 * @param N  Number of Chebyshev coefficients.
 * @param Ta Start of the approximation interval.
 * @param Tb End of the approximation interval.
 * @param Cx Chebyshev coefficients for the x-component.
 * @param Cy Chebyshev coefficients for the y-component.
 * @param Cz Chebyshev coefficients for the z-component.
 * 
 * @return Reference to a Matrix containing the approximated 3D vector at time t.
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);

#endif