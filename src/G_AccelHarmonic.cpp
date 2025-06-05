/**
 * @file G_AccelHarmonic.cpp
 * @brief Implementation of the G_AccelHarmonic function.
 */

#include "../include/G_AccelHarmonic.hpp"

Matrix& G_AccelHarmonic( Matrix& r, Matrix& U, int n_max, int m_max ){
    double d = 1.0;   // Position increment [m]

    Matrix& G = zeros(3,3);
    Matrix& dr = zeros(3);

    // Gradient
    for (int i=1; i<=3; i++){
        // Set offset in i-th component of the position vector
        dr = zeros(3);
        dr(i) = d;

        // Acceleration difference
        Matrix& da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) - AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis
        G.assign_column(i, da/d);      
    }

    return G;
}