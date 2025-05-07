#include "..\include\AccelPointMass.hpp"

Matrix AccelPointMass(Matrix &r, Matrix &s, double GM){
    // Relative position vector of satellite w.r.t. point mass 
    Matrix& d = r - s;

    // Acceleration
    Matrix& a = ( d/pow(norm(d),3) + s/pow(norm(s),3) ) * -GM;

    return a;

}