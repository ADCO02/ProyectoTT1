#include "..\include\AzElPa.hpp"

tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s){ 
    Matrix& dAds = zeros(3);
    Matrix& dEds = zeros(3);

    double pi2 = 2.0*M_PI;

    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    double Az = atan2(s(1),s(2));

    if (Az<0.0){
        Az = Az+pi2;
    }

    double El = atan ( s(3) / rho );

    // Partials
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;

    dEds(1) = -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    dEds(3) = rho;
    dEds = dEds/dot(s,s);

    return tie(Az, El, dAds, dEds);
}