#include "..\include\AzElPa.hpp"

void AzElPa(Matrix &s, double* Az, double* El, Matrix &dAds, Matrix &dEds){
    if (dAds.n_row != 1 || dAds.n_column != 3)
	{
		cout << "AzElPa: dAds must be a size 3 vector\n";
		exit(EXIT_FAILURE);
	}

    if (dEds.n_row != 1 || dEds.n_column != 3)
	{
		cout << "AzElPa: dEds must be a size 3 vector\n";
		exit(EXIT_FAILURE);
	}

    double pi2 = 2.0*M_PI;

    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    *Az = atan2(s(1),s(2));

    if (*Az<0.0){
        *Az = *Az+pi2;
    }

    *El = atan ( s(3) / rho );

    // Partials
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;

    dEds(1) = -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    dEds(3) = rho;
    dEds = dEds/dot(s,s);
}