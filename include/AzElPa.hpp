#ifndef _AZELPA_
#define _AZELPA_
#define _USE_MATH_DEFINES

#include "matrix.hpp"
#include <math.h>

//dads y deds deben inicializarse en 3 antes de pasarlos
void AzElPa(Matrix &s, double* Az, double* El, Matrix &dAds, Matrix &dEds);

#endif