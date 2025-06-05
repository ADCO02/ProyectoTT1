/**
 * @file unit.cpp
 * @brief Implementation of the unit function.
 */

#include "../include/unit.hpp"

Matrix& unit(Matrix& vec){
    double small = 0.000001;
    double magv = norm(vec);

    Matrix& outvec = zeros(3);
    if ( magv > small ){
        for (int i=1;i<=3;i++){
            outvec(i)= vec(i)/magv;
        }
    }else{
        for (int i=1;i<=3;i++){
            outvec(i)= 0.0;
        }
    }

    return outvec;
}