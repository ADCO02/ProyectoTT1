/**
 * @file sign_.cpp
 * @brief Implementation of the sign_ function.
 */

#include "../include/sign_.hpp"

double sign_(double a, double b){
    if (b>=0.0){
        return fabs(a);
    }else{
        return - fabs(a);
    }
}