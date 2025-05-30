#include "..\include\sign_.hpp"

/**
 * @brief Returns the absolute value of a with the sign of b.
 *
 * @param a The value whose magnitude is used.
 * @param b The value whose sign is used.
 * @return The value of a with the sign of b.
 */
double sign_(double a, double b){
    if (b>=0.0){
        return fabs(a);
    }else{
        return - fabs(a);
    }
}