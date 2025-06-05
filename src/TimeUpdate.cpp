/**
 * @file TimeUpdate.cpp
 * @brief Implementation of the TimeUpdate function.
 */

#include "../include/TimeUpdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix& Phi, double Qdt){
    return Phi*P*Phi.transpose() + Qdt;
}
