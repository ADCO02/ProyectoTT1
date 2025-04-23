#include "..\include\TimeUpdate.hpp"

Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt){
    return ((Phi*P)*(Phi.transpose())) + Qdt;
}

Matrix TimeUpdate(Matrix P, Matrix Phi){
    Matrix Qdt = zeros(P.n_row, P.n_column);
    return TimeUpdate(P,Phi,Qdt);
}