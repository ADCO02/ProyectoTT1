#include "..\include\MeasUpdate.hpp"

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n){
    int m = 1;
    double Inv_W = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    Matrix& K = (P*(G.transpose())*inv(G*P*(G.transpose())+Inv_W)).transpose();

    // State update
    x = x + (K*(z-g));

    // Covariance update
    P = (eye(n)-K.transpose()*G)*P;

    return tie(K,x,P);
}