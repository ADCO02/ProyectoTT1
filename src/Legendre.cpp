/**
 * @file Legendre.cpp
 * @brief Implementation of the Legendre function.
 */

#include "../include/Legendre.hpp"

tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi){
    Matrix& pnm = zeros(n+1,m+1);
    Matrix& dpnm = zeros(n+1,m+1);
    pnm(1,1)=1;
    dpnm(1,1)=0;
    pnm(2,2)=sqrt(3)*cos(fi);
    dpnm(2,2)=-sqrt(3)*sin(fi);
    // diagonal coefficients
    for(int i=2; i<=n ; i++){    
        pnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*cos(fi)*pnm(i,i);
    }
    for(int i=2; i<=n ; i++){
        dpnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*((cos(fi)*dpnm(i,i))-(sin(fi)*pnm(i,i)));
    }
    // horizontal first step coefficients
    for(int i=1; i<=n ; i++){
        pnm(i+1,i)= sqrt(2.0*i+1)*sin(fi)*pnm(i,i);
    }
    for(int i=1; i<=n ; i++){
        dpnm(i+1,i)= sqrt(2.0*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
    }
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(1){
        for(int i=k; i<=n;i++){
            pnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*pnm(i,j+1))
                -(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*pnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m) break;
    }
    j = 0;
    k = 2;
    while(1){
        for(int i=k; i<=n;i++){
            dpnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*dpnm(i,j+1))
                +(sqrt(2.0*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m) break;
    }
    return tie(pnm,dpnm);
}