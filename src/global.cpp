/**
 * @file global.cpp
 * @brief Implementation of global variables and functions to load external data files.
 *
 * This file defines global matrices and functions declared in global.hpp, which are
 * used throughout the project to load essential parameters for orbital computations.
 *
 * It implements reading from external data files located in "../data/", including:
 *  - Earth Orientation Parameters (EOP).
 *  - Normalized spherical harmonic gravitational coefficients (GGM03S model).
 *  - Planetary ephemeris Chebyshev coefficients (DE430Coeff).
 *  - Initialization of global auxiliary parameters.
 *
 * The global matrices eopdata, Cnm, Snm, and PC are defined and populated here.
 */

#include "../include/global.hpp"

Matrix eopdata;

void eop19620101(int c){
    eopdata = zeros(13, c);
    FILE *fid = fopen("../data/eop19620101.txt","r");
    if(fid == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    for(int j=1; j<=c; j++){
        fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,j)),&(eopdata(2,j)),
        &(eopdata(3,j)),&(eopdata(4,j)),&(eopdata(5,j)),&(eopdata(6,j)),&(eopdata(7,j)),&(eopdata(8,j)),
        &(eopdata(9,j)),&(eopdata(10,j)),&(eopdata(11,j)),&(eopdata(12,j)),&(eopdata(13,j)));
    }

    fclose(fid);
}

Matrix Cnm;
Matrix Snm;
void GGM03S(int n){    
    Cnm = zeros(n,n);
    Snm = zeros(n,n);
    FILE *fid = fopen("../data/GGM03S.txt","r");
    if(fid == NULL){
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=i; j++){
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",&aux,&aux,&(Cnm(i,j)),&(Snm(i,j)),&aux,&aux);
        }
    }

    fclose(fid);
}

Matrix PC;
void DE430Coeff(int f, int c){
    PC = zeros(f,c);

    FILE *fid = fopen("../data/DE430Coeff.txt","r");
    if(fid == NULL){
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    for(int i=1; i<=f; i++){
        for(int j=1; j<=c; j++){
            fscanf(fid,"%lf",&(PC(i,j)));
        }
    }

    fclose(fid);
}

Param AuxParam;
void initializeAuxParam(){
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.Mjd_TT = 49746.1170623147;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
}