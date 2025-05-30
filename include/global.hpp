/**
 * @file global.hpp
 * 
 * @brief Global data structures and functions to load external data files.
 * 
 * This header declares global matrices and parameters used throughout the
 * project, and functions to initialize/load Earth Orientation Parameters (EOP),
 * gravitational coefficients (GGM03S), planetary ephemeris coefficients (DE430Coeff),
 * and auxiliary parameters.
 * 
 * These functions read data files from the ../data/ directory and fill the global
 * matrices used in orbital computations.
 */

#ifndef _GLOBAL_
#define _GLOBAL_

#include "matrix.hpp"
#include <cstdio>
#include <cstdlib>

/** 
 * @brief Global matrix holding Earth Orientation Parameters loaded from file.
 * 
 * Dimensions: 13 × number_of_records.
 */
extern Matrix eopdata;

/**
 * @brief Loads Earth Orientation Parameters (EOP) from the file "eop19620101.txt".
 * 
 * @param c Number of records (columns) to read.
 * 
 * This function reads EOP data and fills the global matrix `eopdata` accordingly.
 * The data is expected to be stored in "../data/eop19620101.txt".
 */
void eop19620101(int c);

/**
 * @brief Global matrices holding normalized spherical harmonic coefficients.
 * 
 * - Cnm: cosine coefficients matrix
 * - Snm: sine coefficients matrix
 */
extern Matrix Cnm;
extern Matrix Snm;

/**
 * @brief Loads gravitational coefficients (GGM03S model) from file "GGM03S.txt".
 * 
 * @param n Degree and order of coefficients to load.
 * 
 * The coefficients are loaded into the global matrices `Cnm` and `Snm`.
 */
void GGM03S(int n);

/**
 * @brief Global matrix holding planetary ephemeris Chebyshev coefficients.
 */
extern Matrix PC;

/**
 * @brief Loads planetary ephemeris Chebyshev coefficients from "DE430Coeff.txt".
 * 
 * @param f Number of rows to read.
 * @param c Number of columns to read.
 * 
 * Fills the global matrix `PC`.
 */
void DE430Coeff(int f, int c);

/**
 * @brief Struct holding auxiliary parameters used across computations.
 */
struct Param {
    double Mjd_UTC;
    double Mjd_TT;
    int n;
    int m;
    int sun;
    int moon;
    int planets;
};

/**
 * @brief Global auxiliary parameters variable.
 */
extern Param AuxParam;

/**
 * @brief Initializes the global auxiliary parameters with default values.
 */
void initializeAuxParam();

#endif
