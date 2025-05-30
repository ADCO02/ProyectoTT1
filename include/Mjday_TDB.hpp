/**
 * @file Mjday_TDB.hpp
 * 
 * @brief Computes the Modified Julian Date in Barycentric Dynamical Time (TDB).
 *
 * Converts the Modified Julian Date in Terrestrial Time (TT) to 
 * Barycentric Dynamical Time (TDB) using a series approximation.
 */

#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_

#include <cmath>

/**
 * @brief Computes the Modified Julian Date in Barycentric Dynamical Time (TDB).
 *
 * @param Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 * 
 * @return Modified Julian Date in Barycentric Dynamical Time (TDB).
 */
double Mjday_TDB(double Mjd_TT);

#endif