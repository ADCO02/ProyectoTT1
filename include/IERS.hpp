/**
 * @file IERS.hpp
 * 
 * @brief Retrieves Earth orientation parameters (EOP) from IERS data.
 * 
 * Returns polar motion coordinates, UT1-UTC, length of day,
 * nutation corrections, and time offsets for a given Modified Julian Date (UTC).
 * It can perform linear interpolation between data points or return exact values.
 */

#ifndef _IERS_
#define _IERS_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "global.hpp"
#include <cmath>
#include <tuple>

/**
 * @brief Retrieves Earth orientation parameters (EOP) from IERS data.
 * 
 * @param Mjd_UTC Modified Julian Date (UTC) for which EOP values are requested.
 * @param interp  Interpolation flag:
 *                - 'n' : no interpolation (default),
 *                - 'l' : linear interpolation between nearest entries.
 * 
 * @return std::tuple containing:
 *         - x_pole   (rad)   : x coordinate of pole position,
 *         - y_pole   (rad)   : y coordinate of pole position,
 *         - UT1_UTC  (s)     : difference between UT1 and UTC,
 *         - LOD      (s)     : length of day,
 *         - dpsi     (rad)   : nutation in longitude,
 *         - deps     (rad)   : nutation in obliquity,
 *         - dx_pole  (rad)   : correction to x pole coordinate,
 *         - dy_pole  (rad)   : correction to y pole coordinate,
 *         - TAI_UTC  (s)     : difference between TAI and UTC.
 * 
 * @note The function uses the global `eopdata` matrix for Earth orientation parameters.
 *       Exits the program if the MJD is not found in the data.
 */
tuple<double,double,double,double,double,double,double,double,double> IERS(double Mjd_UTC, char interp = 'n');

#endif