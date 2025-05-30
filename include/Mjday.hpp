/**
 * @file Mjday.hpp
 * 
 * @brief Computes the Modified Julian Date from calendar date and time.
 *
 * Converts a calendar date and time (UTC) into the Modified Julian Date (MJD).
 */

#ifndef _MJDAY_
#define _MJDAY_

#include <cmath>

/**
 * @brief Computes the Modified Julian Date from calendar date and time.
 *
 * @param yr  Year.
 * @param mon Month.
 * @param day Day.
 * @param hr  Hour (default is 0).
 * @param min Minute (default is 0).
 * @param sec Second (default is 0).
 * 
 * @return Modified Julian Date.
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0);

#endif