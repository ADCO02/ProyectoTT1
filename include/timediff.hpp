#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <tuple>

using namespace std;

/**
 * @brief Computes several time differences between different time scales.
 * 
 * @param UT1_UTC Difference between UT1 and UTC in seconds.
 * @param TAI_UTC Difference between TAI and UTC in seconds.
 * @return A tuple containing:
 *         - UT1_TAI: difference between UT1 and TAI [s]
 *         - UTC_GPS: difference between UTC and GPS [s]
 *         - UT1_GPS: difference between UT1 and GPS [s]
 *         - TT_UTC:  difference between TT and UTC [s]
 *         - GPS_UTC: difference between GPS and UTC [s]
 */
tuple<double,double,double,double,double> timediff(double UT1_UTC, double TAI_UTC);

#endif