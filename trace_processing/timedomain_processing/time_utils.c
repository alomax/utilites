/*
 * This file is part of the Anthony Lomax C Library.
 *
 * Copyright (C) 2017 Anthony Lomax <anthony@alomax.net www.alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#include "time_utils.h"

/** format two integer second times into pretty hour-min-sec string */

char *timeDiff2string(time_t later_time, time_t earlier_time, char* tdiff_string) {

    double tdiff = difftime(later_time, earlier_time);

    int tdiff_tenths_sec = (int) (0.5 + 10.0 * tdiff);
    int tdiff_hour = 0;
    int tdiff_min = 0;
    while (tdiff_tenths_sec >= 36000) {
        tdiff_tenths_sec -= 36000;
        tdiff_hour++;
    }
    while (tdiff_tenths_sec >= 600) {
        tdiff_tenths_sec -= 600;
        tdiff_min++;
    }
    if (tdiff_hour > 0)
        sprintf(tdiff_string, "%dh %2.2dm %2.2ds", tdiff_hour, tdiff_min, tdiff_tenths_sec / 10);
    else if (tdiff_min > 0)
        sprintf(tdiff_string, "%dm %2.2ds", tdiff_min, tdiff_tenths_sec / 10);
    else if (tdiff_tenths_sec > 0)
        sprintf(tdiff_string, "%ds", tdiff_tenths_sec / 10);
    else
        sprintf(tdiff_string, "0s");

    return (tdiff_string);

}

/** format a time_t time into time string with integer seconds */

char *time2string(double time_dec_sec, char* str) {

    // round origin time to nearest second
    time_t itime_sec = (time_t) (0.5 + time_dec_sec); // 20110930 AJL - modified to avoid times that round to 60s
    struct tm* tm_time = gmtime(&itime_sec);

    sprintf(str, "%4.4d.%2.2d.%2.2d-%2.2d:%2.2d:%2.2d",
            tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, tm_time->tm_sec);

    return (str);
}

/** format a double time (equivalent to time_t time with decimal seconds) into a time string with decimal seconds */

char *timeDecSec2string(double time_dec_sec, char* str, int iformat) {

    time_t itime_sec = (time_t) time_dec_sec;
    struct tm* tm_time = gmtime(&itime_sec);

    double dec_sec = time_dec_sec - itime_sec;

    if (iformat == IRIS_WS_TIME_FORMAT) {   // 2 decimal places for seconds
        sprintf(str, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    } else if (iformat == IRIS_WS_TIME_FORMAT_4) {  // 4 decimal places for seconds
        sprintf(str, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%07.4f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    } else if (iformat == COMMA_DELIMTED_TIME_FORMAT) {
        sprintf(str, "%4.4d,%2.2d,%2.2d,%2.2d,%2.2d,%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    } else {
        sprintf(str, "%4.4d.%2.2d.%2.2d-%2.2d:%2.2d:%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    }

    return (str);
}

/** convert a time string into a double time (equivalent to time_t time with decimal seconds) */

double string2timeDecSec(char *str) {

    //printf("\nDEBUG: string2timeDecSec %s\n", str);

    struct tm tm_time;

    double dec_sec;
    sscanf(str, "%4d.%2d.%2d-%d:%d:%lf", &(tm_time.tm_year), &(tm_time.tm_mon), &(tm_time.tm_mday), &(tm_time.tm_hour), &(tm_time.tm_min), &dec_sec);
    tm_time.tm_year -= 1900;
    tm_time.tm_mon -= 1;
    tm_time.tm_sec = (int) dec_sec;

    //printf("DEBUG: string2timeDecSec %.4d.%.2d.%.2d-%.2d:%.2d:%.2d\n", tm_time.tm_year + 1900, tm_time.tm_mon + 1, tm_time.tm_mday, tm_time.tm_hour, tm_time.tm_min, tm_time.tm_sec);

    time_t time_gmt = timegm(&tm_time);
    //printf("DEBUG: time_gmt %ld\n", time_gmt);

    double dec_time = (double) time_gmt + dec_sec - (int) dec_sec;
    //printf("DEBUG: dec_time %lf\n", dec_time);

    return (dec_time);
}

/** convert a double time (equivalent to time_t time with decimal seconds) into a libmseed hptime_t */

hptime_t timeDecSec2hptime(double time_dec_sec) {

    time_t itime_sec = (time_t) time_dec_sec;
    struct tm* tm_time = gmtime(&itime_sec);

    double dec_sec = time_dec_sec - itime_sec;
    int year = tm_time->tm_year + 1900;
    int month = tm_time->tm_mon + 1;
    int mday = tm_time->tm_mday;
    int jday;
    ms_md2doy(year, month, mday, &jday);
    return (ms_time2hptime(year, jday, tm_time->tm_hour, tm_time->tm_min, tm_time->tm_sec, (int) (0.5 + (dec_sec * (double) 1000000))));

}

/** convert a double time (equivalent to time_t time with decimal seconds) into component date/time elements */

void hptime2dateTimeComponents(hptime_t hptime, int *pyear, int *pjday, int *phour, int *pmin, double *psec) {

    BTime btime;
    ms_hptime2btime(hptime, &btime);

    *pyear = btime.year;
    *pjday = btime.day;
    *phour = btime.hour;
    *pmin = btime.min;
    *psec = (double) btime.sec + ((double) btime.fract / (double) 10000);

}
