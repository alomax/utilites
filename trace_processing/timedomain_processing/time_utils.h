/*
 * This file is part of the Anthony Lomax Java Library.
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

/*
 * File:   time_utils.h
 * Author: anthony
 *
 * Created on 22 December 2017, 13:37
 */

#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>
#include "../miniseed_utils/miniseed_utils.h"

#define DEFAULT_TIME_FORMAT 0
#define IRIS_WS_TIME_FORMAT 1       // 2 decimal places for seconds
#define IRIS_WS_TIME_FORMAT_4 2     // 4 decimal places for seconds
#define COMMA_DELIMTED_TIME_FORMAT 3


char *timeDiff2string(time_t later_time, time_t earlier_time, char* tdiff_string);
char *time2string(double time_dec_sec, char* str);
char *timeDecSec2string(double time_dec_sec, char* str, int iformat);
double string2timeDecSec(char *str);
hptime_t timeDecSec2hptime(double time_dec_sec);
void hptime2dateTimeComponents(hptime_t hptime, int *pyear, int *pjday, int *phour, int *pmin, double *psec);

#ifdef __cplusplus
}
#endif

#endif /* TIME_UTILS_H */

