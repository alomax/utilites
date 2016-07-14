/*
 * This file is part of the Anthony Lomax C Library.
 *
 * Copyright (C) 2006-2010 Anthony Lomax <anthony@alomax.net www.alomax.net>
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


void filter_lp_bu_co_5_n_6_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_lp_bu_co_5_n_6_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_lp_bu_co_5_n_6_40sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_lp_bu_co_5_n_6_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_lp_bu_co_5_n_6_100sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_lp_bu_co_5_n_6_200sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================


// Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 6 -a 2.5000000000e-01 0.0000000000e+00 -l */

#define GAIN_lp_bu_co_5_n_6_20sps   3.379723001e+01

void filter_lp_bu_co_5_n_6_20sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_20sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = (xv[0] + xv[6]) + 6 * (xv[1] + xv[5]) + 15 * (xv[2] + xv[4])
                + 20 * xv[3]
                + (-0.0017509260 * yv[0]) + (-0.0000000000 * yv[1])
                + (-0.1141994251 * yv[2]) + (-0.0000000000 * yv[3])
                + (-0.7776959619 * yv[4]) + (0.0000000000 * yv[5]);
        sampleNew[n] = (float) yv[6];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 6 -a 2.0000000000e-01 0.0000000000e+00 -l */

#define GAIN_lp_bu_co_5_n_6_25sps   9.696617316e+01

void filter_lp_bu_co_5_n_6_25sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_25sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = (xv[0] + xv[6]) + 6 * (xv[1] + xv[5]) + 15 * (xv[2] + xv[4])
                + 20 * xv[3]
                + (-0.0050225266 * yv[0]) + (0.0517530339 * yv[1])
                + (-0.2634693483 * yv[2]) + (0.6743275253 * yv[3])
                + (-1.3052133493 * yv[4]) + (1.1876006802 * yv[5]);
        sampleNew[n] = (float) yv[6];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 6 -a 1.2500000000e-01 0.0000000000e+00 -l */

#define GAIN_lp_bu_co_5_n_6_40sps   9.508895986e+02

void filter_lp_bu_co_5_n_6_40sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_40sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = (xv[0] + xv[6]) + 6 * (xv[1] + xv[5]) + 15 * (xv[2] + xv[4])
                + 20 * xv[3]
                     + ( -0.0433569884 * yv[0]) + (  0.3911172306 * yv[1])
                     + ( -1.5172788447 * yv[2]) + (  3.2597642798 * yv[3])
                     + ( -4.1360809983 * yv[4]) + (  2.9785299261 * yv[5]);
        sampleNew[n] = (float) yv[6];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 6 -a 1.0000000000e-01 0.0000000000e+00 -l */

#define GAIN_lp_bu_co_5_n_6_50sps   2.936532839e+03

void filter_lp_bu_co_5_n_6_50sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_50sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = (xv[0] + xv[6]) + 6 * (xv[1] + xv[5]) + 15 * (xv[2] + xv[4])
                + 20 * xv[3]
                     + ( -0.0837564796 * yv[0]) + (  0.7052741145 * yv[1])
                     + ( -2.5294949058 * yv[2]) + (  4.9654152288 * yv[3])
                     + ( -5.6586671659 * yv[4]) + (  3.5794347983 * yv[5]);
        sampleNew[n] = (float) yv[6];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 6 -a 1.0000000000e-01 0.0000000000e+00 -l */

#define GAIN_lp_bu_co_5_n_6_100sps   1.165969038e+05


void filter_lp_bu_co_5_n_6_100sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_100sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] =   (xv[0] + xv[6]) + 6 * (xv[1] + xv[5]) + 15 * (xv[2] + xv[4])
                     + 20 * xv[3]
                     + ( -0.2951724313 * yv[0]) + (  2.1290387500 * yv[1])
                     + ( -6.4411118810 * yv[2]) + ( 10.4690788930 * yv[3])
                     + ( -9.6495177287 * yv[4]) + (  4.7871354989 * yv[5]);
        sampleNew[n] = (float) yv[6];
    }

}



/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.SG2KCmdRecursionFilter lp_bu_co_5_n_6 200
  Created: Wed Mar 16 18:02:55 CET 2016
*/
#define GAIN_lp_bu_co_5_n_6_200sps  5702375.985737481
void filter_lp_bu_co_5_n_6_200sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = xv[6];
         xv[6] = (double) sample[n] / GAIN_lp_bu_co_5_n_6_200sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = xv[0] + (6.0 * xv[1]) + (15.0 * xv[2]) + (20.0 * xv[3]) + (15.0 * xv[4]) + (6.0 * xv[5]) + xv[6]
                   + (-0.5446010675601198 * yv[0]) + (3.59806353388664 * yv[1])
                   + (-9.923048570770405 * yv[2]) + (14.62378756660761 * yv[3])
                   + (-12.147425170416902 * yv[4]) + (5.393212484861355 * yv[5])
                  ;
        sampleNew[n] = (float) yv[6];
    }
}


/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_lp_bu_co_5_n_6
// _DOC_ =============================

float* filter_lp_bu_co_5_n_6(

        double deltaTime, // dt or timestep of data samples
        float* sample, // array of num_samples data samples
        int num_samples, // the number of samples in array sample

        timedomain_memory** pmem, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        int useMemory // set to TRUE_INT=1 if function is called for packets of data in sequence, FALSE_INT = 0 otherwise

        ) {

    int ierror = 0;

    // _DOC_ =============================
    // _DOC_ apply algoritm

    int memory_size = 6;

    // initialize memory object
    timedomain_memory* mem = *pmem;
    if (mem == NULL) {
        mem = init_timedomain_memory(memory_size, 0.0, memory_size, 0.0, deltaTime);
    }

    // create array for time-series results
    float* sampleNew = NULL;
    sampleNew = calloc(num_samples, sizeof (float));

    // _DOC_ filter
    int i;
    for (i = 0; i < memory_size; i++) {
        xv[i] = mem->input[i];
        yv[i] = mem->output[i];
    }
    if (fabs(deltaTime * 20 - 1.0) < 0.01) // 20 sps
        filter_lp_bu_co_5_n_6_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_lp_bu_co_5_n_6_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_lp_bu_co_5_n_6_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_lp_bu_co_5_n_6_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_lp_bu_co_5_n_6_100sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 200.0 - 1.0) < 0.01) // 200 sps
        filter_lp_bu_co_5_n_6_200sps_impl(sample, num_samples, sampleNew);
    else
        ierror = 1;
    ;

    if (useMemory) {
        for (i = 0; i < memory_size; i++) {
            mem->input[i] = xv[i];
            mem->output[i] = yv[i];
        }
    } else {
        free_timedomain_memory(&mem);
        mem = NULL;
    }
    *pmem = mem;

    if (!ierror) {
        int n;
        for (n = 0; n < num_samples; n++)
            sample[n] = sampleNew[n];
    }

    if (sampleNew != NULL)
        free(sampleNew);

    if (ierror)
        return (NULL);

    return (sample);

}




