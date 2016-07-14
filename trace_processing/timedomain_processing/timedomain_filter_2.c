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


void filter_hp_bu_co_0_005_n_2_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_40sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_80sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_100sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_120sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_hp_bu_co_0_005_n_2_200sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================


// Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Hp -o 2 -a 2.5000000000e-04 0.0000000000e+00 -l */

#define GAIN_hp_bu_co_0_005_n_2_20sps   1.001111338e+00

void filter_hp_bu_co_0_005_n_2_20sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_20sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = (xv[0] + xv[2]) - 2 * xv[1]
                + (-0.9977810241 * yv[0]) + (1.9977785594 * yv[1]);
        sampleNew[n] = (float) yv[2];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Hp -o 2 -a 2.0000000000e-04 0.0000000000e+00 -l */

#define GAIN_hp_bu_co_0_005_n_2_25sps   1.000888972e+00

void filter_hp_bu_co_0_005_n_2_25sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_25sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = (xv[0] + xv[2]) - 2 * xv[1]
                + (-0.9982244250 * yv[0]) + (1.9982228473 * yv[1]);
        sampleNew[n] = (float) yv[2];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Hp -o 2 -a 1.2500000000e-04 0.0000000000e+00 -l */

#define GAIN_hp_bu_co_0_005_n_2_40sps   1.000555515e+00

void filter_hp_bu_co_0_005_n_2_40sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_40sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = (xv[0] + xv[2]) - 2 * xv[1]
                + (-0.9988898959 * yv[0]) + (1.9988892794 * yv[1]);
        sampleNew[n] = (float) yv[2];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Hp -o 2 -a 1.0000000000e-04 0.0000000000e+00 -l */

#define GAIN_hp_bu_co_0_005_n_2_50sps   1.000444387e+00

void filter_hp_bu_co_0_005_n_2_50sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_50sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = (xv[0] + xv[2]) - 2 * xv[1]
                + (-0.9991118181 * yv[0]) + (1.9991114235 * yv[1]);
        sampleNew[n] = (float) yv[2];
    }

}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter hp_bu_co_0_005_n_2 80
  Created: Thu Mar 17 09:14:31 CET 2016
*/
#define GAIN_hp_bu_co_0_005_n_2_80sps  1.0
void filter_hp_bu_co_0_005_n_2_80sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_80sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = xv[0] + (-2.0 * xv[1]) + xv[2]
                   + (-1.0 * yv[0]) + (2.0 * yv[1])
                  ;
        sampleNew[n] = (float) yv[2];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Hp -o 2 -a 5.0000000000e-05 0.0000000000e+00 -l */

#define GAIN_hp_bu_co_0_005_n_2_100sps   1.000222169e+00

void filter_hp_bu_co_0_005_n_2_100sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_100sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] =   (xv[0] + xv[2]) - 2 * xv[1]
                     + ( -0.9995558104 * yv[0]) + (  1.9995557117 * yv[1]);
        sampleNew[n] = (float) yv[2];
    }

}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter hp_bu_co_0_005_n_2 120
  Created: Thu Mar 17 08:39:53 CET 2016
*/
#define GAIN_hp_bu_co_0_005_n_2_120sps  1.0
void filter_hp_bu_co_0_005_n_2_120sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_120sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = xv[0] + (-2.0 * xv[1]) + xv[2]
                   + (-1.0 * yv[0]) + (2.0 * yv[1])
                  ;
        sampleNew[n] = (float) yv[2];
    }
}



/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter hp_bu_co_0_005_n_2 200
  Created: Wed Mar 16 17:59:24 CET 2016
*/
#define GAIN_hp_bu_co_0_005_n_2_200sps  1.0
void filter_hp_bu_co_0_005_n_2_200sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = (double) sample[n] / GAIN_hp_bu_co_0_005_n_2_200sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = xv[0] + (-2.0 * xv[1]) + xv[2]
                   + (-1.0 * yv[0]) + (2.0 * yv[1])
                  ;
        sampleNew[n] = (float) yv[2];
    }
}


/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_hp_bu_co_0_005_n_2
// _DOC_ =============================

float* filter_hp_bu_co_0_005_n_2(

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

    int memory_size = 3;

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
        filter_hp_bu_co_0_005_n_2_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_hp_bu_co_0_005_n_2_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_hp_bu_co_0_005_n_2_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_hp_bu_co_0_005_n_2_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 80.0 - 1.0) < 0.01) // 80 sps
        filter_hp_bu_co_0_005_n_2_80sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_hp_bu_co_0_005_n_2_100sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 120.0 - 1.0) < 0.01) // 120 sps
        filter_hp_bu_co_0_005_n_2_120sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 200.0 - 1.0) < 0.01) // 200 sps
        filter_hp_bu_co_0_005_n_2_200sps_impl(sample, num_samples, sampleNew);
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




