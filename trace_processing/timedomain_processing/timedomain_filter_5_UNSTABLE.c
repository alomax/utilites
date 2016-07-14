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

// WARNING: some of the follwoing may be UNSTABLE:
void filter_bp_bu_co_0_0333_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_0333_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_0333_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew); // !! UNSTABLE
void filter_bp_bu_co_0_0333_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_0333_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================


// Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 1.6665000000e-03 2.5000000000e-01 -l */

#define GAIN_bp_bu_co_0_0333_5_n_4_20sps   1.086231289e+01

void filter_bp_bu_co_0_0333_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = xv[7];
        xv[7] = xv[8];
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_0333_5_n_4_20sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0176778884 * yv[0]) + (0.0753779132 * yv[1])
                + (-0.6015347487 * yv[2]) + (2.0306003928 * yv[3])
                + (-3.9558512236 * yv[4]) + (5.9393164442 * yv[5])
                + (-6.4543248722 * yv[6]) + (3.9840939646 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 1.3332000000e-03 2.0000000000e-01 -l */

#define GAIN_bp_bu_co_0_0333_5_n_4_25sps   2.191021426e+01

void filter_bp_bu_co_0_0333_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = xv[7];
        xv[7] = xv[8];
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_0333_5_n_4_25sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0308310114 * yv[0]) + (0.3100394174 * yv[1])
                + (-1.6125052582 * yv[2]) + (4.7562834793 * yv[3])
                + (-9.0057898810 * yv[4]) + (11.6014060930 * yv[5])
                + (-9.7924684639 * yv[6]) + (4.7738656215 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 8.3325000000e-04 1.2500000000e-01 -l */

#define GAIN_bp_bu_co_0_0333_5_n_4_40sps   1.000186372e+02

void filter_bp_bu_co_0_0333_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = xv[7];
        xv[7] = xv[8];
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_0333_5_n_4_40sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.1222068654 * yv[0]) + (1.2200553791 * yv[1])
                + (-5.4016427918 * yv[2]) + (13.8211691950 * yv[3])
                + (-22.3899045720 * yv[4]) + (23.5246248480 * yv[5])
                + (-15.6208829560 * yv[6]) + (5.9687877634 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 6.6660000000e-04 1.0000000000e-01 -l */

#define GAIN_bp_bu_co_0_0333_5_n_4_50sps   2.117560827e+02

void filter_bp_bu_co_0_0333_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = xv[7];
        xv[7] = xv[8];
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_0333_5_n_4_50sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.1895670169 * yv[0]) + (1.8207519209 * yv[1])
                + (-7.7102418253 * yv[2]) + (18.7972047010 * yv[3])
                + (-28.8657551320 * yv[4]) + (28.5860794080 * yv[5])
                + (-17.8100454900 * yv[6]) + (6.3715734343 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 3.3330000000e-04 5.0000000000e-02 -l */

#define GAIN_bp_bu_co_0_0333_5_n_4_100sps  2.455632197e+03

void filter_bp_bu_co_0_0333_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew) {
    int n;

    for (n = 0; n < num_samples; n++) {

        xv[0] = xv[1];
        xv[1] = xv[2];
        xv[2] = xv[3];
        xv[3] = xv[4];
        xv[4] = xv[5];
        xv[5] = xv[6];
        xv[6] = xv[7];
        xv[7] = xv[8];
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_0333_5_n_4_100sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.4407052408 * yv[0]) + (3.8830497607 * yv[1])
                + (-14.9951533350 * yv[2]) + (33.1479542400 * yv[3])
                + (-45.8765590390 * yv[4]) + (40.7025485750 * yv[5])
                + (-22.6048321970 * yv[6]) + (7.1836972357 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_bp_bu_co_0_0333_5_n_4
// _DOC_ =============================

float* filter_bp_bu_co_0_0333_5_n_4(

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

    int memory_size = 9;

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
        filter_bp_bu_co_0_0333_5_n_4_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_bp_bu_co_0_0333_5_n_4_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_bp_bu_co_0_0333_5_n_4_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_bp_bu_co_0_0333_5_n_4_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_bp_bu_co_0_0333_5_n_4_100sps_impl(sample, num_samples, sampleNew);
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




