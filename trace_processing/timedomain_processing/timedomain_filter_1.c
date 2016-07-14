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


void filter_bp_bu_co_1_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_80sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_120sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_1_5_n_4_200sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================


// Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 5.0000000000e-02 2.5000000000e-01 -l */

#define GAIN_bp_bu_co_1_5_n_4_20sps   2.146699830e+01

void filter_bp_bu_co_1_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_20sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0301188750 * yv[0]) + (0.2202522377 * yv[1])
                + (-0.9262877390 * yv[2]) + (2.4641234718 * yv[3])
                + (-4.3888534987 * yv[4]) + (5.7182667684 * yv[5])
                + (-5.5467354882 * yv[6]) + (3.4743955343 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 5.0000000000e-02 2.5000000000e-01 -l */

#define GAIN_bp_bu_co_1_5_n_4_25sps   4.372432451e+01

void filter_bp_bu_co_1_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_25sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0632116957 * yv[0]) + (0.5497941587 * yv[1])
                + (-2.2946324563 * yv[2]) + (5.8420743415 * yv[3])
                + (-9.9273438837 * yv[4]) + (11.6602832000 * yv[5])
                + (-9.2722033772 * yv[6]) + (4.5020382411 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 2.5000000000e-02 1.2500000000e-01 -l */

#define GAIN_bp_bu_co_1_5_n_4_40sps   2.072730615e+02

void filter_bp_bu_co_1_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_40sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.1873794924 * yv[0]) + (1.6902666989 * yv[1])
                + (-6.8495106788 * yv[2]) + (16.2719675000 * yv[3])
                + (-24.7894815200 * yv[4]) + (24.8026467890 * yv[5])
                + (-15.9059594260 * yv[6]) + (5.9673400561 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 2.0000000000e-02 1.0000000000e-01 -l */

#define GAIN_bp_bu_co_1_5_n_4_50sps   4.474249436e+02

void filter_bp_bu_co_1_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_50sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.2644548164 * yv[0]) + (2.3626908016 * yv[1])
                + (-9.3876754462 * yv[2]) + (21.6613786360 * yv[3])
                + (-31.7461576020 * yv[4]) + (30.2569374390 * yv[5])
                + (-18.3072090400 * yv[6]) + (6.4244688269 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_1_5_n_4 80
  Created: Thu Mar 17 09:13:11 CET 2016
*/
#define GAIN_bp_bu_co_1_5_n_4_80sps  2400.2285970558532
void filter_bp_bu_co_1_5_n_4_80sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_80sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.43826514226197977 * yv[0]) + (3.805268874824978 * yv[1])
                   + (-14.54460702035583 * yv[2]) + (31.963311873423137 * yv[3])
                   + (-44.170549839423266 * yv[4]) + (39.30212514428771 * yv[5])
                   + (-21.986555899199306 * yv[6]) + (7.069271395839914 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 1.0000000000e-02 5.0000000000e-02 -l */

#define GAIN_bp_bu_co_1_5_n_4_100sps  5.457655893e+03

void filter_bp_bu_co_1_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_100sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] =   (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                     + ( -0.5174781998 * yv[0]) + (  4.4348861251 * yv[1])
                     + (-16.6938514710 * yv[2]) + ( 36.0485071710 * yv[3])
                     + (-48.8405001620 * yv[4]) + ( 42.5126474810 * yv[5])
                     + (-23.2155317960 * yv[6]) + (  7.2713207414 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_1_5_n_4 120
  Created: Thu Mar 17 08:38:06 CET 2016
*/
#define GAIN_bp_bu_co_1_5_n_4_120sps  10779.197468842613
void filter_bp_bu_co_1_5_n_4_120sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_120sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.5778338981055613 * yv[0]) + (4.904704853138681 * yv[1])
                   + (-18.263246558410763 * yv[2]) + (38.9650496249689 * yv[3])
                   + (-52.09721037749206 * yv[4]) + (44.69771110851589 * yv[5])
                   + (-24.031147622129552 * yv[6]) + (7.401972842403145 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_1_5_n_4 200
  Created: Wed Mar 16 17:56:38 CET 2016
*/
#define GAIN_bp_bu_co_1_5_n_4_200sps  75217.81988652876
void filter_bp_bu_co_1_5_n_4_200sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_1_5_n_4_200sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.7199103272918717 * yv[0]) + (5.980801257897405 * yv[1])
                   + (-21.758903234813683 * yv[2]) + (45.27863033385741 * yv[3])
                   + (-58.94461344276549 * yv[4]) + (49.15718723254562 * yv[5])
                   + (-25.64597452239443 * yv[6]) + (7.652782702460034 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}


/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_bp_bu_co_1_5_n_4
// _DOC_ =============================

float* filter_bp_bu_co_1_5_n_4(

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
        filter_bp_bu_co_1_5_n_4_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_bp_bu_co_1_5_n_4_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_bp_bu_co_1_5_n_4_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_bp_bu_co_1_5_n_4_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 80.0 - 1.0) < 0.01) // 80 sps
        filter_bp_bu_co_1_5_n_4_80sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_bp_bu_co_1_5_n_4_100sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 120.0 - 1.0) < 0.01) // 120 sps
        filter_bp_bu_co_1_5_n_4_120sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 200.0 - 1.0) < 0.01) // 200 sps
        filter_bp_bu_co_1_5_n_4_200sps_impl(sample, num_samples, sampleNew);
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




