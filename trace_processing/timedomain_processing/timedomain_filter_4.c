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


void filter_bp_bu_co_0_333_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_80sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_120sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_bp_bu_co_0_333_5_n_4_200sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================


// Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 5.0000000000e-02 2.5000000000e-01 -l */

#define GAIN_bp_bu_co_0_333_5_n_4_20sps   1.321029146e+01

void filter_bp_bu_co_0_333_5_n_4_20sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_20sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0189792683 * yv[0]) + (0.1180993002 * yv[1])
                + (-0.6934028616 * yv[2]) + (2.1704960344 * yv[3])
                + (-4.1314930101 * yv[4]) + (5.8886267932 * yv[5])
                + (-6.1697575740 * yv[6]) + (3.8362300590 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 5.0000000000e-02 2.5000000000e-01 -l */

#define GAIN_bp_bu_co_0_333_5_n_4_25sps   2.669753586e+01

void filter_bp_bu_co_0_333_5_n_4_25sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_25sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.0383103754 * yv[0]) + (0.3757335424 * yv[1])
                + (-1.8099117703 * yv[2]) + (5.0952429321 * yv[3])
                + (-9.3212973019 * yv[4]) + (11.6452819910 * yv[5])
                + (-9.6433270076 * yv[6]) + (4.6965506290 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 2.5000000000e-02 1.2500000000e-01 -l */

#define GAIN_bp_bu_co_0_333_5_n_4_40sps   1.231479138e+02

void filter_bp_bu_co_0_333_5_n_4_40sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_40sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.1397228909 * yv[0]) + (1.3523448031 * yv[1])
                + (-5.8240017500 * yv[2]) + (14.5594614240 * yv[3])
                + (-23.1342201310 * yv[4]) + (23.9335612860 * yv[5])
                + (-15.7182523890 * yv[6]) + (5.9708283850 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 2.0000000000e-02 1.0000000000e-01 -l */

#define GAIN_bp_bu_co_0_333_5_n_4_50sps   2.622084615e+02

void filter_bp_bu_co_0_333_5_n_4_50sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_50sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.2103340565 * yv[0]) + (1.9757508770 * yv[1])
                + (-8.2027432078 * yv[2]) + (19.6578196770 * yv[3])
                + (-29.7493361710 * yv[4]) + (29.1088566580 * yv[5])
                + (-17.9694208060 * yv[6]) + (6.3894067848 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_0.333_5_n_4 80
  Created: Thu Mar 17 09:15:23 CET 2016
*/
#define GAIN_bp_bu_co_0_333_5_n_4_80sps  1371.2373652446056
void filter_bp_bu_co_0_333_5_n_4_80sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_80sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.3812906976559329 * yv[0]) + (3.3968676698784934 * yv[1])
                   + (-13.290622557338684 * yv[2]) + (29.828726545975037 * yv[3])
                   + (-41.99987982003667 * yv[4]) + (37.98756304465559 * yv[5])
                   + (-21.549450514950383 * yv[6]) + (7.008086322309175 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}



/** Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 4 -a 1.0000000000e-02 5.0000000000e-02 -l */

#define GAIN_bp_bu_co_0_333_5_n_4_100sps  3.087169919e+03

void filter_bp_bu_co_0_333_5_n_4_100sps_impl(float* sample, int num_samples, float* sampleNew) {
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
        xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_100sps;
        yv[0] = yv[1];
        yv[1] = yv[2];
        yv[2] = yv[3];
        yv[3] = yv[4];
        yv[4] = yv[5];
        yv[5] = yv[6];
        yv[6] = yv[7];
        yv[7] = yv[8];
        yv[8] = (xv[0] + xv[8]) - 4 * (xv[2] + xv[6]) + 6 * xv[4]
                + (-0.4632414443 * yv[0]) + (4.0468645647 * yv[1])
                + (-15.5046457910 * yv[2]) + (34.0262165520 * yv[3])
                + (-46.7818958420 * yv[4]) + (41.2599386970 * yv[5])
                + (-22.7943242500 * yv[6]) + (7.2110875124 * yv[7]);
        sampleNew[n] = (float) yv[8];
    }

}



/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: java net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_0.333_5_n_4 120
  Created: Thu Mar 17 08:41:17 CET 2016
*/
#define GAIN_bp_bu_co_0_333_5_n_4_120sps  6054.6167476774535
void filter_bp_bu_co_0_333_5_n_4_120sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_120sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.5270683373035576 * yv[0]) + (4.542271461202572 * yv[1])
                   + (-17.154676937679593 * yv[2]) + (37.08283412457895 * yv[3])
                   + (-50.1824530431562 * yv[4]) + (43.53148361157432 * yv[5])
                   + (-23.6377102657209 * yv[6]) + (7.345319386183341 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}


/** Digital filter designed by net.alomax.timedom.mkfilter.MakeFilter
  Command line: net.alomax.timedom.mkfilter.SG2KCmdMakeFilter bp_bu_co_0_333_5_n_4 200
  Created: Wed Mar 16 18:14:51 CET 2016
*/
#define GAIN_bp_bu_co_0_333_5_n_4_200sps  41614.975008896996
void filter_bp_bu_co_0_333_5_n_4_200sps_impl(float* sample, int num_samples, float* sampleNew) {
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
         xv[8] = (double) sample[n] / GAIN_bp_bu_co_0_333_5_n_4_200sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = yv[6];
         yv[6] = yv[7];
         yv[7] = yv[8];
         yv[8] = xv[0] + (-4.0 * xv[2]) + (6.0 * xv[4]) + (-4.0 * xv[6]) + xv[8]
                   + (-0.6814485353853469 * yv[0]) + (5.70785879065745 * yv[1])
                   + (-20.92892144587282 * yv[2]) + (43.87684681373103 * yv[3])
                   + (-57.52463640015641 * yv[4]) + (48.294579984087235 * yv[5])
                   + (-25.355043788117698 * yv[6]) + (7.6107645810505 * yv[7])
                  ;
        sampleNew[n] = (float) yv[8];
    }
}


/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_bp_bu_co_0_333_5_n_4
// _DOC_ =============================

float* filter_bp_bu_co_0_333_5_n_4(

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
        filter_bp_bu_co_0_333_5_n_4_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_bp_bu_co_0_333_5_n_4_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_bp_bu_co_0_333_5_n_4_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_bp_bu_co_0_333_5_n_4_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 80.0 - 1.0) < 0.01) // 80 sps
        filter_bp_bu_co_0_333_5_n_4_80sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_bp_bu_co_0_333_5_n_4_100sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 120.0 - 1.0) < 0.01) // 120 sps
        filter_bp_bu_co_0_333_5_n_4_120sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 200.0 - 1.0) < 0.01) // 200 sps
        filter_bp_bu_co_0_333_5_n_4_200sps_impl(sample, num_samples, sampleNew);
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




