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


void filter_WWSSN_SP_DISP_20sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_25sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_40sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_50sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_80sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_100sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_120sps_impl(float* sample, int num_samples, float* sampleNew);
void filter_WWSSN_SP_DISP_200sps_impl(float* sample, int num_samples, float* sampleNew);


/** function to select filter function */
// _DOC_ =============================
// _DOC_ timedomain_filter.c
// _DOC_ =============================



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 20
*/
#define GAIN_WWSSN_SP_DISP_20sps  11.384796745007629  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_20sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_20sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
             yv[5] = 0.0
              + -1.0 * xv[0]
              + 1.0 * xv[1]
              + 2.0 * xv[2]
              + -2.0 * xv[3]
              + -1.0 * xv[4]
              + 1.0 * xv[5]
              + 0.09356489680276334 * yv[0]
              + -0.7883362068655061 * yv[1]
              + 2.5677921900946004 * yv[2]
              + -4.081577851746812 * yv[3]
              + 3.2001821719856034 * yv[4]
         ;
        sampleNew[n] = (float) yv[5];
    }
}



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 25
*/
#define GAIN_WWSSN_SP_DISP_25sps  14.480850029855901  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_25sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_25sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
             yv[5] = 0.0
              + -1.0 * xv[0]
              + 1.0 * xv[1]
              + 2.0 * xv[2]
              + -2.0 * xv[3]
              + -1.0 * xv[4]
              + 1.0 * xv[5]
              + 0.15027625180451618 * yv[0]
              + -1.1306160455905583 * yv[1]
              + 3.335079746441332 * yv[2]
              + -4.846324890355234 * yv[3]
              + 3.488183678115408 * yv[4]
         ;
        sampleNew[n] = (float) yv[5];
    }
}



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 40
*/
#define GAIN_WWSSN_SP_DISP_40sps  26.87498693336043  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_40sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_40sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
             yv[5] = 0.0
              + -1.0 * xv[0]
              + 1.0 * xv[1]
              + 2.0 * xv[2]
              + -2.0 * xv[3]
              + -1.0 * xv[4]
              + 1.0 * xv[5]
              + 0.30588379624093087 * yv[0]
              + -1.959348734774371 * yv[1]
              + 4.986716025966793 * yv[2]
              + -6.311552280946859 * yv[3]
              + 3.9778495112157426 * yv[4]
         ;
        sampleNew[n] = (float) yv[5];
    }
}



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 50
*/
#define GAIN_WWSSN_SP_DISP_50sps  37.592510151554215  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_50sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_50sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
             yv[5] = 0.0
              + -1.0 * xv[0]
              + 1.0 * xv[1]
              + 2.0 * xv[2]
              + -2.0 * xv[3]
              + -1.0 * xv[4]
              + 1.0 * xv[5]
              + 0.3876548101139932 * yv[0]
              + -2.3586431139197312 * yv[1]
              + 5.717175664538988 * yv[2]
              + -6.9056076696279005 * yv[3]
              + 4.159254622810583 * yv[4]
         ;
        sampleNew[n] = (float) yv[5];
    }
}


/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter /Users/anthony/work/mseed_processing/td-filters/WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 80
  Created: Thu Mar 17 09:16:17 CET 2016
*/
#define GAIN_WWSSN_SP_DISP_80sps  81.2516928933165  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_80sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_80sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = (-1.0 * xv[0]) + xv[1] + (2.0 * xv[2]) + (-2.0 * xv[3]) + (-1.0 * xv[4]) + xv[5]
                   + (0.5530676235696055 * yv[0]) + (-3.1210305281574167 * yv[1]) + (7.034838500121374 * yv[2]) + (-7.918218479866248 * yv[3]) + (4.451324125577493 * yv[4]);
        sampleNew[n] = (float) yv[5];
    }
}



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 100
*/
#define GAIN_WWSSN_SP_DISP_100sps  119.88757255809972  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_100sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_100sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
             yv[5] = 0.0
              + -1.0 * xv[0]
              + 1.0 * xv[1]
              + 2.0 * xv[2]
              + -2.0 * xv[3]
              + -1.0 * xv[4]
              + 1.0 * xv[5]
              + 0.6226193139583716 * yv[0]
              + -3.428045704228179 * yv[1]
              + 7.5430154717185225 * yv[2]
              + -8.292092156469884 * yv[3]
              + 4.554496562272316 * yv[4]
         ;
        sampleNew[n] = (float) yv[5];
    }
}


/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter /Users/anthony/work/mseed_processing/td-filters/WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 120
  Created: Thu Mar 17 08:48:58 CET 2016
*/
#define GAIN_WWSSN_SP_DISP_120sps  166.12956493180545  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_120sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_120sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = (-1.0 * xv[0]) + xv[1] + (2.0 * xv[2]) + (-2.0 * xv[3]) + (-1.0 * xv[4]) + xv[5]
                   + (0.6737810988340343 * yv[0]) + (-3.6497741986579535 * yv[1]) + (7.903357975801929 * yv[2]) + (-8.552377148717182 * yv[3]) + (4.625009552082162 * yv[4]);
        sampleNew[n] = (float) yv[5];
    }
}



/** Digital filter designed by net.alomax.timedom.RecursionFilter
  Command line: java net.alomax.timedom.RecursionFilter /Users/anthony/work/mseed_processing/td-filters/WWSSN-SP-DISP.pz WWSSN_SP_DISP 1.2 200
  Created: Wed Mar 16 18:11:23 CET 2016
*/
#define GAIN_WWSSN_SP_DISP_200sps  427.0855747537404  // evaluated at 1.2Hz
void filter_WWSSN_SP_DISP_200sps_impl(float* sample, int num_samples, float* sampleNew) {
     int n;

     for (n = 0; n < num_samples; n++) {

         xv[0] = xv[1];
         xv[1] = xv[2];
         xv[2] = xv[3];
         xv[3] = xv[4];
         xv[4] = xv[5];
         xv[5] = (double) sample[n] / GAIN_WWSSN_SP_DISP_200sps;
         yv[0] = yv[1];
         yv[1] = yv[2];
         yv[2] = yv[3];
         yv[3] = yv[4];
         yv[4] = yv[5];
         yv[5] = (-1.0 * xv[0]) + xv[1] + (2.0 * xv[2]) + (-2.0 * xv[3]) + (-1.0 * xv[4]) + xv[5]
                   + (0.7890623004290417 * yv[0]) + (-4.138368611887681 * yv[1]) + (8.679954077251484 * yv[2]) + (-9.101031602448385 * yv[3]) + (4.770383607967926 * yv[4]);
        sampleNew[n] = (float) yv[5];
    }
}


/** function to filter  */
// _DOC_ =============================
// _DOC_ timedomain_filter_WWSSN_SP_DISP
// _DOC_ =============================

float* filter_WWSSN_SP_DISP(

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
        filter_WWSSN_SP_DISP_20sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 25.0 - 1.0) < 0.01) // 25 sps
        filter_WWSSN_SP_DISP_25sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 40.0 - 1.0) < 0.01) // 40 sps
        filter_WWSSN_SP_DISP_40sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 50.0 - 1.0) < 0.01) // 50 sps
        filter_WWSSN_SP_DISP_50sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 80 - 1.0) < 0.01) // 80 sps
        filter_WWSSN_SP_DISP_80sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 100.0 - 1.0) < 0.01) // 100 sps
        filter_WWSSN_SP_DISP_100sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 120 - 1.0) < 0.01) // 120 sps
        filter_WWSSN_SP_DISP_120sps_impl(sample, num_samples, sampleNew);
    else if (fabs(deltaTime * 200 - 1.0) < 0.01) // 200 sps
        filter_WWSSN_SP_DISP_200sps_impl(sample, num_samples, sampleNew);
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




