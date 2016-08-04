/*
 * This file is part of the Anthony Lomax C Library.
 *
 * Copyright (C) 2006-2009 Anthony Lomax <anthony@alomax.net www.alomax.net>
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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "timedomain_memory.h"
#include "timedomain_processes.h"

#define DEBUG 0

/**  Method to copy data array to newly allocated array */

float* copyData(

        float* sample, // array of num_samples data samples
        int num_samples // the number of samples in array sample

        ) {
    float* new_sample = calloc(num_samples, sizeof (float));
    int n;
    for (n = 0; n < num_samples; n++)
        new_sample[n] = sample[n];

    return (new_sample);
}

/**  Method to calculate square of data values */

float* applySqr(

        float* sample, // array of num_samples data samples
        int num_samples // the number of samples in array sample

        ) {
    int i;

    for (i = 0; i < num_samples; i++) {
        sample[i] = sample[i] * sample[i];
    }

    return (sample);

}

/**  Method to calculate absolute value of data values */

float* applyAbs(

        float* sample, // array of num_samples data samples
        int num_samples // the number of samples in array sample

        ) {
    int i;

    for (i = 0; i < num_samples; i++) {
        sample[i] = fabs(sample[i]);
    }

    return (sample);

}


// 20101223 AJL - Added persistence
static int num_sampleNewArray_alloc = 0;
static float* sampleNewArray = NULL;
static int num_derivArray_alloc = 0;
static double* derivArray = NULL;

/** function to apply triangle smoothing  */

float* applyTriangleSmoothing(

        double deltaTime, // dt or timestep of data samples
        float* sample, // array of num_samples data samples
        int num_samples, // the number of samples in array sample
        double windowHalfWidthSec, // smoothing window half width in seconds

        timedomain_memory** pmem, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        int useMemory // set to TRUE_INT=1 if function is called for packets of data in sequence, FALSE_INT = 0 otherwise

        ) {
    // smoothing window half width in samples
    //int windowHalfWidth = 1 + (int) (windowHalfWidthSec / deltaTime);
    int windowHalfWidth = 1 + (int) (0.5 + windowHalfWidthSec / deltaTime); // 20111004 AJL - avoid integer instabilty due to slight change in deltaTime values

    // initialize memory object
    timedomain_memory* mem = *pmem;
    if (useMemory && mem == NULL) {
        mem = init_timedomain_memory(2 * windowHalfWidth, 0.0, 2 * windowHalfWidth, 0.0, deltaTime);
    }

    // initialize array for time-series results
    if (sampleNewArray == NULL || num_samples > num_sampleNewArray_alloc) {
        sampleNewArray = realloc(sampleNewArray, num_samples * sizeof (float));
        num_sampleNewArray_alloc = num_samples;
    }


    int i1, i2, iwt;
    double sum;
    double weight;

    double *wt = calloc(2 * windowHalfWidth + 1, sizeof (double));
    int n;
    for (n = 0; n < windowHalfWidth; n++) {
        wt[n] = wt[2 * windowHalfWidth - n] = 1.0 - ((double) (windowHalfWidth - n) / (double) windowHalfWidth);
    }
    wt[windowHalfWidth] = 1.0;

    int i;
    for (i = 0; i < num_samples; i++) {

        if (!useMemory) {
            i1 = i - windowHalfWidth;
            if (i1 < 0) // truncated
                i1 = 0;
            i2 = i + windowHalfWidth; // a-causal triangle of width 2 * windowHalfWidth
            if (i2 > num_samples)
                i2 = num_samples;
            iwt = i1 - i + windowHalfWidth; // 0 if not truncated
        } else {
            i1 = i - 2 * windowHalfWidth;
            i2 = i; // causal triangle of width 2 * windowHalfWidth
            iwt = i1 - i + 2 * windowHalfWidth; // 0 if not truncated
        }

        sum = 0.0;
        weight = 0.0;
        float value = 0.0f;
        for (n = i1; n < i2; n++) {
            if (useMemory && n < 0) {
                value = mem->input[2 * windowHalfWidth + n];
            } else {
                value = sample[n];
            }
            sum += value * wt[iwt];
            weight += wt[iwt];
            iwt++;
        }
        if (weight > 0.0)
            sampleNewArray[i] = (float) (sum / weight);
        else
            sampleNewArray[i] = 0.0f;

    }

    if (useMemory) {
        updateInput(mem, sample, num_samples);
    }
    *pmem = mem;

    for (n = 0; n < num_samples; n++)
        sample[n] = sampleNewArray[n];

    free(wt);

    return (sample);

}

/** function to apply boxcar smoothing  */

float* applyBoxcarSmoothing(

        double deltaTime, // dt or timestep of data samples
        float* sample, // array of num_samples data samples
        int num_samples, // the number of samples in array sample
        double windowHalfWidthSec, // smoothing window half width in seconds

        timedomain_memory** pmem, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        int useMemory // set to TRUE_INT=1 if function is called for packets of data in sequence, FALSE_INT = 0 otherwise

        ) {

    // smoothing window half width in samples
    //int windowHalfWidth = 1 + (int) (windowHalfWidthSec / deltaTime); // 20111004 AJL
    int windowHalfWidth = 1 + (int) (0.5 + windowHalfWidthSec / deltaTime); // 20111004 AJL - avoid integer instabilty due to slight change in deltaTime values

    // initialize memory object
    timedomain_memory* mem = *pmem;
    if (useMemory) {
        if (mem == NULL) {
            mem = init_timedomain_memory(2 * windowHalfWidth, 0.0, 2 * windowHalfWidth, 0.0, deltaTime);
        } else if (mem->numInput != 2 * windowHalfWidth) {
            fprintf(stderr, "\nError: timedomain_processes.applyBoxcarSmoothing(): Existing timedomain_memory size = %d not equal to 2 * windowHalfWidth = %d\n", mem->numInput, 2 * windowHalfWidth);
        }
    }

    // initialize array for time-series results
    if (sampleNewArray == NULL || num_samples > num_sampleNewArray_alloc) {
        sampleNewArray = realloc(sampleNewArray, num_samples * sizeof (float));
        num_sampleNewArray_alloc = num_samples;
    }


    int i1, i2, i;

    i = 0;
    if (!useMemory) {
        i1 = i - windowHalfWidth;
        i2 = i + windowHalfWidth; // a-causal boxcar of width 2 * windowHalfWidth
    } else {
        i1 = i - 2 * windowHalfWidth;
        i2 = i; // causal boxcar of width 2 * windowHalfWidth
    }

    double sum = 0.0;
    int icount = 0;

    for (i = 0; i < num_samples; i++) {

        if (!useMemory) {
            if (i1 < 0) // truncated
                i1 = 0;
            if (i2 > num_samples)
                i2 = num_samples;
        }

        float value = 0.0f;
        if (icount == 0) { // first pass, accumulate sum
            int n;
            // 20110919 AJL  for (n = i1; n < i2; n++) {
            for (n = i1; n < i2 + 1; n++) { // 20110919 AJL - bug fix, should inlcude sample i2 in initial sum
                if (useMemory && n < 0) {
                    value = mem->input[2 * windowHalfWidth + n];
                } else {
                    value = sample[n];
                }
                sum += value;
                icount++;
            }
        } else { // later passes, update sum
            if (useMemory && (i1 - 1) < 0) {
                value = mem->input[2 * windowHalfWidth + (i1 - 1)];
            } else {
                value = sample[(i1 - 1)];
            }
            sum -= value;
            if (useMemory && i2 < 0) {
                value = mem->input[2 * windowHalfWidth + i2];
            } else {
                value = sample[i2];
            }
            sum += value;
        }
        if (icount > 0)
            sampleNewArray[i] = (float) (sum / (double) icount);
        else
            sampleNewArray[i] = 0.0f;

        i1++;
        i2++;

    }

    if (useMemory) {
        updateInput(mem, sample, num_samples);
    }
    *pmem = mem;

    int n;
    for (n = 0; n < num_samples; n++)
        sample[n] = sampleNewArray[n];

    return (sample);

}

/*** function to calculate instantaneous period in a fixed window (Tau_c)
 *
 * implements eqs 1-3 in :
 *  Allen, R.M., and H. Kanamori,
 *     The Potential for Earthquake Early Warning in Southern California,
 *     Science, 300 (5620), 786-789, 2003.
 *
 * except use a fixed width window instead of decay function
 */

#define MIN_FLOAT_VAL (1.0e-20)

float* applyInstantPeriodWindowed(

        double deltaTime, // dt or timestep of data samples
        float* sample, // array of num_samples data samples
        int num_samples, // the number of samples in array sample
        double windowWidth, // window width in seconds

        timedomain_memory** pmem, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        timedomain_memory** pmem_dval, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        int useMemory // set to TRUE_INT=1 if function is called for packets of data in sequence, FALSE_INT = 0 otherwise

        ) {

    //int windowLength = 1 + (int) (windowWidth / deltaTime);
    int windowLength = 1 + (int) (0.5 + windowWidth / deltaTime); // 20111004 AJL - avoid integer instabilty due to slight change in deltaTime values
    if (windowLength < 1) {
        windowLength = 1;
    }

    // initialize memory object
    timedomain_memory* memory = *pmem;
    if (useMemory && memory == NULL) {
        memory = init_timedomain_memory(windowLength, 0.0f, 1, 0.0f, deltaTime);
    }
    timedomain_memory* memory_dval = *pmem_dval;
    if (useMemory && memory_dval == NULL) {
        memory_dval = init_timedomain_memory(windowLength, 0.0f, 1, 0.0f, deltaTime);
    }

    // initialize array for time-series results
    if (sampleNewArray == NULL || num_samples > num_sampleNewArray_alloc) {
        sampleNewArray = realloc(sampleNewArray, num_samples * sizeof (float));
        num_sampleNewArray_alloc = num_samples;
    }

    // initialize array for time-series results
    if (derivArray == NULL || num_samples > num_derivArray_alloc) {
        derivArray = realloc(derivArray, num_samples * sizeof (double));
        num_derivArray_alloc = num_samples;
    }

    // use double precision where possible to avoid loss of precision in accumulating xval and dval
    double twopi = 2.0 * PI;
    double sampleLast = 0.0;
    double sample_d = 0.0;
    double deriv_d = 0.0;
    double xval = 0.0;
    double dval = 0.0;
    if (useMemory) { // using stored memory
        sampleLast = (double) memory->input[windowLength - 1];
        xval = (double) memory->output[0];
        dval = (double) memory_dval->output[0];
    }

    int i;
    for (i = 0; i < num_samples; i++) {

        sample_d = (double) sample[i];
        deriv_d = (sample_d - sampleLast) / deltaTime;
        int indexBegin = i - windowLength;
        if (indexBegin >= 0) {
            xval = xval - ((double) sample[indexBegin]) * ((double) sample[indexBegin]) + sample_d * sample_d;
            dval = dval - derivArray[indexBegin] * derivArray[indexBegin] + deriv_d * deriv_d;
        } else {
            int index = i;
            xval = xval - memory->input[index] * memory->input[index] + sample_d * sample_d;
            dval = dval - memory_dval->input[index] * memory_dval->input[index] + deriv_d * deriv_d;
        }
        derivArray[i] = deriv_d;
        sampleLast = sample_d;
        //if (xval > MIN_FLOAT_VAL && dval > MIN_FLOAT_VAL) {
        if (dval > MIN_FLOAT_VAL) {
            sampleNewArray[i] = (float) (twopi * sqrt(xval / dval));
        } else {
            sampleNewArray[i] = 0.0f;
        }

    }

    if (useMemory) {
        memory->output[0] = xval;
        updateInput(memory, sample, num_samples);
        memory_dval->output[0] = dval;
        updateInput_d(memory_dval, derivArray, num_samples);
    }
    *pmem = memory;
    *pmem_dval = memory_dval;

    int n;
    for (n = 0; n < num_samples; n++)
        sample[n] = sampleNewArray[n];

    return (sample);

}

/** function to apply integation  */

float* applyIntegral(

        double deltaTime, // dt or timestep of data samples
        float* sample, // array of num_samples data samples
        int num_samples, // the number of samples in array sample

        timedomain_memory** pmem, // memory structure/object
        // _DOC_ pointer to a memory structure/object is used so that this function can be called repetedly for packets of data in sequence from the same channel.
        // The calling routine is responsible for managing and associating the correct mem structures/objects with each channel.  On first call to this function for each channel set mem = NULL.
        int useMemory // set to TRUE_INT=1 if function is called for packets of data in sequence, FALSE_INT = 0 otherwise

        ) {

    // initialize memory object
    timedomain_memory* mem = *pmem;
    if (useMemory && mem == NULL) {
        mem = init_timedomain_memory(0, 0.0, 1, 0.0, deltaTime);
    }

    double sum = 0.0;
    if (useMemory) {
        sum = mem->output[0];
    }

    int i;
    for (i = 0; i < num_samples; i++) {
        sum += (double) sample[i] * deltaTime;
        sample[i] = (float) sum;
    }

    if (useMemory) {
        mem->output[0] = sum;
    }
    *pmem = mem;

    return (sample);

}

