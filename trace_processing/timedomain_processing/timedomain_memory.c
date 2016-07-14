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

#include "timedomain_memory.h"



/** time domain memory class ***/
// _DOC_ =============================
// _DOC_ timedomain_memory object/structure
// _DOC_ =============================

/** create timedomain_memory object/structure and initialize memory */

timedomain_memory* init_timedomain_memory(
        const int nInput,
        const double inputInitValue,
        const int nOuput,
        const double outputInitValue,
        const double deltaTime
        ) {
    int i;

    timedomain_memory* td_memory = calloc(1, sizeof (timedomain_memory));

    td_memory->deltaTime = deltaTime;

    td_memory->numInput = nInput;
    td_memory->input = calloc(td_memory->numInput, sizeof (double));
    for (i = 0; i < td_memory->numInput; i++)
        td_memory->input[i] = inputInitValue;

    td_memory->numOutput = nOuput;
    td_memory->output = calloc(td_memory->numOutput, sizeof (double));
    for (i = 0; i < td_memory->numOutput; i++)
        td_memory->output[i] = outputInitValue;

    return (td_memory);

}

/** clean up memory */

void free_timedomain_memory(timedomain_memory** ptd_memory) {

    if ((*ptd_memory) == NULL)
        return;

    free((*ptd_memory)->input);
    free((*ptd_memory)->output);

    free(*ptd_memory);
    *ptd_memory = NULL;

}

/** update memory using specifed number of points from end of specified float array */

void update(double* array, int arrayLength, float* sample, int sampleLength) {

    int n;

    if (sampleLength >= arrayLength) {
        int sampleIndex = sampleLength - arrayLength;
        for (n = 0; n < arrayLength; n++) {
            array[n] = sample[sampleIndex];
            sampleIndex++;
        }
    } else { // sample length less than memory length
        // shift data in memory
        for (n = 0; n < arrayLength - sampleLength; n++) {
            array[n] = array[n + sampleLength];
        }
        // append sample data
        int sampleIndex = 0;
        for (n = arrayLength - sampleLength; n < arrayLength; n++) {
            array[n] = sample[sampleIndex];
            sampleIndex++;
        }
    }

}

/** update output memory using specifed number of points from end of specified float array */

void updateOutput(timedomain_memory* td_memory, float* sample, int sampleLength) {

    update(td_memory->output, td_memory->numOutput, sample, sampleLength);

}

/** update input memory using specifed number of points from end of specified float array */

void updateInput(timedomain_memory* td_memory, float* sample, int sampleLength) {

    update(td_memory->input, td_memory->numInput, sample, sampleLength);

}


/** update memory using specifed number of points from end of specified double array */

void update_d(double* array, int arrayLength, double* sample, int sampleLength) {

    int n;

    if (sampleLength >= arrayLength) {
        int sampleIndex = sampleLength - arrayLength;
        for (n = 0; n < arrayLength; n++) {
            array[n] = sample[sampleIndex];
            sampleIndex++;
        }
    } else { // sample length less than memory length
        // shift data in memory
        for (n = 0; n < arrayLength - sampleLength; n++) {
            array[n] = array[n + sampleLength];
        }
        // append sample data
        int sampleIndex = 0;
        for (n = arrayLength - sampleLength; n < arrayLength; n++) {
            array[n] = sample[sampleIndex];
            sampleIndex++;
        }
    }

}


/** update output memory using specifed number of points from end of specified double array */

void updateOutput_d(timedomain_memory* td_memory, double* sample, int sampleLength) {

    update_d(td_memory->output, td_memory->numOutput, sample, sampleLength);

}

/** update input memory using specifed number of points from end of specified double array */

void updateInput_d(timedomain_memory* td_memory, double* sample, int sampleLength) {

    update_d(td_memory->input, td_memory->numInput, sample, sampleLength);

}




