/***************************************************************************
 * timedomain_processing_data.c:
 *
 * TODO: add doc
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2009.02.03
 ***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>

#include "timedomain_processing_data.h"
#include "ttimes.h"
//#include "location.h"

/** timedomain-processing memory class */

/** allocate and initialize a TimedomainProcessingData object
 */

TimedomainProcessingData* init_TimedomainProcessingData(double deltaTime, int flag_do_mwpd, int waveform_export_enable) {

    TimedomainProcessingData* deData = calloc(1, sizeof (TimedomainProcessingData));

    deData->deltaTime = deltaTime;

    deData->pickData = calloc(1, sizeof (PickData));
    strcpy(deData->pick_error_type, "X");
    deData->pick_error = 0.0;
    deData->pick_stream = STREAM_NULL;
    deData->use_for_location = 0;
    deData->can_use_as_location_P = 1;
    deData->use_for_location_twin_data = NULL;
    deData->merged = 0;
    deData->virtual_pick_index = -1;
    deData->sn_pick = -1.0;
    deData->a_ref = -1.0;
    deData->t50 = 0.0;
    deData->flag_complete_t50 = 0; // = 0, not complete; >= 1 completed
    deData->flag_a_ref_not_ok = 0;
    deData->flag_snr_hf_too_low = 0;
    deData->flag_clipped = 0;
    deData->flag_non_contiguous = 0;
    deData->amp_max_in_aref_window = -1.0;
    deData->amp_min_in_aref_window = -1.0;
    // brb
    deData->sn_brb_pick = -1.0;
    deData->sn_brb_signal = -1.0;
    deData->flag_snr_brb_too_low = 0;
    deData->flag_complete_snr_brb = 0;
    deData->sn_brb_int_pick = -1.0;
    deData->sn_brb_int_signal = -1.0;
    deData->flag_snr_brb_int_too_low = 0;
    deData->flag_complete_snr_brb_int = 0;
    // brb bp
    deData->sn_brb_bp_pick = -1.0;
    deData->sn_brb_bp_signal = -1.0;
    deData->flag_snr_brb_bp_too_low = 0;
    deData->flag_complete_snr_brb_bp = 0;
    // BRB HP ground velocity data
    deData->grd_mot = (GrdVel*) calloc(1, sizeof (GrdVel));
    deData->grd_mot->amp_at_pick = 0.0;
    deData->grd_mot->n_amp_max = 2 + (int) (1.0 + (START_ANALYSIS_BEFORE_P_BRB_HP + MAX_GRD_VEL_DUR) / deltaTime);
    deData->grd_mot->n_amp = 0;
    deData->grd_mot->vel = (double*) calloc(deData->grd_mot->n_amp_max, sizeof (double));
    deData->grd_mot->disp = (double*) calloc(deData->grd_mot->n_amp_max, sizeof (double));
    deData->grd_mot->have_used_memory = 0;
    deData->grd_mot->ioffset_pick = -1;
    deData->flag_complete_grd_mot = 0;
    // tauc
    deData->tauc_peak = TAUC_INVALID;
    // mb
    deData->mb = (mBData*) calloc(1, sizeof (mBData));
    deData->mb->amp_at_pick = 0.0;
    deData->mb->int_sum = 0.0;
    deData->mb->polarity = 0;
    deData->mb->peak_index_count = 0;
    deData->mb->n_amplitude_max = 2 + (int) (1.0 + (START_ANALYSIS_BEFORE_P_MB + MAX_MB_DUR) / deltaTime); // 20110331 AJL
    deData->mb->n_amplitude = 0;
    deData->mb->amplitude = (double*) calloc(deData->mb->n_amplitude_max, sizeof (double));
    deData->mb->amplitude[0] = 0.0;
    deData->mb->peak_dur = (double*) calloc(deData->mb->n_amplitude_max, sizeof (double));
    deData->mb->peak_dur[0] = -1;
    deData->mb->mag = MB_INVALID;
    deData->mb->have_used_memory = 0;
    deData->flag_complete_mb = 0;
    // mwp
    deData->mwp = (MwpData*) calloc(1, sizeof (MwpData));
    deData->mwp->amp_at_pick = 0.0;
    deData->mwp->int_sum = 0.0;
    deData->mwp->int_int_sum = 0.0;
    deData->mwp->abs_int_int_sum = 0.0;
    deData->mwp->int_int_sum_mwp_peak = 0.0;
    deData->mwp->polarity = 0;
    deData->mwp->peak_index_count = 0;
    deData->mwp->n_peak_max = 2 + (int) (1.0 + (START_ANALYSIS_BEFORE_P_BRB_HP + MAX_MWP_DUR) / deltaTime);
    deData->mwp->n_peak = 0;
    deData->mwp->peak = (double*) calloc(deData->mwp->n_peak_max, sizeof (double));
    deData->mwp->peak_dur = (double*) calloc(deData->mwp->n_peak_max, sizeof (double));
    deData->mwp->int_int = (double*) calloc(deData->mwp->n_peak_max, sizeof (double));
    deData->mwp->int_int_abs = (double*) calloc(deData->mwp->n_peak_max, sizeof (double));
    //deData->mwp->pulse_first = 0.0;
    //deData->mwp->pulse_first_dur = -1.0;
    deData->mwp->mag = MWP_INVALID;
    deData->mwp->index_mag = 0;
    deData->mwp->have_used_memory = 0;
    deData->flag_complete_mwp = 0;
    // t0
    deData->t0 = (T0Data*) calloc(1, sizeof (T0Data));
    deData->t0->amp_noise = -1.0;
    deData->t0->amp_peak = -1.0;
    deData->t0->index_peak = -1;
    deData->t0->amp_90 = -1.0;
    deData->t0->index_90 = -1;
    deData->t0->amp_80 = -1.0;
    deData->t0->index_80 = -1;
    deData->t0->amp_50 = -1.0;
    deData->t0->index_50 = -1;
    deData->t0->amp_20 = -1.0;
    deData->t0->index_20 = -1;
    deData->t0->duration_raw = T0_INVALID;
    deData->t0->duration_corr = T0_INVALID;
    deData->t0->duration_plot = T0_INVALID;
    deData->t0->have_used_memory = 0;
    deData->flag_complete_t0 = 0;
    // mwpd
    deData->mwpd = NULL;
    if (flag_do_mwpd) {
        deData->mwpd = (MwpdData*) calloc(1, sizeof (MwpdData));
        deData->mwpd->amp_at_pick = 0.0;
        deData->mwpd->int_sum = 0.0;
        // 20110331 AJL - Bug fix, array may have been too small depending on rounding of time shift and duration values
        //int n_int_sum = (int) (1.0 + (START_ANALYSIS_BEFORE_P_BRB_HP + MAX_MWPD_DUR) / deltaTime);
        deData->mwpd->n_int_sum_max = 2 + (int) (1.0 + (START_ANALYSIS_BEFORE_P_BRB_HP + MAX_MWPD_DUR) / deltaTime); // 20110331 AJL
        deData->mwpd->n_int_sum = 0;
        deData->mwpd->int_int_sum_pos = (double*) calloc(deData->mwpd->n_int_sum_max, sizeof (double));
        deData->mwpd->int_int_sum_neg = (double*) calloc(deData->mwpd->n_int_sum_max, sizeof (double));
        deData->mwpd->raw_mag = MWPD_INVALID;
        deData->mwpd->corr_mag = MWPD_INVALID;
        deData->mwpd->have_used_memory = 0;
    }
    deData->flag_complete_mwpd = 0;

    deData->brb_polarity = POLARITY_UNKNOWN;
    deData->brb_polarity_quality = 0.0;
    deData->brb_polarity_delay = 0.0;

    // event association
    //deData->num_assoc = -1;
    deData->is_associated = 0;
    deData->is_full_assoc_loc = -1; // -1 flags uninitialized
    //deData->is_associated_grid_level = 0;
    deData->epicentral_distance = -1.0;
    deData->epicentral_azimuth = -1.0;
    deData->residual = -999.0;
    deData->dist_weight = -1.0;
    deData->station_quality_weight = -1.0;
    deData->loc_weight = -1.0;
    strcpy(deData->phase, "X");
    deData->phase_id = -1;
    deData->take_off_angle_inc = -1.0;
    deData->take_off_angle_az = -1.0;
    // derived values
    deData->ttime_P = VALUE_NOT_SET;
    deData->ttime_S = VALUE_NOT_SET;
    deData->ttime_SminusP = VALUE_NOT_SET;


    // waveform export bookkeeping
    deData->waveform_export = NULL;
    if (waveform_export_enable) {
        deData->waveform_export = (WaveformExport*) calloc(1, sizeof (WaveformExport));
        deData->waveform_export->filename[0] = '\0';
        deData->waveform_export->start_time_written = -1;
        deData->waveform_export->end_time_written = -1;
        deData->waveform_export->hypo_unique_id = -1;
    }

    // amplitude attenuation
    deData->amplitude_error_ratio = -999.0;


    // station corrections  // 20150716 AJL - added to support output of used sta correction in pick csv file
    deData->sta_corr = 0.0;


    return (deData);

}

/** set data values */

TimedomainProcessingData* new_TimedomainProcessingData(char* sladdr, int n_int_tseries, int source_id, char* station, char* location, char* channel,
        char* network, double deltaTime, double lat, double lon, double elev, double station_quality_weight, PickData* pickData, int pick_stream, int init_ellapsed_index_count,
        int year, int month, int mday, int hour, int min, double dsec, int flag_do_mwpd, int waveform_export_enable) {

    TimedomainProcessingData* deData = init_TimedomainProcessingData(deltaTime, flag_do_mwpd, waveform_export_enable);

    deData->sladdr = sladdr;

    deData->n_int_tseries = n_int_tseries;

    deData->source_id = source_id;
    strcpy(deData->station, station);
    strcpy(deData->location, location);
    strcpy(deData->channel, channel);
    strcpy(deData->network, network);
    deData->lat = lat;
    deData->lon = lon;
    deData->elev = elev;
    deData->station_quality_weight = station_quality_weight;
    *(deData->pickData) = *pickData; // make local copy of pickData
    deData->pick_stream = pick_stream;
    deData->virtual_pick_index = init_ellapsed_index_count;
    deData->year = year;
    deData->month = month;
    deData->mday = mday;
    deData->hour = hour;
    deData->min = min;
    deData->dsec = dsec;

    // pick error
    strcpy(deData->pick_error_type, "GAU");
    // set uncertainty to half width between right and left indices
    deData->pick_error = deltaTime * fabs(pickData->indices[1] - pickData->indices[0]);
    deData->pick_error /= 2.0;
    if (deData->pick_error < deltaTime) { // 20101015 AJL - changed 0.0 to deltaTime
        deData->pick_error = deltaTime;
    }

    struct tm tm_time = {0};
    tm_time.tm_year = deData->year - 1900;
    tm_time.tm_mon = deData->month - 1;
    tm_time.tm_mday = deData->mday;
    tm_time.tm_hour = deData->hour;
    tm_time.tm_min = deData->min;
    tm_time.tm_sec = (int) deData->dsec;
    tm_time.tm_isdst = 0;
    tm_time.tm_gmtoff = 0;
    //tm_time.tm_zone = "UTC";
    deData->t_time_t = mktime(&tm_time);
    deData->t_decsec = deData->dsec - (double) ((int) deData->dsec);

    return (deData);

}

/** copy data to new data object */

TimedomainProcessingData* copy_TimedomainProcessingData(TimedomainProcessingData* deData) {

    TimedomainProcessingData* deDataCopy = init_TimedomainProcessingData(deData->deltaTime, deData->mwpd != NULL, deData->waveform_export != NULL);
    *(deDataCopy) = *(deData);

    /*strcpy(deDataCopy->station, deData->station);
    strcpy(deDataCopy->location, deData->location);
    strcpy(deDataCopy->channel, deData->channel);
    strcpy(deDataCopy->network, deData->network);*/

    deDataCopy->pickData = calloc(1, sizeof (PickData));
    *(deDataCopy->pickData) = *(deData->pickData); // make local copy of pickData

    return (deDataCopy);

}

/** clean up data memory */

void free_TimedomainProcessingData(TimedomainProcessingData* deData) {

    if (deData == NULL)
        return;

    free(deData->pickData);
    if (deData->mwpd != NULL) {
        free(deData->mwpd->int_int_sum_pos);
        free(deData->mwpd->int_int_sum_neg);
        free(deData->mwpd);
        deData->mwpd = NULL;
    }
    if (deData->mwp != NULL) {
        free(deData->mwp->peak);
        free(deData->mwp->peak_dur);
        free(deData->mwp->int_int);
        free(deData->mwp->int_int_abs);
        free(deData->mwp);
        deData->mwp = NULL;
    }
    if (deData->mb != NULL) {
        free(deData->mb->amplitude);
        free(deData->mb->peak_dur);
        free(deData->mb);
        deData->mb = NULL;
    }
    if (deData->grd_mot != NULL) {
        free(deData->grd_mot->vel);
        free(deData->grd_mot->disp);
        free(deData->grd_mot);
        deData->grd_mot = NULL;
    }
    free(deData->t0);

    // 20130510 AJL - bug fix (fixed?)
    // make sure other pointers to this data are NULL
    if (deData->use_for_location_twin_data != NULL) {
        deData->use_for_location_twin_data->use_for_location_twin_data = NULL;
        deData->use_for_location_twin_data = NULL;
    }

    if (deData->waveform_export != NULL) {
        free(deData->waveform_export);
    }

    free(deData);

}


/** add a TimedomainProcessingData to a TimedomainProcessingData list
 * list will be sorted by increasing deData->t_time_t
 */

#define SIZE_INCREMENT 16

int addTimedomainProcessingDataToDataList(TimedomainProcessingData* deData, TimedomainProcessingData*** pdata_list, int* pnum_de_data, int check_exists) {


    TimedomainProcessingData** newDataList = NULL;

    if (*pdata_list == NULL) { // list not yet created
        *pdata_list = calloc(SIZE_INCREMENT, sizeof (TimedomainProcessingData*));
        // 20130930 AJL - bug fix
        //} else if (*pnum_de_data != 0 && (*pnum_de_data % SIZE_INCREMENT) == 0) { // list will be too small
    } else if (*pnum_de_data + 1 > sizeof (*pdata_list) / sizeof (TimedomainProcessingData*)) { // list will be too small
        newDataList = calloc(*pnum_de_data + SIZE_INCREMENT, sizeof (TimedomainProcessingData*));
        int n;
        for (n = 0; n < *pnum_de_data; n++)
            newDataList[n] = (*pdata_list)[n];
        free(*pdata_list);
        *pdata_list = newDataList;
    }

    // find first data later than new data
    int ninsert;
    for (ninsert = 0; ninsert < *pnum_de_data; ninsert++)
        if (((*pdata_list)[ninsert]->t_time_t > deData->t_time_t)
                || ((*pdata_list)[ninsert]->t_time_t == deData->t_time_t && (*pdata_list)[ninsert]->t_decsec > deData->t_decsec) // 20101116 AJL - added check on t_decsec for same t_time_t
                )
            break;
    // check that data is not already present in list
    if (check_exists && ninsert > 0) {
        int ntest;
        for (ntest = 0; ntest < *pnum_de_data; ntest++) { // 20101116 AJL - added this loop, before just checked preceeding datum
            TimedomainProcessingData * pdata_test = (*pdata_list)[ntest];
            if (pdata_test->t_time_t < deData->t_time_t - 1) // 20101116 AJL - added   // 20101223 AJL - added -1 (at least +/-1 sec tolerance)
                continue;
            if (pdata_test->t_time_t > deData->t_time_t + 1) // 20101116 AJL - added   // 20101223 AJL - added +1 (at least +/-1 sec tolerance)
                break;
            // 20101223 AJL  if (pdata_test->source_id == deData->source_id && pdata_test->pick_stream == deData->pick_stream && pdata_test->t_decsec == deData->t_decsec)
            if (pdata_test->source_id == deData->source_id && pdata_test->pick_stream == deData->pick_stream) {
                if (0) {
                    printf("DEBUG: pdata_test, deData : \n");
                    fprintf_TimedomainProcessingData(pdata_test, stdout);
                    printf("\n");
                    fprintf_TimedomainProcessingData(deData, stdout);
                    printf("\n");
                }
                return (0);
            }
        }
    }
    // shift all later data to free ninsert position for new data
    if (ninsert < *pnum_de_data) {
        int m;
        for (m = *pnum_de_data - 1; m >= ninsert; m--)
            (*pdata_list)[m + 1] = (*pdata_list)[m];
    }
    // insert new TimedomainProcessingData
    (*pdata_list)[ninsert] = deData;
    (*pnum_de_data)++;

    return (1);

}

/** remove a TimedomainProcessingData from a TimedomainProcessingData list */

void removeTimedomainProcessingDataFromDataList(TimedomainProcessingData* deData, TimedomainProcessingData*** pdata_list, int* pnum_de_data) {

    int i = 0;

    while (i < *pnum_de_data) {
        if ((*pdata_list)[i] == deData)
            break;
        i++;
    }

    if (i == *pnum_de_data) // not found
        return;

    while (i < *pnum_de_data - 1) {
        (*pdata_list)[i] = (*pdata_list)[i + 1];
        i++;
    }

    (*pdata_list)[*pnum_de_data - 1] = NULL;

    (*pnum_de_data)--;

}

/** clean up list memory */

void free_TimedomainProcessingDataList(TimedomainProcessingData** data_list, int num_data) {
    if (data_list == NULL || num_data < 1)
        return;

    int n;
    for (n = 0; n < num_data; n++)
        free_TimedomainProcessingData(*(data_list + n));

    free(data_list);

}

/** check if de data is associated P */

int is_associated_phase(TimedomainProcessingData* deData) {
    return (deData->is_associated && deData->is_associated != NUMBER_ASSOCIATE_IGNORE);
}

/** check if de data is associated P */

int is_unassociated_or_location_P(TimedomainProcessingData* deData) {
    return (!deData->is_associated || is_associated_location_P(deData));
}

/** check if de data is associated P */

int is_associated_location_P(TimedomainProcessingData* deData) {
    // 20130218 AJL - include check for use_for_location
    //return (deData->is_associated && deData->is_associated != NUMBER_ASSOCIATE_IGNORE && is_direct_P(deData->phase_id));
    return (deData->can_use_as_location_P && deData->is_associated && deData->is_associated != NUMBER_ASSOCIATE_IGNORE && is_direct_P(deData->phase_id));
}

/** returns the name for a pick stream */

char* pick_stream_name(TimedomainProcessingData* deData) {

    if (deData->pick_stream == STREAM_HF)
        return (STREAM_HF_NAME);
    if (deData->pick_stream == STREAM_RAW)
        return (STREAM_RAW_NAME);

    return (STREAM_NULL_NAME);
}

/** print data */

int fprintf_TimedomainProcessingData(TimedomainProcessingData* deData, FILE* pfile) {
    if (pfile == NULL || deData == NULL)
        return (0);

    fprintf(pfile, "%d ", deData->source_id); //	int source_id;
    fprintf(pfile, "%s %s %s %s ", deData->station, deData->location, deData->channel, deData->network);
    fprintf(pfile, "%f %f %f ", deData->lat, deData->lon, deData->elev);

    fprintf(pfile, " ");
    fprintf_PickData(deData->pickData, pfile); //  PickData* pickData;   // initial pick
    fprintf(pfile, " ");
    fprintf(pfile, "%d %d %f %f %f %f %d %d %d %d %d %d %d %d %d  %f %ld %f ",
            deData->pick_stream, // int
            deData->virtual_pick_index, // int
            deData->sn_pick, deData->a_ref, deData->t50, deData->tauc_peak, // double
            deData->flag_complete_t50, deData->flag_a_ref_not_ok, deData->flag_snr_hf_too_low, deData->flag_clipped,
            deData->year, deData->month, deData->mday, deData->hour, deData->min, // int
            deData->dsec, // double
            deData->t_time_t, // long int in gcc
            deData->t_decsec // double
            );


    return (1);

}

/** format a double time (equivalent to time_t time with decimal seconds) into a time string with decimal seconds */

char *timeDecSec2NLLocstring(double time_dec_sec, char* str) {

    time_t itime_sec = (time_t) time_dec_sec;
    struct tm* tm_time = gmtime(&itime_sec);
    double dec_sec = time_dec_sec - itime_sec;

    // code data and time integers
    long int idate, ihrmin;
    idate = (tm_time->tm_year + 1900) * 10000 + (tm_time->tm_mon + 1) * 100 + tm_time->tm_mday;
    ihrmin = tm_time->tm_hour * 100 + tm_time->tm_min;

    sprintf(str, "%8.8ld %4.4ld %8.4lf", idate, ihrmin, (double) tm_time->tm_sec + dec_sec);

    return (str);
}

/** print data */

static char tmp_str[128];

int fprintf_NLLoc_TimedomainProcessingData(TimedomainProcessingData* deData, FILE* pfile, int append_evt_ndx_to_phase) {

    if (pfile == NULL || deData == NULL)
        return (0);

    // phase
    char phase[16];
    if (is_associated_phase(deData))
        if (append_evt_ndx_to_phase) {
            // 20130307 AJL - added "_" between phase and event index
            sprintf(phase, "%s_%d", deData->phase, deData->is_associated);
        } else {
            sprintf(phase, "%s", deData->phase);
        } else
        sprintf(phase, "X");

    // location
    char location[16] = "?";
    if (strlen(deData->location) > 0 && strcmp(deData->location, "--")) {
        strcpy(location, deData->location);
    }

    // first motion
    double fmquality = 0.0;
    int fmpolarity = POLARITY_UNKNOWN;
    char fmtype[32] = "Err";
    char first_mot = setPolarity(deData, &fmquality, &fmpolarity, fmtype);

    // error
    char error_type[] = "GAU";
    // set uncertainty to half width between right and left indices
    double error = deData->deltaTime * fabs(deData->pickData->indices[1] - deData->pickData->indices[0]);
    error /= 2.0;
    if (error < 0.0) {
        error = 0.0;
    }

    // misc
    double coda_dur = 0.0;
    double amplitude = deData->pickData->amplitude;
    double period = deData->pickData->period;
    double apriori_weight = 1.0;
    // test simple weighting scheme
    if (1) {
        apriori_weight = 0.0;
        if (deData->is_associated && deData->residual > -9999.9 && deData->residual < 9999.9) {
            apriori_weight = exp(-(fabs(deData->residual) / 1.0));
        }
    }

    // write observation part of NLL_FORMAT_VER_2 phase line
    //sprintf(pick_str,
    //	"%-6s %-4s %-4s %-1s %-6s %-1s %8.8ld %4.4ld %9.4lf %-3s %9.2le %9.2le %9.2le %9.2le %9.4lf",
    // write observation part of orig NLL phase line
    //fprintf(pfile,
    //        "%-6s %-4s %-4s %-1s %-6s %-1s %8.8ld %4.4ld %9.4lf %-3s %9.3le %9.3le %9.3le %9.3le",
    fprintf(pfile,
            "%-6s %-4s %-4s %-1s %-6s %-1c %s %-3s %9.2le %9.2le %9.2le %9.2le %9.4lf",
            deData->station,
            location,
            deData->channel,
            "?",
            phase,
            first_mot,
            timeDecSec2NLLocstring((double) deData->t_time_t + deData->t_decsec, tmp_str),
            error_type, error,
            coda_dur,
            amplitude,
            period,
            apriori_weight // NLL_FORMAT_VER_2
            );

    return (1);

}

/** get time_t structure and seconds as double for the data pick time */

time_t get_time_t(TimedomainProcessingData* deData, double* pdecsec) {

    *pdecsec = deData->t_decsec;

    return (deData->t_time_t);

}

/** set derived values in data object
 *
 * 20111224 AJL - added
 */
void set_derived_values(TimedomainProcessingData* deData, double hypo_depth) {

    deData->ttime_P = VALUE_NOT_SET;
    deData->ttime_S = VALUE_NOT_SET;
    deData->ttime_SminusP = VALUE_NOT_SET; // 20151117 AJL - added so can check in timedomain_processing.c->td_process_timedomain_processing() that s/n window does not pass S time
    // 20140305 AJL - changed condition to set P and S times for all possible waveform export pick data
    //if (is_associated_location_P(deData)) {
    if (is_associated_phase(deData) && (is_first_arrival_P(deData->phase_id) || is_direct_P(deData->phase_id)) && is_count_in_location(deData->phase_id)) {
        int phase_id_P = get_P_phase_index(); // IMPORTANT!: index must match phase order in ttimes.c array phases[]
        deData->ttime_P = get_ttime(phase_id_P, deData->epicentral_distance, hypo_depth);
        int phase_id_S = get_S_phase_index(); // IMPORTANT!: index must match phase order in ttimes.c array phases[]
        deData->ttime_S = get_ttime(phase_id_S, deData->epicentral_distance, hypo_depth);
        if (deData->ttime_P > 0.0 && deData->ttime_S > 0.0) {
            deData->ttime_SminusP = deData->ttime_S - deData->ttime_P;
        }
    }

    double duration_raw = deData->t0->duration_raw;
    double duration_corr = duration_raw;

#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
    double SP_FRACTION_CUTOFF = 0.25;
    // 20120107 AJL if (deData->t0->duration_raw != T0_INVALID && is_associated_location_P(deData)) {
    if (deData->t0->duration_raw != T0_INVALID && is_associated_location_P(deData) && deData->epicentral_distance < 30.0) {
        double ttime_P = deData->ttime_P;
        if (ttime_P > 0.0) {
            double ttime_S = deData->ttime_S;
            if (ttime_S > 0.0) {
                double tSminusP = ttime_S - ttime_P;
                // 20120107 AJL - test  if (ttime_P + duration_raw > ttime_S + tSminusP / 3.0) {
                // 20120221 AJL if (ttime_P + duration_raw > ttime_S + tSminusP) {
                // 20120221 AJL - test
                if (ttime_P + duration_raw > ttime_S + SP_FRACTION_CUTOFF * tSminusP) {

                    duration_corr -= tSminusP;
                    //printf("DEBUG: set_derived_values: ttime_P+duration_raw=%g  ttime_S+(ttime_S-ttime_P)/3.0=%g", ttime_P + duration_raw, ttime_S + tSminusP / 3.0);
                    //printf("  duration_raw=%g -> duration_corr=%g  dist=%g  %s\n", duration_raw, duration_corr, deData->epicentral_distance, (duration_raw == duration_corr ? "" : "CHANGED!"));
                    // 20120107 AJL - added
                    double duration_plot = calculate_duration_plot(deData->t0->duration_raw);
                    deData->t0->duration_plot = duration_plot;
                }
            }
        }
    }
#endif

    deData->t0->duration_corr = duration_corr;

}


#define MAX_P_GRD_MOT_AMP_WINDOW 5.0    // seconds

/** calculate initial P ground motion amplitudes
 *
 * 20140812 AJL - added
 */

double calculate_init_P_grd_mot_amp(TimedomainProcessingData* deData, double snr_brb, double snr_brb_int, int force_polarity, double *pfmquality, int *pfmpolarity, char *fmtype) {

    int index_pick = (int) (0.5 + START_ANALYSIS_BEFORE_P_BRB_HP / deData->deltaTime); // start amplitude check at P

    // set amp evaluation index to maximum available (this will be at the smaller of deData->grd_vel->n_amp, MAX_GRD_VEL_DUR, or deData->t0->duration_raw if set during real-time processing of data)
    int index_max = deData->grd_mot->n_amp;
    // limit to maximum amplitude window
    int max_amp_wind = (int) (0.5 + MAX_P_GRD_MOT_AMP_WINDOW / deData->deltaTime);
    if (index_max > index_pick + max_amp_wind)
        index_max = index_pick + max_amp_wind;
    // limit amp evaluation index to deData->t0->duration_raw
    double duration = deData->t0->duration_raw;
    // 20111230 AJL if (deData->flag_complete_t0 && duration != T0_INVALID) {
    if (duration != T0_INVALID) {
#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
        duration = deData->t0->duration_corr;
        if (duration < 0.0)
            printf("ERROR: initial_brb_bp_P_amplitude: deData->t0->duration_corr not set: this should not happen!\n");
#endif
        // if is associate P, limit amp window to S-P time
        if (is_associated_location_P(deData)) {
            double duration_before = duration;
            double ttime_P = deData->ttime_P;
            if (ttime_P == VALUE_NOT_SET) {
                printf("ERROR: initial_brb_bp_P_amplitude: deData->ttime_P not set: this should not happen! (%s_%s_%s_%s %s)\n",
                        deData->network, deData->station, deData->location, deData->channel, deData->phase);
            }
            if (ttime_P > 0.0) {
                double ttime_S = deData->ttime_S;
                if (ttime_S > 0.0) {
                    if (duration_before > ttime_S - ttime_P)
                        duration = ttime_S - ttime_P;
                    //printf("DEBUG: limit initial_brb_bp_P_amplitude window to S-P time: dist=%g  S-P: %g  duration: %g->%g\n",
                    //        deData->epicentral_distance, ttime_S - ttime_P, duration_before, duration);
                }
            }
        }

        int index_max_test = index_pick + (int) (0.5 + duration / deData->deltaTime); // end magnitude calculation at T0 duration
        if (index_max_test < index_max)
            index_max = index_max_test;
    }

    // find peak amplitude in window
    double amplitude_peak_vel = GRD_MOT_INVALID;
    double amplitude_peak_disp = GRD_MOT_INVALID;
    int pol_peak_disp = 0;
    double amplitude_vel_sum = 0.0;
    int namp_vel = 0;
    double amplitude_disp_sum = 0.0;
    int namp_disp = 0;
    double amp_vel, amp_disp;
    int pol_vel, pol_disp;
    int ndx;
    //int index_min = index_pick; // start amplitude check at P
    int index_min = 0; // start amplitude before P
    for (ndx = index_min; ndx <= index_max; ndx++) {
        if (force_polarity && *pfmpolarity == POLARITY_UNKNOWN)
            continue;
        amp_vel = deData->grd_mot->vel[ndx];
        if (0) { // use ground displacement
            amp_disp = deData->grd_mot->disp[ndx];
        } else { // use integral of ground displacement (moment)
            amp_disp = GRD_MOT_INVALID;
            if (ndx <= deData->mwp->n_peak) {
                amp_disp = deData->mwp->int_int[ndx];
            }
        }
        pol_vel = amp_vel >= 0.0 ? 1 : -1;
        if (!force_polarity || (*pfmpolarity == POLARITY_POS && pol_vel == 1) || (*pfmpolarity == POLARITY_NEG && pol_vel == -1)) {
            if (fabs(amp_vel) > amplitude_peak_vel) {
                amplitude_peak_vel = fabs(amp_vel);
            }
        }
        amplitude_vel_sum += fabs(amp_vel);
        namp_vel++;
        pol_disp = amp_disp >= 0.0 ? 1 : -1;
        if (amp_disp != GRD_MOT_INVALID) {
            double fabs_amp_disp = fabs(amp_disp);
            if (!force_polarity || (*pfmpolarity == POLARITY_POS && pol_disp == 1) || (*pfmpolarity == POLARITY_NEG && pol_disp == -1)) {
                if (fabs_amp_disp > amplitude_peak_disp) {
                    amplitude_peak_disp = fabs_amp_disp;
                    pol_peak_disp = pol_disp;
                }
            }
            amplitude_disp_sum += fabs_amp_disp;
            namp_disp++;
        }
    }

    // 20160316 AJL - gives incorrect polarity at very close distances   if (1 && !force_polarity && amplitude_peak_disp != GRD_MOT_INVALID) {
    if (0 && !force_polarity && amplitude_peak_disp != GRD_MOT_INVALID) {
        // 20140903 AJL - test using waveform polarity as pick polarity
        // use only if force_polarity == 0 !!!
        // TODO: replace Mwp waveform polarity estimate with this?
        // 20151130 AJL  if (snr_brb_int >= 1.1 * SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN) { // minimum quality will be 0.1
        if (snr_brb >= 1.1 * SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN
                && snr_brb_int >= 1.1 * SIGNAL_TO_NOISE_RATIO_BRB_HP_INT_MIN) { // minimum quality will be 0.1
            *pfmpolarity = pol_peak_disp;
            *pfmquality = snr_brb_int >= 2.0 * SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN ? 1.0 : (snr_brb_int - SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN) / SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN;
            sprintf(fmtype, "WMO");
        }
    }

    if (0) { // use mean amplitude
        if (namp_vel > 0) {
            amplitude_peak_vel = amplitude_vel_sum / (double) namp_vel;
        }
        if (namp_disp > 0) {
            amplitude_peak_disp = amplitude_disp_sum / (double) namp_disp;
        }
    }
    //printf("DEBUG: %s_%s, index %d-%d, duration %f, amplitude_peak_vel %e, amplitude_peak_disp %e\n",
    //        deData->network, deData->station, index_min, index_max, duration, amplitude_peak_vel, amplitude_peak_disp);
    deData->grd_mot->peak_amp_vel = amplitude_peak_vel;
    deData->grd_mot->peak_amp_disp = amplitude_peak_disp;

    return (amplitude_peak_vel);

}

/** Applies a magnitude correction to mb
 *
 *     mb_uncorr        uncorrected EE mb
 *
 * From 20150406 proof of :
 *    Bernardi et al. (2015) Early-Est: two years data-set analysis and application to early tsunami warnings.
 * Regression:
 *    the coefficients are:
 *    mb_neic = a * mb_uncorr + b
 *    a = 0.52
 *    b = 2.46
 *
 *    Since regression of Bernardi et al. (2015) only done for M>5:
 *    Only correct downwards, no corrections for mb_uncorr < crossover  between regression line and line mb_neic == mb_uncorr
 *    mb_neic == mb_uncorr for mb_uncorr = -b / (a-1) = 5.125
 */
double calculate_corrected_mB_INGV_EE(double mb_uncorr) {

    if (mb_uncorr < 5.125) {
        return (mb_uncorr);
    }

    double mb_neic = 0.52 * mb_uncorr + 2.46;

    return (mb_neic);
}

/** Applies a distance correction to Mwp
 *
 *     distance        station epicentral distance in degrees
 *
 *     depth           hypocenter depth in km
 *
 *     Return:         distance correction to be subtracted from Mwp
 *
 * 20150210 AJL - first added to code
 * 20150312 AJL - used coefficient from final version of paper.
 *
 * From:
 *    Bernardi et al. (2015) Appraising the Early-est earthquake monitoring system for tsunami alerting at the Italian candidate Tsunami Service Provider
 *    Nat. Hazards Earth Syst. Sci. Discuss., 15, 1–40, 2015 www.nat-hazards-ear th-syst-sci-discuss.net/15/1/2015/ doi:10.5194/nhessd-15-1-2015

 *    for depth<100km
 *
 *  f(d)= -1.32e-6*d^3 +2.40e-4*d^2 -0.0146*d + 0.314
 */
double calculate_Mwp_correction_INGV_EE(double distance, double depth) {

    if (distance < 0.0 || depth > 100.0) {
        return 0.0;
    }

    //double correction = -1.45e-6 * pow(distance, 3) + 2.73e-4 * pow(distance, 2) - 0.0171 * distance + 0.335;   // 20150210
    double correction = -1.32e-6 * pow(distance, 3) + 2.40e-4 * pow(distance, 2) - 0.0146 * distance + 0.314; // 20150312

    return (correction);
}


// 20141201 AJL - added following Mwp distance correction
// from utils_ATWC.c

/***********************************************************************
 *                         calculate_Mwp_correction_ATWC()                             *
 *                                                                     *
 * This function corrects the computed Mwp based on epicentral         *
 * distance.  The correction values are based on Paul Huang's study of *
 * Mwp versus Global CMT over all distance ranges.                     *
 *                                                                     *
 *  Arguments:                                                         *
 *     distance        Station epicentral distance in degrees         *
 *                                                                     *
 *  Return:            Distance correction value to subtract from Mwp *
 *                                                                     *
 ***********************************************************************/
double calculate_Mwp_correction_ATWC(double distance) {

    double dMwp_distance_correction[50] = {
        0.280154895, 0.248159091, 0.192827333, 0.113355523, 0.134518163,
        0.142274094, 0.141015269, 0.225545292, 0.270142832, 0.292807679,
        0.279098314, 0.083464936, 0.032658602, 0.017405582, 0.03411125,
        -0.01300572, -0.087751188, -0.105951149, -0.159169491, -0.109741881,
        -0.040873702, -0.127091897, -0.0291207, -0.051879539, -0.063962226,
        0.000433825, 0.020766625, -0.036926385, -0.042950165, -0.006198259,
        0.013929801, 0.007052744, 0.087026743, -0.007465525, -0.021245992,
        0.039889117, 0.052830981, -0.005707771, -0.022900097, -0.009290172,
        -0.028246945, -0.032688366, -0.020939191, 0.096496123, -0.025355755,
        -0.04532547, -0.048598167, -0.105345188, -0.309972643, -0.299542571
    };

    if (distance < 0.0)
        return 0.0;

    int iIndex; // Index of above array based on distance
    iIndex = (int) (distance / 2);
    if (iIndex < 0 || iIndex >= 50)
        return 0.0;

    return dMwp_distance_correction[iIndex];
}

/** Applies a distance correction to Mwp
 *
 *     distance        Station epicentral distance in degrees
 *
 *     Return:         distance correction to be subtracted from Mwp
 *
 *  * 20140812 AJL - added
 */
double calculate_Mwp_correction(double distance, double depth) {

    static int done = 0;
    if (0 && !done) {
        for (int i = 0; i < 100; i += 2) {
            if (i % 5 == 0) {
                printf("\n");
            }
            printf("%f, ", calculate_Mwp_correction_INGV_EE((double) i, 0.0));
        }
        printf("\n");
    }
    done = 1;

    //return (0.0);
    return (calculate_Mwp_correction_INGV_EE(distance, depth)); // 20150210 AJL - activated
    //return (calculate_Mwp_correction_ATWC(distance));    // used from 20141201 -> 20150210
}



//#define PI 3.14159265359
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define PI M_PI

// 4.213e19 - Tsuboi 1995, 1999
static double MWP_CONST =
        4.0 * PI // 4 PI
        * 3400.0 // rho
        * 7900.0 * 7900.0 * 7900.0 // Pvel**3
        * 2.0 // FP average radiation pattern
        * (10000.0 / 90.0) // distance deg -> km
* 1000.0 // distance km -> meters
;

/** calculate Mwp magnitude  */

double calculate_Mwp_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mwp_distance_correction, int apply_mag_atten_check) {

    // check if amplitude attenuation error ratio too high - gain may be incorrect, ignore magnitude
    if (apply_mag_atten_check) {
        if (deData->amplitude_error_ratio > 0.0
                && (deData->amplitude_error_ratio < AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO
                || deData->amplitude_error_ratio > AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO)) {
            return (MWP_INVALID);
        }
    }

    // set mag evaluation index to maximum available (this will be at the smaller of deData->mwp->n_peak, MAX_MWP_DUR,
    //   or deData->t0->duration_raw if set during real-time processing of data)
    int index_max = deData->mwp->n_peak;
    // limit mag evaluation index to deData->t0->duration_raw
    double duration = deData->t0->duration_raw;
    // 20111230 AJL if (deData->flag_complete_t0 && duration != T0_INVALID) {
    if (duration != T0_INVALID) {
#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
        duration = deData->t0->duration_corr;
        if (duration < 0.0)
            printf("ERROR: calculate_Mwp_Mag: deData->t0->duration_corr not set: this should not happen!\n");
#endif
        if (1) { // 20111224 TEST AJL
            // limit Mwp integration window to S-P time
            double duration_before = duration;
            double ttime_P = deData->ttime_P;
            if (ttime_P == VALUE_NOT_SET) {
                printf("ERROR: calculate_Mwp_Mag: deData->ttime_P not set: this should not happen! (%s_%s_%s_%s %s)\n",
                        deData->network, deData->station, deData->location, deData->channel, deData->phase);
            }
            if (ttime_P > 0.0) {
                double ttime_S = deData->ttime_S;
                if (ttime_S > 0.0) {
                    if (duration_before > ttime_S - ttime_P)
                        duration = ttime_S - ttime_P;
                    //printf("DEBUG: limit Mwpd integration window to S-P time: dist=%g  S-P: %g  duration: %g->%g\n",
                    //        deData->epicentral_distance, ttime_S - ttime_P, duration_before, duration);
                }
            }
        } // END - 20111224 TEST AJL

        int index_max_test = (int) (0.5 + (duration + START_ANALYSIS_BEFORE_P_BRB_HP) / deData->deltaTime); // end magnitude calculation at T0 duration
        if (index_max_test < index_max)
            index_max = index_max_test;
    }
    double peak = deData->mwp->peak[index_max];
    deData->mwp->index_mag = index_max;

    // 20120312 TEST AJL - use Mwp integral for first motion
    deData->brb_polarity = POLARITY_UNKNOWN;
    deData->brb_polarity_delay = 0.0;
    deData->brb_polarity_quality = 0.0;
    if (USE_MWP_INT_FOR_FIRST_MOTION) {
        if (deData->mwp->peak_dur[index_max] >= MWP_POLARITY_PULSE_WIDTH_MIN) {
            // polarity is sign of Mwp integral at time after P = fract *  Mwp duration
            // note that even if Mwp peak is secondary P pulse, polarity is alway taken immediately after P to capture polarity of first main P pulse
            double mwp_polarity_delay = deData->mwp->peak_dur[index_max] * MWP_POLARITY_PULSE_FRACTION;
            if (mwp_polarity_delay > MWP_POLARITY_DELAY_MAX) // make sure delay is not too long (e.g. for very large events)
                mwp_polarity_delay = MWP_POLARITY_DELAY_MAX;
            int polarity_index = (int) (0.5 + (mwp_polarity_delay + START_ANALYSIS_BEFORE_P_BRB_HP) / deData->deltaTime);
            double int_int_abs = deData->mwp->int_int_abs[polarity_index];
            if (int_int_abs > FLT_MIN) {
                double int_int = deData->mwp->int_int[polarity_index];
                deData->brb_polarity = int_int > 0.0 ? 1 : -1;
                deData->brb_polarity_quality = fabs(int_int) / int_int_abs;
                deData->brb_polarity_delay = mwp_polarity_delay;
                //printf("$$>>> %s_%s, pulse %f, delay %f, int %e, int_abs %e, pol %d, qual %f\n",
                //        deData->network, deData->station, deData->mwp->peak_dur[index_max], mwp_polarity_delay, int_int, int_int_abs,
                //        deData->brb_polarity, deData->brb_polarity_quality);
                //if (deData->brb_polarity_quality < 0.0)
                //    deData->brb_polarity_quality = 0.0;
                //else if (deData->brb_polarity_quality > 1.0)
                //    deData->brb_polarity_quality = 1.0;
            }
        }
    }

    // 20120319 AJL - use Mwp duration for tauc
    // DO NOT USE! - much worse results than tauc, mainly captures tsunami earthquakes,
    //      but degrades other tsunamigeneic events and does not filter non-tsunamigenic events
    /*
    if (0) {
        deData->tauc_peak = deData->mwp->peak_dur[index_max];
    } */

    // 20120106 AJL - Test
    /*
    double tauc_max = 0.0;
    int ndx;
    double dur;
    int tauc_index_min = (int) (0.5 + (TIME_DELAY_TAUC_MIN + START_ANALYSIS_BEFORE_P_BRB_HP) / deData->deltaTime);
    if (tauc_index_min > index_max)
        tauc_index_min = index_max;
    int tauc_index_max = (int) (0.5 + (TIME_DELAY_TAUC_MAX + START_ANALYSIS_BEFORE_P_BRB_HP) / deData->deltaTime);
    if (tauc_index_max > index_max)
        tauc_index_max = index_max;
    for (ndx = tauc_index_min; ndx < tauc_index_max; ndx++) {
        dur = deData->mwp->peak_dur[ndx];
        if (dur > 0.0 && dur > tauc_max) {
            tauc_max += dur;
        }
    }
    deData->tauc_peak = tauc_max;
     */

    double moment = MWP_CONST * peak * deData->epicentral_distance;

    double mwp_mag = MWP_INVALID;
    if (moment > FLT_MIN) {
        mwp_mag = (2.0 / 3.0) * (log10(moment) - 9.1);
        if (use_mwp_distance_correction) {
            mwp_mag -= calculate_Mwp_correction(deData->epicentral_distance, hypo_depth);
        }
    }

    deData->mwp->mag = mwp_mag;

    return (mwp_mag);

}

/** set pick polarity and weight based on HF first motion and BRB waveform polarity */

char setPolarity(TimedomainProcessingData* deData, double *pfmquality, int *pfmpolarity, char *fmtype) {

    *pfmquality = 0.0;
    *pfmpolarity = POLARITY_UNKNOWN;
    strcpy(fmtype, "F\0                        ");

    if (deData->pickData->polarity != POLARITY_UNKNOWN) {
        *pfmpolarity = deData->pickData->polarity;
        *pfmquality = deData->pickData->polarityWeight;
        sprintf(fmtype, "F");
    }
    // 20120312 AJL - use Mwp integral for first motion
    if (USE_MWP_INT_FOR_FIRST_MOTION) {
        if (!deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low && deData->brb_polarity != POLARITY_UNKNOWN
                && deData->brb_polarity_quality > *pfmquality
                && deData->brb_polarity_quality > 0.333
                ) {
            *pfmpolarity = deData->brb_polarity;
            *pfmquality = deData->brb_polarity_quality;
            // 201210818 TEST AJL - use snr for first motion quality
            //if (snr_brb < 2.0 * SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN)
            //    *pfmquality = (snr_brb - SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN) / SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN;
            //else
            //    *pfmquality = 1.0;
            sprintf(fmtype, "W%.1f", deData->brb_polarity_delay);
        }
    }

    if (*pfmquality < 0.0)
        *pfmquality = 0.0;

    char first_mot = '?';
    if (*pfmpolarity == POLARITY_POS) {
        first_mot = '+';
    } else if (*pfmpolarity == POLARITY_NEG) {
        first_mot = '-';
    }

    return (first_mot);

}



// calculated following formulas in Lomax & Michelini 2008
/*static double MWPD_CONST =
        1.886E+019 // from net.alomax.seisgram2k.calc.AmplitudeDurationMagnitudeCalculator w/ regression corr from mda6_010_GJI_rev.ods 20080318
 * (10000.0 / 90.0); // distance deg -> km
 */
// since not using geometrical spreading function, better to use Tsuboi constant ??
#define MWPD_CONST MWP_CONST

/** calculate Mwp magnitude
 * calculated following formulas in Lomax & Michelini 2008
 */
double calculate_Raw_Mwpd_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mwp_distance_correction, int apply_mag_atten_check) {

    //printf("MWP_CONST %g  MWPD_CONST %g\n", MWP_CONST, MWPD_CONST);

    double raw_mwpd_mag = MWPD_INVALID;

    // check if amplitude attenuation error ratio too high - gain may be incorrect, ignore magnitude
    if (apply_mag_atten_check) {
        if (deData->amplitude_error_ratio > 0.0
                && (deData->amplitude_error_ratio < AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO
                || deData->amplitude_error_ratio > AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO)) {
            return (MWPD_INVALID);
        }
    }

    int index_duration = 0;
    double duration = deData->t0->duration_raw;
    if (!deData->flag_complete_t0 && duration == T0_INVALID) {
        // calculates Mwpd before TO is completed
        static const int TEST = 1; // 20110908 TEST AJL
        if (TEST) {
            index_duration = deData->mwpd->n_int_sum; // 20110411 AJL - calculates Mwpd before TO is completed
            if (index_duration <= (int) (0.5 + START_ANALYSIS_BEFORE_P_BRB_HP / deData->deltaTime))
                return (MWPD_INVALID);
        } else {
            return (MWPD_INVALID);
        }
        // 20111230 AJL     } else if (deData->flag_complete_t0 && duration != T0_INVALID) {
    } else if (duration != T0_INVALID) {
#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
        duration = deData->t0->duration_corr;
        if (duration < 0.0)
            printf("ERROR: calculate_Raw_Mwpd_Mag: deData->t0->duration_corr not set: this should not happen!\n");
#endif
        if (1) { // 20110908 TEST AJL
            // limit Mwpd integration window to S-P time
            double duration_before = duration;
            double ttime_P = deData->ttime_P;
            if (ttime_P == VALUE_NOT_SET)
                printf("ERROR: calculate_Raw_Mwpd_Mag: deData->ttime_P not set: this should not happen! (%s_%s_%s_%s %s)\n",
                    deData->network, deData->station, deData->location, deData->channel, deData->phase);
            if (ttime_P > 0.0) {
                double ttime_S = deData->ttime_S;
                if (ttime_S > 0.0) {
                    if (duration_before > ttime_S - ttime_P)
                        duration = ttime_S - ttime_P;
                    //printf("DEBUG: limit Mwpd integration window to S-P time: dist=%g  S-P: %g  duration: %g->%g\n",
                    //        deData->epicentral_distance, ttime_S - ttime_P, duration_before, duration);
                }
            }
        } // END - 20110908 TEST AJL
        if (duration > MAX_MWPD_DUR)
            duration = MAX_MWPD_DUR;
        index_duration = (int) (0.5 + (duration + START_ANALYSIS_BEFORE_P_BRB_HP) / deData->deltaTime);
    } else {
        return (MWPD_INVALID);
    }

    if (index_duration >= deData->mwpd->n_int_sum_max) {
        printf("ERROR: calculate_Raw_Mwpd_Mag: dt=%g  Mwpd: index_duration=%d  >= deData->mwpd->n_int_sum_max=%d: this should not happen!\n",
                deData->deltaTime, index_duration, deData->mwpd->n_int_sum_max);
    }

    double amplitudeIntegralPos = deData->mwpd->int_int_sum_pos[index_duration];
    double amplitudeIntegralNeg = deData->mwpd->int_int_sum_neg[index_duration];
    // use greater of pos or neg displacement integral
    double amplitudeIntegral = amplitudeIntegralPos;
    if (amplitudeIntegralNeg > amplitudeIntegralPos)
        amplitudeIntegral = amplitudeIntegralNeg;

    double moment = MWPD_CONST * amplitudeIntegral;
    moment *= deData->epicentral_distance;

    if (moment > FLT_MIN) {
        raw_mwpd_mag = (2.0 / 3.0) * (log10(moment) - 9.1);
        if (use_mwp_distance_correction) {
            raw_mwpd_mag -= calculate_Mwp_correction(deData->epicentral_distance, hypo_depth);
        }
    }
    deData->mwpd->raw_mag = raw_mwpd_mag;

    return (raw_mwpd_mag);

}

// T0 and Td moment corrections (see Mwpd(To,Td) in duration-amplitude spreadsheet, e.g. Mwpd_To_TauC_002.ods->TauC 20101205)
//#define MWPD_MOMENT_CORR_DUR_CUTOFF 100.0    // 20110117 AJL
#define MWPD_MOMENT_CORR_DUR_CUTOFF_LOW 90.0    // 20120418 AJL
#define MWPD_MOMENT_CORR_DUR_CUTOFF_HIGH 110.0    // 20120418 AJL
//#define MWPD_MOMENT_CORR_DUR_CUTOFF 80.0    // 20111023 AJL
#define USE_MPWD_CORR 1
//
// T0 and Td moment corrections (see Mwpd(To,Td) in duration-amplitude spreadsheet, e.g. Mwpd_To_TauC_004.ods->TauC 20110316)
//#define MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC 120.0    // 20110117 AJL
//#define MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC 100.0    // 20110121 AJL
//#define MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC 90.0    // 20110121 AJL
#define MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC 80.0    // 20110121 AJL
#define MWPD_MOMENT_CORR_TAUC_CUTOFF_TAUC 8.0
//#define MWPD_MOMENT_CORR_TAUC_CUTOFF_TAUC 10.0    // 20110121 AJL
#define USE_MPWD_CORR_TAUC 0
//
#define MWPD_MOMENT_CORR_MAG_CUTOFF 7.2
#define MWPD_MOMENT_CORR_POW 0.45
#define USE_MOMENT_CORR 0

/** calculate corrected Mwpd magnitude
 * calculated following formulas in Lomax & Michelini 2008
 */
double calculate_corrected_Mwpd_Mag(double raw_mwpd_mag, double duration, double dominant_period, double depth) {

    //printf("MWP_CONST %g  MWPD_CONST %g\n", MWP_CONST, MWPD_CONST);

    double corr_mwpd_mag = raw_mwpd_mag;

    corr_mwpd_mag += get_depth_corr_mwpd_prem(depth);

    // T0 moment corrections (see Mwpd(To,Td) in duration-amplitude spreadsheet, e.g. Mwpd_To_TauC_002.ods->TauC 20101205)
    if (USE_MPWD_CORR) {
        if (duration > MWPD_MOMENT_CORR_DUR_CUTOFF_LOW) {
            // prevent sharp change in Mwpd at MWPD_MOMENT_CORR_DUR_CUTOFF, introduce linear ramp correction
            double dur_cutoff_ramp_level = (duration - MWPD_MOMENT_CORR_DUR_CUTOFF_LOW) / (MWPD_MOMENT_CORR_DUR_CUTOFF_HIGH - MWPD_MOMENT_CORR_DUR_CUTOFF_LOW);
            if (dur_cutoff_ramp_level > 1.0)
                dur_cutoff_ramp_level = 1.0;
            //printf("DEBUG: dur_cutoff_ramp_level=%f, duration %f, MWPD_MOMENT_CORR_DUR_CUTOFF_LOW %f, MWPD_MOMENT_CORR_DUR_CUTOFF_HIGH %f\n",
            //        dur_cutoff_ramp_level, duration, MWPD_MOMENT_CORR_DUR_CUTOFF_LOW, MWPD_MOMENT_CORR_DUR_CUTOFF_HIGH);
            // moment correction version
            // 20120418 AJL corr_mwpd_mag += (raw_mwpd_mag - MWPD_MOMENT_CORR_MAG_CUTOFF) * MWPD_MOMENT_CORR_POW;
            corr_mwpd_mag += dur_cutoff_ramp_level * (raw_mwpd_mag - MWPD_MOMENT_CORR_MAG_CUTOFF) * MWPD_MOMENT_CORR_POW;
            // duration correction version
            //corr_mwpd_mag += (2.0 / 3.0) * log10(duration / MWPD_MOMENT_CORR_DUR_CUTOFF);
        }
    }

    // T0 and Td moment corrections (see Mwpd(To,Td) in duration-amplitude spreadsheet, e.g. Mwpd_To_TauC_004.ods->TauC 20110316)
    if (USE_MPWD_CORR_TAUC) {
        if (duration > MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC)
            corr_mwpd_mag += (2.0 / 3.0) * log10(duration / MWPD_MOMENT_CORR_DUR_CUTOFF_TAUC);
        if (dominant_period > MWPD_MOMENT_CORR_TAUC_CUTOFF_TAUC)
            corr_mwpd_mag += (2.0 / 3.0) * log10(dominant_period / MWPD_MOMENT_CORR_TAUC_CUTOFF_TAUC);
    }

    // moment correction from Lomax & Michelini 2008
    if (USE_MOMENT_CORR) {

        if (raw_mwpd_mag > MWPD_MOMENT_CORR_MAG_CUTOFF)
            corr_mwpd_mag += (raw_mwpd_mag - MWPD_MOMENT_CORR_MAG_CUTOFF) * MWPD_MOMENT_CORR_POW;
    }

    return (corr_mwpd_mag);

}

/***************************************************************************
 * td_calculateDuration:
 *
 * returns duration calculated following formulas in Lomax & Michelini 2008, Appendix B.
 *
 * adapted from java code net.alomax.seisgram2k.calc.AmplitudeDurationMagnitudeCalculator.calculateDuration()
 *
 ***************************************************************************/

double calculate_duration(TimedomainProcessingData* deData) {

    T0Data* t0 = deData->t0;
    double deltaTime = deData->deltaTime;

    double ttime_P = 0.0; // T0Data index values are relative to P index
    double time_T09 = (double) t0->index_90 * deltaTime;
    double time_T08 = (double) t0->index_80 * deltaTime;
    double time_T05 = (double) t0->index_50 * deltaTime;
    double time_T02 = (double) t0->index_20 * deltaTime;

    // calculate duration
    double duration_estimate = ((time_T05 + time_T08) / 2.0) - ttime_P;
    double factor = (duration_estimate - 20.0) / 40.0; // factor 0->1 for duration 20->60
    if (factor < 0.0) {
        factor = 0.0;
    } else if (factor > 1.0) {
        factor = 1.0;
    }
    double duration = (time_T09 - ttime_P) * (1.0 - factor) + (time_T02 - ttime_P) * factor;
    if (duration < MINIMUN_DURATION)
        duration = MINIMUN_DURATION;

    deData->t0->duration_raw = duration;

    double duration_plot = calculate_duration_plot(deData->t0->duration_raw);
    deData->t0->duration_plot = duration_plot;

    return (duration);

}

/***************************************************************************
 * calculate_duration_plot:
 *
 * returns duration corrected by smoothing window width for plotting.
 * *
 ***************************************************************************/

double calculate_duration_plot(double duration_raw) {

    double smoothing_window = 2.0 * SMOOTHING_WINDOW_HALF_WIDTH_T50;

    // remove smoothing window width for short durations  // 20111222 AJL
    // 20120222 AJL factor = (duration - smoothing_window) / (2.0 * smoothing_window); // factor 0-1 for duration smoothing_window -> 3*smoothing_window
    double factor = (duration_raw - smoothing_window) / (2.0 * smoothing_window); // factor 0->1 for duration= smoothing_window -> 3*smoothing_window
    if (factor < 0.0) {
        factor = 0.0;
    } else if (factor > 1.0) {
        factor = 1.0;
    }
    double duration_plot = duration_raw - smoothing_window * (1.0 - factor);
    if (duration_plot < MINIMUN_DURATION)
        duration_plot = MINIMUN_DURATION;

    return (duration_plot);

}

/** calculate mb magnitude
 *
 * calulated magnitude will be mb if low-pass corner is at 3sec
 * calulated magnitude will be mB if low-pass corner is at 30sec
 *
 * from 2008__Bormann_Saul__The_New_IASPEI_Standard_Broadband__Magnitude_mB__SRL.pdf
 *
 *      mB = log10 (Vmax/2π)+ Q(∆, h)
 *
 *
 * from: 2009__Bormann_et__Summary_of_Magnitude_Working_Group_Recommendations_on_Standard_Procedures_for_Determining_Earthquake_Magnitudes_from_Digital_Data__IASPEI.pdf
 *
 *      mb – short-period body-wave magnitude (5) mb = log10(A/T) + Q(Δ, h)
 *      where,
 *      A = P-wave ground displacement in μm calculated from the maximum trace- amplitude in the entire P-phase train (time spanned by P, pP, sP, and possibly PcP and their codas, and ending preferably before PP);
 *      T = period in seconds, T < 3s;
 *      A and T are measured on output from a vertical-component instrument that is filtered so that the frequency response of the seismograph/filter system replicates that of a WWSSN short-period seismograph;
 *      Q(Δ, h) = attenuation function for PZ (P-waves recorded on vertical component seismographs) established by Gutenberg and Richter (1956);
 *      Δ = epicentral distance in degrees, 21 ̊ ≤ Δ ≤ 100°;
 *      h = focal depth.
 */
double calculate_mb_Mag(TimedomainProcessingData* deData, double hypo_depth, int apply_mag_atten_check) {

    // check if amplitude attenuation error ratio too high - gain may be incorrect, ignore magnitude
    if (apply_mag_atten_check) {
        if (deData->amplitude_error_ratio > 0.0
                && (deData->amplitude_error_ratio < AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO
                || deData->amplitude_error_ratio > AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO)) {
            return (MB_INVALID);
        }
    }

    // set maximum amplitude index to maximum available (this will be at the smaller of deData->mb->n_amplitude, MAX_MB_DUR, or deData->t0->duration_raw if set during real-time processing of data)
    int index_max = deData->mb->n_amplitude;
    // limit amplitude index to deData->t0->duration_raw
    double duration = deData->t0->duration_raw;
    // 20111230 AJL         if (deData->flag_complete_t0 && duration != T0_INVALID) {
    if (duration != T0_INVALID) {
#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
        duration = deData->t0->duration_corr;
        if (duration < 0.0)
            printf("ERROR: calculate_mb_Mag: deData->t0->duration_corr not set: this should not happen!\n");
#endif
        if (1) { // 20111224 TEST AJL
            // limit mb analyis window to S-P time
            double duration_before = duration;
            double ttime_P = deData->ttime_P;
            if (ttime_P == VALUE_NOT_SET)
                printf("ERROR: calculate_mb_Mag: deData->ttime_P not set: this should not happen! (%s_%s_%s_%s %s)\n",
                    deData->network, deData->station, deData->location, deData->channel, deData->phase);
            if (ttime_P > 0.0) {
                double ttime_S = deData->ttime_S;
                if (ttime_S > 0.0) {
                    if (duration_before > ttime_S - ttime_P)
                        duration = ttime_S - ttime_P;
                    //printf("DEBUG: limit Mwpd integration window to S-P time: dist=%g  S-P: %g  duration: %g->%g\n",
                    //        deData->epicentral_distance, ttime_S - ttime_P, duration_before, duration);
                }
            }
        } // END - 20111224 TEST AJL

        int index_max_test = (int) (0.5 + (duration + START_ANALYSIS_BEFORE_P_MB) / deData->deltaTime); // end magnitude calcaulation at T0 duration
        if (index_max_test < index_max)
            index_max = index_max_test;
    }
    double peak = 0.0;
    double peak_dur = 1.0; // 1 sec default
    int peak_index = -1;
    int n;
    //printf("DEBUG: mb->peak_dur[n]: ");
    for (n = 0; n < index_max; n++) {
        //printf(" %g,%g", deData->mb->amplitude[n], deData->mb->peak_dur[n]);
        // check for peak in amplitude and valid peak duration (prevents using amplitudes from incomplete displacement pulse at end of MAX_MB_DUR)
        if (deData->mb->amplitude[n] > peak && deData->mb->peak_dur[n] > 0.0) {
            peak = deData->mb->amplitude[n];
            peak_dur = deData->mb->peak_dur[n];
            peak_index = n;
        }
        //if (strcmp(deData->station, "HGN") == 0)
        //    printf("%d/%g/%g  ", n, deData->mb->amplitude[n], deData->mb->peak_dur[n]);
    }
    //if (strcmp(deData->station, "HGN") == 0)
    //   printf("\n");

    double mb_mag = MB_INVALID;

    //double period = 2.0 * peak_dur; // velocity T is 2 pulses
    double period = peak_dur; // displacement T is 1 pulse width
    double amp_value = peak / period;
    //double Q_value = get_Q_Delta_PV_value(deData->epicentral_distance, hypo_depth);
    double Q_value = get_Q_Depth_Delta_PV_value(deData->epicentral_distance, hypo_depth);
    if (period >= 0.2 && period <= 3.0) {
        if (amp_value > FLT_MIN && Q_value > FLT_MIN) {
            mb_mag = log10(amp_value) + Q_value;
        }
    }

    deData->mb->mag = mb_mag;

    //printf("DEBUG: @@@@  %s_%s : mb_mag: %g = log10(amp_value %g)%g + Q_value %g (distance=%g, depth=%g) : t_max=%gs : peak=%g 2*peak_dur=%g  peak_index=%d, apick=%g\n",
    //        deData->network, deData->station, mb_mag, amp_value, log10(amp_value), Q_value, deData->epicentral_distance, hypo_depth,
    //        (double) index_max * deData->deltaTime, peak, 2.0 * peak_dur, peak_index, deData->mb->amp_at_pick);

    return (mb_mag);

}

/** calculate mB magnitude
 *
 * calulated magnitude will be mb if low-pass corner is at 3sec
 * calulated magnitude will be mB if low-pass corner is at 30sec
 *
 * from 2008__Bormann_Saul__The_New_IASPEI_Standard_Broadband__Magnitude_mB__SRL.pdf
 *
 *      mB = log10 (Vmax/2π)+ Q(∆, h)
 *
 *
 * from: 2009__Bormann_et__Summary_of_Magnitude_Working_Group_Recommendations_on_Standard_Procedures_for_Determining_Earthquake_Magnitudes_from_Digital_Data__IASPEI.pdf
 *
 *      mB(BB) –broadband body-wave magnitude (6)
 *      mBB = log10(A/T)max + Q(Δ, h)
 *      where, (A/T)max = (Vmax/2π),
 *      where Vmax = ground velocity in μm/s associated with the maximum trace-amplitude in the entire P-phase train (time spanned by P, pP, sP,
 *      and possibly PcP and their codas, but ending preferably before PP) as recorded on a vertical-component seismogram
 *      that is proportional to velocity in the period-range 0.2s < T < 30s, and where T should be preserved together with A or Vmax in bulletin data-bases;
 *      Q(Δ, h) = attenuation function for PZ established by Gutenberg and Richter (1956), as discussed above with respect to mb;
 *      Δ = epicentral distance in degrees, 21 ̊ ≤ Δ ≤ 100°;
 *      h = focal depth.
 *
 */
double calculate_mB_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mb_correction, int apply_mag_atten_check) {

    if (apply_mag_atten_check) {
        if (deData->amplitude_error_ratio > 0.0
                && (deData->amplitude_error_ratio < AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO
                || deData->amplitude_error_ratio > AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO)) {
            return (MB_INVALID);
        }
    }

    // set maximum amplitude index to maximum available (this will be at the smaller of deData->mb->n_amplitude, MAX_MB_DUR, or deData->t0->duration_raw if set during real-time processing of data)
    int index_max = deData->mb->n_amplitude;
    // limit amplitude index to deData->t0->duration_raw
    double duration = deData->t0->duration_raw;
    // 20111230 AJL       if (deData->flag_complete_t0 && duration != T0_INVALID) {
    if (duration != T0_INVALID) {
#ifdef USE_CLOSE_S_FOR_DURATION // 20111222 TEST AJL - use S duration
        duration = deData->t0->duration_corr;
        if (duration < 0.0)
            printf("ERROR: calculate_mB_Mag: deData->t0->duration_corr not set: this should not happen!\n");
#endif
        if (1) { // 20111224 TEST AJL
            // limit mb integration window to S-P time
            double duration_before = duration;
            double ttime_P = deData->ttime_P;
            if (ttime_P == VALUE_NOT_SET)
                printf("ERROR: calculate_mB_Mag: deData->ttime_P not set: this should not happen! (%s_%s_%s_%s %s)\n",
                    deData->network, deData->station, deData->location, deData->channel, deData->phase);
            if (ttime_P > 0.0) {
                double ttime_S = deData->ttime_S;
                if (ttime_S > 0.0) {
                    if (duration_before > ttime_S - ttime_P)
                        duration = ttime_S - ttime_P;
                    //printf("DEBUG: limit Mwpd integration window to S-P time: dist=%g  S-P: %g  duration: %g->%g\n",
                    //        deData->epicentral_distance, ttime_S - ttime_P, duration_before, duration);
                }
            }
        } // END - 20111224 TEST AJL

        int index_max_test = (int) (0.5 + (duration + START_ANALYSIS_BEFORE_P_MB) / deData->deltaTime); // end magnitude calcaulation at T0 duration
        if (index_max_test < index_max)
            index_max = index_max_test;
    }
    double peak = 0.0;
    int n;
    //int ierr = 1;
    for (n = 0; n < index_max; n++) {
        if (deData->mb->amplitude[n] > peak) {
            peak = deData->mb->amplitude[n];
            /*if (ierr && peak > 1e10) {
                printf("ERROR: peak %g index %d  sampleRate %g\n", peak, n, 1.0 / deData->deltaTime);
                ierr = 0;
            }*/
        }
    }
    //printf("-----: peak %g index %d  sampleRate %g\n", peak, n, 1.0 / deData->deltaTime);
    double amp_value = peak / (2.0 * PI);

    //double Q_value = get_Q_Delta_PV_value(deData->epicentral_distance, hypo_depth);
    double Q_value = get_Q_Depth_Delta_PV_value(deData->epicentral_distance, hypo_depth);

    double mB_mag = MB_INVALID;
    if (amp_value > FLT_MIN && Q_value > FLT_MIN) {
        mB_mag = log10(amp_value) + Q_value;
        if (use_mb_correction) {
            mB_mag = calculate_corrected_mB_INGV_EE(mB_mag);
        }
    }
    deData->mb->mag = mB_mag;

    //printf("DEBUG: @@@@  %s_%s : mB_mag: %g = log10(amp_value %g)%g + Q_value %g (distance=%g, depth=%g) : peak=%g\n",
    //        deData->network, deData->station, mB_mag, amp_value, log10(amp_value), Q_value, deData->epicentral_distance, hypo_depth, peak);

    return (mB_mag);

}

/** calculate broadband first motion based on Mwp and Mwpd magnitude information
 *
 */
/*
int calculate_BRB_first_motion_polarity(TimedomainProcessingData* deData, double *pfmquality) {

    // default is pick polarity
    int fmpolarity = deData->pickData->polarity;

    // quality is function of signal to noise ratio
    double fmquality = 0.0;
    double snr_hf = 0.0;
    double snr_brb = 0.0;

    if (deData->pick_stream == STREAM_HF) {
        snr_hf = deData->a_ref < 0.0 || deData->sn_pick < FLT_MIN ? 0.0 : deData->a_ref / deData->sn_pick;
        fmquality = (snr_hf - SIGNAL_TO_NOISE_RATIO_HF_MIN) / SIGNAL_TO_NOISE_RATIO_HF_MIN;
    } else {
        snr_brb = deData->mb->pulse_first_polarity_amp < 0.0 || deData->sn_brb_bp_pick < FLT_MIN ? 0.0 : deData->mb->pulse_first_polarity_amp / deData->sn_brb_bp_pick;
        fmquality = (snr_brb - SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN) / SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN;
        // use mb polarity if s/n OK and mb brb polarity defined
        if (!deData->flag_snr_brb_too_low && deData->mb->pulse_first_polarity != 0) {
            // polarity is sign of first mb peak
            fmpolarity = deData->mb->pulse_first_polarity;
        }
    }
    if (fmquality < 0.0)
        fmquality = 0.0;
    else if (fmquality > 1.0)
        fmquality = 1.0;

    printf("%s_%s  snr_brb %g  fmpolarity %d  fmquality %g\n",
            deData->network, deData->station, snr_brb, fmpolarity, fmquality);

 *pfmquality = fmquality;
    return (fmpolarity);

    // polarity is sign of Mwp integral sum for first pulse
    if (deData->mwp->pulse_first > 0.0) {
        fmpolarity = 1;
    } else if (deData->mwp->pulse_first < 0.0) {
        fmpolarity = -1;
    } else {
        fmpolarity = 0;
    }

    printf("%s_%s  snr_brb %g  mwp->pulse_first_dur %g  pulse_first %g  fmpolarity %d  fmquality %g\n",
    deData->network, deData->station, snr_brb, deData->mwp->pulse_first_dur, deData->mwp->pulse_first, fmpolarity, fmquality);

 *pfmquality = fmquality;
    return (fmpolarity);

 }*/

