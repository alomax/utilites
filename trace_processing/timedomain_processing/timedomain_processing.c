/***************************************************************************
 * timedomain_processing.c:
 *
 * TODO: add doc
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2009.02.03
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/stat.h>

#include "../alomax_matrix/alomax_matrix.h"
#include "../alomax_matrix/polarization.h"
#include "../picker/PickData.h"
#include "timedomain_processing_data.h"
#include "ttimes.h"
#include "location.h"
#include "timedomain_processing_report.h"
#include "timedomain_processing.h"
#include "timedomain_memory.h"
#include "timedomain_processes.h"
#include "timedomain_filter.h"
#include "../picker/FilterPicker5_Memory.h"
#include "../picker/FilterPicker5.h"
#include "../response/response_lib.h"
#include "../miniseed_utils/miniseed_utils.h"
#include "../statistics/statistics.h"


// time domain processing memory
static timedomain_memory* pmem_filter_hf[MAX_NUM_SOURCES];
static timedomain_memory* pmem_smooth_a_ref[MAX_NUM_SOURCES];
static timedomain_memory* pmem_smooth_t50[MAX_NUM_SOURCES];
static timedomain_memory* pmem_smooth_brb_hp_rms[MAX_NUM_SOURCES];
static timedomain_memory* pmem_smooth_brb_hp_int_rms[MAX_NUM_SOURCES];
static timedomain_memory* pmem_smooth_wwssn_sp_rms[MAX_NUM_SOURCES];
static timedomain_memory* pmem_filter_hp[MAX_NUM_SOURCES]; // tauc_To, Mwp
static timedomain_memory* pmem_filter_hp_int[MAX_NUM_SOURCES]; // Mwp  // 20110401 AJL - Added
static timedomain_memory* pmem_filter_bp[MAX_NUM_SOURCES]; // mb
//static timedomain_memory* pmem_filter_lp[MAX_NUM_SOURCES]; // tauc_To
static timedomain_memory* pmem_inst_per[MAX_NUM_SOURCES]; // tauc_To
static timedomain_memory* pmem_inst_per_dval[MAX_NUM_SOURCES]; // tauc_To
//static timedomain_memory* pmem_mwp[MAX_NUM_SOURCES]; // mwp
// picker memory
static FilterPicker5_Memory* pmem_pick_hf[MAX_NUM_SOURCES];
static FilterPicker5_Memory* pmem_pick_raw[MAX_NUM_SOURCES];
// local t50 memory
static timedomain_memory* pmem_local_t50[MAX_NUM_SOURCES];
// local a_ref memory
static timedomain_memory* pmem_local_a_ref[MAX_NUM_SOURCES];
static double last_a_ref[MAX_NUM_SOURCES];
static int last_a_ref_ellapsed_index_count[MAX_NUM_SOURCES];
// local brb memory
static timedomain_memory* pmem_local_brb_hp[MAX_NUM_SOURCES]; // 20110401 AJL - Added pmem_local_brb_hp
static timedomain_memory* pmem_local_brb_hp_int_sqr[MAX_NUM_SOURCES]; // 20121115 AJL - Added pmem_local_brb_hp_int_sqr
static timedomain_memory* pmem_local_brb_hp_int_rms[MAX_NUM_SOURCES]; // 20110401 AJL - Added pmem_local_brb_hp_int_rms
static timedomain_memory* pmem_local_brb_hp_sqr[MAX_NUM_SOURCES];
static timedomain_memory* pmem_local_brb_hp_rms[MAX_NUM_SOURCES];
static timedomain_memory* pmem_local_wwssn_sp[MAX_NUM_SOURCES]; // 20110401 AJL - Added pmem_local_wwssn_sp
static timedomain_memory* pmem_local_wwssn_sp_sqr[MAX_NUM_SOURCES];
static timedomain_memory* pmem_local_wwssn_sp_rms[MAX_NUM_SOURCES];

static char pick_str[1024];

static char* geogfile_name;
static char* gainfile_name;
//static time_t geogfile_modtime = 0;

static int num_internet_gain_sources = 0;
#define RESPONSE_NOT_SET -999
#define MWPD_GAIN_FREQUENCY 0.15
static int num_internet_station_sources = 0;

static int num_pick_channels = -1;
static char pick_channel_list[MAX_NUM_PICK_CHANNEL_LIST][32];
static int num_process_orientations = -1;
static char process_orientation_list[MAX_NUM_PROCESS_ORIENTATION_LIST][32];
static PickParams rawPickParams[MAX_NUM_PICK_CHANNEL_LIST];
static PickParams hfPickParams[MAX_NUM_PICK_CHANNEL_LIST];



// 20110401 AJL - added MAX_ACCEPTABLE_PICK_ERROR so that maximum time of pick before begin of data packet can be limited
#define MAX_ACCEPTABLE_PICK_ERROR 5.0

#define HF_MINUS_RAW_PICK_WINDOW 30     // maximum time diffrence between hf and raw pick to use raw pick for location, otherwise, use hf pick

#define USE_TO_FOR_SIGNAL_LEVEL_WINDOW 1 // 20121115 TEST AJL - calculate signal level in window of width To after P, instead of fixed window width

//#define TEST_TO_FILE

//#define USE_RAW_PICKS_UNCONDITIONALLY  // # TEST! 20150930 AJL; seems to work well 20151119 AJL
// 20180823 AJL - IMPORTANT! Disabled above, likely caused false events due to picking strong, long period arrivals from preceding large events without verification of existence of a corresponding high-frequency pick.   Some if these false events had very large Mwp and Mwpd due to large amplitude and long period of picked arrivals.
/* In following, Tonga events are false and have very large Mwp and Mwpd
 1534638104751 	1	85	83	18.6	127	171	-2.10	1.3  	2018.08.19-00:21:46	-15.24	-175.29	   15	307	   28	B	0.0 	 0	35.4 	 1	0.0 	5.7 	 60	7.4 	 54	-1 	 0	7.8 	 22		 Tonga Islands
 1534638073375 	6	94	92	12.8	82	118	-2.74	2.9  	2018.08.19-00:21:13	-17.16	-173.37	   15	20	   12	B	1.6 	 1	22.3 	 2	0.0 	5.9 	 66	7.4 	 58	-1 	 1	7.9 	 22		 Tonga Islands
 1534637975904 	14	609	304	3.8	12	20	-1.49	1.1  	2018.08.19-00:19:39	-18.16	-178.00	   13	576	   14	A	0.4 	 164	4.6 	 30	2.0 	6.2 	 173	7.8 	 151	31 	 164	8.1 	 161		 Fiji Islands Region
 */



/**************************************************************************/
// general use functions
//

/***************************************************************************
 * slp_dtime_curr:
 *
 * Get the current time from the system as Unix/POSIX epoch time with double
 * precision.  This function has microsecond resolution.
 *
 * Return a double precision Unix/POSIX epoch time.
 *
 * 20151116 AJL - adapted from libslink/slplatform.c slp_dtime
 ***************************************************************************/
double
slp_dtime_curr(void) {
    struct timeval tv;
    gettimeofday(&tv, (struct timezone *) 0);
    return ((double) tv.tv_sec + ((double) tv.tv_usec / 1000000.0));
} /* End of slp_dtime() */


/**************************************************************************/
// timedomain processing functions
//

/***************************************************************************
 * calculate_signal_level:
 *
 * Utility function
 ***************************************************************************/
double calculate_signal_level(char *id_str, int source_id, double signal_level_window,
        float *data_float_sqr, timedomain_memory** pmem_local_sqr, int ioffset_pick, int numsamples, double deltaTime) {
    int ioffset_signal_level = ioffset_pick + (int) (0.5 + signal_level_window / deltaTime);
    // check for rounding error near numsamples
    if (ioffset_signal_level == numsamples + 1)
        ioffset_signal_level = numsamples;
    int n;
    double sqr_val_sum = 0.0;
    int icount = 0;
    for (n = ioffset_pick; n < ioffset_signal_level; n++) {
        if (n >= 0 && n < numsamples) { // value in current sqr data array
            sqr_val_sum += data_float_sqr[n];
        } else if (n >= -pmem_local_sqr[source_id]->numInput && n < 0) { // value in memory brb_hp_rms data array
            int index = n + pmem_local_sqr[source_id]->numInput;
            sqr_val_sum += pmem_local_sqr[source_id]->input[index];
        } else if (n >= numsamples) {
            printf("ERROR: calculate_signal_level %s: attempt to access data_float_sqr array of size=%d at invalid index=%d: this should not happen!\n",
                    id_str, numsamples, n);
            return (-1.0); // should never reach here
        } else {
            printf("ERROR: calculate_signal_level %s: attempt to access pmem_local_sqr array of size=%d at invalid index=%d: this should not happen!\n",
                    id_str, pmem_local_sqr[source_id]->numInput, n + pmem_local_sqr[source_id]->numInput);
            return (-1.0); // should never reach here
        }
        icount++;
    }

    double signal_level = 0.0;
    if (icount > 0)
        signal_level = sqrt(sqr_val_sum / (double) icount);

    return (signal_level);
}

/***************************************************************************
 * td_process_timedomain_processing:
 *
 * TODO: add doc.
 ***************************************************************************/
int td_process_timedomain_processing(MSRecord* pmsrecord, char* sladdr, int source_id, char* source_name, double deltaTime, float **pdata_float_in, int numsamples,
        char* network, char* station, char* location, char* channel, int year, int doy, int month, int mday, int hour, int min, double dsec,
        int verbose, int flag_clipped, double data_latency, char *calling_routine,
        int no_aref_level_check,
        int* pnum_new_loc_picks, double count_new_use_loc_picks_cutoff_time
        ) {

    float *data_float_raw = *pdata_float_in;

    // DEBUG
    //if (strcmp(station, "PMSA") == 0 || strcmp(station, "BRVK") == 0)
    //	verbose = 100;

    int return_value = 0;
    int n;


    if (source_id >= MAX_NUM_SOURCES) {
        if (verbose)
            printf("ERROR: source_id >= MAX_NUM_SOURCES: source_id=%d MAX_NUM_SOURCES=%d\n", source_id, MAX_NUM_SOURCES);
        return (-1);
    }

    // check if data contiguous
    // 20150211 AJL - moved this block here to head of function so check is performed
    struct tm tm_time;
    tm_time.tm_year = year - 1900;
    tm_time.tm_mon = month - 1;
    tm_time.tm_mday = mday;
    tm_time.tm_hour = hour;
    tm_time.tm_min = min;
    tm_time.tm_sec = (int) dsec;
    tm_time.tm_isdst = 0;
    tm_time.tm_gmtoff = 0;
    // 20121127 AJL - bug fix
    //double data_start_time = (double) mktime(&tm_time); // assumes mktime returns seconds
    double data_start_time = (double) timegm(&tm_time) + (dsec - (double) ((int) dsec)); // assumes timegm returns seconds
    int flag_non_contiguous = 0; // 20111230 AJL - added
    if (channelParameters[source_id].last_data_end_time > 0.0) {
        double tdiff = data_start_time - channelParameters[source_id].last_data_end_time;
        // 20150211 AJL - separated early and late case - early packets are now ignored in case next packet is contiguous
        if (tdiff < -1.5) { // 1 sec tolerance - early
            if (verbose > 1) {
                printf("Warning: %.2s_%.5s_%.2s_%.3s dt=%f: data packet start time = %lf not contiguous: %fs earlier than last data end time = %lf, ignoring packet\n",
                        network, station, location, channel, deltaTime, data_start_time, tdiff, channelParameters[source_id].last_data_end_time);
            }
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DATA_NON_CONTIGUOUS;
            channelParameters[source_id].count_non_contiguous++;
            return (-1);
        }
        if (tdiff > 1.5) { // 1 sec tolerance - late
            if (verbose > 1) {
                printf("Warning: %.2s_%.5s_%.2s_%.3s dt=%f: data packet start time = %lf not contiguous: %fs later than last data end time = %lf, assuming gap\n",
                        network, station, location, channel, deltaTime, data_start_time, tdiff, channelParameters[source_id].last_data_end_time);
            }
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DATA_NON_CONTIGUOUS;
            channelParameters[source_id].count_non_contiguous++;
            flag_non_contiguous = 1;
        }
    }
    channelParameters[source_id].last_data_end_time = data_start_time + (double) numsamples * deltaTime;

    // waveform export
    if (waveform_export_enable) {
        return_value = mslist_addMsRecord(pmsrecord, &(waveform_export_miniseed_list[source_id]), &(num_waveform_export_miniseed_list[source_id]));
        if (return_value > 0) { // record added correctly
            hptime_t remove_before_time = pmsrecord->starttime - ((long int) waveform_export_memory_sliding_window_length) * HPTMODULUS;
            //printf("DEBUG: remove_before_time=%ld, pmsrecord->starttime=%ld, diff=%ld, waveform_export_memory_sliding_window_length=%ld\n",
            //remove_before_time / HPTMODULUS, pmsrecord->starttime / HPTMODULUS, (pmsrecord->starttime - remove_before_time) / HPTMODULUS, waveform_export_memory_sliding_window_length);
            return_value = mslist_removeMsRecords(remove_before_time, &(waveform_export_miniseed_list[source_id]), &(num_waveform_export_miniseed_list[source_id]));
            if (verbose > 2) {
                char mdtimestr1[128];
                char mdtimestr2[128];
                hptime_t hptime_start = mslist_getStartTime(waveform_export_miniseed_list[source_id], num_waveform_export_miniseed_list[source_id]);
                hptime_t hptime_end = mslist_getEndTime(waveform_export_miniseed_list[source_id], num_waveform_export_miniseed_list[source_id]);
                printf("INFO: %.2s_%.5s_%.2s_%.3s: waveform export memory: start=%s, end=%s, diff=%lds\n", network, station, location, channel,
                        ms_hptime2mdtimestr(hptime_start, mdtimestr1, 1), ms_hptime2mdtimestr(hptime_end, mdtimestr2, 1), (long int) ((hptime_end - hptime_start) / HPTMODULUS));
            }
        }
    }

    // set index of Internet timeseries query structure for this sladdr
    int n_int_tseries = -1;
    if (sladdr != NULL) {
        for (n = 0; n < num_internet_timeseries; n++) {
            if (strcmp(internetTimeseriesQueryParams[n].sladdr, sladdr) == 0) {
                n_int_tseries = n;
                break;
            }
        }
    } else {
        if (num_internet_timeseries > 0)
            n_int_tseries = 0;
    }

    // check have station coordinates set station as active
    //printf("DEBUG: Check ws station: %s %s  time(NULL) %ld >? channelParameters[source_id].internet_station_query_checked_time %ld (diff=%f)\n",
    //        network, station, time(NULL), channelParameters[source_id].internet_station_query_checked_time,
    //        difftime(time(NULL), channelParameters[source_id].internet_station_query_checked_time) / (24 * 3600));
    int icheck_ws_station_coords = difftime(time(NULL), channelParameters[source_id].internet_station_query_checked_time) > (WEB_SERIVCE_METADATA_CHECK_INTERVAL * 24 * 3600);
    if (!channelParameters[source_id].have_coords || icheck_ws_station_coords) {
        td_set_station_coordinates(source_id, network, station, location, channel, year, month, mday, verbose, icheck_ws_station_coords);
        if (!channelParameters[source_id].have_coords) {
            if (verbose)
                printf("ERROR: sensor coordinates not found: %s %s, please add sensor data to file: %s\n", network, station, geogfile_name);
            return (-1);
        }
        // set misc station parameters
        channelParameters[source_id].deltaTime = deltaTime;
        channelParameters[source_id].n_int_tseries = n_int_tseries;
        // find and store references to other channel orientations for this net/sta/chan
        associate3CompChannelSet(channelParameters, num_sources_total, source_id);
    }

    double lat = channelParameters[source_id].lat;
    double lon = channelParameters[source_id].lon;
    double elev = channelParameters[source_id].elev;
    double azimuth = channelParameters[source_id].azimuth;
    double dip = channelParameters[source_id].dip;
    double station_quality_weight = channelParameters[source_id].qualityWeight;
    channelParameters[source_id].staActiveInReportInterval = 1;
    channelParameters[source_id].data_latency = data_latency;

    if (use_station_corrections) {

        // check for station corrections
        int icheck_sta_corr = 0; // TODO: maybe add periodic check for updated station residuals, not necessary if add automatic update of residuals for each event located
        if (!channelParameters[source_id].sta_corr_checked || icheck_sta_corr) {
            td_set_station_corrections(source_id, network, station, location, channel, year, month, mday, verbose, icheck_sta_corr);
        }
    }


    // check have channel gain
    // get response/gain from pole zero source
    int flag_have_gain = 0;
    double gain_meters = 1.0;
    double gain_microns = 1.0; // microns
    //20140223 AJL  if (flag_do_mwp || flag_do_mwpd || flag_do_mb) {
    int icheck_ws_gain = difftime(time(NULL), chan_resp[source_id].internet_gain_query_checked_time) > (WEB_SERIVCE_METADATA_CHECK_INTERVAL * 24 * 3600);
    if (!chan_resp[source_id].have_gain || icheck_ws_gain) {
        td_set_channel_gain(source_id, network, station, location, channel, year, doy, verbose, icheck_ws_gain);
    }
    if (chan_resp[source_id].have_gain && chan_resp[source_id].responseType == DERIVATIVE_TYPE) { // velocity data
        flag_have_gain = 1;
        gain_meters = chan_resp[source_id].gain;
        gain_microns = chan_resp[source_id].gain / 1.0e6; // microns
    }
    //20140223 AJL  }



    // if channel orientation not to be picked and processed, return
    if (channelParameters[source_id].process_this_channel_orientation == 0) {
        // no further processing for this channel, cleanup not needed since no filtering done
        return (0);
    }
    if (channelParameters[source_id].process_this_channel_orientation == -1) {
        // check if channel is in process_orientation_list, if list exists
        char orientation[2];
        strcpy(orientation, channel + 2);
        if (num_process_orientations > 0) {
            for (n = 0; n < num_process_orientations; n++) {
                //printf("DEBUG: n=%d/%d, orientation=%s, in: <%s> \n", n, num_process_orientations, orientation, process_orientation_list[n]);
                if (strstr(process_orientation_list[n], orientation) != NULL) {
                    // orientation is in list, use for pick, processing
                    channelParameters[source_id].process_this_channel_orientation = 1;
                    break;
                }
            }
            if (n >= num_process_orientations) {
                // orientation is not in list, do not use for pick, processing
                //printf("DEBUG: %.2s_%.5s_%.2s_%.3s orient=%s: SKIPPED for pick/process\n", network, station, location, channel, orientation);
                channelParameters[source_id].process_this_channel_orientation = 0;
                // no further processing for this channel, cleanup not needed since no filtering done
                return (0);
            }
            //printf("DEBUG: %.2s_%.5s_%.2s_%.3s orient=%s: USED for pick/process\n", network, station, location, channel, orientation);
        }
    }

    //
    if (flag_clipped || flag_non_contiguous) { // 20111230 AJL - added
        // free all memory for this source to re-initialize time domain processing and avoid false picks
        td_process_free_timedomain_memory(source_id);
    }


    // check if sample rate changed, this will happen if active channel is replaced by channel with different, preferred loc
    /* NOTE: this clean up does not yet work, simpler to leave old active memory remaining until clean-up
     *   see packet_handler() "check if duplicate has preferred location", lines:
     *      // ensures that new source is created, old active memory will remain until clean-up
     *      source_id = num_sources_total;
     *
    if (pmem_filter_hf[source_id] != NULL && pmem_filter_hf[source_id]->deltaTime != deltaTime) {
        if (verbose)
            printf("\nInfo: Existing deltaTime = %f not equal to current deltaTime = %f, OK if active loc replaced by preferred loc.\n",
                pmem_filter_hf[source_id]->deltaTime, deltaTime);
        // clean up memory, this allows change in sample rate
        td_process_free_timedomain_memory(source_id);
        // clear any existing deData for this source
        for (n = num_de_data - 1; n >= 0; n--) {
            TimedomainProcessingData* deData = data_list[n];
            if (deData->source_id == source_id) {
                removeTimedomainProcessingDataFromDataList(deData, &data_list, &num_de_data);
                free_TimedomainProcessingData(deData);
            }
        }
    }
     */

    // declare filtered data arrays before any goto statement, prevents "may be used uninitialized in this function" warnings for these array pointers
    float* data_float_hf = NULL;
    float *data_float_t50 = NULL;
    float *data_float_a_ref = NULL;
    float *data_float_tauc = NULL;
    float *data_float_brb_hp = NULL;
    float *data_float_brb_hp_sqr = NULL;
    float *data_float_brb_hp_rms = NULL;
    float *data_float_brb_hp_int_sqr = NULL;
    float *data_float_brb_hp_int_rms = NULL; // 20120512 AJL - added
    float *data_float_wwssn_sp = NULL;
    float *data_float_wwssn_sp_sqr = NULL;
    float *data_float_wwssn_sp_rms = NULL;


    // filter
    if (measuresEnable || pickOnHFStream) {
        data_float_hf = copyData(data_float_raw, numsamples);
        if (filter_bp_bu_co_1_5_n_4(deltaTime, data_float_hf, numsamples, &pmem_filter_hf[source_id], 1) == NULL) {
            if (verbose)
                printf("ERROR: deltaTime not supported by filter_bp_bu_co_1_5_n_4: deltaTime=%f\n", deltaTime);
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DT_NOT_SUPPORTED_BY_FILTER;
            td_process_free_timedomain_memory(source_id);
            return_value = -1;
            goto cleanup;
            //return (-1); // deltaTime not supported
        }
    }

#ifdef TEST_TO_FILE
    char pickfile[1024];
    sprintf(pickfile, "%s.picks", calling_routine);
    FILE* fp_picks = fopen(pickfile, "a");
#endif

    // set tUpEvent for brb and hf
    double tUpEvent_hf = 4.0;
    double tUpEvent_brb = 2.0;

    // set channel index for pick parameters
    // check if channel is in pick_channel_list, if list exists
    int n_pick_channel = 0;
    if (num_pick_channels > 0) {
        for (n = 0; n < num_pick_channels; n++) {
            //printf("DEBUG: n=%d/%d, channel=%s, in: <%s> \n", n, num_pick_channels, channel, pick_channel_list[n]);
            if (*pick_channel_list[n] == '\0' || strstr(pick_channel_list[n], channel) != NULL) {
                break; // no more channels in list or data channel is specified for this index
            }
        }
        if (n < num_pick_channels)
            n_pick_channel = n;
    }
    //printf("DEBUG: n_pick_channel=%d/%d, channel=%s, in: <%s> \n", n_pick_channel, num_pick_channels, channel, pick_channel_list[n_pick_channel]);


    if (pickOnRawStream) {
        tUpEvent_brb = rawPickParams[n_pick_channel].tUpEvent;
        if (verbose > 2)
            printf("picker_func_test: filp_test filtw %f ltw %f thres1 %f thres2 %f tupevt %f res PICKS\n",
                rawPickParams[n_pick_channel].filterWindow, rawPickParams[n_pick_channel].longTermWindow, rawPickParams[n_pick_channel].threshold1, rawPickParams[n_pick_channel].threshold2, rawPickParams[n_pick_channel].tUpEvent);
        PickData** ppick_list = NULL;
        int num_picks = 0;
        if (verbose > 99)
            printf("%s (%d)  &pmem_pick_raw[source_id] %ld  last_a_ref[source_id] %f\n", station, source_id, (long) & pmem_pick_raw[source_id], last_a_ref[source_id]);
        Pick(deltaTime, data_float_raw, numsamples, rawPickParams[n_pick_channel].filterWindow, rawPickParams[n_pick_channel].longTermWindow,
                rawPickParams[n_pick_channel].threshold1, rawPickParams[n_pick_channel].threshold2, rawPickParams[n_pick_channel].tUpEvent,
                &pmem_pick_raw[source_id], 1, &ppick_list, &num_picks, source_name);
        // save and display picks
        char phase[16];
        int npick;
        for (npick = 0; npick < num_picks; npick++) {
            PickData* ppick = *(ppick_list + npick);
            double dindex = (ppick->indices[0] + ppick->indices[1]) / 2.0;
            // check basic pick quality
            double pick_error = 0.5 * deltaTime * fabs(ppick->indices[1] - ppick->indices[0]);
            if (pick_error > MAX_ACCEPTABLE_PICK_ERROR) {
                if (verbose)
                    printf("Info: new pick not added to deData list, pick error=%g > MAX_ACCEPTABLE_PICK_ERROR=%g: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f\n",
                        pick_error, MAX_ACCEPTABLE_PICK_ERROR, source_name, year, month, mday, hour, min, dsec + dindex * deltaTime);
                continue;
            }
            // save picks to list
            int init_ellapsed_index_count = (int) (0.5 + dindex); // count at beginning of current record
            TimedomainProcessingData* deData = new_TimedomainProcessingData(sladdr, n_int_tseries, source_id, station, location, channel, network, deltaTime, lat, lon, elev, azimuth, dip, station_quality_weight,
                    ppick, STREAM_RAW, init_ellapsed_index_count, year, month, mday, hour, min, dsec + dindex * deltaTime, flag_do_mwpd, waveform_export_enable);
            double time_dec_sec = (double) deData->t_time_t + deData->t_decsec;
            if (!pickOnHFStream) {
                deData->use_for_location = 1;
                if (time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                    (*pnum_new_loc_picks)++;
                }
            }
#ifdef USE_RAW_PICKS_UNCONDITIONALLY
            else {
                deData->use_for_location = 1;
                if (time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                    (*pnum_new_loc_picks)++;
                }
            }
#endif
            int num_added = addTimedomainProcessingDataToDataList(deData, &data_list, &num_de_data, 1, 1);
            if (num_added < 1) {
                if (verbose)
                    printf("Warning: Raw pick: new pick not added to deData list, nearby pick already exists in list: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f (%ld %f)\n",
                        source_name, deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec, deData->t_time_t, deData->t_decsec);
                free_TimedomainProcessingData(deData);
                continue;
            }
            // print pick
            sprintf(phase, "P%d_", npick);
            if (verbose > 1) {
                printf("Pick: %d", npick);
                printf(" %d", ppick->polarity);
                printf(" %f", ppick->polarityWeight);
                printf(" %f", ppick->indices[0]);
                printf(" %f", ppick->indices[1]);
                printf(" %f", ppick->amplitude);
                printf(" %d", ppick->amplitudeUnits);
                printf(" %f\n", ppick->period);
            }

            if (pickOnHFStream) {
                // 20121127 AJL - added in case later HF pick finishes picking before this raw pick
                // if using HF picks, then look for HF picks just after this raw pick; if present, flag this raw pick as use for locate; otherwise, flag HF as use for locate
                {
                    // check for timedomain-processing data for this source_id
                    double diff;
                    double time_dec_sec_hf;
                    int ndata;
                    for (ndata = 0; ndata < num_de_data; ndata++) {
                        TimedomainProcessingData* deData_hf = data_list[ndata];
                        if (deData_hf->source_id != source_id || deData_hf->pick_stream != STREAM_HF
                                || deData_hf->use_for_location_twin_data != NULL
                                )
                            continue;
                        // check if hf pick is within required window after raw pick
                        time_dec_sec_hf = (double) deData_hf->t_time_t + deData_hf->t_decsec;
                        diff = time_dec_sec_hf - time_dec_sec;
                        if (diff > 0.0 && diff <= HF_MINUS_RAW_PICK_WINDOW) { // if diff == 0, use HF since sorts earlier in pick list
                            // do not use HF
                            if (deData_hf->use_for_location && time_dec_sec_hf >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)--;
                            }
                            deData_hf->use_for_location = 0;
                            if (!deData->use_for_location && time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)++;
                            }
                            deData->use_for_location = 1;
                            deData_hf->use_for_location_twin_data = deData;
                            deData->use_for_location_twin_data = deData_hf;
                            break; // 20150423 AJL - bug fix, prevent double twinning and loosing pointer
                        }
#ifdef USE_RAW_PICKS_UNCONDITIONALLY
                        if (diff <= 0.0 && diff >= -HF_MINUS_RAW_PICK_WINDOW) { // if diff <= 0, use HF since sorts earlier in pick list
                            // do not use RAW
                            if (!deData_hf->use_for_location && time_dec_sec_hf >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)++;
                            }
                            deData_hf->use_for_location = 1;
                            if (deData->use_for_location && time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)--;
                            }
                            deData->use_for_location = 0;
                            deData_hf->use_for_location_twin_data = deData;
                            deData->use_for_location_twin_data = deData_hf;
                            break; // 20150423 AJL - bug fix, prevent double twinning and loosing pointer
                        }
#endif
                    }
                }
            }

#ifdef TEST_TO_FILE
            // write pick to <pick_file> in NLLOC_OBS format
            printNlloc(pick_str, ppick, deltaTime, station, location, channel, "?", phase, year, month, mday, hour, min, dsec + dindex * deltaTime);
            fprintf(fp_picks, "%s\n", pick_str);
#endif
            if (verbose > 1)
                printf("%s\n", pick_str);
        }
        free_PickList(ppick_list, num_picks);
    }

    if (pickOnHFStream) {
        // pick
        tUpEvent_hf = hfPickParams[n_pick_channel].tUpEvent;
        if (verbose > 2)
            printf("picker_func_test: filp_test filtw %f ltw %f thres1 %f thres2 %f tupevt %f res PICKS\n",
                hfPickParams[n_pick_channel].filterWindow, hfPickParams[n_pick_channel].longTermWindow, hfPickParams[n_pick_channel].threshold1, hfPickParams[n_pick_channel].threshold2, hfPickParams[n_pick_channel].tUpEvent);
        PickData** ppick_list = NULL;
        int num_picks = 0;
        if (verbose > 99)
            printf("%s (%d)  &pmem_pick_hf[source_id] %ld  last_a_ref[source_id] %f\n", station, source_id, (long) & pmem_pick_hf[source_id], last_a_ref[source_id]);
        Pick(deltaTime, data_float_hf, numsamples, hfPickParams[n_pick_channel].filterWindow, hfPickParams[n_pick_channel].longTermWindow,
                hfPickParams[n_pick_channel].threshold1, hfPickParams[n_pick_channel].threshold2, hfPickParams[n_pick_channel].tUpEvent,
                &pmem_pick_hf[source_id], 1, &ppick_list, &num_picks, source_name);
        // save and display picks
        char phase[16];
        int npick;
        for (npick = 0; npick < num_picks; npick++) {
            PickData* ppick = *(ppick_list + npick);
            double dindex = (ppick->indices[0] + ppick->indices[1]) / 2.0;
            // check basic pick quality
            double pick_error = 0.5 * deltaTime * fabs(ppick->indices[1] - ppick->indices[0]);
            if (pick_error > MAX_ACCEPTABLE_PICK_ERROR) {
                if (verbose)
                    printf("Info: new pick not added to deData list, pick error=%g > MAX_ACCEPTABLE_PICK_ERROR=%g: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f\n",
                        pick_error, MAX_ACCEPTABLE_PICK_ERROR, source_name, year, month, mday, hour, min, dsec + dindex * deltaTime);
                continue;
            }
            // save picks to list
            int init_ellapsed_index_count = (int) (0.5 + dindex); // count at beginning of current record
            TimedomainProcessingData* deData = new_TimedomainProcessingData(sladdr, n_int_tseries, source_id, station, location, channel, network, deltaTime, lat, lon, elev, azimuth, dip, station_quality_weight, ppick, STREAM_HF, init_ellapsed_index_count,
                    year, month, mday, hour, min, dsec + dindex * deltaTime, flag_do_mwpd, waveform_export_enable);
            double time_dec_sec = (double) deData->t_time_t + deData->t_decsec;
            //if (!pickOnRawStream) {
            deData->use_for_location = 1;
            if (time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                (*pnum_new_loc_picks)++;
            }
            // }
            int num_added = addTimedomainProcessingDataToDataList(deData, &data_list, &num_de_data, 1, 1);
            if (num_added < 1) {
                if (verbose)
                    printf("Warning: HF pick: new pick not added to deData list, nearby pick already exists in list: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f\n",
                        source_name, deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec);
                free_TimedomainProcessingData(deData);
                continue;
            }
            // print pick
            sprintf(phase, "P%d_", npick);
            if (verbose > 1) {
                printf("Pick: %d", npick);
                printf(" %d", ppick->polarity);
                printf(" %f", ppick->polarityWeight);
                printf(" %f", ppick->indices[0]);
                printf(" %f", ppick->indices[1]);
                printf(" %f", ppick->amplitude);
                printf(" %d", ppick->amplitudeUnits);
                printf(" %f\n", ppick->period);
            }

            if (pickOnRawStream) {
                // if using raw picks, then look for raw picks just before this HF pick; if present, flag as use for locate; otherwise, flag this HF as use for locate
                {
                    // check for timedomain-processing data for this source_id
                    double diff;
                    double time_dec_sec_raw;
                    int ndata;
                    //for (ndata = 0; ndata < num_de_data; ndata++) { // 20150423 AJL - bug fix, should search backwards from this pick
                    for (ndata = num_de_data - 1; ndata >= 0; ndata--) {
                        TimedomainProcessingData* deData_raw = data_list[ndata];
                        if (deData_raw->source_id != source_id || deData_raw->pick_stream != STREAM_RAW
                                || deData_raw->use_for_location_twin_data != NULL
                                )
                            continue;
                        // check if raw pick is within required window before hf pick
                        time_dec_sec_raw = (double) deData_raw->t_time_t + deData_raw->t_decsec;
                        diff = time_dec_sec - time_dec_sec_raw; // 20101116 AJL - bug fix, before did not check decimal seconds, could allow earlier HF pick to be ignored
                        if (diff >= 0.0 && diff <= HF_MINUS_RAW_PICK_WINDOW) { // if diff == 0, use raw since sorts earlier in pick list
                            // do not use HF
                            //if (diff >= -deltaTime && diff <= HF_MINUS_RAW_PICK_WINDOW) {
                            if (deData->use_for_location && time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)--;
                            }
                            deData->use_for_location = 0;
                            if (!deData_raw->use_for_location && time_dec_sec_raw >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)++;
                            }
                            deData_raw->use_for_location = 1;
                            deData->use_for_location_twin_data = deData_raw;
                            deData_raw->use_for_location_twin_data = deData;
                            break; // 20150423 AJL - bug fix, prevent double twinning and loosing pointer
                        }
#ifdef USE_RAW_PICKS_UNCONDITIONALLY
                        if (diff <= 0.0 && diff >= -HF_MINUS_RAW_PICK_WINDOW) { // if diff <= 0, use HF since sorts earlier in pick list
                            // do not use RAW
                            if (!deData->use_for_location && time_dec_sec >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)++;
                            }
                            deData->use_for_location = 1;
                            if (deData_raw->use_for_location && time_dec_sec_raw >= count_new_use_loc_picks_cutoff_time) {
                                (*pnum_new_loc_picks)--;
                            }
                            deData_raw->use_for_location = 0;
                            // 20160803 AJL - bug fix  deData->use_for_location_twin_data = deData;
                            deData->use_for_location_twin_data = deData_raw;
                            deData_raw->use_for_location_twin_data = deData;
                            break; // 20150423 AJL - bug fix, prevent double twinning and loosing pointer
                        }
#endif
                    }
                }
            }

#ifdef TEST_TO_FILE
            // write pick to <pick_file> in NLLOC_OBS format
            printNlloc(pick_str, ppick, deltaTime, station, location, channel, "?", phase, year, month, mday, hour, min, dsec + dindex * deltaTime);
            fprintf(fp_picks, "%s\n", pick_str);
#endif
            if (verbose > 1)
                printf("%s\n", pick_str);
        }

        free_PickList(ppick_list, num_picks);
    }

#ifdef TEST_TO_FILE
    fclose(fp_picks);
#endif


    //
    if (!measuresEnable) {
        return_value = 0;
        goto cleanup;
    }


    // timedomain-processing time series
    // square value
    applySqr(data_float_hf, numsamples);
    // smooth for T50
    data_float_t50 = copyData(data_float_hf, numsamples);
    applyBoxcarSmoothing(deltaTime, data_float_t50, numsamples, SMOOTHING_WINDOW_HALF_WIDTH_T50, &(pmem_smooth_t50[source_id]), 1);
    // perpare to save t50 to local memory, save enough samples to have t50 values back tUpEvent_hf seconds
    if (pmem_local_t50[source_id] == NULL) {
        double memory_time = tUpEvent_hf + MAX_ACCEPTABLE_PICK_ERROR;
        int num_local_t50_memory = 1 + (int) (memory_time / deltaTime);
        //printf("@@@@@@@@DEBUG: num_local_t50_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_t50_memory, memory_time, deltaTime);
        pmem_local_t50[source_id] = init_timedomain_memory(num_local_t50_memory, -1.0, 0, 0.0, deltaTime);
    }
    // smooth for a_ref
    data_float_a_ref = copyData(data_float_hf, numsamples);
    applyBoxcarSmoothing(deltaTime, data_float_a_ref, numsamples, SMOOTHING_WINDOW_HALF_WIDTH_A_REF, &(pmem_smooth_a_ref[source_id]), 1);
    // perpare to save a_ref to local memory, save enough samples to have a_ref values back tUpEvent and TIME_DELAY_A_REF seconds
    if (pmem_local_a_ref[source_id] == NULL) {
        // 20110401 AJL double memory_time = TIME_DELAY_A_REF > tUpEvent_hf ? TIME_DELAY_A_REF : tUpEvent_hf;
        double memory_time = tUpEvent_hf + MAX_ACCEPTABLE_PICK_ERROR;
        int num_local_a_ref_memory = 1 + (int) (memory_time / deltaTime);
        //printf("@@@@@@@@DEBUG: num_local_a_ref_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_a_ref_memory, memory_time, deltaTime);
        pmem_local_a_ref[source_id] = init_timedomain_memory(num_local_a_ref_memory, -1.0, 0, 0.0, deltaTime);
    }

    // brb processing time series
    if (flag_do_tauc || flag_do_mwp || flag_do_mwpd || flag_do_grd_vel) {
        data_float_brb_hp = copyData(data_float_raw, numsamples);
        // high-pass to remove mean, trend, lp noise
        if (filter_hp_bu_co_0_005_n_2(deltaTime, data_float_brb_hp, numsamples, &pmem_filter_hp[source_id], 1) == NULL) {
            if (verbose)
                printf("ERROR: deltaTime not supported by filter_hp_bu_co_0_005_n_2: deltaTime=%f\n", deltaTime);
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DT_NOT_SUPPORTED_BY_FILTER;
            td_process_free_timedomain_memory(source_id);
            return_value = -1;
            goto cleanup;
            //return (-1); // deltaTime not supported
        }
        // data for BRB S/N
        data_float_brb_hp_sqr = copyData(data_float_brb_hp, numsamples);
        // square value
        applySqr(data_float_brb_hp_sqr, numsamples);
        // data for BRB S/N
        data_float_brb_hp_rms = copyData(data_float_brb_hp_sqr, numsamples);
        // smooth for brb s/n
        applyBoxcarSmoothing(deltaTime, data_float_brb_hp_rms, numsamples, SMOOTHING_WINDOW_HALF_WIDTH_BRB_HP_RMS, &(pmem_smooth_brb_hp_rms[source_id]), 1);
        //int ntest = 10; printf("n=%d, data_float_raw=%f, data_float_brb_hp=%f, data_float_brb_hp_rms=%f\n", ntest, data_float_raw[ntest], data_float_brb_hp[ntest], data_float_brb_hp_rms[ntest]);
        // perpare to save brb_hp_rms to local memory, save enough samples to have brb_hp_rms values back tUpEvent and TIME_DELAY_A_REF seconds
        int num_local_brb_hp_sqr_memory;
        if (pmem_local_brb_hp_sqr[source_id] == NULL) {
            double memory_time = MAX_MWPD_DUR;
            num_local_brb_hp_sqr_memory = 1 + (int) (memory_time / deltaTime);
            //printf("@@@@@@@@DEBUG: num_local_brb_hp_sqr_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_brb_hp_sqr_memory, memory_time, deltaTime);
            pmem_local_brb_hp_sqr[source_id] = init_timedomain_memory(num_local_brb_hp_sqr_memory, -1.0, 0, 0.0, deltaTime);
        } else {
            num_local_brb_hp_sqr_memory = pmem_local_brb_hp_sqr[source_id]->numInput;
        }
        int num_local_brb_hp_rms_memory;
        if (pmem_local_brb_hp_rms[source_id] == NULL) {
            // 20110401 AJL double memory_time = TIME_DELAY_BRB_HP_RMS > tUpEvent_brb ? TIME_DELAY_BRB_HP_RMS : tUpEvent_brb;
            //memory_time += SHIFT_BEFORE_P_BRB_HP_SN;
            double memory_time = SHIFT_BEFORE_P_BRB_HP_SN + tUpEvent_brb + MAX_ACCEPTABLE_PICK_ERROR;
            num_local_brb_hp_rms_memory = 1 + (int) (memory_time / deltaTime);
            //printf("@@@@@@@@DEBUG: num_local_brb_hp_rms_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_brb_hp_rms_memory, memory_time, deltaTime);
            pmem_local_brb_hp_rms[source_id] = init_timedomain_memory(num_local_brb_hp_rms_memory, -1.0, 0, 0.0, deltaTime);
        } else {
            num_local_brb_hp_rms_memory = pmem_local_brb_hp_rms[source_id]->numInput;
        }
        if (flag_do_mwp || flag_do_mwpd) {
            // 20110401 AJL - Added pmem_local_brb_hp
            // prepare to save brb_hp to local memory, save enough samples to have brb_hp values back tUpEvent plus START_ANALYSIS_BEFORE_P_BRB_HP seconds
            int num_local_brb_hp_int_memory;
            if (pmem_local_brb_hp[source_id] == NULL) {
                double memory_time = START_ANALYSIS_BEFORE_P_BRB_HP + tUpEvent_brb + MAX_ACCEPTABLE_PICK_ERROR;
                num_local_brb_hp_int_memory = 1 + (int) (memory_time / deltaTime);
                //printf("@@@@@@@@DEBUG: num_local_brb_hp_int_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_brb_hp_int_memory, memory_time, deltaTime);
                pmem_local_brb_hp[source_id] = init_timedomain_memory(num_local_brb_hp_int_memory, -1.0, 0, 0.0, deltaTime);
            } else {
                num_local_brb_hp_int_memory = pmem_local_brb_hp[source_id]->numInput;
            }
            // 20120512 AJL - added data_float_brb_hp_int_rms, integrated to displacement for rms
            // data for BRB integrated (BRB displacement) rms
            // data for BRB S/N
            data_float_brb_hp_int_sqr = copyData(data_float_brb_hp, numsamples);
            // integrate value
            applyIntegral(deltaTime, data_float_brb_hp_int_sqr, numsamples, &(pmem_filter_hp_int[source_id]), 1);
            // square value
            applySqr(data_float_brb_hp_int_sqr, numsamples);
            data_float_brb_hp_int_rms = copyData(data_float_brb_hp_int_sqr, numsamples);
            // smooth for brb s/n
            applyBoxcarSmoothing(deltaTime, data_float_brb_hp_int_rms, numsamples, SMOOTHING_WINDOW_HALF_WIDTH_BRB_HP_INT_RMS, &(pmem_smooth_brb_hp_int_rms[source_id]), 1);
            //int ntest = 10; printf("n=%d, data_float_raw=%f, data_float_brb_hp=%f, data_float_brb_hp_int_rms=%f\n", ntest, data_float_raw[ntest], data_float_brb_hp[ntest], data_float_brb_hp_int_rms[ntest]);
            // perpare to save brb_hp_int_rms to local memory, save enough samples to have brb_hp_int_rms values back tUpEvent and TIME_DELAY_A_REF seconds
            int num_local_brb_hp_int_sqr_memory;
            if (pmem_local_brb_hp_int_sqr[source_id] == NULL) {
                double memory_time = MAX_MWPD_DUR;
                num_local_brb_hp_int_sqr_memory = 1 + (int) (memory_time / deltaTime);
                //printf("@@@@@@@@DEBUG: num_local_brb_hp_int_sqr_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_brb_hp_int_sqr_memory, memory_time, deltaTime);
                pmem_local_brb_hp_int_sqr[source_id] = init_timedomain_memory(num_local_brb_hp_int_sqr_memory, -1.0, 0, 0.0, deltaTime);
            } else {
                num_local_brb_hp_int_sqr_memory = pmem_local_brb_hp_int_sqr[source_id]->numInput;
            }
            int num_local_brb_hp_int_rms_memory;
            if (pmem_local_brb_hp_int_rms[source_id] == NULL) {
                double memory_time = SHIFT_BEFORE_P_BRB_HP_INT_SN + tUpEvent_brb + MAX_ACCEPTABLE_PICK_ERROR;
                num_local_brb_hp_int_rms_memory = 1 + (int) (memory_time / deltaTime);
                //printf("@@@@@@@@DEBUG: num_local_brb_hp_int_rms_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_brb_hp_int_rms_memory, memory_time, deltaTime);
                pmem_local_brb_hp_int_rms[source_id] = init_timedomain_memory(num_local_brb_hp_int_rms_memory, -1.0, 0, 0.0, deltaTime);
            } else {
                num_local_brb_hp_int_rms_memory = pmem_local_brb_hp_int_rms[source_id]->numInput;
            }
        }
    }
    if (flag_do_mb) {
        data_float_wwssn_sp = copyData(data_float_raw, numsamples);
        float* ifilter_stat = NULL;
        // filter to model WWSSN-SP data
        if (MB_MODE == MB_MODE_mb_DISP || MB_MODE == MB_MODE_mB) // transfer function is WWSSN_SP_DISP, gives velocity->velocity
        {
            ifilter_stat = filter_WWSSN_SP_DISP(deltaTime, data_float_wwssn_sp, numsamples, &pmem_filter_bp[source_id], 1);
        } else if (MB_MODE == MB_MODE_mb_VEL) // transfer function is WWSSN_SP_VEL, gives velocity->disp
        {
            ifilter_stat = filter_WWSSN_SP_VEL(deltaTime, data_float_wwssn_sp, numsamples, &pmem_filter_bp[source_id], 1);
        }
        //if (filter_bp_bu_co_0_333_5_n_4(deltaTime, data_float_wwssn_sp, numsamples, &pmem_filter_bp[source_id], 1) == NULL) {
        // TEST HF filter to better mimick WWSSN SP
        //if (filter_bp_bu_co_1_5_n_4(deltaTime, data_float_wwssn_sp, numsamples, &pmem_filter_bp[source_id], 1) == NULL) {
        if (ifilter_stat == NULL) {
            if (verbose)
                printf("ERROR: deltaTime not supported by filter_bp_bu_co_0_333_5_n_4: deltaTime=%f\n", deltaTime);
            //printf("ERROR: deltaTime not supported by filter_bp_bu_co_1_5_n_4: deltaTime=%f\n", deltaTime);
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DT_NOT_SUPPORTED_BY_FILTER;
            return_value = -1;
            td_process_free_timedomain_memory(source_id);
            goto cleanup;
            //return (-1); // deltaTime not supported
        }
        // 20110401 AJL - Added pmem_local_wwssn_sp
        // perpare to save wwssn_sp to local memory, save enough samples to have wwssn_sp values back tUpEvent and TIME_DELAY_A_REF seconds
        int num_local_wwssn_sp_memory;
        if (pmem_local_wwssn_sp[source_id] == NULL) {
            double memory_time = SHIFT_BEFORE_P_BRB_BP_SN + tUpEvent_brb + MAX_ACCEPTABLE_PICK_ERROR;
            num_local_wwssn_sp_memory = 1 + (int) (memory_time / deltaTime);
            pmem_local_wwssn_sp[source_id] = init_timedomain_memory(num_local_wwssn_sp_memory, -1.0, 0, 0.0, deltaTime);
        } else {
            num_local_wwssn_sp_memory = pmem_local_wwssn_sp[source_id]->numInput;
        }
        // data for BRB S/N
        data_float_wwssn_sp_sqr = copyData(data_float_wwssn_sp, numsamples);
        // square value
        applySqr(data_float_wwssn_sp_sqr, numsamples);
        data_float_wwssn_sp_rms = copyData(data_float_wwssn_sp_sqr, numsamples);
        // smooth for brb s/n
        applyBoxcarSmoothing(deltaTime, data_float_wwssn_sp_rms, numsamples, SMOOTHING_WINDOW_HALF_WIDTH_BRB_BP_RMS, &(pmem_smooth_wwssn_sp_rms[source_id]), 1);
        //int ntest = 10; printf("n=%d, data_float_raw=%f, data_float_wwssn_sp=%f, data_float_wwssn_sp_rms=%f\n", ntest, data_float_raw[ntest], data_float_wwssn_sp[ntest], data_float_wwssn_sp_rms[ntest]);
        // perpare to save wwssn_sp_rms to local memory, save enough samples to have wwssn_sp_rms values back tUpEvent and TIME_DELAY_A_REF seconds
        int num_local_wwssn_sp_sqr_memory;
        if (pmem_local_wwssn_sp_sqr[source_id] == NULL) {
            // 20110401 AJLdouble memory_time = TIME_DELAY_BRB_BP_RMS > tUpEvent_brb ? TIME_DELAY_BRB_BP_RMS : tUpEvent_brb;
            //memory_time += SHIFT_BEFORE_P_BRB_BP_SN;
            double memory_time = MAX_MWPD_DUR;
            num_local_wwssn_sp_sqr_memory = 1 + (int) (memory_time / deltaTime);
            //printf("@@@@@@@@DEBUG: num_local_wwssn_sp_sqr_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_wwssn_sp_sqr_memory, memory_time, deltaTime);
            pmem_local_wwssn_sp_sqr[source_id] = init_timedomain_memory(num_local_wwssn_sp_sqr_memory, -1.0, 0, 0.0, deltaTime);
        } else {
            num_local_wwssn_sp_sqr_memory = pmem_local_wwssn_sp_sqr[source_id]->numInput;
        }
        int num_local_wwssn_sp_rms_memory;
        if (pmem_local_wwssn_sp_rms[source_id] == NULL) {
            // 20110401 AJLdouble memory_time = TIME_DELAY_BRB_BP_RMS > tUpEvent_brb ? TIME_DELAY_BRB_BP_RMS : tUpEvent_brb;
            //memory_time += SHIFT_BEFORE_P_BRB_BP_SN;
            double memory_time = SHIFT_BEFORE_P_BRB_BP_SN + tUpEvent_brb + MAX_ACCEPTABLE_PICK_ERROR;
            num_local_wwssn_sp_rms_memory = 1 + (int) (memory_time / deltaTime);
            //printf("@@@@@@@@DEBUG: num_local_wwssn_sp_rms_memory [%d] = 1 + (int) (memory_time [%g] / deltaTime [%g])\n", num_local_wwssn_sp_rms_memory, memory_time, deltaTime);
            pmem_local_wwssn_sp_rms[source_id] = init_timedomain_memory(num_local_wwssn_sp_rms_memory, -1.0, 0, 0.0, deltaTime);
        } else {
            num_local_wwssn_sp_rms_memory = pmem_local_wwssn_sp_rms[source_id]->numInput;
        }
    }
    // tauc_t0 time series
    if (flag_do_tauc) {
        data_float_tauc = copyData(data_float_brb_hp, numsamples);
        // low-pass to remove hf, aviod (?) problems due to neglecting attenuation
        /*
        if (filter_lp_bu_co_5_n_6(deltaTime, data_float_tauc, numsamples, &pmem_filter_lp[source_id], 1) == NULL) {
            if (verbose)
                printf("ERROR: deltaTime not supported by filter_lp_bu_co_5_n_6: deltaTime=%f\n", deltaTime);
            channelParameters[source_id].error = channelParameters[source_id].error | ERROR_DT_NOT_SUPPORTED_BY_FILTER;
            return_value = -1;
            td_process_free_timedomain_memory(source_id);
            goto cleanup;
            //return (-1); // deltaTime not supported
        }
         */
        applyInstantPeriodWindowed(deltaTime, data_float_tauc, numsamples, INSTANT_PERIOD_WINDOWED_WINDOW_WIDTH,
                &pmem_inst_per[source_id], &pmem_inst_per_dval[source_id], 1);
    }


    // check for a_ref_ok_index where a_ref level below last_a_ref
    //    this prevents starting of new de calculation in "coda" of previous de calculation
    int a_ref_ok_index = 0;
    //printf("a_ref_ok_index %d  last_a_ref[source_id] %f\n", a_ref_ok_index, last_a_ref[source_id]);
    if (last_a_ref[source_id] >= 0.0) {
        a_ref_ok_index = numsamples;
        for (n = 0; n < numsamples; n++) {
            if (sqrt(data_float_a_ref[n]) < last_a_ref[source_id] * A_REF_OK_RATIO) {
                a_ref_ok_index = n;
                last_a_ref[source_id] = -1.0;
                last_a_ref_ellapsed_index_count[source_id] = -1;
                break;
            }
        }
    }

    // check for active (uncompleted) timedomain-processing data for this source_id
    int ndata;
    for (ndata = 0; ndata < num_de_data; ndata++) {

        TimedomainProcessingData* deData = data_list[ndata];

        // 20160811 AJL - debug
        if (deData == NULL) {
            printf("DEBUG: td_process_timedomain_processing(): deData == NULL: deData: %ld, ndata %d, num_de_data %d\n", (long) deData, ndata, num_de_data);
        }

        if (deData->source_id != source_id)
            continue;

        // flag and skip clipped or non_contiguous data     // 20111230 AJL - added
        if (deData->flag_clipped || deData->flag_non_contiguous) {
            //printf("DEBUG: pick flag_clipped %d || deData->flag_non_contiguous %d\n", deData->flag_clipped, deData->flag_non_contiguous);
            continue;
        }

        if ((deData->flag_complete_t50
                && (!(flag_do_tauc || flag_do_mwp || flag_do_mwpd || flag_do_grd_vel) || deData->flag_complete_snr_brb)
                && (!(flag_do_mwp || flag_do_mwpd) || deData->flag_complete_snr_brb_int)
                && (!flag_do_grd_vel || deData->flag_complete_grd_mot)
                && (!flag_do_mb || deData->flag_complete_snr_brb_bp)
                && (!flag_do_mwp || deData->flag_complete_mwp)
                && (!flag_do_t0 || deData->flag_complete_t0)
                && (!flag_do_mwpd || deData->flag_complete_mwpd)
                && (!flag_do_mb || deData->flag_complete_mb))) {
            continue;
        }

        // flag and skip new clipped or non_contiguous data that is not completed
        if (flag_clipped || flag_non_contiguous) { // 20111230 AJL - added
            deData->flag_clipped = flag_clipped;
            deData->flag_non_contiguous = flag_non_contiguous;
            //printf("DEBUG: pick newly flag_clipped %d || deData->flag_non_contiguous %d\n", deData->flag_clipped, deData->flag_non_contiguous);
            continue;
        }

        int ioffset_pick = deData->virtual_pick_index;

        // find virtual pick index of next pick  // 20110310 AJL
        /*int ioffset_pick_hf_next = INT_MAX;
        int nd;
        for (nd = ndata + 1; nd < num_de_data; nd++) {
            TimedomainProcessingData* deDataTest = data_list[nd];
            if (deDataTest->source_id == source_id && deDataTest->pick_stream == STREAM_HF) {
                ioffset_pick_hf_next = deDataTest->virtual_pick_index;
                break;
            }
        }*/


        // check for HF S/N (a_ref value at pick time)
        //printf("DEBUG: pick %d/%d, pick_ndx %d, ioffset_pick %d, sn_pick %f\n", n, num_de_data, ioffset_pick, ioffset_pick, deData->sn_pick);
        if (deData->sn_pick < 0.0) { // do not yet have sn_pick
            if (ioffset_pick >= 0 && ioffset_pick < numsamples) { // value in current a_ref data array
                deData->sn_pick = sqrt(data_float_a_ref[ioffset_pick]);
                if (verbose > 1)
                    printf("a_ref S/N: pick %d/%d, pick_ndx %d, sn_pick %f\n", ndata, num_de_data, ioffset_pick, deData->sn_pick);
            } else if (ioffset_pick >= -pmem_local_a_ref[source_id]->numInput && ioffset_pick < 0) { // value in memory a_ref data array
                int index = ioffset_pick + pmem_local_a_ref[source_id]->numInput;
                deData->sn_pick = sqrt(pmem_local_a_ref[source_id]->input[index]);
                if (verbose > 1)
                    printf("a_ref_memory S/N: pick %d/%d, pick_ndx %d, sn_pick %f\n", ndata, num_de_data, ioffset_pick, deData->sn_pick);
            } else {
                printf("ERROR: HF S/N: attempt to access pmem_local_a_ref array of size=%d at invalid index=%d: this should not happen!\n",
                        pmem_local_a_ref[source_id]->numInput, ioffset_pick + pmem_local_a_ref[source_id]->numInput);
                continue; // should never reach here
            }
        }

        // set A_REF index and offset
        int a_ref_index = (int) (0.5 + TIME_DELAY_A_REF / deltaTime);
        int ioffset_a_ref = a_ref_index + ioffset_pick;

        // check for max and following min of smooth T50 data (SMOOTHING_WINDOW_HALF_WIDTH_T50)
        //   in the window from pick -> A_REF index
        //   this allows identification of smaller evens where HF level drops rapidly after P pick,
        //   to avoid high de values due to S-P ~ 50s or depth_phases-P (pP-P, sP-P, ...) ~ 50s
        if (!deData->flag_a_ref_not_ok && ioffset_a_ref > 0) { // A_REF index at or after beginning of this packet
            int ioffset_begin = ioffset_pick >= 0 ? ioffset_pick : 0;
            int ioffset_end = ioffset_a_ref < numsamples ? ioffset_a_ref + 1 : numsamples;
            double amp_max = deData->amp_max_in_aref_window;
            double amp_min = deData->amp_min_in_aref_window;
            double amp;
            for (n = ioffset_begin; n < ioffset_end; n++) {
                amp = data_float_t50[n];
                if (amp > amp_max) {
                    amp_max = amp;
                    amp_min = amp_max; // reset min amp
                } else if (amp < amp_min) {
                    amp_min = amp;
                }
            }
            deData->amp_max_in_aref_window = amp_max;
            deData->amp_min_in_aref_window = amp_min;
            // if at or past end of A_REF index, check ratio amp_min / amp_max
            if (ioffset_a_ref < numsamples) { // A_REF index at or before end of this packet
                if (sqrt(deData->amp_min_in_aref_window / deData->amp_max_in_aref_window) < A_REF_WINDOW_MIN_TO_MAX_AMP_RATIO_MIN) {
                    deData->t50 = 0.0; // force t50Ex level = 0
                    if (verbose > 1)
                        printf("%s (%d) A_REF Amin/Amax too low: pick %d/%d, pick_ndx %d, t50_index %d, t50 %f, Amin/Amax %f\n",
                            station, source_id, ndata, num_de_data, ioffset_pick, a_ref_index, deData->t50,
                            sqrt(deData->amp_min_in_aref_window / deData->amp_max_in_aref_window));
                    deData->flag_complete_t50 = 1; // flag end of t50 processing for this pick data
                }
            }
        }

        // check for A_REF
        if (ioffset_a_ref >= 0 && ioffset_a_ref < numsamples) {
            // check a_ref index is after a_ref_ok_index
            //printf("ioffset_a_ref %d < a_ref_ok_index %d\n", ioffset_a_ref, a_ref_ok_index);
            if (!no_aref_level_check && ioffset_a_ref < a_ref_ok_index) {
                deData->flag_a_ref_not_ok = 1;
                if (verbose > 99)
                    printf("%s (%d) A_REF:  deData->flag_a_ref_not_ok = 1\n", station, source_id);
            } else {
                // AJL20090519 - test using amp_max_in_aref_window
                deData->a_ref = sqrt(data_float_a_ref[ioffset_a_ref]);
                //deData->a_ref = 0.9 * sqrt(deData->amp_max_in_aref_window);
                // END - AJL20090519
                last_a_ref[source_id] = deData->a_ref;
                last_a_ref_ellapsed_index_count[source_id] = numsamples - ioffset_a_ref;
                deData->a_ref /= gain_microns;
                deData->sn_pick /= gain_microns;
                // check s/n ratio
                if (deData->sn_pick < FLT_MIN || (deData->a_ref / deData->sn_pick) < SIGNAL_TO_NOISE_RATIO_HF_MIN)
                    deData->flag_snr_hf_too_low = 1;
                // re-check for a_ref_ok_index where a_ref level below last_a_ref
                a_ref_ok_index = numsamples;
                for (n = ioffset_a_ref; n < numsamples; n++) {
                    if (sqrt(data_float_a_ref[n]) < last_a_ref[source_id] * A_REF_OK_RATIO) {
                        a_ref_ok_index = n;
                        last_a_ref[source_id] = -1.0;
                        last_a_ref_ellapsed_index_count[source_id] = -1;
                        break;
                    }
                }
                if (verbose > 1)
                    printf("%s (%d) A_REF: pick %d/%d, pick_ndx %d, a_ref_index %d, a_ref %f\n",
                        station, source_id, ndata, num_de_data, ioffset_pick, a_ref_index, deData->a_ref);
            }
        }

        // process t0
        if (flag_do_t0 && !deData->flag_complete_t0) {
            if (deData->t0->amp_noise < 0.0) { // do not yet have noise amplitude
                if (ioffset_pick >= 0 && ioffset_pick < numsamples) { // value in current t50 data array
                    deData->t0->amp_noise = data_float_t50[ioffset_pick];
                    if (verbose > 1)
                        printf("t0 noise amplitude: pick %d/%d, pick_ndx %d, amp %f\n", ndata, num_de_data, ioffset_pick, deData->t0->amp_noise);
                } else if (ioffset_pick >= -pmem_local_t50[source_id]->numInput && ioffset_pick < 0) { // value in memory t50 data array
                    int index = ioffset_pick + pmem_local_t50[source_id]->numInput;
                    deData->t0->amp_noise = pmem_local_t50[source_id]->input[index];
                    if (verbose > 1)
                        printf("t0_memeory noise amplitude: pick %d/%d, pick_ndx %d, amp %f\n", ndata, num_de_data, ioffset_pick, deData->t0->amp_noise);
                }
            }
            int ioffset_t0_min = ioffset_pick;
            // 20110401 AJL
            if (ioffset_t0_min < 0 && deData->t0->have_used_memory)
                ioffset_t0_min = 0;
            else
                deData->t0->have_used_memory = 1;
            int t0_end_index = (int) (0.5 + MAX_MWPD_DUR / deltaTime);
            int ioffset_t0_max = t0_end_index + ioffset_pick;
            if (ioffset_t0_max < numsamples)
                deData->flag_complete_t0 = 1; // will complete
            if (ioffset_t0_max > numsamples - 1)
                ioffset_t0_max = numsamples - 1;
            // check for local and global amplitude extrema for peak and fractions of peak
            double amplitude;
            for (n = ioffset_t0_min; n <= ioffset_t0_max; n++) {
                // 20110401 AJL
                if (n >= 0) {
                    amplitude = data_float_t50[n];
                } else if (n >= -pmem_local_t50[source_id]->numInput) { // value in memory brb_hp data array
                    int index = n + pmem_local_t50[source_id]->numInput;
                    amplitude = pmem_local_t50[source_id]->input[index];
                } else {
                    printf("ERROR: T0: attempt to access pmem_local_t50 array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_t50[source_id]->numInput, n + pmem_local_t50[source_id]->numInput);
                    continue; // should never reach here
                }
                double amp_peak = deData->t0->amp_peak;
                // check if at next pick and amplitude is below 20% level
                // 20110310 AJL - added to avoid erroneously long T0 durations (and Mwpd) when P waves from a later event arrive soon after this pick
                /*if (amplitude < 0.05 * amp_peak && n >= ioffset_pick_hf_next - (int) (0.5 + SMOOTHING_WINDOW_HALF_WIDTH_T50 / deltaTime)) { // at next pick and amplitude is below 20% level
                    deData->flag_complete_t0 = 1; // then declare complete
                    break;
                }*/
                if (amplitude > deData->t0->amp_peak) { // amplitude rises above peak - unset all level amplitudes
                    deData->t0->amp_peak = amplitude;
                    deData->t0->index_peak = n - ioffset_pick;
                    deData->t0->amp_90 = -1.0;
                    deData->t0->amp_80 = -1.0;
                    deData->t0->amp_50 = -1.0;
                    deData->t0->amp_20 = -1.0;
                    //deData->t0->duration_raw = DURATION_INVALID;  // not needed, better to leave last duration estimate if available
                } else {
                    if (deData->t0->amp_90 < 0.0) {
                        if (amplitude <= 0.9 * amp_peak) { // amplitude drops below level - set level amplitude and index
                            deData->t0->amp_90 = amplitude;
                            deData->t0->index_90 = n - ioffset_pick;
                        }
                    } else {
                        if (amplitude > 0.9 * amp_peak) { // amplitude rises above level - unset level amplitudes
                            deData->t0->amp_90 = -1.0;
                            deData->t0->amp_80 = -1.0;
                            deData->t0->amp_50 = -1.0;
                            deData->t0->amp_20 = -1.0;
                        }
                    }
                    if (deData->t0->amp_80 < 0.0) {
                        if (amplitude <= 0.8 * amp_peak) { // amplitude drops below level - set level amplitude and index
                            deData->t0->amp_80 = amplitude;
                            deData->t0->index_80 = n - ioffset_pick;
                        }
                    } else {
                        if (amplitude > 0.8 * amp_peak) { // amplitude rises above level - unset level amplitudes
                            deData->t0->amp_80 = -1.0;
                            deData->t0->amp_50 = -1.0;
                            deData->t0->amp_20 = -1.0;
                        }
                    }
                    if (deData->t0->amp_50 < 0.0) {
                        if (amplitude <= 0.5 * amp_peak) { // amplitude drops below level - set level amplitude and index
                            deData->t0->amp_50 = amplitude;
                            deData->t0->index_50 = n - ioffset_pick;
                        }
                    } else {
                        if (amplitude > 0.5 * amp_peak) { // amplitude rises above level - unset level amplitudes
                            deData->t0->amp_50 = -1.0;
                            deData->t0->amp_20 = -1.0;
                        }
                    }
                    if (deData->t0->amp_20 < 0.0) {
                        if (amplitude <= 0.2 * amp_peak) { // amplitude drops below level - set level amplitude and index
                            deData->t0->amp_20 = amplitude;
                            deData->t0->index_20 = n - ioffset_pick;
                            calculate_duration(deData); // must have all levels, calculate duration
                        }
                    } else {
                        if (amplitude > 0.2 * amp_peak) { // amplitude rises above level - unset level amplitude
                            deData->t0->amp_20 = -1.0;
                        }
                    }
                    if (n - ioffset_pick > (int) (0.5 + MIN_MWPD_DUR / deltaTime) && // passed minimum Mwpd duration
                            ((amplitude / deData->t0->amp_noise < SIGNAL_TO_NOISE_RATIO_T0_END) // amplitude below signal/noise level - terminate t0 calcaulation
                            || (amplitude / deData->t0->amp_peak < SIGNAL_TO_PEAK_RATIO_T0_END)) // amplitude below signal/peak level - terminate t0 calcaulation
                            ) {
                        deData->flag_complete_t0 = 1; // declare complete
                        break;
                    }
                }
                // check if amplitude has been below 50% level for > 2 * (time of 50% level when last duration was set)
                // 20110310 AJL - added to avoid erroneously long T0 durations (and Mwpd) when P waves from a later event arrive soon after this pick
                /*if (amplitude < 0.5 * amp_peak) { // amplitude is below 50% level
                    if (deData->t0->duration_raw != T0_INVALID // a duration has been set
                            && n - ioffset_pick > 2 * deData->t0->index_50_last // current P offset time > 2 * (time of 50% level when last duration was set)
                            && n - ioffset_pick > (int) (0.5 + 2.0 * SMOOTHING_WINDOW_HALF_WIDTH_T50 / deltaTime) // current P offset time > 2 * SMOOTHING_WINDOW_HALF_WIDTH_T50
                            ) {
                        //printf("DEBUG: t0 complete:\n   > pick %s_%s, n - ioffset_pick %d, 2 * deData->t0->index_20 %d\n", deData->network, deData->station, n - ioffset_pick, 2 * deData->t0->index_20);
                        //printf("   > n - ioffset_pick %d, (int) (0.5 + 2.0 * SMOOTHING_WINDOW_HALF_WIDTH_T50 / deltaTime) %d\n", n - ioffset_pick, (int) (0.5 + 2.0 * SMOOTHING_WINDOW_HALF_WIDTH_T50 / deltaTime));
                        deData->flag_complete_t0 = 1; // then declare complete
                        break;
                    }
                } else {    // reset last 50% level index
                    deData->t0->index_50_last = n - ioffset_pick;
                }*/
                // check if time since P is > 2 * duration
                // 20110310 AJL - added to avoid erroneously long T0 durations (and large Mwpd) when P waves from a later event arrive soon after this pick
                if (deData->t0->duration_raw != T0_INVALID // a duration has been set
                        && deData->t0->amp_20 > 0.0 // amplitude is < 20% peak
                        && n - ioffset_pick > (int) (0.5 + 2.0 * deData->t0->duration_raw / deltaTime) // current P offset time > 2 * duration
                        && n - ioffset_pick > (int) (0.5 + 2.0 * SMOOTHING_WINDOW_HALF_WIDTH_T50 / deltaTime) // current P offset time > 2 * SMOOTHING_WINDOW_HALF_WIDTH_T50
                        ) {
                    deData->flag_complete_t0 = 1; // then declare complete
                    break;
                }
            }
        }

        // check for BRB S/N (brb_hp_rms value at pick time)
        //printf("DEBUG: pick %d/%d, pick_ndx %d, ioffset_pick %d, sn_brb_pick %f\n", n, num_de_data, ioffset_pick, ioffset_pick, deData->sn_brb_pick);
        if ((flag_do_tauc || flag_do_mwp || flag_do_mwpd || flag_do_grd_vel)) {
            if (deData->sn_brb_pick < 0.0) { // do not yet have sn_brb_pick
                int ioffset_snr_brb_pick = ioffset_pick - (int) (0.5 + SHIFT_BEFORE_P_BRB_HP_SN / deltaTime); // get rms reference amplitude in advance of P time
                if (ioffset_snr_brb_pick >= 0 && ioffset_snr_brb_pick < numsamples) { // value in current brb_hp_rms data array
                    deData->sn_brb_pick = sqrt(data_float_brb_hp_rms[ioffset_snr_brb_pick]);
                    if (verbose > 1)
                        printf("brb_hp_rms S/N: pick %d/%d, pick_ndx %d, sn_brb_pick %f\n", ndata, num_de_data, ioffset_snr_brb_pick, deData->sn_brb_pick);
                } else if (ioffset_snr_brb_pick >= -pmem_local_brb_hp_rms[source_id]->numInput && ioffset_snr_brb_pick < 0) { // value in memory brb_hp_rms data array
                    int index = ioffset_snr_brb_pick + pmem_local_brb_hp_rms[source_id]->numInput;
                    deData->sn_brb_pick = sqrt(pmem_local_brb_hp_rms[source_id]->input[index]);
                    if (verbose > 1)
                        printf("brb_hp_rms_memory S/N: pick %d/%d, pick_ndx %d, sn_brb_pick %f\n", ndata, num_de_data, ioffset_snr_brb_pick, deData->sn_brb_pick);
                } else {
                    printf("ERROR: BRB HP S/N: %s_%s_%s_%s: attempt to access pmem_local_brb_hp_rms array of size=%d at invalid index=%d: this should not happen!\n",
                            deData->network, deData->station, deData->location, deData->channel, pmem_local_brb_hp_rms[source_id]->numInput, ioffset_snr_brb_pick + pmem_local_brb_hp_rms[source_id]->numInput);
                    continue; // should never reach here
                }
            }
            // check for BRB hp s/n
            if (!deData->flag_complete_snr_brb && deData->sn_brb_pick >= 0.0) {
                if (USE_TO_FOR_SIGNAL_LEVEL_WINDOW) { // 2012115 AJL - added
                    // 20130204 AJL - modified to use maximum window available after P for provisional s/n estimate
                    //          - avoids long wait before declaring s/n too low, which can give temporary, too-high magnitudes
                    double signal_level_window = -1.0;
                    int t0_available = deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID;
                    if (t0_available) {
                        signal_level_window = deData->t0->duration_raw < TIME_DELAY_BRB_HP_RMS ? deData->t0->duration_raw : TIME_DELAY_BRB_HP_RMS;
                        deData->flag_complete_snr_brb = 1; // declare complete
                        //printf("TP 1: win %f  max %f", signal_level_window, TIME_DELAY_BRB_HP_RMS);
                    } else {
                        int snr_brb_index_max = (int) (0.5 + TIME_DELAY_BRB_HP_RMS / deltaTime);
                        int ioffset_snr_brb = snr_brb_index_max + ioffset_pick;
                        if (ioffset_snr_brb <= -snr_brb_index_max) { // current packet starts later than 2 * TIME_DELAY_BRB_HP_RMS after P
                            signal_level_window = TIME_DELAY_BRB_HP_RMS;
                            deData->flag_complete_snr_brb = 1; // declare complete
                            //printf("TP 2: win %f  max %f", signal_level_window, TIME_DELAY_BRB_HP_RMS);
                        } else {
                            // use current maximum window length after P for provisional s/n, do not declare complete
                            signal_level_window = (double) (numsamples - ioffset_pick) * deltaTime;
                            if (signal_level_window > TIME_DELAY_BRB_HP_RMS)
                                signal_level_window = TIME_DELAY_BRB_HP_RMS;
                            //printf("TP 3: win %f  max %f", signal_level_window, TIME_DELAY_BRB_HP_RMS);
                        }
                    }
                    if (signal_level_window > 0.0) {
                        // signal window should not end past S arrival time  // 20151117 AJL added - important to prevent much too high P s/n
                        if (deData->ttime_SminusP > 0.0) {
                            signal_level_window = signal_level_window < deData->ttime_SminusP ? signal_level_window : deData->ttime_SminusP;
                        }
                        deData->flag_snr_brb_too_low = 0;
                        deData->sn_brb_signal = calculate_signal_level("BRB_HP", source_id, signal_level_window, data_float_brb_hp_sqr, pmem_local_brb_hp_sqr,
                                ioffset_pick, numsamples, deltaTime);
                        if (deData->sn_brb_pick < FLT_MIN || (deData->sn_brb_signal / deData->sn_brb_pick) < SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN)
                            deData->flag_snr_brb_too_low = 1;
                        if (verbose > 1)
                            printf("%s (%d) S/N-BRB-HP: pick %d/%d, pick_ndx %d, sn_brb_pick %f, signal_level_window %f, sn_signal %f, snr_brb %f\n",
                                station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_pick, signal_level_window, deData->sn_brb_signal,
                                deData->sn_brb_pick < FLT_MIN ? -1.0 : deData->sn_brb_signal / deData->sn_brb_pick);
                    }
                }
                /*else {
                    int snr_brb_index = (int) (0.5 + TIME_DELAY_BRB_HP_RMS / deltaTime);
                    int ioffset_snr_brb = snr_brb_index + ioffset_pick;
                    if (ioffset_snr_brb >= 0 && ioffset_snr_brb < numsamples) {
                        // check s/n ratio
                        deData->sn_brb_signal = sqrt(data_float_brb_hp_rms[ioffset_snr_brb]);
                        if (deData->sn_brb_pick < FLT_MIN || (deData->sn_brb_signal / deData->sn_brb_pick) < SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN)
                            deData->flag_snr_brb_too_low = 1;
                        deData->flag_complete_snr_brb = 1;
                        if (verbose > 1)
                            printf("%s (%d) S/N-BRB-HP: pick %d/%d, pick_ndx %d, sn_brb_pick %f, snr_brb_index %d, sn_signal %f, snr_brb %f\n",
                                station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_pick, snr_brb_index, deData->sn_brb_signal,
                                deData->sn_brb_pick < FLT_MIN ? -1.0 : deData->sn_brb_signal / deData->sn_brb_pick);
                    }
                }*/
            }
            // check for BRB int (disp) S/N (brb_hp_int_rms value at pick time)
            //printf("DEBUG: pick %d/%d, pick_ndx %d, ioffset_pick %d, sn_brb_int_pick %f\n", n, num_de_data, ioffset_pick, ioffset_pick, deData->sn_brb_int_pick);
            if ((flag_do_mwp || flag_do_mwpd)) {
                if (deData->sn_brb_int_pick < 0.0) { // do not yet have sn_brb_int_pick
                    int ioffset_snr_brb_int_pick = ioffset_pick - (int) (0.5 + SHIFT_BEFORE_P_BRB_HP_INT_SN / deltaTime); // get rms reference amplitude in advance of P time
                    if (ioffset_snr_brb_int_pick >= 0 && ioffset_snr_brb_int_pick < numsamples) { // value in current brb_hp_rms data array
                        deData->sn_brb_int_pick = sqrt(data_float_brb_hp_int_rms[ioffset_snr_brb_int_pick]);
                        if (verbose > 1)
                            printf("brb_hp_rms S/N: pick %d/%d, pick_ndx %d, sn_brb_int_pick %f\n", ndata, num_de_data, ioffset_snr_brb_int_pick, deData->sn_brb_int_pick);
                    } else if (ioffset_snr_brb_int_pick >= -pmem_local_brb_hp_int_rms[source_id]->numInput && ioffset_snr_brb_int_pick < 0) { // value in memory brb_hp_rms data array
                        int index = ioffset_snr_brb_int_pick + pmem_local_brb_hp_int_rms[source_id]->numInput;
                        deData->sn_brb_int_pick = sqrt(pmem_local_brb_hp_int_rms[source_id]->input[index]);
                        if (verbose > 1)
                            printf("brb_hp_rms_memory S/N: pick %d/%d, pick_ndx %d, sn_brb_int_pick %f\n", ndata, num_de_data, ioffset_snr_brb_int_pick, deData->sn_brb_int_pick);
                    } else {
                        printf("ERROR: BRB HP S/N: %s_%s_%s_%s: attempt to access pmem_local_brb_hp_int_rms array of size=%d at invalid index=%d: this should not happen!\n",
                                deData->network, deData->station, deData->location, deData->channel, pmem_local_brb_hp_int_rms[source_id]->numInput, ioffset_snr_brb_int_pick + pmem_local_brb_hp_int_rms[source_id]->numInput);
                        continue; // should never reach here
                    }
                }
                // check for BRB hp int s/n
                if (!deData->flag_complete_snr_brb_int && deData->sn_brb_int_pick >= 0.0) {
                    if (USE_TO_FOR_SIGNAL_LEVEL_WINDOW) { // 2012115 AJL - added
                        // 20130204 AJL - modified to use maximum window available after P for provisional s/n estimate
                        //          - avoids long wait before declaring s/n too low, which can give temporary, too-high magnitudes
                        double signal_level_window = -1.0;
                        int t0_available = deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID;
                        if (t0_available) {
                            signal_level_window = deData->t0->duration_raw < TIME_DELAY_BRB_HP_INT_RMS ? deData->t0->duration_raw : TIME_DELAY_BRB_HP_INT_RMS;
                            deData->flag_complete_snr_brb_int = 1; // declare complete
                        } else {
                            int snr_brb_int_index_max = (int) (0.5 + TIME_DELAY_BRB_HP_INT_RMS / deltaTime);
                            int ioffset_snr_brb_int = snr_brb_int_index_max + ioffset_pick;
                            if (ioffset_snr_brb_int <= -snr_brb_int_index_max) { // current packet starts later than 2 * TIME_DELAY_BRB_HP_INT_RMS after P
                                signal_level_window = TIME_DELAY_BRB_HP_INT_RMS;
                                deData->flag_complete_snr_brb_int = 1; // declare complete
                            } else {
                                // use current maximum window length after P for provisional s/n, do not declare complete
                                signal_level_window = (double) (numsamples - ioffset_pick) * deltaTime;
                                if (signal_level_window > TIME_DELAY_BRB_HP_INT_RMS)
                                    signal_level_window = TIME_DELAY_BRB_HP_INT_RMS;
                            }
                        }
                        if (signal_level_window > 0.0) {
                            // signal window should not end past S arrival time  // 20151117 AJL added - important to prevent much too high P s/n
                            if (deData->ttime_SminusP > 0.0) {
                                signal_level_window = signal_level_window < deData->ttime_SminusP ? signal_level_window : deData->ttime_SminusP;
                            }
                            deData->flag_snr_brb_int_too_low = 0;
                            deData->sn_brb_int_signal = calculate_signal_level("BRB_HP_INT", source_id, signal_level_window, data_float_brb_hp_int_sqr, pmem_local_brb_hp_int_sqr,
                                    ioffset_pick, numsamples, deltaTime);
                            if (deData->sn_brb_int_pick < FLT_MIN || (deData->sn_brb_int_signal / deData->sn_brb_int_pick) < SIGNAL_TO_NOISE_RATIO_BRB_HP_INT_MIN)
                                deData->flag_snr_brb_int_too_low = 1;
                            if (verbose > 1)
                                printf("%s (%d) S/N-BRB-HP-INT: pick %d/%d, pick_ndx %d, sn_brb_int_pick %f, signal_level_window %f, sn_signal %f, snr_brb_int %f\n",
                                    station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_int_pick, signal_level_window, deData->sn_brb_int_signal,
                                    deData->sn_brb_int_pick < FLT_MIN ? -1.0 : deData->sn_brb_int_signal / deData->sn_brb_int_pick);
                        }
                    }
                    /*else {
                        int snr_brb_int_index = (int) (0.5 + TIME_DELAY_BRB_HP_INT_RMS / deltaTime);
                        int ioffset_snr_brb_int = snr_brb_int_index + ioffset_pick;
                        if (ioffset_snr_brb_int >= 0 && ioffset_snr_brb_int < numsamples) {
                            // check s/n ratio
                            deData->sn_brb_int_signal = sqrt(data_float_brb_hp_int_rms[ioffset_snr_brb_int]);
                            if (deData->sn_brb_int_pick < FLT_MIN || (deData->sn_brb_int_signal / deData->sn_brb_int_pick) < SIGNAL_TO_NOISE_RATIO_BRB_HP_INT_MIN)
                                deData->flag_snr_brb_int_too_low = 1;
                            deData->flag_complete_snr_brb_int = 1;
                            if (verbose > 1)
                                printf("%s (%d) S/N-BRB-HP-INT: pick %d/%d, pick_ndx %d, sn_brb_int_pick %f, snr_brb_int_index %d, sn_signal %f, snr_brb_int %f\n",
                                    station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_int_pick, snr_brb_int_index, deData->sn_brb_int_signal,
                                    deData->sn_brb_int_pick < FLT_MIN ? -1.0 : deData->sn_brb_int_signal / deData->sn_brb_int_pick);
                        }
                    }*/
                }
            }
        }

        // check for BRB S/N (wwssn_sp_rms value at pick time)
        //printf("DEBUG: pick %d/%d, pick_ndx %d, ioffset_pick %d, sn_brb_bp_pick %f\n", n, num_de_data, ioffset_pick, ioffset_pick, deData->sn_brb_bp_pick);
        if (flag_do_mb) {
            if (deData->sn_brb_bp_pick < 0.0) { // do not yet have sn_brb_bp_pick
                int ioffset_snr_wwssn_sp_pick = ioffset_pick - (int) (0.5 + SHIFT_BEFORE_P_BRB_BP_SN / deltaTime);
                if (ioffset_snr_wwssn_sp_pick >= 0 && ioffset_snr_wwssn_sp_pick < numsamples) { // value in current wwssn_sp_rms data array
                    deData->sn_brb_bp_pick = sqrt(data_float_wwssn_sp_rms[ioffset_snr_wwssn_sp_pick]);
                    if (verbose > 1)
                        printf("wwssn_sp_rms S/N: pick %d/%d, pick_ndx %d, sn_brb_bp_pick %f\n", ndata, num_de_data, ioffset_snr_wwssn_sp_pick, deData->sn_brb_bp_pick);
                } else if (ioffset_snr_wwssn_sp_pick >= -pmem_local_wwssn_sp_rms[source_id]->numInput && ioffset_snr_wwssn_sp_pick < 0) { // value in memory wwssn_sp_rms data array
                    int index = ioffset_snr_wwssn_sp_pick + pmem_local_wwssn_sp_rms[source_id]->numInput;
                    deData->sn_brb_bp_pick = sqrt(pmem_local_wwssn_sp_rms[source_id]->input[index]);
                    if (verbose > 1)
                        printf("wwssn_sp_rms_memory S/N: pick %d/%d, pick_ndx %d, sn_brb_bp_pick %f\n", ndata, num_de_data, ioffset_snr_wwssn_sp_pick, deData->sn_brb_bp_pick);
                } else {
                    printf("ERROR: BRB BP S/N: attempt to access pmem_local_wwssn_sp_rms array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_wwssn_sp_rms[source_id]->numInput, ioffset_snr_wwssn_sp_pick + pmem_local_wwssn_sp_rms[source_id]->numInput);
                    continue; // should never reach here
                }
            }
            // check for BRB bp s/n
            if (!deData->flag_complete_snr_brb_bp && deData->sn_brb_bp_pick >= 0.0) {
                if (USE_TO_FOR_SIGNAL_LEVEL_WINDOW) { // 2012115 AJL - added
                    // 20130204 AJL - modified to use maximum window available after P for provisional s/n estimate
                    //          - avoids long wait before declaring s/n too low, which can give temporary, too-high magnitudes
                    double signal_level_window = -1.0;
                    int t0_available = deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID;
                    if (t0_available) {
                        signal_level_window = deData->t0->duration_raw < TIME_DELAY_BRB_BP_RMS ? deData->t0->duration_raw : TIME_DELAY_BRB_BP_RMS;
                        deData->flag_complete_snr_brb_bp = 1; // declare complete
                    } else {
                        int snr_brb_bp_index_max = (int) (0.5 + TIME_DELAY_BRB_BP_RMS / deltaTime);
                        int ioffset_snr_brb_bp = snr_brb_bp_index_max + ioffset_pick;
                        if (ioffset_snr_brb_bp <= -snr_brb_bp_index_max) { // current packet starts later than 2 * TIME_DELAY_BRB_BP_RMS after P
                            signal_level_window = TIME_DELAY_BRB_BP_RMS;
                            deData->flag_complete_snr_brb_bp = 1; // declare complete
                        } else {
                            // use current maximum window length after P for provisional s/n, do not declare complete
                            signal_level_window = (double) (numsamples - ioffset_pick) * deltaTime;
                            if (signal_level_window > TIME_DELAY_BRB_BP_RMS)
                                signal_level_window = TIME_DELAY_BRB_BP_RMS;
                        }
                    }
                    if (signal_level_window > 0.0) {
                        // signal window should not end past S arrival time  // 20151117 AJL added - important to prevent much too high P s/n
                        if (deData->ttime_SminusP > 0.0) {
                            signal_level_window = signal_level_window < deData->ttime_SminusP ? signal_level_window : deData->ttime_SminusP;
                        }
                        deData->flag_snr_brb_bp_too_low = 0;
                        deData->sn_brb_bp_signal = calculate_signal_level("BRB_BP", source_id, signal_level_window, data_float_wwssn_sp_sqr, pmem_local_wwssn_sp_sqr,
                                ioffset_pick, numsamples, deltaTime);
                        if (deData->sn_brb_bp_pick < FLT_MIN || (deData->sn_brb_bp_signal / deData->sn_brb_bp_pick) < SIGNAL_TO_NOISE_RATIO_BRB_BP_MIN)
                            deData->flag_snr_brb_bp_too_low = 1;
                        if (verbose > 1)
                            printf("%s (%d) S/N-BRB-HP-INT: pick %d/%d, pick_ndx %d, sn_brb_bp_pick %f, signal_level_window %f, sn_signal %f, snr_brb_bp %f\n",
                                station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_bp_pick, signal_level_window, deData->sn_brb_bp_signal,
                                deData->sn_brb_bp_pick < FLT_MIN ? -1.0 : deData->sn_brb_bp_signal / deData->sn_brb_bp_pick);
                    }
                }
                /*else {
                    int snr_wwssn_sp_index = (int) (0.5 + TIME_DELAY_BRB_BP_RMS / deltaTime);
                    int ioffset_snr_wwssn_sp = snr_wwssn_sp_index + ioffset_pick;
                    if (ioffset_snr_wwssn_sp >= 0 && ioffset_snr_wwssn_sp < numsamples) {
                        // check s/n ratio
                        deData->sn_brb_bp_signal = sqrt(data_float_wwssn_sp_rms[ioffset_snr_wwssn_sp]);
                        if (deData->sn_brb_bp_pick < FLT_MIN || (deData->sn_brb_bp_signal / deData->sn_brb_bp_pick) < SIGNAL_TO_NOISE_RATIO_BRB_BP_MIN)
                            deData->flag_snr_brb_bp_too_low = 1;
                        deData->flag_complete_snr_brb_bp = 1;
                        if (verbose > 1)
                            printf("%s (%d) S/N-BRB-BP: pick %d/%d, pick_ndx %d, sn_brb_bp_pick %f, n_brb_bp_index %d, sn_brb_bp_signal %f, n_brb_bp %f\n",
                                station, source_id, ndata, num_de_data, ioffset_pick, deData->sn_brb_bp_pick, snr_wwssn_sp_index, deData->sn_brb_bp_signal,
                                deData->sn_brb_bp_pick < FLT_MIN ? -1.0 : deData->sn_brb_bp_signal / deData->sn_brb_bp_pick);
                    }
                }*/
            }
        }

        // save BRB HP ground motion for misc processing
        if (flag_do_grd_vel && flag_have_gain && !deData->flag_complete_grd_mot) {
            int ioffset_shift_before_p = (int) (0.5 + START_ANALYSIS_BEFORE_P_BRB_HP / deltaTime); // get reference amplitude in advance of P time
            int ioffset_grd_vel_min = ioffset_pick - ioffset_shift_before_p;
            deData->grd_mot->amp_at_pick = 0.0;
            if (ioffset_grd_vel_min < 0 && deData->grd_mot->have_used_memory)
                ioffset_grd_vel_min = 0;
            else
                deData->grd_mot->have_used_memory = 1;
            int grd_vel_end_index = (int) (0.5 + MAX_GRD_VEL_DUR / deltaTime);
            int ioffset_grd_vel_max = grd_vel_end_index + ioffset_pick;
            if (ioffset_grd_vel_max < numsamples)
                deData->flag_complete_grd_mot = 1; // will complete
            if (ioffset_grd_vel_max > numsamples - 1)
                ioffset_grd_vel_max = numsamples - 1;
            // get amplitude, store vel and disp ground motion
            int index_count_after_analysis_start = deData->grd_mot->n_amp;
            double grd_vel_amp_at_pick = deData->grd_mot->amp_at_pick;
            double amplitude, grd_vel_amplitude;
            for (n = ioffset_grd_vel_min; n <= ioffset_grd_vel_max; n++) {
                index_count_after_analysis_start++;
                if (index_count_after_analysis_start >= deData->grd_mot->n_amp_max)
                    printf("ERROR: grd_vel: index_count_after_analysis_start=%d  >= deData->grd_vel->n_amp_max=%d: this should not happen!\n",
                        index_count_after_analysis_start, deData->grd_mot->n_amp_max);
                if (n >= 0) {
                    amplitude = data_float_brb_hp[n];
                } else if (n >= -pmem_local_brb_hp[source_id]->numInput) { // value in memory brb_hp data array
                    int index = n + pmem_local_brb_hp[source_id]->numInput;
                    amplitude = pmem_local_brb_hp[source_id]->input[index];
                } else {
                    printf("ERROR: grd_vel: attempt to access pmem_local_brb_hp array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_brb_hp[source_id]->numInput, n + pmem_local_brb_hp[source_id]->numInput);
                    continue; // should never reach here
                }
                grd_vel_amplitude = (amplitude - grd_vel_amp_at_pick) / gain_microns;
                deData->grd_mot->vel[index_count_after_analysis_start] = grd_vel_amplitude;
                deData->grd_mot->disp[index_count_after_analysis_start] = grd_vel_amplitude * deltaTime;
            }
            deData->grd_mot->ioffset_pick = ioffset_pick;
            deData->grd_mot->n_amp = index_count_after_analysis_start;
        }

        // process tauc
        if (flag_do_tauc) {
            // find tauc peak value
            int tauc_index_min = (int) (0.5 + TIME_DELAY_TAUC_MIN / deltaTime);
            int ioffset_tauc_min = tauc_index_min + ioffset_pick;
            if (ioffset_tauc_min < 0)
                ioffset_tauc_min = 0;
            int tauc_index_max = (int) (0.5 + TIME_DELAY_TAUC_MAX / deltaTime);
            // 20120109 AJL - added duration limit to tauc determination
            if (deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID && deData->t0->duration_raw < TIME_DELAY_TAUC_MAX)
                tauc_index_max = (int) (0.5 + deData->t0->duration_raw / deltaTime); // end tauc peak determination at T0 duration
            int ioffset_tauc_max = tauc_index_max + ioffset_pick;
            if (ioffset_tauc_max > numsamples - 1)
                ioffset_tauc_max = numsamples - 1;
            for (n = ioffset_tauc_min; n <= ioffset_tauc_max; n++) {
                if (data_float_tauc[n] > deData->tauc_peak)
                    deData->tauc_peak = data_float_tauc[n];
            }
        }

        // process mb
        if (flag_do_mb && flag_have_gain && !deData->flag_complete_mb) {
            // 20110328 AJL - Bug fix:  int ioffset_mb_min = ioffset_pick - START_ANALYSIS_BEFORE_P_MB;
            int ioffset_shift_before_p = (int) (0.5 + START_ANALYSIS_BEFORE_P_MB / deltaTime); // get reference amplitude in advance of P time
            int ioffset_mb_min = ioffset_pick - ioffset_shift_before_p;
            /* 20110630 AJL
            // 20110401 AJL
            if (ioffset_mb_min >= 0 && ioffset_mb_min < numsamples) {
                //deData->mb->amp_at_pick = data_float_wwssn_sp[ioffset_mb_min];
                deData->mb->amp_at_pick = data_float_wwssn_sp[ioffset_mb_min];
            } else if (ioffset_mb_min >= -pmem_local_wwssn_sp[source_id]->numInput && ioffset_mb_min < 0) { // value in memory brb_hp data array
                int index = ioffset_mb_min + pmem_local_wwssn_sp[source_id]->numInput;
                deData->mb->amp_at_pick = pmem_local_wwssn_sp[source_id]->input[index];
            }*/
            // 20110630 AJL - data_float_wwssn_sp is SP, narrow-band, so mean should be zero, and better to use zero than amp_at_pick
            deData->mb->amp_at_pick = 0.0;
            // 20110401 AJL
            //ioffset_mb_min = ioffset_pick; // start integration at P time
            if (ioffset_mb_min < 0 && deData->mb->have_used_memory)
                ioffset_mb_min = 0;
            else
                deData->mb->have_used_memory = 1;
            int mb_end_index = (int) (0.5 + MAX_MB_DUR / deltaTime);
            if (deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID && deData->t0->duration_raw < MAX_MB_DUR)
                mb_end_index = (int) (0.5 + deData->t0->duration_raw / deltaTime); // end magnitude calculation at T0 duration (see Bormann & Saul 2008)
            int ioffset_mb_max = mb_end_index + ioffset_pick;
            if (ioffset_mb_max < numsamples)
                deData->flag_complete_mb = 1; // will complete
            if (ioffset_mb_max > numsamples - 1)
                ioffset_mb_max = numsamples - 1;
            // accumulate amplitudes in window
            int index_count_after_analysis_start = deData->mb->n_amplitude;
            double mb_amp_at_pick = deData->mb->amp_at_pick;
            double mb_int_sum = deData->mb->int_sum;
            double polarity = deData->mb->polarity;
            double peak_index_count = deData->mb->peak_index_count;
            double amp_wwssn_sp;
            double amplitdue_sp;
            double peak_dur;
            for (n = ioffset_mb_min; n <= ioffset_mb_max; n++) {
                index_count_after_analysis_start++;
                if (index_count_after_analysis_start >= deData->mb->n_amplitude_max)
                    printf("ERROR: mb: index_count_after_analysis_start=%d  >= deData->mb->n_amplitude_max=%d: this should not happen!\n",
                        index_count_after_analysis_start, deData->mb->n_amplitude_max);
                // 20110401 AJL
                if (n >= 0) {
                    amp_wwssn_sp = data_float_wwssn_sp[n];
                } else if (n >= -pmem_local_wwssn_sp[source_id]->numInput) { // value in memory wwssn_sp data array
                    int index = n + pmem_local_wwssn_sp[source_id]->numInput;
                    amp_wwssn_sp = pmem_local_wwssn_sp[source_id]->input[index];
                } else {
                    printf("ERROR: mb: attempt to access pmem_local_wwssn_sp array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_wwssn_sp[source_id]->numInput, n + pmem_local_wwssn_sp[source_id]->numInput);
                    continue; // should never reach here
                }
                if (0 && strcmp(deData->station, "DUG") == 0) {
                    printf("TEST: %s  n %d: amp_wwssn_sp %g  mb_amp_at_pick %g  gain %g\n",
                            deData->station, index_count_after_analysis_start, amp_wwssn_sp, mb_amp_at_pick, gain_microns);
                }
                if (MB_MODE == MB_MODE_mb_DISP) {
                    // transfer function is WWSSN_SP_DISP, gives velocity->velocity, integrate to displacement as needed for mb
                    amplitdue_sp = mb_int_sum = mb_int_sum + (amp_wwssn_sp - mb_amp_at_pick) * (deltaTime / gain_microns);
                } else if (MB_MODE == MB_MODE_mB) {
                    // transfer function is WWSSN_SP_DISP, gives velocity->velocity, so already have velocity as needed for mB
                    amplitdue_sp = (amp_wwssn_sp - mb_amp_at_pick) / gain_microns;
                } else if (MB_MODE == MB_MODE_mb_VEL) {
                    // transfer function is WWSSN_SP_VEL, gives velocity->disp, so already have displacement as needed for mb
                    amplitdue_sp = (amp_wwssn_sp - mb_amp_at_pick) / gain_microns;
                }
                //if (strcmp(deData->station, "HGN") == 0)
                //    printf("amplitdue_sp %g = amplitdue_sp + (amp_wwssn_sp %g - mb_amp_at_pick %g) * (deltaTime %g / gain %g)\n",
                //        amplitdue_sp, amp_wwssn_sp, mb_amp_at_pick, deltaTime, gain);
                // mB formulation:  amplitude = (amp_wwssn_sp - mb_amp_at_pick) / gain;
                // check displacement polarity (mb formulation only)
                peak_dur = -1.0;
                if (amplitdue_sp >= 0.0) { // pos
                    // check if past extremum
                    if (polarity < 0) { // passed from neg to pos displacement
                        peak_dur = deltaTime * (double) (peak_index_count);
                        peak_index_count = 0;
                    }
                    polarity = 1;
                } else { // neg
                    // check if past extremum
                    if (polarity > 0) { // passed from pos to neg displacement
                        peak_dur = deltaTime * (double) (peak_index_count);
                        peak_index_count = 0;
                    }
                    polarity = -1;
                }
                peak_index_count++;
                // store abs amplitude
                deData->mb->amplitude[index_count_after_analysis_start] = fabs(amplitdue_sp);
                // mB formulation:  deData->mb->amplitude[index_count_after_analysis_start] = fabs(amplitude);
                // mB formulation:  index_count_after_analysis_start++;
                deData->mb->peak_dur[index_count_after_analysis_start] = peak_dur;
                if (peak_dur > 0.0) {
                    // back-fill peak dur array where not set (under just ended displacement pulse)
                    int ndx;
                    for (ndx = index_count_after_analysis_start - 1; ndx >= 0; ndx--) {
                        if (deData->mb->peak_dur[ndx] >= 0.0)
                            break;
                        deData->mb->peak_dur[ndx] = peak_dur;
                    }
                }
            }
            deData->mb->int_sum = mb_int_sum;
            deData->mb->polarity = polarity;
            deData->mb->peak_index_count = peak_index_count;
            deData->mb->n_amplitude = index_count_after_analysis_start;
        }

        // process mwp
        if (flag_do_mwp && flag_have_gain && !deData->flag_complete_mwp) {
            // 20110328 AJL - Bug fix: int ioffset_mwp_min = ioffset_pick - START_ANALYSIS_BEFORE_P_BRB_HP;
            int ioffset_shift_before_p = (int) (0.5 + START_ANALYSIS_BEFORE_P_BRB_HP / deltaTime); // get reference amplitude in advance of P time
            int ioffset_mwp_min = ioffset_pick - ioffset_shift_before_p;
            // 20110401 AJL
            /*if (ioffset_mwp_min >= 0 && ioffset_mwp_min < numsamples) {
                //deData->mwp->amp_at_pick = data_float_brb_hp[ioffset_mwp_min];
                deData->mwp->amp_at_pick = data_float_brb_hp[ioffset_mwp_min];
            } else if (ioffset_mwp_min >= -pmem_local_brb_hp[source_id]->numInput && ioffset_mwp_min < 0) { // value in memory brb_hp data array
                int index = ioffset_mwp_min + pmem_local_brb_hp[source_id]->numInput;
                deData->mwp->amp_at_pick = pmem_local_brb_hp[source_id]->input[index];
            }*/
            // 20121115 AJL - changed from above
            deData->mwp->amp_at_pick = 0.0;
            // 20110401 AJL
            //ioffset_mwp_min = ioffset_pick; // start integration at P time
            if (ioffset_mwp_min < 0 && deData->mwp->have_used_memory)
                ioffset_mwp_min = 0;
            else
                deData->mwp->have_used_memory = 1;
            int mwp_end_index = (int) (0.5 + MAX_MWP_DUR / deltaTime);
            // should not end Mwp at T0 duration - Mwp is defined to use a 120s window to capture pP or sP pulse
            // but Mwp may be overestimated, especially if there is a P arrival from another event  // 20110310 AJL - try ending Mwp at T0 duration
            if (deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID && deData->t0->duration_raw < MAX_MWP_DUR)
                mwp_end_index = (int) (0.5 + deData->t0->duration_raw / deltaTime); // end magnitude calculation at T0 duration
            //
            int ioffset_mwp_max = mwp_end_index + ioffset_pick;
            if (ioffset_mwp_max < numsamples)
                deData->flag_complete_mwp = 1; // will complete
            if (ioffset_mwp_max > numsamples - 1)
                ioffset_mwp_max = numsamples - 1;
            // apply double integration, check for extrema
            int index_count_after_analysis_start = deData->mwp->n_peak;
            double mwp_amp_at_pick = deData->mwp->amp_at_pick;
            double mwp_int_sum = deData->mwp->int_sum;
            double mwp_int_int_sum = deData->mwp->int_int_sum;
            double mwp_abs_int_int_sum = deData->mwp->abs_int_int_sum;
            double mwp_int_int_sum_mwp_peak = deData->mwp->int_int_sum_mwp_peak;
            double polarity = deData->mwp->polarity;
            double peak_index_count = deData->mwp->peak_index_count;
            double peak = deData->mwp->peak[index_count_after_analysis_start];
            double peak_dur = deData->mwp->peak_dur[index_count_after_analysis_start];
            double amplitude, d_integral;
            for (n = ioffset_mwp_min; n <= ioffset_mwp_max; n++) {
                index_count_after_analysis_start++;
                if (index_count_after_analysis_start >= deData->mwp->n_peak_max)
                    printf("ERROR: Mwp: index_count_after_analysis_start=%d  >= deData->mwp->n_peak_max=%d: this should not happen!\n",
                        index_count_after_analysis_start, deData->mwp->n_peak_max);
                // 20110401 AJL
                if (n >= 0) {
                    amplitude = data_float_brb_hp[n];
                } else if (n >= -pmem_local_brb_hp[source_id]->numInput) { // value in memory brb_hp data array
                    int index = n + pmem_local_brb_hp[source_id]->numInput;
                    amplitude = pmem_local_brb_hp[source_id]->input[index];
                } else {
                    printf("ERROR: Mwp: attempt to access pmem_local_brb_hp array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_brb_hp[source_id]->numInput, n + pmem_local_brb_hp[source_id]->numInput);
                    continue; // should never reach here
                }
                // integrate to displacement
                // 20110401 AJL
                d_integral = (amplitude - mwp_amp_at_pick) * deltaTime / gain_meters;
                mwp_int_sum += d_integral; // displacement
                mwp_int_int_sum += mwp_int_sum * deltaTime; // integral of displacement
                mwp_abs_int_int_sum += fabs(mwp_int_sum * deltaTime); // integral of abs displacement
                // check displacement polarity
                if (mwp_int_sum >= 0.0) { // pos
                    // check if past extremum
                    if (polarity < 0) { // passed from neg to pos displacement
                        mwp_int_int_sum_mwp_peak *= -1.0;
                        if (mwp_int_int_sum_mwp_peak > peak) {
                            peak = mwp_int_int_sum_mwp_peak;
                            peak_dur = deltaTime * (double) (peak_index_count);
                        }
                        /*if (n > ioffset_pick + (int) (0.5 + 1.0 / deltaTime) && deData->mwp->pulse_first_dur < 0.0) {
                            deData->mwp->pulse_first = -mwp_int_int_sum_mwp_peak;
                            deData->mwp->pulse_first_dur = deltaTime * (double) (peak_index_count);
                        }*/
                        mwp_int_int_sum_mwp_peak = 0;
                        peak_index_count = 0;
                    }
                    polarity = 1;
                } else { // neg
                    // check if past extremum
                    if (polarity > 0) { // passed from pos to neg displacement
                        if (mwp_int_int_sum_mwp_peak > peak) {
                            peak = mwp_int_int_sum_mwp_peak;
                            peak_dur = deltaTime * (double) (peak_index_count);
                        }
                        /*if (n > ioffset_pick + (int) (0.5 + 1.0 / deltaTime) && deData->mwp->pulse_first_dur < 0.0) {
                            deData->mwp->pulse_first = mwp_int_int_sum_mwp_peak;
                            deData->mwp->pulse_first_dur = deltaTime * (double) (peak_index_count);
                        }*/
                        mwp_int_int_sum_mwp_peak = 0;
                        peak_index_count = 0;
                    }
                    polarity = -1;
                }
                mwp_int_int_sum_mwp_peak += mwp_int_sum * deltaTime;
                deData->mwp->peak[index_count_after_analysis_start] = peak;
                deData->mwp->peak_dur[index_count_after_analysis_start] = peak_dur;
                deData->mwp->int_int[index_count_after_analysis_start] = mwp_int_int_sum;
                deData->mwp->int_int_abs[index_count_after_analysis_start] = mwp_abs_int_int_sum;
                if (n > ioffset_pick) // 20121115 AJL - bug fix
                    peak_index_count++;
            }
            deData->mwp->int_sum = mwp_int_sum;
            deData->mwp->int_int_sum = mwp_int_int_sum;
            deData->mwp->abs_int_int_sum = mwp_abs_int_int_sum;
            deData->mwp->int_int_sum_mwp_peak = mwp_int_int_sum_mwp_peak;
            deData->mwp->polarity = polarity;
            deData->mwp->peak_index_count = peak_index_count;
            deData->mwp->n_peak = index_count_after_analysis_start;
        }

        // process mwpd
        if (flag_do_mwpd && flag_have_gain && !deData->flag_complete_mwpd) {
            // 20110323 AJL - Bug fix: int ioffset_mwpd_min = ioffset_pick - START_ANALYSIS_BEFORE_P_BRB_HP;
            int ioffset_shift_before_p = (int) (0.5 + START_ANALYSIS_BEFORE_P_BRB_HP / deltaTime); // get reference amplitude in advance of P time
            int ioffset_mwpd_min = ioffset_pick - ioffset_shift_before_p;
            // 20110401 AJL
            /*if (ioffset_mwpd_min >= 0 && ioffset_mwpd_min < numsamples) {
                //deData->mwpd->amp_at_pick = data_float_brb_hp[ioffset_mwpd_min];
                deData->mwpd->amp_at_pick = data_float_brb_hp[ioffset_mwpd_min];
            } else if (ioffset_mwpd_min >= -pmem_local_brb_hp[source_id]->numInput && ioffset_mwpd_min < 0) { // value in memory brb_hp data array
                int index = ioffset_mwpd_min + pmem_local_brb_hp[source_id]->numInput;
                deData->mwpd->amp_at_pick = pmem_local_brb_hp[source_id]->input[index];
            }*/
            // 20121115 AJL - changed from above
            deData->mwpd->amp_at_pick = 0.0;
            // 20110401 AJL
            //ioffset_mwpd_min = ioffset_pick; // start integration at P time
            if (ioffset_mwpd_min < 0 && deData->mwpd->have_used_memory)
                ioffset_mwpd_min = 0;
            else
                deData->mwpd->have_used_memory = 1;

            int mwpd_end_index = (int) (0.5 + MAX_MWPD_DUR / deltaTime);
            if (deData->flag_complete_t0 && deData->t0->duration_raw != T0_INVALID && deData->t0->duration_raw < MAX_MWPD_DUR)
                mwpd_end_index = (int) (0.5 + deData->t0->duration_raw / deltaTime); // end magnitude calcaulation at T0 duration
            int ioffset_mwpd_max = mwpd_end_index + ioffset_pick;
            if (ioffset_mwpd_max < numsamples)
                deData->flag_complete_mwpd = 1; // will complete
            if (ioffset_mwpd_max > numsamples - 1)
                ioffset_mwpd_max = numsamples - 1;
            // apply double integration, check for extrema
            // 20110323 AJL - Bug fix:  int index_count_after_analysis_start = ioffset_mwpd_min - (ioffset_pick - START_ANALYSIS_BEFORE_P_BRB_HP);
            int index_count_after_analysis_start = deData->mwpd->n_int_sum;
            double mwpd_amp_at_pick = deData->mwpd->amp_at_pick;
            double mwpd_int_sum = deData->mwpd->int_sum;
            double mwpd_int_int_sum_pos = deData->mwpd->int_int_sum_pos[index_count_after_analysis_start];
            double mwpd_int_int_sum_neg = deData->mwpd->int_int_sum_neg[index_count_after_analysis_start];
            double amplitude;
            for (n = ioffset_mwpd_min; n <= ioffset_mwpd_max; n++) {
                index_count_after_analysis_start++;
                if (index_count_after_analysis_start >= deData->mwpd->n_int_sum_max)
                    printf("ERROR: Mwpd: index_count_after_analysis_start=%d  >= deData->mwpd->n_int_sum_max=%d: this should not happen!\n",
                        index_count_after_analysis_start, deData->mwpd->n_int_sum_max);
                // 20110401 AJL
                if (n >= 0) {
                    amplitude = data_float_brb_hp[n];
                } else if (n >= -pmem_local_brb_hp[source_id]->numInput) { // value in memory brb_hp data array
                    int index = n + pmem_local_brb_hp[source_id]->numInput;
                    amplitude = pmem_local_brb_hp[source_id]->input[index];
                } else {
                    printf("ERROR: Mwpd: attempt to access pmem_local_brb_hp array of size=%d at invalid index=%d: this should not happen!\n",
                            pmem_local_brb_hp[source_id]->numInput, n + pmem_local_brb_hp[source_id]->numInput);
                    continue; // should never reach here
                }
                // integral to displacement
                // 20110401 AJL
                mwpd_int_sum = mwpd_int_sum + (amplitude - mwpd_amp_at_pick) * deltaTime / gain_meters;
                // check displacement polarity
                if (mwpd_int_sum >= 0.0) { // pos
                    mwpd_int_int_sum_pos += mwpd_int_sum * deltaTime;
                } else { // neg
                    mwpd_int_int_sum_neg += -1.0 * mwpd_int_sum * deltaTime;
                }
                deData->mwpd->int_int_sum_pos[index_count_after_analysis_start] = mwpd_int_int_sum_pos;
                deData->mwpd->int_int_sum_neg[index_count_after_analysis_start] = mwpd_int_int_sum_neg;
            }
            deData->mwpd->n_int_sum = index_count_after_analysis_start;
            deData->mwpd->int_sum = mwpd_int_sum;
            //deData->flag_complete_mwpd = deData->flag_complete_t0; // Mwpd calculate should finish at same time as t0 calculation
        }

        // process T50
        if (!deData->flag_a_ref_not_ok && !deData->flag_complete_t50) {
            int t50_index = (int) (0.5 + TIME_DELAY_T50 / deltaTime);
            int ioffset_t50 = t50_index + ioffset_pick;
            if (verbose > 99)
                printf("%s (%d)  T50: ioffset_t50 %d  numsamples %d\n", station, source_id, ioffset_t50, numsamples);
            if (ioffset_t50 >= 0 && ioffset_t50 < numsamples) {
                deData->t50 = sqrt(data_float_t50[ioffset_t50]);
                deData->t50 /= gain_microns;
                if (verbose > 1)
                    printf("%s (%d) T50: pick %d/%d, pick_ndx %d, t50_index %d, t50 %f, t50/a_ref %f\n",
                        station, source_id, ndata, num_de_data, ioffset_pick, t50_index, deData->t50, deData->t50 / deData->a_ref);
                deData->flag_complete_t50 = 1; // flag end to t50 processing of this pick data
            }
        }

        // increment index count
        //if (!deData->flag_complete_t50 || (flag_do_mwp && !deData->flag_complete_mwp)
        //        || (flag_do_t0 && !deData->flag_complete_t0) || (flag_do_mwpd && !deData->flag_complete_mwpd) || (flag_do_mb && !deData->flag_complete_mb))
        deData->virtual_pick_index -= numsamples;
    }

    // if a_ref valid time passed, reset last_a_ref level
    if (last_a_ref[source_id] >= 0.0) {
        if (last_a_ref_ellapsed_index_count[source_id] > (int) (A_REF_VALID_TIME / deltaTime)) {
            last_a_ref[source_id] = -1.0;
            last_a_ref_ellapsed_index_count[source_id] = -1;
        } else { // 20150408 AJL - bug fix?
            last_a_ref_ellapsed_index_count[source_id] += numsamples;
        }
    }

    // save local memory
    updateInput(pmem_local_t50[source_id], data_float_t50, numsamples);
    updateInput(pmem_local_a_ref[source_id], data_float_a_ref, numsamples);
    if (flag_do_tauc || flag_do_mwp || flag_do_mwpd || flag_do_grd_vel) {
        updateInput(pmem_local_brb_hp[source_id], data_float_brb_hp, numsamples);
        updateInput(pmem_local_brb_hp_sqr[source_id], data_float_brb_hp_sqr, numsamples);
        updateInput(pmem_local_brb_hp_rms[source_id], data_float_brb_hp_rms, numsamples);
        if (flag_do_mwp || flag_do_mwpd) {
            updateInput(pmem_local_brb_hp_int_sqr[source_id], data_float_brb_hp_int_sqr, numsamples);
            updateInput(pmem_local_brb_hp_int_rms[source_id], data_float_brb_hp_int_rms, numsamples);
        }
    }
    if (flag_do_mb) {
        updateInput(pmem_local_wwssn_sp[source_id], data_float_wwssn_sp, numsamples);
        updateInput(pmem_local_wwssn_sp_sqr[source_id], data_float_wwssn_sp_sqr, numsamples);
        updateInput(pmem_local_wwssn_sp_rms[source_id], data_float_wwssn_sp_rms, numsamples);
    }

    // clean up
cleanup:
    ;

    // set returned data
    if (return_value >= 0) {
        // sets stream to return
        float **data_float_return = &data_float_raw;
        //float **data_float_return = &data_float_brb_hp;
        //float **data_float_return = &data_float_wwssn_sp_rms;
        //
        float *data_float_temp = *pdata_float_in; // save ptr to input data stream
        *pdata_float_in = *data_float_return; // puts stream to return into input data stream ptr argument (which is available to calling code)
        *data_float_return = data_float_temp; // make sure non-returned stream memory will be freed below
    }

    if (flag_do_tauc || flag_do_mwp || flag_do_mwpd || flag_do_grd_vel) {
        free(data_float_brb_hp);
        free(data_float_brb_hp_sqr);
        free(data_float_brb_hp_rms);
        if (flag_do_mwp || flag_do_mwpd) {
            free(data_float_brb_hp_int_sqr);
            free(data_float_brb_hp_int_rms);
        }
    }
    if (flag_do_mb) {
        free(data_float_wwssn_sp);
        free(data_float_wwssn_sp_sqr);
        free(data_float_wwssn_sp_rms);
    }
    if (flag_do_tauc)
        free(data_float_tauc);
    free(data_float_t50);
    free(data_float_a_ref);
    free(data_float_hf);

    return (return_value);

} /* End of td_process_timedomain_processing() */

/***************************************************************************
 * td_process_timedomain_processing_init:
 *
 * Do necessary initialization.
 ***************************************************************************/
int td_process_timedomain_processing_init(Settings *settings, char* outpath, char* geogfile, char* gainfile, int verbose) {

    printf("%s\n%s v%s (%s)\n%s %s\n%s\n%s\n%s\n", EARLY_EST_MONITOR_NAME, EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE,
            EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_BANNER_1, EARLY_EST_MONITOR_BANNER_2, EARLY_EST_MONITOR_BANNER_3, EARLY_EST_MONITOR_BANNER_4);

    printf("Info: Velocity model / travel-time grid: %s\n", get_model_ttime_name());

    app_prop_settings = settings;

    int int_param;
    double double_param;
    char out_buf[SETTINGS_MAX_STR_LEN];


    // flag to enable magnitude and tsunami discriminant measures.
    measuresEnable = 1;
    if ((int_param = settings_get_int(app_prop_settings, "TimedomainProcessing", "measures.enable")) != INT_INVALID)
        measuresEnable = int_param;
    printf("Info: property set: [TimedomainProcessing] measures.enable %d\n", measuresEnable);


    // set pick parameters
    // NOTE: do not change defaults!  These values are standard for global EarlyEst.
    // SEE: _DOC_ in FilterPicker5.c for more details on the picker parameters

    // set ordered (0-N), comma separated list channel names for mapping pick parameters
    // default is use first listed parameter for all channel
    if (settings_get(app_prop_settings, "TimedomainProcessing", "pick.channel.list", out_buf, SETTINGS_MAX_STR_LEN)) {
        num_pick_channels = 0;
        printf("Info: property set: [TimedomainProcessing] pick.channel.list %s -> ", out_buf);
        char *str_pos = strtok(out_buf, ", ");
        while (str_pos != NULL && num_pick_channels < MAX_NUM_PICK_CHANNEL_LIST) {
            *pick_channel_list[num_pick_channels] = '\0';
            strcpy(pick_channel_list[num_pick_channels], str_pos);
            printf(" <%s>,", pick_channel_list[num_pick_channels]);
            num_pick_channels++;
            str_pos = strtok(NULL, ",");
        }
        printf("\n"); // DEBUG
    }


    // comma separated list of channel orientations to be picked and processed
    //  default: all channel orientations picked and processed
    //  Note: superseeds channel orientations in pick.channel.list
    if (settings_get(app_prop_settings, "TimedomainProcessing", "process.orientation.list", out_buf, SETTINGS_MAX_STR_LEN)) {
        num_process_orientations = 0;
        printf("Info: property set: [TimedomainProcessing] process.orientation.list %s -> ", out_buf);
        char *str_pos = strtok(out_buf, ", ");
        while (str_pos != NULL && num_process_orientations < MAX_NUM_PROCESS_ORIENTATION_LIST) {
            *process_orientation_list[num_process_orientations] = '\0';
            strcpy(process_orientation_list[num_process_orientations], str_pos);
            printf(" <%s>,", process_orientation_list[num_process_orientations]);
            num_process_orientations++;
            str_pos = strtok(NULL, ",");
        }
        printf("\n"); // DEBUG
    }


    // polarization
    // flag to use P polarization azimuth for location
    polarization_enable = 0;
    if ((int_param = settings_get_int(app_prop_settings, "TimedomainProcessing", "polarization.enable")) != INT_INVALID)
        polarization_enable = int_param;
    printf("Info: property set: [TimedomainProcessing] polarization.enable %d\n", polarization_enable);
    //
    polarization_window_start_delay_after_P = POLARIZATION_START_DELAY_AFTER_P;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "polarization.window.start_delay_after_P")) != DBL_INVALID)
        polarization_window_start_delay_after_P = double_param;
    printf("Info: property set: [TimedomainProcessing] polarization.window.start_delay_after_P %f\n", polarization_window_start_delay_after_P);
    //
    polarization_window_length_max = POLARIZATION_WINDOW__LENGTH_MAX;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "polarization.window.length.max")) != DBL_INVALID)
        polarization_window_length_max = double_param;
    printf("Info: property set: [TimedomainProcessing] polarization.window.length.max %f\n", polarization_window_length_max);
    //
    polarization_window_length_min = POLARIZATION_WINDOW__LENGTH_MIN;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "polarization.window.length.min")) != DBL_INVALID)
        polarization_window_length_min = double_param;
    printf("Info: property set: [TimedomainProcessing] polarization.window.length.min %f\n", polarization_window_length_min);


    // flag to use hf data (high-frequency) for picking.
    pickOnHFStream = 1;
    if ((int_param = settings_get_int(app_prop_settings, "TimedomainProcessing", "pick.hf.enable")) != INT_INVALID)
        pickOnHFStream = int_param;
    printf("Info: property set: [TimedomainProcessing] pick.hf.enable %d\n", pickOnHFStream);
    // flag to use raw data (broadband) for picking.
    pickOnRawStream = 1;
    if ((int_param = settings_get_int(app_prop_settings, "TimedomainProcessing", "pick.raw.enable")) != INT_INVALID)
        pickOnRawStream = int_param;
    printf("Info: property set: [TimedomainProcessing] pick.raw.enable %d\n", pickOnRawStream);

    // raw picker
    // PICKER="filp_5 filtw 100.0 ltw 50.0 thres1 10.0 thres2 10.0 tupevt 2.0 res PICKS"
    // long filter window helps to get TauC Td for early, LP P signal (e.g. for tsunami earthquakes),
    //   compensate with lower MAX_PICK_PERIOD_ASSOC_LOCATE in ttimes_location.h
    rawPickParams[0].filterWindow = 100.0;
    rawPickParams[0].longTermWindow = 50.0;
    rawPickParams[0].threshold1 = 10.0;
    rawPickParams[0].threshold2 = 10.0;
    rawPickParams[0].tUpEvent = 2.0; // must be less than TIME_DELAY_TAUC_MIN
    //*
    double pick_params[MAX_NUM_PICK_CHANNEL_LIST];
    int n_params_returned;
    int n_params, nparam;

    //double settings_get_double_tuple(const Settings *settings, const char *section, const char *key, double *out, unsigned int n_out) {
    // set pick parameters from properties file
    n_params = 1;
    if ((n_params_returned = settings_get_double_tuple(app_prop_settings, "TimedomainProcessing", "pick.raw.filterWindow", pick_params, MAX_NUM_PICK_CHANNEL_LIST)) > 0) {
        for (nparam = 0; nparam < n_params_returned; nparam++) {
            rawPickParams[nparam].filterWindow = pick_params[nparam];
        }
        n_params = n_params_returned;
    }
    printf("Info: property set: [TimedomainProcessing] pick.raw.filterWindow ");
    for (nparam = 0; nparam < n_params; nparam++) {
        printf("%f, ", rawPickParams[nparam].filterWindow);
    }
    printf("\n");
    for (nparam = n_params; nparam < MAX_NUM_PICK_CHANNEL_LIST; nparam++) { // initialize all remaining channel list parameters
        rawPickParams[nparam].filterWindow = pick_params[n_params - 1];
    }
    //
    n_params = 1;
    if ((n_params_returned = settings_get_double_tuple(app_prop_settings, "TimedomainProcessing", "pick.raw.longTermWindow", pick_params, MAX_NUM_PICK_CHANNEL_LIST)) > 0) {
        for (nparam = 0; nparam < n_params_returned; nparam++) {
            rawPickParams[nparam].longTermWindow = pick_params[nparam];
        }
        n_params = n_params_returned;
    }
    printf("Info: property set: [TimedomainProcessing] pick.raw.longTermWindow ");
    for (nparam = 0; nparam < n_params; nparam++) {
        printf("%f, ", rawPickParams[nparam].longTermWindow);
    }
    printf("\n");
    for (nparam = n_params; nparam < MAX_NUM_PICK_CHANNEL_LIST; nparam++) { // initialize all remaining channel list parameters
        rawPickParams[nparam].longTermWindow = pick_params[n_params - 1];
    }
    //
    n_params = 1;
    if ((n_params_returned = settings_get_double_tuple(app_prop_settings, "TimedomainProcessing", "pick.raw.threshold1", pick_params, MAX_NUM_PICK_CHANNEL_LIST)) > 0) {
        for (nparam = 0; nparam < n_params_returned; nparam++) {
            rawPickParams[nparam].threshold1 = pick_params[nparam];
        }
        n_params = n_params_returned;
    }
    printf("Info: property set: [TimedomainProcessing] pick.raw.threshold1 ");
    for (nparam = 0; nparam < n_params; nparam++) {
        printf("%f, ", rawPickParams[nparam].threshold1);
    }
    printf("\n");
    for (nparam = n_params; nparam < MAX_NUM_PICK_CHANNEL_LIST; nparam++) { // initialize all remaining channel list parameters
        rawPickParams[nparam].threshold1 = pick_params[n_params - 1];
    }
    //
    n_params = 1;
    if ((n_params_returned = settings_get_double_tuple(app_prop_settings, "TimedomainProcessing", "pick.raw.threshold2", pick_params, MAX_NUM_PICK_CHANNEL_LIST)) > 0) {
        for (nparam = 0; nparam < n_params_returned; nparam++) {
            rawPickParams[nparam].threshold2 = pick_params[nparam];
        }
        n_params = n_params_returned;
    }
    printf("Info: property set: [TimedomainProcessing] pick.raw.threshold2 ");
    for (nparam = 0; nparam < n_params; nparam++) {
        printf("%f, ", rawPickParams[nparam].threshold2);
    }
    printf("\n");
    for (nparam = n_params; nparam < MAX_NUM_PICK_CHANNEL_LIST; nparam++) { // initialize all remaining channel list parameters
        rawPickParams[nparam].threshold2 = pick_params[n_params - 1];
    }
    //
    n_params = 1;
    if ((n_params_returned = settings_get_double_tuple(app_prop_settings, "TimedomainProcessing", "pick.raw.tUpEvent", pick_params, MAX_NUM_PICK_CHANNEL_LIST)) > 0) {
        for (nparam = 0; nparam < n_params_returned; nparam++) {
            rawPickParams[nparam].tUpEvent = pick_params[nparam];
        }
        n_params = n_params_returned;
    }
    printf("Info: property set: [TimedomainProcessing] pick.raw.tUpEvent ");
    for (nparam = 0; nparam < n_params; nparam++) {
        printf("%f, ", rawPickParams[nparam].tUpEvent);
    }
    printf("\n");
    for (nparam = n_params; nparam < MAX_NUM_PICK_CHANNEL_LIST; nparam++) { // initialize all remaining channel list parameters
        rawPickParams[nparam].tUpEvent = pick_params[n_params - 1];
    }
    //*/
    /*
    // set pick parameters from properties file
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.raw.filterWindow")) != DBL_INVALID)
        rawPickParams[0].filterWindow = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.raw.longTermWindow")) != DBL_INVALID)
        rawPickParams[0].longTermWindow = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.raw.threshold1")) != DBL_INVALID)
        rawPickParams[0].threshold1 = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.raw.threshold2")) != DBL_INVALID)
        rawPickParams[0].threshold2 = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.raw.tUpEvent")) != DBL_INVALID)
        rawPickParams[0].tUpEvent = double_param;
    printf("Info: property set: [TimedomainProcessing] pick.raw. filterWindow %f longTermWindow %f threshold1 %f threshold2 %f tUpEvent %f\n",
            rawPickParams[0].filterWindow, rawPickParams[0].longTermWindow, rawPickParams[0].threshold1, rawPickParams[0].threshold2, rawPickParams[0].tUpEvent);
     * */

    // hf picker
    hfPickParams[0].filterWindow = 10.0;
    hfPickParams[0].longTermWindow = 25.0;
    hfPickParams[0].threshold1 = 10.0;
    hfPickParams[0].threshold2 = 8.0;
    hfPickParams[0].tUpEvent = 2.0; // must be less than TIME_DELAY_TAUC_MIN
    // set pick parameters from properties file
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.hf.filterWindow")) != DBL_INVALID)
        hfPickParams[0].filterWindow = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.hf.longTermWindow")) != DBL_INVALID)
        hfPickParams[0].longTermWindow = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.hf.threshold1")) != DBL_INVALID)
        hfPickParams[0].threshold1 = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.hf.threshold2")) != DBL_INVALID)
        hfPickParams[0].threshold2 = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "TimedomainProcessing", "pick.hf.tUpEvent")) != DBL_INVALID)
        hfPickParams[0].tUpEvent = double_param;
    printf("Info: property set: [TimedomainProcessing] pick.hf. filterWindow %f longTermWindow %f threshold1 %f threshold2 %f tUpEvent %f\n",
            hfPickParams[0].filterWindow, hfPickParams[0].longTermWindow, hfPickParams[0].threshold1, hfPickParams[0].threshold2, hfPickParams[0].tUpEvent);


#ifdef USE_OPENMP
    //
    openmp_parameters.openmp_use = OPENMP_USE_DEFAULT;
    if ((int_param = settings_get_int(app_prop_settings, "ParallelProcessing", "openmp.use")) != INT_INVALID) {
        openmp_parameters.openmp_use = int_param;
    }
    //
#endif


#ifdef USE_RABBITMQ_MESSAGING
    //
    rmq_parameters.rmq_use_rmq = RMQ_USE_RMQ_DEFAULT;
    if ((int_param = settings_get_int(app_prop_settings, "Messaging", "rmq.use_rmq")) != INT_INVALID) {
        rmq_parameters.rmq_use_rmq = int_param;
    }
    //
    strcpy(rmq_parameters.rmq_hostname, RMQ_HOSTNAME_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.hostname", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_hostname, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.hostname %s\n", rmq_parameters.rmq_hostname);
    //
    rmq_parameters.rmq_port = RMQ_PORT_DEFAULT;
    if ((int_param = settings_get_int(app_prop_settings, "Messaging", "rmq.port")) != INT_INVALID) {
        rmq_parameters.rmq_port = int_param;
    }
    printf("Info: property set: [Messaging] rmq.rmq_port %d\n", rmq_parameters.rmq_port);
    //
    strcpy(rmq_parameters.rmq_vhost, RMQ_VHOST_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.vhost", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_vhost, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.vhost %s\n", rmq_parameters.rmq_vhost);
    //
    strcpy(rmq_parameters.rmq_username, RMQ_USERNAME_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.username", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_username, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.username %s\n", rmq_parameters.rmq_username);
    //
    strcpy(rmq_parameters.rmq_password, RMQ_PASSWORD_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.password", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_password, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.password %s\n", rmq_parameters.rmq_password);
    //
    strcpy(rmq_parameters.rmq_exchange, RMQ_EXCHANGE_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.exchange", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_exchange, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.exchange %s\n", rmq_parameters.rmq_exchange);
    //
    strcpy(rmq_parameters.rmq_exchangetype, RMQ_EXCHANGETYPE_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.exchangetype", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_exchangetype, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.exchangetype %s\n", rmq_parameters.rmq_exchangetype);
    //
    strcpy(rmq_parameters.rmq_routingkey, RMQ_ROUTNGKEY_DEFAULT);
    if (settings_get(app_prop_settings, "Messaging", "rmq.routingkey", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(rmq_parameters.rmq_routingkey, out_buf);
    }
    printf("Info: property set: [Messaging] rmq.routingkey %s\n", rmq_parameters.rmq_routingkey);
    //
#endif


    // general initialization

    num_sources_total = 0;
    data_list = NULL;
    num_de_data = 0;
    hypo_list = NULL;
    num_hypocenters = 0;

    int i;

    for (i = 0; i < MAX_NUM_SOURCES; i++) {
        pmem_pick_hf[i] = NULL;
        pmem_pick_raw[i] = NULL;
        pmem_smooth_a_ref[i] = NULL;
        pmem_smooth_t50[i] = NULL;
        pmem_smooth_wwssn_sp_rms[i] = NULL;
        pmem_smooth_brb_hp_rms[i] = NULL;
        pmem_smooth_brb_hp_int_rms[i] = NULL;
        pmem_filter_hf[i] = NULL;
        //pmem_filter_lp[i] = NULL;
        pmem_filter_hp[i] = NULL;
        pmem_filter_hp_int[i] = NULL;
        pmem_filter_bp[i] = NULL;
        pmem_inst_per[i] = NULL;
        pmem_inst_per_dval[i] = NULL;
        //pmem_mwp[i] = NULL;
        pmem_local_t50[i] = NULL;
        pmem_local_a_ref[i] = NULL;
        pmem_local_brb_hp[i] = NULL;
        pmem_local_brb_hp_int_sqr[i] = NULL;
        pmem_local_brb_hp_int_rms[i] = NULL;
        pmem_local_brb_hp_sqr[i] = NULL;
        pmem_local_brb_hp_rms[i] = NULL;
        pmem_local_wwssn_sp[i] = NULL;
        pmem_local_wwssn_sp_sqr[i] = NULL;
        pmem_local_wwssn_sp_rms[i] = NULL;
        last_a_ref[i] = -1.0; // -1.0 indicates a_ref level has dropped below a_ref level for last de measure
        last_a_ref_ellapsed_index_count[i] = -1;
    }

    // must force "local" timezone to be UTC for correct functioning of time_t operations
    setenv("TZ", "UTC", 1);
    /*
     * setenv() on Mac leaks memory intentionally. (refer: man setenv)
     * Valgrind should suppress this on Mac.
     *   https://bugs.kde.org/show_bug.cgi?id=188572
     */

    int nsta;

    // station coordinates, station corrections
    channelParameters = calloc(MAX_NUM_SOURCES, sizeof (ChannelParameters));
    for (nsta = 0; nsta < MAX_NUM_SOURCES; nsta++) {
        initChannelParameters(channelParameters + nsta);
    }
    num_sorted_chan_params = 0;

    // store name
    geogfile_name = geogfile;

    // channel response/gain
    chan_resp = calloc(MAX_NUM_SOURCES, sizeof (ChannelResponse));
    for (nsta = 0; nsta < MAX_NUM_SOURCES; nsta++) {
        chan_resp[nsta].gainfile_checked_time = -1; // set earlier than any possible file mod time to force first check
        chan_resp[nsta].internet_gain_query_checked_time = -1; // set earlier than any possible file mod time to force first check
        chan_resp[nsta].have_gain = 0;
    }
    // store filename
    gainfile_name = gainfile;


    // waveform export memory lists
    waveform_export_miniseed_list = calloc(MAX_NUM_SOURCES, sizeof (MSRecord**));
    num_waveform_export_miniseed_list = calloc(MAX_NUM_SOURCES, sizeof (int));
    for (nsta = 0; nsta < MAX_NUM_SOURCES; nsta++) {
        waveform_export_miniseed_list[nsta] = NULL;
        num_waveform_export_miniseed_list[nsta] = 0;
    }

    int nsource;

    // check and initialize internet pole-zero sources
    if (num_internet_gain_query_host != num_internet_gain_query || num_internet_gain_query_host != num_internet_gain_query_type) {
        fprintf(stderr, "Error - Mismatch in number of -pz-type, -pz-host and -pz parameters.\n");
        return -1;
    }
    //int num_internet_gain = min(min(num_internet_gain_query_host, num_internet_gain_query), num_internet_gain_query_type);
    if (num_internet_station_query_host != num_internet_station_query || num_internet_station_query_host != num_internet_station_query_type) {
        fprintf(stderr, "Error - Mismatch in number of -sta-query-type, -sta-query-host and -sta-query parameters.\n");
        return -1;
    }
    num_internet_gain_sources = 0;
    for (nsource = 0; nsource < num_internet_gain_query; nsource++) {
        if (internet_gain_query_host[nsource] != NULL && internet_gain_query[nsource] != NULL) {
            num_internet_gain_sources++;
        }
    }
    if (num_internet_gain_sources != num_internet_gain_query) {
        if (verbose)
            printf("ERROR: mismatch between -pz-query-host and -pz-query parameters.\n");
        return (-1);
    }

    // check and initialize Internet station parameters sources
    if (num_internet_timeseries_query_hosturl != num_internet_timeseries_query || num_internet_timeseries_query_hosturl != num_internet_timeseries_query_type
            || num_internet_timeseries_query_hosturl != num_internet_timeseries_query_sladdr) {
        fprintf(stderr, "Error - Mismatch in number of -timeseries-query-type, -timeseries-query-host, -timeseries-query and -timeseries-query_sladdr parameters.\n");
        return -1;
    }
    num_internet_timeseries = min(min(num_internet_timeseries_query_hosturl, num_internet_timeseries_query), num_internet_timeseries_query_type);
    int n_int_tseries = 0;
    for (n_int_tseries = 0; n_int_tseries < num_internet_timeseries; n_int_tseries++) {
        strcpy(internetTimeseriesQueryParams[n_int_tseries].query, internet_timeseries_query[n_int_tseries]);
        strcpy(internetTimeseriesQueryParams[n_int_tseries].hosturl, internet_timeseries_query_hosturl[n_int_tseries]);
        strcpy(internetTimeseriesQueryParams[n_int_tseries].sladdr, internet_timeseries_query_sladdr[n_int_tseries]);
        if (strcmp(internet_timeseries_query_type[n_int_tseries], "IRIS_WS_TIMESERIES") == 0) {
            internetTimeseriesQueryParams[n_int_tseries].type = IRIS_WS_TIMESERIES;
        }
    }
    //int num_internet_station = min(min(num_internet_station_query_host, num_internet_station_query), num_internet_station_query_type);
    num_internet_station_sources = 0;
    for (nsource = 0; nsource < num_internet_station_query; nsource++) {
        if (internet_station_query_host[nsource] != NULL && internet_station_query[nsource] != NULL) {
            num_internet_station_sources++;
        }
    }
    if (num_internet_station_sources != num_internet_station_query) {
        if (verbose)
            printf("ERROR: mismatch between -sta-query-host and -sta-query parameters.\n");
        return (-1);
    }


#ifdef TEST_TO_FILE
    // clear pick file
    FILE* fp_picks = fopen("test.picks", "w");
    fclose(fp_picks);
#endif

    if (td_timedomainProcessingReport_init(outpath) < 0)

        return (-1);

    return (0);

}

/***************************************************************************
 * td_set_station_coordinates:
 *
 * Find coordinates for a specified station.
 * Sets lat/lon/elev/azimuth/dip if station found and returns 0.  Returns -1 if station not found.
 ***************************************************************************/
int td_set_station_coordinates(int source_id, char* network, char* station, char* location, char* channel, int data_year, int data_month, int data_mday,
        int verbose, int icheck_ws_station_coords) {

    char net[32];
    char sta[32];
    char loc[32] = "N/A";
    char chan[32] = "N/A";
    double lat = 0.0;
    double lon = 0.0;
    double elev = 0.0;
    double azimuth = -1.0;
    double dip = -999.0;

    // check if geog file modified later than station checked time, and then if station is available
    if (geogfile_name != NULL) {
        // check if file modified
        struct stat file_stats;
        stat(geogfile_name, &file_stats);
        //printf("DEBUG: Check geogfile: %s %s  channelParameters[source_id].geogfile_checked_time %ld <? file_stats.st_mtime %ld\n",
        //        network, station, channelParameters[source_id].geogfile_checked_time, file_stats.st_mtime);
        if (channelParameters[source_id].geogfile_checked_time < file_stats.st_mtime) { // not yet checked or was modified since last check
            //printf("DEBUG: Check geogfile: time: %ld, file: %ld > sta: %ld,  reading from modified file: %s for %s\n",
            //        time(NULL), file_stats.st_mtime, channelParameters[source_id].geogfile_checked_time, geogfile_name, station);
            FILE* fp_geog = fopen(geogfile_name, "r");
            if (fp_geog != NULL) {
                static char line[1024];
                while ((fgets(line, sizeof (line), fp_geog)) != NULL) {
                    // 20160805 AJL - added azimuth and dip
                    int nscan;
                    if ((nscan = sscanf(line, "%s %s %lf %lf %lf %*d %*d %*d %s %s %lf %lf",
                            net, sta, &lat, &lon, &elev, loc, chan, &azimuth, &dip)) >= 5) {
                        //printf("DEBUG: %s[%s] %s[%s] %s[%s] %s[%s] sensor coordinates (%lf %lf %lf), az/dip (%lf %lf) read from geog file: %s\n", net, network, sta, station, loc, location, chan, channel, lat, lon, elev, azimuth, dip, geogfile_name);
                        // 20180205 AJL - added check for "$$" network location channel (from timedomain_processing_data.c->get_next_pick_nll
                        if ((strcmp(net, network) == 0 || strcmp("$$", network) == 0)
                                && (strcmp(sta, station) == 0)) { // found
                            if (nscan <= 5 || ((strcmp(loc, location) == 0 || strcmp("$$", location) == 0)
                                    && (strcmp(chan, channel) == 0 || strcmp("$$", channel) == 0))) { // check loc and chan if available
                                if (verbose)
                                    printf("Info: %s %s %s %s sensor coordinates (%lf %lf %lf), az/dip (%lf %lf) read from geog file: %s\n",
                                        network, station, loc, chan, lat, lon, elev, azimuth, dip, geogfile_name);
                                strcpy(channelParameters[source_id].network, network);
                                strcpy(channelParameters[source_id].station, station);
                                strcpy(channelParameters[source_id].location, location);
                                strcpy(channelParameters[source_id].channel, channel);
                                channelParameters[source_id].lat = lat;
                                channelParameters[source_id].lon = lon;
                                channelParameters[source_id].elev = elev;
                                channelParameters[source_id].azimuth = azimuth;
                                channelParameters[source_id].dip = dip;
                                channelParameters[source_id].have_coords = 1;
                                // add to sorted station params list
                                addChannelParametersToSortedList(channelParameters + source_id, &sorted_chan_params_list, &num_sorted_chan_params);
                                // 20140106 AJL - bug fix
                                //fclose(fp_geog);
                                //return (0);
                                icheck_ws_station_coords = 0; // probably not necessary
                                channelParameters[source_id].internet_station_query_checked_time = file_stats.st_mtime;
                                //printf("DEBUG: Found geogfile: %s %s  lat %f  channelParameters[source_id].lat %f  channelParameters[source_id].lat %f  channelParameters[source_id].lat %f\n",
                                //        network, station, lat, channelParameters[source_id].lat, channelParameters[source_id].lat, channelParameters[source_id].lat);
                                //break;    // 20140107 AJL - no break so last (latest) matching entry in file will be used
                            }
                        }
                    }
                }
                fclose(fp_geog);
            }
            channelParameters[source_id].geogfile_checked_time = file_stats.st_mtime;
        }
    }


    // check if available on internet query service
    if (icheck_ws_station_coords) {
        double lat = channelParameters[source_id].lat;
        double lon = channelParameters[source_id].lon;
        double elev = channelParameters[source_id].elev; // meters
        double azimuth = channelParameters[source_id].azimuth;
        double dip = channelParameters[source_id].dip;
        channelParameters[source_id].internet_station_query_checked_time = time(NULL);
        if (td_get_sta_coords_internet_station_query(source_id, network, station, location, channel, data_year, data_month, data_mday, verbose) == 0) {
            //printf("DEBUG: write geogfile_name %s: %s %s %lf %lf %lf\n", geogfile_name, network, station, channelParameters[source_id].lat, channelParameters[source_id].lon, channelParameters[source_id].elev);
            // found, add to file if values changed
            if (geogfile_name != NULL && (
                    fabs(lat - channelParameters[source_id].lat) > FLT_EPSILON * channelParameters[source_id].lat
                    || fabs(lon - channelParameters[source_id].lon) > FLT_EPSILON * channelParameters[source_id].lon
                    || fabs(elev - channelParameters[source_id].elev) > FLT_EPSILON * channelParameters[source_id].elev
                    || fabs(azimuth - channelParameters[source_id].azimuth) > FLT_EPSILON * channelParameters[source_id].azimuth
                    || fabs(dip - channelParameters[source_id].dip) > FLT_EPSILON * channelParameters[source_id].dip
                    )) {
                //printf("DEBUG: Found internet: %s %s  lat %f  channelParameters[source_id].lat %f  channelParameters[source_id].lat %f  channelParameters[source_id].lat %f\n",
                //        network, station, lat, channelParameters[source_id].lat, channelParameters[source_id].lat, channelParameters[source_id].lat);
                FILE* fp_geog = fopen(geogfile_name, "a");
                if (fp_geog != NULL) {

                    // 20160805 AJL - added azimuth and dip
                    fprintf(fp_geog, "%s %s %lf %lf %lf %d %d %d %s %s %lf %lf\n",
                            network, station, channelParameters[source_id].lat, channelParameters[source_id].lon, channelParameters[source_id].elev,
                            data_year, data_month, data_mday, location, channel, channelParameters[source_id].azimuth, channelParameters[source_id].dip);
                    fclose(fp_geog);
                }
            }
            return (0);
        }
    }


    return (0);
}

/***************************************************************************
 * td_set_station_corrections:
 *
 * Find station correction stats for a specified station.
 * Sets station correction stats if station found and returns 0.  Returns -1 if station not found.
 ***************************************************************************/
int td_set_station_corrections(int source_id, char* network, char* station, char* location, char* channel, int data_year, int data_month, int data_mday,
        int verbose, int icheck_sta_corr) {

    char net[32];
    char sta[32];
    char loc[32];
    char chan[32];
    char phase[32];

    double dist_min, dist_max;
    int num;
    double mean;
    double std_dev;
    double poly[4];

    int year, month, day;

    int ifound = -1;

    // flag that check for available sta corr has been done
    channelParameters[source_id].sta_corr_checked = 1;

    // check if sta corr file modified later than station checked time, and then if station is available
    //if (sta_corr_filename != NULL && strlen(sta_corr_filename) > 0) {
    if (strlen(sta_corr_filename) > 0) {
        // check if file modified
        struct stat file_stats;
        stat(sta_corr_filename, &file_stats);
        //printf("DEBUG: Check sta_corr_file: %s %s  channelParameters[source_id].geogfile_checked_time %ld <? file_stats.st_mtime %ld\n",
        //        network, station, channelParameters[source_id].geogfile_checked_time, file_stats.st_mtime);
        //if (channelParameters[source_id].geogfile_checked_time < file_stats.st_mtime) { // not yet checked or was modified since last check
        //printf("DEBUG: Check sta_corr_file: time: %ld, file: %ld > sta: %ld,  reading from modified file: %s for %s\n",
        //        time(NULL), file_stats.st_mtime, channelParameters[source_id].geogfile_checked_time, sta_corr_filename, station);
        FILE* fp_sta_corr = fopen(sta_corr_filename, "r");
        if (fp_sta_corr != NULL) {
            static char line[1024];
            fgets(line, sizeof (line), fp_sta_corr); // header
            while ((fgets(line, sizeof (line), fp_sta_corr)) != NULL) {
                // network station location channel phase n_data mean std_dev year month day
                // AF SVMA -- BHZ   P   3 0.99 0.556507   2015 01 08
                if (sscanf(line, "%s %s %s %s", net, sta, loc, chan) == 4) {
                    if ((strcmp(net, network) == 0) && (strcmp(sta, station) == 0)
                            && (strcmp(loc, location) == 0) && (strcmp(chan, channel) == 0)) { // found
                        ifound = 0;
                        if (sscanf(line, "%*s %*s %*s %*s %s", phase) == 1) {
                            //if (strcmp(phase, "P") == 0) { // TODO: add other phases (PKPxx ??)
                            if (sscanf(line, "%*s %*s %*s %*s %*s %lf %lf %d %lf %lf %lf %lf %lf %lf %d %d %d",
                                    &dist_min, &dist_max, &num, &mean, &std_dev, &poly[0], &poly[1], &poly[2], &poly[3], &year, &month, &day) == 12) {
                                if (verbose)
                                    printf("Info: %s %s %s %s %s sta corr (%d %lf %lf %.4d %.2d %.2d) read from sta corr file: %s\n",
                                        network, station, location, channel, phase, num, mean, std_dev, year, month, day, sta_corr_filename);
                                int phase_id = phase_id_for_name(phase);
                                if (phase_id >= 0) {
                                    StaCorrections *sta_corr = &(channelParameters[source_id].sta_corr[phase_id]);
                                    sta_corr->dist_min = dist_min;
                                    sta_corr->dist_max = dist_max; // 20150508 AJL added - polynomial only fit and used between dist max and min
                                    sta_corr->num = num;
                                    sta_corr->mean = mean;
                                    sta_corr->std_dev = std_dev;
                                    sta_corr->poly[0] = poly[0];
                                    sta_corr->poly[1] = poly[1];
                                    sta_corr->poly[2] = poly[2];
                                    sta_corr->poly[3] = poly[3];
                                    if (num >= sta_corr_min_num_to_use) {
                                        sta_corr->valid = 1;
                                    }
                                    //break;    // 20140107 AJL - no break so last (latest) matching entry in file will be used
                                }
                            } else {
                                fprintf(stderr, "ERROR: reading sta corr: %s %s %s %s %s sta corr (%d %lf %lf %.4d %.2d %.2d) read from sta corr file: %s\n",
                                        network, station, location, channel, phase, num, mean, std_dev, year, month, day, sta_corr_filename);
                            }
                            //}
                        }
                    }
                }
            }
            fclose(fp_sta_corr);
        } else {
            printf("ERROR: opening station corrections file: %s\n", sta_corr_filename);

            return (-1);
        }
        //channelParameters[source_id].geogfile_checked_time = file_stats.st_mtime;
        //}
    }

    return (ifound);
}

/***************************************************************************
 * td_set_channel_gain:
 *
 * Get and set gain data for a specified channel.
 * Sets gain data if channel response found and returns 0.  Returns -1 if not found or error.
 ***************************************************************************/
int td_set_channel_gain(int source_id, char* network, char* station, char* location, char* channel, int data_year, int data_doy, int verbose, int icheck_ws_gain) {

    char net[32];
    char sta[32];
    char loc[32];
    char chan[32];
    int year = 0;
    int doy = 0;
    double gain = 0.0;
    double freq = 0.0;
    int type = 0;

    // check if gain file modified later than channel checked time, and then if channel is available
    if (gainfile_name != NULL) {
        // check if file modified
        struct stat file_stats;
        stat(gainfile_name, &file_stats);
        //printf("DEBUG: Check gainfile: %s %s %s %s  chan_resp[source_id].gainfile_checked_time %ld <? file_stats.st_mtime %ld\n",
        //        network, station, location, channel, chan_resp[source_id].gainfile_checked_time, file_stats.st_mtime);
        if (chan_resp[source_id].gainfile_checked_time < file_stats.st_mtime) { // not yet checked or was modified since last check
            //printf("DEBUG: Check gainfile: %s time: %ld, file: %ld > sta: %ld,  reading from modified file: %s for %s\n",
            //        gainfile_name, time(NULL), file_stats.st_mtime, chan_resp[source_id].gainfile_checked_time, gainfile_name, station);
            FILE* fp_gain = fopen(gainfile_name, "r");
            if (fp_gain != NULL) {
                static char line[1024];
                while ((fgets(line, sizeof (line), fp_gain)) != NULL) {
                    if (sscanf(line, "%s %s %s %s %d %d %lf %lf %d", net, sta, loc, chan, &year, &doy, &gain, &freq, &type) == 9) {
                        if (gain > FLT_MIN // 20150210 AJL - make sure gain is finite and positive
                                && (strcmp(net, network) == 0) && (strcmp(sta, station) == 0) && (strcmp(loc, location) == 0) && (strcmp(chan, channel) == 0)
                                && (year <= data_year)) { // found - NOTE: effectively no check of year or doy for realtime data
                            if (verbose)
                                printf("Info: %s %s %s %s %d %d  channel gain (%f %f %d) read from gain file: %s\n", net, sta, loc, chan, year, doy, gain, freq, type, gainfile_name);
                            strcpy(chan_resp[source_id].network, net);
                            strcpy(chan_resp[source_id].station, sta);
                            strcpy(chan_resp[source_id].location, loc);
                            strcpy(chan_resp[source_id].channel, chan);
                            chan_resp[source_id].year = year;
                            chan_resp[source_id].doy = doy;
                            chan_resp[source_id].gain = gain;
                            chan_resp[source_id].frequency = freq;
                            chan_resp[source_id].responseType = type;
                            chan_resp[source_id].have_gain = 1;
                            // 20140106 AJL - bug fix
                            //fclose(fp_gain);
                            //return (0);
                            icheck_ws_gain = 0; // probably not necessary
                            chan_resp[source_id].internet_gain_query_checked_time = file_stats.st_mtime;
                            //printf("DEBUG: Found gainfile: %s %s  chan_resp[source_id].gain %f  chan_resp[source_id].frequency %f  chan_resp[source_id].responseType %d\n",
                            //        network, station, chan_resp[source_id].gain, chan_resp[source_id].frequency, chan_resp[source_id].responseType);
                            //break;    // 20140107 AJL - no break so last (latest) matching entry in file will be used
                        }
                    }
                }
                fclose(fp_gain);
            }
            chan_resp[source_id].gainfile_checked_time = file_stats.st_mtime;
        }
    }

    // check if available on internet query service
    if (icheck_ws_gain) {
        double frequency = chan_resp[source_id].frequency;
        double gain = chan_resp[source_id].gain;
        int responseType = chan_resp[source_id].responseType;
        chan_resp[source_id].internet_gain_query_checked_time = time(NULL);
        int data_month;
        int data_day;
        ms_doy2md(data_year, data_doy, &data_month, &data_day);
        if (td_get_gain_internet_query(source_id, network, station, location, channel, data_year, data_month, data_day, verbose) == 0) {
            // found, add to file if values changed
            if (gainfile_name != NULL && (
                    responseType != chan_resp[source_id].responseType
                    || fabs(frequency - chan_resp[source_id].frequency) > FLT_EPSILON * chan_resp[source_id].frequency
                    || fabs(gain - chan_resp[source_id].gain) > FLT_EPSILON * chan_resp[source_id].gain
                    )) {
                //printf("DEBUG: Found internet: %s %s  chan_resp[source_id].gain %f  chan_resp[source_id].frequency %f  chan_resp[source_id].responseType %d\n",
                //        network, station, chan_resp[source_id].gain, chan_resp[source_id].frequency, chan_resp[source_id].responseType);
                FILE* fp_gain = fopen(gainfile_name, "a");
                if (fp_gain != NULL) {

                    fprintf(fp_gain, "%s %s %s %s %d %d %e %f %d\n", network, station, location, channel, data_year, data_doy,
                            chan_resp[source_id].gain, chan_resp[source_id].frequency, chan_resp[source_id].responseType);
                    fclose(fp_gain);
                }
            }
            return (0);
        }
    }

    // not found
    return (-1);

}

/***************************************************************************
 * td_get_gain_internet_gain_query:
 *
 * Get and gain data for a specified channel using an internet query service.
 * Sets gain data if channel response found and returns 0.  Returns -1 if not found or error.
 ***************************************************************************/
int td_get_gain_internet_query(int source_id, char* net, char* sta, char* loc, char* chan, int year, int month, int day, int verbose) {

    int n;
    for (n = 0; n < num_internet_gain_sources; n++) {

        double frequency = MWPD_GAIN_FREQUENCY;
        int responseType;

        // TODO: implement different internet_gain_query_type's

        double gain = -1.0;
        if (strcmp(internet_gain_query_type[n], "IRIS_WS_RESP") == 0) {
            gain = get_gain_internet_gain_query_IRIS_WS_RESP(internet_gain_query_host[n], internet_gain_query[n],
                    net, sta, loc, chan, year, month, day, frequency, &responseType);
        } else if (strcmp(internet_station_query_type[n], "FDSN_WS_STATION") == 0) {
            gain = get_gain_internet_station_query_FDSN_WS_STATION(internet_gain_query_host[n], internet_gain_query[n],
                    net, sta, loc, chan, year, month, day, &frequency, &responseType);
        }

        // 20150210 AJL  if (gain < 0.0)
        if (gain <= FLT_MIN) // 20150210 AJL - make sure gain is finite and positive
            continue;

        if (verbose) {
            fprintf(stdout, "%s: Gain=%e,  Freq=%f", internet_gain_query_host[n], gain, frequency);
            fprintf(stdout, ", Response type=%d (%d=DERIVATIVE, %d=DOUBLE_DERIVATIVE, %d=ERROR)\n", responseType, DERIVATIVE_TYPE, DOUBLE_DERIVATIVE_TYPE, ERROR_TYPE);
        }
        strcpy(chan_resp[source_id].network, net);
        strcpy(chan_resp[source_id].station, sta);
        strcpy(chan_resp[source_id].location, loc);
        strcpy(chan_resp[source_id].channel, chan);
        chan_resp[source_id].gain = gain;
        chan_resp[source_id].frequency = frequency;
        chan_resp[source_id].responseType = responseType;
        chan_resp[source_id].have_gain = 1;

        if (responseType != ERROR_TYPE)

            return (0);

        chan_resp[source_id].have_gain = 0;
        chan_resp[source_id].gain = -1.0;

    }

    return (-1);

}

/***************************************************************************
 * td_get_sta_coords_internet_station_query:
 *
 * Get station data for a specified channel using an internet query service.
 * Sets station data if channel response found and returns 0.  Returns -1 if not found or error.
 ***************************************************************************/
int td_get_sta_coords_internet_station_query(int source_id, char* net, char* sta, char* loc, char* chan, int year, int month, int mday, int verbose) {

    int n;
    for (n = 0; n < num_internet_station_sources; n++) {

        double lat = 0.0;
        double lon = 0.0;
        double elev = 0.0;
        double azimuth = -1.0;
        double dip = -999.0;

        // TODO: implement different internet_station_query_type's

        int return_value = -99;
        if (strcmp(internet_station_query_type[n], "IRIS_WS_STATION") == 0) { // discontinued!
            return_value = get_station_coords_internet_station_query_IRIS_WS_STATION(internet_station_query_host[n], internet_station_query[n],
                    net, sta, year, month, mday, &lat, &lon, &elev);
        } else if (strcmp(internet_station_query_type[n], "FDSN_WS_STATION") == 0) {
            return_value = get_station_coords_internet_station_query_FDSN_WS_STATION(internet_station_query_host[n], internet_station_query[n],
                    net, sta, loc, chan, year, month, mday, &lat, &lon, &elev, &azimuth, &dip);
        }
        // 20140310 AJL - bug fix
        //if (return_value < 0)
        //    return (return_value);
        if (return_value < 0)
            continue;

        if (verbose) {
            fprintf(stdout, "%s: lat=%g lon=%g elev=%g\n", internet_station_query_host[n], lat, lon, elev);
        }
        strcpy(channelParameters[source_id].network, net);
        strcpy(channelParameters[source_id].station, sta);
        strcpy(channelParameters[source_id].location, loc);
        strcpy(channelParameters[source_id].channel, chan);
        channelParameters[source_id].lat = lat;
        channelParameters[source_id].lon = lon;
        channelParameters[source_id].elev = elev;
        channelParameters[source_id].azimuth = azimuth;
        channelParameters[source_id].dip = dip;
        channelParameters[source_id].have_coords = 1;

        // add to sorted station params list
        addChannelParametersToSortedList(channelParameters + source_id, &sorted_chan_params_list, &num_sorted_chan_params);

        return (0);
    }

    return (-1);

}

/***************************************************************************
 * td_process_timedomain_processing_cleanup:
 *
 * Do necessary cleanup.
 ***************************************************************************/
void td_process_timedomain_processing_cleanup() {

    free_TimedomainProcessingDataList(data_list, num_de_data);
    free_HypocenterDescList(&hypo_list, &num_hypocenters);
    free_ChannelParametersList(&sorted_chan_params_list, &num_sorted_chan_params);

    int i;

    for (i = 0; i < MAX_NUM_SOURCES; i++) {

        td_process_free_timedomain_memory(i);
    }

    // unset: force "local" timezone to be UTC for correct function of time_t operations
    unsetenv("TZ");

    for (int nsta = 0; nsta < MAX_NUM_SOURCES; nsta++) {
        freeChannelParameters(channelParameters + nsta);
    }
    free(channelParameters);
    free(chan_resp);

    // waveform export memory lists
    int nsta;
    for (nsta = 0; nsta < MAX_NUM_SOURCES; nsta++) {
        mslist_free(&(waveform_export_miniseed_list[nsta]), &(num_waveform_export_miniseed_list[nsta]));
    }
    if (waveform_export_miniseed_list != NULL)
        free(waveform_export_miniseed_list);

    if (num_waveform_export_miniseed_list != NULL)
        free(num_waveform_export_miniseed_list);


    td_timedomainProcessingReport_cleanup();
}

/***************************************************************************
 * td_process_free_timedomain_memory:
 *
 * Do necessary cleanup.
 ***************************************************************************/
void td_process_free_timedomain_memory(int nsource) {

    free_FilterPicker5_Memory(&(pmem_pick_hf[nsource]));
    free_FilterPicker5_Memory(&(pmem_pick_raw[nsource]));
    free_timedomain_memory(&(pmem_smooth_a_ref[nsource]));
    free_timedomain_memory(&(pmem_smooth_t50[nsource]));
    free_timedomain_memory(&(pmem_smooth_wwssn_sp_rms[nsource]));
    free_timedomain_memory(&(pmem_smooth_brb_hp_rms[nsource]));
    free_timedomain_memory(&(pmem_smooth_brb_hp_int_rms[nsource]));
    free_timedomain_memory(&(pmem_filter_hf[nsource]));
    //free_timedomain_memory(&(pmem_filter_lp[nsource]));
    free_timedomain_memory(&(pmem_filter_hp[nsource]));
    free_timedomain_memory(&(pmem_filter_hp_int[nsource]));
    free_timedomain_memory(&(pmem_filter_bp[nsource]));
    free_timedomain_memory(&(pmem_inst_per[nsource]));
    free_timedomain_memory(&(pmem_inst_per_dval[nsource]));
    //free_timedomain_memory(&(pmem_mwp[nsource]));
    free_timedomain_memory(&(pmem_local_t50[nsource]));
    free_timedomain_memory(&(pmem_local_a_ref[nsource]));
    free_timedomain_memory(&(pmem_local_wwssn_sp[nsource]));
    free_timedomain_memory(&(pmem_local_wwssn_sp_sqr[nsource]));
    free_timedomain_memory(&(pmem_local_wwssn_sp_rms[nsource]));
    free_timedomain_memory(&(pmem_local_brb_hp[nsource]));
    free_timedomain_memory(&(pmem_local_brb_hp_sqr[nsource]));
    free_timedomain_memory(&(pmem_local_brb_hp_rms[nsource]));
    free_timedomain_memory(&(pmem_local_brb_hp_int_sqr[nsource]));
    free_timedomain_memory(&(pmem_local_brb_hp_int_rms[nsource]));

}

/***************************************************************************
 * td_process_free_timedomain_processing_data:
 *
 * Free a all TimedomainProcessingData for a specified source.
 ***************************************************************************/
void td_process_free_timedomain_processing_data(int source_id) {

    int ndata;
    for (ndata = 0; ndata < num_de_data; ndata++) {
        TimedomainProcessingData* deData = data_list[ndata];
        if (deData->source_id == source_id) {

            removeTimedomainProcessingDataFromDataList(deData, &data_list, &num_de_data);
            free_TimedomainProcessingData(deData);
            //data_list[ndata] = NULL; // 20160802 AJL - memory bug fix?
        }
    }

}

/***************************************************************************
 * td_getTimedomainProcessingDataList:
 *
 * returns timedomain-processing data list and its size.
 ***************************************************************************/
void td_getTimedomainProcessingDataList(TimedomainProcessingData*** pdata_list, int* pnum_de_data) {

    *pdata_list = data_list;
    *pnum_de_data = num_de_data;

}

/** check if de data is associated P */

#define DIP_TOLERANCE 2.5 // degrees
#define AZIMUTH_TOLERANCE 5.0 // degrees

int td_doPolarizationAnalysis(TimedomainProcessingData* deData, int ndata) {

    int DEBUG = 0;
    //DEBUG = !strcmp(deData->station, "BALB");
    //DEBUG = !strcmp(deData->station, "ELL");

    deData->polarization.n_analysis_tries++;

    if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: entry: %d %s_%s_%s_%s dt=%f\n",
            ndata + 1, deData->network, deData->station, deData->location, deData->channel, deData->deltaTime);


    // check if analysis already done or not available or not applicable
    if (deData->polarization.status == POL_DONE || deData->polarization.status == POL_NA) {
        if (DEBUG) {
            if (deData->polarization.status == POL_DONE) {
                printf("DEBUG: td_doPolarizationAnalysis: return 1 POL_DONE: az %f+/-%f, dip %f+/-%f, n %d\n",
                        deData->polarization.azimuth, deData->polarization.azimuth_unc, deData->polarization.dip, deData->polarization.dip_unc, deData->polarization.nvalues);
            } else {
                printf("DEBUG: td_doPolarizationAnalysis: return 1 POL_NA: n %d\n", deData->polarization.nvalues);
            }
        }
        return (0);
    }

    int source_id = deData->source_id;
    ChannelParameters* chan_params = channelParameters + source_id;

    // set polarization not available or not applicable
    if (chan_params->channel[2] != 'Z') {
        // only apply polarization analysis for reporting in Z comp channel
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 2  chan_params->channel[2] %c\n", chan_params->channel[2]);
        return (0);
    }
    if (chan_params->channel_set[1] < 0) {
        // channel_set[1] always set last, so must be available to have 3 comp associated
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 3  chan_params->channel_set[1] %d\n", chan_params->channel_set[1]);
        return (0);
    }

    // have vertical comp data with polarization not yet set
    int polarity_z = 1;
    if (fabs(chan_params->dip - (-90.0)) < AZIMUTH_TOLERANCE) { // up
        ;
    } else if (fabs(chan_params->dip - 90.0) < AZIMUTH_TOLERANCE) { // down
        polarity_z = -1;
    } else {
        // not near  vertical
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return not vertical\n");
        return (0);
    }

    // check vertical component is up and near vertical

    // determine x, y channels
    ChannelParameters* chan_params_0 = channelParameters + chan_params->channel_set[0];
    ChannelParameters* chan_params_1 = channelParameters + chan_params->channel_set[1];
    // check both horizontal
    if (fabs(chan_params_0->dip - 0.0) > DIP_TOLERANCE || fabs(chan_params_1->dip - 0.0) > DIP_TOLERANCE) {
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 4: chan_params_0->dip %f\n", chan_params_0->dip);
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 4: chan_params_0->channel %s\n", chan_params_0->channel);
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 4: chan_params_1->channel %s\n", chan_params_1->channel);
        return (0);
    }
    // check azimuth difference near 90 deg
    int source_id_x; // E like orientation
    int source_id_y; // N like orientation
    double az_min, az_max;
    //printf("DEBUG: td_doPolarizationAnalysis: chan_params_0->azimuth %f, chan_params_1->azimuth %f\n", chan_params_0->azimuth, chan_params_1->azimuth);
    if (chan_params_0->azimuth <= chan_params_1->azimuth) {
        az_min = chan_params_0->azimuth;
        az_max = chan_params_1->azimuth;
        source_id_y = chan_params->channel_set[0];
        source_id_x = chan_params->channel_set[1];
    } else {
        az_min = chan_params_1->azimuth;
        az_max = chan_params_0->azimuth;
        source_id_y = chan_params->channel_set[1];
        source_id_x = chan_params->channel_set[0];
    }
    double azimuth_y = az_min;
    if (fabs(az_max - az_min - 90.0) < AZIMUTH_TOLERANCE) { // az_max ~ az_min + 90deg
        ;
    } else if (fabs(az_max - az_min - 270.0) < AZIMUTH_TOLERANCE) { // az_max ~ az_min + 270deg
        // x channel has larger azimuth
        int id_tmp = source_id_x;
        source_id_x = source_id_y;
        source_id_y = id_tmp;
        azimuth_y = az_max;
    } else {
        // not orthogonal
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return not orthogonal\n");
        return (0);
    }

    // check if any required data not available
    if (waveform_export_miniseed_list[source_id] == NULL
            || waveform_export_miniseed_list[chan_params->channel_set[0]] == NULL
            || waveform_export_miniseed_list[chan_params->channel_set[1]] == NULL) {
        deData->polarization.status = POL_NA;
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 5\n");
        return (0);
    }

    double *data_z = NULL;
    double *data_x = NULL;
    double *data_y = NULL;

    // get gain corrections
    // if all zxy gains not available, use 1.0 as approximation (assumes gains are nearly the same for each channel)
    double gain_z = 1.0;
    double gain_x = 1.0;
    double gain_y = 1.0;
    if (chan_resp[source_id].have_gain && chan_resp[source_id_x].have_gain && chan_resp[source_id_y].have_gain) {
        // have all gains, use them
        gain_z = chan_resp[source_id].gain;
        gain_x = chan_resp[source_id_x].gain;
        gain_y = chan_resp[source_id_y].gain;
    }

    // assemble P window data for 3-comp channel set
    double pick_time = (double) deData->t_time_t + deData->t_decsec;
    hptime_t start_time = pick_time * (double) HPTMODULUS;
    hptime_t end_time = (pick_time + polarization_window_start_delay_after_P + polarization_window_length_max + 2.0 * chan_params->deltaTime) * (double) HPTMODULUS;
    hptime_t start_time_return;
    hptime_t end_time_return;
    int num_samples_return;
    int pick_index = 0;
    int z_complete = 0;
    int x_complete = 0;
    int y_complete = 0;
    //int verbose = 99; // DEBUG
    int verbose = 0; // DEBUG
    // Z, vertical component
    data_z = mslist_getDataWindowAsDouble(start_time, end_time, waveform_export_miniseed_list[source_id], num_waveform_export_miniseed_list[source_id],
            verbose, &start_time_return, &end_time_return, &num_samples_return);
    if (end_time_return >= end_time) {
        z_complete = 1;
    }
    // remove offset given by amp at pick
    double pick_amp = data_z[pick_index];
    for (int n = 0; n < num_samples_return; n++) {
        data_z[n] = (data_z[n] - pick_amp) / gain_z;
    }
    int num_samples_max = num_samples_return;
    // X, (N like) component
    data_x = mslist_getDataWindowAsDouble(start_time, end_time, waveform_export_miniseed_list[source_id_x], num_waveform_export_miniseed_list[source_id_x],
            verbose, &start_time_return, &end_time_return, &num_samples_return);
    if (end_time_return >= end_time) {
        x_complete = 1;
    }
    // remove offset given by amp at pick
    pick_amp = data_x[pick_index];
    for (int n = 0; n < num_samples_return; n++) {
        data_x[n] = (data_x[n] - pick_amp) / gain_x;
    }
    num_samples_max = min(num_samples_max, num_samples_return);
    // Y, (E like) component
    data_y = mslist_getDataWindowAsDouble(start_time, end_time, waveform_export_miniseed_list[source_id_y], num_waveform_export_miniseed_list[source_id_y],
            verbose, &start_time_return, &end_time_return, &num_samples_return);
    if (end_time_return >= end_time) {
        y_complete = 1;
    }
    // remove offset given by amp at pick
    pick_amp = data_y[pick_index];
    for (int n = 0; n < num_samples_return; n++) {
        data_y[n] = (data_y[n] - pick_amp) / gain_y;
    }
    num_samples_max = min(num_samples_max, num_samples_return);

    // apply polarization analysis
    // find max number of windows
    double window_width = polarization_window_length_min;
    int nvalues_max = 0;
    while (window_width < polarization_window_length_max + chan_params->deltaTime) {
        nvalues_max++;
        window_width *= 2.0;
    }
    double *az_array = calloc(nvalues_max, sizeof (double));
    double *wt_array = calloc(nvalues_max, sizeof (double));
    // loop through windows
    window_width = polarization_window_length_min;
    double window_width_max = 0.0;
    int istart = (int) polarization_window_start_delay_after_P / chan_params->deltaTime;
    int nvalues_stat = 0;
    //int verbose = DEBUG;
    do {
        int nsamp = 1 + (int) (0.5 + window_width / chan_params->deltaTime);
        if (istart + nsamp < num_samples_max) {
            PolarizationAnalysis *polar_analysis = polarization_covariance(data_z, data_x, data_y, istart, nsamp, polarity_z, azimuth_y, verbose);
            if (DEBUG) printf(
                    "DEBUG: td_doPolarizationAnalysis: %d %s_%s_%s_%s try=%d, dt=%f, w=%f, istart=%d, nsamp=%d: az %f, dip %f, lin %f, plan %f, wt %.2g, zxy %d%d%d %.2g,%.2g,%.2g",
                    ndata + 1, deData->network, deData->station, deData->location, deData->channel,
                    deData->polarization.n_analysis_tries, deData->deltaTime, window_width,
                    istart, nsamp,
                    polar_analysis->azimuth, polar_analysis->dip, polar_analysis->degree_of_linearity, polar_analysis->degree_of_planarity,
                    polar_analysis->mean_vect_amp, z_complete, x_complete, y_complete, gain_z, gain_x, gain_y);
            if (polar_analysis->degree_of_linearity >= POLARIZATION_MIN_DEG_LINEARITY_TO_USE) {
                az_array[nvalues_stat] = polar_analysis->azimuth;
                wt_array[nvalues_stat] = polar_analysis->mean_vect_amp;
                nvalues_stat++;
                if (DEBUG) printf(" ***");
            }
            free(polar_analysis);
            if (DEBUG) printf("\n");
            window_width_max = window_width;
        }
        window_width *= 2.0;
    } while (window_width < polarization_window_length_max + chan_params->deltaTime);

    if (nvalues_stat > 1) {

        // unwrap stat values based on first value
        double az_ref = az_array[0];
        for (int i = 1; i < nvalues_stat; i++) {
            if (az_array[i] - az_ref > 180.0) {
                az_array[i] -= 360.0;
            } else if (az_array[i] - az_ref < -180.0) {
                az_array[i] += 360.0;
            }
        }

        // get stats
        double az_mean;
        // unweighted
        //double az_var = calcVarianceMean(az_array, nvalues_stat, &az_mean);
        // weighted
        double az_var = calcVarianceMeanWeighted(az_array, wt_array, nvalues_stat, &az_mean);
        double az_std_dev = 0.0;
        if (az_var > FLT_MIN) {
            az_std_dev = sqrt(az_var);
        }
        // 20180112 AJL - bug fix, moved outside if statement
        //free(az_array);
        //free(wt_array);
        // constrain az to 0-360deg
        if (az_mean >= 360.0) {
            az_mean -= 360.0;
        } else
            if (az_mean < 0.0) {
            az_mean += 360.0;
        }
        if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: %d %s_%s_%s_%s nd=%d az=%f +/- %f\n",
                ndata + 1, deData->network, deData->station, deData->location, deData->channel, nvalues_stat, az_mean, az_std_dev);

        // successfully set polarization
        if (deData->polarization.n_analysis_tries > POLARIZATION_MAX_NUM_ANALYSIS_TRIES
                || window_width_max >= polarization_window_length_max || (z_complete && x_complete && y_complete)) {
            // have all needed data, polarization analysis done
            deData->polarization.status = POL_DONE;
        } else {
            deData->polarization.status = POL_SET;
        }
        deData->polarization.nvalues = nvalues_stat;
        deData->polarization.azimuth = az_mean;
        deData->polarization.azimuth_unc = az_std_dev;

    } else if (deData->polarization.n_analysis_tries > POLARIZATION_MAX_NUM_ANALYSIS_TRIES
            || window_width_max >= polarization_window_length_max || (z_complete && x_complete && y_complete)) {
        // have all needed data, polarization analysis not available
        deData->polarization.status = POL_NA;
    }

    // 20180112 AJL - bug fix, moved from if statement above
    free(az_array);
    free(wt_array);

    free(data_z);
    free(data_x);
    free(data_y);

    //if (DEBUG) printf("DEBUG: td_doPolarizationAnalysis: return 99\n");

    return (1);
}

