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
#include <errno.h>
#include <sys/resource.h>
#include <dirent.h>
#include <ftw.h>
#include <fcntl.h>

#define EXTERN_MODE

//#include "../geometry/geometry.h"
//#include "../alomax_matrix/alomax_matrix.h"
//#include "../matrix_statistics/matrix_statistics.h"
#include "../statistics/statistics.h"
#include "../picker/PickData.h"
#include "timedomain_processing_data.h"
#include "ttimes.h"
#include "location.h"
#include "timedomain_processing.h"
#include "timedomain_processing_report.h"
#include "timedomain_memory.h"
#include "timedomain_processes.h"
#include "timedomain_filter.h"
#include "../picker/FilterPicker5_Memory.h"
#include "../picker/FilterPicker5.h"
#include "../feregion/feregion.h"
#include "../miniseed_utils/miniseed_utils.h"
#include "loc2xml.h"
#include "../response/response_lib.h"

// 20160513 - following line added per suggestion by mehmety@boun.edu.tr
#define _XOPEN_SOURCE 500

#define IMIN(a,b) (a<b?a:b)

#define DEBUG 1

#undef USE_MWP_LEVEL_ARRAY



// locals set from from properties file
//
static int num_phase_names_use = -1;
static char phase_names_use[MAX_NUM_TTIME_PHASES][32];
static char channel_names_use[MAX_NUM_TTIME_PHASES][128];
static char time_delay_use[MAX_NUM_TTIME_PHASES][64]; // 20150410 AJL - added
//
static double depth_step_full;
static double depth_min_full; // range defines grid cell limits
static double depth_max_full; // range defines grid cell limits
static double lat_step_full;
static double lat_min_full; // range defines grid cell limits
static double lat_max_full; // range defines grid cell limits
static double lon_step_smallest_full; // lon step is nominal for lat = 0 (on equator), will  be adjusted by 1/cos(lat) in assoc/location routine
static double lon_min_full; // range defines grid cell limits
static double lon_max_full; // range defines grid cell limits
//
static double nomimal_critical_oct_tree_node_size; // nominal (approximate) size of oct tree node that must be reached before accepting a location
//    may be reduced automatically in location.c -> octtreeGlobalAssociationLocation()
static double min_critical_oct_tree_node_size; // minimum size of nomimal_critical_oct_tree_node_size
static double nominal_oct_tree_min_node_size; // nominal (approximate) minimum size of oct tree node to be used for association/location
//    may be reduced automatically in location.c -> octtreeGlobalAssociationLocation()
//
static double min_weight_sum_assoc;
static double min_time_delay_between_picks_for_location;
//
static double gap_primary_critical;
static double gap_secondary_critical;

//
// END - locals set from from properties file


static time_t report_time_max = LONG_MAX;
static time_t earliest_time = LONG_MAX;
static char tmp_str[STANDARD_STRLEN];
static char tmp_str_2[STANDARD_STRLEN];
static HypocenterDesc existing_hypo_desc;

static char outname[1024];
static char xmlWriterUri[1024];

#define FEREGION_STR_SIZE 4096
static char feregion_str[STANDARD_STRLEN];

#define MAX_NUM_HYPO 100
static HypocenterDesc **hyp_assoc_loc = NULL;
static int num_hypocenters_associated = 0;
static int hyp_persistent[MAX_NUM_HYPO]; // indices of persistent hypocenters

// reference full search oct-tree for alignment of reduced search volume
static Tree3D *pOctTree = NULL;


// 20140627 AJL - following declarations added to support persistent hypocenters
static float **assoc_scatter_sample = NULL;
static int n_alloc_scatter_sample[MAX_NUM_HYPO];
static int n_assoc_scatter_sample[MAX_NUM_HYPO];
static double global_max_nassociated_P_lat_lon[MAX_NUM_HYPO];


static int numPhaseTypesUse = -1;
static int phaseTypesUse[MAX_NUM_TTIME_PHASES]; // int phase code to use
static char channelNamesUse[MAX_NUM_TTIME_PHASES][128]; // string containing allowable channel names to use
static double timeDelayUse[MAX_NUM_TTIME_PHASES][2]; // min [0] and max [1] allowable time delays after previously defined phase to use // 20150410 AJL - added
// reference phase ttime error for phase weighting
// 20130307 AJL - added to allow use of non P phases for location
static double reference_phase_ttime_error = DBL_MAX;


// 20141212 AJL - made static global
static int nstaHasBeenActive = -1; // number of stations for which data has been received in past
static int nstaIsActive = -1; // number of stations for which data has been received and data_latency < report_interval


static int n_open_files = 0;
static int max_n_open_files = 0;

FILE* fopen_counter(const char *filename, const char *mode) {

    FILE* fptr = fopen(filename, mode);
    if (fptr == NULL) {
        printf("Info: fopen_counter(): cannot open file: %s\n", filename);
    } else {
        n_open_files++;
        if (n_open_files > max_n_open_files)
            max_n_open_files = n_open_files;
    }
    //printf("DEBUG: fopen_counter(): nfiles open: %d\n", n_open_files);

    return (fptr);

}

int fclose_counter(FILE *stream) {

    int retval = fclose(stream);
    if (retval != 0) {
        printf("ERROR: fclose_counter(): closing output file.\n");
    } else {
        n_open_files--;
    }
    //printf("DEBUG: fclose_counter(): nfiles open: %d\n", n_open_files);

    return (retval);

}

int ignoreData(TimedomainProcessingData* deData) {
    // 20130128 AJL - use flag_snr_brb_int_too_low to allow mwp, mwpd, etc., but do not use for ignore tests (e.g. ignore determined by flag_snr_brb_too_low)
    //return (
    //        deData->flag_clipped || deData->flag_non_contiguous || deData->flag_a_ref_not_ok
    //        || (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low && deData->flag_snr_brb_int_too_low)
    //        );
    // 20131022 AJL - try using all picks for location, regardless of HF S/N
    //return (
    //        deData->flag_clipped || deData->flag_non_contiguous || deData->flag_a_ref_not_ok
    //        || (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low)
    //        );
    return (
            deData->flag_clipped || deData->flag_non_contiguous || (!USE_AREF_NOT_OK_PICKS_FOR_LOCATION && deData->flag_a_ref_not_ok)
            // 20150402 AJL  || (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low && !(USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION && deData->is_associated))
            || (!USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION && deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low)
            );
    // END 20131022 AJL - try using all picks for location, regardless of HF S/N

}

/** format two integer second times into pretty hour-min-sec string */

char *timeDiff2string(time_t later_time, time_t earlier_time, char* tdiff_string) {

    double tdiff = difftime(later_time, earlier_time);

    int tdiff_tenths_sec = (int) (0.5 + 10.0 * tdiff);
    int tdiff_hour = 0;
    int tdiff_min = 0;
    while (tdiff_tenths_sec >= 36000) {
        tdiff_tenths_sec -= 36000;
        tdiff_hour++;
    }
    while (tdiff_tenths_sec >= 600) {
        tdiff_tenths_sec -= 600;
        tdiff_min++;
    }
    if (tdiff_hour > 0)
        sprintf(tdiff_string, "%dh %2.2dm %2.2ds", tdiff_hour, tdiff_min, tdiff_tenths_sec / 10);
    else if (tdiff_min > 0)
        sprintf(tdiff_string, "%dm %2.2ds", tdiff_min, tdiff_tenths_sec / 10);
    else if (tdiff_tenths_sec > 0)
        sprintf(tdiff_string, "%ds", tdiff_tenths_sec / 10);
    else
        sprintf(tdiff_string, "0s");

    return (tdiff_string);

}

/** format a time_t time into time string with integer seconds */

char *time2string(double time_dec_sec, char* str) {

    // round origin time to nearest second
    time_t itime_sec = (time_t) (0.5 + time_dec_sec); // 20110930 AJL - modified to avoid times that round to 60s
    struct tm* tm_time = gmtime(&itime_sec);

    sprintf(str, "%4.4d.%2.2d.%2.2d-%2.2d:%2.2d:%2.2d",
            tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, tm_time->tm_sec);

    return (str);
}


/** format a double time (equivalent to time_t time with decimal seconds) into a time string with decimal seconds */

#define DEFAULT_TIME_FORMAT 0
#define IRIS_WS_TIME_FORMAT 1
#define COMMA_DELIMTED_TIME_FORMAT 2

char *timeDecSec2string(double time_dec_sec, char* str, int iformat) {

    time_t itime_sec = (time_t) time_dec_sec;
    struct tm* tm_time = gmtime(&itime_sec);

    double dec_sec = time_dec_sec - itime_sec;

    if (iformat == IRIS_WS_TIME_FORMAT) {
        sprintf(str, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    } else if (iformat == COMMA_DELIMTED_TIME_FORMAT) {
        sprintf(str, "%4.4d,%2.2d,%2.2d,%2.2d,%2.2d,%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    } else {
        sprintf(str, "%4.4d.%2.2d.%2.2d-%2.2d:%2.2d:%05.2f",
                tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);
    }

    return (str);
}

/** convert a time string into a double time (equivalent to time_t time with decimal seconds) */

double string2timeDecSec(char *str) {

    //printf("\nDEBUG: string2timeDecSec %s\n", str);

    struct tm tm_time;

    double dec_sec;
    sscanf(str, "%4d.%2d.%2d-%d:%d:%lf", &(tm_time.tm_year), &(tm_time.tm_mon), &(tm_time.tm_mday), &(tm_time.tm_hour), &(tm_time.tm_min), &dec_sec);
    tm_time.tm_year -= 1900;
    tm_time.tm_mon -= 1;
    tm_time.tm_sec = (int) dec_sec;

    //printf("DEBUG: string2timeDecSec %.4d.%.2d.%.2d-%.2d:%.2d:%.2d\n", tm_time.tm_year + 1900, tm_time.tm_mon + 1, tm_time.tm_mday, tm_time.tm_hour, tm_time.tm_min, tm_time.tm_sec);

    time_t time_gmt = timegm(&tm_time);
    //printf("DEBUG: time_gmt %ld\n", time_gmt);

    double dec_time = (double) time_gmt + dec_sec - (int) dec_sec;
    //printf("DEBUG: dec_time %lf\n", dec_time);

    return (dec_time);
}

/** convert a double time (equivalent to time_t time with decimal seconds) into a libmseed hptime_t */

hptime_t timeDecSec2hptime(double time_dec_sec) {

    time_t itime_sec = (time_t) time_dec_sec;
    struct tm* tm_time = gmtime(&itime_sec);

    double dec_sec = time_dec_sec - itime_sec;
    int year = tm_time->tm_year + 1900;
    int month = tm_time->tm_mon + 1;
    int mday = tm_time->tm_mday;
    int jday;
    ms_md2doy(year, month, mday, &jday);
    return (ms_time2hptime(year, jday, tm_time->tm_hour, tm_time->tm_min, tm_time->tm_sec, (int) (0.5 + (dec_sec * (double) 1000000))));

}

/** convert a double time (equivalent to time_t time with decimal seconds) into component date/time elements */

void hptime2dateTimeComponents(hptime_t hptime, int *pyear, int *pjday, int *phour, int *pmin, double *psec) {

    BTime btime;
    ms_hptime2btime(hptime, &btime);

    *pyear = btime.year;
    *pjday = btime.day;
    *phour = btime.hour;
    *pmin = btime.min;
    *psec = (double) btime.sec + ((double) btime.fract / (double) 10000);

}

/** helper function to remove a file or directory */

int remove_fn(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf) {

    //printf("DEBUG: remove_fn(): remove file: %s\n", fpath);
    int ireturn = remove(fpath);
    if (ireturn) {
        printf("ERROR: remove_fn(): return value = %d, path = %s\n", ireturn, fpath);
    }
    return ireturn;
}

/** function to remove a waveform directory and all files within the directory if the directory name is an event id older
 *      the maximum waveform export file age (waveform_export_file_archive_age_max) relative to the end of the current report interval
 */

static char filepath[STANDARD_STRLEN];

int remove_waveform_export_directories(const char *dirpath) {

    //printf("DEBUG: remove_waveform_export_directories(): root: %s\n", dirpath);

    int ireturn = 0;

    DIR* pDir = opendir(dirpath);
    struct dirent *pFile = NULL;

    while ((pFile = readdir(pDir))) {

        sprintf(filepath, "%s/%s", dirpath, pFile->d_name);
        //printf("DEBUG: remove_waveform_export_directories(): processing %s\n", filepath);

        if (!strcmp(pFile->d_name, ".") || !strcmp(pFile->d_name, "..")) {
            continue;
        }

        // check if old event directory
        time_t archive_earliesttime = report_time_max - waveform_export_file_archive_age_max;
        long hypo_unique_id = atol(pFile->d_name);
        //printf("DEBUG: remove_waveform_export_directories(): hypo_unique_id/1000 %ld, archive_earliesttime %ld, diff %ld\n",
        //        (hypo_unique_id / 1000), archive_earliesttime, (hypo_unique_id / 1000) - archive_earliesttime);
        if (hypo_unique_id > 0 && hypo_unique_id != LONG_MIN && hypo_unique_id != LONG_MAX && (hypo_unique_id / 1000) < archive_earliesttime) {
            ireturn = nftw(filepath, remove_fn, 16, FTW_DEPTH);
            if (ireturn) {
                printf("ERROR: remove_waveform_export_directories(): removing files: return value = %d, path = %s\n", ireturn, filepath);
            }
        }

    }

    closedir(pDir);

    return (ireturn);
}


/** create links if possible for web service to display plot of raw and filtered channel data */

static char hf_str[STANDARD_STRLEN];
static char brb_str[STANDARD_STRLEN];
static char str_channel_linked[STANDARD_STRLEN];

char *create_channel_links(char *network, char *station, char *location, char *channel, int pick_stream, char* stream_name, int n_int_tseries,
        double start_time, int duration, char* str_channel, char* str_stream) {

    //char* stream_name = pick_stream_name(deData);
    strcpy(str_stream, stream_name);
    sprintf(str_channel, "%s_%s_%s_%s", network, station, location, channel);

    // check if internet timeseries info available for this data
    //int n_int_tseries = deData->n_int_tseries;
    if (n_int_tseries >= 0) {

        if (internetTimeseriesQueryParams[n_int_tseries].type == IRIS_WS_TIMESERIES) {

            // channel link
            // assume if IRIS_WS_TIMESERIES, then IRIS DMC MetaData Aggregator available
            sprintf(str_channel_linked, "<a href=\"%s/%s/%s\" "
                    "title=\"Show station meta-data\" target=%s>%s_%s</a>_%s_%s",
                    IRIS_METADATA_BASE_URL, network, station, str_channel, network, station, location, channel
                    );

            // stream link
            /*
            // pick towards beginning of plot window
            int before_pick_lead = 60 * (time_max - deData->t_time_t) / (5 * 60); // 1 min lead time for each 5 min available after pick
            before_pick_lead = 60 * (before_pick_lead / 60); // round to minute
            before_pick_lead += 60;
            time_t start_time = deData->t_time_t - before_pick_lead;
            int duration = time_max - start_time + 60;
            if (duration > 3600)
                duration = 3600;
             */
            //struct tm* tm_gmt = gmtime(&start_time);
            sprintf(brb_str, "<a href=\"%s/%s?net=%s&sta=%s&loc=%s&cha=%s&start=%s&dur=%d&output=plot\" "
                    "title=\"View %s stream\" target=%s_%s>%s</a>",
                    internetTimeseriesQueryParams[n_int_tseries].hosturl, internetTimeseriesQueryParams[n_int_tseries].query,
                    network, station, location, channel,
                    timeDecSec2string(start_time, tmp_str, IRIS_WS_TIME_FORMAT), duration,
                    STREAM_RAW_NAME, str_channel, STREAM_RAW_NAME, STREAM_RAW_NAME
                    );
            sprintf(hf_str, "<a href=\"%s/%s?net=%s&sta=%s&loc=%s&cha=%s&start=%s&dur=%d&output=plot&bpfilter=%s\" "
                    "title=\"View %s stream\" target=%s_%s>%s</a>",
                    internetTimeseriesQueryParams[n_int_tseries].hosturl, internetTimeseriesQueryParams[n_int_tseries].query,
                    network, station, location, channel,
                    timeDecSec2string(start_time, tmp_str, IRIS_WS_TIME_FORMAT), duration, STREAM_HF_BPFILTER,
                    STREAM_HF_NAME, str_channel, STREAM_HF_NAME, STREAM_HF_NAME
                    );
            if (pick_stream == STREAM_HF) {
                sprintf(str_stream, "%s&#8212;<font size=\"-2\">(<i>%s</i>)</font>", hf_str, brb_str);
            } else if (pick_stream == STREAM_RAW) {
                sprintf(str_stream, "%s&#8212;<font size=\"-2\">(<i>%s</i>)</font>", brb_str, hf_str);
            } else {
                sprintf(str_stream, "%s&#8212;%s", brb_str, hf_str);
            }
            strcpy(str_channel, str_channel_linked);

        }
    }

    return (str_channel);
}



static char *qualityCodeBackgroundColor[] = {
    // A B C D
    "\"#00FF00\"", "\"#FFFF00\"", "\"#FFA500\"", "\"#FF0000\"",
};
static char quality_code_list[] = "ABCD";

// 20160307 AJL - added

char *getQualityCodeBackgroundColorStringHtml(char *quality_code) {

    int ndx = strstr(quality_code_list, quality_code) - quality_code_list;
    if (ndx >= 0 && ndx < strlen(quality_code_list)) {
        return (qualityCodeBackgroundColor[ndx]);
    }

    return ("\"#FFFFFF\"");

}


/*
HYPO_COLOR[0]=31/31/191
HYPO_COLOR[1]=159/31/159
HYPO_COLOR[2]=31/127/0
HYPO_COLOR[3]=159/95/0
 */
static char *hypoBackgroundColor[] = {
    "bgcolor=\"#CCCCFF\"", "bgcolor=\"#FFCCFF\"", "bgcolor=\"#99CC99\"", "bgcolor=\"#CCCC99\"",
};
static int hypoBackgroundColorModulo = 4;
static char levelStringHtml[WARNING_LEVEL_STRING_LEN];
static char eventBackgroundColorString[WARNING_LEVEL_STRING_LEN];
static char hypoDataString[STANDARD_STRLEN];
static char hypoMessageHtmlString[STANDARD_STRLEN];

void setLevelString(int numLevelMax, statistic_level *plevelStatistics, char *levelString, int min_num_levels,
        double red_cutoff, double yellow_cutoff, double invalid, int colors_show) {
    //
    if (numLevelMax < min_num_levels) {
        strcpy(levelString, "GREY");
    } else if (!colors_show) {
        strcpy(levelString, "WHITE");
    } else if (plevelStatistics->centralValue >= red_cutoff) {
        strcpy(levelString, "RED");
    } else if (plevelStatistics->centralValue >= yellow_cutoff) {
        strcpy(levelString, "YELLOW");
    } else if (plevelStatistics->centralValue != invalid) {
        strcpy(levelString, "GREEN");
    } else {
        strcpy(levelString, "GREY");
    }

}

int setLevelStringHtml(int numLevelMax, double centralValue, char *prefix, char *rowBackground, char *colorString, int min_num_levels,
        double red_cutoff, double yellow_cutoff, double invalid, int colors_show) {

    int highlight = 0;
    //
    if (numLevelMax < min_num_levels) {
        sprintf(colorString, "%s", rowBackground);
    } else if (!colors_show) {
        //sprintf(colorString, "%s%s", prefix, "\"#FFFFFF\"");
        sprintf(colorString, "%s", rowBackground);
        highlight = 1;
    } else if (centralValue >= red_cutoff) {
        sprintf(colorString, "%s%s", prefix, "\"#FFAAAA\"");
        highlight = 1;
    } else if (centralValue >= yellow_cutoff) {
        sprintf(colorString, "%s%s", prefix, "\"#FFFFBB\"");
        highlight = 1;
    } else if (centralValue != invalid) {
        sprintf(colorString, "%s%s", prefix, "\"#CCFFCC\"");
        highlight = 1;
    } else {
        sprintf(colorString, "%s", rowBackground);
    }

    return (highlight);

}

int setEventBackgroundColorStringHtml(int num_levels, double level, char *prefix, char *colorString, char *qualityCode,
        int min_num_levels, double level_small, double level_intermediate, double level_high, char *loc_quality_acceptable) {


    // greyscale
    int col_small_low = 207;
    int col_high = 255;
    int col_base = 199;
    // grey
    int col_grey = 199;

    //if (!warning_colors_show) { // yellow-scale
    // grey yellow
    int col_small_lowr = 255;
    int col_small_lowg = 255;
    int col_small_lowb = 239;
    // yellow
    int col_intermed_lowr = 255;
    int col_intermed_lowg = 255;
    int col_intermed_lowb = 127;
    // red
    int col_large_r = 255;
    int col_large_g = 127;
    int col_large_b = 127;
    col_base = 239;
    col_grey = 199;
    //}

    //
    if (strstr(loc_quality_acceptable, qualityCode) == NULL) { // 20160307 AJL - added
        sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_grey, col_grey, col_grey);
        return (0);
    } else if (num_levels < min_num_levels) {
        sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_base, col_base, col_base);
        return (0);
    } else if (level <= level_small) {
        sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_base, col_base, col_base);
    } else if (level >= level_high) {
        if (warning_colors_show) { // greyscale
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_high, col_high, col_high);
        } else { // large
            //sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, 255, 255, col_high);
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_large_r, col_large_g, col_large_b);
        }
    } else if (level >= level_intermediate) {
        if (warning_colors_show) { // greyscale
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col_high, col_high, col_high);
        } else { // intermediate-large scale
            double factor = (level - level_intermediate) / (level_high - level_intermediate);
            int colr = (int) ((double) col_intermed_lowr + (double) (col_large_r - col_intermed_lowr) * factor * factor);
            int colg = (int) ((double) col_intermed_lowg + (double) (col_large_g - col_intermed_lowg) * factor * factor);
            int colb = (int) ((double) col_intermed_lowb + (double) (col_large_b - col_intermed_lowb) * factor * factor);
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, colr, colg, colb);
        }
    } else {
        double factor = (level - level_small) / (level_intermediate - level_small);
        if (warning_colors_show) { // greyscale
            int col = (int) ((double) col_small_low + (double) (col_high - col_small_low) * factor * factor);
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, col, col, col);
        } else { // small-intermediate scale
            int colr = (int) ((double) col_small_lowr + (double) (col_intermed_lowr - col_small_lowr) * factor * factor);
            int colg = (int) ((double) col_small_lowg + (double) (col_intermed_lowg - col_small_lowg) * factor * factor);
            int colb = (int) ((double) col_small_lowb + (double) (col_intermed_lowb - col_small_lowb) * factor * factor);
            sprintf(colorString, "%s\"#%02X%02X%02X\"", prefix, colr, colg, colb);
        }
    }
    return (1);

}

/** create compact hypo information string */

void create_compact_hypo_info_string(HypocenterDesc* phypo, char* hypo_info_string) {

    // set FE-region
    feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE);
    sprintf(hypo_info_string, "%s  %s  lat:%.2f lon:%.2f  depth:%.0fkm  Q:%s mb:%.1f[%d]  Mwp:%.1f[%d]  T0:%.0fs[%d]  Mwpd:%.1f[%d]  TdT50Ex:%.1f",
            feregion_str,
            time2string(phypo->otime, tmp_str),
            phypo->lat, phypo->lon, phypo->depth, phypo->qualityIndicators.quality_code,
            phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mbLevelStatistics.centralValue : MB_INVALID, phypo->mbLevelStatistics.numLevel,
            phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpLevelStatistics.centralValue : MWP_INVALID, phypo->mwpLevelStatistics.numLevel,
            phypo->t0LevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->t0LevelStatistics.centralValue : T0_INVALID, phypo->t0LevelStatistics.numLevel,
            phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mwpLevelStatistics.centralValue >= MIN_MWP_FOR_VALID_MWPD
            ? phypo->mwpdLevelStatistics.centralValue : MWPD_INVALID, phypo->mwpdLevelStatistics.numLevel,
            phypo->tdT50ExLevelStatistics.centralValue
            // 20140306 AJL  , phypo->warningLevelString
            );
}

static char map_fileroot[STANDARD_STRLEN] = "\0";
static char map_link_str[STANDARD_STRLEN] = "\0";
static char js_file_str[STANDARD_STRLEN] = "\0";
static char js_link_str[STANDARD_STRLEN] = "\0";
static char image_tag_str[STANDARD_STRLEN] = "\0";
static char event_url_str[STANDARD_STRLEN] = "\0";
static char event_link_str[STANDARD_STRLEN] = "\0";
static char hypo_description_str[STANDARD_STRLEN] = "\0";
static char hypo_description_str2[STANDARD_STRLEN] = "\0";

/** create event page link */

void create_event_link(char* link_root, long event_id, char* event_url, char* event_link, char* feregion_str) {

    sprintf(event_url, "%s/events/hypo.%ld.html", link_root, event_id);
    sprintf(event_link, "<a href=\"%s\" title=\"View event information page - %s\" target=event_%ld>",
            event_url, feregion_str, event_id);

}

/** create mechanism image tag */

void create_mech_image_tag(char* link_root, long event_id, char* image_tag, char* feregion_str) {

    create_event_link(link_root, event_id, event_url_str, event_link_str, feregion_str);
    sprintf(image_tag, "%s<img src=\"%s/events/hypo.%ld.mechanism.thumb.jpg\" height=\"17\"></a>",
            event_link_str, link_root, event_id);

}

#define MAP_LINK_GLOBAL_ZOOM 2
#define MAP_LINK_MED_ZOOM 6
#define MAP_LINK_BIG_ZOOM MAP_LINK_MED_ZOOM
#define MAP_LINK_DEFAULT_ZOOM MAP_LINK_MED_ZOOM

/** create map link */

void create_map_link(char* link_root, long event_id, char* map_link, int zoom_level) {

    if (event_id > 0) {
        sprintf(map_fileroot, "hypo.%ld", event_id);
    } else {
        strcpy(map_fileroot, "station");
    }


    if (zoom_level == MAP_LINK_DEFAULT_ZOOM || zoom_level == MAP_LINK_GLOBAL_ZOOM) {
        sprintf(map_link, "%s/%s%s.map.html", link_root, event_id > 0 ? "events/" : "", map_fileroot);
    } else {
        sprintf(map_link, "%s/%s%s.map.zoom%d.html", link_root, event_id > 0 ? "events/" : "", map_fileroot, zoom_level);
    }

}

/** create map javascript file name */

void create_map_javascript_file_link(char* link_root, long event_id, char* js_file, char* js_link, int zoom_level) {

    if (event_id > 0) {
        sprintf(map_fileroot, "hypo.%ld", event_id);
    } else {
        strcpy(map_fileroot, "station");
    }

    if (zoom_level == MAP_LINK_DEFAULT_ZOOM || zoom_level == MAP_LINK_GLOBAL_ZOOM) {
        sprintf(js_file, "./%s.map.js", map_fileroot);
        sprintf(js_link, "%s/%s%s.map.js", link_root, event_id > 0 ? "events/" : "", map_fileroot);
    } else {
        sprintf(js_file, "./%s.map.zoom%d.js", map_fileroot, zoom_level);
        sprintf(js_link, "%s/%s%s.map.zoom%d.js", link_root, event_id > 0 ? "events/" : "", map_fileroot, zoom_level);
    }

}

/** function to convert lat/long to rectangular km coord */

int latlon2rect_simple(double dlat, double dlong, double* pxrect, double* pyrect, double map_orig_lat, double map_orig_long) {

    double xtemp, ytemp;

    xtemp = dlong - map_orig_long;
    if (xtemp > 180.0)
        xtemp -= 360.0;
    if (xtemp < -180.0)
        xtemp += 360.0;
    xtemp = xtemp * DEG2KM * cos(DE2RA * dlat);
    ytemp = (dlat - map_orig_lat) * DEG2KM;
    *pxrect = xtemp;
    *pyrect = ytemp;

    return (0);

}

/** function to convert rectangular km coord to lat/long */

int rect2latlon_simple(double xrect, double yrect, double* pdlat, double* pdlong, double map_orig_lat, double map_orig_long) {

    double xtemp, ytemp;

    xtemp = xrect;
    ytemp = yrect;
    *pdlat = map_orig_lat + ytemp / DEG2KM;
    *pdlong = map_orig_long + xtemp / (DEG2KM * cos(DE2RA * *pdlat));
    // prevent longitude outside of -180 -> 180 deg range
    if (*pdlong < -180.0)
        *pdlong += 360.0;
    else if (*pdlong > 180.0)
        *pdlong -= 360.0;

    return (0);

}

/** create Google maps page with Google Maps JavaScript API https://developers.google.com/maps/documentation/javascript/reference */
// 20150909 AJL - added plotting of un-associated stations, station parameters and station health map

#define MAX_NUM_SCATTTER_POINTS_MAP 1000

void create_map_html_page(char *outnameroot, HypocenterDesc* phypo_map, time_t time_min, time_t time_max, int zoom_level) {

    if (phypo_map != NULL) { // creating map for associated/located event
        create_map_link(outnameroot, phypo_map->unique_id, map_link_str, zoom_level);
        create_compact_hypo_info_string(phypo_map, hypo_description_str);
        create_map_javascript_file_link(outnameroot, phypo_map->unique_id, js_file_str, js_link_str, zoom_level);
    } else { // creating map for station status/health
        create_map_link(outnameroot, -1, map_link_str, zoom_level);
        strcpy(hypo_description_str, "Station Map");
        create_map_javascript_file_link(outnameroot, -1, js_file_str, js_link_str, zoom_level);
    }
    //printf("DEBUG: js_file_str %s\n js_link_str %s\n map_link_str %s\n", js_file_str, js_link_str, map_link_str);

    // write Google javascript map html page
    if (verbose > 2)
        printf("Opening output file: %s\n", map_link_str);
    FILE * mapHtmlStream = fopen_counter(map_link_str, "w");
    if (mapHtmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", map_link_str);
        perror(tmp_str);
        return;
    }

    // start html and scripts
    fprintf(mapHtmlStream,
            "<!DOCTYPE html>\n"
            "<html>\n"
            "<!-- Automatically generated by %s [v%s %s] www.alomax.net -->\n"
            "<head>\n"
            "<meta name='viewport' content='initial - scale = 1.0, user - scalable = no'>\n"
            "<meta charset='utf-8'>\n"
            "<title>%s - %s</title>\n"
            "<style type='text/css'>\n"
            "  html { height: 100%% }\n"
            "  body { height: 100%%; margin: 0; padding: 0 }\n"
            "  #map-canvas { height: 97%% }\n"
            "h1 { font-family: sans-serif; font-size: medium; text-align: center; border: 0px solid #ffffff; background-color: #ffffff; }\n"
            "</style>\n"
            "<script src='https://maps.googleapis.com/maps/api/js?v=3.exp'></script>\n"
            "<script src='%s'></script>\n"
            "</head>\n"
            "<body>\n"
            "<div id='map-canvas'></div>\n"
            "<h1>%s &nbsp;&nbsp; [%s v%s %s] &nbsp;&nbsp; (Plate Bdys & Faults: USGS, EU-SHARE)</h1>"
            "</body>\n"
            "</html>\n",
            EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE,
            EARLY_EST_MONITOR_SHORT_NAME, hypo_description_str,
            js_file_str,
            hypo_description_str, EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE
            );

    fclose_counter(mapHtmlStream);

    // write Google javascript map html page
    if (verbose > 2)
        printf("Opening output file: %s\n", js_link_str);
    FILE * mapJavascriptStream = fopen_counter(js_link_str, "w");
    if (mapJavascriptStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", map_link_str);
        perror(tmp_str);
        return;
    }

    double lat = 0.0;
    double lon = 150.0;
    if (phypo_map != NULL) { // creating map for associated/located event
        lat = phypo_map->lat;
        lon = phypo_map->lon;
    }

    fprintf(mapJavascriptStream,
            "// Automatically generated by %s [v%s %s] www.alomax.net\n\n"
            "function initialize() {\n"
            "   var myLatlng = new google.maps.LatLng(%f,%f);\n"
            "   var mapOptions = {\n"
            "      zoom: %d,\n"
            "      center: myLatlng,\n"
            "      mapTypeId: google.maps.MapTypeId.HYBRID,\n"
            "      scaleControl: true,\n"
            "      panControlOptions: {\n"
            "         position: google.maps.ControlPosition.LEFT_BOTTOM\n"
            "      },\n"
            "      zoomControlOptions: {\n"
            "          position: google.maps.ControlPosition.LEFT_BOTTOM\n"
            "      },\n"
            "      overviewMapControl: true,\n"
            "      overviewMapControlOptions: {\n"
            "         opened: true\n"
            "      }\n"
            "   };\n"
            "var map = new google.maps.Map(document.getElementById('map-canvas'), mapOptions);\n"
            //
            "//\n"
            "var loc = window.location.href;\n"
            "var path = loc.substring(0, loc.lastIndexOf('/'));\n"
            "//\n"
            "// kml layer elements from http://earthquake.usgs.gov/learn/plate-boundaries.kmz\n"
            "var kmlUrl1 = path + '/../map_data/PlateBoundaries.kmz';\n"
            //var kmlUrl = path + '/../map_data/EarthsTectonicPlates.kmz';\n"
            "var kmlLayer1 = new google.maps.KmlLayer(kmlUrl1, {\n"
            "  suppressInfoWindows: true,\n"
            "  preserveViewport: true,\n"
            "});\n"
            "kmlLayer1.setMap(map);\n"
            "//\n"
            //
            "//\n"
            "// kml layer elements from http://diss.rm.ingv.it/share-edsf\n"
            "var kmlUrl2 = path + '/../map_data/Europe_Crustal_fault_sources_TOP.kmz';\n"
            "var kmlLayer2 = new google.maps.KmlLayer(kmlUrl2, {\n"
            "  suppressInfoWindows: false,\n"
            "  preserveViewport: true,\n"
            "});\n"
            "kmlLayer2.setMap(map);\n"
            //
            "//\n"
            "// kml layer elements from http://earthquake.usgs.gov/hazards/qfaults/KML/Historic.kmz\n"
            "var kmlUrl3 = path + '/../map_data/Historic.kmz';\n"
            "var kmlLayer3 = new google.maps.KmlLayer(kmlUrl3, {\n"
            "  suppressInfoWindows: false,\n"
            "  preserveViewport: true,\n"
            "});\n"
            "kmlLayer3.setMap(map);\n"
            //
            "//\n"
            "// kml layer elements from http://earthquake.usgs.gov/hazards/qfaults/KML/Holocene_LatestPleistocene.kmz\n"
            "// Too big! > 3MB?  var kmlUrl4 = path + '/../map_data/Holocene_LatestPleistocene.kmz';\n"
            "// kml layer elements from genarc.shp\n"
            "var kmlUrl4 = path + '/../map_data/CalifFaults.kmz';\n"
            "var kmlLayer4 = new google.maps.KmlLayer(kmlUrl4, {\n"
            "  suppressInfoWindows: false,\n"
            "  preserveViewport: true,\n"
            "});\n"
            "kmlLayer4.setMap(map);\n",
            EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE,
            lat, lon, zoom_level
            );

    if (phypo_map != NULL) { // creating map for associated/located event

        // main map hypocenter marker
        create_event_link("..", phypo_map->unique_id, event_url_str, event_link_str, feregion(phypo_map->lat, phypo_map->lon, feregion_str, FEREGION_STR_SIZE));
        fprintf(mapJavascriptStream,
                "//\n"
                "var hypoStar = {\n"
                "   path: 'M 0,-100 30,-30 100,-30 40,10 60,80 0,35 -60,80 -40,10 -100,-30 -30,-30 z',\n"
                "   fillColor: 'yellow',\n"
                "   fillOpacity: 1.0,\n"
                "   scale: 0.17,\n"
                "   strokeColor: 'red',\n"
                "   strokeWeight: 2\n"
                "};\n"
                "var hypoMarker = new google.maps.Marker({\n"
                "   zIndex: google.maps.Marker.MAX_ZINDEX + 1000,\n"
                "   position: myLatlng,\n"
                "   icon: hypoStar,\n"
                "   map: map,\n"
                "   title: \'%s\'\n"
                "});\n"
                "google.maps.event.addListener(hypoMarker, 'click', function() {\n"
                "  map.setCenter(hypoMarker.getPosition());\n"
                "  window.open('%s', '_blank');\n"
                "});\n",
                hypo_description_str, event_url_str
                );

        // add hypocenter scatter sample
        int nscatter = phypo_map->nscatter_sample;
        if (nscatter > 0) {
            fprintf(mapJavascriptStream,
                    "// Create an object containing LatLng each scatter sample.\n"
                    "var samp = {};\n"
                    );
            int ns;
            int istep = nscatter / MAX_NUM_SCATTTER_POINTS_MAP;
            if (istep < 1)
                istep = 1;
            float *sample = phypo_map->scatter_sample;
            double lat, lon;
            int nsamp = 0;
            for (ns = 0; ns < nscatter; ns += istep) {
                sample += (istep - 1) * 4; // skip samples if istep > 1
                lon = *sample;
                sample++;
                lat = *sample;
                sample++;
                sample++; // depth
                sample++; // value
                fprintf(mapJavascriptStream,
                        "samp[%d] = {\n"
                        "   center: new google.maps.LatLng(%f, %f)\n"
                        "};\n",
                        nsamp++, lat, lon
                        );
            }
            /* 20150818 AJL - no longer works in Google Maps - Circle radius is scaled to true distance
               fprintf(mapJavascriptStream,
                    "// Construct the circle for each value in samp.\n"
                    "for (var nsamp in samp) {\n"
                    "   var sampleOptions = {\n"
                    "      zIndex: 10,\n"
                    "      strokeColor: 'red',\n"
                    "      strokeOpacity: 1.0,\n"
                    "      strokeWeight: 3,\n"
                    "      fillColor: 'red',\n"
                    "      fillOpacity: 1.0,\n"
                    "      map: map,\n"
                    "      center: samp[nsamp].center,\n"
                    "      radius: 3\n"
                    "   };\n"
                    "   // Add the circle to the map.\n"
                    "   sampCircle = new google.maps.Circle(sampleOptions);\n"
                    "}\n"
                    );*/
            // 20150818 AJL - added Marker version of Circle to allow constant pixel scaling
            fprintf(mapJavascriptStream,
                    "// Construct the circle for each value in samp.\n"
                    "for (var nsamp in samp) {\n"
                    "   var sampleOptions = {\n"
                    "      map: map,\n"
                    "      zIndex: google.maps.Marker.MAX_ZINDEX + 5,\n"
                    "      position: samp[nsamp].center,\n"
                    "      clickable: false,\n"
                    "      icon: {\n"
                    "         path: google.maps.SymbolPath.CIRCLE,\n"
                    "         scale: 1, //pixels\n"
                    "         strokeColor: 'red',\n"
                    "         strokeOpacity: 1.0,\n"
                    "         strokeWeight: 0,\n"
                    "         fillColor: 'red',\n"
                    "         fillOpacity: 1.0\n"
                    "      }\n"
                    "   };\n"
                    "   // Add the circle to the map.\n"
                    "   sampCircle = new google.maps.Marker(sampleOptions);\n"
                    "}\n"
                    );
        }

        // add ellipse
        fprintf(mapJavascriptStream,
                "// Define the LatLng coordinates for the ellipse polygon's path\n"
                "var ellipseCoords = [\n"
                );
        // ellipse is centered on expectation hypo
        double expect_lat = phypo_map->expect.y;
        double expect_lon = phypo_map->expect.x;
        double xrect_hyp, yrect_hyp;
        latlon2rect_simple(expect_lat, expect_lon, &xrect_hyp, &yrect_hyp, expect_lat, expect_lon);
        double lat, lon, xrect, yrect;
        double ell_az_major = DE2RA * (90.0 - (phypo_map->ellipse.az1 + 90.0));
        double ell_len_major = phypo_map->ellipse.len2;
        double ell_len_minor = phypo_map->ellipse.len1;
        int npt = 72;
        double daz = DE2RA * 360.0 / (double) npt;
        double az = 0.0;
        int n;
        for (n = 0; n < npt; n++) {
            xrect = xrect_hyp + ell_len_major * cos(az) * cos(ell_az_major) - ell_len_minor * sin(az) * sin(ell_az_major);
            yrect = yrect_hyp + ell_len_major * cos(az) * sin(ell_az_major) + ell_len_minor * sin(az) * cos(ell_az_major);
            rect2latlon_simple(xrect, yrect, &lat, &lon, expect_lat, expect_lon);
            if (n > 0) { // IE8 (and below) will parse trailing commas in array and object literals incorrectly.
                fprintf(mapJavascriptStream, ",\n");
            }
            fprintf(mapJavascriptStream,
                    "   new google.maps.LatLng(%f, %f)",
                    lat, lon);
            az += daz;
        }
        fprintf(mapJavascriptStream,
                "];\n"
                "// Construct the polygon.\n"
                "ellipsePolygon = new google.maps.Polygon({\n"
                "   clickable: false,\n"
                "   zIndex: google.maps.Marker.MAX_ZINDEX + 10,\n"
                "   paths : ellipseCoords,\n"
                "   strokeColor : 'orange',\n"
                "   strokeOpacity : 0.7,\n"
                "   strokeWeight : 2,\n"
                "   fillColor : 'orange',\n"
                "   fillOpacity : 0.0\n"
                "});\n"
                "ellipsePolygon.setMap(map);\n"
                );

    }


    // add background hypocenter symbol markers
    int nhyp;
    for (nhyp = 0; nhyp < num_hypocenters; nhyp++) {
        HypocenterDesc* phypo = hypo_list[nhyp];
        if (phypo_map != NULL && phypo->unique_id == phypo_map->unique_id) { // do not add marker for current map hypocenter
            continue;
        }
        create_compact_hypo_info_string(phypo, hypo_description_str2);
        double prefMag = 0.0;
        if (phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE) {
            prefMag = phypo->mwpLevelStatistics.centralValue;
        } else if (phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE) {
            prefMag = phypo->mbLevelStatistics.centralValue;
        }
        double iconScale = fmax(prefMag * prefMag / 3.0, 4.0);
        int strokeWeight = 1;
        static char strokeColor[32];
        strcpy(strokeColor, "white");
        if (prefMag > 6.95) {
            strokeWeight = 3;
            strcpy(strokeColor, "red");
        }
        if (prefMag > 5.95) {
            strokeWeight = 2;
            strcpy(strokeColor, "yellow");
        }
        fprintf(mapJavascriptStream,
                "// Define the hypocenter markers\n"
                "var hypo%d = new google.maps.Marker({\n"
                "   zIndex: google.maps.Marker.MAX_ZINDEX + 20,\n"
                "   position: new google.maps.LatLng(%f,%f),\n"
                "   icon: {\n"
                "      path: google.maps.SymbolPath.CIRCLE,\n"
                "      scale: %f,\n"
                "      strokeColor: '%s',\n"
                "      strokeWeight: %d,\n"
                "   },\n"
                "   map: map,\n"
                "   title: \'%s\'\n"
                "});\n",
                nhyp, phypo->lat, phypo->lon, iconScale, strokeColor, strokeWeight, hypo_description_str2
                );
        create_event_link("..", phypo->unique_id, event_url_str, event_link_str, feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE));
        fprintf(mapJavascriptStream,
                "google.maps.event.addListener(hypo%d, 'click', function() {\n"
                "  map.setCenter(hypo%d.getPosition());\n"
                "  window.open('%s', '_blank');\n"
                "});\n",
                nhyp, nhyp, event_url_str
                );
    }

    // add station symbol markers
    fprintf(mapJavascriptStream,
            "//\n"
            "var staTriangleAssoc = {\n"
            "   path: 'M -100,100 0,-100 100,100 z',\n"
            "   fillColor: 'blue',\n"
            "   fillOpacity: 1.0,\n"
            "   scale: 0.08,\n"
            "   strokeColor: 'red',\n"
            "   strokeWeight: 2\n"
            "};\n"
            );
    fprintf(mapJavascriptStream,
            "//\n"
            "var staTriangleUnAssoc = {\n"
            "   path: 'M -100,100 0,-100 100,100 z',\n"
            "   fillColor: 'gray',\n"
            "   fillOpacity: 1.0,\n"
            "   scale: 0.08,\n"
            "   strokeColor: 'white',\n"
            "   strokeWeight: 1\n"
            "};\n"
            );
    fprintf(mapJavascriptStream,
            "//\n"
            "var staTriangleUnAssocLatencyYellow = {\n"
            "   path: 'M -100,100 0,-100 100,100 z',\n"
            "   fillColor: '#FFFFCC',\n"
            "   fillOpacity: 1.0,\n"
            "   scale: 0.08,\n"
            "   strokeColor: 'white',\n"
            "   strokeWeight: 1\n"
            "};\n"
            );
    fprintf(mapJavascriptStream,
            "//\n"
            "var staTriangleUnAssocLatencyRed = {\n"
            "   path: 'M -100,100 0,-100 100,100 z',\n"
            "   fillColor: '#FFBBBB',\n"
            "   fillOpacity: 1.0,\n"
            "   scale: 0.08,\n"
            "   strokeColor: 'white',\n"
            "   strokeWeight: 1\n"
            "};\n"
            );
    fprintf(mapJavascriptStream,
            "var info_open = null;\n"
            "google.maps.event.addListener(map, 'click', function() {\n"
            "   if (info_open !== null)\n"
            "      info_open.close();\n"
            "   info_open = null;\n"
            "});\n"
            );
    ChannelParameters* staParams;
    int duration;
    double start_time;
    int nmarker = 0;
    // loop over stations, add markers for all stations
    int nsta;
    for (nsta = 0; nsta < num_sources_total; nsta++) {
        staParams = &(channelParameters[nsta]);
        if (!staParams->process_this_channel_orientation) {
            continue;
        }
        fprintf(mapJavascriptStream,
                "// Define the station markers\n"
                "var marker%d = new google.maps.Marker({\n"
                "   zIndex: google.maps.Marker.MAX_ZINDEX + 30,\n"
                "   position: new google.maps.LatLng(%f,%f),\n"
                "   icon: %s,\n"
                "   map: map,\n"
                "   title: \'%s_%s\'\n"
                "});\n",
                nmarker, staParams->lat, staParams->lon,
                staParams->data_latency >= LATENCY_RED_CUTOFF ? "staTriangleUnAssocLatencyRed"
                : staParams->data_latency >= LATENCY_YELLOW_CUTOFF ? "staTriangleUnAssocLatencyYellow" : "staTriangleUnAssoc",
                staParams->network, staParams->station
                );
        // otime near middle of plot window
        double gcd = 0.0;
        double center_time = (time_max + time_min) / 2;
        if (phypo_map != NULL) { // creating map for associated/located event
            gcd = GCDistance(phypo_map->lat, phypo_map->lon, staParams->lat, staParams->lon);
            center_time = phypo_map->otime + get_ttime(get_P_phase_index(), gcd, phypo_map->depth); // predicted P arrival time
        }
        duration = 2 * (time_max - (int) center_time);
        duration += 60; // try to make sure latest data is displayed
        if (duration > 3600)
            duration = 3600;
        else if (duration < 300)
            duration = 300;
        start_time = center_time - (double) duration / 2.0;
        if (phypo_map != NULL) {
            fprintf(mapJavascriptStream,
                    "var info%d = new google.maps.InfoWindow({\n"
                    "   content: \'<strong>"
                    "%s&nbsp;%s<br>"
                    "dist=%.1f&deg;&nbsp; az=%.0f&deg;<br>"
                    "no pick or not associated"
                    "</strong>\'\n"
                    "});\n",
                    nmarker, create_channel_links(staParams->network, staParams->station, staParams->location, staParams->channel,
                    STREAM_NULL, "unassociated", staParams->n_int_tseries, start_time, duration,
                    tmp_str, tmp_str_2), tmp_str_2,
                    gcd, phypo_map != NULL ? GCAzimuth(phypo_map->lat, phypo_map->lon, staParams->lat, staParams->lon) : 0.0
                    );
        } else {
            fprintf(mapJavascriptStream,
                    "var info%d = new google.maps.InfoWindow({\n"
                    "   content: \'<strong>"
                    "%s&nbsp;%s<br>"
                    "latency=%s&nbsp;  qualityWt=%.1f"
                    "</strong>\'\n"
                    "});\n",
                    nmarker, create_channel_links(staParams->network, staParams->station, staParams->location, staParams->channel,
                    STREAM_NULL, "unassociated", staParams->n_int_tseries, start_time, duration,
                    tmp_str, tmp_str_2), tmp_str_2,
                    staParams->data_latency_str, staParams->qualityWeight
                    );
        }
        fprintf(mapJavascriptStream,
                "google.maps.event.addListener(marker%d, 'click', function() {\n"
                "   if (info_open !== null)\n"
                "      info_open.close();\n"
                "   info%d.open(marker%d.get('map'), marker%d);\n"
                "   info_open = info%d;\n"
                "});\n",
                nmarker, nmarker, nmarker, nmarker, nmarker
                );
        nmarker++;
    }
    if (phypo_map != NULL) { // creating map for associated/located event
        // loop over data, find associated data
        double fmquality = 0.0;
        int fmpolarity = POLARITY_UNKNOWN;
        char fmtype[32] = "Err";
        int ndata;
        //for (ndata = 0; ndata < num_de_data; ndata++) {
        for (ndata = num_de_data - 1; ndata >= 0; ndata--) { // 20150515 reverse time order so that fist association for each station is priority marker
            TimedomainProcessingData* deData = data_list[ndata];
            if (deData->is_associated == phypo_map->hyp_assoc_index + 1) {
                staParams = &(channelParameters[deData->source_id]);
                fprintf(mapJavascriptStream,
                        "// Define the station markers\n"
                        "var marker%d = new google.maps.Marker({\n"
                        "   zIndex: google.maps.Marker.MAX_ZINDEX + 30,\n"
                        "   position: new google.maps.LatLng(%f,%f),\n"
                        "   icon: staTriangleAssoc,\n"
                        "   map: map,\n"
                        "   title: \'%s_%s\'\n"
                        "});\n",
                        nmarker, staParams->lat, staParams->lon, staParams->network, staParams->station
                        );
                // pick near middle of plot window
                duration = 2 * (time_max - deData->t_time_t);
                duration += 60; // try to make sure latest data is displayed
                if (duration > 3600)
                    duration = 3600;
                start_time = (double) (deData->t_time_t) + deData->t_decsec - (double) duration / 2.0;
                fprintf(mapJavascriptStream,
                        "var info%d = new google.maps.InfoWindow({\n"
                        "   content: \'<strong>"
                        "%s&nbsp;%s<br>"
                        "dist=%.1f&deg;&nbsp; az=%.0f&deg;<br>"
                        "%s%c&nbsp; res=%.1fsec&nbsp; wt=%.2f"
                        "</strong>\'\n"
                        "});\n",
                        nmarker, create_channel_links(staParams->network, staParams->station, staParams->location, staParams->channel,
                        deData->pick_stream, pick_stream_name(deData), staParams->n_int_tseries, start_time, duration,
                        tmp_str, tmp_str_2), tmp_str_2,
                        deData->epicentral_distance, deData->epicentral_azimuth,
                        deData->phase, setPolarity(deData, &fmquality, &fmpolarity, fmtype),
                        deData->residual, deData->loc_weight
                        );

                fprintf(mapJavascriptStream,
                        "google.maps.event.addListener(marker%d, 'click', function() {\n"
                        "   if (info_open !== null)\n"
                        "      info_open.close();\n"
                        "   info%d.open(marker%d.get('map'), marker%d);\n"
                        "   info_open = info%d;\n"
                        "});\n",
                        nmarker, nmarker, nmarker, nmarker, nmarker
                        );
            }
            nmarker++;
        }
    }

    // close script and write html
    fprintf(mapJavascriptStream,
            "   }\n"
            "   google.maps.event.addDomListener(window, 'load', initialize);\n"
            );

    fclose_counter(mapJavascriptStream);

}

void printHypoDataHeaderString(char *hypoDataHeaderString) {

    sprintf(hypoDataHeaderString, "event_id assoc_ndx loc_seq_num ph_assoc ph_used dmin(deg) gap1(deg) gap2(deg) atten sigma_otime(sec) otime(UTC) lat(deg) lon(deg) errH(km) depth(km)"
            " errZ(km) Q T50Ex n Td(sec) n TdT50Ex WL_col mb n Mwp n T0(sec) n Mwpd n region n_sta_tot n_sta_active assoc_latency");

}

/** read hypo data string in csv format
 */
int readHypoDataString(HypocenterDesc* phypo, char *hypoDataString, long *pfirst_assoc_latency) {

    char time_str[64];


#ifdef USE_MWP_MO_POS_NEG
    static char format[] = "%ld %d %d %d %d %lf %lf %lf %lf %lf %s %lf %lf %lf %lf %lf %s %lf %d %lf %d %lf %s %lf %d %lf %d %lf %d %lf %d %*s %d %d %ld %lf %d";
#else
    static char format[] = "%ld %d %d %d %d %lf %lf %lf %lf %lf %s %lf %lf %lf %lf %lf %s %lf %d %lf %d %lf %s %lf %d %lf %d %lf %d %lf %d %*s %d %d %ld";
#endif

    int istat = sscanf(hypoDataString,
            format,
            &phypo->unique_id, &phypo->hyp_assoc_index, &phypo->loc_seq_num, &phypo->nassoc, &phypo->nassoc_P,
            &phypo->dist_min, &phypo->gap_primary, &phypo->gap_secondary,
            &phypo->linRegressPower.power, // amplitude attenuation
            &phypo->ot_std_dev, time_str,
            &phypo->lat, &phypo->lon, &phypo->errh, &phypo->depth, &phypo->errz, phypo->qualityIndicators.quality_code,
            &phypo->t50ExLevelStatistics.centralValue, &phypo->t50ExLevelStatistics.numLevel,
            &phypo->taucLevelStatistics.centralValue, &phypo->taucLevelStatistics.numLevel,
            &phypo->tdT50ExLevelStatistics.centralValue, phypo->warningLevelString,
            &phypo->mbLevelStatistics.centralValue, &phypo->mbLevelStatistics.numLevel,
            &phypo->mwpLevelStatistics.centralValue, &phypo->mwpLevelStatistics.numLevel,
            &phypo->t0LevelStatistics.centralValue, &phypo->t0LevelStatistics.numLevel,
            &phypo->mwpdLevelStatistics.centralValue, &phypo->mwpdLevelStatistics.numLevel,
            // region ignored
            &phypo->nstaHasBeenActive, &phypo->nstaIsActive, pfirst_assoc_latency
#ifdef USE_MWP_MO_POS_NEG
            , &phypo->mwpdMoPosNegLevelStatistics.centralValue, &phypo->mwpdMoPosNegLevelStatistics.numLevel
#endif
            );

    phypo->otime = string2timeDecSec(time_str);

    return (istat);

}

/** print hypo data string in csv format
 */
void printHypoDataString(HypocenterDesc* phypo, char *hypoDataString, int iDecPrecision) {

    feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE);
    char *c = feregion_str;
    while ((c = strchr(c, ' ')) != NULL)
        *c = '_';

    long first_assoc_latency = phypo->first_assoc_time - (long) (0.5 + phypo->otime);
    first_assoc_latency = first_assoc_latency < 0 ? -1 : first_assoc_latency;

    // hypo data csv string
    if (iDecPrecision) {
        sprintf(hypoDataString, "%ld %d %d %d %d %.3f %.3f %.3f ",
                phypo->unique_id, phypo->hyp_assoc_index + 1, phypo->loc_seq_num, phypo->nassoc, phypo->nassoc_P, phypo->dist_min, phypo->gap_primary, phypo->gap_secondary
                );
        // amplitude attenuation
        sprintf(hypoDataString + strlen(hypoDataString), "%.3f ",
                phypo->linRegressPower.power
                );
        sprintf(hypoDataString + strlen(hypoDataString), "%.3f %s %.3f %.3f %.2f %.2f %.2f %s %.3f %d %.3f %d %.3f %s %.3f %d %.3f %d %.2f %d %.3f %d %s",
                phypo->ot_std_dev,
                timeDecSec2string(phypo->otime, tmp_str, DEFAULT_TIME_FORMAT),
                phypo->lat, phypo->lon, phypo->errh, phypo->depth, phypo->errz, phypo->qualityIndicators.quality_code,
                phypo->t50ExLevelStatistics.centralValue, phypo->t50ExLevelStatistics.numLevel,
                phypo->taucLevelStatistics.centralValue, phypo->taucLevelStatistics.numLevel,
                phypo->tdT50ExLevelStatistics.centralValue, phypo->warningLevelString,
                phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mbLevelStatistics.centralValue : MB_INVALID, phypo->mbLevelStatistics.numLevel,
                phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpLevelStatistics.centralValue : MWP_INVALID, phypo->mwpLevelStatistics.numLevel,
                phypo->t0LevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->t0LevelStatistics.centralValue : T0_INVALID, phypo->t0LevelStatistics.numLevel,
                phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpdLevelStatistics.centralValue : MWPD_INVALID, phypo->mwpdLevelStatistics.numLevel,
                feregion_str
                );
        // available station count
        // 20141212 AJL - added to enable later statistical analysis of events (e.g. mag vs. n_sta/n_sta_available)
        sprintf(hypoDataString + strlen(hypoDataString), " %d %d %ld",
                phypo->nstaHasBeenActive, phypo->nstaIsActive, first_assoc_latency
                );
#ifdef USE_MWP_MO_POS_NEG
        sprintf(hypoDataString + strlen(hypoDataString), " %.3f %d",
                phypo->mwpdMoPosNegLevelStatistics.centralValue, phypo->mwpdMoPosNegLevelStatistics.numLevel
                );
#endif
    } else {

        sprintf(hypoDataString, "%ld %d %d %d %d %.1f %.1f %.1f ",
                phypo->unique_id, phypo->hyp_assoc_index + 1, phypo->loc_seq_num, phypo->nassoc, phypo->nassoc_P, phypo->dist_min, phypo->gap_primary, phypo->gap_secondary
                );
        // amplitude attenuation
        sprintf(hypoDataString + strlen(hypoDataString), "%.1f ",
                phypo->linRegressPower.power
                );
        sprintf(hypoDataString + strlen(hypoDataString), "%.1f %s %.1f %.1f %.0f %.0f %.0f %s %.1f %d %.1f %d %.1f %s %.1f %d %.1f %d %.0f %d %.1f %d %s",
                phypo->ot_std_dev,
                time2string(phypo->otime, tmp_str),
                phypo->lat, phypo->lon, phypo->errh, phypo->depth, phypo->errz, phypo->qualityIndicators.quality_code,
                phypo->t50ExLevelStatistics.centralValue, phypo->t50ExLevelStatistics.numLevel,
                phypo->taucLevelStatistics.centralValue, phypo->taucLevelStatistics.numLevel,
                phypo->tdT50ExLevelStatistics.centralValue, phypo->warningLevelString,
                phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mbLevelStatistics.centralValue : MB_INVALID, phypo->mbLevelStatistics.numLevel,
                phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpLevelStatistics.centralValue : MWP_INVALID, phypo->mwpLevelStatistics.numLevel,
                phypo->t0LevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->t0LevelStatistics.centralValue : T0_INVALID, phypo->t0LevelStatistics.numLevel,
                phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpdLevelStatistics.centralValue : MWPD_INVALID, phypo->mwpdLevelStatistics.numLevel,
                feregion_str
                );
        // available station counts, elapsed time after origin of first association/location
        // 20141212 AJL - added to enable later statistical analysis of events (e.g. mag vs. n_sta/n_sta_available)
        sprintf(hypoDataString + strlen(hypoDataString), " %d %d %ld",
                phypo->nstaHasBeenActive, phypo->nstaIsActive, first_assoc_latency
                );
#ifdef USE_MWP_MO_POS_NEG

        sprintf(hypoDataString + strlen(hypoDataString), " %.1f %d",
                phypo->mwpdMoPosNegLevelStatistics.centralValue, phypo->mwpdMoPosNegLevelStatistics.numLevel
                );
#endif
    }
}

void printHypoMessageHtmlHeaderString(char *messageString) {

    sprintf(messageString, "<tr align=right bgcolor=\"#BBBBBB\">\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;event id</th><th>locSeq</th><th>ph&nbsp;&nbsp;<br>assoc</th><th>ph&nbsp;&nbsp;<br>used</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;dmin<br>&nbsp;(deg)</th></th><th>&nbsp;gap1<br>&nbsp;(deg)</th></th><th>&nbsp;gap2<br>&nbsp;(deg)</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;atten</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;&sigma; otime<br>&nbsp;(sec)</th><th>&nbsp;otime<br>&nbsp;(UTC)</th><th>&nbsp;lat&nbsp;<br>&nbsp;(deg)</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;lon&nbsp;<br>&nbsp;(deg)</th><th align=left>&nbsp;&nbsp;errH<br>&nbsp;(km)</th>\n");
    sprintf(messageString + strlen(messageString), "<th>depth<br>&nbsp;(km)</th><th align=left>&nbsp;&nbsp;errZ<br>&nbsp;(km)&nbsp;</th>\n");
    sprintf(messageString + strlen(messageString), "<th>Q</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;T50Ex</th><th align=left>&nbsp;n</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;Td&nbsp;<br>&nbsp;(sec)</th><th align=left>&nbsp;n</th><th>&nbsp;TdT50Ex</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;mb</th><th align=left>&nbsp;n</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;Mwp</th><th align=left>&nbsp;n</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;T0&nbsp;<br>&nbsp;(sec)</th><th align=left>&nbsp;n</th>\n");
    sprintf(messageString + strlen(messageString), "<th>&nbsp;Mwpd</th><th align=left>&nbsp;n</th>\n");
    sprintf(messageString + strlen(messageString), "<th align=left>fm</th>\n");
    sprintf(messageString + strlen(messageString), "<th align=left>region</th>\n");
    sprintf(messageString + strlen(messageString), "</tr>\n");

}

/** print hypo data string in html format
 */

void printHypoMessageHtmlString(HypocenterDesc* phypo, char *messageString, char* hypoBackground, char* rowBackground, long event_index, long event_id) {

    // set feregion string
    feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE);

    // hypo message html
    sprintf(messageString, "<tr align=right %s>", hypoBackground);
    create_event_link(".", event_id, event_url_str, event_link_str, feregion_str);
    sprintf(messageString + strlen(messageString), "<td>%s&nbsp;%ld&nbsp;</a></td>", event_link_str, event_index);
    // location sequence number
    sprintf(messageString + strlen(messageString), "<td>%d</td>", phypo->loc_seq_num);
    sprintf(messageString + strlen(messageString), "<td>%d</td><td>%d</td><td>%.1f</td></td><td>%.0f</td></td><td>%.0f</td><td>%.2f</td><td>%.1f&nbsp;&nbsp;</td>",
            phypo->nassoc, phypo->nassoc_P, phypo->dist_min, phypo->gap_primary, phypo->gap_secondary, phypo->linRegressPower.power, phypo->ot_std_dev);
    sprintf(messageString + strlen(messageString), "<td><strong>%s</strong></td><td><strong>%.2f</strong></td><td><strong>%.2f</strong></td>",
            time2string(phypo->otime, tmp_str), phypo->lat, phypo->lon);
    // errh, with map link
    create_map_link("./", phypo->unique_id, map_link_str, MAP_LINK_BIG_ZOOM);
    sprintf(messageString + strlen(messageString), "<td align=left>&nbsp;&nbsp;&nbsp;<a href=\"%s\" title=\"View epicenter error in Google Maps\" target=map_%ld>%.0f</a></td>",
            map_link_str, phypo->unique_id, phypo->errh);
    //
    sprintf(messageString + strlen(messageString), "<td><strong>%.0f</strong></td><td align=left>&nbsp;&nbsp;&nbsp;%.0f</td>",
            phypo->depth, phypo->errz);
    sprintf(messageString + strlen(messageString), "<td align=\"center\" bgcolor=%s><strong>%s</strong></td>",
            getQualityCodeBackgroundColorStringHtml(phypo->qualityIndicators.quality_code), phypo->qualityIndicators.quality_code);
    int highlight = 0;
    highlight = setLevelStringHtml(phypo->t50ExLevelStatistics.numLevel, phypo->t50ExLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, T50EX_RED_CUTOFF, T50EX_YELLOW_CUTOFF, -1.0, warning_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->t50ExLevelStatistics.centralValue, highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->t50ExLevelStatistics.numLevel);
    highlight = setLevelStringHtml(phypo->taucLevelStatistics.numLevel, phypo->taucLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, TAUC_RED_CUTOFF, TAUC_YELLOW_CUTOFF, -1.0, warning_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->taucLevelStatistics.centralValue, highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->taucLevelStatistics.numLevel);
    highlight = setLevelStringHtml(IMIN(phypo->t50ExLevelStatistics.numLevel, phypo->taucLevelStatistics.numLevel), phypo->tdT50ExLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, TDT50EX_RED_CUTOFF, TDT50EX_YELLOW_CUTOFF, -1.0, warning_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->tdT50ExLevelStatistics.centralValue, highlight ? "&nbsp;</strong>" : "&nbsp;");
    highlight = setLevelStringHtml(phypo->mbLevelStatistics.numLevel, phypo->mbLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, MB_RED_CUTOFF, MB_YELLOW_CUTOFF, MB_INVALID, magnitude_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mbLevelStatistics.centralValue : MB_INVALID,
            highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->mbLevelStatistics.numLevel);
    highlight = setLevelStringHtml(phypo->mwpLevelStatistics.numLevel, phypo->mwpLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, MWP_RED_CUTOFF, MWP_YELLOW_CUTOFF, MWP_INVALID, magnitude_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->mwpLevelStatistics.centralValue : MWP_INVALID,
            highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->mwpLevelStatistics.numLevel);
    highlight = setLevelStringHtml(phypo->t0LevelStatistics.numLevel, phypo->t0LevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, T0_RED_CUTOFF, T0_YELLOW_CUTOFF, T0_INVALID, warning_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.0f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->t0LevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE ? phypo->t0LevelStatistics.centralValue : T0_INVALID,
            highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->t0LevelStatistics.numLevel);
    int nLevels = phypo->mwpdLevelStatistics.numLevel;
    // check if Mwpd below minimum value or Mwp too low

    if (phypo->mwpdLevelStatistics.centralValue < MWPD_MIN_VALUE_USE || phypo->mwpLevelStatistics.centralValue < MIN_MWP_FOR_VALID_MWPD)
        nLevels = -1; // force grey color
    highlight = setLevelStringHtml(nLevels, phypo->mwpdLevelStatistics.centralValue, "bgcolor=", rowBackground, levelStringHtml,
            MIN_NUMBER_VALUES_USE, MWPD_RED_CUTOFF, MWPD_YELLOW_CUTOFF, MWPD_INVALID, magnitude_colors_show);
    sprintf(messageString + strlen(messageString), "<td %s>%s%.1f%s</td><td %s align=left>&nbsp;%d</td>",
            levelStringHtml, highlight ? "<strong>" : "", phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE
            ? phypo->mwpdLevelStatistics.centralValue : MWPD_INVALID, highlight ? "&nbsp;</strong>" : "&nbsp;", rowBackground, phypo->mwpdLevelStatistics.numLevel);
    // focal mechanism image
    create_mech_image_tag(".", event_id, image_tag_str, feregion_str);
    sprintf(messageString + strlen(messageString), "<td align=center>%s</td>", image_tag_str);
    // region, with map link
    create_map_link("./", phypo->unique_id, map_link_str, MAP_LINK_MED_ZOOM);
    sprintf(messageString + strlen(messageString), "<td align=left>&nbsp;<a href=\"%s\" title=\"View epicenter in Google Maps\" target=map_%ld>%s</a></td>",
            map_link_str, phypo->unique_id, feregion_str);
    //
    sprintf(messageString + strlen(messageString), "</tr>");

}

/** write a hypocenter kml file */

int write_hypocenter_kml(char *outnameroot, HypocenterDesc* phypo, int verbose) {

    char outname[1024];
    sprintf(outname, "%s/hypo.%d.kml", outnameroot, phypo->hyp_assoc_index);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * kmlStream = fopen_counter(outname, "w");
    if (kmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }

    fprintf(kmlStream,
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            "<kml xmlns=\"http://earth.google.com/kml/2.1\">"
            "	<Document>"
            "		<name>Sismicite</name>"
            "		<description><![CDATA[<B>SISMOS à l'École - Sismicité</B>]]></description>"
            "		<LookAt id=\"khLookAt704\">"
            "			<longitude>7.0</longitude>"
            "			<latitude>44.0</latitude>"
            "			<altitude>0</altitude>"
            "			<range>250000</range>"
            "			<tilt>0</tilt>"
            "			<heading>0</heading>"
            "		</LookAt>"
            "		<Style id=\"myDefaultStyles\">"
            "			<BalloonStyle>"
            "				<!-- a background color for the balloon -->"
            "				<bgColor>C5D8F5</bgColor>"
            "				<!-- styling of the balloon text -->"
            "				<text><![CDATA["
            "					       <strong><font color=\"#CC0000\" size=\"+3\">$[name]</font></strong>"
            "					       <br/><br/>"
            "					       <font color=\"#000033\" size=\"+2\">$[description]</font>"
            "					       <br/><br/>"
            "					       \"SISMOS à l'École\""
            "					       <br/><br/>"
            "					       <a href=\"http://aster.unice.fr/EduSeisExplorer\">http://aster.unice.fr/EduSeisExplorer</a>"
            "					       ]]></text>"
            "			       </BalloonStyle>"
            "			<IconStyle>"
            "				<color>ff00ffff</color>"
            "				<scale>0.3</scale>"
            "				<Icon>"
            "					<href>/windows/E/www/free/projects/google/cataseism/images/dot.gif</href>"
            "				</Icon>"
            "				<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"
            "			</IconStyle>"
            "			<LabelStyle>"
            "				<color>ccccccff</color>"
            "				<scale>0.75</scale>"
            "			</LabelStyle>"
            "			<LineStyle>"
            "				<color>ffffffff</color>"
            "				<width>15</width>"
            "			</LineStyle>"
            "			<PolyStyle>"
            "				<color>ffffff</color>"
            "				<colorMode>random</colorMode>"
            "			</PolyStyle>"
            "		</Style>"
            "  <Placemark><styleUrl>#myDefaultStyles</styleUrl><Point><coordinates>007.3750,43.8010,-002500</coordinates></Point>  <description><![CDATA[<B> 2000  1019 063254.1  7,3750 43,8010  2,5   M= 1,1]]></description> </Placemark>"
            "	</Document>"
            "</kml>"
            );

    fclose_counter(kmlStream);

    return (0);
}


/** Send a hypocenter mail message */
static char send_mail_params_copy[STANDARD_STRLEN] = "\0";
static char mail_file[STANDARD_STRLEN] = "\0";
static char info_link[STANDARD_STRLEN] = "\0";
static char map_link[STANDARD_STRLEN] = "\0";
static char event_url[STANDARD_STRLEN] = "\0";
static char event_link[STANDARD_STRLEN] = "\0";
static char mail_from[STANDARD_STRLEN] = "\0";
static char mail_to[16384] = "\0";
static char mail_trigger_string[STANDARD_STRLEN] = "\0";
static char mb_str[16] = "\0";
static char mwp_str[16] = "\0";
static char t0_str[16] = "\0";
static char mwpd_str[16] = "\0";
static char alarm_str[16] = "\0";
static char warningLevelString[16] = "\0";
static char alarm_color_open[1024] = "\0";
static char alarm_color_close[1024] = "\0";
static char hypo_info_string[STANDARD_STRLEN] = "\0";
static char latency_string[STANDARD_STRLEN] = "\0";
static char tsunami_subject[STANDARD_STRLEN] = "\0";
static char tsunami_info_string[STANDARD_STRLEN] = "\0";
static char mail_subject[STANDARD_STRLEN] = "\0";
static char sys_command[16384] = "\0";
static message_trigger_theshold messageTriggerThresholdInitial;
//
static char alert_file[STANDARD_STRLEN] = "\0";

int send_hypocenter_alert(HypocenterDesc* phypo, int is_new_hypocenter, HypocenterDesc* pexisting_hypo, char *sendMailParams, char *outnameroot, char *agencyId, int verbose) {

    // time processing
    time_t current_time = time(&current_time);

    // check if new or existing hypocenter and if passed a magnitude cutoff for issuing a mail
    int send_alert = 0;
    int trigger_time = 0;
    int trigger_mb = 0;
    int trigger_mwp = 0;
    int trigger_mwpd = 0;
    int trigger_alarm = 0;
    int trigger_alarm_dropped = 0;
    if (is_new_hypocenter || pexisting_hypo == NULL) {
        phypo->messageTriggerThreshold = messageTriggerThresholdInitial;
    } else {
        //printf("DEBUG: pexisting_hypo->messageTriggerThreshold %g %g %g %g\n",
        //        pexisting_hypo->messageTriggerThreshold.mb, pexisting_hypo->messageTriggerThreshold.mwp, pexisting_hypo->messageTriggerThreshold.mwpd, pexisting_hypo->messageTriggerThreshold.alarm);
        phypo->messageTriggerThreshold = pexisting_hypo->messageTriggerThreshold;
    }
    //printf("DEBUG: phypo->messageTriggerThreshold %g %g %g %g\n", phypo->messageTriggerThreshold.mb, phypo->messageTriggerThreshold.mwp, phypo->messageTriggerThreshold.mwpd, phypo->messageTriggerThreshold.alarm);


    // 20160304 AJL - only send alert for acceptable (A or B) quality events
    if (strstr(LOC_QUALITY_ACCEPTABLE, phypo->qualityIndicators.quality_code) == NULL) {
        return (0);
    }


    // new, check if Mwp, mb, Mwpd or TdT50Ex greater than cutoff
    if (phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mbLevelStatistics.centralValue >= phypo->messageTriggerThreshold.mb) {
        trigger_mb = 1;
        send_alert = 1;
        phypo->messageTriggerThreshold.mb = phypo->mbLevelStatistics.centralValue + MB_MIN_MAIL_INCREMENT;
    }
    if (phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mwpLevelStatistics.centralValue >= phypo->messageTriggerThreshold.mwp) {
        trigger_mwp = 1;
        send_alert = 1;
        phypo->messageTriggerThreshold.mwp = phypo->mwpLevelStatistics.centralValue + MWP_MIN_MAIL_INCREMENT;
    }
    if ((phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mwpLevelStatistics.centralValue >= MIN_MWP_FOR_VALID_MWPD) && // Mwpd requires Mwp above 7
            (phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mwpdLevelStatistics.centralValue >= phypo->messageTriggerThreshold.mwpd)) {
        trigger_mwpd = 1;
        send_alert = 1;
        //phypo->messageTriggerThreshold.mwpd = phypo->mwpLevelStatistics.centralValue + MWPD_MIN_MAIL_INCREMENT;   // 20110623 AJL - Bug fix: mwp should be mwpd
        phypo->messageTriggerThreshold.mwpd = phypo->mwpdLevelStatistics.centralValue + MWPD_MIN_MAIL_INCREMENT;
    }
    if (phypo->tdT50ExLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE) {
        if (phypo->tdT50ExLevelStatistics.centralValue >= phypo->messageTriggerThreshold.alarm) { // value rises above trigger threshold
            trigger_alarm = 1;
            send_alert = 1;
            phypo->messageTriggerThreshold.alarm = phypo->tdT50ExLevelStatistics.centralValue + TDT50EX_MIN_MAIL_INCREMENT;
        } else {
            double previous_trigger_level = phypo->messageTriggerThreshold.alarm - TDT50EX_MIN_MAIL_INCREMENT;
            if (previous_trigger_level > TDT50EX_RED_CUTOFF && phypo->tdT50ExLevelStatistics.centralValue < TDT50EX_YELLOW_CUTOFF) { // alarm level dropped from above RED to below YELLOW
                trigger_alarm_dropped = 1;
                send_alert = 1;
                phypo->messageTriggerThreshold.alarm = phypo->tdT50ExLevelStatistics.centralValue + TDT50EX_MIN_MAIL_INCREMENT;
            }
        }
    }

    // if an alert has been sent, resend after specified time delay after hypo origin time
    double time_since_otime = difftime(current_time, (time_t) phypo->otime);
    //printf("DEBUG: ALERT: if (phypo->alert_sent_count [%d] > 0 && !phypo->alert_sent_time [%d] && time_since_otime [%f] > ALERT_RESEND_TIME_DELAY [%d] (diff=%f)\n",
    //        phypo->alert_sent_count, phypo->alert_sent_time, time_since_otime, phypo->messageTriggerThreshold.resend_time_delay, time_since_otime - phypo->messageTriggerThreshold.resend_time_delay);
    if (phypo->messageTriggerThreshold.resend_time_delay > 0
            && phypo->alert_sent_count > 0 && !phypo->alert_sent_time
            && time_since_otime > phypo->messageTriggerThreshold.resend_time_delay) {
        phypo->alert_sent_time = 1; // flags not to send again alert on time delay
        trigger_time = 1;
        send_alert = 1;
    }

    if (!send_alert)
        return (0);

    phypo->alert_sent_count++;

    // mails trigger string
    sprintf(mail_trigger_string, " ");
    if (trigger_time == 1)
        sprintf(mail_trigger_string, " Time_since_origin>%ds ", phypo->messageTriggerThreshold.resend_time_delay);
    if (trigger_mb == 1)
        strcat(mail_trigger_string, "mb ");
    if (trigger_mwp == 1)
        strcat(mail_trigger_string, "Mwp ");
    if (trigger_mwpd == 1)
        strcat(mail_trigger_string, "Mwpd ");
    if (trigger_alarm == 1)
        strcat(mail_trigger_string, "Td*T50Ex ");
    if (trigger_alarm_dropped == 1)
        strcat(mail_trigger_string, "Td*T50Ex_Decreased ");

    // set magnitudes
    if (phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE)
        sprintf(mb_str, "%.1f", phypo->mbLevelStatistics.centralValue);
    else strcpy(mb_str, "N/A");
    if (phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE)
        sprintf(mwp_str, "%.1f", phypo->mwpLevelStatistics.centralValue);
    else strcpy(mwp_str, "N/A");
    if (phypo->t0LevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE)
        sprintf(t0_str, "%.0f", phypo->t0LevelStatistics.centralValue);
    else strcpy(t0_str, "N/A");
    if (phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE && phypo->mwpLevelStatistics.centralValue >= MIN_MWP_FOR_VALID_MWPD)
        sprintf(mwpd_str, "%.1f", phypo->mwpdLevelStatistics.centralValue);
    else strcpy(mwpd_str, "N/A");
    if (phypo->tdT50ExLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE)
        sprintf(alarm_str, "%.1f", phypo->tdT50ExLevelStatistics.centralValue);
    else strcpy(alarm_str, "N/A");

    // check Td*T50Ex level
    strcpy(warningLevelString, phypo->warningLevelString);
    strcpy(alarm_color_close, "</font>");
    if (!warning_colors_show) {
        warningLevelString[0] = (char) 0;
        strcpy(alarm_color_open, "");
        strcpy(alarm_color_close, "");
    }
    if ((phypo->t50ExLevelStatistics.numLevel < MIN_NUMBER_VALUES_USE) || (phypo->taucLevelStatistics.numLevel < MIN_NUMBER_VALUES_USE)) {
        if (warning_colors_show)
            strcpy(alarm_color_open, "<font style=\"background-color: #BDBDBD\">");
        sprintf(tsunami_info_string, "Tsunami evaluation:\n<br>&nbsp;&nbsp;&nbsp;not available (Td*T50Ex Level: %s%s%s%s)",
                alarm_str, alarm_color_open, warningLevelString, alarm_color_close);
        tsunami_subject[0] = (char) 0;
    } else if (phypo->tdT50ExLevelStatistics.centralValue >= TDT50EX_RED_CUTOFF) {
        if (warning_colors_show)
            strcpy(alarm_color_open, "<font style=\"background-color: #FF5B5B\">");
        sprintf(tsunami_info_string, "Tsunami evaluation:\n<br>&nbsp;&nbsp;&nbsp;<strong>LIKELY TSUNAMIGENIC EVENT if shallow, non-strike-slip, oceanic event (Td*T50Ex Level: %s%s%s%s)</strong>",
                alarm_str, alarm_color_open, warningLevelString, alarm_color_close);
        sprintf(tsunami_subject, "LIKELY TSUNAMIGENIC: ");
    } else if (phypo->tdT50ExLevelStatistics.centralValue >= TDT50EX_YELLOW_CUTOFF) {
        if (warning_colors_show)
            strcpy(alarm_color_open, "<font style=\"background-color: #FFFF5B\">");
        sprintf(tsunami_info_string, "Tsunami evaluation:\n<br>&nbsp;&nbsp;&nbsp;<strong>Possible tsunamigenic event if shallow, non-strike-slip, oceanic event (Td*T50Ex Level: %s%s%s%s)</strong>",
                alarm_str, alarm_color_open, warningLevelString, alarm_color_close);
        sprintf(tsunami_subject, "Possible tsunamigenic: ");
    } else {
        if (warning_colors_show)
            strcpy(alarm_color_open, "<font style=\"background-color: #5BFF5B\">");
        sprintf(tsunami_info_string, "Tsunami evaluation:\n<br>&nbsp;&nbsp;&nbsp;unlikely tsunamigenic event (Td*T50Ex Level: %s%s%s%s)",
                alarm_str, alarm_color_open, warningLevelString, alarm_color_close);
        tsunami_subject[0] = (char) 0;
    }
    if (!tsunami_evaluation_write) {
        tsunami_info_string[0] = (char) 0;
        tsunami_subject[0] = (char) 0;
    }
    if (tsunami_evaluation_write) {
        strcat(tsunami_info_string,
                "\n<br>&nbsp;&nbsp;&nbsp;(DISCLAIMER: Tsunamigenic evaluation is based on the value of the Td*T50Ex Level."
                "\n<br>&nbsp;&nbsp;&nbsp;&nbsp;This evaluation only concerns the likelihood that this earthquake generated a regional or larger scale tsunami."
                "\n<br>&nbsp;&nbsp;&nbsp;&nbsp;This evaluation DOES NOT concern the size and effects of a possible tsunami, which depend on details of the"
                "\n<br>&nbsp;&nbsp;&nbsp;&nbsp;earthquake source, ocean bathymetry, coastal distances and population density, and many other factors.)"
                "\n<br>&nbsp;&nbsp;&nbsp;&nbsp;This evaluation DOES NOT apply to and DOES NOT EXCLUDE the possibility that this earthquake generated a local tsunami,"
                "\n<br>&nbsp;&nbsp;&nbsp;&nbsp;which can be destructive along coasts within a 100 km or more from the earthquake epicenter."
                );
    }

    // set FE-region
    feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE);

    // print hypocenter information string
    if (warning_colors_show) {
        sprintf(hypo_info_string, "%s %s mb%s Mwp%s T0%s Mwpd%s TdT50Ex%s%s%s",
                feregion_str,
                time2string(phypo->otime, tmp_str),
                mb_str, mwp_str, t0_str, mwpd_str,
                alarm_str, strlen(warningLevelString) > 0 ? "|" : "", warningLevelString
                );
    } else {
        sprintf(hypo_info_string, "%s %s mb%s Mwp%s T0%s Mwpd%s TdT50Ex%s",
                feregion_str,
                time2string(phypo->otime, tmp_str),
                mb_str, mwp_str, t0_str, mwpd_str,
                alarm_str
                );
    }
    sprintf(mail_subject, "[Early-est] %s%s", tsunami_subject, hypo_info_string);

    // print multi-line hypocenter information string
    double azMaxHorUnc = phypo->ellipse.az1 + 90.0;
    if (azMaxHorUnc >= 360.0)
        azMaxHorUnc -= 360.0;
    if (azMaxHorUnc >= 180.0)
        azMaxHorUnc -= 180.0;
    sprintf(hypo_info_string,
            "&nbsp;&nbsp;&nbsp;<strong>mb: %s</strong> [%d] &nbsp;&nbsp;&nbsp;<strong>Mwp: %s</strong> [%d] &nbsp;&nbsp;&nbsp;<strong>"
            "T0: %ss</strong> [%d] &nbsp;&nbsp;&nbsp;<strong>Mwpd: %s</strong> [%d] &nbsp;&nbsp;&nbsp;",
            mb_str, phypo->mbLevelStatistics.numLevel, mwp_str, phypo->mwpLevelStatistics.numLevel, t0_str, phypo->t0LevelStatistics.numLevel, mwpd_str, phypo->mwpdLevelStatistics.numLevel
            );
    if (warning_colors_show) {
        sprintf(hypo_info_string + strlen(hypo_info_string),
                "<strong>Td*T50Ex Level: %s</strong> [%d]%s%s%s\n",
                alarm_str, phypo->tdT50ExLevelStatistics.numLevel, alarm_color_open, warningLevelString, alarm_color_close
                );
    } else {
        sprintf(hypo_info_string + strlen(hypo_info_string),
                "<strong>Td*T50Ex Level: %s</strong> [%d]\n",
                alarm_str, phypo->tdT50ExLevelStatistics.numLevel
                );
    }
    sprintf(hypo_info_string + strlen(hypo_info_string),
            "<br>"
            "&nbsp;&nbsp;&nbsp;<strong>%s</strong>\n"
            "<br>"
            "&nbsp;&nbsp;&nbsp;<strong>%s UTC</strong>\n"
            "<br>"
            "&nbsp;&nbsp;&nbsp;<strong>lat: %.1fdeg &nbsp;&nbsp;lon: %.1fdeg</strong> &nbsp;&nbsp;[+/-%.0fkm]  &nbsp;&nbsp;&nbsp;<strong>depth: %.0fkm</strong>[+/-%.0fkm]\n"
            "  &nbsp;&nbsp;&nbsp;<strong>Q: %s</strong>\n"
            "<br>"
            "&nbsp;&nbsp;&nbsp;stdErr: %.1fs &nbsp;&nbsp;assocPhases: %d &nbsp;&nbsp;usedPhases: %d"
            " &nbsp;&nbsp;minDist: %.1fdeg &nbsp;&nbsp;maxDist: %.1fdeg &nbsp;&nbsp;AzGap: %.0fdeg &nbsp;&nbsp;2ndAzGap: %.0fdeg\n",
            feregion_str,
            timeDecSec2string(phypo->otime, tmp_str, DEFAULT_TIME_FORMAT),
            phypo->lat, phypo->lon, phypo->errh, phypo->depth, phypo->errz, phypo->qualityIndicators.quality_code,
            phypo->ot_std_dev, phypo->nassoc, phypo->nassoc_P, phypo->dist_min, phypo->dist_max, phypo->gap_primary, phypo->gap_secondary
            );
    sprintf(hypo_info_string + strlen(hypo_info_string),
            " &nbsp;&nbsp;AmpAtten: %.2f",
            phypo->linRegressPower.power
            );
    sprintf(hypo_info_string + strlen(hypo_info_string),
            "<br>"
            "&nbsp;&nbsp;&nbsp;68%%ConfEllipse: &nbsp;&nbsp;minHorUnc: %.1fkm &nbsp;&nbsp;maxHorUnc: %.1fkm &nbsp;&nbsp;azMaxHorUnc: %.0fdeg",
            phypo->ellipse.len1, phypo->ellipse.len2, azMaxHorUnc
            );

    // print latency string
    timeDiff2string(current_time, (time_t) phypo->otime, latency_string);

    sprintf(mail_file, "%s/mail_%ld.txt", outnameroot, phypo->unique_id);
    FILE * mailTmpStream = fopen_counter(mail_file, "w");
    if (mailTmpStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", mail_file);
        perror(tmp_str);
        return (0);
    }

    // parse and check sendmail parameters
    int can_send_mail = 0;
    strcpy(info_link, "ERROR");
    strcpy(mail_from, "ERROR");
    strcpy(mail_to, "ERROR");
    strncpy(send_mail_params_copy, sendMailParams, STANDARD_STRLEN);
    while (1) {
        // info link
        char *str_pos = strtok(send_mail_params_copy, ",");
        if (str_pos == NULL) break;
        strcpy(info_link, str_pos);
        // mail from
        str_pos = strtok(NULL, ",");
        if (str_pos == NULL) break;
        strcpy(mail_from, str_pos);
        // mail to
        str_pos = strtok(NULL, ",");
        if (str_pos == NULL) break;
        strcpy(mail_to, str_pos);
        //
        can_send_mail = 1;
        // additional mail to
        while (1) {
            str_pos = strtok(NULL, ",");
            if (str_pos == NULL) break;
            strcat(mail_to, ",");
            strcat(mail_to, str_pos);
        }
        //printf("===========>DEBUG: mail_to: %s\n", mail_to);
        break;
    }

    // create Google map link
    create_map_link(info_link, phypo->unique_id, map_link, MAP_LINK_MED_ZOOM);

    fprintf(mailTmpStream, "Subject: %s\n", mail_subject);
    fprintf(mailTmpStream, "Content-Type: text/html; charset=ISO-8859-1\n");
    fprintf(mailTmpStream, "\n");
    fprintf(mailTmpStream, "\n<html><body style=\"font-family:sans-serif;font-size:small\">\n");
    fprintf(mailTmpStream, "---------------------------------------------------------------------------\n<br>");
    fprintf(mailTmpStream, "%s    %s\n<br>", EARLY_EST_MONITOR_NAME, EARLY_EST_MONITOR_AGENCY);
    fprintf(mailTmpStream, "---------------------------------------------------------------------------\n<br>");
    fprintf(mailTmpStream, "\n<br>---- ALERT MESSAGE ----\n<br>");
    fprintf(mailTmpStream, "\n<br><strong>Automatic solution - this event has not been reviewed by a seismologist</strong>\n<br>");
    fprintf(mailTmpStream, "<strong>Automatically determined event information may be erroneous!</strong>\n<br>");
    fprintf(mailTmpStream, "\n<br>Elapsed time since event origin: %s\n<br>", latency_string);
    fprintf(mailTmpStream, "Alert trigger:%s\n<br>", mail_trigger_string);
    fprintf(mailTmpStream, "\n<br>Earthquake details:\n<br>%s\n<br>", hypo_info_string);
    if (tsunami_evaluation_write) {
        fprintf(mailTmpStream, "\n<br>%s\n<br>", tsunami_info_string);
    }
    fprintf(mailTmpStream, "\n<br>Real-time information: <a href='%s'>%s real-time display</a>", info_link, EARLY_EST_MONITOR_SHORT_NAME);
    create_event_link(info_link, phypo->unique_id, event_url, event_link, feregion_str);
    fprintf(mailTmpStream, "\n<br>Event page: %s%s event information</a>", event_link, EARLY_EST_MONITOR_SHORT_NAME);
    fprintf(mailTmpStream, "\n<br>Epicenter in Google Maps: ");
    fprintf(mailTmpStream, " <a href='%s'>Google Map</a>\n<br>", map_link);
    printHypoDataHeaderString(hypoDataString);
    printHypoDataString(phypo, hypo_info_string, 1);
    fprintf(mailTmpStream, "\n\n<br>%s\n<br>%s\n<br>", hypoDataString, hypo_info_string);
    fprintf(mailTmpStream, "Elapsed time after origin of first association/location: %s\n\n<br>", timeDiff2string((time_t) phypo->first_assoc_time, (time_t) phypo->otime, latency_string));
    fprintf(mailTmpStream, "\n<br>Reporting agency: %s\n<br>", agencyId);
    fprintf(mailTmpStream, "Message time (UTC): %s\n<br>", asctime(gmtime(&current_time)));
    fprintf(mailTmpStream, "%s software version: %s (%s)\n<br>", EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE);
    fprintf(mailTmpStream, "</body></html>\n");
    fprintf(mailTmpStream, "\n.\n");
    fclose_counter(mailTmpStream);

    if (can_send_mail) {
        sprintf(sys_command, "sendmail -f %s %s < %s &", mail_from, mail_to, mail_file);
        if (verbose > 1)
            printf("Running command: %s\n", sys_command);
        int process_status = system(sys_command);
        if (verbose > 1)
            printf("Return %d from running command: %s\n", process_status, sys_command);
    }

    if (verbose) {
        if (can_send_mail)
            printf("Info: Hypocenter alert e-mail sent: %s\n", mail_subject);
        else
            printf("Warning: Hypocenter alert e-mail not sent: error in parsing -sendmail parameters:  %s  ->  info_link: %s  mail_from:%s  mail_to:%s\n",
                sendMailParams, info_link, mail_from, mail_to);
        printf("Info: Hypocenter alert e-mail file: %s\n", mail_file);
    }

    // run external alert command
    if (1) {
        // save hypocenter to alert file
        sprintf(alert_file, "%s/alert_%ld.txt", outnameroot, phypo->unique_id);
        FILE * fp_alert = fopen_counter(alert_file, "w");
        if (fp_alert == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", alert_file);
            perror(tmp_str);
        } else {

            printHypoDataString(phypo, hypoDataString, 0);
            fprintf(fp_alert, "%s\n", hypoDataString);
            fclose_counter(fp_alert);
            // run alert script
            sprintf(sys_command, "./run_action__alert_sent.bash %s %s &", alert_file, event_url);
            //DEBUG if (verbose > 1)
            printf("Running command: %s\n", sys_command);
            int process_status = system(sys_command);
            //DEBUG if (verbose > 1)
            printf("Return %d from running command: %s\n", process_status, sys_command);
        }
    }


    // write kml
    // TODO: create correct KML...
    //write_hypocenter_kml(outnameroot, phypo, verbose);

    return (1);

}

/** calculate and write station health information */
// 20150909 AJL - added stations map

void updateStaHealthInformation(char *outnameroot, FILE* staHealthHtmlStream, double report_time_window, double report_interval, double time_since_last_report, time_t time_min, time_t time_max) {

    int n;

    // initialize work fields
    for (n = 0; n < num_sources_total; n++) {
        channelParameters[n].numData = 0;
        channelParameters[n].numDataAssoc = 0;
        channelParameters[n].qualityWeight = 1.0;
    }

    // count station deData
    TimedomainProcessingData* deData = NULL;
    int source_id = -1;
    for (n = 0; n < num_de_data; n++) {
        deData = data_list[n];
        source_id = deData->source_id;
        // 20120316 AJL - add check on flag_complete_t50 and only count assoc data when data counted
        //      avoids premature downweight of station qualityWeight from data that may later be ignored
        if (deData->flag_complete_t50 && !ignoreData(deData)) {
            channelParameters[source_id].numData++;
            if (deData->is_associated)
                channelParameters[source_id].numDataAssoc++;
        }
        // END - 20120316 AJL
    }

    ChannelParameters* chan_params = NULL;
    if (num_sources_total != num_sorted_chan_params)
        printf("Warning: num_sources_total %d != num_sorted_chan_params %d: OK if data from some sources is earlier than report window.\n", num_sources_total, num_sorted_chan_params);

    // accumulate statistics for stations
    nstaHasBeenActive = 0;
    nstaIsActive = 0;
    double data_latency_sum = 0.0;
    double data_latency_sum_sqr = 0.0;
    double data_latency;
    int num_data_unassoc_sum = 0;
    int num_data_unassoc_sum_sqr = 0;
    int num_data_unassoc = 0;
    for (n = 0; n < num_sorted_chan_params; n++) {
        chan_params = sorted_chan_params_list[n];
        // 20160906 AJL  if (chan_params->inactive_duplicate || !chan_params->process_this_channel_orientation) {
        if (chan_params->inactive_duplicate) {
            chan_params->data_latency = -1.0;
            chan_params->staActiveInReportInterval = 0;
            continue;
        }
        /* 20151117 AJL
        if (!chan_params->staActiveInReportInterval) // new station data did not arrive in report interval, increase latency by report_interval
            chan_params->data_latency += report_interval;*/
        if (!chan_params->staActiveInReportInterval) // new station data did not arrive in report interval, increase latency by report_interval
            chan_params->data_latency += time_since_last_report;
        chan_params->staActiveInReportInterval = 0;
        if (!chan_params->process_this_channel_orientation) {
            continue;
        }
        nstaHasBeenActive++;
        data_latency = chan_params->data_latency;
        data_latency_sum += data_latency;
        data_latency_sum_sqr += data_latency * data_latency;
        num_data_unassoc = chan_params->numData - chan_params->numDataAssoc;
        num_data_unassoc_sum += num_data_unassoc;
        num_data_unassoc_sum_sqr += num_data_unassoc * num_data_unassoc;
        if (data_latency < report_interval) {
            nstaIsActive++;
        }
    }
    double data_latency_mean = 0.0;
    double data_latency_variance = 999.0 * 999.0;
    double num_data_unassoc_mean = 0.0;
    double num_data_unassoc_variance = 999.0 * 999.0;
    if (nstaHasBeenActive > 0) {
        data_latency_mean = data_latency_sum / (double) nstaHasBeenActive;
        num_data_unassoc_mean = num_data_unassoc_sum / (double) nstaHasBeenActive;
        if (nstaHasBeenActive > 1) {
            data_latency_variance = (data_latency_sum_sqr - (double) nstaHasBeenActive * data_latency_mean * data_latency_mean) / ((double) nstaHasBeenActive - 1.0);
            num_data_unassoc_variance = (num_data_unassoc_sum_sqr - (double) nstaHasBeenActive * num_data_unassoc_mean * num_data_unassoc_mean) / ((double) nstaHasBeenActive - 1.0);
        }
    }
    double data_latency_std_dev = sqrt(data_latency_variance);
    double num_data_unassoc_std_dev = sqrt(num_data_unassoc_variance);

    // write output

    fprintf(staHealthHtmlStream, "<html>\n<head>\n<title>Station Health - %s</title>\n", EARLY_EST_MONITOR_NAME);
    fprintf(staHealthHtmlStream, "<meta http-equiv=\"refresh\" content=\"30\">\n</head>\n<body>\n");
    fprintf(staHealthHtmlStream, "</head>\n");
    fprintf(staHealthHtmlStream, "<body style=\"font-family:sans-serif;font-size:small\">\n");

    //fprintf(staHealthHtmlStream, "\n<table border=0 cellpadding=0 cellspacing=2>\n<tbody>\n");
    fprintf(staHealthHtmlStream, "\n<table border=0 cellpadding=1 frame=box rules=rows>\n<tbody>\n");
    fprintf(staHealthHtmlStream, "<tr align=right bgcolor=\"#BBBBBB\">\n");
    fprintf(staHealthHtmlStream, "<th>&nbsp;Stations:</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;total</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;DataLatency:</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;n<%gs</th>", report_interval);
    fprintf(staHealthHtmlStream, "<th>&nbsp;mean&nbsp;<br>&nbsp;(sec)</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;std_dev&nbsp;<br>&nbsp;(sec)</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;DataUnassoc:</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;mean</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;std_dev</th>");
    fprintf(staHealthHtmlStream, "</tr>\n");
    fprintf(staHealthHtmlStream, "<tr align=right bgcolor=\"#FFFFFF\">\n");
    fprintf(staHealthHtmlStream, "<td></td>");
    // 20141212 AJL fprintf(staHealthHtmlStream, "<td>%d</td>", num_sources_total);
    fprintf(staHealthHtmlStream, "<td>%d</td>", nstaHasBeenActive);
    fprintf(staHealthHtmlStream, "<td></td>");
    fprintf(staHealthHtmlStream, "<td>%d</td>", nstaIsActive);
    fprintf(staHealthHtmlStream, "<td>%.1f</td>", data_latency_mean);
    fprintf(staHealthHtmlStream, "<td>%.1f</td>", data_latency_std_dev);
    fprintf(staHealthHtmlStream, "<td></td>");
    fprintf(staHealthHtmlStream, "<td>%.2f</td>", num_data_unassoc_mean);
    fprintf(staHealthHtmlStream, "<td>%.2f</td>", num_data_unassoc_std_dev);
    fprintf(staHealthHtmlStream, "</tr>\n");
    fprintf(staHealthHtmlStream, "\n</tbody>\n</table>\n");

    fprintf(staHealthHtmlStream, "<br>\n");

    // create map page after station health values updated
    //create_map_html_page(outnameroot, NULL, time_min, time_max, MAP_LINK_GLOBAL_ZOOM);
    create_map_link("./", -1, map_link_str, MAP_LINK_GLOBAL_ZOOM);
    fprintf(staHealthHtmlStream, "\n<table border=0 cellpadding=1 frame=box rules=rows>\n<tbody>\n");
    fprintf(staHealthHtmlStream, "<tr align=right bgcolor=\"#FFFFFF\">\n");
    fprintf(staHealthHtmlStream, "<th><a href=\"%s\" title=\"View statopms in Google Maps\" target=map_stations>Station map</a></td>", map_link_str);
    fprintf(staHealthHtmlStream, "</tr>\n");
    fprintf(staHealthHtmlStream, "\n</tbody>\n</table>\n");

    fprintf(staHealthHtmlStream, "<br>\n");

    //fprintf(staHealthHtmlStream, "\n<table border=0 cellpadding=0 cellspacing=2>\n<tbody>\n");
    fprintf(staHealthHtmlStream, "\n<table border=0 cellpadding=1 frame=box rules=rows>\n<tbody>\n");
    fprintf(staHealthHtmlStream, "<tr bgcolor=\"#BBBBBB\">\n");
    fprintf(staHealthHtmlStream, "<th>net</th><th>sta</th><th>loc</th><th>chan</th><th>streams</th><th>&nbsp;Hz</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;lat&nbsp;<br>&nbsp;(deg)</th><th>&nbsp;lon&nbsp;<br>&nbsp;(deg)</th><th>&nbsp;elev&nbsp;<br>(m)</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;cmpAz&nbsp;<br>&nbsp;(deg)</th><th>&nbsp;cmpDip&nbsp;<br>&nbsp;(deg)</th><th>&nbsp;gain&nbsp;<br>(cts/m/s)</th><th>&nbsp;staCorr&nbsp;<br>(s)</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;latency</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;nUnassoc</th><th>&nbsp;qualityWt</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;notContig&nbsp;<br>&nbsp;level</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;conflctDt&nbsp;<br>&nbsp;level</th>");
    fprintf(staHealthHtmlStream, "<th>&nbsp;notes</th>");
    fprintf(staHealthHtmlStream, "</tr>\n");
    int ierror;
    double value;
    int latency_hour, latency_min, latency_tenths_sec;
    static char bgColor[64];
    // 20151117 AJL  double decay_factor = (report_time_window - report_interval) / report_time_window;
    double decay_factor = (report_time_window - time_since_last_report) / report_time_window;
    //double std_dev_cutoff = NUM_DATA_UNASSOC_STD_CUTOFF * num_data_unassoc_std_dev;
    int duration = 3600; // 1
    double start_time = (double) (time_max - duration);
    for (n = 0; n < num_sorted_chan_params; n++) {
        chan_params = sorted_chan_params_list[n];
        if (chan_params->inactive_duplicate || !chan_params->process_this_channel_orientation)
            strcpy(bgColor, "bgcolor=\"#EEEEEE\"");
        else
            strcpy(bgColor, "bgcolor=\"#FFFFFF\"");
        fprintf(staHealthHtmlStream, "<tr align=right %s>\n", bgColor);
        // create links for channel data
        create_channel_links(chan_params->network, chan_params->station, chan_params->location, chan_params->channel,
                STREAM_RAW, STREAM_RAW_NAME, chan_params->n_int_tseries, start_time, duration,
                tmp_str, tmp_str_2);
        fprintf(staHealthHtmlStream,
                "<td align=left>%s</td><td align=left>%s</td><td align=left>%s</td><td align=left>%s</td><td align=left>%s</td><td>%g</td><td>&nbsp;&nbsp;%.3f</td><td>%.3f</td><td>%.0f</td>",
                chan_params->network, chan_params->station, chan_params->location,
                tmp_str, tmp_str_2, 1.0 / chan_params->deltaTime,
                chan_params->lat, chan_params->lon, chan_params->elev);
        int source_id = chan_params - channelParameters;
        fprintf(staHealthHtmlStream,
                "<td>%.1f</td><td>%.1f</td><td>&nbsp;&nbsp;%g</td><td>%.2f</td>",
                chan_params->azimuth, chan_params->dip,
                chan_resp[source_id].gain, chan_params->sta_corr->mean);
        // latency
        value = chan_params->data_latency;
        if (chan_params->inactive_duplicate)
            ;
        else if (value >= LATENCY_RED_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFBBBB\"");
        } else if (value >= LATENCY_YELLOW_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFFFCC\"");
        } else {
            strcpy(bgColor, "bgcolor=\"#DDFFDD\"");
        }
        latency_tenths_sec = (int) (0.5 + 10.0 * value);
        latency_hour = latency_min = 0;
        while (latency_tenths_sec >= 36000) {
            latency_tenths_sec -= 36000;
            latency_hour++;
        }
        while (latency_tenths_sec >= 600) {
            latency_tenths_sec -= 600;
            latency_min++;
        }
        if (latency_hour > 0)
            sprintf(chan_params->data_latency_str, "%dh%2.2dm%2.2ds", latency_hour, latency_min, latency_tenths_sec / 10);
        else if (latency_min > 0)
            sprintf(chan_params->data_latency_str, "%dm%2.2ds", latency_min, latency_tenths_sec / 10);
        else if (latency_tenths_sec > 0)
            sprintf(chan_params->data_latency_str, "%ds", latency_tenths_sec / 10);
        else
            sprintf(chan_params->data_latency_str, "0s");
        fprintf(staHealthHtmlStream, "<td %s>%s</td>", bgColor, chan_params->data_latency_str);
        // quality
        /*
        double num_data_unassoc_value = (double) (chan_params->numData - chan_params->numDataAssoc - DATA_UNASSOC_CUTOFF);
        if (num_data_unassoc_value <= 0.0)
            chan_params->qualityWeight = 1.0;
        else
            chan_params->qualityWeight = 1.0 / (num_data_unassoc_value + 1.0);
         */
        //
        if (chan_params->numData - DATA_UNASSOC_CUTOFF > 0) {
            chan_params->qualityWeight = (double) (chan_params->numDataAssoc + 1) / (double) (chan_params->numData - DATA_UNASSOC_CUTOFF + 1);
            //printf("DEBUG: chan_params->qualityWeight %f = (double) (chan_params->numDataAssoc %d + 1) / (double) (chan_params->numData %d - DATA_UNASSOC_CUTOFF + 1)\n",
            //        chan_params->qualityWeight, chan_params->numDataAssoc, chan_params->numData);
            if (chan_params->qualityWeight > 1.0)
                chan_params->qualityWeight = 1.0;
        } else {
            chan_params->qualityWeight = 1.0;
        }
        value = chan_params->qualityWeight;
        //
        if (chan_params->inactive_duplicate || !chan_params->process_this_channel_orientation)
            ;
        else if (value < DATA_UNASSOC_WT_RED_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFBBBB\"");
        } else if (value < DATA_UNASSOC_WT_YELLOW_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFFFCC\"");
        } else {
            strcpy(bgColor, "bgcolor=\"#DDFFDD\"");
        }
        fprintf(staHealthHtmlStream, "<td %s>%d</td><td %s>%.2f</td>",
                bgColor, chan_params->numData - chan_params->numDataAssoc, bgColor, chan_params->qualityWeight);
        // errors
        double level_non_contiguous = chan_params->level_non_contiguous * decay_factor;
        level_non_contiguous += chan_params->count_non_contiguous;
        if (chan_params->inactive_duplicate)
            ;
        else if (level_non_contiguous >= LEVEL_NON_CONTIGUOUS_RED_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFBBBB\"");
        } else if (level_non_contiguous >= LEVEL_NON_CONTIGUOUS_YELLOW_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFFFCC\"");
        } else {
            strcpy(bgColor, "bgcolor=\"#DDFFDD\"");
        }
        fprintf(staHealthHtmlStream, "<td %s>%.1f</td>", bgColor, level_non_contiguous);
        chan_params->level_non_contiguous = level_non_contiguous;
        //
        double level_conflicting_dt = chan_params->level_conflicting_dt * decay_factor;
        level_conflicting_dt += chan_params->count_conflicting_dt;
        if (chan_params->inactive_duplicate)
            ;
        else if (level_conflicting_dt >= LEVEL_CONFLICTING_DT_RED_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFBBBB\"");
        } else if (level_conflicting_dt >= LEVEL_CONFLICTING_DT_YELLOW_CUTOFF) {
            strcpy(bgColor, "bgcolor=\"#FFFFCC\"");
        } else {
            strcpy(bgColor, "bgcolor=\"#DDFFDD\"");
        }
        fprintf(staHealthHtmlStream, "<td %s>%.1f</td>", bgColor, level_conflicting_dt);
        chan_params->level_conflicting_dt = level_conflicting_dt;
        //
        ierror = chan_params->error;
        if (ierror > 0 || chan_params->inactive_duplicate) {
            if (chan_params->inactive_duplicate) {
                fprintf(staHealthHtmlStream, "<td %s>inactive_duplicate", bgColor);
                //} else if (!chan_params->process_this_channel_orientation) {
                //    fprintf(staHealthHtmlStream, "<td %s>other_orientation ", bgColor);
            } else {
                strcpy(bgColor, "bgcolor=\"#FFBBBB\"");
                fprintf(staHealthHtmlStream, "<td %s>", bgColor);
            }
            if (ierror & ERROR_DIFFERENT_SAMPLE_RATES) {
                fprintf(staHealthHtmlStream, "conflicting_dt");
            }
            if (ierror & ERROR_DATA_NON_CONTIGUOUS) {
                fprintf(staHealthHtmlStream, "data_not_contig");
            }
            if (ierror & ERROR_DT_NOT_SUPPORTED_BY_FILTER)
                fprintf(staHealthHtmlStream, "dt_not_supported_by_filter");
            fprintf(staHealthHtmlStream, "</td>");
        } else {

            strcpy(bgColor, "bgcolor=\"#DDFFDD\"");
            fprintf(staHealthHtmlStream, "<td %s>-</td>", bgColor);
        }
        //
        fprintf(staHealthHtmlStream, "</tr>\n");
        // re-initilize station health counters
        chan_params->error = 0;
        chan_params->count_non_contiguous = 0;
        chan_params->count_conflicting_dt = 0;
    }

    fprintf(staHealthHtmlStream, "\n</tbody>\n</table>\n");


    fprintf(staHealthHtmlStream, "</body>\n</html>\n");

    // create map page after station health values updated
    create_map_html_page(outnameroot, NULL, time_min, time_max, MAP_LINK_GLOBAL_ZOOM);

}

/** evaluate and set level statistics */

int setStatistics(char *levelName, double **statisticsArray, int numLevel, statistic_level* pLevelStatistics, int verbose) {

    int n, m;

    double **statisticsArray_tmp;
    double* percentiles;

    pLevelStatistics->numLevel = numLevel;

    // trim before weighting
    int allocated_tmp = 0;
    if (TRIM_BEFORE_WEIGHTING && numLevel > 2) {
        percentiles = vector_percentiles(statisticsArray[0], numLevel, 'a');
        statisticsArray_tmp = calloc(3, sizeof (double*));
        for (m = 0; m < 3; m++) {
            statisticsArray_tmp[m] = calloc(numLevel, sizeof (double));
        }
        allocated_tmp = 1;
        int ntrim = 0;
        double value;
        for (n = 0; n < numLevel; n++) {
            value = statisticsArray[0][n];
            if (value >= percentiles[10] && value <= percentiles[90]) {
                statisticsArray_tmp[0][ntrim] = value;
                statisticsArray_tmp[1][ntrim] = statisticsArray[1][n];
                statisticsArray_tmp[2][ntrim] = statisticsArray[2][n];
                ntrim++;
            }
        }
        numLevel = ntrim;
        free(percentiles);
    } else {
        statisticsArray_tmp = statisticsArray;
    }

    // station distribution weighting
    if (USE_DISTRIBUTION_WEIGHTING) {
        double* weights = setDistributionWeights(statisticsArray_tmp[1], statisticsArray_tmp[2], numLevel);
        if (verbose) {
            int im;
            fprintf(stdout, "= %s levelStatistics ====================\n", levelName);
            for (im = 0; im < numLevel; im++) {
                fprintf(stdout, "stats.addValue(%f, new Position(%f, %f), w=%f);\n",
                        statisticsArray_tmp[0][im], statisticsArray_tmp[1][im], statisticsArray_tmp[2][im], weights[im]);
            }
            fprintf(stdout, "========================\n");
        }
        percentiles = vector_percentiles_weighted(statisticsArray_tmp[0], weights, numLevel);
        free(weights);
    } else {
        percentiles = vector_percentiles(statisticsArray_tmp[0], numLevel, 'a');
    }
    // use all values
    pLevelStatistics->centralValue = percentiles[50];
    pLevelStatistics->upperBound = percentiles[84]; // median plus 1-std-dev
    pLevelStatistics->lowerBound = percentiles[16]; // median minus 1-std-dev

    free(percentiles);
    if (allocated_tmp) {
        for (m = 0; m < 3; m++) {

            free(statisticsArray_tmp[m]);
        }
        free(statisticsArray_tmp);
    }

    return (0);

}

/** evaluate and set level statistics */
/*
int setStatisticsGeometric(char *levelName, double **statisticsArray, int numLevel, statistic_level* pLevelStatistics, int verbose) {

    int n, m;

    double **statisticsArray_tmp = calloc(3, sizeof (double*));
    for (m = 0; m < 3; m++) {
        statisticsArray_tmp[m] = calloc(numLevel, sizeof (double));
    }

    // set logarithmic values
    for (n = 0; n < numLevel; n++) {
        statisticsArray_tmp[0][n] = log10(statisticsArray[0][n]);
        statisticsArray_tmp[1][n] = statisticsArray[1][n];
        statisticsArray_tmp[2][n] = statisticsArray[2][n];
    }

    setStatistics(levelName, statisticsArray_tmp, numLevel, pLevelStatistics, verbose);

    for (m = 0; m < 3; m++) {
        free(statisticsArray_tmp[m]);
    }
    free(statisticsArray_tmp);

    // unset logarithmic values
    pLevelStatistics->centralValue = pow(10.0, pLevelStatistics->centralValue);
    pLevelStatistics->upperBound = pow(10.0, pLevelStatistics->upperBound); // median plus 1-std-dev
    pLevelStatistics->lowerBound = pow(10.0, pLevelStatistics->lowerBound); // median minus 1-std-dev

    return (0);

}
 */

/** calculated the vector sum of station locations in a flat distance-azimuth geometry for all stations less than a specified distance from a specified point
 *
 * returns number of stations used
 */

int calcStationVectorSum(ChannelParameters* sta_coordinates, int num_sta_coordinates, double lat, double lon, double distance_max,
        double *pvector_distance, double *pvector_azimuth) {

    double x_vector_sum = 0.0; // x accumulation for vector sum of epicenter to station vectors
    double y_vector_sum = 0.0; // y accumulation for vector sum of epicenter to station vectors
    int n_vector_sum = 0;

    double distance, azimuth;
    int nsta;
    for (nsta = 0; nsta < num_sta_coordinates; nsta++) {
        // 20160902 AJL - bug fix  if (sta_coordinates[nsta].have_coords && !channelParameters[nsta].inactive_duplicate && channelParameters->process_this_channel_orientation) {
        if (sta_coordinates[nsta].have_coords && !channelParameters[nsta].inactive_duplicate && channelParameters[nsta].process_this_channel_orientation) {
            ChannelParameters* coords = sta_coordinates + nsta;
            distance = GCDistance(lat, lon, coords->lat, coords->lon);
            if (distance <= distance_max) {
                azimuth = GCAzimuth(lat, lon, coords->lat, coords->lon);
                x_vector_sum += distance * sin(azimuth * DE2RA);
                y_vector_sum += distance * cos(azimuth * DE2RA);
                n_vector_sum++;
            }
        }
    }

    // set mean vector of epicenter to station vectors
    *pvector_azimuth = atan2(x_vector_sum, y_vector_sum) * RA2DE;
    if (*pvector_azimuth < 0.0)
        *pvector_azimuth += 360.0;
    *pvector_distance = sqrt(x_vector_sum * x_vector_sum + y_vector_sum * y_vector_sum);

    return (n_vector_sum);

}



// from http://stackoverflow.com/questions/5377450/maximum-number-of-open-filehandles-per-process-on-osx-and-how-to-increase

struct rlimit limit;

/*
 * Set max number of files.
 */
void setLimit(int lim) {
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {
        printf("getrlimit() failed with errno=%d\n", errno);
        exit(1);
    }
    limit.rlim_cur = lim;
    //limit.rlim_max = lim;
    printf("Info: Attempting to set maximum number of open files for this process to: soft limit: %llu, hard limit: %llu\n", limit.rlim_cur, limit.rlim_max);
    if (setrlimit(RLIMIT_NOFILE, &limit) != 0) {

        printf("setrlimit() failed with errno=%d\n", errno);
        exit(1);
    }
}

/*
 * Get max number of files.
 */
void getLimit() {
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {

        printf("getrlimit() failed with errno=%d\n", errno);
        exit(1);
    }
    printf("Info: Maximum number of open files for this process set to: soft limit: %llu, hard limit: %llu\n", limit.rlim_cur, limit.rlim_max);
}

/*
 * Calculate T50Ex level
 */
double getT50Level(TimedomainProcessingData* deData) {

    double t50Level = T50EX_LEVEL_MIN;
    if (deData->flag_complete_t50)
        if (!isnan(deData->t50) && !isnan(deData->a_ref) && deData->a_ref > FLT_MIN) // avoid divide by zero
            t50Level = deData->t50 / deData->a_ref;

    // 20140407 AJL - bug fix.
    if (isnan(t50Level)) {

        t50Level = T50EX_LEVEL_MIN;
    }

    return (t50Level);

}

/*
 * Determine if T0 should be used for reporting
 */
int useT0Report(TimedomainProcessingData* deData) {

    double distance = deData->epicentral_distance;

    if (distance > MIN_EPICENTRAL_DISTANCE_T0) {
        if (distance < MAX_EPICENTRAL_DISTANCE_T0_CONDITIONAL) {
            double t50Level = getT50Level(deData);
            if (t50Level >= T0_CONDITIONAL_T50_LIMIT && !deData->flag_snr_brb_too_low && deData->tauc_peak != TAUC_INVALID && deData->tauc_peak > T0_CONDITIONAL_TAUC_LIMIT)
                return (1);
        } else if (distance < MAX_EPICENTRAL_DISTANCE_T0) {

            return (1);
        }
    }

    return (0);
}

/***************************************************************************
 * td_timedomainProcessingReport_init:
 *
 * Do necessary initialization.
 ***************************************************************************/
int td_timedomainProcessingReport_init(
        char* outnameroot_archive
        ) {

    // if many events associated/located will open many files, increase file number limit
    /* does not seem to work on OS X 10.6
    getLimit();
    setLimit(4096);
    getLimit();
     */

    // set relevant program properties
    char out_buf[STANDARD_STRLEN];


    // get application properties settings
    double double_param;
    int int_param;

    // set phase names to use in location/association
    // default is use all (see ttimes.c)
    if (settings_get(app_prop_settings, "AssociateLocate", "phase.names.use", out_buf, STANDARD_STRLEN)) {
        num_phase_names_use = 0;
        printf("Info: property set: [AssociateLocate] phase.names.use %s -> ", out_buf);
        char *str_pos = strtok(out_buf, ", ");
        while (str_pos != NULL && num_phase_names_use < get_num_ttime_phases()) {
            *channel_names_use[num_phase_names_use] = '\0';
            *time_delay_use[num_phase_names_use] = '\0'; // 20150410 AJL - added
            char *str_pos_pha_chan = strchr(str_pos, ':');
            if (str_pos_pha_chan != NULL) { // there is a ':' separator character, should be a channel string following
                char *str_pos_pha_tdelay = strchr(str_pos_pha_chan + 1, ':');
                if (str_pos_pha_tdelay != NULL) { // there is a ':' separator character, should be a time delay string following // 20150410 AJL - added
                    strcpy(time_delay_use[num_phase_names_use], str_pos_pha_tdelay + 1);
                    *str_pos_pha_tdelay = '\0';
                }
                strcpy(channel_names_use[num_phase_names_use], str_pos_pha_chan + 1);
                *str_pos_pha_chan = '\0';
            }
            strcpy(phase_names_use[num_phase_names_use], str_pos);
            printf(" <%s:%s:%s>,", phase_names_use[num_phase_names_use], channel_names_use[num_phase_names_use], time_delay_use[num_phase_names_use]);
            num_phase_names_use++;
            str_pos = strtok(NULL, ",");
        }
        printf("\n"); // DEBUG
    }

    // associate_locate_octtree parameters
    depth_step_full = get_dist_time_depth_max();
    depth_min_full = 0.0; // range defines grid cell limits
    depth_max_full = get_dist_time_depth_max(); // range defines grid cell limits
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.depth.step")) != DBL_INVALID)
        depth_step_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.depth.min")) != DBL_INVALID)
        depth_min_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.depth.max")) != DBL_INVALID)
        depth_max_full = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.depth. step %f min %f max %f\n", depth_step_full, depth_min_full, depth_max_full);
    // check
    if (depth_max_full > get_dist_time_depth_max()) {
        printf("ERROR: Property AssociateLocate->assoc_loc.depth.max (%f km) > maximum depth in travel time tables (%f km).\n",
                depth_max_full, get_dist_time_depth_max());
        return (-1);
    }
    // 20111110 AJL double lat_step = 5.0;
    lat_step_full = 7.2; // 20111110 AJL
    lat_min_full = -90.0; // range defines grid cell limits
    lat_max_full = 90.0; // range defines grid cell limits
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lat.step")) != DBL_INVALID)
        lat_step_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lat.min")) != DBL_INVALID)
        lat_min_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lat.max")) != DBL_INVALID)
        lat_max_full = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.lat. step %f min %f max %f\n", lat_step_full, lat_min_full, lat_max_full);
    // 20111110 AJL double lon_step_smallest = 5.0; // lon step is nominal for lat = 0 (on equator), will  be adjusted by 1/cos(lat) in assoc/location routine
    lon_step_smallest_full = 7.2; // 20111110 AJL // lon step is nominal for lat = 0 (on equator), will  be adjusted by 1/cos(lat) in assoc/location routine
    lon_min_full = -180.0; // range defines grid cell limits
    lon_max_full = 180.0; // range defines grid cell limits
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lon.step")) != DBL_INVALID)
        lon_step_smallest_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lon.min")) != DBL_INVALID)
        lon_min_full = double_param;
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.lon.max")) != DBL_INVALID)
        lon_max_full = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.lon. step %f min %f max %f\n", lon_step_smallest_full, lon_min_full, lon_max_full);
    //
    nomimal_critical_oct_tree_node_size = 50.0; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.nomimal_critical_oct_tree_node_size")) != DBL_INVALID)
        nomimal_critical_oct_tree_node_size = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.nomimal_critical_oct_tree_node_size %f\n", nomimal_critical_oct_tree_node_size);
    //
    min_critical_oct_tree_node_size = 1.0; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.min_critical_oct_tree_node_size")) != DBL_INVALID)
        min_critical_oct_tree_node_size = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.min_critical_oct_tree_node_size %f\n", min_critical_oct_tree_node_size);
    //
    nominal_oct_tree_min_node_size = 5.0; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.nominal_oct_tree_min_node_size")) != DBL_INVALID)
        nominal_oct_tree_min_node_size = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.nominal_oct_tree_min_node_size %f\n", nominal_oct_tree_min_node_size);
    //
    min_weight_sum_assoc = 4.5; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.min_weight_sum_associate")) != DBL_INVALID)
        min_weight_sum_assoc = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.min_weight_sum_associate %f\n", min_weight_sum_assoc);
    //
    min_time_delay_between_picks_for_location = 25.0; // same as aref, for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.min_time_delay_between_picks_for_location")) != DBL_INVALID)
        min_time_delay_between_picks_for_location = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.min_time_delay_between_picks_for_location %f\n", min_time_delay_between_picks_for_location);
    //
    // 20160906 AJL  gap_primary_critical = 225; // for global / teleseismic monitor mode
    gap_primary_critical = 270; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.gap_primary_critical")) != DBL_INVALID)
        gap_primary_critical = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.gap_primary_critical %f\n", gap_primary_critical);
    // 20160906 AJL  gap_secondary_critical = 270; // for global / teleseismic monitor mode
    gap_secondary_critical = 315; // for global / teleseismic monitor mode
    if ((double_param = settings_get_double(app_prop_settings, "AssociateLocate", "assoc_loc.gap_secondary_critical")) != DBL_INVALID)
        gap_secondary_critical = double_param;
    printf("Info: property set: [AssociateLocate] assoc_loc.gap_secondary_critical %f\n", gap_secondary_critical);

    // associate/locate station corrections (used here and in timedomain_processing.c)
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "sta_corr.nmin", &sta_corr_min_num_to_use, MIN_NUM_STATION_CORRECTIONS_USE,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_helper(app_prop_settings,
            "AssociateLocate", "sta_corr.filename", sta_corr_filename, sizeof (sta_corr_filename), "",
            verbose) == 0) {
        ; // handle error
    }
    use_station_corrections = 0;
    if (strlen(sta_corr_filename) > 0) {
        use_station_corrections = 1;
    }

    // associate/locate station misc settings
    //
    // upweight_picks_sn_cutoff declared in location.h
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.upweight_picks_sn_cutoff", &upweight_picks_sn_cutoff, LOCATION_UPWEIGHT_HIGH_SN_PICKS_SN_CUTOFF_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    // upweight_picks_dist_max declared in location.h
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.upweight_picks_dist_full", &upweight_picks_dist_full, LOCATION_UPWEIGHT_HIGH_SN_PICKS_DIST_MIN_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    // declared in location.h
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.upweight_picks_dist_max", &upweight_picks_dist_max, LOCATION_UPWEIGHT_HIGH_SN_PICKS_DIST_MAX_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    // use_amplitude_attenuation declared in location.h
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.use_amplitude_attenuation", &use_amplitude_attenuation, LOCATION_USE_AMPLITUDE_ATTENUATION_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    //
    initAssociateLocateParameters(upweight_picks_sn_cutoff, upweight_picks_dist_max, upweight_picks_dist_full, use_amplitude_attenuation);

    // event persistence
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.persistence.min_num_def_phases", &event_persistence_min_num_defining_phases, EVENT_PERSISTENCE_MIN_NUM_DEFINING_PHASES_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    use_event_persistence = (event_persistence_min_num_defining_phases > 0);
    //
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.persistence.frac_poss_assoc_cutoff", &event_persistence_frac_poss_assoc_cutoff, EVENT_PERSISTENCE_FRAC_POSS_ASSOC_CUTOFF_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.persistence.tt_err_fact", &event_persistence_tt_err_fact, EVENT_PERSISTENCE_TT_ERR_FACTOR_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }


    // existing event association location  // 20160915 AJL - added
    //
    if (settings_get_double_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.existing_event_assoc.delay_otime.min", &existing_assoc_delay_otime_min, EXISTING_EVENT_ASSOC_DELAY_OTIME_MIN_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.existing_event_assoc.min_num_def_phases", &existing_assoc_min_num_defining_phases, EXISTING_EVENT_ASSOC_MIN_NUM_DEFINING_PHASES,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    use_existing_assoc = (existing_assoc_min_num_defining_phases > 0);

    //
    // use_reference_phase_ttime_error
    int use_reference_phase_ttime_error = 1; // if = 1 use reference_phase_ttime_error, use 0 to disable reference_phase_ttime_error weighting
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.use_reference_phase_ttime_error", &use_reference_phase_ttime_error, LOCATION_USE_REFERENCE_PHASE_TIME_ERROR_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    if (!use_reference_phase_ttime_error) {
        reference_phase_ttime_error = -1.0;
    }


    // magnitude corrections and checks
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.use_mwp_distance_correction", &use_mwp_distance_correction, LOCATION_USE_MWP_DISTANCE_CORRECTION_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.use_mb_correction", &use_mb_correction, LOCATION_USE_MB_CORRECTION_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "AssociateLocate", "assoc_loc.use_magnitude_amp_atten_check", &use_magnitude_amp_atten_check, LOCATION_USE_MAGNITUDE_AMP_ATTEN_CHECK_DEFAULT,
            verbose
            ) == INT_INVALID) {
        ; // handle error
    }


    // enable and configure waveform export (miniseed data files) for specified time window around associated P picks.
    waveform_export_enable = 0;
    if ((int_param = settings_get_int(app_prop_settings, "WaveformExport", "enable")) != INT_INVALID)
        waveform_export_enable = int_param;
    printf("Info: property set: [WaveformExport] enable %d\n", waveform_export_enable);
    //
    if (settings_get_int_helper(app_prop_settings,
            "WaveformExport", "file_archive.age_max", &waveform_export_file_archive_age_max, WAVEFORM_EXPORT_FILE_ARCHIVE_AGE_MAX,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "WaveformExport", "memory.sliding_window.length", &waveform_export_memory_sliding_window_length, WAVEFORM_EXPORT_MEMORY_SLIDING_WINDOW_LENGTH,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "WaveformExport", "window.start.before.P", &waveform_export_window_start_before_P, WAVEFORM_EXPORT_WINDOW_START_BEFORE_P,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_int_helper(app_prop_settings,
            "WaveformExport", "window.end.after.S", &waveform_export_window_end_after_S, WAVEFORM_EXPORT_WINDOW_END_AFTER_S,
            verbose) == INT_INVALID) {
        ; // handle error
    }

    //

    // report
    if (settings_get_int_helper(app_prop_settings,
            "Report", "warning.colors.show", &warning_colors_show, WARNING_COLORS_SHOW_DEFAULT,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    if (settings_get_int_helper(app_prop_settings,
            "Report", "tsunami.evaluation.write", &tsunami_evaluation_write, TSUNAMI_EVALUATION_WRITE_DEFAULT,
            verbose) == INT_INVALID) {
        ; // handle error
    }
    /*
    if (settings_get_int_helper(app_prop_settings,
            "Report", "magnitude.colors.show", &magnitude_colors_show, MAGNITUDE_COLORS_SHOW_DEFAULT,
            verbose) == INT_INVALID) {
        ; // handle error
    }*/
    magnitude_colors_show = warning_colors_show;
    //
    if (settings_get_double_helper(app_prop_settings,
            "Report", "alert.mb_min", &mb_min_mail, MB_MIN_MAIL_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_double_helper(app_prop_settings,
            "Report", "alert.mwp_min", &mwp_min_mail, MWP_MIN_MAIL_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_double_helper(app_prop_settings,
            "Report", "alert.mwpd_min", &mwpd_min_mail, MWPD_MIN_MAIL_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_double_helper(app_prop_settings,
            "Report", "alert.tdt50ex_min", &tdt50ex_min_mail, TDT50EX_MIN_MAIL_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    //
    if (settings_get_double_helper(app_prop_settings,
            "Report", "alert.resend_time_delay", &alert_resend_time_delay_mail, ALERT_RESEND_TIME_DELAY_DEFAULT,
            verbose
            ) == DBL_INVALID) {
        ; // handle error
    }
    messageTriggerThresholdInitial.mb = mb_min_mail;
    messageTriggerThresholdInitial.mwp = mwp_min_mail;
    messageTriggerThresholdInitial.mwpd = mwpd_min_mail;
    messageTriggerThresholdInitial.alarm = tdt50ex_min_mail;
    messageTriggerThresholdInitial.resend_time_delay = alert_resend_time_delay_mail;

    // initialize list of id's for phase types to use in location/association, set reference phase ttime error
    if (numPhaseTypesUse < 0) { //  not initialized
        numPhaseTypesUse = 0;
        //printf("DEBUG: numPhaseTypesUse %d: ", numPhaseTypesUse);
        //printf("DEBUG: num_phase_names_use %d: ", num_phase_names_use);
        int j;
        if (num_phase_names_use > 0) {
            for (j = 0; j < num_phase_names_use; j++) {
                int phase_id_use = phase_id_for_name(phase_names_use[j]);
                //printf("DEBUG: phase.names.use: <%s>->%d:", phase_names_use[j], phase_id_use);
                if (phase_id_use >= 0 && phase_id_use < get_num_ttime_phases()) {
                    phaseTypesUse[numPhaseTypesUse] = phase_id_use;
                    strcpy(channelNamesUse[numPhaseTypesUse], channel_names_use[j]);
                    // 20150410 AJL - added time delay check
                    timeDelayUse[numPhaseTypesUse][0] = timeDelayUse[numPhaseTypesUse][1] = -1.0;
                    char *str_pos = NULL;
                    if (strlen(time_delay_use[j]) > 2 && (str_pos = strchr(time_delay_use[j], '$')) != NULL) { // have min/max time delay
                        timeDelayUse[numPhaseTypesUse][0] = strtod(time_delay_use[j], NULL); // min delay
                        timeDelayUse[numPhaseTypesUse][1] = strtod(str_pos + 1, NULL); // max delay
                    }
                    //printf(" %s ch:%s tdelay:<%s>%f->%f\n", phase_names_use[j], channelNamesUse[numPhaseTypesUse], // DEBUG
                    //        time_delay_use[j],
                    //        timeDelayUse[numPhaseTypesUse][0], timeDelayUse[numPhaseTypesUse][1]);
                    numPhaseTypesUse++;
                    if (use_reference_phase_ttime_error && get_phase_ttime_error(phase_id_use) < reference_phase_ttime_error)
                        reference_phase_ttime_error = get_phase_ttime_error(phase_id_use);
                } else {
                    fprintf(stderr, "ERROR: phase name to use in location not found: %s\n", phase_names_use[j]);
                }
            }
        } else { // set reference phase ttime error for all phases
            for (j = 0; j < get_num_ttime_phases(); j++) {
                if (use_reference_phase_ttime_error && get_phase_ttime_error(j) < reference_phase_ttime_error)
                    reference_phase_ttime_error = get_phase_ttime_error(j);
            }
        }
        //printf("\n"); // DEBUG
    }

    int nhyp;

    // allocate work hypocenters for associate/locate
    hyp_assoc_loc = calloc(MAX_NUM_HYPO, sizeof (HypocenterDesc*));
    for (nhyp = 0; nhyp < MAX_NUM_HYPO; nhyp++) {
        hyp_assoc_loc[nhyp] = new_HypocenterDesc();
        // initialize hypocenter arrays for scatter sample
        hyp_assoc_loc[nhyp]->scatter_sample = (float*) calloc(4 * NUM_LOC_SCATTER_SAMPLE, sizeof (float));
        //printf("DEBUG: calloc hyp_assoc_loc[nhyp]->scatter_sample: %X\n", hyp_assoc_loc[nhyp]->scatter_sample);
        hyp_assoc_loc[nhyp]->nscatter_sample = 0;
    }
    // allocate scatter sample arrays for associate/locate
    assoc_scatter_sample = calloc(MAX_NUM_HYPO, sizeof (float*));
    for (nhyp = 0; nhyp < MAX_NUM_HYPO; nhyp++) {
        assoc_scatter_sample[nhyp] = NULL;
        n_alloc_scatter_sample[nhyp] = 0;
        n_assoc_scatter_sample[nhyp] = 0;
        global_max_nassociated_P_lat_lon[nhyp] = -1.0;
    }


    // read previously archived events from existing hypolist.csv file
    char outname[1024];
    sprintf(outname, "%s/hypolist.csv", outnameroot_archive);
    FILE* hypoListStream = fopen_counter(outname, "r"); // NOTE: opened for reading
    if (hypoListStream == NULL) {
        sprintf(tmp_str, "Warning: opening output file: %s", outname);
        perror(tmp_str);
        printf("Info: file may not yet have been created: %s\n", outname);
        return (0);
    }

    // read header line
    fgets(tmp_str, STANDARD_STRLEN - 1, hypoListStream);
    //"event_id assoc_ndx ph_assoc ph_used dmin(deg) gap1(deg) gap2(deg) atten sigma_otime(sec) otime(UTC) lat(deg) lon(deg) errH(km) depth(km)"
    //        " errZ(km) Q T50Ex n Td(sec) n TdT50Ex WL_col mb n Mwp n T0(sec) n Mwpd n region n_sta_tot n_sta_active assoc_latency"
    HypocenterDesc *phypo = new_HypocenterDesc();
    int icheck_duplicates = 1;
    //char time_str[64];
    long first_assoc_latency = -1;
    int istat = 0;
    while (1) {
        if (fgets(tmp_str, STANDARD_STRLEN - 1, hypoListStream) == NULL)
            break;
        istat = readHypoDataString(phypo, tmp_str, &first_assoc_latency);
        //if (istat < 29) { // format error
        if (istat < 31) { // format error  TODO: use this line when all files have phypo->nstaHasBeenActive and phypo->nstaIsActive
            //if (istat < 32) { // format error  TODO: [20160905] use this line when all files have phypo->loc_seq_num
#ifdef USE_MWP_MO_POS_NEG
            //if (istat < 34) { // format error  TODO: [20161227] use this line when all files have phypo->mwpdMoPosNegLevelStatistics
#endif
            //printf("\nDEBUG: hypo read error: istat=%d\n", istat);
            // try to read to end of line
            if (fgets(tmp_str, STANDARD_STRLEN - 1, hypoListStream) == NULL)
                break;
            continue;
        }
        phypo->hyp_assoc_index = -1; // is not actively associated if read from hypolist
        //phypo->otime = string2timeDecSec(time_str);
        if (first_assoc_latency > 0) {
            phypo->first_assoc_time = first_assoc_latency + (long) (0.5 + phypo->otime);
        } else {
            phypo->first_assoc_time = LONG_MIN / 2;
        }
        phypo->messageTriggerThreshold = messageTriggerThresholdInitial;
        phypo->messageTriggerThreshold.mb = phypo->mbLevelStatistics.centralValue + MB_MIN_MAIL_INCREMENT;
        phypo->messageTriggerThreshold.mwp = phypo->mwpLevelStatistics.centralValue + MWP_MIN_MAIL_INCREMENT;
        phypo->messageTriggerThreshold.mwpd = phypo->mwpLevelStatistics.centralValue + MWPD_MIN_MAIL_INCREMENT;
        phypo->messageTriggerThreshold.alarm = phypo->tdT50ExLevelStatistics.centralValue + TDT50EX_MIN_MAIL_INCREMENT;
        // DEBUG
        //printHypoDataString(phypo, hypoDataString, 1);
        //printf("\nDEBUG: hypo read: %s\n", hypoDataString);
        //
        HypocenterDesc* phypocenter_desc_inserted = NULL;
        addHypocenterDescToHypoList(phypo, &hypo_list, &num_hypocenters, icheck_duplicates, NULL, &phypocenter_desc_inserted); // hypocenter unique_id is set here
    }
    free(phypo);

    fclose_counter(hypoListStream);

    return (0);

}

/***************************************************************************
 * td_timedomainProcessingReport_cleanup:
 *
 * Do necessary cleanup.
 ***************************************************************************/
void td_timedomainProcessingReport_cleanup() {

    // free work hypocenters for associate/locate
    int nhyp;
    for (nhyp = 0; nhyp < MAX_NUM_HYPO; nhyp++) {
        //printf("DEBUG: free hyp_assoc_loc[nhyp]->scatter_sample: %X\n", hyp_assoc_loc[nhyp]->scatter_sample);

        free(hyp_assoc_loc[nhyp]->scatter_sample);
        free(hyp_assoc_loc[nhyp]);
    }
    free(hyp_assoc_loc);

    // free reference full search oct-tree
    if (pOctTree != NULL) {

        int freeDataPointer = 1;
        freeTree3D(pOctTree, freeDataPointer);
        pOctTree = NULL;
    }

    /*
    int ilat, n;
    for (ilat = 0; ilat < nlat_alloc_coarse; ilat++) {
        for (n = 0; n < nlon_alloc_coarse; n++) {
            free((weight_count_grid_coarse)[ilat][n]);
        }
        free((weight_count_grid_coarse)[ilat]);
    }
    free(weight_count_grid_coarse);

    for (ilat = 0; ilat < nlat_alloc; ilat++) {
        for (n = 0; n < nlon_alloc; n++) {

            free((weight_count_grid)[ilat][n]);
        }
        free((weight_count_grid)[ilat]);
    }
    free(weight_count_grid);
     */

}

/***************************************************************************
 * associate_locate_octtree:
 *
 * Do association and location using an oct-tree search
 ***************************************************************************/
int associate_locate_octtree(int num_pass, char *outnameroot, HypocenterDesc **hyp_assoc_loc, int num_hypocenters_associated, ChannelParameters* channelParameters,
        int reassociate_only, time_t time_min, time_t time_max, double lat_min, double lat_max, double lat_step, double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step, int is_local_search, int verbose) {


    // do association/location search

    HypocenterDesc* hyp_work = hyp_assoc_loc[num_hypocenters_associated];
    double wt_count_assoc = 0.0;
    int i_get_assoc_scatter_sample = 1;

    if (reassociate_only) {
        i_get_assoc_scatter_sample = 0;
    } else {
        hyp_work->linRegressPower.power = -9.9;
        if (!is_local_search) { // 20160919 AJL
            hyp_work->ot_std_dev = 0.0;
            hyp_work->otime = 0.0;
        }
        hyp_work->lat = 0.0;
        hyp_work->lon = 0.0;
        hyp_work->depth = 0.0;
        hyp_work->nassoc_P = 0;
        hyp_work->nassoc = 0;
        hyp_work->dist_min = -1.0;
    }
    wt_count_assoc = octtreeGlobalAssociationLocation(num_pass, min_weight_sum_assoc, MAX_NUM_OCTTREE_NODES,
            nomimal_critical_oct_tree_node_size, min_critical_oct_tree_node_size, nominal_oct_tree_min_node_size,
            gap_primary_critical, gap_secondary_critical,
            lat_min, lat_max, lat_step, lon_min, lon_max, lon_step_smallest,
            depth_min, depth_max, depth_step, is_local_search,
            numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
            reference_phase_ttime_error,
            data_list, num_de_data,
            hyp_work,
            &assoc_scatter_sample[num_hypocenters_associated], &n_alloc_scatter_sample[num_hypocenters_associated], i_get_assoc_scatter_sample,
            &n_assoc_scatter_sample[num_hypocenters_associated], &global_max_nassociated_P_lat_lon[num_hypocenters_associated],
            channelParameters, reassociate_only, time_min, time_max);
    printf(" -> %d n=%d/%d/%.1f/%g ",
            num_pass, hyp_work->nassoc, hyp_work->nassoc_P, wt_count_assoc, hyp_work->prob);
    printf("a=%1f ",
            hyp_work->linRegressPower.power
            );
    printf("ot_sd=%.1f ot=%s lat/lon=%.1f/%.1f+/-%.0f depth=%.0f+/-%.0f dist=%.1f->%.1f gap=%.1f/%.1f -> %s",
            hyp_work->ot_std_dev, timeDecSec2string(hyp_work->otime, tmp_str, DEFAULT_TIME_FORMAT),
            hyp_work->lat, hyp_work->lon, hyp_work->errh, hyp_work->depth, hyp_work->errz,
            hyp_work->dist_min, hyp_work->dist_max, hyp_work->gap_primary, hyp_work->gap_secondary,
            wt_count_assoc < min_weight_sum_assoc ? "NOT ASSOCIATED" : feregion(hyp_work->lat, hyp_work->lon, feregion_str, FEREGION_STR_SIZE)
            );
    // 20150904 AJL - added ABCD hypo quality level
    double critical_errh = 10.0 * nominal_oct_tree_min_node_size; // TODO: add specific program property for critical_errh
    double critical_errz = 20.0 * nominal_oct_tree_min_node_size; // TODO: add specific program property for critical_errz
    setHypocenterQuality(hyp_work, min_weight_sum_assoc, critical_errh, critical_errz);
    printf("\n      (qual=aso%.2f/otv%.2f/erh%.2f/erz%.2f/att%.2f/gap%.2f/dcl%.2f/dfr%.2f -> %.2f %s)",
            //printf("\n      (qual=aso%.2f/otv%.2f/erh%.2f/erz%.2f/dep%.2f/att%.2f/gap%.2f/dcl%.2f/dfr%.2f -> %.2f %s)",   // 20160509 AJL - added
            hyp_work->qualityIndicators.wt_count_assoc_weight, hyp_work->qualityIndicators.ot_variance_weight,
            hyp_work->qualityIndicators.errh_weight, hyp_work->qualityIndicators.errz_weight,
            //hyp_work->qualityIndicators.depth_weight,   // 20160509 AJL - added
            hyp_work->qualityIndicators.amp_att_weight, hyp_work->qualityIndicators.gap_weight,
            hyp_work->qualityIndicators.distanceClose_weight, hyp_work->qualityIndicators.distanceFar_weight,
            hyp_work->qualityIndicators.quality_weight, hyp_work->qualityIndicators.quality_code);
    printf("\n");
    //printf("   --->    nlat=%d nlon=%d  nlat_alloc=%d nlon_alloc=%d\n", nlat, nlon, nlat_alloc, nlon_alloc);

    // 20140721 AJL      if (wt_count_assoc < min_weight_sum_assoc) {
    if (!reassociate_only && wt_count_assoc < min_weight_sum_assoc) {
        int n;
        for (n = 0; n < num_de_data; n++) {
            TimedomainProcessingData* deData = data_list[n];
            if (deData->is_associated == num_pass)
                deData->is_associated = 0;
        }
        return (0);
    }

    // accepted

    // plot search association weight count grid
    //double min_plot = 0.90 * global_max_nassociated_P_lat_lon[num_hypocenters_associated] - 1.0;
    double min_plot = 0.1 * global_max_nassociated_P_lat_lon[num_hypocenters_associated];
    if (min_plot < 0.0)
        min_plot = 0.0;
    int PLOT_ALL_SAMPLED_GRID_NODES = 0; // DEBUG: shows all visited nodes
    if (PLOT_ALL_SAMPLED_GRID_NODES)
        min_plot = -0.5; // -1.0 flags unused weight count grid node
    char outname[1024];
    // plotting files
    sprintf(outname, "%s/plot", outnameroot);
    mkdir(outname, 0755);
    sprintf(outname, "%s/plot/hypo.%d.grid.xy", outnameroot, num_pass - 1);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypocenterGridStream = fopen_counter(outname, "w");
    if (hypocenterGridStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (0);
    }
    int i, index;
    double sym_size;
    for (i = 0; i < n_assoc_scatter_sample[num_hypocenters_associated]; i++) {
        index = 4 * i;
        double weight_count = assoc_scatter_sample[num_hypocenters_associated][index + 3];
        if (weight_count > min_plot) { // -1.0 flags unused weight count grid node
            if (PLOT_ALL_SAMPLED_GRID_NODES) {
                sym_size = 0.1;
            } else {
                if (weight_count > global_max_nassociated_P_lat_lon[num_hypocenters_associated])
                    weight_count = global_max_nassociated_P_lat_lon[num_hypocenters_associated];
                //sym_size = 0.05 + 0.05 * (weight_count - min_plot) / (global_max_nassociated_P_lat_lon[num_hypocenters_associated] - min_plot);
                //sym_size *= 2.0; // AJL 20100915
                //sym_size = 0.025; // TEST/DEBUG
                sym_size = 0.0125; // 20160916 AJL - GMT5
            }
            fprintf(hypocenterGridStream, "%f %f %f\n", assoc_scatter_sample[num_hypocenters_associated][index + 1], assoc_scatter_sample[num_hypocenters_associated][index], sym_size);
        }
    }
    fclose_counter(hypocenterGridStream);

    // accepted

    return (1);

}

/***************************************************************************
 * setReducedAssocLocSearchVolume:
 *
 * Set reduced search volume around specified hypocenter
 ***************************************************************************/
void setReducedAssocLocSearchVolume(HypocenterDesc *phypo, double *potime_hypo, double *pot_std_dev_hypo, double *plat_min_hypo, double *plat_max_hypo, double *plat_step_hypo,
        double *plon_min_hypo, double *plon_max_hypo, double *plon_step_smallest_hypo,
        double *pdepth_min_hypo, double *pdepth_max_hypo, double *pdepth_step_hypo) {

    // initialize reference full search oct-tree if necessary
    if (pOctTree == NULL) {
        pOctTree = setUpOctTree(lat_min_full, lat_max_full, lat_step_full,
                lon_min_full, lon_max_full, lon_step_smallest_full,
                depth_min_full, depth_max_full, depth_step_full);
    }

    // get limits of full search oct-tree
    Vect3D coords;
    coords.x = phypo->lon;
    coords.y = phypo->lat;
    coords.z = phypo->depth;
    //printf("DEBUG: setReducedAssocLocSearchVolume: full: %f-%f %f-%f %f-%f", lat_min_full, lat_max_full, lon_min_full, lon_max_full, depth_min_full, depth_max_full);
    //printf("DEBUG: setReducedAssocLocSearchVolume: coords: %f %f %f", coords.y, coords.x, coords.z);
    // check ranges
    if (coords.y < lat_min_full) {
        coords.y = lat_min_full;
    }
    if (coords.y >= lat_max_full) {
        coords.y = lat_max_full;
    }
    int lon_wrapped = fabs(lon_max_full - lon_min_full - 360.0) < 0.001 ? 1 : 0;
    if (coords.x < lon_min_full) {
        if (!lon_wrapped) {
            coords.x = lon_min_full;
        } else {
            coords.x += 360.0;
        }
    }
    if (coords.x >= lon_max_full) {
        if (!lon_wrapped) {
            coords.x = lon_max_full;
        } else {
            coords.x -= 360.0;
        }
    }
    if (coords.z < depth_min_full) {
        coords.z = depth_min_full;
    }
    if (coords.z >= depth_max_full) {
        coords.z = depth_max_full;
    }
    //
    OctNode* poct_node = getTreeNodeContaining(pOctTree, coords);
    //printf(" -> %f %f %f, poct_node %ld\n", coords.y, coords.x, coords.z, poct_node);

#define N_FULL_CELLS_HALF_WIDTH 1.5
    *plat_min_hypo = poct_node->center.y - poct_node->ds.y * N_FULL_CELLS_HALF_WIDTH;
    if (*plat_min_hypo < lat_min_full) {
        *plat_min_hypo = lat_min_full;
    }
    *plat_max_hypo = poct_node->center.y + poct_node->ds.y * N_FULL_CELLS_HALF_WIDTH;
    if (*plat_max_hypo > lat_max_full) {
        *plat_max_hypo = lat_max_full;
    }
    *plat_step_hypo = poct_node->ds.y;

    *plon_min_hypo = poct_node->center.x - poct_node->ds.x * N_FULL_CELLS_HALF_WIDTH;
    if (!lon_wrapped && *plon_min_hypo < lon_min_full) {
        *plon_min_hypo = lon_min_full;
    }
    *plon_max_hypo = poct_node->center.x + poct_node->ds.x * N_FULL_CELLS_HALF_WIDTH;
    if (!lon_wrapped && *plon_max_hypo > lon_max_full) {
        *plon_max_hypo = lon_max_full;
    }
    *plon_step_smallest_hypo = lon_step_smallest_full;

#define N_FULL_CELLS_HALF_DEPTH 0.5
    *pdepth_min_hypo = poct_node->center.z - poct_node->ds.z * N_FULL_CELLS_HALF_DEPTH;
    if (*pdepth_min_hypo < depth_min_full) {
        *pdepth_min_hypo = depth_min_full;
        *pdepth_max_hypo = *pdepth_min_hypo + 2.0 * depth_step_full * N_FULL_CELLS_HALF_DEPTH;
    } else {
        *pdepth_max_hypo = phypo->depth + depth_step_full * N_FULL_CELLS_HALF_DEPTH;
    }
    *pdepth_max_hypo = poct_node->center.z + poct_node->ds.z * N_FULL_CELLS_HALF_DEPTH;
    if (*pdepth_max_hypo > depth_max_full) {

        *pdepth_max_hypo = depth_max_full;
    }
    *pdepth_step_hypo = depth_step_full;


    /*
     *plat_min_hypo = phypo->lat - lat_step_full * N_FULL_CELLS_HALF_WIDTH;
        if (*plat_min_hypo < lat_min_full) {
     *plat_min_hypo = lat_min_full;
        }
     *plat_max_hypo = phypo->lat + lat_step_full * N_FULL_CELLS_HALF_WIDTH;
        if (*plat_max_hypo > lat_max_full) {
     *plat_max_hypo = lat_max_full;
        }
     *plat_step_hypo = *plat_max_hypo - *plat_min_hypo;
     *plat_step_hypo /= (2.0 * N_FULL_CELLS_HALF_WIDTH);

        int lon_wrapped = fabs(lon_max_full - lon_min_full - 360.0) < 0.001 ? 1 : 0;
     *plon_min_hypo = phypo->lon - lon_step_smallest_full * N_FULL_CELLS_HALF_WIDTH;
        if (!lon_wrapped && *plon_min_hypo < lon_min_full) {
     *plon_min_hypo = lon_min_full;
        }
     *plon_max_hypo = phypo->lon + lon_step_smallest_full * N_FULL_CELLS_HALF_WIDTH;
        if (!lon_wrapped && *plon_max_hypo > lon_max_full) {
     *plon_max_hypo = lon_max_full;
        }
     *plon_step_smallest_hypo = *plon_max_hypo - *plon_min_hypo;
     *plon_step_smallest_hypo /= (2.0 * N_FULL_CELLS_HALF_WIDTH);

    #define N_FULL_CELLS_HALF_DEPTH 0.5

     *pdepth_min_hypo = phypo->depth - depth_step_full * N_FULL_CELLS_HALF_DEPTH;
        if (*pdepth_min_hypo < depth_min_full) {
     *pdepth_min_hypo = depth_min_full;
     *pdepth_max_hypo = *pdepth_min_hypo + 2.0 * depth_step_full * N_FULL_CELLS_HALF_DEPTH;
        } else {
     *pdepth_max_hypo = phypo->depth + depth_step_full * N_FULL_CELLS_HALF_DEPTH;
        }
        if (*pdepth_max_hypo > depth_max_full) {
     *pdepth_max_hypo = depth_max_full;
        }
     *pdepth_step_hypo = *pdepth_max_hypo - *pdepth_min_hypo;
     *pdepth_step_hypo /= (2.0 * N_FULL_CELLS_HALF_DEPTH);
     */


    // set otime and otime range (in *pot_std_dev_hypo)
    *potime_hypo = phypo->otime;
    //*pot_std_dev_hypo = setRefOtimeDeviation(phypo, phypo); // same range as used for location.c->isSameEvent())
    *pot_std_dev_hypo = (poct_node->ds.y * N_FULL_CELLS_HALF_WIDTH * DEG2KM) / get_velocity_model_property('P', -1.0, phypo->depth);
    *pot_std_dev_hypo /= 4.0;


    /*printf("DEBUG: setReducedAssocLocSearchVolume: otime: %f+/-%f [+/-%f], lat: %f %f %f [%f/%f], lon: %f %f %f [%f/%f], depth: %f %f %f [%f/%f]\n",
     *potime_hypo, *pot_std_dev_hypo, phypo->ot_std_dev, *plat_min_hypo, phypo->lat, *plat_max_hypo, *plat_step_hypo, lat_step_full,
     *plon_min_hypo, phypo->lon, *plon_max_hypo, *plon_step_smallest_hypo, lon_step_smallest_full,
     *pdepth_min_hypo, phypo->depth, *pdepth_max_hypo, *pdepth_step_hypo, depth_step_full
            );//*/

}


/***************************************************************************
 * td_writeTimedomainProcessingReport:
 *
 * writes timedomain-processing report to gridDataStrem and csvStream
 ***************************************************************************/

static double time_of_last_report = -1.0;

int td_writeTimedomainProcessingReport(char* outnameroot_archive, char* outnameroot,
        time_t time_min, time_t time_max,
        int only_check_for_new_event, // if = 1, use event persistence for all existing events; check for new event association,
        //    if new event do full report, otherwise do nothing (use for real-time, not for running archive data)
        int cut_at_time_max, // if = 1, remove late timedomain-processing data from list (use for running archive data, not for real-time)
        int verbose, int report_interval, int sendMail, char *sendMailParams, char *agencyId
        ) {

    if (!use_event_persistence) {
        if (only_check_for_new_event) {
            printf("ERROR: Event persistence (USE_EVENT_PERSISTENCE) required to use only_check_for_new_event.\n");
            return (-1); // event persistence required to use only_check_for_new_event
        }
    }


    int n, m;
    int nhyp;


    // set global time_max
    report_time_max = time_max;

    // set approximate time of earliest real data
    if (time_max < earliest_time)
        earliest_time = time_max;

    // remove late timedomain-processing data from list (use for running archive data, not for real-time)
    if (cut_at_time_max) {
        for (n = num_de_data - 1; n >= 0; n--) {
            TimedomainProcessingData* deData = data_list[n];
            if (deData->t_time_t > time_max) {
                removeTimedomainProcessingDataFromDataList(deData, &data_list, &num_de_data);
                free_TimedomainProcessingData(deData);
                //data_list[n] = NULL; // 20160802 AJL - memory bug fix?
            }
        }
    }
    //


    // if not using aref level check to filter picks close in time, then check for use_for_location picks too soon after another use_for_location pick for same source_id
    // 20130218 AJL - added
    if (no_aref_level_check) {
        int ndata0;
        for (ndata0 = 0; ndata0 < num_de_data; ndata0++) {
            TimedomainProcessingData* deData = data_list[ndata0];
            deData->can_use_as_location_P = 1; // provisionally set as can be used a associated location P
            int source_id = deData->source_id;
            double time_dec_sec = (double) deData->t_time_t + deData->t_decsec;
            // check other timedomain-processing data for this source_id
            int ndata;
            for (ndata = 0; ndata < ndata0; ndata++) { // data_list sorted in time, so can stop at ndata0
                TimedomainProcessingData* deData_test = data_list[ndata];
                if (deData_test->source_id != source_id || !deData_test->use_for_location)
                    continue;
                // check if other pick is within minimum window before this pick
                double time_dec_sec_test = (double) deData_test->t_time_t + deData_test->t_decsec;
                double diff = time_dec_sec - time_dec_sec_test;
                if (diff < min_time_delay_between_picks_for_location) {
                    deData->use_for_location = 0;
                    deData->can_use_as_location_P = 0;
                    break;
                }
            }
        }
    }


    // update polarization analysis for timedomain-processing data from list
    if (polarization_enable) {
        for (n = 0; n < num_de_data; n++) {
            TimedomainProcessingData* deData = data_list[n];
            //printf("DEBUG: doPolarizationAnalysis: test: %d %s_%s_%s_%s dt=%f stream=%d use=%d\n", n + 1, deData->network, deData->station, deData->location, deData->channel, deData->deltaTime, deData->pick_stream, deData->use_for_location);
            if (!deData->use_for_location && (deData->use_for_location_twin_data == NULL || !deData->use_for_location_twin_data->use_for_location))
                continue;
            td_doPolarizationAnalysis(deData, n);
        }
    }
    //




    // determine time span if not specified
    if (time_min < 0) {
        for (n = 0; n < num_de_data; n++) {
            TimedomainProcessingData* deData = data_list[n];
            double dsec;
            time_t pick_time = get_time_t(deData, &dsec);
            if (time_min < 0 || pick_time < time_min)
                time_min = pick_time;
            if (time_max < 0 || pick_time > time_max)
                time_max = pick_time;
        }
    }
    if (verbose > 2) {
        printf("time_min: %s", asctime(gmtime(&time_min))); // use gmtime since times are already UTC
        printf("time_max: %s", asctime(gmtime(&time_max))); // use gmtime since times are already UTC
        printf("time_max-time_min: %d\n", (int) (time_max - time_min));
    }

    // set plot limits
    time_t plot_time_min = time_min - 61;
    //time_t plot_time_max = time_max + 121;
    time_t plot_time_max = time_max; // for real-time plots, do not go into future!
    if (verbose > 2) {
        printf("plot_time_min: %s", asctime(gmtime(&plot_time_min))); // use gmtime since times are already UTC
        printf("plot_time_max: %s", asctime(gmtime(&plot_time_max))); // use gmtime since times are already UTC
        printf("plot_time_max-plot_time_min: %d\n", (int) (plot_time_max - plot_time_min));
    }


    int hyp_full_assoc_loc_start_index = 0;
    if (use_event_persistence) {
        if (associate_data) { // associating
            // ============================================================================================================
            // check for event persistence - keep previous location results for events with insufficient possible new defining picks
            //
            int nhyp_assoc;
            // remove any not actively associated hypocenters from hyp_assoc_loc list
            // 20151208 AJL - need to update here since events that were archived in last report (associated events with otime before analysis window)
            //      are still in hyp_assoc_loc list but associated data was removed, cannot re-associated such events
            for (nhyp_assoc = num_hypocenters_associated - 1; nhyp_assoc >= 0; nhyp_assoc--) {
                if (hyp_assoc_loc[nhyp_assoc]->hyp_assoc_index < 0) { // not associated
                    // remove hypo and pack list
                    for (nhyp = nhyp_assoc; nhyp < num_hypocenters_associated - 1; nhyp++) {
                        hyp_assoc_loc[nhyp] = hyp_assoc_loc[nhyp + 1];
                    }
                    num_hypocenters_associated--;
                }
            }
            // clear array of persistent hypocenter indices
            for (nhyp_assoc = 0; nhyp_assoc < num_hypocenters_associated; nhyp_assoc++) {
                hyp_persistent[nhyp_assoc] = 0;
            }
            // set persistent events based on type of association location
            if (only_check_for_new_event) {
                // all existing hypocenters in analysis window flagged as persistent
                for (nhyp_assoc = 0; nhyp_assoc < num_hypocenters_associated; nhyp_assoc++) {
                    double otime = hyp_assoc_loc[nhyp_assoc]->otime;
                    if ((time_t) otime <= time_min + report_interval) { // 20160913 AJL
                        // otime is before analysis window, do not preserve event to force full reassociation before event is removed
                        continue;
                    }
                    hyp_persistent[nhyp_assoc] = 1;
                    hyp_assoc_loc[nhyp_assoc]->n_possible_assoc_P = -1;
                    //printf("DEBUG: Mode USE_EVENT_PERSISTENCE+only_check_for_new_event -> ev:%d %s\n", nhyp_assoc, time2string(hyp_assoc_loc[nhyp_assoc]->otime, tmp_str));
                }
            } else {
                // check each associated hypocenter for persistence
                for (nhyp_assoc = 0; nhyp_assoc < num_hypocenters_associated; nhyp_assoc++) {
                    double otime = hyp_assoc_loc[nhyp_assoc]->otime;
                    // 20160910 AJL  if ((time_t) otime <= time_min) {
                    if ((time_t) otime <= time_min + report_interval) { // 20160910 AJL
                        // otime is before analysis window, do not preserve event to force full reassociation before event is removed
                        continue;
                    }
                    if (hyp_assoc_loc[nhyp_assoc]->nassoc_P < event_persistence_min_num_defining_phases) {
                        // too few defining phases, do not preserve event
                        continue;
                    }
                    double lat = hyp_assoc_loc[nhyp_assoc]->lat;
                    double lon = hyp_assoc_loc[nhyp_assoc]->lon;
                    double depth = hyp_assoc_loc[nhyp_assoc]->depth;
                    //printf("DEBUG: check for persistence: ev:%d %s\n", nhyp_assoc, time2string(otime, tmp_str));
                    // find new picks that could be defining picks for this event
                    int n_poss_assoc = 0;
                    double wt_assoc_full = 0.0;
                    for (n = 0; n < num_de_data; n++) {
                        TimedomainProcessingData* deData = data_list[n];
                        // check if pick is skipped for location (location.c->skipData())
                        if (skipData(deData, -1)) {
                            continue;
                        }
                        // check if pick is already associated with full association location
                        if ((deData->is_associated == nhyp_assoc + 1 && deData->is_full_assoc_loc == 1)) {
                            wt_assoc_full += deData->loc_weight;
                            continue;
                        }
                        // check each phase type
                        int phase_id;
                        for (phase_id = 0; phase_id < get_num_ttime_phases(); phase_id++) {
                            // skip if phase not counted in location
                            if (!count_in_location(phase_id, 999.9, deData->use_for_location))
                                continue;
                            // check if phase is in phase type to use list, if list exists
                            if (numPhaseTypesUse > 0) {
                                if (!isPhaseTypeToUse(deData, phase_id, numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                                        data_list, n)) {
                                    continue;
                                }
                            }
                            // get travel-time from hypocenter to station for this phase
                            double gcd = GCDistance(lat, lon, deData->lat, deData->lon);
                            double time_test = get_ttime(phase_id, gcd, depth);
                            if (time_test < 0.0) // station to hypocenter distance outside phase distance range
                                continue;
                            // check if phase arrival time minus travel-time is near event origin time
                            double residual = (double) deData->t_time_t + deData->t_decsec - time_test - otime;
                            /* DEBUG
                            printf("DEBUG: unnassoc pick residual: ev:%d %s -> %s_%s_%s_%s %s %s %fs",
                                    nhyp_assoc, time2string(otime, tmp_str),
                                    deData->network, deData->station, deData->location, deData->channel,
                                    timeDecSec2string((double) deData->t_time_t + deData->t_decsec, tmp_str_2, DEFAULT_TIME_FORMAT),
                                    phase_name_for_id(phase_id), residual);
                            // DEBUG */
                            double tt_err = get_phase_ttime_error(phase_id);
                            double tolerance = event_persistence_tt_err_fact * tt_err;
                            if (fabs(residual) <= tolerance) {
                                //printf(" <= %fs -> may assoc!\n", tolerance); // DEBUG
                                n_poss_assoc++;
                                break; // found one possible association, stop checking phases to prevent double counting of data
                            } else {
                                //printf(" > %fs\n", tolerance); // DEBUG
                            }
                        }
                    }
                    // get ratio of [n_poss_assoc = number possible defining phases not already used for full assoc/loc] to [wt_assoc_full = total weight of fully associated defining phases]
                    //    NOTE: this ratio is conservative (will overestimate true ratio of weights, since number and not eventual weights used for possible associated phases)
                    double frac_poss_assoc = (double) n_poss_assoc / wt_assoc_full;
                    //printf("DEBUG: Mode USE_EVENT_PERSISTENCE -> ev:%d %s  n_poss_assoc: %d  wt_assoc_full: %f frac_poss_assoc %f\n", nhyp_assoc, time2string(otime, tmp_str), n_poss_assoc, wt_assoc_full, frac_poss_assoc);
                    // if ratio is less than cutoff, declare hypocenter as persistent, will be re-associated without full association/location
                    if (frac_poss_assoc < event_persistence_frac_poss_assoc_cutoff) {
                        hyp_persistent[nhyp_assoc] = 1;
                    }
                    hyp_assoc_loc[nhyp_assoc]->n_possible_assoc_P = n_poss_assoc;
                }
            }
            // set persistent events
            int last_saved_persistent = -1;
            for (nhyp_assoc = 0; nhyp_assoc < num_hypocenters_associated; nhyp_assoc++) {
                if (hyp_persistent[nhyp_assoc]) {
                    last_saved_persistent++;
                    if (last_saved_persistent != nhyp_assoc) {
                        // re-index hypocenters
                        HypocenterDesc *temp = hyp_assoc_loc[last_saved_persistent];
                        hyp_assoc_loc[last_saved_persistent] = hyp_assoc_loc[nhyp_assoc]; // 20140717 AJL - shift pointers only
                        hyp_assoc_loc[nhyp_assoc] = temp;
                        // re-index data from to new persistent hypo index
                        for (n = 0; n < num_de_data; n++) {
                            TimedomainProcessingData* deData = data_list[n];
                            if (deData->is_associated == nhyp_assoc + 1) {
                                deData->is_associated = last_saved_persistent + 1; // attach to preserved hypo
                            } else if (deData->is_associated == last_saved_persistent + 1) {
                                deData->is_associated = nhyp_assoc + 1;
                            }
                        }
                    }
                }
            }
            hyp_full_assoc_loc_start_index = last_saved_persistent + 1;
            //printf("DEBUG: hyp_full_assoc_loc_start_index:%d\n", hyp_full_assoc_loc_start_index);
        }
        printf("INFO: Mode USE_EVENT_PERSISTENCE active!\n");
    }


    // ============================================================================================================
    // find sets of data with largest number of associations
    //
    int use_associated_data = 0;
    int num_hypocenters_associated_previously = num_hypocenters_associated;
#if 0
    // Test only with miniseed_process!
    printf("======\n======\n======\n======\n");
    printf("WARNING: ======> Test \"relocate existing\", use only with miniseed_process\n");
    printf("======\n======\n======\n======\n");
    num_hypocenters_associated_previously = 0;
    HypocenterDesc* phypo_test;
#if 0
    //  -> 1 n=83/72/39.4/0.479319 a=-1.912196 ot_sd=1.4 ot=2007.08.02-03:21:44.33 lat/lon=51.1/-180.0+/-12 depth=20+/-16 dist=4.0->89.1 gap=88.3/129.8 -> Andreanof Islands, Aleutian Is.
    // 2007-08-02-mw67-andreanof-islands-aleutian-is_ZNE.1186024907
    phypo_test = hyp_assoc_loc[num_hypocenters_associated_previously];
    phypo_test->lat = 52;
    phypo_test->lon = 180.0;
    phypo_test->depth = 100;
    phypo_test->unique_id = 1186024904370;
    phypo_test->first_assoc_time = phypo_test->unique_id / 1000;
    num_hypocenters_associated_previously++;
#endif
#if 1
    //  -> 1 n=74/61/51.2/0.615973 a=-1.251785 ot_sd=2.0 ot=2016.09.01-16:37:59.98 lat/lon=-37.3/178.9+/-13 depth=10+/-20 dist=1.7->49.8 gap=94.6/120.7 -> Off E. Coast of N. Island, N.Z.
    // 2016-09-01-mww71-off-e-coast-of-n-island-nz_ZNE.1472747877
    phypo_test = hyp_assoc_loc[num_hypocenters_associated_previously];
    phypo_test->lat = -37.25;
    phypo_test->lon = 179.14;
    phypo_test->depth = 17.0;
    /*phypo_test->lat = -35.3;
    phypo_test->lon = 176.9;
    phypo_test->depth = 200.0;*/
    phypo_test->unique_id = 1472747879746;
    phypo_test->first_assoc_time = phypo_test->unique_id / 1000;
    phypo_test->otime = phypo_test->unique_id / 1000.0;
    phypo_test->errh = 19;
    phypo_test->ot_std_dev = 1.9;
    num_hypocenters_associated_previously++;
#endif
#if 1
    //  -> 1 n=101/68/34.1/0.276551 a=-1.672663 ot_sd=2.0 ot=2011.03.11-05:46:21.79 lat/lon=38.1/142.7+/-15 depth=10+/-14 dist=2.8->154.4 gap=54.4/70.2 -> Near East Coast of Honshu, Japan
    // Honshu_ZNE_2011.1299822381
    phypo_test = hyp_assoc_loc[num_hypocenters_associated_previously];
    phypo_test->lat = 38.07;
    phypo_test->lon = 142.91;
    phypo_test->depth = 10.0;
    /*phypo_test->lat = 37.1;
    phypo_test->lon = 143.7;
    phypo_test->depth = 30.0;*/
    phypo_test->unique_id = 1299822379229;
    phypo_test->first_assoc_time = phypo_test->unique_id / 1000;
    num_hypocenters_associated_previously++;

#endif
#endif
    num_hypocenters_associated = 0;
    /* 20160912 AJL  for (nhyp = hyp_full_assoc_loc_start_index; nhyp < MAX_NUM_HYPO; nhyp++) {
        init_HypocenterDesc(hyp_assoc_loc[nhyp]);
    }*/
    if (associate_data) { // associate / locate
        if (verbose > 0) {
            printf("INFO: ======> Associate: %s, num_de_data=%d\n", outnameroot, num_de_data);
        }
        long associate_start_time_total = clock();
        int associated;
        double lat_min_hypo, lat_max_hypo, lat_step_hypo,
                lon_min_hypo, lon_max_hypo, lon_step_smallest_hypo,
                depth_min_hypo, depth_max_hypo, depth_step_hypo,
                otime_hypo, ot_std_dev_hypo;
        // === there are three cases for association: reassociate only (hypo preserved and fixed), relocate existing (search locally around hypo), full association location
        // 20160912 AJL - added case of relocate existing (search locally around hypo) to avoid loosing previously associated events due to failed association in early pass
        int num_pass = 0;
        // === reassociate only (hypo preserved and fixed)
        int reassociate_only = 1;
        int is_local_search = 0;
        while (num_pass < hyp_full_assoc_loc_start_index) {
            num_pass++;
            long associate_start_tim_event = clock();
            associated = associate_locate_octtree(num_pass, outnameroot, hyp_assoc_loc, num_hypocenters_associated, channelParameters,
                    reassociate_only, time_min, time_max, lat_min_full, lat_max_full, lat_step_full, lon_min_full, lon_max_full, lon_step_smallest_full,
                    depth_min_full, depth_max_full, depth_step_full, is_local_search, verbose);
            if (verbose > 0) {
                printf("INFO: event reassociate only: time = %.2f\n", (double) (clock() - associate_start_tim_event) / CLOCKS_PER_SEC);
            }
            if (associated) {
                num_hypocenters_associated++;
                hyp_assoc_loc[num_pass - 1]->loc_type = LOC_TYPE_REASSOC_ONLY;
            } else { // with reassociate_only, should always be associated, otherwise some error
                num_pass--;
                break;
            }
        }
        reassociate_only = 0;
        // === relocate existing (search locally around hypo) // 20160912 AJL - added
        is_local_search = 1;
        int nhypo_test = num_pass;
        static HypocenterDesc hypo_tmp;
        while (nhypo_test < num_hypocenters_associated_previously) {
            // if existing event time since otime less than cutoff and not enough defining phases, force full association location
            if (hyp_assoc_loc[nhypo_test]->loc_seq_num >= 0) {
                double time_since_otime = difftime(time_max, hyp_assoc_loc[nhypo_test]->otime);
                //printf("DEBUG: relocate existing %d: time_since_otime %.1f, <? existing_assoc_delay_otime_min %.1f && nassoc_P %d < existing_assoc_min_num_defining_phases %d\n",
                //        nhypo_test, time_since_otime, existing_assoc_delay_otime_min, hyp_assoc_loc[nhypo_test]->nassoc_P, existing_assoc_min_num_defining_phases);
                if (time_since_otime < existing_assoc_delay_otime_min
                        && hyp_assoc_loc[nhypo_test]->nassoc_P < existing_assoc_min_num_defining_phases) {
                    // too few defining phases, do not relocate existing event
                    nhypo_test++;
                    continue;
                }
            }
            num_pass++;
            long associate_start_tim_event = clock();
            hypo_tmp = *(hyp_assoc_loc[nhypo_test]); // used after assoc/loc to check if same event
            // set reduced associate/locate search volume centered on existing hypo
            setReducedAssocLocSearchVolume(hyp_assoc_loc[nhypo_test], &otime_hypo, &ot_std_dev_hypo, &lat_min_hypo, &lat_max_hypo, &lat_step_hypo,
                    &lon_min_hypo, &lon_max_hypo, &lon_step_smallest_hypo, &depth_min_hypo, &depth_max_hypo, &depth_step_hypo);
            long unique_id = hyp_assoc_loc[nhypo_test]->unique_id;
            // preserve persistent fields
            long first_assoc_time = hyp_assoc_loc[nhypo_test]->first_assoc_time;
            int loc_seq_num = hyp_assoc_loc[nhypo_test]->loc_seq_num;
            int has_valid_magnitude = hyp_assoc_loc[nhypo_test]->has_valid_magnitude;
            int alert_sent_count = hyp_assoc_loc[nhypo_test]->alert_sent_count;
            int alert_sent_time = hyp_assoc_loc[nhypo_test]->alert_sent_time;
            init_HypocenterDesc(hyp_assoc_loc[num_pass - 1]);
            hyp_assoc_loc[num_pass - 1]->otime = otime_hypo;
            hyp_assoc_loc[num_pass - 1]->ot_std_dev = ot_std_dev_hypo;
            associated = associate_locate_octtree(num_pass, outnameroot, hyp_assoc_loc, num_hypocenters_associated, channelParameters,
                    reassociate_only, time_min, time_max, lat_min_hypo, lat_max_hypo, lat_step_hypo, lon_min_hypo, lon_max_hypo, lon_step_smallest_hypo,
                    depth_min_hypo, depth_max_hypo, depth_step_hypo, is_local_search, verbose);
            if (verbose > 0) {
                printf("INFO: existing event associate: time = %.2f\n", (double) (clock() - associate_start_tim_event) / CLOCKS_PER_SEC);
            }
            int is_same_event = 0;
            if (associated) {
                // check that same event has been relocated (may be new event near existing event)
                hyp_assoc_loc[num_pass - 1]->unique_id = -1;
                is_same_event = isSameEvent(&hypo_tmp, hyp_assoc_loc[num_pass - 1]);
            }
            if (is_same_event) {
                num_hypocenters_associated++;
                hyp_assoc_loc[num_pass - 1]->loc_type = LOC_TYPE_RELOC_EXISTING;
                // preserve event unique id
                hyp_assoc_loc[num_pass - 1]->unique_id = unique_id;
                // preserve persistent fields
                hyp_assoc_loc[num_pass - 1]->first_assoc_time = first_assoc_time;
                hyp_assoc_loc[num_pass - 1]->loc_seq_num = loc_seq_num;
                hyp_assoc_loc[num_pass - 1]->has_valid_magnitude = has_valid_magnitude;
                hyp_assoc_loc[num_pass - 1]->alert_sent_count = alert_sent_count;
                hyp_assoc_loc[num_pass - 1]->alert_sent_time = alert_sent_time;
            } else { // with relocate existing, should usually be associated, otherwise new event found, some instability or problem, just break and use full association location
                num_pass--;
                break;
            }
            nhypo_test++;
        }
        // === full association location
        is_local_search = 0;
        while (num_pass < MAX_NUM_HYPO) {
            num_pass++;
            long associate_start_tim_event = clock();
            init_HypocenterDesc(hyp_assoc_loc[num_pass - 1]);
            associated = associate_locate_octtree(num_pass, outnameroot, hyp_assoc_loc, num_hypocenters_associated, channelParameters,
                    reassociate_only, time_min, time_max, lat_min_full, lat_max_full, lat_step_full, lon_min_full, lon_max_full, lon_step_smallest_full,
                    depth_min_full, depth_max_full, depth_step_full, is_local_search, verbose);
            if (verbose > 0) {
                printf("INFO: event full associate: time = %.2f\n", (double) (clock() - associate_start_tim_event) / CLOCKS_PER_SEC);
            }
            if (associated) {
                num_hypocenters_associated++;
                hyp_assoc_loc[num_pass - 1]->loc_type = LOC_TYPE_FULL;
            } else { // with reassociate only, non associated indicates no more successful associations, or error
                break;
            }
        }
        if (verbose > 0) {
            printf("INFO: total associate time = %.2f\n", (double) (clock() - associate_start_time_total) / CLOCKS_PER_SEC);
        }
        if (num_hypocenters_associated > 0)
            use_associated_data = 1;
    } else { // no association
        num_hypocenters_associated = 1; // TODO: test this!
    }

    // return if only_check_for_new_event and no new events associated
    /*if (only_check_for_new_event) {
        printf("DEBUG: Mode only_check_for_new_event -> nassoc:%d <=? ndxfull %d new? %d\n",
                num_hypocenters_associated, hyp_full_assoc_loc_start_index, num_hypocenters_associated <= hyp_full_assoc_loc_start_index);
    }*/
    if (only_check_for_new_event && num_hypocenters_associated <= hyp_full_assoc_loc_start_index) {
        // remove output files
        int ireturn = nftw(outnameroot, remove_fn, 16, FTW_DEPTH);
        if (ireturn) {
            printf("ERROR: removing files: return value = %d, path = %s\n", ireturn, outnameroot);
        }
        return (0);
    }


    // check for duplicate events, remove event with fewer associated P phases
    // 20110408 AJL - added to avoid processing of same event twice and possible loss of instance with greater associated P phases
    if (1) {
        for (n = 0; n < num_de_data; n++) { // 20140722 AJL - bug fix
            TimedomainProcessingData* deData = data_list[n];
            deData->merged = 0;
        }
        int nhyp_assoc, nhyp_assoc_test;
        int nhyp_preserve, nhyp_remove;
        for (nhyp_assoc = num_hypocenters_associated - 1; nhyp_assoc > 0; nhyp_assoc--) {
            for (nhyp_assoc_test = nhyp_assoc - 1; nhyp_assoc_test >= 0; nhyp_assoc_test--) {
                // check if duplicate
                if (isSameEvent(hyp_assoc_loc[nhyp_assoc], hyp_assoc_loc[nhyp_assoc_test])) {
                    // find hypo with most associated P phases
                    if (hyp_assoc_loc[nhyp_assoc]->nassoc_P > hyp_assoc_loc[nhyp_assoc_test]->nassoc_P) {
                        nhyp_preserve = nhyp_assoc;
                        nhyp_remove = nhyp_assoc_test;
                    } else {
                        nhyp_preserve = nhyp_assoc_test;
                        nhyp_remove = nhyp_assoc;
                    }
                    if (verbose > 0) {
                        printf("   ---> Info: REMOVED DUPLICATE HYPOCENTER!  Associated pick data attached to preserved hypocenter.\n");
                        printHypoDataString(hyp_assoc_loc[nhyp_remove], hypoDataString, 0);
                        printf("   ---> Removed:   %s\n", hypoDataString);
                        printHypoDataString(hyp_assoc_loc[nhyp_preserve], hypoDataString, 0);
                        printf("   ---> Preserved: %s\n", hypoDataString);
                    }
                    // attach data from removed hypo to preserved hypo as zero weight associated data  // 20110411 AJL
                    for (n = 0; n < num_de_data; n++) {
                        TimedomainProcessingData* deData = data_list[n];
                        if (deData->is_associated == nhyp_remove + 1) {
                            deData->is_associated = nhyp_preserve + 1; // attach to preserved hypo
                            hyp_assoc_loc[nhyp_preserve]->nassoc++; // increment number associated for preserved hypo
                            deData->merged = 1;
                            // 20140626 AJL  deData->use_for_location = 0; // 20130218 AJL
                            //deData->is_associated_grid_level = 0;
                            deData->epicentral_distance = GCDistance(hyp_assoc_loc[nhyp_preserve]->lat, hyp_assoc_loc[nhyp_preserve]->lon, deData->lat, deData->lon);
                            deData->epicentral_azimuth = GCAzimuth(hyp_assoc_loc[nhyp_preserve]->lat, hyp_assoc_loc[nhyp_preserve]->lon, deData->lat, deData->lon);
                            deData->residual = (double) deData->t_time_t + deData->t_decsec
                                    - (get_ttime(deData->phase_id, deData->epicentral_distance, hyp_assoc_loc[nhyp_preserve]->depth) + hyp_assoc_loc[nhyp_preserve]->otime);
                            deData->dist_weight = 0.0;
                            deData->polarization.weight = 0.0;
                            deData->loc_weight = 0.0;
                            deData->take_off_angle_inc = get_take_off_angle(deData->phase_id, deData->epicentral_distance, hyp_assoc_loc[nhyp_preserve]->depth);
                            deData->take_off_angle_az = deData->epicentral_azimuth;
                        }
                        if (deData->is_associated > nhyp_remove + 1)
                            deData->is_associated--;
                    }
                    // remove hypocenter from hypo list and shift any remaining hypocenters
                    HypocenterDesc *temp = hyp_assoc_loc[nhyp_remove];
                    for (nhyp = nhyp_remove; nhyp < num_hypocenters_associated - 1; nhyp++) {
                        // 20111109 AJL hyp_assoc_loc[nhyp] = hyp_assoc_loc[nhyp + 1];
                        hyp_assoc_loc[nhyp] = hyp_assoc_loc[nhyp + 1]; // 20111109 AJL - Bug fix        // 20140717 AJL - shift pointers only
                        //(hyp_assoc_loc[nhyp]->hyp_assoc_index)--; // 20140716 AJL - reset later
                    }
                    hyp_assoc_loc[num_hypocenters_associated - 1] = temp;
                    num_hypocenters_associated--;
                    // restart loops
                    nhyp_assoc = num_hypocenters_associated;
                    break;
                }
            }
        }
    }



    // ============================================================================================================
    // do report
    //

    // set internal timing information, time_since_last_report
    // IMPORTANT! - time_since_last_report based on time_max, will only have 1 sec precision.  TODO: make sub-second precision?
    double time_since_last_report;
    if (time_of_last_report > 0.0) {
        time_since_last_report = (double) time_max - time_of_last_report;
    } else {
        time_since_last_report = (double) report_interval; // backward compatibility and support for non real-time processing
    }
    time_of_last_report = (double) time_max;



    // ============================================================================================================
    // create plot in plot time range using 1-sec intervals, write asc file

    double dsec;
    // AJL20100118 - support independent results for each event

    // allocate and initialize arrays and variables ===============================================================

    //T50 level
    int numWarningLevelTotal = 1 + (int) ((T50EX_LEVEL_MAX - T50EX_LEVEL_MIN) / T50EX_LEVEL_STEP);
    char* t50ExArray = calloc(numWarningLevelTotal, sizeof (char));
    int** t50ExHistogram = calloc(num_hypocenters_associated, sizeof (int*));
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        t50ExHistogram[nhyp] = calloc(numWarningLevelTotal, sizeof (int));
        for (m = 0; m < numWarningLevelTotal; m++)
            t50ExHistogram[nhyp][m] = 0;
    }
    double ***t50ExStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        t50ExStatisticsArray[nhyp] = calloc(3, sizeof (double*));
        for (m = 0; m < 3; m++) {
            t50ExStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
        }
    }
    int *numT50ExLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numT50ExLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* t50ExLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char t50ExLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        t50ExLevelStatistics[nhyp].centralValue = 0.0;
        t50ExLevelStatistics[nhyp].upperBound = 0.0;
        t50ExLevelStatistics[nhyp].lowerBound = 0.0;
        t50ExLevelStatistics[nhyp].numLevel = 0;
        strcpy(t50ExLevelString[nhyp], "NONE");
    }
    //tauC
    int numTaucLevelTotal = 1 + (int) ((TAUC_LEVEL_MAX - TAUC_LEVEL_MIN) / TAUC_LEVEL_STEP);
    char* taucArray = calloc(numTaucLevelTotal, sizeof (char));
    int** taucHistogram = calloc(num_hypocenters_associated, sizeof (int*));
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            taucHistogram[nhyp] = calloc(numTaucLevelTotal, sizeof (int));
            for (m = 0; m < numTaucLevelTotal; m++)
                taucHistogram[nhyp][m] = 0;
        }
    }
    double ***taucStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            taucStatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                taucStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
        }
    }
    int *numTaucLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numTaucLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* taucLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char taucLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            taucLevelStatistics[nhyp].centralValue = 0.0;
            taucLevelStatistics[nhyp].upperBound = 0.0;
            taucLevelStatistics[nhyp].lowerBound = 0.0;
            taucLevelStatistics[nhyp].numLevel = 0;
            strcpy(taucLevelString[nhyp], "NONE");
        }
    }
    // TdT50Ex
    int *numAlarmLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numAlarmLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* tdT50ExLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char warningLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            tdT50ExLevelStatistics[nhyp].centralValue = 0.0;
            tdT50ExLevelStatistics[nhyp].upperBound = 0.0;
            tdT50ExLevelStatistics[nhyp].lowerBound = 0.0;
            tdT50ExLevelStatistics[nhyp].numLevel = 0;
            strcpy(warningLevelString[nhyp], "NONE");
        }
    }
    //mwp
    int numMwpLevelTotal = 1 + (int) ((MWP_LEVEL_MAX - MWP_LEVEL_MIN) / MWP_LEVEL_STEP);
#ifdef USE_MWP_LEVEL_ARRAY
    char* mwpArray = calloc(numMwpLevelTotal, sizeof (char));
#endif
    int** mwpHistogram = calloc(num_hypocenters_associated, sizeof (int*));
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpHistogram[nhyp] = calloc(numMwpLevelTotal, sizeof (int));
            for (m = 0; m < numMwpLevelTotal; m++)
                mwpHistogram[nhyp][m] = 0;
        }
    }
    double ***mwpStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpStatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                mwpStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
        }
    }
    int *numMwpLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numMwpLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* mwpLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char mwpLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpLevelStatistics[nhyp].centralValue = 0.0;
            mwpLevelStatistics[nhyp].upperBound = 0.0;
            mwpLevelStatistics[nhyp].lowerBound = 0.0;
            mwpLevelStatistics[nhyp].numLevel = 0;
            strcpy(mwpLevelString[nhyp], "NONE");
        }
    }
    // mb
    int numMbLevelTotal = 1 + (int) ((MB_LEVEL_MAX - MB_LEVEL_MIN) / MB_LEVEL_STEP);
    int** mbHistogram = calloc(num_hypocenters_associated, sizeof (int*));
    if (flag_do_mb) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mbHistogram[nhyp] = calloc(numMbLevelTotal, sizeof (int));
            for (m = 0; m < numMbLevelTotal; m++)
                mbHistogram[nhyp][m] = 0;
        }
    }
    double ***mbStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
    if (flag_do_mb) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mbStatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                mbStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
        }
    }
    int *numMbLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numMbLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* mbLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char mbLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    if (flag_do_mb) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mbLevelStatistics[nhyp].centralValue = 0.0;
            mbLevelStatistics[nhyp].upperBound = 0.0;
            mbLevelStatistics[nhyp].lowerBound = 0.0;
            mbLevelStatistics[nhyp].numLevel = 0;
            strcpy(mbLevelString[nhyp], "NONE");
        }
    }
    // t0
    int numT0LevelTotal = 1 + (int) ((T0_LEVEL_MAX - T0_LEVEL_MIN) / T0_LEVEL_STEP);
    int** t0Histogram = calloc(num_hypocenters_associated, sizeof (int*));
    if (flag_do_t0) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            t0Histogram[nhyp] = calloc(numT0LevelTotal, sizeof (int));
            for (m = 0; m < numT0LevelTotal; m++)
                t0Histogram[nhyp][m] = 0;
        }
    }
    double ***t0StatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
    if (flag_do_t0) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            t0StatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                t0StatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
        }
    }
    int *numT0Level = calloc(num_hypocenters_associated, sizeof (int));
    int *numT0LevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* t0LevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char t0LevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
    if (flag_do_t0) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            t0LevelStatistics[nhyp].centralValue = 0.0;
            t0LevelStatistics[nhyp].upperBound = 0.0;
            t0LevelStatistics[nhyp].lowerBound = 0.0;
            t0LevelStatistics[nhyp].numLevel = 0;
            strcpy(t0LevelString[nhyp], "NONE");
        }
    }
    // mwpd
    int numMwpdLevelTotal = 1 + (int) ((MWPD_LEVEL_MAX - MWPD_LEVEL_MIN) / MWPD_LEVEL_STEP);
    int** mwpdHistogram = calloc(num_hypocenters_associated, sizeof (int*));
    if (flag_do_mwpd) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpdHistogram[nhyp] = calloc(numMwpdLevelTotal, sizeof (int));
            for (m = 0; m < numMwpdLevelTotal; m++)
                mwpdHistogram[nhyp][m] = 0;
        }
    }
    double ***mwpdStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
#ifdef USE_MWP_MO_POS_NEG
    double ***mwpdMoPosNegStatisticsArray = calloc(num_hypocenters_associated, sizeof (double**));
#endif
    if (flag_do_mwpd) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpdStatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                mwpdStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
#ifdef USE_MWP_MO_POS_NEG
            mwpdMoPosNegStatisticsArray[nhyp] = calloc(3, sizeof (double*));
            for (m = 0; m < 3; m++) {
                mwpdMoPosNegStatisticsArray[nhyp][m] = calloc(num_de_data, sizeof (double));
            }
#endif
        }
    }
    int *numMwpdLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numMwpdLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* mwpdLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
    char mwpdLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];
#ifdef USE_MWP_MO_POS_NEG
    int *numMwpdMoPosNegLevel = calloc(num_hypocenters_associated, sizeof (int));
    int *numMwpdMoPosNegLevelMax = calloc(num_hypocenters_associated, sizeof (int));
    statistic_level* mwpdMoPosNegLevelStatistics = calloc(num_hypocenters_associated, sizeof (statistic_level));
#endif
    if (flag_do_mwpd) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpdLevelStatistics[nhyp].centralValue = 0.0;
            mwpdLevelStatistics[nhyp].upperBound = 0.0;
            mwpdLevelStatistics[nhyp].lowerBound = 0.0;
            mwpdLevelStatistics[nhyp].numLevel = 0;
            strcpy(mwpdLevelString[nhyp], "NONE");
        }
#ifdef USE_MWP_MO_POS_NEG
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            mwpdMoPosNegLevelStatistics[nhyp].centralValue = 0.0;
            mwpdMoPosNegLevelStatistics[nhyp].upperBound = 0.0;
            mwpdMoPosNegLevelStatistics[nhyp].lowerBound = 0.0;
            mwpdMoPosNegLevelStatistics[nhyp].numLevel = 0;
        }
#endif
    }



    // ============================================================================================================
    // open output files
    //
    // agency
    //
    sprintf(outname, "%s/agency_id.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * agencyStream = fopen_counter(outname, "w");
    if (agencyStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(agencyStream, "%s\n", agencyId);
    fclose_counter(agencyStream);
    // program version
    //
    sprintf(outname, "%s/version.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * versionStream = fopen_counter(outname, "w");
    if (versionStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(versionStream, "v%s (%s)\n", EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE);
    fclose_counter(versionStream);
    // temp directory
    sprintf(outname, "%s", EE_TEMP_DIR);
    mkdir(outname, 0755);
    //
    // event directory
    sprintf(outname, "%s/events", outnameroot);
    mkdir(outname, 0755);
    //
    // plotting files
    sprintf(outname, "%s/plot", outnameroot);
    mkdir(outname, 0755);
    // GMT T50 grd data
    sprintf(outname, "%s/plot/t50.grid.bin", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * t50ExGridDataStream = fopen_counter(outname, "w");
    if (t50ExGridDataStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.grid.sta.code.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * t50ExStaCodeGridStream = fopen_counter(outname, "w");
    if (t50ExStaCodeGridStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    //
    // tauc
    FILE* taucGridDataStrem = NULL;
    FILE* taucStaCodeGridStream = NULL;
    if (flag_do_tauc) {
        // GMT tauc grd data
        sprintf(outname, "%s/plot/tauc.grid.bin", outnameroot);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        taucGridDataStrem = fopen_counter(outname, "w");
        if (taucGridDataStrem == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/tauc.grid.sta.code.txt", outnameroot);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        taucStaCodeGridStream = fopen_counter(outname, "w");
        if (taucStaCodeGridStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
    }
    //
    // mwp
#ifdef USE_MWP_LEVEL_ARRAY
    FILE* mwpGridDataStrem = NULL;
    FILE* mwpStaCodeGridStream = NULL;
    if (flag_do_mwp) {
        // GMT mwp grd data
        sprintf(outname, "%s/plot/mwp.grid.bin", outnameroot);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        mwpGridDataStrem = fopen_counter(outname, "w");
        if (mwpGridDataStrem == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/mwp.grid.sta.code.txt", outnameroot);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        mwpStaCodeGridStream = fopen_counter(outname, "w");
        if (mwpStaCodeGridStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
    }
#endif
    //
    //
    sprintf(outname, "%s/plot/pick.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * pickStream = fopen_counter(outname, "w");
    if (pickStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }

    FILE * pickStreamAssoc[num_hypocenters_associated];
    FILE * staAssociatedStream[num_hypocenters_associated];
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        //
        sprintf(outname, "%s/plot/pick.assoc.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        pickStreamAssoc[nhyp] = fopen_counter(outname, "w");
        if (pickStreamAssoc[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }

        // GMT style
        sprintf(outname, "%s/plot/sta.associated.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        staAssociatedStream[nhyp] = fopen_counter(outname, "w");
        if (staAssociatedStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
    }

    FILE * t50ExCentralValueStream[num_hypocenters_associated];
    FILE * t50ExUpperBoundStream[num_hypocenters_associated];
    FILE * t50ExLowerBoundStream[num_hypocenters_associated];
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        // GMT style
        sprintf(outname, "%s/plot/t50.value.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        t50ExCentralValueStream[nhyp] = fopen_counter(outname, "w");
        if (t50ExCentralValueStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        fprintf(t50ExCentralValueStream[nhyp], ">\n");
        // GMT style
        sprintf(outname, "%s/plot/t50.upper.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        t50ExUpperBoundStream[nhyp] = fopen_counter(outname, "w");
        if (t50ExUpperBoundStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        fprintf(t50ExUpperBoundStream[nhyp], ">\n");
        // GMT style
        sprintf(outname, "%s/plot/t50.lower.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        t50ExLowerBoundStream[nhyp] = fopen_counter(outname, "w");
        if (t50ExLowerBoundStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        fprintf(t50ExLowerBoundStream[nhyp], ">\n");
    }
    FILE * taucCentralValueStream[num_hypocenters_associated];
    FILE * taucUpperBoundStream[num_hypocenters_associated];
    FILE * taucLowerBoundStream[num_hypocenters_associated];
    FILE * tdT50ExCentralValueStream[num_hypocenters_associated];
    FILE * tdT50ExUpperBoundStream[num_hypocenters_associated];
    FILE * tdT50ExLowerBoundStream[num_hypocenters_associated];
#ifdef USE_MWP_LEVEL_ARRAY
    FILE * mwpCentralValueStream[num_hypocenters_associated];
    FILE * mwpUpperBoundStream[num_hypocenters_associated];
    FILE * mwpLowerBoundStream[num_hypocenters_associated];
#endif
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            // GMT style
            sprintf(outname, "%s/plot/tauc.value.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            taucCentralValueStream[nhyp] = fopen_counter(outname, "w");
            if (taucCentralValueStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(taucCentralValueStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/tauc.upper.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            taucUpperBoundStream[nhyp] = fopen_counter(outname, "w");
            if (taucUpperBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(taucUpperBoundStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/tauc.lower.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            taucLowerBoundStream[nhyp] = fopen_counter(outname, "w");
            if (taucLowerBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(taucLowerBoundStream[nhyp], ">\n");
        }
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            // GMT style
            sprintf(outname, "%s/plot/alarm.value.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            tdT50ExCentralValueStream[nhyp] = fopen_counter(outname, "w");
            if (tdT50ExCentralValueStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(tdT50ExCentralValueStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/alarm.upper.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            tdT50ExUpperBoundStream[nhyp] = fopen_counter(outname, "w");
            if (tdT50ExUpperBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(tdT50ExUpperBoundStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/alarm.lower.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            tdT50ExLowerBoundStream[nhyp] = fopen_counter(outname, "w");
            if (tdT50ExLowerBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(tdT50ExLowerBoundStream[nhyp], ">\n");
        }
    }
#ifdef USE_MWP_LEVEL_ARRAY
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            // GMT style
            sprintf(outname, "%s/plot/mwp.value.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mwpCentralValueStream[nhyp] = fopen_counter(outname, "w");
            if (mwpCentralValueStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(mwpCentralValueStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/mwp.upper.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mwpUpperBoundStream[nhyp] = fopen_counter(outname, "w");
            if (mwpUpperBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(mwpUpperBoundStream[nhyp], ">\n");
            // GMT style
            sprintf(outname, "%s/plot/mwp.lower.%d.xy", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mwpLowerBoundStream[nhyp] = fopen_counter(outname, "w");
            if (mwpLowerBoundStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            fprintf(mwpLowerBoundStream[nhyp], ">\n");
        }
    }
#endif
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.red.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50RedStream = fopen_counter(outname, "w");
    if (staT50RedStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.yellow.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50YellowStream = fopen_counter(outname, "w");
    if (staT50YellowStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.green.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50GreenStream = fopen_counter(outname, "w");
    if (staT50GreenStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.ltred.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50LtRedStream = fopen_counter(outname, "w");
    if (staT50LtRedStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.ltyellow.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50LtYellowStream = fopen_counter(outname, "w");
    if (staT50LtYellowStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/t50.sta.ltgreen.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staT50LtGreenStream = fopen_counter(outname, "w");
    if (staT50LtGreenStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }

    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.red.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucRedStream = fopen_counter(outname, "w");
    if (staTaucRedStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.yellow.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucYellowStream = fopen_counter(outname, "w");
    if (staTaucYellowStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.green.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucGreenStream = fopen_counter(outname, "w");
    if (staTaucGreenStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.ltred.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucLtRedStream = fopen_counter(outname, "w");
    if (staTaucLtRedStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.ltyellow.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucLtYellowStream = fopen_counter(outname, "w");
    if (staTaucLtYellowStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/tauc.sta.ltgreen.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staTaucLtGreenStream = fopen_counter(outname, "w");
    if (staTaucLtGreenStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }


    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        numT50ExLevelMax[nhyp] = 0;
    }
    //tauc
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numTaucLevelMax[nhyp] = 0;
        }
    }
    //mwp
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numMwpLevelMax[nhyp] = 0;
        }
    }
    //mb
    if (flag_do_mb) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numMbLevelMax[nhyp] = 0;
        }
    }
    //t0
    if (flag_do_t0) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numT0LevelMax[nhyp] = 0;
        }
    }
    //mwpd
    if (flag_do_mwpd) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numMwpdLevelMax[nhyp] = 0;
        }
#ifdef USE_MWP_MO_POS_NEG
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numMwpdMoPosNegLevelMax[nhyp] = 0;
        }
#endif
    }


    // loop over all time steps in plot time, loop over all data, load level arrays, accumulate statistics ===============================================================

    time_t curr_time = plot_time_min;
    while (curr_time <= plot_time_max) {
        //T50 Level
        for (nhyp = 0; nhyp < numWarningLevelTotal; nhyp++)
            t50ExArray[nhyp] = 0;
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            numT50ExLevel[nhyp] = 0;
        }
        //tauc
        if (flag_do_tauc) {
            for (nhyp = 0; nhyp < numTaucLevelTotal; nhyp++)
                taucArray[nhyp] = 0;
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numTaucLevel[nhyp] = 0;
            }
        }
        //mwp
        if (flag_do_mwp) {
#ifdef USE_MWP_LEVEL_ARRAY
            for (nhyp = 0; nhyp < numMwpLevelTotal; nhyp++)
                mwpArray[nhyp] = 0;
#endif
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numMwpLevel[nhyp] = 0;
            }
        }
        //mb
        if (flag_do_mb) {
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numMbLevel[nhyp] = 0;
            }
        }
        //t0
        if (flag_do_t0) {
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numT0Level[nhyp] = 0;
            }
        }
        //mwpd
        if (flag_do_mwpd) {
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numMwpdLevel[nhyp] = 0;
            }
#ifdef USE_MWP_MO_POS_NEG
            for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
                numMwpdMoPosNegLevel[nhyp] = 0;
            }
#endif
        }
        // loop over data
        int ndata;
        for (ndata = 0; ndata < num_de_data; ndata++) {
            TimedomainProcessingData* deData = data_list[ndata];
            // skip if possible clip in data time span or if HF or BRB s/n ratio too low
            // AJL 20100224
            //if (deData->flag_clipped || deData->flag_snr_hf_too_low || deData->flag_a_ref_not_ok)
            if (deData->flag_clipped || deData->flag_non_contiguous || (!USE_AREF_NOT_OK_PICKS_FOR_LOCATION && deData->flag_a_ref_not_ok))
                continue;
            // skip if data not used for location and not associated
            // AJL 20100630
            if (!deData->use_for_location && !(use_associated_data && is_associated_phase(deData)))
                continue;
            // if time step includes data pick time, update calculated pick values
            time_t pick_time = get_time_t(deData, &dsec);
            if (pick_time >= curr_time && pick_time < curr_time + TIME_STEP) {
                double dtime = difftime(curr_time, plot_time_max) / 60.0;
                if (flag_do_mwpd) {
                    deData->mwp->mag = MWP_INVALID;
                }
                if (flag_do_mb) {
                    deData->mb->mag = MB_INVALID;
                }
                double depth = 0.0;
                if (deData->is_associated)
                    depth = hyp_assoc_loc[deData->is_associated - 1]->depth;
                set_derived_values(deData, depth); // 20111222 TEST AJL - use S duration
                if (use_associated_data && is_associated_phase(deData)) {
                    fprintf(pickStreamAssoc[deData->is_associated - 1], ">\n%f %f\n%f %f\n", dtime, PICK_PLOT_LEVEL_MIN, dtime, PICK_PLOT_LEVEL_MAX);
                    fprintf(staAssociatedStream[deData->is_associated - 1], "%f %f\n", deData->lat, deData->lon);
                    //20110119 AJL if (is_associated_location_P(deData) && flag_do_mwp && deData->flag_complete_mwp) {
                    // 20140808 AJL - exclude Pdiff  if (is_associated_location_P(deData) && flag_do_mwp) {
                    if (is_associated_location_P(deData) && flag_do_mwp && is_P(deData->phase_id)) {
                        calculate_Mwp_Mag(deData, hyp_assoc_loc[deData->is_associated - 1]->depth,
                                use_mwp_distance_correction, use_amplitude_attenuation && use_magnitude_amp_atten_check);
#ifdef USE_MWP_LEVEL_ARRAY
                        // 20121119 AJL if (!deData->flag_snr_brb_too_low && is_unassociated_or_location_P(deData))
                        if (!deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low && is_unassociated_or_location_P(deData))
                            fprintf(mwpStaCodeGridStream, "%f %f %s\n", difftime(pick_time, plot_time_max) / 60.0, deData->mwp->mag, deData->station);
#endif
                    }
                    //20110119 AJL if (is_associated_location_P(deData) && flag_do_mb && deData->flag_complete_mb) {
                    // 20140808 AJL - exclude Pdiff  if (is_associated_location_P(deData) && flag_do_mb) {
                    if (is_associated_location_P(deData) && flag_do_mb && is_P(deData->phase_id)) {
                        if (MB_MODE == MB_MODE_mB)
                            calculate_mB_Mag(deData, hyp_assoc_loc[deData->is_associated - 1]->depth,
                                use_mb_correction, use_amplitude_attenuation && use_magnitude_amp_atten_check);
                        else
                            calculate_mb_Mag(deData, hyp_assoc_loc[deData->is_associated - 1]->depth,
                                use_amplitude_attenuation && use_magnitude_amp_atten_check); // 20110528 AJL
                        //printf("DEBUG: %s S/N-BRB-BP: sn_brb_bp_pick %f, sn_brb_bp_signal %f, snr_brb_bp %f\n",
                        //        deData->station, deData->sn_brb_bp_pick, deData->sn_brb_bp_signal,
                        //        deData->sn_brb_bp_pick < FLT_MIN ? -1.0 : deData->sn_brb_bp_signal / deData->sn_brb_bp_pick);
                    }
                    // 20140808 AJL - exclude Pdiff  if (is_associated_location_P(deData) && flag_do_mwpd) {
                    if (is_associated_location_P(deData) && flag_do_mwpd && is_P(deData->phase_id)) {
                        calculate_Raw_Mwpd_Mag(deData, hyp_assoc_loc[deData->is_associated - 1]->depth,
                                use_mwp_distance_correction, use_amplitude_attenuation && use_magnitude_amp_atten_check);
                    }
                }
                if (!deData->flag_snr_hf_too_low && is_unassociated_or_location_P(deData)) {
                    double wlevel = getT50Level(deData);
                    //wlevel = T50EX_LEVEL_MIN + T50EX_LEVEL_STEP * (int) ((wlevel - T50EX_LEVEL_MIN) / T50EX_LEVEL_STEP);
                    fprintf(t50ExStaCodeGridStream, "%f %f %s\n", difftime(pick_time, plot_time_max) / 60.0, wlevel, deData->station);
                }
                if (flag_do_tauc) {
                    if (!deData->flag_snr_brb_too_low && is_unassociated_or_location_P(deData))
                        fprintf(taucStaCodeGridStream, "%f %f %s\n", difftime(pick_time, plot_time_max) / 60.0, deData->tauc_peak, deData->station);
                }
                // 20130128 AJL - use flag_snr_brb_int_too_low to allow mwp, mwpd, etc., but do not use for ignore tests (e.g. ignore determined by flag_snr_brb_too_low)
                //if (!deData->flag_snr_hf_too_low || !deData->flag_snr_brb_too_low || !deData->flag_snr_brb_int_too_low)
                if (!deData->flag_snr_hf_too_low || !deData->flag_snr_brb_too_low)
                    fprintf(pickStream, ">\n%f %f\n%f %f\n", dtime, PICK_PLOT_LEVEL_MIN, dtime, PICK_PLOT_LEVEL_MAX);
            }
            // skip if not completed
            if (!deData->flag_complete_t50)
                continue;
            // include data discriminant levels in grid and level statistics
            time_t plot_time = pick_time; // plot at pick time
            // set plot length to specified window after OT or default length if this data has no associated hypocenter
            double plot_time_end = plot_time + LEVEL_PLOT_WINDOW_LENGTH_DEFAULT;
            if (use_associated_data && is_associated_location_P(deData))
                plot_time_end = (time_t) hyp_assoc_loc[deData->is_associated - 1]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED;
            if (curr_time >= plot_time && curr_time <= plot_time_end) {
                //T50 Level
                double t50Level = getT50Level(deData);
                int indexT50 = 0;
                int associated_ok = use_associated_data && is_associated_location_P(deData)
                        && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_WARNING && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_WARNING;
                if (!deData->flag_snr_hf_too_low
                        && (!associate_data || !deData->is_associated || associated_ok)) {
                    if (t50Level >= T50EX_LEVEL_MAX - T50EX_LEVEL_STEP)
                        indexT50 = numWarningLevelTotal - 1;
                    else if (t50Level <= T50EX_LEVEL_MIN)
                        indexT50 = 0;
                    else
                        indexT50 = (int) ((t50Level - T50EX_LEVEL_MIN) / T50EX_LEVEL_STEP);
                    if (t50Level >= T50EX_RED_CUTOFF) {
                        if (t50ExArray[indexT50] < 1)
                            t50ExArray[indexT50] += T50EX_MODULO;
                        if (t50ExArray[indexT50] < 127)
                            t50ExArray[indexT50]++;
                    } else if (t50Level >= T50EX_YELLOW_CUTOFF) {
                        if (t50ExArray[indexT50] < T50EX_MODULO)
                            t50ExArray[indexT50]++;
                    } else {
                        //printf("DEBUG: numWarningLevelTotal %d  indexT50 %d  t50Level %f  T50EX_LEVEL_MIN %f  T50EX_LEVEL_STEP %f\n",
                        //        numWarningLevelTotal, indexT50, t50Level, T50EX_LEVEL_MIN, T50EX_LEVEL_STEP);
                        if (t50ExArray[indexT50] > -128)
                            t50ExArray[indexT50]--;
                    }
                    // add to T50 level statistics array - check associated and distance cutoff
                    if (!associate_data || associated_ok) {
                        t50ExStatisticsArray[(deData->is_associated - 1)][0][numT50ExLevel[(deData->is_associated - 1)]] = t50Level;
                        t50ExStatisticsArray[(deData->is_associated - 1)][1][numT50ExLevel[(deData->is_associated - 1)]] = deData->lat;
                        t50ExStatisticsArray[(deData->is_associated - 1)][2][numT50ExLevel[(deData->is_associated - 1)]] = deData->lon;
                        numT50ExLevel[(deData->is_associated - 1)]++;
                        if (numT50ExLevel[(deData->is_associated - 1)] > numT50ExLevelMax[(deData->is_associated - 1)])
                            numT50ExLevelMax[(deData->is_associated - 1)] = numT50ExLevel[(deData->is_associated - 1)];
                    }
                }
                //tauc
                int indexTauc = 0;
                associated_ok = use_associated_data && is_associated_location_P(deData)
                        && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_TAUC && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_TAUC;
                if (flag_do_tauc && deData->tauc_peak != TAUC_INVALID && !deData->flag_snr_brb_too_low
                        && (!associate_data || !deData->is_associated || associated_ok)) {
                    double taucLevel = deData->tauc_peak;
                    if (taucLevel >= TAUC_LEVEL_MAX - TAUC_LEVEL_STEP)
                        indexTauc = numTaucLevelTotal - 1;
                    else if (taucLevel <= TAUC_LEVEL_MIN)
                        indexTauc = 0;
                    else
                        indexTauc = (int) ((taucLevel - TAUC_LEVEL_MIN) / TAUC_LEVEL_STEP);
                    if (taucLevel >= TAUC_RED_CUTOFF) {
                        if (taucArray[indexTauc] < 1)
                            taucArray[indexTauc] += T50EX_MODULO;
                        if (taucArray[indexTauc] < 127)
                            taucArray[indexTauc]++;
                    } else if (taucLevel >= TAUC_YELLOW_CUTOFF) {
                        if (taucArray[indexTauc] < T50EX_MODULO)
                            taucArray[indexTauc]++;
                    } else {
                        if (taucArray[indexTauc] > -128)
                            taucArray[indexTauc]--;
                    }
                    // add tauc to tauc level statistics array - check associated and distance cutoff
                    if (!associate_data || associated_ok) {
                        taucStatisticsArray[(deData->is_associated - 1)][0][numTaucLevel[(deData->is_associated - 1)]] = taucLevel;
                        taucStatisticsArray[(deData->is_associated - 1)][1][numTaucLevel[(deData->is_associated - 1)]] = deData->lat;
                        taucStatisticsArray[(deData->is_associated - 1)][2][numTaucLevel[(deData->is_associated - 1)]] = deData->lon;
                        numTaucLevel[(deData->is_associated - 1)]++;
                        if (numTaucLevel[(deData->is_associated - 1)] > numTaucLevelMax[(deData->is_associated - 1)])
                            numTaucLevelMax[(deData->is_associated - 1)] = numTaucLevel[(deData->is_associated - 1)];
                    }
                }
                //mwp
                int indexMwp = 0;
                associated_ok = use_associated_data && is_associated_location_P(deData)
                        && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MWP && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MWP;
                // 20121119 AJL if (flag_do_mwp && deData->mwp->mag != MWP_INVALID && !deData->flag_snr_brb_too_low
                if (flag_do_mwp && deData->mwp->mag != MWP_INVALID && !deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low
                        && (!associate_data || !deData->is_associated || associated_ok)) {
                    double mwp_mag = deData->mwp->mag;
                    if (mwp_mag >= MWP_LEVEL_MAX - MWP_LEVEL_STEP)
                        indexMwp = numMwpLevelTotal - 1;
                    else if (mwp_mag <= MWP_LEVEL_MIN)
                        indexMwp = 0;
                    else
                        indexMwp = (int) ((mwp_mag - MWP_LEVEL_MIN) / MWP_LEVEL_STEP);
#ifdef USE_MWP_LEVEL_ARRAY
                    if (mwp_mag >= MWP_RED_CUTOFF) {
                        if (mwpArray[indexMwp] < 1)
                            mwpArray[indexMwp] += T50EX_MODULO;
                        if (mwpArray[indexMwp] < 127)
                            mwpArray[indexMwp]++;
                    } else if (mwp_mag >= MWP_YELLOW_CUTOFF) {
                        if (mwpArray[indexMwp] < T50EX_MODULO)
                            mwpArray[indexMwp]++;
                    } else {
                        if (mwpArray[indexMwp] > -128)
                            mwpArray[indexMwp]--;
                    }
#endif
                    // add mwp to mwp level statistics array - check associated and distance cutoff
                    if (!associate_data || associated_ok) {
                        mwpStatisticsArray[(deData->is_associated - 1)][0][numMwpLevel[(deData->is_associated - 1)]] = mwp_mag;
                        mwpStatisticsArray[(deData->is_associated - 1)][1][numMwpLevel[(deData->is_associated - 1)]] = deData->lat;
                        mwpStatisticsArray[(deData->is_associated - 1)][2][numMwpLevel[(deData->is_associated - 1)]] = deData->lon;
                        numMwpLevel[(deData->is_associated - 1)]++;
                        if (numMwpLevel[(deData->is_associated - 1)] > numMwpLevelMax[(deData->is_associated - 1)])
                            numMwpLevelMax[(deData->is_associated - 1)] = numMwpLevel[(deData->is_associated - 1)];
                    }
                }
                // process parameters without results array and histograms and station colors
                if (curr_time + TIME_STEP > plot_time_end || curr_time + TIME_STEP > plot_time_max) { // last time step for this deData
                    // parameters without results array
                    //mb
                    int indexMb = 0;
                    associated_ok = use_associated_data && is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MB && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MB;
                    if (flag_do_mb && deData->mb->mag != MB_INVALID && !deData->flag_snr_brb_bp_too_low
                            && (!associate_data || !deData->is_associated || associated_ok)) {
                        double mb_mag = deData->mb->mag;
                        if (mb_mag >= MB_LEVEL_MAX - MB_LEVEL_STEP)
                            indexMb = numMbLevelTotal - 1;
                        else if (mb_mag <= MB_LEVEL_MIN)
                            indexMb = 0;
                        else
                            indexMb = (int) ((mb_mag - MB_LEVEL_MIN) / MB_LEVEL_STEP);
                        // add mb to mb level statistics array - check associated and distance cutoff
                        if (!associate_data || associated_ok) {
                            mbStatisticsArray[(deData->is_associated - 1)][0][numMbLevel[(deData->is_associated - 1)]] = mb_mag;
                            mbStatisticsArray[(deData->is_associated - 1)][1][numMbLevel[(deData->is_associated - 1)]] = deData->lat;
                            mbStatisticsArray[(deData->is_associated - 1)][2][numMbLevel[(deData->is_associated - 1)]] = deData->lon;
                            numMbLevel[(deData->is_associated - 1)]++;
                            if (numMbLevel[(deData->is_associated - 1)] > numMbLevelMax[(deData->is_associated - 1)])
                                numMbLevelMax[(deData->is_associated - 1)] = numMbLevel[(deData->is_associated - 1)];
                        }
                    }
                    //t0
                    int indexT0 = 0;
                    associated_ok = use_associated_data && is_associated_location_P(deData)
                            && useT0Report(deData);
                    if (flag_do_t0 && deData->t0->duration_raw != T0_INVALID && !deData->flag_snr_hf_too_low
                            && (!associate_data || !deData->is_associated || associated_ok)) {
                        //double t0_dur = deData->t0->duration_raw;
                        //double depth = 0.0;
                        //if (deData->is_associated)
                        //    depth = hyp_assoc_loc[deData->is_associated - 1]->depth;
                        //double t0_dur = calculate_corrected_duration(deData, depth); // 20111222 TEST AJL - use S duration
                        double t0_dur = deData->t0->duration_plot;
                        //if (t0_dur < 0.0)
                        //    printf("ERROR: deData->t0->duration_plot not set: this should not happen!\n");
                        if (t0_dur >= T0_LEVEL_MAX - T0_LEVEL_STEP)
                            indexT0 = numT0LevelTotal - 1;
                        else if (t0_dur <= T0_LEVEL_MIN)
                            indexT0 = 0;
                        else
                            indexT0 = (int) ((t0_dur - T0_LEVEL_MIN) / T0_LEVEL_STEP);
                        // add t0 to t0 level statistics array - check associated and distance cutoff
                        if (!associate_data || associated_ok) {
                            t0StatisticsArray[(deData->is_associated - 1)][0][numT0Level[(deData->is_associated - 1)]] = t0_dur;
                            t0StatisticsArray[(deData->is_associated - 1)][1][numT0Level[(deData->is_associated - 1)]] = deData->lat;
                            t0StatisticsArray[(deData->is_associated - 1)][2][numT0Level[(deData->is_associated - 1)]] = deData->lon;
                            numT0Level[(deData->is_associated - 1)]++;
                            if (numT0Level[(deData->is_associated - 1)] > numT0LevelMax[(deData->is_associated - 1)])
                                numT0LevelMax[(deData->is_associated - 1)] = numT0Level[(deData->is_associated - 1)];
                        }
                    }
                    // histograms and station colors
                    // tauc
                    if (flag_do_tauc && deData->tauc_peak != TAUC_INVALID && !deData->flag_snr_brb_too_low) {
                        int flag_active_tauc = !associate_data || (use_associated_data && is_associated_location_P(deData)
                                && (deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_TAUC && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_TAUC));
                        if (flag_active_tauc) {
                            (taucHistogram[(deData->is_associated - 1)][indexTauc])++;
                        }
                        int not_associated = !deData->is_associated || deData->is_associated == NUMBER_ASSOCIATE_IGNORE;
                        if (curr_time + TIME_STEP > plot_time_max) { // last time step for this deData overlaps plot time max
                            double taucLevel = deData->tauc_peak;
                            if (taucLevel >= TAUC_RED_CUTOFF) {
                                if (flag_active_tauc)
                                    fprintf(staTaucRedStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staTaucLtRedStream, "%f %f\n", deData->lat, deData->lon);
                            } else if (taucLevel >= TAUC_YELLOW_CUTOFF) {
                                if (flag_active_tauc)
                                    fprintf(staTaucYellowStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staTaucLtYellowStream, "%f %f\n", deData->lat, deData->lon);
                            } else {
                                if (flag_active_tauc)
                                    fprintf(staTaucGreenStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staTaucLtGreenStream, "%f %f\n", deData->lat, deData->lon);
                            }
                        }
                    }
                    // 20121119 AJL if (flag_do_mwp && deData->mwp->mag != MWP_INVALID && !deData->flag_snr_brb_too_low) {
                    if (flag_do_mwp && deData->mwp->mag != MWP_INVALID && !deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low) {
                        int flag_active_mwp = !associate_data || (use_associated_data && is_associated_location_P(deData)
                                && (deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MWP && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MWP));
                        if (flag_active_mwp) {
                            (mwpHistogram[(deData->is_associated - 1)][indexMwp])++;
                        }
                    }
                    if (flag_do_mb && deData->mb->mag != MB_INVALID && !deData->flag_snr_brb_bp_too_low) {
                        int flag_active_mb = !associate_data || (use_associated_data && is_associated_location_P(deData)
                                && (deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MB && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MB));
                        if (flag_active_mb) {
                            (mbHistogram[(deData->is_associated - 1)][indexMb])++;
                        }
                    }
                    if (flag_do_t0 && deData->t0->duration_raw != T0_INVALID && !deData->flag_snr_hf_too_low) {
                        int flag_active_t0 = !associate_data || (use_associated_data && is_associated_location_P(deData)
                                && useT0Report(deData));
                        if (flag_active_t0) {
                            (t0Histogram[(deData->is_associated - 1)][indexT0])++;
                        }
                    }
                    if (!deData->flag_snr_hf_too_low) {
                        int flag_active_t50Ex = !associate_data || (use_associated_data && is_associated_location_P(deData)
                                && (deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_WARNING && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_WARNING));
                        if (flag_active_t50Ex) {
                            (t50ExHistogram[(deData->is_associated - 1)][indexT50])++;
                        }
                        int not_associated = !deData->is_associated || deData->is_associated == NUMBER_ASSOCIATE_IGNORE;
                        if (curr_time + TIME_STEP > plot_time_max) { // last time step for this deData overlaps plot time max
                            if (t50Level >= T50EX_RED_CUTOFF) {
                                if (flag_active_t50Ex)
                                    fprintf(staT50RedStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staT50LtRedStream, "%f %f\n", deData->lat, deData->lon);
                            } else if (t50Level >= T50EX_YELLOW_CUTOFF) {
                                if (flag_active_t50Ex)
                                    fprintf(staT50YellowStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staT50LtYellowStream, "%f %f\n", deData->lat, deData->lon);
                            } else {
                                if (flag_active_t50Ex)
                                    fprintf(staT50GreenStream, "%f %f\n", deData->lat, deData->lon);
                                else if (not_associated)
                                    fprintf(staT50LtGreenStream, "%f %f\n", deData->lat, deData->lon);
                            }
                        }
                    }

                }
            }
        }
        fwrite(t50ExArray, sizeof (char), numWarningLevelTotal, t50ExGridDataStream);
        if (flag_do_tauc) {
            fwrite(taucArray, sizeof (char), numTaucLevelTotal, taucGridDataStrem);
        }
#ifdef USE_MWP_LEVEL_ARRAY
        if (flag_do_mwp) {
            fwrite(mwpArray, sizeof (char), numMwpLevelTotal, mwpGridDataStrem);
        }
#endif


        // loop over associated hypocenters, develop statistics ===============================================================

        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {


            //T50 Level
            if (numT50ExLevel[nhyp] > 0) {
                // evaluate and set level statistics
                setStatistics("T50Ex", t50ExStatisticsArray[nhyp], numT50ExLevel[nhyp], &(t50ExLevelStatistics[nhyp]),
                        0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                        && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
            }
            //else {
            //    t50ExLevelStatistics[n].centralValue = t50ExLevelStatistics[n].upperBound = t50ExLevelStatistics[n].lowerBound = 0.0;
            //}
            double t50ExCentralValueWrite = t50ExLevelStatistics[nhyp].centralValue > T50EX_LEVEL_MAX_PLOT ? T50EX_LEVEL_MAX_PLOT : t50ExLevelStatistics[nhyp].centralValue;
            t50ExCentralValueWrite = t50ExCentralValueWrite < T50EX_LEVEL_MIN ? T50EX_LEVEL_MIN : t50ExCentralValueWrite;
            //
            double t50ExUpperBoundWrite = t50ExLevelStatistics[nhyp].upperBound > T50EX_LEVEL_MAX_PLOT ? T50EX_LEVEL_MAX_PLOT : t50ExLevelStatistics[nhyp].upperBound;
            t50ExUpperBoundWrite = t50ExUpperBoundWrite < T50EX_LEVEL_MIN ? T50EX_LEVEL_MIN : t50ExUpperBoundWrite;
            //
            double t50ExLowerBoundWrite = t50ExLevelStatistics[nhyp].lowerBound > T50EX_LEVEL_MAX_PLOT ? T50EX_LEVEL_MAX_PLOT : t50ExLevelStatistics[nhyp].lowerBound;
            t50ExLowerBoundWrite = t50ExLowerBoundWrite < T50EX_LEVEL_MIN ? T50EX_LEVEL_MIN : t50ExLowerBoundWrite;
            // check for enough T50 levels to plot statistics
            if (numT50ExLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                fprintf(t50ExCentralValueStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, t50ExCentralValueWrite);
                fprintf(t50ExUpperBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, t50ExUpperBoundWrite);
                fprintf(t50ExLowerBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, t50ExLowerBoundWrite);
            } else {
                fprintf(t50ExCentralValueStream[nhyp], ">\n");
                fprintf(t50ExUpperBoundStream[nhyp], ">\n");
                fprintf(t50ExLowerBoundStream[nhyp], ">\n");
            }
            setLevelString(numT50ExLevelMax[nhyp], &t50ExLevelStatistics[nhyp], t50ExLevelString[nhyp],
                    MIN_NUMBER_VALUES_USE, T50EX_RED_CUTOFF, T50EX_YELLOW_CUTOFF, -1.0, warning_colors_show);
            //
            //tauc
            if (flag_do_tauc) {
                if (numTaucLevel[nhyp] > 0) {
                    // evaluate and set level statistics
                    setStatistics("Tauc", taucStatisticsArray[nhyp], numTaucLevel[nhyp], &(taucLevelStatistics[nhyp]),
                            0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                            && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
                }
                //else {
                //    taucLevelStatistics[n].centralValue = taucLevelStatistics[n].upperBound = taucLevelStatistics[n].lowerBound = 0.0;
                //}
                double taucCentralValueWrite = taucLevelStatistics[nhyp].centralValue > TAUC_LEVEL_MAX_PLOT ? TAUC_LEVEL_MAX_PLOT : taucLevelStatistics[nhyp].centralValue;
                taucCentralValueWrite = taucCentralValueWrite < TAUC_LEVEL_MIN ? TAUC_LEVEL_MIN : taucCentralValueWrite;
                //
                double taucUpperBoundWrite = taucLevelStatistics[nhyp].upperBound > TAUC_LEVEL_MAX_PLOT ? TAUC_LEVEL_MAX_PLOT : taucLevelStatistics[nhyp].upperBound;
                taucUpperBoundWrite = taucUpperBoundWrite < TAUC_LEVEL_MIN ? TAUC_LEVEL_MIN : taucUpperBoundWrite;
                //
                double taucLowerBoundWrite = taucLevelStatistics[nhyp].lowerBound > TAUC_LEVEL_MAX_PLOT ? TAUC_LEVEL_MAX_PLOT : taucLevelStatistics[nhyp].lowerBound;
                taucLowerBoundWrite = taucLowerBoundWrite < TAUC_LEVEL_MIN ? TAUC_LEVEL_MIN : taucLowerBoundWrite;
                // check for enough tauc levels to plot statistics
                if (numTaucLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                    fprintf(taucCentralValueStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, taucCentralValueWrite);
                    fprintf(taucUpperBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, taucUpperBoundWrite);
                    fprintf(taucLowerBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, taucLowerBoundWrite);
                } else {
                    fprintf(taucCentralValueStream[nhyp], ">\n");
                    fprintf(taucUpperBoundStream[nhyp], ">\n");
                    fprintf(taucLowerBoundStream[nhyp], ">\n");
                }
                setLevelString(numTaucLevelMax[nhyp], &taucLevelStatistics[nhyp], taucLevelString[nhyp],
                        MIN_NUMBER_VALUES_USE, TAUC_RED_CUTOFF, TAUC_YELLOW_CUTOFF, -1.0, warning_colors_show);
                //
                // TdT50Ex level
                // NOTE: number of alarm levels is lessor of number t50 or tauc levels, this value has no statistical meaning with regards to alarm level.
                numAlarmLevel[nhyp] = IMIN(numT50ExLevel[nhyp], numTaucLevel[nhyp]);
                numAlarmLevelMax[nhyp] = IMIN(numT50ExLevelMax[nhyp], numTaucLevelMax[nhyp]);
                // TODO: Is it correct statistics to simply multiply central and bound values?
                /*if (numAlarmLevel[nhyp] >= MIN_NUMBER_VALUES_USE && numTaucLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                    tdT50ExLevelStatistics[nhyp].centralValue = t50ExLevelStatistics[nhyp].centralValue * taucLevelStatistics[nhyp].centralValue;
                    tdT50ExLevelStatistics[nhyp].lowerBound = t50ExLevelStatistics[nhyp].lowerBound * taucLevelStatistics[nhyp].lowerBound;
                    tdT50ExLevelStatistics[nhyp].upperBound = t50ExLevelStatistics[nhyp].upperBound * taucLevelStatistics[nhyp].upperBound;
                }*/
                if (numAlarmLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                    tdT50ExLevelStatistics[nhyp].centralValue = t50ExLevelStatistics[nhyp].centralValue * taucLevelStatistics[nhyp].centralValue;
                    tdT50ExLevelStatistics[nhyp].lowerBound = tdT50ExLevelStatistics[nhyp].centralValue
                            - stdDevProductNormalDist(t50ExLevelStatistics[nhyp].centralValue, t50ExLevelStatistics[nhyp].centralValue - t50ExLevelStatistics[nhyp].lowerBound,
                            taucLevelStatistics[nhyp].centralValue, taucLevelStatistics[nhyp].centralValue - taucLevelStatistics[nhyp].lowerBound);
                    tdT50ExLevelStatistics[nhyp].upperBound = tdT50ExLevelStatistics[nhyp].centralValue
                            + stdDevProductNormalDist(t50ExLevelStatistics[nhyp].centralValue, t50ExLevelStatistics[nhyp].upperBound - t50ExLevelStatistics[nhyp].centralValue,
                            taucLevelStatistics[nhyp].centralValue, taucLevelStatistics[nhyp].upperBound - taucLevelStatistics[nhyp].centralValue);
                    tdT50ExLevelStatistics[nhyp].numLevel = numAlarmLevel[nhyp];
                }
                //
                double tdT50ExCentralValueWrite = tdT50ExLevelStatistics[nhyp].centralValue > TDT50EX_LEVEL_MAX_PLOT ? TDT50EX_LEVEL_MAX_PLOT : tdT50ExLevelStatistics[nhyp].centralValue;
                tdT50ExCentralValueWrite = tdT50ExCentralValueWrite < TDT50EX_LEVEL_MIN ? TDT50EX_LEVEL_MIN : tdT50ExCentralValueWrite;
                //
                double tdT50ExUpperBoundWrite = tdT50ExLevelStatistics[nhyp].upperBound > TDT50EX_LEVEL_MAX_PLOT ? TDT50EX_LEVEL_MAX_PLOT : tdT50ExLevelStatistics[nhyp].upperBound;
                tdT50ExUpperBoundWrite = tdT50ExUpperBoundWrite < TDT50EX_LEVEL_MIN ? TDT50EX_LEVEL_MIN : tdT50ExUpperBoundWrite;
                //
                double tdT50ExLowerBoundWrite = tdT50ExLevelStatistics[nhyp].lowerBound > TDT50EX_LEVEL_MAX_PLOT ? TDT50EX_LEVEL_MAX_PLOT : tdT50ExLevelStatistics[nhyp].lowerBound;
                tdT50ExLowerBoundWrite = tdT50ExLowerBoundWrite < TDT50EX_LEVEL_MIN ? TDT50EX_LEVEL_MIN : tdT50ExLowerBoundWrite;
                // check for enough alarm levels to plot statistics
                if (numAlarmLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                    fprintf(tdT50ExCentralValueStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, tdT50ExCentralValueWrite);
                    fprintf(tdT50ExUpperBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, tdT50ExUpperBoundWrite);
                    fprintf(tdT50ExLowerBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, tdT50ExLowerBoundWrite);
                } else {
                    fprintf(tdT50ExCentralValueStream[nhyp], ">\n");
                    fprintf(tdT50ExUpperBoundStream[nhyp], ">\n");
                    fprintf(tdT50ExLowerBoundStream[nhyp], ">\n");
                }
                setLevelString(numAlarmLevelMax[nhyp], &tdT50ExLevelStatistics[nhyp], warningLevelString[nhyp],
                        MIN_NUMBER_VALUES_USE, TDT50EX_RED_CUTOFF, TDT50EX_YELLOW_CUTOFF, -1.0, warning_colors_show);
            }
            //
            if (flag_do_mwp) {
                //
                //mwp
                if (numMwpLevel[nhyp] > 0) {
                    // evaluate and set level statistics
                    setStatistics("Mwp", mwpStatisticsArray[nhyp], numMwpLevel[nhyp], &(mwpLevelStatistics[nhyp]),
                            0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                            && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
                }
                //else {
                //    mwpLevelStatistics[n].centralValue = mwpLevelStatistics[n].upperBound = mwpLevelStatistics[n].lowerBound = 0.0;
                //}
                double mwpCentralValueWrite = mwpLevelStatistics[nhyp].centralValue > MWP_LEVEL_MAX_PLOT ? MWP_LEVEL_MAX_PLOT : mwpLevelStatistics[nhyp].centralValue;
                mwpCentralValueWrite = mwpCentralValueWrite < MWP_LEVEL_MIN ? MWP_LEVEL_MIN : mwpCentralValueWrite;
                //
                double mwpUpperBoundWrite = mwpLevelStatistics[nhyp].upperBound > MWP_LEVEL_MAX_PLOT ? MWP_LEVEL_MAX_PLOT : mwpLevelStatistics[nhyp].upperBound;
                mwpUpperBoundWrite = mwpUpperBoundWrite < MWP_LEVEL_MIN ? MWP_LEVEL_MIN : mwpUpperBoundWrite;
                //
                double mwpLowerBoundWrite = mwpLevelStatistics[nhyp].lowerBound > MWP_LEVEL_MAX_PLOT ? MWP_LEVEL_MAX_PLOT : mwpLevelStatistics[nhyp].lowerBound;
                mwpLowerBoundWrite = mwpLowerBoundWrite < MWP_LEVEL_MIN ? MWP_LEVEL_MIN : mwpLowerBoundWrite;
#ifdef USE_MWP_LEVEL_ARRAY
                // check for enough mwp levels to plot statistics
                if (numMwpLevel[nhyp] >= MIN_NUMBER_VALUES_USE) {
                    fprintf(mwpCentralValueStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, mwpCentralValueWrite);
                    fprintf(mwpUpperBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, mwpUpperBoundWrite);
                    fprintf(mwpLowerBoundStream[nhyp], "%f %f\n", difftime(curr_time, plot_time_max) / 60.0, mwpLowerBoundWrite);
                } else {
                    fprintf(mwpCentralValueStream[nhyp], ">\n");
                    fprintf(mwpUpperBoundStream[nhyp], ">\n");
                    fprintf(mwpLowerBoundStream[nhyp], ">\n");
                }
#endif
                //
                setLevelString(numMwpLevelMax[nhyp], &mwpLevelStatistics[nhyp], mwpLevelString[nhyp],
                        MIN_NUMBER_VALUES_USE, MWP_RED_CUTOFF, MWP_YELLOW_CUTOFF, MWP_INVALID, magnitude_colors_show);
                //
            }
        }
        curr_time += TIME_STEP;

    }
    //
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        fclose_counter(pickStreamAssoc[nhyp]);
        fclose_counter(staAssociatedStream[nhyp]);
    }
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        fclose_counter(t50ExCentralValueStream[nhyp]);
        fclose_counter(t50ExUpperBoundStream[nhyp]);
        fclose_counter(t50ExLowerBoundStream[nhyp]);
    }
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            fclose_counter(taucCentralValueStream[nhyp]);
            fclose_counter(taucUpperBoundStream[nhyp]);
            fclose_counter(taucLowerBoundStream[nhyp]);
            fclose_counter(tdT50ExCentralValueStream[nhyp]);
            fclose_counter(tdT50ExUpperBoundStream[nhyp]);
            fclose_counter(tdT50ExLowerBoundStream[nhyp]);
        }
    }
#ifdef USE_MWP_LEVEL_ARRAY
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            fclose_counter(mwpCentralValueStream[nhyp]);
            fclose_counter(mwpUpperBoundStream[nhyp]);
            fclose_counter(mwpLowerBoundStream[nhyp]);
        }
    }
#endif
    fclose_counter(staT50RedStream);
    fclose_counter(staT50YellowStream);
    fclose_counter(staT50GreenStream);
    fclose_counter(staT50LtRedStream);
    fclose_counter(staT50LtYellowStream);
    fclose_counter(staT50LtGreenStream);
    fclose_counter(staTaucRedStream);
    fclose_counter(staTaucYellowStream);
    fclose_counter(staTaucGreenStream);
    fclose_counter(staTaucLtRedStream);
    fclose_counter(staTaucLtYellowStream);
    fclose_counter(staTaucLtGreenStream);




    // loop over associated hypocenters, develop statistics ===============================================================

    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {

        //
        if (flag_do_mb) {
            //
            //mb
            if (numMbLevelMax[nhyp] > 0) {
                // evaluate and set level statistics
                setStatistics("Mb", mbStatisticsArray[nhyp], numMbLevelMax[nhyp], &(mbLevelStatistics[nhyp]),
                        0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                        && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
            }
            setLevelString(numMbLevelMax[nhyp], &mbLevelStatistics[nhyp], mbLevelString[nhyp],
                    MIN_NUMBER_VALUES_USE, MB_RED_CUTOFF, MB_YELLOW_CUTOFF, MB_INVALID, magnitude_colors_show);
        }
        //
        if (flag_do_t0) {
            //
            //t0
            if (numT0LevelMax[nhyp] > 0) {
                // evaluate and set level statistics
                setStatistics("T0", t0StatisticsArray[nhyp], numT0LevelMax[nhyp], &(t0LevelStatistics[nhyp]),
                        0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                        && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
            }
            setLevelString(numT0LevelMax[nhyp], &t0LevelStatistics[nhyp], t0LevelString[nhyp],
                    MIN_NUMBER_VALUES_USE, T0_RED_CUTOFF, T0_YELLOW_CUTOFF, T0_INVALID, warning_colors_show);
            //
        }
        //

        // do Mwpd processing - must  be done outside main loops above because depends on T0 statistics
        if (flag_do_mwpd) {

            // loop over data
            int ndata;
            for (ndata = 0; ndata < num_de_data; ndata++) {
                // skip if not associated with this hypocenter
                TimedomainProcessingData* deData = data_list[ndata];
                if (deData->is_associated - 1 != nhyp)
                    continue;
                // skip if possible clip in data time span or if HF or BRB s/n ratio too low
                // AJL 20100224
                //if (deData->flag_clipped || deData->flag_snr_hf_too_low || deData->flag_a_ref_not_ok)
                if (deData->flag_clipped || deData->flag_non_contiguous || (!USE_AREF_NOT_OK_PICKS_FOR_LOCATION && deData->flag_a_ref_not_ok))
                    continue;
                // skip if data not used for location and not associated
                // AJL 20100630
                if (!deData->use_for_location && !(use_associated_data && is_associated_phase(deData)))
                    continue;
                // skip if not completed
                if (!deData->flag_complete_t50)
                    continue;
                //mwpd
                int indexMwpd = 0;
                int associated_ok = use_associated_data && is_associated_location_P(deData);
                if (flag_do_mwpd && deData->mwpd->raw_mag != MWPD_INVALID && !deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low) {// 20120612 AJL - changed s/n check from brb vel to brb disp (brb int)
                    if (!associate_data || !deData->is_associated || associated_ok) {
                        // 20110316 AJL - added taucLevelStatistics, checks against MIN_NUMBER_VALUES_USE
                        // 20121218 AJL - bug fix
                        /*double mwpd_corr_mag = calculate_corrected_Mwpd_Mag(deData->mwpd->raw_mag,
                                t0LevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? t0LevelStatistics[nhyp].centralValue : T0_INVALID,
                                taucLevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? taucLevelStatistics[nhyp].centralValue : T0_INVALID,
                                hyp_assoc_loc[nhyp]->depth);*/
                        double mwpd_corr_mag = calculate_corrected_Mwpd_Mag(deData->mwpd->raw_mag,
                                t0LevelStatistics[nhyp].numLevel >= MIN_NUMBER_VALUES_USE ? t0LevelStatistics[nhyp].centralValue : T0_INVALID,
                                taucLevelStatistics[nhyp].numLevel >= MIN_NUMBER_VALUES_USE ? taucLevelStatistics[nhyp].centralValue : TAUC_INVALID,
                                hyp_assoc_loc[nhyp]->depth);
                        deData->mwpd->corr_mag = mwpd_corr_mag;
                        if (deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MWPD && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MWPD
                                && useT0Report(deData) // 20120416 AJL
                                ) {
                            if (mwpd_corr_mag >= MWPD_LEVEL_MAX - MWPD_LEVEL_STEP)
                                indexMwpd = numMwpdLevelTotal - 1;
                            else if (mwpd_corr_mag <= MWPD_LEVEL_MIN)
                                indexMwpd = 0;
                            else
                                indexMwpd = (int) ((mwpd_corr_mag - MWPD_LEVEL_MIN) / MWPD_LEVEL_STEP);
                            // add mwpd to mwpd level statistics array - check associated and distance cutoff
                            if (!associate_data || associated_ok) {
                                mwpdStatisticsArray[(deData->is_associated - 1)][0][numMwpdLevel[(deData->is_associated - 1)]] = mwpd_corr_mag;
                                // 20110322 AJL - TEST making Mwpd corr after averaging - seems to make no difference to result
                                //mwpdStatisticsArray[(deData->is_associated - 1)][0][numMwpdLevel[(deData->is_associated - 1)]] = deData->mwpd->raw_mag;
                                mwpdStatisticsArray[(deData->is_associated - 1)][1][numMwpdLevel[(deData->is_associated - 1)]] = deData->lat;
                                mwpdStatisticsArray[(deData->is_associated - 1)][2][numMwpdLevel[(deData->is_associated - 1)]] = deData->lon;
                                numMwpdLevel[(deData->is_associated - 1)]++;
                                if (numMwpdLevel[(deData->is_associated - 1)] > numMwpdLevelMax[(deData->is_associated - 1)])
                                    numMwpdLevelMax[(deData->is_associated - 1)] = numMwpdLevel[(deData->is_associated - 1)];
#ifdef USE_MWP_MO_POS_NEG
                                if (deData->mwpd->mo_pos_neg_ratio > 0.0) {
                                    mwpdMoPosNegStatisticsArray[(deData->is_associated - 1)][0][numMwpdMoPosNegLevel[(deData->is_associated - 1)]] = deData->mwpd->mo_pos_neg_ratio;
                                    mwpdMoPosNegStatisticsArray[(deData->is_associated - 1)][1][numMwpdMoPosNegLevel[(deData->is_associated - 1)]] = deData->lat;
                                    mwpdMoPosNegStatisticsArray[(deData->is_associated - 1)][2][numMwpdMoPosNegLevel[(deData->is_associated - 1)]] = deData->lon;
                                    numMwpdMoPosNegLevel[(deData->is_associated - 1)]++;
                                    if (numMwpdMoPosNegLevel[(deData->is_associated - 1)] > numMwpdMoPosNegLevelMax[(deData->is_associated - 1)])
                                        numMwpdMoPosNegLevelMax[(deData->is_associated - 1)] = numMwpdMoPosNegLevel[(deData->is_associated - 1)];
                                }
#endif
                            }
                            (mwpdHistogram[(deData->is_associated - 1)][indexMwpd])++;
                        }
                    }
                    //int flag_active_mwpd = !associate_data || associated_ok;
                    //if (flag_active_mwpd) {
                    //    (mwpdHistogram[(deData->is_associated - 1)][indexMwpd])++;
                    //}
                }
            }
            //
            //mwpd
            if (numMwpdLevelMax[nhyp] > 0) {
                // evaluate and set level statistics
                setStatistics("Mwpd", mwpdStatisticsArray[nhyp], numMwpdLevelMax[nhyp], &(mwpdLevelStatistics[nhyp]),
                        0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                        && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
                // 20110322 AJL - TEST making Mwpd corr after averaging - seems to make no difference to result
                /*
                mwpdLevelStatistics[nhyp].centralValue = calculate_corrected_Mwpd_Mag(mwpdLevelStatistics[nhyp].centralValue,
                        t0LevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? t0LevelStatistics[nhyp].centralValue : T0_INVALID,
                        taucLevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? taucLevelStatistics[nhyp].centralValue : T0_INVALID,
                        hyp_assoc_loc[nhyp]->depth);
                mwpdLevelStatistics[nhyp].lowerBound = calculate_corrected_Mwpd_Mag(mwpdLevelStatistics[nhyp].lowerBound,
                        t0LevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? t0LevelStatistics[nhyp].centralValue : T0_INVALID,
                        taucLevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? taucLevelStatistics[nhyp].centralValue : T0_INVALID,
                        hyp_assoc_loc[nhyp]->depth);
                mwpdLevelStatistics[nhyp].upperBound = calculate_corrected_Mwpd_Mag(mwpdLevelStatistics[nhyp].upperBound,
                        t0LevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? t0LevelStatistics[nhyp].centralValue : T0_INVALID,
                        taucLevelStatistics[nhyp].centralValue >= MIN_NUMBER_VALUES_USE ? taucLevelStatistics[nhyp].centralValue : T0_INVALID,
                        hyp_assoc_loc[nhyp]->depth);*/

            }
            int nLevels = numMwpdLevelMax[nhyp];
            // check if Mwpd below minimum value
            if (mwpdLevelStatistics[nhyp].centralValue < MWPD_MIN_VALUE_USE)
                nLevels = -1; // force grey color
            setLevelString(nLevels, &mwpdLevelStatistics[nhyp], mwpdLevelString[nhyp],
                    MIN_NUMBER_VALUES_USE, MWPD_RED_CUTOFF, MWPD_YELLOW_CUTOFF, MWPD_INVALID, magnitude_colors_show);
            //
#ifdef USE_MWP_MO_POS_NEG
            if (numMwpdMoPosNegLevelMax[nhyp] > 0) {
                // evaluate and set level statistics
                setStatistics("MwpdMoPosNeg", mwpdMoPosNegStatisticsArray[nhyp], numMwpdMoPosNegLevelMax[nhyp], &(mwpdMoPosNegLevelStatistics[nhyp]),
                        0 && DEBUG && curr_time <= ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED)
                        && (curr_time + TIME_STEP) > ((time_t) hyp_assoc_loc[nhyp]->otime + LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED));
            }
#endif
        }


        // fill histograms

        FILE * t50ExHistorgramStream[num_hypocenters_associated];
        FILE * taucHistorgramStream[num_hypocenters_associated];
        FILE * mwpHistorgramStream[num_hypocenters_associated];
        FILE * mbHistorgramStream[num_hypocenters_associated];
        FILE * t0HistorgramStream[num_hypocenters_associated];
        FILE * mwpdHistorgramStream[num_hypocenters_associated];
        // GMT style
        sprintf(outname, "%s/plot/t50.%d.hist", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        t50ExHistorgramStream[nhyp] = fopen_counter(outname, "w");
        if (t50ExHistorgramStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        if (flag_do_tauc) {
            // histogram
            sprintf(outname, "%s/plot/tauc.%d.hist", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            taucHistorgramStream[nhyp] = fopen_counter(outname, "w");
            if (taucHistorgramStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
        }
        if (flag_do_mwp) {
            // histogram
            sprintf(outname, "%s/plot/mwp.%d.hist", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mwpHistorgramStream[nhyp] = fopen_counter(outname, "w");
            if (mwpHistorgramStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
        }
        if (flag_do_mb) {
            // histogram
            sprintf(outname, "%s/plot/mb.%d.hist", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mbHistorgramStream[nhyp] = fopen_counter(outname, "w");
            if (mbHistorgramStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
        }
        if (flag_do_t0) {
            // histogram
            sprintf(outname, "%s/plot/t0.%d.hist", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            t0HistorgramStream[nhyp] = fopen_counter(outname, "w");
            if (t0HistorgramStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
        }
        if (flag_do_mwpd) {
            // histogram
            sprintf(outname, "%s/plot/mwpd.%d.hist", outnameroot, nhyp);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            mwpdHistorgramStream[nhyp] = fopen_counter(outname, "w");
            if (mwpdHistorgramStream[nhyp] == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
        }

        // histograms
        double level = T50EX_LEVEL_MIN;
        for (m = 0; m < numWarningLevelTotal; m++) {
            int num = t50ExHistogram[nhyp][m];
            int m;
            for (m = 0; m < num; m++)
                fprintf(t50ExHistorgramStream[nhyp], " %f\n", level + 0.001);
            level += T50EX_LEVEL_STEP;
        }
        if (flag_do_tauc) {
            level = TAUC_LEVEL_MIN;
            for (m = 0; m < numTaucLevelTotal; m++) {
                int num = taucHistogram[nhyp][m];
                int m;
                for (m = 0; m < num; m++)
                    fprintf(taucHistorgramStream[nhyp], " %f\n", level + 0.001);
                level += TAUC_LEVEL_STEP;
            }
        }
        if (flag_do_mwp) {
            level = MWP_LEVEL_MIN;
            for (m = 0; m < numMwpLevelTotal; m++) {
                int num = mwpHistogram[nhyp][m];
                int m;
                for (m = 0; m < num; m++)
                    fprintf(mwpHistorgramStream[nhyp], " %f\n", level + 0.001);
                level += MWP_LEVEL_STEP;
            }
        }
        if (flag_do_mb) {
            level = MB_LEVEL_MIN;
            for (m = 0; m < numMbLevelTotal; m++) {
                int num = mbHistogram[nhyp][m];
                int m;
                for (m = 0; m < num; m++)
                    fprintf(mbHistorgramStream[nhyp], " %f\n", level + 0.001);
                level += MB_LEVEL_STEP;
            }
        }
        if (flag_do_t0) {
            level = T0_LEVEL_MIN;
            for (m = 0; m < numT0LevelTotal; m++) {
                int num = t0Histogram[nhyp][m];
                int m;
                for (m = 0; m < num; m++)
                    fprintf(t0HistorgramStream[nhyp], " %f\n", level + 0.001);
                level += T0_LEVEL_STEP;
            }
        }
        if (flag_do_mwpd) {
            level = MWPD_LEVEL_MIN;
            for (m = 0; m < numMwpdLevelTotal; m++) {
                int num = mwpdHistogram[nhyp][m];
                int m;
                for (m = 0; m < num; m++)
                    fprintf(mwpdHistorgramStream[nhyp], " %f\n", level + 0.001);
                level += MWPD_LEVEL_STEP;
            }
        }
        fclose_counter(t50ExHistorgramStream[nhyp]);
        if (flag_do_tauc)
            fclose_counter(taucHistorgramStream[nhyp]);
        if (flag_do_mwp)
            fclose_counter(mwpHistorgramStream[nhyp]);
        if (flag_do_mb)
            fclose_counter(mbHistorgramStream[nhyp]);
        if (flag_do_t0)
            fclose_counter(t0HistorgramStream[nhyp]);
        if (flag_do_mwpd)
            fclose_counter(mwpdHistorgramStream[nhyp]);
    }

    // station latency and health
    // 20141212 AJL - moved here from after line 5250: "// miscellaneous output ====" so health info is available for printHypoDataString())
    sprintf(outname, "%s/sta.health.html", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staHealthHtmlStream = fopen_counter(outname, "w");
    if (staHealthHtmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    updateStaHealthInformation(outnameroot, staHealthHtmlStream, difftime(time_max, time_min), (double) report_interval, time_since_last_report, time_min, time_max);
    fclose_counter(staHealthHtmlStream);

    // manage and write hypocenter data ===============================================================
    //
    // update hypocenter list
    // clear flags
    for (nhyp = 0; nhyp < num_hypocenters; nhyp++) {
        HypocenterDesc* phypo = hypo_list[nhyp];
        phypo->hyp_assoc_index = -1; // not associated
    }
    // update list with currently associated hypocenter
    int icheck_duplicates = 1;
    int is_new_hypocenter = 0;
    int have_new_hypocenter = 0;
    int have_new_hypocenter_with_mag = 0;
    int hypocenter_mail_sent = 0;
    // 20110407 AJL - Bug fix, reverse loop order to favor first associated hypocenters for case where same event is associated twice
    //for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
    for (nhyp = num_hypocenters_associated - 1; nhyp >= 0; nhyp--) {
        HypocenterDesc* phypo = hyp_assoc_loc[nhyp];
        // station health
        phypo->nstaHasBeenActive = nstaHasBeenActive;
        phypo->nstaIsActive = nstaIsActive;
        // magnitude and discriminants
        phypo->t50ExLevelStatistics = t50ExLevelStatistics[nhyp];
        phypo->taucLevelStatistics = taucLevelStatistics[nhyp];
        phypo->tdT50ExLevelStatistics = tdT50ExLevelStatistics[nhyp];
        strcpy(phypo->warningLevelString, warningLevelString[nhyp]);
        phypo->mwpLevelStatistics = mwpLevelStatistics[nhyp];
        phypo->mbLevelStatistics = mbLevelStatistics[nhyp];
        phypo->t0LevelStatistics = t0LevelStatistics[nhyp];
        phypo->mwpdLevelStatistics = mwpdLevelStatistics[nhyp];
#ifdef USE_MWP_MO_POS_NEG
        phypo->mwpdMoPosNegLevelStatistics = mwpdMoPosNegLevelStatistics[nhyp];
#endif
        phypo->hyp_assoc_index = nhyp; //  hypo_list is ordered by time, but association index is by association sum weight (lower index -> more phases associated), if re-associated, arbitrary, otherwise
        HypocenterDesc* phypocenter_desc_inserted = NULL;
        if ((is_new_hypocenter = addHypocenterDescToHypoList(phypo, &hypo_list, &num_hypocenters, icheck_duplicates, &existing_hypo_desc, &phypocenter_desc_inserted))) { // hypocenter unique_id is set here
            have_new_hypocenter++;
            phypocenter_desc_inserted->loc_seq_num = 0; // 20160905 AJL - added
        } else {
            if (phypocenter_desc_inserted->loc_type == LOC_TYPE_FULL || phypocenter_desc_inserted->loc_type == LOC_TYPE_RELOC_EXISTING) {
                (phypocenter_desc_inserted->loc_seq_num)++; // 20160905 AJL - added
            }
        }
        // copy inserted hypocenter into working hypocenter since may contain modified fields
        *phypo = *phypocenter_desc_inserted;
        if (!phypo->has_valid_magnitude) {
            if (phypo->mbLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE || phypo->mwpLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE
                    || phypo->mwpdLevelStatistics.numLevel >= MIN_NUMBER_VALUES_USE) {
                have_new_hypocenter_with_mag++;
                phypocenter_desc_inserted->has_valid_magnitude = phypo->has_valid_magnitude = 1;
            }
        }
        // create event map html pages
        create_map_html_page(outnameroot, phypo, time_min, time_max, MAP_LINK_MED_ZOOM);
        if (MAP_LINK_BIG_ZOOM != MAP_LINK_MED_ZOOM) {
            create_map_html_page(outnameroot, phypo, time_min, time_max, MAP_LINK_BIG_ZOOM);
        }
        // check for sending alert
        if (sendMail) {
            // pass phypocenter_desc_inserted here, since send_hypocenter_alert modifies fields in inserted hypocenter
            hypocenter_mail_sent += send_hypocenter_alert(phypocenter_desc_inserted, is_new_hypocenter, &existing_hypo_desc, sendMailParams, outnameroot, agencyId, verbose);
        }
    }
    // act on new hypocenters
    //
    // check for alarm notification
    sprintf(outname, "%s/alarm_notification.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * alarmNotificationStream = fopen_counter(outname, "w");
    if (alarmNotificationStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
    }
    int notification_level = 0;
    if (alarmNotification && (have_new_hypocenter || have_new_hypocenter_with_mag || hypocenter_mail_sent)) {
        fprintf(stdout, "%c%c%c", (char) 7, (char) 7, (char) 7);
        if (hypocenter_mail_sent) { // mail sent, important event
            fprintf(stdout, "%c", (char) 7);
            notification_level = 3;
        } else if (have_new_hypocenter_with_mag) {
            notification_level = 2;
        } else if (have_new_hypocenter) {
            notification_level = 1;
        }
    }
    if (alarmNotificationStream != NULL)
        fprintf(alarmNotificationStream, "%d", notification_level);
    fclose_counter(alarmNotificationStream);
    //
    FILE * hypocenterStream[num_hypocenters_associated];
    FILE * hypocenterOtimeStream[num_hypocenters_associated];
    FILE * hypocenterPwaveFrontStream[num_hypocenters_associated];
    FILE * hypocenterSwaveFrontStream[num_hypocenters_associated];
    FILE * hypocenterDistLimitsStream[num_hypocenters_associated];
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        // GMT style
        sprintf(outname, "%s/plot/hypo.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterStream[nhyp] = fopen_counter(outname, "w");
        if (hypocenterStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/hypo.%d.otime.dat", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterOtimeStream[nhyp] = fopen_counter(outname, "w");
        if (hypocenterOtimeStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/hypo.%d.pwavefront.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterPwaveFrontStream[nhyp] = fopen_counter(outname, "w");
        if (hypocenterPwaveFrontStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/hypo.%d.swavefront.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterSwaveFrontStream[nhyp] = fopen_counter(outname, "w");
        if (hypocenterSwaveFrontStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        // GMT style
        sprintf(outname, "%s/plot/hypo.%d.dist_limits.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterDistLimitsStream[nhyp] = fopen_counter(outname, "w");
        if (hypocenterDistLimitsStream[nhyp] == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
    }
    //
    // hypocenter
    FILE* hypocenterTextStream = NULL;
    if (num_hypocenters_associated > 0) {
        sprintf(outname, "%s/plot/hypo.txt", outnameroot);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        hypocenterTextStream = fopen_counter(outname, "w");
        if (hypocenterTextStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
    }
    //
    // 20160525 AJL - correction: added missing header for hypo csv files
    sprintf(outname, "%s/hypos.csv.hdr", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypocentersCsvHeaderStream = fopen_counter(outname, "w");
    if (hypocentersCsvHeaderStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    printHypoDataHeaderString(hypoDataString);
    fprintf(hypocentersCsvHeaderStream, "%s\n", hypoDataString); // 20160525 AJL - correction: added missing header for hypo csv files
    fclose_counter(hypocentersCsvHeaderStream);
    //
    sprintf(outname, "%s/hypos.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypocentersCsvStream = fopen_counter(outname, "w");
    if (hypocentersCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    //
    sprintf(outname, "%s/hypos_pretty.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypocentersCsvPrettyStream = fopen_counter(outname, "w");
    if (hypocentersCsvPrettyStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // hypoMessageHtmlStream
    char hypomessage_filename[] = "hypomessage.html";
    sprintf(outname, "%s/%s", outnameroot, hypomessage_filename);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypoMessageHtmlStream = fopen_counter(outname, "w");
    if (hypoMessageHtmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(hypoMessageHtmlStream, "<html>\n<body style=\"font-family:sans-serif;font-size:small\">\n<table border=0 cellpadding=1 frame=box rules=rows width=100%%>\n<tbody>\n");
    printHypoMessageHtmlHeaderString(hypoMessageHtmlString);
    fprintf(hypoMessageHtmlStream, "%s\n", hypoMessageHtmlString);
    //
    for (nhyp = num_hypocenters - 1; nhyp >= 0; nhyp--) { // reverse time order, since may remove hypocenters from list, and to give most recent at top of file lists.
        HypocenterDesc* phypo = hypo_list[nhyp];
        if (phypo->hyp_assoc_index >= 0) { // associated

            // hypo data csv string
            printHypoDataString(phypo, hypoDataString, 1);
            fprintf(hypocentersCsvStream, "%s\n", hypoDataString);
            printHypoDataString(phypo, hypoDataString, 0);
            fprintf(hypocentersCsvPrettyStream, "%s\n", hypoDataString);
            //
            printHypoMessageHtmlString(phypo, hypoMessageHtmlString, hypoBackgroundColor[phypo->hyp_assoc_index % hypoBackgroundColorModulo],
                    "bgcolor=#BBBBBB", phypo->hyp_assoc_index + 1, phypo->unique_id);
            fprintf(hypoMessageHtmlStream, "%s\n", hypoMessageHtmlString);
            // misc hypo streams
            fprintf(hypocenterStream[phypo->hyp_assoc_index], "%f %f\n", phypo->lat, phypo->lon);
            fprintf(hypocenterTextStream, "%f %f %d\n", phypo->lat, phypo->lon, phypo->hyp_assoc_index + 1);
            double dtime = difftime((time_t) phypo->otime, plot_time_max) / 60.0;
            //fprintf(hypocenterOtimeStream[phypo->hyp_assoc_index], ">\n%f %f\n%f %f\n", dtime, PICK_PLOT_LEVEL_MIN, dtime, PICK_PLOT_LEVEL_MAX);
            fprintf(hypocenterOtimeStream[phypo->hyp_assoc_index], "%f\n", dtime);

            // generate P and S wave fronts
            int p_phase_id = -1;
            int s_phase_id = -1;
            double time_since_origin = difftime(plot_time_max, (time_t) phypo->otime);
            double p_wave_front_dist = simple_P_distance(time_since_origin, phypo->depth, &p_phase_id); // distance in degrees
            double s_wave_front_dist = simple_S_distance(time_since_origin, phypo->depth, &s_phase_id); // distance in degrees
            //printf("   --->    plot_time_max=%ld hyp.ot=%ld  diff=%.1f p_wave_front_dist: %s = %.1f s_wave_front_dist: %s = %.1f\n",
            //        plot_time_max, (time_t) phypo->otime, time_since_origin,
            //        phase_name_for_id(p_phase_id), p_wave_front_dist, phase_name_for_id(s_phase_id), s_wave_front_dist);
            // P wavefront
            if (p_wave_front_dist > 1.0 && p_wave_front_dist < get_dist_time_dist_max()) {
                //double pfront_az_step = 1.0;  // 20110318 AJL
                double pfront_az_step = 30.0 / p_wave_front_dist;
                //printf("   --->    pfront_az_step=%.1f\n", pfront_az_step);
                double pfront_az = 91.0; // 20110203 AJL - do not start at N or S, and avoid plotting at exactly 0, 90, 180, 270 - may cause problem in GMT ???
                double pf_lat, pf_lon;
                fprintf(hypocenterPwaveFrontStream[phypo->hyp_assoc_index], "> %s\n", phase_name_for_id(p_phase_id));
                while (pfront_az < 91.0 + 360.0 + pfront_az_step + pfront_az_step / 2.0) {
                    PointAtGCDistanceAzimuth(phypo->lat, phypo->lon, p_wave_front_dist, pfront_az, &pf_lat, &pf_lon);
                    if (fabs(pf_lat) < 89.5) // 20110103 AJL - avoid potential problems with latitudes near poles
                        fprintf(hypocenterPwaveFrontStream[phypo->hyp_assoc_index], "%f %f\n", pf_lat, pf_lon);
                    pfront_az += pfront_az_step;
                }
            }
            // S wavefront
            if (s_wave_front_dist > 1.0 && s_wave_front_dist < get_dist_time_dist_max()) {
                //double sfront_az_step = 1.0;  // 20110318 AJL
                double sfront_az_step = 30.0 / s_wave_front_dist;
                //printf("   --->    sfront_az_step=%.1f\n", sfront_az_step);
                double sfront_az = 91.0; // 20110203 AJL - do not start at N or S, and avoid plotting at exactly 0, 90, 180, 270 - may cause problem in GMT ???
                double sf_lat, sf_lon;
                fprintf(hypocenterSwaveFrontStream[phypo->hyp_assoc_index], "> %s\n", phase_name_for_id(s_phase_id));
                while (sfront_az < 91.0 + 360.0 + sfront_az_step + sfront_az_step / 2.0) {
                    PointAtGCDistanceAzimuth(phypo->lat, phypo->lon, s_wave_front_dist, sfront_az, &sf_lat, &sf_lon);
                    if (fabs(sf_lat) < 89.5) // 20110103 AJL - avoid potential problems with latitudes near poles
                        fprintf(hypocenterSwaveFrontStream[phypo->hyp_assoc_index], "%f %f\n", sf_lat, sf_lon);
                    sfront_az += sfront_az_step;
                }
            }

            // calculate ratio of nsta associated P to nsta available within P wavefront
            // 20140818 AJL - added to test check for false large events, this ratio should be large for large events
            int phase_id_found = -1;
            double dist_max = simple_distance(time_since_origin, phypo->depth, "P", &phase_id_found);
            if (dist_max < 0.0)
                dist_max = get_dist_time_dist_max();
            int nStaAvailable = countStationsAvailable(phypo->lat, phypo->lon, dist_max, channelParameters);
            double ratioNasscP2NStaAvailable = -1.0;
            if (nStaAvailable > 0) {
                ratioNasscP2NStaAvailable = (double) phypo->nassoc_P / (double) nStaAvailable;
            }
            printf("INFO TEST: Pdist=%.1f, phase=%s, nassocP=%d, nStaAvailable=%d, ratio_nassocP_2_nStaAvailable=%.3f,\n",
                    dist_max, phase_name_for_id(phase_id_found), phypo->nassoc_P, nStaAvailable, ratioNasscP2NStaAvailable);

            // distance limits T50Ex, Td, Mwp
            /*  2011222 replaced with 5 deg limit only, see below
            for (n = 0; n < 4; n++) {
                double dl_dist;
                char dl_name[64];
                if (n == 0) {
                    dl_dist = MIN_EPICENTRAL_DISTANCE_WARNING;
                    sprintf(dl_name, "%d-T50-min", (int) (0.5 + MIN_EPICENTRAL_DISTANCE_WARNING));
                } else if (n == 1) {
                    dl_dist = MAX_EPICENTRAL_DISTANCE_WARNING;
                    sprintf(dl_name, "%d-T50-Mwp-max", (int) (0.5 + MAX_EPICENTRAL_DISTANCE_WARNING)); // as long as MAX_EPICENTRAL_DISTANCE_WARNING == MAX_EPICENTRAL_DISTANCE_MWP
                } else if (n == 2) {
                    dl_dist = MIN_EPICENTRAL_DISTANCE_TAUC;
                    sprintf(dl_name, "%d-Td-Mwp-min", (int) (0.5 + MIN_EPICENTRAL_DISTANCE_TAUC)); // as long as MIN_EPICENTRAL_DISTANCE_TAUC == MIN_EPICENTRAL_DISTANCE_MWP
                } else if (n == 3) {
                    dl_dist = MAX_EPICENTRAL_DISTANCE_TAUC;
                    sprintf(dl_name, "%d-Td-max", (int) (0.5 + MAX_EPICENTRAL_DISTANCE_TAUC));
                }
                //double dl_az_step = 1.0;  // 20110318 AJL
                double dl_az_step = 30.0 / dl_dist;
                double dl_az = 91.0; // 20110203 AJL - do not start at N or S, and avoid plotting at exactly 0, 90, 180, 270 - may cause problem in GMT ???
                double dl_lat, dl_lon;
                fprintf(hypocenterDistLimitsStream[phypo->hyp_assoc_index], "> %s\n", dl_name);
                while (dl_az < 91.0 + 360.0 + dl_az_step / 2.0) {
                    PointAtGCDistanceAzimuth(phypo->lat, phypo->lon, dl_dist, dl_az, &dl_lat, &dl_lon);
                    //if (dl_lon < -180.0)
                    //    dl_lon += 360.0;
                    //if (dl_lon > 180.0)
                    //    dl_lon -= 360.0; // 20110103
                    if (fabs(dl_lat) < 89.5) // 20110103 AJL - avoid potential problems with latitudes near poles
                        fprintf(hypocenterDistLimitsStream[phypo->hyp_assoc_index], "%f %f\n", dl_lat, dl_lon);
                    dl_az += dl_az_step;
                }
            }*/
            // distance limit at 5 deg only
            if (1) {
                double dl_dist;
                char dl_name[64];
                dl_dist = 5.0;
                sprintf(dl_name, "5deg");
                //double dl_az_step = 1.0;  // 20110318 AJL
                double dl_az_step = 30.0 / dl_dist;
                double dl_az = 91.0; // 20110203 AJL - do not start at N or S, and avoid plotting at exactly 0, 90, 180, 270 - may cause problem in GMT ???
                double dl_lat, dl_lon;
                fprintf(hypocenterDistLimitsStream[phypo->hyp_assoc_index], "> %s\n", dl_name);
                while (dl_az < 91.0 + 360.0 + dl_az_step / 2.0) {
                    PointAtGCDistanceAzimuth(phypo->lat, phypo->lon, dl_dist, dl_az, &dl_lat, &dl_lon);
                    /*if (dl_lon < -180.0)
                        dl_lon += 360.0;
                    if (dl_lon > 180.0)
                        dl_lon -= 360.0;*/ // 20110103
                    if (fabs(dl_lat) < 89.5) // 20110103 AJL - avoid potential problems with latitudes near poles
                        fprintf(hypocenterDistLimitsStream[phypo->hyp_assoc_index], "%f %f\n", dl_lat, dl_lon);
                    dl_az += dl_az_step;
                }
            }
        }
    }
    //
    fclose_counter(hypocentersCsvStream);
    fclose_counter(hypocentersCsvPrettyStream);
    fprintf(hypoMessageHtmlStream, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose_counter(hypoMessageHtmlStream);
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        fclose_counter(hypocenterStream[nhyp]);
        fclose_counter(hypocenterOtimeStream[nhyp]);
        fclose_counter(hypocenterPwaveFrontStream[nhyp]);
        fclose_counter(hypocenterSwaveFrontStream[nhyp]);
        fclose_counter(hypocenterDistLimitsStream[nhyp]);
    }
    if (num_hypocenters_associated > 0)
        fclose_counter(hypocenterTextStream);
    //
    // xml hypocenters messsage
    sprintf(xmlWriterUri, "%s/monitor.xml", outnameroot);
    int iWriteArrivals = 1;
    int iWriteUnAssociatedPicks = 1;
    writeLocXML(xmlWriterUri, time_max, agencyId, hypo_list, num_hypocenters, data_list, num_de_data, iWriteArrivals, iWriteUnAssociatedPicks, printIgnoredData);


    // loop over associated hypocenters, develop epicenter statistics ===============================================================
    // 20150812 AJL - added

    // open epicenter statistics information csv file
    FILE * epicenterDiffCsvStream = NULL;
    sprintf(outname, "%s/epicenter.diff.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    epicenterDiffCsvStream = fopen_counter(outname, "w");
    if (epicenterDiffCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(epicenterDiffCsvStream, "epicenter.diff");
    fprintf(epicenterDiffCsvStream, " %lf", EPI_DIFF_CRITICAL_VALUE);
    fprintf(epicenterDiffCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(epicenterDiffCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(epicenterDiffCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(epicenterDiffCsvStream, " %d", report_interval);
    fprintf(epicenterDiffCsvStream, " %d", nstaIsActive);
    fprintf(epicenterDiffCsvStream, " %d", nstaHasBeenActive);
    fprintf(epicenterDiffCsvStream, " %d", num_hypocenters_associated);
    fprintf(epicenterDiffCsvStream, "\n");
    char epiDiffLevelString[num_hypocenters_associated][WARNING_LEVEL_STRING_LEN];

    // hyp temp directories
    static char hyp_dir_name[1024];
    sprintf(hyp_dir_name, "%s/hyp", EE_TEMP_DIR); // where persistent hyp info is kept
    mkdir(hyp_dir_name, 0755);
    static char hyp_dir_scratch_name[1024];
    sprintf(hyp_dir_scratch_name, "%s/hyp_scratch", EE_TEMP_DIR); // scratch/working dir for active hypos
    mkdir(hyp_dir_scratch_name, 0755);
    static char outname_tmp[1024];

    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        HypocenterDesc* phypo = hyp_assoc_loc[nhyp];
        //
        // save current hypocenter
        sprintf(outname, "%s/hypocenter.history.value.%ld.xy", hyp_dir_name, phypo->unique_id);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        FILE * hypocenterHistoryStream = fopen_counter(outname, "a");
        if (hypocenterHistoryStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        fprintf(hypocenterHistoryStream, "%ld %lf %lf %lf %lf %lf\n", plot_time_max,
                phypo->lat, phypo->lon, phypo->errh, phypo->depth, phypo->otime);
        fclose_counter(hypocenterHistoryStream);
        // move hyp history file to hyp_tmp dir
        sprintf(outname_tmp, "%s/hypocenter.history.value.%ld.xy", hyp_dir_scratch_name, phypo->unique_id);
        rename(outname, outname_tmp);
        //
        // reconstruct epicenter diff relative to hypocenter history
        // re-open hypocenter history file to read
        hypocenterHistoryStream = fopen_counter(outname_tmp, "r");
        if (hypocenterHistoryStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname_tmp);
            perror(tmp_str);
            return (-1);
        }
        // open epicenter diff files
        sprintf(outname, "%s/plot/epicenter.diff.value.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        FILE * epicenterDiffValueStream = fopen_counter(outname, "w");
        if (epicenterDiffValueStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        sprintf(outname, "%s/plot/epicenter.diff.upper.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        FILE * epicenterDiffUpperStream = fopen_counter(outname, "w");
        if (epicenterDiffUpperStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        sprintf(outname, "%s/plot/epicenter.diff.lower.%d.xy", outnameroot, nhyp);
        if (verbose > 2)
            printf("Opening output file: %s\n", outname);
        FILE * epicenterDiffLowerStream = fopen_counter(outname, "w");
        if (epicenterDiffLowerStream == NULL) {
            sprintf(tmp_str, "ERROR: opening output file: %s", outname);
            perror(tmp_str);
            return (-1);
        }
        fprintf(epicenterDiffValueStream, ">\n");
        fprintf(epicenterDiffUpperStream, ">\n");
        fprintf(epicenterDiffLowerStream, ">\n");
        int istat;
        double histtime = 0.0, lat = 0.0, lon = 0.0, errh = 0.0, depth = 0.0, otime = 0.0, depicenter = 0.0;
        while (1) {
            if (fgets(tmp_str, STANDARD_STRLEN - 1, hypocenterHistoryStream) == NULL)
                break;
            istat = sscanf(tmp_str, "%lf %lf %lf %lf %lf %lf ", &histtime, &lat, &lon, &errh, &depth, &otime);
            depicenter = DEG2KM * GCDistance(lat, lon, phypo->lat, phypo->lon);
            fprintf(epicenterDiffValueStream, "%f %f\n", difftime(histtime, plot_time_max) / 60.0, depicenter);
            fprintf(epicenterDiffUpperStream, "%f %f\n", difftime(histtime, plot_time_max) / 60.0, depicenter + errh);
            fprintf(epicenterDiffLowerStream, "%f %f\n", difftime(histtime, plot_time_max) / 60.0, depicenter - errh);
        }
        fclose_counter(hypocenterHistoryStream);
        fclose_counter(epicenterDiffValueStream);
        fclose_counter(epicenterDiffUpperStream);
        fclose_counter(epicenterDiffLowerStream);
        //
        // epicenter statistics information csv file
        strcpy(epiDiffLevelString[nhyp], "NONE");
        double level_plot_epi = depicenter < EPI_DIFF_LEVEL_MAX ? depicenter : EPI_DIFF_LEVEL_MAX;
        level_plot_epi = level_plot_epi > EPI_DIFF_LEVEL_MIN ? level_plot_epi : EPI_DIFF_LEVEL_MIN;
        double level_plot_epi_lower = depicenter - errh < EPI_DIFF_LEVEL_MAX ? depicenter - errh : EPI_DIFF_LEVEL_MAX;
        level_plot_epi_lower = level_plot_epi_lower > EPI_DIFF_LEVEL_MIN ? level_plot_epi_lower : EPI_DIFF_LEVEL_MIN;
        double level_plot_epi_upper = depicenter + errh < EPI_DIFF_LEVEL_MAX ? depicenter + errh : EPI_DIFF_LEVEL_MAX;
        level_plot_epi_upper = level_plot_epi_upper > EPI_DIFF_LEVEL_MIN ? level_plot_epi_upper : EPI_DIFF_LEVEL_MIN;
        fprintf(epicenterDiffCsvStream, " %.1f", level_plot_epi);
        fprintf(epicenterDiffCsvStream, " %.1f", depicenter);
        fprintf(epicenterDiffCsvStream, " %.1f", level_plot_epi_lower);
        fprintf(epicenterDiffCsvStream, " %.1f", level_plot_epi_upper);
        fprintf(epicenterDiffCsvStream, " %s", epiDiffLevelString[nhyp]);
        fprintf(epicenterDiffCsvStream, " %d", phypo->nassoc_P);
        fprintf(epicenterDiffCsvStream, "\n");
    }
    fclose_counter(epicenterDiffCsvStream);

    // clean up old temp files
    int ireturn = nftw(hyp_dir_name, remove_fn, 16, FTW_DEPTH);
    if (ireturn) {
        printf("ERROR: removing files: return value = %d, path = %s\n", ireturn, hyp_dir_name);
    }
    // move new temp files back to persistent temp dir
    rename(hyp_dir_scratch_name, hyp_dir_name);



    // loop over pick data, write output ===============================================================

    //
    sprintf(outname, "%s/pickdata_nlloc.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * pickDataNLLStream = fopen_counter(outname, "w");
    if (pickDataNLLStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // pickMessageHtmlStream
    sprintf(outname, "%s/pickmessage.html", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * pickMessageHtmlStream = fopen_counter(outname, "w");
    if (pickMessageHtmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(pickMessageHtmlStream, "<html>\n<body style=\"font-family:sans-serif;font-size:small\">\n<table border=0 cellpadding=1 frame=box rules=rows width=100%%>\n<tbody>\n");
    fprintf(pickMessageHtmlStream, "<tr align=right bgcolor=\"#BBBBBB\">\n<th>&nbsp;n</th>");
    if (associate_data) {
        fprintf(pickMessageHtmlStream, "<th>&nbsp;evt</th><th>&nbsp;d&nbsp;<br>&nbsp;(deg)</th><th>&nbsp;az&nbsp;<br>&nbsp;(deg)</th>");
    }
    fprintf(pickMessageHtmlStream, "<th>&nbsp;channel</th><th>&nbsp;stream</th><th>&nbsp;loc</th><th>&nbsp;time<br>&nbsp;(UTC)</th><th>&nbsp;unc<br>&nbsp;(sec)</th><th>&nbsp;pol</th><th>&nbsp;_ty</th><th>&nbsp;_wt</th><th>&nbsp;toang<br>&nbsp;(deg)</th><th>&nbsp;paz<br>&nbsp;(deg)</th><th>&nbsp;_unc<br>&nbsp;(deg)</th><th>&nbsp;_calc<br>&nbsp;(deg)</th><th>&nbsp;_wt</th>");
    if (associate_data) {
        fprintf(pickMessageHtmlStream, "<th>&nbsp;phase</th><th>&nbsp;res&nbsp;<br>&nbsp;(sec)</th><th>&nbsp;tot_wt</th><th>&nbsp;dist_wt</th><th>&nbsp;st_q_wt</th>");
    }
    fprintf(pickMessageHtmlStream, "<th>&nbsp;T50</th><th>&nbsp;Aref</th><th>&nbsp;Aerr</th><th>&nbsp;T50Ex</th><th>&nbsp;Td&nbsp;<br>&nbsp;(sec)</th>");
    fprintf(pickMessageHtmlStream, "<th>&nbsp;Avel</th><th>&nbsp;Adisp</th>");
    fprintf(pickMessageHtmlStream, "<th>&nbsp;s/n_HF</th>");
    fprintf(pickMessageHtmlStream, "<th>&nbsp;s/n_BRBV</th><th>&nbsp;s/n_BRBD</th><th>&nbsp;mb</th><th>&nbsp;Mwp</th><th>&nbsp;T0&nbsp;<br>&nbsp;(sec)</th><th>&nbsp;Mwpd</th><th>&nbsp;status</th>");
    fprintf(pickMessageHtmlStream, "</tr>\n");
    // picksCsvStream
    sprintf(outname, "%s/picks.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * picksCsvStream = fopen_counter(outname, "w");
    if (picksCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    if (associate_data) {
        fprintf(picksCsvStream, "event_id n event dist az channel stream loc time unc pol pol_type pol_wt toang paz paz_unc paz_calc paz_wt ");
        fprintf(picksCsvStream, "phase residual tot_wt dist_wt st_q_wt T50 Aref Aerr T50Ex Tdom ");
        fprintf(picksCsvStream, "Avel Adisp ");
        fprintf(picksCsvStream, "s/n_HF s/n_BRBV s/n_BRBD mb Mwp T0 Mwpd status sta_corr\n");
    } else {
        fprintf(picksCsvStream, "n channel stream time T50 Aref T50Ex Tdom Avel s/n_HF s/n_BRBV s/n_BRBD status sta_corr\n");
    }
    // write pick data to file
    int ndata;
    // NLL pick data in time order
    for (ndata = 0; ndata < num_de_data; ndata++) {
        TimedomainProcessingData* deData = data_list[ndata];
        if (printIgnoredData || !ignoreData(deData)) {
            if (deData->t_time_t >= time_min) { // 20130407 AJL - added to prevent writing of picks from interval time_min-report_interval->time_min.
                fprintf_NLLoc_TimedomainProcessingData(deData, pickDataNLLStream, 1);
                if (deData->is_associated > 0) {
                    fprintf(pickDataNLLStream, " %ld\n", hyp_assoc_loc[deData->is_associated - 1]->unique_id);
                } else {
                    fprintf(pickDataNLLStream, " %d\n", -1);
                }
            }
        }
    }
    // general pick data in reverse time order
    static char rowBackgroundColor[64];
    static char t50ExBackgroundColor[64];
    static char taucBackgroundColor[64];
    static char mbBackgroundColor[64];
    static char mwpBackgroundColor[64];
    static char t0BackgroundColor[64];
    static char mwpdBackgroundColor[64];
    int num_pick = 0;
    int num_pick_assoc = 0;
    int num_pick_not_completed = 0;
    int num_pick_clipped = 0;
    int num_pick_non_contiguous = 0;
    int num_pick_sn_low = 0;
    int num_pick_sn_brb_low = 0;
    int num_pick_sn_brb_int_low = 0;
    int num_pick_in_prev_coda = 0;
    int num_pick_ok = 0;
    int num_event = 0;
    // 2010116 AJL for (ndata = 0; ndata < num_de_data; ndata++) {
    for (ndata = num_de_data - 1; ndata >= 0; ndata--) { // reverse time order
        num_pick++;
        int pick_ok = 1;
        TimedomainProcessingData* deData = data_list[ndata];
        if (deData->flag_clipped) {
            num_pick_clipped++;
            pick_ok = 0;
        }
        if (deData->flag_non_contiguous) {
            num_pick_non_contiguous++;
            pick_ok = 0;
        }
        if (deData->flag_snr_hf_too_low) {
            num_pick_sn_low++;
        }
        if (deData->flag_snr_brb_too_low) {
            num_pick_sn_brb_low++;
        }
        if (deData->flag_snr_brb_int_too_low) {
            num_pick_sn_brb_int_low++;
        }
        // 20130128 AJL - use flag_snr_brb_int_too_low to enable mwp, mwpd, etc., but do not use for ignore tests (e.g. ignore determined by flag_snr_brb_too_low)
        //if (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low && deData->flag_snr_brb_int_too_low) {
        // 20131022 AJL - try using all picks for location, regardless of HF S/N
        // 20150408 AJL  if (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low && !(USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION && deData->is_associated)) {
        if (!USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION && deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low) {
            pick_ok = 0;
        }
        if (deData->flag_a_ref_not_ok) {
            num_pick_in_prev_coda++;
            if (!USE_AREF_NOT_OK_PICKS_FOR_LOCATION) {
                pick_ok = 0;
            }
        }
        if (!deData->flag_complete_t50) {
            num_pick_not_completed++;
            pick_ok = 0;
        }
        num_pick_ok += pick_ok;
        if (printIgnoredData || !ignoreData(deData)) {
            double deLevel = -1.0;
            char deLevelStr[64];
            sprintf(deLevelStr, "-1");
            double tauc_peak = TAUC_INVALID;
            char tauc_peakStr[64];
            sprintf(tauc_peakStr, "-1");
            double mb_mag = MB_INVALID;
            double mwp_mag = MWP_INVALID;
            double t0_dur = T0_INVALID;
            double mwpd_corr_mag = MWPD_INVALID;
            // readable form
            strcpy(rowBackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(t50ExBackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(taucBackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(mbBackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(mwpBackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(t0BackgroundColor, "bgcolor=\"#FFFFFF\"");
            strcpy(mwpdBackgroundColor, "bgcolor=\"#FFFFFF\"");
            int ignored = 0;
            if (ignoreData(deData)) {
                ignored = 1;
                strcpy(rowBackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(t50ExBackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(taucBackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(mbBackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(mwpBackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(t0BackgroundColor, "bgcolor=\"#EEEEEE\"");
                strcpy(mwpdBackgroundColor, "bgcolor=\"#EEEEEE\"");
            } else if (!deData->flag_complete_t50) {
                strcpy(rowBackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(t50ExBackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(taucBackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(mbBackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(mwpBackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(t0BackgroundColor, "bgcolor=\"#CCCCCC\"");
                strcpy(mwpdBackgroundColor, "bgcolor=\"#CCCCCC\"");
            } else {
                int not_associated = !deData->is_associated || deData->is_associated == NUMBER_ASSOCIATE_IGNORE;
                deLevel = getT50Level(deData);
                sprintf(deLevelStr, "%.2f", deLevel);
                if (!deData->flag_snr_hf_too_low) {
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_WARNING && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_WARNING);
                    if (!associated_ok) {
                        strcpy(t50ExBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (deLevel >= T50EX_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (deLevel >= T50EX_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (deLevel >= T50EX_LEVEL_MIN) {
                        if (not_associated)
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(t50ExBackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(t50ExBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(t50ExBackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
                if (flag_do_tauc && !deData->flag_snr_brb_too_low) {
                    tauc_peak = deData->tauc_peak;
                    sprintf(tauc_peakStr, "%.2f", tauc_peak);
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_TAUC && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_TAUC);
                    if (!associated_ok) {
                        strcpy(taucBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (tauc_peak >= TAUC_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(taucBackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(taucBackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (tauc_peak >= TAUC_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(taucBackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(taucBackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (tauc_peak >= TAUC_LEVEL_MIN) {
                        if (not_associated)
                            strcpy(taucBackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(taucBackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(taucBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(taucBackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
                if (flag_do_mb && !deData->flag_snr_brb_bp_too_low) {
                    mb_mag = deData->mb->mag;
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MB && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MB);
                    if (!associated_ok) {
                        strcpy(mbBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (mb_mag >= MB_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(mbBackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(mbBackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (mb_mag >= MB_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(mbBackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(mbBackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (mb_mag != MB_INVALID) {
                        if (not_associated)
                            strcpy(mbBackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(mbBackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(mbBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(mbBackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
                // 20121119 AJL if (flag_do_mwp && !deData->flag_snr_brb_too_low) {
                if (flag_do_mwp && !deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low) {
                    mwp_mag = deData->mwp->mag;
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MWP && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MWP);
                    if (!associated_ok) {
                        strcpy(mwpBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (mwp_mag >= MWP_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(mwpBackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(mwpBackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (mwp_mag >= MWP_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(mwpBackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(mwpBackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (mwp_mag != MWP_INVALID) {
                        if (not_associated)
                            strcpy(mwpBackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(mwpBackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(mwpBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(mwpBackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
                if (flag_do_t0 && !deData->flag_snr_hf_too_low) {
                    //t0_dur = deData->t0->duration_raw;
                    //double depth = 0.0;
                    //if (deData->is_associated)
                    //    depth = hyp_assoc_loc[deData->is_associated - 1]->depth;
                    //t0_dur = calculate_corrected_duration(deData, depth); // 20111222 TEST AJL - use S duration
                    t0_dur = deData->t0->duration_plot;
                    //if (t0_dur < 0.0)
                    //    printf("ERROR: deData->t0->duration_plot not set: this should not happen!\n");
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && useT0Report(deData));
                    if (!associated_ok) {
                        strcpy(t0BackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (t0_dur >= T0_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(t0BackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(t0BackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (t0_dur >= T0_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(t0BackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(t0BackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (t0_dur != T0_INVALID) {
                        if (not_associated)
                            strcpy(t0BackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(t0BackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(t0BackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(t0BackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
                // 20121119 AJL if (flag_do_mwpd && !deData->flag_snr_hf_too_low && !deData->flag_snr_brb_int_too_low) { // 20120612 AJL - changed s/n check from brb vel to brb disp (brb int)
                if (flag_do_mwpd && !deData->flag_snr_brb_too_low && !deData->flag_snr_brb_int_too_low) {
                    mwpd_corr_mag = deData->mwpd->corr_mag;
                    int associated_ok = !use_associated_data || not_associated || (is_associated_location_P(deData)
                            && deData->epicentral_distance > MIN_EPICENTRAL_DISTANCE_MWPD && deData->epicentral_distance < MAX_EPICENTRAL_DISTANCE_MWPD
                            && useT0Report(deData) // 20120416 AJL
                            );
                    // 20130208 AJL  if (!associated_ok || mwpd_corr_mag < MWPD_MIN_VALUE_USE) {
                    if (!associated_ok) {
                        strcpy(mwpdBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    } else if (mwpd_corr_mag >= MWPD_RED_CUTOFF) {
                        if (not_associated)
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#FFBBBB\"");
                        else
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#FF5555\"");
                    } else if (mwpd_corr_mag >= MWPD_YELLOW_CUTOFF) {
                        if (not_associated)
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#FFFFCC\"");
                        else
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#FFFF66\"");
                    } else if (mwpd_corr_mag != MWPD_INVALID) {
                        if (not_associated)
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#DDFFDD\"");
                        else
                            strcpy(mwpdBackgroundColor, "bgcolor=\"#77FF77\"");
                    } else {
                        strcpy(mwpdBackgroundColor, "bgcolor=\"#CCCCCC\"");
                    }
                } else {
                    strcpy(mwpdBackgroundColor, "bgcolor=\"#EEEEEE\"");
                }
            }
            if (associate_data) {
                if (use_associated_data && is_associated_phase(deData)) {
                    num_pick_assoc++;
                    if (deData->is_associated > num_event)
                        num_event = deData->is_associated;
                    fprintf(pickMessageHtmlStream, "<tr align=right %s><td>%d</td>", hypoBackgroundColor[(deData->is_associated - 1) % hypoBackgroundColorModulo], ndata + 1);
                    create_event_link(".", hyp_assoc_loc[deData->is_associated - 1]->unique_id, event_url_str, event_link_str,
                            feregion(hyp_assoc_loc[deData->is_associated - 1]->lat, hyp_assoc_loc[deData->is_associated - 1]->lon,
                            feregion_str, FEREGION_STR_SIZE));
                    fprintf(pickMessageHtmlStream, "<td>%s%d</a></td><td>%.1f</td><td>%.0f</td>",
                            event_link_str, deData->is_associated, deData->epicentral_distance, deData->epicentral_azimuth);
                    fprintf(picksCsvStream, "%ld %d ", hyp_assoc_loc[deData->is_associated - 1]->unique_id, ndata + 1);
                    fprintf(picksCsvStream, "%d %.2f %.1f ", deData->is_associated, deData->epicentral_distance, deData->epicentral_azimuth);
                } else {
                    fprintf(pickMessageHtmlStream, "<tr align=right %s><td>%d</td>", rowBackgroundColor, ndata + 1);
                    fprintf(pickMessageHtmlStream, "<td></td><td></td><td></td>");
                    fprintf(picksCsvStream, "-1 %d ", ndata + 1);
                    fprintf(picksCsvStream, "-1 -1 -1 ");
                }
            }
            // channel identification
            char loc_chr = '-';
            if (deData->use_for_location) {
                if (deData->merged)
                    loc_chr = 'M';
                else
                    loc_chr = 'L';
            }
            // pick near middle of plot window
            int duration = 2 * (time_max - deData->t_time_t);
            duration += 60; // try to make sure latest data is displayed
            if (duration > 3600)
                duration = 3600;
            double start_time = (double) (deData->t_time_t) + deData->t_decsec - (double) duration / 2.0;
            fprintf(pickMessageHtmlStream, "<td>%s</td><td>%s</td><td>%c</td>",
                    create_channel_links(deData->network, deData->station, deData->location, deData->channel,
                    deData->pick_stream, pick_stream_name(deData), deData->n_int_tseries, start_time, duration,
                    tmp_str, tmp_str_2), tmp_str_2, loc_chr);
            fprintf(picksCsvStream, "%s_%s_%s_%s %s %c ",
                    deData->network, deData->station, deData->location, deData->channel, pick_stream_name(deData), loc_chr);
            // set s/n ratios
            double snr_hf = deData->a_ref < 0.0 || deData->sn_pick < FLT_MIN ? 0.0 : deData->a_ref / deData->sn_pick;
            double snr_brb = deData->sn_brb_signal < 0.0 || deData->sn_brb_pick < FLT_MIN ? 0.0 : deData->sn_brb_signal / deData->sn_brb_pick;
            double snr_brb_int = deData->sn_brb_int_signal < 0.0 || deData->sn_brb_int_pick < FLT_MIN ? 0.0 : deData->sn_brb_int_signal / deData->sn_brb_int_pick;
            // pick time, error, polarity
            timeDecSec2string((double) deData->t_time_t + deData->t_decsec, tmp_str, DEFAULT_TIME_FORMAT);
            // 20121019 AJL - added pick polarityWeight
            double fmquality = 0.0;
            int fmpolarity = POLARITY_UNKNOWN;
            char fmtype[32] = "Err";
            setPolarity(deData, &fmquality, &fmpolarity, fmtype);
            //
            double grd_vel_peak_amp = GRD_MOT_INVALID;
            double grd_disp_peak_amp = GRD_MOT_INVALID;
            if (is_associated_location_P(deData) && flag_do_grd_vel && is_P(deData->phase_id)) {
                // initialize and calculate brb hp ground motion peaks after pick
                // 20140801 AJL - added
                deData->grd_mot->peak_amp_vel = GRD_MOT_INVALID;
                deData->grd_mot->peak_amp_disp = GRD_MOT_INVALID;
                if (!deData->flag_snr_brb_too_low || !deData->flag_snr_brb_int_too_low) {
                    calculate_init_P_grd_mot_amp(deData, snr_brb, snr_brb_int, 0, &fmquality, &fmpolarity, fmtype);
                    if (!deData->flag_snr_brb_too_low) {
                        grd_vel_peak_amp = deData->grd_mot->peak_amp_vel;
                    }
                    if (!deData->flag_snr_brb_int_too_low) {
                        grd_disp_peak_amp = deData->grd_mot->peak_amp_disp;
                    }
                }
            }
            fprintf(pickMessageHtmlStream, "<td>%s</td><td>%.3f</td><td>%d</td><td>%s&nbsp;</td><td>%.2f</td><td>%.1f</td>",
                    tmp_str, deData->pick_error, fmpolarity, fmtype, fmquality, deData->take_off_angle_inc);
            fprintf(picksCsvStream, "%s %.3f %d %s %.2f %.1f ",
                    tmp_str, deData->pick_error, fmpolarity, fmtype, fmquality, deData->take_off_angle_inc);
            // waveform onset polarization azimuth (e.g. P polarization azimuth)
            fprintf(pickMessageHtmlStream, "<td>%.0f</td><td>%.0f</td><td>%.0f</td>",
                    deData->polarization.azimuth, deData->polarization.azimuth_unc, deData->polarization.azimuth_calc);
            fprintf(picksCsvStream, "%.1f %.1f %.1f ",
                    deData->polarization.azimuth, deData->polarization.azimuth_unc, deData->polarization.azimuth_calc);
            // phase association
            if (associate_data) {
                if (use_associated_data && is_associated_phase(deData)) {
                    if (deData->polarization.weight >= 0.0) {
                        fprintf(pickMessageHtmlStream, "<td>&nbsp;%.2f</td>", deData->polarization.weight);
                    } else {
                        fprintf(pickMessageHtmlStream, "<td>-1</td>");
                    }
                    fprintf(pickMessageHtmlStream, "<td>%s</td>", deData->phase);
                    fprintf(pickMessageHtmlStream, "<td>%.1f</td>", deData->residual);
                    if (deData->loc_weight > 0.001) {
                        fprintf(pickMessageHtmlStream, "<td>%.2f</td>", deData->loc_weight);
                    } else {
                        fprintf(pickMessageHtmlStream, "<td>0</td>");
                    }
                    fprintf(pickMessageHtmlStream, "<td>%.2f</td>", deData->dist_weight);
                    fprintf(pickMessageHtmlStream, "<td>%.2f</td>", deData->station_quality_weight);
                    fprintf(picksCsvStream, "%.2f %s %.2f %.3f %.3f %.3f ", deData->polarization.weight, deData->phase, deData->residual, deData->loc_weight, deData->dist_weight, deData->station_quality_weight);
                } else {
                    fprintf(pickMessageHtmlStream, "<td></td><td></td><td></td><td></td><td></td><td>%.2f</td>", deData->station_quality_weight);
                    fprintf(picksCsvStream, "-1 -1 -1 -1 -1 %.3f ", deData->station_quality_weight);
                }
            }
            double t50_a_ref_have_gain_flag = 1.0; // t50 and a_ref values are positive if not using amplitude attenuation
            // 20140120 AJL - test of amplitude attenuation with distance
            if (!(chan_resp[deData->source_id].have_gain && chan_resp[deData->source_id].responseType == DERIVATIVE_TYPE)) {
                t50_a_ref_have_gain_flag = -1.0; // t50 and a_ref values are negative if not corrected for gain
            }
            fprintf(pickMessageHtmlStream, "<td>%.3g</td><td>%.3g</td>",
                    deData->t50 * t50_a_ref_have_gain_flag, deData->a_ref * t50_a_ref_have_gain_flag);
            if (deData->amplitude_error_ratio > 0.0
                    && deData->amplitude_error_ratio >= AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO
                    && deData->amplitude_error_ratio <= AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO) {
                fprintf(pickMessageHtmlStream, "<td>%.3g</td>", deData->amplitude_error_ratio);
            } else {
                fprintf(pickMessageHtmlStream, "<td bgcolor=\"#EEEEEE\"><em>%.3g</em></td>", deData->amplitude_error_ratio);
            }
            fprintf(pickMessageHtmlStream, "<td %s>%s</td><td %s>%s</td>",
                    t50ExBackgroundColor, deLevelStr, taucBackgroundColor, tauc_peakStr);
            fprintf(pickMessageHtmlStream, "<td>%.3g</td><td>%.3g</td>",
                    grd_vel_peak_amp, grd_disp_peak_amp);
            fprintf(pickMessageHtmlStream, "<td>%.2f</td><td>%.2f</td><td>%.2f</td>",
                    snr_hf, snr_brb, snr_brb_int);
            if (mb_mag != MB_INVALID)
                fprintf(pickMessageHtmlStream, "<td %s>%.1f</td>", mbBackgroundColor, mb_mag);
            else
                fprintf(pickMessageHtmlStream, "<td %s>-9</td>", mbBackgroundColor);
            if (mwp_mag != MWP_INVALID)
                fprintf(pickMessageHtmlStream, "<td %s>%.1f</td>", mwpBackgroundColor, mwp_mag);
            else
                fprintf(pickMessageHtmlStream, "<td %s>-9</td>", mwpBackgroundColor);
            if (t0_dur != T0_INVALID)
                fprintf(pickMessageHtmlStream, "<td %s>%.0f</td>", t0BackgroundColor, t0_dur);
            else
                fprintf(pickMessageHtmlStream, "<td %s>-9</td>", t0BackgroundColor);
            if (mwpd_corr_mag != MWPD_INVALID)
                fprintf(pickMessageHtmlStream, "<td %s>%.1f</td>", mwpdBackgroundColor, mwpd_corr_mag);
            else
                fprintf(pickMessageHtmlStream, "<td %s>-9</td>", mwpdBackgroundColor);
            fprintf(picksCsvStream, "%.3g %.3g %.3g %.3f %.3f ",
                    deData->t50 * t50_a_ref_have_gain_flag, deData->a_ref * t50_a_ref_have_gain_flag, deData->amplitude_error_ratio,
                    deLevel, tauc_peak);
            fprintf(picksCsvStream, "%.3g %.3g ",
                    grd_vel_peak_amp, grd_disp_peak_amp);
            fprintf(picksCsvStream, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f ",
                    deData->a_ref / deData->sn_pick, snr_brb, snr_brb_int, mb_mag, mwp_mag, t0_dur, mwpd_corr_mag);
            fprintf(pickMessageHtmlStream, "<td>");
            if (deData->flag_clipped) {
                fprintf(pickMessageHtmlStream, "%sCLIPPED", ignored ? "Ignored:" : "Loc:");
                fprintf(picksCsvStream, "%sCLIPPED", ignored ? "Ignored:" : "Loc:");
            } else if (deData->flag_non_contiguous) {
                fprintf(pickMessageHtmlStream, "%sNON_CONTIG", ignored ? "Ignored:" : "Loc:");
                fprintf(picksCsvStream, "%sCLIPPED", ignored ? "Ignored:" : "Loc:");
            } else if (deData->flag_a_ref_not_ok) {
                fprintf(pickMessageHtmlStream, "%sin_prev_coda", ignored ? "Ignored:" : "Loc:");
                fprintf(picksCsvStream, "%sin_prev_coda", ignored ? "Ignored:" : "Loc:");
            } else if (!deData->flag_complete_t50) {
                fprintf(pickMessageHtmlStream, "incomplete...");
                fprintf(picksCsvStream, "incomplete");
                // 20130128 AJL - use flag_snr_brb_int_too_low to allow mwp, mwpd, etc., but do not use for ignore tests (e.g. ignore determined by flag_snr_brb_too_low)
                //} else if (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low && deData->flag_snr_brb_int_too_low) {
            } else if (deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low) {
                fprintf(pickMessageHtmlStream, "%ss/n_HF_BRBV_low", ignored ? "Ignored:" : "Loc:");
                fprintf(picksCsvStream, "%ss/n_HF_BRB_low", ignored ? "Ignored:" : "Loc:");
            } else if (deData->flag_snr_hf_too_low) {
                fprintf(pickMessageHtmlStream, "s/n_HF_low");
                fprintf(picksCsvStream, "s/n_HF_low");
            } else if (deData->flag_snr_brb_too_low) {
                fprintf(pickMessageHtmlStream, "s/n_BRBV_low");
                fprintf(picksCsvStream, "s/n_BRBV_low");
            } else if (deData->flag_snr_brb_int_too_low) {
                fprintf(pickMessageHtmlStream, "s/n_BRBD_low");
                fprintf(picksCsvStream, "s/n_BRBD_low");
            } else {
                fprintf(pickMessageHtmlStream, "OK_HF_BRB");
                fprintf(picksCsvStream, "OK_HF_BRB");
            }
            fprintf(pickMessageHtmlStream, "</td>");
            fprintf(pickMessageHtmlStream, "</tr>\n");
            // 20150716 AJL - add station correction to pick csv file
            fprintf(picksCsvStream, " %.3f ", deData->sta_corr);
            fprintf(picksCsvStream, "\n");
        }
    }
    fclose_counter(pickDataNLLStream);
    fprintf(pickMessageHtmlStream, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose_counter(pickMessageHtmlStream);
    fclose_counter(picksCsvStream);

    // write pick summary log
    // pickMessageLogStream
    sprintf(outname, "%s/picks.log", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * pickMessageLogStream = fopen_counter(outname, "w");
    if (pickMessageLogStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    fprintf(pickMessageLogStream, "%ld %d %d %d %d %d %d %d %d %d %d %d\n",
            plot_time_max, num_pick, num_pick_not_completed, num_pick_clipped, num_pick_non_contiguous, num_pick_sn_low, num_pick_sn_brb_low, num_pick_sn_brb_int_low,
            num_pick_in_prev_coda, num_pick_ok, num_pick_assoc, num_event);
    fclose_counter(pickMessageLogStream);


    // miscellaneous output ===============================================================

    // GMT style
    sprintf(outname, "%s/plot/sta.healthy.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staHealthyStream = fopen_counter(outname, "w");
    if (staHealthyStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/sta.requested.xy", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staRequestedStream = fopen_counter(outname, "w");
    if (staRequestedStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // GMT style
    sprintf(outname, "%s/plot/sta.code.txt", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * staCodeStream = fopen_counter(outname, "w");
    if (staCodeStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // write all recently active stations to sta streams
    //int nstaIsActive = 0;
    //int nstaHasBeenActive = 0;
    for (n = 0; n < num_sources_total; n++) {
        if (channelParameters[n].have_coords && !channelParameters[n].inactive_duplicate && channelParameters[n].process_this_channel_orientation) {
            fprintf(staCodeStream, "%f %f %s\n", channelParameters[n].lat, channelParameters[n].lon, channelParameters[n].station);
            //nstaHasBeenActive++;
            //if (channelParameters[n].data_latency < report_interval)
            //    nstaIsActive++;
            if ((channelParameters[n].data_latency < LATENCY_YELLOW_CUTOFF) && (channelParameters[n].qualityWeight > DATA_UNASSOC_WT_YELLOW_CUTOFF) && !(channelParameters[n].error)) {
                fprintf(staHealthyStream, "%f %f\n", channelParameters[n].lat, channelParameters[n].lon);
            } else {
                fprintf(staRequestedStream, "%f %f\n", channelParameters[n].lat, channelParameters[n].lon);
            }
        }
    }
    //
    fclose_counter(staHealthyStream);
    fclose_counter(staRequestedStream);
    fclose_counter(staCodeStream);



    // write alarm csv data ===============================================================

    double level_plot = 0.0;

    // T50 Level
    sprintf(outname, "%s/t50.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * t50ExCsvStream = fopen_counter(outname, "w");
    if (t50ExCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // t50 10.0
    fprintf(t50ExCsvStream, "t50");
    fprintf(t50ExCsvStream, " %lf", T50EX_CRITICAL_VALUE);
    fprintf(t50ExCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(t50ExCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(t50ExCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(t50ExCsvStream, " %d", report_interval);
    fprintf(t50ExCsvStream, " %d", nstaIsActive);
    fprintf(t50ExCsvStream, " %d", nstaHasBeenActive);
    fprintf(t50ExCsvStream, " %d", num_hypocenters_associated);
    fprintf(t50ExCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = t50ExLevelStatistics[nhyp].centralValue < T50EX_LEVEL_MAX ? t50ExLevelStatistics[nhyp].centralValue : T50EX_LEVEL_MAX;
        level_plot = level_plot > T50EX_LEVEL_MIN ? level_plot : T50EX_LEVEL_MIN;
        fprintf(t50ExCsvStream, " %.1f", level_plot);
        fprintf(t50ExCsvStream, " %.1f", t50ExLevelStatistics[nhyp].centralValue);
        fprintf(t50ExCsvStream, " %.1f", t50ExLevelStatistics[nhyp].lowerBound);
        fprintf(t50ExCsvStream, " %.1f", t50ExLevelStatistics[nhyp].upperBound);
        fprintf(t50ExCsvStream, " %s", t50ExLevelString[nhyp]);
        fprintf(t50ExCsvStream, " %d", t50ExLevelStatistics[nhyp].numLevel);
        fprintf(t50ExCsvStream, "\n");
    }
    fclose_counter(t50ExCsvStream);
    //
    FILE* taucCsvStream = NULL;
    sprintf(outname, "%s/tauc.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    taucCsvStream = fopen_counter(outname, "w");
    if (taucCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // tauc 10.0
    fprintf(taucCsvStream, "tauc");
    fprintf(taucCsvStream, " %lf", TAUC_CRITICAL_VALUE);
    fprintf(taucCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(taucCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(taucCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(taucCsvStream, " %d", report_interval);
    fprintf(taucCsvStream, " %d", nstaIsActive);
    fprintf(taucCsvStream, " %d", nstaHasBeenActive);
    fprintf(taucCsvStream, " %d", num_hypocenters_associated);
    fprintf(taucCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = taucLevelStatistics[nhyp].centralValue < TAUC_LEVEL_MAX ? taucLevelStatistics[nhyp].centralValue : TAUC_LEVEL_MAX;
        level_plot = level_plot > TAUC_LEVEL_MIN ? level_plot : TAUC_LEVEL_MIN;
        fprintf(taucCsvStream, " %.1f", level_plot);
        fprintf(taucCsvStream, " %.1f", taucLevelStatistics[nhyp].centralValue);
        fprintf(taucCsvStream, " %.1f", taucLevelStatistics[nhyp].lowerBound);
        fprintf(taucCsvStream, " %.1f", taucLevelStatistics[nhyp].upperBound);
        fprintf(taucCsvStream, " %s", taucLevelString[nhyp]);
        fprintf(taucCsvStream, " %d", taucLevelStatistics[nhyp].numLevel);
        fprintf(taucCsvStream, "\n");
    }
    fclose_counter(taucCsvStream);
    //
    FILE* tdT50ExCsvStream = NULL;
    // TdT50Ex
    sprintf(outname, "%s/alarm.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    tdT50ExCsvStream = fopen_counter(outname, "w");
    if (tdT50ExCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // alarm 10.0
    fprintf(tdT50ExCsvStream, "alarm");
    fprintf(tdT50ExCsvStream, " %lf", TDT50EX_CRITICAL_VALUE);
    fprintf(tdT50ExCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(tdT50ExCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(tdT50ExCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(tdT50ExCsvStream, " %d", report_interval);
    fprintf(tdT50ExCsvStream, " %d", nstaIsActive);
    fprintf(tdT50ExCsvStream, " %d", nstaHasBeenActive);
    fprintf(tdT50ExCsvStream, " %d", num_hypocenters_associated);
    fprintf(tdT50ExCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = tdT50ExLevelStatistics[nhyp].centralValue < TDT50EX_LEVEL_MAX ? tdT50ExLevelStatistics[nhyp].centralValue : TDT50EX_LEVEL_MAX;
        level_plot = level_plot > TDT50EX_LEVEL_MIN ? level_plot : TDT50EX_LEVEL_MIN;
        fprintf(tdT50ExCsvStream, " %.1f", level_plot);
        fprintf(tdT50ExCsvStream, " %.1f", tdT50ExLevelStatistics[nhyp].centralValue);
        fprintf(tdT50ExCsvStream, " %.1f", tdT50ExLevelStatistics[nhyp].lowerBound);
        fprintf(tdT50ExCsvStream, " %.1f", tdT50ExLevelStatistics[nhyp].upperBound);
        fprintf(tdT50ExCsvStream, " %s", warningLevelString[nhyp]);
        fprintf(tdT50ExCsvStream, " %d", tdT50ExLevelStatistics[nhyp].numLevel);
        fprintf(tdT50ExCsvStream, "\n");
    }
    fclose_counter(tdT50ExCsvStream);
    //
    FILE* mwpCsvStream = NULL;
    sprintf(outname, "%s/mwp.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    mwpCsvStream = fopen_counter(outname, "w");
    if (mwpCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // mwp 10.0
    fprintf(mwpCsvStream, "mwp");
    fprintf(mwpCsvStream, " %lf", MWP_CRITICAL_VALUE);
    fprintf(mwpCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(mwpCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(mwpCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(mwpCsvStream, " %d", report_interval);
    fprintf(mwpCsvStream, " %d", nstaIsActive);
    fprintf(mwpCsvStream, " %d", nstaHasBeenActive);
    fprintf(mwpCsvStream, " %d", num_hypocenters_associated);
    fprintf(mwpCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = mwpLevelStatistics[nhyp].centralValue < MWP_LEVEL_MAX ? mwpLevelStatistics[nhyp].centralValue : MWP_LEVEL_MAX;
        level_plot = level_plot > MWP_LEVEL_MIN ? level_plot : MWP_LEVEL_MIN;
        fprintf(mwpCsvStream, " %.1f", level_plot);
        fprintf(mwpCsvStream, " %.1f", mwpLevelStatistics[nhyp].centralValue);
        fprintf(mwpCsvStream, " %.1f", mwpLevelStatistics[nhyp].lowerBound);
        fprintf(mwpCsvStream, " %.1f", mwpLevelStatistics[nhyp].upperBound);
        fprintf(mwpCsvStream, " %s", mwpLevelString[nhyp]);
        fprintf(mwpCsvStream, " %d", mwpLevelStatistics[nhyp].numLevel);
        fprintf(mwpCsvStream, "\n");
    }
    fclose_counter(mwpCsvStream);
    //
    FILE* mbCsvStream = NULL;
    sprintf(outname, "%s/mb.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    mbCsvStream = fopen_counter(outname, "w");
    if (mbCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // mb 10.0
    fprintf(mbCsvStream, "mb");
    fprintf(mbCsvStream, " %lf", MB_CRITICAL_VALUE);
    fprintf(mbCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(mbCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(mbCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(mbCsvStream, " %d", report_interval);
    fprintf(mbCsvStream, " %d", nstaIsActive);
    fprintf(mbCsvStream, " %d", nstaHasBeenActive);
    fprintf(mbCsvStream, " %d", num_hypocenters_associated);
    fprintf(mbCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = mbLevelStatistics[nhyp].centralValue < MB_LEVEL_MAX ? mbLevelStatistics[nhyp].centralValue : MB_LEVEL_MAX;
        level_plot = level_plot > MB_LEVEL_MIN ? level_plot : MB_LEVEL_MIN;
        fprintf(mbCsvStream, " %.1f", level_plot);
        fprintf(mbCsvStream, " %.1f", mbLevelStatistics[nhyp].centralValue);
        fprintf(mbCsvStream, " %.1f", mbLevelStatistics[nhyp].lowerBound);
        fprintf(mbCsvStream, " %.1f", mbLevelStatistics[nhyp].upperBound);
        fprintf(mbCsvStream, " %s", mbLevelString[nhyp]);
        fprintf(mbCsvStream, " %d", mbLevelStatistics[nhyp].numLevel);
        fprintf(mbCsvStream, "\n");
    }
    fclose_counter(mbCsvStream);
    //
    FILE* t0CsvStream = NULL;
    sprintf(outname, "%s/t0.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    t0CsvStream = fopen_counter(outname, "w");
    if (t0CsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // t0 10.0
    fprintf(t0CsvStream, "t0");
    fprintf(t0CsvStream, " %lf", T0_CRITICAL_VALUE);
    fprintf(t0CsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(t0CsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(t0CsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(t0CsvStream, " %d", report_interval);
    fprintf(t0CsvStream, " %d", nstaIsActive);
    fprintf(t0CsvStream, " %d", nstaHasBeenActive);
    fprintf(t0CsvStream, " %d", num_hypocenters_associated);
    fprintf(t0CsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = t0LevelStatistics[nhyp].centralValue < T0_LEVEL_MAX ? t0LevelStatistics[nhyp].centralValue : T0_LEVEL_MAX;
        level_plot = level_plot > T0_LEVEL_MIN ? level_plot : T0_LEVEL_MIN;
        fprintf(t0CsvStream, " %.1f", level_plot);
        fprintf(t0CsvStream, " %.0f", t0LevelStatistics[nhyp].centralValue);
        fprintf(t0CsvStream, " %.0f", t0LevelStatistics[nhyp].lowerBound);
        fprintf(t0CsvStream, " %.0f", t0LevelStatistics[nhyp].upperBound);
        fprintf(t0CsvStream, " %s", t0LevelString[nhyp]);
        fprintf(t0CsvStream, " %d", t0LevelStatistics[nhyp].numLevel);
        fprintf(t0CsvStream, "\n");
    }
    fclose_counter(t0CsvStream);
    //
    FILE* mwpdCsvStream = NULL;
    sprintf(outname, "%s/mwpd.csv", outnameroot);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    mwpdCsvStream = fopen_counter(outname, "w");
    if (mwpdCsvStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    // mwpd 10.0
    fprintf(mwpdCsvStream, "mwpd");
    fprintf(mwpdCsvStream, " %lf", MWPD_CRITICAL_VALUE);
    fprintf(mwpdCsvStream, " %lf", difftime(plot_time_max, plot_time_min) / 60.0);
    fprintf(mwpdCsvStream, " %s", time2string(plot_time_min, tmp_str));
    fprintf(mwpdCsvStream, " %s", time2string(plot_time_max, tmp_str));
    fprintf(mwpdCsvStream, " %d", report_interval);
    fprintf(mwpdCsvStream, " %d", nstaIsActive);
    fprintf(mwpdCsvStream, " %d", nstaHasBeenActive);
    fprintf(mwpdCsvStream, " %d", num_hypocenters_associated);
    fprintf(mwpdCsvStream, "\n");
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        level_plot = mwpdLevelStatistics[nhyp].centralValue < MWPD_LEVEL_MAX ? mwpdLevelStatistics[nhyp].centralValue : MWPD_LEVEL_MAX;
        level_plot = level_plot > MWPD_LEVEL_MIN ? level_plot : MWPD_LEVEL_MIN;
        fprintf(mwpdCsvStream, " %.1f", level_plot);
        fprintf(mwpdCsvStream, " %.1f", mwpdLevelStatistics[nhyp].centralValue);
        fprintf(mwpdCsvStream, " %.1f", mwpdLevelStatistics[nhyp].lowerBound);
        fprintf(mwpdCsvStream, " %.1f", mwpdLevelStatistics[nhyp].upperBound);
        fprintf(mwpdCsvStream, " %s", mwpdLevelString[nhyp]);
        fprintf(mwpdCsvStream, " %d", mwpdLevelStatistics[nhyp].numLevel);
        fprintf(mwpdCsvStream, "\n");
    }
    fclose_counter(mwpdCsvStream);


    // perform update of waveform export for associated P ===============================================================

    if (waveform_export_enable) {
        // initialize waveform export flags for each hypo
        static int waveform_exported[MAX_NUM_HYPO];
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            waveform_exported[nhyp] = 0;
        }
        // set waveform files root path
        static char waveforms_root[STANDARD_STRLEN];
        sprintf(waveforms_root, "%s/waveforms/", outnameroot_archive); // directory
        mkdir(waveforms_root, 0755);
        static char hdr_filename[STANDARD_STRLEN];
        // update waveform export for all associated P data
        //int i = 0;
        //printf("TP %d\n", i++);
        for (n = num_de_data - 1; n >= 0; n--) {
            //i = 1;
            //printf("TP %d num_de_data=%d\n", i++, n);
            TimedomainProcessingData* deData = data_list[n];
            hptime_t start_time_written, end_time_written;
            // 20140417 AJL - make selection less strict to include all first arrival or direct P phases which can be counted in location
            //if (is_associated_location_P(deData) && is_count_in_location(deData->phase_id)) {
            if (is_associated_phase(deData) && (is_first_arrival_P(deData->phase_id) || is_direct_P(deData->phase_id)) && is_count_in_location(deData->phase_id)) {
                //printf("TP %d\n", i++);
                // int waveform_export_memory_sliding_window_length;   // seconds
                // int waveform_export_window_start_before_P;   // seconds
                // int waveform_export_window_end_after_S;   // seconds
                //printf("TP %d waveform_export_window_start_before_P=%ld, waveform_export_window_end_after_S=%ld\n", i++,
                //waveform_export_window_start_before_P, waveform_export_window_end_after_S);
                HypocenterDesc* phypo = hyp_assoc_loc[deData->is_associated - 1];
                hptime_t origin_time = timeDecSec2hptime(phypo->otime);
                if (deData->ttime_P == VALUE_NOT_SET)
                    printf("ERROR: waveform export: %s_%s_%s_%s %s: deData->ttime_P not set: this should not happen!\n",
                        deData->network, deData->station, deData->location, deData->channel, deData->phase);
                hptime_t start_time = origin_time + (deData->ttime_P - waveform_export_window_start_before_P) * (double) HPTMODULUS;
                if (deData->ttime_S == VALUE_NOT_SET)
                    printf("ERROR: waveform export: %s_%s_%s_%s %s: deData->ttime_S not set: this should not happen!\n",
                        deData->network, deData->station, deData->location, deData->channel, deData->phase);
                hptime_t end_time = origin_time + (deData->ttime_S + waveform_export_window_end_after_S) * (double) HPTMODULUS;
                //printf("TP %d ttime_P=%f, ttime_S=%f\n", i++, deData->ttime_P, deData->ttime_S);
                //printf("TP %d start_time=%ld, origin_time=%ld, end_time=%ld\n", i++, start_time/HPTMODULUS, origin_time/HPTMODULUS, end_time/HPTMODULUS);
                //printf("TP %d start_time-origin_time=%ld, end_time-origin_time=%ld\n", i++, (start_time-origin_time)/HPTMODULUS, (end_time - origin_time)/HPTMODULUS);
                //printf("TP %d mslist_getStartTime=%ld, mslist_getEndTime=%ld\n", i++,
                //mslist_getStartTime(waveform_export_miniseed_list[deData->source_id], num_waveform_export_miniseed_list[deData->source_id])/HPTMODULUS,
                //mslist_getEndTime(waveform_export_miniseed_list[deData->source_id], num_waveform_export_miniseed_list[deData->source_id])/HPTMODULUS);

                // 20160808 AJL - add support for 3-comp write
                int source_id_write = deData->source_id;
                ChannelParameters * chan_params = channelParameters + source_id_write;
                for (int ncomp = 0; ncomp < 3; ncomp++) {
                    ChannelParameters* chan_params_write = chan_params; // ncomp = 0
                    if (ncomp > 0) { // check for other channel orientations
                        if (chan_params->channel_set[ncomp - 1] >= 0) {
                            source_id_write = chan_params->channel_set[ncomp - 1];
                            chan_params_write = channelParameters + source_id_write;
                        } else {
                            continue;
                        }
                    }
                    if (
                            (
                            deData->waveform_export[ncomp].start_time_written < 0 || // data not yet written
                            (start_time < deData->waveform_export[ncomp].start_time_written // required start time earlier than start time of written data
                            && deData->waveform_export[ncomp].start_time_written > // AND more data available at beginning (impossible with real-time SeedLink data)
                            mslist_getStartTime(waveform_export_miniseed_list[source_id_write], num_waveform_export_miniseed_list[source_id_write]))
                            )
                            ||
                            (
                            deData->waveform_export[ncomp].end_time_written < 0 || // data not yet written
                            (end_time > deData->waveform_export[ncomp].end_time_written // required end time later than end time of written data
                            && deData->waveform_export[ncomp].end_time_written < // AND more data available at end
                            mslist_getEndTime(waveform_export_miniseed_list[source_id_write], num_waveform_export_miniseed_list[source_id_write]))
                            )
                            ) {
                        sprintf(outname, "%s/%ld/", waveforms_root, phypo->unique_id); // directory
                        mkdir(outname, 0755);
                        char* filename = deData->waveform_export[ncomp].filename;
                        if (filename[0] == '\0' || deData->waveform_export[ncomp].hypo_unique_id != phypo->unique_id) {
                            strcpy(filename, outname);
                            strcat(filename, chan_params_write->network);
                            strcat(filename, ".");
                            strcat(filename, chan_params_write->station);
                            strcat(filename, ".");
                            strcat(filename, chan_params_write->location);
                            strcat(filename, ".");
                            strcat(filename, chan_params_write->channel);
                            strcat(filename, ".mseed");
                        }
                        //printf("TP %d filename=%s\n", i++, filename);
                        int sampleLength = mslist_writeToFile(filename, start_time, end_time,
                                waveform_export_miniseed_list[source_id_write], num_waveform_export_miniseed_list[source_id_write], verbose,
                                &start_time_written, &end_time_written);
                        //printf("DEBUG: mslist_writeToFile: nwritten %d: %s time:%ld->%ld twritten:%ld->%ld diff:%ld->%ld\n",
                        //        nrecords_written, filename, start_time, end_time, start_time_written, end_time_written, start_time_written - start_time, end_time_written - end_time);
                        deData->waveform_export[ncomp].start_time_written = start_time_written;
                        deData->waveform_export[ncomp].end_time_written = end_time_written;
                        deData->waveform_export[ncomp].hypo_unique_id = phypo->unique_id;
                        waveform_exported[deData->is_associated - 1] = 1;
                        // write waveform header file
                        strcpy(hdr_filename, deData->waveform_export[ncomp].filename);
                        strcat(hdr_filename, ".sg2k");
                        sprintf(tmp_str, "%ld", phypo->unique_id);
                        int iyear, ijday, ihour, imin;
                        double sec;
                        double comp_azimuth = chan_params_write->azimuth;
                        // SAC / SG2K: inclination = Component incident angle (degrees from vertical up).   eg.    0, 90, 180
                        // FDSN/SEED: dip = Dip of the instrument in degrees, down from horizontal          e.g. -90,  0,  90
                        double comp_inclination = chan_params_write->dip + 90.0;
                        double baz = GCAzimuth(chan_params_write->lat, chan_params_write->lon, phypo->lat, phypo->lon);
                        char component[4];
                        if (strlen(chan_params_write->channel) >= 3) {
                            strcpy(component, chan_params_write->channel + 2);
                        } else {
                            strcpy(component, chan_params_write->channel);
                        }
                        char* loc;
                        if (strcmp(chan_params_write->location, "--") == 0) {
                            loc = NULL;
                        } else {
                            loc = chan_params_write->location;
                        }
                        // 20140123 AJL
                        //double gain = chan_resp[source_id_write].have_gain ? chan_resp[source_id_write].gain : -1.0;
                        double gain = chan_resp[source_id_write].have_gain && chan_resp[source_id_write].responseType == DERIVATIVE_TYPE
                                ? chan_resp[source_id_write].gain : -1.0;
                        //printf("DEBUG: %s_%s_%s have_gain=%d gain=%e->%e \n", chan_params_write->network, chan_params_write->station, chan_params_write->channel,
                        //        chan_resp[source_id_write].have_gain, chan_resp[source_id_write].gain, gain);
                        hptime2dateTimeComponents(start_time_written, &iyear, &ijday, &ihour, &imin, &sec);
                        writeSG2Kheader(hdr_filename, tmp_str,
                                iyear, ijday, ihour, imin, sec, 0.0,
                                chan_params_write->network, chan_params_write->station,
                                NULL, chan_params_write->channel, loc, chan_params_write->channel,
                                comp_azimuth, comp_inclination,
                                chan_params_write->lat, chan_params_write->lon,
                                0.0, chan_params_write->elev / 1000.0, gain,
                                sampleLength, 1.0 / chan_params_write->deltaTime, "counts", "sec",
                                phypo->lat, phypo->lon, phypo->depth,
                                timeDecSec2string(phypo->otime, tmp_str_2, COMMA_DELIMTED_TIME_FORMAT),
                                -9.9, phypo->mbLevelStatistics.centralValue, phypo->mwpLevelStatistics.centralValue, -9.9, -9.9, phypo->mwpLevelStatistics.centralValue,
                                deData->epicentral_distance, deData->epicentral_distance * DEG2KM, deData->epicentral_azimuth, baz,
                                1
                                );
                    }
                }
            }
        }
        // clean up old waveform data
        if (remove_waveform_export_directories(waveforms_root) != 0) {
            printf("ERROR: waveform export: while cleaning up old waveform data.\n");
        }
        // write hypocenter info and picks for each associated hypocenter for which waveforms were exported
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            if (waveform_exported[nhyp]) {
                HypocenterDesc* phypo = hyp_assoc_loc[nhyp];
                // hypocenter  csv strings
                sprintf(outname, "%s/%ld/%s", waveforms_root, phypo->unique_id, "hypo.csv");
                FILE * fpout = fopen_counter(outname, "w");
                if (fpout == NULL) {
                    sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                    perror(tmp_str);
                    continue;
                }
                printHypoDataHeaderString(hypoDataString);
                fprintf(fpout, "%s\n", hypoDataString);
                printHypoDataString(phypo, hypoDataString, 1);
                fprintf(fpout, "%s\n", hypoDataString);
                fclose_counter(fpout);
                sprintf(outname, "%s/%ld/%s", waveforms_root, phypo->unique_id, "hypo_pretty.csv");
                fpout = fopen_counter(outname, "w");
                if (fpout == NULL) {
                    sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                    perror(tmp_str);
                    continue;
                }
                printHypoDataHeaderString(hypoDataString);
                fprintf(fpout, "%s\n", hypoDataString);
                printHypoDataString(phypo, hypoDataString, 0);
                fprintf(fpout, "%s\n", hypoDataString);
                fclose_counter(fpout);
                // picks
                sprintf(outname, "%s/%ld/%s", waveforms_root, phypo->unique_id, "pickdata_nlloc.txt");
                fpout = fopen_counter(outname, "w");
                if (fpout == NULL) {
                    sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                    perror(tmp_str);
                    continue;
                }
                for (n = 0; n < num_de_data; n++) {
                    TimedomainProcessingData* deData = data_list[n];
                    if (deData->is_associated - 1 == nhyp) {
                        fprintf_NLLoc_TimedomainProcessingData(deData, fpout, 0);
                        fprintf(fpout, " %ld\n", phypo->unique_id);
                    }
                }
                fclose_counter(fpout);
            }
        }
    }


    // perform update of waveform export for all data if requested (file "export_waveforms_all" exists in working directory) ====================

    if (waveform_export_enable) {

        // if file "export_waveforms_all.flag" exists in working directory, then all waveforms for past hour will be exported

        int export_flag = 0;
        sprintf(outname, "%s", "export_waveforms_all");
        int fdtest = open(outname, O_RDONLY);
        if (fdtest > 0) {
            // the file exists
            export_flag = 1;
            close(fdtest);
        }
        printf("Info: export_waveforms_all: flag_file_name:%s fdtest:%d export_flag:%d\n", outname, fdtest, export_flag);
        remove(outname);

        if (export_flag) {

            // set waveform files root path
            static char waveforms_root[STANDARD_STRLEN];
            sprintf(waveforms_root, "%s/waveforms_all/", outnameroot_archive); // directory
            mkdir(waveforms_root, 0755);
            static char hdr_filename[STANDARD_STRLEN];
            // waveform export for all data
            //int i = 0;
            //printf("TP %d\n", i++);
            for (n = num_de_data - 1; n >= 0; n--) {
                TimedomainProcessingData* deData = data_list[n];
                hptime_t start_time = time_min * (double) HPTMODULUS;
                hptime_t end_time = time_max * (double) HPTMODULUS;
                sprintf(outname, "%s/%ld/", waveforms_root, time_min); // directory
                mkdir(outname, 0755);
                strcat(outname, deData->network);
                strcat(outname, ".");
                strcat(outname, deData->station);
                strcat(outname, ".");
                strcat(outname, deData->location);
                strcat(outname, ".");
                strcat(outname, deData->channel);
                strcat(outname, ".mseed");
                hptime_t start_time_written, end_time_written;
                int sampleLength = mslist_writeToFile(outname, start_time, end_time,
                        waveform_export_miniseed_list[deData->source_id], num_waveform_export_miniseed_list[deData->source_id], verbose,
                        &start_time_written, &end_time_written);
                printf("Info: export_waveforms_all: mslist_writeToFile: nsamp_written %d: %s %s->%s time:%ld->%ld twritten:%ld->%ld diff:%ld->%ld\n",
                        sampleLength, outname, asctime(gmtime(&time_min)), asctime(gmtime(&time_max)),
                        (long) start_time, (long) end_time, (long) start_time_written, (long) end_time_written,
                        (long) (start_time_written - start_time), (long) (end_time_written - end_time));
                // write waveform header file
                strcpy(hdr_filename, outname);
                strcat(hdr_filename, ".sg2k");
                sprintf(tmp_str, "%ld", time_min);
                int iyear, ijday, ihour, imin;
                double sec;
                double comp_azimuth = -1.0;
                double comp_inclination = -1.0;
                char component[4];
                if (strlen(deData->channel) >= 3) {
                    strcpy(component, deData->channel + 2);
                } else {
                    strcpy(component, deData->channel);
                }
                char* loc;
                if (strcmp(deData->location, "--") == 0) {
                    loc = NULL;
                } else {
                    loc = deData->location;
                }
                // 20140123 AJL
                //double gain = chan_resp[deData->source_id].have_gain ? chan_resp[deData->source_id].gain : -1.0;
                double gain = chan_resp[deData->source_id].have_gain && chan_resp[deData->source_id].responseType == DERIVATIVE_TYPE
                        ? chan_resp[deData->source_id].gain : -1.0;
                //printf("DEBUG: %s_%s_%s have_gain=%d gain=%e->%e \n", deData->network, deData->station, deData->channel,
                //        chan_resp[deData->source_id].have_gain, chan_resp[deData->source_id].gain, gain);
                hptime2dateTimeComponents(start_time_written, &iyear, &ijday, &ihour, &imin, &sec);
                writeSG2Kheader(hdr_filename, tmp_str,
                        iyear, ijday, ihour, imin, sec, 0.0,
                        deData->network, deData->station,
                        NULL, deData->channel, loc, deData->channel,
                        comp_azimuth, comp_inclination,
                        channelParameters[deData->source_id].lat, channelParameters[deData->source_id].lon,
                        0.0, channelParameters[deData->source_id].elev / 1000.0, gain,
                        sampleLength, 1.0 / deData->deltaTime, "counts", "sec",
                        DBL_INVALID, DBL_INVALID, DBL_INVALID,
                        NULL,
                        DBL_INVALID, DBL_INVALID, DBL_INVALID, DBL_INVALID, DBL_INVALID, DBL_INVALID,
                        DBL_INVALID, DBL_INVALID, DBL_INVALID, DBL_INVALID,
                        0
                        );
            }
        }
    }


    // hypocenter management ===============================================================

    // hypo list archive stream
    sprintf(outname, "%s/hypolist.csv", outnameroot_archive);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypoListStream = fopen_counter(outname, "w");
    if (hypoListStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    printHypoDataHeaderString(hypoDataString);
    fprintf(hypoListStream, "%s\n", hypoDataString);
    // hypo list archive stream html
    sprintf(outname, "%s/hypolist.html", outnameroot_archive);
    if (verbose > 2)
        printf("Opening output file: %s\n", outname);
    FILE * hypoListHtmlStream = fopen_counter(outname, "w");
    if (hypoListHtmlStream == NULL) {
        sprintf(tmp_str, "ERROR: opening output file: %s", outname);
        perror(tmp_str);
        return (-1);
    }
    //double archive_days = (double) MAX_HYPO_ARCHIVE_WINDOW / (double) 3600 / (double) 24;
    int archive_days = MAX_HYPO_ARCHIVE_WINDOW / 3600 / 24;
    fprintf(hypoListHtmlStream, "<html>\n<head>\n<title>Latest Earthquakes - Past %d days - %s</title>\n", archive_days, EARLY_EST_MONITOR_NAME);
    fprintf(hypoListHtmlStream, "<meta http-equiv=\"refresh\" content=\"30\">\n</head>\n<body style=\"font-family:sans-serif;font-size:small\">\n");
    //fprintf(hypoListHtmlStream, "<table border=0 cellpadding=0 cellspacing=2 width=100%%>\n<tbody>\n");
    fprintf(hypoListHtmlStream, "<table border=0 cellpadding=1 frame=box rules=rows width=100%%>\n<tbody>\n");
    printHypoMessageHtmlHeaderString(hypoMessageHtmlString);
    fprintf(hypoListHtmlStream, "%s\n", hypoMessageHtmlString);
    // check for and process hypocenters with otime before archive or before analysis window
    for (nhyp = num_hypocenters - 1; nhyp >= 0; nhyp--) { // reverse time order, since may remove hypocenters from list.
        HypocenterDesc* phypo = hypo_list[nhyp];
        //printf("process: hypo.ot %f\n", phypo->otime);
        // check if hypocenter should be removed or placed in archive
        if ((time_t) phypo->otime <= time_min - MAX_HYPO_ARCHIVE_WINDOW) {
            // otime is before archive window, remove hypocenter
            if (verbose > 0) {
                printf("INFO: otime is before archive window, remove hypocenter: %s, phypo->otime %ld, time_min - MAX_HYPO_ARCHIVE_WINDOW %ld\n",
                        timeDecSec2string(phypo->otime, tmp_str, DEFAULT_TIME_FORMAT), (time_t) phypo->otime, time_min - MAX_HYPO_ARCHIVE_WINDOW);
            }
            removeHypocenterDescFromHypoList(phypo, &hypo_list, &num_hypocenters);
            // free phypo if not from hyp_assoc_loc array, which is freed later
            for (n = 0; n < MAX_NUM_HYPO; n++) {
                if (phypo == hyp_assoc_loc[n])
                    break;
            }
            if (nhyp >= MAX_NUM_HYPO) // not found
                free(phypo);
            // 20160913 AJL - changed to remove if after analysis window + report_interval, to avoid conflicts if using only_check_for_new_event within next report_interval
        } else if (phypo->hyp_assoc_index < 0 && (time_t) phypo->otime > time_min + report_interval) {
            // event not associated and otime is in analysis window (orphan event?), remove hypocenter
            if (verbose > 0) {
                printf("INFO: Event not associated and otime is in analysis window (orphan event?), remove hypocenter: phypo->hyp_assoc_index %d, %s, phypo->otime %ld, time_mintime_min + report_interval %ld\n",
                        phypo->hyp_assoc_index, timeDecSec2string(phypo->otime, tmp_str, DEFAULT_TIME_FORMAT), (time_t) phypo->otime, time_min + report_interval);
            }
            if (num_hypocenters_associated < MAX_NUM_HYPO) { // do not remove events if num_hypocenters_associated == MAX_NUM_HYPO, may be real event
                if ((time_t) phypo->otime > earliest_time) { // try to avoid removing recent events on startup (BUGGY?)
                    if (verbose > 0) {
                        printf("INFO:                  : (time_t) phypo->otime %ld, earliest_time %ld\n", (time_t) phypo->otime, earliest_time);
                    }
                    removeHypocenterDescFromHypoList(phypo, &hypo_list, &num_hypocenters);
                    // free phypo if not from hyp_assoc_loc array, which is freed later
                    for (n = 0; n < MAX_NUM_HYPO; n++) {
                        if (phypo == hyp_assoc_loc[n])
                            break;
                    }
                    if (n >= MAX_NUM_HYPO) // not found
                        free(phypo);
                }
            }
            //} else if (phypo->hyp_assoc_index >= 0 && (time_t) phypo->otime <= time_min) {
            // event associated and otime is before analysis window, archive hypocenter and remove associated data
            // 20160910 AJL - changed to remove if before analysis window + report_interval, to avoid conflicts if using only_check_for_new_event within next report_interval
        } else if (phypo->hyp_assoc_index >= 0 && (time_t) phypo->otime <= time_min + report_interval) {
            if (verbose > 0) {
                printf("INFO: Event associated and otime is before analysis window: phypo->hyp_assoc_index %d, %s, phypo->otime %ld, time_mintime_min + report_interval %ld\n",
                        phypo->hyp_assoc_index, timeDecSec2string(phypo->otime, tmp_str, DEFAULT_TIME_FORMAT), (time_t) phypo->otime, time_min + report_interval);
            }
            // write hypo data csv string to hypo list persistent archive
            // 20141211 AJL - added
            // hypo list persistent archive stream (TODO: !!! ever growing file!)
            sprintf(outname, "%s/hypolist_persistent.csv", outnameroot_archive);
            if (verbose > 2)
                printf("Opening output file: %s\n", outname);
            FILE * hypoListPersistentStream = fopen_counter(outname, "a");
            if (hypoListPersistentStream == NULL) {
                sprintf(tmp_str, "ERROR: opening output file: %s", outname);
                perror(tmp_str);
                return (-1);
            }
            // write header if file empty
            struct stat stbuf;
            if ((fstat(fileno(hypoListPersistentStream), &stbuf) != 0) || (!S_ISREG(stbuf.st_mode))) {
                sprintf(tmp_str, "ERROR: calling fstat() on output file: %s", outname);
                perror(tmp_str);
            } else {
                if (stbuf.st_size < 1) {
                    printHypoDataHeaderString(hypoDataString);
                    fprintf(hypoListPersistentStream, "%s\n", hypoDataString);
                }
            }
            printHypoDataString(phypo, hypoDataString, 1);
            fprintf(hypoListPersistentStream, "%s\n", hypoDataString);
            fclose_counter(hypoListPersistentStream);
            // remove deData for this hypocenter to avoid phantom locations from later phases
            for (n = num_de_data - 1; n >= 0; n--) {
                TimedomainProcessingData* deData = data_list[n];
                if (deData->is_associated && deData->is_associated == (phypo->hyp_assoc_index + 1)) {
                    removeTimedomainProcessingDataFromDataList(deData, &data_list, &num_de_data);
                    free_TimedomainProcessingData(deData);
                    //data_list[n] = NULL; // 20160802 AJL - memory bug fix?
                }
            }
            // flag event as not actively associated, will insure that event is not persistent or relocated at a later time
            phypo->hyp_assoc_index = -1;
        }
    }
    // write hypocenter archives
    for (nhyp = num_hypocenters - 1; nhyp >= 0; nhyp--) { // reverse time order
        HypocenterDesc* phypo = hypo_list[nhyp];
        // write hypocenter to archive files
        // hypo data csv string
        printHypoDataString(phypo, hypoDataString, 1);
        fprintf(hypoListStream, "%s\n", hypoDataString);
        //printf("DEBUG: =========> hypoDataString: %s\n", hypoDataString);
        // hypo data html string
        // set event background color based on Mwp
        if (!setEventBackgroundColorStringHtml(phypo->mwpLevelStatistics.numLevel, phypo->mwpLevelStatistics.centralValue, "bgcolor=", eventBackgroundColorString,
                phypo->qualityIndicators.quality_code,
                MIN_NUMBER_VALUES_USE, MAG_MIN_HIGHLIGHT_CUTOFF, MAG_MAX_HIGHLIGHT_CUTOFF, MAG_HIGH_CUTOFF, LOC_QUALITY_ACCEPTABLE)) {
            // not enough Mwp readings, try mb
            setEventBackgroundColorStringHtml(phypo->mbLevelStatistics.numLevel, phypo->mbLevelStatistics.centralValue, "bgcolor=", eventBackgroundColorString,
                    phypo->qualityIndicators.quality_code,
                    MIN_NUMBER_VALUES_USE, MAG_MIN_HIGHLIGHT_CUTOFF, MAG_MAX_HIGHLIGHT_CUTOFF, MAG_HIGH_CUTOFF, LOC_QUALITY_ACCEPTABLE);
        }
        printHypoMessageHtmlString(phypo, hypoMessageHtmlString, eventBackgroundColorString, eventBackgroundColorString, phypo->unique_id, phypo->unique_id);
        fprintf(hypoListHtmlStream, "%s\n", hypoMessageHtmlString);
        //printf("DEBUG: =========> hypoMessageHtmlString: %s\n", hypoMessageHtmlString);
    }
    fclose_counter(hypoListStream);
    fprintf(hypoListHtmlStream, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose_counter(hypoListHtmlStream);


    // hypocenter timedomain-processing data ===============================================================

    // 20101115 AJL - Bug fix: moved this loop to here from beginning of function.
    //    Before, data from active hypocenter could be removed since time_min is later next time function is entered; this can cause a second or changed location of same event.
    // remove old timedomain-processing data from list
    for (n = num_de_data - 1; n >= 0; n--) {
        TimedomainProcessingData* deData = data_list[n];
        //if (deData->t_time_t < time_min) {
        if (deData->t_time_t < time_min && !deData->is_associated) { // 20150507 AJL - leave associated data, should be removed when hypocenter archived
            removeTimedomainProcessingDataFromDataList(deData, &data_list, &num_de_data);
            free_TimedomainProcessingData(deData);
            //data_list[n] = NULL; // 20160802 AJL - memory bug fix?
        }
    }



    // clean up ===============================================================

    free(t50ExLevelStatistics);
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
        for (m = 0; m < 3; m++) {
            free(t50ExStatisticsArray[nhyp][m]);
        }
        free(t50ExStatisticsArray[nhyp]);
    }
    free(t50ExStatisticsArray);
    free(numT50ExLevel);
    free(numT50ExLevelMax);
    free(t50ExArray);
    for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++)
        free(t50ExHistogram[nhyp]);
    free(t50ExHistogram);
    fclose_counter(pickStream);
    fclose_counter(t50ExGridDataStream);
    fclose_counter(t50ExStaCodeGridStream);
    //
    if (flag_do_tauc) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++)
            free(taucHistogram[nhyp]);
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            for (m = 0; m < 3; m++) {
                free(taucStatisticsArray[nhyp][m]);
            }
            free(taucStatisticsArray[nhyp]);
        }
        fclose_counter(taucGridDataStrem);
        fclose_counter(taucStaCodeGridStream);
    }
    free(taucLevelStatistics);
    free(taucArray);
    free(taucHistogram);
    free(taucStatisticsArray);
    free(numTaucLevel);
    free(numTaucLevelMax);
    free(tdT50ExLevelStatistics);
    free(numAlarmLevel);
    free(numAlarmLevelMax);
    //
    if (flag_do_mwp) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            free(mwpHistogram[nhyp]);
            for (m = 0; m < 3; m++) {
                free(mwpStatisticsArray[nhyp][m]);
            }
            free(mwpStatisticsArray[nhyp]);
        }
#ifdef USE_MWP_LEVEL_ARRAY
        fclose_counter(mwpGridDataStrem);
        fclose_counter(mwpStaCodeGridStream);
#endif
    }
    free(mwpLevelStatistics);
#ifdef USE_MWP_LEVEL_ARRAY
    free(mwpArray);
#endif
    free(mwpHistogram);
    free(mwpStatisticsArray);
    free(numMwpLevel);
    free(numMwpLevelMax);
    //
    if (flag_do_mb) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            free(mbHistogram[nhyp]);
            for (m = 0; m < 3; m++) {
                free(mbStatisticsArray[nhyp][m]);
            }
            free(mbStatisticsArray[nhyp]);
        }
    }
    free(mbLevelStatistics);
    free(mbHistogram);
    free(mbStatisticsArray);
    free(numMbLevel);
    free(numMbLevelMax);
    //
    if (flag_do_t0) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            free(t0Histogram[nhyp]);
            for (m = 0; m < 3; m++) {
                free(t0StatisticsArray[nhyp][m]);
            }
            free(t0StatisticsArray[nhyp]);
        }
    }
    free(t0LevelStatistics);
    free(t0Histogram);
    free(t0StatisticsArray);
    free(numT0Level);
    free(numT0LevelMax);
    //
    if (flag_do_mwpd) {
        for (nhyp = 0; nhyp < num_hypocenters_associated; nhyp++) {
            free(mwpdHistogram[nhyp]);
            for (m = 0; m < 3; m++) {
                free(mwpdStatisticsArray[nhyp][m]);
            }
            free(mwpdStatisticsArray[nhyp]);
#ifdef USE_MWP_MO_POS_NEG
            for (m = 0; m < 3; m++) {
                free(mwpdMoPosNegStatisticsArray[nhyp][m]);
            }
            free(mwpdMoPosNegStatisticsArray[nhyp]);
#endif
        }
    }
    free(mwpdLevelStatistics);
    free(mwpdHistogram);
    free(mwpdStatisticsArray);
    free(numMwpdLevel);
    free(numMwpdLevelMax);
#ifdef USE_MWP_MO_POS_NEG
    free(mwpdMoPosNegLevelStatistics);
    free(mwpdMoPosNegStatisticsArray);
    free(numMwpdMoPosNegLevel);
    free(numMwpdMoPosNegLevelMax);
#endif

    // info messages
    printf("Info: td_writeTimedomainProcessingReport(): maximum number files opened: %d\n", max_n_open_files);

    return (num_hypocenters_associated);


}
