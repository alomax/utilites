/** *************************************************************************
 * usage_common.c
 *
 * Common routines to perform timedomain-processing processing.
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * created: 2013.02.15
 ***************************************************************************/


#include <stdio.h>
#include <sys/stat.h>

#define EXTERN_MODE
#include "app_lib.h"

static Settings *settings = NULL;

/***************************************************************************
 * init_common_globals():
 * Initialize common program globals.
 ***************************************************************************/
void init_common_globals() {

    verbose = 0;
    details = 0;

    strcpy(agencyId, "AGENCY");

    strcpy(data_input_filter_component_accept, DEFAULT_DATA_INPUT_COMPONENT_ACCEPT);
    data_input_filter_ignore_duplicate_net_sta = DEFAULT_DATA_INPUT_IGNORE_DUPLICATE_NET_STA;
    data_input_filter_ignore_duplicate_sta = DEFAULT_DATA_INPUT_IGNORE_DUPLICATE_STA;

    sliding_window_length = DEFAULT_SLIDING_WINDOW_LENGTH;
    report_interval = DEFAULT_REPORT_INTERVAL;
    report_trigger_pick_window = DEFAULT_REPORT_TRIGGER_PICK_WINDOW;
    report_trigger_min_num_picks = DEFAULT_REPORT_TRIGGER_MIN_NUM_PICKS;
    report_trigger_min_delay_time = DEFAULT_REPORT_TRIGGER_MIN_DELAY_TIME;

    num_sourcename_ignore = 0;

    geogfile = NULL;
    gainfile = NULL;

    num_internet_gain_query = 0;
    num_internet_gain_query_host = 0;
    num_internet_gain_query_type = 0;

    num_internet_station_query = 0;
    num_internet_station_query_host = 0;
    num_internet_station_query_type = 0;

    num_internet_timeseries_query = 0;
    num_internet_timeseries_query_hosturl = 0;
    num_internet_timeseries_query_type = 0;
    num_internet_timeseries_query_sladdr = 0;

    flag_do_grd_vel = 0; // 2014080 AJL - added to allow testing of ground velocity based measures (e.g. fmamp)
    flag_do_grd_vel = 1; // 2014080 AJL - added to allow testing of ground velocity based measures (e.g. fmamp)
    flag_do_mwp = 0;
    flag_do_mwpd = 0;
    flag_do_mb = 0;

    associate_data = 0;
    printIgnoredData = 1;
    alarmNotification = 0;

    no_aref_level_check = 0;


}

/***************************************************************************
 * init_common_properties():
 * Initialize properties from properties file
 ***************************************************************************/
int init_properties(char *propfile) {

    // read properties file
    FILE *fp_prop = fopen(propfile, "r");
    if (fp_prop == NULL) {
        fprintf(stdout, "Info: Cannot open application properties file: %s\n", propfile);
        settings = NULL;
        return (0);
    }
    settings = settings_open(fp_prop);
    fclose(fp_prop);
    if (settings == NULL) {
        fprintf(stdout, "ERROR: Reading application properties file: %s\n", propfile);
        return (-1);
    }

    return (0);

}

/***************************************************************************
 * init_common_processing():
 * Initialize common processing items.
 ***************************************************************************/
int init_common_processing(char *outpath) {

    if (td_process_timedomain_processing_init(settings, outpath, geogfile, gainfile, verbose) < 0)
        return -1;

    app_prop_settings = settings;

    int int_param;
    //double double_param;
    char out_buf[SETTINGS_MAX_STR_LEN];

    // data input
    if (settings_get(app_prop_settings, "DataInput", "filter.component.accept", out_buf, SETTINGS_MAX_STR_LEN)) {
        strcpy(data_input_filter_component_accept, out_buf);
    }
    printf("Info: property set: [DataInput] filter.component.accept %s\n", data_input_filter_component_accept);
    //
    if ((int_param = settings_get_int(app_prop_settings, "DataInput", "filter.ignore_duplicate_net_sta")) != INT_INVALID)
        data_input_filter_ignore_duplicate_net_sta = int_param;
    printf("Info: property set: [DataInput] filter.ignore_duplicate_net_sta %d\n", data_input_filter_ignore_duplicate_net_sta);
    //
    if ((int_param = settings_get_int(app_prop_settings, "DataInput", "filter.ignore_duplicate_sta")) != INT_INVALID)
        data_input_filter_ignore_duplicate_sta = int_param;
    printf("Info: property set: [DataInput] filter.ignore_duplicate_sta %d\n", data_input_filter_ignore_duplicate_sta);

    // report
    if ((int_param = settings_get_int(app_prop_settings, "Report", "report.sliding_window.length")) != INT_INVALID)
        sliding_window_length = int_param;
    printf("Info: property set: [Report] report.sliding_window.length %d\n", sliding_window_length);
    //
    if ((int_param = settings_get_int(app_prop_settings, "Report", "report.interval")) != INT_INVALID)
        report_interval = int_param;
    printf("Info: property set: [Report] report.interval %d\n", report_interval);
    //
    if ((int_param = settings_get_int(app_prop_settings, "Report", "report.trigger.pick_window")) != INT_INVALID)
        report_trigger_pick_window = int_param;
    printf("Info: property set: [Report] report.trigger.pick_window %d\n", report_trigger_pick_window);
    //
    if ((int_param = settings_get_int(app_prop_settings, "Report", "report.trigger.min_num_picks")) != INT_INVALID)
        report_trigger_min_num_picks = int_param;
    printf("Info: property set: [Report] report.trigger.min_num_picks %d\n", report_trigger_min_num_picks);
    //
    if ((int_param = settings_get_int(app_prop_settings, "Report", "report.trigger.min_delay_time")) != INT_INVALID)
        report_trigger_min_delay_time = int_param;
    printf("Info: property set: [Report] report.trigger.min_delay_time %d\n", report_trigger_min_delay_time);

    return (0);

}

/***************************************************************************
 * usage_common():
 * Print the usage message for parameters common to all applications.
 ***************************************************************************/
void usage_common() {

    fprintf(stdout,
            " ## General program options ##\n"
            "\n"
            " -V             Report program version\n"
            " -h             Show this usage message\n"
            " -v             Be more verbose, multiple flags can be used\n"
            "\n"
            " -p             Print details of mseed header, multiple flags can be used\n"
            " -agency-id     public id for originatin agency in XML messages (e.g. net.alomax)\n"

            " -i template    Specify source names to ignore with comma separated list of NET_STA format prefixes 'NET[_STA[_LOC[_CHAN]]]' (e.g. TA_,AU_W)\n"
            "\n"
            " -report-sliding-window-length  report sliding window length (3600 sec)\n"
            " -report-interval  report interval (60 sec)\n"
            "\n"
            " -g geogfile    file with geographic coordinates of stations (Network Station Latitude Longitude Elevation(m))\n"
            " -sta-query-type type of query station web service (e.g. FDSN_WS_STATION, IRIS_WS_STATION)\n"
            " -sta-query-host host name for station web service (e.g. www.iris.edu)\n"
            " -sta-query      query root for station web service (e.g. ws/station/query)\n"
            "                 -sta-query-type, -sta-query-host and -sta-query can be repeated in corresponding triplets\n"
            "\n"
            " -pz gainfile   file with gain info for channels (Network Station Location Channel Year, DayOfYear, Gain, Freq, Type)\n"
            " -pz-query-type type of pole-zero response (SEED-RESP format) web service (e.g. FDSN_WS_STATION, IRIS_WS_RESP)\n"
            " -pz-query-host host name for pole-zero response (SEED-RESP format) web service (e.g. www.iris.edu)\n"
            " -pz-query      query root for pole-zero response (SEED-RESP format) web service (e.g. ws/resp/query)\n"
            "                 -pz-query-type, -pz-query-host and -pz-query can be repeated in corresponding triplets\n"
            "\n"
            " -timeseries-query-type type of timeseries web service (e.g. IRIS_WS_TIMESERIES)\n"
            " -timeseries-query-host host name for timeseries web service (e.g. www.iris.edu)\n"
            " -timeseries-query      query root for timeseries web service (e.g. ws/timeseries/query)\n"
            " -timeseries-query-sladdr      host:port to associate data channel with timeseries web service (e.g. rtserve.iris.washington.edu:18000)\n"
            "                 -timeseries-query-type, -timeseries-query-host, -timeseries-query and -timeseries-query-sladdr, can be repeated in corresponding quadruples\n"
            "\n"
            " -mwp           calculate Mwp magnitude for located events (requires -a and -pz, and optionally -pz-query-host and -pz-query)\n"
            " -mwpd          calculate Mwpd magnitude for located events (requires -a and -pz, and optionally -pz-query-host and -pz-query)\n"
            " -mb            calculate mB magnitude for located events (requires -a and -pz, and optionally -pz-query-host and -pz-query)\n"
            " -a             associate picks and locate events based on P travel-time; use associated event data for event parameters statistics\n"
            " -si            suppress output of ignored pick data\n"
            " -alarm         sound audible alarm on new event location\n"
            "\n"
            " -no_aref_level_check            no check of aref level (use to support microseismicity monitoring (e.g. ignore quality of tsunami discriminants)\n"
            "\n");

} /* End of usage_common() */


static char sourcename_ignore_tokens[MAX_NUM_SOURCENAME_IGNORE * 64];

/***************************************************************************
 * parameter_proc():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
int parameter_proc_common(int *poptind, char **argvec, char *PACKAGE, char *VERSION, char *VERSION_DATE, void (*usage)()) {

    // try to process command line argument

    if (strcmp(argvec[*poptind], "-V") == 0) {
        fprintf(stderr, "%s version: %s (%s)\n", PACKAGE, VERSION, VERSION_DATE);
        exit(0);
    } else if (strcmp(argvec[*poptind], "-h") == 0) {
        (*usage)();
        exit(0);
    } else if (strcmp(argvec[*poptind], "-agency-id") == 0) {
        strcpy(agencyId, argvec[++(*poptind)]);
    } else if (strcmp(argvec[*poptind], "-o") == 0) {
        outpath_param = argvec[++(*poptind)];
    } else if (strcmp(argvec[*poptind], "-report-sliding-window-length") == 0) {
        sliding_window_length = atoi(argvec[++(*poptind)]);
    } else if (strcmp(argvec[*poptind], "-report-interval") == 0) {
        report_interval = atoi(argvec[++(*poptind)]);
    } else if (strcmp(argvec[*poptind], "-i") == 0) {
        strcpy(sourcename_ignore_tokens, argvec[++(*poptind)]);
        // break sourcename string into tokens
        num_sourcename_ignore = 0;
        char *str_pos = strtok(sourcename_ignore_tokens, ",");
        if (str_pos == NULL) {
            strcpy(sourcename_ignore[num_sourcename_ignore], sourcename_ignore_tokens);
            //printf("DEBUG: sourcename_ignore[%d]: %s\n", num_sourcename_ignore, sourcename_ignore[num_sourcename_ignore]);
            num_sourcename_ignore++;
        } else {
            while (1) {
                strcpy(sourcename_ignore[num_sourcename_ignore], str_pos);
                //printf("DEBUG: sourcename_ignore[%d]: %s\n", num_sourcename_ignore, sourcename_ignore[num_sourcename_ignore]);
                num_sourcename_ignore++;
                if (num_sourcename_ignore >= MAX_NUM_SOURCENAME_IGNORE) break;
                str_pos = strtok(NULL, ",");
                if (str_pos == NULL) break;
            }
        }
    } else if (strcmp(argvec[*poptind], "-g") == 0) {
        geogfile = argvec[++(*poptind)];
    } else if (strcmp(argvec[*poptind], "-sta-query") == 0) {
        if (num_internet_station_query < MAX_NUM_INTERNET_QUERY) {
            internet_station_query[num_internet_station_query] = argvec[++(*poptind)];
            num_internet_station_query++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-sta-query-host") == 0) {
        if (num_internet_station_query_host < MAX_NUM_INTERNET_QUERY) {
            internet_station_query_host[num_internet_station_query_host] = argvec[++(*poptind)];
            num_internet_station_query_host++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-sta-query-type") == 0) {
        if (num_internet_station_query_type < MAX_NUM_INTERNET_QUERY) {
            internet_station_query_type[num_internet_station_query_type] = argvec[++(*poptind)];
            num_internet_station_query_type++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-pz") == 0) {
        gainfile = argvec[++(*poptind)];
    } else if (strcmp(argvec[*poptind], "-pz-query") == 0) {
        if (num_internet_gain_query < MAX_NUM_INTERNET_QUERY) {
            internet_gain_query[num_internet_gain_query] = argvec[++(*poptind)];
            num_internet_gain_query++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-pz-query-host") == 0) {
        if (num_internet_gain_query_host < MAX_NUM_INTERNET_QUERY) {
            internet_gain_query_host[num_internet_gain_query_host] = argvec[++(*poptind)];
            num_internet_gain_query_host++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-pz-query-type") == 0) {
        if (num_internet_gain_query_type < MAX_NUM_INTERNET_QUERY) {
            internet_gain_query_type[num_internet_gain_query_type] = argvec[++(*poptind)];
            num_internet_gain_query_type++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-timeseries-query") == 0) {
        if (num_internet_timeseries_query < MAX_NUM_INTERNET_QUERY) {
            internet_timeseries_query[num_internet_timeseries_query] = argvec[++*poptind];
            num_internet_timeseries_query++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-timeseries-query-hosturl") == 0) {
        if (num_internet_timeseries_query_hosturl < MAX_NUM_INTERNET_QUERY) {
            internet_timeseries_query_hosturl[num_internet_timeseries_query_hosturl] = argvec[++*poptind];
            num_internet_timeseries_query_hosturl++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-timeseries-query-type") == 0) {
        if (num_internet_timeseries_query_type < MAX_NUM_INTERNET_QUERY) {
            internet_timeseries_query_type[num_internet_timeseries_query_type] = argvec[++*poptind];
            num_internet_timeseries_query_type++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-timeseries-query-sladdr") == 0) {
        if (num_internet_timeseries_query_sladdr < MAX_NUM_INTERNET_QUERY) {
            internet_timeseries_query_sladdr[num_internet_timeseries_query_sladdr] = argvec[++*poptind];
            num_internet_timeseries_query_sladdr++;
        } else {
            fprintf(stderr, "Too many occurrences of option: %s\n", argvec[*poptind]);
            exit(1);
        }
    } else if (strcmp(argvec[*poptind], "-mwp") == 0) {
        flag_do_mwp = 1;
    } else if (strcmp(argvec[*poptind], "-mwpd") == 0) {
        flag_do_mwpd = 1;
        flag_do_t0 = 1; // T0 required for Mwpd
        flag_do_tauc = 1; // also do TauC
    } else if (strcmp(argvec[*poptind], "-mb") == 0) {
        flag_do_mb = 1;
    } else if (strcmp(argvec[*poptind], "-a") == 0) {
        associate_data = 1;
    } else if (strcmp(argvec[*poptind], "-no_aref_level_check") == 0) {
        no_aref_level_check = 1;
    } else if (strcmp(argvec[*poptind], "-si") == 0) {
        printIgnoredData = 0;
    } else if (strcmp(argvec[*poptind], "-alarm") == 0) {
        alarmNotification = 1;
    } else if (strncmp(argvec[*poptind], "-v", 2) == 0) { // must be last argument checked to avoid conflict with other arguments starting with same letter
        verbose += strspn(&argvec[*poptind][1], "v");
    } else if (strncmp(argvec[*poptind], "-p", 2) == 0) { // must be last argument checked to avoid conflict with other arguments starting with same letter
        details += strspn(&argvec[*poptind][1], "p");
    } else {
        return (0);
    }

    return (1);

} /* End of parameter_proc() */

/** *************************************************************************
 * timet2timestr()
 *
 * Returns a time string representation of a struct tm
 ** *********************************************************************** **/

char* tm2timestring(struct tm* timep, double secfrac, char* timestring) {

    snprintf(timestring, 20, "%04d.%02d.%02d-%02d:%02d:%02d.%01.0f", timep->tm_year + 1900, timep->tm_mon + 1, timep->tm_mday, timep->tm_hour, timep->tm_min, timep->tm_sec, secfrac);

    return (timestring);

}

/** *************************************************************************
 * timet2path()
 *
 * Makes a directory path based on a path prefix and the year/month/day/time of a struct tm and
 * returns a string representation of the path
 ** *********************************************************************** **/

char* tm2path(char* prefix, struct tm* timep, char* pathstring) {

    char dirname[STANDARD_STRLEN];

    sprintf(dirname, "%s/%04d", prefix, timep->tm_year + 1900);
    mkdir(dirname, 0755);
    sprintf(dirname, "%s/%04d/%02d", prefix, timep->tm_year + 1900, timep->tm_mon + 1);
    mkdir(dirname, 0755);
    sprintf(dirname, "%s/%04d/%02d/%02d", prefix, timep->tm_year + 1900, timep->tm_mon + 1, timep->tm_mday);
    mkdir(dirname, 0755);

    sprintf(pathstring, "%s/%04d/%02d/%02d/%04d.%02d.%02d.%02d.%02d.%02d", prefix, timep->tm_year + 1900, timep->tm_mon + 1, timep->tm_mday, timep->tm_year + 1900, timep->tm_mon + 1, timep->tm_mday, timep->tm_hour, timep->tm_min, timep->tm_sec);
    mkdir(pathstring, 0755);

    return (pathstring);

}

/***************************************************************************
 * td_process_timedomain_processing_cleanup:
 *
 * Do necessary cleanup.
 ***************************************************************************/
void app_lib_cleanup() {

    settings_delete(settings);

}
