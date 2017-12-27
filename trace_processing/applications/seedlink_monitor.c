/** *************************************************************************
 * seedlink_monitor.c
 *
 * Generic routines to perform timedomain-processing processing on Mini-SEED records from a SeedLink server.
 *
 * Opens a user specified seedlink connection, parses the returned Mini-SEED records and processes
 * each record.
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * created: 2009.03.03
 ************************************************************************* **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef WIN32
#include <signal.h>
static void term_handler(int sig);
#endif

#include <libslink.h>
#include <libmseed.h>

#include "../ran1/ran1.h"
#include "../net/net_lib.h"
#include "../response/response_lib.h"
#include "app_lib.h"
#include "seedlink_monitor.h"

#define VERSION EARLY_EST_MONITOR_VERSION
#define PACKAGE "seedlink_monitor"
#define VERSION_DATE EARLY_EST_MONITOR_VERSION_DATE

//static int reclen = -1;
static char outpath[STANDARD_STRLEN];
static char *command = NULL;


static char begintime[STANDARD_STRLEN];
static int ibackfill = 0;

static int sendMail = 0;
static char *sendMailParams = NULL;

static int num_new_loc_picks = 0;
static double count_new_use_loc_picks_cutoff_time = -FLT_MAX;



//static int write_raw = 0;

//static char sourcename[50];
//static char sourcename_ignore[50];

static int parameter_proc(int argcount, char **argvec);
static void usage(void);
int seedlink_monitor_cleanup();

static char sourcename_list[MAX_NUM_SOURCES][64];
static double source_samplerate_list[MAX_NUM_SOURCES];

static char srcnames_duplicate[MAX_NUM_SOURCES * 64];

static int source_packet_error_reported[MAX_NUM_SOURCES];

//static FILE* outstream_RAW[MAX_NUM_SOURCES];
//static FILE* outstream_HF[MAX_NUM_SOURCES];
/*
static FILE** outstream;
static void record_handler (char *record, int write_reclen, void *srcname) {
// set source id
        int n = 0;
        for (n = 0; n < num_sources; n++) {
                if (strcmp(sourcename_list[n], srcname) == 0) {
                        break;
                }
        }
        if ( n >= MAX_NUM_SOURCES || fwrite(record, write_reclen, 1, outstream[n]) != 1 )
        {
                sl_log(2, 0, "Error writing %s to output file.\n", (char *) srcname);
        }
}
 */

//static TimedomainProcessingData** local_data_list = NULL;
//static int local_num_de_data = 0;


#ifndef WIN32
#include <signal.h>

static void term_handler(int sig);
#else
#define strcasecmp _stricmp
#endif

// 20200507 AJL - Added multiple slconn's to support connections to multiple SeedLink servers.
#define MAX_NUM_SLCON 99
static int num_slconn = 0;
static SLCD* slconn[MAX_NUM_SLCON]; // connection parameters
static int slconn_terminated[MAX_NUM_SLCON]; // flag of terminated slconn's
static int requested_sl_terminate = 0;
static char statefile[MAX_NUM_SLCON][STANDARD_STRLEN]; // state file for saving/restoring state
static char streamfile[MAX_NUM_SLCON][STANDARD_STRLEN]; // stream list file for configuring streams
static char multiselect[MAX_NUM_SLCON][10 * STANDARD_STRLEN];
static char selectors[MAX_NUM_SLCON][STANDARD_STRLEN];
static char location_list[MAX_NUM_SLCON][STANDARD_STRLEN];

static int check_all_sl_servers = 1;
static time_t time_checked_all = -1;

static void packet_handler(char *msrecord, int packet_type, int seqnum, int packet_size, SLCD* slconn_packet, char* location_list_slconn);
static int parameter_proc(int argcount, char **argvec);
static void usage(void);

int main(int argc, char **argv) {

    int nslconn;

    SLpacket* slpack;
    int seqnum;
    int ptype;

#ifndef WIN32
    // Signal handling, use POSIX calls with standardized semantics
    struct sigaction sa;

    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;

    sa.sa_handler = term_handler;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGQUIT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    sa.sa_handler = SIG_IGN;
    sigaction(SIGHUP, &sa, NULL);
    sigaction(SIGPIPE, &sa, NULL);
#endif


    // initialize common globals
    init_common_globals();

    // read properties file
    char propfile[1024];
    sprintf(propfile, "%s.prop", PACKAGE);
    if (init_properties(propfile) < 0) {
        return -1;
    }

    // Process given parameters (command line and parameter file)
    if (parameter_proc(argc, argv) < 0) {
        fprintf(stderr, "ERROR: Parameter processing failed\n\n");
        fprintf(stderr, "Try '-h' for detailed help\n");
        return -1;
    }

    int RandomNumSeed = 9589;
    SRAND_FUNC(RandomNumSeed);
    if (0)
        test_rand_int();

    //printf("outpath <%s>\n", outpath);
    mkdir(outpath, 0755);

    // initialize time domain processing
    if (init_common_processing(outpath) < 0) {
        fprintf(stderr, "ERROR: Initialization of processing items (init_common_processing()) failed\n\n");
        fprintf(stderr, "Try '-h' for detailed help\n");
        return -1;
    }

    // initialize local statics
    strcpy(srcnames_duplicate, "\0");
    // reset last_data_before_report_error list so warning message above will be produced at least once in a report interval
    for (int n = 0; n < MAX_NUM_SOURCES; n++) {
        source_packet_error_reported[n] = 0;
    }

#if 0
    // TEST - force Segmentation fault, used to test external program restart functionality
    double *pxdf = (double *) 999999999999;
    double xcv = *pxdf / 0.0;
    printf("double xcv = %g\n", xcv);
#endif

    // Loop with the connection manager
    int n_slconn_terminated = 0;
    while (n_slconn_terminated < num_slconn) {
        // update ellapsed time and check all servers flag
        time_t time_now = (time_t) slp_dtime_curr();
        if (time_now > time_checked_all + report_interval) {
            time_checked_all = time_now;
            check_all_sl_servers = 1;
        }
        n_slconn_terminated = 0;
        int packet_returned = 0;
        int return_val;
        int iloop_count = 0;
        // check for data packets from each slconn
        for (nslconn = 0; nslconn < num_slconn; nslconn++) {
            // skip slconn's that are down if not first pass
            // 20110629 AJL - added this check to help prevent cumulating large delays (sleep of 0.5s) in sl_collect_nb if a SL conn is down.
            if (!check_all_sl_servers && slconn[nslconn]->stat->sl_state == SL_DOWN && !requested_sl_terminate)
                continue;
            packet_returned = 1;
            iloop_count = 0;
            // 20110324 AJL - added this loop so many packets can be received from each SL conn before checking next SL conn,
            //    this helps to prevent cumulating large delays (sleep of 0.5s) in sl_collect_nb if a SL conn is down.
            while (iloop_count++ < 100 && packet_returned) {
                packet_returned = 0;
                if ((return_val = sl_collect_nb(slconn[nslconn], &slpack) == SLPACKET)) {
                    packet_returned = 1;
                    ptype = sl_packettype(slpack);
                    seqnum = sl_sequence(slpack);
                    packet_handler((char *) &(slpack->msrecord), ptype, seqnum, SLRECSIZE, slconn[nslconn], location_list[nslconn]);
                    // It would be possible to send an in-line INFO request here with sl_request_info().
                } else if (requested_sl_terminate && return_val == SLTERMINATE) {
                    slconn_terminated[nslconn] = 1;
                }
            }
            n_slconn_terminated += slconn_terminated[nslconn];
        }
        check_all_sl_servers = 0;
        if (!packet_returned) { // if no data packets returned, then sleep
            // 20110324 AJL usleep(1000);
            // 20111107 AJL usleep(10);
            if (usleep(10000) < 0)
                fprintf(stderr, "ERROR: seedlink_monitor.main(): usleep(): %s\n", strerror(errno));
        }
    }

    seedlink_monitor_cleanup();

    return 0;

} // End of main()

/***************************************************************************
 * seedlink_monitor_cleanup:
 *
 * Do necessary cleanup.
 ***************************************************************************/
int seedlink_monitor_cleanup() {

    if (verbose) {
        printf("Info: Performing cleanup...\n");
    }

    // Make sure everything is shut down and save the state file
    int nslconn;
    for (nslconn = 0; nslconn < num_slconn; nslconn++) {
        if (slconn[nslconn]->link != -1)
            sl_disconnect(slconn[nslconn]);
        if (statefile[nslconn][0] != '\0')
            sl_savestate(slconn[nslconn], statefile[nslconn]);
        //sl_freeslcd(slconn);  // gives Invalid free() error in valgrind
    }

    // make sure time domain processing is cleaned up
    /*
    int n;
    for (n = 0; n < num_sources; n++) {
            if (write_raw_brb)
                    fclose(outstream_BRB[n]);
            else
                    fclose(outstream_HF[n]);
    }
     */
    td_process_timedomain_processing_cleanup();
    /*
    // remove local timedomain-processing data from list
    while (local_num_de_data > 0) {
            TimedomainProcessingData* deData = local_data_list[local_num_de_data - 1];
            //printf("DEBUG: FREE: deData %ld\n", deData);
            removeTimedomainProcessingDataFromDataList(deData, &local_data_list, &local_num_de_data);
            //printf("DEBUG: FREE: deData %ld  local_num_de_data %d\n", deData, local_num_de_data);
            //printf("DEBUG: FREE: deData->pickData %ld\n", deData->pickData);
            free_TimedomainProcessingData(deData);
    }
     */
    app_lib_cleanup();

    return 0;

}

/** *************************************************************************
 * packet_handler():
 * Process a received packet based on packet type.
 ************************************************************************* **/

static void packet_handler(char *msrecord, int packet_type, int seqnum, int packet_size, SLCD* slconn_packet, char* location_list_slconn) {

    static SLMSrecord* sl_msr = NULL;
    static char outpath_report[STANDARD_STRLEN];
    static char command_report[4 * STANDARD_STRLEN];

    // data time range
    static time_t time_min = -1;
    static time_t time_max = -1;
    static time_t last_report_time = -1;
    static time_t next_report_time = -1;
    static time_t last_report_command_time = -1;
    static time_t last_time_check_for_new_event = INT_MAX / 2;

    static int totalrecs = 0;
    static int totalsamps = 0;

    static char color[64];

    static char srcname[64];
    static char msr_net[4], msr_sta[7];
    static char msr_loc[4], msr_chan[5];
    static char net_sta[50];
    static char timestr[64];

    static char srcname_copy[64];
    static char duplicate_name_copy[64];

    int source_id = 0;


    //  The following is dependent on the packet type values in libslink.h
    char *type[] = {"Data", "Detection", "Calibration", "Timing", "Message", "General", "Request", "Info", "Info (terminated)", "KeepAlive"};


    // system time
    double dtime; // Epoch time
    double secfrac; // Fractional part of epoch time
    time_t isystem_time_t; // Integer part of epoch time
    static char timestamp[64];
    struct tm *timep;
    // Build a current local time string
    dtime = slp_dtime_curr();
    secfrac = (double) ((double) dtime - (int) dtime);
    isystem_time_t = (time_t) dtime;
    timep = localtime(&isystem_time_t);
    tm2timestring(timep, secfrac, timestamp);

    // Process waveform data
    if (verbose > 2)
        printf("\nDEBUG: type %s  packet_type %d  SLDATA %d\n", type[packet_type], packet_type, SLDATA);
    if (packet_type == SLDATA) {
        sl_log(0, 2, "%s, seq %d, Received %s blockette:\n", timestamp, seqnum, type[packet_type]);

        // parse without unpacking since we may not use this data
        if (sl_msr_parse(slconn_packet->log, msrecord, &sl_msr, 1, 0) == NULL) {
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }

        if (verbose > 2 || details)
            sl_msr_print(slconn_packet->log, sl_msr, details);



        // Check source name string
        // Generate clean identifier strings */
        sl_strncpclean(msr_net, sl_msr->fsdh.network, 2);
        sl_strncpclean(msr_sta, sl_msr->fsdh.station, 5);
        sl_strncpclean(msr_loc, sl_msr->fsdh.location, 2);
        sl_strncpclean(msr_chan, sl_msr->fsdh.channel, 3);

        // clean up location code
        char loc_str[3] = "\0\0\0";
        strcpy(loc_str, msr_loc);
        if (isspace(loc_str[0]))
            loc_str[0] = '-';
        if (isspace(loc_str[1]))
            loc_str[1] = '-';
        if (strlen(loc_str) < 2)
            strcpy(loc_str, "--");
        strcpy(msr_loc, loc_str);
        /*
        if (msr_net[0] != '\0')
            strcat(msr_net, "_");
        if (msr_sta[0] != '\0')
            strcat(msr_sta, "_");
        if (prtloc[0] == '\0')
            strcpy(prtloc, "--");
         * */
        // Build the source name string
        sprintf(srcname, "%.2s_%.5s_%.2s_%.3s", msr_net, msr_sta, msr_loc, msr_chan);

        // skip this record if no data samples
        if (verbose > 2)
            printf("DEBUG: sl_msr->fsdh.num_samples %d  sl_msr->numsamples %d  sl_msr->unpackerr %d\n",
                sl_msr->fsdh.num_samples, sl_msr->numsamples, sl_msr->unpackerr);
        if (sl_msr->fsdh.num_samples < 1) {
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }

        /* not needed since selectors can be used
        // skip this record if source name does not begin with requested source name
                        if (sourcename[0] != '\0' && strstr(srcname, sourcename) != srcname)
                        continue;
        // skip this record if source name begins with sourcename_ignore
                        if (sourcename_ignore[0] != '\0' && strstr(srcname, sourcename_ignore) == srcname)
                        continue;
         */
        // skip this record if source name begins with sourcename_ignore
        // 20130128 AJL - added
        if (num_sourcename_ignore > 0) {
            int n_ignore;
            for (n_ignore = 0; n_ignore < num_sourcename_ignore; n_ignore++) {
                //printf("DEBUG: strstr [%d]: %s | %s\n", n_ignore, srcname, sourcename_ignore[n_ignore]);
                if (strstr(srcname, sourcename_ignore[n_ignore]) == srcname) {
                    //printf("DEBUG: IGNORED!\n");
                    sl_msr_free(&sl_msr);
                    sl_msr = NULL;
                    return;
                }
            }
        }

        // skip if not requested component orientation // 20160802 AJL - modified
        // 20160802 AJL - changed to comma separated list of channels to accept
        if (strrchr(data_input_filter_component_accept, msr_chan[2]) == NULL) {
            if (verbose > 1)
                sl_log(0, 1, "Info: Channel orientation in %s (%s) not in %s - Ignored!\n", srcname, msr_chan, data_input_filter_component_accept);
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }


        // check for sourcename in list
        int isource;
        if (verbose > 2)
            printf("DEBUG: num_sources %d\n", num_sources_total);
        int duplicate = 0;
        if (strstr(srcnames_duplicate, srcname) != NULL) { // already flagged as duplicate
            duplicate = 3;
        } else {
            for (isource = 0; isource < num_sources_total; isource++) {
                // if perfect match, OK
                if (strcmp(sourcename_list[isource], srcname) == 0) {
                    // confirm agreement of sample rates
                    // 20160608 AJL  if (source_samplerate_list[isource] == sl_msr_dnomsamprate(sl_msr)) {
                    if (fabs(source_samplerate_list[isource] - sl_msr_dnomsamprate(sl_msr)) / source_samplerate_list[isource] < 0.01) { // 20160608 AJL - changed to 1% accuracy
                        break;
                    } else {
                        // different sample rate !!! NOTE: this should not happen but does,
                        //    try: java -Xmx768m net.alomax.seisgram2k.SeisGram2K -seedlink "rtserve.iris.washington.edu:18000#AT_SVW2:BHZ#7200"
                        sl_log(0, 1, "Warning: conflicting sample rates: %f in ignored %s vs. %f in active %s: this should not happen!\n",
                                sl_msr_dnomsamprate(sl_msr), srcname, source_samplerate_list[isource], sourcename_list[isource]);
                        //duplicate = 4;
                        channelParameters[isource].error = channelParameters[isource].error | ERROR_DIFFERENT_SAMPLE_RATES;
                        channelParameters[isource].count_conflicting_dt++;
                        sl_msr_free(&sl_msr);
                        sl_msr = NULL;
                        return;
                        //break;
                    }
                }
                // 20160802 AJL - only apply net, sta check if same orientation (last char of srcname)
                // check on orientation is done by filter.component.accept properties file parameter, see data_input_filter_component_accept above
                if (sourcename_list[isource][strlen(sourcename_list[isource]) - 1] == srcname[strlen(srcname) - 1]) {
                    if (data_input_filter_ignore_duplicate_net_sta) { // 20130529 AJL - can disable these check with properties file ([DataInput] filter.ignore_duplicate_net_sta)
                        // check for match of network, station only
                        sprintf(net_sta, "%.2s_%.5s_", msr_net, msr_sta);
                        //strcat(strcat(strcpy(net_sta, msr_net), "_"), msr_sta);
                        if (verbose > 2)
                            printf("DEBUG: srcname %s,  net_sta %s,  n %d, sourcename_list[isource] %s,\n", srcname, net_sta, isource, sourcename_list[isource]);
                        if (strstr(sourcename_list[isource], net_sta) != NULL) {
                            duplicate = 1;
                            break;
                        }
                    }
                    if (data_input_filter_ignore_duplicate_sta) { // 20130529 AJL - can disable these check with properties file ([DataInput] filter.ignore_duplicate_sta)
                        // check for possible match of station only
                        sprintf(net_sta, "_%.5s_", msr_sta);
                        if (strstr(sourcename_list[isource], net_sta) != NULL) {
                            duplicate = 2;
                            break;
                        }
                    }
                }
            }
            source_id = isource;
        }
        // act on duplicate net_sta
        if (duplicate) {
            int replace = 0;
            if (duplicate == 3) { // already flagged as duplicate
                ;
            } else if (duplicate == 4) { // different sample rate !!!
                ;
            } else {
                if (duplicate == 1) { // not yet encountered
                    sl_log(0, 1, "Duplicate net_sta \"%s\" (in %s, active is %s)", net_sta, srcname, sourcename_list[source_id]);
                    // check if duplicate has preferred location (set by -locs program parameter)
                    if (location_list_slconn[0] != '\0') {
                        strcpy(duplicate_name_copy, srcname);
                        strtok(duplicate_name_copy, "_");
                        strtok(NULL, "_");
                        char* duplicate_loc = strtok(NULL, "_");
                        strcpy(srcname_copy, sourcename_list[source_id]);
                        strtok(srcname_copy, "_");
                        strtok(NULL, "_");
                        char* srcname_loc = strtok(NULL, "_");
                        char* idup = strstr(location_list_slconn, duplicate_loc);
                        char* icurr = strstr(location_list_slconn, srcname_loc);
                        if (
                                (idup != NULL && icurr == NULL) // duplicate location is in preferred list, active location is not
                                || (idup != NULL && icurr != NULL && icurr - idup > 0) // duplicate location is preferred over active location
                                ) {
                            sl_log(0, 1, " -> active loc \"%s\" replaced by preferred loc \"%s\"\n", srcname_loc, duplicate_loc);
                            // replace active with duplicate, which will be new source
                            replace = 1;
                            // put current source name in duplicates list, later occurrences will be skipped
                            strcat(srcnames_duplicate, "$");
                            strcat(srcnames_duplicate, sourcename_list[source_id]);
                            // make active source name in sourcename_list invalid
                            strcpy(sourcename_list[source_id], "XXX");
                            // flag source as inactive_duplicate
                            channelParameters[source_id].inactive_duplicate = 1;
                            // ensures that new source is created, old active memory will remain until clean-up
                            source_id = num_sources_total;
                        } else {
                            sl_log(0, 1, " -> Ignored.\n");
                        }
                    }
                } else if (duplicate == 2) { // not yet encountered
                    sl_log(0, 1, "Possible duplicate sta \"%s\" (in %s, active is %s) -> Ignored.\n", net_sta, srcname, sourcename_list[source_id]);
                }
                if (!replace) {
                    // put source name in duplicates
                    strcat(srcnames_duplicate, "$");
                    strcat(srcnames_duplicate, srcname);
                }
            }
            if (!replace) {
                sl_msr_free(&sl_msr);
                sl_msr = NULL;
                return;
            }
        }
        // not found, new source name
        if (source_id == num_sources_total) {
            if (num_sources_total < MAX_NUM_SOURCES) { // sourcename not in list
                strcpy(sourcename_list[source_id], srcname);
                source_samplerate_list[source_id] = sl_msr_dnomsamprate(sl_msr);
                num_sources_total++;
                sl_log(0, 1, "New source added: %s: n=%d dt=%f\n", srcname, source_id, 1.0 / sl_msr_dnomsamprate(sl_msr));
                // open output file for this source
                /*
                    if (write_raw) {
                        sprintf(outname, "%s/%s.RAW.mseed", outpath, srcname);
                        outstream_RAW[source_id] = fopen(outname, "w");
                    } else {
                        sprintf(outname, "%s/%s.HF.mseed", outpath, srcname);
                        outstream_HF[source_id] = fopen(outname, "w");
                    }
                 */
            } else { // too many sources
                sl_log(2, 0, "Too many sources (%d max), skipping %s\n", MAX_NUM_SOURCES, srcname);
                sl_msr_free(&sl_msr);
                sl_msr = NULL;
                return;
            }
        }

        // generate time string
        double dstart_time = sl_msr_depochstime(sl_msr);
        double dsec = (double) ((double) dstart_time - (int) dstart_time);
        time_t start_time_t = (time_t) dstart_time;
        time_t end_time_t = start_time_t + (int) ((double) sl_msr->fsdh.num_samples / sl_msr_dnomsamprate(sl_msr));
        //if (strstr(net_sta, "SANT") != NULL)
        //    printf("DEBUG: sl_msr->fsdh.num_samples %d  sl_msr_dnomsamprate(sl_msr) %f  start_time_t %ld  end_time_t %ld\n", sl_msr->fsdh.num_samples, sl_msr_dnomsamprate(sl_msr), start_time_t, end_time_t);
        tm2timestring(gmtime(&start_time_t), dsec, timestr);

        // check time range
        if (time_min < 0 || start_time_t < time_min) {
            time_min = start_time_t;
            if (verbose > 2) {
                printf("DEBUG: time_min found %s\n", timestr);
                sl_msr_print(slconn_packet->log, sl_msr, details);
            }
        }

        //double data_latency = difftime(time_max, end_time_t);
        double data_latency = difftime(isystem_time_t, end_time_t); // 20110105 AJL
        //if (strstr(net_sta, "SANT") != NULL)
        //    printf("DEBUG:   --> %s: data_latency %lf = difftime(isystem_time_t %ld, end_time_t %ld\n", net_sta, data_latency, isystem_time_t, end_time_t);

        // check if data is before report interval -> SL server re-connect after long delay (?) // 20110414 AJL - prevent creation of duplicate events after event OT passed before report window
        if (end_time_t <= last_report_time - sliding_window_length + report_interval) {
            if (!source_packet_error_reported[source_id]) { //  have  not yet logged this srcname in this report interval
                sl_log(0, 1, "Warning: Data end time (%d) earlier than beginning of report window + report interval (%d), skipping %s (This warning is reported only once for each report interval.)\n",
                        end_time_t, last_report_time - sliding_window_length + report_interval, srcname);
                source_packet_error_reported[source_id] = 1;
            }
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            // flag as active, update latency (normally done in timedomain_processing.c->td_process_timedomain_processing())
            channelParameters[source_id].staActiveInReportInterval = 1;
            channelParameters[source_id].data_latency = data_latency;
            return;
        }
        // check if data is in future -> bad data timing (?)
        if (end_time_t > isystem_time_t + report_interval / 2) {
            if (!source_packet_error_reported[source_id]) { //  have  not yet logged this srcname in this report interval
                sl_log(0, 1, "Warning: Data end time (%d) more than %ds after current system time (%ld), skipping %s (This warning is reported only once for each report interval.)\n",
                        end_time_t, report_interval / 2, isystem_time_t, srcname);
                source_packet_error_reported[source_id] = 1;
            }
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }
        // set maximum time, should be close to current realtime
        if (end_time_t > time_max) // first set time_max to last data
            time_max = end_time_t;
        if (time_max > isystem_time_t - report_interval) // if last data within or after last report interval, have reached realtime, set time_max to realtime
            time_max = isystem_time_t;

        // parse with unpacking since we will use this data
        if (sl_msr_parse(slconn_packet->log, msrecord, &sl_msr, 1, 1) == NULL) {
            sl_log(2, 0, "sl_msr_parse() failed, skipping %s\n", srcname);
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }
        if (verbose > 2)
            printf("DEBUG: sl_msr->fsdh.num_samples %d  sl_msr->numsamples %d  sl_msr->unpackerr %d\n",
                sl_msr->fsdh.num_samples, sl_msr->numsamples, sl_msr->unpackerr);


        // increment record and sample count
        totalrecs++;
        totalsamps += sl_msr->numsamples;


        // Pack copy of the record(s) and write to BRB outstream (see record_handler)
        /*
                        MSRecord* dup_msr = msr_duplicate(sl_msr, 1);
                        outstream = outstream_RAW;
                        precords = msr_pack (dup_msr, record_handler, srcname, &psamples, 1, verbose - 2);
                        if (verbose > 2)
                        ms_log (0, "Packed %d samples into %d records, BRB.\n", psamples, precords);
         */

        // check flags
        if (verbose > 1 && sl_msr->fsdh.dq_flags != 0) {
            char b = sl_msr->fsdh.dq_flags;
            sl_log(0, 1, "Warning: Non-zero data quality flag for %s:", srcname);
            if (b & 0x01) sl_log(0, 1, " [Bit 0] Amplifier saturation detected\n");
            if (b & 0x02) sl_log(0, 1, " [Bit 1] Digitizer clipping detected\n");
            if (b & 0x04) sl_log(0, 1, " [Bit 2] Spikes detected\n");
            if (b & 0x08) sl_log(0, 1, " [Bit 3] Glitches detected\n");
            if (b & 0x10) sl_log(0, 1, " [Bit 4] Missing/padded data present\n");
            if (b & 0x20) sl_log(0, 1, " [Bit 5] Telemetry synchronization error\n");
            if (b & 0x40) sl_log(0, 1, " [Bit 6] A digital filter may be charging\n");
            if (b & 0x80) sl_log(0, 1, " [Bit 7] Time tag is questionable\n");
        }

        // process record
        if (slp_process(sl_msr, source_id, msr_net, msr_sta, msr_loc, msr_chan, details, verbose, data_latency, slconn_packet) < 0) {
            sl_log(2, 0, "slp_process() failed, skipping %s\n", srcname);
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }

        //  AJL 20100215 removed, not needed?
        //strlcat(prtloc, "X", 1);

        //strcpy(sl_msr->channel, "X");

        //msr_normalize_header(sl_msr, verbose - 2);

        // Pack the record(s) and write to HF outstream (see record_handler)
        /*
        if (!write_raw) {
                outstream = outstream_HF;
                precords = msr_pack (sl_msr, record_handler, srcname, &psamples, 1, verbose - 2);
                if (verbose > 2)
                        ms_log (0, "Packed %d samples into %d records, HF.\n", psamples, precords);
        }
         */

        sl_msr_free(&sl_msr);
        sl_msr = NULL;

        TimedomainProcessingData** de_data_list;
        int n_de_data;
        td_getTimedomainProcessingDataList(&de_data_list, &n_de_data);

        // check for completed timedomain-processing data
        for (int n = 0; n < n_de_data; n++) {
            TimedomainProcessingData* deData = de_data_list[n];
            if (deData->flag_complete_t50 == 1) {
                // flag that data has been initially processed here
                deData->flag_complete_t50 = 2;
                // save copy of deData
                /*
                TimedomainProcessingData* deDataCopy =  copy_TimedomainProcessingData(deData);
                int num_added = addTimedomainProcessingDataToDataList(deDataCopy, &local_data_list, &local_num_de_data, 1);
                if (num_added < 1) {
                        if (verbose)
                                printf("Warning: completed deData not added to deData list, already exists in list: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f\n", deDataCopy->station, deDataCopy->year, deDataCopy->month, deDataCopy->mday, deDataCopy->hour, deDataCopy->min, deDataCopy->dsec);
                        free_TimedomainProcessingData(deDataCopy);
                        continue;
                }
                 */
                // display results
                double deLevel = deData->t50 / deData->a_ref;
                if (deLevel >= 1.0)
                    strcpy(color, "RED   ");
                else if (deLevel >= 0.7)
                    strcpy(color, "Yellow");
                else if (deLevel >= 0.0)
                    strcpy(color, "green ");
                else
                    strcpy(color, "       ");
                if (verbose)
                    printf("Info: pick %d/%d, %s, elapsed %d, %.4d%.2d%.2d-%.2d:%.2d:%.4f, t50 %.2f, a_ref %.2f, t50/a_ref %.2f, S/N %.2f  %s",
                        n + 1, n_de_data, sourcename_list[deData->source_id], deData->virtual_pick_index,
                        deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec, deData->t50, deData->a_ref, deLevel, deData->a_ref / deData->sn_pick, color);
                if (deData->flag_clipped)
                    printf(" CLIPPED");
                else if (deData->flag_snr_hf_too_low)
                    printf(" S/N_HF_low");
                else if (deData->flag_snr_brb_too_low)
                    printf(" S/N_BRB_low");
                else if (deData->flag_a_ref_not_ok)
                    printf(" in_prev_coda");
                else
                    printf(" OK-------------");
                printf("\n");

            }
            if (deData->flag_complete_mwp == 1) {
                // flag that data has had initial processing here
                deData->flag_complete_mwp = 2;
                // display results
                if (verbose > 1) {
                    printf("Info: Mwp: %d/%d, %s, elapsed %d, %.4d%.2d%.2d-%.2d:%.2d:%.4f, aP %.2f, Mwp_int %.2e, Mwp_int_int %.2e, Pk1 %.2e, T1 %.2f",
                            n + 1, n_de_data, sourcename_list[deData->source_id], deData->virtual_pick_index,
                            deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec,
                            deData->mwp->amp_at_pick, deData->mwp->int_sum, deData->mwp->int_int_sum_mwp_peak,
                            deData->mwp->peak[deData->mwp->index_mag], deData->mwp->peak_dur[deData->mwp->index_mag]
                            );
                    if (deData->flag_clipped)
                        printf(" CLIPPED");
                    else if (deData->flag_snr_hf_too_low)
                        printf(" S/N_HF_low");
                    else if (deData->flag_snr_brb_too_low)
                        printf(" S/N_BRB_low");
                    else if (deData->flag_a_ref_not_ok)
                        printf(" in_prev_coda");
                    else
                        printf(" OK-------------");
                    printf("\n");
                }

            }
        }


        sl_log(0, 2, "Records: %d, Samples: %d\n", totalrecs, totalsamps);

        // create timedomain-processing report
        time_t report_latency = isystem_time_t - time_max; // 20111109 AJL
        int only_check_for_new_event =
                num_new_loc_picks >= report_trigger_min_num_picks // have enough new picks that can contribute to a location
                && time_max >= last_time_check_for_new_event + report_trigger_min_delay_time // last data is more than min delay time after last report time
                && time_max < next_report_time - report_trigger_min_delay_time // last data is earlier than min delay time before next report time
                && report_latency < report_interval / 2; // reporting latency is small
        if (next_report_time < 0 ||
                only_check_for_new_event ||
                (time_max >= next_report_time
                && report_latency < 2 * report_interval) // 20111109 AJL - skip report if processing latency large
                ) {

            /*printf("DEBUG: Mode only_check_for_new_event ->\n"
                    " num_new_loc_picks:%d >= %d,\n"
                    " time_max %ld last_time_check_for_new_event %ld\n"
                    " time_max-last_time_check_for_new_event %ld >= %d,\n"
                    " next_report_time-time_max %ld > %d,\n"
                    " report_latency %ld < %d,\n"
                    " only_check_for_new_event %d\n",
                    num_new_loc_picks, report_trigger_min_num_picks,
                    time_max, last_time_check_for_new_event,
                    time_max - last_time_check_for_new_event, report_trigger_min_delay_time,
                    next_report_time - time_max, report_trigger_min_delay_time,
                    report_latency, report_interval / 2,
                    only_check_for_new_event);*/

            long report_start_time_total = clock();

            if (1) { // DEBUG
                int nslconn;
                for (nslconn = 0; nslconn < num_slconn; nslconn++) {
                    printf("SLCONN: %s sl_state: %d\n", slconn[nslconn]->sladdr, slconn[nslconn]->stat->sl_state);
                }
            }

            tm2path(outpath, gmtime(&time_max), outpath_report);
            int cut_at_time_max = 0;
            int num_hypocenters_associated = td_writeTimedomainProcessingReport(outpath, outpath_report, time_max - sliding_window_length, time_max,
                    only_check_for_new_event, cut_at_time_max,
                    verbose, report_interval,
                    sendMail, sendMailParams, agencyId);
            if (only_check_for_new_event) {
                last_time_check_for_new_event += report_trigger_min_delay_time;
                last_time_check_for_new_event = time_max > last_time_check_for_new_event ? time_max : last_time_check_for_new_event; // 20160407 AJL - bug fix
            }
            if (!only_check_for_new_event || num_hypocenters_associated > 0) {
                sl_log(0, 1, "New report generated: %s\n", outpath_report);
                if (!only_check_for_new_event) {
                    last_report_time = time_max;
                    next_report_time = last_report_time + report_interval;
                    last_time_check_for_new_event = last_report_time;
                    // adjust next report time to an even time increment
                    if (report_interval == 60) {
                        // adjust to whole minute
                        struct tm* tm_gmt = gmtime(&next_report_time);
                        int idelay = tm_gmt->tm_sec % 60;
                        if (idelay != 0) {
                            if (idelay > 10) // adjust to next minute
                                next_report_time += (time_t) (60 - idelay);
                            else // adjust to start of this minute
                                next_report_time -= (time_t) idelay;
                        }
                    }
                }
                // reset last_data_before_report_error list so warning message above will be produced at least once in a report interval
                for (int n = 0; n < MAX_NUM_SOURCES; n++) {
                    source_packet_error_reported[n] = 0;
                }
                // reset new loc picks since last associate/locate counter
                num_new_loc_picks = 0;
                count_new_use_loc_picks_cutoff_time = time_max - report_trigger_pick_window;

                // run command
                //if (command != NULL && (difftime(time(NULL), last_report_command_time) > report_interval / 2)) {    // skip command if running too frequently (e.g. at startup)
                if (command != NULL && (difftime(time(NULL), last_report_command_time) > report_trigger_min_delay_time / 2)) { // skip command if running too frequently (e.g. at startup)
                    char* pchr = strrchr(outpath_report, '/');
                    if (pchr != NULL) {
                        strcpy(command_report, command);
                        strcat(command_report, " ");
                        strcat(command_report, outpath);
                        strcat(command_report, " ");
                        strncat(command_report, outpath_report, pchr - outpath_report);
                        strcat(command_report, " ");
                        strcat(command_report, pchr + 1);
                    } else {
                        sprintf(command_report, "%s %s %s", command, outpath, outpath_report);
                    }
                    strcat(command_report, " &"); // run in background so cannot hang trace processing
                    if (verbose > 1)
                        printf("Running command: %s\n", command_report);
                    //printf("DEBUG: Running command: %s\n", command_report);
                    int process_status = system(command_report);
                    if (verbose > 1)
                        printf("Return %d from running command: %s\n", process_status, command_report);
                    last_report_command_time = time(NULL);
                }

                printf("INFO: total report time = %.2f\n", (double) (clock() - report_start_time_total) / CLOCKS_PER_SEC);
            }
        }


    } else if (packet_type == SLKEEP) {
        sl_log(0, 2, "Keep alive packet received\n");
    } else {
        sl_log(0, 1, "%s, seq %d, Received %s blockette:\n",
                timestamp, seqnum, type[packet_type]);
        // parse without unpacking since we will not use this data
        if (sl_msr_parse(slconn_packet->log, msrecord, &sl_msr, 1, 0) == NULL) {
            sl_msr_free(&sl_msr);
            sl_msr = NULL;
            return;
        }
        if (verbose > 2)
            sl_msr_print(slconn_packet->log, sl_msr, details);
    }


} // End of packet_handler()

/**
 * slp_process:
 *
 * TODO: add doc.
 */
int slp_process(SLMSrecord* sl_msr, int source_id, char* msr_net, char* msr_sta, char* msr_loc, char* msr_chan,
        int details, int verbose, double data_latency, SLCD* slconn_packet) {

    int i;
    int *data_int = NULL;
    float *data_float = NULL;
    //double *data_double = NULL;

    int return_value = 0;

    // simple clip checking
    int flag_clipped = 0;

    if (!sl_msr)
        return (-1);

    if (details)
        sl_msr_print(slconn_packet->log, sl_msr, details);

    double deltaTime = 1.0 / sl_msr_dnomsamprate(sl_msr);

    char sampletype = 'i';
    //sampletype = sl_msr->sampletype'	// libslink has only Unpacked 32-bit data samples (?)

    // make sure data is in a float array
    if (sampletype == 'f') {
        data_float = (float *) sl_msr->datasamples;
    } else if (sampletype == 'i') {
        data_int = (int32_t *) sl_msr->datasamples;
        data_float = calloc(sl_msr->numsamples, sizeof (float));
        long idata_min = LONG_MAX;
        long idata_max = LONG_MIN;
        long idata_last = 0L;
        int iclip_counter = 0;
        for (i = 0; i < sl_msr->numsamples; i++) {
            long idata = data_int[i];
            data_float[i] = (float) data_int[i];
            // clip test
            if (idata < idata_min)
                idata_min = idata;
            if (idata > idata_max)
                idata_max = idata;
            if (idata_last != 0L && idata == idata_last && (idata == idata_min || idata == idata_max)) {
                iclip_counter++;
                if (iclip_counter > 5) {
                    flag_clipped = 1;
                    iclip_counter = 0;
                }
            } else {
                iclip_counter = 0;
            }
            idata_last = idata;
        }
    }

    // set channel id
    // set channel identifiers
    /*
    char chan_code[2];
    sprintf(chan_code, "%c", msr_net[2]);
    char inst[2];
    sprintf(inst, "?");
     * */
    // decode record start time
    struct sl_btime_s btime = sl_msr->fsdh.start_time;
    int month;
    int mday;
    if (sl_doy2md(btime.year, btime.day, &month, &mday) < 0) {
        sl_log(2, 0, "slp_process.sl_doy2md(): cannot convert doy to month/day: ch: %d:  %s %s %s %s  %.4d %.3d->%.2d%.2d  %.2d:%.2d:%2d\n",
                source_id, msr_net, msr_sta, msr_loc, msr_chan, btime.year, btime.day, month, mday, btime.hour, btime.min, btime.sec);
        return_value = -1;
        goto cleanup;
    }
    double dsec = (double) btime.sec + (double) btime.fract / 10000.0;
    // include microsecond offset in blockette 1001 if it was included. // 20120305 AJL - addded
    if (sl_msr->Blkt1001) {
        dsec += (double) sl_msr->Blkt1001->usec / 1000000.0;
        if (verbose) {
            if ((double) sl_msr->Blkt1001->usec / 1000000.0 > 0.1)
                printf("DEBUG: sl_msr->Blkt1001->usec = %d (%lfs)\n", sl_msr->Blkt1001->usec, (double) sl_msr->Blkt1001->usec / 1000000.0);
        }
    }

    // apply processing
    /*
    char network[16];
    char *pchr = strchr(msr_net, ' ');
    if (pchr != NULL)
        strncpy(network, msr_net, pchr - msr_net);
    else
        strncpy(network, msr_net, 15);
    char station[16];
    pchr = strchr(msr_sta, ' ');
    if (pchr != NULL)
        strncpy(station, msr_sta, pchr - msr_sta);
    else
        strncpy(station, msr_sta, 15);
     * */
    MSRecord* pmsrecord = msr_init(NULL);
    msr_unpack((char *) sl_msr->msrecord, SLRECSIZE, &pmsrecord, 1, verbose - 2);
    if (td_process_timedomain_processing(pmsrecord, slconn_packet->sladdr, source_id, sourcename_list[source_id], deltaTime, &data_float, (int) sl_msr->numsamples,
            msr_net, msr_sta, msr_loc, msr_chan, btime.year, btime.day, month, mday, btime.hour, btime.min, dsec,
            (int) verbose, flag_clipped, data_latency, PACKAGE,
            no_aref_level_check,
            &num_new_loc_picks, count_new_use_loc_picks_cutoff_time
            ) < 0) {
        sl_log(2, 0, "applying timedomain-processing processing: ch: %d:  %s %s %s %s  %.4d%.2d%.2d  %.2d:%.2d:%f\n",
                source_id, msr_net, msr_sta, msr_loc, msr_chan, btime.year, month, mday, btime.hour, btime.min, dsec);
        return_value = -1;
    }
    msr_free(&pmsrecord);

cleanup:

    // make sure data is in original array type
    if (sampletype == 'f') {
        sl_msr->datasamples = (void *) data_float;
    } else if (sampletype == 'i') {
        for (i = 0; i < sl_msr->numsamples; i++)
            data_int[i] = (int32_t) data_float[i];
        free(data_float);
        sl_msr->datasamples = (void *) data_int;
    }

    return (return_value);

} // End of slp_process()

/** *************************************************************************
 * parameter_proc:
 *
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ************************************************************************* **/

static int parameter_proc(int argcount, char **argvec) {

    int optind = 1;
    int error = 0;

    outpath_param = NULL;

    if (argcount <= 1)
        error++;

    // Process all command line arguments
    for (optind = 1; optind < argcount; optind++) {
        // check first for common application parameters
        if (parameter_proc_common(&optind, argvec, PACKAGE, VERSION, VERSION_DATE, usage) > 0) {
            continue;
        }
        if (strcmp(argvec[optind], "-nt") == 0) {
            slconn[num_slconn - 1]->netto = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-nd") == 0) {
            slconn[num_slconn - 1]->netdly = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-k") == 0) {
            slconn[num_slconn - 1]->keepalive = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-l") == 0) {
            strcpy(streamfile[num_slconn - 1], argvec[++optind]);
        } else if (strcmp(argvec[optind], "-s") == 0) {
            strcpy(selectors[num_slconn - 1], argvec[++optind]);
        } else if (strcmp(argvec[optind], "-S") == 0) {
            strcpy(multiselect[num_slconn - 1], argvec[++optind]);
        } else if (strcmp(argvec[optind], "-locs") == 0) {
            strcpy(location_list[num_slconn - 1], argvec[++optind]);
        } else if (strcmp(argvec[optind], "-x") == 0) {
            strcpy(statefile[num_slconn - 1], argvec[++optind]);
        } else if (strcmp(argvec[optind], "-t") == 0) {
            strcpy(begintime, argvec[++optind]);
            // check for relative begin time
            if (begintime[0] == '-') {
                sscanf(begintime, "%d", &ibackfill);
                ibackfill *= -1;
                time_t btime = time(NULL) - ibackfill;
                struct tm* tm_gmt = gmtime(&btime);
                sprintf(begintime, "%4.4d,%2.2d,%2.2d,%2.2d,%2.2d,%2.2d", tm_gmt->tm_year + 1900, tm_gmt->tm_mon + 1, tm_gmt->tm_mday, tm_gmt->tm_hour, tm_gmt->tm_min, tm_gmt->tm_sec);
            }
            if (verbose)
                printf("Begin time: %s\n", begintime);
            slconn[num_slconn - 1]->begin_time = begintime;
        } else if (strcmp(argvec[optind], "-e") == 0) {
            slconn[num_slconn - 1]->end_time = argvec[++optind];
        } else if (strcmp(argvec[optind], "-o") == 0) {
            outpath_param = argvec[++optind];
        } else if (strcmp(argvec[optind], "-c") == 0) {
            command = argvec[++optind];
        } else if (strcmp(argvec[optind], "-sendmail") == 0) {
            sendMail = 1;
            sendMailParams = argvec[++optind];
        } else if (strncmp(argvec[optind], "-v", 2) == 0) {
            verbose += strspn(&argvec[optind][1], "v");
        } else if (strncmp(argvec[optind], "-", 1) == 0) {
            fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
            exit(1);
        } // check for common application parameters
        else if (parameter_proc_common(&optind, argvec, PACKAGE, VERSION, VERSION_DATE, usage) > 0) {
            ;
        } else if (num_slconn < MAX_NUM_SLCON) {
            // initialize next slconn
            statefile[num_slconn][0] = '\0'; // state file for saving/restoring state
            slconn_terminated[num_slconn] = 0;
            slconn[num_slconn] = sl_newslcd(); // Allocate and initialize a new connection description
            slconn[num_slconn]->sladdr = argvec[optind];
            // preserve settings from last slconn
            if (num_slconn > 0) {
                slconn[num_slconn]->netto = slconn[num_slconn - 1]->netto;
                slconn[num_slconn]->netdly = slconn[num_slconn - 1]->netdly;
                slconn[num_slconn]->keepalive = slconn[num_slconn - 1]->keepalive;
                slconn[num_slconn]->begin_time = slconn[num_slconn - 1]->begin_time;
                slconn[num_slconn]->end_time = slconn[num_slconn - 1]->end_time;
                strcpy(streamfile[num_slconn], streamfile[num_slconn - 1]);
                strcpy(selectors[num_slconn], selectors[num_slconn - 1]);
                strcpy(multiselect[num_slconn], multiselect[num_slconn - 1]);
                strcpy(location_list[num_slconn], location_list[num_slconn - 1]);
                //statefile[num_slconn] = statefile[num_slconn - 1];
            } else {
                streamfile[num_slconn][0] = '\0';
                selectors[num_slconn][0] = '\0';
                multiselect[num_slconn][0] = '\0';
                location_list[num_slconn][0] = '\0';
                //statefile[num_slconn] = statefile[num_slconn - 1];
            }
            num_slconn++;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
            exit(1);
        }
    }

    // Make sure a server was specified
    if (num_slconn < 1) {
        fprintf(stderr, "No SeedLink server specified\n\n");
        fprintf(stderr, "Usage: %s [options] [host][:port]\n", PACKAGE);
        fprintf(stderr, "Try '-h' for detailed help\n");
        exit(1);
    }

    // Check if an outpath was specified
    if (outpath_param == NULL) {
        strcpy(outpath, "./seedlink_out");
        sl_log(0, 1, "No outpath (-o outpath) was specified, using %s\n", outpath);
    } else {
        strcpy(outpath, outpath_param);
    }

    // Initialize the verbosity for the sl_log function
    if (verbose > 2)
        printf("DEBUG: verbose %d\n", verbose);
    sl_loginit(verbose, NULL, NULL, NULL, NULL);

    // Report the program version
    sl_log(0, 1, "%s version: %s (%s)\n", PACKAGE, VERSION, VERSION_DATE);

    // If errors then report the usage message and quit
    if (error) {
        usage();
        exit(1);
    }

    // If verbosity is 2 or greater print detailed packet info
    if (verbose >= 2 && !details)
        details = 1;


    //
    int nslconn;
    for (nslconn = 0; nslconn < num_slconn; nslconn++) {

        // Load the stream list from a file if specified
        if (streamfile[nslconn][0] != '\0') {
            sl_read_streamlist(slconn[nslconn], streamfile[nslconn], selectors[nslconn]);
        }

        // Parse the 'multiselect' string following '-S'
        if (multiselect[nslconn][0] != '\0') {
            if (sl_parse_streamlist(slconn[nslconn], multiselect[nslconn], selectors[nslconn]) == -1)
                return -1;
        } else if (streamfile[nslconn][0] == '\0') { // No 'streams' array, assuming uni-station mode
            sl_setuniparams(slconn[nslconn], selectors[nslconn], -1, 0);
        }

        // Attempt to recover sequence numbers from state file
        if (statefile[nslconn][0] != '\0') {
            if (sl_recoverstate(slconn[nslconn], statefile[nslconn]) < 0) {
                sl_log(2, 0, "state recovery failed\n");
            }
        }
    }

    return 0;
} // End of parameter_proc()

/** *************************************************************************
 * usage:
 * Print the usage message and exit.
 ************************************************************************* **/

static void usage(void) {
    fprintf(stdout, "\nUsage: %s  [host][:port] [options] [...] [options]\n\n", PACKAGE);
    fprintf(stdout,
            " [host][:port]  Address of the SeedLink server in host:port format\n"
            "                  if host is omitted (i.e. ':18000'), localhost is assumed\n"
            "                  if :port is omitted (i.e. 'localhost'), 18000 is assumed\n"
            "                  multiple host:port SeedLink server addresses may be specified, see below for specifying SeedLink server options\n"
            "\n"
            "\n"
            " ## SeedLink server options ##\n"
            "\n"
            "                    all SeedLink server options must follow immediately after the corresponding SeedLink server address,\n"
            "                    otherwise the last options set for a previous SeedLink server address will be propagated to subsequent server addresses\n"
            " -nd delay      network re-connect delay (seconds), default 30\n"
            " -nt timeout    network timeout (seconds), re-establish connection if no\n"
            "                  data/keepalives are received in this time, default 600\n"
            " -k interval    send keepalive (heartbeat) packets this often (seconds)\n"
            " -t begintime   sets a beginning time for the initiation of data transmission (year,month,day,hour,minute,second) or (-N to begin N seconds before current time)\n"
            " -e endtime     sets an end time for windowed data transmission  (year,month,day,hour,minute,second)\n"
            " -l listfile    read a stream list from this file for multi-station mode\n"
            " -s selectors   selectors for uni-station or default for multi-station\n"
            " -S streams     select streams for multi-station (requires SeedLink >= 2.5)\n"
            "   'streams' = 'stream1[:selectors1],stream2[:selectors2],...'\n"
            "        'stream' is in NET_STA format, for example:\n"
            "        -S \"IU_KONO:BHE BHN,GE_WLF,MN_AQU:HH?.D\"\n"
            " -locs locs     comma separated list of preferred channel location codes (e.g. -locs --,00,10), use -- for none, if multiple locations are available the highest in this list will be used\n"
            " -x statefile   save/restore stream state information to this file (if used, must be specified for each host:port SeedLink server address)\n"
            "\n"
            "\n"
            " ## Application specific options ##\n"
            "\n"
            " -o outpath     path root for output report files, this path must exist\n"
            "                  if outpath is omitted './seedlink.out/' is used\n"
            " -c command     command to execute each time report files written (path root (outpath), path to results directory, and results directory name are appended after command)\n"
            "                  command is run in background, so must complete within report interval (%d seconds) to avoid collision with next invocation of command\n"
            " -sendmail      send e-mail messsage using sendmail alarm on new event location (<base web url>,<from>,<to>[,<to>]...})\n"
            "\n"
            "\n"
            , report_interval
            );
    usage_common();

} // End of usage()



#ifndef WIN32

/** *************************************************************************
 * term_handler:
 * Signal handler routine.
 ************************************************************************* **/
static void term_handler(int sig) {

    int nslconn;
    for (nslconn = 0; nslconn < num_slconn; nslconn++) {
        sl_terminate(slconn[nslconn]);
    }
    requested_sl_terminate = 1;

}
#endif



