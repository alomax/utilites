/** *************************************************************************
 * miniseed_process.c
 *
 * Generic routines to perform timedomain-processing processing on Mini-SEED records from a miniseed file.
 *
 * Opens a user specified file, parses the Mini-SEED records and processes
 * each record.
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * created: 2009.02.03
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>


#ifndef WIN32
#include <signal.h>
static void term_handler(int sig);
#endif

#include <libmseed.h>

#include "../ran1/ran1.h"
#include "app_lib.h"
#include "../miniseed_utils/miniseed_utils.h"
#include "miniseed_process.h"
#include "mseedrtstream_lib.h"

#define VERSION EARLY_EST_MONITOR_VERSION
#define PACKAGE "miniseed_process"
#define VERSION_DATE EARLY_EST_MONITOR_VERSION_DATE

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#define USE_MSEEDRTSTREAM 1     // 20130404 AJL - test mseedrtstream_lib functionality
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

static short int basicsum = 1;
static int reclen = -1;

static int n_write_trace_data = 128;
static int write_traces = 0;
static int write_traces_before_processing = 0;

#ifndef USE_MSEEDRTSTREAM
#define MAX_NUM_INPUT_FILES 1024
static int num_input_files = 0;
static char *inputfile[MAX_NUM_INPUT_FILES];
#endif

static char sourcename[16384];
static char sourcename_list[MAX_NUM_SOURCES][64];
static char srcnames_duplicate[MAX_NUM_SOURCES * 64];
static double source_samplerate_list[MAX_NUM_SOURCES];

static double source_mean[MAX_NUM_SOURCES];

static int report_start = -1;
static int report_end_delay = 0;

static FILE* outstream_trace_data[MAX_NUM_SOURCES];
static FILE** outstream;

static int parameter_proc(int argcount, char **argvec);
static void usage(void);

static void record_handler(char *record, int write_reclen, void *srcname) {
    // set source id
    int n = 0;
    for (n = 0; n < num_sources_total; n++) {
        if (strcmp(sourcename_list[n], srcname) == 0) {
            break;
        }
    }
    if (n >= MAX_NUM_SOURCES || fwrite(record, write_reclen, 1, outstream[n]) != 1) {
        ms_log(2, "Error writing %s to output file.\n", (char *) srcname);
    }
}


//static TimedomainProcessingData** local_data_list = NULL;
//static int local_num_de_data = 0;

int main(int argc, char **argv) {

    MSRecord *msr = NULL;
    int64_t psamples;
    int precords;
    char srcname[50];
    char net_sta[50];
    char srcname_last[50];
    char timestr[25];

    char outname[1024];
    char outnamepath[1024];
    char outpath_report[STANDARD_STRLEN];

    int sendMail = 0; // set to 1 to test mail sending
    //char sendMailParams[] = "http://my.domain.name/early-est,me@my.domain.name,me@my.domain.name,someone.else@domain.name";
    //char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net,alomax@free.fr,anthony@alomax.net";
    char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net,alomax@free.fr";
    //char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net";

    int dataflag = 1;
    int totalrecs = 0;
    int totalsamps = 0;
    int retcode;

    int source_id = 0;
    char color[64];


#ifndef WIN32
    /* Signal handling, use POSIX calls with standardized semantics */
    struct sigaction sa;

    sa.sa_flags = SA_RESTART;
    sigemptyset(&sa.sa_mask);

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

    // initialize globals
    sliding_window_length = -1; // disables sliding window report generation if not specified in app_lib->parameter_proc_common()

    // read properties file
    char propfile[1024];
    sprintf(propfile, "%s.prop", PACKAGE);
    if (init_properties(propfile) < 0) {
        return -1;
    }

    // process given parameters (command line and parameter file)
    strcpy(agencyId, "alomax.net");
    if (parameter_proc(argc, argv) < 0)
        return -1;

#ifdef USE_MSEEDRTSTREAM
    // initialize mseedrtstream parameters
    if (processparam_init_mseedrtstream() < 0)
        return -1;
    RecordMap recmap;
    recmap.recordcnt = 0;
    recmap.first = NULL;
    recmap.last = NULL;
    if (verbose > 0)
        ms_log(1, "Info: mseedrtstream: Reading input files...\n");
    // read and process all files specified on the command line
    if (readfiles(&recmap))
        return -1;
    if (verbose > 0)
        ms_log(1, "Info: mseedrtstream: Sorting record list...\n");
    // sort the record map into time order
    if (sortrecmap(&recmap)) {
        ms_log(2, "ERROR: mseedrtstream: sorting record list\n");
        return -1;
    }
    Filelink* filelist = getFilelist();
    // Check if an outpath was specified
    if (outpath_param == NULL) {
        sprintf(outnamepath, "%s.out", filelist->infilename); // outnamepath is first read data file name
        ms_log(1, "No outpath (-o outpath) was specified, using %s\n", outnamepath);
    } else {
        strcpy(outnamepath, outpath_param);
        ms_log(1, "Outpath (-o outpath) is %s\n", outnamepath);
    }
    mkdir(outnamepath, 0755);
#else
    sprintf(outnamepath, "%s.out", inputfile[0]); // outnamepath is first read data file name
    mkdir(outnamepath, 0755);
#endif

    int RandomNumSeed = 9589;
    SRAND_FUNC(RandomNumSeed);
    if (0)
        test_rand_int();

    // initialize time domain processing
    if (init_common_processing(outnamepath) < 0) {
        fprintf(stderr, "Initialization of processing items (init_common_processing()) failed\n\n");
        fprintf(stderr, "Try '-h' for detailed help\n");
        return -1;
    }

    // initialize local statics
    strcpy(srcnames_duplicate, "\0");

    // set report options
    int do_sliding_window_reports = 1;
    if (sliding_window_length < 0) { // sliding window reports not requested
        do_sliding_window_reports = 0;
        sliding_window_length = DEFAULT_SLIDING_WINDOW_LENGTH;
        report_interval = DEFAULT_REPORT_INTERVAL;
    }

    // data time range
    time_t time_min = report_start;
    time_t time_max = -1;
    if (report_start >= 0)
        time_max = time_min + abs(report_end_delay);
    time_t last_report_time = -1;
    time_t next_report_time = -1;
    //time_t last_report_command_time = -1;

    // misc init
    strcpy(srcname_last, "$$$@@@&&&");

    // loop over input files/records
    char inputfile_name[STANDARD_STRLEN] = "";
#ifdef USE_MSEEDRTSTREAM
    // Loop through record list and process records
    Record *rec = recmap.first;
    Record *rectmp;
    off_t fpos = 0;
    while (rec != NULL) {
        // check if not using same file as last record
        int is_new_file = strcmp(inputfile_name, rec->flp->infilename);
        if (is_new_file) {
            //printf("DEBUG: clean-up! %s %s\n", inputfile_name, rec->flp->infilename);
            // cleanup
            msr_free(&msr);
            if ((retcode = ms_readmsr(&msr, NULL, reclen, &fpos, NULL, 1, dataflag, verbose - 2)) != MS_NOERROR) {
                if (retcode != MS_ENDOFFILE)
                    ms_log(2, "Cannot do cleanup call to ms_readmsr(): %s\n", ms_errorstr(retcode));
            }
        }
        strncpy(inputfile_name, rec->flp->infilename, STANDARD_STRLEN - 1);
        fpos = rec->offset;
        /* From fileutils.c -> ms_readmsr_main() :
         * If *fpos is not NULL it will be updated to reflect the file
         * position (offset from the beginning in bytes) from where the
         * returned record was read.  As a special case, if *fpos is not NULL
         * and the value it points to is less than 0 this will be interpreted
         * as a (positive) starting offset from which to begin reading data;
         * this feature does not work with packed files.
         */
        fpos *= -1;
        //printf("DEBUG: reading file: %s, offset=%ld, time=%ld->%ld", inputfile_name, rec->offset, rec->starttime, rec->endtime);
        rectmp = rec;
        rec = rec->next;
        // 20141126 AJL - added rectmp and following line to prevent memory leak (valgrind)
        free(rectmp); // rec malloc' in trace_processing/applications/mseedrtstream_lib.c
        msr_free(&msr);
        if ((retcode = ms_readmsr(&msr, inputfile_name, reclen, &fpos, NULL, 1, dataflag, verbose - 2)) == MS_NOERROR) {
            //printf("    -> : %s %s %s %s, time=%ld->%ld\n", msr->network, msr->station, msr->channel, msr->location, msr->starttime, msr_endtime(msr));
            //if ((retcode = ms_readmsr(&msr, inputfile_name, reclen, NULL, NULL, 1, dataflag, verbose - 2)) == MS_NOERROR) {
#else
    // loop over the input files and process records
    int ninputfile = 0;
    for (ninputfile = 0; ninputfile < num_input_files; ninputfile++) {
        strncpy(inputfile_name, inputfile[ninputfile], STANDARD_STRLEN - 1);
        while ((retcode = ms_readmsr(&msr, inputfile_name, reclen, NULL, NULL, 1, dataflag, verbose - 2)) == MS_NOERROR) {
#endif
            // Check source name string
            srcname[0] = '\0';
            msr_srcname(msr, srcname, 0);
            /* msr_srcname:
             *
             * Generate a source name string for a specified MSRecord in the
             * format: 'NET_STA_LOC_CHAN' or, if the quality flag is true:
             * 'NET_STA_LOC_CHAN_QUAL'.  The passed srcname must have enough room
             * for the resulting string.
             */

            // skip this record if no data samples
            if (msr->numsamples < 1) {
                msr_free(&msr);
                continue;
            }

            // skip this record if source name does not begin with a requested source name, if specified
            static char last_srcename_accept[1024];
            if (sourcename[0] != '\0') {
                if (strcmp(srcname, last_srcename_accept) != 0) {
                    //ms_log(1, "Info: Checking for sourcename: %s: index: %d\n", srcname, strstr(sourcename, srcname) - sourcename);
                    if (strstr(sourcename, srcname) == NULL) {
                        msr_free(&msr);
                        continue;
                    }
                    strcpy(last_srcename_accept, srcname);
                }
            }

            // skip this record if source name begins with sourcename_ignore
            if (num_sourcename_ignore > 0) {
                int ignore = 0;
                int n_ignore;
                for (n_ignore = 0; n_ignore < num_sourcename_ignore; n_ignore++) {
                    //printf("DEBUG: strstr [%d]: %s | %s\n", n_ignore, srcname, sourcename_ignore[n_ignore]);
                    if (strstr(srcname, sourcename_ignore[n_ignore]) == srcname) {
                        //printf("DEBUG: IGNORED!\n");
                        ignore = 1;
                    }
                }
                if (ignore) {
                    msr_free(&msr);
                    continue;
                }
            }


            // skip if not requested component orientation // 20160802 AJL - modified
            // skip if not vertical component  // 20120130 AJL - added
            //if (strrchr(msr->channel, 'Z') != msr->channel + 2) {
            //if (data_input_filter_component_accept && strrchr(msr->channel, 'Z') != msr->channel + 2) { // 20130509 AJL - can disable this check
            // 20160802 AJL - changed to comma separated list of channels to accept
            if (strrchr(data_input_filter_component_accept, msr->channel[2]) == NULL) {
                if (verbose > 1)
                    ms_log(1, "Info: Channel orientation in %s (%s) not in %s - Ignored!\n", srcname, msr->channel, data_input_filter_component_accept);
                msr_free(&msr);
                continue;
            }

            // 20130419 AJL - TEST: skip if not BHZ
            if (0 && strcmp(msr->channel, "BHZ") != 0) {
                if (verbose > 1)
                    ms_log(1, "Info: Channel in %s (%s) not BHZ - Ignored!\n", srcname, msr->channel);
                msr_free(&msr);
                continue;
            }

            int isource;

            // 20160803 AJL - following if moved outside loop
            if (strstr(srcnames_duplicate, srcname) != NULL) { // already flagged as duplicate
                isource = -1;
            } else {
                // check for sourcename in list
                for (isource = 0; isource < num_sources_total; isource++) {
                    // if perfect match, OK
                    if (strcmp(sourcename_list[isource], srcname) == 0) {
                        // confirm agreement of sample rates
                        // 20160608 AJL  if (source_samplerate_list[isource] == msr_samprate(msr)) {
                        if (fabs(source_samplerate_list[isource] - msr_samprate(msr)) / source_samplerate_list[isource] < 0.01) { // 20160608 AJL - changed to 1% accuracy
                            break;
                        } else {
                            // different sample rate !!! NOTE: this should not happen but does,
                            //    try: java -Xmx768m net.alomax.seisgram2k.SeisGram2K -seedlink "rtserve.iris.washington.edu:18000#AT_SVW2:BHZ#7200"
                            ms_log(1, "Warning: conflicting sample rates: %f in ignored %s vs. %f in active %s: this should not happen!\n",
                                    msr_samprate(msr), srcname, source_samplerate_list[isource], sourcename_list[isource]);
                            //duplicate = 4;
                            channelParameters[isource].error = channelParameters[isource].error | ERROR_DIFFERENT_SAMPLE_RATES;
                            channelParameters[isource].count_conflicting_dt++;
                            msr_free(&msr);
                            msr = NULL;
                            isource = -1;
                            break;
                        }
                    }
                    /* 20160803 AJL - moved outside this loop
                    if (strstr(srcnames_duplicate, srcname) != NULL) { // already flagged as duplicate
                        isource = -1;
                        break;
                    } */
                    // 20160802 AJL - only apply net, sta check if same orientation (last char of srcname)
                    // check on orientation is done by filter.component.accept properties file parameter, see data_input_filter_component_accept above
                    if (sourcename_list[isource][strlen(sourcename_list[isource]) - 1] == srcname[strlen(srcname) - 1]) {
                        if (data_input_filter_ignore_duplicate_net_sta) { // 20130529 AJL - can disable these check with properties file ([DataInput] filter.ignore_duplicate_net_sta)
                            // check for match of network, station only
                            strcat(strcat(strcat(strcpy(net_sta, msr->network), "_"), msr->station), "_");
                            if (strstr(sourcename_list[isource], net_sta) != NULL) {
                                if (strcmp(srcname_last, srcname) != 0)
                                    ms_log(1, "Info: Duplicate net_sta %s in %s - Ignored!\n", net_sta, srcname);
                                strcpy(srcname_last, srcname);
                                isource = -1;
                                strcat(srcnames_duplicate, "$");
                                strcat(srcnames_duplicate, srcname);
                                break;
                            }
                        }
                        if (data_input_filter_ignore_duplicate_sta) { // 20130529 AJL - can disable these check with properties file ([DataInput] filter.ignore_duplicate_sta)
                            // check for match of station only
                            //strcat(strcpy(net_sta, "_"), msr->station);
                            strcat(strcat(strcpy(net_sta, "_"), msr->station), "_");
                            if (strstr(sourcename_list[isource], net_sta) != NULL) {
                                if (strcmp(srcname_last, srcname) != 0)
                                    ms_log(1, "Info: Duplicate sta %s in %s - Ignored!\n", msr->station, srcname);
                                strcpy(srcname_last, srcname);
                                isource = -1;
                                strcat(srcnames_duplicate, "$");
                                strcat(srcnames_duplicate, srcname);
                                break;
                            }
                        }
                    }
                }
            }
            if (isource < 0) { // duplicate net_sta or sample rate mismatch
                msr_free(&msr);
                continue;
            }
            source_id = isource;
            // not found, new source name
            int is_new_source = 0;
            if (source_id == num_sources_total) {
                if (num_sources_total < MAX_NUM_SOURCES) { // sourcename not in list
                    strcpy(sourcename_list[isource], srcname);
                    source_samplerate_list[isource] = msr_samprate(msr);
                    num_sources_total++;
                    source_id = num_sources_total - 1;
                    // open output file for this source
                    if (write_traces && source_id < n_write_trace_data) {
                        sprintf(outname, "%s/%s.mseed", outnamepath, srcname);
                        outstream_trace_data[source_id] = fopen(outname, "w");
                    }
                    // flag as new source
                    // 20140401 AJL - added
                    is_new_source = 1;
                } else { // too many sources
                    ms_log(2, "Too many sources (%d max), skipping %s\n", MAX_NUM_SOURCES, srcname);
                    msr_free(&msr);
                    continue;
                }
            }

            // check time range
            double dsec;
            time_t start_time_t = hptime_t2time_t(msr_starttime(msr), &dsec);
            //if (time_min < 0 || start_time_t < time_min) {
            if (time_min < 0 // if report_start < 0 use strictly first packet start time as time min // 20110411 AJL - to prevent processing past report end, allows realistic simulation of T0, Mwpd integrals, etc.
                    || (report_start < 0 && report_end_delay > 0 && start_time_t < time_min)) {
                time_min = start_time_t;
                ms_hptime2seedtimestr(msr->starttime, timestr, 1);
                //printf("DEBUG: time_min found %s\n", timestr);
                msr_print(msr, 0);
                if (report_end_delay < 0) {
                    time_max = time_min + abs(report_end_delay); // to plot examples shortly after OT
                }
            }
            if (time_max < 0 || (report_start < 0 && report_end_delay >= 0 && start_time_t > time_max))
                time_max = start_time_t;

            // Generate a start time string
            ms_hptime2seedtimestr(msr->starttime, timestr, 1);

            // increment record and sample count
            totalrecs++;
            totalsamps += msr->samplecnt;

            // skip record if past time max // 20110411 AJL - to prevent processing past report end, allows realistic simulation of T0, Mwpd integrals, etc.
            // 20160816 AJL - change from rejecting if packet_end_time > time_max to truncating exactly to time_max
            /* time_t packet_end_time = start_time_t + (int) (msr->numsamples / msr->samprate);
            if ((report_start >= 0 || report_end_delay < 0) && packet_end_time > time_max) {
                msr_free(&msr);
                continue;
            }*/
            int msr_numsamples_use = msr->numsamples;
            if (report_start >= 0 || report_end_delay < 0) {
                if (start_time_t >= time_max) {
                    msr_free(&msr);
                    continue;
                }
                time_t packet_end_time = start_time_t + (int) (msr->numsamples / msr->samprate);
                if (packet_end_time > time_max) {
                    msr_numsamples_use = (int) (time_max - start_time_t) * msr->samprate;
                }
            }

            // write trace data
            if (write_traces && write_traces_before_processing && source_id < n_write_trace_data) {
                // Pack copy of the record(s) and write to BRB outstream (see record_handler)
                MSRecord* dup_msr = msr_duplicate(msr, 1);
                outstream = outstream_trace_data;
                precords = msr_pack(dup_msr, record_handler, srcname, &psamples, 1, verbose - 2);
                // extern int msr_pack (MSRecord *msr, void (*record_handler) (char *, int, void *), void *handlerdata, int64_t *packedsamples, flag flush, flag verbose );
                if (verbose > 2)
                    ms_log(0, "Info: Packed %d samples into %d records.\n", psamples, precords);
                msr_free(&dup_msr);
            }

            // process record
            if (msr->fsdh->dq_flags != 0) {
                if (verbose > 2)
                    printf("DEBUG: msr->fsdh_s.dq_flags %d\n", msr->fsdh->dq_flags);
            }
            // clean up location code
            char loc_str[3] = "\0\0\0";
            strcpy(loc_str, msr->location);
            if (isspace(loc_str[0]))
                loc_str[0] = '-';
            if (isspace(loc_str[1]))
                loc_str[1] = '-';
            if (strlen(loc_str) < 2)
                strcpy(loc_str, "--");
            //printf("DEBUG: %s record length = %lds\n", net_sta, hptime_t2time_t(msr_endtime(msr), &dsec) - hptime_t2time_t(msr_starttime(msr), &dsec));
            msp_process(msr, source_id, is_new_source, msr->network, msr->station, loc_str, msr->channel, msr_numsamples_use, (flag) details, (flag) verbose);
            strcat(msr->location, "X");
            //strcpy(msr->channel, "X");
            msr_normalize_header(msr, verbose - 2);

            // write trace data
            if (write_traces && !write_traces_before_processing && source_id < n_write_trace_data) {
                // Pack copy of the record(s) and write to BRB outstream (see record_handler)
                MSRecord* dup_msr = msr_duplicate(msr, 1);
                outstream = outstream_trace_data;
                precords = msr_pack(dup_msr, record_handler, srcname, &psamples, 1, verbose - 2);
                if (verbose > 2)
                    ms_log(0, "Info: Packed %d samples into %d records.\n", psamples, precords);
                msr_free(&dup_msr);
            }

            msr_free(&msr);

            // check for completed timedomain-processing data
            TimedomainProcessingData** de_data_list;
            int n_de_data;
            td_getTimedomainProcessingDataList(&de_data_list, &n_de_data);
            if (verbose) {
                for (int n = 0; n < n_de_data; n++) {
                    TimedomainProcessingData* deData = de_data_list[n];
                    if (deData->flag_complete_t50 == 1) {
                        // flag that data has had initial processing here
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
                        printf("Info: pick %d/%d, %s, elapsed %d, %.4d%.2d%.2d-%.2d:%.2d:%.4f, t50 %.2f, a_ref %.2f, t50/a_ref %.2f, S/N %.2f  %s",
                                n + 1, n_de_data, sourcename_list[deData->source_id], deData->virtual_pick_index,
                                deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec, deData->t50, deData->a_ref, deLevel, deData->a_ref / deData->sn_pick, color);
                        if (deData->flag_clipped)
                            printf(" Ignored:CLIPPED");
                        else if (deData->flag_snr_hf_too_low)
                            printf(" Ignored:S/N_HF_low");
                        else if (deData->flag_snr_brb_too_low)
                            printf(" Ignored:S/N_BRB_low");
                        else if (deData->flag_a_ref_not_ok)
                            printf(" Ignored:in_prev_coda");
                        else
                            printf(" OK-------------");
                        printf("\n");

                    }
                    if (deData->flag_complete_mwp == 1) {
                        // flag that data has been initially processed here
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
                                printf(" Ignored:CLIPPED");
                            else if (deData->flag_snr_hf_too_low)
                                printf(" Ignored:S/N_HF_low");
                            else if (deData->flag_snr_brb_too_low)
                                printf(" Ignored:S/N_BRB_low");
                            else if (deData->flag_a_ref_not_ok)
                                printf(" Ignored:in_prev_coda");
                            else
                                printf(" OK-------------");
                            printf("\n");

                        }
                    }
                }
            }

            // create timedomain-processing report

            if (do_sliding_window_reports && (next_report_time < 0 || (time_max >= next_report_time))) {

                long report_start_time_total = clock();

                tm2path(outnamepath, gmtime(&time_max), outpath_report);
                int cut_at_time_max = 0;
                int only_check_for_new_event = 0;
                td_writeTimedomainProcessingReport(outnamepath, outpath_report, time_max - sliding_window_length, time_max,
                        only_check_for_new_event, cut_at_time_max, verbose, report_interval, sendMail, sendMailParams, agencyId);
                printf("New report generated: %s\n", outpath_report);
                last_report_time = time_max;
                next_report_time = last_report_time + report_interval;
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

                // run command
                /*
                if (command != NULL && (difftime(time(NULL), last_report_command_time) > report_interval / 2)) {
                    char* pchr = strrchr(outpath_report, '/');
                    if (pchr != NULL) {
                        strcpy(command_report, command);
                                strcat(command_report, " ");
                                strcat(command_report, outnamepath);
                                strcat(command_report, " ");
                                strncat(command_report, outpath_report, pchr - outpath_report);
                                strcat(command_report, " ");
                                strcat(command_report, pchr + 1);
                    } else {
                        sprintf(command_report, "%s %s %s", command, outnamepath, outpath_report);
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
                 */

                printf("INFO: total report time = %.2f\n", (double) (clock() - report_start_time_total) / CLOCKS_PER_SEC);

            }

        }

#ifndef USE_MSEEDRTSTREAM
        if (retcode != MS_ENDOFFILE)
            ms_log(2, "Cannot read %s: %s\n", inputfile_name, ms_errorstr(retcode));
        if (basicsum)
            ms_log(1, "Info: Records: %d, Samples: %d\n", totalrecs, totalsamps);
#endif

    }

#ifdef USE_MSEEDRTSTREAM
    // cleanup
    if ((retcode = ms_readmsr(&msr, NULL, reclen, &fpos, NULL, 1, dataflag, verbose - 2)) != MS_NOERROR) {
        if (retcode != MS_ENDOFFILE)
            ms_log(2, "Cannot do cleanup call to ms_readmsr(): %s\n", ms_errorstr(retcode));
    }
    if (basicsum)
        ms_log(1, "Info: Records: %d, Samples: %d\n", totalrecs, totalsamps);
#endif


    // create timedomain-processing report

    //#define DEBUG_EVENT_PERSISTENCE 1

    if (next_report_time < 0 || last_report_time < time_max) {

        long report_start_time_total = clock();

        //time_max = time_min + 740;  // to plot examples shortly after OT

        if (abs(report_end_delay) > 0) {
            time_max = time_min + abs(report_end_delay); // to plot examples shortly after OT
        }
        int cut_at_time_max = 1;
#ifdef DEBUG_EVENT_PERSISTENCE
        time_max -= 0;
        cut_at_time_max = 0;
#endif
        int only_check_for_new_event = 0;
        if (do_sliding_window_reports) {
            tm2path(outnamepath, gmtime(&time_max), outpath_report);
        } else {
            strcpy(outpath_report, outnamepath);
        }
        td_writeTimedomainProcessingReport(outnamepath, outpath_report, time_max - sliding_window_length, time_max,
                only_check_for_new_event, cut_at_time_max,
                verbose, report_interval, sendMail, sendMailParams, agencyId);
        printf("New report generated: %s\n", outpath_report);

        printf("INFO: total report time = %.2f\n", (double) (clock() - report_start_time_total) / CLOCKS_PER_SEC);


#ifdef DEBUG_EVENT_PERSISTENCE

        report_start_time_total = clock();

        time_max += 120;
        cut_at_time_max = 1;

        if (do_sliding_window_reports) {
            tm2path(outnamepath, gmtime(&time_max), outpath_report);
        } else {
            strcpy(outpath_report, outnamepath);
        }
        td_writeTimedomainProcessingReport(outnamepath, outpath_report, time_max - sliding_window_length, time_max, cut_at_time_max,
                verbose, report_interval, sendMail, sendMailParams, agencyId);
        printf("New report generated: %s\n", outpath_report);

        printf("DEBUG: total report time = %.2f\n", (double) (clock() - report_start_time_total) / CLOCKS_PER_SEC);

#endif

    }

    // make sure everything is cleaned up
    for (int n = 0; n < num_sources_total; n++) {
        if (write_traces && n < n_write_trace_data)
            fclose(outstream_trace_data[n]);
    }
    ms_readmsr(&msr, NULL, 0, NULL, NULL, 0, 0, 0);
    td_process_timedomain_processing_cleanup();
    // remove local timedomain-processing data from list
    /*
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

} /* End of main() */

/***************************************************************************
 * msp_process:
 *
 * TODO: add doc.
 ***************************************************************************/
void msp_process(MSRecord *msr, int source_id, int is_new_source, char* msr_net, char* msr_sta, char* msr_loc, char* msr_chan, int msr_numsamples_use, flag details, flag verbose) {

    static char source_name[64];
    int i;
    int *data_int = NULL;
    double *data_double = NULL; // 20140407 AJL - added 'd' type
    float *data_float = NULL;
    //double *data_double = NULL;

    // simple clip checking
    int flag_clipped = 0;

    if (!msr)
        return;

    if (details)
        msr_print(msr, details);

    double deltaTime = 1.0 / msr->samprate;

    // get msr data in a float array, float array is newly allocated, must be freed later
    data_float = msu_getDataAsFloat(msr, &flag_clipped);
    //printf("DEBUG: msr->sampletype %c\n", msr->sampletype);
    /*if (msr->sampletype == 'f') {
        data_float = (float *) msr->datasamples;
    } else if (msr->sampletype == 'd') { // 20140407 AJL - added 'd' type
        data_double = (double *) msr->datasamples;
        data_float = calloc(msr->numsamples, sizeof (float));
        for (i = 0; i < msr->numsamples; i++) {
            data_float[i] = (float) data_double[i];
        }
    } else if (msr->sampletype == 'i') {
        data_int = (int32_t *) msr->datasamples;
        data_float = calloc(msr->numsamples, sizeof (float));
        long idata_min = LONG_MAX;
        long idata_max = LONG_MIN;
        long idata_last = 0L;
        int iclip_counter = 0;
        for (i = 0; i < msr->numsamples; i++) {
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
    }*/
    // remove mean of first packet from data
    // 20140401 AJL - added to avoid large data offsets giving very low S/N for integrated traces
    if (is_new_source) {
        source_mean[source_id] = 0.0;
        if (msr->numsamples > 0) {
            double dmean = 0.0;
            for (i = 0; i < msr->numsamples; i++) {
                dmean += data_float[i];
            }
            source_mean[source_id] = dmean / (double) msr->numsamples;
        }
    }
    for (i = 0; i < msr->numsamples; i++) {
        data_float[i] -= source_mean[source_id];
    }


    // set channel id
    // set channel identifiers
    char chan_code[2];
    sprintf(chan_code, "%c", msr->channel[2]);
    char inst[2];
    sprintf(inst, "?");
    // decode record start time
    BTime btime;
    ms_hptime2btime(msr->starttime, &btime);
    int month;
    int mday;
    ms_doy2md(btime.year, btime.day, &month, &mday);
    double dsec = (double) btime.sec + (double) btime.fract / 10000.0;

    // apply processing
    //printf("DEBUG: msr->station %s\n", msr->station);
    sprintf(source_name, "%s_%s_%s_%s", msr->network, msr->station, msr->location, msr->channel);
    double data_latency = 0.0;
    int num_new_loc_picks = 0;
    double count_new_use_loc_picks_cutoff_time = -FLT_MAX;
    int numsamples = msr_numsamples_use;
    if (td_process_timedomain_processing(msr, NULL, source_id, source_name, deltaTime, &data_float, numsamples,
            msr_net, msr_sta, msr_loc, msr_chan, btime.year, btime.day, month, mday, btime.hour, btime.min, dsec,
            (int) verbose, flag_clipped,
            data_latency, PACKAGE,
            no_aref_level_check,
            &num_new_loc_picks, count_new_use_loc_picks_cutoff_time
            ) < 0)
        ms_log(2, "Error applying timedomain-processing processing: ch: %d:  %s %s %s  %.4d%.2d%.2d  %.2d:%.2d:%f\n",
            source_id, msr->station, "?", chan_code, btime.year, month, mday, btime.hour, btime.min, dsec);


    // make sure data is in original array type
    if (msr->sampletype == 'f') {
        free(msr->datasamples);
        msr->datasamples = (void *) data_float;
    } else if (msr->sampletype == 'd') { // 20140407 AJL - added 'd' type
        data_double = (double *) msr->datasamples;
        for (i = 0; i < msr->numsamples; i++)
            data_double[i] = (double) data_float[i];
        // 20160804 AJL - bug fix  data_double[i] = (int) data_float[i];
        free(data_float);
        msr->datasamples = (void *) data_double;
    } else if (msr->sampletype == 'i') {
        data_int = (int32_t *) msr->datasamples;
        for (i = 0; i < msr->numsamples; i++)
            data_int[i] = (int) data_float[i];
        free(data_float);
        msr->datasamples = (void *) data_int;
    }


} /* End of msp_process() */

/***************************************************************************
 * get_start_time_t():
 *
 * Returns start time of an MSRecord as a time_t type
 ***************************************************************************/
time_t hptime_t2time_t(hptime_t hptime, double* pdsec) {

    BTime btime;
    ms_hptime2btime(hptime, &btime);

    struct tm tm_time;
    tm_time.tm_year = btime.year - 1900;
    int month;
    int mday;
    ms_doy2md(btime.year, btime.day, &month, &mday);
    tm_time.tm_mon = month - 1;
    tm_time.tm_mday = mday;
    tm_time.tm_hour = btime.hour;
    tm_time.tm_min = btime.min;
    tm_time.tm_sec = btime.sec;
    tm_time.tm_isdst = 0;
    tm_time.tm_gmtoff = 0;
    tm_time.tm_zone = "UTC";

    time_t t_time_t = mktime(&tm_time);
    * pdsec = (double) btime.sec + (double) btime.fract / 10000.0;

    return (t_time_t);

}

/***************************************************************************
 * parameter_proc():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int parameter_proc(int argcount, char **argvec) {
    int optind;
    int sourcename_found = 0;
    outpath_param = NULL;

    // Process all command line arguments
    for (optind = 1; optind < argcount; optind++) {
        if (strcmp(argvec[optind], "-s") == 0) {
            basicsum = 1;
        } else if (strcmp(argvec[optind], "-r") == 0) {
            reclen = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-report-start") == 0) {
            report_start = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-report-delay") == 0) {
            report_end_delay = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-n") == 0) {
            strcpy(sourcename, argvec[++optind]);
            sourcename_found = 1;
        }// check for common application parameters
        else if (parameter_proc_common(&optind, argvec, PACKAGE, VERSION, VERSION_DATE, usage) > 0) {
            ;
        } else if (strncmp(argvec[optind], "-", 1) == 0 && strlen(argvec[optind]) > 1) {
            if (processparam_flag_mseedrtstream(argcount, argvec, &optind) < 0)
                ms_log(2, "Unknown option: %s\n", argvec[optind]);
            exit(1);
        } else {
#ifdef USE_MSEEDRTSTREAM
            processparam_files_mseedrtstream(argcount, argvec, optind);
#else
            if (num_input_files < MAX_NUM_INPUT_FILES) {
                inputfile[num_input_files] = argvec[optind];
                num_input_files++;
            }
#endif
        }
    }

#ifndef USE_MSEEDRTSTREAM
    // Make sure an input file was specified
    if (num_input_files < 1) {
        ms_log(2, "No input file was specified\n");
        ms_log(1, "%s version %s\n", PACKAGE, VERSION);
        ms_log(1, "Try %s -h for usage\n", PACKAGE);
        exit(1);
    }
#endif

    // Check if a sourcename was specified
    if (!sourcename_found) {
        //ms_log (2, "No source name (-n sourcename) was specified\n\n");
        //ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
        //ms_log (1, "Try %s -h for usage\n", PACKAGE);
        sourcename[0] = '\0';
    }

    // Report the program version
    if (verbose)
        ms_log(1, "%s version: %s\n", PACKAGE, VERSION);

    return 0;

} /* End of parameter_proc() */

/***************************************************************************
 * usage():
 * Print the usage message and exit.
 ***************************************************************************/
static void usage(void) {

    fprintf(stdout, "%s version: %s\n\n", PACKAGE, VERSION);
    fprintf(stdout, "Usage: %s file [options] file1 [file2] [file3] ...\n\n", PACKAGE);
    fprintf(stdout,
            " files           Files(s) of Mini-SEED records (output written in directory <file1>.out)\n"
            "\n"
            "\n"
            " ## Application specific options ##\n"
            "\n"
            " -o outpath     path root for output report files, this path must exist\n"
            "                  if outpath is omitted './<first_mseed_file_name>.out/' is used\n"
            " -s             Print a basic summary after processing a file\n"
            " -n sourcename  Specify source names in format 'NET[_STA[_LOC[_CHAN]]]'\n"
            " -r bytes       Specify record length in bytes, required if no Blockette 1000\n"
            "\n"
            " -report-start  report start time (sec since 19700101)\n"
            " -report-delay  report end delay (sec) after -report-start, if specified, otherwise delay after earliest data\n"
            "\n"
            "\n"
            );
    usage_mseedrtstream();
    usage_common();

} /* End of usage() */



#ifndef WIN32

/***************************************************************************
 * term_handler:
 * Signal handler routine.
 ***************************************************************************/
static void
term_handler(int sig) {
    exit(0);
}
#endif

