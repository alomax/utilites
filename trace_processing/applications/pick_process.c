/** *************************************************************************
 * pick_process.c
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
#include "pick_process.h"
#include "mseedrtstream_lib.h"

#define VERSION EARLY_EST_MONITOR_VERSION
#define PACKAGE "pick_process"
#define VERSION_DATE EARLY_EST_MONITOR_VERSION_DATE

static short int basicsum = 1;

#define MAX_NUM_INPUT_FILES 1024
static int num_input_files = 0;
static char *inputfile[MAX_NUM_INPUT_FILES];

static char sourcename[16384];
static char sourcename_list[MAX_NUM_SOURCES][64];
static char srcnames_duplicate[MAX_NUM_SOURCES * 64];

static char phases_count_in_location[1024];

static int report_start = -1;
static int report_end_delay = 0;

static double pick_error_maximum = FLT_MAX;


// forward function declarations
int parameter_proc(int argcount, char **argvec);
void usage(void);
int check_phases_count_in_location();

int main(int argc, char **argv) {

    char srcname[50];
    char net_sta[50];
    char srcname_last[50];

    char outnamepath[1024];
    char outpath_report[STANDARD_STRLEN];

    int sendMail = 0; // set to 1 to test mail sending
    //char sendMailParams[] = "http://my.domain.name/early-est,me@my.domain.name,me@my.domain.name,someone.else@domain.name";
    //char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net,alomax@free.fr,anthony@alomax.net";
    char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net,alomax@free.fr";
    //char sendMailParams[] = "http://early-est.alomax.net,anthony@alomax.net";
    char agencyId[] = "alomax.net";

    int source_id = 0;


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
    flag_do_grd_vel = 0;
    // initialize globals
    sliding_window_length = -1; // disables sliding window report generation if not specified in app_lib->parameter_proc_common()

    // read properties file
    char propfile[1024];
    sprintf(propfile, "%s.prop", PACKAGE);
    if (init_properties(propfile) < 0) {
        return -1;
    }

    // process given parameters (command line and parameter file)
    strcpy(phases_count_in_location, "");
    if (parameter_proc(argc, argv) < 0)
        return -1;
    // check for phases_count_in_location
    check_phases_count_in_location();

    strcpy(outnamepath, outpath_param);
    mkdir(outnamepath, 0755);

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
    // loop over the input files and process records
    TimedomainProcessingData *deData;
    int ninputfile = 0;
    for (ninputfile = 0; ninputfile < num_input_files; ninputfile++) {
        strncpy(inputfile_name, inputfile[ninputfile], STANDARD_STRLEN - 1);
        FILE *fp_pick = fopen(inputfile_name, "r");
        if (fp_pick == NULL) {
            fprintf(stderr, "ERROR: Cannot open pick file: %s\n", inputfile_name);
            continue;
        }
        int pick_format = PICK_FORMAT_NLL;
        while ((deData = get_next_pick(fp_pick, pick_format, verbose - 2)) != NULL) {

            // generate source name string
            sprintf(srcname, "%s_%s_%s_%s", deData->network, deData->station, deData->location, deData->channel);

            // skip this pick if source name does not begin with a requested source name, if specified
            static char last_srcename_accept[1024];
            if (sourcename[0] != '\0') {
                if (strcmp(srcname, last_srcename_accept) != 0) {
                    //printf("Info: Checking for sourcename: %s: index: %d\n", srcname, strstr(sourcename, srcname) - sourcename);
                    if (strstr(sourcename, srcname) == NULL) {
                        free_TimedomainProcessingData(deData);
                        continue;
                    }
                    strcpy(last_srcename_accept, srcname);
                }
            }

            // skip this pick if source name begins with sourcename_ignore
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
                    free_TimedomainProcessingData(deData);
                    continue;
                }
            }


            // skip if not requested component orientation // 20160802 AJL - modified
            if (strlen(deData->channel) >= 3) {
                if (strrchr(data_input_filter_component_accept, deData->channel[2]) == NULL) {
                    if (verbose > 1)
                        printf("Info: Channel orientation in %s (%s) not in %s - Ignored!\n", srcname, deData->channel, data_input_filter_component_accept);
                    free_TimedomainProcessingData(deData);
                    continue;
                }
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
                        break;
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
                            strcat(strcat(strcat(strcpy(net_sta, deData->network), "_"), deData->station), "_");
                            if (strstr(sourcename_list[isource], net_sta) != NULL) {
                                if (strcmp(srcname_last, srcname) != 0)
                                    printf("Info: Duplicate net_sta %s in %s - Ignored!\n", net_sta, srcname);
                                strcpy(srcname_last, srcname);
                                isource = -1;
                                strcat(srcnames_duplicate, "$");
                                strcat(srcnames_duplicate, srcname);
                                break;
                            }
                        }
                        if (data_input_filter_ignore_duplicate_sta) { // 20130529 AJL - can disable these check with properties file ([DataInput] filter.ignore_duplicate_sta)
                            // check for match of station only
                            //strcat(strcpy(net_sta, "_"), deData->station);
                            strcat(strcat(strcpy(net_sta, "_"), deData->station), "_");
                            if (strstr(sourcename_list[isource], net_sta) != NULL) {
                                if (strcmp(srcname_last, srcname) != 0)
                                    printf("Info: Duplicate sta %s in %s - Ignored!\n", deData->station, srcname);
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
                free_TimedomainProcessingData(deData);
                continue;
            }
            source_id = isource;
            // not found, new source name
            int is_new_source = 0;
            if (source_id == num_sources_total) {
                if (num_sources_total < MAX_NUM_SOURCES) { // sourcename not in list
                    strcpy(sourcename_list[isource], srcname);
                    num_sources_total++;
                    source_id = num_sources_total - 1;
                    // open output file for this source
                    // flag as new source
                    // 20140401 AJL - added
                    is_new_source = 1;
                } else { // too many sources
                    printf("ERROR: Too many sources (%d max), skipping %s\n", MAX_NUM_SOURCES, srcname);
                    free_TimedomainProcessingData(deData);
                    continue;
                }
            }

            // check time range
            time_t start_time_t = deData->t_time_t;
            if (time_min < 0 || start_time_t < time_min) {
                time_min = start_time_t;
            }
            if (start_time_t > time_max) {
                time_max = start_time_t;
            }

            // process pick
            if (pick_process(deData, source_id, is_new_source, verbose) < 0) {
                free_TimedomainProcessingData(deData);
                continue;
            }

            int num_added = addTimedomainProcessingDataToDataList(deData, &data_list, &num_de_data, 1);
            if (num_added < 1) {
                if (verbose)
                    printf("Warning: pick: new pick not added to deData list, nearby pick already exists in list: %s, %.4d%.2d%.2d-%.2d:%.2d:%.4f\n",
                        srcname, deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec);
                free_TimedomainProcessingData(deData);
                continue;
            }


            // create timedomain-processing report

            if (do_sliding_window_reports && (next_report_time < 0 || (time_max >= next_report_time))) {

                long report_start_time_total = clock();

                tm2path(outnamepath, gmtime(&time_max), outpath_report);
                int cut_at_time_max = 0;
                int only_check_for_new_event = 0;
                td_writeTimedomainProcessingReport(outnamepath, outpath_report, time_max - sliding_window_length, time_max,
                        only_check_for_new_event, cut_at_time_max, verbose, report_interval, sendMail, sendMailParams, agencyId);
                printf("New sliding window report generated: %s\n", outpath_report);
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

                printf("INFO: total report time = %.2f\n", (double) (clock() - report_start_time_total) / CLOCKS_PER_SEC);

            }

        }

        fclose(fp_pick);

    }


    // create timedomain-processing report

    if (next_report_time < 0 || last_report_time < time_max) {

        long report_start_time_total = clock();

        //time_max = time_min + 740;  // to plot examples shortly after OT

        if (abs(report_end_delay) > 0) {
            time_max = time_min + abs(report_end_delay); // to plot examples shortly after OT
        }
        int cut_at_time_max = 1;
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


    }

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

int check_phases_count_in_location() {

    char *str_pos = strtok(phases_count_in_location, ",");
    int nphases = 0;
    while (str_pos != NULL) {
        add_phase_type_flag_to_phase_name(str_pos, COUNT_IN_LOCATION);
        printf("INFO: phases_count_in_location: added: %s\n", str_pos);
        nphases++;
        str_pos = strtok(NULL, ",");
    }

    return(nphases);

}

/***************************************************************************
 * pick_process:
 *
 * TODO: add doc.
 ***************************************************************************/
int pick_process(TimedomainProcessingData* deData, int source_id, int is_new_source, int verbose) {

    // check have station coordinates set station as active
    //printf("DEBUG: Check ws station: %s %s  time(NULL) %ld >? channelParameters[source_id].internet_station_query_checked_time %ld (diff=%f)\n",
    //        network, station, time(NULL), channelParameters[source_id].internet_station_query_checked_time,
    //        difftime(time(NULL), channelParameters[source_id].internet_station_query_checked_time) / (24 * 3600));
    int icheck_ws_station_coords = difftime(time(NULL), channelParameters[source_id].internet_station_query_checked_time) > (WEB_SERIVCE_METADATA_CHECK_INTERVAL * 24 * 3600);
    if (!channelParameters[source_id].have_coords || icheck_ws_station_coords) {
        td_set_station_coordinates(source_id, deData->network, deData->station, deData->location, deData->channel, deData->year, deData->month, deData->mday, verbose, icheck_ws_station_coords);
        if (!channelParameters[source_id].have_coords) {
            if (verbose)
                printf("ERROR: sensor coordinates not found: %s %s, please add sensor data to file: %s\n", deData->network, deData->station, geogfile);
            return (-1);
        }
        // find and store references to other channel orientations for this net/sta/chan
        // TODO: activate this?  associate3CompChannelSet(channelParameters, num_sources_total, source_id);
    }

    if (use_station_corrections) {

        // check for station corrections
        int icheck_sta_corr = 0; // TODO: maybe add periodic check for updated station residuals, not necessary if add automatic update of residuals for each event located
        if (!channelParameters[source_id].sta_corr_checked || icheck_sta_corr) {
            td_set_station_corrections(source_id, deData->network, deData->station, deData->location, deData->channel, deData->year, deData->month, deData->mday, verbose, icheck_sta_corr);
        }
    }


    // stuff to add or ignore (from timedomain_processing/timedomain_processing_data.c->new_TimedomainProcessingData())
    //deData->sladdr = sladdr;
    //deData->n_int_tseries = n_int_tseries;
    deData->source_id = source_id;
    deData->lat = channelParameters[source_id].lat;
    deData->lon = channelParameters[source_id].lon;
    deData->elev = channelParameters[source_id].elev;
    deData->azimuth = channelParameters[source_id].azimuth;
    deData->dip = channelParameters[source_id].dip;
    //*(deData->pickData) = *pickData; // make local copy of pickData
    deData->pick_stream = STREAM_NULL;
    //deData->virtual_pick_index = init_ellapsed_index_count;

    deData->use_for_location = 1;
    deData->can_use_as_location_P = 1;
    deData->use_for_location_twin_data = NULL;
    deData->flag_complete_t50 = 1;

    // check pick err
    if (deData->pick_error > pick_error_maximum) {
        deData->pick_error = pick_error_maximum;
    }


    return (0);

} /* End of pick_process() */

/***************************************************************************
 * parameter_proc():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
int parameter_proc(int argcount, char **argvec) {
    int optind;
    int sourcename_found = 0;
    outpath_param = NULL;

    // Process all command line arguments
    for (optind = 1; optind < argcount; optind++) {
        if (strcmp(argvec[optind], "-s") == 0) {
            basicsum = 1;
        } else if (strcmp(argvec[optind], "-report-start") == 0) {
            report_start = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-report-delay") == 0) {
            report_end_delay = atoi(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-pick-err-max") == 0) {
            pick_error_maximum = atof(argvec[++optind]);
        } else if (strcmp(argvec[optind], "-phases-count-in-location") == 0) {
            strcpy(phases_count_in_location, argvec[++optind]);
        } else if (strcmp(argvec[optind], "-n") == 0) {
            strcpy(sourcename, argvec[++optind]);
            sourcename_found = 1;
        }// check for common application parameters
        else if (parameter_proc_common(&optind, argvec, PACKAGE, VERSION, VERSION_DATE, usage) > 0) {
            ;
        } else if (strncmp(argvec[optind], "-", 1) == 0 && strlen(argvec[optind]) > 1) {
            printf("ERROR: Unknown option: %s\n", argvec[optind]);
            exit(1);
        } else {
            if (num_input_files < MAX_NUM_INPUT_FILES) {
                inputfile[num_input_files] = argvec[optind];
                num_input_files++;
            }
        }
    }

    // Make sure an input file was specified
    if (num_input_files < 1) {
        printf("ERROR: No input file was specified\n");
        printf("%s version %s\n", PACKAGE, VERSION);
        printf("Try %s -h for usage\n", PACKAGE);
        exit(1);
    }

    // Check if a sourcename was specified
    if (!sourcename_found) {
        sourcename[0] = '\0';
    }

    // Report the program version
    if (verbose)
        printf("%s version: %s\n", PACKAGE, VERSION);

    return 0;

} /* End of parameter_proc() */

/***************************************************************************
 * usage():
 * Print the usage message and exit.
 ***************************************************************************/
void usage(void) {

    fprintf(stdout, "%s version: %s\n\n", PACKAGE, VERSION);
    fprintf(stdout, "Usage: %s file [options] file1 [file2] [file3] ...\n\n", PACKAGE);
    fprintf(stdout,
            " files           Pick files(s) in format: NonLinLoc, ...\n"
            "\n"
            "\n"
            " ## Application specific options ##\n"
            "\n"
            " -o outpath     path root for output report files, this path must exist\n"
            "                  if outpath is omitted './<first_pick_file_name>.out/' is used\n"
            " -s             Print a basic summary after processing a file\n"
            " -n sourcename  Specify source names in format 'NET[_STA[_LOC[_CHAN]]]'\n"
            "\n"
            " -report-start  report start time (sec since 19700101)\n"
            " -report-delay  report end delay (sec) after -report-start, if specified, otherwise delay after earliest data\n"
            " -pick-err_max  maximum pick error to use, error will be reduced to this value for picks with larger error\n"
            " -phases-count-in-location  comma separated list of phases to count in location (e.g. S,Sn,Sg) \n"
            "\n"
            "\n"
            );
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

