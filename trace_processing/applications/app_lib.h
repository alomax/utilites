/*
 * File:   app_lib.h
 * Author: anthony
 *
 * Created on February 15, 2013, 8:32 AM
 */

#ifdef EXTERN_MODE
#define EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif


//#include "../geometry/geometry.h"
//#include "../alomax_matrix/alomax_matrix.h"
//#include "../matrix_statistics/matrix_statistics.h"
//#include "../picker/PickData.h"
#include "../timedomain_processing/timedomain_processing_data.h"
//#include "../timedomain_processing/ttimes.h"
#include "../timedomain_processing/location.h"
#include "../timedomain_processing/timedomain_processing.h"
#include "../timedomain_processing/timedomain_processing_report.h"

#ifndef APP_LIB_H
#define APP_LIB_H

#ifdef __cplusplus
extern "C" {
#endif

    EXTERN_TXT int details;

#define MAX_NUM_SOURCENAME_IGNORE 1024

    EXTERN_TXT int num_sourcename_ignore;
    EXTERN_TXT char sourcename_ignore[MAX_NUM_SOURCENAME_IGNORE][64];

#define DEFAULT_DATA_INPUT_COMPONENT_ACCEPT "Z" // accept Z component (only)
    EXTERN_TXT char data_input_filter_component_accept[SETTINGS_MAX_STR_LEN];
#define DEFAULT_DATA_INPUT_IGNORE_DUPLICATE_NET_STA 1 // ignore duplicate net_sta
    EXTERN_TXT int data_input_filter_ignore_duplicate_net_sta;
#define DEFAULT_DATA_INPUT_IGNORE_DUPLICATE_STA 1 // ignore duplicate net_sta
    EXTERN_TXT int data_input_filter_ignore_duplicate_sta;


#define DEFAULT_SLIDING_WINDOW_LENGTH 3600 // 1 hour
    EXTERN_TXT int sliding_window_length;
#define DEFAULT_REPORT_INTERVAL 60 // 1 minutes
    EXTERN_TXT int report_interval;
#define DEFAULT_REPORT_TRIGGER_PICK_WINDOW 600 // 10 minutes
    EXTERN_TXT int report_trigger_pick_window;
#define DEFAULT_REPORT_TRIGGER_MIN_NUM_PICKS 5 // 5 picks to trigger early assoc/location
    EXTERN_TXT int report_trigger_min_num_picks;
#define DEFAULT_REPORT_TRIGGER_MIN_DELAY_TIME 10  // 10 sec minimum delay after previous report
    EXTERN_TXT int report_trigger_min_delay_time;

    EXTERN_TXT char *outpath_param;


    void init_common_globals();
    int init_properties(char *propfile);
    int init_common_processing(char *outpath);
    void usage_common();
    int parameter_proc_common(int *poptind, char **argvec, char *PACKAGE, char *VERSION, char *VERSION_DATE, void (*usage)());
    char* tm2timestring(struct tm* timep, double secfrac, char* timestring);
    char* tm2path(char* prefix, struct tm* timep, char* pathstring);
    void app_lib_cleanup();


#ifdef __cplusplus
}
#endif

#endif /* APP_LIB_H */

