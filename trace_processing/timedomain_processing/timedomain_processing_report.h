/***************************************************************************
 * timedomain_processing_report.h:
 *
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * 2009->
 *
 ***************************************************************************/


#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

// 20150913 AJL - added to support persistence of event hypocenter information for location statistics
#define EE_TEMP_DIR "./ee_tmp"

// check for event persistence - keep previous location results for events with insufficient possible new defining picks
EXTERN_TXT int use_event_persistence;   // set to 1 if event_persistence_min_num_defining_phases > 0
EXTERN_TXT int event_persistence_min_num_defining_phases;
#define EVENT_PERSISTENCE_MIN_NUM_DEFINING_PHASES_DEFAULT 20
EXTERN_TXT double event_persistence_frac_poss_assoc_cutoff;
#define EVENT_PERSISTENCE_FRAC_POSS_ASSOC_CUTOFF_DEFAULT (1.0 / 20.0)
EXTERN_TXT double event_persistence_tt_err_fact;
#define EVENT_PERSISTENCE_TT_ERR_FACTOR_DEFAULT 2.0

// existing event association location
EXTERN_TXT int use_existing_assoc;   // set to 1 if existing_assoc_min_num_defining_phases > 0
EXTERN_TXT double existing_assoc_delay_otime_min;
#define EXISTING_EVENT_ASSOC_DELAY_OTIME_MIN_DEFAULT (60*12) // 12min
EXTERN_TXT int existing_assoc_min_num_defining_phases;
#define EXISTING_EVENT_ASSOC_MIN_NUM_DEFINING_PHASES 20

#define LOCATION_USE_REFERENCE_PHASE_TIME_ERROR_DEFAULT 1

// magnitude corrections
EXTERN_TXT int use_mwp_distance_correction;
#define LOCATION_USE_MWP_DISTANCE_CORRECTION_DEFAULT 1
EXTERN_TXT int use_mb_correction;
#define LOCATION_USE_MB_CORRECTION_DEFAULT 1
EXTERN_TXT int use_magnitude_amp_atten_check;
#define LOCATION_USE_MAGNITUDE_AMP_ATTEN_CHECK_DEFAULT 1    // requires that location.h->use_amplitude_attenuation = 1


// globals set in app_lib.c
//
EXTERN_TXT int associate_data;
EXTERN_TXT int printIgnoredData;
EXTERN_TXT int alarmNotification;

// waveform export parameters
EXTERN_TXT int waveform_export_enable;
#define WAVEFORM_EXPORT_MEMORY_SLIDING_WINDOW_LENGTH 3600
EXTERN_TXT int waveform_export_memory_sliding_window_length;   // seconds
#define WAVEFORM_EXPORT_WINDOW_START_BEFORE_P 300
EXTERN_TXT int waveform_export_window_start_before_P;   // seconds
#define WAVEFORM_EXPORT_WINDOW_END_AFTER_S 600
EXTERN_TXT int waveform_export_window_end_after_S;   // seconds
#define WAVEFORM_EXPORT_FILE_ARCHIVE_AGE_MAX (10*24*3600)     // 10 days
EXTERN_TXT int waveform_export_file_archive_age_max;   // seconds

// other properties file parameters
#define MAGNITUDE_COLORS_SHOW_DEFAULT 0     // no
EXTERN_TXT int magnitude_colors_show;
#define WARNING_COLORS_SHOW_DEFAULT 0     // no
EXTERN_TXT int warning_colors_show;
#define TSUNAMI_EVALUATION_WRITE_DEFAULT 0     // no
EXTERN_TXT int tsunami_evaluation_write;




#define LEVEL_PLOT_WINDOW_LENGTH_ASSOCIATED 1200    // 20 minutes (approx P travel time at 90 deg + 7 min)
#define LEVEL_PLOT_WINDOW_LENGTH_DEFAULT 360        // 6 minutes (long enough to allow that stations for real event are visible on map before association)
#define TIME_STEP 1

#define EPI_DIFF_LEVEL_MIN 0.0
#define EPI_DIFF_LEVEL_MAX 51.0
#define EPI_DIFF_CRITICAL_VALUE 10.0

#define USE_DISTRIBUTION_WEIGHTING 1
#define TRIM_BEFORE_WEIGHTING 1
#define T50EX_MODULO 100

#define T50EX_LEVEL_MIN 0.0
#define T50EX_LEVEL_MAX 5.199
#define T50EX_LEVEL_STEP 0.1
#define T50EX_LEVEL_MAX_PLOT (T50EX_LEVEL_MAX-2.0*T50EX_LEVEL_STEP)
//#define MIN_EPICENTRAL_DISTANCE_WARNING 10.0
#define MIN_EPICENTRAL_DISTANCE_WARNING 5.0 // 20111222 TEST AJL (5deg gives S-P=57s > TIME_DELAY_T50=55s)
#define MAX_EPICENTRAL_DISTANCE_WARNING 90.0
#define T50EX_CRITICAL_VALUE 1.0
#define T50EX_RED_CUTOFF (1.1*T50EX_CRITICAL_VALUE)
#define T50EX_YELLOW_CUTOFF (0.9*T50EX_CRITICAL_VALUE)

#define TAUC_LEVEL_MIN 0.0
#define TAUC_LEVEL_MAX 25.995
#define TAUC_LEVEL_STEP 0.5
#define TAUC_LEVEL_MAX_PLOT (TAUC_LEVEL_MAX-2.0*TAUC_LEVEL_STEP)
#define MIN_EPICENTRAL_DISTANCE_TAUC 5.0    // (5deg gives S-P=57s > TIME_DELAY_TAUC_MAX=TIME_DELAY_T50=55s)
//#define MAX_EPICENTRAL_DISTANCE_TAUC 30.0
// AJL 20100727
#define MAX_EPICENTRAL_DISTANCE_TAUC 40.0
// AJL 20100204 TauC_To_FP5simpleBRB_noLP_pk5-48_sta-dist_noTrim.ods
#define TAUC_CRITICAL_VALUE 8.0
//#define TAUC_CRITICAL_VALUE 9.0
#define TAUC_RED_CUTOFF (1.1*TAUC_CRITICAL_VALUE)
#define TAUC_YELLOW_CUTOFF (0.9*TAUC_CRITICAL_VALUE)

#define TDT50EX_LEVEL_MIN 0.0
#define TDT50EX_LEVEL_MAX 25.995
#define TDT50EX_LEVEL_STEP 0.5
#define TDT50EX_LEVEL_MAX_PLOT (TDT50EX_LEVEL_MAX-2.0*TDT50EX_LEVEL_STEP)
#define TDT50EX_CRITICAL_VALUE TAUC_CRITICAL_VALUE
#define TDT50EX_RED_CUTOFF (1.1*TDT50EX_CRITICAL_VALUE)
#define TDT50EX_YELLOW_CUTOFF (0.9*TDT50EX_CRITICAL_VALUE)
#define TDT50EX_MIN_MAIL_DEFAULT TDT50EX_YELLOW_CUTOFF
EXTERN_TXT double tdt50ex_min_mail;   // properties file parameter
#define TDT50EX_MIN_MAIL_INCREMENT TDT50EX_CRITICAL_VALUE
//#define TDT50EX_MAIL_CRITICAL_LATENCY 600.0   // wait 10min before issuing "TSUNAMI LIKELY" messages

#define MWP_LEVEL_MIN 3.0
#define MWP_LEVEL_MAX 10.199
#define MWP_LEVEL_STEP 0.1
#define MWP_LEVEL_MAX_PLOT (MWP_LEVEL_MAX-2.0*MWP_LEVEL_STEP)
// 20150316 AJL  #define MIN_EPICENTRAL_DISTANCE_MWP 5.0
#define MIN_EPICENTRAL_DISTANCE_MWP 1.0 // 20150316 AJL - TEST Mwp from 1deg !!!
#define MAX_EPICENTRAL_DISTANCE_MWP 90.0
#define MWP_CRITICAL_VALUE 7.5
#define MWP_RED_CUTOFF 7.65
#define MWP_YELLOW_CUTOFF 7.35
#define MWP_MIN_MAIL_DEFAULT 5.8
EXTERN_TXT double mwp_min_mail;   // properties file parameter
#define MWP_MIN_MAIL_INCREMENT 0.5

#define MAG_MIN_HIGHLIGHT_CUTOFF 4.95   // rounds to 5.0
//#define MAG_MAX_HIGHLIGHT_CUTOFF 6.95   // rounds to 7.0
#define MAG_MAX_HIGHLIGHT_CUTOFF 5.95   // rounds to 6.0    // 20160415 AJL - added
#define MAG_HIGH_CUTOFF 6.95   // rounds to 7.0    // 20160415 AJL - added


#define MB_LEVEL_MIN 3.0
#define MB_LEVEL_MAX 10.199
#define MB_LEVEL_STEP 0.1
#define MIN_EPICENTRAL_DISTANCE_MB 5.0
#define MAX_EPICENTRAL_DISTANCE_MB 105.0
#define MB_CRITICAL_VALUE 7.5
#define MB_RED_CUTOFF 7.65
#define MB_YELLOW_CUTOFF 7.35
//#define MB_MIN_MAIL 5.8
#define MB_MIN_MAIL_DEFAULT 5.5 // changed from 5.8 after discussion with Fabrizio Bernardi, INGV
EXTERN_TXT double mb_min_mail;   // properties file parameter
#define MB_MIN_MAIL_INCREMENT 0.5

#define T0_LEVEL_MIN 0.0
#define T0_LEVEL_MAX 259.95
#define T0_LEVEL_STEP 5.0
#ifdef USE_CLOSE_S_FOR_DURATION
#define MIN_EPICENTRAL_DISTANCE_T0 5.0 // 20111222 TEST AJL - use S duration (5deg gives S-P=57s > TIME_DELAY_T50=55s)
#define MAX_EPICENTRAL_DISTANCE_T0_CONDITIONAL 20.0 // 20120410 AJL - only use close duration for report if T50Ex > 1.0 and tauc > LIMIT
#define T0_CONDITIONAL_T50_LIMIT 1.0 // 20120410 AJL
#define T0_CONDITIONAL_TAUC_LIMIT 5.0 // 20120410 AJL
#else
//#define MIN_EPICENTRAL_DISTANCE_T0 30.0   // 20110908 TEST AJL
#define MIN_EPICENTRAL_DISTANCE_T0 15.0   // 20110908 TEST AJL
#define MAX_EPICENTRAL_DISTANCE_T0_CONDITIONAL -1.0
#endif
#define MAX_EPICENTRAL_DISTANCE_T0 90.0
#define T0_CRITICAL_VALUE 55.0
#define T0_RED_CUTOFF 60.0
#define T0_YELLOW_CUTOFF 50.0

#define MWPD_LEVEL_MIN 3.0
#define MWPD_LEVEL_MAX 10.199
#define MWPD_LEVEL_STEP 0.1
#define MIN_EPICENTRAL_DISTANCE_MWPD MIN_EPICENTRAL_DISTANCE_T0
#define MAX_EPICENTRAL_DISTANCE_MWPD MAX_EPICENTRAL_DISTANCE_T0
#define MWPD_CRITICAL_VALUE 7.5
#define MWPD_RED_CUTOFF 7.65
#define MWPD_YELLOW_CUTOFF 7.35
#define MWPD_MIN_VALUE_USE 6.95   // rounds to 7.0
#define MWPD_MIN_MAIL_DEFAULT 6.95   // rounds to 7.0
EXTERN_TXT double mwpd_min_mail;   // properties file parameter
#define MWPD_MIN_MAIL_INCREMENT 0.5
#define MIN_MWP_FOR_VALID_MWPD 6.95   // rounds to 7.0

#define PICK_PLOT_LEVEL_MIN -999.9
#define PICK_PLOT_LEVEL_MAX 999.9

typedef struct {
    int mb;
    int mwp;
    int mwpd;
    int t0;
    int tdT50Ex;
    int t50Ex;
    int tauc;
}
report_min_number_values_use;
#define MIN_NUM_VALUES_USE_DEFAULT 4

typedef struct {
    double mb;
    double mwp;
    double mwpd;
}
report_preferred_min_value;
#define PREFERRED_MIN_VALUE_MB_DEFAULT -9
#define PREFERRED_MIN_VALUE_MWP_DEFAULT 5.5
#define PREFERRED_MIN_VALUE_MWPD_DEFAULT 8.0


// 20160307 AJL - acceptable location quality levels to use for alerts, hypo display row colors, ...
#define LOC_QUALITY_ACCEPTABLE "AB"

#define ALERT_RESEND_TIME_DELAY_DEFAULT -1     // disabled
//#define ALERT_RESEND_TIME_DELAY_DEFAULT 600     // 10min
EXTERN_TXT double alert_resend_time_delay_mail;   // properties file parameter

//#define GAP_SECONDARY_MAX_ACCEPT_HYPO 270
//#define NUM_STA_WITHIN_DIST_RATIO_MIN_ACCEPT_HYPO 0.6



// ========== ASSOCIATE_LOCATE_OCTTREE ==================================
//#ifdef TTIMES_REGIONAL
#ifdef XXX
#define MAX_NUM_OCTTREE_NODES 20000
#define NOMINAL_CRITICAL_OCT_TREE_NODE_SIZE 50.0
#define MIN_CRITICAL_OCT_TREE_NODE_SIZE 1.0
#define NOMINAL_OCT_TREE_MIN_NODE_SIZE 1.0
#else
//#define MUX_NUM_OCTTREE_NODES 3500
//#define MUX_NUM_OCTTREE_NODES 5000
//#define CRITICAL_OCT_TREE_NODE_SIZE 50.0
#define MAX_NUM_OCTTREE_NODES 10000
//#define CRITICAL_OCT_TREE_NODE_SIZE 50.0
//#define OCTREE_MIN_NODE_SIZE 2.0
#define NOMINAL_CRITICAL_OCT_TREE_NODE_SIZE 50.0
#define MIN_CRITICAL_OCT_TREE_NODE_SIZE 1.0
#define NOMINAL_OCT_TREE_MIN_NODE_SIZE 5.0
#endif
// ========== END - ASSOCIATE_LOCATE_OCTTREE ==================================


#ifdef TTIMES_LONGMENSHAN_REGIONAL
#define MAX_HYPO_ARCHIVE_WINDOW (365*24*3600)     // 365 days
#else
#define MAX_HYPO_ARCHIVE_WINDOW (10*24*3600)     // 10 days
//#define MAX_HYPO_ARCHIVE_WINDOW (2*24*3600)     // testing
#endif

#define LATENCY_RED_CUTOFF 300.0
#define LATENCY_YELLOW_CUTOFF 60.0
//#define NUM_DATA_UNASSOC_STD_CUTOFF 1.0
// 20120911 AJL #define DATA_UNASSOC_CUTOFF 4
#define DATA_UNASSOC_CUTOFF 2  // 20120911 AJL
#define DATA_UNASSOC_WT_RED_CUTOFF 0.3
#define DATA_UNASSOC_WT_YELLOW_CUTOFF 0.6
#define LEVEL_NON_CONTIGUOUS_RED_CUTOFF 0.5
#define LEVEL_NON_CONTIGUOUS_YELLOW_CUTOFF 0.1
#define LEVEL_CONFLICTING_DT_RED_CUTOFF 0.5
#define LEVEL_CONFLICTING_DT_YELLOW_CUTOFF 0.1

// simpleGlobalLocation
/*
static double ***weight_count_grid_coarse = NULL;
static int nlat_alloc_coarse, nlon_alloc_coarse;
static double ***weight_count_grid = NULL;
static int nlat_alloc, nlon_alloc;
*/


int td_timedomainProcessingReport_init(char* outnameroot_archive);
void td_timedomainProcessingReport_cleanup();
int td_writeTimedomainProcessingReport(char* outnameroot_archive, char* outnameroot, time_t time_min, time_t time_max,
        int check_for_new_event, int cut_at_time_max,
        int verbose, int report_interval, int sendMail, char *sendMailParams, char *agencyId);
