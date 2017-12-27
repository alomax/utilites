/***************************************************************************
 * timedomain_processing.h:
 *
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * 2009->
 *
 ***************************************************************************/

#include <libmseed.h>
#include "../settings/settings.h"

#define EARLY_EST_MONITOR_NAME "Early-est - EArthquake Rapid Location sYstem with EStimation of Tsunamigenesis"
#define EARLY_EST_MONITOR_SHORT_NAME "Early-est"
#ifdef ALPHA_VERSION
#define EARLY_EST_MONITOR_VERSION "1.2.2xALPHA"   // use "N.N.NxDEV" for development version
#else
//#define EARLY_EST_MONITOR_VERSION "1.2.1"   // use "N.N.NxDEV" for development version
#define EARLY_EST_MONITOR_VERSION "1.2.2xDEV"   // use "N.N.NxDEV" for development version
#endif
#define EARLY_EST_MONITOR_VERSION_DATE "2017.12.26"
#define EARLY_EST_MONITOR_BANNER_1 "comes with ABSOLUTELY NO WARRANTY."
#define EARLY_EST_MONITOR_BANNER_2 "WARNING: DISCLAIMER: This is prototype software, it is not fully validated for use in continuous, real-time systems."
#define EARLY_EST_MONITOR_BANNER_3 "                     This software produces automatic earthquake information that has not been reviewed by a seismologist."
#define EARLY_EST_MONITOR_BANNER_4 "                     Automatically determined earthquake information may be erroneous!"
// 20140306 AJL  #define EARLY_EST_MONITOR_AGENCY "INGV, ALomax scientific"
#define EARLY_EST_MONITOR_AGENCY ""


// boolean
#define BOOLEAN_INT int
#define FALSE_INT 0
#define TRUE_INT 1

// 20141126 AJL  #define MAX_NUM_SOURCES 16384   // should be much greater than maximum number of data sources (channels) that may be processed
#define MAX_NUM_SOURCES 4096   // should be much greater than maximum number of data sources (channels) that may be processed

#define STANDARD_STRLEN 4096

#define min(x,y) ((x)<(y)?(x):(y))


#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif


EXTERN_TXT Settings *app_prop_settings;

// globals set from command line in app_lib.c and/or from properties file
//
EXTERN_TXT int verbose;
//
EXTERN_TXT char agencyId[STANDARD_STRLEN];
//
EXTERN_TXT char *geogfile;
EXTERN_TXT char *gainfile;
//
EXTERN_TXT int flag_do_grd_vel; // 2014080 AJL - added to allow testing of ground velocity based measures (e.g. fmamp)
EXTERN_TXT int flag_do_mb;
EXTERN_TXT int flag_do_mwp;
EXTERN_TXT int flag_do_mwpd;
EXTERN_TXT int flag_do_tauc;
EXTERN_TXT int flag_do_t0;
//
// 20130218 AJL - added to support microseismicity monitoring (e.g. ignore quality of tsunami discriminants)
EXTERN_TXT int no_aref_level_check;
//
// flag to enable magnitude and tsunami discriminant measures.
EXTERN_TXT int measuresEnable;
//
// flag to use hf data (high-frequency) for picking.
EXTERN_TXT int pickOnHFStream;
// flag to use raw data (broadband) for picking.
EXTERN_TXT int pickOnRawStream;

// associate/locate station corrections (used here and in timedomain_processing_report.c)
EXTERN_TXT int sta_corr_min_num_to_use;
// associate/locate station corrections (used here and in timedomain_processing_report.c)
EXTERN_TXT char sta_corr_filename[STANDARD_STRLEN];

EXTERN_TXT  int use_station_corrections;

// polarization parameters
EXTERN_TXT int polarization_enable;
#define POLARIZATION_START_DELAY_AFTER_P 0.0
EXTERN_TXT double polarization_window_start_delay_after_P;   // seconds
#define POLARIZATION_WINDOW__LENGTH_MAX 32.0
EXTERN_TXT double polarization_window_length_max;   // seconds
#define POLARIZATION_WINDOW__LENGTH_MIN 0.25
EXTERN_TXT double polarization_window_length_min;   // seconds




EXTERN_TXT int num_sources_total;
EXTERN_TXT TimedomainProcessingData** data_list;
EXTERN_TXT int num_de_data;

EXTERN_TXT HypocenterDesc** hypo_list;
EXTERN_TXT int num_hypocenters;


EXTERN_TXT ChannelResponse* chan_resp;

EXTERN_TXT ChannelParameters* channelParameters;

EXTERN_TXT ChannelParameters** sorted_chan_params_list;
EXTERN_TXT int num_sorted_chan_params;

EXTERN_TXT MSRecord*** waveform_export_miniseed_list;
EXTERN_TXT int* num_waveform_export_miniseed_list;


#define MAX_NUM_INTERNET_QUERY 10

EXTERN_TXT char* internet_gain_query[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_gain_query;
EXTERN_TXT char* internet_gain_query_host[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_gain_query_host;
EXTERN_TXT char* internet_gain_query_type[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_gain_query_type;

EXTERN_TXT char* internet_station_query[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_station_query;
EXTERN_TXT char* internet_station_query_host[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_station_query_host;
EXTERN_TXT char* internet_station_query_type[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_station_query_type;

EXTERN_TXT char* internet_timeseries_query[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_timeseries_query;
EXTERN_TXT char* internet_timeseries_query_hosturl[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_timeseries_query_hosturl;
EXTERN_TXT char* internet_timeseries_query_type[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_timeseries_query_type;
EXTERN_TXT char* internet_timeseries_query_sladdr[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_timeseries_query_sladdr;

#define STREAM_HF_BPFILTER "1-5"    // used for setting filter in IRIS_WS_TIMESERIES web service (always 4-pole), must correspons to hf filter used (e.g. filter_bp_bu_co_1_5_n_4)

#define WEB_SERIVCE_METADATA_CHECK_INTERVAL 7 // days
#define IRIS_WS_TIMESERIES 1
#define IRIS_METADATA_BASE_URL "http://www.iris.edu/mda"

// convenience structure to hold Internet query parameters

typedef struct {
    char query[STANDARD_STRLEN];
    char hosturl[STANDARD_STRLEN];
    int type;
    char sladdr[STANDARD_STRLEN];
}
InternetQuery;
EXTERN_TXT InternetQuery internetTimeseriesQueryParams[MAX_NUM_INTERNET_QUERY];
EXTERN_TXT int num_internet_timeseries;


#define MAX_NUM_PICK_CHANNEL_LIST 64
#define MAX_NUM_PROCESS_ORIENTATION_LIST 64

// convenience structure to hold FilterPicker parameters

typedef struct {
    double filterWindow; // FilterPicker5 filter window
    double longTermWindow; // FilterPicker5 long term window
    double threshold1; // FilterPicker5 threshold1
    double threshold2; // FilterPicker5 threshold1
    double tUpEvent; // FilterPicker5 tUpEvent
}
PickParams;

#ifdef USE_RABBITMQ_MESSAGING
typedef struct {
    int rmq_use_rmq;
    char rmq_hostname[STANDARD_STRLEN];
    int rmq_port;
    char rmq_exchange[STANDARD_STRLEN];
    char rmq_routingkey[STANDARD_STRLEN];
}
rmq_params;
EXTERN_TXT rmq_params rmq_parameters;
#define RMQ_USE_RMQ_DEFAULT 0
#define RMQ_HOSTNAME_DEFAULT "localhost"
#define RMQ_PORT_DEFAULT 5672
#define RMQ_EXCHANGE_DEFAULT "amq.direct"
#define RMQ_ROUTNGKEY_DEFAULT "ee_test"
#endif



double slp_dtime_curr(void);
int td_process_timedomain_processing(MSRecord* pmsrecord, char* sladdr, int source_id, char* source_name, double deltaTime, float **pdata_float_in, int numsamples,
        char* network, char* station, char* location, char* channel, int year, int doy, int month, int mday, int hour, int min, double dsec,
        int verbose, int flag_clipped, double data_latency, char *calling_routine,
        int no_aref_level_check,
        int* pnum_new_use_loc_picks, double count_new_use_loc_picks_cutoff_time
        );
int td_process_timedomain_processing_init(Settings *settings, char* outpath, char* geogfile, char* gainfile, int verbose);
int td_set_station_coordinates(int source_id, char* net, char* sta, char* loc, char* chan, int year, int month, int mday, int verbose, int icheck_ws_station_coords);
int td_set_station_corrections(int source_id, char* network, char* station, char* location, char* channel, int data_year, int data_month, int data_mday,
        int verbose, int icheck_mean_res);
int td_set_channel_gain(int source_id, char* net, char* sta, char* loc, char* chan, int year, int doy, int verbose, int icheck_ws_gain);
int td_get_gain_internet_query(int source_id, char* prtnet, char* prtsta, char* prtloc, char* prtchan, int year, int month, int day, int verbose);
int td_get_sta_coords_internet_station_query(int source_id, char* net, char* sta, char* loc, char* chan, int year, int month, int mday, int verbose);
void td_process_timedomain_processing_cleanup();
void td_process_free_timedomain_memory(int nsource);
void td_process_free_timedomain_processing_data(int source_id);
void td_getTimedomainProcessingDataList(TimedomainProcessingData*** pdata_list, int* pnum_de_data);
char *time2sring(time_t time, char* str);
double td_calculateDuration(T0Data* t0, double deltaTime, double smoothing_window);
int td_doPolarizationAnalysis(TimedomainProcessingData* deData, int ndata);

