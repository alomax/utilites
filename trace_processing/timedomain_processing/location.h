// travel times in the ak135 model
//
// iaspei-tau (http://www.iris.edu/pub/programs/iaspei-tau/)
// modified AJL, use ttimes_ajl.f
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "../geometry/geometry.h"
#include "../matrix_statistics/matrix_statistics.h"
#include "../timedomain_processing/ttimes.h"
#include "../octtree/octtree.h"



// hypocenter =====================================================================================

#define WARNING_LEVEL_STRING_LEN 64
#define NUM_LOC_SCATTER_SAMPLE 2048

typedef struct {
    double centralValue;
    double upperBound;
    double lowerBound;
    int numLevel;
}
statistic_level;

typedef struct {
    double mb; // minimum value to send alert/email message
    double mwp; // minimum value to send alert/email message
    double mwpd; // minimum value to send alert/email message
    double alarm; // minimum value to send alert/email message
    int resend_time_delay; // delay in seconds after event otime to send a single update alert/email message
}
message_trigger_theshold;

// linear regression power relation =====================================================================================
// regression values (http://en.wikipedia.org/wiki/Simple_linear_regression))

typedef struct {
    int nvalues; // number of data points used for regression
    double power; // slope of regression in log-log space
    double constant; // y-intercept of regression in log-log space
    double r_square; // square of sample correlation coefficient
    double power_dev; // 1 std deviation of power estimate
    double constant_dev; // 1 std deviation of constant estimate
}
LinRegressPower;

// hypocenter quality =====================================================================================

typedef struct {
    //  calculated during association
    double ot_variance_weight;
    double weight_sum_assoc_unique; // weight sum corrected for duplicate station/phase associations
    double amp_att_weight; // amplitude attenuation (linear regression)
    double gap_weight;
    double distanceClose_weight;
    double distanceFar_weight;
    // set later through call to setHypocenterQuality()
    double wt_count_assoc_weight;
    double errh_weight;
    double errz_weight;
    double depth_weight;
    double quality_weight;
    char quality_code[8];
}
QualityIndicators;

// structure for working storage global best hypocenter location values =====================================================================================

typedef struct {
    double node_prob;
    double node_ot_variance;
    double best_ot_mean;
    double best_ot_variance;
    double best_lat;
    double best_lon;
    double best_depth;
    double node_ot_variance_weight;
    int best_nassociated_P;
    int best_nassociated;
    double best_weight_sum;
    double best_prob;
    double effective_cell_size;
    int nCountClose;
    int numCountDistanceNodeSizeTest;
    int nCountFar;
    double distanceClose;
    double distanceNodeSizeTest;
    double distanceFar;
    // amplitude attenuation
    LinRegressPower linRegressPower;
    QualityIndicators qualityIndicators;
}
GlobalBestValues;


// NOTE: if any pointers to allocated objects are added here, then must modify creation and free functions for this structure

typedef struct {
    // location
    double weight_sum;
    double ot_std_dev;
    double prob;
    //time_t ot; // 1 sec precision
    double otime; // better than 1/10,000 sec precision
    double lat;
    double lon;
    double errh;
    double depth;
    double errz;
    double depth_step;
    int nassoc_P; // number of picks that count with weight > 0 to location ("defining" picks, may be other than P picks)
    int nassoc; // number of picks of any type associated to this hypocenter
    double dist_min; // minimum distance of associated phase counted as nassoc_P
    double dist_max; // maximum distance of associated phase counted as nassoc_P
    double gap_primary; // maximum azimuth gap
    double gap_secondary; // secondary azimuth gap - largest azimuth gap filled by a single station
    double vector_azimuth; // azimuth of sum of vectors from epicenter to stations
    double vector_distance; // 1/n_stations * distance of sum of vectors from epicenter to stations
    // statistics
    Vect3D expect; // "traditional" expectation
    Mtrx3D cov; // "traditional" covariance matrix - units always in km
    Ellipsoid3D ellipsoid; // error ellipsoid description
    Ellipse2D ellipse; // horizontal (X/Y) error ellipse description
    float *scatter_sample; // scatter sample (must be explicity allocated and freed, not handled by HypocenterDesc methods)
    int nscatter_sample; // number of valid samples in scatter_sample array
    // event persistence fields
    // global best solution oct-tree node
    OctNode global_best_oct_node;
    double global_best_critical_node_size_km;
    int n_possible_assoc_P; // cumulative number of picks that may contribute with weight > 0 to location but were attached with reassociate_only
    // warning statistics
    statistic_level t50ExLevelStatistics;
    statistic_level taucLevelStatistics;
    statistic_level tdT50ExLevelStatistics;
    char warningLevelString[WARNING_LEVEL_STRING_LEN];
    statistic_level mwpLevelStatistics;
    statistic_level mbLevelStatistics;
    statistic_level t0LevelStatistics;
    statistic_level mwpdLevelStatistics;
#ifdef USE_MWP_MO_POS_NEG
    statistic_level mwpdMoPosNegLevelStatistics;
#endif
    message_trigger_theshold messageTriggerThreshold;
    // flags
    int hyp_assoc_index; // 0...n index when associated,   event not associated if hyp_assoc_index<0
    // id
    long unique_id;
    int alert_sent_count;
    int alert_sent_time;
    // persistent fields, must be explicitly preserved when copying hypo desc
    long first_assoc_time;
    int has_valid_magnitude; // at least one mag has nreadings > MIN_NUMBER_WARNING_LEVELS_ALARM
    // amplitude attenuation
    LinRegressPower linRegressPower;
    // hypocenter quality
    QualityIndicators qualityIndicators;
    // station health
    // 20141212 AJL - added
    int nstaHasBeenActive; // number of stations for which data has been received in past
    int nstaIsActive; // number of stations for which data has been received and data_latency < report_interval
    // 20160902 AJL - added
    int loc_seq_num;    // sequence number of (re-)association/locations for this event // 20160905 AJL - added
    // 20160914 AJL - added
    int loc_type;    // type of association/location // 20160905 AJL - added
}
HypocenterDesc;
// === there are three cases for association: reassociate only (hypo preserved and fixed), relocate existing (search locally around hypo), full association location
#define LOC_TYPE_UNDEF 0
#define LOC_TYPE_FULL 1
#define LOC_TYPE_RELOC_EXISTING 2
#define LOC_TYPE_REASSOC_ONLY 3


//
void init_HypocenterDesc(HypocenterDesc *phypo);
HypocenterDesc* new_HypocenterDesc();
int addHypocenterDescToHypoList(HypocenterDesc* pnew_hypo_desc, HypocenterDesc*** phypo_list, int* pnum_hypocenters, int icheck_duplicates, HypocenterDesc* pexisting_hypo_desc,
        HypocenterDesc** phypocenter_desc__inserted);
void removeHypocenterDescFromHypoList(HypocenterDesc* hypo, HypocenterDesc*** phypo_list, int* pnum_de_data);
void free_HypocenterDescList(HypocenterDesc*** phypo_list, int* pnum_data);
int isSameEvent(HypocenterDesc* hypo1, HypocenterDesc* hypo2);

typedef struct {
    char network[32];
    char station[32];
    char location[32];
    char channel[32];
    int year;
    int doy;
    double frequency;
    double gain;    // counts/(m/s)
    int responseType;
    time_t gainfile_checked_time;
    time_t internet_gain_query_checked_time;
    int have_gain;
}
ChannelResponse;

#define ERROR_DATA_NON_CONTIGUOUS 1
#define ERROR_DT_NOT_SUPPORTED_BY_FILTER 2
#define ERROR_DIFFERENT_SAMPLE_RATES 4


#define MIN_NUM_STATION_CORRECTIONS_USE 10

typedef struct {
    int checked;
    int valid;
    double dist_min;
    double dist_max;
    int num;
    double mean;
    double std_dev;
    double poly[4];
}
StaCorrections;

typedef struct {
    char network[32];
    char station[32];
    char location[32];
    char channel[32];
    double deltaTime;
    double lat;
    double lon;
    double elev; // meters
    double azimuth; // Azimuth of the sensor in degrees from north, clockwise
    double dip; // Dip of the instrument in degrees, down from horizontal
    int inactive_duplicate;
    int staActiveInReportInterval; // station data has been received in last report interval REPORT_INTERVAL
    double data_latency;
    char data_latency_str[32];
    double last_data_end_time;
    int numData;
    int numDataAssoc;
    double qualityWeight; // quality weight, e.g. function of number of unassociated data
    time_t geogfile_checked_time;
    time_t internet_station_query_checked_time;
    int n_int_tseries; // index of internet timeseries query structure for this data
    int have_coords;
    int error; // error on processing last packet
    int count_non_contiguous;
    double level_non_contiguous;
    int count_conflicting_dt;
    double level_conflicting_dt;
    // 20160601 AJL - sta corr changed from P only to any phase
    int sta_corr_checked;
    StaCorrections *sta_corr;
    int process_this_channel_orientation; // flag indicating if this channel orientation will be processed (-1= not set, 0= no processing, 1= process)
    int channel_set[2]; // source_id's of other two sources that make a 3-comp channel set with orthogonal orientations (e.g. ZNE, Z12)
}
ChannelParameters;

typedef struct Otime_Limit {
    int data_id;
    int pick_stream; // pick stream, e.g. STREAM_HF or STREAM_RAW
    int use_for_location; // flags data to be used for location - allows picking of hf and raw data without using all picks for location
    double azimuth; // degrees CW from N
    double dist; // great-circle distance
    double dist_weight; // distance weight
    double polarization_weight; // polarization weight (polarization analysis azimuth)
    double polarization_azimuth_calc;   // calculated/predicted ray take-off azimuth
    double quality_weight; // quality weight, e.g. function of number of unassociated data
    double total_weight; // total weight for association
    double time; // time of limit
    double otime; // estimated origin time
    int polarity; // = +1 if start time of otime limit, -1 if end time
    int phase_id; // phase id
    int assoc; // currently associated
    struct Otime_Limit *pair; // paired limit
    int index; // index of limit corresponding to initilaization order
    double dist_range; // effective distance range corresponding to otime limits for this data (oct-tree search)
    double time_range_uncertainty; // time range (difference in otime limits) for this data (oct-tree search)
    // derived values
    int count_in_loc; // 1 if phase_id is a phase that counts for location origin time stacking/statistics
    int direct_P_within_dist_wt_dist_max; // flags picks to be used for distance, azimuth, attenuation filtering of locations
    // 20140121 AJL - added so that all data fields are available for association processing
    TimedomainProcessingData* deData;
    // amplitude attenuation
    double log_dist;
    double log_a_ref;
    int have_amp_atten_values;
    // station correction
    double sta_corr;

}
OtimeLimit;

// phases, travel times ===========================================================================

// 20131022 AJL - try using all picks for location, regardless of HF S/N
// 20131029 AJL - produces more events but also many more picks
#ifdef USE_USE_ALL_PICKS_FOR_LOCATION
#define USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION 1 // 20150401 AJL - try again
#else
#define USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION 0
#endif


//#define USE_AREF_NOT_OK_PICKS_FOR_LOCATION 1  // 20150402 AJL
#define USE_AREF_NOT_OK_PICKS_FOR_LOCATION 0  // 20150408 AJL

#define MAX_PICK_PERIOD_ASSOC_LOCATE 25.0

// values for location distance weighting
// global-scale values
// 20150225 AJL - changed from 30 to 20 deg to possibly help avoid distant associations for smaller events
// 20150326 AJL - gave occasional missed events, set back to 30deg
#define DISTANCE_WEIGHT_DIST_MIN 30.0
//#define DISTANCE_WEIGHT_DIST_MAX 90.0
#define DISTANCE_WEIGHT_DIST_MAX 100.0  // 20101124 AJL
//#define DISTANCE_WEIGHT_AT_DIST_MAX 0.5 // 20101228
#define DISTANCE_WEIGHT_AT_DIST_MAX 0.2 // 20130417 AJL - geometrical decay

// values for additional accept/reject hypocenter tests (grid)
//#define GAP_PRIMARY_MAX_ACCEPT 180
//#define GAP_SECONDARY_MAX_ACCEPT 270
//#define RATIO_NUM_STA_WITHIN_DIST_MIN_ACCEPT 0.6
//#define RATIO_NUM_STA_WITHIN_DIST_REJECT 0.1    // 20100609 AJL - was 0.2, may have been too high for smaller events in sparsely covered regions
//#define RATIO_NUM_STA_WITHIN_DIST_REJECT 0.333    // 20100617 AJL - was 0.1, may have been too low, did not remove flase PKP events
//#define RATIO_NUM_STA_WITHIN_DIST_REJECT 0.1    // 20100701 AJL - use 0.1 since changed to use of HF guided RAW picks
//#define RATIO_NUM_STA_WITHIN_DIST_REJECT 0.2    // 20101224 AJL
//#define DISTANCE_RATIO_CRITICAL 2.0

// values for additional accept/reject hypocenter tests (oct-tree)
// 20130503 AJL - moved to properties file  #define GAP_PRIMARY_CRITICAL 225  //20110524 AJL
// 20130503 AJL - moved to properties file  #define GAP_SECONDARY_CRITICAL 270  //20110524 AJL
#define RATIO_NUM_STA_WITHIN_CLOSE_DIST_CRITICAL 3.0  // 20111123 AJL
#define RATIO_NUM_STA_WITHIN_FAR_DIST_CRITICAL 0.3

#define MAX_NUM_DIVIDE_EXTRA 4096

#define OCTTREE_UNDEF_VALUE -VERY_SMALL_DOUBLE
#define OCTTREE_TYPE 1

Tree3D *setUpOctTree(double lat_min, double lat_max, double lat_step,
        double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step);
double octtreeGlobalAssociationLocation(int num_pass, double min_weight_sum_assoc, int max_num_nodes,
        double nominal_critical_node_size_km, double min_critical_node_size_km, double nominal_min_node_size_km,
        double gap_primary_critical, double gap_secondary_critical,
        double lat_min, double lat_max, double lat_step, double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step, int is_local_search,
        int numPhaseTypesUse, int phaseTypesUse[get_num_ttime_phases()], char channelNamesUse[get_num_ttime_phases()][128], double timeDelayUse[MAX_NUM_TTIME_PHASES][2],
        double reference_phase_ttime_error,
        TimedomainProcessingData** pdata_list, int num_de_data,
        HypocenterDesc *hypocenter,
        float **pscatter_sample, int *p_n_alloc_scatter_sample, int i_get_scatter_sample, int *pn_scatter_sample, double *p_global_max_nassociated_P_lat_lon,
        ChannelParameters * channelParameters,
        int reassociate_only, time_t time_min, time_t time_max
        );
void setOtimeLimit(OtimeLimit* otime_limit, int use_for_location, int data_id, int pick_stream, double az, double dist,
        double dist_weight, double polarization_weight, double polarization_azimuth_calc, double quality_weight, double total_weight,
        double time, double otime, int polarity, int phase_id, int index, double dist_range, double time_range_uncertainty, double depth, TimedomainProcessingData* deData,
        double sta_corr);
void addOtimeLimitToList(OtimeLimit* otimeLimit, OtimeLimit*** potime_limit_list, int* pnum_otime_limit);
void free_OtimeLimitList(OtimeLimit*** otime_limit_list, int* pnum_otime_limit);
void clear_OtimeLimitList(int* pnum_otime_limit);
void free_OtimeLimit(OtimeLimit* otime_limit);
int countStationsAvailable(double lat, double lon, double distance_max, ChannelParameters* channelParameters);
int count_in_location(int phase_id, double weight, int use_for_location);
int skipData(TimedomainProcessingData *deData, int num_pass);
int isPhaseTypeToUse(TimedomainProcessingData* deData, int phase_id, int numPhaseTypesUse, int phaseTypesUse[], char channelNamesUse[][128],
        double timeDelayUse[MAX_NUM_TTIME_PHASES][2], TimedomainProcessingData** pdata_list, int ndata_test);

double setRefOtimeDeviation(HypocenterDesc* hypo1, HypocenterDesc* hypo2);
int isSameEvent(HypocenterDesc* hypo1, HypocenterDesc* hypo2);

void freeChannelParameters(ChannelParameters* pchan_params);
void initChannelParameters(ChannelParameters* pchan_params);
int associate3CompChannelSet(ChannelParameters* channel_params, int n_sources, int source_id);
void initAssociateLocateParameters(double upweight_picks_sn_cutoff_init, double upweight_picks_dist_max, double upweight_picks_dist_full, int use_amplitude_attenuation_init);
int addChannelParametersToSortedList(ChannelParameters* pnew_chan_params, ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params);
void removeChannelParametersFromSortedList(ChannelParameters* pchan_params, ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params);
void free_ChannelParametersList(ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params);
void setHypocenterQuality(HypocenterDesc* phypo, double min_weight_sum_assoc, double critical_errh, double critical_errz);


#define MAG_NONE_PREFERRED 0
#define MAG_MB_PREFERRED 1
#define MAG_MWP_PREFERRED 2
#define MAG_MWPD_PREFERRED 3
int getPreferredMagnitude(HypocenterDesc* phypo);

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

//
EXTERN_TXT report_min_number_values_use reportMinNumberValuesUse;
EXTERN_TXT report_preferred_min_value reportPreferredMinValue;



// values (e.g. distance, azimuth =====================================================================================

typedef struct {
    double value;
    // amplitude attenuation
    double log_dist; // log of value, if needed
    double log_amp; // log of additional value, if needed
    TimedomainProcessingData* deData;
}
ValueDesc;

int addValueDescToValueList(ValueDesc* pnew_value_desc, ValueDesc*** pvalue_list, int* pnum_values);
void removeValueDescFromValueList(ValueDesc* pvalue_desc, ValueDesc*** pvalue_list, int* pnum_values);
void free_ValueDescList(ValueDesc*** pvalue_list, int* pnum_values);


