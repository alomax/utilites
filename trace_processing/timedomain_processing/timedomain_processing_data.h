

/** timedomain-processing data structure/class */


#include <time.h>

#include <libmseed.h>

#include "../picker/PickData.h"


/*
enum AmpUnits { NO_AMP_UNITS, DATA_AMP_UNITS, CHAR_FUNCT_AMP_UNITS, INDEP_VAR_UNITS };

#define POLARITY_POS 		1
#define POLARITY_UNKNOWN 	0
#define POLARITY_NEG 		-1
 */

// test of use of Mwpd positive and negative moment ratios as depth & tsunamigenic indicator
// 20161226 AJL - added
// TEST ONLY, DOES NOT HELP
//#define USE_MWP_MO_POS_NEG


// amplitude attenuation
#define AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO 5.0
#define AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO (1.0 / AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO)

// use Mwp integral for first motion
#define USE_MWP_INT_FOR_FIRST_MOTION 1 // 20120312 TEST AJL - use Mwp integral for first motion
#define MWP_POLARITY_PULSE_FRACTION 0.5     // fraction of peak Mwp pulse duration to use for estimating Mwp first motion polarity
#define MWP_POLARITY_DELAY_MAX 20.0     // sec, maximum delay after P for estimating Mwp first motion polarity
#define MWP_POLARITY_PULSE_WIDTH_MIN 1.0   // sec, min Mwp pulse width for estimating Mwp first motion polarity

// use S duration
#define USE_CLOSE_S_FOR_DURATION 1 // 20111222 TEST AJL - use S duration
#define MINIMUN_DURATION 0.2   // highest frequency in HF stream (filter_bp_bu_co_1_5_n_4)

#define VALUE_NOT_SET -999.9

#define TAUC_INVALID -9.0   // 20120106 AJL - added
#define MWP_INVALID -9.0
#define MWPD_INVALID -9.0
#define MB_INVALID -9.0
#define T0_INVALID -1.0
#define GRD_MOT_INVALID -1.0

#define MB_MODE_mb_DISP 0   // classic mb, transfer function is WWSSN_SP_DISP, gives velocity->velocity, integrate to displacement as needed for mb
#define MB_MODE_mb_VEL 1    // classic mb, transfer function is WWSSN_SP_VEL, gives velocity->disp, so already have displacement as needed for mb
#define MB_MODE_mB 2        // mB of 2008__Bormann_Saul w/ WWSSN-SP response, transfer function is WWSSN_SP_DISP, gives velocity->velocity, so alreay have velocity as needed for mB
#define MB_MODE MB_MODE_mB

#define START_ANALYSIS_BEFORE_P_BRB_HP 1.0  // set zero level and start analyses before P to avoid setting zero level after start of strong LP P signal.

#define MAX_GRD_VEL_DUR 120.0   // maximum time to save BRB HP ground velocity data

#define MAX_MWP_DUR 120.0
#define MIN_MWPD_DUR 20.0
#define MAX_MWPD_DUR 1200.0     // 20 min - assumed longer than any real earthquake duration
//#define MAX_MB_DUR 120.0
#define MAX_MB_DUR 30.0    // 20110627 AJL - 30s used by GFZ for mb

// 20110601 AJL - separated from above and changed to shorter shift
//#define START_ANALYSIS_BEFORE_P_MB 1.0  // set zero level and start integration before P to avoid setting zero level after start of strong LP P signal.
// 20110630 AJL - changed to 0 - caused large integration offsets for mb
#define START_ANALYSIS_BEFORE_P_MB 0.0  // set zero level and start integration before P to avoid setting zero level after start of strong LP P signal.

#define NUMBER_ASSOCIATE_IGNORE -999

#define STREAM_NULL -1
#define STREAM_HF 0
#define STREAM_RAW 1
#define STREAM_NULL_NAME "NULL"
#define STREAM_HF_NAME "HF"
#define STREAM_RAW_NAME "BRB"


// 25 sec window (P -> P+25s)
#define TIME_DELAY_A_REF 25.0
#define SMOOTHING_WINDOW_HALF_WIDTH_A_REF 12.5
// 10 sec window (P+50s -> P+60s)
// 20090416 AJL test
//#define TIME_DELAY_T50 60.0
#define TIME_DELAY_T50 55.0
#define SMOOTHING_WINDOW_HALF_WIDTH_T50 5.0
#define A_REF_OK_RATIO 0.5
#define A_REF_VALID_TIME 120.0
// AJL20090519 - important change !!!
#define SIGNAL_TO_NOISE_RATIO_HF_MIN 3.0
// AJL20131021 - TEST !!!
//#define SIGNAL_TO_NOISE_RATIO_HF_MIN 2.0
//#define SIGNAL_TO_NOISE_RATIO_HF_MIN 5.0
// END - AJL20090519
//#define SIGNAL_TO_NOISE_RATIO_HF_MIN 5.0
//#define SIGNAL_TO_NOISE_RATIO_HF_MIN 10.0
#define A_REF_WINDOW_MIN_TO_MAX_AMP_RATIO_MIN 0.3333

// 20180823 AJL - Bug fix: SHIFT_BEFORE_P values were too large, does not capture increase in amplitude just before pick when pick in LP noise of previous event, likely causing false, high magnitudes on strong LP arrivals from preceeding large events
/*
 * In following, Fiji has false, too high Mwp7.3 (note only 4 Mwp readings!)
 1534995679202 	25	220	175	3.7	41	46	-0.91	1.3  	2018.08.23-03:41:18	-18.19	-178.13	   15	559	   17	A	0.0 	 56	4.4 	 2	0.0 	5.1 	 96	7.3 	 4	1 	 55	-9.0 	 2		 Fiji Islands Region <- False Mwp7.3
 1534995310291 	11	345	290	12.6	78	84	-0.71	1.2  	2018.08.23-03:35:22	51.59	-177.93	   13	85	   22	A	0.7 	 198	4.3 	 35	2.9 	5.4 	 249	6.4 	 205	40 	 192	6.5 	 193		 Andreanof Islands, Aleutian Is.
*/
// brb
// Mwp brb hp int  // 20120612 AJL - added
// 20151117 AJL #define SHIFT_BEFORE_P_BRB_HP_INT_SN 10.0  // S/N N window end
// 20151117 AJL #define TIME_DELAY_BRB_HP_INT_RMS 200.0    // S/N S & N window
// 20180823 AJL #define SHIFT_BEFORE_P_BRB_HP_INT_SN 20.0  // S/N N window end
#define SHIFT_BEFORE_P_BRB_HP_INT_SN 5.0  // S/N N window end before P  // 20180823 AJL - changed to 1/4 S/N window TIME_DELAY_*
#define TIME_DELAY_BRB_HP_INT_RMS 20.0    // S/N S & N window
#define SMOOTHING_WINDOW_HALF_WIDTH_BRB_HP_INT_RMS (TIME_DELAY_BRB_HP_INT_RMS/2.0)
#define SIGNAL_TO_NOISE_RATIO_BRB_HP_INT_MIN 3.0
// brb hp
// 20151117 AJL #define SHIFT_BEFORE_P_BRB_HP_SN 10.0  // S/N N window end
// 20151117 AJL #define TIME_DELAY_BRB_HP_RMS 50.0    // S/N S & N window
// 20180823 AJL #define SHIFT_BEFORE_P_BRB_HP_SN 10.0  // S/N N window end
#define SHIFT_BEFORE_P_BRB_HP_SN 2.5  // S/N N window end before P  // 20180823 AJL - changed to 1/4 S/N window TIME_DELAY_*
#define TIME_DELAY_BRB_HP_RMS 10.0    // S/N S & N window
#define SMOOTHING_WINDOW_HALF_WIDTH_BRB_HP_RMS (TIME_DELAY_BRB_HP_RMS/2.0)
//#define SIGNAL_TO_NOISE_RATIO_BRB_MIN 3.0
//#define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN 2.5  // 20100412 AJL
#define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN 3.0  // 20101229 AJL
// mB bp
// 20151117 AJL #define SHIFT_BEFORE_P_BRB_BP_SN 4.0  // S/N N window end
// 20151117 AJL #define TIME_DELAY_BRB_BP_RMS 20.0    // S/N S & N window
// 20180823 AJL #define SHIFT_BEFORE_P_BRB_BP_SN 10.0  // S/N N window end
#define SHIFT_BEFORE_P_BRB_BP_SN 2.5  // S/N N window end before P  // 20180823 AJL - changed to 1/4 S/N window TIME_DELAY_*
#define TIME_DELAY_BRB_BP_RMS 10.0    // S/N S & N window
//#define TIME_DELAY_BRB_BP_RMS 50.0  // 20110610 AJL
#define SMOOTHING_WINDOW_HALF_WIDTH_BRB_BP_RMS (TIME_DELAY_BRB_BP_RMS/2.0)
//#define SMOOTHING_WINDOW_HALF_WIDTH_BRB_BP_RMS 25.0  // 20110610 AJL
#define SIGNAL_TO_NOISE_RATIO_BRB_BP_MIN 1.5  // 20101226 AJL
//#define SIGNAL_TO_NOISE_RATIO_BRB_BP_MIN 3.0  // 20110628 AJL
// t0 Mwpd
//#define SIGNAL_TO_NOISE_RATIO_T0_END 1.5
#define SIGNAL_TO_NOISE_RATIO_T0_END 3.0    // 20110117 AJL
#define SIGNAL_TO_PEAK_RATIO_T0_END 0.025
// tauc
#define INSTANT_PERIOD_WINDOWED_WINDOW_WIDTH  5.0
//#define INSTANT_PERIOD_WINDOWED_WINDOW_WIDTH  (TIME_DELAY_T50-5.0)      // 20130304 AJL - TEST!    -> tsunami_warning_008_Td50s_15m.ods
//#define TIME_DELAY_TAUC_MIN 7.0
#define TIME_DELAY_TAUC_MIN INSTANT_PERIOD_WINDOWED_WINDOW_WIDTH
// 20160503 AJL - TEST!!!
//#define TIME_DELAY_TAUC_MAX 25.0  // 20160503 AJL - TEST!!!
#define TIME_DELAY_TAUC_MAX TIME_DELAY_T50


// support for reading data from pick files
#define PICK_FORMAT_NLL 0


// BRB HP ground velocity data

typedef struct {
    double amp_at_pick; // amplitude on brb hp at pick
    int n_amp_max; // maximum index for amp data
    int n_amp; // current index for amp data
    double* vel; // ground velocity amplitude data
    double peak_amp_vel; // peak of amplitude in velocity
    double* disp; // ground velocity amplitude data
    double peak_amp_disp; // peak of amplitude in displacement
    int have_used_memory;
    int ioffset_pick;
}
GrdVel;

// mB

typedef struct {
    double amp_at_pick; // amplitude on mb brb at pick
    double int_sum; // amplitude on mb first integration at end of packet
    int polarity; // polarity of current int_sum displacement pulse
    int peak_index_count; // count of index length of current int_sum displacement pulse
    int n_amplitude_max; // maximum index for amplitdue accumulation
    int n_amplitude; // current index for amplitude accumulatoin
    double* amplitude; // amplitude
    double* peak_dur; // duration of displacement pulse at each amplitude point (back-filled from duration calcualted at end of pulse)
    double mag;
    int have_used_memory;
}
mBData;

// mwp

typedef struct {
    double amp_at_pick; // amplitude on mwp brb at pick
    double int_sum; // amplitude on mwp first integration at end of packet
    double int_int_sum; //  mwp second integration of amplitude at end of packet
    double abs_int_int_sum; //  mwp integral of abs integral at end of packet
    double int_int_sum_mwp_peak; // amplitude on mwp second integration at end of packet
    int polarity; // polarity of current int_int_sum_mwp_peak sum
    int peak_index_count; // count of index length of current int_int_sum_mwp_peak sum
    int n_peak_max; // maximum index for peak accumulation
    int n_peak; // current index for peak accumulation
    double* peak; // peak values of double integral signal
    double* peak_dur; // time delay after pick of peak values of double integral signal
    double* int_int; // mwp second integration of amplitude
    double* int_int_abs; //  mwp integral of abs integral
    //double pulse_first; // first pulse value of double integral signal
    //double pulse_first_dur; // time delay after pick of first pulse value of double integral signal
    double mag;
    int index_mag; // index used for peak for mag calculation
    int have_used_memory;
}
MwpData;

// t0

typedef struct {
    double amp_noise; // amplitude of noise before pick
    double amp_peak; // amplitude of peak
    int index_peak; //index since pick of peak
    double amp_90; // amplitude at percentage of peak
    int index_90; //index since pick of percentage of peak
    double amp_80; // amplitude at percentage of peak
    int index_80; //index since pick of percentage of peak
    double amp_50; // amplitude at percentage of peak
    int index_50; //index since pick of percentage of peak
    double amp_20; // amplitude at percentage of peak
    int index_20; //index since pick of percentage of peak
    double duration_raw; // HF P duration following Lomax & Michelini 2008, etc.
    double duration_corr; // HF duration using S duration if HF P duration extends past S arrival
    double duration_plot; // HF P duration following Lomax & Michelini 2008, etc.
    int have_used_memory;
}
T0Data;

// mwpd

typedef struct {
    double amp_at_pick; // amplitude on mwpd brb at pick
    double int_sum; // amplitude on mwpd first integration at end of packet
    int n_int_sum_max; // allocated size of int_int_sum_* arrays
    int n_int_sum; // highest intialized index of int_int_sum_* arrays
    double* int_int_sum_pos; // amplitude on mwpd second integration at end of packet for positive displacement
    double* int_int_sum_neg; // amplitude on mwpd second integration at end of packet for negative displacement
    double raw_mag;
    double corr_mag;
    int have_used_memory;
#ifdef USE_MWP_MO_POS_NEG
    // test of use of Mwpd positive and negative moment ratios as depth & tsunamigenic indicator
    // 20161226 AJL - added
    double mo_pos_neg_ratio;
#endif
}
MwpdData;

// waveform export

typedef struct {
    char filename[1024];
    hptime_t start_time_written;
    hptime_t end_time_written;
    long hypo_unique_id;

}
WaveformExport;

// polarization
// 20160810 AJL - added
#define POL_DONE 2   // set and done (maximum data window lengths were available and processed for all zyx)
#define POL_SET 1   // set
#define POL_NA 0    // not set
#define POL_UNK -1  // not available or not applicable

typedef struct {
    int status; // flag indicating polarization state: POL_DONE (set and done), POL_SET (set), POL_UNK (not set), POL_NA (not available or not applicable)
    int n_analysis_tries; // number of times polarization analysis (usually with increasing data windows) was done
    int nvalues; // number of values used for polarization statistics
    double azimuth; // mean azimuth (0 -> 360deg E of N)
    double azimuth_unc; // 1 std-dev azimuth uncertainty
    double azimuth_calc; // calculated/predicted ray take-off azimuth
    double weight; // polarization weight (polarization analysis azimuth)
    double dip; // mean dip (-90deg - down -> 90deg up, following Vidale, 1986)  // NOTE: Not currently set or used.
    double dip_unc; // 1 std-dev dip uncertainty  // NOTE: Not currently set or used.
    double dip_calc; // calculated/predicted ray dip  // NOTE: Not currently set or used.

}
Polarization;
#define POLARIZATION_MIN_DEG_LINEARITY_TO_USE 0.8
#define POLARIZATION_BASELINE_UNCERTAINTY 15.0  // baseline (minimum) uncertainty in degrees of a polarization reading (e.g. due to unmodelled 3D structure, ...)
//#define POLARIZATION_BASELINE_UNCERTAINTY 30.0  // baseline (minimum) uncertainty in degrees of a polarization reading (e.g. due to unmodelled 3D structure, ...)
#define POLARIZATION_MAX_NUM_ANALYSIS_TRIES 10  // maximum number of times to apply polarization analysis
// polarization application thresholds
#define POLARIZATION_MAX_DISTANCE_USE 90.0  // maximum GC arc in degrees to use a polarization reading
#define POLARIZATION_DISTANCE_WEIGHT_DIST_MIN 10.0   // distance at which weight decreases with 1/dist, at closer distance weight = 1
//#define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN_POLARIZATION 10.0  // 20160812 AJL
#define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN_POLARIZATION 6.0  // 20161006 AJL
// 20160912 AJL  #define POLARIZATION_MAX_DISTANCE_USE 30.0  // maximum GC arc in degrees to use a polarization reading
// 20160912 AJL  #define POLARIZATION_DISTANCE_WEIGHT_DIST_MIN 5.0   // distance at which weight decreases with 1/dist, at closer distance weight = 1
// 20160912 AJL  #define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN_POLARIZATION 20.0  // 20160812 AJL
// 20160911 AJL  #define POLARIZATION_MAX_DISTANCE_USE 90.0  // maximum GC arc in degrees to use a polarization reading
// 20160911 AJL  #define POLARIZATION_DISTANCE_WEIGHT_DIST_MIN 10.0   // distance at which weight decreases with 1/dist, at closer distance weight = 1
// 20160911 AJL  #define SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN_POLARIZATION 6.0  // 20160812 AJL

// data

//typedef struct timedomainProcessingData* TimedomainProcessingDataPtr;

typedef struct timedomainProcessingData {
    char* sladdr;
    int n_int_tseries; // index of internet timeseries query structure for this data
    int source_id;
    char station[16];
    char location[16];
    char channel[16];
    char network[16];
    double deltaTime; // dt
    double lat; // latitude -90/90 deg
    double lon; // latitude -180/180 deg
    double elev; // elevation (m)
    double azimuth; // Azimuth of the sensor in degrees from north, clockwise
    double dip; // Dip of the instrument in degrees, down from horizontal
    double station_quality_weight;
    PickData* pickData; // initial pick
    int pick_stream; // pick stream, e.g. STREAM_HF or STREAM_RAW
    int use_for_location; // flags data to be used for location - allows picking of hf and raw streams without using picks from both streams for location
    int can_use_as_location_P; // flags data can be used as definitive location P for statistics, magnitudes, duration, etc.
    struct timedomainProcessingData* use_for_location_twin_data; // other data with related pick
    int merged; // flags that data was merged into another location
    int virtual_pick_index; // elapsed index count since initial pick time
    double sn_pick; // HF noise rms amplitude before pick
    double a_ref; // a_ref rms amplitude
    double t50; // t50 rms amplitude
    int flag_complete_t50; // = 1 if a_ref and t50 values have been set
    int flag_a_ref_not_ok; // = 1 if a_ref less than last of last completed timedomain-processing measure
    int flag_snr_hf_too_low; // = 1 if HF signal-to-noise-ratio too low or noise level <= 0.0
    int flag_clipped; // = 1 if data is or may be clipped
    int flag_non_contiguous; // = 1 if data is non contiguous (gap or overlap)     // 20111230 AJL - added
    double amp_max_in_aref_window; // maximum T50 smoothed amplitude in A_REF window
    double amp_min_in_aref_window; // minimum T50 smoothed amplitude in A_REF window following the maximum T50 smoothed amplitude
    // brb
    double sn_brb_pick; // BRB noise rms amplitude before pick
    double sn_brb_signal; // BRB noise rms amplitude after pick
    double flag_snr_brb_too_low; // = 1 if BRB signal-to-noise-ratio too low or noise level <= 0.0
    int flag_complete_snr_brb; // = 1 if BRB signal-to-noise-ratio value has been set
    double sn_brb_int_pick; // BRB noise rms amplitude before pick
    double sn_brb_int_signal; // BRB noise rms amplitude after pick
    double flag_snr_brb_int_too_low; // = 1 if BRB signal-to-noise-ratio too low or noise level <= 0.0
    int flag_complete_snr_brb_int; // = 1 if BRB int signal-to-noise-ratio value has been set
    double sn_brb_bp_pick; // BRB bp noise rms amplitude before pick
    double sn_brb_bp_signal; // BRB bp noise rms amplitude after pick
    double flag_snr_brb_bp_too_low; // = 1 if BRB bp signal-to-noise-ratio too low or noise level <= 0.0
    int flag_complete_snr_brb_bp; // = 1 if BRB bp signal-to-noise-ratio value has been set
    // BRB HP ground velocity data
    GrdVel* grd_mot;
    int flag_complete_grd_mot; // = 1 if all grd motion amp values have been set
    // tauc
    double tauc_peak; // tauc peak amplitude
    // Mwp
    MwpData* mwp;
    int flag_complete_mwp; // = 1 if Mwp values have been set
    // T0
    T0Data* t0;
    int flag_complete_t0; // = 1 if T0 values have been set
    // Mwpd
    MwpdData* mwpd;
    int flag_complete_mwpd; // = 1 if Mwpd values have been set
    // mB
    mBData* mb;
    int flag_complete_mb; // = 1 if mB values have been set
    // time of pickdata pick
    int year;
    int month;
    int mday;
    int hour;
    int min;
    double dsec;
    char pick_error_type[16];
    double pick_error;
    time_t t_time_t; // pick epoch time integer seconds
    double t_decsec; // pick epoch time decimal seconds remainder
    int brb_polarity; // broad-band polarity measure
    double brb_polarity_quality; // broad-band polarity measure
    double brb_polarity_delay; // peak duration used for polarity measure
    // event association
    //int num_assoc;
    int is_associated; // flag indicating if this data is associated to an event (0 = NO; n>0 = YES, n is event index starting from 1)
    int is_full_assoc_loc; // flag indicating how this data was associated to an event:
    //          0 = with association only w/o re-location (location.c->octtreeGlobalAssociationLocation_full())
    //          1 = with full association location (location.c->octtreeGlobalAssociationLocation_reassociateOnly())
    //          -1 = not associated
    //int is_associated_grid_level;
    double epicentral_distance; // degrees
    double epicentral_azimuth; // degrees CW from N
    double residual;
    double dist_weight;
    double loc_weight;
    char phase[16];
    int phase_id;
    double take_off_angle_inc; // degrees (0/down->180/up)
    double take_off_angle_az; // degrees CW from N
    // derived values
    double ttime_P;
    double ttime_S;
    double ttime_SminusP; // 20151117 AJL - added so can check in timedomain_processing.c->td_process_timedomain_processing() that s/n window does not pass S time
    WaveformExport* waveform_export;
    // amplitude attenuation
    double amplitude_error_ratio; // ratio between a_ref rms amplitude and that predicted for event (see location.c -> calcLinearRegressionPowerRelation)
    // station corrections  // 20150716 AJL - added to support output of used sta correction in pick csv file
    double sta_corr;
    // polarization     // 20160805 AJL - added
    Polarization polarization;
}
TimedomainProcessingData;


/* functions */

TimedomainProcessingData* init_TimedomainProcessingData(double deltaTime, int flag_do_mwpd, int waveform_export_enable);
TimedomainProcessingData* new_TimedomainProcessingData(char* sladdr, int n_int_tseries, int source_id, char* station, char* location, char* channel,
        char* network, double deltaTime, double lat, double lon, double elev, double azimuth, double dip, double station_quality_weight, PickData* pickData, int pick_stream, int init_ellapsed_index_count,
        int year, int month, int mday, int hour, int min, double dsec, int flag_do_mwpd, int waveform_export_enable);
TimedomainProcessingData* copy_TimedomainProcessingData(TimedomainProcessingData* deData);
void free_TimedomainProcessingData(TimedomainProcessingData* deData);
int addTimedomainProcessingDataToDataList(TimedomainProcessingData* deData, TimedomainProcessingData*** pdata_list, int* pnum_data,
        int check_exists, int send_message);
void removeTimedomainProcessingDataFromDataList(TimedomainProcessingData* deData, TimedomainProcessingData*** pdata_list, int* pnum_data);
void free_TimedomainProcessingDataList(TimedomainProcessingData** data_list, int num_data);
int fprintf_TimedomainProcessingData(TimedomainProcessingData* deData, FILE* pfile);
int fprintf_NLLoc_TimedomainProcessingData(TimedomainProcessingData* deData, FILE* pfile, int append_evt_ndx_to_phase, double hypo_depth);
char setPolarity(TimedomainProcessingData* deData, double *pfmquality, int *pfmpolarity, char *fmtype);
time_t get_time_t(TimedomainProcessingData* deData, double* dsec);
double calculate_init_P_grd_mot_amp(TimedomainProcessingData* deData, double snr_brb, double snr_brb_int, int force_polarity, double *pfmquality, int *pfmpolarity, char *fmtype);
double calculate_Mwp_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mwp_distance_correction, int apply_mag_atten_check);
double calculate_mb_Mag(TimedomainProcessingData* deData, double hypo_depth, int apply_mag_atten_check);
double calculate_mB_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mb_correction, int apply_mag_atten_check);
double calculate_Raw_Mwpd_Mag(TimedomainProcessingData* deData, double hypo_depth, int use_mwp_distance_correction, int apply_mag_atten_check);
double calculate_corrected_Mwpd_Mag(double raw_mwpd_mag, double duration, double dominant_period, double depth);
double calculate_duration(TimedomainProcessingData* deData);
double calculate_duration_plot(double duration_raw);
void set_derived_values(TimedomainProcessingData* deData, double hypo_depth);
//int calculate_BRB_first_motion_polarity(TimedomainProcessingData* deData, double *pfmquality);


int is_unassociated_or_location_P(TimedomainProcessingData* deData);
int is_associated_location_P(TimedomainProcessingData* deData);
int is_associated_phase(TimedomainProcessingData* deData);

char* pick_stream_name(TimedomainProcessingData* deData);

// 20161006 AJL - following added to support running of Early-est with pick files
TimedomainProcessingData * get_next_pick(FILE *fp_in, int pick_format, int verbose);
TimedomainProcessingData *get_next_pick_nll(FILE *fp_in, int verbose);
