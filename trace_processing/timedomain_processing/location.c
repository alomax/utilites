// phase arrivals in the ak135 model

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#include "../ran1/ran1.h"
//#include "../geometry/geometry.h"
#include "../alomax_matrix/alomax_matrix.h"
#include "../alomax_matrix/alomax_matrix_svd.h"
//#include "../matrix_statistics/matrix_statistics.h"
#include "../statistics/statistics.h"
#include "../picker/PickData.h"
#include "timedomain_processing_data.h"
#include "ttimes.h"
#include "location.h"
//#include "../ran1/ran1.h"
#include "../response/response_lib.h"

#include "../statistics/student_t_table_1000_68.h"

// 20150728 AJL - performs classic oct-tree search over all levels simultaneously
#define PURE_OCTREE

// FOR DEBUG ONLY
#include "timedomain_processing.h"

#define INVALID_DISTANCE (-1.0e30)
#define INVALID_AZIMUTH (-1.0e30)
#define INVALID_WEIGHT (-1.0e30)

// use library exp function
#define EXP_FUNC exp
// use fast, approximate exp function
//#define EXP_FUNC EXP_FAST_APPROX  // DO NOT USE - makes no difference in total time, but error in exp() is 2-3%

/* fast, approximate exp function
 *  from: Nicol N. Schraudolph, (1999). A Fast, Compact Approximation of the Exponential Function,
 *  Neural Computation 11, 853–862
 */
/*
static union {
    double d;
    struct {
#ifdef LITTLE_ENDIAN
        int j, i;
#else
        int i, j;
#endif
    } n;
} _eco;
#define EXP_A (1048576/M_LN2) // use 1512775 for integer version
#define EXP_C 60801	// see text for choice of c values
#define EXP_FAST_APPROX(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)
 */
/*
 * End - fast, approximate exp function
 */

/** returns 1 if phase_id is a phase that counts for location origin time stacking/statistics
 */

// structure to hold hypocenter values for best location found
static GlobalBestValues globalBestValues;

/** determine if a pick phase should be counted in location weighting */

int count_in_location(int phase_id, double weight, int use_for_location) {

    // do not count if pick stream not same as specified locate stream
    if (!use_for_location) {
        return (0);
    }

    return (weight > 0.0 && (get_phase_type_flag(phase_id) & COUNT_IN_LOCATION));

}

/** initialize a HypocenterDesc object
 */

static Vect3D expect_NULL = {0.0, 0.0, 0.0};
static Mtrx3D cov_NULL = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static Ellipsoid3D ellipsoid_NULL = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static Ellipse2D ellipse_NULL = {0.0, 0.0, 0.0};

void init_HypocenterDesc(HypocenterDesc *phypo) {

    phypo->otime = 0.0;
    phypo->errh = 0.0;
    phypo->errz = 0.0;

    // flags
    phypo->hyp_assoc_index = -1;
    phypo->has_valid_magnitude = 0;

    // id
    phypo->unique_id = -1;
    phypo->first_assoc_time = LONG_MIN / 2; // force negative value when differenced with otime
    phypo->alert_sent_count = 0;
    phypo->alert_sent_time = 0;

    phypo->expect = expect_NULL;
    phypo->cov = cov_NULL;
    phypo->ellipsoid = ellipsoid_NULL;
    phypo->ellipse = ellipse_NULL;

    // amplitude attenuation
    phypo->linRegressPower.nvalues = -1.0;

    // event persistence fields
    phypo->n_possible_assoc_P = 0;

    // station health
    phypo->nstaHasBeenActive = -1;
    phypo->nstaIsActive = -1;

    // location sequence number
    phypo->loc_seq_num = -1;

    // location sequence number
    phypo->loc_type = LOC_TYPE_UNDEF;

}

/** allocate and initialize a HypocenterDesc object
 */

HypocenterDesc* new_HypocenterDesc() {

    HypocenterDesc* phypo = calloc(1, sizeof (HypocenterDesc));

    init_HypocenterDesc(phypo);

    return (phypo);

}

/** add a HypocenterDesc to a HypocenterDesc list
 * will remote any existing HypocenterDesc that is identified as same event
 * list will be sorted by increasing hypo->otime
 */

#define SIZE_INCREMENT_HYPO 16

int addHypocenterDescToHypoList(HypocenterDesc* pnew_hypo_desc, HypocenterDesc*** phypo_list, int* pnum_hypocenters, int icheck_duplicates, HypocenterDesc* pexisting_hypo_desc,
        HypocenterDesc** phypocenterDescInserted) {

    int new_hypocenter = 0;

    if (*phypo_list == NULL) { // list not yet created
        *phypo_list = calloc(SIZE_INCREMENT_HYPO, sizeof (HypocenterDesc*));
        // 20130930 AJL - bug fix
        //} else if (*pnum_hypocenters != 0 && (*pnum_hypocenters % SIZE_INCREMENT_HYPO) == 0) { // list will be too small
    } else if (*pnum_hypocenters + 1 > (int) (sizeof (*phypo_list) / sizeof (HypocenterDesc*))) { // list will be too small
        HypocenterDesc** newHypoList = NULL;
        newHypoList = calloc(*pnum_hypocenters + SIZE_INCREMENT_HYPO, sizeof (HypocenterDesc*));
        int n;
        for (n = 0; n < *pnum_hypocenters; n++)
            newHypoList[n] = (*phypo_list)[n];
        free(*phypo_list);
        *phypo_list = newHypoList;
    }

    // initialize unique_id unset
    long unique_id = pnew_hypo_desc->unique_id;
    long first_assoc_time = pnew_hypo_desc->first_assoc_time;
    long loc_seq_num = pnew_hypo_desc->loc_seq_num;
    int has_valid_magnitude = pnew_hypo_desc->has_valid_magnitude;
    int alert_sent_count = pnew_hypo_desc->alert_sent_count;
    int alert_sent_time = pnew_hypo_desc->alert_sent_time;

    // check if event is already present in list, if so remove existing HypocenterDesc
    if (icheck_duplicates) {
        int ncheck;
        for (ncheck = 0; ncheck < *pnum_hypocenters; ncheck++) {
            HypocenterDesc* phypo = (*phypo_list)[ncheck];
            if (isSameEvent(phypo, pnew_hypo_desc)) {
                // save copy of existing hypo
                if (pexisting_hypo_desc != NULL) {
                    *pexisting_hypo_desc = *phypo; // assumes HypocenterDesc has no fields using allocated memory
                    //printf("isSameEvent: phypo->messageTriggerThreshold %g %g %g %g\n", phypo->messageTriggerThreshold.mb, phypo->messageTriggerThreshold.mwp, phypo->messageTriggerThreshold.mwpd, phypo->messageTriggerThreshold.alarm);
                    //printf("isSameEvent: pexisting_hypo_desc->messageTriggerThreshold %g %g %g %g\n",
                    //        pexisting_hypo_desc->messageTriggerThreshold.mb, pexisting_hypo_desc->messageTriggerThreshold.mwp, pexisting_hypo_desc->messageTriggerThreshold.mwpd, pexisting_hypo_desc->messageTriggerThreshold.alarm);
                }
                unique_id = phypo->unique_id;
                // preserve persistent fields
                first_assoc_time = phypo->first_assoc_time;
                loc_seq_num = phypo->loc_seq_num;
                has_valid_magnitude = phypo->has_valid_magnitude;
                alert_sent_count = phypo->alert_sent_count;
                alert_sent_time = phypo->alert_sent_time;
                removeHypocenterDescFromHypoList(phypo, phypo_list, pnum_hypocenters);
                free(phypo);
            }
        }
    }

    if (unique_id <= 0) { // new hypocenter
        unique_id = (long) (1000.0 * pnew_hypo_desc->otime); // 1/1000 sec precision
        time_t current_time = time(&current_time);
        first_assoc_time = (long) current_time;
        new_hypocenter = 1;
    }
    pnew_hypo_desc->unique_id = unique_id;
    pnew_hypo_desc->first_assoc_time = first_assoc_time;
    pnew_hypo_desc->loc_seq_num = loc_seq_num;
    pnew_hypo_desc->has_valid_magnitude = has_valid_magnitude;
    pnew_hypo_desc->alert_sent_count = alert_sent_count;
    pnew_hypo_desc->alert_sent_time = alert_sent_time;

    // find first hypo later than new hypo
    int ninsert = 0;
    for (ninsert = 0; ninsert < *pnum_hypocenters; ninsert++) {
        if ((*phypo_list)[ninsert]->otime > pnew_hypo_desc->otime)
            break;
    }
    // shift later hypos
    if (ninsert < *pnum_hypocenters) {
        int m;
        for (m = *pnum_hypocenters - 1; m >= ninsert; m--)
            (*phypo_list)[m + 1] = (*phypo_list)[m];
    }
    // insert copy of new HypocenterDesc
    *phypocenterDescInserted = new_HypocenterDesc();
    //printf("DEBUG: ADD: %d sizeof(HypocenterDesc) %ld\n", *pnum_hypocenters, sizeof(HypocenterDesc));
    **phypocenterDescInserted = *pnew_hypo_desc; // simple copy is sufficient as long as HypocenterDesc contains no pointers to allocated objects
    (*phypo_list)[ninsert] = *phypocenterDescInserted;
    (*pnum_hypocenters)++;

    return (new_hypocenter);

}

/** remove a HypocenterDesc from a HypocenterDesc list */

void removeHypocenterDescFromHypoList(HypocenterDesc* phypo, HypocenterDesc*** phypo_list, int* pnum_hypocenters) {

    int i = 0;

    while (i < *pnum_hypocenters) {
        if ((*phypo_list)[i] == phypo)
            break;
        i++;
    }

    if (i == *pnum_hypocenters) // not found
        return;

    while (i < *pnum_hypocenters - 1) {
        (*phypo_list)[i] = (*phypo_list)[i + 1];
        i++;
    }

    (*phypo_list)[*pnum_hypocenters - 1] = NULL;

    //printf("DEBUG: REMOVE: %d sizeof(HypocenterDesc) %ld\n", *pnum_hypocenters, sizeof(HypocenterDesc));
    (*pnum_hypocenters)--;

}

/** clean up list memory */

void free_HypocenterDescList(HypocenterDesc*** phypo_list, int* pnum_hypocenters) {

    //printf("DEBUG: free_HypocenterDescList: %d sizeof(HypocenterDesc) %ld\n", *pnum_hypocenters, sizeof(HypocenterDesc));
    //if (*phypo_list == NULL || *pnum_hypocenters < 1)  // 20111004 AJL
    if (*phypo_list == NULL)
        return;

    int n;
    for (n = 0; n < *pnum_hypocenters; n++)
        free(*(*phypo_list + n)); // simple free is sufficient as long as HypocenterDesc contains no pointers to allocated objects
    //free_HypocenterDesc(*(hypo_list + n));
    *pnum_hypocenters = 0;

    free(*phypo_list);
    *phypo_list = NULL;

}

/** calculate and set hypocenter quality */

// 20150904 AJL - added

void setHypocenterQuality(HypocenterDesc* phypo, double min_weight_sum_assoc, double critical_errh, double critical_errz) {

    phypo->qualityIndicators.wt_count_assoc_weight = 1.0 - exp(-(phypo->qualityIndicators.weight_sum_assoc_unique - 3) / (2.0 * min_weight_sum_assoc));
    //phypo->qualityIndicators.wt_count_assoc_weight = 1.0 - exp(-(phypo->weight_sum - 3) / (2.0 * min_weight_sum_assoc));
    //phypo->qualityIndicators.wt_count_assoc_weight = 1.0 - exp(-(double) (phypo->nassoc_P - 3) / 6.0);
    phypo->qualityIndicators.errh_weight = exp(-phypo->errh / critical_errh);
    phypo->qualityIndicators.errz_weight = exp(-phypo->errz / critical_errz);
    double quality_weight =
            phypo->qualityIndicators.wt_count_assoc_weight
            * phypo->qualityIndicators.ot_variance_weight
            * phypo->qualityIndicators.errh_weight
            * phypo->qualityIndicators.errz_weight
            * phypo->qualityIndicators.amp_att_weight
            * phypo->qualityIndicators.gap_weight
            * phypo->qualityIndicators.distanceClose_weight
            * phypo->qualityIndicators.distanceFar_weight;
    quality_weight = pow(quality_weight, 1.0 / 8.0); // geometrical mean

    char *quality_code = phypo->qualityIndicators.quality_code;
    if (quality_code != NULL && sizeof (quality_code) > 1) {
        if (quality_weight > 0.9) {
            strcpy(quality_code, "A");
        } else if (quality_weight > 0.8) {
            strcpy(quality_code, "B");
            // 20160408 AJL  } else if (quality_weight > 0.5) {
        } else if (quality_weight > 0.6) { // 20160408 AJL
            strcpy(quality_code, "C");
        } else {
            strcpy(quality_code, "D");
        }
    }

    phypo->qualityIndicators.quality_weight = quality_weight;

}

/** initialize / clear deData association related fields
 */
void clear_deData_assoc(TimedomainProcessingData* deDataTmp) {

    deDataTmp->is_associated = 0;
    deDataTmp->epicentral_distance = -1.0;
    deDataTmp->epicentral_azimuth = -1.0;
    deDataTmp->residual = -999.0;
    deDataTmp->dist_weight = -1.0;
    deDataTmp->polarization.weight = -1.0;
    deDataTmp->polarization.azimuth_calc = -1.0;
    //deDataTmp->station_quality_weight = -1.0;
    deDataTmp->loc_weight = -1.0;
    strcpy(deDataTmp->phase, "X");
    deDataTmp->phase_id = -1;
    deDataTmp->take_off_angle_inc = -1.0;
    deDataTmp->take_off_angle_az = -1.0;
    // amplitude attenuation
    deDataTmp->amplitude_error_ratio = -1.0;
    deDataTmp->sta_corr = 0.0;

}

/** simple great-circle distance calculation on a sphere
 * returns distance in degrees
 *
 * Adapted from:
 *    Geographic Distance and Azimuth Calculations
 *    by Andy McGovern
 *    http://www.codeguru.com/algorithms/GeoCalc.html
 */

double GCDistance_local(double latA, double lonA, double latB, double lonB) {

    double d;

    if (latA == latB && lonA == lonB)
        return (0.0);

    latA *= DE2RA;
    lonA *= DE2RA;
    latB *= DE2RA;
    lonB *= DE2RA;

    d = sin(latA) * sin(latB) + cos(latA) * cos(latB) * cos(lonA - lonB);
    // 20120928 AJL return (AVG_ERAD * acos(d) * KM2DEG);
    return (acos(d) * RA2DE);

}

/** simple great-circle distance and azimuth calculation on a sphere
 *
 * returns distance and azimuth in degrees for great-circle from latA/lonA to latB/lonB
 *
 */

double GCDistanceAzimuth_local(double latA, double lonA, double latB, double lonB, double *pazimuth) {

    lonA *= DE2RA;
    latA *= DE2RA;
    lonB *= DE2RA;
    latB *= DE2RA;

    // distance
    double dist = sin(latA) * sin(latB) + cos(latA) * cos(latB) * cos(lonA - lonB);
    dist = acos(dist);

    // 20141106 AJL - added following check, to prevent div by 0 of sin(dist))
    if (dist < FLT_MIN) {
        *pazimuth = 0.0;
        return (dist * RA2DE);
    }

    // azimuth
    double cosAzimuth =
            (cos(latA) * sin(latB)
            - sin(latA) * cos(latB)
            * cos((lonB - lonA)))
            / sin(dist);
    double sinAzimuth =
            cos(latB) * sin((lonB - lonA)) / sin(dist);
    double az = atan2(sinAzimuth, cosAzimuth) * RA2DE;

    if (isnan(az) && fabs(lonB - lonA) < 0.000001) {
        if (latA > latB) {
            az = 180.0;
        } else {
            az = 0.0;
        }
    }

    if (az < 0.0) {
        az += 360.0;
    }

    *pazimuth = az;
    return (dist * RA2DE);

}


/** add a ValueDesc to a ValueDesc list
 * list will be sorted by increasing ValueDesc->value
 */

// 20111004 AJL - Added

#define SIZE_INCREMENT_VALUES 64

int addValueDescToValueList(ValueDesc* pnew_value_desc, ValueDesc*** pvalue_list, int* pnum_values) {

    if (*pvalue_list == NULL) { // list not yet created
        *pvalue_list = calloc(SIZE_INCREMENT_VALUES, sizeof (ValueDesc*));
        // 20130930 AJL - bug fix
        //} else if (*pnum_values != 0 && (*pnum_values % SIZE_INCREMENT_VALUES) == 0) { // list will be too small
    } else if (*pnum_values + 1 > (int) (sizeof (*pvalue_list) / sizeof (ValueDesc*))) { // list will be too small
        ValueDesc** newValueList = NULL;
        newValueList = calloc(*pnum_values + SIZE_INCREMENT_VALUES, sizeof (ValueDesc*));
        int n;
        for (n = 0; n < *pnum_values; n++)
            newValueList[n] = (*pvalue_list)[n];
        free(*pvalue_list);
        *pvalue_list = newValueList;
    }

    // find first value greater than new value
    int ninsert;
    for (ninsert = 0; ninsert < *pnum_values; ninsert++) {
        if ((*pvalue_list)[ninsert]->value > pnew_value_desc->value)
            break;
    }
    // shift newer value
    if (ninsert < *pnum_values) {
        int m;
        for (m = *pnum_values - 1; m >= ninsert; m--)
            (*pvalue_list)[m + 1] = (*pvalue_list)[m];
    }
    // insert new ValueDesc
    (*pvalue_list)[ninsert] = pnew_value_desc;
    (*pnum_values)++;

    return (*pnum_values);

}

/** remove a ValueDesc from a ValueDesc list */

// 20111004 AJL - Added

void removeValueDescFromValueList(ValueDesc* pvalue_desc, ValueDesc*** pvalue_list, int* pnum_values) {

    int i = 0;

    while (i < *pnum_values) {
        if ((*pvalue_list)[i] == pvalue_desc)
            break;
        i++;
    }

    if (i == *pnum_values) { // not found
        printf("ERROR: removeValueDescFromValueList: ValueDesc array element not found: this should not happen!\n");
        return;
    }

    while (i < *pnum_values - 1) {
        (*pvalue_list)[i] = (*pvalue_list)[i + 1];
        i++;
    }

    (*pvalue_list)[*pnum_values - 1] = NULL;

    (*pnum_values)--;

}

/** clean up list memory */

// 20111004 AJL - Added

void free_ValueDescList(ValueDesc*** pvalue_list, int* pnum_values) {

    if (*pvalue_list == NULL)
        return;

    free(*pvalue_list);
    *pnum_values = 0;
    *pvalue_list = NULL;

}

/** set deviation (range) of allowed origin time difference for comparing two event hypocenters
 */
#define OTIME_DEVIATION_MIN 1.0 // sec
double setRefOtimeDeviation(HypocenterDesc* hypo1, HypocenterDesc* hypo2) {

    double ref_ot_dev = hypo1->ot_std_dev + hypo2->ot_std_dev;
    ref_ot_dev *= 4.0;
    if (ref_ot_dev < OTIME_DEVIATION_MIN)
        ref_ot_dev = OTIME_DEVIATION_MIN;

    return(ref_ot_dev);

}

/** compare two HypocenterDesc's to determine if they refer to same event */
//20110412 AJL
#define EPICENTER_DEVIATION_MIN 10.0 // km
int isSameEvent(HypocenterDesc* hypo1, HypocenterDesc* hypo2) {

    // check unique_id
    if (hypo1->unique_id > 0 && (hypo1->unique_id == hypo2->unique_id))
        return (1);

    // check otime and epicenter, both must match

    // otime deviation
    double ot_dev = fabs(hypo1->otime - hypo2->otime);
    double ref_ot_dev = setRefOtimeDeviation(hypo1, hypo2);
    printf("DEBUG: ot_dev %f  ref_ot_dev %f\n", ot_dev, ref_ot_dev);
    // lat/lon deviation
    double dist_dev = GCDistance_local(hypo1->lat, hypo1->lon, hypo2->lat, hypo2->lon) * DEG2KM;
    double ref_epi_dev = hypo1->errh + hypo2->errh;
    printf("DEBUG: dist_dev %f  ref_epi_dev %f\n", dist_dev, ref_epi_dev);
    ref_epi_dev *= 4.0;
    if (ref_epi_dev < EPICENTER_DEVIATION_MIN)
        ref_epi_dev = EPICENTER_DEVIATION_MIN;
    // check for match within sum of allowable deviation ratios
    printf("DEBUG: ot_dev / ref_ot_dev %f  dist_dev / ref_epi_dev %f\n", ot_dev / ref_ot_dev, dist_dev / ref_epi_dev);
    if (
            // 20130308 AJL - change logic from sum to AND
            //ot_dev / ref_ot_dev + dist_dev / ref_epi_dev < 2.0
            ot_dev < ref_ot_dev && dist_dev < ref_epi_dev
            )
        return (1);

    return (0);

}

/* 20110412 AJL - replaced with simplified and more accurate version above
#define LAT_DEVIATION 2.0 // deg
#define LON_DEVIATION 2.0 // deg

int isSameEvent(HypocenterDesc* hypo1, HypocenterDesc* hypo2) {

    // check unique_id
    if (hypo1->unique_id > 0 && hypo1->unique_id == hypo2->unique_id)
        return (1);

    // check otime and epicenter, both must match
    double ot_dev = hypo1->ot_std_dev + hypo2->ot_std_dev;
    ot_dev *= 2.0;
    // correct allowable longitude deviation for convergence of meridiens
    double lon_dev = LON_DEVIATION;
    double cos_lat = cos(DE2RA * (hypo1->lat + hypo2->lat)) / 2.0;
    if (cos_lat > FLT_MIN)
        lon_dev /= cos_lat;
    else
        lon_dev = 999.9; // location at pole
    // correct longitude difference for wrap at +/- 180 deg great meridien
    double lon_diff = fabs(hypo1->lon - hypo2->lon);
    if (lon_diff > 180.0)
        lon_diff -= 360.0;
    // check for match within allowable deviations
    if (
            fabs(hypo1->otime - hypo2->otime) < ot_dev
            && fabs(hypo1->lat - hypo2->lat) < LAT_DEVIATION
            && fabs(lon_diff) < lon_dev
            )
        return (1);

    return (0);

}
 */






static OtimeLimit** OtimeLimitList = NULL;
static int NumOtimeLimit = 0;

/** otime limit class */

void setOtimeLimit(OtimeLimit* otime_limit, int use_for_location, int data_id, int pick_stream, double az, double dist,
        double dist_weight, double polarization_weight, double polarization_azimuth_calc, double quality_weight, double total_weight,
        double time, double otime, int polarity, int phase_id, int index, double dist_range, double time_range_uncertainty, double depth, TimedomainProcessingData* deData,
        double sta_corr) {

    otime_limit->data_id = data_id;
    otime_limit->use_for_location = use_for_location;
    otime_limit->pick_stream = pick_stream;
    otime_limit->azimuth = az;
    otime_limit->dist = dist;
    otime_limit->dist_weight = dist_weight;
    otime_limit->polarization_weight = polarization_weight;
    otime_limit->polarization_azimuth_calc = polarization_azimuth_calc;
    otime_limit->quality_weight = quality_weight;
    otime_limit->total_weight = total_weight;
    otime_limit->time = time;
    otime_limit->otime = otime;
    otime_limit->polarity = polarity;
    otime_limit->phase_id = phase_id;
    otime_limit->assoc = 0;
    otime_limit->index = index;
    otime_limit->dist_range = dist_range;
    otime_limit->time_range_uncertainty = time_range_uncertainty;

    otime_limit->count_in_loc = count_in_location(phase_id, total_weight, use_for_location);

    // 20140128 AJL - bug fix, added direct P check
    // 20140128 AJL  otime_limit->direct_P_within_dist_wt_dist_max = is_direct_P(phase_id) && dist_weight >= DISTANCE_WEIGHT_AT_DIST_MAX * DISTANCE_WEIGHT_AT_DIST_MAX;
    otime_limit->direct_P_within_dist_wt_dist_max = is_direct_P(phase_id) && dist <= DISTANCE_WEIGHT_DIST_MAX;

    // 20140121 AJL - added so that all data fields are available for association processing
    otime_limit->deData = deData;

    otime_limit->log_dist = -DBL_MAX;
    otime_limit->log_a_ref = -DBL_MAX;
    otime_limit->have_amp_atten_values = 0;

    if (use_amplitude_attenuation) {
        if (otime_limit->direct_P_within_dist_wt_dist_max && chan_resp[deData->source_id].have_gain
                && chan_resp[deData->source_id].responseType == DERIVATIVE_TYPE && deData->a_ref > 0.0) {
            if (dist > 0.0) {
                double st_line_dist = sqrt(dist * dist + depth * depth * KM2DEG * KM2DEG); // straight line distance
                // make sur st_line_dist is >> oct-tree cell size and thus station is well outside cell,
                //    otherwise amplitude attenuation estimate can be very erroneous
                if ((st_line_dist * DEG2KM) / dist_range > 2.0) {
                    otime_limit->log_dist = log(st_line_dist); // straight line distance
                    otime_limit->log_a_ref = log(deData->a_ref);
                    otime_limit->have_amp_atten_values = 1;
                }
            }
        }
    }

    otime_limit->sta_corr = sta_corr;


}


/** add an origin time limit to an OtimeLimit list
 * list will be sorted by increasing otime->time
 */

#define SIZE_INCREMENT_OTIME_LIMIT_LIST 128

void addOtimeLimitToList(OtimeLimit* otimeLimit, OtimeLimit*** potime_limit_list, int* pnum_otime_limit) {


    OtimeLimit** newOtimeLimitList = NULL;

    if (*potime_limit_list == NULL) { // list not yet created
        *potime_limit_list = calloc(SIZE_INCREMENT_OTIME_LIMIT_LIST, sizeof (OtimeLimit*));
        *pnum_otime_limit = 0;
        // 20130930 AJL - bug fix
        //} else if (*pnum_otime_limit != 0 && (*pnum_otime_limit % SIZE_INCREMENT_OTIME_LIMIT_LIST) == 0) { // list will be too small
    } else if (*pnum_otime_limit + 1 > (int) (sizeof (*potime_limit_list) / sizeof (OtimeLimit*))) { // list will be too small
        newOtimeLimitList = calloc(*pnum_otime_limit + SIZE_INCREMENT_OTIME_LIMIT_LIST, sizeof (OtimeLimit*));
        int n;
        for (n = 0; n < *pnum_otime_limit; n++)
            newOtimeLimitList[n] = (*potime_limit_list)[n];
        free(*potime_limit_list);
        *potime_limit_list = newOtimeLimitList;
    }

    // find first limit later than new limit
    int ninsert;
    for (ninsert = 0; ninsert < *pnum_otime_limit; ninsert++)
        if ((*potime_limit_list)[ninsert]->time > otimeLimit->time)
            break;
    // shift later data
    if (ninsert < *pnum_otime_limit) {
        int m;
        for (m = *pnum_otime_limit - 1; m >= ninsert; m--)
            (*potime_limit_list)[m + 1] = (*potime_limit_list)[m];
    }
    // insert new OtimeLimit
    (*potime_limit_list)[ninsert] = otimeLimit;
    (*pnum_otime_limit)++;

}

/** remove an origin time limit to an OtimeLimit list
 */

int removeOtimeLimitFromList(OtimeLimit* otimeLimit, OtimeLimit*** potime_limit_list, int* pnum_otime_limit) {

    if (otimeLimit == NULL)
        return (0);

    if (*potime_limit_list == NULL) // list not yet created
        return (0);

    // find first limit later than limit to remove
    int nremove;
    for (nremove = 0; nremove < *pnum_otime_limit; nremove++)
        if ((*potime_limit_list)[nremove] == otimeLimit)
            break;

    // if found, shift later data
    if (nremove < *pnum_otime_limit) {
        int m;
        for (m = nremove; m < *pnum_otime_limit - 1; m++)
            (*potime_limit_list)[m] = (*potime_limit_list)[m + 1];
        (*pnum_otime_limit)--;
        return (0);
    } else {
        printf("ERROR: removeOtimeLimitFromList: OtimeLimit %d not found in list!\n", otimeLimit->index);
        return (-1);
    }

}

/** clean up data memory */

void free_OtimeLimit(OtimeLimit* otime_limit) {

    if (otime_limit == NULL)
        return;

    free(otime_limit);

}

/** clean up list memory */

void free_OtimeLimitList(OtimeLimit*** potime_limit_list, int* pnum_otime_limit) {

    if (*potime_limit_list == NULL)
        return;

    int n;
    for (n = 0; n < *pnum_otime_limit; n++)
        free_OtimeLimit(*(*potime_limit_list + n));

    free(*potime_limit_list);
    *potime_limit_list = NULL;

    *pnum_otime_limit = 0;

}

/** clear list memory */

void clear_OtimeLimitList(int* pnum_otime_limit) {

    *pnum_otime_limit = 0;

}

/** check for and resolve overlao between otime limits for count_in_loc OtimeLimits for this ndata */

void checkOverlapCountInLocation(OtimeLimit* otimeLimitMin, OtimeLimit* otimeLimitMax, OtimeLimit*** potime_limit_list, int* pnum_otime_limit) {

    if (!otimeLimitMin->count_in_loc)
        return;

    int data_id = otimeLimitMin->data_id;
    OtimeLimit* otime_limit_test;
    double otime_limit_test_min_time, otime_limit_test_max_time;

    // check all start limits with same data_id and count_in_loc for overlap
    int nlimit;
    for (nlimit = 0; nlimit < *pnum_otime_limit; nlimit++) {
        otime_limit_test = (*potime_limit_list)[nlimit];
        if (otime_limit_test->data_id == data_id && otime_limit_test->count_in_loc && otime_limit_test->polarity > 0 && otime_limit_test != otimeLimitMin) {
            otime_limit_test_min_time = otime_limit_test->time;
            otime_limit_test_max_time = otime_limit_test->pair->time;
            if (otimeLimitMax->time < otime_limit_test_min_time) { // target limits end before start of test limits
                // no overlap
                continue;
            }
            if (otimeLimitMin->time > otime_limit_test_max_time) { // target limits start after after end of test limits
                // no overlap, past time of target
                return;
            }
            // overlap
            // resolve overlap by truncating limits of target (lower priority?) phase
            if (otimeLimitMin->time < otime_limit_test_min_time && otimeLimitMax->time >= otime_limit_test_min_time) { // target limits overlap start of test limits
                removeOtimeLimitFromList(otimeLimitMax, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                otimeLimitMax->time = otime_limit_test_min_time - DBL_MIN; // truncate target to start of test limits
                addOtimeLimitToList(otimeLimitMax, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                continue;
            }
            if (otimeLimitMin->time <= otime_limit_test_max_time && otimeLimitMax->time > otime_limit_test_max_time) { // target limits overlap end of test limits
                removeOtimeLimitFromList(otimeLimitMin, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                otimeLimitMin->time = otime_limit_test_max_time + DBL_MIN; // truncate target to end of test limits
                addOtimeLimitToList(otimeLimitMin, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                continue;
            }
            // resolve overlap by making target limit not count_in_loc
            otimeLimitMin->count_in_loc = otimeLimitMax->count_in_loc = 0;
            return;
        }
    }

}

/** tests to check if data should be skipped for association
 */
int skipData(TimedomainProcessingData *deData, int num_pass) {

    // skip if associated in previous pass
    if (deData->is_associated && deData->is_associated < num_pass)
        return (1);

    // AJL 20100224 - changed to reflect BRB picking
    // skip if possible clip in data time span, if HF s/n ratio too low or if in previous coda
    //if (deData->flag_clipped || deData->flag_snr_hf_too_low || deData->flag_a_ref_not_ok)
    //    return(1);
    // skip if HF s/n ratio too low
    // 20131022 AJL - try using all picks for location, regardless of HF S/N
    if (!USE_SNR_HF_TOO_LOW_PICKS_FOR_LOCATION && deData->flag_snr_hf_too_low)
        return (1);
    // END 20131022 AJL - try using all picks for location, regardless of HF S/N
    // skip if in previous coda - NOTE: comment out this check to locate sub-events soon after initial rupture (e.g. 2009.09.29 17.48 Samoa Islands)
    // 20150407 AJL if (deData->flag_a_ref_not_ok)
    if (!USE_AREF_NOT_OK_PICKS_FOR_LOCATION && deData->flag_a_ref_not_ok)
        return (1);
    // skip if possible clip in data time span
    if (deData->flag_clipped)
        return (1);
    // 20121127 AJL - bug fix
    // skip if non contiguous in data time span
    if (deData->flag_non_contiguous)
        return (1);

    return (0);

}

/**
 * Returns the mean of the maximum distance between each point on the sphere and all other points
 */
double getMeanDistanceMax(TimedomainProcessingData** pdata_list, int num_de_data) {

    int n, m;
    double distance;
    TimedomainProcessingData *deData0, *deData;

    // find mean of maximum distances
    double distanceMaxMean = 0.0;
    int icount = 0;
    for (n = 0; n < num_de_data; n++) {
        deData0 = pdata_list[n];
        if (!deData0->use_for_location)
            continue;
        double distanceMax = -FLT_MAX;
        for (m = 0; m < num_de_data; m++) {
            if (m != n) {
                deData = pdata_list[m];
                if (!deData->use_for_location)
                    continue;
                distance = GCDistance_local(deData0->lat, deData0->lon, deData->lat, deData->lon);
                if (distance > distanceMax) {
                    distanceMax = distance;
                }
            }
        }
        distanceMaxMean += distanceMax;
        icount++;
    }

    if (icount > 0)
        distanceMaxMean /= (double) icount;

    return (distanceMaxMean);

}

/** function to determine if total weight of unique associated stations/phases is greater than min_weight_sum_assoc.
 *
 */

#define MAX_NUM_LIMIT_CHECK_SAVE 1024 // should be much larger than min number of associated phases to give min_weight_sum_assoc
static OtimeLimit* OtimeLimitCheckSaved[MAX_NUM_LIMIT_CHECK_SAVE];

int checkWtSumUniqueStationsPhases(double min_weight_sum_assoc, TimedomainProcessingData** pdata_list, double *pweight_sum_assoc_unique, int *pnassociated_P_work) {

    int num_saved = 0;
    double weight_sum_assoc = 0.0;
    OtimeLimit *otimeLimit_assoc, *otimeLimit_saved;
    int j, n;
    for (j = 0; j < NumOtimeLimit; j++) {
        otimeLimit_assoc = *(OtimeLimitList + j);
        if (otimeLimit_assoc->polarity > 0 // only use lower limit
                && otimeLimit_assoc->assoc // if currently associated
                && otimeLimit_assoc->count_in_loc
                ) {
            double time_dec_sec = (double) otimeLimit_assoc->deData->t_time_t + otimeLimit_assoc->deData->t_decsec;
            for (n = 0; n < num_saved; n++) {
                otimeLimit_saved = OtimeLimitCheckSaved[n];
                // check if new limit has same station and phase as saved weight
                if (
                        otimeLimit_assoc->phase_id == otimeLimit_saved->phase_id
                        && strcmp(pdata_list[otimeLimit_assoc->data_id]->station, pdata_list[otimeLimit_saved->data_id]->station) == 0
                        ) {
                    // check if new limit has higher weight than saved weight
                    //if (otimeLimit_assoc->total_weight > otimeLimit_saved->total_weight) {
                    // check if new limit is earlier than saved limit   // 20160816 AJL - changed check so that earliest assoc P used, not highest weight
                    if (time_dec_sec < ((double) otimeLimit_saved->deData->t_time_t + otimeLimit_saved->deData->t_decsec)) {
                        weight_sum_assoc += otimeLimit_assoc->total_weight;
                        weight_sum_assoc -= otimeLimit_saved->total_weight;
                        otimeLimit_saved->total_weight = 0.0; // 20160814 AJL - added to flag duplicate sta/phases not used in weight_sum
                        otimeLimit_saved->assoc = 0; // 20160814 AJL - Bug fix, added
                        otimeLimit_saved->count_in_loc = 0; // 20160814 AJL - Bug fix, added
                        (*pnassociated_P_work)--; // 20160814 AJL - Bug fix, added
                        //if (weight_sum_assoc >= min_weight_sum_assoc)
                        //    return (1);
                        OtimeLimitCheckSaved[n] = otimeLimit_assoc; // replace
                    } else { // has lower weight
                        weight_sum_assoc -= otimeLimit_assoc->total_weight;
                        otimeLimit_assoc->total_weight = 0.0; // 20160814 AJL - added to flag duplicate sta/phases not used in weight_sum
                        otimeLimit_assoc->assoc = 0; // 20160814 AJL - Bug fix, added
                        otimeLimit_assoc->count_in_loc = 0; // 20160814 AJL - Bug fix, added
                        (*pnassociated_P_work)--; // 20160814 AJL - Bug fix, added
                    }
                    break; // found
                }
            }
            if (n >= num_saved) {
                // not found, save
                weight_sum_assoc += otimeLimit_assoc->total_weight;
                OtimeLimitCheckSaved[num_saved] = otimeLimit_assoc; // replace
                num_saved++;
                if (num_saved > MAX_NUM_LIMIT_CHECK_SAVE) {
                    printf("ERROR: checkWtSumUniqueStationsPhases: num_saved > MAX_NUM_LIMIT_CHECK_SAVE: this should not happen!\n");
                    return (1);
                }
            }
        }
    }

    *pweight_sum_assoc_unique = weight_sum_assoc;
    if (weight_sum_assoc >= min_weight_sum_assoc)
        return (1);
    //printf("DEBUG: checkWtSumUniqueStationsPhases: weight_sum_assoc %g < min_weight_sum_assoc %g\n", weight_sum_assoc, min_weight_sum_assoc);
    return (0);

}

/** function to determine primary and secondary azimuth gap as seen from epicenter (deg.
 *
 *     primary gap - largest azimuth separation between stations
 *     secondary gap - largest azimuth separation between stations filled by a single station
 *
 */

double calcAzimuthGap(ValueDesc** azimuths, int naz, double *pgap_secondary) {

    if (naz < 2) {
        *pgap_secondary = 360.0;
        return (360.0);
    }

    double az_last, az_last2, az;
    double gap_primary, gap_primary_max = -1.0;
    double gap_secondary, gap_secondary_max = -1.0;

    // find largest gap and secondary gap
    az_last2 = azimuths[naz - 2]->value - 360.0;
    az_last = azimuths[naz - 1]->value - 360.0;
    int narr;
    for (narr = 0; narr < naz; narr++) {
        az = azimuths[narr]->value;
        gap_primary = az - az_last;
        if (gap_primary > gap_primary_max)
            gap_primary_max = gap_primary;
        gap_secondary = az - az_last2;
        if (gap_secondary > gap_secondary_max)
            gap_secondary_max = gap_secondary;
        az_last2 = az_last;
        az_last = az;
    }

    *pgap_secondary = gap_secondary_max;
    return (gap_primary_max);

}


static ValueDesc StationDistances[MAX_NUM_SOURCES];
static ValueDesc** sta_dist_list = NULL; // 20111004 AJL - Added
static int num_sta_dist = 0;
static double stationDistances_lat = -999.0;
static double stationDistances_lon = -999.0;

/** set station distances in static array if lat/lon changed from last set
 *
 * returns number of station distances set
 */

// 20111004 AJL - Modified to handle ValueDesc list

int setStationDistances(double lat, double lon, ChannelParameters* channelParameters) {

    if (lat == stationDistances_lat && lon == stationDistances_lon)
        return (num_sta_dist);

    free_ValueDescList(&sta_dist_list, &num_sta_dist);

    int nsta;
    for (nsta = 0; nsta < num_sources_total; nsta++) {
        if (channelParameters[nsta].have_coords && !channelParameters[nsta].inactive_duplicate) {
            ChannelParameters* coords = channelParameters + nsta;
            StationDistances[nsta].value = GCDistance_local(lat, lon, coords->lat, coords->lon);
            addValueDescToValueList(StationDistances + nsta, &sta_dist_list, &num_sta_dist);
        } else {
            StationDistances[nsta].value = -1.0;
        }
    }

    stationDistances_lat = lat;
    stationDistances_lon = lon;

    return (num_sta_dist);

}

/** count the number of has been active stations less than a specified distance from a specified point
 *
 * returns number of stations
 */

// 20111004 AJL - Modified to handle ValueDesc list

int countStationsAvailable(double lat, double lon, double distance_max, ChannelParameters* channelParameters) {

    setStationDistances(lat, lon, channelParameters);

    int icount = 0;

    if (num_sta_dist > 0) {
        int narr;
        icount = num_sta_dist;
        for (narr = num_sta_dist - 1; narr >= 0; narr--) {
            icount = narr + 1;
            if (sta_dist_list[narr]->value < distance_max)
                break;
            icount--;
        }
    }

    return (icount);

}


#define COUNT_FAR 0  // count number avialable closer than secondary maximum distance in distance array
#define COUNT_CLOSE 1   // count number avialable closer than minimum distance in distance array
#define COUNT_ARITHMETIC_MEAN 2  // count number avialable closer than mean distance in distance array

/** count the number of has been active stations less than a specified distance from a specified point using specified counting function
 *
 * returns number of stations
 */

// 20111004 AJL - Modified to handle ValueDesc list

int countAllStationsAvailable(int mode, double limit_distance, ValueDesc** distances, int ndistances, double lat, double lon,
        ChannelParameters* channelParameters, double *pdistance_used, int *pNumWithinLimit) {


    // find required distance (far, close or mean)
    // count number in distance_array relative to specified limit distance and number available within limit distance
    int nWithinLimit = 0;
    int narr;
    if (mode == COUNT_FAR) { // count number associated at >= limit distance
        nWithinLimit = 0;
        double distance_max_secondary = limit_distance;
        if (ndistances > 1) {
            for (narr = ndistances - 2; narr >= 0; narr--) { // ignore most distant station, may be outlier
                if (distances[narr]->value < limit_distance)
                    break;
                nWithinLimit++;
            }
            distance_max_secondary = distances[ndistances - 2]->value;
        }
        *pdistance_used = distance_max_secondary;
        *pNumWithinLimit = nWithinLimit;
        return (countStationsAvailable(lat, lon, distance_max_secondary, channelParameters)); // count number at <= distance of second to farthest station associated
    } else if (mode == COUNT_CLOSE) { // count number associated at <= limit distance
        nWithinLimit = 0;
        if (ndistances > 1) {
            if (distances[ndistances - 1]->value <= limit_distance) { // all stations within limit distance
                nWithinLimit = ndistances;
            } else {
                for (narr = 0; narr < ndistances; narr++) {
                    if (distances[narr]->value > limit_distance)
                        break;
                    nWithinLimit++;
                }
            }
        }
        double distance_min = 0.0;
        if (ndistances > 0)
            distance_min = distances[0]->value;
        *pdistance_used = distance_min;
        *pNumWithinLimit = nWithinLimit;
        return (countStationsAvailable(lat, lon, distance_min, channelParameters)); // count number at <= distance of closest station associated
    } else if (mode == COUNT_ARITHMETIC_MEAN) {
        double distance_arithmetic_mean = limit_distance;
        double distance_sum = 0.0;
        nWithinLimit = 0;
        if (ndistances > 1) {
            for (narr = ndistances - 2; narr >= 0; narr--) { // ignore most distant station, may be outlier
                if (distances[narr]->value < limit_distance)
                    break;
                distance_sum += distances[narr]->value;
                nWithinLimit++;
            }
        }
        if (nWithinLimit > 0)
            distance_arithmetic_mean = distance_sum / (double) nWithinLimit;
        for (narr = ndistances - 2; narr >= 0; narr--) {
            if (distances[narr]->value <= distance_arithmetic_mean)
                nWithinLimit++;
        }
        *pNumWithinLimit = nWithinLimit;
        *pdistance_used = distance_arithmetic_mean;
        return (countStationsAvailable(lat, lon, distance_arithmetic_mean, channelParameters));
    }

    return (-999);

}

/** count the number of associated stations within a specified distance of hypocenter
 *
 * returns number of stations
 */

// 20111004 AJL - Added

int countStationsWithinDistance(double limit_distance, ValueDesc** distances, int ndistances) {


    // count number in distance_array within limit distance
    int nWithinLimit = 0;
    int narr;

    nWithinLimit = 0;
    if (ndistances > 1) {
        if (distances[ndistances - 1]->value <= limit_distance) { // all stations within limit distance
            nWithinLimit = ndistances;
        } else {
            for (narr = 0; narr < ndistances; narr++) {
                if (distances[narr]->value > limit_distance)
                    break;
                nWithinLimit++;
            }
        }
    }

    return (nWithinLimit);

}


/** find the distance from hypocenter within which there is the specified count of associated stations
 *
 * returns distance or -999 if few associated statoins than station_count
 */

// 20111004 AJL - Added

double findDistanceForStationCount(int station_count, ValueDesc** distances, int ndistances) {


    if (ndistances >= station_count) {
        return (distances[station_count - 1]->value);
    }

    return (-1.0);

}

/** Weighted counts outliers for a simple linear regression of log_a_ref vs. log_dist for specified value_desc
 *
 * returns number of values less[greater] than the specified min[max] ratios between observed and predicted values
 */

double rateOutliersLinearRegressionPowerRelation(ValueDesc** value_desc, int nvalue_desc, double min_dist_use, LinRegressPower *plinRegPow, double min_ratio, double max_ratio,
        int *pnvalues, int *pnoutliers) {

    double weight_sum = 0.0;
    double outlier_wt_sum = 0.0;
    double weight;

    int nvalues = 0;
    int noutliers = 0;

    for (int narr = 0; narr < nvalue_desc; narr++) {
        // check for valid a_ref value
        if (value_desc[narr]->log_dist > -DBL_MAX && value_desc[narr]->value >= min_dist_use) {
            double dist = value_desc[narr]->value;
            if (dist > FLT_MIN) {
                nvalues++;
                double a_ref_calc = plinRegPow->constant * pow(dist, plinRegPow->power);
                double amplitude_error_ratio = exp(value_desc[narr]->log_amp) / a_ref_calc;
                // weight by 1/dist so (fewer) closer readings count more
                //weight = 1.0 / dist;  // 20160114 AJL - may over-weight close stations which typically are amp outliers
                weight = 1.0;
                weight_sum += weight;
                if (amplitude_error_ratio < min_ratio || amplitude_error_ratio > max_ratio) {
                    outlier_wt_sum += weight;
                    noutliers++;
                }
            }
        }
    }

    double outlier_rate = -1.0;
    if (weight_sum > FLT_MIN) {
        outlier_rate = outlier_wt_sum / weight_sum;
    }

    *pnvalues = nvalues;
    *pnoutliers = noutliers;

    return (outlier_rate);

}

/** Performs a simple linear regression of log_a_ref vs. log_dist for specified value_desc array
 *
 * returns number of values with valid a_ref amplitudes found and used in the regression
 */

int calcLinearRegressionPowerRelation(ValueDesc** value_desc, int nvalue_desc, double min_dist_use, LinRegressPower *plinRegPow) {

    int nvalues = 0;
    plinRegPow->nvalues = nvalues;
    plinRegPow->power = -9.9;
    plinRegPow->constant = 0.0;
    plinRegPow->r_square = -1.0;
    plinRegPow->power_dev = 0.0;
    plinRegPow->constant_dev = 0.0;

    if (nvalue_desc < 3) {
        return (-1);
    }

    double x_sum = 0.0;
    double y_sum = 0.0;
    double xy_sum = 0.0;
    double xx_sum = 0.0;
    double yy_sum = 0.0;

    double logx, logy;

    // loop over distance values and accumulate regression sums
    for (int narr = 0; narr < nvalue_desc; narr++) {
        // check for valid a_ref value
        if (value_desc[narr]->log_dist > -DBL_MAX && value_desc[narr]->value >= min_dist_use) {
            logx = value_desc[narr]->log_dist;
            logy = value_desc[narr]->log_amp;
            x_sum += logx;
            y_sum += logy;
            xy_sum += logx * logy;
            xx_sum += logx * logx;
            yy_sum += logy * logy;
            nvalues++;
        }
    }

    // calculate regression values (http://en.wikipedia.org/wiki/Simple_linear_regression))
    double slope, intercept;
    if (nvalues > 3) {
        double nv = (double) nvalues;
        // estimate of slope and intercept
        slope = (xy_sum - x_sum * y_sum / nv) / (xx_sum - x_sum * x_sum / nv);
        plinRegPow->power = slope;
        intercept = (y_sum - slope * x_sum) / nv;
        //double intercept = (y_sum * xx_sum - x_sum * xy_sum) / (nv * xx_sum - x_sum * x_sum);
        plinRegPow->constant = exp(intercept);
        double correl_coeff = (xy_sum - x_sum * y_sum / nv) / sqrt((xx_sum - x_sum * x_sum / nv) * (yy_sum - y_sum * y_sum / nv));
        plinRegPow->r_square = correl_coeff * correl_coeff;
        // confidence interval
        double s_e_2 = (1.0 / (nv * (nv - 2))) * (nv * yy_sum - y_sum * y_sum - slope * slope * (nv * xx_sum - x_sum * x_sum));
        double s_slope = sqrt(nv * s_e_2 / (nv * xx_sum - x_sum * x_sum));
        int ndeg_freedom = nvalues - 2;
        if (ndeg_freedom > STUDENT_T_NUM_DEG_FREEDOM_MAX)
            ndeg_freedom = STUDENT_T_NUM_DEG_FREEDOM_MAX;
        double student_t_quantile = student_t_value[ndeg_freedom];
        //double student_t_quantile = 0.0;
        plinRegPow->power_dev = student_t_quantile * s_slope;
    } else {
        return (-1);
    }

    plinRegPow->nvalues = nvalues;

    return (nvalues);

}

/** returns variance and mean of an array of double values, checking for invalid values
 */

double calcVariance(double* darray, double* wtarray, int nvalues, double value_invalid, double *pmean, double *pdistance_min, double *pdistance_max) {

    double distance_min = FLT_MAX;
    double distance_max = -FLT_MAX;

    double sum = 0.0;
    double sum_sqr = 0.0;

    int nused = 0;
    double weight;
    double weight_sum = 0.0;

    int i;
    double dist;
    for (i = 0; i < nvalues; i++) {
        dist = darray[i];
        if (dist == value_invalid)
            continue;
        nused++;
        weight = wtarray[i];
        weight_sum += weight;
        sum += dist * weight;
        sum_sqr += dist * dist;
        if (dist < distance_min)
            distance_min = dist;
        if (dist > distance_max)
            distance_max = dist;
    }

    double variance = 0.0;

    if (nused < 1) {
        *pmean = 0.0;
        return (variance);
    } else if (nused == 1) {
        *pmean = sum / weight_sum;
        return (variance);
    }

    double mean = sum / weight_sum;
    // 20140226 AJL - bug fix
    //variance = (sum_sqr - mean * mean) / weight_sum;
    variance = (sum_sqr / weight_sum) - mean * mean;

    *pdistance_min = distance_min;
    *pdistance_max = distance_max;
    *pmean = mean;
    return (variance);

}



// global search association and location data

static int sizeOtimeLimitPool = 0;
static OtimeLimit* OtimeLimitPool = NULL;
static int sizeDataTmpPool = 0;
static TimedomainProcessingData *DataTmpPool = NULL;
// initial regular grid distance and azimuth arrays
#define USE_SAVED_INITIAL_GRIDS 1
static int sizeSavedInitialGridDistances = 0;
static double *SavedInitialGridDistances = NULL;
static int sizeSavedInitialGridAzimuths = 0;
static double *SavedInitialGridAzimuths = NULL;
static int sizeSavedInitialGridTravelTimeLatLon = 0;
static int sizeSavedInitialGridTravelTimeDepth = 0;
static double *SavedInitialGridTravelTimes = NULL;
// DEBUG
static int DEBUG_n_not_available = 0;
static int DEBUG_n_set = 0;
static int DEBUG_n_read = 0;
// distances array
static int sizeDistanceArray = 0;
static ValueDesc* DistanceArray = NULL;
static ValueDesc** distance_list = NULL; // 20111004 AJL - Added
static int num_distances = 0;
// azimuths array
static int sizeAzimuthArray = 0;
static ValueDesc* AzimuthArray = NULL;
static ValueDesc** azimuth_list = NULL; // 20111004 AJL - Added
static int num_azimuths = 0;
// weights array
//static int sizeWeightArray = 0;
//static double* WeightArray = NULL;

/** function to generate sample (scatter) of OctTree results */

int octtree_GenEventScatter(ResultTreeNode* resultTreeRoot, int value_type, int num_scatter, double oct_node_value_ref, float* fscatterdata,
        double integral, double *poct_tree_scatter_volume, int level_min, int level_max, char *sample_name) {

    int tot_npoints;
    int fdata_index;


    *poct_tree_scatter_volume = 0.0;

    /* return if no scatter samples requested */
    if (num_scatter < 1)
        return (0);

    printf("Info: Generating event scatter sample for %s, levels: %d-%d...", sample_name, level_min, level_max);


    /* generate scatter points at uniformly-randomly chosen locations in each leaf node */

    tot_npoints = 0;
    fdata_index = 0;
    tot_npoints = getScatterSampleResultTreeAtLevels(resultTreeRoot, value_type, num_scatter, integral,
            fscatterdata, tot_npoints, &fdata_index, oct_node_value_ref, poct_tree_scatter_volume, level_min, level_max);

    printf("   %d points generated, %d points requested, oct_tree_scatter_volume= %le\n",
            tot_npoints, num_scatter, *poct_tree_scatter_volume);

    return (tot_npoints);

}

/** function to calculate minimum and maximum travel-time in a grid cell */

int GetTTminMax(int phase_id, double dist_test_min, double dist_test_max, double depth, double grid_depth_test, double *ptt_min, double *ptt_max) {

    double time_test1, time_test2;

    // min
    time_test1 = get_ttime(phase_id, dist_test_min, depth + grid_depth_test);
    if (time_test1 < 0.0) // station to grid cell distance outside phase distance range
        return (-1);
    time_test2 = get_ttime(phase_id, dist_test_min, depth - grid_depth_test);
    if (time_test2 < 0.0) // station to grid cell distance outside phase distance range
        return (-1);
    *ptt_min = time_test1 < time_test2 ? time_test1 : time_test2;
    //tt_min -= regional_ttime_error;
    // max
    time_test1 = get_ttime(phase_id, dist_test_max, depth + grid_depth_test);
    if (time_test1 < 0.0) // station to grid cell distance outside phase distance range
        return (-1);
    time_test2 = get_ttime(phase_id, dist_test_max, depth - grid_depth_test);
    if (time_test2 < 0.0) // station to grid cell distance outside phase distance range
        return (-1);
    *ptt_max = time_test1 > time_test2 ? time_test1 : time_test2;

    return (0);

}

int isPhaseTypeToUse(TimedomainProcessingData* deData, int phase_id, int numPhaseTypesUse, int phaseTypesUse[], char channelNamesUse[][128],
        double timeDelayUse[MAX_NUM_TTIME_PHASES][2], TimedomainProcessingData** pdata_list, int ndata_test) {

    double time_dec_sec = (double) deData->t_time_t + deData->t_decsec;

    int n;
    for (n = 0; n < numPhaseTypesUse; n++) {
        if (phaseTypesUse[n] == phase_id) { // phase is type to use
            if (*channelNamesUse[n] == '\0' || strstr(channelNamesUse[n], deData->channel) != NULL) {
                // there is no channel to use specified, or data channel is specified to use
                //printf("DEBUG: isPhaseTypeToUse: %s  ndata_test %d  nphstype %d  delay:%f->%f\n", phase_name_for_id(phase_id), ndata_test, n, timeDelayUse[n][0], timeDelayUse[n][1]);
                // check time delay // 20150410 AJL - added
                if (timeDelayUse[n][0] >= 0.0 && timeDelayUse[n][1] >= 0.0) {
                    if (n > 0 && ndata_test > 0) { // need previous data to do this test
                        char *target_chans = channelNamesUse[n - 1];
                        double diff;
                        double time_dec_sec_prev;
                        int ndata;
                        for (ndata = ndata_test - 1; ndata >= 0; ndata--) { // data_list is in time order, work backwards from test data
                            TimedomainProcessingData* deData_prev = pdata_list[ndata];
                            //printf("  station: %s_%s %s_%s  phase_id: %d %d  chans: %s %s\n", deData_prev->network, deData->network, deData_prev->station, deData->station, deData_prev->phase_id, target_phase_id, target_chans, deData_prev->channel);
                            if (strcmp(deData_prev->network, deData->network) != 0 || strcmp(deData_prev->station, deData->station) != 0
                                    || strstr(target_chans, deData_prev->channel) == NULL
                                    )
                                continue;
                            // check if prev data is within required window before this data
                            time_dec_sec_prev = (double) deData_prev->t_time_t + deData_prev->t_decsec;
                            diff = time_dec_sec - time_dec_sec_prev;
                            //printf("  time diff: %f-%f=%f", time_dec_sec, time_dec_sec_prev, diff);
                            if (diff > timeDelayUse[n][1]) {
                                break; // all remaining differences will be too large
                            } else if (diff >= timeDelayUse[n][0]) {
                                //printf("  OK!\n");
                                return (1); // phase type, channel OK, and required phase/channel found within time delay before this data
                            }
                            //printf("  X\n");
                        }
                        //printf("  X\n");
                    }
                } else { // no time delay check required
                    return (1); // phase type and channel OK
                }
            }
        }
    }

    return (0);

}


ResultTreeNode* ppResultTreeRoot_LAST = NULL;

#define ICOUNT_INCREMENT 99999999
#define PROB_UNSET -FLT_MAX

int octtree_core(ResultTreeNode** ppResultTreeRoot, OctNode* poct_node, int indexLatLonSavedInitialGrid, int indexDepthSavedInitialGrid,
        TimedomainProcessingData** pdata_list, int num_de_data,
        int num_pass, double min_weight_sum_assoc, double critical_node_size_km,
        double gap_primary_critical, double gap_secondary_critical,
        int last_reassociate, double last_reassociate_otime, double last_reassociate_otime_sigma,
        int try_assoc_remaining_definitive, int no_reassociate, double distance_weight_dist_min_current,
        int numPhaseTypesUse, int phaseTypesUse[], char channelNamesUse[][128], double timeDelayUse[MAX_NUM_TTIME_PHASES][2],
        double reference_phase_ttime_error,
        ChannelParameters* channelParameters,
        GlobalBestValues *pglobalBestValues,
        /*double *pnode_prob, double *pnode_ot_variance,
        double *pglobal_best_ot_mean, double *pglobal_best_ot_variance, double *pglobal_best_lat, double *pglobal_best_lon, double *pglobal_best_depth,
        double *p_node_ot_variance_weight, int *pglobal_best_nassociated_P, int *pglobal_best_nassociated, double *pglobal_best_weight_sum,
        double *pglobal_best_prob, double *pglobal_effective_cell_size,
        int *pnCountClose, int *pNumCountDistanceNodeSizeTest, int *pnCountFar,
        double *pDistanceClose, double *pDistanceNodeSizeTest, double *pDistanceFar,
        // amplitude attenuation
        LinRegressPower *plinRegressPower,*/
        // oct-tree node volume
        double *pvolume,
        time_t time_min, time_t time_max
        ) {

    //printf("0 ============> channelParameters addr %ld\n", (long) channelParameters);

    double max_node_prob = -1.0;
    double max_node_ot_variance = -1.0;
    double max_node_ot_variance_weight = 0.0;

    // initialize local best solution variables for this node
    int is_global_best = 0;
    double best_prob = PROB_UNSET;
    double best_weight_sum = -1.0;
    double best_effective_cell_size = -1.0;


    TimedomainProcessingData* deData = NULL;

    // reference time is time of earliest arrival
    //time_t reftime = pdata_list[0]->t_time_t; // 1 sec precision
    double reftime = (double) pdata_list[0]->t_time_t + pdata_list[0]->t_decsec;

    double dist_test_min, dist_test_max;
    double tt_min, tt_max;
    double arrival_time;

    double gcd, gcaz;
    double otime_min, otime_max, otime;
    int n_ot_limit, j;


    // set parent value
    /*double parent_value = 0.0;
    if (poct_node->parent != NULL) {
        if (poct_node->parent->pdata == NULL) { // should not get here
            fprintf(stderr, "ERROR: parent OctTree node exists but has no pdata value!\n");
        } else {
            parent_value = *((double *) poct_node->parent->pdata);
        }
    }*/

    // initialize data
    /*if (poct_node->pdata == NULL)
        poct_node->pdata = (void *) malloc(sizeof (double));
    if (poct_node->pdata != NULL) {
     *((double *) poct_node->pdata) = 0.0;
    } else {
        fprintf(stderr, "ERROR: allocating double storage for OctTree node pdata.\n");
    }*/

    // check if node size <= critical
    //int node_size_le_critical = 0;
    //if (poct_node->ds.y <= critical_node_size_deg)
    //    node_size_le_critical = 1;

    // clear OtimeLimitList
    NumOtimeLimit = 0;

    // calculate approximate true cell volume
    double depth_corr = (AVG_ERAD - poct_node->center.z) / AVG_ERAD;
    double dsx_global = poct_node->ds.x * DEG2KM * cos(DE2RA * poct_node->center.y) * depth_corr;
    double dsy_global = poct_node->ds.y * DEG2KM * depth_corr;
    dsx_global *= depth_corr;
    dsy_global *= depth_corr;
    double volume = dsx_global * dsy_global * poct_node->ds.z;
    *pvolume = volume;

    //double node_diagonal = sqrt(3.0) * poct_node->ds.y * DEG2KM; // diagonal of cube with side node_size_km
    //double critical_node_diagonal = sqrt(3.0) * cos(DE2RA * poct_node->center.y)* critical_node_size_km; // diagonal of cube with side critical_node_size_km
    //double dsx_surface_km = poct_node->ds.x * DEG2KM * cos(DE2RA * poct_node->center.y);
    //double dsy_surface_km = poct_node->ds.y * DEG2KM;
    //double node_diagonal = sqrt(dsx_surface_km * dsx_surface_km + dsy_surface_km * dsy_surface_km); // diagonal
    double critical_node_diagonal = sqrt(3.0) * critical_node_size_km; // diagonal of cube with side critical_node_size_km
    //double critical_node_diagonal = node_diagonal * (critical_node_size_km  / dsy_surface_km); // diagonal of cube with side critical_node_size_km

    // set location
    double lon = poct_node->center.x;
    double lat = poct_node->center.y;
    double depth = poct_node->center.z;

    // set depth weight
    // 20160511 AJL - added to help avoid deep locations where not possible or likely
    double depth_for_weight = depth - poct_node->ds.z / 2.0;
    double depth_weight = get_dist_time_latlon_depth_weight(depth_for_weight, lat, lon);

    //printf("depth_for_weight %f  depth_weight %f\n", depth_for_weight, depth_weight);
    // TODO: replace here in association with true prior depth map

    // use large distance and depth margins: full grid cell size to each corner of central point
    double x_dist = poct_node->ds.x * cos(DE2RA * poct_node->center.y);
    double grid_dist_test = x_dist > poct_node->ds.y ? x_dist : poct_node->ds.y;
    grid_dist_test /= 2.0;
    grid_dist_test *= sqrt(2.0); // diagonal from center point to corner
    double grid_depth_test = poct_node->ds.z;
    grid_depth_test /= 2.0;
    //20110518 AJL - removed, gives too large depth in some cases, which gives travel time failure (get_ttime() return value -1)
    //grid_depth_test *= 1.1; // small extra margin

    // set surface P and S velocities
    double surface_property_P = get_velocity_model_surface_property('P');
    double surface_property_S = get_velocity_model_surface_property('S');

    // set travel times for this grid cell
    int ndata;
    for (ndata = 0; ndata < num_de_data; ndata++) {

        deData = pdata_list[ndata];
        /*if (last_reassociate && strcmp(deData->station, "ICESG") == 0) {
            printf("@@@@@@ %s,last_reassociate %d, deData->is_associated %d\n", deData->station, last_reassociate, deData->is_associated);
        }*/

        if (skipData(deData, num_pass))
            continue;

        // skip if not re-associating, and not associated in this pass
        if (no_reassociate && !(deData->is_associated && deData->is_associated == num_pass))
            continue;

        // get distance and azimuth
        if (USE_SAVED_INITIAL_GRIDS && indexLatLonSavedInitialGrid >= 0) {
            int ioffset = indexLatLonSavedInitialGrid * MAX_NUM_SOURCES + deData->source_id;
            if (*(SavedInitialGridDistances + ioffset) < 0.0) { // not initialized
                gcd = GCDistanceAzimuth_local(lat, lon, deData->lat, deData->lon, &gcaz);
                *(SavedInitialGridDistances + ioffset) = gcd;
                *(SavedInitialGridAzimuths + ioffset) = gcaz;
            } else {
                gcd = *(SavedInitialGridDistances + ioffset);
                gcaz = *(SavedInitialGridAzimuths + ioffset);
            }
        } else {
            gcd = GCDistanceAzimuth_local(lat, lon, deData->lat, deData->lon, &gcaz);
        }

        // set back azimuth invalid/unset (used for polarization analysis)
        double gc_back_az = DBL_INVALID;
        // set polarization distance weight invalid/unset (used for polarization analysis)
        double polarization_dist_wt = DBL_INVALID;

        // set min and max dist
        dist_test_min = gcd - grid_dist_test;
        dist_test_max = gcd + grid_dist_test;
        if (dist_test_min < 0.0) // station inside or near lat/lon cell
            dist_test_min = 0.0;
        // set distance weight
        double dist_weight;
#define TRUE_GEOMETRICAL_DECAY_DIST_WIEGHT 1 // AJL 20160405
#ifdef TRUE_GEOMETRICAL_DECAY_DIST_WIEGHT
        // AJL 20160405
        if (dist_test_min < DISTANCE_WEIGHT_DIST_MIN) { // use dist_test_min
            dist_weight = 1.0;
        } else { //if (dist_test_min < DISTANCE_WEIGHT_DIST_MAX) { // use C/dist^2
            dist_weight = (DISTANCE_WEIGHT_DIST_MIN * DISTANCE_WEIGHT_DIST_MIN) / (dist_test_min * dist_test_min);
            //} else { // use C/dist_max^2
            //    dist_weight = (DISTANCE_WEIGHT_DIST_MIN * DISTANCE_WEIGHT_DIST_MIN) / (DISTANCE_WEIGHT_DIST_MAX * DISTANCE_WEIGHT_DIST_MAX);
        }
#else
        // AJL 20100314
        //%%%        distance_weight_dist_min_current = 10.0; // 20150910 AJL - TEST: up-weight close-in stations
        if (dist_test_min < distance_weight_dist_min_current) { // use dist_test_min
            dist_weight = 1.0;
        } else if (dist_test_min < DISTANCE_WEIGHT_DIST_MAX) { // use dist_test_min
            // 1.0 to 0.2, gives full weight for < distance_weight_dist_min_current, zero weight for > DISTANCE_WEIGHT_DIST_MAX deg
            // 20101213 AJL TEST - maybe better to have higher weight for distant stations
            //dist_weight = 1.0 - 0.8 * (dist_test_min - distance_weight_dist_min_current) / (DISTANCE_WEIGHT_DIST_MAX - distance_weight_dist_min_current); // 20101124 AJL
            // AJL 20101215 reduce distance weighting strength
            dist_weight = 1.0 - (1.0 - DISTANCE_WEIGHT_AT_DIST_MAX) * (dist_test_min - distance_weight_dist_min_current) / (DISTANCE_WEIGHT_DIST_MAX - distance_weight_dist_min_current); // 20101124 AJL
            //printf("DEBUG: %d %s %s dist_weight %f, dist_test_min %f, distance_weight_dist_min_current %f\n",
            //        ndata, deData->network, deData->station, dist_weight, dist_test_min, distance_weight_dist_min_current); // 20101124 AJL
            //dist_weight = 1.0;   // 20101226 !!!
            dist_weight *= dist_weight; // 20130417 AJL - geometrical decay
        } else {
            //dist_weight = 0.2; // PKPdf,PKPab
            //dist_weight = 0.1; // PKPdf,PKPab 20110102
            dist_weight = DISTANCE_WEIGHT_AT_DIST_MAX * DISTANCE_WEIGHT_AT_DIST_MAX; // 20130417 AJL - geometrical decay
        }
#endif

        /*        //#define TRY_PERCENT_OF_TTIME_WT
        #ifdef TRY_PERCENT_OF_TTIME_WT
                // do not use distance weight
                dist_weight = 1.0;
        #endif
         */

        // set s/n ratio BRB for later use
        //double snr_brb = deData->a_ref < 0.0 || deData->sn_pick < FLT_MIN ? 0.0 : deData->a_ref / deData->sn_pick;
        double snr_brb = deData->sn_brb_signal < 0.0 || deData->sn_brb_pick < FLT_MIN ? 0.0 : deData->sn_brb_signal / deData->sn_brb_pick;

        // station quality weight
        double station_quality_weight = deData->station_quality_weight;
        double total_weight_station = dist_weight * station_quality_weight;
        // enable up-weighting of picks with high S/N for location   // 20141203 AJL - added
        if (upweight_picks_sn_cutoff > 0.0
                && dist_test_min <= upweight_picks_dist_max // 20151130 AJL - added
                ) {
            if (snr_brb > upweight_picks_sn_cutoff) {
                double upweight = 1.0 + (snr_brb - upweight_picks_sn_cutoff) / upweight_picks_sn_cutoff;
                if (dist_test_min > upweight_picks_dist_full) { // 20151130 AJL - added
                    upweight *= 1.0 - (dist_test_min - upweight_picks_dist_full) / (upweight_picks_dist_max - upweight_picks_dist_full);
                }
                total_weight_station *= upweight < 2.0 ? upweight : 2.0;
            }
        }
        // set invalid min and max gcd times
        tt_min = tt_max = -1.0;

        // set simple elevation corrections
        double elev_corr_P = (deData->elev / 1000.0) / surface_property_P;
        double elev_corr_S = (deData->elev / 1000.0) / surface_property_S;
        if (elev_corr_P < -10.0 || elev_corr_P > 10.0)
            fprintf(stderr, "Warning: unusual elevation correction: deData->elev=%g  elev_corr_P=%g   Vp=%g\n",
                deData->elev, elev_corr_P, surface_property_P);
        if (elev_corr_S < -20.0 || elev_corr_S > 20.0)
            fprintf(stderr, "Warning: unusual elevation correction: deData->elev=%g  elev_corr_S=%g  Vs=%g\n",
                deData->elev, elev_corr_S, surface_property_S);

        // loop over phases
        int phase_id;
        int have_direct_P = 0;
        double direct_P_otime = 0.0;
        OtimeLimit *otimeLimitMin_P = NULL, *otimeLimitMax_P = NULL;
        int have_direct_S = 0;
        double direct_S_otime = 0.0;
        OtimeLimit *otimeLimitMin_S = NULL, *otimeLimitMax_S = NULL;
        for (phase_id = 0; phase_id < get_num_ttime_phases(); phase_id++) {

            // skip if re-associating and before last reassociation pass and phase not counted in location
            if (!last_reassociate && !count_in_location(phase_id, 999.9, deData->use_for_location))
                continue;
            // skip if pick period too large
            /* 20110205 AJL - disabled this check - causes PKP phases to be skipped on large events (Mindanao_2010_MED)
            if (0 && !is_direct_P(phase_id) && deData->pickData->period > MAX_PICK_PERIOD_ASSOC_LOCATE) {
                //if (strcmp(phases[phase_id], "S") == 0)
                //    printf("=====================> DEBUG: skip pick %s_%s %d:%d:%f deData->pickData->period %f > MAX_PICK_PERIOD_ASSOC_LOCATE %f\n",
                //        deData->network, deData->station, deData->hour, deData->min, deData->dsec, deData->pickData->period, MAX_PICK_PERIOD_ASSOC_LOCATE);
                continue;
            }
             */

            // check if phase is in phase type to use list, if list exists
            if (numPhaseTypesUse > 0) {
                if (!isPhaseTypeToUse(deData, phase_id, numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse, pdata_list, ndata)) {
                    continue;
                }
            }

            // set min and max gcd times
            if (USE_SAVED_INITIAL_GRIDS && is_P(phase_id) && indexLatLonSavedInitialGrid >= 0 && indexDepthSavedInitialGrid >= 0) {
                //double *ttvalue = SavedInitialGridTravelTimes[indexLatLonSavedInitialGrid][indexDepthSavedInitialGrid][deData->source_id];
                double *ttvalue = SavedInitialGridTravelTimes + indexLatLonSavedInitialGrid * sizeSavedInitialGridTravelTimeDepth * MAX_NUM_SOURCES * 2
                        + indexDepthSavedInitialGrid * MAX_NUM_SOURCES * 2 + deData->source_id * 2;
                if (*ttvalue == DBL_INVALID) { // not available
                    DEBUG_n_not_available++;
                    continue;
                }
                if (*ttvalue < 0.0) { // not initialized
                    if (GetTTminMax(phase_id, dist_test_min, dist_test_max, depth, grid_depth_test, &tt_min, &tt_max) < 0) {
                        // not available
                        *ttvalue = DBL_INVALID;
                        DEBUG_n_not_available++;
                        continue;
                    }
                    *ttvalue = tt_min;
                    *(ttvalue + 1) = tt_max;
                    DEBUG_n_set++;
                } else {
                    tt_min = *ttvalue;
                    tt_max = *(ttvalue + 1);
                    DEBUG_n_read++;
                }
            } else {
                if (GetTTminMax(phase_id, dist_test_min, dist_test_max, depth, grid_depth_test, &tt_min, &tt_max) < 0) {
                    continue;
                }
            }

            // time_range_variance
            double time_range_uncertainty = 0.0;

            // get velocity at node center
            double velocity = is_P_at_source(phase_id) ? get_velocity_model_property('P', gcd, depth) : get_velocity_model_property('S', gcd, depth);
            if (velocity <= 0.0)
                fprintf(stderr, "ERROR: invalid velocity: gcd=%g depth=%g velocity=%g\n", gcd, depth, velocity);
            // 20101212 keeps minimum effective cell size large enough to associate most likely related phases for an event
            double time_range_min = critical_node_diagonal / velocity;
            double time_range_test = tt_max - tt_min;
            if (time_range_test < time_range_min) { // adjust time range to equal time_range_min
                tt_min = (tt_max + tt_min - time_range_min) / 2.0;
                tt_max = (tt_max + tt_min + time_range_min) / 2.0;
            }

            // set time range for ot variance test
            //double time_range_uncertainty = time_range_min / 2.0; // 20101219
            //double time_range_uncertainty = time_range_min; // 20101225
            //double time_range_uncertainty = current_min_node_size_km / velocity / 2.0; // 20101224
            //double time_range_uncertainty = node_diagonal / velocity / 2.0; // 20101226
            //time_range_uncertainty = (tt_max - tt_min) / 2.0; // 20101226
            time_range_uncertainty = (tt_max - tt_min) / 6.0; // 20150727 - seems to do better (than 2.0 or 4.0) for avoiding false events
            /*if (last_reassociate) {
                printf("DEBUG: time_range_uncertainty %g, grid_dist_test %g, grid_depth_test %g)\n", time_range_uncertainty, grid_dist_test, grid_depth_test);
            }*/

            double total_weight_phase = total_weight_station;

            // AJL 20100312 - added tt error
            double ttime_error = get_phase_ttime_error(phase_id);
            /*if (last_reassociate && strcmp(deData->station, "ICESG") == 0) {
                printf("@@@@@@ %s,last_reassociate %d, deData->is_associated %d\n", deData->station, last_reassociate, deData->is_associated);
            }*/
            //* TEST!
            // 20150910 AJL - added
            // if last re-associate, unassociated, definitive (have non-zero weight for location) picks are assigned increased ttime error (x3).
            // This increase association likelihood and avoid false events due to numerous unassociated picks for large events.
            // 20160114 AJL if (try_assoc_remaining_definitive && !deData->is_associated && count_in_location(phase_id, 999.9, 1)) {
            if (try_assoc_remaining_definitive && !deData->is_associated && count_in_location(phase_id, 999.9, deData->use_for_location)) {
                //printf("DEBUG %s_%s phase_id %d  deData->use_for_location %d  ttime_error/total_weight_phase %f/%f->", channelParameters[deData->source_id].network, channelParameters[deData->source_id].station, phase_id, deData->use_for_location, ttime_error, total_weight_phase);
                ttime_error *= 3.0;
                total_weight_phase = 0.0; // make sure this definitive phase cannot affect the association for this event
                //printf("%f/%f\n", ttime_error, total_weight_phase);
            } // */

            // AJL 20100312 - added tt error
            double time_errors = ttime_error + deData->pick_error;
            /*#ifdef TRY_PERCENT_OF_TTIME_WT
                        // 20150508 AJL - test using percentage of travel-time as ttime uncertainty
                        double tterror = 0.01 * (tt_max + tt_min) / 2.0; // uncertainty is % of mean ttime
                        tterror = fmax(tterror, time_errors); // min uncertainty is specified phase ttime error
                        //tterror = fmax(tterror, 0.5); // min uncertainty
                        //tterror = fmin(tterror, 10.0); // max uncertainty
                        time_errors = tterror;
            #endif*/
            double elev_corr = is_P_at_station(phase_id) ? elev_corr_P : elev_corr_S;
            tt_min = tt_min + elev_corr - time_errors;
            tt_max = tt_max + elev_corr + time_errors;
            time_range_uncertainty += time_errors; // 20101226

            // set minimal time range for calculating otime variance factor
            //double time_range_uncertainty = min_node_size / velocity / 2.0 + time_errors; // 20101222

            // calculate estimated origin time at middle of cell
            arrival_time = (double) deData->t_time_t + deData->t_decsec - reftime;
            // station correction
            // 20160601 AJL - changed from P only to any phase
            double sta_corr = 0.0;
            StaCorrections* psta_corr = &(channelParameters[deData->source_id].sta_corr[phase_id]);
            if (psta_corr->valid) { // based on test: psta_corr->num >= sta_corr_min_num_to_use
                double corr = 0.0;
                // 20150508 AJL added - polynomial only fit and used between dist max and min
                if (gcd >= psta_corr->dist_min && gcd <= psta_corr->dist_max) {
                    corr = psta_corr->poly[0];
                    double dist_pow = 1.0;
                    for (int i = 1; i < 4; i++) {
                        dist_pow *= gcd;
                        corr += dist_pow * psta_corr->poly[i];
                    }
                    /*printf("DEBUG: StaCorrections: %s_%s  %s  corr: %f  gcd: %f  coef: %g %g %g %g\n",
                            channelParameters[deData->source_id].network, channelParameters[deData->source_id].station,
                            phase_name_for_id(phase_id), corr, gcd, psta_corr->poly[0], psta_corr->poly[1], psta_corr->poly[2], psta_corr->poly[3]);*/
                }
                arrival_time -= corr;
                sta_corr = corr;
            }
            otime_max = arrival_time - tt_min;
            otime_min = arrival_time - tt_max;
            otime = (otime_max + otime_min) / 2.0;

            // check if estimated origin time is in analysis window
            // 20150505 AJL - added to speed up association
            // 20150507 AJL - removed: may remove phases from events late in analysis window, causing change in location, many orphaned phases and false events!
            //if ((long) (otime + reftime) < time_min) {
            //    //printf("DEBUG: otime < time_min, (int) otime %ld, time_min %ld \n", (long) (otime + reftime), time_min);
            //    continue;
            //}

            // check if estimated origin time range includes target otime
            // prevents associating to another event at same hypocenter location but with different otime
            // 20140722 AJL - added to support event persistence
            // 20140722 AJL - TODO: may allow small change in otime during association, is this OK???
            // 20160919 AJL if (last_reassociate && last_reassociate_otime > 0.0) {
            if (last_reassociate_otime > 0.0) {
                if (last_reassociate_otime - last_reassociate_otime_sigma - reftime > otime_max
                        || last_reassociate_otime + last_reassociate_otime_sigma - reftime < otime_min)
                    continue;
            }

            // check if earliest direct P and flag if direct arrival
            if (is_direct_P(phase_id) && have_direct_P) {
                // skip phase if already found related times for earlier direct phase
                if (otime > direct_P_otime) // later direct P phase (e.g. Pg)
                    continue;
                removeOtimeLimitFromList(otimeLimitMin_P, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                removeOtimeLimitFromList(otimeLimitMax_P, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
            }
            // check if earliest direct S and flag if direct arrival
            if (is_direct_S(phase_id) && have_direct_S) {
                // skip phase if already found related times for earlier direct phase
                if (otime > direct_S_otime) // later direct S phase (e.g. Sg)
                    continue;
                removeOtimeLimitFromList(otimeLimitMin_S, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
                removeOtimeLimitFromList(otimeLimitMax_S, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
            }

            // set distance range, effective cell linear dimension
            double dist_range = (otime_max - otime_min) * velocity;

            // include ratio of reference phase ttime error to this phase ttime errors for phase weighting
            // 20130307 AJL - added to allow use of non P phases for location
            //#ifdef TTIMES_REGIONAL
            //total_weight_phase = total_weight_phase * (reference_phase_ttime_error / get_phase_ttime_error(phase_id));
            total_weight_phase = total_weight_phase * (reference_phase_ttime_error / time_errors);
            //printf("DEBUG: total_weight_phase=%g total_weight=%g reference_phase_ttime_error=%g time_errors=%g\n", total_weight_phase, total_weight, reference_phase_ttime_error, time_errors);
            //#endif

            // adjust weight for polarization
            // 20160810 AJL - added
            double polarization_weight = -1.0;
            double polarization_azimuth_calc = -1.0;
            if (polarization_enable
                    && (deData->polarization.status == POL_DONE || deData->polarization.status == POL_SET)
                    && is_direct_P(phase_id)) {
                double az_pol = deData->polarization.azimuth;
                // get back azimuth station to epicenter
                // NOTE: assumes 1D model!  TODO: for 3D model must set and use deData->take_off_angle_az before this block
                if (gc_back_az == DBL_INVALID) {
                    GCDistanceAzimuth_local(deData->lat, deData->lon, lat, lon, &gc_back_az);
                }
                polarization_azimuth_calc = gc_back_az + 180.0; // polarization analysis uses azimuth at station away from epicenter
                if (polarization_azimuth_calc >= 360.0) {
                    polarization_azimuth_calc -= 360.0;
                }
                // unwrap GC arc azimuth relative to polarization azimuth
                if (az_pol - polarization_azimuth_calc > 180.0) {
                    polarization_azimuth_calc += 360.0;
                } else if (az_pol - polarization_azimuth_calc < -180.0) {
                    polarization_azimuth_calc -= 360.0;
                }
                // set polarization weight only if will be used for location
                if (total_weight_phase > FLT_MIN && snr_brb >= SIGNAL_TO_NOISE_RATIO_BRB_HP_MIN_POLARIZATION
                        && gcd <= POLARIZATION_MAX_DISTANCE_USE) {

                    // polarization weight based of prob of GC arc azimuth from cell to station +/-base_unc, given polarization azimuth Normal dist
                    double az_unc = deData->polarization.azimuth_unc;
                    // cumulative approach: weight is obs polarization Normal dist pdf: N(az_pol, az_unc) integrated in interval polarization_azimuth_calc-/+base_unc
                    polarization_weight = cumulDistNormal(polarization_azimuth_calc + POLARIZATION_BASELINE_UNCERTAINTY, az_pol, az_unc)
                            - cumulDistNormal(polarization_azimuth_calc - POLARIZATION_BASELINE_UNCERTAINTY, az_pol, az_unc);
                    // relative pdf approach (wt 0->1)
                    //polarization_weight = realtiveDensityNormal(az_pol - polarization_azimuth_calc, 0.0, az_unc + POLARIZATION_BASELINE_UNCERTAINTY);
                    // distance weight accounts for widening of azimuthal arc with distance
                    if (polarization_dist_wt == DBL_INVALID) {
                        polarization_dist_wt = 1.0;
                        if (gcd > POLARIZATION_DISTANCE_WEIGHT_DIST_MIN) {
                            polarization_dist_wt = POLARIZATION_DISTANCE_WEIGHT_DIST_MIN / gcd;
                        }
                    }
                    polarization_weight *= polarization_dist_wt;
                    total_weight_phase = total_weight_phase + polarization_weight * station_quality_weight;
                    // /* TEST! DO NOT USE! */ total_weight_phase = 2.0 * polarization_weight;
                    /*if (last_reassociate && try_assoc_remaining_definitive)
                        printf("DEBUG: %d %s_%s_%s_%s  %s  gcd: %f  az: %f  az_pol: %f  polarization_azimuth_calc: %f  az_unc: %lf,%lf  az_diff: %f  wt_fact: %f  wt: %f->%f\n",
                            ndata + 1, channelParameters[deData->source_id].network, channelParameters[deData->source_id].station,
                            channelParameters[deData->source_id].location, channelParameters[deData->source_id].channel,
                            phase_name_for_id(phase_id), gcd, gcaz, az_pol, polarization_azimuth_calc,
                            az_unc, POLARIZATION_BASELINE_UNCERTAINTY, az_pol - polarization_azimuth_calc,
                            polarization_weight, total_weight_phase - polarization_weight, total_weight_phase);*/
                }
            }

            // set OtimeLimits
            OtimeLimit *otimeLimitMin = OtimeLimitPool + NumOtimeLimit;
            setOtimeLimit(otimeLimitMin, deData->use_for_location, ndata, deData->pick_stream, gcaz, gcd, dist_weight, polarization_weight, polarization_azimuth_calc,
                    station_quality_weight, total_weight_phase, otime_min, otime, 1, phase_id, NumOtimeLimit, dist_range, time_range_uncertainty, depth, deData, sta_corr);
            addOtimeLimitToList(otimeLimitMin, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
            OtimeLimit *otimeLimitMax = OtimeLimitPool + NumOtimeLimit;
            setOtimeLimit(otimeLimitMax, deData->use_for_location, ndata, deData->pick_stream, gcaz, gcd, dist_weight, polarization_weight, polarization_azimuth_calc,
                    station_quality_weight, total_weight_phase, otime_max, otime, -1, phase_id, NumOtimeLimit, dist_range, time_range_uncertainty, depth, deData, sta_corr);
            addOtimeLimitToList(otimeLimitMax, &OtimeLimitList, &NumOtimeLimit); // increments NumOtimeLimit
            otimeLimitMin->pair = otimeLimitMax;
            otimeLimitMax->pair = otimeLimitMin;
            if (doCheckOverlapCountInLocation()) {
                // check for and resolve overlap between otime limits for count_in_loc OtimeLimits for this ndata
                //      prevents overlap of origin times corresponding to different phases for each pick,
                //      preventing multiple phase association for a pick
                checkOverlapCountInLocation(otimeLimitMin, otimeLimitMax, &OtimeLimitList, &NumOtimeLimit);
            }
            // check for direct P or S
            if (is_direct_P(phase_id)) {
                have_direct_P = 1;
                direct_P_otime = otime;
                otimeLimitMin_P = otimeLimitMin;
                otimeLimitMax_P = otimeLimitMax;
            } else if (is_direct_S(phase_id)) {
                have_direct_S = 1;
                direct_S_otime = otime;
                otimeLimitMin_S = otimeLimitMin;
                otimeLimitMax_S = otimeLimitMax;
            }
        }
    }

    // find largest peak of otime limit histogram
    int nassociated_P_work = 0;
    int nassociated_P_within_dist_wt_dist_max = 0;
    double dist_range_sum = 0.0;
    double time_range_sum_sqr = 0.0;
    double ot_sum = 0.0;
    double ot_sum_sqr = 0.0;
    double ot_mean;
    double ot_variance = 0.0;
    double ot_variance_weight = 0.0;
    double weight_sum = 0.0;
    OtimeLimit* otimeLimit;
    double total_wt;
    double gap_primary;
    double gap_secondary;
    int data_id;
    // prepare distance, azimuth weights arrays
    if (sizeDistanceArray < NumOtimeLimit) {
        sizeDistanceArray = NumOtimeLimit;
        DistanceArray = (ValueDesc*) realloc(DistanceArray, sizeDistanceArray * sizeof (ValueDesc));
        if (DistanceArray == NULL)
            fprintf(stderr, "ERROR: realloc of DistanceArray: %s\n", strerror(errno));
    }
    free_ValueDescList(&distance_list, &num_distances);
    if (sizeAzimuthArray < NumOtimeLimit) {
        sizeAzimuthArray = NumOtimeLimit;
        AzimuthArray = (ValueDesc*) realloc(AzimuthArray, sizeAzimuthArray * sizeof (ValueDesc));
        if (AzimuthArray == NULL)
            fprintf(stderr, "ERROR: realloc of AzimuthArray: %s\n", strerror(errno));
    }
    free_ValueDescList(&azimuth_list, &num_azimuths);
    /*
    if (sizeWeightArray < NumOtimeLimit) {
        sizeWeightArray = NumOtimeLimit;
        WeightArray = (double*) realloc(WeightArray, sizeWeightArray * sizeof (double));
        if (WeightArray == NULL)
            fprintf(stderr, "ERROR: realloc of WeightArray: %s\n", strerror(errno));
    }*/
    for (n_ot_limit = 0; n_ot_limit < NumOtimeLimit; n_ot_limit++) {
        DistanceArray[n_ot_limit].value = INVALID_DISTANCE;
        AzimuthArray[n_ot_limit].value = INVALID_AZIMUTH;
        // amplitude attenuation
        DistanceArray[n_ot_limit].deData = NULL;
        DistanceArray[n_ot_limit].log_dist = -DBL_MAX;
        DistanceArray[n_ot_limit].log_amp = -DBL_MAX;
        AzimuthArray[n_ot_limit].deData = NULL;
        //WeightArray[n_ot_limit] = INVALID_WEIGHT;
    }
    // parse otime limits in time order to scan stack (histogram) of weighted origin times
    int location_count_changed = 0;
    for (n_ot_limit = 0; n_ot_limit < NumOtimeLimit; n_ot_limit++) {
        otimeLimit = *(OtimeLimitList + n_ot_limit);
        data_id = otimeLimit->data_id;
        otime = otimeLimit->otime;
        total_wt = otimeLimit->total_weight;
        location_count_changed = 0;
        if (otimeLimit->polarity > 0) { // enter otime limit for this datum
            otimeLimit->assoc = otimeLimit->pair->assoc = 1;
            // count in location statistics only if is first and direct arrival (better than initial?) P phase, helps to prevent false locations with all or mostly later phases
            if (otimeLimit->count_in_loc) {
                location_count_changed = 1;
                ot_sum += total_wt * otime;
                ot_sum_sqr += total_wt * otime * otime;
                weight_sum += total_wt;
                dist_range_sum += total_wt * otimeLimit->dist_range;
                time_range_sum_sqr += total_wt * otimeLimit->time_range_uncertainty * otimeLimit->time_range_uncertainty;
                nassociated_P_work++;
                // 20140122 AJL - bug fix?
                //if (total_wt >= DISTANCE_WEIGHT_AT_DIST_MAX) { // only use direct P at distance < DISTANCE_WEIGHT_DIST_MAX for later distance and azimuth checking
                if (otimeLimit->direct_P_within_dist_wt_dist_max) { // only use direct P at distance < DISTANCE_WEIGHT_DIST_MAX for later distance, azimuth and attenuation checking
                    nassociated_P_within_dist_wt_dist_max++;
                    DistanceArray[otimeLimit->index].value = otimeLimit->dist;
                    if (use_amplitude_attenuation) {
                        if (otimeLimit->have_amp_atten_values) {
                            DistanceArray[otimeLimit->index].log_dist = otimeLimit->log_dist;
                            DistanceArray[otimeLimit->index].log_amp = otimeLimit->log_a_ref;
                        }
                        //    printf("DEBUG: is_direct_P(phase_id):%d, otimeLimit->dist_weight=%f, %f/%f -> %f/%f\n",
                        //            is_direct_P(otimeLimit->phase_id), otimeLimit->dist_weight,
                        //            otimeLimit->log_dist, otimeLimit->log_a_ref,
                        //            DistanceArray[otimeLimit->index].log_dist, DistanceArray[otimeLimit->index].log_a_ref);
                        DistanceArray[otimeLimit->index].deData = otimeLimit->deData; // all corresponding pick data is available with dist
                    }
                    addValueDescToValueList(DistanceArray + otimeLimit->index, &distance_list, &num_distances);
                    AzimuthArray[otimeLimit->index].value = otimeLimit->azimuth;
                    addValueDescToValueList(AzimuthArray + otimeLimit->index, &azimuth_list, &num_azimuths);
                    //WeightArray[otimeLimit->index] = total_wt;
                }
            }
            /*if (last_reassociate)
                printf("DEBUG: ASSOC+ %s weight_sum %g  nassociated_P_work %d  flag_a_ref_not_ok %d  data_id %d  otimeLimits: %s %s gcd %f dwt %f  ot %f - %f (%f) ot=%f\n",
                    phase_name_for_id(otimeLimit->phase_id),
                    weight_sum, nassociated_P_work, deData->flag_a_ref_not_ok, otimeLimit->data_id,
                    pdata_list[otimeLimit->data_id]->station, phase_name_for_id(otimeLimit->phase_id), otimeLimit->dist, otimeLimit->total_weight,
                    otimeLimit->time, otimeLimit->pair->time, otimeLimit->pair->time - otimeLimit->time, otimeLimit->otime);*/
        } else { // leave otime limit for this datum
            otimeLimit->assoc = otimeLimit->pair->assoc = 0;
            if (otimeLimit->count_in_loc) {
                location_count_changed = 1;
                ot_sum -= total_wt * otime;
                ot_sum_sqr -= total_wt * otime * otime;
                weight_sum -= total_wt;
                dist_range_sum -= total_wt * otimeLimit->dist_range;
                time_range_sum_sqr -= total_wt * otimeLimit->time_range_uncertainty * otimeLimit->time_range_uncertainty;
                nassociated_P_work--;
                // 20140122 AJL - bug fix?
                //if (total_wt >= DISTANCE_WEIGHT_AT_DIST_MAX) { // only use direct P at distance < DISTANCE_WEIGHT_DIST_MAX for later distance and azimuth checking
                if (otimeLimit->direct_P_within_dist_wt_dist_max) { // only use direct P at distance < DISTANCE_WEIGHT_DIST_MAX for later distance, azimuth and attenuation checking
                    nassociated_P_within_dist_wt_dist_max--;
                    DistanceArray[otimeLimit->pair->index].value = INVALID_DISTANCE;
                    // amplitude attenuation
                    DistanceArray[otimeLimit->pair->index].deData = NULL;
                    DistanceArray[otimeLimit->pair->index].log_dist = -DBL_MAX;
                    DistanceArray[otimeLimit->pair->index].log_amp = -DBL_MAX;
                    removeValueDescFromValueList(DistanceArray + otimeLimit->pair->index, &distance_list, &num_distances);
                    AzimuthArray[otimeLimit->pair->index].value = INVALID_AZIMUTH;
                    removeValueDescFromValueList(AzimuthArray + otimeLimit->pair->index, &azimuth_list, &num_azimuths);
                    //WeightArray[otimeLimit->pair->index] = INVALID_WEIGHT;
                }
            }
            /*if (last_reassociate)
                printf("DEBUG: ASSOC- %s weight_sum %g  nassociated_P_work %d  flag_a_ref_not_ok %d  data_id %d  otimeLimits: %s %s gcd %f dwt %f  ot %f - %f (%f) ot=%f\n",
                    phase_name_for_id(otimeLimit->phase_id),
                    weight_sum, nassociated_P_work, deData->flag_a_ref_not_ok, otimeLimit->data_id,
                    pdata_list[otimeLimit->data_id]->station, phase_name_for_id(otimeLimit->phase_id), otimeLimit->dist, otimeLimit->total_weight,
                    otimeLimit->time, otimeLimit->pair->time, otimeLimit->pair->time - otimeLimit->time, otimeLimit->otime);*/
            // check if not enough data remaining to get more stations than best
            // i.e. nassociated_P_work + max_new_possible < best_nassociated_P_work
            //if (nassociated_P_work + (NumOtimeLimit - n_ot_limit - 1) / 2 < best_nassociated_P_work)
            //    break;
        }

        //#define DEBUG_PROB 1
#define DEBUG_ALL 1

        // check quality of stack sum
        //if (nassociated_P_work > 1 && weight_sum > 1.0) {
        double adjusted_weight_sum = weight_sum - 1.0; // do not count first otime in stack for location probability
        //if (location_count_changed && nassociated_P_work > 1 && weight_sum >= min_weight_sum_assoc // 20101217
        //        && adjusted_weight_sum > best_prob
        //        ) {
#ifdef PURE_OCTREE
        if (!location_count_changed || nassociated_P_work <= 1 || weight_sum < min_weight_sum_assoc // 20101217
                //|| adjusted_weight_sum < best_prob
                )
            continue; // cannot be better than current best_prob
#else
        if (!location_count_changed || nassociated_P_work <= 1 || weight_sum < min_weight_sum_assoc // 20101217
                || adjusted_weight_sum < best_prob
                )
            continue; // cannot be better than current best_prob
#endif
#ifdef DEBUG_PROB
        if (DEBUG_ALL || last_reassociate) {
            printf("DEBUG_PROB: ======================================\n");
            printf("DEBUG_PROB: weight_sum [%f], min_weight_sum_assoc [%f]\n",
                    weight_sum, min_weight_sum_assoc);
        }
#endif
#ifdef DEBUG_PROB
        if (DEBUG_ALL || last_reassociate) {
            printf("DEBUG_PROB: adjusted_weight_sum [%f] = weight_sum - 1.0 [%f]\n",
                    adjusted_weight_sum, weight_sum - 1.0);
        }
#endif

        // calculate otime variance factor
        ot_mean = ot_sum / weight_sum;
        //ot_variance = (ot_sum_sqr - weight_sum * ot_mean * ot_mean) / weight_sum;
        //ot_variance = (ot_sum_sqr - weight_sum * ot_mean * ot_mean) / (weight_sum - 2.0); // 20101224 AJL - changed  - 1.0  to  - 2.0
        ot_variance = (ot_sum_sqr - weight_sum * ot_mean * ot_mean) / (weight_sum - min_weight_sum_assoc + 1.0); // 20101226 AJL - changed
        //double time_range_variance = time_range_sum_sqr / (weight_sum - 2.0); // 20101224 AJL - changed  - 1.0  to  - 2.0
        double time_range_variance = time_range_sum_sqr / (weight_sum - 1.0); // 20110128 AJL - changed  - 2.0  to  - 1.0
        //ot_variance_weight = time_range_variance / (ot_variance + time_range_variance);
        ot_variance_weight = EXP_FUNC(-ot_variance / time_range_variance);
        //ot_variance_weight = pow(EXP_FUNC(-ot_variance / time_range_variance), 2);   // 20101226
        /*
        if (last_reassociate) {
            printf("DEBUG: ot_variance_weight %g = exp(-ot_variance %g / time_range_variance %g)\n",
                    ot_variance_weight, ot_variance, time_range_variance);
        }
        //*/
        //if (ot_variance > time_range_variance)
        //    ot_variance_weight = time_range_variance / ot_variance;

        // calculate hypocenter "probability"
        double effective_cell_size = dist_range_sum / weight_sum;
#ifdef PURE_OCTREE
        // prob is function of ot variance, assoc weight and cell size; gives higher weight to small cells
        double prob = ot_variance_weight * adjusted_weight_sum / effective_cell_size;
#ifdef DEBUG_PROB
        if (DEBUG_ALL || last_reassociate) {
            printf("DEBUG_PROB: prob [%f] = ot_variance_weight [%f] * adjusted_weight_sum / effective_cell_size [%f]\n",
                    prob, ot_variance_weight, adjusted_weight_sum / effective_cell_size);
        }
#endif

        //VOLUMEdouble prob = ot_variance_weight * adjusted_weight_sum;
#else
        double prob = ot_variance_weight * adjusted_weight_sum;
#endif
        if (isnan(prob)) { // 20101230 AJL
            printf("DEBUG: ===> prob NaN!!! [ot_variance_weight * adjusted_weight_sum] volume=%lg   x y x = %g %g %g\n",
                    volume, poct_node->center.y, poct_node->center.x, poct_node->center.z);
            continue;
        }

        // include depth weight
        // 20160511 AJL - added to help avoid deep locations where not possible or likely
        prob *= depth_weight;
#ifdef DEBUG_PROB
        if (DEBUG_ALL || last_reassociate) {
            printf("DEBUG_PROB: prob [%f] *= depth_weight [%f]\n",
                    prob, depth_weight);
        }
#endif

        // check if best prob for this node
        //if (prob >= best_prob && prob > min_prob_assoc) {    // 20101227
        if (prob < best_prob) // 20101227
            continue; // cannot be better than current best_prob

        // perform a number of checks to further weight the hypocenter probability

        double distanceClose = 0.0;
        double distanceFar = 0.0;

        /*
        if (0) { // 20110513 AJL - check only distance range and azimuth gaps (much more efficient than using countAllStationsAvailable();

            double dist_mean = 0.0;
            double dist_variance = calcVariance(DistanceArray, WeightArray, NumOtimeLimit, INVALID_DISTANCE, &dist_mean, &distanceClose, &distanceFar);
            double dist_stddev = sqrt(dist_variance);
            double dist_low = dist_mean - dist_stddev;
            if (dist_low < distanceClose)
                dist_low = distanceClose;
            double dist_high = dist_mean + dist_stddev;
            if (dist_high > distanceFar)
                dist_high = distanceFar;
            double distance_ratio = dist_high / dist_low; // will be >= 1 or < 0
            double distance_exponent = 0.0;
            if (distance_ratio > 0.0 && distance_ratio < DISTANCE_RATIO_CRITICAL)
                distance_exponent = (DISTANCE_RATIO_CRITICAL - distance_ratio) / (DISTANCE_RATIO_CRITICAL - 1.0);

            gap_primary = calcAzimuthGapSparseArray(AzimuthArray, NumOtimeLimit, nassociated_P_within_dist_wt_dist_max, &gap_secondary);
            double gap_primary_exponent = 0.0;
            if (gap_primary > gap_primary_critical)
                gap_primary_exponent = (gap_primary - gap_primary_critical) / (360.0 - gap_primary_critical);
            double gap_secondary_exponent = 0.0;
            if (gap_secondary > gap_secondary_critical)
                gap_secondary_exponent = (gap_secondary - gap_secondary_critical) / (360.0 - gap_secondary_critical);
            double gap_exponent = (gap_primary_exponent + gap_secondary_exponent) / 2.0;

            double weight_product = EXP_FUNC(-0.5 * (distance_exponent + gap_exponent));

            if (last_reassociate) {
            //if (last_reassociate || distance_exponent > 0.0) {
                printf("DEBUG: nassociated_P_ %d  dist_mean=%g dist_stddev=%g dist_low/high=%g/%g dist_ratio=%.2f/%.1f exp=%g,  gap_primary=%.1f/%d exp=%g  gap_secondary=%.1f/%d exp=%g  weight_product=%g\n",
                        nassociated_P_within_dist_wt_dist_max, dist_mean, dist_stddev, dist_low, dist_high, distance_ratio, DISTANCE_RATIO_CRITICAL, distance_exponent,
                        gap_primary, gap_primary_critical, gap_primary_exponent, gap_secondary, gap_secondary_critical, gap_secondary_exponent, weight_product);
            }

            if (1) {
                prob *= weight_product; // 20101230 AJL
                if (prob < best_prob) {
                    //if (last_reassociate)
                    //   printf("DEBUG:  REJECTED!\n");
                    continue;
                }
                //if (last_reassociate)
                //    printf("DEBUG:  Accepted\n");
            }
        }*/


        // check if total weight of unique associated stations/phases is greater than min_weight_sum_assoc
        // 20111221 AJL - added
        double weight_sum_assoc_unique;
        if (!checkWtSumUniqueStationsPhases(min_weight_sum_assoc, pdata_list, &weight_sum_assoc_unique, &nassociated_P_work))
            continue;


        // check amplitude attenuation
        // 20140122 AJL -
        LinRegressPower linRegressPower;
        linRegressPower.nvalues = -1;
        linRegressPower.power = -9.9;
        linRegressPower.constant = 0.0;
        linRegressPower.r_square = 0.0;
        linRegressPower.power_dev = 0.0;
        linRegressPower.constant_dev = 0.0;
        double amp_att_weight = 1.0;
        if (use_amplitude_attenuation) {
#define MIN_DIST_FACT_NODE_SIZE_AMP_ATTEN 5.0   // 20160115 AJL - do not asses amp atten if station too close to cell relative to cell size
            int noutliers = -1;
            int noutlier_values = -1;
            double outlier_rate = -9.9;
            int nvalues = calcLinearRegressionPowerRelation(distance_list, num_distances,
                    MIN_DIST_FACT_NODE_SIZE_AMP_ATTEN * critical_node_diagonal * KM2DEG, &linRegressPower);
            // 20160113 AJL  #define ACCEPTABLE_POWER_MAX -0.5
            // 20160113 AJL  #define ACCEPTABLE_POWER_MIN -3.5
#define ACCEPTABLE_POWER_MAX -0.5   // 20160113 AJL  hypolist_persistent_ATTENUATION.ods
#define ACCEPTABLE_POWER_MIN -3.0   // 20160113 AJL  hypolist_persistent_ATTENUATION.ods
#define ACCEPTABLE_OUTLIER_RATE_LOW 0.1
            // 20160128 AJL  #define ACCEPTABLE_OUTLIER_RATE_MAX 0.25
#define ACCEPTABLE_OUTLIER_RATE_MAX 0.5   // 20160128 AJL
            if (nvalues > 0) {
                // calculate weight based on normal cumulative distribution function of linear regression power overlapping with pre-defined acceptable power range
                if (isnan(linRegressPower.power) || isnan(linRegressPower.power_dev)) { // 20150202 AJL
                    ;
                    /*printf("DEBUG: ===> prob NaN!!! [amp_att_weight] volume=%lg   x y x = %g %g %g\n",
                            volume, poct_node->center.y, poct_node->center.x, poct_node->center.z);
                    printf("DEBUG: AMP ATTEN: num_distances %d, nvalues %d, power %f [%f,%f], wt=%f, intercept_1deg %g, r^2 %f\n",
                            num_distances, nvalues, linRegressPower.power, linRegressPower.power - linRegressPower.power_dev, linRegressPower.power + linRegressPower.power_dev,
                            amp_att_weight, linRegressPower.constant, linRegressPower.r_square);*/
                } else {
                    double mean = linRegressPower.power;
                    double deviation = linRegressPower.power_dev;
                    amp_att_weight = cumulDistNormal(ACCEPTABLE_POWER_MAX, mean, deviation) - cumulDistNormal(ACCEPTABLE_POWER_MIN, mean, deviation);
                    // 20160112 TEST
                    amp_att_weight = 0.5 * (1.0 + amp_att_weight); // reduce contribution of amp_att_weight
#if 1
                    // 20160112 TEST
                    outlier_rate = rateOutliersLinearRegressionPowerRelation(distance_list, num_distances,
                            MIN_DIST_FACT_NODE_SIZE_AMP_ATTEN * critical_node_diagonal * KM2DEG, &linRegressPower,
                            AMPLITUDE_ATTENUATION_MIN_ERROR_RATIO, AMPLITUDE_ATTENUATION_MAX_ERROR_RATIO, &noutlier_values, &noutliers);
                    if (outlier_rate > ACCEPTABLE_OUTLIER_RATE_MAX) {
                        amp_att_weight *= 0.5;
                    } else if (outlier_rate > ACCEPTABLE_OUTLIER_RATE_LOW) {
                        // additional weight term is 1.0 at low -> 0.5 at high
                        amp_att_weight *= 1.0 - 0.5 * (outlier_rate - ACCEPTABLE_OUTLIER_RATE_LOW) / (ACCEPTABLE_OUTLIER_RATE_MAX - ACCEPTABLE_OUTLIER_RATE_LOW);
                    }
#endif
                    /*
                                // calculate weight based on overlap of linear regression power with pre-defined acceptable power range
                                double power_min = power - power_dev;
                                double power_max = power + power_dev;
                                if (power_min >= ACCEPTABLE_POWER_MIN && power_max <= ACCEPTABLE_POWER_MAX) { // linear regress power inside acceptable range
                                    amp_att_weight = 1.0;
                                } else if (power_min < ACCEPTABLE_POWER_MIN && power_max > ACCEPTABLE_POWER_MAX) { // linear regress power brackets acceptable range
                                    amp_att_weight = (ACCEPTABLE_POWER_MAX - ACCEPTABLE_POWER_MIN) / (power_max - power_min);
                                } else if (power_min < ACCEPTABLE_POWER_MAX && power_max > ACCEPTABLE_POWER_MIN) { // linear regress power overlaps acceptable range
                                    if (power_min > ACCEPTABLE_POWER_MIN) { // overlap acceptable max boundary
                                        amp_att_weight = (power_max - ACCEPTABLE_POWER_MIN) / (power_max - power_min);
                                    } else { // overlap acceptable min boundary
                                    }
                                    amp_att_weight = (ACCEPTABLE_POWER_MAX - power_min) / (power_max - power_min);
                                } else { // no intersection
                                    amp_att_weight = 0.0;
                                }
                     */
                }
            }
            /*if (1 && last_reassociate && prob > best_prob) {
                printf("DEBUG: AMP ATTEN: num_distances %d, nvalues %d, power %f dev %f [%f,%f], noutliers %d/%d, outlier_rate %f, wt %f, intercept_1deg %g, r^2 %f\n",
                        num_distances, nvalues, linRegressPower.power, linRegressPower.power_dev,
                        linRegressPower.power - linRegressPower.power_dev, linRegressPower.power + linRegressPower.power_dev, noutliers, noutlier_values, outlier_rate,
                        amp_att_weight, linRegressPower.constant, linRegressPower.r_square);
            }*/
            //
            prob *= amp_att_weight; // 20101230 AJL
#ifdef DEBUG_PROB
            if (DEBUG_ALL || last_reassociate) {
                printf("DEBUG_PROB: prob [%f] *= amp_att_weight [%f]\n",
                        prob, amp_att_weight);
            }
#endif
            if (prob < best_prob) {
                //if (last_reassociate)
                //    printf("DEBUG:  REJECTED AMP ATTEN!\n");
                continue; // cannot be better than current best_prob
            }
        }


        // additional weighting of location probability based on azimuth gaps and station counts as a function of distance
        // IMPORTANT: these checks are optimized for global scale / teleseismic location, they should be disabled at other scales

        int nCountClose = 0;
        int nCountFar = 0;
        double gap_weight = 1.0;
        double distanceClose_weight = 1.0;
        double distanceFar_weight = 1.0;

        //#ifdef TTIMES_REGIONAL
        //        if (1) { // DANGER! 20130305 AJL - Regional test only, do not use for global scale / teleseismic location
        //#else
        if (1) {
            //#endif

            //int numCountDistanceNodeSizeTest = 0;
            int numAssociatedClose = 0;
            int numAssociatedFar = 0;

            // check gap - penalty if gap_primary > gap_primary_critical or gap_secondary > gap_secondary_critical
            //gap_primary = calcAzimuthGapSparseArray(AzimuthArray, NumOtimeLimit, nassociated_P_within_dist_wt_dist_max, &gap_secondary);
            gap_primary = calcAzimuthGap(azimuth_list, num_azimuths, &gap_secondary);

            /* 20150722 - remove old gap algorithm here
            int test_gap = 0;
            double gap_primary_exponent = 0.0;
            if (gap_primary > gap_primary_critical) {
                gap_primary_exponent = (gap_primary - gap_primary_critical) / (360.0 - gap_primary_critical);
                test_gap = 1;
            }
            double gap_secondary_exponent = 0.0;
            if (gap_secondary > gap_secondary_critical) {
                gap_secondary_exponent = (gap_secondary - gap_secondary_critical) / (360.0 - gap_secondary_critical);
                test_gap = 1;
            }
            if (test_gap) {
                double gap_exponent = (gap_primary_exponent + gap_secondary_exponent) / 2.0;
                //
                gap_weight = EXP_FUNC(-0.5 * gap_exponent);
                //
                prob *= gap_weight; // 20150202 AJL
                if (prob < best_prob) {
                    //if (last_reassociate)
                    //    printf("DEBUG:  REJECTED GAP!\n");
                    continue; // cannot be better than current best_prob
                }
            }
             */
            // 20150730 - new gap algorithm: weight goes to zero at critical gap angles
#define GAP_BASE 180.0
            if (gap_primary > GAP_BASE || gap_secondary > GAP_BASE) {
                if (gap_primary >= gap_primary_critical || gap_secondary >= gap_secondary_critical) {
                    gap_weight = 0.0;
                } else {
                    double gap_primary_weight = gap_primary > GAP_BASE ? (gap_primary_critical - gap_primary) / (gap_primary_critical - GAP_BASE) : 1.0;
                    double gap_secondary_weight = gap_secondary > GAP_BASE ? (gap_secondary_critical - gap_secondary) / (gap_secondary_critical - GAP_BASE) : 1.0;
                    gap_weight = (gap_primary_weight + gap_secondary_weight) / 2.0; // mean of two gap weights
                }
                // 20160817 AJL - give floor to gap weight, allows small gaps
                gap_weight += 0.5;
                if (gap_weight > 1.0) {
                    gap_weight = 1.0;
                }
                //
                prob *= gap_weight; // 20150202 AJL
#ifdef DEBUG_PROB
                if (DEBUG_ALL || last_reassociate) {
                    printf("DEBUG_PROB: prob [%f] *= gap_weight [%f]\n",
                            prob, gap_weight);
                }
#endif
                if (prob < best_prob) {
                    //if (last_reassociate)
                    //    printf("DEBUG:  REJECTED GAP!\n");
                    continue; // cannot be better than current best_prob
                }
            }


            // check close distance ratio - penalty if (# associated / # available closer than closest associated) < RATIO_NUM_STA_WITHIN_CLOSE_DIST_CRITICAL (3.0)
            double limit_distance = FLT_MAX; // no maximum distance for counting close
            nCountClose = countAllStationsAvailable(COUNT_CLOSE, limit_distance, distance_list, num_distances,
                    lat, lon, channelParameters, &distanceClose, &numAssociatedClose);
            double n_count_close_ratio_exponent = 0.0;
            double n_count_close_ratio = 0.0;
            if (nCountClose > 0) {
                n_count_close_ratio = (double) nassociated_P_within_dist_wt_dist_max / (double) nCountClose;
                if (n_count_close_ratio < RATIO_NUM_STA_WITHIN_CLOSE_DIST_CRITICAL)
                    n_count_close_ratio_exponent = 1.0 - n_count_close_ratio / RATIO_NUM_STA_WITHIN_CLOSE_DIST_CRITICAL;
                //
                distanceClose_weight = EXP_FUNC(-0.5 * n_count_close_ratio_exponent);
                //
                prob *= distanceClose_weight; // 20150202 AJL
#ifdef DEBUG_PROB
                if (DEBUG_ALL || last_reassociate) {
                    printf("DEBUG_PROB: prob [%f] *= distanceClose_weight [%f]\n",
                            prob, distanceClose_weight);
                }
#endif
                if (prob < best_prob) {
                    //if (last_reassociate)
                    //    printf("DEBUG:  REJECTED COUNT_CLOSE!\n");
                    continue; // cannot be better than current best_prob
                }
                //printf("DEBUG: n_count_close_ratio %g = nassociated_P_within_dist_wt_dist_max %d / nCountClose %d\n", n_count_close_ratio, nassociated_P_within_dist_wt_dist_max, nCountClose);
            }

            // check far distance ratio - penalty if (# associated / # available farther than depth) < RATIO_NUM_STA_WITHIN_FAR_DIST_CRITICAL (0.3)
            limit_distance = depth * KM2DEG; // minimum distance is depth for counting far and mean
            nCountFar = countAllStationsAvailable(COUNT_FAR, limit_distance, distance_list, num_distances,
                    lat, lon, channelParameters, &distanceFar, &numAssociatedFar);
            double n_count_far_ratio_exponent = 0.0;
            if (nCountFar > 0) {
                double n_count_far_ratio = (double) numAssociatedFar / (double) nCountFar;
                if (n_count_far_ratio < RATIO_NUM_STA_WITHIN_FAR_DIST_CRITICAL)
                    n_count_far_ratio_exponent = 1.0 - n_count_far_ratio / RATIO_NUM_STA_WITHIN_FAR_DIST_CRITICAL;
                //
                distanceFar_weight = EXP_FUNC(-0.5 * n_count_far_ratio_exponent);
                //
                prob *= distanceFar_weight; // 20101230 AJL
#ifdef DEBUG_PROB
                if (DEBUG_ALL || last_reassociate) {
                    printf("DEBUG_PROB: prob [%f] *= distanceFar_weight [%f]\n",
                            prob, distanceFar_weight);
                }
#endif
                if (prob < best_prob) {
                    //if (last_reassociate)
                    //    printf("DEBUG:  REJECTED COUNT_FAR!\n");
                    continue; // cannot be better than current best_prob
                }
            }
        }


        best_prob = prob;
        best_weight_sum = weight_sum;
        best_effective_cell_size = effective_cell_size;

        if (best_prob > max_node_prob) {
            max_node_prob = best_prob;
            max_node_ot_variance = ot_variance;
            max_node_ot_variance_weight = ot_variance_weight;
        }

        // save raw value
        /*if (poct_node->pdata != NULL) {
         *((double *) poct_node->pdata) = adjusted_weight_sum;
        }*/

        // check if globally best solution
        //if (best_weight_count > *pglobal_best_prob || last_reassociate) {
        //if (!no_reassociate && ((best_prob >= *pglobal_best_prob) || last_reassociate)) {
        if ((best_prob >= pglobalBestValues->best_prob) || last_reassociate) {
#ifdef DEBUG_PROB
            if (DEBUG_ALL || last_reassociate) {
                printf("DEBUG_PROB: if ((best_prob [%f] >= pglobalBestValues->best_prob [%f]) || last_reassociate [%d]\n",
                        best_prob, pglobalBestValues->best_prob, last_reassociate);
                printf("DEBUG_PROB: --------------------------------------\n");
            }
#endif
            //if (weight_sum > *pglobal_best_weight_sum || last_reassociate) {

            //printf("ot_variance_weight %g = time_range_variance %g / ot_variance %g \n", ot_variance_weight, time_range_variance, ot_variance);

            /*static int icount = 0;
            if (icount++ % 1000 == 0) {
                printf("DEBUG: CORE:  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  NumOtimeLimit %d  nassoc_P %d  weight_sum %g  eff_cell_size %g  cumu_weight_sum %g prob %g  ot_var %g  ot_v_f %g crit_size_km %g\n",
                        poct_node->center.y, poct_node->center.x, poct_node->center.z,
                        dsy_global, dsx_global, poct_node->ds.z,
                        NumOtimeLimit, nassociated_P_work, weight_sum, best_effective_cell_size, adjusted_weight_sum, best_prob, ot_variance, ot_variance_weight, critical_node_size_km);
                printf("ot_variance_weight %g = time_error_mean %g * time_error_mean / (ot_variance %g + time_error_mean * time_error_mean)\n", ot_variance_weight, time_range_variance, ot_variance);
                fflush(stdout);
            }*/

            is_global_best = 1;
            int nassociated_P = nassociated_P_work;
            int nassociated = 0;

            if (last_reassociate) {
                // clear associations for this pass and un-associated data
                TimedomainProcessingData* deDataTmp;
                for (j = 0; j < num_de_data; j++) {
                    deData = pdata_list[j];
                    if (skipData(deData, num_pass))
                        continue;
                    // skip if associated in previous pass
                    //if (deData->is_associated && deData->is_associated < num_pass)
                    //    continue;
                    deDataTmp = DataTmpPool + j;
                    clear_deData_assoc(deDataTmp);
                }
                // find currently associated otime limits and set corresponding deDataTmp
                nassociated_P = 0;
                nassociated = 0;
                OtimeLimit *otimeLimit_assoc;
                for (j = 0; j < NumOtimeLimit; j++) {
                    otimeLimit_assoc = *(OtimeLimitList + j);
                    if (otimeLimit_assoc->polarity > 0 // only use lower limit
                            && (otimeLimit_assoc->assoc // if currently associated
                            )) {
                        deDataTmp = DataTmpPool + otimeLimit_assoc->data_id;
                        double residual = otimeLimit_assoc->otime - ot_mean;
                        int count_in_loc = otimeLimit_assoc->count_in_loc;
                        if (deDataTmp->is_associated == num_pass) { // this data already associated for another otime limit in this best update
                            int count_in_loc_existing = count_in_location(deDataTmp->phase_id, deDataTmp->loc_weight, deDataTmp->use_for_location);
                            if (
                                    (!count_in_loc && count_in_loc_existing) // keep existing if this not counted and existing is
                                    || (count_in_loc && count_in_loc_existing && fabs(residual) > fabs(deDataTmp->residual)) // or both counted and existing has smaller residual
                                    || (!count_in_loc && !count_in_loc_existing && fabs(residual) > fabs(deDataTmp->residual)) // or neither counted and existing has smaller residual
                                    )
                                continue;
                            // correct associated counts
                            nassociated--;
                            if (count_in_loc_existing) // was an associate P
                                nassociated_P--;
                        }
                        nassociated++;
                        // will be adding new associated datum
                        if (count_in_loc) {
                            nassociated_P++;
                            deDataTmp->loc_weight = otimeLimit_assoc->total_weight;
                        } else {
                            deDataTmp->loc_weight = 0.0;
                        }
                        //}
                        deDataTmp->use_for_location = otimeLimit_assoc->use_for_location;
                        deDataTmp->is_associated = num_pass;
                        //deDataTmp->is_associated_grid_level = num_grid_level;
                        deDataTmp->epicentral_distance = otimeLimit_assoc->dist;
                        deDataTmp->epicentral_azimuth = otimeLimit_assoc->azimuth;
                        deDataTmp->residual = residual;
                        deDataTmp->dist_weight = otimeLimit_assoc->dist_weight;
                        deDataTmp->polarization.weight = otimeLimit_assoc->polarization_weight;
                        deDataTmp->polarization.azimuth_calc = otimeLimit_assoc->polarization_azimuth_calc;
                        deDataTmp->station_quality_weight = otimeLimit_assoc->quality_weight;
                        strcpy(deDataTmp->phase, phase_name_for_id(otimeLimit_assoc->phase_id));
                        deDataTmp->phase_id = otimeLimit_assoc->phase_id;
                        // set take-off angles  // 20110520 AJL
                        deDataTmp->take_off_angle_inc = get_take_off_angle(deDataTmp->phase_id, deDataTmp->epicentral_distance, depth);
                        deDataTmp->take_off_angle_az = deDataTmp->epicentral_azimuth;
                        /*printf("DEBUG: PHASES %s nassociated_P %d  nassociated %d  flag_a_ref_not_ok %d  data_id %d  residual %g  otimeLimits: %s %s gcd %f dwt %f  ot %f - %f (%f) ot=%f\n",
                                phase_name_for_id(otimeLimit_assoc->phase_id),
                                nassociated_P, nassociated, deData->flag_a_ref_not_ok, otimeLimit_assoc->data_id, residual,
                                pdata_list[otimeLimit_assoc->data_id]->station, phase_name_for_id(otimeLimit_assoc->phase_id), otimeLimit_assoc->dist, otimeLimit_assoc->total_weight,
                                otimeLimit_assoc->time, otimeLimit_assoc->pair->time, otimeLimit_assoc->pair->time - otimeLimit_assoc->time, otimeLimit_assoc->otime);
                         */

                        if (use_amplitude_attenuation) {
                            // set data error ratios
                            deData = otimeLimit_assoc->deData;
                            // check for valid a_ref value
                            if (otimeLimit_assoc->have_amp_atten_values) {
                                double a_ref_calc = linRegressPower.constant * pow(otimeLimit_assoc->dist, linRegressPower.power);
                                if (a_ref_calc > FLT_MIN) {
                                    deDataTmp->amplitude_error_ratio = deData->a_ref / a_ref_calc;
                                }
                            }
                        }
                        deDataTmp->sta_corr = otimeLimit_assoc->sta_corr;


                    }
                }
            }

            // 20110513 AJL - moved from above
            //numCountDistanceNodeSizeTest = countAllStationsAvailable(COUNT_ARITHMETIC_MEAN, limit_distance, distance_list, num_distances,
            //        lat, lon, channelParameters, &distanceNodeSizeTest, &numCountDistanceNodeSizeTest);
            // 20111005 AJL - New node size test algorithm
            int numCountDistanceNodeSizeTest = 2 * (int) (0.5 + min_weight_sum_assoc);
            int num_within_depth = countStationsWithinDistance(depth * KM2DEG, distance_list, num_distances); // number of associated stations within distance = hypo depth
            int station_count = numCountDistanceNodeSizeTest + num_within_depth;
            double distanceNodeSizeTest = DEG2KM * findDistanceForStationCount(station_count, distance_list, num_distances); // distance within which there are 2*MIN_NUMBER_ASSOCIATE_OCT_TREE associated stations
            //printf("DEBUG: >>> distanceNodeSizeTest: numCountDistanceNodeSizeTest %d, num_within_depth %d, station_count %d, distanceNodeSizeTest %g, num_distances %d\n",
            //        numCountDistanceNodeSizeTest, num_within_depth, station_count, distanceNodeSizeTest, num_distances);

            // update global best reference values
            pglobalBestValues->node_prob = max_node_prob;
            pglobalBestValues->node_ot_variance = max_node_ot_variance;
            pglobalBestValues->best_ot_mean = reftime + ot_mean;
            pglobalBestValues->best_ot_variance = ot_variance;
            pglobalBestValues->best_lat = lat;
            pglobalBestValues->best_lon = lon;
            pglobalBestValues->best_depth = depth;
            pglobalBestValues->best_nassociated_P = nassociated_P;
            pglobalBestValues->best_nassociated = nassociated;
            pglobalBestValues->best_weight_sum = weight_sum;
            pglobalBestValues->best_prob = best_prob;
            pglobalBestValues->node_ot_variance_weight = max_node_ot_variance_weight;
            pglobalBestValues->effective_cell_size = best_effective_cell_size;
            pglobalBestValues->nCountClose = nCountClose;
            pglobalBestValues->nCountFar = nCountFar;
            pglobalBestValues->distanceClose = distanceClose;
            pglobalBestValues->distanceFar = distanceFar;
            pglobalBestValues->distanceNodeSizeTest = distanceNodeSizeTest;
            pglobalBestValues->numCountDistanceNodeSizeTest = numCountDistanceNodeSizeTest;
            // amplitude attenuation
            pglobalBestValues->linRegressPower = linRegressPower;
            // qualityTestsWeights
            pglobalBestValues->qualityIndicators.ot_variance_weight = ot_variance_weight;
            pglobalBestValues->qualityIndicators.weight_sum_assoc_unique = weight_sum_assoc_unique;
            pglobalBestValues->qualityIndicators.amp_att_weight = amp_att_weight;
            pglobalBestValues->qualityIndicators.gap_weight = gap_weight;
            pglobalBestValues->qualityIndicators.distanceClose_weight = distanceClose_weight;
            pglobalBestValues->qualityIndicators.distanceFar_weight = distanceFar_weight;

        }
    }
    //printf("1 ============> channelParameters addr %ld\n", (long) channelParameters);

    // clear OtimeLimitList
    NumOtimeLimit = 0;

    // put node results in result tree
    if (best_prob != PROB_UNSET) {
        // update node value and result tree with best found for this node location
        poct_node->value = best_prob; // this is quasi- VALUE_IS_PROB_DENSITY_IN_NODE (pdf) (double prob = ot_variance_weight * adjusted_weight_sum / effective_cell_size;)
        if (!isnan(best_prob)) {
#ifdef PURE_OCTREE
            double result_tree_value;
            // effective_cell_size is based on origin time range (otime_max - otime_min), much larger than truer char length of cell, so correct to get truer pdf
            //result_tree_value = best_prob * best_effective_cell_size; // result tree is sorted on probability in node
            //
            // VALUE_IS_PROB_DENSITY_IN_NODE
            //double scatter_pdf_volume = pow(volume, 1.0 / 3.0); // truer characteristic length of cell, will be used for assoc and loc scatter samples
            // effective_cell_size is based on origin time range (otime_max - otime_min), much larger than truer char length of cell, so correct to get truer pdf
            //poct_node->value = best_prob * (best_effective_cell_size / scatter_pdf_volume); // need quasi- VALUE_IS_PROB_DENSITY_IN_NODE (pdf) for assoc and loc scatter samples
            //
            // VALUE_IS_PROBABILITY_IN_NODE
            //poct_node->value = result_tree_value; // value is probability in node
            //double scatter_pdf_volume = 1.0;    // volume not used for VALUE_IS_PROBABILITY_IN_NODE
            //
            // VALUE_IS_PROBABILITY_IN_NODE
            //poct_node->value = best_prob; // value is probability in node
            //double scatter_pdf_volume = 1.0;    // volume not used for VALUE_IS_PROBABILITY_IN_NODE
            //
            // VALUE_IS_PROBABILITY_IN_NODE
            if (no_reassociate) { // location
                //poct_node->value = max_node_ot_variance_weight * pow(volume, 1.0 / 3.0); // value is probability in node
                // effective_cell_size is based on origin time range (otime_max - otime_min), much larger than truer char length of cell, so correct to get truer pdf
                result_tree_value = best_prob * best_effective_cell_size; // result tree is sorted on probability in node
                //result_tree_value *= best_effective_cell_size; // 20150907 AJL - slow down associator, more thorough search
                //poct_node->value = best_prob * best_effective_cell_size * pow(volume, 1.0 / 3.0); // value is probability in node
                poct_node->value = best_prob * best_effective_cell_size; // value is probability in node
            } else { // association
                poct_node->value = best_prob * best_effective_cell_size; // value is probability in node
                // effective_cell_size is based on origin time range (otime_max - otime_min), much larger than truer char length of cell, so correct to get truer pdf
                result_tree_value = best_prob * best_effective_cell_size; // result tree is sorted on probability in node
                //result_tree_value *= best_effective_cell_size; // 20150907 AJL - slow down associator, more thorough search
                //result_tree_value = best_prob * pow(volume, 1.0 / 3.0); // 20150907 AJL - slow down associator, more thorough search
            }
            double scatter_pdf_volume = pow(volume, 1.0 / 3.0); // volume not used for VALUE_IS_PROBABILITY_IN_NODE
            //
            *ppResultTreeRoot = addResult(*ppResultTreeRoot, result_tree_value, scatter_pdf_volume, poct_node);
            //VOLUME*ppResultTreeRoot = addResult(*ppResultTreeRoot, best_prob, volume, poct_node);
#else
            *ppResultTreeRoot = addResult(*ppResultTreeRoot, best_prob, volume, poct_node);
#endif
        } else {
            printf("===================== addResult NaN!!! poct_node->value=%lg volume=%lg   x y x = %g %g %g\n",
                    poct_node->value, volume, poct_node->center.y, poct_node->center.x, poct_node->center.z);
        }
        /*static int icount = 0;
        if ((long) channelParameters == 0 || icount % 1 == 0) {
            printf("==========poct_node->value=%lg volume=%lg log_value_volume=%g  x y z = %g %g %g\n",
                    poct_node->value, volume, log_value_volume, poct_node->center.y, poct_node->center.x, poct_node->center.z);
        }
        icount++;*/
        /*printf("DEBUG: addResult:  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassoc_P %d  best_weight_sum %g   best_prob %g  ot_var %g  ot_v_f %g crit_size_km %g\n",
                poct_node->center.y, poct_node->center.x, poct_node->center.z,
                dsy_global, dsx_global, poct_node->ds.z,
                nassociated_P_work, best_weight_sum, best_prob, ot_variance, ot_variance_weight, critical_node_size_km);
        fflush(stdout);*/
    }

    return (is_global_best);

}


#define ASSOC_DIST_HISTOGRAM_BIN_MIN 0.0
#define ASSOC_DIST_HISTOGRAM_BIN_WIDTH 1.0
#define ASSOC_DIST_HISTOGRAM_BIN_NUM 90
static int assoc_dist_histogram_count[ASSOC_DIST_HISTOGRAM_BIN_NUM]; // bins of 10deg width from 0-90deg
#define CRITICAL_GAP_MIN 30     // 30deg
#define CRITICAL_GAP_MAX 30     // 30deg

/** check for and un-associate outlier associated stations (e.g. isolated, very distant stations associated with small, local event.)
 */
// 20121126 AJL - added

int unassociateOutlierStations(TimedomainProcessingData** pdata_list, int num_de_data, int num_pass, double best_depth, double best_ot_mean) {

    TimedomainProcessingData* deData;

    // clear bins
    int nbin;
    for (nbin = 0; nbin < ASSOC_DIST_HISTOGRAM_BIN_NUM; nbin++)
        assoc_dist_histogram_count[nbin] = 0;

    TimedomainProcessingData* deDataTmp;

    // loop through data associated for this pass and accumulate distance bin count of P arrivals
    int j;
    for (j = 0; j < num_de_data; j++) {
        deData = pdata_list[j];
        if (skipData(deData, num_pass))
            continue;
        deDataTmp = DataTmpPool + j;
        if (deDataTmp->is_associated == num_pass && is_P(deDataTmp->phase_id)) {
            nbin = (int) ((deDataTmp->epicentral_distance - ASSOC_DIST_HISTOGRAM_BIN_MIN) / ASSOC_DIST_HISTOGRAM_BIN_WIDTH);
            if (nbin < ASSOC_DIST_HISTOGRAM_BIN_NUM)
                assoc_dist_histogram_count[nbin]++;
            //printf("DEBUG: >>>> unassociateOutlierStations: ADD P: num_pass %d, nbin %d, %s %s %s %s, dist %f\n",
            //        num_pass, nbin, pdata_list[j]->network, pdata_list[j]->station, pdata_list[j]->channel, deDataTmp->phase, deDataTmp->epicentral_distance);
        }
    }

    // check for gap larger than critical gap
    double distance_limit = 0.0;
    int nempty = 0;
    int critical_gap = CRITICAL_GAP_MIN;
    for (nbin = 0; nbin < ASSOC_DIST_HISTOGRAM_BIN_NUM; nbin++) {
        if (assoc_dist_histogram_count[nbin] == 0) {
            nempty++;
            if (nempty >= critical_gap) {
                break;
            }
        } else {
            nempty = 0;
            distance_limit = ASSOC_DIST_HISTOGRAM_BIN_MIN + (double) (nbin + 1) * ASSOC_DIST_HISTOGRAM_BIN_WIDTH;
            /*critical_gap = nbin / 2; // critical gap is 1/2 current largest distance with assoc P
            if (critical_gap < CRITICAL_GAP_MIN)
                critical_gap = CRITICAL_GAP_MIN;
            else if (critical_gap > CRITICAL_GAP_MAX)
                critical_gap = CRITICAL_GAP_MAX;*/
        }
    }

    // set distance and time limits
    //distance_limit = ASSOC_DIST_HISTOGRAM_BIN_MIN + (double) (nbin + 1) * ASSOC_DIST_HISTOGRAM_BIN_WIDTH;
    double stime_limit = get_ttime(get_S_phase_index(), distance_limit, best_depth);
    double time_limit = best_ot_mean + stime_limit;
    //printf("DEBUG: >>>> unassociateOutlierStations: num_pass %d, nempty %d, nbin_end_gap %d, distance_limit %f, stime_limit %f, best_ot_mean %f, time_limit %f\n",
    //        num_pass, nempty, nbin, distance_limit, stime_limit, best_ot_mean, time_limit);

    // check if critical gap found, un-associate data for distance or time beyond limit
    int num_removed = 0;
    if (nbin < ASSOC_DIST_HISTOGRAM_BIN_NUM && stime_limit > 0.0) {
        for (j = 0; j < num_de_data; j++) {
            deData = pdata_list[j];
            if (skipData(deData, num_pass))
                continue;
            deDataTmp = DataTmpPool + j;
            if (deDataTmp->is_associated == num_pass
                    && !count_in_location(deDataTmp->phase_id, 999.9, 1) // do not use phases used for location (also prevents changing location probability, number associated)
                    && (deDataTmp->epicentral_distance > distance_limit || ((double) pdata_list[j]->t_time_t + pdata_list[j]->t_decsec) > time_limit)
                    ) {
                //printf("DEBUG: >>>> unassociateOutlierStations: REMOVE: %s %s %s %s, dist %f, time %f, wt %f\n",
                //        pdata_list[j]->network, pdata_list[j]->station, pdata_list[j]->channel, deDataTmp->phase, deDataTmp->epicentral_distance,
                //        ((double) pdata_list[j]->t_time_t + pdata_list[j]->t_decsec) - best_ot_mean, deDataTmp->loc_weight);
                num_removed++;
                clear_deData_assoc(deDataTmp);
            }
        }
    }

    return (num_removed);

}

/** free oct-tree memory
 */
void freeOctTreeStructures(ResultTreeNode* pResultTreeRoot, Tree3D * pOctTree) {

    int freeDataPointer = 1;
    freeResultTree(pResultTreeRoot);
    pResultTreeRoot = NULL;
    freeTree3D(pOctTree, freeDataPointer);
    pOctTree = NULL;
}

void setBestHypocenter(HypocenterDesc *hypocenter, double x_vector_sum, // x accumulation for vector sum of epicenter to station vectors
        double y_vector_sum, // y accumulation for vector sum of epicenter to station vectors
        double best_depth_step,
        double gcd_min, double gcd_max,
        OctNode* poct_node,
        double current_critical_node_size_km,
        GlobalBestValues *pglobalBestValues
        ) {

    // set mean vector of epicenter to station vectors
    hypocenter->vector_azimuth = atan2(x_vector_sum, y_vector_sum) * RA2DE;
    if (hypocenter->vector_azimuth < 0.0)
        hypocenter->vector_azimuth += 360.0;
    hypocenter->vector_distance = sqrt(x_vector_sum * x_vector_sum + y_vector_sum * y_vector_sum);

    // set azimuth gap
    hypocenter->gap_primary = calcAzimuthGap(azimuth_list, num_azimuths, &(hypocenter->gap_secondary));

    //hypocenter->ot = (time_t) ((long) best_ot_mean);
    hypocenter->otime = pglobalBestValues->best_ot_mean;
    hypocenter->ot_std_dev = pglobalBestValues->best_ot_variance < 0.0 ? 0.0 : sqrt(pglobalBestValues->best_ot_variance);
    hypocenter->prob = pglobalBestValues->best_prob;
    hypocenter->weight_sum = pglobalBestValues->best_weight_sum;
    hypocenter->lat = pglobalBestValues->best_lat;
    if (pglobalBestValues->best_lon < -180.0)
        pglobalBestValues->best_lon += 360.0;
    else if (pglobalBestValues->best_lon > 180.0)
        pglobalBestValues->best_lon -= 360.0;
    hypocenter->lon = pglobalBestValues->best_lon;
    hypocenter->depth = pglobalBestValues->best_depth;
    hypocenter->depth_step = best_depth_step;
    hypocenter->nassoc_P = pglobalBestValues->best_nassociated_P;
    hypocenter->nassoc = pglobalBestValues->best_nassociated;
    hypocenter->dist_min = gcd_min;
    hypocenter->dist_max = gcd_max;
    // amplitude attenuation
    hypocenter->linRegressPower = pglobalBestValues->linRegressPower;
    // hypocenter quality
    hypocenter->qualityIndicators = pglobalBestValues->qualityIndicators;
    // global best solution oct-tree node
    hypocenter->global_best_oct_node = *poct_node;
    hypocenter->global_best_oct_node.pdata = NULL; // IMPORTANT: pdata not used, if used in future must handle copying of pdata and setting pointer here!
    hypocenter->global_best_critical_node_size_km = current_critical_node_size_km;

}

void setDataAssociationInformation(HypocenterDesc *hypocenter, TimedomainProcessingData** pdata_list, int num_de_data, int num_pass,
        int definitive, int reAssociateOnly, double best_depth_step,
        OctNode* poct_node, double current_critical_node_size_km) {

    TimedomainProcessingData* deData;
    TimedomainProcessingData* deDataTmp;

    double gcd_min = 999.9;
    double gcd_max = -1.0;
    // prepare azimuth array
    if (definitive) {
        if (sizeAzimuthArray < num_de_data) {
            sizeAzimuthArray = num_de_data;
            AzimuthArray = (ValueDesc*) realloc(AzimuthArray, sizeAzimuthArray * sizeof (ValueDesc));
            if (AzimuthArray == NULL)
                fprintf(stderr, "ERROR: realloc of azimuth_array: %s\n", strerror(errno));
        }
    }
    double x_vector_sum = 0.0; // x accumulation for vector sum of epicenter to station vectors
    double y_vector_sum = 0.0; // y accumulation for vector sum of epicenter to station vectors
    int n_vector_sum = 0;
    int num_az_gap = 0;
    if (definitive) {
        free_ValueDescList(&azimuth_list, &num_azimuths);
    }
    int j;
    for (j = 0; j < num_de_data; j++) {
        deData = pdata_list[j];
        if (skipData(deData, num_pass))
            continue;
        // skip if associated in previous pass
        //if (deData->is_associated && deData->is_associated < num_pass)
        //   continue;
        // skip if not associated in previous gridding level of this pass
        //if (num_grid_level > 9999 && deData->is_associated != num_pass)
        //    continue;
        deDataTmp = DataTmpPool + j;
        if (deDataTmp->is_associated == num_pass) {
            // set association for this pass and grid level
            deData->is_associated = deDataTmp->is_associated;
            if (reAssociateOnly) {
                if (deData->is_full_assoc_loc == -1) { // not set, avoid flagging data previously associated with full assoc/loc as associated with re-association only
                    deData->is_full_assoc_loc = 0;
                }
            } else {
                deData->is_full_assoc_loc = 1;
            }
            deData->epicentral_distance = deDataTmp->epicentral_distance;
            deData->epicentral_azimuth = deDataTmp->epicentral_azimuth;
            deData->residual = deDataTmp->residual;
            deData->dist_weight = deDataTmp->dist_weight;
            deData->polarization.weight = deDataTmp->polarization.weight;
            deData->polarization.azimuth_calc = deDataTmp->polarization.azimuth_calc;
            if (deData->is_full_assoc_loc == 1) {
                deData->loc_weight = deDataTmp->loc_weight; // location weight non-zero only if data previously associated with full assoc/loc
            } else if (reAssociateOnly) {
                deData->loc_weight = 0.0; // location weight zero if associated by re-association only
            }
            strcpy(deData->phase, deDataTmp->phase);
            deData->phase_id = deDataTmp->phase_id;
            deData->take_off_angle_inc = deDataTmp->take_off_angle_inc;
            deData->take_off_angle_az = deDataTmp->take_off_angle_az;
            if (count_in_location(deData->phase_id, deData->dist_weight, deData->use_for_location)) {
                // save min distance
                if (deData->epicentral_distance < gcd_min)
                    gcd_min = deData->epicentral_distance;
                if (deData->epicentral_distance > gcd_max)
                    gcd_max = deData->epicentral_distance;
                x_vector_sum += deData->epicentral_distance * sin(deData->epicentral_azimuth * DE2RA);
                y_vector_sum += deData->epicentral_distance * cos(deData->epicentral_azimuth * DE2RA);
                n_vector_sum++;
                // save azimuth info
                if (definitive) {
                    AzimuthArray[num_az_gap].value = deData->epicentral_azimuth;
                    addValueDescToValueList(AzimuthArray + num_az_gap, &azimuth_list, &num_azimuths);
                    num_az_gap++;
                }
            }
            // amplitude attenuation
            deData->amplitude_error_ratio = deDataTmp->amplitude_error_ratio;
            deData->sta_corr = deDataTmp->sta_corr;
            //} else if (deData->is_associated == num_pass && deData->is_associated_grid_level < num_grid_level) {
        } else {
            // clear association for this and later passes
            clear_deData_assoc(deData);
            if (!reAssociateOnly) {
                deData->is_full_assoc_loc = -1; // remaining data is later passes and necessarily will be full assoc/loc or unassociated
            }

        }
    }

    if (definitive) {
        // set fields in Hypocenter for best location
        setBestHypocenter(hypocenter, x_vector_sum, // x accumulation for vector sum of epicenter to station vectors
                y_vector_sum, // y accumulation for vector sum of epicenter to station vectors
                best_depth_step,
                gcd_min, gcd_max,
                poct_node,
                current_critical_node_size_km,
                &globalBestValues
                );
    }
}

Tree3D *setUpOctTree(double lat_min, double lat_max, double lat_step,
        double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step) {

    int num_lat = (int) ((lat_max - lat_min + lat_step * 0.001) / lat_step);
    int num_lon_largest = (int) ((lon_max - lon_min + lon_step_smallest * 0.001) / lon_step_smallest);
    int num_depth = (int) ((depth_max - depth_min + depth_step * 0.001) / depth_step);

    // set up oct-tree x, y, z grid
    Tree3D *pOctTree;
    double integral = 0.0;
    void *pdata = NULL;
    pOctTree = newTree3D_spherical(OCTTREE_TYPE, num_lon_largest,
            num_lat, num_depth,
            lon_min, lat_min, depth_min,
            lon_step_smallest, lat_step, depth_step, OCTTREE_UNDEF_VALUE, integral, pdata);

    return (pOctTree);
}

/** oct-tree global search association and location
 *
 * returns weighted count of associated stations
 *
 */

//#define NUM_ASSOC_SCATTER_SAMPLE 2048
#define NUM_ASSOC_SCATTER_SAMPLE 1024   // 20150401 AJL

double octtreeGlobalAssociationLocation_full(int num_pass, double min_weight_sum_assoc, int max_num_nodes,
        double nominal_critical_node_size_km, double min_critical_node_size_km, double nominal_min_node_size_km,
        double gap_primary_critical, double gap_secondary_critical,
        double lat_min, double lat_max, double lat_step, double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step, int is_local_search,
        int numPhaseTypesUse, int phaseTypesUse[get_num_ttime_phases()], char channelNamesUse[get_num_ttime_phases()][128], double timeDelayUse[MAX_NUM_TTIME_PHASES][2],
        double reference_phase_ttime_error,
        TimedomainProcessingData** pdata_list, int num_de_data,
        HypocenterDesc *hypocenter,
        float **passoc_scatter_sample, int *p_n_alloc_assoc_scatter_sample, int i_get_scatter_sample, int *pn_scatter_sample, double *p_global_max_nassociated_P_lat_lon,
        ChannelParameters * channelParameters,
        time_t time_min, time_t time_max
        ) {


    // exp funct test
    /*
    if (1) {
        printf("EXP FUNC TEST:\n");
        double yarg;
        for (yarg = pow(2.0, 10); yarg > pow(2.0, -10); yarg /= 2.0)
            printf("   exp(%g) = %g\n", yarg, EXP_FUNC(yarg));
        for (yarg = pow(2.0, -10); yarg < pow(2.0, 10); yarg *= 2.0)
            printf("   exp(%g) = %g\n", -yarg, EXP_FUNC(-yarg));
    }*/

    int last_reassociate = 0; // if = 1, will check all phases for association, not just P etc. counted in location.  should use only at optimal hypocenter...
    double last_reassociate_otime = -1.0; // set to a value > 0.0 to only uses phases that give predicted otime range that includes hypocenter->otime
    double last_reassociate_ot_std_dev = -1.0; // must be >= 0.0 if last_reassociate_otime > 0.0
    int try_assoc_remaining_definitive = 0; // if = 1, will check definitive phases (e.g. P, PKP) for association using increased ttime error
    int no_reassociate = 0; // if = 1 will use last reassociate phases only, and will only calculate solution quality and update ResultTree

    // 20160919 AJL - for local search around existing hypo, set allowed otime range around hypocenter->otime
    if (is_local_search) {
        last_reassociate_otime = hypocenter->otime; // set to a value > 0.0 to only uses phases that give predicted otime range that includes hypocenter->otime
        last_reassociate_ot_std_dev = hypocenter->ot_std_dev; // must be >= 0.0 if last_reassociate_otime > 0.0
        //printf("DEBUG: is_local_search: otime %.1f +- %.2f\n",   last_reassociate_otime, last_reassociate_ot_std_dev);
    }

    // need at least 2 stations for location
    if (num_de_data < 2) {
        return (-1);
    }

    int n_ot_limit, j;
    TimedomainProcessingData* deData;

    if (depth_min < 0.0)
        depth_min = 0.0;
    if (depth_max > get_dist_time_depth_max())
        depth_max = get_dist_time_depth_max();

    if (lat_min < -90.0)
        lat_min = -90.0;
    if (lat_max > 90.0)
        lat_max = 90.0;

    /* 20160913 AJL  int num_lat = (int) ((lat_max - lat_min + lat_step * 0.001) / lat_step);
    int num_lon_largest = (int) ((lon_max - lon_min + lon_step_smallest * 0.001) / lon_step_smallest);
    int num_depth = (int) ((depth_max - depth_min + depth_step * 0.001) / depth_step);
    // set up oct-tree x, y, z grid
    Tree3D *pOctTree;
    double integral = 0.0;
    void *pdata = NULL;
    pOctTree = newTree3D_spherical(OCTTREE_TYPE, num_lon_largest,
            num_lat, num_depth,
            lon_min, lat_min, depth_min,
            lon_step_smallest, lat_step, depth_step, OCTTREE_UNDEF_VALUE, integral, pdata);
     */
    Tree3D *pOctTree = setUpOctTree(lat_min, lat_max, lat_step,
            lon_min, lon_max, lon_step_smallest,
            depth_min, depth_max, depth_step);

    // initialize array for scatter sample
    if (i_get_scatter_sample && *passoc_scatter_sample != NULL && (*p_n_alloc_assoc_scatter_sample < NUM_ASSOC_SCATTER_SAMPLE)) {
        free(*passoc_scatter_sample);
        *passoc_scatter_sample = NULL;
    }
    if (i_get_scatter_sample && *passoc_scatter_sample == NULL) {
        *passoc_scatter_sample = (float*) calloc(4 * NUM_ASSOC_SCATTER_SAMPLE, sizeof (float));
        *p_n_alloc_assoc_scatter_sample = NUM_ASSOC_SCATTER_SAMPLE;
    }
    *pn_scatter_sample = 0;


    // allocate OtimeLimit pool array
    int size = 2 * num_de_data * get_num_ttime_phases();
    if (sizeOtimeLimitPool < size) {
        OtimeLimitPool = (OtimeLimit*) realloc(OtimeLimitPool, size * sizeof (OtimeLimit));
        if (OtimeLimitPool == NULL)
            fprintf(stderr, "ERROR: realloc of OtimeLimitPool: %s\n", strerror(errno));
        sizeOtimeLimitPool = size;
    }

    // allocate deDataTmp pool array
    TimedomainProcessingData* deDataTmp = NULL;
    if (sizeDataTmpPool < num_de_data) {
        DataTmpPool = (TimedomainProcessingData*) realloc(DataTmpPool, num_de_data * sizeof (TimedomainProcessingData));
        if (DataTmpPool == NULL)
            fprintf(stderr, "ERROR: realloc of DataTmpPool: %s\n", strerror(errno));
        sizeDataTmpPool = num_de_data;
        // initialize data
        for (j = 0; j < num_de_data; j++) {
            deDataTmp = DataTmpPool + j;
            clear_deData_assoc(deDataTmp);
        }
    }

    // initialize global best solution variables
    globalBestValues.best_ot_variance = -1.0;
    globalBestValues.best_ot_mean = 0;
    globalBestValues.best_lat = 0.0;
    globalBestValues.best_lon = 0.0;
    globalBestValues.best_depth = 0.0;
    globalBestValues.node_ot_variance_weight = 0;
    globalBestValues.best_nassociated_P = 0;
    globalBestValues.best_nassociated = 0;
    globalBestValues.best_weight_sum = FLT_MIN; // non-zero to prevent divide by zero
    globalBestValues.best_prob = -FLT_MAX;
    globalBestValues.effective_cell_size = -1.0;
    globalBestValues.nCountClose = 0;
    globalBestValues.numCountDistanceNodeSizeTest = 0;
    globalBestValues.nCountFar = 0;
    globalBestValues.distanceClose = 0.0;
    globalBestValues.distanceNodeSizeTest = 0.0;
    globalBestValues.distanceFar = 0.0;


    // set min cutoff below which distance weight = 1.0
#ifdef COMPILE_OLD_VERSION
    // 20150225 AJL - remove this algorithm, it changes weighting as function of available picks, always varying and may be somewhat arbitrary
    // use mean of maximum distance between each station and all other stations, helps to stabilize location when only a few nearby station available
    double distance_weight_dist_min_current = getMeanDistanceMax(pdata_list, num_de_data);
    // check that this distance does not exceed DISTANCE_WEIGHT_DIST_MIN
    printf("=====================> 1 DEBUG: distance_weight_dist_min_current %f\n", distance_weight_dist_min_current);
    if (distance_weight_dist_min_current > DISTANCE_WEIGHT_DIST_MIN)
        distance_weight_dist_min_current = DISTANCE_WEIGHT_DIST_MIN;
    printf("=====================> 2 DEBUG: distance_weight_dist_min_current %f\n", distance_weight_dist_min_current);
#else
    // 20150225 AJL - set constant dist min
    double distance_weight_dist_min_current = DISTANCE_WEIGHT_DIST_MIN;
#endif


    // clear data association and perform checks if first pass
    if (num_pass == 1) {
        for (n_ot_limit = 0; n_ot_limit < num_de_data; n_ot_limit++) {
            // clear data association fields
            deData = pdata_list[n_ot_limit];
            clear_deData_assoc(deData);
            // if flagged as use_for_location but will be skipped, check if twin stream pick can  be used for location
            // 20121130 AJL - added
            if (deData->use_for_location && skipData(deData, -999)) {
                TimedomainProcessingData* deDataTwin = deData->use_for_location_twin_data;
                if (deDataTwin != NULL && !skipData(deDataTwin, -999)) {
                    deDataTwin->use_for_location = 1;
                    deData->use_for_location = 0;
                    //printf("DEBUG: >>>> set twin to use for loc: num_pass %d: skipped: %s %s %s %d:%d:%f %s %s -> use: %s %s %s %d:%d:%f %s %s\n",
                    //        num_pass, deData->network, deData->station, deData->channel, deData->hour, deData->min, deData->dsec, pick_stream_name(deData), deData->phase,
                    //        deDataTwin->network, deDataTwin->station, deDataTwin->channel, deData->hour, deData->min, deData->dsec, pick_stream_name(deDataTwin), deDataTwin->phase
                    //        );
                }
            }
            // amplitude attenuation
            deData->amplitude_error_ratio = -1.0;
            deData->sta_corr = 0.0;
        }
    }

    *p_global_max_nassociated_P_lat_lon = 0.0;

    ResultTreeNode* pResultTreeRoot = NULL; // Octtree likelihood*volume results tree root node
    OctNode* poct_node = NULL;

    double smallest_node_size_divide_y = -FLT_MAX;

    // get node size in km
    poct_node = pOctTree->nodeArray[0][0][0];
    smallest_node_size_divide_y = poct_node->ds.y * DEG2KM;
    /*
    double critical_node_size_km = poct_node->ds.y * DEG2KM;
    while (critical_node_size_km > nominal_critical_node_size_km)
        critical_node_size_km /= 2.0;
    critical_node_size_km += FLT_MIN;
     */

    int nSamples = 0;
    int last_info_nsamples = 0;
    long info_nsamples_start_time = clock();

    int is_global_best = 0;
    int found_global_best_during_divide = 0;

    double min_prob_assoc = -FLT_MAX;
    double min_prob_assoc_definative = min_weight_sum_assoc - 1.0;
    //min_prob_assoc = 2.0; // 20101226
    min_prob_assoc = min_prob_assoc_definative; // 20101226B

    globalBestValues.node_prob = 0.0;
    globalBestValues.node_ot_variance = 0.0;
    double oct_node_volume;



    // Stage 1: evaluate each cell in Tree3D =================================================================

    // allocate initial regular grid distance, azimuth and travel-time arrays
    // required initial regular grid distances, azimuths and travel-time arrays are the same for all passes of associate locate
    if (!is_local_search) {
        int initial_grid_lat_lon_size = pOctTree->numx * pOctTree->numy;
        if (USE_SAVED_INITIAL_GRIDS && SavedInitialGridDistances == NULL) {
            printf("Info: Allocation of SavedInitialGridDistances: size: %.1fkB, initial_grid_lat_lon_size: %d\n",
                    (double) (initial_grid_lat_lon_size * sizeof (double)) / 1.0e3, initial_grid_lat_lon_size);
            SavedInitialGridDistances = (double*) calloc(initial_grid_lat_lon_size * MAX_NUM_SOURCES, sizeof (double));
            if (SavedInitialGridDistances == NULL)
                fprintf(stderr, "ERROR: calloc of SavedInitialGridDistances: %s\n", strerror(errno));
            int n, m;
            for (n = 0; n < initial_grid_lat_lon_size; n++) {
                for (m = 0; m < MAX_NUM_SOURCES; m++) {
                    *(SavedInitialGridDistances + n * MAX_NUM_SOURCES + m) = -1.0;
                }
            }
            sizeSavedInitialGridDistances = initial_grid_lat_lon_size;
        }
        if (USE_SAVED_INITIAL_GRIDS && SavedInitialGridAzimuths == NULL) {
            printf("Info: Allocation of SavedInitialGridAzimuths: size: %.1fkB, initial_grid_lat_lon_size: %d\n",
                    (double) (initial_grid_lat_lon_size * sizeof (double)) / 1.0e3, initial_grid_lat_lon_size);
            SavedInitialGridAzimuths = (double*) calloc(initial_grid_lat_lon_size * MAX_NUM_SOURCES, sizeof (double));
            if (SavedInitialGridAzimuths == NULL)
                fprintf(stderr, "ERROR: calloc of SavedInitialGridAzimuths: %s\n", strerror(errno));
            int n, m;
            for (n = 0; n < initial_grid_lat_lon_size; n++) {
                for (m = 0; m < MAX_NUM_SOURCES; m++) {
                    *(SavedInitialGridAzimuths + n * MAX_NUM_SOURCES + m) = -1.0;
                }
            }
            sizeSavedInitialGridAzimuths = initial_grid_lat_lon_size;
        }
        int initial_grid_depth_size = pOctTree->numz;
        if (USE_SAVED_INITIAL_GRIDS && SavedInitialGridTravelTimes == NULL) {
            printf("Info: Allocation of SavedInitialGridTravelTimes: size: %.1fkB, initial_grid_lat_lon_size: %d, initial_grid_depth_size: %d\n",
                    (double) (initial_grid_lat_lon_size * initial_grid_depth_size * MAX_NUM_SOURCES * 2 * sizeof (double)) / 1.0e3,
                    initial_grid_lat_lon_size, initial_grid_depth_size);
            SavedInitialGridTravelTimes = (double*) calloc(initial_grid_lat_lon_size * initial_grid_depth_size * MAX_NUM_SOURCES * 2, sizeof (double***));
            if (SavedInitialGridTravelTimes == NULL)
                fprintf(stderr, "ERROR: calloc of SavedInitialGridTravelTimes: %s\n", strerror(errno));
            int n, i, m;
            for (n = 0; n < initial_grid_lat_lon_size; n++) {
                for (i = 0; i < initial_grid_depth_size; i++) {
                    for (m = 0; m < MAX_NUM_SOURCES; m++) {
                        *(SavedInitialGridTravelTimes + n * initial_grid_depth_size * MAX_NUM_SOURCES * 2 + i * MAX_NUM_SOURCES * 2 + m * 2) = -1.0;
                    }
                }
            }
            sizeSavedInitialGridTravelTimeLatLon = initial_grid_lat_lon_size;
            sizeSavedInitialGridTravelTimeDepth = initial_grid_depth_size;
        }
    }

    // DEBUG
    DEBUG_n_not_available = 0;
    DEBUG_n_set = 0;
    DEBUG_n_read = 0;

    int ilat, ilon, idepth;
    int indexLatLonSavedInitialGrid = 0;
    int indexDepthSavedInitialGrid = 0;
    if (is_local_search) {
        // 20160912 AJL - for local search around existing hypo, search grid will be shifted from full, initial grid
        indexLatLonSavedInitialGrid = -1;
        indexDepthSavedInitialGrid = -1;
    }

    globalBestValues.best_prob = -FLT_MAX;

    for (ilon = 0; ilon < pOctTree->numx; ilon++) {
        for (ilat = 0; ilat < pOctTree->numy; ilat++) {
            indexDepthSavedInitialGrid = 0;
            for (idepth = 0; idepth < pOctTree->numz; idepth++) {

                poct_node = pOctTree->nodeArray[ilon][ilat][idepth];
                if (poct_node == NULL) // case of Tree3D_spherical
                    continue;

                // $$$ NOTE: this block must be identical to block $$$ below
                is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                        num_pass, min_weight_sum_assoc, nominal_critical_node_size_km,
                        gap_primary_critical, gap_secondary_critical,
                        last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                        numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                        reference_phase_ttime_error,
                        channelParameters,
                        &globalBestValues,
                        // oct-tree node volume
                        &oct_node_volume,
                        time_min, time_max
                        );
                // END - this block must be identical to block $$$ below

                //
                /*
                if (is_global_best || nSamples % 1000 == 0) {
                    // get node size in km
                    double depth_corr = (AVG_ERAD - poct_node->center.z) / AVG_ERAD;
                    double dsx_global = poct_node->ds.x * DEG2KM * cos(DE2RA * poct_node->center.y) * depth_corr;
                    double dsy_global = poct_node->ds.y * DEG2KM * depth_corr;
                    printf("DEBUG: OctTree INITIAL: nsamples=%d/%d  is_global_best %d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  crit_size_km %g\n",
                            nSamples, max_num_nodes, is_global_best, poct_node->center.y, poct_node->center.x, poct_node->center.z,
                            dsy_global, dsx_global, poct_node->ds.z,
                            best_nassociated_P, best_nassociated, best_weight_sum, best_prob, best_ot_variance, nominal_critical_node_size_km);
                    fflush(stdout);
                }
                //*/
                nSamples++;

                // save indicator of solution probability
                if (0 && i_get_scatter_sample) {
                    //if (node_prob - min_prob_assoc_definative > 0.0) {
                    //double value = node_prob - min_prob_assoc_definative; // size
                    double value = 1.0; // size
                    //value = pow(value, min_weight_sum_assoc);
                    int index = *pn_scatter_sample * 4;
                    (*passoc_scatter_sample)[index++] = poct_node->center.x;
                    (*passoc_scatter_sample)[index++] = poct_node->center.y;
                    (*passoc_scatter_sample)[index++] = poct_node->center.z;
                    (*passoc_scatter_sample)[index++] = value;
                    (*pn_scatter_sample)++;
                    if (value > *p_global_max_nassociated_P_lat_lon)
                        *p_global_max_nassociated_P_lat_lon = value;
                    //}
                }
                if (!is_local_search) {
                    indexDepthSavedInitialGrid++;
                    if (USE_SAVED_INITIAL_GRIDS && indexDepthSavedInitialGrid > sizeSavedInitialGridTravelTimeDepth) {
                        printf("ERROR: indexDepthSavedInitialGrid > sizeSavedInitialGridTravelTimeDepth: this should not happen!\n");
                    }
                }

            }
            if (!is_local_search) {
                indexLatLonSavedInitialGrid++;
                if (USE_SAVED_INITIAL_GRIDS && indexLatLonSavedInitialGrid > sizeSavedInitialGridDistances) {
                    printf("ERROR: indexLatLonSavedInitialGrid > sizeSavedInitialGridDistances: this should not happen!\n");
                }
                if (USE_SAVED_INITIAL_GRIDS && indexLatLonSavedInitialGrid > sizeSavedInitialGridAzimuths) {
                    printf("ERROR: indexLatLonSavedInitialGrid > sizeSavedInitialGridAzimuths: this should not happen!\n");
                }
                if (USE_SAVED_INITIAL_GRIDS && indexDepthSavedInitialGrid > sizeSavedInitialGridTravelTimeLatLon) {
                    printf("ERROR: indexDepthSavedInitialGrid > sizeSavedInitialGridTravelTimeLatLon: this should not happen!\n");
                }
            }
            /*
            if (icount_sta) {
                (*pstation_count)[ilat][ilon][0] = max_nassociated_P_lat_lon;
                (*pstation_count)[ilat][ilon][1] = lat;
                (*pstation_count)[ilat][ilon][2] = lon;
                if (max_nassociated_P_lat_lon > *p_global_max_nassociated_P_lat_lon)
             *p_global_max_nassociated_P_lat_lon = max_nassociated_P_lat_lon;
            }
             */
        }
    }
    // no further use of initial regular grid distance and azimuth arrays
    indexLatLonSavedInitialGrid = -1;
    indexDepthSavedInitialGrid = -1;
    // DEBUG
    //printf("DEBUG: n_not_available %d, n_set, %d, n_read %d\n", DEBUG_n_not_available, DEBUG_n_set, DEBUG_n_read);


    printf("Info: Completed initial grid search (level %d; node_size_y=%gkm),   nSamples=%d/%d/%d, dt=%.3fs\n",
            0, smallest_node_size_divide_y, nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
    last_info_nsamples = nSamples;
    info_nsamples_start_time = clock();

    if (0 && i_get_scatter_sample) {
        int n_scatter_samples_initial = *pn_scatter_sample;
        // normalize initial scatter sample
        double value_max = *p_global_max_nassociated_P_lat_lon;
        if (value_max > FLT_MIN) {
            int index;
            for (index = 0; index < n_scatter_samples_initial; index++) {
                (*passoc_scatter_sample)[4 * index + 3] /= value_max * 10.0;
            }
        }
    }

    // check if minimum weight sum to associate reached
    if (globalBestValues.best_weight_sum < min_weight_sum_assoc) {
        printf("Info: Minimum weight sum associate not reached during initial search, rejecting location.\n");
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);
    }


    // Stage 2: loop over oct-tree nodes =================================================================

    //nScatterSaved = 0;
    //ipos = 0;

    // 20140228 AJL - Try to avoid patchy, poorly sampled location pdf's
#define SUBDIVIDE_NEIGHBORS 1
#ifdef SUBDIVIDE_NEIGHBORS
    int n_neigh;
    Vect3D coords;
    int nneighbor_tested = 0;
    int noutside = 0;
    int nsubdivided = 0;
#endif
    OctNode* neighbor_node;

    Vect3D smallest_node_size_best;
    smallest_node_size_best.x = FLT_MAX;
    smallest_node_size_best.y = FLT_MAX;
    smallest_node_size_best.z = FLT_MAX;

    double smallest_node_size_y_km = FLT_MAX;
    double current_min_node_size_km = nominal_min_node_size_km;
#ifndef PURE_OCTREE
    int current_level_divide = 0;
    printf("DEBUG: PURE_OCTREE NOT DEFINED! <=======\n");
#endif
    int highest_level_divided_w_assoc = 0;
    // 20150215 AJL - added oct_node_value_max[] for correct normalization of oct_tree_scatter_volume]
#define MAX_NUM_LEVEL_DIVIDES 999   // set very high so no checking needed
    static double oct_node_value_max[MAX_NUM_LEVEL_DIVIDES + 1];
    int n;
    for (n = 0; n < MAX_NUM_LEVEL_DIVIDES + 1; n++) {
        oct_node_value_max[n] = -1.0;
    }

    //int stop_on_min_node_size = 0;

    OctNode* poct_node_global_best = NULL;

    ResultTreeNode* presult_node;
    OctNode* pparent_oct_node;

    int num_divide_this_node_size = 0;

    *p_global_max_nassociated_P_lat_lon = 0;

#ifdef PURE_OCTREE
    //double min_node_size_km_assoc = (nominal_critical_node_size_km + nominal_min_node_size_km) / 2.0;
    double min_node_size_km_assoc = nominal_critical_node_size_km;
    // set maximum number levels to divide
    int max_num_level_assoc = 0;
    double test = lat_step * DEG2KM;
    while (test > min_node_size_km_assoc) {
        test /= 2.0;
        max_num_level_assoc++;
    }
    max_num_level_assoc++; // make sure max assoc level is larger than critical node size level
#else
    // set maximum number levels to divide
    int max_num_level_assoc = 0;
    double test = lat_step * DEG2KM;
    while (test > nominal_min_node_size_km) {
        test /= 2.0;
        max_num_level_assoc++;
    }
#endif

    // set critical node size level (must be reached in Oct-tree divides to accept association)
    int nominal_critical_node_size_level = 0;
    test = lat_step * DEG2KM;
#ifdef PURE_OCTREE
    nominal_critical_node_size_km /= 1.0;
#endif
    while (test > nominal_critical_node_size_km) {
        test /= 2.0;
        nominal_critical_node_size_level++;
    }
    double current_critical_node_size_km = nominal_critical_node_size_km;
    int critical_node_size_level = nominal_critical_node_size_level;
    //printf("DEBUG: critical_node_size_level %d\n", critical_node_size_level);

#ifdef PURE_OCTREE
#else
    int max_num_divide_single_node_level = (max_num_nodes - nSamples) / max_num_level_assoc;
    max_num_divide_single_node_level /= 2; // 20110203 AJL - to preserve available divides if critical node size reduced
    //max_num_divide_single_node_level = 999999;
    //printf("lat_step %g  nominal_min_node_size_km %g  max_num_level_assoc %d  (max_num_nodes - nSamples) %d  max_num_divide_single_node_level %d\n",
    //        lat_step * DEG2KM, nominal_min_node_size_km, max_num_level_assoc, (max_num_nodes - nSamples), max_num_divide_single_node_level);
    printf("Info: Maximum number divide each node level: %d\n", max_num_divide_single_node_level);
    //#define PROB_FACTOR 0.25
#define PROB_FACTOR_CONVERGENT 0.6  // 20111110 AJL
#define PROB_FACTOR_SAMPLE 0.4  // 20111110 AJL
    double prob_factor = PROB_FACTOR_CONVERGENT; // 20111110 AJL
#endif

#ifdef PURE_OCTREE
    globalBestValues.best_prob = -FLT_MAX;
#endif
    double best_prob_test_assoc = -FLT_MAX;

    while (nSamples < max_num_nodes) {
        // 20150907 AJL - TEST while (nSamples < (3 * max_num_nodes) / 4) {

        /*if (stop_on_min_node_size) {
            presult_node = getHighestLeafValue(pResultTreeRoot);
            if (presult_node == NULL) {
                printf("Info: No nodes with finite probability, terminating OctTree search.\n");
                break;
            }
        } else {
            presult_node = getHighestLeafValueMinSize(pResultTreeRoot, 0.0, min_node_size_deg, 0.0);
            // check if null node
            if (presult_node == NULL) {
                printf("Info: No nodes with finite probability or no more nodes larger than min_node_size, terminating OctTree search.\n");
                break;
            }
        }*/
        presult_node = NULL;
#ifdef PURE_OCTREE

        // check if can reduce critical node size
        //if (current_level_divide > critical_node_size_level // only allow change in critical_node_size at or after critical node level reached
        //        && best_prob >= min_prob_assoc) { // and good probability level
        //printf("DEBUG: distanceNodeSizeTest >? 0.0 %g && highest_level_divided_w_assoc %d >=? critical_node_size_level %d\n", distanceNodeSizeTest, highest_level_divided_w_assoc, critical_node_size_level);
        if (globalBestValues.distanceNodeSizeTest > 0.0
                && highest_level_divided_w_assoc >= critical_node_size_level) { // only allow change in critical_node_size at or after critical node level reached
            double new_critical_node_size_km = nominal_critical_node_size_km;
            double new_min_node_size_km = current_min_node_size_km;
            int new_critical_node_size_level = nominal_critical_node_size_level;
            int new_max_num_levels = max_num_level_assoc;
            // 20111005 AJL double ratioTest = globalBestValues.distanceNodeSizeTest / (double) numCountDistanceNodeSizeTest;
            // 20111005 AJL - new critical node size test algorithm
            double distance_node_size_test_scale = globalBestValues.distanceNodeSizeTest / 40.0;
            // compare current node size with 1/40th distance containing (2*min_weight_sum_assoc + number within hypo depth) P associated stations
            //    the latter is the distance within which there are 2*min_weight_sum_assoc associated stations at distance > hypo depth
            double ratioTest = current_critical_node_size_km / distance_node_size_test_scale;
            //printf("DEBUG: >>> Critical node size ratioTest: %g = current_critical_node_size_km %g / distance_node_size_test_scale %g\n", ratioTest, current_critical_node_size_km, distance_node_size_test_scale);
            //while (new_critical_node_size_level <= current_level_divide // only allow reduction to current level
            //        && ratioTest < 1.0) {
            while (ratioTest > 1.0 && new_critical_node_size_km / 2.0 >= min_critical_node_size_km) {
                ratioTest /= 2.0;
                new_critical_node_size_km /= 2.0;
                new_min_node_size_km /= 2.0;
                new_critical_node_size_level++;
                new_max_num_levels++;
            }
            if (new_critical_node_size_km < current_critical_node_size_km) { // only allow reduction in critical_node_size
                double last_critical_node_size_km = current_critical_node_size_km;
                current_critical_node_size_km = new_critical_node_size_km;
                double last_min_node_size_km = current_min_node_size_km;
                current_min_node_size_km = new_min_node_size_km;
                if (current_critical_node_size_km != last_critical_node_size_km) {
                    printf("Info: Critical node size changed: from %g to %g, level from %d to %d, min node size from %g to %g.   nCountClose/Test/All=%d/%d/%d  distClose/Test/All=%g/%g/%g  currSize/distTest=%g\n",
                            last_critical_node_size_km, current_critical_node_size_km, critical_node_size_level, new_critical_node_size_level,
                            last_min_node_size_km, current_min_node_size_km,
                            globalBestValues.nCountClose, globalBestValues.numCountDistanceNodeSizeTest, globalBestValues.nCountFar, globalBestValues.distanceClose,
                            globalBestValues.distanceNodeSizeTest, globalBestValues.distanceFar, current_critical_node_size_km / distance_node_size_test_scale);
                    fflush(stdout);
                }
                critical_node_size_level = new_critical_node_size_level;
                max_num_level_assoc = new_max_num_levels;
                //if (critical_node_size_level + 1 > max_num_level_assoc)
                //    max_num_level_assoc = critical_node_size_level + 1;
            }
        }

        // get best leaf node for subdivision
        presult_node = getHighestLeafValue(pResultTreeRoot);
        if (presult_node != NULL && presult_node->level >= max_num_level_assoc) {
            printf("Info: Oct-tree assoc: Min node size reached, terminating OctTree search (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    highest_level_divided_w_assoc + 1, smallest_node_size_y_km, nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
            last_info_nsamples = nSamples;
            info_nsamples_start_time = clock();
            break;
        }
        /*
        presult_node = getHighestLeafValueLESpecifiedLevel(pResultTreeRoot, critical_node_size_level + 1);
        // stop on first node at level critical_node_size_level + 1
        if (presult_node != NULL && presult_node->level > critical_node_size_level) {
            printf("Info: Terminating Oct-tree association: presult_node->level=%d > critical_node_size_level=%d\n", presult_node->level, critical_node_size_level);
            break;
        }*/
#else
        if (num_divide_this_node_size < max_num_divide_single_node_level) {
            presult_node = getHighestLeafValueAtSpecifiedLevel(pResultTreeRoot, current_level_divide);
            if (presult_node == NULL) {
                ;
                //printf("Info: No nodes with finite probability or no more nodes dividing level %d.\n", current_level_divide);
                //if (presult_node != NULL && best_prob > 0.0 && presult_node->value < best_prob * prob_factor)
                //} else if (best_prob > 0.0 && presult_node->value < best_prob - min_prob_assoc) {
                //    printf("Info: No nodes with prob > best_prob-min_prob_assoc=%g.\n", best_prob - min_prob_assoc);
            } else if (globalBestValues.best_prob > 0.0 && presult_node->value < globalBestValues.best_prob * prob_factor) { // 20110125 AJL - need to check against relatively low prob to avoid early exclusion of better solution cells
                //printf("Info: No nodes with prob > best_prob * prob_factor=%g.\n", best_prob * prob_factor);
                presult_node = NULL;
            }
            // 20150326 AJL  if (presult_node != NULL)
            // 20150326 AJL      highest_level_divided_w_assoc = current_level_divide;
        } else {
            ;
            //printf("Info: At max number nodes for single node size.\n");
        }
        // check if null node
        if (presult_node == NULL) {
            printf("DEBUG: critical_node_size_level %d\n", critical_node_size_level);
            // no more nodes of specified size
            printf("Info: No more nodes available dividing level %d (size=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    current_level_divide, smallest_node_size_divide_y, nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
            last_info_nsamples = nSamples;
            info_nsamples_start_time = clock();
            // check if prob to associate reached
            if (current_level_divide <= critical_node_size_level && globalBestValues.best_weight_sum < min_weight_sum_assoc) {
                printf("Info: Minimum weight sum associate not reached during divide at current level %d, rejecting location.\n", current_level_divide);
                // free oct-tree structures - IMPORTANT!
                freeOctTreeStructures(pResultTreeRoot, pOctTree);
                return (0.0);
            }
            // check if at critical node size level
            if (current_level_divide == critical_node_size_level) { // 20111110 AJL
                prob_factor = PROB_FACTOR_SAMPLE;
                //min_prob_assoc = min_prob_assoc_definative; // enables prob checking in oct tree core and below
            }
            current_level_divide++;
            // check if minimum node size reached
            //if (stop_on_min_node_size && smallest_node_size.y < current_min_node_size_km) {
            if (current_level_divide >= max_num_level_assoc) {
                printf("Info: Min node size reached, terminating OctTree search at level %d (size=%gkm): min_node_size_y_km=%g\n",
                        current_level_divide, smallest_node_size_y_km, current_min_node_size_km);
                break;
            }

            // check if can reduce critical node size
            //if (current_level_divide > critical_node_size_level // only allow change in critical_node_size at or after critical node level reached
            //        && best_prob >= min_prob_assoc) { // and good probability level
            if (globalBestValues.distanceNodeSizeTest > 0.0
                    && current_level_divide > critical_node_size_level) { // only allow change in critical_node_size at or after critical node level reached
                double new_critical_node_size_km = nominal_critical_node_size_km;
                double new_min_node_size_km = nominal_min_node_size_km;
                int new_critical_node_size_level = nominal_critical_node_size_level;
                int new_max_num_levels = max_num_level_assoc;
                // 20111005 AJL double ratioTest = globalBestValues.distanceNodeSizeTest / (double) numCountDistanceNodeSizeTest;
                // 20111005 AJL - new critical node size test algorithm
                double distance_node_size_test_scale = globalBestValues.distanceNodeSizeTest / 40.0;
                // compare current node size with 1/40th distance containing (2*min_weight_sum_assoc + number within hypo depth) P associated stations
                //    the latter is the distance within which there are 2*min_weight_sum_assoc associated stations at distance > hypo depth
                double ratioTest = current_critical_node_size_km / distance_node_size_test_scale;
                //printf("DEBUG: >>> Critical node size ratioTest: %g = current_critical_node_size_km %g / distance_node_size_test_scale %g\n", ratioTest, current_critical_node_size_km, distance_node_size_test_scale);
                //while (new_critical_node_size_level <= current_level_divide // only allow reduction to current level
                //        && ratioTest < 1.0) {
                while (ratioTest > 1.0 && new_critical_node_size_km / 2.0 >= min_critical_node_size_km) {
                    ratioTest /= 2.0;
                    new_critical_node_size_km /= 2.0;
                    new_min_node_size_km /= 2.0;
                    new_critical_node_size_level++;
                    new_max_num_levels++;
                }
                if (new_critical_node_size_km < current_critical_node_size_km) { // only allow reduction in critical_node_size
                    double last_critical_node_size_km = current_critical_node_size_km;
                    current_critical_node_size_km = new_critical_node_size_km;
                    double last_min_node_size_km = current_min_node_size_km;
                    current_min_node_size_km = new_min_node_size_km;
                    if (current_critical_node_size_km != last_critical_node_size_km) {
                        printf("Info: Critical node size changed: from %g to %g, min node size from %g to %g, level from %d to %d.   nCountClose/Test/All=%d/%d/%d  distClose/Test/All=%g/%g/%g  currSize/distTest=%g\n",
                                last_critical_node_size_km, current_critical_node_size_km,
                                last_min_node_size_km, current_min_node_size_km, critical_node_size_level, new_critical_node_size_level,
                                globalBestValues.nCountClose, globalBestValues.numCountDistanceNodeSizeTest, globalBestValues.nCountFar, globalBestValues.distanceClose, globalBestValues.distanceNodeSizeTest, globalBestValues.distanceFar, current_critical_node_size_km / distance_node_size_test_scale);
                        fflush(stdout);
                    }
                    critical_node_size_level = new_critical_node_size_level;
                    max_num_level_assoc = new_max_num_levels;
                    //if (critical_node_size_level + 1 > max_num_level_assoc)
                    //    max_num_level_assoc = critical_node_size_level + 1;
                }
            }

            presult_node = getHighestLeafValueAtSpecifiedLevel(pResultTreeRoot, current_level_divide);

            if (presult_node != NULL) {
                // 20140627 AJL - bug fix?  smallest_node_size_divide_y = poct_node->ds.y * DEG2KM;
                smallest_node_size_divide_y = presult_node->pnode->ds.y * DEG2KM;
                //printf("Info: Node at level %d (size=%gkm) found.   nSamples=%d/%d/%d, dt=%.3fs\n",
                //        current_level_divide, smallest_node_size_divide_y, nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
                //last_info_nsamples = nSamples;
                //info_nsamples_start_time = clock();
                //if (current_level_divide <= critical_node_size_level + 1) // 20101226
                globalBestValues.best_prob = -FLT_MAX;
                num_divide_this_node_size = 0;
            }
        }
#endif

        // check if null node
        if (presult_node == NULL) {
#ifdef PURE_OCTREE
            printf("Info: Oct-tree assoc: No nodes with finite probability or no more nodes dividing (smallest_node_size_divide_y=%g)\n", smallest_node_size_divide_y);
#else
            printf("Info: No nodes with finite probability or no more nodes dividing level %d (smallest_node_size_divide_y=%g)\n", current_level_divide, smallest_node_size_divide_y);
#endif
            break;
        }

        pparent_oct_node = presult_node->pnode;

#ifdef PURE_OCTREE
        smallest_node_size_divide_y = presult_node->pnode->ds.y * DEG2KM;
        if (presult_node->level > highest_level_divided_w_assoc) {
            highest_level_divided_w_assoc = presult_node->level;
        }
#else
        highest_level_divided_w_assoc = current_level_divide; // 20150326 AJL - moved here, before may have missed last divide level
#endif

        // subdivide node and evaluate solution at each child
#ifndef SUBDIVIDE_NEIGHBORS
        subdivide(pparent_oct_node, OCTTREE_UNDEF_VALUE, NULL);
        neighbor_node = pparent_oct_node;
#endif

#ifdef SUBDIVIDE_NEIGHBORS

        // subdivide all HighestLeafValue neighbors

        double fraction_of_node_size = pparent_oct_node->ds.x / 10.0;
        for (n_neigh = 0; n_neigh < 7; n_neigh++) {

            if (n_neigh == 0) {
                neighbor_node = pparent_oct_node;
            } else {
                coords.x = pparent_oct_node->center.x;
                coords.y = pparent_oct_node->center.y;
                coords.z = pparent_oct_node->center.z;
                if (n_neigh == 1) {
                    coords.x = pparent_oct_node->center.x
                            + (pparent_oct_node->ds.x + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 2) {
                    coords.x = pparent_oct_node->center.x
                            - (pparent_oct_node->ds.x + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 3) {
                    coords.y = pparent_oct_node->center.y
                            + (pparent_oct_node->ds.y + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 4) {
                    coords.y = pparent_oct_node->center.y
                            - (pparent_oct_node->ds.y + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 5) {
                    coords.z = pparent_oct_node->center.z
                            + (pparent_oct_node->ds.z + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 6) {
                    coords.z = pparent_oct_node->center.z
                            - (pparent_oct_node->ds.z + fraction_of_node_size) / 2.0;
                }

                // find neighbor node
                // 20160912 AJL - Bug fix: added check before and after removing wrap-around, to support regional or reduced volume searches that straddle +/-180deg
                //printf("getLeafNodeContaining A: parent: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) n_neigh %d %.3f,%.3f,%.3f", pparent_oct_node->center.x, pparent_oct_node->center.y, pparent_oct_node->center.z, pparent_oct_node->ds.x, pparent_oct_node->ds.y, pparent_oct_node->ds.z, n_neigh, coords.x, coords.y, coords.z);
                neighbor_node = getLeafNodeContaining(pOctTree, coords);
                if (neighbor_node == NULL) {
                    // check for strike wrap-around
                    if (coords.x > 180.0)
                        coords.x -= 360.0;
                    else if (coords.x < -180.0)
                        coords.x += 360.0;
                    neighbor_node = getLeafNodeContaining(pOctTree, coords);
                }
                //printf(" -> %.3f,%.3f,%.3f", coords.x, coords.y, coords.z);

                nneighbor_tested++;
                // outside of octTree volume
                if (neighbor_node == NULL) {
                    //printf(" -> outside\n");
                    noutside++;
                    continue;
                }
                // already subdivided
                if (neighbor_node->ds.z < 0.99 * pparent_oct_node->ds.z) {
                    //printf(" -> subdivided\n");
                    nsubdivided++;
                    continue;
                }
                //printf(" -> OK\n");
            }


            // subdivide node and evaluate solution at each child
            subdivide(neighbor_node, OCTTREE_UNDEF_VALUE, NULL);

#endif

            int ix, iy, iz;
            for (ix = 0; ix < 2; ix++) {
                for (iy = 0; iy < 2; iy++) {
                    for (iz = 0; iz < 2; iz++) {

                        num_divide_this_node_size++;

                        poct_node = neighbor_node->child[ix][iy][iz];

                        //if (poct_node->ds.x < min_node_size || poct_node->ds.y < min_node_size || poct_node->ds.z < min_node_size)
                        //fprintf(stderr, "\nnode size too small!! %lf %lf %lf\n", poct_node->ds.x, poct_node->ds.y, poct_node->ds.z);

                        // get node size in km
                        double depth_corr = (AVG_ERAD - poct_node->center.z) / AVG_ERAD;
                        double dsx_global = poct_node->ds.x * DEG2KM * cos(DE2RA * poct_node->center.y) * depth_corr;
                        double dsy_global_surface = poct_node->ds.y * DEG2KM;
                        double dsy_global = dsy_global_surface * depth_corr;
                        if (dsy_global_surface < smallest_node_size_y_km)
                            smallest_node_size_y_km = dsy_global_surface;

                        // $$$ NOTE: this block must be identical to block $$$ above
                        is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                                num_pass, min_weight_sum_assoc, current_critical_node_size_km,
                                gap_primary_critical, gap_secondary_critical,
                                last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                                reference_phase_ttime_error,
                                channelParameters,
                                &globalBestValues,
                                // oct-tree node volume
                                &oct_node_volume,
                                time_min, time_max
                                );
                        // END - this block must be identical to block $$$ above
                        if (globalBestValues.best_prob > oct_node_value_max[poct_node->level]) { // 20150215 AJL - added oct_node_value_max[]
                            oct_node_value_max[poct_node->level] = globalBestValues.best_prob;
                        }
                        nSamples++;


                        // if best node so far
#ifdef PURE_OCTREE
                        if (is_global_best && presult_node->level >= critical_node_size_level) {
                            best_prob_test_assoc = globalBestValues.best_prob * globalBestValues.effective_cell_size;
                            //printf("DEBUG: best_prob_test_assoc [%f] = globalBestValues.best_prob [%f] * globalBestValues.effective_cell_size [%f]\n",
                            //        best_prob_test_assoc, globalBestValues.best_prob, globalBestValues.effective_cell_size);
                            //VOLUMEbest_prob_test_assoc = best_prob;
#else
                        if (is_global_best) {
                            best_prob_test_assoc = globalBestValues.best_prob;
#endif
                            // save node and node size
                            //num_divide_this_node_size = 0;
                            found_global_best_during_divide = 1;
                            poct_node_global_best = poct_node;
                            // check for smallest of global best
                            if (dsx_global < smallest_node_size_best.x)
                                smallest_node_size_best.x = dsx_global;
                            if (dsy_global < smallest_node_size_best.y)
                                smallest_node_size_best.y = dsy_global;
                            if (poct_node->ds.z < smallest_node_size_best.z)
                                smallest_node_size_best.z = poct_node->ds.z;

                        }

                        //
                        /*
                        if (is_global_best || nSamples % 1000 == 0) {
                            printf("DEBUG: OctTree SUBDIVIDE: nsamples=%d/%d  is_global_best %d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  best_prob_test_assoc %g  crit_size_km %g\n",
                                    nSamples, max_num_nodes, is_global_best, poct_node->center.y, poct_node->center.x, poct_node->center.z,
                                    dsy_global, dsx_global, poct_node->ds.z,
                                    best_nassociated_P, best_nassociated, best_weight_sum, best_prob, best_prob_test_assoc, best_ot_variance, current_critical_node_size_km);
                            fflush(stdout);
                        }
                        //*/

                    } // end triple loop over node children
                }
            }

#ifdef SUBDIVIDE_NEIGHBORS
        } // end loop over HighestLeafValue neighbors
#endif

        if (nSamples >= max_num_nodes) {
            // 20150907 AJL - TEST if (nSamples >= (3 * max_num_nodes) / 4) {
#ifdef PURE_OCTREE
            printf("Info: Oct-tree assoc: Reached maximum number nodes (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    highest_level_divided_w_assoc, smallest_node_size_divide_y, nSamples, nSamples + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
#else
            printf("Info: Oct-tree assoc: Reached maximum number nodes dividing level %d (size=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    highest_level_divided_w_assoc, smallest_node_size_divide_y, nSamples, nSamples + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
#endif
            last_info_nsamples += nSamples;
            info_nsamples_start_time = clock();
        }
    } // end while (nSamples < max_num_nodes)


    // check if global best not found during divide
    if (!found_global_best_during_divide) {
        printf("Info: Oct-tree assoc: Global best not found during divide, rejecting location.   nSamples=%d/%d/%d, dt=%.3fs\n",
                nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);
    }
#ifndef PURE_OCTREE
    // check if critical oct-tree cell size reached
    if (smallest_node_size_best.y > nominal_critical_node_size_km) {
        printf("Info: Critical node size dy not reached for best solution, rejecting location.  smallest best node: dy=%g, [dx=%g, dz=%g], critical_cell_size=%gkm.   nSamples=%d/%d/%d, dt=%.3fs\n",
                smallest_node_size_best.y, smallest_node_size_best.x, smallest_node_size_best.z, nominal_critical_node_size_km,
                nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        last_info_nsamples = nSamples;
        info_nsamples_start_time = clock();
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);
    }
#endif

#ifdef PURE_OCTREE
    // check if probability to associate reached
    if (best_prob_test_assoc < min_prob_assoc) {
        printf("Info: Oct-tree assoc: Minimum probability associate not reached not reached during divide, rejecting location.   nSamples=%d/%d/%d, dt=%.3fs\n",
                nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);
    }
#else
    // check if probability to associate reached
    if (best_prob_test_assoc < min_prob_assoc) {
        printf("Info: Minimum probability associate not reached not reached during divide, rejecting location.   nSamples=%d/%d/%d, dt=%.3fs\n",
                nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        last_info_nsamples = nSamples;
        info_nsamples_start_time = clock();
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);
    }
#endif


    // re-evaluate global best solution to set association information
    if (poct_node_global_best != NULL) {
        poct_node = poct_node_global_best;
        last_reassociate = 1;
        //best_weight_sum = FLT_MIN;
        //best_prob = -FLT_MAX;

        // run octtree core to set deDataTmp values for best solution
        // $$$ NOTE: this block must be identical to block $$$ above
        is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                num_pass, min_weight_sum_assoc, current_critical_node_size_km,
                gap_primary_critical, gap_secondary_critical,
                last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                reference_phase_ttime_error,
                channelParameters,
                &globalBestValues,
                // oct-tree node volume
                &oct_node_volume,
                time_min, time_max
                );
        // END - this block must be identical to block $$$ above
        nSamples++;

        //printf("DEBUG: OctTree BEST provisional: nsamples=%d/%d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  crit_size_km %g\n",
        //        nSamples, max_num_nodes, poct_node->center.y, poct_node->center.x, poct_node->center.z,
        //        dsy_global, dsx_global, poct_node->ds.z,
        //        best_nassociated_P, best_nassociated, best_weight_sum, best_prob, best_ot_variance, current_critical_node_size_km);
        //fflush(stdout);
    }

    // check for and un-associate outlier associated stations (e.g. isolated, very distant stations associated with small, local event.)
    // 20121126 AJL - added
    if (1) {
        globalBestValues.best_nassociated -= unassociateOutlierStations(pdata_list, num_de_data, num_pass, globalBestValues.best_depth, globalBestValues.best_ot_mean);
    }

    // set data association information (provisional))
    // needed for location statistics / scatter sample ?
    int definitive = 0;
    int reAssociateOnly = 0;
    double best_depth_step = -1.0; // not used
    setDataAssociationInformation(hypocenter, pdata_list, num_de_data, num_pass, definitive, reAssociateOnly, best_depth_step, poct_node, current_critical_node_size_km);

    //*/


    /*
    // set mean vector of epicenter to station vectors
    hypocenter->vector_azimuth = atan2(x_vector_sum, y_vector_sum) * RA2DE;
    if (hypocenter->vector_azimuth < 0.0)
        hypocenter->vector_azimuth += 360.0;
    hypocenter->vector_distance = sqrt(x_vector_sum * x_vector_sum + y_vector_sum * y_vector_sum);

    // set azimuth gap
    hypocenter->gap_primary = calcAzimuthGap(AzimuthArray, num_az_gap, &(hypocenter->gap_secondary));

    //hypocenter->ot = (time_t) ((long) best_ot_mean);
    hypocenter->otime = best_ot_mean;
    hypocenter->ot_std_dev = best_ot_variance < 0.0 ? 0.0 : sqrt(best_ot_variance);
    hypocenter->lat = best_lat;
    if (best_lon < -180.0)
        best_lon += 360.0;
    else if (best_lon > 180.0)
        best_lon -= 360.0;
    hypocenter->lon = best_lon;
    hypocenter->depth = best_depth;
    hypocenter->depth_step = smallest_node_size_best.z;
    hypocenter->nassoc_P = best_nassociated_P;
    hypocenter->nassoc = best_nassociated;
    hypocenter->dist_min = gcd_min;
    hypocenter->dist_max = gcd_max;
     */



    // Stage 3: loop over oct-tree nodes without re-association =================================================================
    // get extra samples at next divide level past highest_level_divided_w_assoc without re-association to get more precise hypocenter and to use for location statistics

    last_reassociate = 0;
    no_reassociate = 1;
    globalBestValues.best_prob = -FLT_MAX;
    int nSamplesExtra = 0;
    //last_info_nsamples = nSamplesExtra;

#ifdef PURE_OCTREE
    // set maximum number levels to divide
    int max_num_level_loc = 0;
    test = lat_step * DEG2KM;
    while (test > current_min_node_size_km) {
        test /= 2.0;
        max_num_level_loc++;
    }
    int max_num_divide_extra = (max_num_nodes - nSamples);
    int highest_level_divided_w_loc = -1;
    //printf("DEBUG: Oct-tree loc: (level %d; node_size_y=%gkm), nSamplesExtra %d, max_num_divide_extra %d, max_num_nodes %d, nSamples %d\n",
    //        highest_level_divided_w_loc, smallest_node_size_y_km, nSamplesExtra, max_num_divide_extra, max_num_nodes, nSamples);
#else
    //int max_num_divide_extra = MAX_NUM_DIVIDE_EXTRA;
    int max_num_divide_extra = 8 * max_num_divide_single_node_level; // 20111114 AJL
    //int max_num_divide_extra = 16 * max_num_divide_single_node_level; // 20150221 AJL
    current_level_divide = highest_level_divided_w_assoc + 1;
    printf("DEBUG: Oct-tree loc: (level %d; node_size_y=%gkm), nSamplesExtra %d, max_num_divide_extra %d, max_num_nodes %d, nSamples %d\n",
            current_level_divide, smallest_node_size_y_km, nSamplesExtra, max_num_divide_extra, max_num_nodes, nSamples);
#endif

    while (nSamplesExtra < max_num_divide_extra) {

        presult_node = NULL;

#ifdef PURE_OCTREE
        presult_node = getHighestLeafValueGESpecifiedLevel(pResultTreeRoot, highest_level_divided_w_assoc + 1);
        //printf("DEBUG: getHighestLeafValue (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n", presult_node->level, presult_node->volume, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        // check if null node
        if (presult_node == NULL) {
            printf("Info: Oct-tree loc: : No nodes with finite probability or no more nodes (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    highest_level_divided_w_loc, smallest_node_size_y_km, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
            last_info_nsamples += nSamplesExtra;
            info_nsamples_start_time = clock();
            break;
        }
        if (presult_node->level > max_num_level_loc) {
            printf("Info: Oct-tree loc: Min node size reached, terminating OctTree search (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    highest_level_divided_w_loc + 1, smallest_node_size_y_km, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
            last_info_nsamples += nSamplesExtra;
            info_nsamples_start_time = clock();
            break;
        }
#else
        presult_node = getHighestLeafValueAtSpecifiedLevel(pResultTreeRoot, current_level_divide);
        // check if null node
        if (presult_node == NULL) {
            printf("Info: Oct-tree loc: No nodes with finite probability or no more nodes dividing level %d (size=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                    current_level_divide, smallest_node_size_divide_y, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
            last_info_nsamples += nSamplesExtra;
            info_nsamples_start_time = clock();
            break;
        }
#endif
        //else if (nSamplesExtra == 0) {
        // 20140627 AJL - bug fix?  smallest_node_size_divide_y = poct_node->ds.y * DEG2KM;
        smallest_node_size_divide_y = presult_node->pnode->ds.y * DEG2KM;
        //    printf("Info: Oct-tree loc: Node at level %d (size=%gkm) found.\n",
        //            current_level_divide, smallest_node_size_divide_y);
        //}

        pparent_oct_node = presult_node->pnode;

#ifdef PURE_OCTREE
        if (presult_node->level > highest_level_divided_w_loc) {
            highest_level_divided_w_loc = presult_node->level;
        }
#endif

        // subdivide node and evaluate solution at each child
#ifndef SUBDIVIDE_NEIGHBORS
        subdivide(pparent_oct_node, OCTTREE_UNDEF_VALUE, NULL);
        neighbor_node = pparent_oct_node;
#endif

#ifdef SUBDIVIDE_NEIGHBORS

        // subdivide all HighestLeafValue neighbors

        double fraction_of_node_size = pparent_oct_node->ds.x / 10.0;
        for (n_neigh = 0; n_neigh < 7; n_neigh++) {

            if (n_neigh == 0) {
                neighbor_node = pparent_oct_node;
            } else {
                coords.x = pparent_oct_node->center.x;
                coords.y = pparent_oct_node->center.y;
                coords.z = pparent_oct_node->center.z;
                if (n_neigh == 1) {
                    coords.x = pparent_oct_node->center.x
                            + (pparent_oct_node->ds.x + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 2) {
                    coords.x = pparent_oct_node->center.x
                            - (pparent_oct_node->ds.x + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 3) {
                    coords.y = pparent_oct_node->center.y
                            + (pparent_oct_node->ds.y + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 4) {
                    coords.y = pparent_oct_node->center.y
                            - (pparent_oct_node->ds.y + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 5) {
                    coords.z = pparent_oct_node->center.z
                            + (pparent_oct_node->ds.z + fraction_of_node_size) / 2.0;
                } else if (n_neigh == 6) {
                    coords.z = pparent_oct_node->center.z
                            - (pparent_oct_node->ds.z + fraction_of_node_size) / 2.0;
                }

                // find neighbor node
                // 20160912 AJL - Bug fix: added check before and after removing wrap-around, to support regional or reduced volume searches that straddle +/-180deg
                //printf("getLeafNodeContaining A: parent: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) n_neigh %d %.3f,%.3f,%.3f", pparent_oct_node->center.x, pparent_oct_node->center.y, pparent_oct_node->center.z, pparent_oct_node->ds.x, pparent_oct_node->ds.y, pparent_oct_node->ds.z, n_neigh, coords.x, coords.y, coords.z);
                neighbor_node = getLeafNodeContaining(pOctTree, coords);
                if (neighbor_node == NULL) {
                    // check for strike wrap-around
                    if (coords.x > 180.0)
                        coords.x -= 360.0;
                    else if (coords.x < -180.0)
                        coords.x += 360.0;
                    neighbor_node = getLeafNodeContaining(pOctTree, coords);
                }
                //printf(" -> %.3f,%.3f,%.3f", coords.x, coords.y, coords.z);

                nneighbor_tested++;
                // outside of octTree volume
                if (neighbor_node == NULL) {
                    //printf(" -> outside\n");
                    noutside++;
                    continue;
                }
                // already subdivided
                if (neighbor_node->ds.z < 0.99 * pparent_oct_node->ds.z) {
                    //printf(" -> subdivided\n");
                    nsubdivided++;
                    continue;
                }
                //printf(" -> OK\n");
            }


            // subdivide node and evaluate solution at each child
            subdivide(neighbor_node, OCTTREE_UNDEF_VALUE, NULL);

#endif


            int ix, iy, iz;
            for (ix = 0; ix < 2; ix++) {
                for (iy = 0; iy < 2; iy++) {
                    for (iz = 0; iz < 2; iz++) {

                        //num_divide_this_node_size++;

                        poct_node = neighbor_node->child[ix][iy][iz];

                        //if (poct_node->ds.x < min_node_size || poct_node->ds.y < min_node_size || poct_node->ds.z < min_node_size)
                        //fprintf(stderr, "\nnode size too small!! %lf %lf %lf\n", poct_node->ds.x, poct_node->ds.y, poct_node->ds.z);

                        // get node size in km
                        double depth_corr = (AVG_ERAD - poct_node->center.z) / AVG_ERAD;
                        double dsx_global = poct_node->ds.x * DEG2KM * cos(DE2RA * poct_node->center.y) * depth_corr;
                        double dsy_global_surface = poct_node->ds.y * DEG2KM;
                        double dsy_global = dsy_global_surface * depth_corr;
                        if (dsy_global_surface < smallest_node_size_y_km)
                            smallest_node_size_y_km = dsy_global_surface;

                        // $$$ NOTE: this block must be identical to block $$$ above
                        is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                                num_pass, min_weight_sum_assoc, current_critical_node_size_km,
                                gap_primary_critical, gap_secondary_critical,
                                last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                                reference_phase_ttime_error,
                                channelParameters,
                                &globalBestValues,
                                // oct-tree node volume
                                &oct_node_volume,
                                time_min, time_max
                                );
                        // END - this block must be identical to block $$$ above
                        if (globalBestValues.best_prob > oct_node_value_max[poct_node->level]) { // 20150215 AJL - added oct_node_value_max[]
                            oct_node_value_max[poct_node->level] = globalBestValues.best_prob;
                        }
                        nSamplesExtra++;

                        // if best node so far
                        if (is_global_best) {

                            // save node and node size
                            //num_divide_this_node_size = 0;
                            //found_global_best_during_divide = 1;
                            poct_node_global_best = poct_node;
                            // check for smallest of global best
                            if (dsx_global < smallest_node_size_best.x)
                                smallest_node_size_best.x = dsx_global;
                            if (dsy_global < smallest_node_size_best.y)
                                smallest_node_size_best.y = dsy_global;
                            if (poct_node->ds.z < smallest_node_size_best.z)
                                smallest_node_size_best.z = poct_node->ds.z;

                        }

                        if (0 && (is_global_best || nSamplesExtra % 1000 == 0)) {
                            printf("DEBUG: OctTree SUBDIVIDE Extra: nSamplesExtra=%d/%d  level %d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  crit_size_km %g\n",
                                    nSamplesExtra, max_num_divide_extra, poct_node->level, poct_node->center.y, poct_node->center.x, poct_node->center.z,
                                    dsy_global, dsx_global, poct_node->ds.z,
                                    globalBestValues.best_nassociated_P, globalBestValues.best_nassociated, globalBestValues.best_weight_sum, globalBestValues.best_prob, globalBestValues.best_ot_variance, current_critical_node_size_km);
                            fflush(stdout);
                        }

                    } // end triple loop over node children
                }
            }
#ifdef SUBDIVIDE_NEIGHBORS
        } // end loop over HighestLeafValue neighbors
#endif

    } // end while (nSamplesExtra < max_num_divide_extra)

    if (nSamplesExtra >= max_num_divide_extra) {
#ifdef PURE_OCTREE
        printf("Info: Oct-tree loc: Reached maximum number nodes (level %d; node_size_y=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                highest_level_divided_w_loc, smallest_node_size_divide_y, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
#else
        printf("Info: Oct-tree loc: Reached maximum number nodes dividing level %d (size=%gkm).   nSamples=%d/%d/%d, dt=%.3fs\n",
                current_level_divide, smallest_node_size_divide_y, nSamplesExtra, nSamplesExtra + last_info_nsamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
#endif
        last_info_nsamples += nSamplesExtra;
        info_nsamples_start_time = clock();
    }
    nSamples += nSamplesExtra;


    // re-evaluate global best solution to get secondary phases
    if (poct_node_global_best != NULL) {
        poct_node = poct_node_global_best;
        last_reassociate = 1;
        no_reassociate = 0; // if = 1 will use last reassociate phases only, and will only calculate solution quality and update ResultTree
        //best_weight_sum = FLT_MIN;
        //best_prob = -FLT_MAX;

        // run octtree core to set deDataTmp values for best solution
        // $$$ NOTE: this block must be identical to block $$$ above
        is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                num_pass, min_weight_sum_assoc, current_critical_node_size_km,
                gap_primary_critical, gap_secondary_critical,
                last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                reference_phase_ttime_error,
                channelParameters,
                &globalBestValues,
                // oct-tree node volume
                &oct_node_volume,
                time_min, time_max
                );
        // END - this block must be identical to block $$$ above
        nSamples++;

        // set data association information (provisional))
        int definitive = 0;
        int reAssociateOnly = 0;
        double best_depth_step = -1.0; // not used
        setDataAssociationInformation(hypocenter, pdata_list, num_de_data, num_pass, definitive, reAssociateOnly, best_depth_step, poct_node, current_critical_node_size_km);

        // run octtree core again to check for association of unassociated, definitive phase with increased ttime error to avoid outliers
        try_assoc_remaining_definitive = 1;
        // $$$ NOTE: this block must be identical to block $$$ above
        is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
                num_pass, min_weight_sum_assoc, current_critical_node_size_km,
                gap_primary_critical, gap_secondary_critical,
                last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                reference_phase_ttime_error,
                channelParameters,
                &globalBestValues,
                // oct-tree node volume
                &oct_node_volume,
                time_min, time_max
                );
        // END - this block must be identical to block $$$ above
        if (globalBestValues.best_prob > oct_node_value_max[poct_node->level]) { // 20150215 AJL - added oct_node_value_max[]
            oct_node_value_max[poct_node->level] = globalBestValues.best_prob;
        }
        nSamples++;
        try_assoc_remaining_definitive = 0;
        //printf("DEBUG: OctTree BEST provisional: nsamples=%d/%d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  crit_size_km %g\n",
        //        nSamples, max_num_nodes, poct_node->center.y, poct_node->center.x, poct_node->center.z,
        //        dsy_global, dsx_global, poct_node->ds.z,
        //        best_nassociated_P, best_nassociated, best_weight_sum, best_prob, best_ot_variance, current_critical_node_size_km);
        //fflush(stdout);

        // check for and un-associate outlier associated stations (e.g. isolated, very distant stations associated with small, local event.)
        // 20121126 AJL - added
        if (1) {
            globalBestValues.best_nassociated -= unassociateOutlierStations(pdata_list, num_de_data, num_pass, globalBestValues.best_depth, globalBestValues.best_ot_mean);
        }

        // set data association information (provisional))
        // needed for location statistics / scatter sample ?
        definitive = 1;
        reAssociateOnly = 0;
        best_depth_step = -1.0; // not used
        setDataAssociationInformation(hypocenter, pdata_list, num_de_data, num_pass, definitive, reAssociateOnly, best_depth_step, poct_node, current_critical_node_size_km);

    } else {
        printf("Error: Oct-tree assoc: Null global best node, rejecting location: this should not happen!   nSamples=%d/%d/%d, dt=%.3fs\n",
                nSamples - last_info_nsamples, nSamples, max_num_nodes, (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        // free oct-tree structures - IMPORTANT!
        freeOctTreeStructures(pResultTreeRoot, pOctTree);
        return (0.0);

    }


    // Stage 4: get location statistics =================================================================
    // calculate "traditional" statistics for lon/lat/depth/value data in scatter_sample, excluding initial sample (n_scatter_samples_initial)

    if (i_get_scatter_sample) {

        // re-seed random number generator so that results will be identical to previous if no data changed.
        int RandomNumSeed = 9589;
        SRAND_FUNC(RandomNumSeed);

#ifdef COMPILE_OLD_VERSION
        int node_value_type = VALUE_IS_PROBABILITY_IN_NODE; // used before 20150218
#else
#ifdef PURE_OCTREE
        //int node_value_type = VALUE_IS_PROB_DENSITY_IN_NODE;
        int node_value_type = VALUE_IS_PROBABILITY_IN_NODE;
#else
        int node_value_type = VALUE_IS_PROB_DENSITY_IN_NODE; // used after 20150218
#endif
#endif
        Vect3D expect = expect_NULL;
        Mtrx3D cov = cov_NULL;
        Ellipsoid3D ellipsoid = ellipsoid_NULL;
        Ellipse2D ellipse = ellipse_NULL;
        int npoints = 0;
        double oct_tree_scatter_volume = 0.0;
        double oct_node_value_ref = 1.0; // cell probabilities are divided by this value, enables normalization to a specified max value

        if (hypocenter->scatter_sample != NULL) {
            // 20150417 AJL
#ifdef PURE_OCTREE
            //int level_stats_min = nominal_critical_node_size_level + 2;
            int level_stats_min = highest_level_divided_w_assoc + 2;
            //int level_stats_min = critical_node_size_level - 1; // 20150907 AJL - first level whose dist is not directly influenced by assoc dist
            //int level_stats_min = 2 + (nominal_critical_node_size_level + highest_level_divided_w_assoc) / 2;
            int level_stats_max = 999;
#else
            int level_stats_min = highest_level_divided_w_assoc + 2;
            int level_stats_max = 999;
#endif
            // determine integral of all oct-tree leaf node pdf values
            oct_node_value_ref = -FLT_MAX;
            for (int n = level_stats_min; n <= level_stats_max; n++) { // 20150417 AJL
                oct_node_value_ref = fmax(oct_node_value_ref, oct_node_value_max[n]);
            }
            double oct_tree_integral = integrateResultTreeAtLevels(pResultTreeRoot, node_value_type, 0.0, oct_node_value_ref, level_stats_min, level_stats_max);
            //printf("DEBUG: level_stats_min %d  level_stats_max %d  oct_tree_integral %g \n", level_stats_min, level_stats_max, oct_tree_integral);
            //oct_node_value_ref = 1.0; // TEST!
            // get scatter sample for location statistics
            int num_scatter = NUM_LOC_SCATTER_SAMPLE;
            npoints = octtree_GenEventScatter(pResultTreeRoot, node_value_type, num_scatter, oct_node_value_ref,
                    hypocenter->scatter_sample, oct_tree_integral, &oct_tree_scatter_volume, level_stats_min, level_stats_max, "location");
            /* 20150218 AJL
            if (npoints > 64) // min num samples to accept, otherwise decrease level
                break;
            //printf("DEBUG: decrease level: level_stats_min %d  level_stats_max %d\n", level_stats_min, level_stats_max);
            level_stats_max--;
            level_stats_min--;
             */
            // 20150218 AJL  }
            if (npoints < 16 || level_stats_min == 0) {
                printf("ERROR: Too few (%d) scatter samples available at level>=%d!\n", npoints, level_stats_min);
                hypocenter->nscatter_sample = 0;
            } else {
                // calculate location statistics
                hypocenter->nscatter_sample = npoints;
                expect = CalcExpectationSamplesGlobal(hypocenter->scatter_sample, npoints, globalBestValues.best_lon);
                cov = CalcCovarianceSamplesGlobal(hypocenter->scatter_sample, npoints, &expect);
                ellipsoid = CalcErrorEllipsoid(&cov, DELTA_CHI_SQR_68_3);
                ellipse = CalcHorizontalErrorEllipse(&cov, DELTA_CHI_SQR_68_2);
            }
        }
        // get association scatter sample for plotting
#ifdef PURE_OCTREE
        int level_min = critical_node_size_level - 1;
#else
        int level_min = critical_node_size_level;
#endif
        int level_max = highest_level_divided_w_assoc + 1;
        oct_node_value_ref = -FLT_MAX;
        for (int n = level_min; n <= level_max; n++) { // 20150417 AJL
            oct_node_value_ref = fmax(oct_node_value_ref, oct_node_value_max[n]);
        }
        double oct_tree_integral = integrateResultTreeAtLevels(pResultTreeRoot, node_value_type, 0.0, oct_node_value_ref, level_min, level_max);
        //printf("DEBUG: level_min %d  level_max %d  oct_tree_integral %g \n", level_min, level_max, oct_tree_integral);
        //oct_node_value_ref = 1.0; // TEST!
        npoints = octtree_GenEventScatter(pResultTreeRoot, node_value_type, *p_n_alloc_assoc_scatter_sample, oct_node_value_ref,
                *passoc_scatter_sample, oct_tree_integral, &oct_tree_scatter_volume, level_min, level_max, "association");
        *pn_scatter_sample = npoints;
        if (0) { // DEBUG calc and output
            Vect3D expectAssoc = CalcExpectationSamplesGlobal(*passoc_scatter_sample, npoints, globalBestValues.best_lon);
            Mtrx3D covAssoc = CalcCovarianceSamplesGlobal(*passoc_scatter_sample, npoints, &expectAssoc);
            Ellipsoid3D ellipsoidAssoc = CalcErrorEllipsoid(&covAssoc, DELTA_CHI_SQR_68_3);
            fprintf(stdout, "DEBUG: Assoc: EllVol  %le\n", ellipsoidAssoc.len3 * ellipsoidAssoc.len2 * ellipsoidAssoc.len1);
        }
        printf("Info: Finished generating event scatter samples.   dt=%.3fs\n", (double) (clock() - info_nsamples_start_time) / CLOCKS_PER_SEC);
        info_nsamples_start_time = clock();


        // NLL STATISTICS
        fprintf(stdout,
                "Info: Hypo Statistics:");
        fprintf(stdout,
                "  ExpectX %lg Y %lg Z %lg",
                expect.x, expect.y, expect.z);
        fprintf(stdout, "  CovXX %lg XY %lg XZ %lg",
                cov.xx, cov.xy, cov.xz);
        fprintf(stdout, " YY %lg YZ %lg",
                cov.yy, cov.yz);
        fprintf(stdout, " ZZ %lg",
                cov.zz);
        fprintf(stdout, " EllAz1  %lg Dip1  %lg Len1  %lg",
                ellipsoid.az1, ellipsoid.dip1,
                ellipsoid.len1);
        fprintf(stdout, " Az2  %lg Dip2  %lg Len2  %lg",
                ellipsoid.az2, ellipsoid.dip2,
                ellipsoid.len2);
        fprintf(stdout, " Len3  %le\n", ellipsoid.len3);
        fprintf(stdout, " EllVol  %le\n", ellipsoid.len3 * ellipsoid.len2 * ellipsoid.len1);

        /*  20110617 AJL - Bug fix?  TODO: the constant should be DELTA_CHI_SQR_68_3 for 3-D ???
        double errx = sqrt(DELTA_CHI_SQR_68_2 * cov.xx);
        double erry = sqrt(DELTA_CHI_SQR_68_2 * cov.yy);
        double errz = sqrt(DELTA_CHI_SQR_68_2 * cov.zz);
         */
        double errx = sqrt(DELTA_CHI_SQR_68_3 * cov.xx);
        double erry = sqrt(DELTA_CHI_SQR_68_3 * cov.yy);
        double errz = sqrt(DELTA_CHI_SQR_68_3 * cov.zz);

        double azMaxHorUnc = ellipse.az1 + 90.0;
        if (azMaxHorUnc >= 360.0)
            azMaxHorUnc -= 360.0;
        if (azMaxHorUnc >= 180.0)
            azMaxHorUnc -= 180.0;
        fprintf(stdout,
                "  OriginUncertainty (68%%)  latUnc %lg  lonUnc %lg  minHorUnc %lg  maxHorUnc %lg  azMaxHorUnc %lg  vertUnc %g\n",
                erry, errx, ellipse.len1, ellipse.len2, azMaxHorUnc, errz
                );

        hypocenter->errh = ellipse.len2;
        hypocenter->errz = errz;
        hypocenter->expect = expect;
        hypocenter->cov = cov;
        hypocenter->ellipsoid = ellipsoid;
        hypocenter->ellipse = ellipse;

    }

    // clean up
    clear_OtimeLimitList(&NumOtimeLimit);
    //free(AzimuthArray);

    // free oct-tree structures - IMPORTANT!
    freeOctTreeStructures(pResultTreeRoot, pOctTree);

    return (globalBestValues.best_weight_sum);

}

double octtreeGlobalAssociationLocation_reassociateOnly(int num_pass, double min_weight_sum_assoc,
        double gap_primary_critical, double gap_secondary_critical,
        int numPhaseTypesUse, int phaseTypesUse[get_num_ttime_phases()], char channelNamesUse[get_num_ttime_phases()][128], double timeDelayUse[MAX_NUM_TTIME_PHASES][2],
        double reference_phase_ttime_error,
        TimedomainProcessingData** pdata_list, int num_de_data,
        HypocenterDesc *hypocenter,
        ChannelParameters * channelParameters,
        time_t time_min, time_t time_max
        ) {

    OctNode* poct_node_global_best = &(hypocenter->global_best_oct_node);
    if (poct_node_global_best == NULL) {
        return (-1);
    }


    // need at least 2 stations for location
    if (num_de_data < 2) {
        return (-1);
    }

    int n_ot_limit, j;
    TimedomainProcessingData* deData;

    // allocate OtimeLimit pool array
    int size = 2 * num_de_data * get_num_ttime_phases();
    if (sizeOtimeLimitPool < size) {
        OtimeLimitPool = (OtimeLimit*) realloc(OtimeLimitPool, size * sizeof (OtimeLimit));
        if (OtimeLimitPool == NULL)
            fprintf(stderr, "ERROR: realloc of OtimeLimitPool: %s\n", strerror(errno));
        sizeOtimeLimitPool = size;
    }

    // allocate deDataTmp pool array
    TimedomainProcessingData* deDataTmp = NULL;
    if (sizeDataTmpPool < num_de_data) {
        DataTmpPool = (TimedomainProcessingData*) realloc(DataTmpPool, num_de_data * sizeof (TimedomainProcessingData));
        if (DataTmpPool == NULL)
            fprintf(stderr, "ERROR: realloc of DataTmpPool: %s\n", strerror(errno));
        sizeDataTmpPool = num_de_data;
        // initialize data
        for (j = 0; j < num_de_data; j++) {
            deDataTmp = DataTmpPool + j;
            clear_deData_assoc(deDataTmp);
        }
    }

    // initialize global best solution variables
    globalBestValues.best_ot_variance = -1.0;
    globalBestValues.best_ot_mean = 0;
    globalBestValues.best_lat = 0.0;
    globalBestValues.best_lon = 0.0;
    globalBestValues.best_depth = 0.0;
    globalBestValues.node_ot_variance_weight = 0;
    globalBestValues.best_nassociated_P = 0;
    globalBestValues.best_nassociated = 0;
    globalBestValues.best_weight_sum = FLT_MIN; // non-zero to prevent divide by zero
    globalBestValues.best_prob = -FLT_MAX;
    globalBestValues.effective_cell_size = -1.0;
    globalBestValues.nCountClose = 0;
    globalBestValues.numCountDistanceNodeSizeTest = 0;
    globalBestValues.nCountFar = 0;
    globalBestValues.distanceClose = 0.0;
    globalBestValues.distanceNodeSizeTest = 0.0;
    globalBestValues.distanceFar = 0.0;


    // set min cutoff below which distance weight = 1.0
#ifdef COMPILE_OLD_VERSION
    // 20150225 AJL - remove this algorithm, it changes weighting as function of available picks, always varying and may be somewhat arbitrary
    // use mean of maximum distance between each station and all other stations, helps to stabilize location when only a few nearby station available
    double distance_weight_dist_min_current = getMeanDistanceMax(pdata_list, num_de_data);
    // check that this distance does not exceed DISTANCE_WEIGHT_DIST_MIN
    printf("=====================> 1 DEBUG: distance_weight_dist_min_current %f\n", distance_weight_dist_min_current);
    if (distance_weight_dist_min_current > DISTANCE_WEIGHT_DIST_MIN)
        distance_weight_dist_min_current = DISTANCE_WEIGHT_DIST_MIN;
    printf("=====================> 2 DEBUG: distance_weight_dist_min_current %f\n", distance_weight_dist_min_current);
#else
    // 20150225 AJL - set constant dist min, set to 30 deg
    double distance_weight_dist_min_current = DISTANCE_WEIGHT_DIST_MIN;
#endif


    // clear data association and perform checks if first pass
    if (num_pass == 1) {
        for (n_ot_limit = 0; n_ot_limit < num_de_data; n_ot_limit++) {
            // clear data association fields
            deData = pdata_list[n_ot_limit];
            clear_deData_assoc(deData);
            // if flagged as use_for_location but will be skipped, check if twin stream pick can  be used for location
            // 20121130 AJL - added
            if (deData->use_for_location && skipData(deData, -999)) {
                TimedomainProcessingData* deDataTwin = deData->use_for_location_twin_data;
                if (deDataTwin != NULL && !skipData(deDataTwin, -999)) {
                    deDataTwin->use_for_location = 1;
                    deData->use_for_location = 0;
                    //printf("DEBUG: >>>> set twin to use for loc: num_pass %d: skipped: %s %s %s %d:%d:%f %s %s -> use: %s %s %s %d:%d:%f %s %s\n",
                    //        num_pass, deData->network, deData->station, deData->channel, deData->hour, deData->min, deData->dsec, pick_stream_name(deData), deData->phase,
                    //        deDataTwin->network, deDataTwin->station, deDataTwin->channel, deData->hour, deData->min, deData->dsec, pick_stream_name(deDataTwin), deDataTwin->phase
                    //        );
                }
            }
            // amplitude attenuation
            deData->amplitude_error_ratio = -1.0;
            deData->sta_corr = 0.0;
        }
    }

    //*p_global_max_nassociated_P_lat_lon = 0.0;

    ResultTreeNode* pResultTreeRoot = NULL; // Octtree likelihood*volume results tree root node

    int nSamples = 0;
    int is_global_best = 0;

    double min_prob_assoc = -FLT_MAX;
    double min_prob_assoc_definative = min_weight_sum_assoc - 1.0;
    //min_prob_assoc = 2.0; // 20101226
    min_prob_assoc = min_prob_assoc_definative; // 20101226B

    globalBestValues.node_prob = 0.0;
    globalBestValues.node_ot_variance = 0.0;
    double oct_node_volume;

    int last_reassociate = 1; // if = 1, will check all phases for association, not just P etc. counted in location.  should use only at optimal hypocenter...
    double last_reassociate_otime = hypocenter->otime; // set to a value > 0.0 to only uses phases that give predicted otime range that includes hypocenter->otime
    double last_reassociate_ot_std_dev = hypocenter->ot_std_dev; // must be >= 0.0 if last_reassociate_otime > 0.0
    int try_assoc_remaining_definitive = 0; // if = 1, will check definitive phases (e.g. P, PKP) for association using increased ttime error
    int no_reassociate = 0; // if = 1 will use last reassociate phases only, and will only calculate solution quality and update ResultTree

    // re-evaluate global best solution to get secondary phases
    OctNode* poct_node = poct_node_global_best;
    // no use of initial regular grid distance and azimuth arrays
    int indexLatLonSavedInitialGrid = -1;
    int indexDepthSavedInitialGrid = -1;
    double current_critical_node_size_km = hypocenter->global_best_critical_node_size_km;
    //best_weight_sum = FLT_MIN;
    //best_prob = -FLT_MAX;

    // run octtree core to set deDataTmp values for best solution
    // $$$ NOTE: this block must be identical to block $$$ above
    is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
            num_pass, min_weight_sum_assoc, current_critical_node_size_km,
            gap_primary_critical, gap_secondary_critical,
            last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
            numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
            reference_phase_ttime_error,
            channelParameters,
            &globalBestValues,
            // oct-tree node volume
            &oct_node_volume,
            time_min, time_max
            );
    // END - this block must be identical to block $$$ above
    nSamples++;

    // set data association information (provisional))
    int definitive = 0;
    int reAssociateOnly = 1;
    double best_depth_step = -1.0; // not used
    setDataAssociationInformation(hypocenter, pdata_list, num_de_data, num_pass, definitive, reAssociateOnly, best_depth_step, poct_node, current_critical_node_size_km);

    // run octtree core again to check for association of unassociated, definitive phase with increased ttime error to avoid outliers
    try_assoc_remaining_definitive = 1;
    // $$$ NOTE: this block must be identical to block $$$ above
    is_global_best = octtree_core(&pResultTreeRoot, poct_node, indexLatLonSavedInitialGrid, indexDepthSavedInitialGrid, pdata_list, num_de_data,
            num_pass, min_weight_sum_assoc, current_critical_node_size_km,
            gap_primary_critical, gap_secondary_critical,
            last_reassociate, last_reassociate_otime, last_reassociate_ot_std_dev, try_assoc_remaining_definitive, no_reassociate, distance_weight_dist_min_current,
            numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
            reference_phase_ttime_error,
            channelParameters,
            &globalBestValues,
            // oct-tree node volume
            &oct_node_volume,
            time_min, time_max
            );
    // END - this block must be identical to block $$$ above
    nSamples++;
    try_assoc_remaining_definitive = 0;
    //printf("DEBUG: OctTree BEST provisional: nsamples=%d/%d  lat %g lon %g z %g  dlat/lon/z %g/%g/%g  nassociated_P %d  nassociated %d  weight_sum %g  prob %g  ot_var %g  crit_size_km %g\n",
    //        nSamples, max_num_nodes, poct_node->center.y, poct_node->center.x, poct_node->center.z,
    //        dsy_global, dsx_global, poct_node->ds.z,
    //        best_nassociated_P, best_nassociated, best_weight_sum, best_prob, best_ot_variance, current_critical_node_size_km);
    //fflush(stdout);

    // check for and un-associate outlier associated stations (e.g. isolated, very distant stations associated with small, local event.)
    // 20121126 AJL - added
    if (1) {
        globalBestValues.best_nassociated -= unassociateOutlierStations(pdata_list, num_de_data, num_pass, globalBestValues.best_depth, globalBestValues.best_ot_mean);
    }

    // set data association information (provisional))
    // needed for location statistics / scatter sample ?
    definitive = 1;
    reAssociateOnly = 1;
    best_depth_step = -1.0; // not used
    setDataAssociationInformation(hypocenter, pdata_list, num_de_data, num_pass, definitive, reAssociateOnly, best_depth_step, poct_node, current_critical_node_size_km);




    // clean up
    clear_OtimeLimitList(&NumOtimeLimit);
    //free(AzimuthArray);

    // free oct-tree structures - IMPORTANT!
    Tree3D *pOctTree = NULL;
    freeOctTreeStructures(pResultTreeRoot, pOctTree);

    return (globalBestValues.best_weight_sum);

}

double octtreeGlobalAssociationLocation(int num_pass, double min_weight_sum_assoc, int max_num_nodes,
        double nominal_critical_node_size_km, double min_critical_node_size_km, double nominal_min_node_size_km,
        double gap_primary_critical, double gap_secondary_critical,
        double lat_min, double lat_max, double lat_step, double lon_min, double lon_max, double lon_step_smallest,
        double depth_min, double depth_max, double depth_step, int is_local_search,
        int numPhaseTypesUse, int phaseTypesUse[get_num_ttime_phases()], char channelNamesUse[get_num_ttime_phases()][128], double timeDelayUse[MAX_NUM_TTIME_PHASES][2],
        double reference_phase_ttime_error,
        TimedomainProcessingData** pdata_list, int num_de_data,
        HypocenterDesc *hypocenter,
        float **passoc_scatter_sample, int *p_n_alloc_assoc_scatter_sample, int i_get_scatter_sample, int *pn_scatter_sample, double *p_global_max_nassociated_P_lat_lon,
        ChannelParameters * channelParameters,
        int reassociate_only, time_t time_min, time_t time_max
        ) {


    if (reassociate_only) {
        return (
                octtreeGlobalAssociationLocation_reassociateOnly(num_pass, min_weight_sum_assoc,
                gap_primary_critical, gap_secondary_critical,
                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                reference_phase_ttime_error,
                pdata_list, num_de_data,
                hypocenter,
                channelParameters, time_min, time_max)
                );
    } else {
        return (
                octtreeGlobalAssociationLocation_full(num_pass, min_weight_sum_assoc, max_num_nodes,
                nominal_critical_node_size_km, min_critical_node_size_km, nominal_min_node_size_km,
                gap_primary_critical, gap_secondary_critical,
                lat_min, lat_max, lat_step, lon_min, lon_max, lon_step_smallest,
                depth_min, depth_max, depth_step, is_local_search,
                numPhaseTypesUse, phaseTypesUse, channelNamesUse, timeDelayUse,
                reference_phase_ttime_error,
                pdata_list, num_de_data,
                hypocenter,
                passoc_scatter_sample, p_n_alloc_assoc_scatter_sample, i_get_scatter_sample, pn_scatter_sample, p_global_max_nassociated_P_lat_lon,
                channelParameters, time_min, time_max)
                );
    }
}

void freeChannelParameters(ChannelParameters * pchan_params) {

    free(pchan_params->sta_corr);

}

void initChannelParameters(ChannelParameters * pchan_params) {

    strcpy(pchan_params->network, "?");
    strcpy(pchan_params->station, "?");
    strcpy(pchan_params->location, "?");
    strcpy(pchan_params->channel, "?");
    pchan_params->lat = -999.9;
    pchan_params->lon = -999.9;
    pchan_params->elev = 0.0;
    pchan_params->inactive_duplicate = 0;
    pchan_params->staActiveInReportInterval = 0;
    pchan_params->data_latency = 0.0;
    pchan_params->last_data_end_time = -999.0;
    pchan_params->numData = 0;
    pchan_params->numDataAssoc = 0;
    pchan_params->qualityWeight = 1.0;
    pchan_params->geogfile_checked_time = -1; // set earlier than any possible file mod time to force first check
    pchan_params->internet_station_query_checked_time = -1;
    pchan_params->have_coords = 0;
    pchan_params->error = 0;
    pchan_params->count_non_contiguous = 0;
    pchan_params->level_non_contiguous = 0.0;
    pchan_params->count_conflicting_dt = 0;
    pchan_params->level_conflicting_dt = 0.0;

    // 20160601 AJL - sta corr changed from P only to any phase
    pchan_params->sta_corr = calloc(get_num_ttime_phases(), sizeof (StaCorrections));
    for (int n = 0; n < get_num_ttime_phases(); n++) {
        pchan_params->sta_corr[n].valid = 0;
    }
    pchan_params->sta_corr_checked = 0;
    pchan_params->process_this_channel_orientation = -1;
    pchan_params->channel_set[0] = pchan_params->channel_set[1] = -1;
}

/** find and associate source_id for other channel orientations matching this net/sta/chan
 *
 *  returns number of matching channel found and associated
 */

int associate3CompChannelSet(ChannelParameters* channel_params, int n_sources, int source_id) {

    char *network = channel_params[source_id].network;
    char *station = channel_params[source_id].station;
    char *location = channel_params[source_id].location;
    char *channel = channel_params[source_id].channel;

    int first_found = -1;
    int second_found = -1;
    int nfound = 0;

    for (int n = 0; n < n_sources; n++) {
        if (n != source_id) {
            if (strcmp(network, channel_params[n].network) == 0
                    && strcmp(station, channel_params[n].station) == 0
                    && strcmp(location, channel_params[n].location) == 0
                    && strncmp(channel, channel_params[n].channel, 2) == 0) {
                if (first_found < 0) {
                    first_found = n;
                    channel_params[first_found].channel_set[0] = source_id;
                    channel_params[source_id].channel_set[0] = first_found;
                    channel_params[first_found].channel_set[1] = -1;
                    channel_params[source_id].channel_set[1] = -1;
                    nfound++;
                    printf("DEBUG: associate3CompChannelSet %d %s_%s_%s_%s: first_found: %d %s_%s_%s_%s\n",
                            source_id, network, station, location, channel,
                            n, channel_params[n].network, channel_params[n].station, channel_params[n].location, channel_params[n].channel);
                } else {
                    second_found = n;
                    channel_params[first_found].channel_set[1] = second_found;
                    channel_params[source_id].channel_set[1] = second_found;
                    channel_params[second_found].channel_set[0] = first_found;
                    channel_params[second_found].channel_set[1] = source_id;
                    nfound++;
                    printf("DEBUG: associate3CompChannelSet %d %s_%s_%s_%s: second_found: %d %s_%s_%s_%s\n",
                            source_id, network, station, location, channel,
                            n, channel_params[n].network, channel_params[n].station, channel_params[n].location, channel_params[n].channel);
                    break;
                }
            }
        }
    }

    printf("DEBUG: associate3CompChannelSet %d %s_%s_%s_%s: nfound: %d\n",
            source_id, network, station, location, channel, nfound);
    return (nfound);

}

void initAssociateLocateParameters(
        double upweight_picks_sn_cutoff_init,
        double upweight_picks_dist_max_init, double upweight_picks_dist_min_init,
        int use_amplitude_attenuation_init
        ) {

    upweight_picks_sn_cutoff = upweight_picks_sn_cutoff_init;
    upweight_picks_dist_max = upweight_picks_dist_max_init;
    upweight_picks_dist_full = upweight_picks_dist_min_init;
    use_amplitude_attenuation = use_amplitude_attenuation_init;

}


/** add a ChannelParameters to a sorted ChannelParameters list
 * list will be sorted by network and then station
 */

#define SIZE_INCREMENT_STA_PARAMS 64

int addChannelParametersToSortedList(ChannelParameters* pnew_chan_params, ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params) {

    if (*psorted_chan_params_list == NULL) { // list not yet created
        *psorted_chan_params_list = calloc(SIZE_INCREMENT_STA_PARAMS, sizeof (ChannelParameters*));
        // 20130930 AJL - bug fix
        //} else if (*pnum_sorted_chan_params != 0 && (*pnum_sorted_chan_params % SIZE_INCREMENT_STA_PARAMS) == 0) { // list will be too small
    } else if (*pnum_sorted_chan_params + 1 > (int) (sizeof (*psorted_chan_params_list) / sizeof (ChannelParameters*))) { // list will be too small
        ChannelParameters** new_sorted_list = NULL;
        new_sorted_list = calloc(*pnum_sorted_chan_params + SIZE_INCREMENT_STA_PARAMS, sizeof (ChannelParameters*));
        int n;
        for (n = 0; n < *pnum_sorted_chan_params; n++)
            new_sorted_list[n] = (*psorted_chan_params_list)[n];
        free(*psorted_chan_params_list);
        *psorted_chan_params_list = new_sorted_list;
    }

    // find sort position
    int nsort = -1;
    int ninsert;
    for (ninsert = 0; ninsert < *pnum_sorted_chan_params; ninsert++) {
        if ((nsort = strcmp((*psorted_chan_params_list)[ninsert]->network, pnew_chan_params->network)) > 0)
            break;
        if (nsort == 0) {
            if ((nsort = strcmp((*psorted_chan_params_list)[ninsert]->station, pnew_chan_params->station)) > 0)
                break;
            if (nsort == 0) {
                if ((nsort = strcmp((*psorted_chan_params_list)[ninsert]->location, pnew_chan_params->location)) > 0)
                    break;
                if (nsort == 0) {
                    if ((nsort = strcmp((*psorted_chan_params_list)[ninsert]->channel, pnew_chan_params->channel)) >= 0)
                        break;
                }
            }
        }
    }

    if (nsort != 0) { // not identical match, will insert    // 20140114 AJL - bug fix
        // shift higher sort sta params
        if (ninsert < *pnum_sorted_chan_params) {
            int m;
            for (m = *pnum_sorted_chan_params - 1; m >= ninsert; m--)
                (*psorted_chan_params_list)[m + 1] = (*psorted_chan_params_list)[m];
        }
    }
    // insert reference to ChannelParameters
    (*psorted_chan_params_list)[ninsert] = pnew_chan_params;
    if (nsort != 0) { // not identical match, will insert    // 20140114 AJL - bug fix
        (

                *pnum_sorted_chan_params)++;
    }

    return (0);

}

/** remove a ChannelParameters from a sorted ChannelParameters list */

void removeChannelParametersFromSortedList(ChannelParameters* pchan_params, ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params) {

    int i = 0;

    while (i < *pnum_sorted_chan_params) {
        if ((*psorted_chan_params_list)[i] == pchan_params)
            break;
        i++;
    }

    if (i == *pnum_sorted_chan_params) // not found
        return;

    while (i < *pnum_sorted_chan_params - 1) {
        (

                *psorted_chan_params_list)[i] = (*psorted_chan_params_list)[i + 1];
        i++;
    }

    (*psorted_chan_params_list)[*pnum_sorted_chan_params - 1] = NULL;

    (*pnum_sorted_chan_params)--;

}

/** clean up sorted ChannelParameters list memory */

void free_ChannelParametersList(ChannelParameters*** psorted_chan_params_list, int* pnum_sorted_chan_params) {

    if (*psorted_chan_params_list == NULL || *pnum_sorted_chan_params < 1)
        return;

    free(*psorted_chan_params_list);
    *psorted_chan_params_list = NULL;

    *pnum_sorted_chan_params = 0;

}

