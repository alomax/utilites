// phase arrivals in the ak135 model

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#include "../statistics/statistics.h"
#include "../picker/PickData.h"
#include "timedomain_processing_data.h"
#include "ttimes.h"

// conditional compiling of different velocity models / travel-time grids / phase types, etc.
#ifdef TTIMES_AK135_REGIONAL
#include "ttimes/ttimes_ak135_regional.h"
#include "ttimes/ttimes_ak135_regional_times_phases.h"
#include "ttimes/ttimes_ak135_regional_toang_phases.h"
#include "ttimes/tvel_ak135.h"
#else
#ifdef TTIMES_SOCALSN_HK_2010
#include "ttimes/ttimes_SoCalSN_HK_2010.h"
#include "ttimes/ttimes_SoCalSN_HK_2010_times_phases.h"
#include "ttimes/ttimes_SoCalSN_HK_2010_toang_phases.h"
#include "ttimes/tvel_scsn_hk.h"
#else
#ifdef TTIMES_MARSITE1_REGIONAL
#include "ttimes/ttimes_marsite1_regional.h"
#include "ttimes/ttimes_marsite1_regional_times_phases.h"
#include "ttimes/ttimes_marsite1_regional_toang_phases.h"
//#include "ttimes/ttimes_marsite1_regional_tvel.h"     // 20140320 AJL - marsite1 regional tvel produces large errors, rms, etc.  Effect of very low vel in shallow crsut?
#include "ttimes/tvel_ak135.h"
#else
#ifdef TTIMES_MARSITE2_REGIONAL
#include "ttimes/ttimes_marsite2_regional.h"
#include "ttimes/ttimes_marsite2_regional_times_phases.h"
#include "ttimes/ttimes_marsite2_regional_toang_phases.h"
#include "ttimes/tvel_ak135.h"
#else
#ifdef TTIMES_AK135_GLOBAL_TAUP
#include "ttimes/ttimes_ak135_0-800_10_TauP.h"
#include "ttimes/ttimes_ak135_0-800_10_TauP_times_phases.h"
#include "ttimes/ttimes_ak135_0-800_10_TauP_toang_phases.h"
#include "ttimes/tvel_ak135.h"
#else
// default
#ifndef TTIMES_AK135_GLOBAL
#define TTIMES_AK135_GLOBAL
#endif
#endif
#endif
#endif
#endif
#endif
// default
#ifdef TTIMES_AK135_GLOBAL
#include "ttimes/ttimes_ak135_0-800_10.h"
#include "ttimes/ttimes_ak135_0-800_10_times_phases.h"
#include "ttimes/ttimes_ak135_0-800_10_toang_phases.h"
#include "ttimes/tvel_ak135.h"
#include "ttimes/latlon_deep.h"
// 20160509 AJL - added to help avoid deep locations where not possible or likely:
#define LAT_LON_DEEP_AVAILABLE
#define DIST_TIME_DEPTH_LATLON_MIN_DEPTH 100.0
#endif


#include "depth_corr_mwpd.h"
#include "table_Q_Delta_PV.h"
#include "table_Q_Depth_Delta_PV.h"
#include "../geometry/geometry.h"
#include "../ran1/ran1.h"
#include "../octtree/octtree.h"

/** returns 1 if phase_id is P phase
 */

int is_P(int phase_id) {

    return (phase_id == P_PHASE_INDEX);

}

/** returns 1 if phase_id is S phase
 */

int is_S(int phase_id) {

    return (phase_id == S_PHASE_INDEX);

}

/** returns 1 if phase_id is a P phase leaving the source
 */

int is_P_at_source(int phase_id) {

    return (phase_type_flags[phase_id] & P_AT_SOURCE);

}

/** returns 1 if phase_id is a P phase arriving at the station
 */

int is_P_at_station(int phase_id) {

    return (phase_type_flags[phase_id] & P_AT_STATION);

}

/** returns 1 if phase_id is a direct P phase
 */

int is_direct_P(int phase_id) {

    return (phase_type_flags[phase_id] & DIRECT_P);

}

/** returns 1 if phase_id is a first arriving P phase
 */

int is_first_arrival_P(int phase_id) {

    return (phase_type_flags[phase_id] & FIRST_ARRIVAL_P);

}

/** returns 1 if phase_id is a S phase leaving the source
 */

int is_S_at_source(int phase_id) {

    return (phase_type_flags[phase_id] & S_AT_SOURCE);

}

/** returns 1 if phase_id is a S phase arriving at the station
 */

int is_S_at_station(int phase_id) {

    return (phase_type_flags[phase_id] & S_AT_STATION);

}

/** returns 1 if phase_id is a direct S phase
 */

int is_direct_S(int phase_id) {

    return (phase_type_flags[phase_id] & DIRECT_S);

}

/** returns 1 if phase_id is a first arriving S phase
 */

int is_first_arrival_S(int phase_id) {

    return (phase_type_flags[phase_id] & FIRST_ARRIVAL_S);

}

/** returns 1 if phase_id can count in location
 */

int is_count_in_location(int phase_id) {

    return (phase_type_flags[phase_id] & COUNT_IN_LOCATION);

}

/** returns 1 if should check to prevent overlap of origin times corresponding to different phases for each pick,
 *      preventing multiple phase association for a pick
 */

int doCheckOverlapCountInLocation() {

    return (CHECK_OVERLAP_COUNT_IN_LOCATION);

}

/** returns number of travel-time phase types defined for this model and time grid
 */

char *get_model_ttime_name() {

    return (MODEL_TTIME_NAME);

}

/** returns number of travel-time phase types defined for this model and time grid
 */

int get_num_ttime_phases() {

    return (NUM_TTIME_PHASES);

}

/** returns maximum epicentral distance defined for this model and time grid
 */

double get_dist_time_dist_max() {

    return (DIST_TIME_DIST_MAX);

}

/** returns maximum depth defined for this model and time grid
 */

double get_dist_time_depth_max() {

    return (DIST_TIME_DEPTH_MAX);

}

/** returns depth weight for this model and time grid
 *
 * weight is linear from 1.0 at DIST_TIME_DEPTH_MAX_FULL_WEIGHT to 0.0 at DIST_TIME_DEPTH_MAX
 *
 */
// 20160509 AJL - added to help avoid deep locations where not possible or likely

/*double get_dist_time_depth_weight(double depth) {

#ifdef DIST_TIME_DEPTH_MAX_FULL_WEIGHT
    if (depth <= DIST_TIME_DEPTH_MAX_FULL_WEIGHT) {
        return (1.0);
    } else if (depth > DIST_TIME_DEPTH_MAX) {
        return (0.0);
    }
    double weight = (DIST_TIME_DEPTH_MAX - depth) / (DIST_TIME_DEPTH_MAX - DIST_TIME_DEPTH_MAX_FULL_WEIGHT);
    return (weight);
#else
    return (1.0);
#endif

}*/

/** returns depth weight for this model and time grid
 *
 * weight = 1.0 if depth <= DIST_TIME_DEPTH_LATLON_MIN_DEPTH or previous seismicty depth > 100km
 *   otherwise weight = 0.5
 *
 */
// 20160509 AJL - added to help avoid deep locations where not possible or likely

double get_dist_time_latlon_depth_weight(double depth, double lat, double lon) {

#ifdef LAT_LON_DEEP_AVAILABLE
    if (depth <= DIST_TIME_DEPTH_LATLON_MIN_DEPTH) {
        return (1.0);
    }
    int ilat = (int) (lat + 90.0); // 0-179
    if (ilat > 179) {
        ilat = 179;
    } else {
        if (ilat < 0) {
            ilat = 0;
        }
    }
    int ilon = (int) (lon + 180.0); // 0-359
    if (ilon > 359) {
        ilon -= 360;
    } else {
        if (ilon < 0) {
            ilon += 360;
        }
    }
    int ideep = lat_lon_deep[ilat][ilon];
    //printf("lat %f, lon %f, depth %f, ilat %d, ilon %d, ideep %d\n", lat, lon, depth, ilat, ilon, ideep);
    if (ideep) {
        return (1.0);
    } else {
        return (0.5);
    }
#else
    return (1.0);
#endif

}

/** returns phase type index for P
 */

int get_P_phase_index() {

    return (P_PHASE_INDEX);

}

/** returns phase type index for P
 */

int get_S_phase_index() {

    return (S_PHASE_INDEX);

}

/** velocity
 * returns P or S velocity or density (D) from tvel table
 *
 * NOTE: Assumes model has constant velocity layers!  Gives (small?) error if gradient layers!
 */

double get_velocity_model_property(char property, double distance, double depth) {

    int index = NUM_TVEL_DEPTH - 1;
    while (index >= 0 && depth_Vp_Vs_rho[index][0] > depth)
        index--;

    if (index < 0) // check if depth index valid
        return (-1.0);

    if (property == 'P')
        return (depth_Vp_Vs_rho[index][1]);
    else if (property == 'S')
        return (depth_Vp_Vs_rho[index][2]);
    else if (property == 'D')
        return (depth_Vp_Vs_rho[index][3]);

    return (-1.0);

}

/** velocity
 * returns P or S velocity or density (D) from tvel table
 */

double get_velocity_model_surface_property(char property) {

    int index = 0;

    if (property == 'P')
        return (depth_Vp_Vs_rho[index][1]);
    else if (property == 'S')
        return (depth_Vp_Vs_rho[index][2]);
    else if (property == 'D')
        return (depth_Vp_Vs_rho[index][3]);

    return (-1.0);

}

/** depth_corr_mwpd_prem
 * returns Mwpd depth correction from depth_corr_mwpd_prem table
 */

double get_depth_corr_mwpd_prem(double depth) {

    int index = NUM_DEPTH_CORR_MWPD_PREM_DEPTH - 1;
    while (index >= 0 && depth_corr_mwpd_prem[index][0] > depth)
        index--;

    if (index < 0) // check if depth index valid
        return (-999.0);

    return (depth_corr_mwpd_prem[index][1]);

    return (-999.0);

}

/** returns value of Q(∆) for P waves for shallow shocks (h < 70 km) according to Gutenberg and Richter (1956a) if the ground amplitude is given in μm.
 */

double get_Q_Delta_PV_value(double delta, double depth) {

    if (depth > 70.0)
        return (-1.0);

    int index = NUM_Q_DELTA_DISTANCE - 1;
    while (index >= 0 && Q_Delta_PV[index][0] > delta)
        index--;

    if (index < 0) // check if distance index valid
        return (-1.0);

    return (Q_Delta_PV[index][1]);

    return (-1.0);

}


#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

/** returns value of Q(h,∆) for P waves according to Gutenberg and Richter (1956a) if the ground amplitude is given in μm.
 *  Q interpolation algorithm from Joachim Saul 20101223
 */

double get_Q_Depth_Delta_PV_value(double delta, double depth) {

    if (delta < Q_DEPTH_DELTA__DELTA_MIN_USE || delta > Q_DEPTH_DELTA__DELTA_MAX)
        return (-1.0);

    if (depth < 0.0)
        depth = 0.0;
    if (depth > Q_DEPTH_DELTA__DEPTH_MAX)
        return (-1.0);

    int idepth = 0;
    double s1 = 0.0;
    if (depth < 100.0) {
        idepth = (int) (depth / 25.0) + 2;
        s1 = 0.04 * (depth - 25.0 * ((double) (idepth - 2)));
    } else {
        idepth = min((int) (depth / 50.) + 4, 17);
        s1 = 0.02 * (depth - 50.0 * ((double) (idepth - 4)));
    }
    int idelta = max(min((int) (delta), 108), 2);
    double s2 = delta - (double) idelta;
    double q1 = Q_Depth_Delta_PV[idepth - 2][idelta - 2] + s1 * (Q_Depth_Delta_PV[idepth - 1][idelta - 2] - Q_Depth_Delta_PV[idepth - 2][idelta - 2]);
    double Q_value = q1 + s2 * ((Q_Depth_Delta_PV[idepth - 2][idelta - 1] + s1 * (Q_Depth_Delta_PV[idepth - 1][idelta - 1] - Q_Depth_Delta_PV[idepth - 2][idelta - 1])) - q1);

    return (Q_value);

}

/** returns string phase name for specified phase id
 */

char *phase_name_for_id(int phase_id) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return ("X");
    return (phases[phase_id]);

}

/** returns phase id for specified  string phase name
 */

int phase_id_for_name(char *phase_name) {

    int phase_id = 0;
    while (phase_id < NUM_TTIME_PHASES) {
        if (strcmp(phases[phase_id], phase_name) == 0)
            return (phase_id);
        phase_id++;
    }

    return (-1);

}

/** simple travel-time estimate
 * returns time from table for phase n_phasse
 */

double phases_dist_time(int phase_id, int ndist, int ndepth) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return (-1.0);

    if (ndist < 0 || ndist >= NUM_DIST_TIME)
        return (-1.0);

    if (ndepth < 0 || ndepth >= NUM_DIST_TIME_DEPTH)
        return (-1.0);

    if (!dist_time[phase_id][ndist][1]) // check if times available for this depth
        return (-1.0);

    return (dist_time[phase_id][ndist][ndepth + 2]);


}

/** accurate travel-time interpolated from table using 2D Lagrange interpolation
 * returns time in seconds or -1.0 if dist < 0.0 or greater than maximum distance in dist_time table
 * if depth < 0.0 depth is trucated to depth=0.0
 * if depth > max depth is trucated to max depth
 */

double get_ttime(int phase_id, double dist, double depth) {

    if (dist < 0.0 || dist > DIST_TIME_DIST_MAX) {
        //printf("get_time -1 dist: phase_id %d  dist %g  depth %g\n", phase_id, dist, depth);
        return (-1.0);
    }
    int ndist = (int) (dist / DIST_TIME_DIST_STEP);
    double ndist1;
    if (ndist < 0) {
        ndist = ndist1 = 0;
        dist = 0.0;
    } else if (ndist >= NUM_DIST_TIME - 1) {
        ndist = ndist1 = NUM_DIST_TIME - 1;
        dist = DIST_TIME_DIST_MAX;
    } else {
        ndist1 = ndist + 1;
    }
    double dist_diff = (dist - (double) ndist * DIST_TIME_DIST_STEP) / DIST_TIME_DIST_STEP;

    int ndepth = (int) (depth / DIST_TIME_DEPTH_STEP);
    double ndepth1;
    if (ndepth < 0) {
        ndepth = ndepth1 = 0;
        depth = 0.0;
    } else if (ndepth >= NUM_DIST_TIME_DEPTH - 1) {
        ndepth = ndepth1 = NUM_DIST_TIME_DEPTH - 1;
        depth = DIST_TIME_DEPTH_MAX;
    } else {
        ndepth1 = ndepth + 1;
    }
    double depth_diff = (depth - (double) ndepth * DIST_TIME_DEPTH_STEP) / DIST_TIME_DEPTH_STEP;

    double tval00, tval01, tval10, tval11;
    tval00 = tval01 = tval10 = tval11 = 0.0;
    if (dist_diff < 1.0 && depth_diff < 1.0) {
        tval00 = phases_dist_time(phase_id, ndist, ndepth);
        if (tval00 < 0.0) {
            //if (phase_id == 0) printf("get_time -1 tval00: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff < 1.0 && depth_diff > 0.0) {
        tval01 = phases_dist_time(phase_id, ndist, ndepth1);
        if (tval01 < 0.0) {
            //if (phase_id == 0) printf("get_time -1 tval01: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff > 0.0 && depth_diff < 1.0) {
        tval10 = phases_dist_time(phase_id, ndist1, ndepth);
        if (tval10 < 0.0) {
            //if (phase_id == 0) printf("get_time -1 tval10: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff > 0.0 && depth_diff > 0.0) {
        tval11 = phases_dist_time(phase_id, ndist1, ndepth1);
        if (tval11 < 0.0) {
            //if (phase_id == 0) printf("get_time -1 tval01: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }

    // 2D Lagrange interpolation
    double ttime =
            tval00 * (1.0 - dist_diff) * (1.0 - depth_diff)
            + tval01 * (1.0 - dist_diff) * depth_diff
            + tval10 * dist_diff * (1.0 - depth_diff)
            + tval11 * dist_diff * depth_diff;

    //printf("phase_id %d  dist %g  depth %g  ttval %g %g %g %g  depth_diff %g  depth_diff %g  ttime %g\n", phase_id, dist, depth, tval00, tval01, tval10, tval11, dist_diff, depth_diff, ttime);

    return (ttime);

}

/** simple take-off angle estimate
 * returns take-off angle from table for phase n_phasse
 *
 * IMPORTANT! = dist_toang table must have identical dimensions and structure as dist_time table
 */

double phases_dist_toang(int phase_id, int ndist, int ndepth) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return (-1.0);

    if (ndist < 0 || ndist >= NUM_DIST_TIME)
        return (-1.0);

    if (ndepth < 0 || ndepth >= NUM_DIST_TIME_DEPTH)
        return (-1.0);

    if (!dist_toang[phase_id][ndist][1]) // check if take-off angle available for this depth
        return (-1.0);

    return (dist_toang[phase_id][ndist][ndepth + 2]);


}

/** accurate take-off angle interpolated from table using 2D Lagrange interpolation
 * returns take-off angle in deg (0/down->180/up) or -1.0 if dist < 0.0 or greater than maximum distance in dist_toang table
 * if depth < 0.0 depth is trucated to depth=0.0
 * if depth > max depth is trucated to max depth
 *
 * IMPORTANT! = dist_toang table must have identical dimensions and structure as dist_time table
 */

double get_take_off_angle(int phase_id, double dist, double depth) {

    if (dist < 0.0 || dist > DIST_TIME_DIST_MAX) {
        //printf("get_toang -1 dist: phase_id %d  dist %g  depth %g\n", phase_id, dist, depth);
        return (-1.0);
    }
    int ndist = (int) (dist / DIST_TIME_DIST_STEP);
    double ndist1;
    if (ndist < 0) {
        ndist = ndist1 = 0;
        dist = 0.0;
    } else if (ndist >= NUM_DIST_TIME - 1) {
        ndist = ndist1 = NUM_DIST_TIME - 1;
        dist = DIST_TIME_DIST_MAX;
    } else {
        ndist1 = ndist + 1;
    }
    double dist_diff = (dist - (double) ndist * DIST_TIME_DIST_STEP) / DIST_TIME_DIST_STEP;

    int ndepth = (int) (depth / DIST_TIME_DEPTH_STEP);
    double ndepth1;
    if (ndepth < 0) {
        ndepth = ndepth1 = 0;
        depth = 0.0;
    } else if (ndepth >= NUM_DIST_TIME_DEPTH - 1) {
        ndepth = ndepth1 = NUM_DIST_TIME_DEPTH - 1;
        depth = DIST_TIME_DEPTH_MAX;
    } else {
        ndepth1 = ndepth + 1;
    }
    double depth_diff = (depth - (double) ndepth * DIST_TIME_DEPTH_STEP) / DIST_TIME_DEPTH_STEP;

    double tval00, tval01, tval10, tval11;
    tval00 = tval01 = tval10 = tval11 = 0.0;
    if (dist_diff < 1.0 && depth_diff < 1.0) {
        tval00 = phases_dist_toang(phase_id, ndist, ndepth);
        if (tval00 < 0.0) {
            //if (phase_id == 0) printf("get_toang -1 tval00: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff < 1.0 && depth_diff > 0.0) {
        tval01 = phases_dist_toang(phase_id, ndist, ndepth1);
        if (tval01 < 0.0) {
            //if (phase_id == 0) printf("get_toang -1 tval01: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff > 0.0 && depth_diff < 1.0) {
        tval10 = phases_dist_toang(phase_id, ndist1, ndepth);
        if (tval10 < 0.0) {
            //if (phase_id == 0) printf("get_toang -1 tval10: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }
    if (dist_diff > 0.0 && depth_diff > 0.0) {
        tval11 = phases_dist_toang(phase_id, ndist1, ndepth1);
        if (tval11 < 0.0) {
            //if (phase_id == 0) printf("get_toang -1 tval01: phase_id %d  dist %g  depth %g  ndist %d  ndepth %d\n", phase_id, dist, depth, ndist, ndepth);
            return (-1.0);
        }
    }

    // 2D Lagrange interpolation
    double toang =
            tval00 * (1.0 - dist_diff) * (1.0 - depth_diff)
            + tval01 * (1.0 - dist_diff) * depth_diff
            + tval10 * dist_diff * (1.0 - depth_diff)
            + tval11 * dist_diff * depth_diff;

    //printf("phase_id %d  dist %g  depth %g  ttval %g %g %g %g  depth_diff %g  depth_diff %g  toang %g\n", phase_id, dist, depth, tval00, tval01, tval10, tval11, dist_diff, depth_diff, toang);

    return (toang);

}

/** returns travel time error for specified phase
 */

double get_phase_ttime_error(int phase_id) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return (INVALID_ERROR);
    return (phase_ttime_error[phase_id]);

}

/** returns phase type flags for specified phase
 */

unsigned get_phase_type_flag(int phase_id) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return (INVALID_TYPE);
    return (phase_type_flags[phase_id]);

}

/** add a phase type flag to specified phase
 */

unsigned add_phase_type_flag_to_phase_name(char *phase_name, int iflag) {

    int phase_id = phase_id_for_name(phase_name);

    return (add_phase_type_flag(phase_id, iflag));

}

/** add a phase type flag to specified phase
 */

unsigned add_phase_type_flag(int phase_id, int iflag) {

    if (phase_id < 0 || phase_id >= NUM_TTIME_PHASES)
        return (INVALID_TYPE);

    phase_type_flags[phase_id] |= iflag;

    return (phase_type_flags[phase_id]);

}

/** simple distance estimate for first arrival P and S
 * returns approximate distance in degrees for a specified travel time, source depth and first arrival phase P or S
 */

double simple_distance(double ttime, double depth, char *phase_name, int *pphase_id_found) {

    int ndepth = (int) (0.5 + depth / DIST_TIME_DEPTH_STEP);
    if (ndepth < 0)
        ndepth = 0;
    if (ndepth >= NUM_DIST_TIME_DEPTH)
        ndepth = NUM_DIST_TIME_DEPTH - 1;

    int ndist_max = -1;

    int phase_id;
    for (phase_id = 0; phase_id < NUM_TTIME_PHASES; phase_id++) {

        if (!((strcmp(phase_name, phases[phase_id]) == 0)
                || (strcmp(phase_name, "$FA_P") == 0 && is_first_arrival_P(phase_id))
                || (strcmp(phase_name, "$FA_S") == 0 && is_first_arrival_S(phase_id))))
            continue;

        int ndist = NUM_DIST_TIME - 1;
        double tt_prev = phases_dist_time(phase_id, ndist, ndepth);
        ndist--;
        double tt_next = 0.0;
        while (ndist >= 0 && ndist > ndist_max) {
            tt_next = phases_dist_time(phase_id, ndist, ndepth);
            // check if requested ttime falls between adjacent times for phase_name
            if (tt_prev > 0.0 && tt_next > 0.0 && ((ttime < tt_prev && ttime >= tt_next) || (ttime >= tt_prev && ttime < tt_next)))
                break;
            tt_prev = tt_next;
            ndist--;
        }

        if (ndist >= 0 && ndist > ndist_max) {
            ndist_max = ndist;
            *pphase_id_found = phase_id;
        }
    }

    if (ndist_max < 0) {
        *pphase_id_found = -1;
        return (-1.0);
    }

    return ((double) ndist_max * DIST_TIME_DIST_STEP);

}

/** simple distance estimate for first arriving P
 * returns approximate distance in degrees for a specified travel time and source depth
 */

double simple_P_distance(double ttime, double depth, int *pphase_id_found) {

    return (simple_distance(ttime, depth, "$FA_P", pphase_id_found));

}

/** simple distance estimate for first arriving S
 * returns approximate distance in degrees for a specified travel time and source depth
 */

double simple_S_distance(double ttime, double depth, int *pphase_id_found) {

    return (simple_distance(ttime, depth, "$FA_S", pphase_id_found));

}

/** simple great-circle distance calculation on a sphere
 * returns distance in degrees
 *
 * Adapted from:
 *    Geographic Distance and Azimuth Calculations
 *    by Andy McGovern
 *    http://www.codeguru.com/algorithms/GeoCalc.html
 */

double GCDistance(double lat1, double lon1, double lat2, double lon2) {

    double d;

    if (lat1 == lat2 && lon1 == lon2)
        return (0.0);

    lat1 *= DE2RA;
    lon1 *= DE2RA;
    lat2 *= DE2RA;
    lon2 *= DE2RA;

    d = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2);
    // 20120928 AJL return (AVG_ERAD * acos(d) * KM2DEG);
    return (acos(d) * RA2DE);

}

/** simple great-circle azimuth calculation on a sphere
 *
 * returns azimuth in degrees for great-circle from lat1/lon1 to lat2/lon2
 *
 */

double GCAzimuth(double lat1, double lon1, double lat2, double lon2) {

    double lonA = lon1 * DE2RA;
    double latA = lat1 * DE2RA;
    double lonB = lon2 * DE2RA;
    double latB = lat2 * DE2RA;

    // distance
    double dist =
            acos(
            sin(latA) * sin(latB)
            + cos(latA) * cos(latB) * cos((lonB - lonA)));

    // 20141106 AJL - added following check, to prevent div by 0 of sin(dist))
    if (dist < FLT_MIN) {
        return (0.0);
    }

    // azimuth
    double cosAzimuth =
            (cos(latA) * sin(latB)
            - sin(latA) * cos(latB)
            * cos((lonB - lonA)))
            / sin(dist);
    double sinAzimuth =
            cos(latB) * sin((lonB - lonA)) / sin(dist);
    double az = atan2(sinAzimuth, cosAzimuth) / DE2RA;

    if (isnan(az) && fabs(lon2 - lon1) < 0.000001) {
        if (lat1 > lat2) {
            az = 180.0;
        } else {
            az = 0.0;
        }
    }

    if (az < 0.0) {
        az += 360.0;
    }

    return (az);
}


/** calculate end point (latitude/longitude) given a starting point, distance (deg), and azimuth (deg).
 *
 * Adapted from:
 *    Geographic Distance and Azimuth Calculations
 *    by Andy McGovern
 *    http://www.codeguru.com/algorithms/GeoCalc.html
 */
/*
 * Problem 1C. Calculate end point (latitude/longitude) given a starting point, distance, and azimuth
 * Given {lat1, lon1, distance, azimuth} calculate {lat2, lon2}. First, work backwards (relative to Problem 1A)
 * and find b from the distance by dividing by the Earth radius. b = distance / (Earth Radius) making sure distance
 * and (Earth Radius) are the same units so that we end up with b in radians. Knowing b, calculate a using
 * a = arccos(cos(b)*cos(90 - lat1) + sin(90 - lat1)*sin(b)*cos(azimuth))—basically taking the arc cosine of
 * the Law of Cosines for a. From a, we can get lat2, so the only item remaining is to figure lon2;
 * we can get that if we know B. Calculate B using B = arcsin(sin(b)*sin(azimuth)/sin(a)).
 * Then finally, lat2 = 90 - a and lon2 = B + lon1. Essentially, we just worked Problem 1A backwards.
void GCDistanceAzimuth(double lat1, double lon1, double dist, double az, double* lat2, double* lon2)
{
        double b = DEG2KM * dist / AVG_ERAD;
        double sinb = sin(b);
        double cosb = cos(b);
        double sinc = sin(DE2RA * (90.0 - lat1));
        double cosc = cos(DE2RA * (90.0 - lat1));
        double azrad = DE2RA * az;

        double a = acos(cosb*cosc + sinc*sinb*cos(azrad));
        double B = asin(sinb*sin(azrad)/sin(a));

 *lat2 = RA2DE * ((PI/2.0) - a);
 *lon2 = RA2DE * B + lon1;
}
 */

/** calculate end point (latitude/longitude) given a starting point, distance (deg), and azimuth (deg).
 *
 * Adapted from:
 *    http://openmap.bbn.com/doc/api/com/bbn/openmap/proj/GreatCircle.html
 *
 * Calculate point at azimuth and distance from another point.
 * <p>
 * Returns a LatLonPoint at arc distance `c' in direction `Az'
 * from start point.
 * <p>
 *
 * @param phi1 latitude in radians of start point
 * @param lambda0 longitude in radians of start point
 * @param c arc radius in radians (0 &lt; c &lt;= PI)
 * @param Az azimuth (direction) east of north (-PI &lt;= Az &lt;
 *        PI)
 * @return LatLonPoint
 *
 */
void PointAtGCDistanceAzimuth(double lat1, double lon1, double dist, double az, double* lat2, double* lon2) {

    double phi1 = DE2RA * lat1;
    double lambda0 = DE2RA * lon1;
    double c = DE2RA * dist;
    double Az = DE2RA * az;

    double cosphi1 = cos(phi1);
    double sinphi1 = sin(phi1);
    double cosAz = cos(Az);
    double sinAz = sin(Az);
    double sinc = sin(c);
    double cosc = cos(c);

    *lat2 = RA2DE * asin(sinphi1 * cosc + cosphi1 * sinc * cosAz);
    *lon2 = RA2DE * (atan2(sinc * sinAz, cosphi1 * cosc - sinphi1 * sinc * cosAz) + lambda0);

}
