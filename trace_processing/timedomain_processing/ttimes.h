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

// phase_type_flags
#define INVALID_TYPE ((unsigned)0)
#define DIRECT_P ((unsigned)1)
#define FIRST_ARRIVAL_P ((unsigned)2)
#define P_AT_SOURCE ((unsigned)4)
#define P_AT_STATION ((unsigned)8)
#define DIRECT_S ((unsigned)16)
#define FIRST_ARRIVAL_S ((unsigned)32)
#define S_AT_SOURCE ((unsigned)64)
#define S_AT_STATION ((unsigned)128)
#define COUNT_IN_LOCATION ((unsigned)256)

// following must be larger than NUM_TTIME_PHASES set in ttimes_* velocity models / travel-time grids / phase types headers
#define MAX_NUM_TTIME_PHASES 100

int is_P(int phase_id);
int is_P_at_source(int phase_id);
int is_P_at_station(int phase_id);
int is_direct_P(int n_phase);
int is_first_arrival_P(int n_phase);
int is_S_at_source(int phase_id);
int is_S_at_station(int phase_id);
int is_direct_S(int n_phase);
int is_first_arrival_S(int n_phase);
int is_count_in_location(int phase_id);
int doCheckOverlapCountInLocation();
char *get_model_ttime_name();
int get_num_ttime_phases();
double get_dist_time_dist_max();
double get_dist_time_depth_max();
double get_dist_time_depth_weight(double depth);
double get_dist_time_latlon_depth_weight(double depth, double lat, double lon);
int get_P_phase_index();
int get_S_phase_index();
double get_velocity_model_property(char property, double distance, double depth);
double get_velocity_model_surface_property(char property);
double get_depth_corr_mwpd_prem(double depth);
double get_Q_Delta_PV_value(double distance, double depth);
double get_Q_Depth_Delta_PV_value(double delta, double depth);
char *phase_name_for_id(int n_phase);
int phase_id_for_name(char *phase_name);
double get_phase_ttime_error(int phase_id);
unsigned get_phase_type_flag(int phase_id);
unsigned add_phase_type_flag_to_phase_name(char *phase_name, int iflag);
unsigned add_phase_type_flag(int phase_id, int iflag);
double get_ttime(int phase_id, double dist, double depth);
double get_take_off_angle(int phase_id, double dist, double depth);
double simple_distance(double ttime, double depth, char *phase_name, int *pn_phase_found);
double simple_P_distance(double ttime, double depth, int *pn_phase_found);
double simple_S_distance(double ttime, double depth, int *pn_phase_found);

//#define PI 3.14159265359
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#define PI M_PI
//#define RA2DE 57.2957795129
#define RA2DE (180.0/PI)
//#define DE2RA 0.01745329252
#define DE2RA (PI/180.0)
#define AVG_ERAD 6371.0
#define KM2DEG (90.0/10000.0)
#define DEG2KM (10000.0/90.0)
double GCDistance(double lat1, double lon1, double lat2, double lon2);
double GCAzimuth(double lat1, double lon1, double lat2, double lon2);
void PointAtGCDistanceAzimuth(double lat1, double lon1, double dist, double az, double* lat2, double* lon2);



