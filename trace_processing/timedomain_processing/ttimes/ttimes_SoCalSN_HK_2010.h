// travel times in the ak135 model
//
// iaspei-tau (http://www.iris.edu/pub/programs/iaspei-tau/)
// modified AJL, use ttimes_ajl.f
//

#define MODEL_TTIME_NAME "SoCalSN_HK_2010 REGIONAL (SoCalSN_HK_2010)"

// values here must match ttimes_regional_times_phases
#define NUM_DIST_TIME 2001
#define DIST_TIME_DIST_STEP 0.01     // degrees
#define DIST_TIME_DIST_MAX ((double)(NUM_DIST_TIME-1)*DIST_TIME_DIST_STEP)
#define NUM_DIST_TIME_DEPTH 81
#define DIST_TIME_DEPTH_STEP 1.0   // km
#define DIST_TIME_DEPTH_MAX ((double)(NUM_DIST_TIME_DEPTH-1)*DIST_TIME_DEPTH_STEP)

// number and order here must match ttimes_regional_times_phases
#define NUM_TTIME_PHASES 2

// set indices of important phases
// IMPORTANT!: index must match phase order in ttimes.c array phases[]
#define P_PHASE_INDEX 0
#define S_PHASE_INDEX 1



// phase_ttime_error
#define ERROR_FACTOR 1.0
//#define ERROR_FACTOR 2.0
//
#define INVALID_ERROR -1.0
#define P_ERROR ERROR_FACTOR
#define S_ERROR 2.0*ERROR_FACTOR

// order here must match ttimes_regional_times_phases
static double phase_ttime_error[NUM_TTIME_PHASES] = {
    P_ERROR, //P",
    S_ERROR, //S",
};


// IMPORTANT!: Must set set indices of important phases in ttimes.h: P_PHASE_INDEX, S_PHASE_INDEX, ...
static char phases[NUM_TTIME_PHASES][32] = {
    "P",
    "S",
};

// order here must match ttimes_regional_times_phases
// FIRST_ARRIVAL_P from diagrams in 2005__Kennett__Seismological_Tables_ak135__Report
static unsigned phase_type_flags[NUM_TTIME_PHASES] = {
    COUNT_IN_LOCATION | DIRECT_P | FIRST_ARRIVAL_P | P_AT_SOURCE | P_AT_STATION, //P",
    COUNT_IN_LOCATION | DIRECT_S | FIRST_ARRIVAL_S | S_AT_SOURCE | S_AT_STATION, //S",
};
// IMPORTANT: set following to 1 if multiple phases (e.g. P and S) have COUNT_IN_LOCATION flag and network scale is local/regional, 0 otherwise
//      This check will prevent overlap of origin times corresponding to different phases for each pick, preventing multiple phase association for a pick
#define CHECK_OVERLAP_COUNT_IN_LOCATION 1
