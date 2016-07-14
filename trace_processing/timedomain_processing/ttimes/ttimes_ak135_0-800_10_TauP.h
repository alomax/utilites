// travel times in the ak135 model
//
// generated with edu.sc.seis.TauP Java
//

#define MODEL_TTIME_NAME "AK135 GLOBAL (ak135_0-800_10_TauP)"

// values here must match ttimes_ak135_0-800_10_TauP_times_phases
#define NUM_DIST_TIME 1801
#define DIST_TIME_DIST_STEP 0.1     // degrees
#define DIST_TIME_DIST_MAX ((double)(NUM_DIST_TIME-1)*DIST_TIME_DIST_STEP)
#define NUM_DIST_TIME_DEPTH 81
#define DIST_TIME_DEPTH_STEP 10.0   // km
#define DIST_TIME_DEPTH_MAX ((double)(NUM_DIST_TIME_DEPTH-1)*DIST_TIME_DEPTH_STEP)

// number and order here must match ttimes_ak135_0-800_10_TauP_times_phases
#define NUM_TTIME_PHASES 21

// set indices of important phases
// IMPORTANT!: index must match phase order in ttimes.c array phases[]
#define P_PHASE_INDEX 0
#define S_PHASE_INDEX 9



// phase_ttime_error
// 20111123 AJL - added factor to help avoid false events and missed phase associations
#define ERROR_FACTOR 1.0
//#define ERROR_FACTOR 2.0
//
#define INVALID_ERROR -1.0
#define SMALL_ERROR 4.0*ERROR_FACTOR
#define MED_ERROR 6.0*ERROR_FACTOR
// 20120209 AJL #define LARGE_ERROR 10.0*ERROR_FACTOR  // 20111215 AJL added to avoid later phases being used as P to give location of false events
//#define LARGE_ERROR MED_ERROR  // 20120209 AJL disable large since often hides P phases from real events
#define LARGE_ERROR 10.0*ERROR_FACTOR  // 20120209 AJL re-enable disable large since now only applies to Pdiff and Sdiff
//#define P_ERROR 1.0
#define P_ERROR 2.0     // 20110130 AJL changed to help associate more real P arrivals
//#define P_ERROR 4.0     // 20111228 AJL changed to help associate regional P arrivals
// set S errors large to avoid false events from poor S picks at regional distance
// set PP errors large to avoid false events from poor S picks at regional distance
// set Pdiff, Sdiff errors very large since LP signal - HF picks can be delayed (?)

// order here must match ttimes_ak135_0-800_10_TauP_times_phases
// P,p,Pg,Pn Pdiff Pg PKiKP PKP PKIKP PKKP PKPPKP PP S,s,Sg,Sn Sdiff Sg SKS SKP pP sP PcP ScP ScS PKKS SKKS
static double phase_ttime_error[NUM_TTIME_PHASES] = {
    P_ERROR, //P,
    LARGE_ERROR, //Pdiff,
    MED_ERROR, //Pg,    // 20111228 AJL changed to help associate regional P arrivals
    SMALL_ERROR, // PKiKP,
    SMALL_ERROR, //PKP,
    SMALL_ERROR, //PKIKP,
    SMALL_ERROR, //"PKKP,
    SMALL_ERROR, //"PKPPKP,
    SMALL_ERROR, //PP,
    SMALL_ERROR, //S,
    LARGE_ERROR, //Sdiff,
    SMALL_ERROR, //Sg",
    SMALL_ERROR, //SKS,
    SMALL_ERROR, //SKP,
    SMALL_ERROR, //pP",
    SMALL_ERROR, //sP",
    SMALL_ERROR, //PcP",
    SMALL_ERROR, //ScP",
    MED_ERROR, //ScS",
    MED_ERROR, //PKKS,
    MED_ERROR, //SKKS,
};


// IMPORTANT!: Must set set indices of important phases in ttimes.h: P_PHASE_INDEX, S_PHASE_INDEX, ...
static char phases[NUM_TTIME_PHASES][32] = {
// P,p,Pg,Pn Pdiff Pg PKiKP PKP PKIKP PKKP PKPPKP PP S,s,Sg,Sn Sdiff Sg SKS SKP pP sP PcP ScP ScS PKKS SKKS
    "P",
    "Pdiff",
    "Pg",
    "PKiKP",
    "PKP",
    "PKIKP",
    "PKKP",
    "PKPPKP",
    "PP",
    "S",
    "Sdiff",
    "Sg",
    "SKS",
    "SKP",
    "pP",
    "sP",
    "PcP",
    "ScP",
    "ScS",
    "PKKS",
    "SKKS",
};

// order here must match ttimes_ak135_0-800_10_TauP_times_phases
// FIRST_ARRIVAL_P from diagrams in 2005__Kennett__Seismological_Tables_ak135__Report
static unsigned phase_type_flags[NUM_TTIME_PHASES] = {
// P,p,Pg,Pn Pdiff Pg PKiKP PKP PKIKP PKKP PKPPKP PP S,s,Sg,Sn Sdiff Sg SKS SKP pP sP PcP ScP ScS PKKS SKKS
    COUNT_IN_LOCATION | DIRECT_P | FIRST_ARRIVAL_P | P_AT_SOURCE | P_AT_STATION, //P,
    DIRECT_P | FIRST_ARRIVAL_P | P_AT_SOURCE | P_AT_STATION, //Pdiff,  // 20121127 AJL - removed COUNT_IN_LOCATION flag
    COUNT_IN_LOCATION | DIRECT_P | P_AT_SOURCE | P_AT_STATION, // Pg,
    P_AT_SOURCE | P_AT_STATION, // PKiKP,
    COUNT_IN_LOCATION | FIRST_ARRIVAL_P | P_AT_SOURCE | P_AT_STATION, //PKP,
    COUNT_IN_LOCATION | FIRST_ARRIVAL_P | P_AT_SOURCE | P_AT_STATION, //PKIKP,
    P_AT_SOURCE | P_AT_STATION, //PKKP,
    P_AT_SOURCE | P_AT_STATION, //PKPPKP,
    P_AT_SOURCE | P_AT_STATION, //PP,
    DIRECT_S | FIRST_ARRIVAL_S | S_AT_SOURCE | S_AT_STATION, //S",
    DIRECT_S | FIRST_ARRIVAL_S | S_AT_SOURCE | S_AT_STATION, //Sdiff",
    S_AT_SOURCE | S_AT_STATION, //Sg",
    FIRST_ARRIVAL_S | S_AT_SOURCE | S_AT_STATION, //SKS,
    0, //SKP,
    P_AT_SOURCE | P_AT_STATION, //pP,
    0, //sP,
    P_AT_SOURCE | P_AT_STATION, //PcP,
    0, //ScP,
    S_AT_SOURCE | S_AT_STATION, //ScS,
    0, //PKKS,
    S_AT_SOURCE | S_AT_STATION, //SKKS,
};
// IMPORTANT: set following to 1 if multiple phases (e.g. P and S) have COUNT_IN_LOCATION flag and network scale is local/regional, 0 otherwise
//      This check will prevent overlap of origin times corresponding to different phases for each pick, preventing multiple phase association for a pick
#define CHECK_OVERLAP_COUNT_IN_LOCATION 0
