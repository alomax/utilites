/***************************************************************************
 * mseedrtstream_lib.h - Mini-SEED data sorted to simulate a real-time stream.
 *
 * Opens one or more user specified files, applys filtering criteria
 * and outputs any matched data sorted into a time order that
 * simulates a real-time stream.
 *
 * In general critical error messages are prefixed with "ERROR:" and
 * the return code will be 1.  On successful operation the return
 * code will be 0.
 *
 * Written by Chad Trabant, IRIS Data Management Center.
 *
 * modified 2013.087
 *
 * modified 2013.04.04 by AJL - Anthony Lomax (ALomax Scientific)
 *      - convert to library
 *      - removed libdali / dlconn support
 *
 ***************************************************************************/

/* ToDo? Restamp record start times to simulate current data flow */

/* Input/output file information containers */
typedef struct Filelink_s {
    char *infilename; /* Input file name */
    FILE *infp; /* Input file descriptor */
    struct Filelink_s *next;
} Filelink;

/* Mini-SEED record information structures */
typedef struct Record_s {
    struct Filelink_s *flp;
    off_t offset;
    int reclen;
    hptime_t starttime;
    hptime_t endtime;
    struct Record_s *prev;
    struct Record_s *next;
} Record;

/* Record map, holds Record structures for a given MSTrace */
typedef struct RecordMap_s {
    long long int recordcnt;
    struct Record_s *first;
    struct Record_s *last;
} RecordMap;


Filelink* getFilelist();
int readfiles(RecordMap *recmap);
int writerecords(RecordMap *recmap);
int sendrecord(char *recbuf, Record *rec);

int sortrecmap(RecordMap *recmap);
int recordcmp(Record *rec1, Record *rec2);

int processparam_flag_mseedrtstream(int argcount, char **argvec, int *poptind);
int processparam_files_mseedrtstream(int argcount, char **argvec, int optind);
int processparam_init_mseedrtstream();
char *getoptval(int argcount, char **argvec, int argopt);
hptime_t gethptime(void);
int setofilelimit(int limit);
int addfile(char *filename);
int addlistfile(char *filename);
int readregexfile(char *regexfile, char **pppattern);
void usage_mseedrtstream(void);
