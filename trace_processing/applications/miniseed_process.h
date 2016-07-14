/***************************************************************************
 * msprocess_test.h:
 *
 * Generic routines to perform timedomain-processing processing on Mini-SEED records from a miniseed file.
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2009.02.03
 ***************************************************************************/

void msp_process(MSRecord *msr, int source_id, int is_new_source, char* msr_net, char* msr_sta, char* msr_loc, char* msr_chan, flag details, flag verbose);
time_t hptime_t2time_t(hptime_t hptime, double* pdsec);

