/***************************************************************************
 * seedlink_test.h:
 *
 * Generic routines to perform timedomain-processing processing on Mini-SEED records from a SeedLink server.
 *
 * Written by Anthony Lomax
 *   ALomax Scientific www.alomax.net
 *
 * modified: 2009.02.03
 ***************************************************************************/

int slp_process(SLMSrecord* msr, int source_id, char* msr_net, char* msr_sta, char* msr_loc, char* msr_chan,
        int details, int verbose, double data_latency, SLCD* slconn_packet);
