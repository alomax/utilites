/*
 * File:   loc2xml.h
 * Author: anthony
 *
 * Created on January 24, 2011, 11:37 AM
 */


// loc
int writeLocXML(char *xmlWriterUri, time_t time_max, char* agencyId, HypocenterDesc** hypo_list, int num_hypocenters, TimedomainProcessingData** data_list, int num_de_data,
        int iWriteArrivals, int iWriteUnAssociatedPicks, int printIgnoredData);
