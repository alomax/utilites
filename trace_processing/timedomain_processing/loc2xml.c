/*
 * Copyright (C) 1999-2011 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/**   loc2xml.c

        Program to convert location to XML

 */

/**-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/**
        history:

        ver 01    24JAN2011  AJL  Original version developed from testWriter.c
                  12Mar2014  AJL  Updated to conform to QuakeML standard

 */


/** References */

/**
 * section: xmlWriter
 * synopsis: use various APIs for the xmlWriter
 * purpose: tests a number of APIs for the xmlWriter, especially
 *          the various methods to write to a filename, to a memory
 *          buffer, to a new document, or to a subtree. It shows how to
 *          do encoding string conversions too. The resulting
 *          documents are then serialized.
 * usage: testWriter
 * test: testWriter ; for i in 1 2 3 4 ; do diff writer.xml writer$$i.res ; done ; rm writer*.res
 * author: Alfred Mickautsch

 * copy: see Copyright for the status of this software.

The MIT License

Copyright (c) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright © 2006 by the Open Source Initiative
The contents of this website are licensed under the Open Software License 2.1 or Academic Free License 2.1

OSI is a registered non-profit with 501(c)(3) status. Donating to OSI is one way to show your support.

 */


#include <stdio.h>
#include <string.h>
#include <time.h>

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

//#include "../geometry/geometry.h"
#include "../alomax_matrix/alomax_matrix.h"
//#include "../matrix_statistics/matrix_statistics.h"
#include "../picker/PickData.h"
#include "timedomain_processing_data.h"
#include "ttimes.h"
#include "location.h"
#include "timedomain_processing.h"
#include "../feregion/feregion.h"
#include "../response/response_lib.h"

#include "loc2xml.h"


// libxml, check
#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

static char resource_id_root[4096];
static char resource_id_event_parameters[4096];
static char resource_id_event[4096];
static char resource_id_origin[4096];
static char resource_id_tmp[4096];
static char datetime_tmp[128];

#define DOUBLE_NULL (-FLT_MAX)

#define FEREGION_STR_SIZE 4096
static char feregion_str[FEREGION_STR_SIZE];
static char comment_str[4096];

#define MAX_NUM_ATTRIBUTES 10
#define ATTRIBUTES_STR_LEN 128
//static char attNames[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN];
//static char attValues[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN];

// libxml
//#define MY_ENCODING "ISO-8859-1"
#define MY_ENCODING "utf-8"
int initLibXml(void);
int cleanUpLibXml(void);
xmlChar *ConvertInput(const char *in, const char *encoding);

// XMML
xmlTextWriterPtr startXMLDocument(const char *xmlWriterUri);
void endXMLDocument(xmlTextWriterPtr writer);
void startXMLEventElement(xmlTextWriterPtr writer, const char* agencyId, const char* publicIDChar, const time_t creation_time);
void endXMLElement(xmlTextWriterPtr writer);


// loc XMML
void writeXMLOrigin(xmlTextWriterPtr writer,
        char* publicIDChar, HypocenterDesc* phypo, int nhypo, TimedomainProcessingData** data_list, int num_de_data,
        int iWriteArrivals);
void writeXMLAssociatedArrivals(xmlTextWriterPtr writer, int nassoc, HypocenterDesc* phypo, TimedomainProcessingData** data_list, int num_de_data);
void writeXMLUnAssociatedPicks(xmlTextWriterPtr writer, TimedomainProcessingData** data_list, int num_de_data, int printIgnoredData);
void writeXMLArrival(xmlTextWriterPtr writer, HypocenterDesc* phypo, TimedomainProcessingData* deData, int ndata);
void writeXMLPick(xmlTextWriterPtr writer, char *name, TimedomainProcessingData* deData, int ndata);
void writeXMLwaveformID(xmlTextWriterPtr writer, TimedomainProcessingData* deData);

void writeXMLCreationInfo(xmlTextWriterPtr writer, char *name, const char* agencyID, const time_t creation_time);
void writeXMLDateTime(xmlTextWriterPtr writer, char *name, int year, int month, int day, int hour, int min, double sec);
void writeXMLTimeQuantity(xmlTextWriterPtr writer, char *name, int year, int month, int day, int hour, int min, double sec, double uncertainty, double confidenceLevel);
void writeXMLOriginUncertainty(xmlTextWriterPtr writer, HypocenterDesc* phypo);
void writeXMLConfidenceEllipsoid(xmlTextWriterPtr writer, HypocenterDesc* phypo);
void writeXMLMagnitude(xmlTextWriterPtr writer, char* originID, char* name, double dvalue, double lower_uncertainty, double upper_uncertainty, int num_values);

void writeXMLRealQuantity(xmlTextWriterPtr writer, char* name, double dvalue, double uncertainty, double lower_uncertainty, double upper_uncertainty, double confidence_level);

void writeXMLDValue(xmlTextWriterPtr writer, char* name, double dvalue);
void writeXMLDValueWithAttributes(xmlTextWriterPtr writer, char* name, double dvalue,
        char attNames[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN], char attValues[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN], const int numAttributes);
void writeXMLIValue(xmlTextWriterPtr writer, char* name, int ivalue);
void writeXMLCValue(xmlTextWriterPtr writer, const char* name, const char* cvalue);

/** Early-est extensions
 *
 * xmlns:ee="http://net.alomax/earlyest/xmlns/ee
 *
 */
void writeXML_EarlyEst_Discriminant(xmlTextWriterPtr writer, char* originID, char* name, double dvalue, double lower_uncertainty, double upper_uncertainty, int num_values);
void writeXML_EarlyEst_qualityIndicators(xmlTextWriterPtr writer, HypocenterDesc* phypo);

/**
 * writeLocXML
 *
 * @param xmlWriterUri
 * @param publicID
 * @param phypo
 * @param nhypo
 * @param data_list
 * @param num_de_data
 * @param iWriteArrivals
 * @return
 */
int writeLocXML(char *xmlWriterUri, time_t time_stamp, char* agencyId, HypocenterDesc** hypo_list, int num_hypocenters, TimedomainProcessingData** data_list, int num_de_data,
        int iWriteArrivals, int iWriteUnAssociatedPicks, int printIgnoredData) {

    int rc;

    // this initialize the library
    initLibXml();

    // set public_id root
    sprintf(resource_id_root, "smi:%s/ee", agencyId);
    //resource_id_root[0] = '\0';

    // open new XMML document
    xmlTextWriterPtr writer = startXMLDocument(xmlWriterUri);

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "eventParameters");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
    }

    sprintf(resource_id_event_parameters, "%s/event_parameters/%ld", resource_id_root, time_stamp);

    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST resource_id_event_parameters);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return (-1);
    }

    time_t current_time = time(&current_time);
    writeXMLCreationInfo(writer, "creationInfo", agencyId, current_time);

    int nhyp;
    for (nhyp = 0; nhyp < num_hypocenters; nhyp++) {
        HypocenterDesc* phypo = hypo_list[nhyp];
        if (phypo->hyp_assoc_index >= 0) { // associated

            // start new XML event
            sprintf(resource_id_event, "%s/event/%ld", resource_id_root, phypo->unique_id);
            startXMLEventElement(writer, agencyId, resource_id_event, phypo->first_assoc_time);

            // QuakeML
            sprintf(resource_id_origin, "%s/origin/%ld", resource_id_root, phypo->unique_id);
            //T50Ex n Td n TdT50Ex WL Mwp n mb n T0 n Mwpd n
            writeXMLMagnitude(writer, resource_id_origin, "mb", phypo->mbLevelStatistics.centralValue, phypo->mbLevelStatistics.lowerBound,
                    phypo->mbLevelStatistics.upperBound, phypo->mbLevelStatistics.numLevel);
            writeXMLMagnitude(writer, resource_id_origin, "Mwp", phypo->mwpLevelStatistics.centralValue, phypo->mwpLevelStatistics.lowerBound,
                    phypo->mwpLevelStatistics.upperBound, phypo->mwpLevelStatistics.numLevel);
            writeXMLMagnitude(writer, resource_id_origin, "Mwpd", phypo->mwpdLevelStatistics.centralValue, phypo->mwpdLevelStatistics.lowerBound,
                    phypo->mwpdLevelStatistics.upperBound, phypo->mwpdLevelStatistics.numLevel);
            writeXMLOrigin(writer, resource_id_origin, phypo, nhyp, data_list, num_de_data, iWriteArrivals);

            // Early-est
            writeXML_EarlyEst_Discriminant(writer, resource_id_origin, "T50Ex", phypo->t50ExLevelStatistics.centralValue, phypo->t50ExLevelStatistics.lowerBound,
                    phypo->t50ExLevelStatistics.upperBound, phypo->t50ExLevelStatistics.numLevel);
            writeXML_EarlyEst_Discriminant(writer, resource_id_origin, "Td", phypo->taucLevelStatistics.centralValue, phypo->taucLevelStatistics.lowerBound,
                    phypo->taucLevelStatistics.upperBound, phypo->taucLevelStatistics.numLevel);
            writeXML_EarlyEst_Discriminant(writer, resource_id_origin, "TdT50Ex", phypo->tdT50ExLevelStatistics.centralValue, phypo->tdT50ExLevelStatistics.lowerBound,
                    phypo->tdT50ExLevelStatistics.upperBound, phypo->tdT50ExLevelStatistics.numLevel);
            writeXML_EarlyEst_Discriminant(writer, resource_id_origin, "T0", phypo->t0LevelStatistics.centralValue, phypo->t0LevelStatistics.lowerBound,
                    phypo->t0LevelStatistics.upperBound, phypo->t0LevelStatistics.numLevel);

            endXMLElement(writer);
        }
    }

    // Close the element. - eventParameters
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
    }

    // write un-associated picks
    if (iWriteUnAssociatedPicks) {

        // Start an element.
        /*
        rc = xmlTextWriterStartElement(writer, BAD_CAST "unAssociatedPicks");
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterStartElement\n");
        }*/

        writeXMLUnAssociatedPicks(writer, data_list, num_de_data, printIgnoredData);

        // Close the element.
        /*
        rc = xmlTextWriterEndElement(writer);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterEndElement\n");
        }*/

    }


    // end XML document
    endXMLDocument(writer);

    // Cleanup function for the XML library.
    cleanUpLibXml();

    /*
     * this is to debug memory for regression tests
     */
    xmlMemoryDump();

    return 0;
}

/**
 * writeXMLOrigin:
 *
 */
void writeXMLOrigin(xmlTextWriterPtr writer,
        char* publicIDChar, HypocenterDesc* phypo, int nhyp, TimedomainProcessingData** data_list, int num_de_data,
        int iWriteArrivals) {

    int rc;

    // set region
    feregion(phypo->lat, phypo->lon, feregion_str, FEREGION_STR_SIZE);

    // prepare origin time values
    time_t itime_sec = (time_t) phypo->otime;
    struct tm* tm_time = gmtime(&itime_sec);
    double dec_sec = phypo->otime - itime_sec;

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "origin");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "main", BAD_CAST "true");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }*/

    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST publicIDChar);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }

    xmlChar *tmp_str = NULL;
    sprintf(comment_str, "%s;  %4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%6.3f;  %s",
            publicIDChar, tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec, feregion_str);
    strcat(comment_str, " ");
    /* Write a comment_str as child of current element.
     * Please observe, that the input to the xmlTextWriter functions
     * HAS to be in UTF-8, even if the output XML is encoded in iso-8859-1 */
    tmp_str = ConvertInput(comment_str, MY_ENCODING);
    rc = xmlTextWriterWriteComment(writer, tmp_str);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteComment\n");
        return;
    }
    if (tmp_str != NULL) xmlFree(tmp_str);

    // write location children

    writeXMLTimeQuantity(writer, "time", tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec,
            phypo->ot_std_dev, 68.0);
    writeXMLRealQuantity(writer, "latitude", phypo->lat, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL);
    writeXMLRealQuantity(writer, "longitude", phypo->lon, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL);
    writeXMLRealQuantity(writer, "depth", phypo->depth, phypo->errz, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL);

    writeXMLCValue(writer, "region", feregion_str);

    //event_id n nP phs dist_min dist_max gap1 gap2 sigma_ot origin-time lat lon errH depth errZ

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "quality");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    writeXMLIValue(writer, "associatedPhaseCount", phypo->nassoc);
    writeXMLIValue(writer, "usedPhaseCount", phypo->nassoc_P);
    writeXMLDValue(writer, "standardError", phypo->ot_std_dev);
    writeXMLDValue(writer, "azimuthalGap", phypo->gap_primary);
    writeXMLDValue(writer, "secondaryAzimuthalGap", phypo->gap_secondary);
    writeXMLDValue(writer, "minimumDistance", phypo->dist_min);
    writeXMLDValue(writer, "maximumDistance", phypo->dist_max);
    // 20140217 AJL - added as provisional, non-standard QML to report amplitude attenuation regression
    writeXMLRealQuantity(writer, "ee:amplitudeAttenuationPower", phypo->linRegressPower.power, phypo->linRegressPower.power_dev, DOUBLE_NULL, DOUBLE_NULL, 68.0);
    writeXMLRealQuantity(writer, "ee:amplitudeAttenuationIntercept", phypo->linRegressPower.constant, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL, DOUBLE_NULL);
    // 20151120 AJL - added ee quality information
    writeXML_EarlyEst_qualityIndicators(writer, phypo);
    // 20160526 AJL - added ee status indicators
    writeXMLIValue(writer, "ee:nstaHasBeenActive", phypo->nstaHasBeenActive);
    writeXMLIValue(writer, "ee:nstaIsActive", phypo->nstaIsActive);
    long first_assoc_latency = phypo->first_assoc_time - (long) (0.5 + phypo->otime);
    first_assoc_latency = first_assoc_latency < 0 ? -1 : first_assoc_latency;
    writeXMLIValue(writer, "ee:firstAssocLatency", first_assoc_latency);
    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

    writeXMLOriginUncertainty(writer, phypo);

    // write arrivals
    if (iWriteArrivals) {
        writeXMLAssociatedArrivals(writer, phypo->hyp_assoc_index + 1, phypo, data_list, num_de_data);
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLAssociatedArrivals:
 *
 */
void writeXMLAssociatedArrivals(xmlTextWriterPtr writer, int nassoc, HypocenterDesc* phypo, TimedomainProcessingData** data_list, int num_de_data) {


    int ndata;
    for (ndata = 0; ndata < num_de_data; ndata++) {
        TimedomainProcessingData* deData = data_list[ndata];
        if (deData->is_associated == nassoc) {
            writeXMLArrival(writer, phypo, deData, ndata);
        }
    }

}

/**
 * writeXMLUnAssociatedPicks:
 *
 */
void writeXMLUnAssociatedPicks(xmlTextWriterPtr writer, TimedomainProcessingData** data_list, int num_de_data, int printIgnoredData) {

    int rc;

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "qrt:eventParameters");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    xmlChar *tmp_str = NULL;
    sprintf(comment_str, "Unassociated picks");
    strcat(comment_str, " ");
    /* Write a comment_str as child of current element.
     * Please observe, that the input to the xmlTextWriter functions
     * HAS to be in UTF-8, even if the output XML is encoded in iso-8859-1 */
    tmp_str = ConvertInput(comment_str, MY_ENCODING);
    rc = xmlTextWriterWriteComment(writer, tmp_str);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteComment\n");
        return;
    }
    if (tmp_str != NULL) xmlFree(tmp_str);

    int ndata;
    for (ndata = 0; ndata < num_de_data; ndata++) {
        TimedomainProcessingData* deData = data_list[ndata];
        if (!deData->is_associated
                && (printIgnoredData || (!deData->flag_clipped && !(deData->flag_snr_hf_too_low && deData->flag_snr_brb_too_low) && !deData->flag_a_ref_not_ok))
                ) {
            writeXMLPick(writer, "pick", deData, ndata);
        }
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }


}

/**
 * writeXMLArrival:
 *
 */
void writeXMLArrival(xmlTextWriterPtr writer, HypocenterDesc* phypo, TimedomainProcessingData* deData, int ndata) {

    int rc;


    // event_id n event dist az channel stream loc time unc pol phase residual tot_wt dist_wt st_q_wt T50 Aref T50/Aref Tdom s/n_HF s/n_BRB Mwp mb T0 Mwpd status

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "arrival");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    // Add an attribute.
    sprintf(resource_id_tmp, "%s/origin/%ld/arrival/%ld", resource_id_root, phypo->unique_id, (long) ndata);
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST resource_id_tmp);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }

    // QuakeML
    writeXMLCValue(writer, "phase", deData->phase);
    writeXMLDValue(writer, "azimuth", deData->epicentral_azimuth);
    writeXMLDValue(writer, "distance", deData->epicentral_distance);
    writeXMLDValue(writer, "timeResidual", deData->residual);
    writeXMLDValue(writer, "timeWeight", deData->loc_weight);
    writeXMLDValue(writer, "ee:takeOffAngleInc", deData->take_off_angle_inc); // 20160527 AJL - added

    // Pick
    writeXMLPick(writer, "pick", deData, ndata);

    // Early-est
    writeXMLDValue(writer, "ee:totalWeight", deData->loc_weight);
    writeXMLDValue(writer, "ee:distanceWeight", deData->dist_weight);
    writeXMLDValue(writer, "ee:staQualityWeight", deData->station_quality_weight);
    writeXMLDValue(writer, "ee:Aerr", deData->amplitude_error_ratio);
    writeXMLDValue(writer, "ee:mb", deData->mb->mag);
    writeXMLDValue(writer, "ee:Mwp", deData->mwp->mag);
    if (deData->mwpd != NULL)
        writeXMLDValue(writer, "ee:MwpdCorr", deData->mwpd->corr_mag);

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLPick:
 *
 */
void writeXMLPick(xmlTextWriterPtr writer, char* name, TimedomainProcessingData* deData, int ndata) {

    int rc;


    // event_id n event dist az channel stream loc time unc pol phase residual tot_wt dist_wt st_q_wt T50 Aref T50/Aref Tdom s/n_HF s/n_BRB Mwp mb T0 Mwpd status

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    // Add an attribute.
    sprintf(resource_id_tmp, "%s/pick/%s_%s_%s_%s/%ld", resource_id_root,
            deData->network, deData->station, deData->location, deData->channel,
            (long) (((double) deData->t_time_t + deData->t_decsec) * 1000.0)); // 1/1000 sec precision
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST resource_id_tmp);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }

    // QuakeML
    writeXMLwaveformID(writer, deData);
    // 20160527 AJL - added/modified (from timedomain_processing/timedomain_processing_report.c)
    double fmquality = 0.0;
    int fmpolarity = POLARITY_UNKNOWN;
    char fmtype[32] = "Err";
    setPolarity(deData, &fmquality, &fmpolarity, fmtype);
    writeXMLIValue(writer, "polarity", fmpolarity);
    writeXMLCValue(writer, "ee:polarityType", fmtype);
    writeXMLDValue(writer, "ee:polarityWeight", fmquality);
    // END - 20160527 AJL - added/modified
    writeXMLTimeQuantity(writer, "time", deData->year, deData->month, deData->mday, deData->hour, deData->min, deData->dsec, deData->pick_error, 68.0);
    //int numAttributes = 1;
    //strcpy(attNames[0], "type");
    //strcpy(attValues[0], deData->pick_error_type);
    //writeXMLDValueWithAttributes(writer, "uncertainty", deData->pick_error, attNames, attValues, numAttributes);

    // Early-est
    writeXMLCValue(writer, "ee:stream", pick_stream_name(deData));
    writeXMLDValue(writer, "ee:deltaTime", deData->deltaTime);
    writeXMLCValue(writer, "ee:useForLoc", deData->use_for_location ? "true" : "false");
    writeXMLDValue(writer, "ee:stationQualityWeight", deData->station_quality_weight);
    // 20160527 AJL - added/modified (from timedomain_processing/timedomain_processing_report.c)
    double t50_a_ref_have_gain_flag = 1.0; // t50 and a_ref values are positive if not using amplitude attenuation
    if (!(chan_resp[deData->source_id].have_gain && chan_resp[deData->source_id].responseType == DERIVATIVE_TYPE)) {
        t50_a_ref_have_gain_flag = -1.0; // t50 and a_ref values are negative if not corrected for gain
    }
    writeXMLDValue(writer, "ee:T50", deData->t50 * t50_a_ref_have_gain_flag);
    writeXMLDValue(writer, "ee:Aref", deData->a_ref * t50_a_ref_have_gain_flag);
    writeXMLDValue(writer, "ee:Aerr", deData->amplitude_error_ratio);
    double deLevel = 0.0;
    if (!isnan(deData->t50) && !isnan(deData->a_ref) && deData->a_ref > FLT_MIN) // avoid divide by zero
        deLevel = deData->t50 / deData->a_ref;
    writeXMLDValue(writer, "ee:T50Ex", deLevel);
    writeXMLDValue(writer, "ee:Td", deData->tauc_peak);
    writeXMLDValue(writer, "ee:T0", deData->t0->duration_plot);
    // 20160527 AJL - added/modified (from timedomain_processing/timedomain_processing_report.c)
    double grd_vel_peak_amp = GRD_MOT_INVALID;
    double grd_disp_peak_amp = GRD_MOT_INVALID;
    if (!deData->flag_snr_brb_too_low || !deData->flag_snr_brb_int_too_low) {
        if (!deData->flag_snr_brb_too_low) {
            grd_vel_peak_amp = deData->grd_mot->peak_amp_vel;
        }
        if (!deData->flag_snr_brb_int_too_low) {
            grd_disp_peak_amp = deData->grd_mot->peak_amp_disp;
        }
    }
    writeXMLDValue(writer, "ee:Avel", grd_vel_peak_amp);
    writeXMLDValue(writer, "ee:Adisp", grd_disp_peak_amp);
    // END -  20160527 AJL - added/modified
    double snr_hf = deData->a_ref < 0.0 || deData->sn_pick < FLT_MIN ? 0.0 : deData->a_ref / deData->sn_pick;
    writeXMLDValue(writer, "ee:snHF", snr_hf);
    double snr_brb = deData->sn_brb_signal < 0.0 || deData->sn_brb_pick < FLT_MIN ? 0.0 : deData->sn_brb_signal / deData->sn_brb_pick;
    writeXMLDValue(writer, "ee:snBRB", snr_brb);
    double snr_brb_int = deData->sn_brb_int_signal < 0.0 || deData->sn_brb_int_pick < FLT_MIN ? 0.0 : deData->sn_brb_int_signal / deData->sn_brb_int_pick;
    writeXMLDValue(writer, "ee:snBRBD", snr_brb_int);
    writeXMLDValue(writer, "ee:staCorr", deData->sta_corr); // 20160527 AJL - added

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLPick:
 *
 */
void writeXMLwaveformID(xmlTextWriterPtr writer, TimedomainProcessingData* deData) {

    int rc;


    // event_id n event dist az channel stream loc time unc pol phase residual tot_wt dist_wt st_q_wt T50 Aref T50/Aref Tdom s/n_HF s/n_BRB Mwp mb T0 Mwpd status

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "waveformID");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    writeXMLCValue(writer, "networkCode", deData->network);
    writeXMLCValue(writer, "stationCode", deData->station);
    writeXMLCValue(writer, "locationCode", deData->location);
    writeXMLCValue(writer, "channelCode", deData->channel);

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLTimeQuantity:
 *
 */
void writeXMLTimeQuantity(xmlTextWriterPtr writer,
        char *name,
        int year, int month, int day, int hour, int min, double sec,
        double uncertainty, double confidenceLevel
        ) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "data_timezone", BAD_CAST data_timezone);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    // <value>2007-10-10T14:40:39.055</value>
    sprintf(datetime_tmp, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%06.3f", year, month, day, hour, min, sec);

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "value", "%s", datetime_tmp);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    if (uncertainty != DOUBLE_NULL)
        writeXMLDValue(writer, "uncertainty", uncertainty);
    if (confidenceLevel != DOUBLE_NULL)
        writeXMLDValue(writer, "confidenceLevel", confidenceLevel);


    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLDateTime:
 *
 */
void writeXMLDateTime(xmlTextWriterPtr writer, char *name, int year, int month, int day, int hour, int min, double sec) {

    int rc;


    // <value>2007-10-10T14:40:39.055</value>
    sprintf(datetime_tmp, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%06.3f", year, month, day, hour, min, sec);

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST name, "%s", datetime_tmp);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

}

/**
 * writeXMLOriginUncertainty:
 *
 */
void writeXMLOriginUncertainty(xmlTextWriterPtr writer, HypocenterDesc* phypo) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "originUncertainty");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    /*horizontalUncertainty: float [0..1]
     * +	minHorizontalUncertainty: float [0..1]
     * +	maxHorizontalUncertainty: float [0..1]
     * +        azimuthMaxHorizontalUncertainty: float [0..1]
     * +	confidenceEllipsoid: ConfidenceEllipsoid [0..1]
     * +	preferredDescription: OriginUncertaintyDescription [0..1]
     * +	confidenceLevel: float [0..1]*/

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "minHorizontalUncertainty", "%g", phypo->ellipse.len1);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "maxHorizontalUncertainty", "%g", phypo->ellipse.len2);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    double azMaxHorUnc = phypo->ellipse.az1 + 90.0;
    if (azMaxHorUnc >= 360.0)
        azMaxHorUnc -= 360.0;
    if (azMaxHorUnc >= 180.0)
        azMaxHorUnc -= 180.0;
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "azimuthMaxHorizontalUncertainty", "%g", azMaxHorUnc);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "confidenceLevel", "%g", CONFIDENCE_LEVEL);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    writeXMLConfidenceEllipsoid(writer, phypo);


    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }


}

/**
 * writeXMLCreationInfo:
 *
 */
void writeXMLCreationInfo(xmlTextWriterPtr writer, char *name, const char* agencyID, const time_t creation_time) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    struct tm* tm_time = gmtime(&creation_time);
    double dec_sec = 0.0;
    writeXMLDateTime(writer, "creationTime", tm_time->tm_year + 1900, tm_time->tm_mon + 1, tm_time->tm_mday, tm_time->tm_hour, tm_time->tm_min, (double) tm_time->tm_sec + dec_sec);

    writeXMLCValue(writer, "agencyID", agencyID);

    sprintf(comment_str, "%s %s (%s)", EARLY_EST_MONITOR_SHORT_NAME, EARLY_EST_MONITOR_VERSION, EARLY_EST_MONITOR_VERSION_DATE);
    writeXMLCValue(writer, "version", comment_str);

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLConfidenceEllipsoid:
 *
 */
void writeXMLConfidenceEllipsoid(xmlTextWriterPtr writer, HypocenterDesc* phypo) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "ee:confidenceEllipsoidNLL");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    double semiMajorAxisLength, majorAxisPlunge, majorAxisAzimuth,
            semiIntermediateAxisLength, intermediateAxisPlunge, intermediateAxisAzimuth,
            semiMinorAxisLength;

    nllEllipsiod2XMLConfidenceEllipsoid(&(phypo->ellipsoid),
            &semiMajorAxisLength, &majorAxisPlunge, &majorAxisAzimuth,
            &semiIntermediateAxisLength, &intermediateAxisPlunge, &intermediateAxisAzimuth,
            &semiMinorAxisLength);


    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "semiMajorAxisLength", "%g", semiMajorAxisLength);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "majorAxisAzimuth", "%g", majorAxisAzimuth);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "majorAxisPlunge", "%g", majorAxisPlunge);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "semiIntermediateAxisLength", "%g", semiIntermediateAxisLength);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "intermediateAxisAzimuth", "%g", intermediateAxisAzimuth);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "intermediateAxisPlunge", "%g", intermediateAxisPlunge);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "semiMinorAxisLength", "%g", semiMinorAxisLength);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "confidenceLevel", "%g", CONFIDENCE_LEVEL);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }


    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * startXMLDocument:
 * @uri: the output URI
 *
 */
xmlTextWriterPtr startXMLDocument(const char *xmlWriterUri) {
    int rc;
    xmlTextWriterPtr writer;

    // Create a new XmlWriter for xmlWriterUri, with no compression.
    writer = xmlNewTextWriterFilename(xmlWriterUri, 0);
    if (writer == NULL) {
        printf("loc2xml: Error creating the xml writer\n");
        return (writer);
    }

    // Set indenting on
    rc = xmlTextWriterSetIndent(writer, 1);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterSetIndent\n");
        return (writer);
    }

    /* Start the document with the xml default for the version,
     * encoding ISO 8859-1 and the default for the standalone
     * declaration. */
    rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartDocument\n");
        return (writer);
    }

    /* QuakeML
    // Start a namespace element.
    rc = xmlTextWriterStartElementNS(writer, BAD_CAST "q", BAD_CAST "quakeml", BAD_CAST "http://quakeml.org/xmlns/quakeml/1.2");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return (writer);
    }
    // Start a namespace element.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns", BAD_CAST "http://quakeml.org/xmlns/bed/1.2");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return (writer);
    }
     */

    // Start a namespace element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "q:quakeml");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return (writer);
    }

    /* Add an attribute. */
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns", BAD_CAST "http://quakeml.org/xmlns/bed/1.2");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return (writer);
    }

    /* Add an attribute. */
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns:q", BAD_CAST "http://quakeml.org/xmlns/quakeml/1.2");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return (writer);
    }

    /* Add an attribute. */
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns:qrt", BAD_CAST "http://quakeml.org/xmlns/quakeml-rt/1.2");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return (writer);
    }

    /* Add an attribute. */
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns:ee", BAD_CAST "http://net.alomax/earlyest/xmlns/ee");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return (writer);
    }

    return (writer);

}

/**
 * loc2xml:
 * @uri: the output URI
 *
 * test the xmlWriter interface when writing to a new file
 */
void startXMLEventElement(xmlTextWriterPtr writer, const char* agencyId, const char* publicIDChar, const time_t creation_time) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "event");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST publicIDChar);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }

    writeXMLCreationInfo(writer, "creationInfo", agencyId, creation_time);

}

/**
 * writeXMLDValue:
 *
 */
void writeXMLDValue(xmlTextWriterPtr writer, char* name, double dvalue) {
    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }*/

    // Write formatted CDATA.
    rc = xmlTextWriterWriteFormatRaw(writer, "%g", dvalue);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatRaw\n");
        return;
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLRealQuantity:
 *
 */
void writeXMLMagnitude(xmlTextWriterPtr writer, char* originID, char* type, double dvalue, double lower_uncertainty, double upper_uncertainty, int num_values) {

    // 20140530 AJL - do not write this element if no data
    if (num_values < 1) {
        return;
    }

    int rc;

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "magnitude");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    // Add an attribute.
    sprintf(resource_id_tmp, "%s/origin/magnitude/%s", resource_id_root, type);
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "publicID", BAD_CAST resource_id_tmp);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }


    writeXMLCValue(writer, "originID", originID);

    writeXMLRealQuantity(writer, "mag", dvalue, DOUBLE_NULL, lower_uncertainty, upper_uncertainty, DOUBLE_NULL);

    writeXMLCValue(writer, "type", type);

    writeXMLIValue(writer, "stationCount", num_values);



    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLRealQuantity:
 *
 */
void writeXMLAmplitude(xmlTextWriterPtr writer,
        char* type,
        double genericAmplitude, double lower_uncertainty, double upper_uncertainty,
        char *category,
        char *unit,
        char *methodID
        ) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "amplitude");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    writeXMLRealQuantity(writer, "genericAmplitude", genericAmplitude, DOUBLE_NULL, lower_uncertainty, upper_uncertainty, DOUBLE_NULL);

    writeXMLCValue(writer, "type", type);
    writeXMLCValue(writer, "category", category);
    writeXMLCValue(writer, "unit", unit);
    writeXMLCValue(writer, "methodID", methodID);


    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLRealQuantity:
 *
 */
void writeXMLRealQuantity(xmlTextWriterPtr writer, char* name, double dvalue, double uncertainty, double lower_uncertainty, double upper_uncertainty, double confidence_level) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "value", "%g", dvalue);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }

    // Write an element as child.
    if (uncertainty != DOUBLE_NULL) {
        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "uncertainty", "%g", uncertainty);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
            return;
        }
    }

    // Write an element as child.
    if (lower_uncertainty != DOUBLE_NULL) {
        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "lowerUncertainty", "%g", lower_uncertainty);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
            return;
        }
    }

    // Write an element as child.
    if (upper_uncertainty != DOUBLE_NULL) {
        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "upperUncertainty", "%g", upper_uncertainty);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
            return;
        }
    }


    // Write an element as child.
    if (confidence_level != DOUBLE_NULL) {
        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "confidenceLevel", "%g", confidence_level);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
            return;
        }
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLDValueWithAttributes:
 *
 */
void writeXMLDValueWithAttributes(xmlTextWriterPtr writer, char* name, double dvalue,
        char attNames[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN], char attValues[MAX_NUM_ATTRIBUTES][ATTRIBUTES_STR_LEN], const int numAttributes) {

    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }

    int n;
    for (n = 0; n < numAttributes; n++) {
        // Add an attribute.
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST attNames[n], BAD_CAST attValues[n]);
        if (rc < 0) {
            printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
            return;
        }
    }

    // Write formatted CDATA.
    rc = xmlTextWriterWriteFormatRaw(writer, "%g", dvalue);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatRaw\n");
        return;
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLIValue:
 *
 */
void writeXMLIValue(xmlTextWriterPtr writer, char* name, int ivalue) {
    int rc;


    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST name);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    // Write formatted CDATA.
    rc = xmlTextWriterWriteFormatRaw(writer, "%d", ivalue);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatRaw\n");
        return;
    }

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXMLCValue:
 *
 */
void writeXMLCValue(xmlTextWriterPtr writer, const char* name, const char* cvalue) {
    int rc;


    // Write an element as child.
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST name, "%s", cvalue);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteFormatElement\n");
        return;
    }


}

/**
 * endXMLElement:
 *
 */
void endXMLElement(xmlTextWriterPtr writer) {
    int rc;

    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * endXMLDocument:
 *
 */
void endXMLDocument(xmlTextWriterPtr writer) {
    int rc;

    /* Here we could close remaining open elements using the
     * function xmlTextWriterEndElement, but since we do not want to
     * write any other elements, we simply call xmlTextWriterEndDocument,
     * which will do all the work. */
    rc = xmlTextWriterEndDocument(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndDocument\n");
        return;
    }

    xmlFreeTextWriter(writer);

}



/* ----------------------------------------------------------------------------------- */

/** Early-est extensions
 *
 * xmlns:ee="http://net.alomax/earlyest/xmlns/ee
 *
 */

/**
 * writeXML_EarlyEst_Discriminant:
 *
 */
void writeXML_EarlyEst_Discriminant(xmlTextWriterPtr writer, char* originID, char* name, double dvalue, double lower_uncertainty, double upper_uncertainty, int num_values) {

    // 20140530 AJL - do not write this element if no data
    if (num_values < 1) {
        return;
    }

    int rc;

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "ee:discriminant");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    writeXMLCValue(writer, "originID", originID);

    writeXMLRealQuantity(writer, "ee:disc", dvalue, DOUBLE_NULL, lower_uncertainty, upper_uncertainty, DOUBLE_NULL);

    writeXMLCValue(writer, "type", name);

    writeXMLIValue(writer, "stationCount", num_values);



    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}

/**
 * writeXML_EarlyEst_Discriminant:
 *
 */
void writeXML_EarlyEst_qualityIndicators(xmlTextWriterPtr writer, HypocenterDesc* phypo) {

    int rc;

    // Start an element.
    rc = xmlTextWriterStartElement(writer, BAD_CAST "ee:qualityIndicators");
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterStartElement\n");
        return;
    }
    /*
    // Add an attribute.
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterWriteAttribute\n");
        return;
    }
     */

    writeXMLDValue(writer, "ee:wtCountAssocWeight", phypo->qualityIndicators.wt_count_assoc_weight);
    writeXMLDValue(writer, "ee:errhWeight", phypo->qualityIndicators.errh_weight);
    writeXMLDValue(writer, "ee:errzWeight", phypo->qualityIndicators.errz_weight);
    //writeXMLDValue(writer, "ee:depthWeight", phypo->qualityIndicators.depth_weight);
    writeXMLDValue(writer, "ee:otVarianceWeight", phypo->qualityIndicators.ot_variance_weight);
    writeXMLDValue(writer, "ee:ampAttWeight", phypo->qualityIndicators.amp_att_weight);
    writeXMLDValue(writer, "ee:gapWeight", phypo->qualityIndicators.gap_weight);
    writeXMLDValue(writer, "ee:distanceCloseWeight", phypo->qualityIndicators.distanceClose_weight);
    writeXMLDValue(writer, "ee:distanceFarWeight", phypo->qualityIndicators.distanceFar_weight);
    writeXMLCValue(writer, "ee:qualityCode", phypo->qualityIndicators.quality_code);



    // Close the element.
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf("loc2xml: Error at xmlTextWriterEndElement\n");
        return;
    }

}





/* ----------------------------------------------------------------------------------- */

/**
 * this initialize the library and check potential ABI mismatches
 * between the version it was compiled for and the actual shared
 * library used.
 */
int
initLibXml(void) {

    LIBXML_TEST_VERSION

    return 0;

}

/**
 * Cleanup function for the XML library.
 */
int
cleanUpLibXml(void) {

    /*
     * Cleanup function for the XML library.
     */
    xmlCleanupParser();

    return 0;
}

/**
 * ConvertInput:
 * @in: string in a given encoding
 * @encoding: the encoding used
 *
 * Converts @in into UTF-8 for processing with libxml2 APIs
 *
 * Returns the converted UTF-8 string, or NULL in case of error.
 */
xmlChar *
ConvertInput(const char *in, const char *encoding) {
    xmlChar *out;
    int ret;
    int size;
    int out_size;
    int temp;
    xmlCharEncodingHandlerPtr handler;

    if (in == 0)
        return 0;

    handler = xmlFindCharEncodingHandler(encoding);

    if (!handler) {
        printf("ConvertInput: no encoding handler found for '%s'\n",
                encoding ? encoding : "");
        return 0;
    }

    size = (int) strlen(in) + 1;
    out_size = size * 2 - 1;
    out = (unsigned char *) xmlMalloc((size_t) out_size);

    if (out != 0) {
        temp = size - 1;
        ret = handler->input(out, &out_size, (const xmlChar *) in, &temp);
        if ((ret < 0) || (temp - size + 1)) {
            if (ret < 0) {
                printf("ConvertInput: conversion wasn't successful.\n");
            } else {
                printf
                        ("ConvertInput: conversion wasn't successful. converted: %i octets.\n",
                        temp);
            }

            xmlFree(out);
            out = 0;
        } else {
            out = (unsigned char *) xmlRealloc(out, out_size + 1);
            out[out_size] = 0; /*null terminating out */
        }
    } else {
        printf("ConvertInput: no mem\n");
    }

    return out;
}



#else

int main(void) {
    fprintf(stderr, "ERROR - libxml writer or output support not compiled in.\n");
    exit(1);
}
#endif
