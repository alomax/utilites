/*  -*- C -*-  */
/*
 * ======================================================================
 *
 *                           Brad T. Aagaard
 *                        U.S. Geological Survey
 *
 * {LicenseText}
 *
 * ======================================================================
 */

/* Application demonstrating how to do queries from C++ code.
 */

/* Switch used to select whether all values or selected ones are queried.
 */
 
 
/*

!!!! MODIFIED 2005.11.25 by A Lomax to dump P vel values to NLL grid buffer file.
*/


#define ALLVALS


#include "cencalvm/query/cvmquery.h"
#include "cencalvm/query/cvmerror.h"

#include <stdlib.h> /* USES exit() */
#include <unistd.h> /* USES getopt() */
#include <stdio.h> /* USES fopen(), fclose(), fprintf() */
#include <string.h> /* USES strcpy() */

#include <assert.h> /* USES assert() */

/* ------------------------------------------------------------------- */
/* Dump usage to stderr */
void
usage(void)
{ /* usage */
  fprintf(stderr,
	  "usage: cencalvmcquery [-h] -i fileIn -o fileOut -d dbfile\n"
	  "       [-l logfile] [-t queryType] [-r res] [-e dbextfile]\n"
	  "\n"
	  "  -i fileIn   File containing list of locations: 'lon lat elev'.\n"
	  "  -o fileOut  Output file with locations and material properties.\n"
	  "  -d dbfile   Etree database file to query.\n"
	  "  -e dbextfile  Etree extended database file to query.\n"
	  "  -t queryType  Type of query {'maxres', 'fixedres', 'waveres'}\n"
	  "  -r res        Resolution for query (not needed for maxres queries\n"
	  "  -h          Display usage and exit.\n"
	  "  -l logfile  Log file for messages.\n");
} /* usage */

/* ------------------------------------------------------------------- */
/* Parse command line arguments */
void
parseArgs(char* filenameIn,
	  char* filenameOut,
	  char* filenameDB,
	  char* filenameDBExt,
	  char* filenameLog,
	  char* queryType,
	  double* pQueryRes,
	  int argc,
	  char** argv)
{ // parseArgs
  assert(0 != filenameIn);
  assert(0 != filenameOut);
  assert(0 != filenameDB);
  assert(0 != filenameDBExt);
  assert(0 != filenameLog);
  assert(0 != queryType);
  assert(0 != pQueryRes);

  extern char* optarg;

  int nparsed = 1;
  int c = EOF;
  while ( (c = getopt(argc, argv, "hi:l:o:d:e:r:t:") ) != EOF) {
    switch (c)
      { /* switch */
      case 'i' : /* process -i option */
	strcpy(filenameIn, optarg);
	nparsed += 2;
	break;
      case 'o' : /* process -o option */
	strcpy(filenameOut, optarg);
	nparsed += 2;
	break;
      case 'd' : /* process -d option */
	strcpy(filenameDB, optarg);
	nparsed += 2;
	break;
      case 'e' : /* process -e option */
	strcpy(filenameDBExt, optarg);
	nparsed += 2;
	break;
      case 'h' : /* process -h option */
	nparsed += 1;
	usage();
	exit(0);
	break;
      case 'l' : /* process -l option */
	strcpy(filenameLog, optarg);
	nparsed += 2;
	break;
      case 't' : // process -t option
	strcpy(queryType, optarg);
	nparsed += 2;
	break;
      case 'r': // process -r option
	*pQueryRes = atof(optarg);
	nparsed += 2;
	break;
      default :
	usage();
      } /* switch */
  } /* while */
  if (nparsed != argc || 
      0 == strlen(filenameIn) ||
      0 == strlen(filenameOut) ||
      0 == strlen(filenameDB)) {
    usage();
    exit(1);
  } /* if */
} /* parseArgs */

/* ------------------------------------------------------------------- */
/* main */
int
main(int argc,
     char* argv[])
{ /* main */
  char filenameIn[256];
  char filenameOut[256];
  char filenameDB[256];
  char filenameDBExt[256];
  char filenameLog[256];
  char queryType[256];

  strcpy(filenameIn, "");
  strcpy(filenameOut, "");
  strcpy(filenameDB, "");
  strcpy(filenameDBExt, "");
  strcpy(filenameLog, "");
  strcpy(queryType, "maxres");
  double queryRes = 0.0;
  
  fprintf(stderr, "XXX!!!! Special version: MODIFIED 2005.11.25 by A Lomax to dump P vel values to NLL grid buffer file.\n");
  

  
  /* Parse command line arguments */
  parseArgs(filenameIn, filenameOut, filenameDB, filenameDBExt,
	    filenameLog, queryType, &queryRes,
	    argc, argv);

  /* Create query */
  void* query = cencalvm_createQuery();
  if (0 == query) {
    fprintf(stderr, "Could not create query.\n");
    return 1;
  } /* if */

  /* Get handle to error handler */
  void* errHandler = cencalvm_errorHandler(query);
  if (0 == errHandler) {
    fprintf(stderr, "Could not get handle to error handler.\n");
    return 1;
  } /* if */

  /* If log filename has been set, set log filename in error handler */
  if (strlen(filenameLog) > 0)
    cencalvm_error_logFilename(errHandler, filenameLog);

  /* Set database filename */
  if (0 != cencalvm_filename(query, filenameDB)) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return 1;
  } /* if */

  /* Set extended database filename if given */
  if (0 != strcmp("", filenameDBExt)) {
    if (0 != cencalvm_filenameExt(query, filenameDBExt)) {
      fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
      return 1;
    } /* if */
  } /* if */

  /* Set values to be returned in queries (or not) */
#if !defined(ALLVALS)
  int numVals = 2;
  char* pValNames[] = { "FaultBlock", "Zone" };
  if (0 != cencalvm_queryVals(query, pValNames, numVals)) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return 1;
  } /* if */
#else
  int numVals = 8;
#endif

  /* Open database for querying */
  if (0 != cencalvm_open(query)) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return 1;
  } /* if */
  
  /* Set query type and resolution */
  if (0 == strcasecmp(queryType, "maxres"))
    cencalvm_queryType(query, 0);
  else {
    if (queryRes < 0.0) {
      fprintf(stderr, "Query resolution must be a positive value.");
      usage();
      return 1;
    } // if
    cencalvm_queryRes(query, queryRes);
    if (0 == strcasecmp(queryType, "fixedres"))
      cencalvm_queryType(query, 1);
    else if (0 == strcasecmp(queryType, "waveres"))
      cencalvm_queryType(query, 2);
    else {
      fprintf(stderr, "Could not parse query string '%s' "
	      " into a known type of query.", queryType);
      usage();
      return 1;
    } // else
  } // else

  /* Open input file to read locations */
  FILE* fileIn = fopen(filenameIn, "r");
  if (0 == fileIn) {
    fprintf(stderr, "Could not open file '%s' to read query locations.\n",
	    filenameIn);
    return 1;
  } /* if */
  
  /* Open output file to accept data */
  FILE* fileOut = fopen(filenameOut, "w");
  if (0 == fileOut) {
    fprintf(stderr, "Could not open file '%s' to write query data.\n",
	    filenameIn);
    return 1;
  } /* if */

  /* Create array to hold values returned in queries */
  double* pVals = (double*) malloc(sizeof(double)*numVals);

  /* Read location from input file */
  double lon = 0.0;
  double lat = 0.0;
  double elev = 0.0;
  fscanf(fileIn, "%lf %lf %lf", &lon, &lat, &elev);

  float fval = 0.0;
  
  /* Continue operating on locations until end of file */
  while (!feof(fileIn)) {
    /* Query database */
    if (0 != cencalvm_query(query, &pVals, numVals, lon, lat, elev)) {
      fprintf(stderr, "%s", cencalvm_error_message(errHandler));
      /* If query generated an error, then bail out, otherwise reset status */
      if (2 == cencalvm_error_status(errHandler))
	return 1;
      cencalvm_error_resetStatus(errHandler);
    } /* if */

    /* Write values returned by query to output file */
    //AJL fprintf(fileOut, "%9.4f%8.4f%9.1f", lon, lat, elev);
#if !defined(ALLVALS)
    //AJL fprintf(fileOut, "%5d%5d\n", (int)pVals[0], (int)pVals[1]);
#else
    //AJL fprintf(fileOut, "%8.1f%8.1f%8.1f%9.1f%9.1f%9.1f%5d%5d\n",
	//AJL    pVals[0], pVals[1], pVals[2], pVals[3], pVals[4], pVals[5],
	//AJL    (int) pVals[6], (int) pVals[7]);
	fval = (float) (pVals[0] / 1000.0);
	// convert to slowness, 1km grid :
	fval = 1.0 / fval;
	fwrite(&fval, sizeof(float), 1, fileOut);
//fprintf(stderr, "%f ", fval);

#endif
      
    /* Read in next location from input file */
    fscanf(fileIn, "%lf %lf %lf", &lon, &lat, &elev);
  } /* while */

  /* Close database */
  cencalvm_close(query);

  /* Close input and output files */
  fclose(fileIn);
  fclose(fileOut);

  /* If an error was generated, write error message and bail out. */
  if (0 != cencalvm_error_status(errHandler)) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return 1;
  } /* if */

  /* Destroy query handle */
  cencalvm_destroyQuery(query);

  return 0;
} // main

// version
// $Id: cencalvmcquery.c,v 1.8 2005/10/29 18:52:20 brad Exp $

// End of file 
