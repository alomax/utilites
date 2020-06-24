/*
 * domain.h - data types and basic routines for storeing and retrieve
 *            domain information
 *
 * Copyright (c) 2003 Tiankai Tu, David R. O'Hallaron
 * All rights reserved.  May not be used, modified, or copied 
 * without permission.
 *
 * Contact:
 * Tiankai Tu
 * Computer Science Department
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * tutk@cs.cmu.edu
 *
 */

#ifndef DOMAIN_H
#define DOMAIN_H

#define DIST1LAT 110922.0    /* distance(meters) of 1 degree latitude   */
#define DIST1LON 92382.0     /* distance(meters) of 1 degree longitude  */


/*
 * domain_t: Domain specification 
 *
 */
typedef struct domain_t {
    float lonfrom;           /* degree                                  */
    float latfrom;           /* degree                                  */
    float depfrom;           /* meter                                   */

    float xsize;             /* along latitude, from west to east       */
    float ysize;             /* along longitude, from south to north    */
    float zsize;             /* along depth, from surface to kernel     */
} domain_t;



/*
 * property_t: Standard payload for CVM model 
 *
 */
typedef struct property_t {
    float Vp;
    float Vs;
    float density;
} property_t;


int storedomain(domain_t domain, const char *pathname);
int loaddomain(const char *pathname, domain_t * domain);

#endif /* DOMAIN_H */
