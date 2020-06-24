/*
 * domain.c - Routines for store and retrieve domain information
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

#include "stdio.h"

#include "domain.h"


/*
 * storedomain: store the domain specification in a readable text file
 *
 * return 0 if OK, -1 on error
 */
int storedomain(domain_t domain, const char *pathname)
{
    FILE *domainfp;

    if ((domainfp = fopen(pathname, "w+")) == NULL) {
        perror("fopen");
        return -1;
    }

    fprintf(domainfp, "Origin of the domain:\n");
    fprintf(domainfp, "Latitude (Degree): %f\n", domain.latfrom);
    fprintf(domainfp, "Longitude (Degree): %f\n", domain.lonfrom);
    fprintf(domainfp, "Depth (Meter): %f\n", domain.depfrom);
    fprintf(domainfp, "\n");

    fprintf(domainfp, "Size of the domain:\n");
    fprintf(domainfp, "X - along latitude (Meters): %f\n", domain.xsize);
    fprintf(domainfp, "Y - along longitude (Meters): %f\n", domain.ysize);
    fprintf(domainfp, "Z - along depth (Meters): %f\n", domain.zsize);
    fprintf(domainfp, "\n");
    
    
    if (fclose(domainfp) != 0) {
        perror("close");
        return -1;
    }
     
    return 0;
}


/*
 * loaddomain: Read the domain information from a spec file to a domain
 *             structure
 *
 * return 0 if OK, -1 on error
 *
 */
int loaddomain(const char *pathname, domain_t *domain)
{
    FILE *domainfp;

    /* load the domain info */
    if ((domainfp = fopen(pathname, "r")) == NULL) {
        perror("fopen");
        return -1;
    }

    fscanf(domainfp, "Origin of the domain:\n");
    if ((fscanf(domainfp, "Latitude (Degree): %f\n", &domain->latfrom) != 1) 
        ||
        (fscanf(domainfp, "Longitude (Degree): %f\n", &domain->lonfrom) != 1) 
        ||
        (fscanf(domainfp, "Depth (Meter): %f\n", &domain->depfrom) != 1)) {
        perror("fscanf origin");
        return -1;
    }

    fscanf(domainfp, "Size of the domain:\n");
    if ((fscanf(domainfp, "X - along latitude (Meters): %f\n", 
                &domain->xsize) != 1)
        ||
        (fscanf(domainfp, "Y - along longitude (Meters): %f\n", 
                &domain->ysize) != 1)
        ||
        (fscanf(domainfp, "Z - along depth (Meters): %f\n", 
                &domain->zsize) != 1)) {
        perror("fscanf sizes");
        return -1;
    }

    if (fclose(domainfp) != 0) {
        perror("fclose");
        return -1;
    }

    return 0;
}

