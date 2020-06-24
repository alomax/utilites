/*
 * definedomain.c - interactive interface for the user to specify the 
 *                  origin and size of the domain
 *                
 * Tiankai Tu
 * Computer Science Department
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * tutk@cs.cmu.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>

#include "domain.h"

#define MAXBUF 128

int main(int argc, char** argv)
{
    char buf[MAXBUF];
    char domainfile[MAXBUF];
    domain_t domain;


    printf("\n");
    printf(" This utility help you to define the domain of a CVM etree\n\n");

    printf(" Origin of the domain\n");
    printf(" Latitude (Degree): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.latfrom) != 1) {
        perror("read latitude of origin");
        exit(-1);
    }

    printf(" Longitude (Degree): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.lonfrom) != 1) {
        perror("read longitude of origin");
        exit(-1);
    }

    printf(" Depth (Meter): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.depfrom) != 1) {
        perror("read depth of origin");
        exit(-1);
    }
    
    printf("\n");
    printf(" Please specify the size of the domain\n");
    printf(" X - along latitude (Meters): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.xsize) != 1) {
        perror("read size along latitude");
        exit(-1);
    }

    printf(" Y - along longitude (Meters): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.ysize) != 1) {
        perror("read size along longitude");
        exit(-1);
    }

    printf(" Z - along depth (Meters): ");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%f", &domain.zsize) != 1) {
        perror("read size along depth");
        exit(-1);
    }
    
    /* store the meta data to the metadatafile */
    printf("\n");
    printf(" Pathname of this domain specification:");
    fgets(buf, MAXBUF, stdin);
    if (sscanf(buf, "%s", domainfile) != 1) {
        perror("read file name");
        exit(-1);
    }


    if (storedomain(domain, domainfile) != 0) {
        fprintf(stderr, "Cannot store domain spec in %s\n", domainfile);
        exit(-1);
    }

    printf("\n Specification stored in %s:\n\n", domainfile);
    return 0;
}
    

    
