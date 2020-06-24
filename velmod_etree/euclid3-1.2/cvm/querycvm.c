/*
 * querycvm - Query cvm etree for a cross section or volume
 *
 * Copyright (c) 2003 Tiankai Tu, Julio Lopez and David O'Hallaron
 * All rights reserved.  May not be used, modified, or copied 
 * without permission.
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
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "etree.h"
#include "domain.h"


#ifndef CVMBUFFER
#define CVMBUFFER 100
#endif

#define TOTALGRIDS ((etree_tick_t)1 << ETREE_MAXLEVEL)

static void usage();

static int 
query3d(etree_t *ep, const char *fieldname, etree_addr_t origin, 
        etree_tick_t spacing_ticks, int imax, int jmax, int kmax, FILE *outfp);

static int 
query2dnormal(etree_t *ep, const char *fieldname, etree_addr_t origin, 
              etree_tick_t spacing_ticks, int imax, int jmax, int kmax, 
              FILE *outfp);

static int 
query2dslant(etree_t *ep, const char *fieldname, etree_addr_t origin, 
             etree_tick_t spacing_ticks, int nmax, int kmax, double alpha, 
             FILE *outfp);


static int 
query1d(etree_t *ep, const char *fieldname, etree_addr_t origin,
        etree_tick_t spacing_ticks, int nmax, double cosalpha,
        double cosbeta, double cosgamma, FILE *outfp);


static int 
query0d(etree_t *ep, const char *fieldname, etree_addr_t origin, FILE *outfp);

/* performance measurement */
#ifdef PERF
static int queries;
static struct timeval start_t, end_t;
static double querytime, throughput;
#endif

int main(int argc, char **argv)
{
    char *cvmetree, *fieldname, *outfile;
    float p1_lon, p1_lat, p1_dep, p2_lon, p2_lat, p2_dep;
    int dim;
    float spacing;
    etree_tick_t spacing_ticks;
    double slant2d, slant3d;

    /* maximum index offset along X, Y, Z coordinate */
    int imax = 0, jmax = 0, kmax = 0;
    
    /* maximum index offset along the axis perpendicular to the Z
       axis, nmax and alpha are only useful when a vertical slice is
       specified and its normal is neither X or Y axis */
    int nmax = 0;
    double alpha = 0;
    
    /* direction vector of a line specified by (p1, p2)*/
    double cosalpha = 0, cosbeta = 0, cosgamma = 0;

    /* the lower-left corner and upper-right corner of the mininum bounding
       box of the query */
    float mbb_ll_lon, mbb_ll_lat, mbb_ll_dep;
    float mbb_ur_lon, mbb_ur_lat, mbb_ur_dep; 


    /* the origin of the domain of query */
    float origin_lon, origin_lat, origin_dep;
    etree_addr_t originaddr;

    etree_t *ep;
    domain_t domain;
    double cvm_rootsize;
    double cvm_ticksize;
    char *appmeta;

    FILE *outfp;

    /* read commandline arguments */
    if (argc != 12) {
        usage();
        exit(-1);
    }

    cvmetree = argv[1];
    if ((sscanf(argv[2], "%f", &p1_lon) != 1) ||
        (sscanf(argv[3], "%f", &p1_lat) != 1) ||
        (sscanf(argv[4], "%f", &p1_dep) != 1) ||
        (sscanf(argv[5], "%f", &p2_lon) != 1) ||
        (sscanf(argv[6], "%f", &p2_lat) != 1) ||
        (sscanf(argv[7], "%f", &p2_dep) != 1) ||
        (sscanf(argv[8], "%d", &dim) != 1)  ||
        (sscanf(argv[9], "%f", &spacing) != 1)) {
        usage();
        exit(-1);
    }
    fieldname = argv[10];
    if (strcmp(fieldname, "all") == 0) 
        strcpy(fieldname, "*");

    outfile = argv[11];

    /* set maximum index offset properly */
    switch (dim) {
    case(3) :
        imax = (int)(fabs((p1_lon - p2_lon) * DIST1LON) / spacing) + 1;
        jmax = (int)(fabs((p1_lat - p2_lat) * DIST1LAT) / spacing) + 1;
        kmax = (int)(fabs(p1_dep - p2_dep) / spacing) + 1;
        break;
    case(2):
        if (p1_dep == p2_dep) { 
            /* horizontal slice */
            imax = (int)(fabs((p1_lon - p2_lon) * DIST1LON) / spacing) + 1;
            jmax = (int)(fabs((p1_lat - p2_lat) * DIST1LAT) / spacing) + 1;
            kmax = 0;

        } else if (p1_lon == p2_lon) {
            /* vertical slice, with normal along X */
            imax = 0;
            jmax = (int)(fabs((p1_lat - p2_lat) * DIST1LAT) / spacing) + 1;
            kmax = (int)(fabs(p1_dep - p2_dep) / spacing) + 1;

        } else if (p1_lat == p2_lat) {
            /* vertical slice, with normal along Y */
            imax = (int)(fabs((p1_lon - p2_lon) * DIST1LON) / spacing) + 1;
            jmax = 0;
            kmax = (int)(fabs(p1_dep - p2_dep) / spacing) + 1;
        } else {
            /* vertical slice, with a slant normal, use nmax here */

            slant2d = sqrt(pow((p1_lon - p2_lon) * DIST1LON, 2) +
                           pow((p1_lat - p2_lat) * DIST1LAT, 2));

            nmax = (int)(slant2d / spacing) + 1;
            kmax = (int)(fabs(p1_dep - p2_dep) / spacing) + 1;
            
            alpha = atan(((p1_lat - p2_lat) * DIST1LAT) /
                         ((p1_lon - p2_lon) * DIST1LON));
        }

        break;

    case(1):
        slant3d = sqrt(pow((p1_lon - p2_lon) * DIST1LON, 2) +
                       pow((p1_lat - p2_lat) * DIST1LAT, 2) +
                       pow((p1_dep - p2_dep), 2));

        nmax = (int)(slant3d / spacing) + 1;
        cosalpha = ((p2_lon - p1_lon) * DIST1LON / slant3d);
        cosbeta = ((p2_lat - p1_lat) * DIST1LAT / slant3d);
        cosgamma = ((p2_dep - p1_dep) / slant3d);

    }
            

    /* set the minimum bounding box */
    if (p1_lon < p2_lon) {
        mbb_ll_lon = p1_lon;
        mbb_ur_lon = p2_lon;
    } else {
        mbb_ll_lon = p2_lon;
        mbb_ur_lon = p1_lon;
    }
 
    if (p1_lat < p2_lat) {
        mbb_ll_lat = p1_lat;
        mbb_ur_lat = p2_lat;
    } else {
        mbb_ll_lat = p2_lat;
        mbb_ur_lat = p1_lat;
    }

    if (p1_dep < p2_dep) {
        mbb_ll_dep = p1_dep;
        mbb_ur_dep = p2_dep;
    } else {
        mbb_ll_dep = p2_dep;
        mbb_ur_dep = p1_dep;
    }


    /* set the origin of query */
    if (dim == 3) {
        origin_lon = mbb_ll_lon;
        origin_lat = mbb_ll_lat;
        origin_dep = mbb_ll_dep;
    } else if (dim == 2) {
        if (alpha == 0) {
            /* not a vertical slant slice */
        
            origin_lon = mbb_ll_lon;
            origin_lat = mbb_ll_lat;
            origin_dep = mbb_ll_dep;
        }  else {

            if (p1_lon < p2_lon) {
                origin_lon = p1_lon;
                origin_lat = p1_lat;
            } else {
                origin_lon = p2_lon;
                origin_lat = p2_lat;
            }
            origin_dep = mbb_ll_dep;
        }
    } else {
        /* (dim == 1) || (dim == 0) */
        origin_lon = p1_lon;
        origin_lat = p1_lat;
        origin_dep = p1_dep;
    } 

    /* open the cvm etree for query */
    if ((ep = etree_open(cvmetree, O_RDONLY, CVMBUFFER, 0 , 0)) == NULL) {
        fprintf(stderr, "Cannot open CVM etree database %s\n", cvmetree);
        exit(-1);
    }

    /* obtain the material database application meta data */
    appmeta = etree_getappmeta(ep);
    if (appmeta == NULL) {
        fprintf(stderr, "CVM etree has no meta data\n");
        exit(-1);
    }
    
    if (sscanf(appmeta, 
               "LABASE latfrom %f lonfrom %f depfrom %f xsize %f ysize %f zsize %f", 
               &domain.latfrom, &domain.lonfrom, &domain.depfrom, 
               &domain.xsize, &domain.ysize, &domain.zsize) != 6) {
        fprintf(stderr, "Error load CVM etree meta data\n");
        exit(-1);
    }
    free(appmeta); /* remember to release memory of the text string */


    /* check the range of query */
    if ((mbb_ll_lon < domain.lonfrom) ||
        (mbb_ll_lat < domain.latfrom) ||
        (mbb_ll_dep < domain.depfrom) ||
        (mbb_ur_lon > (domain.lonfrom + domain.xsize / DIST1LON)) ||
        (mbb_ur_lat > (domain.latfrom + domain.ysize / DIST1LAT)) ||
        (mbb_ur_dep > (domain.depfrom + domain.zsize))) {
        fprintf(stderr, "Query domain out of the domain of CVM etree\n");
        exit(-1);
    }

    /* remember the real-world size of the CVM model root */
    cvm_rootsize = (domain.xsize > domain.ysize) ? domain.xsize : domain.ysize;
    cvm_rootsize = (cvm_rootsize > domain.zsize) ? cvm_rootsize : domain.zsize;
    cvm_ticksize = cvm_rootsize / TOTALGRIDS;

    /* set the address for the origin of query in the etree address space
       of the CVM model */
    originaddr.level = ETREE_MAXLEVEL;
    originaddr.type = ETREE_LEAF;
    originaddr.x = (origin_lon - domain.lonfrom) * DIST1LON / cvm_ticksize;
    originaddr.y = (origin_lat - domain.latfrom) * DIST1LAT / cvm_ticksize;
    originaddr.z = (origin_dep - domain.depfrom) / cvm_ticksize;

    /* ticks between grid point */
    spacing_ticks = spacing / cvm_ticksize;

    /* create/open ouput ascii file */
    if ((outfp = fopen(outfile, "w+")) == NULL) {
        perror("fopen");
        exit(-1);
    }

#ifdef PERF
    gettimeofday(&start_t, NULL);
#endif

    /* call proper routine to satisfy the query */
    if (dim == 3) {
#ifdef PERF
        queries = (imax ) * (jmax ) * (kmax );
#endif 
        if (query3d(ep, fieldname, originaddr, spacing_ticks, imax, jmax, 
                    kmax, outfp) != 0) {
            fprintf(stderr, "Query 3D region failed\n");
            exit(-1);
        } 
    } else if (dim == 2) {
        if (alpha == 0) {
#ifdef PERF
            queries = (imax ) * (jmax ) * (kmax );
#endif
            if (query2dnormal(ep, fieldname, originaddr, spacing_ticks, imax,
                              jmax, kmax, outfp) != 0) {
                fprintf(stderr, "Query 2D normal slice failed\n");
                exit(-1);
            }
        } else {
            /* a slant slice */
#ifdef PERF
            queries = (nmax ) * (kmax );
#endif
            if (query2dslant(ep, fieldname, originaddr, spacing_ticks,
                             nmax, kmax, alpha, outfp) != 0) {
                fprintf(stderr, "Query 2D slant slice failed\n");
                exit(-1);
            }
        }
    } else if (dim == 1) {
#ifdef PERF
        queries = nmax ;
#endif
        if (query1d(ep, fieldname, originaddr, spacing_ticks, nmax, 
                    cosalpha, cosbeta, cosgamma, outfp) != 0) {
            fprintf(stderr, "Query 1D line failed\n");
            exit(-1);
        }
    } else {
#ifdef PERF
        queries = 1;
#endif
        if (query0d(ep, fieldname, originaddr, outfp) != 0) {
            fprintf(stderr, "Query 0D point failed\n");
            exit(-1);
        }
    }

#ifdef PERF
    gettimeofday(&end_t, NULL);
#endif

    /* close output file */
    if (fclose(outfp) != 0) {
        perror("fclose");
        exit(-1);
    }

    /* close etree */
    if (etree_close(ep) != 0) {
        perror("etree_close");
        exit(-1);
    }

#ifdef PERF
    querytime = (end_t.tv_sec - start_t.tv_sec) * 1000.0 +
        (end_t.tv_usec - start_t.tv_usec) / 1000.0;
    throughput = (queries / querytime) * 1000;

    printf("%dD - Queries: %d Querytime(msec):%f Throughput(queries/sec):%f\n",
           dim, queries, querytime, throughput);
#endif    
           
    return 0;
}


/*
 * usage
 */
void usage()
{
    printf("usage: querycvm cvmetree p1_lon p1_lat p1_dep ");
    printf("p2_lon p2_lat p3_dep dim spacing fieldname outfile\n\n");
    
    printf("cvmetrree: pathname to CVM etree database\n");
    printf("p1_lon: the longitude (degree) of the first end point\n");
    printf("p1_lat: the latitude (degree) of the first end point\n");
    printf("p1_dep: the depth (meter) of the first end point\n");
    printf("p2_lon: the longitude (degree) of the second end point\n");
    printf("p2_lat: the latitude (degree) of the second end point\n");
    printf("p2_dep: the depth (meter) of the second end point\n");
    printf("dim: dimension of the query; 0, 1, 2 or 3\n");
    printf("spacing: spacing between grid point (meters)\n");
    printf("fieldname: name of the field to retrieve\n");
    printf("outfile: output file name\n");
    printf("\n");
    
    return;
}



/*
 * query3d: Service the query of a volume
 *
 * return 0 if OK, -1 on error
 *
 */
int query3d(etree_t *ep, const char *fieldname, 
          etree_addr_t origin, etree_tick_t spacing_ticks,
          int imax, int jmax, int kmax, 
          FILE *outfp)
{
    int i, j, k;
    etree_addr_t pointaddr;
    property_t property;
    float floatvalue;
    void *value;
    int IsSingleField, res;

    /* determine which field(s) should be retrieved */
    if (strcmp(fieldname, "*") == 0) {
        IsSingleField = 0;
        value = &property;
    } else {
        IsSingleField = 1;
        value = &floatvalue;
    }
        
    /* these two fields are never changed */
    pointaddr.level = ETREE_MAXLEVEL;
    pointaddr.type = ETREE_LEAF;

    pointaddr.z = origin.z;

    if (fprintf(outfp, "%d X %d X %d\n", imax , jmax , kmax ) != 3) {
        perror("fprintf indices");
        exit(-1);
    }


    for (k = 0; k < kmax; k++) {

        pointaddr.y = origin.y;
        for (j = 0; j < jmax; j++) {

            pointaddr.x = origin.x;
            for (i = 0; i < imax; i++) {
                
                /* we are not interest in the hit address, set it to NULL */
                res = etree_search(ep, pointaddr, NULL, fieldname, value);
                if (res != 0) {
                    fprintf(stderr, "query3d: %s\n", 
                            etree_strerror(etree_errno(ep)));
                    return -1;
                }
                
                /* write the result to the output file in ASCII format */
                if (fprintf(outfp, "%11d%11d%11d", i, j, k) != 3) {
                    fprintf(stderr, "Output index %d %d %d:", i, j, k);
                    perror("fprintf");
                    exit(-1);
                }

                if (IsSingleField) {
                    if (fprintf(outfp, "%9.2f\n", floatvalue) != 1) {
                        fprintf(stderr, "Output value for %d %d %d:",
                                i, j, k);
                        perror("fprintf");
                        exit(-1);
                    }
                }
                else {
                    if (fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, 
                                property.Vs, property.density) != 3) {
                        fprintf(stderr, "Output values for %d %d %d",
                                i, j, k);
                        perror("fprintf");
                        exit(-1);
                    }
                }

                /* update the X coordinate  */
                pointaddr.x += spacing_ticks;
            }
            
            /* update the Y coordinate  */
            pointaddr.y += spacing_ticks;
        }

        /* update the Z coordinate */
        pointaddr.z += spacing_ticks;
    }

    return 0;
}

                
/*
 * query2dnormal : serive the query of a 2D slice 
 *
 * return 0 if OK, -1 on error
 *
 */
int query2dnormal(etree_t *ep, const char *fieldname, 
                  etree_addr_t origin, etree_tick_t spacing_ticks,
                  int imax, int jmax, int kmax, 
                  FILE *outfp)
{
    int i, j, k;
    etree_addr_t pointaddr;
    property_t property;
    float floatvalue;
    void *value;
    int IsSingleField, res;

    /* determine which field(s) should be retrieved */
    if (strcmp(fieldname, "*") == 0) {
        IsSingleField = 0;
        value = &property;
    } else {
        IsSingleField = 1;
        value = &floatvalue;
    }
    
    pointaddr.level = ETREE_MAXLEVEL;
    pointaddr.type = ETREE_LEAF;
    
    if (kmax == 0) {
        /* horizontal slice, with Z as its normal */
        pointaddr.z = origin.z; 

        if (fprintf(outfp, "%d X %d\n", imax , jmax ) != 3) {
            perror("fprintf indices");
            exit(-1);
        }

        pointaddr.y = origin.y;
        for (j = 0; j < jmax; j++) {

            pointaddr.x = origin.x;
            for (i = 0; i < imax; i++) {
                
                /* we are not interest in the hit address, set it to NULL */
                res = etree_search(ep, pointaddr, NULL, fieldname, value);
                if (res != 0) {
                    fprintf(stderr, "query2d: %s\n", 
                            etree_strerror(etree_errno(ep)));
                    return -1;
                }
                
                /* write the result to the output file in ASCII format */
                if (fprintf(outfp, "%11d%11d", i, j) != 2) {
                    fprintf(stderr, "Ouput index %d %d: ", i, j);
                    perror("");
                    exit(-1);
                }
                if (IsSingleField) 
                    fprintf(outfp, "%9.2f\n", floatvalue);
                else 
                    fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, property.Vs,
                            property.density);

                /* update the X coordinate  */
                pointaddr.x += spacing_ticks;
            }
            
            /* update the Y coordinate  */
            pointaddr.y += spacing_ticks;
        }
    } else if (imax == 0) {
        /* Vertical slice, with X as its normal */
        pointaddr.x = origin.x;
        
        fprintf(outfp, "%d X %d\n", jmax , kmax );

        pointaddr.z = origin.z;
        for (k = 0; k < kmax; k++) {
            
            pointaddr.y = origin.y;
            for (j = 0; j < jmax; j++) {
                
                res = etree_search(ep, pointaddr, NULL, fieldname, value);
                if (res != 0) {
                    fprintf(stderr, "query2d: %s\n", 
                            etree_strerror(etree_errno(ep)));
                    return -1;
                }
                
                /* write the result to the output file in ASCII format */
                fprintf(outfp, "%11d%11d", j, k);
                if (IsSingleField) 
                    fprintf(outfp, "%9.2f\n", floatvalue);
                else 
                    fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, property.Vs,
                            property.density);

                /* update the Y coordinate  */
                pointaddr.y += spacing_ticks;
            }
            
            /* update the Z coordinate  */
            pointaddr.z += spacing_ticks;
        }
    } else if (jmax == 0) {
        /* Vertical slice, with Y as its normal */
        pointaddr.y = origin.y;

        fprintf(outfp, "%d X %d\n", imax , kmax );

        pointaddr.z = origin.z;
        for (k = 0; k < kmax; k++) {
            
            pointaddr.x = origin.x;
            for (i = 0; i < imax; i++) {
                
                res = etree_search(ep, pointaddr, NULL, fieldname, value);
                if (res != 0) {
                    fprintf(stderr, "query2d: %s\n", 
                            etree_strerror(etree_errno(ep)));
                    return -1;
                }
                
                /* write the result to the output file in ASCII format */
                fprintf(outfp, "%11d%11d", i, k);
                if (IsSingleField) 
                    fprintf(outfp, "%9.2f\n", floatvalue);
                else 
                    fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, property.Vs,
                            property.density);

                /* update the X coordinate  */
                pointaddr.x += spacing_ticks;
            }
            
            /* update the Z coordinate  */
            pointaddr.z += spacing_ticks;
        }
    } 

    return 0;
}


/*
 * query2dslant:
 *
 *
 */
int query2dslant(etree_t *ep, const char *fieldname, 
                 etree_addr_t origin, etree_tick_t spacing_ticks,
                 int nmax, int kmax, double alpha, FILE *outfp)
{
    int n, k;
    etree_addr_t pointaddr;
    property_t property;
    float floatvalue;
    void *value;
    int IsSingleField, res;
    etree_tick_t spacing_ticks_x, spacing_ticks_y;

    /* determine which field(s) should be retrieved */
    if (strcmp(fieldname, "*") == 0) {
        IsSingleField = 0;
        value = &property;
    } else {
        IsSingleField = 1;
        value = &floatvalue;
    }
    
    /* adjust the spacing_ticks along X and Y axis */
    spacing_ticks_x = (etree_tick_t)(spacing_ticks * cos(alpha));
    spacing_ticks_y = (etree_tick_t)(spacing_ticks * sin(alpha));

    /* these two fields never change */
    pointaddr.level = ETREE_MAXLEVEL;
    pointaddr.type = ETREE_LEAF;

    fprintf(outfp, "%d X %d\n", nmax , kmax );
    
    pointaddr.z = origin.z;
    for (k = 0; k < kmax; k++) {
        
        pointaddr.y = origin.y;
        pointaddr.x = origin.x;
        for (n = 0; n < nmax; n++) {

                res = etree_search(ep, pointaddr, NULL, fieldname, value);
                if (res != 0) {
                    fprintf(stderr, "query2d: %s\n", 
                            etree_strerror(etree_errno(ep)));
                    return -1;
                }
                
                /* write the result to the output file in ASCII format */
                fprintf(outfp, "%11d%11d", n, k);
                if (IsSingleField) 
                    fprintf(outfp, "%9.2f\n", floatvalue);
                else 
                    fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, 
                            property.Vs, property.density);

                /* update the X and Y coordinate  */
                pointaddr.x += spacing_ticks_x;
                pointaddr.y += spacing_ticks_y;
            }
        
        /* update Z coordinate */
        pointaddr.z += spacing_ticks;

    }
    return 0;
}


/*
 * query1d: query the grid point on a line whose direction vector is
 *          specified by cosalpha, cosbeta, cosgamma
 *
 * return 0 if OK, -1 on error
 *
 */
int query1d(etree_t *ep, const char *fieldname, etree_addr_t origin,
            etree_tick_t spacing_ticks, int nmax, double cosalpha,
            double cosbeta, double cosgamma, FILE *outfp)
{
    int n;
    etree_addr_t pointaddr;
    property_t property;
    float floatvalue;
    void *value;
    int IsSingleField, res;
    etree_tick_t spacing_ticks_x, spacing_ticks_y, spacing_ticks_z;

    /* determine which field(s) should be retrieved */
    if (strcmp(fieldname, "*") == 0) {
        IsSingleField = 0;
        value = &property;
    } else {
        IsSingleField = 1;
        value = &floatvalue;
    }

    /* adjust the spacing_ticks along X and Y axis */
    spacing_ticks_x = (etree_tick_t)(spacing_ticks * cosalpha);
    spacing_ticks_y = (etree_tick_t)(spacing_ticks * cosbeta);
    spacing_ticks_z = (etree_tick_t)(spacing_ticks * cosgamma);

    /* these two fields never change */
    pointaddr.level = ETREE_MAXLEVEL;
    pointaddr.type = ETREE_LEAF;
    
    pointaddr.x = origin.x;
    pointaddr.y = origin.y;
    pointaddr.z = origin.z;

    fprintf(outfp, "%d\n", nmax + 1);

    for (n = 0; n <= nmax; n++) {
        res = etree_search(ep, pointaddr, NULL, fieldname, value);
        if (res != 0) {
            fprintf(stderr, "query1d: %s\n", 
                    etree_strerror(etree_errno(ep)));
            return -1;
        }
                
        /* write the result to the output file in ASCII format */
        fprintf(outfp, "%11d", n);
        if (IsSingleField) 
            fprintf(outfp, "%9.2f\n", floatvalue);
        else 
            fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, 
                    property.Vs, property.density);
        
        /* update the X and Y coordinate  */
        pointaddr.x += spacing_ticks_x;
        pointaddr.y += spacing_ticks_y;
        pointaddr.z += spacing_ticks_z;
    }
        
    return 0;
}


/*
 * query0d: search for a particular point 
 *
 * return 0 if OK, -1 on error
 */
int query0d(etree_t *ep, const char *fieldname, etree_addr_t origin,
            FILE *outfp)
{
    etree_addr_t pointaddr;
    property_t property;
    float floatvalue;
    void *value;
    int IsSingleField, res;

    /* determine which field(s) should be retrieved */
    if (strcmp(fieldname, "*") == 0) {
        IsSingleField = 0;
        value = &property;
    } else {
        IsSingleField = 1;
        value = &floatvalue;
    }

    pointaddr = origin;

    res = etree_search(ep, pointaddr, NULL, fieldname, value);
    if (res != 0) {
        fprintf(stderr, "query0d: %s\n", etree_strerror(etree_errno(ep)));
        return -1;
    }

    if (IsSingleField)
        fprintf(outfp, "%9.2f\n", floatvalue);
        else 
            fprintf(outfp, "%9.2f%9.2f%9.2f\n", property.Vp, 
                    property.Vs, property.density);
    
    return 0;
}
