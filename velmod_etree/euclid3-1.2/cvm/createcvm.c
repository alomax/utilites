/*
 * createcvm.c: create an etree database by repeatedly sampling the 
 *              SCEC 3D community velocity
 * 
 * Tiankai Tu, David O'Hallaron
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/wait.h> 
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "etree.h"
#include "domain.h"

/* specificaition of sampling layers */
typedef struct layer_t {
    etree_tick_t ixmin, ixmax;
    etree_tick_t iymin, iymax;
    etree_tick_t izmin, izmax;
    int level;
}layer_t;

#define DIST1LAT 110922.0    /* distance(meters) of 1 degree latitude   */
#define DIST1LON 92382.0     /* distance(meters) of 1 degree longitude  */

/*same as ibig in in.h of SCEC model (in SCEC_cvm/ ) */
static const int SCEC_ibig = (1<<21); 

static char *domainspec, *layerspec, *destdir, *etreename;
static domain_t domain;
static double domain_rootsize;
static layer_t *layer;
static int samplelayers;


/*
 * Local functions 
 *
 *
 */
static void initlayers(const char *layerspec);
static void samplelayer(int index, etree_t *ep);
static void runFortran();



/*
 * main
 *
 * - obtain the query layer information from an external file specified on 
 *   the command line
 * - change directory to the TARGET dir, where the SCEC Fortran program
 *   should have been copied to.
 * - for each layer, assemble the btestin file, invoke the Fortran cvm 
 *   program to produce the btestout, then load the data into an (unpacked)
 *   cvm etree
 *
 */
int main(int argc, char **argv)
{
    etree_t *ep;
    char appmeta[1024];
    int i;

    /* read command line argument */
    if (argc != 5) {
        fprintf(stderr, "usage: createcvm domainspec layerspec destdir");
        fprintf(stderr, " etreename\n");
        exit(-1);
    }
    
    /* set the pathnames */
    domainspec = argv[1];
    layerspec = argv[2];
    destdir = argv[3];
    etreename = argv[4];

    /* read domain info */
    if (loaddomain(domainspec, &domain) != 0) {
        fprintf(stderr, "Cannot load domain specfication from %s.\n", 
                domainspec);
        exit(-1);
    }

    /* mapping between real-world address space and etree address space */
    domain_rootsize = (domain.xsize > domain.ysize) ? 
        domain.xsize : domain.ysize;
    domain_rootsize = (domain_rootsize > domain.zsize) ? 
        domain_rootsize : domain.zsize;

    /* read sample layer info */
    initlayers(layerspec);

    /* change to the dest directory */
    if (chdir(destdir) != 0) {
        perror("chdir to destdir");
        exit(-1);
    }
        
    /* create an unpacked etree */
    ep = etree_open(etreename, O_CREAT|O_TRUNC|O_RDWR, 0, 0, 3);
    if (ep == NULL) {
        fprintf(stderr, "Fail to create output (unpacked) etree\n");
        exit(-1);
    }
    if (etree_registerschema(ep, "float Vp; float Vs; float density;") != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
        exit(-1);
    }

    /* change to model dir */
    if (chdir("model") != 0) {
        perror("chdir to model directory");
        exit(-1);
    }

    fprintf(stderr, "\n\tQuerying SCEC cvm model (Version 3.0)\n\n");
    for (i = 0 ; i < samplelayers; i++) {
        fprintf(stderr, "Layer %d : ", i);
        fprintf(stderr, "%u %u %u %u %u %u %d\n",
                layer[i].ixmin, layer[i].ixmax, layer[i].iymin,
                layer[i].iymax, layer[i].izmin, layer[i].izmax,
                layer[i].level);

        samplelayer(i, ep);
    }
    fprintf(stderr, "\nFinish querying SCEC cvm model\n");
    
    /* releases the memory used by layer structure */
    free(layer);

    /* write application meta data to the etree before closing */
    sprintf(appmeta, 
            "LABASE latfrom %f lonfrom %f depfrom %f xsize %f ysize %f zsize %f",
            domain.latfrom, domain.lonfrom, domain.depfrom, 
            domain.xsize, domain.ysize, domain.zsize);

    if (etree_setappmeta(ep, appmeta) != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
        exit(-1);
    }

    if (etree_close(ep) != 0) {
        fprintf(stderr, "Error closing etree\n");
        exit(-1);
    }
        

    fprintf(stderr, "CVM etree (unpacked) stored in %s/%s\n", destdir,
            etreename);

    return 0;
}

/*
 * initlayers - read sample layer information from layerspec file 
 *
 * return if OK , exit -1 on error
 *
 */
void initlayers(const char *layerspec)
{
    int i;
    FILE *layerfp;

    /* load sampling layer info */
    if ((layerfp = fopen(layerspec, "r+")) == NULL) {
        perror("fopen layerfile");
        exit(-1);
    }

    if (fscanf(layerfp, "%d\n", &samplelayers) != 1) {
        perror("fscanf number of sample layers");
        exit(-1);
    }

    layer = (layer_t *)malloc(samplelayers * sizeof(layer_t));
    if (layer == NULL) {
        perror("malloc sample layers");
        exit(-1);
    }
    
    for (i = 0 ; i < samplelayers; i++) 
        if (fscanf(layerfp, "%u %u %u %u %u %u %d\n", 
                   &layer[i].ixmin, &layer[i].ixmax, &layer[i].iymin,
                   &layer[i].iymax, &layer[i].izmin, &layer[i].izmax,
                   &layer[i].level) != 7) {
            fprintf(stderr, "Fail to read in sample layer at line %d\n",
                    i + 1);
            exit(-1);
        }
    fclose(layerfp);
    return;
}


/*
 * samplelayer
 *
 * - obtain layer specification from static array
 * - create input file for CVM
 * - run CVM Fortran program in a child process
 * - read output of CVM and store into material database (etree)
 *
 */
void samplelayer(int index, etree_t *ep)
{
    int totalnum, i, j, k, level;
    FILE *infp, *outfp;
    double edgesize, lon, lat, dep;
    property_t property;
    etree_addr_t addr;
    etree_tick_t edgetics;
    
    /* establish relationship between the real-world coodinate system
       and etree address space */

    edgesize = domain_rootsize / ((etree_tick_t)1 << layer[index].level);
    
    level = layer[index].level;
    edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - level);

    /* determine the number of querying/samping points */
    totalnum = (layer[index].ixmax - layer[index].ixmin + 1) * 
        (layer[index].iymax - layer[index].iymin + 1) * 
        (layer[index].izmax - layer[index].izmin + 1);
    assert(totalnum <= SCEC_ibig);
    
    /* creaet input file to CVM */
    infp = fopen("btestin","w+");
    if (infp == NULL){
        perror("fopen btestin");
        exit(-1);
    }    

    fprintf(infp, "%d\n", totalnum);   /* needed by CVM */
    for (j = layer[index].iymin; j <= layer[index].iymax; j++) {
        lat = domain.latfrom + (j + 0.5) * edgesize / DIST1LAT;

        for (i = layer[index].ixmin; i <= layer[index].ixmax; i++) {
            lon = domain.lonfrom + (i + 0.5) * edgesize / DIST1LON;

            for (k = layer[index].izmin; k <= layer[index].izmax; k++) {
                /* sample the top surface center point */
                /* NOTE: for labase.2. The result database may contain
                   very small Vs values */
                dep = domain.depfrom + (k) * edgesize ;  
                fprintf(infp, "%.5f %.5f %.5f\n", lat, lon, dep);
            }
        }
    }
    fclose(infp);
    

    /* fork child process to run CVM Fortran program */
    runFortran();
    

    /* read output and store in etree */
    if ((outfp = fopen("btestout", "r+")) == NULL) {
        perror("fopen btestout");
        exit(-1);
    }
    
    addr.level = level;
    addr.type = ETREE_LEAF;

    for (j = layer[index].iymin; j <= layer[index].iymax; j++) {
        addr.y = edgetics * j;

        for (i = layer[index].ixmin; i <= layer[index].ixmax; i++) {
            addr.x = edgetics * i;

            for (k = layer[index].izmin; k <= layer[index].izmax; k++) {
                addr.z = edgetics * k;

                if (fscanf(outfp, "%*f %*f %*f %f %f %f\n", 
                           &property.Vp, &property.Vs, &property.density) != 3){
                    perror("read from btestout");
                    exit(-1);
                }

                if (etree_insert(ep, addr, &property) != 0) {
                    etree_error_t err = etree_errno(ep);
                    fprintf(stderr, "%s\n", etree_strerror(err));
                    exit(-1);
                }
            }
        }
    }
    fclose(outfp);

    return;
}
        
        

/*
 * runFortran - wrapper function to call cvm  Fortan program
 *
 * return if OK, exit -1 on error
 *
 */
void runFortran()
{
    pid_t pid;

    pid = fork();
    if (pid == -1) {
        perror("fork");
        exit(-1);
    } else if (pid == 0) {
        unlink("btestout");
        fprintf(stderr, "running SCEC cvm program...\n");
        execl("version3.0", "version3.0", (char *)0);
        perror("excel"); /* shall never get to here */
        exit(-1);
    } else {
        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status)) {
            fprintf(stderr, "done\n");
            return;
        }
        else {
            fprintf(stderr, "cvm exits abnormally\n");
            exit(-1);
        }
    }
}


