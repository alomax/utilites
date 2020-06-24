/*
 * bin2etree4d - Load the output produced by Kim Olsen's finite difference
 *               code into a 4D etree
 *
 * Tiankai Tu
 * CSD, CMU
 *
 */


#include <stdio.h>

#include "etree.h"

typedef struct value_t {
    float Vx, Vy, Vz;
} value_t;


static char In3DFile[1024];

static char XVelocityFile[1024];
static char YVelocityFile[1024];
static char ZVelocityFile[1024];

static char etreepath[1024];
static int etreebufsize; /* in Megabyes */

/* 
 *
 * readarg - read in command line argument 
 *
 * return 0 if OK, -1 on error
 *
 */
int readarg(int argc, char **argv)
{
    if (argc != 7) 
        return -1;
    
    sscanf(argv[1], "%s", In3DFile);
    sscanf(argv[2], "%s", XVelocityFile);
    sscanf(argv[3], "%s", YVelocityFile);
    sscanf(argv[4], "%s", ZVelocityFile);
    sscanf(argv[5], "%s", etreepath);
    sscanf(argv[6], "%d", &etreebufsize);

    return 0;
}


/*
 *
 * usage()
 *
 */
void usage()
{
    printf("usage: bin2etree4d In3DFile XVelocityFile YVelocityFile ");
    printf("ZVelocityFile etreepath etreebufsize\n\n");

    printf("In3DFile: the ASCII file to specify the x, y, z, and t dimensions\n");
    printf("XVelocityFile: the binary file containing the X component\n");
    printf("YVelocityFile: the binary file containing the Y component\n");
    printf("ZVelocityFile: the binary file containing the Z component\n");
    printf("etreepath: the full path name of the etree to store the 4D dataset\n");
    printf("etreebufsize: the buffer size in terms of Megabytes to be used by the etree library\n");
    
    return;
}



int main(int argc, char ** argv)
{
    /* coordinate (should start from 1 ) */
    int i, j, k;
    int nbgz, nedz, nskpz;
    int nbgy, nedy, nskpy;
    int nbgx, nedx, nskpx;
    
    /* timestep */
    int nt;
    int timestep;
    int ntiskp;

    /* file descriptors to IN3D, SSX3D, SSY3D, SSZ3D */
    FILE *in3dfp;
    FILE *xfp, *yfp, *zfp;

    /* etree-related variables */
    int flags;
    etree_t *ep;
    etree_addr_t addr;


    /* read in command line arguments */
    if (readarg(argc, argv) != 0) {
        usage();
        exit(-1);
    }

    /*
     * 
     * load the range information into static variables 
     *
     */

    if ((in3dfp = fopen(In3DFile, "r")) == NULL) {
        perror("fopen In3D file");
        exit(-1);
    }
    /* Total timesteps */
    fscanf(in3dfp, "%d\n", &nt);
    fscanf(in3dfp, "%d\n", &ntiskp);

    /* X ranges */
    fscanf(in3dfp, "%d\n", &nbgx);
    fscanf(in3dfp, "%d\n", &nedx);
    fscanf(in3dfp, "%d\n", &nskpx);

    /* Y ranges */
    fscanf(in3dfp, "%d\n", &nbgy);
    fscanf(in3dfp, "%d\n", &nedy);
    fscanf(in3dfp, "%d\n", &nskpy);

    /* Z ranges */
    fscanf(in3dfp, "%d\n", &nbgz);
    fscanf(in3dfp, "%d\n", &nedz);
    fscanf(in3dfp, "%d\n", &nskpz);

    fclose(in3dfp);

    /* 
     *
     * open files (in binary format) containing the velocity along
     * X, Y, and X dimension 
     *
     */
    if (((xfp = fopen(XVelocityFile, "r")) == NULL) ||
        ((yfp = fopen(YVelocityFile, "r")) == NULL) ||
        ((zfp = fopen(ZVelocityFile, "r")) == NULL)) {
        perror("fopen binary velocity file");
        exit(-1);
    }
    
    /* 
     *
     * create a new etree to store the data points 
     *
     */
    flags = O_CREAT|O_TRUNC|O_RDWR;
    if ((ep = etree_open(etreepath, flags, etreebufsize, sizeof(value_t), 4)) 
        == NULL){
        fprintf(stderr, "Unable to create the output etree.\n");
        exit(-1);
    }
    
    /* 
     *
     * read each velocity (3D) from the three files for each grid point
     * and insert into an etree 
     *
     */

    /* these two fields of etree_addr_t never changes */
    addr.type = ETREE_LEAF;
    addr.level = ETREE_MAXLEVEL;

    for (timestep = 1; timestep <=nt; timestep += ntiskp) {

        /* same timestep value for all this iteration */
        addr.t = timestep;
        
        for (k = nbgz; k <= nedz; k+= nskpz) {

            /* k never changes inside this iteration */
            addr.z = k;

            for (j = nbgy; j <= nedy; j+= nskpy ) {
                
                /* j never changes inside this iteration */
                addr.y = j;
                
                for (i = nbgx; i <= nedx; i += nskpx) {
                    value_t value;

                    addr.x = i;

                    /* read binary input */
                    if ((fread(&value.Vx, sizeof(float), 1, xfp) != 1) ||
                        (fread(&value.Vy, sizeof(float), 1, yfp) != 1) ||
                        (fread(&value.Vz, sizeof(float), 1, zfp) != 1)) {

                        perror("read binary data");
                        exit(-1);
                    }

                    /* read formatted text input, for debug */
                    /*
                    if ((fscanf(xfp, "%f\n", &value.Vx) != 1) ||
                        (fscanf(yfp, "%f\n", &value.Vy) != 1) ||
                        (fscanf(zfp, "%f\n", &value.Vz) != 1)) {

                        perror("read ascii input");
                        exit(-1);
                    }
                    */

                    if (etree_insert(ep, addr, &value) != 0) {
                        etree_error_t err;

                        err = etree_errno(ep);
                        fprintf(stderr, "%s\n", etree_strerror(err));
                        exit(-1);
                    }
                } /* for x */
            } /* for y */
        } /* for z */
    } /* for timestep */     

    /* close the etree to flush it to disk. This is important! Otherwise, 
       data will get lost */
    etree_close(ep);

    /* close the three binary velocity files */
    fclose(xfp);
    fclose(yfp);
    fclose(zfp);

    return 0;
}
    


                    
                    

    
