/*
 * packcvm - compactly store an unpacked CVM etree to an packed CVM etree
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "etree.h"
#include "domain.h"

/*
 * main 
 *
 */
int main(int argc, char ** argv)
{
    char *unpackedetree, *packedetree;
    etree_t *unpackedep, *packedep;
    char *appmeta;
    etree_error_t err;
    property_t property;
    etree_addr_t addr;

    /* read command line arguments */
    if (argc != 3) {
        fprintf(stderr, "Usage: pack unpackedetree packedetree\n");
        exit(-1);
    }
    unpackedetree = argv[1];
    packedetree = argv[2];

    /* open the unpacked etree and check the schema */
    if ((unpackedep = etree_open(unpackedetree, O_RDONLY, 0, 0, 0)) == NULL) {
        fprintf(stderr, "Cannot open unpacked etree.\n");
        exit(-1);
    }
    
    printf("Replicating schema\n");
    
    /* create the packed etree to store the database compactly */
    if ((packedep = etree_open(packedetree, O_CREAT|O_TRUNC|O_RDWR, 0, 0, 3))
        == NULL) {
        fprintf(stderr, "Cannot create packed etree.\n");
        exit(-1);
    }

    if (etree_registerschema(packedep, "float Vp; float Vs; float density;") 
        != 0) {
        fprintf(stderr, "packedetree: %s\n", 
                etree_strerror(etree_errno(packedep)));
        exit(-1);
    }

    printf("Storing data compactly to %s\n", packedetree);

    /* Get the initial cursor */
    addr.x = addr.y = addr.z = addr.level = 0;
    if (etree_initcursor(unpackedep, addr) != 0) {
        fprintf(stderr, "unpackedetree: %s\n", 
                etree_strerror(etree_errno(unpackedep)));
        exit(-1);
    }
    
    /* prepare for appending */
    if (etree_beginappend(packedep, 1) != 0) {
        fprintf(stderr, "packedetree: %s\n", 
                etree_strerror(etree_errno(packedep)));
        exit(-1);
    }

    do {
        if (etree_getcursor(unpackedep, &addr, "*", &property) != 0) {
            fprintf(stderr, "unpackedetree: %s\n", 
                    etree_strerror(etree_errno(unpackedep)));
            exit(-1);
        }
        
        if (etree_append(packedep, addr, &property) != 0) {
            fprintf(stderr, "packedetree: %s\n",
                    etree_strerror(etree_errno(packedep)));
            exit(-1);
        }
    } while (etree_advcursor(unpackedep) == 0);

    /* check termination condition */
    if ((err = etree_errno(unpackedep)) != ET_END_OF_TREE) {
        fprintf(stderr, "unpackedetree: %s\n", etree_strerror(err));
        exit(-1);
    }

    /* end append */
    etree_endappend(packedep);

    
    /* copy application metadata */
    printf("Replicate application metadata\n");
    if ((appmeta = etree_getappmeta(unpackedep)) != NULL) {
        if (etree_setappmeta(packedep, appmeta) != 0) {
            err = etree_errno(packedep);
            fprintf(stderr, "packedetree: %s\n", etree_strerror(err));
            exit(-1);
        }
    }

    
    etree_close(unpackedep);
    etree_close(packedep);

    return 0;
}	

