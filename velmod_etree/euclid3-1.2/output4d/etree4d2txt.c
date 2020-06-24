/*
 * etree4d2txt - Convert an etree to a text file
 * 
 * Each line in the output text file has the form:
 *
 *     <x> <y> <z> <level> <type> <value>
 * 
 * where <x>, <y>, and <z> are unsigned integer coordinates, <level>
 * is the octant level [0,..,31], <type> is either 0 (interior) or 1
 * (leaf), and <value> is the integer value of the octant.
 */

#include <stdio.h>
#include "etree.h"

typedef struct value_t {
    float Vx, Vy, Vz;
} value_t;

int traverse(etree_t *ep);

int main(int argc, char **argv)
{
    etree_t *ep;
    etree_addr_t addr;
    etree_error_t err;
    int count=0;
    char *filename;
    value_t value;
    char buf[ETREE_MAXBUF];

    if (argc != 2) {
        fprintf(stderr, "usage: %s etree\n", argv[0]);
        exit(0);
    }
    filename = argv[1];

    /* Open the etree for reading */
    ep = etree_open(filename, O_RDONLY, 0, sizeof(value_t), 3); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", filename);
        exit(1);
    }
    
    /* traverse the 4D dataset */
    addr.x = addr.y = addr.z = addr.t = addr.level = 0;

    if (etree_initcursor(ep, addr) != 0) {
        err = etree_errno(ep);
        fprintf(stderr, "main: %s\n", etree_strerror(err));
        exit(-1);
    }

    count = 0;
    do {
        if (etree_getcursor(ep, &addr, &value) != 0) 
            break;
        count++;

        printf("Visited %s = %f %f %f\n", etree_straddr(ep, buf, addr), 
               value.Vx, value.Vy, value.Vz);
    } while (etree_advcursor(ep) == 0);

    if ((err = etree_errno(ep)) != ET_END_OF_TREE) {
        /* make sure that the cursor terminates correctly */
        fprintf(stderr, " %s\n", etree_strerror(err));
        exit(-1);
    }        

    etree_close(ep);
    exit(0);
}




