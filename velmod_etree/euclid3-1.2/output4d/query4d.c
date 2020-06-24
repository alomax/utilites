/*
 * query4d - Example program that queries an octree
 *  
 * Input a list of octant addresses on stdin and print the
 * corresponding record values on stdout. Each octant address is a
 * text line of the form:
 *
 * <x> <y> <z> <level> <type>
 * 
 * where <x>, <y>, and <z> are unsigned integer coordinates, <level>
 * is the octant level [0,..,31], <type> is either 0 (interior node)
 * or 1 (leaf node).
 */

/* $begin query */
#include <stdio.h>
#include "etree.h"

typedef struct value_t {
    float Vx, Vy, Vz;
} value_t;


int main(int argc, char **argv)
{
    etree_t *ep;
    etree_addr_t addr, res_addr;
    value_t res_val;
    char buf[ETREE_MAXBUF];

    /* $end query */
    if (argc != 2) {
        fprintf(stderr, "usage: %s etree\n", argv[0]);
        exit(0);
    }
    /* $begin query */
    /* Open the etree for reading */
    ep = etree_open(argv[1], O_RDONLY, 0, sizeof(value_t), 3); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", argv[1]);
        exit(1);
    }

    /* Query each octant specified on stdin */
    printf("Input  (x y z t): ");
    addr.level = ETREE_MAXLEVEL;

    while (fgets(buf, ETREE_MAXBUF, stdin)) {
        sscanf(buf, "%u %u %u %u", &addr.x, &addr.y, &addr.z, &addr.t);
	/* $end query */
        /* Ignore comment or blank lines */
        if ((buf[0] == '#') || 
            (buf[0] == '\n') || 
            (buf[0] == '\t')) {
            continue;
        }
	/* $begin query */
        if (etree_search(ep, addr, &res_addr, &res_val) != 0) {
            fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
            printf("Input  (x y z t): ");
	    continue;
        }

	printf("Output (x y z t): %s = %f %f %f\n\n", 
	       etree_straddr(ep, buf, res_addr), res_val.Vx, res_val.Vy,
           res_val.Vz);
	printf("Input  (x y z t): ");
    }

    etree_close(ep);
    exit(0);
}
/* $end query */







