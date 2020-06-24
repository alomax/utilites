/*
 * query - Example program that queries an octree
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


#include <stdio.h>
#include <stdlib.h>




/* include etree.h and other headers, code skipped */

/* $begin query */
#include "etree.h"

int main(int argc, char **argv)
{
    etree_t *ep;
    etree_addr_t addr, res_addr;
    int res_val;
    char buf[ETREE_MAXBUF], *metadata, *schema;

    /* $end query */
    if (argc != 2) {
        fprintf(stderr, "usage: %s etree\n", 
                argv[0]);
        exit(0);
    }
    /* $begin query */
    /* Open the etree for reading */
    ep = etree_open(argv[1], O_RDONLY, 0, 0, 3); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", argv[1]);  
        exit(1);
    }

    /* Print application meta data if it exists. */
    metadata = etree_getappmeta(ep);
    if (metadata != NULL) {
        fprintf(stderr, "Etree meta data: %s\n", metadata);  
        free(metadata);
    }

    /* Print schema if it exists. */
    schema = etree_getschema(ep);
    if (schema != NULL) {
        fprintf(stderr, "Etree schema: %s\n", schema);
        free(schema);
    }
    
    /* Query each octant specified on stdin */
    printf("Input  (x y z l): ");
    while (fgets(buf, ETREE_MAXBUF, stdin)) {
        sscanf(buf, "%u %u %u %d", &addr.x, &addr.y, &addr.z, &addr.level);
	/* $end query */
        /* Ignore comment or blank lines */
        if ((buf[0] == '#') || 
            (buf[0] == '\n') || 
            (buf[0] == '\t')) {
            continue;
        }
    /* $begin query */
        if (etree_search(ep, addr, &res_addr, "val", &res_val) != 0) {
            fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
            printf("Input  (x y z l): ");
            continue;
        }
        /* echo the search result and prompt for input, code skipped */
        /* $end query */

	printf("Output (x y z l): %s = %d\n\n", 
           etree_straddr(ep, buf, res_addr), res_val);
	printf("Input  (x y z l): ");
    /* $begin query */
    }

    etree_close(ep);
    return 0;
}
/* $end query */







