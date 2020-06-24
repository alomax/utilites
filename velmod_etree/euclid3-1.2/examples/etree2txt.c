/*
 * etree2txt - Convert an etree to a text file
 * 
 * Each line in the output text file has the form:
 *
 *     <x> <y> <z> <level> <type> <value.tag><value.val>
 * 
 * where <x>, <y>, and <z> are unsigned integer coordinates, <level>
 * is the octant level [0,..,31], <type> is either 0 (interior) or 1
 * (leaf), and <value> is the integer value of the octant.
 */

#include <stdio.h>
#include <stdlib.h>
#include "etree.h"

int traverse(etree_t *ep);

int main(int argc, char **argv)
{
    etree_t *ep;
    int count = 0;
    char *filename, *metadata, *schema;

    if (argc != 2) {
        fprintf(stderr, "usage: %s etree\n", argv[0]);
        exit(0);
    }
    filename = argv[1];

    /* Open the etree for reading */
    ep = etree_open(filename, O_RDONLY, 0, 0, 3); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", filename);
        exit(1);
    }
    
    /* Print application meta data if it exists. The metadata string is 
       allocated by the etree library. Application must release the 
       memory explicitly
    */
    metadata = etree_getappmeta(ep);
    if (metadata != NULL) {
        fprintf(stderr, "Etree meta data: %s\n", metadata);
        free(metadata);
    }

    /* Print schema if it exists. The schema string is allocated by the
       etree library. Application must release the memory explicitly
    */
    schema = etree_getschema(ep);
    if (schema != NULL) {
        fprintf(stderr, "Etree schema: %s\n", schema);
        free(schema);
    } else {
	fputs ("No schema defined\n", stderr);
    }


    count = traverse(ep);
    fprintf(stderr, "Dumped %d octants\n", count);

    etree_close(ep);
    exit(0);
}




