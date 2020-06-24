/* 
 * tiny - The first etree program 
 */


#include <stdio.h>
#include <stdlib.h>

#include "etree.h"

typedef struct _test_payload_t {
    int32_t val;
    char    tag;
} test_payload_t;

/* main */
int main(int argc, char **argv)
{
    etree_t*	 ep;
    etree_addr_t parent, child, res_addr;
    etree_tick_t edge_len, level;
    int i, j, k, len, count;
    test_payload_t parent_val, res_val;
    char buf[ETREE_MAXBUF];
  
    if (argc != 2) {
	fprintf(stderr, "usage: %s etree\n", argv[0]);
	exit(0);
    } /* if */
  
    /* Create an empty 3d etree where each record is an int */
    ep = etree_open (argv[1], O_RDWR|O_CREAT|O_TRUNC, 0, 
		     sizeof (test_payload_t), 3);
    if (ep == NULL) {
	fprintf(stderr, "Could not open %s\n", argv[1]);  
	exit(1);
    } /* if */

    /* Register schema for the etree to make it portable */
    const char* schema =   
	"int32_t val; "
	"char tag; ";

    if (etree_registerschema(ep, schema) != 0) {
	fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep))); 
	exit(1);
    } /* if */

    level    = 2;
    edge_len = 0x80000000 >> level;
    count    = 0;
    len	     = 1 << level; /* 2^level */

    for (k = 0; k < len; k++) {

	for (j = 0; j < len; j++) {

	    for (i = 0; i < len; i++) {

		/* Insert a parent octant that spans 1 / (len) of the domain */
		parent.level = level;
		parent.x     = i * edge_len;
		parent.y     = j * edge_len;
		parent.z     = k * edge_len;
		parent.type  = ETREE_LEAF;
		parent_val.val = count + 1;
		parent_val.tag = 'A' + count;

		if (etree_insert(ep, parent, &parent_val) != 0) {
		    fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep))); 
		    exit(1);
		} /* if */

		count++;
	    }
	}
    }

    printf ("created %d octants\n", count);


    /* Search for non-existent child of the octant we just inserted */
    child.x = child.y = child.z = 0; 
    child.level = ETREE_MAXLEVEL;

    if (etree_search(ep, child, &res_addr, "*", &res_val) != 0) {
	fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep))); 
	exit(1);
    } /* if */

    printf("Query : %s\n", etree_straddr(ep, buf, child));
    printf("Result: %s =\n", etree_straddr(ep, buf, res_addr));
    printf ("\tval=%d, tag=%c\n", res_val.val, res_val.tag);

    etree_close(ep);

    return 0;
} /* main */
