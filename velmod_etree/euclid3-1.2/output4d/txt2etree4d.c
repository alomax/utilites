/*
 * txt2etree4d - Convert a list of pixels to an etree
 * 
 * Each pixel octant address is a text line of the form:
 *
 * <x> <y> <z> <t> <value>
 * 
 * where <x>, <y>,  <z>, and <t>are unsigned integer coordinates;  the <level>
 * of all the pixels are at ETREE_MAXLEVEL, <type> is set to be ETREE_LEAF
 * for all of them. The value of the pixels are of integer type
 *
 */

/* $begin txt2etree4d */
#include <stdio.h>
#include "etree.h"

int main(int argc, char **argv)
{
    etree_t *ep;
    etree_addr_t addr;
    int value, count=0;
    char buf[ETREE_MAXBUF];

    /* $end txt2etree4d */
    if (argc != 2) {
	fprintf(stderr, "usage: %s etree < input4d.txt\n", argv[0]);
	exit(0);
    }
    /* $begin txt2etree4d */
    /* Create an empty 4d etree */
    ep = etree_open(argv[1], O_RDWR|O_CREAT|O_TRUNC, 1, sizeof(int), 4); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", argv[1]);
        exit(1);
    }

    /* Insert each octant into the tree */
    addr.type = ETREE_LEAF;
    addr.level = ETREE_MAXLEVEL;

    while (fgets(buf, ETREE_MAXBUF, stdin)) {
        sscanf(buf, "%u %u %u %u %d", 
	       &addr.x, &addr.y, &addr.z, &addr.t, &value);

        /* Ignore comment or blank lines */
        if ((buf[0] == '#') || 
            (buf[0] == '\n') || 
            (buf[0] == '\t')) {
            continue;
        }

        count++;
        if (etree_insert(ep, addr, &value) != 0) {
            fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
            exit(1);
        }
	printf("Loaded %s = %d\n", etree_straddr(ep, buf, addr), value);
    }

    printf("Loaded %d octants\n", count);
    etree_close(ep);
    exit(0);
}
/* $end txt2etree4d */


