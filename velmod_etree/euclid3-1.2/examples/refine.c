/*
 * refine.c - Recursively build an octree
 */

#include <stdio.h>
#include <stdlib.h>
#include "etree.h"

typedef struct payload_t {
    int32_t val;
    char tag;
} payload_t;


void refine(etree_t *ep, etree_addr_t root);
int refine_pred(etree_addr_t addr, int payload);
int traverse(etree_t *ep); 

int main(int argc, char **argv)
{
    etree_t *ep;
    etree_addr_t root;
    payload_t payload;
    int count;
    char *filename;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s etree\n", argv[0]);
        exit(0);
    }
    filename = argv[1];

    /* Create an empty 3d etree with an integer record */
    ep = etree_open(filename, O_RDWR|O_CREAT|O_TRUNC, 1, sizeof(payload_t), 3); 
    if (ep == NULL) {
        fprintf(stderr, "Could not open etree file %s\n", filename);
        exit(1);
    }

    /* Register the schema for the etree to make it portable */
    if (etree_registerschema(ep, "int32_t val; char tag;") != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
        exit(1);
    }

    /* Initialize the tree with a single root octant */
    /* $begin refineinit */
    root.x = root.y = root.z = 0;
    root.level = ETREE_MAXLEVEL - 2;
    root.type = ETREE_LEAF;
    payload.tag = 'A';
    payload.val = 0;

    if (etree_insert(ep, root, &payload) != 0) {
        fprintf(stderr, "%s", etree_strerror(etree_errno(ep)));
        exit(1);
    }
    /* $end refineinit */

    /* Recursively refine the root octant */
    /* $begin callrefine */
    refine(ep, root);
    count = traverse(ep);
    fprintf(stderr, "The tree has %d octants\n", count);
    /* $end callrefine */

    /* Clean up */
    etree_close(ep);
    exit(0);
}

/*
 * refine - Depth-first z-order refinement of octant root
 */
/* $begin refine */
void refine(etree_t *ep, etree_addr_t root) 
{
    int i, j, k, incr;
    etree_addr_t child;
    static int val = 0;
    payload_t payload;

    /* $end refine */
    /* The root must be a leaf node */
    if (root.type != ETREE_LEAF) {
        fprintf(stderr, "Tried to refine an interior node\n");
        exit(1);
    }

    /* $begin refine */
    /* Check the terminating conditions for the recursion */
    if ((root.level >= ETREE_MAXLEVEL) || (!refine_pred(root, val)))
        return;

    /* Make the root an interior node */
    if (etree_delete(ep, root) != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
        exit(1);
    }
    root.type = ETREE_INTERIOR; 
    payload.val = val; 
    payload.tag = 'A' + (root.level - (ETREE_MAXLEVEL - 2));
    etree_insert(ep, root, &payload);

    /* Edge of octant at level (ETREE_MAXLEVEL-k) is 2^k ticks */
    child.level = root.level + 1;     
    child.type = ETREE_LEAF;
    incr = 1 << (ETREE_MAXLEVEL - child.level);

    /* Recursively expand the children in z-order */
    for (k = 0; k <= incr; k += incr) {
        child.z = root.z + k;
        for (j = 0; j <= incr; j += incr) {
            child.y = root.y + j;
            for (i = 0; i <= incr; i += incr) {
                child.x = root.x + i;
                payload.val = ++val;
                payload.tag = 'A' + (child.level - (ETREE_MAXLEVEL - 2));
                etree_insert(ep, child, &payload);

                refine(ep, child);
            }	
        }
    }
    return;
}
/* $end refine */

/*
 * refine_pred - Application specific refinement predicate
 */
/* $begin refinepred */
int refine_pred(etree_addr_t addr, int val) {
    if ((val == 0) || (val == 3))
	return 1;
    else
	return 0;
}
/* $end refinepred */
