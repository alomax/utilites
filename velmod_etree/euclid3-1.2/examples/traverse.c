#include <stdio.h>
#include <stdlib.h>
#include "etree.h"
#include "schemax.h"


/**
 * eh_etree_straddr - Format a string representation of an etree address
 */
char* eh_etree_straddr(int dim, char* buf, const etree_addr_t* addr)
{
    if (dim == 3) { /* 3d case */
	sprintf(buf, "[0x%08X 0x%08X 0x%08X %02x-%c]",
		addr->x, addr->y, addr->z, addr->level,
		addr->type == ETREE_INTERIOR ? 'I' : 'L');
    }
  
    else { /* 4d case */
	sprintf(buf, "[0x%08X 0x%08X 0x%08X 0x%08X %02x-%c]",
		addr->x, addr->y, addr->z, addr->t, addr->level,
		addr->type == ETREE_INTERIOR ? 'I' : 'L');
    }

  return buf;
}


/* 
 * traverse - Walk the tree in z-order (pre-order)
 */
/* $begin traverse */
int traverse(etree_t *ep)
{
    etree_addr_t addr;
    void* payload;
    int count;
    char buf[ETREE_MAXBUF];
    etree_error_t err;
    const char* field_desc = NULL;
    char* schema;

    schema = etree_getschema(ep);

    if (schema != NULL) {
	field_desc = "*";
	free (schema);
    }

    /* Get the initial cursor */
    addr.x = addr.y = addr.z = addr.t = addr.level = 0;
    if (etree_initcursor(ep, addr) != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
        exit(-1);
    }

    payload = malloc (etree_getpayloadsize(ep));

    /* Interatively traverse the tree using the cursor mechanism */
    count = 0;
    do {
        if (etree_getcursor(ep, &addr, field_desc, payload) != 0) {
            fprintf(stderr, "%s\n", etree_strerror(etree_errno(ep)));
            exit(-1);
        } 
        count++;

        printf("Visited %s = ", eh_etree_straddr(ep->dimensions, buf, &addr));
	schemax_printpayload(ep, payload, stdout);
	putchar('\n');
    } while (etree_advcursor(ep) == 0);


    free (payload);

    if ((err = etree_errno(ep)) != ET_END_OF_TREE) {
        /* make sure that the cursor terminates correctly */
        fprintf(stderr, "%s\n", etree_strerror(err));
        exit(-1);
    }        

    return count;
}
/* $end traverse */
