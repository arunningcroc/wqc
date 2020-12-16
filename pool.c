/* This file handles memory allocations in bulk.
 * The strategy is to simply allocate one large bulk of memory
 * and then never release it until the calculation is complete
 * This is a conscious tradeoff, due to several different reasons:
 *
 * 1. This reduces the number of machine code instructions in the resulting
 * compiled binary, thus reducing execution time
 * 2. It makes memory management easy - no worries about memory leaks
 * 3. It's simple to write a pool in this way, whereas a more complicated
 * pool could be a fairly nontrivial project.
 *
 * The tradeoff is excess memory use, but on a modern system it doesn't
 * really matter if save a couple of dozen megabytes.
 *
 * One of the things wrong with this is also that it returns unaligned memory,
 * resulting to some loss of performance. I can't be bothered to fiddle bits.
 */
#include "pool.h"
#include <stdio.h>
pool *
begin_calculation(size_t size)
{
    pool *p = malloc(size + sizeof(pool));
    p->next = (char *) &p[1];
    p->end = p->next + size;
    return p;
}
void
end_calculation(pool *p)
{
    free(p);
}

size_t
pool_available(pool *p)
{
    return p->end - p->next;
}

void *
pool_alloc(pool *p, size_t size)
{
    if(pool_available(p) < size){
    printf("Not enough memory allocated");
    return NULL;
    }
    void *mem = (void *) p->next;
    p->next += size;
    return mem;
}
