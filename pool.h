#ifndef POOL
#define POOL

#include <stdlib.h>
typedef struct pool {
    char *next;
    char *end;
} pool;
pool *begin_calculation(size_t size);
void end_calculation(pool *p);
size_t pool_available(pool *p);
void *pool_alloc(pool *p, size_t size);
#endif
