#ifndef __mem_alloc_h
#define __mem_alloc_h
#include <stddef.h>


void *memalign(size_t align, size_t size);
void *malloc(size_t size);
void *realloc(void *p, size_t size);
void free(void *p);
int posix_memalign(void **p, size_t align, size_t size);
void *calloc(size_t n, size_t size);

void *malloc_aligned(size_t size);

#endif
