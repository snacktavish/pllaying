#include "mem_alloc.h"

#define PREFIX(X)   llalloc##X

void *PREFIX(memalign)(size_t align, size_t size);
void *PREFIX(malloc)(size_t size);
void *PREFIX(realloc)(void *p, size_t size);
void PREFIX(free)(void *p);


// forward to the actual malloc implementation


void *memalign(size_t align, size_t size) {
  return PREFIX(memalign)(align, size);
}

void *malloc( size_t size ) {
  return PREFIX(malloc)(size);
}
void *realloc( void *p, size_t size ) {
  return PREFIX(realloc)(p, size);
}


void free(void *p) {
  PREFIX(free)(p);
}

int posix_memalign(void **p, size_t align, size_t size) {
  return PREFIX(posix_memalign)(p, align, size);
}
void *calloc(size_t n, size_t size) {
  return PREFIX(calloc)(n,size);
}

void *malloc_aligned(size_t size) 
{
  const size_t BYTE_ALIGNMENT = 32;
  return memalign(BYTE_ALIGNMENT, size);
  
}
