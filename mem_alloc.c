#include "mem_alloc.h"

#include <stdio.h>

#define PREFIX(X)   llalloc##X

void *PREFIX(memalign)(size_t align, size_t size);
void *PREFIX(malloc)(size_t size);
void *PREFIX(realloc)(void *p, size_t size);
int PREFIX(posix_memalign)(void **p, size_t align, size_t size);
void *PREFIX(calloc)(size_t n, size_t size);
void PREFIX(free)(void *p);


// forward to the actual malloc implementation


void *rax_memalign(size_t align, size_t size) {
  return PREFIX(memalign)(align, size);
}

void *rax_malloc( size_t size ) {
  return PREFIX(malloc)(size);
}
void *rax_realloc( void *p, size_t size ) {
  return PREFIX(realloc)(p, size);
}


void rax_free(void *p) {
  PREFIX(free)(p);
}

int rax_posix_memalign(void **p, size_t align, size_t size) {
  return PREFIX(posix_memalign)(p, align, size);
}
void *rax_calloc(size_t n, size_t size) {
  return PREFIX(calloc)(n,size);
}

void *rax_malloc_aligned(size_t size) 
{
  const size_t BYTE_ALIGNMENT = 32;
  return rax_memalign(BYTE_ALIGNMENT, size);
  
}

#if 0 // this is impossible as the libc uses it internally
// make everyone suffer for using the standard allocator
void *memalign(size_t align, size_t size) {
//   fprintf( stderr, "using forbidden memalign\n" );
  abort();
  return 0;
}

void *malloc( size_t size ) {
//   fprintf( stderr, "using forbidden malloc\n" );
  abort();
  return 0;
}
void *realloc( void *p, size_t size ) {
//   fprintf( stderr, "using forbidden realloc\n" );
  abort();
  return 0;
}


void free(void *p) {
//   fprintf( stderr, "using forbidden free\n" );
  abort();
  
}

int posix_memalign(void **p, size_t align, size_t size) {
//   fprintf( stderr, "using forbidden posix_memalign\n" );
  abort();
  return 0;
}
void *calloc(size_t n, size_t size) {
//   fprintf( stderr, "using forbidden calloc\n" );
  abort();
  return 0;
}

void *malloc_aligned(size_t size) 
{
//   fprintf( stderr, "using forbidden malloc_align\n" );
  abort();
  return 0;
  
}
#endif