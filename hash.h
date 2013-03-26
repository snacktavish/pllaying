#ifndef __pll_HASH__
#define __pll_HASH__

struct pllHashItem
{
  char * str;
  void * data;
  struct pllHashItem * next;
};

struct pllHashTable
{
  unsigned int size;
  struct pllHashItem ** Items;
};

unsigned int pllHashString (const char * s, unsigned int size);
int pllHashAdd  (struct pllHashTable * hTable, const char * s, void * item);
struct pllHashTable * pllHashInit (int n);
int pllHashSearch (struct pllHashTable * hTable, char * s, void ** item);
void pllHashDestroy (struct pllHashTable ** hTable);
#endif
