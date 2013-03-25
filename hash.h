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

#endif
