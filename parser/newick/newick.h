#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "../../stack.h"

typedef struct
{
  int nodes;
  int tips;
  struct pllStack * tree;
} pllNewickTree;

struct item_t
{
  int depth;
  char * name;
  char * branch;
  int leaf;
  int rank;
};


#endif
