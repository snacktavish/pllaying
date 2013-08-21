#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "../../mem_alloc.h"
#include "../../stack.h"
#include "../../lexer.h"
#include "../../axml.h"

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


pllNewickTree * pllNewickParseString (char * newick);
pllNewickTree * pllNewickParseFile (const char * filename);
int pllValidateNewick (pllNewickTree *);
void pllNewickParseDestroy (pllNewickTree **);
char * pllReadFile (const char *, int *);
#endif
