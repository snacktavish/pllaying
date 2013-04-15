#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "../../mem_alloc.h"
#include "../../stack.h"
#include "../../lexer.h"
#include "../../axml.h"

struct pllNewickTree
{
  int nodes;
  int tips;
  struct pllStack * tree;
};

struct item_t
{
  int depth;
  char * name;
  char * branch;
  int leaf;
  int rank;
};


struct pllNewickTree * pllNewickParseString (char * newick);
struct pllNewickTree * pllNewickParseFile (const char * filename);
int pllValidateNewick (struct pllNewickTree *);
tree * pllTreeCreateNewick (struct pllNewickTree * stack);
void pllNewickParseDestroy (struct pllNewickTree ** tree);
void  pllTreeDestroy (tree * t);
#endif
