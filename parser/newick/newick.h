#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "../../mem_alloc.h"
#include "../../axml.h"
#include "../../stack.h"
#include "../../lexer.h"

struct pllNewickTree
{
  int nodes;
  int tips;
  struct pllStack * tree;
};

struct pllNewickTree * pllNewickParseString (char * newick);
struct pllNewickTree * pllNewickParseFile (const char * filename);
int pllValidateNewick (struct pllNewickTree *);
tree * pllTreeCreateNewick (struct pllNewickTree * stack);
void pllNewickParseDestroy (struct pllNewickTree ** tree);
void  pllTreeDestroy (tree * t);
#endif
