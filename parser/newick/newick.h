#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "../../mem_alloc.h"
#include "../../axml.h"
#include "../../stack.h"
#include "../../lexer.h"

int pllNewickParseString (char * newick, struct pllStack ** tree, int * nodes, int * leaves);
int pllNewickParseFile (const char * filename, struct pllStack ** tree, int * nodes, int * leaves);
int pllValidateNewick (struct pllStack * tree, int nodes, int leaves);
tree * pllTreeCreateNewick (struct pllStack * stack, int nodes, int tips);
void pllNewickParseDestroy (struct pllStack ** tree);
void  pllTreeDestroy (tree * t);
#endif
