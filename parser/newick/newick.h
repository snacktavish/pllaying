#ifndef NEWICK_H
#define NEWICK_H
#include "../../axml.h"

struct pllStack
{
  void * item;
  struct pllStack * next;
};

tree * pllTreeCreateNewick (struct pllStack * stack, int nodes, int tips);
int pllNewickParseFile (const char * filename, struct pllStack ** tree, int * nodes, int * leaves);
void pllTreeDestroy (tree * t);
void pllNewickParseDestroy (struct pllStack ** tree);
int pllNewickParseString (char * newick, struct pllStack ** tree, int * nodes, int * leaves);
int pllValidateNewick (struct pllStack * tree, int nodes, int leaves);

/* stack stuff */
void  pllStackClear (struct pllStack ** stack);
void * pllStackPop (struct pllStack ** head);
int pllStackPush (struct pllStack ** head, void * item);
int pllStackSize (struct pllStack ** stack);
#endif
