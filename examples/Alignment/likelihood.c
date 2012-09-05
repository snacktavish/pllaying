#include <stdio.h>
#include <stdlib.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"



void read_msa (tree * tr, char * filename);

int main(int argc, char * argv[])
{

  tree        * tr;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s [binary-alignment-file]\n", argv[0]);
     return (1);
   }
  tr = (tree *)malloc(sizeof(tree));

  /* read the binary input, setup tree, initialize model with alignment */
  read_msa(tr,argv[1]);
  tr->randomNumberSeed = 665;
  makeRandomTree(tr);
  printf("Number of taxa: %d\n", tr->mxtips);
  printf("Number of partitions: %d\n", tr->NumberOfModels);


  /* compute the LH of the full tree */
  printf ("Virtual root: %d\n", tr->start->number);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood: %f\n", tr->likelihood);

  /* 8 rounds of branch length optimization */
  smoothTree(tr, 8);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood after branch length optimization: %f\n", tr->likelihood);
  
  return (0);
}
