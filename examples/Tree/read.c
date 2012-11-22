#include <stdio.h>
#include <stdlib.h>
#include "axml.h"

void getStartingTree(tree *tr);
void makeRandomTree(tree *tr);

#ifdef __cplusplus
extern "C" {
boolean setupTree (tree *tr);
}
#else
boolean setupTree (tree *tr);
#endif

int main(int argc, char * argv[])
{
  tree        * tr;
  
  /* Do we want to print branch lengths? */
  int           printBranchLengths = FALSE;
  
  /* Get a starting tree from input */
  /*
  if (argc != 2)
   {
     fprintf(stderr, " usage: %s [TREE-FILE]\n", argv[0]);
     return(EXIT_FAILURE);
   }
  getStartingTree(NULL); */

  /* Set the minimum required info for the tree structure */
  tr = (tree *)malloc(sizeof(tree));
  tr->mxtips           = 6;
  tr->randomNumberSeed = 345;

  /* Setup some default values 
     TODO: The minimal initialization can be substantially smaller than what is
     described in axml.c 
  */
  setupTree(tr);

  /* Generate a random tree according to the seed given in tr->randomNumberSeed */
  makeRandomTree(tr);

  /* Print the tree */
  printTopology(tr, FALSE);
  printTopology(tr, TRUE);



  return(EXIT_SUCCESS);
}
