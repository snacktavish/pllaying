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
  tr->mxtips           = 20;
  tr->randomNumberSeed = 666;

  /* Setup some default values 
     TODO: The minimal initialization can be substantially smaller than what is
     described in axml.c 
  */
  setupTree(tr);

  /* Generate a random tree according to the seed given in tr->randomNumberSeed */
  makeRandomTree(tr);

  /* Print the tree */
  Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, 0, 0, 0, 0, SUMMARIZE_LH, 0,0);
  fprintf(stderr, "%s", tr->tree_string);



  return(EXIT_SUCCESS);
}


