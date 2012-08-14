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

int  printBranchLengths;

void do_NNI(tree * tr, nodeptr p, int swap)
{
  nodeptr       q;
  nodeptr       tmp;

  q = p->back;
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == 1)
   {
     tmp = p->next->back;
     hookup(p->next, q->next->back, q->next->z, tr->numBranches);
     hookup(q->next, tmp,           p->next->z, tr->numBranches);
   }
  else
   {
      tmp = p->next->next->back;
      hookup(p->next->next, q->next->back, q->next->z,       tr->numBranches);
      hookup(q->next,       tmp,           p->next->next->z, tr->numBranches);
   }
}


int main(int argc, char * argv[])
{
  tree        * tr;
  nodeptr       p;

  
  /* Do we want to print branch lengths? */
  printBranchLengths = FALSE;
  
  /* Set the minimum required info for the tree structure */
  tr = (tree *)malloc(sizeof(tree));
  tr->mxtips           = 5;
  tr->randomNumberSeed = 666;

  /* Setup some default values 
     TODO: The minimal initialization can be substantially smaller than what is
     described in axml.c 
  */
  setupTree(tr);

  /* Generate a random tree according to the seed given in tr->randomNumberSeed */
  makeRandomTree(tr);

  printf ( "Newick notation BEFORE NNI\n" );
  Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
  fprintf(stderr, "%s\n", tr->tree_string);

  
  p = tr->nodep[tr->mxtips + 1];
  /* perform the NNI move */
  do_NNI(tr, p, 1);

  printf ( "Newick notation AFTER NNI\n" );
  Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
  fprintf(stderr, "%s", tr->tree_string);

  

  return(EXIT_SUCCESS);
}


