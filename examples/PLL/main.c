#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../axml.h"
#include "../../utils.h"
#include "../../lexer.h"
#include "../../hash.h"

#include "../../parser/alignment/alignment.h"
#include "../../parser/alignment/phylip.h"

static nodeptr pickMyRandomSubtree(pllInstance *tr)
{
  nodeptr p;
  //do
  {
    /* select a random inner node */
    p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];

    /* select a random orientation */
    int exitDirection = rand() % 3;
    switch(exitDirection)
    {
      case 0:
        break;
      case 1:
        p = p->next;
        break;
      case 2:
        p = p->next->next;
        break;
      default:
        assert(0);
    }
  }
  //while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));
  return p;
}


int main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  pllInstance * tr;
  pllNewickTree * newick;
  partitionList * partitions;
  struct pllQueue * parts;
  pllListSPR * bestList;
  Wrapper * wrapper;
  wrapper = (Wrapper *) malloc (sizeof (Wrapper));
  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Create a PLL tree */
  tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

  /* Parse a PHYLIP file */
  alignmentData = pllParsePHYLIP (argv[1]);


  /* Parse a FASTA file */
  //alignmentData = pllParseFASTA (argv[1]);

  if (!alignmentData)
   {
     fprintf (stderr, "Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }

  /* Parse a NEWICK file */
  newick = pllNewickParseFile (argv[2]);
  if (!newick)
   {
     fprintf (stderr, "Error while parsing newick file %s\n", argv[2]);
     return (EXIT_FAILURE);
   }
  if (!pllValidateNewick (newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
   {
     fprintf (stderr, "Invalid phylogenetic tree\n");
     return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  parts = pllPartitionParse (argv[3]);
  
  /* Validate the partitions */
  if (!pllPartitionsValidate (parts, alignmentData))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (parts, alignmentData);
  
  /* destroy the  intermedia partition queue structure */
  pllQueuePartitionsDestroy (&parts);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllPhylipRemoveDuplicate (alignmentData, partitions);


  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
  printf ("!! ??  In our main file tr->nodep : %p\n\n\n\n\n", tr->nodep);

  pllDummy (wrapper);
  printf ("Address of wrapper->Item1 in main() %p\n", wrapper->Item1);
  printf ("Address of wrapper->Item2 in main() %p\n", wrapper->Item2);


  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels); */

//  /* Connect the alignment with the tree structure */
//  if (!pllLoadAlignment (tr, alignmentData, partitions, PLL_DEEP_COPY))
//   {
//     fprintf (stderr, "Incompatible tree/alignment combination\n");
//     return (EXIT_FAILURE);
//   }
//  
//  /* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
//  //pllInitModel(tr, TRUE, alignmentData, partitions);
//  printf ("!! thread id %d\n", tr->threadID);
//  printf ("\tAddress of tr before calling pllInitModel : %p \n", tr);
//  pllInitModel(tr, alignmentData, partitions);
//  printf ("\tAddress of tr after calling pllInitModel : %p \n", tr);
//  printf ("%p \n", tr);
//  printf ("!! thread id %d\n", tr->threadID);
//  
//  if (tr->start) 
//    printf ("After execution of pllInitModel : tr->start != NULL\n");
//  else
//    printf ("After execution of pllInitModel : tr->start == NULL\n");
//
//  /* TODO: evaluate likelihood, create interface calls */
//  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//  printf ("Likelihood: %f\n", tr->likelihood);
//  Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
//  printf ("Tree: %s\n", tr->tree_string);
//
//  /* another eval*/
//  double computed_lh = tr->likelihood;
//  evaluateGeneric (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
//  assert(computed_lh == tr->likelihood);
//
//  bestList = pllComputeSPR (tr, partitions, tr->nodep[tr->mxtips + 1], 1, 20, 20);
//
//  printf ("Number of SPRs: %d\n", bestList->entries);
//
//  pllDestroyListSPR (&bestList);
//
//  //printf ("Likelihood: %f\n", tr->likelihood);
//  
////  /* optimize BL */
////  {
////    double computed_lh = tr->likelihood;
////    treeEvaluate(tr, partitions, 32);
////    evaluateGeneric (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
////    assert(computed_lh < tr->likelihood);
////    printf ("Likelihood after BL opt: %f\n", tr->likelihood);
////  }
////  
////  /* do some simple SPR to improve your topology */
////  {
////    int i;
////    int max_radius = 15;
////    int min_radius = 1;
////    int num_iterations = 200;
////    tr->startLH = tr->endLH = tr->likelihood;
////    for(i=0; i<num_iterations; i++)
////    {
////      nodeptr p = pickMyRandomSubtree(tr);
////      /* make sure starting and end likelihood are the same */
////      tr->startLH = tr->endLH = tr->likelihood;
////      /* explore a neighbourhood of possible re-insertions */
////      rearrangeBIG(tr, partitions, p, min_radius, max_radius);
////      /* if one of the insertions was better, keep it as a best tree */
////      if(tr->startLH < tr->endLH)
////      {
////        restoreTreeFast(tr, partitions);
////        printf ("new Tree at iter %d: %s\n", i,  tr->tree_string);
////      }
////      /* show the tree we have right now (most of the time did not change)*/
////      Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
////      evaluateGeneric (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
////      if(i % (num_iterations/10) == 0)
////      {
////        modOpt(tr, partitions, 5.0);
////        printf("log lh: after %d iterations: %f \n",i, tr->likelihood);
////      }
////    }
////  }
////  /* Print resulting tree */
////  Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
////  printf ("Tree: %s\n", tr->tree_string);
////  Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
////  printf ("Tree: %s\n", tr->tree_string);
////  printf("Final log lh: %f \n", tr->likelihood);
//
//
//  /* Do some cleanup */
//  pllAlignmentDataDestroy (alignmentData);
//  pllNewickParseDestroy (&newick);
//
//  pllPartitionsDestroy (&partitions, tr->mxtips);
//  pllTreeDestroy (tr);
//

  return (EXIT_SUCCESS);
}
