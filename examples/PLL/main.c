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
  int i;

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


  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tr, alignmentData->sequenceCount, alignmentData->sequenceLabels); */

  /* Connect the alignment with the tree structure */
  if (!pllLoadAlignment (tr, alignmentData, partitions, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  
  /* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
  //pllInitModel(tr, TRUE, alignmentData, partitions);
  pllInitModel(tr, alignmentData, partitions);
  
  /* TODO: evaluate likelihood, create interface calls */
  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("Likelihood: %f\n\n", tr->likelihood);
  //Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
  //printf ("Tree: %s\n", tr->tree_string);

  /* another eval*/
  double computed_lh = tr->likelihood;
  evaluateGeneric (tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
  assert(computed_lh == tr->likelihood);
  int numBranches = partitions->perGeneBranchLengths ? partitions->numberOfPartitions : 1;

  tr->thoroughInsertion = 1;
  printf ("Computing the best 20 SPRs in a radius (1,20)     [thoroughInsertion = enabled]\n");
  bestList = pllComputeSPR (tr, partitions, tr->nodep[tr->mxtips + 1], 1, 20, 20);

  printf ("Number of SPRs computed : %d\n", bestList->entries);

  for (i = 0; i < bestList->entries; ++ i)
   {
     printf ("\t bestList->sprInfo[%2d].likelihood     = %f\n", i, bestList->sprInfo[i].likelihood);
   }


  printf ("Committing bestList->sprInfo[0]                   [thoroughInsertion = enabled]\n");
  tr->thoroughInsertion = 1;
  pllCommitSPR (tr, partitions, &(bestList->sprInfo[0]), PLL_TRUE);

  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  tr->thoroughInsertion = PLL_FALSE;
  pllDestroyListSPR (&bestList);
  printf ("Computing the best 20 SPRs in a radius (1,30)     [thoroughInsertion = enabled]\n");
  bestList = pllComputeSPR (tr, partitions, tr->nodep[tr->mxtips + 1], 1, 30, 20);

  printf ("Number of SPRs computed : %d\n", bestList->entries);
  for (i = 0; i < bestList->entries; ++ i)
   {
     printf ("\t bestList->sprInfo[2%d].likelihood     = %f\n", i, bestList->sprInfo[i].likelihood);
   }

  printf ("Committing bestList->sprInfo[0]                   [thoroughInsertion = false]\n");
  tr->thoroughInsertion = 0;
  pllCommitSPR (tr, partitions, &(bestList->sprInfo[0]), PLL_TRUE);

  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  printf ("Rolling back...\n");
  pllRollbackSPR (tr, numBranches);
  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  printf ("Rolling back...\n");
  pllRollbackSPR (tr, numBranches);
  evaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);


  pllDestroyListSPR (&bestList);

  /* Do some cleanup */
  pllAlignmentDataDestroy (alignmentData);
  pllNewickParseDestroy (&newick);

  pllPartitionsDestroy (&partitions, tr->mxtips);
  pllTreeDestroy (tr);


  return (EXIT_SUCCESS);
}
