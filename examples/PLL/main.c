#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../hash.h"
#include "pll.h"

int main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  pllInstance * tr;
  pllNewickTree * newick;
  partitionList * partitions;
  struct pllQueue * parts;
  int i;
  pllInstanceAttr attr;
  pllRearrangeList * bestList;

#ifdef _FINE_GRAIN_MPI
  pllInitMPI (&argc, &argv);
#endif

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Set the PLL instance attributes */
  attr.rateHetModel     = PLL_GAMMA;
  attr.fastScaling      = PLL_FALSE;
  attr.saveMemory       = PLL_FALSE;
  attr.useRecom         = PLL_FALSE;
  attr.randomNumberSeed = 12345;
  attr.numberOfThreads  = 8;            /* This only affects the pthreads version */

  /* Create a PLL tree */
  tr = pllCreateInstance (&attr);

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
  
  /* Initialize the model */
  pllInitModel(tr, partitions, alignmentData);


  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  

  bestList = pllCreateRearrangeList (20);


  tr->thoroughInsertion = 1;
  printf ("Computing the best 20 SPR and NNIs in a radius (1,20)     [thoroughInsertion = enabled]\n");
  pllRearrangeSearch (tr, partitions, PLL_REARRANGE_SPR, tr->nodep[tr->mxtips + 1], 1, 20, bestList);
  pllRearrangeSearch (tr, partitions, PLL_REARRANGE_NNI, tr->nodep[tr->mxtips + 1], 1, 20, bestList);

  printf ("Number of computed rearrangements: %d\n", bestList->entries);
  printf ("------------------------------------\n");

  for (i = 0; i < bestList->entries; ++ i)
   {
     printf ("%2d  Type: %s  Likelihood: %f\n", i, bestList->rearr[i].rearrangeType == PLL_REARRANGE_SPR ? "SPR" : "NNI", bestList->rearr[i].likelihood);
   }


  printf ("Committing bestList->rearr[0]                   [thoroughInsertion = disabled]\n");
  tr->thoroughInsertion = 0;
  pllRearrangeCommit(tr, partitions, &(bestList->rearr[0]), PLL_TRUE);

  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  tr->thoroughInsertion = PLL_FALSE;
  pllDestroyRearrangeList (&bestList);
  bestList = pllCreateRearrangeList (20);
  printf ("Computing the best 20 SPRs in a radius (1,30)     [thoroughInsertion = enabled]\n");
  pllRearrangeSearch (tr, partitions, PLL_REARRANGE_SPR, tr->nodep[tr->mxtips + 1], 1, 30, bestList);

  printf ("Number of SPRs computed : %d\n", bestList->entries);
  for (i = 0; i < bestList->entries; ++ i)
   {
     printf ("%2d  Type: %s  Likelihood: %f\n", i, bestList->rearr[i].rearrangeType == PLL_REARRANGE_SPR ? "SPR" : "NNI", bestList->rearr[i].likelihood);
   }

  printf ("Committing bestList->rearr[0]                   [thoroughInsertion = false]\n");
  tr->thoroughInsertion = 0;
  pllRearrangeCommit (tr, partitions, &(bestList->rearr[0]), PLL_TRUE);

  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  printf ("Rolling back...\n");
  pllRearrangeRollback (tr, partitions);
  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  printf ("Rolling back...\n");
  pllRearrangeRollback (tr, partitions);
  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
  printf ("New likelihood: %f\n\n", tr->likelihood);

  
  pllDestroyRearrangeList (&bestList);

  /* Do some cleanup */
  pllAlignmentDataDestroy (alignmentData);
  pllNewickParseDestroy (&newick);
  pllPartitionsDestroy (tr, &partitions);
  pllDestroyInstance (tr);
  
  return (EXIT_SUCCESS);
}
