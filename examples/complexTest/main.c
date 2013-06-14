#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../axml.h"
#include "../../utils.h"
#include "../../lexer.h"
#include "../../hash.h"
#include "../../parser/phylip/phylip.h"
#include "../../parser/newick/newick.h"
#include "../../parser/partition/part.h"
#include "../../globalVariables.h"

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
  struct pllPhylip * phylip, *phylip2;
  pllInstance * tr, *tr2;
  struct pllNewickTree * newick;
  partitionList * partitions, *partitions2;
  struct pllQueue * parts;
  int i;

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Create a PLL tree */
  tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
  tr2 = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

  /* Parse a PHYLIP file */
  phylip = pllPhylipParse (argv[1]);
  phylip2 = pllPhylipParse (argv[1]);

  if (!phylip)
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
  if (!pllPartitionsValidate (parts, phylip))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }

  /* commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (parts, phylip);
  partitions2 =  pllPartitionsCommit (parts, phylip2);

  /* destroy the  intermedia partition queue structure */
  pllQueuePartitionsDestroy (&parts);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllPhylipRemoveDuplicate (phylip, partitions);
  pllPhylipRemoveDuplicate (phylip2, partitions2);


  /* Set the topology of the PLL tree from a parsed newick tree */
  //pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
  /* Or instead of the previous function use the next commented line to create
     a random tree topology
  pllTreeInitTopologyRandom (tr, phylip->nTaxa, phylip->label); */

  pllTreeInitTopologyForAlignment(tr, phylip);

  /* Connect the alignment with the tree structure */
  if (!pllLoadAlignment (tr, phylip, partitions, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }

  /* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
 pllInitModel(tr, PLL_TRUE, phylip, partitions);

  /* TODO transform into pll functions !*/

  /*
     allocateParsimonyDataStructures(tr, partitions);
     makeParsimonyTreeFast(tr, partitions);
     freeParsimonyDataStructures(tr, partitions);
  */

  pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
  Tree2String (tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
  printf ("Tree: %s %d\n", tr->tree_string, tr->start->number);
  evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

  double
    firstTree = tr->likelihood;

  printf("%f \n", tr->likelihood);
  //computeBIGRAPID_Test(tr, partitions, PLL_TRUE);
  printf("final like %f\n", tr->likelihood);
  //pllInitModel(tr, PLL_TRUE, phylip, partitions);

  pllTreeInitTopologyNewick (tr2, newick, PLL_TRUE);
  if (!pllLoadAlignment (tr2, phylip2, partitions2, PLL_DEEP_COPY))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  pllInitModel(tr2, PLL_TRUE, phylip2, partitions2);

  Tree2String (tr2->tree_string, tr2, partitions2, tr2->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
  printf ("Tree: %s %d\n", tr2->tree_string, tr2->start->number);
  evaluateGeneric(tr2, partitions2, tr2->start, PLL_TRUE, PLL_FALSE);

  printf("%f \n", tr2->likelihood);

  double
    secondTree = tr2->likelihood;

  assert(firstTree == secondTree);

  pllOptimizeModelParameters(tr2, partitions2, 10.0);

  printf("%f \n", tr2->likelihood);

  pllPhylipDestroy (phylip);
  pllNewickParseDestroy (&newick);

  pllPartitionsDestroy (&partitions, partitions->numberOfPartitions, tr->mxtips);
  pllTreeDestroy (tr);

  pllPhylipDestroy (phylip2); 
  pllPartitionsDestroy (&partitions2, partitions2->numberOfPartitions, tr2->mxtips);
  pllTreeDestroy (tr2);



  for(i = 0; i < 4; i++)
    {
      FILE *f = fopen("dummy", "w");
      
      fprintf(f, "DNA, p1 = 1-200\n");
      fprintf(f, "DNA, p1 = 201-400\n");
      fprintf(f, "DNA, p1 = 401-705\n");
      
      fclose(f);
      
      tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
      
      phylip = pllPhylipParse (argv[1]);
      
      newick = pllNewickParseFile (argv[2]);
      
      parts = pllPartitionParse ("dummy");
      
      /* Validate the partitions */
      if (!pllPartitionsValidate (parts, phylip))
	{
	  fprintf (stderr, "Error: Partitions do not cover all sites\n");
	  return (EXIT_FAILURE);
	}
      
      /* commit the partitions and build a partitions structure */
      partitions = pllPartitionsCommit (parts, phylip);
      
      /* destroy the  intermedia partition queue structure */
      pllQueuePartitionsDestroy (&parts);
      
      /* eliminate duplicate sites from the alignment and update weights vector */
      pllPhylipRemoveDuplicate (phylip, partitions);
      
      pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
      if (!pllLoadAlignment (tr, phylip, partitions, PLL_DEEP_COPY))
	{
	  fprintf (stderr, "Incompatible tree/alignment combination\n");
	  return (EXIT_FAILURE);
	}
      pllInitModel(tr, PLL_TRUE, phylip, partitions);
      
      switch(i)
	{
	case 0:
	  //link params in one way 
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);	
	  break;
	case 1:
	  //link params in another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  break;
	case 2:
	  //link params in yet another way 
	  
	  pllLinkAlphaParameters("0,0,0", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,0", partitions);    	
	  break;

	case 3:
	  //also fiddle around with the Q matrices, make them to be non-GTR, but simpler
	  
	  pllLinkAlphaParameters("0,1,2", partitions);
	  pllLinkFrequencies("0,1,2", partitions);
	  pllLinkRates("0,1,2", partitions);    
	  
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 0);
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,5", partitions, 1);
	  pllSetSubstitutionRateMatrixSymmetries("0,1,2,3,4,0", partitions, 2);
	  break;	
	default:
	  assert(0);
	}
      
      evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
      printf("%f \n", tr->likelihood);
      pllOptimizeModelParameters(tr, partitions, 10.0);
      printf("%f \n", tr->likelihood); 
      //cleanup
      pllPhylipDestroy (phylip);
      pllNewickParseDestroy (&newick);
      
      pllPartitionsDestroy (&partitions, partitions->numberOfPartitions, tr->mxtips);
      pllTreeDestroy (tr);      
    }
  

  return (EXIT_SUCCESS);
}
