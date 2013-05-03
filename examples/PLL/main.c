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
  struct pllPhylip * phylip;
  double ** empiricalFrequencies;
  pllInstance * tr;
  struct pllNewickTree * newick;
  partitionList * partitions;
  struct pllQueue * parts;
  int i;

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Create a PLL tree */
  tr = pllCreateInstance (GAMMA, FALSE, FALSE, FALSE);

  /* Parse a PHYLIP file */
  phylip = pllPhylipParse (argv[1]);
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
  
  /* destroy the  intermedia partition queue structure */
  pllQueuePartitionsDestroy (&parts);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllPhylipRemoveDuplicate (phylip, partitions);

  /* Do the base substitution (from A,C,G....  ->   0,1,2,3....)*/
  pllBaseSubstitute (phylip, partitions);

  /* Compute the empirical frequencies */
  empiricalFrequencies = pllBaseFrequenciesGTR (partitions, phylip);

  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (tr, newick);
  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tr, phylip->nTaxa, phylip->label); */

  /* Connect the alignment with the tree structure */
  if (!pllLoadAlignment (tr, phylip))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  
  /* TODO: initialize partitions and models, create interface calls*/
  initializePartitions (tr, tr, partitions, 0, 0);              
  initModel (tr, empiricalFrequencies, partitions);
  resetBranches(tr);

  /* TODO: evaluate likelihood, create interface calls */
  evaluateGeneric (tr, partitions, tr->start, TRUE, FALSE);
  printf ("Likelihood: %f\n", tr->likelihood);
  Tree2String (tr->tree_string, tr, partitions, tr->start->back, TRUE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
  printf ("Tree: %s\n", tr->tree_string);

  /* another eval*/
  double computed_lh = tr->likelihood;
  evaluateGeneric (tr, partitions, tr->start, TRUE, FALSE);
  assert(computed_lh == tr->likelihood);
  //printf ("Likelihood: %f\n", tr->likelihood);
  
  /* optimize BL */
  {
    double computed_lh = tr->likelihood;
    treeEvaluate(tr, partitions, 32);
    evaluateGeneric (tr, partitions, tr->start, TRUE, FALSE);
    assert(computed_lh < tr->likelihood);
    printf ("Likelihood after BL opt: %f\n", tr->likelihood);
  }

  /* do some simple SPR to improve your topology */
  {
    int i;
    int max_radius = 15;
    int min_radius = 1;
    int num_iterations = 200;
    tr->startLH = tr->endLH = tr->likelihood;
    for(i=0; i<num_iterations; i++)
    {
      nodeptr p = pickMyRandomSubtree(tr);
      /* make sure starting and end likelihood are the same */
      tr->startLH = tr->endLH = tr->likelihood;
      /* explore a neighbourhood of possible re-insertions */
      rearrangeBIG(tr, partitions, p, min_radius, max_radius);
      /* if one of the insertions was better, keep it as a best tree */
      if(tr->startLH < tr->endLH)
      {
        restoreTreeFast(tr, partitions);
        printf ("new Tree at iter %d: %s\n", i,  tr->tree_string);
      }
      /* show the tree we have right now (most of the time did not change)*/
      Tree2String (tr->tree_string, tr, partitions, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
      evaluateGeneric (tr, partitions, tr->start, TRUE, FALSE);
      if(i % (num_iterations/10) == 0)
      {
        modOpt(tr, partitions, 5.0);
        printf("log lh: after %d iterations: %f \n",i, tr->likelihood);
      }
    }
  }
  /* Print resulting tree */
  Tree2String (tr->tree_string, tr, partitions, tr->start->back, TRUE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
  printf ("Tree: %s\n", tr->tree_string);
  Tree2String (tr->tree_string, tr, partitions, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);
  printf ("Tree: %s\n", tr->tree_string);
  printf("Final log lh: %f \n", tr->likelihood);


  /* Do some cleanup */
  pllPhylipDestroy (phylip);
  pllNewickParseDestroy (&newick);
  pllEmpiricalFrequenciesDestroy (&empiricalFrequencies, partitions->numberOfPartitions);

  pllPartitionsDestroy (&partitions, partitions->numberOfPartitions, tr->mxtips);
  pllTreeDestroy (tr);


  return (EXIT_SUCCESS);
}
