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

int main (int argc, char * argv[])
{
  struct pllPhylip * phylip;
  double ** empiricalFrequencies;
  tree * tr;
  struct pllNewickTree * newick;
  partitionList * partitions;
  struct pllQueue * parts;

  if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

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

  /* TODO: do the base substitution */
  pllPhylipSubst (phylip, PHYLIP_DNA_DATA);

  /* Compute the empirical frequencies */
  empiricalFrequencies = pllBaseFrequenciesGTR (partitions, phylip);

  /* Create a PLL tree */
  tr = pllCreateTree (GAMMA, FALSE, FALSE, FALSE);

  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (tr, newick);
  /* Or instead of the previous function use the next commented line to create
     a random tree topology 
  pllTreeInitTopologyRandom (tr, phylip->nTaxa, phylip->label); */

  /* Connect the alignment with the tree structure */
  if (!pllTreeConnectAlignment (tr, phylip))
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


  /* Do some cleanup */
  pllPhylipDestroy (phylip);
  pllNewickParseDestroy (&newick);
  pllEmpiricalFrequenciesDestroy (&empiricalFrequencies, partitions->numberOfPartitions);

  pllPartitionsDestroy (&partitions, partitions->numberOfPartitions, tr->mxtips);
  pllTreeDestroy (tr);


  return (EXIT_SUCCESS);
}
