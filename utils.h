#ifndef __pll_UTILS__
#define __pll_UTILS__

#include "axml.h"
#include "parser/phylip/phylip.h"
#include "parser/newick/newick.h"
#include "queue.h"


void read_msa(tree *tr, const char *filename);
void makeParsimonyTree(tree *tr);
void pllPartitionsDestroy (partitionList **, int, int);
int pllPartitionsValidate (struct pllQueue * parts, struct pllPhylip * phylip);
partitionList * pllPartitionsCommit (struct pllQueue * parts, struct pllPhylip * phylip);
void pllPhylipRemoveDuplicate (struct pllPhylip * phylip, partitionList * pl);
double ** pllBaseFrequenciesGTR (partitionList * pl, struct pllPhylip * phylip);
void pllTreeInitTopologyNewick (tree * tr, struct pllNewickTree * nt);
int pllTreeConnectAlignment (tree * tr, struct pllPhylip * phylip);
void pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models);
void pllTreeInitTopologyRandom (tree * tr, int tips, char ** nameList);
void pllBaseSubstitute (struct pllPhylip * phylip, partitionList * partitions);

#endif /* UTILS_H_ */
