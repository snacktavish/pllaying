#ifndef __pll_UTILS__
#define __pll_UTILS__

#include "axml.h"
#include "parser/phylip/phylip.h"
#include "parser/newick/newick.h"
#include "queue.h"


//void read_msa(pllInstance *tr, const char *filename);
void makeParsimonyTree(pllInstance *tr);
void pllPartitionsDestroy (partitionList **, int, int);
int pllPartitionsValidate (struct pllQueue * parts, struct pllPhylip * phylip);
partitionList * pllPartitionsCommit (struct pllQueue * parts, struct pllPhylip * phylip);
void pllPhylipRemoveDuplicate (struct pllPhylip * phylip, partitionList * pl);
double ** pllBaseFrequenciesGTR (partitionList * pl, struct pllPhylip * phylip);
void pllTreeInitTopologyNewick (pllInstance * tr, struct pllNewickTree * nt, int bUseDefaultZ);
int pllLoadAlignment (pllInstance * tr, struct pllPhylip * phylip, partitionList *, int);
void pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models);
void pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList);
void pllTreeInitTopologyForAlignment (pllInstance * tr, struct pllPhylip * phylip);
void pllBaseSubstitute (struct pllPhylip * phylip, partitionList * partitions);
void  pllTreeDestroy (pllInstance * t);
pllInstance * pllCreateInstance (int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed);
void pllInitModel (pllInstance *, int, struct pllPhylip *, partitionList *);
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions);
#endif /* UTILS_H_ */
