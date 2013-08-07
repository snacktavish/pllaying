#ifndef __pll_UTILS__
#define __pll_UTILS__

#include "axml.h"
//#include "parser/phylip/phylip.h"
#include "parser/alignment/alignment.h"
#include "parser/alignment/phylip.h"
#include "parser/newick/newick.h"
#include "parser/partition/part.h"
#include "queue.h"

linkageList* initLinkageList(int *linkList, partitionList *pr);
//void freeLinkageList( linkageList* ll);

int pllLinkAlphaParameters(char *string, partitionList *pr);
int pllLinkFrequencies(char *string, partitionList *pr);
int pllLinkRates(char *string, partitionList *pr);
int pllSetSubstitutionRateMatrixSymmetries(char *string, partitionList * pr, int model);

void pllSetFixedAlpha(double alpha, int model, partitionList * pr, pllInstance *tr);
void pllSetFixedBaseFrequencies(double *f, int length, int model, partitionList * pr, pllInstance *tr);
int  pllSetOptimizeBaseFrequencies(int model, partitionList * pr, pllInstance *tr);
void pllSetFixedSubstitutionMatrix(double *q, int length, int model, partitionList * pr,  pllInstance *tr);

nodeptr pllGetRandomSubtree(pllInstance *);
void makeParsimonyTree(pllInstance *tr);
void pllPartitionsDestroy (partitionList **, int);
int pllPartitionsValidate (struct pllQueue * parts, pllAlignmentData * alignmentData);
partitionList * pllPartitionsCommit (struct pllQueue * parts, pllAlignmentData * alignmentData);
void pllPhylipRemoveDuplicate (pllAlignmentData * alignmentData, partitionList * pl);
double ** pllBaseFrequenciesGTR (partitionList * pl, pllAlignmentData * alignmentData);
void pllTreeInitTopologyNewick (pllInstance * tr, struct pllNewickTree * nt, int bUseDefaultZ);
int pllLoadAlignment (pllInstance * tr, pllAlignmentData * alignmentData, partitionList *, int);
void pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models);
void pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList);
void pllTreeInitTopologyForAlignment (pllInstance * tr, pllAlignmentData * alignmentData);
void pllBaseSubstitute (pllAlignmentData * alignmentData, partitionList * partitions);
void  pllTreeDestroy (pllInstance * t);
pllInstance * pllCreateInstance (int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed);
int pllInitModel (pllInstance *, pllAlignmentData *, partitionList *);
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions);
int pllOptimizeModelParameters(pllInstance *tr, partitionList *pr, double likelihoodEpsilon);
#endif /* UTILS_H_ */
