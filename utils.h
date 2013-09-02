#ifndef __pll_UTILS__
#define __pll_UTILS__

#include "axml.h"
//#include "parser/phylip/phylip.h"
#include "parser/alignment/alignment.h"
#include "parser/alignment/phylip.h"
#include "parser/newick/newick.h"
#include "parser/partition/part.h"
#include "queue.h"

typedef struct {
  double * zp;
  double * zpn;
  double * zpnn;
  double * zqr;
  nodeptr pn;
  nodeptr pnn;
  nodeptr r;
  nodeptr p;
  nodeptr q;
} sprInfoRollback;

typedef struct
 {
   nodeptr removeNode;
   nodeptr insertNode;
   double likelihood;
   double zqr[NUM_BRANCHES];
 } pllInfoSPR;

typedef struct
 {
   int max_entries;
   int entries;
   pllInfoSPR * sprInfo;
 } pllListSPR;

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
void pllPartitionsDestroy (pllInstance *, partitionList **);
int pllPartitionsValidate (struct pllQueue * parts, pllAlignmentData * alignmentData);
partitionList * pllPartitionsCommit (struct pllQueue * parts, pllAlignmentData * alignmentData);
void pllPhylipRemoveDuplicate (pllAlignmentData * alignmentData, partitionList * pl);
double ** pllBaseFrequenciesGTR (partitionList * pl, pllAlignmentData * alignmentData);
void pllTreeInitTopologyNewick (pllInstance *, pllNewickTree *, int);
int pllLoadAlignment (pllInstance * tr, pllAlignmentData * alignmentData, partitionList *, int);
void pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models);
void pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList);
void pllTreeInitTopologyForAlignment (pllInstance * tr, pllAlignmentData * alignmentData);
void pllBaseSubstitute (pllAlignmentData * alignmentData, partitionList * partitions);
void  pllDestroyInstance (pllInstance *);
pllInstance * pllCreateInstance (pllInstanceAttr *);
int pllInitModel (pllInstance *, partitionList *, pllAlignmentData *);
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions);
int pllOptimizeModelParameters(pllInstance *tr, partitionList *pr, double likelihoodEpsilon);

void pllInitListSPR (pllListSPR ** bestListSPR, int max);
void pllDestroyListSPR (pllListSPR ** bestListSPR);
pllListSPR * pllComputeSPR (pllInstance * tr, partitionList * pr, nodeptr p, int mintrav, int maxtrav, int max);
void pllCommitSPR (pllInstance * tr, partitionList * pr, pllInfoSPR * sprInfo, int saveRollbackInfo);
int pllRollbackSPR (pllInstance * tr, partitionList * pr);
void pllClearSprHistory (pllInstance * tr);
#endif /* UTILS_H_ */
