#ifndef __pll_ALIGNMENT__
#define __pll_ALIGNMENT__

#include "../common.h"

typedef struct
 {
   int              sequenceCount;
   int              sequenceLength;
   char          ** sequenceLabels;
   unsigned char ** sequenceData;
   int            * siteWeights;
 } pllAlignmentData;

void pllAlignmentDataDestroy (pllAlignmentData *);
void pllAlignmentDataDump (pllAlignmentData *);
pllAlignmentData * pllInitAlignmentData (int, int);

#endif
