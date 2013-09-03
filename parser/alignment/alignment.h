#ifndef __pll_ALIGNMENT__
#define __pll_ALIGNMENT__

typedef struct
 {
   int              sequenceCount;
   int              sequenceLength;
   char          ** sequenceLabels;
   unsigned char ** sequenceData;
   int            * siteWeights;
 } pllAlignmentData;

#endif
