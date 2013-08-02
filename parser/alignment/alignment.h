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

char * __pllReadFile (const char *, int *);
void pllAlignmentDataDestroy (pllAlignmentData *);
void pllAlignmentDataDump (pllAlignmentData *);
pllAlignmentData * pllInitAlignmentData (int, int);

#endif
