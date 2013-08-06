#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"

int 
main (int argc, char * argv[])
{
  pllAlignmentData * alignmentData;
  int i;

  if (argc != 2)
   {
     return (EXIT_FAILURE);
   }
  
  alignmentData = pllParseFASTA (argv[1]);
  if (!alignmentData) 
   {
     printf ("Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }

//  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  
  printf ("!! %d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     printf ("%s\n%s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
   }
  
  pllAlignmentDataDestroy (alignmentData);


  return (EXIT_SUCCESS);
}
