#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylip.h"

int 
main (int argc, char * argv[])
{
  struct pllPhylip * phylip;

  if (argc != 2)
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }
  
  phylip = pllPhylipParse (argv[1]);
  if (!phylip) 
   {
     printf ("Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }

//  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  
  
  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  pllPhylipDump (phylip);
  pllPhylipRemoveDuplicate (phylip);
  
  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  int i;
  for (i = 0; i < phylip->seqLen; ++ i)
   {
     printf ("(%d,%d) ", i, phylip->weights[i]);
   }
  printf ("\n");
  pllPhylipDestroy (phylip);


  return (EXIT_SUCCESS);
}
