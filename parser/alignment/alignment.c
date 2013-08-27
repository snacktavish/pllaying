#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "alignment.h"
#include "../../lexer.h"
//#include "ssort.h"
#include "../../mem_alloc.h"

static int
printTokens (int input)
{
  pllLexToken token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.tokenType)
      {
        case PLL_TOKEN_NUMBER:
          printf ("PLL_TOKEN_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_STRING:
          printf ("PLL_TOKEN_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_EOF:
          printf ("PLL_TOKEN_EOF\n");
          break;
        case PLL_TOKEN_WHITESPACE:
          printf ("PLL_TOKEN_WHITESPACE\n");
          break;
        case PLL_TOKEN_NEWLINE:
          printf ("PLL_TOKEN_NEWLINE\n");
          break;
        case PLL_TOKEN_UNKNOWN:
          printf ("PLL_TOKEN_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        default:
          break;
      }
     /* end of parser */


   }
  while (token.tokenType != PLL_TOKEN_EOF && token.tokenType != PLL_TOKEN_UNKNOWN);

  if (token.tokenType == PLL_TOKEN_UNKNOWN) return (0);

  return (1);
}

pllAlignmentData *
pllInitAlignmentData (int sequenceCount, int sequenceLength)
 {
   int i;
   pllAlignmentData * alignmentData;
   void * mem;
   
   /** TODO */
   alignmentData               =  (pllAlignmentData *) malloc (sizeof (pllAlignmentData));
   alignmentData->sequenceData = (unsigned char **) malloc ((sequenceCount + 1) * sizeof (unsigned char *));
   mem = (void *) malloc (sizeof (unsigned char) * (sequenceLength + 1) * sequenceCount);
   for (i = 1; i <= sequenceCount; ++i)
    {
      alignmentData->sequenceData[i]                 = (unsigned char *) (mem + (i - 1) * (sequenceLength + 1) * sizeof (unsigned char));
      alignmentData->sequenceData[i][sequenceLength] = 0;
    }
   alignmentData->sequenceData[0] = NULL;
    
   alignmentData->sequenceLabels = (char **) calloc ((sequenceCount + 1), sizeof (char *));

   alignmentData->sequenceCount  = sequenceCount;
   alignmentData->sequenceLength = sequenceLength;
   /** TODO: remove siteWeights from alignment */
   alignmentData->siteWeights    = NULL;

   return (alignmentData);
 }

void
pllAlignmentDataDestroy (pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     rax_free (alignmentData->sequenceLabels[i]);
   }
  rax_free (alignmentData->sequenceLabels);
  rax_free (alignmentData->sequenceData[1]);
  rax_free (alignmentData->sequenceData);
  rax_free (alignmentData->siteWeights);
  rax_free (alignmentData);
}


static void 
printAlignmentData (pllAlignmentData * alignmentData)
 {
   int i;

   printf ("=> Dumping AlignmentData\n");
   printf ("%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   for (i = 1; i <= alignmentData->sequenceCount; ++ i)
    {
      printf ("|%s| |%s|\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
    }
 }

void
pllAlignmentDataDump (pllAlignmentData * alignmentData)
{
  FILE * fp;
  int i;

  printAlignmentData (alignmentData);
  fp = fopen ("output.seq.phy","w");
  fprintf (fp, "%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     fprintf (fp, "%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
   }

  fclose (fp);
}
