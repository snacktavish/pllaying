#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../../pll.h"

/** @file  alignment.c

    @brief Collection of routines for reading alignments

    Auxiliary functions for storing alignments read from predefined file formats
*/

/** @defgroup alignmentGroup Reading and parsing alignments
    
    This set of functions handles the reading and parsing of several file formats that describe multiple sequence alignments. They are also responsible for storing the alignment in an internal structure
*/

#ifdef __DEBUGGING_MODE
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
#endif

/** @ingroup alignmentGroup
    @brief Initialize alignment structure fields

    Allocates memory for the data structure that will hold the alignment and
    initializes it. It requires the number of sequences \a sequenceCount and
    the length of sequences \a sequenceLength. It returns a pointer to the
    initialized data structure.

    @param sequenceCount
      Number of sequences in the alignment
    
    @param sequenceLength
      Length of the sequences

    @param 
      Initialized alignment data structured
*/
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

/** @ingroup alignmentGroup
    @brief Deallocates the memory associated with the alignment data structure
    
    Deallocates the memory associated with the alignment data structure \a alignmentData.

    @param alignmentData
      The alignment data structure
*/
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


/** @ingroup alignmentGroup
    @brief Prints the alignment to the console

    @param alignmentData
      The alignment data structure
*/
static void 
printAlignmentData (pllAlignmentData * alignmentData)
 {
   int i;

   printf ("%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   for (i = 1; i <= alignmentData->sequenceCount; ++ i)
    {
      printf ("%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
    }
 }

/** @ingroup alignmentGroup
    @brief Dump the alignment to a file in PHYLIP sequential format

    Dumps the alignment contained in \a alignmentData to file \a filename.

    @note If \a filename exists, all contents will be erased

    @param alignmentData
      Alignment data structure

    @param filename
      Output filename

    @return
      Returns \b PLL_TRUE on success, otherwise \b PLL_FALSE.
*/
int
pllAlignmentDataDumpPHYLIP (pllAlignmentData * alignmentData, const char * filename)
{
  FILE * fp;
  int i;

  printAlignmentData (alignmentData);
  fp = fopen (filename,"w");
  if (!fp) return (PLL_FALSE);

  fprintf (fp, "%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     fprintf (fp, "%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
   }

  fclose (fp);
  return (PLL_TRUE);
}

/** @ingroup alignmentGroup
    @brief Dump the alignment to a file in FASTA format

    Dumps the alignment contained in \a alignmentData to file \a filename.

    @note If \a filename exists, all contents will be erased

    @param alignmentData
      Alignment data structure

    @param filename
      Output filename

    @return
      Returns \b PLL_TRUE on success, otherwise \b PLL_FALSE.
*/
int
pllAlignmentDataDumpFASTA (pllAlignmentData * alignmentData, const char * filename)
{
  FILE * fp;
  int i;

  printAlignmentData (alignmentData);
  fp = fopen (filename,"w");
  if (!fp) return (PLL_FALSE);

  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     fprintf (fp, ">%s\n%s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
   }

  fclose (fp);
  return (PLL_TRUE);
}
