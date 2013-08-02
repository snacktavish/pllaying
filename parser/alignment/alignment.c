#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "alignment.h"
#include "../../lexer.h"
//#include "ssort.h"

char * 
__pllReadFile (const char * filename, int * filesize)
{
  FILE * fp;
  char * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *filesize = ftell (fp);
  if (*filesize == -1) return (NULL);
  rewind (fp);

  /* allocate buffer and read file contents */
  rawdata = (char *) rax_malloc ((*filesize) * sizeof (char));
  if (rawdata) 
   {
     if (fread (rawdata, sizeof (char), *filesize, fp) != *filesize) 
      {
        rax_free (rawdata);
        rawdata = NULL;
      }
   }

  fclose (fp);

  return (rawdata);
}

static int
printTokens (int input)
{
  struct ltoken_t token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.class)
      {
        case LEX_NUMBER:
          printf ("LEX_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_STRING:
          printf ("LEX_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_EOF:
          printf ("LEX_EOF\n");
          break;
        case LEX_WHITESPACE:
          printf ("LEX_WHITESPACE\n");
          break;
        case LEX_NEWLINE:
          printf ("LEX_NEWLINE\n");
          break;
        case LEX_UNKNOWN:
          printf ("LEX_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
      }
     /* end of parser */


   }
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN);

  if (token.class == LEX_UNKNOWN) return (0);

  return (1);
}

//static struct pllPhylip *
//alloc_phylip_struct (int nTaxa, int seqLen)
pllAlignmentData *
pllInitAlignmentData (int sequenceCount, int sequenceLength)
 {
   int i;
   pllAlignmentData * alignmentData;
   void * mem;
   
   /** TODO */
   alignmentData               =  (pllAlignmentData *) rax_malloc (sizeof (pllAlignmentData));
   alignmentData->sequenceData = (unsigned char **) rax_malloc ((sequenceCount + 1) * sizeof (unsigned char *));
   mem = (void *) rax_malloc (sizeof (unsigned char) * (sequenceLength + 1) * sequenceCount);
   for (i = 1; i <= sequenceCount; ++i)
    {
      alignmentData->sequenceData[i]                 = (unsigned char *) (mem + (i - 1) * (sequenceLength + 1) * sizeof (unsigned char));
      alignmentData->sequenceData[i][sequenceLength] = 0;
    }
   alignmentData->sequenceData[0] = NULL;
    
   alignmentData->sequenceLabels = (char **) rax_calloc ((sequenceCount + 1), sizeof (char *));

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
