#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fasta.h"
#include "../../mem_alloc.h"



/* only check whether it is a valid alignment in fasta format */
static int
getFastaAlignmentInfo (int * inp, int * seqCount, int * seqLen)
{
  pllLexToken token;
  int input;

  input = *inp;

  *seqCount = *seqLen = 0;

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_STRING) return (0);

  while (1)
   {
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (1);

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          if (token.len < 2 || token.lexeme[0] != '>') return (0);
          break;
        default:
          return (0);
      }
     
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

     /* read second token (sequence) */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (0);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          if (!*seqLen)
            *seqLen = token.len;
          else
           {
             if (*seqLen != token.len) return (0);
           }
          break;
        default:
          return (0);
      }
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     ++ (*seqCount);
   }
}

static int
parseFastaAlignment (pllAlignmentData * alignmentData, int input)
{
  pllLexToken token;
  int i;

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_STRING) return (0);

  i = 1;
  while (1)
   {
     /* first parse the sequence label */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (1);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          alignmentData->sequenceLabels[i] = strndup (token.lexeme + 1, token.len - 1);
          break;
        default:
          return (0);
      }
     
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

     /* now parse the sequence itself */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (0);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          memmove (alignmentData->sequenceData[i], token.lexeme, token.len);
          break;
        default:
          return (0);
      }
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     ++ i;
   }
}


pllAlignmentData *
pllParseFASTA (const char * filename)
{
  int
    i,
    filesize,
    seqLen,
    seqCount,
    input;

  char * rawdata;
  pllAlignmentData * alignmentData;

  rawdata = pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (NULL);
   }

  lex_table_amend_fasta ();
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol ();


  if (!getFastaAlignmentInfo (&input, &seqCount, &seqLen))
   {
     fprintf (stderr, "Finished with error in parsing...\n");
     lex_table_restore ();
     rax_free (rawdata);
     return (NULL);
   }
  
  alignmentData = pllInitAlignmentData (seqCount, seqLen);
  
  printf ("\n---------------\n\n");

  init_lexan (rawdata, filesize);
  input = get_next_symbol ();

  if (!parseFastaAlignment (alignmentData, input))
   {
     printf ("Finished with error in parsing ...\n");
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
     rax_free(rawdata);
     return (NULL);
   }

  /* allocate alignment structure */


  lex_table_restore ();
  rax_free (rawdata);

  alignmentData->siteWeights = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i)
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}
