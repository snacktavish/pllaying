#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylip.h"
#include "../../mem_alloc.h"
static int
read_phylip_header (int * inp, int * sequenceCount, int * sequenceLength)
{
  pllLexToken token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceCount = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceLength = atoi (token.lexeme);

  *inp = input;

  return (*sequenceCount && *sequenceLength);
}

static inline int
parsedOk (int * actLen, int sequenceCount, int sequenceLength)
{
  int i;

  for (i = 1; i <= sequenceCount; ++ i)
   {
     if (actLen[i] != sequenceLength) return (0);
   }
  
  return (1);
}


static int
parse_phylip (pllAlignmentData * alignmentData, int input)
{
  int i,j;
  pllLexToken token;
  int * sequenceLength;
  int rc;

  sequenceLength = (int *) rax_calloc (alignmentData->sequenceCount + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % alignmentData->sequenceCount;
    if (i < alignmentData->sequenceCount) 
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
          rax_free (sequenceLength);
          return (0);
        }

       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)


       if (token.tokenType != PLL_TOKEN_STRING && token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_FLOAT)
        {
          rax_free (sequenceLength);
          return (0);
        }
       alignmentData->sequenceLabels[i + 1] = strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     }
    
    while (1)
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
         rax_free (sequenceLength);
         return (0);
        }
       
       if (token.tokenType == PLL_TOKEN_NEWLINE) break;

       if (token.tokenType != PLL_TOKEN_STRING)
        {
          rax_free (sequenceLength);
          return (0);
        }

       if (sequenceLength[j + 1] + token.len > alignmentData->sequenceLength) 
        {
          fprintf (stderr, "Sequence %d is larger than specified\n", j + 1);
          rax_free (sequenceLength);
          return (0);
        }
       memmove (alignmentData->sequenceData[j + 1] + sequenceLength[j + 1], token.lexeme, token.len);
       sequenceLength[j + 1] += token.len;

       NEXT_TOKEN
       CONSUME (PLL_TOKEN_WHITESPACE)
     }
    CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE);
  }
}

/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL
 *     ENDL -> PLL_TOKEN_WHITESPACE PLL_TOKEN_NEWLINE | PLL_TOKEN_NEWLINE
 *     DATA -> PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA | 
 *             PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF
 */
//struct pllPhylip *
//pllPhylipParse (const char * filename)
pllAlignmentData *
pllParsePHYLIP (const char * filename)
{
  int 
    i, filesize, input, sequenceCount, sequenceLength;
  char * rawdata;
  pllAlignmentData * alignmentData;

  rawdata = pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (NULL);
   }
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (&input, &sequenceCount, &sequenceLength))
   {
     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     return (NULL);
   }

  lex_table_amend_phylip();

  /* allocate alignment structure */
  alignmentData = pllInitAlignmentData (sequenceCount, sequenceLength);

  if (! parse_phylip (alignmentData, input))
   {
     printf ("Finished with error in parsing ...\n");
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
     rax_free (rawdata);
     return (NULL);
   }
  
  lex_table_restore();
  rax_free (rawdata);

  alignmentData->siteWeights  = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i) 
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}

