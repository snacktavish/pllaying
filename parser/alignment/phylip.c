#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylip.h"

static int
read_phylip_header (char * rawdata, int * inp, int * sequenceCount, int * sequenceLength)
{
  struct ltoken_t token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

  if (token.class != LEX_NUMBER) return (0);

  *sequenceCount = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  if (token.class != LEX_NUMBER) return (0);

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
parse_phylip (char * rawdata, pllAlignmentData * alignmentData, int input)
{
  int i,j;
  struct ltoken_t token;
  int * sequenceLength;
  int rc;

  sequenceLength = (int *) rax_calloc (alignmentData->sequenceCount + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % alignmentData->sequenceCount;
    if (i < alignmentData->sequenceCount) 
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
          rax_free (sequenceLength);
          return (0);
        }

       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)


       if (token.class != LEX_STRING && token.class != LEX_NUMBER && token.class != LEX_FLOAT)
        {
          rax_free (sequenceLength);
          return (0);
        }
       alignmentData->sequenceLabels[i + 1] = strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
     }
    
    while (1)
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
         rax_free (sequenceLength);
         return (0);
        }
       
       if (token.class == LEX_NEWLINE) break;

       if (token.class != LEX_STRING)
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
       CONSUME (LEX_WHITESPACE)
     }
    CONSUME(LEX_WHITESPACE | LEX_NEWLINE);
  }
}

/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL |
 *             LEX_WHITESPACE LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL
 *     ENDL -> LEX_WHITESPACE LEX_NEWLINE | LEX_NEWLINE
 *     DATA -> LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA | 
 *             LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF
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

  rawdata = __pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (NULL);
   }
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (rawdata, &input, &sequenceCount, &sequenceLength))
   {
     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     return (NULL);
   }

  lex_table_amend_phylip();

  /* allocate alignment structure */
  alignmentData = pllInitAlignmentData (sequenceCount, sequenceLength);

  if (! parse_phylip (rawdata, alignmentData, input))
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

