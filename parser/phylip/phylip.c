#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "phylip.h"
#include "ssort.h"

//struct rawdata
// {
//   unsigned char ** oa;         /* original alignment */
//   unsigned char ** pats;       /* unique site patterns */
// };

static char * 
readFile (const char * filename, int * n)
{
  FILE * fp;
  char * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *n = ftell (fp);
  if (*n == -1) return (NULL);
  rewind (fp);

  rawdata = (char *) rax_malloc ((*n) * sizeof (char));
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

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

static struct pllPhylip *
alloc_phylip_struct (int nTaxa, int seqLen)
 {
   int i;
   struct pllPhylip * phylip;
   void * mem;
   
   /** TODO */
   phylip = (struct pllPhylip *) rax_malloc (sizeof (struct pllPhylip));
   phylip->seq = (unsigned char **) rax_malloc ((nTaxa + 1) * sizeof (unsigned char *));
   mem = rax_malloc (sizeof (unsigned char) * (seqLen + 1) * nTaxa);
   for (i = 1; i <= nTaxa; ++i)
    {
      phylip->seq[i] = (unsigned char *) (mem + (i - 1) * (seqLen + 1) * sizeof (unsigned char));
      phylip->seq[i][seqLen] = 0;
    }
   phylip->seq[0] = NULL;
    
   phylip->label = (char **) rax_calloc ((nTaxa + 1), sizeof (char *));

   phylip->nTaxa   = nTaxa;
   phylip->seqLen  = seqLen;
   phylip->weights = NULL;

   return (phylip);
 }

void
pllPhylipDestroy (struct pllPhylip * phylip)
{
  int i;

  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     rax_free (phylip->label[i]);
   }
  rax_free (phylip->label);
  rax_free (phylip->seq[1]);
  rax_free (phylip->seq);
  rax_free (phylip->weights);
  rax_free (phylip);
}


static void 
dump_struct (struct pllPhylip * pd)
 {
   int i;

   printf ("=> Dumping pllPhylip\n");
   printf ("%d %d\n", pd->nTaxa, pd->seqLen);
   for (i = 0; i < pd->nTaxa; ++ i)
    {
      printf ("|%s| |%s|\n", pd->label[i], pd->seq[i]);
    }
 }

static int
read_phylip_header (char * rawdata, int * inp, int * nTaxa, int * seqLen)
{
  struct ltoken_t token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

  if (token.class != LEX_NUMBER) return (0);

  *nTaxa = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  if (token.class != LEX_NUMBER) return (0);

  *seqLen = atoi (token.lexeme);

  *inp = input;

  return (*nTaxa && *seqLen);
}

static inline int
parsedOk (int * actLen, int nTaxa, int seqLen  )
{
  int i;

  for (i = 1; i <= nTaxa; ++ i)
   {
     if (actLen[i] != seqLen) return (0);
   }
  
  return (1);
}


static int
parse_phylip (char * rawdata, struct pllPhylip * phylip, int input)
{
  int i,j;
  struct ltoken_t token;
  int * seqLen;
  int rc;

  seqLen = (int *) rax_calloc (phylip->nTaxa + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % phylip->nTaxa;
    if (i < phylip->nTaxa) 
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (seqLen, phylip->nTaxa, phylip->seqLen);
          rax_free (seqLen);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
          rax_free (seqLen);
          return (0);
        }

       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)


       if (token.class != LEX_STRING && token.class != LEX_NUMBER && token.class != LEX_FLOAT)
        {
          rax_free (seqLen);
          return (0);
        }
       phylip->label[i + 1] = strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
     }
    
    while (1)
     {
       if (token.class == LEX_EOF)
        {
          rc = parsedOk (seqLen, phylip->nTaxa, phylip->seqLen);
          rax_free (seqLen);
          return (rc);
        }

       if (token.class == LEX_UNKNOWN)
        {
         rax_free (seqLen);
         return (0);
        }
       
       if (token.class == LEX_NEWLINE) break;

       if (token.class != LEX_STRING)
        {
          rax_free (seqLen);
          return (0);
        }

       if (seqLen[j + 1] + token.len > phylip->seqLen) 
        {
          fprintf (stderr, "Sequence %d is larger than specified\n", j + 1);
          rax_free (seqLen);
          return (0);
        }
       memmove (phylip->seq[j + 1] + seqLen[j + 1], token.lexeme, token.len);
       seqLen[j + 1] += token.len;

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
struct pllPhylip *
pllPhylipParse (const char * filename)
{
  int i, n, input, nTaxa, seqLen;
  char * rawdata;
  struct pllPhylip * phylip;

  rawdata = readFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }
  
  init_lexan (rawdata, n);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (rawdata, &input, &nTaxa, &seqLen))
   {
     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     return (0);
   }

  lex_table_amend_phylip();
  /* allocate the phylip structure */
  phylip = alloc_phylip_struct (nTaxa, seqLen);

  if (! parse_phylip (rawdata, phylip, input))
   {
     printf ("Finished with error in parsing ...\n");
     pllPhylipDestroy (phylip);
     lex_table_restore();
     rax_free (rawdata);
     return (0);
   }
  
  lex_table_restore();
  rax_free (rawdata);

  phylip->weights  = (int *) rax_malloc (phylip->seqLen * sizeof (int));
  for (i = 0; i < phylip->seqLen; ++ i) phylip->weights[i] = 1;
  return (phylip);
}

void
pllPhylipDump (struct pllPhylip * phylip)
{
  FILE * fp;
  int i;

  fp = fopen ("output.seq.phy","w");
  fprintf (fp, "%d %d\n", phylip->nTaxa, phylip->seqLen);
  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     fprintf (fp, "%s %s\n", phylip->label[i], phylip->seq[i]);
   }

  fclose (fp);
}
void 
usage (const char * cmd_name)
{
  fprintf (stderr, "Usage: %s [PHYLIP-FILE]\n", cmd_name);
}

