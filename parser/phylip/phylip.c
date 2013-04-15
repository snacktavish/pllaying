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

void pllPhylipEF (struct pllPhylip * phylip, double ** ef)
{
  int i, j;
  int a, c, g, t, und;
  int pa, pc, pg, pt, pund;
  double total;
  a = c = g = t = und = 0;

  total = 0;



  for (i = 0; i < phylip->seqLen; ++ i)
  {
    for (j = 1; j <= phylip->nTaxa; ++ j)
    {
      //assert (phylip->seq[i][j] == 'A' || phylip->seq[i][j] == 'C' || phylip->seq[i][j] == 'G' || phylip->seq[i][j] == 'T' || phylip->seq[i][j] == '-');
      pa = pc = pg = pt = pund = 0;

      switch (phylip->seq[j][i])
       {
         case 'a':
         case 'A':
           ++ pa;
           break;

         case 'c':
         case 'C':
           ++ pc;
           break;

         case 'g':
         case 'G':
           ++ pg;
           break;

         case 't':
         case 'T':
           ++ pt;
           break;
         
         case '-':
         case '?':
         case 'n':
         case 'N':
           ++ pund;
           break;
         
         default:
           fprintf (stderr, "Error, unidentified base %c\n", phylip->seq[j][i]);
           exit(1);
       }
      a += pa * phylip->weights[i];
      c += pc * phylip->weights[i];
      g += pg * phylip->weights[i];
      t += pt * phylip->weights[i];
      und += pund * phylip->weights[i];
    }
  }
  total = a + c + g + t;
  ef[0][0] = a / total; ef[0][1] = c / total; ef[0][2] = g / total; ef[0][3] = t / total;
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
pllPhylipSubst (struct pllPhylip * phylip, int type)
{
  unsigned char meaningDNA[256];
  unsigned char  meaningAA[256];
  unsigned char * d;
  int i, j;

  for (i = 0; i < 256; ++ i)
   {
     meaningDNA[i] = -1;
     meaningAA[i]  = -1;
   }

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['Y'] = 10;
  meaningDNA['a'] =  1;
  meaningDNA['b'] = 14;
  meaningDNA['c'] =  2;
  meaningDNA['d'] = 13;
  meaningDNA['g'] =  4;
  meaningDNA['h'] = 11;
  meaningDNA['k'] = 12;
  meaningDNA['m'] =  3;
  meaningDNA['r'] =  5;
  meaningDNA['s'] =  6;
  meaningDNA['t'] =  8;
  meaningDNA['u'] =  8;
  meaningDNA['v'] =  7;
  meaningDNA['w'] =  9;
  meaningDNA['y'] = 10;

  meaningDNA['N'] =
  meaningDNA['n'] =
  meaningDNA['O'] =
  meaningDNA['o'] =
  meaningDNA['X'] =
  meaningDNA['x'] =
  meaningDNA['-'] =
  meaningDNA['?'] = 15;
 
  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['a'] =  0;  /* alanine */
  meaningAA['r'] =  1;  /* arginine */
  meaningAA['n'] =  2;  /*  asparagine*/
  meaningAA['d'] =  3;  /* aspartic */
  meaningAA['c'] =  4;  /* cysteine */
  meaningAA['q'] =  5;  /* glutamine */
  meaningAA['e'] =  6;  /* glutamic */
  meaningAA['g'] =  7;  /* glycine */
  meaningAA['h'] =  8;  /* histidine */
  meaningAA['i'] =  9;  /* isoleucine */
  meaningAA['l'] =  10; /* leucine */
  meaningAA['k'] =  11; /* lysine */
  meaningAA['m'] =  12; /* methionine */
  meaningAA['f'] =  13; /* phenylalanine */
  meaningAA['p'] =  14; /* proline */
  meaningAA['s'] =  15; /* serine */
  meaningAA['t'] =  16; /* threonine */
  meaningAA['w'] =  17; /* tryptophan */
  meaningAA['y'] =  18; /* tyrosine */
  meaningAA['v'] =  19; /* valine */
  meaningAA['b'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
  meaningAA['x'] = 
  meaningAA['?'] = 
  meaningAA['*'] = 
  meaningAA['-'] = 22;

  d = (type == PHYLIP_DNA_DATA) ? meaningDNA : meaningAA; 

  for (i = 1; i <= phylip->nTaxa; ++ i)
   {
     for (j = 0; j < phylip->seqLen; ++ j)
      {
        phylip->seq[i][j] = d[phylip->seq[i][j]];
      }
   }
}

void 
usage (const char * cmd_name)
{
  fprintf (stderr, "Usage: %s [PHYLIP-FILE]\n", cmd_name);
}

