#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"
#include "xalloc.h"
#include "msa_sites.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

#define CONSUME(x)         while (token.class & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

struct rawdata
 {
   unsigned char ** oa;         /* original alignment */
   unsigned char ** pats;       /* unique site patterns */
 };


static char * 
read_file (const char * filename, int * n)
{
  FILE       * fp;
  char       * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *n = ftell (fp);
  if (*n == -1) return (NULL);
  rewind (fp);

  rawdata = (char *) xmalloc ((*n) * sizeof (char));
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static int
print_tokens (int input)
{
  struct lex_token token;

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

struct phylip_data *
alloc_phylip_struct (int taxa, int seqlen)
 {
   int i;
   struct phylip_data * pd;

   pd = (struct phylip_data *) malloc (sizeof (struct phylip_data));

   pd->taxa   = taxa;
   pd->seqlen = seqlen;
   pd->label  = (char **) calloc (taxa, sizeof (char *));
   pd->seq    = (char **) calloc (taxa, sizeof (char *));

   for (i = 0; i < taxa; ++ i)
    {
      pd->seq[i] = (char *) calloc ((seqlen + 1), sizeof (char));
    }

   return (pd);
 }

void
free_phylip_struct (struct phylip_data * pd)
{
  int i;

  for (i = 0; i < pd->taxa; ++ i)
   {
     free (pd->label[i]);
     free (pd->seq[i]);
   }
  free (pd->label);
  free (pd->seq);
  free (pd);
}


void 
dump_struct (struct phylip_data * pd)
 {
   int i;

   printf ("=> Dumping phylip_data\n");
   printf ("%d %d\n", pd->taxa, pd->seqlen);
   for (i = 0; i < pd->taxa; ++ i)
    {
      printf ("|%s| |%s|\n", pd->label[i], pd->seq[i]);
    }
 }

static struct phylip_data *
read_phylip_header (char * rawdata, int * inp)
{
  struct lex_token token;
  struct phylip_data * pd;
  int taxa, seqlen, input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

  if (token.class != LEX_NUMBER) return (NULL);

  taxa = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  if (token.class != LEX_NUMBER) return (NULL);

  seqlen = atoi (token.lexeme);
//  NEXT_TOKEN
//  CONSUME(LEX_WHITESPACE)
//  CONSUME(LEX_NEWLINE)

  /* allocate memory for the alignment */

  if (!taxa || !seqlen)
   {
     return (NULL);
   }

  pd   = alloc_phylip_struct (taxa, seqlen);
  *inp = input;

  return (pd);
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
static int
parse_phylip_sequential (char * rawdata, struct phylip_data * pd, int input)
{
  int i;
  struct lex_token token;
  int * seq_size;

  seq_size = (int *) calloc (pd->taxa, sizeof (int));

  NEXT_TOKEN

  /* first parse the TAXA-LABEL SEQUENCE format */
  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN && i < pd->taxa)
   {
     CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

     /* do something with the label */
     if (token.class != LEX_STRING && token.class != LEX_NUMBER) 
      {
        FREE (2, rawdata, seq_size);
        return (0);
      }
     pd->label[i] = strndup(token.lexeme, token.len);

     NEXT_TOKEN
     do 
      {
        CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

        /* do something with the sequence(s) */
        if (token.class != LEX_STRING)
         {
           FREE (2, rawdata, seq_size);
           return (0);
         }
        seq_size[i] += token.len;
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }

        NEXT_TOKEN
      } while (seq_size[i] < pd->seqlen);

     if (seq_size[i] > pd->seqlen) break;
     ++i;
   }
 if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
  {
    FREE (2, rawdata, seq_size);
    return (0);
  }
 
 FREE (2, seq_size, rawdata);

  return (1);
}

static int
parse_phylip_interleaved (char * rawdata, struct phylip_data * pd, int input)
{
  int i;
  struct lex_token token;
  int * seq_size;

  seq_size = (int *) calloc (pd->taxa, sizeof (int));

  NEXT_TOKEN

  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN && i < pd->taxa)
   {
     CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
     
     if (token.class != LEX_STRING && token.class != LEX_NUMBER)
      {
        FREE (2, rawdata, seq_size);
        return (0);
      }
     pd->label[i] = strndup(token.lexeme, token.len);

     NEXT_TOKEN

     while (token.class == LEX_WHITESPACE)
      {
        CONSUME(LEX_WHITESPACE)
        if (token.class != LEX_STRING)
         {
           FREE (2, rawdata, seq_size);
           return (0);
         }
        seq_size[i] += token.len;
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }
        NEXT_TOKEN
      }
     if (seq_size[i] > pd->seqlen) break;
     ++ i;
   }
  if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
   {
     FREE (2, rawdata, seq_size);
     return (0);
   }

  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN)
   {
     if (i % pd->taxa == 0) i = 0;

     CONSUME(LEX_WHITESPACE | LEX_NEWLINE);

     do
      {
        CONSUME(LEX_WHITESPACE)
        if (token.class != LEX_STRING)
         {
           FREE (2, rawdata, seq_size);
           return (0);
         }
        seq_size[i] += token.len;
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }
        NEXT_TOKEN
      } while (token.class == LEX_WHITESPACE);
     if (seq_size[i] > pd->seqlen) break;
     ++i;
     CONSUME(LEX_NEWLINE);
   }
  if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
   {
     FREE (2, rawdata, seq_size);
     return (0);
   }

  for (i = 0; i < pd->taxa; ++ i)
   {
     if (seq_size[i] != pd->seqlen) break;
   }

  FREE (2, seq_size, rawdata);

  if (i < pd->taxa) return (0);

  return (1);

}

struct phylip_data *
parse_phylip (const char * phyfile, int type)
{
  int n, input;
  char * rawdata;
  struct phylip_data * pd;

  rawdata = read_file (phyfile, &n);
  if (!rawdata)
   {
     return (0);
   }

//  printf ("=> Raw data\n");
//  printf ("%s\n", rawdata);

  init_lexan (rawdata, n);
  input = get_next_symbol();

//    return (print_tokens(input));
  
  pd = read_phylip_header (rawdata, &input);
  if (!pd)
   {
     free (rawdata);
     return (0);
   }
  
  switch (type)
   {
     case PHYLIP_SEQUENTIAL:
       if (!parse_phylip_sequential(rawdata, pd, input))
        {
          free_phylip_struct (pd);
          return (NULL);
        }
       break;
     case PHYLIP_INTERLEAVED:
       if (!parse_phylip_interleaved(rawdata, pd, input))
        {
          free_phylip_struct (pd);
          return (NULL);
        }
       break;
     default:
       break;
   }

  return (pd);
}

void 
usage (const char * cmd_name)
{
  fprintf (stderr, "Usage: %s [PHYLIP-FILE] [PHYLIP_INTERLEAVED | PHYLIP_SEQUENTIAL]\n", cmd_name);
}

int 
main (int argc, char * argv[])
{
  struct phylip_data * pd;
  struct msa_sites * ms;
  int format;

  if (argc != 3)
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }
  
  if (strcmp (argv[2], "PHYLIP_INTERLEAVED") && strcmp (argv[2], "PHYLIP_SEQUENTIAL"))
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }

  format = strcmp (argv[2], "PHYLIP_INTERLEAVED") ? PHYLIP_SEQUENTIAL : PHYLIP_INTERLEAVED;

  pd = parse_phylip (argv[1], format);
  if (!pd) 
   {
     return (EXIT_FAILURE);
   }
  
//  dump_struct (pd);

//  printf ("=> Sorting\n");
  ms = construct_msa_sites (pd, SITES_CREATE | SITES_COMPACT | SITES_SORTED);
  dump_sites (ms);

  free_phylip_struct (pd);
  free_sites_struct (ms);



  return (EXIT_SUCCESS);
}
