#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "part.h"

#define GLOBAL_VARIABLES_DEFINITION
#include "../../globalVariables.h"

static struct pllHashTable * hashTable;

static void destroy_model_names(void)
{
  pllHashDestroy (&hashTable, PLL_TRUE);
}

static void init_model_names (void)
{
  int i;
  int * item;

  hashTable = pllHashInit (NUM_PROT_MODELS);

  for (i = 0; i < NUM_PROT_MODELS; ++ i)
   {
     item  = (int *) rax_malloc (sizeof (int));
     *item = i;
     pllHashAdd (hashTable, protModels[i], (void *) item);
   }
}

void
pllQueuePartitionsDestroy (struct pllQueue ** partitions)
{
  struct pllPartitionInfo * pi;
  struct pllPartitionRegion * region;

  while (pllQueueRemove (*partitions, (void **)&pi))
   {
     while (pllQueueRemove (pi->regionList, (void **) &region))
      {
        rax_free (region);
      }
     rax_free (pi->regionList);
     rax_free (pi->partitionName);
     rax_free (pi->partitionModel);
     rax_free (pi);
   }
  rax_free (*partitions);
}

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

  rawdata = (char *) rax_malloc (((*n)  + 1)* sizeof (char));
  rawdata[*n] = 0;
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static struct pllQueue *
parse_partition (char * rawdata, int * inp)
{
  int input, i;
  struct ltoken_t token;
  int lines = 0;
  struct pllQueue * partitions;
  struct pllPartitionInfo * pi;
  struct pllPartitionRegion * region;
  int * item;

  input  = *inp;

  NEXT_TOKEN

  pllQueueInit (&partitions);
  while (token.class != LEX_EOF)
  {
    ++ lines;
    pi = (struct pllPartitionInfo *) rax_calloc (1, sizeof (struct pllPartitionInfo));
    pllQueueInit (&(pi->regionList));
    pllQueueAppend (partitions, (void *)pi);
    CONSUME (LEX_WHITESPACE | LEX_NEWLINE)


    /* read partition type */
    if (token.class != LEX_STRING) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    //pi->partitionModel = strndup (token.lexeme, token.len);
    pi->partitionModel = (char *) rax_malloc ((token.len + 1) * sizeof (char));
    strncpy (pi->partitionModel, token.lexeme, token.len);
    pi->partitionModel[token.len] = 0;
    for (i = 0; i < token.len; ++i) pi->partitionModel[i] = toupper(pi->partitionModel[i]);
    // check partition model
    pi->protModels = 0;
    if (strcmp (pi->partitionModel, "DNA"))
     {
       if (! pllHashSearch (hashTable, pi->partitionModel, (void **)&item))
        {
          pllQueuePartitionsDestroy (&partitions);
          return (0);
        }
       pi->protModels = *item;
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    if (token.class != LEX_COMMA) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read partition name */
    if (token.class != LEX_STRING) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    //pi->partitionName = strndup (token.lexeme, token.len);
    pi->partitionName = (char *) rax_malloc ((token.len + 1) * sizeof (char));
    strncpy (pi->partitionName, token.lexeme, token.len);
    pi->partitionName[token.len] = 0;

    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read equal sign */
    if (token.class != LEX_EQUAL)
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read rhs */
    while (1)
    {
      region = (struct pllPartitionRegion *) rax_malloc (sizeof (struct pllPartitionRegion));
      if (token.class != LEX_NUMBER) 
       {
         pllQueuePartitionsDestroy (&partitions);
         return (0);
       }
      region->start  = region->end = atoi (token.lexeme);  
      region->stride = 1;
      NEXT_TOKEN
      CONSUME(LEX_WHITESPACE)
      
      if  (token.class == LEX_DASH)
       {
         NEXT_TOKEN
         CONSUME(LEX_WHITESPACE)
         if (token.class != LEX_NUMBER) 
          {
            pllQueuePartitionsDestroy (&partitions);
            return (0);
          }
         region->end = atoi (token.lexeme);
         if (region->end < region->start)
          {
            pllQueuePartitionsDestroy (&partitions);
            return (0);
          }
         NEXT_TOKEN
         CONSUME(LEX_WHITESPACE)
         if (token.class == LEX_SLASH)
          {
            NEXT_TOKEN
            CONSUME(LEX_WHITESPACE)
            if (token.class != LEX_NUMBER) 
             {
               pllQueuePartitionsDestroy (&partitions);
               return (0);
             }
            region->stride = atoi (token.lexeme);
            NEXT_TOKEN
          }
         CONSUME(LEX_WHITESPACE)
       }
       pllQueueAppend (pi->regionList, (void *)region);
      
      if (token.class != LEX_COMMA) break;
      NEXT_TOKEN
      CONSUME(LEX_WHITESPACE)
    }
   CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  }
 
 return (partitions);
} 

void 
pllPartitionDump (struct pllQueue * partitions)
{
   struct pllQueueItem * elm;
   struct pllQueueItem * regionList;
   struct pllPartitionInfo * pi;
   struct pllPartitionRegion * region;

   elm = partitions->head;

   while (elm)
    {
      pi  = (struct pllPartitionInfo *) elm->item;
      printf ("%s, %s = ", pi->partitionModel, pi->partitionName);
      regionList = pi->regionList->head;
      while (regionList)
       {
         region = (struct pllPartitionRegion *) regionList->item;
         printf ("%d", region->start);
         if (region->start != region->end)
          {
            printf ("-%d", region->end);
            if (region->stride != 1) printf ("/%d", region->stride);
          }
         regionList = regionList->next;
         if (regionList) printf (", ");
       }
      printf ("\n");

      elm = elm->next;
    }
}

struct pllQueue *
pllPartitionParse (const char * filename)
{
  int n;
  char * rawdata;
  int input;
  struct pllQueue * partitions;

  rawdata = readFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  printf ("%s\n\n", rawdata);

  n = strlen (rawdata);

  init_lexan (rawdata, n);
  input = get_next_symbol();

  init_model_names();
  partitions = parse_partition (rawdata, &input);
  destroy_model_names();
  
  rax_free (rawdata);
  return (partitions);
}
