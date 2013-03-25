#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lexer.h"
#include <math.h>
#include "../../axml.h"
#include "../../queue.h"
#include "../../mem_alloc.h"
#include "part.h"

#define CONSUME(x)         while (token.class & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);


void
pllPartitionsDestroy (struct pllQueue ** partitions)
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
     free (pi->partitionName);
     free (pi->partitionModel);
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
  int input;
  struct ltoken_t token;
  int lines = 0;
  struct pllQueue * partitions;
  struct pllPartitionInfo * pi;
  struct pllPartitionRegion * region;

  input  = *inp;

  NEXT_TOKEN

  pllQueueInit (&partitions);
  while (token.class != LEX_EOF)
  {
    ++ lines;
    pi = (struct pllPartitionInfo *) rax_malloc (sizeof (struct pllPartitionInfo));
    pllQueueInit (&(pi->regionList));
    CONSUME (LEX_WHITESPACE | LEX_NEWLINE)


    /* read partition type */
    if (token.class != LEX_STRING) return (0);
    pi->partitionModel = strndup (token.lexeme, token.len);
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)


    if (token.class != LEX_COMMA) return (0);
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read partition name */
    if (token.class != LEX_STRING) return (0);
    pi->partitionName = strndup (token.lexeme, token.len);
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read equal sign */
    if (token.class != LEX_EQUAL)
     {
       return (0);
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE)

    /* read rhs */
    while (1)
    {
      region = (struct pllPartitionRegion *) rax_malloc (sizeof (struct pllPartitionRegion));
      if (token.class != LEX_NUMBER) return (0);
      region->start  = region->end = atoi (token.lexeme);  
      region->stride = 1;
      NEXT_TOKEN
      CONSUME(LEX_WHITESPACE)
      
      if  (token.class == LEX_DASH)
       {
         NEXT_TOKEN
         CONSUME(LEX_WHITESPACE)
         if (token.class != LEX_NUMBER) return (0);
         region->end = atoi (token.lexeme);
         NEXT_TOKEN
         CONSUME(LEX_WHITESPACE)
         if (token.class == LEX_SLASH)
          {
            NEXT_TOKEN
            CONSUME(LEX_WHITESPACE)
            if (token.class != LEX_NUMBER) return (0);
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
   pllQueueAppend (partitions, (void *)pi);
   CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  }
 
 return (partitions);
} 

void pllPartitionDump (struct pllQueue * partitions)
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

int
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

  partitions = parse_partition (rawdata, &input);
  if (partitions)
   {
     printf ("Parsed successfully...\n");
     pllPartitionDump (partitions);
     pllPartitionsDestroy (&partitions);
   }
  else
   {
     printf ("Error while parsing...\n");
   }
  
  rax_free (rawdata);
  return (1);
}


int main (int argc, char * argv[])
{
  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s FILENAME\n", argv[0]);
     return (EXIT_FAILURE);
   }

  pllPartitionParse (argv[1]);

  return (EXIT_SUCCESS);
}
