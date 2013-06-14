#include <stdio.h>
#include <stdlib.h>
#include "part.h"

int main (int argc, char * argv[])
{
  struct pllQueue * partitions;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s FILENAME\n", argv[0]);
     return (EXIT_FAILURE);
   }

  partitions = pllPartitionParse (argv[1]);
  if (partitions)
   {
     printf ("Parsed successfully...\n");
     pllPartitionDump (partitions);
//     pllPartitionsDestroy (&partitions);
   }
  else
   {
     printf ("Error while parsing...\n");
   }

  return (EXIT_SUCCESS);
}
