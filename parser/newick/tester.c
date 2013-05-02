#include <stdio.h>
#include <stdlib.h>
#include "newick.h"

int main (int argc, char * argv[])
{
  int nodes, leaves;
  struct pllStack * stack = NULL;
  //pllInstance * t;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s FILENAME\n", argv[0]);
     return (EXIT_FAILURE);
   }


  if (pllNewickParseFile (argv[1], &stack, &nodes, &leaves))
   {
     printf ("Parsing successful...\n\n");

     //if (pllValidateNewick (stack, nodes, leaves))
     if (pllValidateNewick (stack, nodes, leaves))
      {
        printf ("Valid phylogenetic tree\n");
      }
     else
       printf ("Not a valid phylogenetic tree\n");

     stack_dump(&stack);

     //t = pllTreeCreateNewick (stack, nodes, leaves);
     pllNewickParseDestroy (&stack);
     //pllTreeDestroy (t);
   }
  else
    printf ("Error while parsing newick tree...\n");


  return (EXIT_SUCCESS);
}
