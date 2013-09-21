#include <stdio.h>
#include <stdlib.h>
#include "../../pll.h"

extern void stack_dump(pllStack ** stack);

int main (int argc, char * argv[])
{
  pllNewickTree * tree;
  //pllInstance * t;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s FILENAME\n", argv[0]);
     return (EXIT_FAILURE);
   }


  if ((tree = pllNewickParseFile (argv[1])))
   {
     printf ("Parsing successful...\n\n");

     //if (pllValidateNewick (stack, nodes, leaves))
     if (pllValidateNewick (tree))
      {
        printf ("Valid phylogenetic tree\n");
      }
     else
       printf ("Not a valid phylogenetic tree\n");

     stack_dump(&(tree->tree));

     //t = pllTreeCreateNewick (stack, nodes, leaves);
     pllNewickParseDestroy (&tree);
   //pllTreeDestroy (t);
   }
  else
    printf ("Error while parsing newick tree...\n");


  return (EXIT_SUCCESS);
}
