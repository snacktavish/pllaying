#include <stdio.h>
#include "../../axml.h"
#include "../../utils.h"

int main (int argc, char * argv [])
{
  pllInstance * tr;
  
  tr = pllCreateInstance (GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

  printf ("Address of tr in main() %p\n", tr);
  printf ("Address of tr->nodep in main() %p\n", tr->nodep);
  printf ("Address of tr->start in main() %p\n\n", tr->start);

  pllDummy (tr);

  printf ("\nAddress of tr in main() %p\n", tr);
  printf ("Address of tr->nodep in main() %p\n", tr->nodep);
  printf ("Address of tr->start in main() %p\n", tr->start);
}
