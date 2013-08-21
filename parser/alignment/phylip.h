#ifndef __pll_PHYLIP__
#define __pll_PHYLIP__

#include "alignment.h"
#include "../../lexer.h"

pllAlignmentData * pllParsePHYLIP (const char *);
char * pllReadFile (const char *, int *);

#endif
