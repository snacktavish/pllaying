#ifndef __pll_FASTA__
#define __pll_FASTA__

#include "alignment.h"
#include "../../lexer.h"

pllAlignmentData * pllParseFASTA (const char *);
char * pllReadFile (const char *, int *);

#endif
