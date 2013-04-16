#ifndef __pll_PART__
#define __pll_PART__
#include "../../lexer.h"
#include "../../axml.h"
#include "../../queue.h"
#include "../../mem_alloc.h"
#include "../../hash.h"

struct pllPartitionRegion
{
  int start;
  int end;
  int stride;
};

struct pllPartitionInfo
{
  char * partitionName;
  char * partitionModel;
  struct pllQueue * regionList;
};

void  pllQueuePartitionsDestroy (struct pllQueue ** partitions);
struct pllQueue * pllPartitionParse (const char * filename);
void pllPartitionDump (struct pllQueue * partitions);
#endif
