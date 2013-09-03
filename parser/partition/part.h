#ifndef __pll_PART__
#define __pll_PART__
#include "../../queue.h"

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
  int protModels;
  int protFreqs;
  int dataType;
  int optimizeBaseFrequencies;
  struct pllQueue * regionList;
};
#endif
