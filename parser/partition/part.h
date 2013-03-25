#ifndef __pll_PART__
#define __pll_PART__

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

#endif
