#ifndef PHYLIP_H
#define PHYLIP_H

#define PHYLIP_SEQUENTIAL       0x00000001
#define PHYLIP_INTERLEAVED      0x00000002

struct phylip_data
 {
   int          taxa;
   int          seqlen;
   char      ** label;
   char      ** seq;
   int        * weight;
 };

struct phylip_data * alloc_phylip_struct (int taxa, int seqlen);

#endif
