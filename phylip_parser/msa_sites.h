#ifndef MSA_SITES
#define MSA_SITES

#define SITES_CREATE  0x00000000
#define SITES_COMPACT 0x00000001
#define SITES_SORTED  0x00000002

struct msa_sites
 {
   int          taxa;
   int          seqlen;
   char      ** label;
   char      ** unique_site;
   char      ** site;
 };

void dump_sites (struct msa_sites * ms);

struct msa_sites * alloc_sites_struct (int taxa, int seqlen);
void free_sites_struct (struct msa_sites * ms);
struct msa_sites * construct_msa_sites (struct phylip_data * pd, int flags);

#endif
