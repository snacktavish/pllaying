#include <stdio.h>
#include <string.h>
#include "hash.h"
#include "mem_alloc.h"

/** @brief Generate the hash value for a string 

    Generates the hash value of a string \a s.

    @param s
      The string for which the hash is computed

    @param size
      Size of the hash table

  @return
      Hash of string \a s, i.e. index at hash table
*/
unsigned int 
pllHashString (char * s, unsigned int size)
{
  unsigned int hash = 0;

  for (; *s; ++s) hash = (hash << 5) - hash + (unsigned int )*s;

  return (hash % size);
}

/** @brief Add a string and its data to a hashtable
    
    Add an \a item associated with string \a s to hashtable \a hTable.

    @param hTable
      Hashtable

    @param s
      String

    @param item
      Data associated with \a s

    @return
      Returns \b 1 if added with success, otherwise \b 0
*/
int
pllHashAdd  (struct pllHashTable * hTable, char * s, void * item)
{
  unsigned int pos;
  struct pllHashItem * hItem;

  pos   = pllHashString (s, hTable->size);
  hItem = hTable->Items[pos];

  for (; hItem; hItem = hItem->next)
   {
     if (!strcmp (s, hItem->str)) return (0);
   }

  hItem = (struct pllHashItem *) rax_malloc (sizeof (struct pllHashItem));

  hItem->str = (char *) rax_malloc ((strlen(s) + 1) * sizeof (char));
  strcpy (hItem->str, s);
  hItem->data = item;

  hItem->next = hTable->Items[pos];
  hTable->Items[pos] = hItem;

  return (1);
}

       
struct pllHashTable *
pllHashInit (int n)
{ 
  struct pllHashTable * hTable;
  unsigned int i;
  unsigned int primeTableLength;
       
  static const unsigned int initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
                                             196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
                                             50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
       
  hTable = (struct pllHashTable *) rax_malloc (sizeof (struct pllHashTable));
  
  primeTableLength = sizeof (initTable) / sizeof(initTable[0]);

  i = 0;
 
  while (initTable[i] < n && i < primeTableLength) ++ i;
 
  n = initTable[i];  
 
  hTable->Items = (struct pllHashItem **) rax_calloc (n, sizeof (struct pllHashItem *));
  hTable->size  = n;
 
  return (hTable);
}

int
pllHashSearch (struct pllHashTable * hTable, char * s, void ** item)
{
  unsigned int pos;
  struct pllHashItem * hItem;

  pos   = pllHashString (s, hTable->size);
  hItem = hTable->Items[pos];

  for (; hItem; hItem = hItem->next)
   {
     if (!strcmp (s, hItem->str))
      {
        *item = hItem->data;
        return (1);
      }
   }

  return (0);
}

void 
pllHashDestroy (struct pllHashTable ** hTable)
{
  int i;
  struct pllHashItem * hItem;
  struct pllHashItem * tmp;

  for (i = 0; i < (*hTable)->size; ++ i)
  {
    hItem = (*hTable)->Items[i];
    while (hItem)
     {
       tmp   = hItem;
       hItem = hItem->next;
       rax_free (tmp->str);
       rax_free (tmp);
     }
  }
  rax_free ((*hTable)->Items);
  rax_free (*hTable);
  *hTable = NULL;
}
