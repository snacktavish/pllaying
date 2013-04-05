#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "newick.h"

struct item_t
{
  int depth;
  char * name;
  char * branch;
  int leaf;
  int rank;
};

static char * 
readFile (const char * filename, int * n)
{
  FILE * fp;
  char * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *n = ftell (fp);
  if (*n == -1) return (NULL);
  rewind (fp);

  rawdata = (char *) rax_malloc (((*n)  + 1)* sizeof (char));
  rawdata[*n] = 0;
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static int
parse_newick (char * rawdata, struct pllStack ** stack, int * inp)
{
  struct item_t * item = NULL;
  int item_active = 0;
  struct ltoken_t token;
  int input;
  struct ltoken_t prev_token;
  int nop = 0;          /* number of open parentheses */
  int depth = 0;

  prev_token.class = LEX_UNKNOWN;

  input = *inp;

  NEXT_TOKEN
  
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN)
  {
    switch (token.class)
     {
       case LEX_OPAREN:
       //printf ("LEX_OPAREN\n");
        ++nop;
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));
        ++depth;
        break;

       case LEX_CPAREN:
       //printf ("LEX_CPAREN\n");
        if (prev_token.class != LEX_CPAREN  &&
            prev_token.class != LEX_UNKNOWN &&
            prev_token.class != LEX_STRING  &&
            prev_token.class != LEX_NUMBER  &&
            prev_token.class != LEX_FLOAT) return (0);

        if (!nop) return (0);
        --nop;
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));

        /* push to the stack */
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t)); // possibly not nec
        //if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("INTERNAL_NODE") + 1) * sizeof (char));
           strcpy (item->name, "INTERNAL_NODE");
         }

        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        item->depth = depth;
        pllStackPush (stack, item);
        item_active  = 1;       /* active = 1 */
        item = NULL;
        --depth;
        break;

       case LEX_STRING:
       //printf ("LEX_STRING\n");
        if (prev_token.class != LEX_OPAREN &&
            prev_token.class != LEX_CPAREN &&
            prev_token.class != LEX_UNKNOWN &&
            prev_token.class != LEX_COMMA) return (0);
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
        //item->name = strndup (token.lexeme, token.len);
        item->name = (char *) rax_malloc ((token.len + 1) * sizeof (char));
        strncpy (item->name, token.lexeme, token.len);
        item->name[token.len] = 0;

        item_active = 1;
        item->depth = depth;
        if (prev_token.class == LEX_COMMA  ||
            prev_token.class == LEX_OPAREN ||
            prev_token.class == LEX_UNKNOWN) item->leaf = 1;
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));
        break;

       case LEX_FLOAT:
       case LEX_NUMBER:
       //if (token.class == LEX_FLOAT) printf ("LEX_FLOAT\n"); else printf ("LEX_NUMBER\n");
         if  (prev_token.class != LEX_OPAREN &&
              prev_token.class != LEX_CPAREN &&
              prev_token.class != LEX_COLON  &&
              prev_token.class != LEX_UNKNOWN &&
              prev_token.class != LEX_COMMA) return (0);
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
        if (prev_token.class == LEX_COLON)
         {
           //item->branch = strndup (token.lexeme, token.len);
           item->branch = (char *) rax_malloc ((token.len + 1) * sizeof (char));
           strncpy (item->branch, token.lexeme, token.len);
           item->branch[token.len] = 0;
         }
        else
         {
           if (prev_token.class == LEX_COMMA  ||
               prev_token.class == LEX_OPAREN ||
               prev_token.class == LEX_UNKNOWN) item->leaf = 1;
           //if (prev_token.class != LEX_UNKNOWN) ++ indent;
           //item->name = strndup (token.lexeme, token.len);
           item->name = (char *) rax_malloc ((token.len + 1) * sizeof (char));
           strncpy (item->name, token.lexeme, token.len);
           item->name[token.len] = 0;
         }
        item_active = 1;
        item->depth = depth;
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));
        break;

       case LEX_COLON:
       //printf ("LEX_COLON\n");
        if (prev_token.class != LEX_CPAREN &&
            prev_token.class != LEX_STRING &&
            prev_token.class != LEX_FLOAT  &&
            prev_token.class != LEX_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));
        break;

       case LEX_COMMA:
       //printf ("LEX_COMMA\n");
        if (prev_token.class != LEX_CPAREN &&
             prev_token.class != LEX_STRING &&
             prev_token.class != LEX_FLOAT && 
             prev_token.class != LEX_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (struct ltoken_t));
        
        /* push to the stack */
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t)); // possibly not nece
        //if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("INTERNAL_NODE") + 1) * sizeof (char));
           strcpy (item->name, "INTERNAL_NODE");
         }
        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        item->depth = depth;
        pllStackPush (stack, item);
        item_active  = 0;
        item = NULL;
        break;

       case LEX_SEMICOLON:
        //printf ("LEX_SEMICOLON\n");
        /* push to the stack */
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
        //if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("ROOT_NODE") + 1) * sizeof (char));
           strcpy (item->name, "ROOT_NODE");
         }
        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        pllStackPush (stack, item);
        item_active  = 0;
        item = NULL;
        break;
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE | LEX_NEWLINE);
  }
  if (item_active)
   {
     if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
     //if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
     if (item->name == NULL) 
      {
        item->name = (char *) rax_malloc ((strlen("ROOT_NODE") + 1) * sizeof (char));
        strcpy (item->name, "ROOT_NODE");
      }
     //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
     if (item->branch == NULL) 
      {
        item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
        strcpy (item->branch, "0.000000");
      }
     pllStackPush (stack, item);
     item_active  = 0;
   }

  if (nop) return (0);
  return (1);
}

void stack_dump(struct pllStack ** stack)
{
  struct item_t * item;
  struct pllStack * head;
  int i;

  head = *stack;
  while (head)
   {
     item = (struct item_t *) head->item;

     for (i = 0; i < item->depth; ++ i) printf ("\t");

     printf ("%s:%s\n", item->name, item->branch);

     head = head->next;
   }
}

static void
assign_ranks (struct pllStack * stack, int * nodes, int * leaves)
{
  struct pllStack * head;
  struct item_t * item, * tmp;
  struct pllStack * preorder = NULL;
  int children;
  int depth;

  *nodes = *leaves = 0;


  head = stack;
  while (head)
  {
    assert (head->item);
    item = (struct item_t *) head->item;
    
    if (item->leaf)  ++ (*leaves);

    if (preorder)
     {
       tmp = (struct item_t *) preorder->item;
       children = 0;
       while (item->depth < tmp->depth)
        {
          children = 1;
          depth = tmp->depth;
          pllStackPop (&preorder);
          tmp = preorder->item;
          while (tmp->depth == depth)
           {
             ++ children;
             pllStackPop (&preorder);
             tmp = (struct item_t *)preorder->item;
           }
          tmp->rank += children;
        }
     }
    
    ++ (*nodes);
    head = head->next;

    if (item->leaf)
     {
       if (!preorder) return;

       children = 1;
       tmp = preorder->item;
       while (tmp->depth == item->depth)
        {
          ++ children;
          pllStackPop (&preorder);
          assert (preorder);
          tmp = (struct item_t *)preorder->item;
        }
       tmp->rank += children;
     }
    else
     {
       pllStackPush (&preorder, item);
     }
  }
  
  while (preorder->item != stack->item)
  {
    item = (struct item_t *)pllStackPop (&preorder);
    tmp  = (struct item_t *) preorder->item;
    children = 1;

    while (tmp->depth == item->depth)
     {
       ++ children;
       item = (struct item_t *) pllStackPop (&preorder);
       tmp  = (struct item_t *) preorder->item;
     }
    tmp->rank += children;
    children = 0;
  }
 assert (preorder->item == stack->item);
 
 pllStackClear (&preorder);
}

/** @brief Validate if a newick tree is a valid phylogenetic tree

    A valid tree is one where the root node is binary or ternary
    and all other internal nodes are binary. In case the root
    is ternary then the tree must contain at least another internal
    node and the total number of nodes must be equal to 
    \f$ 2l - 2\f$, where \f$l\f$ is the number of leaves. If the
    root is binary, then the total number of nodes must be equal
    to \f$2l - 1\f$.

    @param tree
      Stack structure where the tree will be parsed in

    @param nodes
      Number of nodes will be returned in this variable

    @param leaves
      Number of leaves will be returned in this variable

    @return
      Returns \b 1 in case of success, otherwise \b 0
*/
int
pllValidateNewick (struct pllStack * tree, int nodes, int leaves)
{
  struct pllStack * head;
  struct item_t * item;
 
 item = tree->item;
 if (item->rank != 2 && item->rank != 3) return (0);
 head = tree->next;
 while (head)
 {
   item = head->item;
   if (item->rank != 2 &&  item->rank != 0) 
    {
      return (0);
    }
   head = head->next;
 }
 
 item = tree->item;

 if (item->rank == 2) return (nodes == 2 * leaves -1);

 return ((nodes == 2 * leaves - 2) && nodes != 4);
}

/** @brief Parse a newick tree string
  
    Parse a newick string and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children and depth.

    @param newick
      String containing the newick tree

    @param tree
      Stack structure where the tree will be parsed in

    @param nodes
      Number of nodes will be returned in this variable

    @param leaves
      Number of leaves will be returned in this variable

    @return
      Returns \b 1 in case of success, otherwise \b 0
*/
int
pllNewickParseString (char * newick, struct pllStack ** tree, int * nodes, int * leaves)
{
  int n, input, rc;

  n = strlen (newick);

  init_lexan (newick, n);
  input = get_next_symbol();

  rc = parse_newick (newick, tree, &input);
  
  assign_ranks (*tree, nodes, leaves);

  return (rc);
}

/** @brief Deallocate newick parser stack structure

    Deallocates the newick parser stack structure that represents the parsed tree. It
    also frees all memory allocated by elements of the stack structure.

    @param tree
      The tree stack structure
*/
void pllNewickParseDestroy (struct pllStack ** tree)
{
  struct item_t *  item;

  while ((item = (struct item_t *)pllStackPop (tree)))
   {
     rax_free (item->name);
     rax_free (item->branch);
     rax_free (item);
   }
}

/** @brief Deallocate PLL tree structure

    Deallocates the tree structure associated with the library and all
    its elements. 

    @param tr
      The tree structure
*/
void
pllTreeDestroy (tree * tr)
{
  int i;
  for (i = 1; i <= tr->mxtips; ++ i)
    rax_free (tr->nameList[i]);
  
  pllHashDestroy (&(tr->nameHash));
  if (tr->yVector)
   {
     if (tr->yVector[0]) rax_free (tr->yVector[0]);
     rax_free (tr->yVector);
   }
  rax_free (tr->aliaswgt);
  rax_free (tr->rateCategory);
  rax_free (tr->patrat);
  rax_free (tr->patratStored);
  rax_free (tr->lhs);
  rax_free (tr->td[0].parameterValues);
  rax_free (tr->td[0].executeModel);
  rax_free (tr->td[0].ti);
  rax_free (tr->nameList);
  rax_free (tr->nodep);
  rax_free (tr->nodeBaseAddress);
  rax_free (tr->tree_string);
  rax_free (tr->tree0);
  rax_free (tr->tree1);
  rax_free (tr);
}

/** @brief Initialize PLL tree structure with default values
    
    Initialize PLL tree structure with default values and allocate 
    memory for its elements.

    @todo
      STILL NOT FINISHED

*/
void pllTreeInitDefaults (tree * t, int nodes, int tips)
{
  nodeptr p0, p, q;
  int i, j, k;
  int inner;

  

  /* TODO: make a proper static setupTree function */

  inner = tips - 1;

  t->mxtips = tips;

  t->bigCutoff = FALSE;
  t->treeStringLength = t->mxtips * (nmlngth + 128) + 256 + t->mxtips * 2;
  t->tree_string = (char *) rax_calloc ( t->treeStringLength, sizeof(char));
  t->tree0 = (char*)rax_calloc((size_t)t->treeStringLength, sizeof(char));
  t->tree1 = (char*)rax_calloc((size_t)t->treeStringLength, sizeof(char));

  
  p0 = (nodeptr) rax_malloc ((tips + 3 * inner) * sizeof (node));
  assert (p0);

  t->nodeBaseAddress  = p0;

  t->nameList         = (char **)   rax_malloc ((tips + 1) * sizeof (char *));
  t->nodep            = (nodeptr *) rax_malloc ((2 * tips) * sizeof (nodeptr));
  assert (t->nameList && t->nodep);

  t->nodep[0] = NULL;          


  /* TODO: FIX THIS! */
  t->fracchange = -1;

  for (i = 1; i <= tips; ++ i)
   {
     p = p0++;

     //p->hash      = KISS32();     
     p->x         = 0;
     p->xBips     = 0;
     p->number    = i;
     p->next      = p;
     p->back      = NULL;
     p->bInf      = NULL;
     t->nodep[i]  = p;
   }

  for (i = tips + 1; i <= tips + inner; ++i)
   {
     q = NULL;
     for (j = 1; j <= 3; ++ j)
     {
       p = p0++;
       if (j == 1)
        {
          p->xBips = 1;
          p->x = 1; //p->x     = 1;
        }
       else
        {
          p->xBips = 0;
          p->x     = 0;
        }
       p->number = i;
       p->next   = q;
       p->bInf   = NULL;
       p->back   = NULL;
       p->hash   = 0;
       q         = p;
     }
    p->next->next->next = p;
    t->nodep[i]         = p;
   }

  t->likelihood  = unlikely;
  t->start       = NULL;
  t->ntips       = 0;
  t->nextnode    = 0;

  for (i = 0; i < NUM_BRANCHES; ++ i) t->partitionSmoothed[i] = FALSE;

  t->bitVectors = NULL;
  t->vLength    = 0;
  t->h          = NULL;

  /* TODO: Fix hash type */
  t->nameHash   = pllHashInit (10 * t->mxtips);

  
}

/** @brief Construct a PLL tree structure from a parsed newick tree

    Allocate and initialize a PLL tree structure from a parsed newick
    tree. 

    @param stack
      A stack structure containing the parsed newick tree

    @param nodes
      Total number of nodes of the parsed tree

    @param tips
      Number of tips of the parsed tree

    @todo
      Construct a new data structure called pllNewickTree which will
      consist of three elements: stack, nodes, tips. Also code a proper
      SetupTree function.

    @return
      In case of success returns a pointer to the constructed PLL tree,
      otherwise \b NULL.
*/
tree *
pllTreeCreateNewick (struct pllStack * stack, int nodes, int tips)
{
  tree * t;
  struct pllStack * nodeStack = NULL;
  struct pllStack * head;
  struct item_t * item;
  int i, j, k;
  
/*
  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     partitions->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));
     partitions->partitionData[i]->partitionContribution = -1.0;
     partitions->partitionData[i]->partitionLH           =  0.0;
     partitions->partitionData[i]->fracchange            =  1.0;
   }
*/
  
  t = (tree *) rax_calloc (1, sizeof (tree));
  assert (t);

  pllTreeInitDefaults (t, nodes, tips);

  i = tips + 1;
  j = 1;
  nodeptr v;
  
  
  for (head = stack; head; head = head->next)
  {
    item = (struct item_t *) head->item;
    if (!nodeStack)
     {
       pllStackPush (&nodeStack, t->nodep[i]);
       pllStackPush (&nodeStack, t->nodep[i]->next);
       pllStackPush (&nodeStack, t->nodep[i]->next->next);
       ++i;
     }
    else
     {
       v = (nodeptr) pllStackPop (&nodeStack);
       if (item->rank)  /* internal node */
        {
          v->back           = t->nodep[i];
          t->nodep[i]->back = v; //t->nodep[v->number]
          pllStackPush (&nodeStack, t->nodep[i]->next);
          pllStackPush (&nodeStack, t->nodep[i]->next->next);
          double z = exp((-1 * atof(item->branch))/t->fracchange);
          if(z < zmin) z = zmin;
          if(z > zmax) z = zmax;
          for (k = 0; k < NUM_BRANCHES; ++ k)
             v->z[k] = t->nodep[i]->z[k] = z;

          ++ i;
        }
       else             /* leaf */
        {
          v->back           = t->nodep[j];
          t->nodep[j]->back = v; //t->nodep[v->number];

          double z = exp((-1 * atof(item->branch))/t->fracchange);
          if(z < zmin) z = zmin;
          if(z > zmax) z = zmax;
          for (k = 0; k < NUM_BRANCHES; ++ k)
            v->z[k] = t->nodep[j]->z[k] = z;
            
          //t->nameList[j] = strdup (item->name);
          t->nameList[j] = (char *) rax_malloc ((strlen (item->name) + 1) * sizeof (char));
          strcpy (t->nameList[j], item->name);
          pllHashAdd (t->nameHash, t->nameList[j], (void *) (t->nodep[j]));
          ++ j;
        }
     }
  }
  
  t->start = t->nodep[1];
  
  pllStackClear (&nodeStack);
  return (t);
}


/** @brief Parse a newick tree file
  
    Parse a newick file and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children (rank) and depth.

    @param filename
      Filename containing the newick tree

    @param tree
      Stack structure where the tree will be parsed in

    @param nodes
      Number of nodes will be returned in this variable

    @param leaves
      Number of leaves will be returned in this variable

    @return
      Returns \b 1 in case of success, otherwise \b 0
*/
int
pllNewickParseFile (const char * filename, struct pllStack ** tree, int * nodes, int * leaves)
{
  int n, rc;
  char * rawdata;

  rawdata = readFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  printf ("%s\n\n", rawdata);

  rc = pllNewickParseString (rawdata, tree, nodes, leaves);

  rax_free (rawdata);
  return (rc);
}

