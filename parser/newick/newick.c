#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "newick.h"

static int
parse_newick (struct pllStack ** stack, int * inp)
{
  struct item_t * item = NULL;
  int item_active = 0;
  pllLexToken token;
  int input;
  pllLexToken prev_token;
  int nop = 0;          /* number of open parentheses */
  int depth = 0;

  prev_token.tokenType = PLL_TOKEN_UNKNOWN;

  input = *inp;

  NEXT_TOKEN
  
  while (token.tokenType != PLL_TOKEN_EOF && token.tokenType != PLL_TOKEN_UNKNOWN)
  {
    switch (token.tokenType)
     {
       case PLL_TOKEN_OPAREN:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_OPAREN\n");
#endif
        ++nop;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        ++depth;
        break;

       case PLL_TOKEN_CPAREN:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_CPAREN\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN  &&
            prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
            prev_token.tokenType != PLL_TOKEN_STRING  &&
            prev_token.tokenType != PLL_TOKEN_NUMBER  &&
            prev_token.tokenType != PLL_TOKEN_FLOAT) return (0);

        if (!nop) return (0);
        --nop;
        memcpy (&prev_token, &token, sizeof (pllLexToken));

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

       case PLL_TOKEN_STRING:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_STRING      %.*s\n", token.len, token.lexeme);
#endif
        if (prev_token.tokenType != PLL_TOKEN_OPAREN &&
            prev_token.tokenType != PLL_TOKEN_CPAREN &&
            prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
            prev_token.tokenType != PLL_TOKEN_COMMA) return (0);
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
        //item->name = strndup (token.lexeme, token.len);
        item->name = (char *) rax_malloc ((token.len + 1) * sizeof (char));
        strncpy (item->name, token.lexeme, token.len);
        item->name[token.len] = 0;

        item_active = 1;
        item->depth = depth;
        if (prev_token.tokenType == PLL_TOKEN_COMMA  ||
            prev_token.tokenType == PLL_TOKEN_OPAREN ||
            prev_token.tokenType == PLL_TOKEN_UNKNOWN) item->leaf = 1;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_FLOAT:
       case PLL_TOKEN_NUMBER:
#ifdef PLLDEBUG
       if (token.tokenType == PLL_TOKEN_FLOAT) printf ("PLL_TOKEN_FLOAT\n"); else printf ("PLL_TOKEN_NUMBER\n");
#endif
         if  (prev_token.tokenType != PLL_TOKEN_OPAREN &&
              prev_token.tokenType != PLL_TOKEN_CPAREN &&
              prev_token.tokenType != PLL_TOKEN_COLON  &&
              prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
              prev_token.tokenType != PLL_TOKEN_COMMA) return (0);
        if (!item) item = (struct item_t *) rax_calloc (1, sizeof (struct item_t));
        if (prev_token.tokenType == PLL_TOKEN_COLON)
         {
           //item->branch = strndup (token.lexeme, token.len);
           item->branch = (char *) rax_malloc ((token.len + 1) * sizeof (char));
           strncpy (item->branch, token.lexeme, token.len);
           item->branch[token.len] = 0;
         }
        else
         {
           if (prev_token.tokenType == PLL_TOKEN_COMMA  ||
               prev_token.tokenType == PLL_TOKEN_OPAREN ||
               prev_token.tokenType == PLL_TOKEN_UNKNOWN) item->leaf = 1;
           //if (prev_token.tokenType != PLL_TOKEN_UNKNOWN) ++ indent;
           //item->name = strndup (token.lexeme, token.len);
           item->name = (char *) rax_malloc ((token.len + 1) * sizeof (char));
           strncpy (item->name, token.lexeme, token.len);
           item->name[token.len] = 0;
         }
        item_active = 1;
        item->depth = depth;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_COLON:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_COLON\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN &&
            prev_token.tokenType != PLL_TOKEN_STRING &&
            prev_token.tokenType != PLL_TOKEN_FLOAT  &&
            prev_token.tokenType != PLL_TOKEN_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_COMMA:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_COMMA\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN &&
             prev_token.tokenType != PLL_TOKEN_STRING &&
             prev_token.tokenType != PLL_TOKEN_FLOAT && 
             prev_token.tokenType != PLL_TOKEN_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        
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

       case PLL_TOKEN_SEMICOLON:
#ifdef PLLDEBUG
        printf ("PLL_TOKEN_SEMICOLON\n");
#endif
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
       default:
#ifdef PLLDEBUG
         printf ("Unknown token: %d\n", token.tokenType);
#endif
       // TODO: Finish this part and add error codes
        break;
     }
    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE);
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

  if (nop || token.tokenType == PLL_TOKEN_UNKNOWN) 
   {
     return (0);
   }

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
      Newick tree wrapper structure which contains the stack representation of the parsed newick tree

    @return
      Returns \b 1 in case of success, otherwise \b 0
*/
int
pllValidateNewick (pllNewickTree * t)
{
  struct pllStack * head;
  struct item_t * item;
 
 item = t->tree->item;
 if (item->rank != 2 && item->rank != 3) return (0);
 head = t->tree->next;
 while (head)
 {
   item = head->item;
   if (item->rank != 2 &&  item->rank != 0) 
    {
      return (0);
    }
   head = head->next;
 }
 
 item = t->tree->item;

 if (item->rank == 2) return (t->nodes == 2 * t->tips -1);

 return ((t->nodes == 2 * t->tips - 2) && t->nodes != 4);
}

/** @brief Parse a newick tree string
  
    Parse a newick string and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children and depth. The
    stack structure is finally wrapped in a \a pllNewickTree structure which
    also contains the number of nodes and leaves.

    @param newick
      String containing the newick tree

    @return
      Returns a pointer to the created \a pllNewickTree structure in case of success, otherwise \b NULL
*/
pllNewickTree *
pllNewickParseString (char * newick)
{
  int n, input, rc;
  pllNewickTree * t;
  int nodes, leaves;
  
  t = (pllNewickTree *) calloc (1, sizeof (pllNewickTree));

  n = strlen (newick);

  init_lexan (newick, n);
  input = get_next_symbol();

  rc = parse_newick (&(t->tree), &input);
  if (!rc)
   {
     /* TODO: properly clean t->tree */
     rax_free (t);
     t = NULL;
   }
  else
   {
     assign_ranks (t->tree, &nodes, &leaves);
     t->nodes = nodes;
     t->tips  = leaves;
   }

  return (t);
}

/** @brief Deallocate newick parser stack structure

    Deallocates the newick parser stack structure that represents the parsed tree. It
    also frees all memory allocated by elements of the stack structure.

    @param tree
      The tree stack structure
*/
void pllNewickParseDestroy (pllNewickTree ** t)
{
  struct item_t *  item;

  while ((item = (struct item_t *)pllStackPop (&((*t)->tree))))
   {
     rax_free (item->name);
     rax_free (item->branch);
     rax_free (item);
   }
  rax_free (*t);
  (*t) = NULL;
}

/** @brief Parse a newick tree file
  
    Parse a newick file and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children (rank) and depth. The
    stack structure is finally wrapped in a \a pllNewickTree structure which
    also contains the number of nodes and leaves.

    @param filename
      Filename containing the newick tree

    @return
      Returns a pointer to the created \a pllNewickTree structure in case of success, otherwise \b NULL
*/
pllNewickTree *
pllNewickParseFile (const char * filename)
{
  long n;
  char * rawdata;
  pllNewickTree * t;

  rawdata = pllReadFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  printf ("%s\n\n", rawdata);

  t = pllNewickParseString (rawdata);

  rax_free (rawdata);

  return (t);
}

