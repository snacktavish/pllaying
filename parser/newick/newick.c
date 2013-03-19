#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lexer.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

#define CONSUME(x)         while (token.class & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

struct item_t
{
  int indent;
  char * name;
  char * branch;
  int leaf;
  int children;
};

struct stack_t
{
  void * item;
  struct stack_t * next;
};

int 
stack_push (struct stack_t ** head, void * item)
{
  struct stack_t * new;

  new = (struct stack_t *) malloc (sizeof (struct stack_t));
  if (!new) return (0);

  new->item = item;
  new->next = *head;
  *head     = new;

  return (1);
}

void * stack_pop (struct stack_t ** head)
{
  struct item_t * item;
  struct stack_t * tmp;
  if (!*head) return (NULL);

  tmp     = (*head);
  item    = (*head)->item;
  (*head) = (*head)->next;
  free (tmp);

  return (item);
}

inline void 
stack_clear (struct stack_t ** stack)
{
  while (*stack) stack_pop (stack);
}

//struct stack_t * stack = NULL;

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

  rawdata = (char *) malloc (((*n)  + 1)* sizeof (char));
  rawdata[*n] = 0;
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static int
parse_newick (char * rawdata, struct stack_t ** stack, int * inp)
{
  struct item_t * item = NULL;
  int item_active = 0;
  struct ltoken_t token;
  int input;
  struct ltoken_t prev_token;
  int nop = 0;
  int indent = 0;

  printf ("%s\n", rawdata);

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
        ++indent;
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
        if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t)); // possibly not nec
        if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        item->indent = indent;
        stack_push (stack, item);
        item_active  = 1;       /* active = 1 */
        item = NULL;
        --indent;
        break;

       case LEX_STRING:
       //printf ("LEX_STRING\n");
        if (prev_token.class != LEX_OPAREN &&
            prev_token.class != LEX_CPAREN &&
            prev_token.class != LEX_UNKNOWN &&
            prev_token.class != LEX_COMMA) return (0);
        if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t));
        item->name = strndup (token.lexeme, token.len);
        item_active = 1;
        item->indent = indent;
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
        if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t));
        if (prev_token.class == LEX_COLON)
         {
           item->branch = strndup (token.lexeme, token.len);
         }
        else
         {
           if (prev_token.class == LEX_COMMA  ||
               prev_token.class == LEX_OPAREN ||
               prev_token.class == LEX_UNKNOWN) item->leaf = 1;
           //if (prev_token.class != LEX_UNKNOWN) ++ indent;
           item->name = strndup (token.lexeme, token.len);
         }
        item_active = 1;
        item->indent = indent;
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
        if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t)); // possibly not nece
        if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        stack_push (stack, item);
        item_active  = 0;
        item = NULL;
        break;

       case LEX_SEMICOLON:
        //printf ("LEX_SEMICOLON\n");
        /* push to the stack */
        if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t));
        if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
        if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        stack_push (stack, item);
        item_active  = 0;
        item = NULL;
        break;
     }
    NEXT_TOKEN
    CONSUME(LEX_WHITESPACE | LEX_NEWLINE);
  }
  if (item_active)
   {
     if (!item) item = (struct item_t *) calloc (1, sizeof (struct item_t));
     if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
     if (item->branch == NULL) item->branch = strdup ("0.000000"); 
     stack_push (stack, item);
     item_active  = 0;
   }

  if (nop) return (0);
  return (1);
}

void stack_dump(struct stack_t ** stack)
{
  struct item_t * item;
  int i;

  while ((item = (struct item_t *)stack_pop(stack)))
   {
     for (i = 0; i < item->indent; ++ i)
       printf ("|\t");
     printf ("%s:%s %d", item->name, item->branch, item->children);
     if (item->leaf) printf (" *\n"); else printf ("\n");
     free (item->name);
     free (item->branch);
     free (item);
   }
}

static void
assign_ranks (struct stack_t * stack, int * nodes, int * leaves)
{
  struct stack_t * head;
  struct item_t * item, * tmp;
  struct stack_t * numbers = NULL;
  int children;
  int indent;

  *nodes = *leaves = 0;


  head = stack;
  while (head)
  {
    assert (head->item);
    item = (struct item_t *) head->item;
    
    if (item->leaf)  ++ (*leaves);

    if (numbers)
     {
       tmp = (struct item_t *)numbers->item;
       children = 0;
       while (item->indent < tmp->indent)
        {
          children = 1;
          indent = tmp->indent;
          stack_pop (&numbers);
          tmp = numbers->item;
          while (tmp->indent == indent)
           {
             ++ children;
             stack_pop (&numbers);
             tmp = (struct item_t *)numbers->item;
           }
          tmp->children += children;
        }
     }
    
    ++ (*nodes);
    head = head->next;

    if (item->leaf)
     {
       if (!numbers) return;

       children = 1;
       tmp = numbers->item;
       while (tmp->indent == item->indent)
        {
          ++ children;
          stack_pop (&numbers);
          assert (numbers);
          tmp = (struct item_t *)numbers->item;
        }
       tmp->children += children;
     }
    else
     {
       stack_push (&numbers, item);
     }
  }
  
  while (numbers->item != stack->item)
  {
    item = (struct item_t *)stack_pop (&numbers);
    tmp  = (struct item_t *) numbers->item;
    children = 1;

    while (tmp->indent == item->indent)
     {
       ++ children;
       item = (struct item_t *)stack_pop (&numbers);
       tmp  = (struct item_t *) numbers->item;
     }
    tmp->children += children;
    children = 0;
  }
 assert (numbers->item == stack->item);
 
 stack_clear (&numbers);
}

int
pllValidateNewick (struct stack_t * stack, int nodes, int leaves)
{
  struct stack_t * head;
  struct item_t * item;
 
 item = stack->item;
 if (item->children != 2 && item->children != 3) return (0);
 head = stack->next;
 while (head)
 {
   item = head->item;
   if (item->children != 2 &&  item->children != 0) 
    {
      return (0);
    }
   head = head->next;
 }
 
 item = stack->item;

 if (item->children == 2) return (nodes == 2 * leaves -1);

 return ((nodes == 2 * leaves - 2) && nodes != 4);
}

int
pllNewickParse (const char * filename, struct stack_t ** tree, int * nodes, int * leaves)
{
  int n, input, rc;
  char * rawdata;

  rawdata = readFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }


  init_lexan (rawdata, n);
  input = get_next_symbol();

  rc = parse_newick (rawdata, tree, &input);

//  assign_ranks (stack, nodes, leaves);
  assign_ranks (*tree, nodes, leaves);

  free (rawdata);
  return (rc);
}


int main (int argc, char * argv[])
{
  int nodes, leaves;
  struct stack_t * tree = NULL;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s FILENAME\n", argv[0]);
     return (EXIT_FAILURE);
   }


  if (pllNewickParse (argv[1], &tree, &nodes, &leaves))
   {
     printf ("Parsing successful...\n\n");

     //if (pllValidateNewick (stack, nodes, leaves))
     if (pllValidateNewick (tree, nodes, leaves))
      {
        printf ("Valid phylogenetic tree\n");
      }
     else
       printf ("Not a valid phylogenetic tree\n");

     stack_dump(&tree);
   }
  else
    printf ("Error while parsing newick tree...\n");

  return (EXIT_SUCCESS);
}

