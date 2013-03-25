#include <stdio.h>
#include "stack.h"
#include "mem_alloc.h"

int pllStackSize (struct pllStack ** stack)
{
  struct pllStack * top;
  int size = 0;
  top = *stack;
 
  while (top)
  {
    ++ size;
    top = top->next;
  }
  
  return (size);
}

int 
pllStackPush (struct pllStack ** head, void * item)
{
  struct pllStack * new;
 
  new = (struct pllStack *) rax_malloc (sizeof (struct pllStack));
  if (!new) return (0);
 
  new->item = item;
  new->next = *head;
  *head     = new;
 
  return (1);
}

void * pllStackPop (struct pllStack ** head)
{
  struct item_t * item;
  struct pllStack * tmp;
  if (!*head) return (NULL);
 
  tmp     = (*head);
  item    = (*head)->item;
  (*head) = (*head)->next;
  rax_free (tmp);
 
  return (item);
}
 
void 
pllStackClear (struct pllStack ** stack)
{
  while (*stack) pllStackPop (stack);
}

