#ifndef __pll_STACK__
#define __pll_STACK__

struct pllStack
{
  void * item;
  struct pllStack * next;
};

void  pllStackClear (struct pllStack ** stack);
void * pllStackPop (struct pllStack ** head);
int pllStackPush (struct pllStack ** head, void * item);
int pllStackSize (struct pllStack ** stack);

#endif
