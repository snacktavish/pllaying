#ifndef __pll_QUEUE__
#define __pll_QUEUE__

struct pllQueueItem
{  
  void * item;
  struct pllQueueItem * next;
}; 
   
struct pllQueue
{  
  struct pllQueueItem * head;
  struct pllQueueItem * tail;
}; 

int pllQueueInit (struct pllQueue ** q);
int pllQueueSize (struct pllQueue * q);
int pllQueueRemove (struct pllQueue * q, void ** item);
int pllQueueAppend (struct pllQueue * q, void * item);
#endif
