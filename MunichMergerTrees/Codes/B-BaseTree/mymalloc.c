#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



void *mymalloc(size_t n)
{
  void *p;

  if(!(p = malloc(n)))
    {
      if(n)
	{
	  printf("Failed to allocate memory for %u bytes.\n", (int) n);
	  exit(2);
	}
    }

  return p;
}


void myfree(void *ptr)
{
  free(ptr);
}
