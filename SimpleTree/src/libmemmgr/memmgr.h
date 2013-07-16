#ifndef INC_MEMMGR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdint.h>

#define memmgr_max_str 1024

struct memmgr_obj_struct
{
  char name[memmgr_max_str];
  uint64_t size;
};
struct memmgr_obj_wrapper
{
  uint32_t n_obj;
  uint64_t sum_memmgr_obj_size;
  struct memmgr_obj_struct *memmgr_obj;
};

extern struct memmgr_obj_wrapper memngr_central;

extern void init_memmgr();
extern void *memmgr_malloc(size_t n, char name[memmgr_max_str] );
extern void *memmgr_realloc(void *ptr, size_t new, size_t old, char name[memmgr_max_str] );
extern char *memmgr_printsize(uint64_t size);
extern void memmgr_free(void *ptr, size_t n, char name[memmgr_max_str] );
extern void memmgr_printdetails();
#define INC_MEMMGR_H
#endif
