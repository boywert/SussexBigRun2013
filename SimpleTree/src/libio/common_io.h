#ifndef INC_COMMON_IO_H
#define INC_COMMON_IO_H
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

extern int   ReadULong  (FILE *fptr,uint64_t *n,int swap);        
extern int   ReadInt  (FILE *fptr,int32_t *n,int swap);          
extern int   ReadUInt  (FILE *fptr,uint32_t *n,int swap);         
extern int   ReadFloat (FILE *fptr,float *n,int swap);         

#define FALSE 0
#define TRUE 1

#endif
