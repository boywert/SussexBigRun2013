#ifndef READCONFIG_H
#define READCONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

extern long long unsigned npart_box;
extern double  boxsize;
extern int  domain_per_dim,chunk_per_dim,chunk_mpi;
extern char INPUTDIR[1024];
extern char OUTPUTDIR[1024];
extern char CHUNKDIR[1024];
extern void readconfig();
#endif
