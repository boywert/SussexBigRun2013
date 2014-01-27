#ifndef READCONFIG_H
#define READCONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

extern long long unsigned param_npart_box;
extern double  param_boxsize,param_buffer_size,param_fixed_padding;
extern int  param_domain_per_dim,param_chunk_per_dim,param_chunk_mpi,param_domain_per_chunk;
extern char param_INPUTDIR[1024];
extern char param_OUTPUTDIR[1024];
extern char param_CHUNKDIR[1024];
extern char param_CUBEP3MOUT[1024];
extern void readconfig(char* configfile);
#endif
