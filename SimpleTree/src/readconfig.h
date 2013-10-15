#ifndef READCONFIG_H
#define READCONFIG_H

#include <stdio.h>
#include <string.h>

extern long long unsigned npart_box;
extern double  boxsize;
extern int  domain_per_dim;
extern char INPUTDIR[1024];
extern char OUTPUTDIR[1024];

extern void readconfig();
#endif
