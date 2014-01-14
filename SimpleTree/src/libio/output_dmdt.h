#ifndef OUTPUT_DMDT_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "common_io.h"
#include "../libhash/hash.h"
#include "../libmemmgr/memmgr.h"
#include "../treeformat.h"
#include "../common.h"
#include "mpi.h"
extern void sussexbigrun_dm_outputs( m_halo_wrapper_t* haloB, char* outputfolder, int domainid);

#define OUTPUT_DMDT_H
#endif
