#ifndef INC_MAIN_H
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "treeformat.h"
#include "libmemmgr/memmgr.h"
#include "libio/readsussexbigrun.h"
#include "libio/common_io.h"
#include "libio/output_dmdt.h"
#include "cosmology.h"
#include "MPI_distribute.h"
#include <math.h>
#include "mpi.h"
#include "descendant.h"
#include "common.h"


/* Force trigger some options */

#ifdef OUTPUTDMDT
#define TOPLEVELONLY
#endif



#define INC_MAIN_H


#endif
