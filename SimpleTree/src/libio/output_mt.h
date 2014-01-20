#ifndef OUTPUT_MT_H

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

/* convert to lgal */
typedef struct clgal_aux_data
{
  hid_t globalRefID;
  hid_t FirstFOF;
  hid_t NextFOF;
  uint32_t nprogs;
  hid_t *proglist;
  struct Lgalaxy_halo_data lgal_halo_data;
  float redshift;
} clgal_aux_data_t;

typedef struct clgal_aux_data_wrapper
{
  uint32_t already_read; //boolean toggle: 1 if readed
  float redshift;
  uint32_t snapid;
  uint32_t domainid;
  hid_t nHalos;
  clgal_aux_data_t *lgal_aux_halos;
} clgal_aux_data_wrapper_t;

extern void sussexbigrun_dm_outputs(m_halo_wrapper_t* haloB, char* outputfolder, int domainid);
extern void internalaux_outputs(m_halo_wrapper_t* haloB, char* outputfolder, int domainid);
extern void generate_lgal_output(char* outputfolder, float *snaplist,int nSnaps, int totaldomains);
#define OUTPUT_MT_H
#endif
