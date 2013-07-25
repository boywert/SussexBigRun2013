#ifndef INC_READSUSSEXBIGRUN_H
#define INC_READSUSSEXBIGRUN_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include "common_io.h"
#include "../libhash/hash.h"
#include "../libmemmgr/memmgr.h"
#include "../treeformat.h"
#include "../common.h"

#define SUSSEXBIGRUN


double MassCut; // haloes below this mass will not be written to ASCII

//===================================================================
// structures
//===================================================================


typedef struct halo {
  uint64_t ID;
  uint64_t hostHalo;
  uint32_t numSubStruct;
  float   Mvir;
  uint32_t npart;
  float   Xc;
  float   Yc;
  float   Zc;
  float   VXc;
  float   VYc;
  float   VZc;
  float   Rvir;
  float   Rmax;
  float   r2;
  float   mbp_offset;
  float   com_offset;
  float   Vmax;
  float   v_esc;
  float   sigV;
  float   lambda;
  float   lambdaE;
  float   Lx;
  float   Ly;
  float   Lz;
  float   b;
  float   c;
  float   Eax;
  float   Eay;
  float   Eaz;
  float   Ebx;
  float   Eby;
  float   Ebz;
  float   Ecx;
  float   Ecy;
  float   Ecz;
  float   ovdens;
  uint32_t nbins;
  float   fMhires;
  float   Ekin;
  float   Epot;
  float   SurfP;
  float   Phi0;
  float   cNFW;
} halo_t;
#define halo_t_size 180


typedef struct halo_profile {
  uint32_t      nbins;
  uint32_t      numColumns;
  float        *r;
  uint32_t      *npart;
  float        *M_in_r;
  float        *ovdens;
  float        *dens;
  float        *vcirc;
  float        *vesc;
  float        *sigv;
  float        *Lx;
  float        *Ly;
  float        *Lz;
  float        *b;
  float        *c;
  float        *Eax;
  float        *Eay;
  float        *Eaz;
  float        *Ebx;
  float        *Eby;
  float        *Ebz;
  float        *Ecx;
  float        *Ecy;
  float        *Ecz;
  float        *Ekin;
  float        *Epot;
  float        *M_gas;
  float        *M_star;
  float        *u_gas;
  float        *Z_gas_sh;
  float        *Z_star_sh;
} halo_profile_t;

struct particle_buffer
{
  uint64_t ID;
  float energy;
};
//===================================================================
// PROTOTYPES
//===================================================================

/* extern void     alloc_profile                     (uint64_t, halo_profile_t *); */
/* extern void     free_profile                      (halo_profile_t *); */
extern m_halo_wrapper_t  sussexbigrun_load_halo_catalogue_binary(char *folder, float redshift, int tot_domain );
extern m_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary_single_domain(char *folder, float redshift, int domain );

extern m_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(char *folder, float redshift, int domain, int domain_per_dim, double domain_width, double dx);
extern m_halo_wrapper_t sussexbigrun_read_AHF_binary(FILE *fphalo, FILE *fppart, int domain, m_halo_wrapper_t mhalo);
extern m_halo_wrapper_t sussexbigrun_filterhalos_and_particles(m_halo_wrapper_t mhalo);
extern void sussexbigrun_makestruct_tree(m_halo_wrapper_t mhalo);
#endif
