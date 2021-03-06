#ifndef INC_TREEFORMAT_H
#define INC_TREEFORMAT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "libmemmgr/memmgr.h"

typedef uint64_t ptid_t;
typedef uint64_t hid_t;

typedef struct key_sort
{
  ptid_t ID;
  ptid_t order;
} key_sort_t;

typedef struct m_particle
{
  ptid_t ID;
  hid_t haloID;
} m_particle_t;

typedef struct m_merit
{
  hid_t haloID;
  ptid_t NsharedPIDs;
  float Mvir;
  double merit_delucia2007;
  double merit_knollman2009;
  double merit_NFW;
} merit_t;

typedef struct m_merit_embed
{
  double merit_delucia2007;
  double merit_knollman2009;
  double merit_NFW;
  double accumulated_mass;
} merit_embed_t;

typedef struct m_particle_wrapper
{
  ptid_t npart;
  m_particle_t *mparticle;
} m_particle_wrapper_t;


typedef struct particlelist
{
  ptid_t ID;
  float energy;
} particlelist_t;

struct Lgalaxy_halo_data
{
  /* merger tree pointers */
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  /* properties of halo */
  int Len;
  float M_Mean200, M_Crit200, M_TopHat;
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  long long MostBoundID;

  /* original position in subfind output */
  int SnapNum;
  int FileNr;
  int SubhaloIndex;
  float SubHalfMass;
};

typedef struct m_halo 
{
  hid_t ID;
  hid_t oriID;
  hid_t refID;
  hid_t globalRefID;
  float redshift;
#ifdef CUBEP3M
  uint32_t domainID;
#endif
  uint32_t npart;
  float Mvir;
  float Rvir;
  float Xc;
  float Yc;
  float Zc;
  float VXc;
  float VYc;
  float VZc;
  float SpinX;
  float SpinY;
  float SpinZ;
  float Vmax;
  float sigV;
  particlelist_t *Particles;
  /* ####################### */
  hid_t main_progenitor;
  hid_t next_progenitor;
  hid_t descendant;
  hid_t host_halo;
  uint32_t nprogs;
  hid_t *proglist;
  merit_t merit_embed;
  /* ####################### */
  hid_t FirstDownHalo;
  hid_t NextHalo;
  hid_t UpHalo;
  /* ####################### */
  double dm_dt;
  int used;
} m_halo_t;



typedef struct m_halo_wrapper
{
  hid_t nHalos;
  float redshift;
  int snapid;
  m_halo_t *mhalos;
} m_halo_wrapper_t;


extern void allocate_proglist(m_halo_wrapper_t* ptr);
extern void free_m_halo_wrapper(m_halo_wrapper_t* ptr);
#endif
