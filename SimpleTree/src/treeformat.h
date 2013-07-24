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

typedef struct m_particle_wrapper
{
  ptid_t npart;
  m_particle_t *mparticle;
} m_particle_wrapper_t;


typedef struct particlelist
{
  ptid_t ID;
  //float energy;
} particlelist_t;



typedef struct m_halo 
{
  ptid_t ID;
  ptid_t oriID;
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
  particlelist_t *Particles;
  /* ####################### */
  hid_t main_progenitor; 
  hid_t next_progenitor;
  hid_t descendant;
  /* ####################### */
  hid_t main_subhalo;
  hid_t next_subhalo;
  hid_t host_halo;
} m_halo_t;

typedef struct m_halo_wrapper
{
  hid_t nHalos;
  float redshift;
  m_halo_t *mhalos;
} m_halo_wrapper_t;

extern void free_m_halo_wrapper(m_halo_wrapper_t* ptr);
#endif
