#ifndef INC_TREEFORMAT_H
#define INC_TREEFORMAT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdint.h>
#include "libmemmgr/memmgr.h"
typedef uint64_t ptid_t;
typedef uint64_t hid_t;


typedef struct m_particle
{
  ptid_t ID;
  hid_t haloID;
} m_particle_t;

typedef struct particlelist
{
  ptid_t ID;
} particlelist_t;

typedef struct m_halo 
{
  ptid_t ID;
#ifdef CUBEP3M
  uint32_t domainID;
#endif
  uint32_t snapID;
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
  m_halo_t *mhalos;
} m_halo_wrapper_t;

extern void free_m_halo_wrapper(m_halo_wrapper_t* ptr);
#endif
