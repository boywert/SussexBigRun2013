#include "descendant.h"

/* m_particle_wrapper_t* build_m_particle_t_for_search(m_halo_t* mhalo, int allow_dup) */
/* { */
/*   m_particle_wrapper_t* partlist; */
/*   hid_t ihalo; */
/*   ptid_t ipart; */
/*   char memmgr_buff[memmgr_max_str]; */
/*   uint64_t old,new; */
/*   sprintf(memmgr_buff,"Particle Wrapper: Hash") */
/*   partlist = memmgr_malloc(sizeof(m_particle_wrapper_t), memmgr_buff); */
/*   partlist.npart = 0; */
/*   sprintf(memmgr_buff,"Particle inside wrapper: Hash"); */
/*   partlist.mparticle = memmgr_malloc(0,memmgr_buff); */
/*   for(ihalo=0; ihalo<mhalo.nHalos; ihalo++) */
/*     { */
/*       for(ipart=0; ipart<mhalo.mhalos[ihalo].npart; ipart++) */
/* 	{ */
/* 	  pid = mhalo.mhalos[ihalo].Particles[ipart].ID; */
/* 	  if pid is in part list */
/* 	    { */
/* 	      if(allow_dup) */
/* 		{ */
		  
/* 		} */
/* 	      else */
/* 		{ */
		  
/* 		} */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      old = partlist.npart*sizeof(m_particle_t); */
/* 	      partlist.npart++; */
/* 	      new = partlist.npart*sizeof(m_particle_t); */
/* 	      partlist.mparticle = memmgr_realloc(particle.mparticle, new, old, memmgr_buff); */
/* 	      partlist.mparticle[partlist.npart-1].ID = pid; */
/* 	      partlist.mparticle[partlist.npart-1].ID = mhalo.mhalos[ihalo].ID; */
/* 	      /\* one iteration of insert/bubble backward sort  *\/ */
/*  	    } */
/* 	} */
/*     } */
/*   return partlist; */
/* } */
