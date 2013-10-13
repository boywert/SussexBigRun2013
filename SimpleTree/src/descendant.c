#include "descendant.h"


void make_link_AB(m_halo_wrapper_t* haloA, m_halo_wrapper_t* haloB, double dt)
{
  m_particle_wrapper_t *tmppart;
  ptid_t ipart,countpart,ref,curpart;
  hid_t ihalo,jhalo,ihid,max_id;
  double max_Mvir;
  uint64_t old,new;
  merit_t *merit,*merit_prog;
  char memmgr_buff[memmgr_max_str];
  //printf("make link AB\n");
  sprintf(memmgr_buff,"Particle Wrapper: Hash");
  tmppart = memmgr_malloc(sizeof(m_particle_wrapper_t),memmgr_buff);
  tmppart[0].npart = 0;
  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  tmppart[0].mparticle = memmgr_malloc(0,memmgr_buff);
  qsort(haloB->mhalos,haloB->nHalos, sizeof(m_halo_t),compare_m_halo_t_by_Mvir);
  countpart = 0;
  //printf("Start loop for haloB\n");
  for(ihalo=0;ihalo < haloB->nHalos; ihalo++)
    {
      haloB->mhalos[ihalo].ID = ihalo;
      haloB->mhalos[ihalo].main_progenitor = NULLPOINT;
      tmppart[0].npart += haloB->mhalos[ihalo].npart;
      tmppart[0].mparticle = memmgr_realloc(tmppart[0].mparticle,sizeof(m_particle_t)*tmppart[0].npart,sizeof(m_particle_t)*(tmppart[0].npart-haloB->mhalos[ihalo].npart),memmgr_buff);
      for(ipart=0;ipart<haloB->mhalos[ihalo].npart;ipart++)
  	{
    	  tmppart[0].mparticle[countpart].ID =  haloB->mhalos[ihalo].Particles[ipart].ID;
  	  tmppart[0].mparticle[countpart].haloID = haloB->mhalos[ihalo].ID;
  	  countpart++;
    	}
  
    }
  qsort(tmppart[0].mparticle,tmppart[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID);
  ref = NULLPOINT;
  countpart = 0;
  sprintf(memmgr_buff,"TMP particles: Hash");
  //printf("loop to remove dup\n");
  for(ipart=0;ipart<tmppart[0].npart;ipart++)
    {
      if(tmppart[0].mparticle[ipart].ID == ref)
  	{
  	  tmppart[0].mparticle[ipart].ID = NULLPOINT;
  	}
      else
  	{
	  countpart++;
  	  ref = tmppart[0].mparticle[ipart].ID;
  	}
   
    }
  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  qsort(tmppart[0].mparticle, tmppart[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID);
  old = tmppart[0].npart*sizeof(m_particle_t);
  new = countpart*sizeof(m_particle_t);
  tmppart[0].mparticle = memmgr_realloc(tmppart[0].mparticle,new,old,memmgr_buff);
  tmppart[0].npart = countpart;
  /* for(ipart=0; ipart<tmppart[0].npart; ipart++ ) */
  /*   { */
  /*     printf("ipart: %llu => %llu\n",ipart,tmppart[0].mparticle[ipart].ID); */
  /*   } */
  /* exit(0); */
  merit = malloc(haloB->nHalos*sizeof(merit_t));
  for(ihalo = 0; ihalo < haloA->nHalos; ihalo++)
    {
      for(jhalo=0;jhalo<haloB->nHalos;jhalo++)
	{
	  merit[jhalo].haloID = jhalo;
	  merit[jhalo].merit_delucia2007 = 0.;
	}
      for(ipart=0; ipart < haloA->mhalos[ihalo].npart; ipart++)
	{
	  curpart = haloA->mhalos[ihalo].Particles[ipart].ID;
	  //printf("outside: search for %llu\n",curpart);
	  ihid =  search_m_particle_t_for_ID(curpart,tmppart[0].npart,&(tmppart[0].mparticle[0]) );
	  if(ihid < NULLPOINT) 
	    merit[ihid].merit_delucia2007 += pow((double)ipart,-2./3);
	}
      qsort(merit,haloB->nHalos,sizeof(merit_t),compare_merit_t_by_merit_delucia2007);
      if(merit[haloB->nHalos-1].merit_delucia2007 > 2.5)
	haloA->mhalos[ihalo].descendant = merit[haloB->nHalos-1].haloID;
      else
	haloA->mhalos[ihalo].descendant = NULLPOINT;
      //printf("descendant => %llu: delta M \n",haloA->mhalos[ihalo].descendant, haloA->mhalos[ihalo].Mvir-haloB->mhalos[]);
    }
  free(merit)

  qsort(haloA->mhalos,haloA->nHalos, sizeof(m_halo_t),compare_m_halo_t_by_descendant);
  for(ihalo=0; ihalo < haloA->nHalos; ihalo++)
    {
      haloA->mhalos[ihalo].ID = ihalo;
    }
  ihid = NULLPOINT;
  max_id = NULLPOINT;
  max_Mvir = 0.;
  for(ihalo = 0; ihalo < haloA->nHalos; ihalo++)
    {
      if(haloA->mhalos[ihalo].descendant == NULLPOINT)
	break;
      if(haloA->mhalos[ihalo].descendant == ihid)
	{
	  if(haloA->mhalos[ihalo].Mvir > max_Mvir)
	    {
	      //max_Mvir = MAX(mhalos[ihalo].Mvir,max_Mvir);
	      max_id = ihalo;
	      max_Mvir = haloA->mhalos[ihalo].Mvir;
	    }
	}
      else
	{
	  ihid = haloA->mhalos[ihalo].descendant;
	  if(ihalo > 0)
	    {
	      haloB->mhalos[haloA->mhalos[ihalo-1].descendant].main_progenitor = ihalo;
	      max_Mvir = haloA->mhalos[ihalo].Mvir;
	    }
	}
    }
  for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
    {
      if(haloB->mhalos[ihalo].main_progenitor < NULLPOINT)
  	{
	  haloB->mhalos[ihalo].dm_dt = (haloB->mhalos[ihalo].Mvir-haloA->mhalos[haloB->mhalos[ihalo].main_progenitor].Mvir)/dt;
	  // printf("halo %llu<=%llu dm = %lf\n",ihalo,haloB->mhalos[ihalo].main_progenitor,haloB->mhalos[ihalo].Mvir-haloA->mhalos[haloB->mhalos[ihalo].main_progenitor].Mvir);
  	}
      else
  	{
	  haloB->mhalos[ihalo].dm_dt = haloB->mhalos[ihalo].Mvir/dt;
  	  //printf("halo %llu<=%llu dm = %lf\n",ihalo,NULLPOINT,haloB->mhalos[ihalo].Mvir);
  	}
    }
}


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
