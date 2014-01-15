#include "descendant.h"

/* Make link A -> B, timewise */
void make_link_AB(m_halo_wrapper_t* haloA, m_halo_wrapper_t* haloB, double dt)
{
  m_particle_wrapper_t *tmppartB;
  ptid_t ipart,countpart,ref,curpart;
  hid_t ihalo,jhalo,ihid,max_id,count_progs,previous_id,iprog,next_id;
  double max_Mvir;
  uint64_t old,new;
  merit_t *merit,*merit_prog;
  char memmgr_buff[memmgr_max_str];
  hid_t loophalo;
  /* [Boyd] Make the catalogue B exclusive table */

  
  sprintf(memmgr_buff,"Particle Wrapper: Hash");
  tmppartB = memmgr_malloc(sizeof(m_particle_wrapper_t),memmgr_buff);
  tmppartB[0].npart = 0;
  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  tmppartB[0].mparticle = memmgr_malloc(0,memmgr_buff);
  if(haloB->nHalos > 1)
    qsort(haloB->mhalos,haloB->nHalos, sizeof(m_halo_t),compare_m_halo_t_by_Mvir);
  countpart = 0;


  for(ihalo=0;ihalo < haloB->nHalos; ihalo++)
    {
      /* Redefine ID to the order by Mvir */
      haloB->mhalos[ihalo].ID = ihalo;
      /* Set main_progenitor to -1 */
      haloB->mhalos[ihalo].main_progenitor = NULLPOINT;
      /* Increase the particle in tmppartB */
      tmppartB[0].npart += haloB->mhalos[ihalo].npart;
      tmppartB[0].mparticle = memmgr_realloc(tmppartB[0].mparticle,sizeof(m_particle_t)*tmppartB[0].npart,sizeof(m_particle_t)*(tmppartB[0].npart-haloB->mhalos[ihalo].npart),memmgr_buff);
      for(ipart=0;ipart<haloB->mhalos[ihalo].npart;ipart++)
  	{
    	  tmppartB[0].mparticle[countpart].ID =  haloB->mhalos[ihalo].Particles[ipart].ID;
  	  tmppartB[0].mparticle[countpart].haloID = haloB->mhalos[ihalo].ID;
  	  countpart++;
    	}
  
    }
  /* Sort particles by PID  */
  qsort(tmppartB[0].mparticle,tmppartB[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID);
  ref = NULLPOINT;
  countpart = 0;
  sprintf(memmgr_buff,"TMP particles: Hash");
  
  /* Remove all duplicated PID, use only the first one in the set (from the lowest Mvir) */
  for(ipart=0;ipart<tmppartB[0].npart;ipart++)
    {
      if(tmppartB[0].mparticle[ipart].ID == ref)
  	{
  	  tmppartB[0].mparticle[ipart].ID = NULLPOINT;
  	}
      else
  	{
  	  countpart++;
  	  ref = tmppartB[0].mparticle[ipart].ID;
  	}
   
    }
  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  qsort(tmppartB[0].mparticle, tmppartB[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID);
  old = tmppartB[0].npart*sizeof(m_particle_t);
  new = countpart*sizeof(m_particle_t);
  tmppartB[0].mparticle = memmgr_realloc(tmppartB[0].mparticle,new,old,memmgr_buff);
  tmppartB[0].npart = countpart;

  /* Make exclusive PIDs lists and correct npart */

  /* Finish making catalogue B exclusive table */

  merit = malloc(haloB->nHalos*sizeof(merit_t));
  for(ihalo = 0; ihalo < haloA->nHalos; ihalo++)
    {
      for(jhalo=0;jhalo<haloB->nHalos;jhalo++)
  	{
  	  merit[jhalo].haloID = jhalo;
  	  merit[jhalo].merit_delucia2007 = 0.;
  	  merit[jhalo].NsharedPIDs = 0;
  	}
      for(ipart=0; ipart < haloA->mhalos[ihalo].npart; ipart++)
  	{
  	  curpart = haloA->mhalos[ihalo].Particles[ipart].ID;
  	  //printf("outside: search for %llu\n",curpart);
  	  ihid =  search_m_particle_t_for_ID(curpart,tmppartB[0].npart,&(tmppartB[0].mparticle[0]) );
  	  if(ihid < NULLPOINT)
  	    {
  	      merit[ihid].merit_delucia2007 += pow((double)ipart,-2./3);
  	      merit[ihid].NsharedPIDs += 1;
  	    }
  	}
       for(jhalo=0;jhalo<haloB->nHalos;jhalo++)
  	{
  	  merit[jhalo].merit_knollman2009 = (double)merit[jhalo].NsharedPIDs*(double)merit[jhalo].NsharedPIDs/(double)haloA->mhalos[ihalo].npart/haloB->mhalos[jhalo].npart;
  	}
      qsort(merit,haloB->nHalos,sizeof(merit_t),compare_merit_t_by_merit_delucia2007);
      if(merit[haloB->nHalos-1].merit_delucia2007 > 2.5)
  	{
  	  haloA->mhalos[ihalo].descendant = merit[haloB->nHalos-1].haloID;
  	  haloA->mhalos[ihalo].merit_embed.merit_delucia2007 = merit[haloB->nHalos-1].merit_delucia2007;
  	  haloA->mhalos[ihalo].merit_embed.merit_knollman2009 = merit[haloB->nHalos-1].merit_knollman2009;
  	}
      else
  	haloA->mhalos[ihalo].descendant = NULLPOINT;
      //printf("descendant => %llu: delta M \n",haloA->mhalos[ihalo].descendant, haloA->mhalos[ihalo].Mvir-haloB->mhalos[]);
    }
  free(merit);


  memmgr_free(tmppartB[0].mparticle,tmppartB[0].npart*sizeof(m_particle_t),"Particle inside wrapper: Hash");
  memmgr_free(tmppartB,sizeof(m_particle_wrapper_t),"Particle Wrapper: Hash");


  qsort(haloA->mhalos,haloA->nHalos, sizeof(m_halo_t),compare_m_halo_t_by_descendant);

  for(ihalo=0; ihalo < haloA->nHalos; ihalo++)
    {
      haloA->mhalos[ihalo].ID = ihalo;
    }


  ihid = NULLPOINT;
  max_id = NULLPOINT;
  max_Mvir = 0.;
  merit_prog = malloc(0);
  for(ihalo = 0; ihalo < haloA->nHalos; ihalo++)
    {
      if(haloA->mhalos[ihalo].descendant == NULLPOINT)
  	break;
      if(haloA->mhalos[ihalo].descendant == ihid)
  	{
  	  haloB->mhalos[ihid].nprogs += 1;
  	  merit_prog = realloc(merit_prog,haloB->mhalos[ihid].nprogs*sizeof(merit_t));
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].haloID = ihalo;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].merit_delucia2007 = haloA->mhalos[ihalo].merit_embed.merit_delucia2007;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].merit_knollman2009 = haloA->mhalos[ihalo].merit_embed.merit_knollman2009;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].Mvir = haloA->mhalos[ihalo].Mvir;
  	}
      else
  	{
  	  if(ihid < NULLPOINT)
  	    {
  	      qsort(merit_prog,haloB->mhalos[ihid].nprogs,sizeof(merit_t),compare_merit_t_by_Mvir);
  	      haloB->mhalos[ihid].main_progenitor = merit_prog[haloB->mhalos[ihid].nprogs-1].haloID;
  	      if(merit_prog,haloB->mhalos[ihid].nprogs > 1)
  		{
		  printf("print halo list\n");
		  for(loophalo=0;loophalo<haloA->nHalos;loophalo++)
		    {
		      printf("TEST=> ID:%llu -> %llu\n",loophalo,haloA->mhalos[loophalo].oriID);
		    }
  		  for(iprog = haloB->mhalos[ihid].nprogs-2; iprog>0; iprog--)
  		    {
  		      haloA->mhalos[merit_prog[iprog+1].haloID].next_progenitor = merit_prog[iprog].haloID;
  		    }
  		  haloA->mhalos[merit_prog[1].haloID].next_progenitor = merit_prog[0].haloID;
  		  haloA->mhalos[merit_prog[0].haloID].next_progenitor = NULLPOINT;
  		}
  	    }
  	  ihid = haloA->mhalos[ihalo].descendant;
  	  haloB->mhalos[ihid].nprogs = 1;
  	  merit_prog = realloc(merit_prog,haloB->mhalos[ihid].nprogs*sizeof(merit_t));
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].haloID = ihalo;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].merit_delucia2007 = haloA->mhalos[ihalo].merit_embed.merit_delucia2007;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].merit_knollman2009 = haloA->mhalos[ihalo].merit_embed.merit_knollman2009;
  	  merit_prog[haloB->mhalos[ihid].nprogs-1].Mvir = haloA->mhalos[ihalo].Mvir;
  	}
    }
  free(merit_prog);

  /* Calculate dM/dt  */
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
  for(ihalo=0;ihalo < haloB->nHalos; ihalo++)
    {
      printf("%llu -> %llu",haloB->mhalos[ihalo].ID,haloB->mhalos[ihalo].main_progenitor);
      next_id = haloB->mhalos[ihalo].main_progenitor;
      /* while(next_id < NULLPOINT) */
      /* 	{ */
      /* 	  //printf(":%llu:%d",haloA->mhalos[next_id].oriID,haloA->mhalos[next_id].npart); */
      /* 	  next_id = haloA->mhalos[next_id].next_progenitor; */
      /* 	  //printf(" -> %llu",next_id); */
      /* 	} */
      printf("\n");
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
