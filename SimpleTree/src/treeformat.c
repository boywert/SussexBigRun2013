#include "treeformat.h"

void allocate_proglist(m_halo_wrapper_t* ptr)
{
  hid_t i,j;
  for(j=0;j<1;j++)
    {
      for(i=0;i<ptr[j].nHalos;i++)
	{
	  ptr[j].mhalos[i].proglist = malloc(0);
	}

    }
}
void free_m_halo_wrapper(m_halo_wrapper_t* ptr)
{
  hid_t i,j;
  char buff[memmgr_max_str];
  for(j=0;j<1;j++)
    {
      sprintf(buff,"Particle: Halo Array");
      for(i=0;i<ptr[j].nHalos;i++)
	{
	  printf("free particle[%llu]\n",i);
	  memmgr_free(ptr[j].mhalos[i].Particles,ptr[j].mhalos[i].npart*sizeof(particlelist_t),buff);
	  printf("free proglist[%llu]\n",i);
	  free(ptr[j].mhalos[i].proglist);
	}
      sprintf(buff,"Halo Array");
      memmgr_free(ptr[j].mhalos,ptr[j].nHalos*sizeof(m_halo_t),buff);
    }
  sprintf(buff,"Halo wrapper");
  memmgr_free(ptr,sizeof(m_halo_wrapper_t),buff);
}
