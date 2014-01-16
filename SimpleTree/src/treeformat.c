#include "treeformat.h"

void free_m_halo_wrapper(m_halo_wrapper_t* ptr)
{
  hid_t i,j;
  char buff[memmgr_max_str];
  for(j=0;j<1;j++)
    {
      
      sprintf(buff,"Particle: Halo Array");
      for(i=0;i<ptr[j].nHalos;i++)
	{
	  memmgr_free(ptr[j].mhalos[i].Particles,ptr[j].mhalos[i].npart*sizeof(particlelist_t),buff);
	  malloc(ptr[j].mhalos[i].proglist,0);
	}
      sprintf(buff,"Halo Array");
      memmgr_free(ptr[j].mhalos,ptr[j].nHalos*sizeof(m_halo_t),buff);
)
    }
  sprintf(buff,"Halo wrapper");
  memmgr_free(ptr,sizeof(m_halo_wrapper_t),buff);
}
