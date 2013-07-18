#include "MPI_distribute.h"

void MPI_distribute_remove_duplicate_with_structure_level(m_halo_wrapper_t* mhalo,int level)
{
  /* Level 0 only for now */
  hid_t ihalo,ref;
  qsort(mhalo[0].mhalos,mhalo[0].nHalos, sizeof(m_halo_t), compare_m_halo_t_by_Mvir_reverse);
  ref = mhalo[0].mhalos[0].ID;
  for(ihalo=1;ihalo<mhalo[0].nHalos;ihalo++)
    {
      if(ref == mhalo[0].mhalos[ihalo].ID)
	{
	  printf("duplicated %ld \nCannot perform binary search\n",mhalo[0].mhalos[ihalo].ID);
	}
      ref = mhalo[0].mhalos[ihalo].ID;
      /* printf("%ld => %f\n",mhalo[0].mhalos[ihalo].ID,mhalo[0].mhalos[ihalo].Mvir); */
    }
}
