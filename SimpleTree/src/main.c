#include "main.h"

int main(int argc,char **argv)

{   
  char folder[1024];
  m_halo_wrapper_t *halocatA,*halocatB;
  char memmgr_buff[memmgr_max_str];
  double dt;
  hid_t ihalo;
  int i,j,k,l,nodes_return,rank;
  MPI_Init ( &argc, &argv );
  MPI_Comm_size ( MPI_COMM_WORLD, &nodes_return );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  init_memmgr();
  sprintf(folder,"/ccc/cont005/home/ra1089/schneida/scratch/AHF_haloes/cubepm_130315_6_1728_47Mpc_ext2/results");
  sprintf(memmgr_buff,"Halo wrapper");
  //sprintf(folder,"/mnt/lustre/scratch/cs390/testcurie");  
  //halocat = sussexbigrun_load_halo_catalogue_binary(folder,6.000,10*10*10);
  halocatA = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  halocatB = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  dt = get_delta_t_in_hubble_unit(6.354,6.418);
  for(l=0;l<6*6*6;l++)
    {
      if(rank==l)
	{
	  halocatA[0] = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,6.418,l);
	  halocatB[0] = sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(folder, 6.354, l, 6, 47.0/6, speed_of_light*dt*max_part_speed_in_c);
	  //halocatB[0] = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,6.354,l);
	}
    }
  // halocatB[0] = sussexbigrun_load_halo_catalogue_binary(folder,6.418,6);
  /* memmgr_printdetails(); */
  free_m_halo_wrapper(halocatA);
  free_m_halo_wrapper(halocatB);
  //memmgr_printdetails();
  /* for(ihalo=0;ihalo<halocatA[0].nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",halocatA[0].mhalos[ihalo].ID,halocatA[0].mhalos[ihalo].Mvir); */
  /*   } */
  /* exit(0); */

  /* printf(" delta t = %lf, ligt = %lf\n", dt,dt*speed_of_light*max_part_speed_in_c); */
  //(void) MPI_distribute_remove_duplicate_with_structure_level(halocatA,0);
  MPI_Finalize();
  return 0;
}
  
