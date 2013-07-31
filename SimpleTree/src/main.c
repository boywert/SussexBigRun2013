#include "main.h"

int main(int argc,char **argv)

{   
  char folder[1024],outputfolder[1024],snaplistFile[1024];
  m_halo_wrapper_t *halocatA,*halocatB;
  char memmgr_buff[memmgr_max_str];
  double dt,snap1,snap2;
  hid_t ihalo;
  int i,j,k,l;
  int domain_per_dim;
  double boxsize;
  float snaplist[1024];
  FILE *fp;
  
  initialise_MPI(&argc, &argv);
  
  init_memmgr();
  sprintf(folder,"/ccc/cont005/home/ra1089/schneida/scratch/AHF_haloes/cubepm_130315_6_1728_47Mpc_ext2/results");
  sprintf(outputfolder,"/ccc/cont005/home/ra1089/srisawac/scratch/cubepm_130315_6_1728_47Mpc_ext2");
  sprintf(snaplistFile,"halofinds");
  fp = fopen(snaplistFile,"r");
  if(fp == NULL) 
    {
      printf("Cannot open %s\n Exiting...\n",snaplistFile);
    }
  for(i=0;i<1024;i++)
    {
      snaplist[i] = -1.;
    }
  i=0;
  while(fscanf(fp,"%g",&(snaplist[i])) > 0)
    {
      printf("snap:%d  %f\n",i,snaplist[i]);
      i++;
    }
  fclose(fp);
  snap1 = snaplist[20];
  snap2 = snaplist[21];
  sprintf(memmgr_buff,"Halo wrapper");
  //sprintf(folder,"/mnt/lustre/scratch/cs390/testcurie");  
  //halocat = sussexbigrun_load_halo_catalogue_binary(folder,6.000,10*10*10);
  halocatA = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  halocatB = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
  dt = get_delta_t_in_hubble_unit(snap2,snap1);
  for(l=0;l<6*6*6;l++)
    {
      if(mpi_rank==l)
	{
	  halocatA[0] = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,snap1,l);
	  halocatB[0] = sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(folder, snap2, l, 6, 47000.0/6, speed_of_light*dt*max_part_speed_in_c);
	}
    }
  //exit(0);
  make_link_AB(&(halocatA[0]),&(halocatB[0]), dt*kpc2m);
  MPI_Barrier(MPI_COMM_WORLD);

  free_m_halo_wrapper(halocatA);

  sussexbigrun_dm_outputs(&(halocatB[0]),outputfolder);
  
  free_m_halo_wrapper(halocatB);
  //memmgr_printdetails();
  /* for(ihalo=0;ihalo<halocatA[0].nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",halocatA[0].mhalos[ihalo].ID,halocatA[0].mhalos[ihalo].Mvir); */
  /*   } */
  /* exit(0); */

  /* printf(" delta t = %lf, ligt = %lf\n", dt,dt*speed_of_light*max_part_speed_in_c); */
  //(void) MPI_distribute_remove_duplicate_with_structure_level(halocatA,0);
  finalise_MPI();
  return 0;
}
  
