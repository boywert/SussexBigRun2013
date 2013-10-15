#include "main.h"

int main(int argc,char **argv)

{   
  char folder[1024],outputfolder[1024],snaplistFile[1024];
  char command[1024],stringbuff[10];
  m_halo_wrapper_t *halocatA,*halocatB;
  char memmgr_buff[memmgr_max_str];
  double dt,snap1,snap2;
  hid_t ihalo;
  int i,j,k,l,tot_Snap;
  //int domain_per_dim;
  int start_snap;
  //double boxsize;
  float snaplist[1024];
  FILE *fp;
  
  initialise_MPI(&argc, &argv);
  
  init_memmgr();
  //sprintf(folder,"/ccc/cont005/home/ra1089/srisawac/scratch/AHF/cubepm_130315_6_1728_47Mpc_ext2/results");
  //sprintf(outputfolder,"/ccc/cont005/home/ra1089/srisawac/scratch/cubepm_130315_6_1728_47Mpc_ext2");
  readconfig();
  sprintf(folder, INPUTDIR);
  sprintf(outputfolder,OUTPUTDIR);
  sprintf(snaplistFile,"halofinds");
  
  if(mpi_rank==0)
    {
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
      while(fscanf(fp,"%g",&(snaplist[i])) != EOF)
	{
	  //	  printf("snap:%d  %f\n",i,snaplist[i]);
	  i++;
	}
      fclose(fp);
      tot_Snap = i;
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(snaplist, 1024, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&tot_Snap, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(mpi_rank==0)
    {
      fp = fopen("status","r");
      if(fp == NULL)
	start_snap = 1;
      else
	{
	  if (fgets(stringbuff , 10 , fp) != NULL )
	    sscanf(stringbuff,"%d",&start_snap);
	  fclose (fp);
	}
      if(start_snap < 1) start_snap = 1;
      printf("Start making merger trees from Snapshot:%d z=%3.3f\n",start_snap,snaplist[start_snap]);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&start_snap, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=start_snap;i<=tot_Snap;i++)
    {
      snap1 = snaplist[i-1];
      snap2 = snaplist[i];
      sprintf(memmgr_buff,"Halo wrapper");
      //sprintf(folder,"/mnt/lustre/scratch/cs390/testcurie");  
      //halocat = sussexbigrun_load_halo_catalogue_binary(folder,6.000,10*10*10);
      dt = get_delta_t_in_hubble_unit(snap2,snap1);
      if(mpi_rank==0) printf("Making link AB: %3.3f=>%3.3f step %d\n",snap1,snap2,i);
      if(mpi_rank==0)
	{
	  sprintf(command,"echo %d > status",i);
	  system(command);
	}
      for(l=0;l<domain_per_dim*domain_per_dim*domain_per_dim;l++)
	{
	  if(l%mpi_nodes == mpi_rank)
	    {
	      printf("\treading domain %d by rank:%d\n",l,mpi_rank);
	      //printf("\tNode %d is making link AB: %3.3f=>%3.3f in domain %d\n",mpi_rank,snap1,snap2,l);
	      halocatA = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
	      halocatB = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);	      
	      halocatB[0] = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,snap2,l);
	      halocatA[0] = sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(folder, snap1, l, domain_per_dim, boxsize/domain_per_dim, speed_of_light*dt*max_part_speed_in_c);
	      //make_link_AB(&(halocatA[0]),&(halocatB[0]), dt*kpc2m);

	      
	      free_m_halo_wrapper(halocatA);

	      //if(mpi_rank==0) printf("Saving ASCII outputs z = %3.3f\n",halocatB[0].redshift);
	      //sussexbigrun_dm_outputs(&(halocatB[0]),outputfolder,l);
	      free_m_halo_wrapper(halocatB);
	    }
	}
      MPI_Barrier(MPI_COMM_WORLD);
      if(mpi_rank==0)
	{
	  sprintf(command,"echo %d > status",i+1);
	  system(command);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
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
  
