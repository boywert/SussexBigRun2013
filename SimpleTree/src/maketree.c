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
  int start_snap;
  float snaplist[1024];
  FILE *fp;
  float redshiftused[2];

  /* [Boyd] initialise MPI and some aux variables needed */
  initialise_MPI(&argc, &argv);
  
  /* [Boyd] initialise memory menager (don't need to use this actually, but some parts of the codes 
     are still using this method to track the memory used) */
  init_memmgr();

  /* [Boyd] Read config files == TODO need to specify config file later */
  readconfig();

  sprintf(folder, param_INPUTDIR);
  sprintf(outputfolder,param_OUTPUTDIR);

  /* [Boyd] TODO = need to specify this in config file */
  sprintf(snaplistFile,"halofinds");
  
  /* [Boyd] Use rank 0 to read the config/snapshot files and broadcast */
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
	  i++;
	}
      fclose(fp);
      tot_Snap = i;
    }
  /* Broadcast snapshot information to other nodes */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(snaplist, 1024, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&tot_Snap, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  /* Read status file == TODO> specify this in config file */
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

  /* Main part for linking snapshots */
  for(i=start_snap;i<=tot_Snap;i++)
    {
      redshiftused[0] = snaplist[i-1];
      redshiftused[1] = snaplist[i];
      sprintf(memmgr_buff,"Halo wrapper");

      /* calculate time difference */
      dt = get_delta_t_in_hubble_unit(redshiftused[1],redshiftused[0]);

      if(mpi_rank==0) printf("Making link AB: %3.3f=>%3.3f step %d\n",redshiftused[0],redshiftused[1],i);

      /* update status file */
      if(mpi_rank==0)
	{
	  sprintf(command,"echo %d > status",i);
	  system(command);
	}

      /* Use all nodes to link particles in their specified domains */
      for(l=0;l<pow3(param_domain_per_dim);l++)
	{
	  if(l%mpi_nodes == mpi_rank)
	    {
	      //printf("\treading domain %d by rank:%d\n",l,mpi_rank);
	      //printf("\tNode %d is making link AB: %3.3f=>%3.3f in domain %d\n",mpi_rank,snap1,snap2,l);

	      /* Allocate memory for halo catalogues*/
	      halocatA = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);
	      halocatB = memmgr_malloc(1*sizeof(m_halo_wrapper_t),memmgr_buff);	      
	      halocatB[0].snapid = i;
	      halocatA[0].snapid = i-1;
	      halocatB[0] = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,redshiftused[1],l);
	      for(ihalo=0;ihalo<halocatB[0].nHalos;ihalo++)
		{
		  printf("in main before readA: %llu npart:%llu\n",halocatB[0].mhalos[ihalo].ID,halocatB[0].mhalos[ihalo].npart);
		}
	      halocatA[0] = sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(folder, redshiftused[0], l, param_domain_per_dim, param_boxsize/param_domain_per_dim, speed_of_light*dt*max_part_speed_in_c);
	      for(ihalo=0;ihalo<halocatB[0].nHalos;ihalo++)
		{
		  printf("in main after readA: %llu npart:%llu\n",halocatB[0].mhalos[ihalo].ID,halocatB[0].mhalos[ihalo].npart);
		}
	      //make_link_AB(&(halocatA[0]),&(halocatB[0]), dt*kpc2m);

	      free_m_halo_wrapper(halocatA);

#ifdef OUTPUTDMDT
	      if(mpi_rank==0) printf("Saving dM/dt ASCII outputs z = %3.3f\n",halocatB[0].redshift);
	      sussexbigrun_dm_outputs(&(halocatB[0]),outputfolder,l);
#endif
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

  finalise_MPI();
  return 0;
}
  
