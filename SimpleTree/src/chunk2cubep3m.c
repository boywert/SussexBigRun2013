#include "main.h"

int main(int argc,char **argv)

{   
  char folder[1024],outputfolder[1024],snaplistFile[1024];
  char command[1024],stringbuff[10];
  make_catalogue_halo_wrapper_t *halocatA;
  char memmgr_buff[memmgr_max_str];
  double dt,snap1,snapid1;
  hid_t ihalo;
  int i,j,k,l,tot_Snap;
  //int domain_per_dim;
  int start_snap;
  //double boxsize;
  float snaplist[1024];
  FILE *fp;
  
  initialise_MPI(&argc, &argv);
  
  init_memmgr();
  readconfig();
  sprintf(folder, param_CHUNKDIR);
  sprintf(outputfolder,param_OUTPUTDIR);
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
      snap1 = snaplist[i];
      snapid1 = i;
      sprintf(memmgr_buff,"Halo wrapper");
      if(mpi_rank==0)
	{
	  sprintf(command,"echo %d > status",i);
	  system(command);
	}
      for(l=0;l<pow3(param_chunk_per_dim);l++)
	{
	  if(mpi_rank == 0) printf("Redshift: %f\n", snap1);
	  if(l%mpi_nodes == mpi_rank)
	    {
	      printf("\treading chunk %d by rank:%d\n",l,mpi_rank);
	      halocatA = memmgr_malloc(1*sizeof(make_catalogue_halo_wrapper_t),memmgr_buff);      
	      halocatA[0] = sussexbigrun_load_halo_catalogue_binary_single_chunk(folder,snap1,snapid1,l);
	      sussexbigrun_output_cubep3m(halocatA[0],l);
	      free_make_catalogue_halo_wrapper(halocatA);
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
  
