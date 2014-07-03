/** @file main.c
 * @brief The file containing the main() function for L-Galaxies.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"


/**@file main.c
 * @brief Controlling function of L-Galaxies plus Construct Galaxies,
 *        Join Galaxies of progenitors, Evolve Galaxies and check
 *        Makefile options.
 * */
int main(int argc, char **argv) {
  int filenr, tree, halonr,filenrid,i,j;
  struct stat filestatus;
  FILE *fd, *fa;
  char buf[1000];
  time_t initial, final, previous;
  //int NFOFs=767;
  int *FileList, *TreeList;
  int MyNtrees, MytotNHalos, *MyTreeNHalos, NFOFs;
  float weight;
  long long haloid;
  int FileNumber=409024;

  /*Reads the parameter file, given as an argument at run time.*/
  read_parameter_file(argv[1]);

  sprintf(buf, "./sample_nh_%d.dat",FileNumber);
  if(!(fd = fopen(buf, "r")))
     {
       printf("can't open file place 1 `%s'\n", buf);
       exit(1);
     }

  fscanf(fd, "%d \n", &MyNtrees);
  printf("NTrees=%d\n", MyNtrees);


  FileList = malloc(sizeof(int) * MyNtrees);
  TreeList = malloc(sizeof(int) * MyNtrees);

  for(i=0;i<MyNtrees;i++)
  {
	  fscanf(fd, "%lld %d %d %f\n", &haloid, &TreeList[i], &FileList[i], &weight);
	  //printf("i=%d list=%d \n",i,FileList[i]);
  }
  fclose(fd);

  MyTreeNHalos = malloc(sizeof(int) * MyNtrees);
  MytotNHalos=0;
  //MyNtrees=0;
 printf("\n\nReading Started... -> output file %d\n\n\n",FileNumber);

 //tree number 3 numbers for the redshift + Nsampling
 sprintf(buf, "/afs/mpa/data/bmh20/workspace/performance_MCMC_Scaling_backup/MergerTrees/MergerTrees_1/treedata/trees_063.%d",FileNumber);
	  if(!(fd = fopen(buf, "wb")))
	     {
	       printf("can't open file %d\n", __LINE__);
	       exit(1);
	     }

 sprintf(buf, "/afs/mpa/data/bmh20/workspace/performance_MCMC_Scaling_backup/MergerTrees/MergerTrees_1/treedata/tree_dbids_063.%d", FileNumber);
 if(!(fa = fopen(buf, "wb")))
 {
	 printf("can't open file %d\n", __LINE__);
	 exit(1);
 }


//jump header in output file
 /*for (filenr = 0; filenr < NFOFs; filenr++)
   		  	  if(filenr < 0 || ((FileList[filenr]+TreeList[filenr])!= (FileList[filenr-1]+TreeList[filenr-1])))
   		  	      ++MyNtrees;*/


	  fseek(fd, sizeof(int) * (2 + MyNtrees),0);

  	  for (filenr = 0; filenr < MyNtrees; filenr++)
  	  	  {
  		  	   		  load_tree_table(FileList[filenr]);
  		  		  	  load_tree(FileList[filenr], TreeList[filenr]);

  		  		  	  MytotNHalos+=TreeNHalos[TreeList[filenr]];
  		  		  	  MyTreeNHalos[filenr]=TreeNHalos[TreeList[filenr]];

  		  		     printf("Treenr=%d (of %d), Tree Number of Halos=%d\n",filenr, MyNtrees, TreeNHalos[TreeList[filenr]]);

  		  		      if(fwrite(HaloIDs, sizeof(struct halo_ids_data), MyTreeNHalos[filenr], fa) !=  MyTreeNHalos[filenr])
  		  		    		  		  		  printf("error in fwrite\n");

  		  		  	  if(fwrite(Halo, sizeof(struct halo_data), MyTreeNHalos[filenr], fd) !=  MyTreeNHalos[filenr])
  		  		  		  printf("error in fwrite\n");


  		  		  	  free(HaloAux);
  		  		  	  free(Halo);
  		  		  	  free(HaloIDs);

  		  		  	  for(i = NOUT - 1; i >= 0; i--)
  		  		  		  free(TreeNgals[i]);
  		  		  	  free(TreeFirstHalo);
  		  		  	  free(TreeNHalos);

  		  		  	  //printf("done file %d\n", filenr);
  		  	  	 // }
  	  	  }
  	  fclose(fa);
//write header
  	  printf("\n\nRewinding...please wait...\n\n");
      rewind(fd);
  	  fwrite(&MyNtrees, sizeof(int), 1, fd);
      fwrite(&MytotNHalos, sizeof(int), 1, fd);
      fwrite(MyTreeNHalos, sizeof(int), MyNtrees, fd);
      free(MyTreeNHalos);
      free(FileList);
      free(TreeList);
      fclose(fd);


      printf("\n\ndone");

}





