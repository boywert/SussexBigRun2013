#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

void load_tree_table(int filenr, int simulation_number)
{
  int i, n, totNHalos, Ntrees;
  char buf[1000];
  FILE *fd;

  if(simulation_number==1)
  	sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir_MR, LastSnapShotNr, filenr);

  if(simulation_number==2)
  	sprintf(buf, "%s/treedata/trees_sf1_%03d.%d", SimulationDir_MRII, LastSnapShotNr, filenr);

  if(!(fd = fopen(buf, "rb")))
    {
      printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
      exit(1);
    }

  //read header on trees_** file
  fread(&Ntrees, sizeof(int), 1, fd);
  fread(&totNHalos, sizeof(int), 1, fd);

  TreeNHalos = malloc(sizeof(int) * Ntrees);
  TreeFirstHalo = malloc(sizeof(int) * Ntrees);

  fread(TreeNHalos, sizeof(int), Ntrees, fd);

  fclose(fd);

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  /*Define a variable containing the number you have to jump to
   * get from one firshalo to the next. */
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];



}


void load_tree(int filenr, int treenhalos, int treefirsthalo, int simulation_number)
{
	int Ntrees;
	FILE *fd;
  char buf[1000];

  //READ Tree
  if(simulation_number==1)
  	sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir_MR, LastSnapShotNr, filenr);

  if(simulation_number==2)
  	sprintf(buf, "%s/treedata/trees_sf1_%03d.%d", SimulationDir_MRII, LastSnapShotNr, filenr);


  if(!(fd = fopen(buf, "rb")))
    {
      printf("can't open file '%s' in io_tree.c line %d\n",buf, __LINE__);
      exit(1);
    }

  fread(&Ntrees, sizeof(int), 1, fd);
  //printf("load trees Ntrees=%d\n",Ntrees);
  fseek(fd, sizeof(int) * (1 + Ntrees), SEEK_CUR);
  fseek(fd, sizeof(struct halo_data) * treefirsthalo, SEEK_CUR);

  Halo = malloc(sizeof(struct halo_data) * treenhalos);
  fread(Halo, sizeof(struct halo_data), treenhalos, fd);



  //printf("FIRSTPROG=%d firsthalo=%d\n",Halo[0].FirstProgenitor,Halo[0].FirstHaloInFOFgroup);

  fclose(fd);

}


void load_tree_IDs(int filenr, int treenhalos, int treefirsthalo, int simulation_number)
{
	FILE *fd;
  char buf[1000];

  //READ IDS
  if(simulation_number==1)
  	sprintf(buf, "%s/treedata/tree_dbids_%03d.%d", SimulationDir_MR, LastSnapShotNr, filenr);

  if(simulation_number==2)
  	sprintf(buf, "%s/treedata/tree_sf1_dbids_%03d.%d", SimulationDir_MRII, LastSnapShotNr, filenr);

  if(!(fd = fopen(buf, "rb")))
     {
       printf("can't open file '%s' in io_tree.c line %d\n",buf, __LINE__);
       exit(1);
     }

  if(simulation_number==1)
  {
  HaloIDs = malloc(sizeof(struct halo_ids_data) * treenhalos);
  fseek(fd, sizeof(struct halo_ids_data) * treefirsthalo, SEEK_CUR);
  fread(HaloIDs, sizeof(struct halo_ids_data), treenhalos, fd);
  }

  if(simulation_number==2)
   {
   HaloIDs_MRII = malloc(sizeof(struct halo_ids_data_MRII) * treenhalos);
   fseek(fd, sizeof(struct halo_ids_data_MRII) * treefirsthalo, SEEK_CUR);
   fread(HaloIDs_MRII, sizeof(struct halo_ids_data_MRII), treenhalos, fd);
   }


  fclose(fd);

}

void load_aux_header(int filenr, int simulation_number)
{
	int aux_trees;
	char buf[1000];
  FILE *fd;
  struct stat filestatus;

  if(simulation_number==1)
  	sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir_MR, LastSnapShotNr, filenr);

  if(simulation_number==2)
  	sprintf(buf, "%s/treedata/treeaux_sf1_%03d.%d", SimulationDir_MRII, LastSnapShotNr, filenr);


  if(stat(buf, &filestatus) != 0)                  /* seems not to exist */
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
    }

  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
    }

  fread(&NtotHalos, 1, sizeof(int), fd);
  fread(&TotIds, 1, sizeof(int), fd);
  fread(&aux_trees, 1, sizeof(int), fd);
  fread(&TotSnaps, 1, sizeof(int), fd);

  fseek(fd, -4*sizeof(int), SEEK_CUR);
  size_t bytes = 4 * sizeof(int)+2 * TotSnaps * sizeof(int) + 2 * TotSnaps * aux_trees * sizeof(int) + 2 * NtotHalos * sizeof(int);
  TreeAuxData = malloc(bytes);
  fread(TreeAuxData, 1, bytes, fd);

  fclose(fd);


}




void load_all_auxdata(int filenr, int simulation_number)
{
	char buf[1000];
  FILE *fd;
  struct stat filestatus;
  int aNtotHalos, aTotIds, aAUXNtrees, aTotSnaps;
  int acount_snap[64], aoff_snap[64], *acount_tree, *aoff_tree, *acount_halo, *aoff_halo;

  if(simulation_number==1)
  	sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir_MR, LastSnapShotNr, filenr);

  if(simulation_number==2)
  	sprintf(buf, "%s/treedata/treeaux_sf1_%03d.%d", SimulationDir_MRII, LastSnapShotNr, filenr);


  if(stat(buf, &filestatus) != 0)                  /* seems not to exist */
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
    }

  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file `%s'\n", buf);
    }

  size_t bytes = filestatus.st_size;

  TreeAuxData = malloc(bytes);
  fread(TreeAuxData, 1, bytes, fd);

  fclose(fd);


 /* int *header = TreeAuxData;

   NtotHalos = header[0];
   TotIds = header[1];
   Ntrees = header[2];
   TotSnaps = header[3];

   CountIDs_snaptree = header + 4 + 2 * TotSnaps;
   OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * Ntrees;
   CountIDs_halo = OffsetIDs_snaptree + TotSnaps * Ntrees;
   OffsetIDs_halo = CountIDs_halo + NtotHalos;*/
}

