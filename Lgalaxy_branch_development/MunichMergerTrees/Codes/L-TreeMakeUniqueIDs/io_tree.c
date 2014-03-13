#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"



void load_tree_table(int filenr)
{
  int i;
  char buf[1000];
  FILE *fd;


  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }

  fread(&Ntrees, 1, sizeof(int), fd);
  fread(&TotHalos, 1, sizeof(int), fd);

  TreeNHalos = mymalloc(sizeof(int) * Ntrees);
  TreeFirstHalo = mymalloc(sizeof(int) * Ntrees);
  Nunique = mymalloc(sizeof(int) * TotHalos);

  fread(TreeNHalos, Ntrees, sizeof(int), fd);

  fclose(fd);

  /*
  printf("Ntrees=%d TotHaloes=%d\n", Ntrees, TotHalos);
  */

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

}




void free_tree_table(void)
{
  myfree(Nunique);
  myfree(TreeFirstHalo);
  myfree(TreeNHalos);
}




void load_tree(int filenr, int nr)
{
  char buf[1000];
  FILE *fd;


  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir, LastSnapShotNr, filenr);

  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }

  fseek(fd, sizeof(int) * (2 + Ntrees), SEEK_CUR);
  fseek(fd, sizeof(struct halo_data) * TreeFirstHalo[nr], SEEK_CUR);

  Halo = mymalloc(sizeof(struct halo_data) * TreeNHalos[nr]);

  fread(Halo, TreeNHalos[nr], sizeof(struct halo_data), fd);

  fclose(fd);

  /*
  printf("tree %d loaded\n", nr);
  */

  Done = malloc(sizeof(int) * TreeNHalos[nr]);
  HaloAux = malloc(sizeof(struct haloaux_data) * TreeNHalos[nr]);
}




void free_tree(void)
{
  free(HaloAux);
  free(Done);
  myfree(Halo);
}
