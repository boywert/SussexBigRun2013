#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

long long Totlen;
int Totunique;
int Count, CountUnique;

long long *TempList;


struct bigtab_data
{
  int TreeNr, HaloNr, SnapNum;
  long long ID;
} *BigTab;





int main(int argc, char **argv)
{
  int filenr, tree;

  if(argc != 2)
    {
      printf("\n  usage: L-TreeAddIDTab <parameterfile>\n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);

  printf("\n");



  for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
      load_tree_table(filenr);

      BigTab = mymalloc((LastSnapShotNr + 1) * TotHalos * sizeof(struct bigtab_data));

      Totlen = Totunique = 0;

      for(tree = 0; tree < Ntrees; tree++)
	{
	  load_tree(filenr, tree);

	  process_tree(tree);
 
	  free_tree();
	}
 
      printf("Filenr=%d \tTotlen=%lld \tTotunique=%d \tfrac=%g\n", filenr, Totlen, Totunique,
	     (float) Totunique / Totlen);

      printf("Sorting...\n");

      qsort(BigTab, Totunique, sizeof(struct bigtab_data), compare_tab);

      printf("done\n");

      save(filenr);

      myfree(BigTab);

      free_tree_table();

    }

  return 0;
}


void save(int filenr)
{
  FILE *fdout;
  char buf[1000];
  int TotSnaps;
  float pos[3];
  int i, j, *OffsetID_Halo, *OffsetID_SnapTree, *CountID_SnapTree;
  int *CountID_Snap, *OffsetID_Snap;


  sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
  if(!(fdout = fopen(buf, "w+")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }

  TotSnaps = LastSnapShotNr + 1;

  fwrite(&TotHalos, sizeof(int), 1, fdout);
  fwrite(&Totunique, sizeof(int), 1, fdout);
  fwrite(&Ntrees, sizeof(int), 1, fdout);
  fwrite(&TotSnaps, sizeof(int), 1, fdout);

  CountID_Snap = mymalloc(TotSnaps * sizeof(int));
  OffsetID_Snap = mymalloc(TotSnaps * sizeof(int));

  CountID_SnapTree = mymalloc(TotSnaps * Ntrees * sizeof(int));
  OffsetID_SnapTree = mymalloc(TotSnaps * Ntrees * sizeof(int));
  OffsetID_Halo = mymalloc(TotHalos * sizeof(int));

  for(i = 0; i < TotSnaps; i++)
    CountID_Snap[i] = OffsetID_Snap[i] = 0;

  for(i = 0; i < TotSnaps; i++)
    for(j = 0; j < Ntrees; j++)
      CountID_SnapTree[i * Ntrees + j] = OffsetID_SnapTree[i * Ntrees + j] = 0;

  for(i = 0; i < TotHalos; i++)
    OffsetID_Halo[i] = -1;



  for(i = 0; i < Totunique; i++)
    {
      CountID_Snap[BigTab[i].SnapNum]++;

      CountID_SnapTree[BigTab[i].SnapNum * Ntrees + BigTab[i].TreeNr]++;

      if(OffsetID_Halo[BigTab[i].HaloNr] == -1)
	OffsetID_Halo[BigTab[i].HaloNr] = i;
    }


  for(i = 1; i < TotSnaps * Ntrees; i++)
    OffsetID_SnapTree[i] = OffsetID_SnapTree[i - 1] + CountID_SnapTree[i - 1];

  for(i = 1; i < TotSnaps; i++)
    OffsetID_Snap[i] = OffsetID_Snap[i - 1] + CountID_Snap[i - 1];

  fwrite(CountID_Snap, sizeof(int), TotSnaps, fdout);
  fwrite(OffsetID_Snap, sizeof(int), TotSnaps, fdout);

  fwrite(CountID_SnapTree, sizeof(int), TotSnaps * Ntrees, fdout);
  fwrite(OffsetID_SnapTree, sizeof(int), TotSnaps * Ntrees, fdout);

  fwrite(Nunique, sizeof(int), TotHalos, fdout);
  fwrite(OffsetID_Halo, sizeof(int), TotHalos, fdout);




  myfree(OffsetID_Halo);
  myfree(OffsetID_SnapTree);
  myfree(CountID_SnapTree);
  myfree(OffsetID_Snap);
  myfree(CountID_Snap);




  for(i = 0; i < Totunique; i++)
    fwrite(&BigTab[i].ID, sizeof(long long), 1, fdout);


  /* append zero coordinates */

  pos[0] = pos[1] = pos[2] = 0;

  for(i = 0; i < Totunique; i++)
    fwrite(pos, sizeof(float), 3, fdout);

  /* append zero velocities */

  for(i = 0; i < Totunique; i++)
    fwrite(pos, sizeof(float), 3, fdout);

  fclose(fdout);
}


void process_tree(int tree)
{
  int i, halonr;


  TempList = mymalloc(sizeof(long long) * (LastSnapShotNr + 1) * TreeNHalos[tree]);


  for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
    {
      Count = CountUnique = 0;

      find_progenitors(tree, halonr, Halo[halonr].FirstProgenitor);
      
      qsort(TempList, Count, sizeof(long long), compare_id);

      for(i = 0; i < Count; i++)
	{
	  if(CountUnique > 0)
	    {
	      if(TempList[i] == BigTab[Totunique + CountUnique - 1].ID)
		continue;
	    }

	  BigTab[Totunique + CountUnique].ID = TempList[i];
	  BigTab[Totunique + CountUnique].TreeNr = tree;
	  BigTab[Totunique + CountUnique].HaloNr = TreeFirstHalo[tree] + halonr;
	  BigTab[Totunique + CountUnique].SnapNum = Halo[halonr].SnapNum;

	  CountUnique++;
	}

      Nunique[TreeFirstHalo[tree] + halonr] = CountUnique;

      Totlen += Halo[halonr].Len;
      Totunique += CountUnique;
    }

  myfree(TempList);
}


void find_progenitors(int tree, int halonr, int prog)
{
  if(prog >= 0)
    {
      TempList[Count++] = Halo[prog].MostBoundID;
     
      if(Halo[prog].FirstProgenitor >= 0)
	find_progenitors(tree, halonr, Halo[prog].FirstProgenitor);

      if(Halo[prog].NextProgenitor >= 0)
	find_progenitors(tree, halonr, Halo[prog].NextProgenitor);
    }
}


int compare_id(const void *a, const void *b)
{
  if(*(long long *) a < *(long long *) b)
    return -1;

  if(*(long long *) a > *(long long *) b)
    return +1;

  return 0;
}


int compare_tab(const void *a, const void *b)
{
  if(((struct bigtab_data *) a)->SnapNum < ((struct bigtab_data *) b)->SnapNum)
    return -1;

  if(((struct bigtab_data *) a)->SnapNum > ((struct bigtab_data *) b)->SnapNum)
    return +1;

  if(((struct bigtab_data *) a)->TreeNr < ((struct bigtab_data *) b)->TreeNr)
    return -1;

  if(((struct bigtab_data *) a)->TreeNr > ((struct bigtab_data *) b)->TreeNr)
    return +1;

  if(((struct bigtab_data *) a)->HaloNr < ((struct bigtab_data *) b)->HaloNr)
    return -1;

  if(((struct bigtab_data *) a)->HaloNr > ((struct bigtab_data *) b)->HaloNr)
    return +1;

  if(((struct bigtab_data *) a)->ID < ((struct bigtab_data *) b)->ID)
    return -1;

  if(((struct bigtab_data *) a)->ID > ((struct bigtab_data *) b)->ID)
    return +1;

  return 0;
}
