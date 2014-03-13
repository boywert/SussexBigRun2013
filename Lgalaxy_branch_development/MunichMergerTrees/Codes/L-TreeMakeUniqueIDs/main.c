#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"


long long HaloCount;
FILE *FdOut;

#define DOUBLE_to_HASHBITS(y) ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - Hashbits)))

int main(int argc, char **argv)
{
  int filenr, tree;

  if(argc != 2)
    {
      printf("\n  usage: L-TreeMakeUniqueIDs <parameterfile>\n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);

  read_snap_list();

  for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
      load_tree_table(filenr);

      prepare_tree_dbids_output(filenr);

      for(tree = 0; tree < Ntrees; tree++)
	{
          HaloCount = (((filenr * (long long)1000000) + tree) * (long long)1000000);

	  load_tree(filenr, tree);

          process_tree(filenr, tree);
          
	  save_tree_dbids(tree);

	  free_tree();
	}

      close_tree_dbids_output();

      free_tree_table();

      printf("Filenr=%d  HaloCount=%9d%09d\n", filenr, (int) (HaloCount / 1000000000), (int) (HaloCount % 1000000000));
    }
  return 0;
}


int walk(int nr)
{
  int last;

  last = nr;

  if(Done[nr] == 0)
    {
      Done[nr] = 1;
      HaloAux[nr].HaloID = HaloCount++;

      HaloAux[nr].FirstProgenitor = Halo[nr].FirstProgenitor;
      HaloAux[nr].Descendant = Halo[nr].Descendant;
      HaloAux[nr].NextProgenitor = Halo[nr].NextProgenitor;
      HaloAux[nr].FirstHaloInFOFgroup = Halo[nr].FirstHaloInFOFgroup;
      HaloAux[nr].NextHaloInFOFgroup = Halo[nr].NextHaloInFOFgroup;
      
      HaloAux[nr].Redshift = ZZ[Halo[nr].SnapNum];

      if(Halo[nr].FirstProgenitor >= 0)
        {
          last = walk(Halo[nr].FirstProgenitor);
        }

      HaloAux[nr].LastProgenitor = last;

      if(Halo[nr].NextProgenitor >= 0)
        {
          last = walk(Halo[nr].NextProgenitor);
        }
    }

  return last;
}


void process_tree(int filenr, int tree)
{
  int halonr, i, j, k, num;
  double scalefac, xx, yy, zz;

  scalefac = 1.0 / BoxSize;
  
  
  for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
    {
      Done[halonr] = 0;
      HaloAux[halonr].FirstProgenitor = -1;
      HaloAux[halonr].LastProgenitor = -1;
      HaloAux[halonr].NextProgenitor = -1;
      HaloAux[halonr].Descendant = -1;
      HaloAux[halonr].FirstHaloInFOFgroup = -1;
      HaloAux[halonr].NextHaloInFOFgroup = -1;
      HaloAux[halonr].FileTreeNr = (((filenr * (long long) 1000000) + tree) * (long long) 1000000);
  
      xx = Halo[halonr].Pos[0] * scalefac + 1.0;
      yy = Halo[halonr].Pos[1] * scalefac + 1.0;
      zz = Halo[halonr].Pos[2] * scalefac + 1.0;
      
      i = DOUBLE_to_HASHBITS(xx);
      j = DOUBLE_to_HASHBITS(yy);
      k = DOUBLE_to_HASHBITS(zz);

      HaloAux[halonr].PeanoKey = peano_hilbert_key(i, j, k, Hashbits);
  }

  for(num = LastSnapShotNr; num >= 0; num--)
    for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
      {
        if(Halo[halonr].SnapNum == num)
          if(Done[halonr] == 0)
            {
              walk(halonr);
            }
      }

  for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
    {
      if(HaloAux[halonr].FirstProgenitor >= 0)
	HaloAux[halonr].FirstProgenitor = HaloAux[HaloAux[halonr].FirstProgenitor].HaloID;

      if(HaloAux[halonr].LastProgenitor >= 0)
	HaloAux[halonr].LastProgenitor = HaloAux[HaloAux[halonr].LastProgenitor].HaloID;

      if(HaloAux[halonr].NextProgenitor >= 0)
	HaloAux[halonr].NextProgenitor = HaloAux[HaloAux[halonr].NextProgenitor].HaloID;
 
      if(HaloAux[halonr].Descendant >= 0)
	HaloAux[halonr].Descendant = HaloAux[HaloAux[halonr].Descendant].HaloID;

      if(HaloAux[halonr].FirstHaloInFOFgroup >= 0)
	HaloAux[halonr].FirstHaloInFOFgroup = HaloAux[HaloAux[halonr].FirstHaloInFOFgroup].HaloID;

      if(HaloAux[halonr].NextHaloInFOFgroup >= 0)
	HaloAux[halonr].NextHaloInFOFgroup = HaloAux[HaloAux[halonr].NextHaloInFOFgroup].HaloID;


      if(HaloAux[halonr].FirstProgenitor < 0)
        HaloAux[halonr].LastProgenitor = HaloAux[halonr].HaloID;

    }

}



void save_tree_dbids(int tree)
{
  fwrite(HaloAux, sizeof(struct haloaux_data), TreeNHalos[tree], FdOut);
}



void prepare_tree_dbids_output(int filenr)
{
  char buf[1000];

  sprintf(buf, "%s/treedata/tree_dbids_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
  if(!(FdOut = fopen(buf, "w")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }
}

void close_tree_dbids_output(void)
{
  fclose(FdOut);
}


void read_snap_list(void)
{
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithSnapList);

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      exit(1);
    }

  Snaplistlen = 0;
  do
    {
      if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
        {
          ZZ[Snaplistlen] = 1 / AA[Snaplistlen] - 1;
          Snaplistlen++;
        }
      else
	break;
    }
  while(Snaplistlen < MAXSNAPS);

  fclose(fd);

  printf("found %d defined times in snaplist.\n", Snaplistlen);
}

