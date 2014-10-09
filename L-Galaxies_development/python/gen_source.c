#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "L-Galaxies.h"


int main()
{
  struct LGalaxy *lgal;
  FILE *fp;
  double *Sfr;
  float boxsize = 47.0;
  int grid = 306;
  float gridsize = boxsize/grid;
  long long cell;
  int firstfile = 0;
  int lastfile = 127;
  char *basename = "/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/SA_z";
  char *zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt";
  char zlist_string[100][1024];
  char filename[2048];
  char buff[1000];
  int nGals;
  int dummy,*dummyarray;
  int i,j,k;
  int nSnaps;
  
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[i]) != EOF) {
    i++;
  }
  nSnaps = i;
  fclose(fp);
  printf("Total snapshot : %d\n",nSnaps);
  for(j=0;j<nSnaps;j++) {			
    Sfr = calloc(grid*grid*grid,sizeof(double));
    for (i=firstfile;i<=lastfile;i++) {
      sprintf(filename, "%s%s_%d",basename,zlist_string[j],i);
      printf("Reading %s\n",filename);
      fp = fopen(filename,"rb");
      fread(dummy, sizeof(int), 1, fp);
      fread(nGals, sizeof(int),1, fp);
      lgal = malloc(sizeof(struct LGalaxy)*nGals);
      fseek(fp, dummy*sizeof(int), SEEK_CUR); // skip nGalsperTree
      fread(lgal,sizeof(struct LGalaxy),nGals,fp);
      for(k=0;k<nGals;k++) {
	cell = lgal[k].Pos[0]/gridsize + lgal[k].Pos[1]/gridsize*grid + lgal[k].Pos[2]/gridsize*grid*grid;
	Sfr[cell] += lgal[k].Sfr;
      }
      free(lgal)
    }
    free(Sfr);
  }
  return 0;
}
