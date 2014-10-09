#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "L-Galaxies.h"


int main()
{
  struct LGalaxy *tree;
  FILE *fp;
  int firstfile, lastfile;
  char *basename = "/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/SA_z";
  char *zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt";
  char zlist_string[100][1024];
  char filename[2048];
  char buff[1000];
  int nGals;
  int dummy,*dummyarray;
  int i,j;
  int nSnaps;
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[i]) != EOF) {
    // printf("%s%s_\n",basename,zlist_string[i]);
    i++;
  }
  nSnaps = i;
  printf("Total snapshot : %d\n",nSnaps);
  for(j=0;j<nSnaps;j++) {
    for (i=firstfile;i<=lastfile;i++) {
      sprintf(filename, "%s%s_%d",basename,zlist_string[j],i);
      printf(filename,"Reading %s\n",filename);
    }
  }
  return 0;
}
