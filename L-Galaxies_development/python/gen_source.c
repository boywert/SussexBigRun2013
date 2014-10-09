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
  int i;
  int nSnaps;
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[ii]) != EOF) {
    printf("%s%s_\n",basename,zlist_string[ii]);
    i++;
  }
  nSnaps = i;
  

  for (i=firstfile;i<=lastfile;i++) {
    sprintf(filename,"s");
  }
  
  return 0;
}
