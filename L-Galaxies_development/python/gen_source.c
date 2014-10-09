#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "L-Galaxies.h"


int main(int argc, char **argv)
{
  struct LGalaxy *lgal;
  FILE *fp;
  double *Sfr;
  int cubep3m_p = 1728;
  float boxsize = 47.0;
  int grid = 306;
  float gridsize = boxsize/grid;
  int cubep3m_cell = cubep3m_p*2;
  long long cell;
  int firstfile = 0;
  int lastfile = 127;
  char *basename = "/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/SA_z";
  char *outputfolder = "/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/sources/";
  char *zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt";
  char zlist_string[100][1024];
  char filename[2048];
  char buff[1000];
  int nGals;
  int dummy,*dummyarray;
  int i,j,k;
  int nSnaps,selected_snap;
  float gridmass;

  float Mpc2m = 3.08567758e22;
  float m2Mpc = 1./Mpc2m;
  float Msun2kg = 1.98855e30;
  float kg2Msun = 1./Msun2kg;
  float m2km = 0.001;
  float omegam = 0.27;
  float G = 6.674e-11;   // SI
  float h = 0.7;
  float H0 = 100.0;        // km/s / (Mpc/h)
  float pi = 4.0*atan(1.0);
  float rho_crit_0;
  float gridmass_c; //to convert msun to gridmass
  printf("argc:%d argv 0:%s 1:%\n",argc,argv[0],argv[1]);
  exit(1);
  if(argc == 1)
    sscanf(argv[0],"%d",&selected_snap);

  G *= (m2km*m2km) * (m2Mpc) / (kg2Msun); //  (Mpc/h) (km/s)^2 / (Msun/h)
  rho_crit_0 = 3.* (H0*H0)/ (8.*pi*G); //  # (1e10 Msun/h)/(Mpc/h)^3
  gridmass = omegam*rho_crit_0*(boxsize*boxsize*boxsize)/(cubep3m_cell*cubep3m_cell*cubep3m_cell)/h; // Msun
  gridmass_c = 1./gridmass;
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[i]) != EOF) {
    i++;
  }
  nSnaps = i;
  fclose(fp);
  printf("Total snapshot : %d\n",nSnaps);
  for(j=0;j<nSnaps;j++) {
    if(argc == 0)
      selected_snap = j;
    if(j == selected_snap) {
      Sfr = calloc(grid*grid*grid,sizeof(double));
      for (i=firstfile;i<=lastfile;i++) {
	sprintf(filename, "%s%s_%d",basename,zlist_string[j],i);
	printf("Reading %s\n",filename);
	fp = fopen(filename,"rb");
	fread(&dummy, sizeof(int), 1, fp);
	fread(&nGals, sizeof(int),1, fp);
	lgal = malloc(sizeof(struct LGalaxy)*nGals);
	fseek(fp, dummy*sizeof(int), SEEK_CUR); // skip nGalsperTree
	fread(lgal,sizeof(struct LGalaxy),nGals,fp);
	for(k=0;k<nGals;k++) {
	  cell = lgal[k].Pos[0]/gridsize + lgal[k].Pos[1]/gridsize*grid + lgal[k].Pos[2]/gridsize*grid*grid;
	  Sfr[cell] += lgal[k].Sfr*gridmass_c;
	}
	free(lgal);
	fclose(fp);
      }
      for(i=0;i<100;i++){
	printf("%d %lf\n",i,Sfr[i]); 
      }
      free(Sfr);
    }
  }
  return 0;
}
