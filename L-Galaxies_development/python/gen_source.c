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
  double gridmass;
  double total_cell,total_p;
  const float Mpc2m = 3.08567758e22;
  float m2Mpc = 1./Mpc2m;
  const float Msun2kg = 1.98855e30;
  float kg2Msun = 1./Msun2kg;
  const float m2km = 0.001;
  const float omegam = 0.27;
  float G = 6.674e-11;   // SI
  float h = 0.7;
  float H0 = 100.0;        // km/s / (Mpc/h)
  float pi = 4.0*atan(1.0);
  double rho_crit_0;
  double gridmass_c; //to convert msun to gridmass
  if(argc == 2)
    sscanf(argv[1],"%d",&selected_snap);
  printf("argc = %d , argv[0] = %s, argv[1] = %s\n",argc,argv[0],argv[1]);

  G = G*(m2km*m2km) * (m2Mpc) / (kg2Msun); //  (Mpc/h) (km/s)^2 / (Msun/h)
  rho_crit_0 = 3.* (H0*H0)/ (8.*pi*G); //  # (Msun/h)/(Mpc/h)^3
  printf("rho_crit_0 = %g\n",rho_crit_0);
  total_p = (double)cubep3m_p* (double)cubep3m_p* (double)cubep3m_p;
  total_cell = (double)cubep3m_cell*(double)cubep3m_cell*(double)cubep3m_cell;
  printf("total cell = %g\n",total_cell);
  gridmass = (double)omegam*(double)rho_crit_0*(double)(boxsize*boxsize*boxsize)/total_cell; // /(double)h; // Msun
  printf("G = %g, gridmass = %lf\n",G,gridmass);
  gridmass_c = 1./gridmass;
  printf("gridmass_c = %lg\n",gridmass_c);
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[i]) != EOF) {
    i++;
  }
  nSnaps = i;
  fclose(fp);
  printf("Total snapshot : %d\n",nSnaps);
  for(j=0;j<nSnaps;j++) {
    if(argc == 1)
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
	  cell = (int)(lgal[k].Pos[0]/gridsize) + (int)(lgal[k].Pos[1]/gridsize)*grid + (int)(lgal[k].Pos[2]/gridsize)*grid*grid;
	  if(cell == 0)
	    printf("Pos %f:%f:%f cell = %d Srf =%f\n", lgal[k].Pos[0], lgal[k].Pos[0], lgal[k].Pos[0],cell,lgal[k].Sfr);
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
