#include "main.h"
#include "treeformat.h"

int main()
{
  struct Lgalaxy_halo_data sample;
  int i,j,firstsnap,nhalopertree,ntrees,nhalos;
  FILE *fp;
  firstsnap = 0;
  nhalopertree = 61-firstsnap+1;
  ntrees = 1;
  nhalos = ntrees*nhalopertree;
  sample.Len = 1000;
  sample.M_Mean200 = 0.0;
  sample.M_Crit200 = 100.;
  sample.M_TopHat = 0.0;
  sample.Pos[0] = 1.0;
  sample.Pos[1] = 1.0;
  sample.Pos[2] = 1.0;
  sample.Vel[0] = 0.0;
  sample.Vel[1] = 0.0;
  sample.Vel[2] = 0.0;
  sample.VelDisp = 100.0;
  sample.Vmax = 100.0;
  sample.Spin[0] = 10.0;
  sample.Spin[1] = 0.0;
  sample.Spin[1] = 0.0;
  sample.MostBoundID = 0;
  fp = fopen("test","wb+");
  fwrite (&(ntrees),1, sizeof(int), fp);
  fwrite (&(nhalos),1, sizeof(int), fp);
  /* for(j=0;j<ntrees;j++) */
  /*   { */
  /*     for(i=61;i>=firstsnap;i--) */
  /* 	{ */
  /* 	  sample.Descendant = i-61-1; */
  /* 	  sample.FirstProgenitor = i-61+1; */
  /* 	  sample.NextProgenitor = -1; */
  /* 	  sample.FirstHaloInFOFgroup = i-61; */
  /* 	  sample.NextHaloInFOFgroup = -1; */
  /* 	  sample.SnapNum = i; */
  /* 	  sample.FileNr = 0; */
  /* 	  sample.SubhaloIndex = 0; */
  /* 	  sample.SubHalfMass = 0.1; */
  /* 	  fwrite(&(sample),sizeof(struct Lgalaxy_halo_data),1, fp); */
  /* 	} */
  /*   } */
  fclose(fp);
  return 0;
}
