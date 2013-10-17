#include "main.h"
#include "treeformat.h"

int main()
{
  struct Lgalaxy_halo_data sample;
  int i,j,firstsnap,nhalopertree,ntrees,nhalos;
  FILE *fp;
  firstsnap = 40;
  nhalopertree = 61-firstsnap+1;
  ntrees = 10;
  nhalos = ntrees*nhalopertree;
  sample.Len = 1000;
  sample.M_Mean200 = 0.0;
  sample.M_Crit200 = 100.;
  sample.M_TopHat = 0.0;
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
  printf("struct's size = %d\n",sizeof(struct Lgalaxy_halo_data));
  fp = fopen("treedata/trees_061.0","wb+");
  fwrite (&ntrees,1, sizeof(int), fp);
  printf("Ntrees = %d\n",ntrees);
  fwrite (&nhalos,1, sizeof(int), fp);
  printf("Nhalos = %d\n",nhalos);
  for(j=0;j<ntrees;j++)
    {
      fwrite (&nhalopertree,1, sizeof(int), fp);
      printf("Tree:%d Nhalos = %d\n",j,nhalopertree);
    }
  for(j=0;j<ntrees+2;j++)
    {
      fwrite (&nhalopertree,1, sizeof(int), fp);
    }
  for(j=0;j<ntrees;j++)
    {
      sample.Pos[0] = 1.0*j;
      for(i=61;i>=firstsnap;i--)
  	{
  	  sample.Descendant = 61-i-1+j*nhalopertree;
  	  sample.FirstProgenitor = 61-i+1+j*nhalopertree;
  	  sample.NextProgenitor = -1;
  	  sample.FirstHaloInFOFgroup = 61-i+j*nhalopertree;
	  if(i==firstsnap) sample.FirstProgenitor = -1;
  	  sample.NextHaloInFOFgroup = -1;
  	  sample.SnapNum = i;
  	  sample.FileNr = 0;
  	  sample.SubhaloIndex = 0;
  	  sample.SubHalfMass = 0.1;
  	  fwrite(&(sample),sizeof(struct Lgalaxy_halo_data),1, fp);
	  printf("%d, %d, %d\n",i,sample.FirstProgenitor,sample.Descendant);
  	}
    }
  fclose(fp);
  return 0;
}
