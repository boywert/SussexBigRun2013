#include "main.h"
#include "treeformat.h"

int main()
{
  struct Lgalaxy_halo_data sample;
  int i,j,firstsnap,nhalopertree,ntrees,nhalos,number;
  FILE *fp;
  firstsnap = 40;
  nhalopertree = 61-firstsnap+1;
  ntrees = 2;
  nhalos = ntrees*nhalopertree;
  sample.Len = 1000;
  sample.M_Mean200 = 100.0;
  sample.M_Crit200 = 100.;
  sample.M_TopHat = 100.0;
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
  number = 1;
  fwrite (&number,1, sizeof(int), fp);
  printf("Ntrees = %d\n",ntrees);
  number = 44;
  fwrite (&number,1, sizeof(int), fp);
  printf("Nhalos = %d\n",nhalos);
  for(j=0;j<ntrees;j++)
    {
      //fwrite (&nhalopertree,1, sizeof(int), fp);
      printf("Tree:%d Nhalos = %d\n",j,nhalopertree);
    }
  number = nhalopertree*2;
  fwrite (&number,1, sizeof(int), fp);
  for(j=0;j<ntrees+2;j++)
    {
      //fwrite (&j, 1, sizeof(int), fp);
    }
  for(j=0;j<ntrees;j++)
    {
      sample.Pos[0] = 0.5*j;
      for(i=61;i>=firstsnap;i--)
  	{
  	  sample.Descendant = 61-i-1+j*nhalopertree;
	  if(i==61) sample.Descendant = -1;
  	  sample.FirstProgenitor = 61-i+1+j*nhalopertree;
  	  sample.NextProgenitor = -1;
  	  sample.FirstHaloInFOFgroup = 61-i;
	  if(j==1 && i==61) sample.FirstHaloInFOFgroup = 0;
	  if(i==firstsnap) sample.FirstProgenitor = -1;
  	  sample.NextHaloInFOFgroup = -1;
	  if(j==0 && i==61) sample.NextHaloInFOFgroup = nhalopertree;
  	  sample.SnapNum = i;
  	  sample.FileNr = 0;
  	  sample.SubhaloIndex = sample.FirstHaloInFOFgroup;
  	  sample.SubHalfMass = 0.1;
	  /* fwrite(&(sample.Descendant),sizeof(int),1, fp); */
	  /* fwrite(&(sample.FirstProgenitor),sizeof(int),1, fp); */
	  /* fwrite(&(sample.NextProgenitor),sizeof(int),1, fp); */
	  /* fwrite(&(sample.FirstHaloInFOFgroup),sizeof(int),1, fp); */
	  /* fwrite(&(sample.NextHaloInFOFgroup),sizeof(int),1, fp); */
	  /* fwrite(&(sample.Len),sizeof(int),1, fp); */
	  /* fwrite(&(sample.M_Mean200),sizeof(float),1, fp); */
	  /* fwrite(&(sample.M_Crit200),sizeof(float),1, fp); */
	  /* fwrite(&(sample.M_TopHat),sizeof(float),1, fp); */
	  /* fwrite(&(sample.Pos[0]),sizeof(float),3, fp); */
	  /* fwrite(&(sample.Vel[0]),sizeof(float),3, fp); */
	  /* fwrite(&(sample.VelDisp),sizeof(float),1, fp); */
	  /* fwrite(&(sample.Vmax),sizeof(float),1, fp); */
	  /* fwrite(&(sample.Spin[0]),sizeof(float),3, fp); */
	  /* fwrite(&(sample.MostBoundID),sizeof(long long),1, fp); */
	  /* fwrite(&(sample.SnapNum),sizeof(int),1, fp); */
	  /* fwrite(&(sample.FileNr),sizeof(int),1, fp); */
	  /* fwrite(&(sample.SubhaloIndex),sizeof(int),1, fp); */
	  /* fwrite(&(sample.SubHalfMass),sizeof(float),1, fp); */
  

  	  fwrite(&(sample),sizeof(struct Lgalaxy_halo_data),1, fp);
	  //printf("%d, %d, %d\n",i,sample.FirstProgenitor,sample.Descendant);
  	}
    }
  fclose(fp);
  return 0;
}
