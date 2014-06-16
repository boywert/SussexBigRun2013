#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


/*  @file read_xfrac.c
    @brief Read in inonization fraction from C2Ray  
*/

void load_xfrac(int snapnr)
{
  FILE* fp;
  char buf[1024],buf2[1024];
  float redshift;
  int i,j,k,il,cell;
  int dummy,mesh[3];
  if(ThisTask == 0)
    {
      printf("Reading Xfrac data\n\n");
    }

  /* Read mesh from LastSnap */
  redshift = ZZ[MAXSNAPS-1];
  sprintf(buf, "%s/xfrac3d_%2.3f.bin", XfracDir, redshift);
  if((fp = fopen(buf,"r")) == NULL)
    {
      char sbuf[1000];
      printf("can't open file `%s'\n", buf);
      terminate(sbuf);
    }
  fread(&dummy, 1, sizeof(int),fp);
  myfread(XfracMesh, 3, sizeof(int),fp);
  fread(&dummy, 1, sizeof(int),fp);
  fclose(fp);

  il = snapnr;
  redshift = ZZ[il];
  sprintf(buf, "%s/xfrac3d_%2.3f.bin", XfracDir, redshift);
  if((fp = fopen(buf,"r")) == NULL)
    {
      char sbuf[1000];
      printf("can't open file `%s': SKIP\n", buf);
      XfracDataDone[il] = 1;
      XfracData[il] = calloc(XfracMesh[0]*XfracMesh[1]*XfracMesh[2],sizeof(double));
    }
  else
    {
      fread(&dummy, 1, sizeof(int),fp);
      fread(&mesh, 3, sizeof(int),fp);
      fread(&dummy, 1, sizeof(int),fp);
      printf("z=%2.3f, Mesh: %d %d %d\n",redshift,XfracMesh[0],XfracMesh[1],XfracMesh[2]);
      sprintf(buf2,"XfracData[il]",il);
      XfracData[il] = malloc(sizeof(double)*XfracMesh[0]*XfracMesh[1]*XfracMesh[2]);
      fread(&dummy, 1, sizeof(int),fp);
      fread(XfracData[il], XfracMesh[0]*XfracMesh[1]*XfracMesh[2], sizeof(double),fp);
      fread(&dummy, 1, sizeof(int),fp);

      fclose(fp);
      XfracDataDone[il] = 1;
    }
    
}

void free_xfrac(int snapnr)
{

  free(XfracData[snapnr]);
 
}
