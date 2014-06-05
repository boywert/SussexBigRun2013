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

void load_xfrac()
{
  FILE* fp;
  char buf[1024],buf2[1024];
  float redshift;
  int i,j,k,il,cell;
  int dummy;
 
  if(ThisTask == 0)
    {
      printf("Reading Xfrac data\n\n");
    }
  for(il=0;il<MAXSNAPS;il++)
    {
      redshift = ZZ[il];
      sprintf(buf, "%s/xfrac3d_%2.3f.bin", XfracDir, redshift);
      if(!(fp = fopen(buf,"r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s': SKIP\n", buf);
	  XfracDataDone[il] = 0;
	}
      else
	{
	  fread(&dummy, 1, sizeof(int),fp);
	  myfread(XfracMesh, 3, sizeof(int),fp);
	  fread(&dummy, 1, sizeof(int),fp);
	  sprintf(buf2,"XfracData[%d]",il);
	  XfracData[il] = mymalloc(buf2,sizeof(float)*XfracMesh[0]*XfracMesh[1]*XfracMesh[2]);
	  fread(&dummy, 1, sizeof(int),fp);
	  myfread(XfracData[il], XfracMesh[0]*XfracMesh[1]*XfracMesh[2], sizeof(double),fp);
	  fread(&dummy, 1, sizeof(int),fp);

	  /* for(i=0;i<XfracMesh[0];i++) */
	  /*   { */
	  /*     for(j=0;j<XfracMesh[1];j++) */
	  /* 	{ */
	  /* 	  for(k=0;k<XfracMesh[2];k++) */
	  /* 	    { */
	  /* 	      cell = k*XfracMesh[0]*XfracMesh[1]+j*XfracMesh[0]+i; */
	  /* 	      printf("%d: %lf\n",cell,XfracData[il][cell]); */
	  /* 	    } */
	  /* 	} */
	  /*   } */
	  fclose(fp);
	  XfracDataDone[il] = 1;
	}
    }
}

void free_xfrac()
{
  int i;
  for(i=0;i<MAXSNAPS;i++)
    {
      myfree(XfracData[i]);
    }
}
