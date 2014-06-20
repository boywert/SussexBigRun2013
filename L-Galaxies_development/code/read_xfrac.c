#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "allvars.h"
#include "proto.h"

void get_xfrac_mesh()
{
  FILE* fp;
  char buf[1024],sbuf[1024];
  float redshift;
  int dummy;

#ifdef PARALLEL
  if(ThisTask == 0)
    {
#endif
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
#ifdef PARALLEL
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(XfracMesh, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Task:%d => %d %d %d\n",ThisTask,XfracMesh[0],XfracMesh[1],XfracMesh[2]);
#endif

}

/*  @file read_xfrac.c
    @brief Read in inonization fraction from C2Ray  
    return 0 if not found, 1 if success
*/
int read_xfrac(int snapnr, double* xfrac)
{
  FILE* fp;
  char buf[1024],sbuf[1024];
  float redshift;
  int i,j,k,il,cell;
  int dummy,mesh[3];
  int status;

  if(ThisTask == 0)
    {
      printf("Reading Xfrac data: Snap = %d\n",snapnr);
    }
  il = snapnr;

  redshift = ZZ[il];
  sprintf(buf, "%s/xfrac3d_%2.3f.bin", XfracDir, redshift);

  if((fp = fopen(buf,"r")) == NULL)
    {
      char sbuf[1000];
      printf("can't open file `%s': SKIP\n", buf);
      status = 0;
    }
  else
    {
      fread(&dummy, 1, sizeof(int),fp);
      fread(&mesh, 3, sizeof(int),fp);
      fread(&dummy, 1, sizeof(int),fp);
      if(mesh[0] != XfracMesh[0]
	 || mesh[1] != XfracMesh[1]
	 || mesh[2] != XfracMesh[2])
	{
	  sprintf(sbuf,"Meshes do not match\n");
	  terminate(sbuf);
	}
      fread(&dummy, 1, sizeof(int),fp);
      fread(xfrac, XfracMesh[0]*XfracMesh[1]*XfracMesh[2], sizeof(double),fp);
      fread(&dummy, 1, sizeof(int),fp);

      fclose(fp);
      status = 1;
    }


  return status;
}

