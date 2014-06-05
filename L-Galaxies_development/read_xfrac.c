#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


/*  @file read_xfrac.c
    @brief Read in inonization fraction from C2Ray  
*/

void load_xfrac(int SnapNum)
{
  FILE* fp;
  char buf[1024];
  int ncells=306;
  float redshift;
  float *XfracData;
  char *XfracDir;
  XfracDir = "/research/prace/47Mpc_RT/47Mpc_f2_0_306/results/";
  redshift = 6.0;
  sprintf(buf, "%s/xfrac3d_%2.3f.bin", XfracDir,redshift);
  if(!(fp = fopen(buf,"r")))
    {
      /* char sbuf[1000]; */
      /* sprintf(sbuf, "can't open file `%s'\n", buf); */
      /* terminate(sbuf); */
      exit(1);
    }
  // XfracData = mymalloc("XfracData",sizeof(float)*ncells);
  XfracData = malloc(sizeof(float)*ncells*ncells*ncells);
  fread(XfracData, ncells*ncells*ncells, sizeof(float),fp);
  free(XfracData);
  fclose(fp);
}
in main()
{
  load_xfrac(10);
}
