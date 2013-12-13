#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "paramparser.h"


#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300



void read_parameter_file(char *fname)
{
#ifndef PROCCOMPATIBLE
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "SnapshotFileBase");
  addr[nt] = SnapshotFileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "FirstSnapShotNr");
  addr[nt] = &FirstSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "LastSnapShotNr");
  addr[nt] = &LastSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "SnapSkipFac");
  addr[nt] = &SnapSkipFac;
  id[nt++] = INT;

  strcpy(tag[nt], "NumberOfOutputFiles");
  addr[nt] = &NumberOfOutputFiles;
  id[nt++] = INT;

  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  *buf = 0;
	  fgets(buf, 200, fd);
	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case DOUBLE:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
	      printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

      i = strlen(OutputDir);
      if(i > 0)
	if(OutputDir[i - 1] != '/')
	  strcat(OutputDir, "/");
    }
  else
    {
      printf("Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	  printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }
#endif

/*start of the proposed code compatible with ProC*/
#ifdef PROCCOMPATIBLE

  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char *tmp;
  int errorFlag = 0;

  printf("PROPOSED VERSION COMPATIBLE WITH PROC!!! \n");

  void *pParamFile;

  pParamFile = openParams(fname);
  if(pParamFile == 0)
    {
      printf("Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "SnapshotFileBase");
  addr[nt] = SnapshotFileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "FilesPerSnapshot");
  addr[nt] = &FilesPerSnapshot;
  id[nt++] = INT;

  strcpy(tag[nt], "LastSnapShotNr");
  addr[nt] = &LastSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "FirstSnapShotNr");
  addr[nt] = &FirstSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &BoxSize;
  id[nt++] = DOUBLE;

#define PROPOSAL1

#ifdef PROPOSAL1
  /*version using minimal number of functions offered by paramparser */
  for(i = 0; i < nt; i++)
    {
      if(getVarString(pParamFile, tag[i], &tmp) != 0)
	{
	  printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	  break;
	}

      switch (id[i])
	{
	case DOUBLE:
	  *((double *) addr[i]) = atof(tmp);
	  printf("Parameter[%s]= %f \n", tag[i], *((double *) addr[i]));
	  break;
	case STRING:
	  strcpy(addr[i], tmp);
	  printf("Parameter[%s]= %s \n", tag[i], addr[i]);
	  break;
	case INT:
	  *((int *) addr[i]) = atoi(tmp);
	  printf("Parameter[%s]= %d \n", tag[i], *((int *) addr[i]));
	  break;
	}
      free(tmp);
    }
  closeParams(pParamFile);


#endif

#ifndef PROPOSAL1
/* PROPOSEDUSE2 */
/* version using all suitable functions offered by paramparser */

  for(i = 0; i < nt; i++)
    {
      if(!isVar(pParamFile, tag[i]))
	{
	  printf("Error: I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	  break;
	}

      switch (id[i])
	{
	case DOUBLE:
	  getVarDouble(pParamFile, tag[i], ((double *) addr[i]));
	  printf("Parameter[%s]= %f \n", tag[i], *((double *) addr[i]));
	  break;

	case STRING:
	  if(getVarStringProt(pParamFile, tag[i], addr[i], 512) != 0)
	    {
	      printf("Error: Value does not fit", tag[i], fname);
	      errorFlag = 1;
	      break;
	    }
	  printf("Parameter[%s]= %s \n", tag[i], addr[i]);
	  break;

	case INT:
	  getVarInt(pParamFile, tag[i], ((int *) addr[i]));
	  printf("Parameter[%s]= %d \n", tag[i], *((int *) addr[i]));
	  break;
	}
    }
  closeParams(pParamFile);
#endif
#undef PROPOSAL1

#endif
/*end of the proposed code compatible with ProC*/

  if(errorFlag)
    exit(1);

}
