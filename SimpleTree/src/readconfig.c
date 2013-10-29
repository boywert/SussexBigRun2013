#include "readconfig.h"

long long unsigned param_npart_box;
double  param_boxsize,param_buffer_size;
int  param_domain_per_dim,param_chunk_per_dim,param_chunk_mpi;
char param_INPUTDIR[1024];
char param_OUTPUTDIR[1024];
char param_CHUNKDIR[1024];

struct config
{
  char IDENTIFIER[1024];
  int type;
  int used;
  void *pointer; 
};


void readconfig()
{
  FILE* fp;
  char buffer[1024],ident[1024],value[1024];
  char *pch;
  int count,nconf,i,npart;
  struct config config[100];
  
  /* define identifiers */
  nconf=0;

  config[nconf].type = 4;
  config[nconf].pointer = &(param_CHUNKDIR[0]);
  sprintf(config[nconf].IDENTIFIER,"CHUNKDIR");
  config[nconf].used = 0;
  nconf++;
  
  config[nconf].type = 4;
  config[nconf].pointer = &(param_INPUTDIR[0]);
  sprintf(config[nconf].IDENTIFIER,"INPUTDIR");
  config[nconf].used = 0;
  nconf++;

  config[nconf].type = 4;
  config[nconf].pointer = &(param_OUTPUTDIR[0]);
  sprintf(config[nconf].IDENTIFIER,"OUTPUTDIR");
  config[nconf].used = 0;
  nconf++;
  
  config[nconf].type = 1;
  config[nconf].pointer = &param_domain_per_dim;
  sprintf(config[nconf].IDENTIFIER,"CUBEP3MDOMAINS");
  config[nconf].used = 0;
  nconf++;

  config[nconf].type = 1;
  config[nconf].pointer = &param_chunk_per_dim;
  sprintf(config[nconf].IDENTIFIER,"NCHUNKSPERDIM");
  config[nconf].used = 0;
  nconf++;

  config[nconf].type = 1;
  config[nconf].pointer = &param_chunk_mpi;
  sprintf(config[nconf].IDENTIFIER,"CHUNK_MPI");
  config[nconf].used = 0;
  nconf++;

  config[nconf].type = 2;
  config[nconf].pointer = &npart;
  sprintf(config[nconf].IDENTIFIER,"NPARTSPERDIM");
  config[nconf].used = 0;
  nconf++;

  config[nconf].type = 3;
  config[nconf].pointer = &param_boxsize;
  sprintf(config[nconf].IDENTIFIER,"BOXSIZE");
  config[nconf].used = 0;
  nconf++;

  /* end define identifiers */

  fp = fopen("config","r");
  if(fp == NULL) 
    {
      printf("Cannot open config file\nExit(1)\n");
      exit(1);
    }
  while(fgets(buffer, 1024 ,fp) != NULL)
    {
      count = 0;
      pch = strtok (buffer," \t\n");
      while (pch != NULL)
	{
	  count++;
	  //printf("%s\n",pch);
	  if(count == 1) sprintf(ident,"%s",pch);
	  if(count == 2) sprintf(value, "%s", pch);
	  pch = strtok (NULL, " \t\n");
	}
      if(count > 1)
	{
	  for(i=0;i<nconf;i++)
	    {
	      if(strcmp(ident,config[i].IDENTIFIER) == 0)
		{
		  //printf("found %s = %s\n",config[i].IDENTIFIER,value);
		  if(config[i].type == 1) //int
		    {
		      sscanf(value,"%d",config[i].pointer);
		    }
		  else if(config[i].type == 2) //long long unsigned
		    {
		      sscanf(value,"%llu",config[i].pointer);
		    }
		  else if(config[i].type == 3) //double
		    {
		      sscanf(value,"%lf",config[i].pointer);
		    }
		  else if(config[i].type == 4) //string
		    {
		      sprintf(config[i].pointer,"%s",value);
		    }
		  config[i].used = 1;
		}
	    }
	}
      
    }
  param_npart_box = (long long)npart * (long long)npart * (long long)npart;
  printf("npart_box = %llu\n",npart_box);
  printf("boxsize = %lf\n",boxsize);
  printf("domain_per_dim = %d\n",domain_per_dim);
  printf("INPUTDIR = %s\n",INPUTDIR);
  printf("OUTPUTDIR = %s\n",OUTPUTDIR);

  fclose(fp);
}
