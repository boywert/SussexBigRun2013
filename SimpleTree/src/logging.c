#include "logging.h"

FILE* mpi_log_fp;


void logging(char* str)
{
  fprintf(mpi_log_fp,str);
  fprintf(mpi_log_fp,"\n");
}

void init_logging()
{
  char logfile[1024];
  char command[1024];
  printf("start logging\n");
  sprintf(command,"mkdir -p %s",param_logfolder);
  printf("%s\n",command);
  system(command);
  sprintf(logfile,"%s/status_%d.log",param_logfolder,mpi_rank);
  mpi_log_fp = fopen(logfile, "w");
  if(mpi_log_fp == NULL)
    {
      printf("ERROR: cannot open logfile %s/status_%d.log",param_logfolder,mpi_rank);
      exit(1);
    }
  fprintf(mpi_log_fp,"#################################################\n");
  fprintf(mpi_log_fp,"Start SimpleTree\n\n");
  fprintf(mpi_log_fp,"#################################################\n\n");
}



void finalise_logging()
{
  fclose(mpi_log_fp);
}
