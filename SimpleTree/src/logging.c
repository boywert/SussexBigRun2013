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
  sprintf(command,"mkdir -p %s",param_logfolder);
  system(command);
  sprintf(logfile,"%s/status_%d.log",param_logfolder,mpi_rank);
  mpi_log_fp = fopen(logfile, "w+");
  if(mpi_log_fp == NULL)
    {
      printf("ERROR: cannot open logfile %s/status_%d.log",param_logfolder,mpi_rank);
    }
}

void finalise_logging()
{
  fclose(mpi_log_fp);
}
