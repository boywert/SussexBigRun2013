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
  sprintf(logfile,"%s/status_%d.log",param_logfolder,mpi_rank);
  mpi_log_fp = fopen(logfile, "w+");
}

void finalise_logging()
{
  fclose(mpi_log_fp);
}
