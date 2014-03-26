#include "logging.h"

void logging(char* str)
{
  fprintf(mpi_log_fp,str);
  fprintf(mpi_log_fp,"\n");
}

void init_logging()
{
  char logfile[1024];
  sprintf(logfile,"%s/status_%d.log",);
  mpi_log_fp = fopen(logfile, "w+");
}

void finalise_logging()
{
  fclose(mpi_log_fp);
}
