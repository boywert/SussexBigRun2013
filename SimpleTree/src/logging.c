#include "logging.h"

char logfile[1024];
static int SESSION_TRACKER;
FILE *fp;
char* print_time();

char* print_time()
{
    time_t t;
    char *buf;
    
    time(&t);
    buf = (char*)malloc(strlen(ctime(&t))+ 1);
    
    snprintf(buf,strlen(ctime(&t)),"%s ", ctime(&t));
   
    return buf;
}

void log_print(char* filename, int line, char *fmt,...)
{
  va_list         list;
  char            *p, *r;
  int             e;

  if(SESSION_TRACKER > 0)
    fp = fopen (logfile,"a+");
  else
    fp = fopen (logfile,"w");
  if(fp == NULL)
    {
      printf("ERROR: Cannot open log file: %s\n",logfile);
    }
  fprintf(fp,"%s|:  ",print_time());
  va_start( list, fmt );

  for ( p = fmt ; *p ; ++p )
    {
      if ( *p != '%' )//If simple string
        {
	  fputc( *p,fp );
        }
      else
        {
	  switch ( *++p )
            {
	      /* string */
            case 's':
	      {
                r = va_arg( list, char * );

                fprintf(fp,"%s", r);
                continue;
	      }

	      /* integer */
            case 'd':
	      {
                e = va_arg( list, int );

                fprintf(fp,"%d", e);
                continue;
	      }

            default:
	      fputc( *p, fp );
            }
        }
    }
  va_end( list );
  fprintf(fp,"  [%s][line: %d] ",filename,line);
  fputc( '\n', fp );
  SESSION_TRACKER++;
  fclose(fp);
}

void init_logging()
{
  char command[1024];
  printf("start logging\n");
  sprintf(command,"mkdir -p %s",param_logfolder);
  printf("%s\n",command);
  system(command);
  sprintf(logfile,"%s/status_%d.log",param_logfolder,mpi_rank);
  LOG_PRINT("Start SimpleTree");
}


