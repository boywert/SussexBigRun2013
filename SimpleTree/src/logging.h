#ifndef INC_LOGGIN_H
#define INC_LOGGIN_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "common.h"
#include "readconfig.h"

extern FILE* mpi_log_fp;
extern void logging(char* str);
extern void init_logging();
extern void finalise_logging();

#endif
