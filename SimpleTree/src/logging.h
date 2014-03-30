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
#include "mpi.h"
#include <stdarg.h>
#include <time.h>

extern void logging(char* str);
extern void init_logging();
extern void finalise_logging();

#define LOG_PRINT(...) log_print(__FILE__, __LINE__, __VA_ARGS__ )
#endif
