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

#define LOG_PRINT(...) log_print(__FILE__, __LINE__, __VA_ARGS__ )

extern void init_logging();
#endif
