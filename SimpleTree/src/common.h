#ifndef INC_COMMON_H
#define INC_COMMON_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "mpi.h"
#include "readconfig.h"
#define MAXSTRING 1024

extern const double speed_of_light;
extern const uint64_t NULLPOINT;
extern const uint64_t MAXUSEABLE;
extern const double max_part_speed_in_c;
extern const double kpc2m;
extern int mpi_rank;
extern int mpi_nodes;
extern int global_error;
extern void initialise_MPI(int* argc, char ***argv);
extern void finalise_MPI();
#endif
