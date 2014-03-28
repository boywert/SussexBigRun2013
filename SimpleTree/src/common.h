#ifndef INC_COMMON_H
#define INC_COMMON_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include "mpi.h"
#include "readconfig.h"
#define MAXSTRING 1024

#define  pow2(x)   ((x)*(x))
#define  pow3(x)   ((x)*(x)*(x))
#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))


extern const double speed_of_light;
extern const double G;
extern const uint64_t NULLPOINT;
extern const uint64_t MAXUSEABLE;
extern const double max_part_speed_in_c;

/* conversions */
extern const double kpc2m;
extern const double m2Mpc;
extern const double m2km;
extern const double kpc2Mpc;
extern const double Msun2Gadget;
extern const double kpc2Gadget;
extern const double kg2Msun;

extern int mpi_rank;
extern int mpi_nodes;
extern int global_error;


extern void initialise_MPI(int* argc, char ***argv);
extern void finalise_MPI();
#endif
