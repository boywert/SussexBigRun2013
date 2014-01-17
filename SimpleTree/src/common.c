#include "common.h"

const double speed_of_light = 299792458.; /* m/s */
const uint64_t NULLPOINT = (uint64_t)-1;
const uint64_t MAXUSEABLE = (uint64_t)-2;
const double max_part_speed_in_c = 0.01;

/* conversion */
const double kpc2m = 3.08567758e19;
const double kpc2Mpc = 0.001;
const double m2Mpc = 1./3.08567758e22;
const double m2km = 0.001;
const double Msun2Gadget = 1.e-10;
#ifdef GADGETKPC
double kpc2Gadget = 1.0;
#else
double kpc2Gadget = kpc2Mpc;
#endif
const double kg2Msun = 1.989e-30;



int global_error;
int mpi_rank,mpi_nodes;

void initialise_MPI(int* argc, char ***argv )
{
  MPI_Init  (argc,argv);
  MPI_Comm_size ( MPI_COMM_WORLD, &mpi_nodes);
  MPI_Comm_rank ( MPI_COMM_WORLD, &mpi_rank );
}

void finalise_MPI()
{
  MPI_Finalize();
}
