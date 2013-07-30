#include "common.h"

const double speed_of_light = 299792458.; /* m/s */
const uint64_t NULLPOINT = (uint64_t)-1;
const uint64_t MAXUSEABLE = (uint64_t)-2;
const double max_part_speed_in_c = 0.01;

void initialise_MPI(int argc,char **argv)
{
  MPI_Init ( &argc, &argv );
  MPI_Comm_size ( MPI_COMM_WORLD, &mpi_nodes);
  MPI_Comm_rank ( MPI_COMM_WORLD, &mpi_rank );
}

void finalise_MPI()
{
  MPI_Finalize();
}
