#include "common.h"

const double speed_of_light = 299792458.; /* m/s */
//const double G = 6.67384e-11; // m^3/(kgs^2);
//G = G*m2Mpc*m2km**2./(Msun2Gadget*kg2Msun)
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
const double kpc2Gadget = 1.0;
const double G = 43018.9720599; // (kpc/h) (km/s)^2 / (1e10Msun/h)
#else
const double kpc2Gadget = 0.001;
const double G = 43.0189720599; // (Mpc/h) (km/s)^2 / (1e10Msun/h)
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
