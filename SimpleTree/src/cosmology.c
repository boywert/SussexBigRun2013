#include "cosmology.h"

const double h_hubble = 0.70;
const double Omega_0 = 0.73;
const double Omega_m = 0.27;

double get_delta_t_in_hubble_unit(double low_redshift, double high_redshift)
{
  double rieman_sum;
  double delta_z = 0.0001;    
  double z;
  rieman_sum = 0.;
  z= low_redshift;
  while(z < high_redshift)
    {
      rieman_sum += delta_z/sqrt(Omega_m*(1.+z)*(1.+z)*(1.+z)*(1.+z)*(1.+z) + Omega_0*(1.+z)*(1.+z) );
      //printf("sum = %lf\n",rieman_sum);
      z += delta_z;
    }
  rieman_sum /= h_hubble*100.;
  return rieman_sum;
}

