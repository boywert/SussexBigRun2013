#include "cosmology.h"

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
    }
  rieman_sum *= h_hubble*100.;
  return rieman_sum;
}

