#ifndef INC_COSMOLOGY_H
#define INC_COSMOLOGY_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "common.h"

extern const double h_hubble;
extern const double Omega_0;
extern const double Omega_m;

/* unit (kpc)/(m/s) */
extern double get_delta_t_in_hubble_unit(double low_redshift, double high_redshitf);

#endif
