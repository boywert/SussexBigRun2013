#ifndef INC_COSMOLOGY_H
#define INC_COSMOLOGY_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "common.h"

static double h_hubble = 0.70;
static double Omega_0 = 0.73;
static double Omega_m = 0.27;


extern double get_delta_t_in_hubble_unit(double low_redshift, double high_redshitf);

#endif
