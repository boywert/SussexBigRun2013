#include <stdlib.h>

#include <gsl/gsl_math.h>

#include "integrate.h"
#include "network_solver.h"

static const double conv = 1.602177e-12 * 1.0e3 * 6.0221367e23; /* eV2erg * 1.0e3 [keV] * avogadro */

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw)
{
  double sum, xnew;
  int i;

  sum = 0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      sum += x[i];
    }

  if (e)
    {
      for(i = 0; i < nd->nuc_count; i++)
	{
	  xnew = x[i] / sum;
	  *e -= (xnew - x[i]) * nd->nucdata[i].exm * conv;
	  x[i] = xnew;
	}
    }
  else
    {
      for(i = 0; i < nd->nuc_count; i++)
	x[i] /= sum;
    }
}

int network_integrate( double temp, double rho, const double *x, double *dx, double dt, double *dedt, double *drhodt, const struct network_data *nd, struct network_workspace *nw ) {
  double *y;
  double sum;
  int i;

  if (dt == 0 || temp < 1e7) {
    for (i=0; i<nd->nuc_count; i++) dx[i] = 0;
    *dedt = 0;
    if (drhodt) *drhodt = 0;
    return 0;
  }

  /* calculate number densities */
  y = nw->y;
#ifdef NETWORK_SEP_YZ
  y[nd->iYz] = 0.0;
#endif
  for (i=0; i<nd->nuc_count; i++) {
    y[i] = x[i] / nd->nucdata[i].na;
#ifdef NETWORK_SEP_YZ
    y[nd->iYz] += y[i] * gsl_pow_2(nd->nucdata[i].nz);
#endif
  }

#if NETWORK_VARIABLE
  y[nd->iTemp] = temp;
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  y[nd->iRho] = rho;
#endif

  /* run network */
  network_solver_integrate( temp, rho, y, dt, nd, nw );

  /* normalise */
  sum = 0.0;
  for (i=0; i<nd->nuc_count; i++) {
    if (y[i] > 1.0) y[i] = 1.0;
    if (y[i] < 1e-30) y[i] = 1e-30;
    sum += y[i] * nd->nucdata[i].na;
  }
  for (i=0; i<nd->nuc_count; i++) {
    y[i] /= sum;
  }
  /* calculate change of mass fractions and energy release */
  *dedt = 0;
  for (i=0; i<nd->nuc_count; i++) {
    dx[i] = ( y[i] * nd->nucdata[i].na - x[i] ) / dt;
    *dedt -= dx[i] / nd->nucdata[i].na * nd->nucdata[i].exm;
  }
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  if(drhodt) *drhodt = (y[nd->iRho] - rho) / dt;
#else
  if(drhodt) *drhodt = 0.0;
#endif
  *dedt *= conv;

  return 0;
}
