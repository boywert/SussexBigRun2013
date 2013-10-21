#include "network.h"

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw);
int network_integrate( double temp, double rho, const double *x, double *dx, double dt, double *dedt, double *drhodt, const struct network_data *nd, struct network_workspace *nw );
