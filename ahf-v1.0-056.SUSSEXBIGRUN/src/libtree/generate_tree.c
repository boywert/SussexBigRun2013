#include <stddef.h>
#include <stdio.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

void generate_tree(long unsigned npart, partptr fst_part, double Nth_dom)
{
  int NtreeMin;

  NtreeMin = (int)(Nth_dom+0.5);



}
