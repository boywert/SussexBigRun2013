#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        PERIODIC\n"
"        PMGRID=512\n"
"        MULTIPLEDOMAINS=8\n"
"        PEANOHILBERT\n"
"        WALLCLOCK\n"
"        MYSORT\n"
"        DOUBLEPRECISION=2\n"
"        DOUBLEPRECISION_FFTW\n"
"        FOF\n"
"        FOF_PRIMARY_LINK_TYPES=2\n"
"        FOF_GROUP_MIN_LEN=20\n"
"        SUBFIND\n"
"        SUBFIND_SAVE_PARTICLELISTS\n"
"        SO_VEL_DISPERSIONS\n"
"        LINKLENGTH=0.2\n"
"        NO_ISEND_IRECV_IN_DOMAIN\n"
"        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG\n"
"\n");
}
