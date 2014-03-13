
#include "allvars.h"



void generate_trees(void);

void load_subhalo_catalogue(int num);
void count_halos(void);

void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
void set_progenitor_pointers(void);

int peano_hilbert_key(int x, int y, int z, int bits);
int whichfile(float *pos);

size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);

