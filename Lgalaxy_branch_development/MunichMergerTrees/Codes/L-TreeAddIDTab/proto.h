
#include "allvars.h"

void load_tree_table(int filenr);
void load_tree(int filenr, int nr);

void free_tree(void);
void free_tree_table(void);


void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
void process_tree(int tree);

void find_progenitors(int tree, int halonr, int prog);
int compare_id(const void *a, const void *b);
int compare_tab(const void *a, const void *b);
void save(int filenr);
