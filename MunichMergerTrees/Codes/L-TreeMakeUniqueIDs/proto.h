
#include "allvars.h"

void load_tree_table(int filenr);
void load_tree(int filenr, int nr);

void free_tree(void);
void free_tree_table(void);


void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
void process_tree(int filenr, int tree);

void find_progenitors(int tree, int halonr, int prog);
int compare_id(const void *a, const void *b);
int compare_tab(const void *a, const void *b);


void read_snap_list(void);

void save_tree_dbids(int tree);
void prepare_tree_dbids_output(int filenr);
void close_tree_dbids_output(void);
int walk(int nr);
int peano_hilbert_key(int x, int y, int z, int bits);
