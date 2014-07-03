
#include "allvars.h"


void read_parameter_file(char *fname);


void load_tree_table(int filenr, int simulation_number);
void load_tree(int filenr, int treenhalos, int treefirsthalo, int simulation_number);
void load_tree_IDs(int filenr, int treenhalos, int treefirsthalo, int simulation_number);
void load_aux_header(int filenr, int simulation_number);
void load_all_auxdata(int filenr, int simulation_number);

int find_unique_tree_list(int *Unique_TreeList, int *Unique_FileList, int *Unique_TreeCount,
		                      char *Sample_name, int FileNumber, int snap_low, int snap_high, int MAX_UNIQUE_TREES);
