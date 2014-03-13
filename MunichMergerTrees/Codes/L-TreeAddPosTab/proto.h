
#include "allvars.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

int compare_bigtab_id(const void *a, const void *b);
int compare_bigtab_index(const void *a, const void *b);
void count_ids(int snapnum);
void load_ids(int snapnum);
void save_ids(int snapnum);

void swap_Nbyte(char *data, int n, int m);
void swap_header(void);
void find_block(char *label, FILE *fd);
int test_block(char *label, FILE *fd);

int compare_tab_id(const void *a, const void *b);
int compare_snap_id(const void *a, const void *b);
void free_tree_ids(void);
void load_tree_ids(int filenr);
void save_tree_pos(int filenr);
int do_match(int snapnum);
void check_missing(void);
void free_pos_data(void);
void load_pos_data(int num, int filenr);
void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
int count_local_missing(int snapnum);

#ifdef HAVE_HDF5
void read_hdf5_attribute(hid_t file_id, 
			 char *obj_name, char *attr_name,
			 hid_t memtype_id, void *buf);
void read_hdf5_header(hid_t file_id, struct io_header *head);
void read_hdf5_dataset(hid_t file_id, char *dset_name, hid_t dtype_id, void *buf);
void read_all_types(hid_t file_id, struct io_header *head, char *dset_name, hid_t dtype_id, void *buf);
#endif
