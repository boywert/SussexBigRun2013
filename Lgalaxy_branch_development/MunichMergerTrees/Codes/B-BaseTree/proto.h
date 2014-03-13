
#include "allvars.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

void determine_descendants(struct halo_catalogue *catA, struct halo_catalogue *catB, int entry, int snapnum);

int sort_candlist(const void *a, const void *b);
void prepare_index_list(struct halo_catalogue *cat);
void load_subhalo_catalogue(int num, struct halo_catalogue *cat);

void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
#ifdef BACKWARD_CHECKING
void decide_backwards(struct halo_catalogue *catA, struct halo_catalogue *catB);
#endif
void decide_upon_descendant(void);
void count_progenitors(struct halo_catalogue *catA, struct halo_catalogue *catB);

void save_decendant_list(void);

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
void get_id_translation_table(void);
void long_to_str(char *s, long long n);
int sort_IDType(const void *a, const void *b);

int sort_twoids_ord(const void *a, const void *b);
int sort_twoids_id(const void *a, const void *b);
void reassign_ids(MyIDType N, MyIDType *ids); 
void myfree(void *ptr);

#ifdef HAVE_HDF5
void read_hdf5_attribute(hid_t file_id, 
			 char *obj_name, char *attr_name,
			 hid_t memtype_id, void *buf);
void read_hdf5_header(hid_t file_id, struct io_header *head);
void read_hdf5_dataset(hid_t file_id, char *dset_name, hid_t dtype_id, void *buf);


void get_TotNumPart(void);
int attribute_exists(hid_t file_id,
                     char *obj_name, char *attr_name);

#endif

