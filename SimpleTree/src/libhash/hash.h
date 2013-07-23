#ifndef INC_HASH_H
#define INC_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "../common.h"
#include "../treeformat.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
ptid_t m_particle_binary_search_and_insert_element_replace_exist(m_particle_wrapper_t* list, m_particle_t add_element);
extern int compare_m_halo_t_by_Mvir(const void *v1, const void *v2);
extern int compare_m_halo_t_by_Mvir_reverse(const void *v1, const void *v2);
extern int compare_m_halo_t_by_host_halo(const void *v1, const void *v2);
extern int compare_m_halo_t_by_ID(const void *v1, const void *v2);
extern int compare_m_particle_t_by_ID(const void *v1, const void *v2);
extern int compare_uint64_t(const void *v1, const void *v2);
extern uint64_t search_uint64_t_array( uint64_t searchKey, uint64_t n_array ,const void *Array );

#endif
