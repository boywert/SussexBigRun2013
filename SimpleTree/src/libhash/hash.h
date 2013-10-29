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

typedef struct order_uint64 {
  uint64_t id;
  uint64_t ref;
} order_uint64_t;

extern ptid_t m_particle_binary_search_and_insert_element_replace_exist(m_particle_wrapper_t* list, m_particle_t add_element);
extern int compare_m_halo_t_by_Mvir(const void *v1, const void *v2);
extern int compare_m_halo_t_by_Mvir_reverse(const void *v1, const void *v2);
extern int compare_m_halo_t_by_host_halo(const void *v1, const void *v2);
extern int compare_m_halo_t_by_descendant(const void *v1, const void *v2);
extern int compare_m_halo_t_by_host_halo_reverse(const void *v1, const void *v2);
extern int compare_m_halo_t_by_oriID(const void *v1, const void *v2);
extern int compare_key_sort_t_by_ID(const void *v1, const void *v2);
extern int compare_m_halo_t_by_ID(const void *v1, const void *v2);
extern int compare_m_particle_t_by_ID(const void *v1, const void *v2);
extern int compare_uint64_t(const void *v1, const void *v2);
extern int compare_double(const void *v1, const void *v2);
extern int compare_merit_t_by_merit_delucia2007(const void *v1, const void *v2);
extern int compare_merit_t_by_Mvir(const void *v1, const void *v2);
extern int compare_order_uint64_t_by_ref(const void *v1, const void *v2);

extern uint64_t search_uint64_t_array( uint64_t searchKey, uint64_t n_array ,const void *Array );
extern uint64_t search_m_halo_t_array_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array );
extern uint64_t search_m_halo_t_array_for_oriID( uint64_t searchID, uint64_t n_array ,const void *Array );
extern uint64_t search_particlelist_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array );
extern uint64_t search_m_particle_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array );
extern uint64_t search_key_sort_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array );
extern uint64_t search_order_unint64_t_for_ref( uint64_t searchID, uint64_t n_array ,const void *Array );
#endif
