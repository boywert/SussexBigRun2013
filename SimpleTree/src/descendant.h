#ifndef INC_DESCENDANT_H
#define INC_DESCENDANT_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
    
#include "libio/readsussexbigrun.h"
#include "common.h"
#include "treeformat.h"
#include "cosmology.h"
#include "libhash/hash.h"
extern void make_link_AB(m_halo_wrapper_t* haloA, m_halo_wrapper_t* haloB, double dt);

#endif
