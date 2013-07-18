#include "hash.h"


int compare_m_halo_t_by_Mvir(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->Mvir < u2->Mvir)
      ret =  -1;
    else if(u1->Mvir > u2->Mvir)
      ret = 1;
    else if(u1->Mvir == u2->Mvir)
      ret = 0;
    return ret;
}

int compare_m_halo_t_by_Mvir_reverse(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->Mvir < u2->Mvir)
      ret =  1;
    else if(u1->Mvir > u2->Mvir)
      ret = -1;
    else if(u1->Mvir == u2->Mvir)
      ret = 0;
    return ret;
}

int compare_m_halo_t_by_ID(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->ID < u2->ID)
      ret =  -1;
    else if(u1->ID > u2->ID)
      ret = 1;
    else if(u1->ID == u2->ID)
      ret = 0;
    return ret;
}

int compare_m_particle_t_by_ID(const void *v1, const void *v2)
{
    const m_particle_t *u1 = v1;
    const m_particle_t *u2 = v2;
    int ret;
    if(u1->ID < u2->ID)
      ret =  -1;
    else if(u1->ID > u2->ID)
      ret = 1;
    else if(u1->ID == u2->ID)
      ret = 0;
    return ret;
}

int compare_uint64_t(const void *v1, const void *v2)
{
  int ret;
  if (*(uint64_t*)v1 < *(uint64_t *)v2)
    ret = -1;
  else if (*(uint64_t *)v1 > *(uint64_t *)v2)
    ret = 1;
  else if (*(uint64_t *)v1 == *(uint64_t *)v2)
    ret = 0;
  return ret;
}

uint64_t search_uint64_t_array( uint64_t searchKey, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  uint64_t *pool = (uint64_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  if(searchKey < pool[0] || searchKey > pool[n_array-1])
    {
      return NULLPOINT;
    }
  while ( low <= high) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchKey == pool[ middle ] )
	{
	  return middle;
	}
      else if ( searchKey < pool[ middle ] )
	{
	  high = middle - 1;
	}
      else
	{
	  low = middle + 1;
	}
    }
  return NULLPOINT;
}
