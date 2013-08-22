#include "hash.h"

ptid_t m_particle_binary_search_and_insert_element_replace_exist(m_particle_wrapper_t* list, m_particle_t add_element)
{
  ptid_t middle,low,high,target;
  char memmgr_buff[memmgr_max_str];
  
  low = 0;
  high = list.npart-1;

  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  if(list.npart == 0)
    {
      list.npart++;
      list.mparticle = memmgr_realloc(list.mparticle,sizeof(m_particle_t),0,memmgr_buff);
      list.mparticle[0].ID = add_element.ID;
      list.mparticle[0].haloID = add_element.haloID;
      return 0;
    }
  else
    {
      high = list.npart-1;
    }  

  while ( low <= high && high < NULLPOINT && low < NULLPOINT)
    {
      middle = ( low + high ) / 2;
 
      if ( add_element.ID == list.mparticle[middle].ID )
	{
	  list.mparticle[middle].ID = add_element.ID;
	  list.mparticle[middle].haloID = add_element.haloID;
	  return list.npart;
	}
      else if ( add_element.ID < list.mparticle[middle].ID )
	{
	  high = middle - 1;
	}
      else
	{
	  low = middle + 1;
	}
    }
  if(high == NULLPOINT)
    {
      target = 0;
    }
  else
    {
      target = low;
    }
  list.npart++;
  list.mparticle = memmgr_realloc(list.mparticle,sizeof(m_particle_t),0,memmgr_buff);
  for(ipart=list.npart-1;ipart > target; ipart--)
    {
      list.mparticle[ipart].ID = list.mparticle[ipart-1].ID;
      list.mparticle[ipart].haloID  =  list.mparticle[ipart-1].haloID;
    }
  list.mparticle[target].ID = add_element.ID;
  list.mparticle[target].haloID = add_element.haloID;
  return list.npart;
}

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

  /* if(searchKey < pool[0] || searchKey > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
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
