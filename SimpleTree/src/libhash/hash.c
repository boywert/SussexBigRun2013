#include "hash.h"

ptid_t m_particle_binary_search_and_insert_element_replace_exist(m_particle_wrapper_t* list, m_particle_t add_element)
{
  ptid_t middle,low,high,target,ipart;
  char memmgr_buff[memmgr_max_str];
  
  low = 0;
  high = list->npart-1;

  sprintf(memmgr_buff,"Particle inside wrapper: Hash");
  if(list->npart == 0)
    {
      list->npart++;
      list->mparticle = memmgr_realloc(list->mparticle,sizeof(m_particle_t),0,memmgr_buff);
      list->mparticle[0].ID = add_element.ID;
      list->mparticle[0].haloID = add_element.haloID;
      return list->npart;
    }
  else
    {
      high = list->npart-1;
    }  

  while ( low <= high && high < NULLPOINT && low < NULLPOINT)
    {
      middle = ( low + high ) / 2;
 
      if ( add_element.ID == list->mparticle[middle].ID )
	{
	  list->mparticle[middle].ID = add_element.ID;
	  list->mparticle[middle].haloID = add_element.haloID;
	  printf("Duplicate: change\n");
	  return list->npart;
	}
      else if ( add_element.ID < list->mparticle[middle].ID )
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
  list->npart++;
  list->mparticle = memmgr_realloc(list->mparticle,sizeof(m_particle_t)*list->npart,sizeof(m_particle_t)*(list->npart-1),memmgr_buff);
  for(ipart=list->npart-1;ipart > target; ipart--)
    {
      list->mparticle[ipart].ID = list->mparticle[ipart-1].ID;
      list->mparticle[ipart].haloID  =  list->mparticle[ipart-1].haloID;
    }
  list->mparticle[target].ID = add_element.ID;
  list->mparticle[target].haloID = add_element.haloID;
  return list->npart;
}

int compare_order_uint64_t_by_ref(const void *v1, const void *v2)
{
    const order_uint64_t *u1 = v1;
    const order_uint64_t *u2 = v2;
    int ret;
    if(u1->ref < u2->ref)
      ret = -1;
    else if(u1->ref > u2->ref)
      ret = 1;
    else if(u1->ref == u2->ref)
      ret = 0;
    return ret;
}

int compare_m_halo_t_by_host_halo(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->host_halo < u2->host_halo)
      ret =  -1;
    else if(u1->host_halo > u2->host_halo)
      ret = 1;
    else if(u1->host_halo == u2->host_halo)
      ret = 0;
    return ret;
}

int compare_m_halo_t_by_descendant(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->descendant < u2->descendant)
      ret =  -1;
    else if(u1->descendant > u2->descendant)
      ret = 1;
    else if(u1->descendant == u2->descendant)
      ret = 0;
    return ret;
}

int compare_m_halo_t_by_host_halo_reverse(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->host_halo < u2->host_halo)
      ret =  1;
    else if(u1->host_halo > u2->host_halo)
      ret = -1;
    else if(u1->host_halo == u2->host_halo)
      ret = 0;
    return ret;
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
int compare_key_sort_t_by_ID(const void *v1, const void *v2)
{
    const key_sort_t *u1 = v1;
    const key_sort_t *u2 = v2;
    int ret;
    if(u1->ID < u2->ID)
      ret =  -1;
    else if(u1->ID > u2->ID)
      ret = 1;
    else if(u1->ID == u2->ID)
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

int compare_m_halo_t_by_oriID(const void *v1, const void *v2)
{
    const m_halo_t *u1 = v1;
    const m_halo_t *u2 = v2;
    int ret;
    if(u1->oriID < u2->oriID)
      ret =  -1;
    else if(u1->oriID > u2->oriID)
      ret = 1;
    else if(u1->oriID == u2->oriID)
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

int compare_merit_t_by_merit_delucia2007(const void *v1, const void *v2)
{
    const merit_t *u1 = v1;
    const merit_t *u2 = v2;
    int ret;
    if(u1->merit_delucia2007 < u2->merit_delucia2007)
      ret =  -1;
    else if(u1->merit_delucia2007 > u2->merit_delucia2007)
      ret = 1;
    else if(u1->merit_delucia2007 == u2->merit_delucia2007)
      ret = 0;
    return ret;
}

int compare_merit_t_by_Mvir(const void *v1, const void *v2)
{
    const merit_t *u1 = v1;
    const merit_t *u2 = v2;
    int ret;
    if(u1->Mvir < u2->Mvir)
      ret =  -1;
    else if(u1->Mvir > u2->Mvir)
      ret = 1;
    else if(u1->Mvir == u2->Mvir)
      ret = 0;
    return ret;
}

int compare_double(const void *v1, const void *v2)
{
  int ret;
  if (*(double*)v1 < *(double *)v2)
    ret = -1;
  else if (*(double *)v1 > *(double *)v2)
    ret = 1;
  else if (*(double *)v1 == *(double *)v2)
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


uint64_t search_m_halo_t_array_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  m_halo_t *pool = (m_halo_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[0] || searchID > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchID == pool[ middle ].ID )
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].ID )
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

uint64_t search_m_halo_t_array_for_oriID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  m_halo_t *pool = (m_halo_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[0] || searchID > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchID == pool[ middle ].oriID )
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].oriID )
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

uint64_t search_m_particle_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high,ipart;
  m_particle_t *pool = (m_particle_t *) Array;
  /* for(ipart=0; ipart<n_array; ipart++ ) */
  /*   { */
  /*     printf("ipart: %llu => %llu\n",ipart,pool[ipart].ID); */
  /*   } */
  /* exit(0); */
  //printf("start searching for %llu in %llu element array\n",searchID,n_array);
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[0] || searchID > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
      //      printf("%llu : %llu < %llu : %llu < %llu : %llu\n",low,pool[low].ID,middle,pool[middle].ID,high,pool[high].ID);
      if ( searchID == pool[ middle ].ID )
	{
	  return pool[middle].haloID;
	}
      else if ( searchID < pool[ middle ].ID )
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



uint64_t search_particlelist_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  particlelist_t *pool = (particlelist_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[0] || searchID > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchID == pool[ middle ].ID )
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].ID )
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


uint64_t search_key_sort_t_for_ID( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  key_sort_t *pool = (key_sort_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[0] || searchID > pool[n_array-1]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchID == pool[ middle ].ID )
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].ID )
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


uint64_t search_order_unint64_t_for_ref( uint64_t searchID, uint64_t n_array ,const void *Array )
{
  uint64_t middle,low,high;
  order_uint64_t *pool = (order_uint64_t *) Array;
  //printf("start search\n");
  /*
    for(i =0; i< n_array; i++)
    {
    printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  /* if(searchID < pool[low] || searchID > pool[high]) */
  /*   { */
  /*     return NULLPOINT; */
  /*   } */
  while ( low <= high && high < NULLPOINT) 
    {
      middle = ( low + high ) / 2;
 
      if ( searchID == pool[ middle ].ref )
	{
	  return middle;
	}
      else if ( searchID < pool[ middle ].ref )
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
