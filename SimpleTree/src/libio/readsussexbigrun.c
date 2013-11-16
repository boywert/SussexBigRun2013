#include "readsussexbigrun.h"
/* Some private function for this file only */
int compare_make_catalogue_halo_t_by_Mvir_reverse(const void *v1, const void *v2);
void AHF_alloc_profiles( uint32_t nbins, halo_profile_t *prof);
void AHF_free_profiles(halo_profile_t *prof);
/* End private function */
void sussexbigrun_dm_outputs( m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  sprintf(foldername,"%s/%3.3f",outputfolder,haloB->redshift);
  sprintf(command,"mkdir -p %s", foldername);
  system(command);
  sprintf(filename,"%s/%3.3f/dmdt_%d.dat",outputfolder,haloB->redshift,domainid);
  sprintf(command,"rm -f %s",filename);
  system(command);
 
  fp = fopen(filename, "w+");
  for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
    {
      /* if(haloB->mhalos[ihalo].used == 1) */
      /* 	{ */
      fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%llu\t%d\n",
	      haloB->mhalos[ihalo].Xc,
	      haloB->mhalos[ihalo].Yc,
	      haloB->mhalos[ihalo].Zc,
	      haloB->mhalos[ihalo].Mvir,
	      haloB->mhalos[ihalo].Rvir,
	      haloB->mhalos[ihalo].dm_dt,
	      haloB->mhalos[ihalo].oriID,
	      haloB->mhalos[ihalo].domainID);
      /* } */
    }
  fclose(fp);

}

m_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary(char *folder, float redshift, int tot_domain )
{
  FILE *fphalo,*fppart;
  char halofile[MAXSTRING],partfile[MAXSTRING];
  //uint64_t one;
  //uint32_t onep;
  m_halo_wrapper_t mhalo;
  int i;
  hid_t ihalo;
  mhalo.nHalos = 0;
  mhalo.redshift = redshift;
  mhalo.mhalos= memmgr_malloc(0,"Halo Array");
  //tot_domain = 10;
  for (i=0;i<tot_domain;i++)
    {
      sprintf(halofile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_halos.dat_bin",folder,redshift,i);
      sprintf(partfile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_pids.dat_bin",folder,redshift,i);
      fphalo = fopen(halofile,"rb");
      fppart = fopen(partfile,"rb");
      mhalo = sussexbigrun_read_AHF_binary(fphalo, fppart, i, mhalo);
      fclose(fphalo);
      fclose(fppart);
    }
  /* for(ihalo=0;ihalo<mhalo.nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",mhalo.mhalos[ihalo].ID,mhalo.mhalos[ihalo].Mvir); */
  /*   } */
  //memmgr_printdetails();
  mhalo = sussexbigrun_filterhalos_and_particles(mhalo);
  //memmgr_printdetails();
  return mhalo;
}


m_halo_wrapper_t sussexbigrun_add_halo_buffer_binary(char *folder, float redshift, int domain, double domain_width, int domain_per_dim, double buffer_width, int position, m_halo_wrapper_t mhalo_ori)
{
  FILE *fphalo,*fppart;
  char halofile[MAXSTRING],partfile[MAXSTRING];
  uint64_t old,new;
  m_halo_wrapper_t mhalo;
  //uint32_t onep;
  int i,block_x,block_y,block_z;
  hid_t ihalo,tot_halos;
  double min_x,max_x,min_y,max_y,min_z,max_z;
  char memmgr_buff[memmgr_max_str];
  //tot_domain = 10;
  mhalo.nHalos = 0;
  mhalo.redshift = redshift;
  mhalo.mhalos= memmgr_malloc(0,"Halo Array");
  i = domain;
  sprintf(halofile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_halos.dat_bin",folder,redshift,i);
  sprintf(partfile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_pids.dat_bin",folder,redshift,i);
  fphalo = fopen(halofile,"rb");
  fppart = fopen(partfile,"rb");
  if(fppart && fphalo)
    {
      mhalo = sussexbigrun_read_AHF_binary(fphalo, fppart, i, mhalo);
      fclose(fphalo);
      fclose(fppart);
    }
  /* for(ihalo=0;ihalo<mhalo.nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",mhalo.mhalos[ihalo].ID,mhalo.mhalos[ihalo].Mvir); */
  /*   } */
  //memmgr_printdetails();

  mhalo = sussexbigrun_filterhalos_and_particles(mhalo);
  //memmgr_printdetails();
  block_z = (int) (domain/(domain_per_dim*domain_per_dim));
  block_y = (int)((domain - block_z*(domain_per_dim * domain_per_dim))/domain_per_dim);
  block_x = (int)(domain - block_z*(domain_per_dim*domain_per_dim) - block_y*domain_per_dim);
  
  min_x = block_x*domain_width+buffer_width;

  max_x = (block_x+1)*domain_width-buffer_width;
  min_y = block_y*domain_width+buffer_width;
  max_y = (block_y+1)*domain_width-buffer_width;
  min_z = block_z*domain_width+buffer_width;
  max_z = (block_z+1)*domain_width-buffer_width;

  for(ihalo = 0; ihalo < mhalo.nHalos; ihalo++ )
    {
      if(position == 1) //-x
	{
	  if(mhalo.mhalos[ihalo].Xc > min_x)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
      else if(position == 2) //+x
	{
	  if(mhalo.mhalos[ihalo].Xc < max_x)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
      else if(position == 3) //-y
	{
	  if(mhalo.mhalos[ihalo].Yc > min_y)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
      else if(position == 4) //+y
	{
	  if(mhalo.mhalos[ihalo].Yc < max_y)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
      else if(position == 5) //-z
	{
	  if(mhalo.mhalos[ihalo].Zc > min_z)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
      else if(position == 6) //+z
	{
	  if(mhalo.mhalos[ihalo].Zc > max_z)
	    mhalo.mhalos[ihalo].ID = NULLPOINT;
	}
    }
  qsort(mhalo.mhalos,mhalo.nHalos, sizeof(m_halo_t), compare_m_halo_t_by_ID);
  tot_halos = 0;
  for(ihalo=0;ihalo<mhalo.nHalos;ihalo++)
    {
      //printf("ihalo = %llu host =%llu\n",ihalo,mhalo.mhalos[ihalo].host_halo);
      if(mhalo.mhalos[ihalo].ID == NULLPOINT)
	break;
      else
	{
	  tot_halos++;
	}
    }
  sprintf(memmgr_buff,"Particle: Halo Array");
  for(ihalo=tot_halos; ihalo<mhalo.nHalos; ihalo++)
    {
      memmgr_free(mhalo.mhalos[ihalo].Particles, mhalo.mhalos[ihalo].npart*sizeof(particlelist_t), memmgr_buff);
    }
  old = mhalo.nHalos*sizeof(m_halo_t);
  new = tot_halos*sizeof(m_halo_t);
  sprintf(memmgr_buff,"Halo Array");
  mhalo.mhalos = memmgr_realloc(mhalo.mhalos,new,old, memmgr_buff);
  mhalo.nHalos = tot_halos;
  for(ihalo=0;ihalo<mhalo.nHalos;ihalo++)
    {
      mhalo.mhalos[ihalo].used = 0;
    }
  old =  mhalo_ori.nHalos*sizeof(m_halo_t);
  tot_halos = mhalo_ori.nHalos+mhalo.nHalos;
  new =  tot_halos*sizeof(m_halo_t);
  mhalo_ori.mhalos = memmgr_realloc(mhalo_ori.mhalos,new,old, memmgr_buff);
  for(ihalo = mhalo_ori.nHalos; ihalo < tot_halos; ihalo++)
    {
      mhalo_ori.mhalos[ihalo] = mhalo.mhalos[ihalo-mhalo_ori.nHalos];
    }
  mhalo.mhalos = memmgr_realloc(mhalo.mhalos,0,mhalo.nHalos*sizeof(m_halo_t), memmgr_buff);
  return mhalo_ori;
}


m_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary_single_domain_include_buffer(char *folder, float redshift, int domain, int domain_per_dim, double domain_width, double dx)
{
  FILE *fphalo,*fppart;
  char halofile[MAXSTRING],partfile[MAXSTRING];
  //uint64_t one;
  //uint32_t onep;
  m_halo_wrapper_t mhalo;
  int i,j,k,x,y,z,block,position,block_x,block_y,block_z;
  hid_t ihalo;
  double fixed_buffer = 400.0; //400 kpc/h buffer
  
  block_z = (int) (domain/(domain_per_dim*domain_per_dim));
  block_y = (int)((domain - block_z*(domain_per_dim * domain_per_dim))/domain_per_dim);
  block_x = (int)(domain - block_z*(domain_per_dim*domain_per_dim) - block_y*domain_per_dim);
  
  mhalo = sussexbigrun_load_halo_catalogue_binary_single_domain(folder,redshift,domain);
  for(i=-1;i<=1;i++)
    {
      for(j=-1;j<=1;j++)
  	{
  	  for(k=-1;k<=1;k++)
  	    {
  	      x = (block_x + i)%domain_per_dim;
  	      y = (block_y + j)%domain_per_dim;
  	      z = (block_z + k)%domain_per_dim;
  	      if(i != 0 && j != 0 && k != 0)
  		{
  		  block = z*(domain_per_dim*domain_per_dim) + y*domain_per_dim + x;
  		  if(i==1)
  		    position = 1;
  		  if(i==-1)
  		    position = 2;
  		  if(j==1)
  		    position = 3;
  		  if(j==-1)
  		    position = 4;
  		  if(k==1)
  		    position = 5;
  		  if(k==-1)
  		    position = 6;
  		  mhalo = sussexbigrun_add_halo_buffer_binary(folder, redshift, domain, domain_width, domain_per_dim, dx+fixed_buffer, position, mhalo);
  		}
  	    }
  	}
    }
  /* for(ihalo=0;ihalo<mhalo.nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",mhalo.mhalos[ihalo].ID,mhalo.mhalos[ihalo].Mvir); */
  /*   } */
  //memmgr_printdetails();
  mhalo = sussexbigrun_filterhalos_and_particles(mhalo);
  //memmgr_printdetails();
  return mhalo;
}

m_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary_single_domain(char *folder, float redshift, int domain )
{
  FILE *fphalo,*fppart;
  char halofile[MAXSTRING],partfile[MAXSTRING];
  //uint64_t one;
  //uint32_t onep;
  m_halo_wrapper_t mhalo;
  int i;
  hid_t ihalo;
  mhalo.nHalos = 0;
  mhalo.redshift = redshift;
  mhalo.mhalos= memmgr_malloc(0,"Halo Array");
  //tot_domain = 10;
  i = domain;
  sprintf(halofile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_halos.dat_bin",folder,redshift,i);
  sprintf(partfile,"%s/%2.3f_AHF_halos_cubepm_domain_%d_pids.dat_bin",folder,redshift,i);
  fphalo = fopen(halofile,"rb");
  fppart = fopen(partfile,"rb");
  if(fphalo && fphalo)
    {
      mhalo = sussexbigrun_read_AHF_binary(fphalo, fppart, i, mhalo);
      fclose(fphalo);
      fclose(fppart);
    }
  /* for(ihalo=0;ihalo<mhalo.nHalos;ihalo++) */
  /*   { */
  /*     printf("%ld => %f\n",mhalo.mhalos[ihalo].ID,mhalo.mhalos[ihalo].Mvir); */
  /*   } */
  //  memmgr_printdetails();
  mhalo = sussexbigrun_filterhalos_and_particles(mhalo);
  //memmgr_printdetails();
  for(ihalo=0;ihalo<mhalo.nHalos;ihalo++)
    {
      mhalo.mhalos[ihalo].used = 1;
    }
    // memmgr_printdetails();
  return mhalo;
}



m_halo_wrapper_t sussexbigrun_filterhalos_and_particles(m_halo_wrapper_t mhalo)
{
  ptid_t ipart,countpart,ref,p_target,jpart;
  hid_t ihalo,tot_halos,h_target,hostid;
  uint64_t old,new;
  char memmgr_buff[memmgr_max_str];
  // m_particle_wrapper_t *tmp;
  //key_sort_t *key_for_sort;
#ifdef TOPLEVELONLY
  tot_halos = 0;  
  //printf("Filter halos and particles\n");
  qsort(mhalo.mhalos,mhalo.nHalos, sizeof(m_halo_t), compare_m_halo_t_by_host_halo_reverse);
  for(ihalo=0;ihalo<mhalo.nHalos;ihalo++)
    {
      //printf("ihalo = %llu host =%llu\n",ihalo,mhalo.mhalos[ihalo].host_halo);
      if(mhalo.mhalos[ihalo].host_halo < NULLPOINT)
	break;
      else
	{
	  tot_halos++;
	}
    }
  sprintf(memmgr_buff,"Particle: Halo Array");
  for(ihalo=tot_halos; ihalo<mhalo.nHalos; ihalo++)
    {
      memmgr_free(mhalo.mhalos[ihalo].Particles, mhalo.mhalos[ihalo].npart*sizeof(particlelist_t), memmgr_buff);
    }
  old = mhalo.nHalos*sizeof(m_halo_t);
  new = tot_halos*sizeof(m_halo_t);
  sprintf(memmgr_buff,"Halo Array");
  mhalo.mhalos = memmgr_realloc(mhalo.mhalos,new,old, memmgr_buff);
  mhalo.nHalos = tot_halos;
#else //ifndef TOPLEVELONLY => Find hosthalo ID
  qsort(mhalo.mhalos,mhalo.nHalos, sizeof(m_halo_t), compare_m_halo_t_by_oriID);
  for(ihalo=0;ihalo<mhalo.nHalos;ihalo++)
    {
      if(mhalo.mhalos[ihalo].host_halo < NULLPOINT)
	{
	  hostid = search_m_halo_t_array_for_oriID( mhalo.mhalos[ihalo].host_halo,mhalo.nHalos,mhalo.mhalos);
	  if(hostid == NULLPOINT) printf("hostid = %llu <= %llu\n",(hid_t)mhalo.mhalos[ihalo].host_halo, (hid_t)hostid);
	}
    } 
#endif //TOPLEVELONLY


  /* sprintf(memmgr_buff,"Particle Wrapper: Hash"); */
  /* tmp = memmgr_malloc(sizeof(m_particle_wrapper_t),memmgr_buff); */
  /* tmp[0].npart = 0; */
  /* sprintf(memmgr_buff,"Particle inside wrapper: Hash"); */
  /* tmp[0].mparticle = memmgr_malloc(0,memmgr_buff); */
  /* qsort(mhalo.mhalos,mhalo.nHalos, sizeof(m_halo_t),compare_m_halo_t_by_Mvir_reverse); */
  /* countpart = 0; */
  /* for(ihalo=0;ihalo < mhalo.nHalos; ihalo++) */
  /*   { */
  /*     //printf("ihalo: %llu\n",ihalo); */

  /*     tmp[0].npart += mhalo.mhalos[ihalo].npart; */
  /*     //printf("tot npart:%llu\n",tmp[0].npart); */
  /*     tmp[0].mparticle = memmgr_realloc(tmp[0].mparticle,sizeof(m_particle_t)*tmp[0].npart,sizeof(m_particle_t)*(tmp[0].npart-mhalo.mhalos[ihalo].npart),memmgr_buff); */
  /*     for(ipart=0;ipart<mhalo.mhalos[ihalo].npart;ipart++) */
  /* 	{ */
  /* 	  //printf("ihalo: %llu ipart:%llu\n",ihalo,ipart); */
  /* 	  //printf("ipart = %llu, ihalo = %llu\n",ipart,ihalo); */
  /* 	  tmp[0].mparticle[countpart].ID =  mhalo.mhalos[ihalo].Particles[ipart].ID; */
  /* 	  tmp[0].mparticle[countpart].haloID =  mhalo.mhalos[ihalo].ID; */
  /* 	  countpart++; */
  /* 	  //insert.ID = mhalo.mhalos[ihalo].Particles[ipart].ID; */
  /* 	  //printf("pid = %llu/%llu :%llu\n",ipart,mhalo.mhalos[ihalo].npart,ihalo); */
  /* 	  //m_particle_binary_search_and_insert_element_replace_exist(&(tmp[0]), insert); */
  /* 	} */
  /*     //memmgr_printdetails(); */
  /*   } */
  /* qsort(tmp[0].mparticle,tmp[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID); */
  /* ref = tmp[0].mparticle[0].ID; */
  /* countpart = 0; */
  /* qsort(mhalo.mhalos,mhalo.nHalos, sizeof(m_halo_t),compare_m_halo_t_by_ID); */
  /* sprintf(memmgr_buff,"TMP particles: Hash"); */
  /* for(ipart=1;ipart<tmp[0].npart;ipart++) */
  /*   { */
  /*     if(tmp[0].mparticle[ipart].ID == ref) */
  /* 	{ */
  /* 	  /\* ihalo = tmp[0].mparticle[ipart].haloID; *\/ */
  /* 	  /\* h_target = search_m_halo_t_array_for_ID( ihalo, mhalo.nHalos , mhalo.mhalos); *\/ */
  /* 	  /\* //printf("search = %llu\n",h_target); *\/ */
  /* 	  /\* key_for_sort = memmgr_malloc(mhalo.mhalos[h_target].npart*sizeof(key_sort_t),memmgr_buff); *\/ */
  /* 	  /\* for(jpart=0;jpart<mhalo.mhalos[h_target].npart;jpart++) *\/ */
  /* 	  /\*   { *\/ */
  /* 	  /\*     key_for_sort[jpart].ID = mhalo.mhalos[h_target].Particles[jpart].ID; *\/ */
  /* 	  /\*     key_for_sort[jpart].order = jpart; *\/ */
  /* 	  /\*   } *\/ */
  /* 	  /\* qsort(key_for_sort,mhalo.mhalos[h_target].npart, sizeof(key_sort_t),compare_m_halo_t_by_ID); *\/ */
  /* 	  /\* p_target = search_key_sort_t_for_ID( tmp[0].mparticle[ipart].ID, mhalo.mhalos[h_target].npart , key_for_sort); *\/ */
  /* 	  /\* p_target = key_for_sort[p_target].order; *\/ */
  /* 	  /\* memmgr_free(key_for_sort,mhalo.mhalos[h_target].npart*sizeof(key_sort_t),memmgr_buff); *\/ */
  /* 	  /\* //p_target = search_particlelist_t_for_ID(); *\/ */
  /* 	  tmp[0].mparticle[ipart].ID = NULLPOINT; */
  /* 	} */
  /*     else */
  /* 	{ */
  /* 	  ref = tmp[0].mparticle[ipart].ID; */
  /* 	} */
  /* 	//printf("dupplicate pid\n"); */
  /*   } */
  /* //qsort(tmp[0].mparticle,tmp[0].npart, sizeof(m_particle_t),compare_m_particle_t_by_ID); */
  /* printf("total duplicate : %llu\n",countpart); */
  return mhalo;
}

m_halo_wrapper_t sussexbigrun_read_AHF_binary(FILE *fphalo, FILE *fppart, int domain, m_halo_wrapper_t mhalo)
{
  uint64_t numHalos,counthalo,counthalo_local;
  //uint64_t numHaloFromPartFile;
  //uint32_t numColumns;
  uint64_t i,size;
  //int32_t  one;
  int      swap=0,flag;
  halo_t   halo;
  ptid_t ipart,npart;
  //ptid_t id;
  size_t old,new;
  char memmgr_buff[memmgr_max_str];
  struct particle_buffer *pid_buff;;
  if(fphalo == NULL || fppart == NULL) 
    {
      printf("Cannot open file to read\n");
      exit(1);
    }
  // figure out swap status
  /* fread(&one, sizeof(int32_t), 1, fphalo); */
  /* if(one == 1)    */
  /*   swap = 0; */
  /* else */
  /*   swap = 1; */
  /* printf("halo: one = %d\n",one); */
  /* // figure out swap status */
  /* fread(&one, sizeof(int32_t), 1, fppart); */
  /* if(one == 1)    */
  /*    swap = 0; */
  /* else  */
  /*   swap = 1; */
  /* printf("part: one = %d\n",one); */
  swap = 0;
  /* ReadULong(fppart, &numHaloFromPartFile,   swap); */
  /* ReadUInt (fppart, &numColumns, swap); */
  /* printf("particlefile: nhalo = %llu\n",numHaloFromPartFile); */
  /* ReadULong(fphalo, &numHalos,   swap); */
  /* ReadUInt (fphalo, &numColumns, swap);  */
  /* printf("halofile: nhalo = %llu\n",numHalos); */
  /* exit(0); */


  // find the size of the file to calculate the total number of halos
  fseek(fphalo, 0L, SEEK_END);
  size = (uint64_t) ftell(fphalo);
  fseek(fphalo, 0L, SEEK_SET);
  numHalos = size/halo_t_size;
  //printf("size = %f; sizestruct %f; totalhalo %f\n",(float)size,(float)halo_t_size, (float)size/halo_t_size);
  //printf("halo in domain %d :%ld \n",domain,numHalos);
  sprintf(memmgr_buff,"Halo Array");

  counthalo = mhalo.nHalos;
  counthalo_local = 0;
  old = mhalo.nHalos*sizeof(m_halo_t);
  new = old + numHalos*sizeof(m_halo_t);
  mhalo.nHalos += numHalos;
  new = mhalo.nHalos*sizeof(m_halo_t);
  //printf("old:%ld new:%ld totalN:%ld",old,new,mhalo.nHalos);
  mhalo.mhalos = memmgr_realloc(mhalo.mhalos,new,old, memmgr_buff);
  //printf("finish realloc\n");
  // read in halo properties
  //counthalo = 0;
  flag = 0;
  for(i=0; i<numHalos; i++) 
    {
      ReadULong(fphalo, &halo.ID,           swap);    // ID(1)
      ReadULong(fphalo, &halo.hostHalo,     swap);    // hostHalo(2)
      ReadUInt (fphalo, &halo.numSubStruct, swap);    // numSubStruct(3)
      ReadFloat(fphalo, &halo.Mvir,         swap);    // Mvir(4)
      ReadUInt (fphalo, &halo.npart,        swap);    // npart(5)
      ReadFloat(fphalo, &halo.Xc,           swap);    // Xc(6)
      ReadFloat(fphalo, &halo.Yc,           swap);    // Yc(7)
      ReadFloat(fphalo, &halo.Zc,           swap);    // Zc(8)
      ReadFloat(fphalo, &halo.VXc,          swap);    // VXc(9)
      ReadFloat(fphalo, &halo.VYc,          swap);    // VYc(10)
      ReadFloat(fphalo, &halo.VZc,          swap);    // VZc(11)
      ReadFloat(fphalo, &halo.Rvir,         swap);    // Rvir(12)
      ReadFloat(fphalo, &halo.Rmax,         swap);    // Rmax(13)
      ReadFloat(fphalo, &halo.r2,           swap);    // r2(14)
      ReadFloat(fphalo, &halo.mbp_offset,   swap);    // mbp_offset(15)
      ReadFloat(fphalo, &halo.com_offset,   swap);    // com_offset(16)
      ReadFloat(fphalo, &halo.Vmax,         swap);    // Vmax(17)
      ReadFloat(fphalo, &halo.v_esc,        swap);    // v_esc(18)
      ReadFloat(fphalo, &halo.sigV,         swap);    // sigV(19)
      ReadFloat(fphalo, &halo.lambda,       swap);    // lambda(20)
      ReadFloat(fphalo, &halo.lambdaE,      swap);    // lambdaE(21)
      ReadFloat(fphalo, &halo.Lx,           swap);    // Lx(22)
      ReadFloat(fphalo, &halo.Ly,           swap);    // Ly(23)
      ReadFloat(fphalo, &halo.Lz,           swap);    // Lz(24)
      ReadFloat(fphalo, &halo.b,            swap);    // b(25)
      ReadFloat(fphalo, &halo.c,            swap);    // c(26)
      ReadFloat(fphalo, &halo.Eax,          swap);    // Eax(27)
      ReadFloat(fphalo, &halo.Eay,          swap);    // Eay(28)
      ReadFloat(fphalo, &halo.Eaz,          swap);    // Eaz(29)
      ReadFloat(fphalo, &halo.Ebx,          swap);    // Ebx(30)
      ReadFloat(fphalo, &halo.Eby,          swap);    // Eby(31)
      ReadFloat(fphalo, &halo.Ebz,          swap);    // Ebz(32)
      ReadFloat(fphalo, &halo.Ecx,          swap);    // Ecx(33)
      ReadFloat(fphalo, &halo.Ecy,          swap);    // Ecy(34)
      ReadFloat(fphalo, &halo.Ecz,          swap);    // Ecz(35)
      ReadFloat(fphalo, &halo.ovdens,       swap);    // ovdens(36)
      ReadUInt (fphalo, &halo.nbins,        swap);    // nbins(37)
      ReadFloat(fphalo, &halo.fMhires,      swap);    // fMhires(38)
      ReadFloat(fphalo, &halo.Ekin,         swap);    // Ekin(39)
      ReadFloat(fphalo, &halo.Epot,         swap);    // Epot(40)
      ReadFloat(fphalo, &halo.SurfP,        swap);    // SurfP(41)
      ReadFloat(fphalo, &halo.Phi0,         swap);    // Phi0(42)
      ReadFloat(fphalo, &halo.cNFW,         swap);    // cNFW(43)
 
      //=================================================================================
      // write halo properties
      //=================================================================================
      //printf("%ld %f %f %f %f\n",halo.ID,halo.Mvir,halo.Xc,halo.Yc,halo.Zc);

      //
      //printf("ID = %ld\n",halo.ID);
      mhalo.mhalos[counthalo].ID = counthalo;
      mhalo.mhalos[counthalo].refID = pow(10,15)+domain*pow(10,10)+counthalo_local;
      mhalo.mhalos[counthalo].oriID = halo.ID;
      mhalo.mhalos[counthalo].domainID = domain;
      mhalo.mhalos[counthalo].Mvir = halo.Mvir;
      mhalo.mhalos[counthalo].Rvir = halo.Rvir;
      mhalo.mhalos[counthalo].Xc = halo.Xc;
      mhalo.mhalos[counthalo].Yc = halo.Yc;
      mhalo.mhalos[counthalo].Zc = halo.Zc;
      mhalo.mhalos[counthalo].VXc = halo.VXc;
      mhalo.mhalos[counthalo].VYc = halo.VYc;
      mhalo.mhalos[counthalo].VZc = halo.VZc;
      if(halo.hostHalo == 0)
	mhalo.mhalos[counthalo].host_halo = NULLPOINT;
      else
	mhalo.mhalos[counthalo].host_halo = halo.hostHalo;

      ReadULong(fppart, &(npart), swap);
      //memmgr_printdetails();
      //printf("read npart %ld:%ld\n",counthalo,npart);
      mhalo.mhalos[counthalo].npart = npart;
      //memmgr_printdetails();
      //     printf("checking hid:%ld npart:%ld\n",halo.ID,npart);
      if(mhalo.mhalos[counthalo].npart != halo.npart)
	{
	  printf("redshift: %3.3f\n",mhalo.redshift);
	  printf("domain %d\n",domain);
	  printf("haloid:%llu no:%ld local:%ld\n",halo.ID, counthalo, counthalo_local);
	  printf("npart mismatch p:%d, h:%d\n",mhalo.mhalos[counthalo].npart,halo.npart);
	  printf("Xc:%f, Yc:%f, Zc:%f\n",halo.Xc,halo.Yc,halo.Zc);
	  printf("use part from halos\n");
	  mhalo.mhalos[counthalo].npart = halo.npart;
	  flag = 1;
	  //return mhalo;
	  global_error = 1;
	  exit(1);
	  // return mhalo;
	}
      //memmgr_printdetails();
      //printf("counthalo: %ld\n",counthalo);
      sprintf(memmgr_buff,"Particle: Halo Array");
      mhalo.mhalos[counthalo].Particles = memmgr_malloc(mhalo.mhalos[counthalo].npart*sizeof(particlelist_t),memmgr_buff);
      sprintf(memmgr_buff,"Particle: Buffer");
      pid_buff = memmgr_malloc(mhalo.mhalos[counthalo].npart*sizeof(struct particle_buffer),memmgr_buff);
      //memmgr_printdetails();
      //printf("%ld\n",(long unsigned)npart);
      for(ipart=0; ipart<mhalo.mhalos[counthalo].npart; ipart++) 
	{
	  ReadULong(fppart, &(pid_buff[ipart].ID), swap);
	  ReadFloat(fppart, &(pid_buff[ipart].energy), swap);
	  if(flag == 3)
	    {
	      printf("p: %llu \t %f\n",pid_buff[ipart].ID,pid_buff[ipart].energy);
	    }
	}
      /* Do something to sort particles by energy */
      for(ipart=0; ipart<mhalo.mhalos[counthalo].npart; ipart++) 
	{
	  mhalo.mhalos[counthalo].Particles[ipart].ID = pid_buff[ipart].ID;
	}
      memmgr_free(pid_buff,mhalo.mhalos[counthalo].npart*sizeof(struct particle_buffer),memmgr_buff);
      counthalo++;
      counthalo_local++;
      //memmgr_printdetails();
    } // for(numHalos)
  
  return mhalo;
}



make_catalogue_halo_wrapper_t sussexbigrun_load_halo_catalogue_binary_single_chunk(char *folder, float redshift, int snapid, int chunk )
{
  FILE *fphalo,*fppart,*fpprof;
  char halofile[MAXSTRING],partfile[MAXSTRING],proffile[MAXSTRING];
  make_catalogue_halo_wrapper_t chalo;
  uint64_t *maphalo;
  int i;
  hid_t ihalo,jhalo,hostid;
  double dist_sq;
  int indx,indy,indz;
  int chunk_X, chunk_Y, chunk_Z;
  float shift_X, shift_Y, shift_Z;
  float upperbound[3],lowerbound[3];
  chalo.nHalos = 0;
  chalo.redshift = redshift;
  chalo.snapid = snapid;
  chalo.chalos= memmgr_malloc(0,"Halo Array");
  for(i=0;i<param_chunk_mpi;i++)
    {
      sprintf(partfile,"%s/z_%2.3f_178/chunk_%d/%2.3fxv..%04d.z%2.3f.AHF_particles_bin",folder,redshift,chunk,redshift,i,redshift);
      sprintf(halofile,"%s/z_%2.3f_178/chunk_%d/%2.3fxv..%04d.z%2.3f.AHF_halos_bin",folder,redshift,chunk,redshift,i,redshift);
      sprintf(proffile,"%s/z_%2.3f_178/chunk_%d/%2.3fxv..%04d.z%2.3f.AHF_profiles_bin",folder,redshift,chunk,redshift,i,redshift);
      /* printf("reading %s\n",partfile); */
      /* printf("reading %s\n",halofile); */
      /* printf("reading %s\n",proffile); */
      fphalo = fopen(halofile,"rb");
      fppart = fopen(partfile,"rb");  
      fpprof = fopen(proffile,"rb"); 
      if(fphalo && fphalo && fpprof)
	{
	  chalo = sussexbigrun_read_AHF_binary_from_raw(fphalo, fppart, fpprof, chunk, i, chalo);
	  fclose(fphalo);
	  fclose(fppart);
	  fclose(fpprof);
	}
    }
  maphalo = memmgr_malloc(chalo.nHalos*sizeof(uint64_t),"Maphalo");

  /* Check host halo for structure-infant halos */
  /* 1. Sort halos by Mvir  */
  qsort(chalo.chalos,chalo.nHalos,sizeof(make_catalogue_halo_t), compare_make_catalogue_halo_t_by_Mvir_reverse);
  /* 2. Try to map structure-infant to distant < Rvir */
  for(ihalo=0;ihalo<chalo.nHalos;ihalo++)
    {
      maphalo[chalo.chalos[ihalo].refID] = ihalo;
      if(chalo.chalos[ihalo].hostHalo == 0)
	{
	  chalo.chalos[ihalo].hostHalo = NULLPOINT;
	  for(jhalo=ihalo-1;jhalo>=0 && jhalo!=NULLPOINT;jhalo--)
	    {
	      dist_sq = (chalo.chalos[ihalo].Xc-chalo.chalos[jhalo].Xc)*(chalo.chalos[ihalo].Xc-chalo.chalos[jhalo].Xc)
		+(chalo.chalos[ihalo].Yc-chalo.chalos[jhalo].Yc)*(chalo.chalos[ihalo].Yc-chalo.chalos[jhalo].Yc)
		+(chalo.chalos[ihalo].Zc-chalo.chalos[jhalo].Zc)*(chalo.chalos[ihalo].Zc-chalo.chalos[jhalo].Zc);
	      if(sqrt(dist_sq) < chalo.chalos[jhalo].Rvir)
		{
		  chalo.chalos[ihalo].hostHalo = chalo.chalos[jhalo].ID;
		  chalo.chalos[ihalo].UpHalo = chalo.chalos[jhalo].refID;
		  chalo.chalos[ihalo].NextHalo = chalo.chalos[jhalo].FirstDownHalo;
		  chalo.chalos[jhalo].FirstDownHalo = chalo.chalos[ihalo].refID;
		}
	    }
	}
    }
  /* Determine chunk x,y,z */
  for(indz=0;indz<param_chunk_per_dim;indz++)
    {
      for(indy=0;indy<param_chunk_per_dim;indy++)
	{
	  for(indx=0;indx<param_chunk_per_dim;indx++)
	    {
	      if(chunk == indz*pow2(param_chunk_per_dim)+indy*param_chunk_per_dim+indx)
		{
		  chunk_X = indx;
		  chunk_Y = indy;
		  chunk_Z = indz;
		  break;
		}
	    }
	}
    }
  shift_X = chunk_X*param_boxsize/param_chunk_per_dim - 2.*param_buffer_size*param_boxsize/param_domain_per_dim;
  shift_Y = chunk_Y*param_boxsize/param_chunk_per_dim - 2.*param_buffer_size*param_boxsize/param_domain_per_dim;
  shift_Z = chunk_Z*param_boxsize/param_chunk_per_dim - 2.*param_buffer_size*param_boxsize/param_domain_per_dim;
  //printf("chunk: %d,%d,%d\n",chunk_X,chunk_Y,chunk_Z);
  //printf("shift: %f,%f,%f\n",shift_X,shift_Y,shift_Z);

  lowerbound[0] = chunk_X*param_boxsize/param_chunk_per_dim;
  lowerbound[1] = chunk_Y*param_boxsize/param_chunk_per_dim;
  lowerbound[2] = chunk_Z*param_boxsize/param_chunk_per_dim;
  upperbound[0] = (chunk_X+1)*param_boxsize/param_chunk_per_dim;
  upperbound[1] = (chunk_Y+1)*param_boxsize/param_chunk_per_dim;
  upperbound[2] = (chunk_Z+1)*param_boxsize/param_chunk_per_dim;
 
 /* Shift chunk and determine domain */
  for(ihalo=0;ihalo<chalo.nHalos;ihalo++)
    {
      chalo.chalos[ihalo].Xc = chalo.chalos[ihalo].Xc + shift_X;
      chalo.chalos[ihalo].Yc = chalo.chalos[ihalo].Yc + shift_Y;
      chalo.chalos[ihalo].Zc = chalo.chalos[ihalo].Zc + shift_Z;
      if(chalo.chalos[ihalo].hostHalo == NULLPOINT 
	 && chalo.chalos[ihalo].Xc >= lowerbound[0] && chalo.chalos[ihalo].Xc < upperbound[0]
	 && chalo.chalos[ihalo].Yc >= lowerbound[1] && chalo.chalos[ihalo].Yc < upperbound[1]
	 && chalo.chalos[ihalo].Zc >= lowerbound[2] && chalo.chalos[ihalo].Zc < upperbound[2])
	{
	  /* reuse indx,indz,indz - too lazy to define new variables */
	  indx = (int) (chalo.chalos[ihalo].Xc/ (param_boxsize/param_domain_per_dim));
	  indy = (int) (chalo.chalos[ihalo].Yc/ (param_boxsize/param_domain_per_dim));
	  indz = (int) (chalo.chalos[ihalo].Zc/ (param_boxsize/param_domain_per_dim));

	  /* determine domain */
	  chalo.chalos[ihalo].domainid = indz*pow2(param_domain_per_dim)+indy*param_domain_per_dim+indx;
	  //printf("chunk:%d/%d domain:%d = %d,%d,%d\n",chalo.chalos[ihalo].chunkid,chunk,chalo.chalos[ihalo].domainid,indx,indy,indz);
	}
      else if(chalo.chalos[ihalo].hostHalo != NULLPOINT)
	{
	  hostid = maphalo[chalo.chalos[ihalo].UpHalo];
	  chalo.chalos[ihalo].domainid = chalo.chalos[hostid].domainid;
	}
    }
  memmgr_free(maphalo,chalo.nHalos*sizeof(uint64_t),"Maphalo");
  return chalo;
}


/* Use to read from RAW binary AHF */
make_catalogue_halo_wrapper_t sussexbigrun_read_AHF_binary_from_raw(FILE *fphalo, FILE *fppart, FILE *fpprof, int chunk, int partition, make_catalogue_halo_wrapper_t chalo)
{
  uint64_t numHalos,counthalo,counthalo_local,old_nHalos;
  order_uint64_t *maphalo;
  uint64_t numHaloFromPartFile;
  uint64_t numHaloFromProfFile;
  uint32_t numColumns,nbins,ibin;
  uint64_t i,size;
  int32_t  one;
  int      swap=0,flag;
  halo_t   halo;
  ptid_t ipart,npart;

  //ptid_t id;
  size_t old,new;
  char memmgr_buff[memmgr_max_str];
  struct particle_buffer *pid_buff;;
  if(fphalo == NULL || fppart == NULL) 
    {
      printf("Cannot open file to read\n");
      exit(1);
    }

  // figure out swap status Halos file
  fread(&one, sizeof(int32_t), 1, fphalo);
  if(one == 1)
    swap = 0;
  else
    swap = 1;
  
  // figure out swap status Particles file
  fread(&one, sizeof(int32_t), 1, fppart);
  if(one == 1)
     swap = 0;
  else
    swap = 1;

  // figure out swap status Particles file
  fread(&one, sizeof(int32_t), 1, fpprof);
  if(one == 1)
     swap = 0;
  else
    swap = 1;
 
  ReadULong(fpprof, &numHaloFromProfFile,   swap);
  ReadUInt (fpprof, &numColumns, swap);
 
  ReadULong(fppart, &numHaloFromPartFile,   swap);
  ReadUInt (fppart, &numColumns, swap);
  //printf("particlefile: nhalo = %llu\n",numHaloFromPartFile);
  ReadULong(fphalo, &numHalos,   swap);
  ReadUInt (fphalo, &numColumns, swap);
  //printf("halofile: nhalo = %llu\n",numHalos);

  if(numHalos != numHaloFromPartFile ||numHalos != numHaloFromProfFile)
    {
      printf("Number of halos don't match\nExit()\n");
      exit(1);
    }

  sprintf(memmgr_buff,"Halo Array");

  counthalo = chalo.nHalos;
  counthalo_local = 0;
  old_nHalos = chalo.nHalos;
  old = chalo.nHalos*sizeof(make_catalogue_halo_t);
  new = old + numHalos*sizeof(make_catalogue_halo_t);
  chalo.nHalos += numHalos;
  new = chalo.nHalos*sizeof(make_catalogue_halo_t);
  
  chalo.chalos = memmgr_realloc(chalo.chalos,new,old,memmgr_buff);
  flag = 0;

  /* use maphalo to map local IDs to Unique IDs */
  maphalo = memmgr_malloc(numHalos*sizeof(order_uint64_t),"Maphalo");
  for(i=0; i<numHalos; i++) 
    {
      /* Read halo properties from AHF_halos */
      ReadULong(fphalo, &(chalo.chalos[counthalo].ID),           swap);    // ID(1)

      /* maphalo */
      maphalo[i].ref = chalo.chalos[counthalo].ID;
      maphalo[i].id = counthalo;

      /* Make unique ID */
      chalo.chalos[counthalo].ID = chalo.snapid*pow(10,15)+chunk*pow(10,10)+partition*pow(10,7)+i+1;
    
      ReadULong(fphalo, &(chalo.chalos[counthalo].hostHalo),     swap);    // hostHalo(2)
      /* point to (long long unsigned)-1 if hosthalo = 0 */
      if(chalo.chalos[counthalo].hostHalo == 0)
	chalo.chalos[counthalo].hostHalo = NULLPOINT;

      ReadUInt (fphalo, &(chalo.chalos[counthalo].numSubStruct), swap);    // numSubStruct(3)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Mvir),         swap);    // Mvir(4)
      ReadUInt (fphalo, &(chalo.chalos[counthalo].npart),        swap);    // npart(5)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Xc),           swap);    // Xc(6)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Yc),           swap);    // Yc(7)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Zc),           swap);    // Zc(8)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].VXc),          swap);    // VXc(9)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].VYc),          swap);    // VYc(10)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].VZc),          swap);    // VZc(11)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Rvir),         swap);    // Rvir(12)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Rmax),         swap);    // Rmax(13)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].r2),           swap);    // r2(14)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].mbp_offset),   swap);    // mbp_offset(15)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].com_offset),   swap);    // com_offset(16)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Vmax),         swap);    // Vmax(17)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].v_esc),        swap);    // v_esc(18)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].sigV),         swap);    // sigV(19)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].lambda),       swap);    // lambda(20)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].lambdaE),      swap);    // lambdaE(21)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Lx),           swap);    // Lx(22)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ly),           swap);    // Ly(23)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Lz),           swap);    // Lz(24)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].b),            swap);    // b(25)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].c),            swap);    // c(26)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Eax),          swap);    // Eax(27)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Eay),          swap);    // Eay(28)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Eaz),          swap);    // Eaz(29)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ebx),          swap);    // Ebx(30)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Eby),          swap);    // Eby(31)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ebz),          swap);    // Ebz(32)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ecx),          swap);    // Ecx(33)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ecy),          swap);    // Ecy(34)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ecz),          swap);    // Ecz(35)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].ovdens),       swap);    // ovdens(36)
      ReadUInt (fphalo, &(chalo.chalos[counthalo].nbins),        swap);    // nbins(37)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].fMhires),      swap);    // fMhires(38)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Ekin),         swap);    // Ekin(39)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Epot),         swap);    // Epot(40)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].SurfP),        swap);    // SurfP(41)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].Phi0),         swap);    // Phi0(42)
      ReadFloat(fphalo, &(chalo.chalos[counthalo].cNFW),         swap);    // cNFW(43)

      /* Specify other quantities */
      chalo.chalos[counthalo].refID = counthalo;

      chalo.chalos[counthalo].domainid = -1;
      chalo.chalos[counthalo].chunkid = chunk;

 
      /* Set structure tree to default (no relationship) */
      chalo.chalos[counthalo].UpHalo = -1;
      chalo.chalos[counthalo].FirstDownHalo = -1;
      chalo.chalos[counthalo].NextHalo = -1;

      /* Read nparts from AHF_particles */
      ReadULong(fppart, &(npart), swap);

      if(chalo.chalos[counthalo].npart != npart)
	{
	  printf("redshift: %3.3f\n",chalo.redshift);
	  printf("domain %d\n",chunk);
	  printf("haloid:%llu no:%ld local:%ld\n",chalo.chalos[counthalo].ID, counthalo, counthalo_local);
	  printf("npart mismatch p:%d, h:%d\n",chalo.chalos[counthalo].npart,halo.npart);
	  printf("Xc:%f, Yc:%f, Zc:%f\n",chalo.chalos[counthalo].Xc,chalo.chalos[counthalo].Yc,chalo.chalos[counthalo].Zc);
	  flag = 1;
	  global_error = 1;
	  exit(1);
	}

      sprintf(memmgr_buff,"Particle: Halo Array");
      chalo.chalos[counthalo].Particles = memmgr_malloc(chalo.chalos[counthalo].npart*sizeof(particlelist_t),memmgr_buff);
      sprintf(memmgr_buff,"Particle: Buffer");
      pid_buff = memmgr_malloc(chalo.chalos[counthalo].npart*sizeof(struct particle_buffer),memmgr_buff);
      
       /* Read Particle one by one (poor computer,need to read this over and over again) */
      for(ipart=0; ipart<chalo.chalos[counthalo].npart; ipart++) 
	{
	  ReadULong(fppart, &(pid_buff[ipart].ID), swap);
	  ReadFloat(fppart, &(pid_buff[ipart].energy), swap);
	  if(flag == 3)
	    {
	      printf("p: %llu \t %f\n",pid_buff[ipart].ID,pid_buff[ipart].energy);
	    }
	}
      /* Do something to sort particles by energy */
      for(ipart=0; ipart<chalo.chalos[counthalo].npart; ipart++) 
	{
	  chalo.chalos[counthalo].Particles[ipart].ID = pid_buff[ipart].ID;
	}
      memmgr_free(pid_buff,chalo.chalos[counthalo].npart*sizeof(struct particle_buffer),memmgr_buff);
      ReadUInt(fpprof, &nbins, swap);
      if(chalo.chalos[counthalo].nbins != nbins)
	{
	  printf("nbins do not match.\nExit()\n");
	  exit(1);
	}
      AHF_alloc_profiles(chalo.chalos[counthalo].nbins, &(chalo.chalos[counthalo].Profile));
      for(ibin=0; ibin<nbins; ibin++) {
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.r[ibin]),      swap);
	ReadUInt (fpprof, &(chalo.chalos[counthalo].Profile.npart[ibin]),  swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.M_in_r[ibin]), swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.ovdens[ibin]), swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.dens[ibin]),   swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.vcirc[ibin]),  swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.vesc[ibin]),   swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.sigv[ibin]),   swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Lx[ibin]),     swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ly[ibin]),     swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Lz[ibin]),     swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.b[ibin]),      swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.c[ibin]),      swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Eax[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Eay[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Eaz[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ebx[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Eby[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ebz[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ecx[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ecy[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ecz[ibin]),    swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Ekin[ibin]),   swap);
	ReadFloat(fpprof, &(chalo.chalos[counthalo].Profile.Epot[ibin]),   swap);
      }
      
      counthalo++;
      counthalo_local++;
    } // for(numHalos)

  /* Relabel  HostID */
  if(numHalos > 0)
    chalo = sussexbigrun_find_hostHalo(chalo,maphalo,numHalos);
  memmgr_free(maphalo,numHalos*sizeof(order_uint64_t),"Maphalo");
  return chalo;
}

/* This function works only for domain_per_dim%chunk_per_dim = 0 */
/* Use 1 MPI per chunk */
make_catalogue_halo_wrapper_t sussexbigrun_output_cubep3m(make_catalogue_halo_wrapper_t chalo, int chunk)
{
  uint64_t ihalo,count_export[mpi_nodes];
  uint64_t *export_halo[mpi_nodes];
  int domain_to_chunk[pow3(param_domain_per_dim)];
  int i,j,k,inode,jnode,target_chunk;
  int ratio = param_domain_per_dim/param_chunk_per_dim;
  uint64_t send_nhalos,rev_nhalos;
  MPI_Request request[mpi_nodes];
  MPI_Status status[mpi_nodes];
  for(k=0;k<param_domain_per_dim;k++)
    {
      for(j=0;j<param_domain_per_dim;j++)
	{
	  for(i=0;i<param_domain_per_dim;i++)
	    {
	      domain_to_chunk[k*pow2(param_domain_per_dim)+j*param_domain_per_dim+i] =
		k/ratio*pow2(param_chunk_per_dim) + j/ratio*param_chunk_per_dim + i/ratio;
	      //printf("domain:%d => %d\n",k*pow2(param_domain_per_dim)+j*param_domain_per_dim+i,domain_to_chunk[k*pow2(param_domain_per_dim)+j*param_domain_per_dim+i]);
	    }
	}
    }
  for(inode=0;inode<mpi_nodes;inode++)
    {
      count_export[inode] = 0;
      export_halo[inode] = malloc(0);
    }
  for(ihalo=0;ihalo<chalo.nHalos;ihalo++)
    {
      if(chalo.chalos[ihalo].domainid > -1)
	{
	  if(domain_to_chunk[chalo.chalos[ihalo].domainid] != chunk)
	    {
	      //printf("domain: %d is not in chunk %d:%d\n",chalo.chalos[ihalo].domainid,chunk,domain_to_chunk[chalo.chalos[ihalo].domainid]);
	      target_chunk = domain_to_chunk[chalo.chalos[ihalo].domainid];
	      count_export[target_chunk]++;
	      export_halo[target_chunk] = realloc(export_halo[target_chunk], count_export[target_chunk]*sizeof(uint64_t));
	      export_halo[target_chunk][count_export[target_chunk]-1] = ihalo;
	    }
	}
    }

  /* Exchange halos */
  /* Transfer from inode -> jnode */
  rev_nhalos = 0;
  MPI_Barrier(MPI_COMM_WORLD);

  for(inode=0;inode<mpi_nodes;inode++)
    {
      for(jnode=0;jnode<mpi_nodes;jnode++)
  	{
  	  if(mpi_rank == inode)
  	    {
  	      send_nhalos = count_export[jnode];
  	      rev_nhalos = send_nhalos;
  	      MPI_Send(&send_nhalos, 1, MPI_UNSIGNED_LONG_LONG, jnode, mpi_nodes*inode+jnode, MPI_COMM_WORLD);
  	    }
  	  if(mpi_rank == jnode)
  	    {
  	      MPI_Recv(&rev_nhalos, 1, MPI_UNSIGNED_LONG_LONG, inode, mpi_nodes*inode+jnode, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  	    }
  	  MPI_Barrier(MPI_COMM_WORLD);
  	  if(mpi_rank == jnode)
  	    {
  	      chalo.nHalos += rev_nhalos;
  	      chalo.chalos = memmgr_realloc(chalo.chalos,chalo.nHalos*sizeof(make_catalogue_halo_t),(chalo.nHalos-rev_nhalos)*sizeof(make_catalogue_halo_t),"Halo Array");
  	    }
  	  MPI_Barrier(MPI_COMM_WORLD);
	  for(ihalo=0;ihalo<rev_nhalos;ihalo++)
	    {
	      if(mpi_rank == inode)
		{
		  MPI_Send(&(chalo.chalos[export_halo[jnode][ihalo]]), sizeof(make_catalogue_halo_t), MPI_BYTES, jnode, (mpi_nodes*inode+jnode)*rev_nhalos+ihalo, MPI_COMM_WORLD);
		}
	      else if(mpi_rank == jnode)
		{
		  MPI_Recv(&(chalo.chalos[chalo.nHalos-rev_nhalos+ihalo]), sizeof(make_catalogue_halo_t), MPI_BYTES, inode, (mpi_nodes*inode+jnode)*rev_nhalos +ihalo, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  if(mpi_rank == inode || mpi_rank == jnode)
	    rev_nhalos = 0;
  	}
    }
  for(inode=0;inode<mpi_nodes;inode++)
    {
      free(export_halo[inode]);
    }
  return chalo;
}

/* This function will map hostID to the ID we are using */
/* The maphalo needs to be unsorted. I don't want to sort maphalo by id. -Boyd*/
make_catalogue_halo_wrapper_t sussexbigrun_find_hostHalo(make_catalogue_halo_wrapper_t chalo, order_uint64_t *maphalo_unsorted, uint64_t numHalos)
{
  order_uint64_t *maphalo_sorted;
  uint64_t i,hostid_unique_el,startid,stopid;
  startid = maphalo_unsorted[0].id;
  stopid = maphalo_unsorted[numHalos-1].id;
  qsort(maphalo_unsorted, numHalos, sizeof(order_uint64_t), compare_order_uint64_t_by_ref);
  maphalo_sorted = maphalo_unsorted; 	/* This is just to be easy to remember. */
  for(i=startid;i<=stopid;i++)
    { 
 
      if(chalo.chalos[i].hostHalo != NULLPOINT)
	{
	  hostid_unique_el = search_order_unint64_t_for_ref(chalo.chalos[i].hostHalo, numHalos, maphalo_sorted);
	  if(hostid_unique_el != NULLPOINT)
	    {
	      chalo.chalos[i].hostHalo = chalo.chalos[maphalo_sorted[hostid_unique_el].id].ID;
	      chalo.chalos[i].UpHalo = maphalo_sorted[hostid_unique_el].id;
	      chalo.chalos[i].NextHalo = chalo.chalos[maphalo_sorted[hostid_unique_el].id].FirstDownHalo;
	      chalo.chalos[maphalo_sorted[hostid_unique_el].id].FirstDownHalo = chalo.chalos[i].refID;
	    }
	  else
	    { 
	      /* Set hosthalo = 0 for subhalos of mainhalos which are in AHF buffer */
	      chalo.chalos[i].hostHalo = 0;
	    }
	}
    }
  return chalo;
}

void AHF_alloc_profiles( uint32_t nbins, halo_profile_t *prof)
{
  prof->r       = (float *)       calloc(nbins, sizeof(float));
  prof->npart   = (uint32_t*)     calloc(nbins, sizeof(uint32_t));
  prof->M_in_r  = (float *)       calloc(nbins, sizeof(float));
  prof->ovdens  = (float *)       calloc(nbins, sizeof(float));
  prof->dens    = (float *)       calloc(nbins, sizeof(float));
  prof->vcirc   = (float *)       calloc(nbins, sizeof(float));
  prof->vesc    = (float *)       calloc(nbins, sizeof(float));
  prof->sigv    = (float *)       calloc(nbins, sizeof(float));
  prof->Lx      = (float *)       calloc(nbins, sizeof(float));
  prof->Ly      = (float *)       calloc(nbins, sizeof(float));
  prof->Lz      = (float *)       calloc(nbins, sizeof(float));
  prof->b       = (float *)       calloc(nbins, sizeof(float));
  prof->c       = (float *)       calloc(nbins, sizeof(float));
  prof->Eax     = (float *)       calloc(nbins, sizeof(float));
  prof->Eay     = (float *)       calloc(nbins, sizeof(float));
  prof->Eaz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ebx     = (float *)       calloc(nbins, sizeof(float));
  prof->Eby     = (float *)       calloc(nbins, sizeof(float));
  prof->Ebz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecx     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecy     = (float *)       calloc(nbins, sizeof(float));
  prof->Ecz     = (float *)       calloc(nbins, sizeof(float));
  prof->Ekin    = (float *)       calloc(nbins, sizeof(float));
  prof->Epot    = (float *)       calloc(nbins, sizeof(float));
}

void AHF_free_profiles(halo_profile_t *prof)
{
  free(prof->r);
  free(prof->npart);
  free(prof->M_in_r);
  free(prof->ovdens);
  free(prof->dens);
  free(prof->vcirc);
  free(prof->vesc);
  free(prof->sigv);
  free(prof->Lx);
  free(prof->Ly);
  free(prof->Lz);
  free(prof->b);
  free(prof->c);
  free(prof->Eax);
  free(prof->Eay);
  free(prof->Eaz);
  free(prof->Ebx);
  free(prof->Eby);
  free(prof->Ebz);
  free(prof->Ecx);
  free(prof->Ecy);
  free(prof->Ecz);
  free(prof->Ekin);
  free(prof->Epot);

}

void free_make_catalogue_halo_wrapper(make_catalogue_halo_wrapper_t *ptr)
{
  hid_t i,j;
  char buff[memmgr_max_str];
  for(j=0;j<1;j++)
    {
      sprintf(buff,"Particle: Halo Array");
      for(i=0;i<ptr[j].nHalos;i++)
	{
	  memmgr_free(ptr[j].chalos[i].Particles,ptr[j].chalos[i].npart*sizeof(particlelist_t),buff);
	  AHF_free_profiles(&(ptr[j].chalos[i].Profile));
	}
      sprintf(buff,"Halo Array");
      memmgr_free(ptr[j].chalos,ptr[j].nHalos*sizeof(make_catalogue_halo_t),buff);
    }
  sprintf(buff,"Halo wrapper");
  memmgr_free(ptr,sizeof(make_catalogue_halo_t),buff);
}

/* Private functions */
int compare_make_catalogue_halo_t_by_Mvir_reverse(const void *v1, const void *v2)
{
    const make_catalogue_halo_t *u1 = v1;
    const make_catalogue_halo_t *u2 = v2;
    int ret;
    if(u1->Mvir < u2->Mvir)
      ret = 1;
    else if(u1->Mvir > u2->Mvir)
      ret = -1;
    else if(u1->Mvir == u2->Mvir)
      ret = 0;
    return ret;
}

