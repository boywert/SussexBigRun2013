#include "readsussexbigrun.h"

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
/* incomplete */
void sussexbigrun_makestruct_tree(m_halo_wrapper_t mhalo)
{
  ptid_t ip;
  hid_t ih;
  ip = 1;
  ih = 1;
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
  FILE *fphalo,*fppart;
  char halofile[MAXSTRING],partfile[MAXSTRING];
  make_catalogue_halo_wrapper_t chalo;
  int i;
  hid_t ihalo;
  chalo.nHalos = 0;
  chalo.redshift = redshift;
  chalo.snapid = snapid;
  chalo.chalos= memmgr_malloc(0,"Halo Array");
  for(i=0;i<chunk_mpi;i++)
    {
      sprintf(partfile,"%s/z_%2.3f_178/chunk_%d/%2.3fxv..%04d.z%2.3f.AHF_particles_bin",folder,redshift,chunk,redshift,i,redshift);
      sprintf(halofile,"%s/z_%2.3f_178/chunk_%d/%2.3fxv..%04d.z%2.3f.AHF_halos_bin",folder,redshift,chunk,redshift,i,redshift);
      fphalo = fopen(halofile,"rb");
      fppart = fopen(partfile,"rb");  
      if(fphalo && fphalo)
	{
	  chalo = sussexbigrun_read_AHF_binary_from_raw(fphalo, fppart, chunk, i, chalo);
	  fclose(fphalo);
	  fclose(fppart);
	}
    }
  return chalo;
}


/* Use to read from RAW binary AHF */
make_catalogue_halo_wrapper_t sussexbigrun_read_AHF_binary_from_raw(FILE *fphalo, FILE *fppart, int chunk, int partition, make_catalogue_halo_wrapper_t chalo)
{
  uint64_t numHalos,counthalo,counthalo_local,old_nHalos;
  order_uint64_t *maphalo;
  uint64_t numHaloFromPartFile;
  uint32_t numColumns;
  uint64_t i,size;
  int32_t  one;
  int      swap=0,flag;
  halo_t   halo;
  ptid_t ipart,npart;
  //ptid_t id;
  size_t old,new;
  char memmgr_buff[memmgr_max_str];
  float maxx,maxy,maxz;
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
  
  ReadULong(fppart, &numHaloFromPartFile,   swap);
  ReadUInt (fppart, &numColumns, swap);
  //printf("particlefile: nhalo = %llu\n",numHaloFromPartFile);
  ReadULong(fphalo, &numHalos,   swap);
  ReadUInt (fphalo, &numColumns, swap);
  //printf("halofile: nhalo = %llu\n",numHalos);
 
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
  maxx = 0.;
  maxy = 0.;
  maxz = 0.;
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
      //chalo.chalos[counthalo].ID = chalo.snapid*pow(10,15)+chunk*pow(10,10)+partition*pow(10,7)+counthalo+1;
    
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
      chalo.chalos[counthalo].refID = chalo.snapid*pow(10,15)+chunk*pow(10,10)+partition*pow(10,7)+counthalo+1;
      chalo.chalos[counthalo].domainid = -1;
      chalo.chalos[counthalo].chunkid = chunk;

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
      counthalo++;
      counthalo_local++;
    } // for(numHalos)
  //printf("max = %f, %f, %f\n",maxx,maxy,maxz);

  /* Relabel ID and HostID */
  chalo = sussexbigrun_make_treestruct(chalo,maphalo,numHalos);
  
  memmgr_free(maphalo,numHalos*sizeof(order_uint64_t),"Maphalo");
  return chalo;
}

/* This function will map hostID to the ID we are using */
/* The maphalo needs to be unsorted. I don't want to sort maphalo by id. -Boyd*/
make_catalogue_halo_wrapper_t sussexbigrun_make_treestruct(make_catalogue_halo_wrapper_t chalo, order_uint64_t *maphalo_unsorted, uint64_t numHalos)
{
  order_uint64_t *maphalo_sorted;
  uint64_t i,j,count,hostid_unique_el,startid,stopid;
  double dist_sq;
  printf("start: make tree struct\n");
  startid = maphalo_unsorted[0].id;
  stopid = maphalo_unsorted[numHalos-1].id;

  qsort(maphalo_unsorted, numHalos, sizeof(order_uint64_t), compare_order_uint64_t_by_ref);
  maphalo_sorted = maphalo_unsorted; 	/* This is just to be easy to remember. */
  count = 0;
  for(i=startid;i<=stopid;i++)
    {  
      printf("%llu : %llu\n",chalo.chalos[i].hostHalo);
      if(chalo.chalos[i].hostHalo != NULLPOINT)
	{
	  hostid_unique_el = search_order_unint64_t_for_ref(chalo.chalos[i].hostHalo, numHalos, maphalo_sorted);
	  printf("\t%llu:%llu\n",hostid_unique_el,numHalos);
	  if(hostid_unique_el != NULLPOINT)
	    {
	      chalo.chalos[i].hostHalo = chalo.chalos[maphalo_sorted[hostid_unique_el].id].refID;
	    }
	  else
	    {
	      chalo.chalos[i].hostHalo = NULLPOINT;
	    }
	}
      count++;
    }
  printf("End: making tree struct\n");
  return chalo;
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
	}
      sprintf(buff,"Halo Array");
      memmgr_free(ptr[j].chalos,ptr[j].nHalos*sizeof(make_catalogue_halo_t),buff);
    }
  sprintf(buff,"Halo wrapper");
  memmgr_free(ptr,sizeof(make_catalogue_halo_t),buff);
}
