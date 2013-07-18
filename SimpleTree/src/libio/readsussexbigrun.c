#include "readsussexbigrun.h"


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
  printf("halo in domain %d :%ld \n",domain,numHalos);
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
      mhalo.mhalos[counthalo].ID = halo.ID;
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
	  printf("domain %d\n",domain);
	  printf("haloid:%llu no:%ld local:%ld\n",halo.ID, counthalo, counthalo_local);
	  printf("npart mismatch p:%d, h:%d\n",mhalo.mhalos[counthalo].npart,halo.npart);
	  printf("Xc:%f, Yc:%f, Zc:%f\n",halo.Xc,halo.Yc,halo.Zc);
	  printf("use part from halos\n");
	  mhalo.mhalos[counthalo].npart = halo.npart;
	  flag = 1;
	  return mhalo;
	  //exit(1);
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
