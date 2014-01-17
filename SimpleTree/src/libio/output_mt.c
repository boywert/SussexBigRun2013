#include "output_mt.h"

void internalaux_read(m_halo_wrapper_t* haloB, char* outputfolder, int domainid);
void create_subfind_substruct(m_halo_wrapper_t* haloB);

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

/*  */
/* Need haloB->mhalos[x].ID = x */
void create_subfind_substruct(m_halo_wrapper_t* haloB)
{
  hid_t ihalo;
  hid_t upperhost,hosthalo,ref,cursub,curid;
  *haloB = maphalo_to_host_mt(*haloB);
  printf("Generate Subfind structure\n");

  for(ihalo=0;ihalo<haloB->nHalos; ihalo++)
    {
      hosthalo = haloB->mhalos[ihalo].host_halo;
      upperhost = haloB->mhalos[ihalo].host_halo;
      ref = upperhost;
      printf("0 %llu host = %llu\n",ihalo,hosthalo);
      while(upperhost < NULLPOINT)
	{
	  ref = upperhost;
	  upperhost = haloB->mhalos[upperhost].host_halo;
	}
      if(hosthalo != ref)
	{
	  hosthalo = ref;
	}
      haloB->mhalos[ihalo].UpHalo = hosthalo;
      printf("1 %llu hosthalo = %llu\n",ihalo,hosthalo);
      if(haloB->mhalos[ihalo].UpHalo < NULLPOINT)
	{
	  cursub = haloB->mhalos[ihalo].UpHalo;
	  while(cursub < NULLPOINT)
	    {
	      curid = haloB->mhalos[cursub].ID;
	      cursub = haloB->mhalos[cursub].NextHalo;
	    }
	  haloB->mhalos[curid].NextHalo = ihalo;
	  printf("   nexthalo = %llu\n",ihalo);
	}
      printf("2 %llu hosthalo = %llu\n",ihalo,hosthalo);
    }
  printf("finish linking struct\n");
  /* rename everything to globalRefID */
  for(ihalo=0;ihalo< haloB->nHalos; ihalo++)
    {
      if(haloB->mhalos[ihalo].NextHalo < NULLPOINT)
	haloB->mhalos[ihalo].NextHalo = haloB->mhalos[haloB->mhalos[ihalo].NextHalo].globalRefID;
      if(haloB->mhalos[ihalo].UpHalo < NULLPOINT)
	haloB->mhalos[ihalo].UpHalo = haloB->mhalos[haloB->mhalos[ihalo].UpHalo].globalRefID;	
    }
  printf("finish relabel struct\n");
} 

void internalaux_outputs(m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo,whalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  int ih;
  float Pos[3],Vel[3],VelDisp,Vmax,Spin[3];
  long long MostBoundID;


 
  create_subfind_substruct(haloB);

  sprintf(foldername,"%s/%3.3f",outputfolder,haloB->redshift);
  sprintf(command,"mkdir -p %s", foldername);
  system(command);
  sprintf(filename,"%s/%3.3f/mtaux_%d.dat",outputfolder,haloB->redshift,domainid);
  sprintf(command,"rm -f %s",filename);
  system(command);

  fp = fopen(filename, "wb+");
  if(fp != NULL)
    {
      /* write total halo number */
      fwrite(&(haloB->nHalos),sizeof(hid_t),1,fp);
      /* write nprogs */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fwrite(&(haloB->mhalos[ihalo].nprogs),sizeof(uint32_t),1,fp);
	}
      /* write proglist */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  for(whalo=0; whalo < haloB->mhalos[ihalo].nprogs; whalo++)
	    {
	      fwrite(&(haloB->mhalos[ihalo].proglist[whalo]),sizeof(hid_t),1,fp);
	    }
	}

      /* write globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fwrite(&(haloB->mhalos[ihalo].globalRefID),sizeof(hid_t),1,fp);
	}
      /* write FirstFOF globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fwrite(&(haloB->mhalos[ihalo].UpHalo),sizeof(hid_t),1,fp);
	}  
      /* write NextFOF globalRefID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fwrite(&(haloB->mhalos[ihalo].NextHalo),sizeof(hid_t),1,fp);
	}  
      /* write M200 */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fwrite(&(haloB->mhalos[ihalo].Mvir),sizeof(float),1,fp);
	}
      /* write Pos[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  Pos[0] = haloB->mhalos[ihalo].Xc*kpc2Gadget;
	  Pos[1] = haloB->mhalos[ihalo].Yc*kpc2Gadget;
	  Pos[2] = haloB->mhalos[ihalo].Zc*kpc2Gadget;
	  fwrite(Pos,sizeof(float),3,fp);
	}
      /* write Vel[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  Vel[0] = haloB->mhalos[ihalo].VXc;
	  Vel[1] = haloB->mhalos[ihalo].VYc;
	  Vel[2] = haloB->mhalos[ihalo].VZc;
	  fwrite(Vel,sizeof(float),3,fp);
	}
      /* write VelDisp */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  VelDisp = haloB->mhalos[ihalo].sigV;
	  fwrite(&(VelDisp),sizeof(float),1,fp);
	}     
      /* write Vmax */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  Vmax = haloB->mhalos[ihalo].Vmax;
	  fwrite(&(Vmax),sizeof(float),1,fp);
	}    
      /* write Spin[3] */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  Spin[0] = haloB->mhalos[ihalo].SpinX;
	  Spin[1] = haloB->mhalos[ihalo].SpinY;
	  Spin[2] = haloB->mhalos[ihalo].SpinZ;
	  fwrite(Spin,sizeof(float),3,fp);
	}    
      /* write MostBoundID */
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  MostBoundID = haloB->mhalos[ihalo].Particles[0].ID;
	  fwrite(&(MostBoundID),sizeof(float),3,fp);
	}    
      fclose(fp);
    }
  else
    {
      printf("Cannot open file %s\nExiting...\n",filename);
      exit(1);
    }
}


void internalaux_read(m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo,whalo,nhalos;
  FILE *fp;
  char filename[1024];
  int l;
  int ih;
  
  sprintf(filename,"%s/%3.3f/mtaux_%d.dat",outputfolder,haloB->redshift,domainid);

  //printf("start output internalaux\n");
  fp = fopen(filename, "rb");
  if(fp != NULL)
    {
      fread(&(nhalos),sizeof(hid_t),1,fp);
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  fread(&(haloB->mhalos[ihalo].nprogs),sizeof(uint32_t),1,fp);
	}
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  for(whalo=0; whalo < haloB->mhalos[ihalo].nprogs; whalo++)
	    {
	      fread(&(haloB->mhalos[ihalo].proglist[whalo]),sizeof(hid_t),1,fp);
	    }
	}
      fclose(fp);
    }
  else
    {
      printf("Cannot open file %s\nExiting...\n",filename);
      exit(1);
    }
}
