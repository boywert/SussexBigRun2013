#include "output_mt.h"

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

void internalaux_outputs(m_halo_wrapper_t* haloB, char* outputfolder, int domainid)
{
  hid_t ihalo,whalo;
  FILE *fp;
  char filename[1024], foldername[1024];
  char command[1024];
  int l;
  int ih;
  
  sprintf(foldername,"%s/%3.3f",outputfolder,haloB->redshift);
  sprintf(command,"mkdir -p %s", foldername);
  system(command);
  sprintf(filename,"%s/%3.3f/mtaux_%d.dat",outputfolder,haloB->redshift,domainid);
  sprintf(command,"rm -f %s",filename);
  system(command);
 
  fp = fopen(filename, "wb+");
  if(fp != NULL)
    {
      fwrite(&(haloB->nHalos),sizeof(hid_t),1,fp);
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  printf("writing nprog = %d\n",haloB->mhalos[ihalo].nprogs);
	  //fwrite(&(haloB->mhalos[ihalo].nprogs),sizeof(uint32_t),1,fp);
	}
      for(ihalo=0; ihalo < haloB->nHalos; ihalo++)
	{
	  for(whalo=0;whalo<haloB->mhalos[ihalo].nprogs;whalo++)
	    {
	      printf("writing proglist[%llu] = %llu\n",whalo,haloB->mhalos[ihalo].proglist[whalo]);
	      //fwrite(&(haloB->mhalos[ihalo].proglist[whalo]),sizeof(hid_t),1,fp);
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
