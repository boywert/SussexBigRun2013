#include <stdio.h>
#include <stdlib.h>

void main()
{
  char folder[1024];
  char massfname[1024],structfname[1024];
  double pos[3],vel[3];
  FILE *fp;
  int ihalo,jhalo;
  int skip,nhalos,ipart,nparts;
  float mass;
  long long id;
  sprintf(folder,"/scratch/00916/tg459470/clues/4096/reduced/output_00138/fofres/halos");
  sprintf(structfname,"%s/halos_strct_%05d",folder,filter);
  fp = fopen(structfname,"rb");
  
  /* read nHalos */
  fread (&skip,1,sizeof(int),fp);
  fread (&nhalos,1,sizeof(int),fp);
  fread (&skip,1,sizeof(int),fp);

  for(ihalo=0;ihalo<nhalos;ihalo++)
    {
      fread (&skip,1,sizeof(int),fp);
      fread (&nparts,1,sizeof(int),fp);
      fread (&skip,1,sizeof(int),fp);
      
      fread (skip,1,sizeof(int),fp);
      for(ipart=0;ipart<nparts;ipart++)
	{
	  fread (&(pos[0]), 3, sizeof(float),fp);
	}
      fread (skip,1,sizeof(int),fp);

      fread (skip,1,sizeof(int),fp);
      for(ipart=0;ipart<nparts;ipart++)
	{
	  fread (&(vel[0]), 3, sizeof(float),fp);
	}
      fread (skip,1,sizeof(int),fp);

      fread (skip,1,sizeof(int),fp);
      for(ipart=0;ipart<nparts;ipart++)
	{
	  fread (&id, 1, sizeof(long long),fp);
	}
      fread (skip,1,sizeof(int),fp);
    }
}
