#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif


struct halo_data
{
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  /* properties of halo */

  int Len;
  float M_Mean200, M_Crit200, M_TopHat;
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  long long MostBoundID;


  /* original position in subfind output */

  int SnapNum, FileNr, SubhaloIndex;
  float SubhalfMass;
#ifdef SAVE_MASS_TAB
  float SubMassTab[6];
#endif
}
 *Halo, *HaloList;

struct halo_aux_data
{
  int TargetIndex;
  int Origin;
  short int FileNr;
  char UsedFlag;
  char HaloFlag;
}
 *HaloAux;


int CountUsed, CountSumUsed;

void walk_it(int i, int flag)
{
  int p;

  HaloAux[i].UsedFlag = 1;
  HaloAux[i].TargetIndex = CountUsed;
  HaloAux[CountUsed].Origin = i;

  if(flag == 1)
    {
      HaloList[CountUsed] = Halo[i];
      HaloAux[i].FileNr = -1; /* to prevent that it is used again in another file */
    }

  CountUsed++;

  if(Halo[i].Descendant >= 0)
    {
      if(HaloAux[Halo[i].Descendant].UsedFlag == 0)
	walk_it(Halo[i].Descendant, flag);
    }

  p = Halo[i].FirstProgenitor;
  while(p >= 0)
    {
      if(HaloAux[p].UsedFlag == 0)
	walk_it(p, flag);

      p = Halo[p].NextProgenitor;
    }

  p = Halo[i].FirstHaloInFOFgroup;
  if(HaloAux[p].HaloFlag == 0)
    {
      HaloAux[p].HaloFlag = 1;
      while(p >= 0)
	{
	  if(HaloAux[p].UsedFlag == 0)
	    walk_it(p, flag);
	  p = Halo[p].NextHaloInFOFgroup;
	}
    }
}





int main(int argc, char **argv)
{
  int num, count;

  if(argc != 2)
    {
      printf("\n  usage: B-HaloTrees <parameterfile>\n");
      printf("  <parameterfile>    see readparmeterfile.\n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);


  Cats = mymalloc(sizeof(struct halo_catalogue) * (LastSnapShotNr + 1));
  FirstHaloInSnap = mymalloc(sizeof(int) * (LastSnapShotNr + 1));

  count_halos();

  printf("TotHalos = %d\n", TotHalos);
  printf("Want to allocate %g GB\n",
	 ((double) TotHalos) * (sizeof(struct halo_data) +
				sizeof(struct halo_aux_data)) / (1024.0 * 1024.0 * 1024.0));

  Halo = mymalloc(TotHalos * sizeof(struct halo_data));
  HaloAux = mymalloc(TotHalos * sizeof(struct halo_aux_data));


  printf("loading halo catalogues and descendant tree files...\n");
  fflush(stdout);

  for(num = LastSnapShotNr, count = 0; num >= FirstSnapShotNr; num -= SnapSkipFac)
    {
      FirstHaloInSnap[num] = count;

      load_subhalo_catalogue(num);

      count += Cats[num].TotNsubhalos;
    }

  printf("done.\n");


  set_progenitor_pointers();

  printf("progenitor pointers done.\n");
  fflush(stdout);

  generate_trees();

  return 0;
}



void generate_trees(void)
{
  int filenr, i, k, maxhalos, treenr;
  int NtreesPerFile, NhalosPerFile, *npertree;
  char buf[500];
  FILE *fd;
  double drand48(void);

  srand48(42);

  for(i = 0; i < TotHalos; i++)
    HaloAux[i].FileNr = NumberOfOutputFiles * drand48();
 

  for(filenr=0; filenr < NumberOfOutputFiles; filenr++)
    {

      NtreesPerFile = 0;
      NhalosPerFile = 0;
      
      CountSumUsed = 0;
      
      for(i = 0; i < TotHalos; i++)
	HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;
      
      maxhalos = 0;
      
      for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
	{
	  if(HaloAux[i].UsedFlag == 0 && HaloAux[i].FileNr == filenr)
	    {
	      CountUsed = 0;
	      
	      walk_it(i, 0);
	      
	      NtreesPerFile += 1;
	      NhalosPerFile += CountUsed;
	      
	      if(CountUsed > maxhalos)
		maxhalos = CountUsed;
	      
	      CountSumUsed += CountUsed;
	    }
	}

      printf("TotHalos=%d   Used=%d   maxhalos=%d\n", TotHalos, CountSumUsed, maxhalos);
      fflush(stdout);
      
      for(i = 0; i < TotHalos; i++)
	HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;
      

      HaloList = mymalloc(maxhalos * sizeof(struct halo_data));
      
      sprintf(buf, "%s/treedata/trees_%03d.%d", OutputDir, LastSnapShotNr, filenr);
      //sprintf(buf, "%s/treedata/trees_sf%d_%03d.%d", OutputDir, SnapSkipFac, LastSnapShotNr, filenr);
      
      printf("starting: %s\n", buf);
      fflush(stdout);
      
      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      
      my_fwrite(&NtreesPerFile, 1, sizeof(int), fd);
      my_fwrite(&NhalosPerFile, 1, sizeof(int), fd);

      printf("NtreesPerFile=%d  NhalosPerFile=%d\n", NtreesPerFile, NhalosPerFile);
      
      fseek(fd, NtreesPerFile * sizeof(int), SEEK_CUR);
      
      npertree = mymalloc(NtreesPerFile * sizeof(int));
      for(i = 0; i < NtreesPerFile; i++)
	npertree[i] = 0;
      
      treenr = 0;
      
      for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
	{
	  if(HaloAux[i].UsedFlag == 0 && HaloAux[i].FileNr == filenr)
	    {
	      CountUsed = 0;
	      
	      walk_it(i, 1);
	      
	      for(k = 0; k < CountUsed; k++)
		{
		  if(HaloList[k].Descendant >= 0)
		    HaloList[k].Descendant = HaloAux[HaloList[k].Descendant].TargetIndex;
		  
		  if(HaloList[k].FirstProgenitor >= 0)
		    HaloList[k].FirstProgenitor = HaloAux[HaloList[k].FirstProgenitor].TargetIndex;
		  
		  if(HaloList[k].NextProgenitor >= 0)
		    HaloList[k].NextProgenitor = HaloAux[HaloList[k].NextProgenitor].TargetIndex;
		  
		  if(HaloList[k].FirstHaloInFOFgroup >= 0)
		    HaloList[k].FirstHaloInFOFgroup = HaloAux[HaloList[k].FirstHaloInFOFgroup].TargetIndex;
		  
		  if(HaloList[k].NextHaloInFOFgroup >= 0)
		    HaloList[k].NextHaloInFOFgroup = HaloAux[HaloList[k].NextHaloInFOFgroup].TargetIndex;
		}
	      
	      my_fwrite(HaloList, CountUsed, sizeof(struct halo_data), fd);
	      
	      npertree[treenr] = CountUsed;
	      treenr++;
	    }
	  
	}
      
      fclose(fd);
      
      if(!(fd = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      fseek(fd, 2 * sizeof(int), SEEK_SET);
      my_fwrite(npertree, NtreesPerFile, sizeof(int), fd);
      fclose(fd);
      
      myfree(npertree);
      
      
      myfree(HaloList);
      
      printf("Saved=%d\n", CountSumUsed);

    }
}



void count_halos(void)
{
  int num, ngroups, nids, nFiles, nsubhalos;
  long long totNids;
  char buf[1000];
  //FILE *fd;

  TotHalos = 0;

  for(num = LastSnapShotNr; num >= FirstSnapShotNr; num -= SnapSkipFac)
    {
      nFiles = 1;

      sprintf(buf, "%s/groups_%03d/fof_subhalo_tab_%03d.%d.hdf5", OutputDir, num, num, 0);
      hid_t   fd, hd, attr;
      herr_t  ret; 
      fd = H5Fopen(buf, H5F_ACC_RDONLY,  H5P_DEFAULT);
      if(fd < 0) {
	printf("can't open file `%s'\n", buf);
	exit(1);
      }
      hd = H5Gopen (fd, "/Header");
      
      attr = H5Aopen(hd, "Ngroups_ThisFile",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &ngroups);
      ret = H5Aclose(attr);
      //my_fread(&ngroups, sizeof(int), 1, fd);
      // printf("ngroup = %d\n",ngroups);

      attr = H5Aopen(hd, "Ngroups_Total",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &Cats[num].TotNgroups);
      ret = H5Aclose(attr);
      //my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
      // printf("TotNgroup = %d\n",cat.TotNgroups);

      attr = H5Aopen(hd, "Nids_ThisFile",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &nids);
      ret = H5Aclose(attr);
      //my_fread(&nids, sizeof(int), 1, fd);
      //printf("nids = %d\n",nids);

      attr = H5Aopen(hd, "Nids_Total",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_LLONG,  &totNids);
      ret = H5Aclose(attr);
      //my_fread(&cat->TotNids, sizeof(long long), 1, fd);
      //printf("TotNids = %d\n",cat->TotNids);

      attr = H5Aopen(hd, "NumFiles",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &nFiles);
      ret = H5Aclose(attr);    
      //my_fread(&nFiles, sizeof(int), 1, fd);

      attr = H5Aopen(hd, "Nsubgroups_ThisFile",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &nsubhalos);
      ret = H5Aclose(attr);      
      //my_fread(&nsubhalos, sizeof(int), 1, fd);
      //printf("nsubhalos = %d\n",nsubhalos);

      attr = H5Aopen(hd, "Nsubgroups_Total",  H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &Cats[num].TotNsubhalos);
      ret = H5Aclose(attr);     
      //my_fread(&cat->TotNsubhalos, sizeof(int), 1, fd);
      //printf("TotNsubhalos = %d\n",cat.TotNsubhalos);

      ret = H5Gclose(hd);
      ret = H5Fclose(fd);

      
      //my_fread(&ngroups, sizeof(int), 1, fd);
      //my_fread(&Cats[num].TotNgroups, sizeof(int), 1, fd);
      //my_fread(&nids, sizeof(int), 1, fd);
      //my_fread(&totNids, sizeof(long long), 1, fd);
      //my_fread(&nFiles, sizeof(int), 1, fd);
      //my_fread(&nsubhalos, sizeof(int), 1, fd);
      //my_fread(&Cats[num].TotNsubhalos, sizeof(int), 1, fd);

      //fclose(fd);
      
      TotHalos += Cats[num].TotNsubhalos;
    }

  printf("total number of halos=%d\n", TotHalos);

}





void load_subhalo_catalogue(int num)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount;
  int groupcount, filenr, ncount;
  int subgr, gr, nh, sc, gr_nh;
  char buf[1000];
  FILE *fd;
  int *nsubPerHalo, *subLen, *descendant_haloindex, *descendant_snapnum, *filenrOfHalo, *subhaloindex;
  float *halo_M_Mean200, *halo_M_Crit200, *halo_M_TopHat;
  float *subpos, *subvel, *subveldisp, *subvmax, *subspin, *subhalfmass;
  MyIDType *subMostBoundID;

#ifdef SAVE_MASS_TAB
  float *submasstab;
#endif

  printf("Catalogue num=%d\n", num);


  nsubPerHalo = mymalloc(sizeof(int) * Cats[num].TotNgroups);
  subLen = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_haloindex = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_snapnum = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  filenrOfHalo = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  subhaloindex = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  
  halo_M_Mean200 = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  halo_M_Crit200 = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  halo_M_TopHat = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  
  subpos = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subvel = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subveldisp = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
  subvmax = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
  subspin = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subMostBoundID = mymalloc(sizeof(MyIDType) * Cats[num].TotNsubhalos);
  subhalfmass = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
#ifdef SAVE_MASS_TAB
  submasstab = mymalloc(6 * sizeof(float) * Cats[num].TotNsubhalos);
#endif


  subcount = 0;
  groupcount = 0;

  nFiles = 1;

  for(filenr = 0; filenr < nFiles; filenr++) {
    sprintf(buf, "%s/groups_%03d/fof_subhalo_tab_%03d.%d.hdf5", OutputDir, num, num, filenr);
    hid_t   fd, hd, attr, sr, gr, id, dset;
    herr_t  ret;
    fd = H5Fopen(buf, H5F_ACC_RDONLY,  H5P_DEFAULT);
    if(fd < 0) {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }
    /* if(!(fd = fopen(buf, "r"))) */
    /* 	{ */
    /* 	  printf("can't open file `%s'\n", buf); */
    /* 	  exit(1); */
    /* 	} */
      
    long long totNids;

    // Read header
    hd = H5Gopen (fd, "/Header");
      
    attr = H5Aopen(hd, "Ngroups_ThisFile",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &ngroups);
    ret = H5Aclose(attr);
    //my_fread(&ngroups, sizeof(int), 1, fd);
    // printf("ngroup = %d\n",ngroups);

    attr = H5Aopen(hd, "Ngroups_Total",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &Cats[num].TotNgroups);
    ret = H5Aclose(attr);
    //my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
    // printf("TotNgroup = %d\n",cat.TotNgroups);

    attr = H5Aopen(hd, "Nids_ThisFile",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &nids);
    ret = H5Aclose(attr);
    //my_fread(&nids, sizeof(int), 1, fd);
    //printf("nids = %d\n",nids);

    attr = H5Aopen(hd, "Nids_Total",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_LLONG, &totNids);
    ret = H5Aclose(attr);
    //my_fread(&cat->TotNids, sizeof(long long), 1, fd);
    //printf("TotNids = %d\n",cat->TotNids);

    attr = H5Aopen(hd, "NumFiles",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &nFiles);
    ret = H5Aclose(attr);    
    //my_fread(&nFiles, sizeof(int), 1, fd);

    attr = H5Aopen(hd, "Nsubgroups_ThisFile",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &nsubhalos);
    ret = H5Aclose(attr);      
    //my_fread(&nsubhalos, sizeof(int), 1, fd);
    //printf("nsubhalos = %d\n",nsubhalos);

    attr = H5Aopen(hd, "Nsubgroups_Total",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &Cats[num].TotNsubhalos);
    ret = H5Aclose(attr); 
    /* my_fread(&ngroups, sizeof(int), 1, fd); */
    /* my_fread(&Cats[num].TotNgroups, sizeof(int), 1, fd); */
    /* my_fread(&nids, sizeof(int), 1, fd); */
    /* my_fread(&totNids, sizeof(long long), 1, fd); */
    /* my_fread(&nFiles, sizeof(int), 1, fd); */
    /* my_fread(&nsubhalos, sizeof(int), 1, fd); */
    /* my_fread(&Cats[num].TotNsubhalos, sizeof(int), 1, fd); */



    /* fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupLen  *\/ */
    /* fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupOffset  *\/ */
    /* fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  GroupMass  *\/ */
    /* fseek(fd, 3 * sizeof(float) * ngroups, SEEK_CUR);	/\* skip  GroupPos *\/ */


    if(ngroups > 0) {
      gr = H5Gopen (fd, "/Group");
      dset = H5Dopen (gr, "Group_M_Mean200");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &halo_M_Mean200[groupcount]);
      ret = H5Dclose(dset);
      //my_fread(&halo_M_Mean200[groupcount], sizeof(float), ngroups, fd);
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Mean200 */

      dset = H5Dopen (gr, "Group_M_Crit200");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &halo_M_Crit200[groupcount]);
      ret = H5Dclose(dset);
      //my_fread(&halo_M_Crit200[groupcount], sizeof(float), ngroups, fd);
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Crit200 */
      dset = H5Dopen (gr, "Group_M_TopHat200");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &halo_M_TopHat[groupcount]);
      ret = H5Dclose(dset);
      //my_fread(&halo_M_TopHat[groupcount], sizeof(float), ngroups, fd);


      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_TopHat200 */

#ifdef FLAG_GROUP_VELDISP
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Mean200 */
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Crit200 */
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_TopHat200 */
#endif
      //fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupContaminationCount */
      //fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupContaminationMass */



      dset = H5Dopen (gr, "GroupNsubs");
      ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nsubPerHalo[groupcount]);
      ret = H5Dclose(dset);	
      //my_fread(&nsubPerHalo[groupcount], sizeof(int), ngroups, fd);

      //fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupFirstsub */
      ret = H5Gclose(gr);
    }
    if(nsubhalos > 0) {
	
      sr = H5Gopen (fd, "/Subhalo");
      dset = H5Dopen (sr, "SubhaloLen");
      ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subLen[subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subLen[subcount], sizeof(int), nsubhalos, fd);


      //fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  SubOffset */
      //fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  SubParenthalo */

      //fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloMass */
      dset = H5Dopen (sr, "SubhaloPos");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subpos[3*subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subpos[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      dset = H5Dopen (sr, "SubhaloVel");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subvel[3*subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subvel[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      //fseek(fd, 3 * sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloCM */

      dset = H5Dopen (sr, "SubhaloSpin");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subspin[subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subspin[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      dset = H5Dopen (sr, "SubhaloVelDisp");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subveldisp[subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subveldisp[subcount], sizeof(float), nsubhalos, fd);


      dset = H5Dopen (sr, "SubhaloVmax");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subvmax[subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subvmax[subcount], sizeof(float), nsubhalos, fd);

      //fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloVmaxRad */

      dset = H5Dopen (sr, "SubhaloHalfmassRad");
      ret = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subhalfmass[subcount]);
      ret = H5Dclose(dset);
      //my_fread(&subhalfmass[subcount], sizeof(float), nsubhalos, fd);

      dset = H5Dopen (sr, "SubhaloIDMostbound");
      ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &subMostBoundID[subcount]);
      ret = H5Dclose(dset);	
      //my_fread(&subMostBoundID[subcount], sizeof(MyIDType), nsubhalos, fd);
	
      //fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  GrNr */

#ifdef SAVE_MASS_TAB
      printf("SAVE_MASS_TAB cannot be use - too lazy - Boyd\n");
      exit(1);
      my_fread(&submasstab[6 * subcount], 6 * sizeof(float), nsubhalos, fd);
#endif
      //fclose(fd);
      ret = H5Gclose(sr);	  
    }
    for(subgr = 0; subgr < nsubhalos; subgr++)
      filenrOfHalo[subcount + subgr] = filenr;
      
    for(subgr = 0; subgr < nsubhalos; subgr++)
      subhaloindex[subcount + subgr] = subgr;
  
    ret = H5Fclose(fd);
    subcount += nsubhalos;
    groupcount += ngroups;
  }
  if(num < LastSnapShotNr)
    {
      sprintf(buf, "%s/treedata/sub_desc_%03d", OutputDir, num);
      //sprintf(buf, "%s/treedata/sub_desc_sf%d_%03d", OutputDir, SnapSkipFac, num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}

      my_fread(&ncount, sizeof(int), 1, fd);
      my_fread(descendant_haloindex, sizeof(int), Cats[num].TotNsubhalos, fd);
      my_fread(descendant_snapnum, sizeof(int), Cats[num].TotNsubhalos, fd);

      fclose(fd);
    }

  nh = FirstHaloInSnap[num];
  sc = 0;

  for(gr = 0; gr < Cats[num].TotNgroups; gr++)
    {
      for(subgr = 0, gr_nh = nh; subgr < nsubPerHalo[gr]; subgr++, sc++, nh++)
	{
	  Halo[nh].FirstHaloInFOFgroup = gr_nh;
	  if(subgr == nsubPerHalo[gr] - 1)
	    Halo[nh].NextHaloInFOFgroup = -1;
	  else
	    Halo[nh].NextHaloInFOFgroup = nh + 1;

	  if(num < LastSnapShotNr)
	    {
	      if(descendant_haloindex[sc] >= 0)
		Halo[nh].Descendant = FirstHaloInSnap[descendant_snapnum[sc]] + descendant_haloindex[sc];
	      else
		Halo[nh].Descendant = -1;
	    }
	  else
	    Halo[nh].Descendant = -1;

	  Halo[nh].FirstProgenitor = -1;
	  Halo[nh].NextProgenitor = -1;

	  /* assign properties */

	  Halo[nh].Len = subLen[sc];

	  if(subgr == 0)
	    {
	      Halo[nh].M_Mean200 = halo_M_Mean200[gr];
	      Halo[nh].M_Crit200 = halo_M_Crit200[gr];
	      Halo[nh].M_TopHat = halo_M_TopHat[gr];
	    }
	  else
	    {
	      Halo[nh].M_Mean200 = 0;
	      Halo[nh].M_Crit200 = 0;
	      Halo[nh].M_TopHat = 0;
	    }


	  for(i = 0; i < 3; i++)
	    {
	      Halo[nh].Pos[i] = subpos[3 * sc + i];
	      //printf("pos = %f\n", Halo[nh].Pos[i]);
	      Halo[nh].Vel[i] = subvel[3 * sc + i];
	      Halo[nh].Spin[i] = subspin[3 * sc + i];
	    }
	  Halo[nh].VelDisp = subveldisp[sc];
	  Halo[nh].Vmax = subvmax[sc];
	  Halo[nh].MostBoundID = subMostBoundID[sc];


	  /* store position of halo in subfind output */

	  Halo[nh].SnapNum = num;
	  Halo[nh].FileNr = filenrOfHalo[sc];
	  Halo[nh].SubhaloIndex = subhaloindex[sc];
	  Halo[nh].SubhalfMass = subhalfmass[sc];

	  /* auxiliary stuff */

	  HaloAux[nh].UsedFlag = 0;

#ifdef SAVE_MASS_TAB
	  for(i = 0; i < 6; i++)
	    Halo[nh].SubMassTab[i] = submasstab[6 * sc + i];
#endif
	}
    }


  for(gr = 0; gr < nh; gr++)
    {
      if(Halo[gr].NextHaloInFOFgroup == gr)
	{
	  printf("bummer! %d\n", gr);
	}
    }




#ifdef SAVE_MASS_TAB
  myfree(submasstab);
#endif
  myfree(subhalfmass);
  myfree(subMostBoundID);
  myfree(subspin);
  myfree(subvmax);
  myfree(subveldisp);
  myfree(subvel);
  myfree(subpos);

  myfree(halo_M_TopHat);
  myfree(halo_M_Crit200);
  myfree(halo_M_Mean200);

  myfree(subhaloindex);
  myfree(filenrOfHalo);
  myfree(descendant_snapnum);
  myfree(descendant_haloindex);
  myfree(subLen);
  myfree(nsubPerHalo);
}


void set_progenitor_pointers(void)
{
  int i, first, desc;

  for(i = 0; i < TotHalos; i++)
    {
      if((desc = Halo[i].Descendant) >= 0)
	{
	  if((first = Halo[desc].FirstProgenitor) >= 0)
	    {
	      if(Halo[i].Len >= Halo[first].Len)
		{
		  Halo[i].NextProgenitor = first;
		  Halo[desc].FirstProgenitor = i;
		}
	      else
		{
		  Halo[i].NextProgenitor = Halo[first].NextProgenitor;
		  Halo[first].NextProgenitor = i;
		}
	    }
	  else
	    {
	      Halo[desc].FirstProgenitor = i;
	    }
	}
    }
}
