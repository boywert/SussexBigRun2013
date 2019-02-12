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

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "allvars.h"
#include "proto.h"

#ifndef ALPHA
#define ALPHA (2.0/3)
#endif

#define WEIGHT_FAK (3.0)


MyIDType *IdSnapTable;

int main(int argc, char **argv)
{
#ifdef HAVE_HDF5
  H5open();
#endif

  if(argc != 3)
    {
      printf("\n  usage: L-BaseTree <parameterfile>  <outputnum>\n");
      printf("  <parameterfile>    see readparmeterfile.c\n");
      printf("  <outputnum>        snapshot number\n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);
  SnapshotNum = atoi(argv[2]);

#if defined(_OPENMP)
  int nthreads, tid;
#pragma omp parallel private(nthreads, tid)
  {
    /* Obtain thread number */
    tid = omp_get_thread_num();
    printf("Hello from thread = %d\n", tid);

    /* Only master thread does this */
    if(tid == 0) 
    {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n", nthreads);
    }
  }
  fflush(stdout);
#endif

#ifdef IDS_HAVE_GAPS
  get_id_translation_table();
#else
  get_TotNumPart();
#endif

  printf("loading catalogues...\n");
  fflush(stdout);

  load_subhalo_catalogue(SnapshotNum, &CatA);
  printf("CatA loaded\n");
  fflush(stdout);

  load_subhalo_catalogue(SnapshotNum + 1 * SnapSkipFac, &CatB);
  printf("CatB loaded\n");
  fflush(stdout);

  if(SnapshotNum + 2 * SnapSkipFac <= LastSnapShotNr)
    {
      load_subhalo_catalogue(SnapshotNum + 2 * SnapSkipFac, &CatC);
      printf("CatC loaded\n");
      fflush(stdout);
    }

  printf("done.\n");
  fflush(stdout);


  printf("preparing ID-to-halo tables...\n");
  fflush(stdout);
#ifdef BACKWARD_CHECKING
  prepare_index_list(&CatA);
  printf("index A done.\n");
  fflush(stdout);
#endif
  prepare_index_list(&CatB);
  printf("index B done.\n");
  fflush(stdout);


  if(SnapshotNum + 2 * SnapSkipFac <= LastSnapShotNr)
    {
      prepare_index_list(&CatC);
      printf("index C done.\n");
      fflush(stdout);
    }


  printf("determine_descendants...\n");
  fflush(stdout);
  determine_descendants(&CatA, &CatB, 0, SnapshotNum + 1 * SnapSkipFac);
  printf("desc AB done.\n");
  fflush(stdout);

#ifdef BACKWARD_CHECKING
  determine_descendants(&CatB, &CatA, 1, SnapshotNum);
  printf("desc BA done.\n");
  fflush(stdout);
#endif

  if(SnapshotNum + 2 * SnapSkipFac <= LastSnapShotNr)
    {
      determine_descendants(&CatB, &CatC, 0, SnapshotNum + 2 * SnapSkipFac);
      printf("desc BC done.\n");
      fflush(stdout);

      determine_descendants(&CatA, &CatC, 1, SnapshotNum + 2 * SnapSkipFac);	/* secondary descendant */
      printf("desc AC done.\n");
      fflush(stdout);
    }

  printf("descendants done.\n");
  fflush(stdout);



  if(SnapshotNum + 2 * SnapSkipFac <= LastSnapShotNr)
    {
      printf("decide whether we should take secondary descendant...\n");
      fflush(stdout);

      count_progenitors(&CatA, &CatB);
      printf("progcount AB done\n");
      fflush(stdout);

      count_progenitors(&CatB, &CatC);
      printf("progcount BC done\n");
      fflush(stdout);

      decide_upon_descendant();

      printf("decision made\n");
      fflush(stdout);
    }

#ifdef BACKWARD_CHECKING
  printf("Doing Backward decision ...\n");
  fflush(stdout);
  if(SnapshotNum + 2 * SnapSkipFac > LastSnapShotNr)
    {
      count_progenitors(&CatA, &CatB);
      printf("progcount AB done\n");
      fflush(stdout);
    }
  if(SnapshotNum + 2 * SnapSkipFac <= LastSnapShotNr)
    {
      decide_backwards(&CatA, &CatB);
      printf("Backward decision for AB done.\n");
      fflush(stdout);
    }
#endif

  printf("saving descendants...\n");
  fflush(stdout);

  save_decendant_list();

  printf("saving done.\n");
  fflush(stdout);

  return 0;
}


#ifdef BACKWARD_CHECKING

#ifndef HALO_SIZE_INCREASE_FOR_SWITCHING
#define HALO_SIZE_INCREASE_FOR_SWITCHING 1.5
#endif

void decide_backwards(struct halo_catalogue *catA, struct halo_catalogue *catB)
{
  int i, ic = 0, ict = 0, ptA, p, ifound;

  for(i = 0; i < catB->TotNsubhalos; i++) 
    {
      if(catB->CountProgenitors[i] == 0)	/* select halos with no progenitors */
	{
	  ict++;
	  if(catB->Descendant[i].HaloIndex[1] >= 0)	/* But in reality they have one */
	    {
	      ptA = catB->Descendant[i].HaloIndex[1];
	      /* check if the halo without progenitor's most bound
	         id is part of the descendant found at the previous snap. 
	         Only in  this case, continue ... */
	      ifound = 0;
	      for(p = 0; p < catA->SubLen[ptA]; p++)
		if(catB->IdList[catB->SubOffset[i]] == catA->IdList[catA->SubOffset[ptA] + p])
		  ifound++;
	      if(ifound)
		{
		  printf("    But is has a progenitor ... \n");
		  /* now check if the tow descendent found merge in the next step ... */
		  if(catB->Descendant[i].HaloIndex[0] ==
		     catB->Descendant[catA->Descendant[ptA].HaloIndex[0]].HaloIndex[0])
		    {
		      printf("       Both possible decendents remerge at C (%d / %d)!\n", catB->SubLen[i],
			     catB->SubLen[catA->Descendant[ptA].HaloIndex[0]]);
		      /* only redirect if missed descendent has more particles ... */
		      if(catB->SubLen[i] >
			 catB->SubLen[catA->Descendant[ptA].HaloIndex[0]] * HALO_SIZE_INCREASE_FOR_SWITCHING)
			{
			  printf("redirecting sub %d (%lld) from %d (%lld) to %d (%lld) \n",
				 ptA, catA->IdList[ptA],
				 catA->Descendant[ptA].HaloIndex[0], catA->IdList[catA->SubOffset[ptA]],
				 i, catB->IdList[catB->SubOffset[i]]);
			  fflush(stdout);
			  catA->Descendant[ptA].HaloIndex[0] = i;
			  ic++;
			}
		    }
		}
	    }
	}
    }
  printf("Redirected %d of %d descendents ...\n", ic, ict);
  fflush(stdout);

}
#endif



void decide_upon_descendant(void)
{
  int i, index_b, index_c;
  int count_b, count_c, count_w, count_n;
  double sumpart;

#ifdef SKIP_BY_WEIGHT
  int count_s = 0;
#endif

  count_b = count_c = count_w = count_n = 0;
  sumpart = 0.0;


#if defined(_OPENMP)
#pragma omp parallel for private(index_b, index_c) reduction(+:count_b,count_c,count_w,count_n,sumpart) 
#endif
  for(i = 0; i < CatA.TotNsubhalos; i++)
    {
      index_b = CatA.Descendant[i].HaloIndex[0];
      index_c = CatA.Descendant[i].HaloIndex[1];

      if(index_b >= 0)
	count_b++;

      if(index_b >= 0 && index_c >= 0)
	{
	  if(CatB.CountProgenitors[index_b] > 1 && CatC.CountProgenitors[index_c] == 0)
	    {
	      CatB.CountProgenitors[index_b]--;
	      CatC.CountProgenitors[index_c]++;
	      CatA.Descendant[i].HaloIndex[0] = CatA.Descendant[i].HaloIndex[1];
	      CatA.Descendant[i].SnapNum[0] = CatA.Descendant[i].SnapNum[1];
	      count_c++;
	    }
#ifdef SKIP_BY_WEIGHT
	  else
	    {
	      if(CatA.Descendant[i].Weight[1] / WEIGHT_FAK > CatA.Descendant[i].Weight[0])
		{
		  CatB.CountProgenitors[index_b]--;
		  CatC.CountProgenitors[index_c]++;
		  CatA.Descendant[i].HaloIndex[0] = CatA.Descendant[i].HaloIndex[1];
		  CatA.Descendant[i].SnapNum[0] = CatA.Descendant[i].SnapNum[1];
		  count_c++;
		  count_s++;
		}
	    }
#endif
	}

      if(index_b < 0 && index_c >= 0)
	{
	  CatA.Descendant[i].HaloIndex[0] = CatA.Descendant[i].HaloIndex[1];
	  CatA.Descendant[i].SnapNum[0] = CatA.Descendant[i].SnapNum[1];
	  CatC.CountProgenitors[index_c]++;
	  count_w++;
	}

      if(index_b < 0 && index_c < 0)
	{
	  /*
	     printf("len=%d\n", CatA.SubLen[i]);
	   */
	  sumpart += CatA.SubLen[i];
	  count_n++;
	}
    }

  printf("Out of %d primary descendants, %d have been rerouted to the secondary descendant.\n",
	 count_b, count_c);

  printf("Additionally, %d have been pointed to the secondary because they had no primary.\n", count_w);

#ifdef SKIP_BY_WEIGHT
  printf("Additionally, %d have been pointed to the secondary because the primary have had low weights.\n",
	 count_s);
#endif

  printf("This leaves %d without descendant, of average size = %g particles.\n", count_n, sumpart / count_n);
  fflush(stdout);
}



void count_progenitors(struct halo_catalogue *catA, struct halo_catalogue *catB)
{
  int i;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(i = 0; i < catB->TotNsubhalos; i++)
    catB->CountProgenitors[i] = 0;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(i = 0; i < catA->TotNsubhalos; i++)
    {
      if(catA->Descendant[i].HaloIndex[0] >= 0)
	{
	  catB->CountProgenitors[catA->Descendant[i].HaloIndex[0]]++;
	}
    }
}


struct cand_data
{
  int haloindex;
  float weight;
};


void determine_descendants(struct halo_catalogue *catA, struct halo_catalogue *catB, int entry, int snapnum)
{
  int i, j, ndiff, ncand, haloB, prev, maxlen;
  MyIDType id;
  float weightmax;
  int halomax;
  struct cand_data *candlist, *difflist;

  maxlen = 0;

  for(i = 0; i < catA->TotNsubhalos; i++)
    if(catA->SubLen[i] > maxlen)
      maxlen = catA->SubLen[i];


#if defined(_OPENMP)
#pragma omp parallel private(candlist, difflist, ncand, i, j, id, haloB, ndiff, prev, weightmax, halomax) 
#endif
  {
    candlist = mymalloc(maxlen * sizeof(struct cand_data));
    difflist = mymalloc(maxlen * sizeof(struct cand_data));

#if defined(_OPENMP)
#pragma omp for schedule(dynamic) nowait 
#endif
    for(i = 0; i < catA->TotNsubhalos; i++)
      {
	ncand = 0;
	for(j = 0; j < catA->SubLen[i]; j++)
	  {
	    id = catA->IdList[catA->SubOffset[i] + j];

	    if(id >= 0 && id < TotNumPart)
	      {
		haloB = catB->IdToHalo[id];

		if(haloB >= 0)
		  {
		    candlist[ncand].haloindex = haloB;
		    candlist[ncand].weight = 1.0 / pow(j + 1, ALPHA);
		    ncand++;
		  }
	      }
	    /* else */
	    /*   { */
	    /* 	char buf[100]; */

	    /* 	long_to_str(buf, id); */
	    /* 	printf("bummer! i=%d  id=%s TotumPart=%d\n", i, buf, (int)TotNumPart); */
	    /* 	exit(4); */
	    /*   } */
	  }

	qsort(candlist, ncand, sizeof(struct cand_data), sort_candlist);

	for(j = 0, ndiff = 0, prev = -1; j < ncand; j++)
	  {
	    if(candlist[j].haloindex != prev)
	      {
		ndiff++;
		difflist[ndiff - 1].haloindex = candlist[j].haloindex;
		difflist[ndiff - 1].weight = 0;
	      }
	    difflist[ndiff - 1].weight += candlist[j].weight;
	    prev = candlist[j].haloindex;
	  }

	weightmax = 0;
	halomax = -1;

	for(j = 0; j < ndiff; j++)
	  {
	    if(difflist[j].weight > weightmax)
	      {
		weightmax = difflist[j].weight;
		halomax = difflist[j].haloindex;
	      }
	  }


	if(ndiff > 0 && halomax >= 0)
	  {
	    catA->Descendant[i].HaloIndex[entry] = halomax;
	    catA->Descendant[i].SnapNum[entry] = snapnum;
#ifdef SKIP_BY_WEIGHT
	    catA->Descendant[i].Weight[entry] = weightmax;
#endif
	  }
	else
	  {
	    catA->Descendant[i].HaloIndex[entry] = -1;
	    catA->Descendant[i].SnapNum[entry] = -1;
#ifdef SKIP_BY_WEIGHT
	    catA->Descendant[i].Weight[entry] = -1;
#endif
	  }
      }
  }
}



int sort_twoids_id(const void *a, const void *b)
{
  if(((struct twoids *) a)->id < ((struct twoids *) b)->id)
    return -1;

  if(((struct twoids *) a)->id > ((struct twoids *) b)->id)
    return +1;

  return 0;
}

int sort_twoids_ord(const void *a, const void *b)
{
  if(((struct twoids *) a)->ord < ((struct twoids *) b)->ord)
    return -1;

  if(((struct twoids *) a)->ord > ((struct twoids *) b)->ord)
    return +1;

  return 0;
}


int sort_candlist(const void *a, const void *b)
{
  if(((struct cand_data *) a)->haloindex < ((struct cand_data *) b)->haloindex)
    return -1;

  if(((struct cand_data *) a)->haloindex > ((struct cand_data *) b)->haloindex)
    return +1;

  return 0;
}

int sort_IDType(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}

void prepare_index_list(struct halo_catalogue *cat)
{
  MyIDType id;
  signed long long ii;
  int i, j;

  cat->IdToHalo = mymalloc(sizeof(int) * TotNumPart);

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(ii = 0; ii < TotNumPart; ii++)
    cat->IdToHalo[ii] = -1;

#if defined(_OPENMP)
#pragma omp parallel for private(j,id)
#endif
  for(i = 0; i < cat->TotNsubhalos; i++)
    {
      for(j = 0; j < cat->SubLen[i]; j++)
	{
	  id = cat->IdList[cat->SubOffset[i] + j];

	  if(id >= 0 && id < TotNumPart)
	    cat->IdToHalo[id] = i;
	  /* else */
	  /*   { */
	  /*     char buf[100]; */

	  /*     long_to_str(buf, id); */

	  /*     printf("bummer! i=%d j=%d id=%s id=%d TotNumPart=%d)\n", i, j, buf, (int)id, (int)TotNumPart); */
	  /*     exit(1); */
	  /*   } */
	}
    }
}



void load_subhalo_catalogue(int num, struct halo_catalogue *cat)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount;
  MyIDType idcount;
  char buf[1000];
  unsigned int j;

  nFiles = 1;
  subcount = 0;
  idcount = 0;

  for(i = 0; i < nFiles; i++) {
    sprintf(buf, "%s/groups_%03d/fof_subhalo_tab_%03d.%d.hdf5", OutputDir, num, num, i);
    hid_t   fd, hd, attr, sr, id, dset;
    herr_t  ret; 
    fd = H5Fopen(buf, H5F_ACC_RDONLY,  H5P_DEFAULT);
    if(fd
       < 0) {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }	
    /* if(!(fd = fopen(buf, "r"))) */
    /* 	{ */
    /* 	  printf("can't open file `%s'\n", buf); */
    /* 	  exit(1); */
    /* 	} */
    if(i == 0 || i == nFiles - 1)
      printf("reading '%s'\n", buf);
    if(i == 1)
      printf("...to...\n");
    // Read header
    hd = H5Gopen (fd, "/Header");
      
    attr = H5Aopen(hd, "Ngroups_ThisFile",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &ngroups);
    ret = H5Aclose(attr);
    //my_fread(&ngroups, sizeof(int), 1, fd);
    // printf("ngroup = %d\n",ngroups);

    attr = H5Aopen(hd, "Ngroups_Total",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &cat->TotNgroups);
    ret = H5Aclose(attr);
    //my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
    // printf("TotNgroup = %d\n",cat.TotNgroups);

    attr = H5Aopen(hd, "Nids_ThisFile",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &nids);
    ret = H5Aclose(attr);
    //my_fread(&nids, sizeof(int), 1, fd);
    //printf("nids = %d\n",nids);

    attr = H5Aopen(hd, "Nids_Total",  H5P_DEFAULT);
    ret  = H5Aread(attr, H5T_NATIVE_INT, &cat->TotNids);
    ret = H5Aclose(attr);
    //my_fread(&cat->TotNids, sizeof(long long), 1, fd);
    //printf("TotNids = %d\n",cat.TotNids);

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
    ret  = H5Aread(attr, H5T_NATIVE_INT, &cat->TotNsubhalos);
    ret = H5Aclose(attr);     
    //my_fread(&cat->TotNsubhalos, sizeof(int), 1, fd);
    //printf("TotNsubhalos = %d\n",cat.TotNsubhalos);

    ret = H5Gclose(hd);
      
    if(i == 0)
      {
	cat->IdList = mymalloc(sizeof(MyIDType) * cat->TotNids);
	cat->SubLen = mymalloc(sizeof(int) * cat->TotNsubhalos);
	cat->SubOffset = mymalloc(sizeof(MyIDType) * cat->TotNsubhalos);
	cat->SubParentHalo = mymalloc(sizeof(int) * cat->TotNsubhalos);
	cat->Descendant = mymalloc(sizeof(struct descendant_data) * cat->TotNsubhalos);
	cat->CountProgenitors = mymalloc(sizeof(int) * cat->TotNsubhalos);
      }

    /*       fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupLen  *\/ */
    /*       fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupOffset  *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  GroupMass  *\/ */
    /*       fseek(fd, 3 * sizeof(float) * ngroups, SEEK_CUR);	/\* skip  GroupPos *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_M_Mean200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_R_Mean200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_M_Crit200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_R_Crit200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_M_TopHat200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_R_TopHat200 *\/ */
    /* #ifdef FLAG_GROUP_VELDISP */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_Mean200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_Crit200 *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_TopHat200 *\/ */
    /* #endif */
    /*       fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupContaminationCount *\/ */
    /*       fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  GroupContaminationMass *\/ */
    /*       fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupNsubs *\/ */
    /*       fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/\* skip  GroupFirstsub *\/ */

    sr = H5Gopen (fd, "/Subhalo");
    dset = H5Dopen (sr, "SubhaloLen");
    ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cat->SubLen[subcount]);
    ret = H5Dclose(dset);
    //my_fread(&cat->SubLen[subcount], sizeof(int), nsubhalos, fd);

      
    /* I believe we don't need this in HDF5 - Boyd */
    /* tmp = mymalloc(sizeof(int) * nsubhalos); */
      
    /* my_fread(tmp, sizeof(int), nsubhalos, fd); */
    for(j = 0; j < nsubhalos; j++)
      cat->SubOffset[subcount + j] = 0;	/* copy it to 64 bit if needed */

    /* myfree(tmp); */
    //myfree(buf);
      
    dset = H5Dopen (sr, "SubhaloParent");
    ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cat->SubParentHalo[subcount]);
    ret = H5Dclose(dset);
    
    ret = H5Gclose(sr);
    //my_fread(&cat->SubParentHalo[subcount], sizeof(int), nsubhalos, fd);

    //fclose(fd);

    subcount += nsubhalos;

  
    /* sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d", OutputDir, num, num, i); */
    /* if(!(fd = fopen(buf, "r"))) */
    /* 	{ */
    /* 	  printf("can't open file `%s'\n", buf); */
    /* 	  exit(1); */
    /* 	} */
    /* if(i == 0 || i == nFiles - 1) */
    /* 	printf("reading '%s'\n", buf); */
    /* if(i == 1) */
    /* 	printf("...to...\n"); */

    //my_fread(&ngroups, sizeof(int), 1, fd);
    //my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
    //my_fread(&nids, sizeof(int), 1, fd);
    //my_fread(&cat->TotNids, sizeof(long long), 1, fd);
    //my_fread(&nFiles, sizeof(int), 1, fd);

    // I don't think we need this - Boyd
    //my_fread(&offset, sizeof(int), 1, fd);

    id = H5Gopen (fd, "/IDs");
    dset = H5Dopen (id, "ID");
    ret = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cat->IdList[idcount]);
    ret = H5Dclose(dset);

    //my_fread(&cat->IdList[idcount], sizeof(MyIDType), nids, fd);
    //fclose(fd);
    
    ret = H5Gclose(id);

    idcount += nids;
    ret = H5Fclose(fd);
  }

  long_to_str(buf, cat->TotNids);
  printf("cat->TotNsubhalos = %d\n", cat->TotNsubhalos);
  printf("cat->TotNids = %s\n", buf);
  reassign_ids(cat->TotNids, cat->IdList);
}




void save_decendant_list(void)
{
  int i, *data;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/treedata", OutputDir);
  mkdir(buf, 02755);


  sprintf(buf, "%s/treedata/sub_desc_%03d", OutputDir, SnapshotNum);
  //sprintf(buf, "%s/treedata/sub_desc_sf%d_%03d", OutputDir, SnapSkipFac, SnapshotNum);
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }

  fwrite(&CatA.TotNsubhalos, sizeof(int), 1, fd);

  data = mymalloc(sizeof(int) * CatA.TotNsubhalos);

  for(i = 0; i < CatA.TotNsubhalos; i++)
    data[i] = CatA.Descendant[i].HaloIndex[0];

  fwrite(data, sizeof(int), CatA.TotNsubhalos, fd);

  for(i = 0; i < CatA.TotNsubhalos; i++)
    data[i] = CatA.Descendant[i].SnapNum[0];

  fwrite(data, sizeof(int), CatA.TotNsubhalos, fd);

  fclose(fd);
}


#define SKIP   {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}


void get_TotNumPart(void)
{
  FILE *fd;
  char buf[1000], bufA[100];
  int blksize1, blksize2;
  int i;
  if(SnapFormat != 3)
    {
      /* Binary snapshot */
      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, LastSnapShotNr, SnapshotFileBase,
	      LastSnapShotNr, 0);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open '%s'\n", buf);
	  exit(1);
	}

      SKIP;
      my_fread(&header, sizeof(header), 1, fd);
      SKIP2;
      fclose(fd);
    }
  else
    {
#ifdef HAVE_HDF5
      hid_t file_id;

      /* HDF5 snapshot */
      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5",
	      OutputDir, LastSnapShotNr, SnapshotFileBase,
	      LastSnapShotNr, 0);
      if((file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
	{
	  printf("can't open '%s'\n", buf);
	  exit(1);
	}
      read_hdf5_header(file_id, &header);
      H5Fclose(file_id);
#else
      printf("SnapFormat=3 was selected but code was compiled without HDF5.\n");
      exit(1);
#endif
    }
  TotNumPart = 0;
  for(i=0;i<6;i++)
    TotNumPart += header.npartTotal[i] + (((long long) header.npartTotalHighWord[i]) << (long long) 32);
  long_to_str(bufA, TotNumPart);
  printf("TotNumPart=%s\n", bufA);

}



void get_id_translation_table(void)
{
  FILE *fd;
  char buf[1000], bufA[100], bufB[100];
  int filenr, numfiles, blksize1, blksize2;

  MyIDType i, minID, maxID, Nskip = 0;
#ifdef HAVE_HDF5
  hid_t file_id;
#endif

  printf("reading IDs from last snapshot\n");
  fflush(stdout);

  sprintf(buf, "%s/snapdir_%03d/sorted_id_table", OutputDir, LastSnapShotNr);

  if((fd = fopen(buf, "r")))
    {
      printf("ok, I'm reading '%s'\n", buf);
      fflush(stdout);
      my_fread(&TotNumPart, sizeof(TotNumPart), 1, fd);

      long_to_str(bufA, TotNumPart);
      printf("TotNumPart=%s\n", bufA);
      fflush(stdout);

      IdSnapTable = mymalloc(TotNumPart * sizeof(MyIDType));
      my_fread(IdSnapTable, TotNumPart, sizeof(MyIDType), fd);

      fclose(fd);
    }
  else
    {
      numfiles = 1;

      for(filenr = 0; filenr < numfiles; filenr++)
	{
	  if(SnapFormat == 1)
	    {
	      /* Type 1 binary snapshot */
	      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, LastSnapShotNr, SnapshotFileBase,
		      LastSnapShotNr, filenr);
	      if(!(fd = fopen(buf, "r")))
		{
		  printf("can't open '%s'\n", buf);
		  exit(1);
		}

	      if(filenr == 0 || filenr == numfiles - 1)
		{	  
		  printf("reading IDs from '%s'\n", buf);
		  fflush(stdout);
		}
	      if(filenr == 1)
		printf("...to...\n");
	  
	      SKIP;
	      my_fread(&header, sizeof(header), 1, fd);
	      SKIP2;

	      if(filenr == 0)
		{
		  numfiles = header.num_files;

		  TotNumPart =
		    header.npartTotal[1] + (((long long) header.npartTotalHighWord[1]) << (long long) 32);

		  long_to_str(bufA, TotNumPart);

		  IdSnapTable = mymalloc(TotNumPart * sizeof(MyIDType));
		}


	      SKIP;
	      fseek(fd, blksize1, SEEK_CUR);	/* skip positions */
	      SKIP2;

	      SKIP;
	      fseek(fd, blksize1, SEEK_CUR);	/* skip velocities */
	      SKIP2;

	      
	      SKIP;
	      my_fread(&IdSnapTable[Nskip], sizeof(MyIDType), header.npart[1], fd);
	      Nskip += header.npart[1];

	      fseek(fd,
		    sizeof(MyIDType) * (header.npart[2] + header.npart[3] + header.npart[4] + header.npart[5]),
		    SEEK_CUR);
	      SKIP2;

	      if(blksize1 != blksize2)
		{
		  printf("error in blksizes... strange\n");
		  exit(1);
		}
	      fclose(fd);
	    }
	  else if (SnapFormat==3)
	    {
	      /* HDF5 snapshot */
#ifdef HAVE_HDF5
	      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5", OutputDir, LastSnapShotNr, SnapshotFileBase,
		      LastSnapShotNr, filenr);
	   
	      if((file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
		{
		  printf("can't open '%s'\n", buf);
		  exit(1);
		}
	      
	      if(filenr == 0 || filenr == numfiles - 1)
		{	  
		  printf("reading IDs from '%s'\n", buf);
		  fflush(stdout);
		}
	      if(filenr == 1)
		printf("...to...\n");

	      printf("filenr = %d\n", filenr);

	      read_hdf5_header(file_id, &header);

	      if(filenr == 0)
		{
		  numfiles = header.num_files;

		  TotNumPart =
		    header.npartTotal[1] + (((long long) header.npartTotalHighWord[1]) << (long long) 32);

		  long_to_str(bufA, TotNumPart);

		  IdSnapTable = mymalloc(TotNumPart * sizeof(MyIDType));
		}
	      read_hdf5_dataset(file_id, "/PartType1/ParticleIDs", HDF5_ID_TYPE, &(IdSnapTable[Nskip]));
	      Nskip += header.npart[1];
	      H5Fclose(file_id);
#else
	      printf("SnapFormat=3 but code was compiled without HDF5\n");
	      exit(1);
#endif
	    }
	  else
	    {
	      /* Not implemented */
	      printf("Only SnapFormat=1 or 3 is supported\n");
	      exit(1);
	    }
	}

      printf("TotNumPart=%s\n", bufA);

      printf("IDs read.\n");
      fflush(stdout);

      for(i = 1, minID = maxID = IdSnapTable[0]; i < TotNumPart; i++)
	{
	  if(minID > IdSnapTable[i])
	    minID = IdSnapTable[i];

	  if(maxID < IdSnapTable[i])
	    maxID = IdSnapTable[i];
	}

      long_to_str(bufA, minID);
      long_to_str(bufB, maxID);

      printf("min-ID=%s  max-ID=%s\n", bufA, bufB);

      printf("sorting IDs\n");
      fflush(stdout);

      qsort(IdSnapTable, TotNumPart, sizeof(MyIDType), sort_IDType);

      printf("sorting done\n");
      fflush(stdout);

      sprintf(buf, "%s/snapdir_%03d/sorted_id_table", OutputDir, LastSnapShotNr);

      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't write to '%s'\n", buf);
	  exit(1);
	}

      my_fwrite(&TotNumPart, sizeof(TotNumPart), 1, fd);
      my_fwrite(IdSnapTable, TotNumPart, sizeof(MyIDType), fd);

      fclose(fd);
    }
}



void reassign_ids(MyIDType N, MyIDType * ids)
{
#ifdef IDS_HAVE_GAPS

  long long i, j, offset, NN;
  int tid, nthreads;
  struct twoids *TwoIDs;


  printf("reassign IDs...\n");
  fflush(stdout);

#if defined(_OPENMP)
#pragma omp parallel private(tid, nthreads, offset, NN, i, j, TwoIDs) shared(IdSnapTable)
#endif
  {
#if defined(_OPENMP)
    tid = omp_get_thread_num();
    nthreads = omp_get_max_threads();
    
    offset = tid * (N / nthreads);
    NN = (N / nthreads);
    
    if(nthreads > 1 && tid == (nthreads - 1))
      {
	NN = N - offset;
      }
    printf("tid=%d offset=%lld NN=%lld\n", tid, offset, NN); 

#else
    NN = N;
    offset = 0;
    tid = 0;
    nthreads = 1;
#endif

    TwoIDs = mymalloc(NN * sizeof(struct twoids));


    for(i = 0; i < NN; i++)
      {
	TwoIDs[i].id = ids[i + offset];
	TwoIDs[i].ord = i;
      }

    qsort(TwoIDs, NN, sizeof(struct twoids), sort_twoids_id);

    /* now assign */

    j = 0;
    for(i = 0; i < NN; i++)
      {
	while(IdSnapTable[j] < TwoIDs[i].id && j < (TotNumPart - 1))
	  j++;

	if(IdSnapTable[j] != TwoIDs[i].id)
	  {
	    printf("ID mismatch found?\n");
	    exit(1);
	  }

	TwoIDs[i].id = j;
      }


    /* sort back */
    qsort(TwoIDs, NN, sizeof(struct twoids), sort_twoids_ord);


    for(i = 0; i < NN; i++)
      {
	ids[i + offset] = TwoIDs[i].id;
      }
    myfree(TwoIDs);
  }

  printf("done\n");
  fflush(stdout);

#else
  
  signed long long i;
 
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(i=0; i< N; i++)
    ids[i] -= 1;

#endif

}




void long_to_str(char *s, long long n)
{
  if(n >= 1000000000)
    sprintf(s, "%d%09d", (int) (n / 1000000000), (int) (n % 1000000000));
  else
    sprintf(s, "%d", (int) n);
}



