/*==================================================================================================
 *  MergerTree:   Merger Tree AHF_particles files
 *
 *
 *  input:    - how often to perform
 *            - 2x _particles files
 *
 *  output:   - 1x _mtree file
 *
 *
 * it is checked what halos in file2 make up the halos in file1, i.e.
 *
 *   file1   file2
 *
 *    0        0
 *    0       17
 *    0       31    -> halo #0 in file1 shares particles with halos #0,17,31 in file2
 *    1        2
 *    1       12
 *    1        4    -> halo #1 in file1 shares particles with halos #2,12,4  in file2
 *       etc.
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

#define MINCOMMON      10   // we only cross-correlate haloes if they at least share MINCOMMON particles
#define NSATinHOSTfrac 0.5  // at least NSATinHOSTfrac the subhaloes particles need to be inside the prospective host

//#define DEBUG
//#define MTREE_SELF

// writes output that readily allows to find mergers
//#define MERGER_RATIO   0.25

#define CLUES_WM3

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  long   npart;
  long  *Pid;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  long  nhalos;
  long *Hid;
}PARTS;

typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  long unsigned id[2];
  long unsigned npart[2];
  long unsigned common;
} MTREE;

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

HALOptr halos[2];
PARTptr parts[2];
long    nHalos[2];
long    PidMax=-1;
long    PidMin=1234567890;

MTREE *mtree;                // a memory copy of the actual *_mtree file
long unsigned  nlines;       // number of lines in *_mtree file
long unsigned  nstree;       // for how many halos to write *_stree file
long          *ihost;        // assigns one host to each halo

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int  create_mtree           (char HaloFile0[MAXSTRING], char HaloFile1[MAXSTRING], char OutFile[MAXSTRING]);
int  create_stree           (char OutFile[MAXSTRING]);
int  read_particles         (char filename[MAXSTRING], int isimu);
int  particle_halo_mapping  (int  isimu);
int  cross_correlation      (char OutFile[MAXSTRING]);
int  assign_progenitor      (char OutFile[MAXSTRING]);
int  assign_progenitors     (char OutFile[MAXSTRING]);
void read_mtree             (char *infile);
void assign_host            ();
void write_stree            (char *outfile);

/*==================================================================================================
 * main:
 *
 *       simply a wrapper for successive calls to create_mtree()
 *
 *==================================================================================================*/
int main()
{
  int    i, nFiles, isimu, ihalo, ipart;
  char **HaloFile;
  char **OutFile, streeFile[MAXSTRING];
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  printf("=======================================================================\n");
  printf("  construct a cross-correlation between consecutive *_particles files\n");
  printf("=======================================================================\n");
  printf("\nPlease give number of particles files (default=2):      ");
  scanf("%d", &nFiles);
  printf("%d\n",nFiles);
  
  /* allocate memory for nFiles filenames, each of size MAXSTRING */
  for(i=0; i<nFiles; i++)
   {
    HaloFile = (char **) calloc(nFiles, sizeof(char *));
    OutFile  = (char **) calloc(nFiles, sizeof(char *));
   }
  for(i=0; i<nFiles; i++)
   {
    HaloFile[i] = (char *) calloc(MAXSTRING, sizeof(char));
    OutFile[i]  = (char *) calloc(MAXSTRING, sizeof(char));
   }
  
  /* read input filenames from stdin */
  for(i=0; i<nFiles; i++)
   {
    printf("Please give name of %5d. *_particles file:            ",i+1);
    scanf("%s", HaloFile[i]);
    printf("%s\n",HaloFile[i]);
   }
  
  /* read output filenames from stdin */
#ifdef MTREE_SELF
  for(i=0; i<nFiles; i++)
#else
    for(i=0; i<nFiles-1; i++)
#endif
     {
      printf("Please give prefix for %5d. output file:                 ",i+1);
      scanf("%s", OutFile[i]);
      printf("%s\n\n",OutFile[i]);
     }
  
  
  /*======================================================================*
   *  CREATE CROSS-CORRELATION BETWEEN TWO CONSECUTIVE *_PARTICLES FILES  *
   *======================================================================*/
#ifdef MTREE_SELF
  for(i=0; i<nFiles; i++)
#else
    for(i=0; i<nFiles-1; i++)
#endif
     {
      /* mtree[] is filled by assign_progenitor()! */
      mtree = (MTREEptr) calloc(1, sizeof(MTREE));
      
#ifdef MTREE_SELF
      create_mtree(HaloFile[i],HaloFile[i],OutFile[i]);
      create_stree(OutFile[i]);
#else
      create_mtree(HaloFile[i],HaloFile[i+1],OutFile[i]);
      
      /* generate substructure tree? */
      if(strcmp(HaloFile[i],HaloFile[i+1]) == 0)
        create_stree(OutFile[i]);
#endif
      
      
      for(isimu=0; isimu<2; isimu++)
       {
        for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
          if(halos[isimu][ihalo].Pid != NULL) free(halos[isimu][ihalo].Pid);
        for(ipart=0; ipart<PidMax; ipart++)
          if(parts[isimu][ipart].Hid != NULL) free(parts[isimu][ipart].Hid);
        free(halos[isimu]);
        free(parts[isimu]);
       }
      free(mtree);
     }
  
  printf("STOP\n");
  return(1);
}


/*==================================================================================================
 * create_mtree:
 *
 *  generates the cross-correlation between two halo analysis
 *
 *  the correlation is done based upon the individual particle contents of halos
 *  and hence only requires the *_particles files...
 *
 *==================================================================================================*/
int create_mtree(char HaloFile0[MAXSTRING], char HaloFile1[MAXSTRING], char OutFile[MAXSTRING])
{
  long ihalo, jhalo;
  
  /* read in the two *_particle files to be correlated */
  read_particles(HaloFile0, 0);
  read_particles(HaloFile1, 1);
  
  particle_halo_mapping(0);
  particle_halo_mapping(1);
  
  cross_correlation(OutFile);
  
#ifdef MERGER_RATIO
  assign_progenitors(OutFile); //dumps information about progenitors
#else
  assign_progenitor(OutFile);
#endif
  return(1);
}


/*==================================================================================================
 * read_particles:
 *
 * read the file storing the particle IDs for each halo
 *
 *      nHalos = number of halos found in file
 *      Pid    = id's of all those particles
 *
 *==================================================================================================*/
int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE *fpin;
  char  line[MAXSTRING];
  long  nPartInHalo, nPartInUse, ipart, jpart, ihalo, Pid, Ptype, numGoodHalos, haloid;
  long  PidMin_local=1234567890;
  long  PidMax_local=-1;
  
  fprintf(stderr,"o reading file %s ...",filename);
  
  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
   }
  
  /* reset all variables */
  nHalos[isimu] = 0;
  ihalo         = -1;
  halos[isimu]   = NULL;
  
  
  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
  
  /* for AHF_particles files the first line is numGoodHalos which we can happily ignore */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%ld %ld",&haloid,&nPartInHalo) == 1)
    fgets(line,MAXSTRING,fpin);  
  
  do {
    if(strncmp(line,"#",1) != 0)
     {
      
      /* has a haloid been written before the number of particles (true for GalaxiesGoingMAD)*/
      if(sscanf(line,"%ld %ld",&haloid,&nPartInHalo) == 1)
       {
        /* if not, just get the number of particles */
        sscanf(line,"%ld",&nPartInHalo);
       }
      
      /* found yet another halo */
      ihalo++;
      nHalos[isimu] += 1;
      halos[isimu]   = (HALOptr) realloc(halos[isimu], nHalos[isimu]*sizeof(HALOS));
      
      /* store npart and allocate Pid[] accordingly */
      halos[isimu][ihalo].Pid   = (long *) calloc(nPartInHalo, sizeof(long));
      
      /* read all their id's */
      nPartInUse = 0;
      for(ipart=0; ipart<nPartInHalo; ipart++)
       {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read the particle type, too */
        if(sscanf(line,"%ld %ld",&Pid,&Ptype) == 1)
         {
          /* if not, set Ptype to some random number to be ignored */
          sscanf(line,"%ld",&Pid);
          Ptype = -1;
         }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
        if(Ptype == 1)
         {  
           halos[isimu][ihalo].Pid             = (long *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(long));
           halos[isimu][ihalo].Pid[nPartInUse] = Pid;
           
           if(abs(halos[isimu][ihalo].Pid[nPartInUse]) > PidMax) PidMax = abs(halos[isimu][ihalo].Pid[nPartInUse]);
           if(abs(halos[isimu][ihalo].Pid[nPartInUse]) < PidMin) PidMin = abs(halos[isimu][ihalo].Pid[nPartInUse]);
           if(abs(halos[isimu][ihalo].Pid[nPartInUse]) > PidMax_local) PidMax_local = abs(halos[isimu][ihalo].Pid[nPartInUse]);
           if(abs(halos[isimu][ihalo].Pid[nPartInUse]) < PidMin_local) PidMin_local = abs(halos[isimu][ihalo].Pid[nPartInUse]);
           
           nPartInUse++;
         }
       }
      
      /* store number of particles in halo */
      halos[isimu][ihalo].npart = nPartInUse;
      
#ifdef DEBUG
      fprintf(stderr,"  => halo %ld of simu #%ld contains %ld dm+gas particles (expected %ld particles in total)\n",ihalo,isimu,nPartInUse,nPartInHalo);
#endif
     }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
  
  fclose(fpin);
  
  fprintf(stderr," done (full ID range = %ld -> %ld, local ID range = %ld -> %ld)\n",PidMin,PidMax,PidMin_local,PidMax_local);
  
  return(1);
}


/*==================================================================================================
 * particle_halo_mapping:
 *
 *  for each particle remember to which halo(s) it belongs
 *
 *==================================================================================================*/
int particle_halo_mapping(int isimu)
{
  long ihalo, ipart, jpart;
  
  fprintf(stderr,"o creating particle<->halo mapping for file %d ...",isimu);
  
  parts[isimu] = (PARTptr) calloc(PidMax+1, sizeof(PARTS));
  
  /* recording every halo it belongs to */
  for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
   {
    for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
     {
      ipart = abs(halos[isimu][ihalo].Pid[jpart]);
      
      parts[isimu][ipart].nhalos++;
      parts[isimu][ipart].Hid = (long *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(long));
      
      parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
     }
   }
  
  fprintf(stderr," done\n");
  return(1);
}

/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(char OutFile[MAXSTRING])
{
  long  ihalo, jhalo, khalo, ipart, jpart;
  long *common;
  FILE *fpout;
  char outname[MAXSTRING];
  
  fprintf(stderr,"o generating cross-correlation for %ld haloes ...",nHalos[0]);
  
  sprintf(outname,"%s_mtree",OutFile);
  fpout = fopen(outname,"w");
  if(fpout == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",outname);
    exit(0);
   }
  
  for(ihalo=0; ihalo<nHalos[0]; ihalo++)
   {
    /* common records how many particles ihalo(isimu=0) has in common with khalo(isimu=1) */
    common = (long *) calloc(nHalos[1], sizeof(long));
    
    for(jpart=0; jpart<halos[0][ihalo].npart; jpart++)
     {
      ipart = halos[0][ihalo].Pid[jpart];
      
      /* ipart belongs to nhalos halos in isimu=1 */
      for(jhalo=0; jhalo<parts[1][ipart].nhalos; jhalo++)
       {            
         khalo          = parts[1][ipart].Hid[jhalo];
         common[khalo] += 1;
         
#ifdef DEBUG
         if(common[khalo] > halos[1][khalo].npart)
          {
           long kpart;
           int  ifound=0;
           for(kpart=0; kpart<halos[1][khalo].npart; kpart++)
             if(halos[1][khalo].Pid[kpart] == ipart)
              {
               ifound = 1;
               fprintf(stderr,"found: ihalo=%ld khalo=%ld: %ld %ld\n",ihalo,khalo,ipart,halos[1][khalo].Pid[kpart]);
              }
          }
#endif
       }
     }
    
    /* write info to file */
    for(khalo=0; khalo<nHalos[1]; khalo++)
     {
      if(common[khalo] > MINCOMMON)
        fprintf(fpout,"%12ld  %12ld  %12ld  %12ld  %12ld\n",
                ihalo,
                halos[0][ihalo].npart,
                common[khalo],
                khalo,
                halos[1][khalo].npart);
     }
    
    free(common);
    common = NULL;
   }
  
  fclose(fpout);
  
  fprintf(stderr," done\n");
  return(1);
}


/*==================================================================================================
 * assign_progenitor:
 *
 *  assign a unique progenitor to each halo based upon common particles
 *
 *  update 10/10/2007:
 *       each subhalo shares the most particles with its host :-(
 *       -> tried to fix this issue..
 *
 *==================================================================================================*/
int assign_progenitor(char OutFile[MAXSTRING])
{
  FILE   *fpin, *fpout;
  char    OutFile_mtree[MAXSTRING], OutFile_idx[MAXSTRING], line[MAXSTRING];
  long unsigned   id1, npart1, common, id2, npart2, iline;
  double  xcommon, xnpart1, xnpart2, ratio, ratio_max;
  long    prev_id1, iprog;
  
  fprintf(stderr,"o assigning progenitor ...");
  
  strcpy(OutFile_mtree, OutFile);
  strcat(OutFile_mtree, "_mtree");
  strcpy(OutFile_idx, OutFile_mtree);
  strcat(OutFile_idx, "_idx");
  
  /* open output file */
  fpout = fopen(OutFile_idx,"w");
  if(fpout == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_idx);
    exit(0);
   }
  
  /* read *_mtree file */
  read_mtree(OutFile);
  
  prev_id1  = mtree[0].id[0];  
  iprog     = 0;
  ratio_max = -1000.;
  
  for(iline=0; iline<nlines; iline++)
   {
#ifdef DEBUG
    fprintf(stderr,"%ld of %ld\n",iline,nlines);
#endif
    /* copy information from just read *_mtree file */
    id1    = mtree[iline].id[0];
    id2    = mtree[iline].id[1];
    npart1 = mtree[iline].npart[0];
    npart2 = mtree[iline].npart[1];
    common = mtree[iline].common;
    
    /* calculate progenitor criterion */
    /* on some machines (incl. my old PowerMac running Tiger) a long is only 4 bytes :-( */
    xcommon   = (double)common;
    xnpart1   = (double)npart1;
    xnpart2   = (double)npart2;
    ratio     = pow2(xcommon)/(xnpart1*xnpart2);
    
    /* is current line still for the same id1 halo? */
    if(id1 == prev_id1)
     {
      /*yes: is current id2 a more appropriate progenitor? */
      if(ratio > ratio_max)
       {
        ratio_max = ratio;
        iprog     = id2;
       }
     }
    else
     {
      /*no: dump the most appropriate progenitor to file */
      fprintf(fpout,"%ld %ld\n",prev_id1,iprog);
      fflush(fpout);
      
      /* we already read in the first progenitor for the next id1 halo... */
      prev_id1  = id1;
      iprog     = id2;
      
      /* calculate progenitor criterion */
      xcommon   = (double)common;
      xnpart1   = (double)npart1;
      xnpart2   = (double)npart2;
      ratio     = pow2(xcommon)/(xnpart1*xnpart2);
      ratio_max = ratio;
      
      //         fprintf(stderr,"%ld %ld %ld %ld %ld    %g\n",id1,npart1,common,id2,npart2,ratio);
      
     }      
   }
  
  //  fclose(fpin);
  fclose(fpout);
  
  fprintf(stderr," done\n");
  return(1);
}   

#ifdef MERGER_RATIO
/*==================================================================================================
 * assign_progenitors:
 *
 *  get statistics for multiple progenitors (e.g. mass ratios, etc.)
 *
 *  update 10/10/2007:
 *       each subhalo shares the most particles with its host :-(
 *       -> tried to fix this issue..
 *
 *==================================================================================================*/
int assign_progenitors(char OutFile[MAXSTRING])
{
  FILE   *fpin, *fpout, *fpout_merger;
  char    OutFile_mtree[MAXSTRING], OutFile_idx[MAXSTRING], OutFile_merger[MAXSTRING], line[MAXSTRING];
  long unsigned   *id1, *npart1, *common, *id2, *npart2, *idx, iline;
  double  xcommon, xnpart1, xnpart2, *ratio;
  long    prev_id1, nprog, iprog;
  
  fprintf(stderr,"o assigning progenitors ...");
  
  strcpy(OutFile_mtree, OutFile);
  strcat(OutFile_mtree, "_mtree");
  strcpy(OutFile_idx, OutFile_mtree);
  strcat(OutFile_idx, "_progs");
  strcpy(OutFile_merger, OutFile_mtree);
  strcat(OutFile_merger, "_merger");
  
  /* open output file */
  fpout = fopen(OutFile_idx,"w");
  if(fpout == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_idx);
    exit(0);
   }
  fprintf(fpout,"#id(1) Np(2) iprog1(3) Np1(4) ncommon1(5) iprog2(6) Np2(7) ncommon2(8) iprog3(9) Np3(10) ncommon3(11)\n");
  
  /* open output file */
  fpout_merger = fopen(OutFile_merger,"w");
  if(fpout_merger == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_merger);
    exit(0);
   }
  fprintf(fpout_merger,"#id(1) iprog1(2) iprog2(3) common2/ncommon1(4)\n");

  /* read *_mtree file */
  read_mtree(OutFile);
  
  prev_id1  = mtree[0].id[0];  
  nprog     = 0;
  id1    = (long *) calloc(1, sizeof(long));
  id2    = (long *) calloc(1, sizeof(long));
  npart1 = (long *) calloc(1, sizeof(long));
  npart2 = (long *) calloc(1, sizeof(long));
  common = (long *) calloc(1, sizeof(long));
  
  for(iline=0; iline<nlines; iline++)
   {
    if(mtree[iline].id[0] == prev_id1)
     {
      /* make room for one more progenitor */
      nprog++;
      id1    = (long *) realloc(id1,    (nprog)*sizeof(long));
      id2    = (long *) realloc(id2,    (nprog)*sizeof(long));
      npart1 = (long *) realloc(npart1, (nprog)*sizeof(long));
      npart2 = (long *) realloc(npart2, (nprog)*sizeof(long));
      common = (long *) realloc(common, (nprog)*sizeof(long));
      
      /* copy information from mtree[] */
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
     }
    else
     {
      ratio = (double *)        calloc(nprog+1, sizeof(double));
      idx   = (long unsigned *) calloc(nprog+1, sizeof(long unsigned));
      
      for(iprog=0; iprog<nprog; iprog++)
       {
        /* calculate progenitor criterion */
        xcommon          = (double)common[iprog];
        xnpart1          = (double)npart1[iprog];
        xnpart2          = (double)npart2[iprog];
        ratio[iprog]     = pow2(xcommon)/(xnpart1*xnpart2);
       }
      /* sort all progenitor by merit function */
      indexx(nprog, ratio-1, idx-1);
      
      /*-----------------------------------------------------
       * dump information to file
       *-----------------------------------------------------*/
      if(nprog > 2)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1],
                id2[idx[nprog-3]-1],npart2[idx[nprog-3]-1],common[idx[nprog-3]-1]);
        
        /* iprog2 is most credible second progenitor */
        if(id2[idx[nprog-1]-1]<id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* iprog3 is most credible second progenitor */
        else if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-3]-1],
                    (double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* nothing else to do as only one real progenitor exists */
       }
      
      else if (nprog > 1)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1]);            

        /* this indicates that the halo itself is a subhalo */
        if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1])
         {
          /* nothing to do as only one real progenitor exists */
         }
        else
         {
          /* iprog2 is most credible second progenitor */
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
       }
      else
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       -1 -1 -1       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1]);
        
        /* nothing else to do as there are no multiple progenitors */
       }
      
      
      free(ratio);
      free(idx);
      
      /* start a new progenitor list */
      free(id1);
      free(id2);
      free(npart1);
      free(npart2);
      free(common);
      
      nprog  = 1;
      id1    = (long *) calloc(nprog, sizeof(long));
      id2    = (long *) calloc(nprog, sizeof(long));
      npart1 = (long *) calloc(nprog, sizeof(long));
      npart2 = (long *) calloc(nprog, sizeof(long));
      common = (long *) calloc(nprog, sizeof(long));
      
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
      
      prev_id1        = id1[nprog-1];
      
     }
   }
  
  if(id1 != NULL) free(id1);
  if(id2 != NULL) free(id2);
  if(npart1 != NULL) free(npart1);
  if(npart2 != NULL) free(npart2);
  if(common != NULL) free(common);
  
  fclose(fpout);
  fclose(fpout_merger);
  
  fprintf(stderr," done\n");
  
  return(1);
}
#endif

/*==================================================================================================
 * read_mtree:
 *
 *       simply reads in the *_mtree file and 
 *       puts it into the array of structures mtree[iline].XYZ
 *
 * Note: at this stage we just treat these entries as "lines" -> no connection to halos yet!
 *
 *==================================================================================================*/
void read_mtree(char *prefix)
{
  long unsigned iline;
  char          line[MAXSTRING], outname[MAXSTRING];
  FILE         *fpin;
  
  sprintf(outname,"%s_mtree",prefix);
  if((fpin = fopen(outname,"r")) == NULL)
   {
    fprintf(stderr,"cannot open  %s\nEXIT\n",outname);
    exit(0);
   }
  
  // count number of lines
  nlines = 0;
  while(!feof(fpin))
   {
    nlines++;
    fgets(line,MAXSTRING,fpin);
   }
  
  //fprintf(stderr,"o found %ld lines in:  %s\n",nlines,infile);
  
  // allocate memory
  mtree = (MTREEptr) realloc((MTREEptr)mtree, nlines*sizeof(MTREE));
  
  // actually read the file
  rewind(fpin);
  for(iline=0; iline<nlines; iline++)
   {
    // read next line from file
    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%ld %ld %ld %ld %ld",
           &(mtree[iline].id[0]), 
           &(mtree[iline].npart[0]), 
           &(mtree[iline].common), 
           &(mtree[iline].id[1]), 
           &(mtree[iline].npart[1]));
    //fprintf(stderr,"iline=%ld id[0]=%ld\n",iline,mtree[iline].id[0]);
   }
  
  fclose(fpin);
}

/*==================================================================================================
 * assign_host:
 *
 *       loop through *_mtree file (held in memory in array mtree[].XYZ) and
 *       assign one host to each halo in there; host halo "-1" means isolated halo
 *
 *==================================================================================================*/
void assign_host()
{
  long unsigned iline, ihalo;
  long unsigned npart0, npart1, id0, id1, common;
  
  // the largest halo id
  nstree = (mtree[nlines-1].id[0])+1;
  
  // allocate (and initialize) array to hold host id's
  ihost = (long *) realloc((long *)ihost, nstree*sizeof(long));
  for(ihalo=0; ihalo<nstree; ihalo++) ihost[ihalo] = -1;
  
  // fill ihost array for every possible halo
  for(iline=0; iline<nlines; iline++)
   {
    // for simplicity copy *_mtree values over to local variables
    npart0 = mtree[iline].npart[0];
    npart1 = mtree[iline].npart[1];
    id0    = mtree[iline].id[0];
    id1    = mtree[iline].id[1];
    common = mtree[iline].common;
    
    // the current halo for which we want to set the host
    ihalo  = id1;
    
    if(id0 != id1         &&                 // halos cannot be their own hosts
       ihost[ihalo] < 0   &&                 // did we already assign a host
       common >= NSATinHOSTfrac*npart1  &&   // what fraction of particles should be >>inside<< host
       npart1 < npart0                       // the host has to be more massive
       )
     {
      ihost[ihalo] = id0;
     }
   }
}


/*==================================================================================================
 * write_stree:
 *
 *       simply dump ihost[] array to file *_stree
 *
 *==================================================================================================*/
void write_stree(char *outfile)
{
  long unsigned ihalo;
  FILE         *fpout;
  
  fprintf(stderr,"o writing %s\n",outfile);
  
  if((fpout = fopen(outfile,"w")) == NULL)
   {
    fprintf(stderr,"cannot open  $s\nEXIT\n",outfile);
    exit(0);
   }
  
  for(ihalo=0; ihalo<nstree; ihalo++)
   {
    fprintf(fpout,"%16ld   %16ld\n", ihalo, ihost[ihalo]);
   }
  
  fclose(fpout);
}

/*==================================================================================================
 * create_stree:
 *
 *  assign a host to each halo in the catalogue, -1 means isolated halo
 *
 *==================================================================================================*/
int create_stree(char OutFile[MAXSTRING])
{
  char streeFile[MAXSTRING];
  
  strcpy(streeFile, OutFile);
  strcat(streeFile, "_stree");
  
  ihost = (long *) calloc(1, sizeof(long));
  assign_host();
  
  write_stree(streeFile);
  free(ihost);
  
  return(1);
}
