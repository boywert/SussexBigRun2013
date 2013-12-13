#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>


#define  MAXSNAPS 64     /* SIMULATION DEPENDENT!!!  INSERT APPROPRIATE NUMBER */

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
extern int Snaplistlen;

extern struct haloaux_data
{
  long long HaloID;
  long long FileTreeNr;
  long long FirstProgenitor;
  long long LastProgenitor;
  long long NextProgenitor;
  long long Descendant;
  long long FirstHaloInFOFgroup;
  long long NextHaloInFOFgroup;
  double    Redshift;
  int       PeanoKey;
  int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloAux;

extern int *Done;

extern struct halo_data
{
	/* merger tree pointers */
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
	int SnapNum;
	int FileNr;
	int SubhaloIndex;
	float SubHalfMass;
#ifdef SAVE_MASS_TAB
        float SubMassTab[6];
#endif
}
*Halo;


extern int    *Nunique;

extern int    FirstFile;   /* first and last file for processing */
extern int    LastFile;

extern int    Ntrees;     /* number of trees in current file */

extern int    FilesPerSnapshot;
extern int    LastSnapShotNr;

extern char   OutputDir[512];
extern char   SimulationDir[512];
extern char   FileWithSnapList[512];

extern int    TotHalos;

extern int    Hashbits;

extern double BoxSize;

extern int    *FirstHaloInSnap;



extern int    *TreeNHalos;
extern int    *TreeFirstHalo;





#endif





