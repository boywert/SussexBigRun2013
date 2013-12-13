#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>

#ifndef LONGIDS
#define HDF5_ID_TYPE H5T_NATIVE_UINT
#else
#define HDF5_ID_TYPE H5T_NATIVE_ULLONG
#endif

extern struct particle_data
{
  float Pos[3];
  float Vel[3];
  long long ID;
}
*P;


extern struct io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int hashtabsize;
  char fill[84];		/* fills to 256 Bytes */
}
header;


extern int    NumPart;

extern long long *IDs;
extern float *Pos;

extern int    TotIDs;



extern int    *SnapNumHalo;

extern int    FirstFile;   /* first and last file for processing */
extern int    LastFile;

extern int    Ntrees;     /* number of trees in current file */

extern int    FilesPerSnapshot;
extern int    LastSnapShotNr;
extern int    SnapFormat;

extern char   OutputDir[512];
extern char   SimulationDir[512];
extern char   SnapshotFileBase[512];

extern int    TotHalos;

extern int    *FirstHaloInSnap;



extern int    *TreeNHalos;
extern int    *TreeFirstHalo;





#endif





