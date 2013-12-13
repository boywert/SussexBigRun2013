
#include "allvars.h"


struct halo_data
*Halo;

 int    *Nunique;

 int    FirstFile;   /* first and last file for processing */
 int    LastFile;

 int    Ntrees;     /* number of trees in current file */

 int    FilesPerSnapshot;
 int    LastSnapShotNr;

 char   OutputDir[512];
 char   SimulationDir[512];

 int    TotHalos;

 int    *FirstHaloInSnap;

 int    *TreeNHalos;
 int    *TreeFirstHalo;










