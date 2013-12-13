
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>

extern int  LastSnapShotNr;
extern int  FirstSnapShotNr;
extern int  SnapSkipFac;

extern char OutputDir[512];
extern char SnapshotFileBase[512];

extern int TotHalos;

extern int  NumberOfOutputFiles;

extern int    *FirstHaloInSnap;
extern struct halo_catalogue
{
      int TotNsubhalos;
      int TotNgroups;
}
*Cats;

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#endif








