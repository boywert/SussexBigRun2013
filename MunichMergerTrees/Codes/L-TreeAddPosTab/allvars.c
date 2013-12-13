
#include "allvars.h"




struct particle_data *P;


struct io_header header;


int NumPart;

long long *IDs;
float *Pos;

int TotIDs;



int *SnapNumHalo;

int FirstFile;			/* first and last file for processing */
int LastFile;

int Ntrees;			/* number of trees in current file */

int FilesPerSnapshot;
int LastSnapShotNr;
int SnapFormat;

char OutputDir[512];
char SimulationDir[512];
char SnapshotFileBase[512];

int TotHalos;

int *FirstHaloInSnap;



int *TreeNHalos;
int *TreeFirstHalo;
