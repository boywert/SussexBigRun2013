
#include "allvars.h"

int SnapSkipFac;
int SnapshotNum;
long long TotNumPart;

int SnapFormat;

char OutputDir[512];
char SnapshotFileBase[512];

int LastSnapShotNr;

struct halo_catalogue CatA, CatB, CatC;


struct io_header header;
