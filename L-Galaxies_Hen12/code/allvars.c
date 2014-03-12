// TODO add description for all variables that do not have one yet
#include "allvars.h"


struct GALAXY			/* Galaxy data */
 *Gal, *HaloGal;

struct halo_data *Halo, *Halo_Data;

struct halo_aux_data		/* auxiliary halo data */
 *HaloAux;

struct halo_ids_data *HaloIDs, *HaloIDs_Data;


int FirstFile;			/* first and last file for processing */
int LastFile;


double AllocValue_MaxHaloGal;
double AllocValue_MaxGal;
double AllocValue_MaxGalTree;

int Ntrees;			/* number of trees in current file */

int MaxGal;
int NHaloGal, MaxHaloGal;
int NGalTree, MaxGalTree;
int *HaloGalHeap;
int IndexStored;

char PhotDir[512];
char PhotPrefix[50];
char McFile[512];
char FileWithFilterNames[512];
char CoolFunctionsDir[512];
char CosmologyTablesDir[512];
char OutputDir[512];
char FinalOutputDir[512];
char FileNameGalaxies[512];
char SimulationDir[512];
char FileWithOutputRedshifts[512];
char FileWithZList[512];
//variables used to scale to a different cosmology
char FileWithZList_OriginalCosm[512];
double ScalePos, ScaleMass;

#ifdef SPECIFYFILENR
char FileNrDir[512];
int ListInputFilrNr[111];
#endif


int TotHalos;
int TotGalaxies[NOUT];
int *TreeNgals[NOUT];

int LastSnapShotNr;
int LastDarkMatterSnapShot;

int *FirstHaloInSnap;
int *TreeNHalos;
int *TreeFirstHalo;

double MaxMemSize;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

int ThisTask, NTask;


#ifdef GALAXYTREE
int GalCount;
int TotGalCount;
struct galaxy_tree_data *GalTree;
#endif

size_t HighMark;




/* cosmological parameters */
double Omega;
double OmegaLambda;
double Hubble_h;
double PartMass;
double BoxSize;
double Omega_OriginalCosm;
double OmegaLambda_OriginalCosm;
double Hubble_h_OriginalCosm;
double PartMass_OriginalCosm;
double BoxSize_OriginalCosm;
double BaryonFrac;



/* flags */
int StarFormationRecipe;
int FeedbackRecipe;
int EjectionRecipe;
int ReIncorporationRecipe;
int ReionizationOn;
int BlackHoleGrowth;
int AGNrecipeOn;
int DiskRadiusMethod;
int TrackDiskInstability;
int SatelliteRecipe;
int StarBurstRecipe;
int BulgeFormationInMinorMergersOn;
int MetallicityOption;

/* parameters */
double Reionization_z0;
double Reionization_zr;
double SfrEfficiency;
double SfrLawPivotVelocity;
double SfrLawSlope;
double SfrBurstEfficiency;
double SfrBurstSlope;
double Yield;
double RecycleFraction;
double ThreshMajorMerger;
double AgnEfficiency;
double BlackHoleGrowthRate;
double BlackHoleSeedMass;
double BlackHoleAccretionRate;
double BlackHoleCutoffVelocity;
double FeedbackReheatingEpsilon;
double ReheatPreVelocity;
double ReheatSlope;
double FeedbackEjectionEfficiency;
double EjectPreVelocity;
double EjectSlope;
double ReIncorporationFactor;
double ReincZpower;
double ReincVelocitypower;
double FracZtoHot;
double EnergySNcode, EnergySN;
double EtaSNcode, EtaSN;

double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears,
  UnitTime_in_years,
  G,
  Hubble,
  a0, ar;

int ListOutputSnaps[NOUT];
float ListOutputRedshifts[NOUT];


double ZZ[MAXSNAPS];
double AA[MAXSNAPS];
//variable used to scale to a different cosmology
double AA_OriginalCosm[MAXSNAPS];

double Age[MAXSNAPS];

int Zlistlen;

gsl_rng *random_generator;

int NumMergers;

/*  tabulated stuff */

/* fixed-metallicity spectrophotometric model */
/* tables hold magnitues of starburst population as a function of age */

#ifdef STAR_FORMATION_HISTORY
double SFH_t[MAXSNAPS][STEPS][SFH_NBIN];
double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_ibin[MAXSNAPS][STEPS];
#endif //STAR_FORMATION_HISTORY


float AgeTab[ZL_LEN];
float FilterLambda[NMAG + 1];	//wavelength of each filter + 1 To have V-band to calibrate young stars extinction
float RedshiftTab[IZ];
float MagTableZz[NMAG][IZZ][IZ][ZL_LEN];

// dust
long mu_seed;

void *TreeAuxData;

#ifdef UPDATETYPETWO
int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
long long *IdList;
float *PosList, *VelList;
#endif


int Hashbits;
double ScaleFactor;


#ifdef USE_MEMORY_TO_MINIMIZE_IO
char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
size_t offset_auxdata, offset_treedata, offset_dbids;
size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#ifdef OUTPUT_MOMAF_INPUTS
char *ptr_momafdata[NOUT];
size_t offset_momafdata[NOUT], maxstorage_momafdata[NOUT], filled_momafdata[NOUT];
#endif
#endif


#ifdef  GASRECYCLE

/* gas recycle for SF */
float RecGas;
float Frac[ZL_LEN];
#ifdef METALS
struct metals metalcold;
#else
float metalcold;
#endif
#endif

/*H2 formation table */

float Rho[RHO_LEN];
float H2[RHO_LEN][Z_LEN];

/* reionization Okamoto et al. 2008*/
float Reion_z[46],Reion_Mc[46];

FILE *tree_file;
FILE *treeaux_file;
FILE *treedbids_file;
FILE *FdGalTree;
FILE *FdGalDumps[NOUT];



// ~~ predefined metallicities. TBD THIS SHOULD NOT BE HARDCODED HERE
#ifdef M05
float def_zz[IZZ] = { 0.001, 0.01, 0.02, 0.04 };
#else
float def_zz[IZZ] = { 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05 };
#endif

float log10zz[IZZ];
