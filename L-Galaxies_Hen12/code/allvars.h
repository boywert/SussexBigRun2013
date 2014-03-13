/** @file allvars.h
 * @brief TODO add description for all variables that do not have one yet.
 */
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>

#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7


#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))
#define  wrap(x,y) ( (x)>((y)/2.) ? ((x)-(y)) : ((x)<(-(y)/2.)?((x)+(y)):(x)) )
#define  pow2(x)   ((x)*(x))
#define  pow3(x)   ((x)*(x)*(x))

#define  terminate(x) {char termbuf[5000]; sprintf(termbuf, "code termination on task=%d, function %s(), file %s, line %d: %s\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); printf("%s", termbuf); fflush(stdout); endrun(1);}


#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)


#ifdef GALAXYTREE
#define  CORRECTDBFLOAT(x)  ((fabs(x)<(1.e-30) || isnan(x)) ?(0.0):(x))
#else
#define  CORRECTDBFLOAT(x) x
#endif

#define  NMAG 15 		/* Number of filters */

#define  NIDX 10

#ifdef MRII
#define  MAXSNAPS  68     /* Number of snapshots in the dark matter simulation */
#else
//#define  MAXSNAPS  64
#define  MAXSNAPS  62
#endif

#ifdef CUBEP3M
#undef MASSNAPS
#define MAXSNAPS 76
#endif

#define  MAXGALFAC 2.3 /*1.5/2.3 - maximum fraction of satellite without a halo (for memory allocation)*/

#define  STEPS 20		/* Number of integration intervals between two snapshots */

#define  ALLOCPARAMETER 50.  /* new definition !!! THIS HAS TO BE 50 !!! DONT EVER EVER EVER CHANGE !!! */

//To understand the units in the code read through set_units in init.c!!!
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24 // TODO this is read in from input.par
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */

//To understand the units in the code read through set_units in init.c!!!
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifdef GALAXYTREE
#undef  NOUT
#define NOUT MAXSNAPS
#endif


#ifdef STAR_FORMATION_HISTORY
//#define SFH_NMERGE 5
//#define SFH_NBIN 33
#define SFH_NMERGE 3
#define SFH_NBIN 19
/* SFH_ is the reference structure for storing the star formation histories in
 * logarithmic bins. It is computed in init.c generating a binning structure for
 * each snapshot/time step. In the code galaxy structures are adjusted with respect
 * to this structure at each step. */
extern double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present at the lower edge of the bin (code units)
extern double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
extern int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
extern int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
#endif //STAR_FORMATION_HISTORY

#ifdef METALS
struct metals
{
  float type1a;
  float type2;
  float agb;
};
#endif

#ifdef YIELDS
//Individual element histories:
struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif
  float Fe;
};

//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#endif

//Number of interpolated points within the mass ranges for the four types of yield table:
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#define AGB_MASS_NUM 59 //55 //ROB: 59, when going from 0.85 to 7 Msun
#define AGB_Z_NUM 3
#ifdef PORTINARI
#define SNII_MASS_NUM 85  //ROB: 85, from 6 <= M[Msun] <= 120. Change SNII_MIN_MASS and SNII_MAX_MASS for shorter ranges.
#define SNII_Z_NUM 5
#endif
#ifdef CHIEFFI
#define SNII_MASS_NUM 56
#define SNII_Z_NUM 6
#endif
#define SNIA_MASS_NUM 83 //48 //Number increased after extending range to cover M2 masses (07-02-12)

//Mass ranges for the different modes of ejection:
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#define SNIA_MIN_MASS 3.0
#define SNIA_MAX_MASS 8.0 //16.0
#ifdef PORTINARI
#define SNII_MIN_MASS 7.0 //6.0
#define SNII_MAX_MASS 120.0 //Pick from 40.0, 50.0, 60.0, 70.0 or 120.0
#endif
#ifdef CHIEFFI
#define SNII_MIN_MASS 7.0
#define SNII_MAX_MASS 50.0 //40.0
#endif

int ELETOBIGCOUNTA;
int FRACCOUNTA;

//Arrays that yield tables are written to:
float lifetimeMasses[LIFETIME_MASS_NUM];
float lifetimeMetallicities[LIFETIME_Z_NUM];
float lifetimes[LIFETIME_Z_NUM][LIFETIME_MASS_NUM];
float AGBMasses[AGB_MASS_NUM];
float AGBMetallicities[AGB_Z_NUM];
float AGBEjectedMasses[AGB_Z_NUM][AGB_MASS_NUM];
float AGBTotalMetals[AGB_Z_NUM][AGB_MASS_NUM];
float AGBYields[AGB_Z_NUM][11][AGB_MASS_NUM];
float SNIIMasses[SNII_MASS_NUM];
float SNIIMetallicities[SNII_Z_NUM];
float SNIIEjectedMasses[SNII_Z_NUM][SNII_MASS_NUM];
float SNIITotalMetals[SNII_Z_NUM][SNII_MASS_NUM];
float SNIIYields[SNII_Z_NUM][11][SNII_MASS_NUM];
#ifndef DTD
float SNIaMasses[SNIA_MASS_NUM];
float SNIaEjectedMasses[SNIA_MASS_NUM];
float SNIaTotalMetals[SNIA_MASS_NUM];
float SNIaYields[42][SNIA_MASS_NUM];
#else
float SNIaYields[42];
#endif

//Integrated yields arrays:
float NormSNIIMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIIMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIIYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
float NormAGBMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormAGBMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormAGBYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
float NormSNIaMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIaMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIaYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];

//Arrays used to plot SNe rates:
float TheSFH[SFH_NBIN];
float SNIIRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float SNIaRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float AGBRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];

//SNIa parameters:
#define A_FACTOR 0.06 //Fraction of all objects between SNIa_MIN_MASS and SNIA_MAX_MASS that are SN-Ia progenitors
#ifdef DTD
//#define KALPHA 1.4765 //1.59203 //Now set in yield_integrals.c
#define	F316 0.0384 //Integral of the IMF (by number) from 3.0 - 16.0 Msun
#define SNIAEJECMASS 1.2300971 //Total mass (and total metals) ejected by a SNIa explosion in Msun //Value form original yield table (42 elements): 1.3740855. //Value when only considering 11 elements: 1.2300971
#ifdef BIMODALDTD
	#define DTD_NORM 0.903202 //Normalisation constant for Mannucci DTD when using lifetime ranges for stars between 0.85 - 8.0 Msun.
#endif
#ifdef GAUSSIANDTD
	#define TAUCHARAC 1.0 //Characteristic delay time for SNe-Ia (i.e. peak of Gaussian distribution) in Gyrs //default: 2.0
	#define SIGMA_TD 0.2*TAUCHARAC //0.2 for narrow-DTD, 0.5 for wide_DTD
#endif
#ifdef POWERLAWDTD
	#define DTD_NORM 8.30246 //Normalisation constant for Power-law DTD when using lifetime ranges for stars between 0.85 - 8.0 Msun.
	#define DTD_SLOPE -1.12 //Slope of power law, according to Maoz et al. (2012)
#endif
#endif

#endif //YIELDS

/**
 * Galaxy structure for output
 */
#ifdef LIGHT_OUTPUT
struct GALAXY_OUTPUT
{
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
  float DistanceToCentralGal[3];

  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float HotGas; // 10^10/h Msun - Mass in hot gas
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole

  /* magnitudes in various bands */
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif
};
#else
struct GALAXY_OUTPUT
{
#ifdef NO_PROPS_OUTPUTS
#ifdef GALAXYTREE
  long long GalID;		// ID of galaxy, unique within simulation and SAM run.
#endif
#endif
#ifndef NO_PROPS_OUTPUTS
#ifdef GALAXYTREE
  long long GalID; /** ID of galaxy, unique within simulation and SAM run.*/
  long long HaloID; // Unique ID of MPA halo containing this galaxy
#endif
#ifdef MBPID
  long long MostBoundID; // Most bound particle at centre of subhalo last associated with this galaxy.  Put here as want all 8-byte blocks together at top of output record.
#endif
#ifdef GALAXYTREE
  long long FirstProgGal;	// Main progenitor of this galaxy. Also the first progenitor in a linked list representation of the merger tree.
  long long NextProgGal;	// Next progenitor of this galaxy in linked list representation of merger tree
  long long LastProgGal;	// Galaxies with id between this galaxyId and this lastProgenitorId form the merger tree rooted in this galaxy.
  long long FOFCentralGal;	// TODO has been coded, but that code must be tested.
  long long FileTreeNr;
  long long DescendantGal;	// Pointer to the descendant of this galaxy in its merger tree; -1 if there is no descendant
  long long MainLeafId;
  long long TreeRootId;
  long long SubID;
  long long MMSubID; // fofId, the subhaloid of the subhalo at the center of the fof group
  int   PeanoKey; // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
  float Redshift; // redshift of the snapshot where this galaxy resides
#endif
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
#ifndef GALAXYTREE
  int   HaloIndex;
 // long long SubID;
 // long long FirstHaloInFOFgroup;
#endif
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
  float HaloSpin[3];
#endif
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float LookBackTimeToSnap; //The time from a given snapshot to z=0, in years
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float Vel[3]; // km/s - Galaxy Velocities
  int   Len;   
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
  float Vmax; // km/s - Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
  float GasSpin[3]; // Gas Spin
  float StellarSpin[3]; // Stellar Spin
  float InfallVmax; // km/s - Vmax at infall
  int InfallSnap; // Snapnum at infall
  float HotRadius; //Mpc/h - Radius of the hot gas
  /*dynamical friction merger time*/
  float OriMergTime;
  float MergTime;
  float DistanceToCentralGal[3];
  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float HotGas; // 10^10/h Msun - Mass in hot gas
  float EjectedMass; // 10^10/h Msun - Mass in ejected gas
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole
  float BlackHoleGas; // 10^10/h Msun - Mass in BH accretion disk
  /* ICL magnitude and mass*/
  float ICM;            // mass in intra-cluster stars, for type 0,1
#ifdef METALS
  struct metals MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
  struct metals MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
  struct metals MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
  struct metals MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
  struct metals MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
  struct metals MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf; // hot gas metals that come from self
#endif
#else
  float MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
  float MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
  float MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
  float MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
  float MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
  float MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  float MetalsHotGasSelf; // hot gas metals that come from self
#endif
#endif
  /* misc */
  float Sfr;
  float SfrBulge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
  float CosInclination; // cos(angle) between galaxy spin and the z-axis
  int   DisruptOn; // 0: galaxy merged onto merger center; 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center
#ifdef MERGE01
  int   MergeOn;   // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....
#endif  
  float CoolingRadius;  // Q: store this ? (was stored in Delucia20006a)
  float QuasarAccretionRate;
  float RadioAccretionRate;
#endif // NO_PROPS_OUTPUTS

  /* magnitudes in various bands */
#ifdef OUTPUT_REST_MAGS
  float Mag[NMAG]; // rest-frame absolute mags
  float MagBulge[NMAG]; // rest-frame absolute mags for the bulge
  float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
  float MassWeightAge;
#ifdef  POST_PROCESS_MAGS
  float rbandWeightAge;
#endif
  #ifdef ICL
  float MagICL[NMAG];          // rest-frame absolute mags of ICL
#endif
#endif

#ifdef OUTPUT_OBS_MAGS
  float ObsMag[NMAG]; // obs-frame absolute mags
  float ObsMagBulge[NMAG]; // obs-frame absolute mags for the bulge
  float ObsMagDust[NMAG]; // dust-corrected, obs-frame absolute mags
#ifdef ICL
  float ObsMagICL[NMAG];  // observer-frame absolute mags for intra-cluster light
#endif
#ifdef OUTPUT_MOMAF_INPUTS
  float dObsMag[NMAG];
  float dObsMagBulge[NMAG];
  float dObsMagDust[NMAG];
#ifdef ICL
  float dObsMagICL[NMAG];       
#endif	//
#endif	//OUTPUT_MOMAF_INPUTS
#endif	//OUTPUT_OBS_MAGS

#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin currently in use
  float sfh_time[SFH_NBIN]; //time to present at the middle of bin in years.
  float sfh_dt[SFH_NBIN]; //time width of bin in years.
  float sfh_DiskMass[SFH_NBIN];
  float sfh_BulgeMass[SFH_NBIN];
  float sfh_ICM[SFH_NBIN];
#ifdef METALS
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#else
  float sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef YIELDS
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];

  //float DiskMass_elements[ELEMENT_NUM];
  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif
};
#endif //When LIGHT_OUTPUT is not defined

// TODO add documentation, also for all fields
#ifdef OUTPUT_MOMAF_INPUTS
struct MOMAF_INPUTS
{
  long long GalID;
  long long HaloID;
  int   SnapNum;
  float Pos[3];
  float Vel[3];
  float CosInclination;
  float ObsMagBulge[NMAG];
  float ObsMagDust[NMAG];
  float dObsMagBulge[NMAG];
  float dObsMagDust[NMAG];
};
#endif

struct galaxy_tree_data
{
  int HaloGalIndex;
  int IndexStored;
  int SnapNum;
  int GalID;
  int FirstProgGal;
  int NextProgGal;
  int LastProgGal;
  int DescendantGal;
  int MainLeaf;
  int TreeRoot;
  int FOFCentralGal;
  int Done;
}
 *GalTree;

/*Structure with all the data associated with galaxies (this is not the same as the output!)*/
struct GALAXY			/* Galaxy data */
{
  int HeapIndex;
  int GalTreeIndex;
  int NextGalaxy;
#ifdef GALAXYTREE
  int FirstProgGal;
#endif
  int Type;
  int HaloNr;
  long long MostBoundID;
  int SnapNum;
  int CentralGal;
  float CentralMvir;
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3];
  float MergCentralPos[3];
  float Vel[3];
  float Pos_notupdated[3];
  float Vel_notupdated[3];
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
#endif
  float HaloSpin[3];
  float GasSpin[3];
  float StellarSpin[3];
  int   Len;   
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float InfallVmax;
  /*ram pressure*/
  float InfallSnap;
  float CoolingGas;
  float HotRadius;
  /* baryonic reservoirs */
  float ColdGas;
  float BulgeMass;
  float DiskMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float BlackHoleGas;
  /* metals (a la Gabriella) */
#ifdef METALS
  struct metals MetalsColdGas;
  struct metals MetalsBulgeMass;
  struct metals MetalsDiskMass;
  struct metals MetalsHotGas;
  struct metals MetalsEjectedMass;
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf;
#endif
#else
  float MetalsColdGas;
  float MetalsBulgeMass;
  float MetalsDiskMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
#ifdef METALS_SELF
  float MetalsHotGasSelf;
#endif
#endif

  /* misc */
#ifdef SAVE_MEMORY
  float Sfr;
  float SfrBulge;
#else
  float Sfr[NOUT];
  float SfrBulge[NOUT];
#endif
  float StarMerge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
#ifdef GALAXYTREE
  int   DisruptOn;
#endif
  // float halfradius;
  //float periradius;
  float CosInclination; //angle between galaxy spin and the z-axis
  float OriMergTime;
  float MergeSat;
  float MergTime;
  float DistanceToCentralGal[3];
  int MergeOn;
  float CoolingRadius;
  float QuasarAccretionRate;
  float RadioAccretionRate;
  float ICM;
 #ifdef METALS
   struct metals MetalsICM;
 #else
   float MetalsICM;
 #endif
  /* luminosities in various bands */
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
#ifdef ICL
  float ICLLum[NMAG][NOUT];
#endif
  float Lum[NMAG][NOUT];
  float YLum[NMAG][NOUT];
  float LumBulge[NMAG][NOUT];
  float YLumBulge[NMAG][NOUT];
  float LumDust[NMAG][NOUT];
  float MassWeightAge[NOUT];
#endif
#ifdef COMPUTE_OBS_MAGS
  float ObsLum[NMAG][NOUT];
  float ObsYLum[NMAG][NOUT];
  float ObsLumBulge[NMAG][NOUT];
  float ObsYLumBulge[NMAG][NOUT];
  float ObsLumDust[NMAG][NOUT];
#ifdef ICL
  float ObsICL[NMAG][NOUT];
#endif
#ifdef OUTPUT_MOMAF_INPUTS
  float dObsLum[NMAG][NOUT];
  float dObsYLum[NMAG][NOUT];
  float dObsLumBulge[NMAG][NOUT];
  float dObsYLumBulge[NMAG][NOUT];
  float dObsLumDust[NMAG][NOUT];
#endif
#endif
#endif
#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin are currently in use
  double sfh_age; //Time in years of last call to sph_update_bins
  float sfh_dt[SFH_NBIN]; //Size of time interval in units of years
  float sfh_t[SFH_NBIN]; //Time at low-redshift edge of bin in same units
  int sfh_Nbins[SFH_NBIN]; //Number of bins on the time interval
  float sfh_DiskMass[SFH_NBIN]; //Stellar mass in disk, in bin in standard units
  float sfh_BulgeMass[SFH_NBIN]; //Stellar mass in bulge, in bin in standard units
  float sfh_ICM[SFH_NBIN]; //Stellar mass in ICM, in bin in standard units
#ifdef METALS
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#else
  float sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef YIELDS
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];

  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif //YIELDS
} *Gal, *HaloGal;


// TODO add documentation, also for all fields
struct halo_data
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
}
  *Halo, *Halo_Data;


// TODO add documentation, also for all fields
extern struct  halo_ids_data
{
 long long HaloID;
 long long FileTreeNr;
 long long FirstProgenitor;
 long long LastProgenitor;
 long long NextProgenitor;
 long long Descendant;
 long long FirstHaloInFOFgroup;
 long long NextHaloInFOFgroup;
#ifdef MRII
 long long MainLeafID; 
#endif
 double    Redshift;
 int       PeanoKey;
 int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloIDs, *HaloIDs_Data;



// TODO add documentation, also for all fields
struct halo_aux_data  /* auxiliary halo data */
{
	int DoneFlag;
	int HaloFlag;
	int NGalaxies;
	int FirstGalaxy;
	float M_Crit200_Unscaled;
	float Pos_Unscaled[3];
	float Vel_Unscaled[3];
	float Vmax_Unscaled;
	float Spin_Unscaled[3];
}
 *HaloAux;


extern int FirstFile;		/* first and last file for processing */
extern int LastFile;

extern int Ntrees;		/* number of trees in current file */
extern double AllocValue_MaxHaloGal;
extern double AllocValue_MaxGal;
extern double AllocValue_MaxGalTree;

extern int MaxGal;		/* Maximum number of galaxies allowed for Gal[] array */
extern int NHaloGal, MaxHaloGal;
extern int NGalTree, MaxGalTree;
extern int *HaloGalHeap;
extern int IndexStored;

extern int LastSnapShotNr;
extern int LastDarkMatterSnapShot;

extern char PhotDir[512];
extern char PhotPrefix[50];
extern char McFile[512];
extern char FileWithFilterNames[512];
extern char CoolFunctionsDir[512];
extern char CosmologyTablesDir[512];
extern char OutputDir[512];
/* in case a second parameter is given as argument to the code, this will be taken as a
 * temporary outputdir to allow fast I/O. OutputDir will be replaced by this directory
 * and in the end everything will be moved to the FinalOutputDir (original OutputDir
 * given in input.par )*/
extern char FinalOutputDir[512];
extern char FileNameGalaxies[512];
extern char SimulationDir[512];
extern char FileWithOutputRedshifts[512];
extern char FileWithZList[512];
//variable used to scale to a different cosmology
extern char FileWithZList_OriginalCosm[512];
extern double ScalePos;
extern double ScaleMass;

#ifdef SPECIFYFILENR
extern char   FileNrDir[512];
extern int    ListInputFilrNr[111];
#endif

extern int TotHalos;
extern int TotGalaxies[NOUT];
extern int *TreeNgals[NOUT];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

extern void *TreeAuxData;

extern double MaxMemSize;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern int ThisTask, NTask;

#ifdef GALAXYTREE
extern int GalCount;
extern int TotGalCount;
#endif

/* Cosmological parameters */
extern double Omega;
extern double OmegaLambda;
extern double Hubble_h;
extern double PartMass;
extern double BoxSize;
extern double Omega_OriginalCosm;
extern double OmegaLambda_OriginalCosm;
extern double Hubble_h_OriginalCosm;
extern double PartMass_OriginalCosm;
extern double BoxSize_OriginalCosm;
extern double BaryonFrac;

/* flags */
extern int StarFormationRecipe;
extern int FeedbackRecipe;
extern int EjectionRecipe;
extern int ReIncorporationRecipe;
extern int ReionizationOn;
extern int BlackHoleGrowth;
extern int AGNrecipeOn;
extern int DiskRadiusMethod;
extern int TrackDiskInstability;
extern int SatelliteRecipe;
extern int StarBurstRecipe;
extern int BulgeFormationInMinorMergersOn;
extern int MetallicityOption;

/* parameters */
extern double Reionization_z0;
extern double Reionization_zr;
extern double Yield;
extern double RecycleFraction;
extern double ThreshMajorMerger;
extern double SfrEfficiency;
extern double SfrLawPivotVelocity;
extern double SfrLawSlope;
extern double SfrBurstEfficiency;
extern double SfrBurstSlope;
extern double AgnEfficiency;
extern double BlackHoleGrowthRate;
extern double BlackHoleSeedMass;
extern double BlackHoleAccretionRate;
extern double BlackHoleCutoffVelocity;
extern double FeedbackReheatingEpsilon;
extern double ReheatPreVelocity;
extern double ReheatSlope;
extern double FeedbackEjectionEfficiency;
extern double EjectPreVelocity;
extern double EjectSlope;
extern double ReIncorporationFactor;
extern double ReincZpower;
extern double ReincVelocitypower;
extern double FracZtoHot;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

extern double
	UnitLength_in_cm,
	UnitTime_in_s,
	UnitVelocity_in_cm_per_s,
	UnitMass_in_g,
	RhoCrit,
	UnitPressure_in_cgs,
	UnitDensity_in_cgs,
	UnitCoolingRate_in_cgs,
	UnitEnergy_in_cgs,
	UnitTime_in_Megayears, //Using time as stored in the code, this gives Myr/h
	UnitTime_in_years,
	G,
	Hubble,
	a0, ar;

extern int ListOutputSnaps[NOUT];
extern float ListOutputRedshifts[NOUT];

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
//variable used to scale to a different cosmology
extern double AA_OriginalCosm[MAXSNAPS];

extern double Age[MAXSNAPS];

extern int    Zlistlen;

extern gsl_rng *random_generator;


extern int    NumMergers;


/*  tabulated stuff */

/* SSP LOOK UP TABLES
 * tables hold magnitues of starburst population as a function of age */
#ifdef M05
#define ZL_LEN 220		// Age grid of the SSP tables
#define IZZ 4			// Number of Metalicities used
#else
#define ZL_LEN 221		// Age grid of the SSP tables
#define IZZ 6			// Number of Metalicities used
#endif

#ifdef MRII
#define IZ 68			// Number of redshifts used in the tables
#else
#define IZ 64			// Number of redshifts used in the tables
#endif

extern float AgeTab[ZL_LEN];	//table containing the Age grid of the SSP tables
extern float FilterLambda[NMAG + 1];	//wavelength of each filter + 1 for V-band

/*For H2 formation recipe - Not Supported*/
#define RHO_LEN 101
#define Z_LEN 13


/* Variables to hold metallicity and redshift dependent tables */
extern float RedshiftTab[IZ];
extern float MagTableZz[NMAG][IZZ][IZ][ZL_LEN];
extern float IdxTable[2*NIDX][IZZ][ZL_LEN];
extern int   IdxType[NIDX];
extern float Idxw5[NIDX];
extern float Idxw6[NIDX];


#define ExpTauBCBulge 0.5	// constant extinction for young stars in bulges.
#define MUWIDTH  0.2
#define MUCENTER 0.3
extern long mu_seed;

extern size_t HighMark;

#ifdef UPDATETYPETWO
extern int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
extern int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
extern long long *IdList;
extern float *PosList, *VelList;
#endif


extern int Hashbits;
extern double ScaleFactor;	// factor by which to multiply a position to get its ph index (after floring)


#ifdef USE_MEMORY_TO_MINIMIZE_IO
extern char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#ifdef OUTPUT_MOMAF_INPUTS
extern char *ptr_momafdata[NOUT];
extern size_t offset_momafdata[NOUT],maxstorage_momafdata[NOUT],filled_momafdata[NOUT];
#endif
#endif

#endif

#ifdef GASRECYCLE
extern float RecGas;
extern float Frac[ZL_LEN];
#ifdef METALS
extern struct metals metalcold;
#else
extern float metalcold;
#endif
#endif

extern float Rho[RHO_LEN];
extern float H2[RHO_LEN][Z_LEN];

extern float Reion_z[46],Reion_Mc[46];

extern FILE *tree_file;
extern FILE *treeaux_file;
extern FILE *treedbids_file;
extern FILE *FdGalTree;
extern FILE *FdGalDumps[NOUT];



// extract metalicities from init.c::find_interpolated_lum
// and make global. define logs of them
float def_zz[IZZ];
float log10zz[IZZ];
