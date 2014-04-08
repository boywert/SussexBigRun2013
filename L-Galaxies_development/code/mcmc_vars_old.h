
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>



// Variables for the MCMC sampling

int MCMCNpar; //Number of parameters to sample
#define MCMCNConstraints 24  //Nr of Observational Constraints
#define MCMCMaxObsBins 100 //maximum number of bins per observation per redshift

double MCMCConstraintsZZ[NOUT];
#define NgridCosm 50 ;
//Nr of SAM galaxies used to compare with data
int TotMCMCGals[NOUT];
//To allocate structure with SAM Galaxies
#define MCMCAllocFactor 200000
//#define MCMCAllocFactor 1000000
long MCMCseed;
int Nbins[NOUT][MCMCNConstraints]; //bins on each obs test
double lhood1;

int NFofsInSample[NOUT];


#ifdef MR_PLUS_MRII
int NTrees_Switch_MR_MRII;
int Switch_MR_MRII;
#endif
int CurrentMCMCStep;


FILE *FILE_MCMC_LIKELIHOOD;
FILE *FILE_MCMC_PredictionsPerStep[NOUT][MCMCNConstraints];

//READ FROM input.par
char MCMCParameterValues[512];
char MCMCObsConstraints[512];
char ObsConstraintsDir[512];
char MCMCSampleDir[512];
char MCMCSampleFilePrefix[512];
int  MCMCSampleFile;
#ifdef MR_PLUS_MRII
char MCMCSampleFilePrefix_MR[512];
char MCMCSampleFilePrefix_MRII[512];
int  MCMCSampleFile_MR;
int  MCMCSampleFile_MRII;
#endif
int  MCMCSampleFileType;
int  MCMCTreeSampleFile;
char MCMCHaloModelDir[512];
int  ChainLength;
int  Sample_Physical_Parameters;
int  Time_Dependant_PhysPar;
int  Sample_Cosmological_Parameters;
int  MCMCMode;
double MCMC_LogStep_Size;
double MCMC_Initial_Par_Displacement;
double MCMC_Minimum_Obs_Error;


//Structures for the MCMC
struct MCMC_OBSCONSTRAINTS
{
	char Name[1000];
	char TestType[1000];
	//for chi-square like tests
	double Bin[NOUT][MCMCMaxObsBins];
	double Obs[NOUT][MCMCMaxObsBins];
	double Error[NOUT][MCMCMaxObsBins];
  //for binomial like tests
	double ObsUp[NOUT][MCMCMaxObsBins];
	double ObsDown[NOUT][MCMCMaxObsBins];
	int ObsTest_Switch_z[NOUT];
} *MCMC_Obs;

struct MCMC_GALAXY
{
  float StellarMass[NOUT];
  float ColdGas[NOUT];
  float BulgeMass[NOUT];
  float BlackHoleMass[NOUT];
  float Sfr[NOUT];
  float MagU[NOUT];
  float MagB[NOUT];
  float MagV[NOUT];
  float MagJ[NOUT];
  float MagK[NOUT];
  float Magu[NOUT];
  float Magg[NOUT];
  float Magr[NOUT];
  float Magi[NOUT];
  float Magz[NOUT];
  float Weight[NOUT];
#ifdef HALOMODEL
  int fofid[NOUT];
  float fofmass[NOUT];
  float x[NOUT];
  float y[NOUT];
  float z[NOUT];
  int Type[NOUT];
  int ngal[NOUT];
#endif
} *MCMC_GAL;


struct MCMC_PAR
{
    char   Name[1000];
  	double Value[NOUT];
  	double PropValue[NOUT];
  	double PriorMin;
  	double PriorMax;
    char   Type[1000];
    int    Sampling_Switch;
} *MCMC_PAR;

struct MCMC_FOF_struct
{
	//values for the sample of FoF groups
	double Weight[NOUT];
	long long FoFID[NOUT];
	float M_Crit200[NOUT];
	int NGalsInFoF[NOUT];
	int IndexOfCentralGal[NOUT];
} *MCMC_FOF, *MCMC_FOF2;


//Variables to construct the halo model and compute
//the correlation function for a sample of galaxies
#ifdef HALOMODEL

int *HashTable;

int NR;
int massbins;
float minfofmass;
float maxfofmass;

int NTasks,ThisTask;
int* numbymass;
int** indexbymass;

int offset;

double Norm;
double rho_mean,rho_c,Mstar,Rstar,Delta,Delta_invth,delta_c;
double ngal_mean,ncen_mean;
struct M_params { double k; int i; int j; };
struct N_params { int j; };
struct A_params { double R; int i; int j; };
struct mugal_qawo_params { double rs; double alpha; };
double normrho,alpha_ein,rs_ein,beta_ein,gamma_ein;
double alpha1_ein,alpha2_ein,alpha3_ein,alpha4_ein;
double rs1_ein,rs2_ein,rs3_ein,rs4_ein;
double *cutoff_low,*cutoff_high;
double *kPowerTable,*PowerTable;
char rndsuffix[64];
gsl_interp_accel *SigmaAcc,*PowAcc,*nuAcc;
gsl_spline *SigmaSpline,*PowSpline,*nuSpline;
gsl_interp_accel **NgalAcc;
gsl_spline **NgalSpline;
gsl_interp_accel *ellipAcc;
gsl_spline *ellipSpline;
gsl_interp_accel *NewpowAcc;
gsl_spline *NewpowSpline;
gsl_interp_accel *CorrAcc;
gsl_spline *CorrSpline;
const gsl_rng_type *T_rng;
gsl_rng *r_rng;
double fitconst;

double alpharscutoff_low,alpharscutoff_high;
gsl_interp_accel *alphaAcc,*rsAcc;
gsl_spline *alphaSpline,*rsSpline;


#endif //END of variables to construct the halo model







