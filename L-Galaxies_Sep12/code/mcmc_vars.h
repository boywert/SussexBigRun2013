

// Variables for the MCMC sampling

#define Npar 21 //Number of parameters to sample
#define NTests 11  //Nr of Observational Constraints
#define NChiTests 10  //Nr of Chi2 like constraints
#define NBinomTests 1 //Nr of Binomial like constraints

#define NgridCosm 50 ;
//Nr of SAM galaxies used to compare with data
int TotMCMCGals[NOUT];
//To allocate structure with SAM Galaxies
#define MCMCAllocFactor 500000
long MCMCseed;
int Nbins[NOUT][NTests]; //bins on each obs test
double lhood1;
//values for the sample of FoF groups
double *Weight[NOUT];
long long *SampleIDs[NOUT];
int NTreesInSample[NOUT];

//Nr of Obs tests, Nr of chains, Max number of bins in all observations
double **OutputObs[NOUT][NTests];
int CurrentMCMCStep;


//READ FROM input.par
char MCMCParameterValues[512];
char ObsConstraintsDir[512];
char MCMCTreeSampleDir[512];
int  MCMCTreeSampleFile;
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
  //for chi-square like tests
	double Bin[NOUT][NChiTests];
	double Obs[NOUT][NChiTests];
	double Error[NOUT][NChiTests];
  //for binomial like tests
	double ObsUp[NOUT][NBinomTests];
	double ObsDown[NOUT][NBinomTests];
} *MCMC_Obs;

struct MCMC_GALAXY
{
  float StellarMass[NOUT];
  float ColdGas[NOUT];
  float BulgeMass[NOUT];
  float BlackHoleMass[NOUT];
  float Sfr[NOUT];
  float MagB[NOUT];
  float MagV[NOUT];
  float MagK[NOUT];
  float Magu[NOUT];
  float Magg[NOUT];
  float Magr[NOUT];
  float Magi[NOUT];
  float Magz[NOUT];
  float Weight[NOUT];
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


