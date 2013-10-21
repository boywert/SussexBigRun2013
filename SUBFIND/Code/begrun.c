#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"

/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;
#ifdef _OPENMP
  int tid;
#endif
  if(ThisTask == 0)
    {
      /*    printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version()); */
      printf("\nThis is P-Gadget, version %s.\n", GADGETVERSION);
      printf("\nRunning on %d MPI tasks.\n", NTask);
#ifdef _OPENMP
#pragma omp parallel private(tid)
      {
#pragma omp master
	printf("\nUsing %d OpenMP threads\n", omp_get_num_threads());

	tid = omp_get_thread_num();
	/*
	printf("Hello from thread = %d\n", tid);
	*/
      }
#endif
      printf("\nCode was compiled with settings:\n\n");

      output_compile_time_options();

      printf("Size of particle structure       %d  [bytes]\n",(int)sizeof(struct particle_data));
      printf("\nSize of sph particle structure   %d  [bytes]\n",(int)sizeof(struct sph_particle_data));
#ifdef LT_STELLAREVOLUTION
      printf("\nSize of metal particle structure %d  [bytes]\n\n",(int)sizeof(struct met_particle_data));
#endif
#ifdef DETACH_BLACK_HOLES
      printf("\nSize of BH particle structure    %d  [bytes]\n\n",(int)sizeof(struct bh_particle_data));
#endif
    }

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

  mymalloc_init();

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef MODGRAV
  init_modgrav_background_table();
#endif

  set_units();

#ifdef GROWING_DISK_POTENTIAL
  growing_disk_init();
#endif


#ifdef COOLING
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
#ifndef LT_STELLAREVOLUTION
  InitCool();
#endif
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  printf("Initialize chemistry..\n");
  InitChem();
#endif

#if defined(SFR) && !defined(LT_STELLAREVOLUTION)
  init_clouds();
#endif

#ifdef CHEMCOOL
#ifndef CHEMISTRYNETWORK
  terminate("CHEMISTRYNETWORK not specified!\n");
#endif
  All.ChemistryNetwork = CHEMISTRYNETWORK;

  if(ThisTask == 0)
    {
      printf("initializing chemistry...\n");
      fflush(stdout);
    }

  COOLR.tdust = All.InitDustTemp;
  COOLR.deff = All.H2RefDustEff;
  COOLR.abundo = All.OxyAbund;
  COOLR.abundc = All.CarbAbund;
  COOLR.abundsi = All.SiAbund;
  COOLR.abundD = All.DeutAbund;
  COOLR.abundmg = All.MgAbund;
  COOLR.G0 = All.UVField;
  COOLR.phi_pah = All.PhiPAH;
  COOLR.dust_to_gas_ratio = All.DustToGasRatio;
  COOLR.AV_conversion_factor = All.AVConversionFactor;
  COOLR.cosmic_ray_ion_rate = All.CosmicRayIonRate;
  COOLR.redshift = All.InitRedshift;
  COOLR.AV_ext = All.ExternalDustExtinction;
  COOLR.pdv_term = 0.0;
  COOLR.h2_form_ex = All.H2FormEx;
  COOLR.h2_form_kin = All.H2FormKin;
  COOLR.dm_density = 0.0;
  COOLI.iphoto = All.PhotochemApprox;
  COOLI.iflag_mn = All.MNRateFlag;
  COOLI.iflag_ad = All.ADRateFlag;
  COOLI.iflag_atom = All.AtomicFlag;
  COOLI.iflag_3bh2a = All.ThreeBodyFlagA;
  COOLI.iflag_3bh2b = All.ThreeBodyFlagB;
  COOLI.iflag_h3pra = All.H3PlusRateFlag;
  COOLI.idma_mass_option = 0;
  COOLI.no_chem = 0;
  COOLI.irad_heat = All.RadHeatFlag;

  COOLINMO();
  CHEMINMO();
  INIT_TOLERANCES();
  LOAD_H2_TABLE();
  INIT_TEMPERATURE_LOOKUP();

  if(ThisTask == 0)
    {
      printf("initialization of chemistry finished.\n");
      fflush(stdout);
    }
#endif

#ifdef PERIODIC
  ewald_init();
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
  inverse_boxSize = 1. / boxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
  inverse_boxSize_X = 1. / boxSize_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
  inverse_boxSize_Y = 1. / boxSize_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
  inverse_boxSize_Z = 1. / boxSize_Z;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
  All.DecayTime = 1 / All.DecayLength * sqrt((GAMMA - 1) / 2 * GAMMA);
#endif

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */
 
  set_random_numbers();
  
#ifdef PMGRID
  if(RestartFlag != 3 && RestartFlag != 4)
    long_range_init();
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
  long_range_init();
#endif
#endif

#ifdef SUBFIND
  GrNr = -1;
#endif

#ifdef EOS_DEGENERATE
  eos_init(All.EosTable, All.EosSpecies);
#endif

#ifdef NUCLEAR_NETWORK
  network_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses,
	       All.NetworkWeakrates, &All.nd);
  network_workspace_init( &All.nd, &All.nw );
#endif


  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 3 || RestartFlag == 4 || RestartFlag == 5)
    {
      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets
				   all variables in the struct `All'.
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.
				 */

      set_random_numbers();

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MinGasHsmlFractional = all.MinGasHsmlFractional;
#ifdef MODGRAV
      All.MinModGravHsmlFractional_0 = all.MinModGravHsmlFractional_0;
      All.MinModGravHsmlFractional_1 = all.MinModGravHsmlFractional_1;
#endif
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);
      memcpy(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength);

#ifdef TIME_DEP_ART_VISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif

#if defined(MAGNETIC_DISSIPATION) || defined(EULER_DISSIPATION)
      All.ArtMagDispConst = all.ArtMagDispConst;
#ifdef TIME_DEP_MAGN_DISP
      All.ArtMagDispMin = all.ArtMagDispMin;
      All.ArtMagDispSource = all.ArtMagDispSource;
      All.ArtMagDispTime = all.ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
      All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
      All.DivBcleanQ = all.DivBcleanQ;
#endif

#ifdef BLACK_HOLES
      All.BlackHoleMaxAccretionRadius = all.BlackHoleMaxAccretionRadius;
#endif

#ifdef LT_DF_BH
      All.BH_Radiative_Efficiency = all.BH_Radiative_Efficiency;
#endif

#ifdef LT_STELLAREVOLUTION
      All.MaxChemSpreadL = all.MaxChemSpreadL;
#ifdef LT_ADD_GAL_TO_SUB
      strcpy(All.BC_SED_Path,all.BC_SED_Path);
#endif
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.TimebinFile, all.TimebinFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

#ifdef EOS_DEGENERATE
      strcpy(All.EosTable, all.EosTable);
      strcpy(All.EosSpecies, all.EosSpecies);
#endif

#ifdef RELAXOBJECT
      All.RelaxBaseFac = all.RelaxBaseFac;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(All.NetworkRates, all.NetworkRates);
      strcpy(All.NetworkPartFunc, all.NetworkPartFunc);
      strcpy(All.NetworkMasses, all.NetworkMasses);
      strcpy(All.NetworkWeakrates, all.NetworkWeakrates);
      All.nd = all.nd;
      All.nw = all.nw;
#endif

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);
    }

#ifdef LIMIT_HSML
   double meansep;
   //meansep = All.BoxSize / pow(All.TotNumPart, 1.0 / 3.0);
   meansep = All.BoxSize / pow(All.TotN_gas, 1.0 / 3.0);

   MaxGasHsmlFractional = 2.0;

//All.MaxGasHsml = All.MaxGasHsmlFractional * meansep;                                                                             
MaxGasHsml = MaxGasHsmlFractional * meansep;

if(ThisTask == 0)
  printf("Mean gas interparticle separation = %g, MaxGasHsml = %g\n", meansep, MaxGasHsml);
//printf("Mean gas interparticle separation = %g, MaxGasHsml = %g\n", meansep, All.MaxGasHsml);                                    
#endif


  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  open_outputfiles();

#ifdef PMGRID
  long_range_init_regionsize();
#endif

  reconstruct_timebins();

#ifdef COSMIC_RAYS
  int CRpop;

  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
    CR_initialize_beta_tabs(All.CR_Alpha[CRpop], CRpop);
  CR_Tab_Initialize();
#ifdef COSMIC_RAY_TEST
  CR_test_routine();
#endif

#endif

#ifdef TWODIMS
  int i;

  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;
      P[i].Vel[2] = 0;

      P[i].g.GravAccel[2] = 0;

      if(P[i].Type == 0)
	{
	  SphP[i].VelPred[2] = 0;
	  SphP[i].a.HydroAccel[2] = 0;
	}
    }
#endif


#ifdef RADTRANSFER
  if(RestartFlag == 0)
    radtransfer_set_simple_inits();

  rt_get_sigma();
#endif


  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    {
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
    }
  else if(RestartFlag == 1)
    {
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
    }
  else
    {
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);
    }

   All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT
  double coulomb_log;
#endif
#endif
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;
#ifdef TIMEDEPGRAV
  All.Gini = All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

#ifdef DISTORTIONTENSORPS
  /* 5.609589206e23 is the factor to convert from g to GeV/c^2, the rest comes from All.UnitDensity_in_cgs */
  All.UnitDensity_in_Gev_per_cm3 = 5.609589206e23 / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g;
#endif
  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
#ifdef DISTORTIONTENSORPS
      printf("Annihilation radiation units:\n");
      printf("UnitDensity_in_Gev_per_cm3 = %g\n", All.UnitDensity_in_Gev_per_cm3);
#endif

      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifdef LT_STELLAREVOLUTION
  All.MinChemTimeStep /= 86400 * 365 * 1e9;
  /*
     Add here Spreading Length in the case of fixed distance
   */
  All.SnIaEgy /= All.UnitEnergy_in_cgs;
  All.SnIIEgy /= All.UnitEnergy_in_cgs;

#ifdef LT_HOT_EJECTA
  /* All.EgySpecEjecta is supposed to be given in paramfile as Km/sec        */
  /* then we have to transform in erg / g through the 1e10 conversion factor */
  All.EgySpecEjecta = pow(All.EgySpecEjecta, 2) * 1e10 * All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

  Hyd = LT_NMet - 1;
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  init_SN();
#endif


#if defined(SFR) && !defined(LT_STELLAREVOLUTION)
  set_units_sfr();
#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#ifdef NAVIERSTOKES
  /* Braginskii-Spitzer shear viscosity parametrization */
  /* mu = 0.406 * m_p^0.5 * (k_b* T)^(5/2) / e^4 / logLambda  [g/cm/s] */
  /* eta = frac * mu */

  meanweight = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* assuming full ionization */

#ifdef NAVIERSTOKES_CONSTANT
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow(BOLTZMANN * All.ShearViscosityTemperature, 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in cgs units: eta = %g\n", All.NavierStokes_ShearViscosity);

  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#else
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow((meanweight * PROTONMASS * GAMMA_MINUS1), 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */
  /*T = mu*m_p*(gamma-1)/k_b * E * UnitEnergy/UnitMass */

  All.NavierStokes_ShearViscosity *= pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 5. / 2.);	/* now energy can be multiplied later in the internal code units */
  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Variable shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#endif

#ifdef NAVIERSTOKES_BULK
  if(ThisTask == 0)
    printf("Costant bulk viscosity in internal code units: zeta = %g\n", All.NavierStokes_BulkViscosity);
#endif

#ifdef VISCOSITY_SATURATION
  /* calculate ion mean free path assuming complete ionization:
     ion mean free path for hydrogen is similar to that of helium,
     thus we calculate only for hydrogen */
  /* l_i = 3^(3/2)*(k*T)^2 / (4*\pi^(1/2)*ni*(Z*e)^4*lnL) */

  All.IonMeanFreePath = pow(3.0, 1.5) / (4.0 * sqrt(M_PI) * pow(ELECTRONCHARGE, 4) * LOG_LAMBDA);

  All.IonMeanFreePath *= pow(meanweight * PROTONMASS * GAMMA_MINUS1, 2) * pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 2);	/*kT -> u */

  All.IonMeanFreePath /= (HYDROGEN_MASSFRAC / PROTONMASS) *
    (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  /* n_H = rho * Hfr / mp *//* now is cgs units *///changed / to * in front of the unitdensity

  All.IonMeanFreePath *= All.HubbleParam / All.UnitLength_in_cm;
  /* in internal code units */
#endif

#endif

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  /* assuming full ionization */

  coulomb_log = 37.8;
  /* accordin1g to Sarazin's book */

  All.ConductionCoeff *=
    (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif

#endif /* CONDUCTION_CONSTANT */
#endif /* CONDUCTION */

#if defined(CR_DIFFUSION)
  if(All.CR_DiffusionDensZero == 0.0)
    {
      /* Set reference density for CR Diffusion to rhocrit at z=0 */
      All.CR_DiffusionDensZero = 3.0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
    }

  if(All.CR_DiffusionEntropyZero == 0.0)
    {
      All.CR_DiffusionEntropyZero = 1.0e4;
    }

  /* Set reference entropic function to correspond to
     Reference Temperature @ ReferenceDensity */
  if(ThisTask == 0)
    {
      printf("CR Diffusion: T0 = %g\n", All.CR_DiffusionEntropyZero);
    }

  /* convert Temperature value in Kelvin to thermal energy per unit mass
     in internal units, and then to entropy */
  All.CR_DiffusionEntropyZero *=
    BOLTZMANN / (4.0 * PROTONMASS / (3.0 * HYDROGEN_MASSFRAC + 1.0)) *
    All.UnitMass_in_g / All.UnitEnergy_in_cgs / pow(All.CR_DiffusionDensZero, GAMMA_MINUS1);

  /* Change the density scaling, so that the temp scaling is mapped
     onto an entropy scaling that is numerically less expensive */
  All.CR_DiffusionDensScaling += GAMMA_MINUS1 * All.CR_DiffusionEntropyScaling;

  if(ThisTask == 0)
    {
      printf("CR Diffusion: Rho0 = %g -- A0 = %g\n", All.CR_DiffusionDensZero, All.CR_DiffusionEntropyZero);
    }

#endif /* CR_DIFFUSION */


#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);

  /* fac = M200 / Mtot */
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  /* use unsoftened NFW if NFW_Eps=0 */
  if(NFW_Eps > 0.0)
    if(R > Rs * NFW_C)
      R = Rs * NFW_C;

  if(NFW_Eps > 0.0)
    {
      return fac * 4 * M_PI * RhoCrit * Dc *
	(-
	 (Rs * Rs * Rs *
	  (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
	   NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
													NFW_Eps
													* Rs -
													(2 *
													 NFW_Eps
													 -
													 1) *
													(R +
													 Rs) *
													log(R
													    +
													    Rs)
													+
													NFW_Eps
													*
													NFW_Eps
													* (R +
													   Rs)
													*
													log(R
													    +
													    NFW_Eps
													    *
													    Rs)))
	 / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
    }
  else				/* analytic NFW */
    {
      return fac * 4 * M_PI * RhoCrit * Dc *
	(-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
    }
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_details", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(buf, "%sblackhole_details/blackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef CS_FEEDBACK
  sprintf(buf, "%spromotion_%d.txt", All.OutputDir, ThisTask);
  if(!(FdPromotion = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef DISTORTIONTENSORPS
  /* create caustic log file */
  sprintf(buf, "%scaustics_%d.txt", All.OutputDir, ThisTask);
  if(!(FdCaustics = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_TRACK_WINDS
  sprintf(buf, "%s%s%d%s", All.OutputDir, "trckwinds.", ThisTask, ".txt");
  if(!(FdTrackW = fopen(buf, mode)))
    {
      printf("[%d] error in opening file '%s'\n", ThisTask, buf);
      endrun(1);
    }
#endif

#ifdef LT_POPIII_FLAGS
  sprintf(buf, "%s%s%d%s", All.OutputDir, "trckPIII.", ThisTask, ".dat");
  if(!(FdTrackPopIII = fopen(buf, mode)))
    {
      printf("[%d] error in opening file '%s'\n", ThisTask, buf);
      endrun(1);
    }
  PopIIIevents_buffer = (float *) calloc(POPIIIEVENTS_BUFFER_DIM, sizeof(float));
  PopIIIevents_idx = 0;
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimebinFile);
  if(!(FdTimebin = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1],
	  CPU_SymbolImbalance[CPU_TREEWALK1]);
  fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2],
	  CPU_SymbolImbalance[CPU_TREEWALK2]);
  fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1],
	  CPU_SymbolImbalance[CPU_TREEWAIT1]);
  fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2],
	  CPU_SymbolImbalance[CPU_TREEWAIT2]);
  fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND],
	  CPU_SymbolImbalance[CPU_TREESEND]);
  fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV],
	  CPU_SymbolImbalance[CPU_TREERECV]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC],
	  CPU_SymbolImbalance[CPU_TREEMISC]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE],
	  CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
  fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT],
	  CPU_SymbolImbalance[CPU_DENSWAIT]);
  fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC],
	  CPU_SymbolImbalance[CPU_DENSMISC]);
  fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE],
	  CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
  fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT],
	  CPU_SymbolImbalance[CPU_HYDWAIT]);
  fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC],
	  CPU_SymbolImbalance[CPU_HYDMISC]);
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLES],
	  CPU_SymbolImbalance[CPU_BLACKHOLES]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
  fprintf(FdBalance, "\n");

#ifdef SCFPOTENTIAL
  sprintf(buf, "%s%s", All.OutputDir, "scf_coeff.txt");
  if(!(FdSCF = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef RADTRANSFER
  sprintf(buf, "%s%s", All.OutputDir, "radtransfer.txt");
  if(!(FdRad = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "radtransferNew.txt");
  if(!(FdRadNew = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif



#ifdef CS_MODEL
#ifdef CS_TESTS
  sprintf(buf, "%s%s", All.OutputDir, "energy_test.txt");
  if(!(FdEgyTest = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "promotion_test.txt");
  if(!(FdPromTest = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "SN_test.txt");
  if(!(FdSNTest = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#endif


#if defined(CS_SNI) || defined(CS_SNII)
  sprintf(buf, "%s%s", All.OutputDir, "SN.txt");
  if(!(FdSN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef CS_WSS
  printf("\nUsing WSS 08 cooling\n");
#else
  print("\nUsing SD Cooling cooling\n");
#endif

#ifdef CS_SEED
  printf("\nUsing Seeding\n");
#endif

#endif


#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif

#ifdef XXLINFO
  sprintf(buf, "%s%s", All.OutputDir, "xxl.txt");
  if(!(FdXXL = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdXXL, "nstep time ");
#ifdef MAGNETIC
	  fprintf(FdXXL, "<|B|> ");
#ifdef TRACEDIVB
	  fprintf(FdXXL, "max(divB) ");
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "<alpha> ");
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
    }
#endif

#ifdef MODGRAV
  sprintf(buf, "%s%s", All.OutputDir, "modgrav.txt");
  if(!(FdMOG = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
        {
          fprintf(FdMOG, "nstep time H(a) ");
          fprintf(FdMOG, "\n");
          fflush(FdMOG);
        }
    }
#endif


#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "dH dG ");
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
    }
#endif

#ifdef LT_EXTEGY_INFO
  sprintf(buf, "%s%s", All.OutputDir, "extegy.txt");
  if(!(FdExtEgy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_SEv_INFO_DETAILS_onSPREAD
  sprintf(buf, "%s%s", All.OutputDir, "spinfo.txt");
  if(!(FdSPinfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_CharT_INFO
  sprintf(buf, "%s%s", All.OutputDir, "chart.txt");
  if(!(FdCharT = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO

#ifdef LT_SEvDbg
  sprintf(buf, "%s%s", All.OutputDir, "met_sumcheck.txt");
  if(!(FdMetSumCheck = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  sprintf(buf, "%s%s", All.OutputDir, "metals.txt");
  if(!(FdMetals = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn.txt");
  if(!(FdSn = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn_lost.txt");
  if(!(FdSnLost = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef WINDS
  sprintf(buf, "%s%s", All.OutputDir, "winds.txt");
  if(!(FdWinds = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif
#endif


}




/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

#ifdef CAUSTIC_FINDER
  fclose(FdCaustics);		/* needs to be done by everyone */
#endif

#ifdef LT_TRACK_WINDS
  fclose(FdTrackW);
#endif

#ifdef LT_POPIII_FLAGS
  fwrite(PopIIIevents_buffer, PopIIIevents_idx, sizeof(float), FdTrackPopIII);
  fclose(FdTrackPopIII);
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdTimebin);
  fclose(FdBalance);

#ifdef SCFPOTENTIAL
  fclose(FdSCF);
#endif

#ifdef SFR
  fclose(FdSfr);
#endif

#ifdef RADTRANSFER
  fclose(FdRad);
  fclose(FdRadNew);
#endif

#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif

#ifdef XXLINFO
  fclose(FdXXL);
#endif

#ifdef DARKENERGY
  fclose(FdDE);
#endif

#ifdef CS_FEEDBACK
  fclose(FdPromotion);
#endif

#ifdef CS_MODEL
#ifdef CS_TESTS
  fclose(FdEgyTest);
  fclose(FdPromTest);
  fclose(FdSNTest);
#endif
#if defined(CS_SNI) || defined(CS_SNII)
  fclose(FdSN);
#endif
#endif
#ifdef CS_FEEDBACK
  fclose(FdPromotion);
#endif

#ifdef LT_SEv_INFO_DETAILS_onSPREAD
  fclose(FdSPinfo);
#endif

#ifdef LT_CharT_INFO
  fclose(FdCharT);
#endif

#ifdef LT_EXTEGY_INFO
  fclose(FdExtEgy);
#endif

#ifdef LT_CharT_INFO
  fclose(FdCharT);
#endif
#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO
#ifdef LT_SEvDbg
  fclose(FdMetSumCheck);
#endif
  fclose(FdSn);
  fclose(FdSnLost);
  fclose(FdMetals);
#ifdef WINDS
  fclose(FdWinds);
#endif
#endif
#endif


}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

  All.StarformationOn = 0;	/* defaults */

#ifdef COSMIC_RAYS
  char tempBuf[20];

  int CRpop, x;

  double tempAlpha;
#endif


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimebinFile");
      addr[nt] = All.TimebinFile;
      id[nt++] = STRING;
      
      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

#ifdef WINDTUNNEL
      strcpy(tag[nt], "WindGrid");
      addr[nt] = &All.WindGrid;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindVel");
      addr[nt] = &All.WindVel;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindDmean");
      addr[nt] = &All.WindDmean;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEntr");
      addr[nt] = &All.WindEntr;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindDens");
      addr[nt] = &All.WindDens;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindPmass");
      addr[nt] = &All.WindPmass;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindSizeOfInjectionRegion");
      addr[nt] = &All.WindSizeOfInjectionRegion;
      id[nt++] = REAL;

#endif

#ifdef SUBFIND
      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;
#endif

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = REAL;

      //#ifdef LIMIT_HSML
      //strcpy(tag[nt], "MaxGasHsmlFractional");
      //addr[nt] = &All.MaxGasHsmlFractional;
      //id[nt++] = REAL;
      //#endif

#ifdef MODGRAV
      strcpy(tag[nt], "MinModGravHsmlFractional_0");
      addr[nt] = &All.MinModGravHsmlFractional_0;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinModGravHsmlFractional_1");
      addr[nt] = &All.MinModGravHsmlFractional_1;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = REAL;

#ifdef ARTIFICIAL_CONDUCTIVITY
      strcpy(tag[nt], "ArtCondConstant");
      addr[nt] = &All.ArtCondConstant;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;


#ifdef KSPACE_NEUTRINOS
      strcpy(tag[nt], "KspaceNeutrinoSeed");
      addr[nt] = &All.KspaceNeutrinoSeed;
      id[nt++] = INT;

      strcpy(tag[nt], "Nsample");
      addr[nt] = &All.Nsample;
      id[nt++] = INT;

      strcpy(tag[nt], "SphereMode");
      addr[nt] = &All.SphereMode;
      id[nt++] = INT;

      strcpy(tag[nt], "KspaceDirWithTransferfunctions");
      addr[nt] = All.KspaceDirWithTransferfunctions;
      id[nt++] = STRING;

      strcpy(tag[nt], "KspaceBaseNameTransferfunctions");
      addr[nt] = All.KspaceBaseNameTransferfunctions;
      id[nt++] = STRING;

      strcpy(tag[nt], "PrimordialIndex");
      addr[nt] = &All.PrimordialIndex;
      id[nt++] = REAL;

      strcpy(tag[nt], "Sigma8");
      addr[nt] = &All.Sigma8;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaNu");
      addr[nt] = &All.OmegaNu;
      id[nt++] = REAL;

      strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
      addr[nt] = &All.InputSpectrum_UnitLength_in_cm;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

#ifdef START_WITH_EXTRA_NGBDEV
      strcpy(tag[nt], "MaxNumNgbDeviationStart");
      addr[nt] = &All.MaxNumNgbDeviationStart;
      id[nt++] = REAL;
#endif

#ifdef VORONOI_SHAPESCHEME
      strcpy(tag[nt], "VoronoiStiffNess");
      addr[nt] = &All.VoronoiStiffNess;
      id[nt++] = REAL;

      strcpy(tag[nt], "VoronoiRoundNess");
      addr[nt] = &All.VoronoiRoundNess;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

#ifdef DISTORTIONTENSORPS
      strcpy(tag[nt], "DM_velocity_dispersion");
      addr[nt] = &All.DM_velocity_dispersion;
      id[nt++] = REAL;
#endif
#ifdef SCALARFIELD
      strcpy(tag[nt], "ScalarBeta");
      addr[nt] = &All.ScalarBeta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ScalarScreeningLength");
      addr[nt] = &All.ScalarScreeningLength;
      id[nt++] = REAL;
#endif

#ifdef OUTPUTLINEOFSIGHT
      strcpy(tag[nt], "TimeFirstLineOfSight");
      addr[nt] = &All.TimeFirstLineOfSight;
      id[nt++] = REAL;
#endif


#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
      strcpy(tag[nt], "ReferenceGasMass");
      addr[nt] = &All.ReferenceGasMass;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES
      strcpy(tag[nt], "FractionSpitzerViscosity");
      addr[nt] = &All.FractionSpitzerViscosity;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_CONSTANT
      strcpy(tag[nt], "ShearViscosityTemperature");
      addr[nt] = &All.ShearViscosityTemperature;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_BULK
      strcpy(tag[nt], "NavierStokes_BulkViscosity");
      addr[nt] = &All.NavierStokes_BulkViscosity;
      id[nt++] = REAL;
#endif

#ifdef CHEMISTRY
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = REAL;
#endif

#ifdef CONDUCTION
      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeConductionStep");
      addr[nt] = &All.MaxSizeConductionStep;
      id[nt++] = REAL;
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleTimeInterval");
      addr[nt] = &All.BubbleTimeInterval;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "FirstBubbleRedshift");
      addr[nt] = &All.FirstBubbleRedshift;
      id[nt++] = REAL;
#endif

#ifdef MULTI_BUBBLES
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "ClusterMass200");
      addr[nt] = &All.ClusterMass200;
      id[nt++] = REAL;

      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;

#endif

#ifdef BH_BUBBLES
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadioTriggeringFactor");
      addr[nt] = &All.BlackHoleRadioTriggeringFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "DefaultICMDensity");
      addr[nt] = &All.DefaultICMDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioFeedbackFactor");
      addr[nt] = &All.RadioFeedbackFactor;
      id[nt++] = REAL;
#ifdef UNIFIED_FEEDBACK
      strcpy(tag[nt], "RadioThreshold");
      addr[nt] = &All.RadioThreshold;
      id[nt++] = REAL;
#endif
#endif

#ifdef COSMIC_RAYS
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  sprintf(tempBuf, "CR_SpectralIndex_%i", CRpop);
	  strcpy(tag[nt], tempBuf);
	  addr[nt] = &All.CR_Alpha[CRpop];
	  id[nt++] = REAL;
	}

      /* sort the All.CR_Alpha[CRpop] array in assending order */
      if(NUMCRPOP > 1)
	for(x = 0; x < NUMCRPOP - 1; x++)
	  for(CRpop = 0; CRpop < NUMCRPOP - x - 1; CRpop++)
	    if(All.CR_Alpha[CRpop] > All.CR_Alpha[CRpop + 1])
	      {
		tempAlpha = All.CR_Alpha[CRpop];
		All.CR_Alpha[CRpop] = All.CR_Alpha[CRpop + 1];
		All.CR_Alpha[CRpop + 1] = tempAlpha;
	      }

      strcpy(tag[nt], "CR_SupernovaEfficiency");
      addr[nt] = &All.CR_SNEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaSpectralIndex");
      addr[nt] = &All.CR_SNAlpha;
      id[nt++] = REAL;

#if defined(CR_DIFFUSION)
      strcpy(tag[nt], "CR_DiffusionCoeff");
      addr[nt] = &All.CR_DiffusionCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionDensityScaling");
      addr[nt] = &All.CR_DiffusionDensScaling;
      id[nt++] = REAL;

      /* CR Diffusion scaling: reference density rho_0.
         If value is 0, then rho_crit @ z=0 is used. */
      strcpy(tag[nt], "CR_DiffusionReferenceDensity");
      addr[nt] = &All.CR_DiffusionDensZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTemperatureScaling");
      addr[nt] = &All.CR_DiffusionEntropyScaling;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionReferenceTemperature");
      addr[nt] = &All.CR_DiffusionEntropyZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionMaxSizeTimestep");
      addr[nt] = &All.CR_DiffusionMaxSizeTimestep;
      id[nt++] = REAL;
#endif /* CR_DIFFUSION */

#ifdef CR_SHOCK
      strcpy(tag[nt], "CR_ShockEfficiency");
      addr[nt] = &All.CR_ShockEfficiency;
      id[nt++] = REAL;
#if ( CR_SHOCK == 1 )		/* Constant Spectral Index method */
      strcpy(tag[nt], "CR_ShockSpectralIndex");
      addr[nt] = &All.CR_ShockAlpha;
      id[nt++] = REAL;
#else /* Mach-Number - Dependent method */
      strcpy(tag[nt], "CR_ShockCutoffFac");
      addr[nt] = &All.CR_ShockCutoff;
      id[nt++] = REAL;
#endif
#endif /* CR_SHOCK */

#ifdef FIX_QINJ
      strcpy(tag[nt], "Shock_Fix_Qinj");
      addr[nt] = &All.Shock_Fix_Qinj;
      id[nt++] = REAL;
#endif

#ifdef CR_BUBBLES
      strcpy(tag[nt], "CR_AGNEff");
      addr[nt] = &All.CR_AGNEff;
      id[nt++] = REAL;
#endif
#endif /* COSMIC_RAYS */


#ifdef MACHNUM
      strcpy(tag[nt], "Shock_LengthScale");
      addr[nt] = &All.Shock_Length;
      id[nt++] = REAL;

      strcpy(tag[nt], "Shock_DeltaDecayTimeMax");
      addr[nt] = &All.Shock_DeltaDecayTimeMax;
      id[nt++] = REAL;
#endif

#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = REAL;
#endif

#ifdef BLACK_HOLES

#ifdef DETACH_BLACK_HOLES
      strcpy(tag[nt], "BHFormationFactor");
      addr[nt] = &All.BHfactor;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;


      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleNgbFactor");
      addr[nt] = &All.BlackHoleNgbFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
      addr[nt] = &All.BlackHoleMaxAccretionRadius;
      id[nt++] = REAL;

#ifdef FOF
      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

#ifdef LT_DF_BH_BHAR_SWITCH
      strcpy(tag[nt], "BH_radio_treshold");
      addr[nt] = &All.BH_radio_treshold;
      id[nt++] = REAL;
#endif

#ifdef LT_DF_BH
      strcpy(tag[nt], "BH_Radiative_Efficiency");
      addr[nt] = &All.BH_Radiative_Efficiency;
      id[nt++] = REAL;
#endif

#ifdef LT_BH_ACCRETE_SLICES
      strcpy(tag[nt], "BH_N_AccretionSlices");
      addr[nt] = &All.NBHslices;
      id[nt++] = INT;
#endif

#endif  /* BLACK_HOLES */

#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
      /* read the composition from the parameter file */
      strcpy(tag[nt], "START_elec");
      addr[nt] = &All.Startelec;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HI");
      addr[nt] = &All.StartHI;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HII");
      addr[nt] = &All.StartHII;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HM");
      addr[nt] = &All.StartHM;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HeI");
      addr[nt] = &All.StartHeI;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HeII");
      addr[nt] = &All.StartHeII;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HeIII");
      addr[nt] = &All.StartHeIII;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_H2I");
      addr[nt] = &All.StartH2I;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_H2II");
      addr[nt] = &All.StartH2II;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HD");
      addr[nt] = &All.StartHD;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_DI");
      addr[nt] = &All.StartDI;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_DII");
      addr[nt] = &All.StartDII;
      id[nt++] = REAL;

      strcpy(tag[nt], "START_HeHII");
      addr[nt] = &All.StartHeHII;
      id[nt++] = REAL;
#endif

#ifdef SFR
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

#ifndef LT_STELLAREVOLUTION
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;
      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;

#ifdef VARIABLE_WINDS
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "VariableWindSpecMomentum");
      addr[nt] = &All.VariableWindSpecMomentum;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloConcentrationNorm");
      addr[nt] = &All.HaloConcentrationNorm;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloConcentrationSlope");
      addr[nt] = &All.HaloConcentrationSlope;
      id[nt++] = REAL;
#endif

#endif

#if defined(SNIA_HEATING)
      strcpy(tag[nt], "SnIaHeatingRate");
      addr[nt] = &All.SnIaHeatingRate;
      id[nt++] = REAL;
#endif


#ifdef CS_MODEL
      strcpy(tag[nt], "FactorSFR");
      addr[nt] = &All.FactorSFR;
      id[nt++] = REAL;

      strcpy(tag[nt], "DecouplingParam");
      addr[nt] = &All.DecouplingParam;
      id[nt++] = REAL;


      strcpy(tag[nt], "TlifeSNII");
      addr[nt] = &All.TlifeSNII;
      id[nt++] = REAL;

      strcpy(tag[nt], "Raiteri_TlifeSNII");
      addr[nt] = &All.Raiteri_TlifeSNII;
      id[nt++] = INT;

      strcpy(tag[nt], "MinTlifeSNI");
      addr[nt] = &All.MinTlifeSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxTlifeSNI");
      addr[nt] = &All.MaxTlifeSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "RateSNI");
      addr[nt] = &All.RateSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "SN_Energy_cgs");
      addr[nt] = &All.SN_Energy_cgs;
      id[nt++] = REAL;

      strcpy(tag[nt], "Tcrit_Phase");
      addr[nt] = &All.Tcrit_Phase;
      id[nt++] = REAL;

      strcpy(tag[nt], "DensFrac_Phase");
      addr[nt] = &All.DensFrac_Phase;
      id[nt++] = REAL;

      strcpy(tag[nt], "SN_Energy_frac_cold");
      addr[nt] = &All.SN_Energy_frac_cold;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxHotHsmlParam");
      addr[nt] = &All.MaxHotHsmlParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitialHotHsmlFactor");
      addr[nt] = &All.InitialHotHsmlFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxNumHotNgbDeviation");
      addr[nt] = &All.MaxNumHotNgbDeviation;
      id[nt++] = REAL;
#endif

#ifdef SOFTEREQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;
#endif
#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef MODGRAV
      strcpy(tag[nt], "ModifiedGravityFile");
      addr[nt] = All.ModifiedGravityFile;
      id[nt++] = STRING;
#endif

#ifdef TIME_DEP_ART_VISC
      strcpy(tag[nt], "ViscositySourceScaling");
      addr[nt] = &All.ViscSource0;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityDecayLength");
      addr[nt] = &All.DecayLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityAlphaMin");
      addr[nt] = &All.AlphaMin;
      id[nt++] = REAL;
#endif

#if defined(MAGNETIC_DISSIPATION) || defined(EULER_DISSIPATION)
      strcpy(tag[nt], "ArtificialMagneticDissipationConstant");
      addr[nt] = &All.ArtMagDispConst;
      id[nt++] = REAL;

#ifdef TIME_DEP_MAGN_DISP
      strcpy(tag[nt], "ArtificialMagneticDissipationMin");
      addr[nt] = &All.ArtMagDispMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationSource");
      addr[nt] = &All.ArtMagDispSource;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationDecaytime");
      addr[nt] = &All.ArtMagDispTime;
      id[nt++] = REAL;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      strcpy(tag[nt], "DivBcleaningParabolicSigma");
      addr[nt] = &All.DivBcleanParabolicSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
      addr[nt] = &All.DivBcleanHyperbolicSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "DivBcleaningQ");
      addr[nt] = &All.DivBcleanQ;
      id[nt++] = REAL;
#endif

#ifdef MAGNETIC
#ifdef BINISET
      strcpy(tag[nt], "BiniX");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif

#if defined(BSMOOTH)
      strcpy(tag[nt], "BSmoothInt");
      addr[nt] = &All.BSmoothInt;
      id[nt++] = INT;

      strcpy(tag[nt], "BSmoothFrac");
      addr[nt] = &All.BSmoothFrac;
      id[nt++] = REAL;

#ifdef SETMAINTIMESTEPCOUNT
      strcpy(tag[nt], "MainTimestepCount");
      addr[nt] = &All.MainTimestepCountIni;
      id[nt++] = INT;
#endif
#endif

#ifdef MAGNETIC_DIFFUSION
      strcpy(tag[nt], "MagneticEta");
      addr[nt] = &All.MagneticEta;
      id[nt++] = REAL;
#endif
#endif

#ifdef EOS_DEGENERATE
      strcpy(tag[nt], "EosTable");
      addr[nt] = All.EosTable;
      id[nt++] = STRING;

      strcpy(tag[nt], "EosSpecies");
      addr[nt] = All.EosSpecies;
      id[nt++] = STRING;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(tag[nt], "NetworkRates");
      addr[nt] = All.NetworkRates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkPartFunc");
      addr[nt] = All.NetworkPartFunc;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkMasses");
      addr[nt] = All.NetworkMasses;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkWeakrates");
      addr[nt] = All.NetworkWeakrates;
      id[nt++] = STRING;
#endif

#ifdef RELAXOBJECT
      strcpy(tag[nt], "RelaxBaseFac");
      addr[nt] = &All.RelaxBaseFac;
      id[nt++] = REAL;
#endif

#ifdef RADTRANSFER
      strcpy(tag[nt], "IonizingLumPerSolarMass");
      addr[nt] = &All.IonizingLumPerSolarMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "IonizingLumPerSFR");
      addr[nt] = &All.IonizingLumPerSFR;
      id[nt++] = REAL;

#endif

#ifdef SINKS
      strcpy(tag[nt], "SinkHsml");
      addr[nt] = &All.SinkHsml;
      id[nt++] = REAL;

      strcpy(tag[nt], "SinkDensThresh");
      addr[nt] = &All.SinkDensThresh;
      id[nt++] = REAL;
#endif

#ifdef BP_REAL_CRs
      strcpy(tag[nt], "MinCREnergy");
      addr[nt] = &All.ecr_min;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxCREnergy");
      addr[nt] = &All.ecr_max;
      id[nt++] = REAL;

#ifdef BP_SEED_CRs
      strcpy(tag[nt], "CRpInitSlope");
      addr[nt] = &All.pSlope_init;
      id[nt++] = REAL;

      strcpy(tag[nt], "CReInitSlope");
      addr[nt] = &All.eSlope_init;
      id[nt++] = REAL;
#endif
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = REAL;
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_STOP_COOL_BELOW_Z
      strcpy(tag[nt], "CoolStop_redshift");
      addr[nt] = &All.Below_this_redshift_stop_cooling;
      id[nt++] = REAL;
#endif

#ifdef LT_TRACK_CONTRIBUTES
      strcpy(tag[nt], "PackPowerBase");
      addr[nt] = &PowerBase;
      id[nt++] = INT;
#endif

#ifdef LT_STARBURSTS
      strcpy(tag[nt], "StarBurstCondition");
      addr[nt] = &All.StarBurstCondition;
      id[nt++] = INT;

      strcpy(tag[nt], "SB_Density_Thresh");
      addr[nt] = &All.SB_Density_Thresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "SB_DEntropy_Thresh");
      addr[nt] = &All.SB_DEntropy_Thresh;
      id[nt++] = REAL;
#endif

#ifdef LT_SMOOTH_Z
#if defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)
      strcpy(tag[nt], "SmoothRegionSize");
      addr[nt] = &All.SmoothRegionSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "SmoothRegionSizeMax");
      addr[nt] = &All.SmoothRegionSizeMax;
      id[nt++] = REAL;
#endif
#if defined(LT_SMOOTH_NGB) && !defined(LT_SMOOTH_SIZE)
      strcpy(tag[nt], "DesNumNgnSmooth");
      addr[nt] = &All.DesNumNgbSmooth;
      id[nt++] = INT;
#endif
#endif

      strcpy(tag[nt], "TestSuite");
      addr[nt] = &All.TestSuite;
      id[nt++] = INT;

      /* the minimum number of neighbours used to spread
       * metals from Stars.
       */
      strcpy(tag[nt], "InfNeighNum");
      addr[nt] = &All.NeighInfNum;
      id[nt++] = INT;

      /* the desired number of neighbours used to spread
       * metals from Stars.
       */
      strcpy(tag[nt], "DesNumNgbSN");
      addr[nt] = &All.DesNumNgbSN;
      id[nt++] = INT;

      /* the allowed deviation on number of neighbours used
       * to spread metals from Stars.
       */
      strcpy(tag[nt], "SpreadNumNgbDev");
      addr[nt] = &All.SpreadNumNgbDev;
      id[nt++] = INT;

      /* note: this is used just to calculate All.MaxPartMet in
       * read_file(); then, it should give a reasonable "average"
       * of the generations of each specified IMFs.
       *
       * WARNING: this will be overridden by read_imfs() !!!
       */
      strcpy(tag[nt], "Generations");
      addr[nt] = &All.Generations;
      id[nt++] = INT;

      /* the baryon fraction that you expect to end in stars
       */
      strcpy(tag[nt], "SFfactor");
      addr[nt] = &All.SFfactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMFsFileName");
      addr[nt] = All.IMFfilename;
      id[nt++] = STRING;


      strcpy(tag[nt], "SFsFileName");
      addr[nt] = All.SFfilename;
      id[nt++] = STRING;

      /* the minimum mass for CC supernovae (usually 8 Msun) */
      strcpy(tag[nt], "SnII_InfMass");
      addr[nt] = &All.Mup;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinBinSystemMass");
      addr[nt] = &All.MBm;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxBinSystemMass");
      addr[nt] = &All.MBM;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinarySystemFrac");
      addr[nt] = &All.BinFrac;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinStarMass_inBinSystem");
      addr[nt] = &All.MBms;
      id[nt++] = REAL;

#ifdef LT_AVOID_ENRICH_SFGAS
      /* this value is set in Msun/yr ;
       * Gas Particles having a sfr larger that this
       * value will not receive metals from stars.
       */
      strcpy(tag[nt], "SFTh_for_enrich");
      addr[nt] = &All.Enrich_SFGas_Th;
      id[nt++] = REAL;
#endif

#ifdef LT_MOD_EFFM
      strcpy(tag[nt], "ModSEffCrit");
      addr[nt] = &All.Mod_SEff;
      id[nt++] = INT;
#endif

#ifdef LT_HOT_EJECTA
      strcpy(tag[nt], "EgySpecEjecta");
      addr[nt] = &All.EgySpecEjecta;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "SnIaDataFileName");
      addr[nt] = All.SnIaDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnIaDataSetNum");
      addr[nt] = &All.Ia_Nset_ofYields;
      id[nt++] = INT;

      strcpy(tag[nt], "SnIaEnergy");
      addr[nt] = &All.SnIaEgy;
      id[nt++] = REAL;

      strcpy(tag[nt], "LongLiving_Step_Prec");
      addr[nt] = &All.LLv_Step_Prec;
      id[nt++] = REAL;

      /*
         strcpy(tag[nt], "SnIa_Remnant");
         addr[nt] = &All.SnIaRemn;
         id[nt++] = REAL;
       */

      strcpy(tag[nt], "SnIIdataFileName");
      addr[nt] = All.SnIIDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnIIDataSetNum");
      addr[nt] = &All.II_Nset_ofYields;
      id[nt++] = INT;

      strcpy(tag[nt], "SnIIEnergy");
      addr[nt] = &All.SnIIEgy;
      id[nt++] = REAL;

      /*
         strcpy(tag[nt], "ChemTimeStepII");
         addr[nt] = &All.ChemTimeStepII;
         id[nt++] = REAL;

         strcpy(tag[nt], "LongCTStepII");
         addr[nt] = &All.LongChemTimeStepII;
         id[nt++] = REAL;
       */

      strcpy(tag[nt], "metIRAThMass");
      addr[nt] = &All.metIRA_ThMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "egyIRAThMass");
      addr[nt] = &All.egyIRA_ThMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "SnII_Step_Prec");
      addr[nt] = &All.SnII_Step_Prec;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBdataFileName");
      addr[nt] = All.AGBDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "AGBDataSetNum");
      addr[nt] = &All.AGB_Nset_ofYields;
      id[nt++] = INT;

      strcpy(tag[nt], "MinChemTimeStep");
      addr[nt] = &All.MinChemTimeStep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSpreadLength");
      addr[nt] = &All.MinChemSpreadL;
      id[nt++] = REAL;

      strcpy(tag[nt], "Z_toset_SF_DensTh");
      addr[nt] = &All.referenceZ_toset_SF_DensTh;
      id[nt++] = REAL;

      strcpy(tag[nt], "Zdependent_SFTh");
      addr[nt] = &All.SFTh_Zdep;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxSpreadLength");
      addr[nt] = &All.MaxChemSpreadL;
      id[nt++] = REAL;

#endif

#ifdef LT_METAL_COOLING_WAL
      strcpy(tag[nt], "WalCoolTables_Path");
      addr[nt] = All.WalCool_CoolTables_path;
      id[nt++] = STRING;
#endif

#ifdef LT_ADD_GAL_TO_SUB
      strcpy(tag[nt], "BC_SED_Path");
      addr[nt] = All.BC_SED_Path;
      id[nt++] = STRING;
#endif

#ifdef CHEMCOOL
      strcpy(tag[nt], "H2RefDustEff");
      addr[nt] = &All.H2RefDustEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "OxyAbund");
      addr[nt] = &All.OxyAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "CarbAbund");
      addr[nt] = &All.CarbAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SiAbund");
      addr[nt] = &All.SiAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "DeutAbund");
      addr[nt] = &All.DeutAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "MgAbund");
      addr[nt] = &All.MgAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "UVField");
      addr[nt] = &All.UVField;
      id[nt++] = REAL;

      strcpy(tag[nt], "PhiPAH");
      addr[nt] = &All.PhiPAH;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitDustTemp");
      addr[nt] = &All.InitDustTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "DustToGasRatio");
      addr[nt] = &All.DustToGasRatio;
      id[nt++] = REAL;

      strcpy(tag[nt], "AVConversionFactor");
      addr[nt] = &All.AVConversionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "CosmicRayIonRate");
      addr[nt] = &All.CosmicRayIonRate;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitRedshift");
      addr[nt] = &All.InitRedshift;
      id[nt++] = REAL;

      strcpy(tag[nt], "ExternalDustExtinction");
      addr[nt] = &All.ExternalDustExtinction;
      id[nt++] = REAL;

      strcpy(tag[nt], "H2FormEx");
      addr[nt] = &All.H2FormEx;
      id[nt++] = REAL;

      strcpy(tag[nt], "H2FormKin");
      addr[nt] = &All.H2FormKin;
      id[nt++] = REAL;

      strcpy(tag[nt], "PhotochemApprox");
      addr[nt] = &All.PhotochemApprox;
      id[nt++] = INT;

      strcpy(tag[nt], "ADRateFlag");
      addr[nt] = &All.ADRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "MNRateFlag");
      addr[nt] = &All.MNRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "AtomicFlag");
      addr[nt] = &All.AtomicFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagA");
      addr[nt] = &All.ThreeBodyFlagA;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagB");
      addr[nt] = &All.ThreeBodyFlagB;
      id[nt++] = INT;

      strcpy(tag[nt], "H3PlusRateFlag");
      addr[nt] = &All.H3PlusRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "InitMolHydroAbund");
      addr[nt] = &All.InitMolHydroAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitHPlusAbund");
      addr[nt] = &All.InitHPlusAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitDIIAbund");
      addr[nt] = &All.InitDIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitHDAbund");
      addr[nt] = &All.InitHDAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitHeIIAbund");
      addr[nt] = &All.InitHeIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitHeIIIAbund");
      addr[nt] = &All.InitHeIIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitCIIAbund");
      addr[nt] = &All.InitCIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitSiIIAbund");
      addr[nt] = &All.InitSiIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitOIIAbund");
      addr[nt] = &All.InitOIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitCOAbund");
      addr[nt] = &All.InitCOAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitC2Abund");
      addr[nt] = &All.InitC2Abund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitOHAbund");
      addr[nt] = &All.InitOHAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitH2OAbund");
      addr[nt] = &All.InitH2OAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitO2Abund");
      addr[nt] = &All.InitO2Abund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitHCOPlusAbund");
      addr[nt] = &All.InitHCOPlusAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitCHAbund");
      addr[nt] = &All.InitCHAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitCH2Abund");
      addr[nt] = &All.InitCH2Abund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitSiIIIAbund");
      addr[nt] = &All.InitSiIIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitCH3PlusAbund");
      addr[nt] = &All.InitCH3PlusAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitMgPlusAbund");
      addr[nt] = &All.InitMgPlusAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadHeatFlag");
      addr[nt] = &All.RadHeatFlag;
      id[nt++] = INT;
#endif

#ifdef GENERATE_GAS_IN_ICS
#ifdef GENERATE_GAS_TG
      strcpy(tag[nt], "GenGasRefFac");
      addr[nt] = &All.GenGasRefFac;
      id[nt++] = INT;
#endif
#endif

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      printf("Obtaining parameters from file '%s':\n", fname);
	      while(!feof(fd))
		{

		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  fprintf(stdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  fprintf(stdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  fprintf(stdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
#ifdef ALLOWEXTRAPARAMS
		      fprintf(stdout, "WARNING from file %s:   Tag '%s' ignored !\n", fname, buf1);
#else
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
#endif
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      printf("\n");

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

#ifdef INVARIANCETEST
	      i = strlen(All.OutputDir);
	      All.OutputDir[i - 1] = 0;
	      strcat(All.OutputDir, "/run0");
	      i = strlen(All.OutputDir);
	      All.OutputDir[i - 1] += Color;
	      mkdir(All.OutputDir, 02755);
	      strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", All.OutputDir, "logfile.txt");
	      printf("stdout will now appear in the file '%s'\n", buf1);
	      fflush(stdout);
	      freopen(buf1, "w", stdout);
#endif

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}

      /* Counts number of CR populations in parameter file */
      /*
         #ifdef COSMIC_RAYS
         CRpop=0;
         while ((All.CR_Alpha[CRpop] != 0.0 ) && ( CRpop < MaxNumCRpop ))
         CRpop++;
         NUMCRPOP = CRpop;
         #else
         NUMCRPOP = 1;
         #endif
       */
#ifdef COSMIC_RAYS
      printf(" NUMCRPOP = %i \n", NUMCRPOP);
#endif

#ifdef BP_REAL_CRs
      printf(" BP_REAL_CRs = %i \n", BP_REAL_CRs);
#endif

      for(i = 0; i < nt; i++)
	{
#ifdef LT_STELLAREVOLUTION
	  if(*tag[i] &&
	     (strcmp(tag[i], "metIRAThMass") != 0 && strcmp(tag[i], "egyIRAThMass") != 0) &&
	     strcmp(tag[i], "WindEfficiency") != 0 && strcmp(tag[i], "WindEnergyFraction") != 0 &&
	     strcmp(tag[i], "LocalSpreadFactor") != 0 && strcmp(tag[i], "TestSuite"))
#else
	  if(*tag[i])
#endif
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	{
	  errorFlag += read_outputlist(All.OutputListFilename);
	}
      else
	{
	  All.OutputListLength = 0;
	}
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#if defined(VORONOI) && !defined(PERIODIC)
  if(ThisTask == 0)
    {
      printf("Voronoi code requires periodic boundaries.\n");
    }
  endrun(0);
#endif

#ifdef EDDINGTON_TENSOR_BH
#ifndef BLACK_HOLES
  if(ThisTask == 0)
    {
      printf("Code was compiled with EDDINGTON_TENSOR_BH, but not with BLACK_HOLES.\n");
      printf("Switch on BLACK_HOLES.\n");
    }
  endrun(0);
#endif
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef SFR

#ifndef MOREPARAMS
  if(ThisTask == 0)
    {
      printf("Code was compiled with SFR, but not with MOREPARAMS.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif



#ifdef METALS
#ifndef SFR
  if(ThisTask == 0)
    {
      printf("Code was compiled with METALS, but not with SFR.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MOREPARAMS
#ifdef TIME_DEP_ART_VISC
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_ART_VISC, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DARKENERGY, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef TIMEDEPDE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef EULERPOTENTIALS
#ifdef EULER_DISSIPATION
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with EULER_DISSIPATION, but not with EULERPOTENTIALS.\n");
      fprintf(stdout, "This is not sensible.\n");
    }
  endrun(0);
#endif
#endif
#if defined(EULER_DISSIPATION) && defined(MAGNETIC_DISSIPATION)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with EULER_DISSIPATION, and MAGNETIC_DISSIPATION.\n");
      fprintf(stdout, "This is not sensible.\n");
    }
  endrun(0);
#endif

#ifndef MAGNETIC
#ifdef TRACEDIVB
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TRACEDIVB, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DBOUTPUT
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DBOUTPUT, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MAGFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGFORCE, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BSMOOTH
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BSMOOTH, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BFROMROTA
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BFROMROTA, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MU0_UNITY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MU0_UNITY, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef MAGNETIC_DISSIPATION
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGNETIC_DISSIPATION, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef TIME_DEP_MAGN_DISP
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_MAGN_DISP, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef DIVBCLEANING_DEDNER
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBCLEANING_DEDNER, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#endif

#if defined(NOWINDTIMESTEPPING) && defined(MAGNETIC)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with NOWINDTIMESTEPPING and with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed, as it leads to inconsitent MHD for wind particles.\n");
    }
  endrun(0);
#endif

#ifndef MAGFORCE
#if defined(DIVBFORCE) || defined(DIVBFORCE3)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#if defined(DIVBFORCE) && defined(DIVBFORCE3)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE and with DIVBFORCE3.\n");
      fprintf(stdout, "This will lead to no correction at all, better stop.\n");
    }
  endrun(0);
#endif

#ifdef BH_BUBBLES
#ifndef BLACK_HOLES
  if(ThisTask == 0)
    {
      printf("Code was compiled with BH_BUBBLES, but not with BLACK_HOLES.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES) || defined(EBUB_PROPTO_BHAR)
  if(ThisTask == 0)
    {
      printf
	("If the code is compiled with BH_BUBBLES, then BUBBLES, MULTI_BUBBLES or EBUB_PROPTO_BHAR options cannot be used.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#if defined(LT_STELLAREVOLUTION)

#if !defined(LT_PM_LIFETIMES) && !defined(LT_MM_LIFETIMES)
  if(ThisTask == 0)
    printf("LT_STELLAREVOLUTION requires either LT_PM_LIFETIMES or LT_MM_LIFETIMES\n");
  endrun(0);
#endif

#if (defined(LT_WIND_VELOCITY) || defined(LT_HYDROWINDS) || defined(LT_DECOUPLE_POSTWINDS_FROM_SF)) && !defined(WINDS)
  if(ThisTask == 0)
    printf("LT_WIND_VELOCITY, LT_HYDROWINDS and LT_DECOUPLE_POSTWINDS_FROM_SF require WIND\n");
  endrun(0);
#endif


#endif

#if defined(LT_METAL_COOLING)

#if defined(LT_METAL_COOLING_on_SMOOTH_Z) && !defined(LT_SMOOTH_Z)
  if(ThisTask == 0)
    printf("LT_METAL_COOLING_on_SMOOTH_Z requires LT_SMOOTH_Z\n ");
  endrun(0);
#endif

#if defined(LT_METAL_COOLING_on_SMOOTH_Z) && defined(LT_METAL_COOLING_WAL) && !defined(LT_SMOOTH_ALLMETALS)
  if(ThisTask == 0)
    printf("LT_METAL_COOLING_on_SMOOTH_Z requires LT_SMOOTH_ALLMETALS if LT_METAL_COOLING_WAL is defined\n ");
  endrun(0);
#endif

#endif

#ifdef KD_SEED_STAR_MASS_FRACTION
#ifndef KD_CHOOSE_LINKING_LENGTH
  if(ThisTask == 0)
    printf("KD_SEED_STAR_MASS_FRACTION requires KD_CHOOSE_LINKING_LENGTH to be defined\n ");
  endrun(0);
#endif
#endif

#ifdef NUM_THREADS
#ifdef _OPENMP
  if(ThisTask == 0)
    printf("NUM_THREADS is incompatible with enabeling OpenMP in the compiler options \n");
  endrun(0);      
#endif
#endif


#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS


#ifdef COSMIC_RAYS
  if(ThisTask == 0)
    {
      printf("CR SN Efficiency: %g\n", All.CR_SNEff);
    }
#endif
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
	break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
	flag = 1;

      if(count == 1 || count == 2)
	{
	  if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
	    {
	      if(ThisTask == 0)
		printf("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n",
		       (int) MAXLEN_OUTPUTLIST);
	      endrun(13);
	    }

	  All.OutputListFlag[All.OutputListLength] = flag;
	  All.OutputListLength++;
	}
    }

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif
#ifdef CONDUCTION
      All.Conduction_Ti_begstep /= 2;
      All.Conduction_Ti_endstep /= 2;
#endif
#ifdef CR_DIFFUSION
      All.CR_Diffusion_Ti_begstep /= 2;
      All.CR_Diffusion_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_current /= 2;

	  if(P[i].TimeBin > 0)
	    {
	      P[i].TimeBin--;
	      if(P[i].TimeBin <= 0)
		{
		  printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
		  endrun(8765);
		}
	    }
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
