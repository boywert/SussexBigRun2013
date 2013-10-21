#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef BP_REAL_CRs
#include "bp_cosmic_rays/bp_cosmic_rays.h"
#endif
#ifdef CS_MODEL
#include "cs_metals.h"
#endif

#ifdef VORONOI
#include "voronoi.h"
#endif

#ifdef MACHNUM
#ifdef COSMIC_RAYS
#define h  All.HubbleParam
#define cm (h/All.UnitLength_in_cm)
#define s  (h/All.UnitTime_in_s)
#define LightSpeed (2.9979e10*cm/s)
#define c2   ( LightSpeed * LightSpeed )
#endif
#endif



/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3, atime;

#ifdef MAGNETIC
  double a2_fac, gauss2gadget = 1.0;
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#ifdef BP_REAL_CRs
  double Pinit_integral;
  double Pth0;
#endif

#if defined(COSMIC_RAYS) && defined(MACHNUM)
  double Pth1, PCR1[NUMCRPOP], rBeta[NUMCRPOP], C_phys[NUMCRPOP], q_phys[NUMCRPOP];
#endif
#ifdef CR_INITPRESSURE
  double cr_pressure, q_phys, C_phys[NUMCRPOP];
#endif
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif  
#ifdef BLACK_HOLES
  int count_holes = 0;
#endif

#ifdef EULERPOTENTIALS
  double a0, a1, a2;
  double b0, b1, b2;
#endif

#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationMerk;
#endif

#ifdef DISTORTIONTENSORPS
  int i1, i2;
#endif

#ifdef LT_STELLAREVOLUTION
  double this_age, NextChemTime;
  unsigned int IMFi, SFi, I;
  int chem_step, ti_min;
#endif

#ifdef LIMIT_HSML
  double meansep;
#endif
  
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

  if(RestartFlag == 3 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
      endrun(0);
    }

  if(RestartFlag == 4 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if snapshot should be converted\n");
      endrun(0);
    }

  if(RestartFlag == 5 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf
	  ("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");
      endrun(0);
    }


  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
    case 4:
      if(RestartFlag >= 2 && RestartSnapNum >= 0)
	{
	  char fname[1000];

	  if(All.NumFilesPerSnapshot > 1)
	    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
		    RestartSnapNum);
	  else
	    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
	  read_ic(fname);

	}
      else
	{
	  if(All.ICFormat == 4)
	    read_ic_cluster(All.InitCondFile);
	  else
	    read_ic(All.InitCondFile);
	}
      break;

    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();


#ifdef WINDTUNNEL
  All.WindCurrentX = -0.5 * All.WindDmean;
#endif


#ifdef MOL_CLOUDS
 if(ThisTask == 0)
  {
   printf("Init MOL_CLOUDS...\n");
   fflush(stdout);
  }
 init_mol_clouds();
 if(ThisTask == 0)
  {
   printf("done.\n");
   fflush(stdout);
  }
#endif


#ifdef SCFPOTENTIAL
 if(ThisTask == 0)
  {
   printf("Init SCF...\n");
   fflush(stdout);
  }
 SCF_init();
 if(ThisTask == 0)
  {
   printf("Initial random seed = %ld\n", scf_seed);
   printf("done.\n");
   fflush(stdout);
  }
#endif




#ifdef COOLING
#ifdef CS_MODEL
  if(RestartFlag == 0)
    XH = HYDROGEN_MASSFRAC;
#endif
#ifndef LT_METAL_COOLING_WAL  
  IonizeParams();
#endif
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
  InitChem();
#endif

#ifdef BP_REAL_CRs
  init_cr_bound();
  Pinit_integral = pressure_integral(All.ecr_bound[0], All.ecr_bound[BP_REAL_CRs], All.pSlope_init, CNST_MP);
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gauss2gadget *= sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / (All.HubbleParam * All.HubbleParam));
#endif
      a2_fac = (All.Time * All.Time);
#endif
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gauss2gadget *= sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g);
#endif
      a2_fac = 1;
#endif
    }

#ifdef RADTRANSFER
  All.Radiation_Ti_begstep = 0;
#endif


  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;

  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
	{
	  char *underscore = strrchr(All.InitCondFile,'_');
	  if (!underscore) 
	    {
	      char buf[1000];
	      sprintf(buf,"Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n",
		      All.InitCondFile);
	      terminate(buf);
	    } 
	  else 
	    {
	      All.SnapshotFileCount = atoi(underscore+1) + 1;
	    }
	}
      else
	{
	  All.SnapshotFileCount = RestartSnapNum + 1;
	}
    }

#ifdef OUTPUTLINEOFSIGHT
  All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
  if(RestartFlag == 2)
    endrun(78787);
#endif

  All.TotNumOfForces = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
#ifdef SETMAINTIMESTEPCOUNT
  All.MainTimestepCounts = All.MainTimestepCountIni;
#else
  All.MainTimestepCounts = 0;
#endif
#endif

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.7;


#ifdef VORONOI
  Indi.AllocFacNdp = 1.2 * N_gas + MIN_ALLOC_NUMBER;
  Indi.AllocFacNdt = 8.0 * N_gas + MIN_ALLOC_NUMBER;
  Indi.AllocFacNvf = 8.0 * N_gas + MIN_ALLOC_NUMBER;
  Indi.AllocFacNinlist = 1.2 * N_gas + MIN_ALLOC_NUMBER;
  Indi.AllocFacN_DP_Buffer = 0.2 * N_gas + MIN_ALLOC_NUMBER;;
#endif

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
  All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

  for(i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelToTimeBin[i] = 0;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < GRAVCOSTLEVELS; j++)
      P[i].GravCost[j] = 0;

#ifdef BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and first bubble time %g \n", All.TimeBegin, All.TimeOfNextBubble);

  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

#ifdef MULTI_BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and time of the first bubbles %g \n", All.TimeBegin, All.TimeOfNextBubble);
  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }


/* DISTORTION GENERAL SETUP */
#ifdef DISTORTIONTENSORPS
#if defined(SIM_ADAPTIVE_SOFT) || defined(REINIT_AT_TURNAROUND)
#ifdef REINIT_AT_TURNAROUND
  All.SIM_epsilon = REINIT_AT_TURNAROUND;
#endif
#ifdef SIM_ADAPTIVE_SOFT
  All.SIM_epsilon = SIM_ADAPTIVE_SOFT;
#endif
  All.cms_x = 0.0;
  All.cms_y = 0.0;
  All.cms_z = 0.0;
#ifdef COMOVING_DISTORTION
  All.InitialTurnaroundRadius = COMOVING_DISTORTION * pow(All.TimeBegin, 1.0 / (3.0 * All.SIM_epsilon));
  if(ThisTask == 0)
    {
      printf("REINIT_AT_TURNAROUND: comoving initial/final turnaround radius = %g/%g\n",
	     All.InitialTurnaroundRadius, COMOVING_DISTORTION);
      fflush(stdout);
    }
#else
  double turnaround_radius_local = 0.0, turnaround_radius_global = 0.0, v_part, r_part;

  /* get local turnaroundradius */
  for(i = 0; i < NumPart; i++)
    {
      r_part = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
      v_part = (P[i].Pos[0] * P[i].Vel[0] + P[i].Pos[1] * P[i].Vel[1] + P[i].Pos[2] * P[i].Vel[2]);
      if((v_part < 0.0) && (r_part > turnaround_radius_local))
	turnaround_radius_local = r_part;
    }

  /* find global turnaround radius by taking maximum of all CPUs */
  MPI_Allreduce(&turnaround_radius_local, &turnaround_radius_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  All.InitialTurnaroundRadius = turnaround_radius_global;

#endif
  All.CurrentTurnaroundRadius = All.InitialTurnaroundRadius;
  if(ThisTask == 0)
    {
      printf("REINIT_AT_TURNAROUND: current (initial) turnaround radius = %g\n", All.CurrentTurnaroundRadius);
      fflush(stdout);
    }
#endif
#endif


  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = 0;

/* DISTORTION PARTICLE SETUP */
#ifdef DISTORTIONTENSORPS
      /*init tidal tensor for first output (not used for calculation) */
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].tidal_tensorps[i1][i2] = 0.0;

      /* find caustics by sign analysis of configuration space distortion */
      P[i].last_stream_determinant = 1.0;

      /* distortion tensor alrady read from ICs -> get correct last_stream_determinant */
#ifdef DISTORTION_READALL
      MyDouble product_matrix[3][3];

      product_matrix[0][0] = P[i].distortion_tensorps[0][0] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][0];
      product_matrix[0][1] = P[i].distortion_tensorps[0][1] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][1];
      product_matrix[0][2] = P[i].distortion_tensorps[0][2] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][2];
      product_matrix[1][0] = P[i].distortion_tensorps[1][0] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][0];
      product_matrix[1][1] = P[i].distortion_tensorps[1][1] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][1];
      product_matrix[1][2] = P[i].distortion_tensorps[1][2] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][2];
      product_matrix[2][0] = P[i].distortion_tensorps[2][0] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][0];
      product_matrix[2][1] = P[i].distortion_tensorps[2][1] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][1];
      product_matrix[2][2] = P[i].distortion_tensorps[2][2] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][2];

      /* this determinant will change sign when we pass through a caustic -> criterion for caustics */
      P[i].last_stream_determinant =
	((product_matrix[0][0]) * (product_matrix[1][1]) * (product_matrix[2][2]) +
	 (product_matrix[0][1]) * (product_matrix[1][2]) * (product_matrix[2][0]) +
	 (product_matrix[0][2]) * (product_matrix[1][0]) * (product_matrix[2][1]) -
	 (product_matrix[0][2]) * (product_matrix[1][1]) * (product_matrix[2][0]) -
	 (product_matrix[0][0]) * (product_matrix[1][2]) * (product_matrix[2][1]) -
	 (product_matrix[0][1]) * (product_matrix[1][0]) * (product_matrix[2][2]));


#endif

#ifdef REINIT_AT_TURNAROUND
      /* no phase-space analysis for particles that are initially within the turnaround radius */
      if(sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]) <
	 All.CurrentTurnaroundRadius)
#ifdef DISTORTION_READALL
	P[i].turnaround_flag = 1;
#else
	P[i].turnaround_flag = -1;
#endif
      else
	P[i].turnaround_flag = 0;
#endif

#ifdef OUTPUT_LAST_CAUSTIC
      /* all entries zero -> no caustic yet */
      P[i].lc_Time = 0.0;
      P[i].lc_Pos[0] = 0.0;
      P[i].lc_Pos[1] = 0.0;
      P[i].lc_Pos[2] = 0.0;
      P[i].lc_Vel[0] = 0.0;
      P[i].lc_Vel[1] = 0.0;
      P[i].lc_Vel[2] = 0.0;
      P[i].lc_rho_normed_cutoff = 0.0;

      P[i].lc_Dir_x[0] = 0.0;
      P[i].lc_Dir_x[1] = 0.0;
      P[i].lc_Dir_x[2] = 0.0;
      P[i].lc_Dir_y[0] = 0.0;
      P[i].lc_Dir_y[1] = 0.0;
      P[i].lc_Dir_y[2] = 0.0;
      P[i].lc_Dir_z[0] = 0.0;
      P[i].lc_Dir_z[1] = 0.0;
      P[i].lc_Dir_z[2] = 0.0;

      P[i].lc_smear_x = 0.0;
      P[i].lc_smear_y = 0.0;
      P[i].lc_smear_z = 0.0;
#endif


#ifdef PMGRID
      /* long range tidal field init */
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif

#ifndef DISTORTION_READALL
      for(i1 = 0; i1 < 6; i1++)
	for(i2 = 0; i2 < 6; i2++)
	  {
	    if((i1 == i2))
	      P[i].distortion_tensorps[i1][i2] = 1.0;
	    else
	      P[i].distortion_tensorps[i1][i2] = 0.0;
	  }
#endif

#ifdef COMOVING_DISTORTION
      /* for cosmological simulations we do init here, not read from ICs */

      /* Lagrange time of particle */
      P[i].a0 = All.TimeBegin;

#ifndef COMOVING_READIC
      /* no caustic passages in the beginning */
      P[i].caustic_counter = 0.0;

      /* approximation: perfect Hubble Flow -> peculiar sheet orientation is exactly zero */
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].V_matrix[i1][i2] = 0.0;

      /* approximation: initial sream density equals background density */
      P[i].init_density = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif
#endif

      /* annihilation stuff */
      P[i].s_1_last = 1.0;
      P[i].s_2_last = 1.0;
      P[i].s_3_last = 1.0;
      P[i].second_deriv_last = 0.0;
      P[i].rho_normed_cutoff_last = 1.0;

      P[i].s_1_current = 1.0;
      P[i].s_2_current = 1.0;
      P[i].s_3_current = 1.0;
      P[i].second_deriv_current = 0.0;
      P[i].rho_normed_cutoff_current = 1.0;

      P[i].annihilation = 0.0;
      P[i].analytic_caustics = 0.0;
      P[i].analytic_annihilation = 0.0;

#ifdef COMOVING_DISTORTION
      P[i].stream_density = P[i].init_density / (All.TimeBegin * All.TimeBegin * All.TimeBegin);
#else
      P[i].stream_density = P[i].init_density;
#endif

#endif /* DISTORTIONTENSORPS */

#ifdef KEEP_DM_HSML_AS_GUESS
      if(RestartFlag != 1)
	P[i].DM_Hsml = -1;
#endif

#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;

      if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
	P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      P[i].p.Potential = 0;
#endif
#ifdef STELLARAGE
      if(RestartFlag == 0)
	P[i].StellarAge = 0;
#endif

#ifdef METALS
#ifndef CS_MODEL
      if(RestartFlag == 0)
	P[i].Metallicity = 0;
#else
      if(RestartFlag == 0)
	{
	  for(j = 0; j < 12; j++)
	    {
	      P[i].Zm[j] = 0;
	      P[i].ZmReservoir[j] = 0;
	    }

	  P[i].Zm[6] = HYDROGEN_MASSFRAC * (P[i].Mass);
	  P[i].Zm[0] = (1 - HYDROGEN_MASSFRAC) * (P[i].Mass);
	  /*     Order of chemical elements:   He, Carbon,Mg,O,Fe,Si,H,N,Ne,S,Ca,Zn */

	  PPP[i].Hsml = 0;
#ifdef CS_FEEDBACK
	  PPP[i].EnergySN = 0;
	  PPP[i].EnergySNCold = 0;
#endif
	}
#endif

#endif

#ifdef LT_STELLAREVOLUTION
      All.Time_Age = get_age(All.Time);
      
      if(P[i].Type == 4)
	{
	  this_age = get_age(P[i].StellarAge) - All.Time_Age;
	  SFi = get_SF_index(i, &SFi, &IMFi);
	  
          NextChemTime = get_NextChemTime(this_age, SFi, 0x0);
          j = get_chemstep_bin(All.Time, All.Time_Age - NextChemTime, &chem_step, i);

          I = P[i].pt.MetID;
          
          if(All.Ti_Current >= TIMEBASE)
            chem_step = j = 0;
  
          if((TIMEBASE - All.Ti_Current) < chem_step)
            {
              chem_step = TIMEBASE - All.Ti_Current;
              ti_min = TIMEBASE;
              while(ti_min > chem_step)
                ti_min >>= 1;
              chem_step = ti_min;
              j = get_timestep_bin(chem_step);
            }
          
          if(j != MetP[I].ChemTimeBin)
            {
              TimeBinCountStars[MetP[I].ChemTimeBin]--;
              TimeBinCountStars[j]++;
              MetP[I].ChemTimeBin = j;
            }
          MetP[I].LastChemTime = this_age;
	
	  if(this_age < SFs[SFi].ShortLiv_TimeTh)
	    MetP[I].LastChemTime = SFs[SFi].ShortLiv_TimeTh;
	}
#endif

      
#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  count_holes++;

	  if(RestartFlag == 0)
	    BPP(i).BH_Mass = All.SeedBlackHoleMass;
#ifdef BH_BUBBLES
	  if(RestartFlag == 0)
	    {
	      BPP(i).BH_Mass_bubbles = All.SeedBlackHoleMass;
	      BPP(i).BH_Mass_ini = All.SeedBlackHoleMass;
#ifdef UNIFIED_FEEDBACK
	      BPP(i).BH_Mass_radio = All.SeedBlackHoleMass;
#endif
	    }
#endif
	}
#endif
    }

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef CONDUCTION
  All.Conduction_Ti_endstep = All.Conduction_Ti_begstep = 0;
#endif
#ifdef CR_DIFFUSION
  All.CR_Diffusion_Ti_endstep = All.CR_Diffusion_Ti_begstep = 0;
#endif

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      SphP[i].EntropyPred = SphP[i].Entropy;

      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
	}

      SphP[i].e.DtEntropy = 0;

#ifdef EOS_DEGENERATE
      for(j = 0; j < 3; j++)
	SphP[i].xnucPred[j] = SphP[i].xnuc[j];
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif

#ifdef CS_MODEL
      SphP[i].DensityOld = 0;
#endif

#ifdef CS_FEEDBACK
      SphP[i].da.DensityAvg = 0;
      SphP[i].ea.EntropyAvg = 0;
      SphP[i].DensPromotion = 0;
      SphP[i].TempPromotion = 0;
      SphP[i].HotHsml = 0;
#ifdef CS_ATTEMPTS
      SphP[i].PromotionAttempts = 0;
#endif
#endif

      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
	  SphP[i].d.Density = -1;
#ifdef COOLING
#ifndef UM_CHEMISTRY          
	  SphP[i].Ne = 1.0;
#endif
#endif
#ifdef CHEMCOOL
	  SphP[i].TracAbund[IH2] = All.InitMolHydroAbund;
	  SphP[i].TracAbund[IHP] = All.InitHPlusAbund;
	  SphP[i].TracAbund[IDP] = All.InitDIIAbund;
	  SphP[i].TracAbund[IHD] = All.InitHDAbund;
	  SphP[i].TracAbund[IHEP] = All.InitHeIIAbund;
	  SphP[i].TracAbund[IHEPP] = All.InitHeIIIAbund;
#endif
	  SphP[i].v.DivVel = 0;
	}
#ifdef WINDS
      SphP[i].DelayTime = 0;
#ifdef VARIABLE_WINDS
      SphP[i].HostHaloMass = 0;
#endif
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
      SphP[i].elec = All.Startelec;

      SphP[i].HI = All.StartHI;
      SphP[i].HII = All.StartHII;
      SphP[i].HM = All.StartHM;

      SphP[i].HeI = All.StartHeI;
      SphP[i].HeII = All.StartHeII;
      SphP[i].HeIII = All.StartHeIII;

      SphP[i].H2I = All.StartH2I;
      SphP[i].H2II = All.StartH2II;

      SphP[i].HD = All.StartHD;
      SphP[i].DI = All.StartDI;
      SphP[i].DII = All.StartDII;

      SphP[i].HeHII = All.StartHeHII;
#endif      
#ifdef MAGNETIC
#ifdef BINISET
      if (RestartFlag == 0){ /* Set only when starting from ICs */
          SphP[i].b2.BPred[0] = All.BiniX;
          SphP[i].b2.BPred[1] = All.BiniY;
          SphP[i].b2.BPred[2] = All.BiniZ;
       }
#ifdef VECT_POTENTIAL
      if(All.BiniY == 0 && All.BiniZ == 0)
	{
	  SphP[i].APred[0] = 0.;
	  SphP[i].APred[1] = 0.;
	  SphP[i].APred[2] = (All.BiniX * P[i].Pos[1]) * atime * All.HubbleParam;
	  SphP[i].A[0] = SphP[i].APred[0];
	  SphP[i].A[1] = SphP[i].APred[1];
	  SphP[i].A[2] = SphP[i].APred[2];
	  SphP[i].DtA[0] = SphP[i].DtA[1] = SphP[i].DtA[2] = 0.;
	}
      else
	{
	  if(ThisTask == 0)
	    {
	      printf("Can not reconstruct VECTOR potentials from Bini values !\n");
	      endrun(6724);
	    }
	}
#endif
#ifdef EULERPOTENTIALS
      if(All.BiniY == 0 && All.BiniZ == 0)
	{
	  a0 = 0;
	  a1 = sqrt(All.BiniX);
	  a2 = 0;
	  b0 = 0;
	  b1 = sqrt(All.BiniX);
	  b2 = sqrt(All.BiniX);
	}
      else
	{
	  if(All.BiniX != 0 && (All.BiniY != 0 || All.BiniZ != 0))
	    {
	      b0 = -(All.BiniZ + All.BiniY) / All.BiniX;
	      b1 = 1;
	      b2 = 1;
	      a0 = 0;
	      a1 = -sqrt(All.BiniZ) / b0;
	      a2 = sqrt(All.BiniY) / b0;
#ifdef PERIODIC
	      if(ThisTask ==0)
		printf("Initializing EulerPotentials from current Bini settings not yet supported for PERIODIC\n");
              endrun(6519654);
#endif
	    }
	  else
	    {
	      a0 = a1 = a2 = b0 = b1 = b2 = 0;
	      if(ThisTask == 0)
		{
		  printf("Can not reconstruct Euler potentials from Bini values !\n");
		  endrun(6723);
		}
	    }
	}
#ifdef PERIODIC
      /* spread a periode of sin/cos over 10 particles in linear size */
      /* the /fac is needed as euler potentials are gradien in respect to x but we use cos(fac*x) */
      double fac = 2 * M_PI / All.BoxSize * pow(All.TotN_gas,1./3.) / 10; 
      SphP[i].EulerA = (a0 * cos(fac * P[i].Pos[0]) + 
	                a1 * cos(fac * P[i].Pos[1]) + 
                        a2 * cos(fac * P[i].Pos[2])) / fac * sqrt(3);
      SphP[i].EulerB = (b0 * sin(fac * P[i].Pos[0]) + 
 	                b1 * sin(fac * P[i].Pos[1]) + 
                        b2 * sin(fac * P[i].Pos[2])) / fac * sqrt(3);
#else
      SphP[i].EulerA = a0 * P[i].Pos[0] + a1 * P[i].Pos[1] + a2 * P[i].Pos[2];
      SphP[i].EulerB = b0 * P[i].Pos[0] + b1 * P[i].Pos[1] + b2 * P[i].Pos[2];
#endif
      /* a2_fac not needed as this comes in autoatically for bfld from the gradient of a and b */
      SphP[i].EulerA *= sqrt(gauss2gadget) * sqrt(a2_fac);
      SphP[i].EulerB *= sqrt(gauss2gadget) * sqrt(a2_fac);
#endif
#endif /*BINISET*/
#ifndef EULERPOTENTIALS
	for(j = 0; j < 3; j++)
	{
	  SphP[i].b2.BPred[j] *= a2_fac * gauss2gadget;

#ifndef VECT_POTENTIAL
	  SphP[i].DtB[j] = 0;
#endif
	  SphP[i].b1.B[j] = SphP[i].b2.BPred[j];
	}
#endif
#ifdef TIME_DEP_MAGN_DISP
#ifdef HIGH_MAGN_DISP_START
      SphP[i].Balpha = All.ArtMagDispConst;
#else
      SphP[i].Balpha = All.ArtMagDispMin;
#endif
      SphP[i].DtBalpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      if(HIGH_ART_VISC_START == 0)
	SphP[i].alpha = All.ArtBulkViscConst;
      if(HIGH_ART_VISC_START > 0)
	if(P[i].Pos[0] > HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
      if(HIGH_ART_VISC_START < 0)
	if(P[i].Pos[0] < -HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#endif

#ifdef VORONOI_TIME_DEP_ART_VISC
      SphP[i].Dtalpha = 0.0;
      SphP[i].alpha = All.ArtBulkViscConst / 128.0;
#ifdef VORONOI_RELAX_VIA_VISC
      SphP[i].alpha = All.ArtBulkViscConst;
#endif
#endif
#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
      SphP[i].i.Injected_BH_Energy = 0;
#endif

    }


#ifdef MODGRAV 
  for(i = 0; i < NumPart; i++)	/* initialize ModGrav properties */
    {
      for(j = 0; j < 3; j++)
        {
          P[i].mg_grad_phi[j] = 0;
          P[i].mg_grad_phi_tree[j] = 0;
          P[i].mg_grad_phi_PM[j] = 0;
          P[i].ModGravAccel[j] = 0;
        }
      
      if(RestartFlag == 0)
	{
	  P[i].Hsml = 0.0;
	  P[i].mg_Rho_0 = 0.0;
	  P[i].mg_Rho_1 = 0.0;
	  P[i].mg_phi = 0.0;
          P[i].mg_phi_tree = 0.0;
          P[i].mg_phi_PM =0.0;
	  P[i].mg_M_phi = 0.0;

	  if(P[i].Type == 1)        /* dark matter scalar coupling */
	    P[i].mg_beta = 0.0;     /* just a dummy, to be replaced with the actual coupling function */
	  else                      /* baryonic scalar coupling */
	    P[i].mg_beta = 0.0;     /* just a dummy, to be replaced with the actual coupling function */
	}
      else
	{
	  /* read-in Phi and M_phi from restart file */
	}
    }
#endif /* MODGRAV */
  
#ifdef TWODIMS
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

#ifdef ONEDIM
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[1] = P[i].Pos[2] = 0;
      P[i].Vel[1] = P[i].Vel[2] = 0;

      P[i].g.GravAccel[1] = P[i].g.GravAccel[2] = 0;

      if(P[i].Type == 0)
        {
          SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
          SphP[i].a.HydroAccel[1] =SphP[i].a.HydroAccel[2] = 0;
        }
    }
#endif

#ifdef ASSIGN_NEW_IDS
  assign_unique_ids();
#endif

#ifndef NOTEST_FOR_IDUNIQUENESS
  test_id_uniqueness();
#endif


  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;

  domain_Decomposition(0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

#ifdef LIMIT_HSML
  //meansep = All.BoxSize / pow(All.TotNumPart, 1.0 / 3.0);
  meansep = All.BoxSize / pow(All.TotN_gas, 1.0 / 3.0);

  MaxGasHsmlFractional = 2.0;

  //All.MaxGasHsml = All.MaxGasHsmlFractional * meansep;
  MaxGasHsml = MaxGasHsmlFractional * meansep;     

  if(ThisTask == 0)
    printf("Mean gas interparticle separation = %g, MaxGasHsml = %g\n", meansep, MaxGasHsml);
  //printf("Mean gas interparticle separation = %g, MaxGasHsml = %g\n", meansep, All.MaxGasHsml);
#endif

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

#ifdef START_WITH_EXTRA_NGBDEV
  MaxNumNgbDeviationMerk = All.MaxNumNgbDeviation;
  All.MaxNumNgbDeviation = All.MaxNumNgbDeviationStart;
#endif

  if(RestartFlag != 3 && RestartFlag != 5)
    setup_smoothinglengths();

#ifdef VORONOI
  for(i = 0; i < N_gas; i++)
    SphP[i].MaxDelaunayRadius = SphP[i].Hsml;

  voronoi_mesh();
  voronoi_setup_exchange();

#ifdef MESHRELAX_DENSITY_IN_INPUT
  if(RestartFlag == 0)
    for(i = 0; i < N_gas; i++)
      P[i].Mass *= SphP[i].Volume;
#endif

#ifdef VORONOI_MESHRELAX
  double mass, masstot, egy, egytot;

  for(i = 0, mass = 0, egy = 0; i < N_gas; i++)
    {
      mass += P[i].Mass;
      egy += P[i].Mass * SphP[i].Entropy;
      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Entropy * P[i].Mass / SphP[i].Volume;
    }
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&egy, &egytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.MeanMass = masstot / All.TotN_gas;
#ifdef TWODIMS
  All.MeanPressure = GAMMA_MINUS1 * egytot / (boxSize_X * boxSize_Y);
#else
  All.MeanPressure = GAMMA_MINUS1 * egytot / (boxSize_X * boxSize_Y * boxSize_Z);
#endif
#endif

  voronoi_density();

  myfree(List_P);
  myfree(ListExports);

  myfree(DT);
  myfree(DP - 5);
  myfree(VF);
#endif

#ifdef START_WITH_EXTRA_NGBDEV
  All.MaxNumNgbDeviation = MaxNumNgbDeviationMerk;
#endif

#if defined(HEALPIX)

  compute_global_quantities_of_system();

  //this should be readed in the parameterfile
  All.Nside = 32;
  //
  if(ThisTask == 0)
    printf(" First calculation of Healpix %i with %i \n", All.Nside, NSIDE2NPIX(All.Nside));
  // initialize the healpix array (just in case)
  All.healpixmap = (float *) malloc(NSIDE2NPIX(All.Nside) * sizeof(float));
  for(i = 0; i < NSIDE2NPIX(All.Nside); i++)
    All.healpixmap[i] = 0;

  double Minmass, Maxmass;

  All.Maxmass = 0.0;
  Minmass = 1E10;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
	{
	  if(P[i].Mass != 0.0)
	    Minmass = P[i].Mass < Minmass ? P[i].Mass : Minmass;
	  Maxmass = P[i].Mass > Maxmass ? P[i].Mass : Maxmass;
	}
    }
  All.Minmass = Maxmass;

  if(ThisTask == 0)
    printf("Pasamos masas: %f // %f \n", Maxmass, Minmass);

  MPI_Allreduce(&Maxmass, &All.Maxmass, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&Minmass, &All.Minmass, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Pasamos masas x 2: %f // %f \n", All.Maxmass, All.Minmass);

  healpix_halo_cond(All.healpixmap);

#endif /* final of general HEALPIX */

  /* at this point, the entropy variable actually contains the
   * internal energy, read in from the initial conditions file.
   * Once the density has been computed, we can convert to entropy.
   */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{
#ifndef EOS_DEGENERATE

#if !defined(VORONOI_MESHRELAX) && !defined(TRADITIONAL_SPH_FORMULATION) && !defined(ISOTHERM_EQS)

	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");

	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].d.Density / a3, GAMMA_MINUS1);
	  SphP[i].EntropyPred = SphP[i].Entropy;
#endif
#else
	  for(j = 0; j < EOS_NSPECIES; j++)
	    {
	      SphP[i].dxnuc[j] = 0;
	    }

	  SphP[i].Entropy *= All.UnitEnergy_in_cgs;
	  SphP[i].EntropyPred *= All.UnitEnergy_in_cgs;
	  /* call eos with physical units, energy and entropy are always stored in physical units */
	  SphP[i].temp = -1.0;
	  
	  struct eos_result res;
	  eos_calc_egiven(SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnuc, SphP[i].EntropyPred, &SphP[i].temp, &res);
	  SphP[i].Pressure = res.p.v / All.UnitPressure_in_cgs;
	  // Warning: dpdr is in physical units ...
	  SphP[i].dpdr = res.p.drho + res.temp * gsl_pow_2(res.p.dtemp / (SphP[i].d.Density * All.UnitDensity_in_cgs)) / res.e.dtemp;
#endif
	}

      SphP[i].e.DtEntropy = 0;

      SphP[i].v.DivVel = 0;

#ifdef BP_REAL_CRs
      SphP[i].DensityOld = SphP[i].d.Density;
      Pth0 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      init_bp_cr_pressure(i, P[i].Pos[0], Pth0);
      init_bp_cr_arrays(i, Pinit_integral);
#endif

#ifdef MACHNUM
      SphP[i].Shock_MachNumber = 1.0;
#ifdef COSMIC_RAYS
      Pth1 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);

#ifdef CR_IC_PHYSICAL
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] = SphP[i].CR_C0[CRpop];
	  q_phys[CRpop] = SphP[i].CR_q0[CRpop];
	}
#else
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] = SphP[i].CR_C0[CRpop] * pow(SphP[i].d.Density, (All.CR_Alpha[CRpop] - 1.0) / 3.0);
	  q_phys[CRpop] = SphP[i].CR_q0[CRpop] * pow(SphP[i].d.Density, 1.0 / 3.0);
	}
#endif
      SphP[i].PreShock_XCR = 0.0;
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  rBeta[CRpop] = gsl_sf_beta((All.CR_Alpha[CRpop] - 2.0) * 0.5, (3.0 - All.CR_Alpha[CRpop]) * 0.5) *
	    gsl_sf_beta_inc((All.CR_Alpha[CRpop] - 2.0) * 0.5, (3.0 - All.CR_Alpha[CRpop]) * 0.5,
			    1.0 / (1.0 + q_phys[CRpop] * q_phys[CRpop]));

	  PCR1[CRpop] = C_phys[CRpop] * c2 * SphP[i].d.Density * rBeta[CRpop] / 6.0;
	  PCR1[CRpop] *= pow(atime, -3.0 * GAMMA);
	  SphP[i].PreShock_XCR += PCR1[CRpop] / Pth1;
	}

      SphP[i].PreShock_PhysicalDensity = SphP[i].d.Density / a3;
      SphP[i].PreShock_PhysicalEnergy =
	SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density / a3, GAMMA_MINUS1);

      SphP[i].Shock_DensityJump = 1.0001;
      SphP[i].Shock_EnergyJump = 1.0;
#endif /* COSMIC_RAYS */
#ifdef OUTPUT_PRESHOCK_CSND
      Pth1 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      SphP[i].PreShock_PhysicalSoundSpeed =
	sqrt(GAMMA * Pth1 / SphP[i].d.Density) * pow(atime, -3. / 2. * GAMMA_MINUS1);
      SphP[i].PreShock_PhysicalDensity = SphP[i].d.Density / a3;
#endif /* OUTPUT_PRESHOCK_CSND */
#endif /* MACHNUM */

#ifdef REIONIZATION
      All.not_yet_reionized = 1;
#endif

#ifdef CR_IC_PHYSICAL
      /* Scale CR variables so that values from IC file are now the
       * physical values, not the adiabatic invariants
       */
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  SphP[i].CR_C0[CRpop] *= pow(SphP[i].d.Density, (1.0 - All.CR_Alpha[CRpop]) / 3.0);
	  SphP[i].CR_q0[CRpop] *= pow(SphP[i].d.Density, -1.0 / 3.0);
	}
#endif

#ifdef CR_INITPRESSURE

      cr_pressure = CR_INITPRESSURE * SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      SphP[i].Entropy *= (1 - CR_INITPRESSURE);
      q_phys = 1.685;

      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] =
	    cr_pressure / (SphP[i].d.Density / a3 * CR_Tab_Beta(q_phys, CRpop) *
			   (C / All.UnitVelocity_in_cm_per_s) * (C / All.UnitVelocity_in_cm_per_s) / 6.0);

	  SphP[i].CR_C0[CRpop] = C_phys[CRpop] * pow(SphP[i].d.Density, (1.0 - All.CR_Alpha[CRpop]) / 3.0);
	  SphP[i].CR_q0[CRpop] = q_phys * pow(SphP[i].d.Density, -1.0 / 3.0);
	}
#endif

    }

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)

  if(ThisTask == 0)
    {
      printf("Initial abundances (for P[1].ID=%d):\n",P[1].ID);

      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);
      printf("HM=%g, H2I=%g, H2II=%g, elec=%g\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec);

#if defined (UM_CHEMISTRY) && defined(UM_HD_COOLING)
      printf("HD=%g,  DI=%g, DII=%g\n",SphP[1].HD,SphP[1].DI,SphP[1].DII);
      printf("HeHII=%g\n",SphP[1].HeHII);
#endif

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g,\ndensity=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].d.Density, SphP[N_gas - 1].Entropy);



#if defined(UM_METAL_COOLING) && defined (LT_METAL_COOLING)
      printf("::  mass of C=%g, O=%g, Si=%g, Fe=%g, N=%g\n",
	     SphP[N_gas - 1].Metals[Carbon], SphP[N_gas - 1].Metals[Oxygen], SphP[N_gas - 1].Metals[Silicon],
	     SphP[N_gas - 1].Metals[Iron], SphP[N_gas - 1].Metals[Nitrogen]);

#endif
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;


  for(i = 0; i < N_gas; i++)	/* Init Chemistry: */
    {
#ifdef CHEMISTRY
      a_start = All.Time;
      a_end = All.Time + 0.001;	/*  0.001 as an arbitrary value */
      ifunc = compute_abundances(0, i, a_start, a_end);
#endif


#ifdef UM_CHEMISTRY

#if defined (UM_METAL_COOLING) && defined (LT_METAL_COOLING)
      um_ZsPoint = SphP[i].Metals;
      um_mass = P[i].Mass;
#endif
      double um_u;

      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */
      ifunc = compute_abundances(0, i, a_start, a_end, &um_u);

#endif

      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);
    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif

  
  if(RestartFlag == 3)
    {
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
      if(ThisTask == 0)
	printf("SUBFIND_RESHUFFLE_AND_POTENTIAL: Calculating potential energy before reshuffling...\n");
#ifdef PMGRID
      long_range_init_regionsize();
#endif
      compute_potential();
      if(ThisTask == 0)
	printf("potential energy done.\n");

#endif

#ifdef FOF
      fof_fof(RestartSnapNum);
#endif
      endrun(0);
    }

  if(RestartFlag == 5)
    {
      /* calculating powerspec and twopoint function */
#ifndef PM_FROM_L3
#ifdef PMGRID
      long_range_init_regionsize();
#ifdef PERIODIC
      int n, n_type[6];
      long long ntot_type_all[6];
      /* determine global and local particle numbers */
      for(n = 0; n < 6; n++)
	n_type[n] = 0;
      for(n = 0; n < NumPart; n++)
	n_type[P[n].Type]++;
      sumup_large_ints(6, n_type, ntot_type_all);

      calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
      force_treebuild(NumPart, NULL);
      twopoint();
#endif
      endrun(0);
    }

#ifdef CHEMISTRY
  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].d.Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif

  if(RestartFlag == 4)
    {
      All.Time = All.TimeBegin = header.time;
      sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
      if(ThisTask == 0)
        printf("Start writing file %s\n", All.SnapshotFileBase);
      printf("RestartSnapNum %d\n", RestartSnapNum);

      All.TopNodeAllocFactor = 0.008;

      savepositions(RestartSnapNum);
      endrun(0);
    }
}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("masstot = %g\n", masstot);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
#ifdef TIMEDEPGRAV
  omega *= All.Gini / All.G;
#endif

#ifdef KSPACE_NEUTRINOS
  omega += All.OmegaNu;
#endif

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
#ifndef ONEDIM
	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  PPP[i].Hsml = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#else
	  PPP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	  if(All.SofteningTable[0] != 0 && PPP[i].Hsml > 200.0 * All.SofteningTable[0])
	    PPP[i].Hsml = All.SofteningTable[0];
#endif
	}
    }

#ifdef MODGRAV
  if(RestartFlag == 0)
    {
      for(i = 0; i < NumPart; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

	    P[i].Hsml =
	      pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
	  
	    if(P[i].Type == 0)
	      if(P[i].Hsml > 200.0 * All.SofteningTable[0])
		P[i].Hsml = All.SofteningTable[0];
	      
	    if(P[i].Type == 1)
	      if(P[i].Hsml > 200.0 * All.SofteningTable[1])
		P[i].Hsml = All.SofteningTable[1];
	}
    }
  
#endif /* MODGRAV */



#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif

#if defined(RADTRANSFER) || defined(LT_STELLAREVOLUTION)
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 4)
	  PPP[i].Hsml = All.SofteningTable[4];
    }
#endif

#ifdef MODGRAV
  /* compute Rho_0, Rho_1 with the function modgrav_density() in ./modgrav/modgrav_density.c
  
  modgrav_density();

  compute Phi, M_Phi with the function modgrav_fieldsolve in ./modgrav/modgrav_fieldsolve.c

  modgrav_fieldsolve();
  */
#endif

  density();
#ifdef VECT_POTENTIAL
  smoothed_values();
  density();
#endif

#if defined(MAGNETIC) && defined(BFROMROTA)
  if(RestartFlag == 0)
    {
      if(ThisTask == 0)
	printf("Converting: Vector Potential -> Bfield\n");
      rot_a();
    }
#endif

}


void assign_unique_ids(void)
{
  int i, *numpartlist;
  MyIDType idfirst;

  numpartlist = (int *)mymalloc("numpartlist", NTask * sizeof(int));

  MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);

  idfirst = 1;

  for(i = 0; i < ThisTask; i++)
    idfirst += numpartlist[i];

  for(i = 0; i < NumPart; i++)
    {
      P[i].ID = idfirst;
      idfirst++;
    }

  myfree(numpartlist);
}



void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  if(ThisTask == 0)
    {
      printf("Testing ID uniqueness...\n");
      fflush(stdout);
    }

  if(NumPart == 0)
    {
      printf("need at least one particle per cpu\n");
      endrun(8);
    }

  t0 = second();

#ifndef SPH_BND_PARTICLES
  ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < NumPart; i++)
    if(ids[i] == ids[i - 1])
      {
#ifdef LONGIDS
	printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
	       (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);

#else
	printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
	endrun(12);
      }

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    if(ids[NumPart - 1] == ids_first[ThisTask + 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
	endrun(13);
      }

  myfree(ids_first);
  myfree(ids);
#endif

  t1 = second();

  if(ThisTask == 0)
    {
      printf("success.  took=%g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}
