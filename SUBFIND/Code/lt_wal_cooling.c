#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

#ifdef LT_METAL_COOLING_WAL
#include <hdf5.h>

#include "lt_wal_cooling.h"
#include "lt_error_codes.h"


/* :: -------------------------------------------- ::  */
/*    SEGMENT .CODE                                    */



/* :: ---------------------------- PUBLIC ROUTINES ::  */
/*                                                     */
/* +  WalCoolInitialize(void)                          */
/*        initializes tables and memory                */
/*                                                     */
/* +  WalCool_tables_load()                            */
/*        checks whether or not a new table(s) must    */
/*        be loaded and loads them.                    */
/*                                                     */

int DBG=0;

void read_cooling_tables_dummy()
{
  int i;

  ZBins = 1;
  TBins = 1;
  ZMin = ZMax = -4.0;
  TMin = 4.0;
  TMax = 8.0;

  if(ThisTask == 0)
    printf("Setting dummy cooling tables to zero for WAL cooling ...\n");

  CoolTvalue = (double *) mymalloc("CoolTvalue", (TBins + ZBins + TBins * ZBins) * sizeof(double));
  memset(CoolTvalue, 0, (TBins + ZBins + TBins * ZBins) * sizeof(double));

  CoolingTables = (double **) mymalloc("CoolingTables", ZBins * sizeof(double *));
  memset(CoolingTables, 0, ZBins * sizeof(double *));

  CoolZvalue = CoolTvalue + TBins;
  CoolingTables[0] = CoolZvalue + ZBins;
  for(i = 1; i < ZBins; i++)
    CoolingTables[i] = CoolingTables[i - 1] + TBins;

}

void InitCool(void)
{
  read_cooling_tables_dummy();
  if(ThisTask ==0)
    printf("Initializing cooling using WAL tables ...\n");
  WalCool_Initialize();
  return;
}

void WalCool_Initialize()
/* initialize the cooling: loads the common tables and allocate memory to store 2 cooling tables */
{

  if(ThisTask == 0)
    {
      WalCool_get_redshift_table();
      printf("\n[Wiersma et al. Cooling] %d tables found, spanning redshift from %g downto %g %g\n\n", WalCool_CoolTables_num,
             WalCool_CoolTables_redshifts[WalCool_CoolTables_num - 1].redshift, WalCool_CoolTables_redshifts[0].redshift, max_cool_redshift);
    }

  MPI_Bcast(&WalCool_CoolTables_num, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_cool_redshift, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    WalCool_CoolTables_redshifts = mymalloc("WalCool_CoolTables_redshifts", sizeof(struct my_direntry) * WalCool_CoolTables_num);

  MPI_Bcast(WalCool_CoolTables_redshifts, sizeof(struct my_direntry) * WalCool_CoolTables_num, MPI_BYTE, 0, MPI_COMM_WORLD);

                                                                         /* get general informations from headers */
  WalCool_Initialize_get_header();

                                                                         /* allocate space for tables */
  WalCoolMfreeSize = WalCool_n_Hef * WalCool_n_Rho * WalCool_n_T;
  WalCoolSize      = WalCool_n_El * WalCool_n_Rho * WalCool_n_T;
  
  WALCOOLTABLES = (float*)mymalloc("WALCOOLTABLES", sizeof(float) * 2 * (WalCoolMfreeSize + WalCoolSize));

                                                                         /* set pointers to H+He cooling */
  WalMfreeCoolTables = WALCOOLTABLES;

                                                                         /* set pointer to metal cooling tables */
  WalCoolTables = WALCOOLTABLES + 2 * WalCoolMfreeSize;
  

  IDX = (int*) mymalloc("IDX", sizeof(int) * 4);
  dX  = (double*) mymalloc("dX", sizeof(double) * 4);

  WalCool_redshift_index = -1;
  HighZTable = 0;
  LowZTable = 1;
  
  return;
}


void WalCool_tables_load(double redshift)
{
  int change = 0;

  if(redshift > max_cool_redshift)
    {
      if(collisional_table_loaded)
	return;

      WalCool_redshift_index = -1;
      WalCool_get_collis_table();
      collisional_table_loaded = 1;
      if(ThisTask == 0)
	printf("[WalCool] getting collisional table\n"); fflush(stdout);

    }
  else
    {
      if(WalCool_redshift_index == -1)
        {
          WalCool_redshift_index = WalCool_CoolTables_num - 1;
          change = 2;
        }
      else if(WalCool_redshift_index == 0)
        return;
  
      for(; (redshift < WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift) && (WalCool_redshift_index >=0); WalCool_redshift_index--)
	change++;
      
      if(change == 1)
	{

/* 	  memcpy(&WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(HighZTable)], &WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(LowZTable)], sizeof(float) * WalCoolMfreeSize); */
/* 	  memcpy(&WalCoolTables[COOLRATE_z_IDX(HighZTable)], &WalCoolTables[COOLRATE_z_IDX(LowZTable)], sizeof(float) * WalCoolSize); */

/* 	  memcpy(&WalCool_enHS[ENHS_z_IDX(HighZTable)], &WalCool_enHS[ENHS_z_IDX(LowZTable)], sizeof(float) * WalCool_arraysS_size); */
/* 	  memcpy(&WalCool_enH[ENH_z_IDX(HighZTable)], &WalCool_enH[ENH_z_IDX(LowZTable)], sizeof(float) * WalCool_arrays_size); */
/* 	  memcpy(&WalCool_mu[MU_z_IDX(HighZTable)], &WalCool_mu[MU_z_IDX(LowZTable)], sizeof(float) * WalCool_arrays_size); */
/* 	  memcpy(&WalCool_UtoT[UtoT_z_IDX(HighZTable)], &WalCool_UtoT[UtoT_z_IDX(LowZTable)], sizeof(float) * WalCool_arrays_size); */


	  HighZTable ^= 1;
	  LowZTable ^= 1;

	  if(ThisTask == 0)
	    printf("[WalCool] sup = %d, inf = %d, getting %d in inf, z=%f\n", HighZTable, LowZTable, WalCool_redshift_index, WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift);
	  WalCool_get_table(WalCool_redshift_index, LowZTable);

	  MPI_Bcast(&WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(LowZTable)], WalCoolMfreeSize * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&WalCoolTables[COOLRATE_z_IDX(LowZTable)], WalCoolSize * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enHS[ENHS_z_IDX(LowZTable)], WalCool_arraysS_size * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enH[ENH_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_mu[MU_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_UtoT[UtoT_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);

	  return;
	}
      else if(change > 1)
	{
	  if(ThisTask == 0)
	    printf("[WalCool] sup = %d, inf = %d, getting %d(z=%f) in inf, %d(z=%f) in sup\n", HighZTable, LowZTable, WalCool_redshift_index, WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift, WalCool_redshift_index+1, WalCool_CoolTables_redshifts[WalCool_redshift_index+1].redshift);

	  WalCool_get_table(WalCool_redshift_index, LowZTable);
	  WalCool_get_table(WalCool_redshift_index+1, HighZTable);
	}
    }

  MPI_Bcast(WALCOOLTABLES, (2 * WalCoolMfreeSize + WalCoolSize) * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(WalCool_enHS, (3 * 2 * WalCool_arrays_size + 2 * WalCool_arraysS_size) * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
            
  return;
}

int get_cool_redshift(double Redshift, double *DZ)
{
  /* find index for redshift */
  if(WalCool_redshift_index >= 0)
    *DZ = (Redshift - WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift) / (WalCool_CoolTables_redshifts[WalCool_redshift_index + 1].redshift - WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift);
  else
    *DZ = 0;
  return 0;
}

int get_cool_n_el()
{
  return WalCool_n_Ab;
}

int Is_a_Coolant(int i)
{
  return Specie_is_coolant[i];
}


void set_cooltable_index(int idx)
{
  WalCool_redshift_index = idx;
  return;
}


double get_max_cool_redshift()
{
  return max_cool_redshift;
}


double *set_metallicities(int i, double *Metallicities, double dens_to_phys_fact)
/*
 * dens_to_phys_fact must contain all factors needed to convert densities
 * from code units to physical code units (also cosmological factors must be
 * included)
 */
{
  double   nonHydmass;
  MyDouble mass;
  int pos, j, k;

  char buffer[1000];

  dens_to_phys_fact *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  
  mass = P[i].Mass;
  if(P[i].Type != 0)
    {
      printf("[%d] : error : request to calculate metallicities for particle"
             " %d that is not a gas particle\n", ThisTask, i);
      endrun(LT_ERR_WALCOOL_Z_FOR_NOTGAS);
    }
  
  for(pos = 0, j = 0; j < WalCool_n_Ab; j++, pos++)
    if(SpeciesIdx[j] >= 0)
      {
        if(j == myHyd)
          {
            for(nonHydmass = 0, k = 0; k < LT_NMetP; k++)
	      nonHydmass += SphP[i].Metals[k];
            Metallicities[SpeciesPos[j]] = (double)(mass - nonHydmass) / mass * SphP[i].d.Density * dens_to_phys_fact; 
          }
        else
#if defined(LT_METAL_COOLING_on_SMOOTH_Z)	  
	    Metallicities[SpeciesPos[j]] = (double)SphP[i].Zsmooth[SpeciesIdx[j]];
#else
	    Metallicities[SpeciesPos[j]] = (double)SphP[i].Metals[SpeciesIdx[j]] / mass;
#endif
        
	if( j == myHel && UseHeNumberRatio)
	  Metallicities[SpeciesPos[j]] *= SphP[i].d.Density * dens_to_phys_fact;
      }
    else
      printf("%d element %d has no idx\n", ThisTask, j);
  
  Metallicities[HydPos] /= PROTONMASS;
  if(UseHeNumberRatio)
    Metallicities[HelPos] /= (4 * PROTONMASS * Metallicities[HydPos]); /* it would be better to use number ratio, however for    
                                                                        * some unknown reason, the hdf5 lib does not read in the 
                                                                        * Helium_number_ratio_bins data set                      
                                                                        */
  return Metallicities;
}

#ifdef SUBFIND

double *set_metallicities_subfind(int i, double *Metallicities, double dens_to_phys_fact)
/*
 * dens_to_phys_fact must contain all factors needed to convert densities
 * from code units to physical code units (also cosmological factors must be
 * included)
 */
{
  double   nonHydmass;
  MyDouble mass;
  int pos, j, k;

  char buffer[1000];

  dens_to_phys_fact *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  
  mass = P[i].Mass;
  if(P[i].Type != 0)
    {
      printf("[%d] : error : request to calculate metallicities for particle"
             " %d that is not a gas particle\n", ThisTask, i);
      endrun(LT_ERR_WALCOOL_Z_FOR_NOTGAS);
    }
  
  for(pos = 0, j = 0; j < WalCool_n_Ab; j++, pos++)
    if(SpeciesIdx[j] >= 0)
      {
        if(j == myHyd)
          {
            for(nonHydmass = 0, k = 0; k < LT_NMetP; k++)
	      nonHydmass += SphP[P[i].origindex].Metals[k];
            Metallicities[SpeciesPos[j]] = (double)(mass - nonHydmass) / mass * SphP[P[i].origindex].d.Density * dens_to_phys_fact; 
          }
        else
#if defined(LT_METAL_COOLING_on_SMOOTH_Z)	  
	    Metallicities[SpeciesPos[j]] = (double)SphP[P[i].origindex].Zsmooth[SpeciesIdx[j]];
#else
	    Metallicities[SpeciesPos[j]] = (double)SphP[P[i].origindex].Metals[SpeciesIdx[j]] / mass;
#endif
        
	if( j == myHel && UseHeNumberRatio)
	  Metallicities[SpeciesPos[j]] *= SphP[P[i].origindex].d.Density * dens_to_phys_fact;
      }
    else
      printf("%d element %d has no idx\n", ThisTask, j);
  
  Metallicities[HydPos] /= PROTONMASS;
  if(UseHeNumberRatio)
    Metallicities[HelPos] /= (4 * PROTONMASS * Metallicities[HydPos]); /* it would be better to use number ratio, however for    
                                                                        * some unknown reason, the hdf5 lib does not read in the 
                                                                        * Helium_number_ratio_bins data set                      
                                                                        */
  return Metallicities;
}
  
#endif

int set_indexes(double U, double DZ, double *Metallicities)
/*
 * quantities must be given in physical units (Metallicities, for instance,
 * as returned by set_metallicities() )
 */
{
  /* find index for redshift */
  dz = DZ;
  IDX[z_IDX] = WalCool_redshift_index;

  /* if(((dz = DZ)< 0) && (IDX[z_IDX] = WalCool_redshift_index) >= 0) */
  /*   dX[z_IDX] = (Redshift - All.WalCool_redshifts[WalCool_redshift_index]) / (All.WalCool_redshifts[WalCool_redshift_index + 1] - All.WalCool_redshifts[WalCool_redshift_index]);*/
  
  /* find index for hydrogen number density */
  IDX[H_IDX] = find_index(LOG, WalCool_Rho, Rmin, Rrange, WalCool_n_Rho, Metallicities[HydPos], &dX[H_IDX]);
  
  /* find index for the ration between he number density and H number density */
  IDX[He_IDX] = find_index(LIN, WalCool_Hef, WalCool_Hef[0], WalCool_Hef[WalCool_n_Hef-1] - WalCool_Hef[0], WalCool_n_Hef, Metallicities[HelPos], &dX[He_IDX]);
  
  /* find index for internal energy */ 
  IDX[U_IDX] = find_index(LOG, WalCool_U, Umin, Urange, WalCool_n_T, U, &dX[U_IDX]);
  
  return 0;      
}




double GetCoolingTime(double U, double Rho, double Redshift, double DZ, double *Metallicities)
/* quantities must be in physical code units
 */
{
  double Lambda, CTime, ratefact;
                                                                         /* convert in physical cgs units */
  Rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  U   *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  ratefact = Metallicities[HydPos] * Metallicities[HydPos] / Rho;
  
  if((Lambda = CoolingRateFromU(U, Redshift, DZ, Metallicities)) > 0)
    CTime = 0;
  else
    {
      CTime = U / (-ratefact * Lambda);
      CTime *= All.HubbleParam / All.UnitTime_in_s;
    }
    
  return CTime;
}


double convert_u_to_temp(double U, double Redshift, double DZ, double *Metallicities)
/*
 * quantities must be given in physical units (Metallicities, for instance,
 * as returned by set_metallicities)
 */  
{
  double T;

  set_indexes(U, DZ, Metallicities);

  if(WalCool_redshift_index >= 0)
    {      
/*       for(i = 0; i < (1 << 4); i++) */
/*         { */
/*           for(i = 0; i < 4; i++) */
/*             S[i] = j & (1 << i); */
          
/*           V[i] = WalCool_UtoT[UtoT_IDX(WalCool_redshift_index + S[0], */
/*                                        He_idx + S[1], */
/*                                        T_idx + S[2], */
/*                                        R_idx + S[3])]; */
/*         } */

/*       T = get_Ndim_interp(4, V, dC) */


      T = (1 - dz)*(1 - dHe)*(1 - dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU, inH)] +
          (dz)*(1 - dHe)*(1 - dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU, inH)] +
          (1 - dz)*(dHe)*(1 - dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU, inH)] +
          (dz)*(dHe)*(1 - dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU, inH)] +
          (1 - dz)*(1 - dHe)*(dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU + 1, inH)] +
          (dz)*(1 - dHe)*(dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU + 1, inH)] +
          (1 - dz)*(dHe)*(dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU + 1, inH)] +
          (dz)*(dHe)*(dU)*(1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU + 1, inH)] +
          (1 - dz)*(1 - dHe)*(1 - dU)*(dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU, inH + 1)] +
          (dz)*(1 - dHe)*(1 - dU)*(dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU, inH + 1)] +
          (1 - dz)*(dHe)*(1 - dU)*(dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU, inH + 1)] +
          (dz)*(dHe)*(1 - dU)*(dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU, inH + 1)] +
          (1 - dz)*(1 - dHe)*(dU)*(dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU + 1, inH + 1)] +
          (dz)*(1 - dHe)*(dU)*(dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU + 1, inH + 1)] +
          (1 - dz)*(dHe)*(dU)*(dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU + 1, inH + 1)] +
          (dz)*(dHe)*(dU)*(dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU + 1, inH + 1)];

    }
  else
    {
      //iz  = 0;
      //inH = 0;

      T = (1 - dHe)*(1 - dU) * WalCool_UtoT[UtoT_collis_IDX(iHe, iU)] +
          dHe * (1 - dU) * WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU)] +
          (1 - dHe) * dU * WalCool_UtoT[UtoT_collis_IDX(iHe, iU + 1)] +
          dHe * dU * WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU + 1)];

      if(DBG)
	printf("\t\t[cut][%d][%d][%u] %d %d %g %g :: %g %g %g %g :: %g\n",
	       ThisTask, All.NumCurrentTiStep, PID,
	       iHe, iU, dHe, dU,
	       WalCool_UtoT[UtoT_collis_IDX(iHe, iU)], 
	       WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU)],
	       WalCool_UtoT[UtoT_collis_IDX(iHe, iU + 1)],
	       WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU + 1)],
	       T);
	       
    }

  return T;
}



double CoolingRateFromU(double U, double Redshift, double DZ, double *Metallicities)
{
  float Temperature;
  Temperature = convert_u_to_temp(U, Redshift, DZ, Metallicities);

  IDX[U_IDX] = find_index(LOG, WalCool_T, Tmin, Trange, WalCool_n_T, (float)Temperature, &dX[U_IDX]);

  return CoolingRate(Temperature, Redshift, Metallicities);
}


double InnerInterpolation(float *zhigh, float *zlow)
{
  double Lambda;
  double LambdaH, LambdaL;

  LambdaH = (1 - dU)*(1 - dnH) * zhigh[COOLRATE_IIDX(iU, inH)] + 
    (dU)*(1 - dnH) * zhigh[COOLRATE_IIDX(iU + 1, inH)] + 
    (1 - dU)*(dnH) * zhigh[COOLRATE_IIDX(iU, inH + 1)] + 
    (dU)*(dnH) * zhigh[COOLRATE_IIDX(iU + 1, inH + 1)];
    
  LambdaL = (1 - dU)*(1 - dnH) * zlow[COOLRATE_IIDX(iU, inH)] + 
    (dU)*(1 - dnH) * zlow[COOLRATE_IIDX(iU + 1, inH)] + 
    (1 - dU)*(dnH) * zlow[COOLRATE_IIDX(iU, inH + 1)] + 
    (dU)*(dnH) * zlow[COOLRATE_IIDX(iU + 1, inH + 1)];

  Lambda = (1 - dz) * LambdaH + dz * LambdaL;

  return Lambda;
}


double CoolingRate(double Temp, double Redshift, double *Metallicities)
{
  double Lambda1, Lambda2, Lambda, LambdaCmpt, T_cmpt;
  double ne_He1, ne_He2, ne, neS, fact_Z, fact_ne;
  int    el;

  /* H + He cooling first */

  LambdaCmpt = 0;
  if(iz >= 0)
                                                                         /* there is a photo-ionizing background */
    {
      Lambda1 = InnerInterpolation(&WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(HighZTable, iHe)],
                                   &WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(LowZTable, iHe)]);
      Lambda2 = InnerInterpolation(&WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(HighZTable, iHe+1)],
                                   &WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(LowZTable, iHe+1)]);

      ne_He1 = InnerInterpolation(&WalCool_enH[ENH_IDX(HighZTable, iHe, 0, 0)],
                                  &WalCool_enH[ENH_IDX(LowZTable, iHe, 0, 0)]);
      ne_He2 = InnerInterpolation(&WalCool_enH[ENH_IDX(HighZTable, iHe+1, 0, 0)],
                                  &WalCool_enH[ENH_IDX(LowZTable, iHe+1, 0, 0)]);

      ne     = dHe * ne_He1 + (1 - dHe) * ne_He2;

      neS = InnerInterpolation(&WalCool_enHS[ENHS_z_IDX(HighZTable)],
                               &WalCool_enHS[ENHS_z_IDX(LowZTable)]);
    }
  else
                                                                         /* no UV background */
    {
                                                                         // interpolates ne in T separately for each He
      ne_He1 = dU * WalCool_enH[ENH_collis_IDX(iHe, iU)] + (1 - dU) * WalCool_enH[ENH_collis_IDX(iHe, iU+1)];
      ne_He2 = dU * WalCool_enH[ENH_collis_IDX(iHe+1, iU)] + (1 - dU) * WalCool_enH[ENH_collis_IDX(iHe+1, iU+1)];
                                                                         // interpolates ne between the two He
      ne     = dHe * ne_He1 + (1 - dHe) * ne_He2;

      neS    = dU * WalCool_enHS[ENHS_collis_IDX(iU)] + (1 - dU) * WalCool_enHS[ENHS_collis_IDX(iU+1)];

                                                                         // adding inverse Comp. Scatt. off the CMB
      if(All.ComovingIntegrationOn)
        {
	  double redshift_fact;
          
          redshift_fact = (1 + Redshift) * (1 + Redshift);
          redshift_fact *= redshift_fact;

	  LambdaCmpt = 5.65e-36 * ne * (Temp - TCMB(Redshift)) * redshift_fact / Metallicities[HydPos]; // / Metallicities[HydPos];
        }

                                                                         // interpolates lambda separately for each He
      Lambda1 = dU * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU)] +
        (1 - dU) * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU+1)];
      Lambda2 = dU * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe+1, iU)] +
        (1 - dU) * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe+1, iU+1)];

      if(DBG)
	printf("\t\t[CR][%d][%d][%u] %d %d %g %g :: %g %g %g %g :: %g %g %g\n",
	       ThisTask, All.NumCurrentTiStep, PID,
	       iHe, iU, dHe, dU,
	       WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU)],
	       WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU+1)],
	       WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe+1, iU)],
	       WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe+1, iU+1)],
	       LambdaCmpt,
	       Lambda1, Lambda2);
      
    }

                                                                         // interpolates lambda between the two He and adds Compton
  Lambda = LambdaCmpt +
    dHe * Lambda1 + (1 - dHe) * Lambda2;

  fact_ne = ne / neS;
  
  /* add all the elements */
  
  for(el = 0; el < WalCool_n_El; el++)
    {
      fact_Z = fact_ne * Metallicities[el] / WalCool_Ab[el];
      
      if(iz >= 0)
        Lambda += fact_Z * InnerInterpolation(&WalCoolTables[COOLRATE_el_IDX(HighZTable, el)],
					      &WalCoolTables[COOLRATE_el_IDX(LowZTable, el)]) ;
      else
	{
	  Lambda += fact_Z * (dU * WalCoolTables[COOLRATE_collis_IDX(iU, el)] +
			      (1 - dU) * WalCoolTables[COOLRATE_collis_IDX(iU+1, el)]);
	  if(DBG)
	    printf("\t\t[CR][%d][%d][%u] %d %d %s %g %g :: %g %g %g %g %g\n",
		   ThisTask, All.NumCurrentTiStep, PID,\
		   iU, el, WalCool_El_names[el], fact_Z, fact_ne, dU,
		   WalCoolTables[COOLRATE_collis_IDX(iU, el)],
		   WalCoolTables[COOLRATE_collis_IDX(iU+1, el)],
		   fact_Z * (dU * WalCoolTables[COOLRATE_collis_IDX(iU, el)] +
			   (1 - dU) * WalCoolTables[COOLRATE_collis_IDX(iU+1, el)]),
		   Lambda);
	}
    }

  if(DBG)
    DBG = 0;
  return -Lambda;
}


double DoCooling(double U_in, double Rho, double *Metallicities, double Redshift, double DZ, double dt)
/*
 * note: all quantitites must be given in physical code units !
 * note: assumes that the Metallicities array has been given by set_metallicities()
 */
{

  double u, u_bottom, u_top, du;
  double ratefact, Lambda;
  int    iter;
                                                                         /* convert in physical cgs units */
  Rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  U_in*= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt  *= All.UnitTime_in_s / All.HubbleParam;

                                                                         /* in Metallicities H abundance must be in number */
                                                                         /* abundance and He's one in number ratio over H  */
  u = U_in;
  u_bottom = u_top = u;

  ratefact = Metallicities[HydPos] * Metallicities[HydPos] / Rho;

  Lambda = CoolingRateFromU(u, Redshift, DZ, Metallicities);

                                                                         /* cooling/heating are too weak to calculate */
  du = ratefact * Lambda * dt;
  if(fabs(du / u) < 1.0e-6)
    {
      u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;/* to internal units */
      return u;
    }

                                                                         /* bracketing  */
  iter = 0;
  if(u - U_in - ratefact * Lambda * dt < 0)	// heating 
  {
    u_top *= sqrt(1.15);
    u_bottom /= sqrt(1.15);
    while(u_top - U_in - ratefact * CoolingRateFromU(u_top, Redshift, DZ, Metallicities) * dt < 0)
    {
      iter++;
      u_top *= 1.15;
      u_bottom *= 1.15;
      if(iter > 180)
	printf("\t\t[WalCool][Hb][%d][%u][%d]  %g %g %g %g %g %g %g\n", ThisTask, PID, iter, u_bottom, u_top, U_in, Rho, dt, CoolingRateFromU(u_top, Redshift, DZ, Metallicities), ratefact * CoolingRateFromU(u_top, Redshift, DZ, Metallicities) * dt); fflush(stdout);
      if(iter > 200)
	endrun(110011001);
    }    
  }

  iter = 0;
  if(u - U_in - ratefact * Lambda * dt > 0) // cooling
  {
    u_top *= sqrt(1.15);
    u_bottom /= sqrt(1.15);
    while(u_bottom - U_in - ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities) * dt > 0)
    {
      iter++;
      u_top /= 1.15;
      u_bottom /= 1.15;
      if((iter > 180)) // || (PID==80568 && All.NumCurrentTiStep==118))
	{
	  printf("\t\t[WalCool][Cb][%d][%u][%d] %g %g %g %g %g %g %g :: %d %d %d\n",
		 ThisTask, PID, iter, 
		 u_bottom, u_top, U_in, Rho, dt, 
		 CoolingRateFromU(u_top, Redshift, DZ, Metallicities), 
		 ratefact * CoolingRateFromU(u_top, Redshift, DZ, Metallicities) * dt,
		 inH, iHe, iU); fflush(stdout);
	  DBG = 1;
	}
      if(iter > 200)
	endrun(110011002);
    }
  }

  iter = 0;

  do // iterate to convergence
  {
    u = 0.5 * (u_bottom + u_top);
    
    Lambda = CoolingRateFromU(u, Redshift, DZ, Metallicities);
    
    if(u - U_in - ratefact * Lambda * dt > 0)
      {
	u_top = u;
      }
    else
      {
	u_bottom = u;
      }
    
    du = u_top - u_bottom;
    
    iter++;
    
    if(iter >= (MAXITER - 10))
      printf("u= %g\n", u);
  }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);
  
  if(iter >= MAXITER)
    printf("failed to converge in DoCooling()\n");

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;/* to internal units */
	
  return u;
  
}









/* :: -------------------------- INTERNAL ROUTINES ::  */


int WalCool_get_redshift_table(void)
/* this routine gets from the tables' directory what are the redshifts at
   which the tables have been defined */
{

  int n;
  DIR    *WalCool_CoolTables_dir;
  struct dirent *dentry;


                                                                         /* as first, task 0 gets the content of the directory */
                                                                         /* we DO NOT use scandir because it uses the system's */
                                                                         /* malloc routine and we want to stick onto mymalloc  */
                                                                         /* (isnt'it ?)                                        */

  if( (WalCool_CoolTables_dir = opendir(All.WalCool_CoolTables_path)) == NULL)
    {
      printf("The directory %s that should contain the cooling tables from Wiersma et al. does not exist\n", All.WalCool_CoolTables_path);
      exit(9191000);
    }
  
  n = 0;                                                                 /* this first cycle is needed to determine how many   */
                                                                         /* tables are in the directory                        */
  while(dentry = readdir(WalCool_CoolTables_dir))
    if( is_it_a_tablefile(dentry->d_name, -1) )
      n++;

  WalCool_CoolTables_num = n;
  
  closedir(WalCool_CoolTables_dir);
                                                                         /* re-open the directory */
  WalCool_CoolTables_dir = opendir(All.WalCool_CoolTables_path);

  WalCool_CoolTables_redshifts = mymalloc("WalCool_CoolTables_redshifts", sizeof(struct my_direntry) * n); /* allocate memory to store the entries  */

  n = 0;
  while(dentry = readdir(WalCool_CoolTables_dir))
    if( is_it_a_tablefile(dentry->d_name, n) )
      n++;
    
  closedir(WalCool_CoolTables_dir);

  qsort(WalCool_CoolTables_redshifts, WalCool_CoolTables_num, sizeof(struct my_direntry), compare_dir_redshifts);

  max_cool_redshift = WalCool_CoolTables_redshifts[WalCool_CoolTables_num - 1].redshift;
  return WalCool_CoolTables_num;
}



int is_it_a_tablefile(char *name, int n)
/* a file is a table file whether its name is made by the following pattern */
/*      z_[:digit:].[:digit:].hdf5                                          */
/* where [:digit:] means un unknown number od digits                        */
{
  int  pos, dotpos, ndots;
  
  pos = strlen(name);

  if(pos < 10)                                                           /* check for the minimum expected lenght */
    return 0;

  if(strcmp(&name[pos - 5], ".hdf5") != 0x0)                             /* check whether the suffix is ".hdf5" */
    return 0;

  if(strncmp(name, "z_", 2) != 0x0)                                      /* check whether the prefix is "z_" */
    return 0;

  dotpos = pos-5;                                                        /* points to the '.' */

  for(ndots = 0, pos = 2; pos < dotpos; pos++)
    {      
      if(!isdigit(name[pos]))
        {
          if(name[pos] == '.')
            {
              if(ndots == 1)
                return 0;
              else
                ndots = 1;
            }
          else
            return 0;
        }
    }

  if(n >= 0)
    {
      strncpy(WalCool_CoolTables_redshifts[n].redshift_str, &name[2], dotpos-2);
      WalCool_CoolTables_redshifts[n].redshift_str[dotpos-2]='\0';
      WalCool_CoolTables_redshifts[n].redshift = atof(WalCool_CoolTables_redshifts[n].redshift_str);
    }
  
  return 1;
}


void WalCool_get_table(int redshift_index, int tab_index)
/* load the cooling table corresponding to a given entry in the redshift table;
   the table is loaded in the tab_index position */
{
  char   filename[300], dname[300];
  hid_t  fid, dset;
  herr_t hd_report;

  int    element;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_%s.hdf5", All.WalCool_CoolTables_path, WalCool_CoolTables_redshifts[redshift_index].redshift_str);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      
                                                                         /* get cooling from metal-free: H + He */
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, "/Metal_free/Net_Cooling");
#else
      dset      = H5Dopen(fid, "/Metal_free/Net_Cooling", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);
      
      
      
                                                                         /* get cooling for metals */
      for(element = 0; element < WalCool_n_El; element++)
        {
          sprintf(dname, "/%s/Net_Cooling", WalCool_El_names[element]);
#ifdef OLD_HDF5
          dset      = H5Dopen(fid, dname);
#else
          dset      = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
          hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCoolTables[COOLRATE_el_IDX(tab_index, element)]);
          hd_report = H5Dclose(dset);
        }
      
                                                                         /* get ne / nH table */
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dset      = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_enH[ENH_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);


                                                                         /* get ne / nH solar table */
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, "/Solar/Electron_density_over_n_h");
#else
      dset      = H5Dopen(fid, "/Solar/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_enHS[ENHS_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);

                                                                         /* get mu table */
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dset      = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_mu[MU_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);

      
                                                                         /* get U to T table */      
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, "/Metal_free/Temperature/Temperature");
#else
      dset      = H5Dopen(fid, "/Metal_free/Temperature/Temperature", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_UtoT[UtoT_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);

      hd_report = H5Fclose(fid);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  return;  
}


void WalCool_get_collis_table()
/* loads the purely collisional cooling table, to be used at high redshift */
{
  char   filename[300], dname[300];
  hid_t  fid, dset;
  herr_t hd_report;

  float  *rtable;
  int    element, Tidx;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_collis.hdf5", All.WalCool_CoolTables_path);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

                                                                         /* - elecron density over n_h  NO BCKGRND */
#ifdef OLD_HDF5
      dset        = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dset        = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enH);
      hd_report   = H5Dclose(dset);


                                                                         /* - solar elecron density over n_h  NO BCKGRND */
#ifdef OLD_HDF5
      dset        = H5Dopen(fid, "/Solar/Electron_density_over_n_h");
#else
      dset        = H5Dopen(fid, "/Solar/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enHS);
      hd_report   = H5Dclose(dset);


                                                                         /* - mean particle mass  NO BCKGRND */
#ifdef OLD_HDF5
      dset        = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dset        = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_mu);
      hd_report   = H5Dclose(dset);


                                                                         /* - temperature conversion  NO BCKGRND */
#ifdef OLD_HDF5
      dset        = H5Dopen(fid, "/Metal_free/Temperature/Temperature");
#else
      dset        = H5Dopen(fid, "/Metal_free/Temperature/Temperature", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_UtoT);
      hd_report   = H5Dclose(dset);

      
                                                                         /* get cooling from metal-free: H + He */
      sprintf(dname, "/Metal_free/Net_Cooling");
#ifdef OLD_HDF5
      dset      = H5Dopen(fid, dname);
#else
      dset      = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalMfreeCoolTables);
      hd_report = H5Dclose(dset);

      
      rtable = (float*)mymalloc("temporary_table", sizeof(float) * WalCool_n_T);
                                                                         /* get cooling from metals */
      for(element = 0; element < WalCool_n_El; element++)
        {
          sprintf(dname, "/%s/Net_Cooling", WalCool_El_names[element]);
#ifdef OLD_HDF5
          dset      = H5Dopen(fid, dname);
#else
          dset      = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
          hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rtable);
          hd_report = H5Dclose(dset);

          for(Tidx = 0; Tidx < WalCool_n_T; Tidx++)
            WalCoolTables[COOLRATE_collis_IDX(Tidx, element)] = rtable[Tidx];
        }

      myfree(rtable);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  return;  
  
}

void WalCool_Initialize_get_header(void)
/* loads everything is needed to use the tables: Rho, t, He tabs etc. */
{
#define NAME_LEN 20
  char   filename[300];
  hid_t  fid, dataset;
  herr_t hd_report;

  int    dimensions[6], k = 0, i, pos, size;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_%s.hdf5", All.WalCool_CoolTables_path, WalCool_CoolTables_redshifts[WalCool_CoolTables_num-1].redshift_str);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

                                                                         /* get the arrays dimensions */
      
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_density_bins");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_density_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Rho);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Rho;

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_temperature_bins");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_T);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_T;

      /* NOTE :: /Header/Number_of_helium_fractions SEEMS to BE ZERO!
         We Use an alternative way, however.
       
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_helium_fractions");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Hef);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Hef;
      */

      if(UseHeNumberRatio)
#ifdef OLD_HDF5
        dataset = H5Dopen(fid, "/Metal_free/Helium_number_ratio_bins");
#else
        dataset = H5Dopen(fid, "/Metal_free/Helium_number_ratio_bins", H5P_DEFAULT);
#endif
      else
#ifdef OLD_HDF5
        dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins");
#else
        dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins", H5P_DEFAULT);
#endif
      WalCool_n_Hef = H5Dget_storage_size(dataset) / sizeof(float);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Hef;
      
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_species");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_species", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_El);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_El;

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Abundances/Number_of_abundances");
#else
      dataset = H5Dopen(fid, "/Header/Abundances/Number_of_abundances", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Ab);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Ab;
    }

  MPI_Bcast(dimensions, 6, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(ThisTask != 0)
    {
      k = 0;
      WalCool_n_Rho = dimensions[k++];
      WalCool_n_T   = dimensions[k++];
      WalCool_n_Hef = dimensions[k++];
      WalCool_n_El  = dimensions[k++];
      WalCool_n_Ab  = dimensions[k++];
    }
  
  enH_mu_Size = WalCool_n_Rho * WalCool_n_T * WalCool_n_Hef;
  UtoT_Size   = WalCool_n_Rho * WalCool_n_T * WalCool_n_Hef;
  
                                                                         /* allocate memory for array */

  WalCool_arrays_size  = WalCool_n_Hef * WalCool_n_T * WalCool_n_Rho;
  WalCool_arraysS_size = WalCool_n_T * WalCool_n_Rho;

  size = (WalCool_n_Rho   +                                              /* density bins */
          2 * WalCool_n_T +                                              /* T and U bins */
          WalCool_n_Hef   +                                              /* Helium fractions */
          WalCool_n_Ab    +                                              /* solar abundances */
          3 * 2 * WalCool_arrays_size +                                  /* mean molecular mass, electrons over n_H and temp. conversion @ 2 redshifts*/
	  2 * WalCool_arraysS_size);                                     /* ne over nH for solar composition */

  
  WalCool_indexes_and_arrays = (float *) mymalloc("WalCool_indexe_and_arrays", sizeof(float) * size); 

  WalCool_Rho = WalCool_indexes_and_arrays;
  WalCool_T   = WalCool_Rho + WalCool_n_Rho;
  WalCool_U   = WalCool_T   + WalCool_n_T;
  WalCool_Hef = WalCool_U   + WalCool_n_T;
  WalCool_Ab  = WalCool_Hef + WalCool_n_Hef;

  WalCool_enHS = WalCool_Ab  + WalCool_n_Ab;
  WalCool_enH  = WalCool_enHS  + 2 * WalCool_arraysS_size;
  WalCool_mu   = WalCool_enH + 2 * WalCool_arrays_size;
  WalCool_UtoT = WalCool_mu  + 2 * WalCool_arrays_size;

                                                                         /* get the arrays */
  if(ThisTask == 0)
    {
                                                                         /* - solar abundances by mass */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Header/Abundances/Solar_mass_fractions");
#else
      dataset     = H5Dopen(fid, "/Header/Abundances/Solar_mass_fractions", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Ab);
      hd_report   = H5Dclose(dataset);


                                                                         /* - elecron density over n_h */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enH);
      hd_report   = H5Dclose(dataset);
      

                                                                         /* - helium fraction bins */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Hef);
      hd_report   = H5Dclose(dataset);


                                                                         /* - density bins */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Hydrogen_density_bins");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Hydrogen_density_bins", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Rho);
      hd_report   = H5Dclose(dataset);

      
                                                                         /* - mean particle mass */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_mu);
      hd_report   = H5Dclose(dataset);


                                                                         /* - energy density bins */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Temperature/Energy_density_bins");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Temperature/Energy_density_bins", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_U);
      hd_report   = H5Dclose(dataset);

      
                                                                         /* - temperature bins */
#ifdef OLD_HDF5
      dataset     = H5Dopen(fid, "/Metal_free/Temperature_bins");
#else
      dataset     = H5Dopen(fid, "/Metal_free/Temperature_bins", H5P_DEFAULT);
#endif
      hd_report   = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_T);
      hd_report   = H5Dclose(dataset);

      hd_report = H5Fclose(fid);
    }

  MPI_Bcast(&WalCool_indexes_and_arrays[0], size * sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);

  Trange = log10(WalCool_T[WalCool_n_T - 1] / WalCool_T[0]);
  Tmin   = log10(WalCool_T[0]);
  Urange = log10(WalCool_U[WalCool_n_T - 1] / WalCool_U[0]);
  Umin   = log10(WalCool_U[0]);
  Rrange = log10(WalCool_Rho[WalCool_n_Rho - 1] / WalCool_Rho[0]);
  Rmin   = log10(WalCool_Rho[0]);
  Herange= WalCool_Hef[WalCool_n_Hef - 1] - WalCool_Hef[0];
  Hemin  = WalCool_Hef[0];

                                                                         /* allocate memory for names */
  
  WalCool_names_tab = (char**)mymalloc("WalCool_names_tab", (WalCool_n_El + 2 * WalCool_n_Ab) * sizeof(char*));
  WalCool_names     = (char*) mymalloc("WalCool_names",     (WalCool_n_El + 2 * WalCool_n_Ab) * NAME_LEN * sizeof(char));

  WalCool_El_names = WalCool_names_tab;
  for (pos = i = 0; i < WalCool_n_El; i++, pos += NAME_LEN)
    WalCool_El_names[i] = &WalCool_names[pos];

  k = pos;
  WalCool_Ab_names = WalCool_El_names + WalCool_n_El;
  for (i = 0; i < WalCool_n_Ab; i++, pos += NAME_LEN)
    WalCool_Ab_names[i] = &WalCool_names[pos];

  k = pos;
  WalCool_Ab_symbols = WalCool_Ab_names + WalCool_n_Ab;
  for (i = 0; i < WalCool_n_Ab; i++, pos += NAME_LEN)
    WalCool_Ab_symbols[i] = &WalCool_names[pos];

                                                                         /* get the names */
  if(ThisTask == 0)
    {
      char line[200];      
      FILE *myfile;

      sprintf(filename, "%s/Species_names.txt", All.WalCool_CoolTables_path);
      if((myfile = fopen(filename, "r")) == NULL)
	{
	  printf("Can't open file %s\n",filename);
	  endrun(225533);
	}
      for(i = 0; i < WalCool_n_El; i++)
        {
          fgets(WalCool_El_names[i], NAME_LEN, myfile);
          if(WalCool_El_names[i][strlen(WalCool_El_names[i]) - 1] == '\n')
            WalCool_El_names[i][strlen(WalCool_El_names[i]) - 1] = '\0';
        }

      fclose(myfile);
      sprintf(filename, "%s/Abundances_names.txt", All.WalCool_CoolTables_path);
      if((myfile = fopen(filename, "r")) == NULL)
	{
	  printf("Can't open file %s\n",filename);
	  endrun(225544);
	}
      for(i = 0; i < WalCool_n_Ab; i++)
        {
	  fgets(line, 200, myfile);
	  if(sscanf(line, "%s %s\n", WalCool_Ab_names[i], WalCool_Ab_symbols[i]) < 2)
	    {
	      printf("I've got some problem in reading %s/Abundances_names.txt @ line %d : %s\n", All.WalCool_CoolTables_path, i, line);
	      endrun(LT_ERR_WALCOOL_IMPOSSIBLE_TO_READ_ABUNDANCE);
	    }
          if(WalCool_Ab_names[i][strlen(WalCool_Ab_names[i]) - 1] == '\n')
            WalCool_Ab_names[i][strlen(WalCool_Ab_names[i]) - 1] = '\0';
          if(WalCool_Ab_symbols[i][strlen(WalCool_Ab_symbols[i]) - 1] == '\n')
            WalCool_Ab_names[i][strlen(WalCool_Ab_symbols[i]) - 1] = '\0';
        }
      fclose(myfile);
    }

  MPI_Bcast(WalCool_names, (WalCool_n_El + 2 * WalCool_n_Ab) * NAME_LEN, MPI_BYTE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef LT_STELLAREVOLUTION

  SpeciesIdx = (int*)mymalloc("SpeciesIdx", LT_NMet * sizeof(int));      /* defines which position a coolant has in the normal chemical array */
  Specie_is_coolant = (int*)mymalloc("Specie_is_coolant", LT_NMet * sizeof(int)); /* defined whether a specie is a coolant */
  SpeciesPos = (int*)mymalloc("SpeciesPos", WalCool_n_Ab * sizeof(int)); /* defines the ordinal position in the tables for a coolant */

  
  for(i = 0; i < LT_NMet; i++)
    SpeciesIdx[i] = Specie_is_coolant[i] = -1;

  for(pos = 0, i = 0; i < WalCool_n_Ab; i++)
    {
      if(strcmp(WalCool_Ab_names[i], "Hydrogen") == 0)
        {
          myHyd = i;
          HydPos = WalCool_n_El;
          SpeciesPos[i] = HydPos;
	  SpeciesIdx[i] = Hyd;
        }
      else if(strcmp(WalCool_Ab_names[i], "Helium") == 0)
        {
          myHel = i;
          HelPos = WalCool_n_El + 1;
          SpeciesPos[i] = HelPos;
	  SpeciesIdx[i] = Hel;
          Specie_is_coolant[Hel] = SpeciesPos[i];
        }
      else
	{
	  for(k = 0; k < LT_NMet; k++)
	    if(strcmp(MetNames[k], WalCool_Ab_symbols[i]) == 0)
	      break;
	  if(k == LT_NMet)
            {
              if(ThisTask == 0)
                printf("error initializing cooling: %s (%s) is absent in chemical evolution\n", WalCool_Ab_names[i], WalCool_Ab_symbols[i]); fflush(stdout);
              endrun(LT_ERR_WALCOOL_EL_MISSED);
            }
	  SpeciesIdx[i] = k;
	  SpeciesPos[i] = pos++;
          Specie_is_coolant[k] = SpeciesPos[i];
	}

    }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return;
}


int compare_dir_redshifts(const void *A, const void *B)
/* used to sort the tables' redshifts */
{
  if( (*(struct my_direntry*)A).redshift > (*(struct my_direntry*)B).redshift)
    return 1;
  if( (*(struct my_direntry*)A).redshift < (*(struct my_direntry*)B).redshift)
    return -1;
  return 0;
}


int compare_redshifts(const void *A, const void *B)
/* used to sort the tables' redshifts */
{
  if( (*(float*)A) > (*(struct my_direntry*)B).redshift)
    return 1;
  if( (*(float*)A) < (*(struct my_direntry*)B).redshift)
    return -1;
  return 0;
}


int INLINE_FUNC find_index(int mode, float *base, double min, double range, int N, float value, double *d)
/* mode = LOG means logarithmic spacing
   mode = LIN means linear spacing
   where LOG and LIN are defined in the .data segment
*/
{
  int idx;

  
  if(value > base[N - 1])
    {
      *d = 0;
      return N - 2;
    }
  else if(value <= base[0])
    {
      *d = 1;
      return 0;
    }
                                                                         /* guess */
  
  if(mode == LOG)  
    idx = (int)(log10(value) - min) / range * N;
  else
    idx = (int)(value - min) / range * N;

  if(idx < 0)
    idx = 0;
  if(idx > N - 1)
    idx = N -1;
                                                                         /* check the guess */
  for(; value > base[idx] && idx < N - 1; idx++);

  for(; value < base[idx] && idx > 0; idx--);
   
  if(mode == LOG)                                                        /* note: no pathological behaviour should be possible here */
    *d = log10(value / base[idx]) / log10(base[idx + 1] / base[idx]);        /* Hence, we avoid checks in order to speed-up             */
  else
    *d = (value - base[idx]) / (base[idx + 1] - base[idx]); 

  return idx;
}


void WalCool_set_PID(MyIDType i)
{
  PID = i;
  return;
}

/* double binterp(float x1, float x2, float y1, float y2, float vx1y1, float vx2y1, float vx1y2, float vx2y2, float x, float y) */
/* /\* calculate the bilinear interpolation on a square like the one here below *\/ */
/* /\* */
/*         ^ */
/*         |   (Vx1y2)          (Vx2y2) */
/*     y2  |     +---------------+ */
/*         |     |               | */
/*         |     |               | */
/*         |     |               | */
/*         |     |        (x,y)  | */
/*         |     |..........+    | */
/*         |     |          :    | */
/*         |     |          :    | */
/*     y1  |     +---------------+ */
/*         |   (Vx1y1)          (Vx2y1) */
/*         | */
/*         +---------------------------------->     */
/*              x1               x2 */

/*  *\/ */
/* { */
/*   double value; */
/*   double nx, ny; */

/*   if(x1 == x2) */
/*     nx = 0; */
/*   else */
/*     nx = (x - x1) / (x2 - x1); */

/*   if(y1 == y2) */
/*     ny = 0; */
/*   else */
/*     ny = (y - y1) / (y2 - y1); */

/*   value = (vx1y1 * (1 - nx)*(1 - ny) + */
/*            vx2y1 * nx * (1 - ny) + */
/*            vx1y2 * ny * (1 - nx) + */
/*            vx2y2 * nx * ny); */
  
/*   return value; */
/* } */


/* int INLINE_FUNC find_1D_index_in_ND_matrix(int N, int *L, int *gridpoints) */
/* { */
/*   int index, i, partial; */
  
/*   /\* find the index in the uni-dimensional array that represents the */
/*      N-dimensional matrix */

/*      L is an array that stores the length of each dimension of the */
/*      matrix; */
/*      gridpoint is an array that specifies, in N-dimensional coordinates. */
/*      the grid point whose 1-dimensional index we want to find out. */
/*    *\/ */
/*   for(index = 0, i = N-1; i >= 0; i--) */
/*     { */
/*       for(partial = gridpoint[i], k = i-1; k >= 0; k--) */
/*         partial *= L[k]; */
/*       index += partial; */
/*     } */
  
/*   return index; */
/* } */



/* float get_Ndim_interp(int N, double *P, double *Dx) */
/* /\* */
/*   This routine gives the linear interpolation in a N-dimensional */
/*   function defined on a grid. */

/*   --==[ INPUT VALUES ]==-- */
  
/*   N  is the number of dimensions */
/*   P  is a 2^N long array that stores the values of the function */
/*      at the grid points that encompasse the interpolation point */
/*   Dx is an N-dimensional array that stores the normalized */
/*      displacement of the interpolation point with respect to the */
/*      lower grid point: */
/*         Dx_i =  (x_i - x0_i) / (x1_i - x0_i) */

/*         in two dimensions, for instance: */

/*         ^ */
/*         |                             */
/*    x1_1 |     +---------------+ */
/*         |     |               | */
/*         |     |               | */
/*         |     |               | */
/*         |     |      (x_0,x_1)| */
/*         |     |..........+    | */
/*         |     |          :    | */
/*         |     |          :    | */
/*    x0_1 |     +---------------+ */
/*         |                            */
/*         | */
/*         +---------------------------------->     */
/*              x0_0             x1_0 */


/*   --==[ HOW TO FILL THE ARRAYS ]==-- */

/*   As first you must define what is the order of the variables x,i */
/*   (that is irrelevant: interpolation is commutative ini this respect). */

/*   Then you fill the Dx array with that order. */
  
/*   The P array is filled in the same "walking order" of the N-dimensional */
/*   cube - which has 2^N verteces - that encloses the interpolation point. */
             
/*  *\/ */
/* { */
/*   double INTERP = 0, partial; */
/*   int    i, j; */
  
/*   for(j = 0; j < (1<<N); j++) */
/*     { */
/*       for(i = 0; i < N; i++) */
/*         S[i] = j & (1 << i); */

/*       for(partial = 1, i = 0; (i < N) && (partial> 0); i++) */
/*         { */
/*           if(S[i] == 0) */
/*             partial *= (1.0 - Dx[i]); */
/*           else */
/*             partial *= Dx[i]; */
          
/*           partial *= P[i]; */
/*         } */
      
/*       INTERP += partial; */
/*     } */
  
/*   return (float)INTERP; */


/*   /\* */
/*     The above routine is the algorithmical translation of the following */
/*     formula: */
/*                                                _                       _ */
/*                                  ----         /     _____               \                       */
/*                                  \            |      | |      [m,j]     |                       */
/*      interpolated value =        /            |      | |     U      [X] | F(m1 + s1,..,mn+sn)   */
/*                                  ----         \_ j in {1,..n}  s_j     _/                       */
/*                             s1,..,sn in 0,1} */

/*     where: */

/*     n are the dimensions */
/*     m1,..,mn are the grid points */
/*     x1,..,xn are the coordinates of the interpolation point */
/*     s1,..,sn are either 0 or 1 and means the leap from a grid point to */
/*              the subsequent one along the same coordinate */

/*      [m,j] */
/*     U     [X] = (1 - s_j) + (x_j - m_j)(2s_j - 1) */
/*      s_j */

                      
/*    *\/ */
/* } */
#endif
