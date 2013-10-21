#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef RADTRANSFER

void radtransfer_update_chemistry(void)
{
  int i, j;
  double nH, temp, molecular_weight;
  double nHII;
  double dt, dtime, a3inv, c_light;
  double A, B, CC;
  double x;
  double n_gamma;
  double alpha_HII, gamma_HI; 
  double total_nHI, total_V, total_nHI_all, total_V_all;

#ifdef RT_INCLUDE_HE
  double nHe, alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHe, nHeII, nHeIII;
  double y, D, E, F, G, J;
  double total_nHeI, total_nHeI_all;

  total_nHeI = 0;
#endif
      
  total_nHI = total_V = 0;

  dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;
  
  if(All.ComovingIntegrationOn)
    {
      dtime = dt / hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    {
      dtime = dt;
      a3inv = 1.0;
    }

  c_light = C / All.UnitVelocity_in_cm_per_s;

  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
	nH = (HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);
	
	molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].n_elec);
	
	temp = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) *
	  (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) /
	  (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

	/* collisional ionization rate */
	gamma_HI = 5.85e-11 * pow(temp, 0.5) * exp(-157809.1 / temp) / (1.0 + pow(temp / 1e5, 0.5));
	gamma_HI *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	gamma_HI *= All.HubbleParam * All.HubbleParam;
	
	/* alpha_B recombination coefficient */
	alpha_HII = 2.59e-13 * pow(temp / 1e4, -0.7);
	alpha_HII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	alpha_HII *= All.HubbleParam * All.HubbleParam;
	
#ifdef RT_INCLUDE_HE
	nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].d.Density * a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam);   

	/* collisional ionization rate */
	gamma_HeI = 2.38e-11 * pow(temp, 0.5) * exp(-285335.4 / temp) / (1.0 + pow(temp / 1e5, 0.5));
	gamma_HeI *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	gamma_HeI *= All.HubbleParam * All.HubbleParam;
	
	gamma_HeII = 5.68e-12 * pow(temp, 0.5) * exp(-631515 / temp) / (1.0 + pow(temp / 1e5, 0.5));
	gamma_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	gamma_HeII *= All.HubbleParam * All.HubbleParam;
	
	/* alpha_B recombination coefficient */
	alpha_HeII = 1.5e-10 * pow(temp, -0.6353);
	alpha_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	alpha_HeII *= All.HubbleParam * All.HubbleParam;
	
	alpha_HeIII = 3.36e-10 * pow(temp, -0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
	alpha_HeIII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	alpha_HeIII *= All.HubbleParam * All.HubbleParam;
#endif
	
	for(j = 0; j < N_BINS; j++)
	  {
#ifdef RT_INCLUDE_HE
	    x = SphP[i].nHI * nH * rt_sigma_HI[j] / 
	      (SphP[i].nHI * nH * rt_sigma_HI[j] + SphP[i].nHeI * nHe * rt_sigma_HeI[j] + SphP[i].nHeII * nHe * rt_sigma_HeII[j]);
#else
	    x = 1.0;
#endif

	    n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv;

	    /* number of photons should be positive */
	    if(n_gamma < 0 || isnan(n_gamma))
	      {
		printf("NEGATIVE n_gamma: %g %d %d \n", n_gamma, i, ThisTask);
		printf("n_gamma[j] %g mass %g a3inv %g \n", SphP[i].n_gamma[j], P[i].Mass, a3inv);
		endrun(111);
	      }

	    A = dtime * gamma_HI * nH;
	    B = dtime * c_light * rt_sigma_HI[j] * x;
	    CC = dtime * alpha_HII* nH;
	    
	    /* semi-implicit scheme for ionization */
	    nHII =  (SphP[i].nHII + B * n_gamma + A * SphP[i].n_elec) /
	      (1.0 + B * n_gamma + CC * SphP[i].n_elec + A * SphP[i].n_elec);

	    if(nHII < 0 || nHII > 1 || isnan(nHII))
	      {
		printf("ERROR nHII %g \n", nHII);
		endrun(333);
	      }
	    
	    SphP[i].n_elec = nHII;
	    
	    SphP[i].nHII = nHII;
	    
	    SphP[i].nHI = 1.0 - nHII;
	    

#ifdef RT_INCLUDE_HE
	    SphP[i].n_elec += SphP[i].nHeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
	    SphP[i].n_elec += 2.0 * SphP[i].nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

	    y = SphP[i].nHeI * nHe * rt_sigma_HeI[j] / 
	      (SphP[i].nHI * nH * rt_sigma_HI[j] + SphP[i].nHeI * nHe * rt_sigma_HeI[j] + SphP[i].nHeII * nHe * rt_sigma_HeII[j]);
	    
	    D = dtime * gamma_HeII * nH;
	    E = dtime * alpha_HeIII * nH;
	    F = dtime * gamma_HeI * nH;
	    G = dtime * c_light * rt_sigma_HeI[j] * y;
	    J = dtime * alpha_HeII * nH;
	    
	    nHeII =
	      SphP[i].nHeII + F * SphP[i].n_elec + G * n_gamma -
	      ((F * SphP[i].n_elec + G * n_gamma - E * SphP[i].n_elec) / (1.0 +
									     E * SphP[i].n_elec) *
	       SphP[i].nHeIII);
	    
	    nHeII /= 1.0 + F * SphP[i].n_elec + G * n_gamma + D * SphP[i].n_elec + J * SphP[i].n_elec +
	      ((F * SphP[i].n_elec + G * n_gamma - E * SphP[i].n_elec) / (1.0 +
									     E * SphP[i].n_elec) * D *
	       SphP[i].n_elec);
	    
	    if(nHeII < 0 || nHeII > 1 || isnan(nHeII))
	      {
		printf("ERROR neHII %g\n", nHeII);
		endrun(333);
	      }
	    
	    nHeIII = (SphP[i].nHeIII + D * nHeII * SphP[i].n_elec) / (1.0 + E * SphP[i].n_elec);
	    
	    if(nHeIII < 0 || nHeIII > 1)
	      {
		printf("ERROR nHeIII %g\n", nHeIII);
		endrun(222);
	      }
	    
	    SphP[i].n_elec = SphP[i].nHII;
	    SphP[i].n_elec += nHeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
	    SphP[i].n_elec += 2.0 * nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
	    
	    SphP[i].nHeII = nHeII;
	    SphP[i].nHeIII = nHeIII;
	    
	    SphP[i].nHeI = 1.0 - SphP[i].nHeII - SphP[i].nHeIII;
	    
	    if(SphP[i].nHeI < 0 || SphP[i].nHeI > 1 || isnan(SphP[i].nHeI))
	      {
		printf("ERROR nHeI %g\n", SphP[i].nHeI);
		endrun(444);
	      }
#endif
	  }	    


#ifdef RT_INCLUDE_HE
	total_nHeI += SphP[i].nHeI * P[i].Mass / (SphP[i].d.Density * a3inv);
#endif
	total_nHI += SphP[i].nHI * P[i].Mass / (SphP[i].d.Density * a3inv);
	total_V += P[i].Mass / (SphP[i].d.Density * a3inv);
	
      }
  
  MPI_Allreduce(&total_nHI, &total_nHI_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_V, &total_V_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifdef RT_INCLUDE_HE
  MPI_Allreduce(&total_nHeI, &total_nHeI_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif  

  /* output the input of photon density in physical units */
  if(ThisTask == 0)
    {
      fprintf(FdRad, "%g %g ", All.Time, total_nHI_all / total_V_all);
#ifdef RT_INCLUDE_HE
      fprintf(FdRad, "%g\n", total_nHeI_all / total_V_all);
#else
      fprintf(FdRad, "\n");
#endif
      fflush(FdRad);
    }
}

#endif  
