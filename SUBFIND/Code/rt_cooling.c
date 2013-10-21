#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include "cooling.h"

#if defined(RT_COOLING_PHOTOHEATING)

/* rate1 : photoheating for a blackbody spectrum */
/* rate2 : recombination cooling rate */
/* rate3 : collisional ionization cooling rate */
/* rate4 : collisional excitation cooling rate */
/* rate5 : Bremsstrahlung cooling rate */

/* now do the heating (note: we now how many photons we absorbed) */

double rt_DoHeating(int i, double dt_internal)
{
  int j;
  double sigma, a3inv, nH;
  double hubble_a, rate, du, de, c_light;
  double n_gamma, nHI;
  double E;
#ifdef RT_INCLUDE_HE
  double x, y, dE;
  double nHe, nHeI, nHeII;
#endif
  
  if(All.ComovingIntegrationOn)
    {
      a3inv = 1.0 / All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  c_light = C / All.UnitVelocity_in_cm_per_s;
  nH = (HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);
  nHI = SphP[i].nHI * nH;

#ifdef RT_INCLUDE_HE
  nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].d.Density * a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam);  
  nHeI = SphP[i].nHeI * nHe;
  nHeII = SphP[i].nHeII * nHe;
  dE = (end_E - start_E) / (float)N_BINS;
#endif      

  for(j = 0, rate = 0; j < N_BINS; j++)
    {
#ifdef RT_INCLUDE_HE
      x = SphP[i].nHI * nH * rt_sigma_HI[j] / 
        (SphP[i].nHI * nH * rt_sigma_HI[j] + SphP[i].nHeI * nHe * rt_sigma_HeI[j] + SphP[i].nHeII * nHe * rt_sigma_HeII[j]);
      
      y = SphP[i].nHeI * nHe * rt_sigma_HeI[j] / 
	(SphP[i].nHI * nH * rt_sigma_HI[j] + SphP[i].nHeI * nHe * rt_sigma_HeI[j] + SphP[i].nHeII * nHe * rt_sigma_HeII[j]);
      
      E = start_E + (j + 0.5) * dE * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
      
      sigma = rt_sigma_HI[j];
      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv * x;
      rate += nHI * c_light * E * sigma * n_gamma;

      sigma = rt_sigma_HeI[j];
      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv * y;
      rate += nHeI * c_light * E * sigma * n_gamma;

      sigma = rt_sigma_HeII[j];
      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv * (1.0 - x - y);
      rate += nHeII * c_light * E * sigma * n_gamma;
#else
      sigma = 1.49e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv;
      E = 6.4 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
      rate += nHI * c_light * E * sigma * n_gamma;
#endif
    }

  du = rate * dt_internal / hubble_a / (SphP[i].d.Density * a3inv);
  de = du * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
  
  return de / dt_internal;
}

double rt_DoCooling(int i, double dt_internal)
{
  double dtime, a3inv;
  
  if(All.ComovingIntegrationOn)
    {
      dtime = dt_internal / hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    {
      dtime = dt_internal;
      a3inv = 1.0;
    }

  /* do the cooling */
  double lambda = rt_get_cooling_rate(i, SphP[i].Entropy + SphP[i].e.DtEntropy * dt_internal);
  double du = lambda * dtime / (SphP[i].d.Density * a3inv);
  double de = du * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

  if(fabs(de) < 0.2 * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_internal))
    {
      /* cooling is slow, we can do it explicitly */
      return de / dt_internal; 
    }
  else
    {
      /* rapid cooling. Better calculate an implicit solution, which is determined by bisection */
      
      double u_old = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_internal) / 
	GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
      double u_lower = u_old / sqrt(1.1);
      double u_upper = u_old * sqrt(1.1);
      double ratefact = dtime / (SphP[i].d.Density * a3inv);
      int iter = 0;
      
      /* bracketing */
      while(u_lower - u_old - ratefact * 
	    rt_get_cooling_rate(i, u_lower * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1)) > 0)
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
	  
	  if(iter++ >= 1000) //MAXITER)
	    terminate("bracketing failure");
	}
      
      /* bisection */
      double u;
      iter = 0;
      do
	{
	  u = 0.5 * (u_lower + u_upper);
	  
	  if(u - u_old - ratefact * rt_get_cooling_rate(i, u * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1)) > 0)
	    u_upper = u;
	  else
	    u_lower = u;
	  
	  du = u_upper - u_lower;
	  
	  iter++;
	  
	  if(iter >= (MAXITER - 10))
	    printf("u= %g\n", u);
	  
	  if(iter >= MAXITER)
	    terminate("convergence failure");
	}
      while(fabs(du / u) > 1.0e-6);
      
      du = u - u_old;
      
      return du * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / dt_internal;
    }

}
/* returns cooling rate */
double rt_get_cooling_rate(int i, double entropy)
{
  double Lambda;
  double temp, molecular_weight;
  double a3inv;
  double nH;
  double rate2, rate3, rate4, rate5;
  double de2, de3, de4, de5;
#ifdef RT_INCLUDE_HE
  double rateHe1, rateHe2, rateHe3, rateHe4, rateHe5;
  double nHe, deHe1, deHe2, deHe3, deHe4, deHe5;
#endif
  
  if(All.ComovingIntegrationOn)
    a3inv = 1 / All.Time / All.Time / All.Time;
  else
    a3inv = 1;

  nH = (HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);	//physical
  molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].n_elec);

  temp = entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) *
    (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) /
    (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

  /* all rates in erg cm^3 s^-1 in code units */
  /* recombination cooling rate */
  rate2 = 8.7e-27 * pow(temp, 0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
  rate2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate2 *= All.HubbleParam * All.HubbleParam;
  rate2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de2 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate2;

  /* collisional ionization cooling rate */
  rate3 = 1.27e-21 * pow(temp, 0.5) * exp(-157809.1 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rate3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate3 *= All.HubbleParam * All.HubbleParam;
  rate3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de3 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate3;

  /* collisional excitation cooling rate */
  rate4 = 7.5e-19 / (1.0 + pow(temp / 1e5, 0.5)) * exp(-118348 / temp);
  rate4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate4 *= All.HubbleParam * All.HubbleParam;
  rate4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de4 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate4;

  /* Bremsstrahlung cooling rate */
  rate5 = 1.42e-27 * pow(temp, 0.5);
  rate5 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate5 *= All.HubbleParam * All.HubbleParam;
  rate5 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de5 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate5;

  Lambda = de2 + de3 + de4 + de5;

#ifdef RT_INCLUDE_HE
  nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].d.Density * a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam);  

  /* recombination cooling rate */
  rateHe2 = 1.55e-26 * pow(temp, 0.3647);
  rateHe2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe2 *= All.HubbleParam * All.HubbleParam;
  rateHe2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe2 = SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe2;

  rateHe2 = 3.48e-26 * pow(temp, 0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
  rateHe2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe2 *= All.HubbleParam * All.HubbleParam;
  rateHe2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe2 += SphP[i].nHeIII * nHe * SphP[i].n_elec * nH * rateHe2;

  /* collisional ionization cooling rate */
  rateHe3 = 9.38e-22 * pow(temp, 0.5) * exp(-285335.4 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rateHe3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe3 *= All.HubbleParam * All.HubbleParam;
  rateHe3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe3 = SphP[i].nHeI * nHe * SphP[i].n_elec * nH * rateHe3;

  rateHe3 = 4.95e-22 * pow(temp, 0.5) * exp(-631515 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rateHe3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe3 *= All.HubbleParam * All.HubbleParam;
  rateHe3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe3 += SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe3;

  /* collisional excitation cooling rate */
  rateHe4 = 5.54e-17 * pow(temp, -0.397) / (1.0 + pow(temp / 1e5, 0.5)) * exp(-473638 / temp);
  rateHe4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= All.HubbleParam * All.HubbleParam;
  rateHe4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe4 = SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe4;

  rateHe4 = 9.10e-27 * pow(temp, -0.1687) / (1.0 + pow(temp / 1e5, 0.5)) * exp(-13179 / temp);
  rateHe4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= All.HubbleParam * All.HubbleParam * All.HubbleParam * All.HubbleParam * All.HubbleParam;
  rateHe4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe4 += SphP[i].nHeII * nHe * SphP[i].n_elec * nH * SphP[i].n_elec * nH * rateHe4;

  /* Bremsstrahlung cooling rate */
  rateHe5 = 1.42e-27 * pow(temp, 0.5);
  rateHe5 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe5 *= All.HubbleParam * All.HubbleParam;
  rateHe5 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe5 = (SphP[i].nHeII * nHe * SphP[i].n_elec * nH + 4.0 * SphP[i].nHeIII * nHe * SphP[i].n_elec * nH) * rateHe5;

  Lambda += deHe2 + deHe3 + deHe4 + deHe5;
#endif
  
  return - Lambda;
}



#endif
