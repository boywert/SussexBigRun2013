#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

#ifndef LT_STELLAREVOLUTION
#ifdef COOLING

/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */

#ifndef MHM
#ifndef SFR

void cooling_only(void)		/* normal cooling routine when star formation is disabled */
{
  int i;
  double dt, dtime, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, dmax1, dmax2;

#ifdef COSMIC_RAYS
  int CRpop;
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    {
      a3inv = time_hubble_a = hubble_a = 1;
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

	  double fac_entr_to_u = pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	  double uold =  DMAX(All.MinEgySpec, SphP[i].Entropy * fac_entr_to_u);

#ifndef RT_COOLING_PHOTOHEATING
	  unew = DoCooling(uold, SphP[i].d.Density * a3inv, dtime, &ne);
#else
          unew = uold  + dt * fac_entr_to_u * ( rt_DoHeating(i,dt)  + rt_DoCooling(i,dt) );
#endif


#ifdef COSMIC_RAYS
          for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
             unew += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime, CRpop);
#endif

	  SphP[i].Entropy  = unew / fac_entr_to_u;
	  SphP[i].Ne = ne;

          SphP[i].EntropyPred = SphP[i].Entropy;
          SphP[i].Pressure = get_pressure(i);
	}
    }
}


#else

void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, mass_of_star;
  double sum_sm, total_sm, sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate;

#ifdef WINDS
  int j;
  double v;
  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif

#ifdef VARIABLE_WINDS
  double z, r200c, v_esc, c_halo, wind_energy, wind_momentum, wind_mass;
  double rhocrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  
  if(All.ComovingIntegrationOn)
    z = 1 / All.Time - 1;
  else
    z = 0;

  rhocrit *= (All.Omega0 * pow(1 + z, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + z,2) + All.OmegaLambda); /* physical critical density at redshift z */
#endif
#endif

#ifdef METALS
  double w;
#endif
#ifdef COSMIC_RAYS
  int CRpop;

#ifdef CR_SN_INJECTION
  double tinj = 0.0, instant_reheat = 0.0;
  int InjPopulation;
#endif
#endif


#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK) || defined(MODIFIED_SFR)
  double temp, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

#ifdef MODIFIED_SFR

  double SFRTempThresh;

  SFRTempThresh = 5.0e5 / u_to_temp_fac;

#endif


#ifdef FLTROUNDOFFREDUCTION
#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      SphP[i].i.Injected_BH_Energy = FLT(SphP[i].i.dInjected_BH_Energy);
#endif
#endif

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    a3inv = ascale = time_hubble_a = 1;



  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

  for(bits = 0; GENERATIONS > (1 << bits); bits++);

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  /* check whether conditions for star formation are fulfilled.
	   *  
	   * f=1  normal cooling
	   * f=0  star formation
	   */
	  flag = 1;		/* default is normal cooling */

#ifndef MODIFIED_SFR
	  if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
	    flag = 0;
#else
	  if((SphP[i].d.Density * a3inv >= All.PhysDensThresh)
	     && (SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1 <
		 SFRTempThresh))
	    flag = 0;
#endif

	  if(All.ComovingIntegrationOn)
	    if(SphP[i].d.Density < All.OverDensThresh)
	      flag = 1;

#ifdef BLACK_HOLES
	  if(P[i].Mass == 0)
	    flag = 1;
#endif

#ifdef WINDS
	  if(SphP[i].DelayTime > 0)
	    flag = 1;		/* only normal cooling for particles in the wind */

	  if(SphP[i].DelayTime > 0)
	    SphP[i].DelayTime -= dtime;

	  if(SphP[i].DelayTime > 0)
	    if(SphP[i].d.Density * a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh)
	      SphP[i].DelayTime = 0;

	  if(SphP[i].DelayTime < 0)
	    SphP[i].DelayTime = 0;

#endif


#ifdef QUICK_LYALPHA
	  temp = u_to_temp_fac * (SphP[i].Entropy) /
	    GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	  if(SphP[i].d.Density > All.OverDensThresh && temp < 1.0e5)
	    flag = 0;
	  else
	    flag = 1;
#endif
#ifdef MAGNETIC
	  x=0.;
#endif


#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	  if(flag == 1)		/* normal implicit isochoric cooling */
#endif
	    {
	      SphP[i].Sfr = 0;
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_INJECTION)
	      SphP[i].CR_Specific_SupernovaHeatingRate = 0;
#endif
	      ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

	      double fac_entr_to_u = pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	      double uold =  DMAX(All.MinEgySpec, SphP[i].Entropy * fac_entr_to_u);


#ifndef RT_COOLING_PHOTOHEATING
	      unew = DoCooling(uold, SphP[i].d.Density * a3inv, dtime, &ne);
#else
              unew = uold  + dt * fac_entr_to_u * ( rt_DoHeating(i,dt)  + rt_DoCooling(i,dt) );
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
	      if(SphP[i].i.Injected_BH_Energy)
		{
		  if(P[i].Mass == 0)
		    SphP[i].i.Injected_BH_Energy = 0;
		  else
		    unew += SphP[i].i.Injected_BH_Energy / P[i].Mass;

		  temp = u_to_temp_fac * unew;


		  if(temp > 5.0e9)
		    unew = 5.0e9 / u_to_temp_fac;

#ifdef FLTROUNDOFFREDUCTION
		  SphP[i].i.dInjected_BH_Energy = 0;
#else
		  SphP[i].i.Injected_BH_Energy = 0;
#endif
		}
#endif

#ifdef COSMIC_RAYS
	      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	        unew += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime, CRpop);
#endif

	      SphP[i].Entropy  = unew / fac_entr_to_u;
	      SphP[i].Ne = ne;

              SphP[i].EntropyPred = SphP[i].Entropy;
              SphP[i].Pressure = get_pressure(i);
	    }

	  if(flag == 0)		/* active star formation */
	    {
#if !defined(QUICK_LYALPHA)
	      tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

	      factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = SphP[i].Ne;

	      tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);

	      SphP[i].Ne = ne;

	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      cloudmass = x * P[i].Mass;

	      if(tsfr < dtime)
		tsfr = dtime;

	      sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */

	      p = sm / P[i].Mass;

	      sum_sm += P[i].Mass * (1 - exp(-p));


	      if(dt > 0)
		{
		  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		    {
		      trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
		      egycurrent =
			SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;

#ifdef COSMIC_RAYS
#ifdef CR_SN_INJECTION
		      if(All.CR_SNEff > 0)
			{
			  if(NUMCRPOP > 1)
			    InjPopulation = CR_Find_Alpha_to_InjectTo(All.CR_SNAlpha);
			  else
			    InjPopulation = 0;

			  tinj =
			    SphP[i].CR_E0[InjPopulation] / (p * All.FeedbackEnergy * All.CR_SNEff / dtime);

			  instant_reheat =
			    CR_Particle_SupernovaFeedback(&SphP[i], p * All.FeedbackEnergy * All.CR_SNEff,
							  tinj);
			}
		      else
			instant_reheat = 0;

#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_INJECTION)
		      SphP[i].CR_Specific_SupernovaHeatingRate =
			(p * All.FeedbackEnergy * All.CR_SNEff - instant_reheat) / dtime;
#endif
		      egycurrent += instant_reheat;
#endif
		      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
			egycurrent += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime, CRpop);
#endif /* COSMIC_RAYS */


#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
		      if(SphP[i].i.Injected_BH_Energy > 0)
			{
			  egycurrent += SphP[i].i.Injected_BH_Energy / P[i].Mass;

			  temp = u_to_temp_fac * egycurrent;

			  if(temp > 5.0e9)
			    egycurrent = 5.0e9 / u_to_temp_fac;

			  if(egycurrent > egyeff)
			    {
			      tcool = GetCoolingTime(egycurrent, SphP[i].d.Density * a3inv, &ne);

			      if(tcool < trelax && tcool > 0)
				trelax = tcool;
			    }

			  SphP[i].i.Injected_BH_Energy = 0;
			}
#endif
#ifdef MAGNETICSEED
	             SphP[i].MagSeed =  egyhot * factorEVP *  1E-2;   // This is the definition of how much energy we will put in MF (we neglect cooling here, we have to check if thios has sense)
		     egyeff-= SphP[i].MagSeed;                 //Here we also substract that to the feedback

	             SphP[i].MagSeed *= All.UnitMass_in_g / All.UnitEnergy_in_cgs * SphP[i].d.Density * (1.-x) * a3inv *All.Time ;// * a3inv
#endif



#if !defined(NOISMPRESSURE)

#ifdef ALTERNATIVE_SFR
		      double dens=SphP[i].d.Density;
		      if(egycurrent > egyeff)
			{
			  ne = SphP[i].Ne;
			  unew = DoCooling(egycurrent, dens * a3inv, dtime, &ne);

			  if(unew < egyeff)
			    {
#ifdef SLOW_RELAX_TO_EOS
			      tcool = GetCoolingTime(egycurrent, dens * a3inv, &ne);
			      if(tcool < trelax && tcool > 0)
				trelax = tcool;

			      unew = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
#else
			      unew = egyeff;
#endif
			    }
			}
		      else
			unew = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));

                      SphP[i].Entropy = unew  * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
#else                      
		      SphP[i].Entropy =
			(egyeff +
			 (egycurrent -
			  egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
#endif
		      SphP[i].e.DtEntropy = 0;
#endif

		      SphP[i].EntropyPred = SphP[i].Entropy;
		      SphP[i].Pressure = get_pressure(i);
		    }
		}



	      /* the upper bits of the gas particle ID store how man stars this gas
	         particle gas already generated */

	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType)*8 - bits));
	      
	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);


	      SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

	      TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
#ifdef METALS
	      w = get_random_number(P[i].ID);
	      P[i].Metallicity += w * METAL_YIELD * (1 - exp(-p));
#endif

	      prob = P[i].Mass / mass_of_star * (1 - exp(-p));
#else /* belongs to ifndef(QUICK_LYALPHA) */

	      prob = 2.0;	/* this will always cause a star creation event */
	      if(bits == 0)
		number_of_stars_generated = 0;
	      else

		number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType)*8 - bits));

	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);

	      SphP[i].Sfr = 0;

#endif /* ends to QUICK_LYALPHA */

#ifdef DO_NOT_CREATE_STAR_PARTICLES
	      prob = -1;
#endif

	      if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		{
		  if(number_of_stars_generated == (GENERATIONS - 1))
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;

		      sum_mass_stars += P[i].Mass;

		      P[i].Type = 4;
		      TimeBinCountSph[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

#ifdef STELLARAGE
		      P[i].StellarAge = All.Time;
#endif
		    }
		  else
		    {
		      /* here we spawn a new star particle */

		      if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
			  fflush(stdout);
			  endrun(8888);
			}

		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;
#ifdef SNIA_HEATING
		      PPP[NumPart + stars_spawned].Hsml = All.SofteningTable[0];
#endif

		      NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		      FirstActiveParticle = NumPart + stars_spawned;
		      NumForceUpdate++;

		      TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;

		      PrevInTimeBin[NumPart + stars_spawned] = i;
		      NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		      if(NextInTimeBin[i] >= 0)
			PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		      NextInTimeBin[i] = NumPart + stars_spawned;
		      if(LastInTimeBin[P[i].TimeBin] == i)
			LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;

		      P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType)*8 - bits));

		      P[NumPart + stars_spawned].Mass = mass_of_star;
		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
		      sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifdef STELLARAGE
		      P[NumPart + stars_spawned].StellarAge = All.Time;
#endif
		      force_add_star_to_tree(i, NumPart + stars_spawned);

		      stars_spawned++;
		    }
		}

#ifdef METALS
	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		P[i].Metallicity += (1 - w) * METAL_YIELD * (1 - exp(-p));
#endif



#ifdef WINDS
	      /* Here comes the wind model */

	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		{
#ifdef VARIABLE_WINDS
      if (SphP[i].HostHaloMass > 0 && sm > 0)
        {
          r200c = pow(SphP[i].HostHaloMass / (4 * M_PI / 3.0 * 200 *rhocrit), 1.0/3.0);  /* physical r_200,crit value, assuming FoF mass = M_200,crit */
          v_esc = sqrt(All.G * SphP[i].HostHaloMass / r200c);  /* physical circular velocity at r_200,crit */
          c_halo = All.HaloConcentrationNorm * pow(SphP[i].HostHaloMass, All.HaloConcentrationSlope);
          v_esc *= sqrt( 2*c_halo / (log(1+c_halo) - c_halo/(1+c_halo))); /* physical escape velocity of halo */   
          v = All.VariableWindVelFactor * v_esc;  /* physical wind velocity */

          wind_momentum = sm * All.VariableWindSpecMomentum;
          wind_energy = sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);

          wind_mass = (wind_energy + sqrt(wind_energy*wind_energy + v*v*wind_momentum*wind_momentum)) / (v*v); 
	  /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
          p = wind_mass / P[i].Mass;

          /*printf("Mfof %f r200_c %f c_halo %f v_esc %f v %f sm %e wind_mass %e ml %f \n", SphP[i].HostHaloMass, r200c, c_halo, v_esc, v, sm, wind_mass, wind_mass/sm);*/
        }
      else
        {
          v = 0;
          p = 0;
        }
#else
		  p = All.WindEfficiency * sm / P[i].Mass;
#endif

		  prob = 1 - exp(-p);

		  if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		    {
#ifndef VARIABLE_WINDS
		      v =
			sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			     All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
#endif
#ifdef ISOTROPICWINDS
		      theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
		      phi = 2 * M_PI * get_random_number(P[i].ID + 4);

		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
#else
		      dir[0] = P[i].g.GravAccel[1] * P[i].Vel[2] - P[i].g.GravAccel[2] * P[i].Vel[1];
		      dir[1] = P[i].g.GravAccel[2] * P[i].Vel[0] - P[i].g.GravAccel[0] * P[i].Vel[2];
		      dir[2] = P[i].g.GravAccel[0] * P[i].Vel[1] - P[i].g.GravAccel[1] * P[i].Vel[0];
#endif

		      for(j = 0, norm = 0; j < 3; j++)
			norm += dir[j] * dir[j];

		      norm = sqrt(norm);
		      if(get_random_number(P[i].ID + 5) < 0.5)
			norm = -norm;

		      if(norm != 0)
			{
			  for(j = 0; j < 3; j++)
			    dir[j] /= norm;

			  for(j = 0; j < 3; j++)
			    {
			      P[i].Vel[j] += v * ascale * dir[j];
			      SphP[i].VelPred[j] += v * ascale * dir[j];
			    }

        SphP[i].DelayTime = All.WindFreeTravelMaxTimeFactor / hubble_function(All.Time);
			}
#if defined(MAGNETIC)
			x=0.;
#endif
		    }
		}
#endif  /* WINDS */
	    }
#if defined(MAGNETIC)
          SphP[i].XColdCloud = x;
#endif
	} /* End of If Type = 0*/
    }				/* end of main loop over active particles */


  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("SFR: spawned %d stars, converted %d gas particles into stars\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}


      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
    }
}

#if !defined(CS_MODEL)
double get_starformation_rate(int i)
{
  double rateOfSF;
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;



  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  flag = 1;			/* default is normal cooling */

  if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].d.Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].Ne;

  tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);
  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  cloudmass = x * P[i].Mass;

  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}
#endif

#endif /* closes MHM conditional */



#endif /* SFR */

#if defined(SFR) || defined(BLACK_HOLES) || defined(WINDTUNNEL) ||  defined(CA_BH_ACCRETION)
void rearrange_particle_sequence(void)
{
  int i, j, flag = 0, flag_sum;
  struct particle_data psave;

#if defined(BLACK_HOLES) || defined(WINDTUNNEL)
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

#ifdef SFR
  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
      flag = 1;
    }
#endif

#if defined(BLACK_HOLES) || defined(WINDTUNNEL) || defined(CA_BH_ACCRETION)
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	TimeBinCount[P[i].TimeBin]--;

	if(TimeBinActive[P[i].TimeBin])
	  NumForceUpdate--;

	if(P[i].Type == 0)
	  {
	    TimeBinCountSph[P[i].TimeBin]--;

	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];

	    N_gas--;

	    count_gaselim++;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];
	  }

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(count_elim)
    flag = 1;

  if(ThisTask == 0)
    {
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
#if defined(BLACK_HOLES)
  All.TotBHs -= tot_elim - tot_gaselim;
#endif
#endif

  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(flag_sum)
    reconstruct_timebins();
}
#endif /* closing of SFR-conditional */



#if defined(SFR)
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTime(egyhot, dens, &ne);

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
	x / pow(1 - x,
		2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) /
	(All.MaxSfrTimescale * coolrate);

      if(ThisTask == 0)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;

      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;

	  neff = -log(peff) * fac;

	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

#ifdef MODIFIEDBONDI
      All.BlackHoleRefDensity = thresholdStarburst;
      All.BlackHoleRefSoundspeed = sqrt(GAMMA * GAMMA_MINUS1 * egyeff);
#endif


      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      integrate_sfr();

      if(ThisTask == 0)
	{
	  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

	  printf("Isotherm sheet central density: %g   z0=%g\n",
		 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
		 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
	  fflush(stdout);

	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

#ifdef WINDS
      if(All.WindEfficiency > 0)
	if(ThisTask == 0)
	  printf("Windspeed: %g\n",
		 sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) /
		      All.WindEfficiency));
#endif
    }
}

void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2;
  FILE *fd;


  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = All.PhysDensThresh; rho <= 1000 * All.PhysDensThresh; rho *= 1.1)
    {
      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;
      tcool = GetCoolingTime(egyhot, rho, &ne);

      y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  if(ThisTask == 0)
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = All.PhysDensThresh; rho0 <= 10000 * All.PhysDensThresh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > All.PhysDensThresh)
	    {
	      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

	      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;
	      tcool = GetCoolingTime(egyhot, rho, &ne);

	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(All.PhysDensThresh / rho2) * All.MaxSfrTimescale;
	      factorEVP2 = pow(rho2 / All.PhysDensThresh, -0.8) * All.FactorEVP;
	      egyhot2 = All.EgySpecSN / (1 + factorEVP2) + All.EgySpecCold;
	      tcool2 = GetCoolingTime(egyhot2, rho2, &ne);
	      y2 =
		tsfr2 / tcool2 * egyhot2 / (All.FactorSN * All.EgySpecSN -
					    (1 - All.FactorSN) * All.EgySpecCold);
	      x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
	      egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
	      P2 = GAMMA_MINUS1 * rho2 * egyeff2;

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;


	      sigma_u4 += rho * dz;
	    }



	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
	}
    }


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(ThisTask == 0)
    fclose(fd);
}

#endif /* closing of SFR-conditional */


#if defined(SFR)
void set_units_sfr(void)
{
  double meanweight;

#ifdef COSMIC_RAYS
  double feedbackenergyinergs;
#endif

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


#ifdef CS_MODEL
#if defined(CS_SNII)
  All.TlifeSNII /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  All.MinTlifeSNI /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  All.MaxTlifeSNI /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
/* This corresponds to Raiteri estimations of SNII lifetimes for
   different metallicity ranges - note: this is in yr and then to internal*/

/* this is in years!! */
  Raiteri_COEFF_1 = 9.56077e6;	/* = old 0.00682912 in internal (9.8e8yr/h) */
  Raiteri_COEFF_2 = 1.73978e7;	/* 0.0124270 */
  Raiteri_COEFF_3 = 2.29347e7;	/* 0.0163819 */
  Raiteri_COEFF_4 = 2.14039e7;	/* 0.0152885 */
  Raiteri_COEFF_5 = 1.99753e7;	/* 0.0142681 */

  Raiteri_COEFF_1 /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  Raiteri_COEFF_2 /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  Raiteri_COEFF_3 /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  Raiteri_COEFF_4 /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
  Raiteri_COEFF_5 /= (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam);
#endif
#ifdef CS_FEEDBACK
  SN_Energy = All.SN_Energy_cgs / (All.UnitEnergy_in_cgs);	/*conversion to internal 
								   energy */
  if(ThisTask == 0)
    printf("Feedback energy per SN= %g ergs ,   %g internal units\n", All.SN_Energy_cgs, SN_Energy);
#endif
#endif

#ifdef COSMIC_RAYS
  if(All.CR_SNEff < 0.0)
    /* if CR_SNeff < 0.0, then substract CR Feedback energy from thermal
     * feedback energy
     */
    {
      if(ThisTask == 0)
	{
	  printf("%g percent of thermal feedback go into Cosmic Rays.\nRemaining ", -100.0 * All.CR_SNEff);
	}

      All.EgySpecSN *= (1.0 + All.CR_SNEff);
      All.CR_SNEff = -All.CR_SNEff / (1.0 + All.CR_SNEff);

    }

  All.FeedbackEnergy = All.FactorSN / (1 - All.FactorSN) * All.EgySpecSN;

  feedbackenergyinergs = All.FeedbackEnergy / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);

  if(ThisTask == 0)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n", feedbackenergyinergs);
      printf("OverDensThresh= %g\nPhysDensThresh= %g (internal units)\n", All.OverDensThresh,
	     All.PhysDensThresh);
    }
#endif
}


#endif /* closes SFR */

#endif /* closes COOLING */

#endif /* closes LT_STELLAREVOLUTION */
