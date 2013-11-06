#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"


#ifdef MCMC
/** @brief Writes galaxies into a structure to be used by the MCMC */

void save_galaxy_for_mcmc(int gal_index)
{
	//printf("ID=%lld snap=%d sampleID=%d\n", HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup,
	//		HaloGal[gal_index].SnapNum, SampleIDs[treenr]);
	int snap, fof, j, aa;
	float low_mass_limit, high_mass_limit, StellarMass;
	double log10_Hubble_h, PMass;

#ifdef MR_PLUS_MRII
	if(Switch_MR_MRII==1)
	{
		low_mass_limit=9.5;
		high_mass_limit=13.0;
	}
	else
		if(Switch_MR_MRII==2)
		{
			low_mass_limit=6.0;
			high_mass_limit=9.5;
		}
#else
	low_mass_limit=7.27;
	high_mass_limit=13.0;
#ifdef MRII
	low_mass_limit=6.0;
	high_mass_limit=11.27;
#endif
#endif


	log10_Hubble_h=log10(Hubble_h);

	//if(HaloGal[gal_index].SnapNum==53)
	//	printf("id=%lld\n",HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup);
	for(snap=0;snap<NOUT;snap++)
	{
		StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h)+ran3(&MCMCseed)*0.04*(1+MCMCConstraintsZZ[snap]);
		//StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h);

		for(fof=0;fof<NFofsInSample[snap]; fof++)
			if( StellarMass > low_mass_limit &&  StellarMass < high_mass_limit &&
					HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup == MCMC_FOF[fof].FoFID[snap])
			{

				//MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] = log10(1E10 * ((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h));
		  	//if((double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.)>0.9)
		  	//MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] += ran3(&MCMCseed)*0.04*(1+MCMCConstraintsZZ[snap]);
				MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] = StellarMass;
    		MCMC_GAL[TotMCMCGals[snap]].ColdGas[snap] = log10(1E10 * (HaloGal[gal_index].ColdGas*0.54/Hubble_h));
    		MCMC_GAL[TotMCMCGals[snap]].BulgeMass[snap] = log10(1E10 * HaloGal[gal_index].BulgeMass);
    		MCMC_GAL[TotMCMCGals[snap]].BlackHoleMass[snap] = log10(1E10 * HaloGal[gal_index].BlackHoleMass);
    		//in units of Solar Masses per yr
    		MCMC_GAL[TotMCMCGals[snap]].Sfr[snap]
    		= HaloGal[gal_index].Sfr * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;

#ifdef COMPUTE_SPECPHOT_PROPERTIES

    		struct GALAXY_OUTPUT galaxy_output;

    		//in case of postprocess magnitudes they are only calculates here, inside prepare
    		prepare_galaxy_for_output(snap, &HaloGal[gal_index], &galaxy_output);

    		MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = galaxy_output.MagDust[0]-5.*log10_Hubble_h;
    		//printf("MagU=%f ",MCMC_GAL[TotMCMCGals].MagU);
    		MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
    		//printf("MagB=%f ",MCMC_GAL[TotMCMCGals].MagB);
    		MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = galaxy_output.MagDust[2]-5.*log10_Hubble_h;
    		//printf("MagV=%f ",MCMC_GAL[TotMCMCGals].MagV);
    		MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = galaxy_output.MagDust[3]-5.*log10_Hubble_h;
    		//printf("MagJ=%f ",MCMC_GAL[TotMCMCGals].MagJ);
    		MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = galaxy_output.MagDust[4]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].MagK);
    		MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = galaxy_output.MagDust[5]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magu);
    		MCMC_GAL[TotMCMCGals[snap]].Magg[snap] = galaxy_output.MagDust[6]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magg);
    		MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = galaxy_output.MagDust[7]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magr);
    		MCMC_GAL[TotMCMCGals[snap]].Magi[snap] = galaxy_output.MagDust[8]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magi);
    		MCMC_GAL[TotMCMCGals[snap]].Magz[snap] = galaxy_output.MagDust[9]-5.*log10_Hubble_h;
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magz);
#endif

#ifdef HALOMODEL
    		if(snap==0)
    		{
    			MCMC_GAL[TotMCMCGals[snap]].fofid[snap] = fof;
    			MCMC_GAL[TotMCMCGals[snap]].fofmass[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
    			MCMC_GAL[TotMCMCGals[snap]].x[snap] = HaloGal[gal_index].Pos[0];
    			MCMC_GAL[TotMCMCGals[snap]].y[snap] = HaloGal[gal_index].Pos[1];
    			MCMC_GAL[TotMCMCGals[snap]].z[snap] = HaloGal[gal_index].Pos[2];
    			MCMC_GAL[TotMCMCGals[snap]].Type[snap] = HaloGal[gal_index].Type;
    			MCMC_GAL[TotMCMCGals[snap]].ngal[snap] = 0;
    		}
#endif

    		MCMC_GAL[TotMCMCGals[snap]].Weight[snap] = MCMC_FOF[fof].Weight[snap];


    		//#ifdef MR_PLUS_MRII
    		////

    		//NOW GET PROPERTIES FOR FOF GROUPS

    		++MCMC_FOF[fof].NGalsInFoF[snap];

#ifdef HALOMODEL
    		if(HaloGal[gal_index].Type==0)
    		{
    			MCMC_FOF[fof].IndexOfCentralGal[snap]=TotMCMCGals[snap];
    			MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
    		}
#endif

    		++TotMCMCGals[snap];

    		if(TotMCMCGals[snap] > MCMCAllocFactor)
    			terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
    	  }
      }

    return;
}
#endif


