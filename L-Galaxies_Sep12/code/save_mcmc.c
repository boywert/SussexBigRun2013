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
	int snap, tree;

	for(snap=0;snap<NOUT;snap++)
      {
		for(tree=0;tree<NTreesInSample[snap]; tree++)
		  if(HaloGal[gal_index].SnapNum == ListOutputSnaps[snap] &&
    	   log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)) > 7.0 &&
    	   HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup == SampleIDs[snap][tree])
    	  {
    		MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] = log10(1E10 * ((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h));
    		if(snap>0)
    			MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] += ran3(&MCMCseed)*0.25;
    		//printf("Mass=%f \n",MCMC_GAL[TotMCMCGals].StellarMass);

    		MCMC_GAL[TotMCMCGals[snap]].ColdGas[snap] = log10(1E10 * (HaloGal[gal_index].ColdGas/Hubble_h));
    		//printf("ColdMass=%f \n",MCMC_GAL[TotMCMCGals].ColdGas);

    		MCMC_GAL[TotMCMCGals[snap]].BulgeMass[snap] = log10(HaloGal[gal_index].BulgeMass * 1e10);

    		MCMC_GAL[TotMCMCGals[snap]].BlackHoleMass[snap] = log10(HaloGal[gal_index].BlackHoleMass * 1e10);

    		//in units of Solar Masses per Gyr
    		MCMC_GAL[TotMCMCGals[snap]].Sfr[snap]
    		= HaloGal[gal_index].Sfr * 1.e9 * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;

    		MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = lum_to_mag(HaloGal[gal_index].LumDust[0][snap])-5.*log10(Hubble_h);
    		//printf("MagB=%f ",MCMC_GAL[TotMCMCGals].MagB);

    		MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = lum_to_mag(HaloGal[gal_index].LumDust[1][snap])-5.*log10(Hubble_h);
    		//printf("MagV=%f ",MCMC_GAL[TotMCMCGals].MagV);

    		MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = lum_to_mag(HaloGal[gal_index].LumDust[2][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].MagK);

    		MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = lum_to_mag(HaloGal[gal_index].LumDust[3][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magu);

    		MCMC_GAL[TotMCMCGals[snap]].Magg[snap] = lum_to_mag(HaloGal[gal_index].LumDust[4][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magg);

    		MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = lum_to_mag(HaloGal[gal_index].LumDust[5][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magr);

    		MCMC_GAL[TotMCMCGals[snap]].Magi[snap] = lum_to_mag(HaloGal[gal_index].LumDust[6][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magi);

    		MCMC_GAL[TotMCMCGals[snap]].Magz[snap] = lum_to_mag(HaloGal[gal_index].LumDust[7][snap])-5.*log10(Hubble_h);
    		//printf("MagK=%f ",MCMC_GAL[TotMCMCGals].Magz);

    		MCMC_GAL[TotMCMCGals[snap]].Weight[snap] = Weight[snap][tree];
    		++TotMCMCGals[snap];

    		if(TotMCMCGals[snap] > MCMCAllocFactor)
    			terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
    	  }
      }

    return;
}
#endif


