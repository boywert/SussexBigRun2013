/*
 * recipe_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


void update_yields(int p, double dt, int nstep)
{
	int Zi;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp, NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual, NormSNIIMetalEjecRate_actual, NormSNIaMetalEjecRate_actual, NormAGBMetalEjecRate_actual;
	double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
	double MassDiff;
	double timet, sfh_time;
	int n; //Iterator used for loop over NOUT when updating MassWeightedAge
	double DiskSFR, step_width_times_DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units, inverse_DiskMass_physical_units;
	double Disk_total_metallicity, Bulge_total_metallicity;
	double BulgeSFR, step_width_times_BulgeSFR, BulgeSFR_physical_units, step_width_times_BulgeSFR_physical_units, inverse_BulgeMass_physical_units;
	double NormMassEjecRateSumAllTypes;


	double AgeCorrectionDisk[NOUT];
	double AgeCorrectionBulge[NOUT];

	for(n=0;n<NOUT;n++)
	{
		AgeCorrectionDisk[n] = 0.0;
		AgeCorrectionBulge[n] = 0.0;
	}

	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*Gal[p].SnapNum)+nstep;//Bin in Yield tables corresponding to current timestep
	timet = NumToTime(Gal[p].SnapNum) - (nstep + 0.5) * dt; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
  //NB: NumToTime(Gal[p].SnapNum) is the time to z=0 from start of current snapshot
	//    nstep is the number of the current timestep (0-19)
	//    dt is the width of one timestep within current snapshot


	//printf("AgeCorrectionDisk: %f\n", AgeCorrectionDisk);

    int i;
    for (i=0;i<=Gal[p].sfh_ibin;i++)
    {
    	sfh_time=Gal[p].sfh_t[i]+(0.5*Gal[p].sfh_dt[i]);
    //*****************************************
    //ENRICHMENT FROM DISK STARS INTO COLD GAS:
    //*****************************************
    if (Gal[p].sfh_DiskMass[i] > 0.0)
    {
     	//pre-calculations to sped up the code
    	DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_DiskSFR = timestep_width * DiskSFR;
    	DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h);
    	step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units;
    	inverse_DiskMass_physical_units=Hubble_h/(Gal[p].sfh_DiskMass[i]*1.0e10);
    	Disk_total_metallicity=metals_total(Gal[p].sfh_MetalsDiskMass[i])/Gal[p].sfh_DiskMass[i];

    	Zi = find_initial_metallicity(p, i, 1, 1);
    	//interpolate the disk luminosity on the lifetimeMetallicities tables
    	Zi_disp = (Disk_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	/*printf("NormSNIIMassEjecRate_actual = %.11f\n", NormSNIIMassEjecRate_actual);
    	printf("NormAGBMassEjecRate_actual = %.11f\n", NormAGBMassEjecRate_actual);
    	printf("NormSNIIMetalEjecRate_actual = %.11f\n", NormSNIIMetalEjecRate_actual);
    	printf("NormAGBMetalEjecRate_actual = %.11f\n\n", NormAGBMetalEjecRate_actual);*/

    	//pre-calculations to sped up the code
     	NormMassEjecRateSumAllTypes = NormSNIIMassEjecRate_actual + NormSNIaMassEjecRate_actual + NormAGBMassEjecRate_actual;

    	int k;
	    for (k=0;k<NUM_ELEMENTS;k++)
	    {
	    	NormSNIIYieldRate_actual[k] = NormSNIIYieldRate[TimeBin][i][Zi][k] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][k] - NormSNIIYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormSNIaYieldRate_actual[k] = NormSNIaYieldRate[TimeBin][i][Zi][k] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][k] - NormSNIaYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormAGBYieldRate_actual[k] = NormAGBYieldRate[TimeBin][i][Zi][k] + ((NormAGBYieldRate[TimeBin][i][Zi+1][k] - NormAGBYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    }

	    Gal[p].ColdGas += step_width_times_DiskSFR * NormMassEjecRateSumAllTypes;
#ifdef PORTINARI
	    Gal[p].MetalsColdGas.type2 += step_width_times_DiskSFR * (NormSNIIMetalEjecRate_actual + (Disk_total_metallicity * NormSNIIMassEjecRate_actual));
#endif
#ifdef CHIEFFI
	    Gal[p].MetalsColdGas.type2 += step_width_times_DiskSFR * NormSNIIMetalEjecRate_actual;
#endif
	    //Gal[p].MetalsColdGas.type1a += step_width_times_DiskSFR * (NormSNIaMetalEjecRate_actual + (Disk_total_metallicity * NormSNIaMassEjecRate_actual)); //ROB: SHOULD THESE NOT BE W/O THE UNSYNTH COMPONENT, AS THEY USE EJECTED METALS, NOT YIELDS? (19-07-12)
    	Gal[p].MetalsColdGas.type1a += step_width_times_DiskSFR * NormSNIaMetalEjecRate_actual;
    	Gal[p].MetalsColdGas.agb += step_width_times_DiskSFR * (NormAGBMetalEjecRate_actual + (Disk_total_metallicity * NormAGBMassEjecRate_actual)); //ROB: This one should be ok though? (19-07-12)
#ifdef PORTINARI
    	Gal[p].ColdGas_elements.H += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)); //ROB: Why no NormSNIaMassEjecRate_actual here ?! (04-04-12) ANS: SNIa yield tables quote 'ejected element masses' , NOT 'yields'. Therefore, no unsynth component to be calculated separately. (08-05-12)
    	Gal[p].ColdGas_elements.He += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    	Gal[p].ColdGas_elements.Cb += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.N += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.O += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Ne += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Mg += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Si += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.S += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Ca += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Fe += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#else
    	Gal[p].ColdGas_elements.O += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Mg += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Fe += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[p].ColdGas_elements.H += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual)); //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[p].ColdGas_elements.He += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    	Gal[p].ColdGas_elements.Cb += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.N += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.O += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Ne += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Mg += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Si += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.S += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Ca += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Fe += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
#else
    	Gal[p].ColdGas_elements.O += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Mg += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].ColdGas_elements.Fe += step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*(NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#endif //CHIEFFI

#ifdef NORMALIZE
#ifndef MAINELEMENTS
    	//Normalise sum of elements to total mass (only if considering all 11 elements):
    	float HB4;
    	HB4 = Gal[p].ColdGas_elements.H;
    	printf("H before  = %e [Msun]\n", Gal[p].ColdGas_elements.H);
    	MassDiff = ((Gal[p].ColdGas*1.0e10)/Hubble_h) - elements_total(Gal[p].ColdGas_elements); //IN [MSun]
    	printf("MassDiff  = %e [Msun]\n", MassDiff);
    	Gal[p].ColdGas_elements.H += (Gal[p].ColdGas_elements.H/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	printf("H after   = %e [Msun]\n", Gal[p].ColdGas_elements.H);
    	printf("H added   = %e [Msun]\n", (Gal[p].ColdGas_elements.H/elements_total(Gal[p].ColdGas_elements)) * MassDiff);
    	printf("Add check = %e [Msun]\n", Gal[p].ColdGas_elements.H-HB4);
    	Gal[p].ColdGas_elements.He += (Gal[p].ColdGas_elements.He/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Cb += (Gal[p].ColdGas_elements.Cb/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.N += (Gal[p].ColdGas_elements.N/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.O += (Gal[p].ColdGas_elements.O/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Ne += (Gal[p].ColdGas_elements.Ne/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Mg += (Gal[p].ColdGas_elements.Mg/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Si += (Gal[p].ColdGas_elements.Si/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.S += (Gal[p].ColdGas_elements.S/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Ca += (Gal[p].ColdGas_elements.Ca/elements_total(Gal[p].ColdGas_elements)) * MassDiff;
    	Gal[p].ColdGas_elements.Fe += (Gal[p].ColdGas_elements.Fe/elements_total(Gal[p].ColdGas_elements)) * MassDiff;

    	//Normlise total metals to total mass (only if considering all 11 elements):
    	MassDiff = ((metal_elements_total(Gal[p].ColdGas_elements)/1.0e10)*Hubble_h - metals_total(Gal[p].MetalsColdGas)); //IN [MSun*1.0e10/Hubble_h]
    	printf("Sum_metal_elements = %e [Msun]\n", metal_elements_total(Gal[p].ColdGas_elements));
    	printf("Sum_II_Ia_agb      = %e [Msun]\n", (metals_total(Gal[p].MetalsColdGas)*1.0e10)/Hubble_h);
    	printf("Metals_MassDiff    = %e [Msun]\n\n", metal_elements_total(Gal[p].ColdGas_elements) - (metals_total(Gal[p].MetalsColdGas)*1.0e10)/Hubble_h);
    	if (Gal[p].MetalsColdGas.type2 > 0.00008){ Gal[p].MetalsColdGas.type2 += (Gal[p].MetalsColdGas.type2/metals_total(Gal[p].MetalsColdGas)) * MassDiff; } //Cut at 0.00008*(10^9/Hubble_h) Msun is due to resolution issues with using floats.
    	if (Gal[p].MetalsColdGas.type1a > 0.00008){ Gal[p].MetalsColdGas.type1a += (Gal[p].MetalsColdGas.type1a/metals_total(Gal[p].MetalsColdGas)) * MassDiff; }
    	if (Gal[p].MetalsColdGas.agb > 0.00008){ Gal[p].MetalsColdGas.agb += (Gal[p].MetalsColdGas.agb/metals_total(Gal[p].MetalsColdGas)) * MassDiff; }
#endif //NORMALIZE
#endif //MAINELEMENTS

    	/*//if ((elements_total(Gal[p].ColdGas_elements)/(Gal[p].ColdGas/Hubble_h*1.0e10) < 0.0) || (elements_total(Gal[p].HotGas_elements)/(Gal[p].HotGas/Hubble_h*1.0e10) < 0.0) || (elements_total(Gal[p].DiskMass_elements)/(Gal[p].DiskMass/Hubble_h*1.0e10) < 0.0) || (elements_total(Gal[p].BulgeMass_elements)/(Gal[p].BulgeMass/Hubble_h*1.0e10) < 0.0) || (elements_total(Gal[p].ICM_elements)/(Gal[p].ICM/Hubble_h*1.0e10) < 0.0) || (elements_total(Gal[p].EjectedMass_elements)/(Gal[p].EjectedMass/Hubble_h*1.0e10) < 0.0))
    	if (Gal[p].HotGas == 0.0)
    	{
    		printf("------\n");
    		printf("SnapNum = %i, Step = %i\n", Gal[p].SnapNum, nstep);
    		printf("Type = %i\n", Gal[p].Type);
    		printf("ColdGas_elements/ColdGas     = %f\n", elements_total(Gal[p].ColdGas_elements)/(Gal[p].ColdGas/Hubble_h*1.0e10));
    		printf("DiskMass_elements/DiskMass   = %f\n", elements_total(Gal[p].DiskMass_elements)/(Gal[p].DiskMass/Hubble_h*1.0e10));
    		printf("HotGas_elements/HotGas       = %f\n", elements_total(Gal[p].HotGas_elements)/(Gal[p].HotGas/Hubble_h*1.0e10));
    		printf("BulgeMass_elements/BulgeMass = %f\n", elements_total(Gal[p].BulgeMass_elements)/(Gal[p].BulgeMass/Hubble_h*1.0e10));
    		printf("EjecMass_elements/EjecMass   = %f\n", elements_total(Gal[p].EjectedMass_elements)/(Gal[p].EjectedMass/Hubble_h*1.0e10));
    		printf("ICM_elements/ICM             = %f\n", elements_total(Gal[p].ICM_elements)/(Gal[p].ICM/Hubble_h*1.0e10));
    		printf("------\n");
    		printf("HotGas_elements  = %e, HotGas = %e [Msun]\n", elements_total(Gal[p].HotGas_elements), (Gal[p].HotGas/Hubble_h*1.0e10));
    		elements_print("HotGas", Gal[p].HotGas_elements);
    		printf("------\n\n");
    	}*/
/*if (Gal[p].EjectedMass > 0.0 && (elements_total(Gal[p].EjectedMass_elements)/(Gal[p].EjectedMass/Hubble_h*1.0e10) > 1.0001 || elements_total(Gal[p].EjectedMass_elements)/(Gal[p].EjectedMass/Hubble_h*1.0e10) < 0.999)) // && elements_total(Gal[p].EjectedMass_elements) > 0.0)
{
printf("------\n");
printf("ColdGas_elements/ColdGas     = %.11f\n", elements_total(Gal[p].ColdGas_elements)/(Gal[p].ColdGas/Hubble_h*1.0e10));
printf("EjecMass_elements/EjecMass   = %.11f\n", elements_total(Gal[p].EjectedMass_elements)/(Gal[p].EjectedMass/Hubble_h*1.0e10));
elements_print("Ejec", Gal[p].EjectedMass_elements);
printf("------\n\n");
}*/

    	//UPDATE DISK MASS COMPONENTS:
    	Gal[p].DiskMass -= step_width_times_DiskSFR * NormMassEjecRateSumAllTypes;
    	Gal[p].MetalsDiskMass.type2 -= step_width_times_DiskSFR * (Disk_total_metallicity * NormSNIIMassEjecRate_actual);
    	Gal[p].MetalsDiskMass.type1a -= step_width_times_DiskSFR * (Disk_total_metallicity * NormSNIaMassEjecRate_actual);
    	Gal[p].MetalsDiskMass.agb -= step_width_times_DiskSFR * (Disk_total_metallicity * NormAGBMassEjecRate_actual);
    	if (Gal[p].MetalsDiskMass.type2 < 0.0) Gal[p].MetalsDiskMass.type2 = 0.0;
    	if (Gal[p].MetalsDiskMass.type1a < 0.0) Gal[p].MetalsDiskMass.type1a = 0.0;
    	if (Gal[p].MetalsDiskMass.agb < 0.0) Gal[p].MetalsDiskMass.agb = 0.0;

    	Gal[p].DiskMass_elements.H -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].DiskMass_elements.He -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].DiskMass_elements.Cb -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].DiskMass_elements.N -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].DiskMass_elements.O -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].DiskMass_elements.Ne -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].DiskMass_elements.Mg -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].DiskMass_elements.Si -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].DiskMass_elements.S -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].DiskMass_elements.Ca -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].DiskMass_elements.Fe -= step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes);

    	//Update ages:
    	//Gal[p].MassWeightAge[Gal[p].SnapNum] -= (Gal[p].sfh_t[i]-timet)*(step_width_times_DiskSFR * NormMassEjecRateSumAllTypes);
    	//AgeCorrectionDisk += (Gal[p].sfh_t[i]-timet)*(step_width_times_DiskSFR * NormMassEjecRateSumAllTypes);
    	for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionDisk[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(step_width_times_DiskSFR * NormMassEjecRateSumAllTypes);
        	if (AgeCorrectionDisk[n] < 0.0) AgeCorrectionDisk[n] = 0.0;
        }
    	//if (AgeCorrectionDisk < 0.0) AgeCorrectionDisk = 0.0; //To account for small discrepancies between time of LAST SFH bin and time of current timestep (should be exactly the same).
    	//printf("+ %.11f (%i, %i, %i)\n", (Gal[p].sfh_t[i]-timet)*(step_width_times_DiskSFR * NormMassEjecRateSumAllTypes), Gal[p].SnapNum, nstep, i);

    	//NORMALISE (10-04-12):
#ifdef NORMALIZE
#ifndef MAINELEMENTS
    	//Normalise sum of elements to total mass (only if considering all 11 elements):
    	MassDiff = ((Gal[p].DiskMass*1.0e10)/Hubble_h) - elements_total(Gal[p].DiskMass_elements); //IN [MSun]
    	Gal[p].DiskMass_elements.H += (Gal[p].DiskMass_elements.H/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.He += (Gal[p].DiskMass_elements.He/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Cb += (Gal[p].DiskMass_elements.Cb/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.N += (Gal[p].DiskMass_elements.N/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.O += (Gal[p].DiskMass_elements.O/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Ne += (Gal[p].DiskMass_elements.Ne/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Mg += (Gal[p].DiskMass_elements.Mg/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Si += (Gal[p].DiskMass_elements.Si/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.S += (Gal[p].DiskMass_elements.S/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Ca += (Gal[p].DiskMass_elements.Ca/elements_total(Gal[p].DiskMass_elements)) * MassDiff;
    	Gal[p].DiskMass_elements.Fe += (Gal[p].DiskMass_elements.Fe/elements_total(Gal[p].DiskMass_elements)) * MassDiff;

    	//Normlise total metals to total mass (only if considering all 11 elements):
    	MassDiff = ((metal_elements_total(Gal[p].DiskMass_elements)/1.0e10)*Hubble_h - metals_total(Gal[p].MetalsDiskMass)); //IN [MSun*1.0e10/Hubble_h]
    	if (Gal[p].MetalsDiskMass.type2 > 0.00008){ Gal[p].MetalsDiskMass.type2 += (Gal[p].MetalsDiskMass.type2/metals_total(Gal[p].MetalsDiskMass)) * MassDiff; } //Cut at 0.00008*(10^9/Hubble_h) Msun is due to resolution issues with using floats.
    	if (Gal[p].MetalsDiskMass.type1a > 0.00008){ Gal[p].MetalsDiskMass.type1a += (Gal[p].MetalsDiskMass.type1a/metals_total(Gal[p].MetalsDiskMass)) * MassDiff; }
    	if (Gal[p].MetalsDiskMass.agb > 0.00008){ Gal[p].MetalsDiskMass.agb += (Gal[p].MetalsDiskMass.agb/metals_total(Gal[p].MetalsDiskMass)) * MassDiff; }
#endif
#endif
    }



    //*****************************************
    //ENRICHMENT FROM BULGE STARS INTO HOT GAS:
    //*****************************************
    if (Gal[p].sfh_BulgeMass[i] > 0.0)
    {
    	//pre-calculations to sped up the code
    	BulgeSFR = Gal[p].sfh_BulgeMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_BulgeSFR = timestep_width * BulgeSFR;
    	BulgeSFR_physical_units = BulgeSFR * (1.0e10/Hubble_h);
    	step_width_times_BulgeSFR_physical_units = timestep_width * BulgeSFR_physical_units;
    	inverse_BulgeMass_physical_units=Hubble_h/(Gal[p].sfh_BulgeMass[i]*1.0e10);
    	Bulge_total_metallicity=metals_total(Gal[p].sfh_MetalsBulgeMass[i])/Gal[p].sfh_BulgeMass[i];

    	Zi = find_initial_metallicity(p, i, 1, 2);
    	//interpolate the bulge luminosity on the lifetimeMetallicities tables
    	Zi_disp = (Bulge_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);

    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	//pre-calculations to sped up the code
      NormMassEjecRateSumAllTypes = NormSNIIMassEjecRate_actual + NormSNIaMassEjecRate_actual + NormAGBMassEjecRate_actual;

    	int k;
	    for (k=0;k<NUM_ELEMENTS;k++)
	    {
	    	NormSNIIYieldRate_actual[k] = NormSNIIYieldRate[TimeBin][i][Zi][k] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][k] - NormSNIIYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormSNIaYieldRate_actual[k] = NormSNIaYieldRate[TimeBin][i][Zi][k] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][k] - NormSNIaYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormAGBYieldRate_actual[k] = NormAGBYieldRate[TimeBin][i][Zi][k] + ((NormAGBYieldRate[TimeBin][i][Zi+1][k] - NormAGBYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    }

    	//UPDATE HOT GAS COMPONENTS:
    	Gal[p].HotGas += step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes;
    	//Gal[p].MetalsHotGas.type2 += step_width_times_BulgeSFR * (NormSNIIMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIIMassEjecRate_actual));
    	//Gal[p].MetalsHotGas.type1a += step_width_times_BulgeSFR * (NormSNIaMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIaMassEjecRate_actual));
#ifdef PORTINARI
    	Gal[p].MetalsHotGas.type2 += step_width_times_BulgeSFR * (NormSNIIMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIIMassEjecRate_actual));
#endif
#ifdef CHIEFFI
    	Gal[p].MetalsHotGas.type2 += step_width_times_BulgeSFR * NormSNIIMetalEjecRate_actual;
#endif
    	Gal[p].MetalsHotGas.type1a += step_width_times_BulgeSFR * NormSNIaMetalEjecRate_actual;
    	Gal[p].MetalsHotGas.agb += step_width_times_BulgeSFR * (NormAGBMetalEjecRate_actual + (Bulge_total_metallicity * NormAGBMassEjecRate_actual));
#ifdef PORTINARI
    	Gal[p].HotGas_elements.H += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.He += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.N += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.O += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Ne += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Mg += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Si += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.S += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Ca += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Fe += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#else
    	Gal[p].HotGas_elements.O += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Mg += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Fe += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[p].HotGas_elements.H += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[p].HotGas_elements.He += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.N += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.O += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Ne += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Mg += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Si += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.S += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Ca += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Fe += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
#else
    	Gal[p].HotGas_elements.O += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Mg += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
    	Gal[p].HotGas_elements.Fe += step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#endif //CHIEFFI

    	//NORMALISE (10-04-12):
#ifdef NORMALIZE
#ifndef MAINELEMENTS
    	//Normalise sum of elements to total mass (only if considering all 11 elements):
    	MassDiff = ((Gal[p].HotGas*1.0e10)/Hubble_h) - elements_total(Gal[p].HotGas_elements); //IN [MSun]
    	Gal[p].HotGas_elements.H += (Gal[p].HotGas_elements.H/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.He += (Gal[p].HotGas_elements.He/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Cb += (Gal[p].HotGas_elements.Cb/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.N += (Gal[p].HotGas_elements.N/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.O += (Gal[p].HotGas_elements.O/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Ne += (Gal[p].HotGas_elements.Ne/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Mg += (Gal[p].HotGas_elements.Mg/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Si += (Gal[p].HotGas_elements.Si/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.S += (Gal[p].HotGas_elements.S/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Ca += (Gal[p].HotGas_elements.Ca/elements_total(Gal[p].HotGas_elements)) * MassDiff;
    	Gal[p].HotGas_elements.Fe += (Gal[p].HotGas_elements.Fe/elements_total(Gal[p].HotGas_elements)) * MassDiff;

    	//Normlise total metals to total mass (only if considering all 11 elements):
    	MassDiff = ((metal_elements_total(Gal[p].HotGas_elements)/1.0e10)*Hubble_h - metals_total(Gal[p].MetalsHotGas)); //IN [MSun*1.0e10/Hubble_h]
    	if (Gal[p].MetalsHotGas.type2 > 0.0008){ Gal[p].MetalsHotGas.type2 += (Gal[p].MetalsHotGas.type2/metals_total(Gal[p].MetalsHotGas)) * MassDiff; } //Cut at 0.00008*(10^9/Hubble_h) Msun is due to resolution issues with using floats.
    	if (Gal[p].MetalsHotGas.type1a > 0.0008){ Gal[p].MetalsHotGas.type1a += (Gal[p].MetalsHotGas.type1a/metals_total(Gal[p].MetalsHotGas)) * MassDiff; }
    	if (Gal[p].MetalsHotGas.agb > 0.0008){ Gal[p].MetalsHotGas.agb += (Gal[p].MetalsHotGas.agb/metals_total(Gal[p].MetalsHotGas)) * MassDiff; }
#endif
#endif

    	//UPDATE BULGE MASS COMPONENTS:
    	Gal[p].BulgeMass -= step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes;
    	Gal[p].MetalsBulgeMass.type2 -= step_width_times_BulgeSFR * (Bulge_total_metallicity * NormSNIIMassEjecRate_actual);
    	Gal[p].MetalsBulgeMass.type1a -= step_width_times_BulgeSFR * (Bulge_total_metallicity * NormSNIaMassEjecRate_actual);
    	Gal[p].MetalsBulgeMass.agb -= step_width_times_BulgeSFR * (Bulge_total_metallicity * NormAGBMassEjecRate_actual);
    	if (Gal[p].MetalsBulgeMass.type2 < 0.0) Gal[p].MetalsBulgeMass.type2 = 0.0;
    	if (Gal[p].MetalsBulgeMass.type1a < 0.0) Gal[p].MetalsBulgeMass.type1a = 0.0;
    	if (Gal[p].MetalsBulgeMass.agb < 0.0) Gal[p].MetalsBulgeMass.agb = 0.0;

    	Gal[p].BulgeMass_elements.H -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].BulgeMass_elements.He -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Cb -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].BulgeMass_elements.N -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].BulgeMass_elements.O -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Ne -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].BulgeMass_elements.Mg -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Si -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].BulgeMass_elements.S -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
    	Gal[p].BulgeMass_elements.Ca -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);
#endif
    	Gal[p].BulgeMass_elements.Fe -= step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes);

    	//Update ages:
    	//Gal[p].MassWeightAge[Gal[p].SnapNum] -= (Gal[p].sfh_t[i]-timet)*(step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
    	//AgeCorrectionBulge += (Gal[p].sfh_t[i]-timet)*(step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
        for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionBulge[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
        	if (AgeCorrectionBulge[n] < 0.0) AgeCorrectionBulge[n] = 0.0;
        }
    	//if (AgeCorrectionBulge < 0.0) AgeCorrectionBulge = 0.0; //To account for small discrepancies between time of LAST SFH bin and time of current timestep (should be exactly the same).

    	//NORMALISE (10-04-12):
#ifdef NORMALIZE
#ifndef MAINELEMENTS
    	//Normalise sum of elements to total mass (only if considering all 11 elements):
    	MassDiff = ((Gal[p].BulgeMass*1.0e10)/Hubble_h) - elements_total(Gal[p].BulgeMass_elements); //IN [MSun]
    	Gal[p].BulgeMass_elements.H += (Gal[p].BulgeMass_elements.H/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.He += (Gal[p].BulgeMass_elements.He/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Cb += (Gal[p].BulgeMass_elements.Cb/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.N += (Gal[p].BulgeMass_elements.N/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.O += (Gal[p].BulgeMass_elements.O/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Ne += (Gal[p].BulgeMass_elements.Ne/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Mg += (Gal[p].BulgeMass_elements.Mg/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Si += (Gal[p].BulgeMass_elements.Si/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.S += (Gal[p].BulgeMass_elements.S/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Ca += (Gal[p].BulgeMass_elements.Ca/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;
    	Gal[p].BulgeMass_elements.Fe += (Gal[p].BulgeMass_elements.Fe/elements_total(Gal[p].BulgeMass_elements)) * MassDiff;

    	//Normlise total metals to total mass (only if considering all 11 elements):
    	MassDiff = ((metal_elements_total(Gal[p].BulgeMass_elements)/1.0e10)*Hubble_h - metals_total(Gal[p].MetalsBulgeMass)); //IN [MSun*1.0e10/Hubble_h]
    	if (Gal[p].MetalsBulgeMass.type2 > 0.00008){ Gal[p].MetalsBulgeMass.type2 += (Gal[p].MetalsBulgeMass.type2/metals_total(Gal[p].MetalsBulgeMass)) * MassDiff; } //Cut at 0.00008*(10^9/Hubble_h) Msun is due to resolution issues with using floats.
    	if (Gal[p].MetalsBulgeMass.type1a > 0.00008){ Gal[p].MetalsBulgeMass.type1a += (Gal[p].MetalsBulgeMass.type1a/metals_total(Gal[p].MetalsBulgeMass)) * MassDiff; }
    	if (Gal[p].MetalsBulgeMass.agb > 0.00008){ Gal[p].MetalsBulgeMass.agb += (Gal[p].MetalsBulgeMass.agb/metals_total(Gal[p].MetalsBulgeMass)) * MassDiff; }
#endif
#endif
    } //if (Gal[p].sfh_BulgeMass[i] > 0.0)
    } //for (i=0;i<=Gal[p].sfh_ibin;i++)

    //Update MassWeightAge from all SFH bins together:
    //printf("-------------\n");
   //printf("%.11f, %g\n", AgeCorrectionDisk, AgeCorrectionDisk);
    //printf("%.11f\n", Gal[p].MassWeightAge[Gal[p].SnapNum]);
    //Gal[p].MassWeightAge[Gal[p].SnapNum] -= (AgeCorrectionDisk+AgeCorrectionBulge);
    for(n=0;n<NOUT;n++)
    {
    	Gal[p].MassWeightAge[n] -= (AgeCorrectionDisk[n]+AgeCorrectionBulge[n]);
    }
    //printf("%.11f\n\n", Gal[p].MassWeightAge[Gal[p].SnapNum]);

}

int find_initial_metallicity(int p, int sfh_bin, int table_type, int component_type)
{
	if (component_type == 1) //Disk stars
	{
	int i, Zi_bin;
	double initMetals, Z_disk;

	initMetals = metals_total(Gal[p].sfh_MetalsDiskMass[sfh_bin]); //IN [10^10/h Msun]
	Zi_bin = -1;
	i = 0;
	if (initMetals == 0.0 || Gal[p].sfh_DiskMass[sfh_bin] == 0.0)
	{
		Z_disk = 0.0;
	}
	else Z_disk = initMetals/Gal[p].sfh_DiskMass[sfh_bin];

	switch (table_type)
	{
		case 1: //Lifetime metallicity table
			while (Zi_bin == -1)
			{
				if (lifetimeMetallicities[i] < Z_disk)
				{
					i++;
					if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
		case 2: //SN-II metallicity table
			while (Zi_bin == -1)
			{
				if (SNIIMetallicities[i] < Z_disk)
				{
					i++;
					if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
		//case 3 //SNIa yields are NOT metallicity dependent
		case 4: //AGB metallicity table
			while (Zi_bin == -1)
			{
				if (AGBMetallicities[i] < Z_disk)
				{
					i++;
					if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
	}

	if (Zi_bin == 0 ) return Zi_bin;
	else return Zi_bin-1;
	}
	else if (component_type == 2) //Bulge stars
	{
		int i, Zi_bin;
		double initMetals, Z_bulge;

		initMetals = metals_total(Gal[p].sfh_MetalsBulgeMass[sfh_bin]); //IN [10^10/h Msun]
		Zi_bin = -1;
		i = 0;
		if (initMetals == 0.0 || Gal[p].sfh_BulgeMass[sfh_bin] == 0.0)
		{
			Z_bulge = 0.0;
		}
		else Z_bulge = initMetals/Gal[p].sfh_BulgeMass[sfh_bin];

		switch (table_type)
		{
			case 1: //Lifetime metallicity table
				while (Zi_bin == -1)
				{
					if (lifetimeMetallicities[i] < Z_bulge) //Gal[p].sfh_MetalsDiskMass[sfh_bin].type2/Gal[p].sfh_DiskMass[sfh_bin])
					{
						i++;
						if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
			case 2: //SN-II metallicity table
				while (Zi_bin == -1)
				{
					if (SNIIMetallicities[i] < Z_bulge)
					{
						i++;
						if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
			//case 3 //SNIa yields are NOT metallicity dependent
			case 4: //AGB metallicity table
				while (Zi_bin == -1)
				{
					if (AGBMetallicities[i] < Z_bulge)
					{
						i++;
						if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
		}
		if (Zi_bin == 0 ) return Zi_bin;
		else return Zi_bin-1;
	}
	else { printf("Wrong stellar component type for Z_init calculation: Use either 1 (disk) or 2 (bulge)"); exit(1);}
}
