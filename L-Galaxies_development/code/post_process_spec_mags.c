#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* Created on: April 10, 2012
 *         by: Bruno Henriques (bmh20)
 */


/**@file post_process_spec_mags.c
 * @brief post_process_spec_mags.c can be used to compute mags or spectra from
 *        star formation histories. It also applies dust corrections.
 *
 *  When STAR_FORMATION_HISTORY option is ON in the Makefile the code
 *  stores the star formation history in the different components of
 *  the galaxy in log bins. If POST_PROCESS_MAGS is ON these star
 *  formation histories are used in this routine to compute magnitudes
 *  just before output.
 *
 * */
#ifndef NORMALIZEDDB
void post_process_spec_mags(struct GALAXY_OUTPUT *o)
#else
void post_process_spec_mags(struct GALAXY_OUTPUT *o, struct SFH_BIN *sfh_bins)
#endif
{
  int ll, nlum, ibin;
  double diskmass, bulgemass, icmmass, diskmetals, bulgemetals, icmmetals, previous_lum;
  double time, age;
  double Zg;
  //double Age_bin_yr=1.e11, Age_bin, bin_size;
  double Age_bin_yr=1.e7, Age_bin, bin_size;
  int i_AgeBin, MaxBins=1000, N_AgeBin;
// GL: 2014-03-13: ALos use an array of do next to simplify code below.
//     Otherwise would be full of ifdef-s
#ifndef NORMALIZEDDB
  struct SFH_BIN sfh_bins[SFH_NBIN];
  for(ll=0; ll<=o->sfh_ibin; ll++)
  {
	  sfh_bins[ll].sfh_DiskMass = o->sfh_DiskMass[ll];
	  sfh_bins[ll].sfh_BulgeMass = o->sfh_BulgeMass[ll];
	  sfh_bins[ll].sfh_ICM = o->sfh_ICM[ll];
	  sfh_bins[ll].sfh_MetalsDiskMass = o->sfh_MetalsDiskMass[ll];
	  sfh_bins[ll].sfh_MetalsBulgeMass = o->sfh_MetalsBulgeMass[ll];
	  sfh_bins[ll].sfh_MetalsICM = o->sfh_MetalsICM[ll];
  }
#endif


  //used for dust corrections
  Zg = metals_total(o->MetalsColdGas)/o->ColdGas/0.02;
    
  o->rbandWeightAge=0.0;






  for(nlum = 0; nlum < NMAG; nlum++)
    {
	  double LumDisk=0., LumBulge=0., LumICL=0.;
	  double ObsLumDisk=0., ObsLumBulge=0., ObsLumICL=0.;
	  double dObsLumDisk=0., dObsLumBulge=0., dObsLumICL=0.;
	  double dObsLumDisk_forward=0., dObsLumBulge_forward=0., dObsLumICL_forward=0.;

	  double YLumDisk=0., YLumBulge=0., YLumICL=0.;
	  double ObsYLumDisk=0., ObsYLumBulge=0., ObsYLumICL=0.;
	  double dObsYLumDisk=0., dObsYLumBulge=0., dObsYLumICL=0.;
	  double dObsYLumDisk_forward=0., dObsYLumBulge_forward=0., dObsYLumICL_forward=0.;

	  double LumDiskDust=0., LumBulgeDust=0., LumICLDust=0.;
	  double ObsLumDiskDust=0., ObsLumBulgeDust=0., ObsLumICLDust=0.;
	  double dObsLumDiskDust=0., dObsLumBulgeDust=0., dObsLumICLDust=0.;
	  double dObsLumDiskDust_forward=0., dObsLumBulgeDust_forward=0., dObsLumICLDust_forward=0.;

    previous_lum=0.;
    //loop on SFH bins
    ibin = 0;
	  for(ll=0; ll<=o->sfh_ibin; ll++)
	  {

	  	//compute the time at beginning and end of SFH bin, divide into 10Myr steps
	  	//All the SFH bins are divided into this smaller steps for the LUM calculation
	  //	bin_size=o->sfh_dt[ll]* Hubble_h/UnitTime_in_years;
	  //	Age_bin=Age_bin_yr* Hubble_h/UnitTime_in_years;

//#ifdef FULL_SPECTRA
	  /*	N_AgeBin=0;
	  	for (N_AgeBin=1;N_AgeBin<=MaxBins;N_AgeBin++)
	  	{
	  		if(bin_size/((float)N_AgeBin) < Age_bin) break;
	  		if(N_AgeBin==MaxBins)
	  			terminate("MaxBins reached in post_process_mags");
	  	}*/
//#else
	  	N_AgeBin=1;
//#endif

	  //	for(i_AgeBin=0;i_AgeBin<N_AgeBin;i_AgeBin++)
	 // 	{
	  		//sfh_time[ll] has been converted to years so needs to be in code units again
	  		//time at the older side of the bin
	  	//	time = (o->sfh_time[ll]+o->sfh_dt[ll]/2.)* Hubble_h/UnitTime_in_years;
	  		//time at the finer age bin
	  //		time -= (bin_size/((float)N_AgeBin)/2. + i_AgeBin*bin_size/((float)N_AgeBin));
	  		/*printf("ti=%e tf=%e tactual=%e\n",
	  				(o->sfh_time[ll]+o->sfh_dt[ll]/2.)* Hubble_h/UnitTime_in_years,
	  				(o->sfh_time[ll]-o->sfh_dt[ll]/2.)* Hubble_h/UnitTime_in_years,
	  				time);*/

	  	/* The stellar populations tables have magnitudes for all the mass
		   * formed in stars including what will be shortly lost by SNII   */
#ifdef DETAILED_METALS_AND_MASS_RETURN
	  		diskmass = sfh_bins[ll].sfh_DiskMass* 0.1 / Hubble_h / N_AgeBin;
	  		bulgemass = sfh_bins[ll].sfh_BulgeMass* 0.1 / Hubble_h / N_AgeBin;
	  		icmmass = sfh_bins[ll].sfh_ICM* 0.1 / Hubble_h / N_AgeBin;
#else
	  		diskmass = sfh_bins[ll].sfh_DiskMass* 0.1 / (Hubble_h * (1-RecycleFraction) ) / N_AgeBin;
	  		bulgemass = sfh_bins[ll].sfh_BulgeMass* 0.1 / (Hubble_h * (1-RecycleFraction) ) / N_AgeBin;
	  		icmmass = sfh_bins[ll].sfh_ICM* 0.1 / (Hubble_h * (1-RecycleFraction) ) / N_AgeBin;
#endif //DETAILED_METALS_AND_MASS_RETURN

	  		diskmetals = metals_total(sfh_bins[ll].sfh_MetalsDiskMass) / sfh_bins[ll].sfh_DiskMass;
	  		bulgemetals = metals_total(sfh_bins[ll].sfh_MetalsBulgeMass) / sfh_bins[ll].sfh_BulgeMass;
	  		icmmetals = metals_total(sfh_bins[ll].sfh_MetalsICM) / sfh_bins[ll].sfh_ICM;

//	  		 time = (o->sfh_time[ll])* Hubble_h/UnitTime_in_years;
// WRONG?	  		 time = SFH_t[o->SnapNum][ll][0]* Hubble_h/UnitTime_in_years;
	  		time = (SFH_t[o->SnapNum][0][ll]+SFH_dt[o->SnapNum][0][ll]/2.-NumToTime(o->SnapNum));
#ifdef DEBUG
	  		printf("snap=%d ll=%d: %e == %e ?\n",o->SnapNum,ll,time,SFH_t[o->SnapNum][0][ll]);
#endif
	  		//Disk MAGS
	  		if( sfh_bins[ll].sfh_DiskMass > 0.0)
	  			compute_post_process_lum (diskmass, time, diskmetals, o->SnapNum, nlum,
	  					&LumDisk, &ObsLumDisk, &dObsLumDisk, &dObsLumDisk_forward,
	  					&YLumDisk, &ObsYLumDisk, &dObsYLumDisk, &dObsYLumDisk_forward);

	  		//Bulge MAGS
	  		if(sfh_bins[ll].sfh_BulgeMass > 0.0)
	  			compute_post_process_lum (bulgemass, time, bulgemetals, o->SnapNum, nlum,
	  					&LumBulge, &ObsLumBulge, &dObsLumBulge, &dObsLumBulge_forward,
	  					&YLumBulge, &ObsYLumBulge, &dObsYLumBulge, &dObsYLumBulge_forward);

	  		//ICL MAGS
	  		if(sfh_bins[ll].sfh_ICM > 0.0)
	  			compute_post_process_lum (icmmass, time, icmmetals, o->SnapNum, nlum,
	  					&LumICL, &ObsLumICL, &dObsLumICL, &dObsLumICL_forward,
	  					&YLumICL, &ObsYLumICL, &dObsYLumICL, &dObsYLumICL_forward);

	  		//r-band weighted ages and mass weighted
	  		if((o->DiskMass+o->BulgeMass)>0.0 && nlum==17)
	  		{
	  			age = time -NumToTime(o->SnapNum);
	  			o->rbandWeightAge += (age * (LumDisk+LumBulge-previous_lum));
	  			previous_lum=LumDisk+LumBulge;
	  		}

	  //	}//end loop on smaller age bins
	  }//end of loop on sfh bins


#ifdef FULL_SPECTRA
#ifdef OUTPUT_REST_MAGS
	  o->Mag[nlum]=LumDisk+LumBulge;
	  o->MagBulge[nlum]=LumBulge;
#ifdef ICL
	  o->MagICL[nlum]=LumICL;
#endif
#endif
#ifdef COMPUTE_OBS_MAGS
	  o->ObsMag[nlum]=ObsLumDisk+ObsLumBulge;
	  o->ObsMagBulge[nlum]=ObsLumBulge;
#ifdef ICL
	  o->ObsMagICL[nlum]=ObsLumICL;
#endif

#ifdef OUTPUT_MOMAF_INPUTS
	  o->dObsMag[nlum]=dObsLumDisk+dObsLumBulge;
	  o->dObsMagBulge[nlum]=dObsLumBulge;
#ifdef ICL
	  o->dObsMagICL[nlum]=dObsLumICL;
#endif
#ifdef KITZBICHLER
	  o->dObsMag_forward[nlum]=dObsLumDisk_forward+dObsLumBulge_forward;
	  o->dObsMagBulge_forward[nlum]=dObsLumBulge_forward;
#ifdef ICL
	  o->dObsMagICL_forward[nlum]=dObsLumICL_forward;
#endif
#endif
#endif
#endif

#else //#ifndef FULL_SPECTRA

#ifdef OUTPUT_REST_MAGS
	  o->Mag[nlum]=lum_to_mag(LumDisk+LumBulge);
	  o->MagBulge[nlum]=lum_to_mag(LumBulge);
#ifdef ICL
	  o->MagICL[nlum]=lum_to_mag(LumICL);
#endif
#endif
#ifdef COMPUTE_OBS_MAGS
	  o->ObsMag[nlum]=lum_to_mag(ObsLumDisk+ObsLumBulge);
	  o->ObsMagBulge[nlum]=lum_to_mag(ObsLumBulge);
#ifdef ICL
	  o->ObsMagICL[nlum]=lum_to_mag(ObsLumICL);
#endif

#ifdef OUTPUT_MOMAF_INPUTS
	  o->dObsMag[nlum]=lum_to_mag(dObsLumDisk+dObsLumBulge);
	  o->dObsMagBulge[nlum]=lum_to_mag(dObsLumBulge);
#ifdef ICL
	  o->dObsMagICL[nlum]=lum_to_mag(dObsLumICL);
#endif
#ifdef KITZBICHLER
	  o->dObsMag_forward[nlum]=lum_to_mag(dObsLumDisk_forward+dObsLumBulge_forward);
	  o->dObsMagBulge_forward[nlum]=lum_to_mag(dObsLumBulge_forward);
#ifdef ICL
	  o->dObsMagICL_forward[nlum]=lum_to_mag(dObsLumICL_forward);
#endif
#endif
#endif
#endif
#endif //FULL_SPECTRA

#ifndef H2_AND_RINGS
	  o->CosInclination = fabs(o->StellarSpin[2]) /
	  		sqrt(o->StellarSpin[0]*o->StellarSpin[0]+
	  				 o->StellarSpin[1]*o->StellarSpin[1]+
	  		  	 o->StellarSpin[2]*o->StellarSpin[2]);
#else
	  o->CosInclination = fabs(o->DiskSpin[2]) /
	 	  		sqrt(o->DiskSpin[0]*o->DiskSpin[0]+
	 	  				 o->DiskSpin[1]*o->DiskSpin[1]+
	 	  		  	 o->DiskSpin[2]*o->DiskSpin[2]);
#endif

	  //Dust correction for disks (Inter-stellar Medium  + birth clouds)
	  if(o->ColdGas > 0.0)
	  dust_correction_for_post_processing(nlum, o->SnapNum, Zg, o->ColdGas, o->GasDiskRadius, o->CosInclination,
	   			                          LumDisk, ObsLumDisk, dObsLumDisk, dObsLumDisk_forward,
	   			                          YLumDisk, ObsYLumDisk, dObsYLumDisk, dObsYLumDisk_forward,
	   			                          &LumDiskDust, &ObsLumDiskDust, &dObsLumDiskDust, &dObsLumDiskDust_forward);

	  //Dust correction for bulges (remove light from young stars absorbed by birth clouds)
	  LumBulgeDust=LumBulge - YLumBulge * (1. - ExpTauBCBulge);
	  ObsLumBulgeDust=ObsLumBulge - ObsYLumBulge * (1. - ExpTauBCBulge);
	  dObsLumBulgeDust=dObsLumBulge - dObsYLumBulge * (1. - ExpTauBCBulge);
	  dObsLumBulgeDust_forward=dObsLumBulge_forward - dObsYLumBulge_forward * (1. - ExpTauBCBulge);

#ifdef FULL_SPECTRA
#ifdef OUTPUT_REST_MAGS
	  o->MagDust[nlum]=LumDiskDust+LumBulgeDust;
#endif
#ifdef COMPUTE_OBS_MAGS
	  o->ObsMagDust[nlum]=ObsLumDiskDust+ObsLumBulgeDust;
#ifdef OUTPUT_MOMAF_INPUTS
	  o->dObsMagDust[nlum]=dObsLumDiskDust+dObsLumBulgeDust;
#ifdef KITZBICHLER
	  o->dObsMagDust_forward[nlum]=dObsLumDiskDust_forward+dObsLumBulgeDust_forward;
#endif
#endif
#endif
#else //#ifndef FULL_SPECTRA
#ifdef OUTPUT_REST_MAGS
	  o->MagDust[nlum]=lum_to_mag(LumDiskDust+LumBulgeDust);
#endif
#ifdef COMPUTE_OBS_MAGS
	  o->ObsMagDust[nlum]=lum_to_mag(ObsLumDiskDust+ObsLumBulgeDust);
#ifdef OUTPUT_MOMAF_INPUTS
	  o->dObsMagDust[nlum]=lum_to_mag(dObsLumDiskDust+dObsLumBulgeDust);
#ifdef KITZBICHLER
	  o->dObsMagDust_forward[nlum]=lum_to_mag(dObsLumDiskDust_forward+dObsLumBulgeDust_forward);
#endif
#endif
#endif
#endif //#ifdef FULL_SPECTRA

	  if((o->DiskMass+o->BulgeMass)>0.0 && nlum==17)
	  	  {
          //LumDisk & LumBulge are sdss r-band luminosities (nlum==17)
	  			o->rbandWeightAge /= (LumDisk+LumBulge);
	  	  	o->rbandWeightAge = o->rbandWeightAge / 1000. * UnitTime_in_Megayears / Hubble_h; //conversion in age from code units/h -> Gyr
	  	   }

    }//end of loop on bands


}

void compute_post_process_lum (double mass, double time, double metals, int snap, int nlum,
						       double *Lum, double *ObsLum, double *dObsLum, double *dObsLum_forward,
						       double *YLum, double *ObsYLum, double *dObsYLum, double *dObsYLum_forward)
{
	double tbc;
	int metindex, tabindex, zindex;
	double f1, f2, fmet1, fmet2;
	double LumToAdd;

    /* Time bellow which the luminosities are corrected for extinction due to
      * molecular birth clouds.  */
    tbc = 10.0 / UnitTime_in_Megayears * Hubble_h;

	find_interpolated_lum(time,0., metals,&metindex, &tabindex, &f1, &f2, &fmet1, &fmet2);

	// reset met index to use only solar metallicity
	if(MetallicityOption == 0)
	  metindex = 4;

#ifdef OUTPUT_REST_MAGS

	zindex   = 0;
	LumToAdd = mass * (fmet1 * (f1 * LumTables[nlum][metindex][zindex][tabindex] +
			                    f2 * LumTables[nlum][metindex][zindex][tabindex + 1]) +
			           fmet2 * (f1 * LumTables[nlum][metindex + 1][zindex][tabindex] +
				                f2 * LumTables[nlum][metindex + 1][zindex][tabindex + 1]));

	*Lum += LumToAdd;
	if(time<=tbc)
	  *YLum += LumToAdd;


#endif

#ifdef COMPUTE_OBS_MAGS

	double ObsLumToAdd;

	zindex   = (LastDarkMatterSnapShot+1) - 1 - snap;
	ObsLumToAdd = mass * (fmet1 * (f1 * LumTables[nlum][metindex][zindex][tabindex] +
				                f2 * LumTables[nlum][metindex][zindex][tabindex + 1]) +
				       fmet2 * (f1 * LumTables[nlum][metindex + 1][zindex][tabindex] +
					            f2 * LumTables[nlum][metindex + 1][zindex][tabindex + 1]));

	*ObsLum += ObsLumToAdd;
	if(time<=tbc)
	  *ObsYLum += ObsLumToAdd;

#ifdef OUTPUT_MOMAF_INPUTS

	double dObsLumToAdd;

	zindex   = (LastDarkMatterSnapShot+1) - 1 - (snap-1);
	dObsLumToAdd = mass * (fmet1 * (f1 * LumTables[nlum][metindex][zindex][tabindex] +
					            f2 * LumTables[nlum][metindex][zindex][tabindex + 1]) +
					   fmet2 * (f1 * LumTables[nlum][metindex + 1][zindex][tabindex] +
						        f2 * LumTables[nlum][metindex + 1][zindex][tabindex + 1]));

	*dObsLum += dObsLumToAdd;
	if(time<=tbc)
	  *dObsYLum += dObsLumToAdd;
#ifdef KITZBICHLER
	zindex   = (LastDarkMatterSnapShot+1) - 1 - (snap+1);
	if(zindex < 0) zindex = 0; // TODO is this correct, or should dObsLumToAdd be set to 0?
		dObsLumToAdd = mass * (fmet1 * (f1 * LumTables[nlum][metindex][zindex][tabindex] +
					            f2 * LumTables[nlum][metindex][zindex][tabindex + 1]) +
					   fmet2 * (f1 * LumTables[nlum][metindex + 1][zindex][tabindex] +
						        f2 * LumTables[nlum][metindex + 1][zindex][tabindex + 1]));
	*dObsLum_forward += dObsLumToAdd;
	if(time<=tbc)
	  *dObsYLum_forward += dObsLumToAdd;

#endif
#endif
#endif

}



void dust_correction_for_post_processing(int nlum, int snap, double Zg, double ColdGas, double GasDiskRadius, double CosInclination,
		                                 double LumDisk, double ObsLumDisk, double dObsLumDisk, double dObsLumDisk_forward,
		                                 double YLumDisk, double ObsYLumDisk, double dObsYLumDisk, double dObsYLumDisk_forward,
		                                 double *LumDiskDust, double *ObsLumDiskDust, double *dObsLumDiskDust, double *dObsLumDiskDust_forward)
{
  double nh, tau, alam, sec, cosinc;
  double tauv, taubc, tauvbc, mu, VBand_WaveLength=0.55;
  float gasdev(long *idum);

  /* 0.94 = 2.83/3. - 3 to get scale lenght and 2.83 = 1.68^2 */
  nh = ColdGas / (M_PI * pow(GasDiskRadius * 0.94, 2) * 1.4);
  /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
  nh = nh / 3252.37;	// 3252.37 = 10^(3.5122)

  // redshift dependence of dust-to-ColdGas ratio
  nh = nh * pow(1 + ZZ[snap], -0.4);

  cosinc = CosInclination;
  if(cosinc < 0.2)
	cosinc = 0.2;		// minimum inclination ~80 degrees
  sec = 1.0 / cosinc;

  /* mu for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
   * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
  mu = -1.;
  while (mu < 0)
    {
	  mu = gasdev(&mu_seed) * MUWIDTH + MUCENTER;
	  if(mu < 0.1 || mu > 1.0)
		mu = -1.;
    }

  // extinction on Vband used as reference for the BC extinction
  tauv = get_extinction(NMAG, Zg, 0) * nh;

#ifdef OUTPUT_REST_MAGS

  //extinction due to ISM
  tau = get_extinction(nlum, Zg, 0) * nh;
  tau = tau * sec;
  if(tau > 0.0)
	alam = (1.0 - exp(-tau)) / tau;
  else
    alam = 1.;

  *LumDiskDust = LumDisk * alam;

  // now remove light from young stars absorbed by birth clouds
  tauvbc = tauv * (1. / mu - 1.); // TODO can be negative, correct?
  taubc = tauvbc * pow(FilterLambda[nlum] / VBand_WaveLength, -0.7);

  *LumDiskDust -=  (YLumDisk) * alam * (1. - exp(-taubc));

#endif


#ifdef OUTPUT_OBS_MAGS

  //extinction due to ISM
  tau = get_extinction(nlum, Zg, ZZ[snap]) * nh;
  tau = tau * sec;
  if(tau > 0.0)
	alam = (1.0 - exp(-tau)) / tau;
  else
	alam = 1.;

  *ObsLumDiskDust = ObsLumDisk * alam;


  // now remove light from young stars absorbed by birth clouds
  tauvbc = tauv * (1. / mu - 1.); // TODO same as before? no need to recalculate
  taubc = tauvbc * pow((FilterLambda[nlum] * (1. + ZZ[snap])) / VBand_WaveLength, -0.7);

  *ObsLumDiskDust -= (ObsYLumDisk) * alam * (1. - exp(-taubc));



#ifdef OUTPUT_MOMAF_INPUTS   // compute same thing at z + 1

  if(snap < (LastDarkMatterSnapShot+1) - 1)
	tau = get_extinction(nlum, Zg, ZZ[snap + 1]) * nh;
  else
	tau = get_extinction(nlum, Zg, ZZ[snap]) * nh;
  tau = tau * sec;
  if(tau > 0.0)
	alam = (1.0 - exp(-tau)) / tau;
  else
	alam = 1.;

  *dObsLumDiskDust = dObsLumDisk * alam;

  // now remove light from young stars absorbed by birth clouds
  if(snap < (LastDarkMatterSnapShot+1) - 1)
	taubc = tauvbc * pow((FilterLambda[nlum] * (1. + ZZ[snap + 1])) / VBand_WaveLength, -0.7);
  else
	taubc = tauvbc * pow((FilterLambda[nlum] * (1. + ZZ[snap])) / VBand_WaveLength, -0.7);

  *dObsLumDiskDust -= (dObsYLumDisk) * alam * (1. - exp(-taubc));

#ifdef KITZBICHLER
  if(snap > 0)
	tau = get_extinction(nlum, Zg, ZZ[snap - 1]) * nh;
  else
	tau = get_extinction(nlum, Zg, ZZ[snap]) * nh;
  tau = tau * sec;
  if(tau > 0.0)
	alam = (1.0 - exp(-tau)) / tau;
  else
	alam = 1.;

  *dObsLumDiskDust_forward = dObsLumDisk_forward * alam;

  // now remove light from young stars absorbed by birth clouds
  if(snap > 0)
	taubc = tauvbc * pow((FilterLambda[nlum] * (1. + ZZ[snap - 1])) / VBand_WaveLength, -0.7);
  else
	taubc = tauvbc * pow((FilterLambda[nlum] * (1. + ZZ[snap])) / VBand_WaveLength, -0.7);

  *dObsLumDiskDust_forward -= (dObsYLumDisk_forward) * alam * (1. - exp(-taubc));
#endif


#endif //OUTPUT_MOMAF_INPUTS

#endif //OUTPUT_OBS_MAGS

}













