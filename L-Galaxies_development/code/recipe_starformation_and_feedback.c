#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_starformation_and_feedback.c
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars
 *         formed from the cold gas, the amount of gas reheated from cold to hot
 *         and the amount of gas ejected from hot to external.
 *
 * The routine is divided in two parts, star formation and SN feedback, with a
 * number of different implementations controlled by input parameters.
 *
 *
 *  0 -\f$M_{\rm{crit}}=3.8\times 10^9
 *     \left(\frac{V_{\rm{max}}}{200\,\rm{km s}^{-1}}\right)
 *     \left(\frac{r_{\rm{disk}}}{10\,\rm{kpc}}\right)M_{\odot}\f$
 *     (Eq. 16 Guo2010) (StarFormationRecipe = 0), \n
 *        - same as 1 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *

 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackRecipe = 2)
 *     same as FeedbackRecipe = 1 * Vmax dependence.
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive
 * gas from their own satellites when these are outside Rvir of the type 0.
 * */


/** @brief Main recipe, calculates the fraction of cold gas turned into stars due
  *        to star formation; the fraction of mass instantaneously recycled and
  *        returned to the cold gas; the fraction of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   */
void starformation(int p, int centralgal, double time, double dt, int nstep)
{
	/*! Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt*/
  double tdyn, strdot=0., stars, cold_crit, metallicitySF;
#ifdef H2_AND_RINGS
  double strdotr[RNUM], starsr[RNUM], metallicityr[RNUM];
  double sfe, cold_crit_rate, tmp, sigmagas, sigmagasratio;
  //double sfe, vmax, cold_crit_rate, tmp, sigmagas, sigmagasratio, sigmah50=0.0;
  int j;
#endif
  

  if(Gal[p].Type == 0)
  {
  	tdyn = Gal[p].GasDiskRadius / Gal[p].Vmax;
  	cold_crit = SfrColdCrit * Gal[p].Vmax * Gal[p].GasDiskRadius;
  }
  else
  {
  	tdyn = Gal[p].GasDiskRadius / Gal[p].InfallVmax;
  	cold_crit = SfrColdCrit * Gal[p].InfallVmax * Gal[p].GasDiskRadius;
  }

  //standard star formation law (Croton2006, Delucia2007, Guo2010)
  if(StarFormationRecipe == 0)
  {
  	if(Gal[p].ColdGas > cold_crit)
  		strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);
  	else
  		strdot = 0.0;
  }

  //No threshold
  if(StarFormationRecipe == 1)
  {
  	if(Gal[p].ColdGas > 1.e-3)
  		strdot = SfrEfficiency * (Gal[p].ColdGas - 1.e-3) / tdyn * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);
  	else
  		strdot = 0.0;

  		/*	if(Gal[p].ColdGas > 1.e-5)
		  		strdot = SfrEfficiency * (Gal[p].ColdGas) / (2.e-5) * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);

	  			else
					strdot = 0.0;*/
  }

#ifdef H2_AND_RINGS
  update_h2fraction(p);

  sfe=SfrEfficiency/1.0e-4;    // the unit of sfe here is (km/s)/(Mpc/h)=1.022e-3 h^{+1} /Gyr
   if(SFRtdyn==1)
  	 sfe= sfe/1.8/tdyn; // for star formation rate proportional to 1/t_dyn

   for(j=0;j<RNUM;j++)
    {
       if(StarFormationRecipe == 2)
       	{
      	 tmp = sfe * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);
      	 if(Gal[p].Type == 0)
      		 cold_crit_rate = 0.19 * Gal[p].Vmax * Gal[p].GasDiskRadius/Gal[p].ColdGas;
      	  else
      	  	cold_crit_rate = 0.19 * Gal[p].InfallVmax * Gal[p].GasDiskRadius/Gal[p].ColdGas;
      	 if(cold_crit_rate < 1&&cold_crit_rate>=0)
       			strdotr[j] = tmp * Gal[p].ColdGasr[j] * (1 - cold_crit_rate);
       		else strdotr[j] = 0.0;
       	}

       else if(StarFormationRecipe == 3) /*The star formation law in Krumholz et al. 2009*/
       	{
       		if(j==0) sigmagas = Gal[p].ColdGasr[j] / (M_PI* radius[j]*radius[j])/Warmphasefactor*Clumpingfactor;
     	     else sigmagas = Gal[p].ColdGasr[j] / (M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]))/Warmphasefactor*Clumpingfactor;
         /* convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
         sigmagas=sigmagas*0.01*Hubble_h;
         sigmagasratio=pow(sigmagas/85.0,0.33);

  	      if(sigmagasratio<1.0&&sigmagasratio>0.0) strdotr[j] = sfe/sigmagasratio * Gal[p].ColdGasr[j]*Gal[p].H2fractionr[j]/Warmphasefactor ;
  	      	//Only cold H2 component is proportional to star formation rate.
  	      else if(sigmagasratio>=1.0) strdotr[j] = sfe*sigmagasratio * Gal[p].ColdGasr[j]*Gal[p].H2fractionr[j]/Warmphasefactor ;
  	      	else strdotr[j]=0.0;
       	}

       else if(StarFormationRecipe == 4)	/*The star formation law in Fu et al. 2010*/
       	{
       			 if(Gal[p].H2fractionr[j]>=0.0)
  	      	{
  	      		strdotr[j] = sfe * Gal[p].ColdGasr[j]*Gal[p].H2fractionr[j]/Warmphasefactor ; //Only cold H2 component is proportional to star formation rate.
  	      		//if(j==0) sigmah50=Gal[p].ColdGasr[0]/radius[0]/radius[0]/M_PI;
  	      		//  else sigmah50=Gal[p].ColdGasr[j]/(M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]));
  	      	}
  //	        else if(sigmah50>1.0e-7)
  //          {
  //           if(j==0) strdotr[j]=sfe*0.5*Gal[p].ColdGasr[j]/Warmphasefactor;
  //           	//Here 0.5 to keep the continuity of star formation rate at the radius where HI=H_2.
  //	      		else strdotr[j]=sfe*0.5*Gal[p].ColdGasr[j]/Warmphasefactor*(Gal[p].ColdGasr[j]/(M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]))/sigmah50);
  //	      		//The star formation law SFR\propt\Sigma_{gas}^2 in the mass form. m_{new star}\propto M_{gas}*(\Sigma_{gas}/Sigma_{50})
  //	      		//If f_H2<0.5 in ring 0, Sigma_{50}=Sigma_{gas}.
  //          }
  	      	else strdotr[j]=0.0;
       	}
       else  strdotr[j] = 0.0;
    }

    for (j=0,strdot=0;j<RNUM;j++) strdot+=strdotr[j];
#endif



#ifdef H2FORMATION
  /* use H2 formation model from Krumholz */
  Not yet implemented
  H2Mass=cal_H2(p);
  strdot = SfrEfficiency * H2Mass / tdyn;
#endif

  /*TODO - Note that Units of dynamical time are Mpc/Km/s - no conversion on dt needed
   *       be mentioned 3.06e19 to 3.15e19 */
  stars = strdot * dt;
  if(stars < 0.0)
    terminate("***error stars<0.0***\n");

#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    {
   	starsr[j] = strdotr[j] * dt;
   	if(starsr[j] <0.0) starsr[j] =0.0;
    }
#endif


#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if(stars > Gal[p].ColdGas)
  	stars = Gal[p].ColdGas;
#endif

  /*  update the star formation rate */
#ifdef SAVE_MEMORY
  /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
  Gal[p].Sfr += stars / (dt * STEPS);
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++) Gal[p].Sfrr[j] += starsr[j] / (dt * STEPS);
#endif
#else
  int outputbin;
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    if( Gal[p].SnapNum == ListOutputSnaps[outputbin]) {
      Gal[p].Sfr[outputbin] += stars / (dt * STEPS);
#ifdef H2_AND_RINGS
      for(j=0;j<RNUM;j++) Gal[p].Sfrr[j][outputbin] += starsr[j] / (dt * STEPS);
#endif
      break;
    }
  }
#endif
 



  mass_checks("recipe_starform #1",p);
  mass_checks("recipe_starform #1.1",centralgal);

  /* update for star formation
   * updates Mcold, StellarMass, MetalsMcold and MetalsStellarMass
   * in Guo2010 case updates the stellar spin -> hardwired, not an option */

  /* Store the value of the metallicity of the cold phase when SF occurs */
  if (Gal[p].ColdGas > 0.)
  	metallicitySF= metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
  else
    metallicitySF=0.;
 

	if (stars > 0.)
#ifndef H2_AND_RINGS
	update_stars_due_to_reheat(p, centralgal, &stars);
#else
	update_stars_due_to_reheat(p, centralgal, &stars, starsr);
#endif

  mass_checks("recipe_starform #2",p);
  mass_checks("recipe_starform #2.1",centralgal);


  // update_from_star_formation can only be called
   	// after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
   	// (star formation and feedback share the same fraction of cold gas)
   	if (stars > 0.)
 #ifndef H2_AND_RINGS
   		update_from_star_formation(p, stars, false, nstep); // false indicates not a burst
 #else
   	update_from_star_formation(p, stars, starsr, false, nstep); // false indicates not a burst
 #endif

   	update_massweightage(p, stars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
    /* ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN feedback is only called when stars die,
     * inside DETAILED_METALS_AND_MASS_RETURN */
  	if (stars > 0.)
#ifndef H2_AND_RINGS
  		SN_feedback(p, centralgal, stars, "ColdGas");
#else
  	SN_feedback(p, centralgal, stars, starsr, "ColdGas");
#endif

#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  Update the luminosities due to the stars formed */
  if (stars > 0.0)
    add_to_luminosities(p, stars, time, metallicitySF);
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  if(TrackDiskInstability)
    if(Gal[p].DiskMass > 0.0)
      check_disk_instability(p);
#ifndef H2_AND_RINGS
  if (DiskRadiusMethod == 2)
    get_stellar_disk_radius(p);
#endif

}





#ifndef H2_AND_RINGS
void update_stars_due_to_reheat(int p, int centralgal, double *stars)
#else
void update_stars_due_to_reheat(int p, int centralgal, double *stars, double starsr[])
#endif
{
	double MergeCentralVvir=0.;
#ifndef H2_AND_RINGS
	double fac;
	double CentralVvir=0.;
	double reheated_mass=0., ejected_mass=0.;
	/* SN FEEDBACK RECIPES */
#else
	 int j;
	 double reheated_massr[RNUM], facr[RNUM];
#endif

	  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
	   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
	   * for gas flow computations:
	   * If satellite is inside Rvir of main halo, Vvir of main halo used
	   * If it is outside, the Vvir of its central subhalo is used. */

	  //REHEAT
#ifndef H2_AND_RINGS
	 //if (strcmp(feedback_location,"HotGas")==0)
	//	 reheated_mass = 0.;
	// else
	//	 if (strcmp(feedback_location,"ColdGas")==0)
	//	 {
			 CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
			 MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

			 // Feedback depends on the circular velocity of the host halo
			 // Guo2010 - eq 18 & 19
			 if(FeedbackRecipe == 0)
			 {
				 if (Gal[Gal[p].CentralGal].Type == 0)
					 reheated_mass = FeedbackReheatingEpsilon * *stars *
					 (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
				 else
					 reheated_mass = FeedbackReheatingEpsilon * *stars *
					 (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));
			 }

			 //Make sure that the energy used in reheat does not exceed the SN energy (central subhalo Vvir used)
			 if (FeedbackRecipe == 0 || FeedbackRecipe == 1)
			 {
				 if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > *stars * (EtaSNcode * EnergySNcode))
					 reheated_mass = *stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
			 }

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
			 if((*stars + reheated_mass) > Gal[p].ColdGas)
			 {
				 fac = Gal[p].ColdGas / (*stars + reheated_mass);
				 *stars *= fac;
				 reheated_mass *= fac;
			 }
#endif
	//	 }// end if feedback_location

#else //H2_AND_RINGS
	  for(j=0;j<RNUM;j++)
	  {
	  	if(FeedbackRecipe == 0)
	  	{
	  		 if (Gal[Gal[p].CentralGal].Type == 0)
	  			 reheated_massr[j] = FeedbackReheatingEpsilon * starsr[j] *
	  			 (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
	  		 else
	  			 reheated_massr[j] = FeedbackReheatingEpsilon * starsr[j] *
	  			 (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

	  		/* Determine how much of the energy of SN feedback is used to reheat the
	  		 * gas compared to the amount used to eject gas */
	  		if (reheated_massr[j]*MergeCentralVvir*MergeCentralVvir > starsr[j]*(EtaSNcode*EnergySNcode))
	  			reheated_massr[j]=starsr[j]*(EtaSNcode*EnergySNcode)/(MergeCentralVvir*MergeCentralVvir);
	  	}

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
	  	/*  cant use more cold gas than is available! so balance SF and feedback */
	  	if((starsr[j] + reheated_massr[j]) > Gal[p].ColdGasr[j])
	  	{
	  		facr[j] = Gal[p].ColdGasr[j] / (starsr[j] + reheated_massr[j]);
	  		starsr[j] *= facr[j];
	  		reheated_massr[j] *= facr[j];
	  	}
#endif
	  }

	  *stars=0.0;

	  for(j=0;j<RNUM;j++)
	  	*stars+=starsr[j];

#endif //H2_AND_RINGS


}






/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
#ifndef H2_AND_RINGS
void update_from_star_formation(int p, double stars, bool flag_burst, int nstep)
#else
void update_from_star_formation(int p, double stars, double starsr[], bool flag_burst, int nstep)
#endif
{
  int i, j;
  double fraction;
#ifndef DETAILED_METALS_AND_MASS_RETURN
  double stars_nett=0.;
#endif
#ifdef H2_AND_RINGS
  double metallicityr[RNUM], massfrac, diskmass;
#endif
  if(Gal[p].ColdGas <= 0. || stars <= 0.) {
    printf("update_from_star_formation: Gal[p].ColdGas <= 0. || stars <= 0.\n");
    exit(0);
  }

#ifndef H2_AND_RINGS
#ifndef DETAILED_METALS_AND_MASS_RETURN
  stars_nett=(1 - RecycleFraction) * stars; //ROB: No longer an assumed instantaneous recycled fraction. Mass is returned over time via SNe and AGB winds.
#endif
#else
#ifndef DETAILED_METALS_AND_MASS_RETURN
   	for(j=0;j<RNUM;j++)
   		stars_nett+=(1 - RecycleFraction) * starsr[j];
#else
   	for(j=0;j<RNUM;j++)
   		stars+= starsr[j];
#endif
#endif

#ifndef H2_AND_RINGS
  /* Update the Stellar Spin when forming stars */ //ROB: This could be moved to after the yields, total metals and total ejected masses are updated.
#ifdef DETAILED_METALS_AND_MASS_RETURN
  if (Gal[p].DiskMass+stars > 1.e-8)
    for (i = 0; i < 3; i++)
      Gal[p].StellarSpin[i]=((Gal[p].StellarSpin[i])*(Gal[p].DiskMass) + stars*Gal[p].GasSpin[i])/(Gal[p].DiskMass+stars);
#else
  if (Gal[p].DiskMass+stars_nett > 1.e-8)
    for (i = 0; i < 3; i++)
      Gal[p].StellarSpin[i]=((Gal[p].StellarSpin[i])*(Gal[p].DiskMass)+stars_nett*Gal[p].GasSpin[i])/(Gal[p].DiskMass+stars_nett);
#endif //DETAILED_METALS_AND_MASS_RETURN
#else //H2_AND_RINGS
  diskmass=Gal[p].DiskMass+Gal[p].ColdGas;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  massfrac=stars/diskmass;
#else
  massfrac=stars_nett/diskmass;
#endif
  if(massfrac<1.0)
   for (j = 0; j <3 ; j++)
  	 Gal[p].DiskSpin[j]=Gal[p].DiskSpin[j]/(1-massfrac);
#endif
    /*  Update Gas and Metals from star formation */

  mass_checks("update_from_star_formation #0",p);

  //////transfer_gas_to_stars(p,"Disk",p,"Cold",stars_nett/Gal[p].ColdGas);
  //stars_nett=stars; //uncomment FOR DELAYED MASS RETURN
#ifdef DETAILED_METALS_AND_MASS_RETURN
  fraction=stars/Gal[p].ColdGas;
#else
  //fraction=stars_nett/Gal[p].ColdGas; //ROB: Old calculation for fraction of un-recycled ColdGas mass that was formed into stars
  fraction=stars_nett/Gal[p].ColdGas; //ROB: Old calculation for fraction of un-recycled ColdGas mass that was formed into stars
#endif

//STAR FORMATION HISTORY ARRAYS
#ifdef STAR_FORMATION_HISTORY
#ifdef DETAILED_METALS_AND_MASS_RETURN
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars; //ROB: Now, all SF gas is put in SFH array ("recycled' mass will return to gas phase over time)
  Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
  Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin] = elements_add(Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin],Gal[p].ColdGas_elements,fraction);
  //if (Gal[p].SnapNum == 16) {printf("SFH H_added = %f\n", Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin].H + (Gal[p].ColdGas_elements.H)*fraction);}
#endif
#ifdef TRACK_BURST
#ifdef TRACK_BURST_TEST
  Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars;
#else
  if (flag_burst) Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars;
#endif
#endif
#else //DETAILED_METALS_AND_MASS_RETURN
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars_nett; //ROB: Add amount of stars formed to SFH history of the disk. (NOTE: ALL SF OCCURS IN THE DISK. sfh_BulgeMass only increases when stars are transferred to the bulge before they explode)
  Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);
#ifdef TRACK_BURST
#ifdef TRACK_BURST_TEST
  Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars_nett;
#else
  if (flag_burst) Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars_nett;
#endif
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
#endif //STAR_FORMATION_HISTORY


//GLOBAL METALLICITY ARRAYS
#ifndef H2_AND_RINGS
  Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Gal[p].MetalsColdGas,fraction);
  Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsColdGas,-fraction);
#else
 for(j=0;j<RNUM;j++)
 {
	 metallicityr[j]=metals_total(Gal[p].MetalsColdGasr[j])/Gal[p].ColdGasr[j];
#ifdef DETAILED_METALS_AND_MASS_RETURN
	 Gal[p].MetalsColdGasr[j] -= metallicityr[j] * starsr[j];
	 Gal[p].MetalsColdGas -= metallicityr[j] * starsr[j];
	 Gal[p].MetalsDiskMassr[j] += metallicityr[j] * starsr[j];
	 Gal[p].MetalsDiskMass += metallicityr[j] * starsr[j];
#else
	 Gal[p].MetalsColdGasr[j] -= metallicityr[j] * (1 - RecycleFraction) * starsr[j];
	 Gal[p].MetalsColdGas -= metallicityr[j] * (1 - RecycleFraction) * starsr[j];
   Gal[p].MetalsDiskMassr[j] += metallicityr[j] * (1 - RecycleFraction) * starsr[j];
   Gal[p].MetalsDiskMass += metallicityr[j] * (1 - RecycleFraction) * starsr[j];
#endif
 }
#endif

 //GLOBAL PROPERTIES
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS
  Gal[p].DiskMass_elements=elements_add(Gal[p].DiskMass_elements,Gal[p].ColdGas_elements,fraction);
  Gal[p].ColdGas_elements=elements_add(Gal[p].ColdGas_elements,Gal[p].ColdGas_elements,-fraction);
#endif
  Gal[p].DiskMass += stars;
  Gal[p].ColdGas -= stars;
#ifdef H2_AND_RINGS
 for(j=0;j<RNUM;j++)
 {
	 Gal[p].DiskMassr[j] += starsr[j];
	 Gal[p].ColdGasr[j] -= starsr[j];
 }
#endif
#ifdef TRACK_BURST
#ifdef TRACK_BURST_TEST
  Gal[p].BurstMass+=stars;
#else
  if (flag_burst) Gal[p].BurstMass+=stars;
#endif
#endif
#else //DETAILED_METALS_AND_MASS_RETURN
  Gal[p].DiskMass += stars_nett;
  Gal[p].ColdGas -= stars_nett;
#ifdef H2_AND_RINGS
 for(j=0;j<RNUM;j++)
 {
	 Gal[p].DiskMassr[j] += (1 - RecycleFraction) * starsr[j];
	 Gal[p].ColdGasr[j] -= (1 - RecycleFraction) * starsr[j];
 }
#endif
#ifdef TRACK_BURST
#ifdef TRACK_BURST_TEST
  Gal[p].BurstMass+=stars_nett;
#else
  if (flag_burst) Gal[p].BurstMass+=stars_nett;
#endif
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN

  //SIMPLE VERION OF MASS RETURN
 /* if(nstep>-1)
  {
  	double previoustime,newtime, deltaT,time, bin_time, mass, bulge_mass=0., time_present;

  	previoustime = NumToTime(Halo[Gal[p].HaloNr].SnapNum-1);
  	newtime = NumToTime(Halo[Gal[p].HaloNr].SnapNum);
  	// Time between snapshots
  	deltaT = previoustime - newtime;
  	time = previoustime - (nstep + 0.5) * (deltaT / STEPS);

  	// Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars_nett; //ROB: Add amount of stars formed to SFH history of the disk. (NOTE: ALL SF OCCURS IN THE DISK. sfh_BulgeMass only increases when stars are transferred to the bulge before they explode)
  	// Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);

  	time_present=time_to_present(1000.);

  //REMOVE FOR DELAYED MASS RETURN
  	for(i=0;i<=Gal[p].sfh_ibin;i++)
  	{
  		//bin_time=(Gal[p].sfh_t[i]+Gal[p].sfh_dt[i]/2.)*UnitTime_in_years/Hubble_h;
  		bin_time=(Gal[p].sfh_t[i]+Gal[p].sfh_dt[i]/2.-time)*UnitTime_in_years/Hubble_h;

  		if(bin_time>(time_present-time)*UnitTime_in_years/Hubble_h)
  		{
  			printf("current time=%f\n",(time_present-time)*UnitTime_in_years/Hubble_h/1.e9);//time elapsed
  			//time to current present
  			printf("p=%d gal_snap=%d, snap=%d step=%d bin=%d time=%f sfh=%f flag=%d SFH_TIME=%f\n",p,Gal[p].SnapNum,
  					Halo[Gal[p].HaloNr].SnapNum-1, nstep, i, bin_time/1.e9, Gal[p].sfh_DiskMass[i],Gal[p].sfh_flag[i],
  					(SFH_t[Halo[Gal[p].HaloNr].SnapNum-1][nstep][i]+SFH_dt[Halo[Gal[p].HaloNr].SnapNum-1][nstep][i]/2.)*UnitTime_in_years/Hubble_h/1.e9);

  			printf("ola2\n");
  			exit(0);
  		}

  	 //if(i==0)
  		// printf("ola  ");
  	// if(bin_time/1.e9 > 1.0 && bin_time/1.e9<2.0 && Gal[p].sfh_flag[i]==0)
  		if(bin_time/1.e9 > 0.0 && bin_time/1.e9<2.0)
  		{
  		// printf("i=%d ibin=%d\n",i,Gal[p].sfh_ibin);
  			Gal[p].sfh_flag[i]=1;
  			mass=Gal[p].sfh_DiskMass[i]*0.43*(deltaT / STEPS)/2.*UnitTime_in_years/Hubble_h/1.e9;
  		 //printf("mass=%0.2e\n",mass*1.e10);
  			if(mass>Gal[p].DiskMass)
  			{
  				bulge_mass=mass-Gal[p].DiskMass;
  				mass=Gal[p].DiskMass;
  			}
  			Gal[p].sfh_DiskMass[i]-= mass;
  			Gal[p].DiskMass -= mass;
  			Gal[p].ColdGas += mass;
  			fraction=mass/Gal[p].DiskMass;
  			Gal[p].sfh_MetalsDiskMass[i] = metals_add(Gal[p].sfh_MetalsDiskMass[i],Gal[p].sfh_MetalsDiskMass[i],-fraction);
  			Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Gal[p].MetalsDiskMass,-fraction);
  			Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsDiskMass,+fraction);
  			if(bulge_mass>0.)
  			{
  				if(bulge_mass>Gal[p].BulgeMass)
  					bulge_mass=Gal[p].BulgeMass;
  				Gal[p].sfh_BulgeMass[i]-= bulge_mass;
  				Gal[p].BulgeMass -= bulge_mass;
  				Gal[p].ColdGas += bulge_mass;
  				fraction=bulge_mass/Gal[p].BulgeMass;
  				Gal[p].sfh_MetalsDiskMass[i] = metals_add(Gal[p].sfh_MetalsBulgeMass[i],Gal[p].sfh_MetalsBulgeMass[i],-fraction);
  				Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsBulgeMass,Gal[p].MetalsBulgeMass,-fraction);
  				Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsBulgeMass,+fraction);
  			}
  		}



  	}

  }*/


  mass_checks("update_from_star_formation #1",p);

  /* Formation of new metals - instantaneous recycling approximation - only SNII
   * Also recompute the metallicity of the cold phase.
   * TODO In bursts, gas is NOT later transferred from cold to hot - should it be? */
  /* TODO Note that FracZtoHot is broken as there may not be any hot gas.  If metals
   * are given to the hot phase then mass needs to be also. */

  /*DELAYED ENRICHMENT AND MASS RETURN, USING INDIVIDUAL ELEMENT YIELDS: No fixed yield or recycling fraction anymore:*/
  if (FeedbackRecipe == 0 || FeedbackRecipe == 1) {
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef METALS_SELF
	Gal[p].MetalsHotGasSelf.type2 += Yield * FracZtoHot * stars;
#endif
#else //DETAILED_METALS_AND_MASS_RETURN

#ifndef H2_AND_RINGS
    Gal[p].MetalsColdGas += Yield * (1.0 - FracZtoHot) * stars;
#else
    for(j=0;j<RNUM;j++)
           {
           	Gal[p].MetalsColdGasr[j] += Yield * (1.0 - FracZtoHot) * starsr[j];
           	Gal[p].MetalsColdGas += Yield * (1.0 - FracZtoHot) * starsr[j];
           }
#endif
    Gal[Gal[p].CentralGal].MetalsHotGas += Yield * FracZtoHot * stars;// FracZtoHot=0 not used
#ifdef METALS_SELF
    Gal[p].MetalsHotGasSelf += Yield * FracZtoHot * stars;
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
  }

#ifndef H2_AND_RINGS
  if (DiskRadiusMethod == 2)
    get_stellar_disk_radius(p);
#endif

}







/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - onle ejection.*/
#ifndef H2_AND_RINGS
void SN_feedback(int p, int centralgal, double stars, char feedback_location[])
#else
void SN_feedback(int p, int centralgal, double stars, double starsr[], char feedback_location[])
#endif
{
	 double MergeCentralVvir=0., EjectVmax, EjectVvir, SN_Energy, Reheat_Energy;
	 double reheated_mass=0., ejected_mass=0.;
	/* SN FEEDBACK RECIPES */
#ifndef H2_AND_RINGS
	 double CentralVvir;
#else
	 int j;
	 double reheated_massr[RNUM];
#endif

	  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
	   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
	   * for gas flow computations:
	   * If satellite is inside Rvir of main halo, Vvir of main halo used
	   * If it is outside, the Vvir of its central subhalo is used. */

	  //REHEAT
#ifndef H2_AND_RINGS
	 if (strcmp(feedback_location,"HotGas")==0)
		 reheated_mass = 0.;
	 else
		 if (strcmp(feedback_location,"ColdGas")==0)
		 {
			 CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
			 MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

			 mass_checks("recipe_starform #0",p);
			 mass_checks("recipe_starform #0.1",centralgal);

			 // Feedback depends on the circular velocity of the host halo
			 // Guo2010 - eq 18 & 19
			 if(FeedbackRecipe == 0)
			 {
				 if (Gal[Gal[p].CentralGal].Type == 0)
					 reheated_mass = FeedbackReheatingEpsilon * stars *
					 (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
				 else
					 reheated_mass = FeedbackReheatingEpsilon * stars *
					 (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));
			 }

			 //Make sure that the energy used in reheat does not exceed the SN energy (central subhalo Vvir used)
			 if (FeedbackRecipe == 0 || FeedbackRecipe == 1)
			 {
				 if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * EnergySNcode))
					 reheated_mass = stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
			 }


			 if(reheated_mass > Gal[p].ColdGas)
	  		reheated_mass = Gal[p].ColdGas;

		 }// end if feedback_location

#else //H2_AND_RINGS
	  for(j=0;j<RNUM;j++)
	  {
	  	if(FeedbackRecipe == 0)
	  	{
	  		 if (Gal[Gal[p].CentralGal].Type == 0)
	  			 reheated_massr[j] = FeedbackReheatingEpsilon * starsr[j] *
	  			 (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
	  		 else
	  			 reheated_massr[j] = FeedbackReheatingEpsilon * starsr[j] *
	  			 (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

	  		/* Determine how much of the energy of SN feedback is used to reheat the
	  		 * gas compared to the amount used to eject gas */
	  		if (reheated_massr[j]*MergeCentralVvir*MergeCentralVvir > starsr[j]*(EtaSNcode*EnergySNcode))
	  			reheated_massr[j]=starsr[j]*(EtaSNcode*EnergySNcode)/(MergeCentralVvir*MergeCentralVvir);
	  	}


	  	if(reheated_mass > Gal[p].ColdGas)
	  		if((starsr[j] + reheated_massr[j]) > Gal[p].ColdGasr[j])
	  			reheated_massr[j] = Gal[p].ColdGasr[j];

	  }

	  reheated_mass=0.0;
	  stars=0;

	  for(j=0;j<RNUM;j++)
	  {
	  	reheated_mass+=reheated_massr[j];
	  	stars+=starsr[j];
	  }


#endif //H2_AND_RINGS


	  /* Determine ejection (for FeedbackRecipe 2 we have the dependence on Vmax)
	   * Guo2010 - eq 22
	   * Note that satellites can now retain gas and have their own gas cycle*/
	  //EJECT
	  if (Gal[Gal[p].CentralGal].Type == 0)
	    {
		  EjectVmax=Gal[centralgal].Vmax;
		  EjectVvir=Gal[centralgal].Vvir;// main halo Vvir
	    }
	  else
	    {
		  EjectVmax=Gal[Gal[p].CentralGal].InfallVmax;
		  EjectVvir=Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
	    }

	  if(FeedbackRecipe == 0)
	  {
	      ejected_mass =
	    	(FeedbackEjectionEfficiency* (EtaSNcode * EnergySNcode) * stars *
		     min(1./FeedbackEjectionEfficiency, .5+1/pow(EjectVmax/EjectPreVelocity,EjectSlope)) -
		     reheated_mass*EjectVvir*EjectVvir) /(EjectVvir*EjectVvir);
	  }
	  else if(FeedbackRecipe == 1)//the ejected material is assumed to have V_SN
	    {

		  SN_Energy = .5 * stars * (EtaSNcode * EnergySNcode);
		  Reheat_Energy = .5 * reheated_mass * EjectVvir * EjectVvir;

		  ejected_mass = (SN_Energy - Reheat_Energy)/(0.5 * FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode));

		  //if VSN^2<Vvir^2 nothing is ejected
		  if(FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode)<EjectVvir*EjectVvir)
		  		ejected_mass =0.0;
	     }

	  // Finished calculating mass exchanges, so just check that none are negative
	  if (reheated_mass < 0.0) reheated_mass = 0.0;
	  if (ejected_mass < 0.0) ejected_mass = 0.0;

	  /* Update For Feedback */
	   /* update cold, hot, ejected gas fractions and respective metallicities
	    * there are a number of changes introduced by Guo2010 concerning where
	    * the gas ends up */

#ifdef H2_AND_RINGS
	  double lold=0.0, lnew=0.0;
	   if(Gal[p].ColdGas>1.0e-6)
	   	{for(j=0;j<RNUM;j++) lold+=(Gal[p].ColdGasr[j]*radius[j]*Gal[p].Vvir);}
#endif

	  //ejected_mass = 0.01*Gal[centralgal].HotGas;
	  if (reheated_mass + ejected_mass > 0.)
	  {
#ifndef H2_AND_RINGS
	  	update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
#else
	  	update_from_feedback(p, centralgal, reheated_mass, ejected_mass,  reheated_massr);
#endif
#ifdef H2_AND_RINGS
	  	if(Gal[p].ColdGas>1.0e-6)
	  	{for(j=0;j<RNUM;j++) lnew+=(Gal[p].ColdGasr[j]*radius[j]*Gal[p].Vvir); }

	  	       Gal[p].ReheatedGas+=reheated_mass;
	  	       Gal[p].ReheatedL+=(lold-lnew);
#endif
	  }

}


/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
#ifndef H2_AND_RINGS
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass)
#else
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double reheated_massr[])
#endif
{
  double dis=0.;
  double massremain;
  double fraction;
  int merger_centre;
#ifdef H2_AND_RINGS
  double metallicityr[RNUM];
  double coldgas=Gal[p].ColdGas;
  double metallicity=metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
  int j;
#endif

  //mass_checks("update_from_feedback #1",p);

  if(FeedbackRecipe == 0 || FeedbackRecipe == 1)
  {
  	if(Gal[p].ColdGas > 0.)
  	{
  		//REHEAT
  		// if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas being striped,
  		//some of the reheated and ejected masses goes to the type 0 and some stays in the type 1

#ifdef H2_AND_RINGS
  		for(j=0;j<RNUM;j++)
  		{
  			metallicityr[j]=metals_total(Gal[p].MetalsColdGasr[j])/Gal[p].ColdGasr[j];
  			Gal[p].ColdGasr[j] -= reheated_massr[j];
  			Gal[p].ColdGas -= reheated_massr[j];
  			Gal[p].MetalsColdGasr[j] -= metallicityr[j] * reheated_massr[j];
  			Gal[p].MetalsColdGas -= metallicityr[j] * reheated_massr[j];
  		}
#endif


  		if(Gal[p].Type ==0)
  		{
#ifndef H2_AND_RINGS
  			transfer_gas(p,"Hot",p,"Cold",((float)reheated_mass)/Gal[p].ColdGas,"update_from_feedback", __LINE__);
#else
  			//transfer_cold_gas(p,"Hot",p,"Cold",((float)reheated_mass)/Gal[p].ColdGas);
  		  Gal[p].HotGas += reheated_mass;
  		  Gal[p].MetalsHotGas += metallicity * reheated_mass;
#endif
  		}
  		else
  		{
  			if(Gal[p].Type ==1)
  				merger_centre=centralgal;
  			else if(Gal[p].Type ==2)
  				merger_centre=Gal[p].CentralGal;

  			if(HotGasOnType2Galaxies==0) //if no hot gas in type 2's, type 2's share gas between 0 or 1 or 0 and 1
  				dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
  			else if(HotGasOnType2Galaxies==1) //if hot gas in type 2's, type 2's and 1's both share gas between itself or itself and mergercentre
  					dis=separation_gal(merger_centre,p)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

  		  //compute shate of reheated mass
  			if ((dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1 && HotGasOnType2Galaxies==0) ||
  					(dis<Gal[merger_centre].Rvir && HotGasOnType2Galaxies==1))
  			{
  				//mass that remains on type1 (the rest goes to type 0) for reheat - massremain, for eject - ejected mass
  				massremain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
  				ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;

  				if (massremain > reheated_mass)
  					massremain = reheated_mass;
  			}
  			else
  				massremain=reheated_mass;


  			//needed due to precision issues, since we first remove massremain and then (reheated_mass-massremain)
  			//from the satellite into the type 0 and type 1 the fraction might not add up on the second call
  			//since Gal[p].ColdGas is a float and reheated_mass & massremain are doubles
  			if((massremain + reheated_mass)>Gal[p].ColdGas)
  				massremain=Gal[p].ColdGas-reheated_mass;


  			//if(Gal[p].Type !=2)
  			//{
        //transfer massremain
#ifndef H2_AND_RINGS
  			if(HotGasOnType2Galaxies==0) //tranfer to itself if type 1, merger centre if type 2
  				transfer_gas(Gal[p].CentralGal,"Hot",p,"Cold",massremain/Gal[p].ColdGas,"update_from_feedback", __LINE__);
  			else if(HotGasOnType2Galaxies==1) //tranfer to itself
  				transfer_gas(p,"Hot",p,"Cold",massremain/Gal[p].ColdGas,"update_from_feedback", __LINE__);
#else
  			if(HotGasOnType2Galaxies==0) //tranfer to itself if type 1, merger centre if type 2
  			{
  				//transfer_cold_gas(Gal[p].CentralGal,"Hot",p,"Cold",massremain/Gal[p].ColdGas);
  				Gal[Gal[p].CentralGal].HotGas += massremain;
  				Gal[Gal[p].CentralGal].MetalsHotGas += metallicity * massremain;
  			}
  			else if(HotGasOnType2Galaxies==1) //tranfer to itself
  			{
  				//transfer_cold_gas(p,"Hot",p,"Cold",massremain/Gal[p].ColdGas);
  				Gal[p].HotGas += massremain;
  				Gal[p].MetalsHotGas += metallicity * massremain;
  			}
#endif

  			//transfer to the central galaxy
  			if (reheated_mass > massremain)
  				if(Gal[p].ColdGas > 0.) //if the reheat to itself, left cold gas below limit do not reheat to central
  				{
  					//transfer reheated_mass-massremain from galaxy to the type 0
#ifndef H2_AND_RINGS
  					if(HotGasOnType2Galaxies==0) //tranfer to type 0
  						transfer_gas(centralgal,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas,"update_from_feedback", __LINE__);
  					else if(HotGasOnType2Galaxies==1) //tranfer to merger centre
  						transfer_gas(merger_centre,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas,"update_from_feedback", __LINE__);
#else
  					if(HotGasOnType2Galaxies==0) //tranfer to type 0
  					{
  						//transfer_cold_gas(centralgal,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas);
  						Gal[centralgal].HotGas += (reheated_mass-massremain);
  						Gal[centralgal].MetalsHotGas += metallicity * (reheated_mass-massremain);
  					}
  					else if(HotGasOnType2Galaxies==1) //tranfer to merger centre
  					{
  						//transfer_cold_gas(merger_centre,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas);
  						Gal[merger_centre].HotGas += (reheated_mass-massremain);
  						Gal[merger_centre].MetalsHotGas += metallicity * (reheated_mass-massremain);
  					}
#endif
  				}

  		//} //Gal[p].Type !=2

  		}//types

  	}//if(Gal[p].ColdGas > 0.)


    mass_checks("update_from_feedback #2",p);


    //DO EJECTION OF GAS
    if ( (Gal[Gal[p].CentralGal].HotGas > 0. && HotGasOnType2Galaxies==0) ||
    		 (Gal[p].HotGas > 0. && HotGasOnType2Galaxies==1) )
    {

    	if(HotGasOnType2Galaxies==0)
    	{
    		if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
    			ejected_mass = Gal[Gal[p].CentralGal].HotGas;  //either eject own gas or merger_centre gas for ttype 2's

    		fraction=((float)ejected_mass)/Gal[Gal[p].CentralGal].HotGas;
    	}
    	else if(HotGasOnType2Galaxies==1)
    	{
    		if (ejected_mass > Gal[p].HotGas && HotGasOnType2Galaxies==1)
    			ejected_mass = Gal[p].HotGas;  //always eject own gas

    		fraction=((float)ejected_mass)/Gal[p].HotGas;
    	}


    	if (Gal[Gal[p].CentralGal].Type == 1)
    	{
    		/* If type 1, or type 2 orbiting type 1 near type 0 */
    		if (EjectionRecipe == 0)
    		{
#ifndef H2_AND_RINGS
    			if (dis < Gal[centralgal].Rvir)
    				transfer_gas(centralgal,"Hot",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
    			else
    				transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
#else
    			if (dis < Gal[centralgal].Rvir)
    				transfer_hot_gas(centralgal,"Hot",Gal[p].CentralGal,"Hot",fraction);
    			else
    				transfer_hot_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
#endif
    		}
    		else if (EjectionRecipe == 1)
    		{
#ifndef H2_AND_RINGS
    			if(HotGasOnType2Galaxies==0)
    				transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
    			else if(HotGasOnType2Galaxies==1)
    				transfer_gas(Gal[p].CentralGal,"Ejected",p,"Hot",fraction,"update_from_feedback", __LINE__);
#else
    			if(HotGasOnType2Galaxies==0)
    				transfer_hot_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
    			else if(HotGasOnType2Galaxies==1)
    				transfer_hot_gas(Gal[p].CentralGal,"Ejected",p,"Hot",fraction);
#endif
    		}
    	}
    	else // If galaxy type 0 or type 2 merging into type 0
    	{
#ifndef H2_AND_RINGS
    		if(HotGasOnType2Galaxies==0)
    			transfer_gas(centralgal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
    		else if(HotGasOnType2Galaxies==1)
    			transfer_gas(centralgal,"Ejected",p,"Hot",fraction,"update_from_feedback", __LINE__);
#else
    		if(HotGasOnType2Galaxies==0)
    			transfer_hot_gas(centralgal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
    		else if(HotGasOnType2Galaxies==1)
    			transfer_hot_gas(centralgal,"Ejected",p,"Hot",fraction);
#endif
    	}

    	/* if(Gal[p].HotGas < PRECISION_LIMIT)
    	  	Gal[p].HotGas = 0.;
    	 if(Gal[centralgal].HotGas < PRECISION_LIMIT)
    		 Gal[centralgal].HotGas = 0.;
    	 if(Gal[Gal[p].CentralGal].HotGas < PRECISION_LIMIT)
    		 Gal[Gal[p].CentralGal].HotGas = 0.;*/

    }//(Gal[Gal[p].CentralGal].HotGas > 0.)

  } //if(FeedbackRecipe == 0 || FeedbackRecipe == 1)

#ifdef H2_AND_RINGS
  update_h2fraction(p);
#endif


}


#ifdef H2_AND_RINGS
void update_h2fraction(int p)
{
	int j;
	//the central stellar surface density //
	double sigmastar0 = Gal[p].DiskMassr[0]/radius[0]/radius[0]/M_PI*0.01*Hubble_h;

	Gal[p].H2fraction=0;

  for(j=0;j<RNUM;j++)
  {
  	if(H2FractionRecipe==0 || H2FractionRecipe==1)
  	{
  		double sigmahr, metallicityr;

  		metallicityr = metals_total(Gal[p].MetalsColdGasr[j])/Gal[p].ColdGasr[j]/0.02;
			 //if(metallicityr*Clumpingfactor<0.5) metallicityr=0.5/Clumpingfactor;

			if(metallicityr<0.01)
				metallicityr=0.01;
			if(j==0)
				sigmahr = Gal[p].ColdGasr[j] / (M_PI* radius[j]*radius[j])/Warmphasefactor;
			else
				sigmahr = Gal[p].ColdGasr[j] / (M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]))/Warmphasefactor;

			/* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
			sigmahr=sigmahr*0.01*Hubble_h;

			if(metallicityr<1.0)
				sigmahr=sigmahr*Clumpingfactor*pow((1.0/metallicityr),0.7);
			//sigmahr=sigmahr*Clumpingfactor;

			if(H2FractionRecipe==0)
			{
				double tau, khi, s;
				khi=3.1*1.0*(1+3.1*pow(metallicityr,0.365))/4.1;
				tau=0.066*sigmahr*metallicityr;
				s=log(1+0.6*khi+0.01*khi*khi)/0.6/tau;
				if(s<2.0)
					Gal[p].H2fractionr[j]=1-0.75*s/(1+0.25*s);
				else
					Gal[p].H2fractionr[j]=0.0;
				//if(Gal[p].H2fractionr[j]<0.01) Gal[p].H2fractionr[j]=0.01;
			}
			else if(H2FractionRecipe==1)
			{
				/*convert to log10*/
				metallicityr = log10(metallicityr);    sigmahr=log10(sigmahr);
				Gal[p].H2fractionr[j]=H2frac(sigmahr,metallicityr);
			}
  	}
  	else if(H2FractionRecipe==2)
  	{
   		double sigmahr, sigmastarr;
       if(j==0)
       	{
       		sigmahr = Gal[p].ColdGasr[j] / (M_PI* radius[j]*radius[j])/Warmphasefactor;
       		sigmastarr=Gal[p].DiskMassr[j] / (M_PI* radius[j]*radius[j]);
       	}
   	  else
   	  	{
   	  		sigmahr = Gal[p].ColdGasr[j] / (M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]))/Warmphasefactor;
   	  		sigmastarr = Gal[p].DiskMassr[j] / (M_PI*(radius[j]*radius[j]-radius[j-1]*radius[j-1]));
   	  	}
       sigmahr*=(0.01*Hubble_h*1.0);   sigmastarr*=(0.01*Hubble_h);	//from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
       Gal[p].H2fractionr[j]=1.38e-3*pow(sigmahr*(sigmahr+0.1*sqrt(sigmastar0*sigmastarr)),0.92);
       //Gal[p].H2fractionr[j]=6.81e-3*pow(sigmahr*(sigmahr+0.1*sqrt(sigmastar0*sigmastarr)),0.80);
       if(Gal[p].H2fractionr[j]<1.0e-8) Gal[p].H2fractionr[j]=0.0;
       	else Gal[p].H2fractionr[j]=1/(1+1/Gal[p].H2fractionr[j]);
  	}
  	else Gal[p].H2fractionr[j]=0;

  	Gal[p].H2fraction += Gal[p].H2fractionr[j] * Gal[p].ColdGasr[j]/Gal[p].ColdGas;

  }

  if(Gal[p].ColdGas<1.0e-7) Gal[p].H2fraction=0.0;




}
#endif


//Age in Mpc/Km/s/h - code units
void update_massweightage(int p, double stars, double time)
{
	int outputbin;
	double age;

	for(outputbin = 0; outputbin < NOUT; outputbin++)
	  {
	  	age = time - NumToTime(ListOutputSnaps[outputbin]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
	   	Gal[p].MassWeightAge[outputbin] += age * stars;
#else
	   	Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
#endif
	  }
}

/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */

void check_disk_instability(int p)
{

  double Mcrit, fraction, stars, diskmass;
#ifdef H2_AND_RINGS
  double rstar, vmax;
  int j;
#endif
/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used. */

  diskmass = Gal[p].DiskMass;

#ifndef H2_AND_RINGS
  /* check stellar disk -> eq 34 Guo2010*/
  if (Gal[p].Type != 0)
    Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].StellarDiskRadius / G;
  else
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].StellarDiskRadius / G;
#else
   if(diskmass<1.0e-6) rstar=0.5*radius[0];
    else
   	{
   		rstar=0.5*radius[0]*Gal[p].DiskMassr[0];
   		for(j=1;j<RNUM;j++) rstar+=(0.5*(radius[j-1]+radius[j])*Gal[p].DiskMassr[j]);
   		rstar=rstar/diskmass/2.0;      //2.0=mean radius/scale length for exponential disk
   	}
   if (Gal[p].Type != 0)
  	 vmax=Gal[p].InfallVmax;
   	else vmax=Gal[p].Vmax;
     Mcrit = vmax * vmax * rstar / G;
#endif

  stars = diskmass - Mcrit;
  fraction = stars / diskmass;

  /* add excess stars to the bulge */
  if(stars > 0.0) {
    /* to calculate the bulge size */
    update_bulge_from_disk(p,stars);
#ifndef H2_AND_RINGS
    transfer_stars(p,"Bulge",p,"Disk",fraction);
#else
    Gal[p].BulgeMass += stars;
    Gal[p].DiskMass  -= stars;
#endif
#ifndef POST_PROCESS_MAGS
    double Lumdisk;
    int outputbin, j;
#ifdef OUTPUT_REST_MAGS
    for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++)
      {
      	Lumdisk = Gal[p].Lum[j][outputbin]-Gal[p].LumBulge[j][outputbin];
      	Gal[p].LumBulge[j][outputbin] += fraction * Lumdisk;
      	Lumdisk = Gal[p].YLum[j][outputbin]-Gal[p].YLumBulge[j][outputbin];
      	Gal[p].YLumBulge[j][outputbin] += fraction * Lumdisk;
      }
    }
#endif
#ifdef COMPUTE_OBS_MAGS
    for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++)
      {
      	Lumdisk = Gal[p].ObsLum[j][outputbin]-Gal[p].ObsLumBulge[j][outputbin];
      	Gal[p].ObsLumBulge[j][outputbin] += fraction * Lumdisk;
      	Lumdisk = Gal[p].ObsYLum[j][outputbin]-Gal[p].ObsYLumBulge[j][outputbin];
      	Gal[p].ObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#ifdef OUTPUT_MOMAF_INPUTS
      	Lumdisk = Gal[p].dObsLum[j][outputbin]-Gal[p].dObsLumBulge[j][outputbin];
      	Gal[p].dObsLumBulge[j][outputbin] += fraction * Lumdisk;
      	Lumdisk = Gal[p].dObsYLum[j][outputbin]-Gal[p].dObsYLumBulge[j][outputbin];
      	Gal[p].dObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#endif
      }
    }
#endif
#endif

    if ((Gal[p].BulgeMass > 1e-8 && Gal[p].BulgeSize == 0.0)||
    		(Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize >1e-8))
    {
    	char sbuf[1000];
    	sprintf(sbuf, "bulgesize wrong in diskinstablility.c \n");
    	terminate(sbuf);
    }
  }
  
}


/** @brief Introduced in Guo2010 to track the change in size of bulges
 *         after their growth due to disk instabilities. */

void update_bulge_from_disk(int p, double stars)
{      
  double bulgesize, diskmass, fint, massfrac, orisize;
  int  j;
#ifdef H2_AND_RINGS
  double dmass, metallicity;
#endif

#ifndef H2_AND_RINGS

/** @brief Updates bulge from disk instability -> stars represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to the occupied in the disk. */


  orisize=Gal[p].BulgeSize; //remove, not used
  diskmass=(Gal[p].DiskMass);

  /* alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   * the existing and newly formed bulges are concentric)*/
  fint=4.0;

  /* update the stellardisk spin due to the angular momentum transfer
   * from disk to bulge changing the specific angular momentum for disk stars.
   * This should be done on the main routine, as this is update bulge.*/
  massfrac=stars/diskmass;
  for (j = 0; j <3 ; j++)
    Gal[p].StellarSpin[j]=Gal[p].StellarSpin[j]/(1-massfrac);

  /* update disksize done, disk mass is automatically given by total-bulge*/

//GET BULGE SIZE - Eq. 35 in Guo2010
  /* if previous Bulge Mass = 0
     -> bulge size is given directly from newly formed bulge */
  if(Gal[p].BulgeMass <1.e-9) {
    /* size of newly formed bulge, which consists of the stellar mass
     * transfered from the disk. This is calculated using bulge_from_disk
     * which receives Delta_M/DiskMass and returns Rb/Rd. From eq 35 and
     * since DiskMass=2PISigma(Rd)^2 we see that Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd),
     * so function bulge_from_disk avoids calculating the slow "ln" function */
    bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].StellarDiskRadius/3.;
    Gal[p].BulgeSize=bulgesize;

  }      
  else {
    bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].StellarDiskRadius/3.;
    /* combine the old with newly formed bulge and calculate the
     * bulge size assuming energy conservation as for mergers but
     * using alpha=2. - eq 33 */
    Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/
      (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));
  }



  // TODO - check why we need it
  if((Gal[p].BulgeMass + stars > 1.e-8 && Gal[p].BulgeSize == 0.0)
     || (Gal[p].BulgeMass + stars == 0 && Gal[p].BulgeSize > 1.e-8))
    {
      char sbuf[1000];
      sprintf
	(sbuf,
	 "bulgesize wrong in disk instablility. diskmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f masstransfer %f trassize %f, oribulgesize %f\n",
	 Gal[p].DiskMass, Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, Gal[p].GasDiskRadius,
	 Gal[p].StellarDiskRadius, stars, bulgesize, orisize);
      terminate(sbuf);
    }
#else

  diskmass=Gal[p].DiskMass+Gal[p].ColdGas;

  /* alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   * the existing and newly formed bulges are concentric)*/
  fint=4.0;

  /* update the stellardisk spin due to the angular momentum transfer
   * from disk to bulge changing the specific angular momentum for disk stars.
   * This should be done on the main routine, as this is update bulge.*/
  massfrac=stars/diskmass;
  if(massfrac<1.0)
  	for (j = 0; j <3 ; j++)
  		Gal[p].DiskSpin[j]=Gal[p].DiskSpin[j]/(1-massfrac);

  /* update disksize done, disk mass is automatically given by total-bulge*/

  /*size of new formed bulge, which consist of the stellar mass trasfered from the disk*/
  /*combine the old bulge with the new materials and caculate the bulge size assuming energy conservation */
  dmass=stars;
  j=0; //avoid non-definded j if dmass<1e-6
  if(dmass>1.0e-6)
  {
  	for(j=0;j<RNUM;j++)
  	{
  		if(dmass>Gal[p].DiskMassr[j])
  		{
  			dmass-=Gal[p].DiskMassr[j];
  			Gal[p].MetalsBulgeMass+=Gal[p].MetalsDiskMassr[j];
  			Gal[p].MetalsDiskMass-=Gal[p].MetalsDiskMassr[j];
  			Gal[p].DiskMassr[j]=0;
  			Gal[p].MetalsDiskMassr[j]=0;
  		}
  		else break;
  	}
  	if(j==RNUM)
  		bulgesize=radius[RNUM-1];
  	else
  	{
  		if(j==0) bulgesize=dmass/Gal[p].DiskMassr[j]*radius[j];
  		else bulgesize=dmass/Gal[p].DiskMassr[j]*radius[j]+(1-dmass/Gal[p].DiskMassr[j])*radius[j-1];
  	}
  }
  else bulgesize=0.5*radius[0];

  /*update the stellar information in the last ring*/
  metallicity = metals_total(Gal[p].MetalsDiskMassr[j])/Gal[p].DiskMassr[j];
  Gal[p].MetalsBulgeMass+=metallicity*dmass;
  Gal[p].MetalsDiskMass -=metallicity*dmass;
  Gal[p].DiskMassr[j]-=dmass;
  Gal[p].MetalsDiskMassr[j]-=metallicity*dmass;

  if(Gal[p].BulgeMass <1.e-9) Gal[p].BulgeSize=bulgesize;
        else Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/(Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));

    if ((Gal[p].BulgeMass+stars > 1.e-8 && Gal[p].BulgeSize == 0.0)||(Gal[p].BulgeMass+stars == 0 && Gal[p].BulgeSize >1.e-8))
      {
        printf("bulgesize wrong in disk instablility. stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f, masstransfer %f trassize %f, oribulgesize %f\n",Gal[p].DiskMass,Gal[p].BulgeMass, Gal[p].BulgeSize,Gal[p].ColdGas,stars,bulgesize);
        exit(0);
      }

#endif



}


/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. */
double bulge_from_disk(double frac)
{
  double x1,x2,x0,value;
/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. The bulge is assumed
 *         to form with the same size. avoid doing "ln" from eq 35*/
  x1=0.0;
  x2=1.;
  while ((func_size(x2,frac) * func_size(x1,frac))>0) {
    x1=x2;
    x2=x2*2;
  }
  x0=x1+(x2-x1)/2.;
  value=func_size(x0,frac);
  if (value < 0) 
    value = -value;

  while(value>0.00001) {
    if(func_size(x0,frac)*func_size(x2,frac)>0)
      x2=x0;
    else
      x1=x0;
    x0=x1+(x2-x1)/2.;
    value=func_size(x0,frac);
    if (value < 0) 
      value = -value;
  }
    
  return x0;
}


double func_size(double x, double a)
{
  return  exp(-x)*(1+x)-(1-a);
}  

#ifdef H2_AND_RINGS
double H2frac(double logsigmah, double metallicity )
{
	int i,j;
	double logNHtot[LENSIGMAH],lgZ[LENZ],mf,mf1,mf2;
	for ( i=0,logNHtot[0]=-1;i<(LENSIGMAH-1);i++ ) logNHtot[i+1]=logNHtot[i]+0.05;
	for ( j=0,lgZ[0]=-2;j<(LENZ-1);j++ ) lgZ[j+1]=lgZ[j]+0.25;

	if ( logsigmah<logNHtot[0] ) logsigmah=logNHtot[0];
	if ( logsigmah>logNHtot[i-1] ) logsigmah=logNHtot[i-1];
	for ( i=0;logsigmah > logNHtot[i + 1];i++ );

	if ( metallicity<lgZ[0] ) metallicity=lgZ[0];
	if ( metallicity>lgZ[j-1] ) metallicity=lgZ[j-1];
	for ( j=0;metallicity>lgZ[j+1];j++ );

	mf1=h2frac[i][j]+ ( h2frac[i][j+1]-h2frac[i][j] ) * ( metallicity-lgZ[j] ) / ( lgZ[j+1]-lgZ[j] );
	mf2=h2frac[i+1][j]+ ( h2frac[i+1][j+1]-h2frac[i+1][j] ) * ( metallicity-lgZ[j] ) / ( lgZ[j+1]-lgZ[j] );
	mf=mf1+ ( mf2-mf1 ) * ( logsigmah-logNHtot[i] ) / ( logNHtot[i+1]-logNHtot[i] );

	return ( mf );
}
#endif


