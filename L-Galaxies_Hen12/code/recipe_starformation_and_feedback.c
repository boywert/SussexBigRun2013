#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
 * There are 3 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - Delucia 2004 retention
 * or ejection scheme (0),
 *
 * 1 - \f$\Delta m_{\rm reheated}=\epsilon_{\rm disk}\Delta m_{\star}\f$\n
 *     The amount of energy released by supernova during the formation of
 *     \f$\Delta m_{\star}\f$ stars is
 *     \f$\Delta E_{\rm SN}=0.5\,\epsilon_{\rm halo}\,\Delta m_\star V_{\rm SN}^2\f$,
 *     \f$V_{\rm SN}=630$\,km\,s$^{-1}\f$. Any excess energy left over from
 *     reheating the cold gas is used to eject a mass of gas from the galaxy:\n
 *     \f$\Delta m_{\rm ejected}=\left(\epsilon_{\rm halo} \frac{V_{\rm
 *     SN}^2}{V_{\rm vir}^2} -\epsilon_{\rm disk}\right) \Delta m_{\star}\f$,\n
 *     (Eqs. 17, 18, 19 & 20 Croton2006)(FeedbackRecipe = 1) based on moderate
 *     Delucia 2004 ejection scheme.
 *
 * 2 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackRecipe = 2)
 *     same as FeedbackRecipe = 1 * Vmax dependence.
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive
 * gas from their own satellites when these are outside Rvir of the type 0.
 * */

/*TO DO - eliminate duplications by merging CentralVvir and MergeCentralVvir
 *       and centralgal and Gal[p].centragal into single variables with value
 *       defined by a conditions on dis */
/** @brief Main recipe, calculates the fraction of cold gas turned into stars due
  *        to star formation; the fraction of mass instantaneously recycled and
  *        returned to the cold gas; the fraction of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   */
void starformation_and_feedback(int p, int centralgal, double time, double dt, int nstep) 
{
	/*! Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt*/
  double tdyn, strdot, stars, reheated_mass, ejected_mass, fac, cold_crit,
         metallicitySF, CentralVvir, MergeCentralVvir, EjectVmax, EjectVvir,
         SN_Energy, Reheat_Energy;

  
  //standard star formation law (Croton2006, Delucia2007, Guo2010)
  if(StarFormationRecipe == 0)
  {
  	if(Gal[p].Type == 0)
  		tdyn = Gal[p].GasDiskRadius / Gal[p].Vmax;
  	else
  		tdyn = Gal[p].GasDiskRadius / Gal[p].InfallVmax;

  	if(Gal[p].Type == 0)
  		cold_crit = 0.19 * Gal[p].Vmax * Gal[p].GasDiskRadius;
	  else
	  	cold_crit = 0.19 * Gal[p].InfallVmax * Gal[p].GasDiskRadius;

  	if(Gal[p].ColdGas > cold_crit)
  		strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);
  	else
  		strdot = 0.0;
  }

  //No dependce on Tdyn
  if(StarFormationRecipe == 1)
    {

	  if(Gal[p].ColdGas > 1.e-7)
		strdot = SfrEfficiency * (Gal[p].ColdGas) / (2.e-5) * pow(Gal[p].Vvir / SfrLawPivotVelocity, SfrLawSlope);
	  else
		strdot = 0.0;
    }

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

/* SN FEEDBACK RECIPES */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
  MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
   
  mass_checks("recipe_starform #0",p);
  mass_checks("recipe_starform #0.1",centralgal);

  /* Feedback depends on the circular velocity of the host halo
   * Guo2010 - eq 18 & 19*/
  if(FeedbackRecipe == 0) {
    if (Gal[Gal[p].CentralGal].Type == 0)
      reheated_mass = FeedbackReheatingEpsilon * stars*
	(.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
    else
      reheated_mass = FeedbackReheatingEpsilon * stars*
	(.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));
  }

  //Make sure that the energy used in reheat does not exceed the SN energy (central subhalo Vvir used)
  if (FeedbackRecipe == 0 || FeedbackRecipe == 1) {
    if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * EnergySNcode))
      reheated_mass = stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
  }

  /*  cant use more cold gas than is available! so balance SF and feedback */
  if((stars + reheated_mass) > Gal[p].ColdGas) {
    fac = Gal[p].ColdGas / (stars + reheated_mass);
    stars *= fac;
    reheated_mass *= fac;
  }  

  /* Determine ejection (for FeedbackRecipe 2 we have the dependence on Vmax)
   * Guo2010 - eq 22
   * Note that satellites can now retain gas and have their own gas cycle*/

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

  /*  update the star formation rate */
#ifdef SAVE_MEMORY
  /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
  Gal[p].Sfr += stars / (dt * STEPS);
#else
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    if(Halo[halonr].SnapNum == ListOutputSnaps[outputbin]) {
      Gal[p].Sfr[outputbin] += stars / (dt * STEPS);
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
#ifdef METALS
    //metallicitySF=Gal[p].MetalsColdGas.type2/Gal[p].ColdGas;
	metallicitySF= metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
#else
    metallicitySF=Gal[p].MetalsColdGas/Gal[p].ColdGas;
#endif
  else
    metallicitySF=0.;

  /* Metallicity is not a single number when METALS is defined, and it is not needed 
   * in this routine, so I have removed it - PAT */
  //update_from_star_formation(p, time, stars, metallicity);
  if (stars > 0.) 
	update_from_star_formation(p, stars, dt, nstep); //ROB: Added 'dt' and 'nstep' as variables, so that it can be used by 'update_yields', which is called from 'update_from_star_formation'
  
  mass_checks("recipe_starform #2",p);
  mass_checks("recipe_starform #2.1",centralgal);

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */
  if (reheated_mass + ejected_mass > 0.)
  update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
  
#ifndef POST_PROCESS_MAGS
  /*  Update the luminosities due to the stars formed */
  if (stars > 0.0)
    add_to_luminosities(p, stars, time, metallicitySF);
#endif

  if(TrackDiskInstability)
    if(Gal[p].DiskMass > 0.0)
      check_disk_instability(p);

  if (DiskRadiusMethod == 2)
    get_stellar_disk_radius(p);

}


/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
void update_from_star_formation(int p, double stars, double dt, int nstep) // //ROB: Two variables, 'dt' and 'nstep', added, for use in update_yields()
{
  int i;
  double fraction;
#ifndef YIELDS
  double stars_nett;
#endif

  if(Gal[p].ColdGas <= 0. || stars <= 0.) {
    printf("update_from_star_formation: Gal[p].ColdGas <= 0. || stars <= 0.\n");
    exit(0);
  }

#ifndef YIELDS
  stars_nett=(1 - RecycleFraction) * stars; //ROB: No longer an assumed instantaneous recycled fraction. Mass is returned over time via SNe and AGB winds.
#endif

  /* Update the Stellar Spin when forming stars */ //ROB: This could be moved to after the yields, total metals and total ejected masses are updated.
#ifdef YIELDS
  if (Gal[p].DiskMass+stars > 1.e-8)
    for (i = 0; i < 3; i++)
      Gal[p].StellarSpin[i]=((Gal[p].StellarSpin[i])*
			     (Gal[p].DiskMass)+
			     stars*Gal[p].GasSpin[i])/
	(Gal[p].DiskMass+stars);
#else
  if (Gal[p].DiskMass+stars_nett > 1.e-8)
    for (i = 0; i < 3; i++)
      Gal[p].StellarSpin[i]=((Gal[p].StellarSpin[i])*
			     (Gal[p].DiskMass)+
			     stars_nett*Gal[p].GasSpin[i])/
	(Gal[p].DiskMass+stars_nett);
#endif //YIELDS
  
    /*  Update Gas and Metals from star formation */
  
  mass_checks("update_from_star_formation #0",p);

  //transfer_gas_to_stars(p,"Disk",p,"Cold",stars_nett/Gal[p].ColdGas);

#ifdef YIELDS
  fraction=stars/Gal[p].ColdGas;
#else
  fraction=stars_nett/Gal[p].ColdGas; //ROB: Old calculation for fraction of un-recycled ColdGas mass that was formed into stars
#endif

#ifdef STAR_FORMATION_HISTORY
#ifdef YIELDS
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars; //ROB: Now, all SF gas is put in SFH array ("recycled' mass will return to gas phase over time)
  Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);
  Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin] = elements_add(Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin],Gal[p].ColdGas_elements,fraction);
  //if (Gal[p].SnapNum == 16) {printf("SFH H_added = %f\n", Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin].H + (Gal[p].ColdGas_elements.H)*fraction);}
#else
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars_nett; //ROB: Add amount of stars formed to SFH history of the disk. (NOTE: ALL SF OCCURS IN THE DISK. sfh_BulgeMass only increases when stars are transferred to the bulge before they explode)
  Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);
#endif //YIELDS
#endif //STAR_FORMATION_HISTORY

#ifdef YIELDS
  Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Gal[p].MetalsColdGas,fraction);
  Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsColdGas,-fraction);
  Gal[p].DiskMass_elements=elements_add(Gal[p].DiskMass_elements,Gal[p].ColdGas_elements,fraction);
  Gal[p].ColdGas_elements=elements_add(Gal[p].ColdGas_elements,Gal[p].ColdGas_elements,-fraction);
  Gal[p].DiskMass += stars;
  Gal[p].ColdGas -= stars;
#else
  Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Gal[p].MetalsColdGas,fraction);
  Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsColdGas,-fraction);
  Gal[p].DiskMass += stars_nett;
  Gal[p].ColdGas -= stars_nett;
#endif //YIELDS

  mass_checks("update_from_star_formation #1",p);

  /* Formation of new metals - instantaneous recycling approximation - only SNII
   * Also recompute the metallicity of the cold phase.
   * TODO In bursts, gas is NOT later transferred from cold to hot - should it be? */
  /* TODO Note that FracZtoHot is broken as there may not be any hot gas.  If metals
   * are given to the hot phase then mass needs to be also. */

  /*DELAYED ENRICHMENT AND MASS RETURN, USING INDIVIDUAL ELEMENT YIELDS: No fixed yield or recycling fraction anymore:*/
  if (FeedbackRecipe == 0 || FeedbackRecipe == 1) {
#ifdef METALS
#ifdef YIELDS
	update_yields(p, dt, nstep);
#else
    Gal[p].MetalsColdGas.type2 += Yield * (1.0 - FracZtoHot) * stars;
    Gal[Gal[p].CentralGal].MetalsHotGas.type2 += 
      Yield * FracZtoHot * stars;// =0 not used
#endif //YIELDS

#ifdef METALS_SELF
    Gal[p].MetalsHotGasSelf.type2 += Yield * FracZtoHot * stars;
#endif
#else
    Gal[p].MetalsColdGas += Yield * (1.0 - FracZtoHot) * stars;
    Gal[Gal[p].CentralGal].MetalsHotGas += Yield * FracZtoHot * stars;// FracZtoHot=0 not used
#ifdef METALS_SELF
    Gal[p].MetalsHotGasSelf += Yield * FracZtoHot * stars;
#endif
#endif //METALS
  }

  if (DiskRadiusMethod == 2)
    get_stellar_disk_radius(p);

}


/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection.
 * TODO - Should ejected_mass really be modified in reheat section?
 * TODO - Check ejection recipe; it looks wrong */
void update_from_feedback(int p, int centralgal, double reheated_mass, 
			  double ejected_mass)
{
  double dis;
  double massremain;
  double fraction;

  /* Check that there is cold gas: if not then something wrong.
   * This routine should not be called in these circumstances. */
  if (Gal[p].ColdGas < 1.0e-8) {
    printf("Gal[p].ColdGas < 1e-8 in update_form_feedback\n");
    exit(19);
  }

  if(FeedbackRecipe == 0 || FeedbackRecipe == 1) {
	
    dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

    //REHEAT
    if (dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1) {
      // mass that remains on satellite (the rest goes to type 0)
      // for reheat - massremain, for eject - ejected mass
      massremain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
      /**************** Should this line be here? ********************/
      ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;
    }
    else
    	massremain=reheated_mass;

    if (massremain > reheated_mass)
      massremain = reheated_mass;

    transfer_gas(Gal[p].CentralGal,"Hot",p,"Cold",massremain/Gal[p].ColdGas);

    if (reheated_mass > massremain)
      transfer_gas(centralgal,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas);

    //EJECTION
    //In all these schemes, the ejected mass seems to be taken from Gal[Gal[p].CentralGal].HotGas,
    //So test for this being positive and only eject gas if it is
    
    if (Gal[Gal[p].CentralGal].HotGas > 0.)
    {
    	if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
    	{
    		ejected_mass = Gal[Gal[p].CentralGal].HotGas;
    	}

    	fraction=ejected_mass/Gal[Gal[p].CentralGal].HotGas;
	
    	if (Gal[Gal[p].CentralGal].Type == 1)
    	{
    		/* If type 1, or type 2 orbiting type 1 near type 0 */
    		if (dis < Gal[centralgal].Rvir)
    		{
    			if (EjectionRecipe == 0)
    			{
    				// In Delucia2007 gas ejected from satellites goes to hot of type 0
    				transfer_gas(centralgal,"Hot",Gal[p].CentralGal,"Hot",fraction);
    			}
    			else if (EjectionRecipe == 1)
    			{
    				// In Guo2010 ejected gas goes to ejected of type 1
    				transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
    			}
    		}
    		else
    		{
    			// If it is far way from type 0 goes to type 1 in both implementations
    			transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
    		}
    	}
    	else // If central galaxy is a type 0 - it takes everything
    	{
    		transfer_gas(centralgal,"Ejected",Gal[p].CentralGal,"Hot",fraction);
    	}
    }
  } //if(FeedbackRecipe == 0 || FeedbackRecipe == 1)
}


/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */

void check_disk_instability(int p)
{

  double Mcrit, fraction, stars, diskmass;
  double Lumdisk;
  int outputbin, j;

/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used. */


  /* check stellar disk -> eq 34 Guo2010*/
  if (Gal[p].Type != 0)
    Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].StellarDiskRadius / G;
  else
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].StellarDiskRadius / G;
  diskmass = Gal[p].DiskMass;
  stars = diskmass - Mcrit;
  fraction = stars / diskmass;

  /* add excess stars to the bulge */
  if(stars > 0.0) {
    /* to calculate the bulge size */
    update_bulge_from_disk(p,stars);      
    transfer_stars(p,"Bulge",p,"Disk",fraction);
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
    for(outputbin = 0; outputbin < NOUT; outputbin++) {
      for(j = 0; j < NMAG; j++) {
	Lumdisk = Gal[p].Lum[j][outputbin]-Gal[p].LumBulge[j][outputbin];
	Gal[p].LumBulge[j][outputbin] += fraction * Lumdisk;
	Lumdisk = Gal[p].YLum[j][outputbin]-Gal[p].YLumBulge[j][outputbin];
	Gal[p].YLumBulge[j][outputbin] += fraction * Lumdisk;
      }
    }
#endif
#ifdef COMPUTE_OBS_MAGS
    for(outputbin = 0; outputbin < NOUT; outputbin++) {
      for(j = 0; j < NMAG; j++) {
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
	(Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize >1e-8)) {
      printf("bulgesize wrong in diskinstablility.c \n");
      exit(0);
    }
  }
  
}


/** @brief Introduced in Guo2010 to track the change in size of bulges
 *         after their growth due to disk instabilities. */

void update_bulge_from_disk(int p, double stars)
{      
  double bulgesize, diskmass, fint, massfrac, orisize;
  int  j;

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


#ifdef GASRECYCLE

GASRECYLE option not checked

#ifdef GALAXYTREE  
void  cal_gas_recycle(int ngal);
{
  int ind, p, q;
  double ind=0;
  double fraction;

  for(p=0;p<ngal;p++)
    if (Gal[p].Type != 3)
      ind=ind+1;
  for (p = 0; p < ind; p++) {
    RecGas = 0.0; 
    //Note: metalcold is a global variable
    metalcold=metals_init();
    q = p + NumGals -ind ;
    update_from_recyle(q,ptime,previoustime);
    HaloGal[q].MetalsColdGas=metals_add(HaloGal[q].MetalsColdGas,metalcold,1.);
    HaloGal[q].MetalsStellarMass=metals_add(HaloGal[q].MetalsStellarMass,
      metalcold,-1.);
    HaloGal[q].ColdGas = HaloGal[q].ColdGas + RecGas;
    starmass=HaloGal[q].StellarMass;
    HaloGal[q].StellarMass = HaloGal[q].StellarMass - RecGas;
     
    if (HaloGal[q].StellarMass >= 0.) {  
      if (starmass >0.) {
	//I think that this formulation corrects a bug (BulgeMass was being
	//updated before use as a ratio)
	fraction=RecGas/starmass;
	HaloGal[q].BulgeMass-=HaloGal[q].BulgeMass*fraction;
	HaloGal[q].MetalsBulgeMass=metals_add(HaloGal[q].MetalsBulgeMass,
          metalcold,-fraction);
      }
    }
    else {
      printf("wrong in gas recycle . Stellarmass %f, RecGAs %f\n",starmass,RecGas);
      exit(0);
    }
 
  }
}
#endif


#ifdef GALAXYTREE
double update_from_recyle(int p, double time, double previoustime)
{
  int q,metindex, tabindex,outputbin;
  double f1, f2, fmet1, fmet2, T, ptime, pretime;
  float mass;
  ptime=time;
  pretime=previoustime;
  
  q=HaloGal[p].FirstProgGal;
 
  while (q >= 0 ) {
    update_from_recyle(q, ptime,pretime);
    q=HaloGal[q].NextProgGal;
  }
 
  find_interpolation_point(NumToTime(HaloGal[p].SnapNum),previoustime,  &tabindex, 
			   &f1, &f2 );
 
  if(tabindex < 0) {
    printf("Something strange here (recyletable)....%e\t%d snap %d\n", NumToTime(HaloGal[p].SnapNum)-time, tabindex ,HaloGal[p].SnapNum);
    exit(19);
  }
  T = NumToTime(HaloGal[p].SnapNum-1) - NumToTime(HaloGal[p].SnapNum);
   
#ifdef SAVE_MEMORY
  mass=HaloGal[p].Sfr*T*(Frac[tabindex] * f1+Frac[tabindex+1] * f2); 
#else
  mass=HaloGal[p].Sfr[HaloGal[p].SnapNum]*T*(Frac[tabindex] * f1+Frac[tabindex+1] * f2);          
#endif  
     
  find_interpolation_point(NumToTime(HaloGal[p].SnapNum),time,  &tabindex, 
			     &f1, &f2 );
#ifdef SAVE_MEMORY
  mass=HaloGal[p].Sfr*T * (Frac[tabindex] * f1+Frac[tabindex+1] * f2)-mass;
#else
  mass=HaloGal[p].Sfr[HaloGal[p].SnapNum]*T * (Frac[tabindex] * f1+Frac[tabindex+1] * f2)-mass;
#endif
  metalcold=metals_add(metalcold,
		       HaloGal[p].MetalsColdGas,mass/HaloGal[p].ColdGas);
  if (HaloGal[p].FirstProgGal == -1 && mass > HaloGal[p].StellarMass)
    printf("somthing is wrong here(recycle)");
            
  RecGas=RecGas+mass;
  
} 
#endif
#endif
