#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_mergers.c
 *  @brief Calculates the merging time, the central galaxy (for type 1's),
 *         adds galaxies together, calculates SF from bursts and grows
 *         black holes.
 *
 *
 *       <B>set_merger_center</B> - calculates the central galaxy for type 1's,
 *       since if MERGE01=1 (makefile) type 1's can also merge. Therefore,
 *       they need a merger central galaxy and will also have a merger clock
 *       (needed for millennium two, since due to the high resolution, haloes
 *       are very difficult to disrupt and orbit around forever).
 *
 *
 *
 *       <B>estimate_merging_time</B> sets up a merger clock. Originally this
 *       was done only for type 2 galaxies. The positions of galaxies in the
 *       model are normally given by the position of their dark matter halo.
 *       However, when galaxies become satellites, their dark matter haloes
 *       are stripped by the central object to the point where there is none
 *       left. At this point, the dark matter of the satellite becomes part of
 *       the main halo, but the galaxy's position should continue to evolve
 *       due to the dynamical friction force caused by the dark matter around
 *       it. \n
 *       The code does keep track of the position of the most bounded particle
 *       when the satellite's halo was disrupted, but this is not used to track
 *       the galaxy position. Instead a clock is set, giving the time left until
 *       the satellite mergers with the central galaxy. Before, this was only
 *       done for type 2's (satellites that lost a halo). Guo2010 included the
 *       MERGE01 option, that sets the clock also for type 1's (satellites with
 *       a halo), since for the highest resolution millennium 2, their haloes
 *       can be very small and orbit forever around the central companion.\n
 *       This time is computed using the Chandrasekhar's formula for dynamical
 *       friction, as in Binney & Tremaine 1987:
 *
 *       \f$F_{\rm{df}}=
 *               -\frac{4\pi {\rm{G}}^2 m^2_{\rm{sat}} \ln(\Lambda) \rho B(x)}
 *                     {v^2_{\rm{rel}}}.\f$
 *
 *       Which gives (B&T Eq. 7.26):
 *
 *       \f$t_{df}\approx1.17\frac{V_{\rm{vir}}r^2_{\rm{sat}}}
 *                                {{\rm{G}}\,m_{\rm{sat}}\ln(\Lambda)},\f$
 *
 *       that is afterwards multiplied by 2 (after Delucia2007 to fit the
 *       data). When the merging time reaches zero, the satellite is assumed to
 *       merge with the central galaxy.
 *
 *
 *
 *       <B>deal_with_galaxy_merger</B> deals with the process, according to
 *       the mass fraction of the merger. Major if
 *       \f$M_{\rm{sat}}/M_{\rm{central}}>0.3\f$ and minor otherwise. It calls
 *       - add_galaxies_together - Add the cold and stellar phase of the merged
 *       galaxy to the central one. Also form a bulge at the central galaxy
 *       with the stars from the satellite in a minor merger if
 *       BulgeFormationInMinorMergersOn=1 (Major mergers are dealt later).
 *       - Then calls grow_black_hole - Grows black hole through accretion from
 *       cold gas during mergers (due to the instabilities triggered), as in
 *       Kauffmann & Haehnelt (2000). This is commonly referred as the quasar
 *       mode, main responsible for the black hole growth. After Croton2006 this
 *       mode is active even in minor mergers:
 *       \f$\Delta m_{\rm BH,Q}=M_{\rm{BH,min}}
 *          \frac{f_{\rm BH}(m_{\rm sat}/m_{\rm central})\,m_{\rm cold}}
 *           {1+(280\,\mathrm{km\,s}^{-1}/V_{\rm vir})^2}.\f$

 *       - Finally the burst of star formation due to the merger is treated.
 *           - If StarBurstRecipe = 0 (since Croton2006), the Somerville 2001
 *           model of bursts is used collisional_starburst_recipe(). The burst
 *           can happen for both major and minor mergers, with a fraction of
 *           the added cold gas from the satellite and central being consumed:
 *           \f$\dot{m}_{\star}^{\rm{burst}}
 *             = 0.56 \left(\frac{m_{\rm{sat}}}{m_{\rm{central}}}\right)^{0.7}
 *               m_{\rm{gas}}\f$.
 *           SN Feedback from starformation is computed and the sizes of bulge
 *           and disk followed.
 *
 *       - When a major merger occurs, the disk of both merging galaxies is
 *       completely destroyed to form a bulge. In either type of mergers, the
 *       bulge size is updated using Eq. 33 in Guo2010:
 *       \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
 *          C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+\alpha_{\rm{inter}}
 *          \frac{GM_1M_2}{R_1+R_2}\f$*/

/** @brief Calculates the central galaxies for type 1's (needed if MERGE01=1). */

int set_merger_center(int fofhalo)
{
  /** @brief If MERGE01=1 (Guo2010), get id of central galaxy, since type 1's
   *         can have their merger clock started before they become type 2's
   *         if M_star>M_vir. Introduced for millennium 2, where they can have
   *         very small masses, still be followed and never merge. At this
   *         moment the centre is still not known, reason why this function is
   *         needed. Also, if the type 1 merges, all the type 2's associated with
   *         it will need to know the "new" central galaxy they are merging into. */

  int prog, i, first_occupied, type, halonr, currentgal;
  double lenmax;
  i=0;

  halonr = fofhalo;

  //loop on all the haloes in current FOF group - to find a merger centre
  while(halonr >= 0)
  {
	  lenmax = 0;
	  first_occupied = Halo[halonr].FirstProgenitor;
	  prog = Halo[halonr].FirstProgenitor;

	  /* If the main progenitor of the current halo had no galaxies,
	   * set first_occupied to the most massive progenitor. */
	  if(prog >= 0)
	    {
	      	  if(HaloAux[prog].NGalaxies == 0)
	      		  while(prog >= 0)
	      		  {
	      			for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)
	      			  {
	      				type = HaloGal[currentgal].Type;

	      				if(type == 0 || type == 1)
	      				  {
	      					if(Halo[prog].Len > lenmax)
	      					  {
	      						lenmax = Halo[prog].Len;
	      						first_occupied = prog;
	      					  }
	      				  }
	      				currentgal = HaloGal[currentgal].NextGalaxy;
	      			  }
	      			prog = Halo[prog].NextProgenitor;
	      		  }
	    }

      prog = Halo[halonr].FirstProgenitor;

      while(prog >= 0)//loop over all the progenitors
      {
    	  for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)//loop over all the galaxies in a given progenitor
    	  {
    		  type = HaloGal[currentgal].Type;

    		  if(type == 0 || type == 1) // the galaxy is a type 0 or 1?
    			  if(prog == first_occupied) //is the main progenitor?
    				  if(halonr == Halo[halonr].FirstHaloInFOFgroup) //is the main halo?
    					  return currentgal;
    		  currentgal = HaloGal[currentgal].NextGalaxy;
    	  }
    	  prog = Halo[prog].NextProgenitor;
      }

      //if the halo has no galaxies, return 0
      if(i == 0)
    	  if(Halo[halonr].FirstHaloInFOFgroup == halonr)
    		  return i;

      halonr = Halo[halonr].NextHaloInFOFgroup;
  }

  char sbuf[1000];
  sprintf(sbuf, "wrong in finding the central fof %d gal %d\n", fofhalo, currentgal);
  terminate(sbuf);
  return 0;
}

/** @brief Calculates the merging time whenever a galaxy becomes a satellite*/

double estimate_merging_time(int halonr, int mother_halonr, int p)
{
  int central_halonr;
  double coulomb, mergtime, SatelliteMass, SatelliteRadius, MotherHaloRvir;

  /** @brief Binney & Tremaine 1987 - 7.26 merging time for satellites due to
   *         dynamical friction. After Delucia2007 *2, shown to agree with
   *         Kolchin2008 simulations in Delucia2010. This is set when a galaxy
   *         becomes a type 2 or being a type 1 \f$M_{\rm{star}}>M_{\rm{vir}}\f$.
   *         In DeLucia2007 they could only merge into a type 0, now (after
   *         guo2010) they can merge into a type 1. */


  /*  recipe updated for more accurate merging time (see BT eq 7.26),
     now satellite radius at previous timestep is included */
  central_halonr = Halo[Halo[halonr].Descendant].FirstProgenitor;
  if(Gal[p].Type == 1)
    central_halonr=mother_halonr;
  if(central_halonr == halonr)
    {
      terminate("can't be...!\n");
    }


  coulomb = log(Halo[mother_halonr].Len / ((double) Halo[halonr].Len) + 1);

  /*  should include stellar+cold gas in SatelliteMass! */
  SatelliteMass = get_virial_mass(halonr)+(Gal[p].DiskMass+Gal[p].BulgeMass);

  SatelliteRadius = separation_halo(central_halonr,halonr)/(1 + ZZ[Halo[halonr].SnapNum]);

  int j;
  for (j = 0; j < 3; j++)
	Gal[p].DistanceToCentralGal[j] =  wrap(Halo[central_halonr].Pos[j] - Halo[halonr].Pos[j], BoxSize);


  MotherHaloRvir = get_virial_radius(mother_halonr);
  if(SatelliteRadius > MotherHaloRvir)
    SatelliteRadius = MotherHaloRvir;

  if(SatelliteMass > 0.0) {
    mergtime = 1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halonr) /
               (coulomb * G * SatelliteMass); // Binney & Tremaine Eq.7.26

    /* change introduced by Delucia2007 to fit observations */
    mergtime = MergerTimeMultiplier*mergtime;
    //mergtime = 2.*mergtime;
  }
  else
    mergtime = -99999.9;

  return mergtime;

}

/** @brief Deals with all the physics triggered by mergers
 * centralgal is the galaxy at the centre of the group, merger_centralgal
 * is the galaxy onto which galaxy is merging */

void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double deltaT, int nstep)
{

/** @brief Deals with the physics triggered by mergers, according to the mass
 *         fraction of the merger \f$(M_{\rm{sat}}/M_{\rm{central}}><0.3)\f$.
 *         Add the cold and stellar phase of the satellite galaxy to the central
 *         one, form a bulge at the central galaxy with the stars from the
 *         satellite in a minor merger if BulgeFormationInMinorMergersOn=1.
 *         Grows black hole through accretion from cold gas "quasar mode". 
 *         If StarBurstRecipe = 0, the Somerville 2001 model
 *         of bursts is used, SN Feedback from starformation is computed and
 *         the sizes of bulge and disk followed. When a major merger occurs,
 *         the disk of both merging galaxies is completely destroyed to form
 *         a bulge. New stars form of to the bulge*/

  double mi, ma, mass_ratio, Mcstar, Mcgas, Mcbulge, Mpstar, Mpgas,Mpbulge;
  double frac;
#ifdef H2_AND_RINGS
  double rcstar, rcgas; 
  int j;
#endif
#ifdef GALAXYTREE
  int q;

  mass_checks("deal_with_galaxy_merger #0",p);
  mass_checks("deal_with_galaxy_merger #0",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #0",centralgal);

  q = Gal[merger_centralgal].FirstProgGal;
  if(q >= 0)
    {
      while(GalTree[q].NextProgGal >= 0)
        q = GalTree[q].NextProgGal;
        
      GalTree[q].NextProgGal = Gal[p].FirstProgGal;
      
      if(GalTree[q].NextProgGal >= NGalTree)
	{
	  printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		 q, p, GalTree[q].NextProgGal, NGalTree);
	  terminate("problem");
	}
    }

  if(q < 0)
    terminate("may not happen");

  q = GalTree[q].NextProgGal;

  if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
    terminate("inconsistency");

  HaloGal[GalTree[q].HaloGalIndex].MergeOn = 2;

  if(Gal[p].Type == 1)
    HaloGal[GalTree[q].HaloGalIndex].MergeOn = 3;
#endif


  /* flag galaxy as finished */
  Gal[p].Type = 3;

  /*  calculate mass ratio of merging galaxies */
  mi = Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas;
  ma = Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass+Gal[merger_centralgal].ColdGas;
  if(max(mi,ma) > 0.)
    mass_ratio = min(mi,ma) / max(mi,ma);
  else
    mass_ratio = 1.0;

  /* record the gas and stellar component  mass of merger central and satellite
   * galaxies the before merger */
  Mcstar=(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass);
  Mcbulge=Gal[merger_centralgal].BulgeMass;
  Mcgas=Gal[merger_centralgal].ColdGas;
  Mpstar=(Gal[p].DiskMass+Gal[p].BulgeMass);
  Mpbulge=Gal[p].BulgeMass;
  Mpgas=Gal[p].ColdGas;
#ifdef H2_AND_RINGS
  if((Mcstar-Mcbulge)<1.0e-6) rcstar=0.5*radius[0];
   else
  	{
  	  rcstar=0.5*radius[0]*Gal[merger_centralgal].DiskMassr[0];
  		for(j=1;j<RNUM;j++) rcstar+=(0.5*(radius[j-1]+radius[j])*Gal[merger_centralgal].DiskMassr[j]);
  		rcstar=rcstar/(Mcstar-Mcbulge)/2.0;      //2.0=mean radius/scale length for exponential disk
  	}
  
  if(Mcgas<1.0e-6) rcgas=0.5*radius[0];
   else
  	{
  	  rcgas=0.5*radius[0]*Gal[merger_centralgal].ColdGasr[0];
  		for(j=1;j<RNUM;j++) rcgas+=(0.5*(radius[j-1]+radius[j])*Gal[merger_centralgal].ColdGasr[j]);
  		rcgas=rcgas/Mcgas/2.0;      //2.0=mean radius/scale length for exponential disk
  	}  
#endif
  mass_checks("deal_with_galaxy_merger #1",p);
  mass_checks("deal_with_galaxy_merger #1",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #1",centralgal);



  /* Add the cold and stellar phase of the merged galaxy to the central one.
     Also form a bulge if BulgeFormationInMinorMergersOn is set on (transfer stars
     from satellite disk to central bulge). In a major merger (dealt at the
     make_bulge_from_burst) the disk of the central (now made up of central and
     satellite will be moved to the bulge). Any new stars formed will go to the bulge */

  add_galaxies_together(merger_centralgal, p);

  mass_checks("deal_with_galaxy_merger #2",p);
  mass_checks("deal_with_galaxy_merger #2",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #2",centralgal);

  /* grow black hole through accretion from cold disk during mergers, as in
   * Kauffmann & Haehnelt (2000) + minor mergers - Quasar Mode */
  if(AGNRadioModeModel > 0)
    grow_black_hole(merger_centralgal, mass_ratio, deltaT);

  mass_checks("deal_with_galaxy_merger #3",p);
  mass_checks("deal_with_galaxy_merger #3",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #3",centralgal);

  if (StarBurstRecipe == 0) {
    /* Starburst as in Somerville 2001, with feedback computed inside. */
  	/* All star formation happens in the disk, but in a major merger this will then
  	 * be destroyed with everything moved to the bulge. */
    frac=collisional_starburst_recipe(mass_ratio, merger_centralgal, centralgal, time, deltaT);
    bulgesize_from_merger(mass_ratio,merger_centralgal,p,Mcstar,Mcbulge,Mcgas,
			  Mpstar,Mpbulge,Mpgas,frac);

    mass_checks("deal_with_galaxy_merger #3.5",p);
    mass_checks("deal_with_galaxy_merger #3.5",merger_centralgal);
    mass_checks("deal_with_galaxy_merger #3.5",centralgal);

    if(mass_ratio > ThreshMajorMerger)
      make_bulge_from_burst(merger_centralgal);

  }

  mass_checks("deal_with_galaxy_merger #4",p);
  mass_checks("deal_with_galaxy_merger #4",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #4",centralgal);

  /* If we are in the presence of a minor merger, check disk stability (the disk
   * is completely destroyed in major mergers)*/
  if(TrackDiskInstability)
    if(mass_ratio < ThreshMajorMerger &&
       (Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass) > 0.0)
      check_disk_instability(merger_centralgal);

  /* Not supported option to shrink bulge sizes in gas rich mergers */
#ifdef SHRINKINRICHMERGER  
  if (Mcgas+Mcstar+Mpgas+Mpstar > 1.e-8 ) {
    Gal[merger_centralgal].BulgeSize /= 1+pow((Mcgas+Mpgas)/(Mcgas+Mcstar+Mpgas+Mpstar)/0.15,1.);
    if (Gal[merger_centralgal].BulgeSize < 1.e-8)
      Gal[merger_centralgal].BulgeSize = 1.e-8;
  }
#endif

  if ((Gal[merger_centralgal].BulgeMass > 1.e-6 && Gal[merger_centralgal].BulgeSize == 0.0) ||
      (Gal[merger_centralgal].BulgeMass == 0.0 && Gal[merger_centralgal].BulgeSize >1.e-6)) {
  	char sbuf[1000];
  	sprintf(sbuf, "central: stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f \n",
  			(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass),Gal[merger_centralgal].BulgeMass,
  			Gal[merger_centralgal].BulgeSize,Gal[merger_centralgal].ColdGas,Gal[merger_centralgal].GasDiskRadius,
  			Gal[merger_centralgal].StellarDiskRadius);
  	terminate(sbuf);
  }

#ifndef H2_AND_RINGS
  if (DiskRadiusMethod == 2) {
    get_gas_disk_radius(merger_centralgal);
    get_stellar_disk_radius(merger_centralgal);
  }
#endif

  mass_checks("deal_with_galaxy_merger #5",p);
  mass_checks("deal_with_galaxy_merger #5",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #5",centralgal);

}


/** @brief Grows black holes, through accretion from cold gas during mergers,
 *          as in Kauffmann & Haehnelt (2000) - Quasar Mode.  */

void grow_black_hole(int merger_centralgal, double mass_ratio, double deltaT)
{
  double BHaccrete, fraction;

  /** @brief Grows black hole through accretion from cold gas during mergers,
   *         as in Kauffmann & Haehnelt (2000). I have made an addition here -
   *         the black hole can grow during minor mergers but at a reduced rate
   *         and i have included evolution with redshift.
   *         BlackHoleGrowth == 0 gives instantaneous accretion onto the black hole;
   *         BlackHoleGrowth == 1 instead feeds an accretion disk: accretion occurs
   *                              in main.c */

  if(Gal[merger_centralgal].ColdGas > 0.0)
    {
    BHaccrete = BlackHoleGrowthRate * mass_ratio
      / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[merger_centralgal].Vvir))) * Gal[merger_centralgal].ColdGas;
    /* redshift dependent accretion, not published */
    /* BHaccrete = BlackHoleGrowthRate * (1.0 + ZZ[Halo[halonr].SnapNum]) * mass_ratio */

    /* cannot accrete more gas than is available! */
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;
      
    fraction=BHaccrete/Gal[merger_centralgal].ColdGas;
    if (BlackHoleGrowth == 0) {
      Gal[merger_centralgal].BlackHoleMass += BHaccrete;
      Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;
    } else if (BlackHoleGrowth == 1)
      Gal[merger_centralgal].BlackHoleGas += BHaccrete;
    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas=
      metals_add(Gal[merger_centralgal].MetalsColdGas, Gal[merger_centralgal].MetalsColdGas,-fraction);

#ifdef INDIVIDUAL_ELEMENTS
    Gal[merger_centralgal].ColdGas_elements =
      elements_add(Gal[merger_centralgal].ColdGas_elements, Gal[merger_centralgal].ColdGas_elements,-fraction);
#endif

#ifdef H2_AND_RINGS
    int j;
    for(j=0;j<RNUM;j++)
    {
    	Gal[merger_centralgal].ColdGasr[j] *= (1-fraction);
      Gal[merger_centralgal].MetalsColdGasr[j] *= (1-fraction);    
    }		
#endif

  }
}

/** @brief Adds all the components of the satellite galaxy into its
 *         central companion. */

void add_galaxies_together(int t, int p)
{
  /** @brief All the components of the satellite galaxy are added to the
   *         correspondent component of the central galaxy. Cold gas spin
   *         is updated and a bulge is formed at the central galaxy, with
   *         the stars of the satellite if  BulgeFormationInMinorMergersOn=1.
   *         In case of a major merger, everything that was put in the disk of
   *         the central galaxy will be moved into the bulge
   *
   * TODO Even though galaxy p has been set to type 3 (ie a non-galaxy), it would
   * make mass conservation more explicit to zero the properties of galaxy p after
   * the merger.
   * TODO Correct artificial diffusion of metals when BulgeFormationInMinorMergersOn=1. */
  int outputbin, j;
  float tspin[3],tmass,pmass;

  /* t central, p satellite */

  mass_checks("add_galaxies_together #0",p);
  mass_checks("add_galaxies_together #0.1",t);

  /* angular momentum transfer between gas*/
  tmass= Gal[t].ColdGas;
  pmass= Gal[p].ColdGas;
  
  Gal[t].MergeSat +=(Gal[p].DiskMass+Gal[p].BulgeMass);
  Gal[p].MergeSat=0.;

#ifndef H2_AND_RINGS
  transfer_gas(t,"Cold",p,"Cold",1.,"add_galaxies_together", __LINE__);
  transfer_gas(t,"Hot",p,"Hot",1.,"add_galaxies_together", __LINE__);
  transfer_gas(t,"Ejected",p,"Ejected",1.,"add_galaxies_together", __LINE__); //TODO chose move to ejected or hot
#ifdef TRACK_BURST
    /* The whole burst component gets transferred */
  transfer_stars(t,"Burst",p,"Burst",1.);
#endif


  if(BulgeFormationInMinorMergersOn) 
    transfer_stars(t,"Bulge",p,"Disk",1.);
  else
    transfer_stars(t,"Disk",p,"Disk",1.);
  transfer_stars(t,"Bulge",p,"Bulge",1.);
  transfer_stars(t,"ICM",p,"ICM",1.);

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;    
  Gal[p].BlackHoleMass=0.;
  Gal[t].BlackHoleGas += Gal[p].BlackHoleGas;    
  Gal[p].BlackHoleGas=0.;
  Gal[t].StarMerge += Gal[p].StarMerge;
  Gal[p].StarMerge=0.;

#else //ifdef H2_AND_RINGS
  /*satellite cold gas disk directly superpose to the central gas disk*/
//  for(j=0;j<RNUM;j++)
//  {
//  	Gal[t].ColdGasr[j] += Gal[p].ColdGasr[j];
//    Gal[t].MetalsColdGasr[j] += Gal[p].MetalsColdGasr[j];
//  }

/*Satellite cold gas mixed, and then accerated to the central galaxy like the infall process */ 
  double ringtot, ringfracr[RNUM], rd, cosangle;
  cosangle=(Gal[t].DiskSpin[1]*Gal[t].HaloSpin[1]+Gal[t].DiskSpin[2]*Gal[t].HaloSpin[2]+Gal[t].DiskSpin[0]*Gal[t].HaloSpin[0])/
        sqrt(Gal[t].DiskSpin[1]*Gal[t].DiskSpin[1]+Gal[t].DiskSpin[2]*Gal[t].DiskSpin[2]+Gal[t].DiskSpin[0]*Gal[t].DiskSpin[0])/
        sqrt(Gal[t].HaloSpin[1]*Gal[t].HaloSpin[1]+Gal[t].HaloSpin[2]*Gal[t].HaloSpin[2]+Gal[t].HaloSpin[0]*Gal[t].HaloSpin[0]);
  
  /*cosangle=( (Gal[t].GasSpin[0]+Gal[t].StellarSpin[0])*Gal[t].HaloSpin[1] +
  		       (Gal[t].GasSpin[1]+Gal[t].StellarSpin[1])*Gal[t].HaloSpin[2] +
  		       (Gal[t].GasSpin[2]+Gal[t].StellarSpin[2])*Gal[t].HaloSpin[0] ) /
           sqrt( (Gal[t].GasSpin[0]+Gal[t].StellarSpin[0])*(Gal[t].GasSpin[0]+Gal[t].StellarSpin[0]) +
          		   (Gal[t].GasSpin[1]+Gal[t].StellarSpin[1])*(Gal[t].GasSpin[1]+Gal[t].StellarSpin[1]) +
          		   (Gal[t].GasSpin[2]+Gal[t].StellarSpin[2])*(Gal[t].GasSpin[2]+Gal[t].StellarSpin[2]) ) /
           sqrt(Gal[t].HaloSpin[0]*Gal[t].HaloSpin[0]+Gal[t].HaloSpin[1]*Gal[t].HaloSpin[1]+Gal[t].HaloSpin[2]*Gal[t].HaloSpin[2]);*/
 // the cos of the angle between disk and halo
  cosangle=fabs(cosangle);
  if(cosangle>1.0) cosangle=1.0;
  if(cosangle<0.2) cosangle=0.2; //minimum inclination
  
  rd=Gal[t].GasDiskRadius/3.0;

	ringtot=1-(1+radius[RNUM-1]/rd)/exp(radius[RNUM-1]/rd);
	ringfracr[0]=(1-(1+radius[0]/rd)/exp(radius[0]/rd))/ringtot;
	for(j=1; j<RNUM; j++)
	{ringfracr[j]=((1+radius[j-1]/rd)/exp(radius[j-1]/rd)-(1+radius[j]/rd)/exp(radius[j]/rd))/ringtot;}  

	 for(j=0; j<RNUM; j++)
  {
    Gal[t].ColdGasr[j]+= Gal[p].ColdGas*ringfracr[j];
  	Gal[t].MetalsColdGasr[j] += Gal[p].MetalsColdGas*ringfracr[j];  //metal evenly distributed 
  }
  
  Gal[t].ColdGas += Gal[p].ColdGas;
  Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;
  Gal[t].StarMerge +=Gal[p].StarMerge;
  Gal[t].MergeSat +=(Gal[p].DiskMass+Gal[p].BulgeMass);  
  Gal[t].ICM += Gal[p].ICM;
  Gal[t].MetalsICM +=Gal[p].MetalsICM;

  if(BulgeFormationInMinorMergersOn) 
  {
  	  Gal[t].BulgeMass += (Gal[p].DiskMass+Gal[p].BulgeMass);
      Gal[t].MetalsBulgeMass += (Gal[p].MetalsDiskMass+Gal[p].MetalsBulgeMass);
  }
  else
  	{
  		Gal[t].DiskMass += Gal[p].DiskMass;
      Gal[t].MetalsDiskMass += Gal[p].MetalsDiskMass;
      Gal[t].BulgeMass += Gal[p].BulgeMass;
      Gal[t].MetalsBulgeMass += Gal[p].MetalsBulgeMass;
  		for(j=0;j<RNUM;j++)
  		{
  			Gal[t].DiskMassr[j] += Gal[p].DiskMassr[j];
        Gal[t].MetalsDiskMassr[j] += Gal[p].MetalsDiskMassr[j];
#ifdef SAVE_MEMORY
        Gal[t].Sfrr[j] += Gal[p].Sfrr[j];
#else
        for(outputbin=0; outputbin<NOUT; outputbin++) Gal[t].Sfrr[j][outputbin] += Gal[p].Sfrr[j][outputbin];        
#endif
  		}
  	}
  // During merger the angular momentum of reheated gas for satelite will not transfer to	the central galaxy.
  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;
  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;    
  Gal[t].BlackHoleGas += Gal[p].BlackHoleGas;    
#endif //H2_AND_RINGS
  mass_checks("add_galaxies_together #1",p);
  mass_checks("add_galaxies_together #1.1",t);

  /*update the gas spin*/
#ifndef H2_AND_RINGS
  for(j=0;j<3;j++)
    tspin[j]=Gal[t].GasSpin[j]*tmass+Gal[t].HaloSpin[j]*pmass;
  if (Gal[t].ColdGas != 0)
    for (j=0;j<3;j++)
     Gal[t].GasSpin[j]=tspin[j]/(Gal[t].ColdGas); 
#else
  for(j=0;j<3;j++)
     tspin[j]=Gal[t].DiskSpin[j]*tmass+Gal[t].HaloSpin[j]*pmass;
   if ((Gal[t].DiskMass+Gal[t].ColdGas)>0.0)
   {
     for (j=0;j<3;j++) Gal[t].DiskSpin[j]=tspin[j]/(Gal[t].DiskMass+Gal[t].ColdGas);
   }
#endif

#ifdef SAVE_MEMORY
  Gal[t].Sfr += Gal[p].Sfr;
#else
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    Gal[t].Sfr[outputbin] += Gal[p].Sfr[outputbin];
#endif

  if(BulgeFormationInMinorMergersOn) {
#ifdef SAVE_MEMORY
    Gal[t].SfrBulge += Gal[p].Sfr;
#else
    for(outputbin = 0; outputbin < NOUT; outputbin++)
      Gal[t].SfrBulge[outputbin] += Gal[p].Sfr[outputbin];
#endif
  }

  /* Add the black hole accretion rates.  This makes little sense but is not
   * used if the superior BlackHoleGrowth==1 switch is on. */
  Gal[t].QuasarAccretionRate += Gal[p].QuasarAccretionRate;
  Gal[t].RadioAccretionRate += Gal[p].RadioAccretionRate;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
  	 Gal[t].MassWeightAge[outputbin] += Gal[p].MassWeightAge[outputbin];

#ifndef  POST_PROCESS_MAGS
/* Add the luminosities of the satellite and central galaxy */
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {

    for(j = 0; j < NMAG; j++) {
      Gal[t].Lum[j][outputbin] += Gal[p].Lum[j][outputbin];
      Gal[t].YLum[j][outputbin] += Gal[p].YLum[j][outputbin];
#ifdef ICL
      Gal[t].ICLLum[j][outputbin] += Gal[p].ICLLum[j][outputbin];
#endif
    }
    if(BulgeFormationInMinorMergersOn) {
      for(j = 0; j < NMAG; j++) {
	Gal[t].LumBulge[j][outputbin] += Gal[p].Lum[j][outputbin];
	Gal[t].YLumBulge[j][outputbin] += Gal[p].YLum[j][outputbin];
      }
    }
    else {
      for(j = 0; j < NMAG; j++) {
	Gal[t].LumBulge[j][outputbin]  += Gal[p].LumBulge[j][outputbin];
	Gal[t].YLumBulge[j][outputbin] += Gal[p].YLumBulge[j][outputbin];
      }
    }

  }
#endif // OUTPUT_REST_MAGS 

#ifdef COMPUTE_OBS_MAGS 
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++) {
      Gal[t].ObsLum[j][outputbin]   += Gal[p].ObsLum[j][outputbin];
      Gal[t].ObsYLum[j][outputbin]  += Gal[p].ObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].ObsICL[j][outputbin]  += Gal[p].ObsICL[j][outputbin];
#endif

#ifdef OUTPUT_MOMAF_INPUTS
      Gal[t].dObsLum[j][outputbin] += Gal[p].dObsLum[j][outputbin];
      Gal[t].dObsYLum[j][outputbin] += Gal[p].dObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].dObsICL[j][outputbin]  += Gal[p].dObsICL[j][outputbin];
#endif
#endif
    }
    if(BulgeFormationInMinorMergersOn) {
      for(j = 0; j < NMAG; j++) {
	Gal[t].ObsLumBulge[j][outputbin]   += Gal[p].ObsLum[j][outputbin];
	Gal[t].ObsYLumBulge[j][outputbin]  += Gal[p].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	Gal[t].dObsLumBulge[j][outputbin]  += Gal[p].dObsLum[j][outputbin];
	Gal[t].dObsYLumBulge[j][outputbin] += Gal[p].dObsYLum[j][outputbin];
#endif
      }
    }
    else 
    {
      for(j = 0; j < NMAG; j++) {
	Gal[t].ObsLumBulge[j][outputbin]   += Gal[p].ObsLumBulge[j][outputbin];
	Gal[t].ObsYLumBulge[j][outputbin]  += Gal[p].ObsYLumBulge[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	Gal[t].dObsLumBulge[j][outputbin]  += Gal[p].dObsLumBulge[j][outputbin];
	Gal[t].dObsYLumBulge[j][outputbin] += Gal[p].dObsYLumBulge[j][outputbin];
#endif
      }
    }
  }
#endif //COMPUTE_OBS_MAGS 
#endif  //POST_PROCESS_MAGS
}


/** @brief In a major merger, both disks are destroyed and all the mass transferred
 *         to the bulge. The galaxies have already been merged, so all we need to do here
 *         is transfer stars from disk to bulge.
 *  TODO Should we also transfer cold to hot gas? */

void make_bulge_from_burst(int p) 
{
	int outputbin;
  /* generate bulge */
  transfer_stars(p,"Bulge",p,"Disk",1.);

  /*  update the star formation rate */
#ifdef SAVE_MEMORY
  Gal[p].SfrBulge  = Gal[p].Sfr;
#else
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    Gal[p].SfrBulge[outputbin] += Gal[p].Sfr[outputbin];
#endif

#ifndef  POST_PROCESS_MAGS
  int j;
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].LumBulge[j][outputbin]  = Gal[p].Lum[j][outputbin];
      Gal[p].YLumBulge[j][outputbin] = Gal[p].YLum[j][outputbin];
    }
  }
#endif
#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].ObsLumBulge[j][outputbin]   = Gal[p].ObsLum[j][outputbin];
      Gal[p].ObsYLumBulge[j][outputbin]  = Gal[p].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
      Gal[p].dObsLumBulge[j][outputbin]  = Gal[p].dObsLum[j][outputbin];
      Gal[p].dObsYLumBulge[j][outputbin] = Gal[p].dObsYLum[j][outputbin];
#endif
    }
  }
#endif
#endif //POST_PROCESS_MAGS
}

/** @brief Merger burst recipe from Somerville 2001 (used after Croton2006) */

double collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal,
				  double time, double deltaT)
{
  /** @brief If StarBurstRecipe = 1 (since Croton2006), the Somerville 2001
   *         model of bursts is used. The burst can happen for both major
   *         and minor mergers, with a fraction of the added cold gas from
   *         the satellite and central being consumed. SN Feedback from
   *         starformation is computed and the sizes of bulge and disk
   *         followed (not done for the other burst mode).*/

  double mstars, metallicitySF, eburst, Ggas;
#ifdef H2_AND_RINGS
  double mstarsr[RNUM], metallicityr[RNUM];
   int j;
#endif
  /* This is the major and minor merger starburst recipe of Somerville 2001.
   * The coefficients in eburst are taken from TJ Cox's PhD thesis and should
   * be more accurate then previous. */

  Ggas=Gal[merger_centralgal].ColdGas;

  /* the bursting fraction given the mass ratio */
  /* m_dot = 0.56*(m_sat/m_central)^0.7*m_gas */
  eburst = SfrBurstEfficiency * pow(mass_ratio, SfrBurstSlope);
  //eburst = 0.56 * pow(mass_ratio, 0.7);
  mstars = eburst * Gal[merger_centralgal].ColdGas;
  if(mstars < 0.0)
    mstars = 0.0;

#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
  {
  	mstarsr[j] = eburst * Gal[merger_centralgal].ColdGasr[j];
  	if(mstarsr[j] < 0.0) mstarsr[j] = 0.0;
  }
#endif

#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN //otherwise there is another check inside SN_feedback
  if(mstars > Gal[merger_centralgal].ColdGas)
        mstars = Gal[merger_centralgal].ColdGas;
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
  	if(mstarsr[j] > Gal[merger_centralgal].ColdGasr[j])
  	        mstarsr[j] = Gal[merger_centralgal].ColdGasr[j];
#endif
#endif


   
  /*  update the star formation rate */
#ifdef SAVE_MEMORY
  Gal[merger_centralgal].Sfr += mstars / deltaT;
#ifdef H2_AND_RINGS
   for(j=0;j<RNUM;j++)
  		 Gal[merger_centralgal].Sfrr[j] += mstarsr[j] / deltaT;
#endif
#else
  int outputbin;
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    if(Gal[merger_centralgal].SnapNum == ListOutputSnaps[outputbin]) {
      Gal[merger_centralgal].Sfr[outputbin] += mstars / deltaT;
#ifdef H2_AND_RINGS
      for(j=0;j<RNUM;j++) Gal[merger_centralgal].Sfrr[j][outputbin] += mstarsr[j] / deltaT;
#endif
      break;
    }
  }
#endif

  /* Store the value of the metallicity of the cold phase when SF occurs.
   * Used to update luminosities below */
  metallicitySF = metals_total(Gal[merger_centralgal].MetalsColdGas)/Gal[merger_centralgal].ColdGas;

  //feedback is called before starformationto change the value of mstars so that cold gas is
  //shared between star formation and feedback

	if (mstars > 0.)
#ifndef H2_AND_RINGS
	update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars);
#else
	update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars, mstarsr);
#endif

  int nstep=-1;
	if (mstars > 0.)
#ifndef H2_AND_RINGS
		update_from_star_formation(merger_centralgal, mstars, true, nstep); // true indicates starburst
#else
		update_from_star_formation(merger_centralgal, mstars, mstarsr, true, nstep); // true indicates starburst
#endif

	mass_checks("collisional_starburst_recipe #2",merger_centralgal);

	update_massweightage(merger_centralgal, mstars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if (mstars > 0.)
#ifndef H2_AND_RINGS
  	SN_feedback(merger_centralgal, centralgal, mstars, "ColdGas");
#else
  	SN_feedback(merger_centralgal, centralgal, mstars, mstarsr, "ColdGas");
#endif
#endif

    mass_checks("collisional_starburst_recipe #1",merger_centralgal);


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  update the luminosities due to the stars formed */
  if (mstars > 0.0) 
    add_to_luminosities(merger_centralgal, mstars, time, metallicitySF);
#endif //NDEF POST_PROCESS_MAGS
#endif

  if (Ggas > 0.)
    return mstars/Ggas;
  else
    return 0.0;





}

/** @brief Calculates the bulge size after a merger. */
void  bulgesize_from_merger(double mass_ratio,int merger_centralgal,int p,
			    double Mcstar,double Mcbulge,double Mcgas,
			    double Mpstar,double Mpbulge,double Mpgas,double frac)
{
  /** @brief For any type of merger calculates the new bulge size using
   *         Eq. 33 in Guo2010:
   *
   *         \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
   *             C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+		\
   *             alpha_{\rm{inter}}\frac{GM_1M_2}{R_1+R_2}\f$.
   *
   *         This implementation assumed that the new bulge occupies the
   *         same space as the components that formed it. */

  double Mc,Rc;
  double Mp,Rp;
  double fint,c;
#ifdef H2_AND_RINGS
  double rpstar,rpgas;
  int j;
#endif
  fint=0.5;
  c=0.5;

#ifdef H2_AND_RINGS
 Mpgas=Gal[p].ColdGas;
  Mpbulge=Gal[p].BulgeMass;
  Mpstar=Gal[p].DiskMass+Gal[p].BulgeMass;
  
  
  if((Mpstar-Mpbulge)<1.0e-6) rpstar=0.5*radius[0];
   else
  	{
  	  rpstar=0.5*radius[0]*Gal[p].DiskMassr[0];
  		for(j=1;j<RNUM;j++) rpstar+=(0.5*(radius[j-1]+radius[j])*Gal[p].DiskMassr[j]);
  		rpstar=rpstar/(Mpstar-Mpbulge)/2.0;      //2.0=mean radius/scale length for exponential disk
  	}
  	
  if(Mpgas<1.0e-6) rpgas=0.5*radius[0];
   else
  	{
  	  rpgas=0.5*radius[0]*Gal[p].ColdGasr[0];
  		for(j=1;j<RNUM;j++) rpgas+=(0.5*(radius[j-1]+radius[j])*Gal[p].ColdGasr[j]);
  		rpgas=rpgas/Mpgas/2.0;      //2.0=mean radius/scale length for exponential disk
  	}  	
#endif
  /* calculate radius for the object that will form the new bulge - Rc and Rp */
  /* Minor Merger */
  if(mass_ratio < ThreshMajorMerger) {
    /* In a minor merger only the stars of the satellite galaxy are moved
     * to the bulge of the central galaxy, therefore only stellar
     * components are used to compute radius and masses. */
    frac=0.0;
    /* in a minor merger only consider the bulge mass of the central galaxy */
    Mc=Mcbulge;
    Rc=Gal[merger_centralgal].BulgeSize;
    /* and stellarmass of satellite*/
    Mp=Mpstar;
    if (Mp >0.0)
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge)/Mpstar;
    else 
      Rp=0.0;
  }
  /* Major Merger */
  else {
    /* on a major merger both stellar and gas (after a burst) components
     * from the final bulge and need to be considered */
    /* Mc = bulge mass + burst of central*/
    Mc=Mcstar+frac*Mcgas;
    if (Mc > 0.0)
      Rc=(Gal[merger_centralgal].StellarDiskRadius/3.*1.68*(Mcstar-Mcbulge)+Gal[merger_centralgal].BulgeSize*Mcbulge+Gal[merger_centralgal].GasDiskRadius*frac*Mcgas/3.*1.68)/(Mcgas*frac+Mcstar);
    else
      Rc=0.0;
    /* and satellite Mp */
    Mp=Mpstar+frac*Mpgas;
    if (Mp > 0.0)
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge+Gal[p].GasDiskRadius*frac*Mpgas/3.*1.68)/(Mpgas*frac+Mpstar);
    else
      Rp=0.0;
  }

  /* If both original radius are bigger then 0 then this is Eq. 33 in Guo 2010
   * solved for R_new,bulge with all terms divided by G and C. */
  if (Rp >1.e-8 && Rc >1.e-8 )
    Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));
    
  if (Rc >1.e-8 && Rp <1.e-8)
    Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));
  
  if(Rc <1.e-8 && Rp <1.e-8)
    Gal[merger_centralgal].BulgeSize=0.0;
  
  if (Rc <1.e-8 && Rp > 1.e-8)
    Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+fint/c*Mp*Mc/(Rp+Rc));
      

  if ((Mp+Mc > 0.0 && Gal[merger_centralgal].BulgeSize == 0.0 )||(Mp+Mc == 0.0 && Gal[merger_centralgal].BulgeSize> 0.0)) {
  	char sbuf[1000];
  	  	sprintf(sbuf, "bulgesize wrong in merger. bulgemass %f, bulgesize %f, Rp %f, Rc %f,Mp %f,Mc %f, mass ratio %f, halonr %d, merger_centralgal %d\n,   p: stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f \n,   central: stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f \n",
  	  			Gal[merger_centralgal].BulgeMass, Gal[p].BulgeSize, Rp,Rc,Mp,Mc, mass_ratio, Gal[merger_centralgal].HaloNr,
  	  			merger_centralgal, (Gal[p].DiskMass+Gal[p].BulgeMass),Gal[p].BulgeMass, Gal[p].BulgeSize,Gal[p].ColdGas,
  	  			Gal[p].GasDiskRadius,Gal[p].StellarDiskRadius, (Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass),
  	  			Gal[merger_centralgal].BulgeMass,  Gal[merger_centralgal].BulgeSize,Gal[merger_centralgal].ColdGas,
  	  			Gal[merger_centralgal].GasDiskRadius, Gal[merger_centralgal].StellarDiskRadius);
  	  	terminate(sbuf);
   exit(0);
  }

}
